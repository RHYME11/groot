#include <Plugin/GPluginManager.h>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <set>
#include <sstream>

#include <Plugin/GPlugin.h>
#include <Plugin/GPluginEvent.h>
#include <Plugin/GPluginSession.h>
#include <Plugin/GPluginVersion.h>

#include "GPluginLoader.h"
#include "GPluginManifest.h"
#include "GPluginRegistry.h"
#include "GPluginSessionRegistry.h"

namespace {

// ============== AppendSearchPath ==============
// Purpose: Append one unique non-empty plugin search directory.
// Inputs: Candidate path, ordered output, and duplicate tracker.
// Outputs: Updated output containers.
void AppendSearchPath(const std::string& path,
                      std::vector<std::string>& paths,
                      std::set<std::string>& seen) {
  if(path.empty())
    return;
  const std::string normalized =
    std::filesystem::path(path).lexically_normal().string();
  if(seen.insert(normalized).second)
    paths.push_back(normalized);
}

} // namespace

// ============== GPluginManager::Get ==============
// Purpose: Access the process-wide Groot plugin manager.
// Inputs: None.
// Outputs: Singleton manager reference.
GPluginManager& GPluginManager::Get() {
  static GPluginManager* manager = new GPluginManager;
  return *manager;
}

// ============== GPluginManager::GPluginManager ==============
// Purpose: Create an empty process-wide plugin registry.
// Inputs: None.
// Outputs: Uninitialized manager.
GPluginManager::GPluginManager()
  : fRegistry(new GPluginRegistry),
    fSessions(new GPluginSessionRegistry),
    fInitialized(false),
    fShutdown(false) {
}

// ============== GPluginManager::~GPluginManager ==============
// Purpose: Destroy plugin instances and registry metadata at shutdown.
// Inputs: None.
// Outputs: Released manager resources; libraries remain process-loaded.
GPluginManager::~GPluginManager() {
  Shutdown();
  delete fSessions;
  delete fRegistry;
}

// ============== GPluginManager::Shutdown ==============
// Purpose: Close sessions and destroy plugin instances before ROOT cleanup.
// Inputs: None.
// Outputs: Manager enters an inert state; repeated calls are safe.
void GPluginManager::Shutdown() {
  if(fShutdown)
    return;
  fShutdown = true;
  fActionChangedCallback = {};
  fContextProvider = {};
  fStatusCallback = {};
  for(const auto& record : fSessions->TakeAll()) {
    if(record.session)
      record.session->Close();
  }
  fRegistry->Shutdown();
}

// ============== GPluginManager::Initialize ==============
// Purpose: Discover plugin manifests once during application startup.
// Inputs: None.
// Outputs: Populated action registry.
void GPluginManager::Initialize() {
  if(fInitialized || fShutdown)
    return;
  fInitialized = true;
  Rescan();
}

// ============== GPluginManager::Rescan ==============
// Purpose: Rescan configured plugin directories without unloading libraries.
// Inputs: None.
// Outputs: Refreshed action registry and GUI notification.
void GPluginManager::Rescan() {
  if(fShutdown)
    return;
  fRegistry->ClearMetadata();
  for(const auto& path : SearchPaths())
    ScanDirectory(path);
  if(fActionChangedCallback)
    fActionChangedCallback();
}

// ============== GPluginManager::Actions ==============
// Purpose: Return all actions declared by accepted manifests.
// Inputs: None.
// Outputs: Read-only action list.
const std::vector<GPluginAction>& GPluginManager::Actions() const {
  return fRegistry->Actions();
}

// ============== GPluginManager::ExecuteAction ==============
// Purpose: Lazy-load an action owner and execute the selected action.
// Inputs: Global action identifier.
// Outputs: True when the plugin reports success.
bool GPluginManager::ExecuteAction(const std::string& actionId) {
  if(fShutdown)
    return false;
  GPluginRecord* record = fRegistry->FindByAction(actionId);
  if(!record) {
    ReportError("unknown plugin action '" + actionId + "'");
    return false;
  }
  if(record->manifest.apiVersion != kGPluginApiVersion) {
    ReportError("plugin '" + record->manifest.id + "' requires API " +
                std::to_string(record->manifest.apiVersion) + ", Groot provides " +
                std::to_string(kGPluginApiVersion));
    return false;
  }

  std::filesystem::path library(record->manifest.library);
  if(library.is_relative())
    library = std::filesystem::path(record->manifest.directory) / library;

  std::string error;
  if(!GPluginLoader::Load(library.lexically_normal().string(), this,
                          record->loaded, error)) {
    ReportError("plugin '" + record->manifest.id + "': " + error);
    return false;
  }

  GPluginContext context;
  if(fContextProvider)
    context = fContextProvider();
  if(!record->loaded.instance->ExecuteAction(actionId.c_str(), context)) {
    const char* detail = record->loaded.instance->LastError();
    ReportError("plugin action '" + actionId + "' failed" +
                (detail && *detail ? std::string(": ") + detail : ""));
    return false;
  }
  SetStatusMessage(("Plugin action completed: " + actionId).c_str());
  return true;
}

// ============== GPluginManager::CleanArtifacts ==============
// Purpose: Ask loaded plugins to remove artifacts for one ROOT context.
// Inputs: Current canvas, pad, selection, and target.
// Outputs: Cleanup failures are reported without interrupting other plugins.
void GPluginManager::CleanArtifacts(const GPluginContext& context) {
  if(fShutdown)
    return;
  for(GPluginRecord* record : fRegistry->LoadedPlugins()) {
    try {
      if(!record->loaded.instance->CleanArtifacts(context)) {
        const char* detail = record->loaded.instance->LastError();
        ReportError("plugin '" + record->manifest.id +
                    "' artifact cleanup failed" +
                    (detail && *detail ? std::string(": ") + detail : ""));
      }
    } catch(const std::exception& error) {
      ReportError("plugin '" + record->manifest.id +
                  "' artifact cleanup threw: " + error.what());
    } catch(...) {
      ReportError("plugin '" + record->manifest.id +
                  "' artifact cleanup threw an unknown exception");
    }
  }
}

// ============== GPluginManager::SetActionChangedCallback ==============
// Purpose: Attach the GUI action-registry refresh callback.
// Inputs: Optional callback.
// Outputs: Replaced action callback.
void GPluginManager::SetActionChangedCallback(std::function<void()> callback) {
  fActionChangedCallback = std::move(callback);
}

// ============== GPluginManager::SetContextProvider ==============
// Purpose: Attach the host callback that describes the current GUI context.
// Inputs: Optional context provider.
// Outputs: Replaced context provider.
void GPluginManager::SetContextProvider(
  std::function<GPluginContext()> provider) {
  fContextProvider = std::move(provider);
}

// ============== GPluginManager::SetStatusCallback ==============
// Purpose: Attach the optional GUI status-message callback.
// Inputs: Optional status callback.
// Outputs: Replaced status callback.
void GPluginManager::SetStatusCallback(
  std::function<void(const std::string&)> callback) {
  fStatusCallback = std::move(callback);
}

// ============== GPluginManager::SetStatusMessage ==============
// Purpose: Publish a plugin or manager status message.
// Inputs: Null-terminated message.
// Outputs: Message printed and forwarded to the GUI when available.
void GPluginManager::SetStatusMessage(const char* message) {
  const std::string text = message ? message : "";
  if(!text.empty())
    std::printf("groot plugin: %s\n", text.c_str());
  if(fStatusCallback)
    fStatusCallback(text);
}

// ============== GPluginManager::ActivateSession ==============
// Purpose: Register one application-overlay session for a pad and target.
// Inputs: Plugin session and validated context.
// Outputs: True when the pad was inactive and registration succeeded.
bool GPluginManager::ActivateSession(GPluginSession* session,
                                     const GPluginContext& context) {
  if(fShutdown)
    return false;
  if(!session || !context.canvas || !context.pad) {
    ReportError("cannot activate a plugin session without canvas and pad");
    return false;
  }
  GPluginSessionRecord record;
  record.canvas = context.canvas;
  record.pad = context.pad;
  record.target = context.target;
  record.session = session;
  if(!fSessions->Add(record)) {
    std::string reason = "the selected pad already has an active plugin session";
    GPluginSession* owner = nullptr;
    if(auto* padRecord = fSessions->Find(context.pad))
      owner = padRecord->session;
    if(!fSessions->Contains(context.pad) && context.target &&
       fSessions->ContainsTarget(context.target)) {
      reason = "the selected target already has an active plugin session";
      owner = fSessions->FindTarget(context.target);
    }
    if(owner && owner->SessionId() && owner->SessionId()[0] != '\0')
      reason += std::string(": ") + owner->SessionId();
    ReportError(reason);
    return false;
  }
  SetStatusMessage("Plugin session activated");
  return true;
}

// ============== GPluginManager::DeactivateSession ==============
// Purpose: Release one plugin session's pad and target ownership.
// Inputs: Active session identity.
// Outputs: Updated session registry.
void GPluginManager::DeactivateSession(GPluginSession* session) {
  if(fShutdown)
    return;
  fSessions->Remove(session);
  SetStatusMessage("Plugin session deactivated");
}

// ============== GPluginManager::HasActiveSession ==============
// Purpose: Test whether a pad is owned by a plugin session.
// Inputs: Pad identity.
// Outputs: True when Groot interaction must be bypassed.
bool GPluginManager::HasActiveSession(TVirtualPad* pad) const {
  if(fShutdown)
    return false;
  return fSessions->Contains(pad);
}

// ============== GPluginManager::NotifySession ==============
// Purpose: Notify a session after ROOT has processed a native input event.
// Inputs: Normalized event with ROOT's current selection.
// Outputs: True whenever the pad has an application-overlay session.
bool GPluginManager::NotifySession(const GPluginEvent& event) {
  if(fShutdown)
    return false;
  auto* record = fSessions->Find(event.context.pad);
  if(!record || !record->session)
    return false;
  GPluginEvent observed = event;
  observed.context.canvas = record->canvas;
  observed.context.pad = record->pad;
  observed.context.target = record->target;
  record->session->ObserveEvent(observed);
  return true;
}

// ============== GPluginManager::CloseSessionsForPad ==============
// Purpose: Close and detach the session for one pad.
// Inputs: Pad being destroyed or reset.
// Outputs: Closed session callbacks.
void GPluginManager::CloseSessionsForPad(TVirtualPad* pad) {
  if(fShutdown)
    return;
  for(const auto& record : fSessions->TakePad(pad)) {
    if(record.session)
      record.session->Close();
  }
}

// ============== GPluginManager::CloseSessionsForCanvas ==============
// Purpose: Close and detach every session owned by one canvas.
// Inputs: Canvas being closed.
// Outputs: Closed session callbacks.
void GPluginManager::CloseSessionsForCanvas(TCanvas* canvas) {
  if(fShutdown)
    return;
  for(const auto& record : fSessions->TakeCanvas(canvas)) {
    if(record.session)
      record.session->Close();
  }
}

// ============== GPluginManager::SearchPaths ==============
// Purpose: Build the ordered plugin directory search list.
// Inputs: Process environment.
// Outputs: Unique normalized search paths.
std::vector<std::string> GPluginManager::SearchPaths() const {
  std::vector<std::string> paths;
  std::set<std::string> seen;

  if(const char* configured = std::getenv("GROOT_PLUGIN_PATH")) {
    std::stringstream stream(configured);
    std::string item;
    while(std::getline(stream, item, ':'))
      AppendSearchPath(item, paths, seen);
  }
  if(const char* home = std::getenv("HOME"))
    AppendSearchPath(std::string(home) + "/.local/lib/groot/plugins", paths, seen);
  if(const char* gsys = std::getenv("GSYS"))
    AppendSearchPath(std::string(gsys) + "/lib/groot/plugins", paths, seen);
  return paths;
}

// ============== GPluginManager::ScanDirectory ==============
// Purpose: Read sorted .plugin manifests from one directory.
// Inputs: Directory path.
// Outputs: Accepted manifests added to the registry.
void GPluginManager::ScanDirectory(const std::string& path) {
  std::error_code errorCode;
  if(!std::filesystem::is_directory(path, errorCode))
    return;

  std::vector<std::filesystem::path> manifests;
  for(const auto& entry : std::filesystem::directory_iterator(path, errorCode)) {
    if(errorCode)
      break;
    if(entry.is_regular_file() && entry.path().extension() == ".plugin")
      manifests.push_back(entry.path());
  }
  std::sort(manifests.begin(), manifests.end());

  for(const auto& pathEntry : manifests) {
    GPluginManifestData manifest;
    std::string error;
    if(!GPluginManifest::Read(pathEntry.string(), manifest, error)) {
      ReportError("manifest '" + pathEntry.string() + "': " + error);
      continue;
    }
    if(!fRegistry->Add(manifest, error))
      ReportError("manifest '" + pathEntry.string() + "': " + error);
  }
}

// ============== GPluginManager::ReportError ==============
// Purpose: Report a recoverable plugin error without stopping Groot.
// Inputs: Error detail.
// Outputs: Terminal and optional GUI status message.
void GPluginManager::ReportError(const std::string& message) {
  std::fprintf(stderr, "groot plugin ERROR: %s\n", message.c_str());
  if(fStatusCallback)
    fStatusCallback(message);
}
