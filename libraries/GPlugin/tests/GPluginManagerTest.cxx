#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

#include <Plugin/GPluginManager.h>
#include <Plugin/GPluginEvent.h>
#include <Plugin/GPluginSession.h>
#include <Plugin/GPluginVersion.h>

namespace {

int failures = 0;

class MockSession : public GPluginSession {
  public:
    const char* SessionId() const override {
      return "mock.session";
    }
    void ObserveEvent(const GPluginEvent& event) override {
      ++events;
      lastTarget = event.context.target;
    }
    void Close() override {
      closed = true;
    }
    int events = 0;
    bool closed = false;
    TObject* lastTarget = nullptr;
};

// ============== Check ==============
// Purpose: Record one integration-test expectation.
// Inputs: Condition and failure description.
// Outputs: Updated global failure count.
void Check(bool condition, const std::string& message) {
  if(condition)
    return;
  std::cerr << "FAIL: " << message << std::endl;
  ++failures;
}

// ============== WriteManifest ==============
// Purpose: Write one test TEnv plugin manifest.
// Inputs: Path, identity, library, action, and API version.
// Outputs: Manifest file in the test runtime directory.
void WriteManifest(const std::filesystem::path& path,
                   const std::string& id,
                   const std::string& library,
                   const std::string& action,
                   int apiVersion = static_cast<int>(kGPluginApiVersion)) {
  std::ofstream output(path);
  output << "Plugin.Id: " << id << '\n'
         << "Plugin.Name: " << id << '\n'
         << "Plugin.Version: 1.0.0\n"
         << "Plugin.ApiVersion: " << apiVersion << '\n'
         << "Plugin.Library: " << library << '\n'
         << "Plugin.Action.Count: 1\n"
         << "Plugin.Action.0.Id: " << action << '\n'
         << "Plugin.Action.0.Label: " << action << '\n'
         << "Plugin.Action.0.Tooltip: integration test\n";
}

} // namespace

// ============== main ==============
// Purpose: Verify plugin discovery, lazy loading, validation, and isolation.
// Inputs: Mock libraries and writable runtime directory.
// Outputs: Zero on success, non-zero when an expectation fails.
int main(int argc, char** argv) {
  if(argc != 5) {
    std::cerr << "usage: GPluginManagerTest valid wrong-api missing-symbols runtime\n";
    return 2;
  }

  const std::filesystem::path validLibrary = argv[1];
  const std::filesystem::path wrongApiLibrary = argv[2];
  const std::filesystem::path missingSymbolsLibrary = argv[3];
  const std::filesystem::path runtime = argv[4];
  std::filesystem::remove_all(runtime);
  std::filesystem::create_directories(runtime);
  setenv("GROOT_PLUGIN_PATH", runtime.c_str(), 1);
  setenv("HOME", (runtime / "home").c_str(), 1);
  setenv("GSYS", (runtime / "gsys").c_str(), 1);

  auto& manager = GPluginManager::Get();
  std::string lastStatus;
  manager.SetStatusCallback([&lastStatus](const std::string& status) {
    lastStatus = status;
  });
  manager.Initialize();
  Check(manager.Actions().empty(), "empty search path should have no actions");

  MockSession session;
  GPluginContext sessionContext;
  sessionContext.canvas = reinterpret_cast<TCanvas*>(0x1000);
  sessionContext.pad = reinterpret_cast<TVirtualPad*>(0x2000);
  sessionContext.target = reinterpret_cast<TObject*>(0x3000);
  Check(manager.ActivateSession(&session, sessionContext),
        "first pad session should activate");
  Check(manager.HasActiveSession(sessionContext.pad),
        "active pad should be registered");
  MockSession duplicateSession;
  Check(!manager.ActivateSession(&duplicateSession, sessionContext),
        "a pad must reject a second session");
  GPluginContext duplicateTargetContext = sessionContext;
  duplicateTargetContext.pad = reinterpret_cast<TVirtualPad*>(0x4000);
  Check(!manager.ActivateSession(&duplicateSession, duplicateTargetContext),
        "a target must reject a second session in another pad");
  MockSession independentSession;
  GPluginContext independentContext = duplicateTargetContext;
  independentContext.target = reinterpret_cast<TObject*>(0x6000);
  Check(manager.ActivateSession(&independentSession, independentContext),
        "a different pad and target should allow another session");
  GPluginEvent sessionEvent;
  sessionEvent.context = sessionContext;
  sessionEvent.context.target = reinterpret_cast<TObject*>(0x5000);
  Check(manager.NotifySession(sessionEvent) && session.events == 1,
        "active session should observe its event");
  Check(session.lastTarget == sessionContext.target,
        "session notification should retain the locked target");
  manager.CloseSessionsForPad(sessionContext.pad);
  Check(session.closed && !manager.HasActiveSession(sessionContext.pad),
        "closing a pad should close and unregister its session");
  manager.CloseSessionsForPad(independentContext.pad);
  Check(independentSession.closed &&
        !manager.HasActiveSession(independentContext.pad),
        "independent pad session should close cleanly");

  const auto marker = runtime / "mock-created.txt";
  setenv("GROOT_MOCK_LOAD_MARKER", marker.c_str(), 1);
  WriteManifest(runtime / "10-valid.plugin", "mock", validLibrary.string(), "mock.ping");
  manager.Rescan();
  Check(manager.Actions().size() == 1, "valid manifest should register one action");
  Check(!std::filesystem::exists(marker), "plugin should remain unloaded after discovery");
  Check(manager.ExecuteAction("mock.ping"), "valid mock action should execute");
  Check(std::filesystem::exists(marker), "first execution should create plugin instance");

  WriteManifest(runtime / "20-duplicate-plugin.plugin", "mock",
                validLibrary.string(), "duplicate.plugin");
  WriteManifest(runtime / "30-duplicate-action.plugin", "duplicate-action",
                validLibrary.string(), "mock.ping");
  WriteManifest(runtime / "40-missing-library.plugin", "missing-library",
                (runtime / "does-not-exist.so").string(), "missing.library");
  WriteManifest(runtime / "50-wrong-api.plugin", "wrong-api",
                wrongApiLibrary.string(), "wrong.api",
                static_cast<int>(kGPluginApiVersion));
  WriteManifest(runtime / "60-missing-symbols.plugin", "missing-symbols",
                missingSymbolsLibrary.string(), "missing.symbols");
  WriteManifest(runtime / "65-manifest-api.plugin", "manifest-api",
                validLibrary.string(), "manifest.api", 999);
  {
    std::ofstream malformed(runtime / "70-malformed.plugin");
    malformed << "Plugin.Id: malformed\n";
  }

  manager.Rescan();
  Check(manager.Actions().size() == 5,
        "duplicates and malformed manifests should be rejected");
  Check(!manager.ExecuteAction("missing.library"),
        "missing library should fail without terminating");
  Check(!manager.ExecuteAction("wrong.api"),
        "runtime API mismatch should fail without terminating");
  Check(!manager.ExecuteAction("missing.symbols"),
        "missing entry symbols should fail without terminating");
  Check(!manager.ExecuteAction("manifest.api"),
        "manifest API mismatch should fail before loading");
  Check(!manager.ExecuteAction("unknown.action"),
        "unknown action should fail without terminating");
  Check(!lastStatus.empty(), "recoverable failures should reach the status callback");

  std::filesystem::remove_all(runtime);
  return failures == 0 ? 0 : 1;
}
