#include "GPluginRegistry.h"

// ============== GPluginRegistry::~GPluginRegistry ==============
// Purpose: Destroy all plugin instances while their libraries remain loaded.
// Inputs: None.
// Outputs: None.
GPluginRegistry::~GPluginRegistry() {
  Shutdown();
}

// ============== GPluginRegistry::Shutdown ==============
// Purpose: Destroy loaded plugin instances exactly once.
// Inputs: None.
// Outputs: Released plugin instances; libraries remain process-loaded.
void GPluginRegistry::Shutdown() {
  for(auto& entry : fPlugins) {
    auto& plugin = entry.second.loaded;
    if(plugin.instance && plugin.destroy) {
      plugin.destroy(plugin.instance);
      plugin.instance = nullptr;
    }
  }
}

// ============== GPluginRegistry::Add ==============
// Purpose: Register one validated plugin manifest and its actions.
// Inputs: Manifest metadata and error output.
// Outputs: True when the plugin and every action were registered.
bool GPluginRegistry::Add(const GPluginManifestData& manifest,
                          std::string& error) {
  const auto existing = fPlugins.find(manifest.id);
  if(existing != fPlugins.end() && existing->second.loaded.loaded &&
     existing->second.manifest.path == manifest.path) {
    return true;
  }
  if(existing != fPlugins.end()) {
    error = "duplicate plugin id '" + manifest.id + "'";
    return false;
  }
  for(const auto& action : manifest.actions) {
    if(fActionOwners.count(action.id)) {
      error = "duplicate action id '" + action.id + "'";
      return false;
    }
  }

  GPluginRecord record;
  record.manifest = manifest;
  fPlugins.emplace(manifest.id, record);
  for(const auto& action : manifest.actions) {
    fActionOwners.emplace(action.id, manifest.id);
    fActions.push_back(action);
  }
  return true;
}

// ============== GPluginRegistry::FindByAction ==============
// Purpose: Find the owning plugin record for an action identifier.
// Inputs: Global action identifier.
// Outputs: Plugin record pointer, or null when the action is unknown.
GPluginRecord* GPluginRegistry::FindByAction(const std::string& actionId) {
  const auto owner = fActionOwners.find(actionId);
  if(owner == fActionOwners.end())
    return nullptr;
  const auto plugin = fPlugins.find(owner->second);
  return plugin == fPlugins.end() ? nullptr : &plugin->second;
}

// ============== GPluginRegistry::LoadedPlugins ==============
// Purpose: Return records whose shared libraries created plugin instances.
// Inputs: None.
// Outputs: Mutable loaded-plugin records in registry order.
std::vector<GPluginRecord*> GPluginRegistry::LoadedPlugins() {
  std::vector<GPluginRecord*> result;
  for(auto& entry : fPlugins) {
    if(entry.second.loaded.instance)
      result.push_back(&entry.second);
  }
  return result;
}

// ============== GPluginRegistry::Actions ==============
// Purpose: Return actions in manifest discovery order.
// Inputs: None.
// Outputs: Read-only action list.
const std::vector<GPluginAction>& GPluginRegistry::Actions() const {
  return fActions;
}

// ============== GPluginRegistry::ClearMetadata ==============
// Purpose: Clear unloaded metadata before a rescan.
// Inputs: None.
// Outputs: None; loaded plugins are retained for v1 lifetime safety.
void GPluginRegistry::ClearMetadata() {
  std::map<std::string, GPluginRecord> loadedPlugins;
  for(auto& entry : fPlugins) {
    if(entry.second.loaded.loaded)
      loadedPlugins.emplace(entry.first, std::move(entry.second));
  }
  fPlugins = std::move(loadedPlugins);
  fActionOwners.clear();
  fActions.clear();

  for(const auto& entry : fPlugins) {
    for(const auto& action : entry.second.manifest.actions) {
      fActionOwners.emplace(action.id, entry.first);
      fActions.push_back(action);
    }
  }
}
