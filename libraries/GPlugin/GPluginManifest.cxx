#include "GPluginManifest.h"

#include <filesystem>
#include <limits>

#include <TEnv.h>

namespace {

// ============== ReadRequired ==============
// Purpose: Read one required non-empty manifest string.
// Inputs: Environment, key, output value, and error output.
// Outputs: True when the field is present and non-empty.
bool ReadRequired(TEnv& env, const char* key, std::string& value,
                  std::string& error) {
  constexpr const char* missing = "__GROOT_PLUGIN_MISSING__";
  value = env.GetValue(key, missing);
  if(value.empty() || value == missing) {
    error = std::string("missing required field '") + key + "'";
    return false;
  }
  return true;
}

} // namespace

// ============== GPluginManifest::Read ==============
// Purpose: Parse and validate one ROOT TEnv plugin manifest.
// Inputs: Manifest path and output references.
// Outputs: True with populated metadata, or false with an error message.
bool GPluginManifest::Read(const std::string& path,
                           GPluginManifestData& manifest,
                           std::string& error) {
  TEnv env;
  if(env.ReadFile(path.c_str(), kEnvLocal) < 0) {
    error = "unable to read manifest";
    return false;
  }

  manifest = GPluginManifestData{};
  manifest.path = path;
  manifest.directory = std::filesystem::path(path).parent_path().string();

  if(!ReadRequired(env, "Plugin.Id", manifest.id, error) ||
     !ReadRequired(env, "Plugin.Name", manifest.name, error) ||
     !ReadRequired(env, "Plugin.Version", manifest.version, error) ||
     !ReadRequired(env, "Plugin.Library", manifest.library, error)) {
    return false;
  }

  const int invalid = std::numeric_limits<int>::min();
  const int apiVersion = env.GetValue("Plugin.ApiVersion", invalid);
  const int actionCount = env.GetValue("Plugin.Action.Count", invalid);
  if(apiVersion < 0 || apiVersion == invalid) {
    error = "missing or invalid field 'Plugin.ApiVersion'";
    return false;
  }
  if(actionCount < 1 || actionCount == invalid) {
    error = "missing or invalid field 'Plugin.Action.Count'";
    return false;
  }
  manifest.apiVersion = static_cast<std::uint32_t>(apiVersion);

  for(int index = 0; index < actionCount; ++index) {
    GPluginAction action;
    action.pluginId = manifest.id;
    const std::string prefix = "Plugin.Action." + std::to_string(index) + ".";
    if(!ReadRequired(env, (prefix + "Id").c_str(), action.id, error) ||
       !ReadRequired(env, (prefix + "Label").c_str(), action.label, error)) {
      return false;
    }
    action.tooltip = env.GetValue((prefix + "Tooltip").c_str(), "");
    manifest.actions.push_back(action);
  }

  return true;
}
