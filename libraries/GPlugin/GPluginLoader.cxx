#include "GPluginLoader.h"

#include <cstdint>

#include <TSystem.h>

#include <Plugin/GPluginHost.h>
#include <Plugin/GPluginVersion.h>

// ============== GPluginLoader::Load ==============
// Purpose: Load a plugin library and resolve its required C ABI entry points.
// Inputs: Library path, host interface, plugin state, and error output.
// Outputs: True when a compatible plugin instance was created.
bool GPluginLoader::Load(const std::string& library,
                         GPluginHost* host,
                         GLoadedPlugin& plugin,
                         std::string& error) {
  if(plugin.loaded && plugin.instance)
    return true;

  const int loadStatus = gSystem->Load(library.c_str());
  if(loadStatus < 0) {
    error = "ROOT failed to load library '" + library + "'";
    return false;
  }

  auto version = reinterpret_cast<GPluginApiVersionFunction>(
    gSystem->DynFindSymbol(library.c_str(), "GrootPluginApiVersion"));
  auto create = reinterpret_cast<GPluginCreateFunction>(
    gSystem->DynFindSymbol(library.c_str(), "GrootCreatePlugin"));
  auto destroy = reinterpret_cast<GPluginDestroyFunction>(
    gSystem->DynFindSymbol(library.c_str(), "GrootDestroyPlugin"));

  if(!version || !create || !destroy) {
    error = "library '" + library +
      "' is missing GrootPluginApiVersion, GrootCreatePlugin, or GrootDestroyPlugin";
    return false;
  }
  if(version() != kGPluginApiVersion) {
    error = "library '" + library + "' reports plugin API " +
      std::to_string(version()) + ", expected " +
      std::to_string(kGPluginApiVersion);
    return false;
  }

  plugin.instance = create(host);
  if(!plugin.instance) {
    error = "library '" + library + "' returned a null plugin instance";
    return false;
  }
  plugin.destroy = destroy;
  plugin.loaded = true;
  return true;
}
