#ifndef __GPLUGINREGISTRY_H__
#define __GPLUGINREGISTRY_H__

#include <map>
#include <string>
#include <vector>

#include <Plugin/GPluginAction.h>

#include "GPluginLoader.h"
#include "GPluginManifest.h"

struct GPluginRecord {
  GPluginManifestData manifest;
  GLoadedPlugin loaded;
};

class GPluginRegistry {
  public:
    ~GPluginRegistry();

    bool Add(const GPluginManifestData& manifest, std::string& error);
    GPluginRecord* FindByAction(const std::string& actionId);
    const std::vector<GPluginAction>& Actions() const;
    void ClearMetadata();
    void Shutdown();

  private:
    std::map<std::string, GPluginRecord> fPlugins;
    std::map<std::string, std::string> fActionOwners;
    std::vector<GPluginAction> fActions;
};

#endif
