#ifndef __GPLUGINMANIFEST_H__
#define __GPLUGINMANIFEST_H__

#include <cstdint>
#include <string>
#include <vector>

#include <Plugin/GPluginAction.h>

struct GPluginManifestData {
  std::string path;
  std::string directory;
  std::string id;
  std::string name;
  std::string version;
  std::uint32_t apiVersion = 0;
  std::string library;
  std::vector<GPluginAction> actions;
};

class GPluginManifest {
  public:
    static bool Read(const std::string& path,
                     GPluginManifestData& manifest,
                     std::string& error);
};

#endif
