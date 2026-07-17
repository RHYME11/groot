#ifndef __GPLUGINLOADER_H__
#define __GPLUGINLOADER_H__

#include <string>

#include <Plugin/GPlugin.h>

class GPluginHost;

struct GLoadedPlugin {
  GPlugin* instance = nullptr;
  GPluginDestroyFunction destroy = nullptr;
  bool loaded = false;
};

class GPluginLoader {
  public:
    static bool Load(const std::string& library,
                     GPluginHost* host,
                     GLoadedPlugin& plugin,
                     std::string& error);
};

#endif
