#ifndef __GPLUGIN_H__
#define __GPLUGIN_H__

#include <cstdint>

#include <Plugin/GPluginContext.h>

class GPlugin;
class GPluginHost;

class GPlugin {
  public:
    virtual ~GPlugin() = default;
    virtual bool ExecuteAction(const char* actionId,
                               const GPluginContext& context) = 0;
    virtual bool CleanArtifacts(const GPluginContext&) {
      return true;
    }
    virtual const char* LastError() const = 0;
};

using GPluginApiVersionFunction = std::uint32_t (*)();
using GPluginCreateFunction = GPlugin* (*)(GPluginHost*);
using GPluginDestroyFunction = void (*)(GPlugin*);

#endif
