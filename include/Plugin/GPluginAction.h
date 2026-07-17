#ifndef __GPLUGINACTION_H__
#define __GPLUGINACTION_H__

#include <string>

struct GPluginAction {
  std::string pluginId;
  std::string id;
  std::string label;
  std::string tooltip;
};

#endif
