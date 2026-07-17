#ifndef __GPLUGINEVENT_H__
#define __GPLUGINEVENT_H__

#include <Plugin/GPluginContext.h>

struct GPluginEvent {
  GPluginContext context;
  int type = 0;
  int code = 0;
  int state = 0;
  int px = 0;
  int py = 0;
  double x = 0.0;
  double y = 0.0;
};

#endif
