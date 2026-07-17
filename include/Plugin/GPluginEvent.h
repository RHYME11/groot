#ifndef __GPLUGINEVENT_H__
#define __GPLUGINEVENT_H__

struct GPluginEvent {
  int type = 0;
  int code = 0;
  int state = 0;
  int px = 0;
  int py = 0;
};

#endif
