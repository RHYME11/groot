#ifndef __GPLUGINSESSION_H__
#define __GPLUGINSESSION_H__

#include <Plugin/GPluginEvent.h>

class GPluginSession {
  public:
    virtual ~GPluginSession() = default;
    virtual const char* SessionId() const = 0;
    virtual void ObserveEvent(const GPluginEvent& event) = 0;
    virtual void Close() = 0;
};

#endif
