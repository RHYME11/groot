#ifndef __GPLUGINSESSION_H__
#define __GPLUGINSESSION_H__

class GPluginSession {
  public:
    virtual ~GPluginSession() = default;
    virtual void Close() = 0;
};

#endif
