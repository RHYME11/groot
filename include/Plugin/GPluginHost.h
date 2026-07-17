#ifndef __GPLUGINHOST_H__
#define __GPLUGINHOST_H__

class GPluginHost {
  public:
    virtual ~GPluginHost() = default;
    virtual void SetStatusMessage(const char* message) = 0;
};

#endif
