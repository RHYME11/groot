#ifndef __GPLUGINHOST_H__
#define __GPLUGINHOST_H__

struct GPluginContext;
class GPluginSession;

class GPluginHost {
  public:
    virtual ~GPluginHost() = default;
    virtual void SetStatusMessage(const char* message) = 0;
    virtual bool ActivateSession(GPluginSession* session,
                                 const GPluginContext& context) = 0;
    virtual void DeactivateSession(GPluginSession* session) = 0;
};

#endif
