#ifndef __GPLUGINMANAGER_H__
#define __GPLUGINMANAGER_H__

#include <functional>
#include <string>
#include <vector>

#include <Plugin/GPluginAction.h>
#include <Plugin/GPluginContext.h>
#include <Plugin/GPluginHost.h>

class GPluginRegistry;
class GPluginSession;
class GPluginSessionRegistry;
struct GPluginEvent;
class TCanvas;
class TVirtualPad;

class GPluginManager : public GPluginHost {
  public:
    static GPluginManager& Get();
    ~GPluginManager() override;

    void Initialize();
    void Shutdown();
    void Rescan();

    const std::vector<GPluginAction>& Actions() const;
    bool ExecuteAction(const std::string& actionId);
    void CleanArtifacts(const GPluginContext& context);

    void SetActionChangedCallback(std::function<void()> callback);
    void SetContextProvider(std::function<GPluginContext()> provider);
    void SetStatusCallback(std::function<void(const std::string&)> callback);
    void SetStatusMessage(const char* message) override;
    bool ActivateSession(GPluginSession* session,
                         const GPluginContext& context) override;
    void DeactivateSession(GPluginSession* session) override;

    bool HasActiveSession(TVirtualPad* pad) const;
    bool NotifySession(const GPluginEvent& event);
    void CloseSessionsForPad(TVirtualPad* pad);
    void CloseSessionsForCanvas(TCanvas* canvas);

  private:
    GPluginManager();
    GPluginManager(const GPluginManager&) = delete;
    GPluginManager& operator=(const GPluginManager&) = delete;

    std::vector<std::string> SearchPaths() const;
    void ScanDirectory(const std::string& path);
    void ReportError(const std::string& message);

    GPluginRegistry* fRegistry;
    GPluginSessionRegistry* fSessions;
    bool fInitialized;
    bool fShutdown;
    std::function<void()> fActionChangedCallback;
    std::function<GPluginContext()> fContextProvider;
    std::function<void(const std::string&)> fStatusCallback;
};

#endif
