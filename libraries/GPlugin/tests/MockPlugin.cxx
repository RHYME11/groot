#include <cstdlib>
#include <fstream>
#include <string>

#include <Plugin/GPlugin.h>
#include <Plugin/GPluginHost.h>
#include <Plugin/GPluginVersion.h>

namespace {

class GPluginMock : public GPlugin {
  public:
    explicit GPluginMock(GPluginHost* host) : fHost(host) {
      if(const char* marker = std::getenv("GROOT_MOCK_LOAD_MARKER")) {
        std::ofstream output(marker);
        output << "created\n";
      }
    }

    bool ExecuteAction(const char* actionId,
                       const GPluginContext& context) override {
      if(!actionId || std::string(actionId) != "mock.ping") {
        fError = "unexpected action id";
        return false;
      }
      if(fHost)
        fHost->SetStatusMessage(context.IsValid() ? "mock context" : "mock ping");
      return true;
    }

    const char* LastError() const override {
      return fError.c_str();
    }

  private:
    GPluginHost* fHost;
    std::string fError;
};

} // namespace

// ============== GrootPluginApiVersion ==============
// Purpose: Report the mock plugin API version.
// Inputs: None.
// Outputs: Compatible or intentionally incompatible API version.
extern "C" std::uint32_t GrootPluginApiVersion() {
#ifdef MOCK_API_VERSION
  return MOCK_API_VERSION;
#else
  return kGPluginApiVersion;
#endif
}

// ============== GrootCreatePlugin ==============
// Purpose: Create the mock plugin instance.
// Inputs: Groot host interface.
// Outputs: New mock plugin owned by the caller.
extern "C" GPlugin* GrootCreatePlugin(GPluginHost* host) {
  return new GPluginMock(host);
}

// ============== GrootDestroyPlugin ==============
// Purpose: Destroy one mock plugin instance.
// Inputs: Plugin created by GrootCreatePlugin.
// Outputs: Released plugin instance.
extern "C" void GrootDestroyPlugin(GPlugin* plugin) {
  delete plugin;
}
