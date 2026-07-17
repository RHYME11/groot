#ifndef __GPLUGINSESSIONREGISTRY_H__
#define __GPLUGINSESSIONREGISTRY_H__

#include <map>
#include <vector>

class GPlugin;
class GPluginSession;
class TCanvas;
class TVirtualPad;

struct GPluginSessionRecord {
  TCanvas* canvas = nullptr;
  TVirtualPad* pad = nullptr;
  GPlugin* plugin = nullptr;
  GPluginSession* session = nullptr;
};

class GPluginSessionRegistry {
  public:
    bool Add(const GPluginSessionRecord& record);
    GPluginSessionRecord* Find(TVirtualPad* pad);
    bool Contains(TVirtualPad* pad) const;
    void Remove(GPluginSession* session);
    std::vector<GPluginSessionRecord> TakePad(TVirtualPad* pad);
    std::vector<GPluginSessionRecord> TakeCanvas(TCanvas* canvas);
    std::vector<GPluginSessionRecord> TakeAll();

  private:
    std::map<TVirtualPad*, GPluginSessionRecord> fSessions;
};

#endif
