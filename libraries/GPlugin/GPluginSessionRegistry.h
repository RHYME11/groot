#ifndef __GPLUGINSESSIONREGISTRY_H__
#define __GPLUGINSESSIONREGISTRY_H__

#include <map>
#include <vector>

class GPlugin;
class GPluginSession;
class TCanvas;
class TObject;
class TVirtualPad;

struct GPluginSessionRecord {
  TCanvas* canvas = nullptr;
  TVirtualPad* pad = nullptr;
  TObject* target = nullptr;
  GPlugin* plugin = nullptr;
  GPluginSession* session = nullptr;
};

class GPluginSessionRegistry {
  public:
    bool Add(const GPluginSessionRecord& record);
    GPluginSessionRecord* Find(TVirtualPad* pad);
    GPluginSession* FindTarget(TObject* target) const;
    bool Contains(TVirtualPad* pad) const;
    bool ContainsTarget(TObject* target) const;
    void Remove(GPluginSession* session);
    std::vector<GPluginSessionRecord> TakePad(TVirtualPad* pad);
    std::vector<GPluginSessionRecord> TakeCanvas(TCanvas* canvas);
    std::vector<GPluginSessionRecord> TakeAll();

  private:
    std::map<TVirtualPad*, GPluginSessionRecord> fSessions;
    std::map<TObject*, GPluginSession*> fTargets;
};

#endif
