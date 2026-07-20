#include "GPluginSessionRegistry.h"

// ============== GPluginSessionRegistry::Add ==============
// Purpose: Register one application-overlay session for a pad and target.
// Inputs: Session ownership record.
// Outputs: True when the pad was previously inactive.
bool GPluginSessionRegistry::Add(const GPluginSessionRecord& record) {
  if(!record.pad || !record.session || Contains(record.pad) ||
     ContainsTarget(record.target))
    return false;
  fSessions.emplace(record.pad, record);
  if(record.target)
    fTargets.emplace(record.target, record.session);
  return true;
}

// ============== GPluginSessionRegistry::Find ==============
// Purpose: Find the active session for one pad.
// Inputs: Pad identity.
// Outputs: Mutable record or null.
GPluginSessionRecord* GPluginSessionRegistry::Find(TVirtualPad* pad) {
  const auto found = fSessions.find(pad);
  return found == fSessions.end() ? nullptr : &found->second;
}

// ============== GPluginSessionRegistry::FindTarget ==============
// Purpose: Find the session locking one target object.
// Inputs: Target object identity.
// Outputs: Session pointer or null.
GPluginSession* GPluginSessionRegistry::FindTarget(TObject* target) const {
  const auto found = fTargets.find(target);
  return found == fTargets.end() ? nullptr : found->second;
}

// ============== GPluginSessionRegistry::Contains ==============
// Purpose: Test whether one pad has an application-overlay session.
// Inputs: Pad identity.
// Outputs: True when registered.
bool GPluginSessionRegistry::Contains(TVirtualPad* pad) const {
  return pad && fSessions.count(pad) != 0;
}

// ============== GPluginSessionRegistry::ContainsTarget ==============
// Purpose: Test whether an object is locked by an interactive session.
// Inputs: Target object identity.
// Outputs: True when the non-null target is already registered.
bool GPluginSessionRegistry::ContainsTarget(TObject* target) const {
  return target && fTargets.count(target) != 0;
}

// ============== GPluginSessionRegistry::Remove ==============
// Purpose: Unregister one session without closing it.
// Inputs: Session identity.
// Outputs: Updated registry.
void GPluginSessionRegistry::Remove(GPluginSession* session) {
  for(auto it = fSessions.begin(); it != fSessions.end(); ++it) {
    if(it->second.session == session) {
      if(it->second.target)
        fTargets.erase(it->second.target);
      fSessions.erase(it);
      return;
    }
  }
}

// ============== GPluginSessionRegistry::TakePad ==============
// Purpose: Detach the session owned by one pad.
// Inputs: Pad identity.
// Outputs: Detached records.
std::vector<GPluginSessionRecord> GPluginSessionRegistry::TakePad(TVirtualPad* pad) {
  std::vector<GPluginSessionRecord> result;
  const auto found = fSessions.find(pad);
  if(found != fSessions.end()) {
    result.push_back(found->second);
    if(found->second.target)
      fTargets.erase(found->second.target);
    fSessions.erase(found);
  }
  return result;
}

// ============== GPluginSessionRegistry::TakeCanvas ==============
// Purpose: Detach every session below one canvas.
// Inputs: Canvas identity.
// Outputs: Detached records.
std::vector<GPluginSessionRecord> GPluginSessionRegistry::TakeCanvas(TCanvas* canvas) {
  std::vector<GPluginSessionRecord> result;
  for(auto it = fSessions.begin(); it != fSessions.end();) {
    if(it->second.canvas == canvas) {
      result.push_back(it->second);
      if(it->second.target)
        fTargets.erase(it->second.target);
      it = fSessions.erase(it);
    } else {
      ++it;
    }
  }
  return result;
}

// ============== GPluginSessionRegistry::TakeAll ==============
// Purpose: Detach every active session.
// Inputs: None.
// Outputs: Detached records.
std::vector<GPluginSessionRecord> GPluginSessionRegistry::TakeAll() {
  std::vector<GPluginSessionRecord> result;
  for(const auto& entry : fSessions)
    result.push_back(entry.second);
  fSessions.clear();
  fTargets.clear();
  return result;
}
