#include "GPluginSessionRegistry.h"

// ============== GPluginSessionRegistry::Add ==============
// Purpose: Register one exclusive pad session.
// Inputs: Session ownership record.
// Outputs: True when the pad was previously inactive.
bool GPluginSessionRegistry::Add(const GPluginSessionRecord& record) {
  if(!record.pad || !record.session || Contains(record.pad))
    return false;
  fSessions.emplace(record.pad, record);
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

// ============== GPluginSessionRegistry::Contains ==============
// Purpose: Test whether one pad has an exclusive session.
// Inputs: Pad identity.
// Outputs: True when registered.
bool GPluginSessionRegistry::Contains(TVirtualPad* pad) const {
  return pad && fSessions.count(pad) != 0;
}

// ============== GPluginSessionRegistry::Remove ==============
// Purpose: Unregister one session without closing it.
// Inputs: Session identity.
// Outputs: Updated registry.
void GPluginSessionRegistry::Remove(GPluginSession* session) {
  for(auto it = fSessions.begin(); it != fSessions.end(); ++it) {
    if(it->second.session == session) {
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
  return result;
}
