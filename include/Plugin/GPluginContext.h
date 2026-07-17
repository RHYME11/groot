#ifndef __GPLUGINCONTEXT_H__
#define __GPLUGINCONTEXT_H__

class TCanvas;
class TObject;
class TVirtualPad;

struct GPluginContext {
  TCanvas* canvas = nullptr;
  TVirtualPad* pad = nullptr;
  TObject* selected = nullptr;

  bool IsValid() const {
    return canvas && pad;
  }
};

#endif
