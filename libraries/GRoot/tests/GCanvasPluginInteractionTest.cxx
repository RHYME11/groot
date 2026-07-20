#include <cstdio>

#include <GCanvas.h>
#include <Plugin/GPluginContext.h>
#include <Plugin/GPluginEvent.h>
#include <Plugin/GPluginManager.h>
#include <Plugin/GPluginSession.h>

#include <Buttons.h>
#include <TBox.h>
#include <TROOT.h>

namespace {

class ObserverSession : public GPluginSession {
  public:
    const char* SessionId() const override {
      return "groot.test.observer";
    }

    void ObserveEvent(const GPluginEvent& event) override {
      ++events;
      target = event.context.target;
    }

    void Close() override {
      closed = true;
    }

    int events = 0;
    bool closed = false;
    TObject* target = nullptr;
};

// ============== Fail ==============
// Purpose: Report one interaction-test failure.
// Inputs: Failure description.
// Outputs: Non-zero process result.
int Fail(const char* message) {
  std::fprintf(stderr, "GCanvasPluginInteractionTest: %s\n", message);
  return 1;
}

} // namespace

// ============== main ==============
// Purpose: Verify ROOT native handling precedes plugin observation.
// Inputs: None.
// Outputs: Process success or failure.
int main() {
  gROOT->SetBatch(true);
  GCanvas canvas(false);
  canvas.SetName("plugin_native_test");
  canvas.SetCanvasSize(600, 400);
  canvas.Range(0.0, 0.0, 1.0, 1.0);
  canvas.cd();
  TBox probe(0.2, 0.2, 0.8, 0.8);
  probe.Draw();
  canvas.Update();
  canvas.SetSelected(&probe);
  canvas.SetSelectedPad(&canvas);

  ObserverSession session;
  GPluginContext context;
  context.canvas = &canvas;
  context.pad = &canvas;
  context.target = &probe;
  auto& manager = GPluginManager::Get();
  if(!manager.ActivateSession(&session, context))
    return Fail("session activation failed");

  const int px = canvas.XtoAbsPixel(0.5);
  const int py = canvas.YtoAbsPixel(0.5);
  canvas.HandleInput(kButton1Down, px, py);
  if(session.events != 1)
    return Fail("plugin session did not observe exactly one event");
  if(session.target != &probe)
    return Fail("plugin session lost its locked target");

  manager.CloseSessionsForCanvas(&canvas);
  if(!session.closed || manager.HasActiveSession(&canvas))
    return Fail("canvas close path did not release the session");
  manager.Shutdown();
  manager.CloseSessionsForCanvas(&canvas);
  return 0;
}
