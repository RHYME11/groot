
#include <GEventTimer.h>

#include <TList.h>
#include <TMethod.h>
#include <TClass.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TVirtualPad.h>


#include <GObjectManager.h>
#include <GGlobals.h>
#include <Histomatic.h>

TList *GetSignals(const char *classname) { 
  TMethod *method     = 0;
  const char *comment = 0;
  TList *siglist = new TList();
  siglist->Clear();
  TList *methodsList = TClass::GetClass(classname)->GetListOfMethods();
  TIter nextMethod(methodsList);
  while((method = (TMethod*)nextMethod())) {
    comment = method->GetCommentString();
    if(comment && strlen(comment) && strstr(comment,"*SIGNAL"))
      siglist->Add(method);
  }
  return siglist;
}

bool GEventTimer::Notify() {
  static unsigned long lastUpdateCount = 0;
  //printf("\n gHistomatic = 0x%08x \n",gHistomatic); fflush(stdout); 
  if(gHistomatic) { 
    gHistomatic->UpdateInfoPanel();
    unsigned long updateCount = GObjectManager::UpdateCount();
    if(updateCount != lastUpdateCount) {
      TIter next(gROOT->GetListOfCanvases());
      TCanvas *canvas = nullptr;
      while((canvas = static_cast<TCanvas*>(next()))) {
        canvas->Modified();
        canvas->Update();
      }
      lastUpdateCount = updateCount;
    }
    //gHistomatic->SetStatusText("gPad",0);
    //if(gPad) 
    //  gHistomatic->SetStatusText(gPad->GetTitle(),1);  

  }
  return TTimer::Notify(); 
}
