
#include <Histomatic.h>
//#include <utils.h>
#include <vector>

#include "TQObject.h"
#include <TROOT.h>
#include <TVirtualPad.h>
#include <TObjString.h>
#include <TH1.h>
#include <TH2.h>
#include <TFrame.h>
#include <THStack.h>
#include <TGraph.h>
#include <TCutG.h>
#include <TF1.h>
#include <TList.h>
#include <TStyle.h>
#include <TPad.h>
#include <TVirtualX.h>
#include <TGLabel.h>

//#include <GObjectManager.h>
#include <GObjectManager.h>
#include <GCanvas.h>
#include <GH1D.h>
#include <GH2D.h>

#include <GGlobals.h>
#include <GEventTimer.h>

#include <GCommands.h>
#include <Plugin/GPluginContext.h>
#include <Plugin/GPluginManager.h>

#include <GListTree.h>

extern Histomatic *gHistomatic;

////////////////////

Histomatic *Histomatic::fInstance = 0;

////////////////////
////////////////////

GInfoPanel::GInfoPanel(const TGWindow* parent) 
  : TGGroupFrame(parent,"i am a title!") { //,
  //fObject(0),fPosition(0),fBin(0),fCounts(0),fMarker(0),fMode(0) {

    //fObject   = new TGLabel(this, "Object:");
    //fPosition = new TGLabel(this, "Cursor:");
    //fBin      = new TGLabel(this, "Bin:");
    //fCounts   = new TGLabel(this, "Counts:");
    //fMarker   = new TGLabel(this, "Marker:");
    //fMode     = new TGLabel(this, "Mode:");

    //AddFrame(fObject,   new TGLayoutHints(kLHintsExpandX, 4,4,2,2));
    //AddFrame(fPosition, new TGLayoutHints(kLHintsExpandX, 4,4,2,2));
    //AddFrame(fBin,      new TGLayoutHints(kLHintsExpandX, 4,4,2,2));
    //AddFrame(fCounts,   new TGLayoutHints(kLHintsExpandX, 4,4,2,2));
    //AddFrame(fMarker,   new TGLayoutHints(kLHintsExpandX, 4,4,2,2));
    //AddFrame(fMode,     new TGLayoutHints(kLHintsExpandX, 4,4,2,2));

 AddRow("Object",  "");
  AddRow("Cursor",  "");
  AddRow("Bin",     "");
  AddRow("Counts",  "");
  AddRow("Marker",  "");
  AddRow("Mode",    "");


  }

GInfoPanel::~GInfoPanel() { }

void GInfoPanel::AddRow(const std::string &key, const std::string &value) {
  auto* row = new TGHorizontalFrame(this);

  auto* keyLabel = new TGLabel(row,Form("%s:",key.c_str()));
  auto* valLabel = new TGLabel(row,value.c_str());

  keyLabel->SetTextJustify(kTextLeft);
  valLabel->SetTextJustify(kTextRight);

  row->AddFrame(keyLabel, new TGLayoutHints(kLHintsLeft, 4, 10, 2, 2));
  row->AddFrame(valLabel, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

  AddFrame(row, new TGLayoutHints(kLHintsExpandX));

  fRows[key] = valLabel;
  return;
}


void GInfoPanel::SetRow(const std::string &key, const std::string &value) { 
  auto it = fRows.find(key);

  if(it == fRows.end()) {
    AddRow(key, value);
  } else {
    it->second->SetText(value.c_str());
  }

  Layout();
}


void GInfoPanel::Update(const GInteractionInfo &info) {

  if(!info.IsValid()) {
    SetRow("Cursor", "");
    SetRow("Counts", "");
    SetRow("Marker", "Primary");
    SetRow("Mode", "click=marker, Ctrl-click=ignore max");
    return;
  }
  SetRow("Object", info.targetName);
  SetRow("Cursor", Form("x = %.3f", info.x));
  SetRow("Counts", Form("%.3f", info.counts));
  SetRow("Marker", "Primary");
  SetRow("Mode", "click=marker, Ctrl-click=ignore max");

}




////////////////////
////////////////////

//Histomatic::Histomatic() : TGMainFrame(gClient->GetRoot(),200,600) {   
Histomatic::Histomatic() : TGMainFrame(gClient->GetRoot(),350,780), fVf(0) {

  int width  = 350;
  int height = 780;

  GPluginManager::Get().SetActionChangedCallback([this]() {
    RefreshPluginActions();
  });
  GPluginManager::Get().SetContextProvider([]() {
    GPluginContext context;
    context.pad = gPad;
    context.canvas = context.pad ? context.pad->GetCanvas() : nullptr;
    context.selected = context.canvas ? context.canvas->GetSelected() : nullptr;
    return context;
  });
  GPluginManager::Get().SetStatusCallback([this](const std::string& message) {
    SetStatusText(message, 0);
  });

  CreateWindow();
  this->SetWindowName("hist-o-matic");

  int dh = gClient->GetDisplayHeight();
  int dw = gClient->GetDisplayWidth();
  //printf("dh = %i\ndw = %i\n\n",dh,dw);

  Show(width,height);

  fEventTimer = new GEventTimer();
  fEventTimer->Start();

  //this->Move(dw-width*1.1,50);
}

Histomatic *Histomatic::Get() {
  if(!fInstance) fInstance = new Histomatic();
  gHistomatic = fInstance;
  return fInstance;
}

void Histomatic::Show(int width, int height) {
  //MoveResize(1,1,200,600);
  this->Resize(width,height);
  //SetWMPosition(1,1);
  //this->Resize(200,600);
  this->MapWindow();
}

Histomatic::~Histomatic() {
  GPluginManager::Get().SetActionChangedCallback({});
  GPluginManager::Get().SetContextProvider({});
  GPluginManager::Get().SetStatusCallback({});

  if(gHistomatic == this)
    gHistomatic = 0;

  if(fEventTimer) {
    fEventTimer->Stop();
    delete fEventTimer;
  }

  delete fLH0;
  delete fLH1;
  delete fLH2;

  //TGMenuBar *fMenuBar;
  //TGHorizontalFrame *fTopMenuFrame;               
  //TGHorizontalFrame *fPreMenuFrame; 

  //TGPopupMenu       *fMenuFile;
  //TGPopupMenu       *fMenuHelp;     

  for(auto* button : fPluginButtons)
    delete button;
  delete fButton4; 
  delete fButton8; 

  delete fDrawComboBox;
  delete fDrawColz;
  delete fDrawNormalized;

  delete fGListTree;
  delete fGListTreeCanvas;

  delete fInfoPanel;
  //delete fStatusBar;

  delete fButtonRow1;
  delete fButtonRow2;
  delete fButtonContainer;
  delete fDrawOptionContainer;

  delete fVf;
}

void Histomatic::CreateWindow() {

  if(fVf) {
    Show();
    return;
  }


  fVf = new TGVerticalFrame(this,100,100);

  fLH0 = new TGLayoutHints(kLHintsExpandX,0,0,0,0);
  fLH1 = new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0);
  fLH2 = new TGLayoutHints(kLHintsCenterX,5,5,3,4);

  //fTopMenuFrame = new TGHorizontalFrame(fVf,100,20);

  //fPreMenuFrame = new TGHorizontalFrame(fTopMenuFrame,0,20,kRaisedFrame);
  //fMenuBar  = new TGMenuBar(fPreMenuFrame,10,10,kHorizontalFrame);
  //fMenuFile = new TGPopupMenu(gClient->GetDefaultRoot());
  //fMenuFile->AddEntry("thing1",kThing1);
  //fMenuFile->AddEntry("thing2",kThing2);

  //fMenuHelp = new TGPopupMenu(fClient->GetRoot());
  //fMenuHelp->AddEntry("Send Help", kSendHelp);

  //fPreMenuFrame->AddFrame(fMenuBar, fLH1);
  //fTopMenuFrame->AddFrame(fPreMenuFrame, fLH0);

  fButtonContainer = new TGVerticalFrame(fVf,200,200);

  fButtonRow1 = new TGHorizontalFrame(fButtonContainer,10,10);
  fButton4 = new TGTextButton(fButtonRow1,"CloseAll");
  fButton4->Connect("Clicked()","Histomatic",this,"closeAllCanvases()");
  fButtonRow1->AddFrame(fButton4,fLH1);

  fButtonRow2 = new TGHorizontalFrame(fButtonContainer,10,10);
  fButton8 = new TGTextButton(fButtonRow2,"do Draw");
  fButton8->Connect("Clicked()","Histomatic",this,"doDraw()");

  fButtonRow2->AddFrame(fButton8,fLH1);

  fButtonContainer->AddFrame(fButtonRow1,fLH1);  
  fButtonContainer->AddFrame(fButtonRow2,fLH1);  
  RefreshPluginActions();

  fDrawOptionContainer = new TGHorizontalFrame(fVf,100,40);

  fDrawComboBox = new TGComboBox(fDrawOptionContainer,100);
  fDrawComboBox->AddEntry("new canvas ",EDrawOption::eDrawNew);
  fDrawComboBox->AddEntry("same canvas",EDrawOption::eDrawSame);
  fDrawComboBox->AddEntry("stacked    ",EDrawOption::eDrawStacked);
  fDrawComboBox->Select(0);

  fDrawNormalized = new TGCheckButton(fDrawOptionContainer,"normalized",1);
  fDrawNormalized->SetState(kButtonUp);

  fDrawColz = new TGCheckButton(fDrawOptionContainer,"colz",1);
  fDrawColz->SetState(kButtonDown);

  fLockPads = new TGCheckButton(fDrawOptionContainer,"lock pads",1);
  fLockPads->SetState(kButtonUp);
  fLockPads->Connect("Clicked()","Histomatic",this,"doLockPads(TPad*)");
  //Connect("doLockPads()","Histomatic",this,"doLockPads()");
  //Connect("TCanvas","Picked()",this,"doLockPads()");
  TQObject::Connect("TCanvas","Picked(TPad*,TObject*,Int_t)","Histomatic",this,"doLockPads(TPad*)");
  /*
     Connect (const char *sender_class, 
     const char *signal, 
     const char *receiver_class, 
     void *receiver, c
     onst char *slot) */



  fDrawOptionContainer->AddFrame(fDrawComboBox,fLH1);
  fDrawOptionContainer->AddFrame(fDrawNormalized,fLH1);
  fDrawOptionContainer->AddFrame(fDrawColz,fLH0);
  fDrawOptionContainer->AddFrame(fLockPads,fLH0);

  fGListTreeCanvas = new GListTreeCanvas(fVf,10,10);
  fGListTree = new GListTree(fGListTreeCanvas); 

  fInfoPanel = new GInfoPanel(fVf);

  fStatusBar = new TGStatusBar(fVf,100,50);
  fStatusBar->SetParts(4);
  {
    //fStatusBar->SetText("I AM A STATUS BAR");
    SetStatusText("some junk",0);
    fStatusBar->SetBackgroundColor(fStatusBar->GetBlackPixel());
    fStatusBar->SetForegroundColor(fStatusBar->GetWhitePixel());
  }

  fVf->AddFrame(fButtonContainer,fLH0);
  fVf->AddFrame(fDrawOptionContainer,fLH0);
  fVf->AddFrame(fGListTreeCanvas,fLH1);
  fVf->AddFrame(fInfoPanel,new TGLayoutHints(kLHintsExpandX,2,2,4,4));
  fVf->AddFrame(fStatusBar,fLH0);


  this->AddFrame(fVf,fLH1);
  this->MapSubwindows();
  this->Resize(this->GetDefaultSize());
}

void Histomatic::doLockPads(TPad *pad) {
  if(!gPad)
    return;
  if(!gPad->GetCanvas()->InheritsFrom(GCanvas::Class())) 
    return;

  //printf("pad = 0x%p\n",pad);
  //printf("gPad = 0x%p\n",gPad);


  if(pad!=gPad) { // in a new pad - set the button
    if(pad && pad->GetCanvas()->InheritsFrom(GCanvas::Class())) {
      fLockPads->SetState(((GCanvas*)pad->GetCanvas())->GetLockPads() ? kButtonDown : kButtonUp);
    }
  }

  //TODO -the feedback below is just broken.  need to fix
  return;


  //the pad and the gPad are always the same...(?) 
  if(((GCanvas*)gPad->GetCanvas())->GetLockPads()) {
    fLockPads->SetState(kButtonDown);
  } else {
    fLockPads->SetState(kButtonUp);
  }

  //if(fLockPads->GetState() == kButtonDown) {
  //  ((GCanvas*)gPad->GetCanvas())->SetLockPads(true);
  //} else {
  //  ((GCanvas*)gPad->GetCanvas())->SetLockPads(false);
  //}
}


void Histomatic::CloseWindow() {
  gVirtualX->UnmapWindow(fId);
  //fMainWindow->CloseWindow();  
  //  DeleteWindow();
}


// ============== Histomatic::RefreshPluginActions ==============
// Purpose: Rebuild plugin action buttons from the manager registry.
// Inputs: Current plugin action registry.
// Outputs: Updated dynamic button rows.
void Histomatic::RefreshPluginActions() {
  if(!fButtonRow1 || !fButtonRow2)
    return;

  for(auto* button : fPluginButtons) {
    if(button->GetParent() == fButtonRow1)
      fButtonRow1->RemoveFrame(button);
    else
      fButtonRow2->RemoveFrame(button);
    delete button;
  }
  fPluginButtons.clear();
  fPluginButtonActions.clear();

  constexpr int firstPluginButtonId = 1000;
  int index = 0;
  for(const auto& action : GPluginManager::Get().Actions()) {
    TGHorizontalFrame* row = index % 2 == 0 ? fButtonRow1 : fButtonRow2;
    const int buttonId = firstPluginButtonId + index;
    auto* button = new TGTextButton(row, action.label.c_str(), buttonId);
    button->Associate(this);
    if(!action.tooltip.empty())
      button->SetToolTipText(action.tooltip.c_str());
    row->AddFrame(button, fLH1);
    fPluginButtons.push_back(button);
    fPluginButtonActions.emplace(buttonId, action.id);
    ++index;
  }

  fButtonRow1->MapSubwindows();
  fButtonRow2->MapSubwindows();
  Layout();
}

// ============== Histomatic::ProcessMessage ==============
// Purpose: Dispatch dynamic plugin button messages to the plugin manager.
// Inputs: ROOT GUI message and button identifier.
// Outputs: True after the message is handled or delegated.
Bool_t Histomatic::ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2) {
  if(GET_MSG(msg) == kC_COMMAND && GET_SUBMSG(msg) == kCM_BUTTON) {
    const auto action = fPluginButtonActions.find(static_cast<int>(parm1));
    if(action != fPluginButtonActions.end()) {
      GPluginManager::Get().ExecuteAction(action->second);
      return kTRUE;
    }
  }
  return TGMainFrame::ProcessMessage(msg, parm1, parm2);
}

void Histomatic::doDraw() {
  doDraw(fGListTree->GetSelected(),"");
}

void Histomatic::doDraw(std::vector<TGListTreeItem*> selected, Option_t *opt) const {


  TList hList;

  std::vector<TH1*>    hists1D;
  std::vector<TH2*>    hists2D;
  std::vector<TGraph*> graphs;
  std::vector<TF1*>    functions;

  int drawables = 0;
  for(auto item=selected.begin();item!=selected.end();item++) {

    TKey          *key = fGListTree->GetKey(*item);  
    if(!key && !GObjectManager::HasSource(fGListTree->GetFileName(*item))) continue;
    std::string fullPath = fGListTree->GetFullPath(*item);

    TObject *obj = fGListTree->GetObject(*item);
    if(!obj) continue;

    if(obj->InheritsFrom(TH2::Class())) {
      hists2D.push_back(static_cast<TH2*>(obj));
      drawables++;
    } else if(obj->InheritsFrom(TH1::Class())) {
      hists1D.push_back(static_cast<TH1*>(obj));
      drawables++;
    } else if(obj->InheritsFrom(TGraph::Class())) {
      graphs.push_back(static_cast<TGraph*>(obj));
      drawables++;
    } else if(obj->InheritsFrom(TF1::Class())) {
      functions.push_back(static_cast<TF1*>(obj));
      drawables++;
    } else {
      //pass;
    }
  }

  //printf("\n\nnum drawables: %i\n",drawables);
  //printf("found  %lu TH1 to draw.\n",    hists1D.size());
  //printf("found  %lu TH2 to draw.\n",    hists2D.size());
  //printf("found  %lu TGraph to draw.\n", graphs.size());
  //printf("found  %lu TF1    to draw.\n", functions.size());
  //printf("\n\n");

  if(drawables<1) return;

  //TODO - do something with the opt....

  //if(hList.GetSize()<1) return;

  //printf("doDraw\n");
  //printf("draw option: %i\n",fDrawComboBox->GetSelected());
  //printf("normalized:  %i\n",fDrawNormalized->GetState());
  //printf("colz:        %i\n",fDrawColz->GetState());

  TCanvas *g = 0;

  if(hists1D.size()>0 || hists2D.size()>0) {
    switch (fDrawComboBox->GetSelected()) {
      case EDrawOption::eDrawNew:
        g = gROOT->MakeDefCanvas();
        break;
      case EDrawOption::eDrawSame:
        if(gPad)
          g = gPad->GetCanvas();
        else
          g = gROOT->MakeDefCanvas();  //new GCanvas;
      case EDrawOption::eDrawStacked:
        if(gPad)
          g = gPad->GetCanvas();
        else
          g = gROOT->MakeDefCanvas();  //new GCanvas;
        break;
    }
  }

  if(hists1D.size()>0) {

    for(auto it=hists1D.begin();it!=hists1D.end();it++) {
      if((*it)->InheritsFrom(GH1D::Class())) {
        GH1D *current = static_cast<GH1D*>(*it);
        if(fDrawNormalized->GetState() == kButtonDown) {
          if(!current->IsNormalized()) {
            current->Normalize();
          }
        } else {
          if(current->IsNormalized()) {
            current->Normalize(); // this will un-normalize the histogram.
          }
        }
      }
    }
    drawHists(hists1D,g);

    /*
       if(hists1D.size()==1) {
       hists1D.at(0)->Draw();
       } else {
       if(hists1D.size()<=5) {
       g->Divide(1,hists1D.size(),0.01,0);
       int padN=1;
       for(auto it=hists1D.begin();it!=hists1D.end();it++) {
       g->cd(padN++);
       (*it)->Draw();
       }
       } else {
       THStack hs;
       for(auto it=hists1D.begin();it!=hists1D.end();it++) 
       hs.Add(*it);

       if(fDrawComboBox->GetSelected() == EDrawOption::eDrawSame)
       hs.Draw();
    //hs.DrawNormalized();
    else
    hs.Draw("pads");
    //hs.DrawNormalized("pads");

    }
    }
     */
  }
  if(hists2D.size()==1) {
    GH2D *current = static_cast<GH2D*>(hists2D[0]);
    if(fDrawComboBox->GetSelected() == EDrawOption::eDrawNew) {
      if(!g)
        g = gROOT->MakeDefCanvas();//new GCanvas;
      g->cd();

      current->Draw();
    }else if(fDrawComboBox->GetSelected() == EDrawOption::eDrawSame) {
      g = g->GetCanvas(); // let's make sure we are not starting from a subpad.

      std::map<int,std::vector<TObject*> > fPadMap;

      int nPad = 0;
      int mPad = 0;
      TVirtualPad *cpad = g->cd(0);
      TVirtualPad *lastPad = cpad;
      do {
        TList *prim = cpad->GetListOfPrimitives();
        for(int i=0;i<prim->GetSize();i++) {
          TObject *obj = prim->At(i);
          if(obj  && (obj->InheritsFrom(TH1::Class()) || obj->InheritsFrom(TGraph::Class()))) {
            printf("found obj %s\n",obj->GetName());
            fPadMap[mPad].push_back(obj);
          }
        }
        nPad++;
        printf("fPadMap[%i].size() = %lu\n",mPad,fPadMap[mPad].size());
        if(mPad==0 || fPadMap[mPad].size()>0)
          mPad++;
        lastPad = cpad;
        cpad = g->cd(nPad);
      } while(cpad!=lastPad);

      nPad=mPad;

      if(nPad==1 && fPadMap[0].size()>0) {
        g->Clear();
        g->DivideSquare(2);
        g->cd(1); 
        for(auto obj : fPadMap[0]) 
          obj->Draw("same");
        g->cd(2);
      } else if(nPad>1) {
        g->Clear();
        g->DivideSquare(nPad);
        for(int i=1;i<=nPad;i++) {
          g->cd(i); 
          for(size_t j=0;j<fPadMap[i].size();j++) { 
            fPadMap[i][j]->Draw("same");
            //printf("pad[%i]: %s\n",i,fPadMap[i][j]->GetName());
          }
        }
        g->cd(nPad+1);
      }

      current->Draw();
      g->cd();
    }
  } else if(hists2D.size()>1) {

    //ignoring combo box for the momenet.
    if(!g)
      g = gROOT->MakeDefCanvas(); //new GCanvas;
    g->DivideSquare(hists2D.size());
    for(int i=0;i<(int)hists2D.size();i++) {
      g->cd(i+1);
      hists2D[i]->Draw();
    }
    g->cd(1);

    //if(fDrawComboBox->GetSelected() != EDrawOption::eDrawNew)
    //  g = gROOT->MakeDefCanvas(); //new GCanvas;
    //THStack hs;
    //for(auto it=hists1D.begin();it!=hists1D.end();it++) 
    //  hs.Add(*it);
    //hs.Draw("pads");
  }

  if(graphs.size()>0) {
    for(int i=0;i<(int)graphs.size();i++) {
      if(graphs[i]->InheritsFrom(TCutG::Class()))
        graphs[i]->Draw("same");
      else 
        graphs[i]->Draw("A*");
    }
  }

  doUpdate(); 
}

void Histomatic::drawHists(std::vector<TH1*> hists, TCanvas *g) const {

  if(fDrawComboBox->GetSelected() == EDrawOption::eDrawSame ||
      fDrawComboBox->GetSelected() == EDrawOption::eDrawStacked) {
    TList *found = GrabHists(gPad->GetCanvas());
    if(found->GetEntries() > 0 ) { 
      TIter iter(found);
      while(TH1 *hist = (TH1*)iter.Next()) {
        hists.push_back(hist);
      }
      g->Clear();
    }
  }

  if(hists.size()==1) {   
    if(!g) 
      g = gROOT->MakeDefCanvas(); //new GCanvas;
    hists.at(0)->Draw();
    return;
  } 

  int ic;
  switch(fDrawComboBox->GetSelected()) {
    case EDrawOption::eDrawNew:
      for (auto it=hists.begin();it!=hists.end();it++) {
        drawHists(std::vector<TH1*>(it,it+1),0);
      } 
      break;
    case EDrawOption::eDrawStacked:
      //case EDrawOption::eDrawSame:
      //THStack *hs = new THStack("hs","");
      //ic = gStyle->GetColorPalette(0); 
      for(auto it=hists.begin();it!=hists.end();it++) {
        //gPad->IncrementPaletteColor(ic,"plc");
        (*it)->SetLineColor(ic++);
        if(GrabHist()) 
          (*it)->Draw("same");
        else 
          (*it)->Draw("");
      }
      gPad->BuildLegend(0.75,0.75,0.95,0.95,"");
      break;  
    case EDrawOption::eDrawSame:
      //case EDrawOption::eDrawStacked:
      if(hists.size()<=5) {
        g->Divide(1,hists.size(),0.01,0);
        int padN=1;
        for(auto it=hists.begin();it!=hists.end();it++) {
          g->cd(padN++);
          (*it)->Draw();
        }
        g->cd(1);
      } else {
        THStack *hs = new THStack(Form("hs_%s",gPad->GetTitle()),"");
        for(auto it=hists.begin();it!=hists.end();it++) 
          hs->Add(*it);
        //printf("hs.GetNhists() = %i\n",hs.GetNhists());
        hs->Draw("pads");
        g->cd(1);
      }
      break;
    default:
      break;
  };
}


void Histomatic::closeAllCanvases() {
  TSeqCollection *l = gROOT->GetListOfCanvases();
  TIter iter(l);
  while(TObject *obj = iter.Next()) {
    TCanvas *c = (TCanvas*)obj;
    c->Close();
  }
  return;
}


void Histomatic::doUpdate() const {
  if(gPad) {
    gPad->Modified();
    gPad->Update();
  }
}



//TList *Histomatic::GetAllActive() { return fGListTree->GetAllActive(); } 
