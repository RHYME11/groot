
#include <GCommands.h>

#include<TROOT.h>
#include<TVirtualPad.h>
#include<TList.h>
#include<TH1.h>
#include<TGraph.h>
#include<TF1.h>
#include<THStack.h>
#include<TCanvas.h>
#include<TSpectrum.h>
#include<TPolyMarker.h>
#include<TText.h>

#include<GGaus.h>
#include<GPeak.h>


GGaus *GausFit(TH1 *hist,double xlow, double xhigh,Option_t *opt) {
  if(!hist)
    return 0;
  if(xlow>xhigh)
    std::swap(xlow,xhigh);

  GGaus *mypeak= new GGaus(xlow,xhigh);
  std::string options = opt;
  options.append("Q+");
  mypeak->Fit(hist,options.c_str());
  //mypeak->Background()->Draw("SAME");
  //TF1 *bg = new TF1(*mypeak->Background());
  //hist->GetListOfFunctions()->Add(bg);

  double chi2 = GetChi2(hist,mypeak);
  printf("Cal chi2 = %.03f\n",chi2);
  
  return mypeak;
}

GPeak *PhotoPeakFit(TH1 *hist,double xlow, double xhigh,Option_t *opt) {
  if(!hist)
    return 0;
  if(xlow>xhigh)
    std::swap(xlow,xhigh);

  GPeak *mypeak= new GPeak((xlow+xhigh)/2.0,xlow,xhigh);
  std::string options = opt;
  options.append("Q+");
  mypeak->Fit(hist,options.c_str());
  //mypeak->Background()->Draw("SAME");
  //TF1 *bg = new TF1(*mypeak->Background());
  //hist->GetListOfFunctions()->Add(bg);

  double chi2 = GetChi2(hist,mypeak);
  printf("Cal chi2 = %.03f\n",chi2);

  return mypeak;
}


TH1 *GrabHist(int i)  {
  //return the histogram from the current canvas, pad i.
  TH1 *hist = 0;
  if(!gPad)
    return hist;
  TIter iter(gPad->GetListOfPrimitives());
  int j=0;
  while(TObject *obj = iter.Next()) {
    if(obj->InheritsFrom(TH1::Class())) {
      if(j==i) {
        hist = (TH1*)obj;
        break;
      }
      j++;
    } else if(obj->InheritsFrom(THStack::Class())) {
      hist = ((THStack*)obj)->GetHistogram();
      break;
    }
  }
  return hist;
}




TList *GrabHists(TVirtualPad *p) {
  //return all histograms on a canvas or pad. (default is the gPad);
  TList *histList = new TList;
  if(!p) p = gPad;
  if(!p) return histList;

  TVirtualPad *current = gPad;
  //TCanvas *c = p->GetCanvas();
  
  TIter nextp(p->GetListOfPrimitives());
  while(TObject *obj = nextp()) {
    //printf("obj->GetName() = %s\n",obj->GetName());
    if(obj->InheritsFrom(TVirtualPad::Class())) {
     TList *temp = GrabHists((TVirtualPad*)obj);
     TIter nextp2(temp);
     while(TObject *obj2 = nextp2()) histList->Add(obj2);
    }else if(obj->InheritsFrom(TH1::Class())) {
      histList->Add(obj);
    }
  }  
  //printf("found %i histograms...\n",histList->GetEntries());
  return histList;
}

TF1 *GrabFit(int i)  {
  //return the histogram from the current canvas, pad i.
  TH1 *hist = 0;
  TF1 *fit = 0;
  if(!gPad)
    return fit;
  TIter iter(gPad->GetListOfPrimitives());
  int j=0;
  while(TObject *obj = iter.Next()) {
    if(obj->InheritsFrom(TH1::Class())) {
      hist = (TH1*)obj;
      TIter iter2(hist->GetListOfFunctions());
      while(TObject *obj2 = iter2.Next()){
	if(obj2->InheritsFrom(TF1::Class())){
	  if(j==i) {
	    fit=(TF1*)obj2;
	    return fit;
	  }
	  j++;
	}
      }
    }
  }
  return fit;
}



double GetChi2(TObject *obj,TF1 *f=0) {
  if(obj->InheritsFrom(TGraph::Class())) {
    /*
    TGraph *gr = (TGraph*)obj;
    if(f==0) 
      if(gr->GetListOfFunctions()->GetEntries()) 
        f = (TF1*)gr->GetListOfFunctions()->Last();
    if(f==0)
      return sqrt(-1);

    double low,high;
    f->GetRange(low,high);
    double *x = gr->GetX();
    double *y = gr->GetY();
    double chi2 = 0;
    for(int i=0;i<11;i++){            //help 
      double yf = fx->Eval(x[i]);
      chi2 += pow((yf-y[i]),2)/y[i];
    }

    printf("chi2 = %f\n",chi2/10); // 10 is number of dat
    */
    //FIXME
    return sqrt(-1);

  } else if(obj->InheritsFrom(TH1::Class())) {
    TH1 *hist = (TH1*)obj;
    if(hist->GetDimension()>1)
      return sqrt(-1);

    if(f==0) 
      if(hist->GetListOfFunctions()->GetEntries()) 
        f = (TF1*)hist->GetListOfFunctions()->Last();
    if(f==0)
      return sqrt(-1);

    double low,high;
    f->GetRange(low,high);
    int binLow = hist->FindBin(low);
    int binHigh = hist->FindBin(high);

    //double *x = gr->GetX();
    //double *y = gr->GetY();
    double chi2 = 0;
    int pzero = 0;
    for(int i=binLow;i<=binHigh;i++){            //help 
      double obs = hist->GetBinContent(i);
      double cal = f->Eval(hist->GetBinCenter(i));
      if(cal==0){
        pzero++;
        continue;
      }
      chi2 += pow((obs-cal),2)/cal;
    }
    int NDF = binHigh-binLow+1-pzero;
    //printf("chi2 = %f\n",chi2/(NDF-1)); // 10 is number of dat
    
    return chi2/(NDF-1);
  }
  return sqrt(-1);
}

int LabelPeaks(TH1 *hist,double sigma,double thresh,Option_t *opt) {
  TSpectrum::StaticSearch(hist,sigma,"Qnodraw",thresh);
  TPolyMarker *pm = (TPolyMarker*)hist->GetListOfFunctions()->FindObject("TPolyMarker");
  if(!pm) {
    //something has gone wrong....
    return 0;
  }
  TObjArray *array = (TObjArray*)hist->GetListOfFunctions()->FindObject("PeakLabels");
  if(array) {
    hist->GetListOfFunctions()->Remove((TObject*)array);
    array->Delete();
  }
  array = new TObjArray();
  array->SetName("PeakLabels");
  int n = pm->GetN();
  if(n==0)
    return n;
  TText *text;
  double *x = pm->GetX();
  //  double *y = pm->GetY();
  for(int i=0;i<n;i++) {
    //y[i] += y[i]*0.15;
    double y = 0;
    for(int i_x = x[i]-3;i_x<x[i]+3;i_x++){
      if((hist->GetBinContent(hist->GetXaxis()->FindBin(i_x)))>y){
        y = hist->GetBinContent(hist->GetXaxis()->FindBin(i_x));
      }
    }
    y+=y*0.1;
    text = new TText(x[i],y,Form("%.1f",x[i]));
    text->SetTextSize(0.025);
    text->SetTextAngle(90);
    text->SetTextAlign(12);
    text->SetTextFont(42);
    text->SetTextColor(hist->GetLineColor());
    array->Add(text);
  }
  hist->GetListOfFunctions()->Remove(pm);
  pm->Delete();
  hist->GetListOfFunctions()->Add(array);
  return n;
}



bool RemovePeaks(TH1 *hists) {
  bool flag = false;
  //for(unsigned int i=0;i<nhists;i++) {
    if(TObject *obj=hists->GetListOfFunctions()->FindObject("PeakLabels")) {
      hists->GetListOfFunctions()->Remove(obj);
      ((TObjArray*)obj)->Delete();
      flag = true;
    }
  //}
  return flag;
}




