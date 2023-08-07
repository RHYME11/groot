#ifndef __GFUNCTIONS_H__
#define __GFUNCTIONS_H__

#include "TMath.h"

namespace GFunctions {

  Double_t LinFit(Double_t *dim, Double_t *par);
  Double_t QuadFit(Double_t *dim, Double_t *par);
  
  Double_t PolyBg(Double_t *dim, Double_t *par,Int_t order);
  
  Double_t StepBG(Double_t *dim, Double_t *par);
  Double_t StepFunction(Double_t *dim, Double_t *par);
  Double_t PhotoPeak(Double_t *dim, Double_t *par);
  Double_t PhotoPeakBG(Double_t *dim, Double_t *par);
  Double_t PhotoPeakBGExcludeRegion(Double_t *dim,Double_t *par);
  Double_t Gaus(Double_t *dim, Double_t *par);
  Double_t DoubleGaus(Double_t *dim, Double_t *par);
  Double_t SkewedGaus(Double_t *dim, Double_t *par);
  Double_t Efficiency(Double_t *dim, Double_t *par);

  Double_t GausExpo(Double_t *dim,Double_t *par);
  
  Double_t LanGaus(Double_t *dim,Double_t *par);
  Double_t LanGausHighRes(Double_t *dim,Double_t *par);

  Double_t GammaEff(Double_t *dim,Double_t *par);
  Double_t AlignedAD(Double_t *x,Double_t *par);
  Double_t AlignedAD_Norm(Double_t *x,Double_t *par);
  Double_t AlignedADPol_Norm(Double_t *x,Double_t *par);

  Double_t ComptonEnergy(Double_t *x,Double_t *par);
  Double_t ComptonAngle(Double_t *x,Double_t *par);
  Double_t ComptonRatio(Double_t *x,Double_t *par);
  Double_t AnalyzingPower(Double_t *x,Double_t *par);
  Double_t Polarization(Double_t *x,Double_t *par);
  Double_t PolarizationAsymmetry(Double_t *x,Double_t *par);
  Double_t KN_unpol(Double_t *x,Double_t *par);
  Double_t KN_unpol_norm(Double_t *x,Double_t *par);
  Double_t KN_pol(Double_t *x,Double_t *par);


}

#endif
