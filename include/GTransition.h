<<<<<<< HEAD
#ifndef __GTRANSITION_H__
#define __GTRANSITION_H__

#include <TObject.h>

#include <cmath>
#include <cstdio>
#include <string>

class GTransition : public TObject {
  public:
    GTransition();
    virtual ~GTransition();

    Bool_t IsSortable() const override { return kTRUE; }

    int Compare(const TObject* obj) const override;
    int CompareIntensity(const TObject* obj) const;

    void SetEnergy(double energy) { fEnergy = energy; }
    void SetEnergyUncertainty(double error) { fEngUncertainty = error; }
    void SetIntensity(double intensity) { fIntensity = intensity; }
    void SetIntensityUncertainty(double error) { fIntUncertainty = error; }

    double GetEnergy() const { return fEnergy; }
    double GetEnergyUncertainty() const { return fEngUncertainty; }
    double GetIntensity() const { return fIntensity; }
    double GetIntensityUncertainty() const { return fIntUncertainty; }

    void Clear(Option_t* opt = "") override;
    void Print(Option_t* opt = "") const override;

    std::string PrintToString() const;

  private:
    double fEnergy;
    double fEngUncertainty;
    double fIntensity;
    double fIntUncertainty;

  ClassDefOverride(GTransition, 0)
};

#endif
