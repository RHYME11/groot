#include <GTransition.h>

#include <TString.h>

ClassImp(GTransition)

GTransition::GTransition() {
  Clear();
}

GTransition::~GTransition() {
}

void GTransition::Clear(Option_t* opt) {
  TObject::Clear(opt);

  fEnergy = 0;
  fEngUncertainty = 0;
  fIntensity = 0;
  fIntUncertainty = 0;
}

void GTransition::Print(Option_t* opt) const {
  if(!std::isnan(fEngUncertainty))
    printf("Energy: %.02f +/- %.02f", fEnergy, fEngUncertainty);
  else
    printf("Energy: %.02f", fEnergy);

  if(!std::isnan(fIntensity)) {
    if(!std::isnan(fIntUncertainty))
      printf("\tIntensity: %.02f +/- %.02f\n", fIntensity, fIntUncertainty);
    else
      printf("\tIntensity: %.02f\n", fIntensity);
  } else {
    printf("\n");
  }
}

std::string GTransition::PrintToString() const {
  return Form("%f\t%f\t%f\t%f",
              fEnergy,
              fEngUncertainty,
              fIntensity,
              fIntUncertainty);
}

int GTransition::CompareIntensity(const TObject* obj) const {
  const auto* other = static_cast<const GTransition*>(obj);

  if(fIntensity > other->fIntensity)
    return -1;
  if(fIntensity == other->fIntensity)
    return 0;
  return 1;
}

int GTransition::Compare(const TObject* obj) const {
 return CompareIntensity(obj);

 // const auto* other = static_cast<const GTransition*>(obj);

 // if(fEnergy < other->fEnergy)
 //   return -1;
 // if(fEnergy == other->fEnergy)
 //   return 0;
 // return 1;
}
