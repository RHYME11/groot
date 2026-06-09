#ifndef __GNUCLEUS_H__
#define __GNUCLEUS_H__

#include <TList.h>
#include <TNamed.h>

#include <string>

#include <GTransition.h>

class GNucleus : public TNamed {
  public:
    GNucleus();
    GNucleus(const char* symbol);
    GNucleus(int z, int n, double mass, const char* symbol);
    GNucleus(int z, int n, const char* massFile = nullptr);
    virtual ~GNucleus();

    static std::string SortName(const char* name);

    void SetZ(int z);
    void SetN(int n);
    void SetMassExcess(double massExcess);
    void SetMass(double mass);
    void SetMass();
    void SetSymbol(const char* symbol);

    void AddTransition(double energy,
                       double intensity,
                       double energyUncertainty = 0.0,
                       double intensityUncertainty = 0.0);
    void AddTransition(GTransition* transition);

    int GetZ() const { return fZ; }
    int GetN() const { return fN; }
    int GetA() const { return fA; }

    double GetMassExcess() const { return fMassExcess; }
    double GetMass() const { return fMass; }

    const char* GetSymbol() const { return fSymbol.c_str(); }

    double GetEnergyFromBeta(double beta) const;
    double GetBetaFromEnergy(double energyMeV) const;

    GTransition* GetTransition(int index) const;

    int NTransitions() const { return fTransitionList.GetSize(); }
    int GetNTransitions() const { return fTransitionList.GetSize(); }

    double GetRadius() const;
    int GetZfromSymbol(const char* symbol);

    void Print(Option_t* opt = "") const override;
    void WriteSourceFile(std::string outfilename = "") const;

    TList* GetTransitionList() { return &fTransitionList; }
    const TList* GetTransitionList() const { return &fTransitionList; }

  private:
    static std::string SourceDataDirectory();
    static std::string MassFile();
    static bool ParseName(const char* name, int& massNumber, std::string& symbol);

    void Init();
    void SetNameFromSymbol();
    bool LoadMassFile(const char* massFile);
    bool LoadTransitionFile();

    int fA;
    int fN;
    int fZ;
    double fMass;
    double fMassExcess;
    std::string fSymbol;

    TList fTransitionList;

  ClassDefOverride(GNucleus, 1)
};

#endif
