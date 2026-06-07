#ifndef __GCALIBRATION_H__
#define __GCALIBRATION_H__

#include <TObject.h>

#include <string>
#include <vector>

struct GCalibrationPoint {
  double channel = 0;
  double channelErr = 0;
  double energy = 0;
  double energyErr = 0;
  bool enabled = true;
};

class GCalibration : public TObject {
  public:
    GCalibration();
    virtual ~GCalibration();

    void Clear(Option_t* opt="") override;

    bool LoadPoints(const char* pointsFile);
    bool Fit(int order=1, const char* unit="keV");
    bool Save(const char* calFile) const;
    bool Load(const char* calFile);

    void AddPoint(double channel,
                  double energy,
                  double channelErr=0,
                  double energyErr=0,
                  bool enabled=true);

    double Eval(double channel) const;
    double Residual(size_t index) const;

    int GetOrder() const { return fOrder; }
    const char* GetUnit() const { return fUnit.c_str(); }

    double GetC0() const { return fC0; }
    double GetC1() const { return fC1; }
    double GetC2() const { return fC2; }

    const std::vector<GCalibrationPoint>& GetPoints() const { return fPoints; }
    const std::vector<double>& GetResiduals() const { return fResiduals; }

  private:
    std::vector<GCalibrationPoint> fPoints;
    std::vector<double> fResiduals;

    int fOrder;
    double fC0;
    double fC1;
    double fC2;
    std::string fUnit;

  ClassDefOverride(GCalibration, 1)
};

#endif

