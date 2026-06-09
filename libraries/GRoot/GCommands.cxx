
#include <GCommands.h>

#include <TGClient.h>
#include <TGMsgBox.h>

#include<TROOT.h>
#include<TVirtualPad.h>
#include<TList.h>
#include<TH1.h>
#include<TGraph.h>
#include<TF1.h>
#include<THStack.h>
#include<KeySymbols.h>
#include<TFile.h>
#include <TCutG.h>
#include <TGraphErrors.h>

#include <TObjString.h>
#include <TParameter.h>
#include <TAxis.h>

#include <TSpectrum.h>
#include <algorithm>
#include <cmath>
#include <utility>
#include <fstream>
#include <vector>
#include <limits>
#include <map>
#include <numeric>
#include <string>

#include<GCanvas.h>
#include<GGaus.h>
#include<GPeak.h>
#include<GMarker.h>
#include<GROI.h>
#include<GCalibration.h>
#include<GNucleus.h>
#include<GTransition.h>

#include<GH1D.h>
#include<GH2D.h>


namespace {
  GFitResultCallback gFitResultCallback = nullptr;

  void NotifyFitResult(TH1* hist, GGaus* fit) {
    if(!gFitResultCallback || !hist || !fit) return;

    double xlow = 0;
    double xhigh = 0;
    fit->GetRange(xlow, xhigh);

    GFitResult result;
    result.histName = hist->GetName();
    result.fitName = fit->GetName();
    result.centroid = fit->GetCentroid();
    result.fwhm = fit->GetFWHM();
    result.area = fit->GetArea();
    result.areaErr = fit->GetAreaErr();
    result.chi2 = fit->GetChi2();
    result.ndf = fit->GetNdf();
    result.xlow = xlow;
    result.xhigh = xhigh;

    gFitResultCallback(result);
  }

///// something very stupid. if it doesn't work replace with MatchSourceLines(....)  

  struct SourceLineCandidate {
    double energy = 0;
    double energyUncertainty = 0;
    double intensity = 0;
    double intensityUncertainty = 0;
  };

  struct PeakCandidate {
    double centroid = 0;
    double centroidErr = 0;
    double area = 0;
    double areaErr = 0;
    double fwhm = 0;
    double chi2Ndf = 0;
    double height = 0;
  };

  struct CalibrationAssignment {
    PeakCandidate peak;
    SourceLineCandidate line;
    double calibratedEnergy = 0;
    double residual = 0;
  };

  struct CalibrationSolution {
    bool valid = false;
    bool ambiguous = false;
    std::string rejectionReason;

    double c0 = 0;
    double c1 = 1;
    double c2 = 0;

    double rmsResidual = std::numeric_limits<double>::max();
    double score = std::numeric_limits<double>::max();
    double secondBestScore = std::numeric_limits<double>::max();

    std::vector<CalibrationAssignment> assignments;
  };

  bool ComparePeakByCentroid(const PeakCandidate& lhs,
                             const PeakCandidate& rhs) {
    return lhs.centroid < rhs.centroid;
  }

  bool ComparePeakByArea(const PeakCandidate& lhs,
                         const PeakCandidate& rhs) {
    return lhs.area > rhs.area;
  }

  bool CompareSourceLineByEnergy(const SourceLineCandidate& lhs,
                                 const SourceLineCandidate& rhs) {
    return lhs.energy < rhs.energy;
  }

  bool CompareSourceLineByIntensity(const SourceLineCandidate& lhs,
                                    const SourceLineCandidate& rhs) {
    return lhs.intensity > rhs.intensity;
  }

  void BuildIndexCombinations(int n,
                              int k,
                              int start,
                              std::vector<int>& current,
                              std::vector<std::vector<int>>& output) {
    if(static_cast<int>(current.size()) == k) {
      output.push_back(current);
      return;
    }

    for(int i = start; i < n; ++i) {
      current.push_back(i);
      BuildIndexCombinations(n, k, i + 1, current, output);
      current.pop_back();
    }
  }

  std::vector<std::vector<int>> BuildIndexCombinations(int n, int k) {
    std::vector<std::vector<int>> output;
    std::vector<int> current;

    if(k <= 0 || n <= 0 || k > n)
      return output;

    BuildIndexCombinations(n, k, 0, current, output);
    return output;
  }

  bool FitLinearCalibration(const std::vector<PeakCandidate>& peaks,
                            const std::vector<SourceLineCandidate>& lines,
                            double& c0,
                            double& c1) {
    if(peaks.size() != lines.size() || peaks.size() < 2)
      return false;

    double sumX = 0;
    double sumY = 0;
    double sumXX = 0;
    double sumXY = 0;

    for(size_t i = 0; i < peaks.size(); ++i) {
      const double x = peaks[i].centroid;
      const double y = lines[i].energy;

      sumX += x;
      sumY += y;
      sumXX += x * x;
      sumXY += x * y;
    }

    const double n = static_cast<double>(peaks.size());
    const double denominator = n * sumXX - sumX * sumX;

    if(std::abs(denominator) < 1e-12)
      return false;

    c1 = (n * sumXY - sumX * sumY) / denominator;
    c0 = (sumY - c1 * sumX) / n;

    return std::isfinite(c0) && std::isfinite(c1);
  }

  double CalibrationRmsResidual(const std::vector<PeakCandidate>& peaks,
                                const std::vector<SourceLineCandidate>& lines,
                                double c0,
                                double c1) {
    if(peaks.empty() || peaks.size() != lines.size())
      return std::numeric_limits<double>::max();

    double sumResidual2 = 0;

    for(size_t i = 0; i < peaks.size(); ++i) {
      const double calibrated = c0 + c1 * peaks[i].centroid;
      const double residual = calibrated - lines[i].energy;
      sumResidual2 += residual * residual;
    }

    return std::sqrt(sumResidual2 / static_cast<double>(peaks.size()));
  }

  bool HasReasonableLinearCalibration(const std::vector<PeakCandidate>& peaks,
                                      const std::vector<SourceLineCandidate>& lines,
                                      double c0,
                                      double c1) {
    if(peaks.size() < 2 || lines.size() < 2)
      return false;

    if(!std::isfinite(c0) || !std::isfinite(c1))
      return false;

    if(c1 <= 0)
      return false;

    const auto peakMinMax = std::minmax_element(peaks.begin(), peaks.end(),
                                                ComparePeakByCentroid);
    const auto lineMinMax = std::minmax_element(lines.begin(), lines.end(),
                                                CompareSourceLineByEnergy);

    const double peakSpan = peakMinMax.second->centroid - peakMinMax.first->centroid;
    const double lineSpan = lineMinMax.second->energy - lineMinMax.first->energy;

    if(peakSpan <= 0 || lineSpan <= 0)
      return false;

    const double roughSlope = lineSpan / peakSpan;

    if(c1 < roughSlope * 0.25 || c1 > roughSlope * 4.0)
      return false;

    const double lowCal = c0 + c1 * peakMinMax.first->centroid;
    const double highCal = c0 + c1 * peakMinMax.second->centroid;

    if(lowCal > highCal)
      return false;

    if(highCal < 0)
      return false;

    return true;
  }

  CalibrationSolution EvaluateLinearAssignment(std::vector<PeakCandidate> peaks,
                                               std::vector<SourceLineCandidate> lines) {
    CalibrationSolution solution;

    if(peaks.size() != lines.size() || peaks.size() < 2) {
      solution.rejectionReason = "assignment has too few matched points";
      return solution;
    }

    std::sort(peaks.begin(), peaks.end(), ComparePeakByCentroid);
    std::sort(lines.begin(), lines.end(), CompareSourceLineByEnergy);

    if(!FitLinearCalibration(peaks, lines, solution.c0, solution.c1)) {
      solution.rejectionReason = "linear fit failed";
      return solution;
    }

    if(!HasReasonableLinearCalibration(peaks, lines, solution.c0, solution.c1)) {
      solution.rejectionReason = "linear fit failed reasonableness checks";
      return solution;
    }

    solution.rmsResidual = CalibrationRmsResidual(peaks, lines,
                                                  solution.c0,
                                                  solution.c1);

    double fitPenalty = 0;
    for(const auto& peak : peaks) {
      if(peak.chi2Ndf > 0)
        fitPenalty += std::min(peak.chi2Ndf, 50.0) * 0.02;
    }

    
    solution.score = solution.rmsResidual
                   + fitPenalty
                   - 3.0 * static_cast<double>(peaks.size());


    solution.assignments.clear();
    for(size_t i = 0; i < peaks.size(); ++i) {
      CalibrationAssignment assignment;
      assignment.peak = peaks[i];
      assignment.line = lines[i];
      assignment.calibratedEnergy = solution.c0 + solution.c1 * peaks[i].centroid;
      assignment.residual = assignment.calibratedEnergy - lines[i].energy;
      solution.assignments.push_back(assignment);
    }

    solution.valid = true;
    return solution;
  }

  CalibrationSolution FindBestLinearSourceCalibration(std::vector<PeakCandidate> peaks,
                                                      std::vector<SourceLineCandidate> lines,
                                                      int minMatches) {
    CalibrationSolution best;
    CalibrationSolution secondBest;

    std::sort(peaks.begin(), peaks.end(), ComparePeakByArea);
    std::sort(lines.begin(), lines.end(), CompareSourceLineByIntensity);

    const int maxPeaksToTry = std::min(static_cast<int>(peaks.size()), 8);
    const int maxLinesToTry = std::min(static_cast<int>(lines.size()), 12);

    peaks.resize(maxPeaksToTry);
    lines.resize(maxLinesToTry);

    std::sort(peaks.begin(), peaks.end(), ComparePeakByCentroid);
    std::sort(lines.begin(), lines.end(), CompareSourceLineByEnergy);

    const int maxMatches = std::min(maxPeaksToTry, maxLinesToTry);

    for(int matchCount = maxMatches; matchCount >= minMatches; --matchCount) {
      const auto peakCombos = BuildIndexCombinations(maxPeaksToTry, matchCount);
      const auto lineCombos = BuildIndexCombinations(maxLinesToTry, matchCount);

      for(const auto& peakCombo : peakCombos) {
        std::vector<PeakCandidate> selectedPeaks;
        selectedPeaks.reserve(matchCount);

        for(int index : peakCombo)
          selectedPeaks.push_back(peaks[index]);

        for(const auto& lineCombo : lineCombos) {
          std::vector<SourceLineCandidate> selectedLines;
          selectedLines.reserve(matchCount);

          for(int index : lineCombo)
            selectedLines.push_back(lines[index]);

          auto candidate = EvaluateLinearAssignment(selectedPeaks,
                                                    selectedLines);

          if(!candidate.valid)
            continue;

          if(!best.valid || candidate.score < best.score) {
            secondBest = best;
            best = candidate;
          } else if(!secondBest.valid || candidate.score < secondBest.score) {
            secondBest = candidate;
          }
        }
      }

      if(best.valid)
        break;
    }

    if(best.valid && secondBest.valid) {
      best.secondBestScore = secondBest.score;

      const bool sameMatchCount =
        best.assignments.size() == secondBest.assignments.size();

      const bool closeScore =
        secondBest.score < best.score + std::max(2.0, best.rmsResidual);

      best.ambiguous = sameMatchCount && closeScore;
    }

    if(best.valid) {
      const double maxAllowedRms = std::max(5.0, best.assignments.back().line.energy * 0.005);

      if(best.rmsResidual > maxAllowedRms) {
        best.valid = false;
        best.rejectionReason = Form("RMS residual %.3f exceeds limit %.3f",
                                    best.rmsResidual,
                                    maxAllowedRms);
      } else if(best.ambiguous) {
        best.valid = false;
        best.rejectionReason = "best calibration is ambiguous against second-best solution";
      }
    }

    return best;
  }

  std::map<double, double> MatchSourceLines(std::vector<double> peaks,
                                            std::vector<double> sourceLines) {
    std::vector<PeakCandidate> peakCandidates;
    peakCandidates.reserve(peaks.size());

    for(double peak : peaks) {
      PeakCandidate candidate;
      candidate.centroid = peak;
      candidate.area = 1.0;
      candidate.fwhm = 1.0;
      peakCandidates.push_back(candidate);
    }

    std::vector<SourceLineCandidate> lineCandidates;
    lineCandidates.reserve(sourceLines.size());

    for(double line : sourceLines) {
      SourceLineCandidate candidate;
      candidate.energy = line;
      candidate.intensity = 1.0;
      lineCandidates.push_back(candidate);
    }

    std::map<double, double> matched;
    auto solution = FindBestLinearSourceCalibration(peakCandidates,
                                                    lineCandidates,
                                                    2);

    if(!solution.valid)
      return matched;

    for(const auto& assignment : solution.assignments)
      matched[assignment.peak.centroid] = assignment.line.energy;

    return matched;
  }

}

void SetFitResultCallback(GFitResultCallback callback) {
  gFitResultCallback = callback;
}

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
  NotifyFitResult(hist, mypeak);

  return mypeak;
}

///////

int AutoFitPeaks(TH1* hist, double threshold, double sigma, int maxPeaks) {
  if(!hist)
    hist = GrabHist();

  if(!hist || hist->GetDimension() != 1) {
    std::cout << "AutoFitPeaks needs a 1D histogram." << std::endl;
    return 0;
  }

  TSpectrum spectrum(maxPeaks * 4);
  int found = spectrum.Search(hist, sigma, "goff nodraw", threshold);
  if(found <= 0) {
    std::cout << "AutoFitPeaks found no peaks." << std::endl;
    return 0;
  }

  struct PeakCandidate {
    double x;
    double y;
  };

  std::vector<PeakCandidate> candidates;
  double* positions = spectrum.GetPositionX();

  for(int i = 0; i < found; ++i) {
    int bin = hist->GetXaxis()->FindBin(positions[i]);
    candidates.push_back({positions[i], hist->GetBinContent(bin)});
  }

  std::sort(candidates.begin(), candidates.end(),
            [](const PeakCandidate& a, const PeakCandidate& b) {
              return a.y > b.y;
            });

  std::vector<std::pair<double,double>> acceptedRanges;
  int fitted = 0;

  double axisLow = hist->GetXaxis()->GetXmin();
  double axisHigh = hist->GetXaxis()->GetXmax();
  double axisWidth = axisHigh - axisLow;

  for(const auto& candidate : candidates) {
    if(fitted >= maxPeaks)
      break;

    int bin = hist->GetXaxis()->FindBin(candidate.x);
    double binWidth = hist->GetXaxis()->GetBinWidth(bin);
    double halfWidth = std::max(3.0 * sigma * binWidth, 0.01 * axisWidth);

    double xlow = std::max(axisLow, candidate.x - halfWidth);
    double xhigh = std::min(axisHigh, candidate.x + halfWidth);

    bool overlaps = false;
    for(const auto& range : acceptedRanges) {
      if(xlow < range.second && xhigh > range.first) {
        overlaps = true;
        break;
      }
    }

    if(overlaps) {
      std::cout << "\tskipped overlapping peak near " << candidate.x << std::endl;
      continue;
    }

    GGaus* fit = new GGaus(xlow, xhigh);
    fit->SetName(Form("auto_gaus_%d", fitted));
    fit->Fit(hist, "Q+no-print");

    double reducedChi2 = fit->GetNdf() != 0 ? fit->GetChi2() / fit->GetNdf() : 0;
    bool goodFit = fit->GetCentroid() >= xlow &&
                   fit->GetCentroid() <= xhigh &&
                   fit->GetArea() > 0 &&
                   fit->GetFWHM() > 0 &&
                   fit->GetNdf() > 0 &&
                   reducedChi2 < 25.0;

    if(!goodFit) {
      std::cout << "\trejected peak near " << candidate.x
                << " chi2/NDF=" << reducedChi2 << std::endl;
      delete fit;
      continue;
    }

    acceptedRanges.push_back({xlow, xhigh});
    NotifyFitResult(hist, fit);

    std::cout << "\tauto fit " << fit->GetName()
              << " centroid=" << fit->GetCentroid()
              << " area=" << fit->GetArea()
              << " FWHM=" << fit->GetFWHM()
              << " chi2/NDF=" << reducedChi2
              << std::endl;

    ++fitted;
  }

  if(gPad) {
    gPad->Modified();
    gPad->Update();
  }

  std::cout << "AutoFitPeaks fit " << fitted << " of " << found << " found peaks." << std::endl;
  return fitted;
}

///////

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

// making calibration files

bool MakeSourceCalibration(TH1* hist,
                           const char* source,
                           const char* calFile,
                           int order,
                           double sigma,
                           double threshold,
                           const char* unit) {
  if(!hist)
    hist = GrabHist();

  if(!hist || hist->GetDimension() != 1) {
    std::cout << "MakeSourceCalibration needs a 1D histogram." << std::endl;
    return false;
  }

  if(!source || std::strlen(source) == 0) {
    std::cout << "MakeSourceCalibration needs a source name, for example 60Co or 152Eu."
              << std::endl;
    return false;
  }

  if(order != 1) {
    std::cout << "Robust source calibration currently supports linear order=1 only."
              << std::endl;
    return false;
  }

  GNucleus nucleus(source);
  if(nucleus.GetNTransitions() == 0) {
    std::cout << "No source transitions found for " << source << std::endl;
    return false;
  }

  std::vector<SourceLineCandidate> sourceLines;
  sourceLines.reserve(nucleus.GetNTransitions());

  for(int i = 0; i < nucleus.GetNTransitions(); ++i) {
    auto* transition = nucleus.GetTransition(i);
    if(!transition)
      continue;

    if(transition->GetEnergy() <= 0)
      continue;

    SourceLineCandidate line;
    line.energy = transition->GetEnergy();
    line.energyUncertainty = transition->GetEnergyUncertainty();
    line.intensity = transition->GetIntensity();
    line.intensityUncertainty = transition->GetIntensityUncertainty();
    sourceLines.push_back(line);
  }

  if(sourceLines.empty()) {
    std::cout << "No usable source energies found for " << source << std::endl;
    return false;
  }

  TSpectrum spectrum(100);
  const int found = spectrum.Search(hist, sigma, "goff nodraw", threshold);
  if(found <= 0) {
    std::cout << "No peaks found in histogram " << hist->GetName() << std::endl;
    return false;
  }

  std::vector<PeakCandidate> peaks;
  peaks.reserve(found);

  double* positions = spectrum.GetPositionX();

  const double axisLow = hist->GetXaxis()->GetXmin();
  const double axisHigh = hist->GetXaxis()->GetXmax();
  const double axisWidth = axisHigh - axisLow;

  std::vector<std::pair<double, double>> acceptedRanges;

  for(int i = 0; i < found; ++i) {
    const double position = positions[i];
    const int bin = hist->GetXaxis()->FindBin(position);
    const double binWidth = hist->GetXaxis()->GetBinWidth(bin);

    const double halfWidth = std::max(3.0 * sigma * binWidth,
                                      0.01 * axisWidth);

    const double xlow = std::max(axisLow, position - halfWidth);
    const double xhigh = std::min(axisHigh, position + halfWidth);

    bool overlaps = false;
    for(const auto& range : acceptedRanges) {
      if(xlow < range.second && xhigh > range.first) {
        overlaps = true;
        break;
      }
    }

    if(overlaps)
      continue;

    GGaus* fit = new GGaus(xlow, xhigh);
    fit->SetName(Form("source_cal_peak_%d", static_cast<int>(peaks.size())));
    fit->Fit(hist, "Q+no-print");

    const double ndf = fit->GetNdf();
    const double chi2Ndf = ndf > 0 ? fit->GetChi2() / ndf : 0;
    const double centroid = fit->GetCentroid();
    const double area = fit->GetArea();
    const double fwhm = fit->GetFWHM();

    const bool goodFit =
      centroid >= xlow &&
      centroid <= xhigh &&
      area > 0 &&
      fwhm > 0 &&
      ndf > 0 &&
      chi2Ndf < 25.0;

    if(!goodFit) {
      std::cout << "Rejected peak near " << position
                << " chi2/NDF=" << chi2Ndf
                << " area=" << area
                << " FWHM=" << fwhm
                << std::endl;
      delete fit;
      continue;
    }

    PeakCandidate peak;
    peak.centroid = centroid;
    peak.centroidErr = fit->GetCentroidErr();
    peak.area = area;
    peak.areaErr = fit->GetAreaErr();
    peak.fwhm = fwhm;
    peak.chi2Ndf = chi2Ndf;
    peak.height = hist->GetBinContent(bin);

    peaks.push_back(peak);
    acceptedRanges.push_back({xlow, xhigh});

    std::cout << "Accepted peak centroid=" << peak.centroid
              << " area=" << peak.area
              << " FWHM=" << peak.fwhm
              << " chi2/NDF=" << peak.chi2Ndf
              << std::endl;

    delete fit;
  }

  const int minMatches = order == 2 ? 3 : 2;
  if(static_cast<int>(peaks.size()) < minMatches) {
    std::cout << "Only accepted " << peaks.size()
              << " peak fits; need at least " << minMatches
              << " for calibration."
              << std::endl;
    return false;
  }

  auto solution = FindBestLinearSourceCalibration(peaks, sourceLines, minMatches);

  if(!solution.valid) {
    std::cout << "Source calibration failed: "
              << solution.rejectionReason
              << std::endl;
    return false;
  }

  GCalibration calibration;

  for(const auto& assignment : solution.assignments) {
    calibration.AddPoint(assignment.peak.centroid,
                         assignment.line.energy);
  }
///
   if(!calibration.Fit(order, unit))
    return false;

  const bool coefficientsFinite =
    std::isfinite(calibration.GetC0()) &&
    std::isfinite(calibration.GetC1()) &&
    std::isfinite(calibration.GetC2());

  const bool coefficientsAllZero =
    calibration.GetC0() == 0 &&
    calibration.GetC1() == 0 &&
    calibration.GetC2() == 0;

  if(!coefficientsFinite) {
    std::cout << "Source calibration failed: non-finite calibration coefficients."
              << std::endl;
    return false;
  }

  if(calibration.GetC1() <= 0) {
    std::cout << "Source calibration failed: non-positive calibration slope C1="
              << calibration.GetC1()
              << std::endl;
    return false;
  }

  if(coefficientsAllZero) {
    std::cout << "Source calibration failed: calibration coefficients are all zero."
              << std::endl;
    return false;
  }


  const auto& points = calibration.GetPoints();
  const auto& residuals = calibration.GetResiduals();

  std::cout << "Source calibration report for " << source << std::endl;
  std::cout << "histogram: " << hist->GetName() << std::endl;
  std::cout << "candidate peaks accepted: " << peaks.size() << std::endl;
  std::cout << "matched lines: " << solution.assignments.size() << std::endl;
  std::cout << "RMS residual: " << solution.rmsResidual << " " << unit << std::endl;
  std::cout << "score: " << solution.score;
  if(solution.secondBestScore < std::numeric_limits<double>::max())
    std::cout << " second-best score: " << solution.secondBestScore;
  std::cout << std::endl;

  std::cout << "\tcentroid\tcalibrated\tsource\tresidual\tintensity\tarea\tFWHM\tchi2/NDF"
            << std::endl;

  for(size_t i = 0; i < solution.assignments.size(); ++i) {
    const auto& assignment = solution.assignments[i];
    const double calibrated = calibration.Eval(assignment.peak.centroid);
    const double residual = i < residuals.size()
                          ? residuals[i]
                          : calibrated - assignment.line.energy;

    std::cout << "\t"
              << assignment.peak.centroid << "\t"
              << calibrated << "\t"
              << assignment.line.energy << "\t"
              << residual << "\t"
              << assignment.line.intensity << "\t"
              << assignment.peak.area << "\t"
              << assignment.peak.fwhm << "\t"
              << assignment.peak.chi2Ndf
              << std::endl;
  }

  std::cout << "Calibration coefficients:" << std::endl;
  std::cout << "C0=" << calibration.GetC0()
            << ", C1=" << calibration.GetC1()
            << ", C2=" << calibration.GetC2()
            << ", unit=" << calibration.GetUnit()
            << std::endl;

  if(!calibration.Save(calFile))
    return false;

  std::cout << "Wrote source calibration file: " << calFile << std::endl;
  return true;
}



bool MakeCalibration(const char* pointsFile,
                     const char* calFile,
                     int order,
                     const char* unit) {
  GCalibration calibration;

  if(!calibration.LoadPoints(pointsFile))
    return false;

  if(!calibration.Fit(order, unit))
    return false;

  if(!calibration.Save(calFile))
    return false;

  std::cout << "Wrote calibration file: " << calFile << std::endl;
  std::cout << "C0=" << calibration.GetC0()
            << ", C1=" << calibration.GetC1()
            << ", C2=" << calibration.GetC2()
            << ", unit=" << calibration.GetUnit()
            << std::endl;

  const auto& residuals = calibration.GetResiduals();
  const auto& points = calibration.GetPoints();

  std::cout << "Calibration residuals:" << std::endl;
  for(size_t i = 0; i < points.size(); ++i) {
    std::cout << "\tchannel=" << points[i].channel
              << ", energy=" << points[i].energy
              << ", residual=" << residuals[i]
              << std::endl;
  }

  return true;
}

TGraphErrors* DrawCalibrationResiduals(const char* pointsFile,
                                       int order,
                                       const char* unit) {
  GCalibration calibration;

  if(!calibration.LoadPoints(pointsFile))
    return nullptr;

  if(!calibration.Fit(order, unit))
    return nullptr;

  const auto& points = calibration.GetPoints();
  const auto& residuals = calibration.GetResiduals();

  TGraphErrors* graph = new TGraphErrors(points.size());
  graph->SetName("calibration_residuals");
  graph->SetTitle(Form("Calibration Residuals;Channel;Residual (%s)",
                       calibration.GetUnit()));

  for(size_t i = 0; i < points.size(); ++i) {
    graph->SetPoint(i, points[i].channel, residuals[i]);
    graph->SetPointError(i, points[i].channelErr, points[i].energyErr);
  }

  graph->SetMarkerStyle(20);
  graph->SetMarkerSize(1.1);
  graph->Draw("AP");

  std::cout << "Drew calibration residuals for "
            << pointsFile << std::endl;

  return graph;
}


TH1D* ApplyCalibration(TH1* hist, const char* calFile, const char* name) {
  if(!hist)
    hist = GrabHist();

  if(!hist || hist->GetDimension() != 1) {
    std::cout << "ApplyCalibration needs a 1D histogram." << std::endl;
    return nullptr;
  }

  if(hist->GetListOfFunctions() &&
   (hist->GetListOfFunctions()->FindObject("C0") ||
    hist->GetListOfFunctions()->FindObject("C1") ||
    hist->GetListOfFunctions()->FindObject("C2"))) {
  std::cout << "Histogram " << hist->GetName()
            << " already appears to have calibration metadata. "
            << "Skipping to avoid applying calibration twice."
            << std::endl;
  return nullptr;
}

  GCalibration calibration;
  if(!calibration.Load(calFile))
    return nullptr;

  int nBins = hist->GetNbinsX();
  std::vector<double> binEdges;
  binEdges.reserve(nBins + 1);

  for(int i = 1; i <= nBins + 1; ++i) {
    double channelEdge = hist->GetXaxis()->GetBinLowEdge(i);
    binEdges.push_back(calibration.Eval(channelEdge));
  }

  TString histName = name ? name : Form("%s_calibrated", hist->GetName());
  TH1D* calibrated = new TH1D(histName.Data(),
                              hist->GetTitle(),
                              nBins,
                              binEdges.data());

  calibrated->SetDirectory(nullptr);
  calibrated->GetXaxis()->SetTitle(calibration.GetUnit());
  calibrated->GetYaxis()->SetTitle(hist->GetYaxis()->GetTitle());

  for(int bin = 1; bin <= nBins; ++bin) {
    calibrated->SetBinContent(bin, hist->GetBinContent(bin));
    calibrated->SetBinError(bin, hist->GetBinError(bin));
  }

  calibrated->GetListOfFunctions()->Add(new TParameter<double>("C0", calibration.GetC0()));
  calibrated->GetListOfFunctions()->Add(new TParameter<double>("C1", calibration.GetC1()));
  calibrated->GetListOfFunctions()->Add(new TParameter<double>("C2", calibration.GetC2()));
  calibrated->GetListOfFunctions()->Add(new TObjString(Form("unit=%s", calibration.GetUnit())));

  calibrated->Draw();

  std::cout << "Applied calibration from " << calFile
            << " to " << hist->GetName()
            << " -> " << calibrated->GetName()
            << std::endl;

  return calibrated;
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

TObject *GrabPlottable(int i) { 
  TObject *object = 0;
  if(!gPad)
    return object;
  TIter iter(gPad->GetListOfPrimitives());
  int j=0;
  while(TObject *obj = iter.Next()) {
    if(obj->InheritsFrom(TH1::Class())) {
      if(j==i) {
        object = obj;
        break;
      }
      j++;
    } else if(obj->InheritsFrom(TGraph::Class())) {
      if(j==i) {
        object = obj;
        break;
      }
      j++;
    }
  }
  return object;
}



void ls(int n) {
  if(gROOT->GetListOfFiles()->GetEntries()>n) {
    ((TFile*)(gROOT->GetListOfFiles()->At(n)))->ls();
  }
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

void SaveAllCuts(TH1 *hist,const char* fname,Option_t *opt) { 
  if(!hist || !hist->GetListOfFunctions()) return;
  std::vector<TObject*> cuts;
  TIter iter(hist->GetListOfFunctions());
  while(TObject *obj = iter.Next()) {
    if(obj->InheritsFrom(TCutG::Class())) 
      cuts.push_back(obj);
  }
  if(cuts.size()) {
    TDirectory *current = gDirectory;
    printf("opening file: %s with option %s\n",fname,opt);
    TFile *f = TFile::Open(fname,opt);
    for(size_t i=0;i<cuts.size();i++) {
      printf("\twriting %s\n",cuts.at(i)->GetName());
      cuts.at(i)->Write();
    }
    current->cd();
  }
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
      if(obs==0){
        pzero++;
        continue;
      }
      double cal = f->Eval(hist->GetBinCenter(i));
      chi2 += pow((obs-cal),2)/obs;
    }
    int NDF = binHigh-binLow+1-pzero;
    //printf("chi2 = %f\n",chi2/(NDF-1)); // 10 is number of dat

    return chi2/(NDF-1);
  }
  return sqrt(-1);
}


TH1D  *ResidualHist(TH1* hist, TF1* fit) { 
  if(!hist) return nullptr;
  if(!fit)  fit = GrabFit();
  if(!fit)  return nullptr;

  double xlow, xhigh;
  fit->GetRange(xlow, xhigh);

  auto* res = new TH1D(
    Form("%s_residuals", fit->GetName()),
    "Residuals;X;Data-Fit",
    hist->GetNbinsX(),
    hist->GetXaxis()->GetXmin(),
    hist->GetXaxis()->GetXmax()
  );

  res->SetDirectory(nullptr);

  for(int i = 1; i <= hist->GetNbinsX(); ++i) {
    const double x = hist->GetBinCenter(i);

    if(x < xlow || x > xhigh)
      continue;

    const double data = hist->GetBinContent(i);
    const double model = fit->Eval(x);
    const double err = hist->GetBinError(i);

    double r = data - model;
    //if(normalized && err > 0)
    //  r /= err;
    res->SetBinContent(i, r);
  }
  return res;
}

void   DrawResiduals(TH1* hist, TF1* fit,bool normalized) { }


GInteractionInfo BuildInteractionInfo()  {
  GInteractionInfo info;
  info.pad = gPad;
  if(!info.pad) return info;
  info.selected = info.pad->GetSelected();
  info.target   = GrabPlottable();
  info.event    = info.pad->GetEvent();
  info.px       = info.pad->GetEventX();
  info.py       = info.pad->GetEventY();
  info.x        = info.pad->PadtoX(info.pad->AbsPixeltoX(info.px));
  info.y        = info.pad->PadtoY(info.pad->AbsPixeltoY(info.py));
  return info;
}

namespace {
  GInteractionInfo gLastInteractionInfo;
}

const GInteractionInfo& GetLastInteractionInfo() {
  return gLastInteractionInfo;
}


// the below is meant to be added to a pad to make the 
// object interactable. This is currently being handled 
// automatically in the Creation of a GCanvas - it is 
// also passed to subpads using the GCanvas::Divide Method.
// need to be void to prevent useless printing.  
void GRootInteract() {

  GInteractionInfo info = BuildInteractionInfo();
  gLastInteractionInfo = info;

  if(info.pad && gPad && (gPad != info.pad->GetSelectedPad())) 
    return;

  DispatchInteraction(info);

  if(info.modified && info.pad) {
    info.pad->Modified();
    info.pad->Update();
  }
  gLastInteractionInfo = info;

  return;
}

bool DispatchInteraction(GInteractionInfo &info) {
  bool temp = false;
  if(auto* h = dynamic_cast<TH1*>(info.target))
    temp = GRootInteractHist(h,info);
  else if(auto* gr = dynamic_cast<TGraph*>(info.target))
    temp = GRootInteractGraph(gr,info);
  return temp;
}


bool GRootInteractGraph(TGraph *current,GInteractionInfo &info) {return true;}
bool GRootInteractGraphMouseButton(TGraph* current,GInteractionInfo &info) {return true;}
bool GRootInteractGraphKeyPress(TGraph* current,GInteractionInfo &info) {return true;}

bool GRootInteractHist(TH1 *current,GInteractionInfo &info) {
  switch(info.event) {
    case kNoEvent:           
      break;
    case kButton1Down:       
    case kButton2Down:       
    case kButton3Down:       
      GRootInteractHistMouseButton(current,info);
      break;
    case kKeyDown:           
    case kWheelUp:           
    case kWheelDown:         
      break;
    case kButton1Shift:      
      GRootInteractHistMouseButton(current,info);
      break;
    case kButton1ShiftMotion:
    case kButton1Up:         
    case kButton2Up:         
    case kButton3Up:         
    case kKeyUp:             
      break;
    case kButton1Motion:     
    case kButton2Motion:     
    case kButton3Motion:     
      break;
    case kKeyPress:          
      GRootInteractHistKeyPress(current,info);
      break;
    case kArrowKeyPress:     
    case kArrowKeyRelease:  
      break;
    case kButton1Locate:     
    case kButton2Locate:     
    case kButton3Locate:     
      break;
    case kESC:               
      break;
    case kMouseMotion:       
    case kMouseEnter:        
    case kMouseLeave:        
      break;
    case kButton1Double:     
    case kButton2Double:     
    case kButton3Double:     
      break;
    default:
      break;
  }
  return true;
}

bool GRootInteractHistMouseButton(TH1* currentHist,GInteractionInfo &info) {
  if(!currentHist) return false;


  if(info.selected && info.selected->InheritsFrom(GMarker::Class()))
    return false; //should allow GMarker::ExecuteEvent deal with this...

  switch(info.event) {
    case kButton1Down:       
      //if(GCanvas::GetCurrentEvent().fState & kKeyControlMask) {
      //if(GCanvas::GetCurrentEvent().fState && currentHist) {
      {
        const bool ctrlPressed = GCanvas::GetCurrentEvent().fState & kKeyControlMask;

        GMarker *marker = new GMarker();
        marker->SetType(GMarkerType::kPrimary);
        marker->AddTo(currentHist,info.x,info.y,ctrlPressed);
        //gPad->Modified();
        info.modified = true;
      //} else {
        //printf("button1\n");
      }
      break;
    case kButton2Down:       
      //printf("button2\n");
      break;
    case kButton3Down:       
      //printf("button3\n");
      break;
    case kButton1Shift:
      //printf("shiftbutton1\n");
      break;
  }

  return true;
}


void ShowKeyboardShortcutHelp() {
  const char* help =
    "g      Gaussian fit / 2D gate from markers\n"
    "b      Show background\n"
    "B      Toggle background\n"
    "s      Show peaks\n"
    "x      Project 2D histogram onto X\n"
    "y      Project 2D histogram onto Y\n"
    "p      Project selected marked range\n"
    "w      Rebin by 2\n"
    "q      Unbin by 2\n"
    "m      Remove markers\n"
    "o      Unzoom axes\n"
    "z/l    Toggle log scale\n"
    "?      Show this help\n"
    "r      Create ROI from two markers\n"
    "R      Delete all ROIs on current histogram\n"
    "a      Auto-fit peaks\n";

  new TGMsgBox(gClient->GetRoot(), gClient->GetRoot(),
               "Keyboard Shortcuts", help, kMBIconAsterisk, kMBOk);
}



bool GRootInteractHistKeyPress(TH1 *currentHist,GInteractionInfo &info) {
  //printf("key:  %i\t%i\t%i\n",event,px,py);
  std::vector<GMarker*> markers = GMarker::Get(currentHist,GMarkerType::kPrimary);
  switch(info.py) {
    case kKey_b:
      if(currentHist && currentHist->InheritsFrom(GH1D::Class())) {
        dynamic_cast<GH1D*>(currentHist)->ShowBackground();
        //gPad->Modified();
        info.modified = true;
      }
      break;

    case kKey_Question:
	ShowKeyboardShortcutHelp();
      break;

    case kKey_r:
      if(currentHist && currentHist->GetDimension() == 1 && markers.size() > 1) {
        if(GROI::CreateFromMarkers(currentHist)) {
          GMarker::RemoveAll(currentHist);
          info.modified = true;
        }
      }
      break;

    case kKey_R:
      if(currentHist) {
        GROI::RemoveAll(currentHist);
        info.modified = true;
      }
      break;

    case kKey_a:
      if(currentHist && currentHist->GetDimension() == 1) {
	AutoFitPeaks(currentHist);
	info.modified = true;
      }
      break;
    case kKey_B:
      if(currentHist && currentHist->InheritsFrom(GH1D::Class())) {
        dynamic_cast<GH1D*>(currentHist)->ToggleBackground();
        info.modified = true;
        //gPad->Modified();
      }
      break;
    case kKey_c:
      if(markers.size()>1 && currentHist->GetDimension()==1) {
        markers.at(0)->SetType(GMarkerType::kBackground);
        markers.at(1)->SetType(GMarkerType::kBackground);
        info.modified = true;
        //gPad->Modified();
      }
      break;
    case kKey_e:
      if(markers.size()>1) {
        double xlow  = markers.at(0)->X();
        double xhigh = markers.at(1)->X();
        if(xlow>xhigh) std::swap(xlow,xhigh);
        currentHist->GetXaxis()->SetRangeUser(xlow,xhigh);
        if(currentHist->GetDimension()==2) {
          double ylow  = markers.at(0)->Y();
          double yhigh = markers.at(1)->Y();
          if(ylow>yhigh) std::swap(ylow,yhigh);
          currentHist->GetYaxis()->SetRangeUser(ylow,yhigh);
        }
        GMarker::RemoveAll(currentHist);
        //gPad->Modified();
        info.modified = true;
      }
      break;
    case kKey_g:
      if(currentHist->GetDimension()==1 && markers.size()>1) {
        GausFit(currentHist,markers.at(0)->X(),markers.at(1)->X());
        //gPad->Modified();
        info.modified = true;
        //GMarker::RemoveAll(currentHist);
      } else if(currentHist->GetDimension()==2 && markers.size()>1) {
        static int gGateCounter = 0;
        TCutG *cut = GMarker::MakeTCutG(currentHist);    
        cut->SetName(Form("cut%i",gGateCounter++));

        currentHist->GetListOfFunctions()->Add(cut);
        GMarker::RemoveAll(currentHist);

        info.modified = true;

      }
      break;
    case kKey_m:
      GMarker::RemoveAll(currentHist);
      //gPad->Modified();
      info.modified = true;
      break;
    case kKey_n:
      if(currentHist) {
        TListIter iter(currentHist->GetListOfFunctions());
        std::vector<TF1*> funcs;
        while(TObject *obj=iter.Next()) {
          if(obj->InheritsFrom(TF1::Class()))
            funcs.push_back(((TF1*)obj));
        }
        for(auto i=funcs.begin();i!=funcs.end();i++)
          currentHist->GetListOfFunctions()->Remove(*i);
        GMarker::RemoveAll(currentHist);
        currentHist->Sumw2(false);
        if(currentHist->InheritsFrom(GH1D::Class())) {
          GH1D *gcurrentHist = dynamic_cast<GH1D*>(currentHist);
          gcurrentHist->RemovePeaks();
        }
        info.modified = true;
        //gPad->Modified();
      }
      break;
    case kKey_o:
      currentHist->GetXaxis()->UnZoom();
      if(currentHist->GetDimension()==2) 
        currentHist->GetYaxis()->UnZoom();
      //gPad->Modified();
      info.modified = true;
      break;
    case kKey_p: 
      if(currentHist->InheritsFrom(GH1D::Class())) {
        GH1D *ghist = dynamic_cast<GH1D*>(currentHist);
        std::vector<GMarker*> bgmarkers = GMarker::Get(currentHist,GMarkerType::kBackground);
        if(markers.size()==2) {
          double xlow  = markers.at(0)->X();
          double xhigh = markers.at(1)->X();
          if(xlow>xhigh) std::swap(xlow,xhigh);
          GH2D *parent = dynamic_cast<GH2D*>(ghist->GetParent());
          if(parent) {
            GH1D *proj =0;
            //printf("markers.size() = %i\n",int(markers.size()));
            //printf("bgmarkers.size() = %i\n",int(bgmarkers.size()));
            if(bgmarkers.size()==2) {
              double bgxlow  = bgmarkers.at(0)->X();
              double bgxhigh = bgmarkers.at(1)->X();
              if(bgxlow>bgxhigh) std::swap(bgxlow,bgxhigh);
              if(currentHist->TestBits(GH1D::kProjectionX))
                proj = parent->ProjectionY(xlow,xhigh,bgxlow,bgxhigh);
              else
                proj = parent->ProjectionX(xlow,xhigh,bgxlow,bgxhigh);
            } else {
              if(currentHist->TestBits(GH1D::kProjectionX))
                proj = parent->ProjectionY(xlow,xhigh);
              else
                proj = parent->ProjectionX(xlow,xhigh);
            }
            //new GCanvas;
            proj->Draw();
          }
        }
      }
      break;
    case kKey_w:
      if(currentHist->InheritsFrom(GH1D::Class())) {
        double xlow = currentHist->GetXaxis()->GetBinLowEdge(currentHist->GetXaxis()->GetFirst());
        double xup  = currentHist->GetXaxis()->GetBinUpEdge(currentHist->GetXaxis()->GetLast());
        currentHist->Rebin(2);
        currentHist->GetXaxis()->SetRangeUser(xlow,xup);
      }
      //gPad->Modified();
      info.modified = true;
      break;   
    case kKey_q:
      if(currentHist->InheritsFrom(GH1D::Class())) {
        GH1D *gcurrentHist = dynamic_cast<GH1D*>(currentHist);
        double xlow = gcurrentHist->GetXaxis()->GetBinLowEdge(gcurrentHist->GetXaxis()->GetFirst());
        double xup  = gcurrentHist->GetXaxis()->GetBinUpEdge(gcurrentHist->GetXaxis()->GetLast());
        gcurrentHist->Unbin(2);
        gcurrentHist->GetXaxis()->SetRangeUser(xlow,xup);
      }        
      //gPad->Modified();
      info.modified = true;
      break;   
    case kKey_s:
      if(currentHist->InheritsFrom(GH1D::Class())) {
        GH1D *gcurrentHist = dynamic_cast<GH1D*>(currentHist);
        gcurrentHist->ShowPeaks();
        //gPad->Modified();
        info.modified = true;
      }
      break;
    case kKey_x:
      if(currentHist->InheritsFrom(GH2D::Class())) {
        double xlow = currentHist->GetXaxis()->GetBinLowEdge(currentHist->GetXaxis()->GetFirst());
        double xup  = currentHist->GetXaxis()->GetBinUpEdge(currentHist->GetXaxis()->GetLast());
        currentHist->GetXaxis()->UnZoom();
        GH1D *px = dynamic_cast<GH2D*>(currentHist)->ProjectionX();
        px->SetBit(GH1D::kProjectionX,1);
        new GCanvas;
        px->Draw();
      }
      break;
    case kKey_y:
      if(currentHist->InheritsFrom(GH2D::Class())) {
        double ylow = currentHist->GetYaxis()->GetBinLowEdge(currentHist->GetYaxis()->GetFirst());
        double yup  = currentHist->GetYaxis()->GetBinUpEdge(currentHist->GetYaxis()->GetLast());
        currentHist->GetYaxis()->UnZoom();
        GH1D *py = dynamic_cast<GH2D*>(currentHist)->ProjectionY();
        py->SetBit(GH1D::kProjectionX,0);
        new GCanvas;
        py->Draw();
      }
      break;
    case kKey_l:
    case kKey_z:
      if(currentHist->GetDimension()==2) { 
        if(info.pad->GetLogz()) {
          currentHist->GetZaxis()->UnZoom();
          info.pad->SetLogz(0);
        } else {
          if(currentHist->GetMinimum()<0) 
            currentHist->GetZaxis()->SetRangeUser(0,currentHist->GetMaximum());
          info.pad->SetLogz(1);
        }
      } else {
        //printf("GetUymax:  %.2f\n",gPad->GetUxmax());
        if(info.pad->GetLogy()) {
          info.pad->SetLogy(0);
          //currentHist->GetYaxis()->SetRangeUser(0,gPad->GetUymax());
          currentHist->GetYaxis()->UnZoom();
        } else {
          if(info.pad->GetUymin()<0) 
            currentHist->GetYaxis()->SetRangeUser(0,gPad->GetUymax());
          info.pad->SetLogy(1);
        }
      }
      info.modified = true;
      //gPad->Modified();
    default:
      break;
  }

  return true;
}



