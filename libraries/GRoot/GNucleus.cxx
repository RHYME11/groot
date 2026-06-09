#include <GNucleus.h>

#include <TString.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

ClassImp(GNucleus)

namespace {
  constexpr double kAtomicMassUnitMeV = 931.494043;

  std::string Trim(const std::string& input) {
    const auto first = input.find_first_not_of(" \t\n\r");
    if(first == std::string::npos)
      return "";

    const auto last = input.find_last_not_of(" \t\n\r");
    return input.substr(first, last - first + 1);
  }

  std::string ToLower(std::string input) {
    std::transform(input.begin(), input.end(), input.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    return input;
  }

  std::string CanonicalSymbol(std::string symbol) {
    if(symbol.empty())
      return symbol;

    symbol = ToLower(symbol);
    symbol[0] = std::toupper(static_cast<unsigned char>(symbol[0]));
    return symbol;
  }

  bool EqualsIgnoreCase(const std::string& lhs, const std::string& rhs) {
    return ToLower(lhs) == ToLower(rhs);
  }
}

GNucleus::GNucleus() {
  Init();
}

GNucleus::GNucleus(const char* symbol) {
  Init();

  int massNumber = 0;
  std::string parsedSymbol;
  if(!ParseName(symbol, massNumber, parsedSymbol)) {
    std::cout << "Could not parse nucleus name: "
              << (symbol ? symbol : "") << std::endl;
    return;
  }

  fA = massNumber;
  fSymbol = parsedSymbol;
  LoadMassFile(MassFile().c_str());
  SetNameFromSymbol();
  LoadTransitionFile();
}

GNucleus::GNucleus(int z, int n, double mass, const char* symbol) {
  Init();

  fZ = z;
  fN = n;
  fA = fZ + fN;
  fMass = mass;
  fSymbol = CanonicalSymbol(symbol ? symbol : "");

  SetNameFromSymbol();
  LoadTransitionFile();
}

GNucleus::GNucleus(int z, int n, const char* massFile) {
  Init();

  fZ = z;
  fN = n;
  fA = fZ + fN;

  LoadMassFile(massFile ? massFile : MassFile().c_str());
  SetNameFromSymbol();
  LoadTransitionFile();
}

GNucleus::~GNucleus() {
  fTransitionList.Delete();
}

void GNucleus::Init() {
  fA = 0;
  fN = 0;
  fZ = 0;
  fMass = 0;
  fMassExcess = 0;
  fSymbol.clear();

  fTransitionList.SetOwner(kTRUE);
}

std::string GNucleus::SourceDataDirectory() {
  const char* gsys = std::getenv("GSYS");
  if(!gsys || std::strlen(gsys) == 0) {
    std::cout << "GSYS is not set. Cannot locate SourceData directory." << std::endl;
    return "";
  }

  return std::string(gsys) + "/libraries/SourceData";
}

std::string GNucleus::MassFile() {
  const std::string sourceDir = SourceDataDirectory();
  if(sourceDir.empty())
    return "";

  return sourceDir + "/mass.dat";
}


bool GNucleus::ParseName(const char* name, int& massNumber, std::string& symbol) {
  massNumber = 0;
  symbol.clear();

  if(!name)
    return false;

  std::string input = name;
  input.erase(std::remove_if(input.begin(), input.end(),
                             [](unsigned char c) { return std::isspace(c); }),
              input.end());

  if(input.empty())
    return false;

  if(input.length() < 2) {
    switch(std::tolower(static_cast<unsigned char>(input[0]))) {
      case 'p':
        input = "1H";
        break;
      case 'd':
        input = "2H";
        break;
      case 't':
        input = "3H";
        break;
      case 'a':
        input = "4He";
        break;
      default:
        return false;
    }
  }

  const auto firstDigit = input.find_first_of("0123456789");
  const auto firstLetter = input.find_first_not_of("0123456789");

  if(firstDigit == std::string::npos || firstLetter == std::string::npos)
    return false;

  if(firstDigit < firstLetter) {
    massNumber = std::atoi(input.substr(firstDigit, firstLetter - firstDigit).c_str());
    symbol = input.substr(firstLetter);
  } else {
    massNumber = std::atoi(input.substr(firstDigit).c_str());
    symbol = input.substr(firstLetter, firstDigit - firstLetter);
  }

  symbol = CanonicalSymbol(symbol);
  return massNumber > 0 && !symbol.empty();
}

std::string GNucleus::SortName(const char* name) {
  int massNumber = 0;
  std::string symbol;

  if(!ParseName(name, massNumber, symbol))
    return "";

  return Form("%d%s", massNumber, symbol.c_str());
}

void GNucleus::SetNameFromSymbol() {
  if(fSymbol.empty() || GetA() <= 0) {
    TNamed::SetName("");
    return;
  }

  TNamed::SetName(Form("%s%d", fSymbol.c_str(), GetA()));
}

void GNucleus::SetZ(int z) {
  fZ = z;
  fA = fZ + fN;
}

void GNucleus::SetN(int n) {
  fN = n;
  fA = fZ + fN;
}

void GNucleus::SetMassExcess(double massExcess) {
  fMassExcess = massExcess;
}

void GNucleus::SetMass(double mass) {
  fMass = mass;
}

void GNucleus::SetMass() {
  fMass = kAtomicMassUnitMeV * GetA() + GetMassExcess();
}

void GNucleus::SetSymbol(const char* symbol) {
  fSymbol = CanonicalSymbol(symbol ? symbol : "");
  SetNameFromSymbol();
}

bool GNucleus::LoadMassFile(const char* massFile) {
  if(!massFile || std::strlen(massFile) == 0)
    return false;

  std::ifstream input(massFile);
  if(!input.is_open()) {
    std::cout << "Could not open mass file: " << massFile << std::endl;
    return false;
  }

  std::string line;
  while(std::getline(input, line)) {
    line = Trim(line);

    if(line.empty())
      continue;

    if(line.compare(0, 1, "#") == 0 || line.compare(0, 2, "//") == 0)
      continue;

    std::stringstream ss(line);

    int n = 0;
    int z = 0;
    std::string tableName;
    double massExcessKeV = 0;

    if(!(ss >> n >> z >> tableName >> massExcessKeV))
      continue;

    int tableA = 0;
    std::string tableSymbol;
    if(!ParseName(tableName.c_str(), tableA, tableSymbol))
      continue;

    const bool matchesZN = (fZ > 0 || fN > 0) && z == fZ && n == fN;
    const bool matchesName = fA == tableA && EqualsIgnoreCase(fSymbol, tableSymbol);

    if(!matchesZN && !matchesName)
      continue;

    fZ = z;
    fN = n;
    fA = fZ + fN;
    fSymbol = tableSymbol;
    fMassExcess = massExcessKeV / 1000.0;
    SetMass();

    return true;
  }

  std::cout << "Could not find nucleus in mass file: "
            << (fSymbol.empty() ? "" : fSymbol.c_str())
            << fA
            << std::endl;

  return false;
}

bool GNucleus::LoadTransitionFile() {
  if(fTransitionList.GetSize() > 0)
    return false;

  const std::string sourceDir = SourceDataDirectory();
  if(sourceDir.empty())
    return false;

  if(fSymbol.empty() || GetA() <= 0)
    return false;

  std::string filename = sourceDir + "/";
  filename += ToLower(fSymbol);
  filename += Form("%d.sou", GetA());

  std::ifstream input(filename.c_str());
  if(!input.is_open()) {
    return false;
  }

  std::string line;
  while(std::getline(input, line)) {
    line = Trim(line);

    if(line.empty())
      continue;

    if(line.compare(0, 1, "#") == 0 || line.compare(0, 2, "//") == 0)
      continue;

    std::stringstream ss(line);

    double energy = 0;
    double energyUncertainty = 0;
    double intensity = 0;
    double intensityUncertainty = 0;

    if(!(ss >> energy))
      continue;

    ss >> energyUncertainty >> intensity >> intensityUncertainty;

    AddTransition(energy, intensity, energyUncertainty, intensityUncertainty);
  }

  return true;
}

void GNucleus::AddTransition(double energy,
                             double intensity,
                             double energyUncertainty,
                             double intensityUncertainty) {
  auto* transition = new GTransition;
  transition->SetEnergy(energy);
  transition->SetEnergyUncertainty(energyUncertainty);
  transition->SetIntensity(intensity);
  transition->SetIntensityUncertainty(intensityUncertainty);

  AddTransition(transition);
}

void GNucleus::AddTransition(GTransition* transition) {
  if(!transition)
    return;

  fTransitionList.Add(transition);
}

GTransition* GNucleus::GetTransition(int index) const {
  auto* transition = static_cast<GTransition*>(fTransitionList.At(index));
  if(!transition)
    std::cout << "Transition index out of range: " << index << std::endl;

  return transition;
}

double GNucleus::GetRadius() const {
  if(GetA() <= 0)
    return 0;

  return 1.12 * std::pow(GetA(), 1.0 / 3.0)
       - 0.94 * std::pow(GetA(), -1.0 / 3.0);
}

int GNucleus::GetZfromSymbol(const char* symbol) {
  if(!symbol)
    return 0;

  static const char* symbols[] = {
    "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne",
    "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca",
    "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
    "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"
  };

  const std::string query = CanonicalSymbol(symbol);
  for(size_t i = 0; i < sizeof(symbols) / sizeof(symbols[0]); ++i) {
    if(query == symbols[i]) {
      SetZ(static_cast<int>(i + 1));
      return fZ;
    }
  }

  SetZ(0);
  return 0;
}

double GNucleus::GetEnergyFromBeta(double beta) const {
  if(fMass <= 0 || beta <= 0 || beta >= 1)
    return 0;

  const double gamma = 1.0 / std::sqrt(1.0 - beta * beta);
  return fMass * (gamma - 1.0);
}

double GNucleus::GetBetaFromEnergy(double energyMeV) const {
  if(fMass <= 0 || energyMeV <= 0)
    return 0;

  const double gamma = energyMeV / fMass + 1.0;
  return std::sqrt(1.0 - 1.0 / (gamma * gamma));
}

void GNucleus::Print(Option_t* opt) const {
  std::cout << "Nucleus: " << GetName()
            << " Z=" << GetZ()
            << " N=" << GetN()
            << " A=" << GetA()
            << " mass=" << GetMass()
            << " MeV"
            << std::endl;

  for(int i = 0; i < fTransitionList.GetSize(); ++i) {
    auto* transition = static_cast<GTransition*>(fTransitionList.At(i));
    if(!transition)
      continue;

    std::cout << "\t" << i << "\t";
    transition->Print();
  }
}

void GNucleus::WriteSourceFile(std::string outfilename) const {
  if(outfilename.empty())
    return;

  std::ofstream output(outfilename.c_str());
  if(!output.is_open()) {
    std::cout << "Could not write source file: " << outfilename << std::endl;
    return;
  }

  for(int i = 0; i < fTransitionList.GetSize(); ++i) {
    auto* transition = static_cast<GTransition*>(fTransitionList.At(i));
    if(!transition)
      continue;

    output << transition->PrintToString() << std::endl;
  }
}

