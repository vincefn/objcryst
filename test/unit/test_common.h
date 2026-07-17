// Shared helpers for the ObjCryst unit test suite.
//
// This header is included by several independent test binaries
// (api_spacegroup, api_crystal, api_molecule, api_scattering,
// api_powderpattern, api_cif, api_optimization, api_indexing,
// ground_truth). Every function here is declared 'inline' since the
// header may be included by multiple translation units that are each
// linked into their own standalone executable.
#ifndef OBJCRYST_TEST_COMMON_H
#define OBJCRYST_TEST_COMMON_H

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "ObjCryst/ObjCryst/Atom.h"
#include "ObjCryst/ObjCryst/Crystal.h"
#include "ObjCryst/ObjCryst/DiffractionDataSingleCrystal.h"
#include "ObjCryst/ObjCryst/PowderPattern.h"
#include "ObjCryst/ObjCryst/ReflectionProfile.h"
#include "ObjCryst/ObjCryst/ScatteringPower.h"

namespace objcryst_test
{

inline void Check(const bool ok, const char* const msg)
{
   if(!ok)
   {
      std::cerr << msg << std::endl;
      std::exit(1);
   }
}

inline void CheckNearAbs(const REAL actual, const REAL expected, const REAL tol, const std::string& msg)
{
   if(std::fabs(actual - expected) > tol)
   {
      std::ostringstream os;
      os << msg << " (actual=" << actual << ", expected=" << expected << ", tol=" << tol << ")";
      std::cerr << os.str() << std::endl;
      std::exit(1);
   }
}

inline void CheckNearAbsRel(const REAL actual, const REAL expected, const REAL absTol, const REAL relTol,
                            const std::string& msg)
{
   const REAL diff = std::fabs(actual - expected);
   const REAL scale = std::max<REAL>(1.0f, std::fabs(expected));
   const REAL tol = std::max(absTol, relTol * scale);
   if(diff > tol)
   {
      std::ostringstream os;
      os << msg << " (actual=" << actual << ", expected=" << expected
         << ", absTol=" << absTol << ", relTol=" << relTol << ")";
      std::cerr << os.str() << std::endl;
      std::exit(1);
   }
}

inline bool IsTraceEnabled()
{
   const char* const trace = std::getenv("OBJCRYST_TEST_TRACE");
   if(nullptr == trace || '\0' == trace[0]) return false;
   const std::string value(trace);
   return value != "0" && value != "false" && value != "FALSE";
}

inline void Trace(const std::string& message)
{
   if(IsTraceEnabled()) std::cout << "[trace] " << message << std::endl;
}

/// Build the reference PbSO4 crystal structure (Pnma) shared by most tests.
inline ObjCryst::Crystal MakePbso4Crystal()
{
   using namespace ObjCryst;
   Crystal c(8.516, 5.399, 6.989, "Pnma");
   auto* pb = new ScatteringPowerAtom("Pb", "Pb", 1.48);
   auto* s = new ScatteringPowerAtom("S", "S", 0.74);
   auto* o = new ScatteringPowerAtom("O", "O", 1.2);
   c.AddScatteringPower(pb);
   c.AddScatteringPower(s);
   c.AddScatteringPower(o);
   c.AddScatterer(new Atom(.1882, .250, .167, "Pb1", pb, 1.));
   c.AddScatterer(new Atom(.4370, .750, .186, "S1", s, 1.));
   c.AddScatterer(new Atom(.5950, .750, .100, "O1", o, 1.));
   c.AddScatterer(new Atom(.3190, .750, .043, "O2", o, 1.));
   c.AddScatterer(new Atom(.4150, .974, .306, "O3", o, 1.));
   return c;
}

/// Build the reference anisotropic pseudo-Voigt profile used for the
/// PbSO4 ground-truth tests (asymA0=0.9 keeps the asymmetry factor away
/// from the degenerate asym=0 case - see ReflectionProfile.cpp clamp).
inline ObjCryst::ReflectionProfilePseudoVoigtAnisotropic MakePbso4AnisotropicProfile()
{
   using namespace ObjCryst;
   ReflectionProfilePseudoVoigtAnisotropic profile;
   profile.SetProfilePar(.02f * DEG2RAD * DEG2RAD, .003f * DEG2RAD * DEG2RAD, 0, 0,
                         0.002f * DEG2RAD, 0.002f * DEG2RAD, 0, 0, 0, 0, 0, 0, 0.9f, 0, 0, 0, 0);
   return profile;
}

inline CrystVector_REAL MakePowderPatternXGrid(const REAL xmin, const REAL step, const unsigned long nbPoint)
{
   CrystVector_REAL x(nbPoint);
   for(unsigned long i = 0; i < nbPoint; ++i) x(i) = xmin + step * static_cast<REAL>(i);
   return x;
}

inline REAL StolToTwoTheta(const REAL stol, const REAL wavelength)
{
   const REAL x = stol * wavelength;
   if(std::fabs(x) >= 1.0) return 2 * M_PI;
   return 2 * asin(x);
}

struct SingleCrystalPoint
{
   long h;
   long k;
   long l;
   REAL f2;
};

struct PowderPoint
{
   unsigned long idx;
   REAL intensity;
};

inline std::vector<SingleCrystalPoint> LoadSingleCrystalGroundTruth(const std::string& path)
{
   std::ifstream in(path.c_str());
   Check(in.good(), ("Unable to open single-crystal ground truth file: " + path).c_str());

   std::vector<SingleCrystalPoint> out;
   std::string line;
   while(std::getline(in, line))
   {
      if(line.empty() || line[0] == '#') continue;
      std::istringstream is(line);
      SingleCrystalPoint p{};
      is >> p.h >> p.k >> p.l >> p.f2;
      Check(!is.fail(), ("Invalid line in single-crystal ground truth file: " + line).c_str());
      out.push_back(p);
   }
   Check(!out.empty(), ("No data in single-crystal ground truth file: " + path).c_str());
   return out;
}

inline std::vector<PowderPoint> LoadPowderGroundTruth(const std::string& path)
{
   std::ifstream in(path.c_str());
   Check(in.good(), ("Unable to open powder ground truth file: " + path).c_str());

   std::vector<PowderPoint> out;
   std::string line;
   while(std::getline(in, line))
   {
      if(line.empty() || line[0] == '#') continue;
      std::istringstream is(line);
      PowderPoint p{};
      is >> p.idx >> p.intensity;
      Check(!is.fail(), ("Invalid line in powder ground truth file: " + line).c_str());
      out.push_back(p);
   }
   Check(!out.empty(), ("No data in powder ground truth file: " + path).c_str());
   return out;
}

inline void CompareSingleCrystalToGroundTruth(const ObjCryst::DiffractionDataSingleCrystal& sc,
                                              const std::string& filePath,
                                              const REAL absTol,
                                              const REAL relTol)
{
   const std::vector<SingleCrystalPoint> gt = LoadSingleCrystalGroundTruth(filePath);
   const auto& h = sc.GetH();
   const auto& k = sc.GetK();
   const auto& l = sc.GetL();
   const auto& f2 = sc.GetFhklCalcSq();

   Check(static_cast<long>(gt.size()) <= sc.GetNbRefl(), "Ground-truth has more entries than generated reflections");

   for(size_t i = 0; i < gt.size(); ++i)
   {
      Check(static_cast<long>(h(i)) == gt[i].h, "Single-crystal H mismatch against ground truth");
      Check(static_cast<long>(k(i)) == gt[i].k, "Single-crystal K mismatch against ground truth");
      Check(static_cast<long>(l(i)) == gt[i].l, "Single-crystal L mismatch against ground truth");
      CheckNearAbsRel(f2(i), gt[i].f2, absTol, relTol, "Single-crystal F^2 mismatch against ground truth");
   }
}

inline void CompareSingleCrystalSimulationToGroundTruth(const ObjCryst::RadiationType radiation,
                                                        const REAL wavelength,
                                                        const std::string& filePath,
                                                        const REAL absTol,
                                                        const REAL relTol)
{
   using namespace ObjCryst;
   Crystal c = MakePbso4Crystal();
   DiffractionDataSingleCrystal sc(c, false);
   sc.SetRadiationType(radiation);
   sc.SetWavelength(wavelength);
   sc.GenHKLFullSpace2(0.32, true);
   CompareSingleCrystalToGroundTruth(sc, filePath, absTol, relTol);
}

inline void ComparePowderSimulationToGroundTruth(const ObjCryst::RadiationType radiation,
                                                 const REAL wavelength,
                                                 ObjCryst::ReflectionProfile* profile,
                                                 const char* const testLabel,
                                                 const std::string& filePath,
                                                 const REAL absTol,
                                                 const REAL relTol)
{
   using namespace ObjCryst;
   Trace(std::string(testLabel) + ": creating crystal");
   Crystal c = MakePbso4Crystal();

   Trace(std::string(testLabel) + ": preparing powder pattern");
   PowderPattern p;
   p.SetRadiationType(radiation);
   p.SetWavelength(wavelength);
   p.SetPowderPatternPar(10 * DEG2RAD, 0.05 * DEG2RAD, 240);
   CrystVector_REAL obs(240);
   obs = 1;
   p.SetPowderPatternObs(obs);

   auto* phase = new PowderPatternDiffraction;
   phase->SetCrystal(c);
   phase->SetProfile(profile);
   p.AddPowderPatternComponent(*phase);
   Trace(std::string(testLabel) + ": running PowderPattern::Prepare()");
   p.Prepare();
   Trace(std::string(testLabel) + ": PowderPattern::Prepare() done");
   Trace(std::string(testLabel) + ": loading powder ground-truth file");

   const std::vector<PowderPoint> gt = LoadPowderGroundTruth(filePath);
   Trace(std::string(testLabel) + ": powder ground-truth loaded");
   Trace(std::string(testLabel) + ": getting calculated powder pattern");
   const auto& calc = p.GetPowderPatternCalc();
   Trace(std::string(testLabel) + ": calculated powder pattern retrieved");
   Trace(std::string(testLabel) + ": getting number of points");
   const unsigned long nbPoint = p.GetNbPoint();
   Trace(std::string(testLabel) + ": number of points retrieved");
   {
      std::ostringstream os;
      os << testLabel << ": validating " << gt.size() << " ground-truth points";
      Trace(os.str());
   }
   for(size_t i = 0; i < gt.size(); ++i)
   {
      const auto& pt = gt[i];
      Check(pt.idx < nbPoint, "Powder ground-truth index out of range");
      const REAL actual = calc(pt.idx);
      CheckNearAbsRel(actual, pt.intensity, absTol, relTol, "Powder intensity mismatch against ground truth");
      if(IsTraceEnabled() && (i == 0 || i + 1 == gt.size() || (i % 10) == 0))
      {
         std::ostringstream os;
         os << testLabel << ": compared point#" << i
            << " (idx=" << pt.idx
            << ", actual=" << actual
            << ", expected=" << pt.intensity << ")";
         Trace(os.str());
      }
   }
   Trace(std::string(testLabel) + ": comparison done");
}

} // namespace objcryst_test

#endif // OBJCRYST_TEST_COMMON_H
