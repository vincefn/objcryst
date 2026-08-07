// Standalone utility to (re)generate the ground-truth fixtures under
// test/data/ground_truth/. Run with the "dump-ground-truth" argument and
// redirect stdout to inspect/update the corresponding fixture file(s).
//
// This intentionally duplicates the same object construction used by the
// ground-truth regression tests (see api_scattering.cpp / api_powderpattern.cpp)
// so that the printed values are always in sync with what those tests expect.
#include <iomanip>
#include <iostream>
#include <string>
#include <algorithm>
#include <vector>

#include "ObjCryst/ObjCryst/Crystal.h"
#include "ObjCryst/ObjCryst/DiffractionDataSingleCrystal.h"
#include "ObjCryst/ObjCryst/PowderPattern.h"
#include "ObjCryst/ObjCryst/ReflectionProfile.h"

#include "test_common.h"

using namespace objcryst_test;

namespace
{

std::vector<unsigned long> SelectPowderSampleIndices(const CrystVector_REAL& calc,
                                                     const unsigned long targetCount)
{
   const long n = calc.numElements();
   const REAL minUsefulIntensity = static_cast<REAL>(1e-8);
   std::vector<unsigned long> peaks;
   peaks.reserve(static_cast<size_t>(n));
   for(long i = 1; i + 1 < n; ++i)
   {
      const REAL c = calc(i);
      if(c <= minUsefulIntensity) continue;
      if(c >= calc(i - 1) && c >= calc(i + 1)) peaks.push_back(static_cast<unsigned long>(i));
   }
   std::sort(peaks.begin(), peaks.end(), [&](const unsigned long a, const unsigned long b)
   {
      return calc(a) > calc(b);
   });

   std::vector<unsigned long> selected;
   selected.reserve(targetCount);
   std::vector<char> used(static_cast<size_t>(n), 0);
   auto addIndex = [&](const long idx)
   {
      if(idx < 0 || idx >= n) return;
      const unsigned long uidx = static_cast<unsigned long>(idx);
      if(used[uidx] != 0) return;
      if(calc(uidx) <= minUsefulIntensity) return;
      used[uidx] = 1;
      selected.push_back(uidx);
   };

   for(const unsigned long peak : peaks)
   {
      for(long delta = -8; delta <= 8; ++delta) addIndex(static_cast<long>(peak) + delta);
      if(selected.size() >= targetCount) break;
   }

   if(selected.size() < targetCount)
   {
      std::vector<unsigned long> nonZero;
      nonZero.reserve(static_cast<size_t>(n));
      for(long i = 0; i < n; ++i)
         if(calc(i) > minUsefulIntensity) nonZero.push_back(static_cast<unsigned long>(i));
      std::sort(nonZero.begin(), nonZero.end(), [&](const unsigned long a, const unsigned long b)
      {
         return calc(a) > calc(b);
      });
      for(const unsigned long idx : nonZero)
      {
         addIndex(static_cast<long>(idx));
         if(selected.size() >= targetCount) break;
      }
   }

   std::sort(selected.begin(), selected.end());
   return selected;
}

void DumpPowderGroundTruth(const char* const label, const ObjCryst::PowderPattern& pattern)
{
   const auto& calc = pattern.GetPowderPatternCalc();
   const std::vector<unsigned long> indices = SelectPowderSampleIndices(calc, 80);
   std::cout << label << "\n";
   for(const unsigned long idx : indices) std::cout << idx << " " << calc(idx) << "\n";
}

void DumpGroundTruthData()
{
   using namespace ObjCryst;
   std::cout << std::setprecision(9) << std::fixed;

   Crystal cx = MakePbso4Crystal();
   DiffractionDataSingleCrystal scX(cx, false);
   scX.SetRadiationType(RAD_XRAY);
   scX.SetWavelength(1.54056f);
   scX.GenHKLFullSpace2(0.32, true);
   std::cout << "singlecrystal-xray\n";
   for(int i = 0; i < 20; ++i)
   {
      std::cout << static_cast<long>(scX.GetH()(i)) << " "
                << static_cast<long>(scX.GetK()(i)) << " "
                << static_cast<long>(scX.GetL()(i)) << " "
                << scX.GetFhklCalcSq()(i) << "\n";
   }

   Crystal cn = MakePbso4Crystal();
   DiffractionDataSingleCrystal scN(cn, false);
   scN.SetRadiationType(RAD_NEUTRON);
   scN.SetWavelength(1.54056f);
   scN.GenHKLFullSpace2(0.32, true);
   std::cout << "singlecrystal-neutron\n";
   for(int i = 0; i < 20; ++i)
   {
      std::cout << static_cast<long>(scN.GetH()(i)) << " "
                << static_cast<long>(scN.GetK()(i)) << " "
                << static_cast<long>(scN.GetL()(i)) << " "
                << scN.GetFhklCalcSq()(i) << "\n";
   }

   Crystal cXg = MakePbso4Crystal();
   PowderPattern pXg;
   pXg.SetRadiationType(RAD_XRAY);
   pXg.SetWavelength(1.54056f);
   pXg.SetPowderPatternPar(10 * DEG2RAD, 0.05 * DEG2RAD, 240);
   CrystVector_REAL obsXg(240); obsXg = 1; pXg.SetPowderPatternObs(obsXg);
   auto* phaseXg = new PowderPatternDiffraction;
   phaseXg->SetCrystal(cXg);
   auto* pvGx = new ReflectionProfilePseudoVoigt();
   pvGx->SetProfilePar(.03f * DEG2RAD * DEG2RAD, 0, 0, 0.0f, 0);
   phaseXg->SetProfile(pvGx);
   pXg.AddPowderPatternComponent(*phaseXg);
   pXg.Prepare();
   DumpPowderGroundTruth("powder-xray-pseudo-voigt-gaussian", pXg);

   Crystal cXl = MakePbso4Crystal();
   PowderPattern pXl;
   pXl.SetRadiationType(RAD_XRAY);
   pXl.SetWavelength(1.54056f);
   pXl.SetPowderPatternPar(10 * DEG2RAD, 0.05 * DEG2RAD, 240);
   CrystVector_REAL obsXl(240); obsXl = 1; pXl.SetPowderPatternObs(obsXl);
   auto* phaseXl = new PowderPatternDiffraction;
   phaseXl->SetCrystal(cXl);
   auto* pvLx = new ReflectionProfilePseudoVoigt();
   pvLx->SetProfilePar(.03f * DEG2RAD * DEG2RAD, 0, 0, 1.0f, 0);
   phaseXl->SetProfile(pvLx);
   pXl.AddPowderPatternComponent(*phaseXl);
   pXl.Prepare();
   DumpPowderGroundTruth("powder-xray-pseudo-voigt-lorentzian", pXl);

   Crystal cXa = MakePbso4Crystal();
   PowderPattern pXa;
   pXa.SetRadiationType(RAD_XRAY);
   pXa.SetWavelength(1.54056f);
   pXa.SetPowderPatternPar(10 * DEG2RAD, 0.05 * DEG2RAD, 240);
   CrystVector_REAL obsXa(240); obsXa = 1; pXa.SetPowderPatternObs(obsXa);
   auto* phaseXa = new PowderPatternDiffraction;
   phaseXa->SetCrystal(cXa);
   auto* anisoX = new ReflectionProfilePseudoVoigtAnisotropic(MakePbso4AnisotropicProfile());
   phaseXa->SetProfile(anisoX);
   pXa.AddPowderPatternComponent(*phaseXa);
   pXa.Prepare();
   DumpPowderGroundTruth("powder-xray-anisotropic", pXa);

   Crystal cNg = MakePbso4Crystal();
   PowderPattern pNg;
   pNg.SetRadiationType(RAD_NEUTRON);
   pNg.SetWavelength(1.54056f);
   pNg.SetPowderPatternPar(10 * DEG2RAD, 0.05 * DEG2RAD, 240);
   CrystVector_REAL obsNg(240); obsNg = 1; pNg.SetPowderPatternObs(obsNg);
   auto* phaseNg = new PowderPatternDiffraction;
   phaseNg->SetCrystal(cNg);
   auto* pvGn = new ReflectionProfilePseudoVoigt();
   pvGn->SetProfilePar(.03f * DEG2RAD * DEG2RAD, 0, 0, 0.0f, 0);
   phaseNg->SetProfile(pvGn);
   pNg.AddPowderPatternComponent(*phaseNg);
   pNg.Prepare();
   DumpPowderGroundTruth("powder-neutron-pseudo-voigt-gaussian", pNg);

   Crystal cNl = MakePbso4Crystal();
   PowderPattern pNl;
   pNl.SetRadiationType(RAD_NEUTRON);
   pNl.SetWavelength(1.54056f);
   pNl.SetPowderPatternPar(10 * DEG2RAD, 0.05 * DEG2RAD, 240);
   CrystVector_REAL obsNl(240); obsNl = 1; pNl.SetPowderPatternObs(obsNl);
   auto* phaseNl = new PowderPatternDiffraction;
   phaseNl->SetCrystal(cNl);
   auto* pvLn = new ReflectionProfilePseudoVoigt();
   pvLn->SetProfilePar(.03f * DEG2RAD * DEG2RAD, 0, 0, 1.0f, 0);
   phaseNl->SetProfile(pvLn);
   pNl.AddPowderPatternComponent(*phaseNl);
   pNl.Prepare();
   DumpPowderGroundTruth("powder-neutron-pseudo-voigt-lorentzian", pNl);
}

} // namespace

int main(int argc, char* argv[])
{
   if(argc != 2)
   {
      std::cerr << "Usage: ground_truth <test-case>" << std::endl;
      return 2;
   }

   const std::string testName = argv[1];

   if(testName == "dump-ground-truth") DumpGroundTruthData();
   else
   {
      std::cerr << "Unknown test case: " << testName << std::endl;
      return 2;
   }

   return 0;
}
