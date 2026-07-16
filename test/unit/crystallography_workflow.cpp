#include <cstdlib>
#include <fstream>
#include <iostream>

#include "ObjCryst/ObjCryst/Atom.h"
#include "ObjCryst/ObjCryst/CIF.h"
#include "ObjCryst/ObjCryst/Crystal.h"
#include "ObjCryst/ObjCryst/DiffractionDataSingleCrystal.h"
#include "ObjCryst/ObjCryst/Indexing.h"
#include "ObjCryst/ObjCryst/PowderPattern.h"

namespace
{
void Check(const bool ok, const char* const message)
{
   if(!ok)
   {
      std::cerr << message << std::endl;
      std::exit(1);
   }
}
}

int main()
{
   using namespace ObjCryst;

   // Crystal structure creation
   Crystal crystal(8.516, 5.399, 6.989, "Pnma");
   auto* pb = new ScatteringPowerAtom("Pb", "Pb", 1.48);
   auto* s = new ScatteringPowerAtom("S", "S", 0.74);
   crystal.AddScatteringPower(pb);
   crystal.AddScatteringPower(s);
   crystal.AddScatterer(new Atom(.1882, .250, .167, "Pb1", pb, 1.));
   crystal.AddScatterer(new Atom(.437, .750, .186, "S1", s, 1.));
   Check(crystal.GetNbScatterer() == 2, "Crystal creation test failed");

   // Single crystal diffraction
   DiffractionDataSingleCrystal scDiff(crystal, false);
   scDiff.SetRadiationType(RAD_XRAY);
   scDiff.SetWavelength(1.54056);
   scDiff.GenHKLFullSpace2(0.35, true);
   Check(scDiff.GetNbRefl() > 50, "Single-crystal reflection generation failed");
   const CrystVector_REAL& fcalc = scDiff.GetFhklCalcSq();
   Check(fcalc.numElements() == scDiff.GetNbRefl(), "Single-crystal |F|^2 array size mismatch");
   Check(fcalc.max() > 0, "Single-crystal |F|^2 is empty");

   // Powder diffraction calculation
   PowderPattern powderCalc;
   powderCalc.SetRadiationType(RAD_XRAY);
   powderCalc.SetWavelength("CuA1");
   powderCalc.SetPowderPatternPar(0, 0.01 * DEG2RAD, 2000);
   CrystVector_REAL obs(2000);
   obs = 1;
   powderCalc.SetPowderPatternObs(obs);
   auto* phase = new PowderPatternDiffraction;
   phase->SetCrystal(crystal);
   phase->SetReflectionProfilePar(PROFILE_PSEUDO_VOIGT, .03 * DEG2RAD * DEG2RAD, 0, 0, 0.3, 0);
   powderCalc.AddPowderPatternComponent(*phase);
   powderCalc.Prepare();
   const CrystVector_REAL& pcalc = powderCalc.GetPowderPatternCalc();
   Check(pcalc.numElements() == powderCalc.GetNbPoint(), "Powder pattern calculation size mismatch");
   Check(pcalc.max() > 0, "Powder diffraction calculation failed");

   // Powder pattern data import
   PowderPattern powderImported;
   powderImported.SetRadiationType(RAD_XRAY);
   powderImported.SetWavelength("CuA1");
   powderImported.ImportPowderPatternFullprof("../../ObjCryst/example/pbso4/xray-pattern.dat");
   Check(powderImported.GetNbPoint() > 100, "Powder import did not load enough points");
   const CrystVector_REAL& iobs = powderImported.GetPowderPatternObs();
   Check(iobs.numElements() == powderImported.GetNbPoint(), "Imported powder data size mismatch");

   // CIF parsing and crystal creation from COD data
   std::ifstream cifIn("../../test/data/cif/PbSO4-COD-1528837.cif");
   Check(cifIn.good(), "Unable to open COD CIF test file");
   CIF cif(cifIn, true, false);
   Crystal cifCrystal;
   Crystal* cifResult = CreateCrystalFromCIF(cif, false, true, false, false, &cifCrystal);
   Check(cifResult == &cifCrystal, "CIF crystal creation returned unexpected pointer");
   Check(cifCrystal.GetSpaceGroup().GetSpaceGroupNumber() == 62, "CIF spacegroup mismatch");
   Check(cifCrystal.GetNbScatterer() >= 5, "CIF crystal has too few scatterers");

   // Indexing workflow
   PeakList peakList;
   const float volume = peakList.Simulate(0, 8.516f, 5.399f, 6.989f, 90.f, 90.f, 90.f, true, 15, 0, 0, 0, false);
   Check(volume > 250.f, "PeakList simulation produced invalid cell volume");
   const float estimatedVolume = EstimateCellVolume(peakList.GetPeakList().back().dobs,
                                                    peakList.GetPeakList().front().dobs / 10.f,
                                                    static_cast<float>(peakList.GetPeakList().size()),
                                                    ORTHOROMBIC, LATTICE_P, 1.f);
   Check(estimatedVolume > 0, "Cell volume estimation failed");

   RecUnitCell cell(0.f, 1.f / 8.516f, 1.f / 5.399f, 1.f / 6.989f, 0.f, 0.f, 0.f,
                    ORTHOROMBIC, LATTICE_P, 0);
   const float indexingScore = Score(peakList, cell, 0, false, false, false);
   Check(indexingScore > 0, "Indexing score computation failed");

   return 0;
}
