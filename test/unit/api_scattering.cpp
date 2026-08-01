// Extended unit tests for the ScatteringData and DiffractionDataSingleCrystal classes.
// Ground-truth regression checks are included here for radiation types
// exercised via single-crystal simulation (X-ray, neutron).
#include <iostream>
#include <string>

#include "ObjCryst/ObjCryst/Crystal.h"
#include "ObjCryst/ObjCryst/DiffractionDataSingleCrystal.h"

#include "test_common.h"

using namespace objcryst_test;

namespace
{

void TestScatteringDataAndSingleCrystal()
{
   using namespace ObjCryst;
   Crystal c = MakePbso4Crystal();
   DiffractionDataSingleCrystal sc(c, false);
   sc.SetRadiationType(RAD_XRAY);
   sc.SetWavelength(1.54056);
   sc.GenHKLFullSpace2(0.35, true);
   Check(sc.GetNbRefl() > 50, "Single-crystal reflection generation failed");
   Check(sc.GetH().numElements() == sc.GetNbRefl(), "H array size mismatch");
   Check(sc.GetK().numElements() == sc.GetNbRefl(), "K array size mismatch");
   Check(sc.GetL().numElements() == sc.GetNbRefl(), "L array size mismatch");
   Check(sc.GetSinThetaOverLambda().numElements() == sc.GetNbRefl(), "sin(theta)/lambda array size mismatch");
   Check(sc.GetTheta().numElements() == sc.GetNbRefl(), "theta array size mismatch");
   Check(sc.GetFhklCalcSq().max() > 0, "Calculated |F|^2 should have positive values");
}

void TestScatteringDataRadiationTypes()
{
   using namespace ObjCryst;
   // X-rays: form factors decrease with sin(theta)/lambda, so |F(000)| should
   // exceed |F| at higher scattering angle for the same reflection family.
   Crystal cx = MakePbso4Crystal();
   DiffractionDataSingleCrystal scX(cx, false);
   scX.SetRadiationType(RAD_XRAY);
   scX.SetWavelength(1.54056);
   scX.GenHKLFullSpace2(0.35, true);
   Check(scX.GetRadiationType() == RAD_XRAY, "Radiation type should be RAD_XRAY");

   // Neutrons: scattering lengths are essentially independent of scattering angle,
   // and can be negative for some elements - just check the simulation runs and
   // produces finite, non-trivial intensities.
   Crystal cn = MakePbso4Crystal();
   DiffractionDataSingleCrystal scN(cn, false);
   scN.SetRadiationType(RAD_NEUTRON);
   scN.SetWavelength(1.54056);
   scN.GenHKLFullSpace2(0.35, true);
   Check(scN.GetRadiationType() == RAD_NEUTRON, "Radiation type should be RAD_NEUTRON");
   Check(scN.GetFhklCalcSq().max() > 0, "Neutron |F|^2 should have positive values");

   // Electrons: exercised mainly to confirm the radiation type is accepted end-to-end.
   Crystal ce = MakePbso4Crystal();
   DiffractionDataSingleCrystal scE(ce, false);
   scE.SetRadiationType(RAD_ELECTRON);
   scE.SetWavelength(0.025);
   scE.GenHKLFullSpace2(0.35, true);
   Check(scE.GetRadiationType() == RAD_ELECTRON, "Radiation type should be RAD_ELECTRON");
   Check(scE.GetNbRefl() > 0, "Electron diffraction reflection generation failed");
}

void TestDiffractionDataSingleCrystalObservedData()
{
   using namespace ObjCryst;
   Crystal c = MakePbso4Crystal();
   DiffractionDataSingleCrystal sc(c, false);
   sc.SetRadiationType(RAD_XRAY);
   sc.SetWavelength(1.54056);
   sc.GenHKLFullSpace2(0.35, true);

   sc.SetIobsToIcalc();
   Check(sc.GetIobs().numElements() == sc.GetNbRefl(), "Iobs array size mismatch after SetIobsToIcalc");
   CheckNearAbsRel(sc.GetIobs()(0), sc.GetIcalc()(0), 1e-6, 1e-8,
                   "SetIobsToIcalc should make Iobs match Icalc");

   CrystVector_REAL sigma(sc.GetNbRefl());
   sigma = 1.0;
   sc.SetSigma(sigma);
   Check(sc.GetSigma().numElements() == sc.GetNbRefl(), "Sigma array size mismatch");

   CrystVector_REAL weight(sc.GetNbRefl());
   weight = 1.0;
   sc.SetWeight(weight);
   Check(sc.GetWeight().numElements() == sc.GetNbRefl(), "Weight array size mismatch");
}

void TestSingleCrystalGroundTruthXray()
{
   CompareSingleCrystalSimulationToGroundTruth(ObjCryst::RAD_XRAY, 1.54056f,
                                               "../../test/data/ground_truth/singlecrystal_xray_pbso4.txt", 1e-2f, 1e-5f);
}

void TestSingleCrystalGroundTruthNeutron()
{
   CompareSingleCrystalSimulationToGroundTruth(ObjCryst::RAD_NEUTRON, 1.54056f,
                                               "../../test/data/ground_truth/singlecrystal_neutron_pbso4.txt", 1e-2f, 1e-5f);
}

// Regression test for the crash in WXDiffractionSingleCrystal::OnMenuShowGraph
// when called after OnMenuSimulate.
//
// OnMenuSimulate calls GenHKLFullSpace(), which generates reflections and
// populates Icalc, but leaves Iobs empty.  OnMenuShowGraph then calls
// CrystUpdate(), which iterates over Icalc.numElements() entries and reads
// Iobs(i) at each step — a fatal out-of-bounds access when Iobs is empty.
//
// The fix is for GenHKLFullSpace (or a post-simulate step) to initialise Iobs
// to a consistent size.  This test captures that contract: after
// GenHKLFullSpace, GetIcalc() and GetIobs() must have the same number of
// elements so that graph-display code can safely read both arrays in lockstep.
void TestSingleCrystalSimulateThenShowGraph()
{
   using namespace ObjCryst;
   Crystal c = MakePbso4Crystal();
   DiffractionDataSingleCrystal sc(c, false);
   sc.SetRadiationType(RAD_XRAY);
   sc.SetWavelength(1.54056);

   // Reproduce the OnMenuSimulate action: generate reflections up to 50 degrees
   // theta (all reflections, no Friedel merging), then check that there are
   // reflections (otherwise the test is vacuous).
   sc.GenHKLFullSpace(50.0 * DEG2RAD, false);
   Check(sc.GetNbRefl() > 0, "GenHKLFullSpace produced no reflections");

   // GetIcalc() triggers CalcIcalc() and returns a non-empty array.
   const long nbCalc = sc.GetIcalc().numElements();
   Check(nbCalc > 0, "GetIcalc() is empty after GenHKLFullSpace");

   // GetIobs() must return an array of the same size so that code iterating
   // over both in lockstep (as CrystUpdate does when building the graph) does
   // not read out of bounds.
   const long nbObs = sc.GetIobs().numElements();
   Check(nbObs == nbCalc,
         "GetIobs() size differs from GetIcalc() size after GenHKLFullSpace: "
         "crash would occur in WXDiffractionSingleCrystal::CrystUpdate");
}

} // namespace

int main(int argc, char* argv[])
{
   if(argc != 2)
   {
      std::cerr << "Usage: api_scattering <test-case>" << std::endl;
      return 2;
   }

   const std::string testName = argv[1];

   if(testName == "scatteringdata-singlecrystal") TestScatteringDataAndSingleCrystal();
   else if(testName == "scatteringdata-radiation-types") TestScatteringDataRadiationTypes();
   else if(testName == "diffractiondata-observed") TestDiffractionDataSingleCrystalObservedData();
   else if(testName == "singlecrystal-groundtruth-xray") TestSingleCrystalGroundTruthXray();
   else if(testName == "singlecrystal-groundtruth-neutron") TestSingleCrystalGroundTruthNeutron();
   else if(testName == "singlecrystal-simulate-show-graph") TestSingleCrystalSimulateThenShowGraph();
   else
   {
      std::cerr << "Unknown test case: " << testName << std::endl;
      return 2;
   }

   return 0;
}
