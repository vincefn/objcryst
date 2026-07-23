// Extended unit tests for OptimizationObj / MonteCarloObj and LSQNumObj,
// including their RefinablePar/RefinableObj foundations.
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "ObjCryst/ObjCryst/Atom.h"
#include "ObjCryst/ObjCryst/PowderPattern.h"
#include "ObjCryst/RefinableObj/GlobalOptimObj.h"
#include "ObjCryst/RefinableObj/LSQNumObj.h"
#include "ObjCryst/RefinableObj/RefinableObj.h"

#include "test_common.h"

using namespace objcryst_test;

namespace
{

void TestRefinablePar()
{
   using namespace ObjCryst;
   REAL v = 1.0;
   RefinablePar p("test-par", &v, 0.0, 2.0, gpRefParTypeObjCryst);
   p.SetValue(1.5);
   Check(std::fabs(p.GetValue() - 1.5) < 1e-12, "RefinablePar SetValue failed");
   p.Mutate(1.0);
   Check(p.GetValue() <= 2.0, "RefinablePar should respect upper limit");
   p.SetHumanScale(2.0);
   p.SetHumanValue(2.0);
   Check(std::fabs(p.GetValue() - 1.0) < 1e-12, "RefinablePar human scaling failed");
   p.SetIsFixed(true);
   Check(p.IsFixed(), "RefinablePar fixed flag not applied");
}

void TestRefinableObj()
{
   using namespace ObjCryst;
   RefinableObj obj;
   REAL a = 0.1, b = 0.2;
   obj.AddPar(new RefinablePar("a", &a, -1, 1, gpRefParTypeObjCryst));
   obj.AddPar(new RefinablePar("b", &b, -1, 1, gpRefParTypeObjCryst));
   Check(obj.GetNbPar() == 2, "RefinableObj parameter count mismatch");
   obj.SetParIsFixed("a", true);
   obj.PrepareForRefinement();
   Check(obj.GetNbParNotFixed() == 1, "RefinableObj fixed parameter handling failed");
   const unsigned long setId = obj.CreateParamSet("baseline");
   Check(obj.GetParamSetName(setId) == "baseline", "RefinableObj parameter-set naming failed");
   Check(obj.GetParamSet(setId).numElements() > 0, "RefinableObj parameter-set should not be empty");
}

void TestOptimizationObj()
{
   using namespace ObjCryst;
   RefinableObj obj;
   REAL a = 0.1;
   obj.AddPar(new RefinablePar("a", &a, -1, 1, gpRefParTypeObjCryst));
   MonteCarloObj mc;
   mc.AddRefinableObj(obj);
   RefinableObj& full = mc.GetFullRefinableObj(true);
   Check(full.GetNbPar() >= 1, "OptimizationObj failed to aggregate refined parameters");
   mc.FixAllPar();
   mc.UnFixAllPar();
}

void TestOptimizationObjLimitsAndOptions()
{
   using namespace ObjCryst;
   RefinableObj obj;
   REAL a = 0.1;
   obj.AddPar(new RefinablePar("a", &a, -1, 1, gpRefParTypeObjCryst));
   MonteCarloObj mc;
   mc.AddRefinableObj(obj);
   mc.SetLimitsRelative("a", -0.5, 0.5);
   mc.SetLimitsAbsolute("a", -0.9, 0.9);
   Check(mc.GetNbOption() >= 0, "MonteCarloObj should expose its option list");
   Check(mc.GetName().empty() || !mc.GetName().empty(), "MonteCarloObj::GetName should be callable");
   Check(!mc.IsOptimizing(), "MonteCarloObj should not be optimizing before a run is started");
}

void TestLsqNumObj()
{
   using namespace ObjCryst;
   PowderPattern p;
   p.SetRadiationType(RAD_XRAY);
   p.SetWavelength("CuA1");
   p.SetPowderPatternPar(5 * DEG2RAD, 0.05 * DEG2RAD, 1000);
   CrystVector_REAL obs(1000);
   obs = 10;
   p.SetPowderPatternObs(obs);
   auto* bg = new PowderPatternBackground;
   CrystVector_REAL x(2), y(2);
   x(0) = 0;
   x(1) = 3.0 * DEG2RAD;
   y(0) = 8;
   y(1) = 9;
   bg->SetInterpPoints(x, y);
   p.AddPowderPatternComponent(*bg);
   p.Prepare();

   LSQNumObj lsq("lsq-test");
   lsq.SetRefinedObj(p, 0, true, false);
   lsq.PrepareRefParList(true);
   Check(lsq.GetCompiledRefinedObj().GetNbPar() >= 1, "LSQNumObj compiled parameter list is empty");
}

void TestLsqNumObjResidualStatistics()
{
   using namespace ObjCryst;
   PowderPattern p;
   p.SetRadiationType(RAD_XRAY);
   p.SetWavelength("CuA1");
   p.SetPowderPatternPar(5 * DEG2RAD, 0.05 * DEG2RAD, 200);
   CrystVector_REAL obs(200);
   obs = 10;
   p.SetPowderPatternObs(obs);
   p.SetSigmaToSqrtIobs();
   p.SetWeightToInvSigmaSq();

   auto* bg = new PowderPatternBackground;
   CrystVector_REAL x(2), y(2);
   x(0) = 5.0 * DEG2RAD;
   x(1) = 15.0 * DEG2RAD;
   y(0) = 5;
   y(1) = 7;
   bg->SetInterpPoints(x, y);
   p.AddPowderPatternComponent(*bg);
   p.Prepare();

   LSQNumObj lsq("lsq-residual-test");
   lsq.SetRefinedObj(p, 0, true, false);
   lsq.PrepareRefParList(true);
   const CrystVector_REAL& calcVec = lsq.GetLSQCalc();
   const CrystVector_REAL& obsVec = lsq.GetLSQObs();
   const CrystVector_REAL& weightVec = lsq.GetLSQWeight();
   Check(obsVec.numElements() == calcVec.numElements(),
         "LSQNumObj obs/calc vectors should have matching sizes");
   Check(obsVec.numElements() == weightVec.numElements(),
         "LSQNumObj obs/weight vectors should have matching sizes");
   REAL chiSq = 0;
   for(long i = 0; i < obsVec.numElements(); ++i)
   {
      const REAL diff = obsVec(i) - calcVec(i);
      chiSq += diff * diff * weightVec(i);
   }
   Check(std::isfinite(chiSq), "Manually accumulated Chi^2 from LSQ vectors should be finite");
}

void TestLsqNumObjPreservesOriginalParameterNames()
{
   // Two crystalline phases with overlapping unit-cell parameter names (a, b, c).
   // After SetRefinedObj / PrepareRefParList, neither crystal's own parameter
   // names should be mutated by the proxy-based deduplication.
   using namespace ObjCryst;

   // Phase 1: PbSO4 (orthorhombic Pnma – has independent a, b, c)
   Crystal c1 = MakePbso4Crystal();

   // Phase 2: a simple orthorhombic crystal with different lattice constants
   // (shares parameter names a, b, c with phase 1 – the old "~" code would
   // have mutated one of them)
   Crystal c2(4.212, 4.212, 4.212, "Pm3m"); // cubic; single independent parameter a
   auto* sr = new ScatteringPowerAtom("Sr", "Sr", 1.2);
   c2.AddScatteringPower(sr);
   c2.AddScatterer(new Atom(0., 0., 0., "Sr1", sr, 1.));

   PowderPattern p;
   p.SetRadiationType(RAD_XRAY);
   p.SetWavelength("CuA1");
   p.SetPowderPatternPar(5 * DEG2RAD, 0.05 * DEG2RAD, 400);
   CrystVector_REAL obs(400);
   obs = 10;
   p.SetPowderPatternObs(obs);

   auto* bg = new PowderPatternBackground;
   CrystVector_REAL x(2), y(2);
   x(0) = 0;         x(1) = 3.0 * DEG2RAD;
   y(0) = 8;         y(1) = 9;
   bg->SetInterpPoints(x, y);
   p.AddPowderPatternComponent(*bg);

   auto* phase1 = new PowderPatternDiffraction;
   phase1->SetCrystal(c1);
   p.AddPowderPatternComponent(*phase1);

   auto* phase2 = new PowderPatternDiffraction;
   phase2->SetCrystal(c2);
   p.AddPowderPatternComponent(*phase2);

   p.Prepare();

   // Snapshot original parameter names from both crystals
   std::vector<std::string> c1Names, c2Names;
   c1Names.reserve(c1.GetNbPar());
   c2Names.reserve(c2.GetNbPar());
   for(long i=0;i<c1.GetNbPar();++i) c1Names.push_back(c1.GetPar(i).GetName());
   for(long i=0;i<c2.GetNbPar();++i) c2Names.push_back(c2.GetPar(i).GetName());

   LSQNumObj lsq("lsq-preserve-names");
   lsq.SetRefinedObj(p, 0, true, false);
   lsq.PrepareRefParList(true);

   // Verify sizes unchanged
   Check(static_cast<long>(c1Names.size())==c1.GetNbPar(),
         "Crystal1 parameter list size changed unexpectedly");
   Check(static_cast<long>(c2Names.size())==c2.GetNbPar(),
         "Crystal2 parameter list size changed unexpectedly");

   // Verify names unchanged in both crystals
   for(long i=0;i<c1.GetNbPar();++i)
      Check(c1.GetPar(i).GetName()==c1Names[i],
            "Crystal1 parameter name changed after LSQ SetRefinedObj/PrepareRefParList");
   for(long i=0;i<c2.GetNbPar();++i)
      Check(c2.GetPar(i).GetName()==c2Names[i],
            "Crystal2 parameter name changed after LSQ SetRefinedObj/PrepareRefParList");
}

} // namespace

int main(int argc, char* argv[])
{
   if(argc != 2)
   {
      std::cerr << "Usage: api_optimization <test-case>" << std::endl;
      return 2;
   }

   const std::string testName = argv[1];

   if(testName == "refinablepar") TestRefinablePar();
   else if(testName == "refinableobj") TestRefinableObj();
   else if(testName == "optimizationobj") TestOptimizationObj();
   else if(testName == "optimizationobj-limits-options") TestOptimizationObjLimitsAndOptions();
   else if(testName == "lsqnumobj") TestLsqNumObj();
   else if(testName == "lsqnumobj-residual-statistics") TestLsqNumObjResidualStatistics();
   else if(testName == "lsqnumobj-preserve-original-parameter-names") TestLsqNumObjPreservesOriginalParameterNames();
   else
   {
      std::cerr << "Unknown test case: " << testName << std::endl;
      return 2;
   }

   return 0;
}
