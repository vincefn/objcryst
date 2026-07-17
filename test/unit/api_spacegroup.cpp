// Extended unit tests for the ObjCryst::SpaceGroup class.
#include <cmath>
#include <iostream>
#include <string>

#include "ObjCryst/ObjCryst/SpaceGroup.h"

#include "test_common.h"

using namespace objcryst_test;

namespace
{

void TestSpaceGroupBasics()
{
   using namespace ObjCryst;
   SpaceGroup sg("Pnma");
   Check(sg.GetSpaceGroupNumber() == 62, "SpaceGroup number mismatch for Pnma");
   Check(sg.IsCentrosymmetric(), "Pnma should be centrosymmetric");
   Check(sg.GetNbTranslationVectors() == 1, "Pnma should have one translation vector");
   Check(sg.GetNbSymmetrics(false, false) >= 4, "Pnma should expose at least 4 symmetrics");
   Check(sg.GetName().empty() == false, "SpaceGroup name should not be empty");
}

void TestSpaceGroupAlternateSettings()
{
   using namespace ObjCryst;
   SpaceGroup p1("P1");
   Check(p1.GetSpaceGroupNumber() == 1, "Unexpected spacegroup number for P1");
   Check(!p1.IsCentrosymmetric(), "P1 should not be centrosymmetric");

   SpaceGroup pminus1("P-1");
   Check(pminus1.GetSpaceGroupNumber() == 2, "Unexpected spacegroup number for P-1");
   Check(pminus1.IsCentrosymmetric(), "P-1 should be centrosymmetric");

   SpaceGroup fm3m("F m -3 m");
   Check(fm3m.GetSpaceGroupNumber() == 225, "Unexpected spacegroup number for Fm-3m");
   Check(fm3m.GetNbTranslationVectors() == 4, "Fm-3m should have 4 translation vectors (F-centering)");

   SpaceGroup c2c("C 2/c");
   Check(c2c.GetSpaceGroupNumber() == 15, "Unexpected spacegroup number for C2/c");
   Check(c2c.GetUniqueAxis() <= 2, "C2/c unique axis should be a,b or c (index 0..2)");
}

void TestSpaceGroupReflectionProperties()
{
   using namespace ObjCryst;
   // Pnma systematic absences: 0kl absent unless k+l=2n ; hk0 absent unless h=2n ; h00 absent unless h=2n
   SpaceGroup sg("Pnma");
   Check(sg.IsReflSystematicAbsent(1, 0, 0), "Pnma h00 with h odd should be systematically absent");
   Check(!sg.IsReflSystematicAbsent(2, 0, 0), "Pnma h00 with h even should not be systematically absent");
   Check(!sg.IsReflSystematicAbsent(1, 1, 1), "Pnma general reflection should not be systematically absent");

   const unsigned int equiv = sg.AreReflEquiv(1, 1, 1, -1, -1, -1);
   Check(equiv != 0, "Pnma (1,1,1) and (-1,-1,-1) should be equivalent (centrosymmetric)");

   const unsigned int factor = sg.GetExpectedIntensityFactor(1, 1, 1);
   Check(factor >= 1, "Expected intensity factor should be at least 1");

   CrystMatrix_REAL equivRefl = sg.GetAllEquivRefl(1, 1, 1);
   Check(equivRefl.rows() >= 1, "GetAllEquivRefl should return at least the input reflection");
   Check(equivRefl.cols() == 5, "GetAllEquivRefl should return 5 columns (h,k,l,Re,Im)");
}

void TestSpaceGroupSymmetryOperations()
{
   using namespace ObjCryst;
   SpaceGroup sg("Pnma");
   const std::vector<SpaceGroup::SMx>& ops = sg.GetSymmetryOperations();
   Check(!ops.empty(), "Pnma should expose symmetry operations");
   const std::vector<SpaceGroup::TRx>& trans = sg.GetTranslationVectors();
   Check(trans.size() == static_cast<size_t>(sg.GetNbTranslationVectors()), "Translation vector count mismatch");
   Check(trans[0].tr[0] == 0 && trans[0].tr[1] == 0 && trans[0].tr[2] == 0,
        "First translation vector should be (0,0,0)");

   CrystMatrix_REAL symmetrics = sg.GetAllSymmetrics(0.1, 0.2, 0.3);
   Check(symmetrics.rows() == static_cast<long>(sg.GetNbSymmetrics(false, false)),
        "GetAllSymmetrics row count should match GetNbSymmetrics");
   Check(symmetrics.cols() == 3, "GetAllSymmetrics should return 3 columns (x,y,z)");

   REAL x = 0.1, y = 0.2, z = 0.3;
   sg.GetSymmetric(0, x, y, z);
   // GetSymmetric(i,...) is not guaranteed to be the identity for i=0, but the
   // result must be one of the positions already reported by GetAllSymmetrics().
   bool foundMatch = false;
   for(long i = 0; i < symmetrics.rows(); ++i)
   {
      if(std::fabs(symmetrics(i, 0) - x) < 1e-5
         && std::fabs(symmetrics(i, 1) - y) < 1e-5
         && std::fabs(symmetrics(i, 2) - z) < 1e-5)
      {
         foundMatch = true;
         break;
      }
   }
   Check(foundMatch, "GetSymmetric(0,...) result should be among the positions from GetAllSymmetrics()");
}

void TestSpaceGroupAsymmetricUnit()
{
   using namespace ObjCryst;
   SpaceGroup sg("Pnma");
   const AsymmetricUnit& asu = sg.GetAsymUnit();
   Check(asu.Xmax() > 0 && asu.Xmax() <= 1.0, "AsymmetricUnit Xmax should be within (0,1]");
   Check(asu.Ymax() > 0 && asu.Ymax() <= 1.0, "AsymmetricUnit Ymax should be within (0,1]");
   Check(asu.Zmax() > 0 && asu.Zmax() <= 1.0, "AsymmetricUnit Zmax should be within (0,1]");

   // SpaceGroup::IsInAsymmetricUnit() simply forwards to AsymmetricUnit::IsInAsymmetricUnit();
   // verify that forwarding is consistent for a handful of sample positions.
   const REAL testX[] = {0.0, 0.1, 0.3, 0.5};
   for(int i = 0; i < 4; ++i)
   {
      const bool viaAsu = asu.IsInAsymmetricUnit(testX[i], testX[i], testX[i]);
      const bool viaSg = sg.IsInAsymmetricUnit(testX[i], testX[i], testX[i]);
      Check(viaAsu == viaSg, "SpaceGroup::IsInAsymmetricUnit should forward to AsymmetricUnit::IsInAsymmetricUnit");
   }
}

} // namespace

int main(int argc, char* argv[])
{
   if(argc != 2)
   {
      std::cerr << "Usage: api_spacegroup <test-case>" << std::endl;
      return 2;
   }

   const std::string testName = argv[1];

   if(testName == "spacegroup") TestSpaceGroupBasics();
   else if(testName == "spacegroup-alternate-settings") TestSpaceGroupAlternateSettings();
   else if(testName == "spacegroup-reflection-properties") TestSpaceGroupReflectionProperties();
   else if(testName == "spacegroup-symmetry-operations") TestSpaceGroupSymmetryOperations();
   else if(testName == "spacegroup-asymmetric-unit") TestSpaceGroupAsymmetricUnit();
   else
   {
      std::cerr << "Unknown test case: " << testName << std::endl;
      return 2;
   }

   return 0;
}
