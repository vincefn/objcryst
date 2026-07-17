// Extended unit tests for the ObjCryst::Molecule class.
#include <cmath>
#include <iostream>
#include <string>

#include "ObjCryst/ObjCryst/Crystal.h"
#include "ObjCryst/ObjCryst/Molecule.h"
#include "ObjCryst/ObjCryst/ScatteringPower.h"

#include "test_common.h"

using namespace objcryst_test;

namespace
{

ObjCryst::Molecule MakeWaterLikeMolecule(ObjCryst::Crystal& c)
{
   using namespace ObjCryst;
   auto* spO = new ScatteringPowerAtom("O", "O", 1.0);
   c.AddScatteringPower(spO);
   Molecule m(c, "water-like");
   m.AddAtom(0.0, 0.0, 0.0, spO, "O1", false);
   m.AddAtom(1.0, 0.0, 0.0, spO, "O2", false);
   m.AddAtom(0.0, 1.0, 0.0, spO, "O3", false);
   return m;
}

void TestMoleculeAtomsAndBonds()
{
   using namespace ObjCryst;
   Crystal c(10, 10, 10, "P1");
   Molecule m = MakeWaterLikeMolecule(c);
   Check(m.GetAtomList().size() == 3, "Molecule atom insertion failed");

   m.AddBond(m.GetAtom(0), m.GetAtom(1), 1.0, 0.01, 0.02, 1.0, false);
   Check(m.GetBondList().size() == 1, "Molecule bond insertion failed");
   Check(m.FindBond(m.GetAtom(0), m.GetAtom(1)) != m.GetBondList().end(), "Molecule bond lookup failed");

   const REAL length = GetBondLength(m.GetAtom(0), m.GetAtom(1));
   CheckNearAbs(length, 1.0, 1e-4, "GetBondLength should return the Euclidean distance between atoms");
}

void TestMoleculeAnglesAndDihedrals()
{
   using namespace ObjCryst;
   Crystal c(10, 10, 10, "P1");
   Molecule m = MakeWaterLikeMolecule(c);
   m.AddAtom(1.0, 1.0, 0.0, &m.GetAtom(0).GetScatteringPower(), "O4", false);

   m.AddBond(m.GetAtom(0), m.GetAtom(1), 1.0, 0.01, 0.02, 1.0, false);
   m.AddBond(m.GetAtom(0), m.GetAtom(2), 1.0, 0.01, 0.02, 1.0, false);
   m.AddBondAngle(m.GetAtom(1), m.GetAtom(0), m.GetAtom(2), 90.0, 1.0, 1.0, false);
   Check(m.GetBondAngleList().size() == 1, "Molecule bond-angle insertion failed");
   Check(m.FindBondAngle(m.GetAtom(1), m.GetAtom(0), m.GetAtom(2)) != m.GetBondAngleList().end(),
        "Molecule bond-angle lookup failed");

   const REAL angle = GetBondAngle(m.GetAtom(1), m.GetAtom(0), m.GetAtom(2));
   CheckNearAbs(angle * RAD2DEG, 90.0, 1e-2, "GetBondAngle should return 90 degrees for orthogonal bonds");

   m.AddBond(m.GetAtom(2), m.GetAtom(3), 1.0, 0.01, 0.02, 1.0, false);
   m.AddDihedralAngle(m.GetAtom(1), m.GetAtom(0), m.GetAtom(2), m.GetAtom(3), 0.0, 1.0, 1.0, false);
   Check(m.GetDihedralAngleList().size() == 1, "Molecule dihedral-angle insertion failed");
}

void TestMoleculeFormulaAndLogLikelihood()
{
   using namespace ObjCryst;
   Crystal c(10, 10, 10, "P1");
   Molecule m = MakeWaterLikeMolecule(c);
   m.AddBond(m.GetAtom(0), m.GetAtom(1), 1.0, 0.01, 0.02, 1.0, false);
   Check(!m.GetFormula().empty(), "Molecule formula should not be empty");
   // With no restraint violation, the log-likelihood (penalty) should be finite and near zero.
   const REAL llk = m.GetLogLikelihood();
   Check(std::isfinite(llk), "Molecule log-likelihood should be finite");
   Check(llk >= 0, "Molecule log-likelihood (restraint cost) should be non-negative");
}

} // namespace

int main(int argc, char* argv[])
{
   if(argc != 2)
   {
      std::cerr << "Usage: api_molecule <test-case>" << std::endl;
      return 2;
   }

   const std::string testName = argv[1];

   if(testName == "molecule-atoms-bonds") TestMoleculeAtomsAndBonds();
   else if(testName == "molecule-angles-dihedrals") TestMoleculeAnglesAndDihedrals();
   else if(testName == "molecule-formula-loglikelihood") TestMoleculeFormulaAndLogLikelihood();
   else
   {
      std::cerr << "Unknown test case: " << testName << std::endl;
      return 2;
   }

   return 0;
}
