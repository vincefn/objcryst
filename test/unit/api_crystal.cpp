// Extended unit tests for the Crystal, UnitCell, Scatterer and Atom classes.
#include <cmath>
#include <iostream>
#include <string>

#include "ObjCryst/ObjCryst/Atom.h"
#include "ObjCryst/ObjCryst/Crystal.h"
#include "ObjCryst/ObjCryst/ScatteringPower.h"
#include "ObjCryst/ObjCryst/UnitCell.h"

#include "test_common.h"

using namespace objcryst_test;

namespace
{

void TestScatteringPowerAtom()
{
   using namespace ObjCryst;
   ScatteringPowerAtom o("O", "O", 1.2);
   Check(o.GetAtomicNumber() == 8, "Atomic number mismatch for oxygen");
   Check(o.GetAtomicWeight() > 10.0, "Atomic weight for oxygen invalid");
   Check(o.GetRadius() > 0.1, "Atomic radius should be positive");
   Check(o.GetCovalentRadius() > 0.1, "Covalent radius should be positive");
   Check(o.GetForwardScatteringFactor(RAD_XRAY) > 0, "Forward scattering factor should be positive");
   Check(o.GetForwardScatteringFactor(RAD_NEUTRON) != 0, "Neutron scattering length should be non-zero");
   Check(o.GetSymbol() == "O", "ScatteringPowerAtom symbol mismatch");
}

void TestCrystalAndAtom()
{
   using namespace ObjCryst;
   Crystal c = MakePbso4Crystal();
   Check(c.GetNbScatterer() == 5, "Crystal should contain 5 scatterers");
   Check(c.GetScatteringPowerRegistry().GetNb() == 3, "Crystal should contain 3 scattering powers");

   Atom& pb = dynamic_cast<Atom&>(c.GetScatt(0));
   pb.SetX(0.2);
   pb.SetY(0.3);
   pb.SetZ(0.4);
   Check(std::fabs(pb.GetX() - 0.2) < 1e-6, "Atom X coordinate update failed");
   Check(pb.GetMass() > 1., "Atom mass should be positive");
   Check(!pb.IsDummy(), "Pb atom should not be dummy");

   // Restore original position so the crystal remains a valid PbSO4 structure.
   pb.SetX(.1882);
   pb.SetY(.250);
   pb.SetZ(.167);

   Check(&c.GetScatt("Pb1") == &c.GetScatt(0), "Crystal::GetScatt(name) should match GetScatt(index)");
}

void TestCrystalScattererManagement()
{
   using namespace ObjCryst;
   Crystal c = MakePbso4Crystal();
   const long nbBefore = c.GetNbScatterer();
   Scatterer& toRemove = const_cast<Scatterer&>(c.GetScatt("O3"));
   c.RemoveScatterer(&toRemove, true);
   Check(c.GetNbScatterer() == nbBefore - 1, "Crystal::RemoveScatterer should decrement scatterer count");

   const CrystMatrix_REAL distTable = c.GetMinDistanceTable(0.1);
   Check(distTable.rows() == c.GetNbScatterer(), "Min-distance table rows should match scatterer count");
   Check(distTable.cols() == c.GetNbScatterer(), "Min-distance table cols should match scatterer count");
}

void TestUnitCellGeometry()
{
   using namespace ObjCryst;
   UnitCell cell(8.516, 5.399, 6.989, "Pnma");
   Check(cell.GetSpaceGroup().GetSpaceGroupNumber() == 62, "UnitCell spacegroup mismatch for Pnma");

   CheckNearAbs(cell.GetLatticePar(0), 8.516, 1e-4, "UnitCell lattice parameter a mismatch");
   CheckNearAbs(cell.GetLatticePar(1), 5.399, 1e-4, "UnitCell lattice parameter b mismatch");
   CheckNearAbs(cell.GetLatticePar(2), 6.989, 1e-4, "UnitCell lattice parameter c mismatch");
   CheckNearAbs(cell.GetLatticePar(3), M_PI / 2, 1e-4, "UnitCell alpha angle should be 90 degrees for Pnma");

   const REAL expectedVolume = 8.516 * 5.399 * 6.989;
   CheckNearAbsRel(cell.GetVolume(), expectedVolume, 1e-2, 1e-4, "UnitCell volume mismatch for orthorhombic cell");

   const CrystMatrix_REAL& bMatrix = cell.GetBMatrix();
   Check(bMatrix.rows() == 3 && bMatrix.cols() == 3, "B matrix should be 3x3");
   const CrystMatrix_REAL& orthMatrix = cell.GetOrthMatrix();
   Check(orthMatrix.rows() == 3 && orthMatrix.cols() == 3, "Orthogonalization matrix should be 3x3");

   REAL x = 0.123, y = 0.234, z = 0.345;
   const REAL x0 = x, y0 = y, z0 = z;
   cell.FractionalToOrthonormalCoords(x, y, z);
   cell.OrthonormalToFractionalCoords(x, y, z);
   CheckNearAbs(x, x0, 1e-5, "Fractional/orthonormal round-trip failed for x");
   CheckNearAbs(y, y0, 1e-5, "Fractional/orthonormal round-trip failed for y");
   CheckNearAbs(z, z0, 1e-5, "Fractional/orthonormal round-trip failed for z");

   REAL mx = 1, my = 0, mz = 0;
   cell.MillerToOrthonormalCoords(mx, my, mz);
   cell.OrthonormalToMillerCoords(mx, my, mz);
   CheckNearAbs(mx, 1, 1e-4, "Miller/orthonormal round-trip failed for h");
}

} // namespace

int main(int argc, char* argv[])
{
   if(argc != 2)
   {
      std::cerr << "Usage: api_crystal <test-case>" << std::endl;
      return 2;
   }

   const std::string testName = argv[1];

   if(testName == "scattering-power-atom") TestScatteringPowerAtom();
   else if(testName == "crystal-atom") TestCrystalAndAtom();
   else if(testName == "crystal-scatterer-management") TestCrystalScattererManagement();
   else if(testName == "unitcell-geometry") TestUnitCellGeometry();
   else
   {
      std::cerr << "Unknown test case: " << testName << std::endl;
      return 2;
   }

   return 0;
}
