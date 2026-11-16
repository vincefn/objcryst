// Extended unit tests for the Crystal, UnitCell, Scatterer and Atom classes.
#include <cmath>
#include <iostream>
#include <string>

#include "ObjCryst/ObjCryst/Atom.h"
#include "ObjCryst/ObjCryst/Crystal.h"
#include "ObjCryst/ObjCryst/Molecule.h"
#include "ObjCryst/ObjCryst/ScatteringPower.h"
#include "ObjCryst/ObjCryst/UnitCell.h"
#include "ObjCryst/RefinableObj/IO.h"

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

std::vector<ObjCryst::RefinableObj*> gReleasedObj;
bool gKeepReleasedObj = false;

bool RecordReleasedObj(ObjCryst::RefinableObj* obj)
{
   gReleasedObj.push_back(obj);
   return gKeepReleasedObj;
}

// Crystal::RemoveScatterer(.., true) must go through the release hook.
void TestCrystalReleaseHookRemoveScatterer()
{
   using namespace ObjCryst;
   Crystal c = MakePbso4Crystal();
   Scatterer* const o3 = &const_cast<Scatterer&>(c.GetScatt("O3"));
   gReleasedObj.clear();
   gKeepReleasedObj = true;
   const RefinableObjReleaseHook previous = SetReleaseObjHook(&RecordReleasedObj);
   Check(previous == 0, "No RefinableObj release hook should be installed by default");
   c.RemoveScatterer(o3, true);
   SetReleaseObjHook(RefinableObjReleaseHook(0));
   Check(gReleasedObj.size() == 1, "RemoveScatterer(del=true) should call the release hook once");
   // Scatterer uses virtual inheritance from RefinableObj: the most-derived
   // object is recovered with dynamic_cast<void*>.
   Check(dynamic_cast<void*>(gReleasedObj[0]) == dynamic_cast<void*>(o3),
         "Release hook should receive the removed scatterer");
   Check(dynamic_cast<Atom*>(gReleasedObj[0]) != 0, "Released scatterer should still be an Atom");
   Check(o3->GetName() == "O3", "Scatterer kept by the release hook should still be alive");
   delete o3;
}

// Crystal::ConnectAtoms() must respect SetDeleteSubObjInDestructor(false),
// and release the replaced atoms through the hook otherwise.
void TestCrystalConnectAtomsRelease(const bool deleteSubObj)
{
   using namespace ObjCryst;
   auto* spC = new ScatteringPowerAtom("C", "C");
   auto* c = new Crystal(10, 10, 10, "P1");
   c->SetDeleteSubObjInDestructor(deleteSubObj);
   c->AddScatteringPower(spC);
   Atom* atoms[3] = {new Atom(0.10, 0.10, 0.10, "C1", spC),
                     new Atom(0.25, 0.10, 0.10, "C2", spC),
                     new Atom(0.25, 0.25, 0.10, "C3", spC)};
   for(Atom* a : atoms) c->AddScatterer(a);

   gReleasedObj.clear();
   gKeepReleasedObj = false;
   SetReleaseObjHook(&RecordReleasedObj);
   c->ConnectAtoms();
   SetReleaseObjHook(RefinableObjReleaseHook(0));

   Check(c->GetNbScatterer() == 1, "ConnectAtoms should replace the 3 atoms by a single Molecule");
   Molecule* const mol = dynamic_cast<Molecule*>(&const_cast<Scatterer&>(c->GetScatt(0)));
   Check(mol != 0, "ConnectAtoms should create a Molecule");
   Check(mol->GetNbComponent() == 3, "ConnectAtoms Molecule should have 3 atoms");
   if(deleteSubObj)
   {
      Check(gReleasedObj.size() == 3, "ConnectAtoms should release the 3 atoms through the hook");
      delete c;
   }
   else
   {
      Check(gReleasedObj.empty(), "ConnectAtoms should not release atoms if SetDeleteSubObjInDestructor(false)");
      for(Atom* a : atoms) Check(a->GetName().size() == 2, "Atom should still be alive");
      delete c;
      delete mol;
      for(Atom* a : atoms) delete a;
      delete spC;
   }
}

// Crystal::XMLInput() on a non-empty Crystal must release the scatterers
// before the scattering powers they reference.
void TestCrystalXMLInputReleaseOrder()
{
   using namespace ObjCryst;
   std::stringstream xml;
   MakePbso4Crystal().XMLOutput(xml, 0);
   const std::string xmlText = xml.str();

   Crystal restored;
   for(int i = 0; i < 2; ++i)
   {
      std::stringstream is(xmlText);
      XMLCrystTag tag(is);
      restored.XMLInput(is, tag);
   }
   Check(restored.GetNbScatterer() == 5, "Crystal::XMLInput should restore the 5 PbSO4 atoms");

   gReleasedObj.clear();
   gKeepReleasedObj = true;
   SetReleaseObjHook(&RecordReleasedObj);
   {
      std::stringstream is(xmlText);
      XMLCrystTag tag(is);
      restored.XMLInput(is, tag);
   }
   SetReleaseObjHook(RefinableObjReleaseHook(0));
   Check(gReleasedObj.size() == 5 + 3, "Crystal::XMLInput should release 5 atoms and 3 scattering powers");
   bool seenPower = false;
   for(RefinableObj* obj : gReleasedObj)
   {
      if(dynamic_cast<ScatteringPower*>(obj) != 0) seenPower = true;
      else Check(!seenPower, "Crystal::XMLInput should release scatterers before scattering powers");
   }
   for(RefinableObj* obj : gReleasedObj) delete obj;
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
   else if(testName == "crystal-release-hook-remove-scatterer") TestCrystalReleaseHookRemoveScatterer();
   else if(testName == "crystal-connectatoms-release") TestCrystalConnectAtomsRelease(true);
   else if(testName == "crystal-connectatoms-no-delete") TestCrystalConnectAtomsRelease(false);
   else if(testName == "crystal-xmlinput-release-order") TestCrystalXMLInputReleaseOrder();
   else
   {
      std::cerr << "Unknown test case: " << testName << std::endl;
      return 2;
   }

   return 0;
}
