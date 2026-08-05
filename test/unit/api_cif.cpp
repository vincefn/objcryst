// Extended unit tests for the CIF import machinery (CIF, CIFData classes).
#include <fstream>
#include <iostream>
#include <string>

#include "ObjCryst/ObjCryst/CIF.h"
#include "ObjCryst/ObjCryst/Crystal.h"

#include "test_common.h"

using namespace objcryst_test;

namespace
{

const char* const kCifPath = "../../test/data/cif/PbSO4-COD-1528837.cif";

void TestCifImport()
{
   using namespace ObjCryst;
   std::ifstream cifIn(kCifPath);
   Check(cifIn.good(), "Unable to open COD CIF test file");
   CIF cif(cifIn, true, false);
   Crystal c;
   Crystal* result = CreateCrystalFromCIF(cif, false, true, false, false, &c);
   Check(result == &c, "CreateCrystalFromCIF returned wrong pointer");
   Check(c.GetSpaceGroup().GetSpaceGroupNumber() == 62, "CIF spacegroup mismatch");
   Check(c.GetNbScatterer() >= 5, "CIF import created too few scatterers");
}

void TestCifDataFields()
{
   using namespace ObjCryst;
   std::ifstream cifIn(kCifPath);
   Check(cifIn.good(), "Unable to open COD CIF test file");
   // Parse without auto-interpretation so we can inspect the raw CIFData first.
   CIF cif(cifIn, false, false);
   Check(!cif.mvData.empty(), "CIF should contain at least one data block");

   CIFData& data = cif.mvData.begin()->second;
   Check(data.mvLatticePar.empty(), "Lattice parameters should not be extracted before ExtractUnitCell()");
   data.ExtractUnitCell(false);
   Check(data.mvLatticePar.size() == 6, "ExtractUnitCell should populate 6 lattice parameters");
   CheckNearAbs(data.mvLatticePar[0], 8.482, 0.05, "CIF lattice parameter a mismatch");
   CheckNearAbs(data.mvLatticePar[1], 5.398, 0.05, "CIF lattice parameter b mismatch");
   CheckNearAbs(data.mvLatticePar[2], 6.959, 0.05, "CIF lattice parameter c mismatch");

   data.ExtractSpacegroup(false);
   Check(data.mSpacegroupNumberIT == "62", "CIF spacegroup IT number mismatch");

   data.ExtractAtomicPositions(false);
   Check(!data.mvAtom.empty(), "CIF should expose parsed atom records");
   Check(!data.mvAtom[0].mLabel.empty(), "CIF atom record should have a non-empty label");
}

void TestCifCoordinateConversion()
{
   using namespace ObjCryst;
   std::ifstream cifIn(kCifPath);
   Check(cifIn.good(), "Unable to open COD CIF test file");
   CIF cif(cifIn, false, false);
   CIFData& data = cif.mvData.begin()->second;
   data.ExtractUnitCell(false);
   data.CalcMatrices(false);

   REAL x = 0.1, y = 0.2, z = 0.3;
   const REAL x0 = x, y0 = y, z0 = z;
   data.f2c(x, y, z);
   data.c2f(x, y, z);
   CheckNearAbs(x, x0, 1e-4, "CIF fractional/cartesian round-trip failed for x");
   CheckNearAbs(y, y0, 1e-4, "CIF fractional/cartesian round-trip failed for y");
   CheckNearAbs(z, z0, 1e-4, "CIF fractional/cartesian round-trip failed for z");
}

void CheckCifParseFailure(const std::string& contents, const char* const message)
{
   using namespace ObjCryst;
   std::istringstream cifIn(contents);
   try
   {
      CIF cif(cifIn, false, false);
   }
   catch(const ObjCrystException&)
   {
      return;
   }
   Check(false, message);
}

void TestCifTruncatedValues()
{
   CheckCifParseFailure("_a\n", "CIF parser accepted a tag with no value");
   CheckCifParseFailure("_a\n# no value", "CIF parser accepted a comment instead of a tag value");
   CheckCifParseFailure("data_x\n_publ_section_title\n;\nunterminated\n",
                        "CIF parser accepted an unterminated semicolon text field");
   CheckCifParseFailure("data_x\n_chemical_name_common 'unterminated",
                        "CIF parser accepted an unterminated quoted string");
}

} // namespace

int main(int argc, char* argv[])
{
   if(argc != 2)
   {
      std::cerr << "Usage: api_cif <test-case>" << std::endl;
      return 2;
   }

   const std::string testName = argv[1];

   if(testName == "cif-import") TestCifImport();
   else if(testName == "cif-data-fields") TestCifDataFields();
   else if(testName == "cif-coordinate-conversion") TestCifCoordinateConversion();
   else if(testName == "cif-truncated-values") TestCifTruncatedValues();
   else
   {
      std::cerr << "Unknown test case: " << testName << std::endl;
      return 2;
   }

   return 0;
}
