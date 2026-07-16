#include <cmath>
#include <iostream>

#include "ObjCryst/ObjCryst/SpaceGroup.h"
#include "ObjCryst/ObjCryst/UnitCell.h"

int main()
{
   using namespace ObjCryst;

   UnitCell cell(8.482, 5.398, 6.959, "Pnma");
   if(cell.GetSpaceGroup().GetSpaceGroupNumber() != 62)
   {
      std::cerr << "Unexpected spacegroup number for Pnma" << std::endl;
      return 1;
   }

   const REAL volume = cell.GetVolume();
   if(volume < 50. || volume > 1000.)
   {
      std::cerr << "Unexpected unit cell volume: " << volume << std::endl;
      return 1;
   }

   REAL x = 0.123, y = 0.234, z = 0.345;
   const REAL x0 = x, y0 = y, z0 = z;
   cell.FractionalToOrthonormalCoords(x, y, z);
   cell.OrthonormalToFractionalCoords(x, y, z);
   if(std::fabs(x - x0) > 1e-5 || std::fabs(y - y0) > 1e-5 || std::fabs(z - z0) > 1e-5)
   {
      std::cerr << "Fractional/orthonormal coordinate round-trip failed" << std::endl;
      return 1;
   }

   SpaceGroup p1("P1");
   if(p1.GetSpaceGroupNumber() != 1 || p1.IsCentrosymmetric())
   {
      std::cerr << "Unexpected properties for P1" << std::endl;
      return 1;
   }

   SpaceGroup pminus1("P-1");
   if(pminus1.GetSpaceGroupNumber() != 2 || !pminus1.IsCentrosymmetric())
   {
      std::cerr << "Unexpected properties for P-1" << std::endl;
      return 1;
   }

   return 0;
}
