/*  ObjCryst++ Object-Oriented Crystallographic Library
    (c) 2000-2002 Vincent Favre-Nicolin vincefn@users.sourceforge.net
        2000-2001 University of Geneva (Switzerland)

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation;  version 2 of the License.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/*   test.h
*  header file for test functions (speed, etc...)
*
*/
#ifndef _OBJCRYST_TEST_H_
#define _OBJCRYST_TEST_H_

#include "ObjCryst/ObjCryst/General.h"

namespace ObjCryst
{
/** Structure to hold the results of a speedtest (see ObjCryst::SpeedTest())
*
*
*/
struct SpeedTestReport
{
   /// Total number of unique atoms in the test structure.
   unsigned int mNbAtom;
   /// Total number of atom types in the test structure.
   int mNbAtomType;
   /// The symbol for the spacegroup
   string mSpacegroup;
   /// The type of radiation used
   RadiationType mRadiation;
   /// The total number of reflections used for the tests
   unsigned long mNbReflections;
   /// dataType: 0= single crystal, 1= powder pattern (1 background + 1 crystal phase)
   unsigned int mDataType;
   /// Million of Reflections-Atoms computed Per Second (considering all atoms in the unit cell)
   REAL mBogoMRAPS;
   /// Million of Reflections-Atoms computed Per Second (considering all atoms in the unit cell,
   /// except those deduced by a center of symmetry or a lattice translation)
   REAL mBogoMRAPS_reduced;
   /// Number of Structures evaluated Per Second
   REAL mBogoSPS;
};

/**
*
* \param nbAtom: total number of unique atoms.
* \param nbAtomType: total number of atom types. Must be smaller than nbAtom.
* \param spacegroup: the symbol or spacegroup number for the spacegroup to be tested.
* \param radiation: the type of radiation to do the test with (neutron, x-ray).
* \param nbReflections: total number of reflections to use.
* \param time: duration of the test (10s is usually enough).
* \param dataType: 0= single crystal, 1= powder pattern (1 background + 1 crystal phase)
*/
SpeedTestReport SpeedTest(const unsigned int nbAtom, const int nbAtomType,const string spacegroup,
                          const RadiationType radiation, const unsigned long nbReflections,
                          const unsigned int dataType,const REAL time);

}
#endif
