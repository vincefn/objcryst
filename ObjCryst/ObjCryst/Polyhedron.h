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
/*   Polyhedron.h
*  header file for the Polyhdron creation
*
*/
#include "ObjCryst/ObjCryst/ScatteringPower.h"
#include "ObjCryst/ObjCryst/Crystal.h"
#include "ObjCryst/ObjCryst/Molecule.h"
namespace ObjCryst
{
Molecule* MakeTetrahedron(Crystal &cryst,const string &name,
                      const ScatteringPower *centralAtom,
                      const ScatteringPower *peripheralAtom,
                      const REAL dist);

Molecule* MakeOctahedron(Crystal &cryst,const string &name,
                     const ScatteringPower *centralAtom,
                     const ScatteringPower *peripheralAtom,
                     const REAL dist);

Molecule* MakeSquarePlane(Crystal &cryst,const string &name,
                     const ScatteringPower *centralAtom,
                     const ScatteringPower *peripheralAtom,
                     const REAL dist);

Molecule* MakeCube(Crystal &cryst,const string &name,
                     const ScatteringPower *centralAtom,
                     const ScatteringPower *peripheralAtom,
                     const REAL dist);

Molecule* MakeAntiPrismTetragonal(Crystal &cryst,const string &name,
                     const ScatteringPower *centralAtom,
                     const ScatteringPower *peripheralAtom,
                     const REAL dist);

Molecule* MakePrismTrigonal(Crystal &cryst,const string &name,
                     const ScatteringPower *centralAtom,
                     const ScatteringPower *peripheralAtom,
                     const REAL dist);

Molecule* MakeIcosahedron(Crystal &cryst,const string &name,
                     const ScatteringPower *centralAtom,
                     const ScatteringPower *peripheralAtom,
                     const REAL dist);

Molecule* MakeTriangle(Crystal &cryst,const string &name,
                     const ScatteringPower *centralAtom,
                     const ScatteringPower *peripheralAtom,
                     const REAL dist);

}//namespace
