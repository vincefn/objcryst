/*  ObjCryst++ Object-Oriented Crystallographic Library
    (c) 2000-2002 Vincent Favre-Nicolin vincefn@users.sourceforge.net
        2000-2001 University of Geneva (Switzerland)

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
//:NOTE: It may be a good idea to use a static array to compute the structure
//factor, so that less time be spent in the construction of the array.
//For example by using mutable member (assuming h,k and l's are always CrystVector_REAL),
//:NOTE: Or the result could be returned in an array given as a parameter

//NOTA BENE : normally, the formatting of the equations is the same as in
//the Int. Tables for X-Ray Crystallography (1969) : one line of equation
//should correspond to one line in the table, for easier check.

#include "CrystVector/CrystVector.h"

namespace ObjCryst
{


void RealGeomStructFactor      (const REAL x,
                                const REAL y,
                                const REAL z,
                                const CrystVector_REAL&h,
                                const CrystVector_REAL&k,
                                const CrystVector_REAL&l,
                                CrystVector_REAL& rsf)
{
};

void ImagGeomStructFactor      (const REAL x,
                                const REAL y,
                                const REAL z,
                                const CrystVector_REAL&h,
                                const CrystVector_REAL&k,
                                const CrystVector_REAL&l,
                                CrystVector_REAL& isf)
{
};

}//namespace
