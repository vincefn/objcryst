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
// This file declares all the functions used to compute the geometrical
//structure factors

//Use pointers for the calculation of geometrical structure factors ?
//Else use the blitz library...
//This is *not* useful for performance reasons, but rather
//for memory /disk space reasons when using the GNU gcc compiler.
//(complex array expressions using blitz take a *huge* space
//on disk and in memory when compiling). Performance is still the best...

#ifndef _VFN_RC_CRISTALLO_GEOM_STRUCT_FACTOR_H_
#define _VFN_RC_CRISTALLO_GEOM_STRUCT_FACTOR_H_

#include "ObjCryst/CrystVector/CrystVector.h"

#include <cmath>

namespace ObjCryst
{

// This is the default fonction, which does *not* use
//the geometrical structure factor, but calculates all the symmetric
//positions of the atom. This is done for SpaceGroups which
//do not yet have a coded Geom Structure factor, or for those
//where it does not make much difference
/// \deprecated
void RealGeomStructFactor      (const REAL x,
                                const REAL y,
                                const REAL z,
                                const CrystVector_REAL&h,
                                const CrystVector_REAL&k,
                                const CrystVector_REAL&l,
                                CrystVector_REAL& rsf);
// Same for the imaginary part
/// \deprecated
void ImagGeomStructFactor      (const REAL x,
                                const REAL y,
                                const REAL z,
                                const CrystVector_REAL&h,
                                const CrystVector_REAL&k,
                                const CrystVector_REAL&l,
                                CrystVector_REAL& rsf);


/// \deprecated
void RealGeomStructFactor_1    (const REAL x,
                                const REAL y,
                                const REAL z,
                                const CrystVector_REAL&h,
                                const CrystVector_REAL&k,
                                const CrystVector_REAL&l,
                                CrystVector_REAL& rsf);

/// \deprecated
void RealGeomStructFactor_2    (const REAL x,
                                const REAL y,
                                const REAL z,
                                const CrystVector_REAL&h,
                                const CrystVector_REAL&k,
                                const CrystVector_REAL&l,
                                CrystVector_REAL& rsf);

/// \deprecated
void RealGeomStructFactor_67   (const REAL x,
                                const REAL y,
                                const REAL z,
                                const CrystVector_REAL&h,
                                const CrystVector_REAL&k,
                                const CrystVector_REAL&l,
                                CrystVector_REAL& rsf);

/// \deprecated
void RealGeomStructFactor_67ba_c(const REAL x,
                                const REAL y,
                                const REAL z,
                                const CrystVector_REAL&h,
                                const CrystVector_REAL&k,
                                const CrystVector_REAL&l,
                                CrystVector_REAL& rsf);

/// \deprecated
void RealGeomStructFactor_67cab(const REAL x,
                                const REAL y,
                                const REAL z,
                                const CrystVector_REAL&h,
                                const CrystVector_REAL&k,
                                const CrystVector_REAL&l,
                                CrystVector_REAL& rsf);

/// \deprecated
void RealGeomStructFactor_67_cba(const REAL x,
                                const REAL y,
                                const REAL z,
                                const CrystVector_REAL&h,
                                const CrystVector_REAL&k,
                                const CrystVector_REAL&l,
                                CrystVector_REAL& rsf);

/// \deprecated
void RealGeomStructFactor_67bca(const REAL x,
                                const REAL y,
                                const REAL z,
                                const CrystVector_REAL&h,
                                const CrystVector_REAL&k,
                                const CrystVector_REAL&l,
                                CrystVector_REAL& rsf);

/// \deprecated
void RealGeomStructFactor_67a_cb(const REAL x,
                                const REAL y,
                                const REAL z,
                                const CrystVector_REAL&h,
                                const CrystVector_REAL&k,
                                const CrystVector_REAL&l,
                                CrystVector_REAL& rsf);

/// \deprecated
void RealGeomStructFactor_97   (const REAL x,
                                const REAL y,
                                const REAL z,
                                const CrystVector_REAL&h,
                                const CrystVector_REAL&k,
                                const CrystVector_REAL&l,
                                CrystVector_REAL& rsf);

/// \deprecated
void RealGeomStructFactor_230  (const REAL x,
                                const REAL y,
                                const REAL z,
                                const CrystVector_REAL&h,
                                const CrystVector_REAL&k,
                                const CrystVector_REAL&l,
                                CrystVector_REAL& rsf);

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

/// \deprecated
void ImagGeomStructFactor_centro(const REAL x,
                                const REAL y,
                                const REAL z,
                                const CrystVector_REAL&h,
                                const CrystVector_REAL&k,
                                const CrystVector_REAL&l,
                                CrystVector_REAL& isf); //do nothing

/// \deprecated
void ImagGeomStructFactor_1    (const REAL x,
                                const REAL y,
                                const REAL z,
                                const CrystVector_REAL&h,
                                const CrystVector_REAL&k,
                                const CrystVector_REAL&l,
                                CrystVector_REAL& isf);

/// \deprecated
void ImagGeomStructFactor_2    (const REAL x,
                                const REAL y,
                                const REAL z,
                                const CrystVector_REAL&h,
                                const CrystVector_REAL&k,
                                const CrystVector_REAL&l,
                                CrystVector_REAL& isf);

/// \deprecated
void ImagGeomStructFactor_67   (const REAL x,
                                const REAL y,
                                const REAL z,
                                const CrystVector_REAL&h,
                                const CrystVector_REAL&k,
                                const CrystVector_REAL&l,
                                CrystVector_REAL& isf);

/// \deprecated
void ImagGeomStructFactor_67ba_c(const REAL x,
                                const REAL y,
                                const REAL z,
                                const CrystVector_REAL&h,
                                const CrystVector_REAL&k,
                                const CrystVector_REAL&l,
                                CrystVector_REAL& rsf);

/// \deprecated
void ImagGeomStructFactor_67cab(const REAL x,
                                const REAL y,
                                const REAL z,
                                const CrystVector_REAL&h,
                                const CrystVector_REAL&k,
                                const CrystVector_REAL&l,
                                CrystVector_REAL& rsf);

/// \deprecated
void ImagGeomStructFactor_67_cba(const REAL x,
                                const REAL y,
                                const REAL z,
                                const CrystVector_REAL&h,
                                const CrystVector_REAL&k,
                                const CrystVector_REAL&l,
                                CrystVector_REAL& rsf);

/// \deprecated
void ImagGeomStructFactor_67bca(const REAL x,
                                const REAL y,
                                const REAL z,
                                const CrystVector_REAL&h,
                                const CrystVector_REAL&k,
                                const CrystVector_REAL&l,
                                CrystVector_REAL& rsf);

/// \deprecated
void ImagGeomStructFactor_67a_cb(const REAL x,
                                const REAL y,
                                const REAL z,
                                const CrystVector_REAL&h,
                                const CrystVector_REAL&k,
                                const CrystVector_REAL&l,
                                CrystVector_REAL& rsf);

/// \deprecated
void ImagGeomStructFactor_97   (const REAL x,
                                const REAL y,
                                const REAL z,
                                const CrystVector_REAL&h,
                                const CrystVector_REAL&k,
                                const CrystVector_REAL&l,
                                CrystVector_REAL& isf);

/// \deprecated
void ImagGeomStructFactor_230  (const REAL x,
                                const REAL y,
                                const REAL z,
                                const CrystVector_REAL&h,
                                const CrystVector_REAL&k,
                                const CrystVector_REAL&l,
                                CrystVector_REAL& isf);

}//namespace

#endif
