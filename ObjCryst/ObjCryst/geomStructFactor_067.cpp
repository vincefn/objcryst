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
#include "CrystVector/CrystVector.h"
#ifdef __VFN_GEOM_STRUCT_FACTOR_USE_POINTERS

#define H (*h)
#define K (*k)
#define L (*l)
#define SF (*sf)

#define __VFN_GEOM_STRUCT_FACTOR_POINTERS_INIT const REAL *h,*k,*l; REAL*sf;\
   h=hh.data();k=kk.data();l=ll.data();sf=sfsf.data(); for(long i=0;i<hh.numElements();i++){
#define __VFN_GEOM_STRUCT_FACTOR_POINTERS_END h++ ; k++ ; l++ ; sf++; };

#else

#define H hh
#define K kk
#define L ll
#define sfsf SF

#define __VFN_GEOM_STRUCT_FACTOR_POINTERS_INIT
#define __VFN_GEOM_STRUCT_FACTOR_POINTERS_END

#endif

namespace ObjCryst
{

void RealGeomStructFactor_67   (const REAL x,
                                const REAL y,
                                const REAL z,
                                const CrystVector_REAL&hh,
                                const CrystVector_REAL&kk,
                                const CrystVector_REAL&ll,
                                CrystVector_REAL& sfsf)
{
__VFN_GEOM_STRUCT_FACTOR_POINTERS_INIT

   SF+=16*pow(cos((H+K)/4),2)*cos(H*x)*cos(K*y+H/4)*cos(L*z-H/4);

__VFN_GEOM_STRUCT_FACTOR_POINTERS_END
};

void RealGeomStructFactor_67ba_c(const REAL x,
                                const REAL y,
                                const REAL z,
                                const CrystVector_REAL&hh,
                                const CrystVector_REAL&kk,
                                const CrystVector_REAL&ll,
                                CrystVector_REAL& sfsf)
{
__VFN_GEOM_STRUCT_FACTOR_POINTERS_INIT

   SF+=16*pow(cos((K+H)/4),2)*cos(K*y)*cos(H*x+K/4)*cos(L*z-K/4);

__VFN_GEOM_STRUCT_FACTOR_POINTERS_END
};

void RealGeomStructFactor_67cab(const REAL x,
                                const REAL y,
                                const REAL z,
                                const CrystVector_REAL&hh,
                                const CrystVector_REAL&kk,
                                const CrystVector_REAL&ll,
                                CrystVector_REAL& sfsf)
{
__VFN_GEOM_STRUCT_FACTOR_POINTERS_INIT

   SF+=16*pow(cos((L+K)/4),2)*cos(K*y)*cos(L*z+K/4)*cos(H*x-K/4);

__VFN_GEOM_STRUCT_FACTOR_POINTERS_END
};

void RealGeomStructFactor_67_cba(const REAL x,
                                const REAL y,
                                const REAL z,
                                const CrystVector_REAL&hh,
                                const CrystVector_REAL&kk,
                                const CrystVector_REAL&ll,
                                CrystVector_REAL& sfsf)
{
__VFN_GEOM_STRUCT_FACTOR_POINTERS_INIT

   SF+=16*pow(cos((-L+K)/4),2)*cos(L*z)*cos(K*y-L/4)*cos(H*x+L/4);

__VFN_GEOM_STRUCT_FACTOR_POINTERS_END
};

void RealGeomStructFactor_67bca(const REAL x,
                                const REAL y,
                                const REAL z,
                                const CrystVector_REAL&hh,
                                const CrystVector_REAL&kk,
                                const CrystVector_REAL&ll,
                                CrystVector_REAL& sfsf)
{
__VFN_GEOM_STRUCT_FACTOR_POINTERS_INIT

   SF+=16*pow(cos((H+L)/4),2)*cos(L*z)*cos(H*x+L/4)*cos(K*y-L/4);

__VFN_GEOM_STRUCT_FACTOR_POINTERS_END
};

void RealGeomStructFactor_67a_cb(const REAL x,
                                const REAL y,
                                const REAL z,
                                const CrystVector_REAL&hh,
                                const CrystVector_REAL&kk,
                                const CrystVector_REAL&ll,
                                CrystVector_REAL& sfsf)
{
__VFN_GEOM_STRUCT_FACTOR_POINTERS_INIT

   SF+=16*pow(cos((H-L)/4),2)*cos(H*x)*cos(L*z+H/4)*cos(K*y-H/4);

__VFN_GEOM_STRUCT_FACTOR_POINTERS_END
};

}//namespace
#undef H
#undef K
#undef L
#undef SF
#undef hh
#undef kk
#undef ll
#undef sfsf
#undef __VFN_GEOM_STRUCT_FACTOR_POINTERS_INIT
#undef __VFN_GEOM_STRUCT_FACTOR_POINTERS_END

#undef __VFN_GEOM_STRUCT_FACTOR_USE_POINTERS
