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

void RealGeomStructFactor_97   (const REAL x,
                                const REAL y,
                                const REAL z,
                                const CrystVector_REAL&hh,
                                const CrystVector_REAL&kk,
                                const CrystVector_REAL&ll,
                                CrystVector_REAL& sfsf)
{
__VFN_GEOM_STRUCT_FACTOR_POINTERS_INIT

      SF+= 8*pow(cos((H+K+L)/4),2)*cos(L*z)*(cos(H*x)*cos(K*y)+cos(K*x)*cos(H*y));

__VFN_GEOM_STRUCT_FACTOR_POINTERS_END

   return;
};

void ImagGeomStructFactor_97   (const REAL x,
                                const REAL y,
                                const REAL z,
                                const CrystVector_REAL&hh,
                                const CrystVector_REAL&kk,
                                const CrystVector_REAL&ll,
                                CrystVector_REAL& sfsf)
{
__VFN_GEOM_STRUCT_FACTOR_POINTERS_INIT

   SF += (-8)*pow(cos((H+K+L)/4),2)*sin(L*z)*(sin(H*x)*sin(K*y)-sin(K*x)*sin(H*y));

__VFN_GEOM_STRUCT_FACTOR_POINTERS_END

   return;
};

}
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
