#include "CrystVector/CrystVector.h"
#ifdef __VFN_GEOM_STRUCT_FACTOR_USE_POINTERS

#define H (*h)
#define K (*k)
#define L (*l)
#define SF (*sf)

#define __VFN_GEOM_STRUCT_FACTOR_POINTERS_INIT const double *h,*k,*l; double*sf;\
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

void RealGeomStructFactor_230  (const double x,
                                const double y,
                                const double z,
                                const CrystVector_double&hh,
                                const CrystVector_double&kk,
                                const CrystVector_double&ll,
                                CrystVector_double& sfsf)
{
__VFN_GEOM_STRUCT_FACTOR_POINTERS_INIT

   SF+=16*cos((H+K+L)/4)*(cos(H*x+L/4)*cos(K*y+H/4)*cos(L*z+K/4)
       +cos(K*x+H/4)*cos(L*y+K/4)*cos(H*z+L/4)+cos(L*x+K/4)*cos(H*y+L/4)*cos(K*z+H/4)
       +cos((H+K+L)/4)*(cos(K*x+L/4)*cos(H*y+K/4)*cos(L*z+H/4)
       +cos(L*x+H/4)*cos(K*y+L/4)*cos(H*z+K/4)+cos(H*x+K/4)*cos(L*y+H/4)*cos(K*z+L/4)));

__VFN_GEOM_STRUCT_FACTOR_POINTERS_END
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
