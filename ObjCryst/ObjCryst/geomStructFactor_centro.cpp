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

void ImagGeomStructFactor_centro(const double x,
                                const double y,
                                const double z,
                                const CrystVector_double&h,
                                const CrystVector_double&k,
                                const CrystVector_double&l,
                                CrystVector_double& isf){return;}; //do nothing

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
