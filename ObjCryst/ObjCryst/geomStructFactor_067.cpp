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

void RealGeomStructFactor_67   (const double x,
                                const double y,
                                const double z,
                                const CrystVector_double&hh,
                                const CrystVector_double&kk,
                                const CrystVector_double&ll,
                                CrystVector_double& sfsf)
{
__VFN_GEOM_STRUCT_FACTOR_POINTERS_INIT

   SF+=16*pow(cos((H+K)/4),2)*cos(H*x)*cos(K*y+H/4)*cos(L*z-H/4);

__VFN_GEOM_STRUCT_FACTOR_POINTERS_END
};

void RealGeomStructFactor_67ba_c(const double x,
                                const double y,
                                const double z,
                                const CrystVector_double&hh,
                                const CrystVector_double&kk,
                                const CrystVector_double&ll,
                                CrystVector_double& sfsf)
{
__VFN_GEOM_STRUCT_FACTOR_POINTERS_INIT

   SF+=16*pow(cos((K+H)/4),2)*cos(K*y)*cos(H*x+K/4)*cos(L*z-K/4);

__VFN_GEOM_STRUCT_FACTOR_POINTERS_END
};

void RealGeomStructFactor_67cab(const double x,
                                const double y,
                                const double z,
                                const CrystVector_double&hh,
                                const CrystVector_double&kk,
                                const CrystVector_double&ll,
                                CrystVector_double& sfsf)
{
__VFN_GEOM_STRUCT_FACTOR_POINTERS_INIT

   SF+=16*pow(cos((L+K)/4),2)*cos(K*y)*cos(L*z+K/4)*cos(H*x-K/4);

__VFN_GEOM_STRUCT_FACTOR_POINTERS_END
};

void RealGeomStructFactor_67_cba(const double x,
                                const double y,
                                const double z,
                                const CrystVector_double&hh,
                                const CrystVector_double&kk,
                                const CrystVector_double&ll,
                                CrystVector_double& sfsf)
{
__VFN_GEOM_STRUCT_FACTOR_POINTERS_INIT

   SF+=16*pow(cos((-L+K)/4),2)*cos(L*z)*cos(K*y-L/4)*cos(H*x+L/4);

__VFN_GEOM_STRUCT_FACTOR_POINTERS_END
};

void RealGeomStructFactor_67bca(const double x,
                                const double y,
                                const double z,
                                const CrystVector_double&hh,
                                const CrystVector_double&kk,
                                const CrystVector_double&ll,
                                CrystVector_double& sfsf)
{
__VFN_GEOM_STRUCT_FACTOR_POINTERS_INIT

   SF+=16*pow(cos((H+L)/4),2)*cos(L*z)*cos(H*x+L/4)*cos(K*y-L/4);

__VFN_GEOM_STRUCT_FACTOR_POINTERS_END
};

void RealGeomStructFactor_67a_cb(const double x,
                                const double y,
                                const double z,
                                const CrystVector_double&hh,
                                const CrystVector_double&kk,
                                const CrystVector_double&ll,
                                CrystVector_double& sfsf)
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
