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

#include "CrystVector/CrystVector.h"

#include <cmath>

namespace ObjCryst
{

// This is the default fonction, which does *not* use
//the geometrical structure factor, but calculates all the symetric
//positions of the atom. This is done for SpaceGroups which
//do not yet have a coded Geom Structure factor, or for those
//where it does not make much difference
void RealGeomStructFactor      (const double x,
                                const double y,
                                const double z,
                                const CrystVector_double&h,
                                const CrystVector_double&k,
                                const CrystVector_double&l,
                                CrystVector_double& rsf);
// Same for the imaginary part
void ImagGeomStructFactor      (const double x,
                                const double y,
                                const double z,
                                const CrystVector_double&h,
                                const CrystVector_double&k,
                                const CrystVector_double&l,
                                CrystVector_double& rsf);


void RealGeomStructFactor_1    (const double x,
                                const double y,
                                const double z,
                                const CrystVector_double&h,
                                const CrystVector_double&k,
                                const CrystVector_double&l,
                                CrystVector_double& rsf);

void RealGeomStructFactor_2    (const double x,
                                const double y,
                                const double z,
                                const CrystVector_double&h,
                                const CrystVector_double&k,
                                const CrystVector_double&l,
                                CrystVector_double& rsf);

void RealGeomStructFactor_67   (const double x,
                                const double y,
                                const double z,
                                const CrystVector_double&h,
                                const CrystVector_double&k,
                                const CrystVector_double&l,
                                CrystVector_double& rsf);

void RealGeomStructFactor_67ba_c(const double x,
                                const double y,
                                const double z,
                                const CrystVector_double&h,
                                const CrystVector_double&k,
                                const CrystVector_double&l,
                                CrystVector_double& rsf);

void RealGeomStructFactor_67cab(const double x,
                                const double y,
                                const double z,
                                const CrystVector_double&h,
                                const CrystVector_double&k,
                                const CrystVector_double&l,
                                CrystVector_double& rsf);

void RealGeomStructFactor_67_cba(const double x,
                                const double y,
                                const double z,
                                const CrystVector_double&h,
                                const CrystVector_double&k,
                                const CrystVector_double&l,
                                CrystVector_double& rsf);

void RealGeomStructFactor_67bca(const double x,
                                const double y,
                                const double z,
                                const CrystVector_double&h,
                                const CrystVector_double&k,
                                const CrystVector_double&l,
                                CrystVector_double& rsf);
                                
void RealGeomStructFactor_67a_cb(const double x,
                                const double y,
                                const double z,
                                const CrystVector_double&h,
                                const CrystVector_double&k,
                                const CrystVector_double&l,
                                CrystVector_double& rsf);
                                
void RealGeomStructFactor_97   (const double x,
                                const double y,
                                const double z,
                                const CrystVector_double&h,
                                const CrystVector_double&k,
                                const CrystVector_double&l,
                                CrystVector_double& rsf);

void RealGeomStructFactor_230  (const double x,
                                const double y,
                                const double z,
                                const CrystVector_double&h,
                                const CrystVector_double&k,
                                const CrystVector_double&l,
                                CrystVector_double& rsf);

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

void ImagGeomStructFactor_centro(const double x,
                                const double y,
                                const double z,
                                const CrystVector_double&h,
                                const CrystVector_double&k,
                                const CrystVector_double&l,
                                CrystVector_double& isf); //do nothing

void ImagGeomStructFactor_1    (const double x,
                                const double y,
                                const double z,
                                const CrystVector_double&h,
                                const CrystVector_double&k,
                                const CrystVector_double&l,
                                CrystVector_double& isf);

void ImagGeomStructFactor_2    (const double x,
                                const double y,
                                const double z,
                                const CrystVector_double&h,
                                const CrystVector_double&k,
                                const CrystVector_double&l,
                                CrystVector_double& isf);

void ImagGeomStructFactor_67   (const double x,
                                const double y,
                                const double z,
                                const CrystVector_double&h,
                                const CrystVector_double&k,
                                const CrystVector_double&l,
                                CrystVector_double& isf);

void ImagGeomStructFactor_67ba_c(const double x,
                                const double y,
                                const double z,
                                const CrystVector_double&h,
                                const CrystVector_double&k,
                                const CrystVector_double&l,
                                CrystVector_double& rsf);

void ImagGeomStructFactor_67cab(const double x,
                                const double y,
                                const double z,
                                const CrystVector_double&h,
                                const CrystVector_double&k,
                                const CrystVector_double&l,
                                CrystVector_double& rsf);

void ImagGeomStructFactor_67_cba(const double x,
                                const double y,
                                const double z,
                                const CrystVector_double&h,
                                const CrystVector_double&k,
                                const CrystVector_double&l,
                                CrystVector_double& rsf);

void ImagGeomStructFactor_67bca(const double x,
                                const double y,
                                const double z,
                                const CrystVector_double&h,
                                const CrystVector_double&k,
                                const CrystVector_double&l,
                                CrystVector_double& rsf);
                                
void ImagGeomStructFactor_67a_cb(const double x,
                                const double y,
                                const double z,
                                const CrystVector_double&h,
                                const CrystVector_double&k,
                                const CrystVector_double&l,
                                CrystVector_double& rsf);
                                
void ImagGeomStructFactor_97   (const double x,
                                const double y,
                                const double z,
                                const CrystVector_double&h,
                                const CrystVector_double&k,
                                const CrystVector_double&l,
                                CrystVector_double& isf);

void ImagGeomStructFactor_230  (const double x,
                                const double y,
                                const double z,
                                const CrystVector_double&h,
                                const CrystVector_double&k,
                                const CrystVector_double&l,
                                CrystVector_double& isf);

}//namespace

#endif
