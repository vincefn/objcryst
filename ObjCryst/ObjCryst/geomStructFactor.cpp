//:NOTE: It may be a good idea to use a static array to compute the structure
//factor, so that less time be spent in the construction of the array.
//For example by using mutable member (assuming h,k and l's are always CrystVector_double),
//:NOTE: Or the result could be returned in an array given as a parameter

//NOTA BENE : normally, the formatting of the equations is the same as in
//the Int. Tables for X-Ray Crystallography (1969) : one line of equation
//should correspond to one line in the table, for easier check.

#include "CrystVector/CrystVector.h"

namespace ObjCryst
{


void RealGeomStructFactor      (const double x,
                                const double y,
                                const double z,
                                const CrystVector_double&h,
                                const CrystVector_double&k,
                                const CrystVector_double&l,
                                CrystVector_double& rsf)
{
};

void ImagGeomStructFactor      (const double x,
                                const double y,
                                const double z,
                                const CrystVector_double&h,
                                const CrystVector_double&k,
                                const CrystVector_double&l,
                                CrystVector_double& isf)
{
};

}//namespace
