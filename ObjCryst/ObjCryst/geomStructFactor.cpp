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
