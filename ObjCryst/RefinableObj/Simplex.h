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
/*   Simplex.h
*  header file for Simplex Algorithm object
*
*/

#ifndef _SIMPLEX_H
#define _SIMPLEX_H

#include "RefinableObj/GlobalOptimObj.h"

namespace ObjCryst
{
/** Conjugate Gradient Algorithm object
*
* currently does not handle parameters hitting limits, and is not very efficient
* (uses numerical derivatives)
*/
class SimplexObj:public OptimizationObj
{
   public:
      /// Constructor
      SimplexObj(const string name="Unnamed Simplex Object");
      virtual void Optimize(long &nbSteps,const bool silent=false,const REAL finalcost=0,
                            const REAL maxTime=-1);
      virtual void MultiRunOptimize(long &nbCycle,long &nbSteps,const bool silent=false,
                                    const REAL finalcost=0,const REAL maxTime=-1);
      virtual void XMLOutput(ostream &os,int indent=0)const;
      virtual void XMLInput(istream &is,const XMLCrystTag &tag);
   private:
      /// Try a new configuration by expanding the worst vertex from the
      /// center by a factor f. If it is better, store it as new worst.
      /// Return the new obtained llk
      REAL GenerateNewSimplexConfiguration(CrystVector_REAL &vLLK,
                                           CrystVector_long &vIndex,
                                           unsigned long worst,
                                           REAL f);
   #ifdef __WX__CRYST__
   public:
      // :TODO: This should not be required !
      virtual WXCrystObjBasic* WXCreate(wxWindow*);
      virtual WXOptimizationObj* WXGet();
      virtual void WXDelete();
      virtual void WXNotifyDelete();
   #endif
};

}//namespace

#endif //_CONJUGATEGRADIENT_H
