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
/*   LSQNumObj.h
*  header file for Least-Squares refinement witn Numerical derivatives Objects
*
*/
#ifndef _LSQOBJNUM_H
#define _LSQOBJNUM_H

#include "CrystVector/CrystVector.h"
#include "RefinableObj/RefinableObj.h"
#include <string>

using namespace std;
using namespace ObjCryst;
/** \brief (Quick & dirty) Least-Squares Refinement Object with Numerical derivatives
*
* This is still highly experimental !
*/
class LSQNumObj
{
   public:
      LSQNumObj(string objName="Unnamed LSQ object");
      ~LSQNumObj();
      /// Fix one parameter
      void SetParIsFixed(const string& parName,const bool fix);
      /// Fix one family of parameters
      void SetParIsFixed(const RefParType *type,const bool fix);
      /// UnFix All parameters
      void UnFixAllPar();
      /// Set a parameter to be used
      void SetParIsUsed(const string& parName,const bool use);
      /// Set a family of parameters to be used
      void SetParIsUsed(const RefParType *type,const bool use);
      
      void Refine (int nbCycle=1,bool useLevenbergMarquardt=false,
                   const bool silent=false);
      CrystVector_REAL Sigma()const;
      CrystMatrix_REAL CorrelMatrix()const;
      REAL Rfactor()const;
      REAL RwFactor()const;
      REAL ChiSquare()const;   //uses the weight if specified
      ///Add an object to refine
      void SetRefinedObj(RefinableObj &obj, const unsigned int LSQFuncIndex=0);
      void SetUseSaveFileOnEachCycle(bool yesOrNo=true);
      void SetSaveFile(string fileName="refine.save");
      void PrintRefResults()const;
      void SetDampingFactor(const REAL newDampFact);
      void PurgeSaveFile();
      void WriteReportToFile()const;
      
      void OptimizeDerivativeSteps();
   protected:
   private:
      /// Prepare mRefParList for the refinement
      void PrepareRefParList();
      // Refined object
         /// The recursive list of all refined sub-objects
         ObjRegistry<RefinableObj> mRecursiveRefinedObjList;
      /// The refinable par list used during refinement. Only a compilation
      /// of the parameters in RefinableObj and its sub-objects
      mutable RefinableObj mRefParList;
      /// Damping factor for the refinement (unused yet...)
      REAL mDampingFactor;
      ///Save result to file after each cycle ?
      bool mSaveReportOnEachCycle;   
      /// Name of the refined object
      string mName;
      /// File name where refinement info is saved
      string mSaveFileName;
      REAL mR,mRw,mChiSq;
      /// Correlation matrix between all refined parameters.
      CrystMatrix_REAL mCorrelMatrix;
      /// Observed values.
      CrystVector_REAL mObs;
      /// Weight corresponding to all observed values.
      CrystVector_REAL mWeight;
      /// Index of the set of saved values for all refinable parameters, before refinement
      /// and before the last cycle.
      int mIndexValuesSetInitial, mIndexValuesSetLast;
      /// If true, then stop at the end of the cycle. Used in multi-threading environment
      bool mStopAfterCycle;
      /// The opitimized object
      RefinableObj *mpRefinedObj;
      /// The index of the LSQ function in the refined object (if there are several...)
      unsigned int mLSQFuncIndex;
};

#endif //_LSQOBJNUM_H
