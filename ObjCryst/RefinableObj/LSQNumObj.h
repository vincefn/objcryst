/*
* LibCryst++ : a Crystallographic computing library in C++
*
*  (c) 2000 Vincent FAVRE-NICOLIN
*           Laboratoire de Cristallographie
*           24, quai Ernest-Ansermet, CH-1211 Geneva 4, Switzerland
*  Contact: Vincent.Favre-Nicolin@cryst.unige.ch
*           Radovan.Cerny@cryst.unige.ch
*
*/
/*   LSQObjNum.h
*  header file for Least-Squares refinement witn Numerical derivatives Objects
*
*/
#ifndef _LSQOBJNUM_H

#include "CrystVector/CrystVector.h"
#include "RefinableObj/RefinableObj.h"
#include <string>

using namespace std;
using namespace ObjCryst;
/** \brief (Quick & dirty) Least-Squares Refinement Object with Numerical derivatives
*
* This is is very quick& dirty (slow), written one day for tests... Do not expect anything
* from it yet... I'm not sure it currently works !
*/
class LSQObjNum
{
	public:
		LSQObjNum(CrystVector_double (*pRefinedFunc)(),
                                  string objName="Unnamed LSQ object");
		~LSQObjNum();
		void Init(CrystVector_double&obs,CrystVector_double&weight);
		void Init(CrystVector_double&obs);//using unit weight
      
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
      
		void Refine (int nbCycle=1,bool useLevenbergMarquardt=false);
		CrystVector_double Sigma()const;
		CrystMatrix_double CorrelMatrix()const;
		double Rfactor()const;
		double RwFactor()const;
		double ChiSquare()const;	//uses the weight if specified
      ///Add a list of refinable parameters
		void AddRefParList(RefinableObj &refParList);
		void SetUseSaveFileOnEachCycle(bool yesOrNo=true);
		void SetSaveFile(string fileName="refine.save");
		void PrintRefResults()const;
		void SetDampingFactor(const double newDampFact);
		void PurgeSaveFile();
		void WriteReportToFile()const;
      
		void OptimizeDerivativeSteps();
	protected:
	private:
      /// Prepare mRefParList for the refinement
      void PrepareRefParList()const;
      /// Array of pointers to the lists of refined parameters
		RefinableObj *mpRefParList[100];
      /// Number of lists of refined parameters
      int mNbRefParList;
      /// The refinable par list used during refinement. Only a condensed version
      /// of the refinable par list *mpRefParList[]
		mutable RefinableObj mRefParList;
		//The minimised function
		//	Obs and Calc values are stored into Obs and Calc arrays
		CrystVector_double (*mpRefinedFunc) ();
      /// Damping factor for the refinement (unused yet...)
		double mDampingFactor;
      ///Save result to file after each cycle ?
		bool mSaveReportOnEachCycle;	
      /// Name of the refined object
		string mName;
      /// File name where refinement info is saved
		string mSaveFileName;
		double mR,mRw,mChiSq;
      /// Correlation matrix between all refined parameters.
		CrystMatrix_double mCorrelMatrix;
      /// Observed values.
		CrystVector_double mObs;
      /// Weight corresponding to all observed values.
		CrystVector_double mWeight;
      /// Index of the set of saved values for all refinable parameters, before refinement
      /// and before the last cycle.
      int mIndexValuesSetInitial, mIndexValuesSetLast;
      /// If true, then stop at the end of the cycle. Used in multi-threading environment
      bool mStopAfterCycle;
};

#endif //_LSQOBJNUM_H
