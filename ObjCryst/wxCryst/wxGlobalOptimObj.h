/*
* ObjCryst++ : a Crystallographic computing library in C++
*
*  (c) 2000 Vincent FAVRE-NICOLIN
*       vincefn@users.sourceforge.net
*
*/

#ifndef _VFN_WX_GLOBALOPTIM_OBJ_H_
#define _VFN_WX_GLOBALOPTIM_OBJ_H_

#include "wxCryst/wxCryst.h"
namespace ObjCryst
{
   class WXOptimizationObj;
	class WXMonteCarloObj;
}
#include "RefinableObj/GlobalOptimObj.h"
namespace ObjCryst
{
class WXGlobalOptimRunThread;

/// WX Class for a Global Optimization objects
class WXOptimizationObj: public WXCrystObj
{
   public:
      WXOptimizationObj(wxWindow *parent, OptimizationObj*);
      virtual void CrystUpdate();
      virtual bool OnChangeName(const int id);
      virtual void OnSave();
      virtual void OnLoad();
      /// From the menu
      virtual void OnAddRefinedObject();
      /// Added by the library
      virtual void AddRefinedObject(RefinableObj &obj);
      /// From the menu
      virtual void OnAddCostFunction();
      /// Added by the library
      virtual void AddCostFunction(RefinableObj &obj,const int costFuncNum);
		/// Launches the optimization run
      virtual void OnRunOptimization()=0;
      virtual void OnStopOptimization();
		virtual OptimizationObj & GetOptimizationObj()=0;
		virtual const OptimizationObj & GetOptimizationObj()const=0;
   protected:
      WXCrystMenuBar* mpMenuBar;
      WXGlobalOptimRunThread *mpGlobalOptimRunThread;
      WXFieldPar<long> *mpWXFieldNbTrial;
};

/// Class for a GlobalOptimization thread
class WXGlobalOptimRunThread: public wxThread
{
   public:
      WXGlobalOptimRunThread(OptimizationObj &globalOptObj,long &nbTrial,const double finalCost=0);
      
      virtual void *Entry();
      virtual void OnExit();
   private:
      OptimizationObj *mpGlobalOptObj;
      ///This points to the mNbTrial member in WXOptimizationObj
      long *mpNbTrial;
		/// The value of the cost below which the optimization should stop
		/// (0 by default) even if the desired number pf trial has not been reached.
		const double mFinalCost;
};

/** Class for Graphical interface to Monte-Carlo objects 
* (Simulated Annealing, Parallel Tempering)
*
*/
class WXMonteCarloObj: public WXOptimizationObj
{
   public:
      WXMonteCarloObj(wxWindow *parent, MonteCarloObj*);
      //virtual void CrystUpdate();
      virtual void OnRunOptimization();
		/// Called during optimization, to show the user something's still going on...
      void UpdateDisplayNbTrial();
		virtual OptimizationObj & GetOptimizationObj();
		virtual const OptimizationObj & GetOptimizationObj()const;
   protected:
		/// The algorithm object
      MonteCarloObj *mpMonteCarloObj;
      /// The number of trials asked by the user
      long mNbTrial;
      WXFieldPar<long> *mpWXFieldNbTrial;
   DECLARE_EVENT_TABLE()
};


} //namespace

#endif //_VFN_WX_GLOBALOPTIM_OBJ_H_
