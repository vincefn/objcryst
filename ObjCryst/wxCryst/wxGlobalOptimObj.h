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
      virtual void OnRemoveRefinedObject();
      /// From the menu
      virtual void OnAddCostFunction();
      /// Added by the library
      virtual void AddCostFunction(RefinableObj &obj,const int costFuncNum);
      /// From the menu
      virtual void OnRemoveCostFunction();
      /// Launches the optimization run
      virtual void OnRunOptimization()=0;
      virtual void OnStopOptimization();
      virtual OptimizationObj & GetOptimizationObj()=0;
      virtual const OptimizationObj & GetOptimizationObj()const=0;
      virtual void OnUpdateUI(wxUpdateUIEvent& event);
   protected:
      WXCrystMenuBar* mpMenuBar;
      WXGlobalOptimRunThread *mpGlobalOptimRunThread;
      WXFieldPar<long> *mpWXFieldNbTrial;
};

/// Class for a GlobalOptimization thread
class WXGlobalOptimRunThread: public wxThread
{
   public:
      WXGlobalOptimRunThread(OptimizationObj &globalOptObj,long &nbTrial,const REAL finalCost=0);
      
      virtual void *Entry();
      virtual void OnExit();
   private:
      OptimizationObj *mpGlobalOptObj;
      ///This points to the mNbTrial member in WXOptimizationObj
      long *mpNbTrial;
      /// The value of the cost below which the optimization should stop
      /// (0 by default) even if the desired number pf trial has not been reached.
      const REAL mFinalCost;
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
