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

#ifndef _VFN_WX_GLOBALOPTIM_OBJ_H_
#define _VFN_WX_GLOBALOPTIM_OBJ_H_

#include "wxCryst/wxCryst.h"
namespace ObjCryst
{
   class WXGlobalOptimObj;
}
#include "RefinableObj/GlobalOptimObj.h"
namespace ObjCryst
{
class WXGlobalOptimRunThread;

/// WX Class for a Global Optimization objects
class WXGlobalOptimObj: public WXCrystObj
{
   public:
      WXGlobalOptimObj(wxWindow *parent, GlobalOptimObj*);
      virtual void CrystUpdate();
      ///Update the number of trials to go, to show the user something's still going on...
      virtual void UpdateDisplayNbTrial();
      virtual bool OnChangeName(const int id);
      virtual void OnSave();
      virtual void OnLoad();
      ///From the menu
      virtual void OnAddRefinedObject();
      ///Added by the library
      virtual void AddRefinedObject(RefinableObj &obj);
      ///From the menu
      virtual void OnAddCostFunction();
      ///Added by the library
      virtual void AddCostFunction(RefinableObj &obj,const int costFuncNum);
      virtual void OnRunOptimization();
      virtual void OnStopOptimization();
   private:
      GlobalOptimObj *mpGlobalOptimObj;
      WXCrystMenuBar* mpMenuBar;
      WXGlobalOptimRunThread *mpGlobalOptimRunThread;
      /// The number of trials asked by the user
      long mNbTrial;
      WXFieldPar<long> *mpWXFieldNbTrial;
   DECLARE_EVENT_TABLE()
};

/// Class for a GlobalOPtimization thread
class WXGlobalOptimRunThread: public wxThread
{
   public:
      WXGlobalOptimRunThread(GlobalOptimObj *globalOptObj,long &nbTrial);
      
      virtual void *Entry();
      virtual void OnExit();
   private:
      GlobalOptimObj *mpGlobalOptObj;
      ///This points to the mNbTrial member in WXGlobalOptimObj
      long *mpNbTrial;
};


} //namespace

#endif //_VFN_WX_GLOBALOPTIM_OBJ_H_
