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

#ifndef _VFN_WX_SCATTERINGPOWER_OBJ_H_
#define _VFN_WX_SCATTERINGPOWER_OBJ_H_

#include "wxCryst/wxRefinableObj.h"
#include "ObjCryst/ScatteringPower.h"

namespace ObjCryst
{
/// wxCryst class for ScatteringPowerAtom
class WXScatteringPowerAtom: public WXRefinableObj
{
   public:
      WXScatteringPowerAtom(wxWindow *parent, ScatteringPowerAtom*);
      virtual void CrystUpdate();
      virtual bool OnChangeName(const int id);
      void OnChangeColour(wxCommandEvent & event);
   protected:
      ScatteringPowerAtom* mpScatteringPowerAtom;
      WXFieldName *mpFieldSymbol;
      WXFieldRefPar* mpFieldBiso;
   DECLARE_EVENT_TABLE()
};

} //namespace

#endif //_VFN_WX_SCATTERINGPOWER_OBJ_H_
