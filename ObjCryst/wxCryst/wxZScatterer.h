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

#ifndef _VFN_WX_ZSCATTERER_H_
#define _VFN_WX_ZSCATTERER_H_

#include "wxCryst/wxScatterer.h"
#include "ObjCryst/ZScatterer.h"

namespace ObjCryst
{

class WXZAtom:public WXCrystObjBasic
{
   public:
      WXZAtom(wxWindow *parent, ZAtom*);
      virtual void CrystUpdate();
      virtual bool Layout();
      void OnChangeScattPow(wxCommandEvent & WXUNUSED(event));
      void OnChangeName(wxCommandEvent & WXUNUSED(event));
      virtual void OnUpdateUI(wxUpdateUIEvent& event);
   private:
      ZAtom *mpZAtom;
      wxBoxSizer *mpSizer;
      WXCrystObjBasicList mList;
      wxTextCtrl *mpFieldName;
      WXFieldChoice* mpFieldScattPower;
      WXFieldRefPar* mpFieldBond;
      WXFieldRefPar* mpFieldAngle;
      WXFieldRefPar* mpFieldDihed;
   DECLARE_EVENT_TABLE()
};

class WXZScatterer: public WXScatterer
{
   public:
      WXZScatterer(wxWindow *parent, ZScatterer*);
      void OnMenuAddZAtom(wxCommandEvent & WXUNUSED(event));
      void OnMenuSetLimits(wxCommandEvent &event);
      void OnMenuChangePivotAtom(wxCommandEvent &WXUNUSED(event));
      void OnMenuImportZMatrix(wxCommandEvent &WXUNUSED(event));
      void OnMenuExportZMatrix(wxCommandEvent &WXUNUSED(event));
   private:
      ZScatterer* mpZScatterer;
      WXRegistry<ZAtom> *mpWXZAtomRegistry;
   DECLARE_EVENT_TABLE()
};

} //namespace

#endif //_VFN_WX_ZSCATTERER_H_
