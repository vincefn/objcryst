/*
* ObjCryst++ : a Crystallographic computing library in C++
*
*  (c) 2000 Vincent FAVRE-NICOLIN
*  vincefn@users.sourceforge.net
*
*/

#ifndef _VFN_WX_DIFFRACTIONSINGLECRYSTAL_H_
#define _VFN_WX_DIFFRACTIONSINGLECRYSTAL_H_

#include "ObjCryst/DiffractionDataSingleCrystal.h"
#include "wxCryst/wxRefinableObj.h"
namespace ObjCryst
{

/** WX Class for DiffractionDataSingleCrystal objects
*
*
*/
class WXDiffractionSingleCrystal: public WXRefinableObj
{
   public:
      WXDiffractionSingleCrystal(wxWindow *parent, DiffractionDataSingleCrystal*);
      virtual void OnUpdateUI(wxUpdateUIEvent& event);
   private:
      void OnMenuSimulate(wxCommandEvent & WXUNUSED(event));
      void OnMenuImport(wxCommandEvent & event);
      void OnMenuSaveHKLIobsIcalc(wxCommandEvent & WXUNUSED(event));
      void OnMenuSetWavelength(wxCommandEvent &event);
      void OnChangeCrystal(wxCommandEvent & WXUNUSED(event));
      WXFieldChoice* mpFieldCrystal;
      DiffractionDataSingleCrystal *mpData;
   DECLARE_EVENT_TABLE()
};

}//namespace
#endif
