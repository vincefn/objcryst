/*
* ObjCryst++ : a Crystallographic computing library in C++
*
*  (c) 2000 Vincent FAVRE-NICOLIN
*   vincefn@users.sourceforge.net
*
*/

#ifndef _VFN_WX_RADIATION_H_
#define _VFN_WX_RADIATION_H_

#include "wxCryst/wxRefinableObj.h"
#include "ObjCryst/ScatteringData.h"
namespace ObjCryst
{

/// WX Class for Radiation
class WXRadiation: public WXCrystObjBasic
{
   public:
      WXRadiation(wxWindow *parent, Radiation*);
      virtual void CrystUpdate();
   private:
      Radiation *mpRadiation;
      WXFieldOption *mpFieldRadType;
      WXFieldRefPar *mpFieldWavelength;
      wxBoxSizer *mpSizer;
   DECLARE_EVENT_TABLE()
};

}//namespace

#endif
