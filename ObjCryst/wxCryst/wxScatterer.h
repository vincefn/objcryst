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
/*
*  header file for the RefinablePar and RefinableObj classes
*
* This is still in early development stages !! Not secure !
*
*/

#ifndef _VFN_WX_SCATTERER_OBJ_H_
#define _VFN_WX_SCATTERER_OBJ_H_

#include "wxCryst/wxRefinableObj.h"
#include "ObjCryst/Scatterer.h"
namespace ObjCryst
{

/// base wxCryst class for Scatterers
class WXScatterer: public WXRefinableObj
{
   public:
      WXScatterer(wxWindow *parent, Scatterer*);
   protected:
      Scatterer* mpScatterer;
      WXFieldRefPar* mpFieldX;
      WXFieldRefPar* mpFieldY;
      WXFieldRefPar* mpFieldZ;
      WXFieldRefPar* mpFieldPopu;
};

} //namespace

#endif //_VFN_WX_SCATTERER_OBJ_H_
