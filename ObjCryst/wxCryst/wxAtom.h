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

#ifndef _VFN_WX_ATOM_H_
#define _VFN_WX_ATOM_H_

#include "wxCryst/wxScatterer.h"
#include "ObjCryst/Atom.h"
namespace ObjCryst
{
/// wxCryst class for Atoms
class WXAtom: public WXScatterer
{
   public:
      WXAtom(wxWindow *parent, Atom*);
      void OnChangeScattPow(wxCommandEvent & WXUNUSED(event));
      virtual void UpdateUI();
   private:
      Atom* mpAtom;
      WXFieldChoice* mpFieldScattPower;
   DECLARE_EVENT_TABLE()
};

} //namespace

#endif //_VFN_WX_ATOM_H_
