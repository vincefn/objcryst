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

#include "wx/glcanvas.h"
extern "C" {
#include "GL/glu.h"
}

#include "wxCryst/wxScatterer.h"
#include "ObjCryst/Atom.h"
namespace ObjCryst
{
/// wxCryst class for Atoms
class WXAtom: public WXScatterer
{
   public:
      WXAtom(wxWindow *parent, Atom*);
      virtual void CrystUpdate();
      void OnChangeScattPow(wxCommandEvent & WXUNUSED(event));
   private:
      Atom* mpAtom;
      WXFieldChoice* mpFieldScattPower;
   DECLARE_EVENT_TABLE()
};

} //namespace

#endif //_VFN_WX_ATOM_H_
