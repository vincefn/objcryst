//#include <sstream> //for stringstream
#include <fstream>

#include "wx/wx.h"
#include "wxCryst/wxAtom.h"
//Fixes for Cygwin; where do those stupid macros come from ? Somewhere in wxMSW headers
#ifdef max
#undef max
#endif
#ifdef min
#undef min
#endif
#ifdef DrawText
#undef DrawText
#endif
 
namespace ObjCryst
{
////////////////////////////////////////////////////////////////////////
//
//    WXAtom
//
////////////////////////////////////////////////////////////////////////
BEGIN_EVENT_TABLE(WXAtom,wxWindow)
   EVT_BUTTON(ID_ATOM_SCATTPOW,     WXAtom::OnChangeScattPow)
END_EVENT_TABLE()

WXAtom::WXAtom(wxWindow* parent, Atom *obj):
WXScatterer(parent,obj),mpAtom(obj)
{
   VFN_DEBUG_MESSAGE("WXAtom::WXAtom()",6)
      //mpFieldScattPower=new WXFieldName(this,"Symbol:",this,-1);
      mpFieldScattPower=new WXFieldChoice(this,ID_ATOM_SCATTPOW,"Scattering Power:");
      mpSizer->Add(mpFieldScattPower,0,wxALIGN_LEFT);
      mList.Add(mpFieldScattPower);
   
   this->CrystUpdate();
   this->Layout();
}

void WXAtom::OnChangeScattPow(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXAtom::OnChangeScattPow()",6)
   int choice;
   const ScatteringPowerAtom *scatt=dynamic_cast<const ScatteringPowerAtom*>
      ( WXDialogChooseFromRegistry
         (mpAtom->GetCrystal().GetScatteringPowerRegistry(),(wxWindow*)this,
         "Choose a new Scattering Power",choice));
   if(0==scatt) return;
   mpAtom->Init(mpAtom->GetX(),mpAtom->GetY(),mpAtom->GetZ(),mpAtom->GetName(),
                scatt,mpAtom->GetOccupancy());
   this->CrystUpdate();
}

void WXAtom::UpdateUI()
{
   mpFieldScattPower->SetValue(mpAtom->GetScatteringPower().GetName());
	this->WXRefinableObj::UpdateUI();
}


}// namespace 

