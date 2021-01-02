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
//#include <sstream> //for stringstream
#include <fstream>

// wx headers, with or without precompilation
#include "wx/wxprec.h"
#ifdef __BORLANDC__
    #pragma hdrstop
#endif
#ifndef WX_PRECOMP
    #include "wx/wx.h"
#endif

#include "ObjCryst/wxCryst/wxAtom.h"
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
static const long ID_ATOM_SCATTPOW=WXCRYST_ID();

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

   this->BottomLayout(0);
   this->CrystUpdate(true,true);
}
WXAtom::~WXAtom()
{
   mpAtom->WXNotifyDelete();
}

void WXAtom::OnChangeScattPow(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXAtom::OnChangeScattPow()",6)
   WXCrystValidateAllUserInput();
   int choice;
   const ScatteringPower *scatt=
       WXDialogChooseFromRegistry
         (mpAtom->GetCrystal().GetScatteringPowerRegistry(),(wxWindow*)this,
         "Choose a new Scattering Power",choice);
   if(0==scatt) return;
   mpAtom->Init(mpAtom->GetX(),mpAtom->GetY(),mpAtom->GetZ(),mpAtom->GetName(),
                scatt,mpAtom->GetOccupancy());
   this->CrystUpdate(true,true);
}

void WXAtom::UpdateUI(const bool lock)
{
   mpFieldScattPower->SetValue(mpAtom->GetScatteringPower().GetName());
   this->WXRefinableObj::UpdateUI();
}


}// namespace
