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

#include "wxCryst/wxScatteringPowerSphere.h"

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
//    WXScatteringPowerSphere
//
////////////////////////////////////////////////////////////////////////
static const long ID_SCATTPOWATOM_MENU_COLOUR=WXCRYST_ID();
static const long ID_SCATTPOWATOM_MENU_COLOUR_SETRGB=WXCRYST_ID();

BEGIN_EVENT_TABLE(WXScatteringPowerSphere, wxWindow)
   EVT_MENU(ID_SCATTPOWATOM_MENU_COLOUR_SETRGB, WXScatteringPowerSphere::OnChangeColour)
   EVT_UPDATE_UI(ID_CRYST_UPDATEUI,             WXRefinableObj::OnUpdateUI)
END_EVENT_TABLE()

WXScatteringPowerSphere::WXScatteringPowerSphere(wxWindow* parent, 
                                                       ScatteringPowerSphere *obj):
WXRefinableObj(parent,(RefinableObj*)obj),mpScatteringPower(obj)
{
   VFN_DEBUG_MESSAGE("WXScatteringPowerSphere::ScatteringPowerSphere()",6)
   mpWXTitle->SetForegroundColour(wxColour(0,200,0));
   //Lattice
      mpMenuBar->AddMenu("Colour",ID_SCATTPOWATOM_MENU_COLOUR);
         mpMenuBar->AddMenuItem(ID_SCATTPOWATOM_MENU_COLOUR,
                                ID_SCATTPOWATOM_MENU_COLOUR_SETRGB,"Set RGB Colour");
      WXCrystObjBasic* pFieldRadius=
         mpScatteringPower->GetPar(&(mpScatteringPower->mRadius)).WXCreate(this);
      WXCrystObjBasic* pFieldBiso=
         mpScatteringPower->GetPar(&(mpScatteringPower->mBiso  )).WXCreate(this);
      
      mpSizer->Add(pFieldRadius,0,wxALIGN_LEFT);
      mList.Add(pFieldRadius);
      mpSizer->Add(pFieldBiso,0,wxALIGN_LEFT);
      mList.Add(pFieldBiso);
   this->CrystUpdate(true);
   this->Layout();
}

bool WXScatteringPowerSphere::OnChangeName(const int id)
{
   VFN_DEBUG_MESSAGE("WXScatteringPowerSphere::OnChangeName()",6)
   if(this->WXRefinableObj::OnChangeName(id)==true) return true;
   return false;
}

void WXScatteringPowerSphere::OnChangeColour(wxCommandEvent & event)
{
   const float* oldColour=mpScatteringPower->GetColourRGB();
   double r,g,b;
   r=oldColour[0];
   g=oldColour[1];
   b=oldColour[2];
   //red
   {
      wxString str;
      str<<r;
      wxTextEntryDialog dialog(this,"Red",
                              "Enter Red component (0.<r<1.)",str,wxOK | wxCANCEL);
      if(wxID_OK!=dialog.ShowModal())
      {
         VFN_DEBUG_EXIT("WXScatteringPowerSphere::OnChangeColour():Cancelled",6)
         return;
      }
      dialog.GetValue().ToDouble(&r);
   }
   //green
   {
      wxString str;
      str<<g;
      wxTextEntryDialog dialog(this,"Green",
                              "Enter Green component (0.<g<1.)",str,wxOK | wxCANCEL);
      if(wxID_OK!=dialog.ShowModal())
      {
         VFN_DEBUG_EXIT("WXScatteringPowerSphere::OnChangeColour():Cancelled",6)
         return;
      }
      dialog.GetValue().ToDouble(&g);
   }
   //Blue
   {
      wxString str;
      str<<b;
      wxTextEntryDialog dialog(this,"Blue",
                              "Enter Blue component (0.<b<1.)",str,wxOK | wxCANCEL);
      if(wxID_OK!=dialog.ShowModal())
      {
         VFN_DEBUG_EXIT("WXScatteringPowerSphere::OnChangeColour():Cancelled",6)
         return;
      }
      dialog.GetValue().ToDouble(&b);
   }
   mpScatteringPower->SetColour(r,g,b);
}

void WXScatteringPowerSphere::UpdateUI(const bool lock)
{
   this->WXRefinableObj::UpdateUI(lock);
}

}// namespace 

