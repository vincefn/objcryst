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

#include "wx/wx.h"
#include "wxCryst/wxScatteringPowerFullerene.h"

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
//    WXScatteringPowerFullerene
//
////////////////////////////////////////////////////////////////////////
BEGIN_EVENT_TABLE(WXScatteringPowerFullerene, wxWindow)
   EVT_MENU(ID_SCATTPOWATOM_MENU_COLOUR_SETRGB, WXScatteringPowerFullerene::OnChangeColour)
   EVT_UPDATE_UI(ID_CRYST_UPDATEUI,             WXRefinableObj::OnUpdateUI)
END_EVENT_TABLE()

WXScatteringPowerFullerene::WXScatteringPowerFullerene(wxWindow* parent, 
                                                       ScatteringPowerFullerene *obj):
WXRefinableObj(parent,(RefinableObj*)obj),mpScatteringPower(obj)
{
   VFN_DEBUG_MESSAGE("WXScatteringPowerFullerene::WXScatteringPowerAtom()",6)
   mpWXTitle->SetForegroundColour(wxColour(0,200,0));
   //Lattice
      mpMenuBar->AddMenu("Colour",ID_SCATTPOWATOM_MENU_COLOUR);
         mpMenuBar->AddMenuItem(ID_SCATTPOWATOM_MENU_COLOUR,
                                ID_SCATTPOWATOM_MENU_COLOUR_SETRGB,"Set RGB Colour");
      #if 1
      WXFieldRefPar* pFieldRadiusX  =new WXFieldRefPar(this,"Radius:",
            &(mpScatteringPower->GetPar(&(mpScatteringPower->mAxisLengthX))) );
      #else //Anisotropic
      WXFieldRefPar* pFieldRadiusX  =new WXFieldRefPar(this,"XRadius:",
            &(mpScatteringPower->GetPar(&(mpScatteringPower->mAxisLengthX))) );
      WXFieldRefPar* pFieldRadiusY  =new WXFieldRefPar(this,"YRadius:",
            &(mpScatteringPower->GetPar(&(mpScatteringPower->mAxisLengthY))) );
      WXFieldRefPar* pFieldRadiusZ  =new WXFieldRefPar(this,"ZRadius:",
            &(mpScatteringPower->GetPar(&(mpScatteringPower->mAxisLengthZ))) );
      WXFieldRefPar* pFieldPhi  =new WXFieldRefPar(this,"Phi:",
            &(mpScatteringPower->GetPar(&(mpScatteringPower->mPhi))) );
      WXFieldRefPar* pFieldChi  =new WXFieldRefPar(this,"Chi:",
            &(mpScatteringPower->GetPar(&(mpScatteringPower->mChi))) );
      WXFieldRefPar* pFieldPsi  =new WXFieldRefPar(this,"Psi:",
            &(mpScatteringPower->GetPar(&(mpScatteringPower->mPsi))) );
      #endif
      WXFieldRefPar* pFieldNbAtom  =new WXFieldRefPar(this,"Number of atoms:",
            &(mpScatteringPower->GetPar(&(mpScatteringPower->mNbAtom))) );
      WXFieldRefPar* pFieldBiso  =new WXFieldRefPar(this,"Biso:",
            &(mpScatteringPower->mpScatteringPower->GetPar(&(mpScatteringPower->mpScatteringPower->GetBiso()))) );
      
      mpSizer->Add(pFieldRadiusX,0,wxALIGN_LEFT);
      mList.Add(pFieldRadiusX);
      #if 0
      mpSizer->Add(pFieldRadiusY,0,wxALIGN_LEFT);
      mList.Add(pFieldRadiusY);
      mpSizer->Add(pFieldRadiusZ,0,wxALIGN_LEFT);
      mList.Add(pFieldRadiusZ);
      mpSizer->Add(pFieldPhi,0,wxALIGN_LEFT);
      mList.Add(pFieldPhi);
      mpSizer->Add(pFieldChi,0,wxALIGN_LEFT);
      mList.Add(pFieldChi);
      mpSizer->Add(pFieldPsi,0,wxALIGN_LEFT);
      mList.Add(pFieldPsi);
      #endif
      mpSizer->Add(pFieldNbAtom,0,wxALIGN_LEFT);
      mList.Add(pFieldNbAtom);
      mpSizer->Add(pFieldBiso,0,wxALIGN_LEFT);
      mList.Add(pFieldBiso);
      this->CrystUpdate();
   this->Layout();
}

bool WXScatteringPowerFullerene::OnChangeName(const int id)
{
   VFN_DEBUG_MESSAGE("WXScatteringPowerFullerene::OnChangeName()",6)
   if(this->WXRefinableObj::OnChangeName(id)==true) return true;
   return false;
}

void WXScatteringPowerFullerene::OnChangeColour(wxCommandEvent & event)
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
         VFN_DEBUG_EXIT("WXScatteringPowerFullerene::OnChangeColour():Cancelled",6)
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
         VFN_DEBUG_EXIT("WXScatteringPowerFullerene::OnChangeColour():Cancelled",6)
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
         VFN_DEBUG_EXIT("WXScatteringPowerFullerene::OnChangeColour():Cancelled",6)
         return;
      }
      dialog.GetValue().ToDouble(&b);
   }
   mpScatteringPower->SetColour(r,g,b);
}

void WXScatteringPowerFullerene::UpdateUI()
{
   this->WXRefinableObj::UpdateUI();
}

}// namespace 

