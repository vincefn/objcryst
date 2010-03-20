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

#include "ObjCryst/wxCryst/wxScatteringPower.h"

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
//    WXScatteringPowerAtom
//
////////////////////////////////////////////////////////////////////////
static const long ID_WXSCATTPOWATOM_SYMBOL          =WXCRYST_ID();
static const long ID_SCATTPOWATOM_MENU_COLOUR       =WXCRYST_ID();
static const long ID_SCATTPOWATOM_MENU_COLOUR_SETRGB=WXCRYST_ID();

BEGIN_EVENT_TABLE(WXScatteringPowerAtom, wxWindow)
   EVT_MENU(ID_SCATTPOWATOM_MENU_COLOUR_SETRGB, WXScatteringPowerAtom::OnChangeColour)
   EVT_UPDATE_UI(ID_CRYST_UPDATEUI,             WXRefinableObj::OnUpdateUI)
END_EVENT_TABLE()

WXScatteringPowerAtom::WXScatteringPowerAtom(wxWindow* parent, ScatteringPowerAtom *obj):
WXRefinableObj(parent,(RefinableObj*)obj),mpScatteringPowerAtom(obj)
{
   VFN_DEBUG_MESSAGE("WXScatteringPowerAtom::WXScatteringPowerAtom()",6)
   mpWXTitle->SetForegroundColour(wxColour(0,200,0));
   //Lattice
      mpFieldSymbol=new WXFieldName
         (this,"Symbol:",this,ID_WXSCATTPOWATOM_SYMBOL);

      mpTopSizer->Add(mpFieldSymbol  ,0,wxALIGN_TOP);
      mList.Add(mpFieldSymbol);
      
   this->BottomLayout(0);
   this->CrystUpdate(true);
}

bool WXScatteringPowerAtom::OnChangeName(const int id)
{
   VFN_DEBUG_MESSAGE("WXScatteringPowerAtom::OnChangeName()",6)
   if(this->WXRefinableObj::OnChangeName(id)==true) return true;
   if(id==ID_WXSCATTPOWATOM_SYMBOL)
   {
      VFN_DEBUG_MESSAGE("WXScatteringPowerAtom::OnChangeName():Changing Symbol",6)
      mpScatteringPowerAtom->Init(mpScatteringPowerAtom->GetName(),
                                  mpFieldSymbol->GetValue(),
                                  mpScatteringPowerAtom->GetBiso());
      return true;
   }
   return false;
}

void WXScatteringPowerAtom::OnChangeColour(wxCommandEvent & event)
{
   VFN_DEBUG_ENTRY("WXScatteringPowerAtom::OnChangeColour()",6)
   const float* oldColour=mpScatteringPowerAtom->GetColourRGB();
   double r,g,b;
   r=oldColour[0];
   g=oldColour[1];
   b=oldColour[2];
   //red
   {
      wxString str;
      str<<r;
      wxTextEntryDialog dialog(this,_T("Red"),
                              _T("Enter Red component (0.<r<1.)"),str,wxOK | wxCANCEL);
      if(wxID_OK!=dialog.ShowModal())
      {
         VFN_DEBUG_EXIT("WXScatteringPowerAtom::OnChangeColour():Cancelled",6)
         return;
      }
      dialog.GetValue().ToDouble(&r);
   }
   //green
   {
      wxString str;
      str<<g;
      wxTextEntryDialog dialog(this,_T("Green"),
                              _T("Enter Green component (0.<g<1.)"),str,wxOK | wxCANCEL);
      if(wxID_OK!=dialog.ShowModal())
      {
         VFN_DEBUG_EXIT("WXScatteringPowerAtom::OnChangeColour():Cancelled",6)
         return;
      }
      dialog.GetValue().ToDouble(&g);
   }
   //Blue
   {
      wxString str;
      str<<b;
      wxTextEntryDialog dialog(this,_T("Blue"),
                              _T("Enter Blue component (0.<b<1.)"),str,wxOK | wxCANCEL);
      if(wxID_OK!=dialog.ShowModal())
      {
         VFN_DEBUG_EXIT("WXScatteringPowerAtom::OnChangeColour():Cancelled",6)
         return;
      }
      dialog.GetValue().ToDouble(&b);
   }
   mpScatteringPowerAtom->SetColour(r,g,b);
   VFN_DEBUG_EXIT("WXScatteringPowerAtom::OnChangeColour()",6)
}
void WXScatteringPowerAtom::UpdateUI(const bool lock)
{
   VFN_DEBUG_ENTRY("WXScatteringPowerAtom::UpdateUI()",3)
   if(lock) mMutex.Lock();
   mpFieldSymbol->SetValue(mpScatteringPowerAtom->GetSymbol());
   this->WXRefinableObj::UpdateUI(false);
   if(lock) mMutex.Unlock();
   VFN_DEBUG_EXIT("WXScatteringPowerAtom::UpdateUI()",3)
}

}// namespace 

