//#include <sstream> //for stringstream
#include <fstream>

#include "wx/wx.h"
#include "wxCryst/wxScatteringPower.h"

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
BEGIN_EVENT_TABLE(WXScatteringPowerAtom, wxWindow)
   EVT_MENU(ID_SCATTPOWATOM_MENU_COLOUR_SETRGB, WXScatteringPowerAtom::OnChangeColour)
   EVT_UPDATE_UI(ID_CRYST_UPDATEUI, 				WXScatteringPowerAtom::OnUpdateUI)
END_EVENT_TABLE()

WXScatteringPowerAtom::WXScatteringPowerAtom(wxWindow* parent, ScatteringPowerAtom *obj):
WXRefinableObj(parent,(RefinableObj*)obj),mpScatteringPowerAtom(obj)
{
   VFN_DEBUG_MESSAGE("WXScatteringPowerAtom::WXScatteringPowerAtom()",6)
   mpWXTitle->SetForegroundColour(wxColour(0,200,0));
   //mpScatteringPowerAtom->Print();
   //mpScatteringPowerAtom->RefinableObj::Print();
   //cout << mpScatteringPowerAtom->GetSymbol()<<endl;
   //Lattice
      wxBoxSizer* sizer=new wxBoxSizer(wxHORIZONTAL);
      mpMenuBar->AddMenu("Colour",ID_SCATTPOWATOM_MENU_COLOUR);
         mpMenuBar->AddMenuItem(ID_SCATTPOWATOM_MENU_COLOUR,
                                ID_SCATTPOWATOM_MENU_COLOUR_SETRGB,"Set RGB Colour");
   
      mpFieldSymbol=new WXFieldName
         (this,"Symbol:",this,ID_WXSCATTPOWATOM_SYMBOL);

      mpFieldBiso  =new WXFieldRefPar(this,"Biso:",
            &(mpScatteringPowerAtom->GetPar(&(mpScatteringPowerAtom->mBiso))) );
            
      
      sizer->Add(mpFieldSymbol  ,0,wxALIGN_CENTER);
      sizer->Add(mpFieldBiso    ,0,wxALIGN_CENTER);
      
      mpSizer->Add(sizer,0,wxALIGN_LEFT);
      mList.Add(mpFieldSymbol);
      mList.Add(mpFieldBiso);
      this->CrystUpdate();
   this->Layout();
}

bool WXScatteringPowerAtom::OnChangeName(const int id)
{
   VFN_DEBUG_MESSAGE("WXScatteringPowerAtom::OnChangeName()",6)
   if(this->WXRefinableObj::OnChangeName(id)==true) return true;
   switch(id)
   {
      case ID_WXSCATTPOWATOM_SYMBOL:
      {
         VFN_DEBUG_MESSAGE("WXScatteringPowerAtom::OnChangeName():Changing Symbol",6)
         mpScatteringPowerAtom->Init(mpScatteringPowerAtom->GetName(),
                                     mpFieldSymbol->GetValue(),
                                     mpScatteringPowerAtom->GetBiso());
         return true;
      }
   }
   return false;
}

void WXScatteringPowerAtom::OnChangeColour(wxCommandEvent & event)
{
   const float* oldColour=mpScatteringPowerAtom->GetColourRGB();
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
         VFN_DEBUG_EXIT("WXScatteringPowerAtom::OnChangeColour():Cancelled",6)
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
         VFN_DEBUG_EXIT("WXScatteringPowerAtom::OnChangeColour():Cancelled",6)
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
         VFN_DEBUG_EXIT("WXScatteringPowerAtom::OnChangeColour():Cancelled",6)
         return;
      }
      dialog.GetValue().ToDouble(&b);
   }
   mpScatteringPowerAtom->SetColour(r,g,b);
}
void WXScatteringPowerAtom::OnUpdateUI(wxUpdateUIEvent& event)
{
   mpFieldSymbol->SetValue(mpScatteringPowerAtom->GetSymbol());
	this->WXRefinableObj::OnUpdateUI(event);
}

}// namespace 

