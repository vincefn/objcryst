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

// wx headers, with or without precompilation
#include "wx/wxprec.h"
#ifdef __BORLANDC__
    #pragma hdrstop
#endif
#ifndef WX_PRECOMP
    #include "wx/wx.h"
#endif

#include <sstream>
#include <algorithm>
#include "wxCryst/wxMolecule.h"
namespace ObjCryst
{

//:TODO: Move this to wxCryst.h
template<class T> T const* WXDialogChooseFromVector(const vector<T*> &reg,wxWindow*parent,
                                                const string &message,int &choice)
{
   wxString *choices = new wxString[reg.size()];
   for(unsigned int i=0;i<reg.size();i++) 
      choices[i]=(reg[i]->GetName()).c_str();
   wxSingleChoiceDialog dialog
         (parent,message.c_str(),"Choose",reg.size(),choices,0,wxOK | wxCANCEL);
   dialog.SetSize(300,300);
   if(wxID_OK!=dialog.ShowModal()) return 0;
   choice=dialog.GetSelection();
   delete[] choices;
   return reg[choice];
}

template<class T> T * WXDialogChooseFromVector(vector<T*> &reg,wxWindow*parent,
                                                const string &message,int &choice)
{
   wxString *choices = new wxString[reg.size()];
   for(unsigned int i=0;i<reg.size();i++) 
      choices[i]=(reg[i]->GetName()).c_str();
   wxSingleChoiceDialog dialog
         (parent,message.c_str(),"Choose",reg.size(),choices,0,wxOK | wxCANCEL);
   dialog.SetSize(300,300);
   if(wxID_OK!=dialog.ShowModal()) return 0;
   choice=dialog.GetSelection();
   delete[] choices;
   return reg[choice];
}

////////////////////////////////////////////////////////////////////////
//
//    WXMolAtom
//
////////////////////////////////////////////////////////////////////////
static const long ID_MOLATOM_SCATTPOW=WXCRYST_ID();
static const long ID_MOLATOM_NAME=WXCRYST_ID();

BEGIN_EVENT_TABLE(WXMolAtom,wxWindow)
   EVT_BUTTON(ID_MOLATOM_SCATTPOW,    WXMolAtom::OnChangeScattPow)
END_EVENT_TABLE()

WXMolAtom::WXMolAtom(wxWindow *parent, MolAtom*obj):
WXCrystObjBasic(parent),mpMolAtom(obj)
{
   VFN_DEBUG_ENTRY("WXMolAtom::WXMolAtom()",6)
   mpSizer=new wxBoxSizer(wxHORIZONTAL);
   wxStaticText* label=new wxStaticText(this,-1,"Atom");
   mpSizer->Add(label);
   {
      mpFieldName=new WXFieldString(this, mpMolAtom->GetName(),ID_MOLATOM_NAME,80,true);
      mpSizer->Add(mpFieldName,0,wxALIGN_CENTER);
      mpFieldScattPower=new WXFieldChoice(this,ID_MOLATOM_SCATTPOW,"Type:",60);
      mpSizer->Add(mpFieldScattPower,0,wxALIGN_CENTER);
      mList.Add(mpFieldScattPower);
      if(mpMolAtom->IsDummy())
         mpFieldScattPower->SetValue("Dummy");
      else
         mpFieldScattPower->SetValue(mpMolAtom->GetScatteringPower().GetName());
   }
   {
      WXFieldRefPar* pFieldX  =new WXFieldRefPar(this,"x",
                               &(mpMolAtom->GetMolecule().GetPar(&(mpMolAtom->X()))) );
      mpSizer->Add(pFieldX,0,wxALIGN_CENTER);
      mList.Add(pFieldX);
   }
   {
      WXFieldRefPar* pFieldY  =new WXFieldRefPar(this,"y",
                               &(mpMolAtom->GetMolecule().GetPar(&(mpMolAtom->Y()))) );
      mpSizer->Add(pFieldY,0,wxALIGN_CENTER);
      mList.Add(pFieldY);
   }
   {
      WXFieldRefPar* pFieldZ  =new WXFieldRefPar(this,"z",
                               &(mpMolAtom->GetMolecule().GetPar(&(mpMolAtom->Z()))) );
      mpSizer->Add(pFieldZ,0,wxALIGN_CENTER);
      mList.Add(pFieldZ);
   }
   
   this->SetSizer(mpSizer);
   this->CrystUpdate();
   this->Layout();
   VFN_DEBUG_EXIT("WXMolAtom::WXMolAtom()",6)
}

void WXMolAtom::CrystUpdate()
{
   VFN_DEBUG_ENTRY("WXMolAtom::CrystUpdate()",5)
   mList.CrystUpdate();
   VFN_DEBUG_EXIT("WXMolAtom::CrystUpdate()",5)
}

void WXMolAtom::UpdateUI()
{
   VFN_DEBUG_ENTRY("WXMolAtom::UpdateUI()",5)
   mList.UpdateUI();
   mpFieldName->SetValue(mpMolAtom->GetName().c_str());
   if(mpMolAtom->IsDummy())
      mpFieldScattPower->SetValue("Dummy");
   else
      mpFieldScattPower->SetValue(mpMolAtom->GetScatteringPower().GetName());
   VFN_DEBUG_EXIT("WXMolAtom::UpdateUI()",5)
}

bool WXMolAtom::Layout()
{
   VFN_DEBUG_ENTRY("WXMolAtom::Layout()",3)
   //:TODO: Cleeeeaaaaannnn thiiiisss !!!
   for(unsigned int i=0;i<mList.GetNb();i++)
      mpSizer->SetItemMinSize(mList.Get(i),
                              mList.Get(i)->GetSize().GetWidth(),
                              mList.Get(i)->GetSize().GetHeight());
      
   mpSizer->Layout();
   mpSizer->Fit(this);
   wxSizer* s=mWXParent->GetSizer();
   if(s != 0)
   {// Need to do it that way, in case  the parent is not a WXCrystObj
    // with an adequate Layout() function
      s->SetItemMinSize(this,this->GetSize().GetWidth(),this->GetSize().GetHeight());
      s->Fit(mWXParent);
   }
   mWXParent->Layout();
   VFN_DEBUG_EXIT("WXMolAtom::Layout()",3)
   return this->wxWindow::Layout();
}
void WXMolAtom::OnChangeScattPow(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_ENTRY("WXMolAtom::OnChangeScattPow()",6)
   WXCrystValidateAllUserInput();
   int choice;
   const ScatteringPower *scatt=WXDialogChooseFromRegistry(
               mpMolAtom->GetMolecule().GetCrystal().GetScatteringPowerRegistry(),
               (wxWindow*)this,"Choose a new Scattering Power",choice);
   if(0==scatt) return;
   mpMolAtom->SetScatteringPower(*scatt);
   this->CrystUpdate();
   this->UpdateUI();
   VFN_DEBUG_EXIT("WXMolAtom::OnChangeScattPow()",6)
}
////////////////////////////////////////////////////////////////////////
//
//    WXMolBond
//
////////////////////////////////////////////////////////////////////////
WXCRYST_ID ID_MOLBOND_ATOM1;
WXCRYST_ID ID_MOLBOND_ATOM2;

BEGIN_EVENT_TABLE(WXMolBond,wxWindow)
   EVT_BUTTON(ID_MOLBOND_ATOM1,    WXMolBond::OnChangeAtom)
   EVT_BUTTON(ID_MOLBOND_ATOM2,    WXMolBond::OnChangeAtom)
END_EVENT_TABLE()

WXMolBond::WXMolBond(wxWindow *parent, MolBond*obj):
WXCrystObjBasic(parent),mpMolBond(obj)
{
   VFN_DEBUG_ENTRY("WXMolBond::WXMolBond()",6)
   mpSizer=new wxBoxSizer(wxHORIZONTAL);

   mpFieldAtom1=new WXFieldChoice(this,ID_MOLBOND_ATOM1,"Bond:",60);
   mpFieldAtom1->SetValue(mpMolBond->GetAtom1().GetName());
   mpSizer->Add(mpFieldAtom1,0,wxALIGN_CENTER);
   mList.Add(mpFieldAtom1);
   
   mpFieldAtom2=new WXFieldChoice(this,ID_MOLBOND_ATOM2,"-",60);
   mpFieldAtom2->SetValue(mpMolBond->GetAtom2().GetName());
   mpSizer->Add(mpFieldAtom2,0,wxALIGN_CENTER);
   mList.Add(mpFieldAtom2);

   WXFieldPar<REAL> *length=
      new WXFieldPar<REAL>(this,"Length:",-1,&(mpMolBond->Length0()));
   mpSizer->Add(length,0,wxALIGN_CENTER);
   mList.Add(length);

   this->SetSizer(mpSizer);
   this->CrystUpdate();
   this->Layout();
   VFN_DEBUG_EXIT("WXMolBond::WXMolBond()",6)
}

void WXMolBond::CrystUpdate()
{
   VFN_DEBUG_ENTRY("WXMolBond::CrystUpdate()",5)
   mList.CrystUpdate();
   VFN_DEBUG_EXIT("WXMolBond::CrystUpdate()",5)
}

void WXMolBond::UpdateUI()
{
   VFN_DEBUG_ENTRY("WXMolBond::UpdateUI()",5)
   mList.UpdateUI();
   mpFieldAtom1->SetValue(mpMolBond->GetAtom1().GetName().c_str());
   mpFieldAtom2->SetValue(mpMolBond->GetAtom2().GetName().c_str());
   VFN_DEBUG_EXIT("WXMolBond::UpdateUI()",5)
}

bool WXMolBond::Layout()
{
   VFN_DEBUG_ENTRY("WXMolBond::Layout()",5)
   //:TODO: Cleeeeaaaaannnn thiiiisss !!!
   for(unsigned int i=0;i<mList.GetNb();i++)
      mpSizer->SetItemMinSize(mList.Get(i),
                              mList.Get(i)->GetSize().GetWidth(),
                              mList.Get(i)->GetSize().GetHeight());
      
   mpSizer->Layout();
   mpSizer->Fit(this);
   wxSizer* s=mWXParent->GetSizer();
   if(s != 0)
   {// Need to do it that way, in case  the parent is not a WXCrystObj
    // with an adequate Layout() function
      s->SetItemMinSize(this,this->GetSize().GetWidth(),this->GetSize().GetHeight());
      s->Fit(mWXParent);
   }
   mWXParent->Layout();
   VFN_DEBUG_EXIT("WXMolBond::Layout()",5)
   return this->wxWindow::Layout();
}

void WXMolBond::OnChangeAtom(wxCommandEvent &event)
{
   VFN_DEBUG_ENTRY("WXMolBond::OnChangeAtom()",6)
   WXCrystValidateAllUserInput();
   int choice;
   const MolAtom *at=WXDialogChooseFromVector(mpMolBond->GetMolecule().GetAtomList(),
                               (wxWindow*)this,"Choose a new Atom",choice);
   if(0==at) return;
   if((long)ID_MOLBOND_ATOM1==event.GetId())
   {
      if(at==&(mpMolBond->GetAtom2()))
      {
         wxMessageDialog dumbUser(this,"The two atoms must be different !",
                                  "Whooops",wxOK|wxICON_EXCLAMATION);
         dumbUser.ShowModal();
         return;
      }
      mpMolBond->SetAtom1(*at);
   }
   if((long)ID_MOLBOND_ATOM2==event.GetId())
   {
      if(at==&(mpMolBond->GetAtom1()))
      {
         wxMessageDialog dumbUser(this,"The two atoms must be different !",
                                  "Whooops",wxOK|wxICON_EXCLAMATION);
         dumbUser.ShowModal();
         return;
      }
      mpMolBond->SetAtom2(*at);
   }
   
   this->CrystUpdate();
   this->UpdateUI();
   VFN_DEBUG_EXIT("WXMolBond::OnChangeScattPow()",6)
}

////////////////////////////////////////////////////////////////////////
//
//    WXMolBondAngle
//
////////////////////////////////////////////////////////////////////////
WXCRYST_ID ID_MOLBONDANGLE_ATOM1;
WXCRYST_ID ID_MOLBONDANGLE_ATOM2;
WXCRYST_ID ID_MOLBONDANGLE_ATOM3;

BEGIN_EVENT_TABLE(WXMolBondAngle,wxWindow)
   EVT_BUTTON(ID_MOLBONDANGLE_ATOM1,    WXMolBondAngle::OnChangeAtom)
   EVT_BUTTON(ID_MOLBONDANGLE_ATOM2,    WXMolBondAngle::OnChangeAtom)
   EVT_BUTTON(ID_MOLBONDANGLE_ATOM3,    WXMolBondAngle::OnChangeAtom)
END_EVENT_TABLE()

WXMolBondAngle::WXMolBondAngle(wxWindow *parent, MolBondAngle*obj):
WXCrystObjBasic(parent),mpMolBondAngle(obj)
{
   VFN_DEBUG_ENTRY("WXMolBondAngle::WXMolBond()",6)
   mpSizer=new wxBoxSizer(wxHORIZONTAL);

   mpFieldAtom1=new WXFieldChoice(this,ID_MOLBONDANGLE_ATOM1,"Bond Angle:",60);
   mpFieldAtom1->SetValue(mpMolBondAngle->GetAtom1().GetName());
   mpSizer->Add(mpFieldAtom1,0,wxALIGN_CENTER);
   mList.Add(mpFieldAtom1);
   
   mpFieldAtom2=new WXFieldChoice(this,ID_MOLBONDANGLE_ATOM2,"-",60);
   mpFieldAtom2->SetValue(mpMolBondAngle->GetAtom2().GetName());
   mpSizer->Add(mpFieldAtom2,0,wxALIGN_CENTER);
   mList.Add(mpFieldAtom2);
   
   mpFieldAtom3=new WXFieldChoice(this,ID_MOLBONDANGLE_ATOM3,"-",60);
   mpFieldAtom3->SetValue(mpMolBondAngle->GetAtom3().GetName());
   mpSizer->Add(mpFieldAtom3,0,wxALIGN_CENTER);
   mList.Add(mpFieldAtom3);

   WXFieldPar<REAL> *angle=
      new WXFieldPar<REAL>(this,"Angle:",WXCRYST_ID(),&(mpMolBondAngle->Angle0()));
   mpSizer->Add(angle,0,wxALIGN_CENTER);
   mList.Add(angle);

   this->SetSizer(mpSizer);
   this->CrystUpdate();
   this->Layout();
   VFN_DEBUG_EXIT("WXMolBondAngle::WXMolBond()",6)
}

void WXMolBondAngle::CrystUpdate()
{
   VFN_DEBUG_ENTRY("WXMolBondAngle::CrystUpdate()",5)
   mList.CrystUpdate();
   VFN_DEBUG_EXIT("WXMolBondAngle::CrystUpdate()",5)
}

void WXMolBondAngle::UpdateUI()
{
   VFN_DEBUG_ENTRY("WXMolBondAngle::UpdateUI()",5)
   mList.UpdateUI();
   mpFieldAtom1->SetValue(mpMolBondAngle->GetAtom1().GetName().c_str());
   mpFieldAtom2->SetValue(mpMolBondAngle->GetAtom2().GetName().c_str());
   mpFieldAtom3->SetValue(mpMolBondAngle->GetAtom3().GetName().c_str());
   VFN_DEBUG_EXIT("WXMolBondAngle::UpdateUI()",5)
}

bool WXMolBondAngle::Layout()
{
   VFN_DEBUG_ENTRY("WXMolBond::Layout()",5)
   //:TODO: Cleeeeaaaaannnn thiiiisss !!!
   for(unsigned int i=0;i<mList.GetNb();i++)
      mpSizer->SetItemMinSize(mList.Get(i),
                              mList.Get(i)->GetSize().GetWidth(),
                              mList.Get(i)->GetSize().GetHeight());
      
   mpSizer->Layout();
   mpSizer->Fit(this);
   wxSizer* s=mWXParent->GetSizer();
   if(s != 0)
   {// Need to do it that way, in case  the parent is not a WXCrystObj
    // with an adequate Layout() function
      s->SetItemMinSize(this,this->GetSize().GetWidth(),this->GetSize().GetHeight());
      s->Fit(mWXParent);
   }
   mWXParent->Layout();
   VFN_DEBUG_EXIT("WXMolBondAngle::Layout()",5)
   return this->wxWindow::Layout();
}

void WXMolBondAngle::OnChangeAtom(wxCommandEvent &event)
{
   VFN_DEBUG_ENTRY("WXMolBondAngle::OnChangeAtom()",6)
   WXCrystValidateAllUserInput();
   int choice;
   const MolAtom *at=WXDialogChooseFromVector(mpMolBondAngle->GetMolecule().GetAtomList(),
                               (wxWindow*)this,"Choose a new Atom",choice);
   if(0==at) return;
   if((long)ID_MOLBONDANGLE_ATOM1==event.GetId())
   {
      if((at==&(mpMolBondAngle->GetAtom2()) )||(at==&(mpMolBondAngle->GetAtom3())) )
      {
         wxMessageDialog dumbUser(this,"The three atoms must be different !",
                                  "Whooops",wxOK|wxICON_EXCLAMATION);
         dumbUser.ShowModal();
         return;
      }
      mpMolBondAngle->SetAtom1(*at);
   }
   if((long)ID_MOLBONDANGLE_ATOM2==event.GetId())
   {
      if((at==&(mpMolBondAngle->GetAtom1()) )||(at==&(mpMolBondAngle->GetAtom3())) )
      {
         wxMessageDialog dumbUser(this,"The three atoms must be different !",
                                  "Whooops",wxOK|wxICON_EXCLAMATION);
         dumbUser.ShowModal();
         return;
      }
      mpMolBondAngle->SetAtom2(*at);
   }
   if((long)ID_MOLBONDANGLE_ATOM3==event.GetId())
   {
      if((at==&(mpMolBondAngle->GetAtom1()) )||(at==&(mpMolBondAngle->GetAtom2())) )
      {
         wxMessageDialog dumbUser(this,"The three atoms must be different !",
                                  "Whooops",wxOK|wxICON_EXCLAMATION);
         dumbUser.ShowModal();
         return;
      }
      mpMolBondAngle->SetAtom3(*at);
   }
   
   this->CrystUpdate();
   this->UpdateUI();
   VFN_DEBUG_EXIT("WXMolBondAngle::OnChangeScattPow()",6)
}
////////////////////////////////////////////////////////////////////////
//
//    WXMolDihedralAngle
//
////////////////////////////////////////////////////////////////////////
WXCRYST_ID ID_MOLDIHEDRALANGLE_ATOM1;
WXCRYST_ID ID_MOLDIHEDRALANGLE_ATOM2;
WXCRYST_ID ID_MOLDIHEDRALANGLE_ATOM3;
WXCRYST_ID ID_MOLDIHEDRALANGLE_ATOM4;

BEGIN_EVENT_TABLE(WXMolDihedralAngle,wxWindow)
   EVT_BUTTON(ID_MOLDIHEDRALANGLE_ATOM1,    WXMolDihedralAngle::OnChangeAtom)
   EVT_BUTTON(ID_MOLDIHEDRALANGLE_ATOM2,    WXMolDihedralAngle::OnChangeAtom)
   EVT_BUTTON(ID_MOLDIHEDRALANGLE_ATOM3,    WXMolDihedralAngle::OnChangeAtom)
   EVT_BUTTON(ID_MOLDIHEDRALANGLE_ATOM4,    WXMolDihedralAngle::OnChangeAtom)
END_EVENT_TABLE()

WXMolDihedralAngle::WXMolDihedralAngle(wxWindow *parent, MolDihedralAngle*obj):
WXCrystObjBasic(parent),mpMolDihedralAngle(obj)
{
   VFN_DEBUG_ENTRY("WXMolDihedralAngle::WXMolBond()",6)
   mpSizer=new wxBoxSizer(wxHORIZONTAL);

   mpFieldAtom1=new WXFieldChoice(this,ID_MOLDIHEDRALANGLE_ATOM1,"Dihedral Angle:",60);
   mpFieldAtom1->SetValue(mpMolDihedralAngle->GetAtom1().GetName());
   mpSizer->Add(mpFieldAtom1,0,wxALIGN_CENTER);
   mList.Add(mpFieldAtom1);
   
   mpFieldAtom2=new WXFieldChoice(this,ID_MOLDIHEDRALANGLE_ATOM2,"-",60);
   mpFieldAtom2->SetValue(mpMolDihedralAngle->GetAtom2().GetName());
   mpSizer->Add(mpFieldAtom2,0,wxALIGN_CENTER);
   mList.Add(mpFieldAtom2);
   
   mpFieldAtom3=new WXFieldChoice(this,ID_MOLDIHEDRALANGLE_ATOM3,"-",60);
   mpFieldAtom3->SetValue(mpMolDihedralAngle->GetAtom3().GetName());
   mpSizer->Add(mpFieldAtom3,0,wxALIGN_CENTER);
   mList.Add(mpFieldAtom3);

   mpFieldAtom4=new WXFieldChoice(this,ID_MOLDIHEDRALANGLE_ATOM4,"-",60);
   mpFieldAtom4->SetValue(mpMolDihedralAngle->GetAtom4().GetName());
   mpSizer->Add(mpFieldAtom4,0,wxALIGN_CENTER);
   mList.Add(mpFieldAtom4);

   WXFieldPar<REAL> *angle=
      new WXFieldPar<REAL>(this,"Angle:",-1,&(mpMolDihedralAngle->Angle0()));
   mpSizer->Add(angle,0,wxALIGN_CENTER);
   mList.Add(angle);

   this->SetSizer(mpSizer);
   this->CrystUpdate();
   this->Layout();
   VFN_DEBUG_EXIT("WXMolDihedralAngle::WXMolBond()",6)
}

void WXMolDihedralAngle::CrystUpdate()
{
   VFN_DEBUG_ENTRY("WXMolDihedralAngle::CrystUpdate()",5)
   mList.CrystUpdate();
   VFN_DEBUG_EXIT("WXMolDihedralAngle::CrystUpdate()",5)
}

void WXMolDihedralAngle::UpdateUI()
{
   VFN_DEBUG_ENTRY("WXMolDihedralAngle::UpdateUI()",5)
   mList.UpdateUI();
   mpFieldAtom1->SetValue(mpMolDihedralAngle->GetAtom1().GetName().c_str());
   mpFieldAtom2->SetValue(mpMolDihedralAngle->GetAtom2().GetName().c_str());
   mpFieldAtom3->SetValue(mpMolDihedralAngle->GetAtom3().GetName().c_str());
   VFN_DEBUG_EXIT("WXMolDihedralAngle::UpdateUI()",5)
}

bool WXMolDihedralAngle::Layout()
{
   VFN_DEBUG_ENTRY("WXMolBond::Layout()",5)
   //:TODO: Cleeeeaaaaannnn thiiiisss !!!
   for(unsigned int i=0;i<mList.GetNb();i++)
      mpSizer->SetItemMinSize(mList.Get(i),
                              mList.Get(i)->GetSize().GetWidth(),
                              mList.Get(i)->GetSize().GetHeight());
      
   mpSizer->Layout();
   mpSizer->Fit(this);
   wxSizer* s=mWXParent->GetSizer();
   if(s != 0)
   {// Need to do it that way, in case  the parent is not a WXCrystObj
    // with an adequate Layout() function
      s->SetItemMinSize(this,this->GetSize().GetWidth(),this->GetSize().GetHeight());
      s->Fit(mWXParent);
   }
   mWXParent->Layout();
   VFN_DEBUG_EXIT("WXMolDihedralAngle::Layout()",5)
   return this->wxWindow::Layout();
}

void WXMolDihedralAngle::OnChangeAtom(wxCommandEvent &event)
{
   VFN_DEBUG_ENTRY("WXMolDihedralAngle::OnChangeAtom()",6)
   WXCrystValidateAllUserInput();
   int choice;
   const MolAtom *at=WXDialogChooseFromVector(mpMolDihedralAngle->GetMolecule().GetAtomList(),
                               (wxWindow*)this,"Choose a new Atom",choice);
   if(0==at) return;
   if((long)ID_MOLDIHEDRALANGLE_ATOM1==event.GetId())
   {
      if(  (at==&(mpMolDihedralAngle->GetAtom2()))
         ||(at==&(mpMolDihedralAngle->GetAtom3()))
         ||(at==&(mpMolDihedralAngle->GetAtom4())) )
      {
         wxMessageDialog dumbUser(this,"The four atoms must be different !",
                                  "Whooops",wxOK|wxICON_EXCLAMATION);
         dumbUser.ShowModal();
         return;
      }
      mpMolDihedralAngle->SetAtom1(*at);
   }
   if((long)ID_MOLDIHEDRALANGLE_ATOM2==event.GetId())
   {
      if(  (at==&(mpMolDihedralAngle->GetAtom1()))
         ||(at==&(mpMolDihedralAngle->GetAtom3()))
         ||(at==&(mpMolDihedralAngle->GetAtom4())) )
      {
         wxMessageDialog dumbUser(this,"The four atoms must be different !",
                                  "Whooops",wxOK|wxICON_EXCLAMATION);
         dumbUser.ShowModal();
         return;
      }
      mpMolDihedralAngle->SetAtom2(*at);
   }
   if((long)ID_MOLDIHEDRALANGLE_ATOM3==event.GetId())
   {
      if(  (at==&(mpMolDihedralAngle->GetAtom1()))
         ||(at==&(mpMolDihedralAngle->GetAtom2()))
         ||(at==&(mpMolDihedralAngle->GetAtom4())) )
      {
         wxMessageDialog dumbUser(this,"The four atoms must be different !",
                                  "Whooops",wxOK|wxICON_EXCLAMATION);
         dumbUser.ShowModal();
         return;
      }
      mpMolDihedralAngle->SetAtom3(*at);
   }
   if((long)ID_MOLDIHEDRALANGLE_ATOM4==event.GetId())
   {
      if(  (at==&(mpMolDihedralAngle->GetAtom1()))
         ||(at==&(mpMolDihedralAngle->GetAtom2()))
         ||(at==&(mpMolDihedralAngle->GetAtom3())) )
      {
         wxMessageDialog dumbUser(this,"The four atoms must be different !",
                                  "Whooops",wxOK|wxICON_EXCLAMATION);
         dumbUser.ShowModal();
         return;
      }
      mpMolDihedralAngle->SetAtom4(*at);
   }
   
   this->CrystUpdate();
   this->UpdateUI();
   VFN_DEBUG_EXIT("WXMolDihedralAngle::OnChangeScattPow()",6)
}
////////////////////////////////////////////////////////////////////////
//
//    WXMolecule
//
////////////////////////////////////////////////////////////////////////

WXCRYST_ID ID_MENU_OPTIMIZECONFORMATION;
WXCRYST_ID ID_MENU_SETLIMITS;
WXCRYST_ID ID_MOLECULE_MENU_FILE;
WXCRYST_ID ID_MOLECULE_MENU_FORMULA;
WXCRYST_ID ID_MOLECULE_MENU_FORMULA_OPTIMIZECONFORMATION;
WXCRYST_ID ID_MOLECULE_MENU_FORMULA_ADD_ATOM;
WXCRYST_ID ID_MOLECULE_MENU_FORMULA_ADD_BOND;
WXCRYST_ID ID_MOLECULE_MENU_FORMULA_ADD_ANGLE;
WXCRYST_ID ID_MOLECULE_MENU_FORMULA_ADD_DIHEDRAL;
WXCRYST_ID ID_MOLECULE_MENU_FORMULA_TEST;
WXCRYST_ID ID_MOLECULE_MENU_FORMULA_REMOVE_ATOM;
WXCRYST_ID ID_MOLECULE_MENU_FORMULA_REMOVE_BOND;
WXCRYST_ID ID_MOLECULE_MENU_FORMULA_REMOVE_ANGLE;
WXCRYST_ID ID_MOLECULE_MENU_FORMULA_REMOVE_DIHEDRAL;

BEGIN_EVENT_TABLE(WXMolecule,wxWindow)
   EVT_BUTTON(ID_WXOBJ_COLLAPSE,                          WXCrystObj::OnToggleCollapse)
   EVT_MENU(ID_REFOBJ_MENU_PAR_FIXALL,                    WXRefinableObj::OnMenuFixAllPar)
   EVT_MENU(ID_REFOBJ_MENU_PAR_UNFIXALL,                  WXRefinableObj::OnMenuUnFixAllPar)
   EVT_MENU(ID_REFOBJ_MENU_PAR_RANDOMIZE,                 WXRefinableObj::OnMenuParRandomize)
   EVT_MENU(ID_MOLECULE_MENU_FORMULA_OPTIMIZECONFORMATION,WXMolecule::OnMenuOptimizeConformation)
   EVT_MENU(ID_MOLECULE_MENU_FORMULA_ADD_ATOM,            WXMolecule::OnMenuAddAtom)
   EVT_MENU(ID_MOLECULE_MENU_FORMULA_ADD_BOND,            WXMolecule::OnMenuAddBond)
   EVT_MENU(ID_MOLECULE_MENU_FORMULA_ADD_ANGLE,           WXMolecule::OnMenuAddAngle)
   EVT_MENU(ID_MOLECULE_MENU_FORMULA_ADD_DIHEDRAL,        WXMolecule::OnMenuAddDihedralAngle)
   EVT_MENU(ID_MOLECULE_MENU_FORMULA_REMOVE_ATOM,         WXMolecule::OnMenuRemoveAtom)
   EVT_MENU(ID_MOLECULE_MENU_FORMULA_REMOVE_BOND,         WXMolecule::OnMenuRemoveBond)
   EVT_MENU(ID_MOLECULE_MENU_FORMULA_REMOVE_ANGLE,        WXMolecule::OnMenuRemoveAngle)
   EVT_MENU(ID_MOLECULE_MENU_FORMULA_REMOVE_DIHEDRAL,     WXMolecule::OnMenuRemoveDihedralAngle)
   EVT_MENU(ID_MOLECULE_MENU_FORMULA_TEST        ,        WXMolecule::OnMenuTest)
   EVT_MENU(ID_MENU_SETLIMITS,                            WXMolecule::OnMenuSetLimits)
END_EVENT_TABLE()

WXMolecule::WXMolecule(wxWindow *parent, Molecule *mol):
WXScatterer(parent,mol),mpMolecule(mol)
{
   VFN_DEBUG_ENTRY("WXMolecule::WXMolecule()",6)
   //Menus
      mpMenuBar->AddMenu("File",ID_MOLECULE_MENU_FILE);
      mpMenuBar->AddMenu("Parameters",ID_REFOBJ_MENU_PAR);
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_PAR,ID_REFOBJ_MENU_PAR_FIXALL,"Fix all");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_PAR,ID_REFOBJ_MENU_PAR_UNFIXALL,"Unfix all");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_PAR,ID_REFOBJ_MENU_PAR_RANDOMIZE,
                                "Randomize Configuration");
      mpMenuBar->AddMenu("Formula",ID_MOLECULE_MENU_FORMULA);
         mpMenuBar->AddMenuItem(ID_MOLECULE_MENU_FORMULA,ID_MOLECULE_MENU_FORMULA_OPTIMIZECONFORMATION,
                                "Optimize Starting Conformation");
         mpMenuBar->AddMenuItem(ID_MOLECULE_MENU_FORMULA,ID_MOLECULE_MENU_FORMULA_ADD_ATOM,
                                "Add an Atom");
         mpMenuBar->AddMenuItem(ID_MOLECULE_MENU_FORMULA,ID_MOLECULE_MENU_FORMULA_ADD_BOND,
                                "Add Bond Restraint");
         mpMenuBar->AddMenuItem(ID_MOLECULE_MENU_FORMULA,ID_MOLECULE_MENU_FORMULA_ADD_ANGLE,
                                "Add Bond Angle Restraint");
         mpMenuBar->AddMenuItem(ID_MOLECULE_MENU_FORMULA,ID_MOLECULE_MENU_FORMULA_ADD_DIHEDRAL,
                                "Add Dihedral Angle Restraint");
         mpMenuBar->GetMenu(ID_MOLECULE_MENU_FORMULA).AppendSeparator();
         mpMenuBar->AddMenuItem(ID_MOLECULE_MENU_FORMULA,ID_MOLECULE_MENU_FORMULA_REMOVE_ATOM,
                                "Remove an Atom");
         mpMenuBar->AddMenuItem(ID_MOLECULE_MENU_FORMULA,ID_MOLECULE_MENU_FORMULA_REMOVE_BOND,
                                "Remove a Bond Restraint");
         mpMenuBar->AddMenuItem(ID_MOLECULE_MENU_FORMULA,ID_MOLECULE_MENU_FORMULA_REMOVE_ANGLE,
                                "Remove a Bond Angle Restraint");
         mpMenuBar->AddMenuItem(ID_MOLECULE_MENU_FORMULA,ID_MOLECULE_MENU_FORMULA_REMOVE_DIHEDRAL,
                                "Remove a Dihedral Angle Restraint");
         mpMenuBar->GetMenu(ID_MOLECULE_MENU_FORMULA).AppendSeparator();
         mpMenuBar->AddMenuItem(ID_MOLECULE_MENU_FORMULA,ID_MOLECULE_MENU_FORMULA_TEST,
                                "Test");
   
   //sizers
   mpSizerAtomList= new wxBoxSizer(wxVERTICAL);
   mpSizerBondList= new wxBoxSizer(wxVERTICAL);
   mpSizerAngleList= new wxBoxSizer(wxVERTICAL);
   mpSizerDihedralAngleList= new wxBoxSizer(wxVERTICAL);

   mpSizer->Add(mpSizerAtomList,0,wxALIGN_LEFT);
   mpSizer->Add(mpSizerBondList,0,wxALIGN_LEFT);
   mpSizer->Add(mpSizerAngleList,0,wxALIGN_LEFT);
   mpSizer->Add(mpSizerDihedralAngleList,0,wxALIGN_LEFT);
   this->CrystUpdate();
   this->Layout();
   VFN_DEBUG_EXIT("WXMolecule::WXMolecule()",6)
}

WXMolecule::~WXMolecule()
{
   VFN_DEBUG_ENTRY("WXMolecule::~WXMolecule()",10)
   VFN_DEBUG_EXIT("WXMolecule::~WXMolecule()",10)
}

void WXMolecule::OnMenuOptimizeConformation(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_ENTRY("WXMolecule::OnMenuOptimizeConformation()",5)
   WXCrystValidateAllUserInput();
   mpMolecule->OptimizeConformation();
   mpMolecule->GetCrystal().UpdateDisplay();
   VFN_DEBUG_EXIT("WXMolecule::OnMenuOptimizeConformation()",5)
}

void WXMolecule::OnMenuAddAtom(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_ENTRY("WXMolecule::OnMenuAddAtom()",6)
   WXCrystValidateAllUserInput();
   int choice;
   const ScatteringPower *scatt=WXDialogChooseFromRegistry(
               mpMolecule->GetCrystal().GetScatteringPowerRegistry(),
               (wxWindow*)this,"Choose a new Scattering Power",choice);
   if(0==scatt) return;
   stringstream st;
   st<<"_"<<mpMolecule->GetAtomList().size();
   mpMolecule->AddAtom(0.,0.,0.,scatt,scatt->GetName()+st.str());
   VFN_DEBUG_EXIT("WXMolecule::OnMenuAddAtom()",6)
}

void WXMolecule::OnMenuAddBond(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_ENTRY("WXMolecule::OnMenuAddBond()",6)
   WXCrystValidateAllUserInput();
   int choice;
   vector<MolAtom*> v=mpMolecule->GetAtomList();
   MolAtom *at1=WXDialogChooseFromVector(v,
                               (wxWindow*)this,"Choose the first Atom",choice);
   if(0==at1) return;
   MolAtom *at2=WXDialogChooseFromVector(v,
                               (wxWindow*)this,"Choose the second Atom",choice);
   if(0==at2) return;
   
   if(at1==at2)
   {
      wxMessageDialog dumbUser(this,"The two atoms must be different !",
                               "Whooops",wxOK|wxICON_EXCLAMATION);
      dumbUser.ShowModal();
      return;
   }
   mpMolecule->AddBond(*at1,*at2,1.5,.01,.05,1.);
   VFN_DEBUG_EXIT("WXMolecule::OnMenuAddBond()",6)
}

void WXMolecule::OnMenuAddAngle(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_ENTRY("WXMolecule::OnMenuAddBond()",6)
   WXCrystValidateAllUserInput();
   int choice;
   vector<MolAtom*> v=mpMolecule->GetAtomList();
   MolAtom *at1=WXDialogChooseFromVector(v,
                               (wxWindow*)this,"Choose the first Atom",choice);
   if(0==at1) return;
   MolAtom *at2=WXDialogChooseFromVector(v,
                               (wxWindow*)this,"Choose the second Atom",choice);
   if(0==at2) return;

   MolAtom *at3=WXDialogChooseFromVector(v,
                               (wxWindow*)this,"Choose the third Atom",choice);
   if(0==at2) return;
   
   if( (at1==at2) || (at1==at3) ||(at2==at3))
   {
      wxMessageDialog dumbUser(this,"The three atoms must be different !",
                               "Whooops",wxOK|wxICON_EXCLAMATION);
      dumbUser.ShowModal();
      return;
   }
   mpMolecule->AddBondAngle(*at1,*at2,*at3,109*DEG2RAD,.01,0.05);
   VFN_DEBUG_EXIT("WXMolecule::OnMenuAddBond()",6)
}

void WXMolecule::OnMenuAddDihedralAngle(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_ENTRY("WXMolecule::OnMenuAddDihedralAngle()",6)
   WXCrystValidateAllUserInput();
   int choice;
   vector<MolAtom*> v=mpMolecule->GetAtomList();
   MolAtom *at1=WXDialogChooseFromVector(v,
                               (wxWindow*)this,"Choose the first Atom",choice);
   if(0==at1) return;
   MolAtom *at2=WXDialogChooseFromVector(v,
                               (wxWindow*)this,"Choose the second Atom",choice);
   if(0==at2) return;

   MolAtom *at3=WXDialogChooseFromVector(v,
                               (wxWindow*)this,"Choose the third Atom",choice);
   if(0==at3) return;

   MolAtom *at4=WXDialogChooseFromVector(v,
                               (wxWindow*)this,"Choose the fourth Atom",choice);
   if(0==at4) return;
   
   if( (at1==at2) || (at1==at3) || (at1==at4) || (at2==at3) || (at2==at4) || (at3==at4))
   {
      wxMessageDialog dumbUser(this,"The atoms must be different !",
                               "Whooops",wxOK|wxICON_EXCLAMATION);
      dumbUser.ShowModal();
      return;
   }
   mpMolecule->AddDihedralAngle(*at1,*at2,*at3,*at4,180*DEG2RAD,.01,.05);
   VFN_DEBUG_EXIT("WXMolecule::OnMenuAddDihedralAngle()",6)
}

void WXMolecule::OnMenuRemoveAtom(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_ENTRY("WXMolecule::OnMenuRemoveAtom()",6)
   vector<MolAtom*> v=mpMolecule->GetAtomList();
   int choice;
   MolAtom *at=WXDialogChooseFromVector(v,
                               (wxWindow*)this,"Choose the Atom to be removed",choice);
   if(0==at) return;
   mpMolecule->RemoveAtom(*at);
   this->CrystUpdate();
   VFN_DEBUG_EXIT("WXMolecule::OnMenuRemoveAtom()",6)
}

void WXMolecule::OnMenuRemoveBond(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_ENTRY("WXMolecule::OnMenuRemoveBond()",6)
   vector<MolBond*> v=mpMolecule->GetBondList();
   int choice;
   MolBond *b=WXDialogChooseFromVector(v,
                               (wxWindow*)this,"Choose the Bond to be removed",choice);
   if(0==b) return;
   mpMolecule->RemoveBond(*b);
   this->CrystUpdate();
   VFN_DEBUG_EXIT("WXMolecule::OnMenuRemoveBond()",6)
}

void WXMolecule::OnMenuRemoveAngle(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_ENTRY("WXMolecule::OnMenuRemoveAngle()",6)
   vector<MolBondAngle*> v=mpMolecule->GetBondAngleList();
   int choice;
   MolBondAngle *a=WXDialogChooseFromVector(v,
                               (wxWindow*)this,"Choose the Bond Angle to be removed",choice);
   if(0==a) return;
   mpMolecule->RemoveBondAngle(*a);
   this->CrystUpdate();
   VFN_DEBUG_EXIT("WXMolecule::OnMenuRemoveAngle()",6)
}

void WXMolecule::OnMenuRemoveDihedralAngle(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_ENTRY("WXMolecule::OnMenuRemoveDihedralAngle()",6)
   vector<MolDihedralAngle*> v=mpMolecule->GetDihedralAngleList();
   int choice;
   MolDihedralAngle *a=WXDialogChooseFromVector(v,
                               (wxWindow*)this,"Choose the Dihedral Angle to be removed",choice);
   if(0==a) return;
   mpMolecule->RemoveDihedralAngle(*a);
   this->CrystUpdate();
   VFN_DEBUG_EXIT("WXMolecule::OnMenuRemoveDihedralAngle()",6)
}

void WXMolecule::OnMenuTest(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_ENTRY("WXMolecule::OnMenuTest()",6)
   mpMolecule->BeginOptimization();
   mpMolecule->EndOptimization();
   VFN_DEBUG_EXIT("WXMolecule::OnMenuTest()",6)
}

void WXMolecule::OnMenuSetLimits(wxCommandEvent &event)
{
}

void WXMolecule::CrystUpdate()
{
   VFN_DEBUG_ENTRY("WXMolecule::CrystUpdate()",6)
   if(false==mpMolecule->IsBeingRefined())
   {
      //Remove any atom, bond, bond angle or dihedral angle that could have been removed
      {
         vector<MolAtom*>::iterator pos;
         for(pos=mvpAtom.begin();pos!=mvpAtom.end();pos++)
         {
            vector<MolAtom*>::const_iterator pos2=find(mpMolecule->GetAtomList().begin(),
                                                       mpMolecule->GetAtomList().end(),*pos);
            if(pos2==mpMolecule->GetAtomList().end())
            {
               mList.Remove((*pos)->WXGet());
               mpSizerAtomList->Remove((*pos)->WXGet());
               (*pos)->WXDelete();
               pos=mvpAtom.erase(pos);
               --pos;
            }
         }
      }
      {
         vector<MolBond*>::iterator pos;
         for(pos=mvpBond.begin();pos!=mvpBond.end();pos++)
         {
            vector<MolBond*>::const_iterator pos2=find(mpMolecule->GetBondList().begin(),
                                                       mpMolecule->GetBondList().end(),*pos);
            if(pos2==mpMolecule->GetBondList().end())
            {
               mList.Remove((*pos)->WXGet());
               mpSizerBondList->Remove((*pos)->WXGet());
               (*pos)->WXDelete();
               pos=mvpBond.erase(pos);
               --pos;
            }
         }
      }
      {
         vector<MolBondAngle*>::iterator pos;
         for(pos=mvpBondAngle.begin();pos!=mvpBondAngle.end();pos++)
         {
            vector<MolBondAngle*>::const_iterator pos2=
                                    find(mpMolecule->GetBondAngleList().begin(),
                                         mpMolecule->GetBondAngleList().end(),*pos);
            if(pos2==mpMolecule->GetBondAngleList().end())
            {
               mList.Remove((*pos)->WXGet());
               mpSizerAngleList->Remove((*pos)->WXGet());
               (*pos)->WXDelete();
               pos=mvpBondAngle.erase(pos);
               --pos;
            }
         }
      }
      {
         vector<MolDihedralAngle*>::iterator pos;
         for(pos=mvpDihedralAngle.begin();pos!=mvpDihedralAngle.end();pos++)
         {
            vector<MolDihedralAngle*>::const_iterator pos2=
                              find(mpMolecule->GetDihedralAngleList().begin(),
                                   mpMolecule->GetDihedralAngleList().end(),*pos);
            if(pos2==mpMolecule->GetDihedralAngleList().end())
            {
               mList.Remove((*pos)->WXGet());
               mpSizerDihedralAngleList->Remove((*pos)->WXGet());
               (*pos)->WXDelete();
               pos=mvpDihedralAngle.erase(pos);
               --pos;
            }
         }
      }
      //Add any atom, bond, bond angle or dihedral angle that could have been added
      {
         vector<MolAtom*>::iterator pos;
         for(pos=mpMolecule->GetAtomList().begin();pos!=mpMolecule->GetAtomList().end();pos++)
         {
            vector<MolAtom*>::const_iterator pos2=find(mvpAtom.begin(),mvpAtom.end(),*pos);
            if(pos2==mvpAtom.end())
            {
               VFN_DEBUG_MESSAGE("WXMolecule::CrystUpdate():Atom not found:"<<(*pos)->GetName(),5)
               WXCrystObjBasic *at=(*pos)->WXCreate(this);
               mpSizerAtomList->Add(at);
               mList.Add(at);
               mpSizerAtomList->Layout();
               mvpAtom.push_back(*pos);
            }
         }
      }
      {
         vector<MolBond*>::iterator pos;
         for(pos=mpMolecule->GetBondList().begin();pos!=mpMolecule->GetBondList().end();pos++)
         {
            vector<MolBond*>::const_iterator pos2=find(mvpBond.begin(),mvpBond.end(),*pos);
            if(pos2==mvpBond.end())
            {
               VFN_DEBUG_MESSAGE("WXMolecule::CrystUpdate():Bond not found",5)
               WXCrystObjBasic *b=(*pos)->WXCreate(this);
               mpSizerBondList->Add(b);
               mList.Add(b);
               mpSizerBondList->Layout();
               mvpBond.push_back(*pos);
            }
         }
      }
      {
         vector<MolBondAngle*>::iterator pos;
         for(pos=mpMolecule->GetBondAngleList().begin();
             pos!=mpMolecule->GetBondAngleList().end();pos++)
         {
            vector<MolBondAngle*>::const_iterator pos2
               =find(mvpBondAngle.begin(),mvpBondAngle.end(),*pos);
            if(pos2==mvpBondAngle.end())
            {
               VFN_DEBUG_MESSAGE("WXMolecule::CrystUpdate():Bond not found",5)
               WXCrystObjBasic *b=(*pos)->WXCreate(this);
               mpSizerAngleList->Add(b);
               mList.Add(b);
               mpSizerAngleList->Layout();
               mvpBondAngle.push_back(*pos);
            }
         }
      }
      {
         vector<MolDihedralAngle*>::iterator pos;
         for(pos=mpMolecule->GetDihedralAngleList().begin();
             pos!=mpMolecule->GetDihedralAngleList().end();pos++)
         {
            vector<MolDihedralAngle*>::const_iterator pos2
               =find(mvpDihedralAngle.begin(),mvpDihedralAngle.end(),*pos);
            if(pos2==mvpDihedralAngle.end())
            {
               VFN_DEBUG_MESSAGE("WXMolecule::CrystUpdate():Bond not found",5)
               WXCrystObjBasic *b=(*pos)->WXCreate(this);
               mpSizerDihedralAngleList->Add(b);
               mList.Add(b);
               mpSizerDihedralAngleList->Layout();
               mvpDihedralAngle.push_back(*pos);
            }
         }
      }
      this->Layout();
   }
   this->WXRefinableObj::CrystUpdate();
   VFN_DEBUG_EXIT("WXMolecule::CrystUpdate()",6)
}
void WXMolecule::UpdateUI()
{
   VFN_DEBUG_ENTRY("WXMolecule::UpdateUI()",5)
   this->WXRefinableObj::UpdateUI();
   VFN_DEBUG_EXIT("WXMolecule::UpdateUI()",5)
}
} //namespace
