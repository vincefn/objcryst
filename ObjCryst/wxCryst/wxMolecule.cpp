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

template<class T> list<T *> WXDialogChooseMultipleFromVector(vector<T*> &reg,wxWindow*parent,
                                                const string &message)
{
   wxString *choices = new wxString[reg.size()];
   for(unsigned int i=0;i<reg.size();i++) 
      choices[i]=(reg[i]->GetName()).c_str();
   wxMultiChoiceDialog dialog
         (parent,message.c_str(),"Choose",reg.size(),choices,wxOK | wxCANCEL);
   dialog.SetSize(300,300);
   dialog.ShowModal();
   wxArrayInt choice=dialog.GetSelections();
   list<T*> vChoice;
   for(unsigned int i=0;i<choice.GetCount();++i) vChoice.push_back(reg[choice.Item(i)]);
   delete[] choices;
   return vChoice;
}
template<class T> list<T const*> WXDialogChooseMultipleFromVector(const vector<T*> &reg,wxWindow*parent,
                                                const string &message)
{
   wxString *choices = new wxString[reg.size()];
   for(unsigned int i=0;i<reg.size();i++) 
      choices[i]=(reg[i]->GetName()).c_str();
   wxMultiChoiceDialog dialog
         (parent,message.c_str(),"Choose",reg.size(),choices,wxOK | wxCANCEL);
   dialog.SetSize(300,300);
   dialog.ShowModal();
   wxArrayInt choice=dialog.GetSelections();
   list<T const*> vChoice;
   for(unsigned int i=0;i<choice.GetCount();++i) vChoice.push_back(reg[choice.Item(i)]);
   delete[] choices;
   return vChoice;
}
////////////////////////////////////////////////////////////////////////
//
//    WXMolScrolledWindow
//
////////////////////////////////////////////////////////////////////////
WXMolScrolledWindow::WXMolScrolledWindow(wxWindow* parent, WXMolecule* pWXMol):
wxScrolledWindow(parent),mpWXMolecule(pWXMol)
{}
WXMolScrolledWindow::~WXMolScrolledWindow()
{
   mpWXMolecule->NotifyDeleteListWin(this);
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
   VFN_DEBUG_ENTRY("WXMolAtom::WXMolAtom():"<<obj->GetName(),6)
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
                               &(mpMolAtom->GetMolecule().GetPar(&(mpMolAtom->X()))),
                               50,false,false);
      mpSizer->Add(pFieldX,0,wxALIGN_CENTER);
      mList.Add(pFieldX);
   }
   {
      WXFieldRefPar* pFieldY  =new WXFieldRefPar(this,"y",
                               &(mpMolAtom->GetMolecule().GetPar(&(mpMolAtom->Y()))),
                               50,false,false);
      mpSizer->Add(pFieldY,0,wxALIGN_CENTER);
      mList.Add(pFieldY);
   }
   {
      WXFieldRefPar* pFieldZ  =new WXFieldRefPar(this,"z",
                               &(mpMolAtom->GetMolecule().GetPar(&(mpMolAtom->Z()))),
                               50,false,false);
      mpSizer->Add(pFieldZ,0,wxALIGN_CENTER);
      mList.Add(pFieldZ);
   }
   
   this->SetSizer(mpSizer);
   this->BottomLayout(0);
   this->CrystUpdate();
   VFN_DEBUG_EXIT("WXMolAtom::WXMolAtom():"<<obj->GetName(),6)
}

WXMolAtom::~WXMolAtom()
{
   mpMolAtom->WXNotifyDelete();
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
WXCRYST_ID ID_MOLBOND_FREEBUTTON;

BEGIN_EVENT_TABLE(WXMolBond,wxWindow)
   EVT_BUTTON(ID_MOLBOND_ATOM1,    WXMolBond::OnChangeAtom)
   EVT_BUTTON(ID_MOLBOND_ATOM2,    WXMolBond::OnChangeAtom)
   EVT_CHECKBOX(ID_MOLBOND_FREEBUTTON,    WXMolBond::OnToggleFree)
END_EVENT_TABLE()

WXMolBond::WXMolBond(wxWindow *parent, MolBond*obj):
WXCrystObjBasic(parent),mpMolBond(obj),mpButtonFree(0)
{
   VFN_DEBUG_ENTRY("WXMolBond::WXMolBond():"<<obj->GetName(),6)
   mpSizer=new wxBoxSizer(wxHORIZONTAL);
   #if 1
   mpButtonFree=new wxCheckBox(this,ID_MOLBOND_FREEBUTTON,"",wxDefaultPosition, wxDefaultSize);
   mpButtonFree->Fit();
   mpSizer->Add(mpButtonFree,0,wxALIGN_CENTER);
   #endif
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
   
   WXFieldPar<REAL> *delta=
      new WXFieldPar<REAL>(this,",delta=",-1,&(mpMolBond->LengthDelta()));
   mpSizer->Add(delta,0,wxALIGN_CENTER);
   mList.Add(delta);

   WXFieldPar<REAL> *sigma=
      new WXFieldPar<REAL>(this,",sigma=",-1,&(mpMolBond->LengthSigma()));
   mpSizer->Add(sigma,0,wxALIGN_CENTER);
   mList.Add(sigma);

   this->SetSizer(mpSizer);
   this->BottomLayout(0);
   this->CrystUpdate();
   VFN_DEBUG_EXIT("WXMolBond::WXMolBond():"<<obj->GetName(),6)
}
WXMolBond::~WXMolBond()
{
   mpMolBond->WXNotifyDelete();
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
   if(0!=mpButtonFree) mpButtonFree->SetValue(mpMolBond->IsFreeTorsion());
   mList.UpdateUI();
   mpFieldAtom1->SetValue(mpMolBond->GetAtom1().GetName().c_str());
   mpFieldAtom2->SetValue(mpMolBond->GetAtom2().GetName().c_str());
   VFN_DEBUG_EXIT("WXMolBond::UpdateUI()",5)
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

void WXMolBond::OnToggleFree(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXMolBond::OnToggleFree()",6)
   if(0!=mpButtonFree) mpMolBond->SetFreeTorsion(mpButtonFree->GetValue()); 
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
   VFN_DEBUG_ENTRY("WXMolBondAngle::WXMolBond():"<<obj->GetName(),6)
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
   angle->SetHumanValueScale(RAD2DEG);
   mpSizer->Add(angle,0,wxALIGN_CENTER);
   mList.Add(angle);

   WXFieldPar<REAL> *delta=
      new WXFieldPar<REAL>(this,",delta=",WXCRYST_ID(),&(mpMolBondAngle->AngleDelta()));
   delta->SetHumanValueScale(RAD2DEG);
   mpSizer->Add(delta,0,wxALIGN_CENTER);
   mList.Add(delta);

   WXFieldPar<REAL> *sigma=
      new WXFieldPar<REAL>(this,",sigma=",WXCRYST_ID(),&(mpMolBondAngle->AngleSigma()));
   sigma->SetHumanValueScale(RAD2DEG);
   mpSizer->Add(sigma,0,wxALIGN_CENTER);
   mList.Add(sigma);

   this->SetSizer(mpSizer);
   this->BottomLayout(0);
   this->CrystUpdate();
   VFN_DEBUG_EXIT("WXMolBondAngle::WXMolBond():"<<obj->GetName(),6)
}

WXMolBondAngle::~WXMolBondAngle()
{
   mpMolBondAngle->WXNotifyDelete();
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
   angle->SetHumanValueScale(RAD2DEG);
   mpSizer->Add(angle,0,wxALIGN_CENTER);
   mList.Add(angle);

   WXFieldPar<REAL> *delta=
      new WXFieldPar<REAL>(this,",delta=",-1,&(mpMolDihedralAngle->AngleDelta()));
   delta->SetHumanValueScale(RAD2DEG);
   mpSizer->Add(delta,0,wxALIGN_CENTER);
   mList.Add(delta);

   WXFieldPar<REAL> *sigma=
      new WXFieldPar<REAL>(this,",sigma=",-1,&(mpMolDihedralAngle->AngleSigma()));
   sigma->SetHumanValueScale(RAD2DEG);
   mpSizer->Add(sigma,0,wxALIGN_CENTER);
   mList.Add(sigma);

   this->SetSizer(mpSizer);
   this->BottomLayout(0);
   this->CrystUpdate();
   VFN_DEBUG_EXIT("WXMolDihedralAngle::WXMolBond():"<<obj->GetName(),6)
}

WXMolDihedralAngle::~WXMolDihedralAngle()
{
   mpMolDihedralAngle->WXNotifyDelete();
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
WXCRYST_ID ID_MOLECULE_MENU_FORMULA_RIGIDIFY_WITH_DIHEDRALANGLES;
WXCRYST_ID ID_MOLECULE_MENU_FORMULA_TEST;
WXCRYST_ID ID_MOLECULE_MENU_FORMULA_REMOVE_ATOM;
WXCRYST_ID ID_MOLECULE_MENU_FORMULA_REMOVE_BOND;
WXCRYST_ID ID_MOLECULE_MENU_FORMULA_REMOVE_ANGLE;
WXCRYST_ID ID_MOLECULE_MENU_FORMULA_REMOVE_DIHEDRAL;
WXCRYST_ID ID_MOLECULE_MENU_FORMULA_SHOW_BONDLIST;
WXCRYST_ID ID_MOLECULE_MENU_FORMULA_SHOW_BONDANGLELIST;
WXCRYST_ID ID_MOLECULE_MENU_FORMULA_SHOW_DIHEDRALANGLELIST;
WXCRYST_ID ID_MOLECULE_MENU_FORMULA_SET_DELTA_SIGMA;

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
   EVT_MENU(ID_MOLECULE_MENU_FORMULA_RIGIDIFY_WITH_DIHEDRALANGLES,WXMolecule::OnMenuRigidfyWithDihedralAngles)
   EVT_MENU(ID_MOLECULE_MENU_FORMULA_REMOVE_ATOM,         WXMolecule::OnMenuRemoveAtom)
   EVT_MENU(ID_MOLECULE_MENU_FORMULA_REMOVE_BOND,         WXMolecule::OnMenuRemoveBond)
   EVT_MENU(ID_MOLECULE_MENU_FORMULA_REMOVE_ANGLE,        WXMolecule::OnMenuRemoveAngle)
   EVT_MENU(ID_MOLECULE_MENU_FORMULA_REMOVE_DIHEDRAL,     WXMolecule::OnMenuRemoveDihedralAngle)
   EVT_MENU(ID_MOLECULE_MENU_FORMULA_TEST        ,        WXMolecule::OnMenuTest)
   EVT_MENU(ID_MENU_SETLIMITS,                            WXMolecule::OnMenuSetLimits)
   EVT_MENU(ID_MOLECULE_MENU_FORMULA_SHOW_BONDLIST,       WXMolecule::OnMenuShowBondList)
   EVT_MENU(ID_MOLECULE_MENU_FORMULA_SHOW_BONDANGLELIST,  WXMolecule::OnMenuShowBondAngleList)
   EVT_MENU(ID_MOLECULE_MENU_FORMULA_SHOW_DIHEDRALANGLELIST,WXMolecule::OnMenuShowDihedralAngleList)
   EVT_MENU(ID_MOLECULE_MENU_FORMULA_SET_DELTA_SIGMA     ,WXMolecule::OnMenuSetDeltaSigma)
END_EVENT_TABLE()

WXMolecule::WXMolecule(wxWindow *parent, Molecule *mol):
WXScatterer(parent,mol),mpMolecule(mol),
mpSizerBondList(0),mpSizerAngleList(0),mpSizerDihedralAngleList(0),
mpBondWin(0),mpAngleWin(0),mpDihedralAngleWin(0)
{
   VFN_DEBUG_ENTRY("WXMolecule::WXMolecule():"<<mol->GetName(),6)
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
         mpMenuBar->AddMenuItem(ID_MOLECULE_MENU_FORMULA,ID_MOLECULE_MENU_FORMULA_SET_DELTA_SIGMA,
                                "Set Restraints delta && sigma for all bonds && angles");
         mpMenuBar->GetMenu(ID_MOLECULE_MENU_FORMULA).AppendSeparator();
         mpMenuBar->AddMenuItem(ID_MOLECULE_MENU_FORMULA,ID_MOLECULE_MENU_FORMULA_ADD_ATOM,
                                "Add an Atom");
         mpMenuBar->AddMenuItem(ID_MOLECULE_MENU_FORMULA,ID_MOLECULE_MENU_FORMULA_ADD_BOND,
                                "Add Bond Restraint");
         mpMenuBar->AddMenuItem(ID_MOLECULE_MENU_FORMULA,ID_MOLECULE_MENU_FORMULA_ADD_ANGLE,
                                "Add Bond Angle Restraint");
         mpMenuBar->AddMenuItem(ID_MOLECULE_MENU_FORMULA,ID_MOLECULE_MENU_FORMULA_ADD_DIHEDRAL,
                                "Add Dihedral Angle Restraint");
         mpMenuBar->AddMenuItem(ID_MOLECULE_MENU_FORMULA,ID_MOLECULE_MENU_FORMULA_RIGIDIFY_WITH_DIHEDRALANGLES,
                                "Rigidify with Dihedral Angles");
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
         mpMenuBar->AddMenuItem(ID_MOLECULE_MENU_FORMULA,ID_MOLECULE_MENU_FORMULA_SHOW_BONDLIST,
                                "Show Bond List");
         mpMenuBar->AddMenuItem(ID_MOLECULE_MENU_FORMULA,ID_MOLECULE_MENU_FORMULA_SHOW_BONDANGLELIST,
                                "Show Bond Angle List");
         mpMenuBar->AddMenuItem(ID_MOLECULE_MENU_FORMULA,ID_MOLECULE_MENU_FORMULA_SHOW_DIHEDRALANGLELIST,
                                "Show Dihedral Angle List");
         #ifdef __DEBUG__
         mpMenuBar->AddMenuItem(ID_MOLECULE_MENU_FORMULA,ID_MOLECULE_MENU_FORMULA_TEST,"Debug Test");
         #endif
      mpSizer->SetItemMinSize(mpMenuBar,
                              mpMenuBar->GetSize().GetWidth(),
                              mpMenuBar->GetSize().GetHeight());
   
   //sizers
   mpSizerAtomList= new wxBoxSizer(wxVERTICAL);

   mpSizer->Add(mpSizerAtomList,0,wxALIGN_LEFT);
   this->BottomLayout(0);
   this->CrystUpdate();
   VFN_DEBUG_EXIT("WXMolecule::WXMolecule():"<<mol->GetName(),6)
}

WXMolecule::~WXMolecule()
{
   VFN_DEBUG_ENTRY("WXMolecule::~WXMolecule()",10)
   if(0!=mpBondWin) mpBondWin->Destroy();
   if(0!=mpAngleWin) mpAngleWin->Destroy();
   if(0!=mpDihedralAngleWin) mpDihedralAngleWin->Destroy();
   VFN_DEBUG_EXIT("WXMolecule::~WXMolecule()",10)
}

void WXMolecule::OnMenuOptimizeConformation(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_ENTRY("WXMolecule::OnMenuOptimizeConformation()",5)
   WXCrystValidateAllUserInput();
   mpMolecule->OptimizeConformation(100000,mpMolecule->GetAtomList().size());
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
   
   double d;
   if((at1->IsDummy())||(at2->IsDummy())) d=1.5;
   else d= (at1->GetScatteringPower().GetRadius()+at2->GetScatteringPower().GetRadius())*0.9;
   stringstream s;
   s<<d;
   string mes="Enter bond distance (Angstroems) for "+at1->GetName()+"-"+at2->GetName();
   wxTextEntryDialog dialog(this,mes.c_str(),
                           "Bond distance",s.str().c_str(),wxOK | wxCANCEL);
   if(wxID_OK!=dialog.ShowModal())
   {
      VFN_DEBUG_EXIT("WXMolecule::OnMenuAddBond():Canceled",6)
      return;
   }
   dialog.GetValue().ToDouble(&d);
   
   mpMolecule->AddBond(*at1,*at2,d,.01,.02,1.);
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
   
   double a=109.5;
   stringstream s;
   s<<a;
   string mes="Enter bond angle (degrees) for "+at1->GetName()
                                           +"-"+at2->GetName()
                                           +"-"+at3->GetName();
   wxTextEntryDialog dialog(this,mes.c_str(),
                           "Bond angle",s.str().c_str(),wxOK | wxCANCEL);
   if(wxID_OK!=dialog.ShowModal())
   {
      VFN_DEBUG_EXIT("WXMolecule::OnMenuAddAngle():Canceled",6)
      return;
   }
   dialog.GetValue().ToDouble(&a);
   
   mpMolecule->AddBondAngle(*at1,*at2,*at3,a*DEG2RAD,.01,0.02);
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

   double a=180;
   stringstream s;
   s<<a;
   string mes="Enter dihedral angle (degrees) for "+at1->GetName()
                                               +"-"+at2->GetName()
                                               +"-"+at3->GetName()
                                               +"-"+at4->GetName();
   wxTextEntryDialog dialog(this,"Enter dihedral angle (degrees)",
                           "Bond angle",s.str().c_str(),wxOK | wxCANCEL);
   if(wxID_OK!=dialog.ShowModal())
   {
      VFN_DEBUG_EXIT("WXMolecule::OnMenuAddDihedralAngle():Canceled",6)
      return;
   }
   dialog.GetValue().ToDouble(&a);

   mpMolecule->AddDihedralAngle(*at1,*at2,*at3,*at4,a*DEG2RAD,.01,.02);
   VFN_DEBUG_EXIT("WXMolecule::OnMenuAddDihedralAngle()",6)
}

void WXMolecule::OnMenuRemoveAtom(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_ENTRY("WXMolecule::OnMenuRemoveAtom()",6)
   vector<MolAtom*> v=mpMolecule->GetAtomList();
   #if 0
   int choice;
   MolAtom *at=WXDialogChooseFromVector(v,
                               (wxWindow*)this,"Choose the Atom to be removed",choice);
   if(0==at) return;
   mpMolecule->RemoveAtom(*at);
   #else
   list<MolAtom*> vAt=WXDialogChooseMultipleFromVector(v,(wxWindow*)this,
                                                       "Choose the Atom(s) to be removed");
   if(0==vAt.size()) return;
   for(list<MolAtom*>::iterator pos=vAt.begin();pos!=vAt.end();++pos) mpMolecule->RemoveAtom(**pos);
   #endif
   this->CrystUpdate();
   VFN_DEBUG_EXIT("WXMolecule::OnMenuRemoveAtom()",6)
}

void WXMolecule::OnMenuRemoveBond(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_ENTRY("WXMolecule::OnMenuRemoveBond()",6)
   vector<MolBond*> v=mpMolecule->GetBondList();
   #if 0
   int choice;
   MolBond *b=WXDialogChooseFromVector(v,
                               (wxWindow*)this,"Choose the Bond to be removed",choice);
   if(0==b) return;
   mpMolecule->RemoveBond(*b);
   #else
   list<MolBond*> vBond=WXDialogChooseMultipleFromVector(v,(wxWindow*)this,
                                                         "Choose the Bond(s) to be removed");
   if(0==vBond.size()) return;
   const int answer =wxMessageBox
                  ("Remove Bond and Dihedral Angles involving the deleted Bond(s) (if any) ?",
                   "Delete related Restraints ?",wxYES_NO, this);
   for(list<MolBond*>::iterator pos=vBond.begin();pos!=vBond.end();++pos)
   {
      if(answer==wxYES)
      {
         const MolAtom *pAtom1= &((*pos)->GetAtom1());
         const MolAtom *pAtom2= &((*pos)->GetAtom2());
         for(vector<MolBondAngle*>::iterator posb=mpMolecule->GetBondAngleList().begin();
             posb!=mpMolecule->GetBondAngleList().end();++posb)
         {
            if(  ( (pAtom1==&((*posb)->GetAtom1())) && (pAtom2==&((*posb)->GetAtom2())) )
               ||( (pAtom1==&((*posb)->GetAtom2())) && (pAtom2==&((*posb)->GetAtom1())) )
               ||( (pAtom1==&((*posb)->GetAtom2())) && (pAtom2==&((*posb)->GetAtom3())) )
               ||( (pAtom1==&((*posb)->GetAtom3())) && (pAtom2==&((*posb)->GetAtom2())) ))
            {
               posb=mpMolecule->RemoveBondAngle(**posb);
               --posb;
            }
         }
         for(vector<MolDihedralAngle*>::iterator posb=mpMolecule->GetDihedralAngleList().begin();
             posb!=mpMolecule->GetDihedralAngleList().end();++posb)
         {
            if(  ( (pAtom1==&((*posb)->GetAtom1())) && (pAtom2==&((*posb)->GetAtom2())) )
               ||( (pAtom1==&((*posb)->GetAtom2())) && (pAtom2==&((*posb)->GetAtom1())) )
               ||( (pAtom1==&((*posb)->GetAtom2())) && (pAtom2==&((*posb)->GetAtom3())) )
               ||( (pAtom1==&((*posb)->GetAtom3())) && (pAtom2==&((*posb)->GetAtom2())) )
               ||( (pAtom1==&((*posb)->GetAtom3())) && (pAtom2==&((*posb)->GetAtom4())) )
               ||( (pAtom1==&((*posb)->GetAtom4())) && (pAtom2==&((*posb)->GetAtom3())) ))
            {
               posb=mpMolecule->RemoveDihedralAngle(**posb);
               --posb;
            }
         }
      }
      mpMolecule->RemoveBond(**pos);
   }
   #endif
   this->CrystUpdate();
   VFN_DEBUG_EXIT("WXMolecule::OnMenuRemoveBond()",6)
}

void WXMolecule::OnMenuRemoveAngle(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_ENTRY("WXMolecule::OnMenuRemoveAngle()",6)
   vector<MolBondAngle*> v=mpMolecule->GetBondAngleList();
   #if 0
   int choice;
   MolBondAngle *a=WXDialogChooseFromVector(v,
                               (wxWindow*)this,"Choose the Bond Angle to be removed",choice);
   if(0==a) return;
   mpMolecule->RemoveBondAngle(*a);
   #else
   list<MolBondAngle*> vAngle=WXDialogChooseMultipleFromVector(v,(wxWindow*)this,
                                                              "Choose the Bond Angle(s) to be removed");
   if(0==vAngle.size()) return;
   for(list<MolBondAngle*>::iterator pos=vAngle.begin();pos!=vAngle.end();++pos)
      mpMolecule->RemoveBondAngle(**pos);
   #endif
   this->CrystUpdate();
   VFN_DEBUG_EXIT("WXMolecule::OnMenuRemoveAngle()",6)
}

void WXMolecule::OnMenuRemoveDihedralAngle(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_ENTRY("WXMolecule::OnMenuRemoveDihedralAngle()",6)
   vector<MolDihedralAngle*> v=mpMolecule->GetDihedralAngleList();
   #if 0
   int choice;
   MolDihedralAngle *a=WXDialogChooseFromVector(v,
                               (wxWindow*)this,"Choose the Dihedral Angle to be removed",choice);
   if(0==a) return;
   mpMolecule->RemoveDihedralAngle(*a);
   #else
   list<MolDihedralAngle*> vAngle=WXDialogChooseMultipleFromVector(v,(wxWindow*)this,
                                     "Choose the Dihedral Angle(s) to be removed");
   if(0==vAngle.size()) return;
   for(list<MolDihedralAngle*>::iterator pos=vAngle.begin();pos!=vAngle.end();++pos)
      mpMolecule->RemoveDihedralAngle(**pos);
   #endif
   this->CrystUpdate();
   VFN_DEBUG_EXIT("WXMolecule::OnMenuRemoveDihedralAngle()",6)
}

void WXMolecule::OnMenuTest(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_ENTRY("WXMolecule::OnMenuTest()",6)
   #if 0
   mpMolecule->BuildRotorGroup();
   mpMolecule->BuildFlipGroup();
   mpMolecule->FlipAtomGroup(*(mpMolecule->mvFlipGroup.begin()));
   mpMolecule->GetCrystal().UpdateDisplay();
   #endif
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
      bool needLayoutAtom=false;
      bool needLayoutBond=false;
      bool needLayoutAngle=false;
      bool needLayoutDihed=false;
      //Remove any atom, bond, bond angle or dihedral angle that could have been removed
      {
         map<MolAtom*,WXCrystObjBasic*>::iterator pos;
         for(pos=mvpAtom.begin();pos!=mvpAtom.end();)
         {
            vector<MolAtom*>::const_iterator pos2=find(mpMolecule->GetAtomList().begin(),
                                                       mpMolecule->GetAtomList().end(),pos->first);
            if(pos2==mpMolecule->GetAtomList().end())
            {
               mList.Remove(pos->second);
               mvpAtom.erase(pos);
               pos=mvpAtom.begin();
               needLayoutAtom=true;
            }else ++pos;
         }
      }
      if(0!=mpBondWin)
      {
         map<MolBond*,WXCrystObjBasic*>::iterator pos;
         for(pos=mvpBond.begin();pos!=mvpBond.end();)
         {
            vector<MolBond*>::const_iterator pos2=find(mpMolecule->GetBondList().begin(),
                                                       mpMolecule->GetBondList().end(),pos->first);
            if(pos2==mpMolecule->GetBondList().end())
            {
               mList.Remove(pos->second);
               mvpBond.erase(pos);
               pos=mvpBond.begin();
               needLayoutBond=true;
            }else ++pos;
         }
      }
      if(0!=mpAngleWin)
      {
         map<MolBondAngle*,WXCrystObjBasic*>::iterator pos;
         for(pos=mvpBondAngle.begin();pos!=mvpBondAngle.end();)
         {
            vector<MolBondAngle*>::const_iterator pos2=
                                    find(mpMolecule->GetBondAngleList().begin(),
                                         mpMolecule->GetBondAngleList().end(),pos->first);
            if(pos2==mpMolecule->GetBondAngleList().end())
            {
               mList.Remove(pos->second);
               mvpBondAngle.erase(pos);
               pos=mvpBondAngle.begin();
               needLayoutAngle=true;
            }else ++pos;
         }
      }
      if(0!=mpDihedralAngleWin)
      {
         map<MolDihedralAngle*,WXCrystObjBasic*>::iterator pos;
         for(pos=mvpDihedralAngle.begin();pos!=mvpDihedralAngle.end();)
         {
            vector<MolDihedralAngle*>::const_iterator pos2=
                              find(mpMolecule->GetDihedralAngleList().begin(),
                                   mpMolecule->GetDihedralAngleList().end(),pos->first);
            if(pos2==mpMolecule->GetDihedralAngleList().end())
            {
               mList.Remove(pos->second);
               mvpDihedralAngle.erase(pos);
               pos=mvpDihedralAngle.begin();
               needLayoutDihed=true;
            }else ++pos;
         }
      }
      //Add any atom, bond, bond angle or dihedral angle that could have been added
      {
         vector<MolAtom*>::iterator pos;
         for(pos=mpMolecule->GetAtomList().begin();pos!=mpMolecule->GetAtomList().end();pos++)
         {
            map<MolAtom*,WXCrystObjBasic*>::const_iterator pos2=mvpAtom.find(*pos);
            if(pos2==mvpAtom.end())
            {
               VFN_DEBUG_MESSAGE("WXMolecule::CrystUpdate():Atom not found:"<<(*pos)->GetName(),5)
               WXCrystObjBasic *at=(*pos)->WXCreate(this);
               mpSizerAtomList->Add(at);
               mList.Add(at);
               mvpAtom.insert(make_pair(*pos,(*pos)->WXGet()));
               needLayoutAtom=true;
            }
         }
      }
      if(0!=mpBondWin)
      {
         vector<MolBond*>::iterator pos;
         for(pos=mpMolecule->GetBondList().begin();pos!=mpMolecule->GetBondList().end();pos++)
         {
            map<MolBond*,WXCrystObjBasic*>::const_iterator pos2=mvpBond.find(*pos);
            if(pos2==mvpBond.end())
            {
               VFN_DEBUG_MESSAGE("WXMolecule::CrystUpdate():Bond not found",10)
               WXCrystObjBasic *b=(*pos)->WXCreate(mpBondWin);
               mpSizerBondList->Add(b);
               mList.Add(b);
               mvpBond.insert(make_pair(*pos,(*pos)->WXGet()));
               needLayoutBond=true;
            }
         }
      }
      if(0!=mpAngleWin)
      {
         vector<MolBondAngle*>::iterator pos;
         for(pos=mpMolecule->GetBondAngleList().begin();
             pos!=mpMolecule->GetBondAngleList().end();pos++)
         {
            map<MolBondAngle*,WXCrystObjBasic*>::const_iterator pos2=mvpBondAngle.find(*pos);
            if(pos2==mvpBondAngle.end())
            {
               VFN_DEBUG_MESSAGE("WXMolecule::CrystUpdate():Bond Angle not found",5)
               WXCrystObjBasic *b=(*pos)->WXCreate(mpAngleWin);
               mpSizerAngleList->Add(b);
               mList.Add(b);
               mvpBondAngle.insert(make_pair(*pos,(*pos)->WXGet()));
               needLayoutAngle=true;
            }
         }
      }
      if(0!=mpDihedralAngleWin)
      {
         vector<MolDihedralAngle*>::iterator pos;
         for(pos=mpMolecule->GetDihedralAngleList().begin();
             pos!=mpMolecule->GetDihedralAngleList().end();pos++)
         {
            map<MolDihedralAngle*,WXCrystObjBasic*>::const_iterator pos2
                                                         =mvpDihedralAngle.find(*pos);
            if(pos2==mvpDihedralAngle.end())
            {
               VFN_DEBUG_MESSAGE("WXMolecule::CrystUpdate():Bond not found",5)
               WXCrystObjBasic *b=(*pos)->WXCreate(mpDihedralAngleWin);
               mpSizerDihedralAngleList->Add(b);
               mList.Add(b);
               mvpDihedralAngle.insert(make_pair(*pos,(*pos)->WXGet()));
               needLayoutDihed=true;
            }
         }
      }
      if(needLayoutAtom)
      {
         mpTopSizer->SetSizeHints(this);
         this->Layout();
      }
      if(needLayoutBond&&(0!=mpBondWin))
      {
         VFN_DEBUG_MESSAGE("WXMolecule::CrystUpdate(): BondWin Layout",6)
         mpBondWin->Scroll(0,0);
         int width0,height0;
         mpBondWin->GetVirtualSize(&width0,&height0);
         mpSizerBondList->Layout();
         wxSize s=mpSizerBondList->GetMinSize();
         if((height0!=s.GetHeight())||(width0!=s.GetWidth()))
            mpBondWin->SetScrollbars(40,40,s.GetWidth()/40+1,s.GetHeight()/40+1);
         mpBondWin->GetParent()->SetSize(-1,s.GetWidth());
      }
      if(needLayoutAngle&&(0!=mpAngleWin))
      {
         VFN_DEBUG_MESSAGE("WXMolecule::CrystUpdate(): AngleWin Layout",6)
         mpAngleWin->Scroll(0,0);
         int width0,height0;
         mpAngleWin->GetVirtualSize(&width0,&height0);
         //mpSizerAngleList->RecalcSizes();
         mpSizerAngleList->Layout();
         wxSize s=mpSizerAngleList->GetMinSize();
         if((height0!=s.GetHeight())||(width0!=s.GetWidth()))
            mpAngleWin->SetScrollbars(40,40,s.GetWidth()/40+1,s.GetHeight()/40+1);
         mpAngleWin->GetParent()->SetSize(-1,s.GetWidth());
      }
      if(needLayoutDihed&&(0!=mpDihedralAngleWin))
      {
         VFN_DEBUG_MESSAGE("WXMolecule::CrystUpdate(): DihedralAngleWin Layout",6)
         mpDihedralAngleWin->Scroll(0,0);
         int width0,height0;
         mpDihedralAngleWin->GetVirtualSize(&width0,&height0);
         mpSizerDihedralAngleList->Layout();
         wxSize s=mpSizerDihedralAngleList->GetMinSize();
         if((height0!=s.GetHeight())||(width0!=s.GetWidth()))
            mpDihedralAngleWin->SetScrollbars(40,40,s.GetWidth()/40+1,s.GetHeight()/40+1);
         mpDihedralAngleWin->GetParent()->SetSize(-1,s.GetWidth());
      }
   }
   this->WXRefinableObj::CrystUpdate();
   VFN_DEBUG_EXIT("WXMolecule::CrystUpdate()",6)
}
void WXMolecule::OnMenuShowBondList(wxCommandEvent &event)
{
   if(0!=mpBondWin) return;
   wxFrame *frame= new wxFrame(this,-1,("Bond List for: "+mpMolecule->GetName()).c_str(),
                               wxDefaultPosition,wxSize(250,200));
   mpBondWin = new WXMolScrolledWindow(frame,this);
   mpSizerBondList= new wxBoxSizer(wxVERTICAL);
   frame->Show(true);
   this->CrystUpdate();
}
void WXMolecule::OnMenuShowBondAngleList(wxCommandEvent &event)
{
   if(0!=mpAngleWin) return;
   wxFrame *frame= new wxFrame(this,-1,("Bond Angles list for: "+mpMolecule->GetName()).c_str(),
                               wxDefaultPosition,wxSize(250,200));
   mpAngleWin = new WXMolScrolledWindow(frame,this);
   mpSizerAngleList= new wxBoxSizer(wxVERTICAL);
   frame->Show(true);
   this->CrystUpdate();
}
void WXMolecule::OnMenuShowDihedralAngleList(wxCommandEvent &event)
{
   if(0!=mpDihedralAngleWin) return;
   wxFrame *frame= new wxFrame(this,-1,("Dihedral Angles list for: "+mpMolecule->GetName()).c_str(),
                               wxDefaultPosition,wxSize(250,200));
   mpDihedralAngleWin = new WXMolScrolledWindow(frame,this);
   mpSizerDihedralAngleList= new wxBoxSizer(wxVERTICAL);
   frame->Show(true);
   this->CrystUpdate();
}

void WXMolecule::OnMenuRigidfyWithDihedralAngles(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_ENTRY("WXMolecule::OnMenuRigidfyWithDihedralAngles()",6)
   wxString msg;
   msg.Printf( _T("This will add all possible dihedral angles,\n")
               _T("in practice making the Molecule rigid\n\n")
               _T("Are you sure you want to proceed ?")
              );

   wxMessageDialog w(this,msg,"Warning !",wxYES_NO|wxNO_DEFAULT|wxICON_EXCLAMATION);
   int result=w.ShowModal();
   if(wxID_YES==result) mpMolecule->RigidifyWithDihedralAngles();
   
   VFN_DEBUG_EXIT("WXMolecule::OnMenuRigidfyWithDihedralAngles()",6)
}

void WXMolecule::OnMenuSetDeltaSigma(wxCommandEvent &event)
{
   VFN_DEBUG_ENTRY("WXMolecule::OnMenuSetDeltaSigma()",6)
   WXCrystValidateAllUserInput();
   double sigma=0.01,delta=0.02;
   {
      stringstream s;
      s<<delta;
      string title="Choose 'delta' value";
      wxString info;
      info.Printf(_T("The 'delta' value is the allowed range \n")
                  _T("(without penalty) around the expected value.\n\n")
                  _T("It is by default equal to 0.02, in Angstroems for bond lengths,\n")
                  _T("and in radians for angles (0.02rad = 1.15°)\n\n")
                  _T("DO NOT TRY TO CHANGE THE DEFAULT VALUE\n")
                  _T("UNLESS YOU REALLY KNOW WHAT YOU ARE DOING\n")
                  _T("Fox has been optimized with the default values...")
                 );
      wxTextEntryDialog dialog(this,info,title.c_str(),s.str().c_str(),wxOK|wxCANCEL);
      if(wxID_OK!=dialog.ShowModal())
      {
         VFN_DEBUG_EXIT("WXMolecule::OnMenuSetDeltaSigma():Canceled",6)
         return;
      }
      dialog.GetValue().ToDouble(&delta);
   }
   {
      stringstream s;
      s<<sigma;
      string title="Choose 'sigma' value";
      wxString info;
      info.Printf(_T("The 'sigma' value is used to compute \n")
                  _T("penalty)around the expected value\n\n")
                  _T("It is by default equal to 0.01, in Angstroems for bond angles,\n")
                  _T("and in radians for angles (0.01rad = 0.57°)\n\n")
                  _T("DO NOT TRY TO CHANGE THE DEFAULT VALUE\n")
                  _T("UNLESS YOU REALLY KNOW WHAT YOU ARE DOING\n")
                  _T("Fox has been optimized with the default values...")
                 );
      wxTextEntryDialog dialog(this,info,title.c_str(),s.str().c_str(),wxOK|wxCANCEL);
      if(wxID_OK!=dialog.ShowModal())
      {
         VFN_DEBUG_EXIT("WXMolecule::OnMenuSetDeltaSigma():Canceled",6)
         return;
      }
      dialog.GetValue().ToDouble(&sigma);
   }
   for(vector<MolBond*>::iterator pos=mpMolecule->GetBondList().begin();
       pos != mpMolecule->GetBondList().end();++pos)
   {
      (*pos)->SetLengthDelta(delta);
      (*pos)->SetLengthSigma(sigma);
   }
   for(vector<MolBondAngle*>::iterator pos=mpMolecule->GetBondAngleList().begin();
       pos != mpMolecule->GetBondAngleList().end();++pos)
   {
      (*pos)->AngleDelta()=delta;
      (*pos)->AngleSigma()=sigma;
   }
   for(vector<MolDihedralAngle*>::iterator pos=mpMolecule->GetDihedralAngleList().begin();
       pos != mpMolecule->GetDihedralAngleList().end();++pos)
   {
      (*pos)->AngleDelta()=delta;
      (*pos)->AngleSigma()=sigma;
   }
   mpMolecule->GetBondListClock().Click();
   this->CrystUpdate();
   VFN_DEBUG_EXIT("WXMolecule::OnMenuSetDeltaSigma()",6)
}

void WXMolecule::NotifyDeleteListWin(WXMolScrolledWindow *win)
{
   VFN_DEBUG_ENTRY("WXMolecule::NotifyDeleteListWin()",6)
   if(win==mpBondWin)
   {
      VFN_DEBUG_MESSAGE("WXMolecule::NotifyDeleteListWin(): Bond List window",5)
      mpBondWin=0;
      for(map<MolBond*,WXCrystObjBasic*>::iterator pos=mvpBond.begin();pos!=mvpBond.end();pos++)
         mList.Remove(pos->second);
      mvpBond.clear();
   }
   if(win==mpAngleWin)
   {
      VFN_DEBUG_MESSAGE("WXMolecule::NotifyDeleteListWin(): Angle List window",5)
      mpAngleWin=0;
      for(map<MolBondAngle*,WXCrystObjBasic*>::iterator pos=mvpBondAngle.begin();
          pos!=mvpBondAngle.end();pos++) mList.Remove(pos->second);
      mvpBondAngle.clear();
   }
   if(win==mpDihedralAngleWin)
   {
      VFN_DEBUG_MESSAGE("WXMolecule::NotifyDeleteListWin(): Dihedral Angle List window",5)
      mpDihedralAngleWin=0;
      for(map<MolDihedralAngle*,WXCrystObjBasic*>::iterator pos=mvpDihedralAngle.begin();
          pos!=mvpDihedralAngle.end();pos++) mList.Remove(pos->second);
      mvpDihedralAngle.clear();
   }
   VFN_DEBUG_EXIT("WXMolecule::NotifyDeleteListWin()",6)
}
void WXMolecule::UpdateUI()
{
   VFN_DEBUG_ENTRY("WXMolecule::UpdateUI()",5)
   this->WXRefinableObj::UpdateUI();
   VFN_DEBUG_EXIT("WXMolecule::UpdateUI()",5)
}
} //namespace
