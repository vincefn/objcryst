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

#include <stdio.h> //for sprintf()
#include <fstream>

// wx headers, with or without precompilation
#include "wx/wxprec.h"
#ifdef __BORLANDC__
    #pragma hdrstop
#endif
#ifndef WX_PRECOMP
    #include "wx/wx.h"
#endif

#include "wxCryst/wxZScatterer.h"
#include "ObjCryst/Molecule.h"

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
/// Conversion from ZScatterer to the newer Molecule object.
Molecule *ZScatterer2Molecule(ZScatterer *scatt);

Molecule *ZScatterer2Molecule(ZScatterer *scatt)
{
   VFN_DEBUG_ENTRY("ZScatterer2Molecule()",6)
   Molecule *mol=new Molecule(scatt->GetCrystal(),scatt->GetName());
   const unsigned long nb=scatt->GetZAtomRegistry().GetNb();
   const ScatteringComponentList *l=&(scatt->GetScatteringComponentList());
   REAL x,y,z;
   REAL x0=0,y0=0,z0=0;
   for(unsigned int i=0;i<nb;++i)
   {
      x=(*l)(i).mX;
      y=(*l)(i).mY;
      z=(*l)(i).mZ;
      x0+=x;
      y0+=y;
      z0+=z;
      scatt->GetCrystal().FractionalToOrthonormalCoords(x,y,z);
      mol->AddAtom(x,y,z,scatt->GetZAtomRegistry().GetObj(i).GetScatteringPower(),
                   scatt->GetComponentName(i));
      
      if(i>0)
      {
         const RefinablePar* pLength=&(scatt->GetPar(&(scatt->GetZAtomRegistry()
                                                       .GetObj(i).GetZBondLength())));
         if(pLength->IsFixed())
            mol->AddBond(mol->GetAtom(i),mol->GetAtom(scatt->GetZBondAtom(i)),
                         pLength->GetValue(),.01,.05);
         else
            if(pLength->IsLimited())
               mol->AddBond(mol->GetAtom(i),mol->GetAtom(scatt->GetZBondAtom(i)),
                            (pLength->GetMin()+pLength->GetMax())/2.,.01,.05);
      }
      if(i>1)
      {
         const RefinablePar* pAngle=&(scatt->GetPar(&(scatt->GetZAtomRegistry()
                                                      .GetObj(i).GetZAngle())));
         if(pAngle->IsFixed())
            mol->AddBondAngle(mol->GetAtom(i),mol->GetAtom(scatt->GetZBondAtom(i)),
                              mol->GetAtom(scatt->GetZAngleAtom(i)),
                              pAngle->GetValue(),.01,.05);
         else
            if(pAngle->IsLimited())
               mol->AddBondAngle(mol->GetAtom(i),mol->GetAtom(scatt->GetZBondAtom(i)),
                                 mol->GetAtom(scatt->GetZAngleAtom(i)),
                                 (pAngle->GetMin()+pAngle->GetMax())/2.,.01,.05);
      }
      if(i>2)
      {
         const RefinablePar* pDihed=&(scatt->GetPar(&(scatt->GetZAtomRegistry()
                                                      .GetObj(i).GetZDihedralAngle())));
         MolAtom *p1=&(mol->GetAtom(i));
         MolAtom *p2=&(mol->GetAtom(scatt->GetZBondAtom(i)));
         MolAtom *p3=&(mol->GetAtom(scatt->GetZAngleAtom(i)));
         MolAtom *p4=&(mol->GetAtom(scatt->GetZDihedralAngleAtom(i)));
         if(  (abs(GetBondAngle(*p1,*p2,*p3)-M_PI)>0.3)
            &&(abs(GetBondAngle(*p1,*p2,*p4)-M_PI)>0.3)
            &&(abs(GetBondAngle(*p1,*p3,*p4)-M_PI)>0.3)
            &&(abs(GetBondAngle(*p2,*p3,*p4)-M_PI)>0.3))
         {
            if(pDihed->IsFixed())
               mol->AddDihedralAngle(*p1,*p2,*p3,*p4,pDihed->GetValue(),.01,.05);
            else
               if(((pDihed->GetMax()-pDihed->GetMax())<0.3)&&(i>2)&&(pDihed->IsLimited()))
                  mol->AddDihedralAngle(*p1,*p2,*p3,*p4,
                                        (pDihed->GetMin()+pDihed->GetMax())/2.,.01,.05);
         }
      }
      mol->GetAtom(i).SetOccupancy(scatt->GetZAtomRegistry().GetObj(i).GetOccupancy());
      
   }
   for(unsigned int i=0;i<nb;++i)
      for(unsigned int j=i+1;j<nb;++j)
      {
         const REAL dist=GetBondLength(mol->GetAtom(i),mol->GetAtom(j));
         const vector<MolBond*>::const_iterator pos=mol->FindBond(mol->GetAtom(i),mol->GetAtom(j));
         
         if(   (dist<(1.10*( mol->GetAtom(i).GetScatteringPower().GetRadius()
                            +mol->GetAtom(j).GetScatteringPower().GetRadius())))
             &&(pos==mol->GetBondList().end()))
            mol->AddBond(mol->GetAtom(i),mol->GetAtom(j),dist,.01,.05);
      }

   mol->SetX(x0/nb);
   mol->SetY(y0/nb);
   mol->SetZ(z0/nb);
   return mol;
   VFN_DEBUG_EXIT("ZScatterer2Molecule()",6)
}

////////////////////////////////////////////////////////////////////////
//
//    WXZAtom
//
////////////////////////////////////////////////////////////////////////
static const long ID_ZATOM_NAME=    WXCRYST_ID();
static const long ID_ZATOM_SCATTPOW=WXCRYST_ID();
static const long ID_ZATOM_BOND=    WXCRYST_ID();
static const long ID_ZATOM_ANGLE=   WXCRYST_ID();
static const long ID_ZATOM_DIHED=   WXCRYST_ID();

BEGIN_EVENT_TABLE(WXZAtom,wxWindow)
   EVT_BUTTON(ID_ZATOM_SCATTPOW,    WXZAtom::OnChangeScattPow)
END_EVENT_TABLE()

WXZAtom::WXZAtom(wxWindow* parent, ZAtom *obj):
WXCrystObjBasic(parent),mpZAtom(obj)
{
   VFN_DEBUG_ENTRY("WXZAtom::WXZAtom()",6)
   mpSizer=new wxBoxSizer(wxHORIZONTAL);
      
   mpFieldName=new WXFieldString(this, mpZAtom->mName,ID_ZATOM_NAME,80,true);
   mpSizer->Add(mpFieldName,0,wxALIGN_LEFT);
   mpFieldScattPower=new WXFieldChoice(this,ID_ZATOM_SCATTPOW,"Type:",60);
   mpSizer->Add(mpFieldScattPower,0,wxALIGN_LEFT);
   mList.Add(mpFieldScattPower);
#if 1
   if(0<mpZAtom->GetZScatterer().GetZAtomRegistry().Find(*mpZAtom))
   {
      WXFieldRefPar* pFieldBond  =new WXFieldRefPar(this,
               "Bond(w/"+
               (mpZAtom->GetZScatterer().
                  GetZAtomRegistry().GetObj(mpZAtom->GetZBondAtom())).GetName()
               +")",
               &(mpZAtom->GetZScatterer().GetPar(&(mpZAtom->mBondLength))));
      mpSizer->Add(pFieldBond,0,wxALIGN_LEFT);
      mList.Add(pFieldBond);
   }
   
   if(1<mpZAtom->GetZScatterer().GetZAtomRegistry().Find(*mpZAtom))
   {
      WXFieldRefPar* pFieldAngle =new WXFieldRefPar(this,
               "Angle(w/"+
               (mpZAtom->GetZScatterer().
                  GetZAtomRegistry().GetObj(mpZAtom->GetZAngleAtom())).GetName()
               +")",
                        &(mpZAtom->GetZScatterer().GetPar(&(mpZAtom->mAngle))));
      mpSizer->Add(pFieldAngle,0,wxALIGN_LEFT);
      mList.Add(pFieldAngle);
   }
   if(2<mpZAtom->GetZScatterer().GetZAtomRegistry().Find(*mpZAtom))
   {
      WXFieldRefPar* pFieldDihed =new WXFieldRefPar(this,
               "DihedralAngle(w/"+
               (mpZAtom->GetZScatterer().
                  GetZAtomRegistry().GetObj(mpZAtom->GetZDihedralAngleAtom())).GetName()
               +")",
                        &(mpZAtom->GetZScatterer().GetPar(&(mpZAtom->mDihed))));
      mpSizer->Add(pFieldDihed,0,wxALIGN_LEFT);
      mList.Add(pFieldDihed);
   }
   {
      WXFieldRefPar* pFieldOccup =new WXFieldRefPar(this,
               "Occup.",&(mpZAtom->GetZScatterer().GetPar(&(mpZAtom->mOccupancy))));
      mpSizer->Add(pFieldOccup,0,wxALIGN_LEFT);
      mList.Add(pFieldOccup);
   }
#else
   if(0<mpZAtom->GetZScatterer().GetZAtomRegistry().Find(*mpZAtom))
   {
      WXCrystObjBasic* pFieldBond 
         =mpZAtom->GetZScatterer().GetPar(&(mpZAtom->mBondLength)).WXCreate(this);
      mpSizer->Add(pFieldBond,0,wxALIGN_LEFT);
      mList.Add(pFieldBond);
   }
   
   if(1<mpZAtom->GetZScatterer().GetZAtomRegistry().Find(*mpZAtom))
   {
       WXCrystObjBasic* pFieldAngle
          =mpZAtom->GetZScatterer().GetPar(&(mpZAtom->mAngle)).WXCreate(this);
      mpSizer->Add(pFieldAngle,0,wxALIGN_LEFT);
      mList.Add(pFieldAngle);
   }
   if(2<mpZAtom->GetZScatterer().GetZAtomRegistry().Find(*mpZAtom))
   {
       WXCrystObjBasic* pFieldDihed 
       =mpZAtom->GetZScatterer().GetPar(&(mpZAtom->mDihed)).WXCreate(this);
      mpSizer->Add(pFieldDihed,0,wxALIGN_LEFT);
      mList.Add(pFieldDihed);
   }
   {
      WXFieldRefPar* pFieldOccup
         =mpZAtom->GetZScatterer().GetPar(&(mpZAtom->mOccupancy)).WXCreate(this);
      mpSizer->Add(pFieldOccup,0,wxALIGN_LEFT);
      mList.Add(pFieldOccup);
   }
#endif   
   this->SetSizer(mpSizer);
   
   this->CrystUpdate();
   this->Layout();
   VFN_DEBUG_EXIT("WXZAtom::WXZAtom()",6)
}

void WXZAtom::CrystUpdate()
{
   VFN_DEBUG_ENTRY("WXZAtom::CrystUpdate()",6)
   mList.CrystUpdate();
   VFN_DEBUG_EXIT("WXZAtom::CrystUpdate()",6)
}

void WXZAtom::UpdateUI()
{
   VFN_DEBUG_ENTRY("WXZAtom::UpdateUI()",6)
   mList.UpdateUI();
   mpFieldName->SetValue(mpZAtom->GetName().c_str());
   if(0!=mpZAtom->GetScatteringPower())
      mpFieldScattPower->SetValue(mpZAtom->GetScatteringPower()->GetName());
   else
      mpFieldScattPower->SetValue("Dummy");
   VFN_DEBUG_EXIT("WXZAtom::UpdateUI()",6)
}

bool WXZAtom::Layout()
{
   VFN_DEBUG_ENTRY("WXZAtom::Layout()",3)
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
   VFN_DEBUG_EXIT("WXZAtom::Layout()",3)
   return this->wxWindow::Layout();
}

void WXZAtom::OnChangeScattPow(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXAtom::OnChangeScattPow()",6)
   WXCrystValidateAllUserInput();
   int choice;
   const ScatteringPower *scatt=WXDialogChooseFromRegistry(
               mpZAtom->GetZScatterer().GetCrystal().GetScatteringPowerRegistry(),
               (wxWindow*)this,"Choose a new Scattering Power",choice);
   if(0==scatt) return;
   mpZAtom->SetScatteringPower(scatt);
   this->CrystUpdate();
   this->UpdateUI();
}
////////////////////////////////////////////////////////////////////////
//
//    WXZScatterer
//
////////////////////////////////////////////////////////////////////////
static const long ID_ZSCATTERER_MENU_ATOM=                  WXCRYST_ID();
static const long ID_ZSCATTERER_MENU_ATOM_ADD=              WXCRYST_ID();
static const long ID_ZSCATTERER_MENU_ATOM_CHANGE_PIVOT=     WXCRYST_ID();
static const long ID_ZSCATTERER_MENU_PAR_LIMITS_RELAT_BOND= WXCRYST_ID();
static const long ID_ZSCATTERER_MENU_PAR_LIMITS_RELAT_ANGLE=WXCRYST_ID();
static const long ID_ZSCATTERER_MENU_PAR_LIMITS_RELAT_DIHED=WXCRYST_ID();
static const long ID_ZSCATTERER_MENU_FILE=                  WXCRYST_ID();
static const long ID_ZSCATTERER_MENU_IMPORT_FHZ=            WXCRYST_ID();
static const long ID_ZSCATTERER_MENU_EXPORT_FHZ=            WXCRYST_ID();
static const long ID_ZSCATTERER_MENU_CONVERT2MOLECULE=     WXCRYST_ID();
BEGIN_EVENT_TABLE(WXZScatterer,wxWindow)
   EVT_BUTTON(ID_WXOBJ_COLLAPSE,                       WXCrystObj::OnToggleCollapse)
   EVT_MENU(ID_REFOBJ_MENU_PAR_FIXALL,                 WXRefinableObj::OnMenuFixAllPar)
   EVT_MENU(ID_REFOBJ_MENU_PAR_UNFIXALL,               WXRefinableObj::OnMenuUnFixAllPar)
   EVT_MENU(ID_REFOBJ_MENU_PAR_RANDOMIZE,              WXRefinableObj::OnMenuParRandomize)
   EVT_MENU(ID_ZSCATTERER_MENU_PAR_LIMITS_RELAT_BOND,  WXZScatterer::OnMenuSetLimits)
   EVT_MENU(ID_ZSCATTERER_MENU_PAR_LIMITS_RELAT_ANGLE, WXZScatterer::OnMenuSetLimits)
   EVT_MENU(ID_ZSCATTERER_MENU_PAR_LIMITS_RELAT_DIHED, WXZScatterer::OnMenuSetLimits)
   EVT_MENU(ID_ZSCATTERER_MENU_ATOM_ADD,               WXZScatterer::OnMenuAddZAtom)
   EVT_MENU(ID_ZSCATTERER_MENU_ATOM_CHANGE_PIVOT,      WXZScatterer::OnMenuChangePivotAtom)
   EVT_MENU(ID_ZSCATTERER_MENU_IMPORT_FHZ,             WXZScatterer::OnMenuImportZMatrix)
   EVT_MENU(ID_ZSCATTERER_MENU_EXPORT_FHZ,             WXZScatterer::OnMenuExportZMatrix)
   EVT_MENU(ID_ZSCATTERER_MENU_CONVERT2MOLECULE,       WXZScatterer::OnMenuConvert2Molecule)
END_EVENT_TABLE()

WXZScatterer::WXZScatterer(wxWindow* parent, ZScatterer *obj):
WXScatterer(parent,obj),mpZScatterer(obj)
{
   VFN_DEBUG_MESSAGE("WXZScatterer::WXZScatterer()",6)
   //Menus
      mpMenuBar->AddMenu("Import/Export",ID_ZSCATTERER_MENU_FILE);
         mpMenuBar->AddMenuItem(ID_ZSCATTERER_MENU_FILE,
                                ID_ZSCATTERER_MENU_IMPORT_FHZ,
                                "Import Fenske-Hall Zmatrix");
         mpMenuBar->AddMenuItem(ID_ZSCATTERER_MENU_FILE,
                                ID_ZSCATTERER_MENU_EXPORT_FHZ,
                                "Save as Fenske-Hall Zmatrix");
         mpMenuBar->AddMenuItem(ID_ZSCATTERER_MENU_FILE,
                                ID_ZSCATTERER_MENU_CONVERT2MOLECULE,
                                "Convert to Molecule");
      mpMenuBar->AddMenu("Parameters",ID_REFOBJ_MENU_PAR);
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_PAR,ID_REFOBJ_MENU_PAR_FIXALL,"Fix all");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_PAR,ID_REFOBJ_MENU_PAR_UNFIXALL,"Unfix all");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_PAR,ID_REFOBJ_MENU_PAR_RANDOMIZE,
                                "Randomize Configuration");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_PAR,ID_ZSCATTERER_MENU_PAR_LIMITS_RELAT_BOND,
                                "Set limits (relative) on all bondlengths");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_PAR,ID_ZSCATTERER_MENU_PAR_LIMITS_RELAT_ANGLE,
                                "Set limits (relative) on all bond angles");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_PAR,ID_ZSCATTERER_MENU_PAR_LIMITS_RELAT_DIHED,
                                "Set limits (relative) on all dihedral angles");
      mpMenuBar->AddMenu("Atom",ID_ZSCATTERER_MENU_ATOM);
         mpMenuBar->AddMenuItem(ID_ZSCATTERER_MENU_ATOM,ID_ZSCATTERER_MENU_ATOM_ADD,
                                "Add an Atom");
         mpMenuBar->AddMenuItem(ID_ZSCATTERER_MENU_ATOM,ID_ZSCATTERER_MENU_ATOM_CHANGE_PIVOT,
                                "Change Pivot Atom");
   //Orientation
      wxBoxSizer* sizer=new wxBoxSizer(wxHORIZONTAL);
      mpScatterer->RefinableObj::Print();
#if 1
      WXFieldRefPar* pFieldPhi    =new WXFieldRefPar(this,"Phi:",
                                     &(mpZScatterer->GetPar(&(mpZScatterer->mPhi))) );
      WXFieldRefPar* pFieldChi    =new WXFieldRefPar(this,"Chi:",
                                     &(mpZScatterer->GetPar(&(mpZScatterer->mChi))) );
      WXFieldRefPar* pFieldPsi    =new WXFieldRefPar(this,"Psi:",
                                     &(mpZScatterer->GetPar(&(mpZScatterer->mPsi))) );
#else
      WXCrystObjBasic* pFieldPhi =mpZScatterer->GetPar(&(mpZScatterer->mPhi)).WXCreate(this);
      WXCrystObjBasic* pFieldChi =mpZScatterer->GetPar(&(mpZScatterer->mChi)).WXCreate(this);
      WXCrystObjBasic* pFieldPsi =mpZScatterer->GetPar(&(mpZScatterer->mPsi)).WXCreate(this);
#endif
      sizer->Add(pFieldPhi    ,0,wxALIGN_CENTER);
      sizer->Add(pFieldChi    ,0,wxALIGN_CENTER);
      sizer->Add(pFieldPsi    ,0,wxALIGN_CENTER);
      
      mpSizer->Add(sizer,0,wxALIGN_LEFT);
      mList.Add(pFieldPhi);
      mList.Add(pFieldChi);
      mList.Add(pFieldPsi);
   //Atoms
      mpWXZAtomRegistry=mpZScatterer->mZAtomRegistry.WXCreate(this);
      mpSizer->Add(mpWXZAtomRegistry,0,wxALIGN_LEFT);
      mList.Add(mpWXZAtomRegistry);
   
   this->CrystUpdate();
   this->Layout();
}

void WXZScatterer::OnMenuAddZAtom(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_ENTRY("WXZScatterer::OnMenuAddZAtom()",6)
   WXCrystValidateAllUserInput();
   
   int choice;
   int bondAtomId=0;
   double bondLength;
   int angleAtomId=0;;
   double angle;
   int dihedAtomId=0;;
   double dihed;
   //Scattering power
      const ScatteringPower *scattPow=WXDialogChooseFromRegistry(
                                 mpZScatterer->GetCrystal().GetScatteringPowerRegistry(),
                                 this,"Choose an atom type (ScatteringPower):",choice);
      if(0==scattPow)
      {
         VFN_DEBUG_EXIT("WXZScatterer::OnMenuAddZAtom():Cancelled",6)
         return;
      }
   //Bond atom 
   if(0<mpZScatterer->GetZAtomRegistry().GetNb())
   {

      const ZAtom *bondAtom=WXDialogChooseFromRegistry(mpZScatterer->GetZAtomRegistry(),
                                                 this,"Choose the bonded atom",bondAtomId);
      if(0==bondAtom)
      {
         VFN_DEBUG_EXIT("WXZScatterer::OnMenuAddZAtom():Cancelled",6)
         return;
      }
      //Bond length
      wxTextEntryDialog bondLengthDialog(this,"Enter bond length (Angstroems)",
                              "Bond length","1.5",wxOK | wxCANCEL);
      if(wxID_OK!=bondLengthDialog.ShowModal())
      {
         VFN_DEBUG_EXIT("WXZScatterer::OnMenuAddZAtom():Cancelled",6)
         return;
      }
      bondLengthDialog.GetValue().ToDouble(&bondLength);
   }
   //angle atom
   if(1<mpZScatterer->GetZAtomRegistry().GetNb())
   {
      const ZAtom *angleAtom=WXDialogChooseFromRegistry(mpZScatterer->GetZAtomRegistry(),
                                                 this,"Angle Atom",angleAtomId);
      if(0==angleAtom)
      {
         VFN_DEBUG_EXIT("WXZScatterer::OnMenuAddZAtom():Cancelled",6)
         return;
      }
      //angle
      wxTextEntryDialog angleDialog(this,"Enter bond angle (degrees)",
                              "Bond Angle","110",wxOK | wxCANCEL);
      if(wxID_OK!=angleDialog.ShowModal())
      {
         VFN_DEBUG_EXIT("WXZScatterer::OnMenuAddZAtom():Cancelled",6)
         return;
      }
      angleDialog.GetValue().ToDouble(&angle);
   }
   //dihedral atom
   if(2<mpZScatterer->GetZAtomRegistry().GetNb())
   {
      const ZAtom *dihedAtom=WXDialogChooseFromRegistry(mpZScatterer->GetZAtomRegistry(),
                                                 this,"Dihedral angle Atom",dihedAtomId);
      if(0==dihedAtom)
      {
         VFN_DEBUG_EXIT("WXZScatterer::OnMenuAddZAtom():Cancelled",6)
         return;
      }
   //dihedral angle
      wxTextEntryDialog dihedDialog(this,"Enter dihedral angle (degrees)",
                              "Dihedral Angle","0",wxOK | wxCANCEL);
      if(wxID_OK!=dihedDialog.ShowModal())
      {
         VFN_DEBUG_EXIT("WXZScatterer::OnMenuAddZAtom():Cancelled",6)
         return;
      }
      dihedDialog.GetValue().ToDouble(&dihed);
   }
   // Check all atoms are different
   if(1<mpZScatterer->GetZAtomRegistry().GetNb())
      if(bondAtomId==angleAtomId)
      {
         wxMessageDialog dumbUser(this,"Bond, angle an dihedral atoms *must* be different",
                                  "Whooops",wxOK|wxICON_EXCLAMATION);
         dumbUser.ShowModal();
         return;
      }
   if(2<mpZScatterer->GetZAtomRegistry().GetNb())
      if((bondAtomId==angleAtomId)||(bondAtomId==dihedAtomId)||(angleAtomId==dihedAtomId))
      {
         wxMessageDialog dumbUser(this,"Bond, angle an dihedral atoms *must* be different",
                                  "Whooops",wxOK|wxICON_EXCLAMATION);
         dumbUser.ShowModal();
         return;
      }
   char buf [5];
   sprintf(buf,"%d",mpZScatterer->GetNbComponent()+1);
   mpZScatterer->AddAtom (scattPow->GetName()+(string)buf,
                           scattPow,
                           bondAtomId,bondLength,
                           angleAtomId,angle*DEG2RAD,
                           dihedAtomId,dihed*DEG2RAD,
                           1.);
   this->CrystUpdate();
   VFN_DEBUG_EXIT("WXZScatterer::OnMenuAddZAtom()",6)
}

void WXZScatterer::OnMenuSetLimits(wxCommandEvent & event)
{//:TODO: Need to 
   VFN_DEBUG_ENTRY("WXZScatterer::OnMenuSetLimits()",6)
   WXCrystValidateAllUserInput();
   if(event.GetId()==ID_ZSCATTERER_MENU_PAR_LIMITS_RELAT_BOND)
   {
      double limit=.1;
      wxTextEntryDialog limitDialog(this,"Enter maximum shift in Angstroems\n The limits are taken symmetrically around current position:\n X0-l < X < X0+l ",
                              "Set limits (relative) for bondlengths",".1",wxOK | wxCANCEL);
      if(wxID_OK!=limitDialog.ShowModal())
      {
         VFN_DEBUG_EXIT("WXZScatterer::OnMenuSetLimits():Cancelled",6)
         return;
      }
      limitDialog.GetValue().ToDouble(&limit);
      if(limit<=0)
      {
         wxMessageDialog dumbUser(this,"Limit must be > 0 !",
                                  "Whooops",wxOK|wxICON_EXCLAMATION);
         dumbUser.ShowModal();
         return;
      }
      mpZScatterer->SetLimitsRelative(gpRefParTypeScattConformBondLength,-limit,limit);
   }
   if(event.GetId()==ID_ZSCATTERER_MENU_PAR_LIMITS_RELAT_ANGLE)
   {
      double limit=5;
      wxTextEntryDialog limitDialog(this,"Enter maximum shift in Degrees\n The limits are taken symmetrically around current position:\n X0-l < X < X0+l ",
                              "Set limits (relative) for bond angles",".1",wxOK | wxCANCEL);
      if(wxID_OK!=limitDialog.ShowModal())
      {
         VFN_DEBUG_EXIT("WXZScatterer::OnMenuSetLimits():Cancelled",6)
         return;
      }
      limitDialog.GetValue().ToDouble(&limit);
      if(limit<=0)
      {
         wxMessageDialog dumbUser(this,"Limit must be > 0 !",
                                  "Whooops",wxOK|wxICON_EXCLAMATION);
         dumbUser.ShowModal();
         return;
      }
      limit *=DEG2RAD;
      mpZScatterer->SetLimitsRelative(gpRefParTypeScattConformBondAngle,-limit,limit);
   }
   if(event.GetId()==ID_ZSCATTERER_MENU_PAR_LIMITS_RELAT_DIHED)
   {
      double limit=5;
      wxTextEntryDialog limitDialog(this,"Enter maximum shift in Degrees\n The limits are taken symmetrically around current position:\n X0-l < X < X0+l ",
                              "Set limits (relative) for dihedral angles",".1",wxOK | wxCANCEL);
      if(wxID_OK!=limitDialog.ShowModal())
      {
         VFN_DEBUG_EXIT("WXZScatterer::OnMenuSetLimits():Cancelled",6)
         return;
      }
      limitDialog.GetValue().ToDouble(&limit);
      if(limit<=0)
      {
         wxMessageDialog dumbUser(this,"Limit must be > 0 !",
                                  "Whooops",wxOK|wxICON_EXCLAMATION);
         dumbUser.ShowModal();
         return;
      }
      limit *=DEG2RAD;
      mpZScatterer->SetLimitsRelative(gpRefParTypeScattConformDihedAngle,-limit,limit);
   }
}

void WXZScatterer::OnMenuChangePivotAtom(wxCommandEvent &WXUNUSED(event))
{
   WXCrystValidateAllUserInput();
   int pivot=mpZScatterer->mCenterAtomIndex;
   const ZAtom *atom=WXDialogChooseFromRegistry(mpZScatterer->GetZAtomRegistry(),
                                              this,"Choose the new Pivot atom",pivot);
   if(0==atom)
   {
      VFN_DEBUG_EXIT("WXZScatterer::OnMenuAddZAtom():Cancelled",6)
      return;
   }
   mpZScatterer->mCenterAtomIndex=pivot;
   mpZScatterer->mClockScatterer.Click();
}

void WXZScatterer::OnMenuImportZMatrix(wxCommandEvent &WXUNUSED(event))
{
   wxFileDialog open(this,"Choose a file","","","*.fhz",
                                        wxOPEN | wxFILE_MUST_EXIST);
   if(open.ShowModal() != wxID_OK) return;
   ifstream fin (open.GetPath().c_str());
   if(!fin)
   {
      throw ObjCrystException("WXZScatterer::OnMenuImportZMatrix() : \
Error opening file for input:"+string(open.GetPath().c_str()));
   }
   mpZScatterer->ImportFenskeHallZMatrix(fin);
   fin.close();
}
void WXZScatterer::OnMenuExportZMatrix(wxCommandEvent &WXUNUSED(event))
{
   wxFileDialog save(this,"Choose a file","","","*.fhz",wxSAVE);
   if(save.ShowModal() != wxID_OK) return;
   ofstream fout (save.GetPath().c_str());
   if(!fout)
   {
      throw ObjCrystException("WXZScatterer::OnMenuExportZMatrix() : \
Error opening file for input:"+string(save.GetPath().c_str()));
   }
   mpZScatterer->ExportFenskeHallZMatrix(fout);
   fout.close();
}

void WXZScatterer::OnMenuConvert2Molecule(wxCommandEvent &WXUNUSED(event))
{
   Molecule *mol=ZScatterer2Molecule(mpZScatterer);
   mpZScatterer->GetCrystal().RemoveScatterer(mpZScatterer);
   mol->GetCrystal().AddScatterer(mol);
   mol->GetCrystal().UpdateDisplay();
}

}// namespace 












