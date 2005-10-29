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
   REAL x0=0,y0=0,z0=0;
   for(unsigned int i=0;i<nb;++i)
   {
      const REAL x=scatt->GetZAtomX(i);
      const REAL y=scatt->GetZAtomY(i);
      const REAL z=scatt->GetZAtomZ(i);
      x0+=x;
      y0+=y;
      z0+=z;
      mol->AddAtom(x,y,z,scatt->GetZAtomRegistry().GetObj(i).GetScatteringPower(),
                   scatt->GetZAtomRegistry().GetObj(i).GetName());
      
      #if 0
      if(i>0)
      {
         const RefinablePar* pLength=&(scatt->GetPar(&(scatt->GetZAtomRegistry()
                                                       .GetObj(i).GetZBondLength())));
         if(pLength->IsFixed())
            mol->AddBond(mol->GetAtom(i),mol->GetAtom(scatt->GetZBondAtom(i)),
                         pLength->GetValue(),.01,.02,false);
         else
            if(pLength->IsLimited())
               mol->AddBond(mol->GetAtom(i),mol->GetAtom(scatt->GetZBondAtom(i)),
                            (pLength->GetMin()+pLength->GetMax())/2.,.01,.02,false);
      }
      if(i>1)
      {
         const RefinablePar* pAngle=&(scatt->GetPar(&(scatt->GetZAtomRegistry()
                                                      .GetObj(i).GetZAngle())));
         if(pAngle->IsFixed())
            mol->AddBondAngle(mol->GetAtom(i),mol->GetAtom(scatt->GetZBondAtom(i)),
                              mol->GetAtom(scatt->GetZAngleAtom(i)),
                              pAngle->GetValue(),.01,.02,false);
         else
            if(pAngle->IsLimited())
               mol->AddBondAngle(mol->GetAtom(i),mol->GetAtom(scatt->GetZBondAtom(i)),
                                 mol->GetAtom(scatt->GetZAngleAtom(i)),
                                 (pAngle->GetMin()+pAngle->GetMax())/2.,.01,.02,false);
      }
      if(i>2)
      {
         const RefinablePar* pDihed=&(scatt->GetPar(&(scatt->GetZAtomRegistry()
                                                      .GetObj(i).GetZDihedralAngle())));
         MolAtom *p1=&(mol->GetAtom(i));
         MolAtom *p2=&(mol->GetAtom(scatt->GetZBondAtom(i)));
         MolAtom *p3=&(mol->GetAtom(scatt->GetZAngleAtom(i)));
         MolAtom *p4=&(mol->GetAtom(scatt->GetZDihedralAngleAtom(i)));
         if(  (fabs(GetBondAngle(*p1,*p2,*p3)-M_PI)>0.3)
            &&(fabs(GetBondAngle(*p1,*p2,*p4)-M_PI)>0.3)
            &&(fabs(GetBondAngle(*p1,*p3,*p4)-M_PI)>0.3)
            &&(fabs(GetBondAngle(*p2,*p3,*p4)-M_PI)>0.3))
         {
            if(pDihed->IsFixed())
               mol->AddDihedralAngle(*p1,*p2,*p3,*p4,pDihed->GetValue(),.01,.02,false);
            else
               if(((pDihed->GetMax()-pDihed->GetMax())<0.3)&&(i>2)&&(pDihed->IsLimited()))
                  mol->AddDihedralAngle(*p1,*p2,*p3,*p4,
                                        (pDihed->GetMin()+pDihed->GetMax())/2.,.01,.02,false);
         }
      }
      #endif
      mol->GetAtom(i).SetOccupancy(scatt->GetZAtomRegistry().GetObj(i).GetOccupancy());
   }

   CrystVector_REAL x(nb),y(nb),z(nb),radius(nb);
   vector<pair<const ScatteringPowerAtom *,long> > scattpow(nb);
   for(unsigned int i=0;i<nb;++i)
   {
      x(i)=mol->GetAtom(i).GetX();
      y(i)=mol->GetAtom(i).GetY();
      z(i)=mol->GetAtom(i).GetZ();
      if(mol->GetAtom(i).IsDummy())
      {
         radius(i)=-1;
         scattpow[i].first=0;
      }
      else
      {
         radius(i)=mol->GetAtom(i).GetScatteringPower().GetRadius();
         scattpow[i].first=dynamic_cast<const ScatteringPowerAtom *> 
                             (&(mol->GetAtom(i).GetScatteringPower()));
         scattpow[i].second=scattpow[i].first->GetAtomicNumber();
      }
   }
   for(unsigned int i=0;i<nb;++i)
   {
      if(scattpow[i].first==0) continue;
      const REAL x1=x(i),y1=y(i),z1=z(i);
      x += -x1;
      y += -y1;
      z += -z1;
      for(unsigned int j=i+1;j<nb;++j)
      {
         if(scattpow[j].first==0) continue;
         const REAL dist=sqrt(x(j)*x(j)+y(j)*y(j)+z(j)*z(j));
         //cout<<"          -> d="<<dist<<"("<<radius(i)<<","<<radius(j)<<"):"<<scattpow[i].second<<","<<scattpow[j].second<<endl;
         if(dist<(1.10*(radius(i)+radius(j))))
         {
            if((1!=scattpow[i].second)||(1!=scattpow[j].second))
            {
                  mol->AddBond(mol->GetAtom(i),mol->GetAtom(j),dist,.01,.02,false);
            }
         }
      }
      x += x1;
      y += y1;
      z += z1;
   }
   mol->BuildConnectivityTable();
   for(map<MolAtom*,set<MolAtom*> >::const_iterator pos=mol->GetConnectivityTable().begin();
       pos!=mol->GetConnectivityTable().end();++pos)
   {
      for(set<MolAtom*>::const_iterator pos1=pos->second.begin();
          pos1!=pos->second.end();++pos1)
      {
         for(set<MolAtom*>::const_iterator pos2=pos1;
          pos2!=pos->second.end();++pos2)
          {
            if(pos2==pos1) continue;
            if(mol->FindBondAngle(**pos1,*(pos->first),**pos2)== mol->GetBondAngleList().end())
               mol->AddBondAngle(**pos1,*(pos->first),**pos2,
                                 GetBondAngle(**pos1,*(pos->first),**pos2),0.01,0.02,false);
          }
      }
   }
   x0 /= nb;
   y0 /= nb;
   z0 /= nb;
   mol->GetCrystal().OrthonormalToFractionalCoords(x0,y0,z0);
   mol->SetX(x0);
   mol->SetY(y0);
   mol->SetZ(z0);
   mol->UpdateDisplay();
   VFN_DEBUG_EXIT("ZScatterer2Molecule()",6)
   return mol;
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

   if(0<mpZAtom->GetZScatterer().GetZAtomRegistry().Find(*mpZAtom))
   {
      RefinablePar *par=&(mpZAtom->GetZScatterer().GetPar(&(mpZAtom->mBondLength)));
      par->SetName("Bond(w/"+(mpZAtom->GetZScatterer().GetZAtomRegistry()
                     .GetObj(mpZAtom->GetZBondAtom())).GetName()+")");
      WXCrystObjBasic* pFieldBond=par->WXCreate(this);
      mpSizer->Add(pFieldBond,0,wxALIGN_LEFT);
      mList.Add(pFieldBond);
   }
   
   if(1<mpZAtom->GetZScatterer().GetZAtomRegistry().Find(*mpZAtom))
   {
      RefinablePar *par=&(mpZAtom->GetZScatterer().GetPar(&(mpZAtom->mAngle)));
      par->SetName("Angle(w/"+(mpZAtom->GetZScatterer().GetZAtomRegistry()
                     .GetObj(mpZAtom->GetZAngleAtom())).GetName()+")");
       WXCrystObjBasic* pFieldAngle
          =par->WXCreate(this);
      mpSizer->Add(pFieldAngle,0,wxALIGN_LEFT);
      mList.Add(pFieldAngle);
   }
   if(2<mpZAtom->GetZScatterer().GetZAtomRegistry().Find(*mpZAtom))
   {
      RefinablePar *par=&(mpZAtom->GetZScatterer().GetPar(&(mpZAtom->mDihed)));
      par->SetName("DihedralAngle(w/"+(mpZAtom->GetZScatterer().GetZAtomRegistry()
                     .GetObj(mpZAtom->GetZDihedralAngleAtom())).GetName()+")");
      WXCrystObjBasic* pFieldDihed=par->WXCreate(this);
      mpSizer->Add(pFieldDihed,0,wxALIGN_LEFT);
      mList.Add(pFieldDihed);
   }
   {
      RefinablePar *par=&(mpZAtom->GetZScatterer().GetPar(&(mpZAtom->mOccupancy)));
      par->SetName("Occup.");
      WXCrystObjBasic* pFieldOccup=par->WXCreate(this);
      mpSizer->Add(pFieldOccup,0,wxALIGN_LEFT);
      mList.Add(pFieldOccup);
   }

   this->SetSizer(mpSizer);
   
   this->BottomLayout(0);
   this->CrystUpdate(true);
   VFN_DEBUG_EXIT("WXZAtom::WXZAtom()",6)
}

void WXZAtom::CrystUpdate(const bool uui,const bool lock)
{
   VFN_DEBUG_ENTRY("WXZAtom::CrystUpdate()",6)
   if(lock) mMutex.Lock();
   mList.CrystUpdate(uui,false);
   if(lock) mMutex.Unlock();
   VFN_DEBUG_EXIT("WXZAtom::CrystUpdate()",6)
}

void WXZAtom::UpdateUI(const bool lock)
{
   VFN_DEBUG_ENTRY("WXZAtom::UpdateUI()",6)
   if(lock) mMutex.Lock();
   mList.UpdateUI(false);
   mpFieldName->SetValue(mpZAtom->GetName().c_str());
   if(0!=mpZAtom->GetScatteringPower())
      mpFieldScattPower->SetValue(mpZAtom->GetScatteringPower()->GetName());
   else
      mpFieldScattPower->SetValue("Dummy");
   if(lock) mMutex.Unlock();
   VFN_DEBUG_EXIT("WXZAtom::UpdateUI()",6)
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
   this->CrystUpdate(true);
   this->UpdateUI(true);
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
      mpSizer->SetItemMinSize(mpMenuBar,
                              mpMenuBar->GetSize().GetWidth(),
                              mpMenuBar->GetSize().GetHeight());
   //Orientation
      wxBoxSizer* sizer=new wxBoxSizer(wxHORIZONTAL);

      WXCrystObjBasic* pFieldPhi =mpZScatterer->GetPar(&(mpZScatterer->mPhi)).WXCreate(this);
      WXCrystObjBasic* pFieldChi =mpZScatterer->GetPar(&(mpZScatterer->mChi)).WXCreate(this);
      WXCrystObjBasic* pFieldPsi =mpZScatterer->GetPar(&(mpZScatterer->mPsi)).WXCreate(this);

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
   
   this->BottomLayout(0);
   this->CrystUpdate(true);
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
   this->CrystUpdate(true);
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
   //mol->GetCrystal().UpdateDisplay();
}

}// namespace 












