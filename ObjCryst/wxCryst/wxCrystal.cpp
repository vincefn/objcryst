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

#include <stdlib.h>

#include "wx/wx.h"
#include "wx/colordlg.h"
#include "wx/progdlg.h"

#include "wxCryst/wxCrystal.h"

#include "ObjCryst/Atom.h"
#include "ObjCryst/ZScatterer.h"
#include "ObjCryst/ScatteringPowerSphere.h"

extern "C" {
#include "GL/glu.h"
#include "wxCryst/trackball.h"
}

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

//#include "ObjCryst/Map.cpp"

namespace ObjCryst
{

////////////////////////////////////////////////////////////////////////
//
//    WXCrystal
//
////////////////////////////////////////////////////////////////////////
BEGIN_EVENT_TABLE(WXCrystal,wxEvtHandler)
   EVT_BUTTON(ID_WXOBJ_COLLAPSE,                      WXCrystObj::OnToggleCollapse)
   EVT_MENU(ID_REFOBJ_MENU_OBJ_SAVE,                  WXRefinableObj::OnMenuSave)
   EVT_MENU(ID_REFOBJ_MENU_OBJ_LOAD,                  WXRefinableObj::OnMenuLoad)
   EVT_MENU(ID_CRYSTAL_MENU_SAVECIF,                  WXCrystal::OnMenuSaveCIF)
   EVT_MENU(ID_CRYSTAL_MENU_SAVETEXT,                 WXCrystal::OnMenuSaveText)
   EVT_MENU(ID_REFOBJ_MENU_PAR_FIXALL,                WXRefinableObj::OnMenuFixAllPar)
   EVT_MENU(ID_REFOBJ_MENU_PAR_UNFIXALL,              WXRefinableObj::OnMenuUnFixAllPar)
   EVT_MENU(ID_REFOBJ_MENU_PAR_RANDOMIZE,             WXRefinableObj::OnMenuParRandomize)
   EVT_MENU(ID_CRYSTAL_MENU_PAR_ADDANTIBUMP,          WXCrystal::OnMenuAddAntiBumpDist)
   EVT_MENU(ID_CRYSTAL_MENU_PAR_SETRELATIVEXYZLIMITS, WXCrystal::OnMenuSetRelativeXYZLimits)
#ifdef OBJCRYST_GL
   EVT_MENU(ID_CRYSTAL_MENU_DISPLAY_3DVIEW,           WXCrystal::OnMenuCrystalGL)
#endif
   EVT_MENU(ID_CRYSTAL_MENU_SCATT_ADDSCATTPOWATOM,    WXCrystal::OnMenuAddScattPowAtom)
   EVT_MENU(ID_CRYSTAL_MENU_SCATT_ADDSCATTPOWSPHERE,  WXCrystal::OnMenuAddScattPowSphere)
   EVT_MENU(ID_CRYSTAL_MENU_SCATT_REMOVESCATTPOW,     WXCrystal::OnMenuRemoveScattPow)
   EVT_MENU(ID_CRYSTAL_MENU_SCATT_ADDATOM,            WXCrystal::OnMenuAddScatterer)
   EVT_MENU(ID_CRYSTAL_MENU_SCATT_ADDZSCATTERER,      WXCrystal::OnMenuAddScatterer)
   EVT_MENU(ID_CRYSTAL_MENU_SCATT_ADDTETRAHEDRON,     WXCrystal::OnMenuAddScatterer)
   EVT_MENU(ID_CRYSTAL_MENU_SCATT_ADDOCTAHEDRON,      WXCrystal::OnMenuAddScatterer)
   EVT_MENU(ID_CRYSTAL_MENU_SCATT_ADDTRIANGLE,        WXCrystal::OnMenuAddScatterer)
   EVT_MENU(ID_CRYSTAL_MENU_SCATT_ADDSQUAREPLANE,     WXCrystal::OnMenuAddScatterer)
   EVT_MENU(ID_CRYSTAL_MENU_SCATT_ADDCUBE,            WXCrystal::OnMenuAddScatterer)
   EVT_MENU(ID_CRYSTAL_MENU_SCATT_ADDANTIPRISMTETRAGONAL,WXCrystal::OnMenuAddScatterer)
   EVT_MENU(ID_CRYSTAL_MENU_SCATT_ADDPRISMTRIGONAL,   WXCrystal::OnMenuAddScatterer)
   EVT_MENU(ID_CRYSTAL_MENU_SCATT_ADDICOSAHEDRON,     WXCrystal::OnMenuAddScatterer)
   EVT_MENU(ID_CRYSTAL_MENU_SCATT_REMOVESCATTERER,    WXCrystal::OnMenuRemoveScatterer)
   EVT_MENU(ID_CRYSTAL_MENU_SCATT_DUPLICSCATTERER,    WXCrystal::OnMenuDuplicateScatterer)
   EVT_UPDATE_UI(ID_CRYST_UPDATEUI,                   WXRefinableObj::OnUpdateUI)
END_EVENT_TABLE()

WXCrystal::WXCrystal(wxWindow* parent, Crystal *obj):
WXRefinableObj(parent,(RefinableObj*)obj),mpCrystal(obj)
#ifdef OBJCRYST_GL
,mCrystalGLDisplayList(0),
/*mCrystalGLDisplayList(gGLDisplayListNb++),
mCrystalGLDisplayList(glGenLists(1)),*/
mCrystalGLDisplayListIsLocked(false),mpCrystalGL(0)
#endif
{
   VFN_DEBUG_MESSAGE("WXCrystal::WXCrystal()",6)
   //this->SetBackgroundColour("Red");
   //mpWXTitle->SetBackgroundColour(wxColour(255,200,200));
   mpWXTitle->SetForegroundColour(wxColour(255,0,0));
   // Menu
      mpMenuBar->AddMenu("File",ID_REFOBJ_MENU_OBJ);
         //mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_REFOBJ_MENU_OBJ_LOAD,"Load");
         //mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_REFOBJ_MENU_OBJ_SAVE,"Save");
         //mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_CRYSTAL_MENU_SAVECIF,"Save (CIF)");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_CRYSTAL_MENU_SAVETEXT,"Save as text");
      mpMenuBar->AddMenu("Parameters",ID_REFOBJ_MENU_PAR);
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_PAR,ID_CRYSTAL_MENU_PAR_ADDANTIBUMP,
                                "Add Antibump distance");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_PAR,ID_REFOBJ_MENU_PAR_FIXALL,"Fix all");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_PAR,ID_REFOBJ_MENU_PAR_UNFIXALL,"Unfix all");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_PAR,ID_REFOBJ_MENU_PAR_RANDOMIZE,
                                "Randomize Configuration");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_PAR,ID_CRYSTAL_MENU_PAR_SETRELATIVEXYZLIMITS,
                                "Set Relative Limits On All XYZ Parameters");
      mpMenuBar->AddMenu("Scatterers",ID_CRYSTAL_MENU_SCATT);
         mpMenuBar->AddMenuItem(ID_CRYSTAL_MENU_SCATT,ID_CRYSTAL_MENU_SCATT_ADDSCATTPOWATOM,
                                "Add Atomic Scattering Power");
         mpMenuBar->AddMenuItem(ID_CRYSTAL_MENU_SCATT,ID_CRYSTAL_MENU_SCATT_ADDSCATTPOWSPHERE,
                                "Add Sphere Scattering Power");
         mpMenuBar->AddMenuItem(ID_CRYSTAL_MENU_SCATT,ID_CRYSTAL_MENU_SCATT_REMOVESCATTPOW,
                                "Remove Scattering Power");
         //mpMenuBar->AppendSeparator();
         mpMenuBar->AddMenuItem(ID_CRYSTAL_MENU_SCATT,ID_CRYSTAL_MENU_SCATT_ADDATOM,
                                "Add Atom");
         mpMenuBar->AddMenuItem(ID_CRYSTAL_MENU_SCATT,ID_CRYSTAL_MENU_SCATT_ADDZSCATTERER,
                                "Add Z-Matrix Scatterer");
         mpMenuBar->AddMenuItem(ID_CRYSTAL_MENU_SCATT,ID_CRYSTAL_MENU_SCATT_ADDTETRAHEDRON,
                                "Add Tetrahedron");
         mpMenuBar->AddMenuItem(ID_CRYSTAL_MENU_SCATT,ID_CRYSTAL_MENU_SCATT_ADDOCTAHEDRON,
                                "Add Octahedron");
         mpMenuBar->AddMenuItem(ID_CRYSTAL_MENU_SCATT,ID_CRYSTAL_MENU_SCATT_ADDTRIANGLE,
                                "Add Triangle Plane");
         mpMenuBar->AddMenuItem(ID_CRYSTAL_MENU_SCATT,ID_CRYSTAL_MENU_SCATT_ADDSQUAREPLANE,
                                "Add Square Plane");
         mpMenuBar->AddMenuItem(ID_CRYSTAL_MENU_SCATT,ID_CRYSTAL_MENU_SCATT_ADDCUBE,
                                "Add Cube");
         mpMenuBar->AddMenuItem(ID_CRYSTAL_MENU_SCATT,
                                ID_CRYSTAL_MENU_SCATT_ADDANTIPRISMTETRAGONAL,
                                "Add Antiprism Tetragonal");
         mpMenuBar->AddMenuItem(ID_CRYSTAL_MENU_SCATT,ID_CRYSTAL_MENU_SCATT_ADDPRISMTRIGONAL,
                                "Add Prism Trigonal");
         mpMenuBar->AddMenuItem(ID_CRYSTAL_MENU_SCATT,ID_CRYSTAL_MENU_SCATT_ADDICOSAHEDRON,
                                "Add Icosahedron");
         //mpMenuBar->AppendSeparator();
         mpMenuBar->AddMenuItem(ID_CRYSTAL_MENU_SCATT,ID_CRYSTAL_MENU_SCATT_REMOVESCATTERER,
                                "Remove Scatterer");
         mpMenuBar->AddMenuItem(ID_CRYSTAL_MENU_SCATT,ID_CRYSTAL_MENU_SCATT_DUPLICSCATTERER,
                                "Duplicate Scatterer");
      mpMenuBar->AddMenu("Display",ID_CRYSTAL_MENU_DISPLAY);
         mpMenuBar->AddMenuItem(ID_CRYSTAL_MENU_DISPLAY,ID_CRYSTAL_MENU_DISPLAY_3DVIEW,
                                "3D Display");
   // Lattice
      wxBoxSizer* lattice=new wxBoxSizer(wxHORIZONTAL);
#if 1
      WXFieldRefPar* pFieldLatticeA    =new WXFieldRefPar(this,"a:",
                                     &(mpCrystal->GetPar(mpCrystal->mCellDim.data()+0)) );
          
      WXFieldRefPar* pFieldLatticeB    =new WXFieldRefPar(this,"b:",
                                     &(mpCrystal->GetPar(mpCrystal->mCellDim.data()+1)) );
          
      WXFieldRefPar* pFieldLatticeC    =new WXFieldRefPar(this,"c:",
                                     &(mpCrystal->GetPar(mpCrystal->mCellDim.data()+2)) );
          
      WXFieldRefPar* pFieldLatticeAlpha=new WXFieldRefPar(this,"alpha:",
                                     &(mpCrystal->GetPar(mpCrystal->mCellDim.data()+3)) );
          
      WXFieldRefPar* pFieldLatticeBeta =new WXFieldRefPar(this,"beta:",
                                     &(mpCrystal->GetPar(mpCrystal->mCellDim.data()+4)) );
          
      WXFieldRefPar* pFieldLatticeGamma=new WXFieldRefPar(this,"gamma:",
                                     &(mpCrystal->GetPar(mpCrystal->mCellDim.data()+5)) );
#else
      WXCrystObjBasic* pFieldLatticeA 
         =mpCrystal->GetPar(mpCrystal->mCellDim.data()+0).WXCreate(this);
      WXCrystObjBasic* pFieldLatticeB 
         =mpCrystal->GetPar(mpCrystal->mCellDim.data()+1).WXCreate(this);
      WXCrystObjBasic* pFieldLatticeC 
         =mpCrystal->GetPar(mpCrystal->mCellDim.data()+2).WXCreate(this);
      WXCrystObjBasic* pFieldLatticeAlpha 
         =mpCrystal->GetPar(mpCrystal->mCellDim.data()+3).WXCreate(this);
      WXCrystObjBasic* pFieldLatticeBeta 
         =mpCrystal->GetPar(mpCrystal->mCellDim.data()+4).WXCreate(this);
      WXCrystObjBasic* pFieldLatticeGamma
         =mpCrystal->GetPar(mpCrystal->mCellDim.data()+5).WXCreate(this);
#endif
      lattice->Add(pFieldLatticeA    ,0,wxALIGN_CENTER);
      lattice->Add(pFieldLatticeB    ,0,wxALIGN_CENTER);
      lattice->Add(pFieldLatticeC    ,0,wxALIGN_CENTER);
      lattice->Add(pFieldLatticeAlpha,0,wxALIGN_CENTER);
      lattice->Add(pFieldLatticeBeta ,0,wxALIGN_CENTER);
      lattice->Add(pFieldLatticeGamma,0,wxALIGN_CENTER);
      lattice->Layout();
      
      mpSizer->Add(lattice,0,wxALIGN_LEFT);
      mList.Add(pFieldLatticeA);
      mList.Add(pFieldLatticeB);
      mList.Add(pFieldLatticeC);
      mList.Add(pFieldLatticeAlpha);
      mList.Add(pFieldLatticeBeta);
      mList.Add(pFieldLatticeGamma);
      
   // SpaceGroup
      mpFieldSpacegroup=new WXFieldName(this,"SpaceGroup:",this,ID_CRYSTAL_SPACEGROUP);
      mpSizer->Add(mpFieldSpacegroup,0,wxALIGN_LEFT);
      mList.Add(mpFieldSpacegroup);
      
   // Scattering Powers
      mpWXScatteringPowerRegistry=mpCrystal
                    ->GetScatteringPowerRegistry().WXCreate(this);
      mpSizer->Add(mpWXScatteringPowerRegistry,0,wxALIGN_LEFT);
      mList.Add(mpWXScatteringPowerRegistry);
      
   // Scatterers
      mpWXScattererRegistry=mpCrystal
                    ->GetScattererRegistry().WXCreate(this);
      mpSizer->Add(mpWXScattererRegistry,0,wxALIGN_LEFT);
      mList.Add(mpWXScattererRegistry);
      
   this->CrystUpdate();
   this->Layout();
   VFN_DEBUG_MESSAGE("WXCrystal::WXCrystal():End",6)
}

void WXCrystal::CrystUpdate()
{
   VFN_DEBUG_ENTRY("WXCrystal::CrystUpdate()",7)
   this->WXRefinableObj::CrystUpdate();
   //mWXParent->Layout();
   #ifdef OBJCRYST_GL
   this->UpdateGL();
   #endif
   VFN_DEBUG_EXIT("WXCrystal::CrystUpdate():End",7)
}

#ifdef OBJCRYST_GL
void WXCrystal::UpdateGL(const bool onlyIndependentAtoms,
                         const REAL xMin,const REAL xMax,
                         const REAL yMin,const REAL yMax,
                         const REAL zMin,const REAL zMax)
{
   VFN_DEBUG_ENTRY("WXCrystal::UpdateGL()",8)
   WXCrystValidateAllUserInput();
   if(mpCrystalGL!=0)
   {
      VFN_DEBUG_MESSAGE("WXCrystal::UpdateGL():mpCrystalGL",7)
      this->GrabCrystalGLDisplayList();
      if(mCrystalGLDisplayList==0) mCrystalGLDisplayList=glGenLists(1);
      
      // During a refinement (multi-threaded)
      // Wait until the display list has been updated by the main thread...
      static bool cont;//:TODO: not static, but mutable member function (if >1 crystal,...)
      if(false==wxThread::IsMain())
      {
         this->ReleaseCrystalGLDisplayList();
         cont=false;
         wxCommandEvent event(wxEVT_COMMAND_MENU_SELECTED,ID_GLCRYSTAL_MENU_UPDATE);
         wxPostEvent(mpCrystalGL,event);
         while(!cont) wxUsleep(10);
         VFN_DEBUG_EXIT("WXCrystal::UpdateGL()-Not in main thread :End",8)
         return;
      }
      mpCrystalGL->SetCurrent();
      glNewList(mCrystalGLDisplayList,GL_COMPILE);
         glPushMatrix();
            mpCrystal->GLInitDisplayList(onlyIndependentAtoms,xMin,xMax,yMin,yMax,zMin,zMax);
            //ScatteringPowerMap map1(mpCrystal->GetScatteringPowerRegistry().GetObj(0),
            //                            *mpCrystal,.02,.05,.05,RAD_XRAY);
            //map1.GLInitDisplayList(xMin,xMax,yMin,yMax,zMin,zMax);
            //UnitCellScattererDensityMap map2(*mpCrystal,21,21,21);
            //cout << map2.GetMap3D()<<endl;
            //map2.GLInitDisplayList(xMin,xMax,yMin,yMax,zMin,zMax);
         glPopMatrix();
      glEndList();
      //#ifdef __WINDOWS__
      cont=true;
      //#endif
      this->ReleaseCrystalGLDisplayList();
      mpCrystalGL->CrystUpdate();
   }
   else
   {
      VFN_DEBUG_MESSAGE("WXCrystal::UpdateGL():No mpCrystalGL",7)
   }
   VFN_DEBUG_EXIT("WXCrystal::UpdateGL():End",8)
}

int WXCrystal::GrabCrystalGLDisplayList()const
{
   VFN_DEBUG_MESSAGE("WXCrystal::GrabCrystalGLDisplayList()",7)
   //:KLUDGE: ? or OK ?
   while(mCrystalGLDisplayListIsLocked) wxUsleep(5);
   mCrystalGLDisplayListIsLocked=true;
   VFN_DEBUG_MESSAGE("WXCrystal::GrabCrystalGLDisplayList():"<<mCrystalGLDisplayList,7)
   return mCrystalGLDisplayList;
}
void WXCrystal::ReleaseCrystalGLDisplayList()const
{
   VFN_DEBUG_MESSAGE("WXCrystal::ReleaseCrystalGLDisplayList()",7)
   mCrystalGLDisplayListIsLocked=false;
}
bool WXCrystal::GLDisplayListIsLocked()const
{
   return mCrystalGLDisplayListIsLocked;
}

void WXCrystal::OnMenuCrystalGL(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXCrystal::OnMenuCrystalGL()",6)
   if(mpCrystalGL!=0) return;
   wxFrame* frame= new wxFrame(this,-1,mpCrystal->GetName().c_str(),
                               wxDefaultPosition,wxSize(400,400));
   mpCrystalGL=new WXGLCrystalCanvas(this,frame,-1);
   frame->Show(true);
   this->UpdateGL();
}
void WXCrystal::NotifyCrystalGLDelete()
{
   VFN_DEBUG_MESSAGE("WXCrystal::NotifyCrystalGLDelete()",7)
   mpCrystalGL=0;
}
WXGLCrystalCanvas * WXCrystal::GetCrystalGL()
{
   VFN_DEBUG_MESSAGE("WXCrystal::GetCrystalGL()",7)
   return mpCrystalGL;
}
#endif

void WXCrystal::OnMenuSaveCIF(wxCommandEvent & WXUNUSED(event))
{
   WXCrystValidateAllUserInput();
   wxFileDialog save(this,"Choose a file","","","*.cif",wxSAVE | wxOVERWRITE_PROMPT);
   if(save.ShowModal() != wxID_OK) return;
   
   ofstream out(save.GetPath().c_str());
   if(!out) return;//:TODO:
   mpCrystal->CIFOutput(out);
   out.close();
}

void WXCrystal::OnMenuSaveText(wxCommandEvent & WXUNUSED(event))
{
   WXCrystValidateAllUserInput();
   wxFileDialog save(this,"Choose a file","","","*.txt",wxSAVE | wxOVERWRITE_PROMPT);
   if(save.ShowModal() != wxID_OK) return;
   
   ofstream out(save.GetPath().c_str());
   if(!out) return;//:TODO:
   mpCrystal->Print(out);
   mpCrystal->PrintMinDistanceTable(.05,out);
   out.close();
}

void WXCrystal::OnMenuAddScattPowAtom(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXCrystal::OnMenuAddScattPowAtom()",6)
   ScatteringPowerAtom *scatt=new ScatteringPowerAtom("Change me","H");
   mpCrystal->AddScatteringPower(scatt);
   VFN_DEBUG_MESSAGE("WXCrystal::OnMenuAddScattPowAtom():End",6)
   this->Layout();
}

void WXCrystal::OnMenuAddScattPowSphere(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_ENTRY("WXCrystal::OnMenuAddScattSphere()",6)
   ScatteringPower *scatt= new ScatteringPowerSphere;
   mpCrystal->AddScatteringPower(scatt);
   this->Layout();
   VFN_DEBUG_EXIT("WXCrystal::OnMenuAddScattPowSphere()",6)
}

void WXCrystal::OnMenuRemoveScattPow(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_ENTRY("WXCrystal::OnButtonRemoveScattPow()",6)
   WXCrystValidateAllUserInput();
   int choice;
   ScatteringPower *scatt=
      WXDialogChooseFromRegistry(mpCrystal->GetScatteringPowerRegistry(),this,
                                 "Choose Scattering Power to remove:",choice);
   if(0==scatt)
   {
      VFN_DEBUG_EXIT("WXCrystal::OnButtonRemoveScattPow():Cancelled",6)
      return;
   }
   mpCrystal->RemoveScatteringPower(scatt);
   mpCrystal->XMLOutput(cout);
   VFN_DEBUG_EXIT("WXCrystal::OnButtonRemoveScattPow()",6)
   this->Layout();
}

void WXCrystal::OnMenuAddScatterer(wxCommandEvent &event)
{
   VFN_DEBUG_ENTRY("WXCrystal::OnMenuAddScatterer()",6)
   WXCrystValidateAllUserInput();
   Scatterer *scatt;
   switch(event.GetId())
   {
      case ID_CRYSTAL_MENU_SCATT_ADDATOM:
      {
         int choice;
         ScatteringPower *scattPow=
            WXDialogChooseFromRegistry(mpCrystal->GetScatteringPowerRegistry(),this,
                                       "Choose an atom type (ScatteringPower):",choice);
         if(0==scattPow)
         {
            VFN_DEBUG_EXIT("WXCrystal::OnMenuAddScatterer():Canceled",6)
            return;
         }
         scatt=new Atom(0,0,0,"Change Me!",scattPow);
         break;
      }
      case ID_CRYSTAL_MENU_SCATT_ADDZSCATTERER:
      {
         scatt=new ZScatterer("Change Me!",*mpCrystal);
         break;
      }
      case ID_CRYSTAL_MENU_SCATT_ADDTETRAHEDRON:
      {
         int choice;
         //Scattering power 1
            const ScatteringPower *scattPow1=WXDialogChooseFromRegistry(
                                       mpCrystal->GetScatteringPowerRegistry(),
                                       this,"Central atom type (ScatteringPower):",choice);
            if(0==scattPow1)
            {
               VFN_DEBUG_EXIT("WXCrystal::OnMenuAddScatterer():Canceled",6)
               return;
            }
         //Scattering power 2
            const ScatteringPower *scattPow2=WXDialogChooseFromRegistry(
                                       mpCrystal->GetScatteringPowerRegistry(),
                                       this,"Corner atom type (ScatteringPower):",choice);
            if(0==scattPow2)
            {
               VFN_DEBUG_EXIT("WXCrystal::OnMenuAddScatterer():Canceled",6)
               return;
            }
         //Bond length
            wxTextEntryDialog bondLengthDialog(this,"Bond length",
                                    "Enter bond length (Angstroems)","1",wxOK | wxCANCEL);
            if(wxID_OK!=bondLengthDialog.ShowModal())
            {
               VFN_DEBUG_EXIT("WXZScatterer::OnMenuAddZAtom():Cancelled",6)
               return;
            }
            double bondLength;
            bondLengthDialog.GetValue().ToDouble(&bondLength);
            
         scatt=new ZPolyhedron(TETRAHEDRON,*mpCrystal,0,0,0,"Change Me!",
                               scattPow1,scattPow2,bondLength);
         break;
      }
      case ID_CRYSTAL_MENU_SCATT_ADDOCTAHEDRON:
      {
         int choice;
         //Scattering power 1
            const ScatteringPower *scattPow1=WXDialogChooseFromRegistry(
                                       mpCrystal->GetScatteringPowerRegistry(),
                                       this,"Central atom type (ScatteringPower):",choice);
            if(0==scattPow1)
            {
               VFN_DEBUG_EXIT("WXCrystal::OnMenuAddScatterer():Canceled",6)
               return;
            }
         //Scattering power 2
            const ScatteringPower *scattPow2=WXDialogChooseFromRegistry(
                                       mpCrystal->GetScatteringPowerRegistry(),
                                       this,"Corner atom type (ScatteringPower):",choice);
            if(0==scattPow2)
            {
               VFN_DEBUG_EXIT("WXCrystal::OnMenuAddScatterer():Canceled",6)
               return;
            }
         //Bond length
            wxTextEntryDialog bondLengthDialog(this,"Bond length",
                                    "Enter bond length (Angstroems)","1",wxOK | wxCANCEL);
            if(wxID_OK!=bondLengthDialog.ShowModal())
            {
               VFN_DEBUG_EXIT("WXZScatterer::OnMenuAddZAtom():Cancelled",6)
               return;
            }
            double bondLength;
            bondLengthDialog.GetValue().ToDouble(&bondLength);
            
         scatt=new ZPolyhedron(OCTAHEDRON,*mpCrystal,0,0,0,"Change Me!",
                               scattPow1,scattPow2,bondLength);
         break;
      }
      case ID_CRYSTAL_MENU_SCATT_ADDTRIANGLE:
      {
         VFN_DEBUG_MESSAGE("WXCrystal::OnMenuAddScatterer():Add triangle plane",6)
         int choice;
         //Scattering power 1
            const ScatteringPower *scattPow1=WXDialogChooseFromRegistry(
                                       mpCrystal->GetScatteringPowerRegistry(),
                                       this,"Central atom type (ScatteringPower):",choice);
            if(0==scattPow1)
            {
               VFN_DEBUG_EXIT("WXCrystal::OnMenuAddScatterer():Canceled",6)
               return;
            }
         //Scattering power 2
            const ScatteringPower *scattPow2=WXDialogChooseFromRegistry(
                                       mpCrystal->GetScatteringPowerRegistry(),
                                       this,"Corner atom type (ScatteringPower):",choice);
            if(0==scattPow2)
            {
               VFN_DEBUG_EXIT("WXCrystal::OnMenuAddScatterer():Canceled",6)
               return;
            }
         //Bond length
            wxTextEntryDialog bondLengthDialog(this,"Bond length",
                                    "Enter bond length (Angstroems)","1",wxOK | wxCANCEL);
            if(wxID_OK!=bondLengthDialog.ShowModal())
            {
               VFN_DEBUG_EXIT("WXZScatterer::OnMenuAddZAtom():Cancelled",6)
               return;
            }
            double bondLength;
            bondLengthDialog.GetValue().ToDouble(&bondLength);
            
         scatt=new ZPolyhedron(TRIANGLE_PLANE,*mpCrystal,0,0,0,"Change Me!",
                               scattPow1,scattPow2,bondLength);
         break;
      }
      case ID_CRYSTAL_MENU_SCATT_ADDSQUAREPLANE:
      {
         VFN_DEBUG_MESSAGE("WXCrystal::OnMenuAddScatterer():Add square plane",6)
         int choice;
         //Scattering power 1
            const ScatteringPower *scattPow1=WXDialogChooseFromRegistry(
                                       mpCrystal->GetScatteringPowerRegistry(),
                                       this,"Central atom type (ScatteringPower):",choice);
            if(0==scattPow1)
            {
               VFN_DEBUG_EXIT("WXCrystal::OnMenuAddScatterer():Canceled",6)
               return;
            }
         //Scattering power 2
            const ScatteringPower *scattPow2=WXDialogChooseFromRegistry(
                                       mpCrystal->GetScatteringPowerRegistry(),
                                       this,"Corner atom type (ScatteringPower):",choice);
            if(0==scattPow2)
            {
               VFN_DEBUG_EXIT("WXCrystal::OnMenuAddScatterer():Canceled",6)
               return;
            }
         //Bond length
            wxTextEntryDialog bondLengthDialog(this,"Bond length",
                                    "Enter bond length (Angstroems)","1",wxOK | wxCANCEL);
            if(wxID_OK!=bondLengthDialog.ShowModal())
            {
               VFN_DEBUG_EXIT("WXZScatterer::OnMenuAddZAtom():Cancelled",6)
               return;
            }
            double bondLength;
            bondLengthDialog.GetValue().ToDouble(&bondLength);
            
         scatt=new ZPolyhedron(SQUARE_PLANE,*mpCrystal,0,0,0,"Change Me!",
                               scattPow1,scattPow2,bondLength);
         break;
      }
      case ID_CRYSTAL_MENU_SCATT_ADDCUBE:
      {
         int choice;
         //Scattering power 1
            const ScatteringPower *scattPow1=WXDialogChooseFromRegistry(
                                       mpCrystal->GetScatteringPowerRegistry(),
                                       this,"Central atom type (ScatteringPower):",choice);
            if(0==scattPow1)
            {
               VFN_DEBUG_EXIT("WXCrystal::OnMenuAddScatterer():Canceled",6)
               return;
            }
         //Scattering power 2
            const ScatteringPower *scattPow2=WXDialogChooseFromRegistry(
                                       mpCrystal->GetScatteringPowerRegistry(),
                                       this,"Corner atom type (ScatteringPower):",choice);
            if(0==scattPow2)
            {
               VFN_DEBUG_EXIT("WXCrystal::OnMenuAddScatterer():Canceled",6)
               return;
            }
         //Bond length
            wxTextEntryDialog bondLengthDialog(this,"Bond length",
                                    "Enter bond length (Angstroems)","1",wxOK | wxCANCEL);
            if(wxID_OK!=bondLengthDialog.ShowModal())
            {
               VFN_DEBUG_EXIT("WXZScatterer::OnMenuAddZAtom():Cancelled",6)
               return;
            }
            double bondLength;
            bondLengthDialog.GetValue().ToDouble(&bondLength);
            
         scatt=new ZPolyhedron(CUBE,*mpCrystal,0,0,0,"Change Me!",
                               scattPow1,scattPow2,bondLength);
         break;
      }
      case ID_CRYSTAL_MENU_SCATT_ADDANTIPRISMTETRAGONAL:
      {
         int choice;
         //Scattering power 1
            const ScatteringPower *scattPow1=WXDialogChooseFromRegistry(
                                       mpCrystal->GetScatteringPowerRegistry(),
                                       this,"Central atom type (ScatteringPower):",choice);
            if(0==scattPow1)
            {
               VFN_DEBUG_EXIT("WXCrystal::OnMenuAddScatterer():Canceled",6)
               return;
            }
         //Scattering power 2
            const ScatteringPower *scattPow2=WXDialogChooseFromRegistry(
                                       mpCrystal->GetScatteringPowerRegistry(),
                                       this,"Corner atom type (ScatteringPower):",choice);
            if(0==scattPow2)
            {
               VFN_DEBUG_EXIT("WXCrystal::OnMenuAddScatterer():Canceled",6)
               return;
            }
         //Bond length
            wxTextEntryDialog bondLengthDialog(this,"Bond length",
                                    "Enter bond length (Angstroems)","1",wxOK | wxCANCEL);
            if(wxID_OK!=bondLengthDialog.ShowModal())
            {
               VFN_DEBUG_EXIT("WXZScatterer::OnMenuAddZAtom():Cancelled",6)
               return;
            }
            double bondLength;
            bondLengthDialog.GetValue().ToDouble(&bondLength);
            
         scatt=new ZPolyhedron(ANTIPRISM_TETRAGONAL,*mpCrystal,0,0,0,"Change Me!",
                               scattPow1,scattPow2,bondLength);
         break;
      }
      case ID_CRYSTAL_MENU_SCATT_ADDPRISMTRIGONAL:
      {
         int choice;
         //Scattering power 1
            const ScatteringPower *scattPow1=WXDialogChooseFromRegistry(
                                       mpCrystal->GetScatteringPowerRegistry(),
                                       this,"Central atom type (ScatteringPower):",choice);
            if(0==scattPow1)
            {
               VFN_DEBUG_EXIT("WXCrystal::OnMenuAddScatterer():Canceled",6)
               return;
            }
         //Scattering power 2
            const ScatteringPower *scattPow2=WXDialogChooseFromRegistry(
                                       mpCrystal->GetScatteringPowerRegistry(),
                                       this,"Corner atom type (ScatteringPower):",choice);
            if(0==scattPow2)
            {
               VFN_DEBUG_EXIT("WXCrystal::OnMenuAddScatterer():Canceled",6)
               return;
            }
         //Bond length
            wxTextEntryDialog bondLengthDialog(this,"Bond length",
                                    "Enter bond length (Angstroems)","1",wxOK | wxCANCEL);
            if(wxID_OK!=bondLengthDialog.ShowModal())
            {
               VFN_DEBUG_EXIT("WXZScatterer::OnMenuAddZAtom():Cancelled",6)
               return;
            }
            double bondLength;
            bondLengthDialog.GetValue().ToDouble(&bondLength);
            
         scatt=new ZPolyhedron(PRISM_TRIGONAL,*mpCrystal,0,0,0,"Change Me!",
                               scattPow1,scattPow2,bondLength);
         break;
      }
      case ID_CRYSTAL_MENU_SCATT_ADDICOSAHEDRON:
      {
         int choice;
         //Scattering power 1
            const ScatteringPower *scattPow1=WXDialogChooseFromRegistry(
                                       mpCrystal->GetScatteringPowerRegistry(),
                                       this,"Central atom type (ScatteringPower):",choice);
            if(0==scattPow1)
            {
               VFN_DEBUG_EXIT("WXCrystal::OnMenuAddScatterer():Canceled",6)
               return;
            }
         //Scattering power 2
            const ScatteringPower *scattPow2=WXDialogChooseFromRegistry(
                                       mpCrystal->GetScatteringPowerRegistry(),
                                       this,"Corner atom type (ScatteringPower):",choice);
            if(0==scattPow2)
            {
               VFN_DEBUG_EXIT("WXCrystal::OnMenuAddScatterer():Canceled",6)
               return;
            }
         //Bond length
            wxTextEntryDialog bondLengthDialog(this,"Bond length",
                                    "Enter bond length (Angstroems)","1",wxOK | wxCANCEL);
            if(wxID_OK!=bondLengthDialog.ShowModal())
            {
               VFN_DEBUG_EXIT("WXZScatterer::OnMenuAddZAtom():Cancelled",6)
               return;
            }
            double bondLength;
            bondLengthDialog.GetValue().ToDouble(&bondLength);
            
         scatt=new ZPolyhedron(ICOSAHEDRON,*mpCrystal,0,0,0,"Change Me!",
                               scattPow1,scattPow2,bondLength);
         break;
      }
   }
   mpCrystal->AddScatterer(scatt);
   //mpCrystal->XMLOutput(cout);
   VFN_DEBUG_MESSAGE("WXCrystal::OnMenuAddScatterer():calling Layout()",6)
   this->CrystUpdate();
   //this->Layout();
   VFN_DEBUG_EXIT("WXCrystal::OnMenuAddScatterer()",6)
}

void WXCrystal::OnMenuRemoveScatterer(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXCrystal::OnButtonRemoveScatterer()",6)
   WXCrystValidateAllUserInput();
   int choice;
   Scatterer *scatt=WXDialogChooseFromRegistry(mpCrystal->GetScattererRegistry(),this,
                                             "Select the Scatterer to remove:",choice);
   if(0==scatt) return;
   mpCrystal->RemoveScatterer(scatt);
   VFN_DEBUG_MESSAGE("WXCrystal::OnButtonRemoveScatterer():End",6)
   this->Layout();
}

void WXCrystal::OnMenuDuplicateScatterer(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_ENTRY("WXCrystal::OnMenuDuplicateScatterer()",6)
   WXCrystValidateAllUserInput();
   int choice;
   Scatterer *scatt=WXDialogChooseFromRegistry(mpCrystal->GetScattererRegistry(),this,
                                             "Select the Scatterer to duplicate:",choice);
   if(0==scatt) return;
   Scatterer *copy=scatt->CreateCopy();
   scatt->SetName(scatt->GetName()+(string)"(copy)");
   mpCrystal->AddScatterer(copy);
   this->Layout();
   VFN_DEBUG_EXIT("WXCrystal::OnMenuDuplicateScatterer():End",6)
}

void WXCrystal::OnMenuAddAntiBumpDist(wxCommandEvent & WXUNUSED(event))
{
   WXCrystValidateAllUserInput();
   int choice;
   //Scattering power 1
      const ScatteringPower *scattPow1=WXDialogChooseFromRegistry(
                                 mpCrystal->GetScatteringPowerRegistry(),
                                 this,"Atom type (ScatteringPower) #1",choice);
      if(0==scattPow1)
      {
         VFN_DEBUG_EXIT("WXCrystal::OnMenuAddAntiBumpDist():Canceled",6)
         return;
      }
   //Scattering power 2
      const ScatteringPower *scattPow2=WXDialogChooseFromRegistry(
                                 mpCrystal->GetScatteringPowerRegistry(),
                                 this,"Atom type (ScatteringPower) #2",choice);
      if(0==scattPow2)
      {
         VFN_DEBUG_EXIT("WXCrystal::OnMenuAddAntiBumpDist():Canceled",6)
         return;
      }
   //Distance
      wxTextEntryDialog bondLengthDialog(this,"Antibump Distance",
                              "Enter antibmup distance (<bond length) (Angstroems)",
                              "1.2",wxOK | wxCANCEL);
      if(wxID_OK!=bondLengthDialog.ShowModal())
      {
         VFN_DEBUG_EXIT("WXCrystal::OnMenuAddAntiBumpDist():Cancelled",6)
         return;
      }
      double bondLength;
      bondLengthDialog.GetValue().ToDouble(&bondLength);
      
   mpCrystal->SetBumpMergeDistance(*scattPow1,*scattPow2,bondLength);
}

void WXCrystal::OnMenuSetRelativeXYZLimits(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_ENTRY("WXCrystal::OnMenuSetRelativeXYZLimits():Cancelled",6)
   WXCrystValidateAllUserInput();
   wxTextEntryDialog limitDialog(this,"Relative limits",
                           "Enter relative limits for x,y,z (Angstroems)",
                           "0.5",wxOK | wxCANCEL);
   if(wxID_OK!=limitDialog.ShowModal())
   {
      VFN_DEBUG_EXIT("WXCrystal::OnMenuSetRelativeXYZLimits():Cancelled",6)
      return;
   }
   double limit;
   limitDialog.GetValue().ToDouble(&limit);
   limit=fabs(limit);
 
   mpCrystal->SetLimitsRelative(gpRefParTypeScattTranslX,
                                -limit/mpCrystal->GetLatticePar(0),
                                limit/mpCrystal->GetLatticePar(0));
   mpCrystal->SetLimitsRelative(gpRefParTypeScattTranslY,
                                -limit/mpCrystal->GetLatticePar(1),
                                limit/mpCrystal->GetLatticePar(1));
   mpCrystal->SetLimitsRelative(gpRefParTypeScattTranslZ,
                                -limit/mpCrystal->GetLatticePar(2),
                                limit/mpCrystal->GetLatticePar(2));
   VFN_DEBUG_EXIT("WXCrystal::OnMenuSetRelativeXYZLimits()",6)
}

bool WXCrystal::OnChangeName(const int id)
{
   VFN_DEBUG_MESSAGE("WXCrystal::OnChangeName()",6)
   if(this->WXRefinableObj::OnChangeName(id)==true) return true;
   switch(id)
   {
      case ID_CRYSTAL_SPACEGROUP:
      {
         VFN_DEBUG_MESSAGE("WXCrystal::OnChangeName():Changing SpaceGroup",6)
         mpCrystal->Init(mpCrystal->GetLatticePar(0),
                         mpCrystal->GetLatticePar(1),
                         mpCrystal->GetLatticePar(2),
                         mpCrystal->GetLatticePar(3),
                         mpCrystal->GetLatticePar(4),
                         mpCrystal->GetLatticePar(5),
                         mpFieldSpacegroup->GetValue(),
                         mpCrystal->GetName());
         mpCrystal->UpdateLatticePar();
         //Update the lattice parameters
         this->CrystUpdate();
         this->Layout();
         return true;
      }
   }
   return false;
}

void WXCrystal::UpdateUI()
{
   VFN_DEBUG_ENTRY("WXCrystal::UpdateUI()",6)
   mpFieldSpacegroup->SetValue(mpCrystal->GetSpaceGroup().GetName());
   #ifdef OBJCRYST_GL
   if(0!=mpCrystalGL) mpCrystalGL->GetParent()->SetTitle(mpCrystal->GetName().c_str());
   #endif
   this->WXRefinableObj::UpdateUI();
   VFN_DEBUG_EXIT("WXCrystal::UpdateUI()",6)
}
Crystal& WXCrystal::GetCrystal(){return *mpCrystal;}
const Crystal& WXCrystal::GetCrystal()const{return *mpCrystal;}

#ifdef OBJCRYST_GL
////////////////////////////////////////////////////////////////////////
//
//    UnitCellMapImport
//
////////////////////////////////////////////////////////////////////////
UnitCellMapImport::UnitCellMapImport(const Crystal&crystal):
mpCrystal(&crystal)
{}
UnitCellMapImport::~UnitCellMapImport(){}
void UnitCellMapImport::GLInitDisplayList(const float minValue,
                                          const REAL xMin,const REAL xMax,
                                          const REAL yMin,const REAL yMax,
                                          const REAL zMin,const REAL zMax) const
{
   VFN_DEBUG_ENTRY("UnitCellMapImport::GLInitDisplayList()",7)
   // Generate triangles
      VFN_DEBUG_MESSAGE("UnitCellMapImport::GLInitDisplayList(): Generate Triangles",7)
      const int nx=mPoints.cols();
      const int ny=mPoints.rows();
      const int nz=mPoints.depth();
      float step[3];
      step[0]=1/(float)nx;
      step[1]=1/(float)nx;
      step[2]=1/(float)nx;
      const int nxMin = (int)(xMin * nx), nxMax = (int)(xMax * nx);
      const int nyMin = (int)(yMin * ny), nyMax = (int)(yMax * ny);
      const int nzMin = (int)(zMin * nz), nzMax = (int)(zMax * nz);
      const int snx = nxMax-nxMin+1, sny = nyMax-nyMin+1, snz = nzMax-nzMin+1;
      const unsigned int ny_nz = ny*nz, sny_snz = sny*snz;
      int i, j, k;
      unsigned int ni, nj, si, sj, sk, sni, snj, sind;
      float x, y, z;

      //create new set of points
      mp4Vector * subPoints = new mp4Vector[snx*sny*snz];
      for(i=nxMin, si=0; i <= nxMax; i++, si++)
      {
         ni = ((nx + i % nx) % nx)*ny_nz;    //this will 'wrap' around any value (negative or positive)
         sni = si*sny_snz;
         for(j=nyMin, sj=0; j <= nyMax; j++, sj++)
         {
            nj = ((ny + j % ny) % ny)*nz;
            snj = sj*snz;
            for(k=nzMin, sk=0; k <= nzMax; k++, sk++)
            {
               sind = sni + snj + sk;
               x = i*step[0]; y = j*step[1]; z = k*step[2];
               mpCrystal->FractionalToOrthonormalCoords(x, y, z);
               subPoints[sind].x = x; subPoints[sind].y = y; subPoints[sind].z = z;
               subPoints[sind].val = mPoints(nz,nj,ni);
            }
         }
      }
      int numOfTriangles;
      const TRIANGLE *pTriangles= MC(snx-1, sny-1, snz-1, step[0], step[1], step[2], minValue, subPoints, numOfTriangles);
   // OpenGL drawing instructions
      VFN_DEBUG_MESSAGE("UnitCellMapImport::GLInitDisplayList(): OpenGL instructions",7)
      /*
      if(showWire) glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
      else glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
      float mccolor[] = {  (float)fcolor.Red()/255.0, 
                    (float)fcolor.Green()/255.0, 
                    (float)fcolor.Blue()/255.0,    1.0};
      glColor4fv(mccolor);
      glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mccolor);
      glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mccolor);
      glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 100.0);
      */
      glBegin(GL_TRIANGLES);
         for(int i=0; i < numOfTriangles; i++)
         {
            for(int j=0; j < 3; j++)
            {
               //VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnPaint():MC1:"<<i<<" "<<j,5)
               glNormal3f(pTriangles[i].norm[j].x, pTriangles[i].norm[j].y, pTriangles[i].norm[j].z);
               glVertex3f(pTriangles[i].p[j].x    ,pTriangles[i].p[j].y    ,pTriangles[i].p[j].z);
            }
         }
      glEnd();

   delete [] subPoints; 
   VFN_DEBUG_EXIT("UnitCellMapImport::GLInitDisplayList()",7)

}
void UnitCellMapImport::ImportGRD(const string&filename)
{
   ifstream ffile(filename.c_str());
   if(!ffile.is_open())
   {     //if file could not be loaded for some reason then exit
      (*fpObjCrystInformUser)("Error opening file: "+filename);
      return;
   }
   //message for reporting errors
   char buff[99];
   ffile.getline(buff, 100);
   float a, b, c, alpha, beta, gamma;
   ffile >>a >>b >>c >>alpha >>beta >>gamma;
   if(!ffile.good()) {  (*fpObjCrystInformUser)("Error reading file: "+filename); return; }
   //compare dimensions with the original crystal and notify the user if not equal
   /*
   float afac = 180/M_PI, limit = 0.0001;
   if((a - mpWXCrystal->GetCrystal().GetLatticePar()(0)) > limit || (b - mpWXCrystal->GetCrystal().GetLatticePar()(1))> limit ||
      (c - mpWXCrystal->GetCrystal().GetLatticePar()(2)) > limit || (alpha - mpWXCrystal->GetCrystal().GetLatticePar()(3)*afac) > limit || 
      (beta - mpWXCrystal->GetCrystal().GetLatticePar()(4)*afac) > limit || (gamma - mpWXCrystal->GetCrystal().GetLatticePar()(5)*afac) > limit )
      if(wxMessageBox(wxString::Format("Cell dimensions in the file do not match those of the crystal loaded:\n\n" +
         wxString("These are the value:\n") + "  Crystal:                     File:\n   a = %f                  a = %f\n" 
         "   b = %f                  b = %f\n   c = %f                   c = %f\n   alpha = %f             alpha = %f\n" +
         "   beta =  %f            beta = %f\n   gamma = %f          gamma = %f\n\nPercent errors are:\n" +
         "   a: %f\n   b: %f\n   c: %f\n   alpha: %f\n   beta:  %f\n   gamma: %f\n\n\n"+ 
         "Continue loading " + filename.c_str() + " ?",
         mpWXCrystal->GetCrystal().GetLatticePar()(0), a,    mpWXCrystal->GetCrystal().GetLatticePar()(1), b, 
         mpWXCrystal->GetCrystal().GetLatticePar()(2), c,    mpWXCrystal->GetCrystal().GetLatticePar()(3)*afac, alpha, 
         mpWXCrystal->GetCrystal().GetLatticePar()(4)*afac, beta,mpWXCrystal->GetCrystal().GetLatticePar()(5)*afac, gamma, 
         fabs(a-mpWXCrystal->GetCrystal().GetLatticePar()(0)) / mpWXCrystal->GetCrystal().GetLatticePar()(0)*100, 
         fabs(b-mpWXCrystal->GetCrystal().GetLatticePar()(1)) / mpWXCrystal->GetCrystal().GetLatticePar()(1)*100, 
         fabs(c-mpWXCrystal->GetCrystal().GetLatticePar()(2)) / mpWXCrystal->GetCrystal().GetLatticePar()(2)*100,
         fabs(alpha-mpWXCrystal->GetCrystal().GetLatticePar()(3)*afac) / mpWXCrystal->GetCrystal().GetLatticePar()(3)*afac*100,
         fabs(beta-mpWXCrystal->GetCrystal().GetLatticePar()(4)*afac ) / mpWXCrystal->GetCrystal().GetLatticePar()(4)*afac*100,
         fabs(gamma-mpWXCrystal->GetCrystal().GetLatticePar()(5)*afac) / mpWXCrystal->GetCrystal().GetLatticePar()(5)*afac*100 ),
         "Cell Dimensions Notice", wxYES_NO | wxCENTRE, (wxWindow*)this) == wxNO) 
       {
         ffile.close();
         return;
       }
   */
   int nx,ny,nz;
   ffile >>nx >>ny >>nz;
   if(!ffile.good()) {  (*fpObjCrystInformUser)("Error reading file: "+filename); return; }
   mPoints.resize(nz,ny,nx);
   for(int i=0; i < nx; i++) {
     for(int j=0; j < ny; j++) {
        for(int k=0; k < nz; k++) {
           ffile >>mPoints(k,j,i);      //reading rhos
        }
     }
   }
   ffile.close();
}
////////////////////////////////////////////////////////////////////////
//
//    UnitCellMapGLList
//
////////////////////////////////////////////////////////////////////////
UnitCellMapGLList::UnitCellMapGLList(const UnitCellMapImport &map,
                                     const float contour,const bool swhowWire,
                                     const float r,const float g,const float b,
                                     const float t):
mpMap(&map),mGLDisplayList(0)
{
   this->SetColour(r,g,b,t);
   this->SetContour(contour);
}

UnitCellMapGLList::~UnitCellMapGLList()
{
}

void UnitCellMapGLList::SetContour(const float value)
{
   mContourValue=value;
   if(0==mGLDisplayList) mGLDisplayList=glGenLists(1);
}

void UnitCellMapGLList::SetColour(const float r,const float g,const float b,
                                 const float t)
{
   mColour[0]=r;
   mColour[1]=g;
   mColour[2]=b;
   mColour[3]=t;
}

void UnitCellMapGLList::Draw()const
{
   if(mShowWire) glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
   else glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
   glColor4fv(mColour);
   glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mColour);
   glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mColour);
   glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 100.0);
   glCallList(mGLDisplayList);
}

////////////////////////////////////////////////////////////////////////
//
//    WXGLCrystalCanvas
//
////////////////////////////////////////////////////////////////////////
BEGIN_EVENT_TABLE(WXGLCrystalCanvas, wxGLCanvas)
   EVT_SIZE             (WXGLCrystalCanvas::OnSize)
   EVT_PAINT            (WXGLCrystalCanvas::OnPaint)
   EVT_ERASE_BACKGROUND (WXGLCrystalCanvas::OnEraseBackground)
   EVT_MOUSE_EVENTS     (WXGLCrystalCanvas::OnMouse)
   EVT_MENU             (ID_GLCRYSTAL_MENU_UPDATE,              WXGLCrystalCanvas::OnUpdate)
   EVT_MENU             (ID_GLCRYSTAL_MENU_CHANGELIMITS,        WXGLCrystalCanvas::OnChangeLimits)
   EVT_MENU             (ID_GLCRYSTAL_MENU_SHOWCRYSTAL,         WXGLCrystalCanvas::OnShowCrystal)     //shows or hides the crystal
   EVT_MENU             (ID_GLCRYSTAL_MENU_LOADFOURIER,         WXGLCrystalCanvas::OnLoadFourier)
   EVT_MENU             (ID_GLCRYSTAL_MENU_CHANGECONTOUR,       WXGLCrystalCanvas::OnChangeContour)
   EVT_MENU             (ID_GLCRYSTAL_MENU_SHOWFOURIER,         WXGLCrystalCanvas::OnShowFourier)
   EVT_MENU             (ID_GLCRYSTAL_MENU_FOURIERCHANGECOLOR,  WXGLCrystalCanvas::OnFourierChangeColor)
   EVT_MENU             (ID_GLCRYSTAL_MENU_SHOWWIRE,            WXGLCrystalCanvas::OnShowWire)
   EVT_MENU             (ID_GLCRYSTAL_MENU_UNLOADFOURIER,       WXGLCrystalCanvas::OnUnloadFourier)   EVT_CHAR             (WXGLCrystalCanvas::OnKeyDown)
   EVT_KEY_DOWN         (WXGLCrystalCanvas::OnKeyDown)
   EVT_KEY_UP           (WXGLCrystalCanvas::OnKeyUp)
   EVT_UPDATE_UI(ID_GLCRYSTAL_UPDATEUI,WXGLCrystalCanvas::OnUpdateUI)
END_EVENT_TABLE()

WXGLCrystalCanvas::WXGLCrystalCanvas(WXCrystal *wxcryst,
                                     wxFrame *parent, wxWindowID id,
                                     const wxPoint &pos,
                                     const wxSize &size):
wxGLCanvas(parent,id,pos,size,wxDEFAULT_FRAME_STYLE),//
mpWXCrystal(wxcryst),mIsGLInit(false),mDist(60),mX0(0),mY0(0),mZ0(0),mViewAngle(15),
mXmin(-.1),mXmax(1.1),mYmin(-.1),mYmax(1.1),mZmin(-.1),mZmax(1.1)
{
   VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::WXGLCrystalCanvas()",3)
   mpPopUpMenu=new wxMenu("Crystal");
   mpPopUpMenu->Append(ID_GLCRYSTAL_MENU_UPDATE, "&Update");
   mpPopUpMenu->Append(ID_GLCRYSTAL_MENU_CHANGELIMITS, "Change display &Limits");
   
   mcXmin = 0.0; mcXmax = 1.0; mcYmin = 0.0; mcYmax = 1.0; mcZmin = 0.0; mcZmax = 1.0;
   minValue = 1.0;
   showFourier = showWireMC = showCrystal = TRUE;
   fcolor.Set(0xFF, 0, 0);
   initMC = FALSE;
   numOfTriangles = 0;
   Triangles = NULL;
   step[0] = step[1] = step[2] = 1.0;
   mpPopUpMenu->Append(ID_GLCRYSTAL_MENU_SHOWCRYSTAL, "Hide Crystal");
   mpPopUpMenu->AppendSeparator();
   mpPopUpMenu->Append(ID_GLCRYSTAL_MENU_LOADFOURIER, "Load Fourier Map");	
   mpPopUpMenu->Append(ID_GLCRYSTAL_MENU_CHANGECONTOUR, "Change Contour Value");
   mpPopUpMenu->Enable(ID_GLCRYSTAL_MENU_CHANGECONTOUR, FALSE);	//disable it for now
   mpPopUpMenu->Append(ID_GLCRYSTAL_MENU_SHOWFOURIER, "Hide Fourier Map");
   mpPopUpMenu->Enable(ID_GLCRYSTAL_MENU_SHOWFOURIER, FALSE);
   mpPopUpMenu->Append(ID_GLCRYSTAL_MENU_FOURIERCHANGECOLOR, "Change Fourier Color");
   mpPopUpMenu->Enable(ID_GLCRYSTAL_MENU_FOURIERCHANGECOLOR, FALSE);
   mpPopUpMenu->Append(ID_GLCRYSTAL_MENU_SHOWWIRE, "Show Filled");
   mpPopUpMenu->Enable(ID_GLCRYSTAL_MENU_SHOWWIRE, FALSE);
   mpPopUpMenu->Append(ID_GLCRYSTAL_MENU_UNLOADFOURIER, "Unload Fourier Map");
   mpPopUpMenu->Enable(ID_GLCRYSTAL_MENU_UNLOADFOURIER, FALSE);
}

WXGLCrystalCanvas::~WXGLCrystalCanvas()
{
   mpWXCrystal->NotifyCrystalGLDelete();
   this->MCCleanUp();
}

void WXGLCrystalCanvas::OnExit(wxCommandEvent &event)
{
   
}

void WXGLCrystalCanvas::OnPaint(wxPaintEvent &event)
{
   VFN_DEBUG_ENTRY("WXGLCrystalCanvas::OnPaint()",7)
   // This means that another update of the display list is being done, so...
   if(true==mpWXCrystal->GLDisplayListIsLocked()) return;
   
   wxPaintDC dc(this);
   PrepareDC(dc);
   this->GetParent()->PrepareDC(dc);

   #ifndef __WXMOTIF__
   if (!GetContext()) return;
   #endif

   this->SetCurrent();
   if(false==mIsGLInit)
   {
      mIsGLInit=true;
      this->InitGL();
   }

   //glMatrixMode( GL_PROJECTION );
   //glLoadIdentity();
   //gluPerspective( 15, .5, 1, 100 );
   glMatrixMode( GL_MODELVIEW );

   //clear
   glClearColor(0.0f, 0.0f, 0.0f, 1.0f);//Black background
   glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

   //Orientation using the trackball
   GLfloat m[4][4];
   glLoadIdentity();
   glTranslatef( 0, 0, -mDist );
   build_rotmatrix( m,mQuat);
   glMultMatrixf( &m[0][0] );
   glTranslatef( mX0, mY0, mZ0 );
   
   //Draw
      // another update of the display list is being done, so...
      if(true==mpWXCrystal->GLDisplayListIsLocked())
      {
         VFN_DEBUG_EXIT("WXGLCrystalCanvas::OnPaint()",7)
         return;
      }
   
   //drawing triangles returned by MC
   glPushMatrix();
   REAL xc=(mXmin+mXmax)/2.;  //this is also done in Crystal
   REAL yc=(mYmin+mYmax)/2.;
   REAL zc=(mZmin+mZmax)/2.;
   mpWXCrystal->GetCrystal().FractionalToOrthonormalCoords(xc, yc, zc);
   glTranslatef(-xc, -yc, -zc);
   if(showFourier)
   {
      //if(Triangles==NULL) break;
      if(showWireMC) glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
      else glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
      float mccolor[] = {  (float)fcolor.Red()/255.0, 
                    (float)fcolor.Green()/255.0, 
                    (float)fcolor.Blue()/255.0,    1.0};
      glColor4fv(mccolor);
      glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mccolor);
      glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mccolor);
      glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 100.0);
      glBegin(GL_TRIANGLES);
         for(int i=0; i < numOfTriangles; i++)
         {
            for(int j=0; j < 3; j++)
            {
               //VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnPaint():MC1:"<<i<<" "<<j,5)
               glNormal3f(Triangles[i].norm[j].x, Triangles[i].norm[j].y, Triangles[i].norm[j].z);
               glVertex3f(Triangles[i].p[j].x,Triangles[i].p[j].y,Triangles[i].p[j].z);
            }
         }
      glEnd();
   }
   if(initMC && cdial->IsShown())
   {
      glLineWidth(6.0);
      float bound_color[] = {1.0, 0.2, 0.2, 1.0};
      glColor4fv(bound_color);
      glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, bound_color);
      glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, bound_color);
      glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
      cdial->DrawBoundingBox();
      glLineWidth(1.0);
   }
   glPopMatrix();
  
   if(showCrystal)
   {
      glCallList(mpWXCrystal->GrabCrystalGLDisplayList());  //Draw Crystal
      mpWXCrystal->ReleaseCrystalGLDisplayList();
   }
   glFlush();
   SwapBuffers();
   VFN_DEBUG_EXIT("WXGLCrystalCanvas::OnPaint():End",7)
}

void WXGLCrystalCanvas::OnSize(wxSizeEvent& event)
{
   VFN_DEBUG_ENTRY("WXGLCrystalCanvas::OnSize()",7)
   int width, height;
   GetClientSize(& width, & height);

   #ifndef __WXMOTIF__
   if (GetContext())
   #endif
   {
      SetCurrent();
      glViewport(0, 0, width, height);
      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();
      VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnSize():"<<mViewAngle<<","<<width<<","<<height<<","<<mDist,2)
      if( (width>0)&&(height>0)) //in case the window is docked...
         gluPerspective(mViewAngle,(float)width/(float)height,1.f,2.*mDist);   
   }
   VFN_DEBUG_EXIT("WXGLCrystalCanvas::OnSize():End",7)
}

void WXGLCrystalCanvas::OnEraseBackground(wxEraseEvent& event)
{
}

void WXGLCrystalCanvas::OnKeyDown(wxKeyEvent& event)
{
   VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnKeyDown()",2)
   switch(event.GetKeyCode())
   {
      case(45):// +
      {
         VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnKeyDown():Bigger",2)
         int width, height;
         GetClientSize(& width, & height);
         mDist *= 1.05;
         SetCurrent();
         glMatrixMode(GL_PROJECTION);
         glLoadIdentity();
         if( (width>0)&&(height>0)) //in case size is null...
            gluPerspective(mViewAngle,(float)width/(float)height,
                        (mDist>100)?(mDist-100):1.,mDist+100);   
         Refresh(FALSE);
         break;
      }
      case(43):// -
      {
         VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnKeyDown():Smaller",2)
         int width, height;
         GetClientSize(& width, & height);
         mDist *= .95;
         SetCurrent();
         glMatrixMode(GL_PROJECTION);
         glLoadIdentity();
         if( (width>0)&&(height>0)) //in case size is null...
            gluPerspective(mViewAngle,(float)width/(float)height,
                        (mDist>100)?(mDist-100):1.,mDist+100);   
         Refresh(FALSE);
         break;
      }
      case(WXK_INSERT): mY0 += 0.1; Refresh(FALSE); break;
      case(WXK_DELETE): mY0 -= 0.1; Refresh(FALSE); break;
      case(WXK_HOME): mX0 -= 0.1; Refresh(FALSE); break;
      case(WXK_END): mX0 += 0.1; Refresh(FALSE); break;
      case(WXK_PRIOR): mZ0 -= 0.1; Refresh(FALSE); break;
      case(WXK_NEXT): mZ0 += 0.1; Refresh(FALSE); break;
      case(52):// 4
      {
         VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnKeyDown():rotate left",2)
         float spin_quat[4];
         trackball(spin_quat,0,0,-.05,0);
         add_quats( spin_quat, mQuat, mQuat );
         Refresh(FALSE);
         break;
      }
      case(54):// 6
      {
         VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnKeyDown():rotate right",2)
         float spin_quat[4];
         trackball(spin_quat,0,0,.05,0);
         add_quats( spin_quat, mQuat, mQuat );
         Refresh(FALSE);
         break;
      }
      case(50):// 2
      {
         VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnKeyDown():rotate down",2)
         float spin_quat[4];
         trackball(spin_quat,0,0,0,-.05);
         add_quats( spin_quat, mQuat, mQuat );
         Refresh(FALSE);
         break;
      }
      case(56):// 8
      {
         VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnKeyDown():rotate up",2)
         float spin_quat[4];
         trackball(spin_quat,0,0,0,.05);
         add_quats( spin_quat, mQuat, mQuat );
         Refresh(FALSE);
         break;
      }
      case(68):// D
      {
         VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnKeyDown():rotate left",2)
         float spin_quat[4];
         trackball(spin_quat,0,0,-.05,0);
         add_quats( spin_quat, mQuat, mQuat );
         Refresh(FALSE);
         break;
      }
      case(70):// F
      {
         VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnKeyDown():rotate right",2)
         float spin_quat[4];
         trackball(spin_quat,0,0,.05,0);
         add_quats( spin_quat, mQuat, mQuat );
         Refresh(FALSE);
         break;
      }
      case(67):// C
      {
         VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnKeyDown():rotate down",2)
         float spin_quat[4];
         trackball(spin_quat,0,0,0,-.05);
         add_quats( spin_quat, mQuat, mQuat );
         Refresh(FALSE);
         break;
      }
      case(82):// R
      {
         VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnKeyDown():rotate up",2)
         float spin_quat[4];
         trackball(spin_quat,0,0,0,.05);
         add_quats( spin_quat, mQuat, mQuat );
         Refresh(FALSE);
         break;
      }
   }
}

void WXGLCrystalCanvas::OnKeyUp(wxKeyEvent& event)
{
   VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnKeyUp():"<<event.GetKeyCode(),2)
}

void WXGLCrystalCanvas::OnEnterWindow( wxMouseEvent& event )
{
   VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnEnterWindow()",10)
}

void WXGLCrystalCanvas::OnMouse( wxMouseEvent& event )
{
   VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnMouse()",7)
   if (event.Dragging())
   {
      int width, height;
      GetClientSize(& width, & height);
      if(event.LeftIsDown())
      {
         if(event.ShiftDown())
         {
            VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnMouse():Dragging Left Button",2)
            const float v= (mTrackBallLastY-event.GetY())/(float)height;
            const float h= (mTrackBallLastX-event.GetX())/(float)width;
            GLfloat m[4][4];
            build_rotmatrix( m,mQuat);
            const float dx=-h*10,dy=v*10;
            mX0 += m[0][0]* dx +m[0][1]*dy;
            mY0 += m[1][0]* dx +m[1][1]*dy;
            mZ0 += m[2][0]* dx +m[2][1]*dy;
            Refresh(FALSE);
         }
         else
         {
            VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnMouse():Dragging Left Button",2)
            // drag in progress, simulate trackball
            float spin_quat[4];
            trackball(spin_quat,
            (2.0*mTrackBallLastX -       width) / (width+.001),  //normalizing from -1 to 1
            (     height - 2.0*mTrackBallLastY) / (height+.001),
            (     2.0*event.GetX() - width) / (width+.001),
            (    height - 2.0*event.GetY()) / (height+.001));

            add_quats( spin_quat, mQuat, mQuat );
            Refresh(FALSE);
         }
      }
      if(event.MiddleIsDown())
      {
         VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnMouse():Dragging Middle Button",2)
         const float v= (mTrackBallLastY-event.GetY())/(float)height;
         const float h= (mTrackBallLastX-event.GetX())/(float)width;
         
         mDist *= (1.+v)/(1.+h);
         mViewAngle *=(1.+h);
         SetCurrent();
         glMatrixMode(GL_PROJECTION);
         glLoadIdentity();
         if( (width>0)&&(height>0)) //in case size is null...
            gluPerspective(mViewAngle,(float)width/(float)height,
                        (mDist>100)?(mDist-100):1.,mDist+100);   
         Refresh(FALSE);
         VFN_DEBUG_MESSAGE(mViewAngle <<" "<<mDist,2)
      }
   }
   if(event.RightIsDown())
   {
      VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnMouse():Right Button",2)
      this->PopupMenu(mpPopUpMenu, event.GetX(), event.GetY() );
   }

   mTrackBallLastX = event.GetX();
   mTrackBallLastY = event.GetY();
}

void WXGLCrystalCanvas::OnUpdate(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_ENTRY("WXGLCrystalCanvas::OnUpdate()",4)
   mpWXCrystal->UpdateGL(false,mXmin,mXmax,mYmin,mYmax,mZmin,mZmax);
   VFN_DEBUG_EXIT("WXGLCrystalCanvas::OnUpdate()",4)
}

void WXGLCrystalCanvas::CrystUpdate()
{
   VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::CrystUpdate()",7)
   wxUpdateUIEvent event(ID_GLCRYSTAL_UPDATEUI);
   wxPostEvent(this,event);
}

void WXGLCrystalCanvas::OnUpdateUI(wxUpdateUIEvent& WXUNUSED(event))
{
   VFN_DEBUG_ENTRY("WXGLCrystalCanvas::OnUpdateUI()",5)
   this->Refresh(false);
   VFN_DEBUG_EXIT("WXGLCrystalCanvas::OnUpdateUI()",5)
}

void WXGLCrystalCanvas::InitGL()
{
   VFN_DEBUG_ENTRY("WXGLCrystalCanvas::InitGL()",8)
   this->SetCurrent();
    
   glEnable(GL_DEPTH_TEST);
      glEnable(GL_LIGHTING);
      
      const GLfloat color_Ambient [] = {0.3, 0.3, 0.3, 1.00}; 
      const GLfloat color_Diffuse [] = {0.6, 0.6, 0.6, 1.00}; 
      const GLfloat color_Specular[] = {0.8, 0.8, 0.8, 1.00}; 
      
      glLightfv( GL_LIGHT0, GL_AMBIENT,  color_Ambient); 
      glLightfv( GL_LIGHT0, GL_DIFFUSE,  color_Diffuse); 
      glLightfv( GL_LIGHT0, GL_SPECULAR, color_Specular); 
      glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, 1.0); 

      glEnable( GL_LIGHT0 ); 

      glEnable(GL_NORMALIZE);
      glHint(GL_PERSPECTIVE_CORRECTION_HINT,GL_FASTEST);
      glHint(GL_POLYGON_SMOOTH_HINT,GL_FASTEST);
   
   //Initialize Trackball
   trackball(mQuat,0.,0.,0.,0.);
   
   wxSizeEvent event;
   this->OnSize(event);
   
   //First display
   this->CrystUpdate();
   VFN_DEBUG_EXIT("WXGLCrystalCanvas::InitGL()",8)
}
void WXGLCrystalCanvas::OnChangeLimits(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnChangeLimits():End",10)
   double xMin,xMax,yMin,yMax,zMin,zMax;
   {
      wxString str;
      str << mXmin;
      wxTextEntryDialog limitDialog(this,"Enter the minimum x value (reduced coords)",
                              "x min",str,wxOK | wxCANCEL);
      if(wxID_OK!=limitDialog.ShowModal())
      {
         VFN_DEBUG_EXIT("WXGLCrystalCanvas::OnChangeLimits():Cancelled",6)
         return;
      }
      limitDialog.GetValue().ToDouble(&xMin);
   }
   {
      wxString str;
      str << mXmax;
      wxTextEntryDialog limitDialog(this,"Enter the maximum x value (reduced coords)",
                              "x max",str,wxOK | wxCANCEL);
      if(wxID_OK!=limitDialog.ShowModal())
      {
         VFN_DEBUG_EXIT("WXGLCrystalCanvas::OnChangeLimits():Cancelled",6)
         return;
      }
      limitDialog.GetValue().ToDouble(&xMax);
   }
   if(xMax<=xMin)
   {
      wxMessageDialog dumbUser(this,"max <= min !!!",
                               "Whooops",wxOK|wxICON_EXCLAMATION);
      dumbUser.ShowModal();
      return;
   }
   {
      wxString str;
      str << mYmin;
      wxTextEntryDialog limitDialog(this,"Enter the minimum y value (reduced coords)",
                              "y min",str,wxOK | wxCANCEL);
      if(wxID_OK!=limitDialog.ShowModal())
      {
         VFN_DEBUG_EXIT("WXGLCrystalCanvas::OnChangeLimits():Cancelled",6)
         return;
      }
      limitDialog.GetValue().ToDouble(&yMin);
   }
   {
      wxString str;
      str << mYmax;
      wxTextEntryDialog limitDialog(this,"Enter the maximum y value (reduced coords)",
                              "y max",str,wxOK | wxCANCEL);
      if(wxID_OK!=limitDialog.ShowModal())
      {
         VFN_DEBUG_EXIT("WXGLCrystalCanvas::OnChangeLimits():Cancelled",6)
         return;
      }
      limitDialog.GetValue().ToDouble(&yMax);
   }
   if(yMax<=yMin)
   {
      wxMessageDialog dumbUser(this,"max <= min !!!",
                               "Whooops",wxOK|wxICON_EXCLAMATION);
      dumbUser.ShowModal();
      return;
   }
   {
      wxString str;
      str << mZmin;
      wxTextEntryDialog limitDialog(this,"Enter the minimum z value (reduced coords)",
                              "z min",str,wxOK | wxCANCEL);
      if(wxID_OK!=limitDialog.ShowModal())
      {
         VFN_DEBUG_EXIT("WXGLCrystalCanvas::OnChangeLimits():Cancelled",6)
         return;
      }
      limitDialog.GetValue().ToDouble(&zMin);
   }
   {
      wxString str;
      str << mZmax;
      wxTextEntryDialog limitDialog(this,"Enter the maximum z value (reduced coords)",
                              "z max",str,wxOK | wxCANCEL);
      if(wxID_OK!=limitDialog.ShowModal())
      {
         VFN_DEBUG_EXIT("WXGLCrystalCanvas::OnChangeLimits():Cancelled",6)
         return;
      }
      limitDialog.GetValue().ToDouble(&zMax);
   }
   if(zMax<=zMin)
   {
      wxMessageDialog dumbUser(this,"max <= min !!!",
                               "Whooops",wxOK|wxICON_EXCLAMATION);
      dumbUser.ShowModal();
      return;
   }
   mXmin=xMin;mXmax=xMax;
   mYmin=yMin;mYmax=yMax;
   mZmin=zMin;mZmax=zMax;
   mpWXCrystal->UpdateGL(false,mXmin,mXmax,mYmin,mYmax,mZmin,mZmax);
   this->CrystUpdate();
}
#endif
void WXGLCrystalCanvas::OnShowCrystal()
{
   if(showCrystal) mpPopUpMenu->SetLabel(ID_GLCRYSTAL_MENU_SHOWCRYSTAL, "Show Crystal");
   else mpPopUpMenu->SetLabel(ID_GLCRYSTAL_MENU_SHOWCRYSTAL, "Hide Crystal");
   showCrystal = !showCrystal;
   this->CrystUpdate();
}

//read fourier map from 0.0 to 1.0 on all axis
void WXGLCrystalCanvas::OnLoadFourier()
{
   wxFileDialog fd((wxWindow*)this, "Choose a file containing a Fourier Map",
           "", "", "Fourier Map files (*.grd)|*.grd", wxOPEN | wxHIDE_READONLY | wxFILE_MUST_EXIST);
   //if okay then read Fourier map, run MC on it and display the triangles
   if(fd.ShowModal() == wxID_OK)
   {
      this->LoadFourier((string)(fd.GetFilename().c_str()));
   }
}
void WXGLCrystalCanvas::LoadFourier(const string&filename)
{
   int err;
   ifstream ffile(filename.c_str());
   if(!ffile.is_open())
   {     //if file could not be loaded for some reason then exit
      wxMessageDialog error_open((wxWindow*)this, "Error opening file " +
                    wxString(filename.c_str()), "File Open Error");
      err = error_open.ShowModal();
      return;
   }
   //message for reporting errors
   wxMessageDialog errmsg((wxWindow*)this, "Error reading    ", 
                    "File Reading Error");
   char buff[99];
   ffile.getline(buff, 100);
   float a, b, c, alpha, beta, gamma, n[3];
   ffile >>a >>b >>c >>alpha >>beta >>gamma;
   if(!ffile.good()) {  err = errmsg.ShowModal(); ffile.close(); return; }
   //compare dimensions with the original crystal and notify the user if not equal
   float afac = 180/M_PI, limit = 0.0001;
   if((a - mpWXCrystal->GetCrystal().GetLatticePar()(0)) > limit || (b - mpWXCrystal->GetCrystal().GetLatticePar()(1))> limit ||
      (c - mpWXCrystal->GetCrystal().GetLatticePar()(2)) > limit || (alpha - mpWXCrystal->GetCrystal().GetLatticePar()(3)*afac) > limit || 
      (beta - mpWXCrystal->GetCrystal().GetLatticePar()(4)*afac) > limit || (gamma - mpWXCrystal->GetCrystal().GetLatticePar()(5)*afac) > limit )
      if(wxMessageBox(wxString::Format("Cell dimensions in the file do not match those of the crystal loaded:\n\n" +
         wxString("These are the value:\n") + "  Crystal:                     File:\n   a = %f                  a = %f\n" 
         "   b = %f                  b = %f\n   c = %f                   c = %f\n   alpha = %f             alpha = %f\n" +
         "   beta =  %f            beta = %f\n   gamma = %f          gamma = %f\n\nPercent errors are:\n" +
         "   a: %f\n   b: %f\n   c: %f\n   alpha: %f\n   beta:  %f\n   gamma: %f\n\n\n"+ 
         "Continue loading " + filename.c_str() + " ?",
         mpWXCrystal->GetCrystal().GetLatticePar()(0), a,    mpWXCrystal->GetCrystal().GetLatticePar()(1), b, 
         mpWXCrystal->GetCrystal().GetLatticePar()(2), c,    mpWXCrystal->GetCrystal().GetLatticePar()(3)*afac, alpha, 
         mpWXCrystal->GetCrystal().GetLatticePar()(4)*afac, beta,mpWXCrystal->GetCrystal().GetLatticePar()(5)*afac, gamma, 
         fabs(a-mpWXCrystal->GetCrystal().GetLatticePar()(0)) / mpWXCrystal->GetCrystal().GetLatticePar()(0)*100, 
         fabs(b-mpWXCrystal->GetCrystal().GetLatticePar()(1)) / mpWXCrystal->GetCrystal().GetLatticePar()(1)*100, 
         fabs(c-mpWXCrystal->GetCrystal().GetLatticePar()(2)) / mpWXCrystal->GetCrystal().GetLatticePar()(2)*100,
         fabs(alpha-mpWXCrystal->GetCrystal().GetLatticePar()(3)*afac) / mpWXCrystal->GetCrystal().GetLatticePar()(3)*afac*100,
         fabs(beta-mpWXCrystal->GetCrystal().GetLatticePar()(4)*afac ) / mpWXCrystal->GetCrystal().GetLatticePar()(4)*afac*100,
         fabs(gamma-mpWXCrystal->GetCrystal().GetLatticePar()(5)*afac) / mpWXCrystal->GetCrystal().GetLatticePar()(5)*afac*100 ),
         "Cell Dimensions Notice", wxYES_NO | wxCENTRE, (wxWindow*)this) == wxNO) 
       {
         ffile.close();
         return;
       }
   ffile >>n[0] >>n[1] >>n[2];
   if(!ffile.good()) {  err = errmsg.ShowModal(); ffile.close(); return; }
   nx = (int)n[0]; ny = (int)n[1]; nz = (int)n[2];
   int all = nx*ny*nz;              
   if(initMC) delete [] mcPoints;                     //free space from last time
   mcPoints = new float[all];
   step[0] = 1/n[0]; step[1] = 1/n[1]; step[2] = 1/n[2]; //init stepsize
   //READ POINTS
   wxProgressDialog prd("Reading Fourier Map from file " + wxString(filename.c_str()),
      "Reading data... ", nx, (wxWindow*)this, 
      wxPD_AUTO_HIDE | wxPD_APP_MODAL | wxPD_ESTIMATED_TIME | wxPD_REMAINING_TIME );
   unsigned int ni, nj;
   for(int i=0; i < nx; i++) {
     ni = i*ny*nz;
     for(int j=0; j < ny; j++) {
        nj = j*nz;
        for(int k=0; k < nz; k++) {
           ffile >>mcPoints[ni + nj + k];      //reading rhos
        }
     }
     prd.Update(i);
   }
   ffile.close();


   //ask the user for contour value
   //wxTextEntryDialog *cted = new wxTextEntryDialog((wxWindow*)this,"Enter value: ",
   //    "Enter contour value for MC ", wxString::Format("%f", minValue), wxOK | wxCENTRE);
   //err = cted->ShowModal();             //err == wxID_CANCEL should not happen: no cancel button
   //minValue = (float)atof(cted->GetValue().c_str());
   //delete cted;

   if(initMC == FALSE) cdial = new ContourDialog(this);  //asks for contour value and runs MC
   cdial->GetContour(FALSE);

   //enable other options in the pop-up menu:
   mpPopUpMenu->Enable(ID_GLCRYSTAL_MENU_CHANGECONTOUR, TRUE);
   mpPopUpMenu->Enable(ID_GLCRYSTAL_MENU_SHOWFOURIER, TRUE);
   mpPopUpMenu->Enable(ID_GLCRYSTAL_MENU_FOURIERCHANGECOLOR, TRUE);
   mpPopUpMenu->Enable(ID_GLCRYSTAL_MENU_UNLOADFOURIER, TRUE);
   mpPopUpMenu->Enable(ID_GLCRYSTAL_MENU_SHOWWIRE, TRUE);
   initMC = TRUE;
   this->CrystUpdate();
}

void WXGLCrystalCanvas::OnChangeContour()
{
   //ask the user for new contour value
   cdial->GetContour(TRUE);
}

void WXGLCrystalCanvas::RunMC()
{
   VFN_DEBUG_ENTRY("WXGLCrystalCanvas::RunMC()",7)
   //free memory -- this caused problems when changing contour values, but now it seems to work...
   if(initMC && Triangles != NULL) { delete [] Triangles; Triangles=0;} 
   VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::RunMC():0",7)
   int nxMin = (int)(mcXmin * nx), nxMax = (int)(mcXmax * nx);
   int nyMin = (int)(mcYmin * ny), nyMax = (int)(mcYmax * ny);
   int nzMin = (int)(mcZmin * nz), nzMax = (int)(mcZmax * nz);
   int snx = nxMax-nxMin+1, sny = nyMax-nyMin+1, snz = nzMax-nzMin+1;
   unsigned int ny_nz = ny*nz, sny_snz = sny*snz;
   int i, j, k;
   unsigned int ni, nj, si, sj, sk, sni, snj, sind;
   float x, y, z;

   //create new set of points
   VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::RunMC():1",7)
   mp4Vector * subPoints = new mp4Vector[snx*sny*snz];
   for(i=nxMin, si=0; i <= nxMax; i++, si++)
   {
      ni = ((nx + i % nx) % nx)*ny_nz;    //this will 'wrap' around any value (negative or positive)
      sni = si*sny_snz;
      for(j=nyMin, sj=0; j <= nyMax; j++, sj++)
      {
         nj = ((ny + j % ny) % ny)*nz;
         snj = sj*snz;
         for(k=nzMin, sk=0; k <= nzMax; k++, sk++)
         {
            sind = sni + snj + sk;
            x = i*step[0]; y = j*step[1]; z = k*step[2];
            mpWXCrystal->GetCrystal().FractionalToOrthonormalCoords(x, y, z);
            subPoints[sind].x = x; subPoints[sind].y = y; subPoints[sind].z = z;
            subPoints[sind].val = mcPoints[ni + nj + (nz + k % nz) % nz];
         }
      }
   }
   Triangles = MC(snx-1, sny-1, snz-1, step[0], step[1], step[2], minValue, subPoints, numOfTriangles);
   VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::RunMC():2",7)
   delete [] subPoints; 
   VFN_DEBUG_EXIT("WXGLCrystalCanvas::RunMC()",7)
   this->CrystUpdate();
}

void WXGLCrystalCanvas::MCCleanUp()
{
   showFourier = showWireMC = showCrystal = TRUE;
   nx = ny = nz = 0;
   numOfTriangles = 0;
   mcXmin = 0.0; mcXmax = 1.0; mcYmin = 0.0; mcYmax = 1.0; mcZmin = 0.0; mcZmax = 1.0;
   minValue = 1.0;
   if(initMC)
   {
      delete [] mcPoints;
      delete [] Triangles;
   }
   initMC = FALSE;
   mcPoints = NULL;
   Triangles = NULL;
   step[0] = step[1] = step[2] = 1.0;
   mpPopUpMenu->SetLabel(ID_GLCRYSTAL_MENU_SHOWCRYSTAL, "Hide Crystal");
   mpPopUpMenu->Enable(ID_GLCRYSTAL_MENU_CHANGECONTOUR, FALSE);      //disable all of these
   mpPopUpMenu->Enable(ID_GLCRYSTAL_MENU_SHOWFOURIER, FALSE);
   mpPopUpMenu->SetLabel(ID_GLCRYSTAL_MENU_SHOWFOURIER, "Hide Fourier Map");
   mpPopUpMenu->Enable(ID_GLCRYSTAL_MENU_UNLOADFOURIER, FALSE);
   mpPopUpMenu->Enable(ID_GLCRYSTAL_MENU_FOURIERCHANGECOLOR, FALSE);
   mpPopUpMenu->SetLabel(ID_GLCRYSTAL_MENU_SHOWWIRE, "Show Filled");
   mpPopUpMenu->Enable(ID_GLCRYSTAL_MENU_SHOWWIRE, FALSE);
   fcolor.Set(0xFF, 0, 0);
}

void WXGLCrystalCanvas::OnShowFourier()
{
   if(showFourier == TRUE) mpPopUpMenu->SetLabel(ID_GLCRYSTAL_MENU_SHOWFOURIER, "Show Fourier Map");
   else mpPopUpMenu->SetLabel(ID_GLCRYSTAL_MENU_SHOWFOURIER, "Hide Fourier Map");
   showFourier = !showFourier;
   this->CrystUpdate();
}

void WXGLCrystalCanvas::OnFourierChangeColor()
{
   wxColor ncolor;
   ncolor = wxGetColourFromUser((wxWindow*)this, fcolor);   
   if(ncolor.Ok()) 
   { //if user pressed OK
      fcolor = ncolor;
      this->CrystUpdate();
   }
}

void WXGLCrystalCanvas::OnUnloadFourier()
{
   wxMessageDialog * msure = new wxMessageDialog((wxWindow*)this,
     "Are you sure you want to unload all Fourier Map Data?", "Unload Fourier Map", wxYES_NO | wxNO_DEFAULT |
     wxICON_QUESTION );
   if(msure->ShowModal() == wxID_YES)
   {
      MCCleanUp();
      this->CrystUpdate();
   }
   delete msure;
}

void WXGLCrystalCanvas::OnShowWire()
{
   if(showWireMC == TRUE) mpPopUpMenu->SetLabel(ID_GLCRYSTAL_MENU_SHOWWIRE, "Show as Wireframe");    
   else mpPopUpMenu->SetLabel(ID_GLCRYSTAL_MENU_SHOWWIRE, "Show Filled");
   showWireMC = !showWireMC;
   this->CrystUpdate();
}

///////////////////////////
//class ContourDialog
///////////////////////////

enum {ID_BOUNDINGVOLUME_TEXT_CHANGE=1000};

BEGIN_EVENT_TABLE(WXGLCrystalCanvas::ContourDialog, wxDialog)
   EVT_BUTTON(wxOK, WXGLCrystalCanvas::ContourDialog::OnOk)
   EVT_BUTTON(wxCANCEL, WXGLCrystalCanvas::ContourDialog::OnCancel)  
   EVT_CLOSE(WXGLCrystalCanvas::ContourDialog::Closing)
   EVT_TEXT(ID_BOUNDINGVOLUME_TEXT_CHANGE, WXGLCrystalCanvas::ContourDialog::BoundingTextChange)
END_EVENT_TABLE()

WXGLCrystalCanvas::ContourDialog::ContourDialog(WXGLCrystalCanvas * parent) : 
   wxDialog((wxWindow*)parent, -1, "Enter Contour Value", wxDefaultPosition, 
   wxSize(300, 300), wxCAPTION   ), shown(FALSE)
{ 
   this->parent = parent;
   stText[0] = new wxStaticText((wxWindow*)this, -1, "Enter contour value", wxPoint(20, 10));
   conValue = new wxTextCtrl((wxWindow*)this, -1, wxString::Format("%f",parent->minValue), wxPoint(20, 30));

   boundBox = new wxStaticBox((wxWindow*)this, -1, "Enter bounding volume values", wxPoint(10, 60), wxSize(275, 160));

   stText[1] = new wxStaticText((wxWindow*)this, -1, "X min:", wxPoint(20, 100));
   bound[0] = new wxTextCtrl((wxWindow*)this,ID_BOUNDINGVOLUME_TEXT_CHANGE,
                    wxString::Format("%f",parent->mcXmin),wxPoint(50,100),wxSize(70,20));
   stText[2] = new wxStaticText((wxWindow*)this, -1, "X max:", wxPoint(130, 100));
   bound[1] = new wxTextCtrl((wxWindow*)this,ID_BOUNDINGVOLUME_TEXT_CHANGE,
                    wxString::Format("%f",parent->mcXmin),wxPoint(165,100),wxSize(70,20));
   stText[3] = new wxStaticText((wxWindow*)this, -1, "Y min:", wxPoint(20, 140));
   bound[2] = new wxTextCtrl((wxWindow*)this,ID_BOUNDINGVOLUME_TEXT_CHANGE,
                    wxString::Format("%f",parent->mcXmin),wxPoint(50,140),wxSize(70,20));
   stText[4] = new wxStaticText((wxWindow*)this, -1, "Y max:", wxPoint(130, 140));
   bound[3] = new wxTextCtrl((wxWindow*)this,ID_BOUNDINGVOLUME_TEXT_CHANGE,
                    wxString::Format("%f",parent->mcXmin),wxPoint(165,140),wxSize(70,20));
   stText[5] = new wxStaticText((wxWindow*)this, -1, "Z min:", wxPoint(20, 180));
   bound[4] = new wxTextCtrl((wxWindow*)this,ID_BOUNDINGVOLUME_TEXT_CHANGE,
                    wxString::Format("%f",parent->mcXmin),wxPoint(50,180),wxSize(70,20));
   stText[6] = new wxStaticText((wxWindow*)this, -1, "Z max:", wxPoint(130, 180));
   bound[5] = new wxTextCtrl((wxWindow*)this,ID_BOUNDINGVOLUME_TEXT_CHANGE,
                    wxString::Format("%f",parent->mcXmin),wxPoint(165,180),wxSize(70,20));

   butOk = new wxButton((wxWindow*)this, wxOK, "OK", wxPoint(40, 230));
   butCancel = new wxButton((wxWindow*)this, wxCANCEL, "Cancel", wxPoint(170, 230));
   butCancel->Enable(FALSE);
   CenterOnParent();
}

//use this instead of ShowModal()
void WXGLCrystalCanvas::ContourDialog::GetContour(bool showCancelFlag)
{
   if(shown) { Show(TRUE); conValue->SetFocus(); return; }
   butCancel->Enable(showCancelFlag);
   conValue->SetValue(wxString::Format("%f", parent->minValue));
   bound[0]->SetValue(wxString::Format("%f", parent->mcXmin));
   bound[1]->SetValue(wxString::Format("%f", parent->mcXmax));
   bound[2]->SetValue(wxString::Format("%f", parent->mcYmin));
   bound[3]->SetValue(wxString::Format("%f", parent->mcYmax));
   bound[4]->SetValue(wxString::Format("%f", parent->mcZmin));
   bound[5]->SetValue(wxString::Format("%f", parent->mcZmax));
   shown = TRUE;
   conValue->SetFocus();                           //select the textbox with minValue
   conValue->SetSelection(0, conValue->GetValue().Length());
   for(int i=0; i < 6; i++)
      bound[i]->SetSelection(0, bound[i]->GetValue().Length());
   //parent->CrystUpdate();
   Show(TRUE);
}

void WXGLCrystalCanvas::ContourDialog::OnOk()
{ 
   float minx = atof(bound[0]->GetValue().c_str()), maxx = atof(bound[1]->GetValue().c_str()),
         miny = atof(bound[2]->GetValue().c_str()), maxy = atof(bound[3]->GetValue().c_str()),
         minz = atof(bound[4]->GetValue().c_str()), maxz = atof(bound[5]->GetValue().c_str());
   if(minx >= maxx || miny >= maxy || minz >= maxz)
   {
      wxMessageBox("Minimum value has to be less than the maximum!", "Bounding volume error");
      return;
   }
   parent->minValue = atof(conValue->GetValue().c_str());
   parent->mcXmin = minx; parent->mcXmax = maxx;
   parent->mcYmin = miny; parent->mcYmax = maxy;
   parent->mcZmin = minz; parent->mcZmax = maxz;
   parent->RunMC();
   //parent->CrystUpdate();
   Show(FALSE);
   shown = FALSE;
}

void WXGLCrystalCanvas::ContourDialog::OnCancel()
{ 
   Show(FALSE);
   shown = FALSE;
   //parent->CrystUpdate();
}

void WXGLCrystalCanvas::ContourDialog::Closing()
{ /*  Does nothing so user cannot kill the window without pushing OK or Cancel  */ }

void WXGLCrystalCanvas::ContourDialog::BoundingTextChange()
{
  //parent->CrystUpdate(); //everytime new value, update the window
}

bool WXGLCrystalCanvas::ContourDialog::IsShown() const
{  return shown;  }

//draws a box
void WXGLCrystalCanvas::ContourDialog::DrawBoundingBox()
{
   float minx = atof(bound[0]->GetValue().c_str()), maxx = atof(bound[1]->GetValue().c_str()),
         miny = atof(bound[2]->GetValue().c_str()), maxy = atof(bound[3]->GetValue().c_str()),
         minz = atof(bound[4]->GetValue().c_str()), maxz = atof(bound[5]->GetValue().c_str());
   float x, y, z;
   glBegin(GL_QUADS);
      x = minx; y = miny; z = minz; parent->mpWXCrystal->GetCrystal().FractionalToOrthonormalCoords(x,y,z); glVertex3f(x,y,z); 
      x = maxx; y = miny; z = minz; parent->mpWXCrystal->GetCrystal().FractionalToOrthonormalCoords(x,y,z); glVertex3f(x,y,z); 
      x = maxx; y = maxy; z = minz; parent->mpWXCrystal->GetCrystal().FractionalToOrthonormalCoords(x,y,z); glVertex3f(x,y,z); 
      x = minx; y = maxy; z = minz; parent->mpWXCrystal->GetCrystal().FractionalToOrthonormalCoords(x,y,z); glVertex3f(x,y,z);      

      x = minx; y = miny; z = maxz; parent->mpWXCrystal->GetCrystal().FractionalToOrthonormalCoords(x,y,z); glVertex3f(x,y,z); 
      x = maxx; y = miny; z = maxz; parent->mpWXCrystal->GetCrystal().FractionalToOrthonormalCoords(x,y,z); glVertex3f(x,y,z); 
      x = maxx; y = maxy; z = maxz; parent->mpWXCrystal->GetCrystal().FractionalToOrthonormalCoords(x,y,z); glVertex3f(x,y,z); 
      x = minx; y = maxy; z = maxz; parent->mpWXCrystal->GetCrystal().FractionalToOrthonormalCoords(x,y,z); glVertex3f(x,y,z);   


      x = minx; y = miny; z = minz; parent->mpWXCrystal->GetCrystal().FractionalToOrthonormalCoords(x,y,z); glVertex3f(x,y,z); 
      x = minx; y = miny; z = maxz; parent->mpWXCrystal->GetCrystal().FractionalToOrthonormalCoords(x,y,z); glVertex3f(x,y,z); 
      x = minx; y = maxy; z = maxz; parent->mpWXCrystal->GetCrystal().FractionalToOrthonormalCoords(x,y,z); glVertex3f(x,y,z); 
      x = minx; y = maxy; z = minz; parent->mpWXCrystal->GetCrystal().FractionalToOrthonormalCoords(x,y,z); glVertex3f(x,y,z); 

      x = maxx; y = miny; z = minz; parent->mpWXCrystal->GetCrystal().FractionalToOrthonormalCoords(x,y,z); glVertex3f(x,y,z); 
      x = maxx; y = miny; z = maxz; parent->mpWXCrystal->GetCrystal().FractionalToOrthonormalCoords(x,y,z); glVertex3f(x,y,z); 
      x = maxx; y = maxy; z = maxz; parent->mpWXCrystal->GetCrystal().FractionalToOrthonormalCoords(x,y,z); glVertex3f(x,y,z); 
      x = maxx; y = maxy; z = minz; parent->mpWXCrystal->GetCrystal().FractionalToOrthonormalCoords(x,y,z); glVertex3f(x,y,z); 


      x = minx; y = miny; z = minz; parent->mpWXCrystal->GetCrystal().FractionalToOrthonormalCoords(x,y,z); glVertex3f(x,y,z); 
      x = maxx; y = miny; z = minz; parent->mpWXCrystal->GetCrystal().FractionalToOrthonormalCoords(x,y,z); glVertex3f(x,y,z); 
      x = maxx; y = miny; z = maxz; parent->mpWXCrystal->GetCrystal().FractionalToOrthonormalCoords(x,y,z); glVertex3f(x,y,z); 
      x = minx; y = miny; z = maxz; parent->mpWXCrystal->GetCrystal().FractionalToOrthonormalCoords(x,y,z); glVertex3f(x,y,z); 

      x = minx; y = maxy; z = minz; parent->mpWXCrystal->GetCrystal().FractionalToOrthonormalCoords(x,y,z); glVertex3f(x,y,z); 
      x = maxx; y = maxy; z = minz; parent->mpWXCrystal->GetCrystal().FractionalToOrthonormalCoords(x,y,z); glVertex3f(x,y,z); 
      x = maxx; y = maxy; z = maxz; parent->mpWXCrystal->GetCrystal().FractionalToOrthonormalCoords(x,y,z); glVertex3f(x,y,z); 
      x = minx; y = maxy; z = maxz; parent->mpWXCrystal->GetCrystal().FractionalToOrthonormalCoords(x,y,z); glVertex3f(x,y,z); 
   glEnd();
}
}// namespace 

