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
#include <sstream> //for stringstream
#include <fstream>

#include <stdlib.h>

// wx headers, with or without precompilation
#include "wx/wxprec.h"
#ifdef __BORLANDC__
    #pragma hdrstop
#endif
#ifndef WX_PRECOMP
    #include "wx/wx.h"
#endif

#include "wxCryst/wxCrystal.h"

#include "wx/colordlg.h"
#include "wx/progdlg.h"
#include "wx/busyinfo.h"

#include "ObjCryst/Atom.h"
#include "ObjCryst/ZScatterer.h"
#include "ObjCryst/Molecule.h"
#include "ObjCryst/ScatteringPowerSphere.h"
#include "ObjCryst/Polyhedron.h"

#ifdef OBJCRYST_GL
   #ifdef __DARWIN__
      #include <OpenGL/glu.h>
   #else
      #include <GL/glu.h>
   #endif

   #ifdef __LINUX__
      #include "GL/glx.h"
      #ifdef HAVE_GLUT
         #include "GL/glut.h"
      #endif
   #endif
   #ifdef __WIN32__
     #include "gl/glaux.h"
   #endif
   #ifdef __DARWIN
     #include "AGL/agl.h"
   #endif
#endif

extern "C" {
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
#ifndef HAVE_GLUT
// This must be changed for each GL world to the correct first display list,
// i.e. in SetCurrent().
static int sFontDisplayListBase=0;
#endif

GLvoid crystGLPrint(const string &s)
{
   #ifdef HAVE_GLUT
   for(unsigned int l=0;l<s.size();l++)
      glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12,*(s.c_str()+l));
   #else
   glPushAttrib(GL_LIST_BIT);
      glListBase(sFontDisplayListBase - 32);
      glCallLists(s.size(), GL_UNSIGNED_BYTE, s.c_str());
   glPopAttrib();
   #endif
}

/// Conversion from ZScatterer to the newer Molecule object. (in WXZScatterer.cpp)
Molecule *ZScatterer2Molecule(ZScatterer *scatt);

// dialog to get a bounding box
  class UserSelectBoundingBox : public wxDialog {
  public:
    UserSelectBoundingBox(wxWindow * parent, char * title, BBox bbox);
    ~UserSelectBoundingBox ();
    BBox GetBBox ();
  private:
      void OnOk (void);
      wxTextCtrl * mpXminCtrl, *mpXmaxCtrl;
      wxTextCtrl * mpYminCtrl, *mpYmaxCtrl;
      wxTextCtrl * mpZminCtrl, *mpZmaxCtrl;
      BBox mbbox;
      DECLARE_EVENT_TABLE()
  };

// dialog to get a triplet of values box
  class UserXYZBox : public wxDialog {
  public:
    UserXYZBox(wxWindow * parent, char * title, Triple xyz);
    ~UserXYZBox ();
    Triple GetXYZ ();
  private:
      void OnOk (void);
      wxTextCtrl * mpXCtrl;
      wxTextCtrl * mpYCtrl;
      wxTextCtrl * mpZCtrl;
      Triple mXYZ;
      DECLARE_EVENT_TABLE()
  };


////////////////////////////////////////////////////////////////////////
//
//    WXCrystal
//
////////////////////////////////////////////////////////////////////////
static const long ID_CRYSTAL_MENU_SAVECIF                       =WXCRYST_ID();
static const long ID_CRYSTAL_MENU_SAVETEXT                      =WXCRYST_ID();
static const long ID_CRYSTAL_MENU_DISPLAY                       =WXCRYST_ID();
static const long ID_CRYSTAL_MENU_DISPLAY_3DVIEW                =WXCRYST_ID();
static const long ID_CRYSTAL_MENU_SCATT                         =WXCRYST_ID();
static const long ID_CRYSTAL_MENU_PAR_ADDANTIBUMP               =WXCRYST_ID();
static const long ID_CRYSTAL_MENU_PAR_SETRELATIVEXYZLIMITS      =WXCRYST_ID();
static const long ID_CRYSTAL_MENU_SCATT_REMOVESCATTPOW          =WXCRYST_ID();
static const long ID_CRYSTAL_MENU_SCATT_ADDSCATTPOWATOM         =WXCRYST_ID();
static const long ID_CRYSTAL_MENU_SCATT_ADDSCATTPOWSPHERE       =WXCRYST_ID();
static const long ID_CRYSTAL_MENU_SCATT_ADDATOM                 =WXCRYST_ID();
static const long ID_CRYSTAL_MENU_SCATT_ADDZSCATTERER           =WXCRYST_ID();
static const long ID_CRYSTAL_MENU_SCATT_ADDMOLECULE             =WXCRYST_ID();
static const long ID_CRYSTAL_MENU_SCATT_IMPORTFENSKEHALLZMATRIX =WXCRYST_ID();
static const long ID_CRYSTAL_MENU_SCATT_ADDTETRAHEDRON          =WXCRYST_ID();
static const long ID_CRYSTAL_MENU_SCATT_ADDOCTAHEDRON           =WXCRYST_ID();
static const long ID_CRYSTAL_MENU_SCATT_ADDTRIANGLE             =WXCRYST_ID();
static const long ID_CRYSTAL_MENU_SCATT_ADDSQUAREPLANE          =WXCRYST_ID();
static const long ID_CRYSTAL_MENU_SCATT_ADDCUBE                 =WXCRYST_ID();
static const long ID_CRYSTAL_MENU_SCATT_ADDANTIPRISMTETRAGONAL  =WXCRYST_ID();
static const long ID_CRYSTAL_MENU_SCATT_ADDPRISMTRIGONAL        =WXCRYST_ID();
static const long ID_CRYSTAL_MENU_SCATT_ADDICOSAHEDRON          =WXCRYST_ID();
static const long ID_CRYSTAL_MENU_SCATT_REMOVESCATTERER         =WXCRYST_ID();
static const long ID_CRYSTAL_MENU_SCATT_DUPLICSCATTERER         =WXCRYST_ID();
static const long ID_CRYSTAL_SPACEGROUP                         =WXCRYST_ID();
static const long ID_GLCRYSTAL_MENU_UPDATE                      =WXCRYST_ID();

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
   EVT_MENU(ID_CRYSTAL_MENU_SCATT_ADDMOLECULE,        WXCrystal::OnMenuAddScatterer)
   EVT_MENU(ID_CRYSTAL_MENU_SCATT_IMPORTFENSKEHALLZMATRIX,WXCrystal::OnMenuImportMoleculeFromFenskeHallZMatrix)
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
,mCrystalGLDisplayList(0),mCrystalGLNameDisplayList(0),
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
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_CRYSTAL_MENU_SAVETEXT,"Save as text");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_CRYSTAL_MENU_SAVECIF,"Save as CIF");
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
         mpMenuBar->GetMenu(ID_CRYSTAL_MENU_SCATT).AppendSeparator();
         mpMenuBar->AddMenuItem(ID_CRYSTAL_MENU_SCATT,ID_CRYSTAL_MENU_SCATT_ADDATOM,
                                "Add Atom");
         mpMenuBar->AddMenuItem(ID_CRYSTAL_MENU_SCATT,ID_CRYSTAL_MENU_SCATT_ADDMOLECULE,
                                "Add Molecule");
         mpMenuBar->AddMenuItem(ID_CRYSTAL_MENU_SCATT,ID_CRYSTAL_MENU_SCATT_IMPORTFENSKEHALLZMATRIX,
                                "Import Molecule from Fenske-Hall Z-Matrix");
         mpMenuBar->AddMenuItem(ID_CRYSTAL_MENU_SCATT,ID_CRYSTAL_MENU_SCATT_REMOVESCATTERER,
                                "Remove Scatterer");
         mpMenuBar->AddMenuItem(ID_CRYSTAL_MENU_SCATT,ID_CRYSTAL_MENU_SCATT_DUPLICSCATTERER,
                                "Duplicate Scatterer");
         mpMenuBar->GetMenu(ID_CRYSTAL_MENU_SCATT).AppendSeparator();
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
      mpMenuBar->AddMenu("Display",ID_CRYSTAL_MENU_DISPLAY);
         mpMenuBar->AddMenuItem(ID_CRYSTAL_MENU_DISPLAY,ID_CRYSTAL_MENU_DISPLAY_3DVIEW,
                                "3D Display");
   // AntiBump-ProMerge cost
      wxBoxSizer* pAntiBumpSizer=new wxBoxSizer(wxHORIZONTAL);
      WXFieldPar<REAL> *pWXFieldBumpMerge=
         new WXFieldPar<REAL>(this,"AntiBump",-1,&(mpCrystal->mBumpMergeCost),100);
      WXFieldPar<REAL> *pAntiBumpScale=
         new WXFieldPar<REAL>(this,"Scale",-1,&(mpCrystal->mBumpMergeScale));
      pAntiBumpSizer->Add(pWXFieldBumpMerge);
      pAntiBumpSizer->Add(pAntiBumpScale);
      mpSizer->Add(pAntiBumpSizer,0,wxALIGN_LEFT);
      mList.Add(pWXFieldBumpMerge);
      mList.Add(pAntiBumpScale);
   // Lattice
      wxBoxSizer* lattice=new wxBoxSizer(wxHORIZONTAL);
#if 1
      WXFieldRefPar* pFieldLatticeA    =new WXFieldRefPar(this,"a:",
                                     &(mpCrystal->GetPar("a")) );
          
      WXFieldRefPar* pFieldLatticeB    =new WXFieldRefPar(this,"b:",
                                     &(mpCrystal->GetPar("b")) );
          
      WXFieldRefPar* pFieldLatticeC    =new WXFieldRefPar(this,"c:",
                                     &(mpCrystal->GetPar("c")) );
          
      WXFieldRefPar* pFieldLatticeAlpha=new WXFieldRefPar(this,"alpha:",
                                     &(mpCrystal->GetPar("alpha")) );
          
      WXFieldRefPar* pFieldLatticeBeta =new WXFieldRefPar(this,"beta:",
                                     &(mpCrystal->GetPar("beta")) );
          
      WXFieldRefPar* pFieldLatticeGamma=new WXFieldRefPar(this,"gamma:",
                                     &(mpCrystal->GetPar("gamma")) );
#else
      WXCrystObjBasic* pFieldLatticeA 
         =mpCrystal->GetPar("a").WXCreate(this);
      WXCrystObjBasic* pFieldLatticeB 
         =mpCrystal->GetPar("b").WXCreate(this);
      WXCrystObjBasic* pFieldLatticeC 
         =mpCrystal->GetPar("c").WXCreate(this);
      WXCrystObjBasic* pFieldLatticeAlpha 
         =mpCrystal->GetPar("alpha").WXCreate(this);
      WXCrystObjBasic* pFieldLatticeBeta 
         =mpCrystal->GetPar("beta").WXCreate(this);
      WXCrystObjBasic* pFieldLatticeGamma
         =mpCrystal->GetPar("gamma").WXCreate(this);
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
   mpCrystal->GetBumpMergeCost();
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
      if(mCrystalGLDisplayList==0)
      {
         mCrystalGLDisplayList=glGenLists(1);
         mCrystalGLNameDisplayList=glGenLists(1);
         VFN_DEBUG_MESSAGE("WXCrystal::UpdateGL():created mCrystalGLDisplayList="<<mCrystalGLDisplayList,7)
      }
      
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
      glNewList(mCrystalGLNameDisplayList,GL_COMPILE);
         glPushMatrix();
            mpCrystal->GLInitDisplayList(onlyIndependentAtoms,xMin,xMax,yMin,yMax,zMin,zMax,true);
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

int WXCrystal::GrabCrystalGLDisplayList(const bool atomName)const
{
   VFN_DEBUG_MESSAGE("WXCrystal::GrabCrystalGLDisplayList()",7)
   //:KLUDGE: ? or OK ?
   while(mCrystalGLDisplayListIsLocked) wxUsleep(5);
   mCrystalGLDisplayListIsLocked=true;
   VFN_DEBUG_MESSAGE("WXCrystal::GrabCrystalGLDisplayList():"<<mCrystalGLDisplayList,7)
   if(atomName) return mCrystalGLNameDisplayList;
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
   #if wxUSE_STATUSBAR
   frame->CreateStatusBar(1);
   frame->SetStatusText(mpCrystal->GetName().c_str());
   #endif
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
   if(event.GetId()== ID_CRYSTAL_MENU_SCATT_ADDATOM)
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
   }
   if(event.GetId()== ID_CRYSTAL_MENU_SCATT_ADDZSCATTERER)
   {
      scatt=new ZScatterer("Change Me!",*mpCrystal);
   }
   if(event.GetId()== ID_CRYSTAL_MENU_SCATT_ADDMOLECULE)
   {
      scatt=new Molecule(*mpCrystal,"Molecule");
   }
   if(event.GetId()== ID_CRYSTAL_MENU_SCATT_ADDTETRAHEDRON)
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
            VFN_DEBUG_EXIT("WXZScatterer))OnMenuAddZAtom())Cancelled",6)
            return;
         }
         double bondLength;
         bondLengthDialog.GetValue().ToDouble(&bondLength);

      Molecule *mol=MakeTetrahedron(*mpCrystal,scattPow1->GetName()+scattPow2->GetName()+"4",
                                           scattPow1,scattPow2,bondLength);
      mol->RestraintStatus(cout);
      scatt=mol;
   }
   if(event.GetId()== ID_CRYSTAL_MENU_SCATT_ADDOCTAHEDRON)
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
                                    this,"Corner atom type (ScatteringPower))",choice);
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
            VFN_DEBUG_EXIT("WXZScatterer))OnMenuAddZAtom())Cancelled",6)
            return;
         }
         double bondLength;
         bondLengthDialog.GetValue().ToDouble(&bondLength);

      Molecule *mol=MakeOctahedron(*mpCrystal,scattPow1->GetName()+scattPow2->GetName()+"6",
                                   scattPow1,scattPow2,bondLength);
      mol->RestraintStatus(cout);
      scatt=mol;
   }
   if(event.GetId()== ID_CRYSTAL_MENU_SCATT_ADDTRIANGLE)
   {
      VFN_DEBUG_MESSAGE("WXCrystal::OnMenuAddScatterer())Add triangle plane",6)
      int choice;
      //Scattering power 1
         const ScatteringPower *scattPow1=WXDialogChooseFromRegistry(
                                    mpCrystal->GetScatteringPowerRegistry(),
                                    this,"Central atom type (ScatteringPower))",choice);
         if(0==scattPow1)
         {
            VFN_DEBUG_EXIT("WXCrystal::OnMenuAddScatterer():Canceled",6)
            return;
         }
      //Scattering power 2
         const ScatteringPower *scattPow2=WXDialogChooseFromRegistry(
                                    mpCrystal->GetScatteringPowerRegistry(),
                                    this,"Corner atom type (ScatteringPower))",choice);
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

      Molecule *mol=MakeTriangle(*mpCrystal,scattPow1->GetName()+scattPow2->GetName()+"3",
                                 scattPow1,scattPow2,bondLength);
      mol->RestraintStatus(cout);
      scatt=mol;
   }
   if(event.GetId()== ID_CRYSTAL_MENU_SCATT_ADDSQUAREPLANE)
   {
      VFN_DEBUG_MESSAGE("WXCrystal::OnMenuAddScatterer():Add square plane",6)
      int choice;
      //Scattering power 1
         const ScatteringPower *scattPow1=WXDialogChooseFromRegistry(
                                    mpCrystal->GetScatteringPowerRegistry(),
                                    this,"Central atom type (ScatteringPower))",choice);
         if(0==scattPow1)
         {
            VFN_DEBUG_EXIT("WXCrystal::OnMenuAddScatterer():Canceled",6)
            return;
         }
      //Scattering power 2
         const ScatteringPower *scattPow2=WXDialogChooseFromRegistry(
                                    mpCrystal->GetScatteringPowerRegistry(),
                                    this,"Corner atom type (ScatteringPower))",choice);
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

      Molecule *mol=MakeSquarePlane(*mpCrystal,scattPow1->GetName()+scattPow2->GetName()+"4",
                                    scattPow1,scattPow2,bondLength);
      mol->RestraintStatus(cout);
      scatt=mol;
   }
   if(event.GetId()== ID_CRYSTAL_MENU_SCATT_ADDCUBE)
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
                                    this,"Corner atom type (ScatteringPower))",choice);
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

      Molecule *mol=MakeCube(*mpCrystal,scattPow1->GetName()+scattPow2->GetName()+"8",
                             scattPow1,scattPow2,bondLength);
      mol->RestraintStatus(cout);
      scatt=mol;
   }
   if(event.GetId()== ID_CRYSTAL_MENU_SCATT_ADDANTIPRISMTETRAGONAL)
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
                                    this,"Corner atom type (ScatteringPower))",choice);
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

      Molecule *mol=MakeAntiPrismTetragonal(*mpCrystal,scattPow1->GetName()+scattPow2->GetName()+"8",
                                            scattPow1,scattPow2,bondLength);
      mol->RestraintStatus(cout);
      scatt=mol;
   }
   if(event.GetId()== ID_CRYSTAL_MENU_SCATT_ADDPRISMTRIGONAL)
   {
      int choice;
      //Scattering power 1
         const ScatteringPower *scattPow1=WXDialogChooseFromRegistry(
                                    mpCrystal->GetScatteringPowerRegistry(),
                                    this,"Central atom type (ScatteringPower))",choice);
         if(0==scattPow1)
         {
            VFN_DEBUG_EXIT("WXCrystal::OnMenuAddScatterer():Canceled",6)
            return;
         }
      //Scattering power 2
         const ScatteringPower *scattPow2=WXDialogChooseFromRegistry(
                                    mpCrystal->GetScatteringPowerRegistry(),
                                    this,"Corner atom type (ScatteringPower))",choice);
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
            VFN_DEBUG_EXIT("WXZScatterer::OnMenuAddZAtom())Cancelled",6)
            return;
         }
         double bondLength;
         bondLengthDialog.GetValue().ToDouble(&bondLength);

      Molecule *mol=MakePrismTrigonal(*mpCrystal,scattPow1->GetName()+scattPow2->GetName()+"6",
                                      scattPow1,scattPow2,bondLength);
      mol->RestraintStatus(cout);
      scatt=mol;
   }
   if(event.GetId()== ID_CRYSTAL_MENU_SCATT_ADDICOSAHEDRON)
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

      Molecule *mol=MakeIcosahedron(*mpCrystal,scattPow1->GetName()+scattPow2->GetName()+"12",
                                    scattPow1,scattPow2,bondLength);
      mol->RestraintStatus(cout);
      scatt=mol;
   }
   mpCrystal->AddScatterer(scatt);
   //mpCrystal->XMLOutput(cout);
   VFN_DEBUG_MESSAGE("WXCrystal::OnMenuAddScatterer():calling Layout()",6)
   //this->CrystUpdate();
   this->Layout();
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
   this->CrystUpdate();
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

Molecule *ZScatterer2Molecule(ZScatterer *scatt);//defined in wxZScatterer.cpp

void WXCrystal::OnMenuImportMoleculeFromFenskeHallZMatrix(wxCommandEvent &event)
{
   VFN_DEBUG_ENTRY("WXCrystal::OnMenuImportFenskeHallZMatrix()",6)
   wxFileDialog open(this,"Choose a file with a Fenske-Hall Z-matrix","","","*.fhz",
                     wxOPEN | wxFILE_MUST_EXIST);
   if(open.ShowModal() != wxID_OK) return;
   ifstream fin (open.GetPath().c_str());
   if(!fin)
   {
      throw ObjCrystException("WXCrystal::OnMenuImportFenskeHallZMatrix() : \
Error opening file for input:"+string(open.GetPath().c_str()));
   }
   string filename=open.GetPath().c_str();
   string shortName;
   {// Use short name
      std::string::size_type idx =filename.rfind("/");
      std::string::size_type idx2=filename.rfind("\\");
      std::string::size_type idx3=filename.rfind(":");
      if(((long)idx2!=(long)string::npos)&&((long)idx2>(long)idx))idx=idx2;
      if(((long)idx3!=(long)string::npos)&&((long)idx3>(long)idx))idx=idx3;
      if(idx==string::npos)
         shortName=filename;
      else
         shortName=filename.substr(idx+1);
   }
   ZScatterer scatt(shortName,*mpCrystal);
   scatt.ImportFenskeHallZMatrix(fin);
   fin.close();
   mpCrystal->AddScatterer(ZScatterer2Molecule(&scatt));
   this->CrystUpdate();
   VFN_DEBUG_EXIT("WXCrystal::OnMenuImportFenskeHallZMatrix()",6)
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
   mpCrystal->UpdateDisplay();
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
   if(id==ID_CRYSTAL_SPACEGROUP)
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
      this->CrystUpdate();
      this->Layout();
      return true;
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
					  WXGLCrystalCanvas * parentCrystal) const
{
   VFN_DEBUG_ENTRY("UnitCellMapImport::GLInitDisplayList()",7)
   // Generate triangles
      VFN_DEBUG_MESSAGE("UnitCellMapImport::GLInitDisplayList(): Generate Triangles",7)

      const int nx=mPoints.cols();
      const int ny=mPoints.rows();
      const int nz=mPoints.depth();
      float step[3];
      step[0]=1/(float)nx;
      step[1]=1/(float)ny;
      step[2]=1/(float)nz;
      int nxMin, nxMax, nyMin, nyMax, nzMin, nzMax;
      BBox mapbbox = parentCrystal->GetMapBBox();
      // use cell bbox if mapbbox has zero volume (default)
      if (mapbbox.xMin == mapbbox.xMax) mapbbox = parentCrystal->GetCellBBox();
      nxMin = (int)(mapbbox.xMin * nx);
      nxMax = (int)(mapbbox.xMax * nx);
      nyMin = (int)(mapbbox.yMin * ny);
      nyMax = (int)(mapbbox.yMax * ny);
      nzMin = (int)(mapbbox.zMin * nz);
      nzMax = (int)(mapbbox.zMax * nz);
      const int snx = nxMax-nxMin+1, sny = nyMax-nyMin+1, snz = nzMax-nzMin+1;
      const unsigned int sny_snz = sny*snz;
      int i, j, k;
      unsigned int ni, nj, si, sj, sk, sni, snj, sind;
      float x, y, z;

      //create new set of points
      mp4Vector * subPoints = new mp4Vector[snx*sny*snz];
      for(i=nxMin, si=0; i <= nxMax; i++, si++)
      {
         ni = ((nx + i % nx) % nx);    //this will 'wrap' around any value (negative or positive)
         sni = si*sny_snz;
         for(j=nyMin, sj=0; j <= nyMax; j++, sj++)
         {
            nj = ((ny + j % ny) % ny);
            snj = sj*snz;
            for(k=nzMin, sk=0; k <= nzMax; k++, sk++)
            {
               sind = sni + snj + sk;
               x = i*step[0]; y = j*step[1]; z = k*step[2];
               mpCrystal->FractionalToOrthonormalCoords(x, y, z);
               subPoints[sind].x = x; subPoints[sind].y = y; subPoints[sind].z = z;
               //cout << ni <<" "<<nj<<" "<<(nz+ k % nz)<<endl;
               subPoints[sind].val = mPoints((nz+ k % nz)% nz,nj,ni);
            }
         }
      }
      int numOfTriangles;
      VFN_DEBUG_MESSAGE("UnitCellMapImport::GLInitDisplayList(): MC, Min Value="<<minValue,10)
      const TRIANGLE *pTriangles= MC(snx-1, sny-1, snz-1, step[0], step[1], step[2], minValue, subPoints, numOfTriangles);
   // OpenGL drawing instructions
      VFN_DEBUG_MESSAGE("UnitCellMapImport::GLInitDisplayList(): OpenGL instructions",7)
      glBegin(GL_TRIANGLES);
         float normx,normy,normz;
         for(int i=0; i < numOfTriangles; i++)
         {
            for(int j=0; j < 3; j++)
            {
               //VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnPaint():MC1:"<<i<<" "<<j,5)
               //:TODO: Fix normals
               normx=pTriangles[i].norm[j].x;
               normy=pTriangles[i].norm[j].y;
               normz=pTriangles[i].norm[j].z;
               //mpCrystal->FractionalToOrthonormalCoords(normx, normy, normz);
               //mpCrystal->OrthonormalToFractionalCoords(normx, normy, normz);
               glNormal3f(normx, normy, normz);
               glVertex3f(pTriangles[i].p[j].x    ,pTriangles[i].p[j].y    ,pTriangles[i].p[j].z);
            }
         }
      glEnd();

   delete [] subPoints; 
   VFN_DEBUG_EXIT("UnitCellMapImport::GLInitDisplayList()",7)

}
int UnitCellMapImport::ImportGRD(const string&filename)
{
   VFN_DEBUG_ENTRY("UnitCellMapImport::ImportGRD()",7)
   ifstream ffile(filename.c_str());
   if(!ffile.is_open())
   {     //if file could not be loaded for some reason then exit
     VFN_DEBUG_MESSAGE("UnitCellMapImport::ImportGRD() error opening "<<filename.c_str(),10)
      (*fpObjCrystInformUser)("Error opening file: "+filename);
      return 0;
   }
   //message for reporting errors
   char buff[99];
   ffile.getline(buff, 100);
   float a, b, c, alpha, beta, gamma;
   ffile >>a >>b >>c >>alpha >>beta >>gamma;
   if(!ffile.good()) {  (*fpObjCrystInformUser)("Error reading file: "+filename); return 0; }
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
   if(!ffile.good()) {  (*fpObjCrystInformUser)("Error reading file: "+filename); return 0; }
   mPoints.resize(nz,ny,nx);
   for(int i=0; i < nx; i++) {
     for(int j=0; j < ny; j++) {
        for(int k=0; k < nz; k++) {
           ffile >>mPoints(k,j,i);      //reading rhos
        }
     }
   }
   ffile.close();
   {// Use short name
      std::string::size_type idx =filename.rfind("/");
      std::string::size_type idx2=filename.rfind("\\");
      std::string::size_type idx3=filename.rfind(":");
      if(((long)idx2!=(long)string::npos)&&((long)idx2>(long)idx))idx=idx2;
      if(((long)idx3!=(long)string::npos)&&((long)idx3>(long)idx))idx=idx3;
      if(idx==string::npos)
         mName=filename;
      else
      {
         cout<<"name="<<filename.substr(idx+1)<<endl;
         mName=filename.substr(idx+1);
      }
   }
   VFN_DEBUG_EXIT("UnitCellMapImport::ImportGRD()",7)
     return 1;
}

const string & UnitCellMapImport::GetName()const
{
   return mName;
}

REAL UnitCellMapImport::GetValue(const REAL x,const REAL y,const REAL z)const
{
   const int nx=mPoints.cols();
   const int ny=mPoints.rows();
   const int nz=mPoints.depth();
   long ix=((long)(x*nx))%nx;
   long iy=((long)(y*ny))%ny;
   long iz=((long)(z*nz))%nz;
   if(ix<0) ix+=nx;
   if(iy<0) iy+=ny;
   if(iz<0) iz+=nz;
   return mPoints(iz,iy,ix);
}

////////////////////////////////////////////////////////////////////////
//
//    UnitCellMapGLList
//
////////////////////////////////////////////////////////////////////////
UnitCellMapGLList::UnitCellMapGLList(const bool showWire,
                                     const float r,const float g,const float b,
                                     const float t):
mGLDisplayList(0),mShowWire(showWire)
{
   VFN_DEBUG_MESSAGE("UnitCellMapGLList::UnitCellMapGLList()",10)
   this->SetColour(r,g,b,t);
}

UnitCellMapGLList::~UnitCellMapGLList()
{
   VFN_DEBUG_MESSAGE("UnitCellMapGLList::~UnitCellMapGLList()",10)
   if(0!=mGLDisplayList) glDeleteLists(mGLDisplayList,1);
}
void UnitCellMapGLList::GenList(const UnitCellMapImport &ucmap,
				WXGLCrystalCanvas * parent,
                                const float contourValue)
{
   VFN_DEBUG_ENTRY("UnitCellMapGLList::GenList()",7)
   if(0==mGLDisplayList) mGLDisplayList=glGenLists(1);
   glNewList(mGLDisplayList,GL_COMPILE);
      glPushMatrix();
      ucmap.GLInitDisplayList(contourValue, parent);
      glPopMatrix();
   glEndList();
   VFN_DEBUG_EXIT("UnitCellMapGLList::GenList()",7)
}

void UnitCellMapGLList::SetColour(const float r,const float g,const float b,
                                 const float t)
{
   mColour[0]=r;
   mColour[1]=g;
   mColour[2]=b;
   mColour[3]=t;
}

const float* UnitCellMapGLList::GetColour()const
{
   return mColour;
}

void UnitCellMapGLList::ToggleShowWire()
{
   mShowWire =! mShowWire;
}

void UnitCellMapGLList::Draw()const
{
   if(0==mGLDisplayList)
   {
      VFN_DEBUG_MESSAGE("UnitCellMapGLList::Draw():No Display list generated !",7)
      return;
   }
   glPushMatrix();
      if(mShowWire) glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
      else glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
      
      glMaterialfv(GL_FRONT, GL_AMBIENT, mColour);
      glMaterialfv(GL_FRONT, GL_DIFFUSE, mColour);
      glMaterialfv(GL_FRONT, GL_SPECULAR, mColour);
      const GLfloat colour0[] = {0.0f, 0.0f, 0.0f, 0.0f}; 
      glMaterialfv(GL_FRONT, GL_EMISSION, colour0); 
      // :TODO: 
      // Disabled Shininess as there is a problem with normals 
      // and non-orthogonal unit cells
      glMaterialf( GL_FRONT, GL_SHININESS, 0.0); 

      const GLfloat colorBack [] = {mColour[0]/3., mColour[1]/3., mColour[2]/3., 1.00}; 
      glMaterialfv(GL_BACK, GL_AMBIENT, colorBack);
      glMaterialfv(GL_BACK, GL_DIFFUSE, colorBack);
      glMaterialfv(GL_BACK, GL_SPECULAR, colorBack);
      glMaterialf( GL_BACK, GL_SHININESS, 0.0); 
      // :TODO: Check display list is not being modified (lock it), useless for now
      // as the map is not dynamically updated.
      glCallList(mGLDisplayList);
   glPopMatrix();
   VFN_DEBUG_EXIT("UnitCellMapGLList::Draw()",7)
}

void UnitCellMapGLList::SetName(const string &name)
{
   mName=name;
}
const string &UnitCellMapGLList::GetName()const
{
   return mName;
}

////////////////////////////////////////////////////////////////////////
//
//    WXGLCrystalCanvas
//
////////////////////////////////////////////////////////////////////////
static const long ID_GLCRYSTAL_MENU_SHOWATOMLABEL=     WXCRYST_ID(); 
static const long ID_GLCRYSTAL_MENU_SHOWCURSOR=        WXCRYST_ID(); 
static const long ID_GLCRYSTAL_MENU_SETCURSOR=        WXCRYST_ID(); 
static const long ID_GLCRYSTAL_UPDATEUI=               WXCRYST_ID(); 
static const long ID_GLCRYSTAL_MENU_CHANGELIMITS=      WXCRYST_ID(); 
static const long ID_GLCRYSTAL_MENU_SHOWCRYSTAL=       WXCRYST_ID(); 
static const long ID_GLCRYSTAL_MENU_LOADFOURIER=       WXCRYST_ID(); 
static const long ID_GLCRYSTAL_MENU_CHANGECONTOUR=     WXCRYST_ID(); 
static const long ID_GLCRYSTAL_MENU_ADDCONTOUR=        WXCRYST_ID(); 
static const long ID_GLCRYSTAL_MENU_SHOWFOURIER=       WXCRYST_ID(); 
static const long ID_GLCRYSTAL_MENU_FOURIERCHANGECOLOR=WXCRYST_ID(); 
static const long ID_GLCRYSTAL_MENU_SHOWWIRE=          WXCRYST_ID(); 
static const long ID_GLCRYSTAL_MENU_UNLOADFOURIER=     WXCRYST_ID(); 
static const long ID_GLCRYSTAL_MENU_FOURIERCHANGEBBOX= WXCRYST_ID(); 

BEGIN_EVENT_TABLE(WXGLCrystalCanvas, wxGLCanvas)
   EVT_SIZE             (WXGLCrystalCanvas::OnSize)
   EVT_PAINT            (WXGLCrystalCanvas::OnPaint)
   EVT_ERASE_BACKGROUND (WXGLCrystalCanvas::OnEraseBackground)
   EVT_MOUSE_EVENTS     (WXGLCrystalCanvas::OnMouse)
   EVT_MENU             (ID_GLCRYSTAL_MENU_UPDATE,              WXGLCrystalCanvas::OnUpdate)
   EVT_MENU             (ID_GLCRYSTAL_MENU_CHANGELIMITS,        WXGLCrystalCanvas::OnChangeLimits)
   EVT_MENU             (ID_GLCRYSTAL_MENU_SHOWCRYSTAL,         WXGLCrystalCanvas::OnShowCrystal)
   EVT_MENU             (ID_GLCRYSTAL_MENU_SHOWATOMLABEL,       WXGLCrystalCanvas::OnShowAtomLabel)
   EVT_MENU             (ID_GLCRYSTAL_MENU_SHOWCURSOR,          WXGLCrystalCanvas::OnShowCursor)
   EVT_MENU             (ID_GLCRYSTAL_MENU_SETCURSOR,           WXGLCrystalCanvas::OnSetCursor)
   EVT_MENU             (ID_GLCRYSTAL_MENU_LOADFOURIER,         WXGLCrystalCanvas::OnLoadFourier)
   EVT_MENU             (ID_GLCRYSTAL_MENU_CHANGECONTOUR,       WXGLCrystalCanvas::OnChangeContour)
   EVT_MENU             (ID_GLCRYSTAL_MENU_ADDCONTOUR,          WXGLCrystalCanvas::OnAddContour)
   EVT_MENU             (ID_GLCRYSTAL_MENU_SHOWFOURIER,         WXGLCrystalCanvas::OnShowFourier)
   EVT_MENU             (ID_GLCRYSTAL_MENU_FOURIERCHANGECOLOR,  WXGLCrystalCanvas::OnFourierChangeColor)
   EVT_MENU             (ID_GLCRYSTAL_MENU_SHOWWIRE,            WXGLCrystalCanvas::OnShowWire)
   EVT_MENU             (ID_GLCRYSTAL_MENU_UNLOADFOURIER,       WXGLCrystalCanvas::OnUnloadFourier)
   EVT_MENU             (ID_GLCRYSTAL_MENU_FOURIERCHANGEBBOX,   WXGLCrystalCanvas::OnFourierChangeBbox)
   EVT_CHAR             (WXGLCrystalCanvas::OnKeyDown)
   EVT_KEY_DOWN         (WXGLCrystalCanvas::OnKeyDown)
   EVT_KEY_UP           (WXGLCrystalCanvas::OnKeyUp)
   EVT_UPDATE_UI(ID_GLCRYSTAL_UPDATEUI,WXGLCrystalCanvas::OnUpdateUI)
END_EVENT_TABLE()

WXGLCrystalCanvas::WXGLCrystalCanvas(WXCrystal *wxcryst,
                                     wxFrame *parent, wxWindowID id,
                                     const wxPoint &pos,
                                     const wxSize &size):
wxGLCanvas(parent,id,pos,size,wxDEFAULT_FRAME_STYLE),mpParentFrame(parent),
mpWXCrystal(wxcryst),mIsGLInit(false),mDist(60),mX0(0),mY0(0),mZ0(0),mViewAngle(15),
mShowFourier(true),mShowCrystal(true),mShowAtomName(true),mShowCursor(false),
mIsGLFontBuilt(false),mGLFontDisplayListBase(0)
{
   VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::WXGLCrystalCanvas()",3)
     // N.B. xMin=xMax so that the previous cell bbox is used for Maps 
     // until mmapbbox is changed
   mcellbbox.xMin = mcellbbox.yMin = mcellbbox.zMin = -0.1;
   mcellbbox.xMax = mcellbbox.yMax = mcellbbox.zMax = 1.1;
   mmapbbox.xMin = mmapbbox.xMax = mmapbbox.yMin = mmapbbox.zMin = 0.;
   mmapbbox.yMax = mmapbbox.zMax = 1.;
   mpPopUpMenu=new wxMenu("Crystal");
   mpPopUpMenu->Append(ID_GLCRYSTAL_MENU_UPDATE, "&Update");
   mpPopUpMenu->Append(ID_GLCRYSTAL_MENU_CHANGELIMITS, "Change display &Limits");
   
   mpPopUpMenu->Append(ID_GLCRYSTAL_MENU_SHOWCRYSTAL, "Hide Crystal");
   mpPopUpMenu->Append(ID_GLCRYSTAL_MENU_SHOWATOMLABEL, "Hide Atom Labels");
   mpPopUpMenu->Append(ID_GLCRYSTAL_MENU_SHOWCURSOR, "Show Cursor");
   mpPopUpMenu->Append(ID_GLCRYSTAL_MENU_SETCURSOR, "Set view cntr and cursor pos.");
   mpPopUpMenu->AppendSeparator();
   mpPopUpMenu->Append(ID_GLCRYSTAL_MENU_LOADFOURIER, "Load Fourier Map");	
   mpPopUpMenu->Append(ID_GLCRYSTAL_MENU_CHANGECONTOUR, "Change Contour Value");
   mpPopUpMenu->Enable(ID_GLCRYSTAL_MENU_CHANGECONTOUR, FALSE);	//disable it for now
   mpPopUpMenu->Append(ID_GLCRYSTAL_MENU_ADDCONTOUR, "Add Contour Value");
   mpPopUpMenu->Enable(ID_GLCRYSTAL_MENU_ADDCONTOUR, FALSE);	//disable it for now
   mpPopUpMenu->Append(ID_GLCRYSTAL_MENU_SHOWFOURIER, "Hide Fourier Map");
   mpPopUpMenu->Enable(ID_GLCRYSTAL_MENU_SHOWFOURIER, FALSE);
   mpPopUpMenu->Append(ID_GLCRYSTAL_MENU_FOURIERCHANGECOLOR, "Change Fourier Color");
   mpPopUpMenu->Enable(ID_GLCRYSTAL_MENU_FOURIERCHANGECOLOR, FALSE);
   mpPopUpMenu->Append(ID_GLCRYSTAL_MENU_SHOWWIRE, "Toggle WireFrame/Filled");
   mpPopUpMenu->Enable(ID_GLCRYSTAL_MENU_SHOWWIRE, FALSE);
   mpPopUpMenu->Append(ID_GLCRYSTAL_MENU_UNLOADFOURIER, "Unload Fourier Map(s)");
   mpPopUpMenu->Enable(ID_GLCRYSTAL_MENU_UNLOADFOURIER, FALSE);
   mpPopUpMenu->Append(ID_GLCRYSTAL_MENU_FOURIERCHANGEBBOX, "Change Fourier Limits");
   mpPopUpMenu->Enable(ID_GLCRYSTAL_MENU_FOURIERCHANGEBBOX, FALSE);	//disable it for now
}

WXGLCrystalCanvas::~WXGLCrystalCanvas()
{
   mpWXCrystal->NotifyCrystalGLDelete();
   {
      vector<UnitCellMapImport*>::iterator pos;
      for(pos=mvpUnitCellMapImport.begin();pos != mvpUnitCellMapImport.end();pos++)
         delete *pos;
      mvpUnitCellMapImport.clear();
   }
   {
      vector<pair<pair<const UnitCellMapImport*,float>,UnitCellMapGLList* > >::iterator pos;
      for(pos=mvpUnitCellMapGLList.begin();pos != mvpUnitCellMapGLList.end();pos++)
         delete pos->second;
      mvpUnitCellMapGLList.clear();
   }
   #ifndef HAVE_GLUT
   this->DeleteGLFont();
   #endif
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
   
   if(mShowFourier)
   {
      glPushMatrix();
         // The display origin is the center of the Crystal BoundingBox, so translate
            BBox cellbbox = this->GetCellBBox();
            REAL xc=(cellbbox.xMin+cellbbox.xMax)/2.;  
            REAL yc=(cellbbox.yMin+cellbbox.yMax)/2.; 
            REAL zc=(cellbbox.zMin+cellbbox.zMax)/2.; 
            mpWXCrystal->GetCrystal().FractionalToOrthonormalCoords(xc, yc, zc);
            glTranslatef(-xc, -yc, -zc);
         // Draw all Fourier maps
         vector<pair<pair<const UnitCellMapImport*,float>,UnitCellMapGLList* > >::
            const_iterator pos;
         for(pos=mvpUnitCellMapGLList.begin();pos != mvpUnitCellMapGLList.end();++pos)
            pos->second->Draw();
      glPopMatrix();
   }
   if(mShowCrystal)
   {
      glCallList(mpWXCrystal->GrabCrystalGLDisplayList());  //Draw Crystal
      mpWXCrystal->ReleaseCrystalGLDisplayList();
      if(mShowAtomName)
      {
         glLoadIdentity();
         glTranslatef( -0.3, 0, -mDist+1. );// Put labels in front of the atom position
         glMultMatrixf( &m[0][0] );
         glTranslatef( mX0, mY0, mZ0 );
         glCallList(mpWXCrystal->GrabCrystalGLDisplayList(true));  //Draw Atom Names
         mpWXCrystal->ReleaseCrystalGLDisplayList();
      }
         
   }
   if(mShowCursor)
   {
      glLoadIdentity();
      glTranslatef( 0, 0, -mDist);
      glMultMatrixf( &m[0][0] );
      const GLfloat colour0 [] = {0.00, 0.00, 0.00, 0.00}; 
      const GLfloat colour1 [] = {1.0f, 1.0f, 1.0f, 1.00}; 
      glMaterialfv(GL_FRONT, GL_AMBIENT,   colour0); 
      glMaterialfv(GL_FRONT, GL_DIFFUSE,   colour0); 
      glMaterialfv(GL_FRONT, GL_SPECULAR,  colour0); 
      glMaterialfv(GL_FRONT, GL_EMISSION,  colour1); 
      glMaterialfv(GL_FRONT, GL_SHININESS, colour0);
      glBegin(GL_LINES);
         glVertex3f(-1.0f, 0.0f, 0.0f);
         glVertex3f( 1.0f, 0.0f, 0.0f);

         glVertex3f( 0.0f,-1.0f, 0.0f);
         glVertex3f( 0.0f, 1.0f, 0.0f);

         glVertex3f( 0.0f, 0.0f,-1.0f);
         glVertex3f( 0.0f, 0.0f, 1.0f);
      glEnd();
   }
   // Print position of center of image, plus intensity of Fourier maps (if any)
   {
      wxString statusText;
      REAL x=mX0;
      REAL y=mY0;
      REAL z=mZ0;
      mpWXCrystal->GetCrystal().OrthonormalToFractionalCoords(x,y,z);
      x=0.5-x;
      y=0.5-y;
      z=0.5-z;
      statusText.sprintf("Center@(%5.3f,%5.3f,%5.3f)",x,y,z);
      for(unsigned int i=0;i<mvpUnitCellMapImport.size();++i)
      {
         wxString tmp;
         tmp=statusText;
         statusText.sprintf("%s, map(%s)=%5.2fe",tmp.c_str(),
                            mvpUnitCellMapImport[i]->GetName().c_str(),
                            mvpUnitCellMapImport[i]->GetValue(x,y,z));
      }
      mpParentFrame->SetStatusText(statusText);
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
   if(event.Leaving()) return;// wxMSW2.4 bug ?
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
            
            REAL vx1=mTrackBallLastX,vy1=mTrackBallLastY,vz1,
                 vx2=event.GetX(),   vy2=event.GetY(),   vz2;
            
            this->UnProject(vx1,vy1,vz1);
            this->UnProject(vx2,vy2,vz2);

            mX0 += vx2-vx1;
            mY0 += vy2-vy1;
            mZ0 += vz2-vz1;

	         VFN_DEBUG_MESSAGE("Origin (ortho) = "<<mX0<<", "<<mY0<<", "<<mZ0,2)
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
                        (mDist>101)?(mDist-100):1.,mDist+100);   
         Refresh(FALSE);
         VFN_DEBUG_MESSAGE(mViewAngle <<" "<<mDist,2)
      }
   }
   if(event.Leaving()) cout<<"Mouse is leaving window!!"<<endl;
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
   mpWXCrystal->UpdateGL(false,
			 mcellbbox.xMin,mcellbbox.xMax,
			 mcellbbox.yMin,mcellbbox.yMax,
			 mcellbbox.zMin,mcellbbox.zMax);
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

void WXGLCrystalCanvas::SetCurrent()
{
   VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::SetCurrent()",4)
   this->wxGLCanvas::SetCurrent();
   #ifndef HAVE_GLUT
   this->BuildGLFont();
   sFontDisplayListBase=mGLFontDisplayListBase;
   #endif
}

void WXGLCrystalCanvas::InitGL()
{
   VFN_DEBUG_ENTRY("WXGLCrystalCanvas::InitGL()",8)
   this->SetCurrent();
    
   glEnable(GL_DEPTH_TEST);
   glEnable(GL_LIGHTING);
   
   const GLfloat colour_Ambient [] = {0.4, 0.4, 0.4, 1.00}; 
   const GLfloat colour_Diffuse [] = {0.6, 0.6, 0.6, 1.00}; 
   const GLfloat colour_Specular[] = {0.2, 0.2, 0.2, 1.00}; 

   glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, 1.0); 

   const GLfloat LightPosition[]= { -10.0f, 10.0f, 10.0f, 0.0f };   
   glLightfv(GL_LIGHT1, GL_AMBIENT,  colour_Ambient); 
   glLightfv(GL_LIGHT1, GL_DIFFUSE,  colour_Diffuse); 
   glLightfv(GL_LIGHT1, GL_SPECULAR, colour_Specular);  
   glLightfv(GL_LIGHT1, GL_SHININESS,colour_Specular);  
   glLightfv(GL_LIGHT1, GL_POSITION,LightPosition);
   glEnable(GL_LIGHT1);  

   glEnable(GL_NORMALIZE);
   glHint(GL_PERSPECTIVE_CORRECTION_HINT,GL_NICEST);//GL_FASTEST
   glHint(GL_POLYGON_SMOOTH_HINT,GL_NICEST);//GL_FASTEST
   
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
  VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnChangeLimits()",10)
  UserSelectBoundingBox *BoxDlg = new UserSelectBoundingBox(this,
	    "Set bounding box for display of\natoms (fractional coordinates)",
			 mcellbbox);
  if (BoxDlg->ShowModal() == wxID_OK ) {
    mcellbbox =  BoxDlg->GetBBox();
    mpWXCrystal->UpdateGL(false,
			 mcellbbox.xMin,mcellbbox.xMax,
			 mcellbbox.yMin,mcellbbox.yMax,
			 mcellbbox.zMin,mcellbbox.zMax);
    this->CrystUpdate();
    VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnChangeLimits (X: " << 
		      mcellbbox.xMin << ", " << mcellbbox.xMax << 
		      " Y: " << 
		      mcellbbox.yMin << ", " << mcellbbox.yMax << 
		      " Z: " << 
		      mcellbbox.zMin << ", " << mcellbbox.zMax << 
		      ")", 10)
  } 
  BoxDlg->Destroy();
  VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnChangeLimits():UserSelectBoundingBox done",10)
}

void WXGLCrystalCanvas::OnShowCrystal()
{
   if(mShowCrystal) mpPopUpMenu->SetLabel(ID_GLCRYSTAL_MENU_SHOWCRYSTAL, "Show Crystal");
   else mpPopUpMenu->SetLabel(ID_GLCRYSTAL_MENU_SHOWCRYSTAL, "Hide Crystal");
   mShowCrystal = !mShowCrystal;
   this->CrystUpdate();
}

void WXGLCrystalCanvas::OnShowAtomLabel()
{
   if(mShowAtomName) mpPopUpMenu->SetLabel(ID_GLCRYSTAL_MENU_SHOWATOMLABEL, "Show Atom Labels");
   else mpPopUpMenu->SetLabel(ID_GLCRYSTAL_MENU_SHOWATOMLABEL, "Hide Atom Labels");
   mShowAtomName= !mShowAtomName;
   this->CrystUpdate();
}

void WXGLCrystalCanvas::OnShowCursor()
{
   if(mShowCursor) mpPopUpMenu->SetLabel(ID_GLCRYSTAL_MENU_SHOWCURSOR, "Show Cursor");
   else mpPopUpMenu->SetLabel(ID_GLCRYSTAL_MENU_SHOWCURSOR, "Hide Cursor");
   mShowCursor= !mShowCursor;
   this->CrystUpdate();
}

void WXGLCrystalCanvas::OnSetCursor()
{
  VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnSetCursor",1)
  REAL x=mX0;
  REAL y=mY0;
  REAL z=mZ0;
  mpWXCrystal->GetCrystal().OrthonormalToFractionalCoords(x,y,z);

  mViewCntr.x = 0.5 - x;
  mViewCntr.y = 0.5 - y;
  mViewCntr.z = 0.5 - z;
  UserXYZBox *BoxDlg = new UserXYZBox(this,
       "Set fractional coordinates for view\ncenter and cursor position",
			 mViewCntr);
  if (BoxDlg->ShowModal() == wxID_OK ) {
     mViewCntr =  BoxDlg->GetXYZ();
     VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnSetCursor (frac) = " <<
		 mViewCntr.x << "," << mViewCntr.y << "," << mViewCntr.z,1)
     mX0 = 0.5 - mViewCntr.x;
     mY0 = 0.5 - mViewCntr.y;
     mZ0 = 0.5 - mViewCntr.z;
     mpWXCrystal->GetCrystal().FractionalToOrthonormalCoords(mX0, mY0, mZ0);
     VFN_DEBUG_MESSAGE("...ortho" << mX0 << "," << mY0 << "," << mZ0,1)
     Refresh(FALSE);
  }
}

void WXGLCrystalCanvas::OnLoadFourier()
{
   wxFileDialog fd((wxWindow*)this, "Choose a file containing a Fourier Map",
           "", "", "Fourier Map files (*.grd)|*.grd", wxOPEN | wxHIDE_READONLY | wxFILE_MUST_EXIST);
   //if okay then read Fourier map, run MC on it and display the triangles
   if(fd.ShowModal() == wxID_OK)
   {
     VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnLoadFourier opening "<<fd.GetPath().c_str(),10)
      this->LoadFourier((string)(fd.GetPath().c_str()));
   }
}

void WXGLCrystalCanvas::LoadFourier(const string&filename)
{//read fourier map from 0.0 to 1.0 on all axis
   wxBusyInfo wait("Processing Fourier Map...");
   {
      //auto_ptr<UnitCellMapImport> ptr(new UnitCellMapImport(mpWXCrystal->GetCrystal()));
      mvpUnitCellMapImport.push_back(new UnitCellMapImport(mpWXCrystal->GetCrystal()));
      // load the map; exit on error
      if (mvpUnitCellMapImport.back()->ImportGRD(filename) == 0) {
	// how can I insert filename into this message?
	wxMessageBox("Error reading Fourier file", "File error", wxOK, this);
	return;
      }
   }
   {
      //auto_ptr<UnitCellMapGLList> ptr(new UnitCellMapGLList);
      mvpUnitCellMapGLList.push_back(make_pair(make_pair(mvpUnitCellMapImport.back(),1.0),
                                               new UnitCellMapGLList) );
      switch(mvpUnitCellMapGLList.size())
      {
         case 1: mvpUnitCellMapGLList.back().second->SetColour(1.,0.,0.,1.);break;
         case 2: mvpUnitCellMapGLList.back().second->SetColour(0.,0.,1.,1.);break;
         default:mvpUnitCellMapGLList.back().second->SetColour(0.,1.,0.,1.);break;
      }
      this->SetCurrent();
      mvpUnitCellMapGLList.back().second->GenList(*(mvpUnitCellMapImport.back()),
		     this, mvpUnitCellMapGLList.back().first.second);
   }

   mpPopUpMenu->Enable(ID_GLCRYSTAL_MENU_CHANGECONTOUR, TRUE);
   mpPopUpMenu->Enable(ID_GLCRYSTAL_MENU_ADDCONTOUR, TRUE);
   mpPopUpMenu->Enable(ID_GLCRYSTAL_MENU_SHOWFOURIER, TRUE);
   mpPopUpMenu->Enable(ID_GLCRYSTAL_MENU_FOURIERCHANGECOLOR, TRUE);
   mpPopUpMenu->Enable(ID_GLCRYSTAL_MENU_UNLOADFOURIER, TRUE);
   mpPopUpMenu->Enable(ID_GLCRYSTAL_MENU_SHOWWIRE, TRUE);
   mpPopUpMenu->Enable(ID_GLCRYSTAL_MENU_FOURIERCHANGEBBOX, TRUE);
   this->CrystUpdate();
}

void WXGLCrystalCanvas::OnChangeContour()
{
   int mapgl=this->UserSelectUnitCellMapGLList();
   
   double contourValue=(double) (mvpUnitCellMapGLList[mapgl].first.second);
   wxString strValue;
   strValue.Printf("%lf",contourValue);
   wxTextEntryDialog contourValueDialog(this,"New contour value",
                           "New contour value",strValue,wxOK | wxCANCEL);
   if(wxID_OK!=contourValueDialog.ShowModal())
   {
      return;
   }
   wxBusyInfo wait("Processing Fourier Map...");
   contourValueDialog.GetValue().ToDouble(&contourValue);
   mvpUnitCellMapGLList[mapgl].first.second = (float) contourValue;
   mvpUnitCellMapGLList[mapgl].second->GenList(*(mvpUnitCellMapImport.back()),
					       this, mvpUnitCellMapGLList[mapgl].first.second);
   this->CrystUpdate();
}

void WXGLCrystalCanvas::OnAddContour()
{
   int map=this->UserSelectUnitCellMapImport();
   
   // Choose contour value
      double contourValue=1.;
      wxString strValue;
      strValue.Printf("%lf",contourValue);
      wxTextEntryDialog contourValueDialog(this,"Add contour value",
                              "Add contour value",strValue,wxOK | wxCANCEL);
      if(wxID_OK!=contourValueDialog.ShowModal())  return;
      contourValueDialog.GetValue().ToDouble(&contourValue);
   // Choose colour
      wxColor ncolor(255,0,0);
      ncolor = wxGetColourFromUser((wxWindow*)this, ncolor);   
      if(!(ncolor.Ok())) return;
   // Add display map
      wxBusyInfo wait("Processing Fourier Map...");
      mvpUnitCellMapGLList.push_back(make_pair(make_pair(mvpUnitCellMapImport.back(),
                                                         (float)contourValue),
                                               new UnitCellMapGLList()) );
      mvpUnitCellMapGLList.back().second->SetName(mvpUnitCellMapImport[map]->GetName());
      mvpUnitCellMapGLList.back().second
         ->SetColour(ncolor.Red()/255.0,ncolor.Green()/255.0,ncolor.Blue()/255.0,0.5);
      this->SetCurrent();
      mvpUnitCellMapGLList.back().second->GenList(*mvpUnitCellMapImport[map],
						  this, mvpUnitCellMapGLList.back().first.second);
      this->CrystUpdate();
}

void WXGLCrystalCanvas::OnShowFourier()
{
   if(mShowFourier == TRUE) mpPopUpMenu->SetLabel(ID_GLCRYSTAL_MENU_SHOWFOURIER, "Show Fourier Map");
   else mpPopUpMenu->SetLabel(ID_GLCRYSTAL_MENU_SHOWFOURIER, "Hide Fourier Map");
   mShowFourier = !mShowFourier;
   this->CrystUpdate();
}

void WXGLCrystalCanvas::OnFourierChangeColor()
{
   int mapgl=this->UserSelectUnitCellMapGLList();
   wxColor ncolor((char)(255*mvpUnitCellMapGLList[mapgl].second->GetColour()[0]),
                  (char)(255*mvpUnitCellMapGLList[mapgl].second->GetColour()[1]),
                  (char)(255*mvpUnitCellMapGLList[mapgl].second->GetColour()[2]));
   ncolor = wxGetColourFromUser((wxWindow*)this, ncolor);   
   if(ncolor.Ok()) 
   {
      mvpUnitCellMapGLList[mapgl].second
         ->SetColour(ncolor.Red()/255.0,ncolor.Green()/255.0,ncolor.Blue()/255.0,0.5);
      this->CrystUpdate();
   }
}

void WXGLCrystalCanvas::OnFourierChangeBbox()
{
  VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnFourierChangeBbox()",10)
  // change Xmax if in default mode
    BBox bbox = mmapbbox;
  if (bbox.xMin == bbox.xMax) bbox.xMax += 1.0;
  UserSelectBoundingBox *BoxDlg = new UserSelectBoundingBox(this,
      "Set bounding box for display of\nFourier map (fractional coordinates)",
       bbox);
  if (BoxDlg->ShowModal() == wxID_OK ) {
    mmapbbox =  BoxDlg->GetBBox();
    vector<pair<pair<const UnitCellMapImport*,float>,UnitCellMapGLList* > >::iterator pos;
    for(pos=mvpUnitCellMapGLList.begin();pos != mvpUnitCellMapGLList.end();pos++)
      pos->second->GenList(*(mvpUnitCellMapImport.back()),
			   this, pos->first.second);

    this->CrystUpdate();
    VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnFourierChangeBbox (X: " << 
		      mmapbbox.xMin << ", " << mmapbbox.xMax << 
		      " Y: " << 
		      mmapbbox.yMin << ", " << mmapbbox.yMax << 
		      " Z: " << 
		      mmapbbox.zMin << ", " << mmapbbox.zMax << 
		      ")", 10)
  }
  BoxDlg->Destroy();
  VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnFourierChangeBbox done",10)
}

void WXGLCrystalCanvas::OnUnloadFourier()
{
   wxMessageDialog * msure = new wxMessageDialog((wxWindow*)this,
     "Are you sure you want to unload all Fourier Map Data?", "Unload Fourier Map", wxYES_NO | wxNO_DEFAULT |
     wxICON_QUESTION );
   if(msure->ShowModal() == wxID_YES)
   {
      {
         vector<UnitCellMapImport*>::iterator pos;
         for(pos=mvpUnitCellMapImport.begin();pos != mvpUnitCellMapImport.end();pos++)
            delete *pos;
         mvpUnitCellMapImport.clear();
      }
      {
         vector<pair<pair<const UnitCellMapImport*,float>,UnitCellMapGLList* > >::iterator pos;
         for(pos=mvpUnitCellMapGLList.begin();pos != mvpUnitCellMapGLList.end();pos++)
            delete pos->second;
         mvpUnitCellMapGLList.clear();
      }
      mpPopUpMenu->SetLabel(ID_GLCRYSTAL_MENU_SHOWCRYSTAL, "Hide Crystal");
      mpPopUpMenu->Enable(ID_GLCRYSTAL_MENU_CHANGECONTOUR, FALSE);      //disable all of these
      mpPopUpMenu->Enable(ID_GLCRYSTAL_MENU_ADDCONTOUR, FALSE);      //disable all of these
      mpPopUpMenu->Enable(ID_GLCRYSTAL_MENU_SHOWFOURIER, FALSE);
      mpPopUpMenu->SetLabel(ID_GLCRYSTAL_MENU_SHOWFOURIER, "Hide Fourier Map");
      mpPopUpMenu->Enable(ID_GLCRYSTAL_MENU_UNLOADFOURIER, FALSE);
      mpPopUpMenu->Enable(ID_GLCRYSTAL_MENU_FOURIERCHANGECOLOR, FALSE);
      mpPopUpMenu->SetLabel(ID_GLCRYSTAL_MENU_SHOWWIRE, "Show Filled");
      mpPopUpMenu->Enable(ID_GLCRYSTAL_MENU_SHOWWIRE, FALSE);
      mpPopUpMenu->Enable(ID_GLCRYSTAL_MENU_FOURIERCHANGEBBOX, FALSE);

      this->CrystUpdate();
   }
   delete msure;
}

void WXGLCrystalCanvas::OnShowWire()
{
   vector<pair<pair<const UnitCellMapImport*,float>,UnitCellMapGLList* > >::iterator pos;
   for(pos=mvpUnitCellMapGLList.begin();pos != mvpUnitCellMapGLList.end();pos++)
      pos->second->ToggleShowWire();

   this->CrystUpdate();
}

int WXGLCrystalCanvas::UserSelectUnitCellMapGLList()const
{
   int mapgl=0;
   if(mvpUnitCellMapGLList.size()>1)
   {
      wxString *pChoices=new wxString[mvpUnitCellMapGLList.size()];
      for(unsigned int i=0;i<mvpUnitCellMapGLList.size();i++) 
         (pChoices+i)->Printf("%s:contour=%5.3f,rgb=(%5.3f,%5.3f,%5.3f)",
                           mvpUnitCellMapGLList[i].second->GetName().c_str(),
                           mvpUnitCellMapGLList[i].first.second,
                           mvpUnitCellMapGLList[i].second->GetColour()[0],
                           mvpUnitCellMapGLList[i].second->GetColour()[1],
                           mvpUnitCellMapGLList[i].second->GetColour()[2]);
      wxSingleChoiceDialog dialog
         ((wxWindow*)this,"Choose displayed map","Choose displayed map",
          mvpUnitCellMapGLList.size(),pChoices,0,wxOK);
      dialog.ShowModal();
      mapgl=dialog.GetSelection();
      delete[] pChoices;
   }
   return mapgl;
}


int WXGLCrystalCanvas::UserSelectUnitCellMapImport()const
{
   int map=0;
   if(1<mvpUnitCellMapImport.size())
   {
      wxString *pChoices=new wxString[mvpUnitCellMapImport.size()];
      for(unsigned int i=0;i<mvpUnitCellMapImport.size();i++) 
         *(pChoices+i) = mvpUnitCellMapImport[i]->GetName().c_str();
      wxSingleChoiceDialog dialog
         ((wxWindow*)this,"Choose map","Choose map",
          mvpUnitCellMapImport.size(),pChoices,0,wxOK);
      dialog.ShowModal();
      map=dialog.GetSelection();
      delete[] pChoices;
   }
   return map;
}

BBox WXGLCrystalCanvas::GetCellBBox() {
  return mcellbbox;
}
BBox WXGLCrystalCanvas::GetMapBBox() {
  return mmapbbox;
}
void WXGLCrystalCanvas::UnProject(REAL &x, REAL &y, REAL &z)
{
   GLdouble vx,vy,vz,junk;
   GLdouble z0;
   GLdouble modelMatrix[16];
   GLdouble projMatrix[16];
   GLint viewport[4];

   this->SetCurrent();
   glMatrixMode( GL_MODELVIEW );
   glLoadIdentity();
   glTranslatef( 0, 0, -mDist );

   glGetDoublev(GL_MODELVIEW_MATRIX,modelMatrix);
   glGetDoublev(GL_PROJECTION_MATRIX,projMatrix);
   glGetIntegerv(GL_VIEWPORT,viewport);

   // First, get the z depth of where we want to translate
   gluProject(0, 0, 0 ,modelMatrix,projMatrix,viewport,&junk,&junk,&z0);
   // Get the orthonormal coordinates 
   gluUnProject(x,y,z0,modelMatrix,projMatrix,viewport,&vx,&vy,&vz);
   vy = -vy;
   // Use Quaternion to get the correct position
   GLfloat m[4][4];
   build_rotmatrix( m,mQuat);

   x= m[0][0]* vx + m[0][1]*vy + m[0][2]*vz -mX0;
   y= m[1][0]* vx + m[1][1]*vy + m[1][2]*vz -mY0;
   z= m[2][0]* vx + m[2][1]*vy + m[2][2]*vz -mZ0;
	VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::UnProject():X Y Z = "<<x<<" , "<<y<<" , "<<z,5)
}
#ifndef HAVE_GLUT
void WXGLCrystalCanvas::BuildGLFont()
{
   if(mIsGLFontBuilt) return;
   VFN_DEBUG_ENTRY("WXGLCrystalCanvas::BuildGLFont()-gldisplay",6)
   #ifdef __LINUX__
      Display *dpy;
      XFontStruct *fontInfo=NULL;

      mGLFontDisplayListBase = glGenLists(96);

      dpy = XOpenDisplay(NULL); 

      fontInfo = XLoadQueryFont(dpy, "-adobe-helvetica-bold-*-r-*-10-*-*-*-*-*-*-*");
      if (fontInfo == NULL)
         fontInfo = XLoadQueryFont(dpy, "-adobe-helvetica-bold-*-*-*-10-*-*-*-*-*-*-*");
      if (fontInfo == NULL)
         fontInfo = XLoadQueryFont(dpy, "-adobe-times-bold-*-r-*-10-*-*-*-*-*-*-*");
      if (fontInfo == NULL)
         fontInfo = XLoadQueryFont(dpy, "-adobe-helvetica-medium-*-*-*-12-*-*-*-*-*-*-*");
      if (fontInfo == NULL)
         fontInfo = XLoadQueryFont(dpy, "-adobe-times-medium-*-*-*-12-*-*-*-*-*-*-*");
      if (fontInfo == NULL)
         fontInfo = XLoadQueryFont(dpy, "-adobe-helvetica-*-*-*-*-12-*-*-*-*-*-*-*");
      if (fontInfo == NULL)
         fontInfo = XLoadQueryFont(dpy, "-adobe-times-*-*-*-*-12-*-*-*-*-*-*-*");
      if (fontInfo == NULL)
         fontInfo = XLoadQueryFont(dpy, "fixed");
	   if (fontInfo == NULL) cout <<"no X font available..."<<endl;

      glXUseXFont(fontInfo->fid, 32, 96, mGLFontDisplayListBase);
      XFreeFont(dpy, fontInfo);
      XCloseDisplay(dpy);
   #endif
   #ifdef __WIN32__
      HFONT   font;
      HFONT   oldfont;
      wxPaintDC dc(this);
      HDC hDC = (HDC)dc.GetHDC();
      mGLFontDisplayListBase = glGenLists(96);
      font = CreateFont(-12,                       // Height of font
                        0,                         // Width of font
                        0,                         // Angle of escapement
                        0,                         // Orientation angle
                        FW_BOLD,                   // Font weight
                        FALSE,                     // Italic
                        FALSE,                     // Underline
                        FALSE,                     // Strikeout
                        ANSI_CHARSET,              // Character set identifier
                        OUT_TT_PRECIS,             // Output precision
                        CLIP_DEFAULT_PRECIS,       // Clipping precision
                        ANTIALIASED_QUALITY,       // Output quality
                        FF_DONTCARE|DEFAULT_PITCH, // Family and pitch
                        "Helvetica");              // Font name

      oldfont = (HFONT)SelectObject(hDC, font);
      wglUseFontBitmaps(hDC, 32, 96, mGLFontDisplayListBase);
      SelectObject(hDC, oldfont);
      DeleteObject(font);
   #endif
   mIsGLFontBuilt=true;
   sFontDisplayListBase=mGLFontDisplayListBase;
   VFN_DEBUG_EXIT("WXGLCrystalCanvas::BuildGLFont()",6)
}

void WXGLCrystalCanvas::DeleteGLFont() const
{
   if(!mIsGLFontBuilt) return;
   glDeleteLists(mGLFontDisplayListBase, 96);
   mIsGLFontBuilt=false;
   mGLFontDisplayListBase=0;
}
#endif


////////////////////////////////////////////////////////////////////////
//
//    UserSelectBoundingBox
//
////////////////////////////////////////////////////////////////////////
BEGIN_EVENT_TABLE(UserSelectBoundingBox, wxDialog)
   EVT_BUTTON(wxID_OK, UserSelectBoundingBox::OnOk)
END_EVENT_TABLE()

  UserSelectBoundingBox::UserSelectBoundingBox (wxWindow *parent, char * title,
					      const BBox bbox)
  : wxDialog((wxWindow *)parent, -1, "Set bounding box", wxDefaultPosition,
  	     wxSize(250, 250), wxDEFAULT_DIALOG_STYLE) 
{
  wxBoxSizer *dialogSizer = new wxBoxSizer(wxVERTICAL);
  wxFlexGridSizer *inputSizer = new wxFlexGridSizer(4, 3, 10, 10);
  // headers
  inputSizer->Add(new wxStaticText(this, -1, ""), 0, wxALIGN_CENTRE_VERTICAL);
  inputSizer->Add(new wxStaticText(this, -1, "minimum"), 0, wxALIGN_CENTER);
  inputSizer->Add(new wxStaticText(this, -1, "maximum"), 0, wxALIGN_CENTER);
  // 1st row
  inputSizer->Add(new wxStaticText(this, -1, "a"), 0, wxALIGN_CENTRE_VERTICAL);
  inputSizer->Add(mpXminCtrl = new wxTextCtrl(this, -1, 
					      wxString::Format("%f",bbox.xMin)), 
					      0, wxALIGN_CENTRE_VERTICAL);
  inputSizer->Add(mpXmaxCtrl = new wxTextCtrl(this, -1, 
					      wxString::Format("%f",bbox.xMax)), 
					      0, wxALIGN_CENTRE_VERTICAL);
  // 2nd row
  inputSizer->Add(new wxStaticText(this, -1, "b"), 0, wxALIGN_CENTRE_VERTICAL);
  inputSizer->Add(mpYminCtrl = new wxTextCtrl(this, -1, 
					      wxString::Format("%f",bbox.yMin)), 
					      0, wxALIGN_CENTRE_VERTICAL);
  inputSizer->Add(mpYmaxCtrl = new wxTextCtrl(this, -1, 
					      wxString::Format("%f",bbox.yMax)), 
					      0, wxALIGN_CENTRE_VERTICAL);
  // 3rd row
  inputSizer->Add(new wxStaticText(this, -1, "c"), 0, wxALIGN_CENTRE_VERTICAL);
  inputSizer->Add(mpZminCtrl = new wxTextCtrl(this, -1, 
					      wxString::Format("%f",bbox.zMin)), 
					      0, wxALIGN_CENTRE_VERTICAL);
  inputSizer->Add(mpZmaxCtrl = new wxTextCtrl(this, -1, 
					      wxString::Format("%f",bbox.zMax)), 
					      0, wxALIGN_CENTRE_VERTICAL);
  // button section
  wxFlexGridSizer *buttonSizer = new wxFlexGridSizer(1, 2, 10, 10);
  buttonSizer->Add(new wxButton(this, wxID_OK, "OK"), 
		   0, wxALIGN_CENTRE_VERTICAL);
  buttonSizer->Add(new wxButton(this, wxID_CANCEL, "Cancel"), 
		   0, wxALIGN_CENTRE_VERTICAL);

  dialogSizer->Add(10, 10);
  dialogSizer->Add(new wxStaticText(this, -1, title), 0, 
		   wxALIGN_CENTER);
  dialogSizer->Add(10, 10);
  dialogSizer->Add(inputSizer, 0, wxALIGN_CENTER);
  dialogSizer->Add(20, 20);
  dialogSizer->Add(buttonSizer, 0, wxALIGN_CENTER);

  SetSizer(dialogSizer);
  SetAutoLayout(TRUE);
  Layout();
}

UserSelectBoundingBox::~UserSelectBoundingBox () {
  // if the sizers must be deleted, put them in the class and delete them here
  //delete dialogSizer;
  //delete inputSizer;
  //delete buttonSizer;
};

void UserSelectBoundingBox::OnOk () {
  char * strptr;
  const char * val;

  val = mpXminCtrl->GetValue().c_str();
  mbbox.xMin = strtod(val, &strptr);
  if (val == strptr) {wxMessageBox("Invalid value for Xmin!", "Bounding volume error", wxOK, this); return;}
  val = mpXmaxCtrl->GetValue().c_str();
  mbbox.xMax = strtod(val, &strptr);
  if (val == strptr) {wxMessageBox("Invalid value for Xmax!", "Bounding volume error", wxOK, this); return;}
  if (mbbox.xMin == mbbox.xMax) {wxMessageBox("Sorry, Xmin must be less than Xmax!", "Zero bounding volume", wxOK, this); return;}
  if (mbbox.xMin > mbbox.xMax) {
    float tmp = mbbox.xMax;
    mbbox.xMax = mbbox.xMin;
    mbbox.xMin = tmp;
  }
  VFN_DEBUG_MESSAGE("Xmin " << mbbox.xMin << " Xmax " << mbbox.xMax,1)

  val = mpYminCtrl->GetValue().c_str();
  mbbox.yMin = strtod(val, &strptr);
  if (val == strptr) {wxMessageBox("Invalid value for Ymin!", "Bounding volume error", wxOK, this); return;}
  val = mpYmaxCtrl->GetValue().c_str();
  mbbox.yMax = strtod(val, &strptr);
  if (val == strptr) {wxMessageBox("Invalid value for Ymax!", "Bounding volume error", wxOK, this); return;}
  if (mbbox.yMin == mbbox.yMax) {wxMessageBox("Sorry, Ymin must be less than Ymax!", "Zero bounding volume", wxOK, this); return;}
  if (mbbox.yMin > mbbox.yMax) {
    float tmp = mbbox.yMax;
    mbbox.yMax = mbbox.yMin;
    mbbox.yMin = tmp;
  }
  VFN_DEBUG_MESSAGE("Ymin " << mbbox.yMin << " Ymax " << mbbox.yMax,1)

  val = mpZminCtrl->GetValue().c_str();
  mbbox.zMin = strtod(val, &strptr);
  if (val == strptr) {wxMessageBox("Invalid value for Zmin!", "Bounding volume error", wxOK, this); return;}
  val = mpZmaxCtrl->GetValue().c_str();
  mbbox.zMax = strtod(val, &strptr);
  if (val == strptr) {wxMessageBox("Invalid value for Zmax!", "Bounding volume error", wxOK, this); return;}
  if (mbbox.zMin == mbbox.zMax) {wxMessageBox("Sorry, Zmin must be less than Zmax!", "Zero bounding volume", wxOK, this); return;}
  if (mbbox.zMin > mbbox.zMax) {
    float tmp = mbbox.zMax;
    mbbox.zMax = mbbox.zMin;
    mbbox.zMin = tmp;
  }
  VFN_DEBUG_MESSAGE("Zmin " << mbbox.zMin << " Zmax " << mbbox.zMax,1)

    // close the dialog
    EndModal(wxID_OK);
}

BBox UserSelectBoundingBox::GetBBox () {
  return mbbox;
}


////////////////////////////////////////////////////////////////////////
//
//    UserXYZBox
//
////////////////////////////////////////////////////////////////////////
BEGIN_EVENT_TABLE(UserXYZBox, wxDialog)
   EVT_BUTTON(wxID_OK, UserXYZBox::OnOk)
END_EVENT_TABLE()

  UserXYZBox::UserXYZBox (wxWindow *parent, char * title,
					      const Triple xyz)
  : wxDialog((wxWindow *)parent, -1, "Set position", wxDefaultPosition,
  	     wxSize(250, 250), wxDEFAULT_DIALOG_STYLE) 
{
  wxBoxSizer *dialogSizer = new wxBoxSizer(wxVERTICAL);
  wxFlexGridSizer *inputSizer = new wxFlexGridSizer(3, 2, 10, 10);
  // 1st row
  inputSizer->Add(new wxStaticText(this, -1, "x"), 0, wxALIGN_CENTRE_VERTICAL);
  inputSizer->Add(mpXCtrl = new wxTextCtrl(this, -1, 
					   wxString::Format("%.3f",xyz.x)), 
					   0, wxALIGN_CENTRE_VERTICAL);
  // 2nd row
  inputSizer->Add(new wxStaticText(this, -1, "y"), 0, wxALIGN_CENTRE_VERTICAL);
  inputSizer->Add(mpYCtrl = new wxTextCtrl(this, -1, 
					   wxString::Format("%.3f",xyz.y)), 
					   0, wxALIGN_CENTRE_VERTICAL);
  // 3rd row
  inputSizer->Add(new wxStaticText(this, -1, "z"), 0, wxALIGN_CENTRE_VERTICAL);
  inputSizer->Add(mpZCtrl = new wxTextCtrl(this, -1, 
					   wxString::Format("%.3f",xyz.z)), 
					   0, wxALIGN_CENTRE_VERTICAL);
  // button section
  wxFlexGridSizer *buttonSizer = new wxFlexGridSizer(1, 2, 10, 10);
  buttonSizer->Add(new wxButton(this, wxID_OK, "OK"), 
		   0, wxALIGN_CENTRE_VERTICAL);
  buttonSizer->Add(new wxButton(this, wxID_CANCEL, "Cancel"), 
		   0, wxALIGN_CENTRE_VERTICAL);

  dialogSizer->Add(10, 10);
  dialogSizer->Add(new wxStaticText(this, -1, title), 0, 
		   wxALIGN_CENTER);
  dialogSizer->Add(10, 10);
  dialogSizer->Add(inputSizer, 0, wxALIGN_CENTER);
  dialogSizer->Add(20, 20);
  dialogSizer->Add(buttonSizer, 0, wxALIGN_CENTER);

  SetSizer(dialogSizer);
  SetAutoLayout(TRUE);
  Layout();
}

UserXYZBox::~UserXYZBox () {
};

void UserXYZBox::OnOk () {
  char * strptr;
  const char * val;

  val = mpXCtrl->GetValue().c_str();
  mXYZ.x = strtod(val, &strptr);
  if (val == strptr) {wxMessageBox("Invalid value for X!", "Position error", wxOK, this); return;}

  val = mpYCtrl->GetValue().c_str();
  mXYZ.y = strtod(val, &strptr);
  if (val == strptr) {wxMessageBox("Invalid value for Y!", "Position error", wxOK, this); return;}

  val = mpZCtrl->GetValue().c_str();
  mXYZ.z = strtod(val, &strptr);
  if (val == strptr) {wxMessageBox("Invalid value for Z!", "Position error", wxOK, this); return;}

    // close the dialog
    EndModal(wxID_OK);
}

Triple UserXYZBox::GetXYZ () {
  return mXYZ;
}


#endif // #ifdef OBJCRYST_GL

}// namespace 
