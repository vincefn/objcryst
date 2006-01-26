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
   #ifdef __WXGTK__
      #include "GL/glu.h"
   #endif
   
   #ifdef __WXMAC__ // For the wxMac version of wxWindows, i.e. with the "Aqua" look
      #include <OpenGL/glu.h>
      #include "AGL/agl.h"
   #endif

   #ifdef __LINUX__
      #include "GL/glx.h"
   #endif
   
   #ifdef __WIN32__
     #include "gl/glaux.h"
   #endif
   
   #ifdef HAVE_GLUT
      #ifdef __WXMAC__ // For the wxMac version of wxWindows, i.e. with the "Aqua" look
         #include <GLUT/glut.h>
      #else
         #include "GL/glut.h"
       #endif
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
      void OnOk (wxCommandEvent & WXUNUSED(event));
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
      void OnOk (wxCommandEvent & WXUNUSED(event));
      wxTextCtrl * mpXCtrl;
      wxTextCtrl * mpYCtrl;
      wxTextCtrl * mpZCtrl;
      Triple mXYZ;
      DECLARE_EVENT_TABLE()
  };

////////////////////////////////////////////////////////////////////////
//
//    WXMolScrolledWindow
//
////////////////////////////////////////////////////////////////////////
WXCrystalScrolledGridWindow::WXCrystalScrolledGridWindow(wxWindow* parent, WXCrystal* p, long id):
wxGrid(parent,id),mpWXCrystal(p)
{}

WXCrystalScrolledGridWindow::~WXCrystalScrolledGridWindow()
{
   mpWXCrystal->NotifyDeleteListWin(this);
}

////////////////////////////////////////////////////////////////////////
//
//    WXCrystal Grid objects
//
////////////////////////////////////////////////////////////////////////
WXCrystal::RowScattPow::RowScattPow():
mName("H"),
mBiso(1.0),mFormalCharge(0.0),mR(1.0),mG(1.0),mB(1.0),mMaximumLikelihoodError(0.0),mNbGhostAtoms(0.0),mNeedUpdateUI(true)
{}

WXCrystal::RowAntiBump::RowAntiBump():
mName(""),mNeedUpdateUI(true)
{}

WXCrystal::RowBondValence::RowBondValence():
mName(""),mNeedUpdateUI(true)
{}

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
static const long ID_CRYSTAL_WIN_SCATTPOW                       =WXCRYST_ID();
static const long ID_CRYSTAL_WIN_ANTIBUMP                       =WXCRYST_ID();
static const long ID_CRYSTAL_WIN_BONDVALENCE                    =WXCRYST_ID();
static const long ID_CRYSTAL_MENU_SHOW_SCATTPOW_WIN             =WXCRYST_ID();

BEGIN_EVENT_TABLE(WXCrystal,wxWindow)
   EVT_BUTTON(ID_WXOBJ_COLLAPSE,                      WXCrystObj::OnToggleCollapse)
   EVT_MENU(ID_REFOBJ_MENU_OBJ_SAVE,                  WXRefinableObj::OnMenuSave)
   EVT_MENU(ID_REFOBJ_MENU_OBJ_LOAD,                  WXRefinableObj::OnMenuLoad)
   EVT_MENU(ID_CRYSTAL_MENU_SAVECIF,                  WXCrystal::OnMenuSaveCIF)
   EVT_MENU(ID_CRYSTAL_MENU_SAVETEXT,                 WXCrystal::OnMenuSaveText)
   EVT_MENU(ID_REFOBJ_MENU_PAR_FIXALL,                WXRefinableObj::OnMenuFixAllPar)
   EVT_MENU(ID_REFOBJ_MENU_PAR_UNFIXALL,              WXRefinableObj::OnMenuUnFixAllPar)
   EVT_MENU(ID_REFOBJ_MENU_PAR_RANDOMIZE,             WXRefinableObj::OnMenuParRandomize)
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
   EVT_MENU(ID_CRYSTAL_MENU_SHOW_SCATTPOW_WIN,        WXCrystal::OnMenuShowScattPowWindow)
   EVT_UPDATE_UI(ID_CRYST_UPDATEUI,                   WXRefinableObj::OnUpdateUI)
   EVT_GRID_CMD_CELL_CHANGE(ID_CRYSTAL_WIN_SCATTPOW,  WXCrystal::OnEditGridScattPow)
   EVT_GRID_CMD_CELL_CHANGE(ID_CRYSTAL_WIN_ANTIBUMP,  WXCrystal::OnEditGridScattPowAntiBump)
   EVT_GRID_CMD_CELL_CHANGE(ID_CRYSTAL_WIN_BONDVALENCE,WXCrystal::OnEditGridScattPowBondValence)
END_EVENT_TABLE()

WXCrystal::WXCrystal(wxWindow* parent, Crystal *obj):
WXRefinableObj(parent,(RefinableObj*)obj),mpCrystal(obj),
mpScattPowWin(0),mpAntiBumpWin(0),mpBondValenceWin(0),
mIsSelfUpdating(false)
#ifdef OBJCRYST_GL
,mCrystalGLDisplayList(0),mCrystalGLNameDisplayList(0),
mpCrystalGL(0)
#endif
,mpConditionGLUpdate(0)
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
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_PAR,ID_REFOBJ_MENU_PAR_FIXALL,"Fix all");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_PAR,ID_REFOBJ_MENU_PAR_UNFIXALL,"Unfix all");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_PAR,ID_REFOBJ_MENU_PAR_RANDOMIZE,
                                "Randomize Configuration");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_PAR,ID_CRYSTAL_MENU_PAR_SETRELATIVEXYZLIMITS,
                                "Set Relative Limits On All XYZ Parameters");
      mpMenuBar->AddMenu("Scatterers",ID_CRYSTAL_MENU_SCATT);
         mpMenuBar->AddMenuItem(ID_CRYSTAL_MENU_SCATT,ID_CRYSTAL_MENU_SHOW_SCATTPOW_WIN,
                                "Show Scattering Powers Parameters Window");
         mpMenuBar->GetMenu(ID_CRYSTAL_MENU_SCATT).AppendSeparator();
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

      mpSizer->SetItemMinSize(mpMenuBar,
                              mpMenuBar->GetSize().GetWidth(),
                              mpMenuBar->GetSize().GetHeight());
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
   // Bond Valence cost
      wxBoxSizer* pBondValenceSizer=new wxBoxSizer(wxHORIZONTAL);
      WXFieldPar<REAL> *pWXFieldBondValence=
         new WXFieldPar<REAL>(this,"Bond Valence Cost",-1,&(mpCrystal->mBondValenceCost),100);
      WXFieldPar<REAL> *pBondValenceScale=
         new WXFieldPar<REAL>(this,"Scale",-1,&(mpCrystal->mBondValenceCostScale));
      pBondValenceSizer->Add(pWXFieldBondValence);
      pBondValenceSizer->Add(pBondValenceScale);
      mpSizer->Add(pBondValenceSizer,0,wxALIGN_LEFT);
      mList.Add(pWXFieldBondValence);
      mList.Add(pBondValenceScale);
   // Lattice
      wxBoxSizer* lattice=new wxBoxSizer(wxHORIZONTAL);

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
   
   this->BottomLayout(0);
   this->CrystUpdate(true);
   VFN_DEBUG_MESSAGE("WXCrystal::WXCrystal():End",6)
}

WXCrystal::~WXCrystal()
{
   VFN_DEBUG_ENTRY("WXCrystal::~WXCrystal()",10)
   if(0!=mpScattPowWin) mpScattPowWin->GetParent()->Destroy();
   VFN_DEBUG_EXIT("WXCrystal::~WXCrystal()",10)
}

void WXCrystal::CrystUpdate(const bool uui,const bool lock)
{
   VFN_DEBUG_ENTRY("WXCrystal::CrystUpdate()",5)
   mpCrystal->GetBumpMergeCost();
   mpCrystal->GetBondValenceCost();
   #ifdef OBJCRYST_GL
   if(mpCrystalGL!=0)
   {
      if(lock) mMutex.Lock();
      BBox box=mpCrystalGL->GetCellBBox();
      if(lock) mMutex.Unlock();
      this->UpdateGL(false,box.xMin,box.xMax,box.yMin,box.yMax,box.zMin,box.zMax);
   }
   #endif
   if(lock) mMutex.Lock();
   if(false==this->GetCrystal().IsBeingRefined())
   {
      // Update map of scatteringPowers, if needs be.
      //if(mvScattPowIndexClock<this->GetCrystal().GetScatteringPowerRegistry().GetRegistryClock())
      {
         mvScattPowIndex.clear();
         mvScattPowRowIndex.clear();
         int ix=0;
         for(int i=0;i<this->GetCrystal().GetScatteringPowerRegistry().GetNb();++i)
         {
            ScatteringPower *s=&(this->GetCrystal().GetScatteringPowerRegistry().GetObj(i));
            if(s->GetClassName()=="ScatteringPowerAtom")
            {
               mvScattPowRowIndex.push_back(dynamic_cast<ScatteringPowerAtom *>(s));
               mvScattPowIndex[mvScattPowRowIndex.back()]=ix;
            }
         }
         mvScattPowIndexClock.Click();
      }
      if(mpScattPowWin!=0)
      {// Add/remove row(s), if necessary
         bool needLayout=false;
         for(unsigned long i=mvpRowScattPow.size();i<mvScattPowIndex.size();++i)
         {
            mpScattPowWin->AppendRows();
            mvpRowScattPow.push_back(RowScattPow());
            needLayout=true;
         }
         for(unsigned long i=mvScattPowIndex.size();i<mvpRowScattPow.size();++i)
         {
            mpScattPowWin->DeleteRows(i,1,false);
            mvpRowScattPow.pop_back();
            needLayout=true;
         }
         if(needLayout) mpScattPowWin->FitInside();
      }
      if(mpAntiBumpWin!=0)
      {// Add AntiBump rows, if necessary
         bool needLayout=false;
         for(unsigned long i=mvpRowAntiBump.size();i<mvScattPowIndex.size();++i)
         {
            mpAntiBumpWin->AppendRows();
            mpAntiBumpWin->AppendCols();
            mvpRowAntiBump.push_back(RowAntiBump());
            needLayout=true;
         }
         for(unsigned long i=mvScattPowIndex.size();i<mvpRowAntiBump.size();++i)
         {
            mpAntiBumpWin->DeleteRows(i,1,false);
            mpAntiBumpWin->DeleteCols(i,1,false);
            mvpRowAntiBump.pop_back();
            needLayout=true;
         }
         if(needLayout) mpAntiBumpWin->FitInside();
      }
      if(mpBondValenceWin!=0)
      {// Add Bond Valence row, if necessary
         bool needLayout=false;
         for(unsigned long i=mvpRowBondValence.size();i<mvScattPowIndex.size();++i)
         {
            mpBondValenceWin->AppendRows();
            mpBondValenceWin->AppendCols();
            mvpRowBondValence.push_back(RowBondValence());
            needLayout=true;
         }
         for(unsigned long i=mvScattPowIndex.size();i<mvpRowBondValence.size();++i)
         {
            mpBondValenceWin->DeleteRows(i,1,false);
            mpBondValenceWin->DeleteCols(i,1,false);
            mvpRowBondValence.pop_back();
            needLayout=true;
         }
         if(needLayout) mpBondValenceWin->FitInside();
      }
      
      if(mpScattPowWin!=0)
      {
         list<RowScattPow>::iterator row=mvpRowScattPow.begin();
         map<ScatteringPowerAtom*,int>::iterator pow=mvScattPowIndex.begin();
         for(;row!=mvpRowScattPow.end();)
         {
            const string name=pow->first->GetName();
            const REAL biso=pow->first->GetBiso();
            const REAL formalCharge=pow->first->GetFormalCharge();
            const REAL *pRGB=pow->first->GetColourRGB();
            const REAL mlerror=pow->first->GetMaximumLikelihoodPositionError();
            const REAL nbghost=pow->first->GetMaximumLikelihoodNbGhostAtom();
            if(  (name   !=row->mName)
               ||(biso!=row->mBiso)
               ||(formalCharge!=row->mFormalCharge)
               ||(pRGB[0]!=row->mR)
               ||(pRGB[1]!=row->mG)
               ||(pRGB[2]!=row->mB)
               ||(mlerror!=row->mMaximumLikelihoodError)
               ||(nbghost!=row->mNbGhostAtoms))
            {
               row->mName=name;
               row->mBiso=biso;
               row->mFormalCharge=formalCharge;
               row->mR=pRGB[0];
               row->mG=pRGB[1];
               row->mB=pRGB[2];
               row->mMaximumLikelihoodError=mlerror;
               row->mNbGhostAtoms=nbghost;
               row->mNeedUpdateUI=true;
            }
            ++row;++pow;
         }
      }
      if(mpAntiBumpWin!=0)
      {
         list<RowAntiBump>::iterator row=mvpRowAntiBump.begin();
         map<ScatteringPowerAtom*,int>::iterator pow=mvScattPowIndex.begin();
         for(;row!=mvpRowAntiBump.end();)
         {
            const string name=pow->first->GetName();
            const Crystal::VBumpMergePar *pMap=&(mpCrystal->GetBumpMergeParList());
            vector<REAL> dist(mvScattPowIndex.size());
            for(unsigned int i=0;i<mvScattPowRowIndex.size();++i)
            {
               Crystal::VBumpMergePar::const_iterator pos;
               if(pow->first<mvScattPowRowIndex[i]) pos=pMap->find(make_pair(pow->first,mvScattPowRowIndex[i]));
               else pos=pMap->find(make_pair(mvScattPowRowIndex[i],pow->first));
               if(pos==pMap->end()) dist[i]=-999;
               else dist[i]=sqrt(pos->second.mDist2);
            }
            if(  (name!=row->mName)
               ||(dist!=row->mvAntiBumpDistance))
            {
               row->mName=name;
               row->mvAntiBumpDistance=dist;
               row->mNeedUpdateUI=true;
            }
            ++row;++pow;
         }
      }
      if(mpBondValenceWin!=0)
      {
         list<RowBondValence>::iterator row=mvpRowBondValence.begin();
         map<ScatteringPowerAtom*,int>::iterator pow=mvScattPowIndex.begin();
         for(;row!=mvpRowBondValence.end();)
         {
            const string name=pow->first->GetName();
            const std::map<pair<const ScatteringPower*,const ScatteringPower*>, REAL> *pMap=&(mpCrystal->GetBondValenceRoList());
            vector<REAL> ro(mvScattPowIndex.size());
            for(unsigned int i=0;i<mvScattPowRowIndex.size();++i)
            {
               map<pair<const ScatteringPower*,const ScatteringPower*>, REAL>::const_iterator pos;
               if(pow->first<mvScattPowRowIndex[i]) pos=pMap->find(make_pair(pow->first,mvScattPowRowIndex[i]));
               else pos=pMap->find(make_pair(mvScattPowRowIndex[i],pow->first));
               if(pos==pMap->end()) ro[i]=-999;
               else ro[i]=pos->second;
            }
            if(  (name!=row->mName)
               ||(ro!=row->mvBondValenceRo))
            {
               row->mName=name;
               row->mvBondValenceRo=ro;
               row->mNeedUpdateUI=true;
            }
            ++row;++pow;
         }
      }
   }
   if(lock) mMutex.Unlock();

   this->WXRefinableObj::CrystUpdate(uui,lock);
   VFN_DEBUG_EXIT("WXCrystal::CrystUpdate():End",5)
}

#ifdef OBJCRYST_GL
void WXCrystal::UpdateGL(const bool onlyIndependentAtoms,
                         const REAL xMin,const REAL xMax,
                         const REAL yMin,const REAL yMax,
                         const REAL zMin,const REAL zMax)
{
   // :KLUDGE: !!! UGLY !!! This should be done in WXGLCrystalCanvas !
   VFN_DEBUG_ENTRY("WXCrystal::UpdateGL()",8)
   WXCrystValidateAllUserInput();
   if(mpCrystalGL!=0)
   {
      VFN_DEBUG_MESSAGE("WXCrystal::UpdateGL():mpCrystalGL",7)
      
      if(false==wxThread::IsMain())
      {
         mpConditionGLUpdate=new wxCondition(mMutexGLUpdate);
         mMutexGLUpdate.Lock();
         wxCommandEvent event(wxEVT_COMMAND_MENU_SELECTED,ID_GLCRYSTAL_MENU_UPDATE);
         wxPostEvent(mpCrystalGL,event);
         mpConditionGLUpdate->Wait();
         mMutexGLUpdate.Unlock();
         delete mpConditionGLUpdate;
         mpConditionGLUpdate=0;
         VFN_DEBUG_EXIT("WXCrystal::UpdateGL()-Not in main thread :End",8)
         return;
      }
      if(mCrystalGLDisplayList==0)
      {
         mCrystalGLDisplayList=glGenLists(1);
         mCrystalGLNameDisplayList=glGenLists(1);
         VFN_DEBUG_MESSAGE("WXCrystal::UpdateGL():created mCrystalGLDisplayList="<<mCrystalGLDisplayList,7)
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
      if(mpConditionGLUpdate!=0)
      {
         wxMutexLocker lock(mMutexGLUpdate);
         mpConditionGLUpdate->Signal();
      }
      mpCrystalGL->CrystUpdate();
   }
   else
   {
      VFN_DEBUG_MESSAGE("WXCrystal::UpdateGL():No mpCrystalGL",7)
   }
   VFN_DEBUG_EXIT("WXCrystal::UpdateGL():End",8)
}

int WXCrystal::GetCrystalGLDisplayList(const bool atomName)const
{
   VFN_DEBUG_MESSAGE("WXCrystal::GetCrystalGLDisplayList()",7)
   if(atomName) return mCrystalGLNameDisplayList;
   return mCrystalGLDisplayList;
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
   if(mpCrystalGL!=0)
   {
      BBox box=mpCrystalGL->GetCellBBox();
      this->UpdateGL(false,box.xMin,box.xMax,box.yMin,box.yMax,box.zMin,box.zMax);
   }
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
   WXCrystValidateAllUserInput();
   ScatteringPowerAtom *scatt=new ScatteringPowerAtom("Change me","H");
   mpCrystal->AddScatteringPower(scatt);
   VFN_DEBUG_MESSAGE("WXCrystal::OnMenuAddScattPowAtom():End",6)
   this->Layout();
   this->CrystUpdate(true,false);
}

void WXCrystal::OnMenuAddScattPowSphere(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_ENTRY("WXCrystal::OnMenuAddScattSphere()",6)
   WXCrystValidateAllUserInput();
   ScatteringPower *scatt= new ScatteringPowerSphere;
   mpCrystal->AddScatteringPower(scatt);
   this->Layout();
   this->CrystUpdate(true,false);
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
   const ScatteringComponentList *pList=&(mpCrystal->GetScatteringComponentList());
   for(long i=0;i<pList->GetNbComponent();++i)
      if((*pList)(i).mpScattPow==scatt)
      {
         wxMessageDialog dumbUser(this,"This Scattering Power is still used !",
                                  "Whooops",wxOK|wxICON_EXCLAMATION);
         dumbUser.ShowModal();
         VFN_DEBUG_EXIT("WXCrystal::OnButtonRemoveScattPow()",6)
         return;
      }
   mpCrystal->RemoveScatteringPower(scatt);
   VFN_DEBUG_EXIT("WXCrystal::OnButtonRemoveScattPow()",6)
   this->Layout();
   this->CrystUpdate(true,false);
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
   this->CrystUpdate(true);
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
   WXCrystValidateAllUserInput();
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
   this->CrystUpdate(true);
   VFN_DEBUG_EXIT("WXCrystal::OnMenuImportFenskeHallZMatrix()",6)
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
      this->CrystUpdate(true);
      this->Layout();
      return true;
   }
   return false;
}

void WXCrystal::UpdateUI(const bool lock)
{
   VFN_DEBUG_ENTRY("WXCrystal::UpdateUI()",6)
   if(!mpCrystal->IsBeingRefined())
   {
      if(lock) mMutex.Lock();
      mpFieldSpacegroup->SetValue(mpCrystal->GetSpaceGroup().GetName());
      #ifdef OBJCRYST_GL
      if(0!=mpCrystalGL) mpCrystalGL->GetParent()->SetTitle(mpCrystal->GetName().c_str());
      #endif
      if(lock) mMutex.Unlock();
   }
   if(lock) mMutex.Lock();
   if(0!=mpScattPowWin)
   {
      unsigned long i=0;
      for(list<RowScattPow>::iterator pos=mvpRowScattPow.begin();pos!=mvpRowScattPow.end();++pos)
      {
         if(pos->mNeedUpdateUI==true)
         {
            mIsSelfUpdating=true;
            mpScattPowWin->SetRowLabelValue(i,pos->mName.c_str());
            wxString tmp;
            tmp.Printf("%f",pos->mBiso);
            mpScattPowWin->SetCellValue(i, 0, tmp);
            tmp.Printf("%f",pos->mFormalCharge);
            mpScattPowWin->SetCellValue(i, 1, tmp);
            tmp.Printf("%f",pos->mR);
            mpScattPowWin->SetCellValue(i, 2, tmp);
            tmp.Printf("%f",pos->mG);
            mpScattPowWin->SetCellValue(i, 3, tmp);
            tmp.Printf("%f",pos->mB);
            mpScattPowWin->SetCellValue(i, 4, tmp);
            tmp.Printf("%f",pos->mMaximumLikelihoodError);
            mpScattPowWin->SetCellValue(i, 5, tmp);
            tmp.Printf("%f",pos->mNbGhostAtoms);
            mpScattPowWin->SetCellValue(i, 6, tmp);
            mIsSelfUpdating=false;
         }
         ++i;
      }
   }
   if(0!=mpAntiBumpWin)
   {
      unsigned long i=0;
      for(list<RowAntiBump>::iterator pos=mvpRowAntiBump.begin();pos!=mvpRowAntiBump.end();++pos)
      {
         if(pos->mNeedUpdateUI==true)
         {
            mIsSelfUpdating=true;
            mpAntiBumpWin->SetRowLabelValue(i,pos->mName.c_str());
            mpAntiBumpWin->SetColLabelValue(i,pos->mName.c_str());
            wxString tmp;
            for(unsigned long j=0;j<pos->mvAntiBumpDistance.size();++j)
            {
               VFN_DEBUG_MESSAGE("WXCrystal::UpdateUI():Antibump("<<mvScattPowRowIndex[i]->GetName()
                                                                  <<","<<mvScattPowRowIndex[j]->GetName()
                                                                  <<")="<<pos->mvAntiBumpDistance[j],3);
               if(pos->mvAntiBumpDistance[j]>-998)
               {
                  tmp.Printf("%f",pos->mvAntiBumpDistance[j]);
                  mpAntiBumpWin->SetCellValue(i,j,tmp);
               } else mpAntiBumpWin->SetCellValue(i,j,"");
            }
            mIsSelfUpdating=false;
         }
         ++i;
      }
   }
   if(0!=mpBondValenceWin)
   {
      unsigned long i=0;
      for(list<RowBondValence>::iterator pos=mvpRowBondValence.begin();pos!=mvpRowBondValence.end();++pos)
      {
         if(pos->mNeedUpdateUI==true)
         {
            mIsSelfUpdating=true;
            mpBondValenceWin->SetRowLabelValue(i,pos->mName.c_str());
            mpBondValenceWin->SetColLabelValue(i,pos->mName.c_str());
            wxString tmp;
            for(unsigned long j=0;j<pos->mvBondValenceRo.size();++j)
            {
               VFN_DEBUG_MESSAGE("WXCrystal::UpdateUI():BondValence("<<mvScattPowRowIndex[i]->GetName()
                                                                <<","<<mvScattPowRowIndex[j]->GetName()
                                                               <<")="<<pos->mvBondValenceRo[j],3);
               if(pos->mvBondValenceRo[j]>-998)
               {
                  tmp.Printf("%f",pos->mvBondValenceRo[j]);
                  mpBondValenceWin->SetCellValue(i,j,tmp);
               } else mpBondValenceWin->SetCellValue(i,j,"");
            }
            mIsSelfUpdating=false;
         }
         ++i;
      }
   }
   if(lock) mMutex.Unlock();
   this->WXRefinableObj::UpdateUI(lock);
   VFN_DEBUG_EXIT("WXCrystal::UpdateUI()",6)
}
Crystal& WXCrystal::GetCrystal(){return *mpCrystal;}
const Crystal& WXCrystal::GetCrystal()const{return *mpCrystal;}

void WXCrystal::OnMenuShowScattPowWindow(wxCommandEvent &event)
{
   VFN_DEBUG_MESSAGE("WXCrystal::OnMenuShowScattPowWindow()",10)
   if(0!=mpScattPowWin) return;
   WXCrystValidateAllUserInput();
   // Frame with notebook
      wxFrame *frame= new wxFrame(this,-1,("Scattering Powers parameters for: "
                                  +this->GetCrystal().GetName()).c_str(),
                                  wxDefaultPosition,wxSize(800,300));

      wxNotebook *notebook = new wxNotebook(frame, -1);
   {// Individual parameters
      mpScattPowWin = new WXCrystalScrolledGridWindow(notebook,this,ID_CRYSTAL_WIN_SCATTPOW);
      notebook->AddPage(mpScattPowWin, "Scattering Powers", true);
      
      mpScattPowWin->SetDefaultRenderer(new wxGridCellFloatRenderer(5,3));
      mpScattPowWin->SetDefaultEditor(new wxGridCellFloatEditor(5,3));
      mpScattPowWin->SetColMinimalAcceptableWidth(150);
      mpScattPowWin->CreateGrid(0,7);
      
      mpScattPowWin->SetColLabelValue(0,"Biso");
      mpScattPowWin->SetColLabelValue(1,"Charge");
      mpScattPowWin->SetColLabelValue(2,"Red");
      mpScattPowWin->SetColLabelValue(3,"Green");
      mpScattPowWin->SetColLabelValue(4,"Blue");
      mpScattPowWin->SetColLabelValue(5,"ML Error");
      mpScattPowWin->SetColLabelValue(6,"#ghost");
      
      mpScattPowWin->AutoSizeRows();
      mpScattPowWin->AutoSizeColumns();
   }
   {// Anti-Bump
      mpAntiBumpWin = new WXCrystalScrolledGridWindow(notebook,this,ID_CRYSTAL_WIN_ANTIBUMP);
      notebook->AddPage(mpAntiBumpWin, "AntiBump", true);
      
      mpAntiBumpWin->SetDefaultRenderer(new wxGridCellFloatRenderer(5,3));
      mpAntiBumpWin->SetDefaultEditor(new wxGridCellFloatEditor(5,3));
      mpAntiBumpWin->SetColMinimalAcceptableWidth(150);
      mpAntiBumpWin->CreateGrid(0,0);
      
      mpAntiBumpWin->AutoSizeRows();
      mpAntiBumpWin->AutoSizeColumns();
   }
   {// Bond Valence
      mpBondValenceWin = new WXCrystalScrolledGridWindow(notebook,this,ID_CRYSTAL_WIN_BONDVALENCE);
      notebook->AddPage(mpBondValenceWin, "BondValence", true);
      
      mpBondValenceWin->SetDefaultRenderer(new wxGridCellFloatRenderer(5,3));
      mpBondValenceWin->SetDefaultEditor(new wxGridCellFloatEditor(5,3));
      mpBondValenceWin->SetColMinimalAcceptableWidth(150);
      mpBondValenceWin->CreateGrid(0,0);
      
      mpBondValenceWin->AutoSizeRows();
      mpBondValenceWin->AutoSizeColumns();
   }
   notebook->SetSelection(0);
   this->CrystUpdate(true);
   frame->Show(true);
   frame->Layout();
}

void WXCrystal::OnEditGridScattPow(wxGridEvent &e)
{
   if(mIsSelfUpdating) return;
   const int r=e.GetRow();
   const int c=e.GetCol();
   ScatteringPowerAtom *const p=mvScattPowRowIndex[r];
   wxString s=mpScattPowWin->GetCellValue(r,c);
   switch(c)
   {
      case 0:
      {
         if(s!="")
         {
            double d;
            s.ToDouble(&d);
            mvScattPowRowIndex[r]->SetBiso(d);
         }
         break;
      }
      case 1:
      {
         if(s!="")
         {
            double d;
            s.ToDouble(&d);
            p->SetFormalCharge(d);
         }
         break;
      }
      case 2:
      {
         if(s!="")
         {
            double d;
            s.ToDouble(&d);
            const REAL gg=p->GetColourRGB()[1];
            const REAL bb=p->GetColourRGB()[2];
            p->SetColour(d,gg,bb);
         }
         break;
      }
      case 3:
      {
         if(s!="")
         {
            double d;
            s.ToDouble(&d);
            const REAL rr=p->GetColourRGB()[0];
            const REAL bb=p->GetColourRGB()[2];
            p->SetColour(rr,d,bb);
         }
         break;
      }
      case 4:
      {
         if(s!="")
         {
            double d;
            s.ToDouble(&d);
            const REAL rr=p->GetColourRGB()[0];
            const REAL gg=p->GetColourRGB()[1];
            p->SetColour(rr,gg,d);
         }
         break;
      }
      case 5:
      {
         if(s!="")
         {
            double d;
            s.ToDouble(&d);
            p->SetMaximumLikelihoodPositionError(d);
         }
         break;
      }
      case 6:
      {
         if(s!="")
         {
            double d;
            s.ToDouble(&d);
            p->SetMaximumLikelihoodNbGhostAtom(d);
         }
         break;
      }
   }
   this->CrystUpdate();
}

void WXCrystal::OnEditGridScattPowAntiBump(wxGridEvent &e)
{
   if(mIsSelfUpdating) return;
   const int r=e.GetRow();
   const int c=e.GetCol();
   const ScatteringPowerAtom *const p1=mvScattPowRowIndex[r];
   const ScatteringPowerAtom *const p2=mvScattPowRowIndex[c];
   wxString s=mpAntiBumpWin->GetCellValue(r,c);
   double d;
   s.ToDouble(&d);
   if(d>0.01) mpCrystal->SetBumpMergeDistance(*p1,*p2,d);
   else mpCrystal->RemoveBumpMergeDistance(*p1,*p2);
   this->CrystUpdate(true,false);
}

void WXCrystal::OnEditGridScattPowBondValence(wxGridEvent &e)
{
   if(mIsSelfUpdating) return;
   const int r=e.GetRow();
   const int c=e.GetCol();
   const ScatteringPowerAtom *const p1=mvScattPowRowIndex[r];
   const ScatteringPowerAtom *const p2=mvScattPowRowIndex[c];
   wxString s=mpBondValenceWin->GetCellValue(r,c);
   double d;
   s.ToDouble(&d);
   if(d>0.01) mpCrystal->AddBondValenceRo(*p1,*p2,d);
   else mpCrystal->RemoveBondValenceRo(*p1,*p2);
   this->CrystUpdate(true,false);
}

void WXCrystal::NotifyDeleteListWin(WXCrystalScrolledGridWindow *win)
{
   if(win==mpScattPowWin)
   {
      mpScattPowWin=0;
      mvpRowScattPow.clear();
   }
   if(win==mpAntiBumpWin)
   {
      mpAntiBumpWin=0;
      mvpRowAntiBump.clear();
   }
   if(win==mpBondValenceWin)
   {
      mpBondValenceWin=0;
      mvpRowBondValence.clear();
   }
}
bool WXCrystal::Enable(bool e)
{
   if(0!=mpScattPowWin)    mpScattPowWin   ->Enable(e);
   if(0!=mpAntiBumpWin)    mpAntiBumpWin   ->Enable(e);
   if(0!=mpBondValenceWin) mpBondValenceWin->Enable(e);
   return this->::wxWindow::Enable(e);
}

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
   cout<<"Generating OpenGL Triangles for Fourier map:"<<mName<<", contour="<<minValue<<endl;
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
            if(minValue>0)
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
            else
               for(int j=2; j >=0; j--)
               {
                  //VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnPaint():MC1:"<<i<<" "<<j,5)
                  //:TODO: Fix normals
                  normx=-pTriangles[i].norm[j].x;
                  normy=-pTriangles[i].norm[j].y;
                  normz=-pTriangles[i].norm[j].z;
                  //mpCrystal->FractionalToOrthonormalCoords(normx, normy, normz);
                  //mpCrystal->OrthonormalToFractionalCoords(normx, normy, normz);
                  glNormal3f(normx, normy, normz);
                  glVertex3f(pTriangles[i].p[j].x    ,pTriangles[i].p[j].y    ,pTriangles[i].p[j].z);
               }
         }
      glEnd();

   delete [] subPoints; 
   delete [] pTriangles; 
   VFN_DEBUG_EXIT("UnitCellMapImport::GLInitDisplayList()",7)
}

void UnitCellMapImport::POVRayDescription(ostream &os,const float minValue,
                                          const CrystalPOVRayOptions &options)const
{// basically the same code asGLInitDisplayList(), but creates cylinders
   VFN_DEBUG_ENTRY("UnitCellMapImport::POVRayDescription()",7)
   // Generate triangles
      VFN_DEBUG_MESSAGE("UnitCellMapImport::POVRayDescription(): Generate Triangles",7)

      const int nx=mPoints.cols();
      const int ny=mPoints.rows();
      const int nz=mPoints.depth();
      float step[3];
      step[0]=1/(float)nx;
      step[1]=1/(float)ny;
      step[2]=1/(float)nz;
      int nxMin, nxMax, nyMin, nyMax, nzMin, nzMax;
      nxMin = (int)(options.mXmin * nx);
      nxMax = (int)(options.mXmax * nx);
      nyMin = (int)(options.mYmin * ny);
      nyMax = (int)(options.mYmax * ny);
      nzMin = (int)(options.mZmin * nz);
      nzMax = (int)(options.mZmax * nz);
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
      VFN_DEBUG_MESSAGE("UnitCellMapImport::POVRayDescription(): MC, Min Value="<<minValue,10)
      const TRIANGLE *pTriangles= MC(snx-1, sny-1, snz-1, step[0], step[1], step[2], minValue, subPoints, numOfTriangles);
   // drawing instructions
      VFN_DEBUG_MESSAGE("UnitCellMapImport::POVRayDescription(): POVRay instructions",7)
      float normx,normy,normz;
      for(int i=0; i < numOfTriangles; i++)
      {
         const float x1=pTriangles[i].p[0].x;
         const float x2=pTriangles[i].p[1].x;
         const float x3=pTriangles[i].p[2].x;
         const float y1=pTriangles[i].p[0].y;
         const float y2=pTriangles[i].p[1].y;
         const float y3=pTriangles[i].p[2].y;
         const float z1=pTriangles[i].p[0].z;
         const float z2=pTriangles[i].p[1].z;
         const float z3=pTriangles[i].p[2].z;
         
         // Avoid null-length cylinders that make POV-Ray choke
         const float d12=abs(x1-x2)+abs(y1-y2)+abs(z1-z2);
         const float d13=abs(x1-x3)+abs(y1-y3)+abs(z1-z3);
         const float d23=abs(x2-x3)+abs(y2-y3)+abs(z2-z3);
         if((d12<0.05)||(d13<0.05)||(d23<0.05)) continue;
         
         //:TODO: Fix normals
         normx=pTriangles[i].norm[j].x;
         normy=pTriangles[i].norm[j].y;
         normz=pTriangles[i].norm[j].z;
         //mpCrystal->FractionalToOrthonormalCoords(normx, normy, normz);
         //mpCrystal->OrthonormalToFractionalCoords(normx, normy, normz);
         os<<"      ObjCrystMeshTriangle("
           <<x1<<","<<y1<<","<<z1<<","
           <<x2<<","<<y2<<","<<z2<<","
           <<x3<<","<<y3<<","<<z3<<","
           <<normx<<","<<normy<<","<<normz<<","
           <<normx<<","<<normy<<","<<normz<<","
           <<normx<<","<<normy<<","<<normz<<")"
           <<endl;
      }

   delete [] subPoints; 
   delete [] pTriangles; 
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
   
   mMean=mPoints.sum()/(REAL)(mPoints.numElements());
   mMin=mPoints.min();
   mMax=mPoints.max();
   {
      mStandardDeviation=0.0;
      const REAL *tmp=mPoints.data();
      for(long i=0;i<mPoints.numElements();i++)
      {
         mStandardDeviation += (*tmp-mMean) * (*tmp-mMean);
         tmp++;
      }
      mStandardDeviation = sqrt(mStandardDeviation/(REAL)(mPoints.numElements()));
   }
   cout << "Min density value="<<mMin<<endl
        << "Max density value="<<mMax<<endl
        << "Mean density="<<mMean<<endl
        << "Standard Deviation="<<mStandardDeviation<<endl;
   
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

/// Byte-swapping function for DSN6 import
void swap2(void *data, unsigned int nb)
{
   char * dataptr = (char *)data;
   char tmp;

   for (unsigned int i=0; i<(nb-1); i+=2)
   {
      tmp = dataptr[i];
      dataptr[i] = dataptr[i+1];
      dataptr[i+1] = tmp;
   }
}

int UnitCellMapImport::ImportDSN6(const string&filename)
{
   VFN_DEBUG_ENTRY("UnitCellMapImport::ImportDSN6()",7)
   FILE *pfile=fopen(filename.c_str(),"rb");
   if(NULL==pfile)
   {     //if file could not be loaded for some reason then exit
     VFN_DEBUG_MESSAGE("UnitCellMapImport::ImportDSN6() error opening "<<filename.c_str(),10)
      (*fpObjCrystInformUser)("Error opening file: "+filename);
      return 0;
   }
   // :KLUDGE: assume sizeof(short int)==2...
   short header[256];
   fread(header, sizeof(short), 256, pfile);
   bool needswap=false;
   if (header[18] == 25600) needswap=true;
   if(needswap) swap2(header, 19*sizeof(short));
   
   const long xstart=header[0];
   const long ystart=header[1];
   const long zstart=header[2];
   const long xextent=header[3];
   const long yextent=header[4];
   const long zextent=header[5];
   const unsigned long xsamplingrate=header[6];
   const unsigned long ysamplingrate=header[7];
   const unsigned long zsamplingrate=header[8];
   const float celledgea=(float)header[ 9]/(float)header[17];
   const float celledgeb=(float)header[10]/(float)header[17];
   const float celledgec=(float)header[11]/(float)header[17];
   const float alpha=(float)header[12]/(float)header[17];
   const float beta=(float)header[13]/(float)header[17];
   const float gamma=(float)header[14]/(float)header[17];
   const float rhoscale=(float)header[15]/(float)header[18];
   const float rhozero =(float)header[16];
   cout <<"xstart="<<xstart<<endl
        <<"ystart="<<ystart<<endl
        <<"zstart="<<zstart<<endl
        <<"xextent="<<xextent<<endl
        <<"yextent="<<yextent<<endl
        <<"zextent="<<zextent<<endl
        <<"xsamplingrate="<<xsamplingrate<<endl
        <<"ysamplingrate="<<ysamplingrate<<endl
        <<"zsamplingrate="<<zsamplingrate<<endl
        <<"celledgea="<<celledgea<<endl
        <<"celledgeb="<<celledgeb<<endl
        <<"celledgec="<<celledgec<<endl
        <<"alpha="<<alpha<<endl
        <<"beta="<<beta<<endl
        <<"gamma="<<gamma<<endl
        <<"rhoscale="<<rhoscale<<endl
        <<"rhozero="<<rhozero<<endl;

   #define BRICKSIZE 512
   #define BRICKEDGE 8
   
   unsigned char Points[BRICKSIZE];
   // :KLUDGE: stored map is the size of the entire unit cell, not of the input map...
   // The result is that if the input map is smaller than the unit cell size,
   // the rest of the cell will be with zero density
   mPoints.resize(xsamplingrate,ysamplingrate,zsamplingrate);
   mPoints=0;

   const unsigned int nxbrick=((xextent)/BRICKEDGE)+(xextent%8 ? 1 : 0);
   const unsigned int nybrick=((yextent)/BRICKEDGE)+(zextent%8 ? 1 : 0);
   const unsigned int nzbrick=((zextent)/BRICKEDGE)+(zextent%8 ? 1 : 0);

   for(unsigned int zbrick=0;zbrick<nzbrick;zbrick++)
      for(unsigned int ybrick=0;ybrick<nybrick;ybrick++)
         for(unsigned int xbrick=0;xbrick<nxbrick;xbrick++)
         {
            fread(&Points, sizeof(unsigned char), BRICKSIZE, pfile);
            // In this SICK format, data is stored as bytes, but which are
            // nevertheless byteswapped as 2-byte words. Bizarre, vous avez dit bizarre...
            if(needswap) swap2((void *)&Points, BRICKSIZE);
            const unsigned char* pPoint=&Points[0];
            for(unsigned int z=0;z<BRICKEDGE;z++)
               for(unsigned int y=0;y<BRICKEDGE;y++)
                  for(unsigned int x=0;x<BRICKEDGE;x++)
                  {
                     if(  ((xbrick*BRICKEDGE+x)<xsamplingrate)
                        &&((ybrick*BRICKEDGE+y)<ysamplingrate)
                        &&((zbrick*BRICKEDGE+z)<zsamplingrate))
                        mPoints(zbrick*BRICKEDGE+z,ybrick*BRICKEDGE+y,xbrick*BRICKEDGE+x)
                           = ((float) *pPoint - rhozero)/ rhoscale ;
                        pPoint++;
                  }
         }
   fclose (pfile);
   
   mMean=mPoints.sum()/(REAL)(mPoints.numElements());
   mMin=mPoints.min();
   mMax=mPoints.max();
   {
      mStandardDeviation=0.0;
      const REAL *tmp=mPoints.data();
      for(long i=0;i<mPoints.numElements();i++)
      {
         mStandardDeviation += (*tmp-mMean) * (*tmp-mMean);
         tmp++;
      }
      mStandardDeviation = sqrt(mStandardDeviation/(REAL)(mPoints.numElements()));
   }
   cout << "Min density value="<<mMin<<endl
        << "Max density value="<<mMax<<endl
        << "Mean density="<<mMean<<endl
        << "Standard Deviation="<<mStandardDeviation<<endl;
   
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
   VFN_DEBUG_EXIT("UnitCellMapImport::ImportDSN6()",7)
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
REAL UnitCellMapImport::Max()const{return mMax;}
REAL UnitCellMapImport::Min()const{return mMin;}
REAL UnitCellMapImport::Mean()const{return mMean;}
REAL UnitCellMapImport::StandardDeviation()const{return mStandardDeviation;}

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

bool UnitCellMapGLList::ShowWire()const
{
   return mShowWire;
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
static const long ID_GLCRYSTAL_MENU_LIMITS_FULLCELL=   WXCRYST_ID(); 
static const long ID_GLCRYSTAL_MENU_LIMITS_ASYMCELL=   WXCRYST_ID(); 
static const long ID_GLCRYSTAL_MENU_SHOWCRYSTAL=       WXCRYST_ID(); 
static const long ID_GLCRYSTAL_MENU_LOADFOURIERGRD=    WXCRYST_ID(); 
static const long ID_GLCRYSTAL_MENU_LOADFOURIERDSN6=   WXCRYST_ID(); 
static const long ID_GLCRYSTAL_MENU_CHANGECONTOUR=     WXCRYST_ID(); 
static const long ID_GLCRYSTAL_MENU_ADDCONTOUR=        WXCRYST_ID(); 
static const long ID_GLCRYSTAL_MENU_SHOWFOURIER=       WXCRYST_ID(); 
static const long ID_GLCRYSTAL_MENU_FOURIERCHANGECOLOR=WXCRYST_ID(); 
static const long ID_GLCRYSTAL_MENU_SHOWWIRE=          WXCRYST_ID(); 
static const long ID_GLCRYSTAL_MENU_UNLOADFOURIER=     WXCRYST_ID(); 
static const long ID_GLCRYSTAL_MENU_FOURIERCHANGEBBOX= WXCRYST_ID(); 
static const long ID_GLCRYSTAL_MENU_POVRAY=            WXCRYST_ID(); 

BEGIN_EVENT_TABLE(WXGLCrystalCanvas, wxGLCanvas)
   EVT_SIZE             (WXGLCrystalCanvas::OnSize)
   EVT_PAINT            (WXGLCrystalCanvas::OnPaint)
   EVT_ERASE_BACKGROUND (WXGLCrystalCanvas::OnEraseBackground)
   EVT_MOUSE_EVENTS     (WXGLCrystalCanvas::OnMouse)
   EVT_MENU             (ID_GLCRYSTAL_MENU_UPDATE,              WXGLCrystalCanvas::OnUpdate)
   EVT_MENU             (ID_GLCRYSTAL_MENU_CHANGELIMITS,        WXGLCrystalCanvas::OnChangeLimits)
   EVT_MENU             (ID_GLCRYSTAL_MENU_LIMITS_FULLCELL,     WXGLCrystalCanvas::OnChangeLimits)
   EVT_MENU             (ID_GLCRYSTAL_MENU_LIMITS_ASYMCELL,     WXGLCrystalCanvas::OnChangeLimits)
   EVT_MENU             (ID_GLCRYSTAL_MENU_SHOWCRYSTAL,         WXGLCrystalCanvas::OnShowCrystal)
   EVT_MENU             (ID_GLCRYSTAL_MENU_SHOWATOMLABEL,       WXGLCrystalCanvas::OnShowAtomLabel)
   EVT_MENU             (ID_GLCRYSTAL_MENU_SHOWCURSOR,          WXGLCrystalCanvas::OnShowCursor)
   EVT_MENU             (ID_GLCRYSTAL_MENU_SETCURSOR,           WXGLCrystalCanvas::OnSetCursor)
   EVT_MENU             (ID_GLCRYSTAL_MENU_LOADFOURIERGRD,      WXGLCrystalCanvas::OnLoadFourierGRD)
   EVT_MENU             (ID_GLCRYSTAL_MENU_LOADFOURIERDSN6,     WXGLCrystalCanvas::OnLoadFourierDSN6)
   EVT_MENU             (ID_GLCRYSTAL_MENU_CHANGECONTOUR,       WXGLCrystalCanvas::OnChangeContour)
   EVT_MENU             (ID_GLCRYSTAL_MENU_ADDCONTOUR,          WXGLCrystalCanvas::OnAddContour)
   EVT_MENU             (ID_GLCRYSTAL_MENU_SHOWFOURIER,         WXGLCrystalCanvas::OnShowFourier)
   EVT_MENU             (ID_GLCRYSTAL_MENU_FOURIERCHANGECOLOR,  WXGLCrystalCanvas::OnFourierChangeColor)
   EVT_MENU             (ID_GLCRYSTAL_MENU_SHOWWIRE,            WXGLCrystalCanvas::OnShowWire)
   EVT_MENU             (ID_GLCRYSTAL_MENU_UNLOADFOURIER,       WXGLCrystalCanvas::OnUnloadFourier)
   EVT_MENU             (ID_GLCRYSTAL_MENU_FOURIERCHANGEBBOX,   WXGLCrystalCanvas::OnFourierChangeBbox)
   EVT_MENU             (ID_GLCRYSTAL_MENU_POVRAY,              WXGLCrystalCanvas::OnPOVRay)
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
   mcellbbox.xMin = mpWXCrystal->GetCrystal().GetSpaceGroup().GetAsymUnit().Xmin()-0.1;
   mcellbbox.yMin = mpWXCrystal->GetCrystal().GetSpaceGroup().GetAsymUnit().Ymin()-0.1;
   mcellbbox.zMin = mpWXCrystal->GetCrystal().GetSpaceGroup().GetAsymUnit().Zmin()-0.1;
   mcellbbox.xMax = mpWXCrystal->GetCrystal().GetSpaceGroup().GetAsymUnit().Xmax()+0.1;
   mcellbbox.yMax = mpWXCrystal->GetCrystal().GetSpaceGroup().GetAsymUnit().Ymax()+0.1;
   mcellbbox.zMax = mpWXCrystal->GetCrystal().GetSpaceGroup().GetAsymUnit().Zmax()+0.1;
     // N.B. xMin=xMax so that the previous cell bbox is used for Maps 
     // until mmapbbox is changed
   mmapbbox.xMin = mmapbbox.xMax = mmapbbox.yMin = mmapbbox.zMin = 0.;
   mmapbbox.yMax = mmapbbox.zMax = 1.;
   mpPopUpMenu=new wxMenu("Crystal");
   mpPopUpMenu->Append(ID_GLCRYSTAL_MENU_UPDATE, "&Update");
   mpPopUpMenu->AppendSeparator();
   mpPopUpMenu->Append(ID_GLCRYSTAL_MENU_CHANGELIMITS, "Change display &Limits");
   mpPopUpMenu->Append(ID_GLCRYSTAL_MENU_LIMITS_FULLCELL, "Show Full Unit Cell +0.1");
   mpPopUpMenu->Append(ID_GLCRYSTAL_MENU_LIMITS_ASYMCELL, "Show Asymmetric Unit Cell +0.1");
   mpPopUpMenu->AppendSeparator();
   mpPopUpMenu->Append(ID_GLCRYSTAL_MENU_SHOWCRYSTAL, "Hide Crystal");
   mpPopUpMenu->Append(ID_GLCRYSTAL_MENU_SHOWATOMLABEL, "Hide Atom Labels");
   mpPopUpMenu->Append(ID_GLCRYSTAL_MENU_SHOWCURSOR, "Show Cursor");
   mpPopUpMenu->Append(ID_GLCRYSTAL_MENU_SETCURSOR, "Set view cntr and cursor pos.");
   mpPopUpMenu->Append(ID_GLCRYSTAL_MENU_POVRAY, "Create POVRay file");
   mpPopUpMenu->AppendSeparator();
   mpPopUpMenu->Append(ID_GLCRYSTAL_MENU_LOADFOURIERGRD, "Load GRD Fourier Map");	
   mpPopUpMenu->Append(ID_GLCRYSTAL_MENU_LOADFOURIERDSN6,"Load DSN6 Fourier Map");	
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
   //Show display limits
   {
      int w, h;
      GetClientSize(& w, & h);
      glPushMatrix();
         glLoadIdentity();
         glMatrixMode(GL_PROJECTION);
         glPushMatrix();
            glLoadIdentity();
            gluOrtho2D(0,w,0,h);
            glColor3f(1.0,1.0,1.0);
            glRasterPos2i(2,h-12);
            char c[128];
            sprintf(c,"%5.3f<x<%5.3f\n",mcellbbox.xMin,mcellbbox.xMax);
            crystGLPrint(c);
            
            glRasterPos2i(2, h-24);
            sprintf(c,"%5.3f<y<%5.3f\n",mcellbbox.yMin,mcellbbox.yMax);
            crystGLPrint(c);
            
            glRasterPos2i(2, h-36);
            sprintf(c,"%5.3f<z<%5.3f\n",mcellbbox.zMin,mcellbbox.zMax);
            crystGLPrint(c);
         glPopMatrix();
         glMatrixMode( GL_MODELVIEW );
      glPopMatrix();
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
      glCallList(mpWXCrystal->GetCrystalGLDisplayList());  //Draw Crystal
      if(mShowAtomName)
      {
         glLoadIdentity();
         glColor3f(1.0,1.0,1.0);
         glTranslatef( -0.3, 0, -mDist+1. );// Put labels in front of the atom position
         glMultMatrixf( &m[0][0] );
         glTranslatef( mX0, mY0, mZ0 );
         glCallList(mpWXCrystal->GetCrystalGLDisplayList(true));  //Draw Atom Names
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
      x=(mcellbbox.xMax+mcellbbox.xMin)/2.-x;
      y=(mcellbbox.yMax+mcellbbox.yMin)/2.-y;
      z=(mcellbbox.zMax+mcellbbox.zMin)/2.-z;
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
   this->Refresh(false);
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
   VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnEnterWindow()",5)
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
      if(mpWXCrystal->GetCrystal().IsBeingRefined())
      {
         mpPopUpMenu->Enable(ID_GLCRYSTAL_MENU_UPDATE, false);
         mpPopUpMenu->Enable(ID_GLCRYSTAL_MENU_POVRAY, false);
         mpPopUpMenu->Enable(ID_GLCRYSTAL_MENU_LOADFOURIERGRD, false);
         mpPopUpMenu->Enable(ID_GLCRYSTAL_MENU_LOADFOURIERDSN6, false);
      }
      else
      {
         mpPopUpMenu->Enable(ID_GLCRYSTAL_MENU_UPDATE, true);
         mpPopUpMenu->Enable(ID_GLCRYSTAL_MENU_POVRAY, true);
         mpPopUpMenu->Enable(ID_GLCRYSTAL_MENU_LOADFOURIERGRD, true);
         mpPopUpMenu->Enable(ID_GLCRYSTAL_MENU_LOADFOURIERDSN6, true);
      }

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
   
   #ifdef HAVE_GLUT
   static bool needglutinit=true;
   if(needglutinit)
   {
      needglutinit=false;
      glutInit(&(wxApp::GetInstance()->argc),wxApp::GetInstance()->argv);
   }
   #endif

   wxSizeEvent event;
   this->OnSize(event);
   
   //First display
   this->CrystUpdate();
   VFN_DEBUG_EXIT("WXGLCrystalCanvas::InitGL()",8)
}
void WXGLCrystalCanvas::OnChangeLimits(wxCommandEvent &event)
{
  VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnChangeLimits()",10)
   if(event.GetId()==ID_GLCRYSTAL_MENU_LIMITS_FULLCELL)
   {
      mcellbbox.xMin = -0.1;
      mcellbbox.yMin = -0.1;
      mcellbbox.zMin = -0.1;
      mcellbbox.xMax =  1.1;
      mcellbbox.yMax =  1.1;
      mcellbbox.zMax =  1.1;
   }
   if(event.GetId()==ID_GLCRYSTAL_MENU_LIMITS_ASYMCELL)
   {
      mcellbbox.xMin = mpWXCrystal->GetCrystal().GetSpaceGroup().GetAsymUnit().Xmin()-0.1;
      mcellbbox.yMin = mpWXCrystal->GetCrystal().GetSpaceGroup().GetAsymUnit().Ymin()-0.1;
      mcellbbox.zMin = mpWXCrystal->GetCrystal().GetSpaceGroup().GetAsymUnit().Zmin()-0.1;
      mcellbbox.xMax = mpWXCrystal->GetCrystal().GetSpaceGroup().GetAsymUnit().Xmax()+0.1;
      mcellbbox.yMax = mpWXCrystal->GetCrystal().GetSpaceGroup().GetAsymUnit().Ymax()+0.1;
      mcellbbox.zMax = mpWXCrystal->GetCrystal().GetSpaceGroup().GetAsymUnit().Zmax()+0.1;
   }
   if(event.GetId()==ID_GLCRYSTAL_MENU_CHANGELIMITS)
   {

      UserSelectBoundingBox *BoxDlg = new UserSelectBoundingBox(this,
	       "Set bounding box for display of\natoms (fractional coordinates)",
			    mcellbbox);
      if (BoxDlg->ShowModal() == wxID_OK )
      {
         mcellbbox =  BoxDlg->GetBBox();
         mpWXCrystal->UpdateGL(false,
			    mcellbbox.xMin,mcellbbox.xMax,
			    mcellbbox.yMin,mcellbbox.yMax,
			    mcellbbox.zMin,mcellbbox.zMax);
         vector<pair<pair<const UnitCellMapImport*,float>,UnitCellMapGLList* > >::iterator pos;
         for(pos=mvpUnitCellMapGLList.begin();pos != mvpUnitCellMapGLList.end();pos++)
         {
            wxBusyInfo wait("Processing Fourier Map...");
            pos->second->GenList(*(pos->first.first),this, pos->first.second);
         }
         VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnChangeLimits (X: " << 
		         mcellbbox.xMin << ", " << mcellbbox.xMax << 
		         " Y: " << 
		         mcellbbox.yMin << ", " << mcellbbox.yMax << 
		         " Z: " << 
		         mcellbbox.zMin << ", " << mcellbbox.zMax << 
		         ")", 10)
      } 
      BoxDlg->Destroy();
   }
   if(!(mpWXCrystal->GetCrystal().IsBeingRefined()))
      mpWXCrystal->UpdateGL(false,
			                   mcellbbox.xMin,mcellbbox.xMax,
			                   mcellbbox.yMin,mcellbbox.yMax,
			                   mcellbbox.zMin,mcellbbox.zMax);

  VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnChangeLimits():UserSelectBoundingBox done",10)
}

void WXGLCrystalCanvas::OnShowCrystal( wxCommandEvent & WXUNUSED(event))
{
   if(mShowCrystal) mpPopUpMenu->SetLabel(ID_GLCRYSTAL_MENU_SHOWCRYSTAL, "Show Crystal");
   else mpPopUpMenu->SetLabel(ID_GLCRYSTAL_MENU_SHOWCRYSTAL, "Hide Crystal");
   mShowCrystal = !mShowCrystal;
   this->CrystUpdate();
}

void WXGLCrystalCanvas::OnShowAtomLabel( wxCommandEvent & WXUNUSED(event))
{
   if(mShowAtomName) mpPopUpMenu->SetLabel(ID_GLCRYSTAL_MENU_SHOWATOMLABEL, "Show Atom Labels");
   else mpPopUpMenu->SetLabel(ID_GLCRYSTAL_MENU_SHOWATOMLABEL, "Hide Atom Labels");
   mShowAtomName= !mShowAtomName;
   this->CrystUpdate();
}

void WXGLCrystalCanvas::OnShowCursor( wxCommandEvent & WXUNUSED(event))
{
   if(mShowCursor) mpPopUpMenu->SetLabel(ID_GLCRYSTAL_MENU_SHOWCURSOR, "Show Cursor");
   else mpPopUpMenu->SetLabel(ID_GLCRYSTAL_MENU_SHOWCURSOR, "Hide Cursor");
   mShowCursor= !mShowCursor;
   this->CrystUpdate();
}

void WXGLCrystalCanvas::OnSetCursor( wxCommandEvent & WXUNUSED(event))
{
  VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnSetCursor",1)
  REAL x=mX0;
  REAL y=mY0;
  REAL z=mZ0;
  mpWXCrystal->GetCrystal().OrthonormalToFractionalCoords(x,y,z);

  mViewCntr.x = (mcellbbox.xMax+mcellbbox.xMin)/2. - x;
  mViewCntr.y = (mcellbbox.yMax+mcellbbox.yMin)/2. - y;
  mViewCntr.z = (mcellbbox.zMax+mcellbbox.zMin)/2. - z;
  UserXYZBox *BoxDlg = new UserXYZBox(this,
       "Set fractional coordinates for view\ncenter and cursor position",
			 mViewCntr);
  if (BoxDlg->ShowModal() == wxID_OK ) {
     mViewCntr =  BoxDlg->GetXYZ();
     VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnSetCursor (frac) = " <<
		 mViewCntr.x << "," << mViewCntr.y << "," << mViewCntr.z,1)
     mX0 = (mcellbbox.xMax+mcellbbox.xMin)/2. - mViewCntr.x;
     mY0 = (mcellbbox.yMax+mcellbbox.yMin)/2. - mViewCntr.y;
     mZ0 = (mcellbbox.zMax+mcellbbox.zMin)/2. - mViewCntr.z;
     mpWXCrystal->GetCrystal().FractionalToOrthonormalCoords(mX0, mY0, mZ0);
     VFN_DEBUG_MESSAGE("...ortho" << mX0 << "," << mY0 << "," << mZ0,1)
     Refresh(FALSE);
  }
}

void WXGLCrystalCanvas::OnLoadFourierGRD( wxCommandEvent & WXUNUSED(event))
{
   wxFileDialog fd((wxWindow*)this, "Choose a file containing a Fourier Map",
           "", "", "Fourier Map files (*.grd)|*.grd", wxOPEN | wxFILE_MUST_EXIST);
   //if okay then read Fourier map, run MC on it and display the triangles
   if(fd.ShowModal() == wxID_OK)
   {
      const string filename=fd.GetPath().c_str();
      UnitCellMapImport *pMap=new UnitCellMapImport(mpWXCrystal->GetCrystal());
      if (pMap->ImportGRD(filename) == 0)
      {
         string tmp="Error reading Fourier file:"+filename;
         wxMessageBox(tmp.c_str(), "File error", wxOK, this);
      return;
      }
      this->AddFourier(pMap);
   }
}

void WXGLCrystalCanvas::OnLoadFourierDSN6( wxCommandEvent & WXUNUSED(event))
{
   wxFileDialog fd((wxWindow*)this, "Choose a file containing a Fourier Map",
           "", "", "Fourier Map files (*.DN6)|*.DN6", wxOPEN | wxFILE_MUST_EXIST);
   //if okay then read Fourier map, run MC on it and display the triangles
   if(fd.ShowModal() == wxID_OK)
   {
      const string filename=fd.GetPath().c_str();
      UnitCellMapImport *pMap=new UnitCellMapImport(mpWXCrystal->GetCrystal());
      if (pMap->ImportDSN6(filename) == 0)
      {
         string tmp="Error reading Fourier file:"+filename;
         wxMessageBox(tmp.c_str(), "File error", wxOK, this);
      return;
      }
      this->AddFourier(pMap);
   }
}

void WXGLCrystalCanvas::AddFourier(UnitCellMapImport *map)
{
   mvpUnitCellMapImport.push_back(map);
   wxBusyInfo wait("Processing Fourier Map...");
   {
      float contour=map->Mean()+2*map->StandardDeviation();
      if(contour>map->Max()) contour=map->Mean()+0.75*(map->Max()-map->Mean());
      mvpUnitCellMapGLList.push_back(make_pair(make_pair(mvpUnitCellMapImport.back(),contour),
                                               new UnitCellMapGLList) );
      switch(mvpUnitCellMapGLList.size())
      {
         case 1: mvpUnitCellMapGLList.back().second->SetColour(1.,0.,0.,1.);break;
         case 2: mvpUnitCellMapGLList.back().second->SetColour(0.,0.,1.,1.);break;
         default:mvpUnitCellMapGLList.back().second->SetColour(0.,1.,0.,1.);break;
      }
      this->SetCurrent();
      mvpUnitCellMapGLList.back().second->GenList(*map,
		     this, mvpUnitCellMapGLList.back().first.second);
      mvpUnitCellMapGLList.back().second->SetName(map->GetName());
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

void WXGLCrystalCanvas::OnChangeContour( wxCommandEvent & WXUNUSED(event))
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
   mvpUnitCellMapGLList[mapgl].second->GenList(*(mvpUnitCellMapGLList[mapgl].first.first),
					       this, mvpUnitCellMapGLList[mapgl].first.second);
   this->CrystUpdate();
}

void WXGLCrystalCanvas::OnAddContour( wxCommandEvent & WXUNUSED(event))
{
   int map=this->UserSelectUnitCellMapImport();
   
   // Choose contour value
      const REAL min  =mvpUnitCellMapImport[map]->Min();
      const REAL max  =mvpUnitCellMapImport[map]->Max();
      const REAL mean =mvpUnitCellMapImport[map]->Mean();
      const REAL sigma=mvpUnitCellMapImport[map]->StandardDeviation();
      double contour=mean+2*sigma;
      if(contour>max) contour=mean+0.75*(max-mean);
      wxString strValue;
      strValue.Printf("%lf",contour);
      wxString info;
      info.Printf(_T("Add a new contour value\n")
                  _T("For this map: min  =%6.3f\n")
                  _T("              max  =%6.3f\n")
                  _T("              mean =%6.3f\n")
                  _T("              sigma=%6.3f\n\n")
                  _T("Recommended values are mean + 2sigma\n")
                  _T("i.e. %6.3f for a Fobs or Fcalc map\n")
                  _T("or %6.3f and %6.3f for a Fobs-Fcalc map\n"),
                  min,max,mean,sigma,
                  mean+2*sigma,mean+2*sigma,mean-2*sigma
                 );
                  
      wxTextEntryDialog contourValueDialog(this,info,
                              "Add contour value",strValue,wxOK | wxCANCEL);
      if(wxID_OK!=contourValueDialog.ShowModal())  return;
      contourValueDialog.GetValue().ToDouble(&contour);
   // Choose colour
      wxColor ncolor(255,0,0);
      ncolor = wxGetColourFromUser((wxWindow*)this, ncolor);   
      if(!(ncolor.Ok())) return;
   // Add display map
      wxBusyInfo wait("Processing Fourier Map...");
      mvpUnitCellMapGLList.push_back(make_pair(make_pair(mvpUnitCellMapImport[map],
                                                         (float)contour),
                                               new UnitCellMapGLList()) );
      mvpUnitCellMapGLList.back().second->SetName(mvpUnitCellMapImport[map]->GetName());
      mvpUnitCellMapGLList.back().second
         ->SetColour(ncolor.Red()/255.0,ncolor.Green()/255.0,ncolor.Blue()/255.0,0.5);
      this->SetCurrent();
      mvpUnitCellMapGLList.back().second->GenList(*mvpUnitCellMapImport[map],
						  this, mvpUnitCellMapGLList.back().first.second);
      this->CrystUpdate();
}

void WXGLCrystalCanvas::OnShowFourier( wxCommandEvent & WXUNUSED(event))
{
   if(mShowFourier == TRUE) mpPopUpMenu->SetLabel(ID_GLCRYSTAL_MENU_SHOWFOURIER, "Show Fourier Map");
   else mpPopUpMenu->SetLabel(ID_GLCRYSTAL_MENU_SHOWFOURIER, "Hide Fourier Map");
   mShowFourier = !mShowFourier;
   this->CrystUpdate();
}

void WXGLCrystalCanvas::OnFourierChangeColor( wxCommandEvent & WXUNUSED(event))
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

void WXGLCrystalCanvas::OnFourierChangeBbox( wxCommandEvent & WXUNUSED(event))
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
         pos->second->GenList(*(pos->first.first),this, pos->first.second);

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

void WXGLCrystalCanvas::OnUnloadFourier( wxCommandEvent & WXUNUSED(event))
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

void WXGLCrystalCanvas::OnShowWire( wxCommandEvent & WXUNUSED(event))
{
   vector<pair<pair<const UnitCellMapImport*,float>,UnitCellMapGLList* > >::iterator pos;
   for(pos=mvpUnitCellMapGLList.begin();pos != mvpUnitCellMapGLList.end();pos++)
      pos->second->ToggleShowWire();

   this->CrystUpdate();
}

void WXGLCrystalCanvas::OnPOVRay( wxCommandEvent & WXUNUSED(event))
{
   WXCrystValidateAllUserInput();
   wxFileDialog save(this,"Choose filename","","","*.pov",wxSAVE | wxOVERWRITE_PROMPT);
   if(save.ShowModal() != wxID_OK) return;
   
   ofstream os(save.GetPath().c_str());
   //ofstream os("test.pov");
   
   os << "// This File was created by FOX/ObjCryst++ (http://objcryst.sf.net)"<<endl
      << "//"<<endl
      << "// You can produce a ray-traced image using POV-Ray, freely available"<<endl
      << "//from http://www.povray.org"<<endl
      << "//"<<endl
      << "// Example command line to produce an anti-aliase 640x480 image: "<<endl
      << "//      povray +Ifile.pov +.pov +W640 +H480 +A +Q11"<<endl
      << "//   You can add '+UA' at the end to have a transparent background"<<endl
      << "//   You can add '+kff10' to generate 10 rotated images for an animation"<<endl
      << "//   (see the 'clock' in the 'OrientPitch' definition below)"<<endl
      << "//"<<endl
      << "// Notes:"<<endl
      << "//  - This POVRay file is written to produce a 4/3 image, e.g. 640x480"<<endl
      << "//    If your image in the FOX 3D view was not 4/3, some parts may be cut"<<endl
      << "//    You can then get the full structure by increasing the 'angle' "<<endl
      << "//    (viewing angle) in the camera settings below."<<endl
      << "//    You can change the aspect ratio (e.g. to produce a square image)"<<endl
      << "//    by changing the 'right  <-1.33,0,0>' statement in the camera definition"<<endl
      << "//  - You can change the orientation of the view by changing the"<<endl
      << "//    OrientRoll, OrientPitch and OrientYaw angles just below (in degrees)"<<endl
      << "//  - You can change the aspects of atoms by altering the macros below."<<endl
      << "//    The radius of atoms is by default 1/3 of their tabulated atomic radius,"<<endl
      << "//    i.e. as in the FOX/ObjCryst++ 3D Crystal view. To modify this you can"<<endl
      << "//    change the second line of the 'ObjCrystAtom' macro to (e.g. for full radius):"<<endl
      << "//    '{ <atomx,atomy,atomz>,atomr*1.0'"<<endl
      << "//  - The colour of atoms, bonds (free and non-free torsions) can be changed"<<endl
      << "//    in the 'GLOBAL DECLARATIONS FOR ATOMS & BONDS' section"<<endl
      << "//  - Just for fun, you can try getting *very close* to one of the atoms,"<<endl
      << "//    and, in the 'ObjCrystAtom' macro at the end of the 'finish'"<<endl
      << "//    statement, change the 'reflection' value to 1.0, "<<endl
      << "//    to get a mirror effect on the atoms..."<<endl
      << "//"<<endl
      << "// See http://povray.org/documentation/ for more options"<<endl
      << "//"<<endl<<endl;
   
   os << "// Description of Crystal :" << mpWXCrystal->GetCrystal().GetName() <<endl;
   os << "global_settings { assumed_gamma 2.2 ambient_light rgb <1,1,1>}"<<endl;
   float m[4][4];
   REAL xcam=0,ycam=0,zcam=mDist;
   build_rotmatrix( m,mQuat);
   
   REAL x=(mcellbbox.xMin+mcellbbox.xMax)/2.;
   REAL y=(mcellbbox.yMin+mcellbbox.yMax)/2.;
   REAL z=(mcellbbox.zMin+mcellbbox.zMax)/2.;
   mpWXCrystal->GetCrystal().FractionalToOrthonormalCoords(x,y,z);
   x-=mX0;
   y-=mY0;
   z-=mZ0;
   {
      const REAL q1=mQuat[0];const REAL q2=mQuat[1];
      const REAL q3=mQuat[2];const REAL q4=mQuat[3];
      
      REAL yaw =(q4*q4 + q1*q1 - q2*q2 - q3*q3);
      if(abs(yaw)>1e-6)  yaw  =atan( 2*(q1*q2+q4*q3) /yaw )*RAD2DEG;
      else { if((q1*q2+q4*q3)>0) yaw =90.; else yaw =-90;}
      
      const REAL pitch=asin(-2*(q1*q3-q4*q2))*RAD2DEG;
      
      REAL roll=(q4*q4 - q1*q1 - q2*q2 + q3*q3);
      if(abs(roll)>1e-6) roll =atan( 2*(q4*q1+q2*q3) /roll)*RAD2DEG;
      else { if((q4*q1+q2*q3)>0) roll=90.; else roll=-90;}
      
      if((q4*q4 + q1*q1 - q2*q2 - q3*q3)<0) yaw  +=180;
      if((q4*q4 - q1*q1 - q2*q2 + q3*q3)<0) roll +=180;
      
      os<<endl;
      os << "#declare OrientRoll="<<roll<<";"<<endl;
      os << "#declare OrientPitch="<<pitch<<"+360*clock;"<<endl;
      os << "#declare OrientYaw="<<yaw<<";"<<endl<<endl;
   }
   
   os << "camera" <<endl;
   os << "{"<<endl;
   os << "    location  <"<<xcam+x<<","<<ycam+y<<","<<zcam+z<<">"<<endl
      << "    look_at   <" << x << "," << y << "," << z <<">"<<endl
      << "    angle   "<< mViewAngle*1.2 <<endl
      << "    right  <-1.33,0,0> //change handedness as in OpenGL, aspect ratio=4/3"<<endl
      << "    translate   <" <<-x << "," <<-y << "," <<-z <<">"<<endl
      << "    rotate  <OrientRoll,0,0>" <<endl
      << "    rotate  <0,OrientPitch,0>" <<endl
      << "    rotate  <0,0,OrientYaw>" <<endl
      << "    translate   <" << x << "," << y << "," << z <<">"<<endl
      << "}"<<endl;

   REAL xlight=-1000,ylight=1000,zlight=1000;
   os << "light_source"<<endl;
   os << "{" <<endl
      << "   <"<<xlight<<","<<ylight<<","<<zlight<<">"<<endl
      << "   colour rgb <1.0,1.0,1.0>" <<endl
      << "   //shadowless" <<endl
      << "    translate   <" <<-x << "," <<-y << "," <<-z <<">"<<endl
      << "    rotate  <OrientRoll,0,0>" <<endl
      << "    rotate  <0,OrientPitch,0>" <<endl
      << "    rotate  <0,0,OrientYaw>" <<endl
      << "    translate   <" << x << "," << y << "," << z <<">"<<endl
      << "}" <<endl<<endl;
   
   os << "background { colour rgb <0.0, 0.0, 0.0> }"<<endl<<endl;
   
   CrystalPOVRayOptions options;
   options.mXmin=mcellbbox.xMin;
   options.mXmax=mcellbbox.xMax;
   options.mYmin=mcellbbox.yMin;
   options.mYmax=mcellbbox.yMax;
   options.mZmin=mcellbbox.zMin;
   options.mZmax=mcellbbox.zMax;
   options.mShowLabel=mShowAtomName;
   if(mShowCrystal)
   {
      mpWXCrystal->GetCrystal().POVRayDescription(os,options);
   }
   if(mShowFourier)
   {
      wxBusyInfo wait("Processing Fourier Map...");
      // use cell bbox if mapbbox has zero volume (default)
      if (mmapbbox.xMin != mmapbbox.xMax)
      {
         options.mXmin=mmapbbox.xMin;
         options.mXmax=mmapbbox.xMax;
         options.mYmin=mmapbbox.yMin;
         options.mYmax=mmapbbox.yMax;
         options.mZmin=mmapbbox.zMin;
         options.mZmax=mmapbbox.zMax;
      }

      os<<"/////////////////// FOURIER MAPS///////////////////////////"<<endl;
      vector<pair<pair<const UnitCellMapImport*,float>,UnitCellMapGLList* > >::
         const_iterator pos;
      for(pos=mvpUnitCellMapGLList.begin();pos != mvpUnitCellMapGLList.end();++pos)
      {
         const float *prgbf=pos->second->GetColour();
         if(pos->second->ShowWire())
         {
            os << "#macro ObjCrystMeshTriangle(x1,y1,z1,x2,y2,z2,x3,y3,z3,"
               << "nx1,ny1,nz1,nx2,ny2,nz2,nx3,ny3,nz3)"<<endl
               << "   cylinder"<<endl
               << "   {  <x1,y1,z1>,"<<endl
               << "      <x2,y2,z2>,"<<endl
               << "      0.01"<<endl
               << "      finish {ambient 0.5 diffuse 0.4}"<<endl
               << "      pigment { colour rgb<"
               <<*(prgbf+0)<<","<<*(prgbf+1)<<","<<*(prgbf+2)<<">}"<<endl
               << "      no_shadow"<<endl
               << "   }"<<endl
               << "   cylinder"<<endl
               << "   {  <x2,y2,z2>,"<<endl
               << "      <x3,y3,z3>,"<<endl
               << "      0.01"<<endl
               << "      finish {ambient 0.5 diffuse 0.4}"<<endl
               << "      pigment { colour rgb<"
               <<*(prgbf+0)<<","<<*(prgbf+1)<<","<<*(prgbf+2)<<">}"<<endl
               << "      no_shadow"<<endl
               << "   }"<<endl
               << "   cylinder"<<endl
               << "   {  <x1,y1,z1>,"<<endl
               << "      <x3,y3,z3>,"<<endl
               << "      0.01"<<endl
               << "      finish {ambient 0.5 diffuse 0.4}"<<endl
               << "      pigment { colour rgb<"
               <<*(prgbf+0)<<","<<*(prgbf+1)<<","<<*(prgbf+2)<<">}"<<endl
               << "      no_shadow"<<endl
               << "   }"<<endl
               << "#end"<<endl<<endl;
            pos->first.first->POVRayDescription(os,pos->first.second,options);
         }
         else
         {
            os << "#macro ObjCrystMeshTriangle(x1,y1,z1,x2,y2,z2,x3,y3,z3,"
               << "nx1,ny1,nz1,nx2,ny2,nz2,nx3,ny3,nz3)"<<endl
               << "      smooth_triangle"<<endl
               <<"       {<x1,y1,z1>,<nx1,ny1,nz1>"<<endl
               <<"        <x2,y2,z2>,<nx2,ny2,nz2>"<<endl
               <<"        <x3,y3,z3>,<nx3,ny3,nz3>}"<<endl
               << "#end"<<endl<<endl;
            os << "   mesh"<<endl
               << "   {"<<endl;
            pos->first.first->POVRayDescription(os,pos->first.second,options);
            os << "      texture"<<endl
               << "      {"<<endl
               << "         finish {ambient 0.5 diffuse 0.4}"<<endl
               << "         pigment { colour rgb<"
               <<*prgbf++<<",";
            os <<*prgbf++<<",";
            os <<*prgbf++<<">}"<<endl
               << "      }"<<endl
               << "      no_shadow"<<endl
               << "   }"<<endl;
         }
      }
      
   }
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

void UserSelectBoundingBox::OnOk (wxCommandEvent & WXUNUSED(event)) {
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

void UserXYZBox::OnOk (wxCommandEvent & WXUNUSED(event)) {
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
