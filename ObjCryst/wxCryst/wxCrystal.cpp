/*  ObjCryst++ Object-Oriented Crystallographic Library
    (c) 2000-2009 Vincent Favre-Nicolin vincefn@users.sourceforge.net
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
#include "wx/notebook.h"
#include "wx/minifram.h"

#include "ObjCryst/wxCryst/wxCrystal.h"

#include "wx/colordlg.h"
#include "wx/progdlg.h"
#include "wx/busyinfo.h"
#include "wx/config.h"

#include "ObjCryst/Quirks/Chronometer.h"
#include "ObjCryst/ObjCryst/Atom.h"
#include "ObjCryst/ObjCryst/ZScatterer.h"
#include "ObjCryst/ObjCryst/Molecule.h"
#include "ObjCryst/ObjCryst/PowderPattern.h"
#include "ObjCryst/ObjCryst/ScatteringPowerSphere.h"
#include "ObjCryst/ObjCryst/Polyhedron.h"

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
     #include "gl/glu.h"
   #endif
   
   #ifdef HAVE_GLUT
      #ifdef __WXMAC__ // For the wxMac version of wxWindows, i.e. with the "Aqua" look
         #include <GLUT/glut.h>
      #else
         #include "GL/glut.h"
       #endif
   #endif
   #ifdef HAVE_FFTW
      #include "fftw3.h"
   #endif
#endif

extern "C" {
#include "ObjCryst/wxCryst/trackball.h"
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

//#include "ObjCryst/ObjCryst/Map.cpp"

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
   if(mpWXCrystal!=0) mpWXCrystal->NotifyDeleteListWin(this);
}

////////////////////////////////////////////////////////////////////////
//
//    WXCrystal Grid objects
//
////////////////////////////////////////////////////////////////////////
WXCrystal::RowScattPow::RowScattPow():
mName("H"),
mBiso(1.0),mFormalCharge(0.0),mR(1.0),mG(1.0),mB(1.0),mMaximumLikelihoodError(0.0),mNbGhostAtoms(0.0),
mNeedUpdateUI(true),mIdx(-1)
{}

////////////////////////////////////////////////////////////////////////
//
//    Convert a list of atoms to one molecule
//
////////////////////////////////////////////////////////////////////////
Molecule* Atoms2Molecule(list<Atom *> &vAtom)
{
   VFN_DEBUG_ENTRY("Atoms2Molecule()",6)
   Molecule *mol=new Molecule((*vAtom.begin())->GetCrystal(),"Molecule");
   const unsigned long nb=vAtom.size();
   REAL x0=0,y0=0,z0=0;
   unsigned int i=0;
   for(list<Atom *>::iterator pos=vAtom.begin();pos!=vAtom.end();++pos)
   {
      REAL x=(*pos)->GetX();
      REAL y=(*pos)->GetY();
      REAL z=(*pos)->GetZ();
      (*pos)->GetCrystal().FractionalToOrthonormalCoords(x,y,z);
      x0+=x;
      y0+=y;
      z0+=z;
      mol->AddAtom(x,y,z,&((*pos)->GetScatteringPower()),(*pos)->GetName());
      mol->GetAtom(i++).SetOccupancy((*pos)->GetOccupancy());
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
//    WXCrystal
//
////////////////////////////////////////////////////////////////////////
static const long ID_CRYSTAL_MENU_SAVECIF                       =WXCRYST_ID();
static const long ID_CRYSTAL_MENU_SAVETEXT                      =WXCRYST_ID();
static const long ID_CRYSTAL_MENU_DISPLAY                       =WXCRYST_ID();
static const long ID_CRYSTAL_MENU_DISPLAY_3DVIEW                =WXCRYST_ID();
static const long ID_CRYSTAL_MENU_SCATT                         =WXCRYST_ID();
static const long ID_CRYSTAL_MENU_PAR_SETRELATIVEXYZLIMITS      =WXCRYST_ID();
static const long ID_CRYSTAL_MENU_PAR_TEST_RANDOM_MOVES         =WXCRYST_ID();
static const long ID_CRYSTAL_MENU_SCATT_REMOVESCATTPOW          =WXCRYST_ID();
static const long ID_CRYSTAL_MENU_SCATT_ADDSCATTPOWATOM         =WXCRYST_ID();
static const long ID_CRYSTAL_MENU_SCATT_ADDSCATTPOWSPHERE       =WXCRYST_ID();
static const long ID_CRYSTAL_MENU_SCATT_ADDATOM                 =WXCRYST_ID();
static const long ID_CRYSTAL_MENU_SCATT_IMPORTATOMLIST          =WXCRYST_ID();
static const long ID_CRYSTAL_MENU_SCATT_ADDZSCATTERER           =WXCRYST_ID();
static const long ID_CRYSTAL_MENU_SCATT_ADDMOLECULE             =WXCRYST_ID();
static const long ID_CRYSTAL_MENU_SCATT_ATOMS2MOLECULE          =WXCRYST_ID();
static const long ID_CRYSTAL_MENU_SCATT_IMPORTFENSKEHALLZMATRIX =WXCRYST_ID();
static const long ID_CRYSTAL_MENU_SCATT_IMPORTNAMEDZMATRIX      =WXCRYST_ID();
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
static const long ID_GLCRYSTAL_WINDOW                           =WXCRYST_ID();
static const long ID_CRYSTAL_WIN_SCATTPOW                       =WXCRYST_ID();
static const long ID_CRYSTAL_WIN_ANTIBUMP                       =WXCRYST_ID();
static const long ID_CRYSTAL_WIN_BONDVALENCE                    =WXCRYST_ID();
static const long ID_CRYSTAL_MENU_SHOW_SCATTPOW_WIN             =WXCRYST_ID();
//static const long ID_CRYSTAL_MENU_SHOW_PDF                      =WXCRYST_ID();

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
   EVT_MENU(ID_CRYSTAL_MENU_PAR_TEST_RANDOM_MOVES,    WXCrystal::OnMenuTestRandomMoves)
#ifdef OBJCRYST_GL
   EVT_MENU(ID_CRYSTAL_MENU_DISPLAY_3DVIEW,           WXCrystal::OnMenuCrystalGL)
#endif
   EVT_MENU(ID_CRYSTAL_MENU_SCATT_ADDSCATTPOWATOM,    WXCrystal::OnMenuAddScattPowAtom)
   EVT_MENU(ID_CRYSTAL_MENU_SCATT_ADDSCATTPOWSPHERE,  WXCrystal::OnMenuAddScattPowSphere)
   EVT_MENU(ID_CRYSTAL_MENU_SCATT_REMOVESCATTPOW,     WXCrystal::OnMenuRemoveScattPow)
   EVT_MENU(ID_CRYSTAL_MENU_SCATT_ADDATOM,            WXCrystal::OnMenuAddScatterer)
   EVT_MENU(ID_CRYSTAL_MENU_SCATT_IMPORTATOMLIST,     WXCrystal::OnMenuAddScatterer)
   EVT_MENU(ID_CRYSTAL_MENU_SCATT_ADDZSCATTERER,      WXCrystal::OnMenuAddScatterer)
   EVT_MENU(ID_CRYSTAL_MENU_SCATT_ADDMOLECULE,        WXCrystal::OnMenuAddScatterer)
   EVT_MENU(ID_CRYSTAL_MENU_SCATT_ATOMS2MOLECULE,     WXCrystal::OnMenuAtoms2Molecule)
   EVT_MENU(ID_CRYSTAL_MENU_SCATT_IMPORTFENSKEHALLZMATRIX,WXCrystal::OnMenuImportMoleculeFromFenskeHallZMatrix)
   EVT_MENU(ID_CRYSTAL_MENU_SCATT_IMPORTNAMEDZMATRIX, WXCrystal::OnMenuImportMoleculeFromFenskeHallZMatrix)
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
//   EVT_MENU(ID_CRYSTAL_MENU_SHOW_PDF,                 WXCrystal::OnMenuPDF)
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
//,mpPDF(0)
{
   VFN_DEBUG_MESSAGE("WXCrystal::WXCrystal()",6)
   //this->SetBackgroundColour("Red");
   //mpWXTitle->SetBackgroundColour(wxColour(255,200,200));
   mpWXTitle->SetForegroundColour(wxColour(255,0,0));
   mpWXTitle->SetSize(400,-1);
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
         mpMenuBar->GetMenu(ID_REFOBJ_MENU_PAR).AppendSeparator();
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_PAR,ID_CRYSTAL_MENU_PAR_TEST_RANDOM_MOVES,
                                "Test Random Moves for 30s");
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
         mpMenuBar->AddMenuItem(ID_CRYSTAL_MENU_SCATT,ID_CRYSTAL_MENU_SCATT_IMPORTATOMLIST,
                                "Import a List of Atoms");
         mpMenuBar->AddMenuItem(ID_CRYSTAL_MENU_SCATT,ID_CRYSTAL_MENU_SCATT_ADDMOLECULE,
                                "Add Molecule");
         mpMenuBar->AddMenuItem(ID_CRYSTAL_MENU_SCATT,ID_CRYSTAL_MENU_SCATT_ATOMS2MOLECULE,
                                "Convert Atoms to a Molecule");
         mpMenuBar->AddMenuItem(ID_CRYSTAL_MENU_SCATT,ID_CRYSTAL_MENU_SCATT_IMPORTFENSKEHALLZMATRIX,
                                "Import Molecule from Fenske-Hall Z-Matrix");
         mpMenuBar->AddMenuItem(ID_CRYSTAL_MENU_SCATT,ID_CRYSTAL_MENU_SCATT_IMPORTNAMEDZMATRIX,
                                "Import Molecule from a named Z-Matrix");
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
         //mpMenuBar->AddMenuItem(ID_CRYSTAL_MENU_DISPLAY,ID_CRYSTAL_MENU_SHOW_PDF,
         //                       "PDF");

      mpSizer->SetItemMinSize(mpMenuBar,
                              mpMenuBar->GetSize().GetWidth(),
                              mpMenuBar->GetSize().GetHeight());

      // KLUDGE : this only works as long as the option order does not change !
      dynamic_cast<WXFieldOption *>(mpCrystal->GetOption(0).WXGet())->SetToolTip(
         _T("Use this option ONLY if you want to use\n")
         _T("a higher symmetry than the one allowed by\n")
         _T("the spacegroup. This can be useful to search\n")
         _T("structures derived from higher symmetries.\n\n")
         _T("This option should almost never be used."));
      
      dynamic_cast<WXFieldOption *>(mpCrystal->GetOption(1).WXGet())->SetToolTip(
         _T("This option allows Fox to automatically adjust\n")
         _T("the occupancy of atoms that are on a special position,\n")
         _T("or overlapping with another (e.g. two oxygens from\n")
         _T("two polyhedra).\n\n")
         _T("Practically you should choose:\n")
         _T("- Yes for inorganic structures\n")
         _T("- No for organic structures\n\n")
         _T("This option increases computing time\n")
         _T("by up to 50%, so only use when necessary\n\n")
         _T("In doubt, choose Yes"));
      
      dynamic_cast<WXFieldOption *>(mpCrystal->GetOption(2).WXGet())->SetToolTip(
         _T("This option only affects the 3D display,\n")
         _T("and is used to display the enantiomer\n")
         _T("of the crystal structure.\n\n")
         _T("This can be used to compare several\n")
         _T("crystal structures."));


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
      pWXFieldBumpMerge->SetFormat(_T("%8.2f"));
      pAntiBumpScale->SetFormat(_T("%8.2f"));
      pWXFieldBumpMerge->SetToolTip(_T("Current anti-bump cost"));
      pAntiBumpScale->SetToolTip(
         _T("Scale (multiplier) for the anti-bump cost.\n")
         _T("If 0, the anti-bump will be ignored and not calculated\n")
         _T("during optimization (saving time)\n\n")
         _T("Use a value larger than 1 to increase the importance\n")
         _T("of the anti-bump relatively to the diffraction data Chi^2\n\n")
         _T("Note that anti-bump should only be used if the diffraction data\n\n")
         _T("is not of good enough quality to ensure finding the correct\n\n")
         _T("structure."));
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
      pWXFieldBondValence->SetFormat(_T("%8.2f"));
      pBondValenceScale->SetFormat(_T("%8.2f"));
      pWXFieldBondValence->SetToolTip(_T("Current bond valence cost"));
      pBondValenceScale->SetToolTip(
         _T("Scale (multiplier) for the bond valence cost.\n")
         _T("If 0, the bond valence will be ignored and not calculated\n")
         _T("during optimization (saving time)\n\n")
         _T("Use a value larger than 1 to increase the importance\n")
         _T("of the bond valence relatively to the diffraction data Chi^2\n\n")
         _T("Note that bond valence should only be used if the diffraction data\n\n")
         _T("is not of good enough quality to ensure finding the correct\n\n")
         _T("structure."));
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
      
      dynamic_cast<WXFieldRefPar *>(pFieldLatticeA)->SetFormat(_T("%8.4f"));
      dynamic_cast<WXFieldRefPar *>(pFieldLatticeB)->SetFormat(_T("%8.4f"));
      dynamic_cast<WXFieldRefPar *>(pFieldLatticeC)->SetFormat(_T("%8.4f"));
      dynamic_cast<WXFieldRefPar *>(pFieldLatticeAlpha)->SetFormat(_T("%8.3f"));
      dynamic_cast<WXFieldRefPar *>(pFieldLatticeBeta)->SetFormat(_T("%8.3f"));
      dynamic_cast<WXFieldRefPar *>(pFieldLatticeGamma)->SetFormat(_T("%8.3f"));
      
      pFieldLatticeA->SetToolTip(_T("Lattice length parameter (in Angstroems)"));
      pFieldLatticeB->SetToolTip(_T("Lattice length parameter (in Angstroems)"));
      pFieldLatticeC->SetToolTip(_T("Lattice length parameter (in Angstroems)"));
      pFieldLatticeAlpha->SetToolTip(_T("Lattice angle parameter (in degrees)"));
      pFieldLatticeBeta->SetToolTip(_T("Lattice angle parameter (in degrees)"));
      pFieldLatticeGamma->SetToolTip(_T("Lattice angle parameter (in degrees)"));
      
   // SpaceGroup
      mpFieldSpacegroup=new WXFieldName(this,"SpaceGroup:",this,ID_CRYSTAL_SPACEGROUP,100);
      mpSizer->Add(mpFieldSpacegroup,0,wxALIGN_LEFT);
      mList.Add(mpFieldSpacegroup);
      
      mpFieldSpacegroup->SetToolTip(_T("Spacegroup Symbol. You can use:\n\n")
                                    _T("- spacegroup number: \"1\" \"62\" ... \"227\",\"230\"\n")
                                    _T("- Hermann-Mauguin symbol: \"P1\" \"Pnma\" ... \"Fd3m\",\"Ia3d\"\n")
                                    _T("- Hall symbol: \"P1\" \"-P 2ac 2n\" ... \"-F 4vw 2vw 3\",\"-I 4bd 2c 3\"\n\n")
                                    _T("ORIGIN CHOICE: for some spacegroups there are several\n")
                                    _T("possible origins - the default is the one on \n")
                                    _T("the center of symmetry (origin 2). You can specify\n")
                                    _T(" the origin by writing \"Fd3m:1\" or \"Fd3m:2\"\n\n")
                                    _T("CELL CHOICE: to specify a rhomboedral or hexagonal unit cell,\n")
                                    _T("append R or H to the symbol:\"R-3:R\"or \"R-3:H\"\n"));

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
   {
      bool val;
      if(!wxConfigBase::Get()->HasEntry(_T("Crystal/BOOL/Automatically open crystal 3D view")))
         wxConfigBase::Get()->Write(_T("Crystal/BOOL/Automatically open crystal 3D view"), false);
      else
      {
         wxConfigBase::Get()->Read(_T("Crystal/BOOL/Automatically open crystal 3D view"), &val);
         if(val)
         {
            wxCommandEvent event(wxEVT_COMMAND_MENU_SELECTED,ID_CRYSTAL_MENU_DISPLAY_3DVIEW);
            wxPostEvent(this,event);
         }
      }
   }
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
   
   wxWakeUpIdle();
   
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
   // Necessary to change the "used" status of unit cell parameters.
   if((false==this->GetCrystal().IsBeingRefined()) && wxThread::IsMain() ) this->GetCrystal().InitRefParList();
   
   if((false==this->GetCrystal().IsBeingRefined()) && wxThread::IsMain() &&(mpScattPowWin!=0)&&(mpAntiBumpWin!=0)&&(mpBondValenceWin!=0))
   {
      //set<ScatteringPowerAtom*> vpRemovedScattPow;
      //set<ScatteringPowerAtom*> vpAddedScattPow;
      
      bool needLayout=false;
      // Delete rows & cols as required
      for(map<ScatteringPowerAtom*,RowScattPow>::iterator pos=mvpRowScattPow.begin();pos!=mvpRowScattPow.end();)
         if(this->GetCrystal().GetScatteringPowerRegistry().Find(pos->second.mName,"ScatteringPowerAtom",true)<0)
         {
            VFN_DEBUG_MESSAGE("WXCrystal::CrystUpdate(): Removing scattering power: "<<pos->second.mName,5)
            mpScattPowWin->DeleteRows(mvpRowScattPow.size()-1,1,false);
            mpAntiBumpWin->DeleteRows(mvpRowScattPow.size()-1,1,false);
            mpBondValenceWin->DeleteRows(mvpRowScattPow.size()-1,1,false);
            mpAntiBumpWin->DeleteCols(mvpRowScattPow.size()-1,1,false);
            mpBondValenceWin->DeleteCols(mvpRowScattPow.size()-1,1,false);
            mvpRowScattPow.erase(pos++);
            needLayout=true;
         }
         else ++pos; // See Josuttis, p.205
      // Add rows & cols as required
      for(int i=0;i<this->GetCrystal().GetScatteringPowerRegistry().GetNb();++i)
      {
         ScatteringPower *s=&(this->GetCrystal().GetScatteringPowerRegistry().GetObj(i));
         if(s->GetClassName()=="ScatteringPowerAtom")
         {
            ScatteringPowerAtom *p=dynamic_cast<ScatteringPowerAtom *>(s);
            if(mvpRowScattPow.find(p)==mvpRowScattPow.end())
            {
               VFN_DEBUG_MESSAGE("WXCrystal::CrystUpdate(): Adding scattering power: "<<s->GetName(),5)
               mpScattPowWin->AppendRows();
               mpAntiBumpWin->AppendRows();
               mpAntiBumpWin->AppendCols();
               mpBondValenceWin->AppendRows();
               mpBondValenceWin->AppendCols();
               mvpRowScattPow.insert(make_pair(p,RowScattPow()));
               needLayout=true;
            }
         }
      }
      // Put the scattering powers in the same order as they have been declared,
      // for user convenience
      {
         int j=0; // :KLUDGE: number of scattering powers that are not ScatteringPowerAtom
         for(int i=0;i<this->GetCrystal().GetScatteringPowerRegistry().GetNb();++i)
         {
            ScatteringPower *s=&(this->GetCrystal().GetScatteringPowerRegistry().GetObj(i));
            if(s->GetClassName()=="ScatteringPowerAtom")
            {
               ScatteringPowerAtom *p=dynamic_cast<ScatteringPowerAtom *>(s);
               if(mvpRowScattPow[p].mIdx!=i-j)
               {
                  mvpRowScattPow[p].mIdx=i-j;
                  mvpRowScattPow[p].mNeedUpdateUI=true;
               }
            }
            else j++;
         }
      }
      if(needLayout)
      {
         mpScattPowWin->FitInside();
         mpAntiBumpWin->FitInside();
         mpBondValenceWin->FitInside();
      }
      // Update windows
      //if(mpScattPowWin!=0)
      {
         map<ScatteringPowerAtom*,RowScattPow>::iterator pos;
         for(pos=mvpRowScattPow.begin();pos!=mvpRowScattPow.end();++pos)
         {
            const string name=pos->first->GetName();
            const REAL biso=pos->first->GetBiso();
            const REAL formalCharge=pos->first->GetFormalCharge();
            const float *pRGB=pos->first->GetColourRGB();
            const REAL mlerror=pos->first->GetMaximumLikelihoodPositionError();
            const REAL nbghost=pos->first->GetMaximumLikelihoodNbGhostAtom();
            if(  (name   !=pos->second.mName)
               ||(biso   !=pos->second.mBiso)
               ||(formalCharge!=pos->second.mFormalCharge)
               ||(pRGB[0]!=pos->second.mR)
               ||(pRGB[1]!=pos->second.mG)
               ||(pRGB[2]!=pos->second.mB)
               ||(mlerror!=pos->second.mMaximumLikelihoodError)
               ||(nbghost!=pos->second.mNbGhostAtoms)
               || pos->second.mNeedUpdateUI)
            {
               pos->second.mName=name;
               pos->second.mBiso=biso;
               pos->second.mFormalCharge=formalCharge;
               pos->second.mR=pRGB[0];
               pos->second.mG=pRGB[1];
               pos->second.mB=pRGB[2];
               pos->second.mMaximumLikelihoodError=mlerror;
               pos->second.mNbGhostAtoms=nbghost;
               pos->second.mNeedUpdateUI=true;
            }
         }
      }
      //if(mpAntiBumpWin!=0)
      {
         map<ScatteringPowerAtom*,RowScattPow>::iterator pos,pos1;
         for(pos=mvpRowScattPow.begin();pos!=mvpRowScattPow.end();++pos)
         {
            const string name=pos->first->GetName();
            const Crystal::VBumpMergePar *pMap=&(mpCrystal->GetBumpMergeParList());
            vector<REAL> dist(mvpRowScattPow.size());
            for(pos1=mvpRowScattPow.begin();pos1!=mvpRowScattPow.end();++pos1)
            {
               Crystal::VBumpMergePar::const_iterator pos2;
               if(pos->first<pos1->first) pos2=pMap->find(make_pair(pos->first,pos1->first));
               else pos2=pMap->find(make_pair(pos1->first,pos->first));
               if(pos2==pMap->end()) dist[pos1->second.mIdx]=-999;
               else dist[pos1->second.mIdx]=sqrt(pos2->second.mDist2);
            }
            if(  (name!=pos->second.mName)
               ||(dist!=pos->second.mvAntiBumpDistance)
               || pos->second.mNeedUpdateUI)
            {
               pos->second.mName=name;
               pos->second.mvAntiBumpDistance=dist;
               pos->second.mNeedUpdateUI=true;
            }
         }
      }
      //if(mpBondValenceWin!=0)
      {
         map<ScatteringPowerAtom*,RowScattPow>::iterator pos,pos1;
         for(pos=mvpRowScattPow.begin();pos!=mvpRowScattPow.end();++pos)
         {
            const string name=pos->first->GetName();
            const std::map<pair<const ScatteringPower*,const ScatteringPower*>, REAL> *pMap=&(mpCrystal->GetBondValenceRoList());
            vector<REAL> ro(mvpRowScattPow.size());
            for(pos1=mvpRowScattPow.begin();pos1!=mvpRowScattPow.end();++pos1)
            {
               map<pair<const ScatteringPower*,const ScatteringPower*>, REAL>::const_iterator pos2;
               if(pos->first<pos1->first) pos2=pMap->find(make_pair(pos->first,pos1->first));
               else pos2=pMap->find(make_pair(pos1->first,pos->first));
               if(pos2==pMap->end()) ro[pos1->second.mIdx]=-999;
               else ro[pos1->second.mIdx]=pos2->second;
            }
            if(  (name!=pos->second.mName)
               ||(ro!=pos->second.mvBondValenceRo)
               || pos->second.mNeedUpdateUI)
            {
               pos->second.mName=name;
               pos->second.mvBondValenceRo=ro;
               pos->second.mNeedUpdateUI=true;
            }
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
         bool ok=mpConditionGLUpdate->IsOk();
         mMutexGLUpdate.Lock();
         wxCommandEvent event(wxEVT_COMMAND_MENU_SELECTED,ID_GLCRYSTAL_MENU_UPDATE);
         wxPostEvent(mpCrystalGL,event);
         wxWakeUpIdle();
         wxThread::This()->Yield();
         int ct=0;
         #ifdef __LINUX__
         while(mpConditionGLUpdate->WaitTimeout(200)!=wxCOND_NO_ERROR)
         {
            cout<<"WXCrystal::UpdateGL():timeout waiting for mpConditionGLUpdate release: #"<<++ct<<":"<<ok<<endl;
            wxWakeUpIdle();
            if(ct>10) break;//and hope for the best...
         }
         #else
         mpConditionGLUpdate->Wait();
         #endif
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
      mpCrystalGL->CrystUpdate();
      if(mpConditionGLUpdate!=0)
      {
         wxMutexLocker lock(mMutexGLUpdate);
         mpConditionGLUpdate->Signal();
      }
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
   wxFrame* frame;
   if(gvWindowPosition.count(ID_GLCRYSTAL_WINDOW))
     frame= new wxFrame(this,ID_GLCRYSTAL_WINDOW, wxString::FromAscii(mpCrystal->GetName().c_str()),
                        gvWindowPosition[ID_GLCRYSTAL_WINDOW].first,
                        gvWindowPosition[ID_GLCRYSTAL_WINDOW].second);
   else
     frame= new wxFrame(this,ID_GLCRYSTAL_WINDOW, wxString::FromAscii(mpCrystal->GetName().c_str()),
                        wxDefaultPosition,wxSize(400,400));

   mpCrystalGL=new WXGLCrystalCanvas(this,frame,-1);
   #if wxUSE_STATUSBAR
   frame->CreateStatusBar(1);
   frame->SetStatusText( wxString::FromAscii(mpCrystal->GetName().c_str()));
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
   wxFileDialog save(this,_T("Choose a file"),_T(""),_T(""),_T("*.cif"),wxFD_SAVE | wxFD_OVERWRITE_PROMPT);
   if(save.ShowModal() != wxID_OK) return;
   
   ofstream out(save.GetPath().ToAscii());
   if(!out) return;//:TODO:
   mpCrystal->CIFOutput(out);
   out.close();
}

void WXCrystal::OnMenuSaveText(wxCommandEvent & WXUNUSED(event))
{
   WXCrystValidateAllUserInput();
   wxFileDialog save(this,_T("Choose a file"),_T(""),_T(""),_T("*.txt"),wxFD_SAVE | wxFD_OVERWRITE_PROMPT);
   if(save.ShowModal() != wxID_OK) return;
   
   ofstream out(save.GetPath().ToAscii());
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
         wxMessageDialog dumbUser(this,_T("This Scattering Power is still used !"),
                                  _T("Whooops"),wxOK|wxICON_EXCLAMATION);
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
   Scatterer *scatt=0;
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
   if(event.GetId()== ID_CRYSTAL_MENU_SCATT_IMPORTATOMLIST)
   {
      wxFileDialog open(this,_T("Choose a file with a list of atoms: Element x y z occup"),_T(""),_T(""),_T("*"),
                        wxFD_OPEN | wxFD_FILE_MUST_EXIST);
      if(open.ShowModal() != wxID_OK) return;
      ifstream fin (open.GetPath().ToAscii());
      if(!fin)
      {
         throw ObjCrystException("WXCrystal::OnMenuAddScatterer() : Error opening file for input:"+string(open.GetPath().ToAscii()));
      }
      string symbol;
      REAL x,y,z,occup;
      int n=1;
      char buf [10];
      int scattPow;
      while(true)
      {
         fin>>symbol;
         if(fin.eof()) break;
         fin>>x>>y>>z>>occup;
         cout<<symbol<<n<<": "<<x<<", "<<y<<", "<<z<<endl;
         scattPow=mpCrystal->GetScatteringPowerRegistry().Find(symbol,"ScatteringPowerAtom",true);
         if(scattPow<0)
         {
            cout<<"Scattering power "<<symbol<<" not found, creating it..."<<endl;
            mpCrystal->AddScatteringPower(new ScatteringPowerAtom(symbol,symbol));
         }
         sprintf(buf,"%d",n++);
         mpCrystal->AddScatterer(new Atom(x,y,z,symbol+(string)buf,&(mpCrystal->GetScatteringPower(symbol)),occup));
         if(fin.eof()) break;
      }
      fin.close();
      scatt=0;
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
         wxTextEntryDialog bondLengthDialog(this,_T("Bond length"),
                                 _T("Enter bond length (Angstroems)"),_T("1"),wxOK | wxCANCEL);
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
         wxTextEntryDialog bondLengthDialog(this,_T("Bond length"),
                                 _T("Enter bond length (Angstroems)"),_T("1"),wxOK | wxCANCEL);
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
         wxTextEntryDialog bondLengthDialog(this,_T("Bond length"),
                                 _T("Enter bond length (Angstroems)"),_T("1"),wxOK | wxCANCEL);
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
         wxTextEntryDialog bondLengthDialog(this,_T("Bond length"),
                                 _T("Enter bond length (Angstroems)"),_T("1"),wxOK | wxCANCEL);
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
         wxTextEntryDialog bondLengthDialog(this,_T("Bond length"),
                                 _T("Enter bond length (Angstroems)"),_T("1"),wxOK | wxCANCEL);
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
         wxTextEntryDialog bondLengthDialog(this,_T("Bond length"),
                                 _T("Enter bond length (Angstroems)"),_T("1"),wxOK | wxCANCEL);
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
         wxTextEntryDialog bondLengthDialog(this,_T("Bond length"),
                                 _T("Enter bond length (Angstroems)"),_T("1"),wxOK | wxCANCEL);
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
         wxTextEntryDialog bondLengthDialog(this,_T("Bond length"),
                                 _T("Enter bond length (Angstroems)"),_T("1"),wxOK | wxCANCEL);
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
   if(scatt!=0) mpCrystal->AddScatterer(scatt);
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

void WXCrystal::OnMenuAtoms2Molecule(wxCommandEvent &event)
{
   vector<Atom*> v;
   for(unsigned int i=0; i<mpCrystal->GetScattererRegistry().GetNb();++i)
   {
      Atom *pAtom=dynamic_cast<Atom *>(&(mpCrystal->GetScattererRegistry().GetObj(i)));
      if(pAtom!=0) v.push_back(pAtom);
   }
   const unsigned int nb=v.size();
   wxString *choices = new wxString[nb];
   for(unsigned int i=0;i<nb;i++) 
      choices[i]= wxString::FromAscii((v[i]->GetName()).c_str());
   #if 0
   wxMultiChoiceDialog dialog (this,_T("Choose the molecule's atoms"),_T("Select Atoms"),nb,choices,wxOK | wxCANCEL);
   dialog.SetSize(300,300);
   #else
   wxMultiChoiceDialog_ListBox dialog(this,_T("Choose the molecule's atoms"),_T("Select Atoms"),nb,choices);
   #endif
   if(wxID_OK!=dialog.ShowModal()) return;
   wxArrayInt choice=dialog.GetSelections();
   if(choice.GetCount()>0)
   {
      list<Atom*> vChoice;
      for(unsigned int i=0;i<choice.GetCount();++i) vChoice.push_back(v[choice.Item(i)]);

      mpCrystal->AddScatterer(Atoms2Molecule(vChoice));
      for(unsigned int i=0;i<choice.GetCount();++i) mpCrystal->RemoveScatterer(v[choice.Item(i)]);
      mpCrystal->UpdateDisplay();
   }
}

void WXCrystal::OnMenuImportMoleculeFromFenskeHallZMatrix(wxCommandEvent &event)
{
   VFN_DEBUG_ENTRY("WXCrystal::OnMenuImportFenskeHallZMatrix()",6)
   WXCrystValidateAllUserInput();
   string tmp("Fenske-Hall z-matrix|*.fhz;*.fh");
   if(event.GetId()==ID_CRYSTAL_MENU_SCATT_IMPORTNAMEDZMATRIX) tmp="Fox z-matrix|*.zmat";
   wxFileDialog open(this,_T("Choose a file with a Fenske-Hall Z-matrix"),_T(""),_T(""), wxString::FromAscii(tmp.c_str()),
                     wxFD_OPEN | wxFD_FILE_MUST_EXIST);
   if(open.ShowModal() != wxID_OK) return;
   ifstream fin ( open.GetPath().ToAscii());
   if(!fin)
   {
      throw ObjCrystException("WXCrystal::OnMenuImportFenskeHallZMatrix() : \
Error opening file for input:"+string(open.GetPath().ToAscii()));
   }
   string filename(open.GetPath().ToAscii());
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
   bool named=false;
   if(event.GetId()==ID_CRYSTAL_MENU_SCATT_IMPORTNAMEDZMATRIX) named=true;
   scatt.ImportFenskeHallZMatrix(fin,named);
   fin.close();
   mpCrystal->AddScatterer(ZScatterer2Molecule(&scatt));
   this->CrystUpdate(true);
   VFN_DEBUG_EXIT("WXCrystal::OnMenuImportFenskeHallZMatrix()",6)
}

void WXCrystal::OnMenuSetRelativeXYZLimits(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_ENTRY("WXCrystal::OnMenuSetRelativeXYZLimits():Cancelled",6)
   WXCrystValidateAllUserInput();
   wxTextEntryDialog limitDialog(this,_T("Relative limits"),
                           _T("Enter relative limits for x,y,z (Angstroems)"),
                           _T("0.5"),wxOK | wxCANCEL);
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

/// Local class for a thread doing random moves to the structure
class TestCrystalThread: public wxThread
{
   public:
      TestCrystalThread(Crystal &cryst,float seconds):
         wxThread(wxTHREAD_DETACHED),mpCryst(&cryst),mSeconds(seconds){};
      virtual void *Entry()
      {
         cout<<endl<<"Entering refinement thread "<<endl<<endl;
         mpCryst->BeginOptimization();
         Chronometer chrono;
         float dt0=chrono.seconds();
         while(chrono.seconds()<30)
         {
            mpCryst->BeginGlobalOptRandomMove();
            mpCryst->GlobalOptRandomMove(0.05,gpRefParTypeObjCryst);
            wxMilliSleep(1);// Slow down display for simple structures
            if((chrono.seconds()-dt0)>0.05) {mpCryst->UpdateDisplay();dt0=chrono.seconds();}
         }
         mpCryst->EndOptimization();
         return NULL;
      };
      virtual void OnExit()
      {
         cout <<endl<<"Exiting refinement thread "<<endl<<endl;
      };
   private:
      /// The molecule to randomly
      Crystal *mpCryst;
      /// Test duration
      float mSeconds;
};

void WXCrystal::OnMenuTestRandomMoves(wxCommandEvent &event)
{
   TestCrystalThread *pTest = new TestCrystalThread(*mpCrystal,30);
   if(pTest->Create() != wxTHREAD_NO_ERROR) 
      wxLogError(_T("Can't create test optimization thread"));
   else pTest->Run();
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
      if(0!=mpCrystalGL) mpCrystalGL->GetParent()->SetLabel( wxString::FromAscii(mpCrystal->GetName().c_str()));
      #endif
      if(lock) mMutex.Unlock();
   }
   if(lock) mMutex.Lock();
   if(0!=mpScattPowWin)
   {
      map<ScatteringPowerAtom*,RowScattPow>::iterator pos;
      for(pos=mvpRowScattPow.begin();pos!=mvpRowScattPow.end();++pos)
      {
         if(pos->second.mNeedUpdateUI==true)
         {
            mIsSelfUpdating=true;
            mpScattPowWin->SetRowLabelValue(pos->second.mIdx, wxString::FromAscii(pos->second.mName.c_str()));
            wxString tmp;
            tmp.Printf(_T("%f"),pos->second.mBiso);
            mpScattPowWin->SetCellValue(pos->second.mIdx, 0, tmp);
            tmp.Printf(_T("%f"),pos->second.mFormalCharge);
            mpScattPowWin->SetCellValue(pos->second.mIdx, 1, tmp);
            tmp.Printf(_T("%f"),pos->second.mR);
            mpScattPowWin->SetCellValue(pos->second.mIdx, 2, tmp);
            tmp.Printf(_T("%f"),pos->second.mG);
            mpScattPowWin->SetCellValue(pos->second.mIdx, 3, tmp);
            tmp.Printf(_T("%f"),pos->second.mB);
            mpScattPowWin->SetCellValue(pos->second.mIdx, 4, tmp);
            tmp.Printf(_T("%f"),pos->second.mMaximumLikelihoodError);
            mpScattPowWin->SetCellValue(pos->second.mIdx, 5, tmp);
            tmp.Printf(_T("%f"),pos->second.mNbGhostAtoms);
            mpScattPowWin->SetCellValue(pos->second.mIdx, 6, tmp);
            mIsSelfUpdating=false;
         }
      }
   }
   if(0!=mpAntiBumpWin)
   {
      map<ScatteringPowerAtom*,RowScattPow>::iterator pos;
      for(pos=mvpRowScattPow.begin();pos!=mvpRowScattPow.end();++pos)
      {
         if(pos->second.mNeedUpdateUI==true)
         {
            mIsSelfUpdating=true;
            mpAntiBumpWin->SetRowLabelValue(pos->second.mIdx, wxString::FromAscii(pos->second.mName.c_str()));
            mpAntiBumpWin->SetColLabelValue(pos->second.mIdx, wxString::FromAscii(pos->second.mName.c_str()));
            wxString tmp;
            for(unsigned long j=0;j<pos->second.mvAntiBumpDistance.size();++j)
            {
               VFN_DEBUG_MESSAGE("WXCrystal::UpdateUI():Antibump("<<pos->first->GetName()
                                                                  <<",?"//<<mvScattPowRowIndex[j]->GetName()
                                                                  <<")="<<pos->second.mvAntiBumpDistance[j],3);
               if(pos->second.mvAntiBumpDistance[j]>-998)
               {
                  tmp.Printf(_T("%f"),pos->second.mvAntiBumpDistance[j]);
                  mpAntiBumpWin->SetCellValue(pos->second.mIdx,j,tmp);
               } else mpAntiBumpWin->SetCellValue(pos->second.mIdx,j,_T(""));
            }
            mIsSelfUpdating=false;
         }
      }
   }
   if(0!=mpBondValenceWin)
   {
      map<ScatteringPowerAtom*,RowScattPow>::iterator pos;
      for(pos=mvpRowScattPow.begin();pos!=mvpRowScattPow.end();++pos)
      {
         if(pos->second.mNeedUpdateUI==true)
         {
            mIsSelfUpdating=true;
            mpBondValenceWin->SetRowLabelValue(pos->second.mIdx, wxString::FromAscii(pos->second.mName.c_str()));
            mpBondValenceWin->SetColLabelValue(pos->second.mIdx, wxString::FromAscii(pos->second.mName.c_str()));
            wxString tmp;
            for(unsigned long j=0;j<pos->second.mvBondValenceRo.size();++j)
            {
               VFN_DEBUG_MESSAGE("WXCrystal::UpdateUI():BondValence("<<pos->first->GetName()
                                                                <<",?"//<<mvScattPowRowIndex[j]->GetName()
                                                               <<")="<<pos->second.mvBondValenceRo[j],3);
               if(pos->second.mvBondValenceRo[j]>-998)
               {
                  tmp.Printf(_T("%f"),pos->second.mvBondValenceRo[j]);
                  mpBondValenceWin->SetCellValue(pos->second.mIdx,j,tmp);
               } else mpBondValenceWin->SetCellValue(pos->second.mIdx,j,_T(""));
            }
            mIsSelfUpdating=false;
         }
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
      wxFrame *frame= new wxFrame(this,-1,_T("Scattering Powers parameters for: ")
                                  + wxString::FromAscii(this->GetCrystal().GetName().c_str()),
                                  wxDefaultPosition,wxSize(800,300));

      wxNotebook *notebook = new wxNotebook(frame, -1);
   {// Individual parameters
      mpScattPowWin = new WXCrystalScrolledGridWindow(notebook,this,ID_CRYSTAL_WIN_SCATTPOW);
      notebook->AddPage(mpScattPowWin, _T("Scattering Powers"), true);
      
      mpScattPowWin->SetDefaultRenderer(new wxGridCellFloatRenderer(5,3));
      mpScattPowWin->SetDefaultEditor(new wxGridCellFloatEditor(5,3));
      mpScattPowWin->SetColMinimalAcceptableWidth(150);
      mpScattPowWin->CreateGrid(0,7);
      
      mpScattPowWin->SetColLabelValue(0,_T("Biso"));
      mpScattPowWin->SetColLabelValue(1,_T("Charge"));
      mpScattPowWin->SetColLabelValue(2,_T("Red"));
      mpScattPowWin->SetColLabelValue(3,_T("Green"));
      mpScattPowWin->SetColLabelValue(4,_T("Blue"));
      mpScattPowWin->SetColLabelValue(5,_T("ML Error"));
      mpScattPowWin->SetColLabelValue(6,_T("#ghost"));
      
      mpScattPowWin->AutoSizeRows();
      mpScattPowWin->AutoSizeColumns();
   }
   {// Anti-Bump
      mpAntiBumpWin = new WXCrystalScrolledGridWindow(notebook,this,ID_CRYSTAL_WIN_ANTIBUMP);
      notebook->AddPage(mpAntiBumpWin, _T("AntiBump"), true);
      
      mpAntiBumpWin->SetDefaultRenderer(new wxGridCellFloatRenderer(5,3));
      mpAntiBumpWin->SetDefaultEditor(new wxGridCellFloatEditor(5,3));
      mpAntiBumpWin->SetColMinimalAcceptableWidth(150);
      mpAntiBumpWin->CreateGrid(0,0);
      
      mpAntiBumpWin->AutoSizeRows();
      mpAntiBumpWin->AutoSizeColumns();
   }
   {// Bond Valence
      mpBondValenceWin = new WXCrystalScrolledGridWindow(notebook,this,ID_CRYSTAL_WIN_BONDVALENCE);
      notebook->AddPage(mpBondValenceWin, _T("BondValence"), true);
      
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
   map<ScatteringPowerAtom*,RowScattPow>::iterator pos=mvpRowScattPow.begin();
   while(pos->second.mIdx!=r)++pos;
   ScatteringPowerAtom *const p=pos->first;
   
   wxString s=mpScattPowWin->GetCellValue(r,c);
   switch(c)
   {
      case 0:
      {
         if(s!=_T(""))
         {
            double d;
            s.ToDouble(&d);
            p->SetBiso(d);
         }
         break;
      }
      case 1:
      {
         if(s!=_T(""))
         {
            double d;
            s.ToDouble(&d);
            p->SetFormalCharge(d);
         }
         break;
      }
      case 2:
      {
         if(s!=_T(""))
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
         if(s!=_T(""))
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
         if(s!=_T(""))
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
         if(s!=_T(""))
         {
            double d;
            s.ToDouble(&d);
            p->SetMaximumLikelihoodPositionError(d);
         }
         break;
      }
      case 6:
      {
         if(s!=_T(""))
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

   map<ScatteringPowerAtom*,RowScattPow>::iterator pos=mvpRowScattPow.begin();
   while(pos->second.mIdx!=r)++pos;
   const ScatteringPowerAtom *const p1=pos->first;
   
   pos=mvpRowScattPow.begin();
   while(pos->second.mIdx!=c)++pos;
   const ScatteringPowerAtom *const p2=pos->first;
   
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

   map<ScatteringPowerAtom*,RowScattPow>::iterator pos=mvpRowScattPow.begin();
   while(pos->second.mIdx!=r)++pos;
   const ScatteringPowerAtom *const p1=pos->first;
   
   pos=mvpRowScattPow.begin();
   while(pos->second.mIdx!=c)++pos;
   const ScatteringPowerAtom *const p2=pos->first;
   
   wxString s=mpBondValenceWin->GetCellValue(r,c);
   double d;
   s.ToDouble(&d);
   if(d>0.01) mpCrystal->AddBondValenceRo(*p1,*p2,d);
   else mpCrystal->RemoveBondValenceRo(*p1,*p2);
   this->CrystUpdate(true,false);
}

void WXCrystal::NotifyDeleteListWin(WXCrystalScrolledGridWindow *win)
{
   if(win==mpScattPowWin) mpScattPowWin=0;
   if(win==mpAntiBumpWin) mpAntiBumpWin=0;
   if(win==mpBondValenceWin) mpBondValenceWin=0;
   // NOTE : all three subwindows should actually be deleted at the *same* time.
   if((mpScattPowWin==0)&&(mpAntiBumpWin==0)&&(mpBondValenceWin==0)) mvpRowScattPow.clear();
}

/*
void WXCrystal::OnMenuPDF(wxCommandEvent &event)
{
   const unsigned int nb=1000;
   // Simulate data
      if(mpPDF!=0) delete mpPDF;
      mpPDF=new PDF();
      CrystVector_REAL r,obs;
      r.resize(nb);obs.resize(nb);
      for(unsigned int i=0;i<nb;++i) r(i)=(i+1)*.02;
      obs=1.0;
      mpPDF->SetPDFObs(r,obs);
      PDFCrystal *pPDFCrystal=new PDFCrystal(*mpPDF,*mpCrystal);
      mpPDF->AddPDFPhase(*pPDFCrystal);
   // WX window
      wxFrame *frame= new wxFrame(this,-1,"PDF",
                                  wxDefaultPosition,wxSize(300,200));
      WXMultiGraph* pGraph =new WXMultiGraph(frame);

      wxSizer *ps=new wxBoxSizer(wxHORIZONTAL);
      ps->Add(pGraph,1,wxEXPAND);
      frame->CreateStatusBar(2);
      frame->SetSizer(ps);
      frame->SetAutoLayout(true);
      frame->Show(true);
      unsigned long id=pGraph->AddGraph("PDF");
      valarray<float> vr(nb),vcalc(nb);
      CrystVector_REAL v2r,v2calc;
      v2r=mpPDF->GetPDFR();
      v2calc=mpPDF->GetPDFCalc();
      for(unsigned int i=0;i<nb;++i)
      {
         vr[i]=v2r(i);
         vcalc[i]=v2calc(i);
      }
      pGraph->SetGraphData(id,vr,vcalc);
      pGraph->UpdateDisplay();
}
*/
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
//    UnitCellMap
//
////////////////////////////////////////////////////////////////////////
UnitCellMap::UnitCellMap(const Crystal&crystal):
mpCrystal(&crystal)
{}
UnitCellMap::~UnitCellMap(){}
void UnitCellMap::GLInitDisplayList(const float minValue,
					  WXGLCrystalCanvas * parentCrystal) const
{
   VFN_DEBUG_ENTRY("UnitCellMap::GLInitDisplayList()",10)
   //cout<<"Generating OpenGL Triangles for Fourier map:"<<mName<<", contour="<<minValue<<endl;
   // Generate triangles
      VFN_DEBUG_MESSAGE("UnitCellMap::GLInitDisplayList(): Generate Triangles",7)

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
      REAL x, y, z;

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
      VFN_DEBUG_MESSAGE("UnitCellMap::GLInitDisplayList(): MC, Min Value="<<minValue,10)
      const TRIANGLE *pTriangles= MC(snx-1, sny-1, snz-1, step[0], step[1], step[2], minValue, subPoints, numOfTriangles);
   // OpenGL drawing instructions
      VFN_DEBUG_MESSAGE("UnitCellMap::GLInitDisplayList(): OpenGL instructions",7)
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
   VFN_DEBUG_EXIT("UnitCellMap::GLInitDisplayList():nb triangles="<<numOfTriangles,10)
}

void UnitCellMap::POVRayDescription(ostream &os,const float minValue,
                                          const CrystalPOVRayOptions &options)const
{// basically the same code asGLInitDisplayList(), but creates cylinders
   VFN_DEBUG_ENTRY("UnitCellMap::POVRayDescription()",7)
   // Generate triangles
      VFN_DEBUG_MESSAGE("UnitCellMap::POVRayDescription(): Generate Triangles",7)

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
      REAL x, y, z;

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
      VFN_DEBUG_MESSAGE("UnitCellMap::POVRayDescription(): MC, Min Value="<<minValue,10)
      const TRIANGLE *pTriangles= MC(snx-1, sny-1, snz-1, step[0], step[1], step[2], minValue, subPoints, numOfTriangles);
   // drawing instructions
      VFN_DEBUG_MESSAGE("UnitCellMap::POVRayDescription(): POVRay instructions",7)
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
   VFN_DEBUG_EXIT("UnitCellMap::GLInitDisplayList()",7)
}

int UnitCellMap::ImportGRD(const string&filename)
{
   VFN_DEBUG_ENTRY("UnitCellMap::ImportGRD()",7)
   ifstream ffile(filename.c_str());
   if(!ffile.is_open())
   {     //if file could not be loaded for some reason then exit
     VFN_DEBUG_MESSAGE("UnitCellMap::ImportGRD() error opening "<<filename.c_str(),10)
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
   /*
   cout << "Min density value="<<mMin<<endl
        << "Max density value="<<mMax<<endl
        << "Mean density="<<mMean<<endl
        << "Standard Deviation="<<mStandardDeviation<<endl;
   */
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
   mType=3;
   VFN_DEBUG_EXIT("UnitCellMap::ImportGRD()",7)
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

int UnitCellMap::ImportDSN6(const string&filename)
{
   VFN_DEBUG_ENTRY("UnitCellMap::ImportDSN6()",7)
   FILE *pfile=fopen(filename.c_str(),"rb");
   if(NULL==pfile)
   {     //if file could not be loaded for some reason then exit
     VFN_DEBUG_MESSAGE("UnitCellMap::ImportDSN6() error opening "<<filename.c_str(),10)
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
   /*
   cout << "Min density value="<<mMin<<endl
        << "Max density value="<<mMax<<endl
        << "Mean density="<<mMean<<endl
        << "Standard Deviation="<<mStandardDeviation<<endl;
   */
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
         //cout<<"name="<<filename.substr(idx+1)<<endl;
         mName=filename.substr(idx+1);
      }
   }
   VFN_DEBUG_EXIT("UnitCellMap::ImportDSN6()",7)
   mType=3;
   return 1;
}
#ifdef HAVE_FFTW

/** Find the closest value toa given integer v which is larger than v and
* under the form 2^n2 * 3^n3 * 3^n5
*
*/
unsigned int closest235(unsigned int v)
{// This is not particularly fast but we need few iterations so...
   unsigned int n2=0,n3=0,n5=0;
   unsigned int v2,v3,v5=1;
   unsigned int bestdiff=10000;
   unsigned int best=0;
   for(unsigned int i5=1;;)
   {
      v3=1;
      for(unsigned int i3=1;;)
      {
         v2=1;
         for(unsigned int i2=1;;)
         {
            const unsigned int n=v2*v3*v5;
            if(n>v)
            {
               if((n-v)<bestdiff)
               {
                  bestdiff=n-v;
                  n2=i2;
                  n3=i3;
                  n5=i5;
                  best=n;
                  if(best==v) return best;
               }
               else break;
            }
            v2*=2;i2+=1;
         }
         v3*=3;i3+=1;
         if((v3*v5)>(v+bestdiff)) break;
      }
      v5*=5;i5+=1;
      if(v5>(v+bestdiff)) break;
   }
   return best;
}

int UnitCellMap::CalcFourierMap(const ScatteringData& data, unsigned int type0, const bool normalized_sf)
{
   mpData=&data;
   const float resolution=0.3;//Approximate resolution in Ansgtroem
   // We need something like 2^n2 * 3^n3 * 5^n5 - just use a power of 2 now
   const unsigned long sizex=closest235((unsigned int)floor(mpCrystal->GetLatticePar(0)/resolution+.5)) ;//int(pow((double)2, (double)ceil(log(mpCrystal->GetLatticePar(0)/resolution)/log(2)))+.00001);
   const unsigned long sizey=closest235((unsigned int)floor(mpCrystal->GetLatticePar(1)/resolution+.5)) ;//int(pow((double)2, (double)ceil(log(mpCrystal->GetLatticePar(1)/resolution)/log(2)))+.00001);
   const unsigned long sizez=closest235((unsigned int)floor(mpCrystal->GetLatticePar(2)/resolution+.5)) ;//int(pow((double)2, (double)ceil(log(mpCrystal->GetLatticePar(2)/resolution)/log(2)))+.00001);
   //cout<<"UnitCellMap::CalcFourierMap():"<<sizex<<","<<sizey<<","<<sizez<<","<<endl;
   fftwf_complex *in= (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * sizex*sizey*sizez);
   fftwf_plan plan=fftwf_plan_dft_3d(sizez, sizey, sizex,in, in,FFTW_FORWARD, FFTW_ESTIMATE);
   
   float *p=(float*)in;
   for(unsigned long i=0;i<sizex*sizey*sizez*2;i++) *p++=0;
   
   const long nb=data.GetNbReflBelowMaxSinThetaOvLambda();
   
   mType=type0;
   if(data.GetFhklObsSq().numElements()==0) mType=1;
   
   CrystVector_REAL norm_sf;
   if(normalized_sf)
   {
      CrystVector_REAL tmp;
      norm_sf.resize(data.GetFhklCalcReal().numElements());
      norm_sf=0;
      const map<const ScatteringPower*,CrystVector_REAL> *pSF=&(data.GetScatteringFactor());
      const ScatteringComponentList *pComp =&(mpCrystal->GetScatteringComponentList());
      REAL norm0=0;// norm_sf normalized to 1 at low angle
      for(unsigned int i=0;i<pComp->GetNbComponent();i++)
      {
         tmp=pSF->find((*pComp)(i).mpScattPow)->second;// safe enough ?
         tmp*=tmp;
         tmp*=  (*pComp)(i).mOccupancy * (*pComp)(i).mDynPopCorr;
         
         const REAL sf0=(*pComp)(i).mpScattPow->GetForwardScatteringFactor(data.GetRadiationType ());
         norm0+=(*pComp)(i).mOccupancy * (*pComp)(i).mDynPopCorr *sf0*sf0;
         
         norm_sf+=tmp;
      }
      REAL *p=norm_sf.data();
      norm0=1/norm0;
      for(unsigned int i=norm_sf.numElements();i>0;i--) {*p=sqrt(*p * norm0);p++;}
   }
   
   // Auto-scale Fcalc and Fobs ?
   REAL scale_fobs=1.0;
   if(mType!=1)
   {
      REAL tmp=0;
      scale_fobs=0;
      for(long i=0;i<nb;++i) {scale_fobs+=data.GetFhklCalcSq()(i); tmp+=data.GetFhklObsSq()(i);}
      scale_fobs=sqrt(scale_fobs/(tmp+1e-10));
      //cout<<__FILE__<<":"<<__LINE__<<" Fourier map obs/calc scale factor:"<<scale_fobs<<endl;
   }

   const REAL v=1/mpCrystal->GetVolume();//(REAL)(size*size*size);// mpCrystal->GetVolume(); (REAL)(size*size*size);
   for(long i=0;i<nb;++i)
   {
      CrystMatrix_REAL m=mpCrystal->GetSpaceGroup().GetAllEquivRefl (data.GetH()(i),data.GetK()(i),data.GetL()(i),
                                                                     false, data.IsIgnoringImagScattFact(),
                                                                     data.GetFhklCalcReal()(i),data.GetFhklCalcImag()(i));
      REAL norm=1.0;
      if(normalized_sf) norm=1/norm_sf(i);
      for(int j=0;j<m.rows();j++)
      {
         int h=int(m(j,0)),k=int(m(j,1)),l=int(m(j,2));
         if((abs(h*2)>sizex)||(abs(k*2)>sizey)||(abs(l*2)>sizez)) continue;
         h=(h+sizex)%sizex;// e.g. h=-1 is at nx-1
         k=(k+sizey)%sizey;
         l=(l+sizez)%sizez;
         /*
         cout <<int(m(j,0))<<" "<<int(m(j,1))<<" "<<int(m(j,2))<<"("
              <<mpCrystal->GetSpaceGroup().IsReflCentric(data.GetH()(i),data.GetK()(i),data.GetL()(i))<<"):"
              <<m(j,3)<<"+"<<m(j,4)<<"i"<<endl;
         */
         if(mType==2)
         {// Obs-Calc
            const REAL fobs=scale_fobs*sqrt(fabs(data.GetFhklObsSq()(i)));
            const REAL rec=m(j,3),imc=m(j,4),fcalc=sqrt(fabs(data.GetFhklCalcSq()(i)));
            in[h+sizex*k+sizex*sizey*l][0]=v*rec*(fobs-fcalc)/sqrt(rec*rec+imc*imc)*norm;
            in[h+sizex*k+sizex*sizey*l][1]=v*imc*(fobs-fcalc)/sqrt(rec*rec+imc*imc)*norm;
         }
         if(mType==1)
         {// Calc
            in[h+sizex*k+sizex*sizey*l][0]=v*m(j,3)*norm;
            in[h+sizex*k+sizex*sizey*l][1]=v*m(j,4)*norm;
         }
         if(mType==0)
         {// Obs
            const REAL iobs=scale_fobs*sqrt(fabs(data.GetFhklObsSq()(i)));
            const REAL rec=m(j,3),imc=m(j,4),icalc=sqrt(fabs(data.GetFhklCalcSq()(i)));
            in[h+sizex*k+sizex*sizey*l][0]=v*rec*iobs/icalc*norm;
            in[h+sizex*k+sizex*sizey*l][1]=v*imc*iobs/icalc*norm;
         }
      }
      //cout<<endl;
   }
   
   if(mType!=2)
   {// F000, for obs & calc fourier maps ?
      const int nbSymmetrics=mpCrystal->GetSpaceGroup().GetNbSymmetrics(false,false);
      const ScatteringComponentList *pScattCompList=&(mpCrystal->GetScatteringComponentList());
      const long nbComp=pScattCompList->GetNbComponent();
      for(long i=0;i<nbComp;i++)
      {
         //TODO: include f" en forward scattering factor ?
         in[0][0]+= (*pScattCompList)(i).mpScattPow->GetForwardScatteringFactor(data.GetRadiationType())
                   *(*pScattCompList)(i).mOccupancy
                   *(*pScattCompList)(i).mDynPopCorr
                   *nbSymmetrics*v;
      }
      //cout<<"F(000)="<<in[0][0]/v<<endl;
   }
   fftwf_execute(plan);
   mPoints.resize(sizez,sizey,sizex);
   REAL *p1=mPoints.data();
   for(unsigned int i=0;i<sizex*sizey*sizez;i++) *p1++ =in[i][0] ;
   /*
   for(unsigned int ix=0;ix<sizex;ix++)
      for(unsigned int iy=0;iy<sizey;iy++)
         for(unsigned int iz=0;iz<sizez;iz++)
            mPoints(iz,iy,ix)=in[ix+sizex*iy+sizex*sizey*iz][0];
   */
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
   /*
   cout << "Min density value="<<mMin<<endl
        << "Max density value="<<mMax<<endl
        << "Mean density="<<mMean<<endl
        << "Standard Deviation="<<mStandardDeviation<<endl;
   */
   fftwf_destroy_plan(plan);
   fftwf_free(in);
   
   mName=data.GetClassName()+":";
   if(data.GetName()=="") mName+="?";
   else mName+=data.GetName();
   if(data.GetClassName()=="PowderPatternDiffraction")
   {
      mName="P:";
      if(data.GetRadiationType()==RAD_XRAY)     mName+="Xray:";
      if(data.GetRadiationType()==RAD_NEUTRON)  mName+="Neut:";
      if(data.GetRadiationType()==RAD_ELECTRON) mName+="Elec:";
      
      char buf[100];
      if(data.GetRadiation().GetWavelengthType()==WAVELENGTH_TOF) mName+="TOF:";
      else
      {
         sprintf(buf,"%6.3fA:",data.GetWavelength()(0));
         mName+=buf;
      }
      const PowderPatternDiffraction* diff=dynamic_cast<const PowderPatternDiffraction *>(&data);
      if(diff!=0) mName+=diff->GetParentPowderPattern().GetName();
   }
   if(data.GetClassName()=="DiffractionDataSingleCrystal")
   {
      mName="S:";
      if(data.GetRadiationType()==RAD_XRAY)     mName+="Xray:";
      if(data.GetRadiationType()==RAD_NEUTRON)  mName+="Neut:";
      if(data.GetRadiationType()==RAD_ELECTRON) mName+="Elec:";
      
      char buf[100];
      if(data.GetRadiation().GetWavelengthType()==WAVELENGTH_TOF) mName+="TOF:";
      else
      {
         sprintf(buf,"%6.3fA:",data.GetWavelength()(0));
         mName+=buf;
      }
      mName+=data.GetName();
      
   }
   if(mType==0) mName="(Fo)"+mName;
   if(mType==1) mName="(Fc)"+mName;
   if(mType==2) mName="(Fo-Fc)"+mName;
   return 1;
}
#endif

const string & UnitCellMap::GetName()const
{
   return mName;
}

REAL UnitCellMap::GetValue(const REAL x,const REAL y,const REAL z)const
{
   const int nx=mPoints.cols();
   const int ny=mPoints.rows();
   const int nz=mPoints.depth();
   long ix=((long)floor(x*nx+.5))%nx;
   long iy=((long)floor(y*ny+.5))%ny;
   long iz=((long)floor(z*nz+.5))%nz;
   if(ix<0) ix+=nx;
   if(iy<0) iy+=ny;
   if(iz<0) iz+=nz;
   return mPoints(iz,iy,ix);
}
REAL UnitCellMap::Max()const{return mMax;}
REAL UnitCellMap::Min()const{return mMin;}
REAL UnitCellMap::Mean()const{return mMean;}
REAL UnitCellMap::StandardDeviation()const{return mStandardDeviation;}
int UnitCellMap::GetType()const{return mType;}
const Crystal &UnitCellMap::GetCrystal()const{return *mpCrystal;}
const ScatteringData *UnitCellMap::GetData()const{return mpData;}

////////////////////////////////////////////////////////////////////////
//
//    UnitCellMapGLList
//
////////////////////////////////////////////////////////////////////////
UnitCellMapGLList::UnitCellMapGLList(const UnitCellMap &ucmap,WXGLCrystalCanvas * parent,
                                     const bool showWire,float contour,
                                     const float r,const float g,const float b,const float t):
mGLDisplayList(0),mShowWire(showWire),mShow(true),mContour(contour),mpUCMap(&ucmap),mpParent(parent)
{
   VFN_DEBUG_MESSAGE("UnitCellMapGLList::UnitCellMapGLList()",10)
   this->SetColour(r,g,b,t);
}

UnitCellMapGLList::~UnitCellMapGLList()
{
   VFN_DEBUG_MESSAGE("UnitCellMapGLList::~UnitCellMapGLList()",10)
   if(0!=mGLDisplayList) glDeleteLists(mGLDisplayList,1);
}
void UnitCellMapGLList::GenList()
{
   VFN_DEBUG_ENTRY("UnitCellMapGLList::GenList()",7)
   if(0==mGLDisplayList) mGLDisplayList=glGenLists(1);
   glNewList(mGLDisplayList,GL_COMPILE);
      glPushMatrix();
      mpUCMap->GLInitDisplayList(mContour, mpParent);
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

      const GLfloat colorBack [] = {mColour[0]/3., mColour[1]/3., mColour[2]/3., mColour[3]};
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
void UnitCellMapGLList::SetShow(bool show) {mShow=show;}
bool UnitCellMapGLList::Show()const {return mShow;}
void UnitCellMapGLList::SetContour(float contour) {mContour=contour;}
float UnitCellMapGLList::GetContour()const {return mContour;}
const UnitCellMap & UnitCellMapGLList::GetMap()const {return *mpUCMap;}
////////////////////////////////////////////////////////////////////////
//
//    WXGLCrystalCanvas::WXFourierMapList
//
////////////////////////////////////////////////////////////////////////
static const long ID_GLCRYSTAL_FOURIER_ADD=            WXCRYST_ID();
static const long ID_GLCRYSTAL_FOURIER_REMOVE=         WXCRYST_ID();
static const long ID_GLCRYSTAL_FOURIER_UPDATE=         WXCRYST_ID();
static const long ID_GLCRYSTAL_FOURIER_WIREFRAME=      WXCRYST_ID();
static const long ID_GLCRYSTAL_FOURIER_SHOW=           WXCRYST_ID();
static const long ID_GLCRYSTAL_FOURIER_SHARPEN=        WXCRYST_ID();
static const long ID_GLCRYSTAL_FOURIER_LISTMAP=        WXCRYST_ID();
static const long ID_GLCRYSTAL_FOURIER_LISTGLMAP=      WXCRYST_ID();
static const long ID_GLCRYSTAL_FOURIER_CONTOUR=        WXCRYST_ID();
static const long ID_GLCRYSTAL_FOURIER_NEWCONTOUR=     WXCRYST_ID();
static const long ID_GLCRYSTAL_FOURIER_COLOURPICKER=   WXCRYST_ID();

WXGLCrystalCanvas::WXFourierMapList::WXFourierMapList(WXGLCrystalCanvas *pGLCrystalCanvas,wxWindow *parent):
wxWindow(parent,-1),mpGLCrystalCanvas(pGLCrystalCanvas),mIsUpdating(false)
{
   this->SetFont(wxFont(8,wxTELETYPE,wxFONTSTYLE_NORMAL,wxFONTWEIGHT_NORMAL));
   wxBoxSizer* pSizer=new wxBoxSizer(wxVERTICAL);
   // Top buttons
      wxBoxSizer* pSizerButtons=new wxBoxSizer(wxHORIZONTAL);
      wxButton *pButtonUpdate=new wxButton(this,ID_GLCRYSTAL_FOURIER_UPDATE,_T("Update 3D View"));
      mpWireFrame=new wxCheckBox(this,ID_GLCRYSTAL_FOURIER_WIREFRAME,_T("Wireframe"));
      mpShowFourier=new wxCheckBox(this,ID_GLCRYSTAL_FOURIER_SHOW,_T("Show Fourier"));
      mpSharpenMap=new wxCheckBox(this,ID_GLCRYSTAL_FOURIER_SHARPEN,_T("Sharpen maps"));
      pSizerButtons->Add(pButtonUpdate,0,wxALIGN_CENTER);
      pSizerButtons->Add(mpWireFrame,0,wxALIGN_CENTER);
      pSizerButtons->Add(mpShowFourier,0,wxALIGN_CENTER);
      pSizerButtons->Add(mpSharpenMap,0,wxALIGN_CENTER);
      pSizer->Add(pSizerButtons,0,wxALIGN_CENTER);
   
   // Map lists
      wxBoxSizer* pSizerMaps=new wxBoxSizer(wxHORIZONTAL);
      
      // Left column - available maps
         wxBoxSizer* pSizerLeft=new wxBoxSizer(wxVERTICAL);
         pSizerMaps->Add(pSizerLeft,0,wxALIGN_TOP);
         
         wxStaticText *mpLabel0=new wxStaticText(this,-1,_T("Available Maps"));
         pSizerLeft->Add(mpLabel0,0,wxALIGN_CENTER);
         mpAvailableMapList=new wxListBox(this,ID_GLCRYSTAL_FOURIER_LISTMAP,wxDefaultPosition,wxSize(400,150));
         pSizerLeft->Add(mpAvailableMapList,0,wxALIGN_CENTER);
         
         mpMapInfo=new wxStaticText(this,-1,_T("min=+00.00 max=+00.00 sigma=00.00"));
         pSizerLeft->Add(mpMapInfo,0,wxALIGN_CENTER);
         
         wxBoxSizer* pSizerLeft2=new wxBoxSizer(wxHORIZONTAL);
         pSizerLeft->Add(pSizerLeft2,0,wxALIGN_CENTER);
         wxStaticText *mpLabel2=new wxStaticText(this,-1,_T("New Contour:"));
         mpNewContourValue=new wxTextCtrl(this,ID_GLCRYSTAL_FOURIER_NEWCONTOUR,_T(""),wxDefaultPosition,wxDefaultSize,wxTE_PROCESS_ENTER);
         pSizerLeft2->Add(mpLabel2,0,wxALIGN_CENTER);
         pSizerLeft2->Add(mpNewContourValue,0,wxALIGN_CENTER);
      
         wxButton *pButtonAdd=new wxButton(this,ID_GLCRYSTAL_FOURIER_ADD,_T("Add"));
         pSizerLeft->Add(pButtonAdd,0,wxALIGN_CENTER);
      
      pSizerMaps->AddSpacer(5);
      // Right column - displayed maps & contours
         wxBoxSizer* pSizerRight=new wxBoxSizer(wxVERTICAL);
         pSizerMaps->Add(pSizerRight,0,wxALIGN_TOP);
         
         wxStaticText *mpLabel0r=new wxStaticText(this,-1,_T("Displayed Maps"));
         pSizerRight->Add(mpLabel0r,0,wxALIGN_CENTER);
         mpDisplayedMapList=new wxListBox(this,ID_GLCRYSTAL_FOURIER_LISTGLMAP,wxDefaultPosition,wxSize(400,150));
         pSizerRight->Add(mpDisplayedMapList,0,wxALIGN_CENTER);
         
         wxBoxSizer* pSizerRight1=new wxBoxSizer(wxHORIZONTAL);
         pSizerRight->Add(pSizerRight1,0,wxALIGN_CENTER);
         wxStaticText *mpLabel3=new wxStaticText(this,-1,_T("Contour:"));
         mpContourValue=new wxTextCtrl(this,ID_GLCRYSTAL_FOURIER_CONTOUR,_T(""),wxDefaultPosition,wxDefaultSize,wxTE_PROCESS_ENTER);
         pSizerRight1->Add(mpLabel3,0,wxALIGN_CENTER);
         pSizerRight1->Add(mpContourValue,0,wxALIGN_CENTER);
         
         mpColourPicker=new wxColourPickerCtrl(this, ID_GLCRYSTAL_FOURIER_COLOURPICKER, *wxRED, wxDefaultPosition, wxDefaultSize,wxCLRP_USE_TEXTCTRL);
         wxButton *pButtonRemove=new wxButton(this,ID_GLCRYSTAL_FOURIER_REMOVE,_T("Remove"));
         pSizerRight->Add(mpColourPicker,0,wxALIGN_CENTER);
         pSizerRight->Add(pButtonRemove,0,wxALIGN_CENTER);
      pSizer->Add(pSizerMaps,0,wxALIGN_CENTER);
   this->SetSizer(pSizer);
   this->SetAutoLayout(true);
   pSizer->SetSizeHints(this);
   pSizer->SetSizeHints(parent);
   this->Layout();
}
WXGLCrystalCanvas::WXFourierMapList::~WXFourierMapList()
{
   mpGLCrystalCanvas->NotifyDeleteFourierWin();
}

struct GLCrystalConfig
{
  GLCrystalConfig(const bool saved=false);
  float mDist;
  REAL mX0, mY0,mZ0;
  float mViewAngle;
  float mQuat [4];
  bool mSaved;
  bool mShowAtomName;
  bool mShowCursor;
  BBox mcellbbox;
  BBox mmapbbox;
  Triple mViewCntr;
};

GLCrystalConfig::GLCrystalConfig(const bool saved)
{
  mSaved=saved;
}

static GLCrystalConfig sGLCrystalConfig;

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
static const long ID_GLCRYSTAL_MENU_FOURIER=           WXCRYST_ID();
static const long ID_GLCRYSTAL_MENU_LOADFOURIERGRD=    WXCRYST_ID();
static const long ID_GLCRYSTAL_MENU_LOADFOURIERDSN6=   WXCRYST_ID();
//static const long ID_GLCRYSTAL_MENU_UNLOADFOURIER=     WXCRYST_ID();
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
//   EVT_MENU             (ID_GLCRYSTAL_MENU_UNLOADFOURIER,       WXGLCrystalCanvas::OnUnloadFourier)
   EVT_MENU             (ID_GLCRYSTAL_MENU_POVRAY,              WXGLCrystalCanvas::OnPOVRay)
   EVT_MENU             (ID_GLCRYSTAL_MENU_FOURIER,             WXGLCrystalCanvas::OnFourier)
   EVT_LISTBOX          (ID_GLCRYSTAL_FOURIER_LISTGLMAP,        WXGLCrystalCanvas::OnFourier)
   EVT_LISTBOX          (ID_GLCRYSTAL_FOURIER_LISTMAP,          WXGLCrystalCanvas::OnFourier)
   EVT_BUTTON           (ID_GLCRYSTAL_FOURIER_ADD,              WXGLCrystalCanvas::OnFourier)
   EVT_BUTTON           (ID_GLCRYSTAL_FOURIER_REMOVE,           WXGLCrystalCanvas::OnFourier)
   EVT_BUTTON           (ID_GLCRYSTAL_FOURIER_UPDATE,           WXGLCrystalCanvas::OnFourier)
   EVT_CHECKBOX         (ID_GLCRYSTAL_FOURIER_WIREFRAME,        WXGLCrystalCanvas::OnFourier)
   EVT_CHECKBOX         (ID_GLCRYSTAL_FOURIER_SHOW,             WXGLCrystalCanvas::OnFourier)
   EVT_CHECKBOX         (ID_GLCRYSTAL_FOURIER_SHARPEN,          WXGLCrystalCanvas::OnFourier)
   EVT_TEXT_ENTER       (ID_GLCRYSTAL_FOURIER_NEWCONTOUR,       WXGLCrystalCanvas::OnFourier)
   EVT_TEXT_ENTER       (ID_GLCRYSTAL_FOURIER_CONTOUR,          WXGLCrystalCanvas::OnFourier)
   EVT_COLOURPICKER_CHANGED(ID_GLCRYSTAL_FOURIER_COLOURPICKER,  WXGLCrystalCanvas::OnFourierChangeColour)
   EVT_CHAR             (WXGLCrystalCanvas::OnKeyDown)
   EVT_KEY_DOWN         (WXGLCrystalCanvas::OnKeyDown)
   EVT_KEY_UP           (WXGLCrystalCanvas::OnKeyUp)
   EVT_UPDATE_UI(ID_GLCRYSTAL_UPDATEUI,WXGLCrystalCanvas::OnUpdateUI)
END_EVENT_TABLE()

int AttribList [] = {WX_GL_RGBA , WX_GL_DOUBLEBUFFER, WX_GL_DEPTH_SIZE, 16,0};

WXGLCrystalCanvas::WXGLCrystalCanvas(WXCrystal *wxcryst,
                                     wxFrame *parent, wxWindowID id,
                                     const wxPoint &pos,
                                     const wxSize &size):
wxGLCanvas(parent, id,AttribList,pos,size,wxDEFAULT_FRAME_STYLE | wxFULL_REPAINT_ON_RESIZE,_T("GLCanvas"),wxNullPalette),
//wxGLCanvas(parent,id,pos,size,wxDEFAULT_FRAME_STYLE,_T("GLCanvas"),AttribList),
mpParentFrame(parent),
mpWXCrystal(wxcryst),mIsGLInit(false),mDist(60),mX0(0),mY0(0),mZ0(0),mViewAngle(15),
mShowFourier(true),mShowCrystal(true),mShowAtomName(true),mShowCursor(false),mSharpenMap(true),
mIsGLFontBuilt(false),mGLFontDisplayListBase(0),mpFourierMapListWin(0)
{
   mpwxGLContext=new wxGLContext(this);
   VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::WXGLCrystalCanvas()",3)
   if(sGLCrystalConfig.mSaved)
   {
      mDist=sGLCrystalConfig.mDist;
      mX0=sGLCrystalConfig.mX0;
      mY0=sGLCrystalConfig.mY0;
      mZ0=sGLCrystalConfig.mZ0;
      mViewAngle=sGLCrystalConfig.mViewAngle;
      for(int i=0;i<4;++i) mQuat[i]=sGLCrystalConfig.mQuat[i];
      mShowAtomName=sGLCrystalConfig.mShowAtomName;
      mShowCursor=sGLCrystalConfig.mShowCursor;
      
      mcellbbox.xMin=sGLCrystalConfig.mcellbbox.xMin;
      mcellbbox.xMax=sGLCrystalConfig.mcellbbox.xMax;
      mcellbbox.yMin=sGLCrystalConfig.mcellbbox.yMin;
      mcellbbox.yMax=sGLCrystalConfig.mcellbbox.yMax;
      mcellbbox.zMin=sGLCrystalConfig.mcellbbox.zMin;
      mcellbbox.zMax=sGLCrystalConfig.mcellbbox.zMax;
      
      mmapbbox.xMin=sGLCrystalConfig.mmapbbox.xMin;
      mmapbbox.xMax=sGLCrystalConfig.mmapbbox.xMax;
      mmapbbox.yMin=sGLCrystalConfig.mmapbbox.yMin;
      mmapbbox.yMax=sGLCrystalConfig.mmapbbox.yMax;
      mmapbbox.zMin=sGLCrystalConfig.mmapbbox.zMin;
      mmapbbox.zMax=sGLCrystalConfig.mmapbbox.zMax;
      
      mViewCntr.x=sGLCrystalConfig.mViewCntr.x;
      mViewCntr.y=sGLCrystalConfig.mViewCntr.y;
      mViewCntr.z=sGLCrystalConfig.mViewCntr.z;
   }
   else
   {
      mcellbbox.xMin = mpWXCrystal->GetCrystal().GetSpaceGroup().GetAsymUnit().Xmin()-0.1;
      mcellbbox.yMin = mpWXCrystal->GetCrystal().GetSpaceGroup().GetAsymUnit().Ymin()-0.1;
      mcellbbox.zMin = mpWXCrystal->GetCrystal().GetSpaceGroup().GetAsymUnit().Zmin()-0.1;
      mcellbbox.xMax = mpWXCrystal->GetCrystal().GetSpaceGroup().GetAsymUnit().Xmax()+0.1;
      mcellbbox.yMax = mpWXCrystal->GetCrystal().GetSpaceGroup().GetAsymUnit().Ymax()+0.1;
      mcellbbox.zMax = mpWXCrystal->GetCrystal().GetSpaceGroup().GetAsymUnit().Zmax()+0.1;
   }
     // N.B. xMin=xMax so that the previous cell bbox is used for Maps 
     // until mmapbbox is changed
   mmapbbox.xMin = mmapbbox.xMax = mmapbbox.yMin = mmapbbox.zMin = 0.;
   mmapbbox.yMax = mmapbbox.zMax = 1.;
   mpPopUpMenu=new wxMenu(_T("Crystal"));
   mpPopUpMenu->Append(ID_GLCRYSTAL_MENU_UPDATE, _T("&Update"));
   mpPopUpMenu->AppendSeparator();
   mpPopUpMenu->Append(ID_GLCRYSTAL_MENU_CHANGELIMITS, _T("Change display &Limits"));
   mpPopUpMenu->Append(ID_GLCRYSTAL_MENU_LIMITS_FULLCELL, _T("Show Full Unit Cell +0.1"));
   mpPopUpMenu->Append(ID_GLCRYSTAL_MENU_LIMITS_ASYMCELL, _T("Show Asymmetric Unit Cell +0.1"));
   mpPopUpMenu->AppendSeparator();
   mpPopUpMenu->Append(ID_GLCRYSTAL_MENU_SHOWCRYSTAL, _T("Hide Crystal"));
   mpPopUpMenu->Append(ID_GLCRYSTAL_MENU_SHOWATOMLABEL, _T("Hide Atom Labels"));
   mpPopUpMenu->Append(ID_GLCRYSTAL_MENU_SHOWCURSOR, _T("Show Cursor"));
   mpPopUpMenu->Append(ID_GLCRYSTAL_MENU_SETCURSOR, _T("Set view cntr and cursor pos."));
   mpPopUpMenu->Append(ID_GLCRYSTAL_MENU_POVRAY, _T("Create POVRay file"));
   mpPopUpMenu->AppendSeparator();
   mpPopUpMenu->Append(ID_GLCRYSTAL_MENU_FOURIER, _T("Fourier Maps"));
   mpPopUpMenu->Append(ID_GLCRYSTAL_MENU_LOADFOURIERGRD, _T("Load GRD Fourier Map"));	
   mpPopUpMenu->Append(ID_GLCRYSTAL_MENU_LOADFOURIERDSN6,_T("Load DSN6 Fourier Map"));	
   /*
   mpPopUpMenu->Append(ID_GLCRYSTAL_MENU_UNLOADFOURIER, "Unload Fourier Map(s)");
   mpPopUpMenu->Enable(ID_GLCRYSTAL_MENU_UNLOADFOURIER, FALSE);
   */
   if(sGLCrystalConfig.mSaved==false)
   {
      if(!wxConfigBase::Get()->HasEntry(_T("Crystal/BOOL/Default-display only asymmetric unit cell in 3D view")))
          wxConfigBase::Get()->Write(_T("Crystal/BOOL/Default-display only asymmetric unit cell in 3D view"), true);
      else
      {
          bool val;
          wxConfigBase::Get()->Read(_T("Crystal/BOOL/Default-display only asymmetric unit cell in 3D view"), &val);
          if(val)
          {
            mcellbbox.xMin = mpWXCrystal->GetCrystal().GetSpaceGroup().GetAsymUnit().Xmin()-0.1;
            mcellbbox.yMin = mpWXCrystal->GetCrystal().GetSpaceGroup().GetAsymUnit().Ymin()-0.1;
            mcellbbox.zMin = mpWXCrystal->GetCrystal().GetSpaceGroup().GetAsymUnit().Zmin()-0.1;
            mcellbbox.xMax = mpWXCrystal->GetCrystal().GetSpaceGroup().GetAsymUnit().Xmax()+0.1;
            mcellbbox.yMax = mpWXCrystal->GetCrystal().GetSpaceGroup().GetAsymUnit().Ymax()+0.1;
            mcellbbox.zMax = mpWXCrystal->GetCrystal().GetSpaceGroup().GetAsymUnit().Zmax()+0.1;
          }
          else
          {
            mcellbbox.xMin = -0.1;
            mcellbbox.yMin = -0.1;
            mcellbbox.zMin = -0.1;
            mcellbbox.xMax =  1.1;
            mcellbbox.yMax =  1.1;
            mcellbbox.zMax =  1.1;
          }
      }
      if(!wxConfigBase::Get()->HasEntry(_T("Crystal/BOOL/Default-display atom names in 3D view")))
          wxConfigBase::Get()->Write(_T("Crystal/BOOL/Default-display atom names in 3D view"), mShowAtomName);
      else
      {
          wxConfigBase::Get()->Read(_T("Crystal/BOOL/Default-display atom names in 3D view"), &mShowAtomName);
      }
   }
   if(mShowAtomName) mpPopUpMenu->SetLabel(ID_GLCRYSTAL_MENU_SHOWATOMLABEL, _T("Hide Atom Labels"));
   else mpPopUpMenu->SetLabel(ID_GLCRYSTAL_MENU_SHOWATOMLABEL, _T("Show Atom Labels"));
}

WXGLCrystalCanvas::~WXGLCrystalCanvas()
{
   mpWXCrystal->NotifyCrystalGLDelete();
   #ifndef HAVE_GLUT
   this->DeleteGLFont();
   #endif
/*   cout<<"WXGLCrystalCanvas:Store window pos&size:"<<this->GetParent()->GetId()
       <<":"<<this->GetParent()->GetPosition().x<<","
       <<":"<<this->GetParent()->GetPosition().y<<"-"
       <<":"<<this->GetParent()->GetSize().x<<","
       <<":"<<this->GetParent()->GetSize().y<<endl;*/
   gvWindowPosition[this->GetParent()->GetId()]=make_pair(this->GetParent()->GetPosition(),this->GetParent()->GetSize());
   sGLCrystalConfig.mDist=mDist;
   sGLCrystalConfig.mX0=mX0;
   sGLCrystalConfig.mY0=mY0;
   sGLCrystalConfig.mZ0=mZ0;
   sGLCrystalConfig.mViewAngle=mViewAngle;
   for(int i=0;i<4;++i) sGLCrystalConfig.mQuat[i]=mQuat[i];
   sGLCrystalConfig.mShowAtomName=mShowAtomName;
   sGLCrystalConfig.mShowCursor=mShowCursor;
    
   sGLCrystalConfig.mcellbbox.xMin=mcellbbox.xMin;
   sGLCrystalConfig.mcellbbox.xMax=mcellbbox.xMax;
   sGLCrystalConfig.mcellbbox.yMin=mcellbbox.yMin;
   sGLCrystalConfig.mcellbbox.yMax=mcellbbox.yMax;
   sGLCrystalConfig.mcellbbox.zMin=mcellbbox.zMin;
   sGLCrystalConfig.mcellbbox.zMax=mcellbbox.zMax;
    
   sGLCrystalConfig.mmapbbox.xMin=mmapbbox.xMin;
   sGLCrystalConfig.mmapbbox.xMax=mmapbbox.xMax;
   sGLCrystalConfig.mmapbbox.yMin=mmapbbox.yMin;
   sGLCrystalConfig.mmapbbox.yMax=mmapbbox.yMax;
   sGLCrystalConfig.mmapbbox.zMin=mmapbbox.zMin;
   sGLCrystalConfig.mmapbbox.zMax=mmapbbox.zMax;
   
   sGLCrystalConfig.mViewCntr.x=mViewCntr.x;
   sGLCrystalConfig.mViewCntr.y=mViewCntr.y;
   sGLCrystalConfig.mViewCntr.z=mViewCntr.z;
   sGLCrystalConfig.mSaved=true;
}

void WXGLCrystalCanvas::OnExit(wxCommandEvent &event)
{
   
}

void WXGLCrystalCanvas::OnPaint(wxPaintEvent &event)
{
   VFN_DEBUG_ENTRY("WXGLCrystalCanvas::OnPaint()",7)
   wxPaintDC dc(this);
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
      const GLfloat colour2 [] = {1.00, 1.00, 1.00, 1.00}; 
      glMaterialfv(GL_FRONT, GL_EMISSION,  colour2); 
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
      statusText.Printf(_T("Center@(%5.3f,%5.3f,%5.3f)"),x,y,z);
      for(unsigned int i=0;i<mvpUnitCellMap.size();++i)
      {
         statusText+=_T(", map(") + wxString::FromAscii(mvpUnitCellMap[i]->GetName().c_str())
                     +wxString::Format(_T(")=%5.2fe"),mvpUnitCellMap[i]->GetValue(x,y,z));
      }
      mpParentFrame->SetStatusText(statusText);
   }
   if(mShowFourier)
   {
      glLoadIdentity();
      glTranslatef( 0, 0, -mDist );
      build_rotmatrix( m,mQuat);
      glMultMatrixf( &m[0][0] );
      glTranslatef( mX0, mY0, mZ0 );
      glPushMatrix();
         // The display origin is the center of the Crystal BoundingBox, so translate
            BBox cellbbox = this->GetCellBBox();
            REAL xc=(cellbbox.xMin+cellbbox.xMax)/2.;
            REAL yc=(cellbbox.yMin+cellbbox.yMax)/2.;
            REAL zc=(cellbbox.zMin+cellbbox.zMax)/2.;
            mpWXCrystal->GetCrystal().FractionalToOrthonormalCoords(xc, yc, zc);
            glTranslatef(-xc, -yc, -zc);
         // Draw all Fourier maps
         vector<boost::shared_ptr<UnitCellMapGLList> >::const_iterator pos;
         for(pos=mvpUnitCellMapGLList.begin();pos != mvpUnitCellMapGLList.end();++pos)
            if((*pos)->Show()) (*pos)->Draw();
     glPopMatrix();
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

   this->SetCurrent();
   glViewport(0, 0, width, height);
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnSize():"<<mViewAngle<<","<<width<<","<<height<<","<<mDist,2)
   if( (width>0)&&(height>0)) //in case the window is docked...
      gluPerspective(mViewAngle,(float)width/(float)height,1.f,2.*mDist);
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
      case(WXK_PAGEUP): mZ0 -= 0.1; Refresh(FALSE); break;
      case(WXK_PAGEDOWN): mZ0 += 0.1; Refresh(FALSE); break;
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
   event.Skip();
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
   VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnMouse()"
                     <<endl<<"IsButton():"<<event.IsButton()
                     <<endl<<"ButtonDown():"<<event.ButtonDown()
                     <<endl<<"Dragging():"<<event.Dragging()
                     <<endl<<"Entering():"<<event.Entering()
                     <<endl<<"Leaving():"<<event.Leaving()
                     <<endl<<"GetButton()"<<event.GetButton()
                     <<endl<<"GetWheelAxis():"<<event.GetWheelAxis()
                     <<endl<<"GetWheelDelta():"<<event.GetWheelDelta()
                     <<endl<<"GetWheelRotation():"<<event.GetWheelRotation()
                     <<endl<<"Moving():"<<event.Moving()
                     <<endl,7)
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
   //else if(event.Leaving()) cout<<"Mouse is leaving window!!"<<endl;
   else if(event.RightIsDown())
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
   else if (event.GetWheelDelta()>0)
   {// Double-touch event on OSX + trackpad
      if(event.ControlDown())
      {
         VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnMouse(): Mouse Wheel / double touch + control (OSX: command)",2)
         //Change zoom / angle
         int width, height;
         GetClientSize(& width, & height);
         const int delta=event.GetWheelDelta();
         const int rotation=event.GetWheelRotation();
         
         if(event.GetWheelAxis()==0) mDist *= (1.+float(rotation)/100.);
         else
         {
            mDist /= (1.+float(rotation)/100.);
            mViewAngle *=(1.+float(rotation)/100.);
         }
         SetCurrent();
         glMatrixMode(GL_PROJECTION);
         glLoadIdentity();
         if( (width>0)&&(height>0)) //in case size is null...
            gluPerspective(mViewAngle,(float)width/(float)height,
                           (mDist>101)?(mDist-100):1.,mDist+100);
         Refresh(FALSE);
         VFN_DEBUG_MESSAGE(mViewAngle <<" "<<mDist,2)
      }
      else
      {
         VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnMouse(): Mouse Wheel / double touch",2)
         // Rotate view with trackball
         int width, height;
         GetClientSize(& width, & height);
         const int delta=event.GetWheelDelta();
         int dx=0,dy=0;
         if(event.GetWheelAxis()==1) dx=-event.GetWheelRotation();
         else dy=event.GetWheelRotation();
         float spin_quat[4];
         trackball(spin_quat,
                   (2.0*mTrackBallLastX -       width) / (width+.001),  //normalizing from -1 to 1
                   (     height - 2.0*mTrackBallLastY) / (height+.001),
                   (2.0*(mTrackBallLastX+dx) -       width) / (width+.001),
                   (     height - 2.0*(mTrackBallLastY+dy)) / (height+.001));
         
         add_quats( spin_quat, mQuat, mQuat );
         Refresh(FALSE);
      }
   }

   mTrackBallLastX = event.GetX();
   mTrackBallLastY = event.GetY();
   event.Skip();
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
   // This can only be called from the main (graphical) thread
   VFN_DEBUG_ENTRY("WXGLCrystalCanvas::CrystUpdate():"<<wxThread::IsMain(),10)
   // Update the list of available maps & update the displayed ones
   if(mpFourierMapListWin!=0) mpFourierMapListWin->mMutex.Lock();
   // Remove maps that cannot be computed any more
   for(vector<boost::shared_ptr<UnitCellMap> >::iterator pos=mvpUnitCellMap.begin();pos!=mvpUnitCellMap.end();)
   {
      bool keep=false;
      if((*pos)->GetType()==3) keep=true;
      if(((*pos)->GetType()==0)||((*pos)->GetType()==2))
      {
         /*
         cout<<"WXGLCrystalCanvas::CrystUpdate()"<<endl<<(*pos)->GetName()<<(*pos)->GetType()<<" "
            <<mpWXCrystal->GetCrystal().GetScatteringComponentList().GetNbComponent()<<","
            <<(*pos)->GetData()->GetFhklObsSq().numElements()<<endl;
         */
         if(mpWXCrystal->GetCrystal().GetClientRegistry().Find((RefinableObj*)(*pos)->GetData())>=0)
            if(mpWXCrystal->GetCrystal().GetScatteringComponentList().GetNbComponent()>0)
               if((*pos)->GetData()->GetFhklObsSq().numElements()>0)
                  keep=true;
      }
      if((*pos)->GetType()==1)
      {
         if(mpWXCrystal->GetCrystal().GetClientRegistry().Find((RefinableObj*)(*pos)->GetData())>=0)
            if(mpWXCrystal->GetCrystal().GetScatteringComponentList().GetNbComponent()>0)
               keep=true;
      }
      VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::CrystUpdate()"<<(*pos)->GetName()<<(*pos)->GetType()<<": keep="<<keep,8)
      if(!keep)
      {
         //erase corresponding gl maps
         for(vector<boost::shared_ptr<UnitCellMapGLList> >::iterator 
               posgl=mvpUnitCellMapGLList.begin();posgl!=mvpUnitCellMapGLList.end();)
         {
            if(&(**pos)==&((*posgl)->GetMap()))
            {
               VFN_DEBUG_MESSAGE("Erasing GL map:"<<(*posgl)->GetName(),8)
               posgl=mvpUnitCellMapGLList.erase(posgl);
            } else ++posgl;
         }
         pos=mvpUnitCellMap.erase(pos);
      }
      else 
      {
         if((*pos)->GetType()!=3)
         {
            #ifdef HAVE_FFTW
            // During optimization, only update Fourier maps if one is displayed or the Fourier win is opened
            if(  (mvpUnitCellMapGLList.size()>0)
               ||(!(mpWXCrystal->GetCrystal().IsBeingRefined()))
               ||(mpFourierMapListWin!=0))
               (*pos)->CalcFourierMap(*((*pos)->GetData()),(*pos)->GetType(),mSharpenMap);
            #endif
         }
         ++pos;
      }
   }
   #ifdef HAVE_FFTW
   // Add newly computable maps
   if(mpWXCrystal->GetCrystal().GetScatteringComponentList().GetNbComponent()>0)
      for(int i=0;i<mpWXCrystal->GetCrystal().GetClientRegistry().GetNb();++i)
      {
         ScatteringData* data=dynamic_cast<ScatteringData *>(&(mpWXCrystal->GetCrystal().GetClientRegistry().GetObj(i)));
         
         if(data!=0)
         {
            // Add if not already listed
            bool addCalcMap=true,addObsDiffMaps=true;
            for(vector<boost::shared_ptr<UnitCellMap> >::iterator pos=mvpUnitCellMap.begin();pos!=mvpUnitCellMap.end();++pos)
               if((*pos)->GetData()==data)
               {
                  if((*pos)->GetType()==1) addCalcMap=false;
                  if((*pos)->GetType()==0) addObsDiffMaps=false;//type==2 will also be there
               }
            //cout<<__FILE__<<":"<<__LINE__<<":WXGLCrystalCanvas::CrystUpdate()"
            //    <<data<<","<<addCalcMap<<","<<addObsDiffMaps<<endl;
            if(addCalcMap&&(data->GetNbRefl()>0))
            {
               mvpUnitCellMap.push_back(boost::shared_ptr<UnitCellMap>(new UnitCellMap(mpWXCrystal->GetCrystal())));
               mvpUnitCellMap.back()->CalcFourierMap(*data,1);
               VFN_DEBUG_MESSAGE("Added GL map:"<<mvpUnitCellMap.back()->GetName(),8)
            }
            if(addObsDiffMaps && (data->GetFhklObsSq().numElements()>0) )
            {
               mvpUnitCellMap.push_back(boost::shared_ptr<UnitCellMap>(new UnitCellMap(mpWXCrystal->GetCrystal())));
               mvpUnitCellMap.back()->CalcFourierMap(*data,0);
               VFN_DEBUG_MESSAGE("Added GL map:"<<mvpUnitCellMap.back()->GetName()<<":"<<data->GetFhklObsSq().numElements(),8)
               mvpUnitCellMap.push_back(boost::shared_ptr<UnitCellMap>(new UnitCellMap(mpWXCrystal->GetCrystal())));
               mvpUnitCellMap.back()->CalcFourierMap(*data,2);
               VFN_DEBUG_MESSAGE("Added GL map:"<<mvpUnitCellMap.back()->GetName()<<":"<<data->GetFhklObsSq().numElements(),8)
            }
         }
      }
   #endif
   //update GL maps
   for(vector<boost::shared_ptr<UnitCellMapGLList> >::iterator 
         pos=mvpUnitCellMapGLList.begin();pos!=mvpUnitCellMapGLList.end();++pos)
   {
      //cout<<"Updating GL map:"<<(*pos)->GetName()<<endl;
      (*pos)->GenList();
   }
   if(mpFourierMapListWin!=0) mpFourierMapListWin->mMutex.Unlock();
   VFN_DEBUG_EXIT("WXGLCrystalCanvas::CrystUpdate()",10)

   wxUpdateUIEvent event(ID_GLCRYSTAL_UPDATEUI);
   wxPostEvent(this,event);
   /* // To make a movie
   if(mpWXCrystal->GetCrystal().IsBeingRefined())
   {
      // Export POV-Ray file to make a movie
      char povFile[40];
      time_t date=time(0);
      strftime(povFile,sizeof(povFile),"pov/%Y%m%d-%Hh%Mm%Ss%Z.pov",gmtime(&date));//%Y-%m-%dT%H:%M:%S%Z
      this->POVRayOutput(povFile);
   }
   */
}

void WXGLCrystalCanvas::OnUpdateUI(wxUpdateUIEvent&event)
{
   VFN_DEBUG_ENTRY("WXGLCrystalCanvas::OnUpdateUI()",5)
   if(mpFourierMapListWin!=0)
   {
      mpFourierMapListWin->mIsUpdating=true;
      mpFourierMapListWin->mMutex.Lock();
      
      wxArrayString maps;
      for(vector<boost::shared_ptr<UnitCellMap> >::iterator 
          pos=mvpUnitCellMap.begin();pos!=mvpUnitCellMap.end();++pos)
            maps.Add( wxString::FromAscii((*pos)->GetName().c_str()));
      if(mpFourierMapListWin->mpAvailableMapList->GetStrings()!=maps)
         mpFourierMapListWin->mpAvailableMapList->Set(maps);
      
      wxArrayString glmaps;
      for(vector<boost::shared_ptr<UnitCellMapGLList> >::iterator 
          pos=mvpUnitCellMapGLList.begin();pos!=mvpUnitCellMapGLList.end();++pos)
            glmaps.Add( wxString::FromAscii((*pos)->GetName().c_str()));
      if(mpFourierMapListWin->mpDisplayedMapList->GetStrings()!=glmaps)
         mpFourierMapListWin->mpDisplayedMapList->Set(glmaps);
      
      if(mpFourierMapListWin->mpAvailableMapList->GetSelection()>=0)
      {
         boost::shared_ptr<ObjCryst::UnitCellMap> pMap=mvpUnitCellMap[mpFourierMapListWin->mpAvailableMapList->GetSelection()];
         mpFourierMapListWin->mpMapInfo->SetLabel(wxString::Format(_T("min=%5.2f max=%5.2f sigma=%5.2f"),
                                                pMap->Min(),pMap->Max(),pMap->StandardDeviation()));
      }
      mpFourierMapListWin->mMutex.Unlock();
      mpFourierMapListWin->mIsUpdating=false;
   }
   this->Refresh(false);
   event.Skip();
   VFN_DEBUG_EXIT("WXGLCrystalCanvas::OnUpdateUI()",5)
}

void WXGLCrystalCanvas::SetCurrent()
{
   VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::SetCurrent()",4)
   this->wxGLCanvas::SetCurrent(*mpwxGLContext);
   #ifndef HAVE_GLUT
   this->BuildGLFont();
   sFontDisplayListBase=mGLFontDisplayListBase;
   #endif
}

void WXGLCrystalCanvas::NotifyDeleteFourierWin()
{
   mpFourierMapListWin=0;
}

void WXGLCrystalCanvas::InitGL()
{
   VFN_DEBUG_ENTRY("WXGLCrystalCanvas::InitGL()",8)
   this->SetCurrent();
   #ifdef HAVE_GLUT
   static bool needglutinit=true;
   if(needglutinit)
   {
      needglutinit=false;
      //glutInit(&(wxApp::GetInstance()->argc),wxApp::GetInstance()->argv);
      char **argv=new char*;
      int argc=0;
      glutInit(&argc,argv);// We cannot pass arguments directly in Unicode mode, so...
   }
   #endif
   
   int width, height;
   GetClientSize(& width, & height);
   glViewport(0, 0, width, height);
      
   glEnable(GL_DEPTH_TEST);
   glEnable(GL_LIGHTING);
   //glEnable (GL_BLEND);
   //glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
   
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
   
   if(sGLCrystalConfig.mSaved==false)
   {
      //Initialize Trackball
      trackball(mQuat,0.,0.,0.,0.);
   }
   wxSizeEvent ev;
   wxPostEvent(this,ev);
   
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
      vector<boost::shared_ptr<UnitCellMapGLList> >::iterator pos;
      for(pos=mvpUnitCellMapGLList.begin();pos != mvpUnitCellMapGLList.end();pos++)
      {
         wxBusyInfo wait(_T("Processing Fourier Map..."));
         (*pos)->GenList();
      }
   }
   if(event.GetId()==ID_GLCRYSTAL_MENU_LIMITS_ASYMCELL)
   {
      mcellbbox.xMin = mpWXCrystal->GetCrystal().GetSpaceGroup().GetAsymUnit().Xmin()-0.1;
      mcellbbox.yMin = mpWXCrystal->GetCrystal().GetSpaceGroup().GetAsymUnit().Ymin()-0.1;
      mcellbbox.zMin = mpWXCrystal->GetCrystal().GetSpaceGroup().GetAsymUnit().Zmin()-0.1;
      mcellbbox.xMax = mpWXCrystal->GetCrystal().GetSpaceGroup().GetAsymUnit().Xmax()+0.1;
      mcellbbox.yMax = mpWXCrystal->GetCrystal().GetSpaceGroup().GetAsymUnit().Ymax()+0.1;
      mcellbbox.zMax = mpWXCrystal->GetCrystal().GetSpaceGroup().GetAsymUnit().Zmax()+0.1;
      vector<boost::shared_ptr<UnitCellMapGLList> >::iterator pos;
      for(pos=mvpUnitCellMapGLList.begin();pos != mvpUnitCellMapGLList.end();pos++)
      {
         wxBusyInfo wait(_T("Processing Fourier Map..."));
         (*pos)->GenList();
      }
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
         vector<boost::shared_ptr<UnitCellMapGLList> >::iterator pos;
         for(pos=mvpUnitCellMapGLList.begin();pos != mvpUnitCellMapGLList.end();pos++)
         {
            wxBusyInfo wait(_T("Processing Fourier Map..."));
            (*pos)->GenList();
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
   if(mShowCrystal) mpPopUpMenu->SetLabel(ID_GLCRYSTAL_MENU_SHOWCRYSTAL, _T("Show Crystal"));
   else mpPopUpMenu->SetLabel(ID_GLCRYSTAL_MENU_SHOWCRYSTAL, _T("Hide Crystal"));
   mShowCrystal = !mShowCrystal;
   if(!(mpWXCrystal->GetCrystal().IsBeingRefined())) this->CrystUpdate();
}

void WXGLCrystalCanvas::OnShowAtomLabel( wxCommandEvent & WXUNUSED(event))
{
   if(mShowAtomName) mpPopUpMenu->SetLabel(ID_GLCRYSTAL_MENU_SHOWATOMLABEL, _T("Show Atom Labels"));
   else mpPopUpMenu->SetLabel(ID_GLCRYSTAL_MENU_SHOWATOMLABEL, _T("Hide Atom Labels"));
   mShowAtomName= !mShowAtomName;
   if(!(mpWXCrystal->GetCrystal().IsBeingRefined())) this->CrystUpdate();
}

void WXGLCrystalCanvas::OnShowCursor( wxCommandEvent & WXUNUSED(event))
{
   if(mShowCursor) mpPopUpMenu->SetLabel(ID_GLCRYSTAL_MENU_SHOWCURSOR, _T("Show Cursor"));
   else mpPopUpMenu->SetLabel(ID_GLCRYSTAL_MENU_SHOWCURSOR, _T("Hide Cursor"));
   mShowCursor= !mShowCursor;
   if(!(mpWXCrystal->GetCrystal().IsBeingRefined())) this->CrystUpdate();
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
void WXGLCrystalCanvas::OnFourier(wxCommandEvent &event)
{
   if(event.GetId()==ID_GLCRYSTAL_MENU_FOURIER)
   {
      if(mpFourierMapListWin!=0) return;
      if(mpWXCrystal->GetCrystal().IsBeingRefined())
      {
         wxMessageBox(_T("The Fourier maps dialog \ncannot be opened during an optimization"), _T("Error"), wxOK, this);
         return;
      }
      wxFrame *frame= new wxMiniFrame(this,-1, wxString::FromAscii(("Available Fourier maps for "+mpWXCrystal->GetCrystal().GetName()).c_str()),
                                          wxDefaultPosition,wxSize(500,500),wxCLOSE_BOX|wxCAPTION|wxSYSTEM_MENU);
      mpFourierMapListWin=new WXFourierMapList(this,frame);
      mpFourierMapListWin->mpWireFrame->SetValue(true);
      mpFourierMapListWin->mpShowFourier->SetValue(mShowFourier);
      mpFourierMapListWin->mpSharpenMap->SetValue(mSharpenMap);
      frame->Show(true);
      mpWXCrystal->GetCrystal().UpdateDisplay();
      return;
   }
   if(mpFourierMapListWin==0) return;
   if(mpFourierMapListWin->mIsUpdating) return;
   if(  (event.GetId()==ID_GLCRYSTAL_FOURIER_UPDATE)
      ||(event.GetId()==ID_GLCRYSTAL_FOURIER_CONTOUR))
   {
      mpFourierMapListWin->mMutex.Lock();
      //Changed colour or contour ?
      unsigned int choice=mpFourierMapListWin->mpDisplayedMapList->GetSelection();
      if(wxNOT_FOUND!=choice)
      {
         double contour;
         mpFourierMapListWin->mpContourValue->GetValue().ToDouble(&contour);
         wxColour col(mpFourierMapListWin->mpColourPicker->GetColour());
         if(abs((float)contour-mvpUnitCellMapGLList[choice]->GetContour())>.0001)
         {
            mvpUnitCellMapGLList[choice]->SetContour((float)contour);
            if(false==mpWXCrystal->GetCrystal().IsBeingRefined())
            {
               wxBusyInfo wait(_T("Processing Fourier Map..."));
               mvpUnitCellMapGLList[choice]->GenList();
            }
         }
         mvpUnitCellMapGLList[choice]->SetColour(col.Red()/255.0,col.Green()/255.0,col.Blue()/255.0,0.5);
      }
      mpFourierMapListWin->mMutex.Unlock();
   }
   
   if(event.GetId()==ID_GLCRYSTAL_FOURIER_LISTMAP)
   {// Selected one map
      mpFourierMapListWin->mMutex.Lock();
      if(mpFourierMapListWin->mpAvailableMapList->GetSelection()>=0)
      {
         boost::shared_ptr<ObjCryst::UnitCellMap> pMap=mvpUnitCellMap[mpFourierMapListWin->mpAvailableMapList->GetSelection()];
         mpFourierMapListWin->mpMapInfo->SetLabel(wxString::Format(_T("min=%5.2f max=%5.2f sigma=%5.2f"),
                                                pMap->Min(),pMap->Max(),pMap->StandardDeviation()));
      }
	  mpFourierMapListWin->mMutex.Unlock();
   }
   if(event.GetId()==ID_GLCRYSTAL_FOURIER_LISTGLMAP)
   {
      mpFourierMapListWin->mMutex.Lock();
      if(mpFourierMapListWin->mpDisplayedMapList->GetSelection()>=0)
      {
         boost::shared_ptr<UnitCellMapGLList> pMap=mvpUnitCellMapGLList[mpFourierMapListWin->mpDisplayedMapList->GetSelection()];
         mpFourierMapListWin->mpContourValue->SetValue(wxString::Format(_T("%5.2f"),pMap->GetContour()));
         mpFourierMapListWin->mpColourPicker->SetColour(wxColour(pMap->GetColour()[0]*255,pMap->GetColour()[1]*255,
                                                                 pMap->GetColour()[2]*255,pMap->GetColour()[3]*255));
      }
	  mpFourierMapListWin->mMutex.Unlock();
   }
   if((event.GetId()==ID_GLCRYSTAL_FOURIER_ADD)||(event.GetId()==ID_GLCRYSTAL_FOURIER_NEWCONTOUR))
   {
      mpFourierMapListWin->mMutex.Lock();
      if(mpFourierMapListWin->mpAvailableMapList->GetSelection()!=wxNOT_FOUND)
      {
         boost::shared_ptr<ObjCryst::UnitCellMap> pMap=mvpUnitCellMap[mpFourierMapListWin->mpAvailableMapList->GetSelection()];
         double contour=0;
         wxString scontour=mpFourierMapListWin->mpNewContourValue->GetValue();
         if(scontour==_T("")) contour=pMap->Min()+pMap->StandardDeviation()*3;
         else scontour.ToDouble(&contour);
         wxColor ncolor(255,0,0);
         ncolor = wxGetColourFromUser((wxWindow*)this, ncolor);
   
         wxBusyInfo wait(_T("Processing Fourier Map..."));
         mvpUnitCellMapGLList.push_back(boost::shared_ptr<UnitCellMapGLList>(new UnitCellMapGLList(*pMap,this,true,(float)contour)));
         mvpUnitCellMapGLList.back()->SetName(pMap->GetName());
         mvpUnitCellMapGLList.back()->SetColour(ncolor.Red()/255.0,ncolor.Green()/255.0,ncolor.Blue()/255.0,0.5);
         mpFourierMapListWin->mMutex.Unlock();
         if(false==mpWXCrystal->GetCrystal().IsBeingRefined())
         {
            this->SetCurrent();
            mvpUnitCellMapGLList.back()->GenList();
         }
      }
   }
   if(event.GetId()==ID_GLCRYSTAL_FOURIER_REMOVE)
   {
      mpFourierMapListWin->mMutex.Lock();
      unsigned int choice=mpFourierMapListWin->mpDisplayedMapList->GetSelection();
      if(wxNOT_FOUND!=choice)
         mvpUnitCellMapGLList.erase(mvpUnitCellMapGLList.begin()+choice);
      mpFourierMapListWin->mMutex.Unlock();
   }
   if(event.GetId()==ID_GLCRYSTAL_FOURIER_SHOW)
   {
      mpFourierMapListWin->mMutex.Lock();
      mShowFourier=mpFourierMapListWin->mpShowFourier->GetValue();
      mpFourierMapListWin->mMutex.Unlock();
   }
   if(event.GetId()==ID_GLCRYSTAL_FOURIER_WIREFRAME)
   {
      mpFourierMapListWin->mMutex.Lock();
      vector<boost::shared_ptr<UnitCellMapGLList> >::iterator pos;
      for(pos=mvpUnitCellMapGLList.begin();pos != mvpUnitCellMapGLList.end();pos++)
         (*pos)->ToggleShowWire();
      mpFourierMapListWin->mMutex.Unlock();
   }
   
   if(event.GetId()==ID_GLCRYSTAL_FOURIER_SHARPEN)
   {
      mpFourierMapListWin->mMutex.Lock();
      mSharpenMap=mpFourierMapListWin->mpSharpenMap->GetValue();
      mpFourierMapListWin->mMutex.Unlock();
   }
   
   // Update - if the crystal is being refined, it will be done at the next display update
   if(false==mpWXCrystal->GetCrystal().IsBeingRefined())
      this->CrystUpdate();
}

void WXGLCrystalCanvas::OnLoadFourierGRD( wxCommandEvent & WXUNUSED(event))
{
   wxFileDialog fd((wxWindow*)this, _T("Choose a file containing a Fourier Map"),
           _T(""), _T(""), _T("Fourier Map files (*.grd)|*.grd"), wxFD_OPEN | wxFD_FILE_MUST_EXIST);
   //if okay then read Fourier map, run MC on it and display the triangles
   if(fd.ShowModal() == wxID_OK)
   {
      const string filename(fd.GetPath().ToAscii());
      UnitCellMap *pMap=new UnitCellMap(mpWXCrystal->GetCrystal());
      if (pMap->ImportGRD(filename) == 0)
      {
         string tmp="Error reading Fourier file:"+filename;
         wxMessageBox( wxString::FromAscii(tmp.c_str()), _T("File error"), wxOK, this);
      return;
      }
      this->AddFourier(pMap);
   }
}

void WXGLCrystalCanvas::OnLoadFourierDSN6( wxCommandEvent & WXUNUSED(event))
{
   wxFileDialog fd((wxWindow*)this, _T("Choose a file containing a Fourier Map"),
           _T(""), _T(""), _T("Fourier Map files (*.DN6)|*.DN6"), wxFD_OPEN | wxFD_FILE_MUST_EXIST);
   //if okay then read Fourier map, run MC on it and display the triangles
   if(fd.ShowModal() == wxID_OK)
   {
      const string filename(fd.GetPath().ToAscii());
      UnitCellMap *pMap=new UnitCellMap(mpWXCrystal->GetCrystal());
      if (pMap->ImportDSN6(filename) == 0)
      {
         string tmp="Error reading Fourier file:"+filename;
         wxMessageBox( wxString::FromAscii(tmp.c_str()), _T("File error"), wxOK, this);
      return;
      }
      this->AddFourier(pMap);
   }
}

void WXGLCrystalCanvas::AddFourier(UnitCellMap *map)
{
   mvpUnitCellMap.push_back(boost::shared_ptr<UnitCellMap>(map));
   wxBusyInfo wait(_T("Processing Fourier Map..."));
   {
      float contour=map->Mean()+2*map->StandardDeviation();
      if(contour>map->Max()) contour=map->Mean()+0.75*(map->Max()-map->Mean());
      mvpUnitCellMapGLList.push_back(boost::shared_ptr<UnitCellMapGLList>(new UnitCellMapGLList(*map,this)));
      switch(mvpUnitCellMapGLList.size())
      {
         case 1: mvpUnitCellMapGLList.back()->SetColour(1.,0.,0.,.5);break;
         case 2: mvpUnitCellMapGLList.back()->SetColour(0.,0.,1.,.5);break;
         default:mvpUnitCellMapGLList.back()->SetColour(0.,1.,0.,.5);break;
      }
      this->SetCurrent();
      mvpUnitCellMapGLList.back()->GenList();
      mvpUnitCellMapGLList.back()->SetName(map->GetName());
   }
   if(!(mpWXCrystal->GetCrystal().IsBeingRefined())) this->CrystUpdate();
}


void WXGLCrystalCanvas::OnFourierChangeColour(wxColourPickerEvent  &event)
{
   mpFourierMapListWin->mMutex.Lock();
   //Changed colour or contour ?
   unsigned int choice=mpFourierMapListWin->mpDisplayedMapList->GetSelection();
   if(wxNOT_FOUND!=choice)
   {
      double contour;
      mpFourierMapListWin->mpContourValue->GetValue().ToDouble(&contour);
      wxColour col(mpFourierMapListWin->mpColourPicker->GetColour());
      if(abs((float)contour-mvpUnitCellMapGLList[choice]->GetContour())>.0001)
      {
         wxBusyInfo wait(_T("Processing Fourier Map..."));
         mvpUnitCellMapGLList[choice]->SetContour((float)contour);
         mvpUnitCellMapGLList[choice]->GenList();
      }
      mvpUnitCellMapGLList[choice]->SetColour(col.Red()/255.0,col.Green()/255.0,col.Blue()/255.0,0.5);
   }
   mpFourierMapListWin->mMutex.Unlock();
   mpWXCrystal->GetCrystal().UpdateDisplay();
}

/*
void WXGLCrystalCanvas::OnUnloadFourier( wxCommandEvent & WXUNUSED(event))
{
   wxMessageDialog * msure = new wxMessageDialog((wxWindow*)this,
     "Are you sure you want to unload all Fourier Map Data?", "Unload Fourier Map", wxYES_NO | wxNO_DEFAULT |
     wxICON_QUESTION );
   if(msure->ShowModal() == wxID_YES)
   {
      mvpUnitCellMap.clear();
      mvpUnitCellMapGLList.clear();
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
*/
void WXGLCrystalCanvas::OnPOVRay( wxCommandEvent & WXUNUSED(event))
{
   WXCrystValidateAllUserInput();
   wxFileDialog save(this,_T("Choose filename"),_T(""),_T(""),_T("*.pov"),wxFD_SAVE | wxFD_OVERWRITE_PROMPT);
   if(save.ShowModal() != wxID_OK) return;
   this->POVRayOutput(string(save.GetPath().char_str()));
}

void WXGLCrystalCanvas::POVRayOutput(const std::string &filename)
{
   ofstream os(filename.c_str());
   
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
      wxBusyInfo wait(_T("Processing Fourier Map..."));
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
      vector<boost::shared_ptr<UnitCellMapGLList> >::const_iterator pos;
      for(pos=mvpUnitCellMapGLList.begin();pos != mvpUnitCellMapGLList.end();++pos)
      {
         const float *prgbf=(*pos)->GetColour();
         if((*pos)->ShowWire())
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
            (*pos)->GetMap().POVRayDescription(os,(*pos)->GetContour(),options);
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
            (*pos)->GetMap().POVRayDescription(os,(*pos)->GetContour(),options);
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
   os.close();
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
                        _T("Helvetica"));          // Font name

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
  : wxDialog((wxWindow *)parent, -1, _T("Set bounding box"), wxDefaultPosition,
  	     wxSize(250, 250), wxDEFAULT_DIALOG_STYLE) 
{
  wxBoxSizer *dialogSizer = new wxBoxSizer(wxVERTICAL);
  wxFlexGridSizer *inputSizer = new wxFlexGridSizer(4, 3, 10, 10);
  // headers
  inputSizer->Add(new wxStaticText(this, -1, _T("")), 0, wxALIGN_CENTRE_VERTICAL);
  inputSizer->Add(new wxStaticText(this, -1, _T("minimum")), 0, wxALIGN_CENTER);
  inputSizer->Add(new wxStaticText(this, -1, _T("maximum")), 0, wxALIGN_CENTER);
  // 1st row
  inputSizer->Add(new wxStaticText(this, -1, _T("a")), 0, wxALIGN_CENTRE_VERTICAL);
  inputSizer->Add(mpXminCtrl = new wxTextCtrl(this, -1, 
					      wxString::Format(_T("%f"),bbox.xMin)), 
					      0, wxALIGN_CENTRE_VERTICAL);
  inputSizer->Add(mpXmaxCtrl = new wxTextCtrl(this, -1, 
					      wxString::Format(_T("%f"),bbox.xMax)), 
					      0, wxALIGN_CENTRE_VERTICAL);
  // 2nd row
  inputSizer->Add(new wxStaticText(this, -1, _T("b")), 0, wxALIGN_CENTRE_VERTICAL);
  inputSizer->Add(mpYminCtrl = new wxTextCtrl(this, -1, 
					      wxString::Format(_T("%f"),bbox.yMin)), 
					      0, wxALIGN_CENTRE_VERTICAL);
  inputSizer->Add(mpYmaxCtrl = new wxTextCtrl(this, -1, 
					      wxString::Format(_T("%f"),bbox.yMax)), 
					      0, wxALIGN_CENTRE_VERTICAL);
  // 3rd row
  inputSizer->Add(new wxStaticText(this, -1, _T("c")), 0, wxALIGN_CENTRE_VERTICAL);
  inputSizer->Add(mpZminCtrl = new wxTextCtrl(this, -1, 
					      wxString::Format(_T("%f"),bbox.zMin)), 
					      0, wxALIGN_CENTRE_VERTICAL);
  inputSizer->Add(mpZmaxCtrl = new wxTextCtrl(this, -1, 
					      wxString::Format(_T("%f"),bbox.zMax)), 
					      0, wxALIGN_CENTRE_VERTICAL);
  // button section
  wxFlexGridSizer *buttonSizer = new wxFlexGridSizer(1, 2, 10, 10);
  buttonSizer->Add(new wxButton(this, wxID_OK, _T("OK")), 
		   0, wxALIGN_CENTRE_VERTICAL);
  buttonSizer->Add(new wxButton(this, wxID_CANCEL, _T("Cancel")), 
		   0, wxALIGN_CENTRE_VERTICAL);

  dialogSizer->Add(10, 10);
  dialogSizer->Add(new wxStaticText(this, -1,  wxString::FromAscii(title)), 0, 
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
  double val;

  mpXminCtrl->GetValue().ToDouble(&val);
  mbbox.xMin = val;
  mpXmaxCtrl->GetValue().ToDouble(&val);
  mbbox.xMax = val;
  if (mbbox.xMin == mbbox.xMax) {wxMessageBox(_T("Sorry, Xmin must be less than Xmax!"), _T("Zero bounding volume"), wxOK, this); return;}
  if (mbbox.xMin > mbbox.xMax) {
    float tmp = mbbox.xMax;
    mbbox.xMax = mbbox.xMin;
    mbbox.xMin = tmp;
  }
  VFN_DEBUG_MESSAGE("Xmin " << mbbox.xMin << " Xmax " << mbbox.xMax,1)

  mpYminCtrl->GetValue().ToDouble(&val);
  mbbox.yMin = val;
  mpYmaxCtrl->GetValue().ToDouble(&val);
  mbbox.yMax = val;
  if (mbbox.yMin == mbbox.yMax) {wxMessageBox(_T("Sorry, Ymin must be less than Ymax!"), _T("Zero bounding volume"), wxOK, this); return;}
  if (mbbox.yMin > mbbox.yMax) {
    float tmp = mbbox.yMax;
    mbbox.yMax = mbbox.yMin;
    mbbox.yMin = tmp;
  }
  VFN_DEBUG_MESSAGE("Ymin " << mbbox.yMin << " Ymax " << mbbox.yMax,1)

  mpZminCtrl->GetValue().ToDouble(&val);
  mbbox.zMin = val;
  mpZmaxCtrl->GetValue().ToDouble(&val);
  mbbox.zMax = val;
  if (mbbox.zMin == mbbox.zMax) {wxMessageBox(_T("Sorry, Zmin must be less than Zmax!"), _T("Zero bounding volume"), wxOK, this); return;}
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
  : wxDialog((wxWindow *)parent, -1, _T("Set position"), wxDefaultPosition,
  	     wxSize(250, 250), wxDEFAULT_DIALOG_STYLE) 
{
  wxBoxSizer *dialogSizer = new wxBoxSizer(wxVERTICAL);
  wxFlexGridSizer *inputSizer = new wxFlexGridSizer(3, 2, 10, 10);
  // 1st row
  inputSizer->Add(new wxStaticText(this, -1, _T("x")), 0, wxALIGN_CENTRE_VERTICAL);
  inputSizer->Add(mpXCtrl = new wxTextCtrl(this, -1, 
					   wxString::Format(_T("%.3f"),xyz.x)), 
					   0, wxALIGN_CENTRE_VERTICAL);
  // 2nd row
  inputSizer->Add(new wxStaticText(this, -1, _T("y")), 0, wxALIGN_CENTRE_VERTICAL);
  inputSizer->Add(mpYCtrl = new wxTextCtrl(this, -1, 
					   wxString::Format(_T("%.3f"),xyz.y)), 
					   0, wxALIGN_CENTRE_VERTICAL);
  // 3rd row
  inputSizer->Add(new wxStaticText(this, -1, _T("z")), 0, wxALIGN_CENTRE_VERTICAL);
  inputSizer->Add(mpZCtrl = new wxTextCtrl(this, -1, 
					   wxString::Format(_T("%.3f"),xyz.z)), 
					   0, wxALIGN_CENTRE_VERTICAL);
  // button section
  wxFlexGridSizer *buttonSizer = new wxFlexGridSizer(1, 2, 10, 10);
  buttonSizer->Add(new wxButton(this, wxID_OK, _T("OK")), 
		   0, wxALIGN_CENTRE_VERTICAL);
  buttonSizer->Add(new wxButton(this, wxID_CANCEL, _T("Cancel")), 
		   0, wxALIGN_CENTRE_VERTICAL);

  dialogSizer->Add(10, 10);
  dialogSizer->Add(new wxStaticText(this, -1,  wxString::FromAscii(title)), 0, 
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

  val = mpXCtrl->GetValue().ToAscii();
  mXYZ.x = strtod(val, &strptr);
  if (val == strptr) {wxMessageBox(_T("Invalid value for X!"), _T("Position error"), wxOK, this); return;}

  val = mpYCtrl->GetValue().ToAscii();
  mXYZ.y = strtod(val, &strptr);
  if (val == strptr) {wxMessageBox(_T("Invalid value for Y!"), _T("Position error"), wxOK, this); return;}

  val = mpZCtrl->GetValue().ToAscii();
  mXYZ.z = strtod(val, &strptr);
  if (val == strptr) {wxMessageBox(_T("Invalid value for Z!"), _T("Position error"), wxOK, this); return;}

    // close the dialog
    EndModal(wxID_OK);
}

Triple UserXYZBox::GetXYZ () {
  return mXYZ;
}


#endif // #ifdef OBJCRYST_GL

}// namespace 
