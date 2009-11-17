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
#include <algorithm>

#include "cctbx/sgtbx/space_group.h"

// wx headers, with or without precompilation
#include "wx/wxprec.h"
#ifdef __BORLANDC__
    #pragma hdrstop
#endif
#ifndef WX_PRECOMP
    #include "wx/wx.h"
#endif
#include "wx/dcbuffer.h"
#include "wx/config.h"
#include "wx/notebook.h"
#include "wx/progdlg.h"
#include "wx/filename.h"

#include "wxCryst/wxPowderPattern.h"
#include "wxCryst/wxRadiation.h"
#include "RefinableObj/Simplex.h"
#include "RefinableObj/LSQNumObj.h"
#include "ObjCryst/PowderPatternBackgroundBayesianMinimiser.h"
#include "Quirks/VFNStreamFormat.h"
#include "Quirks/Chronometer.h"

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

//#define USE_BACKGROUND_MAXLIKE_ERROR

namespace ObjCryst
{
////////////////////////////////////////////////////////////////////////
//
//    WXRadiation
//
////////////////////////////////////////////////////////////////////////
BEGIN_EVENT_TABLE(WXRadiation, wxWindow)
   EVT_UPDATE_UI(ID_CRYST_UPDATEUI, WXRadiation::OnUpdateUI)
END_EVENT_TABLE()

WXRadiation::WXRadiation(wxWindow *parent, Radiation* rad):
WXCrystObjBasic(parent),mpRadiation(rad)
{
   VFN_DEBUG_ENTRY("WXRadiation::WXRadiation()",6)
   // :TODO: Add a choice for the wavlength type, with 'monochromatic', and a list
   // of X-Ray tubes.
   mpSizer=new wxBoxSizer(wxHORIZONTAL);
   
   mpFieldRadType= new WXFieldOption(this,-1,&(mpRadiation->mRadiationType));
   mpSizer->Add(mpFieldRadType,0);
   mList.Add(mpFieldRadType);
   
   mpFieldWavelengthType= new WXFieldOption(this,-1,&(mpRadiation->mWavelengthType));
   mpSizer->Add(mpFieldWavelengthType,0);
   mList.Add(mpFieldWavelengthType);

   WXCrystObjBasic* pFieldWavelength
      =mpRadiation->GetPar(mpRadiation->mWavelength.data()).WXCreate(this);
   mpSizer->Add(pFieldWavelength,0);
   mList.Add(pFieldWavelength);

   WXFieldPar<REAL> *polarRate=new WXFieldPar<REAL>(this,"Linear Polar Rate:",-1,
                                            &(mpRadiation->mLinearPolarRate));
   mpSizer->Add(polarRate,0,wxALIGN_LEFT);
   mList.Add(polarRate);
   
   WXFieldPar<REAL> *xRayTubeDlambda=new WXFieldPar<REAL>(this,"Tube-DeltaLambda:",-1,
                                                &(mpRadiation->mXRayTubeDeltaLambda));
   mpSizer->Add(xRayTubeDlambda,0,wxALIGN_LEFT);
   mList.Add(xRayTubeDlambda);
      
   WXFieldPar<REAL> *xRayTubeAlpha2Alpha1=new WXFieldPar<REAL>(this,"Tube-Alpha2/Alpha1:",-1,
                                            &(mpRadiation->mXRayTubeAlpha2Alpha1Ratio));
   mpSizer->Add(xRayTubeAlpha2Alpha1,0,wxALIGN_LEFT);
   mList.Add(xRayTubeAlpha2Alpha1);
      
   this->CrystUpdate(true);
   this->SetSizer(mpSizer);
   mpSizer->SetSizeHints(this);
   this->Layout();
   VFN_DEBUG_EXIT("WXRadiation::WXRadiation()",6)
}
WXRadiation::~WXRadiation()
{
   mpRadiation->WXNotifyDelete();
}

void WXRadiation::CrystUpdate(const bool uui,const bool lock)
{
   if(lock) mMutex.Lock();
   mList.CrystUpdate(false,false);
   if(lock) mMutex.Unlock();
   if(uui)
   {
      if(true==wxThread::IsMain()) this->UpdateUI(lock);
      else
      {
         wxUpdateUIEvent event(ID_CRYST_UPDATEUI);
         wxPostEvent(this,event);
      }
   }
}
void WXRadiation::UpdateUI(const bool lock)
{
   mList.UpdateUI(lock);
}
void WXRadiation::OnUpdateUI(wxUpdateUIEvent& event)
{
   this->UpdateUI(true);
   event.Skip();
}

//////////////////////////////////////// WXProfileFitting /////////////////////

class WXProfileFitting:public wxWindow
{
   public:
      WXProfileFitting(wxWindow *parent,PowderPattern *pPattern,PowderPatternDiffraction *pDiff=0);
      ~WXProfileFitting();
      /// Start Fitting
      void OnFit(wxCommandEvent &event);
      void OnExploreSpacegroups(wxCommandEvent &event);
   private:
      PowderPattern *mpPattern;
      PowderPatternDiffraction *mpDiff;
      wxCheckListBox *mpFitCheckList;
      wxTextCtrl *mpLog;
      wxListBox *mpList;
      LSQNumObj mLSQ;
      DECLARE_EVENT_TABLE()
};

////////////////////////////////////////////////////////////////////////
//
//    WXPowderPattern
//
////////////////////////////////////////////////////////////////////////
static const long ID_POWDER_MENU_COMP_ADDBACKGD_BAYESIAN=WXCRYST_ID(); 
static const long ID_POWDER_MENU_COMP_ADDBACKGD=       WXCRYST_ID(); 
static const long ID_POWDER_MENU_COMP_ADDCRYST=        WXCRYST_ID(); 
static const long ID_POWDER_MENU_GRAPH=                     WXCRYST_ID(); 
static const long ID_POWDER_MENU_SAVETEXT=                  WXCRYST_ID(); 
static const long ID_POWDER_MENU_SIMULATE=                  WXCRYST_ID(); 
static const long ID_POWDER_MENU_EXPORT=                    WXCRYST_ID(); 
static const long ID_POWDER_MENU_EXPORT_FULLPROF=           WXCRYST_ID(); 
static const long ID_POWDER_MENU_IMPORT_FULLPROF=           WXCRYST_ID(); 
static const long ID_POWDER_MENU_IMPORT_PSI_DMC=            WXCRYST_ID(); 
static const long ID_POWDER_MENU_IMPORT_ILL_D1A5=           WXCRYST_ID(); 
static const long ID_POWDER_MENU_IMPORT_XDD=                WXCRYST_ID(); 
static const long ID_POWDER_MENU_IMPORT_CPI=                WXCRYST_ID(); 
static const long ID_POWDER_MENU_IMPORT_FULLPROF4=          WXCRYST_ID(); 
static const long ID_POWDER_MENU_IMPORT_MULTIDETECTORLLBG42=WXCRYST_ID(); 
static const long ID_POWDER_MENU_IMPORT_2THETAOBSSIGMA=     WXCRYST_ID(); 
static const long ID_POWDER_MENU_IMPORT_2THETAOBS=          WXCRYST_ID(); 
static const long ID_POWDER_MENU_IMPORT_TOFISISXYSIGMA=     WXCRYST_ID(); 
static const long ID_POWDER_MENU_IMPORT_GSAS=               WXCRYST_ID(); 
static const long ID_POWDER_MENU_IMPORT_CIF=                WXCRYST_ID(); 
static const long ID_POWDER_MENU_FITSCALE_R=                WXCRYST_ID(); 
static const long ID_POWDER_MENU_FITSCALE_RW=               WXCRYST_ID(); 
static const long ID_POWDER_MENU_WAVELENGTH=                WXCRYST_ID(); 
static const long ID_POWDER_MENU_WAVELENGTH_XRAY=           WXCRYST_ID(); 
static const long ID_POWDER_MENU_WAVELENGTH_NEUTRON=        WXCRYST_ID(); 
static const long ID_POWDER_MENU_WAVELENGTH_NEUTRON_TOF=    WXCRYST_ID(); 
static const long ID_POWDER_MENU_WAVELENGTH_SET=            WXCRYST_ID(); 
static const long ID_POWDER_MENU_WAVELENGTH_SET_AG=         WXCRYST_ID(); 
static const long ID_POWDER_MENU_WAVELENGTH_SET_MO=         WXCRYST_ID(); 
static const long ID_POWDER_MENU_WAVELENGTH_SET_CU=         WXCRYST_ID(); 
static const long ID_POWDER_MENU_WAVELENGTH_SET_FE=         WXCRYST_ID(); 
static const long ID_POWDER_MENU_WAVELENGTH_SET_CO=         WXCRYST_ID(); 
static const long ID_POWDER_MENU_WAVELENGTH_SET_CR=         WXCRYST_ID(); 
static const long ID_POWDER_MENU_WAVELENGTH_SET_AGA1=       WXCRYST_ID(); 
static const long ID_POWDER_MENU_WAVELENGTH_SET_MOA1=       WXCRYST_ID(); 
static const long ID_POWDER_MENU_WAVELENGTH_SET_CUA1=       WXCRYST_ID(); 
static const long ID_POWDER_MENU_WAVELENGTH_SET_FEA1=       WXCRYST_ID(); 
static const long ID_POWDER_MENU_WAVELENGTH_SET_COA1=       WXCRYST_ID(); 
static const long ID_POWDER_MENU_WAVELENGTH_SET_CRA1=       WXCRYST_ID(); 
static const long ID_POWDER_MENU_ADD_2THETA_EXCLUDE=        WXCRYST_ID(); 
static const long ID_POWDER_MENU_LEBAIL=                    WXCRYST_ID();
static const long ID_POWDERBACKGROUND_IMPORT=               WXCRYST_ID(); 
static const long ID_POWDERBACKGROUND_OPTIMIZEBAYESIAN=     WXCRYST_ID(); 
static const long ID_POWDERDIFF_CRYSTAL=                    WXCRYST_ID(); 
static const long ID_POWDERDIFF_SAVEHKLFCALC=               WXCRYST_ID(); 
static const long ID_POWDER_GRAPH_NEW_PATTERN=              WXCRYST_ID(); 
static const long ID_POWDERTEXTURE_MENU_ADDPHASE=                   WXCRYST_ID(); 
static const long ID_POWDERTEXTURE_MENU_DELETEPHASE=                WXCRYST_ID(); 
static const long ID_POWDERPATTERN_MENU_COMPONENTS=                 WXCRYST_ID(); 
static const long ID_POWDERPATTERN_MENU_PATTERN=                    WXCRYST_ID(); 

static const long ID_POWDERDIFF_PROFILE_DEPV=                  WXCRYST_ID();


BEGIN_EVENT_TABLE(WXPowderPattern, wxWindow)
   EVT_BUTTON(ID_WXOBJ_COLLAPSE,                    WXCrystObj::OnToggleCollapse)                
   EVT_MENU(ID_POWDER_MENU_EXPORT_FULLPROF,         WXPowderPattern::OnMenuExport)                  
   EVT_MENU(ID_REFOBJ_MENU_OBJ_SAVE,                WXRefinableObj::OnMenuSave)                  
   EVT_MENU(ID_REFOBJ_MENU_OBJ_LOAD,                WXRefinableObj::OnMenuLoad)                  
   EVT_MENU(ID_REFOBJ_MENU_PAR_FIXALL,              WXRefinableObj::OnMenuFixAllPar)             
   EVT_MENU(ID_REFOBJ_MENU_PAR_UNFIXALL,            WXRefinableObj::OnMenuUnFixAllPar)           
   EVT_MENU(ID_POWDER_MENU_COMP_ADDBACKGD,          WXPowderPattern::OnMenuAddCompBackgd)        
   EVT_MENU(ID_POWDER_MENU_COMP_ADDBACKGD_BAYESIAN, WXPowderPattern::OnMenuAddCompBackgdBayesian)   
   EVT_MENU(ID_POWDER_MENU_COMP_ADDCRYST,           WXPowderPattern::OnMenuAddCompCryst)         
   EVT_MENU(ID_POWDER_MENU_SAVETEXT,                WXPowderPattern::OnMenuSaveText)             
   EVT_MENU(ID_POWDER_MENU_SIMULATE,                WXPowderPattern::OnMenuSimulate)             
   EVT_MENU(ID_POWDER_MENU_IMPORT_FULLPROF,         WXPowderPattern::OnMenuImportPattern)       
   EVT_MENU(ID_POWDER_MENU_IMPORT_PSI_DMC,          WXPowderPattern::OnMenuImportPattern)            
   EVT_MENU(ID_POWDER_MENU_IMPORT_ILL_D1A5,         WXPowderPattern::OnMenuImportPattern)            
   EVT_MENU(ID_POWDER_MENU_IMPORT_XDD,              WXPowderPattern::OnMenuImportPattern)            
   EVT_MENU(ID_POWDER_MENU_IMPORT_CPI,              WXPowderPattern::OnMenuImportPattern)            
   EVT_MENU(ID_POWDER_MENU_IMPORT_FULLPROF4,        WXPowderPattern::OnMenuImportPattern)      
   EVT_MENU(ID_POWDER_MENU_IMPORT_MULTIDETECTORLLBG42,WXPowderPattern::OnMenuImportPattern)
   EVT_MENU(ID_POWDER_MENU_IMPORT_2THETAOBSSIGMA,   WXPowderPattern::OnMenuImportPattern)
   EVT_MENU(ID_POWDER_MENU_IMPORT_2THETAOBS,        WXPowderPattern::OnMenuImportPattern)    
   EVT_MENU(ID_POWDER_MENU_IMPORT_TOFISISXYSIGMA,   WXPowderPattern::OnMenuImportPattern)    
   EVT_MENU(ID_POWDER_MENU_IMPORT_GSAS,             WXPowderPattern::OnMenuImportPattern)    
   EVT_MENU(ID_POWDER_MENU_IMPORT_CIF,              WXPowderPattern::OnMenuImportPattern)    
   EVT_MENU(ID_POWDER_MENU_WAVELENGTH_SET,          WXPowderPattern::OnMenuSetWavelength)      
   EVT_MENU(ID_POWDER_MENU_WAVELENGTH_XRAY,         WXPowderPattern::OnMenuSetWavelength)      
   EVT_MENU(ID_POWDER_MENU_WAVELENGTH_NEUTRON,      WXPowderPattern::OnMenuSetWavelength)      
   EVT_MENU(ID_POWDER_MENU_WAVELENGTH_NEUTRON_TOF,  WXPowderPattern::OnMenuSetWavelength)      
   EVT_MENU(ID_POWDER_MENU_WAVELENGTH_SET_AG,       WXPowderPattern::OnMenuSetWavelength)      
   EVT_MENU(ID_POWDER_MENU_WAVELENGTH_SET_MO,       WXPowderPattern::OnMenuSetWavelength)      
   EVT_MENU(ID_POWDER_MENU_WAVELENGTH_SET_CU,       WXPowderPattern::OnMenuSetWavelength)      
   EVT_MENU(ID_POWDER_MENU_WAVELENGTH_SET_FE,       WXPowderPattern::OnMenuSetWavelength)      
   EVT_MENU(ID_POWDER_MENU_WAVELENGTH_SET_CO,       WXPowderPattern::OnMenuSetWavelength)      
   EVT_MENU(ID_POWDER_MENU_WAVELENGTH_SET_CR,       WXPowderPattern::OnMenuSetWavelength)      
   EVT_MENU(ID_POWDER_MENU_WAVELENGTH_SET_AGA1,     WXPowderPattern::OnMenuSetWavelength)      
   EVT_MENU(ID_POWDER_MENU_WAVELENGTH_SET_MOA1,     WXPowderPattern::OnMenuSetWavelength)      
   EVT_MENU(ID_POWDER_MENU_WAVELENGTH_SET_CUA1,     WXPowderPattern::OnMenuSetWavelength)      
   EVT_MENU(ID_POWDER_MENU_WAVELENGTH_SET_FEA1,     WXPowderPattern::OnMenuSetWavelength)      
   EVT_MENU(ID_POWDER_MENU_WAVELENGTH_SET_COA1,     WXPowderPattern::OnMenuSetWavelength)      
   EVT_MENU(ID_POWDER_MENU_WAVELENGTH_SET_CRA1,     WXPowderPattern::OnMenuSetWavelength)      
   EVT_MENU(ID_POWDER_MENU_GRAPH,                   WXPowderPattern::OnMenuShowGraph)          
   EVT_MENU(ID_POWDER_MENU_FITSCALE_R,              WXPowderPattern::OnMenuFitScaleForR)       
   EVT_MENU(ID_POWDER_MENU_FITSCALE_RW,             WXPowderPattern::OnMenuFitScaleForRw)      
   EVT_MENU(ID_POWDER_MENU_ADD_2THETA_EXCLUDE,      WXPowderPattern::OnMenuAddExclude)   
   EVT_MENU(ID_POWDER_MENU_LEBAIL,                  WXPowderPattern::OnMenuLeBail)
   EVT_UPDATE_UI(ID_CRYST_UPDATEUI,                 WXRefinableObj::OnUpdateUI)                
END_EVENT_TABLE()

WXPowderPattern::WXPowderPattern(wxWindow *parent, PowderPattern* pow):
WXRefinableObj(parent,pow),mpPowderPattern(pow),mpGraph(0),
mChi2(0.0),mGoF(0.0),mRwp(0.0),mRp(0.0)
{
   VFN_DEBUG_MESSAGE("WXPowderPattern::WXPowderPattern()",6)
   mpWXTitle->SetForegroundColour(wxColour(255,0,0));
   mpWXTitle->SetSize(400,-1);
   // Menu
      mpMenuBar->AddMenu("Data",ID_REFOBJ_MENU_OBJ);
         //:TODO: reactivate & test those menus
         //mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_REFOBJ_MENU_OBJ_SAVE,"Save");
         //mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_REFOBJ_MENU_OBJ_LOAD,"Load");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_POWDER_MENU_SAVETEXT,
                                "Save pattern (text)");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_POWDER_MENU_SIMULATE,
                                "Simulation mode (no obs. pattern)");
         mpMenuBar->GetMenu(ID_REFOBJ_MENU_OBJ).AppendSeparator();
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_POWDER_MENU_IMPORT_CIF,
                                 "Import CIF Powder Data");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_POWDER_MENU_IMPORT_GSAS,
                                 "Import GSAS Data(CONS-ESD,CONS6STD,RALF-ALT)");
         mpMenuBar->GetMenu(ID_REFOBJ_MENU_OBJ).AppendSeparator();
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_POWDER_MENU_IMPORT_2THETAOBSSIGMA,
                                 "Import 2Theta-Obs-Sigma Pattern");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_POWDER_MENU_IMPORT_2THETAOBS,
                                 "Import 2Theta-Obs Pattern");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_POWDER_MENU_IMPORT_FULLPROF,
                                 "Import Fullprof Pattern");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_POWDER_MENU_IMPORT_FULLPROF4,
                                 "Import FullProf format #4");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_POWDER_MENU_IMPORT_XDD,
                                 "Import Xdd Pattern");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_POWDER_MENU_IMPORT_CPI,
                                 "Import Sietronics CPI Pattern");
         mpMenuBar->GetMenu(ID_REFOBJ_MENU_OBJ).AppendSeparator();
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_POWDER_MENU_IMPORT_ILL_D1A5,
                                 "Import Neutron ILL(D1A-D1B) Pattern (D1A5)");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_POWDER_MENU_IMPORT_PSI_DMC,
                                 "Import PSI(DMC) Pattern");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_POWDER_MENU_IMPORT_MULTIDETECTORLLBG42,
                                 "Import Neutron Multi-Detector Format (LLB G42)");
         mpMenuBar->GetMenu(ID_REFOBJ_MENU_OBJ).AppendSeparator();
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_POWDER_MENU_IMPORT_TOFISISXYSIGMA,
                                 "Import Neutron TOF ISIS X Y Sigma");
      mpMenuBar->AddMenu("Export",ID_POWDER_MENU_EXPORT);
         mpMenuBar->AddMenuItem(ID_POWDER_MENU_EXPORT,ID_POWDER_MENU_EXPORT_FULLPROF,
                                "Export to Fullprof");
      mpMenuBar->AddMenu("Parameters",ID_REFOBJ_MENU_PAR);
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_PAR,ID_REFOBJ_MENU_PAR_FIXALL,"Fix all");
         //mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_PAR,ID_REFOBJ_MENU_PAR_UNFIXALL,"Unfix all");
      mpMenuBar->AddMenu("Phases",ID_POWDERPATTERN_MENU_COMPONENTS);
         mpMenuBar->AddMenuItem(ID_POWDERPATTERN_MENU_COMPONENTS,
                                ID_POWDER_MENU_COMP_ADDBACKGD_BAYESIAN,
                                "Add Background (Bayesian, automatic)");
         mpMenuBar->AddMenuItem(ID_POWDERPATTERN_MENU_COMPONENTS,
                                ID_POWDER_MENU_COMP_ADDBACKGD,
                                "Add user-supplied Background ");
         mpMenuBar->AddMenuItem(ID_POWDERPATTERN_MENU_COMPONENTS,ID_POWDER_MENU_COMP_ADDCRYST,
                                "Add Crystalline Phase");
      mpMenuBar->AddMenu("Radiation",ID_POWDER_MENU_WAVELENGTH);
         mpMenuBar->AddMenuItem(ID_POWDER_MENU_WAVELENGTH,
                                ID_POWDER_MENU_WAVELENGTH_NEUTRON,
                                "Neutron");
         mpMenuBar->AddMenuItem(ID_POWDER_MENU_WAVELENGTH,
                                ID_POWDER_MENU_WAVELENGTH_NEUTRON_TOF,
                                "Neutron Time Of Flight");
         mpMenuBar->AddMenuItem(ID_POWDER_MENU_WAVELENGTH,
                                ID_POWDER_MENU_WAVELENGTH_XRAY,
                                "X-Rays");
         mpMenuBar->AddMenuItem(ID_POWDER_MENU_WAVELENGTH,
                                ID_POWDER_MENU_WAVELENGTH_SET,
                                "Monochromatic Wavelength");
         mpMenuBar->AddMenuItem(ID_POWDER_MENU_WAVELENGTH,
                                ID_POWDER_MENU_WAVELENGTH_SET_AG,
                                "X-Ray Tube Ag Ka12");
         mpMenuBar->AddMenuItem(ID_POWDER_MENU_WAVELENGTH,
                                ID_POWDER_MENU_WAVELENGTH_SET_AGA1,
                                "X-Ray Tube Ag Ka1");
         mpMenuBar->AddMenuItem(ID_POWDER_MENU_WAVELENGTH,
                                ID_POWDER_MENU_WAVELENGTH_SET_MO,
                                "X-Ray Tube Mo Ka12");
         mpMenuBar->AddMenuItem(ID_POWDER_MENU_WAVELENGTH,
                                ID_POWDER_MENU_WAVELENGTH_SET_MOA1,
                                "X-Ray Tube Mo Ka1");
         mpMenuBar->AddMenuItem(ID_POWDER_MENU_WAVELENGTH,
                                ID_POWDER_MENU_WAVELENGTH_SET_CU,
                                "X-Ray Tube Cu Ka12");
         mpMenuBar->AddMenuItem(ID_POWDER_MENU_WAVELENGTH,
                                ID_POWDER_MENU_WAVELENGTH_SET_CUA1,
                                "X-Ray Tube Cu Ka1");
         mpMenuBar->AddMenuItem(ID_POWDER_MENU_WAVELENGTH,
                                ID_POWDER_MENU_WAVELENGTH_SET_FE,
                                "X-Ray Tube Fe Ka12");
         mpMenuBar->AddMenuItem(ID_POWDER_MENU_WAVELENGTH,
                                ID_POWDER_MENU_WAVELENGTH_SET_FEA1,
                                "X-Ray Tube Fe Ka1");
         mpMenuBar->AddMenuItem(ID_POWDER_MENU_WAVELENGTH,
                                ID_POWDER_MENU_WAVELENGTH_SET_CO,
                                "X-Ray Tube Co Ka12");
         mpMenuBar->AddMenuItem(ID_POWDER_MENU_WAVELENGTH,
                                ID_POWDER_MENU_WAVELENGTH_SET_COA1,
                                "X-Ray Tube Co Ka1");
         mpMenuBar->AddMenuItem(ID_POWDER_MENU_WAVELENGTH,
                                ID_POWDER_MENU_WAVELENGTH_SET_CR,
                                "X-Ray Tube Cr Ka12");
         mpMenuBar->AddMenuItem(ID_POWDER_MENU_WAVELENGTH,
                                ID_POWDER_MENU_WAVELENGTH_SET_CRA1,
                                "X-Ray Tube Cr Ka1");
      mpMenuBar->AddMenu("Pattern",ID_POWDERPATTERN_MENU_PATTERN);
         mpMenuBar->AddMenuItem(ID_POWDERPATTERN_MENU_PATTERN,ID_POWDER_MENU_GRAPH,
                                "Show Graph");
         mpMenuBar->AddMenuItem(ID_POWDERPATTERN_MENU_PATTERN,ID_POWDER_MENU_FITSCALE_R,
                                "Fit Scale for R");
         mpMenuBar->AddMenuItem(ID_POWDERPATTERN_MENU_PATTERN,ID_POWDER_MENU_FITSCALE_RW,
                                "Fit Scale for Rw");
         mpMenuBar->AddMenuItem(ID_POWDERPATTERN_MENU_PATTERN,
                                ID_POWDER_MENU_ADD_2THETA_EXCLUDE,
                                "Add excluded region");
         mpMenuBar->GetMenu(ID_POWDERPATTERN_MENU_PATTERN).AppendSeparator();
         mpMenuBar->AddMenuItem(ID_POWDERPATTERN_MENU_PATTERN,
                                ID_POWDER_MENU_LEBAIL,
                                "Fit profile + Le Bail extract");
      mpSizer->SetItemMinSize(mpMenuBar,
                              mpMenuBar->GetSize().GetWidth(),
                              mpMenuBar->GetSize().GetHeight());
   //Radiation
      mpSizer->Add(mpPowderPattern->mRadiation.WXCreate(this),0);
      mList.Add(mpPowderPattern->mRadiation.WXGet());
   // Correction to 2Theta
      wxBoxSizer* thetaCorrSizer=new wxBoxSizer(wxHORIZONTAL);

      WXCrystObjBasic* fieldZero    
         =mpPowderPattern->GetPar(&(mpPowderPattern->mXZero)).WXCreate(this);
      fieldZero->SetToolTip(_T("Zero shift of peaks\n")
                            _T("2Theta = 2Theta_Bragg + Zero\n"));
      WXCrystObjBasic* fieldThetaDispl
         =mpPowderPattern->GetPar(&(mpPowderPattern->m2ThetaDisplacement)).WXCreate(this);
      fieldThetaDispl->SetToolTip(_T("Peak shift due to sample displacement:\n")
                                  _T("2Theta = 2Theta_Bragg + Displacement/cos(Theta)"));
      WXCrystObjBasic* fieldThetaTransp
         =mpPowderPattern->GetPar(&(mpPowderPattern->m2ThetaTransparency)).WXCreate(this);
      fieldThetaTransp->SetToolTip(_T("Zero shift of the peak 2theta positions\n")
                                   _T("2Theta = 2Theta_Bragg + Transparency*sin(Theta)"));

      thetaCorrSizer->Add(fieldZero,0);
      thetaCorrSizer->Add(fieldThetaDispl,0);
      thetaCorrSizer->Add(fieldThetaTransp,0);
      mList.Add(fieldZero);
      mList.Add(fieldThetaDispl);
      mList.Add(fieldThetaTransp);
      mpSizer->Add(thetaCorrSizer);
   // Time OF Flight parameters
      wxBoxSizer* tofSizer=new wxBoxSizer(wxHORIZONTAL);
      WXCrystObjBasic* fieldDIFC=mpPowderPattern->GetPar(&(mpPowderPattern->mDIFC)).WXCreate(this);
      WXCrystObjBasic* fieldDIFA=mpPowderPattern->GetPar(&(mpPowderPattern->mDIFA)).WXCreate(this);
      fieldDIFA->SetToolTip(_T("Peak position (time, in microseconds):\n")
                            _T("t = DIFA * d_hkl + DIFC * d_hkl^2 + ZERO"));
      fieldDIFC->SetToolTip(_T("Peak position (time, in microseconds):\n")
                            _T("t = DIFA * d_hkl + DIFC * d_hkl^2 + ZERO"));
      tofSizer->Add(fieldDIFC,0);
      tofSizer->Add(fieldDIFA,0);
      mList.Add(fieldDIFC);
      mList.Add(fieldDIFA);
      mpSizer->Add(tofSizer);
   // Max Sin(theta/Lambda)
      WXFieldPar<REAL> *maxSiThOvLa=
         new WXFieldPar<REAL>(this,"Max Sin(theta)/lambda:",-1,
                              &(mpPowderPattern->mMaxSinThetaOvLambda));
      mpSizer->Add(maxSiThOvLa,0,wxALIGN_LEFT);
      mList.Add(maxSiThOvLa);
      maxSiThOvLa->SetToolTip(_T("Maximum sin(theta)/lambda=1/2d\n")
                              _T("For global optimization, the default value of ")
                              _T("0.25 (2A resolution) should be sufficient.\n")
                              _T("Use larger values if necessary (0.4(1.25A), 0.5(1A))")
                              _T("but keep in mind that the number of reflections (and")
                              _T("therefore the computing time) varies as [sin(theta/lambda)]^3..."));
   // Statistics
      wxBoxSizer* pStats=new wxBoxSizer(wxHORIZONTAL);
      
      WXFieldPar<REAL> *pWXFieldChi2=new WXFieldPar<REAL>(this,"Chi^2",-1,&mChi2,140);
      pStats->Add(pWXFieldChi2    ,0,wxALIGN_CENTER);
      mList.Add(pWXFieldChi2);
      pWXFieldChi2->SetToolTip(_T("Chi^2=SUM[(Obs_i-Calc_i)^2/Sigma_i^2]"));
      dynamic_cast<WXFieldParBase *>(pWXFieldChi2)->SetFormat(_T("%10.2f"));
      
      WXFieldPar<REAL> *pWXFieldGof=new WXFieldPar<REAL>(this,"GoF",-1,&mGoF,90);
      pStats->Add(pWXFieldGof    ,0,wxALIGN_CENTER);
      mList.Add(pWXFieldGof);
      pWXFieldGof->SetToolTip(_T("GoF=Chi^2/NbPoints"));
      dynamic_cast<WXFieldParBase *>(pWXFieldGof)->SetFormat(_T("%8.3f"));
      
      WXFieldPar<REAL> *pWXFieldRwp=new WXFieldPar<REAL>(this,"Rwp",-1,&mRwp,70);
      pStats->Add(pWXFieldRwp    ,0,wxALIGN_CENTER);
      mList.Add(pWXFieldRwp);
      pWXFieldRwp->SetToolTip(_T("Full profile R-factor (weighted)\n")
                              _T("Will use integrated profiles if option is set."));
      dynamic_cast<WXFieldParBase *>(pWXFieldRwp)->SetFormat(_T("%8.4f"));
      
      WXFieldPar<REAL> *pWXFieldRp=new WXFieldPar<REAL>(this,"Rp",-1,&mRp,70);
      pStats->Add(pWXFieldRp    ,0,wxALIGN_CENTER);
      mList.Add(pWXFieldRp);
      pWXFieldRp->SetToolTip(_T("Full profile R-factor (unweighted)\n")
                             _T("Will use integrated profiles if option is set."));
      dynamic_cast<WXFieldParBase *>(pWXFieldRp)->SetFormat(_T("%8.4f"));
      //pStats->SetSizeHints(this);
      //pStats->Layout();
      
      mpSizer->Add(pStats);
   // Components
      mpWXComponent=mpPowderPattern
                    ->mPowderPatternComponentRegistry.WXCreate(this);
      mpSizer->Add(mpWXComponent,0,wxALIGN_LEFT);
      mList.Add(mpWXComponent);
   
   VFN_DEBUG_MESSAGE("WXPowderPattern::WXPowderPattern():1",6)
   this->BottomLayout(0);
   this->CrystUpdate(true);
   {
      if(!wxConfigBase::Get()->HasEntry(_T("PowderPattern/BOOL/Automatically open powder pattern graph")))
         wxConfigBase::Get()->Write(_T("PowderPattern/BOOL/Automatically open powder pattern graph"), false);
      else
      {
         bool val;
         wxConfigBase::Get()->Read(_T("PowderPattern/BOOL/Automatically open powder pattern graph"), &val);
         if(val)
         {
            wxCommandEvent event(wxEVT_COMMAND_MENU_SELECTED,ID_POWDER_MENU_GRAPH);
            wxPostEvent(this,event);
         }
      }
   }
  VFN_DEBUG_MESSAGE("WXPowderPattern::WXPowderPattern():End",6)
}

void WXPowderPattern::CrystUpdate(const bool uui,const bool lock)
{
   VFN_DEBUG_ENTRY("WXPowderPattern::CrystUpdate()",7)
   wxWakeUpIdle();
   if(lock) mMutex.Lock();
   WXCrystValidateAllUserInput();
   
   if(mpPowderPattern->GetNbPoint()<=0)
   {
      if(lock) mMutex.Unlock();
      this->WXRefinableObj::CrystUpdate(uui,lock);
      return;// nothing to display yet
   }
   
   // Will force re-generating reflection list if the wavelength,
   // or lattice par, or the spacegroup has changed.
   VFN_DEBUG_MESSAGE("WXPowderPattern::CrystUpdate()",7)
   mpPowderPattern->Prepare();
   VFN_DEBUG_MESSAGE("WXPowderPattern::CrystUpdate()",7)
   
   mChi2=mpPowderPattern->GetChi2_Option();
   if(mpPowderPattern->mNbPointUsed>0)
      mGoF=mpPowderPattern->GetChi2()/mpPowderPattern->mNbPointUsed;
   else mGoF=0;
   //cout<<"WXPowderPattern::CrystUpdate():"<<mpPowderPattern->GetChi2()<<"/"<<mpPowderPattern->mNbPointUsed<<"="<<mGoF<<endl;
   mRwp=mpPowderPattern->GetRw();
   mRp=mpPowderPattern->GetR();
   
   if(mpGraph!=0)
   {
      CrystVector_REAL tmp;
      mpPowderPattern->CalcPowderPattern();
      tmp=mpPowderPattern->GetPowderPatternVariance();
      for(long i=0;i<tmp.numElements();i++)
      {
         if(tmp(i)<0) tmp(i)=0;
         else tmp(i)=sqrt(tmp(i));
      }
      mpGraph->SetPattern( mpPowderPattern->GetPowderPatternX(),
                           mpPowderPattern->GetPowderPatternObs(),
                           mpPowderPattern->GetPowderPatternCalc(),
                           tmp,
                           mpPowderPattern->GetChi2Cumul());
   }
   if(lock) mMutex.Unlock();
   this->WXRefinableObj::CrystUpdate(uui,lock);
   VFN_DEBUG_EXIT("WXPowderPattern::CrystUpdate()",7)
} 

void WXPowderPattern::OnMenuAddCompBackgd(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXPowderPattern::OnMenuAddCompBackgd()",6)
   WXCrystValidateAllUserInput();
   const unsigned int nb=mpPowderPattern->GetNbPowderPatternComponent();
   bool hasBack=false;
   for(unsigned int i=0;i<nb;i++)
      if(mpPowderPattern->GetPowderPatternComponent(i).GetClassName()=="PowderPatternBackground")
      {
         hasBack=true;
         break;
      }
   if(hasBack)
   {
      wxMessageDialog dialog(this,_T("You already have one background !\n")
                                  _T(" Are you sure you want to add one ?"),
                             _T("Warning : Duplicate Background !"),
                             wxYES_NO|wxICON_HAND|wxNO_DEFAULT);
      if(wxID_YES!=dialog.ShowModal())
         return;
   }
   PowderPatternBackground *backgdData= new PowderPatternBackground;
   mpPowderPattern->AddPowderPatternComponent(*backgdData);
   if(mpGraph!=0) mpPowderPattern->Prepare();//else this will be done when opening the graph
   //this->Layout();
}

void WXPowderPattern::OnMenuAddCompBackgdBayesian(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_ENTRY("WXPowderPattern::OnMenuAddCompBackgdBayesian()",6)
   WXCrystValidateAllUserInput();
   const unsigned int nb=mpPowderPattern->GetNbPowderPatternComponent();
   bool hasBack=false;
   for(unsigned int i=0;i<nb;i++)
      if(mpPowderPattern->GetPowderPatternComponent(i).GetClassName()=="PowderPatternBackground")
      {
         hasBack=true;
         break;
      }
   if(hasBack)
   {
      wxMessageDialog dialog(this,_T("You already have one background !\n")
                                  _T(" Are you sure you want to add one ?"),
                             _T("Warning : Duplicate Background !"),
                             wxYES_NO|wxICON_HAND|wxNO_DEFAULT);
      if(wxID_YES!=dialog.ShowModal())
         return;
   }

   long nbPointSpline=20;
   wxString mes=_T("Number of Interpolation Points");
   wxString s;
   s.Printf(_T("%i"),nbPointSpline);
   wxTextEntryDialog dialog(this,mes,_T("Automatic Bayesian (David-Sivia) Background"),
                            s,wxOK | wxCANCEL);
   if(wxID_OK!=dialog.ShowModal())
   {
      VFN_DEBUG_EXIT("WXPowderPattern::OnMenuAddCompBackgdBayesian():Canceled",6)
      return;
   }
   dialog.GetValue().ToLong(&nbPointSpline);
   if(nbPointSpline<=1)nbPointSpline=2;
   
   wxProgressDialog dlgProgress(_T("Automatic Bayesian Background"),_T("Automatic Background: Initializing..."),
                                4,this,wxPD_AUTO_HIDE|wxPD_ELAPSED_TIME|wxPD_CAN_ABORT);

   PowderPatternBackground *pBckgd= new PowderPatternBackground;
   VFN_DEBUG_MESSAGE("WXPowderPattern::OnMenuAddCompBackgdBayesian()",6)
   mpPowderPattern->AddPowderPatternComponent(*pBckgd);
   VFN_DEBUG_MESSAGE("WXPowderPattern::OnMenuAddCompBackgdBayesian()",6)
   {
      CrystVector_REAL x(nbPointSpline),backgd(nbPointSpline);
      const CrystVector_REAL *pObs=&(pBckgd->GetParentPowderPattern().GetPowderPatternObs());
      const unsigned long nbPoint=pBckgd->GetParentPowderPattern().GetNbPoint();
      const float xmin=pBckgd->GetParentPowderPattern().GetPowderPatternX()(0),
                  xmax=pBckgd->GetParentPowderPattern().GetPowderPatternX()(nbPoint-1);
      for(int i=0;i<nbPointSpline;i++)
      {// xmax is not necessarily > xmin, but in the right order (TOF)
         x(i)=xmin+(xmax-xmin)/(REAL)(nbPointSpline-1)*REAL(i);
         REAL x1=xmin+(xmax-xmin)/(REAL)(nbPointSpline-1)*REAL(i-.2);
         REAL x2=xmin+(xmax-xmin)/(REAL)(nbPointSpline-1)*REAL(i+.2);
         long n1=(long)(pBckgd->GetParentPowderPattern().X2Pixel(x1));
         long n2=(long)(pBckgd->GetParentPowderPattern().X2Pixel(x2));
         if(n1<0) n1=0;
         if(n2>(long)nbPoint)n2=nbPoint;
         backgd(i)=(*pObs)(n1);
         for(long j=n1;j<n2;j++)
            if((*pObs)(j)<backgd(i))backgd(i)=(*pObs)(j);
      }
      pBckgd->SetInterpPoints(x,backgd);
   }
   if(mpGraph!=0) mpPowderPattern->Prepare();//else this will be done when opening the graph
   
   pBckgd->UnFixAllPar();
   pBckgd->GetOption(0).SetChoice(0);//linear
   if(dlgProgress.Update(1,_T("Automatic Background: Optimizing Linear Model..."))==false) return;
   pBckgd->OptimizeBayesianBackground();
   pBckgd->GetOption(0).SetChoice(1);//spline
   if(dlgProgress.Update(2,_T("Automatic Background: Optimizing Spline Model..."))==false) return;
   pBckgd->OptimizeBayesianBackground();
   pBckgd->FixAllPar();

   //this->Layout();
   VFN_DEBUG_EXIT("WXPowderPattern::OnMenuAddCompBackgdBayesian()",6)
}

void WXPowderPattern::OnMenuAddCompCryst(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_ENTRY("WXPowderPattern::OnMenuAddCompCryst()",10)
   WXCrystValidateAllUserInput();
   PowderPatternDiffraction * diffData=new PowderPatternDiffraction;
   int choice;
   Crystal *cryst=dynamic_cast<Crystal*>
      ( WXDialogChooseFromRegistry(gCrystalRegistry,(wxWindow*)this,
         "Choose a Crystal Structure:",choice));
   if(0==cryst) {delete diffData;return;}
   VFN_DEBUG_MESSAGE("WXPowderPattern::OnMenuAddCompCryst()",10)
   diffData->SetCrystal(*cryst);
   VFN_DEBUG_MESSAGE("WXPowderPattern::OnMenuAddCompCryst()",10)
   mpPowderPattern->AddPowderPatternComponent(*diffData);
   VFN_DEBUG_MESSAGE("WXPowderPattern::OnMenuAddCompCryst()",10)
   if(diffData->GetRadiation().GetWavelengthType()==WAVELENGTH_TOF)
   {
   VFN_DEBUG_MESSAGE("WXPowderPattern::OnMenuAddCompCryst()",10)
      //wxCommandEvent event(wxEVT_COMMAND_MENU_SELECTED,ID_POWDERDIFF_PROFILE_DEPV);
      //wxPostEvent(diffData->WXGet(),event);
      ReflectionProfileDoubleExponentialPseudoVoigt *p=
         new ReflectionProfileDoubleExponentialPseudoVoigt
            (diffData->GetCrystal());
      diffData->SetProfile(p);
   }
   VFN_DEBUG_MESSAGE("WXPowderPattern::OnMenuAddCompCryst()",10)
   if(mpGraph!=0) mpPowderPattern->Prepare();//else this will be done when opening the graph
   this->CrystUpdate();
   VFN_DEBUG_EXIT("WXPowderPattern::OnMenuAddCompCryst()",10)
}

void WXPowderPattern::OnMenuShowGraph(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXPowderPattern::OnMenuShowGraph()"<<mpGraph,6)
   if(mpGraph!=0) return;
   if(mpPowderPattern->GetNbPoint()<=0) return;
   WXCrystValidateAllUserInput();
   mpPowderPattern->Prepare();
   wxFrame *frame= new wxFrame(this,-1, wxString::FromAscii(mpPowderPattern->GetName().c_str()),
                               wxDefaultPosition,wxSize(500,300));
   mpGraph = new WXPowderPatternGraph(frame,this);
   
   wxSizer *ps=new wxBoxSizer(wxHORIZONTAL);
   ps->Add(mpGraph,1,wxEXPAND);
   frame->SetSizer(ps);
   frame->SetAutoLayout(true);
   
   frame->CreateStatusBar(2);
   frame->Show(true);
   this->CrystUpdate(true);
   //frame->SetStatusText("");
}

void WXPowderPattern::OnMenuSaveText(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXPowderPattern::OnMenuSaveText()",6)
   WXCrystValidateAllUserInput();
   wxFileDialog save(this,_T("Choose a file"),_T(""),_T(""),_T("*.txt"),wxSAVE | wxOVERWRITE_PROMPT);
   if(save.ShowModal() != wxID_OK) return;
   
   ofstream out(save.GetPath().ToAscii());
   if(!out) return;//:TODO:
   mpPowderPattern->PrintObsCalcData(out);
   out.close();
}

void WXPowderPattern::OnMenuSimulate(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_ENTRY("WXPowderPattern::OnMenuSimulate()",6)
   WXCrystValidateAllUserInput();
   double min=0.,max=120.;
   long nbPoints=6000;
   {
      wxTextEntryDialog dialog(this,_T("2Theta Min"),
                              _T("Enter minimum 2Theta (degrees)"),_T("5"),wxOK | wxCANCEL);
      if(wxID_OK!=dialog.ShowModal())
      {
         VFN_DEBUG_EXIT("WXPowderPattern::OnMenuSimulate():Cancelled",6)
         return;
      }
      dialog.GetValue().ToDouble(&min);
      if(min<1) min=1.0;
   }
   {
      wxTextEntryDialog dialog(this,_T("2Theta Max"),
                              _T("Enter maximum 2Theta (degrees)"),_T("100"),wxOK | wxCANCEL);
      if(wxID_OK!=dialog.ShowModal())
      {
         VFN_DEBUG_EXIT("WXPowderPattern::OnMenuSimulate():Cancelled",6)
         return;
      }
      dialog.GetValue().ToDouble(&max);
   }
   {
      wxTextEntryDialog dialog(this,_T("Number of points"),
                              _T("Enter the number of points"),_T("1000"),wxOK | wxCANCEL);
      if(wxID_OK!=dialog.ShowModal())
      {
         VFN_DEBUG_EXIT("WXPowderPattern::OnMenuSimulate():Cancelled",6)
         return;
      }
      dialog.GetValue().ToLong(&nbPoints);
   }
   CrystVector_REAL newObs(nbPoints);
   mpPowderPattern->SetPowderPatternPar(min*DEG2RAD,(max-min)/(nbPoints-1)*DEG2RAD,nbPoints);
   newObs=1;//we must not have 0 in case a scale factor is fitted...
   newObs(0)=0.01;//Avoid having the same value for ALL points for scaling the graph.
   if(mpPowderPattern->GetNbPowderPatternComponent()>0)
   {
      // Use the calculated pattern, for indexing simulation
      newObs=mpPowderPattern->GetPowderPatternCalc();
      // Add some noise !
      for(long i=0;i<newObs.numElements();++i)
         newObs(i) += sqrt(newObs(i))*(2*rand()/(REAL)RAND_MAX-1);
   }
   mpPowderPattern->SetPowderPatternObs(newObs);
   VFN_DEBUG_EXIT("WXPowderPattern::OnMenuSimulate()",6)
}

void WXPowderPattern::OnMenuImportPattern(wxCommandEvent &event)
{
   VFN_DEBUG_MESSAGE("WXPowderPattern::OnMenuImportPattern()",6)
   wxFileDialog open(this,_T("Choose a file"),_T(""),_T(""),_T("*.*"),wxOPEN | wxFILE_MUST_EXIST);
   if(open.ShowModal() != wxID_OK) return;
   if(event.GetId()==(long)ID_POWDER_MENU_IMPORT_FULLPROF)
      mpPowderPattern->ImportPowderPatternFullprof(string(open.GetPath().ToAscii()));
   if(event.GetId()==(long)ID_POWDER_MENU_IMPORT_PSI_DMC)
      mpPowderPattern->ImportPowderPatternPSI_DMC(string(open.GetPath().ToAscii()));
   if(event.GetId()==(long)ID_POWDER_MENU_IMPORT_ILL_D1A5)
      mpPowderPattern->ImportPowderPatternILL_D1A5(string(open.GetPath().ToAscii()));
   if(event.GetId()==(long)ID_POWDER_MENU_IMPORT_XDD)
      mpPowderPattern->ImportPowderPatternXdd(string(open.GetPath().ToAscii()));
   if(event.GetId()==(long)ID_POWDER_MENU_IMPORT_CPI)
      mpPowderPattern->ImportPowderPatternSietronicsCPI(string(open.GetPath().ToAscii()));
   if(event.GetId()==(long)ID_POWDER_MENU_IMPORT_FULLPROF4)
      mpPowderPattern->ImportPowderPatternFullprof4(string(open.GetPath().ToAscii()));
   if(event.GetId()==(long)ID_POWDER_MENU_IMPORT_MULTIDETECTORLLBG42)
      mpPowderPattern->ImportPowderPatternMultiDetectorLLBG42(string(open.GetPath().ToAscii()));
   if(event.GetId()==(long)ID_POWDER_MENU_IMPORT_2THETAOBSSIGMA)
      mpPowderPattern->ImportPowderPattern2ThetaObsSigma(string(open.GetPath().ToAscii()));
   if(event.GetId()==(long)ID_POWDER_MENU_IMPORT_2THETAOBS)
      mpPowderPattern->ImportPowderPattern2ThetaObs(string(open.GetPath().ToAscii()));
   if(event.GetId()==(long)ID_POWDER_MENU_IMPORT_TOFISISXYSIGMA)
      mpPowderPattern->ImportPowderPatternTOF_ISIS_XYSigma(string(open.GetPath().ToAscii()));
   if(event.GetId()==(long)ID_POWDER_MENU_IMPORT_GSAS)
      mpPowderPattern->ImportPowderPatternGSAS(string(open.GetPath().ToAscii()));
   if(event.GetId()==(long)ID_POWDER_MENU_IMPORT_CIF)
   {
      ifstream fin (open.GetPath().ToAscii());
      if(!fin)
      {
         throw ObjCrystException("WXPowderPattern::OnMenuImportPattern(): Error opening file for input:"+string(open.GetPath().ToAscii()));
      }
      ObjCryst::CIF cif(fin,true,true);
      mpPowderPattern->ImportPowderPatternCIF(cif);
   }
   bool val;
   wxConfigBase::Get()->Read(_T("PowderPattern/BOOL/Automatically open powder pattern graph"), &val);
   if(val)
   {
      wxCommandEvent event(wxEVT_COMMAND_MENU_SELECTED,ID_POWDER_MENU_GRAPH);
      wxPostEvent(this,event);
   }
}

void WXPowderPattern::OnMenuFitScaleForR(wxCommandEvent & WXUNUSED(event))
{
   if(0==mpGraph) return;
   WXCrystValidateAllUserInput();
   mpPowderPattern->FitScaleFactorForR();//FitScaleFactorForIntegratedR
   this->CrystUpdate(true);
}

void WXPowderPattern::OnMenuFitScaleForRw(wxCommandEvent & WXUNUSED(event))
{
   if(0==mpGraph) return;
   WXCrystValidateAllUserInput();
   mpPowderPattern->FitScaleFactorForRw();//FitScaleFactorForIntegratedRw
   this->CrystUpdate(true);
}


void WXPowderPattern::OnMenuSetWavelength(wxCommandEvent & event)
{
   WXCrystValidateAllUserInput();
   // this looks stupid. In fact, if a user changed the wavelength in the
   // corresponding field, this is (unfortunately) not applied to the
   // components automagically. So we need this function to do the job...
   if(event.GetId()== ID_POWDER_MENU_WAVELENGTH_XRAY)
      mpPowderPattern->SetRadiationType(RAD_XRAY);
   if(event.GetId()== ID_POWDER_MENU_WAVELENGTH_NEUTRON)
      mpPowderPattern->SetRadiationType(RAD_NEUTRON);
   if(event.GetId()== ID_POWDER_MENU_WAVELENGTH_NEUTRON_TOF)
   {
      mpPowderPattern->SetRadiationType(RAD_NEUTRON);
      mpPowderPattern->GetRadiation().SetWavelengthType(WAVELENGTH_TOF);
   }
   if(event.GetId()== ID_POWDER_MENU_WAVELENGTH_SET)
   {
      double lambda;
      wxTextEntryDialog dialog(this,_T("new Wavelength)"),
                              _T("Enter new Wavelength (Angstroems)"),_T("1"),wxOK | wxCANCEL);
      if(wxID_OK!=dialog.ShowModal())
      {
         VFN_DEBUG_EXIT("WXPowderPattern))OnMenuSetWavelength())Monochromatic)Cancelled",6)
         return;
      }
      dialog.GetValue().ToDouble(&lambda);
      mpPowderPattern->SetWavelength(lambda);
   }
   if(event.GetId()== ID_POWDER_MENU_WAVELENGTH_SET_AG)
      mpPowderPattern->SetWavelength("Ag");
   if(event.GetId()== ID_POWDER_MENU_WAVELENGTH_SET_MO)
      mpPowderPattern->SetWavelength("Mo");
   if(event.GetId()== ID_POWDER_MENU_WAVELENGTH_SET_CU)
      mpPowderPattern->SetWavelength("Cu");
   if(event.GetId()== ID_POWDER_MENU_WAVELENGTH_SET_FE)
      mpPowderPattern->SetWavelength("Fe");
   if(event.GetId()== ID_POWDER_MENU_WAVELENGTH_SET_CO)
      mpPowderPattern->SetWavelength("Co");
   if(event.GetId()== ID_POWDER_MENU_WAVELENGTH_SET_CR)
      mpPowderPattern->SetWavelength("Cr");
   if(event.GetId()== ID_POWDER_MENU_WAVELENGTH_SET_AGA1)
      mpPowderPattern->SetWavelength("AgA1");
   if(event.GetId()== ID_POWDER_MENU_WAVELENGTH_SET_MOA1)
      mpPowderPattern->SetWavelength("MoA1");
   if(event.GetId()== ID_POWDER_MENU_WAVELENGTH_SET_CUA1)
      mpPowderPattern->SetWavelength("CuA1");
   if(event.GetId()== ID_POWDER_MENU_WAVELENGTH_SET_FEA1)
      mpPowderPattern->SetWavelength("FeA1");
   if(event.GetId()== ID_POWDER_MENU_WAVELENGTH_SET_COA1)
      mpPowderPattern->SetWavelength("CoA1");
   if(event.GetId()== ID_POWDER_MENU_WAVELENGTH_SET_CRA1)
      mpPowderPattern->SetWavelength("CrA1");
   this->CrystUpdate(true);
}

void WXPowderPattern::OnMenuAddExclude(wxCommandEvent & WXUNUSED(event))
{
   WXCrystValidateAllUserInput();
   double min,max;
   //min
   {
      wxString txt=_T("Enter Min 2theta to exclude (degrees):");
      if(mpPowderPattern->GetRadiation().GetWavelengthType()==WAVELENGTH_TOF)
         txt=_T("Enter Min 2theta to exclude (microseconds):");
      wxTextEntryDialog dialog(this,_T("Min"),txt,_T("0"),wxOK | wxCANCEL);
      if(wxID_OK!=dialog.ShowModal())
      {
         VFN_DEBUG_EXIT("WXPowderPattern::OnMenuAddExclude():Cancelled",6)
         return;
      }
      dialog.GetValue().ToDouble(&min);
   }
   //max
   {
      wxString txt=_T("Enter Max 2theta to exclude (degrees):");
      if(mpPowderPattern->GetRadiation().GetWavelengthType()==WAVELENGTH_TOF)
         txt=_T("Enter Max 2theta to exclude (microseconds):");
      wxTextEntryDialog dialog(this,_T("Max"),txt,_T("5"),wxOK | wxCANCEL);
      if(wxID_OK!=dialog.ShowModal())
      {
         VFN_DEBUG_EXIT("WXPowderPattern::OnMenuAddExclude():Cancelled",6)
         return;
      }
      dialog.GetValue().ToDouble(&max);
   }
   if(max<min)
   {
      VFN_DEBUG_EXIT("WXPowderPattern::OnMenuAddExclude():Stupid user.",6)
      return;
   }
   if(mpPowderPattern->GetRadiation().GetWavelengthType()==WAVELENGTH_TOF)
      mpPowderPattern->AddExcludedRegion(min,max);
   else mpPowderPattern->AddExcludedRegion(min*DEG2RAD,max*DEG2RAD);
}

void WXPowderPattern::OnMenuLeBail(wxCommandEvent& event)
{
   #if 0
   wxFrame *pFrame=new wxFrame(this,-1,_T("Profile Fitting"));
   WXProfileFitting *pFit;
   pFit=new WXProfileFitting(pFrame,&(this->GetPowderPattern()));
   pFrame->Show(true);
   #else
   cout<<"Beginning refinement"<<endl;
   LSQNumObj lsq;
   lsq.SetRefinedObj(this->GetPowderPattern());
   lsq.PrepareRefParList(true);
   lsq.SetParIsFixed(gpRefParTypeObjCryst,true);
   lsq.SetParIsFixed(gpRefParTypeScatt,false);
   lsq.SetParIsFixed(gpRefParTypeScattDataScale,false);
   //lsq.SetParIsUsed(gpRefParTypeScattDataProfile,true);
   //lsq.SetParIsUsed(gpRefParTypeScattDataCorrPos,true);
   //lsq.SetParIsUsed(gpRefParTypeScattDataBackground,true);
   //lsq.SetParIsUsed(gpRefParTypeUnitCell,true);
   try {lsq.Refine(10,true,false);}
   catch(const ObjCrystException &except){};
   cout<<"Finishing refinement"<<endl;
   this->GetPowderPattern().UpdateDisplay();
   #endif
}

void WXPowderPattern::OnMenuExport(wxCommandEvent &event)
{
   WXCrystValidateAllUserInput();
   wxFileDialog save(this,_T("Choose a .pcr file"),_T(""),_T(""),_T("*.pcr"),wxSAVE | wxOVERWRITE_PROMPT);
   if(save.ShowModal() != wxID_OK) return;
   
   wxString path,name,ext;
   wxFileName::SplitPath(save.GetPath(), &path, &name, &ext, wxPATH_NATIVE);
   wxString mes;
   wxString pcr=path+wxFileName::GetPathSeparator()+name+_T(".pcr");
   wxString dat=path+wxFileName::GetPathSeparator()+name+_T(".dat");
   mes.Printf(_T("This will create the files:\n   %s\n   %s"),pcr.c_str(),dat.c_str());
   wxMessageDialog mesd(this,mes,_T("Files"),wxOK|wxCANCEL|wxICON_INFORMATION);
   if(mesd.ShowModal()==wxID_OK)
      mpPowderPattern->ExportFullprof(string((path+wxFileName::GetPathSeparator()+name).ToAscii()));
}

void WXPowderPattern::NotifyDeleteGraph() {mpGraph=0;}
const PowderPattern& WXPowderPattern::GetPowderPattern()const
{ return *mpPowderPattern;}

PowderPattern& WXPowderPattern::GetPowderPattern()
{ return *mpPowderPattern;}

void WXPowderPattern::UpdateUI(const bool lock)
{
   if(lock)mMutex.Lock();
   if(mpGraph!=0)
   {
      mpGraph->GetParent()->SetLabel(wxString::FromAscii(mpPowderPattern->GetName().c_str()));
   }
   this->WXRefinableObj::UpdateUI(false);
   if(lock)mMutex.Unlock();
}
////////////////////////////////////////////////////////////////////////
//
//    WXPowderPatternGraph
//
////////////////////////////////////////////////////////////////////////
static const long ID_POWDERGRAPH_MENU_UPDATE=               WXCRYST_ID(); 
static const long ID_POWDERGRAPH_MENU_TOGGLELABEL=          WXCRYST_ID(); 
static const long ID_POWDERGRAPH_MENU_TOGGPEAK=             WXCRYST_ID(); 
static const long ID_POWDERGRAPH_MENU_FINDPEAKS=            WXCRYST_ID(); 
static const long ID_POWDERGRAPH_MENU_LOADPEAKS=            WXCRYST_ID(); 
static const long ID_POWDERGRAPH_MENU_SAVEPEAKS=            WXCRYST_ID(); 
static const long ID_POWDERGRAPH_MENU_ADDPEAK=              WXCRYST_ID(); 
static const long ID_POWDERGRAPH_MENU_REMOVEPEAK=           WXCRYST_ID(); 
static const long ID_POWDERGRAPH_MENU_INDEX=                WXCRYST_ID(); 
static const long ID_POWDERGRAPH_MENU_XSCALE_DATA=          WXCRYST_ID(); 
static const long ID_POWDERGRAPH_MENU_XSCALE_D=             WXCRYST_ID(); 
static const long ID_POWDERGRAPH_MENU_XSCALE_2PID=          WXCRYST_ID(); 
static const long ID_POWDERGRAPH_MENU_YSCALE_LINEAR=        WXCRYST_ID(); 
static const long ID_POWDERGRAPH_MENU_YSCALE_SQRT=          WXCRYST_ID(); 
static const long ID_POWDERGRAPH_MENU_YSCALE_LOG10=         WXCRYST_ID(); 
static const long ID_POWDERGRAPH_MENU_LEBAIL=               WXCRYST_ID(); 

BEGIN_EVENT_TABLE(WXPowderPatternGraph, wxWindow)
   EVT_PAINT(                                   WXPowderPatternGraph::OnPaint)
   EVT_MOUSE_EVENTS(                            WXPowderPatternGraph::OnMouse)
   EVT_MENU(ID_POWDERGRAPH_MENU_UPDATE,         WXPowderPatternGraph::OnUpdate)
   EVT_MENU(ID_POWDERGRAPH_MENU_TOGGLELABEL,    WXPowderPatternGraph::OnToggleLabel)
   EVT_MENU(ID_POWDERGRAPH_MENU_TOGGPEAK,       WXPowderPatternGraph::OnToggleLabel)
   EVT_MENU(ID_POWDERGRAPH_MENU_FINDPEAKS,      WXPowderPatternGraph::OnFindPeaks)
   EVT_MENU(ID_POWDERGRAPH_MENU_LOADPEAKS,      WXPowderPatternGraph::OnLoadPeaks)
   EVT_MENU(ID_POWDERGRAPH_MENU_SAVEPEAKS,      WXPowderPatternGraph::OnSavePeaks)
   EVT_MENU(ID_POWDERGRAPH_MENU_ADDPEAK,        WXPowderPatternGraph::OnChangePeak)
   EVT_MENU(ID_POWDERGRAPH_MENU_REMOVEPEAK,     WXPowderPatternGraph::OnChangePeak)
   EVT_MENU(ID_POWDERGRAPH_MENU_INDEX,          WXPowderPatternGraph::OnIndex)
   EVT_MENU(ID_POWDERGRAPH_MENU_XSCALE_DATA,    WXPowderPatternGraph::OnChangeScale)
   EVT_MENU(ID_POWDERGRAPH_MENU_XSCALE_D,       WXPowderPatternGraph::OnChangeScale)
   EVT_MENU(ID_POWDERGRAPH_MENU_XSCALE_2PID,    WXPowderPatternGraph::OnChangeScale)
   EVT_MENU(ID_POWDERGRAPH_MENU_YSCALE_LINEAR,  WXPowderPatternGraph::OnChangeScale)
   EVT_MENU(ID_POWDERGRAPH_MENU_YSCALE_SQRT,    WXPowderPatternGraph::OnChangeScale)
   EVT_MENU(ID_POWDERGRAPH_MENU_YSCALE_LOG10,   WXPowderPatternGraph::OnChangeScale)
   EVT_MENU(ID_POWDERGRAPH_MENU_LEBAIL,         WXPowderPatternGraph::OnLeBail)
   EVT_UPDATE_UI(ID_POWDER_GRAPH_NEW_PATTERN,   WXPowderPatternGraph::OnRedrawNewPattern)
   EVT_CHAR(                                    WXPowderPatternGraph::OnKeyDown)
   EVT_MOUSEWHEEL(                              WXPowderPatternGraph::OnMouseWheel)
   EVT_SIZE(                                    WXPowderPatternGraph::OnSize)
END_EVENT_TABLE()

WXPowderPatternGraph::WXPowderPatternGraph(wxFrame *frame, WXPowderPattern* parent):
wxWindow(frame,-1,wxPoint(-1,-1),wxSize(-1,-1)),
mpPattern(parent),mMargin(20),mDiffPercentShift(.20),
mMaxIntensity(-1),mMinIntensity(-1),mMinX(-1),mMaxX(-1),
mpParentFrame(frame),
mIsDragging(false),mDisplayLabel(true),mDisplayPeak(true)
{
   mpPopUpMenu=new wxMenu(_T("Powder Pattern"));
   mpPopUpMenu->Append(ID_POWDERGRAPH_MENU_UPDATE, _T("&Update"));
   mpPopUpMenu->Append(ID_POWDERGRAPH_MENU_TOGGLELABEL, _T("&Hide labels"));
   #if 1
   mpPopUpMenu->AppendSeparator();
   mpPopUpMenu->Append(ID_POWDERGRAPH_MENU_FINDPEAKS, _T("&Find peaks"));
   mpPopUpMenu->Append(ID_POWDERGRAPH_MENU_LOADPEAKS, _T("&Load peaks"));
   mpPopUpMenu->Append(ID_POWDERGRAPH_MENU_SAVEPEAKS, _T("&Save peaks"));
   mpPopUpMenu->Append(ID_POWDERGRAPH_MENU_INDEX, _T("&Index !"));
   mpPopUpMenu->Append(ID_POWDERGRAPH_MENU_TOGGPEAK, _T("&Hide peaks"));
   mpPopUpMenu->Append(ID_POWDERGRAPH_MENU_ADDPEAK, _T("&Add peak"));
   mpPopUpMenu->Append(ID_POWDERGRAPH_MENU_REMOVEPEAK, _T("&Remove peak"));
   mpPopUpMenu->Enable(ID_POWDERGRAPH_MENU_SAVEPEAKS, FALSE);
   mpPopUpMenu->Enable(ID_POWDERGRAPH_MENU_TOGGPEAK, FALSE);
   mpPopUpMenu->Enable(ID_POWDERGRAPH_MENU_ADDPEAK, TRUE);
   mpPopUpMenu->Enable(ID_POWDERGRAPH_MENU_REMOVEPEAK, FALSE);
   mpPopUpMenu->Enable(ID_POWDERGRAPH_MENU_INDEX, FALSE);
   #endif
   mpPopUpMenu->AppendSeparator();
   mpPopUpMenu->Append(ID_POWDERGRAPH_MENU_LEBAIL, _T("Fit &Profile + Le Bail extraction"));
   mpPopUpMenu->AppendSeparator();
   mpPopUpMenu->Append(ID_POWDERGRAPH_MENU_XSCALE_DATA, _T("&X scale: 2theta/TOF"));
   mpPopUpMenu->Append(ID_POWDERGRAPH_MENU_XSCALE_D, _T("&X scale: Q=1/d"));
   mpPopUpMenu->Append(ID_POWDERGRAPH_MENU_XSCALE_2PID, _T("&X scale: Q=2pi/d"));
   mpPopUpMenu->Append(ID_POWDERGRAPH_MENU_YSCALE_LINEAR, _T("&Y scale: I"));
   mpPopUpMenu->Append(ID_POWDERGRAPH_MENU_YSCALE_SQRT, _T("&Y scale: sqrt(I)"));
   mpPopUpMenu->Append(ID_POWDERGRAPH_MENU_YSCALE_LOG10, _T("&Y scale: log10(I)"));
   if(!wxConfigBase::Get()->HasEntry(_T("PowderPattern/BOOL/Default-display reflection indices")))
      wxConfigBase::Get()->Write(_T("PowderPattern/BOOL/Default-display reflection indices"), mDisplayLabel);
   else
   {
      wxConfigBase::Get()->Read(_T("PowderPattern/BOOL/Default-display reflection indices"), &mDisplayLabel);
      if(mDisplayLabel) mpPopUpMenu->SetLabel(ID_POWDERGRAPH_MENU_TOGGLELABEL, _T("Hide labels"));
      else mpPopUpMenu->SetLabel(ID_POWDERGRAPH_MENU_TOGGLELABEL, _T("Show labels"));
   }
   
   // Scale used to display graph x coordinates : 0 - experimental ; 1 - 1/d ; 2 - 2pi/d
   if(!wxConfigBase::Get()->HasEntry(_T("PowderPattern/LONG/graph x scale")))
      wxConfigBase::Get()->Write(_T("PowderPattern/LONG/graph x scale"), 0);
   
   // Scale used to display graph y coordinates : 0 - linear ; 1 - square root ; 2 - log10
   if(!wxConfigBase::Get()->HasEntry(_T("PowderPattern/LONG/graph y scale")))
      wxConfigBase::Get()->Write(_T("PowderPattern/LONG/graph y scale"), 0);
   
   wxConfigBase::Get()->Read(_T("PowderPattern/LONG/graph x scale"), &mXScale);
   wxConfigBase::Get()->Read(_T("PowderPattern/LONG/graph y scale"), &mYScale);

   mpPopUpMenu->Enable(ID_POWDERGRAPH_MENU_XSCALE_DATA, TRUE);
   mpPopUpMenu->Enable(ID_POWDERGRAPH_MENU_XSCALE_D, TRUE);
   mpPopUpMenu->Enable(ID_POWDERGRAPH_MENU_XSCALE_2PID, TRUE);
   mpPopUpMenu->Enable(ID_POWDERGRAPH_MENU_YSCALE_LINEAR, TRUE);
   mpPopUpMenu->Enable(ID_POWDERGRAPH_MENU_YSCALE_SQRT, TRUE);
   mpPopUpMenu->Enable(ID_POWDERGRAPH_MENU_YSCALE_LOG10, TRUE);
   
   if(mXScale==0)  mpPopUpMenu->Enable(ID_POWDERGRAPH_MENU_XSCALE_DATA, FALSE);
   if(mXScale==1)     mpPopUpMenu->Enable(ID_POWDERGRAPH_MENU_XSCALE_D, FALSE);
   if(mXScale==2)  mpPopUpMenu->Enable(ID_POWDERGRAPH_MENU_XSCALE_2PID, FALSE);
   if(mYScale==0)mpPopUpMenu->Enable(ID_POWDERGRAPH_MENU_YSCALE_LINEAR, FALSE);
   if(mYScale==1)  mpPopUpMenu->Enable(ID_POWDERGRAPH_MENU_YSCALE_SQRT, FALSE);
   if(mYScale==2) mpPopUpMenu->Enable(ID_POWDERGRAPH_MENU_YSCALE_LOG10, FALSE);

   mpPattern->CrystUpdate(true);
}

WXPowderPatternGraph::~WXPowderPatternGraph()
{
   mpPattern->NotifyDeleteGraph();
}

void WXPowderPatternGraph::OnPaint(wxPaintEvent& WXUNUSED(event))
{
   if((mObs.numElements()<=0)||(mCalc.numElements()<=0)) return;
   VFN_DEBUG_MESSAGE("WXPowderPatternGraph:OnPaint()",5)
   unsigned int ct=0;
   while(wxMUTEX_NO_ERROR!=mMutex.TryLock())
   {
     if(++ct>10) return;
     wxMilliSleep(50);
   }
   
   wxBufferedPaintDC dc(this);
   PrepareDC(dc);
   mpParentFrame->PrepareDC(dc);
   
   dc.BeginDrawing();
   
   dc.DestroyClippingRegion();
   dc.SetBackground(wxBrush(_T("white"), wxSOLID));
   dc.Clear();

   wxString fontInfo;
   #ifdef __WIN32__
   dc.SetFont(*wxNORMAL_FONT);
   #else
   dc.SetFont(*wxSMALL_FONT);
   #endif

   long nbPoint=mX.numElements();
   
   // Get Window Size
   wxCoord width,height;
   this->GetSize(&width, &height);
   //const int margin=mMargin;
   //width -= margin;
   //height -= margin;
   VFN_DEBUG_MESSAGE("WXPowderPatternGraph:OnPaint():1",5)

   //Check pattern is not being updated
   VFN_DEBUG_MESSAGE("WXPowderPatternGraph:OnPaint():2:"<<mObs.numElements(),5)

   VFN_DEBUG_MESSAGE("WXPowderPatternGraph:OnPaint():3:min="
                     <<mMinX<<", max="<<mMaxX<<", width="
                     <<width<<",margin="<<mMargin,5)
   // Draw sigma bars
   {
      dc.SetPen(* wxLIGHT_GREY_PEN);
      wxCoord x,y1,y2;
      for(long i=0;i<nbPoint;i++)
      {
         if((mX(i)>mMinX)&&(mX(i)<mMaxX))
         {
            x=this->Point2ScreenX(i);
            y1=this->Data2ScreenY(mObs(i)-mSigma(i)/2.);
            y2=this->Data2ScreenY(mObs(i)+mSigma(i)/2.);

            dc.DrawLine(x,y1,x,y2);
         }
      }
   }
   // Draw Axis (sort of)
   {
      wxCoord tmpW,tmpH;
      dc.SetPen(* wxBLACK_PEN);
      dc.DrawLine(mMargin*3,height-mMargin,mMargin*3,mMargin);
      dc.DrawLine(mMargin*3,height-mMargin,width,height-mMargin);
      const int nbTick=10;//approximate
      wxCoord xc,yc;
      //Y axis
      {
         xc=(wxCoord)mMargin*3;
         REAL miny=mMinIntensity,maxy=mMaxIntensity;
         if(mYScale==1) {miny=sqrt(miny) ;maxy=sqrt(maxy);}
         if(mYScale==2) {miny=log10(miny);maxy=log10(maxy);}
         REAL yStep=pow((float)10,(float)floor(log10((maxy-miny)/nbTick)));
         yStep *= floor((maxy-miny)/yStep/nbTick);
         for(REAL ys=yStep*ceil(miny/yStep);ys<maxy;ys+=yStep)
         {
            REAL y=ys;
            if(mYScale==1) {y=ys*ys;}
            if(mYScale==2) {y=pow((float)10,(float)ys);}
            yc=this->Data2ScreenY(y);
            dc.DrawLine(xc-3,yc,xc+3,yc);
            fontInfo.Printf(_T("%g"),y);
            dc.GetTextExtent(fontInfo, &tmpW, &tmpH);
            dc.DrawText(fontInfo,xc-tmpW-3,yc-tmpH/2);
         }
      }
      //X axis
      {
         yc=(wxCoord)(height-mMargin);
         
         REAL minx=mMinX,maxx=mMaxX;
         float mind,maxd;// 1/d
         if(mpPattern->GetPowderPattern().GetRadiation().GetWavelengthType()!=WAVELENGTH_TOF)
         {
            mind=2*mpPattern->GetPowderPattern().X2STOL(minx*DEG2RAD);
            maxd=2*mpPattern->GetPowderPattern().X2STOL(maxx*DEG2RAD);
         }
         else
         {
            mind=2*mpPattern->GetPowderPattern().X2STOL(minx);
            maxd=2*mpPattern->GetPowderPattern().X2STOL(maxx);
         }
         if(mXScale==1) {minx=mind;       maxx=maxd;}
         if(mXScale==2) {minx=2*M_PI*mind;maxx=2*M_PI*maxd;}
         
         REAL xStep=pow((float)10,(float)floor(log10((maxx-minx)/nbTick)));
         xStep *= floor((maxx-minx)/xStep/nbTick);
         for(REAL xs=xStep*ceil(minx/xStep);xs<maxx;xs+=xStep)
         {
            REAL x=xs;
            if(mXScale==1) {x=mpPattern->GetPowderPattern().STOL2X(xs/2);}
            if(mXScale==2) {x=mpPattern->GetPowderPattern().STOL2X(xs/(4*M_PI));}
            if(mXScale==0) xc=this->Data2ScreenX(x);
            else
            {
               if(mpPattern->GetPowderPattern().GetRadiation().GetWavelengthType()!=WAVELENGTH_TOF)
                  xc=this->Data2ScreenX(RAD2DEG*x);
               else xc=this->Data2ScreenX(x);
            }
            dc.DrawLine(xc,yc-3,xc,yc+3);
            fontInfo.Printf(_T("%g"),xs);
            dc.GetTextExtent(fontInfo, &tmpW, &tmpH);
            dc.DrawText(fontInfo,xc-tmpW/2,yc+6);
         }
      }
   }
   // Draw cumulated Chi^2, scaled
   {
      dc.SetPen(* wxGREY_PEN);
      wxCoord x1,y1,x2,y2;
      x2=this->Point2ScreenX(0);
      const REAL s=(mMaxIntensity-mMinIntensity)/mChi2Cumul(mpPattern->GetPowderPattern().GetNbPointUsed()-1);
      y2=this->Data2ScreenY(mMinIntensity+mChi2Cumul(0)*s);
      for(unsigned long i=0;i<mpPattern->GetPowderPattern().GetNbPointUsed();i++)
      {
         if((mX(i)>mMinX)&&(mX(i)<mMaxX))
         {
            x1=x2;
            y1=y2;
            x2=this->Point2ScreenX(i);
            y2=this->Data2ScreenY(mMinIntensity+mChi2Cumul(i)*s);
            dc.DrawLine(x1,y1,x2,y2);
         }
      }
   }
   // Draw observed pattern
   VFN_DEBUG_MESSAGE("WXPowderPatternGraph:OnPaint():4:",5)
   {
      dc.SetPen(* wxCYAN_PEN);
      wxCoord x1,y1,x2,y2;
      x2=this->Point2ScreenX(0);
      y2=this->Data2ScreenY(mObs(0));
      for(long i=0;i<nbPoint;i++)
      {
         if((mX(i)>mMinX)&&(mX(i)<mMaxX))
         {
            x1=x2;
            y1=y2;
            x2=this->Point2ScreenX(i);
            y2=this->Data2ScreenY(mObs(i));
            dc.DrawLine(x1,y1,x2,y2);
         }
      }
   }

   // Draw calculated pattern
   VFN_DEBUG_MESSAGE("WXPowderPatternGraph:OnPaint():5:",5)
   {
      dc.SetPen(* wxRED_PEN);
      wxCoord x1,y1,x2,y2;
      x2=this->Point2ScreenX(0);
      y2=this->Data2ScreenY(mCalc(0));
      for(long i=0;i<nbPoint;i++)
      {
         if((mX(i)>mMinX)&&(mX(i)<mMaxX))
         {
            x1=x2;
            y1=y2;
            x2=this->Point2ScreenX(i);
            y2=this->Data2ScreenY(mCalc(i));
            dc.DrawLine(x1,y1,x2,y2);
         }
      }
   }
   
   // Display labels ?
   list<list<pair<const REAL ,const string > > > vLabel;
   if(true==mDisplayLabel) // "vLabel=mvLabelList;" does not work (gcc 4.1.1)
      for(list<list<pair<const REAL ,const string > > >::const_iterator 
          comp=mvLabelList.begin();comp!=mvLabelList.end();++comp) vLabel.push_back(*comp);
   // Show peaks ?
   if((true==mDisplayPeak)&&(mPeakList.GetPeakList().size()>0))
   {
      list<pair<const REAL ,const string > > peakLabels;
      char buf[50];
      unsigned int ix=0;
      for(vector<PeakList::hkl>::const_iterator pos=mPeakList.GetPeakList().begin();pos!=mPeakList.GetPeakList().end();++pos)
      {
         const float x=mpPattern->GetPowderPattern().STOL2X(pos->dobs/2);
         if(pos->isSpurious)
         {
            if(mpPattern->GetPowderPattern().GetRadiation().GetWavelengthType()!=WAVELENGTH_TOF)
               sprintf(buf,"#%2u,x=%6.3f° d=%6.3fA SPURIOUS?",ix,x*RAD2DEG,1/pos->dobs);
            else
               sprintf(buf,"#%2u,x=%6.3f d=%6.3fA SPURIOUS?",ix,x        ,1/pos->dobs);
         }
         else
         {
            if(mpPattern->GetPowderPattern().GetRadiation().GetWavelengthType()!=WAVELENGTH_TOF)
               sprintf(buf,"#%2u,x=%6.2f° d=%6.3fA",ix,x*RAD2DEG,1/pos->dobs);
            else
               sprintf(buf,"#%2u,x=%6.2f d=%6.3fA",ix,x        ,1/pos->dobs);
         }
         ++ix;
         peakLabels.push_back(make_pair(x,buf));
      }
      vLabel.push_back(peakLabels);
      // Do we have a list of predicted HKL positions as well ?
      if(mPeakList.mvPredictedHKL.size()>0)
      {
         peakLabels.clear();
         for(list<PeakList::hkl>::const_iterator pos=mPeakList.mvPredictedHKL.begin();pos!=mPeakList.mvPredictedHKL.end();++pos)
         {
            const float dobs=sqrt(pos->d2calc);
            const float x=mpPattern->GetPowderPattern().STOL2X(dobs/2);
            sprintf(buf,"?(%2d %2d %2d)?",pos->h,pos->k,pos->l);
            peakLabels.push_back(make_pair(x,buf));
         }
         vLabel.push_back(peakLabels);
      }
   }
   // Draw labels
   VFN_DEBUG_MESSAGE("WXPowderPatternGraph:OnPaint():5:",5)
   if(vLabel.size()>0)
   {
      wxCoord x,y;
      wxCoord tmpW,tmpH;
      int loop=1;
      REAL yr;
      list<list<pair<const REAL ,const string > > >::const_iterator comp;
      list<pair<const REAL ,const string > >::const_iterator pos;
      unsigned int pen=0;
      for(comp=vLabel.begin();comp!=vLabel.end();++comp)
      {
         switch(pen++)
         {
            case 0: dc.SetPen(*wxBLACK_PEN);dc.SetTextForeground(*wxBLACK);break;
            case 1: dc.SetPen(*wxCYAN_PEN );dc.SetTextForeground(*wxCYAN );break;
            case 2: dc.SetPen(*wxGREEN_PEN);dc.SetTextForeground(*wxGREEN);break;
            case 3: dc.SetPen(*wxRED_PEN  );dc.SetTextForeground(*wxRED  );break;
            default:dc.SetPen(*wxGREY_PEN );dc.SetTextForeground(*wxLIGHT_GREY );break;
         }
         unsigned long ct=0;
         for(pos=comp->begin();pos!=comp->end();++pos)
         {
            REAL point=pos->first;
            if(mpPattern->GetPowderPattern().GetRadiation().GetWavelengthType()!=WAVELENGTH_TOF)
               point *= RAD2DEG;
            if((point>=mMinX)&&(point<=mMaxX))
            {
               if(++ct>200)
               {
                  cout <<"Too many labels (>100): displaying only first 100 and ticking 100 more..."<<endl;
                  break;
               }
               x=this->Data2ScreenX(point);
               const REAL pixel=mpPattern->GetPowderPattern().X2Pixel(pos->first);
               if(mCalc((long)pixel)>mObs((long)pixel)) yr=mCalc((long)pixel);
               else yr=mObs((long)pixel);
               y=this->Data2ScreenY(yr);
               
               dc.DrawLine(x,y-5,x,y-10);
               if(ct<100)
               {
                  fontInfo.Printf(wxString::FromAscii(pos->second.c_str()));
                  dc.GetTextExtent(fontInfo, &tmpW, &tmpH);
                  dc.DrawText(fontInfo,x-tmpW/2,y-tmpH*(loop++)-10);
               }
               if(loop==5) loop=1;
            }
         }
      }
   }
   
   dc.EndDrawing();
   mMutex.Unlock();
   VFN_DEBUG_MESSAGE("WXPowderPatternGraph:OnPaint():End",5)
}

void WXPowderPatternGraph::OnMouse(wxMouseEvent &event)
{
   if(event.Leaving()) return;// wxMSW2.4 bug ?
   if(wxMUTEX_NO_ERROR!=mMutex.TryLock())
   {
      mIsDragging=false;
      return;
   }
       VFN_DEBUG_MESSAGE("WXPowderPatternGraph:OnMouse()",5)
   // Write mouse pointer coordinates
      wxClientDC dc(this);
      PrepareDC(dc);
      mpParentFrame->PrepareDC(dc);

      wxPoint pos=event.GetPosition();
      const long x= dc.DeviceToLogicalX(pos.x);
      const long y= dc.DeviceToLogicalY(pos.y);
 
      wxCoord width,height;
      this->GetSize(&width, &height);

   if((x>width)||(y>height))
   {
      mMutex.Unlock();
      return;
   }
   //cout <<pos.x<<" "<<pos.y<<" "<<x<<" "<<y<<" "<<width<<" "<<height<<endl;
   const REAL x0=this->Screen2DataX(x);
   const REAL intensity=this->Screen2DataY(y);
   if(mpPattern->GetPowderPattern().GetRadiation().GetWavelengthType()==WAVELENGTH_TOF)
   {
      wxString str;
      const long pixel=
         (long)(mpPattern->GetPowderPattern().X2PixelCorr(x0));
      str.Printf(_T("tof=%6.2f    ,I=%12.2f.   pixel=#%ld"),x0,intensity,pixel);
      mMutex.Unlock();
      mpParentFrame->SetStatusText(str);
      mMutex.Lock();
      str.Printf(_T("d=%6.3fA"),0.5/mpPattern->GetPowderPattern().X2STOL(x0));
      mMutex.Unlock();
      mpParentFrame->SetStatusText(str,1);
      mMutex.Lock();
   }
   else
   {
      const REAL intensity=this->Screen2DataY(y);

      wxString str;
      const long pixel=
         (long)(mpPattern->GetPowderPattern().X2PixelCorr(x0*DEG2RAD));
      str.Printf(_T("2theta=%6.2f    ,I=%12.2f.   pixel=#%ld"),x0,intensity,pixel);
      // SetStatusText() triggers an OnPaint event ? Avoid deadlock by releasing the data mutex..
      mMutex.Unlock();
      mpParentFrame->SetStatusText(str);
      mMutex.Lock();
      str.Printf(_T("d=%6.3fA"),0.5/mpPattern->GetPowderPattern().X2STOL(x0*DEG2RAD));
      mMutex.Unlock();
      mpParentFrame->SetStatusText(str,1);
      mMutex.Lock();
   }
   if (event.Dragging() && event.LeftIsDown() && (!mIsDragging))
   {//Begin zooming
      mIsDragging=true;
      mDraggingX0=x0;
      mDraggingIntensity0=intensity;
      mMutex.Unlock();
      return;
   }
   if(event.LeftUp() && mIsDragging)
   {//Finished zooming !
      VFN_DEBUG_MESSAGE("WXPowderPatternGraph::OnMouse():Finished zooming...",5)
      mIsDragging=false;
      
      if( (fabs(x0-mDraggingX0)<.1) || (fabs(mDraggingIntensity0-intensity)< fabs(mMaxIntensity*.02)) )
      {
         mMutex.Unlock();
         return;
      }
      if(mDraggingIntensity0>intensity)
      {
         if(mDraggingIntensity0<0.)
         {
            mMutex.Unlock();
            return;
         }
         mMinIntensity=intensity;
         mMaxIntensity=mDraggingIntensity0;
      }
      else
      {
         if(intensity<0.)
         {
            mMutex.Unlock();
            return;
         }
         mMinIntensity=mDraggingIntensity0;
         mMaxIntensity=intensity;
      }
      if(mDraggingX0>x0)
      {
         mMinX=x0;
         mMaxX=mDraggingX0;
      }
      else
      {
         mMinX=mDraggingX0;
         mMaxX=x0;
      }
      mClockAxisLimits.Click();
      mMutex.Unlock();
      wxUpdateUIEvent event(ID_POWDER_GRAPH_NEW_PATTERN);
      wxPostEvent(this,event);
      return;
   }

   if(false==event.Dragging()) mIsDragging=false;

   if(event.LeftDClick())
   {//Reset axis range
      mMutex.Unlock();
      this->ResetAxisLimits();
      wxUpdateUIEvent event(ID_POWDER_GRAPH_NEW_PATTERN);
      wxPostEvent(this,event);
      event.Skip();
      return;
   }
   
   if(event.RightIsDown())
   {//popup menu
      mMutex.Unlock();
      if(mpPattern->GetPowderPattern().IsBeingRefined())
         mpPopUpMenu->Enable(ID_POWDERGRAPH_MENU_UPDATE, false);
      else
         mpPopUpMenu->Enable(ID_POWDERGRAPH_MENU_UPDATE, true);
      // Store x coordinate to allow adding/removing peaks manually
      mDraggingX0=x;
      this->PopupMenu(mpPopUpMenu, event.GetX(), event.GetY() );
      return;
   }
   mMutex.Unlock();
   event.Skip();
}
void WXPowderPatternGraph::OnMouseWheel(wxMouseEvent &event)
{
   VFN_DEBUG_ENTRY("WXPowderPatternGraph::OnMouseWheel()",6)
   wxMutexLocker mlock(mMutex);
   const long nbPoint=mX.numElements();
   if(event.GetWheelRotation()>=event.GetWheelDelta())
   {
      const REAL range=mMaxX-mMinX;
      mMaxX += range/8;
      if(mX(nbPoint-1)>mX(0))
      {
         if(mMaxX>=mX(nbPoint-1)) mMaxX=mX(nbPoint-1);
      }
      else
      {
         if(mMaxX>=mX(0)) mMaxX=mX(0);
      }
      mMinX=mMaxX-range;
   }
   if(event.GetWheelRotation()<=(-event.GetWheelDelta()))
   {
      const REAL range=mMaxX-mMinX;
      mMinX -= range/8;
      if(mX(nbPoint-1)>mX(0))
      {
         if(mMinX<mX(0)) mMinX=mX(0);
      }
      else 
      {
         if(mMinX<mX(nbPoint-1)) mMinX=mX(nbPoint-1);
      }
      mMaxX=mMinX+range;
   }
   mClockAxisLimits.Click();
   wxUpdateUIEvent ev(ID_POWDER_GRAPH_NEW_PATTERN);
   wxPostEvent(this,ev);
   VFN_DEBUG_EXIT("WXPowderPatternGraph::OnMouseWheel()",6)
}

void WXPowderPatternGraph::OnUpdate(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXPowderPatternGraph::OnUpdate()",6)
   mpPattern->CrystUpdate(true,true);
}

void WXPowderPatternGraph::OnToggleLabel(wxCommandEvent &event)
{
   VFN_DEBUG_MESSAGE("WXPowderPatternGraph::OnToggleLabel()",6)
   if(event.GetId()==ID_POWDERGRAPH_MENU_TOGGPEAK)
   {
      mDisplayPeak = !mDisplayPeak;
      if(mDisplayPeak) mpPopUpMenu->SetLabel(ID_POWDERGRAPH_MENU_TOGGPEAK, _T("Hide peaks"));
      else mpPopUpMenu->SetLabel(ID_POWDERGRAPH_MENU_TOGGPEAK, _T("Show peaks"));
   }
   if(event.GetId()==ID_POWDERGRAPH_MENU_TOGGLELABEL)
   {
      mDisplayLabel = !mDisplayLabel;
      if(mDisplayLabel) mpPopUpMenu->SetLabel(ID_POWDERGRAPH_MENU_TOGGLELABEL, _T("Hide labels"));
      else mpPopUpMenu->SetLabel(ID_POWDERGRAPH_MENU_TOGGLELABEL, _T("Show labels"));
   }
   this->Refresh(false);
}

void WXPowderPatternGraph::OnFindPeaks(wxCommandEvent& WXUNUSED(event))
{
   float dmin=1.5;
   while(true)
   {
      mPeakList=mpPattern->GetPowderPattern().FindPeaks(dmin,-1,1000);
      dmin*=0.75;
      if((mPeakList.GetPeakList().size()>30)||(dmin<0.3)) break;
   }
   
   const unsigned int nb=mPeakList.GetPeakList().size();
   //if(nb<5) return;
   
   mpPopUpMenu->Enable(ID_POWDERGRAPH_MENU_SAVEPEAKS, TRUE);
   mpPopUpMenu->Enable(ID_POWDERGRAPH_MENU_TOGGPEAK, TRUE);
   mpPopUpMenu->Enable(ID_POWDERGRAPH_MENU_INDEX, TRUE);
   mpPopUpMenu->Enable(ID_POWDERGRAPH_MENU_ADDPEAK, TRUE);
   mpPopUpMenu->Enable(ID_POWDERGRAPH_MENU_REMOVEPEAK, TRUE);
   
   // Keep lowest peaks
   if(mPeakList.GetPeakList().size()>40) mPeakList.GetPeakList().resize(40);
   //mpPattern->CrystUpdate(true,true);
   if(true)
   {// display 2nd derivative as 'calc' (will remain until an update to the pattern is made)
      CrystVector_REAL obsd2;
      obsd2=SavitzkyGolay(mObs,4,2);
      const float norm=-obsd2.min();
      obsd2 /= -norm;
      
      mCalc=obsd2;
      mCalc*=mObs.max();
   }

   this->Refresh(false);
}

void WXPowderPatternGraph::OnLoadPeaks(wxCommandEvent& WXUNUSED(event))
{
   wxFileDialog fn(this,_T("Choose a file"),_T(""),_T(""),_T("*"),wxOPEN);
   if(fn.ShowModal() != wxID_OK) return;
   ifstream f(fn.GetPath().ToAscii());
   if(!f) return;//:TODO:
   mPeakList.GetPeakList().clear();
   mPeakList.ImportDhklDSigmaIntensity(f);
   f.close();
   
   mpPopUpMenu->Enable(ID_POWDERGRAPH_MENU_SAVEPEAKS, TRUE);
   mpPopUpMenu->Enable(ID_POWDERGRAPH_MENU_TOGGPEAK, TRUE);
   mpPopUpMenu->Enable(ID_POWDERGRAPH_MENU_INDEX, TRUE);
   mpPopUpMenu->Enable(ID_POWDERGRAPH_MENU_ADDPEAK, TRUE);
   mpPopUpMenu->Enable(ID_POWDERGRAPH_MENU_REMOVEPEAK, TRUE);
   this->Refresh(false);
}

void WXPowderPatternGraph::OnSavePeaks(wxCommandEvent& WXUNUSED(event))
{
   wxFileDialog save(this,_T("Choose a file"),_T(""),_T(""),_T("*.txt"),wxSAVE | wxOVERWRITE_PROMPT);
   if(save.ShowModal() != wxID_OK) return;
   
   ofstream out(save.GetPath().ToAscii());
   if(!out) return;//:TODO:
   mPeakList.ExportDhklDSigmaIntensity(out);
   out.close();
}


void WXPowderPatternGraph::OnChangePeak(wxCommandEvent& event)
{
   if(event.GetId()==ID_POWDERGRAPH_MENU_REMOVEPEAK)
   {
      unsigned int idx=0;
      const float d=2*mpPattern->GetPowderPattern().X2STOL(this->Screen2DataX((long)mDraggingX0)*DEG2RAD);
      float dist=100.0;
      for(unsigned int i=0;i<mPeakList.GetPeakList().size();++i)
      {
         float x=mpPattern->GetPowderPattern().STOL2X(mPeakList.GetPeakList()[i].dobs/2);
         if(mpPattern->GetPowderPattern().GetRadiation().GetWavelengthType()!=WAVELENGTH_TOF)
               x *= RAD2DEG;
         x=this->Data2ScreenX(x);
         if(abs(x-mDraggingX0)<dist)
         {
            dist=abs(x-mDraggingX0);
            idx=i;
         }
         cout<<__FILE__<<":"<<__LINE__<<": 1/d0="<<d<<" peak #"<<i<<",d="<<1/mPeakList.GetPeakList()[i].dobs
             <<"("<<mDraggingX0<<","<<x<<"),  mindist="<<dist<<endl;
      }
      if(dist>5)       mpParentFrame->SetStatusText(_T("Could not find peak close enough"),1);
      else
      {
         if(mPeakList.GetPeakList().size()<5) mpPopUpMenu->Enable(ID_POWDERGRAPH_MENU_INDEX, FALSE);
         if(mPeakList.GetPeakList().size()==0)
         {
            mpPopUpMenu->Enable(ID_POWDERGRAPH_MENU_TOGGPEAK, FALSE);
            mpPopUpMenu->Enable(ID_POWDERGRAPH_MENU_SAVEPEAKS, FALSE);
            mpPopUpMenu->Enable(ID_POWDERGRAPH_MENU_REMOVEPEAK, FALSE);
         }
         wxString buf;
         buf.Printf(_T("Removing peak at d=%6.3f"),1/mPeakList.GetPeakList()[idx].dobs);
         mpParentFrame->SetStatusText(buf,1);
         mPeakList.RemovePeak(idx);
         this->Refresh(false);
      }
   }
   if(event.GetId()==ID_POWDERGRAPH_MENU_ADDPEAK)
   {
      float d,sig;
      if(mpPattern->GetPowderPattern().GetRadiation().GetWavelengthType()!=WAVELENGTH_TOF)
      {
         d  =    2*mpPattern->GetPowderPattern().X2STOL(this->Screen2DataX((long)mDraggingX0  )*DEG2RAD);
         sig=abs(2*mpPattern->GetPowderPattern().X2STOL(this->Screen2DataX((long)mDraggingX0+1)*DEG2RAD)-d);
      }
      else 
      {
         d  =    2*mpPattern->GetPowderPattern().X2STOL(this->Screen2DataX((long)mDraggingX0  ));
         sig=abs(2*mpPattern->GetPowderPattern().X2STOL(this->Screen2DataX((long)mDraggingX0+1))-d);
      }
      
      mpPopUpMenu->Enable(ID_POWDERGRAPH_MENU_SAVEPEAKS, TRUE);
      mpPopUpMenu->Enable(ID_POWDERGRAPH_MENU_TOGGPEAK, TRUE);
      if(mPeakList.GetPeakList().size()>=5) mpPopUpMenu->Enable(ID_POWDERGRAPH_MENU_INDEX, TRUE);
      mpPopUpMenu->Enable(ID_POWDERGRAPH_MENU_REMOVEPEAK, TRUE);
      
      wxString buf;
      buf.Printf(_T("Added peak at d=%6.3f"),1/d);
      mPeakList.AddPeak(d,1.0,sig);
      mpParentFrame->SetStatusText(buf,1);
      this->Refresh(false);
      mPeakList.Print(cout);
   }
}
//////////////////////////////////////// WXCellExplorer /////////////////////

class WXCellExplorer:public wxWindow
{
   public:
      WXCellExplorer(wxWindow *parent,PeakList &peaklist,WXPowderPatternGraph *graph=NULL);
      ~WXCellExplorer();
      /// Start indexing
      void OnIndex(wxCommandEvent &event);
      /// Select cell
      void OnSelectCell(wxCommandEvent &event);
      /// Apply cell
      void OnApplyCell(wxCommandEvent &event);
      /// Choose crystal to apply cell to
      void OnChooseCrystal(wxCommandEvent &event);
      /// User chose to automatically run Le Bail & profile fitting,
      ///make sure there is a background and a crystal structure selected
      void OnAutoLeBail(wxCommandEvent &event);
   private:
      WXPowderPatternGraph *mpGraph;
      PeakList *mpPeakList;
      CellExplorer *mpCellExplorer;
      wxRadioBox *mpAlgorithm;
      wxRadioBox *mpBravais;
      wxListBox *mpCell;
      wxTextCtrl *mpLog;
      wxTextCtrl *mpLengthMin,*mpLengthMax;
      wxTextCtrl *mpAngleMin,*mpAngleMax;
      wxTextCtrl *mpVolumeMin,*mpVolumeMax;
      wxTextCtrl *mpNbSpurious;
      wxTextCtrl *mpNbPeak;
      wxTextCtrl *mpErrorD;
      WXFieldChoice *mpFieldCrystal;
      wxTextCtrl *mpStopOnScore,*mpStopOnDepth;
      wxTextCtrl *mpReportOnScore,*mpReportOnDepth;
      Crystal *mpCrystal;
      PowderPatternDiffraction *mpDiff;
      wxCheckBox *mpWeakDiffraction;
      wxCheckBox *mpContinueOnSolution;
      wxCheckBox *mpTryCenteredLattice;
      wxCheckBox *mpAutomaticLeBail;
      DECLARE_EVENT_TABLE()
};
static const long ID_CELLEXPLORER_INDEX= WXCRYST_ID();
static const long ID_CELLEXPLORER_INDEX_QUICK= WXCRYST_ID();
static const long ID_CELLEXPLORER_WEAK= WXCRYST_ID();
static const long ID_CELLEXPLORER_SELECTCELL= WXCRYST_ID();
static const long ID_CELLEXPLORER_APPLYCELL= WXCRYST_ID();
static const long ID_CELLEXPLORER_CHOOSECRYSTAL= WXCRYST_ID();
static const long ID_CELLEXPLORER_LEBAIL= WXCRYST_ID();
static const long ID_CELLEXPLORER_CENTERED= WXCRYST_ID();

BEGIN_EVENT_TABLE(WXCellExplorer, wxWindow)
   EVT_BUTTON(ID_CELLEXPLORER_INDEX,             WXCellExplorer::OnIndex)
   EVT_BUTTON(ID_CELLEXPLORER_INDEX_QUICK,       WXCellExplorer::OnIndex)
   EVT_LISTBOX(ID_CELLEXPLORER_SELECTCELL,       WXCellExplorer::OnSelectCell)
   EVT_LISTBOX_DCLICK(ID_CELLEXPLORER_SELECTCELL,WXCellExplorer::OnApplyCell)
   EVT_BUTTON(ID_CELLEXPLORER_APPLYCELL,         WXCellExplorer::OnApplyCell)
   EVT_BUTTON(ID_CELLEXPLORER_CHOOSECRYSTAL,     WXCellExplorer::OnChooseCrystal)
   EVT_CHECKBOX(ID_CELLEXPLORER_LEBAIL,          WXCellExplorer::OnAutoLeBail)
END_EVENT_TABLE()

WXCellExplorer::WXCellExplorer(wxWindow *parent, PeakList &peaklist, WXPowderPatternGraph *graph):
wxWindow(parent,-1),mpGraph(graph),mpPeakList(&peaklist),mpCellExplorer(0),mpCrystal(0),mpDiff(0)
{
   wxBoxSizer *pSizer1=new wxBoxSizer(wxHORIZONTAL);
   this->SetSizer(pSizer1);
   
   wxNotebook *pNotebook = new wxNotebook(this, -1);
   
   pSizer1->Add(pNotebook,0,wxALIGN_TOP);
   
   // Quick interface
      wxWindow *pQuick=new wxWindow(pNotebook,-1);
      pNotebook->AddPage(pQuick,_T("Quick"));
      
      wxStaticBoxSizer *pSizerQuick=new wxStaticBoxSizer(wxVERTICAL,pQuick);
      
      wxButton *pQuickButtonIndex=new wxButton(pQuick,ID_CELLEXPLORER_INDEX_QUICK,_T("Find cell!"));
      pSizerQuick->Add(pQuickButtonIndex,0,wxALIGN_CENTER);
      
      mpWeakDiffraction=new wxCheckBox(pQuick,ID_CELLEXPLORER_WEAK,_T("Weak Diffraction (scan larger volume)"));
      pSizerQuick->Add(mpWeakDiffraction,0,wxALIGN_CENTER);
      
      mpContinueOnSolution=new wxCheckBox(pQuick,ID_CELLEXPLORER_WEAK,_T("Continue exploring after solution"));
      pSizerQuick->Add(mpContinueOnSolution,0,wxALIGN_CENTER);
      
      mpTryCenteredLattice=new wxCheckBox(pQuick,ID_CELLEXPLORER_CENTERED,_T("Try Centered Lattices"));
      pSizerQuick->Add(mpTryCenteredLattice,0,wxALIGN_CENTER);
      
      pQuick->SetSizer(pSizerQuick);
      pSizerQuick->Fit(pQuick);
      pSizerQuick->RecalcSizes();
      pQuick->Layout();
   // Advanced interface
      wxWindow *pAdvanced=new wxWindow(pNotebook,-1);
      
      wxStaticBoxSizer *pSizerAdvanced=new wxStaticBoxSizer(wxVERTICAL,pAdvanced);
      
      wxButton *pButton1=new wxButton(pAdvanced,ID_CELLEXPLORER_INDEX,_T("Find cell!"));
      pSizerAdvanced->Add(pButton1,0,wxALIGN_CENTER);
      
      wxBoxSizer *pLengthSizer=new wxBoxSizer(wxHORIZONTAL);
      wxStaticText *pLengthText=new wxStaticText(pAdvanced,-1,_T("Length min, max (A):"));
      pLengthSizer->Add(pLengthText,0,wxALIGN_CENTER);
      mpLengthMin=new wxTextCtrl(pAdvanced,-1,_T("3"),wxDefaultPosition,wxSize(30,-1),0,
                                 wxTextValidator(wxFILTER_NUMERIC));
      pLengthSizer->Add(mpLengthMin,0,wxALIGN_CENTER);
      mpLengthMax=new wxTextCtrl(pAdvanced,-1,_T("25"),wxDefaultPosition,wxSize(30,-1),0,
                                 wxTextValidator(wxFILTER_NUMERIC));
      pLengthSizer->Add(mpLengthMax,0,wxALIGN_CENTER);
      pSizerAdvanced->Add(pLengthSizer,0,wxALIGN_CENTER);
   
      wxBoxSizer *pAngleSizer=new wxBoxSizer(wxHORIZONTAL);
      wxStaticText *pAngleText=new wxStaticText(pAdvanced,-1,_T("Angle max(90< <179°):"));
      pAngleSizer->Add(pAngleText,0,wxALIGN_CENTER);
      //mpAngleMin=new wxTextCtrl(this,-1,"90",wxDefaultPosition,wxSize(40,-1),0,
      //                          wxTextValidator(wxFILTER_NUMERIC));
      //pAngleSizer->Add(mpAngleMin,0,wxALIGN_CENTER);
      mpAngleMax=new wxTextCtrl(pAdvanced,-1,_T("130"),wxDefaultPosition,wxSize(40,-1),0,
                              wxTextValidator(wxFILTER_NUMERIC));
      pAngleSizer->Add(mpAngleMax,0,wxALIGN_CENTER);
      pSizerAdvanced->Add(pAngleSizer,0,wxALIGN_CENTER);
   
      wxBoxSizer *pVolumeSizer=new wxBoxSizer(wxHORIZONTAL);
      wxStaticText *pVolumeText=new wxStaticText(pAdvanced,-1,_T("Volume min, max (A3):"));
      pVolumeSizer->Add(pVolumeText,0,wxALIGN_CENTER);
      mpVolumeMin=new wxTextCtrl(pAdvanced,-1,_T("10"),wxDefaultPosition,wxSize(50,-1),0,
                                 wxTextValidator(wxFILTER_NUMERIC));
      pVolumeSizer->Add(mpVolumeMin,0,wxALIGN_CENTER);
      mpVolumeMax=new wxTextCtrl(pAdvanced,-1,_T("2500"),wxDefaultPosition,wxSize(50,-1),0,
                                 wxTextValidator(wxFILTER_NUMERIC));
      pVolumeSizer->Add(mpVolumeMax,0,wxALIGN_CENTER);
      pSizerAdvanced->Add(pVolumeSizer,0,wxALIGN_CENTER);
   
      wxBoxSizer *pSpuriousSizer=new wxBoxSizer(wxHORIZONTAL);
      wxStaticText *pSpuriousText=new wxStaticText(pAdvanced,-1,_T("Nb spurious lines:"));
      pSpuriousSizer->Add(pSpuriousText,0,wxALIGN_CENTER);
      mpNbSpurious=new wxTextCtrl(pAdvanced,-1,_T("0"),wxDefaultPosition,wxSize(40,-1),0,
                                 wxTextValidator(wxFILTER_NUMERIC));
      pSpuriousSizer->Add(mpNbSpurious,0,wxALIGN_CENTER);
      pSizerAdvanced->Add(pSpuriousSizer,0,wxALIGN_CENTER);
      
      wxBoxSizer *pNbPeakSizer=new wxBoxSizer(wxHORIZONTAL);
      wxStaticText *pNbPeakText=new wxStaticText(pAdvanced,-1,_T("Use Nb Peaks:"));
      pNbPeakSizer->Add(pNbPeakText,0,wxALIGN_CENTER);
      mpNbPeak=new wxTextCtrl(pAdvanced,-1,_T("20"),wxDefaultPosition,wxSize(40,-1),0,
                                 wxTextValidator(wxFILTER_NUMERIC));
      pNbPeakSizer->Add(mpNbPeak,0,wxALIGN_CENTER);
      pSizerAdvanced->Add(pNbPeakSizer,0,wxALIGN_CENTER);
      
      wxBoxSizer *pStopSizer=new wxBoxSizer(wxHORIZONTAL);
      wxStaticText* pStopOnScoreText=new wxStaticText(pAdvanced,-1,_T("Stop on Score>"));
      pStopSizer->Add(pStopOnScoreText,0,wxALIGN_CENTER);
      mpStopOnScore=new wxTextCtrl(pAdvanced,-1,_T("50"),wxDefaultPosition,wxSize(50,-1),0,
                                 wxTextValidator(wxFILTER_NUMERIC));
      pStopSizer->Add(mpStopOnScore,0,wxALIGN_CENTER);
      
      wxStaticText* pStopOnDepthText=new wxStaticText(pAdvanced,-1,_T("and depth>="));
      pStopSizer->Add(pStopOnDepthText,0,wxALIGN_CENTER);
      mpStopOnDepth=new wxTextCtrl(pAdvanced,-1,_T("6"),wxDefaultPosition,wxSize(30,-1),0,
                                 wxTextValidator(wxFILTER_NUMERIC));
      pStopSizer->Add(mpStopOnDepth,0,wxALIGN_CENTER);
      pSizerAdvanced->Add(pStopSizer,0,wxALIGN_CENTER);
      
      
      wxBoxSizer *pReportSizer=new wxBoxSizer(wxHORIZONTAL);
      wxStaticText* pReportOnScoreText=new wxStaticText(pAdvanced,-1,_T("Report score>"));
      pReportSizer->Add(pReportOnScoreText,0,wxALIGN_CENTER);
      mpReportOnScore=new wxTextCtrl(pAdvanced,-1,_T("10"),wxDefaultPosition,wxSize(50,-1),0,
                                 wxTextValidator(wxFILTER_NUMERIC));
      pReportSizer->Add(mpReportOnScore,0,wxALIGN_CENTER);
      
      wxStaticText* pReportOnDepthText=new wxStaticText(pAdvanced,-1,_T("or depth>="));
      pReportSizer->Add(pReportOnDepthText,0,wxALIGN_CENTER);
      mpReportOnDepth=new wxTextCtrl(pAdvanced,-1,_T("4"),wxDefaultPosition,wxSize(50,-1),0,
                                 wxTextValidator(wxFILTER_NUMERIC));
      pReportSizer->Add(mpReportOnDepth,0,wxALIGN_CENTER);
      pSizerAdvanced->Add(pReportSizer,0,wxALIGN_CENTER);
      
      wxBoxSizer *pErrorSizer=new wxBoxSizer(wxHORIZONTAL);
      wxStaticText* pErrorText=new wxStaticText(pAdvanced,-1,_T("delta(d)/d^2 error:"));
      pErrorSizer->Add(pErrorText,0,wxALIGN_CENTER);
      mpErrorD=new wxTextCtrl(pAdvanced,-1,_T("0"),wxDefaultPosition,wxSize(50,-1),0,
                              wxTextValidator(wxFILTER_NUMERIC));
      pErrorSizer->Add(mpErrorD,0,wxALIGN_CENTER);
      pSizerAdvanced->Add(pErrorSizer,0,wxALIGN_CENTER);
   
      wxArrayString bravaisChoices;
      bravaisChoices.Add(_T("Triclinic"));
      bravaisChoices.Add(_T("Monoclinic"));
      bravaisChoices.Add(_T("Orthorombic"));
      bravaisChoices.Add(_T("Hexagonal"));
      bravaisChoices.Add(_T("Rhomboedral"));
      bravaisChoices.Add(_T("Tetragonal"));
      bravaisChoices.Add(_T("Cubic"));
      mpBravais=new wxRadioBox((wxWindow*)pAdvanced,-1,_T("Crystal System"),wxDefaultPosition,wxDefaultSize,bravaisChoices,0,wxRA_SPECIFY_ROWS);
      mpBravais->SetSelection(2);
      //mpBravais->Enable(0,false);
      pSizerAdvanced->Add(mpBravais,0,wxALIGN_CENTER);
      
      wxArrayString algoChoices;
      algoChoices.Add(_T("DICVOL"));
      //algoChoices.Add(_T("Differential Evolution"));  // :TODO: re-enable after testing
      #if 0
      mpAlgorithm=new wxRadioBox(pAdvanced,-1,_T("Algorithm"),wxDefaultPosition,wxDefaultSize,algoChoices,0,wxRA_SPECIFY_ROWS);
      mpAlgorithm->Enable(1,false);
      pSizerAdvanced->Add(mpAlgorithm,0,wxALIGN_CENTER);
      #endif
      pAdvanced->SetSizer(pSizerAdvanced);
      pSizerAdvanced->Fit(pAdvanced);
      pSizerAdvanced->RecalcSizes();
      pAdvanced->Layout();
      
      pNotebook->AddPage(pAdvanced,_T("Advanced"));
      
   pNotebook->Layout();
   // Solutions & log
      wxBoxSizer *pSizer2=new wxBoxSizer(wxVERTICAL);
      
      //wxButton *pButton2=new wxButton(this,ID_CELLEXPLORER_APPLYCELL,"Apply selected cell");
      //pSizer2->Add(pButton2,0,wxALIGN_CENTER);
   
      mpFieldCrystal=new WXFieldChoice(this,ID_CELLEXPLORER_CHOOSECRYSTAL,"Choose crystal to apply selected cell to:",200);
      pSizer2->Add(mpFieldCrystal,0,wxALIGN_CENTER);
   
      mpAutomaticLeBail=new wxCheckBox(this,ID_CELLEXPLORER_LEBAIL,_T("Automatic Profile Fitting (Le Bail)"));
      pSizer2->Add(mpAutomaticLeBail,0,wxALIGN_CENTER);
      
      wxArrayString cells;
      mpCell=new wxListBox(this,ID_CELLEXPLORER_SELECTCELL,wxDefaultPosition,wxSize(750,400),cells,wxLB_SINGLE);
      mpCell->SetFont(wxFont(9,wxTELETYPE,wxFONTSTYLE_NORMAL,wxFONTWEIGHT_NORMAL));
      pSizer2->Add(mpCell,0,wxALIGN_CENTER);
      
      mpLog =new wxTextCtrl(this,-1,_T(""),wxDefaultPosition,wxSize(750,250),wxTE_MULTILINE|wxTE_READONLY|wxTE_DONTWRAP);
      mpLog->SetFont(wxFont(9,wxTELETYPE,wxFONTSTYLE_NORMAL,wxFONTWEIGHT_NORMAL));
      pSizer2->Add(mpLog,0,wxALIGN_CENTER);
      
      pSizer1->Add(pSizer2,0,wxALIGN_TOP);
   
   this->Layout();
   pSizer1->Fit(this->GetParent());
   pSizer1->SetSizeHints(this);

   // Estimate volume from number of peaks at a given dmin
   // See J. Appl. Cryst. 20 (1987), 161
   unsigned int nb=mpPeakList->GetPeakList().size();
   if(nb>20) nb=20;// Just use 20 - beyond that we probably have a lot of weak peaks missed
   const float dmin=mpPeakList->GetPeakList()[nb-1].dobs;
   const float dmax=mpPeakList->GetPeakList()[0].dobs/10;// /10: assume no peaks at lower resolution
   mpLog->AppendText(wxString::Format(_T("Predicted unit vell volume from %2u observed peaks between: dmax=%6.3f A-> dmin=%6.3fA\n"),nb,1/dmax,1/dmin));
   mpLog->AppendText(wxString::Format(_T("(Assuming observed lines represent 120%% down to 30%% of existing reflections)\n")));
   mpLog->AppendText(wxString::Format(_T("  Cubic P         v=%6.0f -> %6.0f A\n"),EstimateCellVolume(dmin,dmax,nb,CUBIC      ,LATTICE_P,1.2),EstimateCellVolume(dmin,dmax,nb,CUBIC      ,LATTICE_P,0.3)));
   mpLog->AppendText(wxString::Format(_T("  Cubic I         v=%6.0f -> %6.0f A\n"),EstimateCellVolume(dmin,dmax,nb,CUBIC      ,LATTICE_I,1.2),EstimateCellVolume(dmin,dmax,nb,CUBIC      ,LATTICE_I,0.3)));
   mpLog->AppendText(wxString::Format(_T("  Cubic F         v=%6.0f -> %6.0f A\n"),EstimateCellVolume(dmin,dmax,nb,CUBIC      ,LATTICE_F,1.2),EstimateCellVolume(dmin,dmax,nb,CUBIC      ,LATTICE_F,0.3)));
   mpLog->AppendText(wxString::Format(_T("  Tetragonal P    v=%6.0f -> %6.0f A\n"),EstimateCellVolume(dmin,dmax,nb,TETRAGONAL ,LATTICE_P,1.2),EstimateCellVolume(dmin,dmax,nb,TETRAGONAL ,LATTICE_P,0.3)));
   mpLog->AppendText(wxString::Format(_T("  Tetragonal I    v=%6.0f -> %6.0f A\n"),EstimateCellVolume(dmin,dmax,nb,TETRAGONAL ,LATTICE_I,1.2),EstimateCellVolume(dmin,dmax,nb,TETRAGONAL ,LATTICE_I,0.3)));
   mpLog->AppendText(wxString::Format(_T("  Orthorombic P   v=%6.0f -> %6.0f A\n"),EstimateCellVolume(dmin,dmax,nb,ORTHOROMBIC,LATTICE_P,1.2),EstimateCellVolume(dmin,dmax,nb,ORTHOROMBIC,LATTICE_P,0.3)));
   mpLog->AppendText(wxString::Format(_T("  Orthorombic I,C v=%6.0f -> %6.0f A\n"),EstimateCellVolume(dmin,dmax,nb,ORTHOROMBIC,LATTICE_I,1.2),EstimateCellVolume(dmin,dmax,nb,ORTHOROMBIC,LATTICE_I,0.3)));
   mpLog->AppendText(wxString::Format(_T("  Orthorombic F   v=%6.0f -> %6.0f A\n"),EstimateCellVolume(dmin,dmax,nb,ORTHOROMBIC,LATTICE_F,1.2),EstimateCellVolume(dmin,dmax,nb,ORTHOROMBIC,LATTICE_I,0.3)));
   mpLog->AppendText(wxString::Format(_T("  Hexagonal       v=%6.0f -> %6.0f A\n"),EstimateCellVolume(dmin,dmax,nb,HEXAGONAL  ,LATTICE_P,1.2),EstimateCellVolume(dmin,dmax,nb,HEXAGONAL  ,LATTICE_P,0.3)));
   mpLog->AppendText(wxString::Format(_T("  Monoclinic P    v=%6.0f -> %6.0f A\n"),EstimateCellVolume(dmin,dmax,nb,MONOCLINIC ,LATTICE_P,1.2),EstimateCellVolume(dmin,dmax,nb,MONOCLINIC ,LATTICE_P,0.3)));
   mpLog->AppendText(wxString::Format(_T("  Monoclinic C    v=%6.0f -> %6.0f A\n"),EstimateCellVolume(dmin,dmax,nb,MONOCLINIC ,LATTICE_C,1.2),EstimateCellVolume(dmin,dmax,nb,MONOCLINIC ,LATTICE_C,0.3)));
   mpLog->AppendText(wxString::Format(_T("  Triclinic       v=%6.0f -> %6.0f A\n"),EstimateCellVolume(dmin,dmax,nb,TRICLINIC  ,LATTICE_P,1.2),EstimateCellVolume(dmin,dmax,nb,TRICLINIC  ,LATTICE_P,0.3)));
}

WXCellExplorer::~WXCellExplorer()
{
   if(mpDiff!=0)
   {
      mpDiff->SetExtractionMode(false,false);
      mpDiff->UpdateDisplay();
   }
}

void WXCellExplorer::OnIndex(wxCommandEvent &event)
{
   
   Chronometer chrono;
   if(event.GetId()==ID_CELLEXPLORER_INDEX_QUICK)
   {
      // Erase spurious record
      for(vector<PeakList::hkl>::iterator pos=mpPeakList->mvHKL.begin();pos!=mpPeakList->mvHKL.end();++pos)
         pos->isSpurious=false;
      
      // Use at most 30 reflections for indexing
      PeakList peaklist=*mpPeakList;
      if(peaklist.GetPeakList().size()>20) peaklist.GetPeakList().resize(20);
      
      // Estimate volume from number of peaks at a given dmin
      // See J. Appl. Cryst. 20 (1987), 161
      unsigned int nb=mpPeakList->GetPeakList().size();
      if(nb>20) nb=20;// Just use 20 - beyond that we probably have a lot of weak peaks
      float dmin=mpPeakList->GetPeakList()[nb-1].dobs;
      const float dmax=mpPeakList->GetPeakList()[0].dobs/10;//assume there are no peaks at lower resolution
      mpLog->AppendText(wxString::Format(_T("Predicting volumes from %2u peaks between d=%6.3f and d=%6.3f\n"),nb,1/dmax,1/dmin));
      mpLog->AppendText(wxString::Format(_T("Starting indexing using %2u peaks\n"),nb));

      mpCellExplorer = new CellExplorer(peaklist,CUBIC,0);
      mpCellExplorer->SetLengthMinMax(3,25);
      mpCellExplorer->SetAngleMinMax(90*DEG2RAD,140*DEG2RAD);
      mpCellExplorer->SetD2Error(0);
      
      float weak_f=1.0;
      if(mpWeakDiffraction->GetValue()) weak_f=0.5;
      const bool continueOnSolution=mpContinueOnSolution->GetValue();
      
      const bool noCentered=!(mpTryCenteredLattice->GetValue());
      
      const float        stopOnScore=50, reportOnScore=10;
      const unsigned int stopOnDepth=6+int(continueOnSolution),   reportOnDepth=4;
      
      unsigned int nbSpurious=0;
      wxProgressDialog dlgProgress(_T("Indexing..."),_T("Starting Indexing in Quick Mode"),
                                   7,this,wxPD_AUTO_HIDE|wxPD_ELAPSED_TIME|wxPD_CAN_ABORT);
      while(nbSpurious<=3)
      {
         float t0,minv,maxv,lengthmax;
         mpCellExplorer->SetNbSpurious(nbSpurious);
         CrystalCentering cent;
         char centc;
         for(int lat=0;lat<=2;++lat)
         {
            switch(lat)
            {//LATTICE_P,LATTICE_I,LATTICE_A,LATTICE_B,LATTICE_C,LATTICE_F
               case 0:cent=LATTICE_P;centc='P';break;
               case 1:cent=LATTICE_I;centc='I';break;
               case 2:cent=LATTICE_F;centc='F';break;
            }
            minv=EstimateCellVolume(dmin,dmax,nb,CUBIC      ,cent,1.5);
            maxv=EstimateCellVolume(dmin,dmax,nb,CUBIC      ,cent,0.4*weak_f);
            mpCellExplorer->SetVolumeMinMax(minv,maxv);
            lengthmax=pow(maxv,(float)(1/3.0))*3;
            if(lengthmax<25)lengthmax=25;
            mpCellExplorer->SetLengthMinMax(3,lengthmax);
            mpCellExplorer->SetCrystalSystem(CUBIC);
            mpCellExplorer->SetCrystalCentering(cent);
            mpLog->AppendText(wxString::Format(_T("CUBIC %c      : V= %6.0f -> %6.0f A^3, max length=%6.2fA"),centc,minv,maxv,lengthmax));
            t0=chrono.seconds();
            if(dlgProgress.Update(0,wxString::Format(_T("CUBIC %c (%d spurious), V=%6.0f-%6.0f, l<%6.2fA\n")
                                                   _T("Best Score=%6.1f"),centc,
                                                   nbSpurious,minv,maxv,lengthmax,mpCellExplorer->GetBestScore()))==false) break;
            mpCellExplorer->DicVol(reportOnScore,reportOnDepth,stopOnScore,stopOnDepth);
            mpLog->AppendText(wxString::Format(_T(" -> %3u sols in %6.2fs, best score=%6.1f\n"),
                     mpCellExplorer->GetSolutions().size(),chrono.seconds()-t0,mpCellExplorer->GetBestScore()));
            mpLog->Update();
            if(noCentered) break;
         }
         for(int lat=0;lat<=1;++lat)
         if((mpCellExplorer->GetBestScore()<=stopOnScore)||continueOnSolution)
         {
            switch(lat)
            {//LATTICE_P,LATTICE_I,LATTICE_A,LATTICE_B,LATTICE_C,LATTICE_F
               case 0:cent=LATTICE_P;centc='P';break;
               case 1:cent=LATTICE_I;centc='I';break;
            }
            minv=EstimateCellVolume(dmin,dmax,nb,TETRAGONAL,cent,1.5);
            maxv=EstimateCellVolume(dmin,dmax,nb,TETRAGONAL,cent,0.4*weak_f);
            mpCellExplorer->SetVolumeMinMax(minv,maxv);
            float lengthmax=pow(maxv,(float)(1/3.0))*3;
            if(lengthmax<25)lengthmax=25;
            mpCellExplorer->SetLengthMinMax(3,lengthmax);
            mpCellExplorer->SetCrystalSystem(TETRAGONAL);
            mpCellExplorer->SetCrystalCentering(cent);
            mpLog->AppendText(wxString::Format(_T("TETRAGONAL %c : V= %6.0f -> %6.0f A^3, max length=%6.2fA"),centc,minv,maxv,lengthmax));
            t0=chrono.seconds();
            if(dlgProgress.Update(1,wxString::Format(_T("TETRAGONAL %c (%d spurious), V=%6.0f-%6.0f, l<%6.2fA\n")
                                                   _T("Best Score=%6.1f"),centc,
                                                   nbSpurious,minv,maxv,lengthmax,mpCellExplorer->GetBestScore()))==false) break;
            mpCellExplorer->DicVol(reportOnScore,reportOnDepth,stopOnScore,stopOnDepth);
            mpLog->AppendText(wxString::Format(_T(" -> %3u sols in %6.2fs, best score=%6.1f\n"),
                     mpCellExplorer->GetSolutions().size(),chrono.seconds()-t0,mpCellExplorer->GetBestScore()));
            mpLog->Update();
            if(noCentered) break;
         }
         if((mpCellExplorer->GetBestScore()<=stopOnScore)||continueOnSolution)
         {
            minv=EstimateCellVolume(dmin,dmax,nb,RHOMBOEDRAL,LATTICE_P,1.5);
            maxv=EstimateCellVolume(dmin,dmax,nb,RHOMBOEDRAL,LATTICE_P,0.4*weak_f);
            mpCellExplorer->SetVolumeMinMax(minv,maxv);
            lengthmax=pow(maxv,(float)(1/3.0))*3;
            if(lengthmax<25)lengthmax=25;
            mpCellExplorer->SetLengthMinMax(3,lengthmax);
            mpCellExplorer->SetCrystalSystem(RHOMBOEDRAL);
            mpCellExplorer->SetCrystalCentering(LATTICE_P);
            mpLog->AppendText(wxString::Format(_T("RHOMBOEDRAL  : V= %6.0f -> %6.0f A^3, max length=%6.2fA"),minv,maxv,lengthmax));
            t0=chrono.seconds();
            if(dlgProgress.Update(2,wxString::Format(_T("RHOMBOEDRAL (%d spurious), V=%6.0f-%6.0f, l<%6.2fA\n")
                                                   _T("Best Score=%6.1f"),
                                                   nbSpurious,minv,maxv,lengthmax,mpCellExplorer->GetBestScore()))==false) break;
            mpCellExplorer->DicVol(reportOnScore,reportOnDepth,stopOnScore,stopOnDepth);
            mpLog->AppendText(wxString::Format(_T(" -> %3u sols in %6.2fs, best score=%6.1f\n"),
                     mpCellExplorer->GetSolutions().size(),chrono.seconds()-t0,mpCellExplorer->GetBestScore()));
            mpLog->Update();
         }
         if((mpCellExplorer->GetBestScore()<=stopOnScore)||continueOnSolution)
         {
            minv=EstimateCellVolume(dmin,dmax,nb,HEXAGONAL,LATTICE_P,1.5);
            maxv=EstimateCellVolume(dmin,dmax,nb,HEXAGONAL,LATTICE_P,0.4*weak_f);
            mpCellExplorer->SetVolumeMinMax(minv,maxv);
            lengthmax=pow(maxv,(float)(1/3.0))*3;
            if(lengthmax<25)lengthmax=25;
            mpCellExplorer->SetLengthMinMax(3,lengthmax);
            mpCellExplorer->SetCrystalSystem(HEXAGONAL);
            mpCellExplorer->SetCrystalCentering(LATTICE_P);
            mpLog->AppendText(wxString::Format(_T("HEXAGONAL    : V= %6.0f -> %6.0f A^3, max length=%6.2fA"),minv,maxv,lengthmax));
            t0=chrono.seconds();
            if(dlgProgress.Update(3,wxString::Format(_T("HEXAGONAL (%d spurious), V=%6.0f-%6.0f, l<%6.2fA\n")
                                                   _T("Best Score=%6.1f"),
                                                   nbSpurious,minv,maxv,lengthmax,mpCellExplorer->GetBestScore()))==false) break;
            mpCellExplorer->DicVol(reportOnScore,reportOnDepth,stopOnScore,stopOnDepth);
            mpLog->AppendText(wxString::Format(_T(" -> %3u sols in %6.2fs, best score=%6.1f\n"),
                     mpCellExplorer->GetSolutions().size(),chrono.seconds()-t0,mpCellExplorer->GetBestScore()));
            mpLog->Update();
         }
         for(int lat=0;lat<=5;++lat)
         if((mpCellExplorer->GetBestScore()<=stopOnScore)||continueOnSolution)
         {
            switch(lat)
            {//LATTICE_P,LATTICE_I,LATTICE_A,LATTICE_B,LATTICE_C,LATTICE_F
               case 0:cent=LATTICE_P;centc='P';break;
               case 1:cent=LATTICE_I;centc='I';break;
               case 2:cent=LATTICE_A;centc='A';break;
               case 3:cent=LATTICE_B;centc='B';break;
               case 4:cent=LATTICE_C;centc='C';break;
               case 5:cent=LATTICE_F;centc='F';break;
            }
            minv=EstimateCellVolume(dmin,dmax,nb,ORTHOROMBIC,cent,1.5);
            maxv=EstimateCellVolume(dmin,dmax,nb,ORTHOROMBIC,cent,0.4*weak_f);
            mpCellExplorer->SetVolumeMinMax(minv,maxv);
            lengthmax=pow(maxv,(float)(1/3.0))*3;
            if(lengthmax<25)lengthmax=25;
            mpCellExplorer->SetLengthMinMax(3,lengthmax);
            mpCellExplorer->SetCrystalSystem(ORTHOROMBIC);
            mpCellExplorer->SetCrystalCentering(cent);
            mpLog->AppendText(wxString::Format(_T("ORTHOROMBIC %c: V= %6.0f -> %6.0f A^3, max length=%6.2fA"),centc,minv,maxv,lengthmax));
            t0=chrono.seconds();
            if(dlgProgress.Update(4,wxString::Format(_T("ORTHOROMBIC %c (%d spurious), V=%6.0f-%6.0f, l<%6.2fA\n")
                                                   _T("Best Score=%6.1f"),centc,
                                                   nbSpurious,minv,maxv,lengthmax,mpCellExplorer->GetBestScore()))==false) break;
            mpCellExplorer->DicVol(reportOnScore,reportOnDepth,stopOnScore,stopOnDepth);
            mpLog->AppendText(wxString::Format(_T(" -> %3u sols in %6.2fs, best score=%6.1f\n"),
                     mpCellExplorer->GetSolutions().size(),chrono.seconds()-t0,mpCellExplorer->GetBestScore()));
            mpLog->Update();
            if(noCentered) break;
         }
         for(int lat=0;lat<=3;++lat)
         if((mpCellExplorer->GetBestScore()<=stopOnScore)||continueOnSolution)
         {
            switch(lat)
            {//LATTICE_P,LATTICE_I,LATTICE_A,LATTICE_B,LATTICE_C,LATTICE_F
               case 0:cent=LATTICE_P;centc='P';break;
               case 1:cent=LATTICE_C;centc='C';break;
               case 2:cent=LATTICE_I;centc='I';break;
               case 3:cent=LATTICE_A;centc='A';break;
            }
            minv=EstimateCellVolume(dmin,dmax,nb,MONOCLINIC,cent,1.5);
            maxv=EstimateCellVolume(dmin,dmax,nb,MONOCLINIC,cent,0.4*weak_f);
            mpCellExplorer->SetVolumeMinMax(minv,maxv);
            lengthmax=pow(maxv,(float)(1/3.0))*3;
            if(lengthmax<25)lengthmax=25;
            mpCellExplorer->SetLengthMinMax(3,lengthmax);
            mpCellExplorer->SetCrystalSystem(MONOCLINIC);
            mpCellExplorer->SetCrystalCentering(cent);
            mpLog->AppendText(wxString::Format(_T("MONOCLINIC %c : V= %6.0f -> %6.0f A^3, max length=%6.2fA"),centc,minv,maxv,lengthmax));
            t0=chrono.seconds();
            if(dlgProgress.Update(5,wxString::Format(_T("MONOCLINIC %c (%d spurious), V=%6.0f-%6.0f, l<%6.2fA\n")
                                                   _T("Best Score=%6.1f"),centc,
                                                   nbSpurious,minv,maxv,lengthmax,mpCellExplorer->GetBestScore()))==false) break;
            mpCellExplorer->DicVol(reportOnScore,reportOnDepth,stopOnScore,stopOnDepth);
            mpLog->AppendText(wxString::Format(_T(" -> %3u sols in %6.2fs, best score=%6.1f\n"),
                  mpCellExplorer->GetSolutions().size(),chrono.seconds()-t0,mpCellExplorer->GetBestScore()));
            mpLog->Update();
            if(noCentered) break;
         }
         
         nbSpurious+=1;
         if((mpCellExplorer->GetBestScore()>=stopOnScore)||(nbSpurious>3)) break;
         mpLog->AppendText(wxString::Format(_T("\n Trying now with %2u spurious peaks\n"),nbSpurious));
         mpLog->Update();
      }
   }
   else
   {
      // Erase spurious record
      for(vector<PeakList::hkl>::iterator pos=mpPeakList->mvHKL.begin();pos!=mpPeakList->mvHKL.end();++pos)
         pos->isSpurious=false;
      
      wxString s;
      double lmin,lmax,amin=90,amax,vmin,vmax,error,stopOnScore,reportOnScore;
      long nbspurious,nbPeak,stopOnDepth,reportOnDepth;
      s=mpLengthMin->GetValue();s.ToDouble(&lmin);
      s=mpLengthMax->GetValue();s.ToDouble(&lmax);
      //s=mpAngleMin->GetValue();s.ToDouble(&amin);
      s=mpAngleMax->GetValue();s.ToDouble(&amax);
      s=mpVolumeMin->GetValue();s.ToDouble(&vmin);
      s=mpVolumeMax->GetValue();s.ToDouble(&vmax);
      s=mpNbSpurious->GetValue();s.ToLong(&nbspurious);
      s=mpNbPeak->GetValue();s.ToLong(&nbPeak);
      s=mpErrorD->GetValue();s.ToDouble(&error);
      s=mpStopOnScore->GetValue();s.ToDouble(&stopOnScore);
      s=mpStopOnDepth->GetValue();s.ToLong(&stopOnDepth);
      s=mpReportOnScore->GetValue();s.ToDouble(&reportOnScore);
      s=mpReportOnDepth->GetValue();s.ToLong(&reportOnDepth);
      
      // Use at most 30 reflections for indexing
      if(mpPeakList->GetPeakList().size()>nbPeak) mpPeakList->GetPeakList().resize(nbPeak);
      
      mpCellExplorer = new CellExplorer(*mpPeakList,(CrystalSystem)(mpBravais->GetSelection()),0);
      
      mpCellExplorer->SetLengthMinMax((float)lmin,(float)lmax);
      mpCellExplorer->SetAngleMinMax((float)amin*DEG2RAD,(float)amax*DEG2RAD);
      mpCellExplorer->SetVolumeMinMax((float)vmin,(float)vmax);
      mpCellExplorer->SetNbSpurious((unsigned int)nbspurious);
      mpCellExplorer->SetD2Error((float)(error*error));

      mpCellExplorer->SetCrystalCentering(LATTICE_P);

      cout<<lmin<<" "<<lmax<<" "<<amin<<" "<<amax<<" "<<vmin<<" "<<vmax<<" "<<(unsigned int)nbspurious<<" "<<error*error<<endl;
      #if 1
      mpCellExplorer->DicVol(reportOnScore,reportOnDepth,stopOnScore,stopOnDepth);
      #else
      if(mpAlgorithm->GetSelection()==0) mpCellExplorer->DicVol(reportOnScore,reportOnDepth,stopOnScore,stopOnDepth);
      else
      {
         for(unsigned int i=0;i<20;++i)
         {
            mpCellExplorer->Evolution(5000,true,0.7,0.5,50);
            if(mpCellExplorer->GetBestScore()>stopOnScore) break;
         }
      }
      #endif
   }
   mpLog->AppendText(wxString::Format(_T("Finished indexing, bestscore=%6.1f, elapsed time=%6.2fs\n"),
                     mpCellExplorer->GetBestScore(),chrono.seconds()));
   // Merge similar solutions - useful here for solutions from different runs/different systems
   mpCellExplorer->ReduceSolutions();
   if(mpCellExplorer->GetSolutions().size()>0)
   {
      char buf[200];
      wxArrayString sols;
      float bestvol=0;
      for(list<pair<RecUnitCell,float> >::const_iterator pos=mpCellExplorer->GetSolutions().begin();
         pos!=mpCellExplorer->GetSolutions().end();++pos)
      {
         vector<float> uc=pos->first.DirectUnitCell();
         if(pos==mpCellExplorer->GetSolutions().begin()) bestvol=uc[6]*.99999;
         const float relvol=uc[6]/bestvol;
         string sys;
         switch(pos->first.mlattice)
         {
            case TRICLINIC:sys="TRICLINIC"; break;
            case MONOCLINIC:sys="MONOCLINIC"; break;
            case ORTHOROMBIC:sys="ORTHOROMBIC"; break;
            case HEXAGONAL:sys="HEXAGONAL"; break;
            case RHOMBOEDRAL:sys="RHOMBOEDRAL"; break;
            case TETRAGONAL:sys="TETRAGONAL"; break;
            case CUBIC:sys="CUBIC"; break;
         }
         char centc;
         switch(pos->first.mCentering)
         {
            case LATTICE_P:centc='P'; break;
            case LATTICE_I:centc='I'; break;
            case LATTICE_A:centc='A'; break;
            case LATTICE_B:centc='B'; break;
            case LATTICE_C:centc='C'; break;
            case LATTICE_F:centc='F'; break;
         }
         
         sprintf(buf,"Score=%6.1f V=%6.1f(%3.1fV) %6.3f %6.3f %6.3f %6.2f %6.2f %6.2f %s %c",pos->second,
               uc[6],relvol,uc[0],uc[1],uc[2],uc[3]*RAD2DEG,uc[4]*RAD2DEG,uc[5]*RAD2DEG,sys.c_str(),centc);
         //cout<<buf<<endl;
         sols.Add(wxString::FromAscii(buf));
      }
      mpCell->Set(sols);
   }
   
   if(mpGraph!=0) mpGraph->Refresh(FALSE);
}
void WXCellExplorer::OnSelectCell(wxCommandEvent &event)
{
   VFN_DEBUG_ENTRY("WXCellExplorer::OnSelectCell()",7)
   const int choice=mpCell->GetSelection();
   if(choice!=wxNOT_FOUND)
   {
      wxString s;
      long nbspurious;
      s=mpNbSpurious->GetValue();s.ToLong(&nbspurious);
      
      list<pair<RecUnitCell,float> >::const_iterator pos=mpCellExplorer->GetSolutions().begin();
      for(int i=0;i<choice;++i)++pos;// We need a random access ?
      // This will update the hkl in the list and therefore on the graph
      Score(*mpPeakList,pos->first,nbspurious,true,true,true);
      if(mpCrystal!=NULL)
      {
         // Apply crystal structure
         list<pair<RecUnitCell,float> >::const_iterator pos=mpCellExplorer->GetSolutions().begin();
         for(int i=0;i<choice;++i)++pos;// We need a random access ?
         vector<float> uc=pos->first.DirectUnitCell();
         mpCrystal->GetPar("a").SetValue(uc[0]);
         mpCrystal->GetPar("b").SetValue(uc[1]);
         mpCrystal->GetPar("c").SetValue(uc[2]);
         mpCrystal->GetPar("alpha").SetValue(uc[3]);
         mpCrystal->GetPar("beta").SetValue(uc[4]);
         mpCrystal->GetPar("gamma").SetValue(uc[5]);
         switch(pos->first.mlattice)
         {
            case TRICLINIC:mpCrystal->GetSpaceGroup().ChangeSpaceGroup("P-1");break;
            case MONOCLINIC:mpCrystal->GetSpaceGroup().ChangeSpaceGroup("P2/m");break;
            case ORTHOROMBIC:mpCrystal->GetSpaceGroup().ChangeSpaceGroup("Pmmm");break;
            case HEXAGONAL:mpCrystal->GetSpaceGroup().ChangeSpaceGroup("P6/mmm");break;
            case RHOMBOEDRAL:mpCrystal->GetSpaceGroup().ChangeSpaceGroup("R-3m");break;
            case TETRAGONAL:mpCrystal->GetSpaceGroup().ChangeSpaceGroup("P4/mmm");break;
            case CUBIC:mpCrystal->GetSpaceGroup().ChangeSpaceGroup("Pm-3m");break;
         }
         mpCrystal->UpdateDisplay();
      }
      try{
         if(mpAutomaticLeBail->GetValue())
         {
            VFN_DEBUG_MESSAGE("WXCellExplorer::OnSelectCell():auto-Le Bail",7);
            // run Le Bail
            const bool fitzero=true,
                     fitwidth0=true,
                     fitwidth=true,
                     fiteta=true,
                     fitasym=true,
                     fitdispltransp=true,
                     fitbackgd=true,
                     fitcell=true;
            
            wxProgressDialog dlgProgress(_T("Le Bail and Profile Fitting"),_T("Le Bail Fitting, cycle #0/20"),
                                       25,this,wxPD_AUTO_HIDE|wxPD_ELAPSED_TIME|wxPD_CAN_ABORT);
            mpDiff->SetExtractionMode(true,true);
            VFN_DEBUG_MESSAGE("WXCellExplorer::OnSelectCell():auto-Le Bail",7);
            
            mpLog->AppendText(wxString::Format(_T("Starting 20 Le Bail cycles\n")));
            for(int i=0;i<10;++i)
            {
               VFN_DEBUG_MESSAGE("WXCellExplorer::OnSelectCell():auto-Le Bail #"<<i,7);
               mpDiff->ExtractLeBail(2);
               mpDiff->GetParentPowderPattern().FitScaleFactorForRw();
               VFN_DEBUG_MESSAGE("WXCellExplorer::OnSelectCell():auto-Le Bail #"<<i,7);
               mpDiff->GetParentPowderPattern().UpdateDisplay();
               if(dlgProgress.Update(i,wxString::Format(_T("Le Bail Fitting, cycle #%d/20"),i*2))==false) return;
            }
            mpLog->AppendText(wxString::Format(_T("                  => Rwp=%5.3f%%, GoF=%7.3f\n"),
                                             mpDiff->GetParentPowderPattern().GetRw()*100,
                                             mpDiff->GetParentPowderPattern().GetChi2()
                                             /mpDiff->GetParentPowderPattern().GetNbPointUsed()));
            
            LSQNumObj lsqobj("Profile Fitting object");
            lsqobj.SetRefinedObj(mpDiff->GetParentPowderPattern());
            lsqobj.PrepareRefParList(true);
            lsqobj.SetParIsUsed(gpRefParTypeObjCryst,false);
            lsqobj.SetParIsUsed(gpRefParTypeScattDataScale,true);
            lsqobj.SetParIsUsed(gpRefParTypeScattDataProfile,true);
            lsqobj.SetParIsUsed(gpRefParTypeScattDataCorrPos,true);
            lsqobj.SetParIsUsed(gpRefParTypeScattDataBackground,true);
            lsqobj.SetParIsUsed(gpRefParTypeUnitCell,true);
            lsqobj.SetParIsFixed(gpRefParTypeObjCryst,true);
            lsqobj.SetParIsFixed(gpRefParTypeScattDataScale,false);
            
            // :TODO: take car of other profiles than pseudo-voigt (DE-PV)
            if(fitzero) lsqobj.SetParIsFixed("Zero",false);
            if(fitwidth0) lsqobj.SetParIsFixed("W",false);
            if(fitzero||fitwidth0)
            {
               mpLog->AppendText(wxString::Format(_T("Fitting zero shift && constant width\n")));
               if(dlgProgress.Update(11,_T("Fitting zero shift && constant width"))==false) return;
               lsqobj.Refine(5,true,false);
               mpDiff->GetParentPowderPattern().FitScaleFactorForRw();
               mpDiff->GetParentPowderPattern().UpdateDisplay();
               mpLog->AppendText(wxString::Format(_T("                  => Rwp=%6.3f%%, GoF=%7.3f\n"),
                                                mpDiff->GetParentPowderPattern().GetRw()*100,
                                                mpDiff->GetParentPowderPattern().GetChi2()
                                                /mpDiff->GetParentPowderPattern().GetNbPointUsed()));
            }
            if(fitwidth) lsqobj.SetParIsFixed("U",false);
            if(fitwidth) lsqobj.SetParIsFixed("V",false);
            if(fiteta) lsqobj.SetParIsFixed("Eta0",false);
            if(fitwidth||fiteta)
            {
               mpLog->AppendText(wxString::Format(_T("Fitting width and gaussian/lorentzian fixed mix\n")));
               if(dlgProgress.Update(12,_T("Fitting variable width and gaussian/lorentzian fixed mix"))==false) return;
               lsqobj.Refine(5,true,false);
               mpDiff->GetParentPowderPattern().FitScaleFactorForRw();
               mpDiff->GetParentPowderPattern().UpdateDisplay();
               mpLog->AppendText(wxString::Format(_T("                  => Rwp=%6.3f%%, GoF=%7.3f\n"),
                                                mpDiff->GetParentPowderPattern().GetRw()*100,
                                                mpDiff->GetParentPowderPattern().GetChi2()
                                                /mpDiff->GetParentPowderPattern().GetNbPointUsed()));
            }
            
            if(fiteta) lsqobj.SetParIsFixed("Eta1",false);
            if(fiteta)
            {
               mpLog->AppendText(wxString::Format(_T("Fitting variable width and gaussian/lorentzian mix\n")));
               if(dlgProgress.Update(13,_T("Fitting variable width and gaussian/lorentzian mix"))==false) return;
               lsqobj.Refine(5,true,false);
               mpDiff->GetParentPowderPattern().FitScaleFactorForRw();
               mpDiff->GetParentPowderPattern().UpdateDisplay();
               mpLog->AppendText(wxString::Format(_T("                  => Rwp=%6.3f%%, GoF=%7.3f\n"),
                                                mpDiff->GetParentPowderPattern().GetRw()*100,
                                                mpDiff->GetParentPowderPattern().GetChi2()
                                                /mpDiff->GetParentPowderPattern().GetNbPointUsed()));
            }
            
            if(fitasym) lsqobj.SetParIsFixed("Asym0",false);
            if(fitasym) lsqobj.SetParIsFixed("Asym1",false);
            if(fitasym) lsqobj.SetParIsFixed("Asym2",false);
            if(fitdispltransp) lsqobj.SetParIsFixed("2ThetaDispl",false);
            if(fitdispltransp) lsqobj.SetParIsFixed("2ThetaTransp",false);
            if(fitdispltransp||fitasym)
            {
               mpLog->AppendText(wxString::Format(_T("Fitting assymetry and sample displacement/transparency\n")));
               if(dlgProgress.Update(14,_T("Fitting assymetry and sample displacement/transparency"))==false) return;
               lsqobj.Refine(5,true,false);
               mpDiff->GetParentPowderPattern().FitScaleFactorForRw();
               mpDiff->GetParentPowderPattern().UpdateDisplay();
               mpLog->AppendText(wxString::Format(_T("                  => Rwp=%6.3f%%, GoF=%7.3f\n"),
                                                mpDiff->GetParentPowderPattern().GetRw()*100,
                                                mpDiff->GetParentPowderPattern().GetChi2()
                                                /mpDiff->GetParentPowderPattern().GetNbPointUsed()));
            }
            
            if(fitbackgd)
            {
               lsqobj.SetParIsFixed(gpRefParTypeScattDataBackground,false);
               // Make sure points beyond max resolution are not optimized
               const unsigned int nbcomp= mpGraph->GetWXPowderPattern().GetPowderPattern().GetNbPowderPatternComponent();
               for(unsigned int i=0;i<nbcomp;++i)
                  if(mpGraph->GetWXPowderPattern().GetPowderPattern().GetPowderPatternComponent(i).GetClassName()=="PowderPatternBackground")
                  {
                     PowderPatternBackground *pback=dynamic_cast<PowderPatternBackground *>
                        (&(mpGraph->GetWXPowderPattern().GetPowderPattern().GetPowderPatternComponent(i)));
                     pback->FixParametersBeyondMaxresolution(lsqobj.GetCompiledRefinedObj());
                  }
               
               mpLog->AppendText(wxString::Format(_T("Fitting background\n")));
               if(dlgProgress.Update(15,_T("Fitting background"))==false) return;
               lsqobj.Refine(5,true,false);
               mpDiff->GetParentPowderPattern().FitScaleFactorForRw();
               mpDiff->GetParentPowderPattern().UpdateDisplay();
               mpLog->AppendText(wxString::Format(_T("                  => Rwp=%6.3f%%, GoF=%7.3f\n"),
                                                mpDiff->GetParentPowderPattern().GetRw()*100,
                                                mpDiff->GetParentPowderPattern().GetChi2()
                                                /mpDiff->GetParentPowderPattern().GetNbPointUsed()));
            }
            
            if(fitcell) lsqobj.SetParIsFixed(gpRefParTypeUnitCell,false);
            if(fitcell)
            {
               mpLog->AppendText(wxString::Format(_T("Fitting unit cell\n")));
               if(dlgProgress.Update(16,_T("Fitting unit cell"))==false) return;
               lsqobj.Refine(5,true,false);
               mpDiff->GetParentPowderPattern().FitScaleFactorForRw();
               mpDiff->GetParentPowderPattern().UpdateDisplay();
               mpLog->AppendText(wxString::Format(_T("                  => Rwp=%6.3f%%, GoF=%7.3f\n"),
                                                mpDiff->GetParentPowderPattern().GetRw()*100,
                                                mpDiff->GetParentPowderPattern().GetChi2()
                                                /mpDiff->GetParentPowderPattern().GetNbPointUsed()));
            }
            // Run Le Bail again from scratch
            mpDiff->SetExtractionMode(true,true);
            mpLog->AppendText(wxString::Format(_T("Starting 10 Le Bail cycles\n")));
            for(int i=17;i<22;++i)
            {
               if(dlgProgress.Update(i,wxString::Format(_T("Le Bail Fitting, cycle #%d/10"),(i-17)*2))==false) return;
               mpDiff->ExtractLeBail(2);
               mpDiff->GetParentPowderPattern().FitScaleFactorForRw();
               mpDiff->GetParentPowderPattern().UpdateDisplay();
            }
            mpLog->AppendText(wxString::Format(_T("                  => Rwp=%5.3f%%, GoF=%7.3f\n"),
                                             mpDiff->GetParentPowderPattern().GetRw()*100,
                                             mpDiff->GetParentPowderPattern().GetChi2()
                                             /mpDiff->GetParentPowderPattern().GetNbPointUsed()));
            // Last fit
            mpLog->AppendText(wxString::Format(_T("Last fit...\n")));
            if(dlgProgress.Update(23,_T("Last fit..."))==false) return;
            lsqobj.Refine(5,true,false);
            mpDiff->GetParentPowderPattern().FitScaleFactorForRw();
            mpLog->AppendText(wxString::Format(_T("                  => Rwp=%5.3f%%, GoF=%7.3f\n"),
                                             mpDiff->GetParentPowderPattern().GetRw()*100,
                                             mpDiff->GetParentPowderPattern().GetChi2()
                                             /mpDiff->GetParentPowderPattern().GetNbPointUsed()));
            mpDiff->GetParentPowderPattern().UpdateDisplay();
            mpCrystal->UpdateDisplay();
         }
      }
      catch(const ObjCrystException &except)
      {
            mpLog->AppendText(wxString::Format(_T(" OOPS : refinement diverged ! Aborting.")));
      }
      if(mpGraph!=0) mpGraph->Refresh(FALSE);
      //:TODO: store refined cell parameters, display GoF in cell list
   }
   VFN_DEBUG_EXIT("WXCellExplorer::OnSelectCell",7)
}
void WXCellExplorer::OnApplyCell(wxCommandEvent &event)
{
   const int choice=mpCell->GetSelection();
   if((mpCrystal!=0)&&(choice!=wxNOT_FOUND))
   {
      list<pair<RecUnitCell,float> >::const_iterator pos=mpCellExplorer->GetSolutions().begin();
      for(int i=0;i<choice;++i)++pos;// We need a random access ?
      vector<float> uc=pos->first.DirectUnitCell();
      mpCrystal->GetPar("a").SetValue(uc[0]);
      mpCrystal->GetPar("b").SetValue(uc[1]);
      mpCrystal->GetPar("c").SetValue(uc[2]);
      mpCrystal->GetPar("alpha").SetValue(uc[3]);
      mpCrystal->GetPar("beta").SetValue(uc[4]);
      mpCrystal->GetPar("gamma").SetValue(uc[5]);
      mpCrystal->UpdateDisplay();
      if(mpGraph!=0)
      {
         wxCommandEvent ev(ID_POWDERGRAPH_MENU_UPDATE);
         mpGraph->OnUpdate(ev);
      }
   }
}

void WXCellExplorer::OnChooseCrystal(wxCommandEvent &event)
{
   VFN_DEBUG_MESSAGE("WXCellExplorer::OnChooseCrystal()",6)
   WXCrystValidateAllUserInput();
   int choice;
   mpCrystal=dynamic_cast<Crystal*>
      ( WXDialogChooseFromRegistry(gCrystalRegistry,(wxWindow*)this,
         "Choose a Crystal Structure:",choice));
   if(0==mpCrystal) return;
   mpFieldCrystal->SetValue(mpCrystal->GetName());
}

void WXCellExplorer::OnAutoLeBail(wxCommandEvent &event)
{
   if(mpGraph==NULL)
   {
      mpAutomaticLeBail->SetValue(false);
      return;
   }
   VFN_DEBUG_ENTRY("WXCellExplorer::OnAutoLeBail()",7)
   // Check if powder pattern has a background phase
   const unsigned int nbcomp= mpGraph->GetWXPowderPattern().GetPowderPattern().GetNbPowderPatternComponent();
   bool needBackground=true;
   for(unsigned int i=0;i<nbcomp;++i)
      if(mpGraph->GetWXPowderPattern().GetPowderPattern().GetPowderPatternComponent(i).GetClassName()=="PowderPatternBackground")
      {
         needBackground=false;
         break;
      };
   if(needBackground)
   {
      int answer =wxMessageBox(_T("To automatically run profile-fitting\n")
                               _T("and Le Bail extraction, you must have\n")
                               _T("defined a background phase for the pattern\n")
                               _T("and you will need to choose a crystal phase\n\n")
                               _T("Do you want to do that now ?"),
                               _T("Add Background ?"),wxYES_NO|wxICON_QUESTION);
      if(answer==wxNO)
      {
         mpAutomaticLeBail->SetValue(false);
         return;
      }
      wxCommandEvent ev;
      mpGraph->GetWXPowderPattern().OnMenuAddCompBackgdBayesian(ev);
   }
   // Check if a crystal structure has been selected to apply the calculated cell
   if(mpCrystal==NULL)
   {
      if(gCrystalRegistry.GetNb()==0)
      {
         mpCrystal=new Crystal(4,5,6,"P1");
         mpCrystal->SetName("Indexing Result");
      }
      else
      {
         int answer =wxMessageBox(_T("To automatically run profile-fitting\n")
                                 _T("and Le Bail extraction, you must have\n")
                                 _T("defined a Crystal phase to apply the cell to\n\n")
                                 _T("Do you want to use an EXISTING crystal structure ?\n")
                                 _T("(otherwise a new one will be created for you)"),
                                 _T("Select Crystal ?"),wxYES_NO|wxICON_QUESTION);
         if(answer==wxNO)
         {
            mpCrystal=new Crystal(4,5,6,"P1");
            mpCrystal->SetName("Indexing Result");
         }
         else
         {
            wxCommandEvent ev;
            this->OnChooseCrystal(ev);
         }
      }
   }
   // Now make sure this Crystal structure is used by the powder pattern object
   bool needPowderPatternDiffraction=true;
   unsigned int nbPowderPatternDiffraction=0;
   for(unsigned int i=0;i<nbcomp;++i)
      if(mpGraph->GetWXPowderPattern().GetPowderPattern().GetPowderPatternComponent(i).GetClassName()=="PowderPatternDiffraction")
      {
         nbPowderPatternDiffraction++;
         mpDiff=dynamic_cast<PowderPatternDiffraction*>
            (&(mpGraph->GetWXPowderPattern().GetPowderPattern().GetPowderPatternComponent(i)));
         if(&(mpDiff->GetCrystal())==mpCrystal)
         {
            needPowderPatternDiffraction=false;
            break;
         }
      };
   VFN_DEBUG_MESSAGE("WXCellExplorer::OnAutoLeBail():needPowderPatternDiffraction=="<<needPowderPatternDiffraction,7)
   if(needPowderPatternDiffraction)
   {
      if(nbPowderPatternDiffraction>0)
      {
         int answer =wxMessageBox(_T("To automatically run profile-fitting\n")
                                  _T("and Le Bail extraction, you must assign\n")
                                  _T("the Crystal to a Diffraction Component\n\n")
                                  _T("Do you want to use an already existing one ?"),
                                  _T("Use Crystal Phase ?"),wxYES_NO|wxICON_QUESTION);
         if(answer==wxYES)
         {
            //if(nbPowderPatternDiffraction==1) :TODO: handle multiple phases
            {
               for(unsigned int i=0;i<nbcomp;++i)
                  if(mpGraph->GetWXPowderPattern().GetPowderPattern().GetPowderPatternComponent(i).GetClassName()=="PowderPatternDiffraction")
                  {
                     mpDiff=dynamic_cast<PowderPatternDiffraction*>
                        (&(mpGraph->GetWXPowderPattern().GetPowderPattern().GetPowderPatternComponent(i)));
                     mpDiff->SetCrystal(*mpCrystal);
                     return;
                  }
            }
         }
      }
      else
      {// Create one crystalline phase
         VFN_DEBUG_MESSAGE("WXCellExplorer::OnAutoLeBail():Create PowderPatternDiffraction",7)
         mpDiff=new PowderPatternDiffraction;
         mpDiff->SetCrystal(*mpCrystal);
         mpGraph->GetWXPowderPattern().GetPowderPattern().AddPowderPatternComponent(*mpDiff);
         if(mpGraph->GetWXPowderPattern().GetPowderPattern().GetRadiation().GetWavelengthType()==WAVELENGTH_TOF)
         {
            wxCommandEvent event(wxEVT_COMMAND_MENU_SELECTED,ID_POWDERDIFF_PROFILE_DEPV);
            wxPostEvent(mpDiff->WXGet(),event);
         }
         mpGraph->GetWXPowderPattern().GetPowderPattern().Prepare();
         mpGraph->GetWXPowderPattern().CrystUpdate();
      }
   }
   // Limit resolution (:TODO: Take into account density of peask to limit to ~100 reflections)
   mpDiff->GetParentPowderPattern().SetMaxSinThetaOvLambda(0.25);
   // If one cell is already selected, do optimization immediately
   if(mpCell->GetSelection()>=0)
   {
      VFN_DEBUG_MESSAGE("WXCellExplorer::OnAutoLeBail()->OnSelectCell()",7)
      cout<<mpCell->GetSelection()<<endl;
      wxCommandEvent ev;
      this->OnSelectCell(ev);
   }
   VFN_DEBUG_EXIT("WXCellExplorer::OnAutoLeBail()",7)
}

//////////////////////////////////////// END WXCellExplorer /////////////////////

void WXPowderPatternGraph::OnIndex(wxCommandEvent& WXUNUSED(event))
{
   wxFrame *mpFrame=new wxFrame(this,-1,_T("Fox cell Explorer (EXPERIMENTAL)"));
   WXCellExplorer *mpWXCellExplorer;
   mpWXCellExplorer=new WXCellExplorer(mpFrame,mPeakList,this);
   mpFrame->Show(TRUE);
}

void WXPowderPatternGraph::OnChangeScale(wxCommandEvent& event)
{
   VFN_DEBUG_MESSAGE("WXPowderPatternGraph::OnChangeScale()",10)
   if(event.GetId()==ID_POWDERGRAPH_MENU_XSCALE_DATA)  mXScale=0;
   if(event.GetId()==ID_POWDERGRAPH_MENU_XSCALE_D)     mXScale=1;
   if(event.GetId()==ID_POWDERGRAPH_MENU_XSCALE_2PID)  mXScale=2;
   if(event.GetId()==ID_POWDERGRAPH_MENU_YSCALE_LINEAR)mYScale=0;
   if(event.GetId()==ID_POWDERGRAPH_MENU_YSCALE_SQRT)  mYScale=1;
   if(event.GetId()==ID_POWDERGRAPH_MENU_YSCALE_LOG10) mYScale=2;

   mpPopUpMenu->Enable(ID_POWDERGRAPH_MENU_XSCALE_DATA, TRUE);
   mpPopUpMenu->Enable(ID_POWDERGRAPH_MENU_XSCALE_D, TRUE);
   mpPopUpMenu->Enable(ID_POWDERGRAPH_MENU_XSCALE_2PID, TRUE);
   mpPopUpMenu->Enable(ID_POWDERGRAPH_MENU_YSCALE_LINEAR, TRUE);
   mpPopUpMenu->Enable(ID_POWDERGRAPH_MENU_YSCALE_SQRT, TRUE);
   mpPopUpMenu->Enable(ID_POWDERGRAPH_MENU_YSCALE_LOG10, TRUE);
   
   if(mXScale==0)  mpPopUpMenu->Enable(ID_POWDERGRAPH_MENU_XSCALE_DATA, FALSE);
   if(mXScale==1)     mpPopUpMenu->Enable(ID_POWDERGRAPH_MENU_XSCALE_D, FALSE);
   if(mXScale==2)  mpPopUpMenu->Enable(ID_POWDERGRAPH_MENU_XSCALE_2PID, FALSE);
   if(mYScale==0)mpPopUpMenu->Enable(ID_POWDERGRAPH_MENU_YSCALE_LINEAR, FALSE);
   if(mYScale==1)  mpPopUpMenu->Enable(ID_POWDERGRAPH_MENU_YSCALE_SQRT, FALSE);
   if(mYScale==2) mpPopUpMenu->Enable(ID_POWDERGRAPH_MENU_YSCALE_LOG10, FALSE);

   this->Refresh(false);
}

void WXPowderPatternGraph::OnLeBail(wxCommandEvent& event)
{
   wxFrame *pFrame=new wxFrame(this,-1,_T("Profile Fitting"));
   WXProfileFitting *pFit;
   pFit=new WXProfileFitting(pFrame,&(this->GetWXPowderPattern().GetPowderPattern()));
   pFrame->Show(true);
}

void WXPowderPatternGraph::OnKeyDown(wxKeyEvent& event)
{
   wxMutexLocker mlock(mMutex);
   const long nbPoint=mX.numElements();
   
   switch(event.GetKeyCode())
   {
      case(WXK_LEFT):
      {
         const REAL range=mMaxX-mMinX;
         mMinX -= range/8;
         if(mX(nbPoint-1)>mX(0))
         {
            if(mMinX<mX(0)) mMinX=mX(0);
         }
         else 
         {
            if(mMinX<mX(nbPoint-1)) mMinX=mX(nbPoint-1);
         }
         mMaxX=mMinX+range;
         break;
      }
      case(WXK_RIGHT):
      {
         const REAL range=mMaxX-mMinX;
         mMaxX += range/8;
         if(mX(nbPoint-1)>mX(0))
         {
            if(mMaxX>=mX(nbPoint-1)) mMaxX=mX(nbPoint-1);
         }
         else
         {
            if(mMaxX>=mX(0)) mMaxX=mX(0);
         }
         mMinX=mMaxX-range;
         break;
      }
      case(WXK_UP):
      {
         REAL max=mObs.max(),min=mObs.min();
         if(min<1e-6*max)min=1e-6*max;
         if(mYScale==2)
         {
            const REAL range=log10(max)-log10(min);
            mMinIntensity*=pow(10,range/8);
            mMaxIntensity*=pow(10,range/8);
            break;
         }
         const REAL range=mMaxIntensity-mMinIntensity;
         mMinIntensity+=range/8;
         mMaxIntensity+=range/8;
         break;
      }
      case(WXK_DOWN):
      {
         REAL max=mObs.max(),min=mObs.min();
         if(min<1e-6*max)min=1e-6*max;
         if(mYScale==2)
         {
            const REAL range=log10(max)-log10(min);
            mMinIntensity*=pow(10,-range/8);
            mMaxIntensity*=pow(10,-range/8);
            break;
         }
         const REAL range=mMaxIntensity-mMinIntensity;
         mMinIntensity-=range/8;
         if(mMinIntensity<1e-6*max) mMinIntensity=1e-6*max;
         mMaxIntensity=mMinIntensity+range;
         break;
      }
      case(43):// WXK_ADD ?
      {
         const REAL halfrange=(mMaxX-mMinX)/2;
         const REAL middle=(mMaxX+mMinX)/2;
         mMinX= (long)(middle-halfrange*4./5.);
         mMaxX = (long)(middle+halfrange*4./5.);
         break;
      }
      case(45):// WXK_SUBTRACT ?
      {
         const REAL halfrange=(mMaxX-mMinX)/2;
         const REAL middle=(mMaxX+mMinX)/2;
         mMinX= (long)(middle-halfrange*5./4.);
         mMaxX = (long)(middle+halfrange*5./4.);
         if(mX(nbPoint-1)>mX(0))
         {
            if(mMinX<mX(0)) mMinX=mX(0);
            if(mMaxX>mX(nbPoint-1)) mMaxX=mX(nbPoint-1);
         }
         else
         {
            if(mMinX<mX(nbPoint-1)) mMinX=mX(nbPoint-1);
            if(mMaxX>mX(0)) mMaxX=mX(0);
         }
         break;
      }
      case(42):// WXK_MULTIPLY
      {
         const REAL range=mMaxIntensity-mMinIntensity;
         mMaxIntensity=mMinIntensity+range*4./5.;
         break;
      }
      case(47):// WXK_DIVIDE
      {
         const REAL range=mMaxIntensity-mMinIntensity;
         mMaxIntensity=mMinIntensity+range*5./4.;
         break;
      }
      default: 
      {
         VFN_DEBUG_MESSAGE("WXPowderPatternGraph::OnKeyDown(): no command for key #"<<event.GetKeyCode(),5);
         cout<<"WXPowderPatternGraph::OnKeyDown(): no command for key #"<<event.GetKeyCode()<<endl;
      }
   }
   mClockAxisLimits.Click();
   wxUpdateUIEvent ev(ID_POWDER_GRAPH_NEW_PATTERN);
   wxPostEvent(this,ev);
   event.Skip();
}

void WXPowderPatternGraph::OnSize(wxSizeEvent& event)
{
   this->Refresh(false);
}

WXPowderPattern& WXPowderPatternGraph::GetWXPowderPattern(){return *mpPattern;}
const WXPowderPattern& WXPowderPatternGraph::GetWXPowderPattern()const{return *mpPattern;}

void WXPowderPatternGraph::SetPattern(const CrystVector_REAL &x,
                                      const CrystVector_REAL &obs,
                                      const CrystVector_REAL &calc,
                                      const CrystVector_REAL &sigma,
                                      const CrystVector_REAL &chi2Cumul)
{
   VFN_DEBUG_ENTRY("WXPowderPatternGraph::SetPattern(x,obs,calc,sigma)",4)
   mMutex.Lock();
   mX=x;
   if(mpPattern->GetPowderPattern().GetRadiation().GetWavelengthType()!=WAVELENGTH_TOF) mX*=RAD2DEG;
   mCalc=calc;
   //mCalc=SavitzkyGolay(obs,4,2);
   //mCalc *= -obs.max()/mCalc.max();
   mObs=obs;
   mSigma=sigma;
   mChi2Cumul=chi2Cumul;
   // Reset the zoom parameters, only for the first display or if the limits of the
   // full pattern have changed
   if(  (mMaxX<0)
      ||(mpPattern->GetPowderPattern().GetClockPowderPatternPar()>mClockAxisLimits)) 
   {
      mMutex.Unlock();
      this->ResetAxisLimits();
      mMutex.Lock();
   }
   
   mvLabelList.clear();
   for(unsigned int i=0;i<mpPattern->GetPowderPattern().GetNbPowderPatternComponent();++i)
      mvLabelList.push_back(mpPattern->GetPowderPattern()
                              .GetPowderPatternComponent(i).GetPatternLabelList());

   mMutex.Unlock();
   // If we only send an OnPaint event, only the parts which have been erased are redrawn
   // (under windows). SO we must force the complete Refresh of the window... in the
   // main thread of course...
   if(true==wxThread::IsMain())
   {
      this->Refresh(false);
   }
   else
   {
      wxUpdateUIEvent event(ID_POWDER_GRAPH_NEW_PATTERN);
      wxPostEvent(this,event);
   }
   //cout<<FormatVertVector<REAL>(x,obs,calc,sigma)<<endl;
   VFN_DEBUG_EXIT("WXPowderPatternGraph::SetPattern(x,obs,calc,sigma)"<<mX.numElements()<<","<<mCalc.numElements()<<","<<mObs.numElements()<<","<<mSigma.numElements()<<",",4)
}

void WXPowderPatternGraph::OnRedrawNewPattern(wxUpdateUIEvent& WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXPowderPatternGraph::SetPattern()",5)
   this->Refresh(false);
}

void WXPowderPatternGraph::ResetAxisLimits()
{
   wxMutexLocker mlock(mMutex);
   mMaxIntensity=mObs.max();
   mMinIntensity=mObs.min();
   const float max=mObs.max();
   const float min=mObs.min();
   if(max>mMaxIntensity) mMaxIntensity=max;
   if(min<mMinIntensity) mMinIntensity=min;
   if(mMinIntensity<=0) mMinIntensity=max/1e6;
   mMaxX=mX.max();
   mMinX=mX.min();
   mClockAxisLimits.Click();
   VFN_DEBUG_MESSAGE("WXPowderPatternGraph::ResetAxisLimits():"<<mMinIntensity<<","<<mMaxIntensity<<","<<mMinX<<","<<mMaxX,10)
}
long WXPowderPatternGraph::Data2ScreenX(const REAL x)const
{
   wxCoord width,height;
   this->GetSize(&width, &height);
   REAL xs=x,minx=mMinX,maxx=mMaxX;
   if(xs<minx)xs=minx;
   if(xs>maxx)xs=maxx;
   float d,mind,maxd;
   if(mpPattern->GetPowderPattern().GetRadiation().GetWavelengthType()!=WAVELENGTH_TOF)
   {
      d=2*mpPattern->GetPowderPattern().X2STOL(xs*DEG2RAD);
      mind=2*mpPattern->GetPowderPattern().X2STOL(minx*DEG2RAD);
      maxd=2*mpPattern->GetPowderPattern().X2STOL(maxx*DEG2RAD);
   }
   else
   {
      d=2*mpPattern->GetPowderPattern().X2STOL(xs);
      mind=2*mpPattern->GetPowderPattern().X2STOL(minx);
      maxd=2*mpPattern->GetPowderPattern().X2STOL(maxx);
   }
   if(mXScale==1) {xs=d;minx=mind;maxx=maxd;}
   if(mXScale==2) {xs=2*M_PI*d;minx=2*M_PI*mind;maxx=2*M_PI*maxd;}
   return (long)(mMargin*3+(xs-minx)*(width-3*mMargin)/(maxx-minx));
}
long WXPowderPatternGraph::Point2ScreenX(const long x)const
{
   return this->Data2ScreenX(mX(x));
}
long WXPowderPatternGraph::Data2ScreenY(const REAL y)const
{
   wxCoord width,height;
   this->GetSize(&width, &height);
   REAL ys=y,miny=mMinIntensity,maxy=mMaxIntensity;
   if(ys<miny)ys=miny;
   if(ys>maxy)ys=maxy;
   if(mYScale==1) {ys=sqrt(ys) ;miny=sqrt(miny) ;maxy=sqrt(maxy);}
   if(mYScale==2) {ys=log10(ys);miny=log10(miny);maxy=log10(maxy);}
   return (long)(height-mMargin-(ys-miny)*(height-2*mMargin)/(maxy-miny));
}
REAL WXPowderPatternGraph::Screen2DataX(const long x)const
{
   wxCoord width,height;
   this->GetSize(&width, &height);
   REAL minx=mMinX,maxx=mMaxX;
   float mind,maxd;
   if(mpPattern->GetPowderPattern().GetRadiation().GetWavelengthType()!=WAVELENGTH_TOF)
   {
      mind=2*mpPattern->GetPowderPattern().X2STOL(minx*DEG2RAD);
      maxd=2*mpPattern->GetPowderPattern().X2STOL(maxx*DEG2RAD);
   }
   else
   {
      mind=2*mpPattern->GetPowderPattern().X2STOL(minx);
      maxd=2*mpPattern->GetPowderPattern().X2STOL(maxx);
   }
   if(mXScale==1)
   {
      minx=mind;
      maxx=maxd;
      REAL stol=(minx+(x-mMargin*3)*(maxx-minx)/(REAL)(width-3*mMargin))/2;
      if(mpPattern->GetPowderPattern().GetRadiation().GetWavelengthType()!=WAVELENGTH_TOF)
         return mpPattern->GetPowderPattern().STOL2X(stol)*RAD2DEG;
      else
         return mpPattern->GetPowderPattern().STOL2X(stol);
   }
   if(mXScale==2)
   {
      minx=2*M_PI*mind;
      maxx=2*M_PI*maxd;
      REAL stol=(minx+(x-mMargin*3)*(maxx-minx)/(REAL)(width-3*mMargin))/(4*M_PI);
      if(mpPattern->GetPowderPattern().GetRadiation().GetWavelengthType()!=WAVELENGTH_TOF)
         return mpPattern->GetPowderPattern().STOL2X(stol)*RAD2DEG;
      else
         return mpPattern->GetPowderPattern().STOL2X(stol);
   }
   return mMinX+(x-mMargin*3)*(mMaxX-mMinX)/(REAL)(width-3*mMargin);
}
REAL WXPowderPatternGraph::Screen2DataY(const long y)const
{
   wxCoord width,height;
   this->GetSize(&width, &height);
   REAL miny=mMinIntensity,maxy=mMaxIntensity;
   if(mYScale==1) {miny=sqrt(miny) ;maxy=sqrt(maxy);}
   if(mYScale==2) {miny=log10(miny);maxy=log10(maxy);}
   REAL ys=miny+(height-mMargin-y)*(maxy-miny)/(REAL)(height-2*mMargin);
   if(mYScale==1) ys=ys*ys;
   if(mYScale==2) ys=pow((float)10,(float)ys);
   return ys;
}

////////////////////////////////////////////////////////////////////////
//
//    WXPowderPatternBackgound
//
////////////////////////////////////////////////////////////////////////
static const long ID_POWDERBACKGROUND_GRID= WXCRYST_ID(); 
static const long ID_POWDERBACKGROUND_NEWBAYESIAN= WXCRYST_ID(); 

BEGIN_EVENT_TABLE(WXPowderPatternBackground, wxWindow)
   EVT_MENU(ID_POWDERBACKGROUND_IMPORT, 
                     WXPowderPatternBackground::OnMenuImportUserBackground)
   EVT_MENU(ID_POWDERBACKGROUND_OPTIMIZEBAYESIAN, 
                     WXPowderPatternBackground::OnMenuOptimizeBayesianBackground)
   EVT_GRID_CMD_CELL_CHANGE(ID_POWDERBACKGROUND_GRID,
                     WXPowderPatternBackground::OnEditGridBackgroundPoint)
   EVT_MENU(ID_POWDERBACKGROUND_NEWBAYESIAN,
                     WXPowderPatternBackground::OnMenuAutomaticBayesianBackground)
END_EVENT_TABLE()

WXPowderPatternBackground::WXPowderPatternBackground(wxWindow *parent, 
                                                     PowderPatternBackground *b):
WXRefinableObj(parent,b),mpPowderPatternBackground(b),mNeedUpdateUI(false),mIsSelfUpdating(false)
{
   mpWXTitle->SetForegroundColour(wxColour(0,255,0));
   //Menu
      mpMenuBar->AddMenu("Object",ID_REFOBJ_MENU_OBJ);
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_POWDERBACKGROUND_IMPORT,"Import");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_POWDERBACKGROUND_OPTIMIZEBAYESIAN,
         "Bayesian Optimization");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_POWDERBACKGROUND_NEWBAYESIAN,
         "New Automatic Background (Change Number of Points)");
   VFN_DEBUG_MESSAGE(mpMenuBar->GetSize().GetWidth()<<","<<mpMenuBar->GetSize().GetHeight(),10);
   mpSizer->SetItemMinSize(mpMenuBar,
                           mpMenuBar->GetSize().GetWidth(),
                           mpMenuBar->GetSize().GetHeight());
   
   #ifdef USE_BACKGROUND_MAXLIKE_ERROR
   WXCrystObjBasic* pFieldModelSigma=mpPowderPatternBackground
                                       ->GetPar("ML Model Error").wxCreate(this);
   mpSizer->Add(pFieldModelSigma,0,wxALIGN_LEFT);
   mList.Add(pFieldModelSigma);
   #endif
   // List of background points
      wxGridCellAttr* cellAttrFloat = new wxGridCellAttr;
      cellAttrFloat->SetRenderer(new wxGridCellFloatRenderer(10,3));
      cellAttrFloat->SetEditor(new wxGridCellFloatEditor(10,3));

      mpGridBackgroundPoint= new wxGrid(this,ID_POWDERBACKGROUND_GRID);
      mpGridBackgroundPoint->SetSize(400,300);
      mpGridBackgroundPoint->EnableScrolling(true,true);
      mpGridBackgroundPoint->SetSizeHints(-1,300,-1,300);
      mpGridBackgroundPoint->SetColSize(0,150);
      mpGridBackgroundPoint->CreateGrid(0,2);
      mpGridBackgroundPoint->SetColAttr(0,cellAttrFloat);
      mpGridBackgroundPoint->SetColAttr(1,cellAttrFloat);
      mpGridBackgroundPoint->SetColLabelValue(0,_T("Position"));
      mpGridBackgroundPoint->SetColLabelValue(1,_T("Intensity"));
      mpGridBackgroundPoint->AutoSizeRows();
      mpSizer->Add(mpGridBackgroundPoint,0,wxALIGN_LEFT);
   mpTopSizer->SetSizeHints(this);
   this->Layout();
   this->CrystUpdate(true);
}
void WXPowderPatternBackground::OnMenuImportUserBackground(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXPowderPatternBackground::OnMenuImportUserBackground()",6)
   wxFileDialog *open= new wxFileDialog(this,_T("Choose background file with 2Theta Ibackgd"),
                                        _T(""),_T(""),_T("*.*"),wxOPEN | wxFILE_MUST_EXIST);
   if(open->ShowModal() != wxID_OK) return;
   mpPowderPatternBackground->ImportUserBackground(string(open->GetPath().ToAscii()));
   open->Destroy();
}
void WXPowderPatternBackground::OnMenuOptimizeBayesianBackground(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_ENTRY("WXPowderPatternBackground::OnMenuOptimizeBayesianBackground()",6)
   mpPowderPatternBackground->UnFixAllPar();
   mpPowderPatternBackground->OptimizeBayesianBackground();
   mpPowderPatternBackground->FixAllPar();
   this->CrystUpdate();
   VFN_DEBUG_EXIT("WXPowderPatternBackground::OnMenuOptimizeBayesianBackground()",6)
}
void WXPowderPatternBackground::OnMenuAutomaticBayesianBackground(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_ENTRY("WXPowderPatternBackground::OnMenuAutomaticBayesianBackground()",6)
   WXCrystValidateAllUserInput();

   long nbPointSpline=20;
   wxString mes(_T("Number of Interpolation Points"));
   wxString s;
   s.Printf(_T("%d"),nbPointSpline);
   wxTextEntryDialog dialog(this,mes,_T("Automatic Bayesian (David-Sivia) Background"),
                            s,wxOK | wxCANCEL);
   if(wxID_OK!=dialog.ShowModal())
   {
      VFN_DEBUG_EXIT("WXPowderPatternBackground::OnMenuAutomaticBayesianBackground():Canceled",6)
      return;
   }
   dialog.GetValue().ToLong(&nbPointSpline);
   wxProgressDialog dlgProgress(_T("Automatic Bayesian Background"),_T("Automatic Background, Initializing..."),
                                      4,this,wxPD_AUTO_HIDE|wxPD_ELAPSED_TIME|wxPD_CAN_ABORT);
   if(nbPointSpline<2) nbPointSpline=2;
   {
      CrystVector_REAL x(nbPointSpline),backgd(nbPointSpline);
      const CrystVector_REAL *pObs=&(mpPowderPatternBackground->GetParentPowderPattern().GetPowderPatternObs());
      const unsigned long nbPoint=mpPowderPatternBackground->GetParentPowderPattern().GetNbPoint();
      const float xmin=mpPowderPatternBackground->GetParentPowderPattern()
                       .GetPowderPatternX()(0),
                  xmax=mpPowderPatternBackground->GetParentPowderPattern()
                       .GetPowderPatternX()(nbPoint-1);
      for(int i=0;i<nbPointSpline;i++)
      {// xmax is not necessarily > xmin, but in the right order (TOF)
         x(i)=xmin+(xmax-xmin)/(REAL)(nbPointSpline-1)*REAL(i);
         REAL x1=xmin+(xmax-xmin)/(REAL)(nbPointSpline-1)*REAL(i-.2);
         REAL x2=xmin+(xmax-xmin)/(REAL)(nbPointSpline-1)*REAL(i+.2);
         long n1=(long)(mpPowderPatternBackground->GetParentPowderPattern().X2Pixel(x1));
         long n2=(long)(mpPowderPatternBackground->GetParentPowderPattern().X2Pixel(x2));
         if(n1<0) n1=0;
         if(n2>(long)nbPoint)n2=nbPoint;
         backgd(i)=(*pObs)(n1);
         for(long j=n1;j<n2;j++)
            if((*pObs)(j)<backgd(i))backgd(i)=(*pObs)(j);
      }
      mpPowderPatternBackground->SetInterpPoints(x,backgd);
   }
   //mpPowderPatternBackground->GetParentPowderPattern().Prepare();
   mpPowderPatternBackground->UnFixAllPar();
   mpPowderPatternBackground->GetOption(0).SetChoice(0);//linear
   if(dlgProgress.Update(1,_T("Automatic Background: Optimizing Linear Model..."))==false) return;
   mpPowderPatternBackground->OptimizeBayesianBackground();
   mpPowderPatternBackground->GetOption(0).SetChoice(1);//spline
   if(dlgProgress.Update(2,_T("Automatic Background: Optimizing Spline Model..."))==false) return;
   mpPowderPatternBackground->OptimizeBayesianBackground();
   mpPowderPatternBackground->FixAllPar();

   this->CrystUpdate();
   VFN_DEBUG_EXIT("WXPowderPatternBackground::OnMenuAutomaticBayesianBackground()",6)
}
void WXPowderPatternBackground::OnEditGridBackgroundPoint(wxGridEvent &e)
{
   if(mIsSelfUpdating) return;
   VFN_DEBUG_ENTRY("WXPowderPatternBackground::OnEditGridBackgroundPoint():"<<e.GetRow()<<","<<e.GetCol(),10)
   const long r=e.GetRow();
   const long c=e.GetCol();
   wxString s=mpGridBackgroundPoint->GetCellValue(r,c);
   if(s!=_T(""))
   {
      REAL f=1.0;
      if(mpPowderPatternBackground->GetParentPowderPattern().GetRadiation().GetWavelengthType()!=WAVELENGTH_TOF)f=DEG2RAD;
      double d;
      s.ToDouble(&d);
      if(c==0)
      {
         if(d!=mBackgroundInterpPointX(r))
            mBackgroundInterpPointX(r)=d*f;
      }
      else
      {
         if(d!=mBackgroundInterpPointX(r))
            mBackgroundInterpPointIntensity(r)=d;
      }
      
      mpPowderPatternBackground->SetInterpPoints(mBackgroundInterpPointX,
                                                 mBackgroundInterpPointIntensity);
      // The order of the points might have changed
      mBackgroundInterpPointX        =*(mpPowderPatternBackground->GetInterpPoints().first);
      mBackgroundInterpPointIntensity=*(mpPowderPatternBackground->GetInterpPoints().second);
   }
   mNeedUpdateUI=true,
   this->UpdateUI();
   mpPowderPatternBackground->GetParentPowderPattern().UpdateDisplay();
   VFN_DEBUG_EXIT("WXPowderPatternBackground::OnEditGridBackgroundPoint():"<<e.GetRow()<<","<<e.GetCol(),10)
}

void WXPowderPatternBackground::CrystUpdate(const bool uui,const bool lock)
{
   if(lock) mMutex.Lock();
   this->WXRefinableObj::CrystUpdate(uui,false);
   if(false==mpPowderPatternBackground->IsBeingRefined())
   {
      const long diff=mpPowderPatternBackground->GetInterpPoints().first->numElements()
                      -mpGridBackgroundPoint->GetNumberRows();
      if(diff>0)
      {
         mNeedUpdateUI=true;
         mpGridBackgroundPoint->AppendRows(diff);
      }
      if(diff<0)
      {
         mNeedUpdateUI=true;
         mpGridBackgroundPoint->DeleteRows(0,-diff);
      }
      if(diff==0)
         if(  (MaxDifference(mBackgroundInterpPointX        ,
                             *(mpPowderPatternBackground->GetInterpPoints().first )))
            ||(MaxDifference(mBackgroundInterpPointIntensity,
                             *(mpPowderPatternBackground->GetInterpPoints().second))))
            mNeedUpdateUI=true;
      if(mNeedUpdateUI)
      {
         mBackgroundInterpPointX        =*(mpPowderPatternBackground->GetInterpPoints().first);
         mBackgroundInterpPointIntensity=*(mpPowderPatternBackground->GetInterpPoints().second);
      }
   }
   if(lock) mMutex.Unlock();
}

void WXPowderPatternBackground::UpdateUI(const bool lock)
{
   if(lock) mMutex.Lock();
   if(mNeedUpdateUI)
   {
      REAL f=1.0;
      if(mpPowderPatternBackground->GetParentPowderPattern().GetRadiation().GetWavelengthType()!=WAVELENGTH_TOF)f=RAD2DEG;
      const long nb=mBackgroundInterpPointX.numElements();
      mIsSelfUpdating=true;
      for(long i=0;i<nb;++i)
      {
         wxString tmp;
         tmp.Printf(_T("%f"),f*mBackgroundInterpPointX(i));
         mpGridBackgroundPoint->SetCellValue(i,0,tmp);
         tmp.Printf(_T("%f"),mBackgroundInterpPointIntensity(i));
         mpGridBackgroundPoint->SetCellValue(i,1,tmp);
      }
      mIsSelfUpdating=false;
   }
      
   mNeedUpdateUI=false;
   this->WXRefinableObj::UpdateUI(false);
   if(lock) mMutex.Unlock();
}
////////////////////////////////////////////////////////////////////////
//
//    WXTexturePhaseMarchDollase
//
////////////////////////////////////////////////////////////////////////
WXTexturePhaseMarchDollase::WXTexturePhaseMarchDollase(wxWindow *parent, 
                                                       TexturePhaseMarchDollase *pObj,
                                                       TextureMarchDollase* pTex):
WXCrystObjBasic(parent),mpTexturePhaseMarchDollase(pObj)
{
   VFN_DEBUG_ENTRY("WXTexturePhaseMarchDollase::WXTexturePhaseMarchDollase()",5)
   mpSizer=new wxBoxSizer(wxHORIZONTAL);
   this->SetSizer(mpSizer);
   pTex->Print();
   WXCrystObjBasic* pFieldFraction=pTex->GetPar(&(pObj->mFraction)).WXCreate(this);
   mpSizer->Add(pFieldFraction,0,wxALIGN_LEFT);
   mList.Add(pFieldFraction);

   WXCrystObjBasic* pFieldMarch=pTex->GetPar(&(pObj->mMarchCoeff)).WXCreate(this);
   mpSizer->Add(pFieldMarch,0,wxALIGN_LEFT);
   mList.Add(pFieldMarch);

   WXCrystObjBasic* pFieldH=pTex->GetPar(&(pObj->mH)).WXCreate(this);
   mpSizer->Add(pFieldH,0,wxALIGN_LEFT);
   mList.Add(pFieldH);
   
   WXCrystObjBasic* pFieldK=pTex->GetPar(&(pObj->mK)).WXCreate(this);
   mpSizer->Add(pFieldK,0,wxALIGN_LEFT);
   mList.Add(pFieldK);
   
   WXCrystObjBasic* pFieldL=pTex->GetPar(&(pObj->mL)).WXCreate(this);
   mpSizer->Add(pFieldL,0,wxALIGN_LEFT);
   mList.Add(pFieldL);

   this->BottomLayout(0);
   this->CrystUpdate(true);
   VFN_DEBUG_EXIT("WXTexturePhaseMarchDollase::WXTexturePhaseMarchDollase()",5)
}

WXTexturePhaseMarchDollase::~WXTexturePhaseMarchDollase()
{
   mpTexturePhaseMarchDollase->WXNotifyDelete();
}
void WXTexturePhaseMarchDollase::CrystUpdate(const bool uui,const bool lock)
{
   if(lock) mMutex.Lock();
   mList.CrystUpdate(uui,false);
   if(lock) mMutex.Unlock();
}
void WXTexturePhaseMarchDollase::UpdateUI(const bool lock)
{
   if(lock) mMutex.Lock();
   mList.UpdateUI(false);
   if(lock) mMutex.Unlock();
}

////////////////////////////////////////////////////////////////////////
//
//    WXTextureMarchDollase
//
////////////////////////////////////////////////////////////////////////
BEGIN_EVENT_TABLE(WXTextureMarchDollase, wxWindow)
   EVT_MENU(ID_POWDERTEXTURE_MENU_ADDPHASE,   WXTextureMarchDollase::OnAddTexturePhase)
   EVT_MENU(ID_POWDERTEXTURE_MENU_DELETEPHASE,WXTextureMarchDollase::OnDeleteTexturePhase)
   EVT_UPDATE_UI(ID_CRYST_UPDATEUI,           WXRefinableObj::OnUpdateUI)
END_EVENT_TABLE()

WXTextureMarchDollase::WXTextureMarchDollase(wxWindow *parent, TextureMarchDollase*obj):
WXRefinableObj(parent,(RefinableObj*)obj),mpTextureMarchDollase(obj)
{
   VFN_DEBUG_ENTRY("WXTextureMarchDollase::WXTextureMarchDollase()",5)
   // Menu
      mpMenuBar->AddMenu("Phases",ID_REFOBJ_MENU_OBJ);
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_POWDERTEXTURE_MENU_ADDPHASE,
                                "Add Phase");
   //existing phases
      WXRegistry<TexturePhaseMarchDollase> *pWXPhaseRegistry
         =mpTextureMarchDollase->mPhaseRegistry.WXCreate(this);
      mpSizer->Add(pWXPhaseRegistry,0,wxALIGN_LEFT);
      mList.Add(pWXPhaseRegistry);
   this->BottomLayout(0);
   this->CrystUpdate(true);
   this->SetToolTip(_T("Texture for this crystalline phase.\n")
                    _T("You can describe the preferred orientation using ")
                    _T("the March-Dollase model (use the menu).\n\n")
                    _T("Although possible, it is not recommended to enable ")
                    _T("the global optimization of texture parameters, ")
                    _T("as it is *extremely* slow"));
   VFN_DEBUG_EXIT("WXTextureMarchDollase::WXTextureMarchDollase()",5)
}
void WXTextureMarchDollase::OnAddTexturePhase(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_ENTRY("WXTextureMarchDollase::OnAddTexturePhase()",5)
   mpTextureMarchDollase->AddPhase(0.,1.,1,0,0);
   VFN_DEBUG_EXIT("WXTextureMarchDollase::OnAddTexturePhase()",5)
}
void WXTextureMarchDollase::OnDeleteTexturePhase(wxCommandEvent & WXUNUSED(event))
{
}

////////////////////////////////////////////////////////////////////////
//
//    WXPowderPatternDiffraction
//
////////////////////////////////////////////////////////////////////////
static const long ID_POWDERDIFF_PROFILE=                       WXCRYST_ID();
static const long ID_POWDERDIFF_PROFILE_PV=                    WXCRYST_ID();
static const long ID_POWDERDIFF_LEBAIL=                        WXCRYST_ID(); 
static const long ID_POWDERDIFF_PROFILEFITTINGMODE=            WXCRYST_ID(); 

BEGIN_EVENT_TABLE(WXPowderPatternDiffraction, wxWindow)
   EVT_BUTTON(ID_POWDERDIFF_CRYSTAL,WXPowderPatternDiffraction::OnChangeCrystal)
   EVT_MENU(ID_POWDERDIFF_SAVEHKLFCALC, 
                                            WXPowderPatternDiffraction::OnMenuSaveHKLFcalc)
   EVT_MENU(ID_POWDERDIFF_PROFILE_PV,       WXPowderPatternDiffraction::OnChangeProfile)
   EVT_MENU(ID_POWDERDIFF_PROFILE_DEPV,     WXPowderPatternDiffraction::OnChangeProfile)
   EVT_MENU(ID_POWDERDIFF_LEBAIL,           WXPowderPatternDiffraction::OnLeBail)
   EVT_CHECKBOX(ID_POWDERDIFF_PROFILEFITTINGMODE,WXPowderPatternDiffraction::OnLeBail)
   EVT_UPDATE_UI(ID_CRYST_UPDATEUI,         WXRefinableObj::OnUpdateUI)
END_EVENT_TABLE()

WXPowderPatternDiffraction::WXPowderPatternDiffraction(wxWindow *parent,
                                                         PowderPatternDiffraction *p):
WXRefinableObj(parent,p),mpPowderPatternDiffraction(p)
{
   VFN_DEBUG_ENTRY("WXPowderPatternDiffraction::WXPowderPatternDiffraction()",6)
   mpWXTitle->SetForegroundColour(wxColour(0,255,0));
    //Menu
      mpMenuBar->AddMenu("File",ID_REFOBJ_MENU_OBJ);
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_POWDERDIFF_SAVEHKLFCALC,
                                "Save HKL Fcalc");
      mpMenuBar->AddMenu("Profile",ID_POWDERDIFF_PROFILE);
         mpMenuBar->AddMenuItem(ID_POWDERDIFF_PROFILE,ID_POWDERDIFF_PROFILE_PV,
                                "Pseudo-Voigt (X-Ray & monochromatic neutron)");
         mpMenuBar->AddMenuItem(ID_POWDERDIFF_PROFILE,ID_POWDERDIFF_PROFILE_DEPV,
                                "Double-Exponential Pseudo-Voigt (neutron TOF)");
         mpMenuBar->GetMenu(ID_POWDERDIFF_PROFILE).AppendSeparator();
         mpMenuBar->AddMenuItem(ID_POWDERDIFF_PROFILE,ID_POWDERDIFF_LEBAIL,
                                "Profile Fitting + Le Bail Extraction");
      mpSizer->SetItemMinSize(mpMenuBar,
                              mpMenuBar->GetSize().GetWidth(),
                              mpMenuBar->GetSize().GetHeight());
    // Profile Fitting Mode
      mpProfileFittingMode= new wxCheckBox(this,ID_POWDERDIFF_PROFILEFITTINGMODE,
                                           _T("Profile Fitting (Le Bail) Mode"));
      mpSizer->Add(mpProfileFittingMode,0,wxALIGN_LEFT);
    // Crystal Choice
      mpFieldCrystal=new WXFieldChoice(this,ID_POWDERDIFF_CRYSTAL,"Crystal:",300);
      mpSizer->Add(mpFieldCrystal,0,wxALIGN_LEFT);
      mList.Add(mpFieldCrystal);
      mpFieldCrystal->SetToolTip(_T("Crystal structure for this diffraction phase\n")
                                 _T("Click on the button to select another structure"));
   //Global Biso factor
      WXCrystObjBasic* fieldGlobalBiso
         =mpPowderPatternDiffraction->GetPar(&(mpPowderPatternDiffraction->mGlobalBiso))
            .WXCreate(this);
      mList.Add(fieldGlobalBiso);
      mpSizer->Add(fieldGlobalBiso);
   // Texture
      WXTextureMarchDollase* pTex
         =new WXTextureMarchDollase(this,&(mpPowderPatternDiffraction->mCorrTextureMarchDollase));
      mList.Add(pTex);
      mpSizer->Add(pTex);
   // Profile
      
      if(mpPowderPatternDiffraction->mpReflectionProfile!=0)
      {
         VFN_DEBUG_ENTRY("WXPowderPatternDiffraction::WXPowderPatternDiffraction()",6)
         WXCrystObjBasic* pWXProf=mpPowderPatternDiffraction
                                    ->mpReflectionProfile->WXCreate(this);
         mpSizer->Add(pWXProf);
         mList.Add(pWXProf);
         VFN_DEBUG_EXIT("WXPowderPatternDiffraction::WXPowderPatternDiffraction()",6)
      }
      
   this->BottomLayout(0);
   this->CrystUpdate(true);
   VFN_DEBUG_EXIT("WXPowderPatternDiffraction::WXPowderPatternDiffraction()",6)
}

void WXPowderPatternDiffraction::OnChangeCrystal(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXPowderPatternDiffraction::OnChangeCrystal()",6)
   WXCrystValidateAllUserInput();
   int choice;
   Crystal *cryst=dynamic_cast<Crystal*>
      ( WXDialogChooseFromRegistry(gCrystalRegistry,(wxWindow*)this,
         "Choose a Crystal Structure:",choice));
   if(0==cryst) return;
   mpPowderPatternDiffraction->SetCrystal(*cryst);
   this->CrystUpdate(true);
}
void WXPowderPatternDiffraction::OnMenuSaveHKLFcalc(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXPowderPatternDiffraction::OnMenuSaveHKLFcalc()",6)
   WXCrystValidateAllUserInput();
   wxFileDialog save(this,_T("Choose a file to save to"),_T(""),_T(""),_T("*.txt"),wxSAVE | wxOVERWRITE_PROMPT);
   if(save.ShowModal() != wxID_OK) return;
   
   ofstream out(save.GetPath().ToAscii());
   if(!out) return;//:TODO:
   mpPowderPatternDiffraction->PrintFhklCalc(out);
   out.close();
}
void WXPowderPatternDiffraction::UpdateUI(const bool lock)
{
   if(lock) mMutex.Lock();
   mpFieldCrystal->SetValue(mpPowderPatternDiffraction->GetCrystal().GetName());
   mpProfileFittingMode->SetValue(mpPowderPatternDiffraction->GetExtractionMode());
   if(lock) mMutex.Unlock();
   this->WXRefinableObj::UpdateUI(lock);
}
void WXPowderPatternDiffraction::OnChangeProfile(wxCommandEvent & event)
{
   VFN_DEBUG_ENTRY("WXPowderPatternDiffraction::OnChangeProfile()",6)
   bool add=false;
   if(event.GetId()==ID_POWDERDIFF_PROFILE_PV)
   {
      if(mpPowderPatternDiffraction->mpReflectionProfile==0)
      {
         ReflectionProfilePseudoVoigt *p= new ReflectionProfilePseudoVoigt;
         mpPowderPatternDiffraction->SetProfile(p);
         add=true;
      }
      else
         if(mpPowderPatternDiffraction->mpReflectionProfile->GetClassName()
            !="ReflectionProfilePseudoVoigt")
         {
            ReflectionProfilePseudoVoigt *p= new ReflectionProfilePseudoVoigt;
            mpPowderPatternDiffraction->SetProfile(p);
            add=true;
         }
   }
   if(event.GetId()==ID_POWDERDIFF_PROFILE_DEPV)
   {
      if(mpPowderPatternDiffraction->mpReflectionProfile==0)
      {
         ReflectionProfileDoubleExponentialPseudoVoigt *p=
            new ReflectionProfileDoubleExponentialPseudoVoigt
                (mpPowderPatternDiffraction->GetCrystal());
         mpPowderPatternDiffraction->SetProfile(p);
         add=true;
      }
      else
         if(mpPowderPatternDiffraction->mpReflectionProfile->GetClassName()
            !="ReflectionProfileDoubleExponentialPseudoVoigt")
         {
            ReflectionProfileDoubleExponentialPseudoVoigt *p=
               new ReflectionProfileDoubleExponentialPseudoVoigt
                  (mpPowderPatternDiffraction->GetCrystal());
            mpPowderPatternDiffraction->SetProfile(p);
            add=true;
         }
   }
   if(add)
   {
      mpPowderPatternDiffraction->mpReflectionProfile->WXCreate(this);
      mList.Add(mpPowderPatternDiffraction->mpReflectionProfile->WXGet());
      mpSizer->Add(mpPowderPatternDiffraction->mpReflectionProfile->WXGet());
      this->BottomLayout(mpPowderPatternDiffraction->mpReflectionProfile->WXGet());
      mpPowderPatternDiffraction->GetParentPowderPattern().UpdateDisplay();
   }
   VFN_DEBUG_EXIT("WXPowderPatternDiffraction::OnChangeProfile()",6)
}

//////////////////////////////////////// WXProfileFitting /////////////////////
static const long ID_PROFILEFITTING_RUN= WXCRYST_ID();
static const long ID_PROFILEFITTING_RUN_MANUAL= WXCRYST_ID();
static const long ID_PROFILEFITTING_EXPLORE_SPG= WXCRYST_ID();

BEGIN_EVENT_TABLE(WXProfileFitting, wxWindow)
   EVT_BUTTON(ID_PROFILEFITTING_RUN,             WXProfileFitting::OnFit)
   EVT_BUTTON(ID_PROFILEFITTING_RUN_MANUAL,      WXProfileFitting::OnFit)
   EVT_BUTTON(ID_PROFILEFITTING_EXPLORE_SPG,     WXProfileFitting::OnExploreSpacegroups)
END_EVENT_TABLE()

WXProfileFitting::WXProfileFitting(wxWindow *parent,PowderPattern *pPattern,PowderPatternDiffraction *pDiff):
wxWindow(parent,-1),mpPattern(pPattern),mpDiff(pDiff),mLSQ("Profile Fitting object")
{
   wxBoxSizer *pSizer0=new wxBoxSizer(wxVERTICAL);
   this->SetSizer(pSizer0);

   wxNotebook *pNotebook = new wxNotebook(this, -1);
   pSizer0->Add(pNotebook,0,wxALIGN_CENTER);
   // Quick interface
      wxWindow *pQuick=new wxWindow(pNotebook,-1);
      pNotebook->AddPage(pQuick,_T("Quick Fit"),true);
      wxBoxSizer *pSizer=new wxBoxSizer(wxVERTICAL);
      
      wxButton *pButton1=new wxButton(pQuick,ID_PROFILEFITTING_RUN,_T("Le Bail + Fit Profile !"));
      pSizer->Add(pButton1,0,wxALIGN_CENTER);
      
      if(mpDiff==0)
      {
         // List crystal phases
         wxArrayString choices;
         {
            unsigned int nb=mpPattern->GetNbPowderPatternComponent();
            for(unsigned int i=0;i<nb;++i)
               if(mpPattern->GetPowderPatternComponent(i).GetClassName()==string("PowderPatternDiffraction"))
               {
                  pDiff=dynamic_cast<PowderPatternDiffraction*>(&(mpPattern->GetPowderPatternComponent(i)));
                  cout<<"WXProfileFitting::WXProfileFitting():"<<pDiff<<":"<<mpPattern->GetPowderPatternComponent(i).GetName()<<endl;
                  if(pDiff!=0)
                  {
                     choices.Add(wxString::Format(_T("%s, a=%6.3f b=%6.3f c=%6.3f"),
                                                pDiff->GetCrystal().GetName().c_str(),
                                                pDiff->GetCrystal().GetLatticePar(0),
                                                pDiff->GetCrystal().GetLatticePar(1),
                                                pDiff->GetCrystal().GetLatticePar(2)));
                     cout<<"WXProfileFitting::WXProfileFitting():"<<choices[choices.Count()-1]<<","<<pDiff<<endl;
                  }
               }
         }
         if(choices.GetCount()==1)
         {
            cout<<"WXProfileFitting::WXProfileFitting():"<<choices[0]<<","<<pDiff<<endl;
            mpDiff=pDiff;
         }
         else
         {
            wxStaticText *pLabel=new wxStaticText(pQuick,-1,_T("Crystalline Phase to Fit:"));
            pSizer->Add(pLabel,0,wxALIGN_CENTER);
            
            mpList=new wxListBox(pQuick,-1,wxDefaultPosition,wxDefaultSize,choices,wxLB_SINGLE);
            mpList->SetSelection(0);
            pSizer->Add(mpList,0,wxALIGN_CENTER);
         }
      }
      wxArrayString fitChoices;
      if(mpDiff->GetProfile().GetClassName()=="ReflectionProfileDoubleExponentialPseudoVoigt")
      {
         fitChoices.Add(_T("Fit Zero shift"));
         fitChoices.Add(_T("Fit Instrumental Width (alpha1,beta0,beta1)"));
         fitChoices.Add(_T("Fit Size/Strain Broadening (sigma1,gamma2)"));
         fitChoices.Add(_T("Fit Background"));
         fitChoices.Add(_T("Fit Unit Cell"));
      }
      else
      {
         fitChoices.Add(_T("Fit Zero shift"));
         fitChoices.Add(_T("Fit Constant Width"));
         fitChoices.Add(_T("Fit Variable Width"));
         fitChoices.Add(_T("Fit Gaussian-Lorentzian Mixing"));
         fitChoices.Add(_T("Fit Asymmetric parameters"));
         fitChoices.Add(_T("Fit Displacement+Transparency"));
         fitChoices.Add(_T("Fit Background"));
         fitChoices.Add(_T("Fit Unit Cell"));
      }
      mpFitCheckList=new wxCheckListBox(pQuick,-1,wxDefaultPosition,wxDefaultSize,fitChoices);
      if(mpDiff->GetProfile().GetClassName()=="ReflectionProfileDoubleExponentialPseudoVoigt")
      {
         mpFitCheckList->Check(0,false);
         mpFitCheckList->Check(1,false);
         mpFitCheckList->Check(2,true);
         mpFitCheckList->Check(3,true);
         mpFitCheckList->Check(4,true);
      }
      else
      {
         mpFitCheckList->Check(0,true);
         mpFitCheckList->Check(1,true);
         mpFitCheckList->Check(2,true);
         mpFitCheckList->Check(3,true);
         mpFitCheckList->Check(4,true);
         mpFitCheckList->Check(5,true);
         mpFitCheckList->Check(6,true);
         mpFitCheckList->Check(7,true);
      }
      pSizer->Add(mpFitCheckList,0,wxALIGN_CENTER);
      
      pQuick->SetSizer(pSizer);
      pSizer->SetSizeHints(pQuick);
      pQuick->Layout();

   // Manual interface
      // prepare LSQ object
         mLSQ.SetRefinedObj(pDiff->GetParentPowderPattern());
         mLSQ.PrepareRefParList(true);
         //mLSQ.SetParIsUsed(gpRefParTypeObjCryst,false);
         //mLSQ.SetParIsUsed(gpRefParTypeScattDataScale,true);
         //mLSQ.SetParIsUsed(gpRefParTypeScattDataProfile,true);
         //mLSQ.SetParIsUsed(gpRefParTypeScattDataCorrPos,true);
         //mLSQ.SetParIsUsed(gpRefParTypeScattDataBackground,true);
         //mLSQ.SetParIsUsed(gpRefParTypeUnitCell,true);
         mLSQ.SetParIsUsed(gpRefParTypeScatt,false);
         mLSQ.SetParIsUsed(gpRefParTypeScattPow,false);
         //mLSQ.SetParIsFixed(gpRefParTypeObjCryst,true);
         //mLSQ.SetParIsFixed(gpRefParTypeScattDataScale,false);
      // wxLSQ object
         wxScrolledWindow *pManual=new wxScrolledWindow(pNotebook,-1,wxDefaultPosition,
                                                        wxSize(400,250),wxHSCROLL | wxVSCROLL);
         wxBoxSizer *pSizerManual=new wxBoxSizer(wxVERTICAL);

         wxButton *pButton2=new wxButton(pManual,ID_PROFILEFITTING_RUN_MANUAL,_T("Le Bail + Fit Profile !"));
         pSizerManual->Add(pButton2,0,wxALIGN_CENTER);
         //pManual->SetSize(pQuick->GetSize());
         pSizerManual->Add(mLSQ.WXCreate(pManual),0,wxALIGN_CENTER);
         mLSQ.WXGet()->CrystUpdate(true,true);
         pManual->SetScrollRate(10,10);
         pManual->SetSizer(pSizerManual);
         pManual->Layout();
         pSizerManual->SetVirtualSizeHints(pManual);
      
      pNotebook->AddPage(pManual,_T("Manual Fit"),true);
   
   // Spacegroup exploration
      wxScrolledWindow *pSpgExplor=new wxScrolledWindow(pNotebook,-1,wxDefaultPosition,
                                                         wxSize(400,250),wxHSCROLL | wxVSCROLL);
      wxBoxSizer *pSizerSpgExplor=new wxBoxSizer(wxVERTICAL);

      wxButton *pButton3=new wxButton(pSpgExplor,ID_PROFILEFITTING_EXPLORE_SPG,_T("Try all possible spacegroups !"));
      
      pNotebook->AddPage(pSpgExplor,_T("Spacegroup Explorer"),true);
   
   pNotebook->ChangeSelection(0);
   mpLog =new wxTextCtrl(this,-1,_T(""),wxDefaultPosition,wxSize(400,250),wxTE_MULTILINE|wxTE_READONLY|wxTE_DONTWRAP);
   mpLog->SetFont(wxFont(9,wxTELETYPE,wxFONTSTYLE_NORMAL,wxFONTWEIGHT_NORMAL));
   pSizer0->Add(mpLog,0,wxALIGN_CENTER);


   //pSizer0->Layout();
   pSizer0->SetSizeHints(this);
   //pSizer0->Fit(this);
   pSizer0->Fit(this->GetParent());
   this->Layout();
}

WXProfileFitting::~WXProfileFitting()
{
   if(mpDiff!=0) mpDiff->SetExtractionMode(false);
   else
   {
      unsigned int nb=mpPattern->GetNbPowderPatternComponent();
      for(unsigned int i=0;i<nb;++i)
         if(mpPattern->GetPowderPatternComponent(i).GetClassName()==string("PowderPatternDiffraction"))
         {
            PowderPatternDiffraction *pDiff=dynamic_cast<PowderPatternDiffraction*>(&(mpPattern->GetPowderPatternComponent(i)));
            if(pDiff!=0) pDiff->SetExtractionMode(false);
         }
   }
   mpPattern->UpdateDisplay();
}

void WXProfileFitting::OnFit(wxCommandEvent &event)
{
   // First, Le Bail
   PowderPatternDiffraction *pDiff=0;
   if(mpDiff!=0) pDiff=mpDiff;
   else
   {
      unsigned int nb=mpPattern->GetNbPowderPatternComponent();
      unsigned int n=0;
      const unsigned int n0=mpList->GetSelection();
      for(unsigned int i=0;i<nb;++i)
      {
         if(mpPattern->GetPowderPatternComponent(i).GetClassName()==string("PowderPatternDiffraction"))
         {
            pDiff=dynamic_cast<PowderPatternDiffraction*>(&(mpPattern->GetPowderPatternComponent(i)));
            if(pDiff!=0)
            {
               if(n++==n0) break;
            }
         }
      }
   }
   cout<<"Selected PowderPatternDiffraction:"<<pDiff->GetName()<<","<<pDiff->GetCrystal().GetName()<<endl;

   pDiff->SetExtractionMode(true,true);
   
   mpLog->AppendText(wxString::Format(_T("Starting 20 Le Bail cycles\n")));
   wxProgressDialog dlgProgress(_T("Le Bail and Profile Fitting"),_T("Le Bail Fitting, cycle #0/20"),
                                 18,this,wxPD_AUTO_HIDE|wxPD_ELAPSED_TIME|wxPD_CAN_ABORT);
   for(int i=0;i<10;++i)
   {
      pDiff->ExtractLeBail(2);
      pDiff->GetParentPowderPattern().FitScaleFactorForRw();
      pDiff->GetParentPowderPattern().UpdateDisplay();
      if(dlgProgress.Update(i,wxString::Format(_T("Le Bail Fitting, cycle #%d/20"),i*2))==false) return;
   }
   mpLog->AppendText(wxString::Format(_T("                  => Rwp=%5.3f%%, GoF=%7.3f\n"),
                                    pDiff->GetParentPowderPattern().GetRw()*100,
                                    pDiff->GetParentPowderPattern().GetChi2()
                                    /pDiff->GetParentPowderPattern().GetNbPointUsed()));

   if(event.GetId()==ID_PROFILEFITTING_RUN)
   {
      bool fitzero=false,fitwidth0=false,fitwidth=false,fiteta=false,
         fitasym=false,fitdispltransp=false,fitbackgd=false,fitcell=false,
         fitTOFInstWidth=false,fitTOFBroadening=false;
      if(mpDiff->GetProfile().GetClassName()=="ReflectionProfileDoubleExponentialPseudoVoigt")
      {
         mpDiff->GetProfile().Print();
         fitzero=mpFitCheckList->IsChecked(0);
         fitTOFInstWidth=mpFitCheckList->IsChecked(1);
         fitTOFBroadening=mpFitCheckList->IsChecked(2);
         fitbackgd=mpFitCheckList->IsChecked(3);
         fitcell=mpFitCheckList->IsChecked(4);
      }
      else
      {
         fitzero=mpFitCheckList->IsChecked(0);
         fitwidth0=mpFitCheckList->IsChecked(1);
         fitwidth=mpFitCheckList->IsChecked(2);
         fiteta=mpFitCheckList->IsChecked(3);
         fitasym=mpFitCheckList->IsChecked(4);
         fitdispltransp=mpFitCheckList->IsChecked(5);
         fitbackgd=mpFitCheckList->IsChecked(6);
         fitcell=mpFitCheckList->IsChecked(7);
      }
      try{
         
         // :TODO: take car of other profiles than pseudo-voigt (DE-PV)
         if(fitzero) mLSQ.SetParIsFixed("Zero",false);
         if(fitwidth0) mLSQ.SetParIsFixed("W",false);
         if(fitzero||fitwidth0)
         {
            mpLog->AppendText(wxString::Format(_T("Fitting zero shift && constant width\n")));
            if(dlgProgress.Update(11,_T("Fitting zero shift && constant width"))==false) return;
            mLSQ.Refine(5,true,false);
            pDiff->ExtractLeBail(2);
            //pDiff->GetParentPowderPattern().FitScaleFactorForRw();
            pDiff->GetParentPowderPattern().UpdateDisplay();
            mpLog->AppendText(wxString::Format(_T("                  => Rwp=%6.3f%%, GoF=%7.3f\n"),
                                             pDiff->GetParentPowderPattern().GetRw()*100,
                                             pDiff->GetParentPowderPattern().GetChi2()
                                             /pDiff->GetParentPowderPattern().GetNbPointUsed()));
         }
         if(fitwidth) mLSQ.SetParIsFixed("U",false);
         if(fitwidth) mLSQ.SetParIsFixed("V",false);
         if(fiteta) mLSQ.SetParIsFixed("Eta0",false);
         if(fitwidth||fiteta)
         {
            mpLog->AppendText(wxString::Format(_T("Fitting width and gaussian/lorentzian fixed mix\n")));
            if(dlgProgress.Update(12,_T("Fitting variable width and gaussian/lorentzian fixed mix"))==false) return;
            mLSQ.Refine(5,true,false);
            pDiff->ExtractLeBail(2);
            //pDiff->GetParentPowderPattern().FitScaleFactorForRw();
            pDiff->GetParentPowderPattern().UpdateDisplay();
            mpLog->AppendText(wxString::Format(_T("                  => Rwp=%6.3f%%, GoF=%7.3f\n"),
                                             pDiff->GetParentPowderPattern().GetRw()*100,
                                             pDiff->GetParentPowderPattern().GetChi2()
                                             /pDiff->GetParentPowderPattern().GetNbPointUsed()));
         }
         
         if(fitTOFInstWidth)
         {// TOF
            mpDiff->GetProfile().Print();
            mLSQ.SetParIsFixed("Alpha1",false);
            mLSQ.SetParIsFixed("Beta0",false);
            mLSQ.SetParIsFixed("Beta1",false);
            mpLog->AppendText(wxString::Format(_T("Fitting TOF instrumental width (alpha1,beta0,beta1)\n")));
            if(dlgProgress.Update(12,_T("Fitting TOF instrumental width (alpha1,beta0,beta1)"))==false) return;
            mLSQ.Refine(5,true,false);
            pDiff->ExtractLeBail(2);
            //pDiff->GetParentPowderPattern().FitScaleFactorForRw();
            pDiff->GetParentPowderPattern().UpdateDisplay();
            mpLog->AppendText(wxString::Format(_T("                  => Rwp=%6.3f%%, GoF=%7.3f\n"),
                                             pDiff->GetParentPowderPattern().GetRw()*100,
                                             pDiff->GetParentPowderPattern().GetChi2()
                                             /pDiff->GetParentPowderPattern().GetNbPointUsed()));
         }
         if(fitTOFBroadening)
         {// TOF
            mpDiff->GetProfile().Print();
            mLSQ.SetParIsFixed("GaussianSigma1",false);
            //lsqobj.SetParIsFixed("LorentzianGamma2",false);
            //lsqobj.SetLimitsAbsolute("GaussianSigma1",0,1e6);
            //lsqobj.SetLimitsAbsolute("LorentzianGamma2",0,1e6);
            mpLog->AppendText(wxString::Format(_T("Fitting size/strain broadening parameters (sigma1,gamma2)\n")));
            if(dlgProgress.Update(12,_T("Fitting size/strain broadening parameters (sigma1,gamma2)"))==false) return;
            mLSQ.Refine(5,true,false);
            pDiff->ExtractLeBail(2);
            //pDiff->GetParentPowderPattern().FitScaleFactorForRw();
            pDiff->GetParentPowderPattern().UpdateDisplay();
            mpLog->AppendText(wxString::Format(_T("                  => Rwp=%6.3f%%, GoF=%7.3f\n"),
                                             pDiff->GetParentPowderPattern().GetRw()*100,
                                             pDiff->GetParentPowderPattern().GetChi2()
                                             /pDiff->GetParentPowderPattern().GetNbPointUsed()));
            mLSQ.SetParIsFixed("GaussianSigma1",true);
         }
         
         if(fiteta) mLSQ.SetParIsFixed("Eta1",false);
         if(fiteta)
         {
            mpLog->AppendText(wxString::Format(_T("Fitting gaussian/lorentzian mix\n")));
            if(dlgProgress.Update(13,_T("Fitting variable width and gaussian/lorentzian mix"))==false) return;
            mLSQ.Refine(5,true,false);
            pDiff->ExtractLeBail(2);
            //pDiff->GetParentPowderPattern().FitScaleFactorForRw();
            pDiff->GetParentPowderPattern().UpdateDisplay();
            mpLog->AppendText(wxString::Format(_T("                  => Rwp=%6.3f%%, GoF=%7.3f\n"),
                                             pDiff->GetParentPowderPattern().GetRw()*100,
                                             pDiff->GetParentPowderPattern().GetChi2()
                                             /pDiff->GetParentPowderPattern().GetNbPointUsed()));
         }
         
         if(fitasym) mLSQ.SetParIsFixed("Asym0",false);
         if(fitasym) mLSQ.SetParIsFixed("Asym1",false);
         if(fitasym) mLSQ.SetParIsFixed("Asym2",false);
         if(fitdispltransp) mLSQ.SetParIsFixed("2ThetaDispl",false);
         if(fitdispltransp) mLSQ.SetParIsFixed("2ThetaTransp",false);
         if(fitdispltransp||fitasym)
         {
            mpLog->AppendText(wxString::Format(_T("Fitting assymetry and sample displacement/transparency\n")));
            if(dlgProgress.Update(14,_T("Fitting assymetry and sample displacement/transparency"))==false) return;
            mLSQ.Refine(5,true,false);
            pDiff->ExtractLeBail(2);
            //pDiff->GetParentPowderPattern().FitScaleFactorForRw();
            pDiff->GetParentPowderPattern().UpdateDisplay();
            mpLog->AppendText(wxString::Format(_T("                  => Rwp=%6.3f%%, GoF=%7.3f\n"),
                                             pDiff->GetParentPowderPattern().GetRw()*100,
                                             pDiff->GetParentPowderPattern().GetChi2()
                                             /pDiff->GetParentPowderPattern().GetNbPointUsed()));
         }
         
         if(fitbackgd)
         {
            mLSQ.SetParIsFixed(gpRefParTypeScattDataBackground,false);
            // Make sure points beyond max resolution are not optimized
            const unsigned int nbcomp= pDiff->GetParentPowderPattern().GetNbPowderPatternComponent();
            for(unsigned int i=0;i<nbcomp;++i)
               if(pDiff->GetParentPowderPattern().GetPowderPatternComponent(i).GetClassName()=="PowderPatternBackground")
               {
                  PowderPatternBackground *pback=dynamic_cast<PowderPatternBackground *>
                     (&(pDiff->GetParentPowderPattern().GetPowderPatternComponent(i)));
                  pback->FixParametersBeyondMaxresolution(mLSQ.GetCompiledRefinedObj());
               }
   
            mpLog->AppendText(wxString::Format(_T("Fitting background\n")));
            if(dlgProgress.Update(15,_T("Fitting background"))==false) return;
            mLSQ.Refine(5,true,false);
            pDiff->ExtractLeBail(2);
            //pDiff->GetParentPowderPattern().FitScaleFactorForRw();
            pDiff->GetParentPowderPattern().UpdateDisplay();
            mpLog->AppendText(wxString::Format(_T("                  => Rwp=%6.3f%%, GoF=%7.3f\n"),
                                             pDiff->GetParentPowderPattern().GetRw()*100,
                                             pDiff->GetParentPowderPattern().GetChi2()
                                             /pDiff->GetParentPowderPattern().GetNbPointUsed()));
         }
         
         if(fitcell) mLSQ.SetParIsFixed(gpRefParTypeUnitCell,false);
         if(fitcell)
         {
            mpLog->AppendText(wxString::Format(_T("Fitting unit cell\n")));
            if(dlgProgress.Update(16,_T("Fitting unit cell"))==false) return;
            mLSQ.Refine(5,true,false);
            pDiff->ExtractLeBail(2);
            //pDiff->GetParentPowderPattern().FitScaleFactorForRw();
            pDiff->GetParentPowderPattern().UpdateDisplay();
            mpLog->AppendText(wxString::Format(_T("                  => Rwp=%6.3f%%, GoF=%7.3f\n"),
                                             pDiff->GetParentPowderPattern().GetRw()*100,
                                             pDiff->GetParentPowderPattern().GetChi2()
                                             /pDiff->GetParentPowderPattern().GetNbPointUsed()));
         }
      }
      catch(const ObjCrystException &except)
      {
         mpLog->AppendText(wxString::Format(_T(" OOPS : refinement diverged ! Aborting.")));
      }
   }
   else
   {// Manual fit
      try
      {
         mpLog->AppendText(wxString::Format(_T("Profile fitting (manual):\n")));
         mpLog->AppendText(wxString::Format(_T("    Initial values:  Rwp=%6.3f%%, GoF=%7.3f\n"),
                                                pDiff->GetParentPowderPattern().GetRw()*100,
                                                pDiff->GetParentPowderPattern().GetChi2()
                                                /pDiff->GetParentPowderPattern().GetNbPointUsed()));
         mpLog->AppendText(wxString::Format(_T("3 LSQ cycles...\n")));
         mLSQ.Refine(3,true,false);
         mpLog->AppendText(wxString::Format(_T("                  => Rwp=%6.3f%%, GoF=%7.3f\n"),
                                                pDiff->GetParentPowderPattern().GetRw()*100,
                                                pDiff->GetParentPowderPattern().GetChi2()
                                                /pDiff->GetParentPowderPattern().GetNbPointUsed()));
         mpLog->AppendText(wxString::Format(_T("2 Le Bail cycles...\n")));
         if(dlgProgress.Update(13,_T("Manual Le Bail + Profile fitting"))==false) return;
         pDiff->ExtractLeBail(2);
         mpLog->AppendText(wxString::Format(_T("                  => Rwp=%6.3f%%, GoF=%7.3f\n"),
                                                pDiff->GetParentPowderPattern().GetRw()*100,
                                                pDiff->GetParentPowderPattern().GetChi2()
                                                /pDiff->GetParentPowderPattern().GetNbPointUsed()));
         mpLog->AppendText(wxString::Format(_T("3 LSQ cycles...\n")));
         if(dlgProgress.Update(16,_T("Manual Le Bail + Profile fitting"))==false) return;
         mLSQ.Refine(3,true,false);
         if(dlgProgress.Update(19,_T("Manual Le Bail + Profile fitting"))==false) return;
         pDiff->GetParentPowderPattern().UpdateDisplay();
         mpLog->AppendText(wxString::Format(_T("                  => Rwp=%6.3f%%, GoF=%7.3f\n"),
                                                pDiff->GetParentPowderPattern().GetRw()*100,
                                                pDiff->GetParentPowderPattern().GetChi2()
                                                /pDiff->GetParentPowderPattern().GetNbPointUsed()));
      }
      catch(const ObjCrystException &except)
      {
         mpLog->AppendText(wxString::Format(_T(" OOPS : refinement diverged ! Aborting.")));
      }
   }
   mLSQ.WXGet()->CrystUpdate(true,true);
   pDiff->GetCrystal().UpdateDisplay();
}

struct SPGScore
{
   SPGScore(const string &s, const REAL r, const REAL g):
   hm(s),rw(r),gof(g)
   {}
   string hm;
   REAL rw,gof;
};

bool compareSPGScore(const SPGScore &s1, const SPGScore &s2)
{
   return s1.gof < s2.gof;
}


void WXProfileFitting::OnExploreSpacegroups(wxCommandEvent &event)
{
   PowderPatternDiffraction *pDiff=0;
   if(mpDiff!=0) pDiff=mpDiff;
   else
   {
      unsigned int nb=mpPattern->GetNbPowderPatternComponent();
      unsigned int n=0;
      const unsigned int n0=mpList->GetSelection();
      for(unsigned int i=0;i<nb;++i)
      {
         if(mpPattern->GetPowderPatternComponent(i).GetClassName()==string("PowderPatternDiffraction"))
         {
            pDiff=dynamic_cast<PowderPatternDiffraction*>(&(mpPattern->GetPowderPatternComponent(i)));
            if(pDiff!=0)
            {
               if(n++==n0) break;
            }
         }
      }
   }
   cout<<"Selected PowderPatternDiffraction:"<<pDiff->GetName()<<","<<pDiff->GetCrystal().GetName()<<endl;

   pDiff->SetExtractionMode(true,true);
   Crystal *pCrystal=&(pDiff->GetCrystal());
   
   // Keep initial lattice parameters & spg
   const REAL a=pCrystal->GetLatticePar(0),
              b=pCrystal->GetLatticePar(1),
              c=pCrystal->GetLatticePar(2),
              d=pCrystal->GetLatticePar(3),
              e=pCrystal->GetLatticePar(4),
              f=pCrystal->GetLatticePar(5);
   const string spgname=pCrystal->GetSpaceGroup().GetName();
   const string name=pCrystal->GetName();
   
   const cctbx::uctbx::unit_cell uc(scitbx::af::double6(a,b,c,d*RAD2DEG,e*RAD2DEG,f*RAD2DEG));
   
   cctbx::sgtbx::space_group_symbol_iterator it=cctbx::sgtbx::space_group_symbol_iterator();
   
   // First, count compatible spacegroups
   unsigned int nbspg=0;
   //unsigned int hmlen=0;
   for(;;)
   {
      cctbx::sgtbx::space_group_symbols s=it.next();
      if(s.number()==0) break;
      cctbx::sgtbx::space_group spg(s);
      if(spg.is_compatible_unit_cell(uc,0.01,0.1)) nbspg++;
      //if(s.universal_hermann_mauguin().size()>hmlen) hmlen=s.universal_hermann_mauguin().size();
   }
   mpLog->AppendText(wxString::Format(_T("Beginning spacegroup exploration... %d to go...\n"),nbspg));
   //cout<<"Max HM symbol length:"<<hmlen<<endl;
   const unsigned int nbcycle=3;
   
   wxProgressDialog dlgProgress(_T("Trying compatible spacegroups"),_T("Starting........\n......\n......"),
                                 nbspg*nbcycle,this,wxPD_AUTO_HIDE|wxPD_ELAPSED_TIME|wxPD_CAN_ABORT);
   
   list<SPGScore> vSPG;
   // Try & optimize every spacegroup
   it=cctbx::sgtbx::space_group_symbol_iterator();
   for(int i=0;;)
   {
      cctbx::sgtbx::space_group_symbols s=it.next();
      if(s.number()==0) break;
      cctbx::sgtbx::space_group spg(s);
      bool compat=spg.is_compatible_unit_cell(uc,0.01,0.1);
      if(compat)
      {
         i++;
         const string hm=s.universal_hermann_mauguin();
         //cout<<s.number()<<","<<hm.c_str()<<","<<(int)compat<<endl;
         mpLog->AppendText(wxString::Format(_T(" (#%3d) %-14s:"),s.number(),hm.c_str()));
         pCrystal->Init(a,b,c,d,e,f,hm,name);
         pDiff->GetParentPowderPattern().UpdateDisplay();
         for(unsigned int j=0;j<nbcycle;j++)
         {
            // First, Le Bail
            pDiff->SetExtractionMode(true,true);
            pDiff->ExtractLeBail(5);
            pDiff->GetParentPowderPattern().FitScaleFactorForRw();
            //mpLog->AppendText(wxString::Format(_T("/%5.2f"),pDiff->GetParentPowderPattern().GetRw()*100));
            // LSQ refinement
            LSQNumObj lsq;
            lsq.SetRefinedObj(pDiff->GetParentPowderPattern());
            lsq.PrepareRefParList(true);
            lsq.SetParIsFixed(gpRefParTypeObjCryst,true);
            lsq.SetParIsFixed(gpRefParTypeScattDataScale,false);
            
            // Only do the full monty for P1, keep the parameters for other spacegroups
            if(s.number()==1) lsq.SetParIsFixed("Zero",false);
            lsq.SetParIsFixed(gpRefParTypeUnitCell,false);
            lsq.Refine(2,true,false);
            if(s.number()==1) 
            {
               lsq.SetParIsFixed("2ThetaDispl",false);
               lsq.SetParIsFixed("2ThetaTransp",false);
               lsq.Refine(2,true,false);
               lsq.SetParIsFixed(gpRefParTypeScattDataBackground,false);
               // Fix background point beyond optimized domain
               const unsigned int nbcomp= pDiff->GetParentPowderPattern().GetNbPowderPatternComponent();
               for(unsigned int i=0;i<nbcomp;++i)
                  if(pDiff->GetParentPowderPattern().GetPowderPatternComponent(i).GetClassName()=="PowderPatternBackground")
                  {
                     PowderPatternBackground *pback=dynamic_cast<PowderPatternBackground *>
                        (&(pDiff->GetParentPowderPattern().GetPowderPatternComponent(i)));
                     pback->FixParametersBeyondMaxresolution(lsq.GetCompiledRefinedObj());
                  }
               
               lsq.Refine(2,true,false);
            }
            // restart from equal intensities
            pDiff->SetExtractionMode(true,true);
            pDiff->ExtractLeBail(5);
            lsq.Refine(3,true,false);
            //mpLog->AppendText(wxString::Format(_T("%5.2f%%/"),pDiff->GetParentPowderPattern().GetRw()*100));
            pDiff->GetParentPowderPattern().FitScaleFactorForRw();
            pDiff->GetParentPowderPattern().UpdateDisplay();
            const REAL rw=pDiff->GetParentPowderPattern().GetRw()*100;
            const REAL gof=pDiff->GetParentPowderPattern().GetChi2()
                           /(pDiff->GetParentPowderPattern().GetNbPointUsed()-pDiff->GetNbReflBelowMaxSinThetaOvLambda());
            if(dlgProgress.Update(i*nbcycle+j,wxString::Format(_T("%s  (cycle #%u)\n   Rwp=%5.2f%%\n   GoF=%6.3f"),
                                                               hm.c_str(),j,rw,gof))==false) return;
         }
         const REAL rw=pDiff->GetParentPowderPattern().GetRw()*100;
         const REAL gof=pDiff->GetParentPowderPattern().GetChi2()
                        /(pDiff->GetParentPowderPattern().GetNbPointUsed()-pDiff->GetNbReflBelowMaxSinThetaOvLambda());
         vSPG.push_back(SPGScore(hm.c_str(),rw,gof));
         mpLog->AppendText(wxString::Format(_T(" Rwp= %5.2f%%  GoF=%6.3f\n"),rw,gof));
      }
   }
   // sort results by GoF
   vSPG.sort(compareSPGScore);
   mpLog->AppendText(wxString::Format(_T("\n\n BEST Solutions, from min_GoF to 2*min_Gof:\n")));
   for(list<SPGScore>::const_iterator pos=vSPG.begin();pos!=vSPG.end();++pos)
   {
      if(pos->gof>(2*vSPG.begin()->gof)) break;
      mpLog->AppendText(wxString::Format(_T(" %-14s: Rwp= %5.2f%%  GoF=%6.3f\n"),pos->hm.c_str(),pos->rw,pos->gof));
   }
   // Back to original spacegroup
   pCrystal->Init(a,b,c,d,e,f,spgname,name);
   pDiff->GetParentPowderPattern().UpdateDisplay();
   pDiff->SetExtractionMode(true,true);
   pDiff->ExtractLeBail(5);
   pDiff->GetParentPowderPattern().UpdateDisplay();
}

void WXPowderPatternDiffraction::OnLeBail(wxCommandEvent &event)
{
   if((event.GetId()==ID_POWDERDIFF_PROFILEFITTINGMODE)&&(mpProfileFittingMode->GetValue()==false))
   {
      mpPowderPatternDiffraction->SetExtractionMode(false);
      mpPowderPatternDiffraction->GetParentPowderPattern().UpdateDisplay();
      return;
   }
   mpPowderPatternDiffraction->SetExtractionMode(true,true);
   wxFrame *pFrame=new wxFrame(this,-1,_T("Profile Fitting"));
   WXProfileFitting *pFit;
   pFit=new WXProfileFitting(pFrame,&(mpPowderPatternDiffraction->GetParentPowderPattern()),mpPowderPatternDiffraction);
   pFrame->Show(true);
}


////////////////////////////////////////////////////////////////////////
//
//    WXProfilePseudoVoigt
//
////////////////////////////////////////////////////////////////////////
WXProfilePseudoVoigt::WXProfilePseudoVoigt(wxWindow *parent, ReflectionProfilePseudoVoigt *prof):
WXCrystObj(parent),mpProfile(prof)
{
   VFN_DEBUG_ENTRY("WXProfilePseudoVoigt::WXProfilePseudoVoigt()",6)
   mpWXTitle->SetLabel("Pseudo-Voigt profile");
   mpWXTitle->SetForegroundColour(wxColour(0,0,255));
   // Width
      wxBoxSizer* sizer1=new wxBoxSizer(wxHORIZONTAL);
      WXCrystObjBasic* pFieldCagliotiU=mpProfile->GetPar("U").WXCreate(this);
      WXCrystObjBasic* pFieldCagliotiV=mpProfile->GetPar("V").WXCreate(this);
      WXCrystObjBasic* pFieldCagliotiW=mpProfile->GetPar("W").WXCreate(this);;
      sizer1->Add(pFieldCagliotiU,0);
      sizer1->Add(pFieldCagliotiV,0);
      sizer1->Add(pFieldCagliotiW,0);
      mList.Add(pFieldCagliotiU);
      mList.Add(pFieldCagliotiV);
      mList.Add(pFieldCagliotiW);
      mpSizer->Add(sizer1);
      pFieldCagliotiU->SetToolTip(_T("Width Parameters (Caglioti's law):\n")
                                  _T("fwhm=[W+V*tan(theta)+U*tan^2(theta)]^1/2"));
      pFieldCagliotiV->SetToolTip(_T("Width Parameters (Caglioti's law):\n")
                                  _T("fwhm=[W+V*tan(theta)+U*tan^2(theta)]^1/2"));
      pFieldCagliotiW->SetToolTip(_T("Width Parameters (Caglioti's law):\n")
                                  _T("fwhm=[W+V*tan(theta)+U*tan^2(theta)]^1/2"));
   // Mixing parameter
      wxBoxSizer* sizer2=new wxBoxSizer(wxHORIZONTAL);
      WXCrystObjBasic* pFieldEta0=mpProfile->GetPar("Eta0").WXCreate(this);
      WXCrystObjBasic* pFieldEta1=mpProfile->GetPar("Eta1").WXCreate(this);
      sizer2->Add(pFieldEta0,0);
      sizer2->Add(pFieldEta1,0);
      mList.Add(pFieldEta0);
      mList.Add(pFieldEta1);
      mpSizer->Add(sizer2);
      pFieldEta0->SetToolTip(_T("Gaussian/Lorentzian mixing parameters:\n")
                             _T(" PV(x) = eta*L(x) + (1-eta)*G(x)\n\n")
                             _T("eta=Eta0+Eta11*2theta"));
      pFieldEta1->SetToolTip(_T("Gaussian/Lorentzian mixing parameters:\n")
                             _T(" PV(x) = eta*L(x) + (1-eta)*G(x)\n\n")
                             _T("eta=Eta0+Eta11*2theta"));
   // Asymmetry parameter
      wxBoxSizer* sizer3=new wxBoxSizer(wxHORIZONTAL);
      //WXCrystObjBasic* pFieldAsymA0=mpProfile->GetPar("AsymA0").WXCreate(this);
      //WXCrystObjBasic* pFieldAsymA1=mpProfile->GetPar("AsymA1").WXCreate(this);
      //WXCrystObjBasic* pFieldAsymB0=mpProfile->GetPar("AsymB0").WXCreate(this);
      //WXCrystObjBasic* pFieldAsymB1=mpProfile->GetPar("AsymB1").WXCreate(this);
      //sizer3->Add(pFieldAsymA0,0);
      //sizer3->Add(pFieldAsymA1,0);
      //sizer3->Add(pFieldAsymB0,0);
      //sizer3->Add(pFieldAsymB1,0);
      //mList.Add(pFieldAsymA0);
      //mList.Add(pFieldAsymA1);
      //mList.Add(pFieldAsymB0);
      //mList.Add(pFieldAsymB1);
      WXCrystObjBasic* pFieldAsym0=mpProfile->GetPar("Asym0").WXCreate(this);
      WXCrystObjBasic* pFieldAsym1=mpProfile->GetPar("Asym1").WXCreate(this);
      WXCrystObjBasic* pFieldAsym2=mpProfile->GetPar("Asym2").WXCreate(this);
      sizer3->Add(pFieldAsym0,0);
      sizer3->Add(pFieldAsym1,0);
      sizer3->Add(pFieldAsym2,0);
      mList.Add(pFieldAsym0);
      mList.Add(pFieldAsym1);
      mList.Add(pFieldAsym2);
      mpSizer->Add(sizer3);

      pFieldAsym0->SetToolTip(_T("Asymmetry parameters:\n\n")
                              _T("A=A0+A1/sin(2theta)+A2/sin^2(2theta) "));
      pFieldAsym1->SetToolTip(_T("Asymmetry parameters:\n\n")
                              _T("A=A0+A1/sin(2theta)+A2/sin^2(2theta) "));
      pFieldAsym2->SetToolTip(_T("Asymmetry parameters:\n\n")
                              _T("A=A0+A1/sin(2theta)+A2/sin^2(2theta) "));
   
   this->BottomLayout(0);
   this->CrystUpdate(true);
   VFN_DEBUG_EXIT("WXProfilePseudoVoigt::WXProfilePseudoVoigt()",6)
}
WXProfilePseudoVoigt::~WXProfilePseudoVoigt()
{
   mpProfile->WXNotifyDelete();
}
bool WXProfilePseudoVoigt::OnChangeName(const int id)
{
   return false;
}
////////////////////////////////////////////////////////////////////////
//
//    WXProfileDoubleExponentialPseudoVoigt
//
////////////////////////////////////////////////////////////////////////
WXProfileDoubleExponentialPseudoVoigt::WXProfileDoubleExponentialPseudoVoigt
   (wxWindow *parent, ReflectionProfileDoubleExponentialPseudoVoigt *prof):
WXCrystObj(parent),mpProfile(prof)
{
   VFN_DEBUG_ENTRY("WXProfileDoubleExponentialPseudoVoigt::WXProfileDoubleExponentialPseudoVoigt()",6)
   mpWXTitle->SetLabel("Double-Exponential Pseudo-Voigt profile (for neutron TOF)");
   mpWXTitle->SetForegroundColour(wxColour(0,0,255));
   mpWXTitle->BottomLayout(0);
   // Instrumental
      wxBoxSizer* sizer1=new wxBoxSizer(wxHORIZONTAL);
      WXCrystObjBasic* pFieldCagliotiA0=mpProfile->GetPar("Alpha0").WXCreate(this);
      WXCrystObjBasic* pFieldCagliotiA =mpProfile->GetPar("Alpha1").WXCreate(this);
      WXCrystObjBasic* pFieldCagliotiB0=mpProfile->GetPar("Beta0").WXCreate(this);
      WXCrystObjBasic* pFieldCagliotiB1=mpProfile->GetPar("Beta1").WXCreate(this);
      sizer1->Add(pFieldCagliotiA0,0);
      sizer1->Add(pFieldCagliotiA,0);
      sizer1->Add(pFieldCagliotiB0,0);
      sizer1->Add(pFieldCagliotiB1,0);
      mList.Add(pFieldCagliotiA0);
      mList.Add(pFieldCagliotiA);
      mList.Add(pFieldCagliotiB0);
      mList.Add(pFieldCagliotiB1);
      mpSizer->Add(sizer1);
   // Instrumental
      wxBoxSizer* sizer2=new wxBoxSizer(wxHORIZONTAL);
      WXCrystObjBasic* pFieldSigma0=mpProfile->GetPar("GaussianSigma0").WXCreate(this);
      WXCrystObjBasic* pFieldSigma1=mpProfile->GetPar("GaussianSigma1").WXCreate(this);
      WXCrystObjBasic* pFieldSigma2=mpProfile->GetPar("GaussianSigma2").WXCreate(this);
      sizer2->Add(pFieldSigma0,0);
      sizer2->Add(pFieldSigma1,0);
      sizer2->Add(pFieldSigma2,0);
      mList.Add(pFieldSigma0);
      mList.Add(pFieldSigma1);
      mList.Add(pFieldSigma2);
      mpSizer->Add(sizer2);
   // Instrumental
      wxBoxSizer* sizer3=new wxBoxSizer(wxHORIZONTAL);
      WXCrystObjBasic* pFieldGamma0=mpProfile->GetPar("LorentzianGamma0").WXCreate(this);
      WXCrystObjBasic* pFieldGamma1=mpProfile->GetPar("LorentzianGamma1").WXCreate(this);
      WXCrystObjBasic* pFieldGamma2=mpProfile->GetPar("LorentzianGamma2").WXCreate(this);
      sizer3->Add(pFieldGamma0,0);
      sizer3->Add(pFieldGamma1,0);
      sizer3->Add(pFieldGamma2,0);
      mList.Add(pFieldGamma0);
      mList.Add(pFieldGamma1);
      mList.Add(pFieldGamma2);
      mpSizer->Add(sizer3);
   
   this->BottomLayout(0);
   this->CrystUpdate(true);
   VFN_DEBUG_EXIT("WXProfileDoubleExponentialPseudoVoigt::WXProfileDoubleExponentialPseudoVoigt()",6)
}
WXProfileDoubleExponentialPseudoVoigt::~WXProfileDoubleExponentialPseudoVoigt()
{
   mpProfile->WXNotifyDelete();
}
bool WXProfileDoubleExponentialPseudoVoigt::OnChangeName(const int id)
{
   return false;
}

}// namespace 

