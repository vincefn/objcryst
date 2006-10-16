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

#include "wxCryst/wxPowderPattern.h"
#include "wxCryst/wxRadiation.h"
#include "RefinableObj/Simplex.h"
#include "ObjCryst/PowderPatternBackgroundBayesianMinimiser.h"
#include "Quirks/VFNStreamFormat.h"

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
}

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
static const long ID_POWDERBACKGROUND_IMPORT=               WXCRYST_ID(); 
static const long ID_POWDERBACKGROUND_OPTIMIZEBAYESIAN=     WXCRYST_ID(); 
static const long ID_POWDERDIFF_CRYSTAL=                    WXCRYST_ID(); 
static const long ID_POWDERDIFF_SAVEHKLFCALC=               WXCRYST_ID(); 
static const long ID_POWDER_GRAPH_NEW_PATTERN=              WXCRYST_ID(); 
static const long ID_POWDERTEXTURE_MENU_ADDPHASE=                   WXCRYST_ID(); 
static const long ID_POWDERTEXTURE_MENU_DELETEPHASE=                WXCRYST_ID(); 
static const long ID_POWDERPATTERN_MENU_COMPONENTS=                 WXCRYST_ID(); 
static const long ID_POWDERPATTERN_MENU_PATTERN=                    WXCRYST_ID(); 
  
BEGIN_EVENT_TABLE(WXPowderPattern, wxWindow)
   EVT_BUTTON(ID_WXOBJ_COLLAPSE,                    WXCrystObj::OnToggleCollapse)                
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
         //mpMenuBar->AppendSeparator();
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_POWDER_MENU_IMPORT_FULLPROF,
                                 "Import Fullprof Pattern");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_POWDER_MENU_IMPORT_PSI_DMC,
                                 "Import PSI(DMC) Pattern");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_POWDER_MENU_IMPORT_ILL_D1A5,
                                 "Import ILL(D1A-D1B) Pattern (D1A5)");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_POWDER_MENU_IMPORT_XDD,
                                 "Import Xdd Pattern");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_POWDER_MENU_IMPORT_CPI,
                                 "Import Sietronics CPI Pattern");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_POWDER_MENU_IMPORT_FULLPROF4,
                                 "Import FullProf format #4");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_POWDER_MENU_IMPORT_MULTIDETECTORLLBG42,
                                 "Import Multi-Detector Format (LLB G42)");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_POWDER_MENU_IMPORT_2THETAOBSSIGMA,
                                 "Import 2Theta-Obs-Sigma Pattern");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_POWDER_MENU_IMPORT_2THETAOBS,
                                 "Import 2Theta-Obs Pattern");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_POWDER_MENU_IMPORT_TOFISISXYSIGMA,
                                 "Import ISIS TOF X Y Sigma");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_POWDER_MENU_IMPORT_GSAS,
                                 "Import GSAS Powder Data (Constant Wavelength)");
      mpMenuBar->AddMenu("Parameters",ID_REFOBJ_MENU_PAR);
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_PAR,ID_REFOBJ_MENU_PAR_FIXALL,"Fix all");
         //mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_PAR,ID_REFOBJ_MENU_PAR_UNFIXALL,"Unfix all");
      mpMenuBar->AddMenu("Components",ID_POWDERPATTERN_MENU_COMPONENTS);
         mpMenuBar->AddMenuItem(ID_POWDERPATTERN_MENU_COMPONENTS,
                                ID_POWDER_MENU_COMP_ADDBACKGD_BAYESIAN,
                                "Add Bayesian Background (automatic)");
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
      
      WXFieldPar<REAL> *pWXFieldChi2=new WXFieldPar<REAL>(this,"Chi^2",-1,&mChi2,100);
      pStats->Add(pWXFieldChi2    ,0,wxALIGN_CENTER);
      mList.Add(pWXFieldChi2);
      pWXFieldChi2->SetToolTip(_T("Chi^2=SUM[(Obs_i-Calc_i)^2/Sigma_i^2]"));
      
      WXFieldPar<REAL> *pWXFieldGof=new WXFieldPar<REAL>(this,"GoF",-1,&mGoF,70);
      pStats->Add(pWXFieldGof    ,0,wxALIGN_CENTER);
      mList.Add(pWXFieldGof);
      pWXFieldGof->SetToolTip(_T("GoF=Chi^2/NbPoints"));
      
      WXFieldPar<REAL> *pWXFieldRwp=new WXFieldPar<REAL>(this,"Rwp",-1,&mRwp,70);
      pStats->Add(pWXFieldRwp    ,0,wxALIGN_CENTER);
      mList.Add(pWXFieldRwp);
      pWXFieldRwp->SetToolTip(_T("Full profile R-factor (weighted)\n")
                              _T("Will use integrated profiles if option is set."));
      
      WXFieldPar<REAL> *pWXFieldRp=new WXFieldPar<REAL>(this,"Rp",-1,&mRp,70);
      pStats->Add(pWXFieldRp    ,0,wxALIGN_CENTER);
      mList.Add(pWXFieldRp);
      pWXFieldRp->SetToolTip(_T("Full profile R-factor (unweighted)\n")
                             _T("Will use integrated profiles if option is set."));
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
      bool val;
      if(!wxConfigBase::Get()->HasEntry("PowderPattern/BOOL/Automatically open powder pattern graph"))
         wxConfigBase::Get()->Write("PowderPattern/BOOL/Automatically open powder pattern graph", false);
      else
      {
         wxConfigBase::Get()->Read("PowderPattern/BOOL/Automatically open powder pattern graph/", &val);
         if((val)&&(mpPowderPattern->GetNbPoint()>0))
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
   VFN_DEBUG_MESSAGE("WXPowderPattern::CrystUpdate()",6)
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
   mpPowderPattern->Prepare();
   
   mChi2=mpPowderPattern->GetChi2();
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
                             "Warning : Duplicate Background !",
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
                             "Warning : Duplicate Background !",
                             wxYES_NO|wxICON_HAND|wxNO_DEFAULT);
      if(wxID_YES!=dialog.ShowModal())
         return;
   }

   long nbPointSpline=20;
   string mes="Number of Interpolation Points";
   stringstream s;
   s<<nbPointSpline;
   wxTextEntryDialog dialog(this,mes.c_str(),"Automatic Bayesian (David-Sivia) Background",
                            s.str().c_str(),wxOK | wxCANCEL);
   if(wxID_OK!=dialog.ShowModal())
   {
      VFN_DEBUG_EXIT("WXPowderPattern::OnMenuAddCompBackgdBayesian():Canceled",6)
      return;
   }
   dialog.GetValue().ToLong(&nbPointSpline);
   if(nbPointSpline<=1)nbPointSpline=2;
   
   PowderPatternBackground *pBckgd= new PowderPatternBackground;
   VFN_DEBUG_MESSAGE("WXPowderPattern::OnMenuAddCompBackgdBayesian()",6)
   mpPowderPattern->AddPowderPatternComponent(*pBckgd);
   VFN_DEBUG_MESSAGE("WXPowderPattern::OnMenuAddCompBackgdBayesian()",6)
   {
      CrystVector_REAL x(nbPointSpline),backgd(nbPointSpline);
      const CrystVector_REAL *pObs=&(pBckgd->GetParentPowderPattern().GetPowderPatternObs());
      const unsigned long nbPoint=pBckgd->GetParentPowderPattern().GetNbPoint();
      for(int i=0;i<nbPointSpline;i++)
      {
         VFN_DEBUG_MESSAGE("WXPowderPattern::OnMenuAddCompBackgdBayesian():"<<i,6)
         x(i)=pBckgd->GetParentPowderPattern().GetPowderPatternX()(i*(nbPoint-1)/(nbPointSpline-1));
         VFN_DEBUG_MESSAGE("WXPowderPattern::OnMenuAddCompBackgdBayesian():"<<i,6)
         long n1=(long)((REAL)nbPoint/(REAL)nbPointSpline*((REAL)i-0.2));
         long n2=(long)((REAL)nbPoint/(REAL)nbPointSpline*((REAL)i+0.2));
         if(n1<0) n1=0;
         if(n2>(long)nbPoint)n2=nbPoint;
         VFN_DEBUG_MESSAGE("WXPowderPattern::OnMenuAddCompBackgdBayesian():"<<i,6)
         backgd(i)=(*pObs)(n1);
         VFN_DEBUG_MESSAGE("WXPowderPattern::OnMenuAddCompBackgdBayesian():"<<i,6)
         for(long j=n1;j<n2;j++)
            if((*pObs)(j)<backgd(i))backgd(i)=(*pObs)(j);
      }
      pBckgd->SetInterpPoints(x,backgd);
   }
   VFN_DEBUG_MESSAGE("WXPowderPattern::OnMenuAddCompBackgdBayesian()",6)
   if(mpGraph!=0) mpPowderPattern->Prepare();//else this will be done when opening the graph
   VFN_DEBUG_MESSAGE("WXPowderPattern::OnMenuAddCompBackgdBayesian()",6)
   
   pBckgd->UnFixAllPar();
   pBckgd->GetOption(0).SetChoice(0);//linear
   VFN_DEBUG_MESSAGE("WXPowderPattern::OnMenuAddCompBackgdBayesian()",6)
   pBckgd->OptimizeBayesianBackground();
   VFN_DEBUG_MESSAGE("WXPowderPattern::OnMenuAddCompBackgdBayesian()",6)
   pBckgd->GetOption(0).SetChoice(1);//spline
   pBckgd->OptimizeBayesianBackground();
   VFN_DEBUG_MESSAGE("WXPowderPattern::OnMenuAddCompBackgdBayesian()",6)
   pBckgd->FixAllPar();

   //this->Layout();
   VFN_DEBUG_EXIT("WXPowderPattern::OnMenuAddCompBackgdBayesian()",6)
}

void WXPowderPattern::OnMenuAddCompCryst(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXPowderPattern::OnMenuAddCompCryst()",6)
   WXCrystValidateAllUserInput();
   PowderPatternDiffraction * diffData=new PowderPatternDiffraction;
   int choice;
   Crystal *cryst=dynamic_cast<Crystal*>
      ( WXDialogChooseFromRegistry(gCrystalRegistry,(wxWindow*)this,
         "Choose a Crystal Structure:",choice));
   if(0==cryst) return;
   diffData->SetCrystal(*cryst);
   mpPowderPattern->AddPowderPatternComponent(*diffData);
   if(mpGraph!=0) mpPowderPattern->Prepare();//else this will be done when opening the graph
   //this->Layout();
}

void WXPowderPattern::OnMenuShowGraph(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXPowderPattern::OnMenuShowGraph()"<<mpGraph,6)
   if(mpPowderPattern->GetPowderPatternObs().numElements()==0)
   {
      wxMessageDialog dumbUser(this,"Import a powder pattern or use simulation",
                               "No observed pattern !",wxOK|wxICON_EXCLAMATION);
      dumbUser.ShowModal();
      return;
   }
   if(mpGraph!=0) return;
   WXCrystValidateAllUserInput();
   mpPowderPattern->Prepare();
   wxFrame *frame= new wxFrame(this,-1,mpPowderPattern->GetName().c_str(),
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
   wxFileDialog save(this,"Choose a file","","","*.txt",wxSAVE | wxOVERWRITE_PROMPT);
   if(save.ShowModal() != wxID_OK) return;
   
   ofstream out(save.GetPath().c_str());
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
      wxTextEntryDialog dialog(this,"2Theta Min",
                              "Enter minimum 2Theta (degrees)","5",wxOK | wxCANCEL);
      if(wxID_OK!=dialog.ShowModal())
      {
         VFN_DEBUG_EXIT("WXPowderPattern::OnMenuSimulate():Cancelled",6)
         return;
      }
      dialog.GetValue().ToDouble(&min);
      if(min<1) min=1.0;
   }
   {
      wxTextEntryDialog dialog(this,"2Theta Max",
                              "Enter maximum 2Theta (degrees)","100",wxOK | wxCANCEL);
      if(wxID_OK!=dialog.ShowModal())
      {
         VFN_DEBUG_EXIT("WXPowderPattern::OnMenuSimulate():Cancelled",6)
         return;
      }
      dialog.GetValue().ToDouble(&max);
   }
   {
      wxTextEntryDialog dialog(this,"Number of points",
                              "Enter the number of points","1000",wxOK | wxCANCEL);
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
   mpPowderPattern->SetPowderPatternObs(newObs);
   VFN_DEBUG_EXIT("WXPowderPattern::OnMenuSimulate()",6)
}

void WXPowderPattern::OnMenuImportPattern(wxCommandEvent &event)
{
   VFN_DEBUG_MESSAGE("WXPowderPattern::OnMenuImportPattern()",6)
   wxFileDialog open(this,"Choose a file","","","*.*",wxOPEN | wxFILE_MUST_EXIST);
   if(open.ShowModal() != wxID_OK) return;
   if(event.GetId()==(long)ID_POWDER_MENU_IMPORT_FULLPROF)
      mpPowderPattern->ImportPowderPatternFullprof(open.GetPath().c_str());
   if(event.GetId()==(long)ID_POWDER_MENU_IMPORT_PSI_DMC)
      mpPowderPattern->ImportPowderPatternPSI_DMC(open.GetPath().c_str());
   if(event.GetId()==(long)ID_POWDER_MENU_IMPORT_ILL_D1A5)
      mpPowderPattern->ImportPowderPatternILL_D1A5(open.GetPath().c_str());
   if(event.GetId()==(long)ID_POWDER_MENU_IMPORT_XDD)
      mpPowderPattern->ImportPowderPatternXdd(open.GetPath().c_str());
   if(event.GetId()==(long)ID_POWDER_MENU_IMPORT_CPI)
      mpPowderPattern->ImportPowderPatternSietronicsCPI(open.GetPath().c_str());
   if(event.GetId()==(long)ID_POWDER_MENU_IMPORT_FULLPROF4)
      mpPowderPattern->ImportPowderPatternFullprof4(open.GetPath().c_str());
   if(event.GetId()==(long)ID_POWDER_MENU_IMPORT_MULTIDETECTORLLBG42)
      mpPowderPattern->ImportPowderPatternMultiDetectorLLBG42(open.GetPath().c_str());
   if(event.GetId()==(long)ID_POWDER_MENU_IMPORT_2THETAOBSSIGMA)
      mpPowderPattern->ImportPowderPattern2ThetaObsSigma(open.GetPath().c_str());
   if(event.GetId()==(long)ID_POWDER_MENU_IMPORT_2THETAOBS)
      mpPowderPattern->ImportPowderPattern2ThetaObs(open.GetPath().c_str());
   if(event.GetId()==(long)ID_POWDER_MENU_IMPORT_TOFISISXYSIGMA)
      mpPowderPattern->ImportPowderPatternTOF_ISIS_XYSigma(open.GetPath().c_str());

   if(event.GetId()==(long)ID_POWDER_MENU_IMPORT_GSAS)
      mpPowderPattern->ImportPowderPatternGSAS(open.GetPath().c_str());
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
      wxTextEntryDialog dialog(this,"new Wavelength)",
                              "Enter new Wavelength (Angstroems)","1",wxOK | wxCANCEL);
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
      string txt="Enter Min 2theta to exclude (degrees):";
      if(mpPowderPattern->GetRadiation().GetWavelengthType()==WAVELENGTH_TOF)
         txt="Enter Min 2theta to exclude (microseconds):";
      wxTextEntryDialog dialog(this,"Min",txt.c_str(),"0",wxOK | wxCANCEL);
      if(wxID_OK!=dialog.ShowModal())
      {
         VFN_DEBUG_EXIT("WXPowderPattern::OnMenuAddExclude():Cancelled",6)
         return;
      }
      dialog.GetValue().ToDouble(&min);
   }
   //max
   {
      string txt="Enter Max 2theta to exclude (degrees):";
      if(mpPowderPattern->GetRadiation().GetWavelengthType()==WAVELENGTH_TOF)
         txt="Enter Max 2theta to exclude (microseconds):";
      wxTextEntryDialog dialog(this,"Max",txt.c_str(),"5",wxOK | wxCANCEL);
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

void WXPowderPattern::NotifyDeleteGraph() {mpGraph=0;}
const PowderPattern& WXPowderPattern::GetPowderPattern()const
{ return *mpPowderPattern;}
void WXPowderPattern::UpdateUI(const bool lock)
{
   if(lock)mMutex.Lock();
   if(mpGraph!=0)
   {
      mpGraph->GetParent()->SetTitle(mpPowderPattern->GetName().c_str());
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

BEGIN_EVENT_TABLE(WXPowderPatternGraph, wxWindow)
   EVT_PAINT(                                   WXPowderPatternGraph::OnPaint)
   EVT_MOUSE_EVENTS(                            WXPowderPatternGraph::OnMouse)
   EVT_MENU(ID_POWDERGRAPH_MENU_UPDATE, WXPowderPatternGraph::OnUpdate)
   EVT_MENU(ID_POWDERGRAPH_MENU_TOGGLELABEL, WXPowderPatternGraph::OnToggleLabel)
   EVT_UPDATE_UI(ID_POWDER_GRAPH_NEW_PATTERN,WXPowderPatternGraph::OnRedrawNewPattern)
   EVT_CHAR(                                    WXPowderPatternGraph::OnKeyDown)
   EVT_MOUSEWHEEL(                              WXPowderPatternGraph::OnMouseWheel)
   EVT_SIZE(                                    WXPowderPatternGraph::OnSize)
END_EVENT_TABLE()

WXPowderPatternGraph::WXPowderPatternGraph(wxFrame *frame, WXPowderPattern* parent):
wxWindow(frame,-1,wxPoint(-1,-1),wxSize(-1,-1)),
mpPattern(parent),mMargin(20),mDiffPercentShift(.20),
mMaxIntensity(-1),mMinIntensity(-1),mMinX(-1),mMaxX(-1),
mpParentFrame(frame),
mIsDragging(false),mDisplayLabel(true)
{
   mpPopUpMenu=new wxMenu("Powder Pattern");
   mpPopUpMenu->Append(ID_POWDERGRAPH_MENU_UPDATE, "&Update");
   mpPopUpMenu->Append(ID_POWDERGRAPH_MENU_TOGGLELABEL, "&Hide Labels");
   if(!wxConfigBase::Get()->HasEntry("PowderPattern/BOOL/Default-display reflection indices"))
      wxConfigBase::Get()->Write("PowderPattern/BOOL/Default-display reflection indices", mDisplayLabel);
   else
      wxConfigBase::Get()->Read("PowderPattern/BOOL/Default-display reflection indices", &mDisplayLabel);
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
   wxMutexLocker mlock(mMutex);
   wxBufferedPaintDC dc(this);
   PrepareDC(dc);
   mpParentFrame->PrepareDC(dc);
   
   dc.BeginDrawing();
   
   dc.DestroyClippingRegion();
   dc.SetBackground(wxBrush("white", wxSOLID));
   dc.Clear();

   wxString fontInfo;
   dc.SetFont(*wxSMALL_FONT);

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
         xc=(wxCoord)mMargin*3;
         REAL yStep=pow((float)10,(float)floor(log10((mMaxIntensity-mMinIntensity)/nbTick)));
         yStep *= floor((mMaxIntensity-mMinIntensity)/yStep/nbTick);
         for(REAL y=yStep*ceil(mMinIntensity/yStep);y<mMaxIntensity;y+=yStep)
         {
            yc=this->Data2ScreenY(y);
            dc.DrawLine(xc-3,yc,xc+3,yc);
            fontInfo.Printf("%g",y);
            dc.GetTextExtent(fontInfo, &tmpW, &tmpH);
            dc.DrawText(fontInfo,xc-tmpW-3,yc-tmpH/2);
         }
      //X axis
         yc=(wxCoord)(height-mMargin);
         REAL xStep=pow((float)10,(float)floor(log10((mMaxX-mMinX)/nbTick)));
         xStep *= floor((mMaxX-mMinX)/xStep/nbTick);
         for(REAL x=xStep*ceil(mMinX/xStep);x<mMaxX;x+=xStep)
         {
            xc=this->Data2ScreenX(x);
            dc.DrawLine(xc,yc-3,xc,yc+3);
            fontInfo.Printf("%g",x);
            dc.GetTextExtent(fontInfo, &tmpW, &tmpH);
            dc.DrawText(fontInfo,xc-tmpW/2,yc+6);
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
   // Draw labels
   VFN_DEBUG_MESSAGE("WXPowderPatternGraph:OnPaint():5:",5)
   if(true==mDisplayLabel)
   {
      wxCoord x,y;
      wxCoord tmpW,tmpH;
      int loop=1;
      REAL yr;
      list<list<pair<const REAL ,const string > > >::const_iterator comp;
      list<pair<const REAL ,const string > >::const_iterator pos;
      unsigned int pen=0;
      for(comp=mvLabelList.begin();comp!=mvLabelList.end();++comp)
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
                  fontInfo.Printf("%s",pos->second.c_str());
                  dc.GetTextExtent(fontInfo, &tmpW, &tmpH);
                  dc.DrawText(fontInfo,x-tmpW/2,y-tmpH*(loop++)-10);
               }
               if(loop==5) loop=1;
            }
         }
      }
   }
   
   dc.EndDrawing();
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

      const REAL ttheta=this->Screen2DataX(x);
      const REAL intensity=this->Screen2DataY(y);

      wxString str;
      const long pixel=
         (long)(mpPattern->GetPowderPattern().X2PixelCorr(ttheta*DEG2RAD));
      str.Printf("X=%6.2f    ,I=%12.2f.   pixel=#%ld",ttheta,intensity,pixel);
      mpParentFrame->SetStatusText(str);

   if (event.Dragging() && event.LeftIsDown() && (!mIsDragging))
   {//Begin zooming
      mIsDragging=true;
      mDraggingX0=ttheta;
      mDraggingIntensity0=intensity;
      mMutex.Unlock();
      return;
   }
   if(event.LeftUp() && mIsDragging)
   {//Finished zooming !
      VFN_DEBUG_MESSAGE("WXPowderPatternGraph::OnMouse():Finished zooming...",5)
      mIsDragging=false;
      
      if( (fabs(ttheta-mDraggingX0)<.1) || (fabs(mDraggingIntensity0-intensity)< fabs(mMaxIntensity*.02)) )
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
      if(mDraggingX0>ttheta)
      {
         mMinX=ttheta;
         mMaxX=mDraggingX0;
      }
      else
      {
         mMinX=mDraggingX0;
         mMaxX=ttheta;
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
      return;
   }
   
   if(event.RightIsDown())
   {//popup menu
      mMutex.Unlock();
      if(mpPattern->GetPowderPattern().IsBeingRefined())
         mpPopUpMenu->Enable(ID_POWDERGRAPH_MENU_UPDATE, false);
      else
         mpPopUpMenu->Enable(ID_POWDERGRAPH_MENU_UPDATE, true);
      this->PopupMenu(mpPopUpMenu, event.GetX(), event.GetY() );
      return;
   }
   mMutex.Unlock();
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

void WXPowderPatternGraph::OnToggleLabel(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXPowderPatternGraph::OnToggleLabel()",6)
   mDisplayLabel = !mDisplayLabel;
   this->Refresh(false);
   if(mDisplayLabel) mpPopUpMenu->SetLabel(ID_POWDERGRAPH_MENU_TOGGLELABEL, "Hide Labels");
   else mpPopUpMenu->SetLabel(ID_POWDERGRAPH_MENU_TOGGLELABEL, "Show Labels");
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
         const REAL range=mMaxIntensity-mMinIntensity;
         mMinIntensity+=range/8;
         mMaxIntensity+=range/8;
         break;
      }
      case(WXK_DOWN):
      {
         const REAL range=mMaxIntensity-mMinIntensity;
         mMinIntensity-=range/8;
         mMaxIntensity-=range/8;
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
}

void WXPowderPatternGraph::OnSize(wxSizeEvent& event)
{
   this->Refresh(false);
}

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
   mMaxIntensity=mCalc.max();
   mMinIntensity=mCalc.min();
   const float max=mObs.max();
   const float min=mObs.min();
   if(max>mMaxIntensity) mMaxIntensity=max;
   if(min<mMinIntensity) mMinIntensity=min;
   if(mMinIntensity<0) mMinIntensity=0;
   mMaxX=mX.max();
   mMinX=mX.min();
   mClockAxisLimits.Click();
   VFN_DEBUG_MESSAGE("WXPowderPatternGraph::ResetAxisLimits():"<<mMinIntensity<<","<<mMaxIntensity<<","<<mMinX<<","<<mMaxX,10)
}
long WXPowderPatternGraph::Data2ScreenX(const REAL x)const
{
   wxCoord width,height;
   this->GetSize(&width, &height);
   return (long)(mMargin*3+(x-mMinX)*(width-3*mMargin)/(mMaxX-mMinX));
}
long WXPowderPatternGraph::Point2ScreenX(const long x)const
{
   wxCoord width,height;
   this->GetSize(&width, &height);
   return (long)(mMargin*3+(mX(x)-mMinX)*(width-3*mMargin)/(REAL)(mMaxX-mMinX));
}
long WXPowderPatternGraph::Data2ScreenY(const REAL y)const
{
   wxCoord width,height;
   this->GetSize(&width, &height);
   return (long)(height-mMargin-(y-mMinIntensity)*(height-2*mMargin)
                        /(mMaxIntensity-mMinIntensity));
}
REAL WXPowderPatternGraph::Screen2DataX(const long x)const
{
   wxCoord width,height;
   this->GetSize(&width, &height);
   return mMinX+(x-mMargin*3)*(mMaxX-mMinX)/(REAL)(width-3*mMargin);
}
REAL WXPowderPatternGraph::Screen2DataY(const long y)const
{
   wxCoord width,height;
   this->GetSize(&width, &height);
   return mMinIntensity+(height-mMargin-y)*(mMaxIntensity-mMinIntensity)/(REAL)(height-2*mMargin);
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
      mpGridBackgroundPoint->SetColLabelValue(0,"Position");
      mpGridBackgroundPoint->SetColLabelValue(1,"Intensity");
      mpGridBackgroundPoint->AutoSizeRows();
      mpSizer->Add(mpGridBackgroundPoint,0,wxALIGN_LEFT);
   mpTopSizer->SetSizeHints(this);
   this->Layout();
   this->CrystUpdate(true);
}
void WXPowderPatternBackground::OnMenuImportUserBackground(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXPowderPatternBackground::OnMenuImportUserBackground()",6)
   wxFileDialog *open= new wxFileDialog(this,"Choose background file with 2Theta Ibackgd",
                                        "","","*.*",wxOPEN | wxFILE_MUST_EXIST);
   if(open->ShowModal() != wxID_OK) return;
   mpPowderPatternBackground->ImportUserBackground(open->GetPath().c_str());
   open->Destroy();
}
void WXPowderPatternBackground::OnMenuOptimizeBayesianBackground(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_ENTRY("WXPowderPatternBackground::OnMenuOptimizeBayesianBackground()",6)
   mpPowderPatternBackground->UnFixAllPar();
   mpPowderPatternBackground->OptimizeBayesianBackground();
   mpPowderPatternBackground->FixAllPar();
   VFN_DEBUG_EXIT("WXPowderPatternBackground::OnMenuOptimizeBayesianBackground()",6)
}
void WXPowderPatternBackground::OnMenuAutomaticBayesianBackground(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_ENTRY("WXPowderPatternBackground::OnMenuAutomaticBayesianBackground()",6)
   WXCrystValidateAllUserInput();

   long nbPointSpline=20;
   string mes="Number of Interpolation Points";
   stringstream s;
   s<<nbPointSpline;
   wxTextEntryDialog dialog(this,mes.c_str(),"Automatic Bayesian (David-Sivia) Background",
                            s.str().c_str(),wxOK | wxCANCEL);
   if(wxID_OK!=dialog.ShowModal())
   {
      VFN_DEBUG_EXIT("WXPowderPatternBackground::OnMenuAutomaticBayesianBackground():Canceled",6)
      return;
   }
   dialog.GetValue().ToLong(&nbPointSpline);
   if(nbPointSpline<=1) nbPointSpline=1;
   {
      CrystVector_REAL tth(nbPointSpline),backgd(nbPointSpline);
      const CrystVector_REAL *pObs=&(mpPowderPatternBackground->GetParentPowderPattern().GetPowderPatternObs());
      const unsigned long nbPoint=mpPowderPatternBackground->GetParentPowderPattern().GetNbPoint();
      for(int i=0;i<nbPointSpline;i++)
      {
         if(i==(nbPointSpline-1))
            tth(i)=mpPowderPatternBackground->GetParentPowderPattern()
                     .GetPowderPatternX()(nbPoint-1);
         else tth(i)=mpPowderPatternBackground->GetParentPowderPattern()
                     .GetPowderPatternX()((i*nbPoint)/(nbPointSpline-1));
         long n1=(long)((REAL)nbPoint/(REAL)nbPointSpline*((REAL)i-0.2));
         long n2=(long)((REAL)nbPoint/(REAL)nbPointSpline*((REAL)i+0.2));
         if(n1<0) n1=0;
         if(n2>(long)nbPoint)n2=nbPoint;
         backgd(i)=(*pObs)(n1);
         for(long j=n1;j<n2;j++)
            if((*pObs)(j)<backgd(i))backgd(i)=(*pObs)(j);
      }
      mpPowderPatternBackground->SetInterpPoints(tth,backgd);
   }
   //mpPowderPatternBackground->GetParentPowderPattern().Prepare();
   mpPowderPatternBackground->UnFixAllPar();
   mpPowderPatternBackground->GetOption(0).SetChoice(0);//linear
   mpPowderPatternBackground->OptimizeBayesianBackground();
   mpPowderPatternBackground->GetOption(0).SetChoice(1);//spline
   mpPowderPatternBackground->OptimizeBayesianBackground();
   mpPowderPatternBackground->FixAllPar();

   VFN_DEBUG_EXIT("WXPowderPatternBackground::OnMenuAutomaticBayesianBackground()",6)
}
void WXPowderPatternBackground::OnEditGridBackgroundPoint(wxGridEvent &e)
{
   if(mIsSelfUpdating) return;
   VFN_DEBUG_ENTRY("WXPowderPatternBackground::OnEditGridBackgroundPoint():"<<e.GetRow()<<","<<e.GetCol(),10)
   const long r=e.GetRow();
   const long c=e.GetCol();
   wxString s=mpGridBackgroundPoint->GetCellValue(r,c);
   if(s!="")
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
      if(diff!=0)
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
         tmp.Printf("%f",f*mBackgroundInterpPointX(i));
         mpGridBackgroundPoint->SetCellValue(i,0,tmp);
         tmp.Printf("%f",mBackgroundInterpPointIntensity(i));
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
static const long ID_POWDERDIFF_PROFILE_DEPV=                  WXCRYST_ID();

BEGIN_EVENT_TABLE(WXPowderPatternDiffraction, wxWindow)
   EVT_BUTTON(ID_POWDERDIFF_CRYSTAL,WXPowderPatternDiffraction::OnChangeCrystal)
   EVT_MENU(ID_POWDERDIFF_SAVEHKLFCALC, 
                                            WXPowderPatternDiffraction::OnMenuSaveHKLFcalc)
   EVT_MENU(ID_POWDERDIFF_PROFILE_PV,       WXPowderPatternDiffraction::OnChangeProfile)
   EVT_MENU(ID_POWDERDIFF_PROFILE_DEPV,     WXPowderPatternDiffraction::OnChangeProfile)
   EVT_UPDATE_UI(ID_CRYST_UPDATEUI,         WXRefinableObj::OnUpdateUI)
END_EVENT_TABLE()

WXPowderPatternDiffraction::WXPowderPatternDiffraction(wxWindow *parent,
                                                         PowderPatternDiffraction *p):
WXRefinableObj(parent,p),mpPowderPatternDiffraction(p)
{
   VFN_DEBUG_ENTRY("WXPowderPatternDiffraction::WXPowderPatternDiffraction()",6)
   mpWXTitle->SetForegroundColour(wxColour(0,255,0));
    //Menu
      mpMenuBar->AddMenu("Object",ID_REFOBJ_MENU_OBJ);
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_POWDERDIFF_SAVEHKLFCALC,
                                "Save HKL Fcalc");
      mpMenuBar->AddMenu("Profile",ID_POWDERDIFF_PROFILE);
         mpMenuBar->AddMenuItem(ID_POWDERDIFF_PROFILE,ID_POWDERDIFF_PROFILE_PV,
                                "Pseudo-Voigt (X-Ray & monochromatic neutron)");
         mpMenuBar->AddMenuItem(ID_POWDERDIFF_PROFILE,ID_POWDERDIFF_PROFILE_DEPV,
                                "Double-Exponential Pseudo-Voigt (neutron TOF)");
      mpSizer->SetItemMinSize(mpMenuBar,
                              mpMenuBar->GetSize().GetWidth(),
                              mpMenuBar->GetSize().GetHeight());
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
   wxFileDialog save(this,"Choose a file to save to","","","*.txt",wxSAVE | wxOVERWRITE_PROMPT);
   if(save.ShowModal() != wxID_OK) return;
   
   ofstream out(save.GetPath().c_str());
   if(!out) return;//:TODO:
   mpPowderPatternDiffraction->PrintFhklCalc(out);
   out.close();
}
void WXPowderPatternDiffraction::UpdateUI(const bool lock)
{
   if(lock) mMutex.Lock();
   mpFieldCrystal->SetValue(mpPowderPatternDiffraction->GetCrystal().GetName());
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
         mpPowderPatternDiffraction->mpReflectionProfile= new ReflectionProfilePseudoVoigt;
         add=true;
      }
      else
         if(mpPowderPatternDiffraction->mpReflectionProfile->GetClassName()
            !="ReflectionProfilePseudoVoigt")
         {
            delete mpPowderPatternDiffraction->mpReflectionProfile;
            mpPowderPatternDiffraction->mpReflectionProfile= new ReflectionProfilePseudoVoigt;
            add=true;
         }
   }
   if(event.GetId()==ID_POWDERDIFF_PROFILE_DEPV)
   {
      if(mpPowderPatternDiffraction->mpReflectionProfile==0)
      {
         mpPowderPatternDiffraction->mpReflectionProfile
            = new ReflectionProfileDoubleExponentialPseudoVoigt
                     (mpPowderPatternDiffraction->GetCrystal());
         add=true;
      }
      else
         if(mpPowderPatternDiffraction->mpReflectionProfile->GetClassName()
            !="ReflectionProfileDoubleExponentialPseudoVoigt")
         {
            delete mpPowderPatternDiffraction->mpReflectionProfile;
            mpPowderPatternDiffraction->mpReflectionProfile
               = new ReflectionProfileDoubleExponentialPseudoVoigt
                        (mpPowderPatternDiffraction->GetCrystal());
            add=true;
         }
   }
   if(add)
   {
      mpPowderPatternDiffraction->mpReflectionProfile->WXCreate(this);
      mList.Add(mpPowderPatternDiffraction->mpReflectionProfile->WXGet());
      mpSizer->Add(mpPowderPatternDiffraction->mpReflectionProfile->WXGet());
      this->BottomLayout(mpPowderPatternDiffraction->mpReflectionProfile->WXGet());
      this->CrystUpdate(true);
   }
   VFN_DEBUG_EXIT("WXPowderPatternDiffraction::OnChangeProfile()",6)
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

