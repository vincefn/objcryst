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

#include "wx/wx.h"
#include "wxCryst/wxPowderPattern.h"
#include "wxCryst/wxRadiation.h"

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
////////////////////////////////////////////////////////////////////////
//
//    WXRadiation
//
////////////////////////////////////////////////////////////////////////
BEGIN_EVENT_TABLE(WXRadiation, wxWindow)
   EVT_UPDATE_UI(ID_CRYST_UPDATEUI, WXRadiation::UpdateUI)
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
#if 1
   WXFieldRefPar* pFieldWavelength=new WXFieldRefPar(this,"Wavelength:",
                                     &(mpRadiation->GetPar(mpRadiation->mWavelength.data())));
#else
   WXCrystObjBasic* pFieldWavelength
      =mpRadiation->GetPar(mpRadiation->mWavelength.data()).WXCreate(this);
#endif
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
      
   this->CrystUpdate();
   this->SetSizer(mpSizer);
   mpSizer->Layout();
   mpSizer->Fit(this);
   VFN_DEBUG_EXIT("WXRadiation::WXRadiation()",6)
}

void WXRadiation::CrystUpdate()
{
   mList.CrystUpdate();
   if(true==wxThread::IsMain()) this->UpdateUI();
   else
   {
      wxUpdateUIEvent event(ID_CRYST_UPDATEUI);
      wxPostEvent(this,event);
   }
}
void WXRadiation::UpdateUI()
{
   mList.UpdateUI();
}
void WXRadiation::OnUpdateUI(wxUpdateUIEvent& event)
{
   this->UpdateUI();
}

////////////////////////////////////////////////////////////////////////
//
//    WXPowderPattern
//
////////////////////////////////////////////////////////////////////////
BEGIN_EVENT_TABLE(WXPowderPattern, wxWindow)
   EVT_BUTTON(ID_WXOBJ_COLLAPSE,                       WXCrystObj::OnToggleCollapse)
   EVT_MENU(ID_REFOBJ_MENU_OBJ_SAVE,                   WXRefinableObj::OnMenuSave)
   EVT_MENU(ID_REFOBJ_MENU_OBJ_LOAD,                   WXRefinableObj::OnMenuLoad)
   EVT_MENU(ID_REFOBJ_MENU_PAR_FIXALL,                 WXRefinableObj::OnMenuFixAllPar)
   EVT_MENU(ID_REFOBJ_MENU_PAR_UNFIXALL,               WXRefinableObj::OnMenuUnFixAllPar)
   EVT_MENU(ID_POWDERSPECTRUM_MENU_SCATT_ADDCOMPBACKGD,WXPowderPattern::OnMenuAddCompBackgd)
   EVT_MENU(ID_POWDERSPECTRUM_MENU_SCATT_ADDCOMPCRYST, WXPowderPattern::OnMenuAddCompCryst)
   EVT_MENU(ID_POWDERSPECTRUM_MENU_SAVETEXT,           WXPowderPattern::OnMenuSaveText)
   EVT_MENU(ID_POWDERSPECTRUM_MENU_SIMULATE,           WXPowderPattern::OnMenuSimulate)
   EVT_MENU(ID_POWDERSPECTRUM_MENU_IMPORT_FULLPROF,    WXPowderPattern::OnMenuImportFullProf)
   EVT_MENU(ID_POWDERSPECTRUM_MENU_IMPORT_PSI_DMC,     WXPowderPattern::OnMenuImportPSI)
   EVT_MENU(ID_POWDERSPECTRUM_MENU_IMPORT_ILL_D1A5,    WXPowderPattern::OnMenuImportILL)
   EVT_MENU(ID_POWDERSPECTRUM_MENU_IMPORT_XDD,         WXPowderPattern::OnMenuImportXdd)
   EVT_MENU(ID_POWDERSPECTRUM_MENU_IMPORT_CPI,         WXPowderPattern::OnMenuImportCPI)
   EVT_MENU(ID_POWDERSPECTRUM_MENU_IMPORT_2THETAOBSSIGMA,
                                                WXPowderPattern::OnMenuImport2ThetaObsSigma)
   EVT_MENU(ID_POWDERSPECTRUM_MENU_IMPORT_2THETAOBS,   WXPowderPattern::OnMenuImport2ThetaObs)
   EVT_MENU(ID_POWDERSPECTRUM_MENU_WAVELENGTH_SET,     WXPowderPattern::OnMenuSetWavelength)
   EVT_MENU(ID_POWDERSPECTRUM_MENU_WAVELENGTH_XRAY,    WXPowderPattern::OnMenuSetWavelength)
   EVT_MENU(ID_POWDERSPECTRUM_MENU_WAVELENGTH_NEUTRON, WXPowderPattern::OnMenuSetWavelength)
   EVT_MENU(ID_POWDERSPECTRUM_MENU_WAVELENGTH_SET_AG,  WXPowderPattern::OnMenuSetWavelength)
   EVT_MENU(ID_POWDERSPECTRUM_MENU_WAVELENGTH_SET_MO,  WXPowderPattern::OnMenuSetWavelength)
   EVT_MENU(ID_POWDERSPECTRUM_MENU_WAVELENGTH_SET_CU,  WXPowderPattern::OnMenuSetWavelength)
   EVT_MENU(ID_POWDERSPECTRUM_MENU_WAVELENGTH_SET_FE,  WXPowderPattern::OnMenuSetWavelength)
   EVT_MENU(ID_POWDERSPECTRUM_MENU_WAVELENGTH_SET_CR,  WXPowderPattern::OnMenuSetWavelength)
   EVT_MENU(ID_POWDERSPECTRUM_MENU_WAVELENGTH_SET_AGA1,WXPowderPattern::OnMenuSetWavelength)
   EVT_MENU(ID_POWDERSPECTRUM_MENU_WAVELENGTH_SET_MOA1,WXPowderPattern::OnMenuSetWavelength)
   EVT_MENU(ID_POWDERSPECTRUM_MENU_WAVELENGTH_SET_CUA1,WXPowderPattern::OnMenuSetWavelength)
   EVT_MENU(ID_POWDERSPECTRUM_MENU_WAVELENGTH_SET_FEA1,WXPowderPattern::OnMenuSetWavelength)
   EVT_MENU(ID_POWDERSPECTRUM_MENU_WAVELENGTH_SET_CRA1,WXPowderPattern::OnMenuSetWavelength)
   EVT_MENU(ID_POWDERSPECTRUM_MENU_GRAPH,              WXPowderPattern::OnMenuShowGraph)
   EVT_MENU(ID_POWDERSPECTRUM_MENU_FITSCALE_R,         WXPowderPattern::OnMenuFitScaleForR)
   EVT_MENU(ID_POWDERSPECTRUM_MENU_FITSCALE_RW,        WXPowderPattern::OnMenuFitScaleForRw)
   EVT_MENU(ID_POWDERSPECTRUM_MENU_ADD_2THETA_EXCLUDE, WXPowderPattern::OnMenuAdd2ThetaExclude)
   EVT_UPDATE_UI(ID_CRYST_UPDATEUI,                    WXRefinableObj::OnUpdateUI)
END_EVENT_TABLE()

WXPowderPattern::WXPowderPattern(wxWindow *parent, PowderPattern* pow):
WXRefinableObj(parent,pow),mpPowderPattern(pow),mpGraph(0)
{
   VFN_DEBUG_MESSAGE("WXPowderPattern::WXPowderPattern()",6)
   mpWXTitle->SetForegroundColour(wxColour(255,0,0));
   // Menu
      mpMenuBar->AddMenu("Data",ID_REFOBJ_MENU_OBJ);
         //:TODO: reactivate & test those menus
         //mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_REFOBJ_MENU_OBJ_SAVE,"Save");
         //mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_REFOBJ_MENU_OBJ_LOAD,"Load");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_POWDERSPECTRUM_MENU_SAVETEXT,
                                "Save spectrum (text)");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_POWDERSPECTRUM_MENU_SIMULATE,
                                "Simulation mode (no obs. spectrum)");
         //mpMenuBar->AppendSeparator();
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_POWDERSPECTRUM_MENU_IMPORT_FULLPROF,
                                 "Import Fullprof Pattern");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_POWDERSPECTRUM_MENU_IMPORT_PSI_DMC,
                                 "Import PSI(DMC) Pattern");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_POWDERSPECTRUM_MENU_IMPORT_ILL_D1A5,
                                 "Import ILL(D1A-D1B) Pattern (D1A5)");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_POWDERSPECTRUM_MENU_IMPORT_XDD,
                                 "Import Xdd Pattern");
        
mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_POWDERSPECTRUM_MENU_IMPORT_CPI,
                                 "Import Sietronics CPI Pattern");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_POWDERSPECTRUM_MENU_IMPORT_2THETAOBSSIGMA,
                                 "Import 2Theta-Obs-Sigma Pattern");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_POWDERSPECTRUM_MENU_IMPORT_2THETAOBS,
                                 "Import 2Theta-Obs Pattern");
      mpMenuBar->AddMenu("Parameters",ID_REFOBJ_MENU_PAR);
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_PAR,ID_REFOBJ_MENU_PAR_FIXALL,"Fix all");
         //mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_PAR,ID_REFOBJ_MENU_PAR_UNFIXALL,"Unfix all");
      mpMenuBar->AddMenu("Components",ID_CRYSTAL_MENU_SCATT);
         mpMenuBar->AddMenuItem(ID_CRYSTAL_MENU_SCATT,
                                ID_POWDERSPECTRUM_MENU_SCATT_ADDCOMPBACKGD,
                                "Add Interpolated Background");
         mpMenuBar->AddMenuItem(ID_CRYSTAL_MENU_SCATT,ID_POWDERSPECTRUM_MENU_SCATT_ADDCOMPCRYST,
                                "Add Crystalline Phase");
      mpMenuBar->AddMenu("Radiation",ID_POWDERSPECTRUM_MENU_WAVELENGTH);
         mpMenuBar->AddMenuItem(ID_POWDERSPECTRUM_MENU_WAVELENGTH,
                                ID_POWDERSPECTRUM_MENU_WAVELENGTH_NEUTRON,
                                "Neutron");
         mpMenuBar->AddMenuItem(ID_POWDERSPECTRUM_MENU_WAVELENGTH,
                                ID_POWDERSPECTRUM_MENU_WAVELENGTH_XRAY,
                                "X-Rays");
         mpMenuBar->AddMenuItem(ID_POWDERSPECTRUM_MENU_WAVELENGTH,
                                ID_POWDERSPECTRUM_MENU_WAVELENGTH_SET,
                                "Monochromatic Wavelength");
         mpMenuBar->AddMenuItem(ID_POWDERSPECTRUM_MENU_WAVELENGTH,
                                ID_POWDERSPECTRUM_MENU_WAVELENGTH_SET_AG,
                                "X-Ray Tube Ag Ka12");
         mpMenuBar->AddMenuItem(ID_POWDERSPECTRUM_MENU_WAVELENGTH,
                                ID_POWDERSPECTRUM_MENU_WAVELENGTH_SET_AGA1,
                                "X-Ray Tube Ag Ka1");
         mpMenuBar->AddMenuItem(ID_POWDERSPECTRUM_MENU_WAVELENGTH,
                                ID_POWDERSPECTRUM_MENU_WAVELENGTH_SET_MO,
                                "X-Ray Tube Mo Ka12");
         mpMenuBar->AddMenuItem(ID_POWDERSPECTRUM_MENU_WAVELENGTH,
                                ID_POWDERSPECTRUM_MENU_WAVELENGTH_SET_MOA1,
                                "X-Ray Tube Mo Ka1");
         mpMenuBar->AddMenuItem(ID_POWDERSPECTRUM_MENU_WAVELENGTH,
                                ID_POWDERSPECTRUM_MENU_WAVELENGTH_SET_CU,
                                "X-Ray Tube Cu Ka12");
         mpMenuBar->AddMenuItem(ID_POWDERSPECTRUM_MENU_WAVELENGTH,
                                ID_POWDERSPECTRUM_MENU_WAVELENGTH_SET_CUA1,
                                "X-Ray Tube Cu Ka1");
         mpMenuBar->AddMenuItem(ID_POWDERSPECTRUM_MENU_WAVELENGTH,
                                ID_POWDERSPECTRUM_MENU_WAVELENGTH_SET_FE,
                                "X-Ray Tube Fe Ka12");
         mpMenuBar->AddMenuItem(ID_POWDERSPECTRUM_MENU_WAVELENGTH,
                                ID_POWDERSPECTRUM_MENU_WAVELENGTH_SET_FEA1,
                                "X-Ray Tube Fe Ka1");
         mpMenuBar->AddMenuItem(ID_POWDERSPECTRUM_MENU_WAVELENGTH,
                                ID_POWDERSPECTRUM_MENU_WAVELENGTH_SET_CR,
                                "X-Ray Tube Cr Ka12");
         mpMenuBar->AddMenuItem(ID_POWDERSPECTRUM_MENU_WAVELENGTH,
                                ID_POWDERSPECTRUM_MENU_WAVELENGTH_SET_CRA1,
                                "X-Ray Tube Cr Ka1");
      mpMenuBar->AddMenu("Pattern",ID_CRYSTAL_MENU_DISPLAY);
         mpMenuBar->AddMenuItem(ID_CRYSTAL_MENU_DISPLAY,ID_POWDERSPECTRUM_MENU_GRAPH,
                                "Show Graph");
         mpMenuBar->AddMenuItem(ID_CRYSTAL_MENU_DISPLAY,ID_POWDERSPECTRUM_MENU_FITSCALE_R,
                                "Fit Scale for R");
         mpMenuBar->AddMenuItem(ID_CRYSTAL_MENU_DISPLAY,ID_POWDERSPECTRUM_MENU_FITSCALE_RW,
                                "Fit Scale for Rw");
         mpMenuBar->AddMenuItem(ID_CRYSTAL_MENU_DISPLAY,
                                ID_POWDERSPECTRUM_MENU_ADD_2THETA_EXCLUDE,
                                "Add 2Theta excluded region");
   //Radiation
      mpSizer->Add(mpPowderPattern->mRadiation.WXCreate(this),0);
      mList.Add(mpPowderPattern->mRadiation.WXGet());
   // Correction to 2Theta
      wxBoxSizer* thetaCorrSizer=new wxBoxSizer(wxHORIZONTAL);
#if 1
      WXFieldRefPar* fieldThetaZero    =new WXFieldRefPar(this,"2theta zero:",
                                   &(mpPowderPattern
                                     ->GetPar(&(mpPowderPattern->m2ThetaZero))),90 );
      WXFieldRefPar* fieldThetaDispl    =new WXFieldRefPar(this,"2theta displacement:",
                                   &(mpPowderPattern
                                     ->GetPar(&(mpPowderPattern->m2ThetaDisplacement))),90 );
      WXFieldRefPar* fieldThetaTransp    =new WXFieldRefPar(this,"2theta transparency:",
                                   &(mpPowderPattern
                                     ->GetPar(&(mpPowderPattern->m2ThetaTransparency))),90 );
#else
      WXCrystObjBasic* fieldThetaZero    
         =mpPowderPattern->GetPar(&(mpPowderPattern->m2ThetaZero)).WXCreate(this);
      WXCrystObjBasic* fieldThetaDispl
         =mpPowderPattern->GetPar(&(mpPowderPattern->m2ThetaDisplacement)).WXCreate(this);
      WXCrystObjBasic* fieldThetaTransp
         =mpPowderPattern->GetPar(&(mpPowderPattern->m2ThetaTransparency)).WXCreate(this);
#endif
      thetaCorrSizer->Add(fieldThetaZero,0);
      thetaCorrSizer->Add(fieldThetaDispl,0);
      thetaCorrSizer->Add(fieldThetaTransp,0);
      mList.Add(fieldThetaZero);
      mList.Add(fieldThetaDispl);
      mList.Add(fieldThetaTransp);
      mpSizer->Add(thetaCorrSizer);
   // Max Sin(theta/Lambda)
      WXFieldPar<REAL> *maxSiThOvLa=
         new WXFieldPar<REAL>(this,"Max Sin(theta)/lambda:",-1,
                              &(mpPowderPattern->mMaxSinThetaOvLambda));
      mpSizer->Add(maxSiThOvLa,0,wxALIGN_LEFT);
      mList.Add(maxSiThOvLa);
   // Components
      mpWXComponent=mpPowderPattern
                    ->mPowderPatternComponentRegistry.WXCreate(this);
      mpSizer->Add(mpWXComponent,0,wxALIGN_LEFT);
      mList.Add(mpWXComponent);
   
   VFN_DEBUG_MESSAGE("WXPowderPattern::WXPowderPattern():1",6)
   this->CrystUpdate();
   this->Layout();
   VFN_DEBUG_MESSAGE("WXPowderPattern::WXPowderPattern():End",6)
}

void WXPowderPattern::CrystUpdate()
{
   VFN_DEBUG_MESSAGE("WXPowderPattern::CrystUpdate()",6)
   WXCrystValidateAllUserInput();
   this->WXRefinableObj::CrystUpdate();
   
   // Will force re-generating reflection list if the wavelength,
   // or lattice par, or the spacegroup has changed.
   mpPowderPattern->Prepare();
   
   if(mpGraph!=0)
   {
      mpGraph->SetPattern( mpPowderPattern->GetPowderPatternObs(),
                           mpPowderPattern->GetPowderPatternCalc(),
                           mpPowderPattern->Get2ThetaMin(),
                           mpPowderPattern->Get2ThetaStep());
   }
} 

void WXPowderPattern::OnMenuAddCompBackgd(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXPowderPattern::OnMenuAddCompBackgd()",6)
   WXCrystValidateAllUserInput();
   PowderPatternBackground *backgdData= new PowderPatternBackground;
   mpPowderPattern->AddPowderPatternComponent(*backgdData);
   if(mpGraph!=0) mpPowderPattern->Prepare();//else this will be done when opening the graph
   this->Layout();
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
   this->Layout();
}

void WXPowderPattern::OnMenuShowGraph(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXPowderPattern::OnMenuShowGraph()"<<mpGraph,6)
   if(mpGraph!=0) return;
   WXCrystValidateAllUserInput();
   mpPowderPattern->Prepare();
   wxFrame *frame= new wxFrame(this,-1,mpPowderPattern->GetName().c_str());
   mpGraph = new WXPowderPatternGraph(frame,this);
   frame->Show(true);
   frame->CreateStatusBar(2);
   this->CrystUpdate();
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
                              "Enter minimum 2Theta (degrees)","0",wxOK | wxCANCEL);
      if(wxID_OK!=dialog.ShowModal())
      {
         VFN_DEBUG_EXIT("WXPowderPattern::OnMenuSimulate():Cancelled",6)
         return;
      }
      dialog.GetValue().ToDouble(&min);
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

void WXPowderPattern::OnMenuImportFullProf(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXPowderPattern::OnMenuImportFullProf()",6)
   wxFileDialog *open= new wxFileDialog(this,"Choose a file","","","*.*",
                                        wxOPEN | wxFILE_MUST_EXIST);
   if(open->ShowModal() != wxID_OK) return;
   
   mpPowderPattern->ImportPowderPatternFullprof(open->GetPath().c_str());
   open->Destroy();
}

void WXPowderPattern::OnMenuImportPSI(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXPowderPattern::OnMenuImportPSI()",6)
   wxFileDialog *open= new wxFileDialog(this,"Choose a file","","","*.*",
                                        wxOPEN | wxFILE_MUST_EXIST);
   if(open->ShowModal() != wxID_OK) return;
   
   mpPowderPattern->ImportPowderPatternPSI_DMC(open->GetPath().c_str());
   open->Destroy();
}

void WXPowderPattern::OnMenuImportILL(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXPowderPattern::OnMenuImportPSI()",6)
   wxFileDialog *open= new wxFileDialog(this,"Choose a file","","","*.*",
                                        wxOPEN | wxFILE_MUST_EXIST);
   if(open->ShowModal() != wxID_OK) return;
   
   mpPowderPattern->ImportPowderPatternILL_D1A5(open->GetPath().c_str());
   open->Destroy();
}

void WXPowderPattern::OnMenuImportXdd(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXPowderPattern::OnMenuImportXdd()",6)
   wxFileDialog *open= new wxFileDialog(this,"Choose a file","","","*.xdd",
                                        wxOPEN | wxFILE_MUST_EXIST);
   if(open->ShowModal() != wxID_OK) return;
   
   mpPowderPattern->ImportPowderPatternXdd(open->GetPath().c_str());
   open->Destroy();
}

void WXPowderPattern::OnMenuImportCPI(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXPowderPattern::OnMenuImportCPI()",6)
   wxFileDialog *open= new wxFileDialog(this,"Choose a file","","","*.cpi",
                                        wxOPEN | wxFILE_MUST_EXIST);
   if(open->ShowModal() != wxID_OK) return;
   
   mpPowderPattern->ImportPowderPatternSietronicsCPI(open->GetPath().c_str());
   open->Destroy();
}

void WXPowderPattern::OnMenuImport2ThetaObsSigma(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXPowderPattern::OnMenuImport2ThetaObsSigma()",6)
   wxFileDialog *open= new wxFileDialog(this,"Choose a file","","","*.*",
                                        wxOPEN | wxFILE_MUST_EXIST);
   if(open->ShowModal() != wxID_OK) return;
   
   mpPowderPattern->ImportPowderPattern2ThetaObsSigma(open->GetPath().c_str());
   open->Destroy();
}
void WXPowderPattern::OnMenuImport2ThetaObs(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXPowderPattern::OnMenuImport2ThetaObs()",6)
   wxFileDialog *open= new wxFileDialog(this,"Choose a file","","","*.*",
                                        wxOPEN | wxFILE_MUST_EXIST);
   if(open->ShowModal() != wxID_OK) return;
   
   mpPowderPattern->ImportPowderPattern2ThetaObs(open->GetPath().c_str());
   open->Destroy();
}
void WXPowderPattern::OnMenuFitScaleForR(wxCommandEvent & WXUNUSED(event))
{
   if(0==mpGraph) return;
   WXCrystValidateAllUserInput();
   mpPowderPattern->FitScaleFactorForR();//FitScaleFactorForIntegratedR
   this->CrystUpdate();
}

void WXPowderPattern::OnMenuFitScaleForRw(wxCommandEvent & WXUNUSED(event))
{
   if(0==mpGraph) return;
   WXCrystValidateAllUserInput();
   mpPowderPattern->FitScaleFactorForRw();//FitScaleFactorForIntegratedRw
   this->CrystUpdate();
}


void WXPowderPattern::OnMenuSetWavelength(wxCommandEvent & event)
{
   WXCrystValidateAllUserInput();
   // this looks stupid. In fact, if a user changed the wavelength in the
   // corresponding field, this is (unfortunately) not applied to the
   // components automagically. So we need this function to do the job...
   switch(event.GetId())
   {
      case ID_POWDERSPECTRUM_MENU_WAVELENGTH_XRAY:
         mpPowderPattern->SetRadiationType(RAD_XRAY);break;
      case ID_POWDERSPECTRUM_MENU_WAVELENGTH_NEUTRON:
         mpPowderPattern->SetRadiationType(RAD_NEUTRON);break;
      case ID_POWDERSPECTRUM_MENU_WAVELENGTH_SET:
      {
         double lambda;
         wxTextEntryDialog dialog(this,"new Wavelength:",
                                 "Enter new Wavelength (Angstroems)","1",wxOK | wxCANCEL);
         if(wxID_OK!=dialog.ShowModal())
         {
            VFN_DEBUG_EXIT("WXPowderPattern::OnMenuSetWavelength():Monochromatic:Cancelled",6)
            return;
         }
         dialog.GetValue().ToDouble(&lambda);
         mpPowderPattern->SetWavelength(lambda);break;
      }
      case ID_POWDERSPECTRUM_MENU_WAVELENGTH_SET_AG:
         mpPowderPattern->SetWavelength("Ag");break;
      case ID_POWDERSPECTRUM_MENU_WAVELENGTH_SET_MO:
         mpPowderPattern->SetWavelength("Mo");break;
      case ID_POWDERSPECTRUM_MENU_WAVELENGTH_SET_CU:
         mpPowderPattern->SetWavelength("Cu");break;
      case ID_POWDERSPECTRUM_MENU_WAVELENGTH_SET_FE:
         mpPowderPattern->SetWavelength("Fe");break;
      case ID_POWDERSPECTRUM_MENU_WAVELENGTH_SET_CR:
         mpPowderPattern->SetWavelength("Cr");break;
      case ID_POWDERSPECTRUM_MENU_WAVELENGTH_SET_AGA1:
         mpPowderPattern->SetWavelength("AgA1");break;
      case ID_POWDERSPECTRUM_MENU_WAVELENGTH_SET_MOA1:
         mpPowderPattern->SetWavelength("MoA1");break;
      case ID_POWDERSPECTRUM_MENU_WAVELENGTH_SET_CUA1:
         mpPowderPattern->SetWavelength("CuA1");break;
      case ID_POWDERSPECTRUM_MENU_WAVELENGTH_SET_FEA1:
         mpPowderPattern->SetWavelength("FeA1");break;
      case ID_POWDERSPECTRUM_MENU_WAVELENGTH_SET_CRA1:
         mpPowderPattern->SetWavelength("CrA1");break;
   }
   this->CrystUpdate();
}

void WXPowderPattern::OnMenuAdd2ThetaExclude(wxCommandEvent & WXUNUSED(event))
{
   WXCrystValidateAllUserInput();
   double min,max;
   //min
   {
      wxTextEntryDialog dialog(this,"Min 2Theta",
                              "Enter Min 2Theta to exclude (degrees):","0",wxOK | wxCANCEL);
      if(wxID_OK!=dialog.ShowModal())
      {
         VFN_DEBUG_EXIT("WXPowderPattern::OnMenuAdd2ThetaExclude():Cancelled",6)
         return;
      }
      dialog.GetValue().ToDouble(&min);
   }
   //max
   {
      wxTextEntryDialog dialog(this,"Max 2Theta",
                              "Enter Max 2Theta to exclude (degrees):","5",wxOK | wxCANCEL);
      if(wxID_OK!=dialog.ShowModal())
      {
         VFN_DEBUG_EXIT("WXPowderPattern::OnMenuAdd2ThetaExclude():Cancelled",6)
         return;
      }
      dialog.GetValue().ToDouble(&max);
   }
   if(max<min)
   {
      VFN_DEBUG_EXIT("WXPowderPattern::OnMenuAdd2ThetaExclude():Stupid user.",6)
      return;
   }
   mpPowderPattern->Add2ThetaExcludedRegion(min*DEG2RAD,max*DEG2RAD);
}

void WXPowderPattern::NotifyDeleteGraph() {mpGraph=0;}
const PowderPattern& WXPowderPattern::GetPowderPattern()const
{ return *mpPowderPattern;}
void WXPowderPattern::UpdateUI()
{
   if(mpGraph!=0)
   {
      mpGraph->GetParent()->SetTitle(mpPowderPattern->GetName().c_str());
   }
   this->WXRefinableObj::UpdateUI();
}
////////////////////////////////////////////////////////////////////////
//
//    WXPowderPatternGraph
//
////////////////////////////////////////////////////////////////////////
BEGIN_EVENT_TABLE(WXPowderPatternGraph, wxWindow)
   EVT_PAINT(                                   WXPowderPatternGraph::OnPaint)
   EVT_MOUSE_EVENTS(                            WXPowderPatternGraph::OnMouse)
   EVT_MENU(ID_POWDERSPECTRUMGRAPH_MENU_UPDATE, WXPowderPatternGraph::OnUpdate)
   EVT_UPDATE_UI(ID_POWDERSPECTRUM_GRAPH_NEW_PATTERN,WXPowderPatternGraph::OnRedrawNewPattern)
END_EVENT_TABLE()

WXPowderPatternGraph::WXPowderPatternGraph(wxFrame *frame, WXPowderPattern* parent):
wxWindow(frame,-1,wxPoint(-1,-1),wxSize(-1,-1),wxRETAINED),
mpPattern(parent),mMargin(50),mDiffPercentShift(.20),mpParentFrame(frame),
mCalcPatternIsLocked(false),mIsDragging(false)
{
   mpPopUpMenu=new wxMenu("Crystal");
   mpPopUpMenu->Append(ID_POWDERSPECTRUMGRAPH_MENU_UPDATE, "&Update");
   mpPattern->CrystUpdate();
}

WXPowderPatternGraph::~WXPowderPatternGraph()
{
   mpPattern->NotifyDeleteGraph();
}

void WXPowderPatternGraph::OnPaint(wxPaintEvent& WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXPowderPatternGraph:OnPaint()",5)
   if(true==mCalcPatternIsLocked) return;
   wxPaintDC dc(this);
   PrepareDC(dc);
   mpParentFrame->PrepareDC(dc);
   
   dc.BeginDrawing();
   
   dc.DestroyClippingRegion();
   dc.SetBackground(wxBrush("white", wxSOLID));
   dc.Clear();

   wxString fontInfo;
   //fontInfo.Printf("Powder Pattern : %s",mpData->GetName().c_str());
   //dc.DrawText(fontInfo, 5, 5);

   // Get Window Size
   wxCoord width,height;
   this->GetSize(&width, &height);
   //const int margin=mMargin;
   //width -= margin;
   //height -= margin;
   VFN_DEBUG_MESSAGE("WXPowderPatternGraph:OnPaint():1",5)

   //Check spectrum is not being updated
   while(mCalcPatternIsLocked) wxUsleep(10);
   mCalcPatternIsLocked=true;
   VFN_DEBUG_MESSAGE("WXPowderPatternGraph:OnPaint():2:"<<mObs.numElements(),5)

   const long nbPoints=mLast-mFirst+1;
   VFN_DEBUG_MESSAGE("WXPowderPatternGraph:OnPaint():3",5)

   // Draw Axis (sort of)
   {
      wxCoord tmpW,tmpH;
      dc.SetPen(* wxBLACK_PEN);
      dc.DrawLine(mMargin,height-mMargin,mMargin,mMargin);
      dc.DrawLine(mMargin,height-mMargin,width,height-mMargin);
      const int nbTick=10;//approximate
      wxCoord xc,yc;
      //Y axis
         xc=(wxCoord)mMargin;
         REAL yStep=pow(10,floor(log10((mMaxIntensity-mMinIntensity)/nbTick)));
         yStep *= floor((mMaxIntensity-mMinIntensity)/yStep/nbTick);
         for(REAL y=yStep*floor(mMinIntensity/yStep);y<mMaxIntensity;y+=yStep)
         {
            yc=(wxCoord) (height-mMargin-(y-mMinIntensity)*(height-2*mMargin)
                           /(mMaxIntensity-mMinIntensity));
            dc.DrawLine(xc-3,yc,xc+3,yc);
            fontInfo.Printf("%g",y);
            dc.GetTextExtent(fontInfo, &tmpW, &tmpH);
            dc.DrawText(fontInfo,xc-tmpW,yc-tmpH/2);
         }
      //X axis
         yc=(wxCoord)(height-mMargin);
         REAL xStep=pow(10,floor(log10((mMax2Theta-mMin2Theta)/nbTick)));
         xStep *= floor((mMax2Theta-mMin2Theta)/xStep/nbTick);
         for(REAL x=xStep*floor(mMin2Theta/xStep);x<mMax2Theta;x+=xStep)
         {
            xc=(wxCoord)(mMargin+(x-mMin2Theta)*(width-mMargin)/(mMax2Theta-mMin2Theta));
            dc.DrawLine(xc,yc-3,xc,yc+3);
            fontInfo.Printf("%g",x);
            dc.GetTextExtent(fontInfo, &tmpW, &tmpH);
            dc.DrawText(fontInfo,xc-tmpW/2,yc+tmpH);
         }
   }
   // Draw observed spectrum
   {
      dc.SetPen(* wxCYAN_PEN);
      wxCoord x1,y1,x2,y2;
      x2=(wxCoord)(mMargin+ 0. *(width-mMargin)/(REAL)nbPoints);
      y2=(wxCoord)(height-mMargin-(mObs(mFirst)-mMinIntensity)*(height-2*mMargin)
                     /(mMaxIntensity-mMinIntensity));
      for(long i=mFirst+1;i<=mLast;i++)
      {
         x1=x2;
         y1=y2;
         x2=(wxCoord)(mMargin+ (i-mFirst)*(width-mMargin)/(REAL)nbPoints);
         y2=(wxCoord)(height-mMargin-(mObs(i)-mMinIntensity)*(height-2*mMargin)
                        /(mMaxIntensity-mMinIntensity));
         dc.DrawLine(x1,y1,x2,y2);
      }
   }

   // Draw calculated spectrum
   {
      dc.SetPen(* wxRED_PEN);
      wxCoord x1,y1,x2,y2;
      x2=(wxCoord)(mMargin+ 0.     *(width-mMargin)/(REAL)nbPoints);
      y2=(wxCoord)(height-mMargin- (mCalc(mFirst)-mMinIntensity)*(height-2*mMargin)
                     /(mMaxIntensity-mMinIntensity));
      for(long i=mFirst+1;i<=mLast;i++)
      {
         x1=x2;
         y1=y2;
         x2=(wxCoord)(mMargin+ (i-mFirst)*(width-mMargin)/(REAL)nbPoints);
         y2=(wxCoord)(height-mMargin-(mCalc(i)-mMinIntensity)*(height-2*mMargin)
                        /(mMaxIntensity-mMinIntensity));
         dc.DrawLine(x1,y1,x2,y2);
      }
   }
   mCalcPatternIsLocked=false;
   dc.EndDrawing();

   VFN_DEBUG_MESSAGE("WXPowderPatternGraph:OnPaint():End",5)
}
void WXPowderPatternGraph::OnMouse(wxMouseEvent &event)
{
   if(true==mCalcPatternIsLocked)
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

   if((x>width)||(y>height)) return;
   //cout <<pos.x<<" "<<pos.y<<" "<<x<<" "<<y<<" "<<width<<" "<<height<<endl;

      const REAL 
         ttheta=mMin2Theta+(x-mMargin)*(mMax2Theta-mMin2Theta)/(REAL)(width-mMargin);
      const REAL intensity=mMinIntensity+(height-mMargin-y)*(mMaxIntensity-mMinIntensity)
                                             /(REAL)(height-2*mMargin);

      wxString str;
      const long pixel=
         mpPattern->GetPowderPattern().Get2ThetaCorrPixel(ttheta*DEG2RAD);
      str.Printf("2Theta=%6.2f    ,I=%12.2f.   pixel=#%d",ttheta,intensity,pixel);
      mpParentFrame->SetStatusText(str);

   if (event.Dragging() && event.LeftIsDown() && (!mIsDragging))
   {//Begin zooming
      mIsDragging=true;
      mDragging2Theta0=ttheta;
      mDraggingIntensity0=intensity;
      return;
   }
   if(event.LeftUp() && mIsDragging)
   {//Finished zooming !
      VFN_DEBUG_MESSAGE("WXPowderPatternGraph::OnMouse():Finished zooming...",5)
      mIsDragging=false;
      
      if( (abs(ttheta-mDragging2Theta0)<.3) || (abs(mDraggingIntensity0-intensity)< abs(mMaxIntensity*.05)) )
      {
         return;
      }
      if(mDraggingIntensity0>intensity)
      {
         if(mDraggingIntensity0<0.) return;
         mMinIntensity=intensity;
         mMaxIntensity=mDraggingIntensity0;
      }
      else
      {
         if(intensity<0.) return;
         mMinIntensity=mDraggingIntensity0;
         mMaxIntensity=intensity;
      }
      if(mDragging2Theta0>ttheta)
      {
         mMin2Theta=ttheta;
         mMax2Theta=mDragging2Theta0;
      }
      else
      {
         mMin2Theta=mDragging2Theta0;
         mMax2Theta=ttheta;
      }
      const long nbpoints=m2theta.numElements();
      bool flag=true;
      for(long i=0;i<nbpoints;i++)
      {
         if(flag) if(m2theta(i)>=mMin2Theta) {mFirst=i;flag=false;}
         if(m2theta(i)>=mMax2Theta) {mLast=i;break;}
      }
      if(mFirst<0) mFirst=0;
      if(mLast>=nbpoints) mLast=nbpoints-1;
      wxUpdateUIEvent event(ID_POWDERSPECTRUM_GRAPH_NEW_PATTERN);
      wxPostEvent(this,event);
      return;
   }

   if(false==event.Dragging()) mIsDragging=false;

   if(event.LeftDClick())
   {//Reset axis range
      this->ResetAxisLimits();
      wxUpdateUIEvent event(ID_POWDERSPECTRUM_GRAPH_NEW_PATTERN);
      wxPostEvent(this,event);
      return;
   }
   
   if(event.RightIsDown())
   {//popup menu
      this->PopupMenu(mpPopUpMenu, event.GetX(), event.GetY() );
      return;
   }
}

void WXPowderPatternGraph::OnUpdate(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXPowderPatternGraph::OnUpdate()",6)
   mpPattern->CrystUpdate();
}

void WXPowderPatternGraph::SetPattern(const CrystVector_REAL &obs,
                                        const CrystVector_REAL &calc,
                                        const REAL tthetaMin,const REAL tthetaStep)
{
   VFN_DEBUG_MESSAGE("WXPowderPatternGraph::SetPattern()",5)
   //Make sure spectrum is not being used (for drawing)
   while(mCalcPatternIsLocked) wxUsleep(1);
   mCalcPatternIsLocked=true;
   mCalc=calc;
   mCalcPatternIsLocked=false;
   mObs=obs;
   const long nbPoint=mObs.numElements();
   m2theta.resize(nbPoint);
   for(long i=0;i<nbPoint;i++) m2theta(i)=tthetaMin+i*tthetaStep;
   m2theta*=RAD2DEG;
   this->ResetAxisLimits();
   // If we only send an OnPaint event, only the parts which have been erased are redrawn
   // (under windows). SO we must force the complete Refresh of the window... in the
   // main thread of course...
   wxUpdateUIEvent event(ID_POWDERSPECTRUM_GRAPH_NEW_PATTERN);
   wxPostEvent(this,event);
   VFN_DEBUG_MESSAGE("WXPowderPatternGraph::SetPattern():End",5)
}

void WXPowderPatternGraph::OnRedrawNewPattern(wxUpdateUIEvent& WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXPowderPatternGraph::SetPattern()",5)
   this->Refresh(false);
}
void WXPowderPatternGraph::ResetAxisLimits()
{
   mMaxIntensity=mCalc.max();
   mMinIntensity=mCalc.min();
   const float max=mObs.max();
   const float min=mObs.min();
   if(max>mMaxIntensity) mMaxIntensity=max;
   if(min<mMinIntensity) mMinIntensity=min;
   if(mMinIntensity<0) mMinIntensity=0;
   mMax2Theta=m2theta.max();
   mMin2Theta=m2theta.min();
   mFirst=0;
   mLast=m2theta.numElements()-1;
}

////////////////////////////////////////////////////////////////////////
//
//    WXPowderPatternBackgound
//
////////////////////////////////////////////////////////////////////////
BEGIN_EVENT_TABLE(WXPowderPatternBackground, wxWindow)
   EVT_MENU(ID_POWDERSPECTRUMBACKGROUND_IMPORT, 
                     WXPowderPatternBackground::OnMenuImportUserBackground)
END_EVENT_TABLE()

WXPowderPatternBackground::WXPowderPatternBackground(wxWindow *parent, 
                                                     PowderPatternBackground *b):
WXRefinableObj(parent,b),mpPowderPatternBackground(b)
{
   mpWXTitle->SetForegroundColour(wxColour(0,255,0));
   //Menu
      mpMenuBar->AddMenu("Object",ID_REFOBJ_MENU_OBJ);
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_POWDERSPECTRUMBACKGROUND_IMPORT,"Import");
   
   //:TODO: Display points
   this->Layout();
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

////////////////////////////////////////////////////////////////////////
//
//    WXPowderPatternDiffraction
//
////////////////////////////////////////////////////////////////////////
BEGIN_EVENT_TABLE(WXPowderPatternDiffraction, wxWindow)
   EVT_BUTTON(ID_POWDERSPECTRUMDIFF_CRYSTAL,WXPowderPatternDiffraction::OnChangeCrystal)
   EVT_MENU(ID_POWDERSPECTRUMDIFF_SAVEHKLFCALC, 
                                            WXPowderPatternDiffraction::OnMenuSaveHKLFcalc)
   EVT_UPDATE_UI(ID_CRYST_UPDATEUI,         WXRefinableObj::OnUpdateUI)
END_EVENT_TABLE()

WXPowderPatternDiffraction::WXPowderPatternDiffraction(wxWindow *parent,
                                                         PowderPatternDiffraction *p):
WXRefinableObj(parent,p),mpPowderPatternDiffraction(p)
{
   mpWXTitle->SetForegroundColour(wxColour(0,255,0));
    //Menu
      mpMenuBar->AddMenu("Object",ID_REFOBJ_MENU_OBJ);
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_POWDERSPECTRUMDIFF_SAVEHKLFCALC,
                                "Save HKL Fcalc");
    // Crystal Choice
      mpFieldCrystal=new WXFieldChoice(this,ID_POWDERSPECTRUMDIFF_CRYSTAL,"Crystal:",300);
      mpSizer->Add(mpFieldCrystal,0,wxALIGN_LEFT);
      mList.Add(mpFieldCrystal);
   //Profile Parameters
      wxBoxSizer* profileSizer=new wxBoxSizer(wxHORIZONTAL);
#if 1
      WXFieldRefPar* pFieldCagliotiU    =new WXFieldRefPar(this,"U:",
                                   &(mpPowderPatternDiffraction
                                     ->GetPar(&(mpPowderPatternDiffraction->mCagliotiU))),90 );
      WXFieldRefPar* pFieldCagliotiV    =new WXFieldRefPar(this,"V:",
                                   &(mpPowderPatternDiffraction
                                     ->GetPar(&(mpPowderPatternDiffraction->mCagliotiV))),90 );
      WXFieldRefPar* pFieldCagliotiW    =new WXFieldRefPar(this,"W:",
                                   &(mpPowderPatternDiffraction
                                     ->GetPar(&(mpPowderPatternDiffraction->mCagliotiW))),90 );
      WXFieldRefPar* pFieldEta0         =new WXFieldRefPar(this,"Eta0:",
                                 &(mpPowderPatternDiffraction
                                  ->GetPar(&(mpPowderPatternDiffraction->mPseudoVoigtEta0))));
      WXFieldRefPar* pFieldEta1         =new WXFieldRefPar(this,"Eta1:",
                                 &(mpPowderPatternDiffraction
                                  ->GetPar(&(mpPowderPatternDiffraction->mPseudoVoigtEta1))));
#else
      WXCrystObjBasic* pFieldCagliotiU
         =mpPowderPatternDiffraction->GetPar(&(mpPowderPatternDiffraction->mCagliotiU))
            .WXCreate(this);
      WXCrystObjBasic* pFieldCagliotiV
         =mpPowderPatternDiffraction->GetPar(&(mpPowderPatternDiffraction->mCagliotiV))
            .WXCreate(this);
      WXCrystObjBasic* pFieldCagliotiW
         =mpPowderPatternDiffraction->GetPar(&(mpPowderPatternDiffraction->mCagliotiW))
            .WXCreate(this);
      WXCrystObjBasic* pFieldEta0
         =mpPowderPatternDiffraction->GetPar(&(mpPowderPatternDiffraction->mPseudoVoigtEta0))
            .WXCreate(this);
      WXCrystObjBasic* pFieldEta1=
         mpPowderPatternDiffraction->GetPar(&(mpPowderPatternDiffraction->mPseudoVoigtEta1))
            .WXCreate(this);
#endif
      profileSizer->Add(pFieldCagliotiU,0);
      profileSizer->Add(pFieldCagliotiV,0);
      profileSizer->Add(pFieldCagliotiW,0);
      profileSizer->Add(pFieldEta0,0);
      profileSizer->Add(pFieldEta1,0);
      mList.Add(pFieldCagliotiU);
      mList.Add(pFieldCagliotiV);
      mList.Add(pFieldCagliotiW);
      mList.Add(pFieldEta0);
      mList.Add(pFieldEta1);
      mpSizer->Add(profileSizer);
   //Global Biso factor
      WXCrystObjBasic* fieldGlobalBiso
         =mpPowderPatternDiffraction->GetPar(&(mpPowderPatternDiffraction->mGlobalBiso))
            .WXCreate(this);
      mList.Add(fieldGlobalBiso);
      mpSizer->Add(fieldGlobalBiso);
      
   this->CrystUpdate();
   this->Layout();
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
   this->CrystUpdate();
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
void WXPowderPatternDiffraction::UpdateUI()
{
   mpFieldCrystal->SetValue(mpPowderPatternDiffraction->GetCrystal().GetName());
   this->WXRefinableObj::UpdateUI();
}

}// namespace 

