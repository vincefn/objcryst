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

// wx headers, with or without precompilation
#include "wx/wxprec.h"
#ifdef __BORLANDC__
    #pragma hdrstop
#endif
#ifndef WX_PRECOMP
    #include "wx/wx.h"
#endif

#include "ObjCryst/wxCryst/wxDiffractionSingleCrystal.h"
#include "ObjCryst/wxCryst/wxRadiation.h"

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
//    WXDiffractionSingleCrystalGraph
//
////////////////////////////////////////////////////////////////////////
/** Class to display Iobs and Icalc as a function of 1/d
*
*/
class WXDiffractionSingleCrystalGraph:public WXMultiGraph
{
   public:
      WXDiffractionSingleCrystalGraph(wxFrame *frame,WXDiffractionSingleCrystal *parent):
      WXMultiGraph(frame)
      {
         mpParent=parent;
      }
      virtual ~WXDiffractionSingleCrystalGraph()
      {
         mpParent->NotifyDeleteGraph();
      }
   private:
      WXDiffractionSingleCrystal *mpParent;
};
////////////////////////////////////////////////////////////////////////
//
//    WXDiffractionSingleCrystal
//
////////////////////////////////////////////////////////////////////////
static long ID_DIFFSINGLECRYST_MENU_SAVEHKLIOBSICALC=      WXCRYST_ID(); 
static long ID_DIFFSINGLECRYST_MENU_SAVEHKLFCALC=          WXCRYST_ID(); 
static long ID_DIFFSINGLECRYST_MENU_SIMULATE=              WXCRYST_ID(); 
static long ID_DIFFSINGLECRYST_MENU_IMPORT_HKLIOBS=        WXCRYST_ID(); 
static long ID_DIFFSINGLECRYST_MENU_IMPORT_HKLIOBSSIGMA=   WXCRYST_ID(); 
static long ID_DIFFSINGLECRYST_MENU_IMPORT_SHELXHKLF4=     WXCRYST_ID(); 
static long ID_DIFFSINGLECRYST_MENU_IMPORT_CIF=            WXCRYST_ID(); 
static long ID_DIFFSINGLECRYST_MENU_IMPORT_JANAM91=        WXCRYST_ID(); 
static long ID_DIFFSINGLECRYST_MENU_IMPORT_HKLIOBSGROUP=   WXCRYST_ID(); 
static long ID_DIFFSINGLECRYST_MENU_FITSCALE_R=            WXCRYST_ID(); 
static long ID_DIFFSINGLECRYST_MENU_FITSCALE_RW=           WXCRYST_ID(); 
static long ID_DIFFSINGLECRYST_MENU_WAVELENGTH=            WXCRYST_ID(); 
static long ID_DIFFSINGLECRYST_MENU_WAVELENGTH_XRAY=       WXCRYST_ID(); 
static long ID_DIFFSINGLECRYST_MENU_WAVELENGTH_NEUTRON=    WXCRYST_ID(); 
static long ID_DIFFSINGLECRYST_MENU_WAVELENGTH_ELECTRON=   WXCRYST_ID(); 
static long ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET=        WXCRYST_ID(); 
static long ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_AG=     WXCRYST_ID(); 
static long ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_MO=     WXCRYST_ID(); 
static long ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_CU=     WXCRYST_ID(); 
static long ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_FE=     WXCRYST_ID(); 
static long ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_CO=     WXCRYST_ID(); 
static long ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_CR=     WXCRYST_ID(); 
static long ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_AGA1=   WXCRYST_ID(); 
static long ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_MOA1=   WXCRYST_ID(); 
static long ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_CUA1=   WXCRYST_ID(); 
static long ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_FEA1=   WXCRYST_ID(); 
static long ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_COA1=   WXCRYST_ID(); 
static long ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_CRA1=   WXCRYST_ID(); 
static long ID_DIFFSINGLECRYST_CRYSTAL=                    WXCRYST_ID(); 
static long ID_DIFFSINGLECRYST_MENU_DATA=                  WXCRYST_ID();
static long ID_DIFFSINGLECRYST_MENU_DATA_GRAPH=            WXCRYST_ID();
  
BEGIN_EVENT_TABLE(WXDiffractionSingleCrystal, wxWindow)
   EVT_BUTTON(ID_WXOBJ_COLLAPSE,                        WXCrystObj::OnToggleCollapse)
   EVT_MENU(ID_REFOBJ_MENU_OBJ_SAVE,                    WXRefinableObj::OnMenuSave)
   EVT_MENU(ID_REFOBJ_MENU_OBJ_LOAD,                    WXRefinableObj::OnMenuLoad)
   EVT_MENU(ID_REFOBJ_MENU_PAR_FIXALL,                  WXRefinableObj::OnMenuFixAllPar)
   EVT_MENU(ID_REFOBJ_MENU_PAR_UNFIXALL,                WXRefinableObj::OnMenuUnFixAllPar)
   EVT_MENU(ID_DIFFSINGLECRYST_MENU_SAVEHKLIOBSICALC,   WXDiffractionSingleCrystal::OnMenuSaveHKLIobsIcalc)
   EVT_MENU(ID_DIFFSINGLECRYST_MENU_SAVEHKLFCALC,       WXDiffractionSingleCrystal::OnMenuSaveHKLFcalc)
   EVT_MENU(ID_DIFFSINGLECRYST_MENU_SIMULATE,           WXDiffractionSingleCrystal::OnMenuSimulate)
   EVT_MENU(ID_DIFFSINGLECRYST_MENU_IMPORT_HKLIOBS,     WXDiffractionSingleCrystal::OnMenuImport)
   EVT_MENU(ID_DIFFSINGLECRYST_MENU_IMPORT_HKLIOBSSIGMA,WXDiffractionSingleCrystal::OnMenuImport)
   EVT_MENU(ID_DIFFSINGLECRYST_MENU_IMPORT_SHELXHKLF4,  WXDiffractionSingleCrystal::OnMenuImport)
   EVT_MENU(ID_DIFFSINGLECRYST_MENU_IMPORT_CIF       ,  WXDiffractionSingleCrystal::OnMenuImport)
   EVT_MENU(ID_DIFFSINGLECRYST_MENU_IMPORT_JANAM91,     WXDiffractionSingleCrystal::OnMenuImport)
   EVT_MENU(ID_DIFFSINGLECRYST_MENU_IMPORT_HKLIOBSGROUP,WXDiffractionSingleCrystal::OnMenuImport)
   EVT_BUTTON(ID_DIFFSINGLECRYST_CRYSTAL,               WXDiffractionSingleCrystal::OnChangeCrystal)
   //EVT_MENU(ID_DIFFSINGLECRYST_MENU_FITSCALE_R,         WXDiffractionSingleCrystal::OnMenuFitScaleForR)
   //EVT_MENU(ID_DIFFSINGLECRYST_MENU_FITSCALE_RW,        WXDiffractionSingleCrystal::OnMenuFitScaleForRw)
   EVT_MENU(ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET,     WXDiffractionSingleCrystal::OnMenuSetWavelength)
   EVT_MENU(ID_DIFFSINGLECRYST_MENU_WAVELENGTH_XRAY,    WXDiffractionSingleCrystal::OnMenuSetWavelength)
   EVT_MENU(ID_DIFFSINGLECRYST_MENU_WAVELENGTH_NEUTRON, WXDiffractionSingleCrystal::OnMenuSetWavelength)
   EVT_MENU(ID_DIFFSINGLECRYST_MENU_WAVELENGTH_ELECTRON,WXDiffractionSingleCrystal::OnMenuSetWavelength)
   EVT_MENU(ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_AG,  WXDiffractionSingleCrystal::OnMenuSetWavelength)
   EVT_MENU(ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_MO,  WXDiffractionSingleCrystal::OnMenuSetWavelength)
   EVT_MENU(ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_CU,  WXDiffractionSingleCrystal::OnMenuSetWavelength)
   EVT_MENU(ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_FE,  WXDiffractionSingleCrystal::OnMenuSetWavelength)
   EVT_MENU(ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_CO,  WXDiffractionSingleCrystal::OnMenuSetWavelength)
   EVT_MENU(ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_CR,  WXDiffractionSingleCrystal::OnMenuSetWavelength)
   EVT_MENU(ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_AGA1,WXDiffractionSingleCrystal::OnMenuSetWavelength)
   EVT_MENU(ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_MOA1,WXDiffractionSingleCrystal::OnMenuSetWavelength)
   EVT_MENU(ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_CUA1,WXDiffractionSingleCrystal::OnMenuSetWavelength)
   EVT_MENU(ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_FEA1,WXDiffractionSingleCrystal::OnMenuSetWavelength)
   EVT_MENU(ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_COA1,WXDiffractionSingleCrystal::OnMenuSetWavelength)
   EVT_MENU(ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_CRA1,WXDiffractionSingleCrystal::OnMenuSetWavelength)
   EVT_MENU(ID_DIFFSINGLECRYST_MENU_DATA_GRAPH,         WXDiffractionSingleCrystal::OnMenuShowGraph)
   EVT_UPDATE_UI(ID_CRYST_UPDATEUI,                     WXRefinableObj::OnUpdateUI)
END_EVENT_TABLE()

WXDiffractionSingleCrystal::WXDiffractionSingleCrystal(wxWindow *parent,
                                                       DiffractionDataSingleCrystal* data):
WXRefinableObj(parent,data),mpData(data),mpGraph(0),mGrapIdObs(0),mGrapIdCalc(0)
{
   VFN_DEBUG_MESSAGE("WXDiffractionSingleCrystal::WXDiffractionSingleCrystal()",6)
   mpWXTitle->SetForegroundColour(wxColour(255,0,0));
   mpWXTitle->SetSize(400,-1);
   // Menu
      mpMenuBar->AddMenu("File",ID_REFOBJ_MENU_OBJ);
         //:TODO: reactivate & test those menus
         //mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_REFOBJ_MENU_OBJ_SAVE,"Save");
         //mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_REFOBJ_MENU_OBJ_LOAD,"Load");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_DIFFSINGLECRYST_MENU_SAVEHKLIOBSICALC,
                                "Save HKL Iobs Icalc (text)");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_DIFFSINGLECRYST_MENU_SAVEHKLFCALC,
                                "Save HKL Fcalc (text)");
         mpMenuBar->GetMenu(ID_REFOBJ_MENU_OBJ).AppendSeparator();
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_DIFFSINGLECRYST_MENU_SIMULATE,
                                "Simulation mode (generate HKL list)");
         mpMenuBar->GetMenu(ID_REFOBJ_MENU_OBJ).AppendSeparator();
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_DIFFSINGLECRYST_MENU_IMPORT_HKLIOBS,
                                 "Import HKL Iobs");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_DIFFSINGLECRYST_MENU_IMPORT_HKLIOBSSIGMA,
                                 "Import HKL Iobs Sigma (space or tab-separated)");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_DIFFSINGLECRYST_MENU_IMPORT_SHELXHKLF4,
                                 "Import HKL Iobs Sigma (HKLF 4 Shelx format)");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_DIFFSINGLECRYST_MENU_IMPORT_CIF,
                                 "Import CIF single crystal data");
         //mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_DIFFSINGLECRYST_MENU_IMPORT_JANAM91,
         //                        "Import Jana M91");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_DIFFSINGLECRYST_MENU_IMPORT_HKLIOBSGROUP,
                                 "Import Reflections with group intensity");
        
      mpMenuBar->AddMenu("Radiation",ID_DIFFSINGLECRYST_MENU_WAVELENGTH);
         mpMenuBar->AddMenuItem(ID_DIFFSINGLECRYST_MENU_WAVELENGTH,
                                ID_DIFFSINGLECRYST_MENU_WAVELENGTH_NEUTRON,
                                "Neutron");
         mpMenuBar->AddMenuItem(ID_DIFFSINGLECRYST_MENU_WAVELENGTH,
                                ID_DIFFSINGLECRYST_MENU_WAVELENGTH_XRAY,
                                "X-Ray");
         mpMenuBar->AddMenuItem(ID_DIFFSINGLECRYST_MENU_WAVELENGTH,
                                ID_DIFFSINGLECRYST_MENU_WAVELENGTH_ELECTRON,
                                "Electron");
         mpMenuBar->AddMenuItem(ID_DIFFSINGLECRYST_MENU_WAVELENGTH,
                                ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET,
                                "Monochromatic Wavelength");
         mpMenuBar->AddMenuItem(ID_DIFFSINGLECRYST_MENU_WAVELENGTH,
                                ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_AG,
                                "X-Ray Tube Ag Ka12");
         mpMenuBar->AddMenuItem(ID_DIFFSINGLECRYST_MENU_WAVELENGTH,
                                ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_AGA1,
                                "X-Ray Tube Ag Ka1");
         mpMenuBar->AddMenuItem(ID_DIFFSINGLECRYST_MENU_WAVELENGTH,
                                ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_MO,
                                "X-Ray Tube Mo Ka12");
         mpMenuBar->AddMenuItem(ID_DIFFSINGLECRYST_MENU_WAVELENGTH,
                                ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_MOA1,
                                "X-Ray Tube Mo Ka1");
         mpMenuBar->AddMenuItem(ID_DIFFSINGLECRYST_MENU_WAVELENGTH,
                                ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_CU,
                                "X-Ray Tube Cu Ka12");
         mpMenuBar->AddMenuItem(ID_DIFFSINGLECRYST_MENU_WAVELENGTH,
                                ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_CUA1,
                                "X-Ray Tube Cu Ka1");
         mpMenuBar->AddMenuItem(ID_DIFFSINGLECRYST_MENU_WAVELENGTH,
                                ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_FE,
                                "X-Ray Tube Fe Ka12");
         mpMenuBar->AddMenuItem(ID_DIFFSINGLECRYST_MENU_WAVELENGTH,
                                ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_FEA1,
                                "X-Ray Tube Fe Ka1");
         mpMenuBar->AddMenuItem(ID_DIFFSINGLECRYST_MENU_WAVELENGTH,
                                ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_CO,
                                "X-Ray Tube Co Ka12");
         mpMenuBar->AddMenuItem(ID_DIFFSINGLECRYST_MENU_WAVELENGTH,
                                ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_COA1,
                                "X-Ray Tube Co Ka1");
         mpMenuBar->AddMenuItem(ID_DIFFSINGLECRYST_MENU_WAVELENGTH,
                                ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_CR,
                                "X-Ray Tube Cr Ka12");
         mpMenuBar->AddMenuItem(ID_DIFFSINGLECRYST_MENU_WAVELENGTH,
                                ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_CRA1,
                                "X-Ray Tube Cr Ka1");
      //mpMenuBar->AddMenu("Compute",ID_CRYSTAL_MENU_DISPLAY);
      //   mpMenuBar->AddMenuItem(ID_CRYSTAL_MENU_DISPLAY,ID_DIFFSINGLECRYST_MENU_FITSCALE_R,
      //                          "Fit Scale for R");
      //   mpMenuBar->AddMenuItem(ID_CRYSTAL_MENU_DISPLAY,ID_DIFFSINGLECRYST_MENU_FITSCALE_RW,
      //                          "Fit Scale for Rw");
      mpMenuBar->AddMenu("Data",ID_DIFFSINGLECRYST_MENU_DATA);
         mpMenuBar->AddMenuItem(ID_DIFFSINGLECRYST_MENU_DATA,ID_DIFFSINGLECRYST_MENU_DATA_GRAPH,
                                "Show Graph");
      mpSizer->SetItemMinSize(mpMenuBar,
                              mpMenuBar->GetSize().GetWidth(),
                              mpMenuBar->GetSize().GetHeight());
   //Radiation
      mpSizer->Add(mpData->mRadiation.WXCreate(this),0);
      mList.Add(mpData->mRadiation.WXGet());
   // Crystal
      mpFieldCrystal=new WXFieldChoice(this,ID_DIFFSINGLECRYST_CRYSTAL,"Crystal:",300);
      mpSizer->Add(mpFieldCrystal,0,wxALIGN_LEFT);
      mList.Add(mpFieldCrystal);
   // Max Sin(theta/Lambda)
      WXFieldPar<REAL> *maxSiThOvLa=
         new WXFieldPar<REAL>(this,"Max Sin(theta)/lambda:",-1,&(mpData->mMaxSinThetaOvLambda));
      mpSizer->Add(maxSiThOvLa,0,wxALIGN_LEFT);
      mList.Add(maxSiThOvLa);
   // Statistics
      wxBoxSizer* pStats=new wxBoxSizer(wxHORIZONTAL);
      
      WXFieldPar<REAL> *pWXFieldChi2=new WXFieldPar<REAL>(this,"Chi^2",-1,&mChi2,100);
      pStats->Add(pWXFieldChi2    ,0,wxALIGN_CENTER);
      mList.Add(pWXFieldChi2);
      
      WXFieldPar<REAL> *pWXFieldGof=new WXFieldPar<REAL>(this,"GoF",-1,&mGoF,70);
      pStats->Add(pWXFieldGof    ,0,wxALIGN_CENTER);
      mList.Add(pWXFieldGof);
      
      WXFieldPar<REAL> *pWXFieldRwp=new WXFieldPar<REAL>(this,"Rwp",-1,&mRwp,70);
      pStats->Add(pWXFieldRwp    ,0,wxALIGN_CENTER);
      mList.Add(pWXFieldRwp);
      
      WXFieldPar<REAL> *pWXFieldRp=new WXFieldPar<REAL>(this,"Rp",-1,&mRp,70);
      pStats->Add(pWXFieldRp    ,0,wxALIGN_CENTER);
      mList.Add(pWXFieldRp);
      
      mpSizer->Add(pStats);
         
   this->BottomLayout(0);
   this->CrystUpdate(true);
   VFN_DEBUG_MESSAGE("WXDiffractionSingleCrystal::WXDiffractionSingleCrystal():End",6)
}

void WXDiffractionSingleCrystal::CrystUpdate(const bool uui,const bool lock)
{
   VFN_DEBUG_ENTRY("WXDiffractionSingleCrystal::CrystUpdate()",6)
   if(lock) mMutex.Lock();
   WXCrystValidateAllUserInput();

   mChi2=mpData->GetChi2();
   if(0==mpData->GetIobs().numElements()) mGoF=0;
   else mGoF=mpData->GetChi2()/mpData->GetIobs().numElements();
   mRwp=mpData->GetRw();
   mRp=mpData->GetR();
   if(mpGraph!=0)
   {
      const CrystVector_REAL *mpCalc=&(mpData->GetIcalc());
      const CrystVector_REAL *mpObs =&(mpData->GetIobs());
      const CrystVector_REAL *mpSinThetaOverLambda=&(mpData->GetSinThetaOverLambda());
      const unsigned long nb=mpCalc->numElements();
      mX    .resize(nb);
      mIobs .resize(nb);
      mIcalc.resize(nb);
      for(unsigned long i=0;i<nb;i++)
      {
         mX[i]=(*mpSinThetaOverLambda)(i)*2;//1/d
         mIobs[i] =(*mpObs)(i);
         mIcalc[i]=(*mpCalc)(i);
      }
   }
   if(lock) mMutex.Unlock();
   this->WXRefinableObj::CrystUpdate(uui,lock);
   VFN_DEBUG_EXIT("WXDiffractionSingleCrystal::CrystUpdate()",6)
} 

void WXDiffractionSingleCrystal::NotifyDeleteGraph()
{
   mpGraph=0;
}

void WXDiffractionSingleCrystal::OnMenuSimulate(wxCommandEvent & WXUNUSED(event))
{
   WXCrystValidateAllUserInput();
   double theta;
   {
      wxTextEntryDialog dialog(this,_T("Theta Max"),
                              _T("Enter maximum Theta (degrees)"),_T("50"),wxOK | wxCANCEL);
      if(wxID_OK!=dialog.ShowModal())
      {
         VFN_DEBUG_EXIT("WXDiffractionSingleCrystal::OnMenuSimulate():Cancelled",6)
         return;
      }
      dialog.GetValue().ToDouble(&theta);
   }
   wxString choices[2];
   choices[0]=_T("all reflections (keep Friedel&Bijvoet mates)");
   choices[1]=_T("only unique reflections");
   wxSingleChoiceDialog dialog(this,_T("Choose method"),
                               _T("Choose method for Reflection generation"),
                               2,choices,NULL,wxOK | wxCANCEL);
   if(wxID_OK!=dialog.ShowModal()) return;
   const int choice=dialog.GetSelection();
   if(0==choice) mpData->GenHKLFullSpace(theta*DEG2RAD,false);
   else mpData->GenHKLFullSpace(theta*DEG2RAD,true);
}
void WXDiffractionSingleCrystal::OnMenuImport(wxCommandEvent & event)
{
   if(event.GetId()== ID_DIFFSINGLECRYST_MENU_IMPORT_HKLIOBS)
   {
      wxFileDialog open(this,_T("Choose file to import from"),
                                     _T(""),_T(""),_T("*.*"),wxFD_OPEN | wxFD_FILE_MUST_EXIST);
      if(open.ShowModal() != wxID_OK) return;
      long nb=0;
      {
         wxTextEntryDialog dialog(this,_T("Number of reflections"),
                                 _T("Enter The number of reflections to import"),_T("50"),
                                 wxOK | wxCANCEL);
         if(wxID_OK!=dialog.ShowModal())
         {
            VFN_DEBUG_EXIT("WXDiffractionSingleCrystal))OnMenuImport())Cancelled",6)
            return;
         }
         dialog.GetValue().ToLong(&nb);
      }
      mpData->ImportHklIobs(string(open.GetPath().ToAscii()),nb);
      mpData->UpdateDisplay();
      return;
   }
   if(event.GetId()== ID_DIFFSINGLECRYST_MENU_IMPORT_HKLIOBSSIGMA)
   {
      wxFileDialog open(this,_T("Choose file to import from"),
                                     _T(""),_T(""),_T("*.*"),wxFD_OPEN | wxFD_FILE_MUST_EXIST);
      if(open.ShowModal() != wxID_OK) return;
      long nb=0;
      {
         wxTextEntryDialog dialog(this,_T("Number of reflections"),
                                 _T("Enter The number of reflections to import"),_T("50"),
                                 wxOK | wxCANCEL);
         if(wxID_OK!=dialog.ShowModal())
         {
            VFN_DEBUG_EXIT("WXDiffractionSingleCrystal))OnMenuImport())Cancelled",6)
            return;
         }
         dialog.GetValue().ToLong(&nb);
      }
      mpData->ImportHklIobsSigma(string(open.GetPath().ToAscii()),nb);
      mpData->UpdateDisplay();
      return;
   }
   if(event.GetId()== ID_DIFFSINGLECRYST_MENU_IMPORT_SHELXHKLF4)
   {
      wxFileDialog open(this,_T("Choose Shelx file to import from"),
                                     _T(""),_T(""),_T("*.hkl"),wxFD_OPEN | wxFD_FILE_MUST_EXIST);
      if(open.ShowModal() != wxID_OK) return;
      mpData->ImportShelxHKLF4(string(open.GetPath().ToAscii()));
      mpData->UpdateDisplay();
      return;
   }
   if(event.GetId()== ID_DIFFSINGLECRYST_MENU_IMPORT_CIF)
   {
      wxFileDialog open(this,_T("Choose CIF file to import from"),
                                     _T(""),_T(""),_T("*.cif"),wxFD_OPEN | wxFD_FILE_MUST_EXIST);
      if(open.ShowModal() != wxID_OK) return;
      mpData->ImportCIF(string(open.GetPath().ToAscii()));
      mpData->UpdateDisplay();
      return;
   }
   if(event.GetId()== ID_DIFFSINGLECRYST_MENU_IMPORT_HKLIOBSGROUP)
   {
      wxFileDialog open(this,_T("Choose data file"),
                                     _T(""),_T(""),_T("*.*"),wxFD_OPEN | wxFD_FILE_MUST_EXIST);
      if(open.ShowModal() != wxID_OK) return;
      mpData->ImportHklIobsGroup(string(open.GetPath().ToAscii()));
      mpData->UpdateDisplay();
      return;
   }
   if(event.GetId()== ID_DIFFSINGLECRYST_MENU_IMPORT_JANAM91)
   {
      wxFileDialog open(this,_T("Choose data file"),
                                     _T(""),_T(""),_T("*.*"),wxFD_OPEN | wxFD_FILE_MUST_EXIST);
      if(open.ShowModal() != wxID_OK) return;
      mpData->ImportHklIobsSigmaJanaM91(string(open.GetPath().ToAscii()));
      mpData->UpdateDisplay();
      return;
   }
}
void WXDiffractionSingleCrystal::OnMenuSaveHKLIobsIcalc(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXDiffractionSingleCrystal::OnMenuSaveHKLIobsIcalc()",6)
   WXCrystValidateAllUserInput();
   wxFileDialog save(this,_T("Choose a file"),_T(""),_T(""),_T("*.txt"),wxFD_SAVE | wxFD_OVERWRITE_PROMPT);
   if(save.ShowModal() != wxID_OK) return;
   mpData->SaveHKLIobsIcalc(string(save.GetPath().ToAscii()));
}
void WXDiffractionSingleCrystal::OnMenuSaveHKLFcalc(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXDiffractionSingleCrystal::OnMenuSaveHKLFcalc()",6)
   WXCrystValidateAllUserInput();
   wxFileDialog save(this,_T("Choose a file"),_T(""),_T(""),_T("*.txt"),wxFD_SAVE | wxFD_OVERWRITE_PROMPT);
   if(save.ShowModal() != wxID_OK) return;
   ofstream os(save.GetPath().ToAscii());
   mpData->PrintFhklCalcDetail(os);
   os.close();
   mpData->GetCrystal().GetScatteringComponentList().Print();
}
void WXDiffractionSingleCrystal::OnMenuSetWavelength(wxCommandEvent &event)
{
   WXCrystValidateAllUserInput();
   //:TODO: Use wxRadiation instead
   if(event.GetId()== ID_DIFFSINGLECRYST_MENU_WAVELENGTH_XRAY)
      mpData->SetRadiationType(RAD_XRAY);
   if(event.GetId()== ID_DIFFSINGLECRYST_MENU_WAVELENGTH_NEUTRON)
      mpData->SetRadiationType(RAD_NEUTRON);
   if(event.GetId()== ID_DIFFSINGLECRYST_MENU_WAVELENGTH_ELECTRON)
      mpData->SetRadiationType(RAD_ELECTRON);
   if(event.GetId()== ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET)
   {
      double lambda;
      wxTextEntryDialog dialog(this,_T("new Wavelength)"),
                              _T("Enter new Wavelength (Angstroems)"),_T("1"),wxOK | wxCANCEL);
      if(wxID_OK!=dialog.ShowModal())
      {
         VFN_DEBUG_EXIT("WXDiffractionSingleCrystal))OnMenuSetWavelength())Monochromatic)Cancelled",6)
         return;
      }
      dialog.GetValue().ToDouble(&lambda);
      mpData->SetWavelength(lambda);
   }
   if(event.GetId()== ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_AG)
      mpData->SetWavelength("Ag");
   if(event.GetId()== ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_MO)
      mpData->SetWavelength("Mo");
   if(event.GetId()== ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_CU)
      mpData->SetWavelength("Cu");
   if(event.GetId()== ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_FE)
      mpData->SetWavelength("Fe");
   if(event.GetId()== ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_CO)
      mpData->SetWavelength("Co");
   if(event.GetId()== ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_CR)
      mpData->SetWavelength("Cr");
   if(event.GetId()== ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_AGA1)
      mpData->SetWavelength("AgA1");
   if(event.GetId()== ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_MOA1)
      mpData->SetWavelength("MoA1");
   if(event.GetId()== ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_CUA1)
      mpData->SetWavelength("CuA1");
   if(event.GetId()== ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_FEA1)
      mpData->SetWavelength("FeA1");
   if(event.GetId()== ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_COA1)
      mpData->SetWavelength("CoA1");
   if(event.GetId()== ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_CRA1)
      mpData->SetWavelength("CrA1");
   this->CrystUpdate(true,true);
}

void WXDiffractionSingleCrystal::OnMenuShowGraph(wxCommandEvent &event)
{
   VFN_DEBUG_MESSAGE("WXDiffractionSingleCrystal::OnMenuShowGraph()"<<mpGraph,6)
   if(mpGraph!=0) return;
   if(mpData->GetNbRefl()<=0) return;
   WXCrystValidateAllUserInput();
   std::string s=mpData->GetName();
   if(s.size()==0) s=mpData->GetCrystal().GetName();
   s="Single Crystal data:"+s;
   wxFrame *frame= new wxFrame(this,-1,wxString::FromAscii(s.c_str()),
                               wxDefaultPosition,wxSize(500,300));
   mpGraph = new WXDiffractionSingleCrystalGraph(frame,this);
   mpGraph->SetXLabel(_T("1/d (A)"));
   mpGraph->SetYLabel(_T("Intensity"));
   mGrapIdObs =mpGraph->AddGraph("Iobs");
   mGrapIdCalc=mpGraph->AddGraph("Icalc");
   
   wxSizer *ps=new wxBoxSizer(wxHORIZONTAL);
   ps->Add(mpGraph,1,wxEXPAND);
   frame->SetSizer(ps);
   frame->SetAutoLayout(true);
   
   //frame->CreateStatusBar(2);
   frame->Show(true);
   this->CrystUpdate(true);
}

void WXDiffractionSingleCrystal::OnChangeCrystal(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXDiffractionSingleCrystal::OnChangeCrystal()",6)
   WXCrystValidateAllUserInput();
   int choice;
   Crystal *cryst=dynamic_cast<Crystal*>
      (WXDialogChooseFromRegistry(gCrystalRegistry,(wxWindow*)this,
         "Choose a Crystal Structure:",choice));
   if(0==cryst) return;
   mpData->SetCrystal(*cryst);
   this->CrystUpdate(true,true);
}
void WXDiffractionSingleCrystal::UpdateUI(const bool lock)
{
   if(lock) mMutex.Lock();
   if(&(mpData->GetCrystal())!=0) mpFieldCrystal->SetValue(mpData->GetCrystal().GetName());
   if(mpGraph!=0)
   {
      mpGraph->SetGraphData(mGrapIdObs,mX,mIobs);
      mpGraph->SetGraphData(mGrapIdCalc,mX,mIcalc);
      mpGraph->UpdateDisplay();
   }
   this->WXRefinableObj::UpdateUI(false);
   if(lock) mMutex.Unlock();
}

}//namespace
