//#include <sstream> //for stringstream
#include <fstream>

#include "wx/wx.h"
#include "wxCryst/wxDiffractionSingleCrystal.h"
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
//    WXDiffractionSingleCrystal
//
////////////////////////////////////////////////////////////////////////
BEGIN_EVENT_TABLE(WXDiffractionSingleCrystal, wxWindow)
   EVT_BUTTON(ID_WXOBJ_COLLAPSE,                        WXCrystObj::OnToggleCollapse)
   EVT_MENU(ID_REFOBJ_MENU_OBJ_SAVE,                    WXRefinableObj::OnMenuSave)
   EVT_MENU(ID_REFOBJ_MENU_OBJ_LOAD,                    WXRefinableObj::OnMenuLoad)
   EVT_MENU(ID_REFOBJ_MENU_PAR_FIXALL,                  WXRefinableObj::OnMenuFixAllPar)
   EVT_MENU(ID_REFOBJ_MENU_PAR_UNFIXALL,                WXRefinableObj::OnMenuUnFixAllPar)
   EVT_MENU(ID_DIFFSINGLECRYST_MENU_SAVEHKLIOBSICALC,   WXDiffractionSingleCrystal::OnMenuSaveHKLIobsIcalc)
   EVT_MENU(ID_DIFFSINGLECRYST_MENU_SIMULATE,           WXDiffractionSingleCrystal::OnMenuSimulate)
   EVT_MENU(ID_DIFFSINGLECRYST_MENU_IMPORT_HKLIOBS,     WXDiffractionSingleCrystal::OnMenuImport)
   EVT_MENU(ID_DIFFSINGLECRYST_MENU_IMPORT_HKLIOBSSIGMA,WXDiffractionSingleCrystal::OnMenuImport)
   EVT_MENU(ID_DIFFSINGLECRYST_MENU_IMPORT_JANAM91,     WXDiffractionSingleCrystal::OnMenuImport)
   EVT_BUTTON(ID_DIFFSINGLECRYST_CRYSTAL,               WXDiffractionSingleCrystal::OnChangeCrystal)
   //EVT_MENU(ID_DIFFSINGLECRYST_MENU_FITSCALE_R,         WXDiffractionSingleCrystal::OnMenuFitScaleForR)
   //EVT_MENU(ID_DIFFSINGLECRYST_MENU_FITSCALE_RW,        WXDiffractionSingleCrystal::OnMenuFitScaleForRw)
   EVT_MENU(ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET,     WXDiffractionSingleCrystal::OnMenuSetWavelength)
   EVT_MENU(ID_DIFFSINGLECRYST_MENU_WAVELENGTH_XRAY,    WXDiffractionSingleCrystal::OnMenuSetWavelength)
   EVT_MENU(ID_DIFFSINGLECRYST_MENU_WAVELENGTH_NEUTRON, WXDiffractionSingleCrystal::OnMenuSetWavelength)
   EVT_MENU(ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_AG,  WXDiffractionSingleCrystal::OnMenuSetWavelength)
   EVT_MENU(ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_MO,  WXDiffractionSingleCrystal::OnMenuSetWavelength)
   EVT_MENU(ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_CU,  WXDiffractionSingleCrystal::OnMenuSetWavelength)
   EVT_MENU(ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_FE,  WXDiffractionSingleCrystal::OnMenuSetWavelength)
   EVT_MENU(ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_CR,  WXDiffractionSingleCrystal::OnMenuSetWavelength)
   EVT_MENU(ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_AGA1,WXDiffractionSingleCrystal::OnMenuSetWavelength)
   EVT_MENU(ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_MOA1,WXDiffractionSingleCrystal::OnMenuSetWavelength)
   EVT_MENU(ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_CUA1,WXDiffractionSingleCrystal::OnMenuSetWavelength)
   EVT_MENU(ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_FEA1,WXDiffractionSingleCrystal::OnMenuSetWavelength)
   EVT_MENU(ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_CRA1,WXDiffractionSingleCrystal::OnMenuSetWavelength)
   EVT_UPDATE_UI(ID_CRYST_UPDATEUI, 						  WXDiffractionSingleCrystal::OnUpdateUI)
END_EVENT_TABLE()

WXDiffractionSingleCrystal::WXDiffractionSingleCrystal(wxWindow *parent,
																		 DiffractionDataSingleCrystal* data):
WXRefinableObj(parent,data),mpData(data)
{
   VFN_DEBUG_MESSAGE("WXDiffractionSingleCrystal::WXDiffractionSingleCrystal()",6)
   mpWXTitle->SetForegroundColour(wxColour(255,0,0));
   // Menu
      mpMenuBar->AddMenu("File",ID_REFOBJ_MENU_OBJ);
         //:TODO: reactivate & test those menus
         //mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_REFOBJ_MENU_OBJ_SAVE,"Save");
         //mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_REFOBJ_MENU_OBJ_LOAD,"Load");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_DIFFSINGLECRYST_MENU_SAVEHKLIOBSICALC,
                                "Save HKL Iobs Icalc (text)");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_DIFFSINGLECRYST_MENU_SIMULATE,
                                "Simulation mode (generate HKL list)");
         //mpMenuBar->AppendSeparator();
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_DIFFSINGLECRYST_MENU_IMPORT_HKLIOBS,
                                 "Import HKL Iobs");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_DIFFSINGLECRYST_MENU_IMPORT_HKLIOBSSIGMA,
                                 "Import HKL Iobs Sigma");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_DIFFSINGLECRYST_MENU_IMPORT_JANAM91,
                                 "Import Jana M91");
        
      mpMenuBar->AddMenu("Radiation",ID_DIFFSINGLECRYST_MENU_WAVELENGTH);
         mpMenuBar->AddMenuItem(ID_DIFFSINGLECRYST_MENU_WAVELENGTH,
                                ID_DIFFSINGLECRYST_MENU_WAVELENGTH_NEUTRON,
                                "Neutron");
         mpMenuBar->AddMenuItem(ID_DIFFSINGLECRYST_MENU_WAVELENGTH,
                                ID_DIFFSINGLECRYST_MENU_WAVELENGTH_XRAY,
                                "X-Rays");
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
   //Radiation
      //:TODO: when all objects will use the same radiation object
      //mpSizer->Add(mpPowderPattern->mRadiation.WXCreate(this),0);
      //mList.Add(mpPowderPattern->mRadiation.WXGet());
   // Crystal
      mpFieldCrystal=new WXFieldChoice(this,ID_DIFFSINGLECRYST_CRYSTAL,"Crystal:",300);
      mpSizer->Add(mpFieldCrystal,0,wxALIGN_LEFT);
      mList.Add(mpFieldCrystal);
   
   this->CrystUpdate();
   this->Layout();
   VFN_DEBUG_MESSAGE("WXDiffractionSingleCrystal::WXDiffractionSingleCrystal():End",6)
}

void WXDiffractionSingleCrystal::OnMenuSimulate(wxCommandEvent & WXUNUSED(event))
{
	WXCrystValidateAllUserInput();
	double theta;
   {
      wxTextEntryDialog dialog(this,"Theta Max",
                              "Enter maximum Theta (degrees)","50",wxOK | wxCANCEL);
      if(wxID_OK!=dialog.ShowModal())
      {
         VFN_DEBUG_EXIT("WXDiffractionSingleCrystal::OnMenuSimulate():Cancelled",6)
         return;
      }
      dialog.GetValue().ToDouble(&theta);
   }
   wxString choices[2];
   choices[0]="all reflections (keep Friedel&Bijvoet mates)";
   choices[1]="only unique reflections";
   wxSingleChoiceDialog dialog(this,"Choose method",
										 "Choose method for Reflection generation",
				                   2,choices,0,wxOK | wxCANCEL);
   if(wxID_OK!=dialog.ShowModal()) return;
   const int choice=dialog.GetSelection();
	if(0==choice) mpData->GenHKLFullSpace(theta*DEG2RAD,false);
	else mpData->GenHKLFullSpace(theta*DEG2RAD,true);
	
	CrystVector_REAL iobs(mpData->GetNbRefl());
	iobs=100.;
	mpData->SetIobs(iobs);
}
void WXDiffractionSingleCrystal::OnMenuImport(wxCommandEvent & event)
{
	switch(event.GetId())
	{
		case ID_DIFFSINGLECRYST_MENU_IMPORT_HKLIOBS:
		{
   		wxFileDialog open(this,"Choose file to import from",
                                        "","","*.*",wxOPEN | wxFILE_MUST_EXIST);
   		if(open.ShowModal() != wxID_OK) return;
			long nb=0;
   		{
      		wxTextEntryDialog dialog(this,"Number of reflections",
                              		"Enter The number of reflections to import","50",
												wxOK | wxCANCEL);
      		if(wxID_OK!=dialog.ShowModal())
      		{
         		VFN_DEBUG_EXIT("WXDiffractionSingleCrystal::OnMenuImport():Cancelled",6)
         		return;
      		}
      		dialog.GetValue().ToLong(&nb);
   		}
   		mpData->ImportHklIobs(open.GetPath().c_str(),nb);
			return;
		}
		case ID_DIFFSINGLECRYST_MENU_IMPORT_HKLIOBSSIGMA:
		{
   		wxFileDialog open(this,"Choose file to import from",
                                        "","","*.*",wxOPEN | wxFILE_MUST_EXIST);
   		if(open.ShowModal() != wxID_OK) return;
			long nb=0;
   		{
      		wxTextEntryDialog dialog(this,"Number of reflections",
                              		"Enter The number of reflections to import","50",
												wxOK | wxCANCEL);
      		if(wxID_OK!=dialog.ShowModal())
      		{
         		VFN_DEBUG_EXIT("WXDiffractionSingleCrystal::OnMenuImport():Cancelled",6)
         		return;
      		}
      		dialog.GetValue().ToLong(&nb);
   		}
   		mpData->ImportHklIobsSigma(open.GetPath().c_str(),nb);
			return;
		}
		case ID_DIFFSINGLECRYST_MENU_IMPORT_JANAM91:
		{
		}
	}
}
void WXDiffractionSingleCrystal::OnMenuSaveHKLIobsIcalc(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXDiffractionSingleCrystal::OnMenuSaveHKLIobsIcalc()",6)
	WXCrystValidateAllUserInput();
   wxFileDialog save(this,"Choose a file","","","*.txt",wxSAVE | wxOVERWRITE_PROMPT);
   if(save.ShowModal() != wxID_OK) return;
   mpData->SaveHKLIobsIcalc(save.GetPath().c_str());
}
void WXDiffractionSingleCrystal::OnMenuSetWavelength(wxCommandEvent &event)
{
	WXCrystValidateAllUserInput();
	//:TODO: Use wxRadiation instead
   switch(event.GetId())
   {
      case ID_DIFFSINGLECRYST_MENU_WAVELENGTH_XRAY:
         mpData->SetRadiationType(RAD_XRAY);break;
      case ID_DIFFSINGLECRYST_MENU_WAVELENGTH_NEUTRON:
         mpData->SetRadiationType(RAD_NEUTRON);break;
      case ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET:
      {
         double lambda;
         wxTextEntryDialog dialog(this,"new Wavelength:",
                                 "Enter new Wavelength (Angstroems)","1",wxOK | wxCANCEL);
         if(wxID_OK!=dialog.ShowModal())
         {
            VFN_DEBUG_EXIT("WXDiffractionSingleCrystal::OnMenuSetWavelength():Monochromatic:Cancelled",6)
            return;
         }
         dialog.GetValue().ToDouble(&lambda);
         mpData->SetWavelength(lambda);break;
      }
      case ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_AG:
         mpData->SetWavelength("Ag");break;
      case ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_MO:
         mpData->SetWavelength("Mo");break;
      case ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_CU:
         mpData->SetWavelength("Cu");break;
      case ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_FE:
         mpData->SetWavelength("Fe");break;
      case ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_CR:
         mpData->SetWavelength("Cr");break;
      case ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_AGA1:
         mpData->SetWavelength("AgA1");break;
      case ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_MOA1:
         mpData->SetWavelength("MoA1");break;
      case ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_CUA1:
         mpData->SetWavelength("CuA1");break;
      case ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_FEA1:
         mpData->SetWavelength("FeA1");break;
      case ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_CRA1:
         mpData->SetWavelength("CrA1");break;
   }
   this->CrystUpdate();
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
   this->CrystUpdate();
}
void WXDiffractionSingleCrystal::OnUpdateUI(wxUpdateUIEvent& event)
{
	if(&(mpData->GetCrystal())!=0) mpFieldCrystal->SetValue(mpData->GetCrystal().GetName());
	this->WXRefinableObj::OnUpdateUI(event);
}

}//namespace
