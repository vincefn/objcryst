#ifdef __GNUG__
    #pragma implementation "minimal.cpp"
    #pragma interface "minimal.cpp"
#endif

// For compilers that support precompilation, includes "wx/wx.h".
#include "wx/wxprec.h"

#ifdef __BORLANDC__
    #pragma hdrstop
#endif

// for all others, include the necessary headers (this file is usually all you
// need because it includes almost all "standard" wxWindows headers)
#ifndef WX_PRECOMP
    #include "wx/wx.h"
#endif

#include "wx/notebook.h"

#include <locale.h>

#include "ObjCryst/IO.h"
#include "wxCryst/wxCrystal.h"
#include "wxCryst/wxPowderPattern.h"
#include "wxCryst/wxGlobalOptimObj.h"

using namespace ObjCryst;
using namespace std;
// ----------------------------------------------------------------------------
// private classes
// ----------------------------------------------------------------------------

class MyApp : public wxApp
{
public:
    virtual bool OnInit();
};

// main frame
class WXCrystMainFrame : public wxFrame
{
public:
   WXCrystMainFrame(const wxString& title, const wxPoint& pos, const wxSize& size);
   void OnQuit(wxCommandEvent& WXUNUSED(event));
   void OnAbout(wxCommandEvent& WXUNUSED(event));
   void OnLoad(wxCommandEvent& event);
   void OnSave(wxCommandEvent& WXUNUSED(event));
   void OnAddCrystal(wxCommandEvent& WXUNUSED(event));
   void OnAddPowderPattern(wxCommandEvent& WXUNUSED(event));
   void OnAddGlobalOptimObj(wxCommandEvent& WXUNUSED(event));
   void OnSetDebugLevel0(wxCommandEvent& WXUNUSED(event));
   void OnSetDebugLevel1(wxCommandEvent& WXUNUSED(event));
   void OnSetDebugLevel2(wxCommandEvent& WXUNUSED(event));
   void OnSetDebugLevel3(wxCommandEvent& WXUNUSED(event));
   void OnSetDebugLevel4(wxCommandEvent& WXUNUSED(event));
   void OnSetDebugLevel5(wxCommandEvent& WXUNUSED(event));
   void OnSetDebugLevel6(wxCommandEvent& WXUNUSED(event));
   void OnSetDebugLevel7(wxCommandEvent& WXUNUSED(event));
   void OnSetDebugLevel8(wxCommandEvent& WXUNUSED(event));
   void OnSetDebugLevel9(wxCommandEvent& WXUNUSED(event));
   void OnSetDebugLevel10(wxCommandEvent& WXUNUSED(event));
private:
   wxScrolledWindow *mpWin1,*mpWin2,*mpWin3;
    DECLARE_EVENT_TABLE()
};
// ----------------------------------------------------------------------------
// constants
// ----------------------------------------------------------------------------
enum
{
   // menu items
   MENU_FILE_QUIT,
   MENU_HELP_ABOUT,
   MENU_FILE_LOAD,
   MENU_FILE_LOAD_OXY,
   MENU_FILE_SAVE,
   MENU_OBJECT_CREATE_CRYSTAL,
   MENU_OBJECT_CREATE_POWDERSPECTRUM,
   MENU_OBJECT_CREATE_GLOBALOPTOBJ,
   MENU_DEBUG_LEVEL0,
   MENU_DEBUG_LEVEL1,
   MENU_DEBUG_LEVEL2,
   MENU_DEBUG_LEVEL3,
   MENU_DEBUG_LEVEL4,
   MENU_DEBUG_LEVEL5,
   MENU_DEBUG_LEVEL6,
   MENU_DEBUG_LEVEL7,
   MENU_DEBUG_LEVEL8,
   MENU_DEBUG_LEVEL9,
   MENU_DEBUG_LEVEL10
};

// ----------------------------------------------------------------------------
// event tables and other macros for wxWindows
// ----------------------------------------------------------------------------

BEGIN_EVENT_TABLE(WXCrystMainFrame, wxFrame)
   EVT_MENU(MENU_FILE_QUIT,  WXCrystMainFrame::OnQuit)
   EVT_MENU(MENU_HELP_ABOUT, WXCrystMainFrame::OnAbout)
   EVT_MENU(MENU_FILE_LOAD, WXCrystMainFrame::OnLoad)
   EVT_MENU(MENU_FILE_LOAD_OXY, WXCrystMainFrame::OnLoad)
   EVT_MENU(MENU_FILE_SAVE, WXCrystMainFrame::OnSave)
   EVT_MENU(MENU_OBJECT_CREATE_CRYSTAL, WXCrystMainFrame::OnAddCrystal)
   EVT_MENU(MENU_OBJECT_CREATE_POWDERSPECTRUM, WXCrystMainFrame::OnAddPowderPattern)
   EVT_MENU(MENU_OBJECT_CREATE_GLOBALOPTOBJ, WXCrystMainFrame::OnAddGlobalOptimObj)
   EVT_MENU(MENU_DEBUG_LEVEL0, WXCrystMainFrame::OnSetDebugLevel0)
   EVT_MENU(MENU_DEBUG_LEVEL1, WXCrystMainFrame::OnSetDebugLevel1)
   EVT_MENU(MENU_DEBUG_LEVEL2, WXCrystMainFrame::OnSetDebugLevel2)
   EVT_MENU(MENU_DEBUG_LEVEL3, WXCrystMainFrame::OnSetDebugLevel3)
   EVT_MENU(MENU_DEBUG_LEVEL4, WXCrystMainFrame::OnSetDebugLevel4)
   EVT_MENU(MENU_DEBUG_LEVEL5, WXCrystMainFrame::OnSetDebugLevel5)
   EVT_MENU(MENU_DEBUG_LEVEL6, WXCrystMainFrame::OnSetDebugLevel6)
   EVT_MENU(MENU_DEBUG_LEVEL7, WXCrystMainFrame::OnSetDebugLevel7)
   EVT_MENU(MENU_DEBUG_LEVEL8, WXCrystMainFrame::OnSetDebugLevel8)
   EVT_MENU(MENU_DEBUG_LEVEL9, WXCrystMainFrame::OnSetDebugLevel9)
   EVT_MENU(MENU_DEBUG_LEVEL10,WXCrystMainFrame::OnSetDebugLevel10)
END_EVENT_TABLE()

IMPLEMENT_APP(MyApp)

// ============================================================================
// implementation
// ============================================================================

// 'Main program' equivalent: the program execution "starts" here
bool MyApp::OnInit()
{
	//set locale settings to standard
	setlocale(LC_NUMERIC,"C");

   WXCrystMainFrame *frame ;
   
   frame = new WXCrystMainFrame("FOX: Free Objects for Xtal structures v1.0.2 (dev)",
                                 wxPoint(50, 50), wxSize(450, 340));

   return TRUE;
}

// ----------------------------------------------------------------------------
// main frame
// ----------------------------------------------------------------------------

WXCrystMainFrame::WXCrystMainFrame(const wxString& title, const wxPoint& pos, const wxSize& size)
       : wxFrame((wxFrame *)NULL, -1, title, pos, size)
{
#ifdef __WXMAC__
   // we need this in order to allow the about menu relocation, since ABOUT is
   // not the default id of the about menu
   wxApp::s_macAboutMenuItemId = Minimal_About;
#endif

   // create a menu bar
      wxMenu *menuFile = new wxMenu;//
         menuFile->Append(MENU_FILE_LOAD, "Load", "Load some objects");
         //menuFile->Append(MENU_FILE_LOAD_OXY,"Load OLD .OXY","Load using the old .oxy format");
         menuFile->Append(MENU_FILE_SAVE, "Save", "Save Evrything...");
         menuFile->Append(MENU_FILE_QUIT, "E&xit\tAlt-Q", "Quit ");
      
      wxMenu *objectMenu = new wxMenu("", wxMENU_TEAROFF);
         objectMenu->Append(MENU_OBJECT_CREATE_CRYSTAL, "New Crystal",
                           "Add a new Crystal structure");
         objectMenu->Append(MENU_OBJECT_CREATE_POWDERSPECTRUM, "New PowderPattern",
                           "Add a new PowderPattern Object");
         objectMenu->Append(MENU_OBJECT_CREATE_GLOBALOPTOBJ, "New Global Optimization Object",
                           "Add a new Global Optimization Object");
      
      wxMenu *helpMenu = new wxMenu;
         helpMenu->Append(MENU_HELP_ABOUT, "&About...", "About ObjCryst...");

      wxMenuBar *menuBar = new wxMenuBar();
         menuBar->Append(menuFile,  "&File");
         menuBar->Append(objectMenu,"&Objects");
         #ifdef __DEBUG__
         wxMenu *debugMenu = new wxMenu;
            debugMenu->Append(MENU_DEBUG_LEVEL0, "Debug level 0 (lots of messages)");
            debugMenu->Append(MENU_DEBUG_LEVEL1, "Debug level 1");
            debugMenu->Append(MENU_DEBUG_LEVEL2, "Debug level 2");
            debugMenu->Append(MENU_DEBUG_LEVEL3, "Debug level 3");
            debugMenu->Append(MENU_DEBUG_LEVEL4, "Debug level 4");
            debugMenu->Append(MENU_DEBUG_LEVEL5, "Debug level 5");
            debugMenu->Append(MENU_DEBUG_LEVEL6, "Debug level 6");
            debugMenu->Append(MENU_DEBUG_LEVEL7, "Debug level 7");
            debugMenu->Append(MENU_DEBUG_LEVEL8, "Debug level 8");
            debugMenu->Append(MENU_DEBUG_LEVEL9, "Debug level 9");
            debugMenu->Append(MENU_DEBUG_LEVEL10,"Debug level 10 (few messages)");
         menuBar->Append(debugMenu,  "&Debug");
         #endif
         menuBar->Append(helpMenu,  "&Help");

   // ... and attach this menu bar to the frame
   SetMenuBar(menuBar);

#if wxUSE_STATUSBAR
   // create a status bar just for fun (by default with 1 pane only)
   CreateStatusBar(1);
   SetStatusText("Welcome to ObjCryst++!");
#endif // wxUSE_STATUSBAR

   
   /*   
   //crate something to display
      Crystal *cryst=new Crystal(8.482,5.398,6.959,"Pnma");
      cryst->SetName("PbSO4-(reference structure)");
      //Create 'real' PBSO4 structure (for reference only)
	      ScatteringPowerAtom *ScattPowPb=new ScatteringPowerAtom("Pb","Pb",1.48);
	      ScatteringPowerAtom *ScattPowS =new ScatteringPowerAtom("S" ,"S",0.74);
	      ScatteringPowerAtom *ScattPowO1=new ScatteringPowerAtom("O1","O",1.87);
	      ScatteringPowerAtom *ScattPowO2=new ScatteringPowerAtom("O2","O",1.76);
	      ScatteringPowerAtom *ScattPowO3=new ScatteringPowerAtom("O3","O",1.34);
         cryst->AddScatteringPower(ScattPowPb);
         cryst->AddScatteringPower(ScattPowS);
         cryst->AddScatteringPower(ScattPowO1);
         cryst->AddScatteringPower(ScattPowO2);
         cryst->AddScatteringPower(ScattPowO3);
	      ObjCryst::Atom *Pb=new ObjCryst::Atom(.188,.250,.167,"Pb",ScattPowPb   ,.5);
	      ObjCryst::Atom *S=new ObjCryst::Atom (.437,.750,.186,"S" ,ScattPowS    ,.5);
	      ObjCryst::Atom *O1=new ObjCryst::Atom(.595,.750,.100,"O1",ScattPowO1   ,.5);
	      ObjCryst::Atom *O2=new ObjCryst::Atom(.319,.750,.043,"O2",ScattPowO2   ,.5);
	      ObjCryst::Atom *O3=new ObjCryst::Atom(.415,.974,.306,"O3",ScattPowO3   ,1.);
         cryst->AddScatterer(Pb);
         cryst->AddScatterer(S);
         cryst->AddScatterer(O1);
         cryst->AddScatterer(O2);
         cryst->AddScatterer(O3);
      // Create a Global Optim object
         GlobalOptimObj *globalOptObj=new GlobalOptimObj;
         globalOptObj->SetName("PbSO4 optimization");
      // Create a few PowderPattern objects
         //1
            PowderPattern *data1=new PowderPattern;
            data1->SetWavelength("CuA1");
            data1->SetName("PbSO4-RR");
            data1->ImportPowderPatternFullprof("../example/pbso4-xray/roundrobin.dat");
            
            data1->SetSigmaToPoisson();
            data1->SetWeightToInvSigmaSq();
            //Components
               PowderPatternDiffraction * diffData=new PowderPatternDiffraction;
               diffData->SetCrystal(*cryst);
               data1->AddPowderPatternComponent(*diffData);
               diffData->SetName("PbSo4-diffraction");
               diffData->SetIsIgnoringImagScattFact(true);
               //approximate (hand-determined) background
               PowderPatternBackground *backgdData= new PowderPatternBackground;
               //backgdData->ImportUserBackground("pbso4-xray/background.dat");
               backgdData->SetName("PbSo4-background");
               backgdData->ImportUserBackground("../example/pbso4-xray/background.dat");
               data1->AddPowderPatternComponent(*backgdData);

               diffData->SetReflectionProfilePar(PROFILE_PSEUDO_VOIGT,
                                                 .03*DEG2RAD*DEG2RAD,
                                                 0*DEG2RAD*DEG2RAD,
                                                 0*DEG2RAD*DEG2RAD,
                                                 0.3,0);
         //2
            PowderPattern *data2=new PowderPattern;
            data2->SetRadiationType(RAD_NEUTRON);
            data2->SetName("My second, very stupd PowderDiff object !!");
   
   
   //Simple tests
      //WXFieldRefPar *mpFieldX    =new WXFieldRefPar(this,"x:",0,-1,&(Pb->GetPar(1)) );

      //WXFieldName *mpWXTitle = new WXFieldName(this,"name:",0,0,300);
      //mpWXTitle->SetValue("Toto est content");
      
      //ScattPowPb->WXCreate(this);
      //Pb->WXCreate(this);
      //cryst->GetScatteringPowerRegistry().WXCreate(this);
      //cryst->GetScattererRegistry().WXCreate(this);
      cryst->WXCreate(this);

   
   */
   // Create the notebook

      wxNotebook *notebook = new wxNotebook(this, -1);

      wxLayoutConstraints* c = new wxLayoutConstraints;
      c = new wxLayoutConstraints;
      c->left.SameAs(this, wxLeft, 2);
      c->right.SameAs(this, wxRight, 2);
      c->top.SameAs(this, wxTop, 2);
      c->bottom.SameAs(this, wxBottom, 2);

      notebook->SetConstraints(c);

   // First window -Crystals
      wxScrolledWindow *mpWin1 = new wxScrolledWindow(notebook, -1);
      mpWin1->SetScrollbars( 10, 10, 0, 500 );
      //wxBoxSizer * sizer1=new wxBoxSizer(wxVERTICAL);
      //sizer1->Add(gCrystalRegistry.WXCreate(mpWin1));
      //mpWin1->SetSizer(sizer1);
      //mpWin1->SetAutoLayout(true);
      gCrystalRegistry.WXCreate(mpWin1);
      notebook->AddPage(mpWin1, "Crystals", TRUE);

   // Second window - DiffractionData
      wxScrolledWindow *mpWin2 = new wxScrolledWindow(notebook, -1);
      mpWin2->SetScrollbars( 10, 10, 0, 500 );
      //wxBoxSizer * sizer2=new wxBoxSizer(wxVERTICAL);
      //sizer2->Add(gPowderPatternRegistry.WXCreate(mpWin2));
      //mpWin2->SetSizer(sizer2);
      //win2->SetAutoLayout(true);
      gPowderPatternRegistry.WXCreate(mpWin2);
      notebook->AddPage(mpWin2,"Powder Diffraction",true);
      
   // Third window - Global Optimization
      wxScrolledWindow *mpWin3 = new wxScrolledWindow(notebook, -1);
      mpWin3->SetScrollbars( 10, 10, 0, 500 );
      //wxBoxSizer * sizer3=new wxBoxSizer(wxVERTICAL);
      //sizer3->Add(gGlobalOptimObjRegistry.WXCreate(mpWin3));
      //mpWin3->SetSizer(sizer3);
      //mpWin3->SetAutoLayout(true);
      gGlobalOptimObjRegistry.WXCreate(mpWin3);
      notebook->AddPage(mpWin3,"Global Optimization",true);
   this->Show(TRUE);
   this->Layout();
   //Splash Screen
      wxCommandEvent event;
      this->OnAbout(event);
}

void WXCrystMainFrame::OnQuit(wxCommandEvent& WXUNUSED(event))
{
   // TRUE is to force the frame to close
   Close(TRUE);
}

void WXCrystMainFrame::OnAbout(wxCommandEvent& WXUNUSED(event))
{
   wxString msg;
   msg.Printf( _T("F.O.X. - Free Objects for Xtal structures\n")
               _T("Version 1.0.2 - (dev)\n\n")
               _T("(c) 2000-2001 Vincent FAVRE-NICOLIN, vincefn@users.sourceforge.net\n")
               _T("            & Radovan CERNY, University of Geneva\n\n")
               _T("http://objcryst.sourceforge.net\n")
               _T("http://www.ccp14.ac.uk/ccp/web-mirrors/objcryst/ (Mirror)\n\n")
               _T("Note: this is still in development. Please look \n")
               _T("      for newer version on the above web sites,\n")
               _T("      or subscribe, to the mailing-list at :\n")
               _T("      http://lists.sourceforge.net/lists/listinfo/objcryst-foxx\n\n")
               _T("Project supported by the Swiss National Science Foundation (#21-53847.98)")
              );

   wxMessageBox(msg, "About Fox", wxOK | wxICON_INFORMATION, this);
}
void WXCrystMainFrame::OnLoad(wxCommandEvent& event)
{
   wxFileDialog *open;
   switch(event.GetId())
   {
      case MENU_FILE_LOAD_OXY:
      {
         open= new wxFileDialog(this,"Choose File :",
                                              "","","*.oxy",wxOPEN | wxFILE_MUST_EXIST);
         if(open->ShowModal() != wxID_OK) return;
         IOCrystFileLoadAllObject(open->GetPath().c_str());
         break;
      }
      case MENU_FILE_LOAD:
      {
         open= new wxFileDialog(this,"Choose File :",
                                              "","","*.xml",wxOPEN | wxFILE_MUST_EXIST);
         if(open->ShowModal() != wxID_OK) return;
         XMLCrystFileLoadAllObject(open->GetPath().c_str());
         break;
      }
   }
   open->Destroy();
}
void WXCrystMainFrame::OnSave(wxCommandEvent& WXUNUSED(event))
{
   wxFileDialog *open= new wxFileDialog(this,"Choose File to save all objects:",
                                        "","","*.xml", wxSAVE | wxOVERWRITE_PROMPT);
   if(open->ShowModal() != wxID_OK) return;
   XMLCrystFileSaveGlobal(open->GetPath().c_str());
   open->Destroy();
}
void WXCrystMainFrame::OnAddCrystal(wxCommandEvent& WXUNUSED(event))
{
   Crystal* obj;
   obj=new Crystal;
   obj->SetName("Change Me!");
}
void WXCrystMainFrame::OnAddPowderPattern(wxCommandEvent& WXUNUSED(event))
{
   PowderPattern* obj;
   obj=new PowderPattern;
   obj->SetName("Change Me!");
}
void WXCrystMainFrame::OnAddGlobalOptimObj(wxCommandEvent& WXUNUSED(event))
{
   GlobalOptimObj* obj;
   obj=new GlobalOptimObj;
   obj->SetName("Change Me!");
}
void WXCrystMainFrame::OnSetDebugLevel0(wxCommandEvent& WXUNUSED(event))
{
   VFN_DEBUG_GLOBAL_LEVEL(0);
}
void WXCrystMainFrame::OnSetDebugLevel1(wxCommandEvent& WXUNUSED(event))
{
   VFN_DEBUG_GLOBAL_LEVEL(1);
}
void WXCrystMainFrame::OnSetDebugLevel2(wxCommandEvent& WXUNUSED(event))
{
   VFN_DEBUG_GLOBAL_LEVEL(2);
}
void WXCrystMainFrame::OnSetDebugLevel3(wxCommandEvent& WXUNUSED(event))
{
   VFN_DEBUG_GLOBAL_LEVEL(3);
}
void WXCrystMainFrame::OnSetDebugLevel4(wxCommandEvent& WXUNUSED(event))
{
   VFN_DEBUG_GLOBAL_LEVEL(4);
}
void WXCrystMainFrame::OnSetDebugLevel5(wxCommandEvent& WXUNUSED(event))
{
   VFN_DEBUG_GLOBAL_LEVEL(5);
}
void WXCrystMainFrame::OnSetDebugLevel6(wxCommandEvent& WXUNUSED(event))
{
   VFN_DEBUG_GLOBAL_LEVEL(6);
}
void WXCrystMainFrame::OnSetDebugLevel7(wxCommandEvent& WXUNUSED(event))
{
   VFN_DEBUG_GLOBAL_LEVEL(7);
}
void WXCrystMainFrame::OnSetDebugLevel8(wxCommandEvent& WXUNUSED(event))
{
   VFN_DEBUG_GLOBAL_LEVEL(8);
}
void WXCrystMainFrame::OnSetDebugLevel9(wxCommandEvent& WXUNUSED(event))
{
   VFN_DEBUG_GLOBAL_LEVEL(9);
}
void WXCrystMainFrame::OnSetDebugLevel10(wxCommandEvent& WXUNUSED(event))
{
   VFN_DEBUG_GLOBAL_LEVEL(10);
}
