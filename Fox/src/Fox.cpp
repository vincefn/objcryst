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

// For compilers that support precompilation, includes "wx/wx.h".
#ifndef __DARWIN__ // work around MacOSX type_info bug (??)
   #include "wx/wxprec.h"
#endif

#ifdef __BORLANDC__
    #pragma hdrstop
#endif

// for all others, include the necessary headers (this file is usually all you
// need because it includes almost all "standard" wxWindows headers)
#ifndef WX_PRECOMP
    #include "wx/wx.h"
#endif

#include "wx/tooltip.h"
#include "wx/notebook.h"

#include <locale.h>
#include <sstream>
#include <list>

#include "ObjCryst/IO.h"
#include "ObjCryst/Crystal.h"
#include "ObjCryst/PowderPattern.h"
#include "ObjCryst/DiffractionDataSingleCrystal.h"
#include "ObjCryst/Polyhedron.h"
#include "ObjCryst/test.h"
#include "RefinableObj/GlobalOptimObj.h"
#include "Quirks/VFNStreamFormat.h"
#include "wxCryst/wxCrystal.h"

#if defined(__WXGTK__) || defined(__WXMOTIF__) || defined(__WXMAC__) || defined(__WXMGL__) || defined(__WXX11__)
   #include "Fox.xpm"
#endif


using namespace ObjCryst;
using namespace std;
// ----------------------------------------------------------------------------
// private classes
// ----------------------------------------------------------------------------

class MyApp : public wxApp
{
public:
    virtual bool OnInit();
    virtual int OnExit();
};

// WXCrystScr
class WXCrystScrolledWindow:public wxScrolledWindow
{
   public:
      WXCrystScrolledWindow(wxWindow* parent);
      virtual bool Layout();
      void SetChild(const wxWindow* pChild);
   private:
      const wxWindow* mpChild;
      int mHeight,mWidth;
};

// main frame
class WXCrystMainFrame : public wxFrame
{
public:
   WXCrystMainFrame(const wxString& title, const wxPoint& pos, const wxSize& size,
                    const bool splashscreen=true);
   void OnQuit(wxCommandEvent& WXUNUSED(event));
   void OnAbout(wxCommandEvent& WXUNUSED(event));
   void OnLoad(wxCommandEvent& event);
   void OnSave(wxCommandEvent& WXUNUSED(event));
   void OnAddCrystal(wxCommandEvent& WXUNUSED(event));
   void OnAddPowderPattern(wxCommandEvent& WXUNUSED(event));
   void OnAddSingleCrystalData(wxCommandEvent& WXUNUSED(event));
   void OnAddGlobalOptimObj(wxCommandEvent& WXUNUSED(event));
   void OnAddGeneticAlgorithm(wxCommandEvent& WXUNUSED(event));
   void OnDebugTest(wxCommandEvent& event);
   void OnSetDebugLevel(wxCommandEvent& event);
   void OnUpdateUI(wxUpdateUIEvent& event);
   void OnToggleTooltips(wxCommandEvent& event);
private:
    DECLARE_EVENT_TABLE()
};
// ----------------------------------------------------------------------------
// For messaging the user
// ----------------------------------------------------------------------------
wxFrame *pMainFrameForUserMessage;

void WXCrystInformUserStdOut(const string &str)
{
   pMainFrameForUserMessage->SetStatusText((wxString)str.c_str());
}

// ----------------------------------------------------------------------------
// Speed test
// ----------------------------------------------------------------------------
void standardSpeedTest();

// ----------------------------------------------------------------------------
// constants
// ----------------------------------------------------------------------------
static const long MENU_FILE_QUIT=                      WXCRYST_ID();
static const long MENU_HELP_ABOUT=                     WXCRYST_ID();
static const long MENU_HELP_TOGGLETOOLTIP=             WXCRYST_ID();
static const long MENU_FILE_LOAD=                      WXCRYST_ID();
static const long MENU_FILE_LOAD_OXY=                  WXCRYST_ID();
static const long MENU_FILE_SAVE=                      WXCRYST_ID();
static const long MENU_OBJECT_CREATE_CRYSTAL=          WXCRYST_ID();
static const long MENU_OBJECT_CREATE_POWDERSPECTRUM=   WXCRYST_ID();
static const long MENU_OBJECT_CREATE_SINGLECRYSTALDATA=WXCRYST_ID();
static const long MENU_OBJECT_CREATE_GLOBALOPTOBJ=     WXCRYST_ID();
static const long MENU_OBJECT_CREATE_GENETICALGORITHM= WXCRYST_ID();
static const long MENU_DEBUG_LEVEL0=                   WXCRYST_ID();
static const long MENU_DEBUG_LEVEL1=                   WXCRYST_ID();
static const long MENU_DEBUG_LEVEL2=                   WXCRYST_ID();
static const long MENU_DEBUG_LEVEL3=                   WXCRYST_ID();
static const long MENU_DEBUG_LEVEL4=                   WXCRYST_ID();
static const long MENU_DEBUG_LEVEL5=                   WXCRYST_ID();
static const long MENU_DEBUG_LEVEL6=                   WXCRYST_ID();
static const long MENU_DEBUG_LEVEL7=                   WXCRYST_ID();
static const long MENU_DEBUG_LEVEL8=                   WXCRYST_ID();
static const long MENU_DEBUG_LEVEL9=                   WXCRYST_ID();
static const long MENU_DEBUG_LEVEL10=                  WXCRYST_ID();
static const long MENU_DEBUG_TEST1=                    WXCRYST_ID();
static const long MENU_DEBUG_TEST2=                    WXCRYST_ID();
static const long MENU_DEBUG_TEST3=                    WXCRYST_ID();

// ----------------------------------------------------------------------------
// event tables and other macros for wxWindows
// ----------------------------------------------------------------------------

BEGIN_EVENT_TABLE(WXCrystMainFrame, wxFrame)
   EVT_MENU(MENU_FILE_QUIT,  WXCrystMainFrame::OnQuit)
   EVT_MENU(MENU_HELP_ABOUT, WXCrystMainFrame::OnAbout)
   EVT_MENU(MENU_HELP_TOGGLETOOLTIP, WXCrystMainFrame::OnToggleTooltips)
   EVT_MENU(MENU_FILE_LOAD, WXCrystMainFrame::OnLoad)
   EVT_MENU(MENU_FILE_LOAD_OXY, WXCrystMainFrame::OnLoad)
   EVT_MENU(MENU_FILE_SAVE, WXCrystMainFrame::OnSave)
   EVT_MENU(MENU_OBJECT_CREATE_CRYSTAL, WXCrystMainFrame::OnAddCrystal)
   EVT_MENU(MENU_OBJECT_CREATE_POWDERSPECTRUM, WXCrystMainFrame::OnAddPowderPattern)
   EVT_MENU(MENU_OBJECT_CREATE_SINGLECRYSTALDATA, WXCrystMainFrame::OnAddSingleCrystalData)
   EVT_MENU(MENU_OBJECT_CREATE_GLOBALOPTOBJ, WXCrystMainFrame::OnAddGlobalOptimObj)
   EVT_MENU(MENU_OBJECT_CREATE_GENETICALGORITHM, WXCrystMainFrame::OnAddGeneticAlgorithm)
   EVT_MENU(MENU_DEBUG_LEVEL0, WXCrystMainFrame::OnSetDebugLevel)
   EVT_MENU(MENU_DEBUG_LEVEL1, WXCrystMainFrame::OnSetDebugLevel)
   EVT_MENU(MENU_DEBUG_LEVEL2, WXCrystMainFrame::OnSetDebugLevel)
   EVT_MENU(MENU_DEBUG_LEVEL3, WXCrystMainFrame::OnSetDebugLevel)
   EVT_MENU(MENU_DEBUG_LEVEL4, WXCrystMainFrame::OnSetDebugLevel)
   EVT_MENU(MENU_DEBUG_LEVEL5, WXCrystMainFrame::OnSetDebugLevel)
   EVT_MENU(MENU_DEBUG_LEVEL6, WXCrystMainFrame::OnSetDebugLevel)
   EVT_MENU(MENU_DEBUG_LEVEL7, WXCrystMainFrame::OnSetDebugLevel)
   EVT_MENU(MENU_DEBUG_LEVEL8, WXCrystMainFrame::OnSetDebugLevel)
   EVT_MENU(MENU_DEBUG_LEVEL9, WXCrystMainFrame::OnSetDebugLevel)
   EVT_MENU(MENU_DEBUG_LEVEL10,WXCrystMainFrame::OnSetDebugLevel)
   EVT_MENU(MENU_DEBUG_TEST1,  WXCrystMainFrame::OnDebugTest)
   EVT_MENU(MENU_DEBUG_TEST2,  WXCrystMainFrame::OnDebugTest)
   EVT_MENU(MENU_DEBUG_TEST3,  WXCrystMainFrame::OnDebugTest)
   EVT_UPDATE_UI(ID_CRYST_UPDATEUI, WXCrystMainFrame::OnUpdateUI)
END_EVENT_TABLE()

IMPLEMENT_APP(MyApp)

// ============================================================================
// implementation
// ============================================================================

// 'Main program' equivalent: the program execution "starts" here
bool MyApp::OnInit()
{
   TAU_PROFILE_SET_NODE(0); // sequential code 
   TAU_PROFILE("main()","int()",TAU_DEFAULT);
   //set locale settings to standard
   setlocale(LC_NUMERIC,"C");
   
   bool useGUI(true);
   long nbTrial(1000000);
   REAL finalCost=0.;
   bool silent=false;
   string outfilename("Fox-out.xml");
   bool randomize(false);
   bool only3D(false);
   bool loadFourier(false);
   list<string> vFourierFilename;
   for(int i=1;i<this->argc;i++)
   {
      if('-'==this->argv[i][0])
      {
         if(string("--nogui")==string(this->argv[i]))
         {
            useGUI=false;
            cout << "Running Fox without GUI"<<endl;
            continue;  
         }
         if(string("--randomize")==string(this->argv[i]))
         {
            randomize=true;
            cout << "Randomizing parameters before running"<<endl;
            continue;  
         }
         if(string("--silent")==string(this->argv[i]))
         {
            silent=true;
            cout << "Running Fox quietly"<<endl;
            continue;  
         }
         if(string("--finalcost")==string(this->argv[i]))
         {
            ++i;
            stringstream sstr(this->argv[i]);
            sstr >> finalCost;
            cout << "Fox will stop after reaching cost:"<<finalCost<<endl;
            continue;  
         }
         if(string("-n")==string(this->argv[i]))
         {
            ++i;
            stringstream sstr(this->argv[i]);
            sstr >> nbTrial;
            cout << "Fox will run for "<<nbTrial<<" trials"<<endl;
            continue;
         }
         if(string("-i")==string(this->argv[i]))
         {
            ++i;
            XMLCrystFileLoadAllObject(this->argv[i]);
            continue;
         }
         if(string("-o")==string(this->argv[i]))
         {
            ++i;
            outfilename=string(this->argv[i]);
            continue;
         }
         if(string("--loadfouriergrd")==string(this->argv[i]))
         {
            ++i;
            loadFourier=true;
            vFourierFilename.push_back(string(this->argv[i]));
            continue;
         }
         if(string("--only3d")==string(this->argv[i]))
         {
            only3D=true;
            continue;
         }
         if(string("--speedtest")==string(this->argv[i]))
         {
            standardSpeedTest();
            exit(0);
         }
         #ifdef __DEBUG__
         if(string("--debuglevel")==string(this->argv[i]))
         {
            int level;
            ++i;
            stringstream sstr(this->argv[i]);
            sstr >> level;
            VFN_DEBUG_GLOBAL_LEVEL(level);
            continue;
         }
         #endif
         cout <<"command-line arguments:"<<endl
              <<"   -i input.xml: input 'in.xml' file"<<endl
              <<"   --loadfouriergrd map.grd: load and display 'map.grd' fourier map with (first) crystal structure"<<endl
              <<"   --nogui: run without GUI, automatically launches optimization"<<endl
              <<"      options with --nogui:"<<endl
              <<"         -n 10000     : run for 10000 trials at most (default: 1000000)"<<endl
              <<"         -o out.xml   : output in 'out.xml'"<<endl
              <<"         --randomize  : randomize initial configuration"<<endl
              <<"         --silent     : (almost) no text output"<<endl
              <<"         --finalcost 0.15 : run optimization until cost < 0.15"<<endl
              <<"   --speedtest: run speed test (takes about 1mn30s)"<<endl
              <<endl;
         exit(0);  
      }
   }
   
   
   if(randomize)
      for(int i=0;i<gOptimizationObjRegistry.GetNb();i++)
         gOptimizationObjRegistry.GetObj(i).RandomizeStartingConfig();
   
   if(!useGUI)
   {
      if(nbTrial!=0)
         for(int i=0;i<gOptimizationObjRegistry.GetNb();i++)
            for(int j=0;j<5;j++)
               gOptimizationObjRegistry.GetObj(i).Optimize(nbTrial,silent,finalCost);
      XMLCrystFileSaveGlobal(outfilename);
      cout <<"End of Fox execution. Bye !"<<endl;
      exit (0);
   }
   
   WXCrystMainFrame *frame ;
   
   frame = new WXCrystMainFrame("FOX: Free Objects for Xtal structures v1.6",
                                 wxPoint(50, 50), wxSize(550, 400),!loadFourier);
   // Use the main frame status bar to pass messages to the user
      pMainFrameForUserMessage=frame;
      fpObjCrystInformUser=&WXCrystInformUserStdOut;
      
   if(loadFourier)
   {
      //wxFrame *pWXFrame= new wxFrame((wxFrame *)NULL, -1, "FOX", wxPoint(50, 50), wxSize(550, 400));
      //wxScrolledWindow *pWXScWin=new wxScrolledWindow(pWXFrame,-1);
      //WXCrystal * pWXCrystal=new WXCrystal(pWXScWin,&(gCrystalRegistry.GetObj(0)));
      WXCrystal *pWXCryst=dynamic_cast<WXCrystal*> (gCrystalRegistry.GetObj(0).WXGet());
      wxCommandEvent com;
      pWXCryst->OnMenuCrystalGL(com);
      list<string>::iterator pos;
      for(pos=vFourierFilename.begin();pos!=vFourierFilename.end();++pos)
      {
         pWXCryst->GetCrystalGL()->LoadFourier(*pos);
      }
      return true;
   }
   
   return TRUE;
}
int MyApp::OnExit()
{
   TAU_REPORT_STATISTICS();
   return this->wxApp::OnExit();
}

// ----------------------------------------------------------------------------
// WXCrystScrolledWindow
// ----------------------------------------------------------------------------
WXCrystScrolledWindow::WXCrystScrolledWindow(wxWindow* parent):
wxScrolledWindow(parent),mpChild((wxWindow*)0),mHeight(-1),mWidth(-1)
{}

bool WXCrystScrolledWindow::Layout()
{
   if(0==mpChild) return true;
   int width,height;
   mpChild->GetSize(&width,&height);
   this->Scroll(0,0);//workaround wxMSW 2.2.9 bug
   if((mHeight!=height)||(mWidth!=width)) this->SetScrollbars(40,40,width/40+1,height/40+1);
   mHeight=height;
   mWidth=width;
   return true;
}

void WXCrystScrolledWindow::SetChild(const wxWindow* pChild)
{
   mpChild=pChild;
}

// ----------------------------------------------------------------------------
// main frame
// ----------------------------------------------------------------------------

WXCrystMainFrame::WXCrystMainFrame(const wxString& title, const wxPoint& pos, const wxSize& size,
                                   const bool splashscreen)
       : wxFrame((wxFrame *)NULL, -1, title, pos, size)
{
#ifdef __WXMAC__
   // we need this in order to allow the about menu relocation, since ABOUT is
   // not the default id of the about menu
   wxApp::s_macAboutMenuItemId = MENU_HELP_ABOUT;
#endif

   // create a menu bar
      wxMenu *menuFile = new wxMenu;//
         menuFile->Append(MENU_FILE_LOAD, "Load", "Load some objects");
         //menuFile->Append(MENU_FILE_LOAD_OXY,"Load OLD .OXY","Load using the old .oxy format");
         menuFile->Append(MENU_FILE_SAVE, "Save", "Save Everything...");
         menuFile->Append(MENU_FILE_QUIT, "E&xit\tAlt-Q", "Quit ");
      
      wxMenu *objectMenu = new wxMenu("", wxMENU_TEAROFF);
         objectMenu->Append(MENU_OBJECT_CREATE_CRYSTAL, "New Crystal",
                           "Add a new Crystal structure");
         objectMenu->Append(MENU_OBJECT_CREATE_POWDERSPECTRUM, "New PowderPattern",
                           "Add a new PowderPattern Object");
         objectMenu->Append(MENU_OBJECT_CREATE_SINGLECRYSTALDATA, "New Single Crystal Diffraction",
                           "Add a new Single Crystal Diffraction Object");
         objectMenu->Append(MENU_OBJECT_CREATE_GLOBALOPTOBJ, "New Monte-Carlo Object",
                           "Add a new Monte-Carlo Object");
      
      wxMenu *helpMenu = new wxMenu;
         helpMenu->Append(MENU_HELP_ABOUT, "&About...", "About ObjCryst...");
         helpMenu->Append(MENU_HELP_TOGGLETOOLTIP, "Toggle Tooltip", "Set Tooltips on/off");

      wxMenuBar *menuBar = new wxMenuBar();
         menuBar->Append(menuFile,  "&File");
         menuBar->Append(objectMenu,"&Objects");
         #ifdef __DEBUG__
         wxMenu *debugMenu = new wxMenu;
            debugMenu->Append(MENU_DEBUG_TEST1, "Test #1");
            debugMenu->Append(MENU_DEBUG_TEST2, "Test #2");
            debugMenu->Append(MENU_DEBUG_TEST3, "Test #3");
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

   SetMenuBar(menuBar);

#if wxUSE_STATUSBAR
   CreateStatusBar(1);
   SetStatusText("Welcome to ObjCryst++!");
#endif // wxUSE_STATUSBAR

   // Create the notebook

      wxNotebook *notebook = new wxNotebook(this, -1);

      wxLayoutConstraints* c = new wxLayoutConstraints;
      c->left.SameAs(this, wxLeft, 2);
      c->right.SameAs(this, wxRight, 2);
      c->top.SameAs(this, wxTop, 2);
      c->bottom.SameAs(this, wxBottom, 2);

      notebook->SetConstraints(c);

   // First window -Crystals
      WXCrystScrolledWindow *mpWin1 = new WXCrystScrolledWindow(notebook);
      mpWin1->SetChild(gCrystalRegistry.WXCreate(mpWin1));
      mpWin1->Layout();
      notebook->AddPage(mpWin1, "Crystals", TRUE);

   // Second window - PowderPattern
      WXCrystScrolledWindow *mpWin2 = new WXCrystScrolledWindow(notebook);
      mpWin2->SetChild(gPowderPatternRegistry.WXCreate(mpWin2));
      mpWin2->Layout();
      notebook->AddPage(mpWin2,"Powder Diffraction",true);
      
   // Third window - SingleCrystal
      WXCrystScrolledWindow *mpWin3 = new WXCrystScrolledWindow(notebook);
      mpWin3->SetChild(gDiffractionDataSingleCrystalRegistry.WXCreate(mpWin3));
      mpWin3->Layout();
      notebook->AddPage(mpWin3,"Single Crystal Diffraction",true);
      
   // Fourth window - Global Optimization
      WXCrystScrolledWindow *mpWin4 = new WXCrystScrolledWindow(notebook);
      mpWin4->SetChild(gOptimizationObjRegistry.WXCreate(mpWin4));
      mpWin4->Layout();
      notebook->AddPage(mpWin4,"Global Optimization",true);

   this->SetIcon(wxICON(Fox));
   this->Show(TRUE);
   this->Layout();
   notebook->SetSelection(0);
   //Splash Screen
   if(true==splashscreen)
   {
      wxCommandEvent event;
      this->OnAbout(event);
   }
   // Set tooltip delay
   wxToolTip::SetDelay(500);
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
               _T("Version 1.6\n\n")
               _T("(c) 2000-2003 Vincent FAVRE-NICOLIN, vincefn@users.sourceforge.net\n")
               _T("    2000-2001 Radovan CERNY, University of Geneva\n\n")
               _T("http://objcryst.sourceforge.net\n")
               _T("http://www.ccp14.ac.uk/ccp/web-mirrors/objcryst/ (Mirror)\n\n")
               _T("FOX comes with ABSOLUTELY NO WARRANTY. It is free software, and you are\n")
               _T("welcome to redistribute it under certain conditions. \n")
               _T("See the LICENSE file for details.")
              );

   wxMessageBox(msg, "About Fox", wxOK | wxICON_INFORMATION, this);
}
void WXCrystMainFrame::OnLoad(wxCommandEvent& event)
{
   wxFileDialog *open;
   if(event.GetId()==MENU_FILE_LOAD)
   {
      open= new wxFileDialog(this,"Choose File :",
                                           "","","*.xml",wxOPEN | wxFILE_MUST_EXIST);
      if(open->ShowModal() != wxID_OK) return;
      XMLCrystFileLoadAllObject(open->GetPath().c_str());
   }
   open->Destroy();
}
void WXCrystMainFrame::OnSave(wxCommandEvent& WXUNUSED(event))
{
   WXCrystValidateAllUserInput();
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
   obj->SetMaxSinThetaOvLambda(0.4);
   obj->UpdateDisplay();
}

void WXCrystMainFrame::OnAddSingleCrystalData(wxCommandEvent& WXUNUSED(event))
{
   WXCrystValidateAllUserInput();
   int choice;
   Crystal *cryst=dynamic_cast<Crystal*>
      (WXDialogChooseFromRegistry(gCrystalRegistry,(wxWindow*)this,
         "Choose a Crystal Structure:",choice));
   if(0==cryst) return;

   DiffractionDataSingleCrystal* obj;
   obj=new DiffractionDataSingleCrystal(*cryst);
   obj->SetName("Change Me!");
   obj->SetMaxSinThetaOvLambda(0.4);
   obj->UpdateDisplay();
}
void WXCrystMainFrame::OnAddGlobalOptimObj(wxCommandEvent& WXUNUSED(event))
{
   MonteCarloObj* obj;
   obj=new MonteCarloObj((string)"Change Me!");
}
void WXCrystMainFrame::OnAddGeneticAlgorithm(wxCommandEvent& WXUNUSED(event))
{
   //GeneticAlgorithm* obj;
   //obj=new GeneticAlgorithm("Change Me!");
}
void WXCrystMainFrame::OnSetDebugLevel(wxCommandEvent& event)
{
   if(event.GetId()== MENU_DEBUG_LEVEL0 ){VFN_DEBUG_GLOBAL_LEVEL(0);}
   if(event.GetId()== MENU_DEBUG_LEVEL1 ){VFN_DEBUG_GLOBAL_LEVEL(1);}
   if(event.GetId()== MENU_DEBUG_LEVEL2 ){VFN_DEBUG_GLOBAL_LEVEL(2);}
   if(event.GetId()== MENU_DEBUG_LEVEL3 ){VFN_DEBUG_GLOBAL_LEVEL(3);}
   if(event.GetId()== MENU_DEBUG_LEVEL4 ){VFN_DEBUG_GLOBAL_LEVEL(4);}
   if(event.GetId()== MENU_DEBUG_LEVEL5 ){VFN_DEBUG_GLOBAL_LEVEL(5);}
   if(event.GetId()== MENU_DEBUG_LEVEL6 ){VFN_DEBUG_GLOBAL_LEVEL(6);}
   if(event.GetId()== MENU_DEBUG_LEVEL7 ){VFN_DEBUG_GLOBAL_LEVEL(7);}
   if(event.GetId()== MENU_DEBUG_LEVEL8 ){VFN_DEBUG_GLOBAL_LEVEL(8);}
   if(event.GetId()== MENU_DEBUG_LEVEL9 ){VFN_DEBUG_GLOBAL_LEVEL(9);}
   if(event.GetId()== MENU_DEBUG_LEVEL10){VFN_DEBUG_GLOBAL_LEVEL(10);}
}
void WXCrystMainFrame::OnDebugTest(wxCommandEvent& event)
{
   WXCrystValidateAllUserInput();
   static long saveId=-1;
   static long saveId2=-1;
   if(event.GetId()== MENU_DEBUG_TEST1)
   {
      if(saveId==-1) saveId=gScattererRegistry.GetObj(0).CreateParamSet();
      else gScattererRegistry.GetObj(0).SaveParamSet(saveId);
      gScattererRegistry.GetObj(0).GlobalOptRandomMove(1);
      gCrystalRegistry.GetObj(0).UpdateDisplay();
      if(saveId2==-1) saveId2=gScattererRegistry.GetObj(0).CreateParamSet();
      else gScattererRegistry.GetObj(0).SaveParamSet(saveId2);
   }
   if(event.GetId()== MENU_DEBUG_TEST2)
   {
      gScattererRegistry.GetObj(0).RestoreParamSet(saveId);
      gCrystalRegistry.GetObj(0).UpdateDisplay();
   }
   if(event.GetId()== MENU_DEBUG_TEST3)
   {
      Crystal *cryst=new Crystal(25.,30.,35.,"P1");
      ScatteringPowerAtom *ScattPowS=new ScatteringPowerAtom("S" ,"S",0.74);
      ScatteringPowerAtom *ScattPowO=new ScatteringPowerAtom("O","O",1.87);
      cryst->AddScatteringPower(ScattPowS);
      cryst->AddScatteringPower(ScattPowO);
      Molecule *mol;
      mol=MakeTetrahedron(*cryst,"SO4",ScattPowS,ScattPowO,1.5);
      mol->RestraintStatus(cout);cryst->AddScatterer(mol);
      
      mol=MakeOctahedron(*cryst,"SO6",ScattPowS,ScattPowO,1.5);
      mol->RestraintStatus(cout);cryst->AddScatterer(mol);
      
      mol=MakeSquarePlane(*cryst,"SO6",ScattPowS,ScattPowO,1.5);
      mol->RestraintStatus(cout);cryst->AddScatterer(mol);
      
      mol=MakeCube(*cryst,"SO8",ScattPowS,ScattPowO,1.5);
      mol->RestraintStatus(cout);cryst->AddScatterer(mol);
      
      mol=MakePrismTrigonal(*cryst,"SO6",ScattPowS,ScattPowO,1.5);
      mol->RestraintStatus(cout);cryst->AddScatterer(mol);
      
      mol=MakeIcosahedron(*cryst,"SO12",ScattPowS,ScattPowO,1.5);
      mol->RestraintStatus(cout);cryst->AddScatterer(mol);
      
      mol=MakeTriangle(*cryst,"SO3",ScattPowS,ScattPowO,1.5);
      mol->RestraintStatus(cout);cryst->AddScatterer(mol);
      
      mol=MakeAntiPrismTetragonal(*cryst,"SO8",ScattPowS,ScattPowO,1.5);
      mol->RestraintStatus(cout);cryst->AddScatterer(mol);

      cryst->RandomizeConfiguration();
      WXCrystal *pWXCryst=dynamic_cast<WXCrystal*> (gCrystalRegistry.GetObj(0).WXGet());
      wxCommandEvent com;
      pWXCryst->OnMenuCrystalGL(com);
   }
}

void WXCrystMainFrame::OnUpdateUI(wxUpdateUIEvent& event)
{
   VFN_DEBUG_MESSAGE("WXCrystMainFrame::OnUpdateUI(): Uncaught event !!",10)
}
void WXCrystMainFrame::OnToggleTooltips(wxCommandEvent& event)
{
    static bool tooltip_enabled = true;
    tooltip_enabled = !tooltip_enabled;
    wxToolTip::Enable(tooltip_enabled);
    VFN_DEBUG_MESSAGE("WXCrystMainFrame::OnToggleTooltips(): Tooltips= "<<tooltip_enabled,10)
}
///////////////////////////////////////// Speed Test////////////////////
void standardSpeedTest()
{
	cout << " Beginning Speed tests" << endl ;
   
   std::list<SpeedTestReport> vReport;
   vReport.push_back(SpeedTest(10,2,"P1"  ,RAD_NEUTRON,100,0,5.));
   vReport.push_back(SpeedTest(10,2,"P-1" ,RAD_NEUTRON,100,0,5.));
   vReport.push_back(SpeedTest(10,4,"Pnma",RAD_NEUTRON,100,0,5.));
   vReport.push_back(SpeedTest(10,4,"Ia3d",RAD_NEUTRON,100,0,5.));
   vReport.push_back(SpeedTest(10,2,"P1"  ,RAD_XRAY   ,100,0,5.));
   vReport.push_back(SpeedTest(10,2,"P-1" ,RAD_XRAY   ,100,0,5.));
   vReport.push_back(SpeedTest(10,4,"Pnma",RAD_XRAY   ,100,0,5.));
   vReport.push_back(SpeedTest(10,4,"Ia3d",RAD_XRAY   ,100,0,5.));

   vReport.push_back(SpeedTest(10,2,"P1"  ,RAD_NEUTRON,100,1,5.));
   vReport.push_back(SpeedTest(10,2,"P-1" ,RAD_NEUTRON,100,1,5.));
   vReport.push_back(SpeedTest(10,4,"Pnma",RAD_NEUTRON,100,1,5.));
   vReport.push_back(SpeedTest(10,4,"Ia3d",RAD_NEUTRON,100,1,5.));
   vReport.push_back(SpeedTest(10,2,"P1"  ,RAD_XRAY   ,100,1,5.));
   vReport.push_back(SpeedTest(10,2,"P-1" ,RAD_XRAY   ,100,1,5.));
   vReport.push_back(SpeedTest(10,4,"Pnma",RAD_XRAY   ,100,1,5.));
   vReport.push_back(SpeedTest(10,4,"Ia3d",RAD_XRAY   ,100,1,5.));

   // Results from november 2003 on Vincent's Athlon TB 1.4 GHz,
   //with gcc 3.3.1 with -O3 -ffast-math -march=athlon -funroll-all-loops
   
   
   //Spacegroup NbAtoms NbAtType Radiation Type  NbRefl  BogoSPS     BogoMRAPS   BogoMRAPS(n)
   //P1          10       2      neutron Single   100  30359.2832     30.3593     30.3593
   //P-1         10       2      neutron Single   100  32215.5703     64.4311     32.2156
   //Pnma        10       4      neutron Single   100  11115.5381     88.9243     44.4622
   //Ia3d        10       4      neutron Single   100   2775.5906    266.4567     66.6142
   //P1          10       2      X-ray   Single   100  29940.1211     29.9401     29.9401
   //P-1         10       2      X-ray   Single   100  31437.1270     62.8743     31.4371
   //Pnma        10       4      X-ray   Single   100  10914.5127     87.3161     43.6581
   //Ia3d        10       4      X-ray   Single   100   2797.6191    268.5714     67.1429
   //P1          10       2      neutron Powder   100  12155.6895     12.1557     12.1557
   //P-1         10       2      neutron Powder   100  12584.4922     25.1690     12.5845
   //Pnma        10       4      neutron Powder   100   7157.0571     57.2565     28.6282
   //Ia3d        10       4      neutron Powder   100   2470.5884    237.1765     59.2941
   //P1          10       2      X-ray   Powder   100  11749.5029     11.7495     11.7495
   //P-1         10       2      X-ray   Powder   100  12250.9961     24.5020     12.2510
   //Pnma        10       4      X-ray   Powder   100   7097.4150     56.7793     28.3897
   //Ia3d        10       4      X-ray   Powder   100   2485.2070    238.5799     59.6450
   CrystVector_REAL vfnBogoMRAPS_n_200311(16);
   vfnBogoMRAPS_n_200311(0)=30.3593;
   vfnBogoMRAPS_n_200311(1)=32.2156;
   vfnBogoMRAPS_n_200311(2)=44.4622;
   vfnBogoMRAPS_n_200311(3)=66.6142;
   vfnBogoMRAPS_n_200311(4)=29.9401;
   vfnBogoMRAPS_n_200311(5)=31.4371;
   vfnBogoMRAPS_n_200311(6)=43.6581;
   vfnBogoMRAPS_n_200311(7)=67.1429;
   vfnBogoMRAPS_n_200311(8)=12.1557;
   vfnBogoMRAPS_n_200311(9)=12.5845;
   vfnBogoMRAPS_n_200311(10)=28.6282;
   vfnBogoMRAPS_n_200311(11)=59.2941;
   vfnBogoMRAPS_n_200311(12)=11.7495;
   vfnBogoMRAPS_n_200311(13)=12.2510;
   vfnBogoMRAPS_n_200311(14)=28.3897;
   vfnBogoMRAPS_n_200311(15)=59.6450;
   
   cout<<" Spacegroup NbAtoms NbAtType Radiation Type  NbRefl  BogoSPS     BogoMRAPS   BogoMRAPS(n)  relat"<<endl;
   unsigned int i=0;
   REAL vfnCompar=0.;
   for(std::list<SpeedTestReport>::const_iterator pos=vReport.begin();
       pos != vReport.end();++pos)
   {
      cout<<"  "<<FormatString(pos->mSpacegroup,8)<<" "
          <<FormatInt(pos->mNbAtom)<<" "
          <<FormatInt(pos->mNbAtomType,6)<<"     ";
      switch(pos->mRadiation)
      {
         case(RAD_NEUTRON): cout<<FormatString("neutron",7)<<" ";break;
         case(RAD_XRAY): cout<<FormatString("X-ray",7)<<" ";break;
         case(RAD_ELECTRON): cout<<FormatString("electron",7)<<" ";break;
      }
      switch(pos->mDataType)
      {
         case(0): cout<<FormatString("Single",6)<<" ";break;
         case(1): cout<<FormatString("Powder",6)<<" ";break;
      }
      const REAL relat=pos->mBogoMRAPS_reduced/vfnBogoMRAPS_n_200311(i++);
      cout<<FormatInt(pos->mNbReflections)<<" "
          <<FormatFloat(pos->mBogoSPS)<<" "
          <<FormatFloat(pos->mBogoMRAPS)<<" "
          <<FormatFloat(pos->mBogoMRAPS_reduced)<<" "
          <<FormatFloat(relat*100.,8,2)
          <<endl;
      vfnCompar+=relat;
   }
   vfnCompar/=vfnBogoMRAPS_n_200311.numElements();
   cout<<endl<<"Your FOX/ObjCryst++ speed index is "<<FormatFloat(vfnCompar*100.,8,2)
       <<"% (100% = Athlon 1.4GHz with Linux/gcc 3.3.1, november 2003)"<<endl;
	cout << " End of Crystallographic Speeding !" << endl ;
}

