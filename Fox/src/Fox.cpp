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

#include "wx/tooltip.h"
#include "wx/notebook.h"

#include <locale.h>
#include <sstream>
#include <list>

#include "ObjCryst/IO.h"
#include "ObjCryst/Crystal.h"
#include "ObjCryst/PowderPattern.h"
#include "ObjCryst/DiffractionDataSingleCrystal.h"
#include "RefinableObj/GlobalOptimObj.h"
//#include "RefinableObj/GeneticAlgorithm.h"
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
   wxScrolledWindow *mpWin1,*mpWin2,*mpWin3,*mpWin4;
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
              <<endl;
         exit(1);  
      }
   }
   
   
   if(randomize)
      for(int i=0;i<gOptimizationObjRegistry.GetNb();i++)
         gOptimizationObjRegistry.GetObj(i).RandomizeStartingConfig();
   
   if(!useGUI)
   {
      for(int i=0;i<gOptimizationObjRegistry.GetNb();i++)
         for(int j=0;j<5;j++)
            gOptimizationObjRegistry.GetObj(i).Optimize(nbTrial,silent,finalCost);
      XMLCrystFileSaveGlobal(outfilename);
      cout <<"End of Fox execution. Bye !"<<endl;
      exit (1);
   }
   
   WXCrystMainFrame *frame ;
   
   frame = new WXCrystMainFrame("FOX: Free Objects for Xtal structures v1.5CVS",
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
         //objectMenu->Append(MENU_OBJECT_CREATE_GENETICALGORITHM, "New Genetic Algorithm Object",
         //                  "Add a new Genetic Algorithm Object");
      
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
            debugMenu->Append(MENU_DEBUG_TEST3, "Test #2");
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
   //create something to display
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
      mpWin1->SetScrollbars( 10, 10, 0, 1000 );
      //wxBoxSizer * sizer1=new wxBoxSizer(wxVERTICAL);
      //sizer1->Add(gCrystalRegistry.WXCreate(mpWin1));
      //mpWin1->SetSizer(sizer1);
      //mpWin1->SetAutoLayout(true);
      gCrystalRegistry.WXCreate(mpWin1);
      notebook->AddPage(mpWin1, "Crystals", TRUE);

   // Second window - PowderPattern
      wxScrolledWindow *mpWin2 = new wxScrolledWindow(notebook, -1);
      mpWin2->SetScrollbars( 10, 10, 0, 500 );
      //wxBoxSizer * sizer2=new wxBoxSizer(wxVERTICAL);
      //sizer2->Add(gPowderPatternRegistry.WXCreate(mpWin2));
      //mpWin2->SetSizer(sizer2);
      //win2->SetAutoLayout(true);
      gPowderPatternRegistry.WXCreate(mpWin2);
      notebook->AddPage(mpWin2,"Powder Diffraction",true);
      
   // Third window - SingleCrystal
      wxScrolledWindow *mpWin3 = new wxScrolledWindow(notebook, -1);
      mpWin3->SetScrollbars( 10, 10, 0, 500 );
      //wxBoxSizer * sizer3=new wxBoxSizer(wxVERTICAL);
      //sizer3->Add(gPowderPatternRegistry.WXCreate(mpWin3));
      //mpWin3->SetSizer(sizer3);
      //win3->SetAutoLayout(true);
      gDiffractionDataSingleCrystalRegistry.WXCreate(mpWin3);
      notebook->AddPage(mpWin3,"Single Crystal Diffraction",true);
      
   // Fourth window - Global Optimization
      wxScrolledWindow *mpWin4 = new wxScrolledWindow(notebook, -1);
      mpWin4->SetScrollbars( 10, 10, 0, 500 );
      //wxBoxSizer * sizer4=new wxBoxSizer(wxVERTICAL);
      //sizer4->Add(gOptimizationObjRegistry.WXCreate(mpWin4));
      //mpWin4->SetSizer(sizer4);
      //mpWin4->SetAutoLayout(true);
      gOptimizationObjRegistry.WXCreate(mpWin4);
      notebook->AddPage(mpWin4,"Global Optimization",true);
   this->SetIcon(wxICON(Fox));
   this->Show(TRUE);
   this->Layout();
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
               _T("Version 1.5CVS\n\n")
               _T("(c) 2000-2002 Vincent FAVRE-NICOLIN, vincefn@users.sourceforge.net\n")
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
      gScattererRegistry.GetObj(0).RestoreParamSet(saveId2);
      gCrystalRegistry.GetObj(0).UpdateDisplay();
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
