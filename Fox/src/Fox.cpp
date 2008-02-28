/*  ObjCryst++ Object-Oriented Crystallographic Library
    (c) 2000-2008 Vincent Favre-Nicolin vincefn@users.sourceforge.net
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

#ifdef __WX__CRYST__
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
   #include "wx/wfstream.h"
   #include "wx/zstream.h"
   #include "wx/fileconf.h"
   #include "wx/filesys.h"
   #include <wx/fs_inet.h>
   #include <wx/txtstrm.h>
#endif

#include <locale.h>
#include <sstream>
#include <list>

#include "ObjCryst/IO.h"
#include "ObjCryst/Crystal.h"
#include "ObjCryst/PowderPattern.h"
#include "ObjCryst/DiffractionDataSingleCrystal.h"
#include "ObjCryst/Polyhedron.h"
#include "ObjCryst/test.h"
#include "ObjCryst/CIF.h"
#include "RefinableObj/GlobalOptimObj.h"
#include "Quirks/VFNStreamFormat.h"

#ifdef __WX__CRYST__
   #include "wxCryst/wxCrystal.h"

   #if defined(__WXGTK__) || defined(__WXMOTIF__) || defined(__WXMAC__) || defined(__WXMGL__) || defined(__WXX11__)
      #include "Fox.xpm"
   #endif
#endif


using namespace ObjCryst;
using namespace std;

// Rough version number - must be updated at least for every major version or critical update
// This is used to check for updates...
//:TODO: supply __FOXREVISION__ from the command line (at least under Linux)
#define __FOXREVISION__ 964

static std::string foxVersion;

// ----------------------------------------------------------------------------
// Speed test
// ----------------------------------------------------------------------------
void standardSpeedTest();

#ifdef __WX__CRYST__
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
      void SetChild(wxWindow* pChild);
   private:
      wxWindow* mpChild;
      int mHeight,mWidth;
      wxBoxSizer *mpSizer;
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
   void OnMenuClose(wxCommandEvent& event);
   void OnClose(wxCloseEvent& event);
   void SafeClose();
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
   void OnPreferences(wxCommandEvent& event);
   void OnCheckUpdate(wxCommandEvent& event);
private:
    DECLARE_EVENT_TABLE()
    RefinableObjClock mClockLastSave;
    wxNotebook *mpNotebook;
    /// List of available updates
   std::map<unsigned int,pair<int,wxString> > mvUpdates;
   /// Are we during an autocheck ?
   bool mvUpdatesAutoCheck;
};

// ----------------------------------------------------------------------------
// For messaging the user
// ----------------------------------------------------------------------------
wxFrame *pMainFrameForUserMessage;

void WXCrystInformUserStdOut(const string &str)
{
   if(wxThread::IsMain()) pMainFrameForUserMessage->SetStatusText((wxString)str.c_str());
   else cout<<"Message for user (outside main thread):"<<str<<endl;
}

// ----------------------------------------------------------------------------
// constants
// ----------------------------------------------------------------------------
static const long MENU_FILE_QUIT=                      WXCRYST_ID();
static const long MENU_HELP_ABOUT=                     WXCRYST_ID();
static const long MENU_HELP_TOGGLETOOLTIP=             WXCRYST_ID();
static const long MENU_HELP_UPDATE=                    WXCRYST_ID();
static const long MENU_PREFS_PREFERENCES=              WXCRYST_ID();
static const long MENU_FILE_LOAD=                      WXCRYST_ID();
static const long MENU_FILE_CLOSE=                     WXCRYST_ID();
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

/// Separate thread to check for updates////////////////////////////////////////////////////////////
static const long ID_FOX_UPDATES_RESULT=               WXCRYST_ID();

class WXThreadCheckUpdates:public wxThread
{
   public:
      WXThreadCheckUpdates(std::map<unsigned int,pair<int,wxString> > &vUpdates,WXCrystMainFrame &caller):
      wxThread(wxTHREAD_DETACHED),mpvUpdates(&vUpdates),mpCaller(&caller)
      {}
      wxThread::ExitCode Entry()
      {
         //cout<<"WXThreadCheckUpdates:: OnEntry()"<<endl;
         mpvUpdates->clear();
         if(!(wxFileSystem::HasHandlerForPath("http://objcryst.sourceforge.net/FoxUpdates.txt")))
            wxFileSystem::AddHandler(new wxInternetFSHandler);
         wxFileSystem fs;
         wxFSFile *fp= NULL;
         fp= fs.OpenFile("http://objcryst.sourceforge.net/FoxUpdates.txt",wxFS_READ);
         if(fp!=NULL)
         {
            wxInputStream *fstream = fp->GetStream();
            wxTextInputStream txtis(*fstream);
            txtis.ReadLine();//first line
            while(!fstream->Eof())
            {
               unsigned int revisionfix=txtis.Read16();
               unsigned int revisionbug=txtis.Read16();
               unsigned int severity=txtis.Read16();
               wxString reason=txtis.ReadLine();
               if((revisionfix>__FOXREVISION__)&&(__FOXREVISION__>revisionbug))
               {
                  //cout<<"Revision:"<<revisionfix<<", severity="<<severity<<",reason="<<reason<<endl;
                  (*mpvUpdates)[revisionfix]=make_pair(severity,reason);
               }
            }
         }
         return NULL;
      }
      void OnExit()
      {
         //cout<<"WXThreadCheckUpdates:: OnExit()"<<endl;
         wxCommandEvent event(wxEVT_COMMAND_MENU_SELECTED,ID_FOX_UPDATES_RESULT);
         wxPostEvent(mpCaller,event);
      }
   private:
      std::map<unsigned int,pair<int,wxString> > *mpvUpdates;
      WXCrystMainFrame *mpCaller;
};

// ----------------------------------------------------------------------------
// event tables and other macros for wxWindows
// ----------------------------------------------------------------------------

BEGIN_EVENT_TABLE(WXCrystMainFrame, wxFrame)
   EVT_MENU(MENU_FILE_QUIT,                        WXCrystMainFrame::OnQuit)
   EVT_MENU(MENU_HELP_ABOUT,                       WXCrystMainFrame::OnAbout)
   EVT_MENU(MENU_HELP_TOGGLETOOLTIP,               WXCrystMainFrame::OnToggleTooltips)
   EVT_MENU(MENU_HELP_UPDATE,                      WXCrystMainFrame::OnCheckUpdate)
   EVT_MENU(ID_FOX_UPDATES_RESULT,                 WXCrystMainFrame::OnCheckUpdate)
   EVT_MENU(MENU_PREFS_PREFERENCES,                WXCrystMainFrame::OnPreferences)
   EVT_MENU(MENU_FILE_LOAD,                        WXCrystMainFrame::OnLoad)
   EVT_MENU(MENU_FILE_CLOSE,                       WXCrystMainFrame::OnMenuClose)
   EVT_CLOSE(                                      WXCrystMainFrame::OnClose)
   EVT_MENU(MENU_FILE_SAVE,                        WXCrystMainFrame::OnSave)
   EVT_MENU(MENU_OBJECT_CREATE_CRYSTAL,            WXCrystMainFrame::OnAddCrystal)
   EVT_MENU(MENU_OBJECT_CREATE_POWDERSPECTRUM,     WXCrystMainFrame::OnAddPowderPattern)
   EVT_MENU(MENU_OBJECT_CREATE_SINGLECRYSTALDATA,  WXCrystMainFrame::OnAddSingleCrystalData)
   EVT_MENU(MENU_OBJECT_CREATE_GLOBALOPTOBJ,       WXCrystMainFrame::OnAddGlobalOptimObj)
   EVT_MENU(MENU_OBJECT_CREATE_GENETICALGORITHM,   WXCrystMainFrame::OnAddGeneticAlgorithm)
   EVT_MENU(MENU_DEBUG_LEVEL0,                     WXCrystMainFrame::OnSetDebugLevel)
   EVT_MENU(MENU_DEBUG_LEVEL1,                     WXCrystMainFrame::OnSetDebugLevel)
   EVT_MENU(MENU_DEBUG_LEVEL2,                     WXCrystMainFrame::OnSetDebugLevel)
   EVT_MENU(MENU_DEBUG_LEVEL3,                     WXCrystMainFrame::OnSetDebugLevel)
   EVT_MENU(MENU_DEBUG_LEVEL4,                     WXCrystMainFrame::OnSetDebugLevel)
   EVT_MENU(MENU_DEBUG_LEVEL5,                     WXCrystMainFrame::OnSetDebugLevel)
   EVT_MENU(MENU_DEBUG_LEVEL6,                     WXCrystMainFrame::OnSetDebugLevel)
   EVT_MENU(MENU_DEBUG_LEVEL7,                     WXCrystMainFrame::OnSetDebugLevel)
   EVT_MENU(MENU_DEBUG_LEVEL8,                     WXCrystMainFrame::OnSetDebugLevel)
   EVT_MENU(MENU_DEBUG_LEVEL9,                     WXCrystMainFrame::OnSetDebugLevel)
   EVT_MENU(MENU_DEBUG_LEVEL10,                    WXCrystMainFrame::OnSetDebugLevel)
   EVT_MENU(MENU_DEBUG_TEST1,                      WXCrystMainFrame::OnDebugTest)
   EVT_MENU(MENU_DEBUG_TEST2,                      WXCrystMainFrame::OnDebugTest)
   EVT_MENU(MENU_DEBUG_TEST3,                      WXCrystMainFrame::OnDebugTest)
   EVT_UPDATE_UI(ID_CRYST_UPDATEUI,                WXCrystMainFrame::OnUpdateUI)
END_EVENT_TABLE()

IMPLEMENT_APP(MyApp)

// ============================================================================
// implementation
// ============================================================================

// 'Main program' equivalent: the program execution "starts" here
bool MyApp::OnInit()
#else
int main (int argc, char *argv[])
#endif
{
   TAU_PROFILE_SET_NODE(0); // sequential code 
   TAU_PROFILE("main()","int()",TAU_DEFAULT);
   //set locale settings to standard
   setlocale(LC_NUMERIC,"C");
   
   {// Fox version
      char verBuf[200];
      sprintf(verBuf,"1.7.5.2-#%d",__FOXREVISION__);
      foxVersion=verBuf;
   }
   bool useGUI(true);
   long nbTrial(1000000);
   long nbRun(1);
   REAL finalCost=0.;
   bool silent=false;
   string outfilename("Fox-out.xml");
   long filenameInsertCost=-1;
   bool randomize(false);
   bool only3D(false);
   bool loadFourierGRD(false);
   list<string> vFourierFilenameGRD;
   bool loadFourierDSN6(false);
   list<string> vFourierFilenameDSN6;
   bool cif2pattern=false;
   REAL cif2patternWavelength=1.54056;
   REAL cif2patternPeakWidth=0.01;
   REAL cif2patternNbPoint=1000;
   REAL cif2patternMax2Theta=M_PI*.9;
   for(int i=1;i<argc;i++)
   {
      if(string("--nogui")==string(argv[i]))
      {
         useGUI=false;
         cout << "Running Fox without GUI"<<endl;
         continue;  
      }
      if(string("--randomize")==string(argv[i]))
      {
         randomize=true;
         cout << "Randomizing parameters before running"<<endl;
         continue;  
      }
      if(string("--silent")==string(argv[i]))
      {
         silent=true;
         cout << "Running Fox quietly"<<endl;
         continue;  
      }
      if(string("--finalcost")==string(argv[i]))
      {
         ++i;
         stringstream sstr(argv[i]);
         sstr >> finalCost;
         cout << "Fox will stop after reaching cost:"<<finalCost<<endl;
         continue;  
      }
      if(string("-n")==string(argv[i]))
      {
         ++i;
         stringstream sstr(argv[i]);
         sstr >> nbTrial;
         cout << "Fox will run for "<<nbTrial<<" trials"<<endl;
         continue;
      }
      if(string("--nbrun")==string(argv[i]))
      {
         ++i;
         stringstream sstr(argv[i]);
         sstr >> nbRun;
         cout << "Fox will do "<<nbRun<<" runs, randomizing before each run"<<endl;
         continue;
      }
      if(string("--cif2pattern")==string(argv[i]))
      {
         ++i;
         cif2pattern=true;
         {
            stringstream sstr(argv[i]);
            sstr >> cif2patternWavelength;
         }
         ++i;
         {
            stringstream sstr(argv[i]);
            sstr >> cif2patternMax2Theta;
            cif2patternMax2Theta*=DEG2RAD;
            if(cif2patternMax2Theta>M_PI) cif2patternMax2Theta=M_PI;
         }
         ++i;
         {
            stringstream sstr(argv[i]);
            sstr >> cif2patternNbPoint;
         }
         ++i;
         {
            stringstream sstr(argv[i]);
            sstr >> cif2patternPeakWidth;
            cif2patternPeakWidth*=DEG2RAD;
         }
         continue;
      }
      if(string("-i")==string(argv[i]))
      {// Obsolete, just ignore
         ++i;
         continue;
      }
      if(string("-o")==string(argv[i]))
      {
         ++i;
         outfilename=string(argv[i]);
         filenameInsertCost = outfilename.find("#cost",0);
         cout <<"Fox:#cost, pos="<<filenameInsertCost<<","<<string::npos<<endl;
         if((long)(string::npos)==filenameInsertCost) filenameInsertCost=-1;
         continue;
      }
      if(string("--loadfouriergrd")==string(argv[i]))
      {
         ++i;
         loadFourierGRD=true;
         vFourierFilenameGRD.push_back(string(argv[i]));
         continue;
      }
      if(string("--loadfourierdsn6")==string(argv[i]))
      {
         ++i;
         loadFourierDSN6=true;
         vFourierFilenameDSN6.push_back(string(argv[i]));
         continue;
      }
      if(string("--only3d")==string(argv[i]))
      {
         only3D=true;
         continue;
      }
      if(string("--speedtest")==string(argv[i]))
      {
         standardSpeedTest();
         TAU_REPORT_STATISTICS();
         exit(0);
      }
      if(string("--index")==string(argv[i]))
      {
         ++i;
         ifstream f(argv[i]);
         if(!f) 
         {
            cout<<"Cannot find file to index:"<<argv[i]<<endl;
            exit(0);
         }
         PeakList pl;
         pl.ImportDhklDSigmaIntensity(f);
         f.close();
         // Set uncertainty of position lines to 1/4 of sigma ... or 0
         for(vector<PeakList::hkl>::iterator pos=pl.mvHKL.begin();pos!=pl.mvHKL.end();++pos)
         //{   pos->d2obsmin=(3*pos->d2obs+pos->d2obsmin)/4; pos->d2obsmax=(3*pos->d2obs+pos->d2obsmax)/4;}
         { pos->d2obsmin=pos->d2obs; pos->d2obsmax=pos->d2obs;}
         
         CellExplorer cx(pl,TRICLINIC,0);
         cx.SetAngleMinMax((float)90*DEG2RAD,(float)120*DEG2RAD);
         
         // Use at most 20 lines
         if(pl.GetPeakList().size()>25) pl.GetPeakList().resize(25);
         const unsigned int nb=pl.GetPeakList().size();
         
         const float dmin=pl.GetPeakList()[nb-1].dobs;
         const float dmax=pl.GetPeakList()[0].dobs/10;// /10: assume no peaks at lower resolution

         const float vmin=EstimateCellVolume(dmin,dmax,nb,TRICLINIC  ,LATTICE_P,1.2);
         const float vmax=EstimateCellVolume(dmin,dmax,nb,TRICLINIC  ,LATTICE_P,0.3);
         
         float lengthmax=pow(vmax,(float)(1/3.0))*3;
         if(lengthmax<25)lengthmax=25;
         if(lengthmax>(2/pl.GetPeakList()[0].dobs)) lengthmax=2/pl.GetPeakList()[0].dobs;
         cout<<"Indexing using TRICLINIC lattice, latt=3.0->"<<lengthmax<<"A, V="<<vmin<<"->"<<vmax<<"A^3"<<endl;
         
         cx.SetVolumeMinMax(vmin,vmax);
         cx.SetLengthMinMax(3,lengthmax);
         //cx.SetVolumeMinMax(800,1800);
         //cx.SetLengthMinMax(4,25);
         cx.DicVol(10,4,50,4);
         TAU_REPORT_STATISTICS();
         exit(0);
      }
      #ifdef __DEBUG__
      if(string("--debuglevel")==string(argv[i]))
      {
         int level;
         ++i;
         stringstream sstr(argv[i]);
         sstr >> level;
         VFN_DEBUG_GLOBAL_LEVEL(level);
         continue;
      }
      #endif
      if(string(argv[i]).find(string(".xml"))!=string::npos)
      {
         cout<<"Loading: "<<string(argv[i])<<endl;
         #ifdef __WX__CRYST__
         wxString name(argv[i]);
         if(name.Mid(name.size()-4)==wxString(".xml"))
            XMLCrystFileLoadAllObject(name.c_str());
         else
         {//compressed file
            wxFileInputStream is(name.c_str());
            wxZlibInputStream zstream(is);
            stringstream sst;
            while (!zstream.Eof()) sst<<(char)zstream.GetC();
            XMLCrystFileLoadAllObject(sst);
         }
         #else
         XMLCrystFileLoadAllObject(argv[i]);
         #endif
         if(!cif2pattern)continue;
      }
      if(string(argv[i]).find(string(".cif"))!=string::npos)
      {
         cout<<"Loading: "<<string(argv[i])<<endl;
         ifstream in (argv[i]);
         ObjCryst::CIF cif(in,true,true);
         CreateCrystalFromCIF(cif);
         CreatePowderPatternFromCIF(cif);
         if(!cif2pattern)continue;
      }
      if(cif2pattern)
      {
            for(int j=0;j<gCrystalRegistry.GetNb();++j)
            {
               PowderPattern data;
               data.SetRadiationType(RAD_XRAY);
               data.SetWavelength(cif2patternWavelength);
               data.SetPowderPatternPar(0,cif2patternMax2Theta/cif2patternNbPoint,cif2patternNbPoint);
               //add CaF2 as a Crystalline phase
               PowderPatternDiffraction * diffData=new PowderPatternDiffraction;
               diffData->SetCrystal(gCrystalRegistry.GetObj(j));
               diffData->SetReflectionProfilePar(PROFILE_PSEUDO_VOIGT,cif2patternPeakWidth*cif2patternPeakWidth);
               diffData->GetCrystal().SetUseDynPopCorr(true);
               data.SetMaxSinThetaOvLambda(50.0);
               data.AddPowderPatternComponent(*diffData);
               //we don't have data, so just simulate (0->Pi/2)..
               //give a constant 'obs pattern of unit intensity
               CrystVector_REAL obs(cif2patternNbPoint);
               obs=1;
               data.SetPowderPatternObs(obs);
               data.Prepare();
               // Save the powder pattern in text format
               char buf [200];
               sprintf(buf,"%s_%06d.dat",argv[i],j);
               cout<<"Auto-simulating powder pattern:"<<endl
                   <<"   Crystal #"<<j<<": "<<gCrystalRegistry.GetObj(j).GetName()<<endl
                   <<"   Wavelength: "<<cif2patternWavelength<<endl
                   <<"   2theta: 0->"<<cif2patternMax2Theta*RAD2DEG<<"° ("<<cif2patternNbPoint<<" points)"<<endl
                   <<"   peak width: "<<cif2patternPeakWidth*RAD2DEG<<"°"<<endl
                   <<"   to FILE:"<<buf<<endl;
               ofstream out(buf);
               CrystVector_REAL ttheta,icalc;
               icalc=data.GetPowderPatternCalc();
               icalc*=100/icalc.max();
               ttheta=data.GetPowderPatternX();
               if(data.GetRadiation().GetWavelengthType()!=WAVELENGTH_TOF) ttheta *= RAD2DEG;
               out << "#Simulated data for crystal:"<<gCrystalRegistry.GetObj(j).GetName() << endl;
               out << "#    2Theta/TOF    ICalc" << endl;
               out << FormatVertVector<REAL>(ttheta,icalc,12,4);
               out.close();
            }
            // Erase every data in memory
         gOptimizationObjRegistry.DeleteAll();
         gDiffractionDataSingleCrystalRegistry.DeleteAll();
         gPowderPatternRegistry.DeleteAll();
         gCrystalRegistry.DeleteAll();
         continue;
      }
      cout <<"command-line arguments:"<<endl
           <<"   in.xml: input 'in.xml' file"<<endl
           <<"   --loadfouriergrd map.grd: load and display 'map.grd' fourier map with (first) crystal structure"<<endl
           <<"   --loadfourierdsn6 map.DN6: load and display a DSN6 fourier map with (first) crystal structure"<<endl
           <<"   --nogui: run without GUI, automatically launches optimization"<<endl
           <<"      options with --nogui:"<<endl
           <<"         -n 10000     : run for 10000 trials at most (default: 1000000)"<<endl
           <<"         --nbrun 5     : do 5 runs, randomizing before each run (default: 1), use -1 to run indefinitely"<<endl
           <<"         -o out.xml   : output in 'out.xml'"<<endl
           <<"         --randomize  : randomize initial configuration"<<endl
           <<"         --silent     : (almost) no text output"<<endl
           <<"         --finalcost 0.15 : run optimization until cost < 0.15"<<endl
           <<"         --cif2pattern 1.5406 170 5000 .1 : simulate pattern for input crystal, wavelength=1.5406"<<endl
           <<"                                            up to 170° with 500 points and a peak width of 0.1°"<<endl
           <<endl<<endl<<"           EXAMPLES :"<<endl<<endl
           <<"Load file 'silicon.xml' and launch GUI:"<<endl<<endl
           <<"    Fox silicon.xml"<<endl<<endl
           
           <<"Load file 'ktartrate.xml', randomize, then make 1 optimization of "<<endl
           <<"1 million trials, and save the best structure in 'best.xml' :"<<endl<<endl
           <<"    Fox Cimetidine-powder.xml --nogui --randomize -n 1000000 -o best.xml"<<endl<<endl

           <<"Load file 'Cimetidine-powder.xml', then make 10 runs (starting from "<<endl
           <<"a random structure) of 10 million trials (each run saves one xml file)"<<endl
           <<", and save the best structure in 'best.xml' :"<<endl<<endl
           <<"    Fox Cimetidine-powder.xml --nogui --randomize -n 10000000 --nbrun 10 -o best.xml"<<endl<<endl
           
           <<"Load file 'Cimetidine-powder.xml', then make 10 silent runs of 10 million trials"<<endl
           <<" (each run saves one xml file), and save the best structure in 'best.xml'."<<endl
           <<" For each run, the optimization stops if the cost goes below 200000."<<endl<<endl
           <<"    Fox Cimetidine-powder.xml --nogui --silent --randomize -n 10000000 --nbrun 10 --finalcost 200000 -o best.xml"<<endl<<endl
           <<endl;
      exit(0);  
   }
   
   if(cif2pattern)
   {
      #ifdef __WX__CRYST__
      this->OnExit();
      #endif
      exit(0);
   }
   if(randomize)
      for(int i=0;i<gOptimizationObjRegistry.GetNb();i++)
         gOptimizationObjRegistry.GetObj(i).RandomizeStartingConfig();
   
#ifndef __WX__CRYST__
   useGUI=false;
#endif
   if(!useGUI)
   {
      if(nbTrial!=0)
      {
         if(nbRun==1)
         {
            for(int i=0;i<gOptimizationObjRegistry.GetNb();i++)
                  gOptimizationObjRegistry.GetObj(i).Optimize(nbTrial,silent,finalCost);
         }
         else
         {
            for(int i=0;i<gOptimizationObjRegistry.GetNb();i++)
               gOptimizationObjRegistry.GetObj(i).GetXMLAutoSaveOption().SetChoice(5);
            for(int i=0;i<gOptimizationObjRegistry.GetNb();i++)
               gOptimizationObjRegistry.GetObj(i).MultiRunOptimize(nbRun,nbTrial,silent,finalCost);
         }
      }
      string tmpstr=outfilename;
      if(filenameInsertCost>=0)
      {
         char costAsChar[50];
         sprintf(costAsChar,"-Cost-%f",gOptimizationObjRegistry.GetObj(0).GetLogLikelihood());
         string tmpstr2=costAsChar;
         tmpstr.replace(filenameInsertCost,5,tmpstr2,0,tmpstr2.length());
      }
      XMLCrystFileSaveGlobal(tmpstr);
      cout <<"End of Fox execution. Bye !"<<endl;
      #ifdef __WX__CRYST__
      this->OnExit();
      exit(0);
      #endif
   }
#ifdef __WX__CRYST__
   this->SetVendorName(_T("http://objcryst.sf.net/Fox"));
   this->SetAppName(_T("FOX-Free Objects for Crystallography"));
   // Read (and automatically create if necessary) global Fox preferences
   // We explicitely use a wxFileConfig, to avoid the registry under Windows
   wxConfigBase::Set(new wxFileConfig("FOX-Free Objects for Crystallography"));

   if(wxConfigBase::Get()->HasEntry("Fox/BOOL/Enable tooltips"))
   {
       bool tooltip_enabled;
       wxConfigBase::Get()->Read("Fox/BOOL/Enable tooltips", &tooltip_enabled);
       wxToolTip::Enable(tooltip_enabled);
   }
   else wxConfigBase::Get()->Write("Fox/BOOL/Enable tooltips", true);
   
   if(!wxConfigBase::Get()->HasEntry("Fox/BOOL/Ask confirmation before exiting Fox"))
      wxConfigBase::Get()->Write("Fox/BOOL/Ask confirmation before exiting Fox", true);
   
   if(!wxConfigBase::Get()->HasEntry("Fox/BOOL/Use compressed file format (.xml.gz)"))
      wxConfigBase::Get()->Write("Fox/BOOL/Use compressed file format (.xml.gz)", true);

   if(!wxConfigBase::Get()->HasEntry("Fox/BOOL/Check for Fox updates"))
      wxConfigBase::Get()->Write("Fox/BOOL/Check for Fox updates", true);

   WXCrystMainFrame *frame ;
   string title(string("FOX: Free Objects for Xtal structures v")+foxVersion);
   frame = new WXCrystMainFrame(title.c_str(),
                                 wxPoint(50, 50), wxSize(550, 400),
                                 !(loadFourierGRD||loadFourierDSN6));
   // Use the main frame status bar to pass messages to the user
      pMainFrameForUserMessage=frame;
      fpObjCrystInformUser=&WXCrystInformUserStdOut;
   
   WXCrystal *pWXCryst;
   if(loadFourierGRD || loadFourierDSN6)
   {
      //wxFrame *pWXFrame= new wxFrame((wxFrame *)NULL, -1, "FOX", wxPoint(50, 50), wxSize(550, 400));
      //wxScrolledWindow *pWXScWin=new wxScrolledWindow(pWXFrame,-1);
      //WXCrystal * pWXCrystal=new WXCrystal(pWXScWin,&(gCrystalRegistry.GetObj(0)));
      pWXCryst=dynamic_cast<WXCrystal*> (gCrystalRegistry.GetObj(0).WXGet());
      wxCommandEvent com;
      pWXCryst->OnMenuCrystalGL(com);
   }
   if(loadFourierGRD)
   {
      list<string>::iterator pos;
      for(pos=vFourierFilenameGRD.begin();pos!=vFourierFilenameGRD.end();++pos)
      {
         UnitCellMap *pMap=new UnitCellMap(pWXCryst->GetCrystal());
         cout<<"Reading Fourier file:"<<*pos<<endl;
         if (pMap->ImportGRD(*pos) == 0)
         {
            cout<<"Error reading Fourier file:"<< *pos<<endl;
            return FALSE;
         }
         pWXCryst->GetCrystalGL()->AddFourier(pMap);
      }
   }
   if(loadFourierDSN6)
   {
      list<string>::iterator pos;
      for(pos=vFourierFilenameDSN6.begin();pos!=vFourierFilenameDSN6.end();++pos)
      {
         UnitCellMap *pMap=new UnitCellMap(pWXCryst->GetCrystal());
         cout<<"Reading Fourier file:"<<*pos<<endl;
         if (pMap->ImportDSN6(*pos) == 0)
         {
            cout<<"Error reading Fourier file:"<< *pos<<endl;
            return FALSE;
         }
         pWXCryst->GetCrystalGL()->AddFourier(pMap);
      }
   }
   return TRUE;
#else
   return 0;
#endif
}
#ifdef __WX__CRYST__
int MyApp::OnExit()
{
   delete wxConfigBase::Set((wxConfigBase *) NULL);
   TAU_REPORT_STATISTICS();
   return this->wxApp::OnExit();
}

// ----------------------------------------------------------------------------
// WXCrystScrolledWindow
// ----------------------------------------------------------------------------
WXCrystScrolledWindow::WXCrystScrolledWindow(wxWindow* parent):
wxScrolledWindow(parent),mpChild((wxWindow*)0),mHeight(-1),mWidth(-1)
{
   mpSizer=new wxBoxSizer(wxHORIZONTAL);
   this->SetSizer(mpSizer);
   this->FitInside();
   this->SetScrollRate(10,10);
}

bool WXCrystScrolledWindow::Layout()
{
   //this->Scroll(0,0);//workaround bug ?
   return this->wxScrolledWindow::Layout();
}

void WXCrystScrolledWindow::SetChild(wxWindow* pChild)
{
   mpChild=pChild;
   mpSizer->Add(mpChild);
   // Initialize scrollbars
   //this->SetScrollbars(40,40,2,2);
}

// ----------------------------------------------------------------------------
// main frame
// ----------------------------------------------------------------------------

WXCrystMainFrame::WXCrystMainFrame(const wxString& title, const wxPoint& pos, const wxSize& size,
                                   const bool splashscreen)
       : wxFrame((wxFrame *)NULL, -1, title, pos, size),mvUpdatesAutoCheck(false)
{
#ifdef __WXMAC__
   // we need this in order to allow the about menu relocation, since ABOUT is
   // not the default id of the about menu
   wxApp::s_macAboutMenuItemId = MENU_HELP_ABOUT;
#endif

   // create a menu bar
      wxMenu *menuFile = new wxMenu;//
         menuFile->Append(MENU_FILE_LOAD, "&Open .xml or .cif\tCtrl-O", "Open Fox (.xml, .xml.gz) or CIF file");
         menuFile->Append(MENU_FILE_CLOSE, "Close\tCtrl-W", "Close all");
         menuFile->Append(MENU_FILE_SAVE, "&Save\tCtrl-S", "Save Everything...");
         menuFile->Append(MENU_FILE_QUIT, "E&xit\tCtrl-Q", "Quit ");
      
      wxMenu *objectMenu = new wxMenu("", wxMENU_TEAROFF);
         objectMenu->Append(MENU_OBJECT_CREATE_CRYSTAL, "New Crystal",
                           "Add a new Crystal structure");
         objectMenu->Append(MENU_OBJECT_CREATE_POWDERSPECTRUM, "New PowderPattern",
                           "Add a new PowderPattern Object");
         objectMenu->Append(MENU_OBJECT_CREATE_SINGLECRYSTALDATA, "New Single Crystal Diffraction",
                           "Add a new Single Crystal Diffraction Object");
         objectMenu->Append(MENU_OBJECT_CREATE_GLOBALOPTOBJ, "New Monte-Carlo Object",
                           "Add a new Monte-Carlo Object");
      
      wxMenu *prefsMenu = new wxMenu;
         prefsMenu->Append(MENU_PREFS_PREFERENCES, "&Preferences...", "Fox Preferences...");
      
      wxMenu *helpMenu = new wxMenu;
         helpMenu->Append(MENU_HELP_ABOUT, "&About...", "About ObjCryst...");
         helpMenu->Append(MENU_HELP_TOGGLETOOLTIP, "Toggle Tooltips", "Set Tooltips on/off");
         helpMenu->Append(MENU_HELP_UPDATE, "Check for Updates", "Check for a newer version of Fox");
      wxMenuBar *menuBar = new wxMenuBar();
         menuBar->Append(menuFile,  "&File");
         menuBar->Append(objectMenu,"&Objects");
         menuBar->Append(prefsMenu, "&Preferences");
         menuBar->Append(helpMenu,  "&Help");
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

   SetMenuBar(menuBar);

#if wxUSE_STATUSBAR
   CreateStatusBar(1);
   SetStatusText("Welcome to FOX/ObjCryst++!");
#endif // wxUSE_STATUSBAR

   //Splash Screen
   if(true==splashscreen)
   {
      wxCommandEvent event;
      this->OnAbout(event);
   }
   // Create the notebook

      mpNotebook = new wxNotebook(this, -1);

      wxLayoutConstraints* c = new wxLayoutConstraints;
      c->left.SameAs(this, wxLeft, 2);
      c->right.SameAs(this, wxRight, 2);
      c->top.SameAs(this, wxTop, 2);
      c->bottom.SameAs(this, wxBottom, 2);

      mpNotebook->SetConstraints(c);

   // First window -Crystals
      WXCrystScrolledWindow *mpWin1 = new WXCrystScrolledWindow(mpNotebook);
      mpWin1->SetChild(gCrystalRegistry.WXCreate(mpWin1));
      mpWin1->Layout();
      mpNotebook->AddPage(mpWin1, "Crystals", TRUE);

   // Second window - PowderPattern
      WXCrystScrolledWindow *mpWin2 = new WXCrystScrolledWindow(mpNotebook);
      mpWin2->SetChild(gPowderPatternRegistry.WXCreate(mpWin2));
      mpWin2->Layout();
      mpNotebook->AddPage(mpWin2,"Powder Diffraction",true);
      
   // Third window - SingleCrystal
      WXCrystScrolledWindow *mpWin3 = new WXCrystScrolledWindow(mpNotebook);
      mpWin3->SetChild(gDiffractionDataSingleCrystalRegistry.WXCreate(mpWin3));
      mpWin3->Layout();
      mpNotebook->AddPage(mpWin3,"Single Crystal Diffraction",true);
      
   // Fourth window - Global Optimization
      WXCrystScrolledWindow *mpWin4 = new WXCrystScrolledWindow(mpNotebook);
      mpWin4->SetChild(gOptimizationObjRegistry.WXCreate(mpWin4));
      mpWin4->Layout();
      mpNotebook->AddPage(mpWin4,"Global Optimization",true);

   this->SetIcon(wxICON(Fox));
   this->Show(TRUE);
   this->Layout();
   mpNotebook->SetSelection(0);
   
   // Set tooltip delay
   wxToolTip::SetDelay(500);
   // Reset "last save" clock, in the case we loaded an xml file on startup
   mClockLastSave.Click();
   
   // Check for updates in a separate thread
   bool check;
   wxConfigBase::Get()->Read("Fox/BOOL/Check for Fox updates", &check);
   if(check || ((rand()%10)==0))
   {
      WXThreadCheckUpdates *pThreadCheckUpdates = new WXThreadCheckUpdates(mvUpdates,*this);
      if(pThreadCheckUpdates->Create() != wxTHREAD_NO_ERROR) 
         wxLogError("Can't create updates check thread");
      else pThreadCheckUpdates->Run();
      mvUpdatesAutoCheck=true;
   }
}

void WXCrystMainFrame::OnQuit(wxCommandEvent& WXUNUSED(event))
{
   this->SafeClose();
}

void WXCrystMainFrame::OnAbout(wxCommandEvent& WXUNUSED(event))
{
   string msg(string("F.O.X. - Free Objects for Xtallography\n")
              +"Version "+ foxVersion +" \n\n"
              +"(c) 2000-2008 Vincent FAVRE-NICOLIN, vincefn@users.sourceforge.net\n"
              +"    2000-2001 Radovan CERNY, University of Geneva\n\n"
              +"http://objcryst.sourceforge.net\n\n"
              +"FOX comes with ABSOLUTELY NO WARRANTY. It is free software, and you are\n"
              +"welcome to redistribute it under certain conditions. \n"
              +"See the LICENSE file for details.");

   wxMessageDialog ab(this,msg.c_str(), "About Fox", wxOK | wxICON_INFORMATION | wxSTAY_ON_TOP );
   ab.ShowModal();
}
void WXCrystMainFrame::OnLoad(wxCommandEvent& event)
{
   // First record if any object already in memory need to be saved
   bool saved=true;
   for(int i=0;i<gRefinableObjRegistry.GetNb();i++)
      if(gRefinableObjRegistry.GetObj(i).GetClockMaster()>mClockLastSave)
      {
         saved=false;
         break;
      }

   wxFileDialog *open;
   if(event.GetId()==MENU_FILE_LOAD)
   {
      open= new wxFileDialog(this,"Choose File :",
                             "","","FOX files (*.xml,*.xml.gz) or CIF (*.cif)|*.xml;*.xml.gz;*.cif",wxOPEN | wxFILE_MUST_EXIST);
      if(open->ShowModal() != wxID_OK) return;
      wxString name=open->GetPath();
      if(name.Mid(name.size()-4)==wxString(".xml"))
         XMLCrystFileLoadAllObject(name.c_str());
      else
         if(name.Mid(name.size()-4)==wxString(".cif"))
         {
            ifstream in (open->GetPath().c_str());
            ObjCryst::CIF cif(in,true,true);
            CreateCrystalFromCIF(cif);
            CreatePowderPatternFromCIF(cif);
         }
         else
            if(name.size()>6)
               if(name.Mid(name.size()-6)==wxString("xml.gz"))
               {//compressed file
                  wxFileInputStream is(name.c_str());
                  wxZlibInputStream zstream(is);
                  stringstream sst;
                  while (!zstream.Eof()) sst<<(char)zstream.GetC();
                  XMLCrystFileLoadAllObject(sst);
               }
   }
   open->Destroy();
   if(saved) mClockLastSave.Click();
}

void WXCrystMainFrame::OnMenuClose(wxCommandEvent& event)
{
   bool safe=true;
   wxConfigBase::Get()->Read("Fox/BOOL/Ask confirmation before exiting Fox",&safe);
   if(safe)
   {
      bool saved=true;
      for(int i=0;i<gRefinableObjRegistry.GetNb();i++)
         if(gRefinableObjRegistry.GetObj(i).GetClockMaster()>mClockLastSave)
         {
            saved=false;
            break;
         }
      if(!saved)
      {
         wxString msg;
         msg.Printf( _T("Some objects have not been saved\n")
                  _T("Do you really want to close all ?"));
   
         wxMessageDialog d(this,msg, "Really Close ?", wxYES | wxNO);
         if(wxID_YES!=d.ShowModal()) return;
      }
   }
   cout<<"Removing all Optimization objects..."<<endl;
   gOptimizationObjRegistry.DeleteAll();
   cout<<"Removing all DiffractionDataSingleCrystal objects..."<<endl;
   gDiffractionDataSingleCrystalRegistry.DeleteAll();
   cout<<"Removing all PowderPattern objects..."<<endl;
   gPowderPatternRegistry.DeleteAll();
   cout<<"Removing all Crystal objects..."<<endl;
   gCrystalRegistry.DeleteAll();
}
void WXCrystMainFrame::OnClose(wxCloseEvent& event)
{
   this->SafeClose();
}
void WXCrystMainFrame::SafeClose()
{
   bool safe=true;
   wxConfigBase::Get()->Read("Fox/BOOL/Ask confirmation before exiting Fox",&safe);
   if(!safe)
   {
      this->Destroy();
      return;
   }

   bool saved=true;
   for(int i=0;i<gRefinableObjRegistry.GetNb();i++)
      if(gRefinableObjRegistry.GetObj(i).GetClockMaster()>mClockLastSave)
      {
         saved=false;
         break;
      }
   if(!saved)
   {
      wxString msg;
      msg.Printf( _T("Some objects have not been saved\n")
               _T("Do you really want to Exit ?"));

      wxMessageDialog d(this,msg, "Really Exit ?", wxYES | wxNO);
      if(wxID_YES!=d.ShowModal()) return;
   }
   cout<<"Removing all Optimization objects..."<<endl;
   gOptimizationObjRegistry.DeleteAll();
   cout<<"Removing all DiffractionDataSingleCrystal objects..."<<endl;
   gDiffractionDataSingleCrystalRegistry.DeleteAll();
   cout<<"Removing all PowderPattern objects..."<<endl;
   gPowderPatternRegistry.DeleteAll();
   cout<<"Removing all Crystal objects..."<<endl;
   gCrystalRegistry.DeleteAll();
   this->Destroy();
}
void WXCrystMainFrame::OnSave(wxCommandEvent& WXUNUSED(event))
{
   WXCrystValidateAllUserInput();
   bool compressed;
   wxConfigBase::Get()->Read("Fox/BOOL/Use compressed file format (.xml.gz)",&compressed);
   if(compressed)
   {
      wxFileDialog open(this,"Choose File to save all objects:",
                        "","","FOX compressed files (*.xml.gz)|*.xml.gz", wxSAVE | wxOVERWRITE_PROMPT);
      if(open.ShowModal() != wxID_OK) return;
      string name=open.GetPath().c_str();
      if(name.substr(name.size()-7,7)!=".xml.gz")
      {
         cout<<name<<" -> "<<name+".gz"<<endl;
         if(name.substr(name.size()-4,4)==".xml") name=name+".gz";
         else name=name+".xml.gz";
      }
      stringstream sst;
      XMLCrystFileSaveGlobal(sst);
      wxFileOutputStream ostream(name.c_str());
      wxZlibOutputStream zstream(ostream,-1,wxZLIB_GZIP);
      zstream.Write(sst.str().c_str(),sst.str().size());
   }
   else
   {
      wxFileDialog open(this,"Choose File to save all objects:",
                        "","","*.xml", wxSAVE | wxOVERWRITE_PROMPT);
      if(open.ShowModal() != wxID_OK) return;
      string name=open.GetPath().c_str();
      if(name.substr(name.size()-4,4)!=".xml")
      {
         cout<<name<<" -> "<<name+".xml"<<endl;
         name=name+".xml";
      }
      XMLCrystFileSaveGlobal(name);
   }
   mClockLastSave.Click();
}
void WXCrystMainFrame::OnAddCrystal(wxCommandEvent& WXUNUSED(event))
{
   Crystal* obj;
   obj=new Crystal;
   stringstream s;s<<"Crystal #"<<gCrystalRegistry.GetNb();
   obj->SetName(s.str());
   if(!wxConfigBase::Get()->HasEntry("Crystal/BOOL/Default-use Dynamical Occupancy Correction"))
      wxConfigBase::Get()->Write("Crystal/BOOL/Default-use Dynamical Occupancy Correction", true);
   else
   {
      bool doc;
      wxConfigBase::Get()->Read("Crystal/BOOL/Default-use Dynamical Occupancy Correction", &doc);
      if(doc) obj->GetOption(0).SetChoice(0);
      else obj->GetOption(0).SetChoice(1);
   }
   obj->UpdateDisplay();
   mpNotebook->SetSelection(0);
}
void WXCrystMainFrame::OnAddPowderPattern(wxCommandEvent& WXUNUSED(event))
{
   PowderPattern* obj;
   obj=new PowderPattern;
   stringstream s;s<<"Powder Pattern #"<<gPowderPatternRegistry.GetNb();
   obj->SetName(s.str());
   obj->SetMaxSinThetaOvLambda(0.4);
   obj->UpdateDisplay();
   mpNotebook->SetSelection(1);
}

void WXCrystMainFrame::OnAddSingleCrystalData(wxCommandEvent& WXUNUSED(event))
{
   WXCrystValidateAllUserInput();
   Crystal *cryst;
   if(gCrystalRegistry.GetNb()==1) cryst=&(gCrystalRegistry.GetObj(0));
   else
   {
      int choice;
      cryst=dynamic_cast<Crystal*>
         (WXDialogChooseFromRegistry(gCrystalRegistry,(wxWindow*)this,
            "Choose a Crystal Structure:",choice));
      if(0==cryst) return;
   }

   DiffractionDataSingleCrystal* obj;
   obj=new DiffractionDataSingleCrystal(*cryst);
   stringstream s;s<<"Diffraction data for "<<cryst->GetName();
   obj->SetName(s.str());
   obj->SetMaxSinThetaOvLambda(0.4);
   obj->UpdateDisplay();
   mpNotebook->SetSelection(2);
}
void WXCrystMainFrame::OnAddGlobalOptimObj(wxCommandEvent& WXUNUSED(event))
{
   stringstream s;s<<"OptimizationObj #"<<gOptimizationObjRegistry.GetNb();
   MonteCarloObj* obj=new MonteCarloObj(s.str());
   mpNotebook->SetSelection(3);
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
      Crystal *cryst=new Crystal(25.,30.,35.,"P1");
      ScatteringPowerAtom *ScattPowS=new ScatteringPowerAtom("S" ,"S",0.74);
      ScatteringPowerAtom *ScattPowO=new ScatteringPowerAtom("O","O",1.87);
      cryst->AddScatteringPower(ScattPowS);
      cryst->AddScatteringPower(ScattPowO);
      Molecule *mol;
      mol=MakeOctahedron(*cryst,"SO6",ScattPowS,ScattPowO,1.5);
      cryst->AddScatterer(mol);
      mol->CreateCopy();
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
   bool tooltip_enabled;
   wxConfigBase::Get()->Read("Fox/BOOL/Enable tooltips", &tooltip_enabled);
   tooltip_enabled = !tooltip_enabled;
   wxToolTip::Enable(tooltip_enabled);
   wxConfigBase::Get()->Write("Fox/BOOL/Enable tooltips", tooltip_enabled);
   VFN_DEBUG_MESSAGE("WXCrystMainFrame::OnToggleTooltips(): Tooltips= "<<tooltip_enabled,10)
}

void GetRecursiveConfigEntryList(list<pair<wxString,wxString> > &l)
{
   wxString str;
   wxString path=wxConfigBase::Get()->GetPath();
   if(path=="")path="/"+path;
   long entry;
   bool bCont = wxConfigBase::Get()->GetFirstEntry(str, entry);
   while(bCont)
   {
      //cout<<__FILE__<<":"<<__LINE__<<"Entry:"<<path+"/"+str<<endl;
      l.push_back(make_pair(path,str));
      bCont = wxConfigBase::Get()->GetNextEntry(str, entry);
   }
   long group;
   bCont = wxConfigBase::Get()->GetFirstGroup(str, group);
   while(bCont)
   {
      wxConfigBase::Get()->SetPath(path+"/"+str);
      GetRecursiveConfigEntryList(l);
      wxConfigBase::Get()->SetPath("..");
      bCont = wxConfigBase::Get()->GetNextGroup(str, group);
   }
}

enum FOX_PREF_TYPE {PREF_BOOL,PREF_STRING,PREF_LONG} ;
struct FoxPref
{
   FoxPref(wxString &comp,FOX_PREF_TYPE &t,wxString &e,wxWindow *w):
   component(comp),type(t),entry(e),win(w)
   {}
   wxString component;
   FOX_PREF_TYPE type;
   wxString entry;
   wxWindow  *win;
};

class WXFoxPreferences:public wxDialog
{
   public:
      WXFoxPreferences(wxWindow *parent);
      ~WXFoxPreferences();
      void OnClose(wxCloseEvent& event);
   private:
      list<FoxPref> l;
   DECLARE_EVENT_TABLE()
};

BEGIN_EVENT_TABLE(WXFoxPreferences, wxDialog)
   EVT_CLOSE(WXFoxPreferences::OnClose)
END_EVENT_TABLE()

WXFoxPreferences::WXFoxPreferences(wxWindow *parent):
wxDialog(parent,-1,"FOX Preferences: ",wxDefaultPosition,wxSize(400,300),wxDEFAULT_DIALOG_STYLE)
{
   wxScrolledWindow *sw=new wxScrolledWindow(this);
   sw->FitInside();
   sw->SetScrollRate(10,10);
   //sw=this;
   wxBoxSizer *sizer=new wxBoxSizer(wxVERTICAL);
   sw->SetSizer(sizer);
   list<pair<wxString,wxString> > ltmp;
   GetRecursiveConfigEntryList(ltmp);
   
   wxWindow *w;
   list<FoxPref> l2;
   for(list<pair<wxString,wxString> >::const_iterator pos=ltmp.begin();pos!=ltmp.end();++pos)
   {
      wxString component,entry;
      FOX_PREF_TYPE type;
      
      size_t tmp=pos->first.find("/",1);
      component=pos->first.substr(1,tmp-1);
      
      entry=pos->second;
      
      if(pos->first.find("BOOL"  ,1)!=wxString::npos) type=PREF_BOOL;
      if(pos->first.find("STRING",1)!=wxString::npos) type=PREF_STRING;
      if(pos->first.find("LONG"  ,1)!=wxString::npos) type=PREF_LONG;
      
      switch(type)
      {
         case PREF_BOOL:
         {
            wxCheckBox *win=new wxCheckBox(sw,-1,component+":"+entry);
            bool val;
            wxConfigBase::Get()->Read("/"+component+"/BOOL/"+entry,&val);
            win->SetValue(val);
            sizer->Add(win,0,wxLEFT);
            w=win;
            break;
         }
         case PREF_STRING:
         {
            w=new wxWindow(sw,-1);
            wxBoxSizer *s=new wxBoxSizer(wxHORIZONTAL);
            w->SetSizer(s);
            wxStaticText *txt=new wxStaticText(w,-1,component+":"+entry);
            wxString val;
            wxConfigBase::Get()->Read("/"+component+"/STRING/"+entry,&val);
            wxTextCtrl *win=new wxTextCtrl(w,-1,val);
            s->Add(txt,0,wxALIGN_CENTER);
            s->Add(win,0,wxALIGN_CENTER);
            w->Layout();
            sizer->Add(w,0,wxLEFT);
            w=win;
            break;
         }
         case PREF_LONG:
         {
            w=new wxWindow(sw,-1);
            wxBoxSizer *s=new wxBoxSizer(wxHORIZONTAL);
            w->SetSizer(s);
            wxStaticText *txt=new wxStaticText(w,-1,component+":"+entry);
            wxString val;
            wxConfigBase::Get()->Read("/"+component+"/LONG/"+entry,&val);
            wxTextCtrl *win=new wxTextCtrl(w,-1,val,wxDefaultPosition,wxDefaultSize,0,wxTextValidator(wxFILTER_NUMERIC));
            s->Add(txt,0,wxALIGN_CENTER);
            s->Add(win,0,wxALIGN_CENTER);
            w->Layout();
            sizer->Add(w,0,wxLEFT);
            w=win;
            break;
         }
      }
      l.push_back(FoxPref(component,type,entry,w));
   }
   sw->Layout();
   sizer->Fit(sw);
   this->Layout();
}
WXFoxPreferences::~WXFoxPreferences()
{
   cout<<"WXFoxPreferences::~WXFoxPreferences()"<<endl;
}

void WXFoxPreferences::OnClose(wxCloseEvent& event)
{
   cout<<"WXFoxPreferences::OnClose()"<<endl;
   for(list<FoxPref>::const_iterator pos=l.begin();pos!=l.end();++pos)
   {
      switch(pos->type)
      {
         case PREF_BOOL:
         {
            wxString full="/"+pos->component+"/BOOL/"+pos->entry;
            bool val;
            wxConfigBase::Get()->Read(full,&val);
            wxCheckBox *w=dynamic_cast<wxCheckBox *>(pos->win);
            val=w->GetValue();
            wxConfigBase::Get()->Write(full,val);
            break;
         }
         case PREF_STRING:
         {
            wxString full="/"+pos->component+"/STRING/"+pos->entry;
            wxString val;
            wxConfigBase::Get()->Read(full,&val);
            wxTextCtrl *w=dynamic_cast<wxTextCtrl *>(pos->win);
            val=w->GetValue();
            wxConfigBase::Get()->Write(full,val);
            break;
         }
         case PREF_LONG:
         {
            wxString full="/"+pos->component+"/LONG/"+pos->entry;
            wxString s;
            long val;
            wxConfigBase::Get()->Read(full,&val);
            wxTextCtrl *w=dynamic_cast<wxTextCtrl *>(pos->win);
            s=w->GetValue();
            s.ToLong(&val);
            wxConfigBase::Get()->Write(full,val);
            break;
         }
      }
   }
   this->Destroy();
}


void WXCrystMainFrame::OnPreferences(wxCommandEvent& event)
{
   WXFoxPreferences *prefs= new WXFoxPreferences(this);
   prefs->ShowModal();
}

void WXCrystMainFrame::OnCheckUpdate(wxCommandEvent& event)
{
   cout<<"WXCrystMainFrame::OnCheckUpdate"<<endl;
   if(event.GetId()==ID_FOX_UPDATES_RESULT)
   {
      unsigned int nbminorfeature=0,nbmajorfeature=0,nbrelease=0,nbminorbug=0,nbmajorbug=0,nbcritical=0;
      for(map<unsigned int,pair<int,wxString> >::const_iterator pos=mvUpdates.begin();pos!=mvUpdates.end();++pos)
      {//Critical fixes
         if(pos->second.first==0) nbminorfeature++;
         if(pos->second.first==1) nbmajorfeature++;
         if(pos->second.first==2) nbrelease++;
         if(pos->second.first==10) nbminorbug++;
         if(pos->second.first==11) nbmajorbug++;// Major bug but calculations are still OK
         if(pos->second.first==12) nbcritical++;// A mistake was made, giving erroneous results
      }
      if(mvUpdatesAutoCheck)
      {// Autocheck => only critical updates
         if(nbcritical>0)
         {
            wxString msg;
            msg.Printf( _T("A new version of Fox is available, including CRITICAL bug fixes\n")
                        _T("It is strongly recommended to update to a new version\n\n Major changes: \n"));
   
            for(map<unsigned int,pair<int,wxString> >::const_iterator pos=mvUpdates.begin();pos!=mvUpdates.end();++pos)
            {
               if(pos->second.first==12)
                  msg=msg+wxString::Format("#%d (CRITICAL): %s\n",pos->first,pos->second.second.c_str());
            }
            wxMessageDialog d(this,msg, "! CRITICAL update available !", wxOK|wxICON_EXCLAMATION);
            d.ShowModal();
         }
      }
      else
      {// User asked for this, so give the full version
         wxFrame *frame=new wxFrame(this,-1,"FOX Updates",wxDefaultPosition,wxSize(800,250),wxSTAY_ON_TOP | wxRESIZE_BORDER | wxCAPTION | wxCLOSE_BOX | wxCLIP_CHILDREN);
         wxTextCtrl *wUpdates=new wxTextCtrl(frame,-1,"",wxDefaultPosition,wxDefaultSize,wxTE_MULTILINE|wxTE_READONLY|wxTE_DONTWRAP);
         wUpdates->SetFont(wxFont(10,wxROMAN,wxFONTSTYLE_NORMAL,wxFONTWEIGHT_BOLD));
         if(mvUpdates.size()>0)
         {
            if(nbcritical>0)
            {
               cout<<wxString::Format("\n%d CRITICAL updates available:\n",nbcritical)<<endl;
               wUpdates->AppendText(wxString::Format("\n%d CRITICAL updates available:\n",nbcritical));
               for(map<unsigned int,pair<int,wxString> >::const_iterator pos=mvUpdates.begin();pos!=mvUpdates.end();++pos)
               {
                  if(pos->second.first==12)
                  {
                     wxString mess=wxString::Format("  #%d: %s\n",pos->first,pos->second.second.c_str());
                     wUpdates->AppendText(mess);
                  }
               }
            }
            if(nbrelease>0)
            {
               cout<<wxString::Format("\n%d new RELEASE available :\n",nbrelease)<<endl;
               wUpdates->AppendText(wxString::Format("\n%d new RELEASE available :\n",nbrelease));
               for(map<unsigned int,pair<int,wxString> >::const_iterator pos=mvUpdates.begin();pos!=mvUpdates.end();++pos)
               {
                  if(pos->second.first==2)
                     wUpdates->AppendText(wxString::Format("  #%d: %s\n",pos->first,pos->second.second.c_str()));
               }
            }
            if(nbmajorfeature>0)
            {
               cout<<wxString::Format("\n%d major features updates available:\n",nbmajorfeature)<<endl;
               wUpdates->AppendText(wxString::Format("\n%d major features updates available:\n",nbmajorfeature));
               for(map<unsigned int,pair<int,wxString> >::const_iterator pos=mvUpdates.begin();pos!=mvUpdates.end();++pos)
               {
                  if(pos->second.first==1)
                     wUpdates->AppendText(wxString::Format("  #%d: %s\n",pos->first,pos->second.second.c_str()));
               }
            }
            if(nbminorfeature>0)
            {
               cout<<wxString::Format("\n%d minor features updates available:\n",nbminorfeature)<<endl;
               wUpdates->AppendText(wxString::Format("\n%d minor features updates available:\n",nbminorfeature));
               for(map<unsigned int,pair<int,wxString> >::const_iterator pos=mvUpdates.begin();pos!=mvUpdates.end();++pos)
               {
                  if(pos->second.first==0)
                     wUpdates->AppendText(wxString::Format("  #%d: %s\n",pos->first,pos->second.second.c_str()));
               }
            }
         }
         else
         {
            wUpdates->AppendText(wxString::Format("No updates found !\n"));
         }
         frame->Show(true);
      }
   }
   else
   {
      mvUpdatesAutoCheck=false;
      WXThreadCheckUpdates *pThreadCheckUpdates = new WXThreadCheckUpdates(mvUpdates,*this);
      if(pThreadCheckUpdates->Create() != wxTHREAD_NO_ERROR) 
         wxLogError("Can't create updates check thread");
      else pThreadCheckUpdates->Run();
   }
}

#endif
///////////////////////////////////////// Speed Test////////////////////
void standardSpeedTest()
{
	cout << " Beginning Speed tests" << endl ;
   
   std::list<SpeedTestReport> vReport;
   vReport.push_back(SpeedTest(20,4,"P1"  ,RAD_NEUTRON,100,0,5.));
   vReport.push_back(SpeedTest(20,4,"P-1" ,RAD_NEUTRON,100,0,5.));
   vReport.push_back(SpeedTest(20,4,"Pnma",RAD_NEUTRON,100,0,5.));
   vReport.push_back(SpeedTest(20,4,"Ia3d",RAD_NEUTRON,100,0,5.));
   vReport.push_back(SpeedTest(20,4,"P1"  ,RAD_XRAY   ,100,0,5.));
   vReport.push_back(SpeedTest(20,4,"P-1" ,RAD_XRAY   ,100,0,5.));
   vReport.push_back(SpeedTest(20,4,"Pnma",RAD_XRAY   ,100,0,5.));
   vReport.push_back(SpeedTest(20,4,"Ia3d",RAD_XRAY   ,100,0,5.));

   vReport.push_back(SpeedTest(20,4,"P1"  ,RAD_NEUTRON,100,1,5.));
   vReport.push_back(SpeedTest(20,4,"P-1" ,RAD_NEUTRON,100,1,5.));
   vReport.push_back(SpeedTest(20,4,"Pnma",RAD_NEUTRON,100,1,5.));
   vReport.push_back(SpeedTest(20,4,"Ia3d",RAD_NEUTRON,100,1,5.));
   vReport.push_back(SpeedTest(20,4,"P1"  ,RAD_XRAY   ,100,1,5.));
   vReport.push_back(SpeedTest(20,4,"P-1" ,RAD_XRAY   ,100,1,5.));
   vReport.push_back(SpeedTest(20,4,"Pnma",RAD_XRAY   ,100,1,5.));
   vReport.push_back(SpeedTest(20,4,"Ia3d",RAD_XRAY   ,100,1,5.));

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

