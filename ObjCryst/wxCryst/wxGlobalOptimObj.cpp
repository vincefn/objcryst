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

#include "wx/progdlg.h"

#include "wxCryst/wxGlobalOptimObj.h"

#include "ObjCryst/IO.h"
// Next two just to fix some parameters during global optimization
#include "ObjCryst/Crystal.h"
#include "ObjCryst/ScatteringData.h"

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
//    WXOptimizationObj
//
////////////////////////////////////////////////////////////////////////
static long ID_GLOBALOPT_MENU_OBJECTS=               WXCRYST_ID(); 
static long ID_GLOBALOPT_MENU_OBJECTS_ADDOBJ=        WXCRYST_ID(); 
static long ID_GLOBALOPT_MENU_OBJECTS_REMOVEOBJ=     WXCRYST_ID(); 
static long ID_GLOBALOPT_MENU_OBJECTS_ADDCOSTFUNC=   WXCRYST_ID(); 
static long ID_GLOBALOPT_MENU_OBJECTS_REMOVECOSTFUNC=WXCRYST_ID(); 
static long ID_GLOBALOPT_MENU_OPT=                   WXCRYST_ID(); 
static long ID_GLOBALOPT_MENU_OPT_RUN=               WXCRYST_ID(); 
static long ID_GLOBALOPT_MENU_OPT_RUN_MULTIPLE=      WXCRYST_ID(); 
static long ID_GLOBALOPT_MENU_OPT_STOP=              WXCRYST_ID(); 
static long ID_GLOBALOPT_MENU_OPT_LSQ=               WXCRYST_ID(); 
static long ID_GLOBALOPT_MENU_SOLUTIONS=             WXCRYST_ID(); 
static long ID_GLOBALOPT_MENU_SOLUTIONS_BROWSE=      WXCRYST_ID(); 
static long ID_BROWSE_WIN=                           WXCRYST_ID(); 

WXOptimizationObj::WXOptimizationObj(wxWindow* parent, OptimizationObj *obj):
WXCrystObj(parent),mpGlobalOptimRunThread(0)
{
   VFN_DEBUG_ENTRY("WXOptimizationObj::WXOptimizationObj(wxWindow*,GlobalOptimObj*,)",7)
   #ifdef VFN_CRYST_MUTEX
   cout <<"new CrystMutex("<<&mMutex<<")for WXOptimizationObj:"<<obj->GetName()<<endl;
   #endif
   mpWXTitle->SetForegroundColour(wxColour(255,0,0));
   mpWXTitle->SetLabel("Global Optimization:");
   mpWXTitle->SetSize(400,-1);
   // Menu
      mpMenuBar=new WXCrystMenuBar(this,this);
      mpSizer->Add(mpMenuBar);
      mList.Add(mpMenuBar);
      
      //mpMenuBar->AddMenu("Object",ID_REFOBJ_MENU_OBJ);
         //mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_REFOBJ_MENU_OBJ_SAVE,"Save");
         //mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_REFOBJ_MENU_OBJ_LOAD,"Load");
      mpMenuBar->AddMenu("Optimized Objects",ID_GLOBALOPT_MENU_OBJECTS);
         mpMenuBar->AddMenuItem(ID_GLOBALOPT_MENU_OBJECTS,ID_GLOBALOPT_MENU_OBJECTS_ADDOBJ,
                                "Add object to optimize");
      mpMenuBar->AddMenu("Run/Stop",ID_GLOBALOPT_MENU_OPT);
         mpMenuBar->AddMenuItem(ID_GLOBALOPT_MENU_OPT,
                                ID_GLOBALOPT_MENU_OPT_RUN_MULTIPLE,"Multiple Runs");
         mpMenuBar->AddMenuItem(ID_GLOBALOPT_MENU_OPT,
                                ID_GLOBALOPT_MENU_OPT_RUN,"Single Run");
         mpMenuBar->AddMenuItem(ID_GLOBALOPT_MENU_OPT,
                                ID_GLOBALOPT_MENU_OPT_STOP,"Stop Optimization");
         mpMenuBar->GetMenu(ID_GLOBALOPT_MENU_OPT).AppendSeparator();
         mpMenuBar->AddMenuItem(ID_GLOBALOPT_MENU_OPT,
                                ID_GLOBALOPT_MENU_OPT_LSQ,"Least Squares Fit");
      mpMenuBar->AddMenu("Solutions",ID_GLOBALOPT_MENU_SOLUTIONS);
         mpMenuBar->AddMenuItem(ID_GLOBALOPT_MENU_SOLUTIONS,
                                ID_GLOBALOPT_MENU_SOLUTIONS_BROWSE,"Browse Solutions");
      mpMenuBar->Layout();
      mpSizer->SetItemMinSize(mpMenuBar,
                              mpMenuBar->GetSize().GetWidth(),
                              mpMenuBar->GetSize().GetHeight());
    //Refined Objects
   for(int i=0;i<obj->mRefinedObjList.GetNb();i++)
   {
      WXFieldName *refobj=new WXFieldName(this,"Optimized object:",this,-1,300,false);
      refobj->SetValue(obj->mRefinedObjList.GetObj(i).GetClassName()
                       +":"+obj->mRefinedObjList.GetObj(i).GetName());
      mpSizer->Add(refobj);
      mList.Add(refobj);
   }
      
   // This will be done later
   //this->CrystUpdate();
   VFN_DEBUG_EXIT("WXOptimizationObj::WXOptimizationObj(wxWindow*,GlobalOptimObj*,)",7)
}

void WXOptimizationObj::CrystUpdate(const bool uui,const bool lock)
{
   VFN_DEBUG_ENTRY("WXOptimizationObj::CrystUpdate("<<uui<<lock<<")",7)
   this->WXCrystObj::CrystUpdate(uui,lock);
   if(uui)
   {
      if(true==wxThread::IsMain()) this->UpdateUI(lock);
      else
      {
         wxUpdateUIEvent event(ID_CRYST_UPDATEUI);
         wxPostEvent(this,event);
      }
   }
   VFN_DEBUG_EXIT("WXOptimizationObj::CrystUpdate("<<uui<<lock<<")",7)
}

bool WXOptimizationObj::OnChangeName(const int id)
{
   VFN_DEBUG_MESSAGE("WXOptimizationObj::OnChangeName()",6)
   if(id==ID_WXOBJ_NAME)
   {
   VFN_DEBUG_MESSAGE("WXOptimizationObj::OnChangeName():Changing GlobalOptimObj Name",6)
      this->GetOptimizationObj().SetName(mpWXTitle->GetValue());
      return true;
   }
   return false;
}
void WXOptimizationObj::OnSave(){}

void WXOptimizationObj::OnLoad(){}

void WXOptimizationObj::OnAddRefinedObject(wxCommandEvent & WXUNUSED(event))
{
   WXCrystValidateAllUserInput();
   int choice;
   RefinableObj *obj=
      WXDialogChooseFromRegistry(gTopRefinableObjRegistry,this,
                                 "Choose object to optimize:",choice);
   if(0==obj) return;
   this->GetOptimizationObj().AddRefinableObj(*obj);
}

void WXOptimizationObj::AddRefinedObject(RefinableObj &obj)
{
   WXCrystValidateAllUserInput();
   WXFieldName *refobj=new WXFieldName(this,"Optimized object:",this,-1,300,false);
   refobj->SetValue(obj.GetClassName()+":"+obj.GetName());
   mpSizer->Add(refobj);
   mList.Add(refobj);
   this->BottomLayout(0);
   refobj->UpdateUI(true);
}

void WXOptimizationObj::OnRemoveRefinedObject(wxCommandEvent & WXUNUSED(event))
{
}

void WXOptimizationObj::OnStopOptimization(wxCommandEvent & WXUNUSED(event))
{
   this->GetOptimizationObj().StopAfterCycle();
}
void WXOptimizationObj::OnUpdateUI(wxUpdateUIEvent& event)
{
   VFN_DEBUG_ENTRY("WXOptimizationObj::OnUpdateUI()",5)
   this->UpdateUI(true);
   VFN_DEBUG_EXIT("WXOptimizationObj::OnUpdateUI()",5)
}

void WXOptimizationObj::UpdateUI(const bool lock)
{
   VFN_DEBUG_ENTRY("WXOptimizationObj::UpdateUI()",5)
   if(lock) mMutex.Lock();
   mpWXTitle->SetValue(this->GetOptimizationObj().GetName());
   mpWXTitle->UpdateUI(false);
   this->WXCrystObj::UpdateUI(false);
   if(lock) mMutex.Unlock();
   VFN_DEBUG_EXIT("WXOptimizationObj::UpdateUI()",5)
}

void WXOptimizationObj::OnBrowseParamSet(wxCommandEvent & WXUNUSED(event))
{
   if(this->GetOptimizationObj().IsOptimizing())
   {
      wxMessageDialog dumbUser(this,_T("Cannot browse during Optimisation !"),
                               _T("Cannot browse during Optimisation!"),wxOK|wxICON_EXCLAMATION);
      dumbUser.ShowModal();
      return;
   }
   wxFrame *frame= new wxFrame(this,-1,_T("Stored Configurations"),
                               wxDefaultPosition,wxSize(250,200));
   const long nb=this->GetOptimizationObj().mvSavedParamSet.size();
   wxString *choices = new wxString[nb];
   for(int i=0;i<nb;i++)
   {
      wxString tmpname=wxString::FromAscii(this->GetOptimizationObj().mRefParList.GetParamSetName
                                             (this->GetOptimizationObj().mvSavedParamSet[i].first).c_str());
      choices[i].Printf(_T("%d, cost= %f, %s"),i,
                         this->GetOptimizationObj().mvSavedParamSet[i].second,tmpname.c_str());
      //cout<<choices[i]<<endl;
   }
   mpwxParamSetList=new wxListBox(frame, ID_BROWSE_WIN, wxDefaultPosition, 
                                   wxDefaultSize, nb, choices,
                                   wxLB_SINGLE|wxLB_NEEDED_SB, wxDefaultValidator,
                                   _T("listBox"));
   mpwxParamSetList->SetEventHandler(this);
   mClockParamSetWindow.Click();
   frame->Show(true);
}

void WXOptimizationObj::OnSelectParamSet(wxCommandEvent &event)
{
   if(this->GetOptimizationObj().IsOptimizing())
   {
      wxMessageDialog dumbUser(this,_T("Cannot browse during Optimisation !"),
                               _T("Cannot browse during Optimisation!"),wxOK|wxICON_EXCLAMATION);
      dumbUser.ShowModal();
      return;
   }
   const long n=event.GetSelection();
   if(mClockParamSetWindow>this->GetOptimizationObj().mRefParList.GetRefParListClock())
   {
      try
      {
         this->GetOptimizationObj().mRefParList
            .RestoreParamSet(this->GetOptimizationObj().mvSavedParamSet[n].first);
      }
      catch(const ObjCrystException &except)
      {
         wxMessageDialog bad(this,_T("Impossible ! Model has been altered !"),
                                  _T("Impossible ! Model has been altered !"),wxOK|wxICON_EXCLAMATION);
         mpwxParamSetList->GetParent()->Close();
         mpwxParamSetList=0;
      }
      this->GetOptimizationObj().UpdateDisplay();
      cout <<"Param set #"<<this->GetOptimizationObj().mvSavedParamSet[n].first<<", cost="
           <<this->GetOptimizationObj().mvSavedParamSet[n].second
           <<", now cost="<<this->GetOptimizationObj().GetLogLikelihood()<<endl;
   }
   else
   {
      wxMessageDialog bad(this,_T("Impossible ! The list of parameters has been changed !"),
                               _T("Impossible ! The list of parameters has been changed !"),wxOK|wxICON_EXCLAMATION);
      bad.ShowModal();
   }
}

////////////////////////////////////////////////////////////////////////
//
//    WXGlobalOptimRunThread
//
////////////////////////////////////////////////////////////////////////
WXGlobalOptimRunThread::WXGlobalOptimRunThread(OptimizationObj &globalOptObj,long &nbTrial,
                             const REAL finalCost,long &nbRun,const bool multiple):
wxThread(wxTHREAD_DETACHED),mpGlobalOptObj(&globalOptObj),mpNbTrial(&nbTrial),mpNbRun(&nbRun),
mFinalCost(finalCost),mDoMultiple(multiple)
{
}

void *WXGlobalOptimRunThread::Entry()
{
   cout <<endl<<"Entering refinement thread "<<endl<<endl;
   try{
      if(mDoMultiple) mpGlobalOptObj->MultiRunOptimize(*mpNbRun,*mpNbTrial,false,mFinalCost);
      else mpGlobalOptObj->Optimize(*mpNbTrial,false,mFinalCost);
   }
   catch(...){cout<<"Unhandled exception in WXGlobalOptimRunThread::Entry() ?"<<endl;}
   return NULL;
}

void WXGlobalOptimRunThread::OnExit()
{
   cout <<endl<<"Exiting refinement thread "<<endl<<endl;
   XMLCrystFileSaveGlobal("Fox-LastOptimizationStop.xml");
}
////////////////////////////////////////////////////////////////////////
//
//    WXMonteCarloObj
//
////////////////////////////////////////////////////////////////////////
BEGIN_EVENT_TABLE(WXMonteCarloObj, wxWindow)
   EVT_BUTTON(ID_WXOBJ_COLLAPSE,                       WXCrystObj::OnToggleCollapse)
   //EVT_MENU(ID_REFOBJ_MENU_OBJ_SAVE,                   WXOptimizationObj::OnSave)
   //EVT_MENU(ID_REFOBJ_MENU_OBJ_LOAD,                   WXOptimizationObj::OnLoad)
   EVT_MENU(ID_GLOBALOPT_MENU_OBJECTS_ADDOBJ,    WXOptimizationObj::OnAddRefinedObject)
   EVT_MENU(ID_GLOBALOPT_MENU_OPT_RUN,           WXOptimizationObj::OnRunOptimization)
   EVT_MENU(ID_GLOBALOPT_MENU_OPT_RUN_MULTIPLE,  WXOptimizationObj::OnRunOptimization)
   EVT_MENU(ID_GLOBALOPT_MENU_OPT_STOP,          WXOptimizationObj::OnStopOptimization)
   EVT_MENU(ID_GLOBALOPT_MENU_OPT_LSQ,           WXMonteCarloObj::OnLSQRefine)
   EVT_MENU(ID_GLOBALOPT_MENU_SOLUTIONS_BROWSE,  WXOptimizationObj::OnBrowseParamSet)
   EVT_UPDATE_UI(ID_CRYST_UPDATEUI,              WXOptimizationObj::OnUpdateUI)
   EVT_LISTBOX(ID_BROWSE_WIN,                    WXOptimizationObj::OnSelectParamSet)
   EVT_LISTBOX_DCLICK(ID_BROWSE_WIN,             WXOptimizationObj::OnSelectParamSet)
END_EVENT_TABLE()

WXMonteCarloObj::WXMonteCarloObj(wxWindow *parent, MonteCarloObj* obj):
WXOptimizationObj(parent,obj),mpMonteCarloObj(obj),mNbTrial(10000000),mNbRun(-1)
{
   VFN_DEBUG_ENTRY("WXMonteCarloObj::WXMonteCarloObj()",7)
   //options
      WXFieldOption *opt;
      
      opt=new WXFieldOption(this,-1,&(mpMonteCarloObj->mGlobalOptimType));
      mpSizer->Add(opt,0,wxALIGN_LEFT);
      mList.Add(opt);
      
      opt=new WXFieldOption(this,-1,&(mpMonteCarloObj->mAnnealingScheduleTemp));
      mpSizer->Add(opt,0,wxALIGN_LEFT);
      mList.Add(opt);

      wxBoxSizer *sizerTemp=new wxBoxSizer(wxHORIZONTAL);
      WXFieldPar<REAL> *tempMax=
         new WXFieldPar<REAL>(this,"Temperature Max:",-1,&(mpMonteCarloObj->mTemperatureMax));
      WXFieldPar<REAL> *tempMin=
         new WXFieldPar<REAL>(this,"Temperature Min:",-1,&(mpMonteCarloObj->mTemperatureMin));
      WXFieldPar<REAL> *tempGamma=
         new WXFieldPar<REAL>(this,"Gamma:",-1,&(mpMonteCarloObj->mTemperatureGamma));
      sizerTemp->Add(tempMax);
      sizerTemp->Add(tempMin);
      sizerTemp->Add(tempGamma);
      mpSizer->Add(sizerTemp);
      mList.Add(tempMax);
      mList.Add(tempMin);
      mList.Add(tempGamma);
      tempMax->SetFormat(_T("%8f"));
      tempMin->SetFormat(_T("%8f"));
      tempGamma->SetFormat(_T("%8f"));
      
      opt=new WXFieldOption(this,-1,&(mpMonteCarloObj->mAnnealingScheduleMutation));
      mpSizer->Add(opt,0,wxALIGN_LEFT);
      mList.Add(opt);
      
      wxBoxSizer *sizerAmp=new wxBoxSizer(wxHORIZONTAL);
      WXFieldPar<REAL> *ampMax=
         new WXFieldPar<REAL>(this,"Amplitude Max:",-1,&(mpMonteCarloObj->mMutationAmplitudeMax));
      WXFieldPar<REAL> *ampMin=
         new WXFieldPar<REAL>(this,"Amplitude Min:",-1,&(mpMonteCarloObj->mMutationAmplitudeMin));
      WXFieldPar<REAL> *ampGamma=
         new WXFieldPar<REAL>(this,"Gamma:",-1,&(mpMonteCarloObj->mMutationAmplitudeGamma));
      sizerAmp->Add(ampMax);
      sizerAmp->Add(ampMin);
      sizerAmp->Add(ampGamma);
      mpSizer->Add(sizerAmp);
      mList.Add(ampMax);
      mList.Add(ampMin);
      mList.Add(ampGamma);
      ampMax->SetFormat(_T("%8f"));
      ampMin->SetFormat(_T("%8f"));
      ampGamma->SetFormat(_T("%8f"));
      
      opt=new WXFieldOption(this,-1,&(mpMonteCarloObj->mSaveTrackedData));
      mpSizer->Add(opt,0,wxALIGN_LEFT);
      mList.Add(opt);
      opt->SetToolTip(_T("Saved Tracked values (costs, Chi^2, parameters...)\n\n")
                      _T("This is only useful for Test purposes.\n")
                      _T("Data is saved in the file (Name)-Tracker(-Run#).dat"));

      opt=new WXFieldOption(this,-1,&(mpMonteCarloObj->mXMLAutoSave));
      mpSizer->Add(opt,0,wxALIGN_LEFT);
      mList.Add(opt);
      opt->SetToolTip(_T("Periodically save the best configuration\n\n")
                      _T("Recommended choice is : After Each Run\n")
                      _T("File name is: (name)-(date)-(Run#)-cost.xml\n\n")
                      _T("For Multiple Runs, Note that all choices\n")
                      _T("save the *best* configuration overall, except for\n")
                      _T("'After Each Run', for which the configuration\n")
                      _T("saved are the best for each run."));

      opt=new WXFieldOption(this,-1,&(mpMonteCarloObj->mAutoLSQ));
      mpSizer->Add(opt,0,wxALIGN_LEFT);
      mList.Add(opt);
      opt->SetToolTip(_T("Least squares refinement can be run:\n\n")
                      _T(" - at the end of each run\n")
                      _T(" - perdiodically during the optimization\n\n")
                      _T(" This allows to find the global minimum\n")
                      _T("much faster\n")
                      _T(" Note that if a LSQ refinement is run but does\n")
                      _T("not reach the real global minimum, the returned\n")
                      _T("structure can be very distorted, but this is\n")
                      _T("harmless (restraints will bring back a correct \n")
                      _T("conformation after a few thousand tests)"));
   
   // Number of trials to go
      mpWXFieldNbTrial=new WXFieldPar<long>(this,"Number of trials per run:",-1,&mNbTrial,70);
      mpSizer->Add(mpWXFieldNbTrial);
      mList.Add(mpWXFieldNbTrial);
      mpWXFieldNbTrial->SetFormat(_T("%ld"));
      mpWXFieldNbTrial->SetToolTip(_T("Number of triels per run.\n")
             _T("This number will be updated during the optimization.\n\n")
             _T("Using Multiple Runs:\n")
             _T("  For simple problems (e.g. PbSO4), use 200 000\n")
             _T("  For larger problems (e.g. Cimetidine), use 2 000 000\n")
             _T("  For much larger problems, use 10 000 000\n\n")
             _T("For a single run using Parallel Tempering, use a large number (100 000 000:\n"));
   // Number of cycles (-1=run indefinitely)
      WXFieldPar<long> *pWXFieldNbRun=new WXFieldPar<long>(this,"Number of Runs to perform:",-1,&mNbRun,40);
      mpSizer->Add(pWXFieldNbRun);
      mList.Add(pWXFieldNbRun);
      pWXFieldNbRun->SetFormat(_T("%ld"));
      pWXFieldNbRun->SetToolTip(_T("Number of runs to perform (for Multiple Runs).\n")
                                _T("Use -1 (the default) to run an infinite number of Runs.\n\n")
                                _T("The model will be randomized at the beginning of each run.\n"));
   // Best cost so far
      WXFieldPar<REAL> *pWXFieldBestCost=new WXFieldPar<REAL>(this,"Overall Best Cost:",-1,&(mpMonteCarloObj->GetBestCost()),130);
      mpSizer->Add(pWXFieldBestCost);
      mList.Add(pWXFieldBestCost);
      pWXFieldBestCost->SetFormat(_T("%12.2f"));
      pWXFieldBestCost->SetToolTip(_T("Overall (all runs) Best Cost"));
   this->BottomLayout(0);
   this->CrystUpdate(true);
   VFN_DEBUG_EXIT("WXMonteCarloObj::WXMonteCarloObj()",7)
}

void WXMonteCarloObj::OnRunOptimization(wxCommandEvent & event)
{
   VFN_DEBUG_ENTRY("WXMonteCarloObj::OnRunOptimization()",6)
   WXCrystValidateAllUserInput();
   if(true==this->GetOptimizationObj().IsOptimizing())
   {
      wxMessageDialog dumbUser(this,_T("The optimization is already running !"),_T("Huh ?"),
                               wxOK|wxICON_EXCLAMATION);
      dumbUser.ShowModal();
      VFN_DEBUG_EXIT("WXMonteCarloObj::OnRunOptimization()",6)
      return;
   }
   
   //Fix parameters than really should not be global-optimized
      mpMonteCarloObj->SetParIsFixed(gpRefParTypeUnitCell,true);
      mpMonteCarloObj->SetParIsFixed(gpRefParTypeScattDataScale,true);
      mpMonteCarloObj->SetParIsFixed(gpRefParTypeScattDataProfile,true);
      mpMonteCarloObj->SetParIsFixed(gpRefParTypeScattDataCorrIntAbsorp,true);
      mpMonteCarloObj->SetParIsFixed(gpRefParTypeScattDataCorrIntPolar,true);
      mpMonteCarloObj->SetParIsFixed(gpRefParTypeScattDataCorrIntExtinc,true);
      mpMonteCarloObj->SetParIsFixed(gpRefParTypeScattDataCorrPos,true);
      mpMonteCarloObj->SetParIsFixed(gpRefParTypeScattDataBackground,true);
      mpMonteCarloObj->SetParIsFixed(gpRefParTypeRadiation,true);
      mpMonteCarloObj->UpdateDisplay();
   
   double finalCost=0;
   if(mNbTrial<0)
   {
      mNbTrial = - mNbTrial;
      wxTextEntryDialog costDialog(this,_T("Enter desired cost for the optimization to stop"),
                              _T("Goal Cost"),_T(".20"),wxOK | wxCANCEL);
      if(wxID_OK==costDialog.ShowModal()) costDialog.GetValue().ToDouble(&finalCost);
   }
   if(event.GetId()==ID_GLOBALOPT_MENU_OPT_RUN_MULTIPLE)
      mpGlobalOptimRunThread = new WXGlobalOptimRunThread(this->GetOptimizationObj(),
                                                          mNbTrial,finalCost,mNbRun,true);
   else
      mpGlobalOptimRunThread = new WXGlobalOptimRunThread(this->GetOptimizationObj(),
                                                          mNbTrial,finalCost,mNbRun,false);
   // Tracker window
   if(this->GetOptimizationObj().GetMainTracker().WXGet()==0)
   {
      wxFrame *frame= new wxFrame(this,-1,_T("Tracked data"),
                                  wxDefaultPosition,wxSize(300,200));
      wxWindow* pwxTrackerGraph = this->GetOptimizationObj().GetMainTracker().WXCreate(frame);

      wxSizer *ps=new wxBoxSizer(wxHORIZONTAL);
      ps->Add(pwxTrackerGraph,1,wxEXPAND);
      frame->CreateStatusBar(2);
      frame->SetSizer(ps);
      frame->SetAutoLayout(true);
      frame->Show(true);
   }
   if(mpGlobalOptimRunThread->Create() != wxTHREAD_NO_ERROR) 
      wxLogError(_T("Can't create optimization thread"));
   else mpGlobalOptimRunThread->Run();
   
   VFN_DEBUG_EXIT("WXMonteCarloObj::OnRunOptimization()",6)
}

void WXMonteCarloObj::OnLSQRefine(wxCommandEvent &event)
{
   WXCrystValidateAllUserInput();
   if(true==this->GetOptimizationObj().IsOptimizing())
   {
      wxMessageDialog dumbUser(this,_T("An optimization is already running !"),_T("Huh ?"),
                               wxOK|wxICON_EXCLAMATION);
      dumbUser.ShowModal();
      return;
   }
   char buf[200];
   mpMonteCarloObj->PrepareRefParList();
   mpMonteCarloObj->InitLSQ(true);
   
   sprintf(buf,"LSQ: start");
   REAL cost=mpMonteCarloObj->GetLogLikelihood();
   mpMonteCarloObj->mvSavedParamSet.push_back(make_pair(mpMonteCarloObj->mRefParList.CreateParamSet(buf),cost));
   
   wxProgressDialog dlgProgress(_T("Least Squares refinement"),wxString::Format(_T("Least Squares refinement, cycle #%02d/20, Chi^2=%012.2f"),0,cost),
                                 19,this,wxPD_AUTO_HIDE|wxPD_ELAPSED_TIME|wxPD_CAN_ABORT);
   for(unsigned i=0;i<20;++i)
   {
      try {mpMonteCarloObj->GetLSQObj().Refine(1,true,false);}
      catch(const ObjCrystException &except){};
      mpMonteCarloObj->UpdateDisplay();
      sprintf(buf,"LSQ: cycle #%02d",i);
      cost=mpMonteCarloObj->GetLogLikelihood();
      mpMonteCarloObj->mvSavedParamSet.push_back(make_pair(mpMonteCarloObj->mRefParList.CreateParamSet(buf),cost));
      if(dlgProgress.Update(i,wxString::Format(_T("Least Squares refinement, cycle #%02d/20, Chi^2=%012.2f"),i,cost))==false) return;
   }
}

void WXMonteCarloObj::UpdateDisplayNbTrial()
{
   VFN_DEBUG_MESSAGE("WXMonteCarloObj::UpdateDisplayNbTrial()",5)
   mMutex.Lock();
   mList.CrystUpdate(false);
   mMutex.Unlock();
   wxUpdateUIEvent event(ID_CRYST_UPDATEUI);
   wxPostEvent(this,event);
}

OptimizationObj & WXMonteCarloObj::GetOptimizationObj()
{
   VFN_DEBUG_MESSAGE("WXMonteCarloObj::GetOptimizationObj()",2)
   return *mpMonteCarloObj;
}
const OptimizationObj & WXMonteCarloObj::GetOptimizationObj()const
{
   VFN_DEBUG_MESSAGE("WXMonteCarloObj::GetOptimizationObj() const",2)
   return *mpMonteCarloObj;
}


}// namespace 

