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
static long ID_GLOBALOPT_MENU_SOLUTIONS=             WXCRYST_ID(); 
static long ID_GLOBALOPT_MENU_SOLUTIONS_BROWSE=      WXCRYST_ID(); 
static long ID_BROWSE_WIN=                           WXCRYST_ID(); 

WXOptimizationObj::WXOptimizationObj(wxWindow* parent, OptimizationObj *obj):
WXCrystObj(parent),mpGlobalOptimRunThread(0)
{
   VFN_DEBUG_ENTRY("WXOptimizationObj::WXOptimizationObj(wxWindow*,GlobalOptimObj*,)",7)
   mpWXTitle->SetForegroundColour(wxColour(255,0,0));
   mpWXTitle->SetLabel("Global Optimization:");
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

void WXOptimizationObj::CrystUpdate()
{
   this->WXCrystObj::CrystUpdate();
   if(true==wxThread::IsMain()) this->UpdateUI();
   else
   {
      wxUpdateUIEvent event(ID_CRYST_UPDATEUI);
      wxPostEvent(this,event);
   }
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
   refobj->UpdateUI();
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
   this->UpdateUI();
   VFN_DEBUG_EXIT("WXOptimizationObj::OnUpdateUI()",5)
}

void WXOptimizationObj::UpdateUI()
{
   VFN_DEBUG_ENTRY("WXOptimizationObj::UpdateUI()",5)
   mpWXTitle->SetValue(this->GetOptimizationObj().GetName());
   mpWXTitle->UpdateUI();
   this->WXCrystObj::UpdateUI();
   VFN_DEBUG_EXIT("WXOptimizationObj::UpdateUI()",5)
}

void WXOptimizationObj::OnBrowseParamSet(wxCommandEvent & WXUNUSED(event))
{
   if(this->GetOptimizationObj().IsOptimizing())
   {
      wxMessageDialog dumbUser(this,"Cannot browse during Optimisation !",
                               "Cannot browse during Optimisation!",wxOK|wxICON_EXCLAMATION);
      dumbUser.ShowModal();
      return;
   }
   wxFrame *frame= new wxFrame(this,-1,"Stored Configurations",
                               wxDefaultPosition,wxSize(250,200));
   const long nb=this->GetOptimizationObj().mvSavedParamSet.size();
   wxString *choices = new wxString[nb];
   for(int i=0;i<nb;i++)
   {
      choices[i].sprintf("%d, cost= %f, %s",i,
                         this->GetOptimizationObj().mvSavedParamSet[i].second,
                         this->GetOptimizationObj().mRefParList.GetParamSetName
                           (this->GetOptimizationObj().mvSavedParamSet[i].first).c_str());
      //cout<<choices[i]<<endl;
   }
   wxListBox* wxlist=new wxListBox(frame, ID_BROWSE_WIN, wxDefaultPosition, 
                                   wxDefaultSize, nb, choices,
                                   wxLB_SINGLE|wxLB_NEEDED_SB, wxDefaultValidator,
                                   "listBox");
   wxlist->SetEventHandler(this);
   mClockParamSetWindow.Click();
   frame->Show(true);
}

void WXOptimizationObj::OnSelectParamSet(wxCommandEvent &event)
{
   if(this->GetOptimizationObj().IsOptimizing())
   {
      wxMessageDialog dumbUser(this,"Cannot browse during Optimisation !",
                               "Cannot browse during Optimisation!",wxOK|wxICON_EXCLAMATION);
      dumbUser.ShowModal();
      return;
   }
   const long n=event.GetSelection();
   if(mClockParamSetWindow>this->GetOptimizationObj().mRefParList.GetRefParListClock())
   {
      this->GetOptimizationObj().mRefParList
         .RestoreParamSet(this->GetOptimizationObj().mvSavedParamSet[n].first);
      this->GetOptimizationObj().UpdateDisplay();
      cout <<"Param set #"<<this->GetOptimizationObj().mvSavedParamSet[n].first<<", cost="
           <<this->GetOptimizationObj().mvSavedParamSet[n].second
           <<", now cost="<<this->GetOptimizationObj().GetLogLikelihood()<<endl;
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
   if(mDoMultiple) mpGlobalOptObj->MultiRunOptimize(*mpNbRun,*mpNbTrial,false,mFinalCost);
   else mpGlobalOptObj->Optimize(*mpNbTrial,false,mFinalCost);
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
   EVT_MENU(ID_GLOBALOPT_MENU_OBJECTS_ADDOBJ,        WXOptimizationObj::OnAddRefinedObject)
   EVT_MENU(ID_GLOBALOPT_MENU_OPT_RUN,           WXOptimizationObj::OnRunOptimization)
   EVT_MENU(ID_GLOBALOPT_MENU_OPT_RUN_MULTIPLE,  WXOptimizationObj::OnRunOptimization)
   EVT_MENU(ID_GLOBALOPT_MENU_OPT_STOP,          WXOptimizationObj::OnStopOptimization)
   EVT_MENU(ID_GLOBALOPT_MENU_SOLUTIONS_BROWSE,  WXOptimizationObj::OnBrowseParamSet)
   EVT_UPDATE_UI(ID_CRYST_UPDATEUI,                    WXOptimizationObj::OnUpdateUI)
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
   // Number of trials to go
      mpWXFieldNbTrial=new WXFieldPar<long>(this,"Number of trials per run:",-1,&mNbTrial,70);
      mpSizer->Add(mpWXFieldNbTrial);
      mList.Add(mpWXFieldNbTrial);
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
      pWXFieldNbRun->SetToolTip(_T("Number of runs to perform (for Multiple Runs).\n")
                                _T("Use -1 (the default) to run an infinite number of Runs.\n\n")
                                _T("The model will be randomized at the beginning of each run.\n"));
   this->BottomLayout(0);
   this->CrystUpdate();
   VFN_DEBUG_EXIT("WXMonteCarloObj::WXMonteCarloObj()",7)
}

void WXMonteCarloObj::OnRunOptimization(wxCommandEvent & event)
{
   VFN_DEBUG_ENTRY("WXGeneticAlgorithm::OnRunOptimization()",6)
   WXCrystValidateAllUserInput();
   if(true==this->GetOptimizationObj().IsOptimizing())
   {
      wxMessageDialog dumbUser(this,"The optimization is already running !","Huh ?",
                               wxOK|wxICON_EXCLAMATION);
      dumbUser.ShowModal();
      VFN_DEBUG_EXIT("WXGeneticAlgorithm::OnRunOptimization()",6)
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
      wxTextEntryDialog costDialog(this,"Enter desired cost for the optimization to stop",
                              "Goal Cost",".20",wxOK | wxCANCEL);
      if(wxID_OK==costDialog.ShowModal()) costDialog.GetValue().ToDouble(&finalCost);
   }
   if(event.GetId()==ID_GLOBALOPT_MENU_OPT_RUN_MULTIPLE)
      mpGlobalOptimRunThread = new WXGlobalOptimRunThread(this->GetOptimizationObj(),
                                                          mNbTrial,finalCost,mNbRun,true);
   else
      mpGlobalOptimRunThread = new WXGlobalOptimRunThread(this->GetOptimizationObj(),
                                                          mNbTrial,finalCost,mNbRun,false);
   if(mpGlobalOptimRunThread->Create() != wxTHREAD_NO_ERROR) 
      wxLogError("Can't create optimization thread");
   else mpGlobalOptimRunThread->Run();
   
   VFN_DEBUG_EXIT("WXMonteCarloObj::OnRunOptimization()",6)
}
void WXMonteCarloObj::UpdateDisplayNbTrial()
{
   VFN_DEBUG_MESSAGE("WXMonteCarloObj::UpdateDisplayNbTrial()",5)
   mList.CrystUpdate();
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

