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
/*
*  source file for Global Optimization Objects
*
*/
#include <iomanip>

#include "RefinableObj/GlobalOptimObj.h"
#include "Quirks/VFNStreamFormat.h"
#include "Quirks/VFNDebug.h"
#include "Quirks/Chronometer.h"
#include "ObjCryst/IO.h"
#ifdef __WX__CRYST__
   #include "wxCryst/wxRefinableObj.h"
   #undef GetClassName // Conflict from wxMSW headers ? (cygwin)
#endif

//For some reason, with wxWindows this must be placed after wx headers (Borland c++)
#include <fstream>
#include <sstream>
#include <stdio.h>

namespace ObjCryst
{
void CompareWorlds(const CrystVector_long &idx,const CrystVector_long &swap, const RefinableObj &obj)
{
   const long nb=swap.numElements();
   const CrystVector_REAL *pv0=&(obj.GetParamSet(idx(swap(nb-1))));
   for(long i=0;i<idx.numElements();++i)
   {
      REAL d=0.0;
      const CrystVector_REAL *pv1=&(obj.GetParamSet(idx(swap(i))));
      for(long j=0;j<pv0->numElements();++j) d += ((*pv0)(j)-(*pv1)(j))*((*pv0)(j)-(*pv1)(j));
      cout<<"d("<<i<<")="<<sqrt(d)<<endl;
   }
   cout<<endl;
}
//#################################################################################
//
//       OptimizationObj
//
//#################################################################################
ObjRegistry<OptimizationObj> gOptimizationObjRegistry("List of all Optimization objects");


OptimizationObj::OptimizationObj(const string name):
mName(name),mSaveFileName("GlobalOptim.save"),
mNbTrial(0),mBestCost(-1),
mBestParSavedSetIndex(-1),
mContext(0),
mIsOptimizing(false),mStopAfterCycle(false),
mRefinedObjList("OptimizationObj: "+mName+" RefinableObj registry"),
mRecursiveRefinedObjList("OptimizationObj: "+mName+" recursive RefinableObj registry"),
mLastOptimTime(0)
{
   VFN_DEBUG_ENTRY("OptimizationObj::OptimizationObj()",5)
   // This must be done in a real class to avoid calling a pure virtual method
   // if a graphical representation is automatically called upon registration.
   //  gOptimizationObjRegistry.Register(*this);
   
   // We only copy parameters, so do not delete them !
   mRefParList.SetDeleteRefParInDestructor(false);
   VFN_DEBUG_EXIT("OptimizationObj::OptimizationObj()",5)
}

OptimizationObj::~OptimizationObj()
{
   VFN_DEBUG_ENTRY("OptimizationObj::~OptimizationObj()",5)
   gOptimizationObjRegistry.DeRegister(*this);
   VFN_DEBUG_EXIT("OptimizationObj::~OptimizationObj()",5)
}

void OptimizationObj::RandomizeStartingConfig()
{
   VFN_DEBUG_ENTRY("OptimizationObj::RandomizeStartingConfig()",5)
   this->PrepareRefParList();
   for(int j=0;j<mRefParList.GetNbParNotFixed();j++)
   {
      if(true==mRefParList.GetParNotFixed(j).IsLimited())
      {
         const REAL min=mRefParList.GetParNotFixed(j).GetMin();
         const REAL max=mRefParList.GetParNotFixed(j).GetMax();
         mRefParList.GetParNotFixed(j).MutateTo(min+(max-min)*(rand()/(REAL)RAND_MAX) );
      }
      else if(true==mRefParList.GetParNotFixed(j).IsPeriodic())
             mRefParList.GetParNotFixed(j).
                Mutate(mRefParList.GetParNotFixed(j).GetPeriod()*rand()/(REAL)RAND_MAX);
   }
      //else cout << mRefParList.GetParNotFixed(j).Name() <<" Not limited :-(" <<endl;
   VFN_DEBUG_EXIT("OptimizationObj::RandomizeStartingConfig()",5)
}

void OptimizationObj::FixAllPar()
{
   VFN_DEBUG_ENTRY("OptimizationObj::FixAllPar()",5)
   this->BuildRecursiveRefObjList();
   for(int i=0;i<mRecursiveRefinedObjList.GetNb();i++) 
      mRecursiveRefinedObjList.GetObj(i).FixAllPar();
   VFN_DEBUG_EXIT("OptimizationObj::FixAllPar():End",5)
}
void OptimizationObj::SetParIsFixed(const string& parName,const bool fix)
{
   this->BuildRecursiveRefObjList();
   for(int i=0;i<mRecursiveRefinedObjList.GetNb();i++) 
      mRecursiveRefinedObjList.GetObj(i).SetParIsFixed(parName,fix);
}
void OptimizationObj::SetParIsFixed(const RefParType *type,const bool fix)
{
   this->BuildRecursiveRefObjList();
   for(int i=0;i<mRecursiveRefinedObjList.GetNb();i++) 
      mRecursiveRefinedObjList.GetObj(i).SetParIsFixed(type,fix);
}
   
void OptimizationObj::UnFixAllPar()
{
   this->BuildRecursiveRefObjList();
   for(int i=0;i<mRecursiveRefinedObjList.GetNb();i++) 
      mRecursiveRefinedObjList.GetObj(i).UnFixAllPar();
}
   
void OptimizationObj::SetParIsUsed(const string& parName,const bool use)
{
   this->BuildRecursiveRefObjList();
   for(int i=0;i<mRecursiveRefinedObjList.GetNb();i++) 
      mRecursiveRefinedObjList.GetObj(i).SetParIsUsed(parName,use);
}
void OptimizationObj::SetParIsUsed(const RefParType *type,const bool use)
{
   this->BuildRecursiveRefObjList();
   for(int i=0;i<mRecursiveRefinedObjList.GetNb();i++) 
      mRecursiveRefinedObjList.GetObj(i).SetParIsUsed(type,use);
}
void OptimizationObj::SetLimitsRelative(const string &parName,
                                       const REAL min, const REAL max)
{
   this->BuildRecursiveRefObjList();
   for(int i=0;i<mRecursiveRefinedObjList.GetNb();i++) 
      mRecursiveRefinedObjList.GetObj(i).SetLimitsRelative(parName,min,max);
}
void OptimizationObj::SetLimitsRelative(const RefParType *type,
                                       const REAL min, const REAL max)
{
   this->BuildRecursiveRefObjList();
   for(int i=0;i<mRecursiveRefinedObjList.GetNb();i++) 
      mRecursiveRefinedObjList.GetObj(i).SetLimitsRelative(type,min,max);
}
void OptimizationObj::SetLimitsAbsolute(const string &parName,
                                       const REAL min, const REAL max)
{
   this->BuildRecursiveRefObjList();
   for(int i=0;i<mRecursiveRefinedObjList.GetNb();i++) 
      mRecursiveRefinedObjList.GetObj(i).SetLimitsAbsolute(parName,min,max);
}
void OptimizationObj::SetLimitsAbsolute(const RefParType *type,
                                       const REAL min, const REAL max)
{
   this->BuildRecursiveRefObjList();
   for(int i=0;i<mRecursiveRefinedObjList.GetNb();i++) 
      mRecursiveRefinedObjList.GetObj(i).SetLimitsAbsolute(type,min,max);
}

REAL OptimizationObj::GetLogLikelihood() 
{
   TAU_PROFILE("OptimizationObj::GetLogLikelihood()","void ()",TAU_DEFAULT);
   REAL cost =0.;
   for(int i=0;i<mRecursiveRefinedObjList.GetNb();i++)
   {
      const REAL tmp=mRecursiveRefinedObjList.GetObj(i).GetLogLikelihood();
      if(tmp!=0.)
      {
         LogLikelihoodStats* st=&((mvContextObjStats[mContext])
                                    [&(mRecursiveRefinedObjList.GetObj(i))]);
         st->mTotalLogLikelihood += tmp;
         st->mTotalLogLikelihoodDeltaSq +=
            (tmp-st->mLastLogLikelihood)*(tmp-st->mLastLogLikelihood);
         st->mLastLogLikelihood=tmp;
      }
      cost += mvObjWeight[&(mRecursiveRefinedObjList.GetObj(i))].mWeight * tmp;
   }
   return cost;
}
void OptimizationObj::StopAfterCycle() 
{
   VFN_DEBUG_MESSAGE("OptimizationObj::StopAfterCycle()",5)
   if(mIsOptimizing) mStopAfterCycle=true;
}

void OptimizationObj::DisplayReport() 
{
   //:TODO: ask all objects to print their own report ?
}

void OptimizationObj::AddRefinableObj(RefinableObj &obj)
{
   VFN_DEBUG_MESSAGE("OptimizationObj::AddRefinableObj():"<<obj.GetName(),5)
   //in case some object has been modified, to avoid rebuilding the entire list
      this->BuildRecursiveRefObjList();
      
   mRefinedObjList.Register(obj);
   RefObjRegisterRecursive(obj,mRecursiveRefinedObjList);
   #ifdef __WX__CRYST__
   if(0!=this->WXGet()) this->WXGet()->AddRefinedObject(obj);
   #endif
}

RefinableObj& OptimizationObj::GetFullRefinableObj(const bool rebuild)
{
   if(rebuild) this->PrepareRefParList();
   return mRefParList;
}

const string& OptimizationObj::GetName()const { return mName;}
void OptimizationObj::SetName(const string& name) {mName=name;}

const string OptimizationObj::GetClassName()const { return "OptimizationObj";}

void OptimizationObj::Print()const {this->XMLOutput(cout);}

void OptimizationObj::RestoreBestConfiguration()
{
   //:TODO: check list of refinableObj has not changed, and the list of
   // RefPar has not changed in all sub-objects.
   if(mBestParSavedSetIndex>0) mRefParList.RestoreParamSet(mBestParSavedSetIndex);
}

bool OptimizationObj::IsOptimizing()const{return mIsOptimizing;}

void OptimizationObj::TagNewBestConfig()
{
   for(int i=0;i<mRecursiveRefinedObjList.GetNb();i++)
      mRecursiveRefinedObjList.GetObj(i).TagNewBestConfig();
   this->UpdateDisplay();
}

REAL OptimizationObj::GetLastOptimElapsedTime()const
{
   return mLastOptimTime;
}

void OptimizationObj::PrepareRefParList()
{
   VFN_DEBUG_ENTRY("OptimizationObj::PrepareRefParList()",6)
   
   this->BuildRecursiveRefObjList();
   // As any parameter been added in the recursive list of objects ?
   // or has any object been added/removed ?
      RefinableObjClock clock;
      GetRefParListClockRecursive(mRecursiveRefinedObjList,clock);
   if(  (clock>mRefParList.GetRefParListClock())
      ||(mRecursiveRefinedObjList.GetRegistryClock()>mRefParList.GetRefParListClock()) )
   {
      VFN_DEBUG_MESSAGE("OptimizationObj::PrepareRefParList():Rebuild list",6)
      mRefParList.ResetParList();
      mRefParList.EraseAllParamSet();
      for(int i=0;i<mRecursiveRefinedObjList.GetNb();i++)
         mRefParList.AddPar(mRecursiveRefinedObjList.GetObj(i));
      mvSavedParamSet.clear();
      mBestParSavedSetIndex=mRefParList.CreateParamSet("Best Configuration");
      mvSavedParamSet.push_back(make_pair(mBestParSavedSetIndex,mBestCost));
   }
   // Prepare for refinement, ie get the list of not fixed parameters,
   // and prepare the objects...
   mRefParList.PrepareForRefinement();
   for(int i=0;i<mRecursiveRefinedObjList.GetNb();i++)
      mRecursiveRefinedObjList.GetObj(i).PrepareForRefinement();
   VFN_DEBUG_EXIT("OptimizationObj::PrepareRefParList()",6)
}

void OptimizationObj::InitRandomSeedFromTime()const
{
   VFN_DEBUG_MESSAGE("OptimizationObj::InitRandomSeedFromTime()",3)
   time_t junk;
   time(&junk);
   tm *tmp=localtime(&junk);
   srand((unsigned)( (*tmp).tm_sec+60* (*tmp).tm_min));
   //for(int i=0;i<20;i++) cout << rand() <<endl;
}
void OptimizationObj::InitOptions()
{
   VFN_DEBUG_MESSAGE("OptimizationObj::InitOptions()",5)
   static string xmlAutoSaveName;
   static string xmlAutoSaveChoices[5];
   
   static bool needInitNames=true;
   if(true==needInitNames)
   {
      xmlAutoSaveName="Save Best Config Regularly";
      xmlAutoSaveChoices[0]="No";
      xmlAutoSaveChoices[1]="Every day";
      xmlAutoSaveChoices[2]="Every hour";
      xmlAutoSaveChoices[3]="Every 10mn";
      xmlAutoSaveChoices[4]="Every new best config (a lot ! Not Recommended !)";
      
      needInitNames=false;//Only once for the class
   }
   mXMLAutoSave.Init(5,&xmlAutoSaveName,xmlAutoSaveChoices);
   VFN_DEBUG_MESSAGE("OptimizationObj::InitOptions():End",5)
}

void OptimizationObj::UpdateDisplay()
{
   #ifdef __WX__CRYST__
   if(0!=this->WXGet()) this->WXGet()->CrystUpdate();
   #endif
   for(int i=0;i<mRefinedObjList.GetNb();i++) 
      mRefinedObjList.GetObj(i).UpdateDisplay();
}
void OptimizationObj::BuildRecursiveRefObjList()
{
   // First check if anything has changed (ie if a sub-object has been
   // added or removed in the recursive refinable object list)
   RefinableObjClock clock;
   GetSubRefObjListClockRecursive(mRefinedObjList,clock);
   if(clock>mRecursiveRefinedObjList.GetRegistryClock())
   {
      VFN_DEBUG_ENTRY("OptimizationObj::BuildRecursiveRefObjList()",5)
      mRecursiveRefinedObjList.DeRegisterAll();
      for(int i=0;i<mRefinedObjList.GetNb();i++)
         RefObjRegisterRecursive(mRefinedObjList.GetObj(i),mRecursiveRefinedObjList);
      VFN_DEBUG_EXIT("OptimizationObj::BuildRecursiveRefObjList()",5)
   }
}

//#################################################################################
//
//       MonteCarloObj
//
//#################################################################################
MonteCarloObj::MonteCarloObj(const string name):
OptimizationObj(name),
mCurrentCost(-1),
mHistorySaveFileName("GlobalOptim_history.out"),
mTemperatureMax(1e6),mTemperatureMin(.001),mTemperatureGamma(1.0),
mMutationAmplitudeMax(8.),mMutationAmplitudeMin(.125),mMutationAmplitudeGamma(1.0),
mNbTrialRetry(0),mMinCostRetry(0)
#ifdef __WX__CRYST__
,mpWXCrystObj(0)
#endif
{
   VFN_DEBUG_ENTRY("MonteCarloObj::MonteCarloObj()",5)
   this->InitOptions();
   mGlobalOptimType.SetChoice(GLOBAL_OPTIM_PARALLEL_TEMPERING);
   mAnnealingScheduleTemp.SetChoice(ANNEALING_SMART);
   mAnnealingScheduleMutation.SetChoice(ANNEALING_EXPONENTIAL);
   gOptimizationObjRegistry.Register(*this);
   this->InitRandomSeedFromTime();
   VFN_DEBUG_EXIT("MonteCarloObj::MonteCarloObj()",5)
}

MonteCarloObj::MonteCarloObj(const bool internalUseOnly):
OptimizationObj(""),
mCurrentCost(-1),
mHistorySaveFileName("GlobalOptim_history.out"),
mTemperatureMax(.03),mTemperatureMin(.003),
mMutationAmplitudeMax(16.),mMutationAmplitudeMin(.125),
mNbTrialRetry(0),mMinCostRetry(0)
#ifdef __WX__CRYST__
,mpWXCrystObj(0)
#endif
{
   VFN_DEBUG_ENTRY("MonteCarloObj::MonteCarloObj(bool)",5)
   this->InitOptions();
   mGlobalOptimType.SetChoice(GLOBAL_OPTIM_PARALLEL_TEMPERING);
   mAnnealingScheduleTemp.SetChoice(ANNEALING_SMART);
   mAnnealingScheduleMutation.SetChoice(ANNEALING_EXPONENTIAL);
   if(false==internalUseOnly) gOptimizationObjRegistry.Register(*this);
   this->InitRandomSeedFromTime();
   VFN_DEBUG_EXIT("MonteCarloObj::MonteCarloObj(bool)",5)
}

MonteCarloObj::~MonteCarloObj()
{
   VFN_DEBUG_ENTRY("MonteCarloObj::~MonteCarloObj()",5)
   gOptimizationObjRegistry.DeRegister(*this);
   VFN_DEBUG_EXIT ("MonteCarloObj::~MonteCarloObj()",5)
}
void MonteCarloObj::SetAlgorithmSimulAnnealing(const AnnealingSchedule scheduleTemp,
                           const REAL tMax, const REAL tMin,
                           const AnnealingSchedule scheduleMutation,
                           const REAL mutMax, const REAL mutMin,
                           const long nbTrialRetry,const REAL minCostRetry)
{
   VFN_DEBUG_MESSAGE("MonteCarloObj::SetAlgorithmSimulAnnealing()",5)
   mGlobalOptimType.SetChoice(GLOBAL_OPTIM_SIMULATED_ANNEALING);
   mTemperatureMax=tMax;
   mTemperatureMin=tMin;
   mAnnealingScheduleTemp.SetChoice(scheduleTemp);


   mMutationAmplitudeMax=mutMax;
   mMutationAmplitudeMin=mutMin;
   mAnnealingScheduleMutation.SetChoice(scheduleMutation);
   mNbTrialRetry=nbTrialRetry;
   mMinCostRetry=minCostRetry;
   VFN_DEBUG_MESSAGE("MonteCarloObj::SetAlgorithmSimulAnnealing():End",3)
}

void MonteCarloObj::SetAlgorithmParallTempering(const AnnealingSchedule scheduleTemp,
                                 const REAL tMax, const REAL tMin,
                                 const AnnealingSchedule scheduleMutation,
                                 const REAL mutMax, const REAL mutMin)
{
   VFN_DEBUG_MESSAGE("MonteCarloObj::SetAlgorithmParallTempering()",5)
   mGlobalOptimType.SetChoice(GLOBAL_OPTIM_PARALLEL_TEMPERING);
   mTemperatureMax=tMax;
   mTemperatureMin=tMin;
   mAnnealingScheduleTemp.SetChoice(scheduleTemp);

   mMutationAmplitudeMax=mutMax;
   mMutationAmplitudeMin=mutMin;
   mAnnealingScheduleMutation.SetChoice(scheduleMutation);
   //mNbTrialRetry=nbTrialRetry;
   //mMinCostRetry=minCostRetry;
   VFN_DEBUG_MESSAGE("MonteCarloObj::SetAlgorithmParallTempering():End",3)
}
void MonteCarloObj::Optimize(long &nbStep,const bool silent,const REAL finalcost,
                             const REAL maxTime)
{
   //:TODO: Other algorithms !
   TAU_PROFILE("MonteCarloObj::Optimize()","void (long)",TAU_DEFAULT);
   VFN_DEBUG_ENTRY("MonteCarloObj::Optimize()",5)
   for(int i=0;i<mRefinedObjList.GetNb();i++) mRefinedObjList.GetObj(i).BeginOptimization(true);
   this->PrepareRefParList();
   if(!silent) mRefParList.Print();
   mIsOptimizing=true;
   if(mTemperatureGamma<0.1) mTemperatureGamma= 0.1;
   if(mTemperatureGamma>10.0)mTemperatureGamma=10.0;
   if(mMutationAmplitudeGamma<0.1) mMutationAmplitudeGamma= 0.1;
   if(mMutationAmplitudeGamma>10.0)mMutationAmplitudeGamma=10.0;
   // prepare all objects
   this->TagNewBestConfig();
   mCurrentCost=this->GetLogLikelihood();
   mBestCost=mCurrentCost;
   mvObjWeight.clear();
   switch(mGlobalOptimType.GetChoice())
   {
      case GLOBAL_OPTIM_SIMULATED_ANNEALING:
      {
         this->RunSimulatedAnnealing(nbStep,silent,finalcost,maxTime);
         break;
      }//case GLOBAL_OPTIM_SIMULATED_ANNEALING
      case GLOBAL_OPTIM_PARALLEL_TEMPERING:
      {
         this->RunParallelTempering(nbStep,silent,finalcost,maxTime);
         break;
      }//case GLOBAL_OPTIM_PARALLEL_TEMPERING
      case GLOBAL_OPTIM_GENETIC: //:TODO:
      {
      }//case GLOBAL_OPTIM_GENETIC
   }
   mIsOptimizing=false;
   mStopAfterCycle=false;
      
   for(int i=0;i<mRefinedObjList.GetNb();i++) mRefinedObjList.GetObj(i).EndOptimization();
   VFN_DEBUG_EXIT("MonteCarloObj::Optimize()",5)
}
void MonteCarloObj::MultiRunOptimize(long &nbCycle,long &nbStep,const bool silent,
                                     const REAL finalcost,const REAL maxTime)
{
   //:TODO: Other algorithms !
   TAU_PROFILE("MonteCarloObj::MultiRunOptimize()","void (long)",TAU_DEFAULT);
   VFN_DEBUG_ENTRY("MonteCarloObj::MultiRunOptimize()",5)
   //Keep a copy of the total number of steps, and decrement nbStep
   const long nbStep0=nbStep;
   for(int i=0;i<mRefinedObjList.GetNb();i++) mRefinedObjList.GetObj(i).BeginOptimization(true);
   this->PrepareRefParList();
   if(!silent) mRefParList.Print();
   mIsOptimizing=true;
   if(mTemperatureGamma<0.1) mTemperatureGamma= 0.1;
   if(mTemperatureGamma>10.0)mTemperatureGamma=10.0;
   if(mMutationAmplitudeGamma<0.1) mMutationAmplitudeGamma= 0.1;
   if(mMutationAmplitudeGamma>10.0)mMutationAmplitudeGamma=10.0;
   // prepare all objects
   mCurrentCost=this->GetLogLikelihood();
   mBestCost=mCurrentCost;
   this->TagNewBestConfig();
   mvObjWeight.clear();
   while(nbCycle!=0)
   {
      if(!silent) cout <<"MonteCarloObj::MultiRunOptimize: Starting Run#"<<abs(nbCycle)<<endl;
      nbStep=nbStep0;
      for(int i=0;i<mRefinedObjList.GetNb();i++) mRefinedObjList.GetObj(i).RandomizeConfiguration();
      switch(mGlobalOptimType.GetChoice())
      {
         case GLOBAL_OPTIM_SIMULATED_ANNEALING:
         {
            this->RunSimulatedAnnealing(nbStep,silent,finalcost,maxTime);
            break;
         }
         case GLOBAL_OPTIM_PARALLEL_TEMPERING:
         {
            this->RunParallelTempering(nbStep,silent,finalcost,maxTime);
            break;
         }
         case GLOBAL_OPTIM_GENETIC: //:TODO:
         {
         }
      }
      nbStep=nbStep0;
      this->UpdateDisplay();
      stringstream s;
      s<<"Run #"<<abs(nbCycle);
      mvSavedParamSet.push_back(make_pair(mRefParList.CreateParamSet(s.str()),mCurrentCost));
      if(!silent) cout <<"MonteCarloObj::MultiRunOptimize: Finished Run#"
                       <<abs(nbCycle)<<", Run Best Cost:"<<mCurrentCost
                       <<", Overall Best Cost:"<<mBestCost<<endl;
      if((mStopAfterCycle)||(mBestCost<finalcost)) break;
      nbCycle--;
   }
   mIsOptimizing=false;
   mStopAfterCycle=false;

   mRefParList.RestoreParamSet(mBestParSavedSetIndex);
   this->UpdateDisplay();
      
   for(int i=0;i<mRefinedObjList.GetNb();i++) mRefinedObjList.GetObj(i).EndOptimization();
   VFN_DEBUG_EXIT("MonteCarloObj::MultiRunOptimize()",5)
}

void MonteCarloObj::RunSimulatedAnnealing(long &nbStep,const bool silent,
                                          const REAL finalcost,const REAL maxTime)
{
   //Keep a copy of the total number of steps, and decrement nbStep
   const long nbSteps=nbStep;
   unsigned int accept;// 1 if last trial was accepted? 2 if new best config ? else 0
   mNbTrial=0;
   // If using history, open file
      ofstream outHistory;
      if(mSaveDetailledHistory.GetChoice()>0)
      {
         outHistory.open(mHistorySaveFileName.c_str());
         outHistory << "Trial World origWorld Accept OverallCost ";
         for(long j=0;j<mRefParList.GetNbParNotFixed();j++)
            outHistory << mRefParList.GetParNotFixed(j).GetName() << " ";
         outHistory <<endl;
      }
   // time (in seconds) when last autoSave was made (if enabled)
      unsigned long secondsWhenAutoSave=0;

   if(!silent) cout << "Starting Simulated Annealing Optimization for"<<nbSteps<<" trials"<<endl;
   if(!silent) this->DisplayReport();
   REAL runBestCost;
   mCurrentCost=this->GetLogLikelihood();
   runBestCost=mCurrentCost;
   const long lastParSavedSetIndex=mRefParList.CreateParamSet("MonteCarloObj:Last parameters");
   const long runBestIndex=mRefParList.CreateParamSet("best parameters for current run");
   //Report each ... cycles
      const int nbTryReport=3000;
   // Keep record of the number of accepted moves
      long nbAcceptedMoves=0;//since last report
      long nbAcceptedMovesTemp=0;//since last temperature/mutation rate change
   // Number of tries since best configuration found
      long nbTriesSinceBest=0;
   // Change temperature (and mutation) every...
      const int nbTryPerTemp=300;

   mTemperature=sqrt(mTemperatureMin*mTemperatureMax);
   mMutationAmplitude=sqrt(mMutationAmplitudeMin*mMutationAmplitudeMax);

   Chronometer chrono;
   chrono.start();
   for(mNbTrial=1;mNbTrial<=nbSteps;)
   {
      if((mNbTrial % nbTryPerTemp) == 1)
      {
         VFN_DEBUG_MESSAGE("-> Updating temperature and mutation amplitude.",3)
         // Temperature & displacements amplitude
         switch(mAnnealingScheduleTemp.GetChoice())
         {
            case ANNEALING_BOLTZMANN:
               mTemperature=
                  mTemperatureMin*log((REAL)nbSteps)/log((REAL)(mNbTrial+1));break;
            case ANNEALING_CAUCHY:
               mTemperature=mTemperatureMin*nbSteps/mNbTrial;break;
            //case ANNEALING_QUENCHING:
            case ANNEALING_EXPONENTIAL:
               mTemperature=mTemperatureMax
                              *pow(mTemperatureMin/mTemperatureMax,
                                    mNbTrial/(REAL)nbSteps);break;
            case ANNEALING_GAMMA:
               mTemperature=mTemperatureMax+(mTemperatureMin-mTemperatureMax)
                              *pow(mNbTrial/(REAL)nbSteps,mTemperatureGamma);break;
            case ANNEALING_SMART:
            {
               if((nbAcceptedMovesTemp/(REAL)nbTryPerTemp)>0.30)
                  mTemperature/=1.5;
               if((nbAcceptedMovesTemp/(REAL)nbTryPerTemp)<0.10)
                  mTemperature*=1.5;
               if(mTemperature>mTemperatureMax) mTemperature=mTemperatureMax;
               if(mTemperature<mTemperatureMin) mTemperature=mTemperatureMin;
               nbAcceptedMovesTemp=0;
               break;
            }
            default: mTemperature=mTemperatureMin;break;
         }
         switch(mAnnealingScheduleMutation.GetChoice())
         {
            case ANNEALING_BOLTZMANN:
               mMutationAmplitude=
                  mMutationAmplitudeMin*log((REAL)nbSteps)/log((REAL)(mNbTrial+1));
               break;
            case ANNEALING_CAUCHY:
               mMutationAmplitude=mMutationAmplitudeMin*nbSteps/mNbTrial;break;
            //case ANNEALING_QUENCHING:
            case ANNEALING_EXPONENTIAL:
               mMutationAmplitude=mMutationAmplitudeMax
                              *pow(mMutationAmplitudeMin/mMutationAmplitudeMax,
                                    mNbTrial/(REAL)nbSteps);break;
            case ANNEALING_GAMMA:
               mMutationAmplitude=mMutationAmplitudeMax+(mMutationAmplitudeMin-mMutationAmplitudeMax)
                              *pow(mNbTrial/(REAL)nbSteps,mMutationAmplitudeGamma);break;
            case ANNEALING_SMART:
               if((nbAcceptedMovesTemp/(REAL)nbTryPerTemp)>0.3) mMutationAmplitude*=2.;
               if((nbAcceptedMovesTemp/(REAL)nbTryPerTemp)<0.1) mMutationAmplitude/=2.;
               if(mMutationAmplitude>mMutationAmplitudeMax) 
                  mMutationAmplitude=mMutationAmplitudeMax;
               if(mMutationAmplitude<mMutationAmplitudeMin) 
                  mMutationAmplitude=mMutationAmplitudeMax;
               nbAcceptedMovesTemp=0;
               break;
            default: mMutationAmplitude=mMutationAmplitudeMin;break;
         }
      }

      this->NewConfiguration();
      accept=0;
      REAL cost=this->GetLogLikelihood();
      if(cost<mCurrentCost)
      {
         accept=1;
         mCurrentCost=cost;
         mRefParList.SaveParamSet(lastParSavedSetIndex);
         if(mCurrentCost<runBestCost)
         {
            accept=2;
            runBestCost=mCurrentCost;
            this->TagNewBestConfig();
            mRefParList.SaveParamSet(runBestIndex);
            if(runBestCost<mBestCost)
            {
               mBestCost=mCurrentCost;
               mRefParList.SaveParamSet(mBestParSavedSetIndex);
               if(!silent) cout << "Trial :" << mNbTrial 
                             << " Temp="<< mTemperature
                             << " Mutation Ampl.: "<<mMutationAmplitude
                             << " NEW OVERALL Best Cost="<<runBestCost<< endl;
            }
            else if(!silent) cout << "Trial :" << mNbTrial 
                             << " Temp="<< mTemperature
                             << " Mutation Ampl.: "<<mMutationAmplitude
                             << " NEW Run Best Cost="<<runBestCost<< endl;
            nbTriesSinceBest=0;
         }
         nbAcceptedMoves++;
         nbAcceptedMovesTemp++;
      }
      else
      {
         if( log((rand()+1)/(REAL)RAND_MAX) < (-(cost-mCurrentCost)/mTemperature) )
         {
            accept=1;
            mCurrentCost=cost;
            mRefParList.SaveParamSet(lastParSavedSetIndex);
            nbAcceptedMoves++;
            nbAcceptedMovesTemp++;
         }
      }
      switch(mSaveDetailledHistory.GetChoice())
      {
         case 0:break;
         case 1: 
         {
            if(accept==2)
            {
               outHistory << mNbTrial <<" 0 " <<accept<<" "
                          <<this->GetLogLikelihood()<<" ";
               for(long i=0;i<mRefParList.GetNbParNotFixed();i++) 
                  outHistory << " "<< mRefParList.GetParNotFixed(i).GetHumanValue() ;
               outHistory <<endl;
            }
            break;
         }
         case 2: 
         {
            if(accept>0)
            {
               outHistory << mNbTrial <<" 0 "<<accept<<" "
                          <<this->GetLogLikelihood()<<" ";
               for(long i=0;i<mRefParList.GetNbParNotFixed();i++) 
                  outHistory << " "<< mRefParList.GetParNotFixed(i).GetHumanValue() ;
               outHistory <<endl;
            }
            break;
         }
         case 3: 
         {
            outHistory << mNbTrial <<" 0 "<<accept<<" "
                       <<this->GetLogLikelihood()<<" ";
            for(long i=0;i<mRefParList.GetNbParNotFixed();i++) 
               outHistory << " "<< mRefParList.GetParNotFixed(i).GetHumanValue() ;
            outHistory <<endl;
            break;
         }
      }
      if(accept==0) mRefParList.RestoreParamSet(lastParSavedSetIndex);

      if( (mNbTrial % nbTryReport) == 0)
      {
         if(!silent) cout <<"Trial :" << mNbTrial << " Temp="<< mTemperature;
         if(!silent) cout <<" Mutation Ampl.: " <<mMutationAmplitude<< " Best Cost=" << runBestCost 
                          <<" Current Cost=" << mCurrentCost 
                          <<" Accepting "<<(int)((REAL)nbAcceptedMoves/nbTryReport*100)
                          <<"% moves" << endl;
         nbAcceptedMoves=0;
         #ifdef __WX__CRYST__
         if(0!=mpWXCrystObj) mpWXCrystObj->UpdateDisplayNbTrial();
         #endif
      }
      mNbTrial++;nbStep--;

      if((runBestCost<finalcost) || mStopAfterCycle ||( (maxTime>0)&&(chrono.seconds()>maxTime))) 
      {
         if(!silent) cout << endl <<endl << "Refinement Stopped."<<endl;
         break;
      }
      nbTriesSinceBest++;
      if(  ((mXMLAutoSave.GetChoice()==1)&&((chrono.seconds()-secondsWhenAutoSave)>86400))
         ||((mXMLAutoSave.GetChoice()==2)&&((chrono.seconds()-secondsWhenAutoSave)>3600))
         ||((mXMLAutoSave.GetChoice()==3)&&((chrono.seconds()-secondsWhenAutoSave)> 600))
         ||((mXMLAutoSave.GetChoice()==4)&&(accept==2)) )
      {
         secondsWhenAutoSave=(unsigned long)chrono.seconds();
         string saveFileName=this->GetName();
         time_t date=time(0);
         char strDate[40];
         strftime(strDate,sizeof(strDate),"%Y-%m-%d_%H-%M-%S",localtime(&date));//%Y-%m-%dT%H:%M:%S%Z
         char costAsChar[30];
         if(accept!=2) mRefParList.RestoreParamSet(mBestParSavedSetIndex);
         sprintf(costAsChar,"-Cost-%f",this->GetLogLikelihood());
         saveFileName=saveFileName+(string)strDate+(string)costAsChar+(string)".xml";
         XMLCrystFileSaveGlobal(saveFileName);
         if(accept!=2) mRefParList.RestoreParamSet(lastParSavedSetIndex);
      }
   }
   mLastOptimTime=chrono.seconds();
   //Restore Best values
   mRefParList.RestoreParamSet(runBestIndex);
   mRefParList.ClearParamSet(lastParSavedSetIndex);
   mCurrentCost=this->GetLogLikelihood();
   if(mSaveDetailledHistory.GetChoice()>0) outHistory.close();
   if(!silent) this->DisplayReport();
   if(!silent) chrono.print();
}

void MonteCarloObj::RunParallelTempering(long &nbStep,const bool silent,
                                         const REAL finalcost,const REAL maxTime)
{
   //Keep a copy of the total number of steps, and decrement nbStep
   const long nbSteps=nbStep;
   unsigned int accept;// 1 if last trial was accepted? 2 if new best config ? else 0
   mNbTrial=0;
   // If using history, open file
      ofstream outHistory;
      if(mSaveDetailledHistory.GetChoice()>0)
      {
         outHistory.open(mHistorySaveFileName.c_str());
         outHistory << "Trial World origWorld Accept OverallCost ";
         for(long j=0;j<mRefParList.GetNbParNotFixed();j++)
            outHistory << mRefParList.GetParNotFixed(j).GetName() << " ";
         outHistory <<endl;
      }
   // time (in seconds) when last autoSave was made (if enabled)
      unsigned long secondsWhenAutoSave=0;

   if(!silent) cout << "Starting Parallel Tempering Optimization"<<endl;
   //Total number of parallel refinements,each is a 'World'. The most stable
   // world must be i=nbWorld-1, and the most changing World (high mutation,
   // high temperature) is i=0.
      const long nbWorld=30;
      CrystVector_long worldSwapIndex(nbWorld);
      for(int i=0;i<nbWorld;++i) worldSwapIndex(i)=i;
   // Number of successive trials for each World. At the end of these trials
   // a swap is tried with the upper World (eg i-1). This number effectvely sets
   // the rate of swapping.
      const int nbTryPerWorld=10;
   // Initialize the costs
      mCurrentCost=this->GetLogLikelihood();
      REAL runBestCost=mCurrentCost;
      CrystVector_REAL currentCost(nbWorld);
      currentCost=mCurrentCost;
   // Init the different temperatures
      CrystVector_REAL simAnnealTemp(nbWorld);
      for(int i=0;i<nbWorld;i++)
      {
         switch(mAnnealingScheduleTemp.GetChoice())
         {
            case ANNEALING_BOLTZMANN:
               simAnnealTemp(i)=
                  mTemperatureMin*log((REAL)nbWorld)/log((REAL)(i+2));break;
            case ANNEALING_CAUCHY:
               simAnnealTemp(i)=mTemperatureMin*nbWorld/(i+1);break;
            //case ANNEALING_QUENCHING:
            case ANNEALING_EXPONENTIAL:
               simAnnealTemp(i)=mTemperatureMax
                              *pow(mTemperatureMin/mTemperatureMax,
                                    i/(REAL)(nbWorld-1));break;
            case ANNEALING_GAMMA:
               simAnnealTemp(i)=mTemperatureMax+(mTemperatureMin-mTemperatureMax)
                              *pow(i/(REAL)(nbWorld-1),mTemperatureGamma);break;
            case ANNEALING_SMART:
               simAnnealTemp(i)=mCurrentCost/(100.+(REAL)i/(REAL)nbWorld*900.);break;
            default:
               simAnnealTemp(i)=mCurrentCost/(100.+(REAL)i/(REAL)nbWorld*900.);break;
         }
      }
   //Init the different mutation rate parameters
      CrystVector_REAL mutationAmplitude(nbWorld);
      for(int i=0;i<nbWorld;i++)
      {
         switch(mAnnealingScheduleMutation.GetChoice())
         {
            case ANNEALING_BOLTZMANN:
               mutationAmplitude(i)=
                  mMutationAmplitudeMin*log((REAL)(nbWorld-1))/log((REAL)(i+2));
               break;
            case ANNEALING_CAUCHY:
               mutationAmplitude(i)=mMutationAmplitudeMin*(REAL)(nbWorld-1)/(i+1);break;
            //case ANNEALING_QUENCHING:
            case ANNEALING_EXPONENTIAL:
               mutationAmplitude(i)=mMutationAmplitudeMax
                              *pow(mMutationAmplitudeMin/mMutationAmplitudeMax,
                                    i/(REAL)(nbWorld-1));break;
            case ANNEALING_GAMMA:
               mutationAmplitude(i)=mMutationAmplitudeMax+(mMutationAmplitudeMin-mMutationAmplitudeMax)
                              *pow(i/(REAL)(nbWorld-1),mMutationAmplitudeGamma);break;
            case ANNEALING_SMART:
               mutationAmplitude(i)=sqrt(mMutationAmplitudeMin*mMutationAmplitudeMax);break;
            default:
               mutationAmplitude(i)=sqrt(mMutationAmplitudeMin*mMutationAmplitudeMax);break;
         }
      }
   // Init the parameter sets for each World
   // All Worlds start from the same (current) configuration.
      CrystVector_long worldCurrentSetIndex(nbWorld);
      for(int i=0;i<nbWorld;i++)
      {
         if(i!=nbWorld)
            for(int j=0;j<mRecursiveRefinedObjList.GetNb();j++)
               mRecursiveRefinedObjList.GetObj(j).RandomizeConfiguration();
         worldCurrentSetIndex(i)=mRefParList.CreateParamSet();
      }
      //mNbTrial=nbSteps;;
      const long lastParSavedSetIndex=mRefParList.CreateParamSet("MonteCarloObj:Last parameters");
   const long runBestIndex=mRefParList.CreateParamSet("best parameters for current run");
      CrystVector_REAL swapPar;
   //Keep track of how many trials are accepted for each World
      CrystVector_long worldNbAcceptedMoves(nbWorld);
      worldNbAcceptedMoves=0;
   //Do a report each... And check if mutation rate is OK (for annealing_smart)s
      const int nbTrialsReport=3000;
   // TEST : allow GENETIC mating of configurations
      //Get gene groups list :TODO: check for missing groups
         CrystVector_uint refParGeneGroupIndex(mRefParList.GetNbPar());
         unsigned int first=1;
         for(int i=0;i<mRecursiveRefinedObjList.GetNb();i++) 
            mRecursiveRefinedObjList.GetObj(i).GetGeneGroup(mRefParList,refParGeneGroupIndex,first);
         if(!silent) 
            for(int i=0;i<mRefParList.GetNbPar();i++)
            {
               cout << "Gene Group:"<<refParGeneGroupIndex(i)<<" :";
               mRefParList.GetPar(i).Print();
            }
      // number of gene groups
      // to select which gene groups are exchanged in the mating
         //const unsigned int nbGeneGroup=refParGeneGroupIndex.max();
         //CrystVector_int crossoverGroupIndex(nbGeneGroup);
         //const long parSetOffspringA=mRefParList.CreateParamSet("Offspring A");
         //const long parSetOffspringB=mRefParList.CreateParamSet("Offspring B");
   // record the statistical distribution n=f(cost function) for each World
      //CrystMatrix_REAL trialsDensity(100,nbWorld+1);
      //trialsDensity=0;
      //for(int i=0;i<100;i++) trialsDensity(i,0)=i/(float)100;
   //Do the refinement
   bool makeReport=false;
   Chronometer chrono;
   chrono.start();
   for(;mNbTrial<nbSteps;)
   {
      for(int i=0;i<nbWorld;i++)
      {
         mContext=i;
         //mRefParList.RestoreParamSet(worldCurrentSetIndex(i));
         mMutationAmplitude=mutationAmplitude(i);
         mTemperature=simAnnealTemp(i);
         for(int j=0;j<nbTryPerWorld;j++)
         {
            //mRefParList.SaveParamSet(lastParSavedSetIndex);
            mRefParList.RestoreParamSet(worldCurrentSetIndex(i));
            this->NewConfiguration();
            accept=0;
            REAL cost=this->GetLogLikelihood();
            //trialsDensity((long)(cost*100.),i+1)+=1;
            if(cost<currentCost(i))
            {
               accept=1;
               currentCost(i)=cost;
               mRefParList.SaveParamSet(worldCurrentSetIndex(i));
               if(cost<runBestCost)
               {
                  accept=2;
                  runBestCost=currentCost(i);
                  this->TagNewBestConfig();
                  mRefParList.SaveParamSet(runBestIndex);
                  if(runBestCost<mBestCost)
                  {
                     mBestCost=currentCost(i);
                     mRefParList.SaveParamSet(mBestParSavedSetIndex);
                     if(!silent) cout << "->Trial :" << mNbTrial 
                                   << " World="<< worldSwapIndex(i)
                                   << " Temp="<< mTemperature
                                   << " Mutation Ampl.: "<<mMutationAmplitude
                                   << " NEW OVERALL Best Cost="<<mBestCost<< endl;
                  }
                  else if(!silent) cout << "->Trial :" << mNbTrial 
                                   << " World="<< worldSwapIndex(i)
                                   << " Temp="<< mTemperature
                                   << " Mutation Ampl.: "<<mMutationAmplitude
                                   << " NEW RUN Best Cost="<<runBestCost<< endl;
                  if(!silent) this->DisplayReport();
               }
               worldNbAcceptedMoves(i)++;
            }
            else
            {
               if(log((rand()+1)/(REAL)RAND_MAX)<(-(cost-currentCost(i))/mTemperature) )
               {
                  accept=1;
                  currentCost(i)=cost;
                  mRefParList.SaveParamSet(worldCurrentSetIndex(i));
                  worldNbAcceptedMoves(i)++;
               }
            }
            switch(mSaveDetailledHistory.GetChoice())
            {
               case 0:break;
               case 1: 
               {
                  if(accept==2)
                  {
                     outHistory << mNbTrial <<" "<<i<<" "<<worldSwapIndex(i)<<" "<<accept<<" "
                                <<currentCost(i)<<" ";
                     //for(long j=0;j<mRefParList.GetNbParNotFixed();j++) 
                     //   outHistory << " "<< mRefParList.GetParNotFixed(j).GetHumanValue() ;
                     outHistory <<endl;
                  }
                  break;
               }
               case 2: 
               {
                  if(accept>0)
                  {
                     outHistory << mNbTrial <<" "<<i<<" "<<worldSwapIndex(i)<<" "<<accept<<" "
                                <<currentCost(i)<<" ";
                     //for(long j=0;j<mRefParList.GetNbParNotFixed();j++) 
                     //   outHistory << " "<< mRefParList.GetParNotFixed(j).GetHumanValue() ;
                     outHistory <<endl;
                  }
                  break;
               }
               case 3: 
               {
                  outHistory << mNbTrial <<" "<<i<<" "<<worldSwapIndex(i)<<" "<<accept<<" "
                             <<currentCost(i)<<" ";
                  //for(long j=0;j<mRefParList.GetNbParNotFixed();j++) 
                  //   outHistory << " "<< mRefParList.GetParNotFixed(j).GetHumanValue() ;
                  outHistory <<endl;
                  break;
               }
            }
            if(  ((mXMLAutoSave.GetChoice()==1)&&((chrono.seconds()-secondsWhenAutoSave)>86400))
               ||((mXMLAutoSave.GetChoice()==2)&&((chrono.seconds()-secondsWhenAutoSave)>3600))
               ||((mXMLAutoSave.GetChoice()==3)&&((chrono.seconds()-secondsWhenAutoSave)> 600))
               ||((mXMLAutoSave.GetChoice()==4)&&(accept==2)) )
            {
               secondsWhenAutoSave=(unsigned long)chrono.seconds();
               string saveFileName=this->GetName();
               time_t date=time(0);
               char strDate[40];
               strftime(strDate,sizeof(strDate),"%Y-%m-%d_%H-%M-%S",localtime(&date));//%Y-%m-%dT%H:%M:%S%Z
               char costAsChar[30];
               if(accept!=2) mRefParList.RestoreParamSet(mBestParSavedSetIndex);
               sprintf(costAsChar,"-Cost-%f",this->GetLogLikelihood());
               saveFileName=saveFileName+(string)strDate+(string)costAsChar+(string)".xml";
               XMLCrystFileSaveGlobal(saveFileName);
               //if(accept!=2) mRefParList.RestoreParamSet(lastParSavedSetIndex);
            }
            //if(accept==0) mRefParList.RestoreParamSet(lastParSavedSetIndex);
            mNbTrial++;nbStep--;
            if((mNbTrial%nbTrialsReport)==0) makeReport=true;
         }//nbTryPerWorld trials
      }//For each World
      //Try swapping worlds
      for(int i=1;i<nbWorld;i++)
      {
         #if 0
         mRefParList.RestoreParamSet(worldCurrentSetIndex(i));
         mMutationAmplitude=mutationAmplitude(i);
         cout<<i<<":"<<currentCost(i)<<":"<<this->GetLogLikelihood()<<endl;
         #endif
         #if 1
         if( log((rand()+1)/(REAL)RAND_MAX)
                < (-(currentCost(i-1)-currentCost(i))/simAnnealTemp(i)))
         #else
         // Compare World (i-1) and World (i) with the same amplitude,
         // hence the same max likelihood error
         mRefParList.RestoreParamSet(worldCurrentSetIndex(i-1));
         mMutationAmplitude=mutationAmplitude(i);
         if( log((rand()+1)/(REAL)RAND_MAX)
                < (-(this->GetLogLikelihood()-currentCost(i))/simAnnealTemp(i)))
         #endif
         {  
         /*
            if(i>2)
            {
               cout <<"->Swapping Worlds :" << i <<"(cost="<<currentCost(i)<<")"
                    <<" with "<< (i-1) <<"(cost="<< currentCost(i-1)<<")"<<endl;
            }
            */
            swapPar=mRefParList.GetParamSet(worldCurrentSetIndex(i));
            mRefParList.GetParamSet(worldCurrentSetIndex(i))=
               mRefParList.GetParamSet(worldCurrentSetIndex(i-1));
            mRefParList.GetParamSet(worldCurrentSetIndex(i-1))=swapPar;
            const REAL tmp=currentCost(i);
            currentCost(i)=currentCost(i-1);
            currentCost(i-1)=tmp;
            const long tmpIndex=worldSwapIndex(i);
            worldSwapIndex(i)=worldSwapIndex(i-1);
            worldSwapIndex(i-1)=tmpIndex;
            #if 0
            // Compute correct costs in the case we use maximum likelihood
            mRefParList.RestoreParamSet(worldCurrentSetIndex(i));
            mMutationAmplitude=mutationAmplitude(i);
            currentCost(i)=this->GetLogLikelihood();

            mRefParList.RestoreParamSet(worldCurrentSetIndex(i-1));
            mMutationAmplitude=mutationAmplitude(i-1);
            currentCost(i-1)=this->GetLogLikelihood();
            #endif
         }
      }
      #if 0
      //Try mating worlds- NEW !
      TAU_PROFILE_TIMER(timer1,\
               "MonteCarloObj::Optimize (Try mating Worlds)"\
               ,"", TAU_FIELD);
      TAU_PROFILE_START(timer1);
      if( (rand()/(REAL)RAND_MAX)<.1)
      for(int k=nbWorld-1;k>nbWorld/2;k--)
         for(int i=k-nbWorld/3;i<k;i++)
         {
            #if 0
            // Random switching of gene groups
            for(unsigned int j=0;j<nbGeneGroup;j++) 
               crossoverGroupIndex(j)= (int) floor(rand()/((REAL)RAND_MAX-1)*2);
            for(int j=0;j<mRefParList.GetNbPar();j++)
            {
               if(0==crossoverGroupIndex(refParGeneGroupIndex(j)-1))
               {
                  mRefParList.GetParamSet(parSetOffspringA)(j)=
                     mRefParList.GetParamSet(worldCurrentSetIndex(i))(j);
                  mRefParList.GetParamSet(parSetOffspringB)(j)=
                     mRefParList.GetParamSet(worldCurrentSetIndex(k))(j);
               }
               else
               {
                  mRefParList.GetParamSet(parSetOffspringA)(j)=
                     mRefParList.GetParamSet(worldCurrentSetIndex(k))(j);
                  mRefParList.GetParamSet(parSetOffspringB)(j)=
                     mRefParList.GetParamSet(worldCurrentSetIndex(i))(j);
               }
            }
            #endif
            #if 1
            // Switch gene groups in two parts
            unsigned int crossoverPoint1=
               (int)(1+floor(rand()/((REAL)RAND_MAX-1)*(nbGeneGroup)));
            unsigned int crossoverPoint2=
               (int)(1+floor(rand()/((REAL)RAND_MAX-1)*(nbGeneGroup)));
            if(crossoverPoint2<crossoverPoint1)
            {
               int tmp=crossoverPoint1;
               crossoverPoint1=crossoverPoint2;
               crossoverPoint2=tmp;
            }
            if(crossoverPoint1==crossoverPoint2) crossoverPoint2+=1;
            for(int j=0;j<mRefParList.GetNbPar();j++)
            {
               if((refParGeneGroupIndex(j)>crossoverPoint1)&&refParGeneGroupIndex(j)<crossoverPoint2)
               {
                  mRefParList.GetParamSet(parSetOffspringA)(j)=
                     mRefParList.GetParamSet(worldCurrentSetIndex(i))(j);
                  mRefParList.GetParamSet(parSetOffspringB)(j)=
                     mRefParList.GetParamSet(worldCurrentSetIndex(k))(j);
               }
               else
               {
                  mRefParList.GetParamSet(parSetOffspringA)(j)=
                     mRefParList.GetParamSet(worldCurrentSetIndex(k))(j);
                  mRefParList.GetParamSet(parSetOffspringB)(j)=
                     mRefParList.GetParamSet(worldCurrentSetIndex(i))(j);
               }
            }
            #endif
            // Try both offspring
            for(int junk=0;junk<2;junk++)
            {
               if(junk==0) mRefParList.RestoreParamSet(parSetOffspringA);
               else mRefParList.RestoreParamSet(parSetOffspringB);
               REAL cost=this->GetLogLikelihood();
               //if(log((rand()+1)/(REAL)RAND_MAX)
               //    < (-(cost-currentCost(k))/simAnnealTemp(k)))
               if(cost<currentCost(k))
               {
                  // Also exchange genes for higher-temperature World ?
                  //if(junk==0) 
                  //   mRefParList.GetParamSet(worldCurrentSetIndex(i))=
                  //      mRefParList.GetParamSet(parSetOffspringB);
                  //else
                  //   mRefParList.GetParamSet(worldCurrentSetIndex(i))=
                  //      mRefParList.GetParamSet(parSetOffspringA);
                  currentCost(k)=cost;
                  mRefParList.SaveParamSet(worldCurrentSetIndex(k));
                  //worldNbAcceptedMoves(k)++;
                  if(!silent) cout << "Accepted mating :"<<k<<"(with"<<i<<")"
                       <<" (crossoverGene1="<< crossoverPoint1<<","
                       <<" crossoverGene2="<< crossoverPoint2<<")"
                       <<endl; 
                  if(cost<runBestCost)
                  {
                     runBestCost=cost;
                     this->TagNewBestConfig();
                     mRefParList.SaveParamSet(runBestIndex);
                     if(cost<mBestCost)
                     {
                        mBestCost=cost;
                        mRefParList.SaveParamSet(mBestParSavedSetIndex);
                        if(!silent) cout << "->Trial :" << mNbTrial 
                          << " World="<< worldSwapIndex(k)
                          << " Temp="<< simAnnealTemp(k)
                          << " Mutation Ampl.: "<<mMutationAmplitude
                          << " NEW OVERALL Best Cost="<<mBestCost<< "(MATING !)"<<endl;
                     }
                     else if(!silent) cout << "->Trial :" << mNbTrial 
                          << " World="<< worldSwapIndex(k)
                          << " Temp="<< simAnnealTemp(k)
                          << " Mutation Ampl.: "<<mMutationAmplitude
                          << " NEW RUN Best Cost="<<runBestCost<< "(MATING !)"<<endl;
                     bestConfigNb=mNbTrial;
                     if(!silent) this->DisplayReport();
                     //for(int i=0;i<mRefinedObjList.GetNb();i++) 
                     //   mRefinedObjList.GetObj(i).Print();
                  }
                  i=k;//Don't test other Worlds
                  break;
               }
               //mNbTrial++;nbStep--;
               //if((mNbTrial%nbTrialsReport)==0) makeReport=true;
            }
         }
      TAU_PROFILE_STOP(timer1);
      #endif
      if(true==makeReport)
      {
         makeReport=false;
         worldNbAcceptedMoves*=nbWorld;
         if(!silent)
         {
            #if 0
            {// Experimental, dynamical weighting
               REAL max=0.;
               map<const RefinableObj*,REAL> ll,llvar;
               map<const RefinableObj*,LogLikelihoodStats>::iterator pos;
               for(pos=mvContextObjStats[0].begin();pos!=mvContextObjStats[0].end();++pos)
               {
                  ll   [pos->first]=0.;
                  llvar[pos->first]=0.;
               }
               for(int i=0;i<nbWorld;i++)
               {
                  for(pos=mvContextObjStats[0].begin();pos!=mvContextObjStats[0].end();++pos)
                  {
                     ll   [pos->first] += pos->second.mTotalLogLikelihood;
                     llvar[pos->first] += pos->second.mTotalLogLikelihoodDeltaSq;
                  }
               }
               for(pos=mvContextObjStats[0].begin();pos!=mvContextObjStats[0].end();++pos)
               {
                  cout << pos->first->GetName()
                       << " " << llvar[pos->first]
                       << " " << mvObjWeight[pos->first].mWeight
                       << " " << max<<endl;
                  llvar[pos->first] *= mvObjWeight[pos->first].mWeight;
                  if(llvar[pos->first]>max) max=llvar[pos->first];
               }
               map<const RefinableObj*,REAL>::iterator pos2;
               for(pos2=llvar.begin();pos2!=llvar.end();++pos2)
               {
                  const REAL d=pos2->second;
                  if(d<(max/mvObjWeight.size()/10.))
                  {
                     if(d<1) continue;
                     mvObjWeight[pos2->first].mWeight *=2;
                  }
               }
               REAL ll1=0;
               REAL llt=0;
               for(pos2=ll.begin();pos2!=ll.end();++pos2)
               {
                  llt += pos2->second;
                  ll1 += pos2->second * mvObjWeight[pos2->first].mWeight;
               }
               map<const RefinableObj*,DynamicObjWeight>::iterator posw;
               for(posw=mvObjWeight.begin();posw!=mvObjWeight.end();++posw)
               {
                  posw->second.mWeight *= llt/ll1;
               }
            }
            #endif //Experimental dynamical weighting
            #if 1 //def __DEBUG__
            for(int i=0;i<nbWorld;i++)
            {
               cout<<"   World :"<<worldSwapIndex(i)<<":";
               map<const RefinableObj*,LogLikelihoodStats>::iterator pos;
               for(pos=mvContextObjStats[i].begin();pos!=mvContextObjStats[i].end();++pos)
               {
                  cout << pos->first->GetName() 
                       << "(LLK="
                       << pos->second.mLastLogLikelihood
                       //<< "(<LLK>="
                       //<< pos->second.mTotalLogLikelihood/nbTrialsReport
                       //<< ", <delta(LLK)^2>="
                       //<< pos->second.mTotalLogLikelihoodDeltaSq/nbTrialsReport
                       << ", w="<<mvObjWeight[pos->first].mWeight
                       <<")  ";
                  pos->second.mTotalLogLikelihood=0;
                  pos->second.mTotalLogLikelihoodDeltaSq=0;
               }
               cout << endl;
            }
            #endif
            for(int i=0;i<nbWorld;i++)
            {
               //mRefParList.RestoreParamSet(worldCurrentSetIndex(i));
               cout <<"   World :" << worldSwapIndex(i)
                    <<" Temp.: " << simAnnealTemp(i)
                    <<" Mutation Ampl.: " << mutationAmplitude(i)
                    <<" Current Cost=" << currentCost(i)
                    <<" Accepting "
                    << (int)((REAL)worldNbAcceptedMoves(i)/nbTrialsReport*100)
                    <<"% moves " <<endl;
               //     <<"% moves " << mRefParList.GetPar("Pboccup").GetValue()<<endl;
            }
         }
         if(!silent) cout <<"Trial :" << mNbTrial << " Best Cost=" << runBestCost<< " ";
         if(!silent) chrono.print();
         //Change the mutation rate if necessary for each world
         if(ANNEALING_SMART==mAnnealingScheduleMutation.GetChoice())
         {
            for(int i=0;i<nbWorld;i++)
            {
               if((worldNbAcceptedMoves(i)/(REAL)nbTrialsReport)>0.30)
                  mutationAmplitude(i)*=2.;
               if((worldNbAcceptedMoves(i)/(REAL)nbTrialsReport)<0.10)
                  mutationAmplitude(i)/=2.;
               if(mutationAmplitude(i)>mMutationAmplitudeMax) 
                  mutationAmplitude(i)=mMutationAmplitudeMax;
               if(mutationAmplitude(i)<mMutationAmplitudeMin)
                  mutationAmplitude(i)=mMutationAmplitudeMin;
            }
         }
         if(ANNEALING_SMART==mAnnealingScheduleTemp.GetChoice())
         {
            for(int i=0;i<nbWorld;i++)
            {
               if((worldNbAcceptedMoves(i)/(REAL)nbTrialsReport)>0.30)
                  simAnnealTemp(i)/=1.5;
               if((worldNbAcceptedMoves(i)/(REAL)nbTrialsReport)>0.80)
                  simAnnealTemp(i)/=2;
               if((worldNbAcceptedMoves(i)/(REAL)nbTrialsReport)>0.95)
                  simAnnealTemp(i)/=4;

               if((worldNbAcceptedMoves(i)/(REAL)nbTrialsReport)<0.10)
                  simAnnealTemp(i)*=1.5;
               if((worldNbAcceptedMoves(i)/(REAL)nbTrialsReport)<0.04)
                  simAnnealTemp(i)*=2;
               if((worldNbAcceptedMoves(i)/(REAL)nbTrialsReport)<0.01)
                  simAnnealTemp(i)*=4;
               //cout<<"World#"<<i<<":"<<worldNbAcceptedMoves(i)<<":"<<nbTrialsReport<<endl;
               //if(simAnnealTemp(i)>mTemperatureMax) simAnnealTemp(i)=mTemperatureMax;
               //if(simAnnealTemp(i)<mTemperatureMin) simAnnealTemp(i)=mTemperatureMin;
            }
         }
         worldNbAcceptedMoves=0;
         //this->DisplayReport();

         #ifdef __WX__CRYST__
         if(0!=mpWXCrystObj) mpWXCrystObj->UpdateDisplayNbTrial();
         #endif
      }
      if((runBestCost<finalcost) || mStopAfterCycle ||( (maxTime>0)&&(chrono.seconds()>maxTime))) 
      {
         if(!silent) cout << endl <<endl << "Refinement Stopped:"<<mBestCost<<endl;
         break;
      }
   }//Trials
   mLastOptimTime=chrono.seconds();
   //Restore Best values
      //mRefParList.Print();
      if(!silent) this->DisplayReport();
      mRefParList.RestoreParamSet(runBestIndex);
      //for(int i=0;i<mRefinedObjList.GetNb();i++) mRefinedObjList.GetObj(i).Print();
      mCurrentCost=this->GetLogLikelihood();
      if(!silent) cout<<"Run Best Cost:"<<mCurrentCost<<endl;
      if(!silent) chrono.print();
   //Save density of states
      //ofstream out("densityOfStates.txt");
      //out << trialsDensity<<endl;
      //out.close();
   // Clear temporary param set
      for(int i=0;i<nbWorld;i++)
      {
         mRefParList.ClearParamSet(worldCurrentSetIndex(i));
         //mvSavedParamSet.push_back(make_pair(worldCurrentSetIndex(i),currentCost(i)));
      }
      mRefParList.ClearParamSet(lastParSavedSetIndex);
      mRefParList.ClearParamSet(runBestIndex);
   if(mSaveDetailledHistory.GetChoice()>0) outHistory.close();
}

/*
void MonteCarloObj::SaveOptimHistory() const
{
   VFN_DEBUG_MESSAGE("MonteCarloObj::SaveOptimHistory()",5)
   if(mHistoryNb<=0)
   {
      cout << "MonteCarloObj::SaveOptimHistory(): No History yet !!! "<<endl;
      return;
   }
   
   ofstream out(mHistorySaveFileName.c_str());
   if(!out)
   {
      throw ObjCrystException("MonteCarloObj::SaveOptimHistory() : \
Error opening file for output:"+mHistorySaveFileName);
   }
   out << "# Trial Cost";
   for(long j=0;j<mRefParList.GetNbParNotFixed();j++)
      out << mRefParList.GetParNotFixed(j).GetName() << " ";
   out <<endl;
   for(long i=0;i<mHistoryNb;i++)
   {
      const long saveSet=mHistorySavedParamSetIndex(i);
      out << mHistoryTrialNumber(i) << " " << mHistoryCostFunction(i) ;
      for(long j=0;j<mRefParList.GetNbParNotFixed();j++) 
         out << " "<< mRefParList.GetParamSet_ParNotFixedHumanValue(saveSet,j) ;
      out <<endl;
   }
   out.close();
}
*/
void MonteCarloObj::XMLOutput(ostream &os,int indent)const
{
   VFN_DEBUG_ENTRY("MonteCarloObj::XMLOutput():"<<this->GetName(),5)
   for(int i=0;i<indent;i++) os << "  " ;
   XMLCrystTag tag("GlobalOptimObj");
   tag.AddAttribute("Name",this->GetName());
   
   os <<tag<<endl;
   indent++;
   
   mGlobalOptimType.XMLOutput(os,indent);
   os<<endl;

   mAnnealingScheduleTemp.XMLOutput(os,indent);
   os<<endl;

   mXMLAutoSave.XMLOutput(os,indent);
   os<<endl;

   {
      XMLCrystTag tag2("TempMaxMin");
      for(int i=0;i<indent;i++) os << "  " ;
      os<<tag2<<mTemperatureMax << " "<< mTemperatureMin;
      tag2.SetIsEndTag(true);
      os<<tag2<<endl;
   }

   mAnnealingScheduleMutation.XMLOutput(os,indent);
   os<<endl;
   
   {
      XMLCrystTag tag2("MutationMaxMin");
      for(int i=0;i<indent;i++) os << "  " ;
      os<<tag2<<mMutationAmplitudeMax << " "<< mMutationAmplitudeMin;
      tag2.SetIsEndTag(true);
      os<<tag2<<endl;
   }

   for(int j=0;j<mRefinedObjList.GetNb();j++)
   {
      XMLCrystTag tag2("RefinedObject",false,true);
      tag2.AddAttribute("ObjectType",mRefinedObjList.GetObj(j).GetClassName());
      tag2.AddAttribute("ObjectName",mRefinedObjList.GetObj(j).GetName());
      for(int i=0;i<indent;i++) os << "  " ;
      os<<tag2<<endl;
   }
   
   indent--;
   tag.SetIsEndTag(true);
   for(int i=0;i<indent;i++) os << "  " ;
   os <<tag<<endl;
   VFN_DEBUG_EXIT("MonteCarloObj::XMLOutput():"<<this->GetName(),5)
}

void MonteCarloObj::XMLInput(istream &is,const XMLCrystTag &tagg)
{
   VFN_DEBUG_ENTRY("MonteCarloObj::XMLInput():"<<this->GetName(),5)
   for(unsigned int i=0;i<tagg.GetNbAttribute();i++)
   {
      if("Name"==tagg.GetAttributeName(i)) this->SetName(tagg.GetAttributeValue(i));
   }
   while(true)
   {
      XMLCrystTag tag(is);
      if(("GlobalOptimObj"==tag.GetName())&&tag.IsEndTag())
      {
         VFN_DEBUG_EXIT("MonteCarloObj::Exit():"<<this->GetName(),5)
         this->UpdateDisplay();
         return;
      }
      if("Option"==tag.GetName())
      {
         for(unsigned int i=0;i<tag.GetNbAttribute();i++)
            if("Name"==tag.GetAttributeName(i))
            {
               if("Algorithm"==tag.GetAttributeValue(i))
               {
                  mGlobalOptimType.XMLInput(is,tag);
                  break;
               }
               if("Temperature Schedule"==tag.GetAttributeValue(i))
               {
                  mAnnealingScheduleTemp.XMLInput(is,tag);
                  break;
               }
               if("Displacement Amplitude Schedule"==tag.GetAttributeValue(i))
               {
                  mAnnealingScheduleMutation.XMLInput(is,tag);
                  break;
               }
               if("Save Best Config Regularly"==tag.GetAttributeValue(i))
               {
                  mXMLAutoSave.XMLInput(is,tag);
                  break;
               }
            }
         continue;
      }
      if("TempMaxMin"==tag.GetName())
      {
         is>>mTemperatureMax>> mTemperatureMin;
         if(false==tag.IsEmptyTag()) XMLCrystTag junk(is);//:KLUDGE: for first release
         continue;
      }
      if("MutationMaxMin"==tag.GetName())
      {
         is>>mMutationAmplitudeMax>>mMutationAmplitudeMin;
         if(false==tag.IsEmptyTag()) XMLCrystTag junk(is);//:KLUDGE: for first release
         continue;
      }
      if("RefinedObject"==tag.GetName())
      {
         string name,type;
         for(unsigned int i=0;i<tag.GetNbAttribute();i++)
         {
            if("ObjectName"==tag.GetAttributeName(i)) name=tag.GetAttributeValue(i);
            if("ObjectType"==tag.GetAttributeName(i)) type=tag.GetAttributeValue(i);
         }
         RefinableObj* obj=& (gRefinableObjRegistry.GetObj(name,type));
         this->AddRefinableObj(*obj);
         continue;
      }
   }
}

const string MonteCarloObj::GetClassName()const { return "MonteCarloObj";}

REAL MonteCarloObj::GetLogLikelihood()
{
   #if 0
   const REAL mle=0.10;//+0.025*mMutationAmplitude;//*sqrt(mMutationAmplitude);
   mRefParList.GetPar(7).SetValue(mle);
   mRefParList.GetPar(9).SetValue(mle);
   mRefParList.GetPar(11).SetValue(mle);
   //if((rand()%1000)==0) mRefParList.Print();
   #endif
   return OptimizationObj::GetLogLikelihood();
}

void MonteCarloObj::NewConfiguration(const RefParType *type)
{
   VFN_DEBUG_ENTRY("MonteCarloObj::NewConfiguration()",4)
   for(int i=0;i<mRefinedObjList.GetNb();i++)
      mRefinedObjList.GetObj(i).BeginGlobalOptRandomMove();
   for(int i=0;i<mRefinedObjList.GetNb();i++)
      mRefinedObjList.GetObj(i).GlobalOptRandomMove(mMutationAmplitude,type);
   VFN_DEBUG_EXIT("MonteCarloObj::NewConfiguration()",4)
}

void MonteCarloObj::InitOptions()
{
   VFN_DEBUG_MESSAGE("MonteCarloObj::InitOptions()",5)
   this->OptimizationObj::InitOptions();
   static string GlobalOptimTypeName;
   static string GlobalOptimTypeChoices[2];//:TODO: Add Genetic Algorithm
   
   static string AnnealingScheduleChoices[6];
   
   static string AnnealingScheduleTempName;
   static string AnnealingScheduleMutationName;
   
   static string saveDetailledHistoryName;
   static string saveDetailledHistoryChoices[4];

   static bool needInitNames=true;
   if(true==needInitNames)
   {
      GlobalOptimTypeName="Algorithm";
      GlobalOptimTypeChoices[0]="Simulated Annealing";
      GlobalOptimTypeChoices[1]="Parallel Tempering";
      //GlobalOptimTypeChoices[2]="Genetic";

      AnnealingScheduleTempName="Temperature Schedule";
      AnnealingScheduleMutationName="Displacement Amplitude Schedule";
      AnnealingScheduleChoices[0]="Constant";
      AnnealingScheduleChoices[1]="Boltzmann";
      AnnealingScheduleChoices[2]="Cauchy";
      AnnealingScheduleChoices[3]="Exponential";
      AnnealingScheduleChoices[4]="Smart";
      AnnealingScheduleChoices[5]="Gamma";
      
      saveDetailledHistoryName="Save evolution of all parameters (for tests ONLY!!!)";
      saveDetailledHistoryChoices[0]="No (highly recommended)";
      saveDetailledHistoryChoices[1]="Every new best configuration";
      saveDetailledHistoryChoices[2]="Every accepted configuration";
      saveDetailledHistoryChoices[3]="Every configuration";

      needInitNames=false;//Only once for the class
   }
   mGlobalOptimType.Init(2,&GlobalOptimTypeName,GlobalOptimTypeChoices);
   mAnnealingScheduleTemp.Init(6,&AnnealingScheduleTempName,AnnealingScheduleChoices);
   mAnnealingScheduleMutation.Init(6,&AnnealingScheduleMutationName,AnnealingScheduleChoices);
   mSaveDetailledHistory.Init(4,&saveDetailledHistoryName,saveDetailledHistoryChoices);
   VFN_DEBUG_MESSAGE("MonteCarloObj::InitOptions():End",5)
}

#ifdef __WX__CRYST__
WXCrystObjBasic* MonteCarloObj::WXCreate(wxWindow *parent)
{
   mpWXCrystObj=new WXMonteCarloObj (parent,this);
   return mpWXCrystObj;
}
WXOptimizationObj* MonteCarloObj::WXGet()
{
   return mpWXCrystObj;
}
void MonteCarloObj::WXDelete()
{
   if(0!=mpWXCrystObj) delete mpWXCrystObj;
   mpWXCrystObj=0;
}
void MonteCarloObj::WXNotifyDelete()
{
   mpWXCrystObj=0;
}
#endif

}//namespace
