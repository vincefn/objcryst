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

#include "ObjCryst/RefinableObj/GlobalOptimObj.h"
#include "ObjCryst/ObjCryst/Crystal.h"
#include "ObjCryst/Quirks/VFNStreamFormat.h"
#include "ObjCryst/Quirks/VFNDebug.h"
#include "ObjCryst/Quirks/Chronometer.h"
#include "ObjCryst/ObjCryst/IO.h"
#include "ObjCryst/RefinableObj/LSQNumObj.h"

#include "ObjCryst/ObjCryst/Molecule.h"

#ifdef __WX__CRYST__
   #include "ObjCryst/wxCryst/wxRefinableObj.h"
   #undef GetClassName // Conflict from wxMSW headers ? (cygwin)
#endif

//For some reason, with wxWindows this must be placed after wx headers (Borland c++)
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <boost/format.hpp>

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
mNbTrialPerRun(10000000),mNbTrial(0),mBestCost(-1),
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

   static bool need_initRandomSeed=true;
   if(need_initRandomSeed==true)
   {
      srand(time(NULL));
      need_initRandomSeed=false;
   }
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

REAL OptimizationObj::GetLogLikelihood() const
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
   if(mIsOptimizing)
   {
      #ifdef __WX__CRYST__
      wxMutexLocker lock(mMutexStopAfterCycle);
      #endif
      mStopAfterCycle=true;
   }
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
   mMainTracker.AppendValues(mNbTrial);
}

REAL OptimizationObj::GetLastOptimElapsedTime()const
{
   return mLastOptimTime;
}

MainTracker& OptimizationObj::GetMainTracker(){return mMainTracker;}

const MainTracker& OptimizationObj::GetMainTracker()const{return mMainTracker;}

RefObjOpt& OptimizationObj::GetXMLAutoSaveOption() {return mXMLAutoSave;}
const RefObjOpt& OptimizationObj::GetXMLAutoSaveOption()const {return mXMLAutoSave;}

const REAL& OptimizationObj::GetBestCost()const{return mBestCost;}
REAL& OptimizationObj::GetBestCost(){return mBestCost;}

void OptimizationObj::BeginOptimization(const bool allowApproximations, const bool enableRestraints)
{
   mvContextObjStats.clear();
   for(int i=0;i<mRefinedObjList.GetNb();i++)
   {
      mRefinedObjList.GetObj(i).BeginOptimization(allowApproximations,enableRestraints);
   }
}

void OptimizationObj::EndOptimization()
{
   for(int i=0;i<mRefinedObjList.GetNb();i++) mRefinedObjList.GetObj(i).EndOptimization();
}

long& OptimizationObj::NbTrialPerRun() {return mNbTrialPerRun;}

const long& OptimizationObj::NbTrialPerRun() const {return mNbTrialPerRun;}

unsigned int OptimizationObj::GetNbOption()const
{
   return mOptionRegistry.GetNb();
}

ObjRegistry<RefObjOpt>& OptimizationObj::GetOptionList()
{
   return mOptionRegistry;
}

RefObjOpt& OptimizationObj::GetOption(const unsigned int i)
{
   VFN_DEBUG_MESSAGE("RefinableObj::GetOption()"<<i,3)
   //:TODO: Check
   return mOptionRegistry.GetObj(i);
}

RefObjOpt& OptimizationObj::GetOption(const string & name)
{
   VFN_DEBUG_MESSAGE("OptimizationObj::GetOption()"<<name,3)
   const long i=mOptionRegistry.Find(name);
   if(i<0)
   {
      this->Print();
      throw ObjCrystException("OptimizationObj::GetOption(): cannot find option: "+name+" in object:"+this->GetName());
   }
   return mOptionRegistry.GetObj(i);
}

const RefObjOpt& OptimizationObj::GetOption(const unsigned int i)const
{
   VFN_DEBUG_MESSAGE("RefinableObj::GetOption()"<<i,3)
   //:TODO: Check
   return mOptionRegistry.GetObj(i);
}

const RefObjOpt& OptimizationObj::GetOption(const string & name)const
{
   VFN_DEBUG_MESSAGE("OptimizationObj::GetOption()"<<name,3)
   const long i=mOptionRegistry.Find(name);
   if(i<0)
   {
      this->Print();
      throw ObjCrystException("OptimizationObj::GetOption(): cannot find option: "+name+" in object:"+this->GetName());
   }
   return mOptionRegistry.GetObj(i);
}

const ObjRegistry<RefinableObj>& OptimizationObj::GetRefinedObjList() const
{
   return mRefinedObjList;
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
         mRefParList.AddPar(mRecursiveRefinedObjList.GetObj(i),true);
      mvSavedParamSet.clear();
      mBestParSavedSetIndex=mRefParList.CreateParamSet("Best Configuration");
      mvSavedParamSet.push_back(make_pair(mBestParSavedSetIndex,mBestCost));

      mMainTracker.ClearTrackers();

      REAL (OptimizationObj::*fl)() const;
      fl=&OptimizationObj::GetLogLikelihood;
      mMainTracker.AddTracker(new TrackerObject<OptimizationObj>
         (this->GetName()+"::Overall LogLikelihood",*this,fl));

      for(long i=0;i<mRecursiveRefinedObjList.GetNb();i++)
      {
         REAL (RefinableObj::*fp)() const;
         fp=&RefinableObj::GetLogLikelihood;
         mMainTracker.AddTracker(new TrackerObject<RefinableObj>
            (mRecursiveRefinedObjList.GetObj(i).GetName()+"::LogLikelihood",mRecursiveRefinedObjList.GetObj(i),fp));

         if(mRecursiveRefinedObjList.GetObj(i).GetClassName()=="Crystal")
         {
            REAL (Crystal::*fc)() const;
            const Crystal *pCryst=dynamic_cast<const Crystal *>(&(mRecursiveRefinedObjList.GetObj(i)));
            fc=&Crystal::GetBumpMergeCost;
            mMainTracker.AddTracker(new TrackerObject<Crystal>
               (pCryst->GetName()+"::BumpMergeCost",*pCryst,fc));
            fc=&Crystal::GetBondValenceCost;
            mMainTracker.AddTracker(new TrackerObject<Crystal>
               (pCryst->GetName()+"::BondValenceCost",*pCryst,fc));
         }
      }
   }
   // Prepare for refinement, ie get the list of not fixed parameters,
   // and prepare the objects...
   mRefParList.PrepareForRefinement();
   for(int i=0;i<mRecursiveRefinedObjList.GetNb();i++)
      mRecursiveRefinedObjList.GetObj(i).PrepareForRefinement();
   VFN_DEBUG_EXIT("OptimizationObj::PrepareRefParList()",6)
}

void OptimizationObj::InitOptions()
{
   VFN_DEBUG_MESSAGE("OptimizationObj::InitOptions()",5)
   static string xmlAutoSaveName;
   static string xmlAutoSaveChoices[6];

   static bool needInitNames=true;
   if(true==needInitNames)
   {
      xmlAutoSaveName="Save Best Config Regularly";
      xmlAutoSaveChoices[0]="No";
      xmlAutoSaveChoices[1]="Every day";
      xmlAutoSaveChoices[2]="Every hour";
      xmlAutoSaveChoices[3]="Every 10mn";
      xmlAutoSaveChoices[4]="Every new best config (a lot ! Not Recommended !)";
      xmlAutoSaveChoices[5]="Every Run (Recommended)";

      needInitNames=false;//Only once for the class
   }
   mXMLAutoSave.Init(6,&xmlAutoSaveName,xmlAutoSaveChoices);
   this->AddOption(&mXMLAutoSave);
   VFN_DEBUG_MESSAGE("OptimizationObj::InitOptions():End",5)
}

void OptimizationObj::UpdateDisplay()
{
   Chronometer chrono;
   #ifdef __WX__CRYST__
   if(0!=this->WXGet()) this->WXGet()->CrystUpdate(true,true);
   #endif
   for(int i=0;i<mRefinedObjList.GetNb();i++)
      mRefinedObjList.GetObj(i).UpdateDisplay();
   mMainTracker.UpdateDisplay();
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

void OptimizationObj::AddOption(RefObjOpt *opt)
{
   VFN_DEBUG_ENTRY("OptimizationObj::AddOption()",5)
   mOptionRegistry.Register(*opt);
   VFN_DEBUG_EXIT("OptimizationObj::AddOption()",5)
}

//#################################################################################
//
//       MonteCarloObj
//
//#################################################################################
MonteCarloObj::MonteCarloObj(const string name):
OptimizationObj(name),
mCurrentCost(-1),
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
   mXMLAutoSave.SetChoice(5);//Save after each Run
   mAutoLSQ.SetChoice(0);
   gOptimizationObjRegistry.Register(*this);
   VFN_DEBUG_EXIT("MonteCarloObj::MonteCarloObj()",5)
}

MonteCarloObj::MonteCarloObj(const bool internalUseOnly):
OptimizationObj(""),
mCurrentCost(-1),
mTemperatureMax(.03),mTemperatureMin(.003),mTemperatureGamma(1.0),
mMutationAmplitudeMax(16.),mMutationAmplitudeMin(.125),mMutationAmplitudeGamma(1.0),
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
   mXMLAutoSave.SetChoice(5);//Save after each Run
   mAutoLSQ.SetChoice(0);
   if(false==internalUseOnly) gOptimizationObjRegistry.Register(*this);
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
   this->BeginOptimization(true);
   this->PrepareRefParList();

   this->InitLSQ(false);

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
   mMainTracker.ClearValues();
   Chronometer chrono;
   chrono.start();
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
      case GLOBAL_OPTIM_RANDOM_LSQ: //:TODO:
      {
          long cycles = 1;
          this->RunRandomLSQMethod(cycles);
          break;
      }//case GLOBAL_OPTIM_GENETIC
   }
   mIsOptimizing=false;
   #ifdef __WX__CRYST__
   mMutexStopAfterCycle.Lock();
   #endif
   mStopAfterCycle=false;
   #ifdef __WX__CRYST__
   mMutexStopAfterCycle.Unlock();
   #endif

   mRefParList.RestoreParamSet(mBestParSavedSetIndex);
   this->EndOptimization();
   (*fpObjCrystInformUser)((boost::format("Finished Optimization, final cost=%12.2f (dt=%.1fs)") % this->GetLogLikelihood() % chrono.seconds()).str());

   if(mSaveTrackedData.GetChoice()==1)
   {
      ofstream outTracker;
      outTracker.imbue(std::locale::classic());
      const string outTrackerName=this->GetName()+"-Tracker.dat";
      outTracker.open(outTrackerName.c_str());
      mMainTracker.SaveAll(outTracker);
      outTracker.close();
   }

   for(vector<pair<long,REAL> >::iterator pos=mvSavedParamSet.begin();pos!=mvSavedParamSet.end();++pos)
      if(pos->first==mBestParSavedSetIndex)
      {
         if(  (pos->second>mBestCost)
            ||(pos->second<0))
         {
            pos->second=mBestCost;
            break;
         }
      }

   this->UpdateDisplay();

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
   this->BeginOptimization(true);
   this->PrepareRefParList();

   this->InitLSQ(false);

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
   long nbTrialCumul=0;
   const long nbCycle0=nbCycle;
	Chronometer chrono;
   while(nbCycle!=0)
   {
      if(!silent) cout <<"MonteCarloObj::MultiRunOptimize: Starting Run#"<<abs(nbCycle)<<endl;
      nbStep=nbStep0;
      for(int i=0;i<mRefinedObjList.GetNb();i++) mRefinedObjList.GetObj(i).RandomizeConfiguration();
      mMainTracker.ClearValues();
      chrono.start();
      switch(mGlobalOptimType.GetChoice())
      {
         case GLOBAL_OPTIM_SIMULATED_ANNEALING:
         {
            try{this->RunSimulatedAnnealing(nbStep,silent,finalcost,maxTime);}
            catch(...){cout<<"Unhandled exception in MonteCarloObj::MultiRunOptimize() ?"<<endl;}
            break;
         }
         case GLOBAL_OPTIM_PARALLEL_TEMPERING:
         {
            try{this->RunParallelTempering(nbStep,silent,finalcost,maxTime);}
            catch(...){cout<<"Unhandled exception in MonteCarloObj::MultiRunOptimize() ?"<<endl;}
            break;
         }
         case GLOBAL_OPTIM_RANDOM_LSQ:
         {
            try{this->RunRandomLSQMethod(nbCycle);}
            catch(...){cout<<"Unhandled exception in MonteCarloObj::RunRandomLSQMethod() ?"<<endl;}
            //nbCycle=1;
            break;
         }
      }
      nbTrialCumul+=(nbStep0-nbStep);
      if(finalcost>1)
         (*fpObjCrystInformUser)((boost::format("Finished Run #%d, final cost=%12.2f, nbTrial=%d (dt=%.1fs), so far <nbTrial>=%d")
                                  % (nbCycle0-nbCycle) % this->GetLogLikelihood() % (nbStep0-nbStep) % chrono.seconds() % (nbTrialCumul/(nbCycle0-nbCycle+1))).str());
      else
         (*fpObjCrystInformUser)((boost::format("Finished Run #%d, final cost=%12.2f, nbTrial=%d (dt=%.1fs)")
                                  % (nbCycle0-nbCycle) % this->GetLogLikelihood() % (nbStep0-nbStep) % chrono.seconds()).str());


      nbStep=nbStep0;
      if(false==mStopAfterCycle) this->UpdateDisplay();
      stringstream s;
      s<<"Run #"<<abs(nbCycle);
      mvSavedParamSet.push_back(make_pair(mRefParList.CreateParamSet(s.str()),mCurrentCost));
      if(!silent) cout <<"MonteCarloObj::MultiRunOptimize: Finished Run#"
                       <<abs(nbCycle)<<", Run Best Cost:"<<mCurrentCost
                       <<", Overall Best Cost:"<<mBestCost<<endl;
      if(mXMLAutoSave.GetChoice()==5)
      {
         string saveFileName=this->GetName();
         time_t date=time(0);
         char strDate[40];
         strftime(strDate,sizeof(strDate),"%Y-%m-%d_%H-%M-%S",localtime(&date));//%Y-%m-%dT%H:%M:%S%Z
         char costAsChar[30];
         sprintf(costAsChar,"-Run#%ld-Cost-%f",abs(nbCycle),this->GetLogLikelihood());
         saveFileName=saveFileName+(string)strDate+(string)costAsChar+(string)".xml";
         XMLCrystFileSaveGlobal(saveFileName);
      }
      if(mSaveTrackedData.GetChoice()==1)
      {
         ofstream outTracker;
         outTracker.imbue(std::locale::classic());
         char runNum[40];
         sprintf(runNum,"-Tracker-Run#%ld.dat",abs(nbCycle));
         const string outTrackerName=this->GetName()+runNum;
         outTracker.open(outTrackerName.c_str());
         mMainTracker.SaveAll(outTracker);
         outTracker.close();
      }
      nbCycle--;
      #ifdef __WX__CRYST__
      mMutexStopAfterCycle.Lock();
      #endif
      if(mStopAfterCycle)
      {
         #ifdef __WX__CRYST__
         mMutexStopAfterCycle.Unlock();
         #endif
         break;
      }
      #ifdef __WX__CRYST__
      mMutexStopAfterCycle.Unlock();
      #endif
   }
   mIsOptimizing=false;

   mRefParList.RestoreParamSet(mBestParSavedSetIndex);

   for(vector<pair<long,REAL> >::iterator pos=mvSavedParamSet.begin();pos!=mvSavedParamSet.end();++pos)
      if(pos->first==mBestParSavedSetIndex)
      {
         if(  (pos->second>mBestCost)
            ||(pos->second<0))
         {
            pos->second=mBestCost;
            break;
         }
      }

   this->EndOptimization();

   if(false==mStopAfterCycle) this->UpdateDisplay();

   #ifdef __WX__CRYST__
   mMutexStopAfterCycle.Lock();
   #endif
   mStopAfterCycle=false;
   #ifdef __WX__CRYST__
   mMutexStopAfterCycle.Unlock();
   #endif

   if(finalcost>1)
      cout<<endl<<"Finished all runs, number of trials to reach cost="
          <<finalcost<<" : <nbTrial>="<<nbTrialCumul/(nbCycle0-nbCycle)<<endl;
   VFN_DEBUG_EXIT("MonteCarloObj::MultiRunOptimize()",5)
}

void MonteCarloObj::RunSimulatedAnnealing(long &nbStep,const bool silent,
                                          const REAL finalcost,const REAL maxTime)
{
   //Keep a copy of the total number of steps, and decrement nbStep
   const long nbSteps=nbStep;
   unsigned int accept;// 1 if last trial was accepted? 2 if new best config ? else 0
   mNbTrial=0;
   // time (in seconds) when last autoSave was made (if enabled)
      unsigned long secondsWhenAutoSave=0;

   if(!silent) cout << "Starting Simulated Annealing Optimization for"<<nbSteps<<" trials"<<endl;
   if(!silent) this->DisplayReport();
   REAL runBestCost;
   mCurrentCost=this->GetLogLikelihood();
   runBestCost=mCurrentCost;
   const long lastParSavedSetIndex=mRefParList.CreateParamSet("MonteCarloObj:Last parameters (SA)");
   const long runBestIndex=mRefParList.CreateParamSet("Best parameters for current run (SA)");
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

   // Do we need to update the display ?
   bool needUpdateDisplay=false;
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
            needUpdateDisplay=true;
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

      #ifdef __WX__CRYST__
      mMutexStopAfterCycle.Lock();
      #endif
      if((runBestCost<finalcost) || mStopAfterCycle ||( (maxTime>0)&&(chrono.seconds()>maxTime)))
      {
         #ifdef __WX__CRYST__
         mMutexStopAfterCycle.Unlock();
         #endif
         if(!silent) cout << endl <<endl << "Refinement Stopped."<<endl;
         break;
      }
      #ifdef __WX__CRYST__
      mMutexStopAfterCycle.Unlock();
      #endif
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
      if((mNbTrial%300==0)&&needUpdateDisplay)
      {
         this->UpdateDisplay();
         needUpdateDisplay=false;
         mRefParList.RestoreParamSet(lastParSavedSetIndex);
      }

   }
   //cout<<"Beginning final LSQ refinement? ... ";
   if(mAutoLSQ.GetChoice()>0)
   {// LSQ
      if(!silent) cout<<"Beginning final LSQ refinement"<<endl;
      for(int i=0;i<mRefinedObjList.GetNb();i++) mRefinedObjList.GetObj(i).SetApproximationFlag(false);
      mRefParList.RestoreParamSet(runBestIndex);
      mCurrentCost=this->GetLogLikelihood();
      try {mLSQ.Refine(-50,true,true,false,0.001);}
      catch(const ObjCrystException &except){};
      if(!silent) cout<<"LSQ cost: "<<mCurrentCost<<" -> "<<this->GetLogLikelihood()<<endl;

      // Need to go back to optimization with approximations allowed (they are not during LSQ)
      for(int i=0;i<mRefinedObjList.GetNb();i++) mRefinedObjList.GetObj(i).SetApproximationFlag(true);

      REAL cost=this->GetLogLikelihood();
      if(cost<mCurrentCost)
      {
         mCurrentCost=cost;
         mRefParList.SaveParamSet(lastParSavedSetIndex);
         if(mCurrentCost<runBestCost)
         {
            runBestCost=mCurrentCost;
            mRefParList.SaveParamSet(runBestIndex);
            if(runBestCost<mBestCost)
            {
               mBestCost=mCurrentCost;
               mRefParList.SaveParamSet(mBestParSavedSetIndex);
               if(!silent) cout << "LSQ : NEW OVERALL Best Cost="<<runBestCost<< endl;
            }
            else if(!silent) cout << " LSQ : NEW Run Best Cost="<<runBestCost<< endl;
         }
      }
      if(!silent) cout<<"Finished LSQ refinement"<<endl;
   }


   mLastOptimTime=chrono.seconds();
   //Restore Best values
   mRefParList.RestoreParamSet(runBestIndex);
   mRefParList.ClearParamSet(runBestIndex);
   mRefParList.ClearParamSet(lastParSavedSetIndex);
   mCurrentCost=this->GetLogLikelihood();
   if(!silent) this->DisplayReport();
   if(!silent) chrono.print();
}
/*
void MonteCarloObj::RunNondestructiveLSQRefinement(  int nbCycle,bool useLevenbergMarquardt,
                                            const bool silent, const bool callBeginEndOptimization,
                                            const float minChi2var )
{
    float bsigma=-1, bdelta=-1;
    float asigma=-1, adelta=-1;
    //set the sigma values lower - it makes the molecular model more stable for LSQ
    for(int i=0;i<mRefinedObjList.GetNb();i++) {
        if(mRefinedObjList.GetObj(i).GetClassName()=="Crystal") {
            try {
                Crystal * pCryst = dynamic_cast<Crystal *>(&(mRefinedObjList.GetObj(i)));
                for(int s=0;s<pCryst->GetScattererRegistry().GetNb();s++) {
                    Molecule *pMol=dynamic_cast<Molecule*>(&(pCryst->GetScatt(s)));
                    if(pMol==NULL) continue; // not a Molecule
                    for(vector<MolBond*>::iterator pos = pMol->GetBondList().begin(); pos != pMol->GetBondList().end();++pos) {
                        bsigma = (*pos)->GetLengthSigma();
                        bdelta = (*pos)->GetLengthDelta();
                        (*pos)->SetLengthDelta(0.02);
                        (*pos)->SetLengthSigma(0.001);
                    }
                    for(vector<MolBondAngle*>::iterator pos=pMol->GetBondAngleList().begin();pos != pMol->GetBondAngleList().end();++pos)
                    {
                        asigma = (*pos)->GetAngleSigma();
                        adelta = (*pos)->GetAngleDelta();
                        (*pos)->SetAngleDelta(0.2*DEG2RAD);
                        (*pos)->SetAngleSigma(0.01*DEG2RAD);
                    }
                }
            } catch (const std::bad_cast& e) {

            }
        }
    }
    for(int i=0;i<mRefinedObjList.GetNb();i++) mRefinedObjList.GetObj(i).SetApproximationFlag(false);
    try {
        mLSQ.Refine(nbCycle,useLevenbergMarquardt,silent,callBeginEndOptimization,minChi2var);
    }
    catch(const ObjCrystException &except) {

    };
    for(int i=0;i<mRefinedObjList.GetNb();i++) mRefinedObjList.GetObj(i).SetApproximationFlag(true);

    if(bsigma<0 || bdelta<0 || asigma<0 || adelta<0) return;
    //restore the delta and sigma values
    for(int i=0;i<mRefinedObjList.GetNb();i++) {
        if(mRefinedObjList.GetObj(i).GetClassName()=="Crystal") {
            try {
                Crystal * pCryst = dynamic_cast<Crystal *>(&(mRefinedObjList.GetObj(i)));
                for(int s=0;s<pCryst->GetScattererRegistry().GetNb();s++) {
                    Molecule *pMol=dynamic_cast<Molecule*>(&(pCryst->GetScatt(s)));
                    if(pMol==NULL) continue; // not a Molecule
                    for(vector<MolBond*>::iterator pos = pMol->GetBondList().begin(); pos != pMol->GetBondList().end();++pos) {
                        (*pos)->SetLengthDelta(bdelta);
                        (*pos)->SetLengthSigma(bsigma);
                    }
                    for(vector<MolBondAngle*>::iterator pos=pMol->GetBondAngleList().begin();pos != pMol->GetBondAngleList().end();++pos)
                    {
                        (*pos)->SetAngleDelta(adelta);
                        (*pos)->SetAngleSigma(asigma);
                    }
                }
            } catch (const std::bad_cast& e) {

            }
        }
    }
}
*/
void MonteCarloObj::RunRandomLSQMethod(long &nbCycle)
{
    //perform random move
    mMutationAmplitude=mMutationAmplitudeMax;
    float bsigma=-1, bdelta=-1;
    float asigma=-1, adelta=-1;

    //set the delta and sigma values - low values are good for LSQ!
    for(int i=0;i<mRefinedObjList.GetNb();i++) {
        if(mRefinedObjList.GetObj(i).GetClassName()=="Crystal") {
            try {
                Crystal * pCryst = dynamic_cast<Crystal *>(&(mRefinedObjList.GetObj(i)));
                for(int s=0;s<pCryst->GetScattererRegistry().GetNb();s++)
                {
                    Molecule *pMol=dynamic_cast<Molecule*>(&(pCryst->GetScatt(s)));
                    if(pMol==NULL) continue; // not a Molecule
                    for(vector<MolBond*>::iterator pos = pMol->GetBondList().begin(); pos != pMol->GetBondList().end();++pos) {
                        bsigma = (*pos)->GetLengthSigma();
                        bdelta = (*pos)->GetLengthDelta();
                        (*pos)->SetLengthDelta(0.02);
                        (*pos)->SetLengthSigma(0.001);
                    }
                    for(vector<MolBondAngle*>::iterator pos=pMol->GetBondAngleList().begin();pos != pMol->GetBondAngleList().end();++pos)
                    {
                        asigma = (*pos)->GetAngleSigma();
                        adelta = (*pos)->GetAngleDelta();
                        (*pos)->SetAngleDelta(0.2*DEG2RAD);
                        (*pos)->SetAngleSigma(0.01*DEG2RAD);
                    }
                }
            } catch (const std::bad_cast& e){

            }
        }
    }

    const long starting_point=mRefParList.CreateParamSet("MonteCarloObj:Last parameters (RANDOM-LSQ)");
    mRefParList.SaveParamSet(starting_point);
    while(nbCycle!=0) {
        nbCycle--;
        mRefParList.RestoreParamSet(starting_point);
        //this->NewConfiguration();
        for(int i=0;i<mRefinedObjList.GetNb();i++) mRefinedObjList.GetObj(i).RandomizeConfiguration();
        this->UpdateDisplay();

        //perform LSQ
        for(int i=0;i<mRefinedObjList.GetNb();i++) mRefinedObjList.GetObj(i).SetApproximationFlag(false);
        //mCurrentCost=this->GetLogLikelihood();
        try {
            mLSQ.Refine(20,true,true,false,0.001);
        }
        catch(const ObjCrystException &except) {
            //cout<<"Something wrong?"<<endl;
        };
        for(int i=0;i<mRefinedObjList.GetNb();i++) mRefinedObjList.GetObj(i).SetApproximationFlag(true);
        //cout<<"LSQ cost: "<<mCurrentCost<<" -> "<<this->GetLogLikelihood()<<endl;
        REAL lsq_cost=this->GetLogLikelihood();
        mCurrentCost = lsq_cost;
        //mRefParList.SaveParamSet(lsqtParSavedSetIndex);
        if(mCurrentCost<mBestCost)
        {
            mBestCost=mCurrentCost;
            mRefParList.SaveParamSet(mBestParSavedSetIndex);
        }
        this->UpdateDisplay();

        //save it to the file
        string saveFileName=this->GetName();
        time_t date=time(0);
        char strDate[40];
        strftime(strDate,sizeof(strDate),"%Y-%m-%d_%H-%M-%S",localtime(&date));//%Y-%m-%dT%H:%M:%S%Z
        char costAsChar[30];
        sprintf(costAsChar,"#Run%ld-Cost-%f",nbCycle, mCurrentCost);
        saveFileName=saveFileName+(string)strDate+(string)costAsChar+(string)".xml";
        XMLCrystFileSaveGlobal(saveFileName);

         #ifdef __WX__CRYST__
          mMutexStopAfterCycle.Lock();
          #endif
          if(mStopAfterCycle)
          {
             #ifdef __WX__CRYST__
             mMutexStopAfterCycle.Unlock();
             #endif
             break;
          }
          #ifdef __WX__CRYST__
          mMutexStopAfterCycle.Unlock();
          #endif
    }

    if(bsigma<0 || bdelta<0 || asigma<0 || adelta<0) return;
    //restore the delta and sigma values
    for(int i=0;i<mRefinedObjList.GetNb();i++) {
        if(mRefinedObjList.GetObj(i).GetClassName()=="Crystal") {
            try {
                Crystal * pCryst = dynamic_cast<Crystal *>(&(mRefinedObjList.GetObj(i)));
                for(int s=0;s<pCryst->GetScattererRegistry().GetNb();s++)
                {
                    Molecule *pMol=dynamic_cast<Molecule*>(&(pCryst->GetScatt(s)));
                    if(pMol==NULL) continue; // not a Molecule
                    for(vector<MolBond*>::iterator pos = pMol->GetBondList().begin(); pos != pMol->GetBondList().end();++pos) {
                        (*pos)->SetLengthDelta(bdelta);
                        (*pos)->SetLengthSigma(bsigma);
                    }
                    for(vector<MolBondAngle*>::iterator pos=pMol->GetBondAngleList().begin();pos != pMol->GetBondAngleList().end();++pos)
                    {
                        (*pos)->SetAngleDelta(adelta);
                        (*pos)->SetAngleSigma(asigma);
                    }
                }
            } catch (const std::bad_cast& e){

            }
        }
    }
}

void MonteCarloObj::RunParallelTempering(long &nbStep,const bool silent,
                                         const REAL finalcost,const REAL maxTime)
{
   TAU_PROFILE("MonteCarloObj::RunParallelTempering()","void ()",TAU_DEFAULT);
   TAU_PROFILE_TIMER(timer0a,"MonteCarloObj::RunParallelTempering() Begin 1","", TAU_FIELD);
   TAU_PROFILE_TIMER(timer0b,"MonteCarloObj::RunParallelTempering() Begin 2","", TAU_FIELD);
   TAU_PROFILE_TIMER(timer1,"MonteCarloObj::RunParallelTempering() New Config + LLK","", TAU_FIELD);
   TAU_PROFILE_TIMER(timerN,"MonteCarloObj::RunParallelTempering() Finish","", TAU_FIELD);
   TAU_PROFILE_START(timer0a);
   //Keep a copy of the total number of steps, and decrement nbStep
   const long nbSteps=nbStep;
   unsigned int accept;// 1 if last trial was accepted? 2 if new best config ? else 0
   mNbTrial=0;
   // time (in seconds) when last autoSave was made (if enabled)
      unsigned long secondsWhenAutoSave=0;

   // Periodicity of the automatic LSQ refinements (if the option is set)
   const unsigned int autoLSQPeriod=150000;

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
      for(int i=nbWorld-1;i>=0;i--)
      {
         if((i!=(nbWorld-1))&&(i%2==0))
            for(int j=0;j<mRecursiveRefinedObjList.GetNb();j++)
               mRecursiveRefinedObjList.GetObj(j).RandomizeConfiguration();
         worldCurrentSetIndex(i)=mRefParList.CreateParamSet();
         mRefParList.RestoreParamSet(worldCurrentSetIndex(nbWorld-1));
      }
   TAU_PROFILE_STOP(timer0a);
   TAU_PROFILE_START(timer0b);
      //mNbTrial=nbSteps;;
      const long lastParSavedSetIndex=mRefParList.CreateParamSet("MonteCarloObj:Last parameters (PT)");
      const long runBestIndex=mRefParList.CreateParamSet("Best parameters for current run (PT)");
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
         #if 0
         if(!silent)
            for(int i=0;i<mRefParList.GetNbPar();i++)
            {
               cout << "Gene Group:"<<refParGeneGroupIndex(i)<<" :";
               mRefParList.GetPar(i).Print();
            }
         #endif
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
   // Do we need to update the display ?
   bool needUpdateDisplay=false;
   //Do the refinement
   bool makeReport=false;
   Chronometer chrono;
   chrono.start();
   float lastUpdateDisplayTime=chrono.seconds();
   TAU_PROFILE_STOP(timer0b);
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
            TAU_PROFILE_START(timer1);
            mRefParList.RestoreParamSet(worldCurrentSetIndex(i));
            this->NewConfiguration();
            accept=0;
            REAL cost=this->GetLogLikelihood();
            TAU_PROFILE_STOP(timer1);
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
                  needUpdateDisplay=true;

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
            //if(accept==1 && i==(nbWorld-1)){this->UpdateDisplay();}
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

      if(mAutoLSQ.GetChoice()==2)
         if((mNbTrial%autoLSQPeriod)<(nbTryPerWorld*nbWorld))
         {// Try a quick LSQ ?
            for(int i=0;i<mRefinedObjList.GetNb();i++) mRefinedObjList.GetObj(i).SetApproximationFlag(false);
            for(int i=nbWorld-5;i<nbWorld;i++)
            {
               #ifdef __WX__CRYST__
               mMutexStopAfterCycle.Lock();
               if(mStopAfterCycle)
               {
                  mMutexStopAfterCycle.Unlock();
                  break;
               }
               mMutexStopAfterCycle.Unlock();
               #endif

               mRefParList.RestoreParamSet(worldCurrentSetIndex(i));

               #if 0
               // Report GoF values (Chi^2 / nbObs) values for all objects
               for(map<RefinableObj*,unsigned int>::iterator pos=mLSQ.GetRefinedObjMap().begin();pos!=mLSQ.GetRefinedObjMap().end();++pos)
                  if(pos->first->GetNbLSQFunction()>0)
                  {
                     CrystVector_REAL tmp;
                     tmp =pos->first->GetLSQCalc(pos->second);
                     tmp-=pos->first->GetLSQObs (pos->second);
                     tmp*=tmp;
                     tmp*=pos->first->GetLSQWeight(pos->second);
                     cout<<pos->first->GetClassName()<<":"<<pos->first->GetName()<<": GoF="<<tmp.sum()/tmp.numElements();
                  }
               cout<<endl;
               #endif

               const REAL cost0=this->GetLogLikelihood();// cannot use currentCost(i), approximations changed...
               if(!silent) cout<<"LSQ: World="<<worldSwapIndex(i)<<": cost="<<cost0;
               try {mLSQ.Refine(-30,true,true,false,0.001);}
               catch(const ObjCrystException &except){};
               #if 0
               // Report GoF values (Chi^2 / nbObs) values for all objects
               for(map<RefinableObj*,unsigned int>::iterator pos=mLSQ.GetRefinedObjMap().begin();pos!=mLSQ.GetRefinedObjMap().end();++pos)
                  if(pos->first->GetNbLSQFunction()>0)
                  {
                     CrystVector_REAL tmp;
                     tmp =pos->first->GetLSQCalc(pos->second);
                     tmp-=pos->first->GetLSQObs (pos->second);
                     tmp*=tmp;
                     tmp*=pos->first->GetLSQWeight(pos->second);
                     cout<<pos->first->GetClassName()<<":"<<pos->first->GetName()<<": GoF="<<tmp.sum()/tmp.numElements();
                  }
               cout<<endl;
               #endif
               const REAL cost=this->GetLogLikelihood();
               if(!silent) cout<<" -> "<<cost<<endl;
               if(cost<cost0) mRefParList.SaveParamSet(worldCurrentSetIndex(i));
            }
            //  Need to go back to optimization with approximations allowed (they are not during LSQ)
            for(int i=0;i<mRefinedObjList.GetNb();i++) mRefinedObjList.GetObj(i).SetApproximationFlag(true);
            // And recompute LLK - since they will be lower
            for(int i=nbWorld-5;i<nbWorld;i++)
            {
               mRefParList.RestoreParamSet(worldCurrentSetIndex(i));
               const REAL cost=this->GetLogLikelihood();
               if(!silent) cout<<"LSQ2:"<<currentCost(i)<<"->"<<cost<<endl;
               if(cost<currentCost(i))
               {
                  const REAL oldcost=currentCost(i);
                  mRefParList.SaveParamSet(worldCurrentSetIndex(i));
                  currentCost(i)=cost;
                  if(cost<runBestCost)
                  {
                     runBestCost=currentCost(i);
                     this->TagNewBestConfig();
                     needUpdateDisplay=true;

                     mRefParList.SaveParamSet(runBestIndex);
                     if(runBestCost<mBestCost)
                     {
                        mBestCost=currentCost(i);
                        mRefParList.SaveParamSet(mBestParSavedSetIndex);
                        if(!silent) cout << "->Trial :" << mNbTrial
                                       << " World="<< worldSwapIndex(i)
                                       << " LSQ2: NEW OVERALL Best Cost="<<mBestCost<< endl;
                     }
                     else if(!silent) cout << "->Trial :" << mNbTrial
                                       << " World="<< worldSwapIndex(i)
                                       << " LSQ2: NEW RUN Best Cost="<<runBestCost<< endl;
                     if(!silent) this->DisplayReport();
                  }
                  // KLUDGE - after a successful LSQ, we will be close to a minimum,
                  // which will make most successive global optimization trials to
                  // be rejected, until the temperature is increased a lot - this
                  // is a problem as the temperature increases so much that the
                  // benefit of the LSQ is essentially negated.
                  // So we need to use a higher recorded cost, so that successive trials
                  // may be accepted
                  #if 0
                  mMutationAmplitude=mutationAmplitude(i);
                  for(unsigned int ii=0;ii<4;ii++) this->NewConfiguration(gpRefParTypeObjCryst,false);
                  currentCost(i)=(this->GetLogLikelihood()+cost)/2;
                  if(!silent) cout<<"LSQ3: #"<<worldSwapIndex(i)<<":"<<cost<<"->"<<currentCost(i)<<endl;
                  #else
                  currentCost(i)=oldcost;
                  #endif
               }
            }
         }

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
                     needUpdateDisplay=true;
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
                  simAnnealTemp(i)/=1.5;
               if((worldNbAcceptedMoves(i)/(REAL)nbTrialsReport)>0.95)
                  simAnnealTemp(i)/=1.5;

               if((worldNbAcceptedMoves(i)/(REAL)nbTrialsReport)<0.10)
                  simAnnealTemp(i)*=1.5;
               if((worldNbAcceptedMoves(i)/(REAL)nbTrialsReport)<0.04)
                   simAnnealTemp(i)*=1.5;
               //if((worldNbAcceptedMoves(i)/(REAL)nbTrialsReport)<0.01)
               //   simAnnealTemp(i)*=1.5;
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
      if( (needUpdateDisplay&&(lastUpdateDisplayTime<(chrono.seconds()-1)))||(lastUpdateDisplayTime<(chrono.seconds()-10)))
      {
         mRefParList.RestoreParamSet(runBestIndex);
         this->UpdateDisplay();
         needUpdateDisplay=false;
         lastUpdateDisplayTime=chrono.seconds();
      }
      #ifdef __WX__CRYST__
      mMutexStopAfterCycle.Lock();
      #endif
      if((runBestCost<finalcost) || mStopAfterCycle ||( (maxTime>0)&&(chrono.seconds()>maxTime)))
      {
         #ifdef __WX__CRYST__
         mMutexStopAfterCycle.Unlock();
         #endif
         if(!silent) cout << endl <<endl << "Refinement Stopped:"<<mBestCost<<endl;
         break;
      }
      #ifdef __WX__CRYST__
      mMutexStopAfterCycle.Unlock();
      #endif
   }//Trials

   TAU_PROFILE_START(timerN);
   if(mAutoLSQ.GetChoice()>0)
   {// LSQ
      if(!silent) cout<<"Beginning final LSQ refinement"<<endl;
      for(int i=0;i<mRefinedObjList.GetNb();i++) mRefinedObjList.GetObj(i).SetApproximationFlag(false);
      mRefParList.RestoreParamSet(runBestIndex);
      mCurrentCost=this->GetLogLikelihood();
      try {mLSQ.Refine(-50,true,true,false,0.001);}
      catch(const ObjCrystException &except){};
      if(!silent) cout<<"LSQ cost: "<<mCurrentCost<<" -> "<<this->GetLogLikelihood()<<endl;

      // Need to go back to optimization with approximations allowed (they are not during LSQ)
      for(int i=0;i<mRefinedObjList.GetNb();i++) mRefinedObjList.GetObj(i).SetApproximationFlag(true);

      REAL cost=this->GetLogLikelihood();
      if(cost<mCurrentCost)
      {
         mCurrentCost=cost;
         mRefParList.SaveParamSet(lastParSavedSetIndex);
         if(mCurrentCost<runBestCost)
         {
            runBestCost=mCurrentCost;
            mRefParList.SaveParamSet(runBestIndex);
            if(runBestCost<mBestCost)
            {
               mBestCost=mCurrentCost;
               mRefParList.SaveParamSet(mBestParSavedSetIndex);
               if(!silent) cout << "LSQ : NEW OVERALL Best Cost="<<runBestCost<< endl;
            }
            else if(!silent) cout << " LSQ : NEW Run Best Cost="<<runBestCost<< endl;
         }
      }
      if(!silent) cout<<"Finished LSQ refinement"<<endl;
   }

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
   TAU_PROFILE_STOP(timerN);
}

void MonteCarloObj::XMLOutput(ostream &os,int indent)const
{
   VFN_DEBUG_ENTRY("MonteCarloObj::XMLOutput():"<<this->GetName(),5)
   for(int i=0;i<indent;i++) os << "  " ;
   XMLCrystTag tag("GlobalOptimObj");
   tag.AddAttribute("Name",this->GetName());
   tag.AddAttribute("NbTrialPerRun",(boost::format("%d")%(this->NbTrialPerRun())).str());

   os <<tag<<endl;
   indent++;

   mGlobalOptimType.XMLOutput(os,indent);
   os<<endl;

   mAnnealingScheduleTemp.XMLOutput(os,indent);
   os<<endl;

   mXMLAutoSave.XMLOutput(os,indent);
   os<<endl;

   mAutoLSQ.XMLOutput(os,indent);
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

   mSaveTrackedData.XMLOutput(os,indent);
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
      if("NbTrialPerRun"==tagg.GetAttributeName(i))
      {
         stringstream ss(tagg.GetAttributeValue(i));
         long v;
         ss>>v;
         this->NbTrialPerRun()=v;
      }
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
               if("Save Tracked Data"==tag.GetAttributeValue(i))
               {
                  mSaveTrackedData.XMLInput(is,tag);
                  break;
               }
               if("Automatic Least Squares Refinement"==tag.GetAttributeValue(i))
               {
                  mAutoLSQ.XMLInput(is,tag);
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

LSQNumObj& MonteCarloObj::GetLSQObj() {return mLSQ;}

const LSQNumObj& MonteCarloObj::GetLSQObj() const{return mLSQ;}

void MonteCarloObj::NewConfiguration(const RefParType *type)
{
   TAU_PROFILE("MonteCarloObj::NewConfiguration()","void ()",TAU_DEFAULT);
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

   static string runAutoLSQName;
   static string runAutoLSQChoices[3];

   static string saveTrackedDataName;
   static string saveTrackedDataChoices[2];

   static bool needInitNames=true;
   if(true==needInitNames)
   {
      GlobalOptimTypeName="Algorithm";
      GlobalOptimTypeChoices[0]="Simulated Annealing";
      GlobalOptimTypeChoices[1]="Parallel Tempering";
      //GlobalOptimTypeChoices[2]="Random-LSQ";

      AnnealingScheduleTempName="Temperature Schedule";
      AnnealingScheduleMutationName="Displacement Amplitude Schedule";
      AnnealingScheduleChoices[0]="Constant";
      AnnealingScheduleChoices[1]="Boltzmann";
      AnnealingScheduleChoices[2]="Cauchy";
      AnnealingScheduleChoices[3]="Exponential";
      AnnealingScheduleChoices[4]="Smart";
      AnnealingScheduleChoices[5]="Gamma";

      runAutoLSQName="Automatic Least Squares Refinement";
      runAutoLSQChoices[0]="Never";
      runAutoLSQChoices[1]="At the end of each run";
      runAutoLSQChoices[2]="Every 150000 trials, and at the end of each run";

      saveTrackedDataName="Save Tracked Data";
      saveTrackedDataChoices[0]="No (recommended!)";
      saveTrackedDataChoices[1]="Yes (for tests ONLY)";

      needInitNames=false;//Only once for the class
   }
   mGlobalOptimType.Init(2,&GlobalOptimTypeName,GlobalOptimTypeChoices);
   mAnnealingScheduleTemp.Init(6,&AnnealingScheduleTempName,AnnealingScheduleChoices);
   mAnnealingScheduleMutation.Init(6,&AnnealingScheduleMutationName,AnnealingScheduleChoices);
   mSaveTrackedData.Init(2,&saveTrackedDataName,saveTrackedDataChoices);
   mAutoLSQ.Init(3,&runAutoLSQName,runAutoLSQChoices);
   this->AddOption(&mGlobalOptimType);
   this->AddOption(&mAnnealingScheduleTemp);
   this->AddOption(&mAnnealingScheduleMutation);
   this->AddOption(&mSaveTrackedData);
   this->AddOption(&mAutoLSQ);
   VFN_DEBUG_MESSAGE("MonteCarloObj::InitOptions():End",5)
}

void MonteCarloObj::InitLSQ(const bool useFullPowderPatternProfile)
{
   mLSQ.SetRefinedObj(mRecursiveRefinedObjList.GetObj(0),0,true,true);
   for(unsigned int i=1;i<mRefinedObjList.GetNb();++i)
      mLSQ.SetRefinedObj(mRefinedObjList.GetObj(i),0,false,true);

   if(!useFullPowderPatternProfile)
   {// Use LSQ function #1 for powder patterns (integrated patterns - faster !)
      for(map<RefinableObj*,unsigned int>::iterator pos=mLSQ.GetRefinedObjMap().begin();pos!=mLSQ.GetRefinedObjMap().end();++pos)
         if(pos->first->GetClassName()=="PowderPattern") pos->second=1;
   }
   // Only refine structural parameters (excepting parameters already fixed) and scale factor
   mLSQ.PrepareRefParList(true);
   
   // Intensity corrections can be refined
   std::list<RefinablePar*> vIntCorrPar;
   for(int i=0; i<mLSQ.GetCompiledRefinedObj().GetNbPar();i++)
      if(mLSQ.GetCompiledRefinedObj().GetPar(i).GetType()->IsDescendantFromOrSameAs(gpRefParTypeScattDataCorrInt) && mLSQ.GetCompiledRefinedObj().GetPar(i).IsFixed()==false)
         vIntCorrPar.push_back(&mLSQ.GetCompiledRefinedObj().GetPar(i));
   
   mLSQ.SetParIsFixed(gpRefParTypeScattData,true);
   mLSQ.SetParIsFixed(gpRefParTypeScattDataScale,false);
   
   for(std::list<RefinablePar*>::iterator pos=vIntCorrPar.begin();pos!=vIntCorrPar.end();pos++)
      (*pos)->SetIsFixed(false);
   mLSQ.SetParIsFixed(gpRefParTypeUnitCell,true);
   mLSQ.SetParIsFixed(gpRefParTypeScattPow,true);
   mLSQ.SetParIsFixed(gpRefParTypeRadiation,true);
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
