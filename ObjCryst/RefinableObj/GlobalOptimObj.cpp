/*
* LibCryst++ : a Crystallographic computing library in C++
*
*  (c) 2000 Vincent FAVRE-NICOLIN
*           Laboratoire de Cristallographie
*           24, quai Ernest-Ansermet, CH-1211 Geneva 4, Switzerland
*  Contact: Vincent.Favre-Nicolin@cryst.unige.ch
*           Radovan.Cerny@cryst.unige.ch
*
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
#ifdef __WX__CRYST__
   #include "wxCryst/wxRefinableObj.h"
   #undef GetClassName // Conflict from wxMSW headers ? (cygwin)
#endif
#include <fstream>
#include <sstream>

namespace ObjCryst
{
//#################################################################################
//
//       GlobalOptimObj
//
//#################################################################################
ObjRegistry<GlobalOptimObj> gGlobalOptimObjRegistry("List of all Global Optimization objects");


GlobalOptimObj::GlobalOptimObj(const string name):
mName(name),mSaveFileName("GlobalOptim.save"),
mNbTrial(0),mCurrentCost(-1),mBestCost(-1),
mHistoryNb(0),mHistoryTrialNumber(1000),mHistoryCostFunction(1000),
mHistorySavedParamSetIndex(1000),
mHistorySaveAfterEachOptim(true),mHistorySaveFileName("GlobalOptim_history.out"),
mBestParSavedSetIndex(-1),mLastParSavedSetIndex(-1),
mTemperatureMax(.1),mTemperatureMin(.001),
mMutationAmplitudeBegin(1.),mMutationAmplitudeEnd(1.),
mNbTrialRetry(0),mMinCostRetry(0.),mMaxNbTrialSinceBest(0),
mIsOptimizing(false),mStopAfterCycle(false),
mRefinedObjList("GlobalOptimObj: "+mName+" RefinableObj registry"),
mRecursiveRefinedObjList("GlobalOptimObj: "+mName+" recursive RefinableObj registry"),
mNbCostFunction(0),mMaxNbCostFunction(20),
mpCostFunctionId(mMaxNbCostFunction),mCostFunctionWeight(mMaxNbCostFunction)
#ifdef __WX__CRYST__
,mpWXGlobalOptimObj(0)
#endif

{
   VFN_DEBUG_MESSAGE("GlobalOptimObj::GlobalOptimObj()",5)
   this->InitOptions();
   mGlobalOptimType.SetChoice(GLOBAL_OPTIM_PARALLEL_TEMPERING);
   mAnnealingScheduleTemp.SetChoice(ANNEALING_EXPONENTIAL);
   mAnnealingScheduleMutation.SetChoice(ANNEALING_SMART);
   gGlobalOptimObjRegistry.Register(*this);
}

GlobalOptimObj::~GlobalOptimObj()
{
   gGlobalOptimObjRegistry.DeRegister(*this);
}

void GlobalOptimObj::RandomizeStartingConfig()
{
   VFN_DEBUG_MESSAGE("GlobalOptimObj::RandomizeStartingConfig()",5)
   this->PrepareRefParList();
   this->InitRandomSeedFromTime();
   for(int j=0;j<mRefParList.GetNbParNotFixed();j++)
   {
      if(true==mRefParList.GetParNotFixed(j).IsLimited())
      {
         const double min=mRefParList.GetParNotFixed(j).GetMin();
         const double max=mRefParList.GetParNotFixed(j).GetMax();
         mRefParList.GetParNotFixed(j).MutateTo(min+(max-min)*(rand()/(double)RAND_MAX) );
      }
   }
      //else cout << mRefParList.GetParNotFixed(j).Name() <<" Not limited :-(" <<endl;
   VFN_DEBUG_MESSAGE("GlobalOptimObj::RandomizeStartingConfig():Finished",5)
}
void GlobalOptimObj::SetAlgorithmSimulAnnealing(const AnnealingSchedule scheduleTemp,
                           const double tMax, const double tMin,
                           const AnnealingSchedule scheduleMutation,
                           const double mutMax, const double mutMin,
                           const long nbTrialRetry,const double minCostRetry,
                           const long maxNbTrialSinceBest)
{
   VFN_DEBUG_MESSAGE("GlobalOptimObj::SetAlgorithmSimulAnnealing()",5)
   if(mAnnealingScheduleTemp.GetChoice()==ANNEALING_SMART)
      throw ObjCrystException("GlobalOptimObj::SetAlgorithmSimulAnnealing() : \
Cannot use ANNEALING_SMART for the Temperature schedule (yet).");
   mGlobalOptimType.SetChoice(GLOBAL_OPTIM_SIMULATED_ANNEALING);
   mTemperatureMax=tMax;
   mTemperatureMin=tMin;
   mAnnealingScheduleTemp.SetChoice(scheduleTemp);


   mMutationAmplitudeBegin=mutMax;
   mMutationAmplitudeEnd=mutMin;
   mAnnealingScheduleMutation.SetChoice(scheduleMutation);
   mNbTrialRetry=nbTrialRetry;
   mMinCostRetry=minCostRetry;
   mMaxNbTrialSinceBest=maxNbTrialSinceBest;
   VFN_DEBUG_MESSAGE("GlobalOptimObj::SetAlgorithmSimulAnnealing():End",3)
}

void GlobalOptimObj::SetAlgorithmParallTempering(const AnnealingSchedule scheduleTemp,
                                 const double tMax, const double tMin,
                                 const AnnealingSchedule scheduleMutation,
                                 const double mutMax, const double mutMin)
{
   VFN_DEBUG_MESSAGE("GlobalOptimObj::SetAlgorithmParallTempering()",5)
   if(mAnnealingScheduleTemp.GetChoice()==ANNEALING_SMART)
      throw ObjCrystException("GlobalOptimObj::SetAlgorithmParallTempering() : \
Cannot use ANNEALING_SMART for the Temperature schedule (yet).");
   mGlobalOptimType.SetChoice(GLOBAL_OPTIM_PARALLEL_TEMPERING);
   mTemperatureMax=tMax;
   mTemperatureMin=tMin;
   mAnnealingScheduleTemp.SetChoice(scheduleTemp);

   mMutationAmplitudeBegin=mutMax;
   mMutationAmplitudeEnd=mutMin;
   mAnnealingScheduleMutation.SetChoice(scheduleMutation);
   //mNbTrialRetry=nbTrialRetry;
   //mMinCostRetry=minCostRetry;
   //mMaxNbTrialSinceBest=maxNbTrialSinceBest;
   VFN_DEBUG_MESSAGE("GlobalOptimObj::SetAlgorithmParallTempering():End",3)
}

void GlobalOptimObj::Optimize(long &nbStep)
{
   //Keep a copy of the total number of steps, and decrement nbStep
   const long nbSteps=nbStep;
   //:TODO: Other algorithms !
   TAU_PROFILE("GlobalOptimObj::Optimize()","void (long)",TAU_DEFAULT);
   VFN_DEBUG_MESSAGE("GlobalOptimObj::Optimize()",5)
   for(int i=0;i<mRefinedObjList.GetNb();i++) mRefinedObjList.GetObj(i).BeginOptimization();
   this->PrepareRefParList();
   mRefParList.Print();
   this->InitRandomSeedFromTime();
   mNbTrial=0;
   mRefParList.EraseAllParamSet();
   mIsOptimizing=true;
   // prepare all objects
   
   /*
   if(mBestCost<.3)
   {//:KLUDGE: :TODO: Remove this and do it interactively !!
      this->UnFixPar("Width_CagliotiU");
      this->UnFixPar("Width_CagliotiV");
      this->UnFixPar("Width_CagliotiW");
      this->UnFixPar("PseudoVoigt_Eta0");
   }
   */            
   switch(mGlobalOptimType.GetChoice())
   {
      case GLOBAL_OPTIM_SIMULATED_ANNEALING:
      {
         RESTART_OPTIMIZATION:

         //Re-init History
            mHistoryNb=0;
            mHistoryTrialNumber=-1;
            mHistoryCostFunction=-1;
            mHistorySavedParamSetIndex=-1;
         //Save Starting point
            mHistoryTrialNumber(mHistoryNb)=mNbTrial;
            mHistoryCostFunction(mHistoryNb)=this->GetCostFunctionValue();
            mHistorySavedParamSetIndex(mHistoryNb)=mRefParList.CreateParamSet();
            mHistoryNb++;

         mCurrentCost=this->GetCostFunctionValue();

         mBestCost=mCurrentCost;
         mBestParSavedSetIndex=mRefParList.CreateParamSet("GlobalOptim:Best parameters");
         mLastParSavedSetIndex=mRefParList.CreateParamSet("GlobalOptim:Last parameters");

         //Report each ... cycles
            const int nbTryReport=1000;
         // Keep record of the number of accepted moves
            long nbAcceptedMoves=0;//since last report
            long nbAcceptedMovesTemp=0;//since last temperature/mutation rate change
         // Number of tries since best configuration found
            long nbTriesSinceBest=0;
         // Change temperature (and mutation) every...
            const int nbTryPerTemp=100;

         double simAnnealTemp=-1;
         mMutationAmplitude=1.;

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
                     simAnnealTemp=
                        mTemperatureMin*log((double)nbSteps)/log((double)(mNbTrial+1));break;
                  case ANNEALING_CAUCHY:
                     simAnnealTemp=mTemperatureMin*nbSteps/mNbTrial;break;
                  //case ANNEALING_QUENCHING:
                  case ANNEALING_EXPONENTIAL:
                     simAnnealTemp=mTemperatureMax
                                    *pow(mTemperatureMin/mTemperatureMax,
                                          mNbTrial/(double)nbSteps);break;
                  case ANNEALING_SMART:break;//:TODO:
                  default: simAnnealTemp=mTemperatureMin;break;
               }
               switch(mAnnealingScheduleMutation.GetChoice())
               {
                  case ANNEALING_BOLTZMANN:
                     mMutationAmplitude=
                        mMutationAmplitudeEnd*log((double)nbSteps)/log((double)(mNbTrial+1));
                     break;
                  case ANNEALING_CAUCHY:
                     mMutationAmplitude=mMutationAmplitudeEnd*nbSteps/mNbTrial;break;
                  //case ANNEALING_QUENCHING:
                  case ANNEALING_EXPONENTIAL:
                     mMutationAmplitude=mMutationAmplitudeBegin
                                    *pow(mMutationAmplitudeEnd/mMutationAmplitudeBegin,
                                          mNbTrial/(double)nbSteps);break;
                  case ANNEALING_SMART:
                     if((nbAcceptedMovesTemp/(double)nbTryPerTemp)>0.7) mMutationAmplitude*=2.;
                     if((nbAcceptedMovesTemp/(double)nbTryPerTemp)<0.3) mMutationAmplitude/=2.;
                     if(mMutationAmplitude>16.) mMutationAmplitude=16.;
                     if(mMutationAmplitude<(1./16.)) mMutationAmplitude=1./16.;
                     nbAcceptedMovesTemp=0;
                     break;
                  default: mMutationAmplitude=mMutationAmplitudeEnd;break;
               }
            }

            this->NewConfiguration();
            double cost=this->GetCostFunctionValue();
            if(cost<mCurrentCost)
            {
               mCurrentCost=cost;
               mRefParList.SaveParamSet(mLastParSavedSetIndex);
               if(mCurrentCost<mBestCost)
               {
                  mBestCost=mCurrentCost;
                  mRefParList.SaveParamSet(mBestParSavedSetIndex);
                  mHistoryTrialNumber(mHistoryNb)=mNbTrial;
                  mHistoryCostFunction(mHistoryNb)=mCurrentCost;
                  mHistorySavedParamSetIndex(mHistoryNb)=mRefParList.CreateParamSet();
                  mHistoryNb++;
                  cout << "Trial :" << mNbTrial 
                       << " Temp="<< simAnnealTemp
                       << " Mutation Ampl.: "<<mMutationAmplitude
                       << " Best Cost="<<mBestCost<< endl;
                  nbTriesSinceBest=0;
                  this->UpdateDisplay();
               }
               nbAcceptedMoves++;
               nbAcceptedMovesTemp++;
            }
            else
            {
               if( (rand()/(double)RAND_MAX) < exp(-(cost-mCurrentCost)/simAnnealTemp) )
               {
                  mCurrentCost=cost;
                  mRefParList.SaveParamSet(mLastParSavedSetIndex);
                  nbAcceptedMoves++;
                  nbAcceptedMovesTemp++;
               }
               else mRefParList.RestoreParamSet(mLastParSavedSetIndex);
            }
            if( (mNbTrial % nbTryReport) == 0)
            {
               cout <<"Trial :" << mNbTrial << " Temp="<< simAnnealTemp;
               cout <<" Mutation Ampl.: " <<mMutationAmplitude<< " Best Cost=" << mBestCost 
                    <<" Current Cost=" << mCurrentCost 
                    <<" Accepting "<<(int)((double)nbAcceptedMoves/nbTryReport*100)
                    <<"% moves" << endl;
               nbAcceptedMoves=0;
               #ifdef __WX__CRYST__
               if(0!=mpWXGlobalOptimObj) mpWXGlobalOptimObj->UpdateDisplayNbTrial();
               #endif
            }
            mNbTrial++;nbStep--;
            
            if(mStopAfterCycle) 
            {
               cout << endl <<endl << "Refinement Stopped."<<endl;
               break;
            }
            
            if((mMaxNbTrialSinceBest>0) && (nbTriesSinceBest>mMaxNbTrialSinceBest))
            {
               cout << nbTriesSinceBest
                    <<" trials since best configuration. Going back to Best config."<<endl;
               mRefParList.RestoreParamSet(mBestParSavedSetIndex);
               nbTriesSinceBest=0;
            }
            if((mNbTrialRetry>0) && (mCurrentCost>mMinCostRetry) &&(mNbTrial>mNbTrialRetry))
            {
               cout << mNbTrial
                    <<" trials done and still above" << mMinCostRetry
                    <<"Randomizing and restarting"<<endl;
               goto RESTART_OPTIMIZATION;
            }
            nbTriesSinceBest++;
         }
         //Restore Best values
         mRefParList.RestoreParamSet(mBestParSavedSetIndex);
         mCurrentCost=this->GetCostFunctionValue();
         chrono.print();
         //mRefParList.Print();
         if(true==mHistorySaveAfterEachOptim) this->SaveOptimHistory();
      }//case GLOBAL_OPTIM_SIMULATED_ANNEALING
      case GLOBAL_OPTIM_PARALLEL_TEMPERING:
      {
         //Total number of parallel refinements,each is a 'World'. The most stable
         // world must be i=nbWorld-1, and the most changing World (high mutation,
         // high temperature) is i=0.
            const long nbWorld=10;
         // Init the different temperatures
            CrystVector_double simAnnealTemp(nbWorld);
            for(int i=0;i<nbWorld;i++)
            {
               switch(mAnnealingScheduleTemp.GetChoice())
               {
                  case ANNEALING_BOLTZMANN:
                     simAnnealTemp(i)=
                        mTemperatureMin*log((double)nbWorld)/log((double)(i+1));break;
                  case ANNEALING_CAUCHY:
                     simAnnealTemp(i)=mTemperatureMin*nbWorld/(i+1);break;
                  //case ANNEALING_QUENCHING:
                  case ANNEALING_EXPONENTIAL:
                     simAnnealTemp(i)=mTemperatureMax
                                    *pow(mTemperatureMin/mTemperatureMax,
                                          i/(double)(nbWorld-1));break;
                  case ANNEALING_SMART:break;//:TODO:
                  default: simAnnealTemp(i)=mTemperatureMin;break;
               }
            }
         //Init the different mutation rate parameters
            CrystVector_double mutationAmplitude(nbWorld);
            for(int i=0;i<nbWorld;i++)
            {
               switch(mAnnealingScheduleMutation.GetChoice())
               {
                  case ANNEALING_BOLTZMANN:
                     mutationAmplitude(i)=
                        mMutationAmplitudeEnd*log((double)(nbWorld-1))/log((double)(i+1));
                     break;
                  case ANNEALING_CAUCHY:
                     mutationAmplitude(i)=mMutationAmplitudeEnd*(double)(nbWorld-1)/i;break;
                  //case ANNEALING_QUENCHING:
                  case ANNEALING_EXPONENTIAL:
                     mutationAmplitude(i)=mMutationAmplitudeBegin
                                    *pow(mMutationAmplitudeEnd/mMutationAmplitudeBegin,
                                          i/(double)(nbWorld-1));break;
                  case ANNEALING_SMART:mutationAmplitude(i)=1.;//will be updated later
                  default: mMutationAmplitude=mMutationAmplitudeEnd;break;
               }
            }
         // Number of successive trials for each World. At the end of these trials
         // a swap is tried with the upper World (eg i-1). This number effectvely sets
         // the rate of swapping.
            const int nbTryPerWorld=10;
         // Initialize the costs
            mCurrentCost=this->GetCostFunctionValue();
            mBestCost=mCurrentCost;
            CrystVector_double currentCost(nbWorld);
            currentCost=mCurrentCost;
         // Init the parameter sets for each World
         // All Worlds start from the same (current) configuration.
            CrystVector_long worldCurrentSetIndex(nbWorld);
            for(int i=0;i<nbWorld;i++)
               worldCurrentSetIndex(i)=mRefParList.CreateParamSet();
            mBestParSavedSetIndex=mRefParList.CreateParamSet("GlobalOptim:Best parameters");
            mLastParSavedSetIndex=mRefParList.CreateParamSet("GlobalOptim:Last parameters");
            const long swapParSavedSetIndex=mRefParList.CreateParamSet("Tmp Par Set");
         //Keep track of how many trials are accepted for each World
            CrystVector_long worldNbAcceptedMoves(nbWorld);
            worldNbAcceptedMoves=0;
         //Do a report each... And check if mutation rate is OK (for annealing_smart)s
            const int nbTrialsReport=1000;
         //Do the refinement
         bool makeReport=false;
         Chronometer chrono;
         chrono.start();
         for(;mNbTrial<nbSteps;)
         {
            for(int i=0;i<nbWorld;i++)
            {
               mRefParList.RestoreParamSet(worldCurrentSetIndex(i));
               mMutationAmplitude=mutationAmplitude(i);
               for(int j=0;j<nbTryPerWorld;j++)
               {
                  mRefParList.SaveParamSet(mLastParSavedSetIndex);
                  this->NewConfiguration();
                  double cost=this->GetCostFunctionValue();
                  if(cost<currentCost(i))
                  {
                     currentCost(i)=cost;
                     mRefParList.SaveParamSet(worldCurrentSetIndex(i));
                     if(cost<mBestCost)
                     {
                        mBestCost=currentCost(i);
                        mRefParList.SaveParamSet(mBestParSavedSetIndex);
                        cout << "->Trial :" << mNbTrial 
                             << " World="<< i
                             << " Temp="<< simAnnealTemp(i)
                             << " Mutation Ampl.: "<<mMutationAmplitude
                             << " NEW Best Cost="<<mBestCost<< endl;
                        //nbTriesSinceBest=0;
                        this->DisplayReport();
                        this->UpdateDisplay();
                        //for(int i=0;i<mRefinedObjList.GetNb();i++) 
                        //   mRefinedObjList.GetObj(i).Print();
                     }
                     worldNbAcceptedMoves(i)++;
                  }
                  else
                  {
                     if((rand()/(double)RAND_MAX)<exp(-(cost-currentCost(i))/simAnnealTemp(i)) )
                     {
                        currentCost(i)=cost;
                        mRefParList.SaveParamSet(worldCurrentSetIndex(i));
                        worldNbAcceptedMoves(i)++;
                     }
                     else mRefParList.RestoreParamSet(mLastParSavedSetIndex);
                  }
                  mNbTrial++;nbStep--;
                  if((mNbTrial%nbTrialsReport)==0) makeReport=true;
               }//nbTryPerWorld trials
            }//For each World
            
            //Try swapping worlds
            for(int i=1;i<nbWorld;i++)
            {
               
               if((rand()/(double)RAND_MAX)
                      < exp(-(currentCost(i-1)-currentCost(i))/simAnnealTemp(i)))
               {  
               /*
                  if(i>2)
                  {
                     cout <<"->Swapping Worlds :" << i <<"(cost="<<currentCost(i)<<")"
                          <<" with "<< (i-1) <<"(cost="<< currentCost(i-1)<<")"<<endl;
                  }
                  */
                  mRefParList.RestoreParamSet(worldCurrentSetIndex(i));
                  mRefParList.SaveParamSet(swapParSavedSetIndex);
                  mRefParList.RestoreParamSet(worldCurrentSetIndex(i-1));
                  mRefParList.SaveParamSet(worldCurrentSetIndex(i));
                  mRefParList.RestoreParamSet(swapParSavedSetIndex);
                  mRefParList.SaveParamSet(worldCurrentSetIndex(i-1));
                  const double tmp=currentCost(i);
                  currentCost(i)=currentCost(i-1);
                  currentCost(i-1)=tmp;
               }
            }
            
            if(true==makeReport)
            {
               makeReport=false;
               worldNbAcceptedMoves*=nbWorld;
               cout <<"Trial :" << mNbTrial << " Best Cost=" << mBestCost<< " ";
               chrono.print();
               for(int i=0;i<nbWorld;i++)
               {
                  cout <<"   World :" << i
                       <<" Temp.: " << simAnnealTemp(i)
                       <<" Mutation Ampl.: " << mutationAmplitude(i)
                       <<" Current Cost=" << currentCost(i)
                       <<" Accepting "
                       << (int)((double)worldNbAcceptedMoves(i)/nbTrialsReport*100)
                       <<"% moves" << endl;
               }
               //Change the mutation rate if necessary for each world
               if(ANNEALING_SMART==mAnnealingScheduleMutation.GetChoice())
               {
                  for(int i=0;i<nbWorld;i++)
                  {
                     if((worldNbAcceptedMoves(i)/(double)nbTrialsReport)>0.7)
                        mutationAmplitude(i)*=2.;
                     if((worldNbAcceptedMoves(i)/(double)nbTrialsReport)<0.3)
                        mutationAmplitude(i)/=2.;
                     if(mutationAmplitude(i)>65.) mutationAmplitude(i)=64.;
                     if(mutationAmplitude(i)<(1./9.)) mutationAmplitude(i)=1./8.;
                  }
               }
               worldNbAcceptedMoves=0;
               //this->DisplayReport();
               #ifdef __WX__CRYST__
               if(0!=mpWXGlobalOptimObj) mpWXGlobalOptimObj->UpdateDisplayNbTrial();
               #endif
            }
            if(mStopAfterCycle) 
            {
               cout << endl <<endl << "Refinement Stopped:"<<mBestCost<<endl;
               break;
            }
         }//Trials
         //Restore Best values
            //mRefParList.Print();
            this->DisplayReport();
            mRefParList.RestoreParamSet(mBestParSavedSetIndex);
            //for(int i=0;i<mRefinedObjList.GetNb();i++) mRefinedObjList.GetObj(i).Print();
            mCurrentCost=this->GetCostFunctionValue();
            for(unsigned int i=0;i<mNbCostFunction;i++)
               cout <<mpCostFunctionRefinableObj[i]->GetClassName()<<":"
                    <<mpCostFunctionRefinableObj[i]->GetName()<<":"
                    <<mpCostFunctionRefinableObj[i]
                        ->GetCostFunctionName(mpCostFunctionId(i))<<"="
                    <<mpCostFunctionRefinableObj[i]
                        ->GetCostFunctionValue(mpCostFunctionId(i))
                    <<", weight="<<mCostFunctionWeight(i)<<endl;
            cout<<"Overall cost:"<<mCurrentCost<<"("<<mBestCost<<")"<<endl;
            chrono.print();

      }//case GLOBAL_OPTIM_PARALLEL_TEMPERING
      case GLOBAL_OPTIM_GENETIC: //:TODO:
      {
      }//case GLOBAL_OPTIM_GENETIC
   }
   mIsOptimizing=false;
   mStopAfterCycle=false;
   for(int i=0;i<mRefinedObjList.GetNb();i++) mRefinedObjList.GetObj(i).EndOptimization();
   VFN_DEBUG_MESSAGE("GlobalOptimObj::Optimize():End.",5)
}

//RefinablePar& GlobalOptimObj::GetPar(const string& parName)
//{
//   this->PrepareRefParList();//:TODO: Do not reset the list every time...
//   return mRefParList.GetPar(parName);
//}

void GlobalOptimObj::FixAllPar()
{
   VFN_DEBUG_MESSAGE("GlobalOptimObj::FixAllPar()",5)
   for(int i=0;i<mRecursiveRefinedObjList.GetNb();i++) 
      mRecursiveRefinedObjList.GetObj(i).FixAllPar();
   VFN_DEBUG_MESSAGE("GlobalOptimObj::FixAllPar():End",5)
}
void GlobalOptimObj::SetParIsFixed(const string& parName,const bool fix)
{
   for(int i=0;i<mRecursiveRefinedObjList.GetNb();i++) 
      mRecursiveRefinedObjList.GetObj(i).SetParIsFixed(parName,fix);
}
void GlobalOptimObj::SetParIsFixed(const RefParType *type,const bool fix)
{
   for(int i=0;i<mRecursiveRefinedObjList.GetNb();i++) 
      mRecursiveRefinedObjList.GetObj(i).SetParIsFixed(type,fix);
}
   
void GlobalOptimObj::UnFixAllPar()
{
   for(int i=0;i<mRecursiveRefinedObjList.GetNb();i++) 
      mRecursiveRefinedObjList.GetObj(i).UnFixAllPar();}
   
void GlobalOptimObj::SetParIsUsed(const string& parName,const bool use)
{
   for(int i=0;i<mRecursiveRefinedObjList.GetNb();i++) 
      mRecursiveRefinedObjList.GetObj(i).SetParIsUsed(parName,use);
}
void GlobalOptimObj::SetParIsUsed(const RefParType *type,const bool use)
{
   for(int i=0;i<mRecursiveRefinedObjList.GetNb();i++) 
      mRecursiveRefinedObjList.GetObj(i).SetParIsUsed(type,use);
}
void GlobalOptimObj::SetLimitsRelative(const string &parName,
                                       const double min, const double max)
{
   for(int i=0;i<mRecursiveRefinedObjList.GetNb();i++) 
      mRecursiveRefinedObjList.GetObj(i).SetLimitsRelative(parName,min,max);
}
void GlobalOptimObj::SetLimitsRelative(const RefParType *type,
                                       const double min, const double max)
{
   for(int i=0;i<mRecursiveRefinedObjList.GetNb();i++) 
      mRecursiveRefinedObjList.GetObj(i).SetLimitsRelative(type,min,max);
}
void GlobalOptimObj::SetLimitsAbsolute(const string &parName,
                                       const double min, const double max)
{
   for(int i=0;i<mRecursiveRefinedObjList.GetNb();i++) 
      mRecursiveRefinedObjList.GetObj(i).SetLimitsAbsolute(parName,min,max);
}
void GlobalOptimObj::SetLimitsAbsolute(const RefParType *type,
                                       const double min, const double max)
{
   for(int i=0;i<mRecursiveRefinedObjList.GetNb();i++) 
      mRecursiveRefinedObjList.GetObj(i).SetLimitsAbsolute(type,min,max);
}

void GlobalOptimObj::SaveOptimHistory() const
{
   VFN_DEBUG_MESSAGE("GlobalOptimObj::SaveOptimHistory()",5)
   if(mHistoryNb<=0)
   {
      cout << "GlobalOptimObj::SaveOptimHistory(): No History yet !!! "<<endl;
      return;
   }
   
   ofstream out(mHistorySaveFileName.c_str());
   if(!out)
   {
      throw ObjCrystException("GlobalOptimObj::SaveOptimHistory() : \
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

double GlobalOptimObj::GetCostFunctionValue() 
{
   double cost =0.;
   for(unsigned int i=0;i<mNbCostFunction;i++)
		if(mCostFunctionWeight(i)>0)
      	cost += mCostFunctionWeight(i)*mpCostFunctionRefinableObj[i]
						->GetCostFunctionValue(mpCostFunctionId(i));
   return cost;
}
void GlobalOptimObj::StopAfterCycle() 
{
   VFN_DEBUG_MESSAGE("GlobalOptimObj::StopAfterCycle()",5)
   if(mIsOptimizing) mStopAfterCycle=true;
}

void GlobalOptimObj::DisplayReport() 
{
   for(unsigned int i=0;i<mNbCostFunction;i++)
      cout <<mpCostFunctionRefinableObj[i]->GetClassName()<<":"
           <<mpCostFunctionRefinableObj[i]->GetName()<<":"
           <<mpCostFunctionRefinableObj[i]
               ->GetCostFunctionName(mpCostFunctionId(i))<<"="
           <<mpCostFunctionRefinableObj[i]
               ->GetCostFunctionValue(mpCostFunctionId(i))
           <<", weight="<<mCostFunctionWeight(i)<<endl;
}

void GlobalOptimObj::AddRefinableObj(RefinableObj &obj)
{
   VFN_DEBUG_MESSAGE("GlobalOptimObj::AddRefinableObj():"<<obj.GetName(),5)
   mRefinedObjList.Register(obj);
   RefObjRegisterRecursive(obj,mRecursiveRefinedObjList);
   #ifdef __WX__CRYST__
   if(0!=mpWXGlobalOptimObj) mpWXGlobalOptimObj->AddRefinedObject(obj);
   #endif
}

void GlobalOptimObj::AddCostFunction(RefinableObj &obj,const unsigned int id, const double w)
{
   VFN_DEBUG_MESSAGE("GlobalOptimObj::AddGetCostFunctionValue()",5)
   mpCostFunctionRefinableObj[mNbCostFunction]=&obj;
   mpCostFunctionId(mNbCostFunction)=id;
   mCostFunctionWeight(mNbCostFunction++)=w;
   
   //Forces the object to get ready
      obj.BeginOptimization();
      obj.EndOptimization();
   #ifdef __WX__CRYST__
   if(0!=mpWXGlobalOptimObj) mpWXGlobalOptimObj->AddCostFunction(obj,id);
   #endif
}

void GlobalOptimObj::Output(ostream &os,int indent)const
{
   VFN_DEBUG_ENTRY("GlobalOptimObj::Output():"<<this->GetName(),5)
   for(int i=0;i<indent;i++) os << "  " ;
   XMLCrystTag tag("GlobalOptimObj");
   tag.AddAttribute("Name",this->GetName());
   
   os <<tag<<endl;
   indent++;
   
   mGlobalOptimType.Output(os,indent);
   os<<endl;

   mAnnealingScheduleTemp.Output(os,indent);
   os<<endl;
   {
      XMLCrystTag tag2("TempMaxMin");
      for(int i=0;i<indent;i++) os << "  " ;
      os<<tag2<<mTemperatureMax << " "<< mTemperatureMin;
      tag2.SetIsEndTag(true);
      os<<tag2<<endl;
   }

   mAnnealingScheduleMutation.Output(os,indent);
   os<<endl;
   
   {
      XMLCrystTag tag2("MutationMaxMin");
      for(int i=0;i<indent;i++) os << "  " ;
      os<<tag2<<mMutationAmplitudeBegin << " "<< mMutationAmplitudeEnd;
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
   
   for(unsigned int j=0;j<mNbCostFunction;j++)
   {
      XMLCrystTag tag2("CostFunction",false,true);
      tag2.AddAttribute("ObjectType",mpCostFunctionRefinableObj[j]->GetClassName());
      tag2.AddAttribute("ObjectName",mpCostFunctionRefinableObj[j]->GetName());
      stringstream ss;
      ss<<mpCostFunctionId(j);
      tag2.AddAttribute("FunctionId",ss.str());
      tag2.AddAttribute("FunctionName",
                        mpCostFunctionRefinableObj[j]
                           ->GetCostFunctionDescription(mpCostFunctionId(j)));
      stringstream ss2;
      ss2<<mCostFunctionWeight(j);
      tag2.AddAttribute("Weight",ss2.str());
      for(int i=0;i<indent;i++) os << "  " ;
      os<<tag2<<endl;
   }
   
   indent--;
   tag.SetIsEndTag(true);
   for(int i=0;i<indent;i++) os << "  " ;
   os <<tag<<endl;
   VFN_DEBUG_EXIT("GlobalOptimObj::Output():"<<this->GetName(),5)
}

void GlobalOptimObj::Input(istream &is,const XMLCrystTag &tagg)
{
   VFN_DEBUG_ENTRY("GlobalOptimObj::Input():"<<this->GetName(),5)
   for(unsigned int i=0;i<tagg.GetNbAttribute();i++)
   {
      if("Name"==tagg.GetAttributeName(i)) this->SetName(tagg.GetAttributeValue(i));
   }
   while(true)
   {
      XMLCrystTag tag(is);
      if(("GlobalOptimObj"==tag.GetName())&&tag.IsEndTag())
      {
         VFN_DEBUG_EXIT("GlobalOptimObj::Exit():"<<this->GetName(),5)
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
                  mGlobalOptimType.Input(is,tag);
                  break;
               }
               if("Temperature Schedule"==tag.GetAttributeValue(i))
               {
                  mAnnealingScheduleTemp.Input(is,tag);
                  break;
               }
               if("Displacement Amplitude Schedule"==tag.GetAttributeValue(i))
               {
                  mAnnealingScheduleMutation.Input(is,tag);
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
         is>>mMutationAmplitudeBegin>>mMutationAmplitudeEnd;
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
      if("CostFunction"==tag.GetName())
      {
         string name,type;
         float weight;
         int id;
         for(unsigned int i=0;i<tag.GetNbAttribute();i++)
         {
            if("ObjectName"==tag.GetAttributeName(i)) name=tag.GetAttributeValue(i);
            if("ObjectType"==tag.GetAttributeName(i)) type=tag.GetAttributeValue(i);
            if("Weight"==tag.GetAttributeName(i))
            {
               stringstream ss(tag.GetAttributeValue(i));
               ss>>weight;
            }
            if("FunctionId"==tag.GetAttributeName(i))
            {
               stringstream ss(tag.GetAttributeValue(i));
               ss>>id;
            }
         }
         RefinableObj* obj=& (gRefinableObjRegistry.GetObj(name,type));
         this->AddCostFunction(*obj,id,weight);
      }
   }
}

void GlobalOptimObj::InputOld(istream &is,const IOCrystTag &tagg)
{
   VFN_DEBUG_MESSAGE("GlobalOptimObj::Input():"<<this->GetName(),5)
   switch(tagg.GetVersion())
   {
      case 0:
      {
         this->SetName(tagg.GetName());
         do
         {
            IOCrystTag tag(is);
            if(tag.IsClosingTag()==true)
            {
               VFN_DEBUG_MESSAGE("GlobalOptimObj::Input():End",5)
               return;
            }
            if(tag.GetType()=="Param")
            {
               if(tag.GetName()=="Algorithm")
               {
                  string alg;
                  IOCrystExtractNameSpace(is,alg);
                  if(alg=="SimulatedAnnealing")
                  {
                     mGlobalOptimType.SetChoice(GLOBAL_OPTIM_SIMULATED_ANNEALING);
                     continue;
                  }
                  if(alg=="ParallelTempering")
                  {
                     mGlobalOptimType.SetChoice(GLOBAL_OPTIM_PARALLEL_TEMPERING);
                     continue;
                  }
                  if(alg=="Genetic")
                  {
                     mGlobalOptimType.SetChoice(GLOBAL_OPTIM_GENETIC);
                     continue;
                  }
                  continue;
               }
               if(tag.GetName()=="TemperatureScheduleMaxMin")
               {
                  string temp;
                  IOCrystExtractNameSpace(is,temp);
                  is >> mTemperatureMax>>mTemperatureMin;
                  if(temp=="Constant")
                  {
                     mAnnealingScheduleTemp.SetChoice(ANNEALING_CONSTANT);
                     continue;
                  }
                  if(temp=="Boltzmann")
                  {
                     mAnnealingScheduleTemp.SetChoice(ANNEALING_BOLTZMANN);
                     continue;
                  }
                  if(temp=="Cauchy")
                  {
                     mAnnealingScheduleTemp.SetChoice(ANNEALING_CAUCHY);
                     continue;
                  }
                  if(temp=="Exponential")
                  {
                     mAnnealingScheduleTemp.SetChoice(ANNEALING_EXPONENTIAL);
                     continue;
                  }
                  if(temp=="Smart")
                  {
                     mAnnealingScheduleTemp.SetChoice(ANNEALING_SMART);
                     continue;
                  }
                  continue;
               }
               if(tag.GetName()=="AmplitudeScheduleMaxMin")
               {
                  string mut;
                  IOCrystExtractNameSpace(is,mut);
                  is >> mMutationAmplitudeBegin>>mMutationAmplitudeEnd;
                  if(mut=="Constant")
                  {
                     mAnnealingScheduleMutation.SetChoice(ANNEALING_CONSTANT);
                     continue;
                  }
                  if(mut=="Boltzmann")
                  {
                     mAnnealingScheduleMutation.SetChoice(ANNEALING_BOLTZMANN);
                     continue;
                  }
                  if(mut=="Cauchy")
                  {
                     mAnnealingScheduleMutation.SetChoice(ANNEALING_CAUCHY);
                     continue;
                  }
                  if(mut=="Exponential")
                  {
                     mAnnealingScheduleMutation.SetChoice(ANNEALING_EXPONENTIAL);
                     continue;
                  }
                  if(mut=="Smart")
                  {
                     mAnnealingScheduleMutation.SetChoice(ANNEALING_SMART);
                     continue;
                  }
                  continue;
               }
               if(tag.GetName()=="NbTrialRetry")
               {
                  is>>mNbTrialRetry;
               }
               if(tag.GetName()=="MinCostRetry")
               {
                  is>>mMinCostRetry;
               }
               if(tag.GetName()=="MaxNbTrialSinceBest")
               {
                  is>>mMaxNbTrialSinceBest;
               }
               if(tag.GetName()=="RefinedObject")
               {
                  string className,name;
                  IOCrystExtractNameQuoted(is,className);
                  IOCrystExtractNameQuoted(is,name);
                  RefinableObj* obj=& (gRefinableObjRegistry.GetObj(name,className));
                  this->AddRefinableObj(*obj);
               }
               if(tag.GetName()=="CostFunction")
               {
                  string className,name;
                  IOCrystExtractNameQuoted(is,className);
                  IOCrystExtractNameQuoted(is,name);
                  int func;
                  double weight;
                  is >> func>>weight;
                  RefinableObj* obj=& (gRefinableObjRegistry.GetObj(name,className));
                  this->AddCostFunction(*obj,func,weight);
               }
            }//Param
         } while(true);
         break;
      }
      default: cout << "Unknown tag version !"<<endl;
   }
}

const string& GlobalOptimObj::GetName()const { return mName;}
void GlobalOptimObj::SetName(const string& name) {mName=name;}

const string GlobalOptimObj::GetClassName()const { return "GlobalOptimObj";}

void GlobalOptimObj::Print()const {this->Output(cout);}
      
void GlobalOptimObj::PrepareRefParList()
{
   VFN_DEBUG_MESSAGE("GlobalOptimObj::PrepareRefParList()",5)
   
   //:TODO: instead of resetting the list every time, check if it is necessary.
   mRefParList.ResetParList();
   for(int i=0;i<mRecursiveRefinedObjList.GetNb();i++)
      mRefParList.AddPar(mRecursiveRefinedObjList.GetObj(i));
   mRefParList.PrepareForRefinement();
   for(int i=0;i<mRecursiveRefinedObjList.GetNb();i++)
      mRecursiveRefinedObjList.GetObj(i).PrepareForRefinement();
}

void GlobalOptimObj::InitRandomSeedFromTime()const
{
   VFN_DEBUG_MESSAGE("GlobalOptimObj::InitRandomSeedFromTime()",3)
   time_t junk;
   time(&junk);
   tm *tmp=localtime(&junk);
   srand((unsigned)( (*tmp).tm_sec+60* (*tmp).tm_min));
   //for(int i=0;i<20;i++) cout << rand() <<endl;
}

void GlobalOptimObj::NewConfiguration()
{
   for(int i=0;i<mRefinedObjList.GetNb();i++)
      mRefinedObjList.GetObj(i).GlobalOptRandomMove(mMutationAmplitude);
}

void GlobalOptimObj::InitOptions()
{
   VFN_DEBUG_MESSAGE("GlobalOptimObj::InitOptions()",5)
   static string GlobalOptimTypeName;
   static string GlobalOptimTypeChoices[2];//:TODO: Add Genetic Algorithm
   
   static string AnnealingScheduleChoices[5];
   
   static string AnnealingScheduleTempName;
   static string AnnealingScheduleMutationName;
   
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
      
      
      needInitNames=false;//Only once for the class
   }
   mGlobalOptimType.Init(2,&GlobalOptimTypeName,GlobalOptimTypeChoices);
   mAnnealingScheduleTemp.Init(5,&AnnealingScheduleTempName,AnnealingScheduleChoices);
   mAnnealingScheduleMutation.Init(5,&AnnealingScheduleMutationName,AnnealingScheduleChoices);
   VFN_DEBUG_MESSAGE("GlobalOptimObj::InitOptions():End",5)
}

void GlobalOptimObj::UpdateDisplay()
{
   #ifdef __WX__CRYST__
   if(0!=mpWXGlobalOptimObj) mpWXGlobalOptimObj->CrystUpdate();
   #endif
   for(int i=0;i<mRefinedObjList.GetNb();i++) 
      mRefinedObjList.GetObj(i).UpdateDisplay();
}

#ifdef __WX__CRYST__
WXCrystObj* GlobalOptimObj::WXCreate(wxWindow *parent)
{
   mpWXGlobalOptimObj=new WXGlobalOptimObj (parent,this);
   return mpWXGlobalOptimObj;
}
WXCrystObj* GlobalOptimObj::WXGet()
{
   return mpWXGlobalOptimObj;
}
void GlobalOptimObj::WXDelete()
{
   if(0!=mpWXGlobalOptimObj) delete mpWXGlobalOptimObj;
   mpWXGlobalOptimObj=0;
}
void GlobalOptimObj::WXNotifyDelete()
{
   mpWXGlobalOptimObj=0;
}
#endif

}//namespace
