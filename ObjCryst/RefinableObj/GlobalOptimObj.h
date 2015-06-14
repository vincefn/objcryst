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
/*   OptimizationObj.h
*  header file for the Global Optimization object
*
*/
#ifndef _GLOBALOPTIMOBJ_H
#define _GLOBALOPTIMOBJ_H

#include <map>
#include "ObjCryst/ObjCryst/General.h"

namespace ObjCryst
{
   class OptimizationObj;
   class MonteCarloObj;
}

#include "ObjCryst/RefinableObj/RefinableObj.h"
#include "ObjCryst/RefinableObj/LSQNumObj.h"
#include "ObjCryst/RefinableObj/IO.h"
#include "ObjCryst/RefinableObj/Tracker.h"
#include <string>
#include <iostream>
#ifdef __WX__CRYST__
   //#undef GetClassName // Conflict from wxMSW headers ? (cygwin)
#include "ObjCryst/wxCryst/wxGlobalOptimObj.h"
#endif

namespace ObjCryst
{
/** Annealing schedule type. Used to determine the variation of the
* temperature and the mutation amplitude
* 
* With A=Temperature or A=MutationAMplitude, and the corresponding
, min and max values supplied (the latter is ignored for constant,
* Cauchy and Boltzmann), with 'step' being the current step, and
* NbStep the total number of steps. (In the Parallel Tempering
* algorithm, a 'step' denotes one of the parallel refinement).
* \f[ A_{constant} = A_{min}\f]
* \f[ A_{Boltzmann} = A_{min} \frac{\ln(NbStep)}{\ln(step)} \f]
* \f[ A_{Cauchy} = A_{min} \frac{NbStep}{step} \f]
* \f[ A_{exponential} = A_{max} (\frac{T_{min}}{T_{max}})^{\frac{step}{NbStep}} \f]
* \f[ A_{\gamma} = A_{max} +(A_{min}-A_{max}*(\frac{step}{NbStep})^{\gamma} \f]
*
* For the 'smart' schedule, it is only supported so far for the mutation amplitude:
* it is modulated so that for each temperature between 30 and 70% of trial
* configurations are accepted, within the limits for the mutation.
*/
enum AnnealingSchedule
{
   ANNEALING_CONSTANT,
   ANNEALING_BOLTZMANN,
   ANNEALING_CAUCHY,
//   ANNEALING_QUENCHING,
   ANNEALING_EXPONENTIAL,
   ANNEALING_SMART,
   ANNEALING_GAMMA
};

/** Global optimization type. Eventually it would be better to build
* a base Global Optimization  (or even Optimization) object, and to derive
* it in different classes for Simulated Annealing, Parallel Tempering,
* Genetic Algorithm,...
*/
enum GlobalOptimType
{
   GLOBAL_OPTIM_SIMULATED_ANNEALING,
   GLOBAL_OPTIM_PARALLEL_TEMPERING,
   GLOBAL_OPTIM_RANDOM_LSQ,
   GLOBAL_OPTIM_SIMULATED_ANNEALING_MULTI,
   GLOBAL_OPTIM_PARALLEL_TEMPERING_MULTI,   
};

/** \brief Base object for Optimization methods.
*
* This is an abstract base class, derived for Monte-Cralo type algorithms
* (Simulated Annealing & Parallel Tempering), and hopefully soon
* for Genetic Algorithms.
*
* \remarks Instead of keeping a copy of the list of parameters here,
* maybe it would be better to delegate all parameter handling to the refined 
* objects (they would also have to keep in memory the saved parameter sets, so
* that could be difficult to administrate...).
*/
class OptimizationObj
{
   public:
      /// Constructor
      OptimizationObj(const string name="");
      /// Destructor
      virtual ~OptimizationObj();
      
      /** \brief Randomize starting configuration. Only affects limited and periodic parameters.
      */
      virtual void RandomizeStartingConfig();
      /// Launch optimization (a single run) for N steps
      /// \param nbSteps: the number of steps to go. This number is modified (decreases!)
      /// as the refinement goes on.
      /// \param silent : if true, absolutely no message should be printed (except debugging)
      /// \param finalcost: the optimization will stop if overall cost fallse below this value
      /// \param maxTime: the optimization will stop after the given number of seconds has
      /// been spent optimizing (ignored if <0).
      virtual void Optimize(long &nbSteps,const bool silent=false,const REAL finalcost=0,
                            const REAL maxTime=-1)=0;
      /** Launch optimization for multiple runs of N steps
      * \param nbCycle: the number of runs (cycles) to perform. The structure is randomized
      * at the beginning of each cycle. If nbCycle==-1, this will run indefinitely.
      * The nbCycle parameter is decreased after each run.
      * \param nbSteps: the number of steps to go. This number is modified (decreases!)
      * as the refinement goes on.
      * \param silent : if true, absolutely no message should be printed (except debugging)
      * \param finalcost: the optimization will stop if overall cost fallse below this value
      * \param maxTime: the optimization will stop after the given number of seconds has
      * been spent optimizing (ignored if <0).
      */
      virtual void MultiRunOptimize(long &nbCycle,long &nbSteps,const bool silent=false,const REAL finalcost=0,
                                    const REAL maxTime=-1)=0;
   //Set Refinable parameters status
      /// Fix all parameters
      void FixAllPar();
      /// Fix one parameter
      void SetParIsFixed(const string& parName,const bool fix);
      /// Fix one family of parameters
      void SetParIsFixed(const RefParType *type,const bool fix);
      /// UnFix All parameters
      void UnFixAllPar();
      /// Set a parameter to be used
      void SetParIsUsed(const string& parName,const bool use);
      /// Set a family of parameters to be used
      void SetParIsUsed(const RefParType *type,const bool use);
      /// Change the relative limits for a parameter from its name
      void SetLimitsRelative(const string &parName, const REAL min, const REAL max);
      /// Change the relative limits for a family of parameter
      void SetLimitsRelative(const RefParType *type, const REAL min, const REAL max);
      /// Change the absolute limits for a parameter from its name
      void SetLimitsAbsolute(const string &parName, const REAL min, const REAL max);
      /// Change the absolute limits for a family of parameter
      void SetLimitsAbsolute(const RefParType *type, const REAL min, const REAL max);
      
      /** \brief The optimized (minimized, actually) function.
      *
      * This function is the weighted sum of the chosen Cost Functions for
      * the refined objects.
      */
      virtual REAL GetLogLikelihood()const;

      /// Stop after the current cycle. USed for interactive refinement.
      void StopAfterCycle();
      /// Show report to the user during refinement. Used for GUI update.
      virtual void DisplayReport();

      /// Add a refined object. All sub-objects are also added
      void AddRefinableObj(RefinableObj &);
      /// Get the RefinableObj with all the parameters from all refined objects.
      /// If rebuild=true, prepare again the list of objects/parameters.
      RefinableObj& GetFullRefinableObj(const bool rebuild=true);
      /** \brief Output a description of the object in XML format to a stream.
      *
      * This saves the list of refined object and the cost functions, as well as options
      * for the refinement. The refined objects are \b not saved, so this must be done
      * somewhere else (they must be reloaded before this object).
      */
      virtual void XMLOutput(ostream &os,int indent=0)const=0;
      /** \brief Input in XML format from a stream, restoring the set of refined
      * objects and the associated cost functions. Note that the corresponding objects
      * must have been loaded in memory before, else shit happens.
      *
      */
      virtual void XMLInput(istream &is,const XMLCrystTag &tag)=0;
      //virtual void XMLInputOld(istream &is,const IOCrystTag &tag);
      /// Get the name for this object
      const string& GetName()const;
      /// Set the name for this object
      void SetName(const string&);
      /// Get the name for this class type
      virtual const string GetClassName()const;
      /// Print some information about this object
      virtual void Print()const;
      /// Restore the Best configuration
      void RestoreBestConfiguration();
      /// Are we busy optimizing ?
      bool IsOptimizing()const;
      /** During a global optimization, tell all objects that the current config is
      * the latest "best" config.
      */
      void TagNewBestConfig();
      /// Get the elapsed time (in seconds) during the last optimization
      REAL GetLastOptimElapsedTime()const;
      /// Get the MainTracker
      MainTracker& GetMainTracker();
      /// Get the MainTracker
      const MainTracker& GetMainTracker()const;
      RefObjOpt& GetXMLAutoSaveOption();
      const RefObjOpt& GetXMLAutoSaveOption()const;
      /// Access to current best cost
      const REAL& GetBestCost()const;
      /// Access to current best cost
      REAL& GetBestCost();
      /// Begin optimization for all objects
      virtual void BeginOptimization(const bool allowApproximations=false,
                                     const bool enableRestraints=false);
      /// End optimization for all objects
      virtual void EndOptimization();
      /// Number of trial per run
      virtual long& NbTrialPerRun();
      /// Number of trial per run
      virtual const long& NbTrialPerRun() const;
   protected:
      /// \internal Prepare mRefParList for the refinement
      void PrepareRefParList();
      
      /// Initialization of options.
      virtual void InitOptions();
      /// Update Display (if any display is available), when a new 'relevant' configuration
      /// is reached. This calls all RefinableObj::UpdateDisplay()
      virtual void UpdateDisplay();
      /// (Re)build OptimizationObj::mRecursiveRefinedObjList, if an
      /// object has been added or modified. If no object has been
      /// added and no sub-object has been added/removed, then nothing is done.
      void BuildRecursiveRefObjList();
      /// The refinable par list used during refinement. Only a condensed version
      /// of all objects. This is useful to keep an history of modifications, and to
      /// restore previous values.
      /// \remarks maybe this should be completely delegated to the refined objetcs.
      mutable RefinableObj mRefParList;
      /// Name of the GlobalOptimization object
      string mName;
      /// File name where refinement info is saved (NOT USED so far...)
      string mSaveFileName;
      
      //Status of optimization
         /// Number of trial per run, to be saved/restored in XML output
         long mNbTrialPerRun;
         /// Number of trials so far
         long mNbTrial;
         /// Best value of the cost function so far
         REAL mBestCost;
         /// Index of the 'best' saved parameter set
         long mBestParSavedSetIndex;
         /// The current 'context', in the case the optimization is run in different
         /// parallel contexts
         unsigned long mContext;
         /// Statistics about each object contributing to the overall Log(likelihood)
         struct LogLikelihoodStats
         {
            /// Previous log(likelihood)
            REAL mLastLogLikelihood;
            /// Total Log(Likelihood), to compute the average
            REAL mTotalLogLikelihood;
            /// total of (Delta(Log(Likelihood)))^2 between successive trials
            REAL mTotalLogLikelihoodDeltaSq;
         };
         /// Statistics for each context
         /// (mutable for dynamic update during optimization)
         mutable map <unsigned long, map<const RefinableObj*,LogLikelihoodStats> > mvContextObjStats;
      // Dynamic weights (EXPERIMENTAL!)
         struct DynamicObjWeight
         {
            DynamicObjWeight():mWeight(1.){};
            REAL mWeight;
         };
         /// Weights for each objects in each context
         /// (mutable for dynamic update during optimization)
         mutable map<const RefinableObj*,DynamicObjWeight> mvObjWeight;

         /// List of saved parameter sets. This is used to save possible
         /// solutions during the optimization, so that the user can check them
         /// afterwards.
         ///
         /// The first member of each pair is the \e index of the parameter set,
         /// and the second is the overall cost for that set.
         std::vector<pair<long,REAL> > mvSavedParamSet;
         
      /// True if a refinement is being done. For multi-threaded environment
      bool mIsOptimizing;
      /// If true, then stop at the end of the cycle. Used in multi-threaded environment
      bool mStopAfterCycle;
      
      // Refined objects
         /// The refined objects
         ObjRegistry<RefinableObj> mRefinedObjList;
         /// The refined objects, recursively including all sub-objects.
         /// This is mutable, since it is a function of mRefinedObjList only.
         mutable ObjRegistry<RefinableObj> mRecursiveRefinedObjList;
         
      /// Periodic save of complete environment as an xml file
         RefObjOpt mXMLAutoSave;
      
      /// The time elapsed after the last optimization, in seconds
         REAL mLastOptimTime;
      /// MainTracker object to track the evolution of cost functions, likelihood,
      /// and individual parameters.
      MainTracker mMainTracker;
   private:
   #ifdef __WX__CRYST__
   public:
      /// Create a WXCrystObj for this object.
      virtual WXCrystObjBasic* WXCreate(wxWindow*)=0;
      virtual WXOptimizationObj* WXGet()=0;
      virtual void WXDelete()=0;
      virtual void WXNotifyDelete()=0;
   protected:
      //:TODO: remove this !
      friend class ObjCryst::WXOptimizationObj;
      //friend class ObjCryst::WXMonteCarloObj;
      //friend class ObjCryst::WXGeneticAlgorithm;
      /// Mutex used to protect mStopAfterCycle
      wxMutex mMutexStopAfterCycle;
   #endif
};

/** \brief Base object for Monte-Carlo Global Optimization methods.
*
* The algorithm is quite simple, whith two type of optimizations, either
* simulated Annealing or Parallel Tempering, the latter being recommanded
* for most real-world optimizations.
*
*/
class MonteCarloObj:public OptimizationObj
{
   public:
      /// Constructor
      MonteCarloObj(const string name="");
      /// Constructor. Using internalUseOnly=true will avoid registering the
      /// the object to any registry, and thus (for example) no display will be created,
      /// nor will this object be automatically be saved.
      MonteCarloObj(const bool internalUseOnly);
      /// Destructor
      virtual ~MonteCarloObj();
      
      /** \brief  Set the refinement method to simulated Annealing. Note that
      * Parellel Tempering is more efficient to get out of local minima, so you sould
      * rather use that method.
      *
      * The refinement begins at max and finishes at min temperature.
      *
      *\param scheduleTemp: temperature schedule. See AnnealingSchedule.
      *\param tMax,tMin: Max and Min temperatures.
      *\param scheduleMutation: the mutation schedule. For each new configuration, the
      * variation of each refinable parameter is less than its RefinablePar::GlobalOptimStep(),
      * multiplied by the current mutation amplitude. By default this mutation is equal to 1.,
      * but making bigger steps is a good idea at the beginning of the refinement
      * (for higher temperatures). See AnnealingSchedule. See AnnealingSchedule.
      *\param mutMax,mutMin: Max and Min mutation amplitudes.
      *\param minCostRetry, nbTrialRetry: if after nbTrialRetry, the cost function is still
      * above minCostRetry, then start again from a random configuration. No randomization is
      * made if nbTrialRetry <= 0.
      *\param maxNbTrialSinceBest: if more than maxNbTrialSinceBest trials have been made
      * since the best configuration was recorded, then revert to that configuration. This
      * should be large enough to have an ergodic search (the default is never to revert..)
      *
      * \warning do not use the 'smart' option for the temperature schedule, it is not yet 
      * implemented. Later it will be used to set the temperatures as a function of 
      * the amplitude schedule, so that we accept between 30% and 70%  moves.
      * \note this will be removed when we separate the different algorithms in different
      * classes.
      */
      void SetAlgorithmSimulAnnealing(const AnnealingSchedule scheduleTemp,
                                 const REAL tMax, const REAL tMin,
                                 const AnnealingSchedule scheduleMutation=ANNEALING_CONSTANT,
                                 const REAL mutMax=16., const REAL mutMin=.125,
                                 const long nbTrialRetry=0,const REAL minCostRetry=0.);
      /** \brief  Set the refinement method to Parallel Tempering.
      *
      * The refinement begins at max and finishes at min temperature.
      *
      *\param scheduleTemp: temperature schedule
      *\param tMax,tMin: Max and Min temperatures. See AnnealingSchedule.
      *\param scheduleMutation: the mutation schedule. For each new configuration, the
      * variation of each refinable parameter is less than its RefinablePar::GlobalOptimStep(),
      * multiplied by the current mutation amplitude. By default this mutation is equal to 1.,
      * but making bigger steps can be a good idea at the beginning of the refinement. Thus
      * you can choose a schedule for the amplitude, exactly like for the temperature. 
      *\param mutMax,mutMin: Max and Min mutation amplitudes. 
      * \warning do not use the 'smart' option for the temperature schedule, it is not yet 
      * implemented. Later it will be used to set the temperatures as a function of 
      * the amplitude schedule, so that we keep accepted move between 30% and 70%.
      * \note this will be removed when we separate the different algorithms in different
      * classes.
      */
      void SetAlgorithmParallTempering(const AnnealingSchedule scheduleTemp,
                                 const REAL tMax, const REAL tMin,
                                 const AnnealingSchedule scheduleMutation=ANNEALING_CONSTANT,
                                 const REAL mutMax=16., const REAL mutMin=.125);
      
      virtual void Optimize(long &nbSteps,const bool silent=false,const REAL finalcost=0,
                            const REAL maxTime=-1);
      virtual void MultiRunOptimize(long &nbCycle,long &nbSteps,const bool silent=false,const REAL finalcost=0,
                                    const REAL maxTime=-1);
      
      /** \internal Do a single simulated annealing run. This is called by Optimize(...) and
      * MultiRunOptimize(), which must also prepare the optimization (PrepareRefParList(), etc..).
      */
      void RunSimulatedAnnealing(long &nbSteps,const bool silent=false,const REAL finalcost=0,
                                 const REAL maxTime=-1);
      /** \internal Do a single Parallel Tempering run. This is called by Optimize(...) and
      * MultiRunOptimize(), which must also prepare the optimization (PrepareRefParList(), etc..).
      */
      void RunParallelTempering(long &nbSteps,const bool silent=false,const REAL finalcost=0,
                                const REAL maxTime=-1);

      void RunRandomLSQMethod(long &nbCycle);

      void RunNondestructiveLSQRefinement(  int nbCycle=1,bool useLevenbergMarquardt=false, 
                                            const bool silent=false, const bool callBeginEndOptimization=true, 
                                            const float minChi2var=0.01 );

      
      //Parameter Access by name
      //RefinablePar& GetPar(const string& parName);
      
      // Print information about the current state of optimization (parameters value, 
      // characteristic figures...)
      //virtual ostream& operator<<(ostream& os)const;
      virtual void XMLOutput(ostream &os,int indent=0)const;
      virtual void XMLInput(istream &is,const XMLCrystTag &tag);
      //virtual void XMLInputOld(istream &is,const IOCrystTag &tag);
      virtual const string GetClassName()const;
      /// Access to the builtin LSQ optimization object
      LSQNumObj & GetLSQObj();
      /// Access to the builtin LSQ optimization object
      const LSQNumObj & GetLSQObj() const;
      
      /** Prepare mLSQ for least-squares refinement during the global optimization
      *
      * \param useFullPowderPatternProfile: if true, the refinement will use the full
      * profile version of powder patterns, otherwise only the integrated powder pattern
      * will be used (faster).
      */
      virtual void InitLSQ(const bool useFullPowderPatternProfile=true);
   protected:
      
      /** \brief Make a random change in the configuration.
      * 
      * \internal
      *  This just generates a new configuration with random changes (according
      * to current parameters). The new configuration is \e not tested \e in \e this \e function
      * vs temperature: this should be done in the OptimizationObj::Optimize() function,
      * which also chooses whether to revert to the previous configuration.
      *
      * Random moves are made by the objects and not by this function,
      * because the new configuration can be specific
      * (like, for example, permutations between some of the parameters (atoms)).
      * 
      * \param type: can be used to restrict the move to a given category of parameters.
      */
      virtual void NewConfiguration(const RefParType *type=gpRefParTypeObjCryst);
      
      virtual void InitOptions();
            
      /// Method used for the global optimization. Should be removed when we switch
      /// to using several classes for different algorithms.
      RefObjOpt mGlobalOptimType;
      
      //Status of optimization
         /// Current value of the cost function
         REAL mCurrentCost;
         
      // History, for experimental purposes only !
         /// Option to save the evolution of tracked data (cost functions,
         /// likelihhod, individual parameters,...)
         RefObjOpt mSaveTrackedData;
         
      // Annealing parameters
         /// Current temperature for annealing
         REAL mTemperature;
         /// Beginning temperature for annealing
         REAL mTemperatureMax;
         /// Lower temperature
         REAL mTemperatureMin;
         /// Schedule for the annealing
         RefObjOpt mAnnealingScheduleTemp;
         /// Gamma for the 'gamma' temperature schedule
         REAL mTemperatureGamma;
      //Parameters to create new configurations
         /// Mutation amplitude. From .25 to 64. Random moves will have a maximum amplitude
         /// equal to this amplitude multiplied by the Global optimization step defined
         /// for each RefinablePar. Large amplitude should be used at the beginning of the
         /// refinement (high temeratures).
         REAL mMutationAmplitude;
         /// Mutation amplitude at the beginning of the optimization.
         REAL mMutationAmplitudeMax;
         /// Mutation amplitude at the end of the optimization.
         REAL mMutationAmplitudeMin;
         /// Schedule for the annealing
         RefObjOpt mAnnealingScheduleMutation;
         /// Gamma for the 'gamma' Mutation amplitude schedule
         REAL mMutationAmplitudeGamma;
      //Automatic retry 
         /// Number of trials before testing if we are below the given minimum cost.
         /// If <=0, this will be ignored.
         long mNbTrialRetry;
         /// Cost to reach unless an automatic randomization and retry is done
         REAL mMinCostRetry;
      /// Least squares object
      LSQNumObj mLSQ;
      /// Option to run automatic least-squares refinements
      RefObjOpt mAutoLSQ;
   private:
   #ifdef __WX__CRYST__
   public:
      virtual WXCrystObjBasic* WXCreate(wxWindow*);
      virtual WXOptimizationObj* WXGet();
      virtual void WXDelete();
      virtual void WXNotifyDelete();
   protected:
      WXMonteCarloObj* mpWXCrystObj;
      friend class ObjCryst::WXMonteCarloObj;
   #endif
};

/// Global Registry for all OptimizationObj
extern ObjRegistry<OptimizationObj> gOptimizationObjRegistry;

}// namespace

#endif //_GLOBALOPTIMOBJ_H
