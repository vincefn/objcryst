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

#include "ObjCryst/General.h"

namespace ObjCryst
{
   class OptimizationObj;
   class MonteCarloObj;
}

#include "RefinableObj/RefinableObj.h"
#include "RefinableObj/IO.h"
#include <string>
#include <iostream>
#ifdef __WX__CRYST__
   //#undef GetClassName // Conflict from wxMSW headers ? (cygwin)
#include "wxCryst/wxGlobalOptimObj.h"
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
   ANNEALING_SMART
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
   GLOBAL_OPTIM_GENETIC
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
      /// Launch optimization for N steps
      /// \param nbSteps: the number of steps to go. This number is modified (decreases!)
      /// as the refinement goes on.
      /// \param silent : if true, absolutely no message should be printed (except debugging)
      virtual void Optimize(long &nbSteps,const bool silent=false,const REAL finalcost=0)=0;
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
      * the refined objects. All Cost Functions \b must be strictly positive.
      */
      REAL GetCostFunctionValue();

      /// Stop after the current cycle. USed for interactive refinement.
      void StopAfterCycle();
      /// Show report to the user during refinement. Used for GUI update.
      virtual void DisplayReport();

      /// Add a refined object. All sub-objects are also added
      void AddRefinableObj(RefinableObj &);
      /// Add a cost function, with a given weight. This cost function
      /// should be strictly positive, and ideally should behave like a R/Rw function,
      /// ie a value above 0.50 corresponds to a very inadequate configuration,
      /// while 0.05 is excellent.
      void AddCostFunction(RefinableObj &,const unsigned int id, const REAL weight=1.);
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
   protected:
      /// \internal Prepare mRefParList for the refinement
      void PrepareRefParList();
      
      /// \internal Initialize random seed from time
      void InitRandomSeedFromTime()const;
      /// Initialization of options.
      virtual void InitOptions()=0;
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
         /// Number of trials so far
         long mNbTrial;
         /// Best value of the cost function so far
         REAL mBestCost;
         /// Index of the 'best' saved parameter set
         long mBestParSavedSetIndex;
         
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
         /// Number of Cost Functions used
         unsigned int mNbCostFunction;
         /// Max number of Cost Functions (dynamically adjusted)
         unsigned int mMaxNbCostFunction;
         /// The objects with the Cost functions
         RefinableObj *mpCostFunctionRefinableObj[20];
         /// The id of the cost functions in each RefinableObj
         CrystVector_int mpCostFunctionId;
         /// The weight associated with each cost function
         CrystVector_REAL mCostFunctionWeight;
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
                                 const long nbTrialRetry=0,const REAL minCostRetry=0.,
                                 const long maxNbTrialSinceBest=0);
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
      * the amplitude schedule, so that we keep accepeted move between 30% and 70%.
      * \note this will be removed when we separate the different algorithms in different
      * classes.
      */
      void SetAlgorithmParallTempering(const AnnealingSchedule scheduleTemp,
                                 const REAL tMax, const REAL tMin,
                                 const AnnealingSchedule scheduleMutation=ANNEALING_CONSTANT,
                                 const REAL mutMax=16., const REAL mutMin=.125);
      
      virtual void Optimize(long &nbSteps,const bool silent=false,const REAL finalcost=0);
      
      //Parameter Access by name
      //RefinablePar& GetPar(const string& parName);
      
      // Print information about the current state of optimization (parameters value, 
      // characteristic figures...)
      //virtual ostream& operator<<(ostream& os)const;
      /// Save history of the evolution of parameters to a file. Only non-fixed parameters
      /// are saved. This saves a very crude array in which can bve found the value of
      /// all non-fixed parameters for successive "best" configurations.
      void SaveOptimHistory() const;
      virtual void XMLOutput(ostream &os,int indent=0)const;
      virtual void XMLInput(istream &is,const XMLCrystTag &tag);
      //virtual void XMLInputOld(istream &is,const IOCrystTag &tag);
      virtual const string GetClassName()const;
   protected:
      
      /** \brief Make a random change in the configuration.
      * 
      * \internal
      *  This just generates a new configuration with random changes (according
      * to current parameters). The old config is stored in mRefParList as the
      * last config (index mLastParSavedSetIndex). The new one is \e not tested \e in \e this \e function
      * vs temperature: this should be done in the OptimizationObj::Optimize() function,
      * which also chooses whether to revert to the previous configuration.
      *
      * Random moves are made by the objects and not by this function,
      * because the new configuration can be specific
      * (like, for example, permutations between some of the parameters (atoms)).
      */
      virtual void NewConfiguration();
      
      virtual void InitOptions();
            
      /// Method used for the global optimization. Should be removed when we switch
      /// to using several classes for different algorithms.
      RefObjOpt mGlobalOptimType;
      
      //Status of optimization
         /// Current value of the cost function
         REAL mCurrentCost;
         
      //Keep an history with the evolution of optimization
         /// Total number of saved configurations
         long mHistoryNb;
         /// Trials corresponding to each stored values
         CrystVector_long mHistoryTrialNumber;
         /// Evolution of cost function
         CrystVector_REAL mHistoryCostFunction;
         /// Index of saved parameters set in mRefParList for each saved trial
         CrystVector_long mHistorySavedParamSetIndex;
         /// Save the evolution of refined parameters after optimization ?
         bool mHistorySaveAfterEachOptim;
         /// Save the evolution of refined parameters after optimization ?
         string mHistorySaveFileName;
         
         /// Index of the 'last' parameter set
         long mLastParSavedSetIndex;
      
      // Annealing parameters
         /// Beginning temperature for annealing
         REAL mTemperatureMax;
         /// Lower temperature
         REAL mTemperatureMin;
         /// Schedule for the annealing
         RefObjOpt mAnnealingScheduleTemp;
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
      //Automatic retry 
         /// Number of trials before testing if we are below the given minimum cost.
         /// If <=0, this will be ignored.
         long mNbTrialRetry;
         /// Cost to reach unless an automatic randomization and retry is done
         REAL mMinCostRetry;
         /// If more than mMaxNbTrialSinceBest trials have been made since the best
         /// configuration has been found, then revert to the best configuration. If <=0,
         /// then this is ignored. This must be large enough to have an ergodic 
         /// algorithm (more strictly, should not be used ?)
         long mMaxNbTrialSinceBest;
      
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
