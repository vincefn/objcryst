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
/*   GlobalOptimObj.h
*  header file for Global Optimization Objects
*
*/
#ifndef _GLOBALOPTIMOBJ_H
#define _GLOBALOPTIMOBJ_H

#include "ObjCryst/General.h"

#ifdef __WX__CRYST__
namespace ObjCryst
{
	class GlobalOptimObj;
}
   //#undef GetClassName // Conflict from wxMSW headers ? (cygwin)
#include "wxCryst/wxGlobalOptimObj.h"
#endif

#include "RefinableObj/RefinableObj.h"
#include "RefinableObj/IO.h"
#include <string>
#include <iostream>

namespace ObjCryst
{

enum AnnealingSchedule
{
   ANNEALING_CONSTANT,
   ANNEALING_BOLTZMANN,
   ANNEALING_CAUCHY,
//   ANNEALING_QUENCHING,
   ANNEALING_EXPONENTIAL,
   ANNEALING_SMART
};
enum GlobalOptimType
{
   GLOBAL_OPTIM_SIMULATED_ANNEALING,
   GLOBAL_OPTIM_PARALLEL_TEMPERING,
   GLOBAL_OPTIM_GENETIC
};

/** \brief Base object for Global Optimization method. This is still \e crude.
*
*/

class GlobalOptimObj
{
	public:
		GlobalOptimObj(const string name="");
		virtual ~GlobalOptimObj();
      
      /** \brief Randomize starting configuration
      *
      *  Only limited parameters are modified by this function. Others are not
      * affected.
      */
      void RandomizeStartingConfig();
      
      /** \brief  Set the refinement method to simulated Annealing, and
      *
      * The refinement begins at max and finishes at min temperature.
      *
      *\param scheduleTemp: temperature schedule
      *\param tMax,tMin: Max and Min temperatures. The Max temperature will be ignored
      * for constant, Cauchy, and Boltzmann temperature schedules.
      *\param scheduleMutation: the mutation schedule. For each new configuration, the
      * variation of each refinable parameter is less than its RefinablePar::GlobalOptimStep(),
      * multiplied by the current mutation amplitude. By default this mutation is equal to 1.,
      * but making bigger steps can be a good idea at the beginning of the refinement. Thus
      * you can choose a schedule for the amplitude, exactly like for the temperature.
      *\param mutMax,mutMin: Max and Min mutation amplitudes. The Max temperature will 
      * be ignored for constant, Cauchy, and Boltzmann schedules.
      *\param minCostRetry, nbTrialRetry: if after nbTrialRetry, the cost function is still
      * above minCostRetry, then start again from a random configuration. No randomization is
      * made if nbTrialRetry <= 0.
      *\param maxNbTrialSinceBest: if more than maxNbTrialSinceBest trials have been made
      * since the best configuration was recorded, then revert to that configuration. This
      * should be large enough to have an ergodic search (the default is never to revert..)
      */
      void SetAlgorithmSimulAnnealing(const AnnealingSchedule scheduleTemp,
                                 const double tMax, const double tMin,
                                 const AnnealingSchedule scheduleMutation=ANNEALING_CONSTANT,
                                 const double mutMax=1., const double mutMin=1.,
                                 const long nbTrialRetry=0,const double minCostRetry=0.,
                                 const long maxNbTrialSinceBest=0);
      /** \brief  Set the refinement method to Parallel Tempering.
      *
      * The refinement begins at max and finishes at min temperature.
      *
      *\param scheduleTemp: temperature schedule
      *\param tMax,tMin: Max and Min temperatures. The Max temperature will be ignored
      * for constant, Cauchy, and Boltzmann temperature schedules.
      *\param scheduleMutation: the mutation schedule. For each new configuration, the
      * variation of each refinable parameter is less than its RefinablePar::GlobalOptimStep(),
      * multiplied by the current mutation amplitude. By default this mutation is equal to 1.,
      * but making bigger steps can be a good idea at the beginning of the refinement. Thus
      * you can choose a schedule for the amplitude, exactly like for the temperature. 
      *\param mutMax,mutMin: Max and Min mutation amplitudes. The Max will 
      * be ignored for constant, Cauchy, and Boltzmann schedules. Both parameters
      * are ignored for 'smart' schedule.
      */
      void SetAlgorithmParallTempering(const AnnealingSchedule scheduleTemp,
                                 const double tMax, const double tMin,
                                 const AnnealingSchedule scheduleMutation=ANNEALING_CONSTANT,
                                 const double mutMax=1., const double mutMin=1.);
      
      /// Launch optimization for N steps
      /// \param nbSteps: the number of steps to go. This number is modified (decreases!)
      /// as the refinement goes on.
      void Optimize(long &nbSteps);
      
      //Parameter Access by name
      //RefinablePar& GetPar(const string& parName);
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
      void SetLimitsRelative(const string &parName, const double min, const double max);
      void SetLimitsRelative(const RefParType *type, const double min, const double max);
      void SetLimitsAbsolute(const string &parName, const double min, const double max);
      void SetLimitsAbsolute(const RefParType *type, const double min, const double max);

      
      // Print information about the current state of optimization (parameters value, 
      // characteristic figures...)
      //virtual ostream& operator<<(ostream& os)const;
      /// Save history of the evolution of parameters to a file. Only non-fixed parameters
      /// are saved.
      void SaveOptimHistory() const;
      /** \brief The optimized (minimized, actually) function.
      *
      * It \b must be strictly positive.
      */
      virtual double GetCostFunctionValue();
      /// Stop after the current cycle
      void StopAfterCycle();
      /// Show report to the user during refinement. Overloaded for GUI update 8-)
      virtual void DisplayReport();
      /// Add a refined object
      void AddRefinableObj(RefinableObj &);
      /// Add a cost function
      void AddCostFunction(RefinableObj &,const unsigned int id, const double weight=1.);
      /** \brief Output to stream
      *
      */
      virtual void Output(ostream &os,int indent=0)const;
      /** \brief Input From stream
      *
      */
      virtual void Input(istream &is,const XMLCrystTag &tag);
      //virtual void InputOld(istream &is,const IOCrystTag &tag);
      /// Get the name for this object
      const string& GetName()const;
      /// Set the name for this object
      void SetName(const string&);
      /// Get the name for this class type
      const string GetClassName()const;
      /// Print the configuration for this object
      void Print()const;
	protected:
      
      /** \brief Make a random change in the configuration.
      *
      *  This just generates a new configuration with random changes (according
      * to current parameters). The old config is stored in mRefParList as the
      * last config (index mLastParSavedSetIndex). The new one is \e not tested
      * vs temperature: this should be done in the GlobalOptimObj::Optimize() function,
      * which also chooses whether to return to the previous configuration.
      *
      * The reason this configuration generation is not incorporated to
      * GlobalOptimObj::Optimize() is that the latter is a \e general algorythm, valid
      * for any kind of optimized object, while the new configuration can be specific
      * (like, for example, permutations between some of the parameters (atoms)).
      */
      virtual void NewConfiguration();
      
      /// \internal Prepare mRefParList for the refinement
      void PrepareRefParList();
      
      /// \internal Initialize random seed from time
      void InitRandomSeedFromTime()const;
      
      /// Initial initialization of options.
      void InitOptions();
      /// Update Display (if any diplay is available), when a new 'relevant' configuration
      /// is reached. This calls all RefinableObj::UpdateDisplay()
      void UpdateDisplay();
      
      /// The refinable par list used during refinement. Only a condensed version
      /// of all objects. This is useful to keep an history of modifications, and to
      /// restore previous values
		mutable RefinableObj mRefParList;
      ///Name of the refined object
		string mName;
      ///File name where refinement info is saved (NOT USED so far...)
		string mSaveFileName;
      
      
      ///Method used for the global optimization
      RefObjOpt mGlobalOptimType;
      
      //Status of optimization
         /// Number of trials so far
         long mNbTrial;
         /// Current value of the cost function
         double mCurrentCost;
         /// Best value of the cost function so far
         double mBestCost;
         /* :TODO:
         // Current cost array. one row contains the different contributions to the cost,
         // if more than one cost function is used, and each row correspond to each
         // 'world' refinement (one world for simulated annealing, more for parallel
         // tempering).
         // This is used for reporting results.
         CrystMatrix_double mCurrentCostArray;
         */
         
      //Keep an history with the evolution of optimization
         /// Total number of saved configurations
         long mHistoryNb;
         /// Trials corresponding to each stored values
         CrystVector_long mHistoryTrialNumber;
         /// Evolution of cost function
         CrystVector_double mHistoryCostFunction;
         /// Index of saved parameters set in mRefParList for each saved trial
         CrystVector_long mHistorySavedParamSetIndex;
         /// Save the evolution of refined parameters after optimization ?
         bool mHistorySaveAfterEachOptim;
         /// Save the evolution of refined parameters after optimization ?
         string mHistorySaveFileName;
         
         /// Index of the 'best' saved parameter set
         long mBestParSavedSetIndex;
         /// Index of the 'last' parameter set
         long mLastParSavedSetIndex;
      
      //Simulated Annealing parameters
         /// Beginning temperature for annealing
         double mTemperatureMax;
         /// Lower temperature
         double mTemperatureMin;
         /// Schedule for the annealing
         RefObjOpt mAnnealingScheduleTemp;
      //Parameters to create new configurations
         /// Mutation amplitude. From 1 to 100. Random moves will have a maximum amplitude
         /// equal to this amplitude multiplied by the Global optimization step defined
         /// for each RefinablePar. Large amplitude should be used at the beginning of the
         /// refinement.
         double mMutationAmplitude;
         /// Mutation amplitude at the beginning of the optimization.
         double mMutationAmplitudeBegin;
         /// Mutation amplitude at the end of the optimization.
         double mMutationAmplitudeEnd;
         /// Schedule for the annealing
         RefObjOpt mAnnealingScheduleMutation;
      //Automatic retry 
         /// Number of trials before testing if we are below the given minimum cost.
         /// If <=0, this will be ignored
         long mNbTrialRetry;
         /// Cost to reach unless an automatic randomization and retry is done
         double mMinCostRetry;
         /// If more than mMaxNbTrialSinceBest trials have been made since the best
         /// configuration has been found, then revert to the best configuration. If <=0,
         /// then this is ignored. This must be large enough to have an ergodic 
         /// algorithm (more strictly, should not be used ?)
         long mMaxNbTrialSinceBest;
      /// True if a refinement is being done. For multi-thread environment
      bool mIsOptimizing;
      /// If true, then stop at the end of the cycle. Used in multi-thread environment
      bool mStopAfterCycle;
      
      // Refined objects
         /// The refined objects
         ObjRegistry<RefinableObj> mRefinedObjList;
         /// The refined objects, recursively including all sub-objects
         ObjRegistry<RefinableObj> mRecursiveRefinedObjList;
         /// Number of Cost Functions used
         unsigned int mNbCostFunction;
         /// Max number of Cost Functions (dynamically adjusted)
         unsigned int mMaxNbCostFunction;
         /// The objects with the Cost functions
         RefinableObj *mpCostFunctionRefinableObj[20];
         /// The id of the cost functions in each RefinableObj
         CrystVector_int mpCostFunctionId;
         /// The weight associated with each cost function
         CrystVector_double mCostFunctionWeight;
	private:
   #ifdef __WX__CRYST__
   public:
      /// Create a WXCrystObj for this object.
      virtual WXCrystObj* WXCreate(wxWindow*);
      WXCrystObj* WXGet();
      void WXDelete();
      void WXNotifyDelete();
   protected:
      WXGlobalOptimObj *mpWXGlobalOptimObj;
      friend class ObjCryst::WXGlobalOptimObj;
   #endif
};

/// Global Registry for all GlobalOptimObj
extern ObjRegistry<GlobalOptimObj> gGlobalOptimObjRegistry;

}// namespace

#endif //_GLOBALOPTIMOBJ_H
