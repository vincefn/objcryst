/*
* ObjCryst++ : a Crystallographic computing library in C++
*
*  (c) 2000 Vincent FAVRE-NICOLIN
*  	vincefn@users.sourceforge.net
*
*/
/*
*  header file for the RefinablePar and RefinableObj classes
*
*
*/

#ifndef _VFN_REFINABLE_OBJ_H_
#define _VFN_REFINABLE_OBJ_H_

#include <string>
#include <iostream>
#include "CrystVector/CrystVector.h"
#include "ObjCryst/General.h"
#include "RefinableObj/IO.h"

#ifdef __WX__CRYST__
   class wxWindow;
namespace ObjCryst
{
template<class T> class ObjRegistry;
class RefObjOpt;
template<class T> class RefObjOption;
}
#include "wxCryst/wxRefinableObj.h"
#endif

namespace ObjCryst
{
/// How do we compute steps h for numerical derivative calculation : d=f(x+h)-f(x-h)/h/2
/// either h is fixed (absolute), or relative h=x*derivFactor
enum  RefParDerivStepModel
{
   REFPAR_DERIV_STEP_ABSOLUTE,
   REFPAR_DERIV_STEP_RELATIVE
};

/** \brief class of refinable parameter types.
* 
* This is used to categorize all parameters, to be able to modify
* a complete category of parameters (fix/unfix, set limits,...).
* The types are organized as a tree.
*
* Parameters should be allocated globally in the heap, so we only use
* pointers in the interface.
*
* \note when modifying (fixing, changing limits) for a given RefParType (ie a 'family'
* of parameters), it affects all RefinablePar of this type \e and the parameters
* belonging to the children of this RefParType. eg fixing for the type "gpRefParTypeScatt"
* will fix all the derived postionnal, orientationnal, and population parameters
* for the scatterers.
*
* \remarks In the future, why not use a tree with multiple inheritance ?? It would
* be easy to allow multiple parents... But beware of loops...
*/
class RefParType
{
   public:
      /// Create a \e top parameter type. (in ObjCryst, there is
      /// only one, the "ObjCryst" category.
      RefParType(const string &name);
      /// create a children type.
      RefParType(const RefParType *parent,const string &name);
      ///Destructor
      ~RefParType();
      /// Returns true if the parameter is a descendant of 'type'
      bool IsDescendantFromOrSameAs(const RefParType *type) const;
      /// returns true if the two types are the same.
      bool operator==(const RefParType *parent) const;
      /// Get the name for this parameter
      const string& GetName() const;
   private:
      /// the parent for this RefParType (we could easily allow several...)
       const RefParType* mpParent;
      /// The name/description for this parameter type
      const string mName;
};

/// Top RefParType for the ObjCryst++ library.
extern const RefParType *gpRefParTypeObjCryst;

/// We need to record exactly when refinable objects 
/// have been modified for the last time (to avoid re-computation),
/// and to do that we need a precise time. Since the clock() function is not
/// precise enough (and is architecture-dependant), we use a custom time,
/// which records the number of events in the program which uses the library.
/// This is purely internal, so don't worry about it...
///
/// The clock values have nothing to do with 'time' as any normal person undertands it.
class RefinableObjClock
{
   public:
      RefinableObjClock();
      ~RefinableObjClock();
      bool operator< (const RefinableObjClock &rhs)const;
      bool operator<=(const RefinableObjClock &rhs)const;
      bool operator> (const RefinableObjClock &rhs)const;
      bool operator>=(const RefinableObjClock &rhs)const;
      /// Record an event for this clock (generally, the 'time'
      /// an object has been modified, or some computation has been made)
      void Click();
      /// Reset a Clock to 0, to force an update
      void Reset();
      /// Print clock value. Only for debugging purposes.
      void Print()const;
      /// Print current general clock value. Only for debugging purposes.
      void PrintStatic()const;
   private:
      unsigned long mTick0, mTick1;
      static unsigned long msTick0,msTick1;
};

/** Generic class for parameters of refinable objects.
* These must be continuous.
*
* \todo: define parameters using equations between parameters.
* \todo: define some sort of soft constraint, with the possibility
* of involving several parameters (complex).
* \todo: for complex objects with lots of parameters, give the 
* possibility to define vectors of parameters, all with the same 
* properties, to reduce memory usage.
*/
class RefinablePar
{
   public:
      /// \name Destructor & Constructors
      //@{
      /** Default Constructor
      */
      RefinablePar();
      /** \brief Constructor
      *\par name: the name of the parameter
      *\par refPar: the address of the refined parameter
      *\par min,max: the minimum & maximum value for this parameter. Only used
      * if the parameter is limited.
      *\par type: the type (category) of refinable parameter. This is used to (de)activate
      * groups of parameters.
      *\par derivMode: to compute numerical (partial) derivatives, the step used is either
      * an absolute (fixed) value, or it can be proportionnal to the current value of the
      * parameter.
      *\par hasLimits: if true, then the parameter cannot exceed its limits.
      *\par isFixed: if true, the parameter cannot be refined.
      *\par isUsed: if false, then the parameter does not affect in any way the refined object,
      * and thus is simply ignored and should never appear to the user.
      *\par isPeriodic: if true, then when the parameter exceeds one of its limits, it is
      * shifted by the period (=max-min), in order to be back to the allowed [min,max] range.
      *\par humanScale:this is the scale which should be used to display the value to the
      * end program user. This is mostly used for angles: the values are stored in radians, so
      * that a scale equal to 180/pi must be used for a 'human-understandable' value. Use
      * the RefinablePar::HumanValue() in order to get this value. By default it 
      * is equal to 1.0 (no scaling required).
      */
      RefinablePar(  const string &name,
                     double *refPar,
                     const double min,
                     const double max,
                     const RefParType *type,
                     RefParDerivStepModel derivMode=REFPAR_DERIV_STEP_RELATIVE,
                     const bool hasLimits=true,
                     const bool isFixed=false,
                     const bool isUsed=true,
                     const bool isPeriodic=false,
                     const double humanScale=1.,
                     double period=1.);
      ~RefinablePar();
      /** \brief Constructor
      *\par name: the name of the parameter
      *\par refPar: the address of the refined parameter
      *\par min,max: the minimum & maximum value for this parameter. Only used
      * if the parameter is limited.
      *\par type: the type (category) of refinable parameter. This is used to (de)activate
      * groups of parameters.
      *\par derivMode: to compute numerical (partial) derivatives, the step used is either
      * an absolute (fixed) value, or it can be proportionnal to the current value of the
      * parameter.
      *\par hasLimits: if true, then the parameter cannot exceed its limits.
      *\par isFixed: if true, the parameter cannot be refined.
      *\par isUsed: if false, then the parameter does not affect in any way the refined object,
      * and thus is simply ignored and should never appear to the user.
      *\par isPeriodic: if true, then when the parameter exceeds one of its limits, it is
      * shifted by the period (=max-min), in order to be back to the allowed [min,max] range.
      *\par humanScale:this is the scale which should be used to display the value to the
      * end program user. This is mostly used for angles: the values are stored in radians, so
      * that a scale equal to 180/pi must be used for a 'human-understandable' value. Use
      * the RefinablePar::HumanValue() in order to get this value. By default it 
      * is equal to 1.0 (no scaling required).
      */
      void Init(     const string &name,
                     double *refPar,
                     const double min,
                     const double max,
                     const RefParType *type,
                     RefParDerivStepModel derivMode=REFPAR_DERIV_STEP_RELATIVE,
                     const bool hasLimits=true,
                     const bool isFixed=false,
                     const bool isUsed=true,
                     const bool isPeriodic=false,
                     const double humanScale=1.,
                     double period=1.);
                     
      //@}
      
      /// \name Access & change the current value of the parameter
      //@{
         /** of the parameter. Use the The Mutate() and MutateTo() function
         *  to change this value.
         */
         double GetValue()const;

         /** of the parameter. Use the The Mutate() and MutateTo() function
         *  to change this value.
         */
         void SetValue(const double value);

         /** \brief Current value of parameter, scaled if necessary (for angles) to a
         * human-understandable value.
         */
         const double& GetHumanValue() const;

         /** \brief Current value of parameter, scaled if necessary (for angles) to a
         * human-understandable value.
         */
         void SetHumanValue(const double&) ;

         /** \brief Add the given amount to the parameter current value.
         *
         * If limit is hit, set to limit.
         * If the limit is hit \e and the parameter is periodic, shift by period to bring
         * back to allowed values.
         *\warning Will throw an exception if the parameter is defined by an equation.
         */
         void Mutate(const double mutateValue);
         /**Change the current value to the given one.
         *
         * If the limit is hit, then set to the limit (unless the pameter is periodic,
         * then shift by the period amount back to allowed values).
         *\warning Will throw an exception if the parameter is defined by an equation.
         */
         void MutateTo(const double newValue);

         double GetSigma()const;
         double GetHumanSigma()const;
         void SetSigma(const double);
      //@}
      
      /// \name General info 
      //@{
         ///Parameter type
         const RefParType* GetType()const;
         ///Parameter type
         void SetType(const RefParType*);
         ///Get the parameter's name
         string GetName()const;
      
         void Print() const;
         
         bool IsFixed()const;
         void SetIsFixed(const bool);

         bool IsLimited()const;
         void SetIsLimited(const bool);

         /// Is the parameter used (if not, it is simply irrelevant in the model) ?
         bool IsUsed()const;
         /// Is the parameter used (if not, it is simply irrelevant in the model) ?
         void SetIsUsed(const bool);

         bool IsPeriodic()const;
         void SetIsPeriodic(const bool,double period=1);

         /// Human scale for this parameter : for angles, this is equal to 180/pi.
         double GetHumanScale()const;
         /// Human scale for this parameter : for angles, this is equal to 180/pi.
         void SetHumanScale(const double);
      //@}
      

      /// \name Min, max values
      //@{
         /// Minimum value allowed (if limited or periodic)
         double GetMin()const;
         /// Set the Minimum value allowed (if limited)
         void  SetMin(const double);
         
         ///Get the minimum value allowed (if limited)
         double GetHumanMin()const;
         ///Set the minimum value allowed (if limited)
         void  SetHumanMin(const double);

         ///Get the maximum value allowed (if limited)
         double GetMax()const;
         ///Get the maximum value allowed (if limited)
         void  SetMax(const double);
         
         ///Get the maximum value allowed (if limited)
         double GetHumanMax()const;
         ///Get the maximum value allowed (if limited)
         void  SetHumanMax(const double);
			
         ///Get the period (if periodic)
         double GetPeriod()const;
         ///Set the period value (if periodic)
         void  SetPeriod(const double);
      //@}
      
      /// \name Steps during refinement
      //@{
         ///Fixed step to use to compute numerical derivative
         double GetDerivStep()const;
         ///Fixed step to use to compute numerical derivative
         void  SetDerivStep(const double);

         ///Maximum step to use during Global Optimization algorithms
         double GetGlobalOptimStep()const;
         ///Maximum step to use during Global Optimization algorithms
         void  SetGlobalOptimStep(const double);
      //@}

      
      
      /// \name Equations-In development ! ->do not use or even look.
      //@{
         void SetUseEquation(const bool useItOrNot,const double c0=0.);
         void SetUseEquation(const bool useItOrNot,const double c0,
                             const double c1, const RefinablePar &refpar1);

         void SetUseEquation(const bool useItOrNot,const double c0,
                             const double c1, const RefinablePar &refpar1,
                             const double c2, const RefinablePar &refpar2);
         void SetUseEquation(const bool useItOrNot,const double c0,
                             const double c1, const RefinablePar &refpar1,
                             const double c2, const RefinablePar &refpar2,
                             const double c3, const RefinablePar &refpar3);
      //@}
             
      /// \name Parameter's Clock
      //@{
      /// \internal
      /// Assign a clock to this parameter. Any time this parameter is modified,
      /// the clock will be ticked !
      void AssignClock(RefinableObjClock &clock);
      //@}
      
      /// \name Change Limits
      //@{
         /// Change the limits for this object, giving absolute new limits
         void SetLimitsAbsolute(const double min, const double max);
         /// Change the limits for this object, giving relative new limits (eg giving -.1 
         /// and +.1 will set new limits at the current value + min and current value + max)
         /// Thus min should logically be <0 and max >0.
         void SetLimitsRelative(const double min, const double max);
         /// Change the limits for this object, proportionnaly to the current value.
         /// min should be < 1. and max > 1.
         void SetLimitsProportional(const double min, const double max);
      //@}

      /** \brief XMLOutput to stream in well-formed XML 
      *
      * this will save the fixed & limited flags, as well as limits
      * \param name the name to use instead of the RefPar name.
      */
      void XMLOutput(ostream &os,const string &name,int indent=0)const;
      /** \brief XMLOutput to stream in well-formed XML.
      *
      * this will save the fixed & limited flags, as well as limits.
      * In this function the name used is that of the RefPar.
      */
      void XMLOutput(ostream &os,int indent=0)const;
      /** \brief XMLInput From stream
      *
      */
      void XMLInput(istream &is,const XMLCrystTag &tag);
   private:
      /// Click the Clock ! to telle the RefinableObj it has been modified.
      void Click();
      ///name of the refinable parameter
      string mName;
      ///Type of refined variable
      const RefParType *mpRefParType;
      /// Pointer to the refinable value
      double *mValue;
      /// Min value
      double mMin;
      /// Max value
      double mMax;
      /// Does the refinable parameter need limits (min,max) ?
      bool mHasLimits;
      /// is the parameter currently fixed ?
      bool mIsFixed;
      /// Is the parameter currently used ?
      bool mIsUsed;
      /// Is the parameter periodic ? If this is the case, then when using the 
      /// RefinablePar::Mutate() function, if the parameter goes beyond its limits,
      /// it will be shifted by the value of its period.
      bool mIsPeriodic;
      /// Period value (if relevant)
      double mPeriod;
      /// Step to use for global method search (simulated annealing,...)
      double mGlobalOptimStep;
      /// Step to use for numerical derivative calculation
      double mDerivStep;
      /// Model followed for derivation
      RefParDerivStepModel mRefParDerivStepModel;
      /// Calculated sigma on value
      double mSigma;
      /// Scale to be used to display 'human' value. This is for angular parameters: the computer
      /// stores values in radians, whil the user only understands degrees. So a scale
      /// factor of 180/pi is necessary.
      double mHumanScale;
      
      // Parameter defined by equations ?
         ///Is this parameter deined by an equation ? eg: mValue= c0 +c1*par1.Value+...
         bool mUseEquation;
         /// Max number of other ref. parameters involved in the equation evaluation
         static const int mEquationMaxRefPar=10;
         /// number of other ref. parameters involved in the equation evaluation
         int mEquationNbRefPar;
         /// Coefficient Ci used in equation: Value= C0 + C1 * RefPar1 + C2 * RefPar2 +...
         CrystVector_double mEquationCoeff;
         /// Array of pointers to the RefinablePar used in the equation
         const RefinablePar *mEquationRefPar[10];
      
      /// Is there a clock associated with this parameter ? If yes, then it must Click() it
      /// each time it is modified
         bool mHasAssignedClock;
         RefinableObjClock* mpClock;
         
   friend class RefinableObj;
};
/** Base class for options
*  
*/
class RefObjOpt
{
   public:
      /** Constructor for the option
      * \param obj: the 
      */
      RefObjOpt();
      virtual ~RefObjOpt();
      void Init(const int nbChoice,const string *name,
                const string *choiceNames);
      int GetNbChoice()const;
      int GetChoice()const;
      virtual void SetChoice(const int choice);
      void SetChoice(const string &choiceName);
      const string& GetName()const;
      const string& GetClassName()const;
      const string& GetChoiceName(const int i)const;
      const RefinableObjClock& GetClock()const;
      /** \brief XMLOutput to stream in well-formed XML.
      *
      * In this function the name used is that of the Option.
      */
      void XMLOutput(ostream &os,int indent=0)const;
      /** \brief XMLInput From stream
      *
      */
      void XMLInput(istream &is,const XMLCrystTag &tag);
   protected:
      /// Number of different choice possible for this option
      int mNbChoice;
      /// Current value
      int mChoice;
      /// (short) Name for this option. Should be statically stored in the class
      /// using the option
      const string* mpName;
      /// Names corresponding to each possible value of this option (Human-understandable).
      /// Should be statically stored in the class using the option.
      const string* mpChoiceName;
      /// The clock associated to this option
      RefinableObjClock mClock;
   #ifdef __WX__CRYST__
   public:
      WXCrystObjBasic* WXCreate(wxWindow *parent);
      WXCrystObjBasic* WXGet();
      void             WXDelete();
      void             WXNotifyDelete();
   private:
      WXFieldOption * mpWXFieldOption;
   #endif
};


/**
*  Class for options of RefinableObj, templated so that we can warn
* the object that something has been changed. NOT USED SO FAR.
*/
template<class T> class RefObjOption:public RefObjOpt
{
   public:
      /** Constructor for the option
      * \param obj: the 
      */
      RefObjOption(T* obj);
      ~RefObjOption();
      void Init(const int nbChoice,const string *name,
                const string *choiceNames,
                void (T::*fp)(const int));
      virtual void SetChoice(const int choice);
   protected:
   private:
      /// The object which uses this option
      T* mpObj;
      /// The pointer to the member function to be used when the choice is changed, to notify
      /// immediately the object. If null, the value is just recorded and no notification
      /// is done.
      void (T::*mfpSetNewValue)(const int);
};

/** Object Registry
*
*  This class is used to keep a list of all object of a given class at the global
*  level, or inside another object. This is primarily aimed for the derivative
*  of the RefinableObj class but it
*  can be used for \e any class that has GetName() and GetClassName() function.
*
*  \warning the order of the objects in the registry can change (every time an object
*  is de-registered).
*
* \todo (?) create two derived classes with the same interface, one which is a const
* registry (the 'client' registry for RefinableObj), and one which has a non-const
* access to the registered objects (the 'sub-objects' in RefinableObj).
*/
template<class T> class ObjRegistry
{
   public:
      ObjRegistry();
      ObjRegistry(const string &name);
      ~ObjRegistry();
      /// Register a new object. Already registered objects are skipped.
      void Register(T &obj);
      /// De-register an object.
      void DeRegister(T &obj);
      /// De-register an object from its name.
      void DeRegister(const string &objName);
      /// De-register \e all objects from the list.
      void DeRegisterAll();
      /// Delete all objects in the registry.. Use with caution !!
      void DeleteAll();
      /** Get object #i in the registry.
      * \internal
      * Use with caution. The order of the objects
      * changes as objects are added and removed.
      */
      T& GetObj(const unsigned int i);
      /** Get object #i in the registry.
      * \internal
      * Use with caution. The order of the objects
      * changes as objects are added and removed.
      */
      const T& GetObj(const unsigned int i)const;
      /// Get an object from its name in the registry.
      /// The search starts at the *end* of the registry.
      T& GetObj(const string &objName);
      /// Get an object from its name in the registry.
      /// The search starts at the *end* of the registry.
      const T& GetObj(const string &objName)const;
      /// Get an object from its name in the registry.
      /// The search starts at the *end* of the registry.
      /// Also check the class of the object.
      T& GetObj(const string &objName, const string& className);
      /// Get an object from its name in the registry.
      /// The search starts at the *end* of the registry.
      /// Also check the class of the object.
      const T& GetObj(const string &objName, const string& className)const;
      /// Get the index of an object in the registry, from its name
      /// Warning: it can change if an object is removed from the registry
      long GetNb()const;
      void Print()const;
      void SetName(const string &);
      const string& GetName()const;
      /// Find the number of an object in the registry from its name (slow !)
      /// The search starts at the *end* of the registry.
      long Find(const string &objName)const;
      /// Find the number of an object in the registry from its name (slow !)
      /// The search starts at the *end* of the registry.
      /// Also check the class of the object (inheritance...).
      long Find(const string &objName, const string& className)const;
      /// Find the number of an object in the registry
      /// The search starts at the *end* of the registry.
      long Find(const T &obj)const;
      /// Last time an object was added or removed from the registry
		const RefinableObjClock& GetRegistryClock()const;
   private:
      /// The registry
      T** mpRegistry;
      /// Number of registered objects
      unsigned long mNbRegistered;
      /// Max number of registered objects. Dynamically allocated
      unsigned long mMaxNbRegistered;
      /// Name of this registry
      string mName;
		/// Last time an object was added or removed
		RefinableObjClock mListClock;
		
   #ifdef __WX__CRYST__
   public:
      WXRegistry<T>* WXCreate(wxWindow *parent);
      void             WXDelete();
      void             WXNotifyDelete();
   private:
      WXRegistry<T> * mpWXRegistry;
   #endif
};

/// Register a new object in a registry, and recursively
/// include all included (sub)objects. 
template<class T> void RefObjRegisterRecursive(T &obj,ObjRegistry<T> &reg);
/// Get the last time any object was added in the recursive list of objects.
void GetSubRefObjListClockRecursive(ObjRegistry<RefinableObj> &reg,RefinableObjClock &clock);

/** \brief Generic Refinable Object
*
* This is basically a list of refinable parameters, with other basic common properties
* such as a name,. etc...
* This allows optimization/refinement algorithms to access the parameters without
* knowing what the object really is.
*
* \todo Define more clearly which operations are recursive (ie also modify sub-objects).
*/
class RefinableObj
{
   public:
      /// Constructor
      RefinableObj();
      /// Constructor. Using internalUseOnly=true will avoid registering the
		/// the object to any registry, and thus (for example) no display will be created,
		/// nor will this object be automatically be saved.
      RefinableObj(const bool internalUseOnly);
      /// Defined not implemented... Should never be called
		/// (copying the refinable parameters would allow you to modify the 
		/// input object).
		/// Use the default constructor and RefinableObj::AddPar(RefinableObj&) instead.
      RefinableObj(const RefinableObj &old);
      /// Destructor
      virtual ~RefinableObj();
      /// Name for this class ("RefinableObj", "Crystal",...). This is only useful
      /// to distinguish different classes when picking up objects from the
      /// RefinableObj Global Registry 
      virtual const string GetClassName() const;
      /// Name of the object
      virtual const string& GetName() const;
      /// Name of the object
      virtual void SetName(const string &name);
      /** Defined not implemented... Should never be called
		*/
      void operator=(const RefinableObj &old);
      
      /// Find which parameters are used and \b not fixed, for a refinement /optimization.
      /// This \b must be called before any refinement...
      void PrepareForRefinement() const;
      /// Fix All parameters
      void FixAllPar();
      /// UnFix All parameters
      void UnFixAllPar();
      /// Fix/un-fix one parameter from its #.
      void SetParIsFixed(const long parIndex,const bool fix);
      /// Fix/un-fix one parameter from its name.
      void SetParIsFixed(const string& parName,const bool fix);
      /// Fix/un-fix one family of parameters
      void SetParIsFixed(const RefParType *type,const bool fix);
      /// Set whether a parameter is used
      void SetParIsUsed(const string& parName,const bool use);
      /// Set whether a family of parameters is used
      void SetParIsUsed(const RefParType *type,const bool use);
      /// Total number of refinable parameter in the object. Note that some
      /// may be actually fixed or even not used !!
      /// For refinement use PrepareForRefinement(),  NbRefParNotFixed(), and
      /// ParNotFixed(i)
      long GetNbPar()const;
      /// Total number of non-fixed parameters. Is initialized by PrepareForRefinement()
      long GetNbParNotFixed()const;
      
      /// Access all parameters in the order they were inputted
      RefinablePar& GetPar(const long i);
      /// Access all parameters in the order they were inputted
      const RefinablePar& GetPar(const long i) const;
      
      /// Access all parameters from their name
      RefinablePar& GetPar(const string & name);
      ///Access all parameters from their name
      const RefinablePar& GetPar(const string & name) const;
      
      /// Access parameter from its adress
      RefinablePar& GetPar(const double*);
      /// Access parameter from its adress
      const RefinablePar& GetPar(const double*) const;
      
      /// Access all parameters in the order they were inputted,
      /// skipping fixed parameters. Must call PrepareForRefinement() before !
      RefinablePar& GetParNotFixed(const long i);
      /// Access all parameters in the order they were inputed,
      /// skipping fixed parameters. Must call PrepareForRefinement() before !
      const RefinablePar& GetParNotFixed(const long i)const;
		/** Add a refinable parameter. The parameter is copied, so 
		* it need only be allocated temporarily.
		*
		* \deprecated Use the next function, which supplies the parameter as
		* a pointer, and avoids a useless copy.
		*/ 
      void AddPar(const RefinablePar &newRefPar);
		/** Add a refinable parameter. The parameter is \e not copied, so 
		* it should be allocated in the heap.
		*
		*/ 
      void AddPar(RefinablePar *newRefPar);
		/** Add all the parameters in another RefinableObj. Parameters
		* are \not copied, so they should be allocated in the heap.
		* 
		* \warning If a copy of another RefinableObj parameter list is made,
		* such as in the OptimizationObj class, make sure that upon deletion
		* of this object the parameters will not be destroyed. To do this
		* use RefinableObj::SetDeleteRefParInDestructor(false).
		*/
      void AddPar(RefinableObj &newRefParList);
      
      virtual void Print() const;
      
      /** \brief Save the current set of refined values in a new set.
      *
      * \param name : the name associated to this set of values. Names should be unique.
      * \return an number identifying the set of saved values.
      */
      long CreateParamSet(const string name="") const;
      /** \brief Save the current set of refined values over a previously-created set
      *of saved values.
      * \param id the number identifying the set of saved values.
      */
      void SaveParamSet(const long id)const;
      /** \brief Restore a saved set of values.
      *
      * \param id : the number identifying the set.
      * \warning this only affects parameters which are used and not fixed. Others 
      * remain unchanged.
      */
      void RestoreParamSet(const long id);
      /** \brief Access one save refpar set
      *
      * \param setId : the number identifying the set.
      */
      const CrystVector_double& GetParamSet(const long setId)const;
      /** \brief Access one save refpar set
      *
      * \param setId : the number identifying the set.
      */
      CrystVector_double& GetParamSet(const long setId);
      /** \brief Access the (human) value of one refined parameter in a saved set of parameters
      *
      * \internal
      * \param setId : the number identifying the set.
      * \param parNumber : the number identifying the parameter in the list of refined
      *                    parameters
      * \return if parNumber=5 and setId=37, then the returned value will be the value (scaled
      *if it is an angle) value of the 5th not-fixed parameter in the saved set #37.
      */
      double GetParamSet_ParNotFixedHumanValue(const long setId,const long parNumber)const;
      /** \brief Erase all saved refpar sets
      *
      */
      const void EraseAllParamSet();
      
      /// Change the limits for a given parameter, giving absolute new limits
      void SetLimitsAbsolute(const string &parName, const double min, const double max);
      /// Change the limits for a category of parameters, giving absolute new limits
      void SetLimitsAbsolute(const RefParType *type, const double min, const double max);
      /// Change the limits for a given parameter, giving relative new limits (eg giving -.1 
      /// and +.1 will set new limits at the current value + min and current value + max)
      /// Thus min should logically be <0 and max >0.
      void SetLimitsRelative(const string &parName, const double min, const double max);
      /// Change the limits for a category of parameters, giving relative new limits 
      /// (eg giving -.1 and +.1 will set new limits at the current value + min and 
      /// current value + max). Thus min should logically be <0 and max >0.
      void SetLimitsRelative(const RefParType *type, const double min, const double max);
      /// Change the limits for a given parameter, proportionnaly to the current value.
      /// min should be < 1. and max > 1.
      void SetLimitsProportional(const string &parName, const double min, const double max);
      /// Change the limits for a category of parameters, proportionnaly to their current value.
      /// min should be < 1. and max > 1.
      void SetLimitsProportional(const RefParType *type, const double min, const double max);
      
      /// Access to the registry of RefinableObj used by this object
      ObjRegistry<RefinableObj>& GetSubObjRegistry();
      /// Access to the registry of RefinableObj used by this object
      const ObjRegistry<RefinableObj>& GetSubObjRegistry()const;
      
      /// Register a new object using this object
      /// \todo : the clients should be const, but are not... This need to be fixed...
      virtual void RegisterClient(RefinableObj &)const;
      /// Deregister an object (which not any more) using this object
      virtual void DeRegisterClient(RefinableObj &)const;
      
      /// Is the object being refined ? (Can be refined by one algorithm at a time only.)
      bool IsBeingRefined()const;
      
      /** This should be called by any optimization class at the begining of an optimization
      *
      * This will also check that everything is ready, eg call the RefinableObj::Prepare()
      * function. This also affects all sub-objects.
      * \note this may be called several time for some objects which are used by several
      * other objects.
      */
      virtual void BeginOptimization();
      /** This should be called by any optimization class at the begining of an optimization
      *
      * This also affects all sub-objects.
      * \note this may be called several time for some objects which are used by several
      * other objects.
      */
      virtual void EndOptimization();
      
      /// Randomize Configuration (before a global optimization). This
      /// Affects only parameters which are limited and not fixed.
      /// The randomization also affects all sub-objects (recursive).
      virtual void RandomizeConfiguration();
      
      /** Make a random move of the current configuration.
      *
      *  This is for global optimization algorithms. the moves for each
      *  parameter are less than their global optimization step, multiplied
      *  by the mutation amplitude.
      *
      *  \warning: this makes a random move for the parameter declared
      *  for this object, and it is the duty of the object to decide whether
      *  the included objects should be moved and how. (eg an algorithm should
      *  only call for a move with the top object, and this object decides how
      *  he and his sub-objects moves). By default (RefinableObj implementation) 
      *  all included objects are moved recursively.
      *
      *  RefinableObj::
      *  \param mutationAmplitude: multiplier for the maximum move amplitude,
      *  for all parameters
      */
      virtual void GlobalOptRandomMove(const double mutationAmplitude);
      
      //Cost functions
         /// Number of Cost functions
         virtual unsigned int GetNbCostFunction()const;
         /// Get a Cost function name from its id#.
         virtual const string& GetCostFunctionName(const unsigned int)const;
         /// Get the (short) description of a cost function
         virtual const string& GetCostFunctionDescription(const unsigned int)const;
         /// Get the current value of a cost function
         /// this should be const...
         virtual double GetCostFunctionValue(const unsigned int);
      //LSQ functions
         /// Number of LSQ functions
         virtual unsigned int GetNbLSQFunction()const;
         // Get a Cost function name from its id#.
         //virtual const string& GetLSQFunctionName(const unsigned int)const;
         // Get the (short) description of a cost function
         //virtual const string& GetLSQFunctionDescription(const unsigned int)const;
         /// Get the current calculated value for the LSQ function
         virtual const CrystVector_double& GetLSQCalc(const unsigned int) const;
         /// Get the observed values for the LSQ function
         virtual const CrystVector_double& GetLSQObs(const unsigned int) const;
         /// Get the weight values for the LSQ function
         virtual const CrystVector_double& GetLSQWeight(const unsigned int) const;

      /// Re-init the list of refinable parameters, removing all parameters.
      /// This does \e not delete the RefinablePar if 
		/// RefinableObj::mDeleteRefParInDestructor is false
      void ResetParList();
      
      /** \brief Output to stream in well-formed XML 
      *
      * \todo Use inheritance.. as for XMLInputTag()...
      */
      virtual void XMLOutput(ostream &os,int indent=0)const;
      /** \brief Input From stream
      *
      * \todo Add an bool XMLInputTag(is,tag) function to recognize all the tags
      * from the stream. So that each inherited class can use the XMLInputTag function
      * from its parent (ie take advantage of inheritance). The children class 
      * would first try to interpret the tag, then if unsuccessful would pass it to
      * its parent (thus allowing overloading), etc...
      */
      virtual void XMLInput(istream &is,const XMLCrystTag &tag);
      //virtual void XMLInputOld(istream &is,const IOCrystTag &tag);
      /// If there is an interface, this should be automatically be called each
      /// time there is a 'new, significant' configuration to report.
      virtual void UpdateDisplay()const;
      //Options
         /// Number of Options for this object
         unsigned int GetNbOption()const;
         /// Access to the options
         RefObjOpt& GetOption(const unsigned int i);
         /// const access to the options
         const RefObjOpt& GetOption(const unsigned int i)const;
		// Genetic
			/** \brief Get the gene group assigned to each parameter.
			*
			* Each parameter (a \e gene in terms of genetic algorithms)
			* can be assigned to a gene group. Thus when mating two configurations,
			* genes will be exchanged by groups. By default (in the base RefinabeObj class),
			* each parameter is alone in its group. Derived classes can group genes
			* for a better s** life.
			*
			* The number identifying a gene group only has a meaning in a given
			* object. It can also change on subsequent calls, and thus is not unique.
			*
			* \param obj the \RefinableObj, supplied by an algorithm class (OptimizationObj,..),
			* which contains a list of parameters, some of which (but possibly all or none)
			* are parameters belonging to this object.
			* \param groupIndex a vector of unsigned integers, one for each parameter in the
			* input object, giving an unsigned integer value as gene group index.
			* At the beginning this vector should contain only zeros (no group assigned).
			* \param firstGroup this is the number of groups which have already been assigned,
			* plus one. The gene groups returned by this object will start from this
			* value, and increment \b firstGroup for each gene group used, so that
			* different RefinableObj cannot share a gene group.
			* \note this function is not optimized, and should only be called at the beginning
			* of a refinement.
			*/
			virtual void GetGeneGroup(const RefinableObj &obj, 
											  CrystVector_uint & groupIndex,
											  unsigned int &firstGroup) const;
		/** Set this object not to delete its list of parameters when destroyed.
		*
		* This is used for the RefinableObj in algorithms objects (OptimizationObj),
		* which only hold copies of parameters from the refined objects.
		*/
		void SetDeleteRefParInDestructor(const bool b);
		/** What was the last time a RefinablePar was added/removed ?
		*
		*/
		const RefinableObjClock& GetRefParListClock()const;
   protected:
      /// Find a refinable parameter with a given name
      long FindPar(const string &name) const;
      /// Find a refinable parameter from the adress of its value
      long FindPar(const double*) const;
      
      /// \internal Add an object in the registry of used objects.
      void AddSubRefObj(RefinableObj &);
      /// \internal Remove an object in the registry of used objects.
      void RemoveSubRefObj(RefinableObj &);
      
      /// \internal Init the seed for random number generation from current time.
      void InitRandomSeedFromTime()const;
      /// \internal Add an option for this parameter
      void AddOption(RefObjOpt *opt);
      /// \internal Prepare everything (if necessary) for an optimization/calculation.
      virtual void Prepare();
      
      ///Name for this RefinableObject. Should be unique, at least in the same scope.+
      string mName;
      ///Array of refinable parameters
      RefinablePar **mpRefPar;
      ///Number of refinable parameters
      long mNbRefPar;
      ///Maximum number of refinable parameters (array size-dynamically allocated)
      long mMaxNbRefPar;
      
      //Saved sets of parameters
         ///Max number of saved sets (memory is dynamically allocated...)
         static const int mMaxNbSavedSets=1000;
         ///Array of pointers to arrays used to save sets of values for all parameters
         mutable CrystVector_double **mpSavedValuesSet;
         ///Names associated to the saved values sets
         mutable string **mpSavedValuesSetName;
         ///Is the set associated with (id) currently used ?
         mutable CrystVector_bool mSavedValuesSetIsUsed;
      
      // Used during refinements, initialized by PrepareForRefinement()
         /// Total of not-fixed parameters
         mutable long mNbRefParNotFixed;
         /// Index of not-fixed parameters
         mutable CrystVector_long mRefparNotFixedIndex;
         /// Is the object being refined ?
         bool mIsbeingRefined;
      
      /// Registry of RefinableObject needed for this object (owned by this object or not)
         ObjRegistry<RefinableObj> mSubObjRegistry;
      /// Registry of RefinableObject using this object.
      /// This is mutable so that client can modify it (kludge?)
         mutable ObjRegistry<RefinableObj> mClientObjRegistry;
      
      // Options for this object      
         /// List of options for this object. Note that these are just references,
         /// to options allocated by the object, to have a simple global access to 
         /// all options
         ObjRegistry<RefObjOpt> mOptionRegistry;
		/// If true (the default), then all RefinablePar will be deleted when the 
		/// the object is deleted. The opposite option (false) should only be used
		/// in RefinableObj holding 'copies' of other objects, such as in algorithms.
      bool mDeleteRefParInDestructor;
		/// Last time the RefinableParList was modified (a parameter added or removed).
		RefinableObjClock mRefParListClock;
   #ifdef __WX__CRYST__
   public:
      /// Create a WXCrystObj for this object. Only a generic WXCrystObj pointer is kept.
      virtual WXCrystObjBasic* WXCreate(wxWindow*);
      WXCrystObjBasic* WXGet();
      void WXDelete();
      void WXNotifyDelete();
   protected:
      WXCrystObjBasic *mpWXCrystObj;
   #endif
};
/// Get the last time any RefinablePar was added in a recursive list of objects.
void GetRefParListClockRecursive(ObjRegistry<RefinableObj> &reg,RefinableObjClock &clock);

/// Global Registry for all RefinableObj
extern ObjRegistry<RefinableObj> gRefinableObjRegistry;
/// This is a special registry for 'top' object for an optimization. In
/// the ObjCryst++ class, this currently includes Crystal, PowderPattern and
/// DiffractionDataSingleCrystal.
extern ObjRegistry<RefinableObj> gTopRefinableObjRegistry;
}//namespace

#endif// _VFN_REFINABLE_OBJ_H_
