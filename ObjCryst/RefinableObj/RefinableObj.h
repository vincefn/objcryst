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
*  header file for the RefinablePar and RefinableObj classes
*
*
*/

#ifndef _VFN_REFINABLE_OBJ_H_
#define _VFN_REFINABLE_OBJ_H_

#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <set>

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
   class RefinableObj;
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
      /// Get a Unique id (RefParType::mId)
      void InitId();
      /// the parent for this RefParType (we could easily allow several...)
       const RefParType* mpParent;
      /// The name/description for this parameter type
      const string mName;
      /// The unique number identifying this type
      unsigned long mId;
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
      /// Add a 'child' clock. Whenever a child clock is clicked, it will also click its parent.
      /// This function takes care of adding itself to the list of parent in the children clock.
      void AddChild(const RefinableObjClock &);
      /// remove a child clock. This also tells the child clock to remove the parent.
      void RemoveChild(const RefinableObjClock &);
      /// Add a 'parent' clock. Whenever a clock is clicked, all parent clocks
      /// also are.
      void AddParent(RefinableObjClock &)const;
      /// remove a parent clock
      void RemoveParent(RefinableObjClock &)const;
   private:
      bool HasParent(const RefinableObjClock &) const;
      unsigned long mTick0, mTick1;
      static unsigned long msTick0,msTick1;
      /// List of 'child' clocks, which will click this clock whenever they are clicked.
      std::set<const RefinableObjClock*> mvChild;
      /// List of parent clocks, which will be clicked whenever this one is. This
      /// is a mutable list since reporting to parent clocks does not change
      /// the vlaue of the clock.
      mutable std::set<RefinableObjClock*> mvParent;
};

/** Restraint: generic class for a restraint of a given model. This 
* defines only the category (RefParType) of restraint, and the function
* to access the log(likelihood) associated to this restraint and the current model.
*
* By default: -log(likelihood)=0, and the function must be overloaded
* for "real" restraints.
*
*\note the log(likelihood) \e must always include the normalization term,
* so that variances can also optimized during the maximum likelihood
* optimization. e.g.:
* - \f$ P=\frac{1}{\sqrt{2\pi\sigma^2}}e^{\frac{-(calc-expected)^2}{\sigma^2}}\f$
* - \f$ -\log(P)= \log\left(\sqrt{2\pi\sigma^2}\right)
* + \left(\frac{calc-expected}{\sigma} \right)^2\f$
* 
* forgetting the normalization term would result in making the optimization diverge towards
* infinite variances.
*/
class Restraint
{
   public:
      /// Default constructor, sets RefParType to gpRefParTypeObjCryst
      Restraint();
      ///constructor specifying the type
      Restraint(const RefParType *type);
      virtual ~Restraint();
      virtual const RefParType* GetType()const;
      virtual void SetType(const RefParType *type);
      /// Get -ln(likelihood) for this restraint
      virtual REAL GetLogLikelihood()const;
   private:
      const RefParType *mpRefParType;
};

/** Generic class for parameters of refinable objects.
* These must be continuous.
*
* \todo: define parameters using equations between parameters.
* \todo: for complex objects with lots of parameters, give the 
* possibility to define vectors of parameters, all with the same 
* properties, to reduce memory usage.
*/
class RefinablePar:public Restraint
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
                     REAL *refPar,
                     const REAL min,
                     const REAL max,
                     const RefParType *type,
                     RefParDerivStepModel derivMode=REFPAR_DERIV_STEP_RELATIVE,
                     const bool hasLimits=true,
                     const bool isFixed=false,
                     const bool isUsed=true,
                     const bool isPeriodic=false,
                     const REAL humanScale=1.,
                     REAL period=1.);
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
                     REAL *refPar,
                     const REAL min,
                     const REAL max,
                     const RefParType *type,
                     RefParDerivStepModel derivMode=REFPAR_DERIV_STEP_RELATIVE,
                     const bool hasLimits=true,
                     const bool isFixed=false,
                     const bool isUsed=true,
                     const bool isPeriodic=false,
                     const REAL humanScale=1.,
                     REAL period=1.);
                     
      /// Copy all attributes (limits, flags, etc...) from another RefinablePar object.
      /// This is useful in RefinableObj copy constructors. Everything is copied but the
      /// pointer to the value refined, and the pointer to the clock.
      void CopyAttributes(const RefinablePar&);
      //@}
      
      /// \name Access & change the current value of the parameter
      //@{
         /** of the parameter. Use the The Mutate() and MutateTo() function
         *  to change this value.
         */
         REAL GetValue()const;

         /** of the parameter. Use the The Mutate() and MutateTo() function
         *  to change this value.
         */
         void SetValue(const REAL value);

         /** \brief Current value of parameter, scaled if necessary (for angles) to a
         * human-understandable value.
         */
         const REAL& GetHumanValue() const;

         /** \brief Current value of parameter, scaled if necessary (for angles) to a
         * human-understandable value.
         */
         void SetHumanValue(const REAL&) ;

         /** \brief Add the given amount to the parameter current value.
         *
         * If limit is hit, set to limit.
         * If the limit is hit \e and the parameter is periodic, shift by period to bring
         * back to allowed values.
         *\warning Will throw an exception if the parameter is defined by an equation.
         */
         void Mutate(const REAL mutateValue);
         /**Change the current value to the given one.
         *
         * If the limit is hit, then set to the limit (unless the pameter is periodic,
         * then shift by the period amount back to allowed values).
         *\warning Will throw an exception if the parameter is defined by an equation.
         */
         void MutateTo(const REAL newValue);

         REAL GetSigma()const;
         REAL GetHumanSigma()const;
         void SetSigma(const REAL);
      //@}
      
      /// \name General info 
      //@{
         /// Get the parameter's name
         string GetName()const;
         /// Set the name of the parameter. It should be unique in the RefinableObj.
         void SetName(const string&);
      
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
         void SetIsPeriodic(const bool,REAL period=1);

         /// Human scale for this parameter : for angles, this is equal to 180/pi.
         REAL GetHumanScale()const;
         /// Human scale for this parameter : for angles, this is equal to 180/pi.
         void SetHumanScale(const REAL);
      //@}
      

      /// \name Min, max values
      //@{
         /// Minimum value allowed (if limited or periodic)
         REAL GetMin()const;
         /// Set the Minimum value allowed (if limited)
         void  SetMin(const REAL);
         
         ///Get the minimum value allowed (if limited)
         REAL GetHumanMin()const;
         ///Set the minimum value allowed (if limited)
         void  SetHumanMin(const REAL);

         ///Get the maximum value allowed (if limited)
         REAL GetMax()const;
         ///Get the maximum value allowed (if limited)
         void  SetMax(const REAL);
         
         ///Get the maximum value allowed (if limited)
         REAL GetHumanMax()const;
         ///Get the maximum value allowed (if limited)
         void  SetHumanMax(const REAL);
         
         ///Get the period (if periodic)
         REAL GetPeriod()const;
         ///Set the period value (if periodic)
         void  SetPeriod(const REAL);
      //@}
      
      /// \name Steps during refinement
      //@{
         ///Fixed step to use to compute numerical derivative
         REAL GetDerivStep()const;
         ///Fixed step to use to compute numerical derivative
         void  SetDerivStep(const REAL);

         ///Maximum step to use during Global Optimization algorithms
         REAL GetGlobalOptimStep()const;
         ///Maximum step to use during Global Optimization algorithms
         void  SetGlobalOptimStep(const REAL);
      //@}

      
      #if 0
      /// \name Equations-In development ! ->do not use or even look.
      //@{
         void SetUseEquation(const bool useItOrNot,const REAL c0=0.);
         void SetUseEquation(const bool useItOrNot,const REAL c0,
                             const REAL c1, const RefinablePar &refpar1);

         void SetUseEquation(const bool useItOrNot,const REAL c0,
                             const REAL c1, const RefinablePar &refpar1,
                             const REAL c2, const RefinablePar &refpar2);
         void SetUseEquation(const bool useItOrNot,const REAL c0,
                             const REAL c1, const RefinablePar &refpar1,
                             const REAL c2, const RefinablePar &refpar2,
                             const REAL c3, const RefinablePar &refpar3);
      //@}
      #endif
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
         void SetLimitsAbsolute(const REAL min, const REAL max);
         /// Change the limits for this object, giving relative new limits (eg giving -.1 
         /// and +.1 will set new limits at the current value + min and current value + max)
         /// Thus min should logically be <0 and max >0.
         void SetLimitsRelative(const REAL min, const REAL max);
         /// Change the limits for this object, proportionnaly to the current value.
         /// min should be < 1. and max > 1.
         void SetLimitsProportional(const REAL min, const REAL max);
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
      /// Pointer to the refinable value
      REAL *mpValue;
      /// Hard lower and upper limits.
      REAL mMin,mMax;
      /// Does the refinable parameter need limits (min and max) ?
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
      REAL mPeriod;
      /// Step to use for global method search (simulated annealing,...)
      REAL mGlobalOptimStep;
      /// Step to use for numerical derivative calculation
      REAL mDerivStep;
      /// Model followed for derivation
      RefParDerivStepModel mRefParDerivStepModel;
      /// Calculated sigma on value
      REAL mSigma;
      /// Scale to be used to display 'human' value. This is for angular parameters: the computer
      /// stores values in radians, whil the user only understands degrees. So a scale
      /// factor of 180/pi is necessary.
      REAL mHumanScale;
      #if 0
      // Parameter defined by equations ? :TODO:
         ///Is this parameter deined by an equation ? eg: mValue= c0 +c1*par1.Value+...
         bool mUseEquation;
         /// Max number of other ref. parameters involved in the equation evaluation
         static const int mEquationMaxRefPar=10;
         /// number of other ref. parameters involved in the equation evaluation
         int mEquationNbRefPar;
         /// Coefficient Ci used in equation: Value= C0 + C1 * RefPar1 + C2 * RefPar2 +...
         CrystVector_REAL mEquationCoeff;
         /// Array of pointers to the RefinablePar used in the equation
         const RefinablePar *mEquationRefPar[10];
      #endif
      /// Is there a clock associated with this parameter ? If yes, then it must Click() it
      /// each time it is modified
         bool mHasAssignedClock;
         RefinableObjClock* mpClock;
         
   #ifdef __WX__CRYST__
   public:
      /// Create a WXFieldRefPar representation of the parameter.
      WXCrystObjBasic* WXCreate(wxWindow *parent);
      WXCrystObjBasic* WXGet();
      void             WXDelete();
      void             WXNotifyDelete();
   private:
      WXFieldRefPar * mpWXFieldRefPar;
   #endif
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
*  This class now uses a vector<> approach from the STL.
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
      /// use nothrow=true to avoid having an exception thrown if no object
      /// is found (instead the index returned will be -1)
      long Find(const string &objName, const string& className,
                const bool nothrow=false)const;
      /// Find the number of an object in the registry
      /// The search starts at the *end* of the registry.
      long Find(const T &obj)const;
      /// Last time an object was added or removed from the registry
      const RefinableObjClock& GetRegistryClock()const;
   private:
      /// The registry of objects
      vector<T*> mvpRegistry;
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
      virtual const string& GetClassName() const;
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
      RefinablePar& GetPar(const REAL*);
      /// Access parameter from its adress
      const RefinablePar& GetPar(const REAL*) const;
      
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
      /** Remove a refinable parameter. 
      *
      * This returns an iterator to the next parameter in the vector.
      */ 
      vector<RefinablePar *>::iterator RemovePar(RefinablePar *refPar);
      
      virtual void Print() const;
      
      /** \brief Save the current set of refined values in a new set.
      *
      * \param name : the name associated to this set of values. Names should be unique.
      * \return an number identifying the set of saved values.
      *
      * \warning: there is no limit to the number of parameters sets, so try to
      * release them when you don't need them.
      */
      unsigned long CreateParamSet(const string name="") const;
      /** \brief Erase the param set with the given id, releasing memory.
      */
      void ClearParamSet(const unsigned long id)const;
      /** \brief Save the current set of refined values over a previously-created set
      *of saved values.
      * \param id the number identifying the set of saved values.
      */
      void SaveParamSet(const unsigned long id)const;
      /** \brief Restore a saved set of values.
      *
      * \param id : the number identifying the set.
      * \warning this only affects parameters which are used and not fixed. Others 
      * remain unchanged.
      */
      void RestoreParamSet(const unsigned long id);
      /** \brief Access one save refpar set
      *
      * \param setId : the number identifying the set.
      */
      const CrystVector_REAL& GetParamSet(const unsigned long setId)const;
      /** \brief Access one save refpar set
      *
      * \param setId : the number identifying the set.
      */
      CrystVector_REAL& GetParamSet(const unsigned long setId);
      /** \brief Access the (human) value of one refined parameter in a saved set of parameters
      *
      * \internal
      * \param setId : the number identifying the set.
      * \param parNumber : the number identifying the parameter in the list of refined
      *                    parameters
      * \return if parNumber=5 and setId=37, then the returned value will be the value (scaled
      *if it is an angle) value of the 5th not-fixed parameter in the saved set #37.
      */
      REAL GetParamSet_ParNotFixedHumanValue(const unsigned long setId,const long parNumber)const;
      /** \brief Erase all saved refpar sets
      *
      */
      const void EraseAllParamSet();
      
      /// Change the limits for a given parameter, giving absolute new limits
      void SetLimitsAbsolute(const string &parName, const REAL min, const REAL max);
      /// Change the limits for a category of parameters, giving absolute new limits
      void SetLimitsAbsolute(const RefParType *type, const REAL min, const REAL max);
      /// Change the limits for a given parameter, giving relative new limits (eg giving -.1 
      /// and +.1 will set new limits at the current value + min and current value + max)
      /// Thus min should logically be <0 and max >0.
      void SetLimitsRelative(const string &parName, const REAL min, const REAL max);
      /// Change the limits for a category of parameters, giving relative new limits 
      /// (eg giving -.1 and +.1 will set new limits at the current value + min and 
      /// current value + max). Thus min should logically be <0 and max >0.
      void SetLimitsRelative(const RefParType *type, const REAL min, const REAL max);
      /// Change the limits for a given parameter, proportionnaly to the current value.
      /// min should be < 1. and max > 1.
      void SetLimitsProportional(const string &parName, const REAL min, const REAL max);
      /// Change the limits for a category of parameters, proportionnaly to their current value.
      /// min should be < 1. and max > 1.
      void SetLimitsProportional(const RefParType *type, const REAL min, const REAL max);
      ///Change the maximum step to use during Global Optimization algorithms
      void  SetGlobalOptimStep(const RefParType *type, const REAL step);
      
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
      *
      * \param allowApproximations: if true, then the object can use faster
      * but less precise functions during the optimization. This is useful for
      * global optimization not using derivatives.
      * \param enableRestraints: if true, then restrained parameters will be allowed
      * to go beyond theur hard limits. This implies that the algorithm will take
      * into account the cost (penalty) related to the restraints. Objects which do not
      * use restraints will simply ignore this. WARNING: this parameter may be removed 
      * with the new likelihood scheme.
      */
      virtual void BeginOptimization(const bool allowApproximations=false,
                                     const bool enableRestraints=false);
      /** This should be called by any optimization class at the end of an optimization
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
      *  \param type: restrain the change exclusively to parameters of a given
      *  type (same type or descendant from this RefParType).
      */
      virtual void GlobalOptRandomMove(const REAL mutationAmplitude,
                                       const RefParType *type=gpRefParTypeObjCryst);
      /** Raise a flag, to be sure not to make a random change more than once
      * in each RefinableObj. This calls recursively all sub-objects.
      *
      * This is necessary since one object may be included in several others.
      * This must be called before making a random configuration change on
      * a list of objects.
      */
      void BeginGlobalOptRandomMove();
      
      // Likelihood
         /** Get -log(likelihood) of the current configuration for the object.
         *
         *
         * By default (no likelihood evaluation available), this is equal to 0.
         *
         * This call should not be recursive, it is the task of the algorithm to
         * get the sum of likelihoods for all objects invlolved.
         *
         * \note contrary to the old "Cost Function" approach, with log(Likelihood)
         * there is no 'choice' of cost function, so that it is the task of the 
         * object to give the optimized likelihood (possibly with user options).
         *
         * \warning: this is in under heavy development, so expect changes...
         */
         virtual REAL GetLogLikelihood()const;
      //LSQ functions
         /// Number of LSQ functions
         virtual unsigned int GetNbLSQFunction()const;
         // Get a Cost function name from its id#.
         //virtual const string& GetLSQFunctionName(const unsigned int)const;
         // Get the (short) description of a cost function
         //virtual const string& GetLSQFunctionDescription(const unsigned int)const;
         /// Get the current calculated value for the LSQ function
         virtual const CrystVector_REAL& GetLSQCalc(const unsigned int) const;
         /// Get the observed values for the LSQ function
         virtual const CrystVector_REAL& GetLSQObs(const unsigned int) const;
         /// Get the weight values for the LSQ function
         virtual const CrystVector_REAL& GetLSQWeight(const unsigned int) const;
         /** Get the first derivative values for the LSQ function, for a given
         * parameter. Note that the default method in the base RefinableObj
         * class is to use numerical derivatives, so it should be
         * overridden for better precision.
         * 
         * \todo This should be a const method, and the given RefPar should be const too...
         */
         virtual const CrystVector_REAL& GetLSQDeriv(const unsigned int, RefinablePar&);

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
      // Restraints
         /** Get the restraint cost (overall penalty of all restraints)
         *
         * By default this returns 0, so this \e must be overloaded by any
         * object which actually uses restraint.
         * \note Instead, we could return by default the sum of the restraints,
         * but this is dangerous since we \e want to have derived objects fully
         * responsible for handling restraints.
         */
         virtual REAL GetRestraintCost()const;
         /** Add a new restraint
         *
         * This returns an iterator pointing to the next Restraint in the vector.
         */
         void AddRestraint(Restraint *pNewRestraint);
         /** Remove a restraint from the list of known restraints. This does not
         * delete the Restraint object.
         */
         vector<Restraint*>::iterator RemoveRestraint(Restraint *pRestraint);
      /** During a global optimization, tells the object that the current config is
      * the latest "best" config.
      *
      * This can be used by the object to make more intellingent random moves (use with
      * caution: highly experimental !).
      */
      virtual void TagNewBestConfig()const;
      /// This clocks records _any_ change in the object. See refinableObj::mClockMaster
      const RefinableObjClock& GetClockMaster()const;
   protected:
      /// Find a refinable parameter with a given name
      long FindPar(const string &name) const;
      /// Find a refinable parameter from the adress of its value
      long FindPar(const REAL*) const;
      
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
      
      /// Find a parameter set with a given id (and check if it is there)
      map<unsigned long,pair<CrystVector_REAL,string> >::iterator FindParamSet(unsigned long id)const;

      ///Name for this RefinableObject. Should be unique, at least in the same scope.+
      string mName;
      // Parameters
         /// Vector of pointers to the refinable parameters
         vector<RefinablePar *> mvpRefPar;
      // Restraints
         /// Vector of pointers to the restraints for this object. This excludes
         /// all RefinablePar declared in RefinableObj::mpRefPar, which also
         /// are Restraint.
         vector<Restraint*> mvpRestraint;
         
      //Saved sets of parameters
         /// Map of (index,pointers to arrays) used to save sets of values for all parameters.
         /// Currently there is no limit to the number of saved sets.
         ///
         /// This is mutable since creating/storing a param set does not affect the
         /// 'real' part of the object.
         mutable map<unsigned long,pair<CrystVector_REAL,string> >  mvpSavedValuesSet;
      
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
      /// \internal This true is false if RefinableObj::GlobalOptRandomMove() has been called
      /// since RefinableObj::BeginGlobalOptRandomMove() was called.
      bool mRandomMoveIsDone;
      /// Temporary array used to return derivative values of the LSQ function for given
      /// parameters.
      mutable CrystVector_REAL mLSQDeriv;
      
      /// Master clock, which is changed whenever the object has been altered.
      /// It should be parent to all clocks recording changes in derived classes.
      RefinableObjClock mClockMaster;
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
