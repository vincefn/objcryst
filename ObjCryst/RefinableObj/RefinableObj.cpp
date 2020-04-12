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
*  source file for the RefinablePar and RefinableObj classes
*
*/
#include <ctime>
#include <boost/format.hpp>
#include "ObjCryst/RefinableObj/RefinableObj.h"
#include "ObjCryst/Quirks/VFNStreamFormat.h"
#include "ObjCryst/Quirks/VFNDebug.h"
#ifdef __WX__CRYST__
   #include "ObjCryst/wxCryst/wxRefinableObj.h"
   #undef GetClassName // Conflict from wxMSW headers ? (cygwin)
#endif
#include <algorithm>

#define POSSIBLY_UNUSED(expr) (void)(expr)

namespace ObjCryst
{
//######################################################################
//
//      RefParType
//
//######################################################################
RefParType::RefParType(const string &name):
mpParent(0),mName(name),mId(0)
{
   this->InitId();
}

RefParType::RefParType(const RefParType *parent,const string &name):
mpParent(parent),mName(name),mId(0)
{
   this->InitId();
}

RefParType::~RefParType(){};

bool RefParType::IsDescendantFromOrSameAs(const RefParType *type) const
{
   VFN_DEBUG_MESSAGE("RefParType::IsDescendantFromOrSameAs(RefParType*): "<<this<<" : "<<mpParent,1)
   if(type->mId==mId) return true;
   if(0==mpParent) return false;
   return mpParent->IsDescendantFromOrSameAs(type);
}

bool RefParType::operator==(const RefParType *type) const
{
   if(this==type) return true;
   return false;
}
const string& RefParType::GetName() const{ return mName;}

void RefParType::InitId()
{
   static unsigned long nbRefParType=0;
   mId=nbRefParType++;
}

const RefParType *gpRefParTypeObjCryst=0;
long NiftyStaticGlobalObjectsInitializer_RefinableObj::mCount=0;
//######################################################################
//
//
//
//######################################################################

unsigned long RefinableObjClock::msTick0=0;
unsigned long RefinableObjClock::msTick1=0;
RefinableObjClock::RefinableObjClock()
{
   //this->Click();
   mTick0=0;
   mTick1=0;
}
RefinableObjClock::~RefinableObjClock()
{
   // first copy & clear sets to avoid possible loops in RemoveChild()
   set<const RefinableObjClock*> vChild=mvChild;
   set<RefinableObjClock*> vParent=mvParent;
   mvChild.clear();
   mvParent.clear();
   for(std::set<const RefinableObjClock*>::iterator pos=vChild.begin();
       pos!=vChild.end();++pos) (*pos)->RemoveParent(*this);
   for(std::set<RefinableObjClock*>::iterator pos=vParent.begin();
       pos!=vParent.end();++pos) (*pos)->RemoveChild(*this);
}

bool RefinableObjClock::operator< (const RefinableObjClock &rhs)const
{
   if(mTick1<rhs.mTick1) return true;
   if(mTick1==rhs.mTick1) {if(mTick0<rhs.mTick0) return true;}
   return false;
}
bool RefinableObjClock::operator<=(const RefinableObjClock &rhs)const
{
   if(mTick1<rhs.mTick1) return true;
   if(mTick1==rhs.mTick1) if(mTick0<=rhs.mTick0) return true;
   return false;
}
bool RefinableObjClock::operator> (const RefinableObjClock &rhs)const
{
   if(mTick1>rhs.mTick1) return true;
   if(mTick1==rhs.mTick1) if(mTick0>rhs.mTick0) return true;
   return false;
}
bool RefinableObjClock::operator>=(const RefinableObjClock &rhs)const
{
   if(mTick1>rhs.mTick1) return true;
   if(mTick1==rhs.mTick1) if(mTick0>=rhs.mTick0) return true;
   return false;
}
void RefinableObjClock::Click()
{
   //return;
   if(++msTick0==0) ++msTick1;//Update ObjCryst++ static event counter
   mTick0=msTick0;
   mTick1=msTick1;
   for(std::set<RefinableObjClock*>::iterator pos=mvParent.begin();
       pos!=mvParent.end();++pos) (*pos)->Click();
   VFN_DEBUG_MESSAGE("RefinableObjClock::Click():"<<mTick1<<":"<<mTick0<<"(at "<<this<<")",0)
   //this->Print();
}
void RefinableObjClock::Reset()
{
   mTick0=0;
   mTick1=0;
}
void RefinableObjClock::Print()const
{
   cout <<"Clock():"<<mTick1<<":"<<mTick0;
   VFN_DEBUG_MESSAGE_SHORT(" (at "<<this<<")",4)
   cout <<endl;
}
void RefinableObjClock::PrintStatic()const
{
   cout <<"RefinableObj class Clock():"<<msTick1<<":"<<msTick0<<endl;
}
void RefinableObjClock::AddChild(const RefinableObjClock &clock)
{mvChild.insert(&clock);clock.AddParent(*this);this->Click();}
void RefinableObjClock::RemoveChild(const RefinableObjClock &clock)
{
   const unsigned int i = mvChild.erase(&clock);  POSSIBLY_UNUSED(i);
   VFN_DEBUG_MESSAGE("RefinableObjClock::RemoveChild():"<<i,5)
   clock.RemoveParent(*this);
   this->Click();
}
void RefinableObjClock::AddParent(RefinableObjClock &clock)const
{
   // First check for loop
   if(&clock==this)
      throw ObjCrystException("RefinableObjClock::AddParent(..) child == Parent !!");
   if(clock.HasParent(*this)==true)
      throw ObjCrystException("RefinableObjClock::AddParent(..) Loop in clock tree !!");
   mvParent.insert(&clock);
}
void RefinableObjClock::RemoveParent(RefinableObjClock &clock)const
{
   // avoid warnings about unused i when not debugging.
   const unsigned int i = mvParent.erase(&clock);  POSSIBLY_UNUSED(i);
   VFN_DEBUG_MESSAGE("RefinableObjClock::RemoveParent():"<<i,5)
}

void RefinableObjClock::operator=(const RefinableObjClock &rhs)
{
   mTick0=rhs.mTick0;
   mTick1=rhs.mTick1;
   for(std::set<RefinableObjClock*>::iterator pos=mvParent.begin();
       pos!=mvParent.end();++pos) if( (*this) > (**pos) ) **pos = *this;
}

bool RefinableObjClock::HasParent(const RefinableObjClock &clock) const
{
   for(std::set<RefinableObjClock*>::iterator pos=mvParent.begin();
       pos!=mvParent.end();++pos)
   {
      if((*pos)==&clock) return true;
      if((*pos)->HasParent(clock)) return true;
   }
   return false;
}
//######################################################################
//    Restraint
//######################################################################
Restraint::Restraint():
mpRefParType(gpRefParTypeObjCryst)
{}

Restraint::Restraint(const RefParType *type):
mpRefParType(type)
{}

Restraint::~Restraint()
{}

const RefParType* Restraint::GetType()const{return mpRefParType;}

void Restraint::SetType(const RefParType *type){mpRefParType=type;}

REAL Restraint::GetLogLikelihood()const{return 0.;}

//######################################################################
//    RefinablePar
//######################################################################

RefinablePar::RefinablePar():
Restraint(),
mName(""),mpValue(0),mMin(0),mMax(0),
mHasLimits(false),mIsFixed(true),mIsUsed(true),mIsPeriodic(false),
mPeriod(0.),mGlobalOptimStep(1.),mDerivStep(1e-5),mRefParDerivStepModel(REFPAR_DERIV_STEP_ABSOLUTE),
mSigma(0.),mHumanScale(1.),mHasAssignedClock(false),mpClock(0)
#ifdef __WX__CRYST__
,mpWXFieldRefPar(0)
#endif
{}

RefinablePar::RefinablePar(  const string &name,
                     REAL *refPar,
                     const REAL min,
                     const REAL max,
                     const RefParType *type,
                     RefParDerivStepModel derivMode,
                     const bool hasLimits,
                     const bool isFixed,
                     const bool isUsed,
                     const bool isPeriodic,
                     const REAL humanScale,
                     REAL period):
Restraint(type),
mName(name),mpValue(refPar),mMin(min),mMax(max),
mHasLimits(hasLimits),mIsFixed(isFixed),mIsUsed(isUsed),mIsPeriodic(isPeriodic),mPeriod(period),
mGlobalOptimStep((max-min)/100.),mDerivStep(1e-5),mRefParDerivStepModel(derivMode),
mSigma(0.),mHumanScale(humanScale),
#if 0
mUseEquation(false),mEquationNbRefPar(0),mEquationCoeff(0),
#endif
mHasAssignedClock(false),mpClock(0)
#ifdef __WX__CRYST__
,mpWXFieldRefPar(0)
#endif
{}

RefinablePar::~RefinablePar()
{
   #ifdef __WX__CRYST__
   this->WXDelete();
   #endif
}

void RefinablePar::Init(const string &name,
                        REAL *refPar,
                        const REAL min,
                        const REAL max,
                        const RefParType *type,
                        RefParDerivStepModel derivMode,
                        const bool hasLimits,
                        const bool isFixed,
                        const bool isUsed,
                        const bool isPeriodic,
                        const REAL humanScale,
                        REAL period)
{
   mName=name;
   mpValue=refPar;
   mMin=min;
   mMax=max;
   Restraint::SetType(type);
   mHasLimits=hasLimits;
   mIsFixed=isFixed;
   mIsUsed=isUsed;
   mIsPeriodic=isPeriodic;
   mPeriod=period;
   mGlobalOptimStep=(max-min)/100.;
   mDerivStep=1e-5;
   mRefParDerivStepModel=derivMode;
   mSigma=0.;
   mHumanScale=humanScale;
   #if 0
   mUseEquation=false;
   mEquationNbRefPar=0;
   mEquationCoeff=0;
   #endif
   mHasAssignedClock=false;
   mpClock=0;
}
RefinablePar::RefinablePar(const RefinablePar &old):
Restraint(old)
{
   mpValue=old.mpValue;
   #ifdef __WX__CRYST__
   mpWXFieldRefPar=0;
   #endif
   this->CopyAttributes(old);
   mHasAssignedClock=old.mHasAssignedClock;
   mpClock=old.mpClock;
}

void RefinablePar::CopyAttributes(const RefinablePar&old)
{
   mName=old.mName;
   mMin=old.GetMin();
   mMax=old.GetMax();
   mHasLimits=old.mHasLimits;
   mIsFixed=old.mIsFixed;
   mIsUsed=old.mIsUsed;
   mIsPeriodic=old.mIsPeriodic;
   mPeriod=old.mPeriod;
   mGlobalOptimStep=old.mGlobalOptimStep;
   mDerivStep=old.mDerivStep;
   mRefParDerivStepModel=old.mRefParDerivStepModel;
   mSigma=old.mSigma;
   mHumanScale=old.mHumanScale;
   #if 0
   mUseEquation=old.mUseEquation;
   mEquationNbRefPar=old.mEquationNbRefPar;
   mEquationCoeff=old.mEquationCoeff;
   #endif
}

REAL RefinablePar::GetValue()const
{
   #if 0
   if(true==mUseEquation)
   {
      VFN_DEBUG_MESSAGE("RefinablePar::Value():Evaluating Equation",0)
      REAL tmp=mEquationCoeff(0);
      for(int i=0;i<mEquationNbRefPar;i++)
         tmp += mEquationCoeff(i+1) * mEquationRefPar[i]->GetValue();
      *mpValue = tmp;
   }
   #endif
   return *mpValue;
}

const REAL* RefinablePar::GetPointer()const
{
   return mpValue;
}

void RefinablePar::SetValue(const REAL value)
{
   if(*mpValue == value) return;
   this->Click();
   VFN_DEBUG_MESSAGE("RefinablePar::SetValue()",2)
   #if 0
   if(true==mUseEquation)
   {
      cout << "RefinablePar::SetValue(): this parameter is defined by an equation !!" <<endl;
      throw 0;
   }
   #endif
   *mpValue = value;
   /*
   if(this->IsLimited() ==true)
   {
      if(true==this->IsPeriodic())
      {
         if(*mpValue > this->GetMax()) *mpValue -= this->GetMax()-this->GetMin();
         if(*mpValue < this->GetMin()) *mpValue += this->GetMax()-this->GetMin();
      }
      else
      {
         if(*mpValue > this->GetMax()) *mpValue=this->GetMax();
         if(*mpValue < this->GetMin()) *mpValue=this->GetMin();
      }
   }
   */
   if(this->IsLimited() ==true)
   {
      if(*mpValue > this->GetMax()) *mpValue=this->GetMax();
      if(*mpValue < this->GetMin()) *mpValue=this->GetMin();
   }
   else if(true==this->IsPeriodic())
   {
      if(*mpValue > mPeriod) *mpValue -= mPeriod;
      if(*mpValue < 0) *mpValue += mPeriod;
   }
}

const REAL& RefinablePar::GetHumanValue() const
{
   static REAL val;
   val = *mpValue * mHumanScale;
   return val;
}

void RefinablePar::SetHumanValue(const REAL &value)
{
   if(*mpValue == (value/mHumanScale)) return;
   this->Click();
   VFN_DEBUG_MESSAGE("RefinablePar::SetHumanValue()",2)
   #if 0
   if(true==mUseEquation)
   {
      cout << "RefinablePar::SetValue(): this parameter is defined by an equation !!" <<endl;
      throw 0;
   }
   #endif
   *mpValue = value/mHumanScale;
   /*
   if(this->IsLimited() ==true)
   {
      if(true==this->IsPeriodic())
      {
         if(*mpValue > this->GetMax()) *mpValue -= this->GetMax()-this->GetMin();
         if(*mpValue < this->GetMin()) *mpValue += this->GetMax()-this->GetMin();
      }
      else
      {
         if(*mpValue > this->GetMax()) *mpValue=this->GetMax();
         if(*mpValue < this->GetMin()) *mpValue=this->GetMin();
      }
   }
   */
   if(this->IsLimited() ==true)
   {
      if(*mpValue > this->GetMax()) *mpValue=this->GetMax();
      if(*mpValue < this->GetMin()) *mpValue=this->GetMin();
   }
   else if(true==this->IsPeriodic())
   {
      if(*mpValue > mPeriod) *mpValue -= mPeriod;
      if(*mpValue < 0) *mpValue += mPeriod;
   }
}

void RefinablePar::Mutate(const REAL mutateValue)
{
   if(0==mutateValue) return;
   VFN_DEBUG_MESSAGE("RefinablePar::Mutate():"<<this->GetName(),1)
   if(true==mIsFixed) return;
   this->Click();
   #if 0
   if(true==mUseEquation)
   {
      cout << "RefinablePar::Mutate(): this parameter is defined by an equation !!" <<endl;
      throw 0;
   }
   #endif
   *mpValue += mutateValue;
   /*
   if(this->IsLimited() ==true)
   {
      if(true==this->IsPeriodic())
      {
         if(*mpValue > this->GetMax()) *mpValue -= this->GetMax()-this->GetMin();
         if(*mpValue < this->GetMin()) *mpValue += this->GetMax()-this->GetMin();
      }
      else
      {
         if(*mpValue > this->GetMax()) *mpValue=this->GetMax();
         if(*mpValue < this->GetMin()) *mpValue=this->GetMin();
      }
   }
   */
   if(this->IsLimited() ==true)
   {
      if(*mpValue > this->GetMax()) *mpValue=this->GetMax();
      if(*mpValue < this->GetMin()) *mpValue=this->GetMin();
   }
   else if(true==this->IsPeriodic())
   {
      //if(*mpValue > mPeriod) *mpValue -= mPeriod;
      *mpValue=fmod((REAL)*mpValue,(REAL)mPeriod);
      if(*mpValue < 0) *mpValue += mPeriod;
   }
   VFN_DEBUG_MESSAGE("RefinablePar::Mutate():End",0)
}

void RefinablePar::MutateTo(const REAL mutateValue)
{
   VFN_DEBUG_MESSAGE("RefinablePar::MutateTo()",2)
   if(true==mIsFixed) return;
   if(*mpValue == mutateValue)return;
   this->Click();
   #if 0
   if(true==mUseEquation)
   {
      cout << "RefinablePar::Mutate(): this parameter is defined by an equation !!" <<endl;
      throw 0;
   }
   #endif
   *mpValue = mutateValue;
   /*
   if(this->IsLimited() ==true)
   {
      if(true==this->IsPeriodic())
      {
         if(*mpValue > this->GetMax()) *mpValue -= this->GetMax()-this->GetMin();
         if(*mpValue < this->GetMin()) *mpValue += this->GetMax()-this->GetMin();
      }
      else
      {
         if(*mpValue > this->GetMax()) *mpValue=this->GetMax();
         if(*mpValue < this->GetMin()) *mpValue=this->GetMin();
      }
   }
   */
   if(this->IsLimited() ==true)
   {
      if(*mpValue > this->GetMax()) *mpValue=this->GetMax();
      if(*mpValue < this->GetMin()) *mpValue=this->GetMin();
   }
   else if(true==this->IsPeriodic())
   {
      if(*mpValue > mPeriod) *mpValue -= mPeriod;
      if(*mpValue < 0) *mpValue += mPeriod;
   }
}

REAL  RefinablePar::GetSigma()const {return mSigma;}
REAL  RefinablePar::GetHumanSigma()const {return mSigma*mHumanScale;}
void   RefinablePar::SetSigma(const REAL sigma) {mSigma=sigma; this->Click();}

void RefinablePar::Print() const
{
      cout << this->GetName() << " : " << this->GetHumanValue()
           << " Fixed:"<< mIsFixed <<" Periodic:"<<mIsPeriodic<<" Limited:"<<mHasLimits
           << " Min:" << this->GetHumanMin() << " Max:" << this->GetHumanMax()
           << " Step:" <<GetGlobalOptimStep()
      #ifdef __DEBUG__
           << ",HasClock=" << mHasAssignedClock << " at " << mpClock
      #endif
           <<endl;
}

string RefinablePar::GetName()const {return mName;}
void RefinablePar::SetName(const string &name) {mName=name;}

bool RefinablePar::IsFixed()const {return mIsFixed;}
void RefinablePar::SetIsFixed(const bool b)
{
   VFN_DEBUG_MESSAGE("RefinablePar::SetIsFixed():"<<this->GetName(),1)
   mIsFixed=b;
}

bool RefinablePar::IsLimited()const {return mHasLimits;}
void RefinablePar::SetIsLimited(const bool b) {mHasLimits=b;this->Click();}

bool RefinablePar::IsUsed()const {return mIsUsed;}
void RefinablePar::SetIsUsed(const bool b) {mIsUsed=b;this->Click();}

bool RefinablePar::IsPeriodic()const {return mIsPeriodic;}
void RefinablePar::SetIsPeriodic(const bool b,REAL period)
{mIsPeriodic=b;mPeriod=period;this->Click();}


REAL RefinablePar::GetMin()const   {return mMin;}
void  RefinablePar::SetMin(const REAL min) { mMin=min;this->Click();}
REAL RefinablePar::GetHumanMin()const   {return mMin * mHumanScale;}
void  RefinablePar::SetHumanMin(const REAL min) { mMin=min/mHumanScale;this->Click();}

REAL RefinablePar::GetMax()const   {return mMax;}
void  RefinablePar::SetMax(const REAL max) { mMax=max;this->Click();}
REAL RefinablePar::GetHumanMax()const   {return mMax * mHumanScale;}
void  RefinablePar::SetHumanMax(const REAL max) { mMax=max/mHumanScale;this->Click();}

REAL RefinablePar::GetPeriod()const   {return mPeriod;}
void  RefinablePar::SetPeriod(const REAL period)
{ mPeriod=period;this->Click();}

REAL  RefinablePar::GetDerivStep()const
{
   if(REFPAR_DERIV_STEP_ABSOLUTE==mRefParDerivStepModel) return mDerivStep;
   REAL d=mDerivStep* (*mpValue);

   //:KLUDGE: Parameter will probably has a singular value, so it should not matter..
   if(d == 0.) return 1e-8;
   return d;
}

void RefinablePar::SetDerivStep(const REAL step)
{
   this->Click();
   mDerivStep = step;
}

REAL RefinablePar::GetGlobalOptimStep()const {return mGlobalOptimStep;}
void  RefinablePar::SetGlobalOptimStep(const REAL step) {mGlobalOptimStep=step;}

REAL RefinablePar::GetHumanScale()const {return mHumanScale;}
void  RefinablePar::SetHumanScale(const REAL scale) {mHumanScale=scale;}
#if 0
void RefinablePar::SetUseEquation(const bool useItOrNot,const REAL c0)
{
   this->Click();
   mUseEquation=useItOrNot;
   if(true==mUseEquation)
   {
      mEquationCoeff.resize(mEquationMaxRefPar);
      mEquationCoeff(0)=c0;
   }
}

void RefinablePar::SetUseEquation(const bool useItOrNot,const REAL c0,
                                  const REAL c1, const RefinablePar &refpar1)
{
   this->Click();
   mUseEquation=useItOrNot;
   if(true==mUseEquation)
   {
      mEquationCoeff.resize(mEquationMaxRefPar);
      mEquationCoeff(0)=c0;
      mEquationCoeff(1)=c1;
      mEquationRefPar[0]=&refpar1;
   }
}

void RefinablePar::SetUseEquation(const bool useItOrNot,const REAL c0,
                                  const REAL c1, const RefinablePar &refpar1,
                                  const REAL c2, const RefinablePar &refpar2)
{
   this->Click();
   mUseEquation=useItOrNot;
   if(true==mUseEquation)
   {
      mEquationCoeff.resize(mEquationMaxRefPar);
      mEquationCoeff(0)=c0;
      mEquationCoeff(1)=c1;
      mEquationCoeff(2)=c2;
      mEquationRefPar[0]=&refpar1;
      mEquationRefPar[1]=&refpar2;
   }
}

void RefinablePar::SetUseEquation(const bool useItOrNot,const REAL c0,
                                  const REAL c1, const RefinablePar &refpar1,
                                  const REAL c2, const RefinablePar &refpar2,
                                  const REAL c3, const RefinablePar &refpar3)
{
   this->Click();
   mUseEquation=useItOrNot;
   if(true==mUseEquation)
   {
      mEquationCoeff.resize(mEquationMaxRefPar);
      mEquationCoeff(0)=c0;
      mEquationCoeff(1)=c1;
      mEquationCoeff(2)=c2;
      mEquationCoeff(3)=c2;
      mEquationRefPar[0]=&refpar1;
      mEquationRefPar[1]=&refpar2;
      mEquationRefPar[2]=&refpar2;
   }
}
#endif
void RefinablePar::AssignClock(RefinableObjClock &clock)
{
   VFN_DEBUG_MESSAGE("RefinablePar::AssignClock() for "<<this->GetName()<< "at "<<&clock,4)
   mpClock=&clock;
   mHasAssignedClock=true;
}
void RefinablePar::Click()
{
   if(false==mHasAssignedClock) return;
   VFN_DEBUG_MESSAGE("RefinablePar::Click():"<<this->GetName(),1)
   mpClock->Click();
   //mpClock->Print();
   //VFN_DEBUG_MESSAGE("RefinablePar::Click():End",2)
}

void RefinablePar::SetLimitsAbsolute(const REAL min, const REAL max)
{
   //:TODO: check limits
   mMin=min;
   mMax=max;
   mHasLimits=true;
}
void RefinablePar::SetLimitsRelative(const REAL min, const REAL max)
{
   VFN_DEBUG_MESSAGE("RefinablePar::SetLimitsRelative():"<<this->GetName(),1)
   //:TODO: check limits
   mMin=this->GetValue()+min;
   mMax=this->GetValue()+max;
   this->SetIsLimited(true);
}
void RefinablePar::SetLimitsProportional(const REAL min, const REAL max)
{
   //:TODO: check limits
   mMin=this->GetValue()*min;
   mMax=this->GetValue()*max;
   this->SetIsLimited(true);
}
#ifdef __WX__CRYST__
WXCrystObjBasic* RefinablePar::WXCreate(wxWindow *parent)
{
   VFN_DEBUG_MESSAGE("RefinablePar::WXCreate()",8)
   if(mpWXFieldRefPar!=0)
   {
      throw ObjCrystException((string)"RefinablePar::WXCreate():"+this->GetName()+(string)" WXFieldRefPar already exists !");
   }
   mpWXFieldRefPar=new WXFieldRefPar (parent,this->GetName(),this);
   return (WXCrystObjBasic*) mpWXFieldRefPar;
}
WXCrystObjBasic* RefinablePar::WXGet()
{
   return (WXCrystObjBasic*) mpWXFieldRefPar;
}
void RefinablePar::WXDelete()
{
   if(0!=mpWXFieldRefPar)
   {
      VFN_DEBUG_MESSAGE("RefinablePar::WXDelete():"<<mName,5)
      delete mpWXFieldRefPar;
   }
   mpWXFieldRefPar=0;
}
void RefinablePar::WXNotifyDelete()
{
   VFN_DEBUG_MESSAGE("RefinablePar::WXNotifyDelete():"<<mName,5)
   mpWXFieldRefPar=0;
}
#endif

//######################################################################
//    RefObjOpt
//######################################################################
RefObjOpt::RefObjOpt()
{
   #ifdef __WX__CRYST__
   mpWXFieldOption=0;
   #endif
}

RefObjOpt::~RefObjOpt()
{}

void RefObjOpt::Init(const int nbChoice,
                     const string *name,
                     const string *choiceNames)
{
   VFN_DEBUG_MESSAGE("RefObjOpt::Init()"<<*name,5)
   mNbChoice=nbChoice;
   mChoice=0;
   mpName=name;
   mpChoiceName=choiceNames;
}

int RefObjOpt::GetNbChoice()const
{ return mNbChoice;}

int RefObjOpt::GetChoice()const
{ return mChoice;}

void RefObjOpt::SetChoice(const int choice)
{
   if(mChoice==choice)return;
   VFN_DEBUG_MESSAGE("RefObjOpt::SetChoice()"<<this->GetName()<< \
                     " to "<<this->GetChoiceName(choice),5)
   mChoice=choice;
   mClock.Click();
}
void RefObjOpt::SetChoice(const string &choiceName)
{
   int choice;
   for(choice=0;choice<mNbChoice;choice++) if(choiceName==*(mpChoiceName+choice)) break;
   if(choice==mNbChoice) choice=0;
   this->SetChoice(choice);
}

const string& RefObjOpt::GetName()const
{
   return *mpName;
}

const string& RefObjOpt::GetClassName()const
{
   static string className="Option";
   return className;
}

const string& RefObjOpt::GetChoiceName(const int i)const
{
   return *(mpChoiceName+i);
}

const RefinableObjClock& RefObjOpt::GetClock()const{return mClock;}

#ifdef __WX__CRYST__
WXCrystObjBasic* RefObjOpt::WXCreate(wxWindow *parent)
{
   VFN_DEBUG_MESSAGE("RefObjOpt::WXCreate()",8)
   mpWXFieldOption=new WXFieldOption (parent,-1,this);
   return mpWXFieldOption;
}
WXCrystObjBasic* RefObjOpt::WXGet()
{
   return mpWXFieldOption;
}
void RefObjOpt::WXDelete()
{
   if(0!=mpWXFieldOption)
   {
      VFN_DEBUG_MESSAGE("RefObjOpt::WXDelete()",5)
      delete mpWXFieldOption;
   }
   mpWXFieldOption=0;
}
void RefObjOpt::WXNotifyDelete()
{
   VFN_DEBUG_MESSAGE("RefObjOpt::WXNotifyDelete()",5)
   mpWXFieldOption=0;
}
#endif

//######################################################################
//    RefObjOption
//######################################################################
template<class T> RefObjOption<T>::RefObjOption(T* obj):
mpObj(obj)
{}

template<class T> RefObjOption<T>::~RefObjOption()
{}

template<class T> void RefObjOption<T>::SetChoice(const int choice)
{
   if(mChoice==choice)return;
   VFN_DEBUG_MESSAGE("RefObjOption<T>::SetChoice()"<<this->GetName()<< \
                     " to "<<this->GetChoiceName(choice),5)
   mChoice=choice;
   mClock.Click();
   if(mfpSetNewValue !=0) (mpObj->*mfpSetNewValue)(choice);
}

template<class T> void RefObjOption<T>::Init(const int nbChoice,
                                             const string *name,
                                             const string *choiceNames,
                                             void (T::*fp)(const int))
{
   this->RefObjOpt::Init(nbChoice,name,choiceNames);
   mfpSetNewValue=fp;
}

//######################################################################
//    ObjRegistry
//######################################################################
#ifdef __WX__CRYST__
bool operator==(const wxString&wx,const string&str)
{
   return wx==str.c_str();
}
bool operator==(const string&str,const wxString&wx)
{
   return wx==str.c_str();
}
#endif

template<class T> ObjRegistry<T>::ObjRegistry():
mName(""),mAutoUpdateUI(true)
#ifdef __WX__CRYST__
,mpWXRegistry(0)
#endif
{
   VFN_DEBUG_MESSAGE("ObjRegistry::ObjRegistry()",5)
}

template<class T> ObjRegistry<T>::ObjRegistry(const string &name):
mName(name),mAutoUpdateUI(true)
#ifdef __WX__CRYST__
,mpWXRegistry(0)
#endif
{
   VFN_DEBUG_MESSAGE("ObjRegistry::ObjRegistry(name):"<<mName,5)
}

//:TODO: a copy constructor
template<class T> ObjRegistry<T>::~ObjRegistry()
{
   VFN_DEBUG_MESSAGE("ObjRegistry::~ObjRegistry():"<<mName,5)
   #ifdef __WX__CRYST__
   this->WXDelete();
   #endif
}

template<class T> void ObjRegistry<T>::Register(T &obj)
{
   VFN_DEBUG_ENTRY("ObjRegistry("<<mName<<")::Register():"<<obj.GetName(),2)
   typename vector<T*>::iterator pos=find(mvpRegistry.begin(),mvpRegistry.end(),&obj);
   if(pos!=mvpRegistry.end())
   {
      VFN_DEBUG_EXIT("ObjRegistry("<<mName<<")::Register():"<<obj.GetName()<<"Already registered!",2)
      return;
   }
   mvpRegistry.push_back(&obj);
   mListClock.Click();
   #ifdef __WX__CRYST__
   if((0!=mpWXRegistry) && mAutoUpdateUI)
      mpWXRegistry->Add(obj.WXCreate(mpWXRegistry));
   #endif
   //this->Print();
   VFN_DEBUG_EXIT("ObjRegistry("<<mName<<")::Register():"<<obj.GetName(),2)
}

template<class T> void ObjRegistry<T>::DeRegister(T &obj)
{
   VFN_DEBUG_ENTRY("ObjRegistry("<<mName<<")::Deregister(&obj)"<<mvpRegistry.size(),2)
   if (mvpRegistry.size() == 0)
   {// This may happen if an object is deleted several times due to inherited destructors
    // :TODO: make sure it does not happen, while making sure WXGet() below
    //is not pure virtual because the child destructor has already done its job...
      VFN_DEBUG_EXIT("ObjRegistry(" << mName << ")::Deregister(&obj): EMPTY registry", 2)
      return;
   }
   //this->Print();
   typename vector<T*>::iterator pos=find(mvpRegistry.begin(),mvpRegistry.end(),&obj);
   if(pos==mvpRegistry.end())
   {
      VFN_DEBUG_EXIT("ObjRegistry("<<mName<<")::Deregister(&obj):NOT FOUND !!!",2)
      return; //:TODO: throw something ?
   }
   #ifdef __WX__CRYST__
   if(0!=mpWXRegistry) mpWXRegistry->Remove(obj.WXGet());
   #endif
   mvpRegistry.erase(pos);
   mListClock.Click();
   VFN_DEBUG_EXIT("ObjRegistry("<<mName<<")::Deregister(&obj)",2)
}

template<class T> void ObjRegistry<T>::DeRegister(const string &objName)
{
   VFN_DEBUG_ENTRY("ObjRegistry("<<mName<<")::Deregister(name):"<<objName,2)

   const long i=this->Find(objName);
   if(-1==i)
   {
      VFN_DEBUG_EXIT("ObjRegistry("<<mName<<")::Deregister(name): NOT FOUND !!!",2)
      return; //:TODO: throw something ?
   }
   //:KLUDGE: should directly do an iterator search on the name...
   typename vector<T*>::iterator pos=find(mvpRegistry.begin(),mvpRegistry.end(),mvpRegistry[i]);

   #ifdef __WX__CRYST__
   if(0!=mpWXRegistry) mpWXRegistry->Remove((*pos)->WXGet());
   #endif
   mvpRegistry.erase(pos);
   mListClock.Click();
   VFN_DEBUG_EXIT("ObjRegistry("<<mName<<")::Deregister(name):",2)
}

template<class T> void ObjRegistry<T>::DeRegisterAll()
{
   VFN_DEBUG_ENTRY("ObjRegistry("<<mName<<")::DeRegisterAll():",5)
   #ifdef __WX__CRYST__
   if(0!=mpWXRegistry)
   {
      typename vector<T*>::iterator pos;
      for(pos=mvpRegistry.begin();pos!=mvpRegistry.end();++pos)
         mpWXRegistry->Remove((*pos)->WXGet());
   }
   #endif
   mvpRegistry.clear();
   mListClock.Click();
   VFN_DEBUG_EXIT("ObjRegistry("<<mName<<")::DeRegisterAll():",5)
}

template<class T> void ObjRegistry<T>::DeleteAll()
{
   VFN_DEBUG_ENTRY("ObjRegistry("<<mName<<")::DeleteAll():",5)
   vector<T*> reg=mvpRegistry;//mvpRegistry will be modified as objects are deleted, so use a copy
   typename vector<T*>::iterator pos;
   for(pos=reg.begin();pos!=reg.end();++pos) delete *pos;
   mvpRegistry.clear();
   mListClock.Click();
   VFN_DEBUG_EXIT("ObjRegistry("<<mName<<")::DeleteAll():",5)
}

template<class T> T& ObjRegistry<T>::GetObj(const unsigned int i)
{
   if(i>=this->GetNb()) throw ObjCrystException("ObjRegistry<T>::GetObj(i): i >= nb!");
   return *(mvpRegistry[i]);
}

template<class T> const T& ObjRegistry<T>::GetObj(const unsigned int i) const
{
   if(i>=this->GetNb()) throw ObjCrystException("ObjRegistry<T>::GetObj(i): i >= nb!");
   return *(mvpRegistry[i]);
}

template<class T> T& ObjRegistry<T>::GetObj(const string &objName)
{
   const long i=this->Find(objName);
   return *(mvpRegistry[i]);
}

template<class T> const T& ObjRegistry<T>::GetObj(const string &objName) const
{
   const long i=this->Find(objName);
   return *(mvpRegistry[i]);
}

template<class T> T& ObjRegistry<T>::GetObj(const string &objName,
                                                  const string& className)
{
   const long i=this->Find(objName,className);
   return *(mvpRegistry[i]);
}

template<class T> const T& ObjRegistry<T>::GetObj(const string &objName,
                                                        const string& className) const
{
   const long i=this->Find(objName,className);
   return *(mvpRegistry[i]);
}

template<class T> long ObjRegistry<T>::GetNb()const{return (long)mvpRegistry.size();}

template<class T> void ObjRegistry<T>::Print()const
{
   VFN_DEBUG_MESSAGE("ObjRegistry::Print():",2)
   cout <<mName<<" :"<<this->GetNb()<<" object registered:" <<endl;

   for(long i=0;i<this->GetNb();++i)
      cout << boost::format("#%3d:%s(%s)") %i %this->GetObj(i).GetClassName() %this->GetObj(i).GetName()<<endl;
}

template<class T> void ObjRegistry<T>::SetName(const string &name){ mName=name;}

template<class T> const string& ObjRegistry<T>::GetName()const { return mName;}

template<class T> long ObjRegistry<T>::Find(const string &objName) const
{
   VFN_DEBUG_MESSAGE("ObjRegistry::Find(objName)",2)
   long index=-1;
   //bool error=false;
   for(long i=this->GetNb()-1;i>=0;i--)
      if( mvpRegistry[i]->GetName() == objName) return i;
   //      if(-1 != index) error=true ;else index=i;
   //if(true == error)
   //{
   //   cout << "ObjRegistry::Find(name) : ";
   //   cout << "found duplicate name ! This *cannot* be !!" ;
   //   cout << objName <<endl;
   //   this->Print();
   //   throw 0;
   //}
   cout << "ObjRegistry<T>::Find("<<objName<<"): Not found !!"<<endl;
   this->Print();
   throw ObjCrystException("ObjRegistry<T>::Find("+objName+"): Not found !!");
   return index;
}

template<class T> long ObjRegistry<T>::Find(const string &objName,
                                            const string &className,
                                             const bool nothrow) const
{
   VFN_DEBUG_MESSAGE("ObjRegistry::Find(objName,className)",2)
   long index=-1;
   //bool error=false;
   for(long i=this->GetNb()-1;i>=0;i--)
      if( mvpRegistry[i]->GetName() == objName)
         if(className==mvpRegistry[i]->GetClassName()) return i;
   //      if(-1 != index) error=true ;else index=i;
   //if(true == error)
   //{
   //   cout << "ObjRegistry::Find(name) : ";
   //   cout << "found duplicate name ! This *cannot* be !!" ;
   //   cout << objName <<endl;
   //   this->Print();
   //   throw 0;
   //}
   cout << "ObjRegistry<T>::Find("<<objName<<","<<className<<"): Not found !!"<<endl;
   this->Print();
   if(nothrow==false)
      throw ObjCrystException("ObjRegistry<T>::Find("+objName+","+className+"): Not found !!");
   return index;
}

template<class T> long ObjRegistry<T>::Find(const T &obj) const
{
   VFN_DEBUG_MESSAGE("ObjRegistry::Find(&obj)",2)
   for(long i=this->GetNb()-1;i>=0;i--)
      if( mvpRegistry[i]== &obj)  return i;
   //:TODO: throw something
   return -1;
}

template<class T> long ObjRegistry<T>::Find(const T *pobj) const
{
   VFN_DEBUG_MESSAGE("ObjRegistry::Find(&obj)",2)
   for(long i=this->GetNb()-1;i>=0;i--)
      if( mvpRegistry[i]== pobj)  return i;
   //:TODO: throw something
   return -1;
}

template<class T> const RefinableObjClock& ObjRegistry<T>::GetRegistryClock()const{return mListClock;}

template<class T> void ObjRegistry<T>::AutoUpdateUI(const bool autoup)
{
   mAutoUpdateUI=autoup;
}

template<class T> void ObjRegistry<T>::UpdateUI()
{
   #ifdef __WX__CRYST__
   for(unsigned int i=0;i<this->GetNb();i++)
   {
      if((this->GetObj(i).WXGet()==NULL) && (0!=mpWXRegistry))
            mpWXRegistry->Add(this->GetObj(i).WXCreate(mpWXRegistry));
   }
   #endif
}

template<class T> std::size_t ObjRegistry<T>::size() const
{
   return (std::size_t) mvpRegistry.size();
}

template<class T> typename vector<T*>::const_iterator ObjRegistry<T>::begin() const
{
   return mvpRegistry.begin();
}

template<class T> typename vector<T*>::const_iterator ObjRegistry<T>::end() const
{
   return mvpRegistry.end();
}

#ifdef __WX__CRYST__
template<class T> WXRegistry<T>* ObjRegistry<T>::WXCreate(wxWindow *parent)
{
   VFN_DEBUG_MESSAGE("ObjRegistry<T>::WXCreate()",2)
   mpWXRegistry=new WXRegistry<T> (parent,this);
   for(int i=0;i<this->GetNb();i++)
      mpWXRegistry->Add(this->GetObj(i).WXCreate(mpWXRegistry));
   return mpWXRegistry;
}
template<class T> void ObjRegistry<T>::WXDelete()
{
   if(0!=mpWXRegistry)
   {
      VFN_DEBUG_MESSAGE("ObjRegistry<T>::WXDelete()",2)
      delete mpWXRegistry;
   }
   mpWXRegistry=0;
}
template<class T> void ObjRegistry<T>::WXNotifyDelete()
{
   VFN_DEBUG_MESSAGE("ObjRegistry<T>::WXNotifyDelete()",2)
   mpWXRegistry=0;
}
#endif

//######################################################################
//    function RefObjRegisterRecursive
//######################################################################
template<class T> void RefObjRegisterRecursive(T &obj,ObjRegistry<T> &reg)
{
   VFN_DEBUG_MESSAGE("RefObjRegisterRecursive()",3)
   reg.Register(obj);
   ObjRegistry<T> *pObjReg=&(obj.GetSubObjRegistry());
   for(int i=0;i<pObjReg->GetNb();i++)
      RefObjRegisterRecursive(pObjReg->GetObj(i),reg);
   return;
}
//######################################################################
//    function RefObjRegisterRecursive
//######################################################################

void GetSubRefObjListClockRecursive(ObjRegistry<RefinableObj> &reg,RefinableObjClock &clock)
{
   if(reg.GetRegistryClock()>clock) clock=reg.GetRegistryClock();
   for(int i=0;i<reg.GetNb();i++)
      GetSubRefObjListClockRecursive(reg.GetObj(i).GetSubObjRegistry(),clock);
}

//######################################################################
//    RefinableObj
//######################################################################

ObjRegistry<RefinableObj> gRefinableObjRegistry("Global RefinableObj registry");
ObjRegistry<RefinableObj> gTopRefinableObjRegistry("Global Top RefinableObj registry");

RefinableObj::RefinableObj():
mName(""),
mNbRefParNotFixed(-1),mOptimizationDepth(0),mDeleteRefParInDestructor(true)
#ifdef __WX__CRYST__
,mpWXCrystObj(0)
#endif
{
   VFN_DEBUG_MESSAGE("RefinableObj::RefinableObj()",3)
   gRefinableObjRegistry.Register(*this);
   mSubObjRegistry.SetName("Registry for sub-objects");
   mClientObjRegistry.SetName("Registry for Clients");

   VFN_DEBUG_MESSAGE("RefinableObj::RefinableObj():End",2)
}
RefinableObj::RefinableObj(const bool internalUseOnly):
mName(""),
mNbRefParNotFixed(-1),mOptimizationDepth(0),mDeleteRefParInDestructor(true)
#ifdef __WX__CRYST__
,mpWXCrystObj(0)
#endif
{
   VFN_DEBUG_MESSAGE("RefinableObj::RefinableObj(bool)",3)
   if(false==internalUseOnly) gRefinableObjRegistry.Register(*this);
   mSubObjRegistry.SetName("Registry for sub-objects");
   mClientObjRegistry.SetName("Registry for Clients");

   VFN_DEBUG_MESSAGE("RefinableObj::RefinableObj(bool):End",2)
}

RefinableObj::RefinableObj(const RefinableObj &old) {}
/*
RefinableObj::RefinableObj(const RefinableObj &old):
mName(old.mName),mMaxNbRefPar(old.mMaxNbRefPar),mSavedValuesSetIsUsed(mMaxNbSavedSets),
mOptimizationDepth(0),mDeleteRefParInDestructor(true)
#ifdef __WX__CRYST__
,mpWXCrystObj(0)
#endif
{
   VFN_DEBUG_MESSAGE("RefinableObj::RefinableObj(RefinableObj&)",3)
   mpRefPar = new RefinablePar*[mMaxNbRefPar];
   mpSavedValuesSet = new CrystVector_REAL* [mMaxNbSavedSets];
   mpSavedValuesSetName = new string* [mMaxNbSavedSets];
   mSavedValuesSetIsUsed=false;
   *this=old;
   mSubObjRegistry.SetName("Registry for sub-objects of "+mName);
   mClientObjRegistry.SetName("Registry for Clients of "+mName);
   gRefinableObjRegistry.Register(*this);
}
*/
RefinableObj::~RefinableObj()
{
   VFN_DEBUG_MESSAGE("RefinableObj::~RefinableObj():"<<this->GetName(),5)
   if(mvpRefPar.size()>0)
   {
      if(true==mDeleteRefParInDestructor)
      {
         vector<RefinablePar*>::iterator pos;
         for(pos=mvpRefPar.begin();pos!=mvpRefPar.end();pos++) delete *pos;
      }
   }
   gRefinableObjRegistry.DeRegister(*this);
   for(int i=0;i<mSubObjRegistry.GetNb();i++)
      mSubObjRegistry.GetObj(i).DeRegisterClient(*this);
   VFN_DEBUG_MESSAGE("RefinableObj::~RefinableObj():End",4)
   #ifdef __WX__CRYST__
   this->WXDelete();
   #endif
}

const string& RefinableObj::GetClassName() const
{
   const static string className="RefinableObj";
   return className;
}

const string& RefinableObj::GetName() const {return mName;}
void RefinableObj::SetName(const string &name)
{
   VFN_DEBUG_MESSAGE("RefinableObj::SetName()to :"<<name,6)
   mName=name;
   mSubObjRegistry.SetName("Registry for sub-objects of "+mName);
}
/*
void RefinableObj::operator=(const RefinableObj &old)
{
   VFN_DEBUG_MESSAGE("RefinableObj::operator=(RefinableObj&)",3)
   this->ResetParList();
   //this->AddPar(old);
   // Do not copy old saved sets
   //... but erase any that may be stored
   for(long i=0;i<mMaxNbSavedSets;i++)
      if(true==mSavedValuesSetIsUsed(i))
      {
         delete *(mpSavedValuesSetName+i);
         delete *(mpSavedValuesSet+i);
      }
   mSavedValuesSetIsUsed=false;
}
*/
void RefinableObj::PrepareForRefinement() const
{
   VFN_DEBUG_MESSAGE("RefinableObj::PrepareForRefinement()",5)
   mNbRefParNotFixed=0;
   mRefparNotFixedIndex.resize(this->GetNbPar());
   for(long i=0;i<this->GetNbPar();i++)
      if ( (this->GetPar(i).IsFixed() == false) && (this->GetPar(i).IsUsed() == true))
      {
         mRefparNotFixedIndex(mNbRefParNotFixed) = i;
         mNbRefParNotFixed++;
      }
   //this->Print();
   VFN_DEBUG_MESSAGE("RefinableObj::PrepareForRefinement():End",5)
}

void RefinableObj::FixAllPar()
{
   VFN_DEBUG_ENTRY("RefinableObj("<<this->GetClassName()<<":"
                                  <<this->GetName()<<")::FixAllPar()",4)
   for(long i=0;i<this->GetNbPar();i++) this->GetPar(i).SetIsFixed(true);
   for(int i=0;i<this->GetSubObjRegistry().GetNb();i++)
      this->GetSubObjRegistry().GetObj(i).FixAllPar();
   VFN_DEBUG_EXIT("RefinableObj("<<this->GetName()<<")::FixAllPar()",4)
}

void RefinableObj::UnFixAllPar()
{
   for(long i=0;i<this->GetNbPar();i++) this->GetPar(i).SetIsFixed(false);
   for(int i=0;i<this->GetSubObjRegistry().GetNb();i++)
      this->GetSubObjRegistry().GetObj(i).UnFixAllPar();
}

void RefinableObj::SetParIsFixed(const long parIndex,const bool fix)
{
   this->GetPar(parIndex).SetIsFixed(fix);
}

void RefinableObj::SetParIsFixed(const string& name,const bool fix)
{
   for(long i=this->GetNbPar()-1;i>=0;i--)
      if( this->GetPar(i).GetName() == name)
         this->GetPar(i).SetIsFixed(fix);
}

void RefinableObj::SetParIsFixed(const RefParType *type,const bool fix)
{
   for(long i=0;i<this->GetNbPar();i++)
      if( this->GetPar(i).GetType()->IsDescendantFromOrSameAs(type))
      {
         //cout << " Fixing ..." << this->GetPar(i).Name()<<endl;
         this->GetPar(i).SetIsFixed(fix);
      }
   for(int i=0;i<this->GetSubObjRegistry().GetNb();i++)
      this->GetSubObjRegistry().GetObj(i).SetParIsFixed(type,fix);
}

void RefinableObj::SetParIsUsed(const string& name,const bool use)
{
   for(long i=this->GetNbPar()-1;i>=0;i--)
      if( this->GetPar(i).GetName() == name)
         this->GetPar(i).SetIsUsed(use);
}

void RefinableObj::SetParIsUsed(const RefParType *type,const bool use)
{
   for(long i=0;i<this->GetNbPar();i++)
      if( this->GetPar(i).GetType()->IsDescendantFromOrSameAs(type))
      {
         //cout << " Now used (Waow!) : ..." << this->GetPar(i).Name()<<endl;
         this->GetPar(i).SetIsUsed(use);
      }
   for(int i=0;i<this->GetSubObjRegistry().GetNb();i++)
      this->GetSubObjRegistry().GetObj(i).SetParIsUsed(type,use);
}

long RefinableObj::GetNbPar()const { return mvpRefPar.size();}

long RefinableObj::GetNbParNotFixed()const {return mNbRefParNotFixed;}

RefinablePar& RefinableObj::GetPar(const long i)
{
   return *(mvpRefPar[i]);
}

const RefinablePar& RefinableObj::GetPar(const long i) const
{
   return *(mvpRefPar[i]);
}

RefinablePar& RefinableObj::GetPar(const string & name)
{
   const long i=this->FindPar(name);
   if(-1==i)
   {
      this->Print();
      throw ObjCrystException("RefinableObj::GetPar(): cannot find parameter: "+name+" in object:"+this->GetName());
   }
   return *(mvpRefPar[i]);
}

const RefinablePar& RefinableObj::GetPar(const string & name) const
{
   const long i=this->FindPar(name);
   if(-1==i)
   {
      this->Print();
      throw ObjCrystException("RefinableObj::GetPar(): cannot find parameter: "+name+" in object:"+this->GetName());
   }
   return *(mvpRefPar[i]);
}

RefinablePar& RefinableObj::GetPar(const REAL *p)
{
   const long i=this->FindPar(p);
   if(-1==i)
   {
      this->Print();
      throw ObjCrystException("RefinableObj::GetPar(*p): cannot find parameter in object:"+this->GetName());
   }
   return *(mvpRefPar[i]);
}

const RefinablePar& RefinableObj::GetPar(const REAL *p) const
{
   const long i=this->FindPar(p);
   if(-1==i)
   {
      this->Print();
      throw ObjCrystException("RefinableObj::GetPar(*p): cannot find parameter in object:"+this->GetName());
   }
   return *(mvpRefPar[i]);
}

RefinablePar& RefinableObj::GetParNotFixed(const long i)
{
   return *(mvpRefPar[mRefparNotFixedIndex(i)]);
}

const RefinablePar& RefinableObj::GetParNotFixed(const long i) const
{
   return *(mvpRefPar[mRefparNotFixedIndex(i)]);
}

void RefinableObj::AddPar(const RefinablePar &newRefPar)
{
   VFN_DEBUG_MESSAGE("RefinableObj::AddPar(RefPar&)",2)
   string name=newRefPar.GetName();
   long ct=0;
   if(this->FindPar(name)>=0)
      while(this->FindPar(name)!=-1)
      {// KLUDGE ? Extend name if another parameter already exists with the same name
         VFN_DEBUG_MESSAGE("RefinableObj::AddPar(): need to change name ?! -> "<<name<<","<<this->FindPar(name),10)
         name += "~";
         if(++ct==100) break;// KLUDGE, let go and hope for the best...
      }

   mvpRefPar.push_back(new RefinablePar(newRefPar));
   mvpRefPar.back()->SetName(name);
   mRefParListClock.Click();
}

void RefinableObj::AddPar(RefinablePar *newRefPar)
{
   VFN_DEBUG_MESSAGE("RefinableObj::AddPar(RefPar&)",2)
   string name=newRefPar->GetName();
   long ct=0;
   if(this->FindPar(name)>=0)
      while(this->FindPar(name)!=-1)
      {// KLUDGE ? Extend name if another parameter already exists with the same name
         VFN_DEBUG_MESSAGE("RefinableObj::AddPar(): need to change name ?! -> "<<name<<","<<this->FindPar(name),10)
         name += "~";
         if(++ct==100) break;// KLUDGE, let go and hope for the best...
      }
   mvpRefPar.push_back(newRefPar);
   mvpRefPar.back()->SetName(name);
   mRefParListClock.Click();
}

void RefinableObj::AddPar(RefinableObj &newRefParList,const bool copyParam)
{
   VFN_DEBUG_MESSAGE("RefinableObj::AddPar(RefParList&)" <<newRefParList.GetNbPar() ,2)
   RefinablePar *p;
   for(long i=0;i<newRefParList.GetNbPar();i++)
   {
      if(copyParam) p=new RefinablePar(newRefParList.GetPar(i));
      else p=&(newRefParList.GetPar(i));
      this->AddPar(p);
   }
}

vector<RefinablePar *>::iterator RefinableObj::RemovePar(RefinablePar *refPar)
{
   VFN_DEBUG_MESSAGE("RefinableObj::RemovePar(RefPar&)",2)
   vector<RefinablePar *>::iterator pos=find(mvpRefPar.begin(),mvpRefPar.end(),refPar);
   if(pos==mvpRefPar.end())
   {
      throw ObjCrystException("RefinableObj::RemovePar():"+refPar->GetName()
                              +"is not in this object:"+this->GetName());
   }
   return mvpRefPar.erase(pos);
}

void RefinableObj::Print() const
{
   VFN_DEBUG_ENTRY("RefinableObj::Print()",2)
   cout << "Refinable Object:"<<this->GetName()
        <<", with " << this->GetNbPar() << " parameters" <<endl;
   for(int i=0;i<this->GetNbPar();i++)
   {
      if(this->GetPar(i).IsUsed() == false) continue;
      cout << "#"<<i<<"#" << this->GetPar(i).GetName() << ": " ;
      cout << FormatFloat(this->GetPar(i).GetHumanValue(),18,12) << " ";
      if(true == this->GetPar(i).IsFixed()) cout << "Fixed";
      else
         if(true == this->GetPar(i).IsLimited())
         {
            cout << "Limited (" << this->GetPar(i).GetHumanMin()<<","
                 <<this->GetPar(i).GetHumanMax()<<")";
            if(true == this->GetPar(i).IsPeriodic()) cout << ",Periodic" ;
         }
      VFN_DEBUG_MESSAGE_SHORT(" (at "<<this->GetPar(i).mpValue<<")",5)
      if(true == this->GetPar(i).mHasAssignedClock)
      {
         VFN_DEBUG_MESSAGE_SHORT(" (Clock at "<<this->GetPar(i).mpClock<<")",5)
      }
      cout << endl;
   }
   VFN_DEBUG_EXIT("RefinableObj::Print()",2)
}

unsigned long RefinableObj::CreateParamSet(const string name) const
{
   VFN_DEBUG_ENTRY("RefinableObj::CreateParamSet()",3)
   unsigned long id;
   for(id=0;id<=mvpSavedValuesSet.size();id++)
      if(mvpSavedValuesSet.end()==mvpSavedValuesSet.find(id)) break;

   pair< CrystVector_REAL ,string> p;
   p.second=name;
   mvpSavedValuesSet.insert(make_pair(id,p));

   this->SaveParamSet(id);
   VFN_DEBUG_MESSAGE("RefinableObj::CreateParamSet(): new parameter set with id="<<id<<" and name:"<<name,2)
   VFN_DEBUG_EXIT("RefinableObj::CreateParamSet()",3)
   return id;
}

void RefinableObj::ClearParamSet(const unsigned long id)const
{
   VFN_DEBUG_ENTRY("RefinableObj::ClearParamSet()",2)
   mvpSavedValuesSet.erase(this->FindParamSet(id));
   VFN_DEBUG_EXIT("RefinableObj::ClearParamSet()",2)
}

void RefinableObj::SaveParamSet(const unsigned long id)const
{
   VFN_DEBUG_MESSAGE("RefinableObj::SaveRefParSet()",2)
   map<unsigned long,pair<CrystVector_REAL,string> >::iterator pos=this->FindParamSet(id);
   pos->second.first.resize(mvpRefPar.size());
   REAL *p=pos->second.first.data();
   for(long i=0;i<this->GetNbPar();i++) *p++ = this->GetPar(i).GetValue();
}

void RefinableObj::RestoreParamSet(const unsigned long id)
{
   VFN_DEBUG_MESSAGE("RefinableObj::RestoreRefParSet()",2)
   map<unsigned long,pair<CrystVector_REAL,string> >::iterator pos=this->FindParamSet(id);
   REAL *p=pos->second.first.data();
   for(long i=0;i<this->GetNbPar();i++)
   {
      //if( !this->GetPar(i).IsFixed() && this->GetPar(i).IsUsed())
      if(this->GetPar(i).IsUsed())
         this->GetPar(i).SetValue(*p);
      p++;
   }
}

const CrystVector_REAL & RefinableObj::GetParamSet(const unsigned long id)const
{
   VFN_DEBUG_MESSAGE("RefinableObj::GetParamSet() const",2)
   map<unsigned long,pair<CrystVector_REAL,string> >::const_iterator pos=this->FindParamSet(id);
   return pos->second.first;
}

CrystVector_REAL & RefinableObj::GetParamSet(const unsigned long id)
{
   VFN_DEBUG_MESSAGE("RefinableObj::GetParamSet()",2)
   map<unsigned long,pair<CrystVector_REAL,string> >::iterator pos=this->FindParamSet(id);
   return pos->second.first;
}

REAL RefinableObj::GetParamSet_ParNotFixedHumanValue(const unsigned long id,
                                                      const long par)const
{
   VFN_DEBUG_MESSAGE("RefinableObj::RefParSetNotFixedHumanValue()",0)
   map<unsigned long,pair<CrystVector_REAL,string> >::iterator pos=this->FindParamSet(id);
   return pos->second.first(mRefparNotFixedIndex(par));
}

const void RefinableObj::EraseAllParamSet()
{
   mvpSavedValuesSet.clear();
}

const string& RefinableObj::GetParamSetName(const unsigned long id)const
{
   VFN_DEBUG_MESSAGE("RefinableObj::GetParamSetName()",2)
   map<unsigned long,pair<CrystVector_REAL,string> >::const_iterator pos=this->FindParamSet(id);
   return pos->second.second;
}

void RefinableObj::SetLimitsAbsolute(const string &name,const REAL min,const REAL max)
{
   for(long i=this->GetNbPar()-1;i>=0;i--)
      if( this->GetPar(i).GetName() == name)
         this->GetPar(i).SetLimitsAbsolute(min,max);
}
void RefinableObj::SetLimitsAbsolute(const RefParType *type,
                                     const REAL min,const REAL max)
{
   for(long i=0;i<this->GetNbPar();i++)
      if(this->GetPar(i).GetType()->IsDescendantFromOrSameAs(type))
         this->GetPar(i).SetLimitsAbsolute(min,max);
   for(int i=0;i<this->GetSubObjRegistry().GetNb();i++)
      this->GetSubObjRegistry().GetObj(i).SetLimitsAbsolute(type,min,max);
}
void RefinableObj::SetLimitsRelative(const string &name, const REAL min, const REAL max)
{
   for(long i=this->GetNbPar()-1;i>=0;i--)
      if( this->GetPar(i).GetName() == name)
         this->GetPar(i).SetLimitsRelative(min,max);
}
void RefinableObj::SetLimitsRelative(const RefParType *type,
                                     const REAL min, const REAL max)
{
   VFN_DEBUG_MESSAGE("RefinableObj::SetLimitsRelative(RefParType*):"<<this->GetName(),2)
   for(long i=0;i<this->GetNbPar();i++)
   {
      VFN_DEBUG_MESSAGE("RefinableObj::SetLimitsRelative(RefParType*):par #"<<i,2)
      if(this->GetPar(i).GetType()->IsDescendantFromOrSameAs(type))
         this->GetPar(i).SetLimitsRelative(min,max);
   }
   for(int i=0;i<this->GetSubObjRegistry().GetNb();i++)
      this->GetSubObjRegistry().GetObj(i).SetLimitsRelative(type,min,max);
}
void RefinableObj::SetLimitsProportional(const string &name,const REAL min,const REAL max)
{
   for(long i=this->GetNbPar()-1;i>=0;i--)
      if( this->GetPar(i).GetName() == name)
         this->GetPar(i).SetLimitsProportional(min,max);
}
void RefinableObj::SetLimitsProportional(const RefParType *type,
                                         const REAL min, const REAL max)
{
   for(long i=0;i<this->GetNbPar();i++)
      if(this->GetPar(i).GetType()->IsDescendantFromOrSameAs(type))
         this->GetPar(i).SetLimitsProportional(min,max);
   for(int i=0;i<this->GetSubObjRegistry().GetNb();i++)
      this->GetSubObjRegistry().GetObj(i).SetLimitsProportional(type,min,max);
}
void RefinableObj::SetGlobalOptimStep(const RefParType *type, const REAL step)
{
   for(long i=0;i<this->GetNbPar();i++)
      if(this->GetPar(i).GetType()->IsDescendantFromOrSameAs(type))
         this->GetPar(i).SetGlobalOptimStep(step);
   for(int i=0;i<this->GetSubObjRegistry().GetNb();i++)
      this->GetSubObjRegistry().GetObj(i).SetGlobalOptimStep(type,step);
}

ObjRegistry<RefinableObj>& RefinableObj::GetSubObjRegistry()
{return mSubObjRegistry;}

const ObjRegistry<RefinableObj>& RefinableObj::GetSubObjRegistry()const
{return mSubObjRegistry;}

void RefinableObj::RegisterClient(RefinableObj &obj)const
{mClientObjRegistry.Register(obj);}

void RefinableObj::DeRegisterClient(RefinableObj &obj)const
{mClientObjRegistry.DeRegister(obj);}

const ObjRegistry<RefinableObj>& RefinableObj::GetClientRegistry()const{return mClientObjRegistry;}
      ObjRegistry<RefinableObj>& RefinableObj::GetClientRegistry()     {return mClientObjRegistry;}

bool RefinableObj::IsBeingRefined()const {return mOptimizationDepth>0;}

extern const long ID_WXOBJ_ENABLE; //These are defined in wxCryst/wxCryst.cpp
extern const long ID_WXOBJ_DISABLE;
void RefinableObj::BeginOptimization(const bool allowApproximations,
                                     const bool enableRestraints)
{
   mOptimizationDepth++;
   if(mOptimizationDepth>1) return;
   this->Prepare();
   for(int i=0;i<mSubObjRegistry.GetNb();i++)
      mSubObjRegistry.GetObj(i).BeginOptimization(allowApproximations);
   #ifdef __WX__CRYST__
   if(0!=mpWXCrystObj)
   {
      if(true==wxThread::IsMain()) mpWXCrystObj->Enable(false);
      else
      {
         wxUpdateUIEvent event(ID_WXOBJ_DISABLE);
         wxPostEvent(mpWXCrystObj,event);
      }
   }
   #endif
}

void RefinableObj::EndOptimization()
{
   mOptimizationDepth--;
   if(mOptimizationDepth<0) throw ObjCrystException("RefinableObj::EndOptimization(): mOptimizationDepth<0 !!");

   if(mOptimizationDepth>0) return;
   for(int i=0;i<mSubObjRegistry.GetNb();i++)
      mSubObjRegistry.GetObj(i).EndOptimization();
   #ifdef __WX__CRYST__
   if(0!=mpWXCrystObj)
   {
      if(true==wxThread::IsMain()) mpWXCrystObj->Enable(true);
      else
      {
         wxUpdateUIEvent event(ID_WXOBJ_ENABLE);
         wxPostEvent(mpWXCrystObj,event);
      }
   }
   #endif
}

void RefinableObj::SetApproximationFlag(const bool allow)
{
for(int i=0;i<mSubObjRegistry.GetNb();i++)
      mSubObjRegistry.GetObj(i).SetApproximationFlag(allow);
}

void RefinableObj::RandomizeConfiguration()
{
   VFN_DEBUG_ENTRY("RefinableObj::RandomizeConfiguration():"<<mName,5)
   this->PrepareForRefinement();
   for(int j=0;j<this->GetNbParNotFixed();j++)
   {
      if(true==this->GetParNotFixed(j).IsLimited())
      {
         const REAL min=this->GetParNotFixed(j).GetMin();
         const REAL max=this->GetParNotFixed(j).GetMax();
         this->GetParNotFixed(j).MutateTo(min+(max-min)*(rand()/(REAL)RAND_MAX) );
      }
      else
         if(true==this->GetParNotFixed(j).IsPeriodic())
         {

            this->GetParNotFixed(j).MutateTo((rand()/(REAL)RAND_MAX)
                  * this->GetParNotFixed(j).GetPeriod());
         }
   }
   for(int i=0;i<this->GetSubObjRegistry().GetNb();i++)
      this->GetSubObjRegistry().GetObj(i).RandomizeConfiguration();
   VFN_DEBUG_EXIT("RefinableObj::RandomizeConfiguration():Finished",5)
}

void RefinableObj::GlobalOptRandomMove(const REAL mutationAmplitude,
                                       const RefParType *type)
{
   if(mRandomMoveIsDone) return;
   VFN_DEBUG_ENTRY("RefinableObj::GlobalOptRandomMove()",2)
   for(int j=0;j<this->GetNbParNotFixed();j++)
   {
      if(this->GetParNotFixed(j).GetType()->IsDescendantFromOrSameAs(type))
         this->GetParNotFixed(j).Mutate( this->GetParNotFixed(j).GetGlobalOptimStep()
                     *2*(rand()/(REAL)RAND_MAX-0.5)*mutationAmplitude);
   }
   for(int i=0;i<mSubObjRegistry.GetNb();i++)
      mSubObjRegistry.GetObj(i).GlobalOptRandomMove(mutationAmplitude,type);
   mRandomMoveIsDone=true;
   VFN_DEBUG_EXIT("RefinableObj::GlobalOptRandomMove()",2)
}
void RefinableObj::BeginGlobalOptRandomMove()
{
   mRandomMoveIsDone=false;
   for(int i=0;i<mSubObjRegistry.GetNb();i++)
      mSubObjRegistry.GetObj(i).BeginGlobalOptRandomMove();
}

unsigned int RefinableObj::GetNbLSQFunction()const{return 0;}

REAL RefinableObj::GetLogLikelihood()const
{
   VFN_DEBUG_ENTRY("RefinableObj::GetLogLikelihood()",3)
   if(0==mvpRestraint.size())
   {
      VFN_DEBUG_EXIT("RefinableObj::GetLogLikelihood()=0",3)
      return 0;
   }
   VFN_DEBUG_MESSAGE("RefinableObj::GetLogLikelihood()there are restraints...",2)
   REAL loglike=0;
   vector< Restraint * >::const_iterator pos;
   for(pos=mvpRestraint.begin();pos!=mvpRestraint.end();pos++)
   {
      VFN_DEBUG_MESSAGE("RefinableObj::GetLogLikelihood()Restraint: "<<*pos,2)
      loglike+= (*pos)->GetLogLikelihood();
   }
   VFN_DEBUG_EXIT("RefinableObj::GetLogLikelihood()="<<loglike,3)
   return loglike;
}

// std::map<RefinablePar*, REAL>& RefinableObj::GetLogLikelihood_FullDeriv(std::set<RefinablePar *> &vPar)
// {
//    //TODO
//    return mLogLikelihood_FullDeriv;
// }

const CrystVector_REAL& RefinableObj::GetLSQCalc(const unsigned int) const
{
   throw ObjCrystException("Error: called RefinableObj::GetLSQCalc()");
   CrystVector_REAL *noWarning=new CrystVector_REAL;
   return *noWarning;
}

const CrystVector_REAL& RefinableObj::GetLSQObs(const unsigned int) const
{
   throw ObjCrystException("Error: called RefinableObj::GetLSQObs()");
   CrystVector_REAL *noWarning=new CrystVector_REAL;
   return *noWarning;
}

const CrystVector_REAL& RefinableObj::GetLSQWeight(const unsigned int) const
{
   throw ObjCrystException("Error: called RefinableObj::GetLSQWeight()");
   CrystVector_REAL *noWarning=new CrystVector_REAL;
   return *noWarning;
}

const CrystVector_REAL& RefinableObj::GetLSQDeriv(const unsigned int n, RefinablePar&par)
{
   // By default, use numerical derivatives
   const REAL step=par.GetDerivStep();
   par.Mutate(step);
   mLSQDeriv  =this->GetLSQCalc(n);
   par.Mutate(-2*step);
   mLSQDeriv -=this->GetLSQCalc(n);
   par.Mutate(step);
   mLSQDeriv /= step*2;
   return mLSQDeriv;
}

std::map<RefinablePar*, CrystVector_REAL> & RefinableObj::GetLSQ_FullDeriv(const unsigned int n,std::set<RefinablePar *> &vPar)
{
   mLSQ_FullDeriv[n].clear();
   mLSQ_FullDeriv[n][(RefinablePar*)0]=this->GetLSQCalc(n);
   for(std::set<RefinablePar *>::const_iterator pos=vPar.begin();pos!=vPar.end();pos++)
      mLSQ_FullDeriv[n][*pos]=this->GetLSQDeriv(n,**pos);

   return mLSQ_FullDeriv[n];
}

void RefinableObj::ResetParList()
{
   VFN_DEBUG_MESSAGE("RefinableObj::ResetParList()",3)
   if(mvpRefPar.size()>0)
   {
      if(true==mDeleteRefParInDestructor)
      {
         vector<RefinablePar*>::iterator pos;
         for(pos=mvpRefPar.begin();pos!=mvpRefPar.end();pos++) delete *pos;
      }
      mvpRefPar.clear();
   }
   mNbRefParNotFixed=-1;

   VFN_DEBUG_MESSAGE("RefinableObj::ResetParList():Deleting Saved Sets....",2)
   this->EraseAllParamSet();
   mRefParListClock.Click();
   VFN_DEBUG_MESSAGE("RefinableObj::ResetParList():End.",3)
}

ObjRegistry<RefObjOpt>& RefinableObj::GetOptionList()
{
    return mOptionRegistry;
}

unsigned int RefinableObj::GetNbOption()const
{
   return mOptionRegistry.GetNb();
}

RefObjOpt& RefinableObj::GetOption(const unsigned int i)
{
   VFN_DEBUG_MESSAGE("RefinableObj::GetOption()"<<i,3)
   //:TODO: Check
   return mOptionRegistry.GetObj(i);
}

RefObjOpt& RefinableObj::GetOption(const string & name)
{
    VFN_DEBUG_MESSAGE("RefinableObj::GetOption()"<<name,3)
    const long i=mOptionRegistry.Find(name);
    if(i<0)
    {
        this->Print();
        throw ObjCrystException("RefinableObj::GetOption(): cannot find option: "+name+" in object:"+this->GetName());
    }
    return mOptionRegistry.GetObj(i);
}

const RefObjOpt& RefinableObj::GetOption(const unsigned int i)const
{
   //:TODO: Check
   return mOptionRegistry.GetObj(i);
}

const RefObjOpt& RefinableObj::GetOption(const string & name)const
{
    VFN_DEBUG_MESSAGE("RefinableObj::GetOption()"<<name,3)
    const long i=mOptionRegistry.Find(name);
    if(i<0)
    {
        this->Print();
        throw ObjCrystException("RefinableObj::GetOption(): cannot find option: "+name+" in object:"+this->GetName());
    }
    return mOptionRegistry.GetObj(i);
}

void RefinableObj::GetGeneGroup(const RefinableObj &obj,
                                CrystVector_uint & groupIndex,
                                unsigned int &first) const
{
   VFN_DEBUG_MESSAGE("RefinableObj::GetGeneGroup()",4)
   for(long i=0;i<obj.GetNbPar();i++)
      for(long j=0;j<this->GetNbPar();j++)
         if(&(obj.GetPar(i)) == &(this->GetPar(j))) groupIndex(i)= first++;
}
void RefinableObj::SetDeleteRefParInDestructor(const bool b) {mDeleteRefParInDestructor=b;}

const RefinableObjClock& RefinableObj::GetRefParListClock()const{return mRefParListClock;}

REAL  RefinableObj::GetRestraintCost()const
{
   vector<Restraint*>::const_iterator pos;
   REAL cost(0);
   for(pos=mvpRestraint.begin();pos != mvpRestraint.end();++pos)
      cost += (*pos)->GetLogLikelihood();
   return 0;
}

void RefinableObj::AddRestraint(Restraint *pNewRestraint)
{
   VFN_DEBUG_MESSAGE("RefinableObj::AddRestraint(Restraint*)",2)
   mvpRestraint.push_back(pNewRestraint);
}

vector<Restraint*>::iterator RefinableObj::RemoveRestraint(Restraint *pRestraint)
{
   VFN_DEBUG_MESSAGE("RefinableObj::RemoveRestraint(Restraint*)",2)
   vector<Restraint*>::iterator pos=find(mvpRestraint.begin(),mvpRestraint.end(),pRestraint);
   if(mvpRestraint.end() != pos)
   {
      return mvpRestraint.erase(pos);
   }
   else
   {
      cout <<"RefinableObj::RemoveRestraint(..)"
           <<" Whoops... tried to remove a Restraint which does not exist..."<<endl;
   }
   return pos;
}

void RefinableObj::TagNewBestConfig()const
{
}
const RefinableObjClock& RefinableObj::GetClockMaster()const{return mClockMaster;}


void RefinableObj::UpdateDisplay()const
{
   #ifdef __WX__CRYST__
   VFN_DEBUG_ENTRY("RefinableObj::UpdateDisplay()",3)
   if(0!=mpWXCrystObj) mpWXCrystObj->CrystUpdate(true,true);
   VFN_DEBUG_EXIT("RefinableObj::UpdateDisplay()",3)
   #endif
}

long RefinableObj::FindPar(const string &name) const
{
   long index=-1;
   bool warning=false;
   for(long i=this->GetNbPar()-1;i>=0;i--)
      if( this->GetPar(i).GetName() == name)
      {
         if(-1 != index) warning=true ;else index=i;
      }
   if(true == warning)
   {
      throw ObjCrystException("RefinableObj::FindPar("+name+"): found duplicate refinable variable name in object:"+this->GetName());
   }
   return index;
}

long RefinableObj::FindPar(const REAL *p) const
{
   long index=-1;
   bool warning=false;
   for(long i=this->GetNbPar()-1;i>=0;i--)
      if( p == mvpRefPar[i]->mpValue ) //&(this->GetPar(i).GetValue())
      {
         if(-1 != index) warning=true ;else index=i;
      }
   if(true == warning)
   {
      throw ObjCrystException("RefinableObj::FindPar(*p): Found duplicate parameter in object:"+this->GetName());
   }

   return index;
}

void RefinableObj::AddSubRefObj(RefinableObj &obj)
{
   VFN_DEBUG_MESSAGE("RefinableObj::AddSubRefObj()",3)
   mSubObjRegistry.Register(obj);
   mClockMaster.AddChild(obj.GetClockMaster());
}

void RefinableObj::RemoveSubRefObj(RefinableObj &obj)
{
   VFN_DEBUG_MESSAGE("RefinableObj::RemoveSubRefObj()",3)
   mSubObjRegistry.DeRegister(obj);
   mClockMaster.RemoveChild(obj.GetClockMaster());
}

void RefinableObj::AddOption(RefObjOpt *opt)
{
   VFN_DEBUG_MESSAGE("RefinableObj::AddOption()",5)
   //:TODO: automagically resize the option array if necessary
   mOptionRegistry.Register(*opt);
   mClockMaster.AddChild(opt->GetClock());
   VFN_DEBUG_MESSAGE("RefinableObj::AddOption():End",5)
}

void RefinableObj::Prepare()
{
   VFN_DEBUG_MESSAGE("RefinableObj::Prepare()",5)
   for(int i=0;i<this->GetSubObjRegistry().GetNb();i++)
      this->GetSubObjRegistry().GetObj(i).Prepare();
}

map<unsigned long,pair<CrystVector_REAL,string> >::iterator
   RefinableObj:: FindParamSet(const unsigned long id)const
{
   VFN_DEBUG_ENTRY("RefinableObj::FindParamSet()",2)
   map<unsigned long,pair<CrystVector_REAL,string> >::iterator pos;
   pos=mvpSavedValuesSet.find(id);
   if(mvpSavedValuesSet.end() == pos)
   {//throw up
      throw ObjCrystException("RefinableObj::FindParamSet(long): Unknown saved set ! In object:"+this->GetName());
   }
   VFN_DEBUG_EXIT("RefinableObj::FindParamSet()",2)
   return pos;
}

#ifdef __WX__CRYST__
WXCrystObjBasic* RefinableObj::WXCreate(wxWindow *parent)
{
   VFN_DEBUG_MESSAGE("RefinableObj::WXCreate()",8)
   mpWXCrystObj=new WXRefinableObj (parent,this);
   return mpWXCrystObj;
}
WXCrystObjBasic* RefinableObj::WXGet()
{
   return mpWXCrystObj;
}
void RefinableObj::WXDelete()
{
   if(0!=mpWXCrystObj)
   {
      VFN_DEBUG_MESSAGE("RefinableObj::WXDelete()",5)
      delete mpWXCrystObj;
   }
   mpWXCrystObj=0;
}
void RefinableObj::WXNotifyDelete()
{
   VFN_DEBUG_MESSAGE("RefinableObj::WXNotifyDelete():"<<mName,5)
   mpWXCrystObj=0;
}
#endif
//######################################################################
//    function GetRefParListClockRecursive
//######################################################################

void GetRefParListClockRecursive(ObjRegistry<RefinableObj> &reg,RefinableObjClock &clock)
{
   for(int i=0;i<reg.GetNb();i++)
   {
      if(reg.GetObj(i).GetRefParListClock()>clock)
         clock=reg.GetObj(i).GetRefParListClock();
      GetRefParListClockRecursive(reg.GetObj(i).GetSubObjRegistry(),clock);
   }
}

//***********EXPLICIT INSTANTIATION*******************//
template void RefObjRegisterRecursive(RefinableObj &obj,ObjRegistry<RefinableObj> &reg);
}//namespace
#include "ObjCryst/ObjCryst/Crystal.h"
#include "ObjCryst/ObjCryst/Scatterer.h"
#include "ObjCryst/ObjCryst/ScatteringPower.h"
#include "ObjCryst/ObjCryst/ZScatterer.h"
#include "ObjCryst/ObjCryst/PowderPattern.h"
#include "ObjCryst/ObjCryst/DiffractionDataSingleCrystal.h"
#include "ObjCryst/ObjCryst/ScatteringCorr.h"
#include "ObjCryst/RefinableObj/GlobalOptimObj.h"
#include "ObjCryst/RefinableObj/IO.h"
#include "ObjCryst/ObjCryst/ReflectionProfile.h"

namespace ObjCryst {

template class ObjRegistry<RefObjOpt>;
template class ObjRegistry<RefinableObj>;
template class ObjRegistry<Crystal>;
template class ObjRegistry<Scatterer>;
template class ObjRegistry<ScatteringPower>;
template class ObjRegistry<ScatteringPowerAtom>;
template class ObjRegistry<PowderPattern>;
template class ObjRegistry<PowderPatternComponent>;
template class ObjRegistry<DiffractionDataSingleCrystal>;
template class ObjRegistry<OptimizationObj>;
//template class ObjRegistry<IOCrystTag>;//to be removed
template class ObjRegistry<XMLCrystTag>;
template class ObjRegistry<ZAtom>;
template class ObjRegistry<TexturePhaseMarchDollase>;
template class ObjRegistry<TextureEllipsoid>;
template class ObjRegistry<ReflectionProfile>;

template class RefObjOption<RefinableObj>;
template class RefObjOption<Crystal>;
template class RefObjOption<Radiation>;
template class RefObjOption<Scatterer>;
template class RefObjOption<ScatteringPower>;
template class RefObjOption<ScatteringPowerAtom>;
template class RefObjOption<PowderPattern>;
template class RefObjOption<PowderPatternComponent>;
template class RefObjOption<PowderPatternBackground>;
template class RefObjOption<DiffractionDataSingleCrystal>;
template class RefObjOption<PowderPatternDiffraction>;
//template class RefObjOption<GlobalOptimObj>;
//template class RefObjOption<IOCrystTag>;

}   // namespace ObjCryst
