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
*  source file for LibCryst++ PowderPattern class
*
*/

#include <cstdlib>

#include <typeinfo>
#include <stdio.h> //for sprintf()

#include "cctbx/sgtbx/space_group.h" // For fullprof export

#include "ObjCryst/ObjCryst/PowderPattern.h"
#include "ObjCryst/ObjCryst/Molecule.h" // For fullprof export
#include "ObjCryst/ObjCryst/PowderPatternBackgroundBayesianMinimiser.h"
#include "ObjCryst/RefinableObj/Simplex.h"
#include "ObjCryst/Quirks/VFNDebug.h"
#include "ObjCryst/Quirks/VFNStreamFormat.h"
#include "ObjCryst/ObjCryst/CIF.h"
#ifdef __WX__CRYST__
   #include "ObjCryst/wxCryst/wxPowderPattern.h"
#endif

#include <fstream>
#include <iomanip>
#include <sstream>

#ifdef _MSC_VER // MS VC++ predefined macros....
#undef min
#undef max
#endif

//#define USE_BACKGROUND_MAXLIKE_ERROR

namespace ObjCryst
{
////////////////////////////////////////////////////////////////////////
//
//    PowderPatternComponent
//
////////////////////////////////////////////////////////////////////////
ObjRegistry<PowderPatternComponent>
   gPowderPatternComponentRegistry("List of all PowderPattern Components");

PowderPatternComponent::PowderPatternComponent():
mIsScalable(false),mpParentPowderPattern(0)
{
   gPowderPatternComponentRegistry.Register(*this);
   mClockMaster.AddChild(mClockBraggLimits);
}

PowderPatternComponent::PowderPatternComponent(const PowderPatternComponent &old):
mIsScalable(old.mIsScalable),
mpParentPowderPattern(old.mpParentPowderPattern)
{
   mClockMaster.AddChild(mClockBraggLimits);
   if(mpParentPowderPattern!=0)
   {
      mClockMaster.AddChild(mpParentPowderPattern->GetClockPowderPatternPar());
      mClockMaster.AddChild(mpParentPowderPattern->GetClockPowderPatternXCorr());
      mClockMaster.AddChild(mpParentPowderPattern->GetClockPowderPatternRadiation());
   }
}

PowderPatternComponent::~PowderPatternComponent()
{
   gPowderPatternComponentRegistry.DeRegister(*this);
}

const string& PowderPatternComponent::GetClassName() const
{
   const static string className="PowderPatternComponent";
   return className;
}

const PowderPattern& PowderPatternComponent::GetParentPowderPattern()const
{
   return *mpParentPowderPattern;
}

std::map<RefinablePar*,CrystVector_REAL>& PowderPatternComponent::GetPowderPattern_FullDeriv(std::set<RefinablePar *> &vPar)
{
   this->CalcPowderPattern_FullDeriv(vPar);
   return mPowderPattern_FullDeriv;
}

std::map<RefinablePar*,CrystVector_REAL>& PowderPatternComponent::GetPowderPatternIntegrated_FullDeriv(std::set<RefinablePar *> &vPar)
{
   this->CalcPowderPatternIntegrated_FullDeriv(vPar);
   return mPowderPatternIntegrated_FullDeriv;
}

PowderPattern& PowderPatternComponent::GetParentPowderPattern()
{
   return *mpParentPowderPattern;
}

bool PowderPatternComponent::IsScalable()const {return mIsScalable;}

const RefinableObjClock& PowderPatternComponent::GetClockPowderPatternCalc()const
{
   return mClockPowderPatternCalc;
}
const RefinableObjClock& PowderPatternComponent::GetClockBraggLimits()const
{
   return mClockBraggLimits;
}

const list<pair<const REAL, const string> >& PowderPatternComponent::GetPatternLabelList()const
{
   return mvLabel;
}

void PowderPatternComponent::CalcPowderPattern_FullDeriv(std::set<RefinablePar *> &vPar)
{
   TAU_PROFILE("PowderPatternComponent::CalcPowderPattern_FullDeriv()","void ()",TAU_DEFAULT);
   mPowderPattern_FullDeriv.clear();
   for(std::set<RefinablePar *>::iterator par=vPar.begin();par!=vPar.end();par++)
   {
      if(*par==0)
      {
         mPowderPattern_FullDeriv[*par]=this->GetPowderPatternCalc();
         continue;
      }
      const REAL step=(*par)->GetDerivStep();
      (*par)->Mutate(step);
      mPowderPattern_FullDeriv[*par] =this->GetPowderPatternCalc();
      (*par)->Mutate(-2*step);
      mPowderPattern_FullDeriv[*par]-=this->GetPowderPatternCalc();
      (*par)->Mutate(step);
      mPowderPattern_FullDeriv[*par]/= step*2;
      if(MaxAbs(mPowderPattern_FullDeriv[*par])==0)
         mPowderPattern_FullDeriv[*par].resize(0);
   }
}

void PowderPatternComponent::CalcPowderPatternIntegrated_FullDeriv(std::set<RefinablePar *> &vPar)
{
   TAU_PROFILE("PowderPatternComponent::CalcPowderPatternIntegrated_FullDeriv()","void ()",TAU_DEFAULT);
   mPowderPatternIntegrated_FullDeriv.clear();
   for(std::set<RefinablePar *>::iterator par=vPar.begin();par!=vPar.end();par++)
   {
      if(*par==0)
      {
         mPowderPatternIntegrated_FullDeriv[*par]=*(this->GetPowderPatternIntegratedCalc().first);
         continue;
      }
      const REAL step=(*par)->GetDerivStep();
      (*par)->Mutate(step);
      mPowderPatternIntegrated_FullDeriv[*par] =*(this->GetPowderPatternIntegratedCalc().first);
      (*par)->Mutate(-2*step);
      mPowderPatternIntegrated_FullDeriv[*par]-=*(this->GetPowderPatternIntegratedCalc().first);
      (*par)->Mutate(step);
      mPowderPatternIntegrated_FullDeriv[*par]/= step*2;
      if(MaxAbs(mPowderPatternIntegrated_FullDeriv[*par])==0)
         mPowderPatternIntegrated_FullDeriv[*par].resize(0);
   }
}

////////////////////////////////////////////////////////////////////////
//
//        PowderPatternBackground
//
////////////////////////////////////////////////////////////////////////
PowderPatternBackground::PowderPatternBackground():
mBackgroundNbPoint(0),
mMaxSinThetaOvLambda(10),mModelVariance(0)
{
   mClockMaster.AddChild(mClockBackgroundPoint);
   this->InitOptions();
}

PowderPatternBackground::PowderPatternBackground(const  PowderPatternBackground &old):
mBackgroundNbPoint(old.mBackgroundNbPoint),
mBackgroundInterpPointX(old.mBackgroundInterpPointX),
mBackgroundInterpPointIntensity(old.mBackgroundInterpPointIntensity),
mMaxSinThetaOvLambda(10),mModelVariance(0)
{
   mClockMaster.AddChild(mClockBackgroundPoint);
   this->InitOptions();
   mInterpolationModel.SetChoice(old.mInterpolationModel.GetChoice());
}

PowderPatternBackground::~PowderPatternBackground(){}
const string& PowderPatternBackground::GetClassName() const
{
   const static string className="PowderPatternBackground";
   return className;
}

void PowderPatternBackground::SetParentPowderPattern(PowderPattern &s)
{
   if(mpParentPowderPattern!=0)
      mClockMaster.RemoveChild(mpParentPowderPattern->GetIntegratedProfileLimitsClock());
   mpParentPowderPattern = &s;
   mClockMaster.AddChild(mpParentPowderPattern->GetIntegratedProfileLimitsClock());
   mClockMaster.AddChild(mpParentPowderPattern->GetClockPowderPatternPar());
   mClockMaster.AddChild(mpParentPowderPattern->GetClockPowderPatternXCorr());
   mClockMaster.AddChild(mpParentPowderPattern->GetClockPowderPatternRadiation());
}
const CrystVector_REAL& PowderPatternBackground::GetPowderPatternCalc()const
{
   this->CalcPowderPattern();
   return mPowderPatternCalc;
}

pair<const CrystVector_REAL*,const RefinableObjClock*>
   PowderPatternBackground::GetPowderPatternIntegratedCalc()const
{
   VFN_DEBUG_MESSAGE("PowderPatternBackground::GetPowderPatternIntegratedCalc()",3)
   this->CalcPowderPatternIntegrated();
   return make_pair(&mPowderPatternIntegratedCalc,&mClockPowderPatternIntegratedCalc);
}

void PowderPatternBackground::ImportUserBackground(const string &filename)
{
   VFN_DEBUG_MESSAGE("PowderPatternBackground::ImportUserBackground():"<<filename,5)
   ifstream fin (filename.c_str());
   if(!fin)
   {
      throw ObjCrystException("PowderPatternBackground::ImportUserBackground() : \
Error opening file for input:"+filename);
   }
   long max=100;
   CrystVector_REAL bckgd2Theta(max),bckgd(max);

   long nbPoints=0;
   do
   {
      fin >> bckgd2Theta(nbPoints);
      fin >> bckgd(nbPoints);
      if(!fin) break;
      VFN_DEBUG_MESSAGE("Background=" << bckgd(nbPoints)\
         <<" at 2theta="<<bckgd2Theta(nbPoints),3)
      nbPoints++;
      if(nbPoints==max)
      {
         max += 100;
         bckgd2Theta.resizeAndPreserve(max);
         bckgd.resizeAndPreserve(max);
      }
   } while(fin.eof() == false);
   fin.close();
   bckgd2Theta.resizeAndPreserve(nbPoints);
   bckgd.resizeAndPreserve(nbPoints);
   if(mpParentPowderPattern!=0)
   {   if((this->GetParentPowderPattern().GetRadiation().GetWavelengthType()==WAVELENGTH_MONOCHROMATIC)
         ||(this->GetParentPowderPattern().GetRadiation().GetWavelengthType()==WAVELENGTH_ALPHA12))
      bckgd2Theta*= DEG2RAD;
   } else bckgd2Theta*= DEG2RAD;

   this->SetInterpPoints(bckgd2Theta,bckgd);
   this->InitRefParList();
   mClockBackgroundPoint.Click();
   {
      char buf [200];
      sprintf(buf,"Imported %d background points",(int)nbPoints);
      (*fpObjCrystInformUser)((string)buf);
   }
   this->UpdateDisplay();
   VFN_DEBUG_MESSAGE("PowderPatternBackground::ImportUserBackground():finished",5)
}
void PowderPatternBackground::SetInterpPoints(const CrystVector_REAL tth,
                                              const CrystVector_REAL backgd)
{
   VFN_DEBUG_ENTRY("PowderPatternBackground::SetInterpPoints():",5)
   if(  (tth.numElements()!=backgd.numElements())
      ||(tth.numElements()<2))
   {
      throw ObjCrystException("PowderPatternBackground::SetInterpPoints() : \
number of points differ or less than 2 points !");
   }
   mBackgroundNbPoint=tth.numElements();
   mBackgroundInterpPointX.resize(mBackgroundNbPoint);
   mBackgroundInterpPointIntensity.resize(mBackgroundNbPoint);
   // Sort in ascending order, disregarding radiation type.
   CrystVector<long> subs;
   subs=SortSubs(tth);

   for(long i=0;i<mBackgroundNbPoint;++i)
   {
      mBackgroundInterpPointX(i)=tth(subs(i));
      mBackgroundInterpPointIntensity(i)=backgd(subs(i));
   }
   this->InitRefParList();
   mClockBackgroundPoint.Click();
   VFN_DEBUG_EXIT("PowderPatternBackground::SetInterpPoints()",5)
}

const pair<const CrystVector_REAL*,const CrystVector_REAL*> PowderPatternBackground::GetInterpPoints()const
{
   return make_pair(&mBackgroundInterpPointX,&mBackgroundInterpPointIntensity);
}

void PowderPatternBackground::GetGeneGroup(const RefinableObj &obj,
                                CrystVector_uint & groupIndex,
                                unsigned int &first) const
{
   // One group for all background points
   unsigned int index=0;
   VFN_DEBUG_MESSAGE("PowderPatternBackground::GetGeneGroup()",4)
   for(long i=0;i<obj.GetNbPar();i++)
      for(long j=0;j<this->GetNbPar();j++)
         if(&(obj.GetPar(i)) == &(this->GetPar(j)))
         {
            if(index==0) index=first++;
            groupIndex(i)=index;
         }
}

void PowderPatternBackground::BeginOptimization(const bool allowApproximations,
                                                const bool enableRestraints)
{
   this->RefinableObj::BeginOptimization(allowApproximations,enableRestraints);
}

const CrystVector_REAL& PowderPatternBackground::GetPowderPatternCalcVariance()const
{
   this->CalcPowderPattern();
   return mPowderPatternCalcVariance;
}

pair<const CrystVector_REAL*,const RefinableObjClock*>
   PowderPatternBackground::GetPowderPatternIntegratedCalcVariance()const
{
   this->CalcPowderPatternIntegrated();
   return make_pair(&mPowderPatternIntegratedCalcVariance,
                    &mClockPowderPatternIntegratedVarianceCalc);
}

bool PowderPatternBackground::HasPowderPatternCalcVariance()const
{
   #ifdef USE_BACKGROUND_MAXLIKE_ERROR
   return true;
   #else
   return false;
   #endif
}

void PowderPatternBackground::TagNewBestConfig()const
{
}

void PowderPatternBackground::OptimizeBayesianBackground()
{
   VFN_DEBUG_ENTRY("PowderPatternBackground::OptimizeBayesianBackground()",5);
   TAU_PROFILE("PowderPatternBackground::OptimizeBayesianBackground()","void ()",TAU_DEFAULT);
   PowderPatternBackgroundBayesianMinimiser min(*this);
   SimplexObj simplex("Simplex Test");
   simplex.AddRefinableObj(min);

   long nbcycle;
   REAL llk=simplex.GetLogLikelihood();
   long ct=0;
   cout<<"Initial Chi^2(BayesianBackground)="<<llk<<endl;
   this->SetGlobalOptimStep(gpRefParTypeScattDataBackground,
                            mBackgroundInterpPointIntensity.max()/1000.0);

   {
      char buf [200];
      sprintf(buf,"Optimizing Background, Cycle %d, Chi^2(Background)=%f",
              (int)ct,(float)llk);
      (*fpObjCrystInformUser)((string)buf);
   }
   nbcycle=500*mBackgroundNbPoint;
   simplex.Optimize(nbcycle,false);
   llk=simplex.GetLogLikelihood();
   cout<<ct<<", Chi^2(BayesianBackground)="<<llk<<endl;

   this->SetGlobalOptimStep(gpRefParTypeScattDataBackground,10.0);
   {
      char buf [200];
      sprintf(buf,"Done Optimizing Bayesian Background, Chi^2(Background)=%f",(float)llk);
      (*fpObjCrystInformUser)((string)buf);
   }

   LSQNumObj lsq;
   lsq.SetRefinedObj(min,0,true,true);
   lsq.PrepareRefParList(true);
   lsq.SetParIsUsed(gpRefParTypeScattDataBackground,true);
   try {lsq.Refine(10,true,false);}
   catch(const ObjCrystException &except){};

   this->GetParentPowderPattern().UpdateDisplay();
   VFN_DEBUG_EXIT("PowderPatternBackground::OptimizeBayesianBackground()",5);
}

void PowderPatternBackground::FixParametersBeyondMaxresolution(RefinableObj &obj)
{
   //Auto-fix points beyond used range
   unsigned long nbpoint=this->GetParentPowderPattern().GetNbPointUsed();
   for(long j=0;j<mBackgroundNbPoint;j++)
      if(this->GetParentPowderPattern().X2Pixel(mBackgroundInterpPointX(j))>nbpoint)
      {
         obj.GetPar(&mBackgroundInterpPointIntensity(j)).Print();
         obj.GetPar(&mBackgroundInterpPointIntensity(j)).SetIsFixed(true);
      }
}

void PowderPatternBackground::CalcPowderPattern() const
{
   if(mClockPowderPatternCalc>mClockMaster) return;

   //:TODO: This needs serious optimization !
   if(   (mClockPowderPatternCalc>mClockBackgroundPoint)
       &&(mClockPowderPatternCalc>mpParentPowderPattern->GetClockPowderPatternPar())
       &&(mClockPowderPatternCalc>mInterpolationModel.GetClock())) return;
   TAU_PROFILE("PowderPatternBackground::CalcPowderPattern()","void ()",TAU_DEFAULT);
   VFN_DEBUG_MESSAGE("PowderPatternBackground::CalcPowderPattern()",3);

   const unsigned long nb=mpParentPowderPattern->GetNbPoint();
   mPowderPatternCalc.resize(nb);
   if(nb!=0)
      switch(mInterpolationModel.GetChoice())
      {
         case POWDER_BACKGROUND_LINEAR:
         {
            VFN_DEBUG_MESSAGE("PowderPatternBackground::CalcPowderPattern()..Linear",2)
            REAL p1,p2;
            REAL b1,b2;
            if(mBackgroundNbPoint==0)
            {
               mPowderPatternCalc=0;
               break;
            }
            VFN_DEBUG_MESSAGE("PowderPatternBackground::CalcPowderPattern()"<<nb,2)
            this->InitSpline();
            REAL *b=mPowderPatternCalc.data();
            p1=this->GetParentPowderPattern().X2Pixel(mBackgroundInterpPointX(mPointOrder(0)));
            p2=this->GetParentPowderPattern().X2Pixel(mBackgroundInterpPointX(mPointOrder(1)));
            b1=mBackgroundInterpPointIntensity(mPointOrder(0));
            b2=mBackgroundInterpPointIntensity(mPointOrder(1));
            long point=1;
            for(unsigned long i=0;i<nb;i++)
            {
               if(i >= p2)
               {
                  if(point < mBackgroundNbPoint-1)
                  {
                     b1=b2;
                     p1=p2;
                     b2=mBackgroundInterpPointIntensity(mPointOrder(point+1));
                     p2=this->GetParentPowderPattern().X2Pixel(mBackgroundInterpPointX(mPointOrder(point+1)));
                     point++ ;
                  }
               }
               *b = (b1*(p2-i)+b2*(i-p1))/(p2-p1) ;
               b++;
            }
            break;
         }
         case POWDER_BACKGROUND_CUBIC_SPLINE:
         {
            if(mBackgroundNbPoint==0) mPowderPatternCalc=0;
            else
            {
               this->InitSpline();
               mPowderPatternCalc=mvSpline((REAL)0,(REAL)1,nb);
            }
            break;
         }
      }
   VFN_DEBUG_MESSAGE("PowderPatternBackground::CalcPowderPattern()",3);
   #ifdef USE_BACKGROUND_MAXLIKE_ERROR
   {
      mPowderPatternCalcVariance.resize(nb);
      const REAL step=mModelVariance*mModelVariance/(REAL)nbPoint;
      REAL var=0;
      REAL *p=mPowderPatternCalcVariance.data();
      for(long i=0;i<nb;i++) {*p++ = var;var +=step;}
   }
   mClockPowderPatternVarianceCalc.Click();
   #endif
   mClockPowderPatternCalc.Click();
   VFN_DEBUG_MESSAGE("PowderPatternBackground::CalcPowderPattern():End",3);
}

void PowderPatternBackground::CalcPowderPattern_FullDeriv(std::set<RefinablePar*> &vPar)
{
   TAU_PROFILE("PowderPatternBackground::CalcPowderPattern_FullDeriv()","void ()",TAU_DEFAULT);
   const unsigned long nb=mpParentPowderPattern->GetNbPoint();
   mPowderPattern_FullDeriv.clear();
   if((nb==0)||(mBackgroundNbPoint==0)) return;
   for(std::set<RefinablePar*>::iterator par=vPar.begin();par!=vPar.end();++par)
   {
      if((*par)==0) mPowderPattern_FullDeriv[*par]=this->GetPowderPatternCalc();
      else
         for(int j = 0; j < mBackgroundNbPoint; j++)
         {
            if((*par)->GetPointer()!=mBackgroundInterpPointIntensity.data()+j) continue;
            const REAL step=(*par)->GetDerivStep();
            (*par)->Mutate(step);
            mPowderPattern_FullDeriv[*par] =this->GetPowderPatternCalc();
            (*par)->Mutate(-2*step);
            mPowderPattern_FullDeriv[*par]-=this->GetPowderPatternCalc();
            (*par)->Mutate(step);
            mPowderPattern_FullDeriv[*par]/= step*2;
            if(MaxAbs(mPowderPattern_FullDeriv[*par])==0)
               mPowderPattern_FullDeriv[*par].resize(0);
         }
   }
}

void PowderPatternBackground::CalcPowderPatternIntegrated() const
{
   if(mClockPowderPatternCalc>mClockMaster) return;

   this->CalcPowderPattern();// :TODO: Optimize
   if(  (mClockPowderPatternIntegratedCalc>mClockPowderPatternCalc)
      &&(mClockPowderPatternIntegratedCalc>mpParentPowderPattern->GetIntegratedProfileLimitsClock()))
         return;

   VFN_DEBUG_ENTRY("PowderPatternBackground::CalcPowderPatternIntegrated()",3)
   TAU_PROFILE("PowderPatternBackground::CalcPowderPatternIntegrated()","void ()",TAU_DEFAULT);
   const CrystVector_long *pMin=&(mpParentPowderPattern->GetIntegratedProfileMin());
   const CrystVector_long *pMax=&(mpParentPowderPattern->GetIntegratedProfileMax());

   const long numInterval=pMin->numElements();
   mPowderPatternIntegratedCalc.resize(numInterval);
   REAL * RESTRICT p2=mPowderPatternIntegratedCalc.data();
   for(int j=0;j<numInterval;j++)
   {
      const long max=(*pMax)(j);
      const REAL * RESTRICT p1=mPowderPatternCalc.data()+(*pMin)(j);
      *p2=0;
      for(int k=(*pMin)(j);k<=max;k++) *p2 += *p1++;
      p2++;
   }
   #ifdef USE_BACKGROUND_MAXLIKE_ERROR
   mPowderPatternIntegratedCalcVariance.resize(numInterval);
   p2=mPowderPatternIntegratedCalcVariance.data();
   for(int j=0;j<numInterval;j++)
   {
      const long max=(*pMax)(j);
      const REAL *p1=mPowderPatternCalcVariance.data()+(*pMin)(j);
      *p2=0;
      for(int k=(*pMin)(j);k<=max;k++) *p2 += *p1++;
      p2++;
   }
   mClockPowderPatternIntegratedVarianceCalc.Click();
   #endif
   mClockPowderPatternIntegratedCalc.Click();
   VFN_DEBUG_EXIT("PowderPatternBackground::CalcPowderPatternIntegrated()",3)
}

void PowderPatternBackground::CalcPowderPatternIntegrated_FullDeriv(std::set<RefinablePar*> &vPar)
{
   TAU_PROFILE("PowderPatternBackground::CalcPowderPatternIntegrated_FullDeriv()","void ()",TAU_DEFAULT);
   //cout<<"PowderPatternBackground::CalcPowderPatternIntegrated_FullDeriv"<<endl;
   const unsigned long nb=mpParentPowderPattern->GetNbPoint();
   mPowderPatternIntegrated_FullDeriv.clear();
   if((nb==0)||(mBackgroundNbPoint==0)) return;
   for(std::set<RefinablePar*>::iterator par=vPar.begin();par!=vPar.end();++par)
   {
      if((*par)==0) mPowderPattern_FullDeriv[*par]=this->GetPowderPatternCalc();
      else
         for(int j = 0; j < mBackgroundNbPoint; j++)
         {
            if((*par)->GetPointer()!=mBackgroundInterpPointIntensity.data()+j) continue;
            const REAL step=(*par)->GetDerivStep();
            (*par)->Mutate(step);
            mPowderPatternIntegrated_FullDeriv[*par] =*(this->GetPowderPatternIntegratedCalc().first);
            (*par)->Mutate(-2*step);
            mPowderPatternIntegrated_FullDeriv[*par]-=*(this->GetPowderPatternIntegratedCalc().first);
            (*par)->Mutate(step);
            mPowderPatternIntegrated_FullDeriv[*par]/= step*2;
            if(MaxAbs(mPowderPatternIntegrated_FullDeriv[*par])==0)
               mPowderPatternIntegrated_FullDeriv[*par].resize(0);
         }
   }
   #if 0
   std::map<RefinablePar*, CrystVector_REAL> newDeriv=mPowderPatternIntegrated_FullDeriv;
   this->PowderPatternComponent::CalcPowderPatternIntegrated_FullDeriv(vPar);
   std::vector<const CrystVector_REAL*> v;
   int n=0;
   for(std::map<RefinablePar*, CrystVector_REAL>::reverse_iterator pos=mPowderPatternIntegrated_FullDeriv.rbegin();pos!=mPowderPatternIntegrated_FullDeriv.rend();++pos)
   {
      cout<<pos->first->GetName()<<":"<<pos->second.size()<<","<<newDeriv[pos->first].size()<<endl;
      if(pos->second.size()==0) continue;
      v.push_back(&(pos->second));
      v.push_back(&(newDeriv[pos->first]));
      if(++n>8) break;
   }
   if(v.size()>0) cout<<"PowderPatternBackground::CalcPowderPatternIntegrated_FullDeriv():"<<endl<<FormatVertVector<REAL>(v,12,1,20)<<endl;
   //exit(0);
   #endif
}

void PowderPatternBackground::Prepare()
{
}
const CrystVector_long& PowderPatternBackground::GetBraggLimits()const
{
   // no integration interval for the background
   mIntegratedReflLimits.resize(0);
   return mIntegratedReflLimits;
}

void PowderPatternBackground::SetMaxSinThetaOvLambda(const REAL max)
{
   mMaxSinThetaOvLambda=max;
   mClockMaster.Click();
}

void PowderPatternBackground::InitRefParList()
{
   this->ResetParList();
   REAL *p=mBackgroundInterpPointIntensity.data();
   char buf [10];
   string str="Background_Point_";
   //for(int i=0;i<3;i++)
   for(int i=0;i<mBackgroundNbPoint;i++)
   {
      sprintf(buf,"%d",i);
      RefinablePar tmp(str+(string)buf,p++,
                        0.,1000.,gpRefParTypeScattDataBackground,REFPAR_DERIV_STEP_RELATIVE,
                        false,true,true,false,1.);
      tmp.SetGlobalOptimStep(10.);
      tmp.AssignClock(mClockBackgroundPoint);
      tmp.SetDerivStep(1e-3);
      this->AddPar(tmp);
   }
   #ifdef USE_BACKGROUND_MAXLIKE_ERROR
   {
      RefinablePar tmp("ML Model Error",&mModelVariance,
                        0.,100000.,gpRefParTypeObjCryst,REFPAR_DERIV_STEP_RELATIVE,
                        true,true,true,false,1.);
      tmp.AssignClock(mClockBackgroundPoint);
      tmp.SetDerivStep(1e-3);
      //tmp.SetGlobalOptimStep(10.);
      tmp.SetGlobalOptimStep(sqrt(mBackgroundInterpPointIntensity.sum()
                             /mBackgroundInterpPointIntensity.numElements()));
      this->AddPar(tmp);
   }
   #endif
}

void PowderPatternBackground::InitOptions()
{
   VFN_DEBUG_MESSAGE("PowderPatternBackground::InitOptions()",5)
   static string InterpolationModelName;
   static string InterpolationModelChoices[2];

   static bool needInitNames=true;
   if(true==needInitNames)
   {
      InterpolationModelName="Interpolation Model";
      InterpolationModelChoices[0]="Linear";
      InterpolationModelChoices[1]="Spline";
      //InterpolationModelChoices[2]="Chebyshev";

      needInitNames=false;//Only once for the class
   }
   mInterpolationModel.Init(2,&InterpolationModelName,InterpolationModelChoices);
   this->AddOption(&mInterpolationModel);
   mClockMaster.AddChild(mInterpolationModel.GetClock());
   mInterpolationModel.SetChoice(1);
}

void PowderPatternBackground::InitSpline()const
{
   if(  (mClockSpline>mClockBackgroundPoint)
      &&(mClockSpline>mpParentPowderPattern->GetClockPowderPatternPar())
      &&(mClockSpline>this->GetParentPowderPattern().GetRadiation().GetClockWavelength())) return;

   mvSplinePixel.resize(mBackgroundNbPoint);

   // The points must be in ascending order
   // Take care later of neutron TOF, as the powder apttern data may not have been initialized yet.
   CrystVector<long> subs;
   subs=SortSubs(mBackgroundInterpPointX);

   if(this->GetParentPowderPattern().GetRadiation().GetWavelengthType()==WAVELENGTH_TOF)
   {
      mPointOrder.resize(mBackgroundNbPoint);
      for(long i=0;i<mBackgroundNbPoint;++i)
         mPointOrder(i)=subs(mBackgroundNbPoint-1-i);
   }
   else mPointOrder=subs;

   CrystVector_REAL ipixel(mBackgroundNbPoint);

   for(long i=0;i<mBackgroundNbPoint;++i)
   {
      mvSplinePixel(i)=
         this->GetParentPowderPattern().X2Pixel(mBackgroundInterpPointX(mPointOrder(i)));
      ipixel(i)=mBackgroundInterpPointIntensity(mPointOrder(i));
   }

   mvSpline.Init(mvSplinePixel,ipixel);
   mClockSpline.Click();
}

#ifdef __WX__CRYST__
WXCrystObjBasic* PowderPatternBackground::WXCreate(wxWindow* parent)
{
   //:TODO: Check mpWXCrystObj==0
   mpWXCrystObj=new WXPowderPatternBackground(parent,this);
   return mpWXCrystObj;
}
#endif
////////////////////////////////////////////////////////////////////////
//
//        PowderPatternDiffraction
//
////////////////////////////////////////////////////////////////////////
PowderPatternDiffraction::PowderPatternDiffraction():
mpReflectionProfile(0),
mCorrLorentz(*this),mCorrPolar(*this),mCorrSlitAperture(*this),
mCorrTextureMarchDollase(*this),mCorrTextureEllipsoid(*this),mCorrTOF(*this),mExtractionMode(false),
mpLeBailData(0),mFrozenLatticePar(6),mFreezeLatticePar(false),mFrozenBMatrix(3,3)
{
   VFN_DEBUG_MESSAGE("PowderPatternDiffraction::PowderPatternDiffraction()",10)
   mIsScalable=true;
   this->InitOptions();
   this->SetProfile(new ReflectionProfilePseudoVoigt);
   this->SetIsIgnoringImagScattFact(true);
   this->AddSubRefObj(mCorrTextureMarchDollase);
   this->AddSubRefObj(mCorrTextureEllipsoid);
   mClockMaster.AddChild(mClockProfilePar);
   mClockMaster.AddChild(mClockLorentzPolarSlitCorrPar);
   mClockMaster.AddChild(mpReflectionProfile->GetClockMaster());
   for(unsigned int i=0;i<3;++i) mFrozenLatticePar(i)=5;
   for(unsigned int i=3;i<6;++i) mFrozenLatticePar(i)=M_PI/2;
}

PowderPatternDiffraction::PowderPatternDiffraction(const PowderPatternDiffraction &old):
mpReflectionProfile(0),
mCorrLorentz(*this),mCorrPolar(*this),mCorrSlitAperture(*this),
mCorrTextureMarchDollase(*this),mCorrTextureEllipsoid(*this),mCorrTOF(*this),mExtractionMode(false),
mpLeBailData(0),mFrozenLatticePar(6),mFreezeLatticePar(old.FreezeLatticePar()),mFrozenBMatrix(3,3)
{
   this->AddSubRefObj(mCorrTextureMarchDollase);
   this->AddSubRefObj(mCorrTextureEllipsoid);
   this->SetIsIgnoringImagScattFact(true);
   this->SetProfile(old.mpReflectionProfile->CreateCopy());
   #if 0 //:TODO:
   if(old.mpLeBailData!=0)
   {
      mpLeBailData=new DiffractionDataSingleCrystal(false);
      *mpLeBailData = *(old.mpLeBailData);
   }
   #endif
   mClockMaster.AddChild(mClockProfilePar);
   mClockMaster.AddChild(mClockLorentzPolarSlitCorrPar);
   mClockMaster.AddChild(mpReflectionProfile->GetClockMaster());
   for(unsigned int i=0;i<6;++i) mFrozenLatticePar(i)=old.GetFrozenLatticePar(i);
   mFrozenBMatrix=old.GetBMatrix();
}

PowderPatternDiffraction::~PowderPatternDiffraction()
{
   if(mpReflectionProfile!=0)
   {
      this->RemoveSubRefObj(*mpReflectionProfile);
      delete mpReflectionProfile;
   }
}
const string& PowderPatternDiffraction::GetClassName() const
{
   const static string className="PowderPatternDiffraction";
   return className;
}

PowderPatternDiffraction* PowderPatternDiffraction::CreateCopy()const
{
   return new PowderPatternDiffraction(*this);
}

void PowderPatternDiffraction::SetParentPowderPattern(PowderPattern &s)
{
   if(mpParentPowderPattern!=0)
      mClockMaster.RemoveChild(mpParentPowderPattern->GetIntegratedProfileLimitsClock());
   mpParentPowderPattern = &s;
   mClockMaster.AddChild(mpParentPowderPattern->GetIntegratedProfileLimitsClock());
   mClockMaster.AddChild(mpParentPowderPattern->GetClockPowderPatternPar());
   mClockMaster.AddChild(mpParentPowderPattern->GetClockPowderPatternXCorr());
   mClockMaster.AddChild(mpParentPowderPattern->GetClockPowderPatternRadiation());
}

const CrystVector_REAL& PowderPatternDiffraction::GetPowderPatternCalc()const
{
   this->CalcPowderPattern();
   return mPowderPatternCalc;
}

pair<const CrystVector_REAL*,const RefinableObjClock*>
   PowderPatternDiffraction::GetPowderPatternIntegratedCalc()const
{
   this->CalcPowderPatternIntegrated();
   return make_pair(&mPowderPatternIntegratedCalc,&mClockPowderPatternIntegratedCalc);
}

void PowderPatternDiffraction::SetReflectionProfilePar(const ReflectionProfileType prof,
                                                       const REAL w, const REAL u, const REAL v,
                                                       const REAL eta0, const REAL eta1)
{
   VFN_DEBUG_MESSAGE("PowderPatternDiffraction::SetReflectionProfilePar()",5)
   ReflectionProfilePseudoVoigt* p=new ReflectionProfilePseudoVoigt();
   p->SetProfilePar(w,u,v,eta0,eta1);
   this->SetProfile(p);
}

void PowderPatternDiffraction::SetProfile(ReflectionProfile *p)
{
   if(p==mpReflectionProfile) return;
   if(mpReflectionProfile!=0)
   {
      this->RemoveSubRefObj(*mpReflectionProfile);
      delete mpReflectionProfile;
   }
   mpReflectionProfile= p;
   this->AddSubRefObj(*mpReflectionProfile);
   mClockMaster.AddChild(mpReflectionProfile->GetClockMaster());
}

const ReflectionProfile& PowderPatternDiffraction::GetProfile()const
{
   return *mpReflectionProfile;
}

ReflectionProfile& PowderPatternDiffraction::GetProfile()
{
   return *mpReflectionProfile;
}

// Disable the base-class function.
void PowderPatternDiffraction::GenHKLFullSpace(
        const REAL maxTheta, const bool useMultiplicity)
{
   // This should be never called.
   abort();
}

void PowderPatternDiffraction::GenHKLFullSpace()
{
   VFN_DEBUG_ENTRY("PowderPatternDiffraction::GenHKLFullSpace():",5)
   float stol;
   if(this->GetRadiation().GetWavelengthType()==WAVELENGTH_TOF)
      stol=mpParentPowderPattern->X2STOL(mpParentPowderPattern->GetPowderPatternXMin());
   else
      stol=mpParentPowderPattern->X2STOL(mpParentPowderPattern->GetPowderPatternXMax());
   if(stol>1) stol=1; // Do not go beyond 0.5 A resolution (mostly for TOF data)
   this->ScatteringData::GenHKLFullSpace2(stol,true);
   if((mExtractionMode) && (mFhklObsSq.numElements()!=this->GetNbRefl()))
   {// Reflections changed, so ScatteringData::PrepareHKLarrays() probably reseted mFhklObsSq
      VFN_DEBUG_ENTRY("PowderPatternDiffraction::GenHKLFullSpace(): need to reset observed intensities",7)
      mFhklObsSq.resize(this->GetNbRefl());
      mFhklObsSq=100;
   }
   mCorrTextureEllipsoid.InitRefParList();// #TODO: SHould this be here ?
   VFN_DEBUG_EXIT("PowderPatternDiffraction::GenHKLFullSpace():"<<this->GetNbRefl(),5)
}
void PowderPatternDiffraction::BeginOptimization(const bool allowApproximations,
                                                 const bool enableRestraints)
{
   if(mUseFastLessPreciseFunc!=allowApproximations)
   {
      mClockProfileCalc.Reset();
      mClockGeomStructFact.Reset();
      mClockStructFactor.Reset();
      mClockMaster.Click();
   }
   mUseFastLessPreciseFunc=allowApproximations;
   this->GetNbReflBelowMaxSinThetaOvLambda();
   this->RefinableObj::BeginOptimization(allowApproximations,enableRestraints);
}
void PowderPatternDiffraction::EndOptimization()
{
   if(mOptimizationDepth==1)
   {
      if(mUseFastLessPreciseFunc==true)
      {
         mClockProfileCalc.Reset();
         mClockGeomStructFact.Reset();
         mClockStructFactor.Reset();
         mClockMaster.Click();
      }
      mUseFastLessPreciseFunc=false;
      this->GetNbReflBelowMaxSinThetaOvLambda();
   }
   this->RefinableObj::EndOptimization();
}

void PowderPatternDiffraction::SetApproximationFlag(const bool allow)
{
   if(mUseFastLessPreciseFunc!=allow)
   {
      mClockProfileCalc.Reset();
      mClockGeomStructFact.Reset();
      mClockStructFactor.Reset();
      mClockMaster.Click();
   }
   mUseFastLessPreciseFunc=allow;
   this->GetNbReflBelowMaxSinThetaOvLambda();
   this->RefinableObj::SetApproximationFlag(allow);
}

void PowderPatternDiffraction::GetGeneGroup(const RefinableObj &obj,
                                CrystVector_uint & groupIndex,
                                unsigned int &first) const
{
   // One group for all profile parameters
   unsigned int index=0;
   VFN_DEBUG_MESSAGE("PowderPatternDiffraction::GetGeneGroup()",4)
   for(long i=0;i<obj.GetNbPar();i++)
      for(long j=0;j<this->GetNbPar();j++)
         if(&(obj.GetPar(i)) == &(this->GetPar(j)))
         {
            //if(this->GetPar(j).GetType()->IsDescendantFromOrSameAs())
            //{
               if(index==0) index=first++;
               groupIndex(i)=index;
            //}
            //else //no parameters other than unit cell
         }
}

const CrystVector_REAL& PowderPatternDiffraction::GetPowderPatternCalcVariance()const
{

   return mPowderPatternCalcVariance;
}

pair<const CrystVector_REAL*,const RefinableObjClock*>
   PowderPatternDiffraction::GetPowderPatternIntegratedCalcVariance()const
{
   this->CalcPowderPatternIntegrated();
   return make_pair(&mPowderPatternIntegratedCalcVariance,
                    &mClockPowderPatternIntegratedVarianceCalc);
}

bool PowderPatternDiffraction::HasPowderPatternCalcVariance()const
{
   return true;
}

void PowderPatternDiffraction::SetCrystal(Crystal &crystal)
{
   bool reprep=(mpCrystal!=0);
   this->ScatteringData::SetCrystal(crystal);
   // Check if we use DE-PV
   if(mpReflectionProfile!=0)
      if(mpReflectionProfile->GetClassName()=="ReflectionProfileDoubleExponentialPseudoVoigt")
      {
         ReflectionProfileDoubleExponentialPseudoVoigt *p
            =dynamic_cast<ReflectionProfileDoubleExponentialPseudoVoigt*>(mpReflectionProfile);
            p->SetUnitCell((UnitCell)crystal);
      }
   mClockHKL.Reset();
   if(reprep) this->Prepare();
}

const Radiation& PowderPatternDiffraction::GetRadiation()const
{ return mpParentPowderPattern->GetRadiation();}

void PowderPatternDiffraction::SetExtractionMode(const bool extract,const bool init)
{
   VFN_DEBUG_ENTRY("PowderPatternDiffraction::SetExtractionMode(),ExtractionMode="<<mExtractionMode<<", nbrefl="<<this->GetNbRefl(),7)
   mExtractionMode=extract;
   bool needInit=false;
   if(extract)
   {
      this->FreezeLatticePar(false);
      this->Prepare();
      mFhklObsSq.resizeAndPreserve(this->GetNbRefl());
   }
   if(extract && (!init) && (mpLeBailData!=0))
   {
      // Re-use existing Le Bail data, if list of hkl's is consistent
      const long nbrefl=this->GetNbReflBelowMaxSinThetaOvLambda();
      if(nbrefl==mpLeBailData->GetNbRefl())
      {
         for(int i=0;i<nbrefl;++i)
         {
            if(  (mpLeBailData->GetH()(i)==this->GetH()(i))
               &&(mpLeBailData->GetK()(i)==this->GetK()(i))
               &&(mpLeBailData->GetL()(i)==this->GetL()(i)))
            {
               mFhklObsSq(i)=mpLeBailData->GetFhklObsSq()(i);
            }
            else
            {
               needInit=true;
               VFN_DEBUG_MESSAGE("PowderPatternDiffraction::SetExtractionMode():: Forcing initialize, cannot re-use: hkl list differs",10);
               break;
            }

         }
      }
      else
      {
        needInit=true;
        VFN_DEBUG_MESSAGE("PowderPatternDiffraction::SetExtractionMode():: Forcing initialize, cannot re-use: different number of reflections",10);
      }
   }
   if((extract && init) || needInit)
   {
      VFN_DEBUG_MESSAGE("PowderPatternDiffraction::SetExtractionMode():: Initializing intensities to 100",10);
      mFhklObsSq=100;
   }
   if((mExtractionMode==false)&&(mFhklObsSq.numElements()>0))
   {// Leaving extraction mode, so update extracted single crystal data
      VFN_DEBUG_ENTRY("PowderPatternDiffraction::SetExtractionMode(),LEAVING Le Bail Mode",7)
      if(mpLeBailData==0)  mpLeBailData=new DiffractionDataSingleCrystal(this->GetCrystal(),false);
      // Update wavelength & name
      mpLeBailData->SetWavelength(this->GetRadiation().GetWavelength()(0));
      mpLeBailData->SetRadiationType(this->GetRadiation().GetRadiationType());
      char buf[200];
      sprintf(buf,"LeBail (d=%4.2fA):",1/(2*abs(mMaxSinThetaOvLambda)+1e-6));
      mpLeBailData->SetName(string(buf)+this->GetCrystal().GetName());

      const unsigned long nbrefl=this->GetNbReflBelowMaxSinThetaOvLambda();
      CrystVector_REAL iobs(nbrefl),sigma(nbrefl);
      CrystVector_long h(nbrefl),k(nbrefl),l(nbrefl);
      sigma=1;
      for(unsigned long i=0;i<nbrefl;++i)
      {
         h(i)=mIntH(i);
         k(i)=mIntK(i);
         l(i)=mIntL(i);
         iobs(i)=mFhklObsSq(i);
      }
      mpLeBailData->SetHklIobs(h,k,l,iobs,sigma);
      // Erase mFhklObsSq - only used during extraction mode.
      mFhklObsSq.resize(0);
   }
   mClockIhklCalc.Reset();mClockMaster.Reset();
   mClockFhklObsSq.Click();
   VFN_DEBUG_EXIT("PowderPatternDiffraction::SetExtractionMode(),ExtractionMode="<<mExtractionMode<<", nbrefl="<<this->GetNbRefl(),7)
}

bool PowderPatternDiffraction::GetExtractionMode()const{return mExtractionMode;}

void PowderPatternDiffraction::ExtractLeBail(unsigned int nbcycle)
{
   VFN_DEBUG_ENTRY("PowderPatternDiffraction::ExtractLeBail()",7)
   TAU_PROFILE("PowderPatternDiffraction::ExtractLeBail()","void (int)",TAU_DEFAULT);

   if(mExtractionMode==false) this->SetExtractionMode(true,true);// Should not have to do this here !
   if(mFhklObsSq.numElements()!=this->GetNbRefl())
   {//Something went wrong !
      VFN_DEBUG_ENTRY("PowderPatternDiffraction::ExtractLeBail() mFhklObsSq.size() != NbRefl !!!!!!",7)
      mFhklObsSq.resize(this->GetNbRefl());
      mFhklObsSq=100;
   }
   // First get the observed powder pattern, minus the contribution of all other phases.
   CrystVector_REAL obs,iextract,calc;
   iextract=mFhklObsSq;
   mFhklObsSq=0;
   mClockFhklObsSq.Click();
   // Get the observed and calculated powder pattern (excluding this diffraction phase)
   obs=mpParentPowderPattern->GetPowderPatternObs();
   obs-=mpParentPowderPattern->GetPowderPatternCalc();
   mFhklObsSq=iextract;
   mClockFhklObsSq.Click();
   // We take here the reflections which are centered below the max(sin(theta)/lambda)
   // actually more reflections are calculated, but the pattern is only calculated up to
   // max(sin(theta)/lambda).
   const unsigned long nbrefl=this->ScatteringData::GetNbReflBelowMaxSinThetaOvLambda();
   iextract=0;
   for(;nbcycle>0;nbcycle--)
   {
      //cout<<"PowderPatternDiffraction::ExtractLeBail(): cycle #"<<nbcycle<<endl;
      calc=this->GetPowderPatternCalc();
      for(unsigned int k0=0;k0<nbrefl;++k0)
      {
         if(mvReflProfile[k0].profile.numElements()==0) continue; // May happen for reflections near limits ?
         REAL s1=0;
         //cout<<mH(k0)<<" "<<mK(k0)<<" "<<mL(k0)<<" , Iobs=??"<<endl;
         long last=mvReflProfile[k0].last,first;
         if(last>=(long)(mpParentPowderPattern->GetNbPointUsed())) last=mpParentPowderPattern->GetNbPointUsed();
         if(mvReflProfile[k0].first<0)first=0;
         else first=(mvReflProfile[k0].first);
         const REAL *p1=mvReflProfile[k0].profile.data()+(first-mvReflProfile[k0].first);
         const REAL *p2=calc.data()+first;
         const REAL *pobs=obs.data()+first;
         for(long i=first;i<=last;++i)
         {
            const REAL s2=*p2++;
            const REAL tmp=*pobs++ * *p1++;
            if( (s2<1e-8) ) // || (tmp<=0)
            {// Avoid <0 intensities (should not happen, it means profile is <0)
               //cout<<"S2? "<< int(mH(k0))<<" "<<int(mK(k0))<<" "<<int(mL(k0)) <<" calc(i="<<i<<")"<<calc(i)<<" obs(i="<<i<<")="<<obs(i)<<", tmp="<<tmp<<" profile(i)="<<mvReflProfile[k0].profile(i-mvReflProfile[k0].first)<<" "<<mFhklObsSq(k0)<<endl;
               continue ;
            }
            s1 += tmp /s2;
            //cout<<"   "<<s2<<" "<<obs(i)<<" "<<mvReflProfile[k0].profile(i-mvReflProfile[k0].first)<<" "<<mFhklObsSq(k0)<<endl;
         }
         if((s1>1e-8)&&(!ISNAN_OR_INF(s1))) iextract(k0)=s1*mFhklObsSq(k0);
         else iextract(k0)=1e-8;//:KLUDGE: should <0 intensities be allowed ?
         //if(nbcycle==1) cout<<" Le Bail "<<int(mH(k0))<<" "<<int(mK(k0))<<" "<<int(mL(k0))<<" , Iobs="<<iextract(k0)<<endl;
      }
      mFhklObsSq=iextract;
      if(this->GetCrystal().GetScatteringComponentList().GetNbComponent()>0)
      {// Change scale factor if we have some atoms in the structure
         const REAL* p1=this->GetFhklCalcSq() .data();
         const REAL* p2=mFhklObsSq.data();
         REAL tmp1=0,tmp2=0;
         for(long i=nbrefl;i>0;i--)
         {
            tmp1 += (*p1) * (*p2++);
            tmp2 += (*p1) * (*p1);
            p1++;
         }
         //cout<<"SCALING: tmp2="<<tmp2<<",tmp1="<<tmp1<<endl;
         mFhklObsSq*=tmp2/tmp1;
      }
      mClockFhklObsSq.Click();
      //cout<<"PowderPatternDiffraction::ExtractLeBail():results (scale factor="<<mpParentPowderPattern->GetScaleFactor(*this)*1e6<<")" <<endl<< FormatVertVectorHKLFloats<REAL>(mH,mK,mL,this->GetFhklCalcSq(),mFhklObsSq,10,4,nbrefl)<<endl;
   }
   // Store extracted data in a single crystal data object
   if(mpLeBailData==0) mpLeBailData=new DiffractionDataSingleCrystal(*mpCrystal,false);
   {
      VFN_DEBUG_MESSAGE("PowderPatternDiffraction::ExtractLeBail(): creating single crystal extracted data",7)
      CrystVector_REAL iobs(nbrefl),sigma(nbrefl);
      CrystVector_long h(nbrefl),k(nbrefl),l(nbrefl);
      sigma=1;
      for(unsigned long i=0;i<nbrefl;++i)
      {
         h(i)=mIntH(i);
         k(i)=mIntK(i);
         l(i)=mIntL(i);
         iobs(i)=mFhklObsSq(i);
      }
      mpLeBailData->SetHklIobs(h,k,l,iobs,sigma);
   }
   VFN_DEBUG_EXIT("PowderPatternDiffraction::ExtractLeBail()mFhklObsSq.size()=="<<mFhklObsSq.numElements(),7)
}
long PowderPatternDiffraction::GetNbReflBelowMaxSinThetaOvLambda()const
{
   if(this->IsBeingRefined()) return mNbReflUsed;
   VFN_DEBUG_MESSAGE("PowderPatternDiffraction::GetNbReflBelowMaxSinThetaOvLambda(): "<<mNbReflUsed<<"/"<<mNbRefl<<" [max sin(theta)/lambda="<<mMaxSinThetaOvLambda<<"]",4)
   this->CalcPowderReflProfile();
   const long nbpoint=mpParentPowderPattern->GetNbPointUsed();
   if((mNbReflUsed>0)&&(mNbReflUsed<mNbRefl))
   {
      if(  (mvReflProfile[mNbReflUsed  ].first>nbpoint)
         &&(mvReflProfile[mNbReflUsed-1].first<=nbpoint)) return mNbReflUsed;
   }

   if((mNbReflUsed==mNbRefl) && (mvReflProfile[mNbReflUsed-1].profile.numElements()>0))
      if(mvReflProfile[mNbReflUsed-1].first<=nbpoint)return mNbReflUsed;


   long i;
   for(i=0;i<mNbRefl;i++)
   {
      if(mvReflProfile[i].first>nbpoint) break;
   }
   if(i!=mNbReflUsed)
   {
      mNbReflUsed=i;
      mClockNbReflUsed.Click();
      VFN_DEBUG_MESSAGE("->Changed Max sin(theta)/lambda="<<mMaxSinThetaOvLambda\
                        <<" nb refl="<<mNbReflUsed,4)
   }
   return mNbReflUsed;
}

void PowderPatternDiffraction::SetFrozenLatticePar(const unsigned int i, REAL v)
{
   const REAL old=mFrozenLatticePar(i);
   cout<<"PowderPatternDiffraction::SetFrozenLatticePar("<<i<<":"<<v<<")"<<endl;
   if(old==v) return;
   mFrozenLatticePar(i)=v;
   this->CalcFrozenBMatrix();
}

REAL PowderPatternDiffraction::GetFrozenLatticePar(const unsigned int i) const {return mFrozenLatticePar(i);}

void PowderPatternDiffraction::FreezeLatticePar(const bool use)
{
   cout<<"PowderPatternDiffraction::FreezeLatticePar("<<use<<")"<<endl;
   VFN_DEBUG_MESSAGE("PowderPatternDiffraction::FreezeLatticePar("<<use<<")", 10)
   if(use==mFreezeLatticePar) return;
   mFreezeLatticePar=use;
   mFrozenLatticePar=this->GetCrystal().GetLatticePar();
   if(use) this->CalcFrozenBMatrix();
   mClockTheta.Reset();
   this->UpdateDisplay();
}

bool PowderPatternDiffraction::FreezeLatticePar() const {return mFreezeLatticePar;}

unsigned int PowderPatternDiffraction::GetProfileFitNetNbObs()const
{
   unsigned int nb=0;
   unsigned int irefl=0;
   unsigned int ilast=0;
   while(this->GetParentPowderPattern().STOL2Pixel(mSinThetaLambda(irefl))<0)
   {
      irefl++;
      if(irefl>=this->GetNbReflBelowMaxSinThetaOvLambda()) break;
   }
   REAL stol=mSinThetaLambda(irefl);
   while(irefl<this->GetNbReflBelowMaxSinThetaOvLambda())
   {
      while(mSinThetaLambda(irefl)==stol)
      {
         //cout<<int(mH(irefl))<<" "<<int(mK(irefl))<<" "<<int(mL(irefl))<<endl;
         irefl++;
         if(irefl>=this->GetNbReflBelowMaxSinThetaOvLambda()) break;
      }
      const int nbnew =this->GetParentPowderPattern().STOL2Pixel(stol)-ilast;
      if(nbnew>1) nb += nbnew-1;
      //cout<<"     => Added "<< nbnew-1<< "net observed points ("<<this->GetParentPowderPattern().STOL2Pixel(stol)<<"-"<<ilast<<")"<<endl;
      ilast=this->GetParentPowderPattern().STOL2Pixel(stol);
      stol = mSinThetaLambda(irefl);
   }
   //cout<<"Final number of net observed points: "<<nb<<endl;
   return nb;
}

void PowderPatternDiffraction::CalcPowderPattern() const
{
   this->GetNbReflBelowMaxSinThetaOvLambda();
   if(mClockPowderPatternCalc>mClockMaster) return;
   TAU_PROFILE("PowderPatternDiffraction::CalcPowderPattern()-Apply profiles","void (bool)",TAU_DEFAULT);

   VFN_DEBUG_ENTRY("PowderPatternDiffraction::CalcPowderPattern():",3)

   // :TODO: Can't do this as this is non-const
   //if(this->GetCrystal().GetSpaceGroup().GetClockSpaceGroup()>mClockHKL)
   //   this->GenHKLFullSpace();
   //
   // The workaround is to call Prepare() (non-const) before every calculation
   // when a modifictaion may have occured.

   this->CalcIhkl();
   this->CalcPowderReflProfile();

   if(  (mClockPowderPatternCalc>mClockIhklCalc)
      &&(mClockPowderPatternCalc>mClockProfileCalc)) return;

   if(true) //:TODO: false == mUseFastLessPreciseFunc
   {
      const long nbRefl=this->GetNbRefl();
      VFN_DEBUG_MESSAGE("PowderPatternDiffraction::CalcPowderPattern\
Applying profiles for "<<nbRefl<<" reflections",2)
      long step; // number of reflections at the same place and with the same (assumed) profile
      const long  specNbPoints=mpParentPowderPattern->GetNbPoint();
      mPowderPatternCalc.resize(specNbPoints);
      mPowderPatternCalc=0;
      const bool useML= (mIhklCalcVariance.numElements() != 0);
      if(useML)
      {
         mPowderPatternCalcVariance.resize(specNbPoints);
         mPowderPatternCalcVariance=0;
      }
      else mPowderPatternCalcVariance.resize(0);
      VFN_DEBUG_MESSAGE("PowderPatternDiffraction::CalcPowderPattern() Has variance:"<<useML,2)

      for(long i=0;i<mNbRefl;i += step)
      {
         if(mvReflProfile[i].profile.numElements()==0)
         {
            step=1;
            if(i>=mNbReflUsed) break;// After sin(theta)/lambda limit
            else continue; // before beginning of pattern ?
         }

         VFN_DEBUG_MESSAGE("PowderPatternDiffraction::CalcPowderPattern()#"<<i,2)
         REAL intensity=0.;
         REAL var=0.;
         //check if the next reflection is at the same theta. If this is true,
         //Then assume that the profile is exactly the same, unless it is anisotropic
         for(step=0; ;)
         {
            intensity += mIhklCalc(i + step);
            if(useML) var += mIhklCalcVariance(i + step);
            step++;
            if(mpReflectionProfile->IsAnisotropic()) break;// Anisotropic profiles
            if( (i+step) >= nbRefl) break;
            if(mSinThetaLambda(i+step) > (mSinThetaLambda(i)+1e-5) ) break;
         }
         VFN_DEBUG_MESSAGE("Apply profile(Monochromatic)Refl("<<i<<")"\
            <<mIntH(i)<<" "<<mIntK(i)<<" "<<mIntL(i)<<" "\
            <<"  I="<<intensity<<"  stol="<<mSinThetaLambda(i)\
            <<",pixel #"<<mvReflProfile[i].first<<"->"<<mvReflProfile[i].last,2)
         {
            const long first=mvReflProfile[i].first,last=mvReflProfile[i].last;
            const REAL *p2 = mvReflProfile[i].profile.data();
            REAL *p3 = mPowderPatternCalc.data()+first;
            for(long j=first;j<=last;j++) *p3++ += *p2++ * intensity;
            if(useML)
            {
               const REAL *p2 = mvReflProfile[i].profile.data();
               REAL *p3 = mPowderPatternCalcVariance.data()+first;
               for(long j=first;j<=last;j++) *p3++ += *p2++ * var;
            }
         }
      }
   }
   else
   {
      //:TODO:
      throw ObjCrystException("PowderPatternDiffraction::CalcPowderPattern() : \
         FAST option not yet implemented !");
   }
   mClockPowderPatternCalc.Click();
   VFN_DEBUG_EXIT("PowderPatternDiffraction::CalcPowderPattern: End.",3)
}

void PowderPatternDiffraction::CalcPowderPattern_FullDeriv(std::set<RefinablePar*> &vPar)
{
   TAU_PROFILE("PowderPatternDiffraction::CalcPowderPattern_FullDeriv()","void ()",TAU_DEFAULT);
   //cout<<"PowderPatternDiffraction::CalcPowderPattern_FullDeriv"<<endl;
   this->CalcPowderPattern();
   bool notYetDerivProfiles=true;
   mIhkl_FullDeriv.clear();
   mvReflProfile_FullDeriv.clear();
   mPowderPattern_FullDeriv.clear();
   for(std::set<RefinablePar*>::iterator par=vPar.begin();par!=vPar.end();++par)
   {
      if(*par==0) continue;
      if((*par)->IsFixed()) continue;
      if((*par)->IsUsed()==false) continue;
      if(  (*par)->GetType()->IsDescendantFromOrSameAs(gpRefParTypeScatt)
         ||(*par)->GetType()->IsDescendantFromOrSameAs(gpRefParTypeScattPow))
      {
         this->CalcIhkl_FullDeriv(vPar);
      }
      if(notYetDerivProfiles)
      {
         if(  (*par)->GetType()->IsDescendantFromOrSameAs(gpRefParTypeRadiation)
            ||(*par)->GetType()->IsDescendantFromOrSameAs(gpRefParTypeUnitCell)
            ||(*par)->GetType()->IsDescendantFromOrSameAs(gpRefParTypeScattDataCorrPos)
            ||(*par)->GetType()->IsDescendantFromOrSameAs(gpRefParTypeScattDataProfile))
         {
            this->CalcPowderReflProfile_FullDeriv(vPar);
            notYetDerivProfiles=false;
         }
      }
   }

   //this->CalcPowderReflProfile();
   for(std::set<RefinablePar*>::iterator par=vPar.begin();par!=vPar.end();++par)
   {
      if(*par==0) mPowderPattern_FullDeriv[*par]=this->GetPowderPatternCalc();
      else
      {
         if((*par)->IsFixed()) continue;
         if((*par)->IsUsed()==false) continue;
         if(mIhkl_FullDeriv[*par].size()!=0)
         {
            const long nbRefl=this->GetNbRefl();
            long step; // number of reflections at the same place and with the same (assumed) profile
            const long  specNbPoints=mpParentPowderPattern->GetNbPoint();
            mPowderPattern_FullDeriv[*par].resize(specNbPoints);
            mPowderPattern_FullDeriv[*par]=0;

            for(long i=0;i<mNbReflUsed;i += step)
            {
               if(mvReflProfile[i].profile.numElements()==0)
               {
                  step=1;
                  if(i>=mNbReflUsed) break;
                  else continue;
               }

               REAL intensity=0.;
               //check if the next reflection is at the same theta. If this is true,
               //Then assume that the profile is exactly the same, unless it is anisotropic
               for(step=0; ;)
               {
                  intensity += mIhkl_FullDeriv[*par](i + step);
                  step++;
                  if(mpReflectionProfile->IsAnisotropic()) break;// Anisotropic profiles
                  if( (i+step) >= nbRefl) break;
                  if(mSinThetaLambda(i+step) > (mSinThetaLambda(i)+1e-5) ) break;
               }
               {
                  const long first=mvReflProfile[i].first,last=mvReflProfile[i].last;
                  const REAL *p2 = mvReflProfile[i].profile.data();
                  REAL *p3 = mPowderPattern_FullDeriv[*par].data()+first;
                  for(long j=first;j<=last;j++) *p3++ += *p2++ * intensity;
               }
            }
         }
         if(mvReflProfile_FullDeriv[*par].size()!=0)
         {
            const long nbRefl=this->GetNbRefl();
            long step; // number of reflections at the same place and with the same (assumed) profile
            const long  specNbPoints=mpParentPowderPattern->GetNbPoint();
            mPowderPattern_FullDeriv[*par].resize(specNbPoints);
            mPowderPattern_FullDeriv[*par]=0;// :TODO: use only the number of points actually used
            cout<<__FILE__<<":"<<__LINE__<<":PowderPatternDiffraction::CalcPowderPattern_FullDeriv():par="<<(*par)->GetName()<<endl;
            for(long i=0;i<mNbReflUsed;i += step)
            {
               if(mvReflProfile[i].profile.numElements()==0)
               {
                  step=1;
                  if(i>=mNbReflUsed) break;
                  else continue;
               }

               REAL intensity=0.;
               //check if the next reflection is at the same theta. If this is true,
               //Then assume that the profile is exactly the same, unless it is anisotropic
               for(step=0; ;)
               {
                  intensity += mIhklCalc(i + step);
                  step++;
                  if(mpReflectionProfile->IsAnisotropic()) break;// Anisotropic profiles
                  if( (i+step) >= nbRefl) break;
                  if(mSinThetaLambda(i+step) > (mSinThetaLambda(i)+1e-5) ) break;
               }
               if(mvReflProfile_FullDeriv[*par][i].size()>0)// Some profiles may be unaffected by a given parameter
               {
                  const long first=mvReflProfile[i].first,last=mvReflProfile[i].last;
                  const REAL *p2 = mvReflProfile_FullDeriv[*par][i].data();
                  REAL *p3 = mPowderPattern_FullDeriv[*par].data()+first;
                  for(long j=first;j<=last;j++) *p3++ += *p2++ * intensity;
               }
            }
         }
      }
   }
   #if 0
   std::map<RefinablePar*, CrystVector_REAL> newDeriv=mPowderPattern_FullDeriv;
   this->PowderPatternComponent::CalcPowderPattern_FullDeriv(vPar);
   std::vector<const CrystVector_REAL*> v;
   int n=0;
   for(std::map<RefinablePar*, CrystVector_REAL>::reverse_iterator pos=mPowderPattern_FullDeriv.rbegin();pos!=mPowderPattern_FullDeriv.rend();++pos)
   {
      if(pos->first==0) continue;
      if(pos->second.size()==0) continue;
      v.push_back(&(newDeriv[pos->first]));
      v.push_back(&(pos->second));
      cout<<pos->first->GetName()<<":"<<pos->second.size()<<","<<newDeriv[pos->first].size()<<endl;
      if(++n>8) break;
   }
   cout<<FormatVertVector<REAL>(v,16,4,1000)<<endl;
   //exit(0);
   #endif
}

void PowderPatternDiffraction::CalcPowderPatternIntegrated() const
{
   this->GetNbReflBelowMaxSinThetaOvLambda();
   if(mClockPowderPatternIntegratedCalc>mClockMaster) return;
   TAU_PROFILE("PowderPatternDiffraction::CalcPowderPatternIntegrated()","void (bool)",TAU_DEFAULT);
   TAU_PROFILE_TIMER(timer1,"PowderPatternDiffraction::CalcPowderPatternIntegrated()1","", TAU_FIELD);
   TAU_PROFILE_TIMER(timer2,"PowderPatternDiffraction::CalcPowderPatternIntegrated()2","", TAU_FIELD);

   this->CalcIhkl();
   TAU_PROFILE_START(timer1);
   this->PrepareIntegratedProfile();
   TAU_PROFILE_STOP(timer1);

   if(  (mClockPowderPatternIntegratedCalc>mClockIhklCalc)
      &&(mClockPowderPatternIntegratedCalc>mClockIntegratedProfileFactor)
      &&(mClockPowderPatternIntegratedCalc>mpParentPowderPattern->GetIntegratedProfileLimitsClock()))
            return;
   VFN_DEBUG_ENTRY("PowderPatternDiffraction::CalcPowderPatternIntegrated()",3)
   const long nbRefl=this->GetNbRefl();

   const long  nb=mpParentPowderPattern->GetIntegratedProfileMin().numElements();
   mPowderPatternIntegratedCalc.resize(nb);
   mPowderPatternIntegratedCalc=0;
   const bool useML= (mIhklCalcVariance.numElements() != 0);
   if(useML)
   {
      mPowderPatternIntegratedCalcVariance.resize(nb);
      mPowderPatternIntegratedCalcVariance=0;
   }
   else mPowderPatternIntegratedCalcVariance.resize(0);
   const REAL * RESTRICT psith=mSinThetaLambda.data();
   const REAL * RESTRICT pI=mIhklCalc.data();
   const REAL * RESTRICT pIvar=mIhklCalcVariance.data();
   vector< pair<unsigned long, CrystVector_REAL> >::const_iterator pos;
   pos=mIntegratedProfileFactor.begin();
   TAU_PROFILE_START(timer2);
   for(long i=0;i<mNbReflUsed;)
   {
      VFN_DEBUG_MESSAGE("PowderPatternDiffraction::CalcPowderPatternIntegrated():"<<i,2)
      REAL intensity=0.;
      REAL var=0.;
      const REAL thmax=*psith+1e-5;
      //check if the next reflection is at the same theta. If this is true,
      //Then assume that the profile is exactly the same, unless profiles are anisotropic.
      for(;;)
      {
         intensity += *pI++;
         if(useML) var += *pIvar++;
         if( ++i >= nbRefl) break;
         if( *(++psith) > thmax ) break;
         if(mpReflectionProfile->IsAnisotropic()) break;// Anisotropic profile
         ++pos;
      }
      VFN_DEBUG_MESSAGE("PowderPatternDiffraction::CalcPowderPatternIntegrated():"<<i,2)
      REAL * RESTRICT pData=mPowderPatternIntegratedCalc.data()+pos->first;
      const REAL * RESTRICT pFact=pos->second.data();
      const unsigned long nb=pos->second.numElements();
      //cout <<i<<" - "<< intensity<<"*:";
      for(unsigned long j=nb;j>0;j--)
      {
         //cout <<pos->first+j<<"("<<*pFact<<","<<*pData<<") ";
         *pData++ += intensity * *pFact++ ;
      }
      //cout<<endl;

      if(useML)
      {
         const REAL * RESTRICT pFact=pos->second.data();
         REAL * RESTRICT pVar=mPowderPatternIntegratedCalcVariance.data()+pos->first;
         for(unsigned long j=nb;j>0;j--) *pVar++ += var * *pFact++ ;
      }
      ++pos;
   }
   TAU_PROFILE_STOP(timer2);
   #ifdef __DEBUG__
   if(gVFNDebugMessageLevel<3)
   {
      this->CalcPowderPattern();
      CrystVector_REAL integr(nb),min(nb),max(nb),diff(nb),index(nb);
      integr=0;
      for(long i=0;i<nb;i++)
      {
         index(i)=i;
         min(i)=mpParentPowderPattern->GetIntegratedProfileMin()(i);
         max(i)=mpParentPowderPattern->GetIntegratedProfileMax()(i);
         integr(i)=0;
         for(long j=mpParentPowderPattern->GetIntegratedProfileMin()(i);
                 j<=mpParentPowderPattern->GetIntegratedProfileMax()(i);j++)
         {
            integr(i) += mPowderPatternCalc(j);
         }
         diff(i)=1.-mPowderPatternIntegratedCalc(i)/integr(i);
      }
      cout << "Integrated intensities, Component"<<endl
           << FormatVertVectorHKLFloats<REAL>
                (index,min,max,integr,mPowderPatternIntegratedCalc,diff,20,6)<<endl;
   }
   #endif
   mClockPowderPatternIntegratedCalc.Click();
   VFN_DEBUG_EXIT("PowderPatternDiffraction::CalcPowderPatternIntegrated",3)
}

void PowderPatternDiffraction::CalcPowderPatternIntegrated_FullDeriv(std::set<RefinablePar*> &vPar)
{
   TAU_PROFILE("PowderPatternDiffraction::CalcPowderPatternIntegrated_FullDeriv()","void ()",TAU_DEFAULT);
   //cout<<"PowderPatternDiffraction::CalcPowderPatternIntegrated_FullDeriv"<<endl;
   //this->PowderPatternComponent::CalcPowderPatternIntegrated_FullDeriv(vPar);
   //return;

   this->CalcPowderPatternIntegrated();
   this->CalcIhkl_FullDeriv(vPar);
   const long nbRefl=this->GetNbRefl();
   //#define PowderPatternDiffraction_CalcPowderPatternIntegrated_FullDerivDEBUG
   #ifdef PowderPatternDiffraction_CalcPowderPatternIntegrated_FullDerivDEBUG
   this->PowderPatternComponent::CalcPowderPatternIntegrated_FullDeriv(vPar);
   std::map<RefinablePar*, CrystVector_REAL> oldDeriv=mPowderPatternIntegrated_FullDeriv;
   #endif
   const long  nbprof=mpParentPowderPattern->GetIntegratedProfileMin().size();
   long ctpar=0;
   mPowderPatternIntegrated_FullDeriv.clear();
   for(std::set<RefinablePar*>::iterator par=vPar.begin();par!=vPar.end();++par)
   {
      if(*par==0) mPowderPatternIntegrated_FullDeriv[*par]=mPowderPatternIntegratedCalc;
      else
      {
         if(mIhkl_FullDeriv[*par].size()==0) continue;
         if(mPowderPatternIntegrated_FullDeriv[*par].size()==0)
         {
            mPowderPatternIntegrated_FullDeriv[*par].resize(nbprof);
            mPowderPatternIntegrated_FullDeriv[*par]=0;
         }
         const REAL * RESTRICT psith=mSinThetaLambda.data();
         const REAL * RESTRICT pI=mIhkl_FullDeriv[*par].data();
         // :TODO: also handle derivatives of profile parameters ? Though one should not
         // refine profiles using integrated profiles !
         vector< pair<unsigned long, CrystVector_REAL> >::const_iterator pos=mIntegratedProfileFactor.begin();

         for(long i=0;i<mNbReflUsed;)
         {
            REAL intensity=0.;
            const REAL thmax=*psith+1e-5;
            for(;;)
            {
               intensity += *pI++;
               if( ++i >= nbRefl) break;
               if( *(++psith) > thmax ) break;
               if(mpReflectionProfile->IsAnisotropic()) break;// Anisotropic profile
               ++pos;
            }
            REAL * RESTRICT pData=mPowderPatternIntegrated_FullDeriv[*par].data()+pos->first;
            const REAL * RESTRICT pFact=pos->second.data();
            const unsigned long nb=pos->second.numElements();
            #ifdef PowderPatternDiffraction_CalcPowderPatternIntegrated_FullDerivDEBUG
            if((i<5)&&(ctpar<8)) cout<<__FILE__<<":"<<__LINE__<<":"<<(*par)->GetName()<<"i="<<setw(16)<<i<<":I="<<setw(16)<<intensity<<endl;
            #endif
            for(unsigned long j=nb;j>0;j--)
            {
               #ifdef PowderPatternDiffraction_CalcPowderPatternIntegrated_FullDerivDEBUG
               *pData += intensity * *pFact ;
               if((i<5)&&((*par)->GetName()=="Cimetidine_C11_x")&&(pos->first==0)&&(nb==j)) cout<<nb-j<<" SUM1"<<setw(16)<<*pData<<", dI="<<setw(16)<<intensity<<", prof="<<setw(16)<<*pFact<<endl;
               pData++;pFact++ ;
               #else
               *pData++ += intensity * *pFact++ ;
               #endif
            }
            ++pos;
            ctpar++;
         }
      }
   }
   #ifdef PowderPatternDiffraction_CalcPowderPatternIntegrated_FullDerivDEBUG
   std::vector<const CrystVector_REAL*> v;
   int n=0;
   cout<<"PowderPatternDiffraction::CalcPowderPatternIntegrated_FullDeriv():parameters:"<<endl;
   for(std::set<RefinablePar*>::iterator par=vPar.begin();par!=vPar.end();++par)
   {
      if(mPowderPatternIntegrated_FullDeriv[*par].size()==0) continue;
      v.push_back(&(mPowderPatternIntegrated_FullDeriv[*par]));
      v.push_back(&(oldDeriv[*par]));
      cout<<(*par)->GetName()<<":"<<mPowderPatternIntegrated_FullDeriv[*par].size()<<","<<oldDeriv[*par].size()<<endl;
      if(++n>6) break;
   }
   cout<<"PowderPatternDiffraction::CalcPowderPatternIntegrated_FullDeriv():"<<endl<<FormatVertVector<REAL>(v,12,1,20)<<endl;
   //exit(0);
   #endif
}

void PowderPatternDiffraction::CalcPowderReflProfile()const
{
   this->CalcSinThetaLambda();
   mpParentPowderPattern->GetNbPointUsed();//
   if(  (mClockProfileCalc>mClockProfilePar)
      &&(mClockProfileCalc>mpReflectionProfile->GetClockMaster())
      &&(mClockProfileCalc>mClockTheta)
      &&(mClockProfileCalc>this->GetRadiation().GetClockWavelength())
      &&(mClockProfileCalc>mpParentPowderPattern->GetClockPowderPatternXCorr())
      &&(mClockProfileCalc>mClockHKL)
      &&(mClockProfileCalc>mpParentPowderPattern->GetClockNbPointUsed())) return;

   TAU_PROFILE("PowderPatternDiffraction::CalcPowderReflProfile()","void (bool)",TAU_DEFAULT);
   VFN_DEBUG_ENTRY("PowderPatternDiffraction::CalcPowderReflProfile()",5)

   //Calc all profiles
   mvLabel.clear();
   stringstream label;

   unsigned int nbLine=1;
   CrystVector_REAL spectrumDeltaLambdaOvLambda;
   CrystVector_REAL spectrumFactor;//relative weigths of different lines of X-Ray tube
   switch(this->GetRadiation().GetWavelengthType())
   {
      case WAVELENGTH_MONOCHROMATIC:
      {
         spectrumDeltaLambdaOvLambda.resize(1);spectrumDeltaLambdaOvLambda=0.0;
         spectrumFactor.resize(1);spectrumFactor=1.0;
         break;
      }
      case WAVELENGTH_ALPHA12:
      {
         nbLine=2;
         spectrumDeltaLambdaOvLambda.resize(2);
         spectrumDeltaLambdaOvLambda(0)
            =-this->GetRadiation().GetXRayTubeDeltaLambda()
             *this->GetRadiation().GetXRayTubeAlpha2Alpha1Ratio()
             /(1+this->GetRadiation().GetXRayTubeAlpha2Alpha1Ratio())
             /this->GetRadiation().GetWavelength()(0);
         spectrumDeltaLambdaOvLambda(1)
            = this->GetRadiation().GetXRayTubeDeltaLambda()
             /(1+this->GetRadiation().GetXRayTubeAlpha2Alpha1Ratio())
             /this->GetRadiation().GetWavelength()(0);

         spectrumFactor.resize(2);
         spectrumFactor(0)=1./(1.+this->GetRadiation().GetXRayTubeAlpha2Alpha1Ratio());
         spectrumFactor(1)=this->GetRadiation().GetXRayTubeAlpha2Alpha1Ratio()
                           /(1.+this->GetRadiation().GetXRayTubeAlpha2Alpha1Ratio());
         break;
      }
      case WAVELENGTH_TOF:
      {
         spectrumDeltaLambdaOvLambda.resize(1);spectrumDeltaLambdaOvLambda=0.0;
         spectrumFactor.resize(1);spectrumFactor=1.0;
         break;
      }
      default: throw ObjCrystException("PowderPatternDiffraction::PrepareIntegratedProfile():\
Radiation must be either monochromatic, from an X-Ray Tube, or neutron TOF !!");
   }


   VFN_DEBUG_MESSAGE("PowderPatternDiffraction::CalcPowderReflProfile():\
Computing all Profiles",5)
   REAL center,// center of current reflection (depends on line if several)
        x0;    // theoretical (uncorrected for zero's, etc..) position of center of line
   long first,last;// first & last point of the stored profile
   CrystVector_REAL vx,reflProfile,tmpV;
   mvReflProfile.resize(this->GetNbRefl());
   for(unsigned int i=0;i<this->GetNbRefl();i++)
   {
      mvReflProfile[i].first=0;
      mvReflProfile[i].last=0;
      mvReflProfile[i].profile.resize(0);
   }
   VFN_DEBUG_MESSAGE("PowderPatternDiffraction::CalcPowderReflProfile()",5)

   for(unsigned int line=0;line<nbLine;line++)
   {
      for(long i=0;i<this->GetNbRefl();i++)
      {// Only the reflections contributing below the max(sin(theta)/lambda) will be computed
         VFN_DEBUG_ENTRY("PowderPatternDiffraction::CalcPowderReflProfile()#"<<i,5)
         x0=mpParentPowderPattern->STOL2X(mSinThetaLambda(i));

         VFN_DEBUG_MESSAGE("PowderPatternDiffraction::CalcPowderReflProfile()#"<<i,5)
         if(nbLine>1)
         {// we have several lines, not centered on the profile range
            center = mpParentPowderPattern->X2XCorr(
                        x0+2*tan(x0/2.0)*spectrumDeltaLambdaOvLambda(line));
         }
         else center=mpParentPowderPattern->X2XCorr(x0);
         REAL fact=1.0;
         if(!mUseFastLessPreciseFunc) fact=5.0;
         const REAL halfwidth=mpReflectionProfile->GetFullProfileWidth(0.04,center,mH(i),mK(i),mL(i))*fact;
         if(line==0)
         {
            // For an X-Ray tube, label on first (strongest) of reflections lines (Kalpha1)
            label.str("");
            label<<mIntH(i)<<" "<<mIntK(i)<<" "<<mIntL(i);
            mvLabel.push_back(make_pair(center,label.str()));
            REAL spectrumwidth=0.0;
            if(this->GetRadiation().GetWavelengthType()==WAVELENGTH_ALPHA12)
            {// We need to shift the last point to include 2 lines in the profile
               spectrumwidth=2*this->GetRadiation().GetXRayTubeDeltaLambda()
                              /this->GetRadiation().GetWavelength()(0)*tan(x0/2.0);
            }
            first=(long)(mpParentPowderPattern->X2Pixel(center-halfwidth));
            last =(long)(mpParentPowderPattern->X2Pixel(center+halfwidth+spectrumwidth));
            if(this->GetRadiation().GetWavelengthType()==WAVELENGTH_TOF)
            {
               const long f=first;
               first=last;
               last=f;
            }
            if(first>last)
            { // Whoops - should not happen !! Unless there is a strange (dis)order for the x coordinates...
               cout<<"PowderPatternDiffraction::CalcPowderReflProfile(), line"<<__LINE__<<"first>last !! :"<<first<<","<<last<<endl;
               first=(first+last)/2;
               last=first;
            }
            first -=1;
            last+=1;
            VFN_DEBUG_MESSAGE("PowderPatternDiffraction::CalcPowderReflProfile():"<<first<<","<<last<<","<<center,3)
            if(first>last)
            {
               cout<<__FILE__<<__LINE__<<endl;
               exit(0);
            }
            if((last>=0)&&(first<(long)(mpParentPowderPattern->GetNbPoint())))
            {
               if(first<0) first=0;
               if(last>=(long)(mpParentPowderPattern->GetNbPoint()))
                  last=mpParentPowderPattern->GetNbPoint()-1;
               vx.resize(last-first+1);
            }
            else vx.resize(0); // store no profile if reflection out of pattern
            mvReflProfile[i].first=first;
            mvReflProfile[i].last=last;
         }
         else
         {
            first=mvReflProfile[i].first;
            last=mvReflProfile[i].last;
            if((last>=0)&&(first<(long)(mpParentPowderPattern->GetNbPoint())))
               vx.resize(last-first+1);
            else vx.resize(0);
            vx.resize(last-first+1);
         }
         if((last>=0)&&(first<(long)(mpParentPowderPattern->GetNbPoint())))
         {
            {
               const REAL *p0=mpParentPowderPattern->GetPowderPatternX().data()+first;
               REAL *p1=vx.data();
               for(long i=first;i<=last;i++) *p1++ = *p0++;
            }

            VFN_DEBUG_MESSAGE("PowderPatternDiffraction::CalcPowderReflProfile():"<<first<<","<<last<<","<<center,3)
            reflProfile=mpReflectionProfile->GetProfile(vx,center,mH(i),mK(i),mL(i));
            VFN_DEBUG_MESSAGE("PowderPatternDiffraction::CalcPowderReflProfile()",2)
            if(nbLine>1) reflProfile *=spectrumFactor(line);
            if(line==0) mvReflProfile[i].profile = reflProfile;
            else mvReflProfile[i].profile += reflProfile;
         }
         else
         { // reflection is out of pattern, so store no profile
            mvReflProfile[i].profile.resize(0);
         }
         VFN_DEBUG_EXIT("PowderPatternDiffraction::CalcPowderReflProfile():\
Computing all Profiles: Reflection #"<<i,5)
         if(first>(long)(mpParentPowderPattern->GetNbPointUsed())) break;
      }
   }
   mClockProfileCalc.Click();
   VFN_DEBUG_EXIT("PowderPatternDiffraction::CalcPowderReflProfile()",5)
}

void PowderPatternDiffraction::CalcPowderReflProfile_FullDeriv(std::set<RefinablePar *> &vPar)
{
   TAU_PROFILE("PowderPatternDiffraction::CalcPowderReflProfile_FullDeriv()","void (bool)",TAU_DEFAULT);
   cout<<__FILE__<<":"<<__LINE__<<":PowderPatternDiffraction::CalcPowderReflProfile_FullDeriv()"<<endl;
   this->CalcPowderReflProfile();
   unsigned int nbLine=1;
   CrystVector_REAL spectrumDeltaLambdaOvLambda;
   CrystVector_REAL spectrumFactor;//relative weigths of different lines of X-Ray tube
   switch(this->GetRadiation().GetWavelengthType())
   {
      case WAVELENGTH_MONOCHROMATIC:
      {
         spectrumDeltaLambdaOvLambda.resize(1);spectrumDeltaLambdaOvLambda=0.0;
         spectrumFactor.resize(1);spectrumFactor=1.0;
         break;
      }
      case WAVELENGTH_ALPHA12:
      {
         nbLine=2;
         spectrumDeltaLambdaOvLambda.resize(2);
         spectrumDeltaLambdaOvLambda(0)
            =-this->GetRadiation().GetXRayTubeDeltaLambda()
             *this->GetRadiation().GetXRayTubeAlpha2Alpha1Ratio()
             /(1+this->GetRadiation().GetXRayTubeAlpha2Alpha1Ratio())
             /this->GetRadiation().GetWavelength()(0);
         spectrumDeltaLambdaOvLambda(1)
            = this->GetRadiation().GetXRayTubeDeltaLambda()
             /(1+this->GetRadiation().GetXRayTubeAlpha2Alpha1Ratio())
             /this->GetRadiation().GetWavelength()(0);

         spectrumFactor.resize(2);
         spectrumFactor(0)=1./(1.+this->GetRadiation().GetXRayTubeAlpha2Alpha1Ratio());
         spectrumFactor(1)=this->GetRadiation().GetXRayTubeAlpha2Alpha1Ratio()
                           /(1.+this->GetRadiation().GetXRayTubeAlpha2Alpha1Ratio());
         break;
      }
      case WAVELENGTH_TOF:
      {
         spectrumDeltaLambdaOvLambda.resize(1);spectrumDeltaLambdaOvLambda=0.0;
         spectrumFactor.resize(1);spectrumFactor=1.0;
         break;
      }
      default: throw ObjCrystException("PowderPatternDiffraction::CalcPowderReflProfile_FullDeriv():\
Radiation must be either monochromatic, from an X-Ray Tube, or neutron TOF !!");
   }
   REAL center,// center of current reflection (depends on line if several)
        x0;    // theoretical (uncorrected for zero's, etc..) position of center of line
   long first,last;// first & last point of the stored profile
   CrystVector_REAL vx,reflProfile,tmpV;

   // Derivative vs the shift of the reflection center
   vector<CrystVector_REAL> vReflProfile_DerivCenter(mNbReflUsed);

   mvReflProfile_FullDeriv.clear();
   for(std::set<RefinablePar*>::iterator par=vPar.begin();par!=vPar.end();++par)
   {
      if(*par==0) continue;
      if(  (*par)->GetType()->IsDescendantFromOrSameAs(gpRefParTypeRadiation)
         ||(*par)->GetType()->IsDescendantFromOrSameAs(gpRefParTypeUnitCell)
         ||(*par)->GetType()->IsDescendantFromOrSameAs(gpRefParTypeScattDataCorrPos)
         ||(*par)->GetType()->IsDescendantFromOrSameAs(gpRefParTypeScattDataProfile))
      {
         mvReflProfile_FullDeriv[*par].resize(mNbReflUsed);

         for(unsigned int line=0;line<nbLine;line++)
         {
            for(long i=0;i<mNbReflUsed;i++)
            {
               x0=mpParentPowderPattern->STOL2X(mSinThetaLambda(i));

               if(nbLine>1)
               {// we have several lines, not centered on the profile range
                  center = mpParentPowderPattern->X2XCorr(
                              x0+2*tan(x0/2.0)*spectrumDeltaLambdaOvLambda(line));
               }
               else center=mpParentPowderPattern->X2XCorr(x0);
               REAL fact=1.0;
               if(!mUseFastLessPreciseFunc) fact=5.0;

               first=mvReflProfile[i].first;
               last=mvReflProfile[i].last;
               if((last>=0)&&(first<(long)(mpParentPowderPattern->GetNbPoint())))
                  vx.resize(last-first+1);
               else vx.resize(0);
               vx.resize(last-first+1);
               if((last>=0)&&(first<(long)(mpParentPowderPattern->GetNbPoint())))
               {
                  {
                     const REAL *p0=mpParentPowderPattern->GetPowderPatternX().data()+first;
                     REAL *p1=vx.data();
                     for(long i=first;i<=last;i++) *p1++ = *p0++;
                  }

                  if((*par)->GetType()->IsDescendantFromOrSameAs(gpRefParTypeScattDataProfile))
                  {// Parameter only affects profile
                     //if(i==0) cout<<"PowderPatternDiffraction::CalcPowderReflProfile_FullDeriv()par="<<(*par)->GetName()<<":refl #"<<i<<endl;
                     //:TODO: analytical derivatives
                     const REAL step=(*par)->GetDerivStep();
                     (*par)->Mutate(step);
                     reflProfile=mpReflectionProfile->GetProfile(vx,center,mH(i),mK(i),mL(i));
                     (*par)->Mutate(-2*step);
                     reflProfile-=mpReflectionProfile->GetProfile(vx,center,mH(i),mK(i),mL(i));
                     (*par)->Mutate(step);
                     reflProfile/=2*step;
                  }
                  else
                  {// Parameter affects reflection center
                     REAL dcenter=0;
                     {
                        //:TODO: analytical derivatives
                        const REAL step=(*par)->GetDerivStep();
                        (*par)->Mutate(step);
                        REAL x1=mpParentPowderPattern->STOL2X(this->CalcSinThetaLambda(mH(i),mK(i),mL(i)));
                        if(nbLine>1) dcenter = mpParentPowderPattern->X2XCorr(x1+2*tan(x1/2.0)*spectrumDeltaLambdaOvLambda(line));
                        else         dcenter = mpParentPowderPattern->X2XCorr(x1);
                        (*par)->Mutate(-2*step);
                        x1=mpParentPowderPattern->STOL2X(this->CalcSinThetaLambda(mH(i),mK(i),mL(i)));
                        if(nbLine>1) dcenter-= mpParentPowderPattern->X2XCorr(x1+2*tan(x1/2.0)*spectrumDeltaLambdaOvLambda(line));
                        else         dcenter-= mpParentPowderPattern->X2XCorr(x1);
                        (*par)->Mutate(step);
                        dcenter/=2*step;
                     }

                     if(dcenter!=0)
                     {
                        //if(i==0) cout<<"PowderPatternDiffraction::CalcPowderReflProfile_FullDeriv()par="<<(*par)->GetName()<<":refl #"<<i<<", dcenter="<<setw(8)<<dcenter<<endl;
                        if(vReflProfile_DerivCenter[i].size()==0)
                        {
                           const REAL step=1e-4;//:TODO: adapt for TOF
                           vReflProfile_DerivCenter[i] =mpReflectionProfile->GetProfile(vx,center+step,mH(i),mK(i),mL(i));
                           vReflProfile_DerivCenter[i]-=mpReflectionProfile->GetProfile(vx,center-step,mH(i),mK(i),mL(i));
                           vReflProfile_DerivCenter[i]/=2*step;
                        }
                        reflProfile=vReflProfile_DerivCenter[i];
                        reflProfile*=dcenter;
                     }
                     else
                     {
                        //if(i==0) cout<<"PowderPatternDiffraction::CalcPowderReflProfile_FullDeriv()par="<<(*par)->GetName()<<":refl #"<<i<<" => Parameter affects nothing ?"<<endl;
                        reflProfile.resize(0);
                     }
                  }
                  if(reflProfile.size()>0)
                  {
                     if(nbLine>1) reflProfile *=spectrumFactor(line);
                     if(line==0) mvReflProfile_FullDeriv[*par][i] = reflProfile;
                     else mvReflProfile_FullDeriv[*par][i] += reflProfile;
                  }
               }
            }
         }
      }
   }
}

void PowderPatternDiffraction::CalcIntensityCorr()const
{
   bool needRecalc=false;

   this->CalcSinThetaLambda();
   if((mClockIntensityCorr<mClockTheta)||(mClockIntensityCorr<this->GetClockNbReflBelowMaxSinThetaOvLambda())) needRecalc=true;

   const CrystVector_REAL *mpCorr[5] = {0, 0, 0, 0, 0};

   if(this->GetRadiation().GetWavelengthType()==WAVELENGTH_TOF)
   {
      mpCorr[0]=&(mCorrTOF.GetCorr());
      if(mClockIntensityCorr<mCorrTOF.GetClockCorr()) needRecalc=true;
   }
   else
   {
      mpCorr[0]=&(mCorrLorentz.GetCorr());
      if(mClockIntensityCorr<mCorrLorentz.GetClockCorr()) needRecalc=true;

      if(this->GetRadiation().GetRadiationType()==RAD_XRAY)
      {
         mpCorr[1]=&(mCorrPolar.GetCorr());
         if(mClockIntensityCorr<mCorrPolar.GetClockCorr()) needRecalc=true;
      }

      mpCorr[2]=&(mCorrSlitAperture.GetCorr());
      if(mClockIntensityCorr<mCorrSlitAperture.GetClockCorr()) needRecalc=true;
   }

   if(mCorrTextureMarchDollase.GetNbPhase()>0)
   {
      mpCorr[3]=&(mCorrTextureMarchDollase.GetCorr());
      if(mClockIntensityCorr<mCorrTextureMarchDollase.GetClockCorr()) needRecalc=true;
   }
   mpCorr[4]=&(mCorrTextureEllipsoid.GetCorr());
   if(mClockIntensityCorr<mCorrTextureEllipsoid.GetClockCorr()) needRecalc=true;


   if(needRecalc==false) return;

   TAU_PROFILE("PowderPatternDiffraction::CalcIntensityCorr()","void ()",TAU_DEFAULT);
   VFN_DEBUG_MESSAGE("PowderPatternDiffraction::CalcIntensityCorr()",2)
   mIntensityCorr.resize(mNbRefl);
   REAL *pCorr=mIntensityCorr.data();
   const REAL *p=mpCorr[0]->data();
   for(long i=mNbReflUsed;i>0;i--) *pCorr++ = *p++;
   if(this->GetRadiation().GetWavelengthType()!=WAVELENGTH_TOF)
   {
      if(this->GetRadiation().GetRadiationType()==RAD_XRAY)
      {
         pCorr=mIntensityCorr.data();
         p=mpCorr[1]->data();
         const REAL* p2=mpCorr[2]->data();
         for(long i=mNbReflUsed;i>0;i--) *pCorr++ *= *p++ * *p2++;
      }
      else
      {
         pCorr=mIntensityCorr.data();
         p=mpCorr[2]->data();
         for(long i=mNbReflUsed;i>0;i--) *pCorr++ *= *p++;
      }
   }
   if(mCorrTextureMarchDollase.GetNbPhase()>0)
   {
      pCorr=mIntensityCorr.data();
      p=mpCorr[3]->data();
      for(long i=mNbReflUsed;i>0;i--) *pCorr++ *= *p++;
   }
   if(mpCorr[4]->numElements()>0)
   {
      pCorr=mIntensityCorr.data();
      p=mpCorr[4]->data();
      for(long i=mNbReflUsed;i>0;i--) *pCorr++ *= *p++;
   }
   mClockIntensityCorr.Click();
   VFN_DEBUG_MESSAGE("PowderPatternDiffraction::CalcIntensityCorr():finished",2)
}

void PowderPatternDiffraction::CalcIhkl() const
{
   this->CalcIntensityCorr();
   if(mExtractionMode==true)
   {
      VFN_DEBUG_MESSAGE("PowderPatternDiffraction::CalcIhkl():"<<mFhklObsSq.numElements()<<","<<mIntensityCorr.numElements()<<","<<mMultiplicity.numElements(),7);
      mIhklCalc=mFhklObsSq;
      mIhklCalc*=mIntensityCorr;
      mIhklCalc*=mMultiplicity;
      mClockIhklCalc.Click();
      return;
   }
   this->CalcStructFactor();
   if(  (mClockIhklCalc>mClockIntensityCorr)
      &&(mClockIhklCalc>mClockStructFactor)
      &&(mClockIhklCalc>mClockNbReflUsed)) return;

   VFN_DEBUG_MESSAGE("PowderPatternDiffraction::CalcIhkl()",3)
   TAU_PROFILE("PowderPatternDiffraction::CalcIhkl()","void ()",TAU_DEFAULT);
   const REAL * RESTRICT pr,* RESTRICT pi,* RESTRICT pcorr;
   const int * RESTRICT mult;
   REAL * RESTRICT p;

   pr=mFhklCalcReal.data();
   pi=mFhklCalcImag.data();
   pcorr=mIntensityCorr.data();

   mult=mMultiplicity.data();
   mIhklCalc.resize(mNbRefl);
   p=mIhklCalc.data();
   if(mFhklCalcVariance.numElements()>0)
   {
      const REAL * RESTRICT pv=mFhklCalcVariance.data();
      for(long i=mNbReflUsed;i>0;i--)
      {
         *p++ = *mult++ * (*pr * *pr + *pi * *pi + 2 * *pv++) * *pcorr++;
         pr++;
         pi++;
      }
   }
   else
   {
      for(long i=mNbReflUsed;i>0;i--)
      {
         *p++ = *mult++ * (*pr * *pr + *pi * *pi) * *pcorr++;
         pr++;
         pi++;
      }
   }
   if(mFhklCalcVariance.numElements()==0)
   {
      VFN_DEBUG_MESSAGE("PowderPatternDiffraction::CalcIhkl(): No Calc Variance",2)
      mIhklCalcVariance.resize(0);
      VFN_DEBUG_MESSAGE(endl<<
                     FormatVertVectorHKLFloats<REAL>(mH,mK,mL,mSinThetaLambda,
                                                      mFhklCalcReal,
                                                      mFhklCalcImag,
                                                      mIhklCalc,
                                                      mIntensityCorr
                                                      ),2)
   }
   else
   {
      VFN_DEBUG_MESSAGE("PowderPatternDiffraction::CalcIhkl(): Calc Variance",2)
      mIhklCalcVariance.resize(mNbRefl);
      REAL * RESTRICT pVar2=mIhklCalcVariance.data();

      const REAL * RESTRICT pInt=mIhklCalc.data();
      const REAL * RESTRICT pVar=mFhklCalcVariance.data();
      pcorr=mIntensityCorr.data();
      mult=mMultiplicity.data();

      for(long j=mNbReflUsed;j>0;j--)
      {
         *pVar2++ = (4* *mult) * *pcorr * *pVar *(*pInt++ - (*mult * *pcorr) * *pVar);
         pVar++;mult++;pcorr++;
      }
      VFN_DEBUG_MESSAGE(endl<<
                     FormatVertVectorHKLFloats<REAL>(mH,mK,mL,mSinThetaLambda,
                                                   mFhklCalcReal,
                                                   mIhklCalc,
                                                   mExpectedIntensityFactor,
                                                   mIntensityCorr,
                                                   mMultiplicity,
                                                   mvLuzzatiFactor[&(mpCrystal->GetScatteringPowerRegistry().GetObj(0))],
                                                   mFhklCalcVariance,
                                                   mIhklCalcVariance),2);
      VFN_DEBUG_MESSAGE(mNbRefl<<" "<<mNbReflUsed,2)
   }

   //cout <<FormatVertVector<REAL>(mTheta,mIhklCalc,mMultiplicity,
   //                               mFhklCalcReal,mFhklCalcImag,mIntensityCorr);
   mClockIhklCalc.Click();
   VFN_DEBUG_MESSAGE("PowderPatternDiffraction::CalcIhkl():End",3)
}

void PowderPatternDiffraction::CalcIhkl_FullDeriv(std::set<RefinablePar*> &vPar)
{
   TAU_PROFILE("PowderPatternDiffraction::CalcIhkl_FullDeriv()","void ()",TAU_DEFAULT);
   //cout<<"PowderPatternDiffraction::CalcIhkl_FullDeriv()"<<endl;
   this->CalcIntensityCorr();//:TODO: derivatives of intensity corrections (Texture, displacement parameters,...)
   mIhkl_FullDeriv.clear();
   if(mExtractionMode==true)
   {
      //:TODO: handle Pawley refinements of I(hkl)
      return;
   }
   this->CalcStructFactor_FullDeriv(vPar);

   for(std::set<RefinablePar*>::iterator par=vPar.begin();par!=vPar.end();++par)
   {
      if(*par==0) mIhkl_FullDeriv[*par]=mIhklCalc;
      else
      {
         if(mFhklCalcReal_FullDeriv[*par].size()==0) continue;
         const REAL * RESTRICT pr,* RESTRICT pi,* RESTRICT prd,* RESTRICT pid,* RESTRICT pcorr;
         const int * RESTRICT mult;
         REAL * RESTRICT p;

         pr=mFhklCalcReal.data();
         pi=mFhklCalcImag.data();
         prd=mFhklCalcReal_FullDeriv[*par].data();
         pid=mFhklCalcImag_FullDeriv[*par].data();
         pcorr=mIntensityCorr.data();//:TODO: derivatives of intensity corrections (Texture, displacement parameters,...)

         mult=mMultiplicity.data();
         mIhkl_FullDeriv[*par].resize(mNbRefl);
         p=mIhkl_FullDeriv[*par].data();
         if(mFhklCalcImag_FullDeriv[*par].size()==0)
            for(long i=mNbReflUsed;i>0;i--) *p++ = *mult++ * 2 * *pr++ * *prd++ * *pcorr++;
         else
            for(long i=mNbReflUsed;i>0;i--) *p++ = *mult++ * 2 *(*pr++ * *prd++ + *pi++ * *pid++) * *pcorr++;
      }
   }
   #if 0
   std::map<RefinablePar*, CrystVector_REAL> oldDeriv;
   std::vector<const CrystVector_REAL*> v;
   v.push_back(&mH);
   v.push_back(&mK);
   v.push_back(&mL);
   CrystVector_REAL m;
   m=mMultiplicity;
   v.push_back(&m);
   v.push_back(&mIntensityCorr);
   int n=0;
   for(std::set<RefinablePar*>::iterator par=vPar.begin();par!=vPar.end();++par)
   {
      if((*par)==0) continue;
      if(mIhkl_FullDeriv[*par].size()==0) continue;

      const REAL step=(*par)->GetDerivStep();
      (*par)->Mutate(step);
      this->CalcIhkl();
      oldDeriv[*par]=mIhklCalc;
      (*par)->Mutate(-2*step);
      this->CalcIhkl();
      oldDeriv[*par]-=mIhklCalc;
      oldDeriv[*par]/=2*step;
      (*par)->Mutate(step);

      v.push_back(&(mIhkl_FullDeriv[*par]));
      v.push_back(&(oldDeriv[*par]));
      cout<<(*par)->GetName()<<":"<<mIhkl_FullDeriv[*par].size()<<","<<oldDeriv[*par].size()<<",  step="<<setw(16)<<step<<endl;
      if(++n>5) break;
   }
   cout<<"PowderPatternDiffraction::CalcIhkl_FullDeriv():"<<endl<<FormatVertVectorHKLFloats<REAL>(v,12,1,20)<<endl;
   //exit(0);
   #endif
}

void PowderPatternDiffraction::Prepare()
{
   if(  (this->GetCrystal().GetSpaceGroup().GetClockSpaceGroup()>mClockHKL)
      ||(this->GetCrystal().GetClockLatticePar()>mClockHKL)
      ||(this->GetRadiation().GetClockWavelength()>mClockHKL)
      ||(mpParentPowderPattern->GetClockPowderPatternPar()>mClockHKL))
         this->GenHKLFullSpace();
   //if(0==this->GetNbRefl()) this->GenHKLFullSpace();
}
void PowderPatternDiffraction::InitOptions()
{
   VFN_DEBUG_MESSAGE("PowderPatternDiffraction::InitOptions()",5)
   #if 0
   static string ReflectionProfileTypeName;
   static string ReflectionProfileTypeChoices[3];

   static bool needInitNames=true;
   if(true==needInitNames)
   {
      ReflectionProfileTypeName="Profile Type";
      ReflectionProfileTypeChoices[0]="Gaussian";
      ReflectionProfileTypeChoices[1]="Lorentzian";
      ReflectionProfileTypeChoices[2]="Pseudo-Voigt";

      needInitNames=false;//Only once for the class
   }
   mReflectionProfileType.Init(3,&ReflectionProfileTypeName,ReflectionProfileTypeChoices);
   this->AddOption(&mReflectionProfileType);
   #endif
}
const CrystVector_long& PowderPatternDiffraction::GetBraggLimits()const
{
   this->CalcPowderReflProfile();
   if((mClockProfileCalc>mClockBraggLimits)&&(this->GetNbReflBelowMaxSinThetaOvLambda()>0))
   {
      VFN_DEBUG_ENTRY("PowderPatternDiffraction::GetBraggLimits(*min,*max)",3)
      TAU_PROFILE("PowderPatternDiffraction::GetBraggLimits()","void ()",TAU_DEFAULT);
      mIntegratedReflLimits.resize(this->GetNbReflBelowMaxSinThetaOvLambda());
      long i = 0;
      mIntegratedReflLimits(i)=mvReflProfile[0].first;
      for(;i<(this->GetNbReflBelowMaxSinThetaOvLambda()-1);++i)
         mIntegratedReflLimits(i+1)=(mvReflProfile[i].first+mvReflProfile[i].last+mvReflProfile[i+1].first+mvReflProfile[i+1].last)/4;
      mIntegratedReflLimits(i)=mvReflProfile[i].last;
      mClockBraggLimits.Click();
      VFN_DEBUG_EXIT("PowderPatternDiffraction::GetBraggLimits(*min,*max)",3)
   }
   return mIntegratedReflLimits;
}

void PowderPatternDiffraction::SetMaxSinThetaOvLambda(const REAL max)
{this->ScatteringData::SetMaxSinThetaOvLambda(max);}

const CrystMatrix_REAL& PowderPatternDiffraction::GetBMatrix()const
{
   if(mFreezeLatticePar) return mFrozenBMatrix;
   return this->ScatteringData::GetBMatrix();
}

void PowderPatternDiffraction::CalcFrozenBMatrix()const
{
   VFN_DEBUG_MESSAGE("PowderPatternDiffraction::CalcFrozenBMatrix()", 10)
   REAL a,b,c,alpha,beta,gamma;//direct space parameters
   REAL aa,bb,cc,alphaa,betaa,gammaa;//reciprocal space parameters
   REAL v;//volume of the unit cell
   a=mFrozenLatticePar(0);
   b=mFrozenLatticePar(1);
   c=mFrozenLatticePar(2);
   alpha=mFrozenLatticePar(3);
   beta=mFrozenLatticePar(4);
   gamma=mFrozenLatticePar(5);

   v=sqrt(1-cos(alpha)*cos(alpha)-cos(beta)*cos(beta)-cos(gamma)*cos(gamma)
          +2*cos(alpha)*cos(beta)*cos(gamma));

   aa=sin(alpha)/a/v;
   bb=sin(beta )/b/v;
   cc=sin(gamma)/c/v;

   alphaa=acos( (cos(beta )*cos(gamma)-cos(alpha))/sin(beta )/sin(gamma) );
   betaa =acos( (cos(alpha)*cos(gamma)-cos(beta ))/sin(alpha)/sin(gamma) );
   gammaa=acos( (cos(alpha)*cos(beta )-cos(gamma))/sin(alpha)/sin(beta ) );

   mFrozenBMatrix = aa ,  bb*cos(gammaa) , cc*cos(betaa) ,
                   0  , bb*sin(gammaa) ,-cc*sin(betaa)*cos(alpha),
                   0  , 0              ,1/c;
}

void PowderPatternDiffraction::PrepareIntegratedProfile()const
{
   this->CalcPowderReflProfile();

   if(  (mClockIntegratedProfileFactor>mClockProfileCalc)
      &&(mClockIntegratedProfileFactor>mpParentPowderPattern->GetIntegratedProfileLimitsClock())
      &&(mClockIntegratedProfileFactor>mClockNbReflUsed))
   return;
   VFN_DEBUG_ENTRY("PowderPatternDiffraction::PrepareIntegratedProfile()",7)
   TAU_PROFILE("PowderPatternDiffraction::PrepareIntegratedProfile()","void ()",TAU_DEFAULT);
   const CrystVector_long *pMin=&(mpParentPowderPattern->GetIntegratedProfileMin());
   const CrystVector_long *pMax=&(mpParentPowderPattern->GetIntegratedProfileMax());

   const long numInterval=pMin->numElements();

   vector< map<long, REAL> > vIntegratedProfileFactor;
   vIntegratedProfileFactor.resize(mNbReflUsed);
   vector< map<long, REAL> >::iterator pos1;
   pos1=vIntegratedProfileFactor.begin();

   mIntegratedProfileFactor.resize(mNbReflUsed);
   vector< pair<unsigned long, CrystVector_REAL> >::iterator pos2;
   pos2=mIntegratedProfileFactor.begin();
   for(long i=0;i<mNbReflUsed;i++)
   {
      pos1->clear();
      long firstInterval=numInterval;
      for(long j=0;j<numInterval;j++)
      {
         const long first0 = mvReflProfile[i].first;
         const long last0  = mvReflProfile[i].last ;
         const long first= first0>(*pMin)(j) ? first0:(*pMin)(j);
         const long last = last0 <(*pMax)(j) ? last0 :(*pMax)(j);
         if((first<=last) && (mvReflProfile[i].profile.size()>0))
         {
            if(firstInterval>j) firstInterval=j;
            if(pos1->find(j) == pos1->end()) (*pos1)[j]=0.;
            REAL *fact = &((*pos1)[j]);//this creates the 'j' entry if necessary
            const REAL *p2 = mvReflProfile[i].profile.data()+(first-first0);
            //cout << i<<","<<j<<","<<first<<","<<last<<":"<<*fact<<"/"<<mNbReflUsed<<","<<mNbRefl<<endl;
            for(long k=first;k<=last;k++) *fact += *p2++;
         }
      }
      pos2->first=firstInterval;
      pos2->second.resize(pos1->size());
      REAL *pFact=pos2->second.data();
      for(map<long, REAL>::const_iterator pos=pos1->begin();pos!=pos1->end();++pos)
         *pFact++ = pos->second;
      pos1++;
      pos2++;
   }
   mClockIntegratedProfileFactor.Click();
   #ifdef __DEBUG__
   if(gVFNDebugMessageLevel<3)
   {
      unsigned long i=0;
      for(vector< pair<unsigned long, CrystVector_REAL> >::const_iterator
            pos=mIntegratedProfileFactor.begin();
            pos!=mIntegratedProfileFactor.end();++pos)
      {
         cout <<"Integrated profile factors for reflection #"<<i++<<"  ";
         for(int j=0;j<pos->second.numElements();++j)
            cout << j+pos->first<<"("<<pos->second(j)<<")  ";
         cout<<endl;
      }
   }
   #endif

   VFN_DEBUG_EXIT("PowderPatternDiffraction::PrepareIntegratedProfile()",7)
}

#ifdef __WX__CRYST__
WXCrystObjBasic* PowderPatternDiffraction::WXCreate(wxWindow* parent)
{
   //:TODO: Check mpWXCrystObj==0
   mpWXCrystObj=new WXPowderPatternDiffraction(parent,this);
   return mpWXCrystObj;
}
#endif

////////////////////////////////////////////////////////////////////////
//
//        PowderPattern
//
////////////////////////////////////////////////////////////////////////
ObjRegistry<PowderPattern>
   gPowderPatternRegistry("List of all PowderPattern objects");

PowderPattern::PowderPattern():
mIsXAscending(true),mNbPoint(0),
mXZero(0.),m2ThetaDisplacement(0.),m2ThetaTransparency(0.),
mDIFC(48277.14),mDIFA(-6.7),
mScaleFactor(20),mUseFastLessPreciseFunc(false),
mStatisticsExcludeBackground(false),mMaxSinThetaOvLambda(10),mNbPointUsed(0)
{
   mScaleFactor=1;
   mSubObjRegistry.SetName("SubObjRegistry for a PowderPattern object");
   mPowderPatternComponentRegistry.SetName("Powder Pattern Components");
   this->AddSubRefObj(mRadiation);
   this->Init();
   gPowderPatternRegistry.Register(*this);
   gTopRefinableObjRegistry.Register(*this);
   mClockMaster.AddChild(mClockPowderPatternPar);
   mClockMaster.AddChild(mClockNbPointUsed);
   mClockMaster.AddChild(mClockPowderPatternXCorr);
   mClockMaster.AddChild(mClockScaleFactor);
   mClockMaster.AddChild(mClockPowderPatternRadiation);
}

PowderPattern::PowderPattern(const PowderPattern &old):
mIsXAscending(old.mIsXAscending),mNbPoint(old.mNbPoint),
mRadiation(old.mRadiation),
mXZero(old.mXZero),m2ThetaDisplacement(old.m2ThetaDisplacement),
m2ThetaTransparency(old.m2ThetaTransparency),
mDIFC(old.mDIFC),mDIFA(old.mDIFA),
mPowderPatternComponentRegistry(old.mPowderPatternComponentRegistry),
mScaleFactor(old.mScaleFactor),
mUseFastLessPreciseFunc(old.mUseFastLessPreciseFunc),
mStatisticsExcludeBackground(old.mStatisticsExcludeBackground),
mMaxSinThetaOvLambda(old.mMaxSinThetaOvLambda),mNbPointUsed(old.mNbPointUsed)
{
   mX=old.mX;
   this->Init();
   mSubObjRegistry.SetName("SubObjRegistry for a PowderPattern :"+mName);
   gPowderPatternRegistry.Register(*this);
   gTopRefinableObjRegistry.Register(*this);
   this->AddSubRefObj(mRadiation);
   mClockMaster.AddChild(mClockPowderPatternPar);
   mClockMaster.AddChild(mClockNbPointUsed);
   mClockMaster.AddChild(mClockPowderPatternXCorr);
   mClockMaster.AddChild(mClockScaleFactor);
   mClockMaster.AddChild(mClockPowderPatternRadiation);
}

PowderPattern::~PowderPattern()
{
   gPowderPatternRegistry.DeRegister(*this);
   for(int i=0;i<mPowderPatternComponentRegistry.GetNb();i++)
   {
      mPowderPatternComponentRegistry.GetObj(i).DeRegisterClient(*this);
      this->RemoveSubRefObj(mPowderPatternComponentRegistry.GetObj(i));
      delete &(mPowderPatternComponentRegistry.GetObj(i));
   }
   gTopRefinableObjRegistry.DeRegister(*this);
}
const string& PowderPattern::GetClassName() const
{
   const static string className="PowderPattern";
   return className;
}

void PowderPattern::AddPowderPatternComponent(PowderPatternComponent &comp)
{
   VFN_DEBUG_ENTRY("PowderPattern::AddPowderPatternComponent():"<<comp.GetName(),5)
   comp.SetParentPowderPattern(*this);
   this->AddSubRefObj(comp);
   comp.RegisterClient(*this);
   mClockPowderPatternCalc.Reset();
   mClockIntegratedFactorsPrep.Reset();
   mPowderPatternComponentRegistry.Register(comp);
   //:TODO: check if there are enough scale factors
   //mScaleFactor.resizeAndPreserve(mPowderPatternComponentRegistry.GetNb());
   mScaleFactor(mPowderPatternComponentRegistry.GetNb()-1)=1.;
   mClockScaleFactor.Click();
   if(comp.IsScalable())
   {//Init refinable parameter
      RefinablePar tmp("Scale_"+comp.GetName(),mScaleFactor.data()+mPowderPatternComponentRegistry.GetNb()-1,
                        1e-10,1e10,gpRefParTypeScattDataScale,REFPAR_DERIV_STEP_RELATIVE,
                        false,true,true,false,1.);
      tmp.SetGlobalOptimStep(0.);
      tmp.AssignClock(mClockScaleFactor);
      tmp.SetDerivStep(1e-4);
      this->AddPar(tmp);
   }
   //this->UpdateDisplay();
   VFN_DEBUG_EXIT("PowderPattern::AddPowderPatternComponent():"<<comp.GetName(),5)
}

unsigned int PowderPattern::GetNbPowderPatternComponent()const
{
   return mPowderPatternComponentRegistry.GetNb();
}

const PowderPatternComponent& PowderPattern::GetPowderPatternComponent
                                                   (const string &name)const
{
   return mPowderPatternComponentRegistry.GetObj(name);
}

const PowderPatternComponent& PowderPattern::GetPowderPatternComponent
                                                   (const int i) const
{
   return mPowderPatternComponentRegistry.GetObj(i);
}

PowderPatternComponent& PowderPattern::GetPowderPatternComponent
                                                   (const string &name)
{
   return mPowderPatternComponentRegistry.GetObj(name);
}

PowderPatternComponent& PowderPattern::GetPowderPatternComponent
                                                   (const int i)
{
   return mPowderPatternComponentRegistry.GetObj(i);
}

REAL PowderPattern::GetScaleFactor(const int i)const{return mScaleFactor(i);}

REAL PowderPattern::GetScaleFactor(const PowderPatternComponent &comp)const
{
   unsigned int i=0;
   for(;i<mPowderPatternComponentRegistry.GetNb();++i)
   {
      if(&(mPowderPatternComponentRegistry.GetObj(i))==&comp) break;
   }
   if(i==mPowderPatternComponentRegistry.GetNb())
      throw ObjCrystException("PowderPattern::GetScaleFactor(comp) : no such component");
   return mScaleFactor(i);
}
void PowderPattern::SetScaleFactor(const int i, REAL s){ mScaleFactor(i)=s;}

void PowderPattern::SetScaleFactor(const PowderPatternComponent &comp, REAL s)
{
   unsigned int i=0;
   for(;i<mPowderPatternComponentRegistry.GetNb();++i)
   {
      if(&(mPowderPatternComponentRegistry.GetObj(i))==&comp) break;
   }
   if(i==mPowderPatternComponentRegistry.GetNb())
      throw ObjCrystException("PowderPattern::GetScaleFactor(comp) : no such component");
   mScaleFactor(i)=s;
}

void PowderPattern::SetPowderPatternPar(const REAL min,
                                        const REAL step,
                                        unsigned long nbPoint)
{
   VFN_DEBUG_MESSAGE("PowderPattern::SetPowderPatternPar():"<<min<<","<<step<<","<<nbPoint,3)
   mNbPoint=nbPoint;
   mX.resize(mNbPoint);
   for(unsigned long i=0;i<mNbPoint;i++) mX(i)=min+step*i;
   mPowderPatternObs.resizeAndPreserve(mNbPoint);
   mPowderPatternObsSigma.resizeAndPreserve(mNbPoint);
   mPowderPatternWeight.resizeAndPreserve(mNbPoint);
   mClockPowderPatternPar.Click();
}
void PowderPattern::SetPowderPatternX(const CrystVector_REAL &x)
{
   mNbPoint=x.numElements();
   if(&x != &mX) mX=x;
   mPowderPatternObs.resizeAndPreserve(mNbPoint);
   mPowderPatternObsSigma.resizeAndPreserve(mNbPoint);
   mPowderPatternWeight.resizeAndPreserve(mNbPoint);
   mClockPowderPatternPar.Click();
   if(mX(mNbPoint-1)>mX(0))mIsXAscending=true;
   else mIsXAscending=false;
   VFN_DEBUG_MESSAGE("PowderPattern::SetPowderPatternX() is ascending="<<mIsXAscending,5)
}

unsigned long PowderPattern::GetNbPoint()const {return mNbPoint;}

unsigned long PowderPattern::GetNbPointUsed()const
{
   if(!this->IsBeingRefined()) this->CalcNbPointUsed();
   return mNbPointUsed;
}

const RefinableObjClock& PowderPattern::GetClockNbPointUsed()const{return mClockNbPointUsed;}

void PowderPattern::SetRadiation(const Radiation &radiation)
{
   mRadiation=radiation;
   mClockPowderPatternRadiation.Click();
}
const Radiation& PowderPattern::GetRadiation()const {return mRadiation;}

Radiation& PowderPattern::GetRadiation() {return mRadiation;}

void PowderPattern::SetRadiationType(const RadiationType rad)
{
   mRadiation.SetRadiationType(rad);
}

RadiationType PowderPattern::GetRadiationType()const {return mRadiation.GetRadiationType();}
void PowderPattern::SetWavelength(const REAL lambda)
{
   VFN_DEBUG_MESSAGE("PowderPattern::SetWavelength(lambda)",3)
   mRadiation.SetWavelength(lambda);
}

void PowderPattern::SetWavelength(const string &XRayTubeElementName,const REAL alpha12ratio)
{
   VFN_DEBUG_MESSAGE("PowderPattern::SetWavelength(wavelength)",3)
   mRadiation.SetWavelength(XRayTubeElementName,alpha12ratio);
}

REAL PowderPattern::GetWavelength()const{return mRadiation.GetWavelength()(0);}

const CrystVector_REAL& PowderPattern::GetPowderPatternCalc()const
{
   this->CalcPowderPattern();
   return mPowderPatternCalc;
}

std::map<RefinablePar*,CrystVector_REAL>& PowderPattern::GetPowderPattern_FullDeriv(std::set<RefinablePar *> &vPar)
{
   this->CalcPowderPattern_FullDeriv(vPar);
   return mPowderPattern_FullDeriv;
}

const CrystVector_REAL& PowderPattern::GetPowderPatternObs()const
{
   return mPowderPatternObs;
}

const CrystVector_REAL& PowderPattern::GetPowderPatternObsSigma()const
{
   return mPowderPatternObsSigma;
}

const CrystVector_REAL& PowderPattern::GetPowderPatternVariance()const
{
   return mPowderPatternVariance;
}

const CrystVector_REAL& PowderPattern::GetPowderPatternWeight()const
{
   return mPowderPatternWeight;
}

REAL PowderPattern::GetPowderPatternXMin()const
{
   if(mNbPoint==0) return 0;//:KLUDGE: ?
   if(true==mIsXAscending) return mX(0);
   return mX(mNbPoint-1);
}

REAL PowderPattern::GetPowderPatternXStep()const
{
   if(mNbPoint==0) return 0;//:KLUDGE: ?
   return abs((-mX(0)+mX(mNbPoint-1))/(mNbPoint-1));
}

REAL PowderPattern::GetPowderPatternXMax()const
{
   if(mNbPoint==0) return 0;//:KLUDGE: ?
   if(true==mIsXAscending) return mX(mNbPoint-1);
   return mX(0);
}
const CrystVector_REAL& PowderPattern::GetPowderPatternX()const
{
   return mX;
}

const CrystVector_REAL& PowderPattern::GetChi2Cumul(const int m)const
{
   VFN_DEBUG_ENTRY("PowderPattern::GetChi2Cumul()",3)
   mChi2Cumul.resize(mNbPoint);
   int mode = m;
   if((mode!=0) && (mode!=1)) mode = mOptProfileIntegration.GetChoice();
   if(0 == mode)
   {
      this->CalcPowderPatternIntegrated();
      if(mNbIntegrationUsed==0)
      	mChi2Cumul=0;
      else
      {
         const REAL *pObs=mIntegratedObs.data();
         const REAL *pCalc=mPowderPatternIntegratedCalc.data();
         const REAL *pWeight;
         if(mIntegratedWeight.numElements()==0) pWeight=mIntegratedWeightObs.data();
         else pWeight=mIntegratedWeight.data();

         REAL *pC2Cu=mChi2Cumul.data();
         for(int i=0;i<mIntegratedPatternMin(0);i++) *pC2Cu++ = 0;
         REAL chi2cumul=0,tmp;
         for(unsigned long j=1;j<mNbIntegrationUsed;j++)
         {
            tmp=(*pObs++ - *pCalc++) ;
            chi2cumul += *pWeight++ * tmp*tmp;
            for(int i=mIntegratedPatternMin(j-1);i<mIntegratedPatternMin(j);i++) *pC2Cu++ =chi2cumul;

            if(mIntegratedPatternMin(j)>(int)mNbPointUsed)
            {
               for(unsigned int i=mIntegratedPatternMin(j);i<mNbPoint;i++) *pC2Cu++ =chi2cumul;
               break;
            }
         }
         pC2Cu=mChi2Cumul.data()+mIntegratedPatternMin(mNbIntegrationUsed-1);
         for(unsigned int i=mIntegratedPatternMin(mNbIntegrationUsed-1);i<mNbPoint;i++) *pC2Cu++ =chi2cumul;
      }
   }
   else
   {
      this->CalcPowderPattern();
      const REAL *pObs=mPowderPatternObs.data();
      const REAL *pCalc=mPowderPatternCalc.data();
      const REAL *pWeight=mPowderPatternWeight.data();
      REAL *pC2Cu=mChi2Cumul.data();
      REAL chi2cumul=0,tmp;
      for(unsigned int i=0;i<mNbPointUsed;i++)
      {
         tmp = (*pObs++ - *pCalc++) ;
         chi2cumul += *pWeight++ * tmp*tmp;
         *pC2Cu++ = chi2cumul;
      }
   }
   VFN_DEBUG_EXIT("PowderPattern::GetChi2Cumul()",3)
   return mChi2Cumul;
}

const RefinableObjClock& PowderPattern::GetClockPowderPatternCalc()const
{  return mClockPowderPatternCalc;}

const RefinableObjClock& PowderPattern::GetClockPowderPatternPar()const
{  return mClockPowderPatternPar;}

const RefinableObjClock& PowderPattern::GetClockPowderPatternRadiation()const
{  return mClockPowderPatternRadiation;}

const RefinableObjClock& PowderPattern::GetClockPowderPatternXCorr()const
{  return mClockPowderPatternXCorr;}

void PowderPattern::SetXZero(const REAL newZero)
{
   mXZero=newZero;
   mClockPowderPatternPar.Click();
}

void PowderPattern::Set2ThetaDisplacement(const REAL displacement)
{
   m2ThetaDisplacement=displacement;
   mClockPowderPatternPar.Click();
}

void PowderPattern::Set2ThetaTransparency(const REAL transparency)
{
   m2ThetaTransparency=transparency;
   mClockPowderPatternPar.Click();
}

REAL PowderPattern::X2XCorr(const REAL x0)const
{
   REAL x=x0;
   if(  (mRadiation.GetWavelengthType()==WAVELENGTH_MONOCHROMATIC)
      ||(mRadiation.GetWavelengthType()==WAVELENGTH_ALPHA12))
      x += m2ThetaDisplacement*cos(x/2) +m2ThetaTransparency*sin(x);

   return x+mXZero;
}

REAL PowderPattern::X2PixelCorr(const REAL x0)const
{
   return this->X2Pixel(this->X2XCorr(x0));
}

REAL PowderPattern::X2Pixel(const REAL x)const
{
   //:TODO: faster if the step is actually constant.
   // Step may not be constant, so we guess twice before step-search
   REAL pixx;
   if(mIsXAscending==false)
   {
      VFN_DEBUG_MESSAGE("PowderPattern::X2Pixel()",1)
      long pix=(long)(mNbPoint-1-(x-this->GetPowderPatternXMin())/this->GetPowderPatternXStep());
      if((pix>0)&&(pix<((long)mNbPoint-1)))
      {
         // Why floor() and ceil() don't return a bloody integer is beyond me
         const REAL localStep=mX(pix)-mX(pix+1);
         if(localStep>0) pix -= (long)((x-mX(pix))/localStep);
      }
      VFN_DEBUG_MESSAGE("PowderPattern::X2Pixel():"<<x<<","<<pix,1)
      if(pix<1) pix=1;
      if(pix>((long)mNbPoint-2))pix=(long)mNbPoint-2;
      VFN_DEBUG_MESSAGE("PowderPattern::X2Pixel():"<<x<<","<<pix<<","<<mX(pix),1)
      if(mX(pix)<x)
      {
         for(;;pix--)
         {
            VFN_DEBUG_MESSAGE("PowderPattern::X2Pixel():"<<x<<","<<pix<<","<<mX(pix),1)
            if(mX(pix)>=x) break;
            if(pix==0) break;
         }
      }
      else
      {
         for(;;pix++)
         {
            if(mX(pix)<=x) {pix--;break;}
            if(pix==((long)mNbPoint-2)) break;
         }
      }
      // This assumes step is at least localy constant...
      VFN_DEBUG_MESSAGE("PowderPattern::X2Pixel():"<<x<<","<<pix<<","<<mX(pix),1)
      const REAL localStep=mX(pix)-mX(pix+1);
      pixx = (REAL)pix-(x-mX(pix))/localStep;
      VFN_DEBUG_MESSAGE("PowderPattern::X2Pixel():"<<x<<","<<pix<<","<<mX(pix),1)
   }
   else
   {
      VFN_DEBUG_MESSAGE("PowderPattern::X2Pixel():"<<x<<","<<this->GetPowderPatternXMin()<<","<<this->GetPowderPatternXMax(),1)
      long pix=(long)((x-this->GetPowderPatternXMin())/this->GetPowderPatternXStep());
      if((pix>0)&&(pix<((long)mNbPoint-1)))
      {
         // Why floor() and ceil() don't return a bloody integer is beyond me
         const REAL localStep=mX(pix+1)-mX(pix);
         if(localStep>0) pix += (long)((x-mX(pix))/localStep);
      }
      VFN_DEBUG_MESSAGE("PowderPattern::X2Pixel():"<<x<<","<<pix,1)
      if(pix<1) pix=1;
      if(pix>((long)mNbPoint-2))pix=(long)mNbPoint-2;
      VFN_DEBUG_MESSAGE("PowderPattern::X2Pixel():"<<x<<","<<pix<<","<<mX(pix),1)
      if(x<mX(pix))
      {
         for(;;pix--)
         {
            VFN_DEBUG_MESSAGE("PowderPattern::X2Pixel():"<<x<<","<<pix<<","<<mX(pix),1)
            if(mX(pix)<=x) break;
            if(pix==0) break;
         }
      }
      else
      {
         for(;;pix++)
         {
            VFN_DEBUG_MESSAGE("PowderPattern::X2Pixel():"<<x<<","<<pix<<","<<mX(pix),1)
            if(mX(pix)>=x) {pix-- ;break;}
            if(pix==((long)mNbPoint-2)) break;
         }
      }
      VFN_DEBUG_MESSAGE("PowderPattern::X2Pixel():"<<x<<","<<pix<<","<<mX(pix),1)
      if(pix>((long)mNbPoint-2))pix=(long)mNbPoint-2;
      // This assumes step is at least localy constant...
      const REAL localStep=mX(pix+1)-mX(pix);
      VFN_DEBUG_MESSAGE("PowderPattern::X2Pixel():"<<x<<","<<pix<<","<<mX(pix)<<","<<localStep,1)
      pixx = (REAL)pix+(x-mX(pix))/localStep;
   }
   VFN_DEBUG_MESSAGE("PowderPattern::X2Pixel():"<<x<<","<<pixx,1)
   return pixx;
}

void PowderPattern::ImportPowderPatternFullprof(const string &filename)
{
   //15.000   0.030  70.000 LANI4FE#1 REC 800 4JRS
   //2447.   2418.   2384.   2457.   2398.   2374.   2378.   2383.
   //...
   VFN_DEBUG_MESSAGE("PowderPattern::ImportPowderPatternFullprof() : \
from file : "+filename,5)
   ifstream fin(filename.c_str());
   if(!fin)
   {
      throw ObjCrystException("PowderPattern::ImportPowderPatternFullprof() : \
Error opening file for input:"+filename);
   }
   REAL min,max,step;
   fin >> min >> step >> max;
   min  *= DEG2RAD;
   max  *= DEG2RAD;
   step *= DEG2RAD;
   this->SetPowderPatternPar(min,step,(long)((max-min)/step+1.001));
   VFN_DEBUG_MESSAGE("PowderPattern::ImportPowderPatternFullprof() :"\
      << " 2Theta min=" << min*RAD2DEG << " 2Theta max=" << max*RAD2DEG \
      << " NbPoints=" << mNbPoint,5)
   mPowderPatternObs.resize (mNbPoint);
   mPowderPatternObsSigma.resize (mNbPoint);
   mPowderPatternWeight.resize(mNbPoint);

   char tmpComment[200];
   fin.getline(tmpComment,100);
   //if(""==mName) mName.append(tmpComment);

   for(unsigned long i=0;i<mNbPoint;i++) fin >> mPowderPatternObs(i);
   fin.close();
   this->SetSigmaToSqrtIobs();
   this->SetWeightToInvSigmaSq();
   mClockPowderPatternPar.Click();
   this->UpdateDisplay();
   {
      char buf [200];
      sprintf(buf,"Imported powder pattern: %d points, 2theta=%7.3f -> %7.3f, step=%6.3f",
              (int)mNbPoint,min*RAD2DEG,max*RAD2DEG,step*RAD2DEG);
      (*fpObjCrystInformUser)((string)buf);
   }
   VFN_DEBUG_MESSAGE("PowderPattern::ImportFullProfPattern():finished:"<<mNbPoint<<" points",5)
}

void PowderPattern::ImportPowderPatternPSI_DMC(const string &filename)
{
   VFN_DEBUG_MESSAGE("PowderPattern::ImportPowderPatternPSI_DMC() : \
from file : "+filename,5)
   ifstream fin(filename.c_str());
   if(!fin)
   {
      throw ObjCrystException("PowderPattern::ImportPowderPatternPSI_DMC() : \
Error opening file for input:"+filename);
   }
   //Skip the first two lines
      char tmpComment[200];
      fin.getline(tmpComment,190);
      fin.getline(tmpComment,190);
   REAL min,max,step;
   fin >> min >> step >> max;
   min  *= DEG2RAD;
   max  *= DEG2RAD;
   step *= DEG2RAD;
   this->SetPowderPatternPar(min,step,(unsigned long)((max-min)/step+1.001));
   VFN_DEBUG_MESSAGE("PowderPattern::ImportPowderPatternPSI_DMC() :"\
      << " 2Theta min=" << min*RAD2DEG << " 2Theta max=" << max*RAD2DEG \
      << " NbPoints=" << mNbPoint,5)
   mPowderPatternObs.resize (mNbPoint);
   mPowderPatternObsSigma.resize (mNbPoint);
   mPowderPatternWeight.resize(mNbPoint);

   fin.getline(tmpComment,100);
   //if(""==mName) mName.append(tmpComment);

   for(unsigned long i=0;i<mNbPoint;i++) fin >> mPowderPatternObs(i);
   for(unsigned long i=0;i<mNbPoint;i++) fin >> mPowderPatternObsSigma(i);
   fin.close();
   this->SetWeightToInvSigmaSq();
   mClockPowderPatternPar.Click();
   this->UpdateDisplay();
   {
      char buf [200];
      sprintf(buf,"Imported powder pattern: %d points, 2theta=%7.3f -> %7.3f, step=%6.3f",
              (int)mNbPoint,min*RAD2DEG,max*RAD2DEG,step*RAD2DEG);
      (*fpObjCrystInformUser)((string)buf);
   }
   VFN_DEBUG_MESSAGE("PowderPattern::ImportPowderPatternPSI_DMC():finished",5)
}

void PowderPattern::ImportPowderPatternILL_D1A5(const string &filename)
{
   VFN_DEBUG_MESSAGE("PowderPattern::ImportPowderPatternILL_D1AD2B() : \
from file : "+filename,5)
   ifstream fin(filename.c_str());
   if(!fin)
   {
      throw ObjCrystException("PowderPattern::ImportPowderPatternILL_D1AD2B() : \
Error opening file for input:"+filename);
   }
   //Skip the first three lines
      char tmpComment[200];
      fin.getline(tmpComment,190);
      fin.getline(tmpComment,190);
      fin.getline(tmpComment,190);

   fin >> mNbPoint;
   fin.getline(tmpComment,190);
   REAL min,step;
   fin >> min >> step;
   min  *= DEG2RAD;
   step *= DEG2RAD;
   this->SetPowderPatternPar(min,step,mNbPoint);
   VFN_DEBUG_MESSAGE("PowderPattern::ImportPowderPatternILL_D1AD2B() :"\
      << " 2Theta min=" << min*RAD2DEG << " 2Theta max=" << min*RAD2DEG+mNbPoint*step*RAD2DEG \
      << " NbPoints=" << mNbPoint,5)
   mPowderPatternObs.resize (mNbPoint);
   mPowderPatternObsSigma.resize (mNbPoint);
   mPowderPatternWeight.resize(mNbPoint);

   //if(""==mName) mName.append(tmpComment);

   for(unsigned long i=0;i<mNbPoint;i++) fin >> mPowderPatternObs(i);
   for(unsigned long i=0;i<mNbPoint;i++) fin >> mPowderPatternObsSigma(i);
   fin.close();
   this->SetWeightToInvSigmaSq();
   mClockPowderPatternPar.Click();
   this->UpdateDisplay();
   {
      char buf [200];
      sprintf(buf,"Imported powder pattern: %d points, 2theta=%7.3f -> %7.3f, step=%6.3f",
              (int)mNbPoint,min*RAD2DEG,(min+step*(mNbPoint-1))*RAD2DEG,
              step*RAD2DEG);
      (*fpObjCrystInformUser)((string)buf);
   }
   VFN_DEBUG_MESSAGE("PowderPattern::ImportPowderPatternILL_D1AD2B():finished",5)
}

void PowderPattern::ImportPowderPatternXdd(const string &filename)
{
   VFN_DEBUG_MESSAGE("PowderPattern::ImportPowderPatternXdd():from file :" \
                           +filename,5)
   ifstream fin (filename.c_str());
   if(!fin)
   {
      throw ObjCrystException("PowderPattern::ImportPowderPatternXdd() : \
Error opening file for input:"+filename);
   }
   char tmpComment[200];
   fin.getline(tmpComment,100);
   //if(""==mName) mName.append(tmpComment);
   REAL min,max,step,tmp;
   fin >> min >> step >> max;
   min  *= DEG2RAD;
   max *= DEG2RAD;
   step *= DEG2RAD;
   this->SetPowderPatternPar(min,step,(long)((max-min)/step+1.001));
   mPowderPatternObs.resize (mNbPoint);
   mPowderPatternObsSigma.resize(mNbPoint);
   mPowderPatternCalc.resize(mNbPoint);
   mPowderPatternWeight.resize(mNbPoint);
   mPowderPatternWeight=1.;

   fin >> tmp; //Count time
   fin >> tmp; //unused
   fin >> tmp; //unused (wavelength?)

   for(unsigned long i=0;i<mNbPoint;i++) fin >> mPowderPatternObs(i);
   fin.close();
   this->SetSigmaToSqrtIobs();
   this->SetWeightToInvSigmaSq();
   this->UpdateDisplay();
   {
      char buf [200];
      sprintf(buf,"Imported powder pattern: %d points, 2theta=%7.3f -> %7.3f, step=%6.3f",
              (int)mNbPoint,min*RAD2DEG,max*RAD2DEG,step*RAD2DEG);
      (*fpObjCrystInformUser)((string)buf);
   }
   VFN_DEBUG_MESSAGE("DiffractionDataPowder::ImportXddPattern() :finished",5)
}

void PowderPattern::ImportPowderPatternSietronicsCPI(const string &filename)
{
   VFN_DEBUG_ENTRY("PowderPattern::ImportPowderPatternSietronicsCPI():from file :" \
                           +filename,5)
   ifstream fin (filename.c_str());
   if(!fin)
   {
      throw ObjCrystException("PowderPattern::ImportPowderPatternSietronicsCPI() : \
Error opening file for input:"+filename);
   }
   char tmpComment[300];
   fin.getline(tmpComment,100);
   VFN_DEBUG_MESSAGE(" ->Discarded comment :"<<tmpComment,1)
   REAL min,max,step;
   fin >> min >> max >> step;
   min  *= DEG2RAD;
   max *= DEG2RAD;
   step *= DEG2RAD;
   this->SetPowderPatternPar(min,step,(long)((max-min)/step+1.001));
   mPowderPatternObs.resize (mNbPoint);
   mPowderPatternObsSigma.resize(mNbPoint);
   mPowderPatternCalc.resize(mNbPoint);
   mPowderPatternWeight.resize(mNbPoint);
   mPowderPatternWeight=1.;

   //Following lines are ignored (no fixed format ?)
   string str;
   do
   {
      fin>>str;
      VFN_DEBUG_MESSAGE(" ->Read :"<<str,1)
   } while ("SCANDATA"!=str);

   for(unsigned long i=0;i<mNbPoint;i++) fin >> mPowderPatternObs(i);
   fin.close();
   this->SetSigmaToSqrtIobs();
   this->SetWeightToInvSigmaSq();
   mClockPowderPatternPar.Click();
   this->UpdateDisplay();
   {
      char buf [200];
      sprintf(buf,"Imported powder pattern: %d points, 2theta=%7.3f -> %7.3f, step=%6.3f",
              (int)mNbPoint,min*RAD2DEG,max*RAD2DEG,step*RAD2DEG);
      (*fpObjCrystInformUser)((string)buf);
   }
   VFN_DEBUG_EXIT("DiffractionDataPowder::ImportPowderPatternSietronicsCPI()",5)
}

void PowderPattern::ImportPowderPattern2ThetaObsSigma(const string &filename,const int nbSkip)
{
   VFN_DEBUG_MESSAGE("DiffractionDataPowder::ImportPowderPattern2ThetaObsSigma():from:" \
                           +filename,5)
   ifstream fin (filename.c_str());
   if(!fin)
   {
      throw ObjCrystException("PowderPattern::ImportPowderPattern2ThetaObsSigma():\
Error opening file for input:"+filename);
   }
   {//Get rid of first lines
      char tmpComment[200];
      for(int i=0;i<nbSkip;i++) fin.getline(tmpComment,150);
   }
   mPowderPatternObs.resize (500);
   mPowderPatternObsSigma.resize(500);
   mX.resize(500);
   mNbPoint=0;

   do
   {
      fin >> mX              (mNbPoint);
      fin >> mPowderPatternObs     (mNbPoint);
      fin >> mPowderPatternObsSigma(mNbPoint);
      //cout << mX              (mNbPoint)<<" "
      //     << mPowderPatternObs     (mNbPoint)<<" "
      //    << mPowderPatternObsSigma(mNbPoint)<<endl;
      if(!fin) break;
      mNbPoint++;
      if( (mNbPoint%500)==0)
      {
         mX.resizeAndPreserve(mNbPoint+500);
         mPowderPatternObs.resizeAndPreserve(mNbPoint+500);
         mPowderPatternObsSigma.resizeAndPreserve(mNbPoint+500);
      }
   } while(fin.eof() == false);
   fin.close();

   mX.resizeAndPreserve                    (mNbPoint);
   mPowderPatternObs.resizeAndPreserve     (mNbPoint);
   mPowderPatternObsSigma.resizeAndPreserve(mNbPoint);
   mPowderPatternWeight.resize(mNbPoint);

   mX *= DEG2RAD;

   this->SetWeightToInvSigmaSq();
   mClockPowderPatternPar.Click();
   this->UpdateDisplay();
   {
      char buf [200];
      sprintf(buf,"Imported powder pattern: %d points, 2theta=%7.3f -> %7.3f, step=%6.3f",
              (int)mNbPoint,this->GetPowderPatternXMin()*RAD2DEG,
              this->GetPowderPatternXMax()*RAD2DEG,
              this->GetPowderPatternXStep()*RAD2DEG);
      (*fpObjCrystInformUser)((string)buf);
   }
   VFN_DEBUG_MESSAGE("PowderPattern::ImportPowderPattern2ThetaObsSigma()\
:finished: "<<mNbPoint<<" points",5)
}

void PowderPattern::ImportPowderPattern2ThetaObs(const string &filename,const int nbSkip)
{
   VFN_DEBUG_MESSAGE("PowderPattern::ImportPowderPattern2ThetaObs():from:" \
                           +filename,5)
   ifstream fin (filename.c_str());
   if(!fin)
   {
      throw ObjCrystException("PowderPattern::ImportPowderPattern2ThetaObs():\
Error opening file for input:"+filename);
   }
   {//Get rid of first lines
      char tmpComment[200];
      for(int i=0;i<nbSkip;i++) fin.getline(tmpComment,150);
   }
   mPowderPatternObs.resize (500);
   mPowderPatternObsSigma.resize(500);
   mX.resize(500);
   mNbPoint=0;

   do
   {
      fin >> mX(mNbPoint);
      fin >> mPowderPatternObs     (mNbPoint);
      mPowderPatternObsSigma(mNbPoint)
         =sqrt(mPowderPatternObs(mNbPoint));
      if(!fin) break;
      mNbPoint++;
      if( (mNbPoint%500)==0)
      {
         mX.resizeAndPreserve(mNbPoint+500);
         mPowderPatternObs.resizeAndPreserve(mNbPoint+500);
         mPowderPatternObsSigma.resizeAndPreserve(mNbPoint+500);
      }
   } while(fin.eof() == false);
   fin.close();

   mX.resizeAndPreserve              (mNbPoint);
   mPowderPatternObs.resizeAndPreserve     (mNbPoint);
   mPowderPatternObsSigma.resizeAndPreserve(mNbPoint);
   mPowderPatternWeight.resize(mNbPoint);

   mX *= DEG2RAD;

   this->SetSigmaToSqrtIobs();
   this->SetWeightToInvSigmaSq();
   mClockPowderPatternPar.Click();
   this->UpdateDisplay();
   {
      char buf [200];
      sprintf(buf,"Imported powder pattern: %d points, 2theta=%7.3f -> %7.3f, step=%6.3f",
              (int)mNbPoint,this->GetPowderPatternXMin()*RAD2DEG,
              this->GetPowderPatternXMax()*RAD2DEG,
              this->GetPowderPatternXStep()*RAD2DEG);
      (*fpObjCrystInformUser)((string)buf);
   }
   VFN_DEBUG_MESSAGE("PowderPattern::ImportPowderPattern2ThetaObs():finished",5)
}

void PowderPattern::ImportPowderPatternMultiDetectorLLBG42(const string &filename)
{

   //Sample 4: NaY + CF2=CCL2 T=20K, Lambda: 2.343 A.
   //     100       0   0.100      70       0       0
   //   3.000
   // 500000.  12000.    0.00    0.00
   //70  570369  569668  562868  532769  527469  495669  481969  452767  429468  4132
   //68  393269  372067  353769  337068  328268  310270  299469  296470  294668  2780
   //...
   //14  255814  282714  274014  281314  300314  302714  301114  298014  313214  3097
   //14  286914  295714  305214  300714  288311  295511  288511  3024 8  2937 7  2883
   // 7  2905 7  2895 7  2767 7  2777 7  2758 7  2495 7  2507 7  2496 7  2382 7  2329
   // 7  2542 7  2415 7  2049 7  2389 7  2270 7  2157 6  2227 6  2084 3  1875 2  2094
   // 1  1867
   //   -1000
   //  -10000
   VFN_DEBUG_MESSAGE("PowderPattern::ImportPowderPatternMultiDetectorLLBG42() : \
from file : "+filename,5)
   ifstream fin(filename.c_str());
   if(!fin)
   {
      throw ObjCrystException("PowderPattern::ImportPowderPatternMultiDetectorLLBG42() : \
Error opening file for input:"+filename);
   }

   string str;
   getline(fin,str);
   float junk;
   REAL min,step;
   fin >>junk>>junk>>step>>junk>>junk>>junk>>min>>junk>>junk>>junk>>junk;
   min  *= DEG2RAD;
   step *= DEG2RAD;
   VFN_DEBUG_MESSAGE("PowderPattern::ImportPowderPatternMultiDetectorLLBG42() :"\
      << " 2Theta min=" << min*RAD2DEG << " 2Theta step=" <<  step*RAD2DEG,5)
   mPowderPatternObs.resize (500);
   mPowderPatternObsSigma.resize (500);

   getline(fin,str);//finish reading line

   float tmp;
   string sub;
   float ct,iobs;
   mNbPoint=0;
   for(;;)
   {
      getline(fin,str);
      sscanf(str.c_str(),"%f",&tmp);
      if(tmp<0) break;
      const unsigned int nb=str.length()/8;
      for(unsigned int i=0;i<nb;i++)
      {
         if(mNbPoint==(unsigned int)mPowderPatternObs.numElements())
         {
            mPowderPatternObs.resizeAndPreserve(mNbPoint+500);
            mPowderPatternObsSigma.resizeAndPreserve(mNbPoint+500);
         }
         sub=str.substr(i*8,8);
         sscanf(sub.c_str(),"%2f%6f",&ct,&iobs);
         mPowderPatternObs(mNbPoint)=iobs;
         mPowderPatternObsSigma(mNbPoint++)=sqrt(iobs/ct);
      }
   }
   this->SetPowderPatternPar(min,step,mNbPoint);
   //exit(1);
   mPowderPatternObs.resizeAndPreserve (mNbPoint);
   mPowderPatternObsSigma.resizeAndPreserve (mNbPoint);
   mPowderPatternWeight.resizeAndPreserve(mNbPoint);

   fin.close();
   this->SetWeightToInvSigmaSq();
   mClockPowderPatternPar.Click();
   this->UpdateDisplay();
   {
      char buf [200];
      sprintf(buf,"Imported powder pattern: %d points, 2theta=%7.3f -> %7.3f, step=%6.3f",
              (int)mNbPoint,min*RAD2DEG,mX(mNbPoint-1)*RAD2DEG,step*RAD2DEG);
      (*fpObjCrystInformUser)((string)buf);
   }
   VFN_DEBUG_MESSAGE("PowderPattern::ImportPowderPatternMultiDetectorLLBG42():finished:"<<mNbPoint<<" points",5)
}

void PowderPattern::ImportPowderPatternFullprof4(const string &filename)
{
   //1.550   0.005  66.000
   // 213.135 193.243 208.811 185.873 231.607 200.995 196.792 187.516 215.977 199.634
   //  17.402  16.570  12.180  11.491  18.141  11.950  11.824  11.542  17.518  11.909
   // 211.890 185.740 204.610 200.645 199.489 169.549 203.189 178.298 186.241 198.522
   //  12.269  11.487  17.051  11.939  11.905  10.975  16.992  11.255  11.503  11.876
   VFN_DEBUG_MESSAGE("PowderPattern::ImportPowderPatternFullprof4() : \
from file : "+filename,5)
   ifstream fin(filename.c_str());
   if(!fin)
   {
      throw ObjCrystException("PowderPattern::ImportPowderPatternFullprof4() : \
Error opening file for input:"+filename);
   }
   REAL min,step,max;
   fin >> min >> step >> max;
   min *= DEG2RAD;
   max *= DEG2RAD;
   step *= DEG2RAD;
   this->SetPowderPatternPar(min,step,(long)((max-min)/step+1.001));
   VFN_DEBUG_MESSAGE("PowderPattern::ImportPowderPatternFullprof4() :"\
      << " 2Theta min=" << min*RAD2DEG << " 2Theta max=" << max*RAD2DEG \
      << " NbPoints=" << mNbPoint,5)
   mPowderPatternObs.resize (mNbPoint);
   mPowderPatternObsSigma.resize (mNbPoint);
   mPowderPatternWeight.resize(mNbPoint);

   string str;
   getline(fin,str);//read end of first line

   unsigned ct=0;
   unsigned ctSig=0;
   float line[10];
   for(;ct<mNbPoint;)
   {
      getline(fin,str);
      for(unsigned int i=0;i<str.size();i++) if(' '==str[i]) str[i]='0';
      sscanf(str.c_str(),"%8f%8f%8f%8f%8f%8f%8f%8f%8f%8f",
                         line+0,line+1,line+2,line+3,line+4,line+5,line+6,line+7,line+8,line+9);
      for(unsigned int j=0;j<10;j++)
         if(ct<mNbPoint) mPowderPatternObs(ct++)=line[j];
      getline(fin,str);
      for(unsigned int i=0;i<str.size();i++) if(' '==str[i]) str[i]='0';
      sscanf(str.c_str(),"%8f%8f%8f%8f%8f%8f%8f%8f%8f%8f",
                         line+0,line+1,line+2,line+3,line+4,line+5,line+6,line+7,line+8,line+9);
      for(unsigned int j=0;j<10;j++)
         if(ctSig<mNbPoint) mPowderPatternObsSigma(ctSig++)=line[j];
   }
   fin.close();
   this->SetWeightToInvSigmaSq();
   mClockPowderPatternPar.Click();
   this->UpdateDisplay();
   {
      char buf [200];
      sprintf(buf,"Imported powder pattern: %d points, 2theta=%7.3f -> %7.3f, step=%6.3f",
              (int)mNbPoint,min*RAD2DEG,max*RAD2DEG,
              step*RAD2DEG);
      (*fpObjCrystInformUser)((string)buf);
   }
   VFN_DEBUG_MESSAGE("PowderPattern::ImportFullProfPattern4():finished:"<<mNbPoint<<" points",5)
}

void PowderPattern::ImportPowderPatternTOF_ISIS_XYSigma(const string &filename)
{
   VFN_DEBUG_MESSAGE("DiffractionDataPowder::ImportPowderPatternTOF_ISIS_XYSigma():from:" \
                           +filename,5)
   ifstream fin (filename.c_str());
   if(!fin)
   {
      throw ObjCrystException("PowderPattern::ImportPowderPatternTOF_ISIS_XYSigma():\
Error opening file for input:"+filename);
   }
   {//Get rid of first line
      char junk[400];
      fin.getline(junk,150);
   }
   mPowderPatternObs.resize (500);
   mPowderPatternObsSigma.resize(500);
   mX.resize(500);
   mNbPoint=0;

   do
   {
      fin >> mX              (mNbPoint);
      fin >> mPowderPatternObs     (mNbPoint);
      fin >> mPowderPatternObsSigma(mNbPoint);
      //cout << mX              (mNbPoint)<<" "
      //     << mPowderPatternObs     (mNbPoint)<<" "
      //    << mPowderPatternObsSigma(mNbPoint)<<endl;
      if(!fin) break;
      mNbPoint++;
      if( (mNbPoint%500)==0)
      {
         mX.resizeAndPreserve(mNbPoint+500);
         mPowderPatternObs.resizeAndPreserve(mNbPoint+500);
         mPowderPatternObsSigma.resizeAndPreserve(mNbPoint+500);
      }
   } while(fin.eof() == false);
   fin.close();

   mX.resizeAndPreserve                    (mNbPoint);
   mPowderPatternObs.resizeAndPreserve     (mNbPoint);
   mPowderPatternObsSigma.resizeAndPreserve(mNbPoint);
   mPowderPatternWeight.resize(mNbPoint);

   // Reverse order of arrays, so that we are in ascending order of sin(theta)/lambda
   REAL tmp;
   for(unsigned long i=0;i<(mNbPoint/2);i++)
   {
      tmp=mX(i);
      mX(i)=mX(mNbPoint-1-i);
      mX(mNbPoint-1-i)=tmp;

      tmp=mPowderPatternObs(i);
      mPowderPatternObs(i)=mPowderPatternObs(mNbPoint-1-i);
      mPowderPatternObs(mNbPoint-1-i)=tmp;

      tmp=mPowderPatternObsSigma(i);
      mPowderPatternObsSigma(i)=mPowderPatternObsSigma(mNbPoint-1-i);
      mPowderPatternObsSigma(mNbPoint-1-i)=tmp;
   }
   this->SetPowderPatternX(mX);

   this->SetWeightToInvSigmaSq();
   this->SetRadiationType(RAD_NEUTRON);
   this->GetRadiation().SetWavelengthType(WAVELENGTH_TOF);
   mClockPowderPatternPar.Click();
   this->UpdateDisplay();
   {
      char buf [200];
      sprintf(buf,"Imported TOF powder pattern: %d points, TOF=%7.3f -> %7.3f",
              (int)mNbPoint,this->GetPowderPatternXMin(),
              this->GetPowderPatternXMax());
      (*fpObjCrystInformUser)((string)buf);
   }
   VFN_DEBUG_MESSAGE("PowderPattern::ImportPowderPatternTOF_ISIS_XYSigma()\
:finished: "<<mNbPoint<<" points",5)
}

void PowderPattern::ImportPowderPatternGSAS(const string &filename)
{
   VFN_DEBUG_ENTRY("PowderPattern::ImportPowderPatternGSAS():file:"<<filename,5)
   ifstream fin (filename.c_str());
   if(!fin)
   {
      throw ObjCrystException("PowderPattern::ImportPowderPatternGSAS():\
Error opening file for input:"+filename);
   }
   {//Get rid of title
      char title[81];
      fin.read(title,80);
      while(isprint(fin.peek())==false)
      {
         if(fin.eof()) break;
         fin.get();
      }
      cout<<"Title:"<<title<<endl;
      if(this->GetName()=="Change Me!") this->SetName(title);
   }

   //BANK  1    38101    7620 CONST     200.00      0.10 0 0 ESD
   int numBank,nbRecords;
   string binType, type;
   float bcoeff[4];
   //string line;
   char line[81];
   char bank[5];
   do
   {
      fin.getline(line,80);
      while(isprint(fin.peek())==false)
      {
         if(fin.eof()) break;
         fin.get();
      }
      sscanf(line,"%4s",bank);
      if(fin.eof())
         throw ObjCrystException("PowderPattern::ImportPowderPatternGSAS():\
Could not find BANK statement !! In file: "+filename);
   }
   while(string(bank)!=string("BANK"));

   {
      line[80]='\0';
      char binTypeC[20],typeC[20];
      sscanf(line,"%4s%d %ld %d %s %f %f %f %f %s",bank,&numBank,&mNbPoint,&nbRecords,
             binTypeC,&bcoeff[0],&bcoeff[1],&bcoeff[2],&bcoeff[3],typeC);
      binType=binTypeC;
      type=typeC;
   }
   if(binType=="CONST") binType="CONS";
   if((type!="ALT")&&(type!="ESD")) type="STD";

   cout<<"BANK #"<<numBank<<endl;
   cout<<"Number of data points:"<<mNbPoint<<endl;
   cout<<"Number of records:"<<nbRecords<<endl;
   cout<<"BinType:"<<binType<<endl;
   cout<<"BCoeff[1-4]:"<<bcoeff[0]<<","<<bcoeff[1]<<","<<bcoeff[2]<<","<<bcoeff[3]<<endl;
   cout<<"Type:"<<type<<endl;

   mPowderPatternObs.resize (mNbPoint);
   mPowderPatternObsSigma.resize(mNbPoint);
   mX.resize(mNbPoint);
   bool importOK=false;
   if((binType=="CONS") && (type=="ESD"))
   {
      this->SetPowderPatternPar(bcoeff[0]*DEG2RAD/100,bcoeff[1]*DEG2RAD/100,mNbPoint);
      string sub;
      unsigned long point=0;
      REAL iobs,isig;
      string substr;
      for(long i=0;i<nbRecords;i++)
      {
         fin.read(line,80);
         line[80]='\0';
         while(isprint(fin.peek())==false)
         {
            if(fin.eof()) break;
            fin.get();
         }
         for(unsigned int j=0;j<5;j++)
         {
            /*
            substr=string(line).substr(j*16,16);
            sscanf(substr.c_str(),"%8f%8f",&iobs,&isig);
            */
            substr=string(line).substr(j*16+0 ,8);
            istringstream(substr) >> iobs;
            substr=string(line).substr(j*16+8 ,8);
            istringstream(substr) >> isig;

            mPowderPatternObs(point)=iobs;
            mPowderPatternObsSigma(point++)=isig;
            if(point==mNbPoint) break;
         }
         if(point==mNbPoint) break;
      }
      importOK=true;
   }
   if((binType=="CONS") && (type=="STD"))
   {
      this->SetPowderPatternPar(bcoeff[0]*DEG2RAD/100,bcoeff[1]*DEG2RAD/100,mNbPoint);
      unsigned long point=0;
      REAL iobs;
      int nc;
      string substr;
      for(long i=0;i<nbRecords;i++)
      {
         fin.read(line,80);
         line[80]='\0';
         while(isprint(fin.peek())==false)
         {
            if(fin.eof()) break;
            fin.get();
         }
         for(unsigned int j=0;j<10;j++)
         {
            /*
            substr=string(line).substr(j*8,8);
            if(substr.substr(0,2)==string("  "))
            {
               nc=1;
               sscanf(substr.c_str(),"%8f",&iobs);
            }
            else sscanf(substr.c_str(),"%2d%6f",&nc,&iobs);
            */
            substr=string(line).substr(j*8+0 ,2);
            if(substr=="  ") nc=1;
            else sscanf(substr.c_str(),"%d",&nc);
            substr=string(line).substr(j*8+2 ,6);
            istringstream(substr) >> iobs;
            mPowderPatternObs(point)=iobs;
            mPowderPatternObsSigma(point++)=sqrt(iobs)/sqrt((REAL)nc);
            if(point==mNbPoint) break;
         }
         if(point==mNbPoint) break;
      }
      importOK=true;
   }
   if((binType=="RALF") && (type=="ALT"))
   {
      this->SetRadiationType(RAD_NEUTRON);
      this->GetRadiation().SetWavelengthType(WAVELENGTH_TOF);
      mClockPowderPatternPar.Click();

      unsigned long point=0;
      REAL x,iobs,iobssigma;
      string substr;
      for(long i=0;i<nbRecords;i++)
      {
         fin.read(line,80);
         line[80]='\0';
         while(isprint(fin.peek())==false)
         {
            if(fin.eof()) break;
            fin.get();
         }
         for(unsigned int j=0;j<4;j++)
         {//4 records per line
            /* Does not work because sscanf ignores the leading spaces and shifts the reading !
            substr=string(line).substr(j*20,20);
            sscanf(substr.c_str(),"%8f%7f%5f",&x,&iobs,&iobssigma);
            */
            substr=string(line).substr(j*20+0 ,8);
            istringstream(substr) >> x;
            substr=string(line).substr(j*20+8 ,7);
            istringstream(substr) >> iobs;
            substr=string(line).substr(j*20+15,5);
            istringstream(substr) >> iobssigma;
            mPowderPatternObs(point)=iobs;
            mPowderPatternObsSigma(point)=iobssigma;
            mX(point)=x/32;
            if(++point==mNbPoint) break;
         }
         if(point==mNbPoint) break;
      }
      // Reverse order of arrays, so that we are in ascending order of sin(theta)/lambda
      REAL tmp;
      for(unsigned long i=0;i<(mNbPoint/2);i++)
      {
         tmp=mX(i);
         mX(i)=mX(mNbPoint-1-i);
         mX(mNbPoint-1-i)=tmp;

         tmp=mPowderPatternObs(i);
         mPowderPatternObs(i)=mPowderPatternObs(mNbPoint-1-i);
         mPowderPatternObs(mNbPoint-1-i)=tmp;

         tmp=mPowderPatternObsSigma(i);
         mPowderPatternObsSigma(i)=mPowderPatternObsSigma(mNbPoint-1-i);
         mPowderPatternObsSigma(mNbPoint-1-i)=tmp;
      }
      importOK=true;
   }
   fin.close();
   if(!importOK)
   {
      mNbPoint=0;
      mPowderPatternObs.resize (mNbPoint);
      mPowderPatternObsSigma.resize(mNbPoint);
      mX.resize(mNbPoint);
      throw ObjCrystException("PowderPattern::ImportPowderPatternGSAS(): Sorry, \
this type of format is not handled yet (send an example file to the Fox author)!:"+filename);
   }
   mPowderPatternWeight.resize(mNbPoint);
   this->SetPowderPatternX(mX);
   this->SetWeightToInvSigmaSq();

   this->UpdateDisplay();
   if(this->GetRadiation().GetWavelengthType()==WAVELENGTH_TOF)
   {
      char buf [200];
      sprintf(buf,"Imported powder pattern: %d points, tof=%7.3f us-> %7.3f us",
              (int)mNbPoint,this->GetPowderPatternXMin(),
              this->GetPowderPatternXMax());
      (*fpObjCrystInformUser)((string)buf);
   }
   else
   {
      char buf [200];
      sprintf(buf,"Imported powder pattern: %d points, 2theta=%7.3f -> %7.3f, step=%6.3f",
              (int)mNbPoint,this->GetPowderPatternXMin()*RAD2DEG,
              this->GetPowderPatternXMax()*RAD2DEG,
              this->GetPowderPatternXStep()*RAD2DEG);
      (*fpObjCrystInformUser)((string)buf);
   }
   VFN_DEBUG_EXIT("PowderPattern::ImportPowderPatternGSAS():file:"<<filename,5)
}

void PowderPattern::ImportPowderPatternCIF(const CIF &cif)
{
   VFN_DEBUG_ENTRY("PowderPattern::ImportPowderPatternCIF():file:",5)
   for(map<string,CIFData>::const_iterator pos=cif.mvData.begin();pos!=cif.mvData.end();++pos)
      if(pos->second.mPowderPatternObs.size()>10)
      {
         mNbPoint=pos->second.mPowderPatternObs.size();
         mX.resize(mNbPoint);
         mPowderPatternObs.resize(mNbPoint);
         mPowderPatternObsSigma.resize(mNbPoint);
         mPowderPatternWeight.resize(mNbPoint);
         if(pos->second.mDataType==WAVELENGTH_TOF)
         {
            this->SetRadiationType(RAD_NEUTRON);
            this->GetRadiation().SetWavelengthType(WAVELENGTH_TOF);
            mClockPowderPatternPar.Click();
         }
         else this->GetRadiation().SetWavelength(pos->second.mWavelength);
         for(unsigned long i=0;i<mNbPoint;++i)
         {
            mPowderPatternObs(i)=pos->second.mPowderPatternObs[i];
            mX(i)=pos->second.mPowderPatternX[i];
            mPowderPatternObsSigma(i)=pos->second.mPowderPatternSigma[i];
         }
         this->SetWeightToInvSigmaSq();
         this->SetPowderPatternX(mX);
      }
   VFN_DEBUG_EXIT("PowderPattern::ImportPowderPatternCIF():file:",5)
}


void PowderPattern::SetPowderPatternObs(const CrystVector_REAL& obs)
{
   VFN_DEBUG_MESSAGE("PowderPattern::ImportPowderPatternObs()",5)
   if((unsigned long)obs.numElements() != mNbPoint)
   {
      cout << obs.numElements()<<" "<<mNbPoint<<" "<<this<<endl;
      throw(ObjCrystException("PowderPattern::SetPowderPatternObs(vect): The \
supplied vector of observed intensities does not have the expected number of points!"));
   }

   mPowderPatternObs=obs;
   mPowderPatternObsSigma.resize(mPowderPatternObs.numElements());
   mPowderPatternWeight.resize(mNbPoint);

   this->SetSigmaToSqrtIobs();
   this->SetWeightToInvSigmaSq();
   mClockIntegratedFactorsPrep.Reset();
   {
      char buf [200];
      sprintf(buf,"Changed powder pattern: %d points, 2theta=%7.3f -> %7.3f, step=%6.3f",
              (int)mNbPoint,mX(0)*RAD2DEG,mX(mNbPoint-1)*RAD2DEG,
              (mX(mNbPoint-1)-mX(0))/(mNbPoint-1)*RAD2DEG);
      (*fpObjCrystInformUser)((string)buf);
   }
}
void PowderPattern::SavePowderPattern(const string &filename) const
{
   VFN_DEBUG_MESSAGE("PowderPattern::SavePowderPattern",5)
   this->CalcPowderPattern();
   ofstream out(filename.c_str());
   CrystVector_REAL ttheta;
   ttheta=mX;
   if(this->GetRadiation().GetWavelengthType()!=WAVELENGTH_TOF) ttheta *= RAD2DEG;

   CrystVector_REAL diff;
   diff=mPowderPatternObs;
   diff-=mPowderPatternCalc;
   out << "#    2Theta/TOF Iobs       ICalc   Iobs-Icalc    Weight  Comp0" << endl;
   out << FormatVertVector<REAL>(ttheta,
                    mPowderPatternObs,
                    mPowderPatternCalc,
                    diff,mPowderPatternWeight,
                    mPowderPatternComponentRegistry.GetObj(0).mPowderPatternCalc,16,8);
   out.close();
   VFN_DEBUG_MESSAGE("DiffractionDataPowder::SavePowderPattern:End",3)
}

void PowderPattern::PrintObsCalcData(ostream&os)const
{
   VFN_DEBUG_MESSAGE("DiffractionDataPowder::PrintObsCalcData()",5);
   CrystVector_REAL ttheta;
   ttheta=mX;
   if(this->GetRadiation().GetWavelengthType()!=WAVELENGTH_TOF) ttheta *= RAD2DEG;
   os << "PowderPattern : " << mName <<endl;
   os << "      2Theta/TOF  Obs          Sigma        Calc        Weight" <<endl;
   os << FormatVertVector<REAL>(ttheta,mPowderPatternObs,mPowderPatternObsSigma,
               mPowderPatternCalc,mPowderPatternWeight,16,8);
   //            mPowderPatternComponentRegistry.GetObj(0).mPowderPatternCalc,12,4);
}

REAL PowderPattern::GetR()const
{
   if(  (0==this->GetPowderPatternObs().numElements())
      ||(0==GetNbPowderPatternComponent()))
   {
      return 0;
   }
   this->CalcPowderPattern();
   TAU_PROFILE("PowderPattern::GetR()","void ()",TAU_DEFAULT);

   REAL tmp1=0.;
   REAL tmp2=0.;

   unsigned long maxPoints=mNbPointUsed;
   if(  (true==mStatisticsExcludeBackground)
      &&(mPowderPatternBackgroundCalc.numElements()>0))
   {
      const REAL *p1, *p2, *p3;
      p1=mPowderPatternCalc.data();
      p2=mPowderPatternObs.data();
      p3=mPowderPatternBackgroundCalc.data();
      const long nbExclude=mExcludedRegionMinX.numElements();
      if(0==nbExclude)
      {
         VFN_DEBUG_MESSAGE("PowderPattern::GetR():Exclude Backgd",4);
         for(unsigned long i=0;i<maxPoints;i++)
         {
            tmp1 += ((*p1)-(*p2)) * ((*p1)-(*p2));
            tmp2 += ((*p2)-(*p3)) * ((*p2)-(*p3));
            p1++;p2++;p3++;
         }
      }
      else
      {
         VFN_DEBUG_MESSAGE("PowderPattern::GetR():Exclude Backgd,Exclude regions",4);
         unsigned long min,max;
         unsigned long i=0;
         for(int j=0;j<nbExclude;j++)
         {
            min=(unsigned long)floor(this->X2Pixel(mExcludedRegionMinX(j)));
            max=(unsigned long)ceil (this->X2Pixel(mExcludedRegionMaxX(j)));
            if(min>maxPoints) break;
            if(max>maxPoints)max=maxPoints;
            for(;i<min;i++)//! min is the *beginning* of the excluded region !
            {
               tmp1 += ((*p1)-(*p2)) * ((*p1)-(*p2));
               tmp2 += ((*p2)-(*p3)) * ((*p2)-(*p3));
               p1++;p2++;p3++;
            }
            p1 += max-i;
            p2 += max-i;
            p3 += max-i;
            i  += max-i;
         }
         for(;i<maxPoints;i++)
         {
            tmp1 += ((*p1)-(*p2)) * ((*p1)-(*p2));
            tmp2 += ((*p2)-(*p3)) * ((*p2)-(*p3));
            p1++;p2++;p3++;
         }

      }
   } // Exclude Background ?
   else
   {
      const REAL *p1, *p2;
      p1=mPowderPatternCalc.data();
      p2=mPowderPatternObs.data();
      const long nbExclude=mExcludedRegionMinX.numElements();
      if(0==nbExclude)
      {
         VFN_DEBUG_MESSAGE("PowderPattern::GetR()",4);
         for(unsigned long i=0;i<maxPoints;i++)
         {
            tmp1 += ((*p1)-(*p2))*((*p1)-(*p2));
            tmp2 += (*p2) * (*p2);
            //cout <<i<<":"<< tmp1 << " "<<tmp2 << " " << *p1 <<" "<<*p2<<endl;
            p1++;p2++;
         }
      }
      else
      {
         VFN_DEBUG_MESSAGE("PowderPattern::GetR(),Exclude regions",4);
         unsigned long min,max;
         unsigned long i=0;
         for(int j=0;j<nbExclude;j++)
         {
            min=(unsigned long)floor(this->X2Pixel(mExcludedRegionMinX(j)));
            max=(unsigned long)ceil (this->X2Pixel(mExcludedRegionMaxX(j)));
            if(min>maxPoints) break;
            if(max>maxPoints)max=maxPoints;
            for(;i<min;i++)//! min is the *beginning* of the excluded region !
            {
               tmp1 += ((*p1)-(*p2))*((*p1)-(*p2));
               tmp2 += (*p2) * (*p2);
               p1++;p2++;
            }
            p1 += max-i;
            p2 += max-i;
            i  += max-i;
         }
         for(;i<maxPoints;i++)
         {
            tmp1 += ((*p1)-(*p2))*((*p1)-(*p2));
            tmp2 += (*p2) * (*p2);
            p1++;p2++;
         }
      }
   }

   VFN_DEBUG_MESSAGE("PowderPattern::GetR()="<<sqrt(tmp1/tmp2),4);
   //cout << FormatVertVector<REAL>(mPowderPatternCalc,mPowderPatternObs);
   //this->SavePowderPattern("refinedPattern.out");
   //abort();
   return sqrt(tmp1/tmp2);
}

REAL PowderPattern::GetIntegratedR()const
{
   if(  (0==this->GetPowderPatternObs().numElements())
      ||(0==GetNbPowderPatternComponent()))
   {
      return 0;
   }
   this->CalcPowderPattern();
   this->PrepareIntegratedRfactor();
   VFN_DEBUG_ENTRY("PowderPattern::GetIntegratedR()",4);
   TAU_PROFILE("PowderPattern::GetIntegratedR()","void ()",TAU_DEFAULT);

   REAL tmp1=0.;
   REAL tmp2=0.;
   const long numInterval=mIntegratedPatternMin.numElements();
   if(  (true==mStatisticsExcludeBackground)
      &&(mPowderPatternBackgroundCalc.numElements()>0))
   {
      const REAL *p1, *p2, *p3;
      CrystVector_REAL integratedCalc(numInterval);
      integratedCalc=0;
      CrystVector_REAL backgdCalc(numInterval);
      backgdCalc=0;
      REAL *pp1=integratedCalc.data();
      REAL *pp2=backgdCalc.data();
      for(int i=0;i<numInterval;i++)
      {
         const long max=mIntegratedPatternMax(i);
         p1=mPowderPatternCalc.data()+mIntegratedPatternMin(i);
         for(int j=mIntegratedPatternMin(i);j<=max;j++) *pp1 += *p1++;
         pp1++;
         p1=mPowderPatternBackgroundCalc.data()+mIntegratedPatternMin(i);
         for(int j=mIntegratedPatternMin(i);j<=max;j++) *pp2 += *p1++;
         pp2++;
      }

      p1=integratedCalc.data();
      p2=mIntegratedObs.data();
      p3=backgdCalc.data();
      VFN_DEBUG_MESSAGE("PowderPattern::GetIntegratedR():Exclude Backgd",2);
      for(long i=0;i<numInterval;i++)
      {
         tmp1 += ((*p1)-(*p2)) * ((*p1)-(*p2));
         tmp2 += ((*p2)-(*p3)) * ((*p2)-(*p3));
         p1++;p2++;p3++;
      }
   } // Exclude Background ?
   else
   {
      const REAL *p1, *p2;
      CrystVector_REAL integratedCalc(numInterval);
      integratedCalc=0;
      REAL *pp1=integratedCalc.data();
      for(int i=0;i<numInterval;i++)
      {
         const long max=mIntegratedPatternMax(i);
         p1=mPowderPatternCalc.data()+mIntegratedPatternMin(i);
         for(int j=mIntegratedPatternMin(i);j<=max;j++) *pp1 += *p1++;
         pp1++;
      }
      p1=integratedCalc.data();
      p2=mIntegratedObs.data();
      VFN_DEBUG_MESSAGE("PowderPattern::GetIntegratedR()",2);
      for(long i=0;i<numInterval;i++)
      {
         tmp1 += ((*p1)-(*p2))*((*p1)-(*p2));
         tmp2 += (*p2) * (*p2);
         //cout <<i<<":"<< tmp1 << " "<<tmp2 << " " << *p1 <<" "<<*p2<<endl;
         p1++;p2++;
      }
   }

   VFN_DEBUG_EXIT("PowderPattern::GetIntegratedR()="<<sqrt(tmp1/tmp2),4);
   return sqrt(tmp1/tmp2);
}

REAL PowderPattern::GetRw()const
{
   if(  (0==this->GetPowderPatternObs().numElements())
      ||(0==GetNbPowderPatternComponent()))
   {
      return 0;
   }
   this->CalcPowderPattern();
   TAU_PROFILE("PowderPattern::GetRw()","void ()",TAU_DEFAULT);
   VFN_DEBUG_MESSAGE("PowderPattern::GetRw()",3);


   //cout <<FormatVertVector<REAL>(mPowderPatternObs,
   //                               mPowderPatternCalc,
   //                               mPowderPatternWeight);
   REAL tmp1=0.;
   REAL tmp2=0.;

   unsigned long maxPoints=mNbPointUsed;

   if(  (true==mStatisticsExcludeBackground)
      &&(mPowderPatternBackgroundCalc.numElements()>0))
   {
      VFN_DEBUG_MESSAGE("PowderPattern::GetRw():Exclude Backgd",3);
      const REAL *p1, *p2, *p3, *p4;
      p1=mPowderPatternCalc.data();
      p2=mPowderPatternObs.data();
      p3=mPowderPatternBackgroundCalc.data();
      p4=mPowderPatternWeight.data();
      const long nbExclude=mExcludedRegionMinX.numElements();
      if(0==nbExclude)
      {
         for(unsigned long i=0;i<maxPoints;i++)
         {
            tmp1 += *p4   * ((*p1)-(*p2)) * ((*p1)-(*p2));
            tmp2 += *p4++ * ((*p2)-(*p3)) * ((*p2)-(*p3));
            p1++;p2++;p3++;
         }
      }
      else
      {
         unsigned long min,max;
         unsigned long i=0;
         for(int j=0;j<nbExclude;j++)
         {
            min=(unsigned long)floor(this->X2Pixel(mExcludedRegionMinX(j)));
            max=(unsigned long)ceil (this->X2Pixel(mExcludedRegionMaxX(j)));
            if(min>maxPoints) break;
            if(max>maxPoints)max=maxPoints;
            for(;i<min;i++)//! min is the *beginning* of the excluded region !
            {
               tmp1 += *p4   * ((*p1)-(*p2)) * ((*p1)-(*p2));
               tmp2 += *p4++ * ((*p2)-(*p3)) * ((*p2)-(*p3));
               p1++;p2++;p3++;
            }
            p1 += max-i;
            p2 += max-i;
            p3 += max-i;
            p4 += max-i;
            i  += max-i;
         }
         for(;i<maxPoints;i++)
         {
            tmp1 += *p4   * ((*p1)-(*p2)) * ((*p1)-(*p2));
            tmp2 += *p4++ * ((*p2)-(*p3)) * ((*p2)-(*p3));
            p1++;p2++;p3++;
         }

      }
   }
   else
   {
      VFN_DEBUG_MESSAGE("PowderPattern::GetRw()",3);
      const REAL *p1, *p2, *p4;
      p1=mPowderPatternCalc.data();
      p2=mPowderPatternObs.data();
      p4=mPowderPatternWeight.data();
      const long nbExclude=mExcludedRegionMinX.numElements();
      if(0==nbExclude)
      {
         for(unsigned long i=0;i<maxPoints;i++)
         {
            tmp1 += *p4   * ((*p1)-(*p2))*((*p1)-(*p2));
            tmp2 += *p4++ * (*p2) * (*p2);
            p1++;p2++;
         }
      }
      else
      {
         unsigned long min,max;
         unsigned long i=0;
         for(int j=0;j<nbExclude;j++)
         {
            min=(unsigned long)floor(this->X2Pixel(mExcludedRegionMinX(j)));
            max=(unsigned long)ceil (this->X2Pixel(mExcludedRegionMaxX(j)));
            if(min>maxPoints) break;
            if(max>maxPoints)max=maxPoints;
            for(;i<min;i++)//! min is the *beginning* of the excluded region !
            {
               tmp1 += *p4   * ((*p1)-(*p2))*((*p1)-(*p2));
               tmp2 += *p4++ * (*p2) * (*p2);
               p1++;p2++;
            }
            p1 += max-i;
            p2 += max-i;
            p4 += max-i;
            i  += max-i;
         }
         for(;i<maxPoints;i++)
         {
            tmp1 += *p4   * ((*p1)-(*p2))*((*p1)-(*p2));
            tmp2 += *p4++ * (*p2) * (*p2);
            p1++;p2++;
         }
      }
   }
   VFN_DEBUG_MESSAGE("PowderPattern::GetRw()="<<sqrt(tmp1/tmp2),3);
   return sqrt(tmp1/tmp2);
}
REAL PowderPattern::GetIntegratedRw()const
{
   if(  (0==this->GetPowderPatternObs().numElements())
      ||(0==GetNbPowderPatternComponent()))
   {
      return 0;
   }
   this->CalcPowderPattern();
   this->PrepareIntegratedRfactor();
   TAU_PROFILE("PowderPattern::GetIntegratedRw()","void ()",TAU_DEFAULT);

   REAL tmp1=0.;
   REAL tmp2=0.;
   const long numInterval=mIntegratedPatternMin.numElements();
   if(  (true==mStatisticsExcludeBackground)
      &&(mPowderPatternBackgroundCalc.numElements()>0))
   {
      const REAL *p1, *p2, *p3, *p4;
      CrystVector_REAL integratedCalc(numInterval);
      integratedCalc=0;
      CrystVector_REAL backgdCalc(numInterval);
      backgdCalc=0;
      REAL *pp1=integratedCalc.data();
      REAL *pp2=backgdCalc.data();
      for(int i=0;i<numInterval;i++)
      {
         const long max=mIntegratedPatternMax(i);
         p1=mPowderPatternCalc.data()+mIntegratedPatternMin(i);
         for(int j=mIntegratedPatternMin(i);j<=max;j++) *pp1 += *p1++;
         pp1++;
         p1=mPowderPatternBackgroundCalc.data()+mIntegratedPatternMin(i);
         for(int j=mIntegratedPatternMin(i);j<=max;j++) *pp2 += *p1++;
         pp2++;
      }

      p1=integratedCalc.data();
      p2=mIntegratedObs.data();
      p3=backgdCalc.data();
      if(mIntegratedWeight.numElements()==0) p4=mIntegratedWeightObs.data();
      else p4=mIntegratedWeight.data();
      VFN_DEBUG_MESSAGE("PowderPattern::GetIntegratedRw():Exclude Backgd",4);
      for(long i=0;i<numInterval;i++)
      {
         tmp1 += *p4   * ((*p1)-(*p2)) * ((*p1)-(*p2));
         tmp2 += *p4++ * ((*p2)-(*p3)) * ((*p2)-(*p3));
         //cout <<i<<": " <<mIntegratedPatternMin(i)<<"->"<<mIntegratedPatternMax(i)
         //     <<" "<< tmp1 << " "<<tmp2 << " " << *p1 <<" "<<*p2<<" "<<*p3<<" "<<*(p4-1) <<endl;
         p1++;p2++;p3++;
      }
   } // Exclude Background ?
   else
   {
      const REAL *p1, *p2, *p4;
      CrystVector_REAL integratedCalc(numInterval);
      integratedCalc=0;
      REAL *pp1=integratedCalc.data();
      for(int i=0;i<numInterval;i++)
      {
         const long max=mIntegratedPatternMax(i);
         p1=mPowderPatternCalc.data()+mIntegratedPatternMin(i);
         for(int j=mIntegratedPatternMin(i);j<=max;j++) *pp1 += *p1++;
         pp1++;
      }
      p1=integratedCalc.data();
      p2=mIntegratedObs.data();
      if(mIntegratedWeight.numElements()==0) p4=mIntegratedWeightObs.data();
      else p4=mIntegratedWeight.data();
      VFN_DEBUG_MESSAGE("PowderPattern::GetIntegratedRw()",4);
      for(long i=0;i<numInterval;i++)
      {
         tmp1 += *p4   * ((*p1)-(*p2))*((*p1)-(*p2));
         tmp2 += *p4++ * (*p2) * (*p2);
         //cout <<i<<":"<< tmp1 << " "<<tmp2 << " " << *p1 <<" "<<*p2<<endl;
         p1++;p2++;
      }
   }

   VFN_DEBUG_MESSAGE("PowderPattern::GetIntegratedRw()="<<sqrt(tmp1/tmp2),4);
   //cout << FormatVertVector<REAL>(mPowderPatternCalc,mPowderPatternObs);
   //this->SavePowderPattern("refinedPattern.out");
   //abort();
   return sqrt(tmp1/tmp2);
}

REAL PowderPattern::GetChi2()const
{
   if(  (0==this->GetPowderPatternObs().numElements())
      ||(0==GetNbPowderPatternComponent()))
   {
      mChi2=0.;
      return mChi2;
   }
   this->CalcNbPointUsed();
   if(mClockChi2>mClockMaster) return mChi2;

   this->CalcPowderPattern();
   if(  (mClockChi2>mClockPowderPatternPar)
      &&(mClockChi2>mClockScaleFactor)
      &&(mClockChi2>mClockPowderPatternCalc)) return mChi2;
   // We want the best scale factor
   this->FitScaleFactorForRw();

   TAU_PROFILE("PowderPattern::GetChi2()","void ()",TAU_DEFAULT);

   VFN_DEBUG_ENTRY("PowderPattern::GetChi2()",3);

   const unsigned long maxPoints=mNbPointUsed;

   mChi2=0.;
   mChi2LikeNorm=0.;
   VFN_DEBUG_MESSAGE("PowderPattern::GetChi2()Integrated profiles",3);
   const REAL * RESTRICT p1, * RESTRICT p2, * RESTRICT p3;
   p1=mPowderPatternCalc.data();
   p2=mPowderPatternObs.data();
   p3=mPowderPatternWeight.data();
   const long nbExclude=mExcludedRegionMinX.numElements();
   if(0==nbExclude)
   {
      for(unsigned long i=0;i<maxPoints;i++)
      {
         mChi2 += *p3 * ((*p1)-(*p2))*((*p1)-(*p2));
         if(*p3<=0) p3++;
         else mChi2LikeNorm -= log(*p3++);
         p1++;p2++;
      }
   }
   else
   {
      unsigned long min,max;
      unsigned long i=0;
      for(int j=0;j<nbExclude;j++)
      {
         min=(unsigned long)floor(this->X2Pixel(mExcludedRegionMinX(j)));
         max=(unsigned long)ceil (this->X2Pixel(mExcludedRegionMaxX(j)));
         if(min>maxPoints) break;
         if(max>maxPoints)max=maxPoints;
         for(;i<min;i++)//! min is the *beginning* of the excluded region !
         {
            mChi2 += *p3 * ((*p1)-(*p2))*((*p1)-(*p2));
            if(*p3<=0) p3++;
            else mChi2LikeNorm -= log(*p3++);
            p1++;p2++;
         }
         p1 += max-i;
         p2 += max-i;
         p3 += max-i;
         i  += max-i;
      }
      for(;i<maxPoints;i++)
      {
         mChi2 += *p3 * ((*p1)-(*p2))*((*p1)-(*p2));
         if(*p3<=0) p3++;
         else mChi2LikeNorm -= log(*p3++);
         p1++;p2++;
      }
   }
   mChi2LikeNorm/=2;
   VFN_DEBUG_MESSAGE("Chi^2="<<mChi2<<", log(norm)="<<mChi2LikeNorm,3)
   mClockChi2.Click();
   VFN_DEBUG_EXIT("PowderPattern::GetChi2()="<<mChi2,3);
   return mChi2;
}

REAL PowderPattern::GetIntegratedChi2()const
{
   if(  (0==this->GetPowderPatternObs().numElements())
      ||(0==GetNbPowderPatternComponent()))
   {
      mIntegratedChi2=0.;
      return mIntegratedChi2;
   }
   this->CalcNbPointUsed();
   if(mClockIntegratedChi2>mClockMaster) return mIntegratedChi2;

   this->CalcPowderPatternIntegrated();
   if(  (mClockChi2>mClockPowderPatternPar)
      &&(mClockChi2>mClockScaleFactor)
      &&(mClockChi2>mClockPowderPatternIntegratedCalc)) return mIntegratedChi2;

   // We want the best scale factor
   this->FitScaleFactorForIntegratedRw();

   TAU_PROFILE("PowderPattern::GetIntegratedChi2()","void ()",TAU_DEFAULT);

   VFN_DEBUG_ENTRY("PowderPattern::GetIntegratedChi2()",3);

   mIntegratedChi2=0.;
   mIntegratedChi2LikeNorm=0.;
   VFN_DEBUG_MESSAGE("PowderPattern::GetIntegratedChi2()",3);
   const REAL * RESTRICT p1, * RESTRICT p2, * RESTRICT p3;
   p1=mPowderPatternIntegratedCalc.data();
   p2=mIntegratedObs.data();
   if(mIntegratedWeight.numElements()==0) p3=mIntegratedWeightObs.data();
   else p3=mIntegratedWeight.data();
   double weightProd=1.;
   VFN_DEBUG_MESSAGE("PowderPattern::GetIntegratedIntegratedRw()",4);
   for(unsigned long i=0;i<mNbIntegrationUsed;)
   {
      // group weights to avoid computing too many log()
      // group only a limited number to avoid underflow...
      for(unsigned long j=0;j<32;++j)
      {
         mIntegratedChi2 += *p3 * ((*p1)-(*p2))*((*p1)-(*p2));
         if(*p3>0) weightProd *= *p3;
         p1++;p2++;p3++;
         if(++i == mNbIntegrationUsed) break;
      }
      mIntegratedChi2LikeNorm -= log(weightProd);
      weightProd=1.;
   }
   mIntegratedChi2LikeNorm/=2;
   VFN_DEBUG_MESSAGE("Chi^2="<<mIntegratedChi2<<", log(norm)="<<mIntegratedChi2LikeNorm,3)
   mClockIntegratedChi2.Click();
   VFN_DEBUG_EXIT("PowderPattern::GetChi2()="<<mIntegratedChi2,3);
   return mIntegratedChi2;
}

REAL PowderPattern::GetChi2_Option()const
{
   if(0 == mOptProfileIntegration.GetChoice()) return this->GetIntegratedChi2();
   else return this->GetChi2();
}

void PowderPattern::FitScaleFactorForR()const
{
   if(  (0==this->GetPowderPatternObs().numElements())
      ||(0==GetNbPowderPatternComponent()))
   {
      return ;
   }
   this->CalcPowderPattern();
   TAU_PROFILE("PowderPattern::FitScaleFactorForR()","void ()",TAU_DEFAULT);
   VFN_DEBUG_ENTRY("PowderPattern::FitScaleFactorForR()",3);
   // Which components are scalable ?
      mScalableComponentIndex.resize(mPowderPatternComponentRegistry.GetNb());
      int nbScale=0;
      for(int i=0;i<mPowderPatternComponentRegistry.GetNb();i++)
      {
         if(mPowderPatternComponentRegistry.GetObj(i).IsScalable())
            mScalableComponentIndex(nbScale++)=i;
      }
   VFN_DEBUG_MESSAGE("-> Number of Scale Factors:"<<nbScale<<":Index:"<<endl<<mScalableComponentIndex,3);
   if(0==nbScale)
   {
      VFN_DEBUG_EXIT("PowderPattern::FitScaleFactorForR(): No scalable component!",3);
      return;
   }
   mScalableComponentIndex.resizeAndPreserve(nbScale);
   // prepare matrices
      mFitScaleFactorM.resize(nbScale,nbScale);
      mFitScaleFactorB.resize(nbScale,1);
      mFitScaleFactorX.resize(nbScale,1);
   // Build Matrix & Vector for LSQ
   const long nbExclude=mExcludedRegionMinX.numElements();
   if(0==nbExclude)
   {
      for(int i=0;i<nbScale;i++)
      {
         for(int j=i;j<nbScale;j++)
         {
            // Here use a direct access to the powder spectrum, since
            // we know it has just been recomputed
            const REAL *p1=mPowderPatternComponentRegistry.GetObj(mScalableComponentIndex(i))
                              .mPowderPatternCalc.data();
            const REAL *p2=mPowderPatternComponentRegistry.GetObj(mScalableComponentIndex(j))
                              .mPowderPatternCalc.data();
            REAL m=0.;
            for(unsigned long k=0;k<mNbPointUsed;k++) m += *p1++ * *p2++;
            mFitScaleFactorM(i,j)=m;
            mFitScaleFactorM(j,i)=m;
         }
      }
      for(int i=0;i<nbScale;i++)
      {
         const REAL *p1=mPowderPatternObs.data();
         const REAL *p2=mPowderPatternComponentRegistry.GetObj(mScalableComponentIndex(i))
                           .mPowderPatternCalc.data();
         REAL b=0.;
         if(mPowderPatternBackgroundCalc.numElements()<=1)
            for(unsigned long k=0;k<mNbPointUsed;k++) b += *p1++ * *p2++;
         else
         {
            const REAL *p3=mPowderPatternBackgroundCalc.data();
            for(unsigned long k=0;k<mNbPointUsed;k++) b += (*p1++ - *p3++) * *p2++;
         }
         mFitScaleFactorB(i,0) =b;
      }
   }
   else
   {
      unsigned long min,max;
      for(int i=0;i<nbScale;i++)
      {
         for(int j=i;j<nbScale;j++)
         {
            // Here use a direct access to the powder spectrum, since
            // we know it has just been recomputed
            const REAL *p1=mPowderPatternComponentRegistry.GetObj(mScalableComponentIndex(i))
                              .mPowderPatternCalc.data();
            const REAL *p2=mPowderPatternComponentRegistry.GetObj(mScalableComponentIndex(j))
                              .mPowderPatternCalc.data();
            REAL m=0.;
            unsigned long l=0;
            for(int k=0;k<nbExclude;k++)
            {
               min=(unsigned long)floor(this->X2Pixel(mExcludedRegionMinX(j)));
               max=(unsigned long)ceil (this->X2Pixel(mExcludedRegionMaxX(j)));
               if(min>mNbPointUsed) break;
               if(max>mNbPointUsed)max=mNbPointUsed;
               //! min is the *beginning* of the excluded region
               for(;l<min;l++) m += *p1++ * *p2++;
               p1 += max-l;
               p2 += max-l;
               l  = max;
            }
            for(;l<mNbPointUsed;l++) m += *p1++ * *p2++;
            mFitScaleFactorM(i,j)=m;
            mFitScaleFactorM(j,i)=m;
         }
      }
      for(int i=0;i<nbScale;i++)
      {
         const REAL *p1=mPowderPatternObs.data();
         const REAL *p2=mPowderPatternComponentRegistry
                           .GetObj(mScalableComponentIndex(i))
                              .mPowderPatternCalc.data();
         REAL b=0.;
         unsigned long l=0;
         if(mPowderPatternBackgroundCalc.numElements()<=1)
         {
            for(int k=0;k<nbExclude;k++)
            {
               min=(unsigned long)floor(this->X2Pixel(mExcludedRegionMinX(k)));
               max=(unsigned long)ceil (this->X2Pixel(mExcludedRegionMaxX(k)));
               if(min>mNbPointUsed) break;
               if(max>mNbPointUsed)max=mNbPointUsed;
               //! min is the *beginning* of the excluded region
               for(;l<min;l++) b += *p1++ * *p2++;
               p1 += max-l;
               p2 += max-l;
               l  = max;
            }
            for(;l<mNbPointUsed;l++) b += *p1++ * *p2++;
         }
         else
         {
            const REAL *p3=mPowderPatternBackgroundCalc.data();
            for(int k=0;k<nbExclude;k++)
            {
               min=(unsigned long)floor(this->X2Pixel(mExcludedRegionMinX(k)));
               max=(unsigned long)ceil (this->X2Pixel(mExcludedRegionMaxX(k)));
               if(min>mNbPointUsed) break;
               if(max>mNbPointUsed)max=mNbPointUsed;
               //! min is the *beginning* of the excluded region
               for(;l<min;l++) b += (*p1++ - *p3++) * *p2++;
               p1 += max-l;
               p2 += max-l;
               l  = max;
            }
            for(;l<mNbPointUsed;l++) b += (*p1++ - *p3++) * *p2++;
         }
         mFitScaleFactorB(i,0) =b;
      }
   }
   if(1==nbScale) mFitScaleFactorX=mFitScaleFactorB(0)/mFitScaleFactorM(0);
   else
      mFitScaleFactorX=product(InvertMatrix(mFitScaleFactorM),mFitScaleFactorB);
   VFN_DEBUG_MESSAGE("B, M, X"<<endl<<mFitScaleFactorB<<endl<<mFitScaleFactorM<<endl<<mFitScaleFactorX,2)
   for(int i=0;i<nbScale;i++)
   {
      const REAL * p1=mPowderPatternComponentRegistry.GetObj(mScalableComponentIndex(i))
                        .mPowderPatternCalc.data();
      REAL * p0 = mPowderPatternCalc.data();
      const REAL s = mFitScaleFactorX(i)
                       -mScaleFactor(mScalableComponentIndex(i));
      if(ISNAN_OR_INF(s))
      {
         (*fpObjCrystInformUser)("Warning: working around NaN scale factor...");
         continue;
      }
      for(unsigned long j=0;j<mNbPointUsed;j++) *p0++ += s * *p1++;
      VFN_DEBUG_MESSAGE("-> Old:"<<mScaleFactor(mScalableComponentIndex(i)) <<" Change:"<<mFitScaleFactorX(i),2);
      mScaleFactor(mScalableComponentIndex(i)) = mFitScaleFactorX(i);
      mClockScaleFactor.Click();
      mClockPowderPatternCalc.Click();//we *did* correct the spectrum
   }
   VFN_DEBUG_EXIT("PowderPattern::FitScaleFactorForR():End",3);
}

void PowderPattern::FitScaleFactorForIntegratedR()const
{
   if(  (0==this->GetPowderPatternObs().numElements())
      ||(0==GetNbPowderPatternComponent()))
   {
      return ;
   }
   this->CalcPowderPattern();
   this->PrepareIntegratedRfactor();
   VFN_DEBUG_ENTRY("PowderPattern::FitScaleFactorForIntegratedR()",3);
   TAU_PROFILE("PowderPattern::FitScaleFactorForIntegratedR()","void ()",TAU_DEFAULT);
   // Which components are scalable ?
      mScalableComponentIndex.resize(mPowderPatternComponentRegistry.GetNb());
      int nbScale=0;
      for(int i=0;i<mPowderPatternComponentRegistry.GetNb();i++)
      {
         if(mPowderPatternComponentRegistry.GetObj(i).IsScalable())
            mScalableComponentIndex(nbScale++)=i;
      }
   VFN_DEBUG_MESSAGE("-> Number of Scale Factors:"<<nbScale<<":Index:"<<endl<<mScalableComponentIndex,2);
   if(0==nbScale)
   {
      VFN_DEBUG_EXIT("PowderPattern::FitScaleFactorForIntegratedR(): No scalable component!",3);
      return;
   }
   mScalableComponentIndex.resizeAndPreserve(nbScale);
   // prepare matrices
      mFitScaleFactorM.resize(nbScale,nbScale);
      mFitScaleFactorB.resize(nbScale,1);
      mFitScaleFactorX.resize(nbScale,1);
   // Build Matrix & Vector for LSQ
   VFN_DEBUG_MESSAGE("PowderPattern::FitScaleFactorForIntegratedR():1",2);
      const long numInterval=mIntegratedPatternMin.numElements();
      CrystVector_REAL *integratedCalc= new CrystVector_REAL[nbScale];
      for(int i=0;i<nbScale;i++)
      {
         integratedCalc[i].resize(numInterval);

         // Here use a direct access to the powder spectrum, since
         // we know it has just been recomputed
         const REAL *p1=mPowderPatternComponentRegistry.GetObj(mScalableComponentIndex(i))
                          .mPowderPatternCalc.data();

         REAL *p2=integratedCalc[i].data();
         for(int j=0;j<numInterval;j++)
         {
            const long max=mIntegratedPatternMax(j);
            p1=mPowderPatternComponentRegistry.GetObj(mScalableComponentIndex(i))
                          .mPowderPatternCalc.data()+mIntegratedPatternMin(j);
            *p2=0;
            for(int k=mIntegratedPatternMin(j);k<=max;k++) *p2 += *p1++;
            //cout <<"Calc#"<<i<<":"<< mIntegratedPatternMin(j) << " "
            //     <<mIntegratedPatternMax(j)<<" "
            //     << *p2<<endl;
            p2++;
         }
      }
   VFN_DEBUG_MESSAGE("PowderPattern::FitScaleFactorForIntegratedR():2",2);
      CrystVector_REAL backdIntegrated(numInterval);
      if(mPowderPatternBackgroundCalc.numElements()>1)
      {
         const REAL *p1;
         REAL *p2=backdIntegrated.data();
         for(int j=0;j<numInterval;j++)
         {
            const long max=mIntegratedPatternMax(j);
            p1=mPowderPatternBackgroundCalc.data()+mIntegratedPatternMin(j);
            *p2=0;
            for(int k=mIntegratedPatternMin(j);k<=max;k++) *p2 += *p1++;
            //cout <<"Backgd:"<< mIntegratedPatternMin(j) << " "
            //     <<mIntegratedPatternMax(j)<<" "
            //     << *p2<<endl;
            p2++;
         }
      }

   //if(mPowderPatternBackgroundCalc.numElements()<=1)
   //   cout<< FormatVertVector<REAL>(integratedCalc[0],mIntegratedObs,mIntegratedWeight,backdIntegrated)<<endl;
   //else
   //   cout<< FormatVertVector<REAL>(integratedCalc[0],mIntegratedObs,mIntegratedWeight)<<endl;
   VFN_DEBUG_MESSAGE("PowderPattern::FitScaleFactorForIntegratedR():3",2);
      for(int i=0;i<nbScale;i++)
      {
         for(int j=i;j<nbScale;j++)
         {
            const REAL *p1=integratedCalc[i].data();
            const REAL *p2=integratedCalc[j].data();
            REAL m=0.;
            for(long k=0;k<numInterval;k++)
            {
               m += *p1++ * *p2++;
               //cout <<"M:"<< mIntegratedPatternMin(k) << " "<<mIntegratedPatternMax(k)<<" "<<m<<endl;
            }
            mFitScaleFactorM(i,j)=m;
            mFitScaleFactorM(j,i)=m;
         }
      }
   VFN_DEBUG_MESSAGE("PowderPattern::FitScaleFactorForIntegratedR():4",2);
      for(int i=0;i<nbScale;i++)
      {
         const REAL *p1=mIntegratedObs.data();
         const REAL *p2=integratedCalc[i].data();
         REAL b=0.;
         if(mPowderPatternBackgroundCalc.numElements()<=1)
            for(long k=0;k<numInterval;k++)
            {
               b += *p1++ * *p2++;
               //cout<<"B:"<<mIntegratedPatternMin(k)<<" "<<mIntegratedPatternMax(k)<<" "<<b<<endl;
            }
         else
         {
            const REAL *p3=backdIntegrated.data();
            for(long k=0;k<numInterval;k++)
            {
               //cout<<"B(minus backgd):"<<mIntegratedPatternMin(k)<<" "
               //    <<mIntegratedPatternMax(k)<<" "
               //    <<*p1<<" "<<*p2<<" "<<*p3<<" "<<b<<endl;
               b += (*p1++ - *p3++) * *p2++;
            }
         }
         mFitScaleFactorB(i,0) =b;
      }
   VFN_DEBUG_MESSAGE("PowderPattern::FitScaleFactorForIntegratedR():5",2);

   if(1==nbScale) mFitScaleFactorX=mFitScaleFactorB(0)/mFitScaleFactorM(0);
   else
      mFitScaleFactorX=product(InvertMatrix(mFitScaleFactorM),mFitScaleFactorB);
   VFN_DEBUG_MESSAGE("B, M, X"<<endl<<mFitScaleFactorB<<endl<<mFitScaleFactorM<<endl<<mFitScaleFactorX,3)
   for(int i=0;i<nbScale;i++)
   {
      const REAL * p1=mPowderPatternComponentRegistry.GetObj(mScalableComponentIndex(i))
                        .mPowderPatternCalc.data();
      REAL * p0 = mPowderPatternCalc.data();
      const REAL s = mFitScaleFactorX(i)
                       -mScaleFactor(mScalableComponentIndex(i));
      if(ISNAN_OR_INF(s))
      {
         (*fpObjCrystInformUser)("Warning: working around NaN scale factor...");
         continue;
      }
      for(unsigned long j=0;j<mNbPointUsed;j++) *p0++ += s * *p1++;
      VFN_DEBUG_MESSAGE("-> Old:"<<mScaleFactor(mScalableComponentIndex(i)) <<" New:"<<mFitScaleFactorX(i),3);
      mScaleFactor(mScalableComponentIndex(i)) = mFitScaleFactorX(i);
      mClockScaleFactor.Click();
      mClockPowderPatternCalc.Click();//we *did* correct the spectrum
   }
   delete[] integratedCalc;
   VFN_DEBUG_EXIT("PowderPattern::FitScaleFactorForIntegratedR():End",3);
}

void PowderPattern::FitScaleFactorForRw()const
{
   if(  (0==this->GetPowderPatternObs().numElements())
      ||(0==GetNbPowderPatternComponent()))
   {
      return ;
   }
   TAU_PROFILE("PowderPattern::FitScaleFactorForRw()","void ()",TAU_DEFAULT);
   VFN_DEBUG_ENTRY("PowderPattern::FitScaleFactorForRw()",3);
   this->CalcPowderPattern();
   //:TODO: take into account excluded regions...
   // Which components are scalable ?
      mScalableComponentIndex.resize(mPowderPatternComponentRegistry.GetNb());
      int nbScale=0;
      for(int i=0;i<mPowderPatternComponentRegistry.GetNb();i++)
      {
         if(mPowderPatternComponentRegistry.GetObj(i).IsScalable())
            mScalableComponentIndex(nbScale++)=i;
      }
   VFN_DEBUG_MESSAGE("-> Number of Scale Factors:"<<nbScale<<":Index:"<<endl<<mScalableComponentIndex,2);
   if(0==nbScale)
   {
      VFN_DEBUG_EXIT("PowderPattern::FitScaleFactorForRw(): No scalable component!",3);
      return;
   }
   mScalableComponentIndex.resizeAndPreserve(nbScale);
   // prepare matrices
      mFitScaleFactorM.resize(nbScale,nbScale);
      mFitScaleFactorB.resize(nbScale,1);
      mFitScaleFactorX.resize(nbScale,1);
   // Build Matrix & Vector for LSQ
   const long nbExclude=mExcludedRegionMinX.numElements();
   if(0==nbExclude)
   {
      for(int i=0;i<nbScale;i++)
      {
         for(int j=i;j<nbScale;j++)
         {
            // Here use a direct access to the powder spectrum, since
            // we know it has just been recomputed
            const REAL *p1=mPowderPatternComponentRegistry.GetObj(mScalableComponentIndex(i))
                              .mPowderPatternCalc.data();
            const REAL *p2=mPowderPatternComponentRegistry.GetObj(mScalableComponentIndex(j))
                              .mPowderPatternCalc.data();
            const REAL *p3=mPowderPatternWeight.data();
            REAL m=0.;
            for(unsigned long k=0;k<mNbPointUsed;k++) m += *p1++ * *p2++ * *p3++;
            mFitScaleFactorM(i,j)=m;
            mFitScaleFactorM(j,i)=m;
         }
      }
      for(int i=0;i<nbScale;i++)
      {
         const REAL *p1=mPowderPatternObs.data();
         const REAL *p2=mPowderPatternComponentRegistry.GetObj(mScalableComponentIndex(i))
                           .mPowderPatternCalc.data();
         const REAL *p3=mPowderPatternWeight.data();
         REAL b=0.;
         if(mPowderPatternBackgroundCalc.numElements()<=1)
            for(unsigned long k=0;k<mNbPointUsed;k++) b += *p1++ * *p2++ * *p3++;
         else
         {
            const REAL *p4=mPowderPatternBackgroundCalc.data();
            for(unsigned long k=0;k<mNbPointUsed;k++)
               b += (*p1++ - *p4++) * *p2++ * *p3++;
         }
         mFitScaleFactorB(i,0) =b;
      }
   }
   else
   {
      unsigned long min,max;
      for(int i=0;i<nbScale;i++)
      {
         for(int j=i;j<nbScale;j++)
         {
            // Here use a direct access to the powder spectrum, since
            // we know it has just been recomputed
            const REAL *p1=mPowderPatternComponentRegistry.GetObj(mScalableComponentIndex(i))
                              .mPowderPatternCalc.data();
            const REAL *p2=mPowderPatternComponentRegistry.GetObj(mScalableComponentIndex(j))
                              .mPowderPatternCalc.data();
            const REAL *p3=mPowderPatternWeight.data();
            REAL m=0.;
            unsigned long l=0;
            for(int k=0;k<nbExclude;k++)
            {
               min=(unsigned long)floor(this->X2Pixel(mExcludedRegionMinX(j)));
               max=(unsigned long)ceil (this->X2Pixel(mExcludedRegionMaxX(j)));
               if(min>mNbPointUsed) break;
               if(max>mNbPointUsed)max=mNbPointUsed;
               //! min is the *beginning* of the excluded region
               for(;l<min;l++) m += *p1++ * *p2++ * *p3++;
               p1 += max-l;
               p2 += max-l;
               p3 += max-l;
               l  = max;
            }
            for(;l<mNbPointUsed;l++) m += *p1++ * *p2++ * *p3++;
            mFitScaleFactorM(i,j)=m;
            mFitScaleFactorM(j,i)=m;
         }
      }
      for(int i=0;i<nbScale;i++)
      {
         const REAL *p1=mPowderPatternObs.data();
         const REAL *p2=mPowderPatternComponentRegistry
                           .GetObj(mScalableComponentIndex(i))
                              .mPowderPatternCalc.data();
         const REAL *p3=mPowderPatternWeight.data();
         REAL b=0.;
         unsigned long l=0;
         if(mPowderPatternBackgroundCalc.numElements()<=1)
         {
            for(int k=0;k<nbExclude;k++)
            {
               min=(unsigned long)floor(this->X2Pixel(mExcludedRegionMinX(k)));
               max=(unsigned long)ceil (this->X2Pixel(mExcludedRegionMaxX(k)));
               if(min>mNbPointUsed) break;
               if(max>mNbPointUsed)max=mNbPointUsed;
               //! min is the *beginning* of the excluded region
               for(;l<min;l++) b += *p1++ * *p2++ * *p3++;
               p1 += max-l;
               p2 += max-l;
               p3 += max-l;
               l  = max;
            }
            for(;l<mNbPointUsed;l++) b += *p1++ * *p2++ * *p3++;
         }
         else
         {
            const REAL *p4=mPowderPatternBackgroundCalc.data();
            for(int k=0;k<nbExclude;k++)
            {
               min=(unsigned long)floor(this->X2Pixel(mExcludedRegionMinX(k)));
               max=(unsigned long)ceil (this->X2Pixel(mExcludedRegionMaxX(k)));
               if(min>mNbPointUsed) break;
               if(max>mNbPointUsed)max=mNbPointUsed;
               //! min is the *beginning* of the excluded region
               for(;l<min;l++) b += (*p1++ - *p4++) * *p2++ * *p3++;
               p1 += max-l;
               p2 += max-l;
               p3 += max-l;
               l  = max;
            }
            for(;l<mNbPointUsed;l++) b += (*p1++ - *p4++) * *p2++ * *p3++;
         }
         mFitScaleFactorB(i,0) =b;
      }
   }
   if(1==nbScale) mFitScaleFactorX=mFitScaleFactorB(0)/mFitScaleFactorM(0);
   else
      mFitScaleFactorX=product(InvertMatrix(mFitScaleFactorM),mFitScaleFactorB);
   VFN_DEBUG_MESSAGE("B, M, X"<<endl<<mFitScaleFactorB<<endl<<mFitScaleFactorM<<endl<<mFitScaleFactorX,2)
   for(int i=0;i<nbScale;i++)
   {
      const REAL * p1=mPowderPatternComponentRegistry.GetObj(mScalableComponentIndex(i))
                        .mPowderPatternCalc.data();
      REAL * p0 = mPowderPatternCalc.data();
      const REAL s = mFitScaleFactorX(i)
                       -mScaleFactor(mScalableComponentIndex(i));
      if(ISNAN_OR_INF(s))
      {
         (*fpObjCrystInformUser)("Warning: working around NaN scale factor...");
         continue;
      }
      for(unsigned long j=0;j<mNbPointUsed;j++) *p0++ += s * *p1++;
      VFN_DEBUG_MESSAGE("-> Old:"<<mScaleFactor(mScalableComponentIndex(i)) <<" Change:"<<mFitScaleFactorX(i),3);
      mScaleFactor(mScalableComponentIndex(i)) = mFitScaleFactorX(i);
      mClockScaleFactor.Click();
      mClockPowderPatternCalc.Click();//we *did* correct the spectrum
   }
   VFN_DEBUG_EXIT("PowderPattern::FitScaleFactorForRw():End",3);
}

void PowderPattern::FitScaleFactorForIntegratedRw()const
{
   if(  (0==this->GetPowderPatternObs().numElements())
      ||(0==GetNbPowderPatternComponent()))
   {
      return ;
   }
   this->CalcPowderPatternIntegrated();
   if(mClockScaleFactor>mClockPowderPatternIntegratedCalc)return;
   VFN_DEBUG_ENTRY("PowderPattern::FitScaleFactorForIntegratedRw()",3);
   TAU_PROFILE("PowderPattern::FitScaleFactorForIntegratedRw()","void ()",TAU_DEFAULT);
   // Which components are scalable ? Which scalable calculated components have a non-null variance ?
      mScalableComponentIndex.resize(mPowderPatternComponentRegistry.GetNb());
      CrystVector_REAL vScalableVarianceIndex(mPowderPatternComponentRegistry.GetNb());
      int nbScale=0;
      int nbVarCalc=0;
      for(int i=0;i<mPowderPatternComponentRegistry.GetNb();i++)
      {
         if(mPowderPatternComponentRegistry.GetObj(i).IsScalable())
            mScalableComponentIndex(nbScale++)=i;
         if(mPowderPatternComponentRegistry.GetObj(i).HasPowderPatternCalcVariance())
         {
            if(0!=mPowderPatternComponentRegistry.GetObj(i)
                    .GetPowderPatternIntegratedCalcVariance().first->numElements())
                        vScalableVarianceIndex(nbVarCalc++)=i;
         }
      }
   VFN_DEBUG_MESSAGE("-> Number of Scale Factors:"<<nbScale<<"("<<nbVarCalc
                     <<"with variance). Index:"<<endl<<mScalableComponentIndex,2);
   if(0==nbScale)
   {
      VFN_DEBUG_EXIT("PowderPattern::FitScaleFactorForIntegratedRw(): No scalable component!",3);
      return;
   }
   bool again;
   unsigned int ctagain=0;
   RECALC_SCALE_FACTOR_VARIANCE_FitScaleFactorForIntegratedRw:
   if(false)//if((nbScale==1)&&(nbVarCalc==1))
   {// Special case when using ML error, the scale appears also in the weight so we have to
    // use a 2nd-degree equation (ax^2+bx+c=0) to get the scale factor.
      double a=0.,b=0.,c=0.;// possible overflows...
      REAL newscale;
      {
         const REAL * RESTRICT p1=mIntegratedObs.data();
         const REAL * RESTRICT p2=mPowderPatternComponentRegistry.GetObj(mScalableComponentIndex(0))
                              .GetPowderPatternIntegratedCalc().first->data();
         const REAL * RESTRICT p1v=mIntegratedVarianceObs.data();
         const REAL * RESTRICT p2v=mPowderPatternComponentRegistry.GetObj(mScalableComponentIndex(0))
                              .GetPowderPatternIntegratedCalcVariance().first->data();
         if(mPowderPatternBackgroundIntegratedCalc.numElements()<=1)
         {
            for(unsigned long k=mNbIntegrationUsed;k>0;k--)
            {
               a += *p2v * *p1 * *p2;
               b += *p2 * *p2 * *p1v - *p1 * *p1 * *p2v;
               c -= *p1 * *p2 * *p1v;
               p1++;p2++;p1v++;p2v++;
            }
         }
         else
         {
            const REAL * RESTRICT p3=mPowderPatternBackgroundIntegratedCalc.data();
            for(unsigned long k=mNbIntegrationUsed;k>0;k--)
            {
               a += *p2v * (*p1 - *p3) * *p2;
               b += *p2 * *p2 * *p1v - (*p1 - *p3) * (*p1 - *p3) * *p2v;
               c -= (*p1 - *p3) * *p2 * *p1v;
               p1++;p2++;p1v++;p2v++;p3++;
            }
         }
         // Only one >0 solution to the 2nd degree equation
         newscale=(REAL)((-b+sqrt(b*b-4*a*c))/(2*a));
      }
      {// Store new scale, and correct old calculated pattern
         const REAL s = newscale-mScaleFactor(mScalableComponentIndex(0));
         const REAL s2 = newscale*newscale
                         -mScaleFactor(mScalableComponentIndex(0))
                         *mScaleFactor(mScalableComponentIndex(0));
         REAL * RESTRICT p0 = mPowderPatternIntegratedCalc.data();
         const REAL * RESTRICT p1=mPowderPatternComponentRegistry.GetObj(mScalableComponentIndex(0))
                           .GetPowderPatternIntegratedCalc().first->data();
         REAL * RESTRICT p0v = mPowderPatternVarianceIntegrated.data();
         const REAL * RESTRICT p1v=mPowderPatternComponentRegistry.GetObj(mScalableComponentIndex(0))
                           .GetPowderPatternIntegratedCalcVariance().first->data();
         REAL * RESTRICT p0w = mIntegratedWeight.data();
         for(unsigned long j=mNbIntegrationUsed;j>0;j--)
         {
            *p0++  += s * *p1++;
            *p0v += s2 * *p1v++;
            if(*p0v <=0) {*p0w++ =0;p0v++;}else *p0w++ = 1. / *p0v++;
         }
      }
      mScaleFactor(mScalableComponentIndex(0)) = newscale;
      mClockScaleFactor.Click();
      mClockPowderPatternIntegratedCalc.Click();//we *did* correct the spectrum
   }
   else
   {
      again=false;
      mScalableComponentIndex.resizeAndPreserve(nbScale);
      // prepare matrices
         mFitScaleFactorM.resize(nbScale,nbScale);
         mFitScaleFactorB.resize(nbScale,1);
         mFitScaleFactorX.resize(nbScale,1);
      // Build Matrix & Vector for LSQ
      VFN_DEBUG_MESSAGE("PowderPattern::FitScaleFactorForIntegratedRw():1",2);
         vector<const CrystVector_REAL*> integratedCalc;
         for(int i=0;i<nbScale;i++)
         {
            integratedCalc.push_back(mPowderPatternComponentRegistry.
                                        GetObj(mScalableComponentIndex(i))
                                           .GetPowderPatternIntegratedCalc().first);
         }
      VFN_DEBUG_MESSAGE("PowderPattern::FitScaleFactorForIntegratedRw():3",2);
         for(int i=0;i<nbScale;i++)
         {
            for(int j=i;j<nbScale;j++)
            {
               const REAL * RESTRICT p1=integratedCalc[i]->data();
               const REAL * RESTRICT p2=integratedCalc[j]->data();
               const REAL * RESTRICT p3;
               if(mIntegratedWeight.numElements()==0)
               {
                  p3=mIntegratedWeightObs.data();
                  if(ctagain>5) VFN_DEBUG_MESSAGE("ctagain="<<ctagain<<", using mIntegratedWeightObs",5);
               }
               else
               {
                  p3=mIntegratedWeight.data();
                  if(ctagain>5) VFN_DEBUG_MESSAGE("ctagain="<<ctagain<<", using mIntegratedWeight",5);
               }
               REAL m=0.;
               if(j==i)
                  for(unsigned long k=mNbIntegrationUsed;k>0;k--)
                     {m += *p1 * *p1 * *p3++; p1++;}
               else
                  for(unsigned long k=mNbIntegrationUsed;k>0;k--)
                     m += *p1++ * *p2++ * *p3++;
               mFitScaleFactorM(i,j)=m;
               mFitScaleFactorM(j,i)=m;
            }
         }
      VFN_DEBUG_MESSAGE("PowderPattern::FitScaleFactorForIntegratedRw():4",2);
         for(int i=0;i<nbScale;i++)
         {
            const REAL * RESTRICT p1=mIntegratedObs.data();
            const REAL * RESTRICT p2=integratedCalc[i]->data();
            const REAL * RESTRICT p4;
            if(mIntegratedWeight.numElements()==0) p4=mIntegratedWeightObs.data();
            else p4=mIntegratedWeight.data();
            REAL b=0.;
            if(mPowderPatternBackgroundIntegratedCalc.numElements()<=1)
            {
               for(unsigned long k=mNbIntegrationUsed;k>0;k--)
               {
                  b += *p1++ * *p2++ * *p4++;
                  //cout<<"B:"<<mIntegratedPatternMin(k)<<" "<<mIntegratedPatternMax(k)<<" "<<b<<endl;
               }
            }
            else
            {
               const REAL * RESTRICT p3=mPowderPatternBackgroundIntegratedCalc.data();
               for(unsigned long k=mNbIntegrationUsed;k>0;k--)
               {
                  //cout<<"B(minus backgd):"<<mIntegratedPatternMin(k)<<" "
                  //    <<mIntegratedPatternMax(k)<<" "
                  //    <<*p1<<" "<<*p2<<" "<<*p3<<" "<<b<<endl;
                  b += (*p1++ - *p3++) * *p2++ * *p4++;
               }
            }
            mFitScaleFactorB(i,0) =b;
         }
      VFN_DEBUG_MESSAGE("PowderPattern::FitScaleFactorForIntegratedRw():5",2);

      if(1==nbScale) mFitScaleFactorX=mFitScaleFactorB(0)/mFitScaleFactorM(0);
      else if(1<nbScale)
         mFitScaleFactorX=product(InvertMatrix(mFitScaleFactorM),mFitScaleFactorB);
      VFN_DEBUG_MESSAGE("B, M, X"<<endl<<mFitScaleFactorB<<endl<<mFitScaleFactorM<<endl<<mFitScaleFactorX,3)
      //Correct the variance, if necessary
      if(nbVarCalc>0)
      {
         //if(ctagain>0) cout <<"ctgain, sumvar="<<log(mPowderPatternVarianceIntegrated.sum());
         for(int i=0;i<nbScale;i++)
         {
            if(mPowderPatternComponentRegistry.GetObj(mScalableComponentIndex(i)).HasPowderPatternCalcVariance())
            {
               if(0!=mPowderPatternComponentRegistry.GetObj(mScalableComponentIndex(i))
                     .GetPowderPatternIntegratedCalcVariance().first->numElements())
               {
                  const REAL * RESTRICT p1=mPowderPatternComponentRegistry.GetObj(mScalableComponentIndex(i))
                                    .GetPowderPatternIntegratedCalcVariance().first->data();
                  //cout <<",sumvar(i)="<<log(mPowderPatternComponentRegistry.GetObj(mScalableComponentIndex(i))
                  //                    .GetPowderPatternIntegratedCalcVariance().first->sum());
                  REAL * RESTRICT p0 = mPowderPatternVarianceIntegrated.data();
                  const REAL s2 = mFitScaleFactorX(i)*mFitScaleFactorX(i)
                                   -mScaleFactor(mScalableComponentIndex(i))
                                    *mScaleFactor(mScalableComponentIndex(i));
                  for(unsigned long j=mNbIntegrationUsed;j>0;j--) *p0++ += s2 * *p1++;
               }
            }
         }
         //if(ctagain>0) cout <<" ->"<<log(mPowderPatternVarianceIntegrated.sum())
         //                   <<" , sumobsvar="<<log(mIntegratedVarianceObs.sum())<<endl;
         REAL * RESTRICT p0 = mIntegratedWeight.data();
         const REAL * RESTRICT p1=mPowderPatternVarianceIntegrated.data();
         for(unsigned long j=mNbIntegrationUsed;j>0;j--)
            if(*p1 <=0) {*p0++ =0;p1++;}
            else *p0++ = 1. / *p1++;

      }
      // Correct the calculated integrated pattern
      for(int i=0;i<nbScale;i++)
      {
         const REAL * RESTRICT p1=mPowderPatternComponentRegistry.GetObj(mScalableComponentIndex(i))
                           .GetPowderPatternIntegratedCalc().first->data();
         REAL * RESTRICT p0 = mPowderPatternIntegratedCalc.data();
         const REAL s = mFitScaleFactorX(i)
                          -mScaleFactor(mScalableComponentIndex(i));
         if(ISNAN_OR_INF(s))
         {
            (*fpObjCrystInformUser)("Warning: working around NaN scale factor...");
            continue;
         }
         if(nbVarCalc>0)
         {
            if(abs(s/mFitScaleFactorX(i))>0.001)
            {
               again=true;
                  //cout<<"log(scale) :"<<log(mScaleFactor(mScalableComponentIndex(i)))
                  //    <<"->"<<log(mFitScaleFactorX(i))<<" Again="<<ctagain++<<endl;
            }
            if((!again)&&(ctagain>0))
            {
               VFN_DEBUG_MESSAGE("log(scale) :"<<log(mScaleFactor(mScalableComponentIndex(i)))
                   <<"->"<<log(mFitScaleFactorX(i))<<" Again="<<ctagain,5);
            }
         }
         for(unsigned long j=mNbIntegrationUsed;j>0;j--) *p0++ += s * *p1++;
         if(ctagain>5)
         {
           VFN_DEBUG_MESSAGE("->log(scale) Old :"<<log(mScaleFactor(mScalableComponentIndex(i))) <<" New:"<<log(mFitScaleFactorX(i)),10);
         }
         mScaleFactor(mScalableComponentIndex(i)) = mFitScaleFactorX(i);
      }

      mClockScaleFactor.Click();
      mClockPowderPatternIntegratedCalc.Click();//we *did* correct the spectrum
   }
   if(again && (ctagain<20))
   {
     ctagain++;
     if(ctagain>5)
     {
       VFN_DEBUG_MESSAGE("PowderPattern::FitScaleFactorForIntegratedRw(), scaling again #"<<ctagain,10)
     }
     goto RECALC_SCALE_FACTOR_VARIANCE_FitScaleFactorForIntegratedRw;
   }
   VFN_DEBUG_EXIT("PowderPattern::FitScaleFactorForIntegratedRw():End",3);
}

void PowderPattern::SetSigmaToSqrtIobs()
{
   VFN_DEBUG_MESSAGE("PowderPattern::SetSigmaToSqrtIobs()",5);
   for(long i=0;i<mPowderPatternObs.numElements();i++)
      mPowderPatternObsSigma(i)=sqrt(fabs(mPowderPatternObs(i)));
}
void PowderPattern::SetWeightToInvSigmaSq(const REAL minRelatSigma)
{
   VFN_DEBUG_MESSAGE("PowderPattern::SetWeightToInvSigmaSq()",5);
   //:KLUDGE: If less than 1e-4*max, set to 0.... Do not give weight to unobserved points
   const REAL min=MaxAbs(mPowderPatternObsSigma)*minRelatSigma;
   //mPowderPatternWeight.resize(mPowderPatternObsSigma.numElements());
   REAL tmp;
   for(long i=0;i<mPowderPatternObsSigma.numElements();i++)
   {
      tmp = mPowderPatternObsSigma(i);
      if(tmp<min) mPowderPatternWeight(i)= 1./min/min;
      else  mPowderPatternWeight(i) =1./tmp/tmp;
   }
}
void PowderPattern::SetWeightToUnit()
{
   VFN_DEBUG_MESSAGE("PowderPattern::SetWeightToSinTheta()",5);
   //mPowderPatternWeight.resize(mPowderPatternObs.numElements());
   mPowderPatternWeight=1;
}
void PowderPattern::SetWeightPolynomial(const REAL a, const REAL b,
                                        const REAL c,
                                         const REAL minRelatIobs)
{
   VFN_DEBUG_MESSAGE("PowderPattern::SetWeightPolynomial()",5);
   const REAL min=MaxAbs(mPowderPatternObs)*minRelatIobs;
   REAL tmp;
   for(long i=0;i<mPowderPatternWeight.numElements();i++)
   {
      tmp=mPowderPatternObs(i);
      if(tmp<min)
         mPowderPatternWeight(i) =   1./(a+min+b*min*min+c*min*min*min);
      else mPowderPatternWeight(i) = 1./(a+tmp+b*tmp*tmp+c*tmp*tmp*tmp);
   }
}

void PowderPattern::BeginOptimization(const bool allowApproximations,
                                      const bool enableRestraints)
{
   this->Prepare();
   if(0 == mOptProfileIntegration.GetChoice()) this->FitScaleFactorForIntegratedRw();
   else this->FitScaleFactorForRw();
   this->RefinableObj::BeginOptimization(allowApproximations,enableRestraints);
}

//void PowderPattern::SetApproximationFlag(const bool allow)
//{// Do we need this ?
//   this->Prepare();
//   if(0 == mOptProfileIntegration.GetChoice()) this->FitScaleFactorForIntegratedRw();
//   else this->FitScaleFactorForRw();
//   this->RefinableObj::SetApproximationFlag(allow);
//}

void PowderPattern::GlobalOptRandomMove(const REAL mutationAmplitude,
                         const RefParType *type)
{
   if(mRandomMoveIsDone) return;
   this->RefinableObj::GlobalOptRandomMove(mutationAmplitude,type);
}

void PowderPattern::AddExcludedRegion(const REAL min,const REAL max)
{
   VFN_DEBUG_MESSAGE("PowderPattern::AddExcludedRegion()",5)
   const int num=mExcludedRegionMinX.numElements();
   if(num>0)
   {
      mExcludedRegionMinX.resizeAndPreserve(num+1);
      mExcludedRegionMaxX.resizeAndPreserve(num+1);
   }
   else
   {
      mExcludedRegionMinX.resize(1);
      mExcludedRegionMaxX.resize(1);
   }
   mExcludedRegionMinX(num)=min;
   mExcludedRegionMaxX(num)=max;

   //ensure regions are sorted by ascending 2thetamin
   CrystVector_long subs;
   subs=SortSubs(mExcludedRegionMinX);
   CrystVector_REAL tmp1,tmp2;
   tmp1=mExcludedRegionMinX;
   tmp2=mExcludedRegionMaxX;
   for(int i=0;i<mExcludedRegionMinX.numElements();i++)
   {
      mExcludedRegionMinX(i)=tmp1(subs(i));
      mExcludedRegionMaxX(i)=tmp2(subs(i));
   }
   VFN_DEBUG_MESSAGE(FormatVertVector<REAL>(mExcludedRegionMinX,mExcludedRegionMaxX),5)
   VFN_DEBUG_MESSAGE("PowderPattern::Add2ThetaExcludedRegion():End",5)
}

REAL PowderPattern::GetLogLikelihood()const
{
   REAL tmp=this->GetChi2_Option();
   if(mOptProfileIntegration.GetChoice()==0) tmp+=mIntegratedChi2LikeNorm;
   else tmp+=mChi2LikeNorm;
   return tmp;
}

unsigned int PowderPattern::GetNbLSQFunction()const{return 2;}

const CrystVector_REAL&
   PowderPattern::GetLSQCalc(const unsigned int idx) const
{
   TAU_PROFILE("PowderPattern::GetLSQCalc()","void ()",TAU_DEFAULT);
   switch(idx)
   {
      case 1:
      {
         this->CalcPowderPatternIntegrated();
         mPowderPatternUsedCalc=mPowderPatternIntegratedCalc;
         break;
      }
      default:
      {
         mPowderPatternUsedCalc=this->GetPowderPatternCalc();
         mPowderPatternUsedCalc.resizeAndPreserve(mNbPointUsed);
         break;
      }
   }
   return mPowderPatternUsedCalc;
}

const CrystVector_REAL&
   PowderPattern::GetLSQObs(const unsigned int idx) const
{
   TAU_PROFILE("PowderPattern::GetLSQObs()","void ()",TAU_DEFAULT);
   switch(idx)
   {
      case 1:
      {
         this->PrepareIntegratedRfactor();
         mPowderPatternUsedObs=mIntegratedObs;
         break;
      }
      default:
      {
         mPowderPatternUsedObs=this->GetPowderPatternObs();
         mPowderPatternUsedObs.resizeAndPreserve(mNbPointUsed);
         break;
      }
   }
   return mPowderPatternUsedObs;
}

const CrystVector_REAL&
   PowderPattern::GetLSQWeight(const unsigned int idx) const
{
   TAU_PROFILE("PowderPattern::GetLSQWeight()","void ()",TAU_DEFAULT);
   switch(idx)
   {
      case 1:
      {
         this->PrepareIntegratedRfactor();
         // :KLUDGE: When variance is used, mIntegratedWeight will change at each powder pattern calculation,
         // so this might be quite wrong...
         if(mIntegratedWeight.numElements()==0)
            mPowderPatternUsedWeight=mIntegratedWeightObs;
         else mPowderPatternUsedWeight=mIntegratedWeight;
         break;
      }
      default:
      {
         mPowderPatternUsedWeight=this->GetPowderPatternWeight();
         mPowderPatternUsedWeight.resizeAndPreserve(mNbPointUsed);
         break;
      }
   }
   return mPowderPatternUsedWeight;
}

std::map<RefinablePar*, CrystVector_REAL>& PowderPattern::GetLSQ_FullDeriv(const unsigned int idx,std::set<RefinablePar *> &vPar)
{
   TAU_PROFILE("PowderPattern::GetLSQ_FullDeriv()","void ()",TAU_DEFAULT);
   //return this->RefinableObj::GetLSQ_FullDeriv(idx,vPar);
   if(idx==1)
   {
      this->CalcPowderPatternIntegrated_FullDeriv(vPar);
      #if 0
      std::map<RefinablePar*, CrystVector_REAL> fullderiv_old;
      std::vector<const CrystVector_REAL*> v;
      int n=0;
      //cout<<"PowderPattern::GetLSQ_FullDeriv(integrated):scales:"<<mScaleFactor<<endl;
      cout<<"PowderPattern::GetLSQ_FullDeriv(integrated):parameters:"<<endl;
      for(std::set<RefinablePar*>::iterator par=vPar.begin();par!=vPar.end();++par)
      {
         v.push_back(&(mPowderPatternIntegrated_FullDeriv[*par]));
         fullderiv_old[*par]=this->GetLSQDeriv(idx,*(*par));
         v.push_back(&(fullderiv_old[*par]));
         cout<<(*par)->GetName()<<":"<<mPowderPatternIntegrated_FullDeriv[*par].size()<<","<<fullderiv_old[*par].size()<<endl;
         if(++n>8) break;
      }
      cout<<"PowderPattern::GetLSQ_FullDeriv(integrated):"<<endl<<FormatVertVector<REAL>(v,12,1,20)<<endl;
      //exit(0);
      #endif
      return mPowderPatternIntegrated_FullDeriv;
   }
   mPowderPattern_FullDeriv=this->GetPowderPattern_FullDeriv(vPar);
   return mPowderPattern_FullDeriv;
}

void PowderPattern::Prepare()
{
   VFN_DEBUG_MESSAGE("PowderPattern::Prepare()",5);
   for(int i=0;i<mPowderPatternComponentRegistry.GetNb();i++)
   {
      mPowderPatternComponentRegistry.GetObj(i).SetMaxSinThetaOvLambda(mMaxSinThetaOvLambda);
      mPowderPatternComponentRegistry.GetObj(i).Prepare();
   }
}
void PowderPattern::GetGeneGroup(const RefinableObj &obj,
                                CrystVector_uint & groupIndex,
                                unsigned int &first) const
{
   // One group for scales, one for theta error parameters
   unsigned int scaleIndex=0;
   unsigned int thetaIndex=0;
   VFN_DEBUG_MESSAGE("PowderPattern::GetGeneGroup()",4)
   for(long i=0;i<obj.GetNbPar();i++)
      for(long j=0;j<this->GetNbPar();j++)
         if(&(obj.GetPar(i)) == &(this->GetPar(j)))
         {
            if(this->GetPar(j).GetType()->IsDescendantFromOrSameAs(gpRefParTypeScattDataScale))
            {
               if(scaleIndex==0) scaleIndex=first++;
               groupIndex(i)=scaleIndex;
            }
            else //gpRefParTypeScattDataCorrPos
            {
               if(thetaIndex==0) thetaIndex=first++;
               groupIndex(i)=thetaIndex;
            }
         }
}

void PowderPattern::SetMaxSinThetaOvLambda(const REAL max)
{
   mMaxSinThetaOvLambda=max;
   for(int i=0;i<mPowderPatternComponentRegistry.GetNb();i++)
   {
      mPowderPatternComponentRegistry.GetObj(i).SetMaxSinThetaOvLambda(mMaxSinThetaOvLambda);
   }
}

REAL PowderPattern::GetMaxSinThetaOvLambda()const{return mMaxSinThetaOvLambda;}

const CrystVector_long& PowderPattern::GetIntegratedProfileMin()const
{
   this->PrepareIntegratedRfactor();
   return mIntegratedPatternMin;
}

const CrystVector_long& PowderPattern::GetIntegratedProfileMax()const
{
   this->PrepareIntegratedRfactor();
   return mIntegratedPatternMax;
}

const RefinableObjClock& PowderPattern::GetIntegratedProfileLimitsClock()const
{
   return mClockIntegratedFactorsPrep;
}

REAL PowderPattern::STOL2X(const REAL stol)const
{
   REAL x;
   if(this->GetRadiation().GetWavelengthType()==WAVELENGTH_TOF)
   {
      if(stol>0) x = 1.0/(2*stol); else return 0;
      x = mDIFC*x+mDIFA*x*x;
      VFN_DEBUG_MESSAGE("PowderPattern::STOL2X("<<stol<<","<<1.0/(2*stol)<<")="<<x,2)
   }
   else
   {
      x=stol*this->GetWavelength();
      if(abs(x)<1.0) x=2*asin(x); else x=2*M_PI;
   }
   return x;
}

REAL PowderPattern::X2STOL(const REAL x)const
{
   REAL stol;
   if(this->GetRadiation().GetWavelengthType()==WAVELENGTH_TOF)
   {
      if(abs(mDIFA)>abs(mDIFC*1e-6))
      {
         const REAL delta=mDIFC*mDIFC+4.0*mDIFA*x;
         stol = (-mDIFC+sqrt(delta))/(2.0*mDIFA);
         stol = 1/(2.0*stol);
      }
      else stol=mDIFC/(2.0*x);
   }
   else
      stol=sin(x/2.0)/this->GetWavelength();
   VFN_DEBUG_MESSAGE("PowderPattern::X2STOL("<<x<<")="<<stol,2)
   return stol;
}

REAL PowderPattern::STOL2Pixel(const REAL stol)const
{
   return this->X2Pixel(this->STOL2X(stol));
}

PeakList PowderPattern::FindPeaks(const float dmin,const float maxratio,const unsigned int maxpeak)
{
   const long nb=this->GetNbPoint() ;
   // Limit peak detection to 1.5A resolution
   long start,finish;
   if(this->GetRadiation().GetWavelengthType()!=WAVELENGTH_TOF)
   {
      start=1;// do not start at 0, if this is a simulation that really start at theta=0...
      for(finish=0;finish<nb;++finish)
      {
         const REAL d=1/(this->X2STOL(this->GetPowderPatternX()(finish))*2);
         if(d<dmin) break;
      }
   }
   else
   {
      finish=nb-1;
      for(start=nb-1;start>=0;--start)
      {
         const REAL d=1/(this->X2STOL(this->GetPowderPatternX()(start))*2);
         if(d<dmin) break;
      }
   }
   // First evaluate approximate width (in number of pixels) of reflections
   unsigned int width_golay=10;
   {
      CrystVector_REAL obs;
      const int nbwidth=9;
      CrystVector_long width(nbwidth);
      width=0;
      obs=this->GetPowderPatternObs();
      // Zero excluded regions.
      for(long i= 0;i<mExcludedRegionMinX.numElements();i++)
      {
         long min,max;
         min=(long)floor(this->X2Pixel(mExcludedRegionMinX(i)));
         max=(long)ceil (this->X2Pixel(mExcludedRegionMaxX(i)));
         if(min<0) min = 0;
         if(max>=obs.numElements()) max = obs.numElements();
         for(long j=min;j<max;j++) obs(j) = 0;
      }
      const long nb=obs.numElements();
      for(int j=0;j<nbwidth;j++)
      {
         const long imax=obs.imax(nb/10,nb-1);
         const REAL iobs_max=obs(imax);
         REAL thres=iobs_max;
         long i;
         for(i=imax-100;i<(imax+100);++i)
         {
            if(i<0){i=0;continue;}
            if(i>=nb) break;
            if(obs(i)<thres) thres=obs(i);
         }
         thres=(iobs_max+thres)/2;
         i=imax;
         while(obs(i)>=thres)
         {
            cout<<obs(i)<<"   ";
            obs(i--)=0;
            width(j)+=1;
            if(i<0) break;
         }
         i=imax+1;
         while(obs(i)>=thres)
         {
            cout<<obs(i)<<"   ";
            obs(i++)=0;
            width(j)+=1;
            if(i>=nb) break;
         }
         cout<<endl<<" => "<<width(j)<<endl;
         for(i=imax-width(j)*2;i<(imax+width(j)*2);++i)
         {
            if(i<0) continue;
            if(i>=nb) break;
            obs(i)=0;
         }
      }
      cout<<"Width of "<<nbwidth<<" strongest peaks:"<<endl<<width;
      width_golay=width(SortSubs(width)(nbwidth/2));
      cout<<"median width:"<<width_golay<<endl;
      if(width_golay<=4)width_golay=4;
      if(width_golay>=16)width_golay=16;
   }

   // get 2nd derivative
   CrystVector_REAL obsd2;
   obsd2=SavitzkyGolay(this->GetPowderPatternObs(),width_golay,2);
   // Zero excluded regions.
   for(long i= 0;i<mExcludedRegionMinX.numElements();i++)
   {
      long min,max;
      min=(long)floor(this->X2Pixel(mExcludedRegionMinX(i)));
      max=(long)ceil (this->X2Pixel(mExcludedRegionMaxX(i)));
      if(min<0) min = 0;
      if(max>=obsd2.numElements()) max = obsd2.numElements();
      for(long j=min;j<max;j++) obsd2(j) = 0;
   }
   const float norm=-obsd2.min();
   // Normalize, so that the derivative has the same extent as the observed pattern
   obsd2 *= mPowderPatternObs.max()/(-norm);

   REAL min_iobs;
   if(maxratio<0)
   {//Automatic discrimination - get an idea of noise from the distribution of the scond derivative
      CrystVector_REAL tmp;
      tmp=obsd2;
      tmp.resizeAndPreserve(tmp.numElements()/4);// First quarter, avoid too many peaks
      CrystVector<long> sub(tmp.numElements());
      QuickSortSubs(tmp,sub,tmp.numElements()-1,0);
      min_iobs=5*(tmp(tmp.numElements()/2)-tmp(tmp.numElements()/4));
      //cout<<__FILE__<<":"<<__LINE__<<" MIN_IOBS (automatic)="<<min_iobs<<endl;
   }else min_iobs=-1;// This will be set after highest peak is found
   PeakList pl;
   int nbav_min=0;//minimum numerb of points over which the peak is integrated
   while(true)
   {// Start from max
      const long imax=obsd2.imax(start,finish);
      REAL iobs=obsd2(imax);
      if(iobs<=0) break;
      REAL xmax=mX(imax)*iobs;
      long nbav=1;
      long i=imax;
      REAL lastiobs=obsd2(i);
      const REAL iobs_max=lastiobs;
      //cout<<i<<":"<<lastiobs<<":"<<mX(i)*RAD2DEG<<endl;
      while(true)
      {
         if(i<=1) break;
         if(obsd2(--i)>=lastiobs) break;
         lastiobs=obsd2(i);
         obsd2(i)=0;
         iobs+=lastiobs;
         xmax+=mX(i)*lastiobs;nbav++;
         if(lastiobs<=0) break;
         //cout<<i<<":"<<lastiobs<<":"<<mX(i)*RAD2DEG<<endl;
      }
      float dleft=mX(i+1);
      i=imax;
      lastiobs=obsd2(i);
      while(true)
      {
         if(i>=(nb-2)) break;
         if(obsd2(++i)>=lastiobs) break;
         lastiobs=obsd2(i);
         obsd2(i)=0;
         iobs+=lastiobs;
         xmax+=mX(i)*lastiobs;nbav++;
         if(lastiobs<=0) break;
         //cout<<i<<":"<<lastiobs<<":"<<mX(i)*RAD2DEG<<endl;
      }
      float dright=mX(i-1);
      xmax/=iobs;
      obsd2(imax)=0;
      REAL dmax=this->X2STOL(xmax)*2;
      dright=this->X2STOL(dright)*2;
      dleft =this->X2STOL(dleft)*2;
      //TODO : evaluate min intensity ratio from noise ?
      if(min_iobs<0)
      {
         cout<<this->GetPowderPatternObsSigma()(imax)<<","<<this->GetPowderPatternObs()(imax)<<endl;
         min_iobs=iobs/nbav*maxratio;
         //cout<<__FILE__<<":"<<__LINE__<<" MIN_IOBS="<<min_iobs<<endl;
      }
      if(nbav_min==0)
      {
         nbav_min=nbav/2;
         if(nbav_min<3)nbav_min=3;
      }
      if(pl.GetPeakList().size()>maxpeak) break;
      float sigma=this->GetPowderPatternObsSigma()(imax);
      if(this->GetRadiation().GetWavelengthType()!=WAVELENGTH_TOF)
         xmax  *= RAD2DEG;
      cout<<"Peak #"<<pl.GetPeakList().size()<<"imax="<<imax<<", x="<<xmax<<",d="<<1/dmax<<", d2iobs_max="<<iobs_max
          <<", d2Iobs="<<iobs<<", nbav="<<nbav<<", min_iobs="<<min_iobs<<",sigma="<<sigma<<endl;
      if(  ((nbav>=nbav_min)  &&(iobs_max>min_iobs)&&((iobs/nbav)>min_iobs))
         ||((nbav>=nbav_min)  &&(iobs_max>min_iobs)&&((iobs/nbav)>min_iobs*.2)&&((iobs/nbav)>3*sigma))
         ||((nbav>=nbav_min/2)&&(iobs_max>min_iobs)&&((iobs/nbav)>min_iobs*2 )&&((iobs/nbav)>6*sigma)))
      {
         pl.AddPeak(dmax,iobs,abs(dright-dleft)*.25);
         if((pl.GetPeakList().size()==1)&&(maxratio<0)&&(min_iobs<0.005*iobs/nbav)) min_iobs=0.005*iobs/nbav;
      }
   }
   pl.Print(cout);
   return pl;
}

const CrystVector_REAL& PowderPattern::GetScaleFactor() const{return mScaleFactor;}

CrystVector_REAL& PowderPattern::GetScaleFactor(){return mScaleFactor;}

//Local structures to export atoms, bond and angle restraints
struct exportAtom
{
   exportAtom(string n,REAL X, REAL Y, REAL Z,REAL b,REAL o,const ScatteringPower *pow):
   name(n),x(X),y(Y),z(Z),biso(b),occ(o),occMult(1),mpScattPow(pow){}
   string name;
   REAL x,y,z,biso,occ;
   /// Number of overlapping atoms - 1 if no overlapping, otherwise we need to multiply by occMult
   int occMult;
   const ScatteringPower *mpScattPow;
};

struct exportBond
{
   exportBond(const string &a1,const string &a2, REAL d, REAL s):
   at1(a1),at2(a2),dist(d),sigma(s){}
   string at1,at2;
   REAL dist,sigma;
};

struct exportAngle
{
   exportAngle(const string &a1,const string &a2,const string &a3, REAL a, REAL s):
   at1(a1),at2(a2),at3(a3),ang(a),sigma(s){}
   string at1,at2,at3;
   REAL ang,sigma;
};

void PowderPattern::ExportFullprof(const std::string &prefix)const
{
   // Analyze our data - background ? number of crystalline phases ?
   const PowderPatternBackground *pBackground=0;
   vector<const PowderPatternDiffraction*> vDiff;
   for(unsigned int i=0;i<this->GetNbPowderPatternComponent();i++)
   {
      if(this->GetPowderPatternComponent(i).GetClassName()=="PowderPatternBackground")
         pBackground=dynamic_cast<const PowderPatternBackground*>(&(this->GetPowderPatternComponent(i)));
      if(this->GetPowderPatternComponent(i).GetClassName()=="PowderPatternDiffraction")
         vDiff.push_back(dynamic_cast<const PowderPatternDiffraction*>(&(this->GetPowderPatternComponent(i))));
   }
   if((pBackground==0)||vDiff.size()==0) return;

   // Powder data file
   ofstream dat((prefix+".dat").c_str());
   dat<<"XYDATA"<<endl
      <<"INTER 1.0 1.0 0"<<endl<<endl<<endl<<endl<<endl;

   CrystVector_REAL ttheta;
   ttheta=mX;
   if(this->GetRadiation().GetWavelengthType()!=WAVELENGTH_TOF) ttheta *= RAD2DEG;
   dat << FormatVertVector<REAL>(ttheta,mPowderPatternObs,mPowderPatternObsSigma,12,4);
   dat.close();

   // PCR file
   ofstream pcr((prefix+".pcr").c_str());
   // if(!pcr) :TODO:

   // Title
   pcr<<"Fox/ObjCryst exported file:"<<this->GetName()<<endl;
   // Number of patterns
   pcr<<"NPATT 1"<<endl;
   // Weight of each pattern
   pcr<<"W_PAT 1.0"<<endl;
   // Multi-pattern format
   pcr<<"! Nph Dum Ias Nre Cry Opt Aut"<<endl;
   pcr<<"   "<<vDiff.size()
           <<"   0   0   0   0   1   1 "<<endl;
   // For each phase
   {
      int job=0;
      if(this->GetRadiation().GetRadiationType()==RAD_XRAY) job=0;
      if(this->GetRadiation().GetRadiationType()==RAD_NEUTRON) job=1;
      //:TODO: TOF:
      pcr<<"! Job Npr Nba Nex Nsc Nor Iwg Ilo Res Ste Uni Cor Anm"<<endl
         <<"   "<<job
              <<"   5   "<<pBackground->GetInterpPoints().first->numElements()
                       <<"  0   0   1   0   0   0   1   0   0   0"<<endl;
   }
   // Names of data files
   string shortName;
   {// Strip path for data file
      std::string::size_type idx =prefix.rfind("/");
      std::string::size_type idx2=prefix.rfind("\\");
      std::string::size_type idx3=prefix.rfind(":");
      if(((long)idx2!=(long)string::npos)&&((long)idx2>(long)idx))idx=idx2;
      if(((long)idx3!=(long)string::npos)&&((long)idx3>(long)idx))idx=idx3;
      if(idx==string::npos)
         shortName=prefix;
      else
         shortName=prefix.substr(idx+1);
   }
   pcr<<"! File names of data files"<<endl;
   pcr<<shortName<<".dat"<<endl;
   // Output options...
   pcr<<"! Mat Pcr Syo Rpa Sym Sho"<<endl
      <<"   1   1   0   -1  0   0 "<<endl;
   // Output options... For each pattern
   pcr<<"! Ipr Ppl Ioc Ls1 Ls2 Ls3 Prf Ins Hkl Fou Ana"<<endl
      <<"   0   0   0   0   0   0   1   10  0   0   1 "<<endl;
   // Fixed experimental parameters For each 2-theta pattern :TODO: Check !
   int wdt=16;
   if(this->GetRadiation().GetRadiationType()==RAD_NEUTRON) wdt=10;
   pcr<<"!lambda1 lambda2 Ratio Bkpos Wdt Cthm muR AsyLim Rpolarz -> Patt #1"<<endl
      <<this->GetRadiation().GetWavelength()(0)<<" "<<this->GetRadiation().GetWavelength()(0)
      <<                  " 0     0   "<<wdt
                                       <<"  0   0    0     0.95"<<endl;
   // Refinement parameters - changes are damped !!
   pcr<<"!NCY Eps R_at R_an R_pr R_gl"<<endl
      <<"  5  0.2  1.0  1.0  1.0  1.0"<<endl;
   // Refinement parameters & powder data range, for each 2theta pattern
   pcr<<"! Thmin Step Thmax PSD Sent0 -> Patt #1"<<endl
      <<"    0     0    0    0    0"<<endl;
   // Background points
   pcr<<"!2Theta Background for Pattern #1"<<endl;
   for(long i=0;i<pBackground->GetInterpPoints().first->numElements();i++)
      pcr<<(*(pBackground->GetInterpPoints().first))(i)*RAD2DEG<<" "
      <<(*(pBackground->GetInterpPoints().second))(i)<<" 0.0"<<endl;
   // Number of refined parameters - just use one for the scale factor !
   pcr<<"!"<<endl<<"!"<<endl<<"1 !Number of refined parameters"<<endl;
   // Powder data experimental set-up II (refinable parameters)
   pcr<<"! Zero Code Sycos Code Sysin Code Lambda Code More -> Patt #1"<<endl;
   pcr<<" "<<mXZero*RAD2DEG <<" 0.0 "
                   <<m2ThetaDisplacement*RAD2DEG <<" 0.0 "
                                <<m2ThetaTransparency*RAD2DEG <<" 0.0 "
                                         <<"0.000  0.0  0"<<endl;
   // PHASE DESCRIPTIONS
   for(unsigned int i=0;i<vDiff.size();++i)
   {
      pcr<<"!-------------------------------------------------------------------------------"<<endl
         <<"!  Data for PHASE number:   "<<i<<"  ==> Current R_Bragg for Pattern#  1:     0.00    "<<endl
         <<"!-------------------------------------------------------------------------------"<<endl;
      //Phase name
      pcr<<vDiff[i]->GetCrystal().GetName()<<endl;

      // List all atoms, remove overlapping ones
      map<int,exportAtom> vExportAtom;
      list<exportBond> vExportBond;
      list<exportAngle> vExportAngle;
      {
         CrystMatrix_REAL minDistTable;
         minDistTable=vDiff[i]->GetCrystal().GetMinDistanceTable(-1.);
         unsigned long k=0;
         // list0 is the full scattering component list with all atoms except dummies,
         // and a correct mDynPopCorr
         const ScatteringComponentList list0=vDiff[i]->GetCrystal().GetScatteringComponentList();
         for(int s=0;s<vDiff[i]->GetCrystal().GetScattererRegistry().GetNb();s++)
         {
            const ScatteringComponentList list=vDiff[i]->GetCrystal().GetScatt(s).GetScatteringComponentList();

            // If we have a Molecule, remember the names used for the atoms to describe restraints
            // We can't use the original atom names as they might not be unique in the crystal
            const Molecule *pMol=0;
            if(vDiff[i]->GetCrystal().GetScatt(s).GetClassName()=="Molecule")
               pMol=dynamic_cast<const Molecule*>(&(vDiff[i]->GetCrystal().GetScatt(s)));
            map<const MolAtom*,string> vMolAtomName;

            for(int j=0;j<list.GetNbComponent();j++)
            {
               if(0==list(j).mpScattPow) continue;//Can this happen ?
               bool redundant=false;
               for(unsigned long l=0;l<k;++l)
                  if(abs(minDistTable(l,k))<0.5)
                  {
                     map<int,exportAtom>::iterator pos=vExportAtom.find(l);
                     if(pos!=vExportAtom.end()) pos->second.occMult+=1;
                     redundant=true;//-1 means dist > 10A
                  }
               if(!redundant)
               {
                  //:TODO: avoid non-alphanumeric characters in name
                  stringstream name;
                  name<<list(j).mpScattPow->GetName()<<k+1;
                  vExportAtom.insert(make_pair(k,exportAtom(name.str(),
                                                            list(j).mX,list(j).mY,list(j).mZ,
                                                            list(j).mpScattPow->GetBiso(),
                                                            list(j).mOccupancy*list0(k).mDynPopCorr,
                                                            list(j).mpScattPow)));
                  if(pMol!=0) vMolAtomName.insert(make_pair(pMol->GetAtomList()[j],name.str()));
               }
               k++;
            }
            if(pMol!=0)
            {
               for(vector<MolBond*>::const_iterator pos=pMol->GetBondList().begin();
                     pos!=pMol->GetBondList().end();++pos)
               {
                  map<const MolAtom*,string>::const_iterator p1,p2;
                  p1=vMolAtomName.find(&((*pos)->GetAtom1()));
                  p2=vMolAtomName.find(&((*pos)->GetAtom2()));
                  if( (p1!=vMolAtomName.end()) && (p2!=vMolAtomName.end()))
                     vExportBond.push_back(exportBond(p1->second, p2->second,
                                                      (*pos)->GetLength0(),(*pos)->GetLengthSigma()));
               }

               for(vector<MolBondAngle*>::const_iterator pos=pMol->GetBondAngleList().begin();
                     pos!=pMol->GetBondAngleList().end();++pos)
               {
                  map<const MolAtom*,string>::const_iterator p1,p2,p3;
                  p1=vMolAtomName.find(&((*pos)->GetAtom1()));
                  p2=vMolAtomName.find(&((*pos)->GetAtom2()));
                  p3=vMolAtomName.find(&((*pos)->GetAtom3()));
                  if( (p1!=vMolAtomName.end()) && (p2!=vMolAtomName.end()) && (p3!=vMolAtomName.end()))
                     vExportAngle.push_back(exportAngle(p1->second, p2->second,p3->second,
                                                        (*pos)->GetAngle0(),(*pos)->GetAngleSigma()));
               }
            }
         }
         // :TODO: recognize special positions, and move the atoms on them.
         // :TODO: list atoms excluded, commented out
      }

      // Main control codes line for the phase
      //:TODO:  extract distance (Dis) and bond angle (Ang) restraints whenever possible
      //const ScatteringComponentList *pSC=&(vDiff[i]->GetCrystal().GetScatteringComponentList());
      pcr<<"!Nat Dis Ang Jbt Isy Str Furth  ATZ Nvk More"<<endl
         <<  vExportAtom.size() <<"  "<<vExportBond.size()<<"  "<<vExportAngle.size()
                    <<"   0   0   0    0    1.0  0   1"<<endl;
      pcr<<"!Jvi Jdi Hel Sol Mom Ter  N_Domains"<<endl
         <<"  0   3   0   0   0   0      0"<<endl;
      // Contribution to the patterns
      pcr<<"!Contributions (0/1) of this phase to the  patterns"<<endl
         <<" 1"<<endl;
      //
      pcr<<"!Irf Npr Jtyp  Nsp_Ref Ph_Shift for Pattern#"<<i<<endl
         <<"  0   0   0      0      0"<<endl;
      pcr<<"! Pr1    Pr2    Pr3   Brind.   Rmua   Rmub   Rmuc     for Pattern#"<<i<<endl
         <<"  1.0    1.0    1.0    1.0      0.0    0.0    0.0"<<endl;

      // Limits for distance calculations
      pcr<<"! Max_dst(dist) (angles)  Bond-Valence Calc."<<endl
         <<"    2.7000      1.5000        0"<<endl;

      // Space group symbol
      pcr<<vDiff[i]->GetCrystal().GetSpaceGroup().GetCCTbxSpg().match_tabulated_settings().hermann_mauguin()
         <<"                       <- Space Group Symbol"<<endl;
      // Atomic parameters
      pcr<<"!Atom Typ X Y Z Biso Occ In Fin N_t Spc / Codes"<<endl;
      for(map<int,exportAtom>::const_iterator pos=vExportAtom.begin();pos!=vExportAtom.end();++pos)
      {
         pcr<<pos->second.name
            <<" "<<pos->second.mpScattPow->GetSymbol()<<" "
            <<pos->second.x<<" "<<pos->second.y<<" "<<pos->second.z<<" "
            <<pos->second.biso<<" "
            <<pos->second.occ*pos->second.occMult<<" 0  0   0   0"<<endl
            <<"       0 0 0  0    0"<<endl;
      }
      // POWDER DATA-I: PROFILE PARAMETERS FOR EACH PATTERN
      REAL eta0=vDiff[0]->GetProfile().GetPar("Eta0").GetHumanValue();
      if(eta0<.01) eta0=.01;
      else if(eta0>.99) eta0=.99;
      pcr<<"!Scale Shape1 Bov Str1 Str2 Str3 Strain-Model"<<endl
         <<" 1.0 "<<eta0
                   <<" 0.0  0.0  0.0  0.0       0"<<endl
         <<" 1.0     0.0  0.0  0.0  0.0  0.0       0"<<endl;

      // :TODO: make sure the profile used corrseponds to pseudo-Voigt first !
         pcr<<"!     U     V     W     X     Y     GauSiz     LorSiz Size-Model"<<endl;
         // :NOTE: we need to separate the 3 next lines, or the same number is generated
         // three times when compiled with Visual c++ express 2008 (compiler bug ?)
         pcr<<vDiff[i]->GetProfile().GetPar("U").GetHumanValue()<<" ";
         pcr<<vDiff[i]->GetProfile().GetPar("V").GetHumanValue()<<" ";
         pcr<<vDiff[i]->GetProfile().GetPar("W").GetHumanValue()<<" ";
         pcr<<                     "  0.0   0.0      0.0        0.0 "<<endl
            << "    0.0   0.0   0.0   0.0   0.0      0.0        0.0 "<<endl;
      // Cell parameters
      pcr<<"!     a          b         c        alpha      beta       gamma      #Cell Info"<<endl
         <<vDiff[i]->GetCrystal().GetLatticePar(0)<<" "
         <<vDiff[i]->GetCrystal().GetLatticePar(1)<<" "
         <<vDiff[i]->GetCrystal().GetLatticePar(2)<<" "
         <<vDiff[i]->GetCrystal().GetLatticePar(3)*RAD2DEG<<" "
         <<vDiff[i]->GetCrystal().GetLatticePar(4)*RAD2DEG<<" "
         <<vDiff[i]->GetCrystal().GetLatticePar(5)*RAD2DEG<<endl
         <<"    0.0        0.0       0.0        0.0        0.0        0.0"<<endl;
      pcr<<"! Pref1 Pref2 alpha0 beta0 alpha1 beta1 ?"<<endl
         <<"   0.0   0.0    0.0   0.0    0.0   0.0"<<endl
         <<"   0.0   0.0    0.0   0.0    0.0   0.0"<<endl;
      // ??
      //pcr<<"! ??"<<endl
      //   <<"0.00   0.00   0.00000    0.00"<<endl;
      if(vExportBond.size()>0)
      {
         pcr<<"!Soft distance constraints"<<endl;
         for(list<exportBond>::const_iterator pos=vExportBond.begin();pos!=vExportBond.end();++pos)
         {
            pcr<<pos->at1<<" "<<pos->at2<<" 1 0 0 0 "<<pos->dist<<" "<<pos->sigma<<endl;
         }
      }
      if(vExportBond.size()>0)
      {
         pcr<<"!Soft angle constraints"<<endl;
         for(list<exportAngle>::const_iterator pos=vExportAngle.begin();pos!=vExportAngle.end();++pos)
         {
            pcr<<pos->at1<<" "<<pos->at2<<" "<<pos->at3<<" 1 1  0 0 0  0 0 0 "
               <<pos->ang*RAD2DEG<<" "<<pos->sigma*RAD2DEG<<endl;
         }
      }
   }
   pcr.close();
}

void PowderPattern::CalcPowderPattern() const
{
   this->CalcNbPointUsed();
   if(mClockPowderPatternCalc>mClockMaster) return;

   TAU_PROFILE("PowderPattern::CalcPowderPattern()","void ()",TAU_DEFAULT);
   VFN_DEBUG_ENTRY("PowderPattern::CalcPowderPattern()",3);
   if(mPowderPatternComponentRegistry.GetNb()==0)
   {
      mPowderPatternCalc.resize(mNbPoint);
      mPowderPatternCalc=0;

      mPowderPatternVariance.resize(mNbPoint);
      mPowderPatternVariance  = mPowderPatternObsSigma;
      mPowderPatternVariance *= mPowderPatternObsSigma;

      const REAL *p0 = mPowderPatternVariance.data();
      REAL *p1=mPowderPatternWeight.data();
      for(unsigned long j=0;j<mNbPoint;j++)
      {
         if(*p0 <=0) {*p1 =0.;}
         else *p1 = 1. / *p0;
         p0++;p1++;
      }
      VFN_DEBUG_EXIT("PowderPattern::CalcPowderPattern():no components!",3);
      return;
   }
   TAU_PROFILE_TIMER(timer1,"PowderPattern::CalcPowderPattern1()Calc components","", TAU_FIELD);
   TAU_PROFILE_TIMER(timer2,"PowderPattern::CalcPowderPattern2(Add spectrums-scaled)"\
                     ,"", TAU_FIELD);
   TAU_PROFILE_TIMER(timer3,"PowderPattern::CalcPowderPattern2(Add spectrums-backgd1)"\
                     ,"", TAU_FIELD);
   TAU_PROFILE_TIMER(timer4,"PowderPattern::CalcPowderPattern2(Add spectrums-backgd2)"\
                     ,"", TAU_FIELD);
   TAU_PROFILE_TIMER(timer5,"PowderPattern::CalcPowderPattern3(Variance)","", TAU_FIELD);
   TAU_PROFILE_START(timer1);
   for(int i=0;i<mPowderPatternComponentRegistry.GetNb();i++)
      mPowderPatternComponentRegistry.GetObj(i).CalcPowderPattern();
   TAU_PROFILE_STOP(timer1);
   VFN_DEBUG_MESSAGE("PowderPattern::CalcPowderPattern():Calculated components..",3);
   bool b=false;
   if(mClockPowderPatternCalc<mClockScaleFactor)
   {
      b=true;
   }
   else
      for(int i=0;i<mPowderPatternComponentRegistry.GetNb();i++)
         if(mClockPowderPatternCalc<
               mPowderPatternComponentRegistry.GetObj(i).GetClockPowderPatternCalc())
         {
            b=true;
            break;
         }

   if(false==b)
   {
      VFN_DEBUG_EXIT("PowderPattern::CalcPowderPattern():no need to recalc",3);
      return;
   }
   mPowderPatternCalc.resize(mNbPoint);
   int nbBackgd=0;//count number of background phases
   for(int i=0;i<mPowderPatternComponentRegistry.GetNb();i++)
   {//THIS SHOULD GO FASTER (PRE-FETCHING ARRAY DATA?)
      VFN_DEBUG_MESSAGE("PowderPattern::CalcPowderPattern():Adding "<< mPowderPatternComponentRegistry.GetObj(i).GetName(),3);
      if(true==mPowderPatternComponentRegistry.GetObj(i).IsScalable())
      {
         TAU_PROFILE_START(timer2);
         if(0==i)
         {
            const REAL * p1=mPowderPatternComponentRegistry.GetObj(i)
                                 .mPowderPatternCalc.data();
            REAL * p0 = mPowderPatternCalc.data();
            const REAL s = mScaleFactor(i);
            for(unsigned long j=0;j<mNbPointUsed;j++) *p0++ = s * *p1++;
            if(!(this->IsBeingRefined())) for(unsigned long j=mNbPointUsed;j<mNbPoint;j++) *p0++ = 0;
         }
         else
         {
            const REAL * p1=mPowderPatternComponentRegistry.GetObj(i)
                                 .mPowderPatternCalc.data();
            REAL * p0 = mPowderPatternCalc.data();
            const REAL s = mScaleFactor(i);
            for(unsigned long j=0;j<mNbPointUsed;j++) *p0++ += s * *p1++;
         }
          TAU_PROFILE_STOP (timer2);
      }
      else
      {// This is a background phase
         TAU_PROFILE_START(timer3);
         if(0==i)
         {
            const REAL * p1=mPowderPatternComponentRegistry.GetObj(i)
                                 .mPowderPatternCalc.data();
            REAL * p0 = mPowderPatternCalc.data();
            for(unsigned long j=0;j<mNbPointUsed;j++) *p0++ = *p1++;
            if(!(this->IsBeingRefined())) for(unsigned long j=mNbPointUsed;j<mNbPoint;j++) *p0++ = 0;
         }
         else
         {
            const REAL * p1=mPowderPatternComponentRegistry.GetObj(i)
                                 .mPowderPatternCalc.data();
            REAL * p0 = mPowderPatternCalc.data();
            for(unsigned long j=0;j<mNbPointUsed;j++) *p0++ += *p1++;

         }
         TAU_PROFILE_STOP(timer3);
         TAU_PROFILE_START(timer4);
         // The following is useless if there is only one background phase...
         if(0==nbBackgd)
         {
            mPowderPatternBackgroundCalc.resize(mNbPoint);
            REAL *p0 = mPowderPatternBackgroundCalc.data();
            const REAL *p1=mPowderPatternComponentRegistry.GetObj(i)
                              .mPowderPatternCalc.data();
            for(unsigned long j=0;j<mNbPointUsed;j++) *p0++ = *p1++;
            if(!(this->IsBeingRefined())) for(unsigned long j=mNbPointUsed;j<mNbPoint;j++) *p0++ = 0;
         }
         else
         {
            REAL *p0 = mPowderPatternBackgroundCalc.data();
            const REAL *p1=mPowderPatternComponentRegistry.GetObj(i)
                              .mPowderPatternCalc.data();
            for(unsigned long j=0;j<mNbPointUsed;j++) *p0++ += *p1++;
         }
         nbBackgd++;
         TAU_PROFILE_STOP(timer4);
      }
   }
   if(0==nbBackgd) mPowderPatternBackgroundCalc.resize(0);//:KLUDGE:

   TAU_PROFILE_START(timer5);
   // Calc variance
   {
      VFN_DEBUG_MESSAGE("PowderPattern::CalcPowderPattern():variance",2);
      mPowderPatternVariance=mPowderPatternObsSigma;
      mPowderPatternVariance *= mPowderPatternVariance;
      for(int i=0;i<mPowderPatternComponentRegistry.GetNb();i++)
      {
         if(mPowderPatternComponentRegistry.GetObj(i).HasPowderPatternCalcVariance())
         {
            if(0==mPowderPatternComponentRegistry.GetObj(i).GetPowderPatternCalcVariance()
                    .numElements()) break;

            REAL *p0 = mPowderPatternVariance.data();
            const REAL *p1=mPowderPatternComponentRegistry.GetObj(i)
                              .GetPowderPatternCalcVariance().data();

            if(true==mPowderPatternComponentRegistry.GetObj(i).IsScalable())
            {
               const REAL s2 = mScaleFactor(i) * mScaleFactor(i);
               for(unsigned long j=0;j<mNbPointUsed;j++) *p0++ += s2 * *p1++;
            }
            else for(unsigned long j=0;j<mNbPointUsed;j++) *p0++ += *p1++;
         }
      }
      REAL *p0 = mPowderPatternWeight.data();
      const REAL *p1=mPowderPatternVariance.data();
      for(unsigned long j=0;j<mNbPointUsed;j++)
         if(*p1 <=0) {*p0++ =0;p1++;}
         else *p0++ = 1. / *p1++;
   }
   mClockPowderPatternCalc.Click();
   TAU_PROFILE_STOP(timer5);
   VFN_DEBUG_EXIT("PowderPattern::CalcPowderPattern():End",3);
}

void PowderPattern::CalcPowderPattern_FullDeriv(std::set<RefinablePar*> &vPar)
{
   TAU_PROFILE("PowderPattern::CalcPowderPattern_FullDeriv()","void ()",TAU_DEFAULT);
   this->CalcPowderPattern();
   mPowderPattern_FullDeriv.clear();
   if(mPowderPatternComponentRegistry.GetNb()==0) return;
   std::vector<std::map<RefinablePar*,CrystVector_REAL>*> comps;
   for(int i=0;i<mPowderPatternComponentRegistry.GetNb();i++)
      comps.push_back(&(mPowderPatternComponentRegistry.GetObj(i).GetPowderPattern_FullDeriv(vPar)));

   for(std::set<RefinablePar *>::iterator par=vPar.begin();par!=vPar.end();++par)
   {
      if(*par==0) continue;
      for(int i=0;i<mPowderPatternComponentRegistry.GetNb();i++)
      {
         if((*par)->GetPointer()==mScaleFactor.data()+i)
         {
            mPowderPattern_FullDeriv[*par]=mPowderPatternComponentRegistry.GetObj(i).GetPowderPatternCalc();
            continue;
         }
         else
         {
            if((*(comps[i]))[*par].size()==0) continue;
         }
         if(true==mPowderPatternComponentRegistry.GetObj(i).IsScalable())
         {
            if(mPowderPattern_FullDeriv[*par].size()==0)
            {
               mPowderPattern_FullDeriv[*par].resize(mNbPoint);// :TODO: only use mNbPointUsed
               const REAL * p1=(*comps[i])[*par].data();
               REAL * p0 = mPowderPattern_FullDeriv[*par].data();
               const REAL s = mScaleFactor(i);
               for(unsigned long j=0;j<mNbPointUsed;j++) *p0++ = s * *p1++;
               for(unsigned long j=mNbPointUsed;j<mNbPoint;j++) *p0++ = 0;// :TODO: only use mNbPointUsed
            }
            else
            {
               const REAL * p1=(*comps[i])[*par].data();
               REAL * p0 = mPowderPattern_FullDeriv[*par].data();
               const REAL s = mScaleFactor(i);
               for(unsigned long j=0;j<mNbPointUsed;j++) *p0++ += s * *p1++;
            }
         }
         else
         {// This is a background phase
            if(mPowderPattern_FullDeriv[*par].size()==0)
            {
               mPowderPattern_FullDeriv[*par].resize(mNbPoint);// :TODO: only use mNbPointUsed
               const REAL * p1=(*comps[i])[*par].data();
               REAL * p0 = mPowderPattern_FullDeriv[*par].data();
               for(unsigned long j=0;j<mNbPointUsed;j++) *p0++ = *p1++;
               for(unsigned long j=mNbPointUsed;j<mNbPoint;j++) *p0++ = 0;
            }
            else
            {
               const REAL * p1=(*comps[i])[*par].data();
               REAL * p0 = mPowderPattern_FullDeriv[*par].data();
               for(unsigned long j=0;j<mNbPointUsed;j++) *p0++ += *p1++;

            }
         }
      }
   }
   #if 0
   std::map<RefinablePar*, CrystVector_REAL> oldDeriv;
   std::vector<const CrystVector_REAL*> v;
   int n=0;
   for(std::map<RefinablePar*, CrystVector_REAL>::reverse_iterator pos=mPowderPattern_FullDeriv.rbegin();pos!=mPowderPattern_FullDeriv.rend();++pos)
   {
      if(pos->first==0) continue;
      if(pos->second.size()==0) continue;

      const REAL step=pos->first->GetDerivStep();
      pos->first->Mutate(step);
      this->CalcPowderPattern();
      oldDeriv[pos->first]=mPowderPatternCalc;
      pos->first->Mutate(-2*step);
      this->CalcPowderPattern();
      oldDeriv[pos->first]-=mPowderPatternCalc;
      oldDeriv[pos->first]/=2*step;
      pos->first->Mutate(step);

      v.push_back(&(pos->second));
      v.push_back(&(oldDeriv[pos->first]));
      cout<<pos->first->GetName()<<":"<<pos->second.size()<<","<<oldDeriv[pos->first].size()<<endl;
      if(++n>8) break;
   }
   cout<<FormatVertVector<REAL>(v,16)<<endl;
   exit(0);
   #endif
}

void PowderPattern::CalcPowderPatternIntegrated() const
{
   this->CalcNbPointUsed();
   if(mClockPowderPatternIntegratedCalc>mClockMaster) return;

   this->PrepareIntegratedRfactor();
   TAU_PROFILE("PowderPattern::CalcPowderPatternIntegrated()","void ()",TAU_DEFAULT);
   VFN_DEBUG_ENTRY("PowderPattern::CalcPowderPatternIntegrated()",4);
   if(mPowderPatternComponentRegistry.GetNb()==0)
   {
      mPowderPatternIntegratedCalc.resize(0);
      mPowderPatternVarianceIntegrated.resize(0);
      VFN_DEBUG_EXIT("PowderPattern::CalcPowderPatternIntegrated():no components!",4);
      return;
   }
   TAU_PROFILE_TIMER(timer1,"PowderPattern::CalcPowderPatternIntegrated()1:Calc components",\
                     "", TAU_FIELD);
   TAU_PROFILE_TIMER(timer2,"PowderPattern::CalcPowderPatternIntegrated()2:Add comps-scaled"\
                     ,"", TAU_FIELD);
   TAU_PROFILE_TIMER(timer3,"PowderPattern::CalcPowderPatternIntegrated()2:Add backgd1"\
                     ,"", TAU_FIELD);
   TAU_PROFILE_TIMER(timer4,"PowderPattern::CalcPowderPatternIntegrated()2:Add backgd2"\
                     ,"", TAU_FIELD);
   TAU_PROFILE_TIMER(timer5,"PowderPattern::CalcPowderPatternIntegrated()3:Variance"
                     ,"", TAU_FIELD);
   TAU_PROFILE_START(timer1);
   vector< pair<const CrystVector_REAL*,const RefinableObjClock*> > comps;
   for(int i=0;i<mPowderPatternComponentRegistry.GetNb();i++)
   {
      comps.push_back(mPowderPatternComponentRegistry.GetObj(i).
                      GetPowderPatternIntegratedCalc());
   }
   TAU_PROFILE_STOP(timer1);
   bool b=false;
   if(mClockPowderPatternCalc<mClockScaleFactor)
   {
      b=true;
   }
   else
      for(vector< pair<const CrystVector_REAL*,const RefinableObjClock*> >::iterator
          pos=comps.begin();pos!=comps.end();++pos)
         if(mClockPowderPatternCalc < *(pos->second) )
         {
            b=true;
            break;
         }

   if(false==b)
   {
      VFN_DEBUG_EXIT("PowderPattern::CalcPowderPatternIntegrated():no need to recalc",4);
      return;
   }
   VFN_DEBUG_MESSAGE("PowderPattern::CalcPowderPatternIntegrated():Recalc",3);
   mPowderPatternIntegratedCalc.resize(mNbIntegrationUsed);
   int nbBackgd=0;//count number of background phases
   for(int i=0;i<mPowderPatternComponentRegistry.GetNb();i++)
   {
      VFN_DEBUG_MESSAGE("PowderPattern::CalcPowderPatternIntegrated():Adding "\
                        << mPowderPatternComponentRegistry.GetObj(i).GetName(),3);
      if(true==mPowderPatternComponentRegistry.GetObj(i).IsScalable())
      {
         TAU_PROFILE_START(timer2);
         if(0==i)
         {
            const REAL * RESTRICT p1= comps[i].first->data();
            REAL * RESTRICT p0 = mPowderPatternIntegratedCalc.data();
            const REAL s = mScaleFactor(i);
            for(unsigned long j=mNbIntegrationUsed;j>0;j--) *p0++ = s * *p1++;
         }
         else
         {
            const REAL * RESTRICT p1= comps[i].first->data();
            REAL * RESTRICT p0 = mPowderPatternIntegratedCalc.data();
            const REAL s = mScaleFactor(i);
            for(unsigned long j=mNbIntegrationUsed;j>0;j--) *p0++ += s * *p1++;
         }
          TAU_PROFILE_STOP (timer2);
      }
      else
      {// This is a background phase
         TAU_PROFILE_START(timer3);
         if(0==i)
         {
            const REAL * RESTRICT p1= comps[i].first->data();
            REAL * RESTRICT p0 = mPowderPatternIntegratedCalc.data();
            for(unsigned long j=mNbIntegrationUsed;j>0;j--) *p0++ = *p1++;
         }
         else
         {
            const REAL * RESTRICT p1= comps[i].first->data();
            REAL * RESTRICT p0 = mPowderPatternIntegratedCalc.data();
            for(unsigned long j=mNbIntegrationUsed;j>0;j--) *p0++ += *p1++;

         }
         TAU_PROFILE_STOP(timer3);
         TAU_PROFILE_START(timer4);
         // The following is useless if there is only one background phase...
         if(0==nbBackgd)
         {
            mPowderPatternBackgroundIntegratedCalc.resize(mNbIntegrationUsed);
            const REAL * RESTRICT p1= comps[i].first->data();
            REAL * RESTRICT p0 = mPowderPatternBackgroundIntegratedCalc.data();
            for(unsigned long j=mNbIntegrationUsed;j>0;j--) *p0++ = *p1++;
         }
         else
         {
            const REAL * RESTRICT p1= comps[i].first->data();
            REAL * RESTRICT p0 = mPowderPatternBackgroundIntegratedCalc.data();
            for(unsigned long j=mNbIntegrationUsed;j>0;j--) *p0++ += *p1++;
         }
         nbBackgd++;
         TAU_PROFILE_STOP(timer4);
      }
   }
   TAU_PROFILE_START(timer5);
   if(0==nbBackgd) mPowderPatternBackgroundIntegratedCalc.resize(0);
   // Calc variance
   bool useCalcVariance=false;
   for(int i=0;i<mPowderPatternComponentRegistry.GetNb();i++)
      if(mPowderPatternComponentRegistry.GetObj(i).HasPowderPatternCalcVariance())
         if(mPowderPatternComponentRegistry.GetObj(i)
               .GetPowderPatternIntegratedCalcVariance().first->numElements() !=0)
                  useCalcVariance=true;
   if(useCalcVariance)
   {
      VFN_DEBUG_MESSAGE("PowderPattern::CalcPowderPatternIntegrated():variance",3);
      {
         mPowderPatternVarianceIntegrated.resize(mNbIntegrationUsed);
         mIntegratedWeight.resize(mNbIntegrationUsed);
         const REAL * RESTRICT p1= mIntegratedVarianceObs.data();
         REAL * RESTRICT p0 = mPowderPatternVarianceIntegrated.data();
         for(unsigned long j=mNbIntegrationUsed;j>0;j--) *p0++ = *p1++;
      }
      //cout <<"PowderPattern::CalcPowderPatternIntegrated():variance"
      //     <<"obsvarsum="<<log(mIntegratedVarianceObs.sum());
      for(int i=0;i<mPowderPatternComponentRegistry.GetNb();i++)
      {
         if(mPowderPatternComponentRegistry.GetObj(i).HasPowderPatternCalcVariance())
         {
            if(0==mPowderPatternComponentRegistry.GetObj(i)
                    .GetPowderPatternIntegratedCalcVariance().first->numElements()) break;

            const REAL * RESTRICT p1= mPowderPatternComponentRegistry.GetObj(i)
                                .GetPowderPatternIntegratedCalcVariance().first->data();
            //cout <<",sumvar(i)="<<log(mPowderPatternComponentRegistry.GetObj(i)
            //                    .GetPowderPatternIntegratedCalcVariance().first->sum());
            REAL * RESTRICT p0 = mPowderPatternVarianceIntegrated.data();

            if(true==mPowderPatternComponentRegistry.GetObj(i).IsScalable())
            {
               const REAL s2 = mScaleFactor(i) * mScaleFactor(i);
               for(unsigned long j=mNbIntegrationUsed;j>0;j--) *p0++ += s2 * *p1++;
            }
            else for(unsigned long j=mNbIntegrationUsed;j>0;j--) *p0++ += *p1++;

         }
      }
      //cout <<endl;
      REAL *p0 = mIntegratedWeight.data();
      const REAL *p1=mPowderPatternVarianceIntegrated.data();
      for(unsigned long j=mNbIntegrationUsed;j>0;j--)
         if(*p1 <=0) {*p0++ =0;p1++;}
         else *p0++ = 1. / *p1++;
   }
   else mIntegratedWeight.resize(0);
   mClockPowderPatternIntegratedCalc.Click();
   TAU_PROFILE_STOP(timer5);
   /*
   // Compare-DEBUG ONLY
   {
      VFN_DEBUG_MESSAGE("PowderPattern::CalcPowderPatternIntegrated():Check",10);
      this->CalcPowderPattern();
      CrystVector_REAL integr(mNbIntegrationUsed);
      for(int k=0;k<mPowderPatternComponentRegistry.GetNb();k++)
      {
         VFN_DEBUG_MESSAGE("PowderPattern::CalcPowderPatternIntegrated():Check #"<<k,10);
         integr=0;
         const CrystVector_REAL *v
            =&(mPowderPatternComponentRegistry.GetObj(k).GetPowderPatternCalc());
         for(unsigned int i=0;i<mNbIntegrationUsed;i++)
         {
            integr(i)=0;
            for(int j=mIntegratedPatternMin(i);j<=mIntegratedPatternMax(i);j++)
            {
               integr(i) += (*v)(j);
            }
         }
         cout << "Integrated intensities, Component #"<<k<<endl
              << FormatVertVector<REAL> (integr,*(comps[k].first))<<endl;
      }
   }
   */
   VFN_DEBUG_EXIT("PowderPattern::CalcPowderPatternIntegrated():End",4);
}

void PowderPattern::CalcPowderPatternIntegrated_FullDeriv(std::set<RefinablePar *> &vPar)
{
   TAU_PROFILE("PowderPattern::CalcPowderPatternIntegrated_FullDeriv()","void ()",TAU_DEFAULT);
   this->CalcPowderPatternIntegrated();
   this->CalcNbPointUsed();

   this->PrepareIntegratedRfactor();
   mPowderPatternUsed_FullDeriv.clear();
   if(mPowderPatternComponentRegistry.GetNb()==0) return;
   std::vector<map<RefinablePar*,CrystVector_REAL>*> comps;
   for(int i=0;i<mPowderPatternComponentRegistry.GetNb();i++)
   {
      comps.push_back(&(mPowderPatternComponentRegistry.GetObj(i).
                        GetPowderPatternIntegrated_FullDeriv(vPar)));
   }
   //RefinablePar *scalePar;
   //cout<<"PowderPattern::CalcPowderPatternIntegrated_FullDeriv():"<<endl;
   mPowderPatternIntegrated_FullDeriv.clear();
   for(std::set<RefinablePar *>::iterator par=vPar.begin();par!=vPar.end();++par)
   {
      if(*par==0) continue; //:TODO: store the calculated (non-derived) pattern here ?
      for(int i=0;i<mPowderPatternComponentRegistry.GetNb();i++)
      {
         if((*par)->GetPointer()==mScaleFactor.data()+i)
         {
            //scalePar=*par;
            mPowderPatternIntegrated_FullDeriv[*par]=*(mPowderPatternComponentRegistry.GetObj(i).GetPowderPatternIntegratedCalc().first);
            //cout<<"PowderPattern::CalcPowderPatternIntegrated_FullDeriv():scale #"<<i<<":"<<(*par)->GetName()<<":"<<(*par)->GetPointer()<<":"<<mPowderPatternUsed_FullDeriv[*par]<<endl;
            continue;
         }
         else
         {
            if((*(comps[i]))[*par].size()==0) continue;
         }
         if(true==mPowderPatternComponentRegistry.GetObj(i).IsScalable())
         {
            if(mPowderPatternIntegrated_FullDeriv[*par].size()==0)
            {
               mPowderPatternIntegrated_FullDeriv[*par].resize(mNbIntegrationUsed);
               const REAL * RESTRICT p1= (*(comps[i]))[*par].data();
               REAL * RESTRICT p0 = mPowderPatternIntegrated_FullDeriv[*par].data();
               const REAL s = mScaleFactor(i);
               for(unsigned long j=mNbIntegrationUsed;j>0;j--)
               {
                  #if 1
                  *p0++ = s * *p1++;
                  #else
                  *p0 = s * *p1;
                  if((j==mNbIntegrationUsed)&&((*par)->GetName()=="Cimetidine_C11_x")) cout<<__FILE__<<":"<<__LINE__<<":"<<*p0<<","<<s<<"*"<<*p1<<endl;
                  p0++;p1++;
                  #endif
               }
            }
            else
            {
               const REAL * RESTRICT p1= (*(comps[i]))[*par].data();
               REAL * RESTRICT p0 = mPowderPatternIntegrated_FullDeriv[*par].data();
               const REAL s = mScaleFactor(i);
               for(unsigned long j=mNbIntegrationUsed;j>0;j--)
               {
                  #if 1
                  *p0++ += s * *p1++;
                  #else
                  *p0 += s * *p1;
                  if((j==mNbIntegrationUsed)&&((*par)->GetName()=="Cimetidine_C11_x")) cout<<__FILE__<<":"<<__LINE__<<":"<<*p0<<","<<s<<"*"<<*p1<<endl;
                  p0++;p1++;
                  #endif
               }
            }
         }
         else
         {// This is a background phase
            if(mPowderPatternIntegrated_FullDeriv[*par].size()==0)
            {
               mPowderPatternIntegrated_FullDeriv[*par].resize(mNbIntegrationUsed);
               const REAL * RESTRICT p1= (*(comps[i]))[*par].data();
               REAL * RESTRICT p0 = mPowderPatternIntegrated_FullDeriv[*par].data();
               for(unsigned long j=mNbIntegrationUsed;j>0;j--)
               {
                  #if 1
                  *p0++ = *p1++;
                  #else
                  *p0 = *p1;
                  if((j==mNbIntegrationUsed)&&((*par)->GetName()=="Cimetidine_C11_x")) cout<<__FILE__<<":"<<__LINE__<<":"<<*p0<<","<<*p1<<endl;
                  p0++;p1++;
                  #endif
               }
            }
            else
            {
               const REAL * RESTRICT p1= (*(comps[i]))[*par].data();
               REAL * RESTRICT p0 = mPowderPatternIntegrated_FullDeriv[*par].data();
               for(unsigned long j=mNbIntegrationUsed;j>0;j--)
               {
                  #if 1
                  *p0++ += *p1++;
                  #else
                  *p0 += *p1;
                  if((j==mNbIntegrationUsed)&&((*par)->GetName()=="Cimetidine_C11_x")) cout<<__FILE__<<":"<<__LINE__<<":"<<*p0<<","<<*p1<<endl;
                  p0++;p1++;
                  #endif
               }
            }
         }
         #if 0
         if((*par)->GetName()=="Cimetidine_C11_x")
         cout<<"PowderPattern::CalcPowderPatternIntegrated_FullDeriv():"
             <<(*par)->GetName()<<":s="<<mScaleFactor(i)<<", d[0]="<<(*(comps[i]))[*par](0)
             <<", integ="<<mPowderPatternIntegrated_FullDeriv[*par](0)<<endl;
         #endif
      }
   }


   //cout<<"PowderPattern::CalcPowderPatternIntegrated_FullDeriv():scale #"<<1<<":"<<scalePar->GetName()<<":"<<scalePar->GetPointer()<<":"<<mPowderPatternUsed_FullDeriv[scalePar]<<endl;
}

void PowderPattern::Init()
{
   this->InitOptions();
   {
      RefinablePar tmp("Zero",&mXZero,-.05,.05,gpRefParTypeScattDataCorrPos,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,RAD2DEG);
      tmp.AssignClock(mClockPowderPatternXCorr);
      tmp.SetDerivStep(1e-6);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("2ThetaDispl",&m2ThetaDisplacement,-.05,.05,gpRefParTypeScattDataCorrPos,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,RAD2DEG);
      tmp.AssignClock(mClockPowderPatternXCorr);
      tmp.SetDerivStep(1e-6);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("2ThetaTransp",&m2ThetaTransparency,-.05,.05,gpRefParTypeScattDataCorrPos,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,RAD2DEG);
      tmp.AssignClock(mClockPowderPatternXCorr);
      tmp.SetDerivStep(1e-6);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("DIFC",&mDIFC,0,1e6,gpRefParTypeScattDataCorrPos,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,1.0);
      tmp.AssignClock(mClockPowderPatternXCorr);
      tmp.SetDerivStep(1e-2);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("DIFA",&mDIFA,-1e4,1e4,gpRefParTypeScattDataCorrPos,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,1.0);
      tmp.AssignClock(mClockPowderPatternXCorr);
      tmp.SetDerivStep(1e-4);
      this->AddPar(tmp);
   }
}
void PowderPattern::PrepareIntegratedRfactor()const
{
   bool needPrep=false;
   for(int i=0;i<mPowderPatternComponentRegistry.GetNb();i++)
   {
      mPowderPatternComponentRegistry.GetObj(i).GetBraggLimits();
      if(mPowderPatternComponentRegistry.GetObj(i).GetClockBraggLimits()
            >mClockIntegratedFactorsPrep)
      {
         needPrep=true;
         break;
      }
   }

   // If using max sin(theta)/lambda
      this->CalcNbPointUsed();
      if(mClockIntegratedFactorsPrep<mClockNbPointUsed) needPrep=true;

   if(false==needPrep)return;
   VFN_DEBUG_ENTRY("PowderPattern::PrepareIntegratedRfactor()",3);
   TAU_PROFILE("PowderPattern::PrepareIntegratedRfactor()","void ()",TAU_DEFAULT);

   // Aggregate all limiting pixels in a single list
   {
      list<long> vLimits;

      for(int i=0;i<mPowderPatternComponentRegistry.GetNb();i++)
      {
         const CrystVector_long vLim=mPowderPatternComponentRegistry.GetObj(i).GetBraggLimits();
         for(i=0;i<vLim.numElements();i++) vLimits.push_back(vLim(i));
      }
      if(vLimits.size()<2)
      {
         mIntegratedPatternMin.resize(0);
         mIntegratedPatternMax.resize(0);
         mNbIntegrationUsed=0;
         mClockIntegratedFactorsPrep.Click();
         return;
         VFN_DEBUG_EXIT("PowderPattern::PrepareIntegratedRfactor(): no intervals !",3);
      }
      vLimits.sort();
      if(*(vLimits.begin())<0)
      {
         vLimits.push_back(0);
         vLimits.sort();
      }
      for(list<long>::iterator pos=vLimits.begin();pos!=vLimits.end();)
      {
         if( (*pos<0) || (*pos>=long(mNbPointUsed)) ) pos=vLimits.erase(pos);
         else ++pos;
      }

      // Try to avoid too small intervals
      list<long> vLimits2;
      list<long>::iterator pos1=vLimits.begin();
      list<long>::iterator pos2=pos1;pos2++;
      for(;pos2!=vLimits.end();)
      {
         const long pix1=*pos1;
         //cout<<__FILE__<<":"<<__LINE__<<":"<<pix1<<endl;
         vLimits2.push_back(pix1);
         for(;;)
         {
            pos1=pos2++;
            if(pos2==vLimits.end()) break;
            if(*pos2>(pix1+8)) break;
         }
      }
      vLimits2.push_back(*pos1);

      // Try to avoid too small intervals (2nd pass)
      pos1=vLimits2.begin();
      pos2=pos1;pos2++;
      for(;pos2!=vLimits2.end();)
      {
         //cout<<__FILE__<<":"<<__LINE__<<":"<<*pos1<<" -> "<<*pos2<<endl;
         if( *pos2<=((*pos1)+2))
         {
            //cout<<__FILE__<<":"<<__LINE__<<":"<<*pos1<<" -> "<<*pos2<<"...PLONK";
            pos2=vLimits2.erase(pos2);
            //cout<<"->"<<*pos2<<endl;
         }
         else {pos1++;pos2++;}
      }

      // Create min/max pairs
      list<pair<long,long> > vLimits3;
      pos1=vLimits2.begin();
      pos2=pos1;pos2++;
      for(;pos2!=vLimits2.end();)
      {
         if(*pos2!=(long(mNbPointUsed)-1)) vLimits3.push_back(make_pair(*pos1++,*pos2++-1));
         else vLimits3.push_back(make_pair(*pos1++,*pos2++));
         //cout<<__FILE__<<":"<<__LINE__<<":"<<vLimits3.back().first<<" -> "<<vLimits3.back().second<<endl;
      }

      mIntegratedPatternMin.resize(vLimits3.size());
      mIntegratedPatternMax.resize(vLimits3.size());
      unsigned long i=0;
      for(list<pair<long,long> >::iterator pos=vLimits3.begin();pos!=vLimits3.end();++pos)
      {
         mIntegratedPatternMin(i)=pos->first;
         mIntegratedPatternMax(i++)=pos->second;
      }
   }
   long numInterval=mIntegratedPatternMin.numElements();
   CrystVector_bool keep(numInterval);
   keep=true;
   // Take care of excluded regions (change integration areas accordingly)
   // regions are sorted by ascending theta
      const long nbExclude=mExcludedRegionMinX.numElements();
      if(nbExclude>0)
      {
         VFN_DEBUG_MESSAGE("PowderPattern::PrepareIntegratedRfactor():5:Excluded regions("<<nbExclude<<")",3);
         long j=0;
         long minExcl,maxExcl;
         minExcl=(long)floor(this->X2Pixel(mExcludedRegionMinX(j)));
         maxExcl=(long)ceil (this->X2Pixel(mExcludedRegionMaxX(j)));
         for(int i=0;i<nbExclude;i++)
         {
            while(mIntegratedPatternMax(j)<minExcl)
            {
               j++;
               if(j>=numInterval) break;
            }
            if(j>=numInterval) break;
            while(mIntegratedPatternMin(j)<maxExcl)
            {
               if( (mIntegratedPatternMin(j)>minExcl) &&(mIntegratedPatternMax(j)<maxExcl))
                  keep(j)=false;
               if( (mIntegratedPatternMin(j)<minExcl) &&(mIntegratedPatternMax(j)<maxExcl))
                  mIntegratedPatternMax(j)=minExcl;
               if( (mIntegratedPatternMin(j)>minExcl) &&(mIntegratedPatternMax(j)>maxExcl))
                  mIntegratedPatternMin(j)=maxExcl;
               if(j==(numInterval-1)) break;
               j++;
            }
            minExcl=(long)(this->X2Pixel(mExcludedRegionMinX(i)));
            maxExcl=(long)(this->X2Pixel(mExcludedRegionMaxX(i)));
            //go back if one integration segment is concerned by several exclusion zones...
            if(j!=0)
               while(mIntegratedPatternMax(j)>=minExcl)
               {
                  j--;
                  if(j==0) break;
               }
         }
      }
   // Keep only the selected intervals
   VFN_DEBUG_MESSAGE("PowderPattern::PrepareIntegratedRfactor():6",3);
   long j=0;
   for(int i=0;i<numInterval;i++)
   {
      if(keep(i))
      {
         mIntegratedPatternMin(j  )=mIntegratedPatternMin(i);
         mIntegratedPatternMax(j++)=mIntegratedPatternMax(i);
      }
   }
   numInterval=j;
   mIntegratedPatternMax.resizeAndPreserve(numInterval);
   mIntegratedPatternMin.resizeAndPreserve(numInterval);

   VFN_DEBUG_MESSAGE("PowderPattern::PrepareIntegratedRfactor():intervals"<<endl\
                  <<FormatVertVector<long>(mIntegratedPatternMin,mIntegratedPatternMax),2);
   // Integrate Obs and weight arrays
   mIntegratedObs.resize(numInterval);
   mIntegratedVarianceObs.resize(numInterval);
   mIntegratedVarianceObs=0;
   mIntegratedObs=0;
   mIntegratedWeightObs.resize(numInterval);
   for(int i=0;i<numInterval;i++)
   {
      for(int j=mIntegratedPatternMin(i);j<=mIntegratedPatternMax(i);j++)
      {
         mIntegratedObs   (i)+=mPowderPatternObs(j);
         mIntegratedVarianceObs(i)+=mPowderPatternObsSigma(j)*mPowderPatternObsSigma(j);
      }
      if(mIntegratedVarianceObs(i) <= 0) mIntegratedWeightObs(i)=0;
      else mIntegratedWeightObs(i)=1./mIntegratedVarianceObs(i);
   }


   //cout<<FormatVertVector<REAL>(mIntegratedPatternMin,
   //                               mIntegratedPatternMax,
   //                               mIntegratedObs,mIntegratedWeight,12,6)<<endl;
   mNbIntegrationUsed=mIntegratedPatternMin.numElements();
   mClockIntegratedFactorsPrep.Click();
   VFN_DEBUG_EXIT("PowderPattern::PrepareIntegratedRfactor()",3);
}
void PowderPattern::CalcNbPointUsed()const
{
   if(this->IsBeingRefined())return;
   unsigned long tmp;
   // Use the first point of the profile of the first reflection not calculated
   if(this->GetRadiation().GetWavelengthType()==WAVELENGTH_TOF)
   {
      tmp=(unsigned long)(this->X2PixelCorr(this->STOL2X(mMaxSinThetaOvLambda)));
   }
   else
   {
      REAL sinth=mMaxSinThetaOvLambda*this->GetWavelength();
      if(1>fabs(sinth)) tmp=(unsigned long)(this->X2PixelCorr(2*asin(sinth))); else tmp=mNbPoint;
   }
   if(tmp>mNbPoint) tmp= mNbPoint;
   if(tmp !=mNbPointUsed)
   {
      mNbPointUsed=tmp;
      mClockNbPointUsed.Click();
      VFN_DEBUG_MESSAGE("PowderPattern::CalcNbPointUsed():"<<mNbPointUsed<<" max(sin(theta)/lambda)="<<mMaxSinThetaOvLambda, 3)
   }

}

void PowderPattern::InitOptions()
{
   VFN_DEBUG_MESSAGE("PowderPattern::InitOptions()",5)
   static string OptProfileIntegrationName;
   static string OptProfileIntegrationChoices[2];

   static bool needInitNames=true;
   if(true==needInitNames)
   {
      OptProfileIntegrationName="Use Integrated Profiles";
      OptProfileIntegrationChoices[0]="Yes (recommended)";
      OptProfileIntegrationChoices[1]="No";

      needInitNames=false;//Only once for the class
   }
   mOptProfileIntegration.Init(2,&OptProfileIntegrationName,OptProfileIntegrationChoices);
   this->AddOption(&mOptProfileIntegration);
}

#ifdef __WX__CRYST__
WXCrystObjBasic* PowderPattern::WXCreate(wxWindow* parent)
{
   //:TODO: Check mpWXCrystObj==0
   mpWXCrystObj=new WXPowderPattern(parent,this);
   return mpWXCrystObj;
}
#endif

}//namespace ObjCryst
