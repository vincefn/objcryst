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
*  source file for LibCryst++ DiffractionData class
*
*/

#include <cmath>

#include <typeinfo>

#include "ObjCryst/DiffractionDataSingleCrystal.h"
#include "Quirks/VFNDebug.h"
#include "Quirks/VFNStreamFormat.h"

#ifdef __WX__CRYST__
#include "wxCryst/wxDiffractionSingleCrystal.h"
#endif

#include <fstream>
#include <iomanip>
#include <stdio.h> //for sprintf()

//#include <xmmintrin.h>

#ifdef _MSC_VER // MS VC++ predefined macros....
#undef min
#undef max
#endif

namespace ObjCryst
{
//######################################################################
//    DiffractionDataSingleCrystal
//######################################################################
ObjRegistry<DiffractionDataSingleCrystal> 
   gDiffractionDataSingleCrystalRegistry("Global DiffractionDataSingleCrystal Registry");

DiffractionDataSingleCrystal::DiffractionDataSingleCrystal():
mHasObservedData(false),mScaleFactor(1.)
{
   VFN_DEBUG_MESSAGE("DiffractionDataSingleCrystal::DiffractionDataSingleCrystal()",5)
   this->InitRefParList();
   this->InitOptions();
   gDiffractionDataSingleCrystalRegistry.Register(*this);
   gTopRefinableObjRegistry.Register(*this);
   mClockMaster.AddChild(mClockScaleFactor);
   this->AddSubRefObj(mRadiation);
}
DiffractionDataSingleCrystal::DiffractionDataSingleCrystal(Crystal &cryst):
mHasObservedData(false),mScaleFactor(1.)
{
   VFN_DEBUG_MESSAGE("DiffractionDataSingleCrystal::DiffractionDataSingleCrystal()",5)
   this->InitRefParList();
   this->SetCrystal(cryst);
   this->InitOptions();
   gDiffractionDataSingleCrystalRegistry.Register(*this);
   gTopRefinableObjRegistry.Register(*this);
   mClockMaster.AddChild(mClockScaleFactor);
   this->AddSubRefObj(mRadiation);
}

DiffractionDataSingleCrystal::DiffractionDataSingleCrystal(const DiffractionDataSingleCrystal &old):
ScatteringData(old),
mHasObservedData(old.mHasObservedData),mRadiation(old.mRadiation)
{
   mObsIntensity=old.mObsIntensity;
   // Keep a copy as squared F(hkl), to enable fourier maps
   // :TODO: stop using mObsIntensity and just keep mFhklObsSq ?
   mFhklObsSq=mObsIntensity;
   mClockGetFhklObsSq.Click();
   
   mObsSigma=old.mObsSigma;
   mWeight=old.mWeight;
   mCalcIntensity=old.mCalcIntensity;
   mScaleFactor=old.mScaleFactor;
   this->InitOptions();
   mGroupOption.SetChoice(old.mGroupOption.GetChoice());
   gDiffractionDataSingleCrystalRegistry.Register(*this);
   gTopRefinableObjRegistry.Register(*this);
   mClockMaster.AddChild(mClockScaleFactor);
   this->AddSubRefObj(mRadiation);
}

DiffractionDataSingleCrystal::~DiffractionDataSingleCrystal()
{
   VFN_DEBUG_MESSAGE("DiffractionDataSingleCrystal::~DiffractionDataSingleCrystal()",5)
   gDiffractionDataSingleCrystalRegistry.DeRegister(*this);
   gTopRefinableObjRegistry.DeRegister(*this);
}

DiffractionDataSingleCrystal* DiffractionDataSingleCrystal::CreateCopy()const
{
   VFN_DEBUG_MESSAGE("DiffractionDataSingleCrystal::CreateCopy()",5)
   return new DiffractionDataSingleCrystal(*this);
}

const string& DiffractionDataSingleCrystal::GetClassName() const
{
   const static string className="DiffractionDataSingleCrystal";
   return className;
}

void DiffractionDataSingleCrystal::SetHklIobs(const CrystVector_long &h,
                                              const CrystVector_long &k,
                                              const CrystVector_long &l,
                                              const CrystVector_REAL &iObs,
                                              const CrystVector_REAL &sigma)
{
   VFN_DEBUG_ENTRY("DiffractionDataSingleCrystal::SetHklIobs(h,k,l,i,s)",5)
   mNbRefl=h.numElements();
   mH.resize(mNbRefl);
   mK.resize(mNbRefl);
   mL.resize(mNbRefl);
   mObsIntensity.resize(mNbRefl);
   mObsSigma.resize(mNbRefl);
   mWeight.resize(mNbRefl);
   mMultiplicity.resize(mNbRefl);
   
   mH=h;
   mK=k;
   mL=l;
   mObsIntensity=iObs;
   mObsSigma=sigma;
   mWeight=0;
   mMultiplicity=1;
   
   this->PrepareHKLarrays();
   
   this->CalcSinThetaLambda();
   
   mHasObservedData=true;
   
   // Keep a copy as squared F(hkl), to enable fourier maps
   // :TODO: stop using mObsIntensity and just keep mFhklObsSq ?
   mFhklObsSq=mObsIntensity;
   mClockGetFhklObsSq.Click();
   
   {
      char buf [200];
      sprintf(buf,"Changed HKL list, with %d reflections",(int)mNbRefl);
      (*fpObjCrystInformUser)((string)buf);
   }
   VFN_DEBUG_EXIT("DiffractionDataSingleCrystal::SetHklIobs(h,k,l,i,s)",5)
}

const CrystVector_REAL& DiffractionDataSingleCrystal::GetIcalc()const
{
   this->CalcIcalc();
   return mCalcIntensity;
}
const CrystVector_REAL& DiffractionDataSingleCrystal::GetIobs()const
{
   //if(mHasObservedData==false) DoSomething
   return mObsIntensity;
}

void DiffractionDataSingleCrystal::SetIobs(const CrystVector_REAL &obs)
{
   mObsIntensity=obs;
   // Keep a copy as squared F(hkl), to enable fourier maps
   // :TODO: stop using mObsIntensity and just keep mFhklObsSq ?
   mFhklObsSq=mObsIntensity;
   mClockGetFhklObsSq.Click();
}

const CrystVector_REAL& DiffractionDataSingleCrystal::GetSigma()const
{
   //if(mHasObservedData=false) DoSomething
   return mObsSigma;
}

void DiffractionDataSingleCrystal::SetSigma(const CrystVector_REAL& sigma) {mObsSigma=sigma;}

const CrystVector_REAL& DiffractionDataSingleCrystal::GetWeight()const
{
   //if(mHasObservedData=false) DoSomething
   return mWeight;
}

void DiffractionDataSingleCrystal::SetWeight(const CrystVector_REAL& weight)
{
   VFN_DEBUG_MESSAGE("DiffractionDataSingleCrystal::SetWeight(w)",5)
   mWeight=weight;
   mClockMaster.Click();
}

void DiffractionDataSingleCrystal::SetIobsToIcalc()
{
   VFN_DEBUG_MESSAGE("DiffractionDataSingleCrystal::SetIobsToIcalc()",5)
   mObsIntensity=this->GetIcalc();
   mObsSigma.resize(mNbRefl);
   mWeight.resize(mNbRefl);
   mWeight=1;
   mObsSigma=0;
   mHasObservedData=true;
   // Keep a copy as squared F(hkl), to enable fourier maps
   // :TODO: stop using mObsIntensity and just keep mFhklObsSq ?
   mFhklObsSq=mObsIntensity;
   mClockGetFhklObsSq.Click();
   mClockMaster.Click();
}


void DiffractionDataSingleCrystal::ImportHklIobs(const string &fileName,
                                    const long nbRefl,
                                    const int skipLines)
{
   //configure members
      mNbRefl=nbRefl;
      mH.resize(mNbRefl);
      mK.resize(mNbRefl);
      mL.resize(mNbRefl);
      mObsIntensity.resize(mNbRefl);
      mObsSigma.resize(mNbRefl);
   
   //Import data
   {   
      //:TODO: Skip the lines if required !!!
      cout << "inputing reflections from file : "+fileName<<endl;
      ifstream fin (fileName.c_str());
      if(!fin)
      {
         throw ObjCrystException("DiffractionDataSingleCrystal::ImportHklIobs() : \
Error opening file for input:"+fileName);
      }
      cout << "Number of reflections to import : " << nbRefl << endl ;
      for(long i=0;i<nbRefl;i++)
      { // :TODO: A little faster....
         //:TODO: Check for the end of the stream if(!fin.good()) {throw...}
         fin >> mH(i);
         fin >> mK(i);
         fin >> mL(i);
         fin >> mObsIntensity(i);
         mObsSigma(i)=sqrt(mObsIntensity(i));
         //cout << mObsIntensity(i) <<endl;
      }
      cout << "Finished reading data>"<< endl;
      fin.close();
   }
  //Finish
   mWeight.resize(mNbRefl);
   const REAL minIobs=mObsIntensity.max()*1e-6;
   for(int i=0;i<mNbRefl;i++) 
      if(mObsIntensity(i)<minIobs) mWeight(i)=1./minIobs;
      else mWeight(i)=1./mObsIntensity(i);
   mHasObservedData=true;
   
   mMultiplicity.resize(mNbRefl);
   mMultiplicity=1;
   
   // Keep a copy as squared F(hkl), to enable fourier maps
   // :TODO: stop using mObsIntensity and just keep mFhklObsSq ?
   mFhklObsSq=mObsIntensity;
   mClockGetFhklObsSq.Click();
   
   this->PrepareHKLarrays();
   this->SortReflectionBySinThetaOverLambda();
   {
      char buf [200];
      sprintf(buf,"Imported HKLIobs, with %d reflections",(int)mNbRefl);
      (*fpObjCrystInformUser)((string)buf);
   }
   
   cout << "Finished storing data..."<< endl ;
}

void DiffractionDataSingleCrystal::ImportHklIobsSigma(const string &fileName,
                                         const long nbRefl,
                                         const int skipLines)
{
   //configure members
      mNbRefl=nbRefl;
      mH.resize(mNbRefl);
      mK.resize(mNbRefl);
      mL.resize(mNbRefl);
      mObsIntensity.resize(mNbRefl);
      mObsSigma.resize(mNbRefl);
   
   //Import data
   {   
      cout << "inputing reflections from file : "+fileName<<endl;
      ifstream fin (fileName.c_str());
      if(!fin)
      {
         throw ObjCrystException("DiffractionDataSingleCrystal::ImportHklIobsSigma() : \
Error opening file for input:"+fileName);
      }
      if(skipLines>0)
      {//Get rid of first lines if required
         char tmpComment[200];
         for(int i=0;i<skipLines;i++) fin.getline(tmpComment,150);
      }
      cout << "Number of reflections to import : " << nbRefl << endl ;
      for(long i=0;i<nbRefl;i++)
      { // :TODO: A little faster....
         //:TODO: Check for the end of the stream if(!fin.good()) {throw...}
         fin >> mH(i);
         fin >> mK(i);
         fin >> mL(i);
         fin >> mObsIntensity(i);
         fin >> mObsSigma(i);
        //cout << mH(i)<<" "<<mK(i)<<" "<<mL(i)<<" "<<mObsIntensity(i)<<" "<<mObsSigma(i)<<endl;
      }
      cout << "Finished reading data>"<< endl;
      fin.close();
   }
  //Finish
   mWeight.resize(mNbRefl);
   const REAL minSigma=mObsSigma.max()*1e-3;
   for(int i=0;i<mNbRefl;i++) 
      if(mObsSigma(i)<minSigma) mWeight(i)=1./minSigma/minSigma;
      else mWeight(i)=1./mObsSigma(i)/mObsSigma(i);
   mHasObservedData=true;

   mMultiplicity.resize(mNbRefl);
   mMultiplicity=1;

   // Keep a copy as squared F(hkl), to enable fourier maps
   // :TODO: stop using mObsIntensity and just keep mFhklObsSq ?
   mFhklObsSq=mObsIntensity;
   mClockGetFhklObsSq.Click();
   
   this->PrepareHKLarrays();
   this->SortReflectionBySinThetaOverLambda();
   {
      char buf [200];
      sprintf(buf,"Imported HKLIobsSigma, with %d reflections",(int)mNbRefl);
      (*fpObjCrystInformUser)((string)buf);
   }

   cout << "Finished storing data..."<< endl ;

}

void DiffractionDataSingleCrystal::ImportHklIobsSigmaJanaM91(const string &fileName)
{
   //configure members
      mNbRefl=1000;//reasonable beginning value ?
      mH.resize(mNbRefl);
      mK.resize(mNbRefl);
      mL.resize(mNbRefl);
      mObsIntensity.resize(mNbRefl);
      mObsSigma.resize(mNbRefl);
   
   //Import data
   {   
      cout << "inputing reflections from Jana98 file : "+fileName<<endl;
      ifstream fin (fileName.c_str());
      if(!fin)
      {
         throw ObjCrystException("DiffractionDataSingleCrystal::ImportHklIobsSigmaJanaM91() : \
Error opening file for input:"+fileName);
      }
      long i=0;
      REAL tmpH;
      REAL junk;
      fin >> tmpH;
      while(tmpH != 999)
      { // :TODO: A little faster....
         //:TODO: Check for the end of the stream if(!fin.good()) {throw...}
         i++;
         if(i>=mNbRefl)
         {
            cout << mNbRefl << " reflections imported..." << endl;
            mNbRefl+=1000;
            mH.resizeAndPreserve(mNbRefl);
            mK.resizeAndPreserve(mNbRefl);
            mL.resizeAndPreserve(mNbRefl);
            mObsIntensity.resizeAndPreserve(mNbRefl);
            mObsSigma.resizeAndPreserve(mNbRefl);
         }
         mH(i)=tmpH;
         fin >> mK(i);
         fin >> mL(i);
         fin >> mObsIntensity(i);
         fin >> mObsSigma(i);
         fin >> junk;
         fin >> junk;
         fin >> junk;
         fin >> junk;
         fin >> tmpH;
      }
      fin.close();
      mNbRefl=i;
      cout << mNbRefl << " reflections imported." << endl;
      mH.resizeAndPreserve(mNbRefl);
      mK.resizeAndPreserve(mNbRefl);
      mL.resizeAndPreserve(mNbRefl);
      mObsIntensity.resizeAndPreserve(mNbRefl);
      mObsSigma.resizeAndPreserve(mNbRefl);
   }
  //Finish
   mWeight.resize(mNbRefl);
   mWeight=1;
   
   mMultiplicity.resize(mNbRefl);
   mMultiplicity=1;

   // Keep a copy as squared F(hkl), to enable fourier maps
   // :TODO: stop using mObsIntensity and just keep mFhklObsSq ?
   mFhklObsSq=mObsIntensity;
   mClockGetFhklObsSq.Click();
   
   this->PrepareHKLarrays();
   this->SortReflectionBySinThetaOverLambda();
   cout << "Finished storing data..."<< endl ;

   mHasObservedData=true;
   {
      char buf [200];
      sprintf(buf,"Imported HKLIobsSigma from Jana, with %d reflections",(int)mNbRefl);
      (*fpObjCrystInformUser)((string)buf);
   }
}

void DiffractionDataSingleCrystal::ImportHklIobsGroup(const string &fileName,const unsigned int skipLines)
{
   //configure members
      mNbRefl=0;
      mNbGroup=0;
      mH.resize(500);
      mK.resize(500);
      mL.resize(500);
      mObsIntensity.resize(500);
      mObsSigma.resize(500);
      mGroupIndex.resize(500);
      mGroupIobs.resize(500);
      mGroupSigma.resize(500);
      mGroupWeight.resize(500);
   //Import data
   {   
      //:TODO: Skip the lines if required !!!
      cout << "inputing reflections from file : "+fileName<<endl;
      ifstream fin (fileName.c_str());
      if(!fin)
      {
         throw ObjCrystException("DiffractionDataSingleCrystal::ImportHklIobs() : \
Error opening file for input:"+fileName);
      }
      string buffer;
      REAL h,k,l,iobs,sigma;
      while(true)
      {
         getline(fin,buffer);
         const int n=sscanf(buffer.c_str(),"%f %f %f %f %f",&h,&k,&l,&iobs,&sigma);
         if(n<3) break;
         mH(mNbRefl)=h;
         mK(mNbRefl)=k;
         mL(mNbRefl)=l;
         mGroupIndex(mNbRefl)=mNbGroup;
         //cout<<mNbRefl<<" "<<h<<" "<<k<<" "<<l<<"(g="<<mNbGroup<<") n="<<n;
         if(n>=4)
         {
            //cout<<" Iobs="<<iobs;
            mObsIntensity(mNbRefl)=iobs;
            if(n!=5)
               sigma=sqrt(fabs(iobs)+1e-6);
            mObsSigma(mNbRefl)=sigma;
            mGroupIobs(mNbGroup)=iobs;
            mGroupSigma(mNbGroup)=sigma;
            mGroupWeight(mNbGroup)=1/(sigma*sigma+1e-6);
            mNbGroup++;
            if(mNbGroup==mGroupIobs.numElements())
            {
               mGroupIobs.resizeAndPreserve(mNbGroup+500);
               mGroupSigma.resizeAndPreserve(mNbGroup+500);
               mGroupWeight.resizeAndPreserve(mNbGroup+500);
            }
         }
         else
         {
            mObsIntensity(mNbRefl)=0;
            mObsSigma(mNbRefl)=0;
         }
         //cout<<endl;
         mNbRefl++;
         if(mNbRefl==mH.numElements())
         {
            mH.resizeAndPreserve(mNbRefl+500);
            mK.resizeAndPreserve(mNbRefl+500);
            mL.resizeAndPreserve(mNbRefl+500);
            mObsIntensity.resizeAndPreserve(mNbRefl+500);
            mObsSigma.resizeAndPreserve(mNbRefl+500);
            mGroupIndex.resizeAndPreserve(mNbRefl+500);
         }
         if(fin.eof()) break;
      }
      fin.close();
   }
   // This must NOT be changed with this kind of data.
   mGroupOption.SetChoice(2);
   // So de-register the option so that it is hidden from the user's view
   mOptionRegistry.DeRegister(mGroupOption);
   mClockMaster.RemoveChild(mGroupOption.GetClock());
   //Finish
   mH.resizeAndPreserve(mNbRefl);
   mK.resizeAndPreserve(mNbRefl);
   mL.resizeAndPreserve(mNbRefl);
   mObsIntensity.resizeAndPreserve(mNbRefl);
   mObsSigma.resizeAndPreserve(mNbRefl);
   mWeight.resize(mNbRefl);
   mGroupIndex.resizeAndPreserve(mNbRefl);// this will change after sorting reflections
   mGroupIobs.resizeAndPreserve(mNbGroup);
   mGroupWeight.resizeAndPreserve(mNbGroup);
   mGroupSigma.resizeAndPreserve(mNbGroup);

   const REAL minIobs=mObsIntensity.max()*1e-6;
   for(int i=0;i<mNbRefl;i++) 
      if(mObsIntensity(i)<minIobs) mWeight(i)=1./minIobs;
      else mWeight(i)=1./mObsIntensity(i);
   mHasObservedData=true;
   
   mMultiplicity.resize(mNbRefl);
   mMultiplicity=1;
   
   this->PrepareHKLarrays();
   {
      char buf [200];
      sprintf(buf,"Imported HKLIobs, with %d reflections",(int)mNbRefl);
      (*fpObjCrystInformUser)((string)buf);
   }
   this->SortReflectionBySinThetaOverLambda();
   this->CalcIcalc();
}


REAL DiffractionDataSingleCrystal::GetRw()const
{
   TAU_PROFILE("DiffractionData::Rw()"," REAL()",TAU_DEFAULT);
   VFN_DEBUG_MESSAGE("DiffractionData::Rw()",3);
   if(mHasObservedData==false)
   {
      return 0;
   }
   REAL tmp1=0;
   REAL tmp2=0;
   const REAL *p1;
   const REAL *p2;
   const REAL *p3;
   long nb;
   if(mGroupOption.GetChoice()==0)
   {
      p1=mCalcIntensity.data();
      p2=mObsIntensity.data();
      p3=mWeight.data();
      nb=mNbReflUsed;
   }
   else
   {
      p1=mGroupIcalc.data();
      p2=mGroupIobs.data();
      p3=mGroupWeight.data();
      nb=mGroupIobs.numElements();
   }
   
   for(long i=nb;i>0;i--)
   {
      tmp1 += *p3 * ( *p1 - *p2) * ( *p1 - *p2);
      tmp2 += *p3 * *p2 * *p2;
      p1++;p2++;p3++;
   }
   tmp1=sqrt(tmp1/tmp2);
   VFN_DEBUG_MESSAGE("DiffractionData::Rw()="<<tmp1,3);
   return tmp1 ;// /this->GetNbRefl();
}

REAL DiffractionDataSingleCrystal::GetR()const
{
   TAU_PROFILE("DiffractionData::R()"," REAL()",TAU_DEFAULT);
   VFN_DEBUG_MESSAGE("DiffractionData::R()",3);
   if(mHasObservedData==false)
   {
      return 0;
   }
   
   REAL tmp1=0;
   REAL tmp2=0;
   const REAL *p1;
   const REAL *p2;
   long nb;
   if(mGroupOption.GetChoice()==0)
   {
      p1=mCalcIntensity.data();
      p2=mObsIntensity.data();
      nb=mNbReflUsed;
   }
   else
   {
      p1=mGroupIcalc.data();
      p2=mGroupIobs.data();
      nb=mGroupIobs.numElements();
   }
   
   for(long i=nb;i>0;i--)
   {
      tmp1 += ( *p1 - *p2) * ( *p1 - *p2);
      tmp2 += *p2 * *p2;
      p1++;p2++;
   }
   tmp1=sqrt(tmp1/tmp2);
   
   VFN_DEBUG_MESSAGE("DiffractionData::R()="<<tmp1,3);
   return tmp1 ;
}

REAL DiffractionDataSingleCrystal::GetChi2()const
{
   if(mHasObservedData==false)
   {
      mChi2=0;
      return mChi2;
   }
   this->GetNbReflBelowMaxSinThetaOvLambda();
   if(mClockChi2>mClockMaster) return mChi2;
   
   this->CalcIcalc();
   if(mClockChi2>mClockIcalc) return mChi2;
   
   TAU_PROFILE("DiffractionData::Chi2()"," REAL()",TAU_DEFAULT);
   VFN_DEBUG_ENTRY("DiffractionData::Chi2()",3);
   this->FitScaleFactorForRw();
   mChi2=0;
   const REAL *p1;
   const REAL *p2;
   const REAL *p3;
   long nb;
   if(mGroupOption.GetChoice()==0)
   {
      p1=mCalcIntensity.data();
      p2=mObsIntensity.data();
      p3=mWeight.data();
      nb=mNbReflUsed;
   }
   else
   {
      p1=mGroupIcalc.data();
      p2=mGroupIobs.data();
      p3=mGroupWeight.data();
      nb=mGroupIobs.numElements();
   }
   
   for(long i=nb;i>0;i--)
   {
      mChi2 += *p3++ * ( *p1 - *p2) * ( *p1 - *p2);
      p1++;p2++;
   }
   /*
   // SSE code gives about 30% faster code on P3 (tested with 10000 reflections),
   // but scarcely any improvement on athlon-xp
   union sse4f
   {
      __m128 m128;
      struct
      {
         float x, y, z, w;
      };
   } ;
   const long nb=mNbReflUsed/4;
   const float* p1=mCalcIntensity.data();
   const float* p2=mObsIntensity.data();
   const float* p3=mWeight.data();
   for(long i=0;i<nb;i++)
   {
      //segfaults with _mm_load_ps() instead of _mm_loadu_ps? Alignment problem ?
      __m128 a = _mm_sub_ps(_mm_loadu_ps(p1),_mm_loadu_ps(p2));
      union sse4f b;
      b.m128= _mm_mul_ps(_mm_loadu_ps(p3),_mm_mul_ps(a,a));
      mChi2 += b.x+b.y+b.z+b.w;
      p1+=4;p2+=4;p3+=4;
   }
   p1-=4;p2-=4;p3-=4;
   for(long i=nb*4;i<mNbReflUsed;i++)
   {
      mChi2 += *p3++ * ( *p1 - *p2) * ( *p1 - *p2);
      p1++;p2++;
   }
   */
   
   mClockChi2.Click();
   VFN_DEBUG_EXIT("DiffractionData::Chi2()="<<mChi2,3);
   return mChi2;
}

void DiffractionDataSingleCrystal::FitScaleFactorForRw() const
{
   TAU_PROFILE("DiffractionData::FitScaleFactorForRw()","void ()",TAU_DEFAULT);
   VFN_DEBUG_MESSAGE("DiffractionData::FitScaleFactorForRw()",3);
   if(mHasObservedData==false)
   {//throw exception here ?
      return;
      //throw ObjCrystException("DiffractionData::FitScaleFactorForRw() Cannot compute Rw 
      //   or scale factor: there is no observed data !");
   }
   REAL tmp1=0;
   REAL tmp2=0;
   const REAL *p1;
   const REAL *p2;
   const REAL *p3;
   long nb;
   if(mGroupOption.GetChoice()==0)
   {
      p1=mCalcIntensity.data();
      p2=mObsIntensity.data();
      p3=mWeight.data();
      nb=mNbReflUsed;
   }
   else
   {
      p1=mGroupIcalc.data();
      p2=mGroupIobs.data();
      p3=mGroupWeight.data();
      nb=mGroupIobs.numElements();
   }
   
   for(long i=nb;i>0;i--)
   {
      tmp1 += *p3 * (*p1) * (*p2++);
      tmp2 += *p3++ * (*p1) * (*p1);
      p1++;
   }
   mScaleFactor      *= tmp1/tmp2;
   mClockScaleFactor.Click();

   mCalcIntensity *= tmp1/tmp2;
   if(0!=mGroupOption.GetChoice()) mGroupIcalc*= tmp1/tmp2;
   mClockIcalc.Click();
}

void DiffractionDataSingleCrystal::FitScaleFactorForR() const
{
   TAU_PROFILE("DiffractionData::FitScaleFactorForR()","void ()",TAU_DEFAULT);
   VFN_DEBUG_MESSAGE("DiffractionData::FitScaleFactorForR()",3);
   if(mHasObservedData==false)
   {//throw exception here ?
      return;
      //throw ObjCrystException("DiffractionData::FitScaleFactorForR() Cannot compute R 
      //   or scale factor: there is no observed data !");
   }
   REAL tmp1=0;
   REAL tmp2=0;
   const REAL *p1;
   const REAL *p2;
   long nb;
   if(mGroupOption.GetChoice()==0)
   {
      p1=mCalcIntensity.data();
      p2=mObsIntensity.data();
      nb=mNbReflUsed;
   }
   else
   {
      p1=mGroupIcalc.data();
      p2=mGroupIobs.data();
      nb=mGroupIobs.numElements();
   }
   
   for(long i=nb;i>0;i--)
   {
      tmp1 += (*p1) * (*p2++);
      tmp2 += (*p1) * (*p1);
      p1++;
   }
   mScaleFactor      *= tmp1/tmp2;
   mClockScaleFactor.Click();
   
   mCalcIntensity *= tmp1/tmp2;
   if(0!=mGroupOption.GetChoice()) mGroupIcalc*= tmp1/tmp2;
   mClockIcalc.Click();
}

REAL DiffractionDataSingleCrystal::GetBestRFactor() const
{
   TAU_PROFILE("DiffractionData::GetBestRFactor()","void ()",TAU_DEFAULT);
   VFN_DEBUG_MESSAGE("DiffractionData::GetBestRFactor()",3);
   this->FitScaleFactorForR();
   return this->GetR();
}

void DiffractionDataSingleCrystal::SetSigmaToSqrtIobs()
{
   for(long i=0;i<mObsIntensity.numElements();i++) mObsSigma(i)=sqrt(fabs(mObsIntensity(i)));
   if(1==mGroupOption.GetChoice()) mClockPrepareTwinningCorr.Reset();
   // This is not needed for mGroupOption==2
}

void DiffractionDataSingleCrystal::SetWeightToInvSigma2(const REAL minRelatSigma)
{
   //:KLUDGE: If less than 1e-6*max, set to 0.... Do not give weight to unobserved points
   const REAL min=MaxAbs(mObsSigma)*minRelatSigma;
   for(long i=0;i<mObsSigma.numElements();i++)
   {
      if(mObsSigma(i)<min) mWeight(i)=0 ; else  mWeight(i) =1./mObsSigma(i)/mObsSigma(i);
   }
   if(1==mGroupOption.GetChoice()) mClockPrepareTwinningCorr.Reset();
   // This is not needed for mGroupOption==2
}

REAL DiffractionDataSingleCrystal::GetScaleFactor()const {return mScaleFactor;}

//void DiffractionDataSingleCrystal::SetScaleFactor(const REAL s) {mScaleFactor=s;}

void DiffractionDataSingleCrystal::PrintObsData()const
{
   this->CalcSinThetaLambda();
   cout << "DiffractionData : " << mName <<endl;
   cout << "Number of observed reflections : " << mNbRefl << endl;
   cout << "       H        K        L     Iobs        Sigma       sin(theta)/lambda)" <<endl;
   cout << mH.numElements()<<endl;
   cout << mK.numElements()<<endl;
   cout << mL.numElements()<<endl;
   cout << mObsIntensity.numElements()<<endl;
   cout << mObsSigma.numElements()<<endl;
   cout << mSinThetaLambda.numElements()<<endl;
   cout << FormatVertVectorHKLFloats<REAL>
               (mH,mK,mL,mObsIntensity,mObsSigma,mSinThetaLambda,12,4);
}

void DiffractionDataSingleCrystal::PrintObsCalcData()const
{
   this->CalcIcalc();
   CrystVector_REAL tmpTheta=mTheta;
   tmpTheta*= RAD2DEG;
   /*
   CrystVector_REAL tmp=mObsIntensity;
   CrystVector_REAL tmpS=mObsSigma;
   if(true==mHasObservedData)
   {
      cout << "Scale factor : " << mScaleFactor <<endl;
      tmp*= mScaleFactor*mScaleFactor;//the scale factor is computed for F's
      tmpS*=mScaleFactor;
   } else cout << "No observed data. Default scale factor =1." <<endl;
   */
   cout << "DiffractionData : " << mName <<endl;
   cout << " Scale Factor : " << mScaleFactor <<endl;
   cout << "Number of observed reflections : " << mNbRefl << endl;
   
   cout << "       H        K        L     Iobs        Sigma       Icalc  ";
   cout << "      multiplicity     Theta      SiThSL       Re(F)     Im(F)    Weight" <<endl;
   cout << FormatVertVectorHKLFloats<REAL>(mH,mK,mL,
               mObsIntensity,mObsSigma,mCalcIntensity,
               mMultiplicity,tmpTheta,mSinThetaLambda,
               mFhklCalcReal,mFhklCalcImag,mWeight,12,4);
}

void DiffractionDataSingleCrystal::SetUseOnlyLowAngleData(
                     const bool useOnlyLowAngle,const REAL angle)
{
   throw ObjCrystException("DiffractionDataSingleCrystal::SetUseOnlyLowAngleData() :\
 not yet implemented for DiffractionDataSingleCrystal.");
}

void DiffractionDataSingleCrystal::SaveHKLIobsIcalc(const string &filename)
{
   VFN_DEBUG_MESSAGE("DiffractionDataSingleCrystal::SaveHKLIobsIcalc",5)
   this->GetIcalc();
   ofstream out(filename.c_str());
   CrystVector_REAL theta;
   theta=mTheta;
   theta *= RAD2DEG;
   
   if(false == mHasObservedData)
   {
      out << "#    H        K        L      Icalc    theta  sin(theta)/lambda"
          <<"  Re(F)   Im(F)" << endl;
      out << FormatVertVectorHKLFloats<REAL>(mH,mK,mL,mCalcIntensity,
                           theta,mSinThetaLambda,mFhklCalcReal,mFhklCalcImag,12,4);
   }
   else
   {
      out << "#    H        K        L      Iobs   Icalc    theta"
          <<" sin(theta)/lambda  Re(F)   Im(F)" << endl;
      out << FormatVertVectorHKLFloats<REAL>(mH,mK,mL,mObsIntensity,mCalcIntensity,
                           theta,mSinThetaLambda,mFhklCalcReal,mFhklCalcImag,12,4);
   }
   out.close();
   VFN_DEBUG_MESSAGE("DiffractionDataSingleCrystal::SaveHKLIobsIcalc:End",3)
}

void DiffractionDataSingleCrystal::GlobalOptRandomMove(const REAL mutationAmplitude,
                         const RefParType *type)
{
   this->RefinableObj::GlobalOptRandomMove(mutationAmplitude,type);
   //this->FitScaleFactorForRw();
}

REAL DiffractionDataSingleCrystal::GetLogLikelihood()const
{
   return this->GetChi2();
}

void DiffractionDataSingleCrystal::InitRefParList()
{
   VFN_DEBUG_MESSAGE("DiffractionDataSingleCrystal::InitRefParList()",5)
   RefinablePar tmp("Scale factor",&mScaleFactor,
                     1e-10,1e10,gpRefParTypeScattDataScale,REFPAR_DERIV_STEP_RELATIVE,
                     false,true,true,false,1.);
   tmp.SetGlobalOptimStep(0.);
   tmp.AssignClock(mClockScaleFactor);
   tmp.SetDerivStep(1e-4);
   this->AddPar(tmp);
}
unsigned int DiffractionDataSingleCrystal::GetNbLSQFunction()const{return 1;}
const CrystVector_REAL& 
   DiffractionDataSingleCrystal::GetLSQCalc(const unsigned int) const
{return this->GetIcalc();}

const CrystVector_REAL& 
   DiffractionDataSingleCrystal::GetLSQObs(const unsigned int) const
{return this->GetIobs();}

const CrystVector_REAL& 
   DiffractionDataSingleCrystal::GetLSQWeight(const unsigned int) const
{return this->GetWeight();}

const Radiation& DiffractionDataSingleCrystal::GetRadiation()const { return mRadiation;}
Radiation& DiffractionDataSingleCrystal::GetRadiation() { return mRadiation;}
void DiffractionDataSingleCrystal::SetRadiationType(const RadiationType radiation)
{
   VFN_DEBUG_MESSAGE("DiffractionDataSingleCrystal::SetRadiationType():End",5)
   mRadiation.SetRadiationType(radiation);
}

void DiffractionDataSingleCrystal::SetWavelength(const REAL lambda)
{
   VFN_DEBUG_MESSAGE("DiffractionDataSingleCrystal::SetWavelength() to "<<lambda,5)
   this->GetRadiation().SetWavelength(lambda);
}

void DiffractionDataSingleCrystal::SetWavelength(const string &XRayTubeElementName,
                                   const REAL alpha2Alpha2ratio)
{
   VFN_DEBUG_MESSAGE("DiffractionDataSingleCrystal::SetWavelength() to "<<XRayTubeElementName,5)
   this->GetRadiation().SetWavelength(XRayTubeElementName,alpha2Alpha2ratio);
}

void DiffractionDataSingleCrystal::SetEnergy(const REAL energy)
{
   this->SetWavelength(12398.4/energy);
}

void DiffractionDataSingleCrystal::CalcIcalc() const
{
   TAU_PROFILE("DiffractionData::CalcIcalc()","void ()",TAU_DEFAULT);
   VFN_DEBUG_MESSAGE("DiffractionData::CalcIcalc():"<<this->GetName(),3)
   this->GetFhklCalcSq();
   if( (mClockStructFactorSq<mClockIcalc) && (mClockScaleFactor<mClockIcalc)
        && ((0==mGroupOption.GetChoice()) || (mClockPrepareTwinningCorr<mClockIcalc)) ) return;
   
   mCalcIntensity=mFhklCalcSq;
   mCalcIntensity*=mScaleFactor;
   if(0!=mGroupOption.GetChoice())
   {
      if(1==mGroupOption.GetChoice()) this->PrepareTwinningCalc();
      mGroupIcalc.resize(mNbGroup);
      mGroupIcalc=0;
      long first=0;
      for(long i=0;i<mNbGroup;i++)
      {
         //cout<<"Group #"<<i<<":"<<first<<"->"<<mGroupIndex(i)<<endl;
         for(long j=first;j<mGroupIndex(i);j++)
         {
            //cout <<"       " << mIntH(j)<<" "<< mIntK(j)<<" "<< mIntL(j)<<" " << mCalcIntensity(j)<<endl;
            mGroupIcalc(i)+=mCalcIntensity(j);
         }
         first=mGroupIndex(i);
         //cout  << "   => Icalc="<< mGroupIcalc(i) <<" , Iobs="<< mGroupIobs(i)<<endl<<endl;
      }
   }
   mClockIcalc.Click();
}

CrystVector_long DiffractionDataSingleCrystal::SortReflectionBySinThetaOverLambda(const REAL maxSTOL)
{
   TAU_PROFILE("DiffractionDataSingleCrystal::SortReflectionBySinThetaOverLambda()","void ()",TAU_DEFAULT);
   VFN_DEBUG_ENTRY("DiffractionDataSingleCrystal::SortReflectionBySinThetaOverLambda()",5)
   // ScatteringData::SortReflectionBySinThetaOverLambda only sorts H,K,L and multiplicity.
   CrystVector_long index=this->ScatteringData::SortReflectionBySinThetaOverLambda(maxSTOL);
   
   if(mObsIntensity.numElements()==mNbRefl)
   {
      CrystVector_REAL tmpObs,tmpSigma,tmpWeight;
      tmpObs=mObsIntensity;
      tmpSigma=mObsSigma;
      tmpWeight=mWeight;
      for(long i=0;i<mNbRefl;i++)
      {
         mObsIntensity(i)=tmpObs(index(i));
         mObsSigma(i)=tmpSigma(index(i));
         mWeight(i)=tmpWeight(index(i));
      }
      // Keep a copy as squared F(hkl), to enable fourier maps
      // :TODO: stop using mObsIntensity and just keep mFhklObsSq ?
      mFhklObsSq=mObsIntensity;
      mClockGetFhklObsSq.Click();

      if(mGroupOption.GetChoice()==2)
      {
         CrystVector_long tmp;//,oldgroup(mNbGroup);;
         tmp=mGroupIndex;
         for(long i=0;i<mNbRefl;i++) mGroupIndex(i)=tmp(index(i));
         /*
         for(long i=0;i<mNbRefl;i++)
            cout<<mIntH(i)<<" "<<mIntK(i)<<" "<<mIntL(i)<<" "
                <<mObsIntensity(i)<<" "<<mObsSigma(i)<<" "<<mWeight(i)<< ":"<<mGroupIndex(i)<<endl;
         */
         // Now re-index the groups of reflections in
         // ascending order
         {
            index.resize(mNbGroup);
            index=-1;
            long group=0;
            CrystVector_REAL oldGroupIobs, oldGroupWeight,oldGroupSigma;
            oldGroupIobs  =mGroupIobs;
            oldGroupWeight=mGroupWeight;
            oldGroupSigma=mGroupSigma;
            for(long i=0;i<mNbRefl;i++)
            {
               if(index(mGroupIndex(i))==-1)// first reflection of a group ?
               {
                  mGroupIobs(group)=oldGroupIobs(mGroupIndex(i));
                  mGroupSigma(group)=oldGroupSigma(mGroupIndex(i));
                  mGroupWeight(group)=oldGroupWeight(mGroupIndex(i));
                  //oldgroup(group)=mGroupIndex(i);
                  index(mGroupIndex(i))=group++;
               }
               mGroupIndex(i)=index(mGroupIndex(i));
            }
         }
         /*
         cout<<mIntH.numElements()<<","
             <<mIntK.numElements()<<","
             <<mIntL.numElements()<<","
             <<oldgroup.numElements()<<","<<mGroupIndex.numElements()<<endl;
         for(long i=0;i<mNbRefl;i++)
            cout<<mIntH(i)<<" "<<mIntK(i)<<" "<<mIntL(i)<<endl
                <<"             :"<<oldgroup(mGroupIndex(i))<<"->"<<mGroupIndex(i)<<endl;
         */
         // Now re-group the reflections
         index=SortSubs(mGroupIndex);
         {
            CrystVector_long oldH,oldK,oldL,oldMult;
            oldH=mH;
            oldK=mK;
            oldL=mL;
            oldMult=mMultiplicity;
            tmpObs=mObsIntensity;
            tmpSigma=mObsSigma;
            tmpWeight=mWeight;
            for(long i=0;i<mNbRefl;i++)
            {
               const long subs=index(i);
               mH(i)=oldH(subs);
               mK(i)=oldK(subs);
               mL(i)=oldL(subs);
               mMultiplicity(i)=oldMult(subs);
               mObsIntensity(i)=tmpObs(subs);
               mObsSigma(i)=tmpSigma(subs);
               mWeight(i)=tmpWeight(subs);
            }
            mClockHKL.Click();
            this->PrepareHKLarrays();
            this->CalcSinThetaLambda();
         }
         
         // re-write mGroupIndex so that it marks the
         // last reflection of each group.
         index=mGroupIndex;
         mGroupIndex.resize(mNbGroup);
         long group=0;
         for(long i=0;i<mNbRefl;i++)
            if(index(i)!=group)
               mGroupIndex(group++)=i;
         mGroupIndex(mNbGroup-1)=mNbRefl;
      }
   }
   else
   {// if there are no observed values, enter dummy ones
      mObsIntensity.resize(mNbRefl);
      mObsSigma.resize(mNbRefl);
      mWeight.resize(mNbRefl);
      mObsIntensity=100.;
      mObsSigma=1.;
      mWeight=1.;
      mFhklObsSq.resize(0);
      mClockGetFhklObsSq.Click();
   }
   VFN_DEBUG_EXIT("DiffractionDataSingleCrystal::SortReflectionBySinThetaOverLambda()",5)
   return index;
}

void DiffractionDataSingleCrystal::InitOptions()
{
   static string GroupOption;
   static string GroupOptionChoices[3];
   static bool needInitNames=true;
   if(true==needInitNames)
   {
      GroupOption="Group Reflections";
      GroupOptionChoices[0]="No";
      GroupOptionChoices[1]="Sum equally-spaced reflections";
      GroupOptionChoices[2]="Sum according to user data";
      needInitNames=false;
   }
   mGroupOption.Init(3,&GroupOption,GroupOptionChoices);
   mGroupOption.SetChoice(0);
   this->AddOption(&mGroupOption);
}

void DiffractionDataSingleCrystal::PrepareTwinningCalc() const
{
   if(mClockPrepareTwinningCorr>mClockHKL) return;
   VFN_DEBUG_ENTRY("DiffractionDataSingleCrystal::PrepareTwinningCalc()",5)
   // first get the index of reflections which limit each block of summed reflections
   mNbGroup=0;
   {
      const REAL dSiThOvLa=.0001;
      mGroupIndex.resize(mNbReflUsed);
      this->CalcSinThetaLambda();
      REAL sithol0=mSinThetaLambda(0)+dSiThOvLa;
      for(long i=1;i<mNbReflUsed;i++)
      {
         if(mSinThetaLambda(i)>sithol0)
         {
            mGroupIndex(mNbGroup++)=i;
            sithol0=mSinThetaLambda(i)+dSiThOvLa;
         }
      }
      mGroupIndex(mNbGroup++)=mNbReflUsed;
      mGroupIndex.resizeAndPreserve(mNbGroup);
   }
   // Calculate summed Iobs and weight
   {
      mGroupIobs.resize(mNbGroup);
      mGroupIobs=0;
      mGroupWeight.resize(mNbGroup);
      mGroupWeight=0;
      mGroupSigma.resize(mNbGroup);
      mGroupSigma=0;
      long first=0;
      for(long i=0;i<mNbGroup;i++)
      {
         for(long j=first;j<mGroupIndex(i);j++) 
         {
            mGroupIobs(i)+=mObsIntensity(j);
            mGroupSigma(i)+=mObsSigma(j)*mObsSigma(j);
            mGroupWeight(i)+=mObsSigma(j)*mObsSigma(j);
         }
         mGroupWeight(i)=1./mGroupWeight(i);
         mGroupSigma(i)=sqrt(mGroupSigma(i));
         first=mGroupIndex(i);
      }
   }
   mClockPrepareTwinningCorr.Click();
   VFN_DEBUG_EXIT("DiffractionDataSingleCrystal::PrepareTwinningCalc()",5)
}

#ifdef __WX__CRYST__
WXCrystObjBasic* DiffractionDataSingleCrystal::WXCreate(wxWindow* parent)
{
   //:TODO: Check mpWXCrystObj==0
   mpWXCrystObj=new WXDiffractionSingleCrystal(parent,this);
   return mpWXCrystObj;
}
#endif

}//namespace
