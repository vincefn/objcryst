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

namespace ObjCryst
{
//######################################################################
//    DiffractionDataSingleCrystal
//######################################################################
ObjRegistry<DiffractionDataSingleCrystal> 
   gDiffractionDataSingleCrystalRegistry("Global DiffractionDataSingleCrystal Registry");

DiffractionDataSingleCrystal::DiffractionDataSingleCrystal():mScaleFactor(1.)
{
   VFN_DEBUG_MESSAGE("DiffractionDataSingleCrystal::DiffractionDataSingleCrystal()",5)
   this->InitRefParList();
   this->InitOptions();
   gDiffractionDataSingleCrystalRegistry.Register(*this);
   gTopRefinableObjRegistry.Register(*this);
}
DiffractionDataSingleCrystal::DiffractionDataSingleCrystal(Crystal &cryst):mScaleFactor(1.)
{
   VFN_DEBUG_MESSAGE("DiffractionDataSingleCrystal::DiffractionDataSingleCrystal()",5)
   this->InitRefParList();
   this->SetCrystal(cryst);
   this->InitOptions();
   gDiffractionDataSingleCrystalRegistry.Register(*this);
   gTopRefinableObjRegistry.Register(*this);
}

DiffractionDataSingleCrystal::DiffractionDataSingleCrystal(const DiffractionDataSingleCrystal &old):
ScatteringData(old),
mHasObservedData(old.mHasObservedData),mRadiation(old.mRadiation)
{
   mObsIntensity=old.mObsIntensity;
   mObsSigma=old.mObsSigma;
   mWeight=old.mWeight;
   mCalcIntensity=old.mCalcIntensity;
   mScaleFactor=old.mScaleFactor;
   this->InitOptions();
   mTwinningOption.SetChoice(old.mTwinningOption.GetChoice());
   gDiffractionDataSingleCrystalRegistry.Register(*this);
   gTopRefinableObjRegistry.Register(*this);
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

void DiffractionDataSingleCrystal::SetIobs(const CrystVector_REAL &obs) {mObsIntensity=obs;}

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
   
   this->PrepareHKLarrays();
   this->SortReflectionByTheta();
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

   this->PrepareHKLarrays();
   this->SortReflectionByTheta();
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

   this->PrepareHKLarrays();
   this->SortReflectionByTheta();
   cout << "Finished storing data..."<< endl ;

   mHasObservedData=true;
   {
      char buf [200];
      sprintf(buf,"Imported HKLIobsSigma from Jana, with %d reflections",(int)mNbRefl);
      (*fpObjCrystInformUser)((string)buf);
   }
}
REAL DiffractionDataSingleCrystal::GetRw()const
{
   TAU_PROFILE("DiffractionData::Rw()"," REAL()",TAU_DEFAULT);
   VFN_DEBUG_MESSAGE("DiffractionData::Rw()",3);
   if(mHasObservedData==false)
   {
      throw ObjCrystException("DiffractionData::Rw() Cannot compute Rw ! \
         There is no observed data !");
   }
   REAL tmp1=0;
   REAL tmp2=0;
   const REAL *p1;
   const REAL *p2;
   const REAL *p3;
   long nb;
   switch(mTwinningOption.GetChoice())
   {
      case 0:
      {
         p1=mCalcIntensity.data();
         p2=mObsIntensity.data();
         p3=mWeight.data();
         nb=mNbReflUsed;
         break;
      }
      case 1:
      {
         p1=mTwinnedIcalcSum.data();
         p2=mTwinnedIobsSum.data();
         p3=mTwinnedWeight.data();
         nb=mTwinnedIobsSum.numElements();
         break;
      }
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
      throw ObjCrystException("DiffractionData::R() Cannot compute R ! \
         There is no observed data !");
   }
   
   REAL tmp1=0;
   REAL tmp2=0;
   const REAL *p1;
   const REAL *p2;
   long nb;
   switch(mTwinningOption.GetChoice())
   {
      case 0:
      {
         p1=mCalcIntensity.data();
         p2=mObsIntensity.data();
         nb=mNbReflUsed;
         break;
      }
      case 1:
      {
         p1=mTwinnedIcalcSum.data();
         p2=mTwinnedIobsSum.data();
         nb=mTwinnedIobsSum.numElements();
         break;
      }
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
   TAU_PROFILE("DiffractionData::Chi2()"," REAL()",TAU_DEFAULT);
   VFN_DEBUG_MESSAGE("DiffractionData::Chi2()",3);
   if(mHasObservedData==false)
   {
      throw ObjCrystException("DiffractionData::Chi2() Cannot compute Chi^2 ! \
         There is no observed data !");
   }
   REAL tmp1=0;
   const REAL *p1;
   const REAL *p2;
   const REAL *p3;
   long nb;
   switch(mTwinningOption.GetChoice())
   {
      case 0:
      {
         p1=mCalcIntensity.data();
         p2=mObsIntensity.data();
         p3=mWeight.data();
         nb=mNbReflUsed;
         break;
      }
      case 1:
      {
         p1=mTwinnedIcalcSum.data();
         p2=mTwinnedIobsSum.data();
         p3=mTwinnedWeight.data();
         nb=mTwinnedIobsSum.numElements();
         break;
      }
   }
   for(long i=nb;i>0;i--)
   {
      tmp1 += *p3++ * ( *p1 - *p2) * ( *p1 - *p2);
      p1++;p2++;
   }
   VFN_DEBUG_MESSAGE("DiffractionData::Chi2()="<<tmp1,3);
   return tmp1;
}

void DiffractionDataSingleCrystal::FitScaleFactorForRw()
{
   TAU_PROFILE("DiffractionData::FitScaleFactorForRw()","void ()",TAU_DEFAULT);
   VFN_DEBUG_MESSAGE("DiffractionData::FitScaleFactorForRw()",3);
   if(mHasObservedData==false)
   {//throw exception here ?
      throw ObjCrystException("DiffractionData::FitScaleFactorForRw() Cannot compute Rw \
         or scale factor: there is no observed data !");
   }
   REAL tmp1=0;
   REAL tmp2=0;
   const REAL *p1;
   const REAL *p2;
   const REAL *p3;
   long nb;
   switch(mTwinningOption.GetChoice())
   {
      case 0:
      {
         p1=mCalcIntensity.data();
         p2=mObsIntensity.data();
         p3=mWeight.data();
         nb=mNbReflUsed;
         break;
      }
      case 1:
      {
         p1=mTwinnedIcalcSum.data();
         p2=mTwinnedIobsSum.data();
         p3=mTwinnedWeight.data();
         nb=mTwinnedIobsSum.numElements();
         break;
      }
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
   if(1==mTwinningOption.GetChoice()) mTwinnedIcalcSum*= tmp1/tmp2;
   mClockIcalc.Click();
}

void DiffractionDataSingleCrystal::FitScaleFactorForR()
{
   TAU_PROFILE("DiffractionData::FitScaleFactorForR()","void ()",TAU_DEFAULT);
   VFN_DEBUG_MESSAGE("DiffractionData::FitScaleFactorForR()",3);
   if(mHasObservedData==false)
   {//throw exception here ?
      throw ObjCrystException("DiffractionData::FitScaleFactorForR() Cannot compute R \
         or scale factor: there is no observed data !");
   }
   REAL tmp1=0;
   REAL tmp2=0;
   const REAL *p1;
   const REAL *p2;
   long nb;
   switch(mTwinningOption.GetChoice())
   {
      case 0:
      {
         p1=mCalcIntensity.data();
         p2=mObsIntensity.data();
         nb=mNbReflUsed;
         break;
      }
      case 1:
      {
         p1=mTwinnedIcalcSum.data();
         p2=mTwinnedIobsSum.data();
         nb=mTwinnedIobsSum.numElements();
         break;
      }
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
   if(1==mTwinningOption.GetChoice()) mTwinnedIcalcSum*= tmp1/tmp2;
   mClockIcalc.Click();
}

REAL DiffractionDataSingleCrystal::GetBestRFactor()
{
   TAU_PROFILE("DiffractionData::GetBestRFactor()","void ()",TAU_DEFAULT);
   VFN_DEBUG_MESSAGE("DiffractionData::GetBestRFactor()",3);
   if(mHasObservedData==false)
   {
      throw ObjCrystException("DiffractionData::GetBestRFactor() Cannot compute R \
         or scale factor: there is no observed data !");
   }
   this->FitScaleFactorForR();
   return this->GetR();
}

void DiffractionDataSingleCrystal::SetSigmaToSqrtIobs()
{
   for(long i=0;i<mObsIntensity.numElements();i++) mObsSigma(i)=sqrt(fabs(mObsIntensity(i)));
   if(0!=mTwinningOption.GetChoice()) mClockPrepareTwinningCorr.Reset();
}

void DiffractionDataSingleCrystal::SetWeightToInvSigma2(const REAL minRelatSigma)
{
   //:KLUDGE: If less than 1e-6*max, set to 0.... Do not give weight to unobserved points
   const REAL min=MaxAbs(mObsSigma)*minRelatSigma;
   for(long i=0;i<mObsSigma.numElements();i++)
   {
      if(mObsSigma(i)<min) mWeight(i)=0 ; else  mWeight(i) =1./mObsSigma(i)/mObsSigma(i);
   }
   if(0!=mTwinningOption.GetChoice()) mClockPrepareTwinningCorr.Reset();
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

unsigned int DiffractionDataSingleCrystal::GetNbCostFunction()const {return 2;}

const string& DiffractionDataSingleCrystal::GetCostFunctionName(const unsigned int id)const
{
   static string costFunctionName[2];
   if(0==costFunctionName[0].length())
   {
      costFunctionName[0]="Best R()";
      costFunctionName[1]="Best Rw()";
   }
   switch(id)
   {
      case 0: return costFunctionName[0];
      case 1: return costFunctionName[1];
      default:
      {
         cout << "DiffractionDataSingleCrystal::GetCostFunctionName(): Not Found !" <<endl;
         throw 0;
      }
   }
}

const string& DiffractionDataSingleCrystal::GetCostFunctionDescription(const unsigned int id)const
{
   static string costFunctionDescription[2];
   if(0==costFunctionDescription[0].length())
   {
      costFunctionDescription[0]="Crystallographic, unweighted R-factor with best scale";
      costFunctionDescription[1]="Crystallographic, weigthed R-factor with best scale";
   }
   switch(id)
   {
      case 0: return costFunctionDescription[0];
      case 1: return costFunctionDescription[1];
      default:
      {
         cout << "DiffractionDataSingleCrystal::GetCostFunctionDescription(): Not Found !" 
              <<endl;
         throw 0;
      }
   }
}

REAL DiffractionDataSingleCrystal::GetCostFunctionValue(const unsigned int n)
{
   VFN_DEBUG_MESSAGE("DiffractionDataSingleCrystal::GetCostFunctionValue():"<<mName,4)
   this->CalcIcalc();
   switch(n)
   {
      case 0: this->FitScaleFactorForR()  ;return this->GetR();
      case 1: this->FitScaleFactorForRw() ;return this->GetRw();
      default:
      {
         cout << "DiffractionDataSingleCrystal::GetCostFunctionValue(): Not Found !" <<endl;
         throw 0;
      }
   }
}

void DiffractionDataSingleCrystal::InitRefParList()
{
   VFN_DEBUG_MESSAGE("DiffractionDataSingleCrystal::InitRefParList()",5)
   //:TODO:
//   throw ObjCrystException("DiffractionDataSingleCrystal::InitRefParList() :
// not yet implemented !");
   //this->ResetParList();
   cout << "DiffractionDataSingleCrystal::InitRefParList():no parameters !" <<endl;
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
        && ((0==mTwinningOption.GetChoice()) || (mClockPrepareTwinningCorr<mClockIcalc)) ) return;
   
   mCalcIntensity=mFhklCalcSq;
   mCalcIntensity*=mScaleFactor;
   if(1==mTwinningOption.GetChoice())
   {
      this->PrepareTwinningCalc();
      const long nbGroup=mTwinnedGroupIndex.numElements();
      mTwinnedIcalcSum.resize(nbGroup);
      mTwinnedIcalcSum=0;
      long first=0;
      for(long i=0;i<nbGroup;i++)
      {
         for(long j=first;j<mTwinnedGroupIndex(i);j++)
         {
            #if 0 //def __DEBUG__
            cout <<endl << mIntH(j)<<" "<< mIntK(j)<<" "<< mIntL(j)<<" " << mObsIntensity(j);
            #endif
            mTwinnedIcalcSum(i)+=mCalcIntensity(j);
         }
         first=mTwinnedGroupIndex(i);
         #if 0 //def __DEBUG__
         cout  <<endl<< "         "<< mTwinnedIcalcSum(i) <<" "<< mTwinnedIobsSum(i)<<endl;
         #endif
      }
   }
   mClockIcalc.Click();
}

CrystVector_long DiffractionDataSingleCrystal::SortReflectionByTheta(const REAL maxTheta)
{
   TAU_PROFILE("DiffractionDataSingleCrystal::SortReflectionByTheta()","void ()",TAU_DEFAULT);
   VFN_DEBUG_ENTRY("DiffractionDataSingleCrystal::SortReflectionByTheta()",5)
   const CrystVector_long index=this->ScatteringData::SortReflectionByTheta(maxTheta);
   
   if(mObsIntensity.numElements()==mNbRefl)
   {
      CrystVector_REAL tmp;
      tmp=mObsIntensity;
      mObsIntensity.resize(mNbRefl);
      for(long i=0;i<mNbRefl;i++) mObsIntensity(i)=tmp(index(i));

      tmp=mObsSigma;
      mObsSigma.resize(mNbRefl);
      for(long i=0;i<mNbRefl;i++) mObsSigma(i)=tmp(index(i));

      tmp=mWeight;
      mWeight.resize(mNbRefl);
      for(long i=0;i<mNbRefl;i++) mWeight(i)=tmp(index(i));
   }
   else
   {// if there are no observed values, enter dummy ones
      mObsIntensity.resize(mNbRefl);
      mObsSigma.resize(mNbRefl);
      mWeight.resize(mNbRefl);
      mObsIntensity=100.;
      mObsSigma=1.;
      mWeight=1.;
   }
   VFN_DEBUG_EXIT("DiffractionDataSingleCrystal::SortReflectionByTheta()",5)
   return index;
}

void DiffractionDataSingleCrystal::InitOptions()
{
   static string TwinningOption;
   static string TwinningOptionChoices[2];
   static bool needInitNames=true;
   if(true==needInitNames)
   {
      TwinningOption="Twinning correction";
      TwinningOptionChoices[0]="None";
      TwinningOptionChoices[1]="Sum metrically equivalent reflections";
      needInitNames=false;
   }
   mTwinningOption.Init(2,&TwinningOption,TwinningOptionChoices);
   mTwinningOption.SetChoice(0);
   this->AddOption(&mTwinningOption);
}

void DiffractionDataSingleCrystal::PrepareTwinningCalc() const
{
   if(mClockPrepareTwinningCorr>mClockHKL) return;
   VFN_DEBUG_ENTRY("DiffractionDataSingleCrystal::PrepareTwinningCalc()",5)
   // first get the index of reflections which limit each block of summed reflections
   long nbGroup=0;
   {
      const REAL dSiThOvLa=.0001;
      mTwinnedGroupIndex.resize(mNbReflUsed);
      this->CalcSinThetaLambda();
      REAL sithol0=mSinThetaLambda(0)+dSiThOvLa;
      for(long i=1;i<mNbReflUsed;i++)
      {
         if(mSinThetaLambda(i)>sithol0)
         {
            mTwinnedGroupIndex(nbGroup++)=i;
            sithol0=mSinThetaLambda(i)+dSiThOvLa;
         }
      }
      mTwinnedGroupIndex(nbGroup++)=mNbReflUsed;
      mTwinnedGroupIndex.resizeAndPreserve(nbGroup);
   }
   // Calculate summed Iobs and weight
   {
      mTwinnedIobsSum.resize(nbGroup);
      mTwinnedIobsSum=0;
      mTwinnedWeight.resize(nbGroup);
      mTwinnedWeight=0;
      long first=0;
      for(long i=0;i<nbGroup;i++)
      {
         for(long j=first;j<mTwinnedGroupIndex(i);j++) 
         {
            mTwinnedIobsSum(i)+=mObsIntensity(j);
            mTwinnedWeight(i)+=mObsSigma(j)*mObsSigma(j);
         }
         mTwinnedWeight(i)=1./mTwinnedWeight(i);
         first=mTwinnedGroupIndex(i);
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
