/*
* ObjCryst++ : a Crystallographic computing library in C++
*
*  (c) 2000 Vincent FAVRE-NICOLIN
*           Laboratoire de Cristallographie
*           24, quai Ernest-Ansermet, CH-1211 Geneva 4, Switzerland
*  Contact: Vincent.Favre-Nicolin@cryst.unige.ch
*           Radovan.Cerny@cryst.unige.ch
*
*/
/*
*  source file for LibCryst++ PowderPattern class
*
*/

#include <cmath>

#include <typeinfo>
#include <stdio.h> //for sprintf()

#include "ObjCryst/PowderPattern.h"
#include "Quirks/VFNDebug.h"
#include "Quirks/VFNStreamFormat.h"
#ifdef __WX__CRYST__
   #include "wxCryst/wxPowderPattern.h"
#endif

#include <fstream>
#include <iomanip>

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
}

PowderPatternComponent::PowderPatternComponent(const PowderPatternComponent &old):
mIsScalable(old.mIsScalable),
mpParentPowderPattern(old.mpParentPowderPattern)
{}

PowderPatternComponent::~PowderPatternComponent()
{
   gPowderPatternComponentRegistry.DeRegister(*this);
}

const string PowderPatternComponent::GetClassName() const{return "PowderPatternComponent";}

bool PowderPatternComponent::IsScalable()const {return mIsScalable;}

const RefinableObjClock& PowderPatternComponent::GetClockPowderPatternCalc()const
{
   return mClockPowderPatternCalc;
}

////////////////////////////////////////////////////////////////////////
//
//        PowderPatternBackground
//
////////////////////////////////////////////////////////////////////////
PowderPatternBackground::PowderPatternBackground():
mBackgroundType(POWDER_BACKGROUND_LINEAR),mBackgroundNbPoint(0)
{}

PowderPatternBackground::PowderPatternBackground(const  PowderPatternBackground &old):
mBackgroundType(old.mBackgroundType),mBackgroundNbPoint(old.mBackgroundNbPoint),
mBackgroundInterpPoint2Theta(old.mBackgroundInterpPoint2Theta),
mBackgroundInterpPointIntensity(old.mBackgroundInterpPointIntensity)
{}

PowderPatternBackground::~PowderPatternBackground(){}
const string PowderPatternBackground::GetClassName() const{return "PowderPatternBackground";}

void PowderPatternBackground::SetParentPowderPattern(const PowderPattern &s)
{
   mpParentPowderPattern = &s;
}
const CrystVector_double& PowderPatternBackground::GetPowderPatternCalc()const
{
   this->CalcPowderPattern();
   return mPowderPatternCalc;
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
   CrystVector_double bckgd2Theta(max),bckgd(max);
   
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
   bckgd2Theta.resizeAndPreserve(nbPoints);
   bckgd.resizeAndPreserve(nbPoints);
   //cout << bckgd << endl;
   mBackgroundNbPoint=nbPoints;
   mBackgroundType=POWDER_BACKGROUND_LINEAR;
   mBackgroundInterpPoint2Theta=bckgd2Theta;
   mBackgroundInterpPoint2Theta*= DEG2RAD;
   mBackgroundInterpPointIntensity=bckgd;
   //Init refinable parameters
   {
      double *p=mBackgroundInterpPointIntensity.data();
      char buf [10];
      string str="Background_Point_";
      //for(int i=0;i<3;i++)
      for(int i=0;i<mBackgroundNbPoint;i++)
      {
         sprintf(buf,"%d",i);
         RefinablePar tmp(str+(string)buf,p++,
                           0.,1000.,gpRefParTypeScattDataBackground,REFPAR_DERIV_STEP_RELATIVE,
                           false,true,true,false,1.);
         tmp.AssignClock(mClockBackgroundPoint);
         tmp.SetDerivStep(1e-3);
         this->AddPar(tmp);
      }
   }
   
   mClockBackgroundPoint.Click();
   
   VFN_DEBUG_MESSAGE("PowderPatternBackground::ImportUserBackground():finished",5)
}

void PowderPatternBackground::SetUseFastLessPreciseFunc(const bool useItOrNot)
{//:TODO
}

void PowderPatternBackground::SetUseOnlyLowAngleData(const bool useOnlyLowAngle,
                                                      const double angle)
{
   if(  (mUseOnlyLowAngleData != useOnlyLowAngle)
      ||( fabs(mUseOnlyLowAngleDataLimit-angle)>.001))
         mClockPowderPatternCalc.Reset();
   mUseOnlyLowAngleData=useOnlyLowAngle;
   mUseOnlyLowAngleDataLimit=angle;
}

void PowderPatternBackground::CalcPowderPattern() const
{
   //:TODO: This needs serious optimization !
   if(   (mClockPowderPatternCalc>mClockBackgroundPoint)
       &&(mClockPowderPatternCalc>mpParentPowderPattern->GetClockPowderPatternPar())) return;
   TAU_PROFILE("PowderPatternBackground::CalcPowderPattern()","void ()",TAU_DEFAULT);
   VFN_DEBUG_MESSAGE("PowderPatternBackground::CalcPowderPattern()",3);
   
   switch(mBackgroundType)
   {
      case POWDER_BACKGROUND_LINEAR:
      {
         VFN_DEBUG_MESSAGE("PowderPatternBackground::CalcPowderPattern()..Linear",2)
         double t1,t2,b1,b2,t;
         const long nbPoint=mpParentPowderPattern->GetNbPoint();
         mPowderPatternCalc.resize(nbPoint);
         //mPowderPatternCalc=0.;
         double *b=mPowderPatternCalc.data();
         t1=mBackgroundInterpPoint2Theta(0);
         t2=mBackgroundInterpPoint2Theta(1);
         b1=mBackgroundInterpPointIntensity(0);
         b2=mBackgroundInterpPointIntensity(1);
         t=mpParentPowderPattern->Get2ThetaMin();
         long point=1;
         for(long i=0;i<nbPoint;i++)
         {
            if(t >= t2)
            {
               if(point < mBackgroundNbPoint-1)
               {
                  b1=b2;
                  t1=t2;
                  b2=mBackgroundInterpPointIntensity(point+1);
                  t2=mBackgroundInterpPoint2Theta(point+1);
                  point++ ;
               }
            }
            *b = (b1*(t2-t)+b2*(t-t1))/(t2-t1) ;
            b++;
            t+=mpParentPowderPattern->Get2ThetaStep();
         }
         break;
      }
      case POWDER_BACKGROUND_CUBIC_SPLINE:
      {
         //:TODO:
         break;
      }
   }
   mClockPowderPatternCalc.Click();
   VFN_DEBUG_MESSAGE("PowderPatternBackground::CalcPowderPattern():End",3);
}

void PowderPatternBackground::SetRadiation(const Radiation& rad)
{
}
void PowderPatternBackground::Prepare()
{
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
mFullProfileWidthFactor(5.),
mCagliotiU(0),mCagliotiV(0),mCagliotiW(3e-5),
mPseudoVoigtEta0(0.5),mPseudoVoigtEta1(0.),mUseAsymmetricProfile(false),
mNeedLorentzCorr(true),
mNeedPolarCorr(true),
mNeedSlitApertureCorr(true),
mPolarAfactor(1.),
mUseOnlyLowAngleData(false),
mUseOnlyLowAngleDataLimit(0.)
{
   mIsScalable=true;
   this->InitOptions();
   mReflectionProfileType.SetChoice(PROFILE_PSEUDO_VOIGT);
}

PowderPatternDiffraction::PowderPatternDiffraction(const PowderPatternDiffraction &old):
mReflectionProfileType(old.mReflectionProfileType),
mFullProfileWidthFactor(old.mFullProfileWidthFactor),
mCagliotiU(old.mCagliotiU),mCagliotiV(old.mCagliotiV),mCagliotiW(old.mCagliotiW),
mPseudoVoigtEta0(old.mPseudoVoigtEta0),mPseudoVoigtEta1(old.mPseudoVoigtEta1),
mUseAsymmetricProfile(old.mUseAsymmetricProfile),
mNeedLorentzCorr(old.mNeedLorentzCorr),
mNeedPolarCorr(old.mNeedPolarCorr),
mNeedSlitApertureCorr(old.mNeedSlitApertureCorr),
mPolarAfactor(old.mPolarAfactor),
mUseOnlyLowAngleData(old.mUseOnlyLowAngleData),
mUseOnlyLowAngleDataLimit(old.mUseOnlyLowAngleDataLimit),
mUseOnlyLowAngleData_SavedH(old.mUseOnlyLowAngleData_SavedH),
mUseOnlyLowAngleData_SavedK(old.mUseOnlyLowAngleData_SavedK),
mUseOnlyLowAngleData_SavedL(old.mUseOnlyLowAngleData_SavedL)
{}

PowderPatternDiffraction::~PowderPatternDiffraction()
{}
const string PowderPatternDiffraction::GetClassName() const
{return "PowderPatternDiffraction";}

PowderPatternDiffraction* PowderPatternDiffraction::CreateCopy()const
{
   return new PowderPatternDiffraction(*this);
}

void PowderPatternDiffraction::SetParentPowderPattern(const PowderPattern &s)
{
    mpParentPowderPattern = &s;
}

const CrystVector_double& PowderPatternDiffraction::GetPowderPatternCalc()const
{
   this->CalcPowderPattern();
   return mPowderPatternCalc;
}

void PowderPatternDiffraction::SetReflectionProfilePar(const ReflectionProfileType prof,
                          const double fwhmCagliotiW,
                          const double fwhmCagliotiU,
                          const double fwhmCagliotiV,
                          const double eta0,
                          const double eta1)
{
   VFN_DEBUG_MESSAGE("PowderPatternDiffraction::SetReflectionProfilePar()",5)
   mReflectionProfileType.SetChoice(prof);
   mCagliotiU=fwhmCagliotiU;
   mCagliotiV=fwhmCagliotiV;
   mCagliotiW=fwhmCagliotiW;
   mPseudoVoigtEta0=eta0;
   mPseudoVoigtEta1=eta1;
   mClockProfilePar.Click();
}

void PowderPatternDiffraction::SetUseFastLessPreciseFunc(const bool useItOrNot)
{
   this->ScatteringData::SetUseFastLessPreciseFunc(useItOrNot);
}

void PowderPatternDiffraction::SetUseOnlyLowAngleData(const bool useOnlyLowAngle,
                                                       const double angle)
{
   if((mUseOnlyLowAngleData==useOnlyLowAngle) &&
      (fabs(mUseOnlyLowAngleDataLimit-angle) < 1e-4)) return;
      
   VFN_DEBUG_MESSAGE("PowderPatternDiffraction::SetUseOnlyLowAngleData(theta)",5)
   
   bool useSavedArrays=false;
   if((true==useOnlyLowAngle)&&(true==mUseOnlyLowAngleData))
      useSavedArrays=true;
   
   mUseOnlyLowAngleDataLimit=angle;
   mUseOnlyLowAngleData=useOnlyLowAngle;
   
   if((true==mUseOnlyLowAngleData) && (mUseOnlyLowAngleDataLimit <= 0.))
   {//shit happens
      mUseOnlyLowAngleData=false;
      throw ObjCrystException("PowderPatternDiffraction::SetUseOnlyLowAngleData(bool,angle) \
Setting to use only low angle data, but given limit is <= 0 !");
   }
   
   if((true==mUseOnlyLowAngleData) &&  ((2*mUseOnlyLowAngleDataLimit) >=
       mpParentPowderPattern->Get2ThetaMax()))
   {
      mUseOnlyLowAngleData=false;
      cout << "PowderPatternDiffraction::SetUseOnlyLowAngleData(theta)"<<endl;
      cout << "-> The given theta low limit is higher than the spectrum limit !"
           << " Stupid you... Ignoring limit" <<endl;
      return;
   }
   
   if(true==mUseOnlyLowAngleData)
   {
      if(true==useSavedArrays)
      {
         this ->SetHKL( mUseOnlyLowAngleData_SavedH,
                        mUseOnlyLowAngleData_SavedK,
                        mUseOnlyLowAngleData_SavedL);
      }
      else
      {
      //Make copies of arrays for when the flag is un-raised
         mUseOnlyLowAngleData_SavedH=mIntH;
         mUseOnlyLowAngleData_SavedK=mIntK;
         mUseOnlyLowAngleData_SavedL=mIntL;
      }
      this->SortReflectionByTheta(mUseOnlyLowAngleDataLimit);
   }
   else
   {//Restore full data
      this ->SetHKL( mUseOnlyLowAngleData_SavedH,
                     mUseOnlyLowAngleData_SavedK,
                     mUseOnlyLowAngleData_SavedL);
      mUseOnlyLowAngleData_SavedH.resize(0);
      mUseOnlyLowAngleData_SavedK.resize(0);
      mUseOnlyLowAngleData_SavedL.resize(0);
   }
}

void PowderPatternDiffraction::GenHKLFullSpace()
{
   VFN_DEBUG_MESSAGE("PowderPatternDiffraction::GenHKLFullSpace():",3)
   this->ScatteringData::GenHKLFullSpace(mpParentPowderPattern->Get2ThetaMax()/2,true);
}

void PowderPatternDiffraction::CalcPowderPattern() const
{
   VFN_DEBUG_MESSAGE("PowderPatternDiffraction::CalcPowderPattern():",3)
   TAU_PROFILE("PowderPatternDiffraction::CalcPowderPattern()","void (bool)",TAU_DEFAULT);
   TAU_PROFILE_TIMER(timer1,"PowderPatternDiffraction::CalcPowderPattern1(Calc  F(hkl))"\
                     ,"", TAU_FIELD);
   TAU_PROFILE_TIMER(timer2,"PowderPatternDiffraction::CalcPowderPattern2(Calc Profiles)"\
                     ,"", TAU_FIELD);
   TAU_PROFILE_TIMER(timer3,"PowderPatternDiffraction::CalcPowderPattern3(Apply Profiles)"\
                     ,"", TAU_FIELD);

   //:TODO: Synchronize some other way the two wavelengths... (or remove one ?)
   //can't do this since it's not const...
   //
   //if(mpParentPowderPattern->GetRadiation().GetClockWavelength()
   //   >mRadiation.GetClockWavelength())
   //      mRadiation.SetWavelength(mpParentPowderPattern->GetRadiation().GetWavelength()(0));

   TAU_PROFILE_START(timer1);
   this->CalcIhkl();
   TAU_PROFILE_STOP(timer1);
   TAU_PROFILE_START(timer2);
   this->CalcPowderReflProfile();
   TAU_PROFILE_STOP(timer2);
   
   if(  (mClockPowderPatternCalc>mClockIhklCalc)
      &&(mClockPowderPatternCalc>mClockProfileCalc)) return;
   
   //mpCrystal->Print();
   if(true) //:TODO: false == mUseFastLessPreciseFunc
   {
      TAU_PROFILE_START(timer3);
      double theta;
      long thetaPt,first,last,shift;
      const long nbPoints=2*mSavedPowderReflProfileNbPoint+1;
      long step;
      const long nbRefl=this->GetNbRefl();
      VFN_DEBUG_MESSAGE("PowderPatternDiffraction::CalcPowderPattern\
Applying profiles for "<<nbRefl<<" reflections",3)

      const long  specNbPoints=mpParentPowderPattern->GetNbPoint();
      mPowderPatternCalc.resize(specNbPoints);
      mPowderPatternCalc=0;
      
      for(long i=0;i<nbRefl;i += step)
      {
         double intensity=0.;
         //check if the next reflection is at the same theta. If this is true,
         //Then assume that the profile is exactly the same.
         for(step=0; ;)
         {
            intensity += mIhklCalc(i + step++);
            if( (i+step) >= nbRefl) break;
            if(mTheta(i+step) > (mTheta(i)+1e-4) ) break;
         }
         //intensity *= 1e-14;
         switch(mRadiation.GetWavelengthType())
         {
            case WAVELENGTH_MONOCHROMATIC:
            {
               theta =  mpParentPowderPattern->Get2ThetaCorr(2*mTheta(i));
               thetaPt= mpParentPowderPattern->Get2ThetaCorrPixel(2*mTheta(i));
               VFN_DEBUG_MESSAGE("Apply profile(Monochromatic)Refl("<<i<<")"\
                  <<mIntH(i)<<" "<<mIntK(i)<<" "<<mIntL(i)<<" "\
                  <<"  I="<<intensity<<"  2Theta="<<2*mTheta(i)*RAD2DEG\
                  <<",pixel #"<<thetaPt,2)
               
               first=thetaPt-mSavedPowderReflProfileNbPoint;
               
               if( first >= specNbPoints) continue;
               if( first < 0)
               {
                  shift = -first;
                  first =0;
               } else shift =0;
               last=thetaPt+mSavedPowderReflProfileNbPoint;
               if( last > specNbPoints) last=specNbPoints;
               VFN_DEBUG_MESSAGE("first:"<<first<<"  last:"<<last,1)
               {
                  const double *p2 = mSavedPowderReflProfile.data() + i*nbPoints +shift;
                  double *p3 = mPowderPatternCalc.data()+first;
                  for(long j=first;j<last;j++) *p3++ += *p2++ * intensity;
               }
               break;
            }
            case WAVELENGTH_ALPHA12:
            {
               VFN_DEBUG_MESSAGE("Apply profile(X-Ray Tube)Refl="<<i<<\
               "  I="<<intensity<<"  Theta="<<mTheta(i)*RAD2DEG,2)
               //:TODO: Use only ONE profile array for both alpha1&2 (faster)
               {//Alpha1
                  intensity /= (1+mRadiation.GetXRayTubeAlpha2Alpha1Ratio());
                  theta=mTheta(i);
                  theta+=mTanTheta(i)*(
                     -mRadiation.GetXRayTubeDeltaLambda()
                      *mRadiation.GetXRayTubeAlpha2Alpha1Ratio())
                        /(1+mRadiation.GetXRayTubeAlpha2Alpha1Ratio())
                           /mRadiation.GetWavelength()(0);
                  thetaPt= mpParentPowderPattern->Get2ThetaCorrPixel(2*theta);
                  theta=mpParentPowderPattern->Get2ThetaCorr(2*theta);
                  first=thetaPt-mSavedPowderReflProfileNbPoint;
                  if( first >= specNbPoints) continue;
                  if( first < 0)
                  {
                     shift = -first;
                     first =0;
                  } else shift =0;
                  last=thetaPt+mSavedPowderReflProfileNbPoint;
                  if( last > specNbPoints) last=specNbPoints;
                  VFN_DEBUG_MESSAGE("first:"<<first<<"  last:"<<last,1)
                  {
                     const double *p2 = mSavedPowderReflProfile.data() + i*nbPoints +shift;
                     double *p3 = mPowderPatternCalc.data()+first;
                     for(long j=first;j<last;j++) *p3++ += *p2++ * intensity;
                  }
               }
               {//Alpha2
                  intensity *= mRadiation.GetXRayTubeAlpha2Alpha1Ratio();
                  theta=mTheta(i);
                  theta+=mTanTheta(i)*(mRadiation.GetXRayTubeDeltaLambda()
                              /(1+mRadiation.GetXRayTubeAlpha2Alpha1Ratio()))
                           /mRadiation.GetWavelength()(0);
                  thetaPt= mpParentPowderPattern->Get2ThetaCorrPixel(2*theta);
                  theta=mpParentPowderPattern->Get2ThetaCorr(2*theta);
                  first=thetaPt-mSavedPowderReflProfileNbPoint;
                  if( first >= specNbPoints) continue;
                  if( first < 0)
                  {
                     shift = -first;
                     first =0;
                  } else shift =0;
                  last=thetaPt+mSavedPowderReflProfileNbPoint;
                  if( last > specNbPoints) last=specNbPoints;
                  VFN_DEBUG_MESSAGE("first:"<<first<<"  last:"<<last,0)
                  {
                     const double *p2 = mSavedPowderReflProfile.data() + i*nbPoints +shift;
                     double *p3 = mPowderPatternCalc.data()+first;
                     for(long j=first;j<last;j++) *p3++ += *p2++ * intensity;
                  }
               }
               break;
            }
            default: throw ObjCrystException("PowderPatternDiffraction::CalcPowderPattern():\
Beam must either be monochromatic or from an XRay Tube !!");
         }
      }
      TAU_PROFILE_STOP(timer3);
   }
   else
   {
      //:TODO:
      throw ObjCrystException("PowderPatternDiffraction::CalcPowderPattern() : \
         FAST option not yet implemented !");
   }
   mClockPowderPatternCalc.Click();
   VFN_DEBUG_MESSAGE("PowderPatternDiffraction::CalcPowderPattern: End.",3)
}

void PowderPatternDiffraction::CalcPowderReflProfile()const
{
   this->CalcSinThetaLambda();
   //mClockProfileCalc.Print();
   //mRadiation.GetClockWavelength().Print();
   //mRadiation.GetClockRadiation().Print();
   if(  (mClockProfileCalc>mClockProfilePar)
      &&(mClockProfileCalc>mReflectionProfileType.GetClock())
      &&(mClockProfileCalc>mClockTheta)
      &&(mClockProfileCalc>mRadiation.GetClockWavelength())
      &&(mClockProfileCalc>mClockHKL)) return;
   
   TAU_PROFILE("PowderPatternDiffraction::CalcPowderReflProfile()","void (bool)",TAU_DEFAULT);
   VFN_DEBUG_MESSAGE("PowderPatternDiffraction::CalcPowderReflProfile()",5)

   const double specMin     =mpParentPowderPattern->Get2ThetaMin();
   //const double specMax     =mpParentPowderPattern->Get2ThetaMax();
   const double specStep    =mpParentPowderPattern->Get2ThetaStep();
   const long  specNbPoints=mpParentPowderPattern->GetNbPoint();
   
   long thetaPt;
   double fwhm,tmp;
   double *p1,*p2;
   CrystVector_double ttheta,tmp2theta,reflProfile,tmpV;
   //Width of calc profiles
   VFN_DEBUG_MESSAGE("PowderPatternDiffraction::CalcPowderReflProfile():\
Computing Widths",5)
   {
      tmp = specMin+specStep*specNbPoints;
      tmp/= 2;
      fwhm=      mCagliotiW + mCagliotiV*tmp + mCagliotiU*tmp*tmp;
      if(fwhm<1e-10) fwhm=1e-10;
      fwhm=sqrt(fwhm);
                
      mSavedPowderReflProfileNbPoint =(long)(mFullProfileWidthFactor
                                             *fwhm/specStep);
      VFN_DEBUG_MESSAGE("PowderPatternDiffraction::CalcPowderReflProfile():"<<\
      "Profiles half-width="<<mFullProfileWidthFactor*fwhm*RAD2DEG<<" ("<<\
      mSavedPowderReflProfileNbPoint<<" points)",5)
   }
   //Prepare arrays
      ttheta.resize(2*mSavedPowderReflProfileNbPoint+1);
      reflProfile.resize(2*mSavedPowderReflProfileNbPoint+1);
      p2=ttheta.data();
      for(double i=-mSavedPowderReflProfileNbPoint;i<=mSavedPowderReflProfileNbPoint;i++)
         *p2++ = i * specStep;
      mSavedPowderReflProfile.resize(this->GetNbRefl(),2*mSavedPowderReflProfileNbPoint+1);
   //Calc all profiles
   p1=mSavedPowderReflProfile.data();
   VFN_DEBUG_MESSAGE("PowderPatternDiffraction::CalcPowderReflProfile():\
Computing all Profiles",5)
   for(long i=0;i<this->GetNbRefl();i++)
   {
   VFN_DEBUG_MESSAGE("PowderPatternDiffraction::CalcPowderReflProfile():\
Computing all Profiles: Reflection #"<<i,2)
      tmp=mTheta(i);
		fwhm=mCagliotiW + mCagliotiV*tmp + mCagliotiU*tmp*tmp;
		if(fwhm<1e-10) fwhm=1e-10;
      fwhm =sqrt(fwhm);
                 
      double powderAsym=1.;
      //if(true == mUseAsymmetricProfile) 
      //   powderAsym=mPowderAsymA0+mPowderAsymA1/sin(tmp)+mPowderAsymA2/sin(tmp)/sin(tmp);
         
      //:KLUDGE: Need to use the same 'thetaPt' when calculating the complete powder spectrum
      thetaPt =(long) ((2*tmp-specMin)/specStep);
      tmp2theta = ttheta;
      tmp2theta += 2*tmp-(specMin + thetaPt * specStep);
      switch(mReflectionProfileType.GetChoice())
      {
         case PROFILE_GAUSSIAN:
         {
            if(true == mUseAsymmetricProfile)
               reflProfile=PowderProfileGauss(tmp2theta,fwhm,powderAsym);
            else reflProfile=PowderProfileGauss(tmp2theta,fwhm);
            break;
         }
         case PROFILE_LORENTZIAN:
         {
            if(true == mUseAsymmetricProfile)
               reflProfile=PowderProfileLorentz(tmp2theta,fwhm,powderAsym);
            else reflProfile=PowderProfileLorentz(tmp2theta,fwhm);
            break;
         }
         case PROFILE_PSEUDO_VOIGT:
         {
            if(true == mUseAsymmetricProfile)
               reflProfile=PowderProfileGauss(tmp2theta,fwhm,powderAsym);
            else reflProfile=PowderProfileGauss(tmp2theta,fwhm);
            reflProfile *= 1-(mPseudoVoigtEta0+2*mTheta(i)*mPseudoVoigtEta1);
            if(true == mUseAsymmetricProfile)
               tmpV=PowderProfileLorentz(tmp2theta,fwhm,powderAsym);
            else tmpV=PowderProfileLorentz(tmp2theta,fwhm);
            tmpV *= mPseudoVoigtEta0+2*mTheta(i)*mPseudoVoigtEta1;
            reflProfile += tmpV;
            break;
         }
         case PROFILE_PSEUDO_VOIGT_FINGER_COX_JEPHCOAT:
         {
            throw ObjCrystException(
               "PROFILE_PSEUDO_VOIGT_FINGER_COX_JEPHCOAT Not implemented");
            /*
            tmp2theta*=180./M_PI;//Keep degrees:external library used
            fwhm *=180./M_PI;
            double eta=mPowderPseudoVoigtEta0+2*mTheta(i)*mPowderPseudoVoigtEta1;
            double dPRdT, dPRdG, dPRdE , dPRdS, dPRdD;
            bool useAsym=true;
            if( (mPowderAsymSourceWidth < 1e-8) && (mPowderAsymSourceWidth < 1e-8) )
               useAsym=false;
               
            for(long j=0;j<2*mSavedPowderReflProfileNbPoint+1;j++)
               reflProfile(j)=Profval( eta,fwhm,mPowderAsymSourceWidth,
                                       mPowderAsymDetectorWidth,tmp2theta(j) ,
                                       0.,&dPRdT,&dPRdG,&dPRdE,&dPRdS,&dPRdD,
                                       useAsym);
            break;
            */
         }
         case PROFILE_PEARSON_VII: throw ObjCrystException("PEARSON_VII Not implemented yet");

      }
      p2=reflProfile.data();
      for(int j=0;j<2*mSavedPowderReflProfileNbPoint+1;j++) *p1++ = *p2++;
   }
   mClockProfileCalc.Click();
   VFN_DEBUG_MESSAGE("PowderPatternDiffraction::CalcPowderReflProfile():finished",5)
}

void PowderPatternDiffraction::CalcLorentzPolarCorr()const
{
   this->CalcSinThetaLambda();
   
   if(  (mClockLorentzPolarSlitCorrCalc>mClockLorentzPolarSlitCorrPar)
      &&(mClockLorentzPolarSlitCorrCalc>mClockTheta) ) return;
      
   TAU_PROFILE("PowderPatternDiffraction::CalcLorentzPolarCorr()","void ()",TAU_DEFAULT);
   VFN_DEBUG_MESSAGE("PowderPatternDiffraction::CalcLorentzPolarCorr()",3)
   mLorentzPolarSlitCorr.resize(this->GetNbRefl());
   mLorentzPolarSlitCorr=1.;
   if(true==mNeedLorentzCorr)
   {
      VFN_DEBUG_MESSAGE("->Need Lorentz Corr",2)
      for(long i=0;i<this->GetNbRefl();i++)
         mLorentzPolarSlitCorr(i) /=sin(2*mTheta(i));
   }
   if(true==mNeedPolarCorr)
   {
      VFN_DEBUG_MESSAGE("->Need Polarization Corr: A="<<mPolarAfactor,2)
      for(long i=0;i<this->GetNbRefl();i++)
         mLorentzPolarSlitCorr(i) *=(1.+mPolarAfactor*cos(mTheta(i))*cos(mTheta(i)))
                                 /(1.+mPolarAfactor);
   }
   if(true==mNeedSlitApertureCorr)
   {
      VFN_DEBUG_MESSAGE("->Need Slit Aperture Corr",2)
      for(long i=0;i<this->GetNbRefl();i++)
         mLorentzPolarSlitCorr(i) /= sin(mTheta(i));
   }
         
   mClockLorentzPolarSlitCorrCalc.Click();
   VFN_DEBUG_MESSAGE("PowderPatternDiffraction::CalcLorentzPolarCorr():finished",0)
}

void PowderPatternDiffraction::CalcIhkl() const
{
   this->CalcStructFactor();
   this->CalcLorentzPolarCorr();
   if(  (mClockIhklCalc>mClockLorentzPolarSlitCorrCalc)
      &&(mClockIhklCalc>mClockStructFactor)) return;
      
   VFN_DEBUG_MESSAGE("PowderPatternDiffraction::CalcIhkl()",3)
   TAU_PROFILE("PowderPatternDiffraction::CalcIhkl()","void ()",TAU_DEFAULT);
   const double *pr,*pi,*pcorr;
   const int *mult;
   double *p;
   
   pr=mFhklCalcReal.data();
   pi=mFhklCalcImag.data();
   pcorr=mLorentzPolarSlitCorr.data();
   
   mult=mMultiplicity.data();
   mIhklCalc.resize(mNbRefl);
   p=mIhklCalc.data();
   
   for(long i=0;i<mNbRefl;i++)
   {
      *p++ = *mult++ * (*pr * *pr + *pi * *pi) * *pcorr++;
      pr++;
      pi++;
   }
   //cout <<FormatVertVector<double>(mTheta,mIhklCalc,mMultiplicity,
   //                               mFhklCalcReal,mFhklCalcImag,mLorentzPolarSlitCorr);
   mClockIhklCalc.Click();
   VFN_DEBUG_MESSAGE("PowderPatternDiffraction::CalcIhkl():End",3)
}

void PowderPatternDiffraction::SetRadiation(const Radiation& rad)
{
   mRadiation=rad;
}

void PowderPatternDiffraction::Prepare()
{
   if(0==this->GetNbRefl()) this->GenHKLFullSpace();
}
void PowderPatternDiffraction::InitOptions()
{
   VFN_DEBUG_MESSAGE("PowderPatternDiffraction::InitOptions()",5)
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
   
   //Init parameters. This should not be done here !!!
   {
      RefinablePar tmp("U",&mCagliotiU,-1/RAD2DEG/RAD2DEG,1./RAD2DEG/RAD2DEG,
                        gpRefParTypeScattDataProfileWidth,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,RAD2DEG*RAD2DEG);
      tmp.AssignClock(mClockProfilePar);
      tmp.SetDerivStep(1e-5);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("V",&mCagliotiV,-1/RAD2DEG/RAD2DEG,1./RAD2DEG/RAD2DEG,
                        gpRefParTypeScattDataProfileWidth,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,RAD2DEG*RAD2DEG);
      tmp.AssignClock(mClockProfilePar);
      tmp.SetDerivStep(1e-5);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("W",&mCagliotiW,0,1./RAD2DEG/RAD2DEG,
                        gpRefParTypeScattDataProfileWidth,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,RAD2DEG*RAD2DEG);
      tmp.AssignClock(mClockProfilePar);
      tmp.SetDerivStep(1e-5);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("Eta0",&mPseudoVoigtEta0,0,1.,gpRefParTypeScattDataProfileType,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false);
      tmp.AssignClock(mClockProfilePar);
      tmp.SetDerivStep(1e-4);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("Eta1",&mPseudoVoigtEta1,-1,1.,gpRefParTypeScattDataProfileType,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false);
      tmp.AssignClock(mClockProfilePar);
      tmp.SetDerivStep(1e-4);
      this->AddPar(tmp);
   }
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
m2ThetaMin(0),m2ThetaStep(0),mNbPoint(0),mWavelength(1.),
m2ThetaZero(0.),m2ThetaDisplacement(0.),m2ThetaTransparency(0.),
mScaleFactor(20),mUseFastLessPreciseFunc(false),
mStatisticsExcludeBackground(false),mStatisticsUseIntegratedPeak(false),
mUseOnlyLowAngleData(false),mUseOnlyLowAngleDataLimit(0)
{
   mScaleFactor=1;
   mSubObjRegistry.SetName("SubObjRegistry for a PowderPattern object");
   mPowderPatternComponentRegistry.SetName("Powder Pattern Components");
   gPowderPatternRegistry.Register(*this);
   gTopRefinableObjRegistry.Register(*this);
   this->Init();
}

PowderPattern::PowderPattern(const PowderPattern &old):
m2ThetaMin(old.m2ThetaMin),m2ThetaStep(old.m2ThetaStep),mNbPoint(old.mNbPoint),
mWavelength(old.mWavelength),mRadiation(old.mRadiation),
m2ThetaZero(old.m2ThetaZero),m2ThetaDisplacement(old.m2ThetaDisplacement),
m2ThetaTransparency(old.m2ThetaTransparency),
mPowderPatternComponentRegistry(old.mPowderPatternComponentRegistry),
mScaleFactor(old.mScaleFactor),
mUseFastLessPreciseFunc(old.mUseFastLessPreciseFunc),
mStatisticsExcludeBackground(old.mStatisticsExcludeBackground),
mStatisticsUseIntegratedPeak(old.mStatisticsUseIntegratedPeak),
mUseOnlyLowAngleData(old.mUseOnlyLowAngleData),
mUseOnlyLowAngleDataLimit(old.mUseOnlyLowAngleDataLimit)
{
   mSubObjRegistry.SetName("SubObjRegistry for a PowderPattern :"+mName);
   gPowderPatternRegistry.Register(*this);
   gTopRefinableObjRegistry.Register(*this);
   this->Init();
}

PowderPattern::~PowderPattern()
{
   for(int i=0;i<mPowderPatternComponentRegistry.GetNb();i++)
   {
      mPowderPatternComponentRegistry.GetObj(i).DeRegisterClient(*this);
      mSubObjRegistry.DeRegister(mPowderPatternComponentRegistry.GetObj(i));
      delete &(mPowderPatternComponentRegistry.GetObj(i));
   }
   gPowderPatternRegistry.DeRegister(*this);
   gTopRefinableObjRegistry.DeRegister(*this);
}
const string PowderPattern::GetClassName() const {return "PowderPattern";}

void PowderPattern::AddPowderPatternComponent(PowderPatternComponent &comp)
{
   VFN_DEBUG_ENTRY("PowderPattern::AddPowderPatternComponent():"<<comp.GetName(),5)
   comp.SetParentPowderPattern(*this);
   mSubObjRegistry.Register(comp);
   comp.RegisterClient(*this);
   mClockPowderPatternCalc.Reset();
   comp.SetRadiation(mRadiation);
   comp.SetUseFastLessPreciseFunc(mUseFastLessPreciseFunc);
   mPowderPatternComponentRegistry.Register(comp);
	//:TODO: check if there are enough scale factors
   //mScaleFactor.resizeAndPreserve(mPowderPatternComponentRegistry.GetNb());
   mScaleFactor(mPowderPatternComponentRegistry.GetNb()-1)=1.;
	mClockScaleFactor.Click();
	if(comp.IsScalable())
   {//Init refinable parameter
      RefinablePar tmp("Scale_"+comp.GetName(),mScaleFactor.data()+mPowderPatternComponentRegistry.GetNb()-1,
                        0.,0.,gpRefParTypeScattDataScale,REFPAR_DERIV_STEP_RELATIVE,
                        false,true,true,false,1.);
      tmp.AssignClock(mClockScaleFactor);
      tmp.SetDerivStep(1e-4);
      this->AddPar(tmp);
   }
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

void PowderPattern::SetPowderPatternPar(const double tthetaMin,
                                          const double tthetaStep,
                                          unsigned long nbPoint)
{
   m2ThetaMin=tthetaMin;
   m2ThetaStep=tthetaStep;
   mNbPoint=nbPoint;
   mClockPowderPatternPar.Click();
}

unsigned long PowderPattern::GetNbPoint()const {return mNbPoint;}

void PowderPattern::SetRadiation(const Radiation &radiation)
{
   mRadiation=radiation;
   for(int i=0;i<mPowderPatternComponentRegistry.GetNb();i++)
      mPowderPatternComponentRegistry.GetObj(i).SetRadiation(mRadiation);
   mClockPowderPatternRadiation.Click();
}
const Radiation& PowderPattern::GetRadiation()const {return mRadiation;}

void PowderPattern::SetRadiationType(const RadiationType rad)
{
   mRadiation.SetRadiationType(rad);
   for(int i=0;i<mPowderPatternComponentRegistry.GetNb();i++)
      mPowderPatternComponentRegistry.GetObj(i).SetRadiation(mRadiation);
}

RadiationType PowderPattern::GetRadiationType()const {return mRadiation.GetRadiationType();}
void PowderPattern::SetWavelength(const double lambda)
{
   VFN_DEBUG_MESSAGE("PowderPattern::SetWavelength(lambda)",3)
   mRadiation.SetWavelength(lambda);
   for(int i=0;i<mPowderPatternComponentRegistry.GetNb();i++)
      mPowderPatternComponentRegistry.GetObj(i).SetRadiation(mRadiation);
}

void PowderPattern::SetWavelength(const string &XRayTubeElementName,const double alpha12ratio)
{
   VFN_DEBUG_MESSAGE("PowderPattern::SetWavelength(wavelength)",3)
   mRadiation.SetWavelength(XRayTubeElementName,alpha12ratio);
   for(int i=0;i<mPowderPatternComponentRegistry.GetNb();i++)
      mPowderPatternComponentRegistry.GetObj(i).SetRadiation(mRadiation);
}

double PowderPattern::GetWavelength()const{return mRadiation.GetWavelength()(0);}

void PowderPattern::SetUseFastLessPreciseFunc(const bool useItOrNot)
{
   mUseFastLessPreciseFunc=useItOrNot;
   for(int i=0;i<mPowderPatternComponentRegistry.GetNb();i++)
      mPowderPatternComponentRegistry.GetObj(i)
         .SetUseFastLessPreciseFunc(mUseFastLessPreciseFunc);
}

void PowderPattern::SetUseOnlyLowAngleData(const bool useOnlyLowAngle,const double angle)
{
   VFN_DEBUG_MESSAGE("PowderPattern::SetUseOnlyLowAngleData()",3)
   mUseOnlyLowAngleData=useOnlyLowAngle;
   mUseOnlyLowAngleDataLimit=angle;
   this->Prepare();//Will prepare components, if necessary
   for(int i=0;i<mPowderPatternComponentRegistry.GetNb();i++)
      mPowderPatternComponentRegistry.GetObj(i)
         .SetUseOnlyLowAngleData(useOnlyLowAngle,angle);
}

const CrystVector_double& PowderPattern::GetPowderPatternCalc()const
{
   this->CalcPowderPattern();
   return mPowderPatternCalc;
}

const CrystVector_double& PowderPattern::GetPowderPatternObs()const
{
   return mPowderPatternObs;
}

const CrystVector_double& PowderPattern::GetPowderPatternObsSigma()const
{
   return mPowderPatternObsSigma;
}

const CrystVector_double& PowderPattern::GetPowderPatternWeight()const
{
   return mPowderPatternWeight;
}

//CrystVector_double PowderPattern::Get2Theta()const

double PowderPattern::Get2ThetaMin()const {return m2ThetaMin;}

double PowderPattern::Get2ThetaStep()const {return m2ThetaStep;}

double PowderPattern::Get2ThetaMax()const 
{
   VFN_DEBUG_MESSAGE("PowderPattern::Get2ThetaMax()",2)
   return m2ThetaMin+m2ThetaStep*mNbPoint;
}

const RefinableObjClock& PowderPattern::GetClockPowderPatternCalc()const
{  return mClockPowderPatternCalc;}

const RefinableObjClock& PowderPattern::GetClockPowderPatternPar()const
{  return mClockPowderPatternPar;}

const RefinableObjClock& PowderPattern::GetClockPowderPatternRadiation()const
{  return mClockPowderPatternRadiation;}

void PowderPattern::Set2ThetaZero(const double newZero)
{
   m2ThetaZero=newZero;
   mClockPowderPatternPar.Click();
}

void PowderPattern::Set2ThetaDisplacement(const double displacement)
{
   m2ThetaDisplacement=displacement;
   mClockPowderPatternPar.Click();
}

void PowderPattern::Set2ThetaTransparency(const double transparency)
{
   m2ThetaTransparency=transparency;
   mClockPowderPatternPar.Click();
}

double PowderPattern::Get2ThetaCorr(const double ttheta)const
{
   double t=ttheta;
   t +=  m2ThetaZero +m2ThetaDisplacement/cos(t/2) +m2ThetaTransparency*sin(t);
   return t;
}

long PowderPattern::Get2ThetaCorrPixel(const double ttheta)const
{
   double t=ttheta;
   t +=  m2ThetaZero +m2ThetaDisplacement/cos(t/2) +m2ThetaTransparency*sin(t);
   return (long) ((t-m2ThetaMin)/m2ThetaStep);
}
void PowderPattern::ImportPowderPatternFullprof(const string &filename)
{
   //15.000   0.030  70.000 LANI4FE#1 REC 800ø 4JRS                                 
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
   fin >> m2ThetaMin;
   fin >> m2ThetaStep;
   double tmp;
   fin >> tmp;
   mNbPoint=(long)((tmp-m2ThetaMin)/m2ThetaStep+1.001);
   VFN_DEBUG_MESSAGE("PowderPattern::ImportPowderPatternFullprof() :"\
      << " 2Theta min=" << m2ThetaMin << " 2Theta max=" << tmp \
      << " NbPoints=" << mNbPoint,5)
   m2ThetaMin  *= DEG2RAD;
   m2ThetaStep *= DEG2RAD;
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
   
   fin >> m2ThetaMin;
   fin >> m2ThetaStep;
   double tmp;
   fin >> tmp;
   mNbPoint=(long)((tmp-m2ThetaMin+.001)/m2ThetaStep)+1;
   VFN_DEBUG_MESSAGE("PowderPattern::ImportPowderPatternPSI_DMC() :"\
      << " 2Theta min=" << m2ThetaMin << " 2Theta max=" << tmp \
      << " NbPoints=" << mNbPoint,5)
   m2ThetaMin  *= DEG2RAD;
   m2ThetaStep *= DEG2RAD;
   mPowderPatternObs.resize (mNbPoint);
   mPowderPatternObsSigma.resize (mNbPoint);
   mPowderPatternWeight.resize(mNbPoint);
   
   fin.getline(tmpComment,100);
   //if(""==mName) mName.append(tmpComment);
   
   for(unsigned long i=0;i<mNbPoint;i++) fin >> mPowderPatternObs(i);
   for(unsigned long i=0;i<mNbPoint;i++) fin >> mPowderPatternObsSigma(i);
   fin.close();
	this->SetWeightToInvSigmaSq();
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
   fin >> m2ThetaMin;
   fin >> m2ThetaStep;
   VFN_DEBUG_MESSAGE("PowderPattern::ImportPowderPatternILL_D1AD2B() :"\
      << " 2Theta min=" << m2ThetaMin << " 2Theta max=" << m2ThetaMin+mNbPoint*m2ThetaStep \
      << " NbPoints=" << mNbPoint,5)
   m2ThetaMin  *= DEG2RAD;
   m2ThetaStep *= DEG2RAD;
   mPowderPatternObs.resize (mNbPoint);
   mPowderPatternObsSigma.resize (mNbPoint);
   mPowderPatternWeight.resize(mNbPoint);
   
   //if(""==mName) mName.append(tmpComment);
   
   for(unsigned long i=0;i<mNbPoint;i++) fin >> mPowderPatternObs(i);
   for(unsigned long i=0;i<mNbPoint;i++) fin >> mPowderPatternObsSigma(i);
   fin.close();
	this->SetWeightToInvSigmaSq();
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
   
   fin >> m2ThetaMin;
   fin >> m2ThetaStep;
   double tmp;
   fin >> tmp;
   mNbPoint=(long)((tmp-m2ThetaMin)/m2ThetaStep+1);
   m2ThetaMin  *= DEG2RAD;
   m2ThetaStep *= DEG2RAD;
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
   
   fin >> m2ThetaMin;
   double tmp;
   fin >> tmp;//2Theta Max
   fin >> m2ThetaStep;
   mNbPoint=(long)((tmp-m2ThetaMin)/m2ThetaStep+1);
   m2ThetaMin  *= DEG2RAD;
   m2ThetaStep *= DEG2RAD;
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
   CrystVector_double tmp2Theta(500);
   mNbPoint=0;
   
   do
   {
      fin >> tmp2Theta              (mNbPoint);
      fin >> mPowderPatternObs     (mNbPoint);
      fin >> mPowderPatternObsSigma(mNbPoint);
      //cout << tmp2Theta              (mNbPoint)<<" "
      //     << mPowderPatternObs     (mNbPoint)<<" "
      //    << mPowderPatternObsSigma(mNbPoint)<<endl;
      if(!fin) break;
      mNbPoint++;
      if( (mNbPoint%500)==0)
      {
         tmp2Theta.resizeAndPreserve(mNbPoint+500);
         mPowderPatternObs.resizeAndPreserve(mNbPoint+500);
         mPowderPatternObsSigma.resizeAndPreserve(mNbPoint+500);
      }
   } while(fin.eof() == false);
   fin.close();
   
   tmp2Theta.resizeAndPreserve              (mNbPoint);
   mPowderPatternObs.resizeAndPreserve     (mNbPoint);
   mPowderPatternObsSigma.resizeAndPreserve(mNbPoint);
   mPowderPatternWeight.resize(mNbPoint);
   
   m2ThetaMin=tmp2Theta(0);
   m2ThetaStep=(tmp2Theta(mNbPoint-1)-tmp2Theta(0))
                                 /(mNbPoint-1);
   m2ThetaMin  *= DEG2RAD;
   m2ThetaStep *= DEG2RAD;
	
	this->SetWeightToInvSigmaSq();
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
   CrystVector_double tmp2Theta(500);
   mNbPoint=0;
   
   do
   {
      fin >> tmp2Theta              (mNbPoint);
      fin >> mPowderPatternObs     (mNbPoint);
      mPowderPatternObsSigma(mNbPoint)
         =sqrt(mPowderPatternObs(mNbPoint));
      if(!fin) break;
      mNbPoint++;
      if( (mNbPoint%500)==0)
      {
         tmp2Theta.resizeAndPreserve(mNbPoint+500);
         mPowderPatternObs.resizeAndPreserve(mNbPoint+500);
         mPowderPatternObsSigma.resizeAndPreserve(mNbPoint+500);
      }
   } while(fin.eof() == false);
   fin.close();
   
   tmp2Theta.resizeAndPreserve              (mNbPoint);
   mPowderPatternObs.resizeAndPreserve     (mNbPoint);
   mPowderPatternObsSigma.resizeAndPreserve(mNbPoint);
   mPowderPatternWeight.resize(mNbPoint);

   m2ThetaMin=tmp2Theta(0);
   m2ThetaStep=(tmp2Theta(mNbPoint-1)-tmp2Theta(0))/(mNbPoint-1);
   m2ThetaMin  *= DEG2RAD;
   m2ThetaStep *= DEG2RAD;

	this->SetSigmaToSqrtIobs();
	this->SetWeightToInvSigmaSq();

   VFN_DEBUG_MESSAGE("PowderPattern::ImportPowderPattern2ThetaObs():finished",5)
}

void PowderPattern::SetPowderPatternObs(const CrystVector_double& obs)
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
}
void PowderPattern::SavePowderPattern(const string &filename) const
{
   VFN_DEBUG_MESSAGE("PowderPattern::SavePowderPattern",5)
   this->CalcPowderPattern();
   ofstream out(filename.c_str());
   CrystVector_double ttheta(mNbPoint);
   //double *p=ttheta.data();
   for(unsigned long i=0;i<mNbPoint;i++) 
      ttheta(i) = m2ThetaMin+ i * m2ThetaStep;
   ttheta *= RAD2DEG;
   
   CrystVector_double diff;
   diff=mPowderPatternObs;
   diff-=mPowderPatternCalc;
   out << "#    2Theta     Iobs       ICalc   Iobs-Icalc    Weight  Comp0" << endl;
   out << FormatVertVector<double>(ttheta,
                    mPowderPatternObs,
                    mPowderPatternCalc,
                    diff,mPowderPatternWeight,
                    mPowderPatternComponentRegistry.GetObj(0).mPowderPatternCalc,12,4);
   out.close();
   VFN_DEBUG_MESSAGE("DiffractionDataPowder::SavePowderPattern:End",3)
}

void PowderPattern::PrintObsCalcData(ostream&os)const
{
   VFN_DEBUG_MESSAGE("DiffractionDataPowder::PrintObsCalcData()",5);
   CrystVector_double ttheta(mNbPoint);
   //double *p=ttheta.data();
   for(unsigned long i=0;i<mNbPoint;i++) 
      ttheta(i) = m2ThetaMin+ i * m2ThetaStep;
   ttheta *= RAD2DEG;
   os << "PowderPattern : " << mName <<endl;
   os << "      2Theta      Obs          Sigma        Calc        Weight" <<endl;
   os << FormatVertVector<double>(ttheta,mPowderPatternObs,mPowderPatternObsSigma,
               mPowderPatternCalc,mPowderPatternWeight,12,4);
   //            mPowderPatternComponentRegistry.GetObj(0).mPowderPatternCalc,12,4);
}

double PowderPattern::GetR()const
{
   this->CalcPowderPattern();
   TAU_PROFILE("PowderPattern::GetR()","void ()",TAU_DEFAULT);
   
   double tmp1=0.;
   double tmp2=0.;
   
   if(false==mStatisticsUseIntegratedPeak)
   {
      long maxPoints=mNbPoint;
      if(true==mUseOnlyLowAngleData) 
         maxPoints= (long)((2*mUseOnlyLowAngleDataLimit-m2ThetaMin)
                        /m2ThetaStep+1);

      if(  (true==mStatisticsExcludeBackground)
         &&(mPowderPatternBackgroundCalc.numElements()>0))
      {
         const double *p1, *p2, *p3;
         p1=mPowderPatternCalc.data();
         p2=mPowderPatternObs.data();
         p3=mPowderPatternBackgroundCalc.data();
         const long nbExclude=mExcludedRegionMin2Theta.numElements();
         if(0==nbExclude)
         {
            VFN_DEBUG_MESSAGE("PowderPattern::GetR():Exclude Backgd",4);
            for(long i=0;i<maxPoints;i++)
            {
               tmp1 += ((*p1)-(*p2)) * ((*p1)-(*p2));
               tmp2 += ((*p2)-(*p3)) * ((*p2)-(*p3));
               p1++;p2++;p3++;
            }
         }
         else
         {
            VFN_DEBUG_MESSAGE("PowderPattern::GetR():Exclude Backgd,Exclude regions",4);
            long min,max;
            long i=0;
            for(int j=0;j<nbExclude;j++)
            {
               min=(long)floor((mExcludedRegionMin2Theta(j)-m2ThetaMin)
                                    /m2ThetaStep);
               max=(long)ceil((mExcludedRegionMax2Theta(j)-m2ThetaMin)
                                    /m2ThetaStep);
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
         const double *p1, *p2;
         p1=mPowderPatternCalc.data();
         p2=mPowderPatternObs.data();
         const long nbExclude=mExcludedRegionMin2Theta.numElements();
         if(0==nbExclude)
         {
            VFN_DEBUG_MESSAGE("PowderPattern::GetR()",4);
            for(unsigned long i=0;i<mNbPoint;i++)
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
            long min,max;
            long i=0;
            for(int j=0;j<nbExclude;j++)
            {
               min=(long)floor((mExcludedRegionMin2Theta(j)-m2ThetaMin)
                                    /m2ThetaStep);
               max=(long)ceil((mExcludedRegionMax2Theta(j)-m2ThetaMin)
                                    /m2ThetaStep);
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
   }
   else //Use integrated peaks ?
   { 
      throw ObjCrystException("PowderPattern::GetR(): Unimplemented yet !");
      /*
      long nbRefl=this->NbRefl();
      double theta;
      long thetaPt,first,last;
      double fwhm;
      {
         double tmp = m2ThetaMin+m2ThetaStep*mNbPoint;
         tmp/= 2;
         fwhm=sqrt(mPowderCagliotiW + mPowderCagliotiV*tmp + mPowderCagliotiU*tmp*tmp);
      }
      const long nbIntegrPt=(long)(fwhm/m2ThetaStep);
      #ifdef __DEBUG__
      if(  (true==mStatisticsExcludeBackground)
         &&(mPowderPatternBackgroundCalc.numElements()>0))
      {
         VFN_DEBUG_MESSAGE("PowderPattern::GetR() Integrated R (2x"<<nbIntegrPt<<
            "points) & Exclude Backgd",3)
      }
      else
      {
         VFN_DEBUG_MESSAGE("PowderPattern::GetR()Integrated R",3)
      }
      #endif
         
      for(long i=0;i<nbRefl;i++)
      {
         theta=mTheta(i);
         if((true==mUseOnlyLowAngleData) && (theta>mUseOnlyLowAngleDataLimit)) continue;
         theta +=  m2ThetaZero/2.
                  +m2ThetaDisplacement/cos(theta)
                  +m2ThetaTransparency*sin(2*theta);
         thetaPt =(long) ((2*theta-(m2ThetaMin))
                           /m2ThetaStep);
         
         first=thetaPt-nbIntegrPt;
         if( first >= mNbPoint) continue;
         if( first < 0) first =0;
         last=thetaPt+nbIntegrPt;
         if( last >= mNbPoint) last=mNbPoint;
         
         double ttmp1 = 0;
         double ttmp2 = 0;
         double ttmp3 = 0;
         if((true==mStatisticsExcludeBackground)&&(true==mHasPowderPatternBackground))
         {
            const double *p1 = mPowderPatternCalc.data()+first;
            const double *p2 = mPowderPatternObs.data()+first;
            const double *p3 = mPowderPatternBackgroundCalc.data()+first;
            for(long j=first;j<last;j++)
            {
               ttmp1 += *p1++;
               ttmp2 += *p2++;
               ttmp3 += *p3++;
            }
         }
         else
         {
            const double *p1 = mPowderPatternCalc.data()+first;
            const double *p2 = mPowderPatternObs.data()+first;
            for(long j=first;j<last;j++)
            {
               ttmp1 += *p1++;
               ttmp2 += *p2++;
            }
         }
         tmp1 += (ttmp1-ttmp2) * (ttmp1-ttmp2);
         tmp2 += (ttmp2-ttmp3) * (ttmp2-ttmp3);
      }
      */
   }
   
   VFN_DEBUG_MESSAGE("PowderPattern::GetR()="<<sqrt(tmp1/tmp2),4);
   //cout << FormatVertVector<double>(mPowderPatternCalc,mPowderPatternObs);
   //this->SavePowderPattern("refinedPattern.out");
   //abort();
   return sqrt(tmp1/tmp2);
}

double PowderPattern::GetRw()const
{
   this->CalcPowderPattern();
   TAU_PROFILE("PowderPattern::GetRw()","void ()",TAU_DEFAULT);
   VFN_DEBUG_MESSAGE("PowderPattern::GetRw()",3);
   
   
   //cout <<FormatVertVector<double>(mPowderPatternObs,
   //                               mPowderPatternCalc,
   //                               mPowderPatternWeight);
   double tmp1=0.;
   double tmp2=0.;
   
   if(false==mStatisticsUseIntegratedPeak)
   {
      long maxPoints=mNbPoint;
      if(true==mUseOnlyLowAngleData) 
         maxPoints= (long)((2*mUseOnlyLowAngleDataLimit-m2ThetaMin)
                        /m2ThetaStep+1);

      if(  (true==mStatisticsExcludeBackground)
         &&(mPowderPatternBackgroundCalc.numElements()>0))
      {
         VFN_DEBUG_MESSAGE("PowderPattern::GetRw():Exclude Backgd",3);
         const double *p1, *p2, *p3, *p4;
         p1=mPowderPatternCalc.data();
         p2=mPowderPatternObs.data();
         p3=mPowderPatternBackgroundCalc.data();
         p4=mPowderPatternWeight.data();
         const long nbExclude=mExcludedRegionMin2Theta.numElements();
         if(0==nbExclude)
         {
            for(long i=0;i<maxPoints;i++)
            {
               tmp1 += *p4   * ((*p1)-(*p2)) * ((*p1)-(*p2));
               tmp2 += *p4++ * ((*p2)-(*p3)) * ((*p2)-(*p3));
               p1++;p2++;p3++;
            }
         }
         else
         {
            long min,max;
            long i=0;
            for(int j=0;j<nbExclude;j++)
            {
               min=(long)floor((mExcludedRegionMin2Theta(j)-m2ThetaMin)
                                    /m2ThetaStep);
               max=(long)ceil((mExcludedRegionMax2Theta(j)-m2ThetaMin)
                                    /m2ThetaStep);
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
         const double *p1, *p2, *p4;
         p1=mPowderPatternCalc.data();
         p2=mPowderPatternObs.data();
         p4=mPowderPatternWeight.data();
         const long nbExclude=mExcludedRegionMin2Theta.numElements();
         if(0==nbExclude)
         {
            for(unsigned long i=0;i<mNbPoint;i++)
            {
               tmp1 += *p4   * ((*p1)-(*p2))*((*p1)-(*p2));
               tmp2 += *p4++ * (*p2) * (*p2);
               p1++;p2++;
            }
         }
         else
         {
            long min,max;
            long i=0;
            for(int j=0;j<nbExclude;j++)
            {
               min=(long)floor((mExcludedRegionMin2Theta(j)-m2ThetaMin)
                                    /m2ThetaStep);
               max=(long)ceil((mExcludedRegionMax2Theta(j)-m2ThetaMin)
                                    /m2ThetaStep);
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
   }
   else throw ObjCrystException("Use of Integrated R-factor for Rw not implemented yet !!");
   VFN_DEBUG_MESSAGE("PowderPattern::GetRw()="<<sqrt(tmp1/tmp2),3);
   return sqrt(tmp1/tmp2);
}

double PowderPattern::GetChiSq()const
{
   TAU_PROFILE("PowderPattern::GetChiSq()","void ()",TAU_DEFAULT);
   VFN_DEBUG_MESSAGE("PowderPattern::GetChiSq()",3);
   
   long maxPoints=mNbPoint;
   if(true==mUseOnlyLowAngleData) 
      maxPoints= (long)((2*mUseOnlyLowAngleDataLimit-m2ThetaMin)
                     /m2ThetaStep+1);
                     
   double tmp1=0.;
   {
      const double *p1, *p2, *p3;
      p1=mPowderPatternCalc.data();
      p2=mPowderPatternObs.data();
      p3=mPowderPatternWeight.data();
      for(long i=0;i<maxPoints;i++)
      {
         tmp1 += *p3 * ((*p1)-(*p2))*((*p1)-(*p2));
         p1++;p2++;
      }
   }
   VFN_DEBUG_MESSAGE("PowderPattern::GetChiSq()="<<tmp1,3);
   return tmp1;
}

void PowderPattern::FitScaleFactorForR()
{
   this->CalcPowderPattern();
   TAU_PROFILE("PowderPattern::FitScaleFactorForR()","void ()",TAU_DEFAULT);
   VFN_DEBUG_ENTRY("PowderPattern::FitScaleFactorForR()",3);
   //:TODO: take into account excluded regions...
   // Which components are scalable ?
      mScalableComponentIndex.resize(mPowderPatternComponentRegistry.GetNb());
      int nbScale=0;
      for(int i=0;i<mPowderPatternComponentRegistry.GetNb();i++)
		{
         if(mPowderPatternComponentRegistry.GetObj(i).IsScalable())
				mScalableComponentIndex(nbScale++)=i;
		}
      mScalableComponentIndex.resizeAndPreserve(nbScale);
   VFN_DEBUG_MESSAGE("-> Number of Scale Factors:"<<nbScale<<":Index:"<<endl<<mScalableComponentIndex,3);
   // prepare matrices
      //mFitScaleFactorD.resize(mNbPoint,nbScale);
      mFitScaleFactorM.resize(nbScale,nbScale);
      mFitScaleFactorB.resize(nbScale,1);
      mFitScaleFactorX.resize(nbScale,1);
   // Build Matrix & Vector for LSQ
   const long nbExclude=mExcludedRegionMin2Theta.numElements();
   if(0==nbExclude)
   {
   	for(int i=0;i<nbScale;i++)
   	{
      	for(int j=i;j<nbScale;j++)
      	{
         	// Here use a direct access to the powder spectrum, since
         	// we know it has just been recomputed
         	const double *p1=mPowderPatternComponentRegistry.GetObj(mScalableComponentIndex(i))
                           	.mPowderPatternCalc.data();
         	const double *p2=mPowderPatternComponentRegistry.GetObj(mScalableComponentIndex(j))
                           	.mPowderPatternCalc.data();
         	double m=0.;
         	for(unsigned long k=0;k<mNbPoint;k++) m += *p1++ * *p2++;
         	mFitScaleFactorM(i,j)=m;
         	mFitScaleFactorM(j,i)=m;
      	}
   	}
   	for(int i=0;i<nbScale;i++)
   	{
      	const double *p1=mPowderPatternObs.data();
      	const double *p2=mPowderPatternComponentRegistry.GetObj(mScalableComponentIndex(i))
                        	.mPowderPatternCalc.data();
      	double b=0.;
			if(mPowderPatternBackgroundCalc.numElements()<=1)
      		for(unsigned long k=0;k<mNbPoint;k++) b += *p1++ * *p2++;
			else
			{
      		const double *p3=mPowderPatternBackgroundCalc.data();
      		for(unsigned long k=0;k<mNbPoint;k++) b += (*p1++ - *p3++) * *p2++;
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
         	const double *p1=mPowderPatternComponentRegistry.GetObj(mScalableComponentIndex(i))
                           	.mPowderPatternCalc.data();
         	const double *p2=mPowderPatternComponentRegistry.GetObj(mScalableComponentIndex(j))
                           	.mPowderPatternCalc.data();
         	double m=0.;
            unsigned long l=0;
            for(int k=0;k<nbExclude;k++)
            {
               min=(unsigned long)floor((mExcludedRegionMin2Theta(k)-m2ThetaMin)
                                    /m2ThetaStep);
               max=(unsigned long)ceil((mExcludedRegionMax2Theta(k)-m2ThetaMin)
                                    /m2ThetaStep);
               if(min>mNbPoint) break;
               if(max>mNbPoint)max=mNbPoint;
					//! min is the *beginning* of the excluded region
               for(;l<min;l++) m += *p1++ * *p2++;
               p1 += max-l;
               p2 += max-l;
               l  = max;
            }
            for(;l<mNbPoint;l++) m += *p1++ * *p2++;
         	mFitScaleFactorM(i,j)=m;
         	mFitScaleFactorM(j,i)=m;
      	}
   	}
   	for(int i=0;i<nbScale;i++)
   	{
      	const double *p1=mPowderPatternObs.data();
      	const double *p2=mPowderPatternComponentRegistry
									.GetObj(mScalableComponentIndex(i))
                        		.mPowderPatternCalc.data();
      	double b=0.;
         unsigned long l=0;
			if(mPowderPatternBackgroundCalc.numElements()<=1)
			{
         	for(int k=0;k<nbExclude;k++)
         	{
            	min=(long)floor((mExcludedRegionMin2Theta(k)-m2ThetaMin)
                                 	/m2ThetaStep);
            	max=(long)ceil((mExcludedRegionMax2Theta(k)-m2ThetaMin)
                                 	/m2ThetaStep);
            	if(min>mNbPoint) break;
            	if(max>mNbPoint)max=mNbPoint;
					//! min is the *beginning* of the excluded region
            	for(;l<min;l++) b += *p1++ * *p2++;
            	p1 += max-l;
            	p2 += max-l;
            	l  = max;
         	}
         	for(;l<mNbPoint;l++) b += *p1++ * *p2++;
			}
			else
			{
      		const double *p3=mPowderPatternBackgroundCalc.data();
         	for(int k=0;k<nbExclude;k++)
         	{
            	min=(long)floor((mExcludedRegionMin2Theta(k)-m2ThetaMin)
                                 	/m2ThetaStep);
            	max=(long)ceil((mExcludedRegionMax2Theta(k)-m2ThetaMin)
                                 	/m2ThetaStep);
            	if(min>mNbPoint) break;
            	if(max>mNbPoint)max=mNbPoint;
					//! min is the *beginning* of the excluded region
            	for(;l<min;l++) b += (*p1++ - *p3++) * *p2++;
            	p1 += max-l;
            	p2 += max-l;
            	l  = max;
         	}
         	for(;l<mNbPoint;l++) b += (*p1++ - *p3++) * *p2++;
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
      const double * p1=mPowderPatternComponentRegistry.GetObj(mScalableComponentIndex(i))
                        .mPowderPatternCalc.data();
      double * p0 = mPowderPatternCalc.data();
      const double s = mFitScaleFactorX(i)
							  -mScaleFactor(mScalableComponentIndex(i));
      for(unsigned long j=0;j<mNbPoint;j++) *p0++ += s * *p1++;
      VFN_DEBUG_MESSAGE("-> Old:"<<mScaleFactor(mScalableComponentIndex(i)) <<" Change:"<<mFitScaleFactorX(i),2);
      mScaleFactor(mScalableComponentIndex(i)) = mFitScaleFactorX(i);
		mClockScaleFactor.Click();
   }
   VFN_DEBUG_EXIT("PowderPattern::FitScaleFactorForR():End",3);
}

void PowderPattern::FitScaleFactorForRw()
{
   this->CalcPowderPattern();
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
      mScalableComponentIndex.resizeAndPreserve(nbScale);
   VFN_DEBUG_MESSAGE("-> Number of Scale Factors:"<<nbScale<<":Index:"<<endl<<mScalableComponentIndex,2);
   // prepare matrices
      //mFitScaleFactorD.resize(mNbPoint,nbScale);
      mFitScaleFactorM.resize(nbScale,nbScale);
      mFitScaleFactorB.resize(nbScale,1);
      mFitScaleFactorX.resize(nbScale,1);
   // Build Matrix & Vector for LSQ
   const long nbExclude=mExcludedRegionMin2Theta.numElements();
   if(0==nbExclude)
   {
   	for(int i=0;i<nbScale;i++)
   	{
      	for(int j=i;j<nbScale;j++)
      	{
         	// Here use a direct access to the powder spectrum, since
         	// we know it has just been recomputed
         	const double *p1=mPowderPatternComponentRegistry.GetObj(mScalableComponentIndex(i))
                           	.mPowderPatternCalc.data();
         	const double *p2=mPowderPatternComponentRegistry.GetObj(mScalableComponentIndex(j))
                           	.mPowderPatternCalc.data();
         	const double *p3=mPowderPatternWeight.data();
         	double m=0.;
         	for(unsigned long k=0;k<mNbPoint;k++) m += *p1++ * *p2++ * *p3++;
         	mFitScaleFactorM(i,j)=m;
         	mFitScaleFactorM(j,i)=m;
      	}
   	}
   	for(int i=0;i<nbScale;i++)
   	{
      	const double *p1=mPowderPatternObs.data();
      	const double *p2=mPowderPatternComponentRegistry.GetObj(mScalableComponentIndex(i))
                        	.mPowderPatternCalc.data();
      	const double *p3=mPowderPatternWeight.data();
      	double b=0.;
			if(mPowderPatternBackgroundCalc.numElements()<=1)
      		for(unsigned long k=0;k<mNbPoint;k++) b += *p1++ * *p2++ * *p3++;
			else
			{
      		const double *p4=mPowderPatternBackgroundCalc.data();
      		for(unsigned long k=0;k<mNbPoint;k++)
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
         	const double *p1=mPowderPatternComponentRegistry.GetObj(mScalableComponentIndex(i))
                           	.mPowderPatternCalc.data();
         	const double *p2=mPowderPatternComponentRegistry.GetObj(mScalableComponentIndex(j))
                           	.mPowderPatternCalc.data();
         	const double *p3=mPowderPatternWeight.data();
         	double m=0.;
            unsigned long l=0;
            for(int k=0;k<nbExclude;k++)
            {
               min=(unsigned long)floor((mExcludedRegionMin2Theta(k)-m2ThetaMin)
                                    /m2ThetaStep);
               max=(unsigned long)ceil((mExcludedRegionMax2Theta(k)-m2ThetaMin)
                                    /m2ThetaStep);
               if(min>mNbPoint) break;
               if(max>mNbPoint)max=mNbPoint;
					//! min is the *beginning* of the excluded region
               for(;l<min;l++) m += *p1++ * *p2++ * *p3++;
               p1 += max-l;
               p2 += max-l;
               p3 += max-l;
               l  = max;
            }
            for(;l<mNbPoint;l++) m += *p1++ * *p2++ * *p3++;
         	mFitScaleFactorM(i,j)=m;
         	mFitScaleFactorM(j,i)=m;
      	}
   	}
   	for(int i=0;i<nbScale;i++)
   	{
      	const double *p1=mPowderPatternObs.data();
      	const double *p2=mPowderPatternComponentRegistry
									.GetObj(mScalableComponentIndex(i))
                        		.mPowderPatternCalc.data();
         const double *p3=mPowderPatternWeight.data();
      	double b=0.;
         unsigned long l=0;
			if(mPowderPatternBackgroundCalc.numElements()<=1)
			{
         	for(int k=0;k<nbExclude;k++)
         	{
            	min=(long)floor((mExcludedRegionMin2Theta(k)-m2ThetaMin)
                                 	/m2ThetaStep);
            	max=(long)ceil((mExcludedRegionMax2Theta(k)-m2ThetaMin)
                                 	/m2ThetaStep);
            	if(min>mNbPoint) break;
            	if(max>mNbPoint)max=mNbPoint;
					//! min is the *beginning* of the excluded region
            	for(;l<min;l++) b += *p1++ * *p2++ * *p3++;
            	p1 += max-l;
            	p2 += max-l;
            	p3 += max-l;
            	l  = max;
         	}
         	for(;l<mNbPoint;l++) b += *p1++ * *p2++ * *p3++;
			}
			else
			{
      		const double *p4=mPowderPatternBackgroundCalc.data();
         	for(int k=0;k<nbExclude;k++)
         	{
            	min=(long)floor((mExcludedRegionMin2Theta(k)-m2ThetaMin)
                                 	/m2ThetaStep);
            	max=(long)ceil((mExcludedRegionMax2Theta(k)-m2ThetaMin)
                                 	/m2ThetaStep);
            	if(min>mNbPoint) break;
            	if(max>mNbPoint)max=mNbPoint;
					//! min is the *beginning* of the excluded region
            	for(;l<min;l++) b += (*p1++ - *p4++) * *p2++ * *p3++;
            	p1 += max-l;
            	p2 += max-l;
            	p3 += max-l;
            	l  = max;
         	}
         	for(;l<mNbPoint;l++) b += (*p1++ - *p4++) * *p2++ * *p3++;
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
      const double * p1=mPowderPatternComponentRegistry.GetObj(mScalableComponentIndex(i))
                        .mPowderPatternCalc.data();
      double * p0 = mPowderPatternCalc.data();
      const double s = mFitScaleFactorX(i)
							  -mScaleFactor(mScalableComponentIndex(i));
      for(unsigned long j=0;j<mNbPoint;j++) *p0++ += s * *p1++;
      VFN_DEBUG_MESSAGE("-> Old:"<<mScaleFactor(mScalableComponentIndex(i)) <<" Change:"<<mFitScaleFactorX(i),2);
      mScaleFactor(mScalableComponentIndex(i)) = mFitScaleFactorX(i);
		mClockScaleFactor.Click();
   }
   VFN_DEBUG_EXIT("PowderPattern::FitScaleFactorForRw():End",3);
}

void PowderPattern::SetSigmaToSqrtIobs()
{
   VFN_DEBUG_MESSAGE("PowderPattern::SetSigmaToSqrtIobs()",5);
   for(long i=0;i<mPowderPatternObs.numElements();i++)
      mPowderPatternObsSigma(i)=sqrt(fabs(mPowderPatternObs(i)));
}
void PowderPattern::SetWeightToInvSigmaSq(const double minRelatSigma)
{
   VFN_DEBUG_MESSAGE("PowderPattern::SetWeightToInvSigmaSq()",5);
   //:KLUDGE: If less than 1e-4*max, set to 0.... Do not give weight to unobserved points
   const double min=MaxAbs(mPowderPatternObsSigma)*minRelatSigma;
   //mPowderPatternWeight.resize(mPowderPatternObsSigma.numElements());
   double tmp;
   for(long i=0;i<mPowderPatternObsSigma.numElements();i++)
   {
      tmp = mPowderPatternObsSigma(i);
      if(tmp<min) mPowderPatternWeight(i)= 1./min/min; 
      else  mPowderPatternWeight(i) =1./tmp/tmp;
   }
}
void PowderPattern::SetWeightToSinTheta(const double power)
{
   VFN_DEBUG_MESSAGE("PowderPattern::SetWeightToSinTheta()",5);
   //mPowderPatternWeight.resize(mPowderPatternObs.numElements());
   for(long i=0;i<mPowderPatternWeight.numElements();i++) 
      mPowderPatternWeight(i) =pow(sin((m2ThetaMin+ i * m2ThetaStep)/2),power);
}
void PowderPattern::SetWeightToUnit()
{
   VFN_DEBUG_MESSAGE("PowderPattern::SetWeightToSinTheta()",5);
   //mPowderPatternWeight.resize(mPowderPatternObs.numElements());
   mPowderPatternWeight=1;
}
void PowderPattern::SetWeightPolynomial(const double a, const double b,
													 const double c,
												 	 const double minRelatIobs)
{
   VFN_DEBUG_MESSAGE("PowderPattern::SetWeightPolynomial()",5);
   const double min=MaxAbs(mPowderPatternObs)*minRelatIobs;
	double tmp;
   for(long i=0;i<mPowderPatternWeight.numElements();i++)
	{
		tmp=mPowderPatternObs(i);
		if(tmp<min)
			mPowderPatternWeight(i) =   1./(a+min+b*min*min+c*min*min*min);
		else mPowderPatternWeight(i) = 1./(a+tmp+b*tmp*tmp+c*tmp*tmp*tmp);
	}
}

void PowderPattern::Add2ThetaExcludedRegion(const double min2Theta,const double max2Theta)
{
   VFN_DEBUG_MESSAGE("PowderPattern::Add2ThetaExcludedRegion()",5)
   const int num=mExcludedRegionMin2Theta.numElements();
   if(num>0)
   {
      mExcludedRegionMin2Theta.resizeAndPreserve(num+1);
      mExcludedRegionMax2Theta.resizeAndPreserve(num+1);
   }
   else
   {
      mExcludedRegionMin2Theta.resize(1);
      mExcludedRegionMax2Theta.resize(1);
   }
   mExcludedRegionMin2Theta(num)=min2Theta;
   mExcludedRegionMax2Theta(num)=max2Theta;
   
   //ensure regions are sorted by ascending 2thetamin
   CrystVector_long subs;
   subs=SortSubs(mExcludedRegionMin2Theta);
   CrystVector_double tmp1,tmp2;
   tmp1=mExcludedRegionMin2Theta;
   tmp2=mExcludedRegionMax2Theta;
   for(int i=0;i<mExcludedRegionMin2Theta.numElements();i++)
   {
      mExcludedRegionMin2Theta(i)=tmp1(subs(i));
      mExcludedRegionMax2Theta(i)=tmp2(subs(i));
   }
   VFN_DEBUG_MESSAGE(FormatVertVector<double>(mExcludedRegionMin2Theta,mExcludedRegionMax2Theta),10)
   VFN_DEBUG_MESSAGE("PowderPattern::Add2ThetaExcludedRegion():End",5)
}

unsigned int PowderPattern::GetNbCostFunction()const {return 2;}

const string& PowderPattern::GetCostFunctionName(const unsigned int id)const
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
         cout << "RefinableObj::GetCostFunctionName(): Not Found !" <<endl;
         throw 0;
      }
   }
}

const string& PowderPattern::GetCostFunctionDescription(const unsigned int id)const
{
   static string costFunctionDescription[2];
   if(0==costFunctionDescription[0].length())
   {
      costFunctionDescription[0]="unweighted R-factor (best scale)";
      costFunctionDescription[1]="weigthed R-factor (best scale)";
   }
   switch(id)
   {
      case 0: return costFunctionDescription[0];
      case 1: return costFunctionDescription[1];
      default:
      {
         cout << "RefinableObj::GetCostFunctionDescription(): Not Found !" <<endl;
         throw 0;
      }
   }
}

double PowderPattern::GetCostFunctionValue(const unsigned int n)
{
   switch(n)
   {
      case 0:
      {
         this->FitScaleFactorForR()  ;
         return this->GetR();
      }
      case 1: this->FitScaleFactorForRw() ;return this->GetRw();
      default:
      {
         cout << "RefinableObj::GetCostFunctionValue(): Not Found !" <<endl;
         throw 0;
      }
   }
}

void PowderPattern::Prepare()
{
   VFN_DEBUG_MESSAGE("PowderPattern::Prepare()",5);
   for(int i=0;i<mPowderPatternComponentRegistry.GetNb();i++)
      mPowderPatternComponentRegistry.GetObj(i).Prepare();
}

void PowderPattern::CalcPowderPattern() const
{
   TAU_PROFILE("PowderPattern::CalcPowderPattern()","void ()",TAU_DEFAULT);
   VFN_DEBUG_MESSAGE("PowderPattern::CalcPowderPattern()",3);
   if(mPowderPatternComponentRegistry.GetNb()==0)
   {
      mPowderPatternCalc.resize(mNbPoint);
      mPowderPatternCalc=0;
      return;
   }
   TAU_PROFILE_TIMER(timer1,"PowderPattern::CalcPowderPattern1(Calc spectrum components)"\
                     ,"", TAU_FIELD);
   TAU_PROFILE_TIMER(timer2,"PowderPattern::CalcPowderPattern2(Add spectrums)"\
                     ,"", TAU_FIELD);
   TAU_PROFILE_START(timer1);
   for(int i=0;i<mPowderPatternComponentRegistry.GetNb();i++)
      mPowderPatternComponentRegistry.GetObj(i).CalcPowderPattern();
   VFN_DEBUG_MESSAGE("PowderPattern::CalcPowderPattern():Calculated components..",3);
   bool b=false;
	if(mClockPowderPatternCalc<mClockScaleFactor) b=true;
	else
   	for(int i=0;i<mPowderPatternComponentRegistry.GetNb();i++) 
      	if(mClockPowderPatternCalc<
            	mPowderPatternComponentRegistry.GetObj(i).GetClockPowderPatternCalc())
      	{
         	b=true;
         	break;
      	}
   TAU_PROFILE_STOP (timer1);
      
   if(true==b)
   {
      TAU_PROFILE_START(timer2);
      mPowderPatternCalc.resize(mNbPoint);
      mPowderPatternCalc=0;
      mPowderPatternBackgroundCalc.resize(mNbPoint);
      mPowderPatternBackgroundCalc=0;
		int nbBackgd=0;//count number of background phases
      for(int i=0;i<mPowderPatternComponentRegistry.GetNb();i++)
      {
         VFN_DEBUG_MESSAGE("PowderPattern::CalcPowderPattern():Adding "<< mPowderPatternComponentRegistry.GetObj(i).GetName(),3);
         const double * p1=mPowderPatternComponentRegistry.GetObj(i)
                              .mPowderPatternCalc.data();
         double * p0 = mPowderPatternCalc.data();
         if(true==mPowderPatternComponentRegistry.GetObj(i).IsScalable())
         {
            const double s = mScaleFactor(i);
            for(unsigned long j=0;j<mNbPoint;j++) *p0++ += s * *p1++;
         }
         else
			{
				nbBackgd++;
				for(unsigned long j=0;j<mNbPoint;j++) *p0++ += *p1++;
				
         	p0 = mPowderPatternBackgroundCalc.data();
				p1=mPowderPatternComponentRegistry.GetObj(i)
                              .mPowderPatternCalc.data();
				for(unsigned long j=0;j<mNbPoint;j++) *p0++ += *p1++;
			}
      }
		if(0==nbBackgd) mPowderPatternBackgroundCalc.resize(1);//:KLUDGE:
      mClockPowderPatternCalc.Click();
      TAU_PROFILE_STOP (timer2);
   }
   VFN_DEBUG_MESSAGE("PowderPattern::CalcPowderPattern():End",3);
}
void PowderPattern::Init()
{
   {
      RefinablePar tmp("2Theta0",&m2ThetaZero,-.05,.05,gpRefParTypeScattDataCorrPos,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,RAD2DEG);
      tmp.AssignClock(mClockPowderPattern2ThetaCorr);
      tmp.SetDerivStep(1e-4);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("2ThetaDispl",&m2ThetaDisplacement,-.05,.05,gpRefParTypeScattDataCorrPos,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,RAD2DEG);
      tmp.AssignClock(mClockPowderPattern2ThetaCorr);
      tmp.SetDerivStep(1e-4);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("2ThetaTransp",&m2ThetaTransparency,-.05,.05,gpRefParTypeScattDataCorrPos,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,RAD2DEG);
      tmp.AssignClock(mClockPowderPattern2ThetaCorr);
      tmp.SetDerivStep(1e-4);
      this->AddPar(tmp);
   }
}

#ifdef __WX__CRYST__
WXCrystObjBasic* PowderPattern::WXCreate(wxWindow* parent)
{
   //:TODO: Check mpWXCrystObj==0
   mpWXCrystObj=new WXPowderPattern(parent,this);
   return mpWXCrystObj;
}
#endif
//######################################################################
//    PROFILE FUNCTIONS (for powder diffraction)
//######################################################################

CrystVector_double PowderProfileGauss  (const CrystVector_double ttheta,const double fwhm,
                                      const double asymmetryPar)
{
   TAU_PROFILE("PowderProfileGauss()","Vector (Vector,double)",TAU_DEFAULT);
   //:TODO: faster... What if we return to using blitz++ ??
   const long nbPoints=ttheta.numElements();
   CrystVector_double result(nbPoints);
   result=ttheta;
   result *= result;
   double *p=result.data();
   
   if( fabs(asymmetryPar-1.) < 1e-5)
   {
      //reference: IUCr Monographs on Crystallo 5 - The Rietveld Method (ed RA Young)
      result *= -4.*log(2.)/fwhm/fwhm;
   }
   else
   {
      //reference : Toraya J. Appl. Cryst 23(1990),485-491
      
      //Search values <0
      long middlePt;
      const double *pt=ttheta.data();
      for( middlePt=0;middlePt<nbPoints;middlePt++) if( *pt++ > 0) break;
      
      const double c1= -(1.+asymmetryPar)/asymmetryPar*log(2.)/fwhm/fwhm;
      const double c2= -(1.+asymmetryPar)             *log(2.)/fwhm/fwhm;
      for(long i=0;i<middlePt;i++) *p++ *= c1;
      for(long i=middlePt;i<nbPoints;i++) *p++ *= c2;
   }
   p=result.data();
   for(long i=0;i<nbPoints;i++) { *p = exp(*p) ; p++ ;}
   
   result *= 2. / fwhm * sqrt(log(2.)/M_PI);
   return result;
}

CrystVector_double PowderProfileLorentz(const CrystVector_double ttheta,const double fwhm,
                                      const double asymmetryPar)
{
   TAU_PROFILE("PowderProfileLorentz()","Vector (Vector,double)",TAU_DEFAULT);
   //:TODO: faster... What if we return to using blitz++ ??
   //reference: IUCr Monographs on Crystallo 5 - The Rietveld Method (ed RA Young)
   const long nbPoints=ttheta.numElements();
   CrystVector_double result(nbPoints);
   result=ttheta;
   result *= result;
   double *p=result.data();
   if( fabs(asymmetryPar-1.) < 1e-5)
   {
      //reference: IUCr Monographs on Crystallo 5 - The Rietveld Method (ed RA Young)
      result *= 4./fwhm/fwhm;
   }
   else
   {
      //reference : Toraya J. Appl. Cryst 23(1990),485-491
      //Search values <0
      long middlePt;
      const double *pt=ttheta.data();
      for( middlePt=0;middlePt<nbPoints;middlePt++) if( *pt++ > 0) break;
      //cout << nbPoints << " " << middlePt <<endl;
      const double c1= (1+asymmetryPar)/asymmetryPar*(1+asymmetryPar)/asymmetryPar/fwhm/fwhm;
      const double c2= (1+asymmetryPar)*(1+asymmetryPar)/fwhm/fwhm;
      for(long i=0;i<middlePt;i++) *p++ *= c1 ;
      for(long i=middlePt;i<nbPoints;i++)  *p++ *= c2 ;
   }
   p=result.data();
   result += 1. ;
   for(long i=0;i<nbPoints;i++) { *p = 1/(*p) ; p++ ;}
   
   
   result *= 2./M_PI/fwhm;
   return result;
}


}//namespace ObjCryst
