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
*  source file for LibCryst++ ScatteringData class
*
*/

#include <cmath>

#include <typeinfo>

#include "ObjCryst/ScatteringData.h"
#include "Quirks/VFNDebug.h"
#include "Quirks/VFNStreamFormat.h"

#ifdef __WX__CRYST__
   #include "wxCryst/wxPowderPattern.h"
#endif

#include <fstream>
#include <iomanip>

namespace ObjCryst
{

extern "C"
{
//XRay Tubes wavelengths from atominfo
//AtomInfo (c) 1994-96 Ralf W. Grosse-Kunstleve 
#include "atominfo/atominfo.h"
}
const RefParType *gpRefParTypeScattData
   = new RefParType(gpRefParTypeObjCryst,"Scattering Data");
const RefParType *gpRefParTypeScattDataScale
   = new RefParType(gpRefParTypeObjCryst,"Scale Factor");
const RefParType *gpRefParTypeScattDataProfile
   = new RefParType(gpRefParTypeScattData,"Profile");
const RefParType *gpRefParTypeScattDataProfileType
   = new RefParType(gpRefParTypeScattDataProfile,"Type");
const RefParType *gpRefParTypeScattDataProfileWidth
   = new RefParType(gpRefParTypeScattDataProfile,"Width");
const RefParType *gpRefParTypeScattDataProfileAsym
   = new RefParType(gpRefParTypeScattDataProfile,"Asymmetry");
const RefParType *gpRefParTypeScattDataCorr
   = new RefParType(gpRefParTypeScattData,"Correction");
const RefParType *gpRefParTypeScattDataCorrInt
   = new RefParType(gpRefParTypeScattDataCorr,"Intensities");
const RefParType *gpRefParTypeScattDataCorrIntAbsorp
   = new RefParType(gpRefParTypeScattDataCorrInt,"Absorption");
const RefParType *gpRefParTypeScattDataCorrIntPolar
   = new RefParType(gpRefParTypeScattDataCorrInt,"Polarization");
const RefParType *gpRefParTypeScattDataCorrIntExtinc
   = new RefParType(gpRefParTypeScattDataCorrInt,"Extinction");
const RefParType *gpRefParTypeScattDataCorrPos
   = new RefParType(gpRefParTypeScattDataCorr,"Reflections Positions");
const RefParType *gpRefParTypeScattDataBackground
   = new RefParType(gpRefParTypeScattData,"Background");

const RefParType *gpRefParTypeRadiation
   = new RefParType(gpRefParTypeObjCryst,"Radiation");
const RefParType *gpRefParTypeRadiationWavelength
   = new RefParType(gpRefParTypeRadiation,"Wavelength");
//######################################################################
//    Tabulated math functions for faster (&less precise) F(hkl) calculation
//These function are defined and used in cristallo-spacegroup.cpp
//Currently tabulating sine and cosine only
//######################################################################

static bool sLibCrystTabulCosineIsInit=false;
//conversion value
static double sLibCrystTabulCosineRatio;
//initialize tabulated values of cosine
void InitLibCrystTabulCosine();
// Number of tabulated values of cosine between [0;2pi]
// 100 000 is far enough for a model search, yielding a maximum
// error less than .05%
static const long sLibCrystNbTabulSine=100000;
//storage of tabulated values of cosine and sine
static double *spLibCrystTabulCosine;
static double *spLibCrystTabulSine;

void InitLibCrystTabulCosine()
{
   VFN_DEBUG_MESSAGE("InitLibCrystTabulCosine()",10)
   spLibCrystTabulCosine=new double[sLibCrystNbTabulSine];
   spLibCrystTabulSine=new double[sLibCrystNbTabulSine];
   double *tmp=spLibCrystTabulCosine;
   sLibCrystTabulCosineRatio=sLibCrystNbTabulSine/2./M_PI;
   for(double i=0;i<sLibCrystNbTabulSine;i++) *tmp++ = cos(i/sLibCrystTabulCosineRatio);
   tmp=spLibCrystTabulSine;
   for(double i=0;i<sLibCrystNbTabulSine;i++) *tmp++ = sin(i/sLibCrystTabulCosineRatio);
}

//:KLUDGE: The allocated memory for cos and sin table is never freed...
// This should be done after the last ScatteringData object is deleted.

////////////////////////////////////////////////////////////////////////
//
//    Radiation
//
////////////////////////////////////////////////////////////////////////
Radiation::Radiation():
mWavelength(1),mXRayTubeName(""),mXRayTubeDeltaLambda(0.),
mXRayTubeAlpha2Alpha1Ratio(0.5)
{
   mWavelength=1;
   this->InitOptions();
   mRadiationType.SetChoice(RAD_XRAY);
   mWavelengthType.SetChoice(WAVELENGTH_MONOCHROMATIC);
}

Radiation::Radiation(const RadiationType rad,const double wavelength)
{
   this->InitOptions();
   mRadiationType.SetChoice(rad);
   mWavelengthType.SetChoice(WAVELENGTH_MONOCHROMATIC);
   mWavelength=wavelength;
   mXRayTubeName="";
   mXRayTubeDeltaLambda=0.;//useless here
   mXRayTubeAlpha2Alpha1Ratio=0.5;//useless here
}

Radiation::Radiation(const string &XRayTubeElementName,const double alpha2Alpha2ratio)
{
   this->InitOptions();
   this->SetWavelength(XRayTubeElementName,alpha2Alpha2ratio);
}

Radiation::Radiation(const Radiation &old):
mRadiationType(old.mRadiationType),
mWavelengthType(old.mWavelengthType),
mWavelength(old.mWavelength),
mXRayTubeName(old.mXRayTubeName),
mXRayTubeDeltaLambda(old.mXRayTubeDeltaLambda),
mXRayTubeAlpha2Alpha1Ratio(old.mXRayTubeAlpha2Alpha1Ratio)
{
   mClockWavelength.Click();
}

Radiation::~Radiation()
{}

const string Radiation::GetClassName() const {return "Radiation";}

void Radiation::operator=(const Radiation &old)
{
   mRadiationType             =old.mRadiationType;
   mWavelengthType            =old.mWavelengthType;
   mWavelength                =old.mWavelength;
   mXRayTubeName              =old.mXRayTubeName;
   mXRayTubeDeltaLambda       =old.mXRayTubeDeltaLambda;
   mXRayTubeAlpha2Alpha1Ratio =old.mXRayTubeAlpha2Alpha1Ratio;
   mClockWavelength.Click();
   mRadiationType.SetChoice(old.mRadiationType.GetChoice());
}

RadiationType Radiation::GetRadiationType()const 
{return (RadiationType) mRadiationType.GetChoice();}

void Radiation::SetRadiationType(const RadiationType rad)
{
   mRadiationType.SetChoice(rad);
}

WavelengthType Radiation::GetWavelengthType()const
{return (WavelengthType) mWavelengthType.GetChoice();}

const CrystVector_double& Radiation::GetWavelength() const {return mWavelength;}
void Radiation::SetWavelength(const double l)
{
   mWavelength.resize(1);
   mWavelength=l;
   mClockWavelength.Click();
}
void Radiation::SetWavelength(const string &XRayTubeElementName,
                              const double alpha2Alpha2ratio)
{
   mXRayTubeName=XRayTubeElementName;
   mRadiationType.SetChoice(RAD_XRAY);
   mWavelength.resize(1);
   
   if(XRayTubeElementName.length() >=3) //:KLUDGE:
   {
      mWavelengthType.SetChoice(WAVELENGTH_MONOCHROMATIC);
      const T_ChXrayWaveLength *xrayWaveLength;
      xrayWaveLength=ChXrayWaveLengthOf(mXRayTubeName.c_str());
      mWavelength=xrayWaveLength->Length;
   }
   else
   {
      mWavelengthType.SetChoice(WAVELENGTH_ALPHA12);
      mXRayTubeAlpha2Alpha1Ratio=alpha2Alpha2ratio;
      const T_ChXrayWaveLength *xrayWaveLength;
      double lambda1,lambda2;
      xrayWaveLength=ChXrayWaveLengthOf((mXRayTubeName+"A1").c_str());
      lambda1=xrayWaveLength->Length;
      xrayWaveLength=ChXrayWaveLengthOf((mXRayTubeName+"A2").c_str());
      lambda2=xrayWaveLength->Length;
      mXRayTubeDeltaLambda=lambda2-lambda1;
      mWavelength=lambda1
            +mXRayTubeDeltaLambda*mXRayTubeAlpha2Alpha1Ratio/(1.+mXRayTubeAlpha2Alpha1Ratio);
   }
   mClockWavelength.Click();
}

double Radiation::GetXRayTubeDeltaLambda()const {return mXRayTubeDeltaLambda;}

double Radiation::GetXRayTubeAlpha2Alpha1Ratio()const {return mXRayTubeAlpha2Alpha1Ratio;}
      

const RefinableObjClock& Radiation::GetClockWavelength() const {return mClockWavelength;}
const RefinableObjClock& Radiation::GetClockRadiation()const {return mRadiationType.GetClock();}

void Radiation::Print()const
{
   VFN_DEBUG_MESSAGE("Radiation::Print():"<<this->GetName(),5)
   cout << "Radiation:" << " " ;
   
   switch(mRadiationType.GetChoice())
   {
      case RAD_NEUTRON:  cout<< "Neutron,";break;
      case RAD_XRAY:     cout<< "X-Ray,";break;
      case RAD_ELECTRON: cout<< "Electron,";break;
   }
      
   cout << "Wavelength=" <<" ";
   switch(mWavelengthType.GetChoice())
   {
      case WAVELENGTH_MONOCHROMATIC: cout<< "monochromatic:"<<" "<<mWavelength(0) <<endl;break;
      case WAVELENGTH_ALPHA12:     cout  << "tube:"<<" "<<mXRayTubeName<<", Alpha1/Alpha2= "
                                       << mXRayTubeAlpha2Alpha1Ratio<<endl;break;
      case WAVELENGTH_MAD: cout<< "mad"<<" "<<endl;break;
      case WAVELENGTH_DAFS: cout<< "dafs"<<" "<<endl;break;
      case WAVELENGTH_LAUE: cout<< "laue"<<" "<<endl;break;
   }
}

void Radiation::InitOptions()
{
   static string RadiationTypeName;
   static string RadiationTypeChoices[2];
   static string WavelengthTypeName;
   static string WavelengthTypeChoices[5];
   
   static bool needInitNames=true;
   if(true==needInitNames)
   {
      RadiationTypeName="Radiation";
      RadiationTypeChoices[0]="Neutron";
      RadiationTypeChoices[1]="X-Ray";
      
      WavelengthTypeName="Spectrum";
      WavelengthTypeChoices[0]="Monochromatic";
      WavelengthTypeChoices[1]="X-Ray Tube";
      WavelengthTypeChoices[2]="MAD";
      WavelengthTypeChoices[3]="DAFS";
      WavelengthTypeChoices[4]="LAUE";
      
      needInitNames=false;//Only once for the class
   }
   mRadiationType.Init(2,&RadiationTypeName,RadiationTypeChoices);
   mWavelengthType.Init(5,&WavelengthTypeName,WavelengthTypeChoices);
   this->AddOption(&mRadiationType);
   this->AddOption(&mWavelengthType);
   
   {//Fixed by default
      RefinablePar tmp("Wavelength",mWavelength.data(),0.05,20.,
                        gpRefParTypeRadiationWavelength,REFPAR_DERIV_STEP_ABSOLUTE,
                        true,true,true,false,1.0);
      tmp.SetDerivStep(1e-4);
      tmp.AssignClock(mClockWavelength);
      this->AddPar(tmp);
   }
}

#ifdef __WX__CRYST__
WXCrystObjBasic* Radiation::WXCreate(wxWindow* parent)
{
   //:TODO: Check mpWXCrystObj==0
   mpWXCrystObj=new WXRadiation(parent,this);
   return mpWXCrystObj;
}
#endif

////////////////////////////////////////////////////////////////////////
//
//    ScatteringData
//
////////////////////////////////////////////////////////////////////////

ScatteringData::ScatteringData():
mNbRefl(0),
mpCrystal(0),mUseFastLessPreciseFunc(false),
mpTemperatureFactor(0),mpScatteringFactor(0),mpRealGeomSF(0),mpImagGeomSF(0),
mpScattCompList(0),mNbScatteringPower(0),
mIgnoreImagScattFact(false)
{
   VFN_DEBUG_MESSAGE("ScatteringData::ScatteringData()",10)
   
   //:TODO: To be removed when using Clocks
   mAnomalousNeedRecalc=true;
   mThermicNeedRecalc=true;
   mScattFactNeedRecalc=true;
   mGeomFhklCalcNeedRecalc=true;
   mFhklCalcNeedRecalc=true;
}

ScatteringData::ScatteringData(const ScatteringData &old):
mNbRefl(old.mNbRefl),mRadiation(old.mRadiation),
mpCrystal(old.mpCrystal),mUseFastLessPreciseFunc(old.mUseFastLessPreciseFunc),
//Do not copy temporary arrays
mpTemperatureFactor(0),mpScatteringFactor(0),mpRealGeomSF(0),mpImagGeomSF(0),
mClockHKL(old.mClockHKL),
mpScattCompList(0),mNbScatteringPower(0),

mIgnoreImagScattFact(old.mIgnoreImagScattFact)
{
   VFN_DEBUG_MESSAGE("ScatteringData::ScatteringData(&old)",10)
   mClockStructFactor.Reset();
   mClockTheta.Reset();
   mClockScattFactor.Reset();
   mClockScattFactorResonant.Reset();
   mClockThermicFact.Reset();
   this->SetHKL(old.GetH(),old.GetK(),old.GetL());
   VFN_DEBUG_MESSAGE("ScatteringData::ScatteringData(&old):End",5)
   
   //:TODO: To be removed when using Clocks
   mAnomalousNeedRecalc=true;
   mThermicNeedRecalc=true;
   mScattFactNeedRecalc=true;
   mGeomFhklCalcNeedRecalc=true;
   mFhklCalcNeedRecalc=true;
}

ScatteringData::~ScatteringData()
{
   VFN_DEBUG_MESSAGE("ScatteringData::~ScatteringData()",10)
   if(0 != mpTemperatureFactor) delete[] mpTemperatureFactor;
   if(0 != mpScatteringFactor) delete[] mpScatteringFactor;
   if(0 != mpRealGeomSF) delete[] mpRealGeomSF;
   if(0 != mpImagGeomSF) delete[] mpImagGeomSF;
   //if(0 != mpCrystal) mpCrystal->DeRegisterClient(*this);
}

void ScatteringData::SetHKL( CrystVector_double const &h,
                     CrystVector_double const &k,
                     CrystVector_double const &l)
{
   VFN_DEBUG_ENTRY("ScatteringData::SetHKL(h,k,l)",5)
   mNbRefl=h.numElements();
   mH=h;
   mK=k;
   mL=l;
   mClockHKL.Click();
   this->PrepareHKLarrays();
   VFN_DEBUG_EXIT("ScatteringData::SetHKL(h,k,l):End",5)
}

void ScatteringData::GenHKLFullSpace(const double maxTheta,const bool useMultiplicity)
{
   VFN_DEBUG_ENTRY("ScatteringData::GenHKLFullSpace()",5)
   TAU_PROFILE("ScatteringData::GenHKLFullSpace()","void (double,bool)",TAU_DEFAULT);
   if(mRadiation.GetWavelength()(0) <=.01)
   {
      throw ObjCrystException("ScatteringData::GenHKLFullSpace() \
      no wavelength assigned yet to this ScatteringData object.");;
   }
   if(0==mpCrystal)
   {
      throw ObjCrystException("ScatteringData::GenHKLFullSpace() \
      no crystal assigned yet to this ScatteringData object.");;
   }
   VFN_DEBUG_MESSAGE("ScatteringData::GenHKLFullSpace():Max theta="<<maxTheta \
   << " Using Multiplicity : "<<useMultiplicity,3)
   long maxH,maxK,maxL;
   maxH=(int) (sin(maxTheta)/mRadiation.GetWavelength()(0) * mpCrystal->GetLatticePar(0)*2+1);
   maxK=(int) (sin(maxTheta)/mRadiation.GetWavelength()(0) * mpCrystal->GetLatticePar(1)*2+1);
   maxL=(int) (sin(maxTheta)/mRadiation.GetWavelength()(0) * mpCrystal->GetLatticePar(2)*2+1);
   VFN_DEBUG_MESSAGE("->maxH : " << maxH << "  maxK : " << maxK << "maxL : " << maxL,5)
   mNbRefl=(2*maxH+1)*(2*maxK+1)*(2*maxL+1);
   CrystVector_long H(mNbRefl);
   CrystVector_long K(mNbRefl);
   CrystVector_long L(mNbRefl);
   long i=0;
   for(int h = -maxH ; h <= maxH;h++)
      for(int k = -maxK ; k <= maxK;k++)
         for(int l = -maxL ; l<= maxL;l++)
         {
            H(i)=h;
            K(i)=k;
            L(i)=l;
            i++;
         }
   VFN_DEBUG_MESSAGE("ScatteringData::GenHKLFullSpace():Finished setting h, k and l...",3)
   this->SetHKL(H,K,L);
   //this->CalcSinThetaLambda();//calc theta
   this->SortReflectionByTheta(maxTheta);
   if(true==useMultiplicity)
   {
      VFN_DEBUG_MESSAGE("ScatteringData::GenHKLFullSpace():Multiplicity...",3)
      {
         //CrystVector_double tmp;
         //tmp=mTheta;
         //cout << FormatVertVectorHKLFloats<double>(mH,mK,mL,tmp,12,4);
      }
      //generate 20 random atom positions, to check which reflections are equivalent
         int nbTestPositions=20;
         {//Init pseudo-random number generator
            time_t junk;
            time(&junk);
            tm *tmp=localtime(&junk);
            srand((unsigned)( (*tmp).tm_sec+60* (*tmp).tm_min+3600* (*tmp).tm_hour));
         }
         ScatteringComponentList scattList(nbTestPositions);
         CrystVector_long structFactorIndex(nbTestPositions);
         for(int i=0;i<nbTestPositions;i++)
         {
            scattList(i).mX=rand()/(double)RAND_MAX;
            scattList(i).mY=rand()/(double)RAND_MAX;
            scattList(i).mZ=rand()/(double)RAND_MAX;
            scattList(i).mOccupancy=1.;
            scattList(i).mDynPopCorr=1.;
            structFactorIndex(i)=i;
         }
      //calc geometrical struct factor for these positions
         CrystVector_double* realGeomSF=new CrystVector_double[nbTestPositions];
         CrystVector_double* imagGeomSF=new CrystVector_double[nbTestPositions];
         for(int i=0;i<nbTestPositions;i++)
         {
            this->CalcGeomStructFactor(scattList,mpCrystal->GetSpaceGroup(),structFactorIndex,
                  realGeomSF,imagGeomSF);
         }
      //OK, now sort reflections to keep or remove
         long nbKeptRefl=0;
         CrystVector_long subscriptKeptRefl(mNbRefl);
         mMultiplicity.resize(mNbRefl);
         CrystVector_bool treatedRefl(mNbRefl);
         long currentBaseRefl=0,testedRefl=0;
         double currentTheta=0;
         double compare;
         double h,k,l,h1,k1,l1;
         subscriptKeptRefl=0;
         mMultiplicity=0;
         treatedRefl=false;
      VFN_DEBUG_MESSAGE("ScatteringData::GenHKLFullSpace():Multiplicity 1",2)
         do
         {
      		VFN_DEBUG_MESSAGE("...Multiplicity 2",1)
            if(true==treatedRefl(currentBaseRefl)) continue;
            subscriptKeptRefl(nbKeptRefl)=currentBaseRefl;
            mMultiplicity(nbKeptRefl)=1;
            currentTheta=mTheta(currentBaseRefl);
            treatedRefl(currentBaseRefl)=true;
            h=mH(currentBaseRefl)+.001;
            k=mK(currentBaseRefl)+.001;
            l=mL(currentBaseRefl)+.001;
            testedRefl=currentBaseRefl+1;
            if(testedRefl==mNbRefl) break;
            bool test;
            do
            {
      			VFN_DEBUG_MESSAGE("...Multiplicity 3",1)
               compare=0;
               if(true==mIgnoreImagScattFact) //Friedel pairs are equivalent.
                  for(int i=0;i<nbTestPositions;i++) 
                     compare+= 
                        fabs(   fabs( (*(realGeomSF+i))(testedRefl     ) )
                              -fabs( (*(realGeomSF+i))(currentBaseRefl) ))
                       +fabs(   fabs( (*(imagGeomSF+i))(testedRefl     ) )
                              -fabs( (*(imagGeomSF+i))(currentBaseRefl) ));
               else
                  for(int i=0;i<nbTestPositions;i++) 
                     compare+= fabs(  (*(realGeomSF+i))(testedRefl)
                                    -(*(realGeomSF+i))(currentBaseRefl))
                              +fabs(  (*(imagGeomSF+i))(testedRefl)
                                    -(*(imagGeomSF+i))(currentBaseRefl));
      			VFN_DEBUG_MESSAGE("...Multiplicity 4",1)
               if(.001 > compare)
               {
                  mMultiplicity(nbKeptRefl) +=1;
                  treatedRefl(testedRefl)=true;
                  
                  //keep the reflection with 0) max indices positive then 
                  //1)max H, 2)max K and 3) max L
                  h1=mH(testedRefl)+.001;
                  k1=mK(testedRefl)+.001;
                  l1=mL(testedRefl)+.001;
      				VFN_DEBUG_MESSAGE("...Multiplicity 5",1)
      				VFN_DEBUG_MESSAGE(h1<<","<<k1<<","<<l1<<",",1)
      				VFN_DEBUG_MESSAGE(fabs(h1)<<","<<fabs(k1)<<","<<fabs(l1)<<",",1)
                  if( ((int)(h1/fabs(h1)+k1/fabs(k1)+l1/fabs(l1)))
                        > ((int)(h/fabs(h)+k/fabs(k)+l/fabs(l))) )
                  {
      					VFN_DEBUG_MESSAGE("...Multiplicity 6a",1)
                     subscriptKeptRefl(nbKeptRefl)=testedRefl;
                     h=h1;
                     k=k1;
                     l=l1;
                  } else
                  {
      					VFN_DEBUG_MESSAGE("...Multiplicity 6b",1)
                  	if( (int)(h1/fabs(h1)+k1/fabs(k1)+l1/fabs(l1))
                                 == (int)(h/fabs(h)+k/fabs(k)+l/fabs(l)) )
                  	{
                     	if(  (mH(testedRefl) > mH(subscriptKeptRefl(nbKeptRefl)))  ||

                         	((mH(testedRefl) == mH(subscriptKeptRefl(nbKeptRefl))) &&
                          	(mK(testedRefl) > mK(subscriptKeptRefl(nbKeptRefl)))) ||

                         	((mH(testedRefl) == mH(subscriptKeptRefl(nbKeptRefl))) &&
                          	(mK(testedRefl) == mK(subscriptKeptRefl(nbKeptRefl))) &&
                          	(mL(testedRefl) > mL(subscriptKeptRefl(nbKeptRefl)))) )
                     	{
                        	subscriptKeptRefl(nbKeptRefl)=testedRefl;
                        	h=h1;
                        	k=k1;
                        	l=l1;
                     	}
                     }
                  }
                  //cout << currentTheta*RAD2DEG << "  " <<
                  //       mIntH(subscriptKeptRefl(nbKeptRefl))<<"  "<<
                  //       mIntK(subscriptKeptRefl(nbKeptRefl))<<"  "<<
                  //       mIntL(subscriptKeptRefl(nbKeptRefl))<<"  ";
                  VFN_DEBUG_MESSAGE(mIntH(testedRefl)<<"  "<< mIntK(testedRefl)<<"  "<<mIntL(testedRefl),1);
                  //cout << "   " << compare << "  " << h/fabs(h)+k/fabs(k)+l/fabs(l) << endl;
               }
      			VFN_DEBUG_MESSAGE("...Multiplicity 5",1)
               testedRefl++;
               if(testedRefl<mNbRefl)
               {
                  if(fabs(currentTheta-mTheta(testedRefl)) < .002) test=true;
                  else test=false;
               }
               else test=false;
            } while(test);
            nbKeptRefl++;
         } while( ++currentBaseRefl < mNbRefl);
      VFN_DEBUG_MESSAGE("ScatteringData::GenHKLFullSpace():Multiplicity 2",2)
      //Keep only the elected reflections
         mNbRefl=nbKeptRefl;
         {
            CrystVector_double oldH,oldK,oldL;
            CrystVector_double oldWeight;
            long subs;
            
            oldH=mH;
            oldK=mK;
            oldL=mL;
            
            mMultiplicity.resizeAndPreserve(mNbRefl);
            subscriptKeptRefl.resizeAndPreserve(mNbRefl);
            mH.resize(mNbRefl);
            mK.resize(mNbRefl);
            mL.resize(mNbRefl);
            for(long i=0;i<mNbRefl;i++)
            {
               subs=subscriptKeptRefl(i);
               mH(i)=oldH(subs);
               mK(i)=oldK(subs);
               mL(i)=oldL(subs);
            }
         }
         this->PrepareHKLarrays();
         //this->CalcSinThetaLambda(true);
      delete[] realGeomSF;
      delete[] imagGeomSF;
   } //true==useMultiplicity
   else
   {
      mMultiplicity.resize(mNbRefl);
      mMultiplicity=1;
   }
   this->EliminateExtinctReflections();
   VFN_DEBUG_EXIT("ScatteringData::GenHKLFullSpace():End",5)
}

void ScatteringData::SetRadiationType(const RadiationType radiation)
{
   VFN_DEBUG_MESSAGE("ScatteringData::SetRadiationType():End",5)
   mRadiation.SetRadiationType(radiation);
}

RadiationType ScatteringData::GetRadiationType()const {return mRadiation.GetRadiationType();}

void ScatteringData::SetCrystal(Crystal &crystal)
{
   VFN_DEBUG_MESSAGE("ScatteringData::SetCrystal()",5)
   mpCrystal=&crystal;
   mSubObjRegistry.Register(crystal);
   //crystal.RegisterClient(*this);
   mClockGeomStructFact.Reset();
   mClockStructFactor.Reset();
}
const Crystal& ScatteringData::GetCrystal()const {return *mpCrystal;}
Crystal& ScatteringData::GetCrystal() {return *mpCrystal;}

long ScatteringData::GetNbRefl() const {return mNbRefl;}

const CrystVector_double& ScatteringData::GetH() const {return mH;}
const CrystVector_double& ScatteringData::GetK() const {return mK;}
const CrystVector_double& ScatteringData::GetL() const {return mL;}

const CrystVector_double& ScatteringData::GetH2Pi() const {return mH2Pi;}
const CrystVector_double& ScatteringData::GetK2Pi() const {return mK2Pi;}
const CrystVector_double& ScatteringData::GetL2Pi() const {return mH2Pi;}

const CrystVector_double& ScatteringData::GetSinThetaOverLambda()const
{
   VFN_DEBUG_ENTRY("ScatteringData::GetSinThetaOverLambda()",1)
   this->CalcSinThetaLambda();
   VFN_DEBUG_EXIT("ScatteringData::GetSinThetaOverLambda()",1)
   return mSinThetaLambda;
}
const CrystVector_double& ScatteringData::GetFhklCalcSq() const
{
   VFN_DEBUG_ENTRY("ScatteringData::GetFhklCalcSq()",2)
   this->CalcStructFactor();
   if(mClockStructFactorSq>mClockStructFactor) return mFhklCalcSq;
   #ifdef __LIBCRYST_VECTOR_USE_BLITZ__
   mFhklCalcSq=pow2(mFhklCalcReal)+pow2(mFhklCalcImag);
   #else
   const double *pr,*pi;
   double *p;
   pr=mFhklCalcReal.data();
   pi=mFhklCalcImag.data();
   p=mFhklCalcSq.data();
   for(long i=0;i<mNbRefl;i++)
   {
      *p++ = *pr * *pr + *pi * *pi;
      pr++;
      pi++;
   }
   #endif
   mClockStructFactorSq.Click();
   VFN_DEBUG_EXIT("ScatteringData::GetFhklCalcSq()",2)
   return mFhklCalcSq;
}
const CrystVector_double& ScatteringData::GetFhklCalcReal() const
{
   VFN_DEBUG_ENTRY("ScatteringData::GetFhklCalcReal()",2)
   this->CalcStructFactor();
   VFN_DEBUG_EXIT("ScatteringData::GetFhklCalcReal()",2)
   return mFhklCalcReal;
}

const CrystVector_double& ScatteringData::GetFhklCalcImag() const
{
   VFN_DEBUG_ENTRY("ScatteringData::GetFhklCalcImag()",2)
   this->CalcStructFactor();
   VFN_DEBUG_EXIT("ScatteringData::GetFhklCalcImag()",2)
   return mFhklCalcImag;
}

void ScatteringData::SetWavelength(const double lambda)
{
   VFN_DEBUG_MESSAGE("ScatteringData::SetWavelength() to "<<lambda,5)
   mRadiation.SetWavelength(lambda);
}

void ScatteringData::SetWavelength(const string &XRayTubeElementName,
                                   const double alpha2Alpha2ratio)
{
   VFN_DEBUG_MESSAGE("ScatteringData::SetWavelength() to "<<XRayTubeElementName,5)
   mRadiation.SetWavelength(XRayTubeElementName,alpha2Alpha2ratio);
}

void ScatteringData::SetEnergy(const double energy)
{
   this->SetWavelength(12398.4/energy);
}
CrystVector_double ScatteringData::GetWavelength()const {return mRadiation.GetWavelength();}

void ScatteringData::SetUseFastLessPreciseFunc(const bool useItOrNot)
{
   mUseFastLessPreciseFunc=useItOrNot;
   mClockGeomStructFact.Reset();
   mClockStructFactor.Reset();
}
void ScatteringData::SetIsIgnoringImagScattFact(const bool b)
{
   mIgnoreImagScattFact=b;
   mClockGeomStructFact.Reset();
   mClockStructFactor.Reset();
}
bool ScatteringData::IsIgnoringImagScattFact() const {return mIgnoreImagScattFact;}

void ScatteringData::PrintFhklCalc()const
{
   VFN_DEBUG_ENTRY("ScatteringData::PrintFhklCalc()",5)
   this->GetFhklCalcSq();
   CrystVector_double theta;
   theta=mTheta;
   theta *= RAD2DEG;
   cout <<" Number of reflections:"<<mNbRefl<<endl;
   cout <<"       H        K        L     F(hkl)^2     Re(F)         Im(F)";
   cout <<"        Theta       1/2d"<<endl;
   cout << FormatVertVectorHKLFloats<double>
               (mH,mK,mL,mFhklCalcSq,mFhklCalcReal,mFhklCalcImag,theta,mSinThetaLambda,12,4);
   VFN_DEBUG_EXIT("ScatteringData::PrintFhklCalc()",5)
}

void ScatteringData::PrepareHKLarrays()
{
   VFN_DEBUG_ENTRY("ScatteringData::PrepareHKLarrays()",5)
   mFhklCalcReal.resize(mNbRefl);
   mFhklCalcImag.resize(mNbRefl);
   mFhklCalcSq.resize(mNbRefl);
   
   mIntH=mH;
   mIntK=mK;
   mIntL=mL;
   
   mH2Pi=mH;
   mK2Pi=mK;
   mL2Pi=mL;
   mH2Pi*=(2*M_PI);
   mK2Pi*=(2*M_PI);
   mL2Pi*=(2*M_PI);
   mClockHKL.Click();
   VFN_DEBUG_EXIT("ScatteringData::PrepareHKLarrays()",5)
}

CrystVector_long ScatteringData::SortReflectionByTheta(const double maxTheta)
{
   TAU_PROFILE("ScatteringData::SortReflectionByTheta()","void ()",TAU_DEFAULT);
   VFN_DEBUG_ENTRY("ScatteringData::SortReflectionByTheta()",5)
   this->CalcSinThetaLambda();
   CrystVector_long sortedSubs;
   sortedSubs=SortSubs(mSinThetaLambda);
   CrystVector_long oldH,oldK,oldL;
   oldH=mH;
   oldK=mK;
   oldL=mL;
   long subs;
   long shift=0;
   
   //get rid of [0,0,0] reflection
   VFN_DEBUG_MESSAGE("ScatteringData::SortReflectionByTheta() 1",2)
   if(0==mSinThetaLambda(sortedSubs(0)))
   {
      shift=1;
      mNbRefl -= 1;
      mH.resize(mNbRefl);
      mK.resize(mNbRefl);
      mL.resize(mNbRefl);
   }
   VFN_DEBUG_MESSAGE("ScatteringData::SortReflectionByTheta() 2",2)
   for(long i=0;i<mNbRefl;i++)
   {
      subs=sortedSubs(i+shift);
      mH(i)=oldH(subs);
      mK(i)=oldK(subs);
      mL(i)=oldL(subs);
   }
   VFN_DEBUG_MESSAGE("ScatteringData::SortReflectionByTheta() 3",2)
   this->PrepareHKLarrays();
   this->CalcSinThetaLambda();
   
   VFN_DEBUG_MESSAGE("ScatteringData::SortReflectionByTheta() 4",2)
   if(0<maxTheta)
   {
      double maxsithsl=sin(maxTheta)/mRadiation.GetWavelength()(0);
      long maxSubs;
      for(maxSubs=0;(mSinThetaLambda(maxSubs)<maxsithsl) && (maxSubs<mNbRefl) ;maxSubs++);
      if(maxSubs==mNbRefl) return sortedSubs;
      mNbRefl=maxSubs;
      mH.resizeAndPreserve(mNbRefl);
      mK.resizeAndPreserve(mNbRefl);
      mL.resizeAndPreserve(mNbRefl);
      sortedSubs.resizeAndPreserve(mNbRefl);
      this->PrepareHKLarrays();
   }
   VFN_DEBUG_EXIT("ScatteringData::SortReflectionByTheta()",5)
   return sortedSubs;
}

CrystVector_long ScatteringData::EliminateExtinctReflections()
{
   TAU_PROFILE("ScatteringData::EliminateExtinctReflections()","void ()",TAU_DEFAULT);
   VFN_DEBUG_ENTRY("ScatteringData::EliminateExtinctReflections()",7)
   //return;
   //generate 20 random atom positions, to check which reflections are equivalent
      int nbTestPositions=20;
      {//Init pseudo-random number generator
         time_t junk;
         time(&junk);
         tm *tmp=localtime(&junk);
         srand((unsigned)( (*tmp).tm_sec+60* (*tmp).tm_min+3600* (*tmp).tm_hour));
      }
      ScatteringComponentList scattList(nbTestPositions);
      CrystVector_long structFactorIndex(nbTestPositions);
      for(int i=0;i<nbTestPositions;i++)
      {
         scattList(i).mX=rand()/(double)RAND_MAX;
         scattList(i).mY=rand()/(double)RAND_MAX;
         scattList(i).mZ=rand()/(double)RAND_MAX;
         scattList(i).mOccupancy=1.;
         scattList(i).mDynPopCorr=1.;
         structFactorIndex(i)=i;
      }
   VFN_DEBUG_MESSAGE("ScatteringData::EliminateExtinctReflections():2",5)
   //calc geometrical struct factor for these positions
      CrystVector_double* realGeomSF=new CrystVector_double[nbTestPositions];
      CrystVector_double* imagGeomSF=new CrystVector_double[nbTestPositions];
      for(int i=0;i<nbTestPositions;i++)
      {
         this->CalcGeomStructFactor(scattList,mpCrystal->GetSpaceGroup(),structFactorIndex,
               realGeomSF,imagGeomSF);
      }
   VFN_DEBUG_MESSAGE("ScatteringData::EliminateExtinctReflections():3",5)
   //OK, get reflections to keep
      long nbKeptRefl=0;
      CrystVector_long subscriptKeptRefl(mNbRefl);
      double value;
      subscriptKeptRefl=0;
      for(long j=0;j<mNbRefl;j++)
      {
         value=0;
         for(int i=0;i<nbTestPositions;i++) 
            value+= fabs( (*(realGeomSF+i))(j) )+fabs( (*(imagGeomSF+i))(j) );
         if(.01 < value) subscriptKeptRefl(nbKeptRefl++)=j;
      }
   VFN_DEBUG_MESSAGE("ScatteringData::EliminateExtinctReflections():4",5)
   //Keep only the elected reflections
      mNbRefl=nbKeptRefl;
      {
         CrystVector_long oldH,oldK,oldL;
         CrystVector_int oldMulti;
         long subs;

         oldH=mH;
         oldK=mK;
         oldL=mL;
         oldMulti=mMultiplicity;
         
         mMultiplicity.resize(mNbRefl);
         mH.resize(mNbRefl);
         mK.resize(mNbRefl);
         mL.resize(mNbRefl);
         for(long i=0;i<mNbRefl;i++)
         {
            subs=subscriptKeptRefl(i);
            mH(i)=oldH(subs);
            mK(i)=oldK(subs);
            mL(i)=oldL(subs);
            mMultiplicity(i)=oldMulti(subs);
         }
      }
   delete[] realGeomSF;
   delete[] imagGeomSF;
   this->PrepareHKLarrays();
   VFN_DEBUG_EXIT("ScatteringData::EliminateExtinctReflections():End",7)
   return subscriptKeptRefl;
}
void ScatteringData::PrepareCalcStructFactor()const
{
   //:TODO: Optimize using Clocks and not flags !!!
   // This is REALLY UGLY KLUDGE SHREK whatever... 

   VFN_DEBUG_ENTRY("ScatteringData::PrepareCalcStructFactor()"<<this->GetName(),4)
   TAU_PROFILE("ScatteringData::PrepareCalcStructFactor()","void (bool)",TAU_DEFAULT);
   mpScattCompList = &(mpCrystal->GetScatteringComponentList());
   const long nbComponent=mpScattCompList->GetNbComponent();
   
   //Compare this with the old list to know if any scattering component has moved
   VFN_DEBUG_MESSAGE("ScatteringData::PrepareCalcStructFactor():Compare with last list",2)
      if(nbComponent != mLastScattCompList.GetNbComponent())
      {
         VFN_DEBUG_MESSAGE("ScatteringData::PrepareCalcStructFactor():Resize GeomSF arrays",2)
         mGeomFhklCalcNeedRecalc=true;//mFhklCalcNeedRecalc is automatically set to true
      }
      else
      {
         for(long i=0;i<nbComponent;i++)
         {
            if(   ((*mpScattCompList)(i).mX != mLastScattCompList(i).mX)
                ||((*mpScattCompList)(i).mY != mLastScattCompList(i).mY)
                ||((*mpScattCompList)(i).mZ != mLastScattCompList(i).mZ))
            {
               mGeomFhklCalcNeedRecalc=true;
               break;
            }
            if(   ((*mpScattCompList)(i).mOccupancy != mLastScattCompList(i).mOccupancy)
                ||((*mpScattCompList)(i).mpScattPow != mLastScattCompList(i).mpScattPow))
                mFhklCalcNeedRecalc=true;
         }
      }
      
   //Build an index of ScatteringPower as they are stored in ScatteringData
   VFN_DEBUG_MESSAGE("ScatteringData::PrepareCalcStructFactor():Prepare index",2)
      mNbScatteringPower=0;
      CrystVector_long lastIndex;
      lastIndex=mScatteringPowerIndex;
      
      //this is the greatest index for the Scattering powers used for this calculation
      long nbScatteringPower=0;
      long tmp;
      for(long i=0;i<nbComponent;i++)
      {  
         tmp=(*mpScattCompList)(i).mpScattPow->GetScatteringPowerId()+1;
         //cout << "YOYO"<<tmp<<endl;
         if(tmp>nbScatteringPower) nbScatteringPower=tmp;
      }
      
   VFN_DEBUG_MESSAGE("ScatteringData::PrepareCalcStructFactor():which Scatt Pow are used?",2)
      
      // See which are used
         CrystVector_bool scatteringPowerIndexIsUsed(nbScatteringPower);
         scatteringPowerIndexIsUsed=false;
         for(long i=0;i<nbComponent;i++)
            scatteringPowerIndexIsUsed((*mpScattCompList)(i).mpScattPow->GetScatteringPowerId())
               =true;
      //Build the index CECI EST A REECRIRE UN BORDEL IMMONDE & ILLISIBLE....
         mScatteringPowerIndex.resize(nbScatteringPower);
         mScatteringPowerIndex=0;
      VFN_DEBUG_MESSAGE("ScatteringData::PrepareCalcStructFactor():Build index",2)
         for(long i=0;i<nbScatteringPower;i++)
            if(true==scatteringPowerIndexIsUsed(i))
               mScatteringPowerIndex(i)=mNbScatteringPower++;
         mScatteringPowerIndex2.resize(mNbScatteringPower);
         tmp=0;
         for(long i=0;i<nbScatteringPower;i++)
            if(true==scatteringPowerIndexIsUsed(i))
               mScatteringPowerIndex2(tmp++)=i;
      
      VFN_DEBUG_MESSAGE("->mScatteringPowerIndex:"<<endl<<mScatteringPowerIndex,1)
   //Compare this index with the last...
   VFN_DEBUG_MESSAGE("ScatteringData::PrepareCalcStructFactor():Compare with last index",2)
   if(mScatteringPowerIndex.numElements() != lastIndex.numElements())
   {//Recalc everything...
      VFN_DEBUG_MESSAGE("->Number of elements in index has changed!",2)
      mAnomalousNeedRecalc=true;
      mThermicNeedRecalc=true;
      mScattFactNeedRecalc=true;
      
      if(0!=mpRealGeomSF) delete[] mpRealGeomSF;
      if(0!=mpImagGeomSF) delete[] mpImagGeomSF;
      mpRealGeomSF=new CrystVector_double[mNbScatteringPower];
      mpImagGeomSF=new CrystVector_double[mNbScatteringPower];
   }
   else
   {
      if(MaxDifference(mScatteringPowerIndex,lastIndex) >0)
      {//Recalc everything...
         VFN_DEBUG_MESSAGE("->Contents of index has changed!",2)
         cout << "->Contents of index has changed!"<<endl;
         mAnomalousNeedRecalc=true;
         mThermicNeedRecalc=true;
         mScattFactNeedRecalc=true;
      }
      else
      {//check if any ScatteringPower has been modified since last time
         for(int i=0;i<nbScatteringPower;i++)
            if(true==scatteringPowerIndexIsUsed(i))
               if( GetScatteringPower(i).GetLastChangeClock()>mClockScattFactor)
               {
                  VFN_DEBUG_MESSAGE("->At least one scattering power has been modified:"<<GetScatteringPower(i).GetName(),2)
                  //mClockScattFactor.Print();
                  //GetScatteringPower(i).GetLastChangeClock().Print();
                  //cout << clock() <<endl;
                  mAnomalousNeedRecalc=true;
                  mThermicNeedRecalc=true;
                  mScattFactNeedRecalc=true;
                  break;
               }
      }
   }
   
   if(  (mClockScattFactor<mRadiation.GetClockWavelength()) || (mClockScattFactor<mClockHKL)
      ||(mClockScattFactor<mpCrystal->GetClockLatticePar())) 
   {
      VFN_DEBUG_MESSAGE("Wavelength, or HKLs, or lattice has changed",2)
      mScattFactNeedRecalc=true;
   }
   if(  (mClockThermicFact<mRadiation.GetClockWavelength()) || (mClockThermicFact<mClockHKL)
      ||(mClockThermicFact<mpCrystal->GetClockLatticePar())) 
   {
      VFN_DEBUG_MESSAGE("Wavelength, or HKLs, or lattice has changed",2)
      mThermicNeedRecalc=true;
   }
   if(  (mClockScattFactor        <mRadiation.GetClockRadiation())
      ||(mClockScattFactorResonant<mRadiation.GetClockRadiation())) 
   {
      VFN_DEBUG_MESSAGE("Radiation type has changed !",2)
      mAnomalousNeedRecalc=true;
      mScattFactNeedRecalc=true;
   }

   mLastScattCompList=*mpScattCompList;
   mClockScattFactor.Click();//update clock
   VFN_DEBUG_MESSAGE("->mGeomFhklCalcNeedRecalc:"<<mGeomFhklCalcNeedRecalc,2)
   VFN_DEBUG_MESSAGE("->mFhklCalcNeedRecalc    :"<<mFhklCalcNeedRecalc,2)
   VFN_DEBUG_MESSAGE("->mAnomalousNeedRecalc   :"<<mAnomalousNeedRecalc,2)
   VFN_DEBUG_MESSAGE("->mThermicNeedRecalc     :"<<mThermicNeedRecalc,2)
   VFN_DEBUG_MESSAGE("->mScattFactNeedRecalc   :"<<mScattFactNeedRecalc,2)
   VFN_DEBUG_EXIT("ScatteringData::PrepareCalcStructFactor():End",4)
}

void ScatteringData::CalcSinThetaLambda()const
{
   if( 0 == mpCrystal) throw ObjCrystException("ScatteringData::CalcSinThetaLambda() \
      Cannot compute sin(theta)/lambda : there is no crystal affected to this \
      ScatteringData object yet.");

   if( 0 == this->GetNbRefl()) throw ObjCrystException("ScatteringData::CalcSinThetaLambda() \
      Cannot compute sin(theta)/lambda : there are no reflections !");

   if(  (mClockTheta>mRadiation.GetClockWavelength()) && (mClockTheta>mClockHKL)
      &&(mClockTheta>mpCrystal->GetClockLatticePar())) return;
   
   VFN_DEBUG_ENTRY("ScatteringData::CalcSinThetaLambda()",3)
   TAU_PROFILE("ScatteringData::CalcSinThetaLambda()","void (bool)",TAU_DEFAULT);
   mSinThetaLambda.resize(mNbRefl);
   
   const CrystMatrix_double bMatrix= mpCrystal->GetBMatrix();
   CrystMatrix_double xyz(this->GetNbRefl(),3);
   for(int i=0;i<this->GetNbRefl();i++)
   {  //:TODO: faster,nicer (implement a bloody matrix product !)
      // bMatrix(row,column)
      xyz(i,0)=bMatrix(0,0)*mH(i)+bMatrix(0,1)*mK(i)+bMatrix(0,2)*mL(i);
      xyz(i,1)=bMatrix(1,0)*mH(i)+bMatrix(1,1)*mK(i)+bMatrix(1,2)*mL(i);
      xyz(i,2)=bMatrix(2,0)*mH(i)+bMatrix(2,1)*mK(i)+bMatrix(2,2)*mL(i);
   }
   //cout << bMatrix << endl << xyz<<endl;
   for(int i=0;i< (this->GetNbRefl());i++)
      mSinThetaLambda(i)=sqrt(pow(xyz(i,0),2)+pow(xyz(i,1),2)+pow(xyz(i,2),2))/2;
   
   if(mRadiation.GetWavelength()(0) > 0)
   {
      mTheta.resize(mNbRefl);
      mTanTheta.resize(mNbRefl);
      for(int i=0;i< (this->GetNbRefl());i++) 
      {  
      	if( (mSinThetaLambda(i)*mRadiation.GetWavelength()(0))>1)
      	{
      		//:KLUDGE: :TODO:
      		mTheta(i)=M_PI;
      		mTanTheta(i)=1e6;
      		/*
      		ofstream out("log.txt");
      		out << "Error when computing Sin(theta) :"
      			 << "i="<<i<<" ,mSinThetaLambda(i)="<<mSinThetaLambda(i)
      			 << " ,mRadiation.GetWavelength()(0)="
      			 << mRadiation.GetWavelength()(0) 
      			 << " ,H="<<mH(i)
      			 << " ,K="<<mK(i)
      			 << " ,L="<<mL(i)
      			 <<endl;
   			out.close();
      		abort();
      		*/
      	}
         else 
         {
         	mTheta(i)=asin(mSinThetaLambda(i)*mRadiation.GetWavelength()(0));
         	mTanTheta(i)=tan(mTheta(i));
         }
      }
   } else 
   {
      cout << "Wavelength not given in ScatteringData::CalcSinThetaLambda() !" <<endl;
      throw 0;
   }
      
   mClockTheta.Click();
   VFN_DEBUG_EXIT("ScatteringData::CalcSinThetaLambda()",3)
}

void ScatteringData::CalcScattFactor()const
{
   //check if the number of ScatteringPower has changed, and if
   bool force=false;
   if(true == mScattFactNeedRecalc) force=true;
   if(force == false) return;
   TAU_PROFILE("ScatteringData::CalcScattFactor()","void (bool)",TAU_DEFAULT);
   VFN_DEBUG_ENTRY("ScatteringData::CalcScattFactor()",4)
   this->CalcResonantScattFactor();
   
   if(mpScatteringFactor !=0) delete[] mpScatteringFactor;
   mpScatteringFactor=new CrystVector_double[mNbScatteringPower];
   this->CalcResonantScattFactor();
   for(int i=0;i<mNbScatteringPower;i++)
   {
      mpScatteringFactor[i]=
         GetScatteringPower(mScatteringPowerIndex2(i)).GetScatteringFactor(*this);
      //Directly add Fprime
      // if(mIgnoreImagScattFact==true) 
      mpScatteringFactor[i] += this->mFprime(i);
   }
   mScattFactNeedRecalc=false;
   mFhklCalcNeedRecalc=true;
   #if 0 //#ifdef __DEBUG__
   #endif
   #ifdef __DEBUG__
   if(gVFNDebugMessageLevel<=1)
   {
      cout <<"Scattering factors:"<<endl;
      for(int i=0;i<mNbScatteringPower;i++)
         GetScatteringPower(mScatteringPowerIndex2(i)).Print();
      cout << "#    Sin(theta/lambda) " ;
      for(int i=0;i<mNbScatteringPower;i++) 
         cout << "  " <<GetScatteringPower(mScatteringPowerIndex2(i)).GetName();
      cout << endl;
      cout << FormatVertVector<double>(mSinThetaLambda,mpScatteringFactor,mNbScatteringPower);
   }
   #endif
   mClockScattFactor.Click();
   VFN_DEBUG_EXIT("ScatteringData::CalcScattFactor()",4)
}

void ScatteringData::CalcTemperatureFactor()const
{
   bool force=false;
   //check if the number of atoms has changed, and if the thermic factors
   //are unchanged
   if(true == mThermicNeedRecalc) force=true;
   if(force == false) return;
   TAU_PROFILE("ScatteringData::CalcTemperatureFactor()","void (bool)",TAU_DEFAULT);
   VFN_DEBUG_ENTRY("ScatteringData::CalcTemperatureFactor()",4)

   if(mpTemperatureFactor !=0) delete[] mpTemperatureFactor;
   mpTemperatureFactor=new CrystVector_double[mNbScatteringPower];
   
   for(int i=0;i<mNbScatteringPower;i++)
      mpTemperatureFactor[i]=
         GetScatteringPower(mScatteringPowerIndex2(i)).GetTemperatureFactor(*this);
   mThermicNeedRecalc=false;
   mFhklCalcNeedRecalc=true;
   #ifdef __DEBUG__
   if(gVFNDebugMessageLevel<2)
   {
      cout <<"Thermic factors:"<<endl;
      cout << "#    Sin(theta/lambda) " ;
      for(int i=0;i<mNbScatteringPower;i++) 
         cout << "  " <<GetScatteringPower(mScatteringPowerIndex2(i)).GetName();
      cout << endl;
      cout << FormatVertVector<double>(mSinThetaLambda,mpTemperatureFactor,mNbScatteringPower);
   }
   #endif
   mClockThermicFact.Click();
   VFN_DEBUG_EXIT("ScatteringData::CalcTemperatureFactor()",4)
}

void ScatteringData::CalcResonantScattFactor()const
{
   bool force=false;
   //check if the number of atoms has changed, and if the atomic numbers
   //are unchanged (a change of wavelength automatically calls this function, 
   //so it is not checked)
   //:TODO: do not assume monochromatic
   if(true == mAnomalousNeedRecalc) force=true;
   if(force == false) return;
   VFN_DEBUG_ENTRY("ScatteringData::GetResonantScattFactor()",4)
   TAU_PROFILE("ScatteringData::GetResonantScattFactor()","void (bool)",TAU_DEFAULT);
   mFprime.resize(mNbScatteringPower);
   mFsecond.resize(mNbScatteringPower);
   if(mRadiation.GetWavelength()(0) == 0)
   {
      mFprime=0;
      mFsecond=0;
      VFN_DEBUG_MESSAGE("->Lambda=0. fprime=fsecond=0",1)
      return;
   }
   for(int i=0;i<mNbScatteringPower;i++)
   {//probably real slow, but hardly ever executed
      mFprime(i) =
         (GetScatteringPower(mScatteringPowerIndex2(i)).GetResonantScattFactReal(*this))(0);
      mFsecond(i)=
         (GetScatteringPower(mScatteringPowerIndex2(i)).GetResonantScattFactImag(*this))(0);
      VFN_DEBUG_MESSAGE("->(scatt#"<<i<<"(Lambda=0. fprime="<<mFprime(i)<<"fsecond="<<mFsecond(i),1)
   }
   mAnomalousNeedRecalc=false;
   mFhklCalcNeedRecalc=true;
   #ifdef __DEBUG__
   if(gVFNDebugMessageLevel<2)
   {
      cout <<"Resonant Scattering factors:"<<endl;
      cout <<"  f' : "<< mFprime <<endl;
      cout <<"  f' : "<< mFsecond <<endl;
   }
   #endif
   mClockScattFactorResonant.Click();
   VFN_DEBUG_EXIT("ScatteringData::GetResonantScattFactor()",4)
}

void ScatteringData::CalcStructFactor() const
{
   TAU_PROFILE("ScatteringData::CalcStructFactor()","void ()",TAU_DEFAULT);
   VFN_DEBUG_ENTRY("ScatteringData::CalcStructFactor()",3)
   //:TODO: Anisotropic Thermic factors
   TAU_PROFILE_TIMER(timer1,"ScatteringData::CalcStructFactor1:Prepare","", TAU_FIELD);
   TAU_PROFILE_TIMER(timer2,"ScatteringData::CalcStructFactor2:GeomStructFact","", TAU_FIELD);
   TAU_PROFILE_TIMER(timer3,"ScatteringData::CalcStructFactor3:Scatt.Factors","", TAU_FIELD);
   TAU_PROFILE_TIMER(timer4,"ScatteringData::CalcStructFactor4:Finish,DynCorr","", TAU_FIELD);

   TAU_PROFILE_START(timer1);
   this->PrepareCalcStructFactor();
   
   const long nbRefl=this->GetNbRefl();
   const long nbScattComp=mpScattCompList->GetNbComponent();
   this->CalcSinThetaLambda();
   
   //cout <<"CalcStructFactor:"<< endl << FormatVertVector<double>(mH,mK,mL)<<endl ;
   //get atom positions
   //cout << "Getting all atoms coordinates..."<<endl;
   TAU_PROFILE_STOP(timer1);
   //compute structure factor for each atom & reflection
   if(true==mGeomFhklCalcNeedRecalc)
   {
      TAU_PROFILE_START(timer2);
      VFN_DEBUG_MESSAGE("ScatteringData::CalcStructFactor():Getting Geom Struct Fact...",3)
      CrystVector_long structFactorIndex(nbScattComp);
      for(long i=0;i<nbScattComp;i++)
         structFactorIndex(i)=
            mScatteringPowerIndex((*mpScattCompList)(i).mpScattPow->GetScatteringPowerId());
            
      //mpRealGeomSF and mpImagGeomSF are allocated by PrepareCalcStructFactor()
      this->CalcGeomStructFactor(*mpScattCompList,mpCrystal->GetSpaceGroup(),structFactorIndex,
                                              mpRealGeomSF,mpImagGeomSF,
                                              mUseFastLessPreciseFunc);
      //cout << mpRealGeomSF[0]<<endl;
      mFhklCalcNeedRecalc=true;
      #ifdef __DEBUG__
      if(gVFNDebugMessageLevel<2)
      {
         cout <<"Geometrical Structure Factor, real:"<<endl;
         cout << "#    Sin(theta/lambda) " ;
         for(int i=0;i<mNbScatteringPower;i++) 
            cout << "  " <<(*mpScattCompList)(i).mpScattPow->GetName();
         cout << endl;
         cout << FormatVertVector<double>(mSinThetaLambda,mpRealGeomSF,mNbScatteringPower);
         cout <<"Geometrical Structure Factor, imag:"<<endl;
         cout << "#    Sin(theta/lambda) " ;
         for(int i=0;i<mNbScatteringPower;i++) 
            cout << "  " <<(*mpScattCompList)(i).mpScattPow->GetName();
         cout << endl;
         cout << FormatVertVector<double>(mSinThetaLambda,mpImagGeomSF,mNbScatteringPower);
      }
      #endif
      TAU_PROFILE_STOP(timer2);
   }
   // else cout << "ScatteringData::CalcStructFactor(): No need to recalc :-)))"
   

   
   //Recompute scattering factors ?
   TAU_PROFILE_START(timer3);
   this->CalcScattFactor();
   this->CalcResonantScattFactor();
   this->CalcTemperatureFactor();
   TAU_PROFILE_STOP(timer3);
   
   
   //OK, really must recompute SFs?
   if(true==mFhklCalcNeedRecalc)
   {
      TAU_PROFILE_START(timer4);
      VFN_DEBUG_MESSAGE("ScatteringData::CalcStructFactor():Fhkl Recalc...",3)
      //reset Fcalc
         mFhklCalcReal.resize(nbRefl);
         mFhklCalcImag.resize(nbRefl);
         mFhklCalcReal=0;
         mFhklCalcImag=0;
      //Add all contributions
      const double *pGeomR;
      const double *pGeomI;
      const double *pScatt;
      const double *pTemp;
      double *pReal;
      double *pImag;
      double fsecond;
      for(long i=0;i<mNbScatteringPower;i++)
      {
         VFN_DEBUG_MESSAGE("ScatteringData::CalcStructFactor():Fhkl Recalc, comp#"<<i,2)
         pGeomR=(mpRealGeomSF+i)->data();
         pGeomI=(mpImagGeomSF+i)->data();
         pScatt=(mpScatteringFactor+i)->data();
         pTemp=(mpTemperatureFactor+i)->data();
         
         pReal=mFhklCalcReal.data();
         pImag=mFhklCalcImag.data();
         
         VFN_DEBUG_MESSAGE("->(mpRealGeomSF+i) "<<(mpRealGeomSF+i)->numElements()<<"elements",2)
         VFN_DEBUG_MESSAGE("->(mpImagGeomSF+i) "<<(mpImagGeomSF+i)->numElements()<<"elements",2)
         VFN_DEBUG_MESSAGE("->(mpScatteringFactor+i)"<<(mpScatteringFactor+i)->numElements()
            <<"elements",1)
         VFN_DEBUG_MESSAGE("->(mpTemperatureFactor+i)" << (mpTemperatureFactor+i)->numElements()<<"elements",1)
         VFN_DEBUG_MESSAGE("->mFhklCalcReal "<<mFhklCalcReal.numElements()<<"elements",2)
         VFN_DEBUG_MESSAGE("->mFhklCalcImag "<<mFhklCalcImag.numElements()<<"elements",2)
         
         // //Now where are stored the scattering power values ?
         //const long index2=(*mpScattCompList)(i).mpScattPow->GetScatteringPowerId();
         // const long index=mScatteringPowerIndex(index2);

         //VFN_DEBUG_MESSAGE(FormatHorizVector<double>(*(mpScatteringFactor+index)),0);
         //VFN_DEBUG_MESSAGE("->Scattering Power:"<<GetScatteringPower(index2).GetName()<<i,10)
         //cout << i <<":";GetScatteringPower(index2).Print();
         //VFN_DEBUG_MESSAGE(FormatVertVector<double>(mH,mK,mL,
         //                                          *(mpRealGeomSF+i),
         //                                          *(mpImagGeomSF+i)),1;
         VFN_DEBUG_MESSAGE("->   H      K      L   sin(t/l)     Re(F)      Im(F)      scatt      Temp",1)

         VFN_DEBUG_MESSAGE(FormatVertVectorHKLFloats<double>(mH,mK,mL,mSinThetaLambda,
                                                            *(mpRealGeomSF+i),
                                                            *(mpImagGeomSF+i),
                                                            *(mpScatteringFactor+i),
                                                            *(mpTemperatureFactor+i) 
                                                            ),1);

         if(false==mIgnoreImagScattFact)
         {
            fsecond=mFsecond(i);
            VFN_DEBUG_MESSAGE("->fsecond= "<<fsecond,2)
            for(long j=0;j<nbRefl;j++)
            {
               *pReal += *pGeomR  * *pTemp * *pScatt - *pGeomI * *pTemp * fsecond;
               *pImag += *pGeomI  * *pTemp * *pScatt + *pGeomR * *pTemp * fsecond;
            VFN_DEBUG_MESSAGE("-->"<<j<<" "<<*pReal<<" "<<*pImag<<" "<<*pGeomR<<" "<<*pGeomI<<" "<<*pScatt<<" "<<*pTemp,1)
               pGeomR++;pGeomI++;pTemp++;pScatt++;pReal++;pImag++;
            }
         }
         else
         {
            for(long j=0;j<nbRefl;j++)
            {
               *pReal++ += *pGeomR  * *pTemp * *pScatt;
               *pImag++ += *pGeomI  * *pTemp * *pScatt;
               pGeomR++;pGeomI++;pTemp++;pScatt++;
            }
         }
      }
      TAU_PROFILE_STOP(timer4);
   }
   else 
   {
      VFN_DEBUG_EXIT("ScatteringData::CalcStructFactor():Fhkl NOT Recalc...",3)
      return;
   }
   

   mFhklCalcNeedRecalc=false;
   mGeomFhklCalcNeedRecalc=false;
   mClockStructFactor.Click();
   VFN_DEBUG_EXIT("ScatteringData::CalcStructFactor()",3)
}

void ScatteringData::CalcGeomStructFactor(const ScatteringComponentList &scattCompList,
                          const SpaceGroup &spg,
                          const CrystVector_long &structFactorIndex,
                          CrystVector_double* rsf2,
                          CrystVector_double* isf2,
                          bool useFastTabulatedTrigFunctions) const
{
   TAU_PROFILE("ScatteringData::GeomStructFactor()","void (Vx,Vy,Vz,data,M,M,bool)",TAU_DEFAULT);
   VFN_DEBUG_ENTRY("ScatteringData::GeomStructFactor(Vx,Vy,Vz,...)",3)
   VFN_DEBUG_MESSAGE("-->Using fast functions:"<<useFastTabulatedTrigFunctions,2)
   VFN_DEBUG_MESSAGE("-->Number of translation vectors:"<<spg.GetNbTranslationVectors()-1,2)
   VFN_DEBUG_MESSAGE("-->Has an inversion Center:"<<spg.HasInversionCenter(),2)
   VFN_DEBUG_MESSAGE("-->Number of symetry operations (w/o transl&inv cent.):"\
                     <<spg.GetNbSymmetrics(true,true),2)
   VFN_DEBUG_MESSAGE("-->Number of Scattering Components :"<<scattCompList.GetNbComponent(),2)
   VFN_DEBUG_MESSAGE("-->Number of reflections:"<<this->GetNbRefl(),2)
   #ifdef __DEBUG__
   static long counter=0;
   VFN_DEBUG_MESSAGE("-->Number of GeomStructFactor calculations so far:"<<counter++,3)
   #endif
   
   //:TODO: implement for geometrical structure factor calculation
   //bool useGeomStructFactor=mUseGeomStructFactor;
   
   //if((mfpRealGeomStructFactor==0)||(mfpImagGeomStructFactor==0)) useGeomStructFactor=false ;
   
   //if(useGeomStructFactor==true)
   //{
   //   (*mfpRealGeomStructFactor)(x,y,z,data.H2Pi(),data.K2Pi(),data.L2Pi(),rsf);
   //   if(this->IsCentrosymmetric())return;
   //   (*mfpImagGeomStructFactor)(x,y,z,data.H2Pi(),data.K2Pi(),data.L2Pi(),isf);
   //   return;
   //}
   //else
   {  
      const int nbSymmetrics=spg.GetNbSymmetrics(true,true);
      const int nbTranslationVectors=spg.GetNbTranslationVectors();
      const long nbComp=scattCompList.GetNbComponent();
      const int nbRefl=this->GetNbRefl();
                  
      CrystMatrix_double transVect(nbTranslationVectors,3);
      transVect=spg.GetTranslationVectors();
      CrystMatrix_double allCoords(nbSymmetrics,3);
      CrystVector_double tmpVect(nbRefl);
      
      if(useFastTabulatedTrigFunctions==true)
      {
         if(sLibCrystTabulCosineIsInit==false )
         {
            cout << "Init tabulated sin and cosine functions."<<endl;
            InitLibCrystTabulCosine();
            sLibCrystTabulCosineIsInit=true;
         }
      }
      CrystVector_long intVect(nbRefl);//not used if useFastTabulatedTrigFunctions==false
      
      //Get the number of scattering components (and therefore arrays) used
      long nbScattArrays=0;
      for(long i=0;i<nbComp;i++)
         if(structFactorIndex(i)>nbScattArrays) nbScattArrays=structFactorIndex(i);
      nbScattArrays++;
      
      //Resize all arrays and set them to 0
      for(long i=0;i<nbScattArrays;i++)
      {
         (rsf2+i)->resize(nbRefl);
         (isf2+i)->resize(nbRefl);
         *(rsf2+i)=0;
         *(isf2+i)=0;
      }
      
      double x,y,z;
      const register double * hh;
      const register double * kk;
      const register double * ll;
      
      for(long i=0;i<nbComp;i++)
      {
         x=scattCompList(i).mX;
         y=scattCompList(i).mY;
         z=scattCompList(i).mZ;
         
         const double popu= scattCompList(i).mOccupancy * scattCompList(i).mDynPopCorr;
         
         const long index=structFactorIndex(i);
         
         register long *tmpInt;
         long intX,intY,intZ;
         const long *intH,*intK,*intL;
         register double * tmp;
         
         allCoords=spg.GetAllSymmetrics(x,y,z,true,true);
         if((true==spg.HasInversionCenter()) && (false==spg.IsInversionCenterAtOrigin()))
         {
            for(int j=0;j<nbSymmetrics;j++)
            {
               //The phase of the structure factor will be wrong
               //This is fixed a bit further...
               allCoords(j,0) -= ((double)spg.GetSgOps().InvT[0])/STBF/2.;
               allCoords(j,1) -= ((double)spg.GetSgOps().InvT[1])/STBF/2.;
               allCoords(j,2) -= ((double)spg.GetSgOps().InvT[2])/STBF/2.;
            }
         }
         for(int j=0;j<nbSymmetrics;j++)
         {
            
            if(useFastTabulatedTrigFunctions==true)
            {
               double *rrsf,*iisf;
               rrsf=(rsf2+index)->data();
               iisf=(isf2+index)->data();
               
               intX=(long)(allCoords(j,0)*sLibCrystNbTabulSine);
               intY=(long)(allCoords(j,1)*sLibCrystNbTabulSine);
               intZ=(long)(allCoords(j,2)*sLibCrystNbTabulSine);

               intH=mIntH.data();
               intK=mIntK.data();
               intL=mIntL.data();

               tmpInt=intVect.data();

               //the +sLibCrystNbTabulSine*1000 ensures that the resulting integer is >0
               for(int jj=0;jj<nbRefl;jj++) 
                *tmpInt++ = (*intH++ * intX + *intK++ * intY + *intL++ *intZ 
                              +sLibCrystNbTabulSine*1000)% sLibCrystNbTabulSine ;

               //Doing 2 loops for sine and cosine is faster than just 1 (for not-centro)
               tmpInt=intVect.data();
               for(int jj=0;jj<nbRefl;jj++) 
                  *rrsf++ += popu * spLibCrystTabulCosine[*tmpInt++];
               tmpInt=intVect.data();

               if(false==spg.HasInversionCenter()) 
                  for(int jj=0;jj<nbRefl;jj++) 
                     *iisf++ += popu * spLibCrystTabulSine[*tmpInt++];
            }
            else
            {
               x=allCoords(j,0);
               y=allCoords(j,1);
               z=allCoords(j,2);
               hh=mH2Pi.data();
               kk=mK2Pi.data();
               ll=mL2Pi.data();
               tmp=tmpVect.data();
               
               
               for(int jj=0;jj<nbRefl;jj++) *tmp++ = *hh++ * x + *kk++ * y + *ll++ *z;
               
               double *sf=(rsf2+index)->data();
               tmp=tmpVect.data();
               
               for(int jj=0;jj<nbRefl;jj++) *sf++ += popu * cos(*tmp++);
               
               if(false==spg.HasInversionCenter()) 
               {
                  sf=(isf2+index)->data();
                  tmp=tmpVect.data();
                  for(int jj=0;jj<nbRefl;jj++) *sf++ += popu * sin(*tmp++);
               }
            }
         }
      }//for all components...
      
      if(nbTranslationVectors > 1)
      {
         tmpVect=1;
         if( (spg.GetSpaceGroupNumber()>= 143) && (spg.GetSpaceGroupNumber()<= 167))
         {//Special case for trigonal groups R3,...
            double *p1=tmpVect.data();
            hh=mH2Pi.data();
            kk=mK2Pi.data();
            ll=mL2Pi.data();
            for(long j=0;j<nbRefl;j++) *p1++ += 2*cos((*hh++ - *kk++ - *ll++)/3.);
         }
         else
         {
            for(int j=1;j<nbTranslationVectors;j++)
            {
               x=transVect(j,0);
               y=transVect(j,1);
               z=transVect(j,2);
               double *p1=tmpVect.data();
               hh=mH2Pi.data();
               kk=mK2Pi.data();
               ll=mL2Pi.data();
               for(long j=0;j<nbRefl;j++) *p1++ += cos(*hh++ *x + *kk++ *y + *ll++ *z );
            }
         }
         for(long i=0;i<nbScattArrays;i++) *(rsf2+i) *= tmpVect;
         
         if(false==spg.HasInversionCenter())
            for(long i=0;i<nbScattArrays;i++) *(isf2+i) *= tmpVect;
      }
      if(true==spg.HasInversionCenter())
      {
         for(long i=0;i<nbScattArrays;i++) *(rsf2+i) *=2.;
         if(false==spg.IsInversionCenterAtOrigin())
         {
            VFN_DEBUG_MESSAGE("ScatteringData::GeomStructFactor(Vx,Vy,Vz):\
               Inversion Center not at the origin...",2)
            //fix the phase of each reflection when the inversion center is not
            //at the origin, using :
            // Re(F) = RSF*cos(2pi(h*Xc+k*Yc+l*Zc))
            // Re(F) = RSF*sin(2pi(h*Xc+k*Yc+l*Zc))
            //cout << "Glop Glop"<<endl;
            {
               const double xc=((double)spg.GetSgOps().InvT[0])/STBF/2.;
               const double yc=((double)spg.GetSgOps().InvT[1])/STBF/2.;
               const double zc=((double)spg.GetSgOps().InvT[2])/STBF/2.;
               #ifdef __LIBCRYST_VECTOR_USE_BLITZ__
               tmpVect = mH2Pi() * xc + mK2PI() * yc + mL2PI() * zc;
               #else
               {
                  const double *hh=mH2Pi.data();
                  const double *kk=mK2Pi.data();
                  const double *ll=mL2Pi.data();
                  double *ttmpVect=tmpVect.data();
                  for(long ii=0;ii<nbRefl;ii++) 
                     *ttmpVect++ = *hh++ * xc + *kk++ * yc + *ll++ * zc;
               }
               #endif
            }
            CrystVector_double cosTmpVect;
            CrystVector_double sinTmpVect;
            cosTmpVect=cos(tmpVect);
            sinTmpVect=sin(tmpVect);
            
            for(long i=0;i<nbScattArrays;i++)
            {
               *(isf2+i) = *(rsf2+i);
               *(isf2+i) *= sinTmpVect;
            }
            
            for(long i=0;i<nbScattArrays;i++)
               *(rsf2+i) *= cosTmpVect;
         }
      }

   }
   //cout << FormatVertVector<double>(*rsf2,*isf2)<<endl;
   VFN_DEBUG_EXIT("ScatteringData::GeomStructFactor(Vx,Vy,Vz,...)",3)
}

}//namespace ObjCryst
