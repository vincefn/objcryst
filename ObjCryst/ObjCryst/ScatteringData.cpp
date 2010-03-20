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
*  source file for LibCryst++ ScatteringData class
*
*/

#include <cmath>

#include <typeinfo>

#include "cctbx/sgtbx/space_group.h"
#include "cctbx/miller/index_generator.h"
#include "cctbx/miller/sym_equiv.h"
#include "cctbx/eltbx/wavelengths.h"

#include "ObjCryst/ObjCryst/ScatteringData.h"
#include "ObjCryst/Quirks/VFNDebug.h"
#include "ObjCryst/Quirks/VFNStreamFormat.h"
#include "ObjCryst/Quirks/Chronometer.h"

#ifdef __WX__CRYST__
   #include "ObjCryst/wxCryst/wxPowderPattern.h"
   #include "ObjCryst/wxCryst/wxRadiation.h"
#endif

#include <fstream>
#include <iomanip>
#include <stdio.h> //for sprintf()

namespace ObjCryst
{
const RefParType *gpRefParTypeScattData= 0;
const RefParType *gpRefParTypeScattDataScale=0;
const RefParType *gpRefParTypeScattDataProfile=0;
const RefParType *gpRefParTypeScattDataProfileType=0;
const RefParType *gpRefParTypeScattDataProfileWidth=0;
const RefParType *gpRefParTypeScattDataProfileAsym=0;
const RefParType *gpRefParTypeScattDataCorr=0;
const RefParType *gpRefParTypeScattDataCorrInt=0;
const RefParType *gpRefParTypeScattDataCorrIntPO_Direction=0;
const RefParType *gpRefParTypeScattDataCorrIntPO_Fraction=0;
const RefParType *gpRefParTypeScattDataCorrIntPO_Amplitude=0;
const RefParType *gpRefParTypeScattDataCorrIntAbsorp=0;
const RefParType *gpRefParTypeScattDataCorrIntPolar=0;
const RefParType *gpRefParTypeScattDataCorrIntExtinc=0;
const RefParType *gpRefParTypeScattDataCorrPos=0;
const RefParType *gpRefParTypeScattDataBackground=0;

const RefParType *gpRefParTypeRadiation=0;
const RefParType *gpRefParTypeRadiationWavelength=0;

long NiftyStaticGlobalObjectsInitializer_ScatteringData::mCount=0;
//######################################################################
//    Tabulated math functions for faster (&less precise) F(hkl) calculation
//These function are defined and used in cristallo-spacegroup.cpp
//Currently tabulating sine and cosine only
//######################################################################

static bool sLibCrystTabulCosineIsInit=false;
//conversion value
static REAL sLibCrystTabulCosineRatio;
// Number of tabulated values of cosine between [0;2pi]
// 100 000 is far enough for a model search, yielding a maximum
// error less than .05%... 10000 should be enough, too, with (probably) a higher cache hit
#define sLibCrystNbTabulSine 8192
#define sLibCrystNbTabulSineMASK 8191
//storage of tabulated values of cosine, and a table with interlaced csoine/sine values
static REAL *spLibCrystTabulCosine;
static REAL *spLibCrystTabulCosineSine;

void InitLibCrystTabulCosine()
{
   VFN_DEBUG_MESSAGE("InitLibCrystTabulCosine()",10)
   spLibCrystTabulCosine=new REAL[sLibCrystNbTabulSine];
   spLibCrystTabulCosineSine=new REAL[sLibCrystNbTabulSine*2];
   REAL *tmp=spLibCrystTabulCosine;
   sLibCrystTabulCosineRatio=sLibCrystNbTabulSine/2./M_PI;
   for(REAL i=0;i<sLibCrystNbTabulSine;i++) *tmp++ = cos(i/sLibCrystTabulCosineRatio);
   tmp=spLibCrystTabulCosineSine;
   for(REAL i=0;i<sLibCrystNbTabulSine;i++)
   {
      *tmp++ = cos(i/sLibCrystTabulCosineRatio);
      *tmp++ = sin(i/sLibCrystTabulCosineRatio);
   }
}

void DeleteLibCrystTabulCosine()
{
   delete[] spLibCrystTabulCosine;
   delete[] spLibCrystTabulCosineSine;
}

// Same for exponential calculations (used for global temperature factors)
static bool sLibCrystTabulExpIsInit=false;
static const long sLibCrystNbTabulExp=10000;
static const REAL sLibCrystMinTabulExp=-5.;
static const REAL sLibCrystMaxTabulExp=10.;
static REAL *spLibCrystTabulExp;
void InitLibCrystTabulExp()
{
   VFN_DEBUG_MESSAGE("InitLibCrystTabulExp()",10)
   spLibCrystTabulExp=new REAL[sLibCrystNbTabulExp];
   REAL *tmp=spLibCrystTabulExp;
   for(REAL i=0;i<sLibCrystNbTabulExp;i++) 
      *tmp++ = exp(sLibCrystMinTabulExp+i*(sLibCrystMaxTabulExp-sLibCrystMinTabulExp)/sLibCrystNbTabulExp);
   sLibCrystTabulExpIsInit=true;
}

void DeleteLibCrystTabulExp() { delete[] spLibCrystTabulExp;}

//:KLUDGE: The allocated memory for cos and sin table is never freed...
// This should be done after the last ScatteringData object is deleted.

////////////////////////////////////////////////////////////////////////
//
//    Radiation
//
////////////////////////////////////////////////////////////////////////
Radiation::Radiation():
mWavelength(1),mXRayTubeName(""),mXRayTubeDeltaLambda(0.),
mXRayTubeAlpha2Alpha1Ratio(0.5),mLinearPolarRate(0)
{
   mWavelength=1;
   this->InitOptions();
   mRadiationType.SetChoice(RAD_XRAY);
   mWavelengthType.SetChoice(WAVELENGTH_MONOCHROMATIC);
   mClockMaster.AddChild(mClockWavelength);
   mClockMaster.AddChild(mClockRadiation);
}

Radiation::Radiation(const RadiationType rad,const REAL wavelength)
{
   this->InitOptions();
   mRadiationType.SetChoice(rad);
   mWavelengthType.SetChoice(WAVELENGTH_MONOCHROMATIC);
   mWavelength=wavelength;
   mXRayTubeName="";
   mXRayTubeDeltaLambda=0.;//useless here
   mXRayTubeAlpha2Alpha1Ratio=0.5;//useless here
   mLinearPolarRate=0.95;//assume it's synchrotron ?
   mClockMaster.AddChild(mClockWavelength);
   mClockMaster.AddChild(mClockRadiation);
}

Radiation::Radiation(const string &XRayTubeElementName,const REAL alpha2Alpha2ratio)
{
   this->InitOptions();
   this->SetWavelength(XRayTubeElementName,alpha2Alpha2ratio);
   mClockMaster.AddChild(mClockWavelength);
   mClockMaster.AddChild(mClockRadiation);
}

Radiation::Radiation(const Radiation &old):
mRadiationType(old.mRadiationType),
mWavelengthType(old.mWavelengthType),
mWavelength(old.mWavelength),
mXRayTubeName(old.mXRayTubeName),
mXRayTubeDeltaLambda(old.mXRayTubeDeltaLambda),
mXRayTubeAlpha2Alpha1Ratio(old.mXRayTubeAlpha2Alpha1Ratio),
mLinearPolarRate(old.mLinearPolarRate)
{
   mClockWavelength.Click();
   mClockMaster.AddChild(mClockWavelength);
   mClockMaster.AddChild(mClockRadiation);
}

Radiation::~Radiation()
{}

const string& Radiation::GetClassName() const
{
   const static string className="Radiation";
   return className;
}

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
   if(rad == RAD_NEUTRON) mLinearPolarRate=0;
   if(rad == RAD_ELECTRON) mLinearPolarRate=0;
}

void Radiation::SetWavelengthType(const WavelengthType &type)
{
   mWavelengthType.SetChoice((unsigned long) type);
   if(type==WAVELENGTH_TOF) this->SetRadiationType(RAD_NEUTRON);
   if(type==WAVELENGTH_ALPHA12) this->SetRadiationType(RAD_XRAY);
}

WavelengthType Radiation::GetWavelengthType()const
{return (WavelengthType) mWavelengthType.GetChoice();}

const CrystVector_REAL& Radiation::GetWavelength() const {return mWavelength;}
void Radiation::SetWavelength(const REAL l)
{
   mWavelength.resize(1);
   mWavelength=l;
   mClockWavelength.Click();
}
void Radiation::SetWavelength(const string &XRayTubeElementName,
                              const REAL alpha2Alpha2ratio)
{
   VFN_DEBUG_MESSAGE("Radiation::SetWavelength(tubeName,ratio):",5)
   mXRayTubeName=XRayTubeElementName;
   this->SetRadiationType(RAD_XRAY);
   mWavelength.resize(1);
   mLinearPolarRate=0;
   
   if(XRayTubeElementName.length() >=3) //:KLUDGE:
   {
      mWavelengthType.SetChoice(WAVELENGTH_MONOCHROMATIC);
      if(XRayTubeElementName=="CoA1")
      {
         mWavelength=1.78901;
      }
      else
      {
         try
         {
            cctbx::eltbx::wavelengths::characteristic ch(mXRayTubeName);
            if(!ch.is_valid())
            {
               cout << "WARNING: could not interpret X-Ray tube name:"<<XRayTubeElementName<<endl
                    << "         not modifying wavelength !"<<endl;
               return;
            }
            mWavelength=ch.as_angstrom();
         }
         catch(cctbx::error)
         {
            cout << "WARNING: could not interpret X-Ray tube name:"<<XRayTubeElementName<<endl
                 << "         not modifying wavelength !"<<endl;
         }
      }
   }
   else
   {
      mWavelengthType.SetChoice(WAVELENGTH_ALPHA12);
      mXRayTubeAlpha2Alpha1Ratio=alpha2Alpha2ratio;
      REAL lambda1,lambda2;
      if(XRayTubeElementName=="Co")
      {
         lambda1=1.78901;
         lambda2=1.79290;
      }
      else
      {
         try
         {
            cctbx::eltbx::wavelengths::characteristic ch(mXRayTubeName+"A1");
            if(!ch.is_valid())
            {
               cout << "WARNING: could not interpret X-Ray tube name:"<<XRayTubeElementName<<endl
                    << "         not modifying wavelength !"<<endl;
               return;
            }
            lambda1=ch.as_angstrom();
            cctbx::eltbx::wavelengths::characteristic ch2(mXRayTubeName+"A2");
            if(!ch2.is_valid())
            {
               cout << "WARNING: could not interpret X-Ray tube name:"<<XRayTubeElementName<<endl
                    << "         not modifying wavelength !"<<endl;
               return;
            }
            lambda2=ch2.as_angstrom();
         }
         catch(cctbx::error)
         {
            cout << "WARNING: could not interpret X-Ray tube name:"<<XRayTubeElementName<<endl
                 << "         not modifying wavelength !"<<endl;
         }
      }
      mXRayTubeDeltaLambda=lambda2-lambda1;
      mWavelength=lambda1
            +mXRayTubeDeltaLambda*mXRayTubeAlpha2Alpha1Ratio/(1.+mXRayTubeAlpha2Alpha1Ratio);
   }
   mClockWavelength.Click();
}

REAL Radiation::GetXRayTubeDeltaLambda()const {return mXRayTubeDeltaLambda;}

REAL Radiation::GetXRayTubeAlpha2Alpha1Ratio()const {return mXRayTubeAlpha2Alpha1Ratio;}
      

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
      case WAVELENGTH_TOF: cout<< "Time Of Flight"<<" "<<endl;break;
   }
}

REAL Radiation::GetLinearPolarRate()const{return mLinearPolarRate;}

void Radiation::SetLinearPolarRate(const REAL f){mLinearPolarRate=f;}

void Radiation::InitOptions()
{
   static string RadiationTypeName;
   static string RadiationTypeChoices[3];
   static string WavelengthTypeName;
   static string WavelengthTypeChoices[3];
   
   static bool needInitNames=true;
   if(true==needInitNames)
   {
      RadiationTypeName="Radiation";
      RadiationTypeChoices[0]="Neutron";
      RadiationTypeChoices[1]="X-Ray";
      RadiationTypeChoices[2]="Electron";
      
      WavelengthTypeName="Spectrum";
      WavelengthTypeChoices[0]="Monochromatic";
      WavelengthTypeChoices[1]="X-Ray Tube";
      WavelengthTypeChoices[2]="Time Of Flight";
      //WavelengthTypeChoices[2]="MAD";
      //WavelengthTypeChoices[3]="DAFS";
      //WavelengthTypeChoices[4]="LAUE";
      
      needInitNames=false;//Only once for the class
   }
   mRadiationType.Init(3,&RadiationTypeName,RadiationTypeChoices);
   mWavelengthType.Init(3,&WavelengthTypeName,WavelengthTypeChoices);
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
mpCrystal(0),mGlobalBiso(0),mUseFastLessPreciseFunc(false),
mIgnoreImagScattFact(false),mMaxSinThetaOvLambda(10)
{
   VFN_DEBUG_MESSAGE("ScatteringData::ScatteringData()",10)
   {//This should be done elsewhere...
      RefinablePar tmp("Global Biso",&mGlobalBiso,-1.,1.,
                        gpRefParTypeScattPowTemperatureIso,REFPAR_DERIV_STEP_ABSOLUTE,
                        true,true,true,false,1.0);
      tmp.SetDerivStep(1e-4);
      tmp.AssignClock(mClockGlobalBiso);
      this->AddPar(tmp);
   }
   mClockMaster.AddChild(mClockHKL);
   mClockMaster.AddChild(mClockGlobalBiso);
   mClockMaster.AddChild(mClockNbReflUsed);
   mClockMaster.AddChild(mClockFhklObsSq);
}

ScatteringData::ScatteringData(const ScatteringData &old):
mNbRefl(old.mNbRefl),
mpCrystal(old.mpCrystal),mUseFastLessPreciseFunc(old.mUseFastLessPreciseFunc),
//Do not copy temporary arrays
mClockHKL(old.mClockHKL),
mIgnoreImagScattFact(old.mIgnoreImagScattFact),
mMaxSinThetaOvLambda(old.mMaxSinThetaOvLambda)
{
   VFN_DEBUG_MESSAGE("ScatteringData::ScatteringData(&old)",10)
   mClockStructFactor.Reset();
   mClockTheta.Reset();
   mClockScattFactor.Reset();
   mClockScattFactorResonant.Reset();
   mClockThermicFact.Reset();
   this->SetHKL(old.GetH(),old.GetK(),old.GetL());
   VFN_DEBUG_MESSAGE("ScatteringData::ScatteringData(&old):End",5)
   {//This should be done elsewhere...
      RefinablePar tmp("Global Biso",&mGlobalBiso,-1.,1.,
                        gpRefParTypeScattPowTemperatureIso,REFPAR_DERIV_STEP_ABSOLUTE,
                        true,true,true,false,1.0);
      tmp.SetDerivStep(1e-4);
      tmp.AssignClock(mClockGlobalBiso);
      this->AddPar(tmp);
   }
   mClockMaster.AddChild(mClockHKL);
   mClockMaster.AddChild(mClockGlobalBiso);
   mClockMaster.AddChild(mClockNbReflUsed);
   mClockMaster.AddChild(mClockFhklObsSq);
}

ScatteringData::~ScatteringData()
{
   VFN_DEBUG_MESSAGE("ScatteringData::~ScatteringData()",10)
}

void ScatteringData::SetHKL(const CrystVector_REAL &h,
                            const CrystVector_REAL &k,
                            const CrystVector_REAL &l)
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

void ScatteringData::GenHKLFullSpace2(const REAL maxSTOL,const bool unique)
{
   (*fpObjCrystInformUser)("Generating Full HKL list...");
   VFN_DEBUG_ENTRY("ScatteringData::GenHKLFullSpace2()",5)
   TAU_PROFILE("ScatteringData::GenHKLFullSpace2()","void (REAL,bool)",TAU_DEFAULT);
   if(0==mpCrystal)
   {
      throw ObjCrystException("ScatteringData::GenHKLFullSpace2() \
      no crystal assigned yet to this ScatteringData object.");
   }
   cctbx::uctbx::unit_cell uc=cctbx::uctbx::unit_cell(scitbx::af::double6(mpCrystal->GetLatticePar(0),
                                                                          mpCrystal->GetLatticePar(1),
						                          mpCrystal->GetLatticePar(2),
									  mpCrystal->GetLatticePar(3)*RAD2DEG,
									  mpCrystal->GetLatticePar(4)*RAD2DEG,
									  mpCrystal->GetLatticePar(5)*RAD2DEG));
   cctbx::miller::index_generator igen(uc,
                                this->GetCrystal().GetSpaceGroup().GetCCTbxSpg().type(),
                                !(this->IsIgnoringImagScattFact()),
                                1/(2*maxSTOL));
   if(unique)
   {
      mNbRefl=0;
      CrystVector_long H(mNbRefl);
      CrystVector_long K(mNbRefl);
      CrystVector_long L(mNbRefl);
      mMultiplicity.resize(mNbRefl);
      for(;;)
      {
         if(mNbRefl==H.numElements())
         {
            H.resizeAndPreserve(mNbRefl+100);
            K.resizeAndPreserve(mNbRefl+100);
            L.resizeAndPreserve(mNbRefl+100);
            mMultiplicity.resizeAndPreserve(mNbRefl+100);
         }
         cctbx::miller::index<> h = igen.next();
         if (h.is_zero()) break;
         H(mNbRefl)=h[0];
         K(mNbRefl)=h[1];
         L(mNbRefl)=h[2]; 
         cctbx::miller::sym_equiv_indices sei(this->GetCrystal().GetSpaceGroup().GetCCTbxSpg(),h);
         mMultiplicity(mNbRefl)=sei.multiplicity(!(this->IsIgnoringImagScattFact()));
         mNbRefl++;
      }
      H.resizeAndPreserve(mNbRefl);
      K.resizeAndPreserve(mNbRefl);
      L.resizeAndPreserve(mNbRefl);
      mMultiplicity.resizeAndPreserve(mNbRefl);
      this->SetHKL(H,K,L);
      this->SortReflectionBySinThetaOverLambda(maxSTOL);
   }
   else
   {
      mNbRefl=0;
      CrystVector_long H(mNbRefl);
      CrystVector_long K(mNbRefl);
      CrystVector_long L(mNbRefl);
      mMultiplicity.resize(mNbRefl);
      for(;;)
      {
         cctbx::miller::index<> h = igen.next();
         if (h.is_zero()) break;
         cctbx::miller::sym_equiv_indices sei(this->GetCrystal().GetSpaceGroup().GetCCTbxSpg(),h);
         for(int i=0;i<sei.multiplicity(true);i++)
         {
            cctbx::miller::index<> k = sei(i).h();
            if(mNbRefl==H.numElements())
            {
               H.resizeAndPreserve(mNbRefl+100);
               K.resizeAndPreserve(mNbRefl+100);
               L.resizeAndPreserve(mNbRefl+100);
               mMultiplicity.resizeAndPreserve(mNbRefl+100);
            }
            mMultiplicity(mNbRefl)=sei.multiplicity(!(this->IsIgnoringImagScattFact()));
            H(mNbRefl)=k[0];
            K(mNbRefl)=k[1];
            L(mNbRefl++)=k[2];
         }
      }
      H.resizeAndPreserve(mNbRefl);
      K.resizeAndPreserve(mNbRefl);
      L.resizeAndPreserve(mNbRefl);
      mMultiplicity.resizeAndPreserve(mNbRefl);
      this->SetHKL(H,K,L);
      this->SortReflectionBySinThetaOverLambda(maxSTOL);
   }
   mClockHKL.Click();
   {
      char buf [200];
      sprintf(buf,"Generating Full HKL list...Done (kept %d reflections)",(int)mNbRefl);
      (*fpObjCrystInformUser)((string)buf);
   }
   VFN_DEBUG_EXIT("ScatteringData::GenHKLFullSpace2():End",5)
}

void ScatteringData::GenHKLFullSpace(const REAL maxTheta,const bool useMultiplicity)
{
   VFN_DEBUG_ENTRY("ScatteringData::GenHKLFullSpace()",5)
   if(this->GetRadiation().GetWavelength()(0) <=.01)
   {
      throw ObjCrystException("ScatteringData::GenHKLFullSpace() \
      no wavelength assigned yet to this ScatteringData object.");;
   }
   this->GenHKLFullSpace2(sin(maxTheta)/this->GetRadiation().GetWavelength()(0),useMultiplicity);
   VFN_DEBUG_EXIT("ScatteringData::GenHKLFullSpace()",5)
}

RadiationType ScatteringData::GetRadiationType()const {return this->GetRadiation().GetRadiationType();}

void ScatteringData::SetCrystal(Crystal &crystal)
{
   VFN_DEBUG_MESSAGE("ScatteringData::SetCrystal()",5)
   if(mpCrystal!=0) mpCrystal->DeRegisterClient(*this);
   mpCrystal=&crystal;
   this->AddSubRefObj(crystal);
   crystal.RegisterClient(*this);
   mClockMaster.AddChild(mpCrystal->GetClockLatticePar());
   mClockGeomStructFact.Reset();
   mClockStructFactor.Reset();
}
const Crystal& ScatteringData::GetCrystal()const {return *mpCrystal;}
Crystal& ScatteringData::GetCrystal() {return *mpCrystal;}

long ScatteringData::GetNbRefl() const {return mNbRefl;}

const CrystVector_REAL& ScatteringData::GetH() const {return mH;}
const CrystVector_REAL& ScatteringData::GetK() const {return mK;}
const CrystVector_REAL& ScatteringData::GetL() const {return mL;}

const CrystVector_REAL& ScatteringData::GetH2Pi() const {return mH2Pi;}
const CrystVector_REAL& ScatteringData::GetK2Pi() const {return mK2Pi;}
const CrystVector_REAL& ScatteringData::GetL2Pi() const {return mH2Pi;}

const CrystVector_REAL& ScatteringData::GetReflX() const
{
   VFN_DEBUG_ENTRY("ScatteringData::GetReflX()",1)
   this->CalcSinThetaLambda();
   VFN_DEBUG_EXIT("ScatteringData::GetReflX()",1)
   return mX;
}
const CrystVector_REAL& ScatteringData::GetReflY() const
{
   VFN_DEBUG_ENTRY("ScatteringData::GetReflY()",1)
   this->CalcSinThetaLambda();
   VFN_DEBUG_EXIT("ScatteringData::GetReflY()",1)
   return mY;
}
const CrystVector_REAL& ScatteringData::GetReflZ() const
{
   VFN_DEBUG_ENTRY("ScatteringData::GetReflZ()",1)
   this->CalcSinThetaLambda();
   VFN_DEBUG_EXIT("ScatteringData::GetReflZ()",1)
   return mZ;
}

const CrystVector_REAL& ScatteringData::GetSinThetaOverLambda()const
{
   VFN_DEBUG_ENTRY("ScatteringData::GetSinThetaOverLambda()",1)
   this->CalcSinThetaLambda();
   VFN_DEBUG_EXIT("ScatteringData::GetSinThetaOverLambda()",1)
   return mSinThetaLambda;
}

const CrystVector_REAL& ScatteringData::GetTheta()const
{
   VFN_DEBUG_ENTRY("ScatteringData::GetTheta()",1)
   this->CalcSinThetaLambda();
   VFN_DEBUG_EXIT("ScatteringData::GetTheta()",1)
   return mTheta;
}

const RefinableObjClock& ScatteringData::GetClockTheta()const
{
   return mClockTheta;
}

const CrystVector_REAL& ScatteringData::GetFhklCalcSq() const
{
   VFN_DEBUG_ENTRY("ScatteringData::GetFhklCalcSq()",2)
   this->CalcStructFactor();
   if(mClockStructFactorSq>mClockStructFactor) return mFhklCalcSq;
   #ifdef __LIBCRYST_VECTOR_USE_BLITZ__
   mFhklCalcSq=pow2(mFhklCalcReal)+pow2(mFhklCalcImag);
   #else
   const REAL *pr,*pi;
   REAL *p;
   pr=mFhklCalcReal.data();
   pi=mFhklCalcImag.data();
   p=mFhklCalcSq.data();
   for(long i=0;i<mNbReflUsed;i++)
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
const CrystVector_REAL& ScatteringData::GetFhklCalcReal() const
{
   VFN_DEBUG_ENTRY("ScatteringData::GetFhklCalcReal()",2)
   this->CalcStructFactor();
   VFN_DEBUG_EXIT("ScatteringData::GetFhklCalcReal()",2)
   return mFhklCalcReal;
}

const CrystVector_REAL& ScatteringData::GetFhklCalcImag() const
{
   VFN_DEBUG_ENTRY("ScatteringData::GetFhklCalcImag()",2)
   this->CalcStructFactor();
   VFN_DEBUG_EXIT("ScatteringData::GetFhklCalcImag()",2)
   return mFhklCalcImag;
}

const CrystVector_REAL& ScatteringData::GetFhklObsSq() const
{
   return mFhklObsSq;
}

const map<const ScatteringPower*,CrystVector_REAL>& ScatteringData::GetScatteringFactor() const
{
   this->CalcScattFactor();
   return mvScatteringFactor;
}

CrystVector_REAL ScatteringData::GetWavelength()const {return this->GetRadiation().GetWavelength();}
#if 0
void ScatteringData::SetUseFastLessPreciseFunc(const bool useItOrNot)
{
   mUseFastLessPreciseFunc=useItOrNot;
   mClockGeomStructFact.Reset();
   mClockStructFactor.Reset();
}
#endif
void ScatteringData::SetIsIgnoringImagScattFact(const bool b)
{
   VFN_DEBUG_MESSAGE("ScatteringData::SetIsIgnoringImagScattFact():"<<b,10)
   mIgnoreImagScattFact=b;
   mClockGeomStructFact.Reset();
   mClockStructFactor.Reset();
}
bool ScatteringData::IsIgnoringImagScattFact() const {return mIgnoreImagScattFact;}

void ScatteringData::PrintFhklCalc(ostream &os)const
{
   VFN_DEBUG_ENTRY("ScatteringData::PrintFhklCalc()",5)
   this->GetFhklCalcSq();
   CrystVector_REAL theta;
   theta=mTheta;
   theta *= RAD2DEG;
   os <<" Number of reflections:"<<mNbRefl<<endl;
   os <<"       H        K        L     F(hkl)^2     Re(F)         Im(F)";
   os <<"        Theta       1/2d"<<endl;
   os << FormatVertVectorHKLFloats<REAL>
               (mH,mK,mL,mFhklCalcSq,mFhklCalcReal,mFhklCalcImag,theta,mSinThetaLambda,12,4);
   VFN_DEBUG_EXIT("ScatteringData::PrintFhklCalc()",5)
}

void ScatteringData::PrintFhklCalcDetail(ostream &os)const
{
   VFN_DEBUG_ENTRY("ScatteringData::PrintFhklCalcDetail()",5)
   this->GetFhklCalcSq();
   CrystVector_REAL theta;
   theta=mTheta;
   theta *= RAD2DEG;
   vector<const CrystVector_REAL *> v;
   v.push_back(&mH);
   v.push_back(&mK);
   v.push_back(&mL);
   v.push_back(&mSinThetaLambda);
   v.push_back(&theta);
   v.push_back(&mFhklCalcSq);
   v.push_back(&mFhklCalcReal);
   v.push_back(&mFhklCalcImag);
   os <<" Number of reflections:"<<mNbRefl<<endl;
   os <<"       H        K        L       1/2d        Theta       F(hkl)^2";
   os <<"     Re(F)         Im(F)       ";
   vector<CrystVector_REAL> sf;
   sf.resize(mvRealGeomSF.size()*2);
   map<const ScatteringPower*,CrystVector_REAL>::const_iterator pos=mvRealGeomSF.begin();
   long i=0;
   for(map<const ScatteringPower*,CrystVector_REAL>::const_iterator 
         pos=mvRealGeomSF.begin();pos!=mvRealGeomSF.end();++pos)
   {
      os << FormatString("Re(F)_"+pos->first->GetName(),14)
         << FormatString("Im(F)_"+pos->first->GetName(),14);
      cout<<pos->first->GetName()<<":"<<pos->first->GetForwardScatteringFactor(RAD_XRAY)<<endl;
      sf[2*i]  = mvRealGeomSF[pos->first];
      sf[2*i] *= mvScatteringFactor[pos->first];
      sf[2*i] *= mvTemperatureFactor[pos->first];
      sf[2*i+1]  = mvImagGeomSF[pos->first];
      sf[2*i+1] *= mvScatteringFactor[pos->first];
      sf[2*i+1] *= mvTemperatureFactor[pos->first];
      v.push_back(&(sf[2*i]));
      v.push_back(&(sf[2*i+1]));
      //v.push_back(mvRealGeomSF[pos->first]);
      //v.push_back(mvImagGeomSF[pos->first]);
      //v.push_back(mvScatteringFactor[pos->first]);
      i++;
   }
   os<<endl;
   os << FormatVertVectorHKLFloats<REAL>(v,12,4);
   VFN_DEBUG_EXIT("ScatteringData::PrintFhklCalcDetail()",5)
}

void ScatteringData::BeginOptimization(const bool allowApproximations,
                                       const bool enableRestraints)
{
   if(mUseFastLessPreciseFunc!=allowApproximations)
   {
      mClockGeomStructFact.Reset();
      mClockStructFactor.Reset();
      mClockMaster.Click();
   }
   mUseFastLessPreciseFunc=allowApproximations;
   this->RefinableObj::BeginOptimization(allowApproximations,enableRestraints);
}
void ScatteringData::EndOptimization()
{
   if(mOptimizationDepth==1)
   {
      if(mUseFastLessPreciseFunc==true)
      {
         mClockGeomStructFact.Reset();
         mClockStructFactor.Reset();
         mClockMaster.Click();
      }
      mUseFastLessPreciseFunc=false;
   }
   this->RefinableObj::EndOptimization();
}

void ScatteringData::SetApproximationFlag(const bool allow)
{
   if(mUseFastLessPreciseFunc!=allow)
   {
      mClockGeomStructFact.Reset();
      mClockStructFactor.Reset();
      mClockMaster.Click();
   }
   mUseFastLessPreciseFunc=allow;
   this->RefinableObj::SetApproximationFlag(allow);
}

void ScatteringData::PrepareHKLarrays()
{
   VFN_DEBUG_ENTRY("ScatteringData::PrepareHKLarrays()"<<mNbRefl<<" reflections",5)
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
   
   // If we extracted some intensities, try to keep them.
   // Do not do this if the number of reflections changed too much, if there is no crystal
   // structure associated, or if the spacegroup changed.
   bool noSpgChange=false;
   if(mpCrystal!=0) noSpgChange = mpCrystal->GetSpaceGroup().GetClockSpaceGroup()<mClockFhklObsSq;
   if( (mFhklObsSq.numElements()>0) && (abs(mFhklObsSq.numElements()-mNbRefl)<(0.1*mNbRefl)) && (noSpgChange) )
      mFhklObsSq.resizeAndPreserve(mNbRefl);
   else mFhklObsSq.resize(0);
   mClockFhklObsSq.Click();
   
   mNbReflUsed=mNbRefl;
   
   mExpectedIntensityFactor.resize(mNbRefl);
   for(long i=0;i<mNbRefl;i++)
   {
      mExpectedIntensityFactor(i)=
         mpCrystal->GetSpaceGroup().GetExpectedIntensityFactor(mH(i),mK(i),mL(i));
   }
   /*
   {
      mpCrystal->GetSpaceGroup().Print();
      for(long i=0;i<mNbRefl;i++)
      {
         cout<<mIntH(i)<<" "<<mIntK(i)<<" "<<mIntL(i)<<" "<<mExpectedIntensityFactor(i)<<endl;
      }
   }
   */

   mClockHKL.Click();
   VFN_DEBUG_EXIT("ScatteringData::PrepareHKLarrays()"<<mNbRefl<<" reflections",5)
}

void ScatteringData::SetMaxSinThetaOvLambda(const REAL max){mMaxSinThetaOvLambda=max;}
REAL ScatteringData::GetMaxSinThetaOvLambda()const{return mMaxSinThetaOvLambda;}
long ScatteringData::GetNbReflBelowMaxSinThetaOvLambda()const
{
   if(this->IsBeingRefined()) return mNbReflUsed;
   VFN_DEBUG_MESSAGE("ScatteringData::GetNbReflBelowMaxSinThetaOvLambda()",4)
   this->CalcSinThetaLambda();
   if((mNbReflUsed>0)&&(mNbReflUsed<mNbRefl))
   {
      if(  (mSinThetaLambda(mNbReflUsed  )>mMaxSinThetaOvLambda)
         &&(mSinThetaLambda(mNbReflUsed-1)<=mMaxSinThetaOvLambda)) return mNbReflUsed;
   }
   
   if((mNbReflUsed==mNbRefl)&&(mSinThetaLambda(mNbRefl-1)<=mMaxSinThetaOvLambda))
      return mNbReflUsed;
   
   long i;
   for(i=0;i<mNbRefl;i++) if(mSinThetaLambda(i)>mMaxSinThetaOvLambda) break;
   if(i!=mNbReflUsed)
   {
      mNbReflUsed=i;
      mClockNbReflUsed.Click();
      VFN_DEBUG_MESSAGE("->Changed Max sin(theta)/lambda="<<mMaxSinThetaOvLambda\
                        <<" nb refl="<<mNbReflUsed,4)
   }
   return mNbReflUsed;
}
const RefinableObjClock& ScatteringData::GetClockNbReflBelowMaxSinThetaOvLambda()const
{return mClockNbReflUsed;}

CrystVector_long ScatteringData::SortReflectionBySinThetaOverLambda(const REAL maxSTOL)
{
   TAU_PROFILE("ScatteringData::SortReflectionBySinThetaOverLambda()","void ()",TAU_DEFAULT);
   VFN_DEBUG_ENTRY("ScatteringData::SortReflectionBySinThetaOverLambda()",5)
   this->CalcSinThetaLambda();
   CrystVector_long sortedSubs;
   sortedSubs=SortSubs(mSinThetaLambda);
   CrystVector_long oldH,oldK,oldL,oldMult;
   oldH=mH;
   oldK=mK;
   oldL=mL;
   oldMult=mMultiplicity;
   long subs;
   long shift=0;
   
   //get rid of [0,0,0] reflection
   VFN_DEBUG_MESSAGE("ScatteringData::SortReflectionBySinThetaOverLambda() 1",2)
   if(0==mSinThetaLambda(sortedSubs(0)))
   {
      shift=1;
      mNbRefl -= 1;
      mH.resize(mNbRefl);
      mK.resize(mNbRefl);
      mL.resize(mNbRefl);
      mMultiplicity.resize(mNbRefl);
   }
   VFN_DEBUG_MESSAGE("ScatteringData::SortReflectionBySinThetaOverLambda() 2",2)
   for(long i=0;i<mNbRefl;i++)
   {
      subs=sortedSubs(i+shift);
      mH(i)=oldH(subs);
      mK(i)=oldK(subs);
      mL(i)=oldL(subs);
      mMultiplicity(i)=oldMult(subs);
   }
   mClockHKL.Click();
   VFN_DEBUG_MESSAGE("ScatteringData::SortReflectionBySinThetaOverLambda() 3",2)
   this->PrepareHKLarrays();
   this->CalcSinThetaLambda();
   
   VFN_DEBUG_MESSAGE("ScatteringData::SortReflectionBySinThetaOverLambda() 4",2)
   if(0<maxSTOL)
   {
      VFN_DEBUG_MESSAGE("ScatteringData::SortReflectionBySinThetaOverLambda() 5"<<maxSTOL,2)
      long maxSubs=0;
      VFN_DEBUG_MESSAGE("  "<< mIntH(maxSubs)<<" "<< mIntK(maxSubs)<<" "<< mIntL(maxSubs)<<" "<<mSinThetaLambda(maxSubs),1)
      for(maxSubs=0;mSinThetaLambda(maxSubs)<maxSTOL;maxSubs++)
      {
         VFN_DEBUG_MESSAGE("  "<< mIntH(maxSubs)<<" "<< mIntK(maxSubs)<<" "<< mIntL(maxSubs)<<" "<<mSinThetaLambda(maxSubs),1)
         if(maxSubs==(mNbRefl-1))
         {
            maxSubs=mNbRefl;
            break;
         }
      }
      if(maxSubs==mNbRefl)
      {
         VFN_DEBUG_EXIT("ScatteringData::SortReflectionBySinThetaOverLambda():"<<mNbRefl<<" reflections",5)
         return sortedSubs;
      }
      mNbRefl=maxSubs;
      mH.resizeAndPreserve(mNbRefl);
      mK.resizeAndPreserve(mNbRefl);
      mL.resizeAndPreserve(mNbRefl);
      mMultiplicity.resizeAndPreserve(mNbRefl);
      sortedSubs.resizeAndPreserve(mNbRefl);
      mClockHKL.Click();
      this->PrepareHKLarrays();
   }
   VFN_DEBUG_EXIT("ScatteringData::SortReflectionBySinThetaOverLambda():"<<mNbRefl<<" reflections",5)
   return sortedSubs;
}

CrystVector_long ScatteringData::EliminateExtinctReflections()
{
   TAU_PROFILE("ScatteringData::EliminateExtinctReflections()","void ()",TAU_DEFAULT);
   VFN_DEBUG_ENTRY("ScatteringData::EliminateExtinctReflections()",7)
   
   long nbKeptRefl=0;
   CrystVector_long subscriptKeptRefl(mNbRefl);
   subscriptKeptRefl=0;
   for(long j=0;j<mNbRefl;j++)
   {
      if( this->GetCrystal().GetSpaceGroup().IsReflSystematicAbsent(mH(j),mK(j),mL(j))==false )
         subscriptKeptRefl(nbKeptRefl++)=j;
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
   this->PrepareHKLarrays();
   VFN_DEBUG_EXIT("ScatteringData::EliminateExtinctReflections():End",7)
   return subscriptKeptRefl;
}

void ScatteringData::CalcSinThetaLambda()const
{
   if(mClockTheta>mClockMaster) return;
   if( 0 == mpCrystal) throw ObjCrystException("ScatteringData::CalcSinThetaLambda() \
      Cannot compute sin(theta)/lambda : there is no crystal affected to this \
      ScatteringData object yet.");

   if( 0 == this->GetNbRefl()) throw ObjCrystException("ScatteringData::CalcSinThetaLambda() \
      Cannot compute sin(theta)/lambda : there are no reflections !");

   if(  (mClockTheta>this->GetRadiation().GetClockWavelength()) && (mClockTheta>mClockHKL)
      &&(mClockTheta>mpCrystal->GetClockLatticePar())) return;
   
   VFN_DEBUG_ENTRY("ScatteringData::CalcSinThetaLambda()",3)
   TAU_PROFILE("ScatteringData::CalcSinThetaLambda()","void (bool)",TAU_DEFAULT);
   mSinThetaLambda.resize(mNbRefl);
   
   const CrystMatrix_REAL bMatrix= mpCrystal->GetBMatrix();
   mX.resize(this->GetNbRefl());
   mY.resize(this->GetNbRefl());
   mZ.resize(this->GetNbRefl());
   for(int i=0;i<this->GetNbRefl();i++)
   {  //:TODO: faster,nicer 
      mX(i)=bMatrix(0,0)*mH(i)+bMatrix(0,1)*mK(i)+bMatrix(0,2)*mL(i);
      mY(i)=bMatrix(1,0)*mH(i)+bMatrix(1,1)*mK(i)+bMatrix(1,2)*mL(i);
      mZ(i)=bMatrix(2,0)*mH(i)+bMatrix(2,1)*mK(i)+bMatrix(2,2)*mL(i);
   }
   //cout << bMatrix << endl << xyz<<endl;
   for(int i=0;i< (this->GetNbRefl());i++)
      mSinThetaLambda(i)=sqrt(pow(mX(i),2)+pow(mY(i),2)+pow(mZ(i),2))/2;
   if(this->GetRadiation().GetWavelengthType()!=WAVELENGTH_TOF)
   {
      if(this->GetRadiation().GetWavelength()(0) > 0)
      {
         mTheta.resize(mNbRefl);
         for(int i=0;i< (this->GetNbRefl());i++) 
         {  
            if( (mSinThetaLambda(i)*this->GetRadiation().GetWavelength()(0))>1)
            {
               //:KLUDGE: :TODO:
               mTheta(i)=M_PI;
               /*
               ofstream out("log.txt");
               out << "Error when computing Sin(theta) :"
                   << "i="<<i<<" ,mSinThetaLambda(i)="<<mSinThetaLambda(i)
                   << " ,this->GetRadiation().GetWavelength()(0)="
                   << this->GetRadiation().GetWavelength()(0) 
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
               mTheta(i)=asin(mSinThetaLambda(i)*this->GetRadiation().GetWavelength()(0));
            }
         }
      } else 
      {
         cout << "Wavelength not given in ScatteringData::CalcSinThetaLambda() !" <<endl;
         throw 0;
      }
   }
   else mTheta.resize(0);
      
   mClockTheta.Click();
   VFN_DEBUG_EXIT("ScatteringData::CalcSinThetaLambda()",3)
}

void ScatteringData::CalcScattFactor()const
{
   //if(mClockScattFactor>mClockMaster) return;
   if(  (mClockScattFactor>this->GetRadiation().GetClockWavelength()) 
      &&(mClockScattFactor>mClockHKL)
      &&(mClockScattFactor>mClockTheta)
      &&(mClockScattFactor>mpCrystal->GetClockLatticePar())
      &&(mClockThermicFact>mpCrystal->GetMasterClockScatteringPower())) return;
   TAU_PROFILE("ScatteringData::CalcScattFactor()","void (bool)",TAU_DEFAULT);
   VFN_DEBUG_ENTRY("ScatteringData::CalcScattFactor()",4)
   this->CalcResonantScattFactor();
   mvScatteringFactor.clear();
   for(int i=mpCrystal->GetScatteringPowerRegistry().GetNb()-1;i>=0;i--)
   {
      const ScatteringPower *pScattPow=&(mpCrystal->GetScatteringPowerRegistry().GetObj(i));
      mvScatteringFactor[pScattPow]=pScattPow->GetScatteringFactor(*this);
      //Directly add Fprime
      mvScatteringFactor[pScattPow]+= this->mvFprime[pScattPow];
      VFN_DEBUG_MESSAGE("->   H      K      L   sin(t/l)     f0+f'"
                        <<FormatVertVectorHKLFloats<REAL>(mH,mK,mL,mSinThetaLambda,
                                                          mvScatteringFactor[pScattPow]),1);
   }
   mClockScattFactor.Click();
   VFN_DEBUG_EXIT("ScatteringData::CalcScattFactor()",4)
}

void ScatteringData::CalcTemperatureFactor()const
{
   //if(mClockThermicFact>mClockMaster) return;
   if(  (mClockThermicFact>this->GetRadiation().GetClockWavelength())
      &&(mClockThermicFact>mClockHKL)
      &&(mClockThermicFact>mClockTheta)
      &&(mClockThermicFact>mpCrystal->GetClockLatticePar())
      &&(mClockThermicFact>mpCrystal->GetMasterClockScatteringPower())) return;
   TAU_PROFILE("ScatteringData::CalcTemperatureFactor()","void (bool)",TAU_DEFAULT);
   VFN_DEBUG_ENTRY("ScatteringData::CalcTemperatureFactor()",4)
   mvTemperatureFactor.clear();
   for(int i=mpCrystal->GetScatteringPowerRegistry().GetNb()-1;i>=0;i--)
   {
      const ScatteringPower *pScattPow=&(mpCrystal->GetScatteringPowerRegistry().GetObj(i));
      mvTemperatureFactor[pScattPow]=pScattPow->GetTemperatureFactor(*this);
   }
   mClockThermicFact.Click();
   VFN_DEBUG_EXIT("ScatteringData::CalcTemperatureFactor()",4)
}

void ScatteringData::CalcResonantScattFactor()const
{
   if(  (mClockScattFactorResonant>mpCrystal->GetMasterClockScatteringPower())
      &&(mClockScattFactorResonant>this->GetRadiation().GetClockWavelength())) return;
   VFN_DEBUG_ENTRY("ScatteringData::CalcResonantScattFactor()",4)
   TAU_PROFILE("ScatteringData::CalcResonantScattFactor()","void (bool)",TAU_DEFAULT);
   
   mvFprime.clear();
   mvFsecond.clear();
   if(this->GetRadiation().GetWavelength()(0) == 0)
   {
      VFN_DEBUG_EXIT("ScatteringData::CalcResonantScattFactor()->Lambda=0. fprime=fsecond=0",4)
      return;
   }
   else
   {
      for(int i=mpCrystal->GetScatteringPowerRegistry().GetNb()-1;i>=0;i--)
      {
         const ScatteringPower *pScattPow=&(mpCrystal->GetScatteringPowerRegistry().GetObj(i));
         mvFprime [pScattPow]=pScattPow->GetResonantScattFactReal(*this)(0);
         mvFsecond[pScattPow]=pScattPow->GetResonantScattFactImag(*this)(0);
      }
   }
   mClockScattFactorResonant.Click();
   VFN_DEBUG_EXIT("ScatteringData::CalcResonantScattFactor()",4)
}

void ScatteringData::CalcGlobalTemperatureFactor() const
{
   this->GetNbReflBelowMaxSinThetaOvLambda();//update mNbReflUsed, also recalc sin(theta)/lambda
   if(mClockGlobalTemperatureFact>mClockMaster) return;
   if(  (mClockGlobalBiso<mClockGlobalTemperatureFact)
      &&(mClockTheta     <mClockGlobalTemperatureFact)
      &&(mClockHKL       <mClockGlobalTemperatureFact)
      &&(mClockNbReflUsed<mClockGlobalTemperatureFact)) return;
   VFN_DEBUG_MESSAGE("ScatteringData::CalcGlobalTemperatureFactor()",2)
   TAU_PROFILE("ScatteringData::CalcGlobalTemperatureFactor()","void ()",TAU_DEFAULT);
   
   mGlobalTemperatureFactor.resize(mNbRefl);
   //if(true==mUseFastLessPreciseFunc) //:TODO:
   {
   }
   //else
   {
      const REAL *stol=this->GetSinThetaOverLambda().data();
      REAL *fact=mGlobalTemperatureFactor.data();
      for(long i=0;i<mNbReflUsed;i++) {*fact++ = exp(-mGlobalBiso * *stol * *stol);stol++;}
   }
   mClockGlobalTemperatureFact.Click();
}

void ScatteringData::CalcStructFactor() const
{
   this->GetNbReflBelowMaxSinThetaOvLambda();//check mNbReflUsed, also recalc sin(theta)/lambda
   if(mClockStructFactor>mClockMaster) return;
   
   //:TODO: Anisotropic Thermic factors
   //TAU_PROFILE_TIMER(timer1,"ScatteringData::CalcStructFactor1:Prepare","", TAU_FIELD);
   //TAU_PROFILE_TIMER(timer2,"ScatteringData::CalcStructFactor2:GeomStructFact","", TAU_FIELD);
   //TAU_PROFILE_TIMER(timer3,"ScatteringData::CalcStructFactor3:Scatt.Factors","", TAU_FIELD);
   //TAU_PROFILE_TIMER(timer4,"ScatteringData::CalcStructFactor4:Finish,DynCorr","", TAU_FIELD);

   //TAU_PROFILE_START(timer1);
   
   const long nbRefl=this->GetNbRefl();
   this->CalcSinThetaLambda();
   //TAU_PROFILE_STOP(timer1);
   //TAU_PROFILE_START(timer2);
   this->CalcGeomStructFactor();
   //TAU_PROFILE_STOP(timer2);
   //TAU_PROFILE_START(timer3);
   this->CalcScattFactor();
   this->CalcResonantScattFactor();
   this->CalcTemperatureFactor();
   this->CalcGlobalTemperatureFactor();
   this->CalcLuzzatiFactor();
   this->CalcStructFactVariance();
   //TAU_PROFILE_STOP(timer3);
   
   //OK, really must recompute SFs?
   VFN_DEBUG_MESSAGE("ScatteringData::CalcStructFactor():Fhkl Recalc ?"<<endl
      <<"mClockStructFactor<mClockGlobalTemperatureFact"<<(bool)(mClockStructFactor<mClockGlobalTemperatureFact)<<endl
      <<"mClockStructFactor<mClockGeomStructFact"       <<(bool)(mClockStructFactor<mClockGeomStructFact)<<endl
      <<"mClockStructFactor<mClockScattFactorResonant"  <<(bool)(mClockStructFactor<mClockScattFactorResonant)<<endl
      <<"mClockStructFactor<mClockThermicFact"          <<(bool)(mClockStructFactor<mClockThermicFact)<<endl
      <<"mClockStructFactor<mClockLuzzatiFactor"        <<(bool)(mClockStructFactor<mClockLuzzatiFactor)<<endl      
      ,2)
   if(  (mClockStructFactor>mClockGlobalTemperatureFact)
      &&(mClockStructFactor>mClockGeomStructFact)
      &&(mClockStructFactor>mClockScattFactorResonant)
      &&(mClockStructFactor>mClockThermicFact)
      &&(mClockStructFactor>mClockFhklCalcVariance)
      &&(mClockStructFactor>mClockLuzzatiFactor)) return;
   VFN_DEBUG_ENTRY("ScatteringData::CalcStructFactor()",3)
   TAU_PROFILE("ScatteringData::CalcStructFactor()","void ()",TAU_DEFAULT);
   //TAU_PROFILE_START(timer4);
   //reset Fcalc
      mFhklCalcReal.resize(nbRefl);
      mFhklCalcImag.resize(nbRefl);
      mFhklCalcReal=0;
      mFhklCalcImag=0;
   //Add all contributions
   for(map<const ScatteringPower*,CrystVector_REAL>::const_iterator pos=mvRealGeomSF.begin();
       pos!=mvRealGeomSF.end();++pos)
   {
      const ScatteringPower* pScattPow=pos->first;
      VFN_DEBUG_MESSAGE("ScatteringData::CalcStructFactor():Fhkl Recalc, "<<pScattPow->GetName(),2)
      const REAL * RESTRICT pGeomR=mvRealGeomSF[pScattPow].data();
      const REAL * RESTRICT pGeomI=mvImagGeomSF[pScattPow].data();
      const REAL * RESTRICT pScatt=mvScatteringFactor[pScattPow].data();
      const REAL * RESTRICT pTemp=mvTemperatureFactor[pScattPow].data();

      REAL * RESTRICT pReal=mFhklCalcReal.data();
      REAL * RESTRICT pImag=mFhklCalcImag.data();

      VFN_DEBUG_MESSAGE("->mvRealGeomSF[i] "
         <<mvRealGeomSF[pScattPow].numElements()<<"elements",2)
      VFN_DEBUG_MESSAGE("->mvImagGeomSF[i] "
         <<mvImagGeomSF[pScattPow].numElements()<<"elements",2)
      VFN_DEBUG_MESSAGE("->mvScatteringFactor[i]"
         <<mvScatteringFactor[pScattPow].numElements()<<"elements",1)
      VFN_DEBUG_MESSAGE("->mvTemperatureFactor[i]" 
         <<mvTemperatureFactor[pScattPow].numElements()<<"elements",1)
      VFN_DEBUG_MESSAGE("->mFhklCalcReal "<<mFhklCalcReal.numElements()<<"elements",2)
      VFN_DEBUG_MESSAGE("->mFhklCalcImag "<<mFhklCalcImag.numElements()<<"elements",2)
      VFN_DEBUG_MESSAGE("->   H      K      L   sin(t/l)     Re(F)      Im(F)      scatt      Temp->"<<pScattPow->GetName(),1)

      VFN_DEBUG_MESSAGE(FormatVertVectorHKLFloats<REAL>(mH,mK,mL,mSinThetaLambda,
                                                        mvRealGeomSF[pScattPow],
                                                        mvImagGeomSF[pScattPow],
                                                        mvScatteringFactor[pScattPow],
                                                        mvTemperatureFactor[pScattPow] 
                                                        ),1);
      if(mvLuzzatiFactor[pScattPow].numElements()>0)
      {// using maximum likelihood
         const REAL* RESTRICT pLuzzati=mvLuzzatiFactor[pScattPow].data();
         if(false==mIgnoreImagScattFact)
         {
            const REAL fsecond=mvFsecond[pScattPow];
            VFN_DEBUG_MESSAGE("->fsecond= "<<fsecond,10)
            for(long j=mNbReflUsed;j>0;j--)
            {
               VFN_DEBUG_MESSAGE("-->"<<j<<" "<<*pReal<<" "<<*pImag<<" "<<*pGeomR<<" "<<*pGeomI<<" "<<*pScatt<<" "<<*pTemp,1)
               *pReal++ += (*pGeomR   * *pScatt   - *pGeomI   * fsecond)* *pTemp * *pLuzzati;
               *pImag++ += (*pGeomI++ * *pScatt++ + *pGeomR++ * fsecond)* *pTemp++ * *pLuzzati++;
            }
         }
         else
         {
            for(long j=mNbReflUsed;j>0;j--)
            {
               *pReal++ += *pGeomR++  * *pTemp   * *pScatt   * *pLuzzati;
               *pImag++ += *pGeomI++  * *pTemp++ * *pScatt++ * *pLuzzati++;
            }
         }
         VFN_DEBUG_MESSAGE("ScatteringData::CalcStructFactor():"<<mIgnoreImagScattFact
                           <<",f\"="<<mvFsecond[pScattPow]<<endl<<
                           FormatVertVectorHKLFloats<REAL>(mH,mK,mL,mSinThetaLambda,
                                                           mvRealGeomSF[pScattPow],
                                                           mvImagGeomSF[pScattPow],
                                                           mvScatteringFactor[pScattPow],
                                                           mvTemperatureFactor[pScattPow], 
                                                           mvLuzzatiFactor[pScattPow],
                                                           mFhklCalcReal, 
                                                           mFhklCalcImag
                                                           ),2);
      }
      else
      { 
         if(false==mIgnoreImagScattFact)
         {
            const REAL fsecond=mvFsecond[pScattPow];
            VFN_DEBUG_MESSAGE("->fsecond= "<<fsecond,2)
            for(long j=mNbReflUsed;j>0;j--)
            {
               *pReal += (*pGeomR   * *pScatt - *pGeomI * fsecond)* *pTemp;
               *pImag += (*pGeomI   * *pScatt + *pGeomR * fsecond)* *pTemp;
               VFN_DEBUG_MESSAGE("-->"<<j<<" "<<*pReal<<" "<<*pImag<<" "<<*pGeomR<<" "<<*pGeomI<<" "<<*pScatt<<" "<<*pTemp,1)
               pGeomR++;pGeomI++;pTemp++;pScatt++;pReal++;pImag++;
            }
         }
         else
         {
            for(long j=mNbReflUsed;j>0;j--)
            {
               *pReal++ += *pGeomR++  * *pTemp * *pScatt;
               *pImag++ += *pGeomI++  * *pTemp++ * *pScatt++;
            }
         }
         VFN_DEBUG_MESSAGE(FormatVertVectorHKLFloats<REAL>(mH,mK,mL,mSinThetaLambda,
                                                            mvRealGeomSF[pScattPow],
                                                            mvImagGeomSF[pScattPow],
                                                            mvScatteringFactor[pScattPow],
                                                            mvTemperatureFactor[pScattPow], 
                                                            mFhklCalcReal, 
                                                            mFhklCalcImag
                                                            ),2);
      }
   }
   //TAU_PROFILE_STOP(timer4);
   {
      this->CalcGlobalTemperatureFactor();
      if(mGlobalTemperatureFactor.numElements()>0)
      {//else for some reason it's useless
         REAL *pReal=mFhklCalcReal.data();
         REAL *pImag=mFhklCalcImag.data();
         const REAL *pTemp=mGlobalTemperatureFactor.data();
         for(long j=0;j<mNbReflUsed;j++)
         {
            *pReal++ *= *pTemp;
            *pImag++ *= *pTemp++;
         }
      }
   }
   mClockStructFactor.Click();

   VFN_DEBUG_EXIT("ScatteringData::CalcStructFactor()",3)
}

void ScatteringData::CalcGeomStructFactor() const
{
   // This also updates the ScattCompList if necessary.
   const ScatteringComponentList *pScattCompList
      =&(this->GetCrystal().GetScatteringComponentList());
   if(  (mClockGeomStructFact>mpCrystal->GetClockScattCompList())
      &&(mClockGeomStructFact>mClockHKL)
      &&(mClockGeomStructFact>mClockNbReflUsed)
      &&(mClockGeomStructFact>mpCrystal->GetMasterClockScatteringPower())) return;
   TAU_PROFILE("ScatteringData::GeomStructFactor()","void (Vx,Vy,Vz,data,M,M,bool)",TAU_DEFAULT);
   VFN_DEBUG_ENTRY("ScatteringData::GeomStructFactor(Vx,Vy,Vz,...)",3)
   VFN_DEBUG_MESSAGE("-->Using fast functions:"<<mUseFastLessPreciseFunc,2)
   VFN_DEBUG_MESSAGE("-->Number of translation vectors:"
      <<this->GetCrystal().GetSpaceGroup().GetNbTranslationVectors()-1,2)
   VFN_DEBUG_MESSAGE("-->Has an inversion Center:"
      <<this->GetCrystal().GetSpaceGroup().HasInversionCenter(),2)
   VFN_DEBUG_MESSAGE("-->Number of symetry operations (w/o transl&inv cent.):"\
                     <<this->GetCrystal().GetSpaceGroup().GetNbSymmetrics(true,true),2)
   VFN_DEBUG_MESSAGE("-->Number of Scattering Components :"
      <<this->GetCrystal().GetScatteringComponentList().GetNbComponent(),2)
   VFN_DEBUG_MESSAGE("-->Number of reflections:"
      <<this->GetNbRefl()<<" (actually used:"<<mNbReflUsed<<")",2)
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
      const SpaceGroup *pSpg=&(this->GetCrystal().GetSpaceGroup());
      
      const int nbSymmetrics=pSpg->GetNbSymmetrics(true,true);
      const int nbTranslationVectors=pSpg->GetNbTranslationVectors();
      const long nbComp=pScattCompList->GetNbComponent();
      const int nbRefl=this->GetNbRefl();
      const std::vector<SpaceGroup::TRx> *pTransVect=&(pSpg->GetTranslationVectors());
      CrystMatrix_REAL allCoords(nbSymmetrics,3);
      CrystVector_REAL tmpVect(nbRefl);
      
      CrystVector_long intVect(nbRefl);//not used if mUseFastLessPreciseFunc==false
      
      // which scattering powers are actually used ?
      map<const ScatteringPower*,bool> vUsed;
      // Add existing previously used scattering power to the test;
      for(map<const ScatteringPower*,CrystVector_REAL>::const_iterator pos=mvRealGeomSF.begin();pos!=mvRealGeomSF.end();++pos)
         vUsed[pos->first]=false;// this will be changed to true later if they are actually used
      
      for(int i=mpCrystal->GetScatteringPowerRegistry().GetNb()-1;i>=0;i--)
      {// Here we make sure scattering power that only contribute ghost atoms are taken into account
         const ScatteringPower*pow=&(mpCrystal->GetScatteringPowerRegistry().GetObj(i));
         if(pow->GetMaximumLikelihoodNbGhostAtom()>0) vUsed[pow]=true;
         else vUsed[pow]=false;
      }
      for(long i=0;i<nbComp;i++)
         vUsed[(*pScattCompList)(i).mpScattPow]=true;
      //Resize all arrays and set them to 0
      for(map<const ScatteringPower*,bool>::const_iterator pos=vUsed.begin();pos!=vUsed.end();++pos)
      {
         if(pos->second)
         {// this will create the entry if it does not already exist
            mvRealGeomSF[pos->first].resize(nbRefl);
            mvImagGeomSF[pos->first].resize(nbRefl);
            mvRealGeomSF[pos->first]=0;
            mvImagGeomSF[pos->first]=0;
         }
         else
         {// erase entries that are not useful any more (e.g. ScatteringPower that were
          // used but are not any more).
            map<const ScatteringPower*,CrystVector_REAL>::iterator
               poubelle=mvRealGeomSF.find(pos->first);
            if(poubelle!=mvRealGeomSF.end()) mvRealGeomSF.erase(poubelle);
            poubelle=mvImagGeomSF.find(pos->first);
            if(poubelle!=mvImagGeomSF.end()) mvImagGeomSF.erase(poubelle);
         }
      }
      
      for(long i=0;i<nbComp;i++)
      {
         VFN_DEBUG_MESSAGE("ScatteringData::GeomStructFactor(),comp"<<i,3)
         const REAL x=(*pScattCompList)(i).mX;
         const REAL y=(*pScattCompList)(i).mY;
         const REAL z=(*pScattCompList)(i).mZ;
         const ScatteringPower *pScattPow=(*pScattCompList)(i).mpScattPow;
         REAL centrMult=1.0;
         if(true==pSpg->HasInversionCenter()) centrMult=2.0;
         const REAL popu= (*pScattCompList)(i).mOccupancy 
                         *(*pScattCompList)(i).mDynPopCorr
                         *centrMult;
         
         allCoords=pSpg->GetAllSymmetrics(x,y,z,true,true);
         if((true==pSpg->HasInversionCenter()) && (false==pSpg->IsInversionCenterAtOrigin()))
         {
            const REAL STBF=2.*pSpg->GetCCTbxSpg().inv_t().den();
            for(int j=0;j<nbSymmetrics;j++)
            {
               //The phase of the structure factor will be wrong
               //This is fixed a bit further...
               allCoords(j,0) -= ((REAL)pSpg->GetCCTbxSpg().inv_t()[0])/STBF;
               allCoords(j,1) -= ((REAL)pSpg->GetCCTbxSpg().inv_t()[1])/STBF;
               allCoords(j,2) -= ((REAL)pSpg->GetCCTbxSpg().inv_t()[2])/STBF;
            }
         }
         for(int j=0;j<nbSymmetrics;j++)
         {
            VFN_DEBUG_MESSAGE("ScatteringData::GeomStructFactor(),comp #"<<i<<", sym #"<<j,3)
            
            if(mUseFastLessPreciseFunc==true)
            {
               REAL * RESTRICT rrsf=mvRealGeomSF[pScattPow].data();
               REAL * RESTRICT iisf=mvImagGeomSF[pScattPow].data();
               
               const long intX=(long)(allCoords(j,0)*sLibCrystNbTabulSine);
               const long intY=(long)(allCoords(j,1)*sLibCrystNbTabulSine);
               const long intZ=(long)(allCoords(j,2)*sLibCrystNbTabulSine);

               const long * RESTRICT intH=mIntH.data();
               const long * RESTRICT intK=mIntK.data();
               const long * RESTRICT intL=mIntL.data();

               register long * RESTRICT tmpInt=intVect.data();
               // :KLUDGE: using a AND to bring back within [0;sLibCrystNbTabulSine[ may
               // not be portable, depending on the model used to represent signed integers
               // a test should be added to throw up in that case.
               //
               // This work if we are using "2's complement" to represent negative numbers,
               // but not with a "sign magnitude" approach
               for(int jj=mNbReflUsed;jj>0;jj--)
                *tmpInt++ = (*intH++ * intX + *intK++ * intY + *intL++ *intZ)
                              &sLibCrystNbTabulSineMASK;
               if(false==pSpg->HasInversionCenter())
               {

                  tmpInt=intVect.data();
                  for(int jj=mNbReflUsed;jj>0;jj--)
                  {
                     const REAL *pTmp=&spLibCrystTabulCosineSine[*tmpInt++ <<1];
                     *rrsf++ += popu * *pTmp++;
                     *iisf++ += popu * *pTmp;
                  }

               }
               else
               {
                  tmpInt=intVect.data();
                  for(int jj=mNbReflUsed;jj>0;jj--)
                     *rrsf++ += popu * spLibCrystTabulCosine[*tmpInt++];
               }
            }
            else
            {
               const REAL x=allCoords(j,0);
               const REAL y=allCoords(j,1);
               const REAL z=allCoords(j,2);
               const register REAL *hh=mH2Pi.data();
               const register REAL *kk=mK2Pi.data();
               const register REAL *ll=mL2Pi.data();
               register REAL *tmp=tmpVect.data();
               
               
               for(int jj=0;jj<mNbReflUsed;jj++) *tmp++ = *hh++ * x + *kk++ * y + *ll++ *z;
               
               REAL *sf=mvRealGeomSF[pScattPow].data();
               tmp=tmpVect.data();
               
               for(int jj=0;jj<mNbReflUsed;jj++) *sf++ += popu * cos(*tmp++);
               
               if(false==pSpg->HasInversionCenter()) 
               {
                  sf=mvImagGeomSF[pScattPow].data();
                  tmp=tmpVect.data();
                  for(int jj=0;jj<mNbReflUsed;jj++) *sf++ += popu * sin(*tmp++);
               }
            }
         }
      }//for all components...
      
      if(nbTranslationVectors > 1)
      {
         tmpVect=1;
         if( (pSpg->GetSpaceGroupNumber()>= 143) && (pSpg->GetSpaceGroupNumber()<= 167))
         {//Special case for trigonal groups R3,...
            REAL * RESTRICT p1=tmpVect.data();
            const register REAL * RESTRICT hh=mH2Pi.data();
            const register REAL * RESTRICT kk=mK2Pi.data();
            const register REAL * RESTRICT ll=mL2Pi.data();
            for(long j=mNbReflUsed;j>0;j--) *p1++ += 2*cos((*hh++ - *kk++ - *ll++)/3.);
         }
         else
         {
            for(int j=1;j<nbTranslationVectors;j++)
            {
               const REAL x=(*pTransVect)[j].tr[0];
               const REAL y=(*pTransVect)[j].tr[1];
               const REAL z=(*pTransVect)[j].tr[2];
               REAL *p1=tmpVect.data();
               const register REAL *hh=mH2Pi.data();
               const register REAL *kk=mK2Pi.data();
               const register REAL *ll=mL2Pi.data();
               for(long j=mNbReflUsed;j>0;j--) *p1++ += cos(*hh++ *x + *kk++ *y + *ll++ *z );
            }
         }
         for(map<const ScatteringPower*,CrystVector_REAL>::iterator
               pos=mvRealGeomSF.begin();pos!=mvRealGeomSF.end();++pos) 
                  pos->second *= tmpVect;
         
         if(false==pSpg->HasInversionCenter())
            for(map<const ScatteringPower*,CrystVector_REAL>::iterator
                  pos=mvImagGeomSF.begin();pos!=mvImagGeomSF.end();++pos)
                     pos->second *= tmpVect;
      }
      if(true==pSpg->HasInversionCenter())
      {
         // we already multiplied real geom struct factor by 2
         if(false==pSpg->IsInversionCenterAtOrigin())
         {
            VFN_DEBUG_MESSAGE("ScatteringData::GeomStructFactor(Vx,Vy,Vz):\
               Inversion Center not at the origin...",2)
            //fix the phase of each reflection when the inversion center is not
            //at the origin, using :
            // Re(F) = RSF*cos(2pi(h*Xc+k*Yc+l*Zc))
            // Re(F) = RSF*sin(2pi(h*Xc+k*Yc+l*Zc))
            //cout << "Glop Glop"<<endl;
            const REAL STBF=2*pSpg->GetCCTbxSpg().inv_t().den();
            {
               const REAL xc=((REAL)pSpg->GetCCTbxSpg().inv_t()[0])/STBF;
               const REAL yc=((REAL)pSpg->GetCCTbxSpg().inv_t()[1])/STBF;
               const REAL zc=((REAL)pSpg->GetCCTbxSpg().inv_t()[2])/STBF;
               #ifdef __LIBCRYST_VECTOR_USE_BLITZ__
               tmpVect = mH2Pi() * xc + mK2PI() * yc + mL2PI() * zc;
               #else
               {
                  const REAL * RESTRICT hh=mH2Pi.data();
                  const REAL * RESTRICT kk=mK2Pi.data();
                  const REAL * RESTRICT ll=mL2Pi.data();
                  REAL * RESTRICT ttmpVect=tmpVect.data();
                  for(long ii=mNbReflUsed;ii>0;ii--) 
                     *ttmpVect++ = *hh++ * xc + *kk++ * yc + *ll++ * zc;
               }
               #endif
            }
            CrystVector_REAL cosTmpVect;
            CrystVector_REAL sinTmpVect;
            cosTmpVect=cos(tmpVect);
            sinTmpVect=sin(tmpVect);
            
            map<const ScatteringPower*,CrystVector_REAL>::iterator posi=mvImagGeomSF.begin();
            map<const ScatteringPower*,CrystVector_REAL>::iterator posr=mvRealGeomSF.begin();
            for(;posi!=mvImagGeomSF.end();)
            {
               posi->second = posr->second;
               posi->second *= sinTmpVect;
               posr->second *= cosTmpVect;
               posi++;posr++;
            }
         }
      }

   }
   //cout << FormatVertVector<REAL>(*mvRealGeomSF,*mvImagGeomSF)<<endl;
   mClockGeomStructFact.Click();
   VFN_DEBUG_EXIT("ScatteringData::GeomStructFactor(Vx,Vy,Vz,...)",3)
}

void ScatteringData::CalcLuzzatiFactor()const
{
   // Assume this is  called by ScatteringData::CalcStructFactor()
   // and that we already have computed geometrical structure factors
   VFN_DEBUG_ENTRY("ScatteringData::CalcLuzzatiFactor",3)
   bool useLuzzati=false;
   for(map<const ScatteringPower*,CrystVector_REAL>::const_iterator
       pos=mvRealGeomSF.begin();pos!=mvRealGeomSF.end();++pos)
   {
      if(pos->first->GetMaximumLikelihoodPositionError()!=0)
      {
         useLuzzati=true;
         break;
      }
   }
   if(!useLuzzati)
   {
      mvLuzzatiFactor.clear();
      VFN_DEBUG_EXIT("ScatteringData::CalcLuzzatiFactor(): not needed, no positionnal errors",3)
      return;
   }
   bool recalc=false;
   if(  (mClockTheta     >mClockLuzzatiFactor)
      ||(mClockGeomStructFact>mClockLuzzatiFactor)//checks if occupancies changed
      ||(mClockNbReflUsed>mClockLuzzatiFactor))
   {
      //if(mClockTheta     >mClockLuzzatiFactor)cout<<"1"<<endl;
      //if(mClockGeomStructFact>mClockLuzzatiFactor)cout<<"2"<<endl;
      //if(mClockNbReflUsed>mClockLuzzatiFactor)cout<<"3"<<endl;
      recalc=true;
   }
   else
   {
      for(int i=mpCrystal->GetScatteringPowerRegistry().GetNb()-1;i>=0;i--)
      {
         if(mpCrystal->GetScatteringPowerRegistry().GetObj(i)
            .GetMaximumLikelihoodParClock()>mClockLuzzatiFactor)
         {
            recalc=true;
            break;
         }
      }
   }
   
   if(false==recalc)
   {
      VFN_DEBUG_EXIT("ScatteringData::CalcLuzzatiFactor(): no recalc needed",3)
      return;
   }
   TAU_PROFILE("ScatteringData::CalcLuzzatiFactor()","void ()",TAU_DEFAULT);
   
   for(int i=mpCrystal->GetScatteringPowerRegistry().GetNb()-1;i>=0;i--)
   {
      const ScatteringPower* pScattPow=&(mpCrystal->GetScatteringPowerRegistry().GetObj(i));
      if(0 == pScattPow->GetMaximumLikelihoodPositionError())
      {
         mvLuzzatiFactor[pScattPow].resize(0);
      }
      else
      {
         mvLuzzatiFactor[pScattPow].resize(mNbRefl);
         const REAL b=-(8*M_PI*M_PI)* pScattPow->GetMaximumLikelihoodPositionError()
                                    * pScattPow->GetMaximumLikelihoodPositionError();
         const REAL *stol=this->GetSinThetaOverLambda().data();
         REAL *fact=mvLuzzatiFactor[pScattPow].data();
         for(long j=0;j<mNbReflUsed;j++) {*fact++ = exp(b * *stol * *stol);stol++;}
         VFN_DEBUG_MESSAGE("ScatteringData::CalcLuzzatiFactor():"<<pScattPow->GetName()<<endl<<
                           FormatVertVectorHKLFloats<REAL>(mH,mK,mL,mSinThetaLambda,
                           mvRealGeomSF[pScattPow],mvImagGeomSF[pScattPow],
                           mvScatteringFactor[pScattPow],mvLuzzatiFactor[pScattPow]
                           ),2);
      }
   }
   mClockLuzzatiFactor.Click();
   VFN_DEBUG_EXIT("ScatteringData::CalcLuzzatiFactor(): no recalc needed",3)
}

void ScatteringData::CalcStructFactVariance()const
{
   // this is called by CalcStructFactor(), after the calculation of the structure factors,
   // and the recomputation of Luzzati factors has already been asked
   // So we only recompute if these clocks have changed.
   //
   // The Crystal::mMasterClockScatteringPower will tell the last time the number of ghost
   // atoms has been changed in any of the scattpow.
   
   if(  (mClockFhklCalcVariance>mClockLuzzatiFactor)
      &&(mClockFhklCalcVariance>mClockStructFactor)
      &&(mClockFhklCalcVariance>mpCrystal->GetMasterClockScatteringPower())) return;

   bool hasGhostAtoms=false;
   for(map<const ScatteringPower*,CrystVector_REAL>::const_iterator
       pos=mvRealGeomSF.begin();pos!=mvRealGeomSF.end();++pos)
   {
      if(pos->first->GetMaximumLikelihoodNbGhostAtom()!=0)
      {
         hasGhostAtoms=true;
         break;
      }
   }

   if( (0==mvLuzzatiFactor.size())&&(!hasGhostAtoms))
   {
      mFhklCalcVariance.resize(0);
      return;
   }
   VFN_DEBUG_ENTRY("ScatteringData::CalcStructFactVariance()",3)
   TAU_PROFILE("ScatteringData::CalcStructFactVariance()","void ()",TAU_DEFAULT);
   bool needVar=false;
   
   map<const ScatteringPower*,REAL> vComp;
   {
      const ScatteringComponentList *pList= & (this->GetCrystal().GetScatteringComponentList());
      const long nbComp=pList->GetNbComponent();
      const ScatteringComponent *pComp;
      for(long i=0;i<nbComp;i++)
      {
         pComp=&((*pList)(i));
         vComp[pComp->mpScattPow]=0;
      }
      for(long i=0;i<nbComp;i++)
      {
         pComp=&((*pList)(i));
         vComp[pComp->mpScattPow]+= pComp->mOccupancy * pComp->mDynPopCorr;
      }
      for(map<const ScatteringPower*,REAL>::iterator
          pos=vComp.begin();pos!=vComp.end();++pos)
            pos->second *= this->GetCrystal().GetSpaceGroup().GetNbSymmetrics();
   }
   // Ghost atoms
   map<const ScatteringPower*,REAL> vGhost;
   {
      const long nbScattPow=mpCrystal->GetScatteringPowerRegistry().GetNb();
      const long mult=this->GetCrystal().GetSpaceGroup().GetNbSymmetrics();
      for(int i=0;i<nbScattPow;i++)
      {
         const ScatteringPower* pow=&(mpCrystal->GetScatteringPowerRegistry().GetObj(i));
         const REAL nb=pow->GetMaximumLikelihoodNbGhostAtom();
         vGhost[pow]=nb*mult;
      }
   }

   if(mFhklCalcVariance.numElements() == mNbRefl)
   {
      REAL *pVar=mFhklCalcVariance.data();
      for(long j=0;j<mNbReflUsed;j++) *pVar++ = 0;
   }
   
   for(int i=mpCrystal->GetScatteringPowerRegistry().GetNb()-1;i>=0;i--)
   {
      const ScatteringPower* pScattPow=&(mpCrystal->GetScatteringPowerRegistry().GetObj(i));
      if(  (mvLuzzatiFactor[pScattPow].numElements()==0)
         &&(vGhost[pScattPow]==0)) continue;
      needVar=true;
      if(mFhklCalcVariance.numElements() != mNbRefl)
      {
         mFhklCalcVariance.resize(mNbRefl);
         REAL *pVar=mFhklCalcVariance.data();
         for(long j=0;j<mNbReflUsed;j++) *pVar++ = 0;
      }
      // variance on real & imag parts of the structure factor
      const REAL *pScatt=mvScatteringFactor[pScattPow].data();
      const int  *pExp=mExpectedIntensityFactor.data();
      REAL *pVar=mFhklCalcVariance.data();
      if(mvLuzzatiFactor[pScattPow].numElements()==0)
      {
         const REAL nbghost=vGhost[pScattPow];
         for(long j=0;j<mNbReflUsed;j++)
         {
            *pVar++ += *pExp++ * *pScatt * *pScatt * nbghost;
            pScatt++; 
         }
      }
      else
      {
         const REAL *pLuz=mvLuzzatiFactor[pScattPow].data();
         const REAL occ=vComp[pScattPow];
         const REAL nbghost=vGhost[pScattPow];
         for(long j=0;j<mNbReflUsed;j++)
         {
            *pVar++ += *pExp++ * *pScatt * *pScatt * ( occ*(1 - *pLuz * *pLuz) + nbghost);
            pScatt++; pLuz++;
         }
      }
   }
   if(false == needVar) mFhklCalcVariance.resize(0);

   mClockFhklCalcVariance.Click();
   VFN_DEBUG_EXIT("ScatteringData::CalcStructFactVariance()",3)
}

}//namespace ObjCryst
