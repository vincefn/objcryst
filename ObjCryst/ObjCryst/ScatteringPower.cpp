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
#include <cmath>
#include <typeinfo>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <stdio.h> //for sprintf()


#include "cctbx/eltbx/xray_scattering.h"
#include "cctbx/eltbx/tiny_pse.h"
#include "cctbx/eltbx/icsd_radii.h"
#include "cctbx/eltbx/covalent_radii.h"
#include "cctbx/eltbx/henke.h"
#include "cctbx/eltbx/neutron.h"

#include "ObjCryst/ObjCryst/ScatteringPower.h"
#include "ObjCryst/Quirks/VFNStreamFormat.h"
#include "ObjCryst/Quirks/VFNDebug.h"
#include "ObjCryst/ObjCryst/Colours.h" 

#ifdef __WX__CRYST__
   #include "ObjCryst/wxCryst/wxScatteringPower.h"
#endif

namespace ObjCryst
{

const RefParType *gpRefParTypeScattPow=0;
const RefParType *gpRefParTypeScattPowResonant=0;
const RefParType *gpRefParTypeScattPowTemperature=0;
const RefParType *gpRefParTypeScattPowTemperatureIso=0;
const RefParType *gpRefParTypeScattPowTemperatureAniso=0;
long NiftyStaticGlobalObjectsInitializer_ScatteringPower::mCount=0;
//######################################################################
//
//      Bij to Betaij conversion
//
//######################################################################
CrystVector_REAL Bij2Betaij(const CrystVector_REAL &Bij, const UnitCell &cell)
{
   // Willis & Pryor,p 101: (betaij) = 2*pi^2 * transpose(cell.Bmatrix) * (Bij) * cell.Bmatrix
   // :TODO: this needs to be checked before being used
   const REAL B11=Bij(0);
   const REAL B22=Bij(1);
   const REAL B33=Bij(2);
   const REAL B12=Bij(3);
   const REAL B13=Bij(4);
   const REAL B23=Bij(5);
   CrystMatrix_REAL B(3,3);
   B(0,0)=B11;
   B(0,1)=B12;
   B(0,1)=B13;
   B(1,0)=B12;
   B(1,1)=B22;
   B(1,1)=B23;
   B(2,0)=B13;
   B(2,1)=B23;
   B(2,1)=B33;
   CrystMatrix_REAL b(3,3);
   b=cell.GetBMatrix().transpose().Mult(B.Mult(cell.GetBMatrix()));
   b*=2*M_PI*M_PI;
}

//######################################################################
//
//      SCATTERING POWER
//
//######################################################################
ObjRegistry<ScatteringPower> gScatteringPowerRegistry("Global ScatteringPower Registry");

ScatteringPower::ScatteringPower():mDynPopCorrIndex(0),mBiso(1.0),mIsIsotropic(true),
mMaximumLikelihoodNbGhost(0),mFormalCharge(0.0)
{
   VFN_DEBUG_MESSAGE("ScatteringPower::ScatteringPower():"<<mName,5)
   mBeta.resize(6);
   mBeta = 0;
   mB.resize(6);
   mB = 0;
   gScatteringPowerRegistry.Register(*this);
   this->Init();
   mClockMaster.AddChild(mClock);
   mClockMaster.AddChild(mMaximumLikelihoodParClock);
}
ScatteringPower::ScatteringPower(const ScatteringPower& old):
mDynPopCorrIndex(old.mDynPopCorrIndex),mBiso(old.mBiso),mIsIsotropic(old.mIsIsotropic),
mBeta(old.mBeta),mB(old.mB),
mFormalCharge(old.mFormalCharge)
{
   VFN_DEBUG_MESSAGE("ScatteringPower::ScatteringPower(&old):"<<mName,5)
   gScatteringPowerRegistry.Register(*this);
   this->Init();
   mMaximumLikelihoodPositionError=old.mMaximumLikelihoodPositionError;
   mMaximumLikelihoodNbGhost=old.mMaximumLikelihoodNbGhost;
   mClockMaster.AddChild(mClock);
   mClockMaster.AddChild(mMaximumLikelihoodParClock);
}
ScatteringPower::~ScatteringPower()
{
   VFN_DEBUG_MESSAGE("ScatteringPower::~ScatteringPower():"<<mName,5)
   gScatteringPowerRegistry.DeRegister(*this);
}

const string& ScatteringPower::GetClassName()const
{
   const static string className="ScatteringPower";
   return className;
}

void ScatteringPower::operator=(const ScatteringPower& rhs)
{
   VFN_DEBUG_MESSAGE("ScatteringPower::operator=():"<<mName,2)
   mDynPopCorrIndex=rhs.mDynPopCorrIndex;
   mBiso=rhs.mBiso;
   mIsIsotropic=rhs.mIsIsotropic;
   mBeta=rhs.mBeta;
   mB=rhs.mB;
}

bool ScatteringPower::IsScatteringFactorAnisotropic()const{return false;}
bool ScatteringPower::IsTemperatureFactorAnisotropic()const{return false;}
bool ScatteringPower::IsResonantScatteringAnisotropic()const{return false;}

const string& ScatteringPower::GetSymbol() const {return this->GetName();}
REAL ScatteringPower::GetBiso() const {return mBiso;}
REAL& ScatteringPower::GetBiso() {mClock.Click();return mBiso;}
void ScatteringPower::SetBiso(const REAL newB) { mClock.Click();mBiso=newB;mIsIsotropic=true;}
REAL ScatteringPower::GetBij(const size_t &i, const size_t &j) const
{
    size_t idx = 0;
    if(i == j)
    {
        idx = i - 1;
    }
    else
    {
        idx = i + j;
    }
    return this->GetBij(idx);
}
REAL ScatteringPower::GetBij(const size_t &idx) const
{
    return mB(idx);
}
void ScatteringPower::SetBij(const size_t &i, const size_t &j, const REAL newB)
{
    size_t idx = 0;
    if(i == j)
    {
        idx = i - 1;
    }
    else
    {
        idx = i + j;
    }
    this->SetBij(idx, newB);
}
void ScatteringPower::SetBij(const size_t &idx, const REAL newB)
{
    mClock.Click();
    mIsIsotropic=false;
    mB(idx) = newB;
}
bool ScatteringPower::IsIsotropic() const {return mIsIsotropic;}
long ScatteringPower::GetDynPopCorrIndex() const {return mDynPopCorrIndex;}
long ScatteringPower::GetNbScatteringPower()const {return gScatteringPowerRegistry.GetNb();}
const RefinableObjClock& ScatteringPower::GetLastChangeClock()const {return mClock;}

const string& ScatteringPower::GetColourName()const{ return mColourName;}
const float* ScatteringPower::GetColourRGB()const{ return mColourRGB;}
void ScatteringPower::SetColour(const string& colourName)
{
   mColourName=colourName;
   this->InitRGBColour();
}
void ScatteringPower::SetColour(const float r,const float g,const float b)
{
   mColourRGB[0]=r;
   mColourRGB[1]=g;
   mColourRGB[2]=b;
}
void ScatteringPower::GetGeneGroup(const RefinableObj &obj,
                                CrystVector_uint & groupIndex,
                                unsigned int &first) const
{
   // One group for all parameters
   unsigned int index=0;
   VFN_DEBUG_MESSAGE("ScatteringPower::GetGeneGroup()",4)
   for(long i=0;i<obj.GetNbPar();i++)
      for(long j=0;j<this->GetNbPar();j++)
         if(&(obj.GetPar(i)) == &(this->GetPar(j)))
         {
            if(index==0) index=first++;
            groupIndex(i)=index;
         }
}

REAL ScatteringPower::GetMaximumLikelihoodPositionError()const 
{return mMaximumLikelihoodPositionError;}

const RefinableObjClock& ScatteringPower::GetMaximumLikelihoodParClock()const
{return mMaximumLikelihoodParClock;}

void ScatteringPower::SetMaximumLikelihoodPositionError(const REAL mle) 
{
   if(mle!=mMaximumLikelihoodPositionError)
   {
      mMaximumLikelihoodPositionError=mle;
      mMaximumLikelihoodParClock.Click();
   }
}

REAL ScatteringPower::GetMaximumLikelihoodNbGhostAtom()const
{return mMaximumLikelihoodNbGhost;}

void ScatteringPower::SetMaximumLikelihoodNbGhostAtom(const REAL nb)
{
   if(nb!=mMaximumLikelihoodNbGhost)
   {
      mMaximumLikelihoodNbGhost=nb;
      mMaximumLikelihoodParClock.Click();
   }
}

REAL ScatteringPower::GetFormalCharge()const{return mFormalCharge;}
void ScatteringPower::SetFormalCharge(const REAL charge)
{mFormalCharge=charge;}

void ScatteringPower::Init()
{
   VFN_DEBUG_MESSAGE("ScatteringPower::Init():"<<mName,2)
   mColourName="White";
   mMaximumLikelihoodPositionError=0;
   mMaximumLikelihoodNbGhost=0;
   VFN_DEBUG_MESSAGE("ScatteringPower::Init():End",2)
}
void ScatteringPower::InitRGBColour()
{
   VFN_DEBUG_MESSAGE("ScatteringPower::InitRGBColour()",2)
   for(long i=0;;)
   {
      if(gPOVRayColours[i].mName==mColourName)
      {
         mColourRGB[0]=gPOVRayColours[i].mRGB[0];
         mColourRGB[1]=gPOVRayColours[i].mRGB[1];
         mColourRGB[2]=gPOVRayColours[i].mRGB[2];
         break;
      }
      i++;
      if(strncmp(gPOVRayColours[i].mName,"",3)==0)
      {//could not find colour !
         cout << "Could not find colour:"<<mColourName<<" for ScatteringPower "<<mName<<endl;
         mColourRGB[0]=1;
         mColourRGB[1]=1;
         mColourRGB[2]=1;
         break;
      }
   }
   VFN_DEBUG_MESSAGE("->RGBColour:"<<mColourName<<mColourRGB[0]<<" "<<mColourRGB[1]<<" "<<mColourRGB[2],2)
}

//######################################################################
//
//      SCATTERING POWER ATOM
//
//######################################################################
ObjRegistry<ScatteringPowerAtom> 
   gScatteringPowerAtomRegistry("Global ScatteringPowerAtom Registry");

ScatteringPowerAtom::ScatteringPowerAtom():
ScatteringPower(),mSymbol(""),mAtomicNumber(0),mpGaussian(0)
{
   VFN_DEBUG_MESSAGE("ScatteringPowerAtom::ScatteringPowerAtom():"<<mName,5)
   gScatteringPowerAtomRegistry.Register(*this);
   this->InitRefParList();
}

ScatteringPowerAtom::ScatteringPowerAtom(const string &name,
                                         const string &symbol,
                                         const REAL bIso):
mpGaussian(0)
{
   VFN_DEBUG_MESSAGE("ScatteringPowerAtom::ScatteringPowerAtom(n,s,B):"<<name,5)
   gScatteringPowerAtomRegistry.Register(*this);
   this->InitRefParList();
   this->Init(name,symbol,bIso);
}

ScatteringPowerAtom::ScatteringPowerAtom(const ScatteringPowerAtom& old):
mpGaussian(0)
{
   VFN_DEBUG_MESSAGE("ScatteringPowerAtom::ScatteringPowerAtom(&old):"<<old.mName,5)
   gScatteringPowerAtomRegistry.Register(*this);
   this->Init(old.GetName(),old.mSymbol,old.mBiso);
   //this->InitRefParList(); //?? :TODO: Check
}

ScatteringPowerAtom::~ScatteringPowerAtom()
{
   VFN_DEBUG_MESSAGE("ScatteringPowerAtom::~ScatteringPowerAtom():"<<mName,5)
   gScatteringPowerAtomRegistry.DeRegister(*this);
   delete mpGaussian;
}

const string& ScatteringPowerAtom::GetClassName() const
{
   const static string className="ScatteringPowerAtom";
   return className;
}

void ScatteringPowerAtom::Init(const string &name,const string &symbol,const REAL bIso)
{
   VFN_DEBUG_MESSAGE("ScatteringPowerAtom::Init(n,s,b)"<<mName,4)
   this->ScatteringPower::Init();
   this->SetName(name);
   mSymbol=symbol;
   mBiso=bIso;
   mIsIsotropic=true;
   if(mpGaussian!=0) delete mpGaussian;
   try
   {
      cctbx::eltbx::xray_scattering::wk1995 wk95t(mSymbol);
      mpGaussian=new cctbx::eltbx::xray_scattering::gaussian(wk95t.fetch());

      this->InitAtNeutronScattCoeffs();

      cctbx::eltbx::tiny_pse::table tpse(mSymbol);
      mAtomicNumber=tpse.atomic_number();

      cctbx::eltbx::icsd_radii::table ticsd(mSymbol);
      mRadius= ticsd.radius();
      cctbx::eltbx::covalent_radii::table tcov(mSymbol);
      mCovalentRadius=tcov.radius();
   }
   catch(cctbx::error)
   {
      cout << "WARNING: could not interpret Symbol name !"<<mSymbol<<endl
           << "         Reverting to H !"<<endl;
      (*fpObjCrystInformUser)("Symbol not understood:"+mSymbol);
      this->Init(name,"H",bIso);
   }
   
   
   VFN_DEBUG_MESSAGE("ScatteringPowerAtom::Init():/Name="<<this->GetName() \
      <<" /Symbol="<<mSymbol<<" /Atomic Number=" << mAtomicNumber,4)
   
   mDynPopCorrIndex=mAtomicNumber;

   //Init default atom colours for POVRay/GUI
   cctbx::eltbx::tiny_pse::table tpse(mSymbol);
   mColourName= tpse.symbol();
   this->InitRGBColour();
   switch(mAtomicNumber)
   {// Values from OpenBabel element.txt
      case 0: mMaxCovBonds=0;break;break;break;//Xx
      case 1: mMaxCovBonds=1;break;break;break;//H
      case 2: mMaxCovBonds=0;break;break;break;//He
      case 3: mMaxCovBonds=1;break;break;break;//Li
      case 4: mMaxCovBonds=2;break;//Be
      case 5: mMaxCovBonds=4;break;//B
      case 6: mMaxCovBonds=4;break;//C
      case 7: mMaxCovBonds=4;break;//N
      case 8: mMaxCovBonds=2;break;//O
      case 9: mMaxCovBonds=1;break;//F
      case 10: mMaxCovBonds=0;break;//Ne
      case 11: mMaxCovBonds=1;break;//Na
      case 12: mMaxCovBonds=2;break;//Mg
      case 13: mMaxCovBonds=6;break;//Al
      case 14: mMaxCovBonds=6;break;//Si
      case 15: mMaxCovBonds=6;break;//P
      case 16: mMaxCovBonds=6;break;//S
      case 17: mMaxCovBonds=1;break;//Cl
      case 18: mMaxCovBonds=0;break;//Ar
      case 19: mMaxCovBonds=1;break;//K
      case 20: mMaxCovBonds=2;break;//Ca
      case 21: mMaxCovBonds=6;break;//Sc
      case 22: mMaxCovBonds=6;break;//Ti
      case 23: mMaxCovBonds=6;break;//V
      case 24: mMaxCovBonds=6;break;//Cr
      case 25: mMaxCovBonds=8;break;//Mn
      case 26: mMaxCovBonds=6;break;//Fe
      case 27: mMaxCovBonds=6;break;//Co
      case 28: mMaxCovBonds=6;break;//Ni
      case 29: mMaxCovBonds=6;break;//Cu
      case 30: mMaxCovBonds=6;break;//Zn
      case 31: mMaxCovBonds=3;break;//Ga
      case 32: mMaxCovBonds=4;break;//Ge
      case 33: mMaxCovBonds=3;break;//As
      case 34: mMaxCovBonds=2;break;//Se
      case 35: mMaxCovBonds=1;break;//Br
      case 36: mMaxCovBonds=0;break;//Kr
      case 37: mMaxCovBonds=1;break;//Rb
      case 38: mMaxCovBonds=2;break;//Sr
      case 39: mMaxCovBonds=6;break;//Y
      case 40: mMaxCovBonds=6;break;//Zr
      case 41: mMaxCovBonds=6;break;//Nb
      case 42: mMaxCovBonds=6;break;//Mo
      case 43: mMaxCovBonds=6;break;//Tc
      case 44: mMaxCovBonds=6;break;//Ru
      case 45: mMaxCovBonds=6;break;//Rh
      case 46: mMaxCovBonds=6;break;//Pd
      case 47: mMaxCovBonds=6;break;//Ag
      case 48: mMaxCovBonds=6;break;//Cd
      case 49: mMaxCovBonds=3;break;//In
      case 50: mMaxCovBonds=4;break;//Sn
      case 51: mMaxCovBonds=3;break;//Sb
      case 52: mMaxCovBonds=2;break;//Te
      case 53: mMaxCovBonds=1;break;//I
      case 54: mMaxCovBonds=0;break;//Xe
      case 55: mMaxCovBonds=1;break;//Cs
      case 56: mMaxCovBonds=2;break;//Ba
      case 57: mMaxCovBonds=12;break;//La
      case 58: mMaxCovBonds=6;break;//Ce
      case 59: mMaxCovBonds=6;break;//Pr
      case 60: mMaxCovBonds=6;break;//Nd
      case 61: mMaxCovBonds=6;break;//Pm
      case 62: mMaxCovBonds=6;break;//Sm
      case 63: mMaxCovBonds=6;break;//Eu
      case 64: mMaxCovBonds=6;break;//Gd
      case 65: mMaxCovBonds=6;break;//Tb
      case 66: mMaxCovBonds=6;break;//Dy
      case 67: mMaxCovBonds=6;break;//Ho
      case 68: mMaxCovBonds=6;break;//Er
      case 69: mMaxCovBonds=6;break;//Tm
      case 70: mMaxCovBonds=6;break;//Yb
      case 71: mMaxCovBonds=6;break;//Lu
      case 72: mMaxCovBonds=6;break;//Hf
      case 73: mMaxCovBonds=6;break;//Ta
      case 74: mMaxCovBonds=6;break;//W
      case 75: mMaxCovBonds=6;break;//Re
      case 76: mMaxCovBonds=6;break;//Os
      case 77: mMaxCovBonds=6;break;//Ir
      case 78: mMaxCovBonds=6;break;//Pt
      case 79: mMaxCovBonds=6;break;//Au
      case 80: mMaxCovBonds=6;break;//Hg
      case 81: mMaxCovBonds=3;break;//Tl
      case 82: mMaxCovBonds=4;break;//Pb
      case 83: mMaxCovBonds=3;break;//Bi
      case 84: mMaxCovBonds=2;break;//Po
      case 85: mMaxCovBonds=1;break;//At
      case 86: mMaxCovBonds=0;break;//Rn
      case 87: mMaxCovBonds=1;break;//Fr
      case 88: mMaxCovBonds=2;break;//Ra
      default: mMaxCovBonds=6;break;
   }

   VFN_DEBUG_MESSAGE("ScatteringPowerAtom::Init(n,s,b):End",3)
}

CrystVector_REAL ScatteringPowerAtom::GetScatteringFactor(const ScatteringData &data,
                                                            const int spgSymPosIndex) const
{
   VFN_DEBUG_MESSAGE("ScatteringPower::GetScatteringFactor(&data):"<<mName,3)
   CrystVector_REAL sf(data.GetNbRefl());
   switch(data.GetRadiationType())
   {
      case(RAD_NEUTRON):
      {
         VFN_DEBUG_MESSAGE("ScatteringPower::GetScatteringFactor():NEUTRON:"<<mName,3)
         sf=mNeutronScattLengthReal;
         break;
      }
      case(RAD_XRAY):
      {
         VFN_DEBUG_MESSAGE("ScatteringPower::GetScatteringFactor():XRAY:"<<mName,3)
         if(mpGaussian!=0)
         {
            const long nb=data.GetSinThetaOverLambda().numElements();
            const REAL *pstol=data.GetSinThetaOverLambda().data();
            for(long i=0;i<nb;i++)
               sf(i)=mpGaussian->at_stol(*pstol++);
         }
         else sf=1.0;//:KLUDGE:  Should never happen
         break;
      }
      case(RAD_ELECTRON):
      {
         VFN_DEBUG_MESSAGE("ScatteringPower::GetScatteringFactor():ELECTRON:"<<mName,3)
         if(mpGaussian!=0)
         {
            const REAL z=this->GetAtomicNumber();
            const long nb=data.GetSinThetaOverLambda().numElements();
            const REAL *pstol=data.GetSinThetaOverLambda().data();
            for(long i=0;i<nb;i++)
               sf(i)=(z-mpGaussian->at_stol(*pstol))/(*pstol * *pstol);
               pstol++;
         }
         else sf=1.0;//:KLUDGE:  Should never happen
         break;
      }
   }
   VFN_DEBUG_MESSAGE("ScatteringPower::GetScatteringFactor(&data):End",3)
   return sf;
}

REAL ScatteringPowerAtom::GetForwardScatteringFactor(const RadiationType type) const
{
   REAL sf;
   switch(type)
   {
      case(RAD_NEUTRON):
      {
         sf=mNeutronScattLengthReal;
         break;
      }
      case(RAD_XRAY):
      {
         if(mpGaussian!=0)
            sf=mpGaussian->at_stol(0);
         else sf=1.0;
         break;
      }
      case(RAD_ELECTRON):
      {
         const REAL z=this->GetAtomicNumber();
         sf=(z-mpGaussian->at_stol(0.0001))/(.0001 * .0001);
      }
   }
   VFN_DEBUG_MESSAGE("ScatteringPower::GetScatteringFactor(&data):End",3)
   return sf;
}

static bool warnADP=true;
CrystVector_REAL ScatteringPowerAtom::GetTemperatureFactor(const ScatteringData &data,
                                                             const int spgSymPosIndex) const
{
   VFN_DEBUG_MESSAGE("ScatteringPower::GetTemperatureFactor(&data):"<<mName,3)
   CrystVector_REAL sf(data.GetNbRefl());
   if((mIsIsotropic==false) && warnADP)
   {  // Warn once
      cout<<"========================== WARNING ========================="<<endl
          <<"   In ScatteringPowerAtom::GetTemperatureFactor():"<<endl
          <<"   Anisotropic Displacement Parameters are not currently properly handled"<<endl
          <<"   for Debye-Waller calculations (no symmetry handling for ADPs)."<<endl
          <<"   =>The Debye-Waller calculations will instead use only isotropic DPs"<<endl<<endl;
      warnADP=false;
   }

   if(true)//(mIsIsotropic)
   {
      CrystVector_REAL stolsq(data.GetNbRefl());
      const CrystVector_REAL stol=data.GetSinThetaOverLambda();
      stolsq=stol;
      stolsq*=stol;
      
      #ifdef __VFN_VECTOR_USE_BLITZ__
         #define SF sf
         #define STOLSQ stolsq
      #else
         #define SF (*ssf)
         #define STOLSQ (*sstolsq)

         REAL *ssf=sf.data();
         const REAL *sstolsq=stolsq.data();

         for(long ii=0;ii<sf.numElements();ii++)
         {
      #endif
      
         SF=exp(-mBiso*STOLSQ);
         
      #ifdef __VFN_VECTOR_USE_BLITZ__

      #else
         ssf++;
         sstolsq++;
         }
      #endif

      #undef SF
      #undef STOLSQ
   }
   else
   {// :TODO: handle ADP - requires taking into account symmetries... 
      const REAL b11=mBeta(0);
      const REAL b22=mBeta(1);
      const REAL b33=mBeta(2);
      const REAL b12=mBeta(3);
      const REAL b13=mBeta(4);
      const REAL b23=mBeta(5);
      #ifdef __VFN_VECTOR_USE_BLITZ__
         #define HH data.H()
         #define KK data.K()
         #define LL data.L()
         #define SF sf
      #else
         #define HH (*hh)
         #define KK (*kk)
         #define LL (*ll)
         #define SF (*ssf)

         const REAL *hh=(data.GetH()).data();
         const REAL *kk=(data.GetK()).data();
         const REAL *ll=(data.GetL()).data();
         REAL *ssf=sf.data();

         for(long ii=0;ii<sf.numElements();ii++)
         {
      #endif
   
      SF=   exp( -b11*pow(HH,2)
                 -b22*pow(KK,2)
                 -b33*pow(LL,2)
                 -2*b12*HH*KK
                 -2*b13*HH*LL
                 -2*b23*KK*LL);
                 
      #ifdef __VFN_VECTOR_USE_BLITZ__

      #else
         hh++;
         kk++;
         ll++;
         ssf++;
         }
      #endif

      #undef HH
      #undef KK
      #undef LL
      #undef SF
   }
   return sf;
}

CrystMatrix_REAL ScatteringPowerAtom::
   GetResonantScattFactReal(const ScatteringData &data,
                            const int spgSymPosIndex) const
{
   VFN_DEBUG_MESSAGE("ScatteringPower::GetResonantScattFactReal(&data):"<<mName,3)
   CrystMatrix_REAL fprime(1,1);//:TODO: More than one wavelength
   CrystMatrix_REAL fsecond(1,1);
   switch(data.GetRadiationType())
   {
      case(RAD_NEUTRON):
      {
         fprime=0;
         fsecond=mNeutronScattLengthImag;
         break;
      }
      case(RAD_XRAY):
      {
         try
         {
            cctbx::eltbx::henke::table thenke(mSymbol);
            cctbx::eltbx::fp_fdp f=thenke.at_angstrom(data.GetWavelength()(0));

            if(f.is_valid_fp()) fprime(0)=f.fp();
            else fprime(0)=0;
            if(f.is_valid_fdp()) fsecond(0)=f.fdp();
            else fsecond(0)=0;
         }
         catch(cctbx::error)
         {
            fprime(0)=0;
            fsecond(0)=0;
         }
         break;
      }
      case(RAD_ELECTRON):
      {
         fprime=0;
         fsecond=0;
         break;
      }
   }
   return fprime;
}

CrystMatrix_REAL ScatteringPowerAtom::
   GetResonantScattFactImag(const ScatteringData &data,
                            const int spgSymPosIndex) const
{
   VFN_DEBUG_MESSAGE("ScatteringPower::GetResonantScattFactImag():"<<mName,3)
   CrystMatrix_REAL fprime(1,1);//:TODO: More than one wavelength
   CrystMatrix_REAL fsecond(1,1);
   switch(data.GetRadiationType())
   {
      case(RAD_NEUTRON):
      {
         fprime=0;
         fsecond=mNeutronScattLengthImag;
         break;
      }
      case(RAD_XRAY):
      {
         try
         {
            cctbx::eltbx::henke::table thenke(mSymbol);
            cctbx::eltbx::fp_fdp f=thenke.at_angstrom(data.GetWavelength()(0));

            if(f.is_valid_fp()) fprime(0)=f.fp();
            else fprime(0)=0;
            if(f.is_valid_fdp()) fsecond(0)=f.fdp();
            else fsecond(0)=0;
         }
         catch(cctbx::error)
         {
            fprime(0)=0;
            fsecond(0)=0;
         }
         break;
      }
      case(RAD_ELECTRON):
      {
         fprime=0;
         fsecond=0;
         break;
      }
   }
   return fsecond;
}


void ScatteringPowerAtom::SetSymbol(const string& symbol)
{
   VFN_DEBUG_MESSAGE("ScatteringPowerAtom::SetSymbol():"<<mName,5)
   this->Init(this->GetName(),symbol,this->GetBiso());
}
const string& ScatteringPowerAtom::GetSymbol() const
{
   VFN_DEBUG_MESSAGE("ScatteringPowerAtom::GetSymbol():"<<mName,5)
   return mSymbol;
}

string ScatteringPowerAtom::GetElementName() const
{
   VFN_DEBUG_MESSAGE("ScatteringPowerAtom::GetElementName():"<<mName,2)
   try
   {
      cctbx::eltbx::tiny_pse::table tpse(mSymbol);
      return tpse.name();
   }
   catch(cctbx::error)
   {
      cout << "WARNING: could not interpret Symbol:"<<mSymbol<<endl;
   }
   return "Unknown";
}

int ScatteringPowerAtom::GetAtomicNumber() const {return mAtomicNumber;}
REAL ScatteringPowerAtom::GetRadius() const {return mRadius;}
REAL ScatteringPowerAtom::GetCovalentRadius() const {return mCovalentRadius;}
unsigned int ScatteringPowerAtom::GetMaxCovBonds()const{ return mMaxCovBonds;}

void ScatteringPowerAtom::Print()const
{
   VFN_DEBUG_MESSAGE("ScatteringPowerAtom::Print()",1)
   cout << "ScatteringPowerAtom ("<<this->GetName()<<","
        << FormatString(this->GetSymbol(),4) << ") :"
        << FormatFloat(this->GetBiso());
   VFN_DEBUG_MESSAGE_SHORT("at "<<this,10)
   cout << endl;
}

void ScatteringPowerAtom::InitAtNeutronScattCoeffs()
{
   VFN_DEBUG_MESSAGE("ScatteringPowerAtom::InitAtNeutronScattCoeffs():"<<mName,3)
   mClock.Click();
   try
   {
      cctbx::eltbx::neutron::neutron_news_1992_table nn92t(mSymbol);
      mNeutronScattLengthReal=nn92t.bound_coh_scatt_length_real();
      mNeutronScattLengthImag=nn92t.bound_coh_scatt_length_imag();
   }
   catch(cctbx::error)
   {
      cout << "WARNING: could not interpret symbol for neutron coeefs:"<<mSymbol<<endl;
   }
   
   VFN_DEBUG_MESSAGE("ScatteringPowerAtom::InitAtNeutronScattCoeffs():End",3)
}

void ScatteringPowerAtom::InitRefParList()
{
   VFN_DEBUG_MESSAGE("ScatteringPowerAtom::InitRefParList():"<<mName,5)
   {
      RefinablePar tmp("Biso",&mBiso,0.1,5.,
                        gpRefParTypeScattPowTemperatureIso,REFPAR_DERIV_STEP_ABSOLUTE,
                        true,true,true,false);
      tmp.SetDerivStep(1e-3);
      tmp.SetGlobalOptimStep(.5);
      tmp.AssignClock(mClock);
      this->AddPar(tmp);
   }
   {
      REAL* bdata = (REAL*) mB.data();

      RefinablePar B11("B11",&bdata[0],0.1,5.,
              gpRefParTypeScattPowTemperatureAniso,REFPAR_DERIV_STEP_ABSOLUTE,
              true,true,false,false);
      B11.SetDerivStep(1e-3);
      B11.SetGlobalOptimStep(.5);
      B11.AssignClock(mClock);
      this->AddPar(B11);

      RefinablePar B22("B22",&bdata[1],0.1,5.,
              gpRefParTypeScattPowTemperatureAniso,REFPAR_DERIV_STEP_ABSOLUTE,
              true,true,false,false);
      B22.SetDerivStep(1e-3);
      B22.SetGlobalOptimStep(.5);
      B22.AssignClock(mClock);
      this->AddPar(B22);

      RefinablePar B33("B33",&bdata[2],0.1,5.,
              gpRefParTypeScattPowTemperatureAniso,REFPAR_DERIV_STEP_ABSOLUTE,
              true,true,false,false);
      B33.SetDerivStep(1e-3);
      B33.SetGlobalOptimStep(.5);
      B33.AssignClock(mClock);
      this->AddPar(B33);

      RefinablePar B12("B12",&bdata[3],-5,5.,
              gpRefParTypeScattPowTemperatureAniso,REFPAR_DERIV_STEP_ABSOLUTE,
              true,true,false,false);
      B12.SetDerivStep(1e-3);
      B12.SetGlobalOptimStep(.5);
      B12.AssignClock(mClock);
      this->AddPar(B12);

      RefinablePar B13("B13",&bdata[4],-5,5.,
              gpRefParTypeScattPowTemperatureAniso,REFPAR_DERIV_STEP_ABSOLUTE,
              true,true,false,false);
      B13.SetDerivStep(1e-3);
      B13.SetGlobalOptimStep(.5);
      B13.AssignClock(mClock);
      this->AddPar(B13);

      RefinablePar B23("B23",&bdata[5],-5,5.,
              gpRefParTypeScattPowTemperatureAniso,REFPAR_DERIV_STEP_ABSOLUTE,
              true,true,false,false);
      B23.SetDerivStep(1e-3);
      B23.SetGlobalOptimStep(.5);
      B23.AssignClock(mClock);
      this->AddPar(B23);
   }
   {
      RefinablePar tmp("ML Error",&mMaximumLikelihoodPositionError,0.,1.,
                        gpRefParTypeScattPow,REFPAR_DERIV_STEP_ABSOLUTE,
                        false,true,true,false);
      tmp.SetDerivStep(1e-4);
      tmp.SetGlobalOptimStep(.001);
      tmp.AssignClock(mMaximumLikelihoodParClock);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("ML-Nb Ghost Atoms",&mMaximumLikelihoodNbGhost,0.,10.,
                        gpRefParTypeScattPow,REFPAR_DERIV_STEP_ABSOLUTE,
                        true,true,true,false);
      tmp.SetDerivStep(1e-3);
      tmp.SetGlobalOptimStep(.05);
      tmp.AssignClock(mMaximumLikelihoodParClock);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("Formal Charge",&mFormalCharge,-10.,10.,
                        gpRefParTypeScattPow,REFPAR_DERIV_STEP_ABSOLUTE,
                        true,true,true,false);
      tmp.SetDerivStep(1e-3);
      tmp.SetGlobalOptimStep(.05);
      tmp.AssignClock(mClock);
      this->AddPar(tmp);
   }
}
#ifdef __WX__CRYST__
WXCrystObjBasic* ScatteringPowerAtom::WXCreate(wxWindow* parent)
{
   //:TODO: Check mpWXCrystObj==0
   mpWXCrystObj=new WXScatteringPowerAtom(parent,this);
   return mpWXCrystObj;
}
#endif

//######################################################################
//
//      SCATTERING COMPONENT
//
//######################################################################
ScatteringComponent::ScatteringComponent():
mX(0),mY(0),mZ(0),mOccupancy(0),mpScattPow(0),mDynPopCorr(0)
{}
bool ScatteringComponent::operator==(const ScatteringComponent& rhs) const
{
   return ((mX==rhs.mX) && (mY==rhs.mY) && (mZ==rhs.mZ) &&
           (mOccupancy==rhs.mOccupancy) && (mpScattPow==rhs.mpScattPow));
}
bool ScatteringComponent::operator!=(const ScatteringComponent& rhs) const
{
   return ((mX!=rhs.mX) || (mY!=rhs.mY) || (mZ!=rhs.mZ) ||
           (mOccupancy!=rhs.mOccupancy) || (mpScattPow!=rhs.mpScattPow));
}
void ScatteringComponent::Print()const
{
   cout <<mX<<" "<<mY<<" "<<mZ<<" "<<mOccupancy<<" "<<mDynPopCorr<<" "<<mpScattPow;
   if(0!=mpScattPow) cout <<" "<<mpScattPow->GetName();
   cout<<endl;
}

//######################################################################
//
//      SCATTERING COMPONENT LIST
//
//######################################################################

ScatteringComponentList::ScatteringComponentList()
{
}

ScatteringComponentList::ScatteringComponentList(const long nbComponent)
{
   mvScattComp.resize(nbComponent);
}

ScatteringComponentList::ScatteringComponentList(const ScatteringComponentList &old):
mvScattComp(old.mvScattComp)
{
}

ScatteringComponentList::~ScatteringComponentList()
{
}
void ScatteringComponentList::Reset()
{
   mvScattComp.clear();
}

const ScatteringComponent& ScatteringComponentList::operator()(const long i) const
{
   VFN_DEBUG_MESSAGE("ScatteringComponentList::operator()("<<i<<")",1)
   if(i>=this->GetNbComponent())
   {
      this->Print();
      throw ObjCrystException("ScatteringComponentList::operator()(i)::i>mNbComponent!!");
   }
   if(i<0) throw ObjCrystException("ScatteringComponentList::operator()&(i)::i<0!!");
   return mvScattComp[i];
}

ScatteringComponent& ScatteringComponentList::operator()(const long i)
{
   VFN_DEBUG_MESSAGE("ScatteringComponentList::operator()&("<<i<<")",1)
   if(i>=this->GetNbComponent())
   {
      this->Print();
      throw ObjCrystException("ScatteringComponentList::operator()&(i)::i>mNbComponent!!");
   }
   if(i<0) throw ObjCrystException("ScatteringComponentList::operator()&(i):: i<0!!");
   return mvScattComp[i];
}

long ScatteringComponentList::GetNbComponent() const {return mvScattComp.size();}

void ScatteringComponentList::operator=(const ScatteringComponentList &rhs)
{
   VFN_DEBUG_MESSAGE("ScatteringComponentList::operator=()",1)
   mvScattComp=rhs.mvScattComp;
   VFN_DEBUG_MESSAGE("ScatteringComponentList::operator=():End",0)
}
bool ScatteringComponentList::operator==(const ScatteringComponentList &rhs)const
{
   if(rhs.GetNbComponent() != this->GetNbComponent()) return false;
   for(long i=0;i<this->GetNbComponent();i++)
      if( (*this)(i) != rhs(i) ) return false;
   return true;
}

void ScatteringComponentList::operator+=(const ScatteringComponentList &rhs)
{
   for(long i=0;i<rhs.GetNbComponent();i++) 
      mvScattComp.push_back(rhs(i));
}
void ScatteringComponentList::operator+=(const ScatteringComponent &rhs)
{
   VFN_DEBUG_MESSAGE("ScatteringComponentList::operator+=()",1)
   mvScattComp.push_back(rhs);
}

void ScatteringComponentList::operator++()
{
   VFN_DEBUG_MESSAGE("ScatteringComponentList::operator++()",1)
   mvScattComp.resize(this->GetNbComponent()+1);
}

void ScatteringComponentList::operator--()
{
   VFN_DEBUG_MESSAGE("ScatteringComponentList::operator--()",1)
   if(this->GetNbComponent()>0)
      mvScattComp.resize(this->GetNbComponent()-1);
}

void ScatteringComponentList::Print()const
{
   VFN_DEBUG_ENTRY("ScatteringComponentList::Print()",5)
   cout<<"Number of Scattering components:"<<this->GetNbComponent()<<endl;
   for(long i=0;i<this->GetNbComponent();i++)
   {
      cout << i<<":";
      (*this)(i).Print();
   }
   VFN_DEBUG_EXIT("ScatteringComponentList::Print()",5)
}

}//namespace
