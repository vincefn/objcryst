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

#include <stdio.h> //for sprintf()

#include "ObjCryst/ScatteringPower.h"
#include "Quirks/VFNStreamFormat.h"
#include "Quirks/VFNDebug.h"

#include "ObjCryst/Colours.h"

#ifdef __WX__CRYST__
   #include "wxCryst/wxScatteringPower.h"
#endif

#include <iomanip>
#include <fstream>

namespace ObjCryst
{

////////////////////////////////////////////////////////////////////////
//
// Using external functions from 'atominfo' package
//
//-> Coefficients for the analytical approximation of scattering factors.
//AtomInfo (c) 1994-96 Ralf W. Grosse-Kunstleve 
//
////////////////////////////////////////////////////////////////////////

extern "C"
{
#include "atominfo/atominfo.h"
}

const RefParType *gpRefParTypeScattPow
   =new RefParType(gpRefParTypeObjCryst,"Scattering Power");
const RefParType *gpRefParTypeScattPowResonant
   =new RefParType(gpRefParTypeScattPow,"Resonant Scatt.");
const RefParType *gpRefParTypeScattPowTemperature
   =new RefParType(gpRefParTypeScattPow,"Temperature");
const RefParType *gpRefParTypeScattPowTemperatureIso
   =new RefParType(gpRefParTypeScattPowTemperature,"Isotropic");
const RefParType *gpRefParTypeScattPowTemperatureAniso
   =new RefParType(gpRefParTypeScattPowTemperatureIso,"Anisotropic");

//######################################################################
//
//      SCATTERING POWER
//
//######################################################################
ObjRegistry<ScatteringPower> gScatteringPowerRegistry("Global ScatteringPower Registry");

long ScatteringPower::mNbScatteringPower=0;
ScatteringPower* ScatteringPower::mspScatteringPowerGlobalList[1000];
CrystVector_bool ScatteringPower::mspScatteringPowerGlobalListIsUsed(1000);


ScatteringPower::ScatteringPower():mDynPopCorrIndex(0),mBiso(1.0),mIsIsotropic(true),
mScatteringPowerId(mNbScatteringPower)
{
   VFN_DEBUG_MESSAGE("ScatteringPower::ScatteringPower():"<<mName,5)
   if(mNbScatteringPower>1000) throw ObjCrystException("ScatteringPower::ScatteringPower() \
      Reached maximum number of ScatteringPower objects (1000!)");
   mspScatteringPowerGlobalList[mNbScatteringPower] = this;
   if(0==mNbScatteringPower) mspScatteringPowerGlobalListIsUsed=false;
   mspScatteringPowerGlobalListIsUsed(mNbScatteringPower)=true;
   mNbScatteringPower++;
   gScatteringPowerRegistry.Register(*this);
   this->Init();
}
ScatteringPower::ScatteringPower(const ScatteringPower& old):
mDynPopCorrIndex(old.mDynPopCorrIndex),mBiso(old.mBiso),mIsIsotropic(old.mIsIsotropic),
mBeta(old.mBeta),mScatteringPowerId(mNbScatteringPower)
{
   VFN_DEBUG_MESSAGE("ScatteringPower::ScatteringPower(&old):"<<mName,5)
   if(mNbScatteringPower>1000) throw ObjCrystException("ScatteringPower::ScatteringPower() \
      Reached maximum number of ScatteringPower objects (1000!)");
   mspScatteringPowerGlobalList[mNbScatteringPower] = this;
   mspScatteringPowerGlobalListIsUsed(mNbScatteringPower)=true;
   mNbScatteringPower++;
   gScatteringPowerRegistry.Register(*this);
   this->Init();
}
ScatteringPower::~ScatteringPower()
{
   VFN_DEBUG_MESSAGE("ScatteringPower::~ScatteringPower():"<<mName,5)
   mspScatteringPowerGlobalListIsUsed(mScatteringPowerId)=false;
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
}

bool ScatteringPower::IsScatteringFactorAnisotropic()const{return false;}
bool ScatteringPower::IsTemperatureFactorAnisotropic()const{return false;}
bool ScatteringPower::IsResonantScatteringAnisotropic()const{return false;}

const string& ScatteringPower::GetSymbol() const {return this->GetName();}
REAL ScatteringPower::GetBiso() const {return mBiso;}
REAL& ScatteringPower::GetBiso() {mClock.Click();return mBiso;}
void ScatteringPower::SetBiso(const REAL newB) { mClock.Click();mBiso=newB;}
bool ScatteringPower::IsIsotropic() const {return mIsIsotropic;}
long ScatteringPower::GetDynPopCorrIndex() const {return mDynPopCorrIndex;}
long ScatteringPower::GetScatteringPowerId()const {return mScatteringPowerId;}
long ScatteringPower::GetNbScatteringPower()const {return mNbScatteringPower;}
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

void ScatteringPower::Init()
{
   VFN_DEBUG_MESSAGE("ScatteringPower::Init():"<<mName,2)
   mColourName="White";
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
      if(gPOVRayColours[i].mName=="")
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
// GetScatteringPower
//######################################################################
const ScatteringPower& GetScatteringPower(const long i) 
{  
   VFN_DEBUG_MESSAGE("GetScatteringPower(i):"<<i,2)
   if(false==ScatteringPower::mspScatteringPowerGlobalListIsUsed(i))
   {
      cout <<i<<endl;
      cout <<ScatteringPower::mspScatteringPowerGlobalListIsUsed<<endl;
      throw ObjCrystException("GetScatteringPower(i):this Scattering power does not exist !!");
   }

   return *ScatteringPower::mspScatteringPowerGlobalList[i];
}


//######################################################################
//
//      SCATTERING POWER ATOM
//
//######################################################################
ObjRegistry<ScatteringPowerAtom> 
   gScatteringPowerAtomRegistry("Global ScatteringPowerAtom Registry");

ScatteringPowerAtom::ScatteringPowerAtom():
ScatteringPower(),mSymbol(""),mAtomicNumber(0),mScattAi(5),mScattBi(5)
{
   VFN_DEBUG_MESSAGE("ScatteringPowerAtom::ScatteringPowerAtom():"<<mName,5)
   gScatteringPowerAtomRegistry.Register(*this);
   this->InitRefParList();
}

ScatteringPowerAtom::ScatteringPowerAtom(const string &name,
                                         const string &symbol,
                                         const REAL bIso):
mScattAi(5),mScattBi(5)
{
   VFN_DEBUG_MESSAGE("ScatteringPowerAtom::ScatteringPowerAtom(n,s,B):"<<name,5)
   gScatteringPowerAtomRegistry.Register(*this);
   this->InitRefParList();
   this->Init(name,symbol,bIso);
}

ScatteringPowerAtom::ScatteringPowerAtom(const ScatteringPowerAtom& old):
mScattAi(5),mScattBi(5)
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
}

const string& ScatteringPowerAtom::GetClassName() const
{
	const static string className="ScatteringPowerAtom";
	return className;
}

void ScatteringPowerAtom::Init(const string &name,const string &symbol,const REAL bIso)
{
   VFN_DEBUG_MESSAGE("ScatteringPowerAtom::Init(n,s,b)"<<mName<<":"<<mScatteringPowerId,4)
   this->ScatteringPower::Init();
   this->SetName(name);
   mSymbol=symbol;
   mBiso=bIso;
   mIsIsotropic=true;
   this->InitAtScattCoeffsWK95();
   this->InitAtNeutronScattCoeffs();
   
   const T_PSE * tmp;
   tmp=FindInPSE(mSymbol.c_str(),0);
   mAtomicNumber=tmp->Z;
   //delete tmp;
   
   const T_AtomRadius *atomRadius;
   atomRadius=FindAtomRadius(mSymbol.c_str(),0);
   mRadius= atomRadius->Radius;
   //delete atomRadius;
   
   VFN_DEBUG_MESSAGE("ScatteringPowerAtom::Init():/Name="<<this->GetName() \
      <<" /Symbol="<<mSymbol<<" /Atomic Number=" << mAtomicNumber,4)
   
   mDynPopCorrIndex=mAtomicNumber;

   //Init default atom colours for POVRay
   switch(mAtomicNumber)
   {
      case   1: mColourName="White";           break;   //hydrogen    
      case   2: mColourName="White";           break;   //helium      
      case   3: mColourName="White";           break;   //lithium     
      case   4: mColourName="White";           break;   //beryllium   
      case   5: mColourName="White";           break;   //boron       
      case   6: mColourName="Gray50";          break;   //carbon      
      case   7: mColourName="DarkGreen";           break;   //nitrogen    
      case   8: mColourName="Red";             break;   //oxygen      
      case   9: mColourName="White";           break;   //fluorine    
      case  10: mColourName="White";           break;   //neon        
      case  11: mColourName="White";           break;   //sodium      
      case  12: mColourName="White";           break;   //magnesium   
      case  13: mColourName="White";           break;   //aluminium   
      case  14: mColourName="White";           break;   //silicon     
      case  15: mColourName="White";           break;   //phosphorus  
      case  16: mColourName="Yellow";          break;   //sulphur     
      case  17: mColourName="SeaGreen";        break;   //chlorine    
      case  18: mColourName="White";           break;   //argon       
      case  19: mColourName="White";           break;   //potassium   
      case  20: mColourName="White";           break;   //calcium     
      case  21: mColourName="White";           break;   //scandium    
      case  22: mColourName="White";           break;   //titanium    
      case  23: mColourName="White";           break;   //vanadium    
      case  24: mColourName="Gray20";          break;   //chromium    
      case  25: mColourName="MediumWood";      break;   //manganese   
      case  26: mColourName="Orange";          break;   //iron        
      case  27: mColourName="Pink";            break;   //cobalt      
      case  28: mColourName="Green";           break;   //nickel      
      case  29: mColourName="Copper";          break;   //copper      
      case  30: mColourName="Gray80";          break;   //zinc        
      case  31: mColourName="White";           break;   //gallium     
      case  32: mColourName="White";           break;   //germanium   
      case  33: mColourName="White";           break;   //arsenic     
      case  34: mColourName="Flesh";           break;   //selenium    
      case  35: mColourName="White";           break;   //bromine     
      case  36: mColourName="White";           break;   //krypton     
      case  37: mColourName="White";           break;   //rubidium    
      case  38: mColourName="White";           break;   //strontium   
      case  39: mColourName="White";           break;   //yttrium     
      case  40: mColourName="White";           break;   //zirconium   
      case  41: mColourName="White";           break;   //niobium     
      case  42: mColourName="White";           break;   //molybdenum  
      case  43: mColourName="White";           break;   //technetium  
      case  44: mColourName="White";           break;   //ruthenium   
      case  45: mColourName="White";           break;   //rhodium     
      case  46: mColourName="White";           break;   //palladium   
      case  47: mColourName="Silver";          break;   //silver      
      case  48: mColourName="White";           break;   //cadmium     
      case  49: mColourName="White";           break;   //indium      
      case  50: mColourName="White";           break;   //tin         
      case  51: mColourName="White";           break;   //antimony    
      case  52: mColourName="White";           break;   //tellurium   
      case  53: mColourName="OrangeRed";       break;   //iodine      
      case  54: mColourName="White";           break;   //xenon       
      case  55: mColourName="White";           break;   //caesium     
      case  56: mColourName="White";           break;   //barium      
      case  57: mColourName="White";           break;   //lanthanum   
      case  58: mColourName="White";           break;   //cerium      
      case  59: mColourName="White";           break;   //praseodymium
      case  60: mColourName="CornflowerBlue";  break;   //neodymium   
      case  61: mColourName="White";           break;   //promethium  
      case  62: mColourName="White";           break;   //samarium    
      case  63: mColourName="White";           break;   //europium    
      case  64: mColourName="White";           break;   //gadolinium  
      case  65: mColourName="White";           break;   //terbium     
      case  66: mColourName="White";           break;   //dysprosium  
      case  67: mColourName="White";           break;   //holmium     
      case  68: mColourName="SkyBlue";           break;   //erbium      
      case  69: mColourName="White";           break;   //thulium     
      case  70: mColourName="White";           break;   //ytterbium   
      case  71: mColourName="White";           break;   //lutetium    
      case  72: mColourName="White";           break;   //hafnium     
      case  73: mColourName="Blue";            break;   //tantalum    
      case  74: mColourName="White";           break;   //tungsten    
      case  75: mColourName="White";           break;   //rhenium     
      case  76: mColourName="White";           break;   //osmium      
      case  77: mColourName="ForestGreen";     break;   //iridium     
      case  78: mColourName="White";           break;   //platinum    
      case  79: mColourName="Gold";            break;   //gold        
      case  80: mColourName="VioletRed";       break;   //mercury     
      case  81: mColourName="White";           break;   //thallium    
      case  82: mColourName="White";           break;   //lead        
      case  83: mColourName="White";           break;   //bismuth     
      case  84: mColourName="White";           break;   //polonium    
      case  85: mColourName="White";           break;   //astatine    
      case  86: mColourName="White";           break;   //radon       
      case  87: mColourName="White";           break;   //francium    
      case  88: mColourName="White";           break;   //radium      
      case  89: mColourName="White";           break;   //actinium    
      case  90: mColourName="White";           break;   //thorium     
      case  91: mColourName="White";           break;   //protactinium
      case  92: mColourName="White";           break;   //uranium     
      case  93: mColourName="White";           break;   //neptunium   
      case  94: mColourName="White";           break;   //plutonium   
      case  95: mColourName="White";           break;   //americium   
      case  96: mColourName="White";           break;   //curium      
      case  97: mColourName="White";           break;   //berkelium   
      case  98: mColourName="White";           break;   //californium 
      case  99: mColourName="White";           break;   //einsteinium 
      case 100: mColourName="White";           break;   //fermium     
      case 101: mColourName="White";           break;   //mendelevium 
      case 102: mColourName="White";           break;   //nobelium    
      case 103: mColourName="White";           break;   //lawrencium  
      default : mColourName="White";           break;   
   }
   this->InitRGBColour();

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
         CrystVector_REAL stolsq(data.GetNbRefl());
         stolsq=data.GetSinThetaOverLambda();
         stolsq*=data.GetSinThetaOverLambda();
         sf=mScattC;
         REAL a,b;
         for(int i=0;i<5;i++) 
         {
            a=  mScattAi(i);
            b= -mScattBi(i);
            
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
            
               SF+=a*exp(STOLSQ*b);
               
            #ifdef __VFN_VECTOR_USE_BLITZ__

            #else
               ssf++;
               sstolsq++;
               }
            #endif

            #undef SF
            #undef STOLSQ
         }
         break;
      }
      case(RAD_ELECTRON):
      {
         throw ObjCrystException("ScatteringPowerAtom::GetScatteringFactor(data): Scattering factors not implemented for electrons yet !!!:");
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
         sf=  mScattAi.sum();
         break;
      }
      case(RAD_ELECTRON):
      {
         throw ObjCrystException("ScatteringPowerAtom::GetForwardScatteringFactor(data): Scattering factors not implemented for electrons yet !!!:");
      }
   }
   VFN_DEBUG_MESSAGE("ScatteringPower::GetScatteringFactor(&data):End",3)
   return sf;
}

CrystVector_REAL ScatteringPowerAtom::GetTemperatureFactor(const ScatteringData &data,
                                                             const int spgSymPosIndex) const
{
   VFN_DEBUG_MESSAGE("ScatteringPower::GetTemperatureFactor(&data):"<<mName,3)
   CrystVector_REAL sf(data.GetNbRefl());
   if(mIsIsotropic)
   {
      // :NOTE: can't use 'return exp(-mBiso*pow2(diffData.GetSinThetaOverLambda()))'
      //using kcc (OK with gcc)
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
   {
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
         const char *tmpSymbol=mSymbol.c_str();
         float fp,fs;
         Get_fpfs_Henke(tmpSymbol, XRAY_WAVELENGTH_TO_ENERGY/data.GetWavelength()(0),
                           &fp, &fs);
         fprime(0)=fp;
         fsecond(0)=fs;
         //do a check : any value missing ?
         if( (fprime(0)<-9000) || (fsecond(0)<-9000) )
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
         const char *tmpSymbol=mSymbol.c_str();
         float fp,fs;
         Get_fpfs_Henke(tmpSymbol, XRAY_WAVELENGTH_TO_ENERGY/data.GetWavelength()(0),
                           &fp, &fs);
         fprime(0)=fp;
         fsecond(0)=fs;
         //do a check : any value missing ?
         if( (fprime(0)<-9000) || (fsecond(0)<-9000) )
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
   const T_PSE * atomMass;
   atomMass=FindInPSE(mSymbol.c_str(),0);
   return atomMass->Name;
}

int ScatteringPowerAtom::GetAtomicNumber() const {return mAtomicNumber;}
REAL ScatteringPowerAtom::GetRadius() const {return mRadius;}

void ScatteringPowerAtom::Print()const
{
   VFN_DEBUG_MESSAGE("ScatteringPowerAtom::Print()",1)
   cout << "ScatteringPowerAtom ("<<this->GetName()<<","
        << FormatString(this->GetSymbol(),4) << ") :"
        << FormatFloat(this->GetBiso());
   VFN_DEBUG_MESSAGE_SHORT("at "<<this,10)
   cout << endl;
}

void ScatteringPowerAtom::InitAtScattCoeffsWK95()
{
   VFN_DEBUG_MESSAGE("ScatteringPowerAtom::InitAtScattCoeffsWK95():"<<mName,3)
   mClock.Click();
   const char *tmpSymbol;
   tmpSymbol=mSymbol.c_str();
   const T_SF_WK95_CAA *atScattFact;
   atScattFact=FindSF_WK95_CAA( tmpSymbol, 0);
   VFN_DEBUG_MESSAGE("ScatteringPowerAtom::InitAtScattCoeffsWK95():1",2)
   for(int i=0;i<5;i++)
   {
      mScattAi(i)=atScattFact->a[i];
      mScattBi(i)=atScattFact->b[i];
   }
   mScattC=atScattFact->c;
   VFN_DEBUG_MESSAGE("ScatteringPowerAtom::InitAtScattCoeffsWK95():End",2)
}

void ScatteringPowerAtom::InitAtNeutronScattCoeffs()
{
   VFN_DEBUG_MESSAGE("ScatteringPowerAtom::InitAtNeutronScattCoeffs():"<<mName,3)
   mClock.Click();
   const char *tmpSymbol;
   tmpSymbol=mSymbol.c_str();
   const T_NeutronBondSL_NN92 *neutronScattFact;
   neutronScattFact=FindNeutronBondSL_NN92( tmpSymbol, 0);
   mNeutronScattLengthReal=neutronScattFact->BondCohScattLength;
   mNeutronScattLengthImag=neutronScattFact->BondCohScattLengthImag;
   //mNeutronAbsCrossSection=neutronScattFact->AbsCrossSect;
   VFN_DEBUG_MESSAGE("ScatteringPowerAtom::InitAtNeutronScattCoeffs():End",3)
}

void ScatteringPowerAtom::InitRefParList()
{
   VFN_DEBUG_MESSAGE("ScatteringPowerAtom::InitRefParList():"<<mName,5)
   {
      RefinablePar tmp(this->GetName()+"_Biso",&mBiso,0.1,5.,
                        gpRefParTypeScattPowTemperatureIso,REFPAR_DERIV_STEP_ABSOLUTE,
                        true,true,true,false);
      tmp.SetDerivStep(1e-3);
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
   cout <<mX<<" "<<mY<<" "<<mZ<<" "<<mOccupancy<<" "<<mDynPopCorr<<" "<<mpScattPow
        <<" "<<mpScattPow->GetScatteringPowerId()<<" "<<mpScattPow->GetName()<<endl;
}

//######################################################################
//
//      SCATTERING COMPONENT LIST
//
//######################################################################

ScatteringComponentList::ScatteringComponentList():
mNbComponent(0),mpScattComp(0),mMaxNbComponent(0)
{
}

ScatteringComponentList::ScatteringComponentList(const long nbComponent):
mNbComponent(nbComponent),mpScattComp(0),mMaxNbComponent(nbComponent)
{
   mpScattComp=new ScatteringComponent[mMaxNbComponent];
}

ScatteringComponentList::ScatteringComponentList(const ScatteringComponentList &old):
mNbComponent(old.mNbComponent),mMaxNbComponent(old.mMaxNbComponent)
{
   mpScattComp=new ScatteringComponent[mMaxNbComponent];
   for(long i=0;i<mNbComponent;i++) *(mpScattComp+i) = old(i);
}

ScatteringComponentList::~ScatteringComponentList()
{
   if(0 != mpScattComp) delete[] mpScattComp;
}
void ScatteringComponentList::Reset()
{
   mNbComponent=0;
}

const ScatteringComponent& ScatteringComponentList::operator()(const long i) const
{
   VFN_DEBUG_MESSAGE("ScatteringComponentList::operator()("<<i<<")",1)
   if(i>=mNbComponent)
   {
      this->Print();
      throw ObjCrystException("ScatteringComponentList::operator()(i)::\
 i>mNbComponent!!");
   }
   if(i<0) throw ObjCrystException("ScatteringComponentList::operator()&(i)::\
 i<0!!");
   return *(mpScattComp+i);
}

ScatteringComponent& ScatteringComponentList::operator()(const long i)
{
   VFN_DEBUG_MESSAGE("ScatteringComponentList::operator()&("<<i<<")",1)
   if(i>=mNbComponent)
   {
      this->Print();
      if(i>=mNbComponent) throw ObjCrystException("ScatteringComponentList::operator()&(i)::\
 i>mNbComponent!!");
   }
   if(i<0) throw ObjCrystException("ScatteringComponentList::operator()&(i)::\
 i<0!!");
   return *(mpScattComp+i);
}

long ScatteringComponentList::GetNbComponent() const {return mNbComponent;}

void ScatteringComponentList::operator=(const ScatteringComponentList &rhs)
{
   VFN_DEBUG_MESSAGE("ScatteringComponentList::operator=()",1)
   this->Reset();
   const long newNbComp=rhs.GetNbComponent();
   if(newNbComp > mMaxNbComponent)
   {//need to resize
      mMaxNbComponent=newNbComp;
      if(mpScattComp!=0) delete[] mpScattComp;
      //VFN_DEBUG_MESSAGE("ScatteringComponentList::operator=():1",2)
      mpScattComp=new ScatteringComponent[mMaxNbComponent];
      //VFN_DEBUG_MESSAGE("ScatteringComponentList::operator=():2",2)
   }
   mNbComponent=newNbComp;
   for(long i=0;i<mNbComponent;i++)
   {
      VFN_DEBUG_MESSAGE(".."<<i,0)
      *(mpScattComp+i) = rhs(i);
      //*(mpScattComp+i) = rhs(i-mNbComponent);
   }
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
   const long newNbComp=mNbComponent+rhs.GetNbComponent();
   if(newNbComp >= mMaxNbComponent)
   {//need to resize
      ScatteringComponent* pScattComp=mpScattComp;
      //:TODO: Check this works all right (error messages in window$ to this line, among others)
      mpScattComp=new ScatteringComponent[newNbComp];
      mMaxNbComponent=newNbComp;
      if(pScattComp!=0)
      {
         for(long i=0;i<mNbComponent;i++) *(mpScattComp+i) = *(pScattComp+i);
         delete[] pScattComp;
      }
   }
   for(long i=mNbComponent;i<newNbComp;i++) 
      *(mpScattComp+i) = rhs(i-mNbComponent);
   mNbComponent=newNbComp;
}
void ScatteringComponentList::operator+=(const ScatteringComponent &rhs)
{
   VFN_DEBUG_MESSAGE("ScatteringComponentList::operator+=()",1)
   ++(*this);
   *(mpScattComp + mNbComponent-1) = rhs;
}

void ScatteringComponentList::operator++()
{
   VFN_DEBUG_MESSAGE("ScatteringComponentList::operator++()",1)
   const long newNbComp=mNbComponent+1;
   if(newNbComp >= mMaxNbComponent)
   {//need to resize
      ScatteringComponent* pScattComp=mpScattComp;
      mpScattComp=new ScatteringComponent[newNbComp];
      mMaxNbComponent=newNbComp;
      if(pScattComp!=0)
      {
         for(long i=0;i<mNbComponent;i++) *(mpScattComp+i) = *(pScattComp+i);
         delete[] pScattComp;
      }
   }
   mNbComponent++;
}

void ScatteringComponentList::Print()const
{
   VFN_DEBUG_MESSAGE("ScatteringComponentList::Print()",5)
   cout<<"Number of Scattering components:"<<this->GetNbComponent()<<endl;
   for(long i=0;i<this->GetNbComponent();i++)
   {
      cout << i<<":";
      (mpScattComp+i)->Print();
   }
}

void ScatteringComponentList::ChangeMaxNbComponent(const long num)
{
   VFN_DEBUG_MESSAGE("ScatteringComponentList::ChangeMaxNbComponents(num)",2)
   if(num >= mMaxNbComponent)
   {//need to resize
      mMaxNbComponent=num;
      ScatteringComponent* pScattComp=mpScattComp;
      mpScattComp=new ScatteringComponent[mMaxNbComponent];
      if(pScattComp!=0)
      {
         for(long i=0;i<mNbComponent;i++) *(mpScattComp+i) = *(pScattComp+i);
         delete[] pScattComp;
      }
   }
}

}//namespace
