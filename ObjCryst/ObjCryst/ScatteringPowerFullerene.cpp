/*  ObjCryst++ Object-Oriented Crystallographic Library
    (c) 2000-2002 Vincent Favre-Nicolin vincefn@users.sourceforge.net
        2000-2001 University of Geneva (Switzerland)

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; version 2 of the License.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/*   ScatteringPowerFullerene.cpp
*  source file for the Fullerene-type scattering power
*
*/

#include <cmath>

#include "ObjCryst/ScatteringPowerFullerene.h"
#include "Quirks/VFNStreamFormat.h"

#ifdef __WX__CRYST__
#include "wxCryst/wxScatteringPowerFullerene.h"
#endif

namespace ObjCryst
{
////////////////////////////////////////////////////////////////////////
//
//    ScatteringPowerFullerene
//
////////////////////////////////////////////////////////////////////////

ScatteringPowerFullerene::ScatteringPowerFullerene():
mpScatteringPower(0)
{
   cout<<"Creating scattering power for fullerene"<<endl;
   this->Init("C60",60,3.56,3.56,3.56,1.0,"C");
}

ScatteringPowerFullerene::ScatteringPowerFullerene(const string &name,
                                                   const unsigned int nbAtom,
                                                   const REAL axisLengthX,
                                                   const REAL axisLengthY,
                                                   const REAL axisLengthZ,
                                                   const REAL bIso,
                                                   const string symbol)
{
   this->Init(name,nbAtom,axisLengthX,axisLengthY,axisLengthZ,bIso,symbol);
}

ScatteringPowerFullerene::~ScatteringPowerFullerene()
{
   if(mpScatteringPower!=0) delete mpScatteringPower;
}

void ScatteringPowerFullerene::Init(const string &name,
                                    const unsigned int nbAtom,
                                    const REAL axisLengthX,
                                    const REAL axisLengthY,
                                    const REAL axisLengthZ,
                                    const REAL bIso,
                                    const string symbol)
{
   this->SetName(name);
   mNbAtom=nbAtom;
   mAxisLengthX=axisLengthX;
   mAxisLengthY=axisLengthY;
   mAxisLengthZ=axisLengthZ;
   mPhi=0;
   mChi=0;
   mPsi=0;
   if(mpScatteringPower!=0) delete mpScatteringPower;
   mpScatteringPower=new ScatteringPowerAtom(name+"_atom",symbol,bIso);
   mDynPopCorrIndex=nbAtom*mpScatteringPower->GetDynPopCorrIndex();
   mColourRGB[0]=0.2;
   mColourRGB[1]=0.2;
   mColourRGB[2]=0.2;
   if(this->GetNbPar()==0) this->InitRefParList();
}

const string& ScatteringPowerFullerene::GetClassName() const
{
   static string className="ScatteringPowerFullerene";
   return className;
}

CrystVector_REAL ScatteringPowerFullerene::GetScatteringFactor(const ScatteringData &data,
                                       const int spgSymPosIndex) const
{  
   
   CrystVector_REAL sf;
   #if 0
   // :TODO: Support anisotropic form factor. For the moment use only mAxisLengthX
   REAL a,b,c;
   CrystMatrix_REAL phiMatrix(3,3),chiMatrix(3,3),psiMatrix(3,3);
   phiMatrix= cos(mPhi)   , -sin(mPhi)   , 0,
              sin(mPhi)   , cos(mPhi)    , 0,
              0           ,0             ,1;

   chiMatrix= cos(mChi)   ,0             ,-sin(mChi),
              0           ,1             ,0,
              sin(mChi)   ,0             ,cos(mChi);

   psiMatrix= 1           , 0            , 0,
              0           ,cos(mPsi)     ,-sin(mPsi),
              0           ,sin(mPsi)     ,cos(mPsi);

   CrystMatrix_REAL phiChiPsiMatrix=product(chiMatrix,product(phiMatrix,psiMatrix));

   a=phiChiPsiMatrix(0,0)*mAxisLengthX+phiChiPsiMatrix(0,1)*mAxisLengthY+phiChiPsiMatrix(0,2)*mAxisLengthZ;
   b=phiChiPsiMatrix(1,0)*mAxisLengthX+phiChiPsiMatrix(1,1)*mAxisLengthY+phiChiPsiMatrix(1,2)*mAxisLengthZ;
   c=phiChiPsiMatrix(2,0)*mAxisLengthX+phiChiPsiMatrix(2,1)*mAxisLengthY+phiChiPsiMatrix(2,2)*mAxisLengthZ;
   const REAL *pX=data.GetReflX().data();
   const REAL *pY=data.GetReflY().data();
   const REAL *pZ=data.GetReflZ().data();
   sf.resize(data.GetNbRefl());
   REAL *pSF=sf.data();
   for(long i=0;i<data.GetNbReflBelowMaxSinThetaOvLambda();i++)
      *pSF++ = 2*M_PI * ( *pX++ * a +  *pY++ * b + *pZ++ * c );
   pSF=sf.data();
   for(long i=0;i<data.GetNbReflBelowMaxSinThetaOvLambda();i++)
   {
      *pSF = sin(*pSF) / *pSF;
      *pSF++;
   }
   #else
   sf=data.GetSinThetaOverLambda();
   sf *= 2*M_PI*mAxisLengthX;
   REAL *pSF=sf.data();
   for(long i=0;i<data.GetNbReflBelowMaxSinThetaOvLambda();i++)
   {
      *pSF = sin(*pSF) / *pSF;
      *pSF++;
   }
   #endif
   sf *= mpScatteringPower->GetScatteringFactor(data,spgSymPosIndex);
   sf *= mNbAtom;
   cout << FormatVertVectorHKLFloats<REAL>(data.GetH(),data.GetK(),data.GetL(),
                                     data.GetSinThetaOverLambda(),
                                     mpScatteringPower->GetScatteringFactor(data,spgSymPosIndex),
                                     sf)<<endl;
   return sf;
}

REAL ScatteringPowerFullerene::GetForwardScatteringFactor(const RadiationType type) const
{
   return mpScatteringPower->GetForwardScatteringFactor(type)*mNbAtom;
}

CrystVector_REAL ScatteringPowerFullerene::GetTemperatureFactor(const ScatteringData &data,
                                       const int spgSymPosIndex) const
{
   return mpScatteringPower->GetTemperatureFactor(data, spgSymPosIndex);
}

CrystMatrix_REAL ScatteringPowerFullerene::GetResonantScattFactReal(const ScatteringData &data,
                                       const int spgSymPosIndex) const
{
   CrystMatrix_REAL rsfr;
   rsfr=mpScatteringPower->GetResonantScattFactReal(data,spgSymPosIndex);
   rsfr *= mNbAtom;
   return rsfr;
}

CrystMatrix_REAL ScatteringPowerFullerene::GetResonantScattFactImag(const ScatteringData &data,
                                                     const int spgSymPosIndex) const
{
   CrystMatrix_REAL rsfi;
   rsfi=mpScatteringPower->GetResonantScattFactImag(data,spgSymPosIndex);
   rsfi *= mNbAtom;
   return rsfi;
}

REAL ScatteringPowerFullerene::GetAxisLengthX() const{return mAxisLengthX;}
REAL ScatteringPowerFullerene::GetAxisLengthY() const{return mAxisLengthY;}
REAL ScatteringPowerFullerene::GetAxisLengthZ() const{return mAxisLengthZ;}

REAL ScatteringPowerFullerene::GetRadius()const{return mAxisLengthX;}

void ScatteringPowerFullerene::Print()const
{
   
}

void ScatteringPowerFullerene::XMLOutput(ostream &os,int indent)const
{
   VFN_DEBUG_ENTRY("ScatteringPowerFullerene::XMLOutput():"<<this->GetName(),5)
   for(int i=0;i<indent;i++) os << "  " ;
   XMLCrystTag tag("ScatteringPowerFullerene");
   tag.AddAttribute("Name",mName);
   tag.AddAttribute("Symbol",mpScatteringPower->GetSymbol());
   os <<tag<<endl;
   
   for(int i=0;i<=indent;i++) os << "  " ;
      this->GetPar(&mNbAtom).XMLOutput(os,"NbAtom",0);
   os<<endl;
   
   for(int i=0;i<=indent;i++) os << "  " ;
      this->GetPar(&mAxisLengthX).XMLOutput(os,"AxisLengthX",0);
   os<<endl;
   
   for(int i=0;i<=indent;i++) os << "  " ;
      mpScatteringPower
         ->GetPar(&(mpScatteringPower->GetBiso())).XMLOutput(os,"Biso",0);
   os<<endl;
   
   for(int i=0;i<=indent;i++) os << "  " ;
   XMLCrystTag tag2("RGBColour");
   os << tag2
      << mColourRGB[0]<<" "
      << mColourRGB[1]<<" "
      << mColourRGB[2];
   tag2.SetIsEndTag(true);
   os << tag2<<endl;
   
   tag.SetIsEndTag(true);
   for(int i=0;i<indent;i++) os << "  " ;
   os <<tag<<endl;
   VFN_DEBUG_EXIT("ScatteringPowerFullerene::XMLOutput():"<<this->GetName(),5)
}

void ScatteringPowerFullerene::XMLInput(istream &is,const XMLCrystTag &tagg)
{
   VFN_DEBUG_ENTRY("ScatteringPowerFullerene::XMLInput():"<<this->GetName(),5)
   string symbol;
   for(unsigned int i=0;i<tagg.GetNbAttribute();i++)
   {
      if("Name"==tagg.GetAttributeName(i)) this->SetName(tagg.GetAttributeValue(i));
      if("Symbol"==tagg.GetAttributeName(i)) symbol=tagg.GetAttributeValue(i);
   }
   this->Init(mName,60,3.56,3.56,3.56,1.0,symbol);
   while(true)
   {
      XMLCrystTag tag(is);
      if(("ScatteringPowerFullerene"==tag.GetName())&&tag.IsEndTag())
      {
         VFN_DEBUG_EXIT("ScatteringPowerFullerene::Exit():"<<this->GetName(),5)
         return;
      }
      if("RGBColour"==tag.GetName())
      {
         float r,g,b;
         is>>r>>g>>b;
         this->SetColour(r,g,b);
         XMLCrystTag junk(is);
      }
      if("Par"==tag.GetName())
      {
         for(unsigned int i=0;i<tag.GetNbAttribute();i++)
         {
            if("Name"==tag.GetAttributeName(i))
            {  cout <<tag<<endl;
               if("NbAtom"==tag.GetAttributeValue(i)) this->GetPar(&mNbAtom).XMLInput(is,tag);
               if("AxisLengthX"==tag.GetAttributeValue(i)) this->GetPar(&mAxisLengthX).XMLInput(is,tag);
               if("Biso"==tag.GetAttributeValue(i))
                  mpScatteringPower
                     ->GetPar(&(mpScatteringPower->GetBiso())).XMLInput(is,tag);
               break;
            }
         }
         continue;
      }
   }
}
void ScatteringPowerFullerene::InitRefParList()
{
   #if 1
   {
      RefinablePar tmp("Radius",&mAxisLengthX,3.,5.,
                        gpRefParTypeScattPow,REFPAR_DERIV_STEP_RELATIVE,
                        true,true,true,false);
      tmp.SetDerivStep(1e-3);
      tmp.SetGlobalOptimStep(.1);
      tmp.AssignClock(mClock);
      this->AddPar(tmp);
   }
   #else
   // :TODO: Anisotropic support
   {
      RefinablePar tmp("RadiusX",&mAxisLengthX,2.,5.,
                        gpRefParTypeScattPow,REFPAR_DERIV_STEP_RELATIVE,
                        true,true,true,false);
      tmp.SetDerivStep(1e-3);
      tmp.SetGlobalOptimStep(.1);
      tmp.AssignClock(mClock);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("RadiusY",&mAxisLengthY,2.,5.,
                        gpRefParTypeScattPow,REFPAR_DERIV_STEP_RELATIVE,
                        true,true,true,false);
      tmp.SetDerivStep(1e-3);
      tmp.SetGlobalOptimStep(.1);
      tmp.AssignClock(mClock);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("RadiusZ",&mAxisLengthZ,2.,5.,
                        gpRefParTypeScattPow,REFPAR_DERIV_STEP_RELATIVE,
                        true,true,true,false);
      tmp.SetDerivStep(1e-3);
      tmp.SetGlobalOptimStep(.1);
      tmp.AssignClock(mClock);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("Phi",&mPhi,0.,2*M_PI,
                        gpRefParTypeScattPow,REFPAR_DERIV_STEP_ABSOLUTE,
                        false,true,true,true,RAD2DEG,2*M_PI);
      tmp.AssignClock(mClock);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("Chi",&mChi,0.,2*M_PI,
                        gpRefParTypeScattPow,REFPAR_DERIV_STEP_ABSOLUTE,
                        false,true,true,true,RAD2DEG,2*M_PI);
      tmp.AssignClock(mClock);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("Psi",&mPsi,0.,2*M_PI,
                        gpRefParTypeScattPow,REFPAR_DERIV_STEP_ABSOLUTE,
                        false,true,true,true,RAD2DEG,2*M_PI);
      tmp.AssignClock(mClock);
      this->AddPar(tmp);
   }
   #endif
   {
      RefinablePar tmp("NbAtom",&mNbAtom,1.,100.,
                        gpRefParTypeScattPow,REFPAR_DERIV_STEP_RELATIVE,
                        true,true,true,false);
      tmp.SetDerivStep(1e-3);
      tmp.SetGlobalOptimStep(.1);
      tmp.AssignClock(mClock);
      this->AddPar(tmp);
   }
}

#ifdef __WX__CRYST__

WXCrystObjBasic* ScatteringPowerFullerene::WXCreate(wxWindow* parent)
{
   //:TODO: Check mpWXCrystObj==0
   mpWXCrystObj=new WXScatteringPowerFullerene(parent,this);
   return mpWXCrystObj;
}
#endif


}//namespace
