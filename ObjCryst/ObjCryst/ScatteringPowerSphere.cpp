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
/*   ScatteringPowerSphere.cpp
*  source file for the spherical scattering power
*
*/

#include <cmath>

#include "ObjCryst/ObjCryst/ScatteringPowerSphere.h"
#include "ObjCryst/Quirks/VFNStreamFormat.h"

#ifdef __WX__CRYST__
#include "ObjCryst/wxCryst/wxScatteringPowerSphere.h"
#endif

namespace ObjCryst
{
////////////////////////////////////////////////////////////////////////
//
//    ScatteringPowerSphere
//
////////////////////////////////////////////////////////////////////////

ScatteringPowerSphere::ScatteringPowerSphere()
{
   cout<<"Creating scattering power for fullerene"<<endl;
   this->Init("C60",3.56,1.0);
}

ScatteringPowerSphere::ScatteringPowerSphere(const string &name,
                                             const REAL radius,
                                             const REAL bIso)
{
   this->Init(name,radius,bIso);
}

ScatteringPowerSphere::ScatteringPowerSphere(const ScatteringPowerSphere& old) 
{}

ScatteringPowerSphere::~ScatteringPowerSphere()
{
}

void ScatteringPowerSphere::Init(const string &name,
                                             const REAL radius,
                                             const REAL bIso)
{
   this->SetName(name);
   mRadius=radius;
   mColourRGB[0]=0.2;
   mColourRGB[1]=0.2;
   mColourRGB[2]=0.2;
   if(this->GetNbPar()==0) this->InitRefParList();
}

const string& ScatteringPowerSphere::GetClassName() const
{
   const static string className="ScatteringPowerSphere";
   return className;
}

CrystVector_REAL ScatteringPowerSphere::GetScatteringFactor(const ScatteringData &data,
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
   sf *= 2*M_PI*mRadius;
   REAL *pSF=sf.data();
   for(long i=0;i<data.GetNbReflBelowMaxSinThetaOvLambda();i++)
   {
      *pSF = sin(*pSF) / *pSF;
      pSF++;
   }
   #endif
   return sf;
}

REAL ScatteringPowerSphere::GetForwardScatteringFactor(const RadiationType type) const
{
   return 1;
}

CrystVector_REAL ScatteringPowerSphere::GetTemperatureFactor(const ScatteringData &data,
                                       const int spgSymPosIndex) const
{
   VFN_DEBUG_MESSAGE("ScatteringPowerSphere::GetTemperatureFactor(&data):"<<mName,3)
   CrystVector_REAL sf(data.GetNbRefl());
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
   return sf;
}

CrystMatrix_REAL ScatteringPowerSphere::GetResonantScattFactReal(const ScatteringData &data,
                                       const int spgSymPosIndex) const
{
   CrystMatrix_REAL fprime(1,1);//:TODO: more than one lambda
   fprime=0;
   return fprime;
}

CrystMatrix_REAL ScatteringPowerSphere::GetResonantScattFactImag(const ScatteringData &data,
                                                     const int spgSymPosIndex) const
{
   CrystMatrix_REAL fsecond(1,1);//:TODO: more than one lambda
   fsecond=0;
   return fsecond;
}

REAL ScatteringPowerSphere::GetRadius()const{return mRadius;}

void ScatteringPowerSphere::Print()const
{
   
}

void ScatteringPowerSphere::XMLOutput(ostream &os,int indent)const
{
   VFN_DEBUG_ENTRY("ScatteringPowerSphere::XMLOutput():"<<this->GetName(),5)
   for(int i=0;i<indent;i++) os << "  " ;
   XMLCrystTag tag("ScatteringPowerSphere");
   tag.AddAttribute("Name",mName);
   os <<tag<<endl;
      
   for(int i=0;i<=indent;i++) os << "  " ;
      this->GetPar(&mRadius).XMLOutput(os,"Radius",0);
   os<<endl;
   
   for(int i=0;i<=indent;i++) os << "  " ;
      this->GetPar(&mBiso).XMLOutput(os,"Biso",0);
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
   VFN_DEBUG_EXIT("ScatteringPowerSphere::XMLOutput():"<<this->GetName(),5)
}

void ScatteringPowerSphere::XMLInput(istream &is,const XMLCrystTag &tagg)
{
   VFN_DEBUG_ENTRY("ScatteringPowerSphere::XMLInput():"<<this->GetName(),5)
   string symbol;
   for(unsigned int i=0;i<tagg.GetNbAttribute();i++)
   {
      if("Name"==tagg.GetAttributeName(i)) this->SetName(tagg.GetAttributeValue(i));
      if("Symbol"==tagg.GetAttributeName(i)) symbol=tagg.GetAttributeValue(i);
   }
   this->Init(mName,3.56,1.0);
   while(true)
   {
      XMLCrystTag tag(is);
      if(("ScatteringPowerSphere"==tag.GetName())&&tag.IsEndTag())
      {
         VFN_DEBUG_EXIT("ScatteringPowerSphere::Exit():"<<this->GetName(),5)
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
               if("Radius"==tag.GetAttributeValue(i)) this->GetPar(&mRadius).XMLInput(is,tag);
               if("Biso"==tag.GetAttributeValue(i)) this->GetPar(&mBiso).XMLInput(is,tag);
               break;
            }
         }
         continue;
      }
   }
}
void ScatteringPowerSphere::InitRefParList()
{
   {
      RefinablePar tmp("Radius",&mRadius,3.,5.,
                        gpRefParTypeScattPow,REFPAR_DERIV_STEP_RELATIVE,
                        true,true,true,false);
      tmp.SetDerivStep(1e-3);
      tmp.SetGlobalOptimStep(.1);
      tmp.AssignClock(mClock);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("Biso",&mBiso,0.1,5.,
                        gpRefParTypeScattPowTemperatureIso,REFPAR_DERIV_STEP_ABSOLUTE,
                        true,true,true,false);
      tmp.SetDerivStep(1e-3);
      tmp.SetGlobalOptimStep(.5);
      tmp.AssignClock(mClock);
      this->AddPar(tmp);
   }
}

#ifdef __WX__CRYST__

WXCrystObjBasic* ScatteringPowerSphere::WXCreate(wxWindow* parent)
{
   //:TODO: Check mpWXCrystObj==0
   mpWXCrystObj=new WXScatteringPowerSphere(parent,this);
   return mpWXCrystObj;
}
#endif


}//namespace
