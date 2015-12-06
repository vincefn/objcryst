/*  ObjCryst++ Object-Oriented Crystallographic Library
    (c) 2000-2005 Vincent Favre-Nicolin vincefn@users.sourceforge.net
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
*  source file for LibCryst++ ReflectionProfile and derived classes
*
*/
#include <limits>
#include "ObjCryst/ObjCryst/ReflectionProfile.h"
#include "ObjCryst/Quirks/VFNStreamFormat.h"
#ifdef __WX__CRYST__
   #include "ObjCryst/wxCryst/wxPowderPattern.h"
#endif

#ifdef HAVE_SSE_MATHFUN
#include "ObjCryst/Quirks/sse_mathfun.h"
#endif

namespace ObjCryst
{
#if defined(_MSC_VER) || defined(__BORLANDC__)
#undef min // Predefined macros.... (wx?)
#undef max

double erfc(const double x)// in C99, but not in VC++....
{
   if(x<0.0) return 2.0-erfc(-x);
   if(x<3.8)
   { // Series, Abramowitz & Stegun 7.1.6
      double y=x,y0=x;
      for(int i=1;i<=50;i++)
      {
         y0*=2*x*x/(2*i+1.0);
         y+=y0;
      }
      static const double spi=2/sqrt(M_PI);
      return 1-spi*exp(-x*x)*y;
   }
   double y=1.0,y0=1.0;
   for(int i=1;i<=10;i++)
   {// Asymptotic, Abramowitz & Stegun 7.1.23
      y0*=-(2*i-1)/(2*x*x);
      y+=y0;
   }
   static const double invsqrtpi=1.0/sqrt(M_PI);
   return invsqrtpi*exp(-x*x)/x*y;
}

#endif
extern const RefParType *gpRefParTypeScattDataProfile;
extern const RefParType *gpRefParTypeScattDataProfileType;
extern const RefParType *gpRefParTypeScattDataProfileWidth;
extern const RefParType *gpRefParTypeScattDataProfileAsym;

ObjRegistry<ReflectionProfile>
   gReflectionProfileRegistry("List of all ReflectionProfile types");;
////////////////////////////////////////////////////////////////////////
//
//    ReflectionProfile
//
////////////////////////////////////////////////////////////////////////
ReflectionProfile::ReflectionProfile():
RefinableObj()
{}
ReflectionProfile::ReflectionProfile(const ReflectionProfile &old)
{}
ReflectionProfile::~ReflectionProfile()
{}
bool ReflectionProfile::IsAnisotropic()const
{return false;}
////////////////////////////////////////////////////////////////////////
//
//    ReflectionProfilePseudoVoigt
//
////////////////////////////////////////////////////////////////////////
ReflectionProfilePseudoVoigt::ReflectionProfilePseudoVoigt():
ReflectionProfile(),
mCagliotiU(0),mCagliotiV(0),mCagliotiW(.01*DEG2RAD*DEG2RAD),
mPseudoVoigtEta0(0.5),mPseudoVoigtEta1(0.0),
mAsymBerarBaldinozziA0(0.0),mAsymBerarBaldinozziA1(0.0),
mAsymBerarBaldinozziB0(0.0),mAsymBerarBaldinozziB1(0.0),
mAsym0(1.0),mAsym1(0.0),mAsym2(0.0)
{
   this->InitParameters();
}

ReflectionProfilePseudoVoigt::ReflectionProfilePseudoVoigt
   (const ReflectionProfilePseudoVoigt &old):
mCagliotiU(old.mCagliotiU),mCagliotiV(old.mCagliotiV),mCagliotiW(old.mCagliotiW),
mPseudoVoigtEta0(old.mPseudoVoigtEta0),mPseudoVoigtEta1(old.mPseudoVoigtEta1),
mAsymBerarBaldinozziA0(old.mAsymBerarBaldinozziA0),
mAsymBerarBaldinozziA1(old.mAsymBerarBaldinozziA1),
mAsymBerarBaldinozziB0(old.mAsymBerarBaldinozziB0),
mAsymBerarBaldinozziB1(old.mAsymBerarBaldinozziB1),
mAsym0(old.mAsym0),mAsym1(old.mAsym1),mAsym2(old.mAsym2)
{
   this->InitParameters();
}

ReflectionProfilePseudoVoigt::~ReflectionProfilePseudoVoigt()
{
   #ifdef __WX__CRYST__
   if(mpWXCrystObj!=0)
   {
      delete mpWXCrystObj;
      mpWXCrystObj=0;
   }
   #endif
}

ReflectionProfilePseudoVoigt* ReflectionProfilePseudoVoigt::CreateCopy()const
{
   return new ReflectionProfilePseudoVoigt(*this);
}

const string& ReflectionProfilePseudoVoigt::GetClassName()const
{
   static string className="ReflectionProfilePseudoVoigt";
   return className;
}

CrystVector_REAL ReflectionProfilePseudoVoigt::GetProfile(const CrystVector_REAL &x,
                            const REAL center,const REAL h, const REAL k, const REAL l)const
{
   VFN_DEBUG_ENTRY("ReflectionProfilePseudoVoigt::GetProfile(),c="<<center,2)
   REAL fwhm= mCagliotiW
             +mCagliotiV*tan(center/2.0)
             +mCagliotiU*pow(tan(center/2.0),2);
   if(fwhm<=0)
   {
      VFN_DEBUG_MESSAGE("ReflectionProfilePseudoVoigt::GetProfile(): fwhm**2<0 ! "
          <<h<<","<<k<<","<<l<<":"<<center<<","<<mCagliotiU<<","<<mCagliotiV<<","<<","<<mCagliotiW<<":"<<fwhm,10);
      fwhm=1e-6;
   }
   else fwhm=sqrt(fwhm);
   CrystVector_REAL profile,tmpV;
   const REAL asym=mAsym0+mAsym1/sin(center)+mAsym2/pow((REAL)sin(center),(REAL)2.0);
   profile=PowderProfileGauss(x,fwhm,center,asym);
   
   // Eta for gaussian/lorentzian mix. Make sure 0<=eta<=1, else profiles could be <0 !
   REAL eta=mPseudoVoigtEta0+center*mPseudoVoigtEta1;
   if(eta>1) eta=1;
   if(eta<0) eta=0;
   
   profile *= 1-eta;
   tmpV=PowderProfileLorentz(x,fwhm,center,asym);
   tmpV *= eta;
   profile += tmpV;
   //profile *= AsymmetryBerarBaldinozzi(x,fwhm,center,
   //                                    mAsymBerarBaldinozziA0,mAsymBerarBaldinozziA1,
   //                                    mAsymBerarBaldinozziB0,mAsymBerarBaldinozziB1);
   VFN_DEBUG_EXIT("ReflectionProfilePseudoVoigt::GetProfile()",2)
   return profile;
}

void ReflectionProfilePseudoVoigt::SetProfilePar(const REAL fwhmCagliotiW,
                   const REAL fwhmCagliotiU,
                   const REAL fwhmCagliotiV,
                   const REAL eta0,
                   const REAL eta1)
{
   mCagliotiU=fwhmCagliotiU;
   mCagliotiV=fwhmCagliotiV;
   mCagliotiW=fwhmCagliotiW;
   mPseudoVoigtEta0=eta0;
   mPseudoVoigtEta1=eta1;
   mClockMaster.Click();
}

bool ReflectionProfilePseudoVoigt::IsAnisotropic()const{return false;}
REAL ReflectionProfilePseudoVoigt::GetFullProfileWidth(const REAL relativeIntensity,
                            const REAL center,const REAL h, const REAL k, const REAL l)
{
   VFN_DEBUG_ENTRY("ReflectionProfilePseudoVoigt::GetFullProfileWidth()",2)
   const int nb=100;
   const int halfnb=nb/2;
   CrystVector_REAL x(nb);
   REAL n=5.0;
   REAL fwhm= mCagliotiW
             +mCagliotiV*tan(center/2.0)
             +mCagliotiU*pow(tan(center/2.0),2);
   if(fwhm<=0) fwhm=1e-6;
   else fwhm=sqrt(fwhm);
   CrystVector_REAL prof;
   while(true)
   {
      //Create an X array with 100 elements reaching +/- n*FWHM/2
      REAL *p=x.data();
      const REAL tmp=fwhm*n/nb;
      for(int i=0;i<nb;i++) *p++ = tmp*(i-halfnb);
      x+=center;
      
      prof=this->GetProfile(x,center,0,0,0);
      const REAL max=prof.max();
      const REAL test=max*relativeIntensity;
      int n1=0,n2=0;
      if((prof(0)<test)&&(prof(nb-1)<test))
      {
         p=prof.data();
         while(*p<test){ p++; n1++;n2++;}
         n1--;
         while(*p>test){ p++; n2++;}
         VFN_DEBUG_EXIT("ReflectionProfilePseudoVoigt::GetFullProfileWidth():"<<x(n2)-x(n1),10)
         return x(n2)-x(n1);
      }
      VFN_DEBUG_MESSAGE("ReflectionProfilePseudoVoigt::GetFullProfileWidth():"<<relativeIntensity<<","
                        <<fwhm<<","<<center<<","<<h<<","<<k<<","<<l<<","<<max<<","<<test,2)
      VFN_DEBUG_MESSAGE(FormatVertVector<REAL>(x,prof),2)
      n*=2.0;
      //if(n>200) exit(0);
   }
}

void ReflectionProfilePseudoVoigt::InitParameters()
{
   {
      RefinablePar tmp("U",&mCagliotiU,-1/RAD2DEG/RAD2DEG,1./RAD2DEG/RAD2DEG,
                        gpRefParTypeScattDataProfileWidth,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,RAD2DEG*RAD2DEG);
      tmp.AssignClock(mClockMaster);
      tmp.SetDerivStep(1e-9);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("V",&mCagliotiV,-1/RAD2DEG/RAD2DEG,1./RAD2DEG/RAD2DEG,
                        gpRefParTypeScattDataProfileWidth,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,RAD2DEG*RAD2DEG);
      tmp.AssignClock(mClockMaster);
      tmp.SetDerivStep(1e-9);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("W",&mCagliotiW,0,1./RAD2DEG/RAD2DEG,
                        gpRefParTypeScattDataProfileWidth,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,RAD2DEG*RAD2DEG);
      tmp.AssignClock(mClockMaster);
      tmp.SetDerivStep(1e-9);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("Eta0",&mPseudoVoigtEta0,0,1.,gpRefParTypeScattDataProfileType,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false);
      tmp.AssignClock(mClockMaster);
      tmp.SetDerivStep(1e-4);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("Eta1",&mPseudoVoigtEta1,-1,1.,gpRefParTypeScattDataProfileType,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false);
      tmp.AssignClock(mClockMaster);
      tmp.SetDerivStep(1e-4);
      this->AddPar(tmp);
   }
   #if 0
   {
      RefinablePar tmp("AsymA0",&mAsymBerarBaldinozziA0,-0.05,0.05,gpRefParTypeScattDataProfileAsym,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false);
      tmp.AssignClock(mClockMaster);
      tmp.SetDerivStep(1e-4);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("AsymA1",&mAsymBerarBaldinozziA1,-0.05,0.05,gpRefParTypeScattDataProfileAsym,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false);
      tmp.AssignClock(mClockMaster);
      tmp.SetDerivStep(1e-4);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("AsymB0",&mAsymBerarBaldinozziB0,-0.01,0.01,gpRefParTypeScattDataProfileAsym,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false);
      tmp.AssignClock(mClockMaster);
      tmp.SetDerivStep(1e-4);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("AsymB1",&mAsymBerarBaldinozziB1,-0.01,0.01,gpRefParTypeScattDataProfileAsym,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false);
      tmp.AssignClock(mClockMaster);
      tmp.SetDerivStep(1e-4);
      this->AddPar(tmp);
   }
   #endif
   {
      RefinablePar tmp("Asym0",&mAsym0,0.01,10.0,gpRefParTypeScattDataProfileAsym,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false);
      tmp.AssignClock(mClockMaster);
      tmp.SetDerivStep(1e-4);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("Asym1",&mAsym1,-1.0,1.0,gpRefParTypeScattDataProfileAsym,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false);
      tmp.AssignClock(mClockMaster);
      tmp.SetDerivStep(1e-4);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("Asym2",&mAsym2,-1.0,1.0,gpRefParTypeScattDataProfileAsym,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false);
      tmp.AssignClock(mClockMaster);
      tmp.SetDerivStep(1e-4);
      this->AddPar(tmp);
   }
}

void ReflectionProfilePseudoVoigt::XMLOutput(ostream &os,int indent)const
{
   VFN_DEBUG_ENTRY("ReflectionProfilePseudoVoigt::XMLOutput():"<<this->GetName(),5)
   for(int i=0;i<indent;i++) os << "  " ;
   XMLCrystTag tag("ReflectionProfilePseudoVoigt");
   os <<tag<<endl;
   indent++;

   this->GetPar(&mCagliotiU).XMLOutput(os,"U",indent);
   os <<endl;

   this->GetPar(&mCagliotiV).XMLOutput(os,"V",indent);
   os <<endl;

   this->GetPar(&mCagliotiW).XMLOutput(os,"W",indent);
   os <<endl;

   this->GetPar(&mPseudoVoigtEta0).XMLOutput(os,"Eta0",indent);
   os <<endl;

   this->GetPar(&mPseudoVoigtEta1).XMLOutput(os,"Eta1",indent);
   os <<endl;

   this->GetPar(&mAsym0).XMLOutput(os,"Asym0",indent);
   os <<endl;

   this->GetPar(&mAsym1).XMLOutput(os,"Asym1",indent);
   os <<endl;

   this->GetPar(&mAsym2).XMLOutput(os,"Asym2",indent);
   os <<endl;
   #if 0
   this->GetPar(&mAsymBerarBaldinozziA0).XMLOutput(os,"AsymA0",indent);
   os <<endl;

   this->GetPar(&mAsymBerarBaldinozziA1).XMLOutput(os,"AsymA1",indent);
   os <<endl;

   this->GetPar(&mAsymBerarBaldinozziB0).XMLOutput(os,"AsymB0",indent);
   os <<endl;

   this->GetPar(&mAsymBerarBaldinozziB1).XMLOutput(os,"AsymB1",indent);
   os <<endl;
   #endif
   indent--;
   tag.SetIsEndTag(true);
   for(int i=0;i<indent;i++) os << "  " ;
   os <<tag<<endl;
   VFN_DEBUG_EXIT("ReflectionProfilePseudoVoigt::XMLOutput():"<<this->GetName(),5)
}
void ReflectionProfilePseudoVoigt::XMLInput(istream &is,const XMLCrystTag &tagg)
{
   VFN_DEBUG_ENTRY("ReflectionProfilePseudoVoigt::XMLInput():"<<this->GetName(),5)
   for(unsigned int i=0;i<tagg.GetNbAttribute();i++)
   {
      if("Name"==tagg.GetAttributeName(i)) this->SetName(tagg.GetAttributeValue(i));
   }
   while(true)
   {
      XMLCrystTag tag(is);
      if(("ReflectionProfilePseudoVoigt"==tag.GetName())&&tag.IsEndTag())
      {
         this->UpdateDisplay();
         VFN_DEBUG_EXIT("ReflectionProfilePseudoVoigt::Exit():"<<this->GetName(),5)
         return;
      }
      if("Par"==tag.GetName())
      {
         for(unsigned int i=0;i<tag.GetNbAttribute();i++)
         {
            if("Name"==tag.GetAttributeName(i))
            {
               if("U"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&mCagliotiU).XMLInput(is,tag);
                  break;
               }
               if("V"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&mCagliotiV).XMLInput(is,tag);
                  break;
               }
               if("W"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&mCagliotiW).XMLInput(is,tag);
                  break;
               }
               if("Eta0"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&mPseudoVoigtEta0).XMLInput(is,tag);
                  break;
               }
               if("Eta1"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&mPseudoVoigtEta1).XMLInput(is,tag);
                  break;
               }
               if("Asym0"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&mAsym0).XMLInput(is,tag);
                  break;
               }
               if("Asym1"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&mAsym1).XMLInput(is,tag);
                  break;
               }
               if("Asym2"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&mAsym2).XMLInput(is,tag);
                  break;
               }
               #if 0
               if("AsymA0"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&mAsymBerarBaldinozziA0).XMLInput(is,tag);
                  break;
               }
               if("AsymA1"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&mAsymBerarBaldinozziA1).XMLInput(is,tag);
                  break;
               }
               if("AsymB0"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&mAsymBerarBaldinozziB0).XMLInput(is,tag);
                  break;
               }
               if("AsymB1"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&mAsymBerarBaldinozziB1).XMLInput(is,tag);
                  break;
               }
               #endif
            }
         }
         continue;
      }
      if("Option"==tag.GetName())
      {
         for(unsigned int i=0;i<tag.GetNbAttribute();i++)
            if("Name"==tag.GetAttributeName(i))
               mOptionRegistry.GetObj(tag.GetAttributeValue(i)).XMLInput(is,tag);
         continue;
      }
   }
}
#ifdef __WX__CRYST__
WXCrystObjBasic* ReflectionProfilePseudoVoigt::WXCreate(wxWindow* parent)
{
   VFN_DEBUG_ENTRY("ReflectionProfilePseudoVoigt::WXCreate()",6)
   if(mpWXCrystObj==0)
      mpWXCrystObj=new WXProfilePseudoVoigt(parent,this);
   VFN_DEBUG_EXIT("ReflectionProfilePseudoVoigt::WXCreate()",6)
   return mpWXCrystObj;
}
#endif

////////////////////////////////////////////////////////////////////////
//
//    ReflectionProfilePseudoVoigtAnisotropic
//
////////////////////////////////////////////////////////////////////////

ReflectionProfilePseudoVoigtAnisotropic::ReflectionProfilePseudoVoigtAnisotropic():
mCagliotiU(0),mCagliotiV(0),mCagliotiW(.01*DEG2RAD*DEG2RAD),mScherrerP(0),mLorentzX(0),mLorentzY(0),
mLorentzGammaHH(0),mLorentzGammaKK(0),mLorentzGammaLL(0),mLorentzGammaHK(0),mLorentzGammaHL(0),mLorentzGammaKL(0),
mPseudoVoigtEta0(0.5),mPseudoVoigtEta1(0),mAsym0(1.0),mAsym1(0),mAsym2(0)
{
   this->InitParameters();
}

ReflectionProfilePseudoVoigtAnisotropic::ReflectionProfilePseudoVoigtAnisotropic(const ReflectionProfilePseudoVoigtAnisotropic &old):
mCagliotiU(old.mCagliotiU),mCagliotiV(old.mCagliotiV),mCagliotiW(old.mCagliotiW),mScherrerP(old.mScherrerP),mLorentzX(old.mLorentzX),mLorentzY(old.mLorentzY),
mLorentzGammaHH(old.mLorentzGammaHH),mLorentzGammaKK(old.mLorentzGammaKK),mLorentzGammaLL(old.mLorentzGammaLL),mLorentzGammaHK(old.mLorentzGammaHK),mLorentzGammaHL(old.mLorentzGammaHL),mLorentzGammaKL(old.mLorentzGammaKL),
mPseudoVoigtEta0(old.mPseudoVoigtEta0),mPseudoVoigtEta1(old.mPseudoVoigtEta1),mAsym0(old.mAsym0),mAsym1(old.mAsym1),mAsym2(old.mAsym2)
{
   this->InitParameters();
}
ReflectionProfilePseudoVoigtAnisotropic::~ReflectionProfilePseudoVoigtAnisotropic()
{
   #ifdef __WX__CRYST__
   if(mpWXCrystObj!=0)
   {
      delete mpWXCrystObj;
      mpWXCrystObj=0;
   }
   #endif
}
   
ReflectionProfilePseudoVoigtAnisotropic*   ReflectionProfilePseudoVoigtAnisotropic::CreateCopy()const
{
   return new ReflectionProfilePseudoVoigtAnisotropic(*this);
}

const string& ReflectionProfilePseudoVoigtAnisotropic::GetClassName()const
{
   static string className="ReflectionProfilePseudoVoigtAnisotropic";
   return className;
}

CrystVector_REAL ReflectionProfilePseudoVoigtAnisotropic::GetProfile(const CrystVector_REAL &x, const REAL center,
                            const REAL h, const REAL k, const REAL l)const
{
   VFN_DEBUG_ENTRY("ReflectionProfilePseudoVoigtAnisotropic::GetProfile()",2)
   const REAL tantheta=tan(center/2.0);
   const REAL costheta=cos(center/2.0);
   const REAL sintheta=sin(center/2.0);
   const REAL fwhmG=sqrt(abs( mCagliotiW+mCagliotiV*tantheta+mCagliotiU*tantheta*tantheta+mScherrerP/(costheta*costheta)));
   const REAL gam=mLorentzGammaHH*h*h+mLorentzGammaKK*k*k+mLorentzGammaLL*l*l+2*mLorentzGammaHK*h*k+2*mLorentzGammaHL*h*l+2*mLorentzGammaKL*k*l;
   const REAL fwhmL= mLorentzX/costheta+(mLorentzY+gam/(sintheta*sintheta))*tantheta;
   // Eta for gaussian/lorentzian mix. Make sure 0<=eta<=1, else profiles could be <0 !
   REAL eta=mPseudoVoigtEta0+center*mPseudoVoigtEta1;
   if(eta>1) eta=1;
   if(eta<0) eta=0;

   CrystVector_REAL profile(x.numElements()),tmpV(x.numElements());
   const REAL asym=mAsym0+mAsym1/sin(center)+mAsym2/pow((REAL)sin(center),(REAL)2.0);
   VFN_DEBUG_MESSAGE("ReflectionProfilePseudoVoigtAnisotropic::GetProfile():("<<int(h)<<","<<int(k)<<","<<int(l)<<"),fwhmG="<<fwhmG<<",fwhmL="<<fwhmL<<",gam="<<gam<<",asym="<<asym<<",center="<<center<<",eta="<<eta, 2)
   if(fwhmG>0)
   {
      profile=PowderProfileGauss(x,fwhmG,center,asym);
      profile *= 1-eta;
   }
   else profile=0;
   if(fwhmL>0)
   {
      tmpV=PowderProfileLorentz(x,fwhmL,center,asym);
      tmpV *= eta;
      profile += tmpV;
   }
   VFN_DEBUG_MESSAGE(FormatVertVector<REAL>(x,profile),1)
   VFN_DEBUG_EXIT("ReflectionProfilePseudoVoigtAnisotropic::GetProfile()",2)
   return profile;
}

void ReflectionProfilePseudoVoigtAnisotropic::SetProfilePar(const REAL fwhmCagliotiW,
                   const REAL fwhmCagliotiU,
                   const REAL fwhmCagliotiV,
                   const REAL fwhmGaussP,
                   const REAL fwhmLorentzX,
                   const REAL fwhmLorentzY,
                   const REAL fwhmLorentzGammaHH,
                   const REAL fwhmLorentzGammaKK,
                   const REAL fwhmLorentzGammaLL,
                   const REAL fwhmLorentzGammaHK,
                   const REAL fwhmLorentzGammaHL,
                   const REAL fwhmLorentzGammaKL,
                   const REAL pseudoVoigtEta0,
                   const REAL pseudoVoigtEta1,
                   const REAL asymA0,
                   const REAL asymA1,
                   const REAL asymA2
                   )
{
   mCagliotiU=fwhmCagliotiU;
   mCagliotiV=fwhmCagliotiV;
   mCagliotiW=fwhmCagliotiW;
   mLorentzX=fwhmLorentzX;
   mLorentzY=fwhmLorentzY;
   mLorentzGammaHH=fwhmLorentzGammaHH;
   mLorentzGammaKK=fwhmLorentzGammaKK;
   mLorentzGammaLL=fwhmLorentzGammaLL;
   mLorentzGammaHK=fwhmLorentzGammaHK;
   mLorentzGammaHL=fwhmLorentzGammaHL;
   mLorentzGammaKL=fwhmLorentzGammaKL;
   mPseudoVoigtEta0=pseudoVoigtEta0;
   mPseudoVoigtEta1=pseudoVoigtEta1;
   mAsym0=asymA0;
   mAsym1=asymA1;
   mAsym2=asymA2;
   mClockMaster.Click();
}

REAL ReflectionProfilePseudoVoigtAnisotropic::GetFullProfileWidth(const REAL relativeIntensity, const REAL center,
                                 const REAL h, const REAL k, const REAL l)
{
   VFN_DEBUG_ENTRY("ReflectionProfilePseudoVoigt::GetFullProfileWidth()",2)
   const int nb=100;
   const int halfnb=nb/2;
   CrystVector_REAL x(nb);
   REAL n=5.0;
   const REAL tantheta=tan(center/2.0);
   const REAL costheta=cos(center/2.0);
   const REAL sintheta=sin(center/2.0);
   const REAL fwhmG=sqrt(abs( mCagliotiW+mCagliotiV*tantheta+mCagliotiU*tantheta*tantheta+mScherrerP/(costheta*costheta)));
   const REAL gam=mLorentzGammaHH*h*h+mLorentzGammaKK*k*k+mLorentzGammaLL*l*l+2*mLorentzGammaHK*h*k+2*mLorentzGammaHL*h*l+2*mLorentzGammaKL*k*l;
   const REAL fwhmL= mLorentzX/costheta+(mLorentzY+gam/(sintheta*sintheta))*tantheta;
   const REAL eta=mPseudoVoigtEta0+mPseudoVoigtEta1*center;
   // Obviously this is not the REAL FWHM, just a _very_ crude starting approximation
   REAL fwhm=fwhmL*eta+fwhmG*(1-eta);
   if(fwhm<=0) fwhm=1e-3;
   CrystVector_REAL prof;
   while(true)
   {
      //Create an X array with 100 elements reaching +/- n*FWHM/2
      REAL *p=x.data();
      const REAL tmp=fwhm*n/nb;
      for(int i=0;i<nb;i++) *p++ = tmp*(i-halfnb);
      x+=center;
      prof=this->GetProfile(x,center,h,k,l);
      const REAL max=prof.max();
      const REAL test=max*relativeIntensity;
      int n1=0,n2=0;
      if((prof(0)<test)&&(prof(nb-1)<test))
      {
         p=prof.data();
         while(*p<test){ p++; n1++;n2++;}
         n1--;
         while(*p>test){ p++; n2++;}
         VFN_DEBUG_EXIT("ReflectionProfilePseudoVoigtAnisotropic::GetFullProfileWidth():"<<x(n2)-x(n1),2)
         return x(n2)-x(n1);
      }
      VFN_DEBUG_MESSAGE("ReflectionProfilePseudoVoigtAnisotropic::GetFullProfileWidth():"<<relativeIntensity<<","
                        <<fwhm<<","<<center<<","<<h<<","<<k<<","<<l<<","<<max<<","<<test,2)
      VFN_DEBUG_MESSAGE(FormatVertVector<REAL>(x,prof),1)
      n*=2.0;
   }
}

bool ReflectionProfilePseudoVoigtAnisotropic::IsAnisotropic()const
{
   return true;
}

void ReflectionProfilePseudoVoigtAnisotropic::XMLOutput(ostream &os,int indent)const
{
   VFN_DEBUG_ENTRY("ReflectionProfilePseudoVoigtAnisotropic::XMLOutput():"<<this->GetName(),5)
   for(int i=0;i<indent;i++) os << "  " ;
   XMLCrystTag tag("ReflectionProfilePseudoVoigtAnisotropic");
   os <<tag<<endl;
   indent++;
   
   this->GetPar(&mCagliotiU).XMLOutput(os,"U",indent);
   os <<endl;
   
   this->GetPar(&mCagliotiV).XMLOutput(os,"V",indent);
   os <<endl;
   
   this->GetPar(&mCagliotiW).XMLOutput(os,"W",indent);
   os <<endl;

   this->GetPar(&mScherrerP).XMLOutput(os,"P",indent);
   os <<endl;

   this->GetPar(&mLorentzX).XMLOutput(os,"X",indent);
   os <<endl;

   this->GetPar(&mLorentzY).XMLOutput(os,"Y",indent);
   os <<endl;

   this->GetPar(&mLorentzGammaHH).XMLOutput(os,"G_HH",indent);
   os <<endl;
   
   this->GetPar(&mLorentzGammaKK).XMLOutput(os,"G_KK",indent);
   os <<endl;
   
   this->GetPar(&mLorentzGammaLL).XMLOutput(os,"G_LL",indent);
   os <<endl;
   
   this->GetPar(&mLorentzGammaHK).XMLOutput(os,"G_HK",indent);
   os <<endl;
   
   this->GetPar(&mLorentzGammaHL).XMLOutput(os,"G_HL",indent);
   os <<endl;
   
   this->GetPar(&mLorentzGammaKL).XMLOutput(os,"G_KL",indent);
   os <<endl;
   
   this->GetPar(&mPseudoVoigtEta0).XMLOutput(os,"Eta0",indent);
   os <<endl;
   
   this->GetPar(&mPseudoVoigtEta1).XMLOutput(os,"Eta1",indent);
   os <<endl;
   
   this->GetPar(&mAsym0).XMLOutput(os,"Asym0",indent);
   os <<endl;
   
   this->GetPar(&mAsym1).XMLOutput(os,"Asym1",indent);
   os <<endl;
   
   this->GetPar(&mAsym2).XMLOutput(os,"Asym2",indent);
   os <<endl;
   indent--;
   tag.SetIsEndTag(true);
   for(int i=0;i<indent;i++) os << "  " ;
   os <<tag<<endl;
   VFN_DEBUG_EXIT("ReflectionProfilePseudoVoigtAnisotropic::XMLOutput():"<<this->GetName(),5)
}

void ReflectionProfilePseudoVoigtAnisotropic::XMLInput(istream &is,const XMLCrystTag &tagg)
{
   VFN_DEBUG_ENTRY("ReflectionProfilePseudoVoigtAnisotropic::XMLInput():"<<this->GetName(),5)
   for(unsigned int i=0;i<tagg.GetNbAttribute();i++)
   {
      if("Name"==tagg.GetAttributeName(i)) this->SetName(tagg.GetAttributeValue(i));
   }
   while(true)
   {
      XMLCrystTag tag(is);
      if(("ReflectionProfilePseudoVoigtAnisotropic"==tag.GetName())&&tag.IsEndTag())
      {
         this->UpdateDisplay();
         VFN_DEBUG_EXIT("ReflectionProfilePseudoVoigtAnisotropic::Exit():"<<this->GetName(),5)
         return;
      }
      if("Par"==tag.GetName())
      {
         for(unsigned int i=0;i<tag.GetNbAttribute();i++)
         {
            if("Name"==tag.GetAttributeName(i))
            {
               if("U"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&mCagliotiU).XMLInput(is,tag);
                  break;
               }
               if("V"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&mCagliotiV).XMLInput(is,tag);
                  break;
               }
               if("W"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&mCagliotiW).XMLInput(is,tag);
                  break;
               }
               if("P"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&mScherrerP).XMLInput(is,tag);
                  break;
               }
               if("X"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&mLorentzX).XMLInput(is,tag);
                  break;
               }
               if("Y"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&mLorentzY).XMLInput(is,tag);
                  break;
               }
               if("G_HH"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&mLorentzGammaHH).XMLInput(is,tag);
                  break;
               }
               if("G_KK"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&mLorentzGammaKK).XMLInput(is,tag);
                  break;
               }
               if("G_LL"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&mLorentzGammaLL).XMLInput(is,tag);
                  break;
               }
               if("G_HK"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&mLorentzGammaHK).XMLInput(is,tag);
                  break;
               }
               if("G_HL"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&mLorentzGammaHL).XMLInput(is,tag);
                  break;
               }
               if("G_KL"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&mLorentzGammaKL).XMLInput(is,tag);
                  break;
               }
               if("Eta0"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&mPseudoVoigtEta0).XMLInput(is,tag);
                  break;
               }
               if("Eta1"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&mPseudoVoigtEta1).XMLInput(is,tag);
                  break;
               }
               if("Asym0"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&mAsym0).XMLInput(is,tag);
                  break;
               }
               if("Asym1"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&mAsym1).XMLInput(is,tag);
                  break;
               }
               if("Asym2"==tag.GetAttributeValue(i))
               {
                  this->GetPar(&mAsym2).XMLInput(is,tag);
                  break;
               }
            }
         }
         continue;
      }
      if("Option"==tag.GetName())
      {
         for(unsigned int i=0;i<tag.GetNbAttribute();i++)
            if("Name"==tag.GetAttributeName(i))
               mOptionRegistry.GetObj(tag.GetAttributeValue(i)).XMLInput(is,tag);
         continue;
      }
   }
}

void ReflectionProfilePseudoVoigtAnisotropic::InitParameters()
{
   {
      RefinablePar tmp("U",&mCagliotiU,-1/RAD2DEG/RAD2DEG,1./RAD2DEG/RAD2DEG,
                       gpRefParTypeScattDataProfileWidth,
                       REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,RAD2DEG*RAD2DEG);
      tmp.AssignClock(mClockMaster);
      tmp.SetDerivStep(1e-9);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("V",&mCagliotiV,-1/RAD2DEG/RAD2DEG,1./RAD2DEG/RAD2DEG,
                       gpRefParTypeScattDataProfileWidth,
                       REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,RAD2DEG*RAD2DEG);
      tmp.AssignClock(mClockMaster);
      tmp.SetDerivStep(1e-9);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("W",&mCagliotiW,0,1./RAD2DEG/RAD2DEG,
                       gpRefParTypeScattDataProfileWidth,
                       REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,RAD2DEG*RAD2DEG);
      tmp.AssignClock(mClockMaster);
      tmp.SetDerivStep(1e-9);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("P",&mScherrerP,-1./RAD2DEG/RAD2DEG,1./RAD2DEG/RAD2DEG,
                       gpRefParTypeScattDataProfileWidth,
                       REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,RAD2DEG*RAD2DEG);
      tmp.AssignClock(mClockMaster);
      tmp.SetDerivStep(1e-9);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("X",&mLorentzX,0,1./RAD2DEG,
                       gpRefParTypeScattDataProfileWidth,
                       REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,RAD2DEG);
      tmp.AssignClock(mClockMaster);
      tmp.SetDerivStep(1e-9);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("Y",&mLorentzY,-1./RAD2DEG,1./RAD2DEG,
                       gpRefParTypeScattDataProfileWidth,
                       REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,RAD2DEG);
      tmp.AssignClock(mClockMaster);
      tmp.SetDerivStep(1e-9);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("G_HH",&mLorentzGammaHH,-1./RAD2DEG,1./RAD2DEG,
                       gpRefParTypeScattDataProfileWidth,
                       REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,RAD2DEG);
      tmp.AssignClock(mClockMaster);
      tmp.SetDerivStep(1e-9);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("G_KK",&mLorentzGammaKK,-1./RAD2DEG,1./RAD2DEG,
                       gpRefParTypeScattDataProfileWidth,
                       REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,RAD2DEG);
      tmp.AssignClock(mClockMaster);
      tmp.SetDerivStep(1e-9);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("G_LL",&mLorentzGammaLL,-1./RAD2DEG,1./RAD2DEG,
                       gpRefParTypeScattDataProfileWidth,
                       REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,RAD2DEG);
      tmp.AssignClock(mClockMaster);
      tmp.SetDerivStep(1e-9);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("G_HK",&mLorentzGammaHK,-1./RAD2DEG,1./RAD2DEG,
                       gpRefParTypeScattDataProfileWidth,
                       REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,RAD2DEG);
      tmp.AssignClock(mClockMaster);
      tmp.SetDerivStep(1e-9);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("G_HL",&mLorentzGammaHL,-1./RAD2DEG,1./RAD2DEG,
                       gpRefParTypeScattDataProfileWidth,
                       REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,RAD2DEG);
      tmp.AssignClock(mClockMaster);
      tmp.SetDerivStep(1e-9);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("G_KL",&mLorentzGammaKL,-1./RAD2DEG,1./RAD2DEG,
                       gpRefParTypeScattDataProfileWidth,
                       REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,RAD2DEG);
      tmp.AssignClock(mClockMaster);
      tmp.SetDerivStep(1e-9);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("Eta0",&mPseudoVoigtEta0,0,1.,gpRefParTypeScattDataProfileType,
                       REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false);
      tmp.AssignClock(mClockMaster);
      tmp.SetDerivStep(1e-4);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("Eta1",&mPseudoVoigtEta1,-1,1.,gpRefParTypeScattDataProfileType,
                       REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false);
      tmp.AssignClock(mClockMaster);
      tmp.SetDerivStep(1e-4);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("Asym0",&mAsym0,0.01,10.0,gpRefParTypeScattDataProfileAsym,
                       REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false);
      tmp.AssignClock(mClockMaster);
      tmp.SetDerivStep(1e-4);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("Asym1",&mAsym1,-1.0,1.0,gpRefParTypeScattDataProfileAsym,
                       REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false);
      tmp.AssignClock(mClockMaster);
      tmp.SetDerivStep(1e-4);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("Asym2",&mAsym2,-1.0,1.0,gpRefParTypeScattDataProfileAsym,
                       REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false);
      tmp.AssignClock(mClockMaster);
      tmp.SetDerivStep(1e-4);
      this->AddPar(tmp);
   }
}

#ifdef __WX__CRYST__
WXCrystObjBasic* ReflectionProfilePseudoVoigtAnisotropic::WXCreate(wxWindow* parent)
{
   VFN_DEBUG_ENTRY("ReflectionProfilePseudoVoigt::WXCreate()",6)
   if(mpWXCrystObj==0)
      mpWXCrystObj=new WXProfilePseudoVoigtAnisotropic(parent,this);
   VFN_DEBUG_EXIT("ReflectionProfilePseudoVoigt::WXCreate()",6)
   return mpWXCrystObj;
}
#endif

////////////////////////////////////////////////////////////////////////
//
//    ReflectionProfileDoubleExponentialPseudoVoigt
//
////////////////////////////////////////////////////////////////////////
ReflectionProfileDoubleExponentialPseudoVoigt::ReflectionProfileDoubleExponentialPseudoVoigt():
ReflectionProfile(),
mInstrumentAlpha0(0.0),
mInstrumentAlpha1(0.0952),
mInstrumentBeta0(0.0239),
mInstrumentBeta1(0.0043),
mGaussianSigma0(0.0),
mGaussianSigma1(7.0),
mGaussianSigma2(0.0),
mLorentzianGamma0(0.0),
mLorentzianGamma1(0.0),
mLorentzianGamma2(0.414),
mpCell(0)
{
   VFN_DEBUG_MESSAGE("ReflectionProfileDoubleExponentialPseudoVoigt::ReflectionProfileDoubleExponentialPseudoVoigt()",10)
   this->InitParameters();
}

ReflectionProfileDoubleExponentialPseudoVoigt
   ::ReflectionProfileDoubleExponentialPseudoVoigt(const UnitCell &cell):
ReflectionProfile(),
mInstrumentAlpha0(0.0),
mInstrumentAlpha1(0.0952),
mInstrumentBeta0(0.0239),
mInstrumentBeta1(0.0043),
mGaussianSigma0(0.0),
mGaussianSigma1(7.0),
mGaussianSigma2(0.0),
mLorentzianGamma0(0.0),
mLorentzianGamma1(0.0),
mLorentzianGamma2(0.414),
mpCell(&cell)
{
   VFN_DEBUG_MESSAGE("ReflectionProfileDoubleExponentialPseudoVoigt::ReflectionProfileDoubleExponentialPseudoVoigt()",10)
   this->InitParameters();
}

ReflectionProfileDoubleExponentialPseudoVoigt::ReflectionProfileDoubleExponentialPseudoVoigt
   (const ReflectionProfileDoubleExponentialPseudoVoigt &old):
ReflectionProfile(),
mInstrumentAlpha0(old.mInstrumentAlpha0),
mInstrumentAlpha1(old.mInstrumentAlpha1),
mInstrumentBeta0(old.mInstrumentBeta0),
mInstrumentBeta1(old.mInstrumentBeta1),
mGaussianSigma0(old.mGaussianSigma0),
mGaussianSigma1(old.mGaussianSigma1),
mGaussianSigma2(old.mGaussianSigma2),
mLorentzianGamma0(old.mLorentzianGamma0),
mLorentzianGamma1(old.mLorentzianGamma1),
mLorentzianGamma2(old.mLorentzianGamma2),
mpCell(old.mpCell)
{
   VFN_DEBUG_MESSAGE("ReflectionProfileDoubleExponentialPseudoVoigt::ReflectionProfileDoubleExponentialPseudoVoigt()",10)
   this->InitParameters();
}

ReflectionProfileDoubleExponentialPseudoVoigt::~ReflectionProfileDoubleExponentialPseudoVoigt()
{
   #ifdef __WX__CRYST__
   if(mpWXCrystObj!=0)
   {
      delete mpWXCrystObj;
      mpWXCrystObj=0;
   }
   #endif
}

ReflectionProfileDoubleExponentialPseudoVoigt*
   ReflectionProfileDoubleExponentialPseudoVoigt::CreateCopy()const
{
   return new ReflectionProfileDoubleExponentialPseudoVoigt(*this);
}

const string& ReflectionProfileDoubleExponentialPseudoVoigt::GetClassName()const
{
   static string className="ReflectionProfileDoubleExponentialPseudoVoigt";
   return className;
}

CrystVector_REAL ReflectionProfileDoubleExponentialPseudoVoigt
   ::GetProfile(const CrystVector_REAL &x, const REAL center,
                const REAL h, const REAL k, const REAL l)const
{
   VFN_DEBUG_ENTRY("ReflectionProfileDoubleExponentialPseudoVoigt::GetProfile()",4)
   REAL dcenter=0;
   if(mpCell!=0)
   {
      REAL hh=h,kk=k,ll=l;// orthonormal coordinates in reciprocal space
      mpCell->MillerToOrthonormalCoords(hh,kk,ll);
      dcenter=1.0/sqrt(hh*hh+kk*kk+ll*ll);//d_hkl, in Angstroems
   }
   const REAL alpha=mInstrumentAlpha0+mInstrumentAlpha1/dcenter;
   const REAL beta=mInstrumentBeta0+mInstrumentBeta1/pow(dcenter,4);
   const REAL siggauss2= mGaussianSigma0
                        +mGaussianSigma1*pow(dcenter,2)
                        +mGaussianSigma2*pow(dcenter,4);
   static const REAL log2=log(2.0);
   const REAL hg=sqrt(8*siggauss2*log2);
   const REAL hl= mLorentzianGamma0
                 +mLorentzianGamma1*dcenter
                 +mLorentzianGamma2*dcenter*dcenter;
   const REAL hcom=pow(pow(hg,5)+2.69269*pow(hg,4)*hl+2.42843*pow(hg,3)*hl*hl
                       +4.47163*hg*hg*pow(hl,3)+0.07842*hg*pow(hl,4)+pow(hl,5),0.2);
   const REAL sigcom2=hcom*hcom/(8.0*log2);
   const REAL eta=1.36603*hl/hcom-0.47719*pow(hl/hcom,2)+0.11116*pow(hl/hcom,3);
   const long nbPoints=x.numElements();
   VFN_DEBUG_MESSAGE("ReflectionProfileDoubleExponentialPseudoVoigt::GetProfile():alpha="
                     <<alpha<<",beta="<<beta<<",siggauss2="<<siggauss2
                     <<",hg="<<hg<<",hl="<<hl<<",hcom="<<hcom<<",sigcom2="<<sigcom2
                     <<",eta="<<eta,2)
   CrystVector_REAL prof;
   prof=x;
   prof+=-center;
   REAL *pp=prof.data();
   for(long i=0;i<nbPoints;i++)
   {
      const double u=alpha/2*(alpha*sigcom2+2* *pp);
      const double nu=beta/2*(beta *sigcom2-2* *pp);
      const double y=(alpha*sigcom2+*pp)/sqrt(2*sigcom2);
      const double z=(beta *sigcom2-*pp)/sqrt(2*sigcom2);
      const complex<double> p(alpha* *pp,alpha*hcom/2);
      const complex<double> q(-beta* *pp, beta*hcom/2);
      const complex<double> e1p=ExponentialIntegral1_ExpZ(p);
      const complex<double> e1q=ExponentialIntegral1_ExpZ(q);
      VFN_DEBUG_MESSAGE("dt="<<*pp<<",  u="<<u<<",nu="<<nu<<",y="<<y<<",z="<<z
                        <<",p=("<<p.real()<<","<<p.imag()
                        <<"),q=("<<q.real()<<","<<q.imag()
                        <<"),e^p*E1(p)=("<<e1p.real()<<","<<e1p.imag()
                        <<"),e^q*E1(q)=("<<e1q.real()<<","<<e1q.imag(),2)
      REAL expnu_erfcz,expu_erfcy;
      // Use asymptotic value for erfc(x) = 1/(sqrt(pi)*x*exp(x^2)) [A&S 7.1.23]
      if(z>10.0) expnu_erfcz=exp(nu-z*z)/(z*sqrt(M_PI));
      else expnu_erfcz=exp(nu)*erfc(z);

      if(y>10.0) expu_erfcy=exp(u-y*y)/(y*sqrt(M_PI));
      else expu_erfcy=exp(u)*erfc(y);

      #if 0
      double tmp=(1-eta)*alpha*beta/(2*(alpha+beta))*(expu_erfcy+expnu_erfcz)
           -eta*alpha*beta/(M_PI*(alpha+beta))*(e1p.imag()+e1q.imag());
      if(isnan(*pp))// Is this portable ? Test for numeric_limits<REAL>::quiet_NaN()
      {
         cout<<"*pp==numeric_limits<REAL>::quiet_NaN()"<<endl;
         cout<<"ReflectionProfileDoubleExponentialPseudoVoigt::GetProfile():"<<endl
                     <<"   alpha="<<alpha<<",beta="<<beta<<",siggauss2="<<siggauss2
                     <<",hg="<<hg<<",hl="<<hl<<",hcom="<<hcom<<",sigcom2="<<sigcom2
                     <<",eta="<<eta<<endl;
         cout<<"   dt="<<*pp<<",  u="<<u<<",nu="<<nu<<",y="<<y<<",z="<<z
                        <<",e^u*E1(y)="<<expu_erfcy
                        <<",e^nu*E1(z)="<<expnu_erfcz
                        <<endl
                        <<"   p=("<<p.real()<<","<<p.imag()
                        <<"),q=("<<q.real()<<","<<q.imag()
                        <<"),e^p*E1(p)=("<<e1p.real()<<","<<e1p.imag()
                        <<"),e^q*E1(q)=("<<e1q.real()<<","<<e1q.imag()<<endl;
         cout<<(1-eta)*alpha*beta/(2*(alpha+beta))*(expu_erfcy+expnu_erfcz)<<endl
             <<eta*alpha*beta/(M_PI*(alpha+beta))*(e1p.imag()+e1q.imag())<<endl
             << *pp<<endl
             << tmp<<endl;
         exit(0);
      }
      if(abs(*pp)==numeric_limits<REAL>::infinity())
      {
         cout<<"*pp==numeric_limits<REAL>::infinity()"<<endl;
         exit(0);
      }
      //if(*pp>1e30) exit(0);
      #endif
      *pp++=(1-eta)*alpha*beta/(2*(alpha+beta))*(expu_erfcy+expnu_erfcz)
            -eta*alpha*beta/(M_PI*(alpha+beta))*(e1p.imag()+e1q.imag());
   }
   VFN_DEBUG_EXIT("ReflectionProfileDoubleExponentialPseudoVoigt::GetProfile()",4)
   return prof;
}

void ReflectionProfileDoubleExponentialPseudoVoigt
   ::SetProfilePar(const REAL instrumentAlpha0,
                   const REAL instrumentAlpha1,
                   const REAL instrumentBeta0,
                   const REAL instrumentBeta1,
                   const REAL gaussianSigma0,
                   const REAL gaussianSigma1,
                   const REAL gaussianSigma2,
                   const REAL lorentzianGamma0,
                   const REAL lorentzianGamma1,
                   const REAL lorentzianGamma2)
{
   mInstrumentAlpha0=instrumentAlpha0;
   mInstrumentAlpha1=instrumentAlpha1;
   mInstrumentBeta0=instrumentBeta0;
   mInstrumentBeta1=instrumentBeta1;
   mGaussianSigma0=gaussianSigma0;
   mGaussianSigma1=gaussianSigma1;
   mGaussianSigma2=gaussianSigma2;
   mLorentzianGamma0=lorentzianGamma0;
   mLorentzianGamma1=lorentzianGamma1;
   mLorentzianGamma2=lorentzianGamma2;
   mClockMaster.Click();
}

REAL ReflectionProfileDoubleExponentialPseudoVoigt
   ::GetFullProfileWidth(const REAL relativeIntensity, const REAL center,
                         const REAL h, const REAL k, const REAL l)
{
   VFN_DEBUG_ENTRY("ReflectionProfileDoubleExponentialPseudoVoigt::GetFullProfileWidth()",5)
   REAL dcenter=0;
   if(mpCell!=0)
   {
      REAL hh=h,kk=k,ll=l;// orthonormal coordinates in reciprocal space
      VFN_DEBUG_MESSAGE("ReflectionProfileDoubleExponentialPseudoVoigt::GetFullProfileWidth(),"<<dcenter<<","<<mpCell->GetName(),5)
      mpCell->MillerToOrthonormalCoords(hh,kk,ll);
      VFN_DEBUG_MESSAGE("ReflectionProfileDoubleExponentialPseudoVoigt::GetFullProfileWidth(),"<<dcenter,5)
      dcenter=sqrt(hh*hh+kk*kk+ll*ll);//1/d
   }
   VFN_DEBUG_MESSAGE("ReflectionProfileDoubleExponentialPseudoVoigt::GetFullProfileWidth(),"<<dcenter,5)
   const int nb=100;
   const int halfnb=nb/2;
   CrystVector_REAL x(nb);
   REAL n=5.0;
   const REAL siggauss2= mGaussianSigma0
                        +mGaussianSigma1*pow(dcenter,2)
                        +mGaussianSigma2*pow(dcenter,4);
   static const REAL log2=log(2.0);
   const REAL hg=sqrt(8*siggauss2*log2);
   const REAL hl= mLorentzianGamma0
                 +mLorentzianGamma1*dcenter
                 +mLorentzianGamma2*dcenter*dcenter;
   const REAL fwhm=pow(pow(hg,5)+2.69269*pow(hg,4)*hl+2.42843*pow(hg,3)*hl*hl
                       +4.47163*hg*hg*pow(hl,3)+0.07842*hg*pow(hl,4)+pow(hl,5),0.2);
   CrystVector_REAL prof;
   while(true)
   {
      REAL *p=x.data();
      const REAL tmp=fwhm*n/nb;
      for(int i=0;i<nb;i++) *p++ = tmp*(i-halfnb);
      x+=center;
      prof=this->GetProfile(x,center,h,k,l);
      const REAL max=prof.max();
      const REAL test=max*relativeIntensity;
      int n1=0,n2=0;
      if((prof(0)<test)&&(prof(nb-1)<test))
      {
         p=prof.data();
         while(*p<test){ p++; n1++;n2++;}
         n1--;
         while(*p>test){ p++; n2++;}
         VFN_DEBUG_EXIT("ReflectionProfilePseudoVoigt::GetFullProfileWidth():"<<x(n2)-x(n1),5)
         return abs(x(n2)-x(n1));
      }
      VFN_DEBUG_MESSAGE("ReflectionProfilePseudoVoigt::GetFullProfileWidth():"<<max<<","<<test
                        <<endl<<FormatVertVector<REAL>(x,prof),5)
      n*=2.0;
      //if(n>200) exit(0);
   }
   VFN_DEBUG_EXIT("ReflectionProfileDoubleExponentialPseudoVoigt::GetFullProfileWidth()",5)
}

bool ReflectionProfileDoubleExponentialPseudoVoigt
   ::IsAnisotropic()const{return false;}

void ReflectionProfileDoubleExponentialPseudoVoigt
   ::XMLOutput(ostream &os,int indent)const
{
   VFN_DEBUG_ENTRY("ReflectionProfileDoubleExponentialPseudoVoigt::XMLOutput():"<<this->GetName(),5)
   for(int i=0;i<indent;i++) os << "  " ;
   XMLCrystTag tag("ReflectionProfileDoubleExponentialPseudoVoigt");
   os <<tag<<endl;
   indent++;

   this->GetPar(&mInstrumentAlpha0).XMLOutput(os,"Alpha0",indent);
   os <<endl;

   this->GetPar(&mInstrumentAlpha1).XMLOutput(os,"Alpha1",indent);
   os <<endl;

   this->GetPar(&mInstrumentBeta0).XMLOutput(os,"Beta0",indent);
   os <<endl;

   this->GetPar(&mInstrumentBeta1).XMLOutput(os,"Beta1",indent);
   os <<endl;

   this->GetPar(&mGaussianSigma0).XMLOutput(os,"GaussianSigma0",indent);
   os <<endl;

   this->GetPar(&mGaussianSigma1).XMLOutput(os,"GaussianSigma1",indent);
   os <<endl;

   this->GetPar(&mGaussianSigma2).XMLOutput(os,"GaussianSigma2",indent);
   os <<endl;

   this->GetPar(&mLorentzianGamma0).XMLOutput(os,"LorentzianGamma0",indent);
   os <<endl;

   this->GetPar(&mLorentzianGamma1).XMLOutput(os,"LorentzianGamma1",indent);
   os <<endl;

   this->GetPar(&mLorentzianGamma2).XMLOutput(os,"LorentzianGamma2",indent);
   os <<endl;

   indent--;
   tag.SetIsEndTag(true);
   for(int i=0;i<indent;i++) os << "  " ;
   os <<tag<<endl;
   VFN_DEBUG_EXIT("ReflectionProfileDoubleExponentialPseudoVoigt::XMLOutput():"<<this->GetName(),5)
}

void ReflectionProfileDoubleExponentialPseudoVoigt
   ::XMLInput(istream &is,const XMLCrystTag &tagg)
{
   VFN_DEBUG_ENTRY("ReflectionProfileDoubleExponentialPseudoVoigt::XMLInput():"<<this->GetName(),5)
   for(unsigned int i=0;i<tagg.GetNbAttribute();i++)
   {
      if("Name"==tagg.GetAttributeName(i)) this->SetName(tagg.GetAttributeValue(i));
   }
   while(true)
   {
      XMLCrystTag tag(is);
      if(("ReflectionProfileDoubleExponentialPseudoVoigt"==tag.GetName())&&tag.IsEndTag())
      {
         this->UpdateDisplay();
         VFN_DEBUG_EXIT("ReflectionProfileDoubleExponentialPseudoVoigt::Exit():"<<this->GetName(),5)
         return;
      }
      if("Par"==tag.GetName())
      {
         for(unsigned int i=0;i<tag.GetNbAttribute();i++)
         {
            if("Name"==tag.GetAttributeName(i))
            {
               this->GetPar(tag.GetAttributeValue(i)).XMLInput(is,tag);
            }
         }
         continue;
      }
      if("Option"==tag.GetName())
      {
         for(unsigned int i=0;i<tag.GetNbAttribute();i++)
            if("Name"==tag.GetAttributeName(i))
               mOptionRegistry.GetObj(tag.GetAttributeValue(i)).XMLInput(is,tag);
         continue;
      }
   }
}

void ReflectionProfileDoubleExponentialPseudoVoigt::SetUnitCell(const UnitCell &cell)
{
   VFN_DEBUG_MESSAGE("ReflectionProfileDoubleExponentialPseudoVoigt::SetUnitCell()",10)
   mpCell=&cell;
}

void ReflectionProfileDoubleExponentialPseudoVoigt
   ::InitParameters()
{
   {
      RefinablePar tmp("Alpha0",&mInstrumentAlpha0,0,1e6,gpRefParTypeScattDataProfile,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false);
      tmp.AssignClock(mClockMaster);
      tmp.SetDerivStep(1e-4);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("Alpha1",&mInstrumentAlpha1,0,1e6,gpRefParTypeScattDataProfile,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false);
      tmp.AssignClock(mClockMaster);
      tmp.SetDerivStep(1e-6);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("Beta0",&mInstrumentBeta0,0,1e6,gpRefParTypeScattDataProfile,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false);
      tmp.AssignClock(mClockMaster);
      tmp.SetDerivStep(1e-6);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("Beta1",&mInstrumentBeta1,0,1e6,gpRefParTypeScattDataProfile,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false);
      tmp.AssignClock(mClockMaster);
      tmp.SetDerivStep(1e-6);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("GaussianSigma0",&mGaussianSigma0,0,1e6,gpRefParTypeScattDataProfileWidth,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false);
      tmp.AssignClock(mClockMaster);
      tmp.SetDerivStep(1e-4);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("GaussianSigma1",&mGaussianSigma1,0,1e6,gpRefParTypeScattDataProfileWidth,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false);
      tmp.AssignClock(mClockMaster);
      tmp.SetDerivStep(1e-4);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("GaussianSigma2",&mGaussianSigma2,0,1e6,gpRefParTypeScattDataProfileWidth,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false);
      tmp.AssignClock(mClockMaster);
      tmp.SetDerivStep(1e-4);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("LorentzianGamma0",&mLorentzianGamma0,0,1e6,gpRefParTypeScattDataProfileWidth,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false);
      tmp.AssignClock(mClockMaster);
      tmp.SetDerivStep(1e-4);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("LorentzianGamma1",&mLorentzianGamma1,0,1e6,gpRefParTypeScattDataProfileWidth,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false);
      tmp.AssignClock(mClockMaster);
      tmp.SetDerivStep(1e-4);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("LorentzianGamma2",&mLorentzianGamma2,0,1e6,gpRefParTypeScattDataProfileWidth,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false);
      tmp.AssignClock(mClockMaster);
      tmp.SetDerivStep(1e-4);
      this->AddPar(tmp);
   }
}

#ifdef __WX__CRYST__
WXCrystObjBasic* ReflectionProfileDoubleExponentialPseudoVoigt::WXCreate(wxWindow* parent)
{
   if(mpWXCrystObj==0)
      mpWXCrystObj=new WXProfileDoubleExponentialPseudoVoigt(parent,this);
   return mpWXCrystObj;
}
#endif

//######################################################################
//    Basic PROFILE FUNCTIONS
//######################################################################

CrystVector_REAL PowderProfileGauss  (const CrystVector_REAL ttheta,const REAL fw,
                                      const REAL center, const REAL asym)
{
   TAU_PROFILE("PowderProfileGauss()","Vector (Vector,REAL)",TAU_DEFAULT);
   REAL fwhm=fw;
   if(fwhm<=0) fwhm=1e-6;
   const long nbPoints=ttheta.numElements();
   CrystVector_REAL result(nbPoints);
   result=ttheta;
   result+= -center;
   result *= result;
   REAL *p;
   if(false)// fabs(asym-1.) < 1e-5)
   {
      //reference: IUCr Monographs on Crystallo 5 - The Rietveld Method (ed RA Young)
      result *= -4.*log(2.)/fwhm/fwhm;
   }
   else
   {  // Adapted from Toraya J. Appl. Cryst 23(1990),485-491
      const REAL c1= -(1.+asym)/asym*(1.+asym)/asym*log(2.)/fwhm/fwhm;
      const REAL c2= -(1.+asym)     *(1.+asym)     *log(2.)/fwhm/fwhm;
      long i;
      p=result.data();
      const REAL *pt=ttheta.data();
      for(i=0;i<nbPoints;i++){ *p++ *= c1;if(*pt++>center) break;}
      i++;
      for(   ;i<nbPoints;i++)  *p++ *= c2;
   }
   p=result.data();
   #ifdef _MSC_VER
   // Bug from Hell (in MSVC++) !
   // The *last* point ends up sometimes with an arbitrary large value...
   for(long i=0;i<nbPoints;i++) { *p = pow((float)2.71828182846,(float)*p) ; p++ ;}
   #else
   long i=nbPoints;
   for(;i>3;i-=4)
   {
     #ifdef HAVE_SSE_MATHFUN
     v4sf x=_mm_loadu_ps(p);
     _mm_storeu_ps(p,exp_ps(x));
     p+=4;
     #else
     for(unsigned int j=0;j<4;++j)
     {// Fixed-length loop enables vectorization
       *p = exp(*p) ;
       p++ ;
     }
     #endif
   }
   for(;i>0;i--) { *p = exp(*p) ; p++ ;}
   #endif

#if 0
   #if 1 //def _MSC_VER
   // Bug from Hell (in MSVC++) !
   // The *last* point ends up sometimes with an arbitrary large value...
   long i=0;
   for(;i<nbPoints;i+=4)
      for(unsigned int j=0;j<4;++j)
      {// Fixed-length loop enables vectorization
        *p = pow((float)2.71828182846,(float)*p) ;
        p++ ;
      }
   #else
   long i=0;
   for(;i<nbPoints;i+=1)
   {
     //for(unsigned int j=0;j<4;++j)
     {// Fixed-length loop enables vectorization
       *p = exp(*p) ;
       p++ ;
     }
   }
   #endif
#endif
   result *= 2. / fwhm * sqrt(log(2.)/M_PI);
   return result;
}

CrystVector_REAL PowderProfileLorentz(const CrystVector_REAL ttheta,const REAL fw,
                                      const REAL center, const REAL asym)
{
   TAU_PROFILE("PowderProfileLorentz()","Vector (Vector,REAL)",TAU_DEFAULT);
   REAL fwhm=fw;
   if(fwhm<=0) fwhm=1e-6;
   const long nbPoints=ttheta.numElements();
   CrystVector_REAL result(nbPoints);
   result=ttheta;
   result+= -center;
   result *= result;
   REAL *p;
   if(false)// fabs(asym-1.) < 1e-5)
   {
      //reference: IUCr Monographs on Crystallo 5 - The Rietveld Method (ed RA Young)
      result *= 4./fwhm/fwhm;
   }
   else
   {  // Adapted from Toraya J. Appl. Cryst 23(1990),485-491
      const REAL c1= (1+asym)/asym*(1+asym)/asym/fwhm/fwhm;
      const REAL c2= (1+asym)     *(1+asym)     /fwhm/fwhm;
      long i;
      p=result.data();
      const REAL *pt=ttheta.data();
      for(i=0;i<nbPoints;i++){ *p++ *= c1;if(*pt++>center) break;}
      i++;
      for(   ;i<nbPoints;i++)  *p++ *= c2 ;
   }
   p=result.data();
   result += 1. ;
   for(long i=0;i<nbPoints;i++) { *p = 1/(*p) ; p++ ;}
   result *= 2./M_PI/fwhm;
   return result;
}

CrystVector_REAL AsymmetryBerarBaldinozzi(const CrystVector_REAL x,
                                          const REAL fw, const REAL center,
                                          const REAL a0, const REAL a1,
                                          const REAL b0, const REAL b1)
{
   TAU_PROFILE("AsymmetryBerarBaldinozzi()","Vector (Vector,REAL)",TAU_DEFAULT);
   REAL fwhm=fw;
   if(fwhm<=0) fwhm=1e-6;
   const long nbPoints=x.numElements();
   CrystVector_REAL result(nbPoints);
   result=x;
   result+= -center;
   result *= 1/fwhm;
   REAL *p=result.data();
   const REAL a=a0/tan(center/2)+a1/tan(center);
   const REAL b=b0/tan(center/2)+b1/tan(center);
   for(long i=0;i<nbPoints;i++)
   {
      *p = 1+*p * exp(-*p * *p)*(2*a+b*(8* *p * *p-12));
      p++ ;
   }
   return result;
}
/*
from python:
E1(1)= 0.219383934396       (0.219383934396+0j)
E1(1j)= (-0.337403922901-0.624713256428j)
E1(1+1j)= (0.000281624451981-0.179324535039j)
E1(100+1j)= (1.95936883899e-46-3.11904399563e-46j)
E1(10+20j)= (-1.20141500252e-06-1.58298052926e-06j)

this code (REAL=float)
CE1(1.000000000000+0.000000000000j) = 0.219383955002+0.000000000000j
CE1(0.000000000000+1.000000000000j) = -0.337403953075+-0.624713361263j
CE1(1.000000000000+1.000000000000j) = 0.000281602144+-0.179324567318j
CE1(100.000000000000+1.000000000000j) = 0.000000000000+-0.000000000000j
CE1(10.000000000000+20.000000000000j) = -0.000001201415+-0.000001582981j
   {
      complex<REAL>z(1.0,0.0);
      complex<REAL>ce1=ExponentialIntegral1(z);
      cout<<"CE1("<<z.real()<<"+"<<z.imag()<<"j) = "<<ce1.real()<<"+"<<ce1.imag()<<"j"<<endl;
   }
   {
      complex<REAL>z(0.0,1.0);
      complex<REAL>ce1=ExponentialIntegral1(z);
      cout<<"CE1("<<z.real()<<"+"<<z.imag()<<"j) = "<<ce1.real()<<"+"<<ce1.imag()<<"j"<<endl;
   }
   {
      complex<REAL>z(1.0,1.0);
      complex<REAL>ce1=ExponentialIntegral1(z);
      cout<<"CE1("<<z.real()<<"+"<<z.imag()<<"j) = "<<ce1.real()<<"+"<<ce1.imag()<<"j"<<endl;
   }
   {
      complex<REAL>z(100.0,1.0);
      complex<REAL>ce1=ExponentialIntegral1(z);
      cout<<"CE1("<<z.real()<<"+"<<z.imag()<<"j) = "<<ce1.real()<<"+"<<ce1.imag()<<"j"<<endl;
   }
   {
      complex<REAL>z(10.0,20.0);
      complex<REAL>ce1=ExponentialIntegral1(z);
      cout<<"CE1("<<z.real()<<"+"<<z.imag()<<"j) = "<<ce1.real()<<"+"<<ce1.imag()<<"j"<<endl;
   }
   exit(0);
*/

template <class T>std::complex<T>ExponentialIntegral1(const complex<T> z)
{
   return exp(-z)*ExponentialIntegral1_ExpZ(z);
}

template <class T>std::complex<T>ExponentialIntegral1_ExpZ(const complex<T> z)
{
   const T zr=z.real();
   const T zn=abs(z);
   complex<T> ce1;
   if(zn==0.0) return 1e100;// Should return an error ? std::numeric_limits::quiet_NaN() ?
   if((zn<10.0)||((zr<0.0)&&(zn<20.0)))// Abramowitz & Stegun 5.1.11
   {
      ce1=complex<T>(1,0);
      complex<T> y(1,0);
      for(unsigned int i=1;i<=150;i++)
      {
         y=-y*(T)i*z / (T)((i+1)*(i+1));
         ce1+=y;
         if(abs(y)<=abs(ce1)*1e-15) break;
      }
      static const T EulerMascheroni=0.5772156649015328606065120900;
      return exp(z)*(z*ce1-EulerMascheroni-log(z));// Euler-Mascheroni constant
   }
   else// Abramowitz & Stegun 5.1.51
   {
      if(zn>500) return 1.0/z;
      complex<T> y(0.0,0.0);
      for(unsigned int i=120;i>=1;i--) y=(T)i/((T)1+(T)i/(z+y));
      ce1/=(z+y);
      if((zr<0)&&(z.imag()==0)) ce1 -= complex<T>(0.0,M_PI)*exp(z);
      return ce1;
   }
}

}
