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
#include "ObjCryst/ScatteringCorr.h"
#include "Quirks/VFNStreamFormat.h"
#include <stdio.h> //for sprintf()
#ifdef __WX__CRYST__
namespace ObjCryst
{
   class TexturePhaseMarchDollase;
}
#include "wxCryst/wxPowderPattern.h"
#endif
namespace ObjCryst
{
////////////////////////////////////////////////////////////////////////
//
//        ScatteringCorr
//
////////////////////////////////////////////////////////////////////////

ScatteringCorr::ScatteringCorr(const ScatteringData & data):
mpData(&data)
{
   VFN_DEBUG_MESSAGE("ScatteringCorr::ScatteringCorr(&scattData)",5)
}

ScatteringCorr::~ScatteringCorr()
{
   VFN_DEBUG_MESSAGE("ScatteringCorr::~ScatteringCorr()",5)
}

const CrystVector_REAL& ScatteringCorr::GetCorr() const
{
   this->CalcCorr();
   return mCorr;
}

const RefinableObjClock& ScatteringCorr::GetClockCorr()const {return mClockCorrCalc;}
////////////////////////////////////////////////////////////////////////
//
//        LorentzCorr
//
////////////////////////////////////////////////////////////////////////
LorentzCorr::LorentzCorr(const ScatteringData & data):
ScatteringCorr(data)
{}

LorentzCorr::~LorentzCorr()
{}

const string & LorentzCorr::GetName() const
{
   //So far, we do not need a personalized name...
   const static string mName="LorentzCorr";
   return mName;
}

const string & LorentzCorr::GetClassName() const
{
   const static string className="LorentzCorr";
   return className;
}

void LorentzCorr::CalcCorr() const
{
   const CrystVector_REAL *theta=&(mpData->GetTheta());
   if(mpData->GetClockTheta()<mClockCorrCalc) return;
   TAU_PROFILE("LorentzCorr::CalcCorr()","void ()",TAU_DEFAULT);
   mCorr.resize(mpData->GetNbRefl());
   for(long i=0;i<mpData->GetNbRefl();i++)mCorr(i) =1/sin(2*(*theta)(i));
   mClockCorrCalc.Click();
}

////////////////////////////////////////////////////////////////////////
//
//        PolarizationCorr
//
////////////////////////////////////////////////////////////////////////
PolarizationCorr::PolarizationCorr(const ScatteringData & data):
ScatteringCorr(data),mPolarAfactor(1)
{}

PolarizationCorr::~PolarizationCorr()
{}

const string & PolarizationCorr::GetName() const
{
   //So far, we do not need a personalized name...
   const static string mName="PolarizationCorr";
   return mName;
}

const string & PolarizationCorr::GetClassName() const
{
   const static string className="PolarizationCorr";
   return className;
}

void PolarizationCorr::CalcCorr() const
{
   const CrystVector_REAL *theta=&(mpData->GetTheta());
   const REAL f=mpData->GetRadiation().GetLinearPolarRate();
   if(    (mpData->GetClockTheta()<mClockCorrCalc)
       && (fabs(mPolarAfactor-((1-f)/(1+f))) <(mPolarAfactor*.0001)) ) return;
   VFN_DEBUG_MESSAGE("PolarizationCorr::CalcCorr()",10)
   TAU_PROFILE("PolarizationCorr::CalcCorr()","void ()",TAU_DEFAULT);
   mPolarAfactor=((1-f)/(1+f));
   mCorr.resize(mpData->GetNbRefl());
   for(long i=0;i<mpData->GetNbRefl();i++)
      mCorr(i) =(1.+mPolarAfactor*cos((*theta)(i))*cos((*theta)(i)))/(1.+mPolarAfactor);
   mClockCorrCalc.Click();
}


////////////////////////////////////////////////////////////////////////
//
//        PowderSlitApertureCorr
//
////////////////////////////////////////////////////////////////////////
PowderSlitApertureCorr::PowderSlitApertureCorr(const ScatteringData & data):
ScatteringCorr(data)
{}

PowderSlitApertureCorr::~PowderSlitApertureCorr()
{}

const string & PowderSlitApertureCorr::GetName() const
{
   //So far, we do not need a personalized name...
   const static string mName="PowderSlitApertureCorr";
   return mName;
}

const string & PowderSlitApertureCorr::GetClassName() const
{
   const static string className="PowderSlitApertureCorr";
   return className;
}

void PowderSlitApertureCorr::CalcCorr() const
{
   const CrystVector_REAL *theta=&(mpData->GetTheta());
   if(mpData->GetClockTheta()<mClockCorrCalc) return;
   TAU_PROFILE("PowderSlitApertureCorr::CalcCorr()","void ()",TAU_DEFAULT);
   mCorr.resize(mpData->GetNbRefl());
   for(long i=0;i<mpData->GetNbRefl();i++)mCorr(i) =1/sin((*theta)(i));
   mClockCorrCalc.Click();
}
////////////////////////////////////////////////////////////////////////
//
//        TexturePhaseMarchDollase
//
////////////////////////////////////////////////////////////////////////
TexturePhaseMarchDollase::TexturePhaseMarchDollase(const REAL f, 
                                                   const REAL c,
                                                   const REAL h,
                                                   const REAL k, 
                                                   const REAL l,
                                                   TextureMarchDollase &tex):
mFraction(f),mMarchCoeff(c),mH(h),mK(k),mL(l),mpTextureMarchDollase(&tex)
#ifdef __WX__CRYST__
,mpWXCrystObj(0)
#endif
{}

TexturePhaseMarchDollase::~TexturePhaseMarchDollase()
{
   #ifdef __WX__CRYST__
   this->WXDelete();
   #endif
}
const string& TexturePhaseMarchDollase::GetClassName()const
{
   const static string className="March-Dollase Texture Phase";
   return className;
}
const string& TexturePhaseMarchDollase::GetName()const
{
   const static string name="March-Dollase Texture Phase";
   return name;
}

void TexturePhaseMarchDollase::SetPar(const REAL f, const REAL c,const REAL h,const REAL k, const REAL l)
{mFraction=f;mMarchCoeff=c;mH=h;mK=k;mL=l;}
void TexturePhaseMarchDollase::XMLOutput(ostream &os,int indent)const
{
   VFN_DEBUG_ENTRY("TexturePhaseMarchDollase::XMLOutput():"<<this->GetName(),5)
   for(int i=0;i<indent;i++) os << "  " ;
   XMLCrystTag tag("TexturePhaseMarchDollase");
   os <<tag<<endl;
   indent++;
   
   mpTextureMarchDollase->GetPar(&mFraction).XMLOutput(os,"Fraction",indent);
   os <<endl;
   
   mpTextureMarchDollase->GetPar(&mMarchCoeff).XMLOutput(os,"MarchCoeff",indent);
   os <<endl;
   
   mpTextureMarchDollase->GetPar(&mH).XMLOutput(os,"H",indent);
   os <<endl;
   
   mpTextureMarchDollase->GetPar(&mK).XMLOutput(os,"K",indent);
   os <<endl;
   
   mpTextureMarchDollase->GetPar(&mL).XMLOutput(os,"L",indent);
   os <<endl;
   
   indent--;
   tag.SetIsEndTag(true);
   for(int i=0;i<indent;i++) os << "  " ;
   os <<tag<<endl;
   VFN_DEBUG_EXIT("TexturePhaseMarchDollase::XMLOutput():"<<this->GetName(),5)
}

void TexturePhaseMarchDollase::XMLInput(istream &is,const XMLCrystTag &tagg)
{
   VFN_DEBUG_ENTRY("TexturePhaseMarchDollase::XMLInput():",5)
   for(unsigned int i=0;i<tagg.GetNbAttribute();i++)
   {
      //No attribute to read
   }
   while(true)
   {
      XMLCrystTag tag(is);
      if(("TexturePhaseMarchDollase"==tag.GetName())&&tag.IsEndTag())
      {
         VFN_DEBUG_EXIT("TexturePhaseMarchDollase::XMLInput()",5)
         return;
      }
      if("Par"==tag.GetName())
      {
         for(unsigned int i=0;i<tag.GetNbAttribute();i++)
         {
            if("Name"==tag.GetAttributeName(i))
            {
               if("Fraction"==tag.GetAttributeValue(i))
               {
                  mpTextureMarchDollase->GetPar(&mFraction).XMLInput(is,tag);
                  break;
               }
               if("MarchCoeff"==tag.GetAttributeValue(i))
               {
                  mpTextureMarchDollase->GetPar(&mMarchCoeff).XMLInput(is,tag);
                  break;
               }
               if("H"==tag.GetAttributeValue(i))
               {
                  mpTextureMarchDollase->GetPar(&mH).XMLInput(is,tag);
                  break;
               }
               if("K"==tag.GetAttributeValue(i))
               {
                  mpTextureMarchDollase->GetPar(&mK).XMLInput(is,tag);
                  break;
               }
               if("L"==tag.GetAttributeValue(i))
               {
                  mpTextureMarchDollase->GetPar(&mL).XMLInput(is,tag);
                  break;
               }
            }
         }
         continue;
      }
   }
}


#ifdef __WX__CRYST__
WXCrystObjBasic* TexturePhaseMarchDollase::WXCreate(wxWindow* parent)
{
   mpWXCrystObj=new WXTexturePhaseMarchDollase(parent,this,mpTextureMarchDollase);
   return mpWXCrystObj;
}
WXCrystObjBasic* TexturePhaseMarchDollase::WXGet()
{
   return mpWXCrystObj;
}

void TexturePhaseMarchDollase::WXDelete()
{
   if(0!=mpWXCrystObj)
   {
      VFN_DEBUG_MESSAGE("TexturePhaseMarchDollase::WXDelete()",5)
      delete mpWXCrystObj;
   }
}
void TexturePhaseMarchDollase::WXNotifyDelete(){mpWXCrystObj=0;}
#endif
////////////////////////////////////////////////////////////////////////
//
//        TextureMarchDollase
//
////////////////////////////////////////////////////////////////////////
TextureMarchDollase::TextureMarchDollase(const ScatteringData & data):
ScatteringCorr(data),mNbReflUsed(0)
{}

TextureMarchDollase::~TextureMarchDollase()
{
}

const string & TextureMarchDollase::GetName() const
{
   //So far, we do not need a personalized name...
   const static string name="TextureMarchDollase";
   return name;
}

const string & TextureMarchDollase::GetClassName() const
{
   //So far, we do not need a personalized name...
   const static string name="TextureMarchDollase";
   return name;
}

void TextureMarchDollase::AddPhase(const REAL f, const REAL c,
                                   const REAL h,const REAL k, const REAL l)
 
{
   VFN_DEBUG_ENTRY("TextureMarchDollase::AddPhase()",5)
   TexturePhaseMarchDollase* phase=new TexturePhaseMarchDollase(f,c,h,k,l,*this);
   this->Print();
   //Add parameters
   const unsigned int nbPhase=this->GetNbPhase();
   char buf [5];
   sprintf(buf,"%d",nbPhase);
   {
      RefinablePar tmp("Fraction_"+(string)buf,&(phase->mFraction),0.,1.,
                        gpRefParTypeScattDataCorrInt,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,1.);
      tmp.AssignClock(mClockTexturePar);
      tmp.SetDerivStep(1e-7);
      tmp.SetGlobalOptimStep(.05);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("MarchCoeff_"+(string)buf,&(phase->mMarchCoeff),.1,10.,
                        gpRefParTypeScattDataCorrInt,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,1.);
      tmp.AssignClock(mClockTexturePar);
      tmp.SetDerivStep(1e-7);
      tmp.SetGlobalOptimStep(.1);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("H_"+(string)buf,&(phase->mH),-10.,10.,
                        gpRefParTypeScattDataCorrInt,
                        REFPAR_DERIV_STEP_ABSOLUTE,false,true,true,false,1.);
      tmp.AssignClock(mClockTexturePar);
      tmp.SetDerivStep(1e-7);
      tmp.SetGlobalOptimStep(.01);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("K_"+(string)buf,&(phase->mK),-10.,10.,
                        gpRefParTypeScattDataCorrInt,
                        REFPAR_DERIV_STEP_ABSOLUTE,false,true,true,false,1.);
      tmp.AssignClock(mClockTexturePar);
      tmp.SetDerivStep(1e-7);
      tmp.SetGlobalOptimStep(.01);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("L_"+(string)buf,&(phase->mL),-10.,10.,
                        gpRefParTypeScattDataCorrInt,
                        REFPAR_DERIV_STEP_ABSOLUTE,false,true,true,false,1.);
      tmp.AssignClock(mClockTexturePar);
      tmp.SetDerivStep(1e-7);
      tmp.SetGlobalOptimStep(.01);
      this->AddPar(tmp);
   }
   this->Print();

   mPhaseRegistry.Register(*phase);
   VFN_DEBUG_EXIT("TextureMarchDollase::AddPhase()=",5)
}

void TextureMarchDollase::SetPhasePar(const unsigned int i, const REAL f, const REAL c,
                 const REAL h,const REAL k, const REAL l)
{
   mPhaseRegistry.GetObj(i).SetPar(f,c,h,k,l);
}

//void DeletePhase(const unsigned int i);

unsigned int TextureMarchDollase::GetNbPhase() const {return mPhaseRegistry.GetNb();}

REAL TextureMarchDollase::GetFraction(const unsigned int i)const
{return mPhaseRegistry.GetObj(i).mFraction;}

REAL TextureMarchDollase::GetMarchCoeff(const unsigned int i)const
{return mPhaseRegistry.GetObj(i).mMarchCoeff;}

REAL TextureMarchDollase::GetPhaseH(const unsigned int i)const
{return mPhaseRegistry.GetObj(i).mH;}

REAL TextureMarchDollase::GetPhaseK(const unsigned int i)const
{return mPhaseRegistry.GetObj(i).mK;}

REAL TextureMarchDollase::GetPhaseL(const unsigned int i)const
{return mPhaseRegistry.GetObj(i).mL;}

void TextureMarchDollase::GlobalOptRandomMove(const REAL mutationAmplitude,
                                              const RefParType *type)
{
   if(mRandomMoveIsDone) return;
   if(!(gpRefParTypeScattDataCorrInt->IsDescendantFromOrSameAs(type)))
   {
      mRandomMoveIsDone=true;
      return;
   }
   //if((rand()/(REAL)RAND_MAX)<.3)//only 30% proba to make a random move
   {
      VFN_DEBUG_MESSAGE("TextureMarchDollase::GlobalOptRandomMove()",1)
      for(unsigned int i=0;i<this->GetNbPhase();i++)
      {
         // :TODO: Give some probability (1% ?) to invert the March coefficient
         // with a direction perpendicular to the current one ?
         
         RefinablePar *pF=&(this->GetPar(&(mPhaseRegistry.GetObj(i).mFraction)));
         RefinablePar *pM=&(this->GetPar(&(mPhaseRegistry.GetObj(i).mMarchCoeff)));
         RefinablePar *pH=&(this->GetPar(&(mPhaseRegistry.GetObj(i).mH)));
         RefinablePar *pK=&(this->GetPar(&(mPhaseRegistry.GetObj(i).mK)));
         RefinablePar *pL=&(this->GetPar(&(mPhaseRegistry.GetObj(i).mL)));
         if(pF->IsFixed()==false)
         {
            const REAL delta=pF->GetGlobalOptimStep()*mutationAmplitude;
            const REAL sig=4*delta;
            const REAL y0=mPhaseRegistry.GetObj(i).mBiasFraction;
            REAL y,ymin,ymax;
            y=pF->GetValue();
            
            ymax=.5+1/M_PI*atan((y+delta-y0)/(2.*sig));
            ymin=.5+1/M_PI*atan((y-delta-y0)/(2.*sig));
            y=ymin+rand()/(REAL)RAND_MAX*(ymax-ymin);
            y-=.5;
            if(y<-.499)y=-.499;//Should not happen but make sure we remain in [-pi/2;pi/2]
            if(y> .499)y= .499;
            pF->MutateTo(y0+2*sig*tan(M_PI*y));
         }
         if((pH->IsFixed()==false)||(pK->IsFixed()==false)||(pL->IsFixed()==false))
         {
            REAL tx=pH->GetValue();
            REAL ty=pK->GetValue();
            REAL tz=pL->GetValue();
            mpData->GetCrystal().MillerToOrthonormalCoords(tx,ty,tz);
            {
               REAL tx0=mPhaseRegistry.GetObj(i).mBiasH;
               REAL ty0=mPhaseRegistry.GetObj(i).mBiasK;
               REAL tz0=mPhaseRegistry.GetObj(i).mBiasL;
               mpData->GetCrystal().MillerToOrthonormalCoords(tx0,ty0,tz0);
               const REAL delta=pH->GetGlobalOptimStep()*mutationAmplitude*mPhaseRegistry.GetObj(i).mNorm;
               const REAL sig=2*delta;
               REAL y,ymin,ymax;

               ymax=.5+1/M_PI*atan((tx+delta-tx0)/(2.*sig));
               ymin=.5+1/M_PI*atan((tx-delta-tx0)/(2.*sig));
               y=ymin+rand()/(REAL)RAND_MAX*(ymax-ymin);
               y-=.5;
               if(y<-.499)y=-.499;
               if(y> .499)y= .499;
               tx=tx0+2*sig*tan(M_PI*y);

               ymax=.5+1/M_PI*atan((ty+delta-ty0)/(2.*sig));
               ymin=.5+1/M_PI*atan((ty-delta-ty0)/(2.*sig));
               y=ymin+rand()/(REAL)RAND_MAX*(ymax-ymin);
               y-=.5;
               if(y<-.499)y=-.499;
               if(y> .499)y= .499;
               ty=ty0+2*sig*tan(M_PI*y);

               ymax=.5+1/M_PI*atan((tz+delta-tz0)/(2.*sig));
               ymin=.5+1/M_PI*atan((tz-delta-tz0)/(2.*sig));
               y=ymin+rand()/(REAL)RAND_MAX*(ymax-ymin);
               y-=.5;
               if(y<-.499)y=-.499;
               if(y> .499)y= .499;
               tz=tz0+2*sig*tan(M_PI*y);
            }
            const REAL factor=mPhaseRegistry.GetObj(i).mNorm/sqrt(tx*tx+ty*ty+tz*tz);
            tx *= factor;
            ty *= factor;
            tz *= factor;
            mpData->GetCrystal().OrthonormalToMillerCoords(tx,ty,tz);
            pH->MutateTo(tx);
            pK->MutateTo(ty);
            pL->MutateTo(tz);
         }
         if(pM->IsFixed()==false)
         {
            // Given the nature of this param, we use a proportionnal max step
            const REAL delta=pM->GetGlobalOptimStep()*mutationAmplitude;
            const REAL sig=2*delta;
            REAL y,ymin,ymax;
            y=log(pM->GetValue());
            const REAL y0=log(mPhaseRegistry.GetObj(i).mBiasMarchCoeff);

            ymin=.5+1/M_PI*atan((y-delta-y0)/(2.*sig));
            ymax=.5+1/M_PI*atan((y+delta-y0)/(2.*sig));
            y=ymin+rand()/(REAL)RAND_MAX*(ymax-ymin);
               y-=.5;
               if(y<-.499)y=-.499;
               if(y> .499)y= .499;
            pM->MutateTo(exp(y0+2*sig*tan(M_PI*y)));
         }
      }
   }
   //this->RefinableObj::Print();
   mRandomMoveIsDone=true;
}
REAL TextureMarchDollase::GetBiasingCost()const
{
   REAL cost=0;
   REAL tmp;
   for(unsigned int i=0; i<this->GetNbPhase();i++)
   {
      tmp =(mPhaseRegistry.GetObj(i).mBiasFraction-mPhaseRegistry.GetObj(i).mFraction)/.04;
      cost += tmp*tmp;
      
      tmp =log10(mPhaseRegistry.GetObj(i).mBiasMarchCoeff/mPhaseRegistry.GetObj(i).mMarchCoeff)/.04;
      cost += tmp*tmp;
      
      REAL tx=mPhaseRegistry.GetObj(i).mH-mPhaseRegistry.GetObj(i).mBiasH;
      REAL ty=mPhaseRegistry.GetObj(i).mK-mPhaseRegistry.GetObj(i).mBiasK;
      REAL tz=mPhaseRegistry.GetObj(i).mL-mPhaseRegistry.GetObj(i).mBiasL;
      mpData->GetCrystal().MillerToOrthonormalCoords(tx,ty,tz);
      
      cost +=(tx*tx+ty*ty+tz*tz)/mPhaseRegistry.GetObj(i).mNorm/.04;
   }
   VFN_DEBUG_MESSAGE("TextureMarchDollase::GetBiasingCost()="<<cost<<"("<<mName<<")",1)
   return cost;
}
void TextureMarchDollase::XMLOutput(ostream &os,int indent)const
{
   VFN_DEBUG_ENTRY("TextureMarchDollase::XMLOutput():",5)
   for(int i=0;i<indent;i++) os << "  " ;
   XMLCrystTag tag("TextureMarchDollase");
   os <<tag<<endl;
   indent++;
   
   for(int i=0;i<mPhaseRegistry.GetNb();i++) mPhaseRegistry.GetObj(i).XMLOutput(os,indent);
   
   indent--;
   tag.SetIsEndTag(true);
   for(int i=0;i<indent;i++) os << "  " ;
   os <<tag<<endl;
   VFN_DEBUG_EXIT("TextureMarchDollase::XMLOutput():"<<this->GetName(),5)
}

void TextureMarchDollase::XMLInput(istream &is,const XMLCrystTag &tagg)
{
   VFN_DEBUG_ENTRY("TextureMarchDollase::XMLInput():"<<this->GetName(),5)
   for(unsigned int i=0;i<tagg.GetNbAttribute();i++)
   {
      // No attribute to read
   }
   while(true)
   {
      XMLCrystTag tag(is);
      if(("TextureMarchDollase"==tag.GetName())&&tag.IsEndTag())
      {
         VFN_DEBUG_EXIT("TextureMarchDollase::XMLInput()",5)
         return;
      }
      if("TexturePhaseMarchDollase"==tag.GetName())
      {
         const long nb=mPhaseRegistry.GetNb();
         this->AddPhase(0.,1.,1.,0.,0.);
         mPhaseRegistry.GetObj(nb).XMLInput(is,tag);
      }
   }
}
void TextureMarchDollase::BeginOptimization(const bool allowApproximations,
                                            const bool enableRestraints)
{
   for(unsigned int i=0; i<this->GetNbPhase();i++)
   {
      RefinablePar *pH=&(this->GetPar(&(mPhaseRegistry.GetObj(i).mH)));
      RefinablePar *pK=&(this->GetPar(&(mPhaseRegistry.GetObj(i).mK)));
      RefinablePar *pL=&(this->GetPar(&(mPhaseRegistry.GetObj(i).mL)));
      if((pH->IsFixed()==false)||(pK->IsFixed()==false)||(pL->IsFixed()==false))
      {
         pH->SetIsFixed(false);
         pK->SetIsFixed(false);
         pL->SetIsFixed(false);
      }
      REAL tx=pH->GetValue();
      REAL ty=pK->GetValue();
      REAL tz=pL->GetValue();
      mpData->GetCrystal().MillerToOrthonormalCoords(tx,ty,tz);
      mPhaseRegistry.GetObj(i).mNorm=sqrt(tx*tx+ty*ty+tz*tz);
      // Something went wrong, preferred orientation vector has null norm !
      if(mPhaseRegistry.GetObj(i).mNorm<1e-6) mPhaseRegistry.GetObj(i).mNorm=1;
   }
   this->RefinableObj::BeginOptimization(allowApproximations,enableRestraints);
}
void TextureMarchDollase::TagNewBestConfig()const
{
   for(unsigned int i=0; i<this->GetNbPhase();i++)
   {
      mPhaseRegistry.GetObj(i).mBiasFraction  =mPhaseRegistry.GetObj(i).mFraction;
      mPhaseRegistry.GetObj(i).mBiasMarchCoeff=mPhaseRegistry.GetObj(i).mMarchCoeff;
      mPhaseRegistry.GetObj(i).mBiasH         =mPhaseRegistry.GetObj(i).mH;
      mPhaseRegistry.GetObj(i).mBiasK         =mPhaseRegistry.GetObj(i).mK;
      mPhaseRegistry.GetObj(i).mBiasL         =mPhaseRegistry.GetObj(i).mL;
   }
}

void TextureMarchDollase::CalcCorr() const
{
   const long nbReflUsed=mpData->GetNbReflBelowMaxSinThetaOvLambda();
   if(  (mClockTexturePar<mClockCorrCalc)
      &&(mpData->GetClockNbReflBelowMaxSinThetaOvLambda()<mClockCorrCalc)) return;
   VFN_DEBUG_ENTRY("TextureMarchDollase::CalcCorr()",3)
   TAU_PROFILE("TextureMarchDollase::CalcCorr()","void ()",TAU_DEFAULT);
   // normalizer for the sum of fractions, and non-texture fraction
   // (the sum of fractions must be equal to 1, but we cannot
   // modify fractions here since this is a const function)
      REAL fractionNorm=0,nonTexturedFraction;
      for(unsigned int i=0; i<this->GetNbPhase();i++) fractionNorm+=this->GetFraction(i);
      if(fractionNorm<1)
      {
         nonTexturedFraction= 1.-fractionNorm;
         fractionNorm=1.;
      }
      else nonTexturedFraction=0.;
   //compute correction for each phase
      const long nbRefl=mpData->GetNbRefl();
      mCorr.resize(nbRefl);
      mCorr=nonTexturedFraction;
      CrystVector_REAL reflNorm(nbRefl);
      {
         const REAL *xx=mpData->GetReflX().data();
         const REAL *yy=mpData->GetReflY().data();
         const REAL *zz=mpData->GetReflZ().data();
         for(long i=0;i<nbReflUsed;i++)
         {
            reflNorm(i)= sqrt(*xx * *xx + *yy * *yy + *zz * *zz);
            xx++;yy++;zz++;
         }
      }
      CrystMatrix_REAL hkl;
      for(unsigned int i=0; i<this->GetNbPhase();i++)
      {
         // We are using multiplicity for powder diffraction, therefore with only
         // unique reflections. But Equivalent reflections do not have the same
         // texture correction ! So we must use the symmetry oprators, and it is simpler
         // to apply the symmetries to the texture vector than to all reflections
         hkl=mpData->GetCrystal().GetSpaceGroup()
               .GetAllEquivRefl(this->GetPhaseH(i),this->GetPhaseK(i),this->GetPhaseL(i),true);
         //coefficients
            const REAL march=1./(this->GetMarchCoeff(i)+1e-6);
            const REAL march2=this->GetMarchCoeff(i)*this->GetMarchCoeff(i)-march;
            // Normalized by the number of symmetrical reflections
            const REAL frac=this->GetFraction(i)/(fractionNorm+1e-6)/hkl.rows();
         
         for(long j=0;j<hkl.rows();j++)
         {
            //orthonormal coordinates for T (texture) vector
               REAL tx=hkl(j,0),
                    ty=hkl(j,1),
                    tz=hkl(j,2);
            // reflection coordinates
               const REAL *xx=mpData->GetReflX().data();
               const REAL *yy=mpData->GetReflY().data();
               const REAL *zz=mpData->GetReflZ().data();
               const REAL *xyznorm=reflNorm.data();
               {
                  mpData->GetCrystal().MillerToOrthonormalCoords(tx,ty,tz);
                  const REAL norm=sqrt(tx*tx+ty*ty+tz*tz);
                  tx/=(norm+1e-6);
                  ty/=(norm+1e-6);
                  tz/=(norm+1e-6);
               }
            // Calculation
               REAL tmp;
               for(long k=0;k<nbReflUsed;k++)
               {
                  tmp=(tx * (*xx++) + ty * (*yy++) + tz * (*zz++))/ (*xyznorm++);
                  tmp=march+march2*tmp*tmp;
                  if(tmp<0) tmp=0;// rounding errors ?
                  mCorr(k)+=frac*pow((float)tmp,(float)-1.5);
               }
         }
      }
   //if(mIsbeingRefined==false)
   //{
   //   cout <<FormatVertVectorHKLFloats<REAL>(mpData->GetH(),
   //                                          mpData->GetK(),
   //                                          mpData->GetL(),
   //                                          mpData->GetSinThetaOverLambda(),
   //                                          mCorr)<<endl;
   //   this->Print();
   //}
   mClockCorrCalc.Click();
   VFN_DEBUG_EXIT("TextureMarchDollase::CalcCorr()",3)
}
void TextureMarchDollase::DeleteAllPhase()
{
}

#ifdef __WX__CRYST__
WXCrystObjBasic* TextureMarchDollase::WXCreate(wxWindow* parent)
{
   //:TODO: Check mpWXCrystObj==0
   //mpWXCrystObj=new WXTextureMarchDollase(parent,this);
   return mpWXCrystObj;
}
#endif
////////////////////////////////////////////////////////////////////////
//
//        Time Of Flight correction
//
////////////////////////////////////////////////////////////////////////
TOFCorr::TOFCorr(const ScatteringData & data):
ScatteringCorr(data)
{}

TOFCorr::~TOFCorr()
{}

const string & TOFCorr::GetName() const
{
   //So far, we do not need a personalized name...
   const static string mName="TOFCorr";
   return mName;
}

const string & TOFCorr::GetClassName() const
{
   const static string className="TOFCorr";
   return className;
}

void TOFCorr::CalcCorr() const
{
   const REAL *pstol=mpData->GetSinThetaOverLambda().data();
   if(mpData->GetClockTheta()<mClockCorrCalc) return;
   TAU_PROFILE("TOFCorr::CalcCorr()","void ()",TAU_DEFAULT);
   mCorr.resize(mpData->GetNbRefl());
   for(long i=0;i<mpData->GetNbRefl();i++) mCorr(i) = pow((float)(1.0/(2.0* *pstol++)),(int)4);
   mClockCorrCalc.Click();
}

}//namespace
