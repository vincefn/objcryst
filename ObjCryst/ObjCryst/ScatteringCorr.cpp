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
       && ((mPolarAfactor-((1-f)/(1+f))) <(mPolarAfactor*.0001)) ) return;
   VFN_DEBUG_MESSAGE("PolarizationCorr::CalcCorr()",10)
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
void TexturePhaseMarchDollase::XMLOutput(ostream &os,int indent=0)const
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
ScatteringCorr(data)
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
      tmp.SetGlobalOptimStep(.01);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("MarchCoeff_"+(string)buf,&(phase->mMarchCoeff),.01,100.,
                        gpRefParTypeScattDataCorrInt,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,1.);
      tmp.AssignClock(mClockTexturePar);
      tmp.SetDerivStep(1e-7);
      tmp.SetGlobalOptimStep(.01);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("H_"+(string)buf,&(phase->mH),-10.,10.,
                        gpRefParTypeScattDataCorrInt,
                        REFPAR_DERIV_STEP_ABSOLUTE,false,true,true,false,1.);
      tmp.AssignClock(mClockTexturePar);
      tmp.SetDerivStep(1e-7);
      tmp.SetGlobalOptimStep(.1);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("K_"+(string)buf,&(phase->mK),-10.,10.,
                        gpRefParTypeScattDataCorrInt,
                        REFPAR_DERIV_STEP_ABSOLUTE,false,true,true,false,1.);
      tmp.AssignClock(mClockTexturePar);
      tmp.SetDerivStep(1e-7);
      tmp.SetGlobalOptimStep(.1);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("L_"+(string)buf,&(phase->mL),-10.,10.,
                        gpRefParTypeScattDataCorrInt,
                        REFPAR_DERIV_STEP_ABSOLUTE,false,true,true,false,1.);
      tmp.AssignClock(mClockTexturePar);
      tmp.SetDerivStep(1e-7);
      tmp.SetGlobalOptimStep(.1);
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

void TextureMarchDollase::GlobalOptRandomMove(const REAL mutationAmplitude)
{
   if(mRandomMoveIsDone) return;
   if((rand()/(REAL)RAND_MAX)<.3)//only 30% proba to make a random move
   {
      VFN_DEBUG_MESSAGE("TextureMarchDollase::GlobalOptRandomMove()",1)
      for(unsigned int i=0;i<this->GetNbPhase();i++)
      {
         RefinablePar *pF=&(this->GetPar(&(mPhaseRegistry.GetObj(i).mFraction)));
         RefinablePar *pM=&(this->GetPar(&(mPhaseRegistry.GetObj(i).mMarchCoeff)));
         RefinablePar *pH=&(this->GetPar(&(mPhaseRegistry.GetObj(i).mH)));
         RefinablePar *pK=&(this->GetPar(&(mPhaseRegistry.GetObj(i).mK)));
         RefinablePar *pL=&(this->GetPar(&(mPhaseRegistry.GetObj(i).mL)));
         if(pF->IsFixed()==false)
            pF->Mutate(pF->GetGlobalOptimStep()*2*(rand()/(REAL)RAND_MAX-0.5)*mutationAmplitude);
         if((pM->IsFixed()==false) &&
            ((pH->IsFixed()==false)||(pK->IsFixed()==false)||(pL->IsFixed()==false)))
         {
            REAL tx=pH->GetValue();
            REAL ty=pK->GetValue();
            REAL tz=pL->GetValue();
            mpData->GetCrystal().MillerToOrthonormalCoords(tx,ty,tz);

            tx += .01*2*(rand()/(REAL)RAND_MAX-0.5)*mutationAmplitude*mPhaseRegistry.GetObj(i).mNorm;
            ty += .01*2*(rand()/(REAL)RAND_MAX-0.5)*mutationAmplitude*mPhaseRegistry.GetObj(i).mNorm;
            tz += .01*2*(rand()/(REAL)RAND_MAX-0.5)*mutationAmplitude*mPhaseRegistry.GetObj(i).mNorm;

            const REAL factor=mPhaseRegistry.GetObj(i).mNorm/sqrt(tx*tx+ty*ty+tz*tz);
            pM->MutateTo(pM->GetValue()/factor);
            tx *= factor;
            ty *= factor;
            tz *= factor;
            mpData->GetCrystal().OrthonormalToMillerCoords(tx,ty,tz);
            pH->MutateTo(tx);
            pK->MutateTo(ty);
            pL->MutateTo(tz);
         }
         else
         {
            if(pM->IsFixed()==false)
            {
               pM->MutateTo(pM->GetValue()*(1.+.012*(rand()/(REAL)RAND_MAX-0.5)*mutationAmplitude));
            }
            if((pH->IsFixed()==false)||(pK->IsFixed()==false)||(pL->IsFixed()==false))
            {
               REAL tx=pH->GetValue();
               REAL ty=pK->GetValue();
               REAL tz=pL->GetValue();
               mpData->GetCrystal().MillerToOrthonormalCoords(tx,ty,tz);
               REAL norm=sqrt(tx*tx+ty*ty+tz*tz);

               tx += .01*2*(rand()/(REAL)RAND_MAX-0.5)*mutationAmplitude*norm;
               ty += .01*2*(rand()/(REAL)RAND_MAX-0.5)*mutationAmplitude*norm;
               tz += .01*2*(rand()/(REAL)RAND_MAX-0.5)*mutationAmplitude*norm;

               REAL factor=norm/sqrt(tx*tx+ty*ty+tz*tz);
               tx *= factor;
               ty *= factor;
               tz *= factor;
               mpData->GetCrystal().OrthonormalToMillerCoords(tx,ty,tz);
               pH->MutateTo(tx);
               pK->MutateTo(ty);
               pL->MutateTo(tz);
            }
         }
      }
   }
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
      
      cost +=(tx*tx+ty*ty+tz*tz)/.04;
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
   if(mClockTexturePar<mClockCorrCalc) return;
   VFN_DEBUG_ENTRY("TextureMarchDollase::CalcCorr()",3)
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
         for(long i=0;i<nbRefl;i++)
         {
            reflNorm(i)= sqrt(*xx * *xx + *yy * *yy + *zz * *zz);
            xx++;yy++;zz++;
         }
      }
      for(unsigned int i=0; i<this->GetNbPhase();i++)
      {
         //coefficients
            const REAL march=1./this->GetMarchCoeff(i);
            const REAL march2=this->GetMarchCoeff(i)*this->GetMarchCoeff(i)-march;
            const REAL frac=this->GetFraction(i)/fractionNorm;
         // reflection coordinates
            const REAL *xx=mpData->GetReflX().data();
            const REAL *yy=mpData->GetReflY().data();
            const REAL *zz=mpData->GetReflZ().data();
            const REAL *xyznorm=reflNorm.data();
         
         //orthonormal coordinates for T (texture) vector
            REAL tx=this->GetPhaseH(i),
                 ty=this->GetPhaseK(i),
                 tz=this->GetPhaseL(i);
            {
               mpData->GetCrystal().MillerToOrthonormalCoords(tx,ty,tz);
               const REAL norm=sqrt(tx*tx+ty*ty+tz*tz);
               tx/=norm;
               ty/=norm;
               tz/=norm;
            }
         // Calculation
            REAL tmp;
            for(long i=0;i<nbRefl;i++)
            {
               tmp=(tx * (*xx++) + ty * (*yy++) + tz * (*zz++))/ (*xyznorm++);
               mCorr(i)+=frac*powf(march+march2*tmp*tmp,-1.5);
            }
            
      }
   if(mIsbeingRefined==false)
   {
      cout <<FormatVertVectorHKLFloats<REAL>(mpData->GetH(),
                                             mpData->GetK(),
                                             mpData->GetL(),
                                             mpData->GetSinThetaOverLambda(),
                                             mCorr)<<endl;
      this->Print();
   }
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

}//namespace
