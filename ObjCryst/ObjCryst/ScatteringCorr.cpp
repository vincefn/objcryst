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
#include <stdio.h> //for sprintf()

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
       && (mPolarAfactor==((1-f)/(1+f)))) return;
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
//        TextureMarchDollase
//
////////////////////////////////////////////////////////////////////////
TextureMarchDollase::TextureMarchDollase(const ScatteringData & data):
ScatteringCorr(data),mNbPhase(0),mpTexturePhase(0)
{}

TextureMarchDollase::~TextureMarchDollase()
{
   for(unsigned int i=0;i<mNbPhase;i++) delete mpTexturePhase[i];
   delete[] mpTexturePhase;
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
   TexturePhase** tmp=new TexturePhase*[mNbPhase+1];
   if(mNbPhase>0)
   {
      for(unsigned int i=0;i<mNbPhase;i++) tmp[i]=mpTexturePhase[i];
      delete[] mpTexturePhase;
      mpTexturePhase=tmp;
   }
   mpTexturePhase[mNbPhase]=new TexturePhase(f,c,h,k,l);
   
   //Add parameters
   char buf [3];
   sprintf(buf,"%d",mNbPhase);
   {
      RefinablePar tmp("Fraction_"+(string)buf,&(mpTexturePhase[mNbPhase]->mFraction),0,1.,
                        gpRefParTypeScattDataCorrInt,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,1.);
      tmp.AssignClock(mClockTexturePar);
      tmp.SetDerivStep(1e-7);
      tmp.SetGlobalOptimStep(.1);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("MarchCoeff_"+(string)buf,&(mpTexturePhase[mNbPhase]->mFraction),0,1.,
                        gpRefParTypeScattDataCorrInt,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,1.);
      tmp.AssignClock(mClockTexturePar);
      tmp.SetDerivStep(1e-7);
      tmp.SetGlobalOptimStep(.1);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("H_"+(string)buf,&(mpTexturePhase[mNbPhase]->mFraction),0,1.,
                        gpRefParTypeScattDataCorrInt,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,1.);
      tmp.AssignClock(mClockTexturePar);
      tmp.SetDerivStep(1e-7);
      tmp.SetGlobalOptimStep(.1);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("K_"+(string)buf,&(mpTexturePhase[mNbPhase]->mFraction),0,1.,
                        gpRefParTypeScattDataCorrInt,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,1.);
      tmp.AssignClock(mClockTexturePar);
      tmp.SetDerivStep(1e-7);
      tmp.SetGlobalOptimStep(.1);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("L_"+(string)buf,&(mpTexturePhase[mNbPhase]->mFraction),0,1.,
                        gpRefParTypeScattDataCorrInt,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,1.);
      tmp.AssignClock(mClockTexturePar);
      tmp.SetDerivStep(1e-7);
      tmp.SetGlobalOptimStep(.1);
      this->AddPar(tmp);
   }
   mNbPhase++;
}

void TextureMarchDollase::SetPhasePar(const unsigned int i, const REAL f, const REAL c,
                 const REAL h,const REAL k, const REAL l)
{
   mpTexturePhase[i]->SetPar(f,c,h,k,l);
}

unsigned int TextureMarchDollase::GetNbPhase() const {return mNbPhase;}

REAL TextureMarchDollase::GetFraction(const unsigned int i)const
{return mpTexturePhase[i]->mFraction;}

REAL TextureMarchDollase::GetMarchCoeff(const unsigned int i)const
{return mpTexturePhase[i]->mMarchCoeff;}

REAL TextureMarchDollase::GetPhaseH(const unsigned int i)const
{return mpTexturePhase[i]->mH;}

REAL TextureMarchDollase::GetPhaseK(const unsigned int i)const
{return mpTexturePhase[i]->mK;}

REAL TextureMarchDollase::GetPhaseL(const unsigned int i)const
{return mpTexturePhase[i]->mL;}

void TextureMarchDollase::GlobalOptRandomMove(const REAL mutationAmplitude)
{
}
REAL TextureMarchDollase::GetBiasingCost()const
{
   REAL cost=0;
   //for(int i=0;i<mNbRestraint;i++) cost+=mpRestraint[i]->GetRestraintCost();
   for(int i=0;i<this->GetNbPar();i++) cost+=this->GetPar(i).GetBiasingCost();
   VFN_DEBUG_MESSAGE("TextureMarchDollase::GetBiasingCost()="<<cost<<"("<<mName<<")",1)
   return cost;
}
void TextureMarchDollase::XMLOutput(ostream &os,int indent)const
{
}
void TextureMarchDollase::XMLInput(istream &is,const XMLCrystTag &tag)
{
}
void TextureMarchDollase::CalcCorr() const
{
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
