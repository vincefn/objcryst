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
//        
//
////////////////////////////////////////////////////////////////////////


}
