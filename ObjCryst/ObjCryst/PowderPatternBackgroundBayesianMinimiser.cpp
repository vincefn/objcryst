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
/*   PowderPatternBackgroundBayesianMinimiser.h
*  Bayesian estimation of powder pattern background
*
*/

#include "ObjCryst/ObjCryst/PowderPatternBackgroundBayesianMinimiser.h"
namespace ObjCryst
{
PowderPatternBackgroundBayesianMinimiser::
   PowderPatternBackgroundBayesianMinimiser(PowderPatternBackground &backgd):
mpBackground(&backgd)
{
   this->AddSubRefObj(*mpBackground);
}

PowderPatternBackgroundBayesianMinimiser::~PowderPatternBackgroundBayesianMinimiser()
{
   this->RemoveSubRefObj(*mpBackground);
}

const string& PowderPatternBackgroundBayesianMinimiser::GetClassName()const
{
   static string className="PowderPatternBackgroundBayesianMinimiser";
   return className;
}

REAL PowderPatternBackgroundBayesianMinimiser::GetLogLikelihood()const
{
   TAU_PROFILE("PowderPatternBackgroundBayesianMinimiser::GetLogLikelihood()","void ()",TAU_DEFAULT);
   REAL llk=0;
   const long nb=mpBackground->GetPowderPatternCalc().numElements();
   const long step=1;
   {// Calc (obs-calc)/sigma
      const REAL *pBackgd=mpBackground->GetPowderPatternCalc().data();
      const REAL *pObs=mpBackground->GetParentPowderPattern().GetPowderPatternObs().data();
      const REAL *pSigma=mpBackground->GetParentPowderPattern().GetPowderPatternObsSigma().data();
      for(long i=0;i<nb;i+=step)
      {
         if(*pSigma>0)
         {
            llk += PowderPatternBackgroundBayesianMinimiser
                     ::BayesianBackgroundLogLikelihood
                        ((*pObs-*pBackgd) / (1.4142135623730951**pSigma));
         }
         pObs+=step;
         pBackgd+=step;
         pSigma+=step;
      }
   }
   return llk;
}

unsigned int PowderPatternBackgroundBayesianMinimiser::GetNbLSQFunction () const {return 1;}

const CrystVector_REAL& PowderPatternBackgroundBayesianMinimiser::GetLSQCalc (const unsigned int id) const
{
  const long nb=mpBackground->GetPowderPatternCalc().numElements();
  const long step=1;
  if(mBayesianCalc.numElements()!=nb) mBayesianCalc.resize(nb);

  const REAL *pBackgd=mpBackground->GetPowderPatternCalc().data();
  const REAL *pObs=mpBackground->GetParentPowderPattern().GetPowderPatternObs().data();
  const REAL *pSigma=mpBackground->GetParentPowderPattern().GetPowderPatternObsSigma().data();
  REAL *pBayesCalc=mBayesianCalc.data();
  for(long i=0;i<nb;i+=step)
  {
    if(*pSigma>0)
    {
      *pBayesCalc = 1 + PowderPatternBackgroundBayesianMinimiser::BayesianBackgroundLogLikelihood((*pObs-*pBackgd) / (1.4142135623730951**pSigma));
    }
    else *pBayesCalc = 1;
    pObs+=step;
    pBackgd+=step;
    pSigma+=step;
    pBayesCalc+=step;
  }
  return mBayesianCalc;
}

const CrystVector_REAL& PowderPatternBackgroundBayesianMinimiser::GetLSQObs (const unsigned int) const
{
  const long nb=mpBackground->GetPowderPatternCalc().numElements();
  if(mBayesianObs.numElements()!=nb)
  {
    mBayesianObs.resize(nb);
    mBayesianObs=1; // Avoid having all observed values==0, raises issues when computing R and Rw in LSQNumObj
  }
  return mBayesianObs;
}

const CrystVector_REAL& PowderPatternBackgroundBayesianMinimiser::GetLSQWeight (const unsigned int) const
{
  const long nb=mpBackground->GetPowderPatternCalc().numElements();
  if(mBayesianWeight.numElements()!=nb)
  {
     mBayesianWeight.resize(nb);
     mBayesianWeight=1;
  }
  return mBayesianWeight;
}

REAL PowderPatternBackgroundBayesianMinimiser::BayesianBackgroundLogLikelihood(const REAL t)
{
   static const REAL vllk[11]={0.00000000e+00,
                                1e-4,
                                1.77123249e-03,
                                1.00997634e+00,
                                2.89760310e+00,
                                3.61881096e+00,
                                3.93024374e+00,
                                4.16063018e+00,
                                4.34600620e+00,
                                4.50155649e+00,
                                4.63573160e+00};
   static const REAL vt[11]={  0. ,   0.01,  0.1 ,  1.1 ,  2.1 ,  3.1 ,  4.1 ,  5.1 ,  6.1 ,  7.1 ,  8.1};
   static const CubicSpline spline(vt,vllk,11);
   static const REAL s1=spline(8)-log((REAL)8);
   if(t<=0) return 5*t*t;
   if(t<8)return spline(t);
   return s1+log(t);
}


}//namespace
