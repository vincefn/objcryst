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

#include "ObjCryst/PowderPatternBackgroundBayesianMinimiser.h"
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
   const long step=nb/500+1;
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

//integral{epsilon->+inf}[1/u*exp(-(t-u)^2)du]

//    t\epsilon    1.0e-06      1.0e-05      1.0e-04      0.0010       0.01000    0.100000
//    0.00000      2.60443      2.41787      2.18830      1.88979      1.46229     0.702584
//   0.200000      2.59357      2.41288      2.19215      1.90839      1.51004     0.826473
//   0.400000      2.51008      2.33647      2.12621      1.85953      1.49353     0.889141
//   0.600000      2.35717      2.19223      1.99459      1.74786      1.41794     0.895476
//   0.800000      2.13979      1.98568      1.80339      1.58012      1.29034     0.851766
//    1.00000      1.86570      1.72521      1.56168      1.36591      1.12028     0.765916
//    1.20000      1.54687      1.42345      1.28260      1.11844     0.920072     0.647484
//    1.40000      1.20106      1.09842     0.984005     0.854640     0.704420     0.507339
//    1.60000     0.852220     0.773306     0.687608     0.593751     0.488919     0.356793
//    1.80000     0.527629     0.472979     0.415154     0.353694     0.287337     0.206228
//    2.00000     0.249382     0.216138     0.181717     0.145989     0.108731    0.0633822
//    2.20000     0.0271344    0.00959706  -0.00828458   -0.0265562   -0.0438729   -0.0668048
//    2.40000    -0.161050    -0.161093    -0.161512    -0.164514    -0.172413    -0.182920
//    2.60000    -0.276636    -0.276654    -0.276827    -0.278066    -0.281322    -0.285711
//    2.80000    -0.373540    -0.373547    -0.373612    -0.374076    -0.375300    -0.377011
//    3.00000    -0.457712    -0.457719    -0.457779    -0.458069    -0.458523    -0.458850
//    3.20000    -0.533130    -0.533132    -0.533151    -0.533242    -0.533383    -0.533564
//    3.40000    -0.601793    -0.601794    -0.601799    -0.601825    -0.601866    -0.601921
//    3.60000    -0.665178    -0.665178    -0.665180    -0.665187    -0.665199    -0.665227
//    3.80000    -0.724197    -0.724197    -0.724198    -0.724200    -0.724208    -0.724255
//    4.00000    -0.785116    -0.785116    -0.785113    -0.785086    -0.784814    -0.779562
//    4.20000    -0.832013    -0.832013    -0.832013    -0.832013    -0.832013    -0.832011
//    4.40000    -0.881490    -0.881490    -0.881490    -0.881490    -0.881490    -0.881490
//    4.60000    -0.928466    -0.928466    -0.928466    -0.928466    -0.928467    -0.928475
//    4.80000    -0.973167    -0.973167    -0.973167    -0.973168    -0.973172    -0.973201
//    5.00000     -1.01576     -1.01576     -1.01576     -1.01576     -1.01577     -1.01585
//    5.20000     -1.05686     -1.05686     -1.05686     -1.05686     -1.05686     -1.06303
//    5.40000     -1.09609     -1.09609     -1.09609     -1.09609     -1.09609     -1.09609
//    5.60000     -1.13377     -1.13377     -1.13377     -1.13377     -1.13377     -1.13378
//    5.80000     -1.17001     -1.17001     -1.17001     -1.17001     -1.17001     -1.17003
//    6.00000     -1.20487     -1.20487     -1.20487     -1.20487     -1.20487     -1.20493
//    6.20000     -1.23873     -1.23873     -1.23873     -1.23873     -1.24505     -1.23851
//    6.40000     -1.27134     -1.27134     -1.27134     -1.27134     -1.27134     -1.27134
//    6.60000     -1.30289     -1.30289     -1.30289     -1.30289     -1.30289     -1.30289
//    6.80000     -1.33343     -1.33343     -1.33343     -1.33343     -1.33343     -1.33344
//    7.00000     -1.36300     -1.36300     -1.36300     -1.36300     -1.36301     -1.36304
//    7.20000     -1.39160     -1.39160     -1.39160     -1.39160     -1.39161     -1.39170
//    7.40000     -1.41977     -1.41977     -1.41977     -1.41977     -1.41977     -1.41977
//    7.60000     -1.44694     -1.44694     -1.44694     -1.44694     -1.44694     -1.44694
//    7.80000     -1.47337     -1.47337     -1.47337     -1.47337     -1.47337     -1.47337
//    8.00000     -1.49908     -1.49908     -1.49908     -1.49908     -1.49908     -1.49910
REAL PowderPatternBackgroundBayesianMinimiser::BayesianBackgroundLogLikelihood(const REAL t)
{
   // We could probably reduce the number of points to get a speed boost
   static const REAL vllk[41]={2.60443,
                               2.59357,
                               2.51008,
                               2.35717,
                               2.13979,
                               1.86570,
                               1.54687,
                               1.20106,
                               0.852220,
                               0.527629,
                               0.249382,
                               0.0271344,
                              -0.161050, 
                              -0.276636, 
                              -0.373540, 
                              -0.457712, 
                              -0.533130, 
                              -0.601793, 
                              -0.665178, 
                              -0.724197, 
                              -0.785116, 
                              -0.832013, 
                              -0.881490, 
                              -0.928466, 
                              -0.973167, 
                              -1.01576,  
                              -1.05686,  
                              -1.09609,  
                              -1.13377,  
                              -1.17001,  
                              -1.20487,  
                              -1.23873,  
                              -1.27134,  
                              -1.30289,  
                              -1.33343,  
                              -1.36300,  
                              -1.39160,  
                              -1.41977,  
                              -1.44694,  
                              -1.47337,  
                              -1.49908}; 
   static const REAL vt[41]={  0.00000,  
                               0.20000,  
                               0.40000,  
                               0.60000,  
                               0.80000,  
                               1.00000,  
                               1.20000,  
                               1.40000,  
                               1.60000,  
                               1.80000,  
                               2.00000,  
                               2.20000,  
                               2.40000,  
                               2.60000,  
                               2.80000,  
                               3.00000,  
                               3.20000,  
                               3.40000,  
                               3.60000,  
                               3.80000,  
                               4.00000,  
                               4.20000,  
                               4.40000,  
                               4.60000,  
                               4.80000,  
                               5.00000,  
                               5.20000,  
                               5.40000,  
                               5.60000,  
                               5.80000,  
                               6.00000,  
                               6.20000,  
                               6.40000,  
                               6.60000,  
                               6.80000,  
                               7.00000,  
                               7.20000,  
                               7.40000,  
                               7.60000,  
                               7.80000,  
                               8.00000};
   static const CubicSpline spline(vt,vllk,41);
   static const REAL s0=spline(0);
   static const REAL s1=s0-spline(8)-log(8);
   if(t<=0) return t*t;
   if(t<8)return s0-spline(t);
   return s1+log(t);
}


}//namespace
