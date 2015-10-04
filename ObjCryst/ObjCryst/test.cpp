/*  ObjCryst++ Object-Oriented Crystallographic Library
    (c) 2000-2002 Vincent Favre-Nicolin vincefn@users.sourceforge.net
        2000-2001 University of Geneva (Switzerland)

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation;  version 2 of the License.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/*   test.cpp
*  source file for test functions (speed, etc...)
*
*/
#include <stdlib.h>
#include <list>
#include "ObjCryst/ObjCryst/test.h"
#include "ObjCryst/ObjCryst/Crystal.h"
#include "ObjCryst/ObjCryst/Atom.h"
#include "ObjCryst/ObjCryst/DiffractionDataSingleCrystal.h"
#include "ObjCryst/ObjCryst/PowderPattern.h"
#include "ObjCryst/RefinableObj/GlobalOptimObj.h"
#include "ObjCryst/Quirks/VFNStreamFormat.h"

namespace ObjCryst
{

SpeedTestReport SpeedTest(const unsigned int nbAtom, const int nbAtomType,const string spacegroup,
                          const RadiationType radiation, const unsigned long nbReflections,
                          const unsigned int dataType,const REAL time)
{
   Crystal cryst(9,11,15,1.2,1.3,1.7,spacegroup);
   for(int i=0;i<nbAtomType;++i)
   {
      cryst.AddScatteringPower(new ScatteringPowerAtom("O","O",1.5));
   }
   for(unsigned int i = 0; i < nbAtom; ++i)
   {
      cryst.AddScatterer(new Atom(.0,.0,.0,"O",
                                  &(cryst.GetScatteringPowerRegistry().GetObj(i%nbAtomType)),
                                  1.));
   }
   cryst.SetUseDynPopCorr(false);

   RefinableObj *pData = NULL;
   switch(dataType)
   {
      case 0:
      {
         DiffractionDataSingleCrystal *pDataTmp=new DiffractionDataSingleCrystal;
         pDataTmp->SetWavelength(0.25);
         pDataTmp->SetRadiationType(radiation);
         pDataTmp->SetMaxSinThetaOvLambda(100.);
         pDataTmp->SetCrystal(cryst);
         float maxtheta=0.1;
         for(;;)
         {
            pDataTmp->GenHKLFullSpace(maxtheta, true);
            if(pDataTmp->GetNbRefl()>(long)nbReflections) break;
            maxtheta*=1.5;
            if(maxtheta>=M_PI/2.) break;
         }
         CrystVector_REAL hh; hh=pDataTmp->GetH();hh.resizeAndPreserve(nbReflections);hh+=0.0001;
         CrystVector_REAL kk; kk=pDataTmp->GetK();kk.resizeAndPreserve(nbReflections);kk+=0.0001;
         CrystVector_REAL ll; ll=pDataTmp->GetL();ll.resizeAndPreserve(nbReflections);ll+=0.0001;

         CrystVector_long h(nbReflections); h=hh;
         CrystVector_long k(nbReflections); k=kk;
         CrystVector_long l(nbReflections); l=ll;

         CrystVector_REAL iobs(nbReflections);
         for(unsigned int i=0;i<nbReflections;++i) iobs(i)=(REAL)rand();
         CrystVector_REAL sigma(nbReflections);sigma=1;

         pDataTmp->SetHklIobs (h, k, l, iobs, sigma);
         pDataTmp->SetWeightToInvSigma2();

         pData=pDataTmp;
         break;
      }
      case 1:
      {
         PowderPattern *pDataTmp=new PowderPattern;
         pDataTmp->SetWavelength(0.25);
         pDataTmp->SetRadiationType(radiation);
         pDataTmp->SetPowderPatternPar(0.001,.001,3140);
         CrystVector_REAL iobs(3140);
         for(unsigned int i=0;i<3140;++i) iobs(i)=(REAL)rand()+1.;
         pDataTmp->SetPowderPatternObs(iobs);
         pDataTmp->SetMaxSinThetaOvLambda(100.);

         PowderPatternBackground *backgdData= new PowderPatternBackground;
         backgdData->SetName("PbSo4-background");
         {
            CrystVector_REAL tth(2),backgd(2);
            tth(0)=0.;tth(1)=3.14;
            backgd(0)=1.;backgd(1)=9.;
            backgdData->SetInterpPoints(tth,backgd);
         }
         pDataTmp->AddPowderPatternComponent(*backgdData);

         PowderPatternDiffraction * diffData=new PowderPatternDiffraction;
         diffData->SetCrystal(cryst);
         pDataTmp->AddPowderPatternComponent(*diffData);
         diffData->SetName("Crystal phase");
         diffData->SetReflectionProfilePar(PROFILE_PSEUDO_VOIGT,
                                           .03*DEG2RAD*DEG2RAD,0.,0.,0.3,0);
         {
            float maxtheta=0.1;
            for(;;)
            {
               diffData->ScatteringData::GenHKLFullSpace(maxtheta, true);
               if(diffData->GetNbRefl()>(long)nbReflections) break;
               maxtheta*=1.5;
               if(maxtheta>=M_PI/2.) break;
            }
            CrystVector_REAL hh; hh=diffData->GetH();hh.resizeAndPreserve(nbReflections);
            CrystVector_REAL kk; kk=diffData->GetK();kk.resizeAndPreserve(nbReflections);
            CrystVector_REAL ll; ll=diffData->GetL();ll.resizeAndPreserve(nbReflections);

            diffData->SetHKL (hh, kk, ll);
         }
         pData=pDataTmp;
         break;
      }
   }

   //Create the global optimization object
      MonteCarloObj *pGlobalOptObj=new MonteCarloObj;
      pGlobalOptObj->AddRefinableObj(*pData);
      pGlobalOptObj->AddRefinableObj(cryst);

   //Refine only positionnal parameters
      pGlobalOptObj->FixAllPar();
      pGlobalOptObj->SetParIsFixed(gpRefParTypeScattTransl,false);
      pGlobalOptObj->SetParIsFixed(gpRefParTypeScattOrient,false);

   //Don't cheat ;-)
      pGlobalOptObj->RandomizeStartingConfig();

   //Annealing parameters (schedule, Tmax, Tmin, displacement schedule,
      pGlobalOptObj->SetAlgorithmParallTempering(ANNEALING_SMART,1e8,1e-8,
                                               ANNEALING_EXPONENTIAL,8,.125);

   //Global Optimization
      //The real job-first test
      long nbTrial=50000000;
      pGlobalOptObj->Optimize(nbTrial,true,0,time);


   SpeedTestReport report;
   report.mNbAtom=nbAtom;
   report.mNbAtomType=nbAtomType;
   report.mSpacegroup=spacegroup;
   report.mRadiation=radiation;
   report.mNbReflections=nbReflections;
   report.mDataType=dataType;
   report.mBogoMRAPS=(REAL)nbAtom*cryst.GetSpaceGroup().GetNbSymmetrics()*(REAL)nbReflections
                     *(50000000-nbTrial)/pGlobalOptObj->GetLastOptimElapsedTime()/1e6;
   report.mBogoMRAPS_reduced=(REAL)nbAtom*cryst.GetSpaceGroup().GetNbSymmetrics(true,true)
                             *(REAL)nbReflections
                             *(50000000-nbTrial)/pGlobalOptObj->GetLastOptimElapsedTime()/1e6;
   report.mBogoSPS=(50000000-nbTrial)/pGlobalOptObj->GetLastOptimElapsedTime();
   delete pGlobalOptObj;
   delete pData;
   return report;
}
}
