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
/*   Simplex.cpp
*  source file for Conjugate Gradient Algorithm object
*
*/
#include "RefinableObj/Simplex.h"
#include "Quirks/VFNStreamFormat.h"

namespace ObjCryst
{
SimplexObj::SimplexObj(const string name):
OptimizationObj(name)
{
}
void SimplexObj::Optimize(long &nbSteps,const bool silent,const REAL finalcost,
                          const REAL maxTime)
{
   VFN_DEBUG_ENTRY("SimplexObj::Optimize()",10)
   for(int i=0;i<mRefinedObjList.GetNb();i++) mRefinedObjList.GetObj(i).BeginOptimization(true);
   this->PrepareRefParList();
   const unsigned long n=mRefParList.GetNbParNotFixed();
   // Create the n+1 set of parameters (n+1 vertices for the initial simplex)
   // To obtain the n new vertices we just move from the starting point
   // by the amount of the declared global optimization step.
      CrystVector_long vIndex(n+1);
      vIndex(0)=mRefParList.CreateParamSet();
      CrystVector_REAL vLLK(n+1);
      vLLK(0)=this->GetLogLikelihood();
      for(unsigned long i=0;i<n;i++)
      {
         vIndex(0)=mRefParList.CreateParamSet();
         mRefParList.GetParNotFixed(i).
            Mutate(mRefParList.GetParNotFixed(i).GetGlobalOptimStep()*100.0);
         vIndex(i+1)=mRefParList.CreateParamSet();
         vLLK(i+1)=this->GetLogLikelihood();
         mRefParList.RestoreParamSet(vIndex(0));
      }
   unsigned long best=0,worst=0,nextworst=0;
   for(;nbSteps>=0;--nbSteps)
   {
      // determine best, worst and next worst points
         if(vLLK(0)>vLLK(1)){worst=0;nextworst=1;}
         else{worst=1;nextworst=0;}
         for(unsigned long i=1;i<=n;i++)
         {
            if(vLLK(i)<=vLLK(best)) best =i;
            if(vLLK(i)>vLLK(worst))
            {
               nextworst=worst;
               worst=i;
            }
            else if((vLLK(i)>vLLK(nextworst)) && (i!=worst)) nextworst=i;
         }
      {
         CrystVector_REAL center(n);
         center = 0.0;
         for(unsigned long i=0;i<=n;i++)
         {
            if(i==worst) continue;
            center += mRefParList.GetParamSet(vIndex(i));
         }
         center /= (REAL)n;
         CrystVector_REAL worstdiff;
         worstdiff=mRefParList.GetParamSet(vIndex(worst));
         worstdiff-=center;
         for(unsigned long i=0;i<n;i++) worstdiff(i)/=mRefParList.GetParNotFixed(i).GetGlobalOptimStep();
         
         if(!silent) cout<<"Simplex:cycle="<<nbSteps<<", cost="<<vLLK(best)
                         <<",best="<<best<<",worst="<<worst<<",nextworst="<<nextworst
                         <<", Worst diff="<<abs(worstdiff.min())+abs(worstdiff.max())<<endl;
         if((abs(worstdiff.min())+abs(worstdiff.max()))<0.1) break;
      }
      #if 0
      cout<<FormatHorizVector<REAL>(vLLK)<<endl;
      cout<<FormatVertVector<REAL>(mRefParList.GetParamSet(vIndex(best)),
                                   mRefParList.GetParamSet(vIndex(worst)),
                                   mRefParList.GetParamSet(vIndex(nextworst)),
                                   mRefParList.GetParamSet(vIndex(6)),
                                   mRefParList.GetParamSet(vIndex(8)),
                                   mRefParList.GetParamSet(vIndex(10)))<<endl;
      #endif
      mRefParList.RestoreParamSet(vIndex(best));
      REAL llktry=this->GenerateNewSimplexConfiguration(vLLK,vIndex,worst,-1.0);
      if(llktry<=vLLK(best))
         llktry=this->GenerateNewSimplexConfiguration(vLLK,vIndex,worst,2.0);
      else if(llktry>=vLLK(nextworst))
      {
         const REAL llksave=vLLK(worst);
         llktry=this->GenerateNewSimplexConfiguration(vLLK,vIndex,worst,0.5);
         if(llktry>=llksave)
         {
            CrystVector_REAL *p0=&(mRefParList.GetParamSet(vIndex(best)));
            for(unsigned long i=0;i<=n;i++)
            {
               if(i==best) continue;
               CrystVector_REAL *pi=&(mRefParList.GetParamSet(vIndex(i)));
               for(unsigned long j=0;j<n;j++)
                  {(*pi)(j) = 0.5*((*p0)(j) + (*pi)(j));}
               mRefParList.RestoreParamSet(vIndex(i));
               vLLK(i)=this->GetLogLikelihood();
            }
         }
      }
   }
   mRefParList.RestoreParamSet(vIndex(best));
   mRefParList.EraseAllParamSet();
   for(int i=0;i<mRefinedObjList.GetNb();i++) mRefinedObjList.GetObj(i).EndOptimization();
   VFN_DEBUG_EXIT("SimplexObj::Optimize()",10)
}
void SimplexObj::MultiRunOptimize(long &nbCycle,long &nbSteps,const bool silent,
                                    const REAL finalcost,const REAL maxTime)
{
   const long nbStep0=nbSteps;
   while(nbCycle--!=0)
   {
      if(!silent) cout <<"SimplexObj::MultiRunOptimize: Starting Run#"<<abs(nbCycle)<<endl;
      nbSteps=nbStep0;
      this->Optimize(nbSteps,silent,finalcost,maxTime);
      if(!silent) cout <<"SimplexObj::MultiRunOptimize: Finished Run#"<<abs(nbCycle)<<endl;
   }
}
void SimplexObj::XMLOutput(ostream &os,int indent)const
{
}
void SimplexObj::XMLInput(istream &is,const XMLCrystTag &tag)
{
}
REAL SimplexObj::GenerateNewSimplexConfiguration(CrystVector_REAL &vLLK,
                                                 CrystVector_long &vIndex,
                                                 unsigned long worst,
                                                 REAL f)
{
   // Compute center of vertices
   const unsigned long n=mRefParList.GetNbParNotFixed();
   const unsigned long iworst=vIndex(worst);
   CrystVector_REAL center(n);
   center = 0.0;
   for(unsigned long i=0;i<=n;i++)
   {
      if(i==worst) continue;
      center += mRefParList.GetParamSet(vIndex(i));
   }
   center /= (REAL)n;
   for(unsigned long i=0;i<n;i++)
   {
      mRefParList.GetParNotFixed(i)
         .MutateTo(center(i)+f*(mRefParList.GetParamSet(iworst)(i)-center(i)));
   }
   const REAL llk=this->GetLogLikelihood();
   //cout<<"New,f="<<f<<",worst="<<worst<<":"<<llk<<endl;
   if(llk<vLLK(worst))
   {
      vLLK(worst)=llk;
      mRefParList.SaveParamSet(iworst);
   }
   return llk;
}
#ifdef __WX__CRYST__
WXCrystObjBasic* SimplexObj::WXCreate(wxWindow*)
{
   return NULL;
}
WXOptimizationObj* SimplexObj::WXGet()
{
   return NULL;
}
void SimplexObj::WXDelete(){}
void SimplexObj::WXNotifyDelete(){}
#endif
}//namespace
