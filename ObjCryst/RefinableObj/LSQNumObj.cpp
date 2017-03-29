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
#include "ObjCryst/Quirks/VFNStreamFormat.h"

#include "ObjCryst/RefinableObj/LSQNumObj.h"

#ifdef __WX__CRYST__
   #include "ObjCryst/wxCryst/wxLSQ.h"
#endif

#include "newmat/newmatap.h" //for SVD decomposition
#include "newmat/newmatio.h"

#ifdef use_namespace
using namespace NEWMAT;
#endif
using namespace std;

#include <iomanip>

#define POSSIBLY_UNUSED(expr) (void)(expr)

namespace ObjCryst
{

LSQNumObj::LSQNumObj(string objName)
#ifdef __WX__CRYST__
:mpWXCrystObj(0)
#endif
{
   mDampingFactor=1.;
   mSaveReportOnEachCycle=false;
   mName=objName;
   mSaveFileName="LSQrefinement.save";
   mR=0;
   mRw=0;
   mChiSq=0;
   mStopAfterCycle=false;
}

LSQNumObj::~LSQNumObj()
{
   #ifdef __WX__CRYST__
   this->WXDelete();
   #endif
}

void LSQNumObj::SetParIsFixed(const string& parName,const bool fix)
{
   if(mRefParList.GetNbPar()==0) this->PrepareRefParList();
   mRefParList.SetParIsFixed(parName,fix);
}
void LSQNumObj::SetParIsFixed(const RefParType *type,const bool fix)
{
   if(mRefParList.GetNbPar()==0) this->PrepareRefParList();
   mRefParList.SetParIsFixed(type,fix);
}

void LSQNumObj::SetParIsFixed(RefinablePar &par,const bool fix)
{
   if(mRefParList.GetNbPar()==0) this->PrepareRefParList();
   mRefParList.GetPar(par.GetPointer()).SetIsFixed(fix);
}

void LSQNumObj::SetParIsFixed(RefinableObj &obj,const bool fix)

{
   if(mRefParList.GetNbPar()==0) this->PrepareRefParList();
   for(unsigned int i=0;i<obj.GetNbPar();++i)
      this->SetParIsFixed(obj.GetPar(i),fix);
}

void LSQNumObj::UnFixAllPar()
{
   if(mRefParList.GetNbPar()==0) this->PrepareRefParList();
   mRefParList.UnFixAllPar();
}

void LSQNumObj::SetParIsUsed(const string& parName,const bool use)
{
   if(mRefParList.GetNbPar()==0) this->PrepareRefParList();
   mRefParList.SetParIsUsed(parName,use);
}
void LSQNumObj::SetParIsUsed(const RefParType *type,const bool use)
{
   if(mRefParList.GetNbPar()==0) this->PrepareRefParList();
   mRefParList.SetParIsUsed(type,use);
}

void LSQNumObj::Refine (int nbCycle,bool useLevenbergMarquardt,
                        const bool silent, const bool callBeginEndOptimization,
                        const float minChi2var)
{
   TAU_PROFILE("LSQNumObj::Refine()","void ()",TAU_USER);
   TAU_PROFILE_TIMER(timer1,"LSQNumObj::Refine() 1 - Init","", TAU_FIELD);
   TAU_PROFILE_TIMER(timer2,"LSQNumObj::Refine() 2 - LSQ Deriv","", TAU_FIELD);
   TAU_PROFILE_TIMER(timer3,"LSQNumObj::Refine() 3 - LSQ MB","", TAU_FIELD);
   TAU_PROFILE_TIMER(timer4,"LSQNumObj::Refine() 4 - LSQ Singular Values","", TAU_FIELD);
   TAU_PROFILE_TIMER(timer5,"LSQNumObj::Refine() 5 - LSQ Newmat, eigenvalues...","", TAU_FIELD);
   TAU_PROFILE_TIMER(timer6,"LSQNumObj::Refine() 6 - LSQ Apply","", TAU_FIELD);
   TAU_PROFILE_TIMER(timer7,"LSQNumObj::Refine() 7 - LSQ Finish","", TAU_FIELD);
   TAU_PROFILE_START(timer1);
   if(callBeginEndOptimization) this->BeginOptimization();
   mObs=this->GetLSQObs();
   mWeight=this->GetLSQWeight();

   bool terminateOnDeltaChi2=false;
   if(nbCycle<0)
   {
      nbCycle=-nbCycle;
      terminateOnDeltaChi2=true;
   }

   if(!silent) cout << "LSQNumObj::Refine():Beginning "<<endl;
   //Prepare for refinement (get non-fixed parameters)
      if(mRefParList.GetNbPar()==0) this->PrepareRefParList();
      mRefParList.PrepareForRefinement();
      //if(!silent) mRefParList.Print();
      if(mRefParList.GetNbPar()==0) throw ObjCrystException("LSQNumObj::Refine():no parameter to refine !");

   //variables
      long nbVar=mRefParList.GetNbParNotFixed();
      const long nbObs=mObs.numElements();
      CrystVector_REAL calc,calc0,calc1,tmpV1,tmpV2;
      CrystMatrix_REAL M(nbVar,nbVar);
      CrystMatrix_REAL N(nbVar,nbVar);
      CrystVector_REAL B(nbVar);
      CrystMatrix_REAL designMatrix(nbVar,nbObs);
      CrystVector_REAL deltaVar(nbVar);
      long i,j,k;
      REAL R_ini,Rw_ini;  POSSIBLY_UNUSED(R_ini);
      REAL *pTmp1,*pTmp2;

      REAL marquardt=1e-2;
      const REAL marquardtMult=4.;
   //initial Chi^2, needed for Levenberg-Marquardt
   this->CalcChiSquare();
   //store old values
   mIndexValuesSetInitial=mRefParList.CreateParamSet("LSQ Refinement-Initial Values");
   mIndexValuesSetLast=mRefParList.CreateParamSet("LSQ Refinement-Last Cycle Values");
   TAU_PROFILE_STOP(timer1);
   //refine
   for(int cycle=1 ; cycle <=nbCycle;cycle++)
   {
      TAU_PROFILE_START(timer2);
      const REAL ChisSqPreviousCycle=mChiSq;
      mRefParList.SaveParamSet(mIndexValuesSetLast);// end of last cycle
      if(!silent) cout << "LSQNumObj::Refine():Cycle#"<< cycle <<endl;
      //initial value of function
         calc0=this->GetLSQCalc();
         //R
            tmpV1 =  mObs;
            tmpV1 -= calc0;
            tmpV1 *= tmpV1;
            tmpV2 =  mObs;
            tmpV2 *= mObs;
            R_ini=sqrt(tmpV1.sum()/tmpV2.sum());
         //Rw
            tmpV1 *= mWeight;
            tmpV2 *= mWeight;
            Rw_ini=sqrt(tmpV1.sum()/tmpV2.sum());
      //derivatives
      //designMatrix=0.;
      pTmp2=designMatrix.data();
      //cout <<"obs:"<<FormatHorizVector<REAL>(calc0,10,8);
      //cout <<"calc:"<<FormatHorizVector<REAL>(mObs,10,8);
      //cout <<"weight:"<<FormatHorizVector<REAL>(mWeight,10,8);
      #if 1
      for(i=0;i<nbVar;i++)
      {
         //:NOTE: Real design matrix is the transposed of the one computed here
         //if(!silent) cout << "........." << mRefParList.GetParNotFixed(i).GetName() <<endl;

         tmpV1=this->GetLSQDeriv(mRefParList.GetParNotFixed(i));
         pTmp1=tmpV1.data();
         //cout <<"deriv#"<<i<<":"<<FormatHorizVector<REAL>(tmpV1,10,8);
         for(j=0;j<nbObs;j++) *pTmp2++ = *pTmp1++;
      }
      #else
      this->GetLSQ_FullDeriv();
      for(i=0;i<nbVar;i++)
      {
         pTmp1=mLSQ_FullDeriv[&(mRefParList.GetParNotFixed(i))].data();
         //if(i>=(nbVar-2)) cout<<__FILE__<<":"<<__LINE__<<":"<<(mRefParList.GetParNotFixed(i)).GetName()<<"size="<<mLSQ_FullDeriv[&(mRefParList.GetParNotFixed(i))].size()<<":"<<mLSQ_FullDeriv[&(mRefParList.GetParNotFixed(i))]<<endl;
         for(j=0;j<nbObs;j++) *pTmp2++ = *pTmp1++;
      }
      #endif
         //cout << designMatrix;

      TAU_PROFILE_STOP(timer2);
      LSQNumObj_Refine_Restart: //Used in case of singular matrix
      TAU_PROFILE_START(timer3);

      //Calculate M and B matrices
         tmpV1.resize(nbObs);
         tmpV2.resize(nbObs);
         for(i=0;i<nbVar;i++)
         {
            #if 1
            {
               const register REAL * RESTRICT pD=designMatrix.data()+i*designMatrix.cols();
               const register REAL * RESTRICT pW=mWeight.data();
               REAL * RESTRICT p=tmpV1.data();
               for(k=0;k<nbObs;k++) *p++ = *pD++ * *pW++;
            }
            const register REAL * pD=designMatrix.data();
            for(j=0;j<nbVar;j++)
            {
               const register REAL * p1=tmpV1.data();
               REAL v2=0;
               for(k=0;k<nbObs;k++) v2+= *pD++ * *p1++;
               M(j,i)=v2;
            }
            REAL b=0;
            const register REAL * pObs=mObs.data();
            const register REAL * pCalc=calc0.data();
            const register REAL * p1=tmpV1.data();
            for(k=0;k<nbObs;k++) b+= (*pObs++ - *pCalc++)* *p1++;
            B(i)=b;
            #else
            for(k=0;k<nbObs;k++) tmpV1(k)=designMatrix(i,k);
            tmpV1 *= mWeight;
            for(j=0;j<nbVar;j++)
            {
               for(k=0;k<nbObs;k++) tmpV2(k)=designMatrix(j,k);
               tmpV2 *= tmpV1;
               M(j,i)= tmpV2.sum();
                  //M(i,j)=total(designMatrix.row(i)*designMatrix.row(j)*weight);
            }
            tmpV2=mObs;
            tmpV2 -= calc0;
            tmpV2 *= tmpV1;
            B(i)=tmpV2.sum();//total((obs-calc0)*weight*designMatrix.row(i))
            #endif

         }
      TAU_PROFILE_STOP(timer3);
      bool increaseMarquardt=false;
      LSQNumObj_Refine_RestartMarquardt: //Used in case of singular matrix or for Marquardt
      TAU_PROFILE_START(timer4);

       //Apply LevenBerg-Marquardt factor
         if(true==useLevenbergMarquardt)
         {
            const REAL marquardtOLD=marquardt;
            if(increaseMarquardt) marquardt=marquardt*marquardtMult;
            const REAL lmfact=(1+marquardt)/(1+marquardtOLD);
            for(i=0;i<nbVar;i++) M(i,i) *= lmfact;
         }
       // Check for singular values
         for(i=0;i<nbVar;i++)
         {
            //cout<<__FILE__<<":"<<__LINE__<<":"<<i<<":"<<mRefParList.GetParNotFixed(i).GetName()<<":"<<mRefParList.GetParNotFixed(i).GetPointer()<<":"<<M(i,i)<<endl;
            if( (M(i,i) < 1e-20)||(ISNAN_OR_INF(M(i,i)))) //:TODO: Check what value to use as a limit
            {
               if(!silent) cout << "LSQNumObj::Refine() Singular parameter !"
                                    << "(null derivate in all points) : "<<M(i,i)<<":"
                                    << mRefParList.GetParNotFixed(i).GetName() << endl;
               /*
               if(!silent)
               {
                  for(i=0;i<nbVar;i++)
                  {
                     tmpV1=this->GetLSQDeriv(mRefParList.GetParNotFixed(i));
                     cout <<"deriv#"<<i<<":"<<FormatHorizVector<REAL>(tmpV1,10,8);
                  }
               }
               */
               if(!silent) cout << "LSQNumObj::Refine(): Automatically fixing parameter";
               if(!silent) cout << " and re-start cycle..";
               if(!silent) cout << endl;
               mRefParList.GetParNotFixed(i).SetIsFixed(true);
               mRefParList.PrepareForRefinement();
               nbVar=mRefParList.GetNbParNotFixed();
               if(nbVar<=1)
               {
                  mRefParList.RestoreParamSet(mIndexValuesSetInitial);
                  if(callBeginEndOptimization) this->EndOptimization();
                  if(!silent) mRefParList.Print();
                  throw ObjCrystException("LSQNumObj::Refine(): not enough (1) parameters after fixing one...");
               }
               N.resize(nbVar,nbVar);
               deltaVar.resize(nbVar);

               //Just remove the ith line in the design matrix
                        REAL *p1=designMatrix.data();
                  const REAL *p2=designMatrix.data();
                  p1 += i*nbObs;
                  p2 += i*nbObs+nbObs;
                  for(long j=i*nbObs+nbObs;j<nbObs*nbVar;j++) *p1++ = *p2++;
               designMatrix.resizeAndPreserve(nbVar,nbObs);

               //:TODO: Make this work...
               /*
               //Remove ith line &Column in M & ith element in B
                  p1=M.data();
                  p2=M.data();
                  for(long j=0;j<=nbVar;j++)
                  {
                     if( (j>=i) && (j<nbVar) ) B(j)=B(j+1);
                     for(long k=0;k<=nbVar;k++)
                     {
                        if((j==i) || (k==i)) p2++;
                        else *p1++ = *p2++;
                     }
                  }
               M.resizeAndPreserve(nbVar,nbVar);
               B.resizeAndPreserve(nbVar);
               */
               M.resize(nbVar,nbVar);
               B.resize(nbVar);
               TAU_PROFILE_STOP(timer4);
               goto LSQNumObj_Refine_Restart;
            }
         }
      TAU_PROFILE_STOP(timer4);

/*
      //Solve with Singular value Decomposition on design matrix (using newmat)
      {
         cout << "LSQNumObj::Refine():Performing SVD on design matrix..." << endl;
         Matrix newmatA(nbObs,nbVar), newmatU(nbObs,nbVar), newmatV(nbVar,nbVar);
         DiagonalMatrix newmatW(nbVar);
         ColumnVector newmatB(nbObs);

         CrystMatrix_REAL U(nbVar,nbObs), V(nbVar,nbVar);
         CrystVector_REAL W(nbVar),invW(nbVar);
         for(long i=0;i<nbVar;i++)
         {
            for(long j=0;j<nbObs;j++)
               newmatA(j+1,i+1) = designMatrix(i,j);//design matrix (need to transpose)
         }

         for(long j=0;j<nbObs;j++) newmatB(j+1)=mObs(j)-calc0(j);

         SVD(newmatA,newmatW,newmatU,newmatV);

         ColumnVector newmatDelta(nbVar);
         DiagonalMatrix newmatInvW(nbVar);
         REAL max=newmatW.MaximumAbsoluteValue();
         REAL minAllowedValue=1e-16*max;// :TODO: Check if reasonable !
         //Avoid singular values,following 'Numerical Recipes in C'
         for(long i=0;i<nbVar;i++)
            if(newmatW(i+1,i+1) > minAllowedValue) newmatInvW(i+1,i+1) = 1./newmatW(i+1,i+1);
            else
            {
               cout << "LSQNumObj::Refine():SVD: fixing singular value "<< i <<endl;
               newmatInvW(i+1,i+1) = 0.;
            }
         newmatDelta=newmatV * newmatInvW * newmatU.t() * newmatB;

         for(long i=0;i<nbVar;i++)
         {
            W(i)=newmatW(i+1,i+1);
            invW(i)=newmatInvW(i+1,i+1);
            deltaVar(i)=newmatDelta(i+1);
         }
         cout << newmatW << endl;
      }
*/
      TAU_PROFILE_START(timer5);
      //Perform "Eigenvalue Filtering" on normal matrix (using newmat library)
      {
         //if(!silent) cout << "LSQNumObj::Refine():Eigenvalue Filtering..." <<endl;
         CrystMatrix_REAL V(nbVar,nbVar);
         CrystVector_REAL invW(nbVar);//diagonal matrix, in fact
         //if(!silent) cout << "LSQNumObj::Refine():Eigenvalue Filtering...1" <<endl;
         {
            SymmetricMatrix newmatA(nbVar);
            Matrix   newmatV(nbVar,nbVar),
                     newmatN(nbVar,nbVar);
            DiagonalMatrix newmatW(nbVar);
            ColumnVector newmatB(nbVar);
            //'Derivative scaling' matrix
               DiagonalMatrix newmatDscale(nbVar);
               for(long i=0;i<nbVar;i++)
                  newmatDscale(i+1,i+1) = 1./sqrt(M(i,i));
            for(long i=0;i<nbVar;i++)
            {
               newmatB(i+1)=B(i);//NOT scaled
               for(long j=0;j<nbVar;j++)
                  newmatA(i+1,j+1) = M(i,j) * newmatDscale(i+1,i+1) * newmatDscale(j+1,j+1);
            }
         //if(!silent) cout << "LSQNumObj::Refine():Eigenvalue Filtering...2" <<endl;

            //cout << newmatA.SubMatrix(1,3,1,3) <<endl;
         //if(!silent) cout << "LSQNumObj::Refine():Eigenvalue Filtering...3" <<endl;

            //Jacobi(newmatA,newmatW,newmatV);
            try
            {
               EigenValues(newmatA,newmatW,newmatV);
            }
            catch(...)
            {
               cout<<"Caught a Newmat exception :"<<BaseException::what()<<endl;
               cout<<"A:"<<endl<<newmatA<<endl<<"W:"<<endl<<newmatW<<endl<<"V:"<<endl<<newmatV<<endl<<"Dscale:"<<newmatDscale*1e6<<endl;
               cout<<setw(5)<<"B:"<<endl;
               for(unsigned int i=0;i<B.size();i++) cout<<B(i)<<" ";
               cout<<endl<<endl<<"M:"<<endl;
               for(unsigned int i=0;i<M.rows();i++)
               {
                  for(unsigned int j=0;j<M.cols();j++) cout<<M(i,j)<<" ";
                  cout<<endl;
               }
               cout<<endl<<endl<<"D("<<designMatrix.rows()<<"x"<<designMatrix.cols()<<"):"<<endl;
               for(unsigned int i=0;i<designMatrix.rows();i++)
               {
                  for(unsigned int j=0;j<designMatrix.cols();j++) cout<<designMatrix(i,j)<<" ";
                  cout<<endl;
               }
               throw ObjCrystException("LSQNumObj::Refine():caught a newmat exception during Eigenvalues computing !");
            }
               ColumnVector newmatDelta(nbVar);
               DiagonalMatrix newmatInvW(nbVar);
        //if(!silent) cout << "LSQNumObj::Refine():Eigenvalue Filtering...4" <<endl;
            //Avoid singular values
            {
               REAL max=newmatW.MaximumAbsoluteValue();
               REAL minAllowedValue=1e-5*max;// :TODO: Check if reasonable !
               for(long i=0;i<nbVar;i++)
                  if(newmatW(i+1,i+1) > minAllowedValue)
                     newmatInvW(i+1,i+1)= 1./newmatW(i+1,i+1);
                  else
                  {
                     if(!silent) cout << "LSQNumObj::Refine():fixing ill-cond EigenValue "<< i <<endl;
                     newmatInvW(i+1,i+1) = 0.;
                  }

               /*
               REAL max;
               int maxIndex;
               const REAL minRatio=1e-6;
               for(long i=0;i<nbVar;i++)
               {
                  max=fabs(newmatW(1)*newmatV(1,i+1)*newmatV(1,i+1));
                  maxIndex=0;
                  for(long j=1;j<nbVar;j++)
                     if( fabs(newmatW(j+1)*newmatV(j+1,i+1)*newmatV(j+1,i+1)) > max )
                     {
                        max = fabs(newmatW(j+1)*newmatV(j+1,i+1)*newmatV(j+1,i+1));
                        maxIndex=j;
                     }

                  cout << FormatFloat(newmatW(i+1,i+1),36,2) << " " ;
                  cout << FormatFloat(max,36,2);
                  cout << FormatFloat(newmatW(maxIndex+1,maxIndex+1),36,2) << " " ;
                  cout << FormatFloat(newmatV(i+1,maxIndex+1)*100,8,2) ;
                  //cout << endl;

                  if(newmatW(i+1,i+1) > (max*minRatio))
                  {
                     newmatInvW(i+1,i+1) = 1./newmatW(i+1,i+1);
                     cout << endl;
                  }
                  else
                  {
                     cout << "LSQNumObj::Refine():fixing ill-cond EigenValue "<< i <<endl;
                     newmatInvW(i+1,i+1) = 0.;
                  }
               }
               */
            }
            newmatN=newmatV * newmatInvW * newmatV.t();

            //Back 'Derivative Scaling' :TODO:
               newmatN = newmatDscale * newmatN * newmatDscale;
            newmatDelta = newmatN * newmatB;

            for(long i=0;i<nbVar;i++)
            {
               invW(i)=newmatInvW(i+1,i+1);
               deltaVar(i)=newmatDelta(i+1);
               for(long j=0;j<nbVar;j++)
               {
                  N(i,j) = newmatN(i+1,j+1);
                  V(i,j) = newmatV(i+1,j+1);
               }
            }
            //cout << invW <<endl << U << endl << V << endl << B << endl << deltaVar << endl;
            //cout<<setw(10)<<setprecision(4)<< newmatV*100. << endl;
            //cout << newmatA/1e15 << endl;
            //cout << newmatW << endl;
            //cout << newmatDscale*1e20 << endl;
            //cout << newmatB <<endl;
         }
      }//End EigenValue filtering
      TAU_PROFILE_STOP(timer5);
/*
*/
/*
      //Perform Singular value Decomposition normalon normal matrix (using newmat library)
      {
         cout << "LSQNumObj::Refine():Performing SVD on normal matrix" <<endl;
         CrystMatrix_REAL U(nbVar,nbVar),V(nbVar,nbVar);
         CrystVector_REAL invW(nbVar);//diagonal matrix, in fact
         {
            Matrix   newmatA(nbVar,nbVar),
                     newmatU(nbVar,nbVar),
                     newmatV(nbVar,nbVar),
                     newmatN(nbVar,nbVar);
            DiagonalMatrix newmatW(nbVar);
            ColumnVector newmatB(nbVar);
            for(long i=0;i<nbVar;i++)
            {
               newmatB(i+1)=B(i);
               for(long j=0;j<nbVar;j++) newmatA(i+1,j+1) = M(i,j);
            }
            SVD(newmatA,newmatW,newmatU,newmatV);
            ColumnVector newmatDelta(nbVar);
            DiagonalMatrix newmatInvW(nbVar);
            REAL max=newmatW.MaximumAbsoluteValue();
            REAL minAllowedValue=1e-10*max;// :TODO: Check if reasonable !
            //Avoid singular values,following 'Numerical Recipes in C'
            for(long i=0;i<nbVar;i++)
               if(newmatW(i+1,i+1) > minAllowedValue) newmatInvW(i+1,i+1) = 1./newmatW(i+1,i+1);
               else
               {
                  cout << "LSQNumObj::Refine():SVD: fixing ill-cond value "<< i <<endl;
                  newmatInvW(i+1,i+1) = 0.;
               }
            newmatN=newmatV * newmatInvW * newmatU.t();
            newmatDelta = newmatN * newmatB;


            for(long i=0;i<nbVar;i++)
            {
               invW(i)=newmatInvW(i+1,i+1);
               deltaVar(i)=newmatDelta(i+1);
               for(long j=0;j<nbVar;j++)
               {
                  N(i,j) = newmatN(i+1,j+1);
                  U(i,j) = newmatU(i+1,j+1);
                  V(i,j) = newmatV(i+1,j+1);
               }
            }
            //cout << invW <<endl << U << endl << V << endl << B << endl << deltaVar << endl;
         }
        }//End SVD
*/
/*
      //OLD VERSION USING GAUS INVERSION
      //Calculate new values for variables
      cout << "LSQNumObj::Refine():Computing new values for variables" <<endl;
         try { N=InvertMatrix(M);}
         catch(long svix)
         {
            cout << "LSQNumObj::Refine(): WARNING : Singular Matrix !!!";
            cout << "With refinable parameter: " << mRefParList.GetParNotFixed(svix).Name() <<endl;
            cout << "LSQNumObj::Refine(): Automatically fixing parameter and re-start cycle..";              cout << endl;
            {
               CrystMatrix_REAL Mtmp;
               Mtmp=M;
               Mtmp /= 1e14;
               cout << Mtmp << endl << endl;
            }
            mRefParList.GetParNotFixed(svix).IsFixed()=true;
            mRefParList.PrepareForRefinement();
            nbVar=mRefParList.NbGetParNotFixed();
            M.resize(nbVar,nbVar);
            N.resize(nbVar,nbVar);
            B.resize(nbVar);
            designMatrix.resize(nbVar,nbObs);
            deltaVar.resize(nbVar);
            goto LSQNumObj_Refine_Restart;
         }
         deltaVar=0;
         for(i=0;i<nbVar;i++)
            for(j=0;j<nbVar;j++) deltaVar(i) += N(i,j) * B(j);
*/
         /*
         {
            CrystMatrix_REAL Mtmp,Ntmp;
            Mtmp=M;
            Ntmp=N;
            Mtmp /= 1e20;
            Ntmp *= 1e20;
            //cout << Mtmp << endl << endl;
            cout << Ntmp << endl << endl;
            //cout << B <<endl;
            //cout << deltaVar  << endl;
         }
         */
      TAU_PROFILE_START(timer6);
      /// Applying new computed values :TODO: & Check if a limit has been hit
         for(i=0;i<nbVar;i++)
         {
            //const REAL oldvalue=mRefParList.GetParNotFixed(i).GetValue();
            //const REAL expected=oldvalue+deltaVar(i);
            mRefParList.GetParNotFixed(i).Mutate(deltaVar(i));
            //const REAL newvalue=mRefParList.GetParNotFixed(i).GetValue();
         }

      //for statistics...
         //mRefParList.Print();
         calc=this->GetLSQCalc();
         //Chi^2
         {
            REAL oldChiSq=mChiSq;
            tmpV1 = mObs;
            tmpV1 -= calc;
            tmpV1 *= tmpV1;
            tmpV1 *= mWeight;
            mChiSq=tmpV1.sum();
            if(true==useLevenbergMarquardt)
            {
               if(mChiSq > (oldChiSq*1.0001))
               {
                  mRefParList.RestoreParamSet(mIndexValuesSetLast);
                  increaseMarquardt=true;
                  if(!silent)
                  {
                     cout << "LSQNumObj::Refine(Chi^2="<<oldChiSq<<"->"<<mChiSq
                          <<")=>Increasing Levenberg-Marquardt factor :"
                          << FormatFloat(marquardt*marquardtMult,18,14) <<endl;
                  }
                  mChiSq=oldChiSq;
                  if(marquardt>1e4)
                  {
                     // :TODO: Revert to previous parameters. Or initial ?
                     mRefParList.RestoreParamSet(mIndexValuesSetLast);
                     if(callBeginEndOptimization) this->EndOptimization();
                     //if(!silent) mRefParList.Print();
                     return;
                     //throw ObjCrystException("LSQNumObj::Refine():Levenberg-Marquardt diverging !");
                  }
                  TAU_PROFILE_STOP(timer6);
                  goto LSQNumObj_Refine_RestartMarquardt;
               }
               else
               {
                  if(!silent && (marquardt>1e-2))
                  {
                     cout << "LSQNumObj::Refine(Chi^2="<<oldChiSq<<"->"<<mChiSq
                          <<")=>Decreasing Levenberg-Marquardt factor :" << FormatFloat(marquardt/marquardtMult,18,14) <<endl;
                  }
                  marquardt /= marquardtMult;
                  if(marquardt<1e-2) marquardt=1e-2;
               }
            }
         }
      TAU_PROFILE_STOP(timer6);
      TAU_PROFILE_START(timer7);

         //Sigmas
            if(nbObs==nbVar)
               for(i=0;i<nbVar;i++) mRefParList.GetParNotFixed(i).SetSigma(0);
            else
               for(i=0;i<nbVar;i++) mRefParList.GetParNotFixed(i).SetSigma(sqrt(N(i,i)*mChiSq/(nbObs-nbVar)));
         //Correlations :TODO: re-compute M and N if using Levenberg-Marquardt
         mCorrelMatrix.resize(nbVar,nbVar);
         mvVarCovar.clear();

         for(i=0;i<nbVar;i++)
         {
            RefinablePar *pi=&(mRefParList.GetParNotFixed(i));
            for(j=0;j<nbVar;j++)
            {
               RefinablePar *pj=&(mRefParList.GetParNotFixed(j));
               mCorrelMatrix(i,j)=sqrt(N(i,j)*N(i,j)/N(i,i)/N(j,j));
               if(nbObs!=nbVar)
               mvVarCovar[make_pair(pi,pj)]=N(i,j)*mChiSq/(nbObs-nbVar);
            }
         }
         //R-factor
            tmpV1 = mObs;
            tmpV1 -= calc;
            tmpV1 *= tmpV1;
            tmpV2 = mObs;
            tmpV2 *= tmpV2;
            mR=sqrt(tmpV1.sum()/tmpV2.sum());
         //Rw-factor
            tmpV1 *= mWeight;
            tmpV2 *= mWeight;
            mRw=sqrt(tmpV1.sum()/tmpV2.sum());
      //OK, finished
         if(!silent) cout << "finished cycle #"<<cycle <<"/"<<nbCycle <<". Rw="<<Rw_ini<<"->"<<mRw<<",    Chi^2="<<ChisSqPreviousCycle<<"->"<<mChiSq<<endl;
         if (mSaveReportOnEachCycle) this->WriteReportToFile();

      if(!silent) this->PrintRefResults();
      TAU_PROFILE_STOP(timer7);
      if( terminateOnDeltaChi2 && (minChi2var>( (ChisSqPreviousCycle-mChiSq)/abs(ChisSqPreviousCycle+1e-6) ) ) ) break;
   }
   if(callBeginEndOptimization) this->EndOptimization();
}

bool LSQNumObj::SafeRefine(std::list<RefinablePar*> vnewpar, std::list<const RefParType*> vnewpartype,
                                 REAL maxChi2factor,
                                 int nbCycle, bool useLevenbergMarquardt,
                                 const bool silent, const bool callBeginEndOptimization,
                                 const float minChi2var)
{
   if(callBeginEndOptimization) this->BeginOptimization();
   // :TODO: update mObs and mWeight in a centralized way... Not in BeginOptimization() (not always called)
   mObs=this->GetLSQObs();
   mWeight=this->GetLSQWeight();
   
   //Prepare for refinement (get non-fixed parameters)
   if(mRefParList.GetNbPar()==0) this->PrepareRefParList();
   mRefParList.PrepareForRefinement();
   if(mRefParList.GetNbPar()==0) throw ObjCrystException("LSQNumObj::SafeRefine():no parameter to refine !");

   this->CalcChiSquare();
   const REAL chi2_0 = mChiSq;
   for(std::list<RefinablePar*>::iterator pos=vnewpar.begin(); pos!=vnewpar.end(); pos++)
   {
      this->SetParIsFixed(**pos, false);
   }
   for(std::list<const RefParType*>::iterator pos=vnewpartype.begin(); pos!=vnewpartype.end(); pos++)
   {
      this->SetParIsFixed(*pos, false);
   }
   bool diverged = false;
   try
   {
      this->Refine(nbCycle, useLevenbergMarquardt, silent, false, minChi2var);
   }
   catch(const ObjCrystException &except)
   {
      diverged = true;
      cout << "Refinement did not converge !";
   }
   const REAL deltachi2 = (mChiSq-chi2_0)/(chi2_0+1e-6);
   if(callBeginEndOptimization) this->EndOptimization();
   if(deltachi2>maxChi2factor)
   {
      cout << "Refinement did not converge ! Chi2 increase("<<chi2_0<<"->"<<mChiSq<<") by a factor: "<< deltachi2<<endl;
      diverged = true;
   }
   if(diverged)
   {
      mRefParList.RestoreParamSet(mIndexValuesSetInitial);
      for(std::list<RefinablePar*>::iterator pos=vnewpar.begin(); pos!=vnewpar.end(); pos++)
      {
         this->SetParIsFixed(**pos, true);
      }
      for(std::list<const RefParType*>::iterator pos=vnewpartype.begin(); pos!=vnewpartype.end(); pos++)
      {
         this->SetParIsFixed(*pos, true);
      }
      this->CalcRfactor();
      this->CalcRwFactor();
      this->CalcChiSquare();
      cout <<"=> REVERTING to initial parameters values and fixing new parameters"<<endl;
      return false;
   }
   return true;
}

CrystMatrix_REAL LSQNumObj::CorrelMatrix()const{return mCorrelMatrix;};

void LSQNumObj::CalcRfactor()const
{
   CrystVector_REAL calc, tmpV1, tmpV2;
   calc=this->GetLSQCalc();
   tmpV1 = this->GetLSQObs();
   tmpV1 -= calc;
   tmpV1 *= tmpV1;
   tmpV2 = this->GetLSQObs();
   tmpV2 *= tmpV2;
   mR=sqrt(tmpV1.sum()/tmpV2.sum());
}

REAL LSQNumObj::Rfactor()const{return mR;};

void LSQNumObj::CalcRwFactor()const
{
   CrystVector_REAL calc, tmpV1, tmpV2;
   calc=this->GetLSQCalc();
   tmpV1 = this->GetLSQObs();
   tmpV1 -= calc;
   tmpV1 *= tmpV1;
   tmpV2 = this->GetLSQObs();
   tmpV2 *= tmpV2;
   tmpV1 *= mWeight;
   tmpV2 *= mWeight;
   mRw=sqrt(tmpV1.sum()/tmpV2.sum());
}

REAL LSQNumObj::RwFactor()const{return mRw;};

void LSQNumObj::CalcChiSquare()const
{
   CrystVector_REAL calc, tmpV1;
   calc=this->GetLSQCalc();
   tmpV1 = mObs;
   tmpV1 -= calc;
   tmpV1 *= tmpV1;
   tmpV1 *= mWeight;
   mChiSq=tmpV1.sum();
}

REAL LSQNumObj::ChiSquare()const{return mChiSq;};


void RecursiveMapFunc(RefinableObj &obj,map<RefinableObj*,unsigned int> &themap, const unsigned int value)
{
   themap[&obj]=value;
   ObjRegistry<RefinableObj> *pObjReg=&(obj.GetSubObjRegistry());
   for(int i=0;i<pObjReg->GetNb();i++)
      RecursiveMapFunc(pObjReg->GetObj(i),themap,value);
   return;
}

void LSQNumObj::SetRefinedObj(RefinableObj &obj, const unsigned int LSQFuncIndex, const bool init, const bool recursive)

{
   if(init)
   {
      mvRefinedObjMap.clear();
   }
   if(recursive) RecursiveMapFunc(obj,mvRefinedObjMap,LSQFuncIndex);
   else mvRefinedObjMap[&obj]=LSQFuncIndex;
}

//ObjRegistry<RefinableObj> &LSQNumObj::GetRefinedObjList(){return mRecursiveRefinedObjList;}

const map<RefinableObj*,unsigned int>& LSQNumObj::GetRefinedObjMap() const
{
   return mvRefinedObjMap;
}

map<RefinableObj*,unsigned int>& LSQNumObj::GetRefinedObjMap()
{
   return mvRefinedObjMap;
}


RefinableObj& LSQNumObj::GetCompiledRefinedObj(){return mRefParList;}

const RefinableObj& LSQNumObj::GetCompiledRefinedObj()const{return mRefParList;}

void LSQNumObj::SetUseSaveFileOnEachCycle(bool yesOrNo)
{
   mSaveReportOnEachCycle=yesOrNo;
}

void LSQNumObj::SetSaveFile(string fileName)
{
   mSaveFileName=fileName;
}

void LSQNumObj::PrintRefResults() const
{
   //:TODO:
   //this->PrepareRefParList(); activate this when PrepareRefParList() will be more savy
   cout << "Results after last refinement :(" ;
   cout << mRefParList.GetNbParNotFixed()<< " non-fixed parameters)"<<endl;
   cout << "Variable information : Initial, last cycle , current values and sigma"<<endl;
   for (int i=0;i<mRefParList.GetNbPar();i++)
   {
      if( (true==mRefParList.GetPar(i).IsFixed())
            || (false == mRefParList.GetPar(i).IsUsed()) ) continue;
      cout << FormatString(mRefParList.GetPar(i).GetName(),30) << "  " ;
      cout << FormatFloat((mRefParList.GetParamSet(mIndexValuesSetInitial))(i)*mRefParList.GetPar(i).GetHumanScale(),15,8) << "  " ;
      cout << FormatFloat((mRefParList.GetParamSet(mIndexValuesSetLast))(i)*mRefParList.GetPar(i).GetHumanScale(),15,8) << "  " ;
      cout << FormatFloat(mRefParList.GetPar(i).GetHumanValue(),15,8) << "  " ;
      cout << FormatFloat(mRefParList.GetPar(i).GetHumanSigma(),15,8) << "  " ;
      //cout << FormatFloat(mRefParList(i).DerivStep(),16,12) << "  ";
      //cout << varNames[i] << "  " << var0(i) << "  " << varLast(i)
      //cout<< "  " << varCurrent(i)<< "  " << sigmaValues(i)<<endl;
      cout << endl;
   }
   cout << "R-factor  : " << mR<<endl;
   cout << "Rw-factor : " << mRw<<endl;
   cout << "Chi-Square: " << mChiSq<<endl;
   cout << "GoF: " << mChiSq/this->GetLSQWeight().numElements()<<endl;
   cout <<endl;
}

void LSQNumObj::SetDampingFactor(const REAL newDampFact)
{
   mDampingFactor=newDampFact;
}

void LSQNumObj::PurgeSaveFile()
{
   //:TODO:
}

void LSQNumObj::WriteReportToFile() const
{
   //:TODO:
}

void LSQNumObj::OptimizeDerivativeSteps()
{
   //:TODO:
}

const std::map<pair<const RefinablePar*,const RefinablePar*>,REAL > & LSQNumObj::GetVarianceCovarianceMap()const
{ return mvVarCovar;}

void LSQNumObj::PrepareRefParList(const bool copy_param)
{
   mRefParList.ResetParList();
   for(map<RefinableObj*,unsigned int>::iterator pos=mvRefinedObjMap.begin();pos!=mvRefinedObjMap.end();++pos)
   {
      VFN_DEBUG_MESSAGE("LSQNumObj::PrepareRefParList():"<<pos->first->GetName(),4);
      //mRecursiveRefinedObjList.GetObj(i).Print();
      mRefParList.AddPar(*(pos->first),copy_param);
   }
   //mRefParList.Print();
   if(copy_param) mRefParList.SetDeleteRefParInDestructor(true);
   else mRefParList.SetDeleteRefParInDestructor(false);
}

const CrystVector_REAL& LSQNumObj::GetLSQCalc() const
{
   const CrystVector_REAL *pV;
   long nb=0;
   for(map<RefinableObj*,unsigned int>::const_iterator pos=mvRefinedObjMap.begin();pos!=mvRefinedObjMap.end();++pos)
   {
      if(pos->first->GetNbLSQFunction()==0) continue;
      pV=&(pos->first->GetLSQCalc(pos->second));
      const long n2 = pV->numElements();
      if((nb+n2)>mLSQCalc.numElements()) mLSQCalc.resizeAndPreserve(nb+pV->numElements());
      const REAL *p1=pV->data();
      REAL *p2=mLSQCalc.data()+nb;
      for(long j = 0; j < n2; ++j)  *p2++ = *p1++;
      nb+=n2;
   }
   if(mLSQCalc.numElements()>nb) mLSQCalc.resizeAndPreserve(nb);
   return mLSQCalc;
}

const CrystVector_REAL& LSQNumObj::GetLSQObs() const
{
   const CrystVector_REAL *pV;
   long nb=0;
   for(map<RefinableObj*,unsigned int>::const_iterator pos=mvRefinedObjMap.begin();pos!=mvRefinedObjMap.end();++pos)
   {
      if(pos->first->GetNbLSQFunction()==0) continue;
      pV=&(pos->first->GetLSQObs(pos->second));
      const long n2 = pV->numElements();
      mvRefinedObjLSQSize[pos->first]=n2;
      if((nb+n2)>mLSQObs.numElements()) mLSQObs.resizeAndPreserve(nb+pV->numElements());
      const REAL *p1=pV->data();
      REAL *p2=mLSQObs.data()+nb;
      for(long j = 0; j < n2; ++j)  *p2++ = *p1++;
      nb+=n2;
   }
   if(mLSQObs.numElements()>nb) mLSQObs.resizeAndPreserve(nb);
   return mLSQObs;
}

const CrystVector_REAL& LSQNumObj::GetLSQWeight() const
{
   const CrystVector_REAL *pV;
   long nb=0;
   for(map<RefinableObj*,unsigned int>::const_iterator pos=mvRefinedObjMap.begin();pos!=mvRefinedObjMap.end();++pos)
   {
      if(pos->first->GetNbLSQFunction()==0) continue;
      pV=&(pos->first->GetLSQWeight(pos->second));
      const long n2 = pV->numElements();
      if((nb+n2)>mLSQWeight.numElements()) mLSQWeight.resizeAndPreserve(nb+pV->numElements());
      const REAL *p1=pV->data();
      REAL *p2=mLSQWeight.data()+nb;
      for(long j = 0; j < n2; ++j)  *p2++ = *p1++;
      nb+=n2;
   }
   if(mLSQWeight.numElements()>nb) mLSQWeight.resizeAndPreserve(nb);
   return mLSQWeight;
}

const CrystVector_REAL& LSQNumObj::GetLSQDeriv(RefinablePar&par)
{
   const CrystVector_REAL *pV;
   long nb=0;
   for(map<RefinableObj*,unsigned int>::iterator pos=mvRefinedObjMap.begin();pos!=mvRefinedObjMap.end();++pos)
   {
      if(pos->first->GetNbLSQFunction()==0) continue;
      pV=&(pos->first->GetLSQDeriv(pos->second,par));
      const long n2 = pV->numElements();
      if((nb+n2)>mLSQDeriv.numElements()) mLSQDeriv.resizeAndPreserve(nb+pV->numElements());
      const REAL *p1=pV->data();
      REAL *p2=mLSQDeriv.data()+nb;
      for(long j = 0; j < n2; ++j)  *p2++ = *p1++;
      nb+=n2;
   }
   if(mLSQDeriv.numElements()>nb) mLSQDeriv.resizeAndPreserve(nb);
   return mLSQDeriv;
}

const std::map<RefinablePar*,CrystVector_REAL>& LSQNumObj::GetLSQ_FullDeriv()
{
   long nbVar=mRefParList.GetNbParNotFixed();
   std::set<RefinablePar*> vPar;
   for(unsigned int i=0;i<nbVar;++i)
      vPar.insert(&(mRefParList.GetParNotFixed(i)));
   mLSQ_FullDeriv.clear();
   unsigned long nb=0;// full length of derivative vector
   for(map<RefinableObj*,unsigned int>::iterator pos=mvRefinedObjMap.begin();pos!=mvRefinedObjMap.end();++pos)
   {
      if(pos->first->GetNbLSQFunction()==0) continue;
      const unsigned long n2=mvRefinedObjLSQSize[pos->first];
      if(n2==0) continue;//this object does not have an LSQ function

      const std::map<RefinablePar*,CrystVector_REAL> *pvV=&(pos->first->GetLSQ_FullDeriv(pos->second,vPar));
      for(std::map<RefinablePar*,CrystVector_REAL>::const_iterator d=pvV->begin();d!=pvV->end();d++)
      {
         if(mLSQ_FullDeriv[d->first].size()==0) mLSQ_FullDeriv[d->first].resize(mLSQObs.size());
         REAL *p2=mLSQ_FullDeriv[d->first].data()+nb;
         if(d->second.size()==0)
         {  //derivative can be null and then the vector missing
            // But we must still fill in zeros
            cout<<__FILE__<<":"<<__LINE__<<":"<<pos->first->GetClassName()<<":"<<pos->first->GetName()<<":"<<d->first->GetName()<<" (all deriv=0)"<<endl;
            for(unsigned long j=0;j<n2;++j) *p2++ = 0;
         }
         else
         {
            const REAL *p1=d->second.data();
            for(unsigned long j=0;j<n2;++j) *p2++ = *p1++;
         }
      }
      nb+=n2;
   }
   return mLSQ_FullDeriv;
}

void LSQNumObj::BeginOptimization(const bool allowApproximations, const bool enableRestraints)
{
   for(map<RefinableObj*,unsigned int>::iterator pos=mvRefinedObjMap.begin();pos!=mvRefinedObjMap.end();++pos)
      pos->first->BeginOptimization(allowApproximations, enableRestraints);
}

void LSQNumObj::EndOptimization()
{
   for(map<RefinableObj*,unsigned int>::iterator pos=mvRefinedObjMap.begin();pos!=mvRefinedObjMap.end();++pos)
      pos->first->EndOptimization();
}

#ifdef __WX__CRYST__
WXCrystObjBasic* LSQNumObj::WXCreate(wxWindow* parent)
{
   this->WXDelete();
   mpWXCrystObj=new WXLSQ(parent,this);
   return mpWXCrystObj;
}
WXCrystObjBasic* LSQNumObj::WXGet()
{
   return mpWXCrystObj;
}
void LSQNumObj::WXDelete()
{
   if(0!=mpWXCrystObj)
   {
      VFN_DEBUG_MESSAGE("LSQNumObj::WXDelete()",5)
      delete mpWXCrystObj;
   }
   mpWXCrystObj=0;
}
void LSQNumObj::WXNotifyDelete()
{
   VFN_DEBUG_MESSAGE("LSQNumObj::WXNotifyDelete():"<<mName,5)
   mpWXCrystObj=0;
}
#endif

}//namespace
