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
#include "RefinableObj/LSQNumObj.h"
#include "Quirks/VFNStreamFormat.h"

#include "newmat/newmatap.h" //for SVD decomposition
#include "newmat/newmatio.h"

using namespace NEWMAT;
using namespace std;

#include <iomanip>

namespace ObjCryst
{

LSQNumObj::LSQNumObj(string objName)
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
                        const bool silent)
{
   mpRefinedObj->BeginOptimization();
   mObs=mpRefinedObj->GetLSQObs(mLSQFuncIndex);
   mWeight=mpRefinedObj->GetLSQWeight(mLSQFuncIndex);

   //Check if we are ready for the refinement
   //:TODO:
   if(!silent) cout << "LSQNumObj::Refine():Beginning "<<endl;
   //Prepare for refinement (get non-fixed parameters)
      if(mRefParList.GetNbPar()==0) this->PrepareRefParList();
      mRefParList.PrepareForRefinement();
      if(!silent) mRefParList.Print();

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
      REAL R_ini,Rw_ini;
      REAL *pTmp1,*pTmp2;
      
      REAL marquardt=1e-2;
      const REAL marquardtMult=4.;
   //initial Chi^2, needed for Levenberg-Marquardt
   {
      calc=mpRefinedObj->GetLSQCalc(mLSQFuncIndex);
      tmpV1 = mObs;
      tmpV1 -= calc;
      tmpV1 *= tmpV1;
      tmpV1 *= mWeight;
      mChiSq=tmpV1.sum();
   }
   //store old values
   mIndexValuesSetInitial=mRefParList.CreateParamSet("LSQ Refinement-Initial Values");
   mIndexValuesSetLast=mRefParList.CreateParamSet("LSQ Refinement-Last Cycle Values");
   //refine
   for(int cycle=1 ; cycle <=nbCycle;cycle++)
   {
      mRefParList.SaveParamSet(mIndexValuesSetLast);// end of last cycle
      if(!silent) cout << "LSQNumObj::Refine():Cycle#"<< cycle <<endl;
      //initial value of function
         calc0=mpRefinedObj->GetLSQCalc(mLSQFuncIndex);
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
      designMatrix=0.;
      pTmp2=designMatrix.data();
      //cout <<"obs:"<<FormatHorizVector<REAL>(calc0,10,8);
      //cout <<"calc:"<<FormatHorizVector<REAL>(mObs,10,8);
      //cout <<"weight:"<<FormatHorizVector<REAL>(mWeight,10,8);
      for(i=0;i<nbVar;i++)
      {
         //:NOTE: Real design matrix is the transposed of the one computed here
         //if(!silent) cout << "........." << mRefParList.GetParNotFixed(i).GetName() <<endl;
         
         tmpV1=mpRefinedObj->GetLSQDeriv(mLSQFuncIndex,mRefParList.GetParNotFixed(i));
         pTmp1=tmpV1.data();
         //cout <<"deriv#"<<i<<":"<<FormatHorizVector<REAL>(tmpV1,10,8);
         for(j=0;j<nbObs;j++) *pTmp2++ = *pTmp1++;
      }
         //cout << designMatrix;
         
      LSQNumObj_Refine_Restart: //Used in case of singular matrix or for Marquardt
      
      //Calculate M and B matrices
         tmpV1.resize(nbObs);
         tmpV2.resize(nbObs);
         for(i=0;i<nbVar;i++) 
         {
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
         }
         
       //Apply LevenBerg-Marquardt factor
         if(true==useLevenbergMarquardt)
         {
            const REAL lmfact=1.+marquardt;
            for(i=0;i<nbVar;i++) M(i,i) *= lmfact;
         }
        
       // Check for singular values
         for(i=0;i<nbVar;i++)
         {
            if( (M(i,i) < 1e-20)||(ISNAN_OR_INF(M(i,i)))) //:TODO: Check what value to use as a limit
            {  
               if(!silent) cout << "LSQNumObj::Refine() Singular parameter !";
               if(!silent) cout << "(null derivate in all points) : "<<M(i,i)<<":";
               if(!silent) cout << mRefParList.GetParNotFixed(i).GetName() << endl;
               /*
               if(!silent)
               {
                  for(i=0;i<nbVar;i++)
                  {
                     tmpV1=mpRefinedObj->GetLSQDeriv(mLSQFuncIndex,mRefParList.GetParNotFixed(i));
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
                  mpRefinedObj->EndOptimization();
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
               goto LSQNumObj_Refine_Restart;
            }
         }
         
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
            EigenValues(newmatA,newmatW,newmatV);
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
         calc=mpRefinedObj->GetLSQCalc(mLSQFuncIndex);
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
                  marquardt *= marquardtMult;
                  if(!silent)
                  {
                     cout << "LSQNumObj::Refine(Chi^2="<<oldChiSq<<"->"<<mChiSq
                          <<")=>new Levenberg-Marquardt factor :"
                          << FormatFloat(marquardt,18,14) <<endl;
                  }
                  mChiSq=oldChiSq;
                  if(marquardt>1e8)
                  {
                     mRefParList.RestoreParamSet(mIndexValuesSetInitial);
                     mpRefinedObj->EndOptimization();
                     if(!silent) mRefParList.Print();
                     throw ObjCrystException("LSQNumObj::Refine():Levenberg-Marquardt diverging !");
                  }
                  goto LSQNumObj_Refine_Restart;
               }
               else marquardt /= marquardtMult;
               if(marquardt<1e-2) marquardt=1e-2;
               if(!silent) cout << "LSQNumObj::Refine():new Levenberg-Marquardt factor :" ;
               if(!silent) cout << FormatFloat(marquardt,18,14) <<endl;
            }
         }
            
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
         if(!silent) cout << "finished cycle #"<<cycle <<"/"<<nbCycle <<". Rw="<<Rw_ini<<"->"<<mRw<<endl;
         if (mSaveReportOnEachCycle) this->WriteReportToFile();
      
      if(!silent)this->PrintRefResults();
   }
   mpRefinedObj->EndOptimization();
}

CrystMatrix_REAL LSQNumObj::CorrelMatrix()const{return mCorrelMatrix;};

REAL LSQNumObj::Rfactor()const{return mR;};

REAL LSQNumObj::RwFactor()const{return mRw;};

REAL LSQNumObj::ChiSquare()const{return mChiSq;};

void LSQNumObj::SetRefinedObj(RefinableObj &obj, const unsigned int LSQFuncIndex)
{
   mpRefinedObj=&obj;
   mRefParList.ResetParList();
   mLSQFuncIndex=LSQFuncIndex;
   RefObjRegisterRecursive(obj,mRecursiveRefinedObjList);
}

ObjRegistry<RefinableObj> &LSQNumObj::GetRefinedObjList(){return mRecursiveRefinedObjList;}

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
   cout << "R-factor  : " << mR<<endl;
   cout << "Rw-factor : " << mRw<<endl;
   cout << "Chi-Square: " << mChiSq<<endl;
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
   for(int i=0;i<mRecursiveRefinedObjList.GetNb();i++)
   {
      VFN_DEBUG_MESSAGE("LSQNumObj::PrepareRefParList():"<<mRecursiveRefinedObjList.GetObj(i).GetName(),4);
      //mRecursiveRefinedObjList.GetObj(i).Print();
      mRefParList.AddPar(mRecursiveRefinedObjList.GetObj(i),copy_param);
   }
   //mRefParList.Print();
   if(copy_param) mRefParList.SetDeleteRefParInDestructor(true);
   else mRefParList.SetDeleteRefParInDestructor(false);
}

}//namespace
