#include "RefinableObj/LSQNumObj.h"
#include "Quirks/VFNStreamFormat.h"

#include "newmat/newmatap.h" //for SVD decomposition
#include "newmat/newmatio.h"

using namespace NEWMAT;

#include <iomanip>

LSQNumObj::LSQNumObj(CrystVector_double (*pFunc)(),
                     string objName)
{
	mpRefinedFunc=pFunc;
	mDampingFactor=1.;
	mSaveReportOnEachCycle=false;
	mName=objName;
	mSaveFileName="LSQrefinement.save";
   mR=0;
   mRw=0;
   mChiSq=0;
   mStopAfterCycle=false;
	//We only use copies of parameters.
	mRefParList.SetDeleteRefParInDestructor(false);
}

LSQNumObj::~LSQNumObj()
{
}

void LSQNumObj::Init(const CrystVector_double &obs,const CrystVector_double&w)
{
	mObs=obs;
	mWeight=w;
}

void LSQNumObj::Init(const CrystVector_double &obs)
{
	CrystVector_double tmp(obs.numElements());
	tmp=1.;
	this->Init(obs,tmp);
}

void LSQNumObj::SetParIsFixed(const string& parName,const bool fix)
{
   for(int i=0;i<mRecursiveRefinedObjList.GetNb();i++) 
      mRecursiveRefinedObjList.GetObj(i).SetParIsFixed(parName,fix);
}
void LSQNumObj::SetParIsFixed(const RefParType *type,const bool fix)
{
   for(int i=0;i<mRecursiveRefinedObjList.GetNb();i++) 
      mRecursiveRefinedObjList.GetObj(i).SetParIsFixed(type,fix);
}
   
void LSQNumObj::UnFixAllPar()
{
   for(int i=0;i<mRecursiveRefinedObjList.GetNb();i++) 
      mRecursiveRefinedObjList.GetObj(i).UnFixAllPar();
}
   
void LSQNumObj::SetParIsUsed(const string& parName,const bool use)
{
   for(int i=0;i<mRecursiveRefinedObjList.GetNb();i++) 
      mRecursiveRefinedObjList.GetObj(i).SetParIsUsed(parName,use);
}
void LSQNumObj::SetParIsUsed(const RefParType *type,const bool use)
{
   for(int i=0;i<mRecursiveRefinedObjList.GetNb();i++) 
      mRecursiveRefinedObjList.GetObj(i).SetParIsUsed(type,use);
}

void LSQNumObj::Refine (int nbCycle,bool useLevenbergMarquardt)
{

	//Check if we are ready for the refinement
	//:TODO:
   cout << "LSQNumObj::Refine():Beginning "<<endl;
   //Prepare for refinement (get non-fixed parameters)
      this->PrepareRefParList();
      mRefParList.PrepareForRefinement();
      mRefParList.Print();

	//variables
      long nbVar=mRefParList.GetNbParNotFixed();
		const long nbObs=mObs.numElements();
		CrystVector_double calc,calc0,calc1,tmpV1,tmpV2;
		CrystMatrix_double M(nbVar,nbVar);
		CrystMatrix_double N(nbVar,nbVar);
		CrystVector_double B(nbVar);
		CrystMatrix_double designMatrix(nbVar,nbObs);
		CrystVector_double deltaVar(nbVar);
		long i,j,k;
		double R_ini,Rw_ini;
      double *pTmp1,*pTmp2;
      
      double marquardt=1e-3;
      const double marquardtMult=10.;
	//store old values
   mIndexValuesSetInitial=mRefParList.CreateParamSet("LSQ Refinement-Initial Values");
   mIndexValuesSetLast=mRefParList.CreateParamSet("LSQ Refinement-Last Cycle Values");
	//refine
	for(int cycle=1 ; cycle <=nbCycle;cycle++)
	{
		mRefParList.SaveParamSet(mIndexValuesSetLast);// end of last cycle
      cout << "LSQNumObj::Refine():Cycle#"<< cycle <<endl;
      cout << "LSQNumObj::Refine():Computing initial values" <<endl;
		//initial value of function
			calc0=(*mpRefinedFunc)();
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
      cout << "LSQNumObj::Refine():Computing derivatives" <<endl;
      designMatrix=0.;
      pTmp2=designMatrix.data();
		for(i=0;i<nbVar;i++)
		{
         //:NOTE: Real design matrix is the transposed of the one computed here
         
         cout << "........." << mRefParList.GetParNotFixed(i).GetName() <<endl;
         //cout << mRefParList.GetParNotFixed(i).DerivStep() <<endl;
         
         /*
         //Forward derivative
			mRefParList.GetParNotFixed(i).Value() += mRefParList.GetParNotFixed(i).DerivStep();
			tmpV1=(*mpRefinedFunc)();
			mRefParList.GetParNotFixed(i).Value() -= mRefParList.GetParNotFixed(i).DerivStep();
         tmpV1 -= calc0;
         tmpV1 /= mRefParList.GetParNotFixed(i).DerivStep();
         */
         
         //centered derivative
			mRefParList.GetParNotFixed(i).Mutate(mRefParList.GetParNotFixed(i).GetDerivStep());
			calc=(*mpRefinedFunc)();
			mRefParList.GetParNotFixed(i).Mutate(-2*mRefParList.GetParNotFixed(i).GetDerivStep());
			calc1=(*mpRefinedFunc)();
			mRefParList.GetParNotFixed(i).Mutate(mRefParList.GetParNotFixed(i).GetDerivStep());
         tmpV1 =  calc;
         tmpV1 -= calc1;
         tmpV1 /= mRefParList.GetParNotFixed(i).GetDerivStep()/2;
         
         
         pTmp1=tmpV1.data();
         for(j=0;j<nbObs;j++) *pTmp2++ = *pTmp1++;
		}
         //cout << designMatrix;
         
      LSQNumObj_Refine_Restart: //Used in case of singular matrix or for Marquardt
      
      cout << "LSQNumObj::Refine():Computing M and B Matrices" <<endl;
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
            const double lmfact=1.+marquardt;
			   for(i=0;i<nbVar;i++) M(i,i) *= lmfact;
         }
        
       // Check for singular values
         for(i=0;i<nbVar;i++)
         {
            if( M(i,i) < 1e-10) //:TODO: Check what value to use as a limit
            {  
               cout << "LSQNumObj::Refine()Singular parameter";
               cout << "(null derivate in all points) : ";
               cout << mRefParList.GetParNotFixed(i).GetName() << endl;
               cout << "LSQNumObj::Refine(): Automatically fixing parameter";
               cout << " and re-start cycle..";
               cout << endl;
               mRefParList.GetParNotFixed(i).SetIsFixed(true);
               mRefParList.PrepareForRefinement();
               nbVar=mRefParList.GetNbParNotFixed();
		         N.resize(nbVar,nbVar);
		         deltaVar.resize(nbVar);

               //Just remove the ith line in the design matrix
                        double *p1=designMatrix.data();
                  const double *p2=designMatrix.data();
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
         
         CrystMatrix_double U(nbVar,nbObs), V(nbVar,nbVar);
         CrystVector_double W(nbVar),invW(nbVar);
         for(long i=0;i<nbVar;i++)
         {
            for(long j=0;j<nbObs;j++)
               newmatA(j+1,i+1) = designMatrix(i,j);//design matrix (need to transpose)
         }
         
         for(long j=0;j<nbObs;j++) newmatB(j+1)=mObs(j)-calc0(j);
         
         SVD(newmatA,newmatW,newmatU,newmatV);
         
         ColumnVector newmatDelta(nbVar);
         DiagonalMatrix newmatInvW(nbVar);
         double max=newmatW.MaximumAbsoluteValue();
         double minAllowedValue=1e-16*max;// :TODO: Check if reasonable !
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
         cout << "LSQNumObj::Refine():Eigenvalue Filtering..." <<endl;
         CrystMatrix_double V(nbVar,nbVar);
         CrystVector_double invW(nbVar);//diagonal matrix, in fact
         cout << "LSQNumObj::Refine():Eigenvalue Filtering...1" <<endl;
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
         cout << "LSQNumObj::Refine():Eigenvalue Filtering...2" <<endl;
            
            //cout << newmatA.SubMatrix(1,3,1,3) <<endl;
         cout << "LSQNumObj::Refine():Eigenvalue Filtering...3" <<endl;
            
            //Jacobi(newmatA,newmatW,newmatV);
            EigenValues(newmatA,newmatW,newmatV);
            ColumnVector newmatDelta(nbVar);
            DiagonalMatrix newmatInvW(nbVar);
         cout << "LSQNumObj::Refine():Eigenvalue Filtering...4" <<endl;
            //Avoid singular values
            {
               double max=newmatW.MaximumAbsoluteValue();
               double minAllowedValue=1e-5*max;// :TODO: Check if reasonable !
               for(long i=0;i<nbVar;i++) 
                  if(newmatW(i+1,i+1) > minAllowedValue) 
                     newmatInvW(i+1,i+1)= 1./newmatW(i+1,i+1);
                  else 
                  {
                     cout << "LSQNumObj::Refine():fixing ill-cond EigenValue "<< i <<endl;
                     newmatInvW(i+1,i+1) = 0.;
                  }
               
               /*
               double max;
               int maxIndex;
               const double minRatio=1e-6;
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
               cout << "Back-scaling" << endl;
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
         CrystMatrix_double U(nbVar,nbVar),V(nbVar,nbVar);
         CrystVector_double invW(nbVar);//diagonal matrix, in fact
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
            double max=newmatW.MaximumAbsoluteValue();
            double minAllowedValue=1e-10*max;// :TODO: Check if reasonable !
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
               CrystMatrix_double Mtmp;
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
            CrystMatrix_double Mtmp,Ntmp;
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
         cout << "LSQNumObj::Refine():Computing new values for variables" <<endl;
         

			for(i=0;i<nbVar;i++) mRefParList.GetParNotFixed(i).Mutate(deltaVar(i));
         
      cout << "LSQNumObj::Refine():Computing statistics for last cycle..." <<endl;
		//for statistics...
         //mRefParList.Print();
			calc=(*mpRefinedFunc)();
         //Chi^2
      cout << "LSQNumObj::Refine():Computing Chi^2 for last cycle..." <<endl;
         {
            double oldChiSq=mChiSq;
            tmpV1 = mObs;
            tmpV1 -= calc;
            tmpV1 *= mWeight;
            tmpV1 *= tmpV1;
			   mChiSq=tmpV1.sum()/(nbObs-nbVar);
            if(true==useLevenbergMarquardt)
            {
               if(mChiSq > oldChiSq)
               {
                  marquardt *= marquardtMult;
                  cout << "LSQNumObj::Refine():new Levenberg-Marquardt factor :" ;
                  cout << FormatFloat(marquardt,18,14) <<endl;
                  mChiSq=oldChiSq;
                  goto LSQNumObj_Refine_Restart;
               }
               else marquardt /= marquardtMult;
               cout << "LSQNumObj::Refine():new Levenberg-Marquardt factor :" ;
               cout << FormatFloat(marquardt,18,14) <<endl;
            }
         }
            
         //Sigmas
			   for(i=0;i<nbVar;i++) mRefParList.GetParNotFixed(i).SetSigma(sqrt(N(i,i)*mChiSq));
         //Correlations :TODO: re-compute M and N if using Levenberg-Marquardt
         mCorrelMatrix.resize(nbVar,nbVar);
			for(i=0;i<nbVar;i++)
            for(j=0;j<nbVar;j++) mCorrelMatrix(i,j)=sqrt(N(i,j)*N(i,j)/N(i,i)/N(j,j));
         //R-factor
            tmpV1 = mObs;
            tmpV1 -= calc;
            tmpV1 *= tmpV1;
            tmpV2 = mObs;
            tmpV2 *= tmpV2;
			   mR=sqrt(tmpV1.sum()/tmpV2.sum());
         //Rw-factor
            tmpV1 = mObs;
            tmpV1 -= calc;
            tmpV1 *= mWeight;
            tmpV1 *= tmpV1;
            tmpV2 = mObs;
            tmpV2 *= mWeight;
            tmpV2 *= tmpV2;
			   mRw=sqrt(tmpV1.sum()/tmpV2.sum());
		//OK, finished
			cout << "finished cycle #"<<cycle <<"/"<<nbCycle <<". R="<<R_ini<<"->"<<mR<<endl;
			if (mSaveReportOnEachCycle) this->WriteReportToFile();
      
      this->PrintRefResults();
	}
}

CrystMatrix_double LSQNumObj::CorrelMatrix()const{return mCorrelMatrix;};

double LSQNumObj::Rfactor()const{return mR;};

double LSQNumObj::RwFactor()const{return mRw;};

double LSQNumObj::ChiSquare()const{return mChiSq;};

void LSQNumObj::AddRefinableObj(RefinableObj &obj)
{
   mRefinedObjList.Register(obj);
   RefObjRegisterRecursive(obj,mRecursiveRefinedObjList);
}

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
      cout << FormatFloat((mRefParList.GetParamSet(mIndexValuesSetInitial))(i),15,8) << "  " ;
      cout << FormatFloat((mRefParList.GetParamSet(mIndexValuesSetLast))(i),15,8) << "  " ;
      cout << FormatFloat(mRefParList.GetPar(i).GetHumanValue(),15,8) << "  " ;
      cout << FormatFloat(mRefParList.GetPar(i).GetHumanSigma(),15,8) << "  " ;
      //cout << FormatFloat(mRefParList(i).DerivStep(),16,12) << "  ";
		//cout << varNames[i] << "  " << var0(i) << "  " << varLast(i) 
      //cout<< "  " << varCurrent(i)<< "  " << sigmaValues(i)<<endl;
      cout << endl;
	}
	cout <<endl;
}

void LSQNumObj::SetDampingFactor(const double newDampFact)
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

void LSQNumObj::PrepareRefParList()
{
   //:TODO: instead of resetting the list every time, check if it is necessary !
   mRefParList.ResetParList();
   for(int i=0;i<mRecursiveRefinedObjList.GetNb();i++)
   {
      mRecursiveRefinedObjList.GetObj(i).Print();
      mRefParList.AddPar( mRecursiveRefinedObjList.GetObj(i));
   }
   mRefParList.Print();
}


