/*
* ObjCryst++ : a Crystallographic computing library in C++
*
*  (c) 2000 Vincent FAVRE-NICOLIN
*           Laboratoire de Cristallographie
*           24, quai Ernest-Ansermet, CH-1211 Geneva 4, Switzerland
*  Contact: Vincent.Favre-Nicolin@cryst.unige.ch
*           Radovan.Cerny@cryst.unige.ch
*
*/
/*
*  source file LibCryst++ AsymmetricUnit and SpaceGroup classes
*
*/

#include <iomanip>
#include <cmath>
#include <typeinfo>

#include "ObjCryst/SpaceGroup.h"
#include "Quirks/VFNStreamFormat.h" //simple formatting of integers, REALs..
#include "ObjCryst/GeomStructFactor.h" //Geometrical Struct Factor definitions
#include "Quirks/VFNDebug.h"

#include <fstream>

namespace ObjCryst
{

////////////////////////////////////////////////////////////////////////
//
//   AsymmetricUnit
//
////////////////////////////////////////////////////////////////////////
AsymmetricUnit::AsymmetricUnit()
{
   VFN_DEBUG_MESSAGE("AsymmetricUnit::AsymmetricUnit()",5)
   mXmin=0;
   mYmin=0;
   mZmin=0;
   mXmax=1;
   mYmax=1;
   mZmax=1;
}

AsymmetricUnit::AsymmetricUnit(const SpaceGroup &spg)
{
   VFN_DEBUG_MESSAGE("AsymmetricUnit::AsymmetricUnit(SpGroup)",5)
   this->SetSpaceGroup(spg);
}

AsymmetricUnit::~AsymmetricUnit()
{
   VFN_DEBUG_MESSAGE("AsymmetricUnit::~AsymmetricUnit(SpGroup)",5)
}

void AsymmetricUnit::SetSpaceGroup(const SpaceGroup &spg)
{
   VFN_DEBUG_MESSAGE("AsymmetricUnit::SetSpaceGroup(SpGroup)",5)
   TAU_PROFILE("(AsymmetricUnit::SetSpaceGroup)","void (SpaceGroup)",TAU_DEFAULT);
   mXmin=0.;
   mYmin=0.;
   mZmin=0.;
   mXmax=1.;
   mYmax=1.;
   mZmax=1.;
   if(1==spg.GetSpaceGroupNumber()) return;//no need to search an asymmetric unit
   // Test points=reular grid of points inside the unit cell
   // All points must be or have at least a symmetric in the asymmetric unit
   const long nbPoints=13;
   CrystMatrix_REAL testPoints(nbPoints*nbPoints*nbPoints,3);
   {
      long l=0;
      for(long i=0;i<nbPoints;i++)
         for(long j=0;j<nbPoints;j++)
            for(long k=0;k<nbPoints;k++)
            {
               testPoints(l  ,0)=i/(REAL)nbPoints;
               testPoints(l  ,1)=j/(REAL)nbPoints;
               testPoints(l++,2)=k/(REAL)nbPoints;
            }
   }
   testPoints += 0.01;
   
   CrystVector_REAL vert(8);//vertices limits
   vert(0)=1/8.; vert(1)=1/6.; vert(2)=1/4.; vert(3)=1/3.;
   vert(4)=1/2.; vert(5)=2/3.; vert(6)=3/4.; vert(7)=1.;
   
   const int NbStep=vert.numElements();
                
   CrystMatrix_REAL coords;
   
   double junk;
   
   REAL minVolume=1.;
   
   bool allPtsInAsym,tmp;
   for(long nx=0;nx<NbStep;nx++)
      for(long ny=0;ny<NbStep;ny++)
         for(long nz=0;nz<NbStep;nz++)
         {
            if(minVolume<(vert(nx)*vert(ny)*vert(nz)-.0001)) break;
            allPtsInAsym=true;
            for(int i=0;i<testPoints.rows();i++)
            {
               coords=spg.GetAllSymmetrics(testPoints(i,0),testPoints(i,1),testPoints(i,2));
               for(long j=0;j<coords.numElements();j++) coords(j)=modf(coords(j)+10.,&junk) ;
               tmp=false;
               for(long j=0;j<coords.rows();j++)
               {//Test if at least one of the symmetrics is in the parallelepiped
                  if(  (coords(j,0) < vert(nx))
                     &&(coords(j,1) < vert(ny))
                     &&(coords(j,2) < vert(nz)))
                  {
                     //cout  << modf(coords(j,0)+10.,junk) << " "
                     //      << modf(coords(j,1)+10.,junk) << " "
                     //      << modf(coords(j,2)+10.,junk) << endl;
                     tmp=true;
                     break;
                  }
               }
               if(false==tmp)
               {
                  //cout << " Rejected:"<<vert(nx)<<" "<<vert(ny)<<" "<<vert(nz)<<" "<<i<<endl;
                  //cout << coords <<endl;
                  allPtsInAsym=false;
                  break;
               }
            }
            if( (true==allPtsInAsym))
            {
               mXmax=vert(nx);
               mYmax=vert(ny);
               mZmax=vert(nz);
               VFN_DEBUG_MESSAGE("->ACCEPTED:" << mXmax <<" "<< mYmax <<" "<< mZmax <<endl,2)
               //cout << "->ACCEPTED:" << mXmax <<" "<< mYmax <<" "<< mZmax <<endl;
               minVolume=vert(nx)*vert(ny)*vert(nz);
               break;//no need to grow any more along z
            }
         }
   VFN_DEBUG_MESSAGE("->Finished Generating Asymmetric Unit, with:",5)
   VFN_DEBUG_MESSAGE("     0 <= x <= "<< mXmax,10)
   VFN_DEBUG_MESSAGE("     0 <= y <= "<< mYmax,10)
   VFN_DEBUG_MESSAGE("     0 <= z <= "<< mZmax,10)
}

bool AsymmetricUnit::IsInAsymmetricUnit(const REAL x, const REAL y, const REAL z)const
{
   return (  ( x <= mXmin) && ( x >= mXmax)
           &&( y <= mYmin) && ( y >= mYmax)
           &&( z <= mZmin) && ( z >= mZmax));
}
REAL AsymmetricUnit::Xmin() const {return mXmin;}
REAL AsymmetricUnit::Xmax() const {return mXmax;}
REAL AsymmetricUnit::Ymin() const {return mYmin;}
REAL AsymmetricUnit::Ymax() const {return mYmax;}
REAL AsymmetricUnit::Zmin() const {return mZmin;}
REAL AsymmetricUnit::Zmax() const {return mZmax;}

////////////////////////////////////////////////////////////////////////
//
//    SpaceGroup 
//
////////////////////////////////////////////////////////////////////////

SpaceGroup::SpaceGroup():mId("P1")
{
   InitSpaceGroup(mId);
}

SpaceGroup::SpaceGroup(const string &spgId):mId(spgId)
{
   InitSpaceGroup(spgId);
}

SpaceGroup::~SpaceGroup()
{
}

void SpaceGroup::ChangeSpaceGroup(const string &spgId)
{
   VFN_DEBUG_MESSAGE("SpaceGroup::ChangeSpaceGroup():"<<spgId,5)
   this->InitSpaceGroup(spgId);
}

const string& SpaceGroup::GetName()const{return mId;}

bool SpaceGroup::IsInAsymmetricUnit(const REAL x, const REAL y, const REAL z) const
{
   return mAsymmetricUnit.IsInAsymmetricUnit(x,y,z);
}

void SpaceGroup::ChangeToAsymmetricUnit(REAL x, REAL y, REAL z) const
{
   //:TODO:
}

const AsymmetricUnit& SpaceGroup::GetAsymUnit() const {return mAsymmetricUnit;}


/// Id number of the spacegroup
int SpaceGroup::GetSpaceGroupNumber()const
{
   return mHM_as_Hall.SgNumber;
}
      
bool SpaceGroup::IsCentrosymmetric()const
{
   return mHasInversionCenter;
}

int SpaceGroup::GetNbTranslationVectors()const
{
   return mSgOps.nLTr;
}
      
CrystMatrix_REAL SpaceGroup::GetTranslationVectors()const
{
   const int *t1;
   REAL *t2;
   CrystMatrix_REAL transVect(this->GetNbTranslationVectors(),3);
   t1=&(mSgOps.LTr[0].v[0]);
   t2=transVect.data();
   for(int i=0;i<transVect.numElements();i++) *t2++ =  *t1++ / (REAL)STBF ;
   return transVect;
}


CrystMatrix_REAL SpaceGroup::GetAllSymmetrics(const REAL x, const REAL y, const REAL z,
                                const bool noCenter,const bool noTransl,
                                const bool noIdentical)const
{
   TAU_PROFILE("SpaceGroup::GetAllSymmetrics()","Matrix (x,y,z)",TAU_DEFAULT);
   VFN_DEBUG_MESSAGE("SpaceGroup::GetAllSymmetrics()",0)
   int nbMatrix, nbTrans,coeffInvert,i,j,k;
   nbMatrix=mSgOps.nSMx;
   nbTrans=this->GetNbTranslationVectors();
   if(this->IsCentrosymmetric()) coeffInvert=2 ; else coeffInvert=1;
   
   if(noCenter==true) coeffInvert=1;   //skip center of symmetry
   if(noTransl==true) nbTrans=1; //skip translation operations
   CrystMatrix_REAL coords(nbMatrix*nbTrans*coeffInvert,3);
   
   REAL tx,ty,tz;
   const int *t;
   t=&(mSgOps.LTr[0].v[0]);
   const T_RTMx *pMatrix;
   k=0;
   for(i=0;i<nbTrans;i++)
   {
      tx=*t++;
      ty=*t++;
      tz=*t++;
      //if(noTransl==false) cout << nbTrans <<endl;
      //if(noTransl==false) cout << tx <<" "<< ty<<" "<< tz<<" "<<endl;
      pMatrix=&(mSgOps.SMx[0]);
      for(j=0;j<nbMatrix;j++)
      {
         coords(k,0)= (*pMatrix).s.R[0]*x+(*pMatrix).s.R[1]*y+(*pMatrix).s.R[2]*z
                     +(*pMatrix).s.T[0]/(REAL)STBF+tx/(REAL)STBF;
         coords(k,1)= (*pMatrix).s.R[3]*x+(*pMatrix).s.R[4]*y+(*pMatrix).s.R[5]*z
                     +(*pMatrix).s.T[1]/(REAL)STBF+ty/(REAL)STBF;
         coords(k,2)= (*pMatrix).s.R[6]*x+(*pMatrix).s.R[7]*y+(*pMatrix).s.R[8]*z
                     +(*pMatrix).s.T[2]/(REAL)STBF+tz/(REAL)STBF;
         //cout<<(*pMatrix).s.R[0]<<" "<<(*pMatrix).s.R[1]<<" "<<(*pMatrix).s.R[2];
         //cout <<" "<<(*pMatrix).s.T[0]/(REAL)STBF<<endl;
         //cout<<(*pMatrix).s.R[3]<<" "<<(*pMatrix).s.R[4]<<" "<<(*pMatrix).s.R[5];
         //cout <<" "<<(*pMatrix).s.T[1]/(REAL)STBF<<endl;
         //cout<<(*pMatrix).s.R[5]<<" "<<(*pMatrix).s.R[6]<<" "<<(*pMatrix).s.R[7];
         //cout <<" "<<(*pMatrix).s.T[2]/(REAL)STBF<<endl<<endl;
         pMatrix++;
         k++;
      }
   }
   if(coeffInvert==2) //inversion center not in ListSeitzMx, but to be applied
   {
      int shift=nbMatrix*nbTrans;
      const REAL dx=((REAL)mSgOps.InvT[0])/STBF;//inversion not at the origin
      const REAL dy=((REAL)mSgOps.InvT[1])/STBF;
      const REAL dz=((REAL)mSgOps.InvT[2])/STBF;
      for(i=0;i<shift;i++)
      {
         coords(i+shift,0)=dx-coords(i,0);
         coords(i+shift,1)=dy-coords(i,1);
         coords(i+shift,2)=dz-coords(i,2);
      }
   }
   //for(i=0;i<nbTrans*nbMatrix*coeffInvert;i++) 
   //cout <<FormatFloat(coords(0,i))<<FormatFloat(coords(1,i))<<FormatFloat(coords(2,i))<<endl;
   //if(noTransl==false) cout <<coords<<endl;
   
   if(true==noIdentical)
   {
      VFN_DEBUG_MESSAGE("SpaceGroup::GetAllSymmetrics():Removing identical atoms",5)
      //Bring back all coordinates to [0;1[
      REAL *p=coords.data();
      double junk;
      for(long i=0;i<coords.numElements();i++)
      {
         *p = modf(*p,&junk);
         if(*p<0) *p += 1.;
         p++;
      }
      CrystMatrix_REAL newCoords;
      newCoords=coords;
      const REAL eps=1e-5;
      long nbKeep=0;
      for(long i=0;i<coords.rows();i++)
      {
         bool keep=true;
         for(long j=0;j<i;j++)
         {
            if(  ( fabs(coords(i,0)-coords(j,0)) < eps )
               &&( fabs(coords(i,1)-coords(j,1)) < eps )
               &&( fabs(coords(i,2)-coords(j,2)) < eps )) keep=false;
         }
         if(true==keep)
         {
            newCoords(nbKeep  ,0) = coords(i,0);
            newCoords(nbKeep  ,1) = coords(i,1);
            newCoords(nbKeep++,2) = coords(i,2);
         }
      }
      newCoords.resizeAndPreserve(nbKeep,3);
      return newCoords;
   }
   VFN_DEBUG_MESSAGE("SpaceGroup::GetAllSymmetrics():End",0)
   return coords;
}

int SpaceGroup::GetNbSymmetrics(const bool noCenter,const bool noTransl)const
{
   int nbMatrix, nbTrans,coeffInvert;
   nbMatrix=mSgOps.nSMx;
   nbTrans=this->GetNbTranslationVectors();
   if(this->IsCentrosymmetric()) coeffInvert=2 ; else coeffInvert=1;
   
   if(noCenter==true) coeffInvert=1;   //skip center of symmetry
   if(noTransl==true) nbTrans=1; //skip translation operations
   
   return nbMatrix*nbTrans*coeffInvert;
}

void SpaceGroup::Print() const
{
   cout << "SpaceGroup:" <<endl;
   cout << "  Schoenflies symbol = " << mHM_as_Hall.Schoenfl << endl ;
   cout << "  Hermann-Maugin symbol = " << mHM_as_Hall.HM << endl ;
   cout << "  Hall symbol = " << mHM_as_Hall.Hall << endl ;
   cout << "  SgNumber = " <<  mHM_as_Hall.SgNumber << endl ;
   cout << "  Number of Seitz Matrix = " <<  mSgOps.nSMx << endl ;
   cout << "  Number of Translation Vectors = " <<  mSgOps.nLTr << endl ;
   cout << "  List of Seitz Matrices : " << endl ;
   for(int i=0;i<mSgOps.nSMx;i++)
      cout << "    " << RTMx2XYZ(&mSgOps.SMx[i],1,STBF,0,0,1,NULL,NULL,80) <<endl;
   if(true==mHasInversionCenter)
   {
      cout << "  There is an inversion center at "
           << ((REAL)mSgOps.InvT[0])/STBF/2. << " "
           << ((REAL)mSgOps.InvT[1])/STBF/2. << " "
           << ((REAL)mSgOps.InvT[2])/STBF/2. << endl;
   }
   if(mSgOps.nLTr>0)
   {
      cout <<"  List of Translation vectors :"<<endl;
      for(int i=0;i<mSgOps.nLTr;i++)
         cout << "     "<< mSgOps.LTr[i].v[0]/(REAL)STBF<<","
                        << mSgOps.LTr[i].v[1]/(REAL)STBF<<","
                        << mSgOps.LTr[i].v[2]/(REAL)STBF<<endl;
   }
}
bool SpaceGroup::HasInversionCenter() const {return mHasInversionCenter;}
bool SpaceGroup::IsInversionCenterAtOrigin() const {return mIsInversionCenterAtOrigin;}
const T_SgOps& SpaceGroup::GetSgOps()const{return mSgOps;}

const RefinableObjClock& SpaceGroup::GetClockSpaceGroup() const{return mClock;}

const T_HM_as_Hall& SpaceGroup::GetHM_as_Hall()const
{
   return mHM_as_Hall;
}

unsigned int SpaceGroup::GetUniqueAxis()const{return mUniqueAxisId;}

void SpaceGroup::InitSpaceGroup(const string &spgId)
{
   VFN_DEBUG_MESSAGE("SpaceGroup::InitSpaceGroup():"<<spgId,5)
	(*fpObjCrystInformUser)("Initializing spacegroup: "+spgId);
#if 0
   mfpRealGeomStructFactor=0;
   mfpImagGeomStructFactor=0;
#endif
   int match;
   
   ResetSgOps(&mSgOps);
   match=SgSymbolLookup('A',spgId.c_str(),&mHM_as_Hall);
   if( match <= 0)
   {
      match=ParseHallSymbol(spgId.c_str(),&mSgOps,1);
      if( match <= 0)
      {
         cout << "An Error occured ! Cannot build SpaceGroup Info :" ;
         cout << "Cannot understand Spacegroup Symbol !" ;
         //:TODO: throw Exception
         this->InitSpaceGroup(mId);
			(*fpObjCrystInformUser)("Could not understand spacegroup symbol: "+spgId);
         return;
         //throw 0;
      }
      MatchTabulatedSettings(&mSgOps,&mHM_as_Hall);
   }
   else ParseHallSymbol(mHM_as_Hall.Hall,&mSgOps,1);
   mId=spgId;
      
   //Inversion center
   if(mSgOps.fInv == 2)
   {
      mHasInversionCenter=true ;
      if( (mSgOps.InvT[0] !=0) || 
          (mSgOps.InvT[1] !=0) || 
          (mSgOps.InvT[2] !=0)   ) mIsInversionCenterAtOrigin=false;
      else mIsInversionCenterAtOrigin=true;
   }
   else
   {
      mHasInversionCenter=false ;
      mIsInversionCenterAtOrigin=true;
   }   
   
   //initialize asymmetric unit
   mAsymmetricUnit.SetSpaceGroup(*this);
   
   //:TODO: This may be obsolete... Not using Geometrical structure factors any more ?
#if 0
   const char *ext=mHM_as_Hall.Qualif ;
   switch(mHM_as_Hall.SgNumber)
   {
      case 1:  mfpRealGeomStructFactor=RealGeomStructFactor_1        ;break;
      case 2:  mfpRealGeomStructFactor=RealGeomStructFactor_2        ;break;
      case 67: 
         mfpRealGeomStructFactor=RealGeomStructFactor_67;
         if(strcmp(ext,"ba-c")==0)mfpRealGeomStructFactor=RealGeomStructFactor_67ba_c;
         if(strcmp(ext,"cab" )==0)mfpRealGeomStructFactor=RealGeomStructFactor_67cab;
         if(strcmp(ext,"-cba")==0)mfpRealGeomStructFactor=RealGeomStructFactor_67_cba;
         if(strcmp(ext,"bca" )==0)mfpRealGeomStructFactor=RealGeomStructFactor_67bca;
         if(strcmp(ext,"a-cb")==0)mfpRealGeomStructFactor=RealGeomStructFactor_67a_cb;
         break;
      case 97: mfpRealGeomStructFactor=RealGeomStructFactor_97       ;break;
      case 230:mfpRealGeomStructFactor=RealGeomStructFactor_230      ;break;

      //should return an errror there...
      default: break;
   }
   
   switch(mHM_as_Hall.SgNumber)
   {
      case 1:  mfpImagGeomStructFactor=ImagGeomStructFactor_1        ;break;
      case 2:  mfpImagGeomStructFactor=ImagGeomStructFactor_centro   ;break;
      case 67: mfpImagGeomStructFactor=ImagGeomStructFactor_centro   ;break;
      case 97: mfpImagGeomStructFactor=ImagGeomStructFactor_97       ;break;
      case 230:mfpImagGeomStructFactor=ImagGeomStructFactor_centro   ;break;
      
      //should throw an exception here...
      default: break ;
   }
#endif
   if( (mHM_as_Hall.SgNumber >2) && (mHM_as_Hall.SgNumber <16))
   {
      const char * ch=this->GetHM_as_Hall().Hall;
      while(true)
      {
         if(*ch=='x') {mUniqueAxisId=0;break;}
         if(*ch=='y') {mUniqueAxisId=1;break;}
         if(*ch=='z') {mUniqueAxisId=2;break;}
         if(*ch=='\0'){mUniqueAxisId=2;break;}//:TODO: check ??
         ch++;
      }
   }else mUniqueAxisId=0;
   //this->Print();
   mClock.Click();
	(*fpObjCrystInformUser)("Initializing spacegroup: "+spgId+"... Done");
   VFN_DEBUG_MESSAGE("SpaceGroup::InitSpaceGroup():End",4)
}

}//namespace
