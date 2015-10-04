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
/*
*  source file LibCryst++ AsymmetricUnit and SpaceGroup classes
*
*/

#include <iomanip>
#include <cmath>
#include <typeinfo>

#include "cctbx/sgtbx/space_group.h"
#include "cctbx/sgtbx/brick.h"
#include "cctbx/miller/sym_equiv.h"
#include "boost/rational.hpp"

#include "ObjCryst/ObjCryst/SpaceGroup.h"
#include "ObjCryst/Quirks/VFNStreamFormat.h" //simple formatting of integers, REALs..
#include "ObjCryst/ObjCryst/GeomStructFactor.h" //Geometrical Struct Factor definitions
#include "ObjCryst/Quirks/VFNDebug.h"


#include <fstream>

namespace ObjCryst
{

#include "ObjCryst/Quirks/VFNDebug.h"

// We need to force the C locale when using cctbx (when interpreting xyz strings)
tmp_C_Numeric_locale::tmp_C_Numeric_locale()
{
   char *old;
   old=setlocale(LC_NUMERIC,NULL);
   mLocale=old;
   setlocale(LC_NUMERIC,"C");
}

tmp_C_Numeric_locale::~tmp_C_Numeric_locale()
{
   setlocale(LC_NUMERIC,mLocale.c_str());
}

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
   tmp_C_Numeric_locale tmploc;
   # if 0
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
   cout<<"->Finished Generating (pseudo) Asymmetric Unit, with:"<<endl
       <<"     0 <= x <= "<< mXmax<<endl
       <<"     0 <= y <= "<< mYmax<<endl
       <<"     0 <= z <= "<< mZmax<<endl<<endl;
   #else
   const cctbx::sgtbx::brick b(spg.GetCCTbxSpg().type());
   #ifdef __DEBUG__
   cout<<"->>Parallelepipedic Asymmetric Unit, from cctbx::sgtbx::brick:"<<endl
       <<b.as_string()<<endl;
   #endif
   mXmin=boost::rational_cast<REAL,int>(b(0,0).value());
   mYmin=boost::rational_cast<REAL,int>(b(1,0).value());
   mZmin=boost::rational_cast<REAL,int>(b(2,0).value());
   mXmax=boost::rational_cast<REAL,int>(b(0,1).value());
   mYmax=boost::rational_cast<REAL,int>(b(1,1).value());
   mZmax=boost::rational_cast<REAL,int>(b(2,1).value());

   #endif
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

SpaceGroup::SpaceGroup():mId("P1"),mpCCTbxSpaceGroup(0)
{
   InitSpaceGroup(mId);
}

SpaceGroup::SpaceGroup(const string &spgId):mId(spgId),mpCCTbxSpaceGroup(0)
{
   InitSpaceGroup(spgId);
}

SpaceGroup::~SpaceGroup()
{
   if(mpCCTbxSpaceGroup!=0) delete mpCCTbxSpaceGroup;
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
   return mSpgNumber;
}

bool SpaceGroup::IsCentrosymmetric()const
{
   return mHasInversionCenter;
}

int SpaceGroup::GetNbTranslationVectors()const
{
   return mNbTrans;
}

const std::vector<SpaceGroup::TRx>& SpaceGroup::GetTranslationVectors()const
{
   return mvTrans;
}

const std::vector<SpaceGroup::SMx>& SpaceGroup::GetSymmetryOperations()const
{
   return mvSym;
}

CrystMatrix_REAL SpaceGroup::GetAllSymmetrics(const REAL x, const REAL y, const REAL z,
                                const bool noCenter,const bool noTransl,
                                const bool noIdentical)const
{
   //TAU_PROFILE("SpaceGroup::GetAllSymmetrics()","Matrix (x,y,z)",TAU_DEFAULT);
   VFN_DEBUG_MESSAGE("SpaceGroup::GetAllSymmetrics()",0)
   int nbMatrix, nbTrans,coeffInvert,i,j,k;
   nbMatrix=this->GetCCTbxSpg().n_smx();
   nbTrans=this->GetNbTranslationVectors();
   if(this->IsCentrosymmetric()) coeffInvert=2 ; else coeffInvert=1;

   if(noCenter==true) coeffInvert=1;   //skip center of symmetry
   if(noTransl==true) nbTrans=1; //skip translation operations
   CrystMatrix_REAL coords(nbMatrix*nbTrans*coeffInvert,3);

   k=0;
   for(i=0;i<nbTrans;i++)
   {
      const REAL ltr_den=1/(REAL)(this->GetCCTbxSpg().ltr(i).den());
      const REAL tx=this->GetCCTbxSpg().ltr(i)[0]*ltr_den;
      const REAL ty=this->GetCCTbxSpg().ltr(i)[1]*ltr_den;
      const REAL tz=this->GetCCTbxSpg().ltr(i)[2]*ltr_den;
      //if(noTransl==false) cout << nbTrans <<endl;
      //if(noTransl==false) cout << tx <<" "<< ty<<" "<< tz<<" "<<endl;
      for(j=0;j<nbMatrix;j++)
      {
         const cctbx::sgtbx::rt_mx *pMatrix=&(this->GetCCTbxSpg().smx(j));
         const cctbx::sgtbx::rot_mx *pRot=&(pMatrix->r());
         const cctbx::sgtbx::tr_vec *pTrans=&(pMatrix->t());
         const REAL r_den=1/(REAL)(pMatrix->r().den());
         const REAL t_den=1/(REAL)(pMatrix->t().den());
         coords(k,0)= ((*pRot)[0]*x+(*pRot)[1]*y+(*pRot)[2]*z)*r_den+(*pTrans)[0]*t_den+tx;
         coords(k,1)= ((*pRot)[3]*x+(*pRot)[4]*y+(*pRot)[5]*z)*r_den+(*pTrans)[1]*t_den+ty;
         coords(k,2)= ((*pRot)[6]*x+(*pRot)[7]*y+(*pRot)[8]*z)*r_den+(*pTrans)[2]*t_den+tz;
         k++;
      }
   }
   if(coeffInvert==2) //inversion center not in ListSeitzMx, but to be applied
   {
      int shift=nbMatrix*nbTrans;
      const REAL dx=((REAL)this->GetCCTbxSpg().inv_t()[0])/(REAL)this->GetCCTbxSpg().inv_t().den();//inversion not at the origin
      const REAL dy=((REAL)this->GetCCTbxSpg().inv_t()[1])/(REAL)this->GetCCTbxSpg().inv_t().den();
      const REAL dz=((REAL)this->GetCCTbxSpg().inv_t()[2])/(REAL)this->GetCCTbxSpg().inv_t().den();
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
void SpaceGroup::GetSymmetric(unsigned int idx, REAL &x, REAL &y, REAL &z,
                              const bool noCenter,const bool noTransl,
                              const bool derivative) const
{
   int coeffInvert;
   const int nbMatrix=this->GetCCTbxSpg().n_smx();
   int nbTrans=this->GetNbTranslationVectors();
   if(this->IsCentrosymmetric()) coeffInvert=2 ; else coeffInvert=1;

   if(noCenter==true) coeffInvert=1;   //skip center of symmetry
   if(noTransl==true) nbTrans=1; //skip translation operations

   unsigned int idx0=idx;
   const unsigned int mxidx = nbTrans * nbMatrix;
   if(idx > mxidx) idx0 = idx % mxidx;
   const int i=idx/nbMatrix;//translation index
   const int j=idx%nbMatrix;

   const REAL ltr_den=1/(REAL)(this->GetCCTbxSpg().ltr(i).den());
   const REAL tx=this->GetCCTbxSpg().ltr(i)[0]*ltr_den;
   const REAL ty=this->GetCCTbxSpg().ltr(i)[1]*ltr_den;
   const REAL tz=this->GetCCTbxSpg().ltr(i)[2]*ltr_den;
   const cctbx::sgtbx::rt_mx *pMatrix=&(this->GetCCTbxSpg().smx(j));
   const cctbx::sgtbx::rot_mx *pRot=&(pMatrix->r());
   const cctbx::sgtbx::tr_vec *pTrans=&(pMatrix->t());
   const REAL r_den=1/(REAL)(pMatrix->r().den());
   const REAL t_den=1/(REAL)(pMatrix->t().den());
   const REAL x1= ((*pRot)[0]*x+(*pRot)[1]*y+(*pRot)[2]*z)*r_den;
   const REAL y1= ((*pRot)[3]*x+(*pRot)[4]*y+(*pRot)[5]*z)*r_den;
   const REAL z1= ((*pRot)[6]*x+(*pRot)[7]*y+(*pRot)[8]*z)*r_den;
   if(derivative==false)
   {
      x=x1+(*pTrans)[0]*t_den+tx;
      y=y1+(*pTrans)[1]*t_den+ty;
      z=z1+(*pTrans)[2]*t_den+tz;
   }
   else
   {
      x=x1;
      y=y1;
      z=z1;
   }
   if(coeffInvert==2) //inversion center not in ListSeitzMx, but to be applied
   {
      const REAL dx=((REAL)this->GetCCTbxSpg().inv_t()[0])/(REAL)this->GetCCTbxSpg().inv_t().den();//inversion not at the origin
      const REAL dy=((REAL)this->GetCCTbxSpg().inv_t()[1])/(REAL)this->GetCCTbxSpg().inv_t().den();
      const REAL dz=((REAL)this->GetCCTbxSpg().inv_t()[2])/(REAL)this->GetCCTbxSpg().inv_t().den();
      x=dx-x;
      y=dy-y;
      z=dz-z;
   }
}

int SpaceGroup::GetNbSymmetrics(const bool noCenter,const bool noTransl)const
{
   if(noCenter || (!mHasInversionCenter))
   {
      if(noTransl) return mNbSym;
      else return mNbSym*mNbTrans;
   }
   else
   {
      if(noTransl) return 2*mNbSym;
      else return 2*mNbSym*mNbTrans;
   }
   return 2*mNbSym*mNbTrans;
}

void SpaceGroup::Print() const
{
   cout << "SpaceGroup:" <<endl;
   cout << "  Schoenflies symbol = " << this->GetCCTbxSpg().match_tabulated_settings().schoenflies() << endl ;
   cout << "  Hermann-Maugin symbol = " <<  this->GetCCTbxSpg().match_tabulated_settings().hermann_mauguin() << endl ;
   cout << "  Hall symbol = " <<  this->GetCCTbxSpg().match_tabulated_settings().hall() << endl ;
   cout << "  SgNumber = " <<   this->GetCCTbxSpg().match_tabulated_settings().number() << endl ;
   cout << "  Number of Seitz Matrix = " <<  this->GetCCTbxSpg().n_smx() << endl ;
   cout << "  Number of Translation Vectors = " <<   this->GetCCTbxSpg().n_ltr() << endl ;
   cout << "  List of Seitz Matrices : " << endl ;
   for(unsigned int i=0;i<this->GetCCTbxSpg().n_smx();i++)
      cout << "    " <<  this->GetCCTbxSpg().smx(i).as_xyz() <<endl;
   if(true==mHasInversionCenter)
   {
      cout << "  There is an inversion center at "
           << (this->GetCCTbxSpg().inv_t()[0])/(REAL)this->GetCCTbxSpg().inv_t().den()/2. << " "
           << (this->GetCCTbxSpg().inv_t()[1])/(REAL)this->GetCCTbxSpg().inv_t().den()/2. << " "
           << (this->GetCCTbxSpg().inv_t()[2])/(REAL)this->GetCCTbxSpg().inv_t().den()/2. << endl;
   }
   if(this->GetCCTbxSpg().n_ltr()>0)
   {
      cout <<"  List of Translation vectors :"<<endl;
      for(unsigned int i=0;i<this->GetCCTbxSpg().n_ltr();i++)
         cout << "     "<< this->GetCCTbxSpg().ltr(i)[0]/(REAL)this->GetCCTbxSpg().ltr(i).den()<<","
                        << this->GetCCTbxSpg().ltr(i)[1]/(REAL)this->GetCCTbxSpg().ltr(i).den()<<","
                        << this->GetCCTbxSpg().ltr(i)[2]/(REAL)this->GetCCTbxSpg().ltr(i).den()<<endl;
   }
   cout<<"Extension (origin choice, rhomboedral/hexagonal):"<<mExtension<<endl;
}
bool SpaceGroup::HasInversionCenter() const {return mHasInversionCenter;}
bool SpaceGroup::IsInversionCenterAtOrigin() const {return mIsInversionCenterAtOrigin;}
const cctbx::sgtbx::space_group& SpaceGroup::GetCCTbxSpg()const{return *mpCCTbxSpaceGroup;}

const RefinableObjClock& SpaceGroup::GetClockSpaceGroup() const{return mClock;}

unsigned int SpaceGroup::GetUniqueAxis()const{return mUniqueAxisId;}

char SpaceGroup::GetExtension()const{return mExtension;}

CrystVector_REAL SpaceGroup::GetInversionCenter()const {
   CrystVector_REAL center(3);
   center(0) =((REAL)this->GetCCTbxSpg().inv_t()[0])/(REAL)this->GetCCTbxSpg().inv_t().den();//inversion not at the origin
   center(1) =((REAL)this->GetCCTbxSpg().inv_t()[1])/(REAL)this->GetCCTbxSpg().inv_t().den();
   center(2) =((REAL)this->GetCCTbxSpg().inv_t()[2])/(REAL)this->GetCCTbxSpg().inv_t().den();
   return center;
}

unsigned int SpaceGroup::AreReflEquiv(const REAL h1, const REAL k1, const REAL l1,
                                      const REAL h2, const REAL k2, const REAL l2)const
{
   const int ih1=scitbx::math::iround(h1);
   const int ik1=scitbx::math::iround(k1);
   const int il1=scitbx::math::iround(l1);
   cctbx::miller::index<long> h0(scitbx::vec3<long>(ih1,ik1,il1));
   const int ih2=scitbx::math::iround(h2);
   const int ik2=scitbx::math::iround(k2);
   const int il2=scitbx::math::iround(l2);
   cctbx::miller::index<long> k0(scitbx::vec3<long>(ih2,ik2,il2));
   cctbx::miller::sym_equiv_indices sei(this->GetCCTbxSpg(),k0);
   int equiv=0;
   //cout<<h0.as_string()<<" - "<<k0.as_string()<<","<<sei.f_mates(false)<<","<<sei.f_mates(true)<<endl;
   for(size_t i_indices=0;i_indices<sei.indices().size();i_indices++)
   {
      const size_t sfm = sei.f_mates(false);
      for(size_t i_mate = 0; i_mate < sfm; i_mate++)
      {
         cctbx::miller::index<long> k = sei(i_mate, i_indices).h();
         //cout<<" ->("<<i_indices<<","<<i_mate<<")"<<k.as_string()<<endl;
         if(h0==k)
         {
            if(i_mate==0) equiv=1;
            else equiv=2;
            break;
         }
      }
   }
   VFN_DEBUG_MESSAGE("SpaceGroup::AreReflEquiv("<<ih1<<","<<ik1<<","<<il1<<"),("<<ih2<<","<<ik2<<","<<il2<<"):"<<equiv,2)
   return equiv;
}

CrystMatrix_REAL SpaceGroup::GetAllEquivRefl(const REAL h0, const REAL k0, const REAL l0,
                                             const bool excludeFriedelMate,
                                             const bool forceFriedelLaw,
                                             const REAL sf_re,const REAL sf_im) const
{
   VFN_DEBUG_ENTRY("SpaceGroup::GetAllEquivRefl():",5)
   const int ih0=scitbx::math::iround(h0);
   const int ik0=scitbx::math::iround(k0);
   const int il0=scitbx::math::iround(l0);
   cctbx::miller::index<long> h(scitbx::vec3<long>(ih0,ik0,il0));
   cctbx::miller::sym_equiv_indices sei(this->GetCCTbxSpg(),h);
   int nbEquiv;
   if(forceFriedelLaw) nbEquiv=sei.multiplicity(false);
   else nbEquiv=sei.multiplicity(true);
   if( ((this->IsCentrosymmetric())||forceFriedelLaw) && excludeFriedelMate) nbEquiv/=2;
   CrystMatrix_REAL equivReflList(nbEquiv,5);
   complex<double> sf0((double)(sf_re),(double)(sf_im));
   for(int i=0;i<nbEquiv;i+=1)
   {
      cctbx::miller::index<long> k = sei(i).h();
      equivReflList(i,0)=(REAL)k[0];
      equivReflList(i,1)=(REAL)k[1];
      equivReflList(i,2)=(REAL)k[2];
      equivReflList(i,3)=(REAL)sei(i).complex_eq(sf0).real();
      equivReflList(i,4)=(REAL)sei(i).complex_eq(sf0).imag();
   }
   VFN_DEBUG_EXIT("SpaceGroup::GetAllEquivRefl():",5)
   return equivReflList;
}

bool SpaceGroup::IsReflSystematicAbsent(const REAL h0, const REAL k0, const REAL l0)const
{
   const int ih0=scitbx::math::iround(h0);
   const int ik0=scitbx::math::iround(k0);
   const int il0=scitbx::math::iround(l0);
   cctbx::miller::index<long> h(scitbx::vec3<long>(ih0,ik0,il0));
   return this->GetCCTbxSpg().is_sys_absent(h);
}

bool SpaceGroup::IsReflCentric(const REAL h0, const REAL k0, const REAL l0)const
{
   const int ih0=scitbx::math::iround(h0);
   const int ik0=scitbx::math::iround(k0);
   const int il0=scitbx::math::iround(l0);
   cctbx::miller::index<long> h(scitbx::vec3<long>(ih0,ik0,il0));
   return this->GetCCTbxSpg().is_centric(h);
}

unsigned int SpaceGroup::GetExpectedIntensityFactor(const REAL h0,
                                                    const REAL k0,
                                                    const REAL l0)const
{
   const int ih0=scitbx::math::iround(h0);
   const int ik0=scitbx::math::iround(k0);
   const int il0=scitbx::math::iround(l0);
   cctbx::miller::index<long> h(scitbx::vec3<long>(ih0,ik0,il0));
   return this->GetCCTbxSpg().epsilon(h);
}

void SpaceGroup::InitSpaceGroup(const string &spgId)
{
   if((mId==spgId)&&(mpCCTbxSpaceGroup!=0)) return;
   VFN_DEBUG_ENTRY("SpaceGroup::InitSpaceGroup():"<<spgId,8)
   #ifdef __DEBUG__
   (*fpObjCrystInformUser)("Initializing spacegroup: "+spgId);
   #endif
   try
   {
      cctbx::sgtbx::space_group_symbols sgs=cctbx::sgtbx::space_group_symbols(spgId);
      if(mpCCTbxSpaceGroup!=0) delete mpCCTbxSpaceGroup;
      mpCCTbxSpaceGroup=0;
      mpCCTbxSpaceGroup = new cctbx::sgtbx::space_group(sgs);
   }
   catch(exception &ex1)
   {
      try
      {
         (*fpObjCrystInformUser)("Failed lookup symbol ! try Hall symbol ?");
         if(mpCCTbxSpaceGroup!=0) delete mpCCTbxSpaceGroup;
         mpCCTbxSpaceGroup=0;
         mpCCTbxSpaceGroup = new cctbx::sgtbx::space_group(spgId);
      }
      catch(exception &ex2)
      {
         (*fpObjCrystInformUser)("Could not interpret Spacegroup Symbol:"+spgId);
         this->InitSpaceGroup(mId);
         VFN_DEBUG_EXIT("SpaceGroup::InitSpaceGroup() could not interpret spacegroup:"<<spgId<<":"<<ex1.what()<<":"<<ex2.what(),8)
         return;
      }
   }

   try
   {
      //Inversion center
      if(this->GetCCTbxSpg().f_inv() == 2)
      {
         mHasInversionCenter=true ;
         if( (this->GetCCTbxSpg().inv_t()[0] !=0) ||
            (this->GetCCTbxSpg().inv_t()[1] !=0) ||
            (this->GetCCTbxSpg().inv_t()[2] !=0)   ) mIsInversionCenterAtOrigin=false;
         else mIsInversionCenterAtOrigin=true;
      }
      else
      {
         mHasInversionCenter=false ;
         mIsInversionCenterAtOrigin=true;
      }

      //initialize asymmetric unit
      mAsymmetricUnit.SetSpaceGroup(*this);

      mUniqueAxisId=0;
      if(  (this->GetCCTbxSpg().type().number() >2)
         &&(this->GetCCTbxSpg().type().number() <16))
      {
         string ch=this->GetCCTbxSpg().match_tabulated_settings().hall();
         if(ch.find("x")!=std::string::npos) {mUniqueAxisId=0;}
         else
            if(ch.find("y")!=std::string::npos) {mUniqueAxisId=1;}
            else mUniqueAxisId=2;
      }

      mNbSym    =this->GetCCTbxSpg().n_smx();
      mNbTrans  =this->GetCCTbxSpg().n_ltr();
      mSpgNumber=this->GetCCTbxSpg().type().number();

      mExtension='\0'; //this->GetCCTbxSpg().type().extension();
   }
   catch(exception &ex)
   {
      (*fpObjCrystInformUser)("Error initializing spacegroup (Incorrect Hall symbol ?):"+spgId);
      this->InitSpaceGroup(mId);
      VFN_DEBUG_EXIT("SpaceGroup::InitSpaceGroup() could not interpret spacegroup:"<<spgId<<":"<<ex.what(),8)
      return;
   }

   mExtension=this->GetCCTbxSpg().match_tabulated_settings().extension();

   // Force using the H-M symbol
   if(mExtension=='\0') mId=this->GetCCTbxSpg().match_tabulated_settings().hermann_mauguin();
   else                 mId=this->GetCCTbxSpg().match_tabulated_settings().hermann_mauguin()+":"+mExtension;

   // Store rotation matrices & translation vectors
   mvTrans.resize(mNbTrans);
   for(unsigned int i=0;i<mNbTrans;i++)
   {
      for(unsigned int j=0;j<3;j++)
         mvTrans[i].tr[j] = this->GetCCTbxSpg().ltr(i)[j]/(REAL)(this->GetCCTbxSpg().ltr(i).den());
   }
   mvSym.resize(mNbSym);
   for(unsigned int j=0;j<mNbSym;j++)
   {
      const cctbx::sgtbx::rt_mx *pMatrix=&(this->GetCCTbxSpg().smx(j));
      const cctbx::sgtbx::rot_mx *pRot=&(pMatrix->r());
      const cctbx::sgtbx::tr_vec *pTrans=&(pMatrix->t());
      const REAL r_den=1/(REAL)(pMatrix->r().den());
      const REAL t_den=1/(REAL)(pMatrix->t().den());
      for(unsigned int i=0;i<9;++i) mvSym[j].mx[i]=(*pRot)[i]*r_den;
      for(unsigned int i=0;i<3;++i) mvSym[j].tr[i]=(*pTrans)[i]*t_den;
   }
   #ifdef __DEBUG__
   this->Print();
   #endif
   mClock.Click();
   string extension("");
   if(mExtension=='1') extension=" (Using origin choice #1)";
   if(mExtension=='2') extension=" (Using origin choice #2)";
   if(mExtension=='R') extension=" (Using Rhombohedral cell)";
   if(mExtension=='H') extension=" (Using Hexagonal cell)";
   #ifdef __DEBUG__
  (*fpObjCrystInformUser)("Initialized spacegroup, HM: "+spgId+extension+" , Hall:"+this->GetCCTbxSpg().type().hall_symbol());
   #endif
   VFN_DEBUG_EXIT("SpaceGroup::InitSpaceGroup():"<<spgId,8)
}

}//namespace
