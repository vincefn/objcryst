/*  ObjCryst++ Object-Oriented Crystallographic Library
    (c) 2000-2009 Vincent Favre-Nicolin vincefn@users.sourceforge.net
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
/*   Molecule.cpp
*  Source file for the Molecule scatterer
*
*/
#include <sstream>
#include <fstream>
#include <iterator>
#include <algorithm>

#include "ObjCryst/Quirks/VFNStreamFormat.h"
#include "ObjCryst/ObjCryst/Molecule.h"
#include "ObjCryst/RefinableObj/GlobalOptimObj.h"

#ifdef OBJCRYST_GL
   #ifdef __DARWIN__
      #include <OpenGL/glu.h>
   #else
      #include <GL/glu.h>
   #endif
#endif

#ifdef __WX__CRYST__
   #include "ObjCryst/wxCryst/wxMolecule.h"
#endif

//#include <xmmintrin.h>

// Try new approach for rigid bodies ?
#define RIGID_BODY_STRICT_EXPERIMENTAL

// Tighter restraints x**2+x**4+x**6 instead of just X**2
#undef RESTRAINT_X2_X4_X6

using namespace std;


namespace ObjCryst
{
XYZ::XYZ(REAL x0,REAL y0,REAL z0):x(x0),y(y0),z(z0){};
   
REAL GetBondLength(const MolAtom&at1,const MolAtom&at2)
{
   //TAU_PROFILE("GetBondLength()","REAL (...)",TAU_DEFAULT);
   
   /*
   __m128 m128=_mm_set_ps(0.0f,
                          at1.GetZ()-at2.GetZ(),
                          at1.GetY()-at2.GetY(),
                          at1.GetX()-at2.GetX());

   __m128 a = _mm_mul_ps(m128,m128);

   // horizontal add
   __m128 b = _mm_add_ss(_mm_shuffle_ps(a, a, _MM_SHUFFLE(0,0,0,0)),
                         _mm_add_ss(_mm_shuffle_ps(a, a, _MM_SHUFFLE(1,1,1,1)),
                                    _mm_shuffle_ps(a, a, _MM_SHUFFLE(2,2,2,2))));
   union m128_float
   {
      __m128 m128;
      struct
      {
         float x, y, z, pad;
      };
   };
   union m128_float l;
   l.m128 = _mm_sqrt_ss(b);
   
   return l.x;
   */
   return  sqrt( (at1.GetX()-at2.GetX())
                *(at1.GetX()-at2.GetX())
                +(at1.GetY()-at2.GetY())
                *(at1.GetY()-at2.GetY())
                +(at1.GetZ()-at2.GetZ())
                *(at1.GetZ()-at2.GetZ()) );
}
REAL GetBondAngle(const MolAtom &at1,const MolAtom &at2,const MolAtom &at3)
{
   //TAU_PROFILE("GetBondAngle()","REAL (...)",TAU_DEFAULT);
   const REAL x21=at1.GetX()-at2.GetX();
   const REAL y21=at1.GetY()-at2.GetY();
   const REAL z21=at1.GetZ()-at2.GetZ();
   const REAL x23=at3.GetX()-at2.GetX();
   const REAL y23=at3.GetY()-at2.GetY();
   const REAL z23=at3.GetZ()-at2.GetZ();
   const REAL norm21_norm23= sqrt( (x21*x21+y21*y21+z21*z21)
                                  *(x23*x23+y23*y23+z23*z23)+1e-6);
   const REAL angle=(x21*x23+y21*y23+z21*z23)/norm21_norm23;
   if(angle>=1)  return 0;
   if(angle<=-1) return M_PI;
   return acos(angle);
}
REAL GetDihedralAngle(const MolAtom &at1,const MolAtom &at2,const MolAtom &at3,const MolAtom &at4)
{
   //TAU_PROFILE("GetDihedralAngle()","REAL (...)",TAU_DEFAULT);
   const REAL x21=at1.GetX()-at2.GetX();
   const REAL y21=at1.GetY()-at2.GetY();
   const REAL z21=at1.GetZ()-at2.GetZ();
   
   const REAL x34=at4.GetX()-at3.GetX();
   const REAL y34=at4.GetY()-at3.GetY();
   const REAL z34=at4.GetZ()-at3.GetZ();
   
   const REAL x23=at3.GetX()-at2.GetX();
   const REAL y23=at3.GetY()-at2.GetY();
   const REAL z23=at3.GetZ()-at2.GetZ();
   
   // v21 x v23
   const REAL x123= y21*z23-z21*y23;
   const REAL y123= z21*x23-x21*z23;
   const REAL z123= x21*y23-y21*x23;
   
   // v32 x v34 (= -v23 x v34)
   const REAL x234= -(y23*z34-z23*y34);
   const REAL y234= -(z23*x34-x23*z34);
   const REAL z234= -(x23*y34-y23*x34);

   const REAL norm123_norm234= sqrt( (x123*x123+y123*y123+z123*z123)
                                    *(x234*x234+y234*y234+z234*z234)+1e-6);
   
   REAL angle=(x123*x234+y123*y234+z123*z234)/norm123_norm234;
   if(angle>= 1) angle=0;
   else 
   {
      if(angle<=-1) angle=M_PI;
      else angle=acos(angle);
   }
   if((x21*x234 + y21*y234 + z21*z234)<0) return -angle;
   return angle;
}

   
void ExpandAtomGroupRecursive(MolAtom* atom,
                              const map<MolAtom*,set<MolAtom*> > &connect,
                              set<MolAtom*> &atomlist,const MolAtom* finalAtom)
{
   const pair<set<MolAtom*>::iterator,bool> status=atomlist.insert(atom);
   if(false==status.second) return;
   if(finalAtom==atom) return;
   map<MolAtom*,set<MolAtom*> >::const_iterator c=connect.find(atom);
   set<MolAtom*>::const_iterator pos;
   for(pos=c->second.begin();pos!=c->second.end();++pos)
   {
      ExpandAtomGroupRecursive(*pos,connect,atomlist,finalAtom);
   }
}

void ExpandAtomGroupRecursive(MolAtom* atom,
                              const map<MolAtom*,set<MolAtom*> > &connect,
                              map<MolAtom*,unsigned long> &atomlist,const unsigned long maxdepth,unsigned long depth)
{
   if(atomlist.count(atom)>0)
     if(atomlist[atom]<=depth) return;
   atomlist[atom]=depth;
   if(depth==maxdepth) return;//maxdepth reached
   map<MolAtom*,set<MolAtom*> >::const_iterator c=connect.find(atom);
   set<MolAtom*>::const_iterator pos;
   for(pos=c->second.begin();pos!=c->second.end();++pos)
   {
      ExpandAtomGroupRecursive(*pos,connect,atomlist,maxdepth,depth+1);
   }
}

//######################################################################
//
//      MolAtom
//
//######################################################################

MolAtom::MolAtom(const REAL x, const REAL y, const REAL z,
                 const ScatteringPower *pPow, const string &name, Molecule &parent):
mName(name),mX(x),mY(y),mZ(z),mOccupancy(1.),mpScattPow(pPow),mpMol(&parent),mIsInRing(false)
#ifdef __WX__CRYST__
,mpWXCrystObj(0)
#endif
{
   VFN_DEBUG_MESSAGE("MolAtom::MolAtom()",4)
}

MolAtom::~MolAtom()
{
#ifdef __WX__CRYST__
this->WXDelete();
#endif
}

void MolAtom::SetName(const string &name)
{
   mName=name;
   // Set parameter's name in the parent molecule
   // The atom should already be part of a Molecule, but just in case be careful
   if((mName!="") && (mpMol!=0))
   {
      try
      {
         mpMol->GetPar(&mX).SetName(mpMol->GetName()+"_"+mName+"_x");
         mpMol->GetPar(&mY).SetName(mpMol->GetName()+"_"+mName+"_y");
         mpMol->GetPar(&mZ).SetName(mpMol->GetName()+"_"+mName+"_z");
      }
      catch(const ObjCrystException &except)
      {
         cout<<"MolAtom::SetName(): Atom parameters not yet declared in a Molecule ?"<<endl;
      }
   }
}

const string& MolAtom::GetName()const{return mName;}
      string& MolAtom::GetName()     {return mName;}
      
const Molecule& MolAtom::GetMolecule()const{return *mpMol;}
      Molecule& MolAtom::GetMolecule()     {return *mpMol;}

const REAL& MolAtom::X()const{return mX;}
const REAL& MolAtom::Y()const{return mY;}
const REAL& MolAtom::Z()const{return mZ;}

REAL& MolAtom::X(){return mX;}
REAL& MolAtom::Y(){return mY;}
REAL& MolAtom::Z(){return mZ;}

REAL MolAtom::GetX()const{return mX;}
REAL MolAtom::GetY()const{return mY;}
REAL MolAtom::GetZ()const{return mZ;}
REAL MolAtom::GetOccupancy()const{return mOccupancy;}

void MolAtom::SetX(const REAL a)const{ mX=a;mpMol->GetAtomPositionClock().Click();}
void MolAtom::SetY(const REAL a)const{ mY=a;mpMol->GetAtomPositionClock().Click();}
void MolAtom::SetZ(const REAL a)const{ mZ=a;mpMol->GetAtomPositionClock().Click();}
void MolAtom::SetOccupancy(const REAL a){ mOccupancy=a;}

bool MolAtom::IsDummy()const{return mpScattPow==0;}
const ScatteringPower& MolAtom::GetScatteringPower()const{return *mpScattPow;}
void MolAtom::SetScatteringPower(const ScatteringPower& pow){mpScattPow=&pow;}

void MolAtom::XMLOutput(ostream &os,int indent)const
{
   VFN_DEBUG_ENTRY("MolAtom::XMLOutput()",4)
   for(int i=0;i<indent;i++) os << "  " ;
   XMLCrystTag tag("Atom",false,true);
   tag.AddAttribute("Name",this->GetName());
   if(!this->IsDummy())tag.AddAttribute("ScattPow",this->GetScatteringPower().GetName());
   {
      stringstream ss;
      ss.precision(os.precision());
      ss <<mX;
      tag.AddAttribute("x",ss.str());
   }
   {
      stringstream ss;
      ss.precision(os.precision());
      ss <<mY;
      tag.AddAttribute("y",ss.str());
   }
   {
      stringstream ss;
      ss.precision(os.precision());
      ss <<mZ;
      tag.AddAttribute("z",ss.str());
   }
   {
      stringstream ss;
      ss.precision(os.precision());
      ss <<mOccupancy;
      tag.AddAttribute("Occup",ss.str());
   }
   os <<tag<<endl;
   VFN_DEBUG_EXIT("MolAtom::XMLOutput()",4)
}

void MolAtom::XMLInput(istream &is,const XMLCrystTag &tag)
{
   VFN_DEBUG_ENTRY("MolAtom::XMLInput()",7)
   string name;
   for(unsigned int i=0;i<tag.GetNbAttribute();i++)
   {
      if("Name"==tag.GetAttributeName(i))
      {
         name=tag.GetAttributeValue(i);
      }
      if("ScattPow"==tag.GetAttributeName(i))
      {
         mpScattPow=&(mpMol->GetCrystal().GetScatteringPower(tag.GetAttributeValue(i)));
      }
      if("x"==tag.GetAttributeName(i))
      {
         stringstream ss(tag.GetAttributeValue(i));
         ss >>mX;
      }
      if("y"==tag.GetAttributeName(i))
      {
         stringstream ss(tag.GetAttributeValue(i));
         ss >>mY;
      }
      if("z"==tag.GetAttributeName(i))
      {
         stringstream ss(tag.GetAttributeValue(i));
         ss >>mZ;
      }
      if("Occup"==tag.GetAttributeName(i))
      {
         stringstream ss(tag.GetAttributeValue(i));
         ss >>mOccupancy;
      }
   }
   this->SetName(name);
   VFN_DEBUG_EXIT("MolAtom::XMLInput()",7)
}

void MolAtom::SetIsInRing(const bool r)const{mIsInRing=r;}
bool MolAtom::IsInRing()const{return mIsInRing;}

#ifdef __WX__CRYST__
WXCrystObjBasic* MolAtom::WXCreate(wxWindow* parent)
{
   VFN_DEBUG_ENTRY("MolAtom::WXCreate()",5)
   mpWXCrystObj=new WXMolAtom(parent,this);
   VFN_DEBUG_EXIT("MolAtom::WXCreate()",5)
   return mpWXCrystObj;
}
WXCrystObjBasic* MolAtom::WXGet(){return mpWXCrystObj;}
void MolAtom::WXDelete(){if(0!=mpWXCrystObj) delete mpWXCrystObj;mpWXCrystObj=0;}
void MolAtom::WXNotifyDelete(){mpWXCrystObj=0;}
#endif
//######################################################################
//
//      MolBond
//
//######################################################################
MolBond::MolBond(MolAtom &atom1, MolAtom &atom2,
                 const REAL length0, const REAL sigma, const REAL delta,
                 Molecule &parent,const REAL bondOrder):
mAtomPair(make_pair(&atom1,&atom2)),
mLength0(length0),mDelta(delta),mSigma(sigma),
mBondOrder(bondOrder),mIsFreeTorsion(false),mpMol(&parent)
#ifdef __WX__CRYST__
,mpWXCrystObj(0)
#endif
{}

MolBond::~MolBond()
{
#ifdef __WX__CRYST__
this->WXDelete();
#endif
}

const Molecule& MolBond::GetMolecule()const{return *mpMol;}
      Molecule& MolBond::GetMolecule()     {return *mpMol;}

string MolBond::GetName()const
{return this->GetAtom1().GetName()+"-"+this->GetAtom2().GetName();}

void MolBond::XMLOutput(ostream &os,int indent)const
{
   VFN_DEBUG_ENTRY("MolBond::XMLOutput()",4)
   for(int i=0;i<indent;i++) os << "  " ;
   XMLCrystTag tag("Bond",false,true);
   tag.AddAttribute("Atom1",mAtomPair.first->GetName());
   tag.AddAttribute("Atom2",mAtomPair.second->GetName());
   {
      stringstream ss;
      ss.precision(os.precision());
      ss <<mLength0;
      tag.AddAttribute("Length",ss.str());
   }
   {
      stringstream ss;
      ss.precision(os.precision());
      ss <<mDelta;
      tag.AddAttribute("Delta",ss.str());
   }
   {
      stringstream ss;
      ss.precision(os.precision());
      ss <<mSigma;
      tag.AddAttribute("Sigma",ss.str());
   }
   {
      stringstream ss;
      ss.precision(os.precision());
      ss <<mBondOrder;
      tag.AddAttribute("BondOrder",ss.str());
   }
   {
      stringstream ss;
      ss.precision(os.precision());
      ss <<mIsFreeTorsion;
      tag.AddAttribute("FreeTorsion",ss.str());
   }
   os <<tag<<endl;
   VFN_DEBUG_EXIT("MolBond::XMLOutput()",4)
}

void MolBond::XMLInput(istream &is,const XMLCrystTag &tag)
{
   VFN_DEBUG_ENTRY("MolBond::XMLInput():",7)
   for(unsigned int i=0;i<tag.GetNbAttribute();i++)
   {
      if("Atom1"==tag.GetAttributeName(i))
      {
         mAtomPair.first=&(mpMol->GetAtom(tag.GetAttributeValue(i)));
      }
      if("Atom2"==tag.GetAttributeName(i))
      {
         mAtomPair.second=&(mpMol->GetAtom(tag.GetAttributeValue(i)));
      }
      if("Length"==tag.GetAttributeName(i))
      {
         stringstream ss(tag.GetAttributeValue(i));
         ss >>mLength0;
      }
      if("Delta"==tag.GetAttributeName(i))
      {
         stringstream ss(tag.GetAttributeValue(i));
         ss >>mDelta;
      }
      if("Sigma"==tag.GetAttributeName(i))
      {
         stringstream ss(tag.GetAttributeValue(i));
         ss >>mSigma;
      }
      if("BondOrder"==tag.GetAttributeName(i))
      {
         stringstream ss(tag.GetAttributeValue(i));
         ss >>mBondOrder;
      }
      if("FreeTorsion"==tag.GetAttributeName(i))
      {
         stringstream ss(tag.GetAttributeValue(i));
         ss >>mIsFreeTorsion;
      }
   }
   VFN_DEBUG_EXIT("MolBond::XMLInput():",7)
}

REAL MolBond::GetLogLikelihood()const{return this->GetLogLikelihood(false,true);}

REAL MolBond::GetLogLikelihood(const bool calcDeriv, const bool recalc)const
{
   if(!recalc) return mLLK;
   VFN_DEBUG_ENTRY("MolBond::GetLogLikelihood():",2)
   //TAU_PROFILE("MolBond::GetLogLikelihood()","REAL (bool,bool)",TAU_DEFAULT);
   //const REAL length=this->GetLength();
   const REAL x=this->GetAtom2().GetX()-this->GetAtom1().GetX();
   const REAL y=this->GetAtom2().GetY()-this->GetAtom1().GetY();
   const REAL z=this->GetAtom2().GetZ()-this->GetAtom1().GetZ();
   const REAL length=sqrt(abs(x*x+y*y+z*z));
   
   if(calcDeriv)
   {
      const REAL tmp2=1/(length+1e-10);
      mDerivAtom1.x=-x*tmp2;
      mDerivAtom1.y=-y*tmp2;
      mDerivAtom1.z=-z*tmp2;

      mDerivAtom2.x=-mDerivAtom1.x;
      mDerivAtom2.y=-mDerivAtom1.y;
      mDerivAtom2.z=-mDerivAtom1.z;
   }

   if(mSigma<1e-6)
   {
      if(calcDeriv) mDerivLLKCoeff=0;
      mLLK=0;
      return 0;
   }
   mLLK=length-(mLength0+mDelta);
   if(mLLK>0)
   {
      mLLK /= mSigma;
      if(calcDeriv) mDerivLLKCoeff=2*mLLK/mSigma;
      #ifdef RESTRAINT_X2_X4_X6
      const float mLLK2=mLLK*mLLK;
      //if(calcDeriv) mDerivLLKCoeff=(2*mLLK+4*mLLK2+6*mLLK2*mLLK2)/mSigma;
      //mLLK=mLLK2*(1+mLLK2+mLLK2*mLLK2);
      if(calcDeriv) mDerivLLKCoeff=(2*mLLK+4*mLLK2)/mSigma;
      mLLK=mLLK2*(1+mLLK2);
      #else
      if(calcDeriv) mDerivLLKCoeff=2*mLLK/mSigma;
      mLLK *= mLLK;
      #endif
      VFN_DEBUG_EXIT("MolBond::GetLogLikelihood():",2)
      return mLLK;
   }
   mLLK=length-(mLength0-mDelta);
   if(mLLK<0)
   {
      mLLK /= mSigma;
      #ifdef RESTRAINT_X2_X4_X6
      const float mLLK2=mLLK*mLLK;
      //if(calcDeriv) mDerivLLKCoeff=(2*mLLK+4*mLLK2+6*mLLK2*mLLK2)/mSigma;
      //mLLK=mLLK2*(1+mLLK2+mLLK2*mLLK2);
      if(calcDeriv) mDerivLLKCoeff=(2*mLLK+4*mLLK2)/mSigma;
      mLLK=mLLK2*(1+mLLK2);
      #else
      if(calcDeriv) mDerivLLKCoeff=2*mLLK/mSigma;
      mLLK *= mLLK;
      #endif
      VFN_DEBUG_EXIT("MolBond::GetLogLikelihood():",2)
      return mLLK;
   }
   if(calcDeriv) mDerivLLKCoeff=0;
   mLLK=0;
   VFN_DEBUG_EXIT("MolBond::GetLogLikelihood():",2)
   return mLLK;
}

REAL MolBond::GetDeriv(const map<const MolAtom*,XYZ> &m, const bool llk)const
{
   //TAU_PROFILE("MolBond::GetDeriv()","REAL (mak,bool)",TAU_DEFAULT);
   REAL d=0;
   map<const MolAtom*,XYZ>::const_iterator pos;
   pos=m.find(mAtomPair.first);
   if(pos!=m.end())
      d+= pos->second.x*mDerivAtom1.x
         +pos->second.y*mDerivAtom1.y
         +pos->second.z*mDerivAtom1.z;
   pos=m.find(mAtomPair.second);
   if(pos!=m.end())
      d+= pos->second.x*mDerivAtom2.x
         +pos->second.y*mDerivAtom2.y
         +pos->second.z*mDerivAtom2.z;
   if(llk) return mDerivLLKCoeff*d;
   return d;
}

void MolBond::CalcGradient(std::map<MolAtom*,XYZ> &m)const
{
   this->GetLogLikelihood(true,true);
   map<MolAtom*,XYZ>::iterator pos;
   pos=m.find(mAtomPair.first);
   if(pos!=m.end())
   {
      pos->second.x+=mDerivLLKCoeff*mDerivAtom1.x;
      pos->second.y+=mDerivLLKCoeff*mDerivAtom1.y;
      pos->second.z+=mDerivLLKCoeff*mDerivAtom1.z;
   }
   pos=m.find(mAtomPair.second);
   if(pos!=m.end())
   {
      pos->second.x+=mDerivLLKCoeff*mDerivAtom2.x;
      pos->second.y+=mDerivLLKCoeff*mDerivAtom2.y;
      pos->second.z+=mDerivLLKCoeff*mDerivAtom2.z;
   }
   #if 0
   // Display gradient - for tests
   cout<<this->GetName()<<" :LLK"<<endl;
   for(map<MolAtom*,XYZ>::const_iterator pos=m.begin();pos!=m.end();++pos)
   {
      char buf[100];
      sprintf(buf,"%10s Grad LLK= (%8.3f %8.3f %8.3f)",
            pos->first->GetName().c_str(),pos->second.x,pos->second.y,pos->second.z);
      cout<<buf<<endl;
   }
   #endif
}

const MolAtom& MolBond::GetAtom1()const{return *(mAtomPair.first);}
const MolAtom& MolBond::GetAtom2()const{return *(mAtomPair.second);}
MolAtom& MolBond::GetAtom1(){return *(mAtomPair.first);}
MolAtom& MolBond::GetAtom2(){return *(mAtomPair.second);}
void MolBond::SetAtom1(MolAtom &at){mAtomPair.first =&at;}
void MolBond::SetAtom2(MolAtom &at){mAtomPair.second=&at;}
REAL MolBond::GetLength()const
{
   return GetBondLength(GetAtom1(),this->GetAtom2());
}

REAL MolBond::GetLength0()const{return mLength0;}
REAL MolBond::GetLengthDelta()const{return mDelta;}
REAL MolBond::GetLengthSigma()const{return mSigma;}
REAL MolBond::GetBondOrder()const{return mBondOrder;}

REAL& MolBond::Length0(){return mLength0;}
REAL& MolBond::LengthDelta(){return mDelta;}
REAL& MolBond::LengthSigma(){return mSigma;}
REAL& MolBond::BondOrder(){return mBondOrder;}

void MolBond::SetLength0(const REAL a){mLength0=a;}
void MolBond::SetLengthDelta(const REAL a){mDelta=a;}
void MolBond::SetLengthSigma(const REAL a){mSigma=a;}
void MolBond::SetBondOrder(const REAL a){mBondOrder=a;}

bool MolBond::IsFreeTorsion()const{return mIsFreeTorsion;}
void MolBond::SetFreeTorsion(const bool isFreeTorsion)
{
   if(mIsFreeTorsion==isFreeTorsion) return;
   mIsFreeTorsion=isFreeTorsion;
   mpMol->GetBondListClock().Click();
}
#ifdef __WX__CRYST__
WXCrystObjBasic* MolBond::WXCreate(wxWindow* parent)
{
   VFN_DEBUG_ENTRY("MolBond::WXCreate()",5)
   mpWXCrystObj=new WXMolBond(parent,this);
   VFN_DEBUG_EXIT("MolBond::WXCreate()",5)
   return mpWXCrystObj;
}
WXCrystObjBasic* MolBond::WXGet(){return mpWXCrystObj;}
void MolBond::WXDelete(){if(0!=mpWXCrystObj) delete mpWXCrystObj;mpWXCrystObj=0;}
void MolBond::WXNotifyDelete(){mpWXCrystObj=0;}
#endif
//######################################################################
//
//      MolBondAngle
//
//######################################################################
MolBondAngle::MolBondAngle(MolAtom &atom1,MolAtom &atom2,MolAtom &atom3,
                           const REAL angle, const REAL sigma, const REAL delta,
                           Molecule &parent):
mAngle0(angle),mDelta(delta),mSigma(sigma),mpMol(&parent)
#ifdef __WX__CRYST__
,mpWXCrystObj(0)
#endif
{
   mvpAtom.push_back(&atom1);
   mvpAtom.push_back(&atom2);
   mvpAtom.push_back(&atom3);
}

MolBondAngle::~MolBondAngle()
{
#ifdef __WX__CRYST__
this->WXDelete();
#endif
}

const Molecule& MolBondAngle::GetMolecule()const{return *mpMol;}
      Molecule& MolBondAngle::GetMolecule()     {return *mpMol;}

string MolBondAngle::GetName()const
{
   return this->GetAtom1().GetName()+"-"
         +this->GetAtom2().GetName()+"-"
         +this->GetAtom3().GetName();
}

void MolBondAngle::XMLOutput(ostream &os,int indent)const
{
   VFN_DEBUG_ENTRY("MolBondAngle::XMLOutput()",4)
   for(int i=0;i<indent;i++) os << "  " ;
   XMLCrystTag tag("BondAngle",false,true);
   tag.AddAttribute("Atom1",this->GetAtom1().GetName());
   tag.AddAttribute("Atom2",this->GetAtom2().GetName());
   tag.AddAttribute("Atom3",this->GetAtom3().GetName());
   {
      stringstream ss;
      ss.precision(os.precision());
      ss <<mAngle0*RAD2DEG;
      tag.AddAttribute("Angle",ss.str());
   }
   {
      stringstream ss;
      ss.precision(os.precision());
      ss <<mDelta*RAD2DEG;
      tag.AddAttribute("Delta",ss.str());
   }
   {
      stringstream ss;
      ss.precision(os.precision());
      ss <<mSigma*RAD2DEG;
      tag.AddAttribute("Sigma",ss.str());
   }
   os <<tag<<endl;
   VFN_DEBUG_EXIT("MolBondAngle::XMLOutput()",4)
}

void MolBondAngle::XMLInput(istream &is,const XMLCrystTag &tag)
{
   VFN_DEBUG_ENTRY("MolBondAngle::XMLInput():",4)
   mvpAtom.resize(3);
   for(unsigned int i=0;i<tag.GetNbAttribute();i++)
   {
      if("Atom1"==tag.GetAttributeName(i))
      {
         mvpAtom[0]=&(mpMol->GetAtom(tag.GetAttributeValue(i)));
      }
      if("Atom2"==tag.GetAttributeName(i))
      {
         mvpAtom[1]=&(mpMol->GetAtom(tag.GetAttributeValue(i)));
      }
      if("Atom3"==tag.GetAttributeName(i))
      {
         mvpAtom[2]=&(mpMol->GetAtom(tag.GetAttributeValue(i)));
      }
      if("Angle"==tag.GetAttributeName(i))
      {
         stringstream ss(tag.GetAttributeValue(i));
         ss >>mAngle0;
         mAngle0*=DEG2RAD;
      }
      if("Delta"==tag.GetAttributeName(i))
      {
         stringstream ss(tag.GetAttributeValue(i));
         ss >>mDelta;
         mDelta*=DEG2RAD;
      }
      if("Sigma"==tag.GetAttributeName(i))
      {
         stringstream ss(tag.GetAttributeValue(i));
         ss >>mSigma;
         mSigma*=DEG2RAD;
      }
   }
   VFN_DEBUG_EXIT("MolBondAngle::XMLInput():",4)
}
REAL& MolBondAngle::Angle0()
{
   return mAngle0;
}
REAL& MolBondAngle::AngleDelta(){return mDelta;}
REAL& MolBondAngle::AngleSigma(){return mSigma;}

REAL MolBondAngle::GetAngle0()const{return mAngle0;}
REAL MolBondAngle::GetAngleDelta()const{return mDelta;}
REAL MolBondAngle::GetAngleSigma()const{return mSigma;}

void MolBondAngle::SetAngle0(const REAL angle){mAngle0=angle;}
void MolBondAngle::SetAngleDelta(const REAL delta){mDelta=delta;}
void MolBondAngle::SetAngleSigma(const REAL sigma){mSigma=sigma;}

REAL MolBondAngle::GetAngle()const
{
   return GetBondAngle(this->GetAtom1(),this->GetAtom2(),this->GetAtom3());
}

REAL MolBondAngle::GetLogLikelihood()const{return this->GetLogLikelihood(false,true);}

REAL MolBondAngle::GetLogLikelihood(const bool calcDeriv, const bool recalc)const
{
   if(!recalc) return mLLK;
   VFN_DEBUG_ENTRY("MolBondAngle::GetLogLikelihood():",2)
   //TAU_PROFILE("MolBondAngle::GetLogLikelihood()","REAL (bool,bool)",TAU_DEFAULT);
   //const REAL angle=this->GetAngle();
   const REAL x21=this->GetAtom1().GetX()-this->GetAtom2().GetX();
   const REAL y21=this->GetAtom1().GetY()-this->GetAtom2().GetY();
   const REAL z21=this->GetAtom1().GetZ()-this->GetAtom2().GetZ();
   const REAL x23=this->GetAtom3().GetX()-this->GetAtom2().GetX();
   const REAL y23=this->GetAtom3().GetY()-this->GetAtom2().GetY();
   const REAL z23=this->GetAtom3().GetZ()-this->GetAtom2().GetZ();
   
   const REAL n1=sqrt(abs(x21*x21+y21*y21+z21*z21));
   const REAL n3=sqrt(abs(x23*x23+y23*y23+z23*z23));
   const REAL p=x21*x23+y21*y23+z21*z23;
   
   const REAL a0=p/(n1*n3+1e-10);
   REAL angle;
   if(a0>=1)  angle=0;
   else
   {
      if(a0<=-1) angle=M_PI;
      else angle=acos(a0);
   }
   
   if(calcDeriv)
   {
      const REAL s=1/(sqrt(1-a0*a0+1e-6));
      const REAL s0=-s/(n1*n3+1e-10);
      const REAL s1= s*p/(n3*n1*n1*n1+1e-10);
      const REAL s3= s*p/(n1*n3*n3*n3+1e-10);
      mDerivAtom1.x=s0*x23+s1*x21;
      mDerivAtom1.y=s0*y23+s1*y21;
      mDerivAtom1.z=s0*z23+s1*z21;

      mDerivAtom3.x=s0*x21+s3*x23;
      mDerivAtom3.y=s0*y21+s3*y23;
      mDerivAtom3.z=s0*z21+s3*z23;

      mDerivAtom2.x=-(mDerivAtom1.x+mDerivAtom3.x);
      mDerivAtom2.y=-(mDerivAtom1.y+mDerivAtom3.y);
      mDerivAtom2.z=-(mDerivAtom1.z+mDerivAtom3.z);
   }

   if(mSigma<1e-6)
   {
      if(calcDeriv) mDerivLLKCoeff=0;
      mLLK=0;
      return 0;
   }
   
   mLLK=angle-(mAngle0+mDelta);
   if(mLLK>0)
   {
      mLLK/=mSigma;
      #ifdef RESTRAINT_X2_X4_X6
      const float mLLK2=mLLK*mLLK;
      //if(calcDeriv) mDerivLLKCoeff=(2*mLLK+4*mLLK2+6*mLLK2*mLLK2)/mSigma;
      //mLLK=mLLK2*(1+mLLK2+mLLK2*mLLK2);
      if(calcDeriv) mDerivLLKCoeff=(2*mLLK+4*mLLK2)/mSigma;
      mLLK=mLLK2*(1+mLLK2);
      #else
      if(calcDeriv) mDerivLLKCoeff=2*mLLK/mSigma;
      mLLK *= mLLK;
      #endif
      VFN_DEBUG_EXIT("MolBondAngle::GetLogLikelihood():",2)
      return mLLK;
   }
   mLLK=angle-(mAngle0-mDelta);
   if(mLLK<0)
   {
      mLLK/=mSigma;
      #ifdef RESTRAINT_X2_X4_X6
      const float mLLK2=mLLK*mLLK;
      //if(calcDeriv) mDerivLLKCoeff=(2*mLLK+4*mLLK2+6*mLLK2*mLLK2)/mSigma;
      //mLLK=mLLK2*(1+mLLK2+mLLK2*mLLK2);
      if(calcDeriv) mDerivLLKCoeff=(2*mLLK+4*mLLK2)/mSigma;
      mLLK=mLLK2*(1+mLLK2);
      #else
      if(calcDeriv) mDerivLLKCoeff=2*mLLK/mSigma;
      mLLK *= mLLK;
      #endif
      VFN_DEBUG_EXIT("MolBondAngle::GetLogLikelihood():",2)
      return mLLK;
   }
   VFN_DEBUG_EXIT("MolBondAngle::GetLogLikelihood():",2)
   if(calcDeriv) mDerivLLKCoeff=0;
   mLLK=0;
   return mLLK;
}

REAL MolBondAngle::GetDeriv(const std::map<const MolAtom*,XYZ> &m,const bool llk)const
{
   //TAU_PROFILE("MolBondAngle::GetDeriv()","REAL (mak,bool)",TAU_DEFAULT);
   REAL d=0;
   map<const MolAtom*,XYZ>::const_iterator pos;
   pos=m.find(mvpAtom[0]);
   if(pos!=m.end())
      d+= pos->second.x*mDerivAtom1.x
         +pos->second.y*mDerivAtom1.y
         +pos->second.z*mDerivAtom1.z;
   pos=m.find(mvpAtom[1]);
   if(pos!=m.end())
      d+= pos->second.x*mDerivAtom2.x
         +pos->second.y*mDerivAtom2.y
         +pos->second.z*mDerivAtom2.z;
   pos=m.find(mvpAtom[2]);
   if(pos!=m.end())
      d+= pos->second.x*mDerivAtom3.x
         +pos->second.y*mDerivAtom3.y
         +pos->second.z*mDerivAtom3.z;
   if(llk) return mDerivLLKCoeff*d;
   return d;
}

void MolBondAngle::CalcGradient(std::map<MolAtom*,XYZ> &m)const
{
   this->GetLogLikelihood(true,true);
   map<MolAtom*,XYZ>::iterator pos;
   pos=m.find(mvpAtom[0]);
   if(pos!=m.end())
   {
      pos->second.x+=mDerivLLKCoeff*mDerivAtom1.x;
      pos->second.y+=mDerivLLKCoeff*mDerivAtom1.y;
      pos->second.z+=mDerivLLKCoeff*mDerivAtom1.z;
   }
   pos=m.find(mvpAtom[1]);
   if(pos!=m.end())
   {
      pos->second.x+=mDerivLLKCoeff*mDerivAtom2.x;
      pos->second.y+=mDerivLLKCoeff*mDerivAtom2.y;
      pos->second.z+=mDerivLLKCoeff*mDerivAtom2.z;
   }
   pos=m.find(mvpAtom[2]);
   if(pos!=m.end())
   {
      pos->second.x+=mDerivLLKCoeff*mDerivAtom3.x;
      pos->second.y+=mDerivLLKCoeff*mDerivAtom3.y;
      pos->second.z+=mDerivLLKCoeff*mDerivAtom3.z;
   }
   #if 0
   // Display gradient - for tests
   cout<<this->GetName()<<" :LLK"<<endl;
   for(map<MolAtom*,XYZ>::const_iterator pos=m.begin();pos!=m.end();++pos)
   {
      char buf[100];
      sprintf(buf,"%10s Grad LLK= (%8.3f %8.3f %8.3f)",
            pos->first->GetName().c_str(),pos->second.x,pos->second.y,pos->second.z);
      cout<<buf<<endl;
   }
   #endif
}

const MolAtom& MolBondAngle::GetAtom1()const{return *(mvpAtom[0]);}
const MolAtom& MolBondAngle::GetAtom2()const{return *(mvpAtom[1]);}
const MolAtom& MolBondAngle::GetAtom3()const{return *(mvpAtom[2]);}
MolAtom& MolBondAngle::GetAtom1(){return *(mvpAtom[0]);}
MolAtom& MolBondAngle::GetAtom2(){return *(mvpAtom[1]);}
MolAtom& MolBondAngle::GetAtom3(){return *(mvpAtom[2]);}
void MolBondAngle::SetAtom1(MolAtom& at){mvpAtom[0]=&at;}
void MolBondAngle::SetAtom2(MolAtom& at){mvpAtom[1]=&at;}
void MolBondAngle::SetAtom3(MolAtom& at){mvpAtom[2]=&at;}
//MolAtom& MolBondAngle::GetAtom1(){return *(mvpAtom[0]);}
//MolAtom& MolBondAngle::GetAtom2(){return *(mvpAtom[1]);}
//MolAtom& MolBondAngle::GetAtom3(){return *(mvpAtom[2]);}
#ifdef __WX__CRYST__
WXCrystObjBasic* MolBondAngle::WXCreate(wxWindow* parent)
{
   VFN_DEBUG_ENTRY("MolBondAngle::WXCreate()",5)
   mpWXCrystObj=new WXMolBondAngle(parent,this);
   VFN_DEBUG_EXIT("MolBondAngle::WXCreate()",5)
   return mpWXCrystObj;
}
WXCrystObjBasic* MolBondAngle::WXGet(){return mpWXCrystObj;}
void MolBondAngle::WXDelete(){if(0!=mpWXCrystObj) delete mpWXCrystObj;mpWXCrystObj=0;}
void MolBondAngle::WXNotifyDelete(){mpWXCrystObj=0;}
#endif
//######################################################################
//
//      MolDihedralAngle
//
//######################################################################
MolDihedralAngle::MolDihedralAngle(MolAtom &atom1, MolAtom &atom2,
                                   MolAtom &atom3, MolAtom &atom4,
                                   const REAL angle, const REAL sigma, const REAL delta,
                                   Molecule &parent):
mAngle0(angle),mDelta(delta),mSigma(sigma),mpMol(&parent)
#ifdef __WX__CRYST__
,mpWXCrystObj(0)
#endif
{
   VFN_DEBUG_ENTRY("MolDihedralAngle::MolDihedralAngle()",5)
   mvpAtom.push_back(&atom1);
   mvpAtom.push_back(&atom2);
   mvpAtom.push_back(&atom3);
   mvpAtom.push_back(&atom4);
   // We want the angle in [-pi;pi]
   mAngle0=fmod((REAL)mAngle0,(REAL)(2*M_PI));
   if(mAngle0<(-M_PI)) mAngle0+=2*M_PI;
   if(mAngle0>M_PI) mAngle0-=2*M_PI;
   VFN_DEBUG_EXIT("MolDihedralAngle::MolDihedralAngle()",5)
}

MolDihedralAngle::~MolDihedralAngle()
{
#ifdef __WX__CRYST__
this->WXDelete();
#endif
}

const Molecule& MolDihedralAngle::GetMolecule()const{return *mpMol;}
      Molecule& MolDihedralAngle::GetMolecule()     {return *mpMol;}

string MolDihedralAngle::GetName()const
{
   return this->GetAtom1().GetName()+"-"
         +this->GetAtom2().GetName()+"-"
         +this->GetAtom3().GetName()+"-"
         +this->GetAtom4().GetName();
}

void MolDihedralAngle::XMLOutput(ostream &os,int indent)const
{
   VFN_DEBUG_ENTRY("MolDihedralAngle::XMLOutput()",4)
   for(int i=0;i<indent;i++) os << "  " ;
   XMLCrystTag tag("DihedralAngle",false,true);
   tag.AddAttribute("Atom1",this->GetAtom1().GetName());
   tag.AddAttribute("Atom2",this->GetAtom2().GetName());
   tag.AddAttribute("Atom3",this->GetAtom3().GetName());
   tag.AddAttribute("Atom4",this->GetAtom4().GetName());
   {
      stringstream ss;
      ss.precision(os.precision());
      ss <<mAngle0*RAD2DEG;
      tag.AddAttribute("Angle",ss.str());
   }
   {
      stringstream ss;
      ss.precision(os.precision());
      ss <<mDelta*RAD2DEG;
      tag.AddAttribute("Delta",ss.str());
   }
   {
      stringstream ss;
      ss.precision(os.precision());
      ss <<mSigma*RAD2DEG;
      tag.AddAttribute("Sigma",ss.str());
   }
   os <<tag<<endl;
   VFN_DEBUG_EXIT("MolDihedralAngle::XMLOutput()",4)
}

void MolDihedralAngle::XMLInput(istream &is,const XMLCrystTag &tag)
{
   VFN_DEBUG_ENTRY("MolDihedralAngle::XMLInput():",5)
   mvpAtom.resize(4);
   for(unsigned int i=0;i<tag.GetNbAttribute();i++)
   {
      if("Atom1"==tag.GetAttributeName(i))
      {
         mvpAtom[0]=&(mpMol->GetAtom(tag.GetAttributeValue(i)));
      }
      if("Atom2"==tag.GetAttributeName(i))
      {
         mvpAtom[1]=&(mpMol->GetAtom(tag.GetAttributeValue(i)));
      }
      if("Atom3"==tag.GetAttributeName(i))
      {
         mvpAtom[2]=&(mpMol->GetAtom(tag.GetAttributeValue(i)));
      }
      if("Atom4"==tag.GetAttributeName(i))
      {
         mvpAtom[3]=&(mpMol->GetAtom(tag.GetAttributeValue(i)));
      }
      if("Angle"==tag.GetAttributeName(i))
      {
         stringstream ss(tag.GetAttributeValue(i));
         ss >>mAngle0;
         mAngle0*=DEG2RAD;
      }
      if("Delta"==tag.GetAttributeName(i))
      {
         stringstream ss(tag.GetAttributeValue(i));
         ss >>mDelta;
         mDelta*=DEG2RAD;
      }
      if("Sigma"==tag.GetAttributeName(i))
      {
         stringstream ss(tag.GetAttributeValue(i));
         ss >>mSigma;
         mSigma*=DEG2RAD;
      }
   }
   VFN_DEBUG_EXIT("MolDihedralAngle::XMLInput():",5)
}

REAL MolDihedralAngle::GetAngle()const
{
   //Get the angle [2pi] closest to the restraint
   const REAL angle=GetDihedralAngle(this->GetAtom1(),this->GetAtom2(),this->GetAtom3(),this->GetAtom4());
   if((angle-mAngle0)>M_PI) return angle-2*M_PI;
   else if((angle-mAngle0)<(-M_PI)) return angle+2*M_PI;

   return angle;
}

REAL& MolDihedralAngle::Angle0(){return mAngle0;}
REAL& MolDihedralAngle::AngleDelta(){return mDelta;}
REAL& MolDihedralAngle::AngleSigma(){return mSigma;}

REAL MolDihedralAngle::GetAngle0()const{return mAngle0;}
REAL MolDihedralAngle::GetAngleDelta()const{return mDelta;}
REAL MolDihedralAngle::GetAngleSigma()const{return mSigma;}

void MolDihedralAngle::SetAngle0(const REAL angle){mAngle0=angle;}
void MolDihedralAngle::SetAngleDelta(const REAL delta){mDelta=delta;}
void MolDihedralAngle::SetAngleSigma(const REAL sigma){mSigma=sigma;}

REAL MolDihedralAngle::GetLogLikelihood()const{return this->GetLogLikelihood(false,true);}

REAL MolDihedralAngle::GetLogLikelihood(const bool calcDeriv, const bool recalc)const
{
   if(!recalc) return mLLK;
   VFN_DEBUG_ENTRY("MolDihedralAngle::GetLogLikelihood():",2)
   //TAU_PROFILE("MolDihedralAngle::GetLogLikelihood()","REAL (bool,bool)",TAU_DEFAULT);
   const REAL x21=this->GetAtom1().GetX()-this->GetAtom2().GetX();
   const REAL y21=this->GetAtom1().GetY()-this->GetAtom2().GetY();
   const REAL z21=this->GetAtom1().GetZ()-this->GetAtom2().GetZ();
   
   const REAL x34=this->GetAtom4().GetX()-this->GetAtom3().GetX();
   const REAL y34=this->GetAtom4().GetY()-this->GetAtom3().GetY();
   const REAL z34=this->GetAtom4().GetZ()-this->GetAtom3().GetZ();
   
   const REAL x23=this->GetAtom3().GetX()-this->GetAtom2().GetX();
   const REAL y23=this->GetAtom3().GetY()-this->GetAtom2().GetY();
   const REAL z23=this->GetAtom3().GetZ()-this->GetAtom2().GetZ();
   
   // v21 x v23
   const REAL x123= y21*z23-z21*y23;
   const REAL y123= z21*x23-x21*z23;
   const REAL z123= x21*y23-y21*x23;
   
   // v32 x v34 (= -v23 x v34)
   const REAL x234= -(y23*z34-z23*y34);
   const REAL y234= -(z23*x34-x23*z34);
   const REAL z234= -(x23*y34-y23*x34);

   const REAL n123= sqrt(x123*x123+y123*y123+z123*z123+1e-7);
   const REAL n234= sqrt(x234*x234+y234*y234+z234*z234+1e-6);
   
   const REAL p=x123*x234+y123*y234+z123*z234;
   const REAL a0=p/(n123*n234+1e-10);
   REAL angle;
   if(a0>= 1) angle=0;
   else 
   {
      if(a0<=-1) angle=M_PI;
      else angle=acos(a0);
   }
   REAL sgn=1.0;
   if((x21*x34 + y21*y34 + z21*z34)<0) {angle=-angle;sgn=-1;}
   
   
   if(calcDeriv)
   {
      const REAL s=sgn/(sqrt(1-a0*a0+1e-6));
      const REAL s0=-s/(n123*n234+1e-10);
      const REAL s1= s*p/(n234*n123*n123*n123+1e-10);
      const REAL s3= s*p/(n123*n234*n234*n234+1e-10);
      mDerivAtom1.x=s0*(-z23*y234+y23*z234)+s1*(-z23*y123+y23*z123);
      mDerivAtom1.y=s0*(-x23*z234+z23*x234)+s1*(-x23*z123+z23*x123);
      mDerivAtom1.z=s0*(-y23*x234+x23*y234)+s1*(-y23*x123+x23*y123);

      mDerivAtom4.x=s0*(-z23*y123+y23*z123)+s3*(-z23*y234+y23*z234);
      mDerivAtom4.y=s0*(-x23*z123+z23*x123)+s3*(-x23*z234+z23*x234);
      mDerivAtom4.z=s0*(-y23*x123+x23*y123)+s3*(-y23*x234+x23*y234);

      mDerivAtom2.x=s0*((z23-z21)*y234-y123*z34+(y21-y23)*z234+z123*y34)+s1*(y123*(z23-z21)+z123*(y21-y23))+s3*(-y234*z34+z234*y34);
      mDerivAtom2.y=s0*((x23-x21)*z234-z123*x34+(z21-z23)*x234+x123*z34)+s1*(z123*(x23-x21)+x123*(z21-z23))+s3*(-z234*x34+x234*z34);
      mDerivAtom2.z=s0*((y23-y21)*x234-x123*y34+(x21-x23)*y234+y123*x34)+s1*(x123*(y23-y21)+y123*(x21-x23))+s3*(-x234*y34+y234*x34);

      mDerivAtom3.x=-(mDerivAtom1.x+mDerivAtom2.x+mDerivAtom4.x);
      mDerivAtom3.y=-(mDerivAtom1.y+mDerivAtom2.y+mDerivAtom4.y);
      mDerivAtom3.z=-(mDerivAtom1.z+mDerivAtom2.z+mDerivAtom4.z);
   }
   
   if(mSigma<1e-6)
   {
      if(calcDeriv) mDerivLLKCoeff=0;
      mLLK=0;
      return mLLK;
   }
   mLLK=angle-(mAngle0+mDelta);
   if(mLLK<(-M_PI)) mLLK += 2*M_PI;
   if(mLLK>  M_PI ) mLLK -= 2*M_PI;
   if(mLLK>0)
   {
      mLLK/=mSigma;
      #ifdef RESTRAINT_X2_X4_X6
      const float mLLK2=mLLK*mLLK;
      //if(calcDeriv) mDerivLLKCoeff=(2*mLLK+4*mLLK2+6*mLLK2*mLLK2)/mSigma;
      //mLLK=mLLK2*(1+mLLK2+mLLK2*mLLK2);
      if(calcDeriv) mDerivLLKCoeff=(2*mLLK+4*mLLK2)/mSigma;
      mLLK=mLLK2*(1+mLLK2);
      #else
      if(calcDeriv) mDerivLLKCoeff=2*mLLK/mSigma;
      mLLK *= mLLK;
      #endif
      VFN_DEBUG_EXIT("MolDihedralAngle::GetLogLikelihood():",2)
      return mLLK;
   }
   mLLK=angle-(mAngle0-mDelta);
   if(mLLK<(-M_PI)) mLLK += 2*M_PI;
   if(mLLK>  M_PI ) mLLK -= 2*M_PI;
   if(mLLK<0)
   {
      mLLK/=mSigma;
      #ifdef RESTRAINT_X2_X4_X6
      const float mLLK2=mLLK*mLLK;
      //if(calcDeriv) mDerivLLKCoeff=(2*mLLK+4*mLLK2+6*mLLK2*mLLK2)/mSigma;
      //mLLK=mLLK2*(1+mLLK2+mLLK2*mLLK2);
      if(calcDeriv) mDerivLLKCoeff=(2*mLLK+4*mLLK2)/mSigma;
      mLLK=mLLK2*(1+mLLK2);
      #else
      if(calcDeriv) mDerivLLKCoeff=2*mLLK/mSigma;
      mLLK *= mLLK;
      #endif
      VFN_DEBUG_EXIT("MolDihedralAngle::GetLogLikelihood():",2)
      return mLLK;
   }
   VFN_DEBUG_EXIT("MolDihedralAngle::GetLogLikelihood():",2)
   if(calcDeriv) mDerivLLKCoeff=0;
   mLLK=0;
   return 0;
}

REAL MolDihedralAngle::GetDeriv(const std::map<const MolAtom*,XYZ> &m,const bool llk)const
{
   //TAU_PROFILE("MolDihedralAngle::GetDeriv()","REAL (mak,bool)",TAU_DEFAULT);
   REAL d=0;
   map<const MolAtom*,XYZ>::const_iterator pos;
   pos=m.find(mvpAtom[0]);
   if(pos!=m.end())
      d+= pos->second.x*mDerivAtom1.x
         +pos->second.y*mDerivAtom1.y
         +pos->second.z*mDerivAtom1.z;
   pos=m.find(mvpAtom[1]);
   if(pos!=m.end())
      d+= pos->second.x*mDerivAtom2.x
         +pos->second.y*mDerivAtom2.y
         +pos->second.z*mDerivAtom2.z;
   pos=m.find(mvpAtom[2]);
   if(pos!=m.end())
      d+= pos->second.x*mDerivAtom3.x
         +pos->second.y*mDerivAtom3.y
         +pos->second.z*mDerivAtom3.z;
   pos=m.find(mvpAtom[3]);
   if(pos!=m.end())
      d+= pos->second.x*mDerivAtom4.x
         +pos->second.y*mDerivAtom4.y
         +pos->second.z*mDerivAtom4.z;
   if(llk) return mDerivLLKCoeff*d;
   return d;
}

void MolDihedralAngle::CalcGradient(std::map<MolAtom*,XYZ> &m)const
{
   this->GetLogLikelihood(true,true);
   map<MolAtom*,XYZ>::iterator pos;
   pos=m.find(mvpAtom[0]);
   if(pos!=m.end())
   {
      pos->second.x+=mDerivLLKCoeff*mDerivAtom1.x;
      pos->second.y+=mDerivLLKCoeff*mDerivAtom1.y;
      pos->second.z+=mDerivLLKCoeff*mDerivAtom1.z;
   }
   pos=m.find(mvpAtom[1]);
   if(pos!=m.end())
   {
      pos->second.x+=mDerivLLKCoeff*mDerivAtom2.x;
      pos->second.y+=mDerivLLKCoeff*mDerivAtom2.y;
      pos->second.z+=mDerivLLKCoeff*mDerivAtom2.z;
   }
   pos=m.find(mvpAtom[2]);
   if(pos!=m.end())
   {
      pos->second.x+=mDerivLLKCoeff*mDerivAtom3.x;
      pos->second.y+=mDerivLLKCoeff*mDerivAtom3.y;
      pos->second.z+=mDerivLLKCoeff*mDerivAtom3.z;
   }
   pos=m.find(mvpAtom[3]);
   if(pos!=m.end())
   {
      pos->second.x+=mDerivLLKCoeff*mDerivAtom4.x;
      pos->second.y+=mDerivLLKCoeff*mDerivAtom4.y;
      pos->second.z+=mDerivLLKCoeff*mDerivAtom4.z;
   }
   #if 0
   // Display gradient - for tests
   cout<<this->GetName()<<" :LLK"<<endl;
   for(map<MolAtom*,XYZ>::const_iterator pos=m.begin();pos!=m.end();++pos)
   {
      char buf[100];
      sprintf(buf,"%10s Grad LLK= (%8.3f %8.3f %8.3f)",
            pos->first->GetName().c_str(),pos->second.x,pos->second.y,pos->second.z);
      cout<<buf<<endl;
   }
   #endif
}

const MolAtom& MolDihedralAngle::GetAtom1()const{return *(mvpAtom[0]);}
const MolAtom& MolDihedralAngle::GetAtom2()const{return *(mvpAtom[1]);}
const MolAtom& MolDihedralAngle::GetAtom3()const{return *(mvpAtom[2]);}
const MolAtom& MolDihedralAngle::GetAtom4()const{return *(mvpAtom[3]);}
void MolDihedralAngle::SetAtom1(MolAtom& at){mvpAtom[0]=&at;}
void MolDihedralAngle::SetAtom2(MolAtom& at){mvpAtom[1]=&at;}
void MolDihedralAngle::SetAtom3(MolAtom& at){mvpAtom[2]=&at;}
void MolDihedralAngle::SetAtom4(MolAtom& at){mvpAtom[3]=&at;}
MolAtom& MolDihedralAngle::GetAtom1(){return *(mvpAtom[0]);}
MolAtom& MolDihedralAngle::GetAtom2(){return *(mvpAtom[1]);}
MolAtom& MolDihedralAngle::GetAtom3(){return *(mvpAtom[2]);}
MolAtom& MolDihedralAngle::GetAtom4(){return *(mvpAtom[3]);}
#ifdef __WX__CRYST__
WXCrystObjBasic* MolDihedralAngle::WXCreate(wxWindow* parent)
{
   VFN_DEBUG_ENTRY("MolDihedralAngle::WXCreate()",5)
   mpWXCrystObj=new WXMolDihedralAngle(parent,this);
   VFN_DEBUG_EXIT("MolDihedralAngle::WXCreate()",5)
   return mpWXCrystObj;
}
WXCrystObjBasic* MolDihedralAngle::WXGet(){return mpWXCrystObj;}
void MolDihedralAngle::WXDelete(){if(0!=mpWXCrystObj) delete mpWXCrystObj;mpWXCrystObj=0;}
void MolDihedralAngle::WXNotifyDelete(){mpWXCrystObj=0;}
#endif
//######################################################################
//
//      RigidGroup
//
//######################################################################
string RigidGroup::GetName()const
{
   set<MolAtom *>::const_iterator at=this->begin();
   string name=(*at++)->GetName();
   for(;at!=this->end();++at) name+=", "+(*at)->GetName();
   return name;
}
//######################################################################
//
//      MolRing
//
//######################################################################
MolRing::MolRing()
{}

const std::list<MolAtom*>& MolRing::GetAtomList()const
{return mvpAtom;}

std::list<MolAtom*>& MolRing::GetAtomList()
{return mvpAtom;}
//######################################################################
//
//      Quaternion
//
//######################################################################
Quaternion::Quaternion():
mQ0(1),mQ1(0),mQ2(0),mQ3(0),mIsUniQuaternion(true)
{
   VFN_DEBUG_MESSAGE("Quaternion::Quaternion()",5)
}

Quaternion::Quaternion(const REAL q0,
                               const REAL q1,
                               const REAL q2,
                               const REAL q3,
                               bool unit):
mQ0(q0),mQ1(q1),mQ2(q2),mQ3(q3),mIsUniQuaternion(unit)
{
   VFN_DEBUG_MESSAGE("Quaternion::Quaternion()",5)
   if(unit) this->Normalize();
}

Quaternion::~Quaternion()
{
   VFN_DEBUG_MESSAGE("Quaternion::~Quaternion()",5)
}

Quaternion Quaternion::RotationQuaternion(const REAL ang,
                                          const REAL v1,
                                          const REAL v2,
                                          const REAL v3)
{
   VFN_DEBUG_MESSAGE("Quaternion::RotationQuaternion()",4)
   const REAL s=sin(ang/2.)/sqrt(v1*v1+v2*v2+v3*v3+1e-7);
   return Quaternion(cos(ang/2.),s*v1,s*v2,s*v3,
                     true);
}

Quaternion Quaternion::GetConjugate()const
{
   return Quaternion(mQ0,-mQ1,-mQ2,-mQ3);
}
Quaternion Quaternion::operator*(const Quaternion &q)const
{
   // http://www.cs.berkeley.edu/~laura/cs184/quat/quaternion.html
   return Quaternion
      (this->Q0()*q.Q0()-this->Q1()*q.Q1()-this->Q2()*q.Q2()-this->Q3()*q.Q3(),
       this->Q0()*q.Q1()+this->Q1()*q.Q0()+this->Q2()*q.Q3()-this->Q3()*q.Q2(),
       this->Q0()*q.Q2()-this->Q1()*q.Q3()+this->Q2()*q.Q0()+this->Q3()*q.Q1(),
       this->Q0()*q.Q3()+this->Q1()*q.Q2()-this->Q2()*q.Q1()+this->Q3()*q.Q0(),false);
}

void Quaternion::operator*=(const Quaternion &q)
{
   //cout<<"Quaternion::operator*= before:";this->XMLOutput(cout);
   //cout<<"Quaternion::operator*= by    :";q.XMLOutput(cout);
   const REAL q1=this->Q0()*q.Q1()+this->Q1()*q.Q0()+this->Q2()*q.Q3()-this->Q3()*q.Q2();
   const REAL q2=this->Q0()*q.Q2()+this->Q2()*q.Q0()-this->Q1()*q.Q3()+this->Q3()*q.Q1();
   const REAL q3=this->Q0()*q.Q3()+this->Q3()*q.Q0()+this->Q1()*q.Q2()-this->Q2()*q.Q1();
   this->Q0()=   this->Q0()*q.Q0()-this->Q1()*q.Q1()-this->Q2()*q.Q2()-this->Q3()*q.Q3();
   this->Q1()=q1;
   this->Q2()=q2;
   this->Q3()=q3;
   this->Normalize();
   //cout<<"Quaternion::operator*= after :";this->XMLOutput(cout);
}

void Quaternion::XMLOutput(ostream &os,int indent)const
{
   VFN_DEBUG_ENTRY("Quaternion::XMLOutput()",4)
   for(int i=0;i<indent;i++) os << "  " ;
   XMLCrystTag tag("Quaternion",false,true);
   //#error "which atoms for this bond ?"
   {
      stringstream ss;
      ss.precision(os.precision());
      ss <<mQ0;
      tag.AddAttribute("Q0",ss.str());
   }
   {
      stringstream ss;
      ss.precision(os.precision());
      ss <<mQ1;
      tag.AddAttribute("Q1",ss.str());
   }
   {
      stringstream ss;
      ss.precision(os.precision());
      ss <<mQ2;
      tag.AddAttribute("Q2",ss.str());
   }
   {
      stringstream ss;
      ss.precision(os.precision());
      ss <<mQ3;
      tag.AddAttribute("Q3",ss.str());
   }
   {
      stringstream ss;
      ss.precision(os.precision());
      ss <<mIsUniQuaternion;
      tag.AddAttribute("IsUnitQuaternion",ss.str());
   }
   os <<tag<<endl;
   VFN_DEBUG_EXIT("Quaternion::XMLOutput()",4)
}

void Quaternion::XMLInput(istream &is,const XMLCrystTag &tag)
{
   VFN_DEBUG_ENTRY("Quaternion::XMLInput()",5)
   for(unsigned int i=0;i<tag.GetNbAttribute();i++)
   {
      if("Q0"==tag.GetAttributeName(i))
      {
         stringstream ss(tag.GetAttributeValue(i));
         ss >>mQ0;
      }
      if("Q1"==tag.GetAttributeName(i))
      {
         stringstream ss(tag.GetAttributeValue(i));
         ss >>mQ1;
      }
      if("Q2"==tag.GetAttributeName(i))
      {
         stringstream ss(tag.GetAttributeValue(i));
         ss >>mQ2;
      }
      if("Q3"==tag.GetAttributeName(i))
      {
         stringstream ss(tag.GetAttributeValue(i));
         ss >>mQ3;
      }
      if("IsUnitQuaternion"==tag.GetAttributeName(i))
      {
         stringstream ss(tag.GetAttributeValue(i));
         ss >>mIsUniQuaternion;
      }
   }
   if(mIsUniQuaternion) this->Normalize();
   VFN_DEBUG_EXIT("Quaternion::XMLInput()",5)
}

void Quaternion::RotateVector(REAL &v1,REAL &v2, REAL &v3)const
{
   #if 0
   //#error P should not be a _UNIT_ quaternion...
   Quaternion P(0,v1,v2,v3,false);
   //cout<<"RotQuat:(n="<<this->GetNorm()<<")";this->XMLOutput(cout);
   //cout<<"before :(n="<<P.GetNorm()<<")";P.XMLOutput(cout);
   P= (*this)* P * this->GetConjugate();
   //cout<<"rotated:(n="<<P.GetNorm()<<")";P.XMLOutput(cout);
   v1=P.Q1();
   v2=P.Q2();
   v3=P.Q3();
   #endif
   const REAL p0=-mQ1*v1          - mQ2*v2 - mQ3*v3;
   const REAL p1= mQ0*v1          + mQ2*v3 - mQ3*v2;
   const REAL p2= mQ0*v2          - mQ1*v3 + mQ3*v1;
   const REAL p3= mQ0*v3          + mQ1*v2 - mQ2*v1;
   
   v1           =-p0*mQ1 + p1*mQ0 - p2*mQ3 + p3*mQ2;
   v2           =-p0*mQ2 + p2*mQ0 + p1*mQ3 - p3*mQ1;
   v3           =-p0*mQ3 + p3*mQ0 - p1*mQ2 + p2*mQ1;
}

void Quaternion::Normalize()const
{
   const REAL norm=sqrt( this->Q0()*this->Q0()
                        +this->Q1()*this->Q1()
                        +this->Q2()*this->Q2()
                        +this->Q3()*this->Q3());
   mQ0 /= norm;
   mQ1 /= norm;
   mQ2 /= norm;
   mQ3 /= norm;
}
REAL Quaternion::GetNorm()const
{return sqrt( this->Q0()*this->Q0()
             +this->Q1()*this->Q1()
             +this->Q2()*this->Q2()
             +this->Q3()*this->Q3());}

const REAL& Quaternion::Q0()const{return mQ0;}
const REAL& Quaternion::Q1()const{return mQ1;}
const REAL& Quaternion::Q2()const{return mQ2;}
const REAL& Quaternion::Q3()const{return mQ3;}
REAL& Quaternion::Q0(){return mQ0;}
REAL& Quaternion::Q1(){return mQ1;}
REAL& Quaternion::Q2(){return mQ2;}
REAL& Quaternion::Q3(){return mQ3;}
//######################################################################
//
//      Molecule Stretch Modes
//
//######################################################################
StretchMode::~StretchMode(){}

StretchModeBondLength::StretchModeBondLength(MolAtom &at0,MolAtom &at1,
                                             const MolBond *pBond):
mpAtom0(&at0),mpAtom1(&at1),mpBond(pBond)
{
   mBaseAmplitude=0.1;
   mpMol = &(at1.GetMolecule());
}

StretchModeBondLength::~StretchModeBondLength(){}

void StretchModeBondLength::CalcDeriv(const bool derivllk)const
{
   //TAU_PROFILE("StretchModeBondLength::CalcDeriv()","void ()",TAU_DEFAULT);
   // Derivative of the atom positions
      //mDerivXYZ.clear();
      REAL dx=mpAtom1->GetX()-mpAtom0->GetX();
      REAL dy=mpAtom1->GetY()-mpAtom0->GetY();
      REAL dz=mpAtom1->GetZ()-mpAtom0->GetZ();
      const REAL n=sqrt(dx*dx+dy*dy+dz*dz+1e-7);
      if(n<1e-6) return;//:KLUDGE: ?
      dx/=n;
      dy/=n;
      dz/=n;
      for(set<MolAtom *>::const_iterator pos=mvTranslatedAtomList.begin();
          pos!=mvTranslatedAtomList.end();++pos)
      {
         XYZ *const p=&(mDerivXYZ[*pos]);
         p->x=dx;
         p->y=dy;
         p->z=dz;
      }
   // Derivative of the LLK
   if(derivllk)
   {
      mLLKDeriv=0;
      for(map<const MolBond*,REAL>::const_iterator pos=this->mvpBrokenBond.begin();
          pos!=this->mvpBrokenBond.end();++pos) mLLKDeriv += pos->first->GetDeriv(mDerivXYZ,true);
      for(map<const MolBondAngle*,REAL>::const_iterator pos=this->mvpBrokenBondAngle.begin();
          pos!=this->mvpBrokenBondAngle.end();++pos) mLLKDeriv += pos->first->GetDeriv(mDerivXYZ,true);
      for(map<const MolDihedralAngle*,REAL>::const_iterator pos=this->mvpBrokenDihedralAngle.begin();
          pos!=this->mvpBrokenDihedralAngle.end();++pos) mLLKDeriv += pos->first->GetDeriv(mDerivXYZ,true);
   }
}

void StretchModeBondLength::Print(ostream &os,bool full)const
{
   cout<<mpAtom0->GetName()<<"-"<<mpAtom1->GetName();
   if(full)
   {
      cout<<", Translated Atoms:";
      for(set<MolAtom*>::const_iterator atom=mvTranslatedAtomList.begin();
          atom!=mvTranslatedAtomList.end();++atom)
      {
         cout<<(*atom)->GetName()<<",";
      }
   }
}

void StretchModeBondLength::Stretch(const REAL amplitude, 
                                    const bool keepCenter)
{
   REAL dx=mpAtom1->GetX()-mpAtom0->GetX();
   REAL dy=mpAtom1->GetY()-mpAtom0->GetY();
   REAL dz=mpAtom1->GetZ()-mpAtom0->GetZ();
   const REAL l=sqrt(dx*dx+dy*dy+dz*dz+1e-7);
   if(l<1e-6) return;// :KLUDGE:
   const REAL change=amplitude/l;
   dx*=change;
   dy*=change;
   dz*=change;
   mpMol->TranslateAtomGroup(mvTranslatedAtomList,dx,dy,dz,keepCenter);
}

void StretchModeBondLength::RandomStretch(const REAL amplitude, 
                                          const bool keepCenter)
{
   mpMol->BondLengthRandomChange(*this,amplitude,keepCenter);
}

StretchModeBondAngle::StretchModeBondAngle(MolAtom &at0,MolAtom &at1,MolAtom &at2,
                                           const MolBondAngle *pBondAngle):
mpAtom0(&at0),mpAtom1(&at1),mpAtom2(&at2),mpBondAngle(pBondAngle)
{
   mBaseAmplitude=M_PI*0.02;
   mpMol = &(at1.GetMolecule());
}

StretchModeBondAngle::~StretchModeBondAngle(){}

void StretchModeBondAngle::CalcDeriv(const bool derivllk)const
{
   //TAU_PROFILE("StretchModeBondAngle::CalcDeriv()","void ()",TAU_DEFAULT);
   // Derivative of the atomic positions
      const REAL x1=mpAtom1->GetX(),
                 y1=mpAtom1->GetY(),
                 z1=mpAtom1->GetZ();

      const REAL dx10=mpAtom0->GetX()-x1,
                 dy10=mpAtom0->GetY()-y1,
                 dz10=mpAtom0->GetZ()-z1,
                 dx12=mpAtom2->GetX()-x1,
                 dy12=mpAtom2->GetY()-y1,
                 dz12=mpAtom2->GetZ()-z1;

      REAL vx=dy10*dz12-dz10*dy12,
           vy=dz10*dx12-dx10*dz12,
           vz=dx10*dy12-dy10*dx12;
      //const REAL n=sqrt((dx10*dx10+dy10*dy10+dz10*dz10)*(dx12*dx12+dy12*dy12+dz12*dz12))+1e-10;
      const REAL n=sqrt(vx*vx+vy*vy+vz*vz+1e-10);
      vx/=n;
      vy/=n;
      vz/=n;

      if(n<1e-6)
      {
         mDerivXYZ.clear();
         return;//:KLUDGE: ?
      }

      for(set<MolAtom *>::const_iterator pos=mvRotatedAtomList.begin();
          pos!=mvRotatedAtomList.end();++pos)
      {
         XYZ *const p=&(mDerivXYZ[*pos]);
         const REAL x=(*pos)->GetX()-x1,
                    y=(*pos)->GetY()-y1,
                    z=(*pos)->GetZ()-z1;
         p->x=z*vy-y*vz;
         p->y=x*vz-z*vx;
         p->z=y*vx-x*vy;
      }
   // Derivative of the LLK
   if(derivllk)
   {
      mLLKDeriv=0;
      for(map<const MolBond*,REAL>::const_iterator pos=this->mvpBrokenBond.begin();
          pos!=this->mvpBrokenBond.end();++pos) mLLKDeriv += pos->first->GetDeriv(mDerivXYZ,true);
      for(map<const MolBondAngle*,REAL>::const_iterator pos=this->mvpBrokenBondAngle.begin();
          pos!=this->mvpBrokenBondAngle.end();++pos) mLLKDeriv += pos->first->GetDeriv(mDerivXYZ,true);
      for(map<const MolDihedralAngle*,REAL>::const_iterator pos=this->mvpBrokenDihedralAngle.begin();
          pos!=this->mvpBrokenDihedralAngle.end();++pos) mLLKDeriv += pos->first->GetDeriv(mDerivXYZ,true);
   }
}

void StretchModeBondAngle::Print(ostream &os,bool full)const
{
   cout<<mpAtom0->GetName()<<"-"<<mpAtom1->GetName()<<"-"<<mpAtom2->GetName();
   if(full)
   {
      cout<<", Rotated Atoms:";
      for(set<MolAtom*>::const_iterator atom=mvRotatedAtomList.begin();
          atom!=mvRotatedAtomList.end();++atom)
      {
         cout<<(*atom)->GetName()<<",";
      }
   }
}

void StretchModeBondAngle::Stretch(const REAL amplitude, 
                                   const bool keepCenter)
{
   REAL dx10=mpAtom0->GetX()-mpAtom1->GetX();
   REAL dy10=mpAtom0->GetY()-mpAtom1->GetY();
   REAL dz10=mpAtom0->GetZ()-mpAtom1->GetZ();
   REAL dx12=mpAtom2->GetX()-mpAtom1->GetX();
   REAL dy12=mpAtom2->GetY()-mpAtom1->GetY();
   REAL dz12=mpAtom2->GetZ()-mpAtom1->GetZ();
   
   const REAL vx=dy10*dz12-dz10*dy12;
   const REAL vy=dz10*dx12-dx10*dz12;
   const REAL vz=dx10*dy12-dy10*dx12;
   mpMol->RotateAtomGroup(*mpAtom1,vx,vy,vz,mvRotatedAtomList,amplitude,keepCenter);
}

void StretchModeBondAngle::RandomStretch(const REAL amplitude, 
                                         const bool keepCenter)
{
   mpMol->BondAngleRandomChange(*this,amplitude,keepCenter);
}

StretchModeTorsion::StretchModeTorsion(MolAtom &at1,MolAtom &at2,
                                       const MolDihedralAngle *pAngle):
mpAtom1(&at1),mpAtom2(&at2),mpDihedralAngle(pAngle)
{
   mBaseAmplitude=M_PI*0.02;
   mpMol = &(at1.GetMolecule());
}

StretchModeTorsion::~StretchModeTorsion(){}

void StretchModeTorsion::CalcDeriv(const bool derivllk)const
{
   //TAU_PROFILE("StretchModeTorsion::CalcDeriv()","void ()",TAU_DEFAULT);
   // Derivative of the LLK
   //mDerivXYZ.clear();
      const REAL x1=mpAtom1->GetX(),
                 y1=mpAtom1->GetY(),
                 z1=mpAtom1->GetZ();

      REAL vx=mpAtom2->GetX()-x1,
           vy=mpAtom2->GetY()-y1,
           vz=mpAtom2->GetZ()-z1;

      const REAL n=sqrt(vx*vx+vy*vy+vz*vz+1e-10);
      vx/=n;
      vy/=n;
      vz/=n;

      for(set<MolAtom *>::const_iterator pos=mvRotatedAtomList.begin();
          pos!=mvRotatedAtomList.end();++pos)
      {
         XYZ *const p=&(mDerivXYZ[*pos]);
         const REAL x=(*pos)->GetX()-x1,
                    y=(*pos)->GetY()-y1,
                    z=(*pos)->GetZ()-z1;
         p->x=z*vy-y*vz;
         p->y=x*vz-z*vx;
         p->z=y*vx-x*vy;
      }
   // Derivative of the LLK
   if(derivllk)
   {
      mLLKDeriv=0;
      for(map<const MolBond*,REAL>::const_iterator pos=this->mvpBrokenBond.begin();
          pos!=this->mvpBrokenBond.end();++pos) mLLKDeriv += pos->first->GetDeriv(mDerivXYZ,true);
      for(map<const MolBondAngle*,REAL>::const_iterator pos=this->mvpBrokenBondAngle.begin();
          pos!=this->mvpBrokenBondAngle.end();++pos) mLLKDeriv += pos->first->GetDeriv(mDerivXYZ,true);
      for(map<const MolDihedralAngle*,REAL>::const_iterator pos=this->mvpBrokenDihedralAngle.begin();
          pos!=this->mvpBrokenDihedralAngle.end();++pos) mLLKDeriv += pos->first->GetDeriv(mDerivXYZ,true);
   }
}

void StretchModeTorsion::Print(ostream &os,bool full)const
{
   cout<<mpAtom1->GetName()<<"-"<<mpAtom2->GetName();
   if(full)
   {
      cout<<", Rotated Atoms:";
      for(set<MolAtom*>::const_iterator atom=mvRotatedAtomList.begin();
          atom!=mvRotatedAtomList.end();++atom)
      {
         cout<<(*atom)->GetName()<<",";
      }
   }
}

void StretchModeTorsion::Stretch(const REAL amplitude, const bool keepCenter)
{
   mpMol->RotateAtomGroup(*mpAtom1,*mpAtom2,mvRotatedAtomList,amplitude,keepCenter);
}

void StretchModeTorsion::RandomStretch(const REAL amplitude, 
                                       const bool keepCenter)
{
   mpMol->DihedralAngleRandomChange(*this,amplitude,keepCenter);
}


//######################################################################
//
//      StretchModeTwist
//
//######################################################################
StretchModeTwist::StretchModeTwist(MolAtom &at1,MolAtom &at2):
mpAtom1(&at1),mpAtom2(&at2)
{
   mBaseAmplitude=M_PI*0.02;
   mpMol = &(at1.GetMolecule());
}

StretchModeTwist::~StretchModeTwist(){}

void StretchModeTwist::CalcDeriv(const bool derivllk)const
{// Identical to StretchModeTorsion::CalcDeriv()
   // Derivative of the LLK
   //mDerivXYZ.clear();
      const REAL x1=mpAtom1->GetX(),
                 y1=mpAtom1->GetY(),
                 z1=mpAtom1->GetZ();

      REAL vx=mpAtom2->GetX()-x1,
           vy=mpAtom2->GetY()-y1,
           vz=mpAtom2->GetZ()-z1;

      const REAL n=sqrt(vx*vx+vy*vy+vz*vz+1e-10);
      vx/=n;
      vy/=n;
      vz/=n;

      for(set<MolAtom *>::const_iterator pos=mvRotatedAtomList.begin();
          pos!=mvRotatedAtomList.end();++pos)
      {
         XYZ *const p=&(mDerivXYZ[*pos]);
         const REAL x=(*pos)->GetX()-x1,
                    y=(*pos)->GetY()-y1,
                    z=(*pos)->GetZ()-z1;
         p->x=z*vy-y*vz;
         p->y=x*vz-z*vx;
         p->z=y*vx-x*vy;
      }
   // Derivative of the LLK
   if(derivllk)
   {
      mLLKDeriv=0;
      for(map<const MolBond*,REAL>::const_iterator pos=this->mvpBrokenBond.begin();
          pos!=this->mvpBrokenBond.end();++pos) mLLKDeriv += pos->first->GetDeriv(mDerivXYZ,true);
      for(map<const MolBondAngle*,REAL>::const_iterator pos=this->mvpBrokenBondAngle.begin();
          pos!=this->mvpBrokenBondAngle.end();++pos) mLLKDeriv += pos->first->GetDeriv(mDerivXYZ,true);
      for(map<const MolDihedralAngle*,REAL>::const_iterator pos=this->mvpBrokenDihedralAngle.begin();
          pos!=this->mvpBrokenDihedralAngle.end();++pos) mLLKDeriv += pos->first->GetDeriv(mDerivXYZ,true);
   }
}

void StretchModeTwist::Print(ostream &os,bool full)const
{
   os<<mpAtom1->GetName()<<"/"<<mpAtom2->GetName()<<"-"<<mpAtom2->GetName();
   if(full)
   {
      os<<", Rotated Atoms:";
      for(set<MolAtom*>::const_iterator atom=mvRotatedAtomList.begin();
          atom!=mvRotatedAtomList.end();++atom)
      {
         os<<(*atom)->GetName()<<",";
      }
   }
}

void StretchModeTwist::Stretch(const REAL amplitude, const bool keepCenter)
{
   mpMol->RotateAtomGroup(*mpAtom1,*mpAtom2,mvRotatedAtomList,amplitude,keepCenter);
}

void StretchModeTwist::RandomStretch(const REAL amplitude, 
                                     const bool keepCenter)
{
   const REAL dx=mpAtom2->GetX()-mpAtom1->GetX();
   const REAL dy=mpAtom2->GetY()-mpAtom1->GetY();
   const REAL dz=mpAtom2->GetZ()-mpAtom1->GetZ();
   if((abs(dx)+abs(dy)+abs(dz))<1e-6) return;// :KLUDGE:
   const REAL change=(REAL)(2.*rand()-RAND_MAX)/(REAL)RAND_MAX*mBaseAmplitude*amplitude;
   mpMol->RotateAtomGroup(*mpAtom1,*mpAtom2,mvRotatedAtomList,change,keepCenter);
}

//######################################################################
//
//      MDAtomGroup
//
//######################################################################
MDAtomGroup::MDAtomGroup(){};

MDAtomGroup::MDAtomGroup(set<MolAtom*> &vat,
                         set<MolBond*> &vb,
                         set<MolBondAngle*> &va,
                         set<MolDihedralAngle*> &vd):
mvpAtom(vat)
{
   // Use vector instead of sets for MolecularDynamicsEvolve & general
   // storage in molecule compatibility
   for(set<MolBond*>::iterator pos=vb.begin();pos!=vb.end();++pos)
      mvpBond.push_back(*pos);
   for(set<MolBondAngle*>::iterator pos=va.begin();pos!=va.end();++pos)
      mvpBondAngle.push_back(*pos);
   for(set<MolDihedralAngle*>::iterator pos=vd.begin();pos!=vd.end();++pos)
      mvpDihedralAngle.push_back(*pos);
}

void MDAtomGroup::Print(ostream &os,bool full)const
{
   if(full) os<<"MDAtomGroup: ";
   for(set<MolAtom*>::const_iterator pos=mvpAtom.begin();pos!=mvpAtom.end();++pos)
      os<<(*pos)->GetName()<<" ";
   if(full)
   {
      os<<endl<<"  Involving bond restraints:";
      for(vector<MolBond*>::const_iterator pos=mvpBond.begin();pos!=mvpBond.end();++pos)
         os<<(*pos)->GetName()<<"  ";
      os<<endl<<"  Involving bond angle restraints:";
      for(vector<MolBondAngle*>::const_iterator pos=mvpBondAngle.begin();pos!=mvpBondAngle.end();++pos)
         os<<(*pos)->GetName()<<"  ";
      os<<endl<<"  Involving dihedral angle restraints:";
      for(vector<MolDihedralAngle*>::const_iterator pos=mvpDihedralAngle.begin();pos!=mvpDihedralAngle.end();++pos)
         os<<(*pos)->GetName()<<"  ";
      os<<endl;
   }
}

//######################################################################
//
//      Molecule
//
//######################################################################
Molecule::Molecule(Crystal &cryst, const string &name):
mBaseRotationAmplitude(M_PI*0.02),mIsSelfOptimizing(false),mpCenterAtom(0),
mMDMoveFreq(0.0),mMDMoveEnergy(40.),mDeleteSubObjInDestructor(1),mLogLikelihoodScale(1.0)
{
   VFN_DEBUG_MESSAGE("Molecule::Molecule()",5)
   this->SetName(name);
   mpCryst=&cryst;
   {
      RefinablePar tmp(this->GetName()+"_x",&mXYZ(0),0.,1.,
                        gpRefParTypeScattTranslX,
                        REFPAR_DERIV_STEP_ABSOLUTE,false,false,true,true,1.,1.);
      tmp.AssignClock(mClockScatterer);
      tmp.SetDerivStep(1e-5);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp(this->GetName()+"_y",&mXYZ(1),0,1,
                        gpRefParTypeScattTranslY,
                        REFPAR_DERIV_STEP_ABSOLUTE,false,false,true,true,1.,1.);
      tmp.AssignClock(mClockScatterer);
      tmp.SetDerivStep(1e-5);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp(this->GetName()+"_z",&mXYZ(2),0,1,
                        gpRefParTypeScattTranslZ,
                        REFPAR_DERIV_STEP_ABSOLUTE,false,false,true,true,1.,1.);
      tmp.AssignClock(mClockScatterer);
      tmp.SetDerivStep(1e-5);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp(this->GetName()+"_Occ",&mOccupancy,0,1,
                        gpRefParTypeScattOccup,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,1.,1.);
      tmp.AssignClock(mClockScatterer);
      tmp.SetDerivStep(1e-5);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp(this->GetName()+"Q0",&(mQuat.Q0()),0,1,
                        gpRefParTypeScattOrient,
                        REFPAR_DERIV_STEP_ABSOLUTE,false,false,true,false,1.,1.);
      tmp.AssignClock(mClockScatterer);
      tmp.SetGlobalOptimStep(0.04);
      tmp.SetDerivStep(1e-4);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp(this->GetName()+"Q1",&(mQuat.Q1()),0,1,
                        gpRefParTypeScattOrient,
                        REFPAR_DERIV_STEP_ABSOLUTE,false,false,true,false,1.,1.);
      tmp.AssignClock(mClockScatterer);
      tmp.SetGlobalOptimStep(0.04);
      tmp.SetDerivStep(1e-4);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp(this->GetName()+"Q2",&(mQuat.Q2()),0,1,
                        gpRefParTypeScattOrient,
                        REFPAR_DERIV_STEP_ABSOLUTE,false,false,true,false,1.,1.);
      tmp.AssignClock(mClockScatterer);
      tmp.SetGlobalOptimStep(0.04);
      tmp.SetDerivStep(1e-4);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp(this->GetName()+"Q3",&(mQuat.Q3()),0,1,
                        gpRefParTypeScattOrient,
                        REFPAR_DERIV_STEP_ABSOLUTE,false,false,true,false,1.,1.);
      tmp.AssignClock(mClockScatterer);
      tmp.SetGlobalOptimStep(0.04);
      tmp.SetDerivStep(1e-4);
      this->AddPar(tmp);
   }
   mLocalParamSet=this->CreateParamSet("saved parameters for local minimization");
   this->InitOptions();
   mClockScatterer.AddChild(mClockAtomList);
   mClockScatterer.AddChild(mClockBondList);
   mClockScatterer.AddChild(mClockBondAngleList);
   mClockScatterer.AddChild(mClockDihedralAngleList);
   mClockScatterer.AddChild(mClockRingList);
   mClockScatterer.AddChild(mClockRigidGroup);
   mClockScatterer.AddChild(mClockAtomPosition);
   mClockScatterer.AddChild(mClockAtomScattPow);
   mClockScatterer.AddChild(mClockOrientation);
}

Molecule::Molecule(const Molecule &old):
mIsSelfOptimizing(false),mpCenterAtom(0),mDeleteSubObjInDestructor(old.mDeleteSubObjInDestructor)
{
   VFN_DEBUG_ENTRY("Molecule::Molecule(old&)",5)
   // a hack, but const-correct
   mpCryst=&(gCrystalRegistry.GetObj(gCrystalRegistry.Find(old.GetCrystal())));
   {
      RefinablePar tmp(this->GetName()+"_x",&mXYZ(0),0.,1.,
                        gpRefParTypeScattTranslX,
                        REFPAR_DERIV_STEP_ABSOLUTE,false,false,true,true,1.,1.);
      tmp.AssignClock(mClockScatterer);
      tmp.SetDerivStep(1e-5);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp(this->GetName()+"_y",&mXYZ(1),0,1,
                        gpRefParTypeScattTranslY,
                        REFPAR_DERIV_STEP_ABSOLUTE,false,false,true,true,1.,1.);
      tmp.AssignClock(mClockScatterer);
      tmp.SetDerivStep(1e-5);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp(this->GetName()+"_z",&mXYZ(2),0,1,
                        gpRefParTypeScattTranslZ,
                        REFPAR_DERIV_STEP_ABSOLUTE,false,false,true,true,1.,1.);
      tmp.AssignClock(mClockScatterer);
      tmp.SetDerivStep(1e-5);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp(this->GetName()+"_Occ",&mOccupancy,0,1,
                        gpRefParTypeScattOccup,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,1.,1.);
      tmp.AssignClock(mClockScatterer);
      tmp.SetDerivStep(1e-5);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp(this->GetName()+"Q0",&(mQuat.Q0()),0,1,
                        gpRefParTypeScattOrient,
                        REFPAR_DERIV_STEP_ABSOLUTE,false,false,true,false,1.,1.);
      tmp.AssignClock(mClockScatterer);
      tmp.SetGlobalOptimStep(0.04);
      tmp.SetDerivStep(1e-4);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp(this->GetName()+"Q1",&(mQuat.Q1()),0,1,
                        gpRefParTypeScattOrient,
                        REFPAR_DERIV_STEP_ABSOLUTE,false,false,true,false,1.,1.);
      tmp.AssignClock(mClockScatterer);
      tmp.SetGlobalOptimStep(0.04);
      tmp.SetDerivStep(1e-4);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp(this->GetName()+"Q2",&(mQuat.Q2()),0,1,
                        gpRefParTypeScattOrient,
                        REFPAR_DERIV_STEP_ABSOLUTE,false,false,true,false,1.,1.);
      tmp.AssignClock(mClockScatterer);
      tmp.SetGlobalOptimStep(0.04);
      tmp.SetDerivStep(1e-4);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp(this->GetName()+"Q3",&(mQuat.Q3()),0,1,
                        gpRefParTypeScattOrient,
                        REFPAR_DERIV_STEP_ABSOLUTE,false,false,true,false,1.,1.);
      tmp.AssignClock(mClockScatterer);
      tmp.SetGlobalOptimStep(0.04);
      tmp.SetDerivStep(1e-4);
      this->AddPar(tmp);
   }
   mLocalParamSet=this->CreateParamSet("saved parameters for local minimization");
   this->InitOptions();
   mClockScatterer.AddChild(mClockAtomList);
   mClockScatterer.AddChild(mClockBondList);
   mClockScatterer.AddChild(mClockBondAngleList);
   mClockScatterer.AddChild(mClockDihedralAngleList);
   mClockScatterer.AddChild(mClockRingList);
   mClockScatterer.AddChild(mClockRigidGroup);
   mClockScatterer.AddChild(mClockAtomPosition);
   mClockScatterer.AddChild(mClockAtomScattPow);
   mClockScatterer.AddChild(mClockOrientation);
   
   stringstream str;
   old.XMLOutput(str);
   XMLCrystTag tag(str);
   this->XMLInput(str,tag);
   VFN_DEBUG_EXIT("Molecule::Molecule(old&)",5)
}

Molecule::~Molecule()
{
   VFN_DEBUG_ENTRY("Molecule::~Molecule()",5)
   if(mDeleteSubObjInDestructor)
   {
       {
          vector<MolAtom*>::iterator pos;
          for(pos=mvpAtom.begin();pos!=mvpAtom.end();pos++) delete *pos;
       }
       {
          vector<MolBond*>::iterator pos;
          for(pos=mvpBond.begin();pos!=mvpBond.end();pos++) delete *pos;
       }
       {
          vector<MolBondAngle*>::iterator pos;
          for(pos=mvpBondAngle.begin();pos!=mvpBondAngle.end();pos++) delete *pos;
       }
       {
          vector<MolDihedralAngle*>::iterator pos;
          for(pos=mvpDihedralAngle.begin();pos!=mvpDihedralAngle.end();pos++) delete *pos;
       }
   }
   VFN_DEBUG_EXIT("Molecule::~Molecule()",5)
}

Molecule* Molecule::CreateCopy() const
{
   VFN_DEBUG_MESSAGE("Molecule::CreateCopy()",5)
   return new Molecule(*this);
}

const string& Molecule::GetClassName() const
{
   static const string className="Molecule";
   return className;
}

void Molecule::SetName(const string &name)
{
   if(mName==name) return;
   this->RefinableObj::SetName(name);
   // Set parameter's name including the Molecule's name
   this->GetPar(&mXYZ(0)).SetName(mName+"_x");
   this->GetPar(&mXYZ(1)).SetName(mName+"_y");
   this->GetPar(&mXYZ(2)).SetName(mName+"_z");
   this->GetPar(&mOccupancy).SetName(mName+"_Occ");
   this->GetPar(&(mQuat.Q0())).SetName(mName+"Q0");
   this->GetPar(&(mQuat.Q1())).SetName(mName+"Q1");
   this->GetPar(&(mQuat.Q2())).SetName(mName+"Q2");
   this->GetPar(&(mQuat.Q3())).SetName(mName+"Q3");
}

void Molecule::Print()const
{
   VFN_DEBUG_MESSAGE("Molecule::Print()",5)
   this->XMLOutput(cout);
}

void Molecule::XMLOutput(ostream &os,int indent)const
{
   VFN_DEBUG_ENTRY("Molecule::XMLOutput()",4)
   
   // :KLUDGE: this may be dangerous if the molecule is beaing refined !
   this->ResetRigidGroupsPar();
   
   for(int i=0;i<indent;i++) os << "  " ;
   XMLCrystTag tag("Molecule");
   tag.AddAttribute("Name",mName);
   stringstream sst;
   sst<<mMDMoveFreq;
   tag.AddAttribute("MDMoveFreq",sst.str());
   sst.str("");
   sst<<mMDMoveEnergy;
   tag.AddAttribute("MDMoveEnergy",sst.str());
   sst.str("");
   sst<<mLogLikelihoodScale;
   tag.AddAttribute("LogLikelihoodScale",sst.str());
   os <<tag<<endl;
   indent++;
   
   mQuat.Normalize();
   mQuat.XMLOutput(os,indent);
   
   this->GetPar(mXYZ.data()+0).XMLOutput(os,"x",indent);
   os <<endl;
   
   this->GetPar(mXYZ.data()+1).XMLOutput(os,"y",indent);
   os <<endl;
   
   this->GetPar(mXYZ.data()+2).XMLOutput(os,"z",indent);
   os <<endl;
   
   this->GetPar(&mOccupancy).XMLOutput(os,"Occup",indent);
   os <<endl<<endl;
   
   for(unsigned int i=0;i<this->GetNbOption();i++)
   {
      this->GetOption(i).XMLOutput(os,indent);
      os <<endl;
   }
   os <<endl;
   
   {
      vector<MolAtom*>::const_iterator pos;
      for(pos=mvpAtom.begin();pos!=mvpAtom.end();++pos)
         (*pos)->XMLOutput(os,indent);
   }
   {
      vector<MolBond*>::const_iterator pos;
      for(pos=mvpBond.begin();pos!=mvpBond.end();++pos)
         (*pos)->XMLOutput(os,indent);
   }
   {
      vector<MolBondAngle*>::const_iterator pos;
      for(pos=mvpBondAngle.begin();pos!=mvpBondAngle.end();++pos)
         (*pos)->XMLOutput(os,indent);
   }
   {
      vector<MolDihedralAngle*>::const_iterator pos;
      for(pos=mvpDihedralAngle.begin();pos!=mvpDihedralAngle.end();++pos)
         (*pos)->XMLOutput(os,indent);
   }
   {
      vector<RigidGroup*>::const_iterator pos;
      for(pos=mvRigidGroup.begin();pos!=mvRigidGroup.end();++pos)
      {
         XMLCrystTag tagg("RigidGroup",false,true);
         for(set<MolAtom *>::const_iterator at=(*pos)->begin();at!=(*pos)->end();++at)
            tagg.AddAttribute("Atom",(*at)->GetName());
         /*
         tagg.AddAttribute("Q0",(*pos)->mQuat.Q0());
         tagg.AddAttribute("Q1",(*pos)->mQuat.Q1());
         tagg.AddAttribute("Q2",(*pos)->mQuat.Q2());
         tagg.AddAttribute("Q3",(*pos)->mQuat.Q3());
         tagg.AddAttribute("dX",(*pos)->mX);
         tagg.AddAttribute("dY",(*pos)->mY);
         tagg.AddAttribute("dZ",(*pos)->mZ);
         */
         for(int i=0;i<indent;i++) os << "  " ;
         os <<tagg<<endl;
      }
   }
   if(this->GetCenterAtom()!=0)
   {
         
      XMLCrystTag tagg("CenterAtom",false,true);
      tagg.AddAttribute("Name",this->GetCenterAtom()->GetName());
      for(int i=0;i<indent;i++) os << "  " ;
      os <<tagg<<endl;
   }
   
   indent--;
   tag.SetIsEndTag(true);
   for(int i=0;i<indent;i++) os << "  " ;
   os <<tag<<endl;
   VFN_DEBUG_EXIT("Molecule::XMLOutput()",4)
}

void Molecule::XMLInput(istream &is,const XMLCrystTag &tag)
{
   VFN_DEBUG_ENTRY("Molecule::XMLInput()",5)
   for(unsigned int i=0;i<tag.GetNbAttribute();i++)
   {
      if("Name"==tag.GetAttributeName(i))
      {
         mName=tag.GetAttributeValue(i);
      }
      if("MDMoveFreq"==tag.GetAttributeName(i))
      {
         mMDMoveFreq=atof(tag.GetAttributeValue(i).c_str());
      }
      if("MDMoveEnergy"==tag.GetAttributeName(i))
      {
         mMDMoveEnergy=atof(tag.GetAttributeValue(i).c_str());
      }
      if("LogLikelihoodScale"==tag.GetAttributeName(i))
      {
         mLogLikelihoodScale=atof(tag.GetAttributeValue(i).c_str());
      }
   }
   while(true)
   {
      XMLCrystTag tagg(is);
      if(("Molecule"==tagg.GetName())&&tagg.IsEndTag())
      {
         //this->XMLOutput(cout);
         this->UpdateDisplay();
         VFN_DEBUG_EXIT("Molecule::XMLInput():"<<this->GetName(),5)
         return;
      }
      if("Quaternion"==tagg.GetName())
      {
         mQuat.XMLInput(is,tagg);
      }
      if("Atom"==tagg.GetName())
      {
         this->AddAtom(0.,0.,0.,(ScatteringPower *)0,"",false);
         mvpAtom.back()->XMLInput(is,tagg);
      }
      if("Bond"==tagg.GetName())
      {
         this->AddBond(this->GetAtom(0),this->GetAtom(1),1.5,.01,.05,1.,false);
         mvpBond.back()->XMLInput(is,tagg);
      }
      if("BondAngle"==tagg.GetName())
      {
         this->AddBondAngle(this->GetAtom(0),this->GetAtom(1),this->GetAtom(2),1.5,.01,.05,false);
         mvpBondAngle.back()->XMLInput(is,tagg);
      }
      if("DihedralAngle"==tagg.GetName())
      {
         this->AddDihedralAngle(this->GetAtom(0),this->GetAtom(1),
                                this->GetAtom(2),this->GetAtom(3),1.5,.01,.05,false);
         mvpDihedralAngle.back()->XMLInput(is,tagg);
      }
      if("RigidGroup"==tagg.GetName())
      {
         RigidGroup s;
         for(unsigned int i=0;i<tagg.GetNbAttribute();i++)
            if("Atom"==tagg.GetAttributeName(i))
               s.insert(&(this->GetAtom(tagg.GetAttributeValue(i))));
         this->AddRigidGroup(s);
      }
      if("CenterAtom"==tagg.GetName())
      {
         RigidGroup s;
         for(unsigned int i=0;i<tagg.GetNbAttribute();i++)
            if("Name"==tagg.GetAttributeName(i))
               this->SetCenterAtom(this->GetAtom(tagg.GetAttributeValue(i)));
      }
      if("Option"==tagg.GetName())
      {
         for(unsigned int i=0;i<tagg.GetNbAttribute();i++)
            if("Name"==tagg.GetAttributeName(i)) 
               mOptionRegistry.GetObj(tagg.GetAttributeValue(i)).XMLInput(is,tagg);
         continue;
      }
      if("Par"==tagg.GetName())
      {
         for(unsigned int i=0;i<tagg.GetNbAttribute();i++)
         {
            if("Name"==tagg.GetAttributeName(i))
            {
               if("x"==tagg.GetAttributeValue(i))
               {
                  this->GetPar(mXYZ.data()+0).XMLInput(is,tagg);
                  break;
               }
               if("y"==tagg.GetAttributeValue(i))
               {
                  this->GetPar(mXYZ.data()+1).XMLInput(is,tagg);
                  break;
               }
               if("z"==tagg.GetAttributeValue(i))
               {
                  this->GetPar(mXYZ.data()+2).XMLInput(is,tagg);
                  break;
               }
               if("Occup"==tagg.GetAttributeValue(i))
               {
                  this->GetPar(&mOccupancy).XMLInput(is,tagg);
                  break;
               }
            }
         }
      }
   }
   VFN_DEBUG_EXIT("Molecule::XMLInput()",5)
}

void Molecule::UpdateDisplay()const
{
   this->ResetRigidGroupsPar();
   this->RefinableObj::UpdateDisplay();
}

void Molecule::BeginOptimization(const bool allowApproximations,const bool enableRestraints)
{
   if(this->IsBeingRefined())
   {
      // RefinableObj::BeginOptimization() will increase the optimization depth
      this->RefinableObj::BeginOptimization(allowApproximations,enableRestraints);
      return;
   }
   TAU_PROFILE("Molecule::BeginOptimization()","void (bool,bool)",TAU_DEFAULT);
   #if 1 // Is doing this automatically too dangerous ?
   if(  (!mIsSelfOptimizing)
      &&(mAutoOptimizeConformation.GetChoice()==0))
   {
      if(this->GetLogLikelihood()>(mvpRestraint.size()*500))
      {
         (*fpObjCrystInformUser)("Optimizing initial conformation of Molecule:"+this->GetName());
         this->OptimizeConformation(100000,(REAL)(mvpRestraint.size()));
         (*fpObjCrystInformUser)("");
      }
      mAutoOptimizeConformation.SetChoice(1);
   }
   #endif
   
   RefinableObjClock clockConf, clockMode;
   clockConf=mClockAtomList;
   if(clockConf<mClockBondList) clockConf=mClockBondList;
   if(clockConf<mClockBondAngleList) clockConf=mClockBondAngleList;
   if(clockConf<mClockDihedralAngleList) clockConf=mClockDihedralAngleList;
   if(clockConf<mClockRigidGroup) clockConf=mClockRigidGroup;
   if(clockConf<mClockAtomScattPow) clockConf=mClockAtomScattPow;
   
   clockMode=mClockConnectivityTable;
   if(clockMode<mClockRingList) clockMode=mClockRingList;
   if(clockMode<mClockRotorGroup) clockMode=mClockRotorGroup;
   if(clockMode<mClockFlipGroup) clockMode=mClockFlipGroup;
   if(clockMode<mClockStretchModeBondLength) clockMode=mClockStretchModeBondLength;
   if(clockMode<mClockStretchModeBondAngle) clockMode=mClockStretchModeBondAngle;
   if(clockMode<mClockStretchModeTorsion) clockMode=mClockStretchModeTorsion;
   if(clockMode<mClockStretchModeTwist) clockMode=mClockStretchModeTwist;
   if(clockMode<mClockMDAtomGroup) clockMode=mClockMDAtomGroup;


   if( (!mIsSelfOptimizing) && (clockMode<clockConf))
   {
      #if 0
      this->BuildRotorGroup();
      #endif
      if(mFlexModel.GetChoice()!=1)
      {
         this->BuildFlipGroup();
         this->BuildRingList();
         this->BuildStretchModeBondLength();
         this->BuildStretchModeBondAngle();
         this->BuildStretchModeTorsion();
         //this->BuildStretchModeTwist();
         this->TuneGlobalOptimRotationAmplitude();
         this->BuildStretchModeGroups();
         this->BuildMDAtomGroups();
      }
   }
   if(mOptimizeOrientation.GetChoice()==1)
   {
     this->GetPar(&(mQuat.Q0())).SetIsFixed(true);
     this->GetPar(&(mQuat.Q1())).SetIsFixed(true);
     this->GetPar(&(mQuat.Q2())).SetIsFixed(true);
     this->GetPar(&(mQuat.Q3())).SetIsFixed(true);
   }
   else
   {
     this->GetPar(&(mQuat.Q0())).SetIsFixed(false);
     this->GetPar(&(mQuat.Q1())).SetIsFixed(false);
     this->GetPar(&(mQuat.Q2())).SetIsFixed(false);
     this->GetPar(&(mQuat.Q3())).SetIsFixed(false);
   }
   if(1==mFlexModel.GetChoice())
   {// Molecule is a rigid body - fix all individual atomic parameters
      for(vector<MolAtom*>::iterator at=this->GetAtomList().begin();at!=this->GetAtomList().end();++at)
      {
         this->GetPar(&((*at)->X())).SetIsFixed(true);
         this->GetPar(&((*at)->Y())).SetIsFixed(true);
         this->GetPar(&((*at)->Z())).SetIsFixed(true);
      }
   }
   else
   {// Molecule is flexible - rigid groups are handled below
      for(vector<MolAtom*>::iterator at=this->GetAtomList().begin();at!=this->GetAtomList().end();++at)
      {
         this->GetPar(&((*at)->X())).SetIsFixed(false);
         this->GetPar(&((*at)->Y())).SetIsFixed(false);
         this->GetPar(&((*at)->Z())).SetIsFixed(false);
      }
   }
   #ifdef RIGID_BODY_STRICT_EXPERIMENTAL
   // Block individual refinable parameters from atoms in rigid groups
   // And create the index of the atoms
   for(vector<RigidGroup *>::iterator pos=this->GetRigidGroupList().begin();pos!=this->GetRigidGroupList().end();++pos)
   {
      // Init the translation & rotation parameters (ignored outside an optimization)
      (*pos)->mX=0;
      (*pos)->mY=0;
      (*pos)->mZ=0;
      (*pos)->mQuat.Q0()=1;
      (*pos)->mQuat.Q1()=0;
      (*pos)->mQuat.Q2()=0;
      (*pos)->mQuat.Q3()=0;
      (*pos)->mvIdx.clear();
      for(set<MolAtom *>::iterator at=(*pos)->begin();at!=(*pos)->end();++at)
      {
         this->GetPar(&((*at)->X())).SetIsFixed(true);
         this->GetPar(&((*at)->Y())).SetIsFixed(true);
         this->GetPar(&((*at)->Z())).SetIsFixed(true);
         for(unsigned int i=0;i<this->GetNbComponent();++i)
            if(&(this->GetAtom(i))==*at)
            {
              (*pos)->mvIdx.insert(i);
              break;
            }
      }
   }
   #endif
   
   this->RefinableObj::BeginOptimization(allowApproximations,enableRestraints);
   mRandomConformChangeNbTest=0;
   mRandomConformChangeNbAccept=0;
   mRandomConformChangeTemp=1.;//(REAL)this->GetNbComponent();
}

void Molecule::EndOptimization()
{
   /*
   if(mOptimizationDepth>1)
   {
      this->RefinableObj::EndOptimization();
      return;
   }
   */
   #ifdef RIGID_BODY_STRICT_EXPERIMENTAL
   // Un-block individual refinable parameters from atoms in rigid groups
   for(vector<RigidGroup *>::iterator pos=this->GetRigidGroupList().begin();pos!=this->GetRigidGroupList().end();++pos)
   {
      for(set<MolAtom *>::iterator at=(*pos)->begin();at!=(*pos)->end();++at)
      {
         this->GetPar(&((*at)->X())).SetIsFixed(false);
         this->GetPar(&((*at)->Y())).SetIsFixed(false);
         this->GetPar(&((*at)->Z())).SetIsFixed(false);
      }
   }
   // Apply the translations & rotations of the rigid group parameters, and
   // use this as the newly stored atomic coordinates.
   this->ResetRigidGroupsPar();
   #endif
   this->RefinableObj::EndOptimization();
}

void Molecule::RandomizeConfiguration()
{
   TAU_PROFILE("Molecule::RandomizeConfiguration()","void ()",TAU_DEFAULT);
   VFN_DEBUG_ENTRY("Molecule::RandomizeConfiguration()",4)
   
   if(  (!mIsSelfOptimizing)
      &&(this->GetLogLikelihood()>(mvpRestraint.size()*500))
      &&(mAutoOptimizeConformation.GetChoice()==0))
   {// This is only done once, for a newly-created molecule with atoms not conforming to restraints
      (*fpObjCrystInformUser)("Optimizing initial conformation of Molecule:"+this->GetName());
      this->OptimizeConformation(100000,(REAL)(mvpRestraint.size()));
      (*fpObjCrystInformUser)("");
   }
   
   if(   (!(this->IsBeingRefined())) 
      && (mvStretchModeTorsion.size()==0)
      &&(mvStretchModeBondAngle.size()==0)
      &&(mvStretchModeBondLength.size()==0)
      &&(mvStretchModeTwist.size()==0)
      &&(mvMDAtomGroup.size()==0))
   {
      //This will build stretch modes & MD groups
      if(mFlexModel.GetChoice()!=1)
      {
         this->BuildStretchModeTorsion();
         this->TuneGlobalOptimRotationAmplitude();
         //this->BuildStretchModeGroups();
         this->BuildMDAtomGroups();
      }
   }
   
   #if 0
   this->BuildRotorGroup();
   if((mFlexModel.GetChoice()==0)||(mFlexModel.GetChoice()==2))
      for(list<RotorGroup>::const_iterator pos=mvRotorGroupTorsion.begin();
          pos!=mvRotorGroupTorsion.end();++pos)
      {
         const REAL angle=(REAL)rand()*2.*M_PI/(REAL)RAND_MAX;
         this->RotateAtomGroup(*(pos->mpAtom1),*(pos->mpAtom2),
                               pos->mvRotatedAtomList,angle);
      }
   if(mFlexModel.GetChoice()==0)
      for(list<RotorGroup>::const_iterator pos=mvRotorGroupTorsionSingleChain.begin();
          pos!=mvRotorGroupTorsionSingleChain.end();++pos)
      {
         const REAL angle=(REAL)rand()*2.*M_PI/(REAL)RAND_MAX;
         this->RotateAtomGroup(*(pos->mpAtom1),*(pos->mpAtom2),
                               pos->mvRotatedAtomList,angle);
      }
   if(mFlexModel.GetChoice()==0)
      for(list<RotorGroup>::const_iterator pos=mvRotorGroupInternal.begin();
          pos!=mvRotorGroupInternal.end();++pos)
      {
         const REAL angle=(REAL)rand()*2.*M_PI/(REAL)RAND_MAX;
         this->RotateAtomGroup(*(pos->mpAtom1),*(pos->mpAtom2),
                               pos->mvRotatedAtomList,angle);
      }
   #else
   for(list<StretchModeTorsion>::const_iterator
         pos=mvStretchModeTorsion.begin();
       pos!=mvStretchModeTorsion.end();++pos)
   {
      const REAL amp=2*M_PI*rand()/(REAL)RAND_MAX;
      this->DihedralAngleRandomChange(*pos,amp,true);
   }
   // Molecular dynamics moves
   if((mvMDAtomGroup.size()>0)&&(mMDMoveFreq>0)&&(mMDMoveEnergy>0))
   {
      // Random initial speed for all atoms
      map<MolAtom*,XYZ> v0;
      for(vector<MolAtom*>::iterator at=this->GetAtomList().begin();at!=this->GetAtomList().end();++at)
         v0[*at]=XYZ(rand()/(REAL)RAND_MAX+0.5,rand()/(REAL)RAND_MAX+0.5,rand()/(REAL)RAND_MAX+0.5);
      
      const REAL nrj0=mMDMoveEnergy*( this->GetBondList().size()
                                     +this->GetBondAngleList().size()
                                     +this->GetDihedralAngleList().size());
      map<RigidGroup*,std::pair<XYZ,XYZ> > vr;
      this->MolecularDynamicsEvolve(v0, 5000,0.004,
                                    this->GetBondList(),
                                    this->GetBondAngleList(),
                                    this->GetDihedralAngleList(),
                                    vr,nrj0);
   }
   
   #endif
   // this will only change limited parameters i.e. translation
   this->RefinableObj::RandomizeConfiguration();
   #ifdef RIGID_BODY_STRICT_EXPERIMENTAL
   // Init rigid groups translation & rotation parameters to zero
   for(vector<RigidGroup *>::iterator pos=this->GetRigidGroupList().begin();pos!=this->GetRigidGroupList().end();++pos)
   {
      // Init the translation & rotation parameters (ignored outside an optimization
      (*pos)->mX=0;
      (*pos)->mY=0;
      (*pos)->mZ=0;
      (*pos)->mQuat.Q0()=1;
      (*pos)->mQuat.Q1()=0;
      (*pos)->mQuat.Q2()=0;
      (*pos)->mQuat.Q3()=0;
   }
   #endif
   if(mOptimizeOrientation.GetChoice()==0)
   {//Rotate around an arbitrary vector
      const REAL amp=M_PI/RAND_MAX;
      mQuat *= Quaternion::RotationQuaternion
                  ((2.*(REAL)rand()-(REAL)RAND_MAX)*amp,
                   (REAL)rand(),(REAL)rand(),(REAL)rand());
      mQuat.Normalize();
      mClockOrientation.Click();
   }
   VFN_DEBUG_EXIT("Molecule::RandomizeConfiguration()",4)
}

void Molecule::GlobalOptRandomMove(const REAL mutationAmplitude,
                                       const RefParType *type)
{
   if(mRandomMoveIsDone) return;
   if(mIsSelfOptimizing) 
   {
      this->RefinableObj::GlobalOptRandomMove(mutationAmplitude,type);
      mQuat.Normalize();
      VFN_DEBUG_EXIT("Molecule::GlobalOptRandomMove()",4)
      return;
   }
   TAU_PROFILE("Molecule::GlobalOptRandomMove()","void (REAL,RefParType*)",TAU_DEFAULT);
   TAU_PROFILE_TIMER(timer1,"Molecule::GlobalOptRandomMove 1","", TAU_FIELD);
   TAU_PROFILE_TIMER(timer2,"Molecule::GlobalOptRandomMove 2","", TAU_FIELD);
   TAU_PROFILE_TIMER(timer3,"Molecule::GlobalOptRandomMove 3","", TAU_FIELD);
   TAU_PROFILE_TIMER(timer4,"Molecule::GlobalOptRandomMove 4","", TAU_FIELD);
   TAU_PROFILE_TIMER(timer5,"Molecule::GlobalOptRandomMove 5","", TAU_FIELD);
   VFN_DEBUG_ENTRY("Molecule::GlobalOptRandomMove()",4)
   mClockScatterer.Click();
   
   #if 1
   // From time to time, just do one flip
   if(  (mFlexModel.GetChoice()!=1)
      &&(gpRefParTypeScattConform->IsDescendantFromOrSameAs(type))
      &&(mvFlipGroup.size()>0)
      &&(((rand()%100)==0)))
   {
      this->SaveParamSet(mLocalParamSet);
      const REAL llk0=this->GetLogLikelihood()/mLogLikelihoodScale;
      const unsigned long i=rand() % mvFlipGroup.size();
      list<FlipGroup>::iterator pos=mvFlipGroup.begin();
      for(unsigned long j=0;j<i;++j)++pos;
      this->FlipAtomGroup(*pos,true);
      #if 0
      static unsigned long ctflip1=0,ctflip2=0;
      if(pos->mvRotatedChainList.begin()->first==pos->mpAtom0)
      {
         cout <<"TRYING: Flip group from atom "
              <<pos->mpAtom0->GetName()<<",exchanging bonds with "
              <<pos->mpAtom1->GetName()<<" and "
              <<pos->mpAtom2->GetName()<<", resulting in a 180deg rotation of atoms : ";
         for(set<MolAtom*>::iterator pos1=pos->mvRotatedChainList.begin()->second.begin();
             pos1!=pos->mvRotatedChainList.begin()->second.end();++pos1)
            cout<<(*pos1)->GetName()<<"  ";
      }
      else
      {
         cout <<"TRYING: Flip group with respect to: "
              <<pos->mpAtom1->GetName()<<"-"
              <<pos->mpAtom0->GetName()<<"-"
              <<pos->mpAtom2->GetName()<<" : ";
         for(list<pair<const MolAtom *,set<MolAtom*> > >::const_iterator 
             chain=pos->mvRotatedChainList.begin();
             chain!=pos->mvRotatedChainList.end();++chain)
         {
            cout<<"    -"<<chain->first->GetName()<<":";
            for(set<MolAtom*>::const_iterator pos1=chain->second.begin();
                pos1!=chain->second.end();++pos1)
               cout<<(*pos1)->GetName()<<"  ";
         }
      }
      ctflip1++;
      if((this->GetLogLikelihood()/mLogLikelihoodScale-llk0)>.3*llk0)
      {
        cout<<"      FLIP REJECTED ("<<int(float(ctflip2)/ctflip1*100)<<"% accepted): llk="<<llk0<<"  ->  "<<this->GetLogLikelihood()/mLogLikelihoodScale<<endl;
        this->RestoreParamSet(mLocalParamSet);
      }
      else
      {
         ctflip2++;
         cout<<"      FLIP ACCEPTED ("<<float(ctflip2)/ctflip1*100<<"% accepted)"<<endl;
      }
      #else
      if((this->GetLogLikelihood()/mLogLikelihoodScale-llk0)>.3*llk0)
         this->RestoreParamSet(mLocalParamSet);
      #endif
   }
   else
   #endif
   {
      TAU_PROFILE_START(timer1);
      if(mOptimizeOrientation.GetChoice()==0)
      {//Rotate around an arbitrary vector
         static const REAL amp=mBaseRotationAmplitude/(REAL)RAND_MAX;
         REAL mult=1.0;
         if((1==mFlexModel.GetChoice())||(mvRotorGroupTorsion.size()<2)) mult=2.0;
         mQuat *= Quaternion::RotationQuaternion
                     ((2.*(REAL)rand()-(REAL)RAND_MAX)*amp*mutationAmplitude*mult,
                      (REAL)rand(),(REAL)rand(),(REAL)rand());
         mQuat.Normalize();
         mClockOrientation.Click();
      }
      // Occupancy
      if(gpRefParTypeScattOccup->IsDescendantFromOrSameAs(type))
         this->RefinableObj::GlobalOptRandomMove(mutationAmplitude,gpRefParTypeScattOccup);
      mRandomMoveIsDone=false;
      //translation
      if(gpRefParTypeScattTransl->IsDescendantFromOrSameAs(type))
         this->RefinableObj::GlobalOptRandomMove(mutationAmplitude,gpRefParTypeScattTransl);
      mRandomMoveIsDone=false;
      TAU_PROFILE_STOP(timer1);
      if(gpRefParTypeScattConform->IsDescendantFromOrSameAs(type))
      {
         if(mFlexModel.GetChoice()!=1)
         {
            #if 1 // Move as many atoms as possible
            if((mvMDFullAtomGroup.size()>3)&&(rand()<(RAND_MAX*mMDMoveFreq)))
            {
               #if 0
               // Use one center for the position of an impulsion, applied to all atoms with an exponential decrease
               // Determine extent of atom group
               REAL xmin=(*mvMDFullAtomGroup.begin())->GetX();
               REAL xmax=(*mvMDFullAtomGroup.begin())->GetX();
               REAL ymin=(*mvMDFullAtomGroup.begin())->GetY();
               REAL ymax=(*mvMDFullAtomGroup.begin())->GetY();
               REAL zmin=(*mvMDFullAtomGroup.begin())->GetZ();
               REAL zmax=(*mvMDFullAtomGroup.begin())->GetZ();
               for(set<MolAtom*>::iterator at=mvMDFullAtomGroup.begin();at!=mvMDFullAtomGroup.end();++at)
               {
                  if((*at)->GetX()<xmin) xmin=(*at)->GetX();
                  if((*at)->GetX()>xmax) xmax=(*at)->GetX();
                  if((*at)->GetY()<ymin) ymin=(*at)->GetY();
                  if((*at)->GetY()>ymax) ymax=(*at)->GetY();
                  if((*at)->GetZ()<zmin) zmin=(*at)->GetZ();
                  if((*at)->GetZ()>zmax) zmax=(*at)->GetZ();
               }
               // Apply a gaussian impulsion to part of the atom group (FWHM=1/3 of group size)
               REAL dx=(xmax-xmin)/5.,dy=(ymax-ymin)/5.,dz=(zmax-zmin)/5.;
               if(dx<2) dx=2;
               if(dy<2) dy=2;
               if(dz<2) dz=2;
               const REAL xc=xmin+rand()/(REAL)RAND_MAX*(xmax-xmin);
               const REAL yc=ymin+rand()/(REAL)RAND_MAX*(ymax-ymin);
               const REAL zc=zmin+rand()/(REAL)RAND_MAX*(zmax-zmin);
               map<MolAtom*,XYZ> v0;
               const REAL ax=-4.*log(2.)/(dx*dx);
               const REAL ay=-4.*log(2.)/(dy*dy);
               const REAL az=-4.*log(2.)/(dz*dz);
               for(set<MolAtom*>::iterator at=mvMDFullAtomGroup.begin();at!=mvMDFullAtomGroup.end();++at)
                  v0[*at]=XYZ(exp(ax*((*at)->GetX()-xc)*((*at)->GetX()-xc)),
                              exp(ay*((*at)->GetY()-yc)*((*at)->GetY()-yc)),
                              exp(az*((*at)->GetZ()-zc)*((*at)->GetZ()-zc)));
               #else
               // Use one atom for the center of the impulsion, 'push' atoms depending on distance & connectivity table
               map<MolAtom*,XYZ> v0;
               for(set<MolAtom*>::iterator at=this->mvMDFullAtomGroup.begin();at!=this->mvMDFullAtomGroup.end();++at)
                  v0[*at]=XYZ(0,0,0);
               std::map<MolAtom*,unsigned long> pushedAtoms;
               unsigned long idx=rand()%v0.size();
               set<MolAtom*>::iterator at0=this->mvMDFullAtomGroup.begin();
               for(unsigned int i=0;i<idx;i++) at0++;
               const REAL xc=(*at0)->GetX();
               const REAL yc=(*at0)->GetY();
               const REAL zc=(*at0)->GetZ();
               const map<MolAtom *,set<MolAtom *> > *pConnect=&(this-> GetConnectivityTable());
               ExpandAtomGroupRecursive(*at0,*pConnect,pushedAtoms,3);
               REAL ux,uy,uz,n=0;
               while(n<1)
               {
                  ux=REAL(rand()-RAND_MAX/2);
                  uy=REAL(rand()-RAND_MAX/2);
                  uz=REAL(rand()-RAND_MAX/2);
                  n=sqrt(ux*ux+uy*uy+uz*uz);
               }
               ux=ux/n;uy=uy/n;uz=uz/n;
               const REAL a=-4.*log(2.)/(2*2);//FWHM=2 Angstroems
               if(rand()%2==0)
                  for(map<MolAtom*,unsigned long>::iterator at=pushedAtoms.begin() ;at!=pushedAtoms.end();++at)
                     v0[at->first]=XYZ(ux*exp(a*(at->first->GetX()-xc)*(at->first->GetX()-xc)),
                                 uy*exp(a*(at->first->GetY()-yc)*(at->first->GetY()-yc)),
                                 uz*exp(a*(at->first->GetZ()-zc)*(at->first->GetZ()-zc)));
               else
                  for(map<MolAtom*, unsigned long>::iterator at=pushedAtoms.begin() ;at!=pushedAtoms.end();++at)
                     v0[at->first]=XYZ((at->first->GetX()-xc)*ux*exp(a*(at->first->GetX()-xc)*(at->first->GetX()-xc)),
                                       (at->first->GetY()-yc)*uy*exp(a*(at->first->GetY()-yc)*(at->first->GetY()-yc)),
                                       (at->first->GetZ()-zc)*uz*exp(a*(at->first->GetZ()-zc)*(at->first->GetZ()-zc)));
               

               #endif
               const REAL nrj0=mMDMoveEnergy*( this->GetBondList().size()
                                      +this->GetBondAngleList().size()
                                      +this->GetDihedralAngleList().size());
               map<RigidGroup*,std::pair<XYZ,XYZ> > vr;
               this->MolecularDynamicsEvolve(v0, int(100*sqrt(mutationAmplitude)),0.004,
                                             this->GetBondList(),
                                             this->GetBondAngleList(),
                                             this->GetDihedralAngleList(),
                                             vr,nrj0);
            }
            #else // Move atoms belonging to a MD group
            if((mvMDAtomGroup.size()>0)&&(rand()<(RAND_MAX*mMDMoveFreq)))
            {
               const unsigned int n=rand()%mvMDAtomGroup.size();
               list<MDAtomGroup>::iterator pos=mvMDAtomGroup.begin();
               for(unsigned int i=0;i<n;++i)++pos;
               map<MolAtom*,XYZ> v0;
               for(set<MolAtom*>::iterator at=pos->mvpAtom.begin();at!=pos->mvpAtom.end();++at)
                  v0[*at]=XYZ(rand()/(REAL)RAND_MAX+0.5,rand()/(REAL)RAND_MAX+0.5,rand()/(REAL)RAND_MAX+0.5);
               
               const REAL nrj0=mMDMoveEnergy*( pos->mvpBond.size()
                                    +pos->mvpBondAngle.size()
                                    +pos->mvpDihedralAngle.size());
               map<RigidGroup*,std::pair<XYZ,XYZ> > vr;
               float nrjMult=1.0+mutationAmplitude*0.2;
               if((rand()%20)==0) nrjMult=4.0;
               this->MolecularDynamicsEvolve(v0, int(100*sqrt(mutationAmplitude)),0.004,
                                             pos->mvpBond,
                                             pos->mvpBondAngle,
                                             pos->mvpDihedralAngle,
                                             vr,nrj0*nrjMult);
            }
            #endif
            else
            {
            #if 0 // For tests
            mLogLikelihood=0;
            for(vector<MolBond*>::const_iterator pos=mvpBond.begin();pos!=mvpBond.end();++pos)
               mLogLikelihood+=(*pos)->GetLogLikelihood(true,true);
            for(vector<MolBondAngle*>::const_iterator pos=mvpBondAngle.begin();pos!=mvpBondAngle.end();++pos)
               mLogLikelihood+=(*pos)->GetLogLikelihood(true,true);
            for(vector<MolDihedralAngle*>::const_iterator pos=mvpDihedralAngle.begin();pos!=mvpDihedralAngle.end();++pos)
               mLogLikelihood+=(*pos)->GetLogLikelihood(true,true);
            for(list<StretchMode*>::const_iterator mode=mvpStretchModeNotFree.begin();
                mode!=mvpStretchModeNotFree.end();++mode)
            {
               //if((rand()%3)==0)
               {
                  // 2) Get the derivative of the overall LLK for this mode
                  (*mode)->CalcDeriv();
                  REAL llk=0;
                  for(map<const MolBond*,REAL>::const_iterator pos=(*mode)->mvpBrokenBond.begin();
                      pos!=(*mode)->mvpBrokenBond.end();++pos) llk+=pos->first->GetLogLikelihood(false,false);
                  for(map<const MolBondAngle*,REAL>::const_iterator pos=(*mode)->mvpBrokenBondAngle.begin();
                      pos!=(*mode)->mvpBrokenBondAngle.end();++pos) llk+=pos->first->GetLogLikelihood(false,false);
                  for(map<const MolDihedralAngle*,REAL>::const_iterator pos=(*mode)->mvpBrokenDihedralAngle.begin();
                      pos!=(*mode)->mvpBrokenDihedralAngle.end();++pos) llk+=pos->first->GetLogLikelihood(false,false);
                  // 3) Calculate MD move. base step =0.1 A (accelerated moves may go faster)
                  REAL change=(2.*(REAL)rand()-(REAL)RAND_MAX)/(REAL)RAND_MAX;
                  // if llk>100, change has to be in the opposite direction
                  // For a single restraint, sqrt(llk)=dx/sigma, so do not go above 10*sigma
                  if((*mode)->mLLKDeriv>0)
                  {
                     change -= 0.3*sqrt(llk);
                     if(change<-1) change=-1;
                  }
                  else
                  {
                     change += 0.3*sqrt(llk);
                     if(change>1) change=1;
                  }
                  (*mode)->Print(cout);
                  change *= mutationAmplitude * (*mode)->mBaseAmplitude;
                  cout <<"      Change="<<change<<" (dLLK= "<<(*mode)->mLLKDeriv<<"), llk= "<<llk<<"     ?->"<<llk+(*mode)->mLLKDeriv*change<<endl;
                  //change *= mutationAmplitude * (*mode)->mBaseAmplitude;
                  (*mode)->Stretch(change);
                  llk=0;
                  //(*mode)->RandomStretch(change * mutationAmplitude * (*mode)->mBaseAmplitude);
                  for(map<const MolBond*,REAL>::const_iterator pos=(*mode)->mvpBrokenBond.begin();
                      pos!=(*mode)->mvpBrokenBond.end();++pos) 
                  {
                      cout<<"      "<<pos->first->GetName()<<", llk= "<<pos->first->GetLogLikelihood(false,false)
                          <<"   ?->"<<pos->first->GetLogLikelihood(false,false)+pos->first->GetDeriv((*mode)->mDerivXYZ,true)*change
                          <<"?   (deriv="<<pos->first->GetDeriv((*mode)->mDerivXYZ)<<", "<<pos->first->GetDeriv((*mode)->mDerivXYZ,true);
                      cout<<")  ->" <<pos->first->GetLogLikelihood()<<endl;
                     llk+=pos->first->GetLogLikelihood(false,false);
                  }
                  for(map<const MolBondAngle*,REAL>::const_iterator pos=(*mode)->mvpBrokenBondAngle.begin();
                      pos!=(*mode)->mvpBrokenBondAngle.end();++pos) 
                  {
                      cout<<"      "<<pos->first->GetName()<<", llk= "<<pos->first->GetLogLikelihood(false,false)
                          <<"   ?->"<<pos->first->GetLogLikelihood(false,false)+pos->first->GetDeriv((*mode)->mDerivXYZ,true)*change
                          <<"?   (deriv="<<pos->first->GetDeriv((*mode)->mDerivXYZ)<<", "<<pos->first->GetDeriv((*mode)->mDerivXYZ,true);
                      cout<<")  ->" <<pos->first->GetLogLikelihood()<<endl;
                     llk+=pos->first->GetLogLikelihood(false,false);
                  }
                  for(map<const MolDihedralAngle*,REAL>::const_iterator pos=(*mode)->mvpBrokenDihedralAngle.begin();
                      pos!=(*mode)->mvpBrokenDihedralAngle.end();++pos) 
                  {
                      cout<<"      "<<pos->first->GetName()<<", llk= "<<pos->first->GetLogLikelihood(false,false)
                          <<"   ?->"<<pos->first->GetLogLikelihood(false,false)+pos->first->GetDeriv((*mode)->mDerivXYZ,true)*change
                          <<"?   (deriv="<<pos->first->GetDeriv((*mode)->mDerivXYZ)<<", "<<pos->first->GetDeriv((*mode)->mDerivXYZ,true);
                      cout<<")  ->" <<pos->first->GetLogLikelihood()<<endl;
                     llk+=pos->first->GetLogLikelihood(false,false);
                  }
                  cout <<" -> "<<llk<<endl;
               }
            }
            #else
            // First move free Stretch modes
            TAU_PROFILE_START(timer2);
            for(list<StretchMode*>::iterator mode=mvpStretchModeFree.begin();
                mode!=mvpStretchModeFree.end();++mode)
            {
               if((rand()%2)==0) (*mode)->RandomStretch(mutationAmplitude);
            }
            TAU_PROFILE_STOP(timer2);
            if((rand()%3)==0)
            {
               // Now do an hybrid move for other modes, with a smaller amplitude (<=0.5)
               // 1) Calc LLK and derivatives for restraints
               mLogLikelihood=0;
               TAU_PROFILE_START(timer3);
               for(vector<MolBond*>::const_iterator pos=mvpBond.begin();pos!=mvpBond.end();++pos)
                  mLogLikelihood+=(*pos)->GetLogLikelihood(true,true);
               for(vector<MolBondAngle*>::const_iterator pos=mvpBondAngle.begin();pos!=mvpBondAngle.end();++pos)
                  mLogLikelihood+=(*pos)->GetLogLikelihood(true,true);
               for(vector<MolDihedralAngle*>::const_iterator pos=mvpDihedralAngle.begin();pos!=mvpDihedralAngle.end();++pos)
                  mLogLikelihood+=(*pos)->GetLogLikelihood(true,true);
               TAU_PROFILE_STOP(timer3);

               TAU_PROFILE_START(timer4);
               for(list<StretchMode*>::const_iterator mode=mvpStretchModeNotFree.begin();
                   mode!=mvpStretchModeNotFree.end();++mode)
               {
                  // 2) Choose Stretch modes
                  if((rand()%3)==0)
                  {
                     // 2) Get the derivative of the overall LLK for this mode
                     (*mode)->CalcDeriv();
                     REAL llk=0;
                     for(map<const MolBond*,REAL>::const_iterator pos=(*mode)->mvpBrokenBond.begin();
                         pos!=(*mode)->mvpBrokenBond.end();++pos) llk+=pos->first->GetLogLikelihood(false,false);
                     for(map<const MolBondAngle*,REAL>::const_iterator pos=(*mode)->mvpBrokenBondAngle.begin();
                         pos!=(*mode)->mvpBrokenBondAngle.end();++pos) llk+=pos->first->GetLogLikelihood(false,false);
                     for(map<const MolDihedralAngle*,REAL>::const_iterator pos=(*mode)->mvpBrokenDihedralAngle.begin();
                         pos!=(*mode)->mvpBrokenDihedralAngle.end();++pos) llk+=pos->first->GetLogLikelihood(false,false);
                     REAL change=(2.*(REAL)rand()-(REAL)RAND_MAX)/(REAL)RAND_MAX;
                     // if llk>100, change has to be in the direction minimising the llk
                     if((*mode)->mLLKDeriv>0)
                     {
                        change -= 0.01*llk;
                        if(change<-1) change=-1;
                     }
                     else
                     {
                        change += 0.01*llk;
                        if(change>1) change=1;
                     }
                     if(mutationAmplitude<0.5) change *= mutationAmplitude * (*mode)->mBaseAmplitude;
                     else change *= 0.5 * (*mode)->mBaseAmplitude;
                     (*mode)->Stretch(change);
                  }
               }
               // Here we do not take mLogLikelihoodScale into account
               // :TODO: take into account cases where the lllk cannot go down to 0 because of 
               // combined restraints.
               if( ((rand()%100)==0) && (mLogLikelihood>(mvpRestraint.size()*10)))
                  this->OptimizeConformationSteepestDescent(0.02,5);
               TAU_PROFILE_STOP(timer4);
            }
            
            // Perform MD moves if there are MDAtomGroups
            #if 0
            for(list<MDAtomGroup>::iterator pos=mvMDAtomGroup.begin();pos!=mvMDAtomGroup.end();++pos)
            {
               if((rand()%100)==0)
               {
                  map<MolAtom*,XYZ> v0;
                  for(set<MolAtom*>::iterator at=pos->mvpAtom.begin();at!=pos->mvpAtom.end();++at)
                     v0[*at]=XYZ(rand()/(REAL)RAND_MAX+0.5,rand()/(REAL)RAND_MAX+0.5,rand()/(REAL)RAND_MAX+0.5);
                  
                  const REAL nrj0=20*(pos->mvpBond.size()+pos->mvpBondAngle.size()+pos->mvpDihedralAngle.size());
                  map<RigidGroup*,std::pair<XYZ,XYZ> > vr;
                  this->MolecularDynamicsEvolve(v0, int(100*sqrt(mutationAmplitude)),0.002,
                                                (const vector<MolBond*>) (pos->mvpBond),
                                                (const vector<MolBondAngle*>) (pos->mvpBondAngle),
                                                (const vector<MolDihedralAngle*>) (pos->mvpDihedralAngle),
                                                vr,nrj0);
               }
            }
            #endif
            }
            // Do a steepest descent from time to time
            if((rand()%100)==0) this->OptimizeConformationSteepestDescent(0.02,1);

            mClockLogLikelihood.Click();
            #endif
         }
      }
   }
   if((rand()%100)==0)
   {// From time to time, bring back average position to 0
      REAL x0=0,y0=0,z0=0;
      for(vector<MolAtom*>::iterator pos=mvpAtom.begin();pos!=mvpAtom.end();++pos)
      {
         x0 += (*pos)->X();
         y0 += (*pos)->Y();
         z0 += (*pos)->Z();
      }
      x0 /= mvpAtom.size();
      y0 /= mvpAtom.size();
      z0 /= mvpAtom.size();
      for(vector<MolAtom*>::iterator pos=mvpAtom.begin();pos!=mvpAtom.end();++pos)
      {
         (*pos)->X() -= x0;
         (*pos)->Y() -= y0;
         (*pos)->Z() -= z0;
      }
   }
   mRandomMoveIsDone=true;
   VFN_DEBUG_EXIT("Molecule::GlobalOptRandomMove()",4)
}

REAL Molecule::GetLogLikelihood()const
{
   if(  (mClockLogLikelihood>mClockAtomList)
      &&(mClockLogLikelihood>mClockBondList)
      &&(mClockLogLikelihood>mClockBondAngleList)
      &&(mClockLogLikelihood>mClockDihedralAngleList)
      &&(mClockLogLikelihood>mClockAtomPosition)
      &&(mClockLogLikelihood>mClockScatterer)) return mLogLikelihood*mLogLikelihoodScale;
   TAU_PROFILE("Molecule::GetLogLikelihood()","REAL ()",TAU_DEFAULT);
   mLogLikelihood=this->RefinableObj::GetLogLikelihood();
   mClockLogLikelihood.Click();
   return mLogLikelihood*mLogLikelihoodScale;
}

unsigned int Molecule::GetNbLSQFunction()const
{
   return 1;
}

const CrystVector_REAL& Molecule::GetLSQCalc(const unsigned int) const
{
   mLSQCalc.resize(mvpRestraint.size());
   REAL *p=mLSQCalc.data();
   for(vector<MolBond*>::const_iterator pos=this->GetBondList().begin();pos!=this->GetBondList().end();++pos)
      *p++=(*pos)->GetLength();
   for(vector<MolBondAngle*>::const_iterator pos=this->GetBondAngleList().begin();pos!=this->GetBondAngleList().end();++pos)
      *p++=(*pos)->GetAngle();
   for(vector<MolDihedralAngle*>::const_iterator pos=this->GetDihedralAngleList().begin();pos!=this->GetDihedralAngleList().end();++pos)
      *p++=(*pos)->GetAngle();
   return mLSQCalc;
}

const CrystVector_REAL& Molecule::GetLSQObs(const unsigned int) const
{
   mLSQObs.resize(mvpRestraint.size());
   REAL *p=mLSQObs.data();
   for(vector<MolBond*>::const_iterator pos=this->GetBondList().begin();pos!=this->GetBondList().end();++pos)
      *p++=(*pos)->GetLength0();
   for(vector<MolBondAngle*>::const_iterator pos=this->GetBondAngleList().begin();pos!=this->GetBondAngleList().end();++pos)
      *p++=(*pos)->GetAngle0();
   for(vector<MolDihedralAngle*>::const_iterator pos=this->GetDihedralAngleList().begin();pos!=this->GetDihedralAngleList().end();++pos)
      *p++=(*pos)->GetAngle0();
   return mLSQObs;
}

const CrystVector_REAL& Molecule::GetLSQWeight(const unsigned int) const
{
   //:TODO: USe a clock to avoid re-computation
   mLSQWeight.resize(mvpRestraint.size());
   REAL *p=mLSQWeight.data();
   for(vector<MolBond*>::const_iterator pos=this->GetBondList().begin();pos!=this->GetBondList().end();++pos)
      *p++=1/((*pos)->GetLengthSigma()* (*pos)->GetLengthSigma()+1e-6);
   for(vector<MolBondAngle*>::const_iterator pos=this->GetBondAngleList().begin();pos!=this->GetBondAngleList().end();++pos)
      *p++=1/((*pos)->GetAngleSigma()* (*pos)->GetAngleSigma()+1e-6);
   for(vector<MolDihedralAngle*>::const_iterator pos=this->GetDihedralAngleList().begin();pos!=this->GetDihedralAngleList().end();++pos)
      *p++=1/((*pos)->GetAngleSigma()* (*pos)->GetAngleSigma()+1e-6);
   mLSQWeight*=mLogLikelihoodScale;
   return mLSQWeight;
}

const CrystVector_REAL& Molecule::GetLSQDeriv(const unsigned int n, RefinablePar&par)
{
   //:TODO: return analytical derivatives
   return RefinableObj::GetLSQDeriv(n,par);
}

void Molecule::TagNewBestConfig()const
{
   this->ResetRigidGroupsPar();
   #if 0
   cout<<"Molecule::TagNewBestConfig()"<<endl;
   {
      vector<MolBond*>::const_iterator pos;
      for(pos=mvpBond.begin();pos!=mvpBond.end();++pos)
      {
         cout<<"BondLength="<<(*pos)->GetLength();
         (*pos)->XMLOutput(cout);
      }
   }
   {
      vector<MolBondAngle*>::const_iterator pos;
      for(pos=mvpBondAngle.begin();pos!=mvpBondAngle.end();++pos)
      {
         cout<<"BondAngle="<<(*pos)->GetAngle();
         (*pos)->XMLOutput(cout);
      }
   }
   {
      vector<MolDihedralAngle*>::const_iterator pos;
      for(pos=mvpDihedralAngle.begin();pos!=mvpDihedralAngle.end();++pos)
      {
         cout<<"DihedralAngle="<<(*pos)->GetAngle();
         (*pos)->XMLOutput(cout);
      }
   }
   
   for(list<MDAtomGroup>::iterator pos=mvMDAtomGroup.begin();pos!=mvMDAtomGroup.end();++pos)
   {
      char buf[100];
      for(set<MolAtom*>::iterator at1=pos->mvpAtom.begin();at1!=pos->mvpAtom.end();++at1)
      {
         sprintf(buf,"%5s : ",(*at1)->GetName().c_str());
         cout<<buf;
         for(set<MolAtom*>::iterator at2=at1;at2!=pos->mvpAtom.end();++at2)
         {
            if(at1==at2) continue;
            sprintf(buf,"%5s(%6.3f),",(*at2)->GetName().c_str(),GetBondLength(**at1,**at2));
            cout<<buf;
         }
         cout<<endl;
      }
   }
   #endif
}

int Molecule::GetNbComponent() const { return mvpAtom.size();}

const ScatteringComponentList& Molecule::GetScatteringComponentList() const
{
   VFN_DEBUG_ENTRY("Molecule::GetScatteringComponentList()",3)
   this->UpdateScattCompList();
   VFN_DEBUG_EXIT("Molecule::GetScatteringComponentList()",3)
   return mScattCompList;
}

string Molecule::GetComponentName(const int i) const
{
   //if(mvpAtom[i]->IsDummy()) return "Dummy";
   return mvpAtom[i]->GetName();
} 

ostream& Molecule::POVRayDescription(ostream &os,const CrystalPOVRayOptions &options)const
{
   VFN_DEBUG_ENTRY("Molecule::POVRayDescription()",3)
   const REAL xMin=options.mXmin; const REAL xMax=options.mXmax;
   const REAL yMin=options.mYmin; const REAL yMax=options.mYmax;
   const REAL zMin=options.mZmin; const REAL zMax=options.mZmax;
   if(mvpAtom.size()==0)
   {
      VFN_DEBUG_EXIT("Molecule::POVRayDescription():No atom to display !",4)
      return os;
   }
   REAL en=1;
   this->UpdateScattCompList();
   
   const float colour_bondnonfree[]= { 0.3, .3, .3, 1.0 };
   const float colour_bondfree[]= { 0.7, .7, .7, 1.0 };
   const float colour0[] = {0.0f, 0.0f, 0.0f, 0.0f}; 
   
   os << "// Description of Molecule :" << this->GetName() <<endl;
   vector<CrystMatrix_REAL> vXYZCoords;
   {
      this->GetScatteringComponentList();
      REAL x0,y0,z0;
      for(long i=0;i<mScattCompList.GetNbComponent();++i)
      {
         x0=mScattCompList(i).mX;
         y0=mScattCompList(i).mY;
         z0=mScattCompList(i).mZ;
         vXYZCoords.push_back(this->GetCrystal().GetSpaceGroup().
                        GetAllSymmetrics(x0,y0,z0,false,false,false));
      }
   }
   CrystMatrix_int translate(27,3);
   translate=  -1,-1,-1,
               -1,-1, 0,
               -1,-1, 1,
               -1, 0,-1,
               -1, 0, 0,
               -1, 0, 1,
               -1, 1,-1,
               -1, 1, 0,
               -1, 1, 1,
                0,-1,-1,
                0,-1, 0,
                0,-1, 1,
                0, 0,-1,
                0, 0, 0,
                0, 0, 1,
                0, 1,-1,
                0, 1, 0,
                0, 1, 1,
                1,-1,-1,
                1,-1, 0,
                1,-1, 1,
                1, 0,-1,
                1, 0, 0,
                1, 0, 1,
                1, 1,-1,
                1, 1, 0,
                1, 1, 1;
   REAL dx,dy,dz;
   CrystVector_REAL x(mvpAtom.size()),y(mvpAtom.size()),z(mvpAtom.size());
   CrystVector_REAL xSave,ySave,zSave;
   const int nbSymmetrics=vXYZCoords[0].rows();
   unsigned int ct=0;
   for(int i=0;i<nbSymmetrics;i++)
   {
      VFN_DEBUG_ENTRY("Molecule::POVRayDescription():Symmetric#"<<i,3)
      for(unsigned int j=0;j<mvpAtom.size();j++)
      {
         x(j)=vXYZCoords[j](i,0);
         y(j)=vXYZCoords[j](i,1);
         z(j)=vXYZCoords[j](i,2);
      }
      //Bring back central atom in unit cell; move peripheral atoms with the same amount
         dx=x(0);
         dy=y(0);
         dz=z(0);
         x(0) = fmod((float) x(0),(float)1); if(x(0)<0) x(0)+=1.;
         y(0) = fmod((float) y(0),(float)1); if(y(0)<0) y(0)+=1.;
         z(0) = fmod((float) z(0),(float)1); if(z(0)<0) z(0)+=1.;
         dx = x(0)-dx;
         dy = y(0)-dy;
         dz = z(0)-dz;
         for(unsigned int j=1;j<mvpAtom.size();j++)
         {
            x(j) += dx;
            y(j) += dy;
            z(j) += dz;
         }
      //Generate also translated atoms near the unit cell
      xSave=x;
      ySave=y;
      zSave=z;
      for(int j=0;j<translate.rows();j++)
      {
         x += translate(j,0);
         y += translate(j,1);
         z += translate(j,2);
         const REAL tmpxc=x.sum()/(REAL)x.numElements();
         const REAL tmpyc=y.sum()/(REAL)y.numElements();
         const REAL tmpzc=z.sum()/(REAL)z.numElements();
         if(   (tmpxc>xMin) && (tmpxc<xMax)
             &&(tmpyc>yMin) && (tmpyc<yMax)
             &&(tmpzc>zMin) && (tmpzc<zMax))
         {
            os<<"  //Symetric#"<<++ct<<endl;
            for(unsigned int k=0;k<mvpAtom.size();k++)
            {
               this->GetCrystal().FractionalToOrthonormalCoords(x(k),y(k),z(k));
               if(mvpAtom[k]->IsDummy()) continue;
               const float r=mvpAtom[k]->GetScatteringPower().GetColourRGB()[0];
               const float g=mvpAtom[k]->GetScatteringPower().GetColourRGB()[1];
               const float b=mvpAtom[k]->GetScatteringPower().GetColourRGB()[2];
               if(options.mShowLabel)
               {
                  /*
                  GLfloat colourChar [] = {1.0, 1.0, 1.0, 1.0}; 
                  if((r>0.8)&&(g>0.8)&&(b>0.8))
                  {
                     colourChar[0] = 0.5;
                     colourChar[1] = 0.5;
                     colourChar[2] = 0.5;
                  }
                  glMaterialfv(GL_FRONT, GL_AMBIENT,  colour0); 
                  glMaterialfv(GL_FRONT, GL_DIFFUSE,  colour0); 
                  glMaterialfv(GL_FRONT, GL_SPECULAR, colour0); 
                  glMaterialfv(GL_FRONT, GL_EMISSION, colourChar); 
                  glMaterialfv(GL_FRONT, GL_SHININESS,colour0);
                  glRasterPos3f(x(k)*en, y(k), z(k));
                  crystGLPrint(mvpAtom[k]->GetName());
                  */
               }
               os << "    ObjCrystAtom("
                  <<x(k)<<","
                  <<y(k)<<","
                  <<z(k)<<","
                  <<mvpAtom[k]->GetScatteringPower().GetRadius()<<","
                  <<"colour_"+mvpAtom[k]->GetScatteringPower().GetName()
                  <<")"<<endl;
            }
            for(unsigned int k=0;k<mvpBond.size();k++)
            {
               if(  (mvpBond[k]->GetAtom1().IsDummy())
                  ||(mvpBond[k]->GetAtom2().IsDummy()) ) continue;
               unsigned long n1,n2;
               //:KLUDGE: Get the atoms
               for(n1=0;n1<mvpAtom.size();n1++)
                  if(mvpAtom[n1]==&(mvpBond[k]->GetAtom1())) break;
               for(n2=0;n2<mvpAtom.size();n2++)
                  if(mvpAtom[n2]==&(mvpBond[k]->GetAtom2())) break;
               os << "    ObjCrystBond("
                  <<x(n1)<<","<<y(n1)<<","<<z(n1)<< ","
                  <<x(n2)<<","<<y(n2)<<","<<z(n2)<< ","
                  << "0.1,";
               if(mvpBond[k]->IsFreeTorsion()) os<<"colour_freebond)"<<endl;
               else os<<"colour_nonfreebond)"<<endl;

            }
         }//if in limits
         x=xSave;
         y=ySave;
         z=zSave;
      }//for translation
      VFN_DEBUG_EXIT("Molecule::POVRayDescription():Symmetric#"<<i,3)
   }//for symmetrics
   VFN_DEBUG_EXIT("Molecule::POVRayDescription()",3)
   return os;
}

void Molecule::GLInitDisplayList(const bool onlyIndependentAtoms,
                               const REAL xMin,const REAL xMax,
                               const REAL yMin,const REAL yMax,
                               const REAL zMin,const REAL zMax,
                               const bool displayEnantiomer,
                               const bool displayNames)const
{
   #ifdef OBJCRYST_GL
   VFN_DEBUG_ENTRY("Molecule::GLInitDisplayList()",3)
   if(mvpAtom.size()==0)
   {
      VFN_DEBUG_EXIT("Molecule::GLInitDisplayList():No atom to display !",4)
      return;
   }
   bool large=false;
   if(mvpAtom.size()>200) large=true;
   REAL en=1;
   if(displayEnantiomer==true) en=-1;
   this->UpdateScattCompList();
   //this->BuildRingList();
   //this->BuildStretchModeBondLength();
   //this->BuildStretchModeBondAngle();
   //this->BuildStretchModeTorsion();
   
   const GLfloat colour_bondnonfree[]= { 0.2, .2, .2, 1.0 };
   const GLfloat colour_bondrigid[]= { 0.5, .3, .3, 1.0 };
   const GLfloat colour_bondfree[]= { 0.8, .8, .8, 1.0 };
   const GLfloat colour0[] = {0.0f, 0.0f, 0.0f, 0.0f}; 
   
   GLUquadricObj* pQuadric = gluNewQuadric();
   
   if(true==onlyIndependentAtoms)//
   {
      REAL xc=mXYZ(0),yc=mXYZ(1),zc=mXYZ(2);
      this->GetCrystal().FractionalToOrthonormalCoords(xc,yc,zc);
      vector<MolAtom*>::const_iterator pos;
      for(pos=mvpAtom.begin();pos!=mvpAtom.end();pos++)
      {
         
         if((*pos)->IsDummy())continue;
         const float r=(*pos)->GetScatteringPower().GetColourRGB()[0];
         const float g=(*pos)->GetScatteringPower().GetColourRGB()[1];
         const float b=(*pos)->GetScatteringPower().GetColourRGB()[2];
         glPushMatrix();
            if(displayNames)
            {
               GLfloat colourChar [] = {1.0, 1.0, 1.0, 1.0}; 
               GLfloat colourCharRing [] = {1.0, 1.0, 0.8, 1.0}; 
               if((r>0.8)&&(g>0.8)&&(b>0.8))
               {
                  colourChar[0] = 0.5;
                  colourChar[1] = 0.5;
                  colourChar[2] = 0.5;
               }
               glMaterialfv(GL_FRONT, GL_AMBIENT,   colour0); 
               glMaterialfv(GL_FRONT, GL_DIFFUSE,   colour0); 
               glMaterialfv(GL_FRONT, GL_SPECULAR,  colour0);
               if((*pos)->IsInRing())
                  glMaterialfv(GL_FRONT, GL_EMISSION,  colourCharRing); 
               else
                  glMaterialfv(GL_FRONT, GL_EMISSION,  colourChar); 
               glMaterialfv(GL_FRONT, GL_SHININESS, colour0);
               glTranslatef((*pos)->X()*en+xc, (*pos)->Y()+yc, (*pos)->Z()+zc);
               crystGLPrint((*pos)->GetName());
            }
            else
            {
               const GLfloat colourAtom [] = {r, g, b, 1.0}; 
               glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE,   colourAtom); 
               glMaterialfv(GL_FRONT, GL_SPECULAR,              colour0); 
               glMaterialfv(GL_FRONT, GL_EMISSION,              colour0); 
               glMaterialfv(GL_FRONT, GL_SHININESS,             colour0);
               glPolygonMode(GL_FRONT, GL_FILL);
               glTranslatef((*pos)->X()*en+xc, (*pos)->Y()+yc, (*pos)->Z()+zc);
               gluSphere(pQuadric,(*pos)->GetScatteringPower().GetRadius()/3.,20,20);
            }
         glPopMatrix();
      }
   }//Only independent atoms ?
   else
   {
      VFN_DEBUG_ENTRY("Molecule::GLInitDisplayList():Show all symmetrics",3)
      // Reverse index of atoms
      map<const MolAtom*,unsigned long> rix;
         {
            long i=0;
            for(vector<MolAtom*>::const_iterator pos=mvpAtom.begin();pos!=mvpAtom.end();++pos)
               rix[*pos]=i++;
         }
      vector<CrystMatrix_REAL> vXYZCoords;
      {
         this->GetScatteringComponentList();
         REAL x0,y0,z0;
         for(long i=0;i<mScattCompList.GetNbComponent();++i)
         {
            x0=mScattCompList(i).mX;
            y0=mScattCompList(i).mY;
            z0=mScattCompList(i).mZ;
            vXYZCoords.push_back(this->GetCrystal().GetSpaceGroup().
                           GetAllSymmetrics(x0,y0,z0,false,false,false));
         }
      }
      CrystMatrix_int translate(27,3);
      translate=  -1,-1,-1,
                  -1,-1, 0,
                  -1,-1, 1,
                  -1, 0,-1,
                  -1, 0, 0,
                  -1, 0, 1,
                  -1, 1,-1,
                  -1, 1, 0,
                  -1, 1, 1,
                   0,-1,-1,
                   0,-1, 0,
                   0,-1, 1,
                   0, 0,-1,
                   0, 0, 0,
                   0, 0, 1,
                   0, 1,-1,
                   0, 1, 0,
                   0, 1, 1,
                   1,-1,-1,
                   1,-1, 0,
                   1,-1, 1,
                   1, 0,-1,
                   1, 0, 0,
                   1, 0, 1,
                   1, 1,-1,
                   1, 1, 0,
                   1, 1, 1;
      REAL dx,dy,dz;
      CrystVector_REAL x(mvpAtom.size()),y(mvpAtom.size()),z(mvpAtom.size());
      CrystVector_REAL xSave,ySave,zSave;
      const int nbSymmetrics=vXYZCoords[0].rows();
      for(int i=0;i<nbSymmetrics;i++)
      {
         VFN_DEBUG_ENTRY("Molecule::GLInitDisplayList():Symmetric#"<<i,3)
         for(unsigned int j=0;j<mvpAtom.size();j++)
         {
            x(j)=vXYZCoords[j](i,0);
            y(j)=vXYZCoords[j](i,1);
            z(j)=vXYZCoords[j](i,2);
         }
         //Bring back central atom in unit cell; move peripheral atoms with the same amount
            dx=x(0);
            dy=y(0);
            dz=z(0);
            x(0) = fmod((float) x(0),(float)1); if(x(0)<0) x(0)+=1.;
            y(0) = fmod((float) y(0),(float)1); if(y(0)<0) y(0)+=1.;
            z(0) = fmod((float) z(0),(float)1); if(z(0)<0) z(0)+=1.;
            dx = x(0)-dx;
            dy = y(0)-dy;
            dz = z(0)-dz;
            for(unsigned int j=1;j<mvpAtom.size();j++)
            {
               x(j) += dx;
               y(j) += dy;
               z(j) += dz;
            }
         //Generate also translated atoms near the unit cell
         xSave=x;
         ySave=y;
         zSave=z;
         for(int j=0;j<translate.rows();j++)
         {
            x += translate(j,0);
            y += translate(j,1);
            z += translate(j,2);
            const REAL tmpxc=x.sum()/(REAL)x.numElements();
            const REAL tmpyc=y.sum()/(REAL)y.numElements();
            const REAL tmpzc=z.sum()/(REAL)z.numElements();
            if(   (tmpxc>xMin) && (tmpxc<xMax)
                &&(tmpyc>yMin) && (tmpyc<yMax)
                &&(tmpzc>zMin) && (tmpzc<zMax))
            {
               for(unsigned int k=0;k<mvpAtom.size();k++)
               {
                  this->GetCrystal().FractionalToOrthonormalCoords(x(k),y(k),z(k));
                  if(mvpAtom[k]->IsDummy()) continue;
                  glPushMatrix();
                     const float r=mvpAtom[k]->GetScatteringPower().GetColourRGB()[0];
                     const float g=mvpAtom[k]->GetScatteringPower().GetColourRGB()[1];
                     const float b=mvpAtom[k]->GetScatteringPower().GetColourRGB()[2];
                     if(displayNames)
                     {
                        GLfloat colourChar [] = {1.0, 1.0, 1.0, 1.0}; 
                        GLfloat colourCharRing [] = {1.0, 1.0, 0.8, 1.0}; 
                        if((r>0.8)&&(g>0.8)&&(b>0.8))
                        {
                           colourChar[0] = 0.5;
                           colourChar[1] = 0.5;
                           colourChar[2] = 0.5;
                        }
                        glMaterialfv(GL_FRONT, GL_AMBIENT,  colour0); 
                        glMaterialfv(GL_FRONT, GL_DIFFUSE,  colour0); 
                        glMaterialfv(GL_FRONT, GL_SPECULAR, colour0); 
                        if(mvpAtom[k]->IsInRing())
                           glMaterialfv(GL_FRONT, GL_EMISSION,  colourCharRing); 
                        else
                           glMaterialfv(GL_FRONT, GL_EMISSION,  colourChar); 
                        glMaterialfv(GL_FRONT, GL_SHININESS,colour0);
                        glRasterPos3f(x(k)*en, y(k), z(k));
                        crystGLPrint(mvpAtom[k]->GetName());
                     }
                     else
                     {
                        if(!large)
                        {
                           const GLfloat colourAtom [] = {r, g, b, 1.0}; 
                           glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE,colourAtom); 
                           glMaterialfv(GL_FRONT, GL_SPECULAR,           colour0); 
                           glMaterialfv(GL_FRONT, GL_EMISSION,           colour0); 
                           glMaterialfv(GL_FRONT, GL_SHININESS,          colour0);
                           glPolygonMode(GL_FRONT, GL_FILL);
                           glTranslatef(x(k)*en, y(k), z(k));
                           gluSphere(pQuadric,
                              mvpAtom[k]->GetScatteringPower().GetRadius()/3.,20,20);
                        }
                        else
                        {
                           const GLfloat colourAtom [] = {r, g, b, 1.0}; 
                           glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE,colourAtom); 
                           glMaterialfv(GL_FRONT, GL_SPECULAR,           colour0); 
                           glMaterialfv(GL_FRONT, GL_EMISSION,           colour0); 
                           glMaterialfv(GL_FRONT, GL_SHININESS,          colour0);
                           glPolygonMode(GL_FRONT, GL_FILL);
                           glTranslatef(x(k)*en, y(k), z(k));
                           gluSphere(pQuadric,.15,10,10);
                        }
                     }
                  glPopMatrix();
               }
               if(displayNames==false)
               {
                  if(large)
                  {
                     glMaterialfv(GL_FRONT, GL_SPECULAR,  colour0); 
                     glMaterialfv(GL_FRONT, GL_EMISSION,  colour0); 
                     glMaterialfv(GL_FRONT, GL_SHININESS, colour0);
                     //glPolygonMode(GL_FRONT, GL_FILL);
                     for(unsigned int k=0;k<mvpBond.size();k++)
                     {
                        if(  (mvpBond[k]->GetAtom1().IsDummy())
                           ||(mvpBond[k]->GetAtom2().IsDummy()) ) continue;
                        const float r1=mvpBond[k]->GetAtom1().GetScatteringPower().GetColourRGB()[0];
                        const float g1=mvpBond[k]->GetAtom1().GetScatteringPower().GetColourRGB()[1];
                        const float b1=mvpBond[k]->GetAtom1().GetScatteringPower().GetColourRGB()[2];
                        const float r2=mvpBond[k]->GetAtom2().GetScatteringPower().GetColourRGB()[0];
                        const float g2=mvpBond[k]->GetAtom2().GetScatteringPower().GetColourRGB()[1];
                        const float b2=mvpBond[k]->GetAtom2().GetScatteringPower().GetColourRGB()[2];
                        const GLfloat colourAtom1 [] = {r1, g1, b1, 1.0};
                        const GLfloat colourAtom2 [] = {r2, g2, b2, 1.0};
                        const unsigned long n1=rix[&(mvpBond[k]->GetAtom1())],
                                            n2=rix[&(mvpBond[k]->GetAtom2())];
                        #if 0
                        glPushMatrix();
                           glBegin(GL_LINE_STRIP);
                              glMaterialfv(GL_FRONT, GL_SPECULAR,  colourAtom1); 
                              glMaterialfv(GL_FRONT, GL_EMISSION,  colourAtom1); 
                              glMaterialfv(GL_FRONT, GL_SHININESS, colourAtom1);
                              glVertex3f(x(n1)*en,y(n1),z(n1));
                              glVertex3f((x(n1)+x(n2))/2*en,(y(n1)+y(n2))/2,(z(n1)+z(n2))/2);
                              glMaterialfv(GL_FRONT, GL_SPECULAR,  colourAtom2); 
                              glMaterialfv(GL_FRONT, GL_EMISSION,  colourAtom2); 
                              glMaterialfv(GL_FRONT, GL_SHININESS, colourAtom2);
                              glVertex3f(x(n2)*en,y(n2),z(n2));
                           glEnd();
                        glPopMatrix();
                        #else
                        const REAL height= sqrt(abs(  (x(n2)-x(n1))*(x(n2)-x(n1))
                                                     +(y(n2)-y(n1))*(y(n2)-y(n1))
                                                     +(z(n2)-z(n1))*(z(n2)-z(n1))));
                        glMaterialfv(GL_FRONT, GL_SPECULAR,  colour0); 
                        glMaterialfv(GL_FRONT, GL_EMISSION,  colour0); 
                        glMaterialfv(GL_FRONT, GL_SHININESS, colour0);
                        glPolygonMode(GL_FRONT, GL_FILL);
                        glPushMatrix();
                           glTranslatef(x(n1)*en, y(n1), z(n1));
                           GLUquadricObj *quadobj = gluNewQuadric();
                           glRotatef(180,(x(n2)-x(n1))*en,y(n2)-y(n1),z(n2)-z(n1)+height);// ?!?!?!
                           
                           glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,colourAtom1);
                           gluCylinder(quadobj,.1,.1,height/2,10,1 );
                           gluDeleteQuadric(quadobj);
                           
                           glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,colourAtom2);
                           GLUquadricObj *quadobj2 = gluNewQuadric();
                           glTranslatef(0, 0, height/2);
                           gluCylinder(quadobj2,.1,.1,height/2,10,1 );
                           gluDeleteQuadric(quadobj2);
                        glPopMatrix();
                        #endif
                     }
                  }
                  else
                  {
                     for(unsigned int k=0;k<mvpBond.size();k++)
                     {
                        if(  (mvpBond[k]->GetAtom1().IsDummy())
                           ||(mvpBond[k]->GetAtom2().IsDummy()) ) continue;
                        const unsigned long n1=rix[&(mvpBond[k]->GetAtom1())],
                                            n2=rix[&(mvpBond[k]->GetAtom2())];
                        // Is the bond in a rigid group ?
                        bool isRigidGroup=false;
                        for(vector<RigidGroup *>::const_iterator pos=this->GetRigidGroupList().begin();pos!=this->GetRigidGroupList().end();++pos)
                        {
                           if(  ((*pos)->find(&(mvpBond[k]->GetAtom1()))!=(*pos)->end()) 
                              &&((*pos)->find(&(mvpBond[k]->GetAtom2()))!=(*pos)->end()) )
                           {
                              isRigidGroup=true;
                              break;
                           }
                        }
                        if(isRigidGroup)
                           glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,colour_bondrigid);
                        else
                        {
                           if(mvpBond[k]->IsFreeTorsion())
                              glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,colour_bondfree);
                           else
                              glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,colour_bondnonfree);
                        }
                        glMaterialfv(GL_FRONT, GL_SPECULAR,  colour0); 
                        glMaterialfv(GL_FRONT, GL_EMISSION,  colour0); 
                        glMaterialfv(GL_FRONT, GL_SHININESS, colour0);
                        glPolygonMode(GL_FRONT, GL_FILL);
                        glPushMatrix();
                           glTranslatef(x(n1)*en, y(n1), z(n1));
                           GLUquadricObj *quadobj = gluNewQuadric();
                           //glColor4f(1.0f,1.0f,1.0f,1.0);
                           const REAL height= sqrt(abs(  (x(n2)-x(n1))*(x(n2)-x(n1))
                                                        +(y(n2)-y(n1))*(y(n2)-y(n1))
                                                        +(z(n2)-z(n1))*(z(n2)-z(n1))));
                           glRotatef(180,(x(n2)-x(n1))*en,y(n2)-y(n1),z(n2)-z(n1)+height);// ?!?!?!
                           gluCylinder(quadobj,.1,.1,height,10,1 );
                           gluDeleteQuadric(quadobj);
                        glPopMatrix();
                     }
                  }
               }
            }//if in limits
            x=xSave;
            y=ySave;
            z=zSave;
         }//for translation
         VFN_DEBUG_EXIT("Molecule::GLInitDisplayList():Symmetric#"<<i,3)
      }//for symmetrics
      VFN_DEBUG_EXIT("Molecule::GLInitDisplayList():Show all symmetrics",3)
   }//else
   gluDeleteQuadric(pQuadric);
   VFN_DEBUG_EXIT("Molecule::GLInitDisplayList()",3)
   #endif //GLCryst
}

void Molecule::AddAtom(const REAL x, const REAL y, const REAL z,
                       const ScatteringPower *pPow, const string &name,
                       const bool updateDisplay)
{
   VFN_DEBUG_ENTRY("Molecule::AddAtom():"<<name,5)
   string thename=name;
   if(thename==string(""))
   {// This should not be needed, the atom will reset the parameters name when its name is set
      char buf[100];
      if(pPow!=0) sprintf(buf,"%s_%s_%d",this->GetName().c_str(),pPow->GetName().c_str(),mvpAtom.size()+1);
      else sprintf(buf,"%s_X_%d",this->GetName().c_str(),mvpAtom.size()+1);
      thename=buf;
   }
   mvpAtom.push_back(new MolAtom(x,y,z,pPow,thename,*this));
   mClockAtomPosition.Click();
   mClockAtomScattPow.Click();
   ++mScattCompList;
   {
      RefinablePar tmp(thename+"_x",&(mvpAtom.back()->X()),0.,1.,
                        gpRefParTypeScattConformX,
                        REFPAR_DERIV_STEP_ABSOLUTE,false,false,true,false,1.,1.);
      tmp.AssignClock(mClockAtomPosition);
      tmp.SetGlobalOptimStep(0.05);
      tmp.SetDerivStep(1e-4);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp(thename+"_y",&(mvpAtom.back()->Y()),0.,1.,
                        gpRefParTypeScattConformY,
                        REFPAR_DERIV_STEP_ABSOLUTE,false,false,true,false,1.,1.);
      tmp.AssignClock(mClockAtomPosition);
      tmp.SetGlobalOptimStep(0.05);
      tmp.SetDerivStep(1e-4);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp(thename+"_z",&(mvpAtom.back()->Z()),0.,1.,
                        gpRefParTypeScattConformZ,
                        REFPAR_DERIV_STEP_ABSOLUTE,false,false,true,false,1.,1.);
      tmp.AssignClock(mClockAtomPosition);
      tmp.SetGlobalOptimStep(0.05);
      tmp.SetDerivStep(1e-4);
      this->AddPar(tmp);
   }
   
   mClockScatterer.Click();
   if(updateDisplay) this->UpdateDisplay();
   VFN_DEBUG_EXIT("Molecule::AddAtom()",5)
}

vector<MolAtom*>::iterator Molecule::RemoveAtom(MolAtom &atom, const bool del)
{
   VFN_DEBUG_ENTRY("Molecule::RemoveAtom():"<<atom.GetName(),6)
   vector<MolAtom*>::iterator pos=find(mvpAtom.begin(),mvpAtom.end(),&atom);
   if(pos==mvpAtom.end())
   {
      throw ObjCrystException("Molecule::RemoveAtom():"+atom.GetName()
                        +" is not in this Molecule:"+this->GetName());
   }
   // Delete parameters
      this->RemovePar(&(this->GetPar(&(atom.X()))));
      this->RemovePar(&(this->GetPar(&(atom.Y()))));
      this->RemovePar(&(this->GetPar(&(atom.Z()))));
   // Delete relevant bonds, bond angles, dihedral angles...
   for(vector<MolBond*>::iterator posb=mvpBond.begin();posb!=mvpBond.end();)
   {
      if( (&atom==&((*posb)->GetAtom1())) || (&atom==&((*posb)->GetAtom2())) )
      {
         posb=this->RemoveBond(**posb, del);
      }
      else ++posb;
   }
   for(vector<MolBondAngle*>::iterator posb=mvpBondAngle.begin();posb!=mvpBondAngle.end();)
   {
      if(  (&atom==&((*posb)->GetAtom1())) || (&atom==&((*posb)->GetAtom2()))
         ||(&atom==&((*posb)->GetAtom3())))
      {
         posb=this->RemoveBondAngle(**posb, del);
      }
      else ++posb;
   }
   for(vector<MolDihedralAngle*>::iterator posb=mvpDihedralAngle.begin();
       posb!=mvpDihedralAngle.end();)
   {
      if(  (&atom==&((*posb)->GetAtom1())) || (&atom==&((*posb)->GetAtom2()))
         ||(&atom==&((*posb)->GetAtom3())) || (&atom==&((*posb)->GetAtom4())))
         posb=this->RemoveDihedralAngle(**posb, del);
      else ++posb;
  }
   mClockAtomList.Click();
   mClockScatterer.Click();
   
   if(mpCenterAtom==*pos) mpCenterAtom=0;
   if(del) delete *pos;
   pos=mvpAtom.erase(pos);
   --mScattCompList;

   this->UpdateDisplay();
   VFN_DEBUG_EXIT("Molecule::RemoveAtom()",6)
   return pos;
}

void Molecule::AddBond(MolAtom &atom1, MolAtom &atom2,
                       const REAL length, const REAL sigma, const REAL delta,
                       const REAL bondOrder,
                       const bool updateDisplay)
{
   VFN_DEBUG_ENTRY("Molecule::AddBond()",5)
   mvpBond.push_back(new MolBond(atom1,atom2,length,sigma,delta,*this,bondOrder));
   this->AddRestraint(mvpBond.back());
   mClockBondList.Click();
   if(updateDisplay) this->UpdateDisplay();
   VFN_DEBUG_EXIT("Molecule::AddBond()",5)
}

vector<MolBond*>::iterator Molecule::RemoveBond(const MolBond &bond, const bool del)
{
   VFN_DEBUG_ENTRY("Molecule::RemoveBond():"<<bond.GetName(),6)
   vector<MolBond*>::iterator pos=find(mvpBond.begin(),mvpBond.end(),&bond);
   if(pos==mvpBond.end())
   {
      throw ObjCrystException("Molecule::RemoveBond():"+bond.GetAtom1().GetName()
                              +"-"+bond.GetAtom2().GetName()
                              +" is not in this Molecule:"+this->GetName());
   }
   this->RemoveRestraint(*pos);
   mClockBondList.Click();
   if(del) delete *pos;
   pos= mvpBond.erase(pos);
   this->UpdateDisplay();
   VFN_DEBUG_EXIT("Molecule::RemoveBond():",6)
   return pos;
}

vector<MolBond*>::const_iterator Molecule::FindBond(const MolAtom &at1,const MolAtom &at2)const
{
   for(vector<MolBond*>::const_iterator pos=mvpBond.begin();pos!=mvpBond.end();++pos)
   {
      if(  ((&((*pos)->GetAtom1())==&at1)&&(&((*pos)->GetAtom2())==&at2))
         ||((&((*pos)->GetAtom1())==&at2)&&(&((*pos)->GetAtom2())==&at1)))
         return pos;
   }
   return mvpBond.end();
}
vector<MolBond*>::iterator Molecule::FindBond(const MolAtom &at1,const MolAtom &at2)
{
   for(vector<MolBond*>::iterator pos=mvpBond.begin();pos!=mvpBond.end();++pos)
   {
      if(  ((&((*pos)->GetAtom1())==&at1)&&(&((*pos)->GetAtom2())==&at2))
         ||((&((*pos)->GetAtom1())==&at2)&&(&((*pos)->GetAtom2())==&at1)))
         return pos;
   }
   return mvpBond.end();
}

void Molecule::AddBondAngle(MolAtom &atom1, MolAtom &atom2, MolAtom &atom3,
                            const REAL angle, const REAL sigma, const REAL delta,
                            const bool updateDisplay)
{
   VFN_DEBUG_ENTRY("Molecule::AddBondAngle()",5)
   mvpBondAngle.push_back(new MolBondAngle(atom1,atom2,atom3,angle,sigma,delta,*this));
   this->AddRestraint(mvpBondAngle.back());
   mClockBondAngleList.Click();
   if(updateDisplay) this->UpdateDisplay();
   VFN_DEBUG_EXIT("Molecule::AddBondAngle()",5)
}

vector<MolBondAngle*>::iterator Molecule::RemoveBondAngle(const MolBondAngle &angle, const bool del)
{
   VFN_DEBUG_ENTRY("Molecule::RemoveBondAngle():"<<angle.GetName(),6)
   vector<MolBondAngle*>::iterator pos=find(mvpBondAngle.begin(),mvpBondAngle.end(),&angle);
   if(pos==mvpBondAngle.end())
   {
      throw ObjCrystException("Molecule::RemoveBondAngle():"+angle.GetAtom1().GetName()
                              +"-"+angle.GetAtom2().GetName()+"-"+angle.GetAtom3().GetName()
                              +" is not in this Molecule:"+this->GetName());
   }
   this->RemoveRestraint(*pos);
   mClockBondAngleList.Click();
   if(del) delete *pos;
   pos=mvpBondAngle.erase(pos);
   this->UpdateDisplay();
   VFN_DEBUG_EXIT("Molecule::RemoveBondAngle():",6)
   return pos;
}

vector<MolBondAngle*>::const_iterator Molecule::FindBondAngle(const MolAtom &at1,
                                                              const MolAtom &at2,
                                                              const MolAtom &at3)const
{
   for(vector<MolBondAngle*>::const_iterator pos=mvpBondAngle.begin();
       pos!=mvpBondAngle.end();++pos)
   {
      if(  (&((*pos)->GetAtom2())==&at2)
         &&(  ((&((*pos)->GetAtom1())==&at1)&&(&((*pos)->GetAtom3())==&at3))
            ||((&((*pos)->GetAtom1())==&at3)&&(&((*pos)->GetAtom3())==&at1))))
         return pos;
   }
   return mvpBondAngle.end();
}

void Molecule::AddDihedralAngle(MolAtom &atom1, MolAtom &atom2,
                                MolAtom &atom3, MolAtom &atom4,
                                const REAL angle, const REAL sigma, const REAL delta,
                                const bool updateDisplay)
{
   VFN_DEBUG_ENTRY("Molecule::AddDihedralAngle()",5)
   mvpDihedralAngle.push_back(new MolDihedralAngle(atom1,atom2,atom3,atom4,
                                                   angle,sigma,delta,*this));
   this->AddRestraint(mvpDihedralAngle.back());
   mClockDihedralAngleList.Click();
   if(updateDisplay) this->UpdateDisplay();
   VFN_DEBUG_EXIT("Molecule::AddDihedralAngle()",5)
}

vector<MolDihedralAngle*>::iterator Molecule::RemoveDihedralAngle(const MolDihedralAngle &angle, const bool del)
{
   VFN_DEBUG_ENTRY("Molecule::RemoveDihedralAngle():"<<angle.GetName(),6)
   vector<MolDihedralAngle*>::iterator pos=find(mvpDihedralAngle.begin(),
                                                mvpDihedralAngle.end(),&angle);
   if(pos==mvpDihedralAngle.end())
   {
      throw ObjCrystException("Molecule::RemoveDihedralAngle():"+angle.GetAtom1().GetName()
                              +"-"+angle.GetAtom2().GetName()+"-"+angle.GetAtom3().GetName()
                              +"-"+angle.GetAtom4().GetName()
                              +" is not in this Molecule:"+this->GetName());
   }
   this->RemoveRestraint(*pos);
   mClockDihedralAngleList.Click();
   if(del) delete *pos;
   pos=mvpDihedralAngle.erase(pos);
   this->UpdateDisplay();
   VFN_DEBUG_ENTRY("Molecule::RemoveDihedralAngle():",6)
   return pos;
}
vector<MolDihedralAngle*>::const_iterator Molecule::FindDihedralAngle(const MolAtom &at1,
                                                                      const MolAtom &at2,
                                                                      const MolAtom &at3,
                                                                      const MolAtom &at4)const
{
   for(vector<MolDihedralAngle*>::const_iterator pos=mvpDihedralAngle.begin();
       pos!=mvpDihedralAngle.end();++pos)
   {
      if(  (  ((&((*pos)->GetAtom1())==&at1)&&(&((*pos)->GetAtom2())==&at2))
            &&((&((*pos)->GetAtom3())==&at3)&&(&((*pos)->GetAtom4())==&at4)))
         ||(  ((&((*pos)->GetAtom4())==&at1)&&(&((*pos)->GetAtom3())==&at2))
            &&((&((*pos)->GetAtom2())==&at3)&&(&((*pos)->GetAtom1())==&at4))))
         return pos;
   }
   return mvpDihedralAngle.end();
}

void Molecule::AddRigidGroup(const RigidGroup &group,
                             const bool updateDisplay)
{
  mvRigidGroup.push_back(new RigidGroup(group));
  #ifdef RIGID_BODY_STRICT_EXPERIMENTAL
  char buf[50];
  const unsigned int i=mvRigidGroup.size();
  RigidGroup* p=this->GetRigidGroupList().back();
  p->mX=0;
  p->mY=0;
  p->mZ=0;
  p->mQuat.Q0()=1;
  p->mQuat.Q1()=0;
  p->mQuat.Q2()=0;
  p->mQuat.Q3()=0;
  {
    sprintf(buf,"RigidGroup%d_x",i);
    RefinablePar tmp(buf,&(p->mX),0.,1.,
                    gpRefParTypeScattConformX,
                    REFPAR_DERIV_STEP_ABSOLUTE,false,false,true,false,1.,1.);
    tmp.AssignClock(mClockAtomPosition);
    tmp.SetGlobalOptimStep(0.05);
    tmp.SetDerivStep(1e-4);
    this->AddPar(tmp);
  }
  {
    sprintf(buf,"RigidGroup%d_y",i);
    RefinablePar tmp(buf,&(p->mY),0.,1.,
                    gpRefParTypeScattConformY,
                    REFPAR_DERIV_STEP_ABSOLUTE,false,false,true,false,1.,1.);
    tmp.AssignClock(mClockAtomPosition);
    tmp.SetGlobalOptimStep(0.05);
    tmp.SetDerivStep(1e-4);
    this->AddPar(tmp);
  }
  {
    sprintf(buf,"RigidGroup%d_z",i);
    RefinablePar tmp(buf,&(p->mZ),0.,1.,
                    gpRefParTypeScattConformZ,
                    REFPAR_DERIV_STEP_ABSOLUTE,false,false,true,false,1.,1.);
    tmp.AssignClock(mClockAtomPosition);
    tmp.SetGlobalOptimStep(0.05);
    tmp.SetDerivStep(1e-4);
    this->AddPar(tmp);
  }
  {
    sprintf(buf,"RigidGroup%d_Q1",i);
    RefinablePar tmp(buf,&(p->mQuat.Q1()),-1,1.,
                    gpRefParTypeScattConform,
                    REFPAR_DERIV_STEP_ABSOLUTE,true,false,true,false,1.,1.);
    tmp.AssignClock(mClockAtomPosition);
    tmp.SetGlobalOptimStep(0.01);
    tmp.SetDerivStep(1e-4);
    this->AddPar(tmp);
  }
  {
    sprintf(buf,"RigidGroup%d_Q2",i);
    RefinablePar tmp(buf,&(p->mQuat.Q2()),-1,1.,
                    gpRefParTypeScattConform,
                    REFPAR_DERIV_STEP_ABSOLUTE,true,false,true,false,1.,1.);
    tmp.AssignClock(mClockAtomPosition);
    tmp.SetGlobalOptimStep(0.01);
    tmp.SetDerivStep(1e-4);
    this->AddPar(tmp);
  }
  {
    sprintf(buf,"RigidGroup%d_Q3",i);
    RefinablePar tmp(buf,&(p->mQuat.Q3()),-1,1.,
                    gpRefParTypeScattConform,
                    REFPAR_DERIV_STEP_ABSOLUTE,true,false,true,false,1.,1.);
    tmp.AssignClock(mClockAtomPosition);
    tmp.SetGlobalOptimStep(0.01);
    tmp.SetDerivStep(1e-4);
    this->AddPar(tmp);
  }
  #endif
  mClockRigidGroup.Click();
  if(updateDisplay) this->UpdateDisplay();
}

vector<RigidGroup*>::iterator Molecule::RemoveRigidGroup(const RigidGroup &g,const bool updateDisplay, const bool del)
{
   vector<RigidGroup*>::iterator pos=find(mvRigidGroup.begin(),mvRigidGroup.end(),&g);
   if(pos==mvRigidGroup.end()) return pos;
   #ifdef RIGID_BODY_STRICT_EXPERIMENTAL
   // Remove the refinable parameters (even if del==False - used for python delayed deletion)
   // NOTE - this should only be done outside an optimization, since rigid group translationnal 
   // and rotationnal parameters are resetted at the end of the optimization, and the atomic
   // parameters are directly the correct ones (thus deletion of the rigid group does not change
   // the final coordinates).
   this->RemovePar(&(this->GetPar(&((*pos)->mX))));
   this->RemovePar(&(this->GetPar(&((*pos)->mY))));
   this->RemovePar(&(this->GetPar(&((*pos)->mZ))));
   this->RemovePar(&(this->GetPar(&((*pos)->mQuat.Q1()))));
   this->RemovePar(&(this->GetPar(&((*pos)->mQuat.Q2()))));
   this->RemovePar(&(this->GetPar(&((*pos)->mQuat.Q3()))));
   #endif
   if(del) delete *pos;
   pos=mvRigidGroup.erase(pos);
   if(updateDisplay) this->UpdateDisplay();
   return pos;
}

MolAtom &Molecule::GetAtom(unsigned int i){return *mvpAtom[i];}

const MolAtom &Molecule::GetAtom(unsigned int i)const{return *mvpAtom[i];}

MolAtom &Molecule::GetAtom(const string &name){return **(this->FindAtom(name));}

const MolAtom &Molecule::GetAtom(const string &name)const{return **(this->FindAtom(name));}

void Molecule::OptimizeConformation(const long nbTrial,const REAL stopCost)
{
   VFN_DEBUG_ENTRY("Molecule::OptimizeConformation()",5)
   MonteCarloObj globalOptObj(true);
   globalOptObj.AddRefinableObj(*this);
   globalOptObj.SetAlgorithmParallTempering(ANNEALING_EXPONENTIAL,10000.,1.,
                                            ANNEALING_EXPONENTIAL,10,.1);      

   long nb=nbTrial;
   mIsSelfOptimizing=true;
   globalOptObj.Optimize(nb,false,stopCost);
   mIsSelfOptimizing=false;
   // Must rebuild Flip & Rotor group, in case they were tested with an absurd conformation
   mClockFlipGroup.Reset();
   mClockRotorGroup.Reset();
   VFN_DEBUG_EXIT("Molecule::OptimizeConformation()",5)
}

void Molecule::OptimizeConformationSteepestDescent(const REAL maxStep,const unsigned nbStep)
{
   //cout<<"LLK="<<this->GetLogLikelihood()<<endl;
   for(int i=0;i<nbStep;++i)
   {
      // Calc full gradient
      //Calculate & display gradient
      map<MolAtom*,XYZ> grad;
      for(vector<MolAtom*>::iterator pos=this->GetAtomList().begin();
         pos!=this->GetAtomList().end();++pos) grad[*pos]=XYZ(0,0,0);
      
      // :TODO: remove atoms that are in rigid groups ?
      for(vector<MolBond*>::iterator pos=this->GetBondList().begin();
         pos!=this->GetBondList().end();++pos)          (*pos)->CalcGradient(grad);
      
      for(vector<MolBondAngle*>::iterator pos=this->GetBondAngleList().begin();
         pos!=this->GetBondAngleList().end();++pos)     (*pos)->CalcGradient(grad);
      
      for(vector<MolDihedralAngle*>::iterator pos=this->GetDihedralAngleList().begin();
         pos!=this->GetDihedralAngleList().end();++pos) (*pos)->CalcGradient(grad);
      
      #if 0
      // Display gradient - for tests
      for(map<MolAtom*,XYZ>::const_iterator pos=grad.begin();pos!=grad.end();++pos)
      {
         char buf[100];
         sprintf(buf,"%10s Grad LLK= (%8.3f %8.3f %8.3f)",
               pos->first->GetName().c_str(),pos->second.x,pos->second.y,pos->second.z);
         cout<<buf<<endl;
      }
      #endif
      // Find maximum absolute value of gradient
      REAL f=0;
      for(map<MolAtom*,XYZ>::const_iterator pos=grad.begin();pos!=grad.end();++pos)
      {
         if(abs(pos->second.x)>f) f=abs(pos->second.x);
         if(abs(pos->second.y)>f) f=abs(pos->second.y);
         if(abs(pos->second.z)>f) f=abs(pos->second.z);
      }
      if(f>1e-6) f=maxStep/f;
      else break;//nothing to optimize ?
      // Average derivatives inside rigid groups
      //:TODO: still allow rotations of rigid groups ?
      //:TODO: Handle case when one atom belongs to several rigid groups...
      for(vector<RigidGroup *>::const_iterator pos=this->GetRigidGroupList().begin();pos!=this->GetRigidGroupList().end();++pos)
      {
         if((*pos)->size()==0) continue; // Just in case...
         REAL dx=0,dy=0,dz=0;
         for(set<MolAtom *>::const_iterator at=(*pos)->begin();at!=(*pos)->end();++at)
         {
            dx+=(*at)->GetX();
            dy+=(*at)->GetY();
            dz+=(*at)->GetZ();
         }
         dx/=(*pos)->size();
         dy/=(*pos)->size();
         dz/=(*pos)->size();
         for(set<MolAtom *>::const_iterator at=(*pos)->begin();at!=(*pos)->end();++at)
         {
            grad[*at].x=dx;
            grad[*at].y=dy;
            grad[*at].z=dz;
         }
      }
      // Move according to max step to minimize LLK
      for(map<MolAtom*,XYZ>::const_iterator pos=grad.begin();pos!=grad.end();++pos)
      {
         pos->first->SetX(pos->first->GetX()-pos->second.x*f);
         pos->first->SetY(pos->first->GetY()-pos->second.y*f);
         pos->first->SetZ(pos->first->GetZ()-pos->second.z*f);
      }
      //this->RestraintStatus(cout);
      //cout<<"LLK="<<this->GetLogLikelihood()<<endl;
   }
}

void Molecule::MolecularDynamicsEvolve(map<MolAtom*,XYZ> &v0,const unsigned nbStep,const REAL dt,
                                       const vector<MolBond*> &vb,const vector<MolBondAngle*> &va,
                                       const vector<MolDihedralAngle*> &vd,
                                       map<RigidGroup*,std::pair<XYZ,XYZ> > &vr, REAL nrj0)
{
   const vector<MolBond*> *pvb=&vb;
   const vector<MolBondAngle*> *pva=&va;
   const vector<MolDihedralAngle*> *pvd=&vd;
   map<RigidGroup*,std::pair<XYZ,XYZ> > *pvr=&vr;
   if((pvb->size()==0)&&(pva->size()==0)&&(pvd->size()==0))
   {
      pvb=&(this->GetBondList());
      pva=&(this->GetBondAngleList());
      pvd=&(this->GetDihedralAngleList());
      for(vector<RigidGroup *>::iterator pos=this->GetRigidGroupList().begin();pos!=this->GetRigidGroupList().end();++pos)
         (*pvr)[*pos]=make_pair(XYZ(0,0,0),XYZ(0,0,0));
   }
   const REAL m=500;// mass
   const REAL im=1./m;
   
   // Try to keep total energy constant
   REAL e_v,e_k,v_r=1.0;
   for(int i=0;i<nbStep;++i)
   {
      // Calc full gradient
      map<MolAtom*,XYZ> grad;
      for(map<MolAtom*,XYZ>::iterator pos=v0.begin();pos!=v0.end();++pos) grad[pos->first]=XYZ(0,0,0);
      
      // :TODO: handle rigid groups ?
      e_v=0;
      for(vector<MolBond*>::const_iterator pos=pvb->begin();pos!=pvb->end();++pos)
      {
         (*pos)->CalcGradient(grad);
         e_v+=(*pos)->GetLogLikelihood(false,false);
      }
      
      for(vector<MolBondAngle*>::const_iterator pos=pva->begin();pos!=pva->end();++pos)
      {
         (*pos)->CalcGradient(grad);
         e_v+=(*pos)->GetLogLikelihood(false,false);
      }
      
      for(vector<MolDihedralAngle*>::const_iterator pos=pvd->begin();pos!=pvd->end();++pos)
      {
         (*pos)->CalcGradient(grad);
         e_v+=(*pos)->GetLogLikelihood(false,false);
      }
      
      //kinetic energy
      e_k=0;
      for(map<MolAtom*,XYZ>::const_iterator pos=v0.begin();pos!=v0.end();++pos)
         e_k += 0.5*m*(pos->second.x*pos->second.x + pos->second.y*pos->second.y + pos->second.z*pos->second.z);

      if(nrj0==0) nrj0=e_k+e_v;
      else
      {
         // Apply a coefficient to the speed to keep the overall energy constant
         const REAL de=e_k+e_v-nrj0;
         if(de<e_k) v_r=sqrt((e_k-de)/e_k);
         else v_r=0.0;
      }

      #if 0
      char buf[100];
      sprintf(buf,"(i) LLK + Ek = %10.3f + %10.3f =%10.3f (nrj0=%10.3f)",e_v,e_k,e_v+e_k,nrj0);
      cout<<buf<<endl;
      
      // Display gradient & speed - for tests
      for(map<MolAtom*,XYZ>::iterator pos=v0.begin();pos!=v0.end();++pos)
      {
         sprintf(buf,"%10s xyz= (%8.3f %8.3f %8.3f) v= (%8.3f %8.3f %8.3f) m*a= (%8.3f %8.3f %8.3f)",
               pos->first->GetName().c_str(),
               pos->first->GetX(),pos->first->GetY(),pos->first->GetZ(),
                v0[*pos].x    ,v0[*pos].y     ,v0[*pos].z,
               -grad[*pos].x*im,-grad[*pos].y*im,-grad[*pos].z*im);
         //cout<<buf<<endl;
      }
      //cout<<endl<<endl;
      #endif
      // Move according to max step to minimize LLK
      for(map<MolAtom*,XYZ>::const_iterator pos=v0.begin();pos!=v0.end();++pos)
      {
         const XYZ *pa=&(grad[pos->first]);
         pos->first->SetX(pos->first->GetX()+pos->second.x*dt*v_r-0.5*im*dt*dt*pa->x);
         pos->first->SetY(pos->first->GetY()+pos->second.y*dt*v_r-0.5*im*dt*dt*pa->y);
         pos->first->SetZ(pos->first->GetZ()+pos->second.z*dt*v_r-0.5*im*dt*dt*pa->z);
      }
      // Compute new speed
      for(map<MolAtom*,XYZ>::iterator pos=v0.begin();pos!=v0.end();++pos)
      {
         const XYZ *pa=&(grad[pos->first]);
         pos->second.x = v_r*pos->second.x - pa->x*dt*im;
         pos->second.y = v_r*pos->second.y - pa->y*dt*im;
         pos->second.z = v_r*pos->second.z - pa->z*dt*im;
      }
   }
}

const vector<MolAtom*>& Molecule::GetAtomList()const{return mvpAtom;}
const vector<MolBond*>& Molecule::GetBondList()const{return mvpBond;}
const vector<MolBondAngle*>& Molecule::GetBondAngleList()const{return mvpBondAngle;}
const vector<MolDihedralAngle*>& Molecule::GetDihedralAngleList()const{return mvpDihedralAngle;}

vector<MolAtom*>& Molecule::GetAtomList(){return mvpAtom;}
vector<MolBond*>& Molecule::GetBondList(){return mvpBond;}
vector<MolBondAngle*>& Molecule::GetBondAngleList(){return mvpBondAngle;}
vector<MolDihedralAngle*>& Molecule::GetDihedralAngleList(){return mvpDihedralAngle;}

list<StretchModeBondLength>& Molecule::GetStretchModeBondLengthList(){return mvStretchModeBondLength;}
list<StretchModeBondAngle>& Molecule::GetStretchModeBondAngleList(){return mvStretchModeBondAngle;}
list<StretchModeTorsion>& Molecule::GetStretchModeTorsionList(){return mvStretchModeTorsion;}

const list<StretchModeBondLength>& Molecule::GetStretchModeBondLengthList()const{return mvStretchModeBondLength;}
const list<StretchModeBondAngle>& Molecule::GetStretchModeBondAngleList()const{return mvStretchModeBondAngle;}
const list<StretchModeTorsion>& Molecule::GetStretchModeTorsionList()const{return mvStretchModeTorsion;}

const std::vector<RigidGroup*>& Molecule::GetRigidGroupList()const{return mvRigidGroup;}
std::vector<RigidGroup*>& Molecule::GetRigidGroupList(){return mvRigidGroup;}

void Molecule::RotateAtomGroup(const MolAtom &at1,const MolAtom &at2,
                               const set<MolAtom *> &atoms, const REAL angle,
                               const bool keepCenter)
{
   const REAL vx=at2.X()-at1.X();
   const REAL vy=at2.Y()-at1.Y();
   const REAL vz=at2.Z()-at1.Z();
   this->RotateAtomGroup(at1,vx,vy,vz,atoms,angle,keepCenter);
}
void Molecule::RotateAtomGroup(const MolAtom &at,const REAL vx,const REAL vy,const REAL vz,
                               const set<MolAtom *> &atoms, const REAL angle,
                               const bool keepCenter)
{
   TAU_PROFILE("Molecule::RotateAtomGroup(MolAtom&,vx,vy,vz,...)","void (...)",TAU_DEFAULT);
   if(atoms.size()==0) return;
   const REAL x0=at.X();
   const REAL y0=at.Y();
   const REAL z0=at.Z();
   // :KLUDGE: ? Refuse to do anything if vector is not well defined
   if((fabs(vx)+fabs(vy)+fabs(vz))<1e-6) return;
   REAL dx=0.,dy=0.,dz=0.;
   bool keepc=keepCenter;
   if(keepc)
      if(  (this->GetPar(mXYZ.data()  ).IsFixed())
         ||(this->GetPar(mXYZ.data()+1).IsFixed())
         ||(this->GetPar(mXYZ.data()+2).IsFixed())) keepc=false;
   #if 0
   const REAL ca=cos(angle),sa=sin(angle);
   const REAL ca1=1-ca;
   const REAL vnorm=1/sqrt(vx*vx+vy*vy+vz*vz);
   const REAL ux=vx*vnorm,uy=vy*vnorm,uz=vz*vnorm;
   const REAL m00=ca+ux*ux*ca1;// See http://en.wikipedia.org/wiki/Rotation_matrix
   const REAL m01=ux*uy*ca1-uz*sa;
   const REAL m02=ux*uz*ca1+uy*sa;
   const REAL m10=uy*ux*ca1+uz*sa; // :TODO: Check formulas !
   const REAL m11=ca+uy*uy*ca1;
   const REAL m12=uy*uz*ca1-ux*sa;
   const REAL m20=uz*ux*ca1-uy*sa;
   const REAL m21=uz*uy*ca1+ux*sa;
   const REAL m22=ca+uz*uz*ca1;
   for(set<MolAtom *>::const_iterator pos=atoms.begin();pos!=atoms.end();++pos)
   {
      if(keepc)
      {
         dx -= (*pos)->X();
         dy -= (*pos)->Y();
         dz -= (*pos)->Z();
      }
      const REAL x=(*pos)->X() - x0;
      const REAL y=(*pos)->Y() - y0;
      const REAL z=(*pos)->Z() - z0;
      
      (*pos)->X() = m00*x+m01*y+m02*z+x0;
      (*pos)->Y() = m10*x+m11*y+m12*z+y0;
      (*pos)->Z() = m20*x+m21*y+m22*z+z0;
      if(keepc)
      {
         dx += (*pos)->X();
         dy += (*pos)->Y();
         dz += (*pos)->Z();
      }
   }
   #else
   const Quaternion quat=Quaternion::RotationQuaternion(angle,vx,vy,vz);
   for(set<MolAtom *>::const_iterator pos=atoms.begin();pos!=atoms.end();++pos)
   {
      if(keepc)
      {
         dx -= (*pos)->X();
         dy -= (*pos)->Y();
         dz -= (*pos)->Z();
      }
      (*pos)->X() -= x0;
      (*pos)->Y() -= y0;
      (*pos)->Z() -= z0;
      quat.RotateVector((*pos)->X(),(*pos)->Y(),(*pos)->Z());
      (*pos)->X() += x0;
      (*pos)->Y() += y0;
      (*pos)->Z() += z0;
      if(keepc)
      {
         dx += (*pos)->X();
         dy += (*pos)->Y();
         dz += (*pos)->Z();
      }
   }
   #endif
   // (dx,dy,dz) = vector of the translation of the center of the molecule due to the rotation
   if(keepc)
   {
      dx /= (REAL)(this->GetNbComponent());
      dy /= (REAL)(this->GetNbComponent());
      dz /= (REAL)(this->GetNbComponent());
      mQuat.RotateVector(dx,dy,dz);
      this->GetCrystal().OrthonormalToFractionalCoords(dx,dy,dz);
      mXYZ(0) += dx;
      mXYZ(1) += dy;
      mXYZ(2) += dz;
   }
   mClockAtomPosition.Click();
   mClockScatterer.Click();
}
void Molecule::TranslateAtomGroup(const set<MolAtom *> &atoms, 
                                  const REAL dx,const REAL dy,const REAL dz,
                                  const bool keepCenter)
{
   for(set<MolAtom *>::const_iterator pos=atoms.begin();pos!=atoms.end();++pos)
   {
      (*pos)->X() += dx;
      (*pos)->Y() += dy;
      (*pos)->Z() += dz;
   }
   bool keepc=keepCenter;
   if(keepc)
      if(  (this->GetPar(mXYZ.data()  ).IsFixed())
         ||(this->GetPar(mXYZ.data()+1).IsFixed())
         ||(this->GetPar(mXYZ.data()+2).IsFixed())) keepc=false;
   if(keepc)
   {
      const REAL r= (REAL)(atoms.size())/(REAL)(this->GetNbComponent());
      REAL dxc=dx*r,dyc=dy*r,dzc=dz*r;
      this->GetCrystal().OrthonormalToFractionalCoords(dxc,dyc,dzc);
      mXYZ(0) += dxc;
      mXYZ(1) += dyc;
      mXYZ(2) += dzc;
   }
   mClockAtomPosition.Click();
   mClockScatterer.Click();
}
void Molecule::RestraintExport(ostream &os)const
{
   VFN_DEBUG_ENTRY("Molecule::RestraintExport()",5)
   os<<"BondName, IdealLength, Length, log(likelihood)"<<endl;
   for(vector<MolBond*>::const_iterator pos=mvpBond.begin();pos!=mvpBond.end();++pos)
        os <<(*pos)->GetName()
           <<", "<<(*pos)->GetLength0()
           <<", "<<(*pos)->GetLength()
           <<", "<<(*pos)->GetLogLikelihood()<<endl;
   os<<"BondAngle, IdealAngle, Angle, log(likelihood)"<<endl;
   for(vector<MolBondAngle*>::const_iterator pos=mvpBondAngle.begin();
       pos!=mvpBondAngle.end();++pos)
        os <<(*pos)->GetName()
           <<", "<<(*pos)->Angle0()*180/M_PI
           <<", "<<(*pos)->GetAngle()*180/M_PI
           <<", "<<(*pos)->GetLogLikelihood()<<endl;
   os<<"DihedralAngle, IdealAngle, Angle, log(likelihood)"<<endl;
   for(vector<MolDihedralAngle*>::const_iterator pos=mvpDihedralAngle.begin();
       pos!=mvpDihedralAngle.end();++pos)
        os <<(*pos)->GetName()
           <<", "<<(*pos)->Angle0()*180/M_PI
           <<", "<<(*pos)->GetAngle()*180/M_PI
           <<", "<<(*pos)->GetLogLikelihood()<<endl;
   VFN_DEBUG_EXIT("Molecule::RestraintExport()",5)
}
void Molecule::RestraintStatus(ostream &os)const
{
   VFN_DEBUG_ENTRY("Molecule::RestraintStatus()",5)
   for(vector<MolBond*>::const_iterator pos=mvpBond.begin();pos!=mvpBond.end();++pos)
      cout <<"Bond "<<(*pos)->GetName()
           <<", IdealLength="<<(*pos)->GetLength0()
           <<", Length="<<(*pos)->GetLength()
           <<", log(likelihood)="<<(*pos)->GetLogLikelihood()<<endl;
   for(vector<MolBondAngle*>::const_iterator pos=mvpBondAngle.begin();
       pos!=mvpBondAngle.end();++pos)
      cout <<"Bond Angle "<<(*pos)->GetName()
           <<", IdealAngle="<<(*pos)->Angle0()*180/M_PI
           <<", Angle="<<(*pos)->GetAngle()*180/M_PI
           <<", log(likelihood)="<<(*pos)->GetLogLikelihood()<<endl;
   for(vector<MolDihedralAngle*>::const_iterator pos=mvpDihedralAngle.begin();
       pos!=mvpDihedralAngle.end();++pos)
      cout <<"Dihedral Angle "<<(*pos)->GetName()
           <<", IdealAngle="<<(*pos)->Angle0()*180/M_PI
           <<", Angle="<<(*pos)->GetAngle()*180/M_PI
           <<", log(likelihood)="<<(*pos)->GetLogLikelihood()<<endl;
   VFN_DEBUG_EXIT("Molecule::RestraintStatus()",5)
}

const map<MolAtom *,set<MolAtom *> > &Molecule::GetConnectivityTable()
{
   this->BuildConnectivityTable();
   return mConnectivityTable;
}

RefinableObjClock& Molecule::GetBondListClock(){return mClockBondList;}
const RefinableObjClock& Molecule::GetBondListClock()const{return mClockBondList;}

RefinableObjClock& Molecule::GetAtomPositionClock(){return mClockAtomPosition;}
const RefinableObjClock& Molecule::GetAtomPositionClock()const{return mClockAtomPosition;}

      RefinableObjClock& Molecule::GetRigidGroupClock(){return mClockRigidGroup;}
const RefinableObjClock& Molecule::GetRigidGroupClock()const{return mClockRigidGroup;}

void Molecule::RigidifyWithDihedralAngles()
{
   this->BuildConnectivityTable();
   for(vector<MolBond*>::iterator bond=mvpBond.begin();bond!=mvpBond.end();++bond)
   {
      MolAtom* at2=&((*bond)->GetAtom1());
      MolAtom* at3=&((*bond)->GetAtom2());
      for(set<MolAtom *>::const_iterator c2=mConnectivityTable[at2].begin();
          c2!=mConnectivityTable[at2].end();++c2)
      {
         //MolAtom* at1=mvpAtom[*c2];
         if(*c2==at3) continue;
         if(GetBondAngle(**c2,*at2,*at3)<(10 *DEG2RAD)) continue;
         if(GetBondAngle(**c2,*at2,*at3)>(180*DEG2RAD)) continue;
         for(set<MolAtom*>::const_iterator c3=mConnectivityTable[at3].begin();
             c3!=mConnectivityTable[at3].end();++c3)
         {
            //MolAtom* at4=mvpAtom[*c3];
            if((*c3==at2)||(*c3==*c2)) continue;
            if(GetBondAngle(*at2,*at3,**c3)<(10 *DEG2RAD)) continue;
            if(GetBondAngle(*at2,*at3,**c3)>(180*DEG2RAD)) continue;
            if(this->FindDihedralAngle(**c2,*at2,*at3,**c3)==mvpDihedralAngle.end())
            {
               const REAL dihed=GetDihedralAngle(**c2,*at2,*at3,**c3);
               this->AddDihedralAngle(**c2,*at2,*at3,**c3,dihed,.01,.05,false);
            }
         }
      }
   }
   this->UpdateDisplay();
}
#if 0
/** x array for the Integral[0->x](exp(+t*t)dt)
*
*/
static const REAL svGaussianIntX[51]
   ={0. ,  0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8,  0.9,  1. ,  1.1,  1.2,
     1.3,  1.4,  1.5,  1.6,  1.7,  1.8,  1.9,  2. ,  2.1,  2.2,  2.3,
     2.4,  2.5,  2.6,  2.7,  2.8,  2.9,  3. ,  3.1,  3.2,  3.3,  3.4,
     3.5,  3.6,  3.7,  3.8,  3.9,  4. ,  4.1,  4.2,  4.3,  4.4,  4.5,
     4.6,  4.7,  4.8,  4.9,  5. };
/** This is Integral[0->x](exp(+t*t)dt)
*
*/
static const REAL svGaussianIntY[51]
   ={0.00000000e+00,   1.00285768e-01,   2.02498389e-01,   3.08782899e-01,
     4.21537858e-01,   5.43577677e-01,   6.78339751e-01,   8.30161516e-01,
     1.00466359e+00,   1.20929269e+00,   1.45410564e+00,   1.75292008e+00,
     2.12502806e+00,   2.59778353e+00,   3.21056320e+00,   4.02091257e+00,
     5.11421527e+00,   6.61911896e+00,   8.73249677e+00,   1.17604260e+01,
     1.61864569e+01,   2.27870547e+01,   3.28297946e+01,   4.84189023e+01,
     7.31071567e+01,   1.12996748e+02,   1.78751757e+02,   2.89337246e+02,
     4.79080996e+02,   8.11232960e+02,   1.40443983e+03,   2.48531490e+03,
     4.49461494e+03,   8.30539633e+03,   1.56790580e+04,   3.02354014e+04,
     5.95525189e+04,   1.19793234e+05,   2.46080284e+05,   5.16182008e+05,
     1.10556243e+06,   2.41765307e+06,   5.39775929e+06,   1.23033271e+07,
     2.86288375e+07,   6.80050408e+07,   1.64899866e+08,   4.08157783e+08,
     1.03122228e+09,   2.65938808e+09,   7.00012860e+09};
/** Spline giving x=f(y) for the Integral[0->x](exp(+t*t)dt)
*
*/
static const CubicSpline sInvGaussianInt(svGaussianIntY,svGaussianIntX,51);

REAL FlatGaussianProba(const REAL x, const REAL sigma, const REAL delta)
{
   if(abs(x)<delta) return 1.0;
   const REAL z=(abs(x)-delta)/sigma;
   return exp(-z*z);
}

REAL FlatGaussianIntegral(const REAL x1,const REAL x2, const REAL sigma, const REAL delta)
{
   static const REAL SPI2=0.88622692545275794;//sqrt(pi)/2
   if((x1<=-delta)&&(x2<=-delta)) return SPI2*sigma*(erf((x2+delta)/sigma)-erf((x1+delta)/sigma));
   if((x1<=-delta)&&(x2<= delta)) return (x2+delta)-SPI2*sigma*erf((x1+delta)/sigma);
   if(x1<=-delta) return 2*delta+SPI2*sigma*(erf((x2-delta)/sigma)-erf((x1+delta)/sigma));
   if((x1<= delta)&&(x2<= delta)) return x2-x1;
   if(x1<= delta) return SPI2*sigma*erf((x2-delta)/sigma)+(delta-x1);
   return SPI2*sigma*(erf((x2-delta)/sigma)-erf((x1-delta)/sigma));
}
#endif
REAL FlatLorentzianProba(const REAL x, const REAL sigma, const REAL delta)
{
   if(abs(x)<delta) return 1.0;
   const REAL z=(abs(x)-delta)/sigma;
   return 1/(1+z*z);
}

REAL FlatLorentzianIntegral(const REAL x1,const REAL x2, const REAL sigma, const REAL delta)
{
   if((x1<=-delta)&&(x2<=-delta)) return atan((x2+delta)/sigma)-atan((x1+delta)/sigma);
   if((x1<=-delta)&&(x2<= delta)) return (x2+delta)-atan((x1+delta)/sigma);
   if(x1<=-delta) return 2*delta+atan((x2-delta)/sigma)-atan((x1+delta)/sigma);
   if((x1<= delta)&&(x2<= delta)) return x2-x1;
   if(x1<= delta) return atan((x2-delta)/sigma)+(delta-x1);
   return atan((x2-delta)/sigma)-atan((x1-delta)/sigma);
}

ofstream f;

/** Random move respecting a gaussian probability distribution with a flat top.
* i.e. x in [-delta;+delta], P(x)=1
* outside, P(x)=1/(1+(abs(x)-delta)^2/sigma^2)
*
* If sigma<1e-6, it is treated as a step-like probability non-null only in [-delta;+delta]
*/
REAL LorentzianBiasedRandomMove(const REAL x0,const REAL sigma,const REAL delta,const REAL amplitude)
{
   //static const REAL SPI2=0.88622692545275794;//sqrt(pi)/2
   REAL r=(REAL)rand()/(REAL)RAND_MAX;
   if(sigma<1e-6)
   {
      REAL x=x0+amplitude*(2*r-1.0);
      if(x> delta)x= delta;
      if(x<-delta)x=-delta;
      return x;
   }
   // Compute xmin and xmax around x0 so that:
   // (x0-xmin)/(xmax-x0)=Integral[xmin->x0](P(x)dx)/Integral[x0->xmax](P(x)dx)
   REAL xmin;
   REAL xmax;
   #if 1
   {
      REAL Pmax,Pmin;
      xmin=x0-amplitude;
      xmax=x0+amplitude;
      Pmin=FlatLorentzianProba((x0+xmin)/2,sigma,delta);
      Pmax=FlatLorentzianProba((x0+xmax)/2,sigma,delta);
      //Pmin=FlatLorentzianIntegral(xmin,x0,sigma,delta);
      //Pmax=FlatLorentzianIntegral(x0,xmax,sigma,delta);
      if(Pmin>Pmax)
      {
         for(int i=0;i<5;i++)
         {
            const REAL r=Pmax/Pmin;
            xmax=x0+r*amplitude;
            Pmax=FlatLorentzianProba((x0+xmax)/2,sigma,delta);
            //Pmax=FlatLorentzianIntegral(x0,xmax,sigma,delta);
         }
      }
      else
      {
         for(int i=0;i<5;i++)
         {
            const REAL r=Pmin/Pmax;
            xmin=x0-r*amplitude;
            Pmin=FlatLorentzianProba((x0+xmin)/2,sigma,delta);
            //Pmin=FlatLorentzianIntegral(xmin,x0,sigma,delta);
         }
      }
   }
   #else
   if(abs(x0)<=delta)
   {
      xmin=x0-amplitude;
      xmax=x0+amplitude;
   }
   else
   {
      REAL d=(abs(x0)-delta)/sigma;
      d=2*d*exp(-d*d);
      d=(1-amplitude*d/2)/(1+amplitude*d/2);//d<1
      if(x0>0)
      {
         xmin=x0-amplitude*2/(1+d);
         xmax=x0+amplitude*2*d/(1+d);
      }
      else
      {
         xmax=x0+amplitude*2/(1+d);
         xmin=x0-amplitude*2*d/(1+d);
      }
   }
   #endif
   //xmin=x0-amplitude;
   //xmax=x0+amplitude;
   //Now get the biased move...
   if(xmax<=-delta)
   {
      REAL ymin=(abs(xmin)-delta)/sigma;
      ymin=atan(ymin);
      REAL ymax=(abs(xmax)-delta)/sigma;
      ymax=atan(ymax);
      const REAL y=ymin+(ymax-ymin)*r;
      return -tan(y)*sigma-delta;
   }
   if(xmin<=-delta)
   {
      if(xmax<=delta) 
      {
         //probability of being in [xmin;-delta] rather than [-delta;xmax]
         REAL p0=atan((abs(xmin)-delta)/sigma)*sigma;
         REAL p1=xmax+delta;
         const REAL n=p0+p1;
         if(r<p0/n)
         {
            REAL ymin=(abs(xmin)-delta)/sigma;
            ymin=atan(ymin);
            const REAL y=ymin*(REAL)rand()/(REAL)RAND_MAX;
            return -delta-tan(y)*sigma;
         }
         else
         {
            return -delta+(REAL)rand()/(REAL)RAND_MAX*(xmax+delta);
         }
      }
      else //xmax>delta && xmin <= -delta
      {
         REAL p0=atan((abs(xmin)-delta)/sigma)*sigma;//probability of being in [xmin;-delta]
         REAL p1=2*delta;//probability of being in [-delta;delta]
         REAL p2=atan((xmax-delta)/sigma)*sigma;//probability of being in [delta;xmax]
         const REAL n=p0+p1+p2;
         if(r<(p0/n))
         {
            REAL ymin=(abs(xmin)-delta)/sigma;
            ymin=atan(ymin);//exp(ymin*ymin);
            const REAL y=ymin*(REAL)rand()/(REAL)RAND_MAX;
            const REAL x=-delta-tan(y)*sigma;
            return x;
         }
         if(r<(p0+p1)/n)
         {
            const REAL x=-delta+(REAL)rand()/(REAL)RAND_MAX*2*delta;
            return x;
         }

         REAL ymax=(xmax-delta)/sigma;
         ymax=atan(ymax);
         const REAL y=ymax*(REAL)rand()/(REAL)RAND_MAX;
         const REAL x=delta+tan(y)*sigma;
         return x;
      }
   }
   if(xmin<delta)
   {
      if(xmax<delta)
      {
         return xmin+r*(xmax-xmin);
      }

      const REAL p0=delta-xmin;//relative probability of being in [xmin;delta]
      const REAL p1=atan((xmax-delta)/sigma)*sigma;// proba in[delta;xmax]
      if(r<(p0/(p0+p1)))
      {
         return xmin+(REAL)rand()/(REAL)RAND_MAX*(delta-xmin);
      }

      REAL ymax=(xmax-delta)/sigma;
      ymax=atan(ymax);
      const REAL y=ymax*(REAL)rand()/(REAL)RAND_MAX;
      return delta+tan(y)*sigma;
   }
   //xmin>delta
   REAL ymin=(xmin-delta)/sigma;
   ymin=atan(ymin);
   REAL ymax=(xmax-delta)/sigma;
   ymax=atan(ymax);
   const REAL y=ymin+(ymax-ymin)*r;
   return tan(y)*sigma+delta;
}

void TestLorentzianBiasedRandomMove()
{
   srand(time(NULL));
   REAL x=0,sigma=0.1,delta=0.5,amplitude=0.05;
   f.open("test.dat");
   for(long i=0;i<400000;i++)
   {
      f<<x<<endl;
      x=LorentzianBiasedRandomMove(x,sigma,delta,amplitude);
   }
   f.close();
   exit(0);
   //#Histogram in Python
   //from scipy import *
   //
   //def histogram(y,n=100):
   //   ymin=min(y)
   //   ymax=max(y)
   //   step=(ymax-ymin)/n
   //   hx=arange(ymin,ymax+step,step)
   //   hy=arange(ymin,ymax+step,step)
   //   for i in xrange(n-1):
   //      hy[i]=sum((y>hx[i])*(y<(hx[i]+step)))
   //   return hx,hy
   //
   //
   //f=open("test.dat",'r')
   //ll=f.readlines()
   //f.close()
   //y=zeros(len(ll),Float)
   //for i in xrange(len(ll)):
   //   y[i]=float(ll[i])
   //
   //hx,hy=histogram(y)
   //gplt.plot(hx,hy)
}  

REAL Molecule::BondLengthRandomChange(const StretchModeBondLength& mode, const REAL amplitude,
                                      const bool respectRestraint)
{
   REAL dx=mode.mpAtom1->GetX()-mode.mpAtom0->GetX();
   REAL dy=mode.mpAtom1->GetY()-mode.mpAtom0->GetY();
   REAL dz=mode.mpAtom1->GetZ()-mode.mpAtom0->GetZ();
   const REAL l=sqrt(dx*dx+dy*dy+dz*dz+1e-7);
   REAL change=0.0;
   if(l<1e-6) return change;// :KLUDGE:
   if(respectRestraint && mode.mpBond!=0)
   {
      const REAL d0=l-mode.mpBond->GetLength0();
      const REAL sigma=mode.mpBond->GetLengthSigma();
      const REAL delta=mode.mpBond->GetLengthDelta();
      const REAL max=delta+sigma*5.0;
      if(sigma<1e-6)
      {
         REAL d1=d0+(REAL)(2*rand()-RAND_MAX)/(REAL)RAND_MAX*amplitude*0.1;
         if(d1> delta)d1= delta;
         if(d1<-delta)d1=-delta;
         change=d1-d0;
      }
      else change=LorentzianBiasedRandomMove(d0,sigma,delta,amplitude*0.1)-d0;
      if((d0+change)>max) change=max-d0;
      else if((d0+change)<(-max)) change=-max-d0;
      #if 0
      if(rand()%10000==0)
      {
         cout<<"BOND LENGTH change("<<change<<"):"
             <<mode.mpAtom0->GetName()<<"-"
             <<mode.mpAtom1->GetName()
             <<"(Restraint="<<l-d0<<"s"<<sigma<<"d"<<delta<<"):"
             <<l<<"->"<<l+change<<endl;
      }
      #endif
   }
   else change=(2.*(REAL)rand()-(REAL)RAND_MAX)/(REAL)RAND_MAX*amplitude*0.1;
   dx*=change/l;
   dy*=change/l;
   dz*=change/l;
   this->TranslateAtomGroup(mode.mvTranslatedAtomList,dx,dy,dz,true);
   return change;
}

REAL Molecule::BondAngleRandomChange(const StretchModeBondAngle& mode, const REAL amplitude,
                                     const bool respectRestraint)
{
   REAL dx10=mode.mpAtom0->GetX()-mode.mpAtom1->GetX();
   REAL dy10=mode.mpAtom0->GetY()-mode.mpAtom1->GetY();
   REAL dz10=mode.mpAtom0->GetZ()-mode.mpAtom1->GetZ();
   REAL dx12=mode.mpAtom2->GetX()-mode.mpAtom1->GetX();
   REAL dy12=mode.mpAtom2->GetY()-mode.mpAtom1->GetY();
   REAL dz12=mode.mpAtom2->GetZ()-mode.mpAtom1->GetZ();
   
   const REAL vx=dy10*dz12-dz10*dy12;
   const REAL vy=dz10*dx12-dx10*dz12;
   const REAL vz=dx10*dy12-dy10*dx12;
   
   REAL change=0.0;
   if((abs(vx)+abs(vy)+abs(vz))<1e-6) return change;// :KLUDGE:
   REAL angle0;
   if(respectRestraint && mode.mpBondAngle!=0)
   {
      const REAL norm= sqrt( (dx10*dx10+dy10*dy10+dz10*dz10)*(dx12*dx12+dy12*dy12+dz12*dz12)+1e-6);
      angle0=(dx10*dx12+dy10*dy12+dz10*dz12)/norm;
      if(angle0>=1)  angle0=0;
      else 
      {
         if(angle0<=-1) angle0=M_PI;
         else angle0= acos(angle0);
      }
      
      const REAL a0=angle0-mode.mpBondAngle->GetAngle0();
      const REAL sigma=mode.mpBondAngle->GetAngleSigma();
      const REAL delta=mode.mpBondAngle->GetAngleDelta();
      if(sigma<1e-6)
      {
         REAL a1=a0+(REAL)(2*rand()-RAND_MAX)/(REAL)RAND_MAX*amplitude*mode.mBaseAmplitude;
         if(a1> delta)a1= delta;
         if(a1<-delta)a1=-delta;
         change=a1-a0;
      }
      else change=LorentzianBiasedRandomMove(a0,sigma,delta,amplitude*mode.mBaseAmplitude)-a0;
      if((a0+change)>(delta+sigma*5.0))       change= delta+sigma*5.0-a0;
      else if((a0+change)<(-delta-sigma*5.0)) change=-delta-sigma*5.0-a0;
      #if 0
      if(rand()%1==0)
      {
         cout<<"ANGLE change("<<change*RAD2DEG<<"):"
             <<mode.mpAtom0->GetName()<<"-"
             <<mode.mpAtom1->GetName()<<"-"
             <<mode.mpAtom2->GetName()
             <<"(Restraint="<<(angle0-a0)*RAD2DEG<<"s"<<sigma*RAD2DEG<<"d"<<delta*RAD2DEG<<"):"
             <<angle0*RAD2DEG<<"->"<<(angle0+change)*RAD2DEG<<endl;
      }
      #endif
   }
   else change=(2.*(REAL)rand()-(REAL)RAND_MAX)/(REAL)RAND_MAX*mode.mBaseAmplitude*amplitude;
   this->RotateAtomGroup(*(mode.mpAtom1),vx,vy,vz,mode.mvRotatedAtomList,change,true);
   return change;
}
REAL Molecule::DihedralAngleRandomChange(const StretchModeTorsion& mode, const REAL amplitude,
                                         const bool respectRestraint)
{
   const REAL dx=mode.mpAtom2->GetX()-mode.mpAtom1->GetX();
   const REAL dy=mode.mpAtom2->GetY()-mode.mpAtom1->GetY();
   const REAL dz=mode.mpAtom2->GetZ()-mode.mpAtom1->GetZ();
   REAL change=0.0;
   if((abs(dx)+abs(dy)+abs(dz))<1e-6) return change;// :KLUDGE:
   if(respectRestraint && mode.mpDihedralAngle!=0)
   {
      const REAL angle0=mode.mpDihedralAngle->GetAngle();
      const REAL a0=angle0-mode.mpDihedralAngle->GetAngle0();
      const REAL sigma=mode.mpDihedralAngle->GetAngleSigma();
      const REAL delta=mode.mpDihedralAngle->GetAngleDelta();
      if(sigma<1e-6)
      {
         REAL a1=a0+(REAL)(2*rand()-RAND_MAX)/(REAL)RAND_MAX*amplitude*mode.mBaseAmplitude;
         if(a1> delta)a1= delta;
         if(a1<-delta)a1=-delta;
         change=a1-a0;
      }
      else change=LorentzianBiasedRandomMove(a0,sigma,delta,amplitude*mode.mBaseAmplitude)-a0;
      if((a0+change)>(delta+sigma*5.0))       change= delta+sigma*5.0-a0;
      else if((a0+change)<(-delta-sigma*5.0)) change=-delta-sigma*5.0-a0;
      #if 0
      if(rand()%1==0)
      {
         cout<<"TORSION change ("
             <<mode.mpAtom1->GetName()<<"-"<<mode.mpAtom2->GetName()<<"):"<<endl
             <<"     initial angle="<<angle0*RAD2DEG
             <<" (Restraint("<<mode.mpDihedralAngle->GetName()<<")="
             <<(angle0-a0)*RAD2DEG<<"s"<<sigma*RAD2DEG<<"d"<<delta*RAD2DEG<<"):"
             <<endl<<"     New angle:"<<(angle0+change)*RAD2DEG<<", change="<<change*RAD2DEG<<endl;
      }
      #endif
   }
   else change=(REAL)(2.*rand()-RAND_MAX)/(REAL)RAND_MAX*mode.mBaseAmplitude*amplitude;
   this->RotateAtomGroup(*(mode.mpAtom1),*(mode.mpAtom2),mode.mvRotatedAtomList,change,true);
   return change;
}

const MolAtom* Molecule::GetCenterAtom()const
{
   return mpCenterAtom;
}

void Molecule::SetCenterAtom(const MolAtom &at)
{
   mpCenterAtom=&at;
   mClockAtomPosition.Click();
   this->UpdateDisplay();
}

void BuildZMatrixRecursive(long &z,const long curr,
                           const vector<MolAtom*> &vpAtom,
                           const map<MolAtom *, set<MolAtom *> > &connT,
                           vector<MolZAtom> &zmatrix,
                           const map<const MolAtom*,long> &vIndex,
                           vector<long> &vZIndex,
                           vector<long> &vrZIndex)
{
   zmatrix[z].mpPow=&(vpAtom[curr]->GetScatteringPower());
   vZIndex[curr]=z;
   vrZIndex[z]=curr;
   const long n=vpAtom.size();
   // Get the list of connected atoms and sort them
      map<MolAtom *, set<MolAtom *> >::const_iterator pConn=connT.find(vpAtom[curr]);
      const long nc=pConn->second.size();
      vector<long> conn(nc);
      vector<long> zconn(nc);
      vector<long>::iterator pos=conn.begin();
      vector<long>::iterator zpos=zconn.begin();
      for(set<MolAtom *>::const_iterator pos1=pConn->second.begin();pos1!=pConn->second.end();++pos1)
      {
         *pos = vIndex.find(*pos1)->second;
         *zpos = vZIndex[*pos];
         cout<<(*pos1)->GetName()<<"("<<*pos<<","<<*zpos<<")"<<endl;
         zpos++;pos++;
      }
      sort(conn.begin(),conn.end());
      sort(zconn.begin(),zconn.end());
   if(z>0)
   {
      // Use the most recent atom in the z-matrix
      const long b=zconn[nc-1];
      zmatrix[z].mBondAtom=b;
      zmatrix[z].mBondLength=GetBondLength(*vpAtom[vrZIndex[b]],*vpAtom[curr]);
      if(z>1)
      {
         const long a=zmatrix[b].mBondAtom;
         zmatrix[z].mBondAngleAtom=a;
         zmatrix[z].mBondAngle=GetBondAngle(*vpAtom[vrZIndex[a]],*vpAtom[vrZIndex[b]],*vpAtom[curr]);
         if(z>2)
         {
            const long d=zmatrix[b].mBondAngleAtom;
            zmatrix[z].mDihedralAtom=d;
            zmatrix[z].mDihedralAngle=fmod(GetDihedralAngle(*vpAtom[vrZIndex[d]],*vpAtom[vrZIndex[a]],
                                                            *vpAtom[vrZIndex[b]],*vpAtom[curr])+2*M_PI,
                                           2*M_PI);
         }
         else
         {
            zmatrix[z].mDihedralAtom=0;
            zmatrix[z].mDihedralAngle=0;
         }
      }
      else
      {
         zmatrix[z].mBondAngleAtom=0;
         zmatrix[z].mBondAngle=0;
      }
   }
   else
   {
      zmatrix[z].mBondAtom=0;
      zmatrix[z].mBondLength=0;
   }
   z++;
   // Continue filling up the zmatrix, beginning from thz first atoms not already in the zmatrix
   for(pos=conn.begin();pos!=conn.end();++pos)
   {
      if(*pos!=-1)
      {
         if(vZIndex[*pos]==-1)
            BuildZMatrixRecursive(z,*pos,vpAtom,connT,zmatrix,vIndex,vZIndex,vrZIndex);
      }
   }
}

const vector<MolZAtom>& Molecule::AsZMatrix(const bool keeporder)const
{
   this->BuildConnectivityTable();
   const long n=mvpAtom.size();
   // index of the atoms in the list
   map<const MolAtom*,long> vIndex;
   {
      long i=0;
      for(vector<MolAtom*>::const_iterator pos=mvpAtom.begin();pos!=mvpAtom.end();++pos)
         vIndex[*pos]=i++;
   }
   mAsZMatrix.resize(n);
   if(keeporder)
   {
      for(long i=0;i<n;++i)
      {
         mAsZMatrix[i].mpPow=&(mvpAtom[i]->GetScatteringPower());
         if(i>0)
         {
            const set<MolAtom *> *pConn=&(mConnectivityTable.find(mvpAtom[i])->second);
            // Find a connected atom already in the mAsZMatrix, prefereably the most recent
            long b=-1;
            for(set<MolAtom *>::const_iterator pos=pConn->begin();pos!=pConn->end();++pos)
               if((vIndex[*pos]<i)&&(vIndex[*pos]>b)) b=vIndex[*pos];
            // Did not find a connected atom already in the z-matrix ? Take the last one
            if(b==-1) b=i-1;
            mAsZMatrix[i].mBondAtom=b;
            mAsZMatrix[i].mBondLength=GetBondLength(*mvpAtom[b],*mvpAtom[i]);
            if(i>1)
            {
               const long a= (b==0)?1 : mAsZMatrix[b].mBondAtom;
               mAsZMatrix[i].mBondAngleAtom=a;
               mAsZMatrix[i].mBondAngle=GetBondAngle(*mvpAtom[a],*mvpAtom[b],*mvpAtom[i]);
               if(i>2)
               {
                  long d= mAsZMatrix[a].mBondAtom;
                  if(d==b)
                  {// Dihedral atom is already bond atom, find another connected to angle atom
                     d=-1;
                     const set<MolAtom *> *pConn=&(mConnectivityTable.find(mvpAtom[a])->second);
                     // Find a connected atom already in the mAsZMatrix, prefereably the most recent
                     for(set<MolAtom *>::const_iterator pos=pConn->begin();pos!=pConn->end();++pos)
                        if((vIndex[*pos]<i) && (vIndex[*pos]!=b) && (vIndex[*pos]>d)) d=vIndex[*pos];
                  }
                  if(d==-1)
                  {// Can't find an angle connected to angle atom, so find another with bond atom
                     const set<MolAtom *> *pConn=&(mConnectivityTable.find(mvpAtom[b])->second);
                     // Find a connected atom already in the mAsZMatrix, prefereably the most recent
                     for(set<MolAtom *>::const_iterator pos=pConn->begin();pos!=pConn->end();++pos)
                        if((vIndex[*pos]<i) && (vIndex[*pos]!=a) && (vIndex[*pos]>d)) d=vIndex[*pos];
                  }
                  if(d==-1)
                  {// Maybe another connected to this atom ??
                     const set<MolAtom *> *pConn=&(mConnectivityTable.find(mvpAtom[i])->second);
                     // Find a connected atom already in the mAsZMatrix, prefereably the most recent
                     for(set<MolAtom *>::const_iterator pos=pConn->begin();pos!=pConn->end();++pos)
                        if(  (vIndex[*pos]<i)  && (vIndex[*pos]!=a)
                           &&(vIndex[*pos]!=b) && (vIndex[*pos]>d)) d=vIndex[*pos];
                  }
                  if(d==-1)
                  {// OK, pick *any* (can this happen ? Really ?)
                     for(long j=0;j<i;++j)
                        if((j!=a) &&(j!=b) && (j>d)) d=j;
                  }
                  mAsZMatrix[i].mDihedralAtom=d;
                  mAsZMatrix[i].mDihedralAngle=fmod(GetDihedralAngle(*mvpAtom[d],
                                                                     *mvpAtom[a],
                                                                     *mvpAtom[b],
                                                                     *mvpAtom[i])+2*M_PI,
                                                 2*M_PI);
               }
            }
         }
      }
   }
   else
   {
      // vZIndex[i] tells where mvpAtom[i] is in the z-matrix
      vector<long> vZIndex(n);
      // vrZIndex[i] tells where which index in vpAtom is ZAtom #i
      vector<long> vrZIndex(n);
      for(long i=0;i<n;++i)
      {
         vZIndex [i]=-1;
         vrZIndex[i]=-1;
      }
      long z=0;
      BuildZMatrixRecursive(z,0,mvpAtom,mConnectivityTable,mAsZMatrix,vIndex,vZIndex,vrZIndex);
   }
   return mAsZMatrix;
}

void Molecule::InitRefParList()
{
}

/** Find rings, starting from a one atom, and given a connectivity table. The
* ring begins and ends with the given atom.
*
* \param atom: the current atom
* \param connect: the connectivity table
* \returns: the list of atoms to which will be appended the atoms of the ring,
*\b if one is found. Otherwise, an empty set of atoms.
* \param atomlist: the current list of atoms in the list.
*/
void BuildRingRecursive(MolAtom * currentAtom,
                        MolAtom * previousAtom,
                        const map<MolAtom *, set<MolAtom *> > &connect,
                        list<MolAtom *> &atomlist,
                        map<set<MolAtom *>,list<MolAtom *> > &ringlist)
{
   list<MolAtom *>::const_iterator f=find(atomlist.begin(),atomlist.end(),currentAtom);
   if(f!=atomlist.end())
   {// This atom was already in the list ! We have found a ring !
      #ifdef __DEBUG__
      cout<<currentAtom->GetName()<<" was already in the list : ring found !"<<endl;
      for(list<MolAtom *>::const_iterator atom=atomlist.begin();atom!=atomlist.end();++atom)
         cout<<(*atom)->GetName()<<" ";
      cout<<endl;
      #endif
      set<MolAtom *> ring1;
      list<MolAtom *> ring2;
      for(list<MolAtom *>::const_iterator pos=f;pos!=atomlist.end();++pos)
      {
         ring1.insert(*pos);
         ring2.push_back(*pos);
      }
      ringlist.insert(make_pair(ring1,ring2));
   }
   else
   {
      atomlist.push_back(currentAtom);
      map<MolAtom *,set<MolAtom *> >::const_iterator c=connect.find(currentAtom);
      set<MolAtom *>::const_iterator pos;
      for(pos=c->second.begin();pos!=c->second.end();++pos)
      {
         if(*pos==previousAtom) continue;
         BuildRingRecursive(*pos,currentAtom,connect,atomlist,ringlist);
      }
      atomlist.pop_back(); //??
   }
}

void Molecule::BuildRingList()
{
   this->BuildConnectivityTable();
   if(mClockRingList>mClockConnectivityTable) return;
   VFN_DEBUG_ENTRY("Molecule::BuildRingList()",7)
   for(vector<MolAtom*>::const_iterator pos=mvpAtom.begin();pos!=mvpAtom.end();++pos)
      (*pos)->SetIsInRing(false);
   list<MolAtom *> atomlist;
   // Use a map with a set for key to eliminate duplicate rings
   map<set<MolAtom *>,list<MolAtom *> > ringlist;
   for(unsigned long i=0;i<mvpAtom.size();i++)
   {
      atomlist.clear();
      BuildRingRecursive(mvpAtom[i],mvpAtom[i],mConnectivityTable,atomlist,ringlist);
   }
   for(map<set<MolAtom *>,list<MolAtom *> >::const_iterator pos0=ringlist.begin();pos0!=ringlist.end();pos0++)
   {
      mvRing.resize(mvRing.size()+1);
      std::list<MolAtom*> *pList=&(mvRing.back().GetAtomList());
      #if 1//def __DEBUG__
      cout<<"Found ring:";
      #endif
      for(list<MolAtom *>::const_iterator atom=pos0->second.begin();atom!=pos0->second.end();++atom)
      {
         pList->push_back(*atom);
         (*atom)->SetIsInRing(true);
         #if 1//def __DEBUG__
         cout<<(*atom)->GetName()<<" ";
         #endif
      }
      #if 1//def __DEBUG__
      cout<<endl;
      #endif
   }
   
   cout<<"Rings found :"<<ringlist.size()<<", "<<mvRing.size()<<" unique."<<endl;
   mClockRingList.Click();
   VFN_DEBUG_EXIT("Molecule::BuildRingList()",7)
}

void Molecule::BuildConnectivityTable()const
{
   if(  (mClockConnectivityTable>mClockBondList)
      &&(mClockConnectivityTable>mClockAtomList)) return;
   VFN_DEBUG_ENTRY("Molecule::BuildConnectivityTable()",5)
   TAU_PROFILE("Molecule::BuildConnectivityTable()","void ()",TAU_DEFAULT);
   mConnectivityTable.clear();
   for(unsigned long i=0;i<mvpBond.size();++i)
   {
      mConnectivityTable[&(mvpBond[i]->GetAtom1())].insert(&(mvpBond[i]->GetAtom2()));
      mConnectivityTable[&(mvpBond[i]->GetAtom2())].insert(&(mvpBond[i]->GetAtom1()));
   }
   
   #ifdef __DEBUG__
   {
      map<MolAtom *,set<MolAtom *> >::const_iterator pos;
      unsigned long at=0;
      for(pos=mConnectivityTable.begin();pos!=mConnectivityTable.end();++pos)
      {
         cout<<"Atom "<<pos->first->GetName()<<" is connected to atoms: ";
         set<MolAtom *>::const_iterator pos1;
         for(pos1=pos->second.begin();pos1!=pos->second.end();++pos1)
         {
            cout<<(*pos1)->GetName()<<"  ";
         }
         cout<<endl;
         if(pos->second.size()>10) exit(0);
      }
   }
   #endif
   mClockConnectivityTable.Click();
   VFN_DEBUG_EXIT("Molecule::BuildConnectivityTable()",5)
}

Molecule::RotorGroup::RotorGroup(const MolAtom &at1,const MolAtom &at2):
mpAtom1(&at1),mpAtom2(&at2),mBaseRotationAmplitude(M_PI*0.04)
{}

void Molecule::BuildRotorGroup()
{
   if(  (mClockRotorGroup>mClockBondList)
      &&(mClockRotorGroup>mClockAtomList)
      &&(mClockRotorGroup>mClockBondAngleList)
      &&(mClockRotorGroup>mClockDihedralAngleList)) return;
   VFN_DEBUG_ENTRY("Molecule::BuildRotorGroup()",5)
   TAU_PROFILE("Molecule::BuildRotorGroup()","void ()",TAU_DEFAULT);
   this->BuildConnectivityTable();
   mvRotorGroupTorsion.clear();
   mvRotorGroupTorsionSingleChain.clear();
   mvRotorGroupInternal.clear();
   
   // Build Rotation groups around bonds
   for(unsigned long i=0;i<mvpBond.size();++i)
   {
      if((mFlexModel.GetChoice()!=0)&&(false==mvpBond[i]->IsFreeTorsion())) continue;
      MolAtom *const atom1=&(mvpBond[i]->GetAtom1()),
              *const atom2=&(mvpBond[i]->GetAtom2());
      for(unsigned int j=1;j<=2;++j)
      {
         const set<MolAtom*> *pConn;
         if(j==1) pConn=&(mConnectivityTable[atom1]);
         else pConn=&(mConnectivityTable[atom2]);
         
         mvRotorGroupTorsion.push_back(RotorGroup(mvpBond[i]->GetAtom1(),
                                                  mvpBond[i]->GetAtom2()));
         mvRotorGroupTorsion.back().mvRotatedAtomList.insert(atom1);
         mvRotorGroupTorsion.back().mvRotatedAtomList.insert(atom2);
         
         for(set<MolAtom*>::const_iterator pos=pConn->begin();pos!=pConn->end();++pos)
         {
            if((j==1)&&(*pos==atom2)) continue;
            if((j==2)&&(*pos==atom1)) continue;
            ExpandAtomGroupRecursive(*pos,mConnectivityTable,
                                     mvRotorGroupTorsion.back().mvRotatedAtomList);
            if(pConn->size()>2)
            {
               mvRotorGroupTorsionSingleChain.push_back(RotorGroup(mvpBond[i]->GetAtom1(),
                                                                   mvpBond[i]->GetAtom2()));
               mvRotorGroupTorsionSingleChain.back().mvRotatedAtomList.insert(atom1);
               mvRotorGroupTorsionSingleChain.back().mvRotatedAtomList.insert(atom2);
               ExpandAtomGroupRecursive(*pos,mConnectivityTable,
                                        mvRotorGroupTorsionSingleChain.back().mvRotatedAtomList);
               mvRotorGroupTorsionSingleChain.back().mvRotatedAtomList.erase(atom1);
               mvRotorGroupTorsionSingleChain.back().mvRotatedAtomList.erase(atom2);

               if(  (mvRotorGroupTorsionSingleChain.back().mvRotatedAtomList.size()>=((mvpAtom.size()+1)/2))
                  ||(mvRotorGroupTorsionSingleChain.back().mvRotatedAtomList.size()==0))
                  mvRotorGroupTorsionSingleChain.pop_back();
            }
         }
         mvRotorGroupTorsion.back().mvRotatedAtomList.erase(atom1);
         mvRotorGroupTorsion.back().mvRotatedAtomList.erase(atom2);
         if(  (mvRotorGroupTorsion.back().mvRotatedAtomList.size()>=((mvpAtom.size()+1)/2))
            ||(mvRotorGroupTorsion.back().mvRotatedAtomList.size()==0))
           mvRotorGroupTorsion.pop_back();
      }
   }
   #if 1
   // Build 'internal' rotation groups between random atoms
   //:TODO: This should be tried for *random* configuration of free torsion angles...
   if(mFlexModel.GetChoice()==0)
   {
      for(vector<MolAtom*>::const_iterator atom1=this->GetAtomList().begin();
          atom1!=this->GetAtomList().end();++atom1)
      {
         const set<MolAtom*> *pConn=&(mConnectivityTable[*atom1]);
         vector<MolAtom*>::const_iterator atom2=atom1;
         atom2++;
         for(;atom2!=this->GetAtomList().end();++atom2)
         {
            for(set<MolAtom*>::const_iterator pos=pConn->begin();pos!=pConn->end();++pos)
            {
               if(*pos==*atom2) continue;
               mvRotorGroupInternal.push_back(RotorGroup(**atom1,**atom2));
               mvRotorGroupInternal.back().mvRotatedAtomList.insert(*atom1);
               ExpandAtomGroupRecursive(*pos,mConnectivityTable,
                                        mvRotorGroupInternal.back().mvRotatedAtomList,
                                        *atom2);
               //Check if this chains leads to atom2
               set<MolAtom*>::const_iterator check
                     =find(mvRotorGroupInternal.back().mvRotatedAtomList.begin(),
                           mvRotorGroupInternal.back().mvRotatedAtomList.end(),*atom2);
               if(  (check==mvRotorGroupInternal.back().mvRotatedAtomList.end())
                  ||(mvRotorGroupInternal.back().mvRotatedAtomList.size()<3)
                  ||(mvRotorGroupInternal.back().mvRotatedAtomList.size()>=((mvpAtom.size()+1)/2)))
               {
                  mvRotorGroupInternal.pop_back();
               }
               else
               {
                  mvRotorGroupInternal.back().mvRotatedAtomList.erase(*atom1);
                  mvRotorGroupInternal.back().mvRotatedAtomList.erase(*atom2);
               }
            }
         }
      }
   }
   #endif
   
   // Remove identical groups
   for(unsigned int i=1;i<=3;++i)
   {
      list<RotorGroup> *pRotorGroup1;
      switch(i)
      {
         case 1: pRotorGroup1=&mvRotorGroupTorsion;break;
         case 2: pRotorGroup1=&mvRotorGroupTorsionSingleChain;break;
         case 3: pRotorGroup1=&mvRotorGroupInternal;break;
      }
      for(list<RotorGroup>::iterator pos1=pRotorGroup1->begin();
          pos1!=pRotorGroup1->end();++pos1)
      {
         for(unsigned int j=i;j<=3;++j)
         {
            list<RotorGroup> *pRotorGroup2;
            switch(j)
            {
               case 1: pRotorGroup2=&mvRotorGroupTorsion;break;
               case 2: pRotorGroup2=&mvRotorGroupTorsionSingleChain;break;
               case 3: pRotorGroup2=&mvRotorGroupInternal;break;
            }
            for(list<RotorGroup>::iterator pos2=pRotorGroup2->begin();
                pos2!=pRotorGroup2->end();)
            {
               if(pos2==pos1) {++pos2;continue;}
               if((  ((pos1->mpAtom1 == pos2->mpAtom1) && (pos1->mpAtom2 == pos2->mpAtom2))
                   ||((pos1->mpAtom2 == pos2->mpAtom1) && (pos1->mpAtom1 == pos2->mpAtom2))) 
                  &&pos1->mvRotatedAtomList.size() == pos2->mvRotatedAtomList.size())
               {
                  bool ident=true;
                  for(set<MolAtom*>::const_iterator pos=pos1->mvRotatedAtomList.begin();
                      pos!=pos1->mvRotatedAtomList.end();++pos)
                  {
                     set<MolAtom*>::const_iterator tmp=pos2->mvRotatedAtomList.find(*pos);
                     if(tmp == pos2->mvRotatedAtomList.end())
                     {
                        ident=false;
                        break;
                     }
                  }
                  if(ident)
                  {
                     #if 0
                     cout<<"Identical groups:"<<endl;
                     cout<<"    G1:"
                         <<pos1->mpAtom1->GetName()<<"-"
                         <<pos1->mpAtom2->GetName()<<" : ";
                     for(set<MolAtom*>::iterator pos=pos1->mvRotatedAtomList.begin();
                         pos!=pos1->mvRotatedAtomList.end();++pos)
                        cout<<(*pos)->GetName()<<"  ";
                     cout<<endl;
                     cout<<"    G2:"
                         <<pos2->mpAtom1->GetName()<<"-"
                         <<pos2->mpAtom2->GetName()<<" : ";
                     for(set<MolAtom*>::iterator pos=pos2->mvRotatedAtomList.begin();
                         pos!=pos2->mvRotatedAtomList.end();++pos)
                        cout<<(*pos)->GetName()<<"  ";
                     cout<<endl;
                     #endif
                     pos2=pRotorGroup2->erase(pos2);
                  }
                  else ++pos2;
               }
               else ++pos2;
            }
         }
      }
   }
   // Remove all rotations which break restraints and therefore are not "free torsion"
   this->SaveParamSet(mLocalParamSet);
   const REAL llk0=this->GetLogLikelihood();
   for(unsigned int i=1;i<=3;++i)
   {
      list<RotorGroup> *pRotorGroup1;
      switch(i)
      {
         case 1: pRotorGroup1=&mvRotorGroupTorsion;break;
         case 2: pRotorGroup1=&mvRotorGroupTorsionSingleChain;break;
         case 3: pRotorGroup1=&mvRotorGroupInternal;break;
      }
      for(list<RotorGroup>::iterator pos=pRotorGroup1->begin();
          pos!=pRotorGroup1->end();)
      {
         REAL llk=0;
         for(unsigned int j=0;j<36;++j)
         {
            const REAL angle=(REAL)j*M_PI/36.;
            this->RotateAtomGroup(*(pos->mpAtom1),*(pos->mpAtom2),
                                  pos->mvRotatedAtomList,angle);
            // use fabs in case we are not starting from the minimum of a restraint..
            llk += fabs(this->GetLogLikelihood() - llk0);
            this->RestoreParamSet(mLocalParamSet);
         }
         #ifdef __DEBUG__
         switch(i)
         {
            case 1: cout<<"Rotation Group around bond :";break;
            case 2: cout<<"Rotation Group (single chain) around bond :";break;
            case 3: cout<<"Rotation Group (internal) between :";break;
         }
         cout <<pos->mpAtom1->GetName()<<"-"
             <<pos->mpAtom2->GetName()<<" : ";
         for(set<MolAtom *>::iterator pos1=pos->mvRotatedAtomList.begin();
             pos1!=pos->mvRotatedAtomList.end();++pos1)
            cout<<(*pos1)->GetName()<<"  ";
         cout<<"   <d(LLK)>="<< llk/36.;
         #endif
         if((llk/50.)>100.)
         {
            pos = pRotorGroup1->erase(pos);
            //cout <<" -> NOT a free torsion"<<endl;
         }
         else ++pos;
         //else
         //   cout <<" -> free torsion"<<endl;
      }
   }
   //cout<<endl;
   
   // Label free torsions
   for(vector<MolBond*>::iterator pos=mvpBond.begin();pos!=mvpBond.end();++pos)
      (*pos)->SetFreeTorsion(false);
   for(list<RotorGroup>::iterator pos=mvRotorGroupTorsion.begin();
       pos!=mvRotorGroupTorsion.end();++pos)
   {
      vector<MolBond*>::iterator  bd=this->FindBond((*pos->mpAtom1),(*pos->mpAtom2));
      if(bd!=mvpBond.end()) (*bd)->SetFreeTorsion(true);
   }
   
   mClockRotorGroup.Click();
   VFN_DEBUG_EXIT("Molecule::BuildRotorGroup()",5)
}

void Molecule::BuildStretchModeBondLength()
{
   #if 0
   if(  (mClockStretchModeBondLength>mClockBondList)
      &&(mClockStretchModeBondLength>mClockAtomList)
      &&(mClockStretchModeBondLength>mClockBondAngleList)
      &&(mClockStretchModeBondLength>mClockDihedralAngleList)) return;
   #endif
   VFN_DEBUG_ENTRY("Molecule::BuildStretchModeBondLength()",7)
   this->BuildConnectivityTable();
   TAU_PROFILE("Molecule::BuildStretchModeBondLength()","void ()",TAU_DEFAULT);
   TAU_PROFILE_TIMER(timer1,"Molecule::BuildStretchModeBondLength 1","", TAU_FIELD);
   TAU_PROFILE_TIMER(timer2,"Molecule::BuildStretchModeBondLength 2","", TAU_FIELD);
   TAU_PROFILE_TIMER(timer3,"Molecule::BuildStretchModeBondLength 3","", TAU_FIELD);
   TAU_PROFILE_TIMER(timer4,"Molecule::BuildStretchModeBondLength 4","", TAU_FIELD);
   TAU_PROFILE_TIMER(timer5,"Molecule::BuildStretchModeBondLength 5","", TAU_FIELD);
   mvStretchModeBondLength.clear();
   // Build list of atoms moved when stretching a bond length. Only keep the group
   // of atoms on the smaller side.
   TAU_PROFILE_START(timer1);
   for(unsigned long i=0;i<mvpBond.size();++i)
   {
      //if((mFlexModel.GetChoice()!=0)&&(false==mvpBond[i]->IsFreeTorsion())) continue;
      MolAtom* const atom1=&(mvpBond[i]->GetAtom1());
      MolAtom* const atom2=&(mvpBond[i]->GetAtom2());
      for(unsigned int j=1;j<=2;++j)
      {
         const set<MolAtom*> *pConn;
         if(j==1) pConn=&(mConnectivityTable[atom1]);
         else pConn=&(mConnectivityTable[atom2]);
         
         if(j==1)
            mvStretchModeBondLength.push_back(StretchModeBondLength(mvpBond[i]->GetAtom2(),
                                                                    mvpBond[i]->GetAtom1(),
                                                                    mvpBond[i]));
         else
            mvStretchModeBondLength.push_back(StretchModeBondLength(mvpBond[i]->GetAtom1(),
                                                                    mvpBond[i]->GetAtom2(),
                                                                    mvpBond[i]));
         if(j==1) mvStretchModeBondLength.back().mvTranslatedAtomList.insert(atom1);
         if(j==2) mvStretchModeBondLength.back().mvTranslatedAtomList.insert(atom2);
         
         for(set<MolAtom*>::const_iterator pos=pConn->begin();pos!=pConn->end();++pos)
         {
            if((j==1)&&(*pos==atom2)) continue;
            if((j==2)&&(*pos==atom1)) continue;
            ExpandAtomGroupRecursive(*pos,mConnectivityTable,
                                     mvStretchModeBondLength.back().mvTranslatedAtomList);
         }
         const unsigned long ct1=mvStretchModeBondLength.back().mvTranslatedAtomList.count(atom1),
                             ct2=mvStretchModeBondLength.back().mvTranslatedAtomList.count(atom2);
         if( ((j==1)&&(ct2>0)) || ((j==2)&&(ct1>0)) )
         {
            // We have found a ring. No use looking at the other side.
            // :TODO: handle this properly..
            mvStretchModeBondLength.pop_back();
            break;
         }
         if(  (mvStretchModeBondLength.back().mvTranslatedAtomList.size()>((mvpAtom.size()+1)/2))
            ||(mvStretchModeBondLength.back().mvTranslatedAtomList.size()==0))
         {
            #ifdef __DEBUG__
            cout<<"Rejecting StretchModeBondLength ";mvStretchModeBondLength.back().Print(cout);cout<<endl;
            #endif
            mvStretchModeBondLength.pop_back();
         }
         else if(mvStretchModeBondLength.back().mvTranslatedAtomList.size()==((mvpAtom.size()+1)/2))
                  break;//we translate exactly half of the atoms, so skip the other half
      }
   }
   TAU_PROFILE_STOP(timer1);
   
   // find rigid groups broken by each mode
   for(list<StretchModeBondLength>::iterator pos=mvStretchModeBondLength.begin();
       pos!=mvStretchModeBondLength.end();)
   {
      TAU_PROFILE_START(timer5);
      bool keep=true;
      for(vector<RigidGroup*>::const_iterator group=mvRigidGroup.begin();
          group!=mvRigidGroup.end();++group)
      {
         unsigned long ct=0;
         for(set<MolAtom *>::const_iterator at=(*group)->begin();at!=(*group)->end();++at)
            ct += pos->mvTranslatedAtomList.count(*at);
         if((ct>0)&&(ct!=(*group)->size()))
         {
            keep=false;
            #ifdef __DEBUG__
            cout<<"       Breaks Rigid Group:";
            for(set<MolAtom *>::const_iterator at=(*group)->begin();at!=(*group)->end();++at)
               cout<<(*at)->GetName()<<" ";
            cout<<endl;
            #endif
            break;
         }
      }
      if(keep) ++pos;
      else pos=mvStretchModeBondLength.erase(pos);
      TAU_PROFILE_STOP(timer5);
   }
   // Generate 5 completely random atomic positions
      this->SaveParamSet(mLocalParamSet);
      unsigned long paramSetRandom[5];
      for(unsigned long i=0;i<5;++i)
      {
         for(vector<MolAtom*>::iterator pos=mvpAtom.begin();pos!=mvpAtom.end();++pos)
         {
            (*pos)->SetX(100.*rand()/(REAL) RAND_MAX);
            (*pos)->SetY(100.*rand()/(REAL) RAND_MAX);
            (*pos)->SetZ(100.*rand()/(REAL) RAND_MAX);
         }
         paramSetRandom[i]=this->CreateParamSet();
      }
   // find bond lengths broken by each mode
   for(list<StretchModeBondLength>::iterator pos=mvStretchModeBondLength.begin();
       pos!=mvStretchModeBondLength.end();)
   {
      TAU_PROFILE_START(timer2);
      bool keep=true;
      pos->mvpBrokenBond.clear();
      for(vector<MolBond*>::const_iterator r=mvpBond.begin();r!=mvpBond.end();++r)
      {
         unsigned int ct=0;
         for(set<MolAtom *>::const_iterator at=pos->mvTranslatedAtomList.begin();
             at!=pos->mvTranslatedAtomList.end();++at)
         {
            if(*at==&((*r)->GetAtom1())) ct++;
            if(*at==&((*r)->GetAtom2())) ct++;
         }
         // If we moved either both or non of the bond atom, the bond length is unchanged.
         if((ct!=0)&&(ct !=2)) pos->mvpBrokenBond.insert(make_pair(*r,0));
      }
      if(mFlexModel.GetChoice()==2)
      {
         int nb=pos->mvpBrokenBond.size();
         if(pos->mpBond!=0) nb -= 1;
         if(nb>0) keep=false;
      }
      if(keep) ++pos;
      else pos=mvStretchModeBondLength.erase(pos);
      TAU_PROFILE_STOP(timer2);
   }
   // find bond angles broken by each mode
   for(list<StretchModeBondLength>::iterator pos=mvStretchModeBondLength.begin();
       pos!=mvStretchModeBondLength.end();)
   {
      TAU_PROFILE_START(timer3);
      bool keep=true;
      pos->mvpBrokenBondAngle.clear();
      for(vector<MolBondAngle*>::const_iterator r=mvpBondAngle.begin();r!=mvpBondAngle.end();++r)
      {
         unsigned int ct=0;
         for(set<MolAtom *>::const_iterator at=pos->mvTranslatedAtomList.begin();
             at!=pos->mvTranslatedAtomList.end();++at)
         {
            if(*at==&((*r)->GetAtom1())) ct++;
            if(*at==&((*r)->GetAtom2())) ct++;
            if(*at==&((*r)->GetAtom3())) ct++;
         }
         bool broken=true;
         if((ct==0)||(ct==3)) broken=false;
         if(broken)
         {// Make sure with derivatives
            REAL d=0;
            for(unsigned long i=0;i<5;++i)
            {
               this->RestoreParamSet(paramSetRandom[i]);
               pos->CalcDeriv(false);
               (*r)->GetLogLikelihood(true,true);
               d += abs((*r)->GetDeriv(pos->mDerivXYZ));
               if(d>0.01) break;
            }
            if(abs(d)<=0.01) broken=false;
         }
         if(broken) pos->mvpBrokenBondAngle.insert(make_pair(*r,0));
      }
      if(mFlexModel.GetChoice()==2)
      {
         if(pos->mvpBrokenBondAngle.size()>0) keep=false;
      }
      if(keep) ++pos;
      else pos=mvStretchModeBondLength.erase(pos);
      TAU_PROFILE_STOP(timer3);
   }
   // find dihedral angles broken by each mode
   for(list<StretchModeBondLength>::iterator pos=mvStretchModeBondLength.begin();
       pos!=mvStretchModeBondLength.end();)
   {
      TAU_PROFILE_START(timer4);
      bool keep=true;
      pos->mvpBrokenDihedralAngle.clear();
      for(vector<MolDihedralAngle*>::const_iterator r=mvpDihedralAngle.begin();r!=mvpDihedralAngle.end();++r)
      {
         unsigned int ct=0;
         for(set<MolAtom *>::const_iterator at=pos->mvTranslatedAtomList.begin();
             at!=pos->mvTranslatedAtomList.end();++at)
         {
            if(*at==&((*r)->GetAtom1())) ct++;
            if(*at==&((*r)->GetAtom2())) ct++;
            if(*at==&((*r)->GetAtom3())) ct++;
            if(*at==&((*r)->GetAtom4())) ct++;
         }
         bool broken=true;
         if((ct==0)||(ct==4)) broken=false;
         if(broken)
         {// Make sure with derivatives
            REAL d=0;
            for(unsigned long i=0;i<5;++i)
            {
               this->RestoreParamSet(paramSetRandom[i]);
               pos->CalcDeriv(false);
               (*r)->GetLogLikelihood(true,true);
               d += abs((*r)->GetDeriv(pos->mDerivXYZ));
               if(d>0.01) break;
            }
            if(abs(d)<=0.01) broken=false;
         }
         if(broken) pos->mvpBrokenDihedralAngle.insert(make_pair(*r,0));
      }
      if(mFlexModel.GetChoice()==2)
      {
         if(pos->mvpBrokenDihedralAngle.size()>0) keep=false;
      }
      if(keep) ++pos;
      else pos=mvStretchModeBondLength.erase(pos);
      TAU_PROFILE_STOP(timer4);
   }
   this->RestoreParamSet(mLocalParamSet);
   for(unsigned long i=0;i<5;++i) this->ClearParamSet(paramSetRandom[i]);
   #if 1//def __DEBUG__
   cout<<"List of Bond Length stretch modes"<<endl;
   for(list<StretchModeBondLength>::const_iterator pos=mvStretchModeBondLength.begin();
       pos!=mvStretchModeBondLength.end();++pos)
   {
      cout<<"   Bond:"<<pos->mpAtom0->GetName()<<"-"<<pos->mpAtom1->GetName()<<", Translated Atoms:  ";
      for(set<MolAtom*>::const_iterator atom=pos->mvTranslatedAtomList.begin();
          atom!=pos->mvTranslatedAtomList.end();++atom)
      {
         cout<<(*atom)->GetName()<<",";
      }
      if(pos->mpBond!=0) cout<< " ; restrained to length="<<pos->mpBond->GetLength0()
                             <<", sigma="<<pos->mpBond->GetLengthSigma()
                             <<", delta="<<pos->mpBond->GetLengthDelta();
      if(pos->mvpBrokenBond.size()>0)
      {
         cout<<endl<<"       Broken bonds:";
         for(map<const MolBond*,REAL>::const_iterator bond=pos->mvpBrokenBond.begin();
             bond!=pos->mvpBrokenBond.end();++bond)
            cout<<bond->first->GetName()<<", ";
      }
      if(pos->mvpBrokenBondAngle.size()>0)
      {
         cout<<endl<<"       Broken bond angles:";
         for(map<const MolBondAngle*,REAL>::const_iterator angle=pos->mvpBrokenBondAngle.begin();
             angle!=pos->mvpBrokenBondAngle.end();++angle)
            cout<<angle->first->GetName()<<", ";
      }
      if(pos->mvpBrokenDihedralAngle.size()>0)
      {
         cout<<endl<<"       Broken dihedral angles:";
         for(map<const MolDihedralAngle*,REAL>::const_iterator 
             angle=pos->mvpBrokenDihedralAngle.begin();
             angle!=pos->mvpBrokenDihedralAngle.end();++angle)
            cout<<angle->first->GetName()<<", ";
      }
      cout<<endl;
   }
   #endif
   mClockStretchModeBondLength.Click();
   VFN_DEBUG_EXIT("Molecule::BuildStretchModeBondLength()",7)
}

void Molecule::BuildStretchModeBondAngle()
{
   #if 0
   if(  (mClockStretchModeBondAngle>mClockBondList)
      &&(mClockStretchModeBondAngle>mClockAtomList)
      &&(mClockStretchModeBondAngle>mClockBondAngleList)
      &&(mClockStretchModeBondAngle>mClockDihedralAngleList)) return;
   #endif
   VFN_DEBUG_ENTRY("Molecule::BuildStretchModeBondAngle()",10)
   this->BuildConnectivityTable();
   TAU_PROFILE("Molecule::BuildStretchModeBondAngle()","void ()",TAU_DEFAULT);
   TAU_PROFILE_TIMER(timer1,"Molecule::BuildStretchModeBondAngle 1","", TAU_FIELD);
   TAU_PROFILE_TIMER(timer2,"Molecule::BuildStretchModeBondAngle 2","", TAU_FIELD);
   TAU_PROFILE_TIMER(timer3,"Molecule::BuildStretchModeBondAngle 3","", TAU_FIELD);
   TAU_PROFILE_TIMER(timer4,"Molecule::BuildStretchModeBondAngle 4","", TAU_FIELD);
   TAU_PROFILE_TIMER(timer5,"Molecule::BuildStretchModeBondAngle 5","", TAU_FIELD);
   mvStretchModeBondAngle.clear();
   // Build list of atoms moved when stretching a bond angle. Only keep the group
   // of atoms on the smaller side.
   TAU_PROFILE_START(timer1);
   for(unsigned long i=0;i<mvpAtom.size();++i)
   {
      //if((mFlexModel.GetChoice()!=0)&&(false==mvpBond[i]->IsFreeTorsion())) continue;
      set<MolAtom*> *pConn0=&(mConnectivityTable[mvpAtom[i]]);
      if(pConn0->size()<2) continue;
      for(set<MolAtom*>::const_iterator pos1=pConn0->begin();pos1!=pConn0->end();++pos1)
      {
         set<MolAtom*>::const_iterator pos2=pos1;
         pos2++;
         for(;pos2!=pConn0->end();++pos2)
         {
            VFN_DEBUG_MESSAGE("Molecule::BuildStretchModeBondAngle():"<<i<<","<<*pos1<<","<<*pos2,10)
            //Do we have a bond angle restraint corresponding to these atoms ?
            MolBondAngle *pMolBondAngle=0;
            for(vector<MolBondAngle*>::const_iterator pos=mvpBondAngle.begin();pos!=mvpBondAngle.end();++pos)
            {
               if(&((*pos)->GetAtom2())==mvpAtom[i])
               {
                  if(  ((&((*pos)->GetAtom1())==*pos1)&&(&((*pos)->GetAtom3())==*pos2))
                     ||((&((*pos)->GetAtom1())==*pos2)&&(&((*pos)->GetAtom3())==*pos1)))
                  {
                     pMolBondAngle=*pos;
                     break;
                  }
               }
            }
            for(unsigned int j=1;j<=2;++j)
            {
               const set<MolAtom*> *pConn;
               if(j==1)
               {
                  pConn=&(mConnectivityTable[*pos1]);
                  mvStretchModeBondAngle.push_back(StretchModeBondAngle(**pos2,
                                                                        *mvpAtom[i],
                                                                        **pos1,
                                                                        pMolBondAngle));
                  mvStretchModeBondAngle.back().mvRotatedAtomList.insert(*pos1);
               }
               else
               {
                  pConn=&(mConnectivityTable[*pos2]);
                  mvStretchModeBondAngle.push_back(StretchModeBondAngle(**pos1,
                                                                        *mvpAtom[i],
                                                                        **pos2,
                                                                        pMolBondAngle));
                  mvStretchModeBondAngle.back().mvRotatedAtomList.insert(*pos2);
               }
               mvStretchModeBondAngle.back().mvRotatedAtomList.insert(mvpAtom[i]);

               for(set<MolAtom*>::const_iterator pos=pConn->begin();pos!=pConn->end();++pos)
               {
                  if(*pos==mvpAtom[i]) continue;
                  ExpandAtomGroupRecursive(*pos,mConnectivityTable,
                                           mvStretchModeBondAngle.back().mvRotatedAtomList);
               }
               //if(j==1)mvStretchModeBondAngle.back().mvRotatedAtomList.erase(*pos2);
               //if(j==2)mvStretchModeBondAngle.back().mvRotatedAtomList.erase(*pos1);
               mvStretchModeBondAngle.back().mvRotatedAtomList.erase(mvpAtom[i]);
               
               if(  (mvStretchModeBondAngle.back().mvRotatedAtomList.size()>=((mvpAtom.size()+1)/2))
                  ||(mvStretchModeBondAngle.back().mvRotatedAtomList.size()==0))
               {
                  #ifdef __DEBUG__
                  cout<<"Rejecting StretchModeBondAngle ";mvStretchModeBondAngle.back().Print(cout);cout<<endl;
                  #endif
                  mvStretchModeBondAngle.pop_back();
               }
               else
               {
                  if((j==1) &&(mvStretchModeBondAngle.back().mvRotatedAtomList.find(*pos2)
                               !=mvStretchModeBondAngle.back().mvRotatedAtomList.end()))
                  {
                     #if 1//def __DEBUG__
                     cout<<"Rejecting StretchModeBondAngle (ring) ";mvStretchModeBondAngle.back().Print(cout);cout<<endl;
                     #endif
                     mvStretchModeBondAngle.pop_back();
                  }
                  if((j==2) &&(mvStretchModeBondAngle.back().mvRotatedAtomList.find(*pos1)
                               !=mvStretchModeBondAngle.back().mvRotatedAtomList.end()))
                  {
                     #if 1//def __DEBUG__
                     cout<<"Rejecting StretchModeBondAngle (ring) ";mvStretchModeBondAngle.back().Print(cout);cout<<endl;
                     #endif
                     mvStretchModeBondAngle.pop_back();
                  }
               }
               
            }
         }
      }
   }
   TAU_PROFILE_STOP(timer1);
   // find rigid groups broken by each mode
   for(list<StretchModeBondAngle>::iterator pos=mvStretchModeBondAngle.begin();
       pos!=mvStretchModeBondAngle.end();)
   {
      TAU_PROFILE_START(timer5);
      bool keep=true;
      for(vector<RigidGroup*>::const_iterator group=mvRigidGroup.begin();
          group!=mvRigidGroup.end();++group)
      {
         unsigned long ct=0;
         for(set<MolAtom *>::const_iterator at=(*group)->begin();at!=(*group)->end();++at)
            ct += pos->mvRotatedAtomList.count(*at);
         if(ct>0)
         {
            // Add the origin atom, which does not move relatively to the rotated atoms
            ct += (*group)->count(pos->mpAtom1);
            if(ct!=(*group)->size())
            {
               keep=false;
               #ifdef __DEBUG__
               pos->Print(cout);
               cout<<"       Breaks Rigid Group:";
               for(set<MolAtom *>::const_iterator at=(*group)->begin();at!=(*group)->end();++at)
                  cout<<(*at)->GetName()<<" ";
               cout<<endl;
               #endif
               break;
            }
         }
      }
      if(keep) ++pos;
      else pos=mvStretchModeBondAngle.erase(pos);
      TAU_PROFILE_STOP(timer5);
   }
   // Generate 5 completely random atomic positions
      this->SaveParamSet(mLocalParamSet);
      unsigned long paramSetRandom[5];
      for(unsigned long i=0;i<5;++i)
      {
         for(vector<MolAtom*>::iterator pos=mvpAtom.begin();pos!=mvpAtom.end();++pos)
         {
            (*pos)->SetX(100.*rand()/(REAL) RAND_MAX);
            (*pos)->SetY(100.*rand()/(REAL) RAND_MAX);
            (*pos)->SetZ(100.*rand()/(REAL) RAND_MAX);
         }
         paramSetRandom[i]=this->CreateParamSet();
      }
   // find bond lengths broken by each mode
   for(list<StretchModeBondAngle>::iterator pos=mvStretchModeBondAngle.begin();
       pos!=mvStretchModeBondAngle.end();)
   {
      TAU_PROFILE_START(timer2);
      bool keep=true;
      pos->mvpBrokenBond.clear();
      for(vector<MolBond*>::const_iterator r=mvpBond.begin();r!=mvpBond.end();++r)
      {
         unsigned int ct=0;
         for(set<MolAtom *>::const_iterator at=pos->mvRotatedAtomList.begin();
             at!=pos->mvRotatedAtomList.end();++at)
         {
            if(*at==&((*r)->GetAtom1())) ct++;
            if(*at==&((*r)->GetAtom2())) ct++;
         }
         bool broken=true;
         // If we moved either both or non of the bond atom, the bond length is unchanged.
         if((ct==0)||(ct==2)) broken=false;
         if(broken)
         {// Make sure with derivatives
            REAL d=0;
            for(unsigned long i=0;i<5;++i)
            {
               this->RestoreParamSet(paramSetRandom[i]);
               pos->CalcDeriv(false);
               (*r)->GetLogLikelihood(true,true);
               d += abs((*r)->GetDeriv(pos->mDerivXYZ));
               if(d>0.01) break;
            }
            if(abs(d)<=0.01) broken=false;
         }
         if(broken) pos->mvpBrokenBond.insert(make_pair(*r,0));
      }
      if(mFlexModel.GetChoice()==2)
      {
         if(pos->mvpBrokenBond.size()>0) keep=false;
      }
      if(keep) ++pos;
      else pos=mvStretchModeBondAngle.erase(pos);
      TAU_PROFILE_STOP(timer2);
   }
   // find bond angles broken by each mode
   for(list<StretchModeBondAngle>::iterator pos=mvStretchModeBondAngle.begin();
       pos!=mvStretchModeBondAngle.end();)
   {
      TAU_PROFILE_START(timer3);
      bool keep=true;
      pos->mvpBrokenBondAngle.clear();
      for(vector<MolBondAngle*>::const_iterator r=mvpBondAngle.begin();r!=mvpBondAngle.end();++r)
      {
         unsigned int ct=0;
         for(set<MolAtom *>::const_iterator at=pos->mvRotatedAtomList.begin();
             at!=pos->mvRotatedAtomList.end();++at)
         {
            if(*at==&((*r)->GetAtom1())) ct++;
            if(*at==&((*r)->GetAtom2())) ct++;
            if(*at==&((*r)->GetAtom3())) ct++;
         }
         bool broken=true;
         if(ct==0) broken=false;
         if(ct==3) broken=false;
         if(broken)
         {// Make sure with derivatives
            REAL d=0;
            for(unsigned long i=0;i<5;++i)
            {
               this->RestoreParamSet(paramSetRandom[i]);
               pos->CalcDeriv(false);
               (*r)->GetLogLikelihood(true,true);
               d += abs((*r)->GetDeriv(pos->mDerivXYZ));
               if(d>0.01) break;
            }
            if(abs(d)<=0.01) broken=false;
         }
         if(broken) pos->mvpBrokenBondAngle.insert(make_pair(*r,0));
      }
      if(mFlexModel.GetChoice()==2)
      {
         int nb=pos->mvpBrokenBond.size();
         if(pos->mpBondAngle!=0) nb -= 1;
         if(nb>0) keep=false;
      }
      if(keep) ++pos;
      else pos=mvStretchModeBondAngle.erase(pos);
      TAU_PROFILE_STOP(timer3);
   }
   // find dihedral angles broken by each mode
   for(list<StretchModeBondAngle>::iterator pos=mvStretchModeBondAngle.begin();
       pos!=mvStretchModeBondAngle.end();)
   {
      TAU_PROFILE_START(timer4);
      bool keep=true;
      pos->mvpBrokenDihedralAngle.clear();
      for(vector<MolDihedralAngle*>::const_iterator r=mvpDihedralAngle.begin();r!=mvpDihedralAngle.end();++r)
      {
         unsigned int ct=0;
         for(set<MolAtom *>::const_iterator at=pos->mvRotatedAtomList.begin();
             at!=pos->mvRotatedAtomList.end();++at)
         {
            if(*at==&((*r)->GetAtom1())) ct++;
            if(*at==&((*r)->GetAtom2())) ct++;
            if(*at==&((*r)->GetAtom3())) ct++;
            if(*at==&((*r)->GetAtom4())) ct++;
         }
         bool broken=true;
         if(ct==0) broken=false;
         if(ct==4) broken=false;
         if(broken)
         {// Make sure with derivatives
            REAL d=0;
            for(unsigned long i=0;i<5;++i)
            {
               this->RestoreParamSet(paramSetRandom[i]);
               pos->CalcDeriv(false);
               (*r)->GetLogLikelihood(true,true);
               d += abs((*r)->GetDeriv(pos->mDerivXYZ));
               if(d>0.01) break;
            }
            if(abs(d)<=0.01) broken=false;
         }
         if(broken) pos->mvpBrokenDihedralAngle.insert(make_pair(*r,0));
      }
      if(mFlexModel.GetChoice()==2)
      {
         if(pos->mvpBrokenDihedralAngle.size()>0) keep=false;
      }
      if(keep) ++pos;
      else pos=mvStretchModeBondAngle.erase(pos);
      TAU_PROFILE_STOP(timer4);
   }
   this->RestoreParamSet(mLocalParamSet);
   for(unsigned long i=0;i<5;++i) this->ClearParamSet(paramSetRandom[i]);
   #ifdef __DEBUG__
   cout<<"List of Bond Angle stretch modes"<<endl;
   for(list<StretchModeBondAngle>::const_iterator pos=mvStretchModeBondAngle.begin();
       pos!=mvStretchModeBondAngle.end();++pos)
   {
      cout<<"   Angle:"<<pos->mpAtom0->GetName()<<"-"
                       <<pos->mpAtom1->GetName()<<"-"
                       <<pos->mpAtom2->GetName()<<", Rotated Atoms:  ";
      for(set<MolAtom*>::const_iterator atom=pos->mvRotatedAtomList.begin();
          atom!=pos->mvRotatedAtomList.end();++atom)
      {
         cout<<(*atom)->GetName()<<",";
      }
      if(pos->mpBondAngle!=0) cout<< " ; restrained to angle="<<pos->mpBondAngle->GetAngle0()*RAD2DEG
                                  <<"�, sigma="<<pos->mpBondAngle->GetAngleSigma()*RAD2DEG
                                  <<"�, delta="<<pos->mpBondAngle->GetAngleDelta()*RAD2DEG<<"�";
      if(pos->mvpBrokenBond.size()>0)
      {
         cout<<endl<<"       Broken bonds:";
         for(map<const MolBond*,REAL>::const_iterator bond=pos->mvpBrokenBond.begin();
             bond!=pos->mvpBrokenBond.end();++bond)
            cout<<bond->first->GetName()<<", ";
      }
      if(pos->mvpBrokenBondAngle.size()>0)
      {
         cout<<endl<<"       Broken bond angles:";
         for(map<const MolBondAngle*,REAL>::const_iterator angle=pos->mvpBrokenBondAngle.begin();
             angle!=pos->mvpBrokenBondAngle.end();++angle)
            cout<<angle->first->GetName()<<", ";
      }
      if(pos->mvpBrokenDihedralAngle.size()>0)
      {
         cout<<endl<<"       Broken dihedral angles:";
         for(map<const MolDihedralAngle*,REAL>::const_iterator 
             angle=pos->mvpBrokenDihedralAngle.begin();
             angle!=pos->mvpBrokenDihedralAngle.end();++angle)
            cout<<angle->first->GetName()<<", ";
      }
      cout<<endl;
   }
   #endif
   mClockStretchModeBondAngle.Click();
   VFN_DEBUG_EXIT("Molecule::BuildStretchModeBondAngle()",10)
}

void Molecule::BuildStretchModeTorsion()
{
   #if 0
   if(  (mClockStretchModeTorsion>mClockBondList)
      &&(mClockStretchModeTorsion>mClockAtomList)
      &&(mClockStretchModeTorsion>mClockBondAngleList)
      &&(mClockStretchModeTorsion>mClockDihedralAngleList)) return;
   #endif
   VFN_DEBUG_ENTRY("Molecule::BuildStretchModeTorsion()",7)
   TAU_PROFILE("Molecule::BuildStretchModeTorsion()","void ()",TAU_DEFAULT);
   TAU_PROFILE_TIMER(timer1,"Molecule::BuildStretchModeTorsion 1","", TAU_FIELD);
   TAU_PROFILE_TIMER(timer2,"Molecule::BuildStretchModeTorsion 2","", TAU_FIELD);
   TAU_PROFILE_TIMER(timer3,"Molecule::BuildStretchModeTorsion 3","", TAU_FIELD);
   TAU_PROFILE_TIMER(timer4,"Molecule::BuildStretchModeTorsion 4","", TAU_FIELD);
   TAU_PROFILE_TIMER(timer5,"Molecule::BuildStretchModeTorsion 5","", TAU_FIELD);
   this->BuildConnectivityTable();
   mvStretchModeTorsion.clear();
   // Build list of atoms moved when changing the angle. Only keep the group
   // of atoms on the smaller side.
   for(unsigned long i=0;i<mvpBond.size();++i)
   {
      //if((mFlexModel.GetChoice()!=0)&&(false==mvpBond[i]->IsFreeTorsion())) continue;
      MolAtom* const atom1=&(mvpBond[i]->GetAtom1());
      MolAtom* const atom2=&(mvpBond[i]->GetAtom2());
      for(unsigned int j=1;j<=2;++j)
      {
         const set<MolAtom*> *pConn;
         if(j==1) pConn=&(mConnectivityTable[atom1]);
         else pConn=&(mConnectivityTable[atom2]);
         mvStretchModeTorsion.push_back(StretchModeTorsion(*atom1,*atom2,0));
         mvStretchModeTorsion.back().mvRotatedAtomList.insert(atom1);
         mvStretchModeTorsion.back().mvRotatedAtomList.insert(atom2);
         
         for(set<MolAtom*>::const_iterator pos=pConn->begin();pos!=pConn->end();++pos)
         {
            if((*pos==atom2)||(*pos==atom1)) continue;
            ExpandAtomGroupRecursive(*pos,mConnectivityTable,
                                     mvStretchModeTorsion.back().mvRotatedAtomList);
         }
         mvStretchModeTorsion.back().mvRotatedAtomList.erase(atom1);
         mvStretchModeTorsion.back().mvRotatedAtomList.erase(atom2);
         
         for(vector<MolDihedralAngle*>::const_iterator dih=mvpDihedralAngle.begin();dih!=mvpDihedralAngle.end();++dih)
         {
            // :TODO: There are some other weird cases to take into account,
            // for restraints with atoms *not* connected to another
            // More generally, should check the list of atoms rotated.
            if(  ((&((*dih)->GetAtom2())==atom1) && (&((*dih)->GetAtom3())==atom2))
               ||((&((*dih)->GetAtom3())==atom1) && (&((*dih)->GetAtom2())==atom2)))
            {
               const unsigned long ct1=mvStretchModeTorsion.back().mvRotatedAtomList.count(&((*dih)->GetAtom1()));
               const unsigned long ct4=mvStretchModeTorsion.back().mvRotatedAtomList.count(&((*dih)->GetAtom4()));

               if((ct1+ct4)==1)// One of the atom is rotated, not the other
               {
                  mvStretchModeTorsion.back().mpDihedralAngle=*dih;
                  //:TODO: Check sense of rotation !
                  if(ct4==1)
                  {
                     mvStretchModeTorsion.back().mpAtom1=&((*dih)->GetAtom2());
                     mvStretchModeTorsion.back().mpAtom2=&((*dih)->GetAtom3());
                  }
                  else
                  {
                     mvStretchModeTorsion.back().mpAtom1=&((*dih)->GetAtom3());
                     mvStretchModeTorsion.back().mpAtom2=&((*dih)->GetAtom2());
                  }
               }
            }
         }
         
         if(mvStretchModeTorsion.size()>1)
         {//Duplicate ?
            // Does not work with a const_reverse_iterator ?
            // http://gcc.gnu.org/bugzilla/show_bug.cgi?id=11729
            list<StretchModeTorsion>::reverse_iterator mode=mvStretchModeTorsion.rbegin();
            ++mode;
            for(;mode!=mvStretchModeTorsion.rend();++mode)
            {
               if( (  ((mode->mpAtom1==atom1)&&(mode->mpAtom2==atom2)) 
                    ||((mode->mpAtom1==atom2)&&(mode->mpAtom2==atom1)))
                  &&(mode->mvRotatedAtomList==mvStretchModeTorsion.back().mvRotatedAtomList))
               {
                  #ifdef __DEBUG__
                  cout<<"Duplicate StretchModeTorsion ";mvStretchModeTorsion.back().Print(cout);cout<<endl;
                  #endif
                  mvStretchModeTorsion.pop_back();
               }
            }
         }
         if(  (mvStretchModeTorsion.back().mvRotatedAtomList.size()>((mvpAtom.size()+1)/2))
            ||(mvStretchModeTorsion.back().mvRotatedAtomList.size()==0))
         {
            #ifdef __DEBUG__
            cout<<"Rejecting StretchModeTorsion ";mvStretchModeTorsion.back().Print(cout);cout<<endl;
            #endif
            mvStretchModeTorsion.pop_back();
         }
         else if(mvStretchModeTorsion.back().mvRotatedAtomList.size()==((mvpAtom.size()+1)/2))
                  break;//we translate exactly half of the atoms, so skip the other half
      }
   }
   for(unsigned long i=0;i<mvpBond.size();++i)
   {//Single chains
      //if((mFlexModel.GetChoice()!=0)&&(false==mvpBond[i]->IsFreeTorsion())) continue;
      MolAtom* const atom1=&(mvpBond[i]->GetAtom1());
      MolAtom* const atom2=&(mvpBond[i]->GetAtom2());
      for(unsigned int j=1;j<=2;++j)
      {
         const set<MolAtom*> *pConn;
         if(j==1) pConn=&(mConnectivityTable[atom1]);
         else pConn=&(mConnectivityTable[atom2]);
         
         for(set<MolAtom*>::const_iterator pos=pConn->begin();pos!=pConn->end();++pos)
         {
            if((*pos==atom2)||(*pos==atom1)) continue;
            mvStretchModeTorsion.push_back(StretchModeTorsion(*atom1,*atom2,0));
            mvStretchModeTorsion.back().mvRotatedAtomList.insert(atom1);
            mvStretchModeTorsion.back().mvRotatedAtomList.insert(atom2);
            ExpandAtomGroupRecursive(*pos,mConnectivityTable,
                                     mvStretchModeTorsion.back().mvRotatedAtomList);
            mvStretchModeTorsion.back().mvRotatedAtomList.erase(atom1);
            mvStretchModeTorsion.back().mvRotatedAtomList.erase(atom2);
            
            for(vector<MolDihedralAngle*>::const_iterator dih=mvpDihedralAngle.begin();dih!=mvpDihedralAngle.end();++dih)
            {
               // :TODO: There are some other weird cases to take into account,
               // for restraints with atoms *not* connected to another
               // More generally, should check the list of atoms rotated.
               if(  ((&((*dih)->GetAtom2())==atom1) && (&((*dih)->GetAtom3())==atom2))
                  ||((&((*dih)->GetAtom3())==atom1) && (&((*dih)->GetAtom2())==atom2)))
               {
                  const unsigned long ct1=mvStretchModeTorsion.back().mvRotatedAtomList.count(&((*dih)->GetAtom1()));
                  const unsigned long ct4=mvStretchModeTorsion.back().mvRotatedAtomList.count(&((*dih)->GetAtom4()));

                  if((ct1+ct4)==1)// One of the atom is rotated, not the other
                  {
                     mvStretchModeTorsion.back().mpDihedralAngle=*dih;
                     //:TODO: Check sense of rotation !
                     if(ct4==1)
                     {
                        mvStretchModeTorsion.back().mpAtom1=&((*dih)->GetAtom2());
                        mvStretchModeTorsion.back().mpAtom2=&((*dih)->GetAtom3());
                     }
                     else
                     {
                        mvStretchModeTorsion.back().mpAtom1=&((*dih)->GetAtom3());
                        mvStretchModeTorsion.back().mpAtom2=&((*dih)->GetAtom2());
                     }
                  }
               }
            }
            
            if(mvStretchModeTorsion.size()>1)
            {//Duplicate ?
               // Does not work with a const_reverse_iterator ?
               // http://gcc.gnu.org/bugzilla/show_bug.cgi?id=11729
               list<StretchModeTorsion>::reverse_iterator mode=mvStretchModeTorsion.rbegin();
               ++mode;
               for(;mode!=mvStretchModeTorsion.rend();++mode)
               {
                  if( (  ((mode->mpAtom1==atom1)&&(mode->mpAtom2==atom2)) 
                       ||((mode->mpAtom1==atom2)&&(mode->mpAtom2==atom1)))
                     &&(mode->mvRotatedAtomList==mvStretchModeTorsion.back().mvRotatedAtomList))
                  {
                     #ifdef __DEBUG__
                     cout<<"Duplicate StretchModeTorsion ";mvStretchModeTorsion.back().Print(cout);cout<<endl;
                     #endif
                     mvStretchModeTorsion.pop_back();
                     break;
                  }
               }
            }
            if(  (mvStretchModeTorsion.back().mvRotatedAtomList.size()>((mvpAtom.size()+1)/2))
               ||(mvStretchModeTorsion.back().mvRotatedAtomList.size()==0))
            {
               #ifdef __DEBUG__
               cout<<"Rejecting StretchModeTorsion ";mvStretchModeTorsion.back().Print(cout);cout<<endl;
               #endif
               mvStretchModeTorsion.pop_back();
            }
            else if(mvStretchModeTorsion.back().mvRotatedAtomList.size()==((mvpAtom.size()+1)/2))
                     break;//we translate exactly half of the atoms, so skip the other half
         }
      }
   }
   // find rigid groups broken by each mode
   for(list<StretchModeTorsion>::iterator pos=mvStretchModeTorsion.begin();
       pos!=mvStretchModeTorsion.end();)
   {
      TAU_PROFILE_START(timer5);
      bool keep=true;
      for(vector<RigidGroup*>::const_iterator group=mvRigidGroup.begin();
          group!=mvRigidGroup.end();++group)
      {
         unsigned long ct=0;
         for(set<MolAtom *>::const_iterator at=(*group)->begin();at!=(*group)->end();++at)
            ct += pos->mvRotatedAtomList.count(*at);
         if(ct>0)
         {
            // Add the axis atoms, which do not move relatively to the rotated atoms
            ct += (*group)->count(pos->mpAtom1);
            ct += (*group)->count(pos->mpAtom2);
            if(ct!=(*group)->size())
            {
               keep=false;
               #ifdef __DEBUG__
               pos->Print(cout);
               cout<<"       Breaks Rigid Group:";
               for(set<MolAtom *>::const_iterator at=(*group)->begin();at!=(*group)->end();++at)
                  cout<<(*at)->GetName()<<" ";
               cout<<endl;
               #endif
               break;
            }
         }
      }
      if(keep) ++pos;
      else pos=mvStretchModeTorsion.erase(pos);
      TAU_PROFILE_STOP(timer5);
   }
   // Generate 5 completely random atomic positions
      this->SaveParamSet(mLocalParamSet);
      unsigned long paramSetRandom[5];
      for(unsigned long i=0;i<5;++i)
      {
         for(vector<MolAtom*>::iterator pos=mvpAtom.begin();pos!=mvpAtom.end();++pos)
         {
            (*pos)->SetX(100.*rand()/(REAL) RAND_MAX);
            (*pos)->SetY(100.*rand()/(REAL) RAND_MAX);
            (*pos)->SetZ(100.*rand()/(REAL) RAND_MAX);
         }
         paramSetRandom[i]=this->CreateParamSet();
      }
   // find bond lengths broken by each mode
   for(list<StretchModeTorsion>::iterator pos=mvStretchModeTorsion.begin();
       pos!=mvStretchModeTorsion.end();)
   {
      TAU_PROFILE_START(timer2);
      bool keep=true;
      pos->mvpBrokenBond.clear();
      for(vector<MolBond*>::const_iterator r=mvpBond.begin();r!=mvpBond.end();++r)
      {
         unsigned int ct=0;
         for(set<MolAtom *>::const_iterator at=pos->mvRotatedAtomList.begin();
             at!=pos->mvRotatedAtomList.end();++at)
         {
            if(*at==&((*r)->GetAtom1())) ct++;
            if(*at==&((*r)->GetAtom2())) ct++;
         }
         bool broken=true;
         // If we moved either both or non of the bond atom, the bond length is unchanged.
         if((ct==0)||(ct==2)) broken=false;
         if(broken)
         {// Make sure with derivatives
            REAL d=0;
            for(unsigned long i=0;i<5;++i)
            {
               this->RestoreParamSet(paramSetRandom[i]);
               pos->CalcDeriv(false);
               (*r)->GetLogLikelihood(true,true);
               d += abs((*r)->GetDeriv(pos->mDerivXYZ));
               if(d>0.01) break;
            }
            if(abs(d)<=0.01) broken=false;
         }
         if(broken) pos->mvpBrokenBond.insert(make_pair(*r,0));
      }
      if(mFlexModel.GetChoice()==2)
      {
         if(pos->mvpBrokenBond.size()>0) keep=false;
      }
      if(keep) ++pos;
      else pos=mvStretchModeTorsion.erase(pos);
      TAU_PROFILE_STOP(timer2);
   }
   // find bond angles broken by each mode
   for(list<StretchModeTorsion>::iterator pos=mvStretchModeTorsion.begin();
       pos!=mvStretchModeTorsion.end();)
   {
      TAU_PROFILE_START(timer3);
      bool keep=true;
      pos->mvpBrokenBondAngle.clear();
      for(vector<MolBondAngle*>::const_iterator r=mvpBondAngle.begin();r!=mvpBondAngle.end();++r)
      {
         unsigned int ct=0;
         for(set<MolAtom *>::const_iterator at=pos->mvRotatedAtomList.begin();
             at!=pos->mvRotatedAtomList.end();++at)
         {
            if(*at==&((*r)->GetAtom1())) ct++;
            if(*at==&((*r)->GetAtom2())) ct++;
            if(*at==&((*r)->GetAtom3())) ct++;
         }
         bool broken=true;
         if((ct==0)||(ct==3)) broken=false;
         if(broken)
         {// Make sure with derivatives
            REAL d=0;
            for(unsigned long i=0;i<5;++i)
            {
               this->RestoreParamSet(paramSetRandom[i]);
               pos->CalcDeriv(false);
               (*r)->GetLogLikelihood(true,true);
               d += abs((*r)->GetDeriv(pos->mDerivXYZ));
               if(d>0.01) break;
            }
            if(abs(d)<=0.01) broken=false;
         }
         if(broken) pos->mvpBrokenBondAngle.insert(make_pair(*r,0));
      }
      if(mFlexModel.GetChoice()==2)
      {
         if(pos->mvpBrokenBond.size()>0) keep=false;
      }
      if(keep) ++pos;
      else pos=mvStretchModeTorsion.erase(pos);
      TAU_PROFILE_STOP(timer3);
   }
   // find dihedral angles broken by each mode
   for(list<StretchModeTorsion>::iterator pos=mvStretchModeTorsion.begin();
       pos!=mvStretchModeTorsion.end();)
   {
      TAU_PROFILE_START(timer4);
      bool keep=true;
      pos->mvpBrokenDihedralAngle.clear();
      for(vector<MolDihedralAngle*>::const_iterator r=mvpDihedralAngle.begin();r!=mvpDihedralAngle.end();++r)
      {
         unsigned int ct=0;
         for(set<MolAtom *>::const_iterator at=pos->mvRotatedAtomList.begin();
             at!=pos->mvRotatedAtomList.end();++at)
         {
            if(*at==&((*r)->GetAtom1())) ct++;
            if(*at==&((*r)->GetAtom2())) ct++;
            if(*at==&((*r)->GetAtom3())) ct++;
            if(*at==&((*r)->GetAtom4())) ct++;
         }
         bool broken=true;
         if((ct==0)||(ct==4)) broken=false;
         if(broken)
         {// Make sure with derivatives
            REAL d=0;
            for(unsigned long i=0;i<5;++i)
            {
               this->RestoreParamSet(paramSetRandom[i]);
               pos->CalcDeriv(false);
               (*r)->GetLogLikelihood(true,true);
               d += abs((*r)->GetDeriv(pos->mDerivXYZ));
               if(d>0.01) break;
            }
            if(abs(d)<=0.01) broken=false;
         }
         if(broken) pos->mvpBrokenDihedralAngle.insert(make_pair(*r,0));
      }
      if(mFlexModel.GetChoice()==2)
      {
         int nb=pos->mvpBrokenDihedralAngle.size();
         if(pos->mpDihedralAngle!=0) nb -= 1;
         if(nb>0) keep=false;
      }
      if(keep) ++pos;
      else pos=mvStretchModeTorsion.erase(pos);
      TAU_PROFILE_STOP(timer4);
   }

   for(unsigned long i=0;i<5;++i) this->ClearParamSet(paramSetRandom[i]);
   this->RestoreParamSet(mLocalParamSet);
   
   #ifdef __DEBUG__
   cout<<"List of Dihedral Angle stretch modes("<<mvStretchModeTorsion.size()<<")"<<endl;
   for(list<StretchModeTorsion>::const_iterator pos=mvStretchModeTorsion.begin();
       pos!=mvStretchModeTorsion.end();++pos)
   {
      cout<<"   Dihedral Angle:"
          <<pos->mpAtom1->GetName()<<"-"
          <<pos->mpAtom2->GetName()<<", Rotated Atoms:  ";
      for(set<MolAtom*>::const_iterator atom=pos->mvRotatedAtomList.begin();
          atom!=pos->mvRotatedAtomList.end();++atom)
      {
         cout<<(*atom)->GetName()<<",";
      }
      if(pos->mpDihedralAngle!=0) 
         cout<<endl<< "      ->restrained by dihedral angle "<<pos->mpDihedralAngle->GetName()
             <<"to :"<<pos->mpDihedralAngle->GetAngle0()*RAD2DEG
             <<"�, sigma="<<pos->mpDihedralAngle->GetAngleSigma()*RAD2DEG
             <<"�, delta="<<pos->mpDihedralAngle->GetAngleDelta()*RAD2DEG<<"�";
      if(pos->mvpBrokenBond.size()>0)
      {
         cout<<endl<<"       Broken bonds:";
         for(map<const MolBond*,REAL>::const_iterator bond=pos->mvpBrokenBond.begin();
             bond!=pos->mvpBrokenBond.end();++bond)
            cout<<bond->first->GetName()<<", ";
      }
      if(pos->mvpBrokenBondAngle.size()>0)
      {
         cout<<endl<<"       Broken bond angles:";
         for(map<const MolBondAngle*,REAL>::const_iterator angle=pos->mvpBrokenBondAngle.begin();
             angle!=pos->mvpBrokenBondAngle.end();++angle)
            cout<<angle->first->GetName()<<", ";
      }
      if(pos->mvpBrokenDihedralAngle.size()>0)
      {
         cout<<endl<<"       Broken dihedral angles:";
         for(map<const MolDihedralAngle*,REAL>::const_iterator 
             angle=pos->mvpBrokenDihedralAngle.begin();
             angle!=pos->mvpBrokenDihedralAngle.end();++angle)
            cout<<angle->first->GetName()<<", ";
      }
      cout<<endl;
   }
   #endif
   mClockStretchModeTorsion.Click();
   VFN_DEBUG_EXIT("Molecule::BuildStretchModeTorsion()",7)
}

void Molecule::BuildStretchModeTwist()
{
   #if 0
   if(  (mClockStretchModeTwist>mClockBondList)
      &&(mClockStretchModeTwist>mClockAtomList)
      &&(mClockStretchModeTwist>mClockBondAngleList)
      &&(mClockStretchModeTwist>mClockDihedralAngleList)) return;
   #endif
   VFN_DEBUG_ENTRY("Molecule::BuildStretchModeTwist()",7)
   this->BuildConnectivityTable();
   mvStretchModeTwist.clear();
   
   // For each pair of atoms, build an internal chain to twist.
   for(vector<MolAtom*>::const_iterator atom1=this->GetAtomList().begin();
       atom1!=this->GetAtomList().end();++atom1)
   {
      const set<MolAtom*> *pConn=&(mConnectivityTable[*atom1]);
      vector<MolAtom*>::const_iterator atom2=atom1;
      atom2++;
      for(;atom2!=this->GetAtomList().end();++atom2)
      {
         for(set<MolAtom*>::const_iterator pos=pConn->begin();pos!=pConn->end();++pos)
         {// Start from one atom connected to atom1
            if(*pos==*atom2) continue;
            mvStretchModeTwist.push_back(StretchModeTwist(**atom1,**atom2));
            mvStretchModeTwist.back().mvRotatedAtomList.insert(*atom1);
            ExpandAtomGroupRecursive(*pos,mConnectivityTable,
                                     mvStretchModeTwist.back().mvRotatedAtomList,
                                     *atom2);
            //Check if this chains actually leads to atom2
            set<MolAtom*>::const_iterator check
                  =find(mvStretchModeTwist.back().mvRotatedAtomList.begin(),
                        mvStretchModeTwist.back().mvRotatedAtomList.end(),*atom2);
            bool keep =true;
            if(  (check==mvStretchModeTwist.back().mvRotatedAtomList.end())
               ||(mvStretchModeTwist.back().mvRotatedAtomList.size()<3)
               ||(mvStretchModeTwist.back().mvRotatedAtomList.size()>=((mvpAtom.size()+1)/2)))
            {
               keep=false;
            }
            if(keep)
            {
               mvStretchModeTwist.back().mvRotatedAtomList.erase(*atom1);
               mvStretchModeTwist.back().mvRotatedAtomList.erase(*atom2);
               if(  (mvStretchModeTwist.back().mvRotatedAtomList.size()>=(mvpAtom.size()/2))
                  ||(mvStretchModeTwist.back().mvRotatedAtomList.size()==0))
               {
                  #ifdef __DEBUG__
                  cout<<"Rejecting StretchModeTwist ";mvStretchModeTwist.back().Print(cout);cout<<endl;
                  #endif
                  keep=false;
               }
            }
            if(keep)
            {
               if(mvStretchModeTwist.size()>1)
               {//Duplicate ?
                  // Does not work with a const_reverse_iterator ?
                  // http://gcc.gnu.org/bugzilla/show_bug.cgi?id=11729
                  list<StretchModeTwist>::reverse_iterator mode=mvStretchModeTwist.rbegin();
                  ++mode;
                  for(;mode!=mvStretchModeTwist.rend();++mode)
                  {
                     if( (  ((mode->mpAtom1==*atom1)&&(mode->mpAtom2==*atom2)) 
                          ||((mode->mpAtom1==*atom2)&&(mode->mpAtom2==*atom1)))
                        &&(mode->mvRotatedAtomList==mvStretchModeTwist.back().mvRotatedAtomList))
                     {
                        #ifdef __DEBUG__
                        cout<<"Duplicate StretchModeTwist ";mvStretchModeTwist.back().Print(cout);cout<<endl;
                        #endif
                        keep=false;
                     }
                     if(!keep) break;
                  }
               }
            }
            if(keep)
            {
               if(mvStretchModeTorsion.size()>0) // Torsion and twist modes can be identical for cycles.
               {//Duplicate ?
                  // Does not work with a const_reverse_iterator ?
                  // http://gcc.gnu.org/bugzilla/show_bug.cgi?id=11729
                  list<StretchModeTorsion>::reverse_iterator mode;
                  for(mode=mvStretchModeTorsion.rbegin();mode!=mvStretchModeTorsion.rend();++mode)
                  {
                     if( (  ((mode->mpAtom1==*atom1)&&(mode->mpAtom2==*atom2)) 
                          ||((mode->mpAtom1==*atom2)&&(mode->mpAtom2==*atom1)))
                        &&(mode->mvRotatedAtomList==mvStretchModeTwist.back().mvRotatedAtomList))
                     {
                        #ifdef __DEBUG__
                        cout<<"Duplicate StretchModeTwist (with Torsion) ";mvStretchModeTwist.back().Print(cout);cout<<endl;
                        #endif
                        keep=false;
                     }
                     if(!keep) break;
                  }
               }
            }
            if(!keep) mvStretchModeTwist.pop_back();
         }
      }

   }
   // Generate 5 completely random atomic positions
      this->SaveParamSet(mLocalParamSet);
      unsigned long paramSetRandom[5];
      for(unsigned long i=0;i<5;++i)
      {
         for(vector<MolAtom*>::iterator pos=mvpAtom.begin();pos!=mvpAtom.end();++pos)
         {
            (*pos)->SetX(100.*rand()/(REAL) RAND_MAX);
            (*pos)->SetY(100.*rand()/(REAL) RAND_MAX);
            (*pos)->SetZ(100.*rand()/(REAL) RAND_MAX);
         }
         paramSetRandom[i]=this->CreateParamSet();
      }
   // find bond, bond angles and dihedral angles broken by each mode
   for(list<StretchModeTwist>::iterator pos=mvStretchModeTwist.begin();
       pos!=mvStretchModeTwist.end();)
   {
      pos->mvpBrokenBond.clear();
      pos->mvpBrokenBondAngle.clear();
      pos->mvpBrokenDihedralAngle.clear();
      pos->CalcDeriv();
      #ifdef __DEBUG__
      cout<<"   DerivLLK for Twist mode around:"
          <<pos->mpAtom1->GetName()<<"-"
          <<pos->mpAtom2->GetName()<<": Moving atoms:";
      for(set<MolAtom*>::const_iterator atom=pos->mvRotatedAtomList.begin();
          atom!=pos->mvRotatedAtomList.end();++atom)
         cout<<(*atom)->GetName()<<",";
      cout<<endl;
      #endif
      for(vector<MolBond*>::const_iterator r=mvpBond.begin();r!=mvpBond.end();++r)
      {
         (*r)->GetLogLikelihood(true,true);
         REAL d=0;
         for(unsigned long i=0;i<5;++i)
         {
            this->RestoreParamSet(paramSetRandom[i]);
            d += abs((*r)->GetDeriv(pos->mDerivXYZ));
         }
         if(abs(d)>0.1)
         {
            #ifdef __DEBUG__
            cout<<"       Bond "<<(*r)->GetName()
                <<": dLength/da="<<d<<endl;
            #endif
            pos->mvpBrokenBond.insert(make_pair(*r,0.0));
         }
      }
      for(vector<MolBondAngle*>::const_iterator r=mvpBondAngle.begin();r!=mvpBondAngle.end();++r)
      {
         (*r)->GetLogLikelihood(true,true);
         REAL d=0;
         for(unsigned long i=0;i<5;++i)
         {
            this->RestoreParamSet(paramSetRandom[i]);
            d += abs((*r)->GetDeriv(pos->mDerivXYZ));
         }
         if(abs(d)>0.01)
         {
            #ifdef __DEBUG__
            cout<<"       BondAngle:"<<(*r)->GetName()<<": dAngle/da="<<d<<endl;
            #endif
            pos->mvpBrokenBondAngle.insert(make_pair(*r,0.0));
         }
      }
      for(vector<MolDihedralAngle*>::const_iterator r=mvpDihedralAngle.begin();
          r!=mvpDihedralAngle.end();++r)
      {
         //if(*r==pos->mpDihedralAngle) continue;
         (*r)->GetLogLikelihood(true,true);
         REAL d=0;
         for(unsigned long i=0;i<5;++i)
         {
            this->RestoreParamSet(paramSetRandom[i]);
            d += abs((*r)->GetDeriv(pos->mDerivXYZ));
         }
         if(abs(d)>0.01)
         {
            #ifdef __DEBUG__
            cout<<"       DihedralAngle:"<<(*r)->GetName()<<": dAngle/da="<<d<<endl;
            #endif
            pos->mvpBrokenDihedralAngle.insert(make_pair(*r,0.0));
         }
      }
      //Get rid of stretch modes that break rigid groups
      bool keep=true;
      for(vector<RigidGroup*>::const_iterator group=mvRigidGroup.begin();
          group!=mvRigidGroup.end();++group)
      {
         unsigned long ct=0;
         for(set<MolAtom *>::const_iterator at=(*group)->begin();at!=(*group)->end();++at)
            ct += pos->mvRotatedAtomList.count(*at);
         if(ct>0)
         {
            // Add atom1 and atom2 to the count only if they are in the group
            // They do not move in absolute or relatively to the group.
            ct += (*group)->count(pos->mpAtom1);
            ct += (*group)->count(pos->mpAtom2);
            if(ct!=(*group)->size())
            {
               keep=false;
               #ifdef __DEBUG__
               cout<<"       Breaks Rigid Group:"<<ct<<"!="<<(*group)->size()<<":";
               for(set<MolAtom *>::const_iterator at=(*group)->begin();at!=(*group)->end();++at)
                  cout<<(*at)->GetName()<<" ";
               cout<<endl;
               #endif
               break;
            }
         }
      }
      if(mFlexModel.GetChoice()==2)
      {
         if(pos->mvpBrokenBond.size()+pos->mvpBrokenBondAngle.size()+pos->mvpBrokenDihedralAngle.size()) keep=false;
      }
      if(keep) ++pos;
      else pos=mvStretchModeTwist.erase(pos);
   }
   for(unsigned long i=0;i<5;++i) this->ClearParamSet(paramSetRandom[i]);
   this->RestoreParamSet(mLocalParamSet);
   
   #ifdef __DEBUG__
   cout<<"List of Twist stretch modes("<<mvStretchModeTwist.size()<<")"<<endl;
   for(list<StretchModeTwist>::const_iterator pos=mvStretchModeTwist.begin();
       pos!=mvStretchModeTwist.end();++pos)
   {
      cout<<"   Twist mode:"
          <<pos->mpAtom1->GetName()<<"-"
          <<pos->mpAtom2->GetName()<<", Rotated Atoms:  ";
      for(set<MolAtom*>::const_iterator atom=pos->mvRotatedAtomList.begin();
          atom!=pos->mvRotatedAtomList.end();++atom)
      {
         cout<<(*atom)->GetName()<<",";
      }
      if(pos->mvpBrokenBond.size()>0)
      {
         cout<<endl<<"       Broken bonds:";
         for(map<const MolBond*,REAL>::const_iterator bond=pos->mvpBrokenBond.begin();
             bond!=pos->mvpBrokenBond.end();++bond)
            cout<<bond->first->GetName()<<", ";
      }
      if(pos->mvpBrokenBondAngle.size()>0)
      {
         cout<<endl<<"       Broken bond angles:";
         for(map<const MolBondAngle*,REAL>::const_iterator angle=pos->mvpBrokenBondAngle.begin();
             angle!=pos->mvpBrokenBondAngle.end();++angle)
            cout<<angle->first->GetName()<<", ";
      }
      if(pos->mvpBrokenDihedralAngle.size()>0)
      {
         cout<<endl<<"       Broken dihedral angles:";
         for(map<const MolDihedralAngle*,REAL>::const_iterator 
             angle=pos->mvpBrokenDihedralAngle.begin();
             angle!=pos->mvpBrokenDihedralAngle.end();++angle)
            cout<<angle->first->GetName()<<", ";
      }
      cout<<endl;
   }
   #endif
   mClockStretchModeTwist.Click();
   VFN_DEBUG_EXIT("Molecule::BuildStretchModeTwist()",7)
}

void Molecule::TuneGlobalOptimRotationAmplitude()
{
   VFN_DEBUG_ENTRY("Molecule::TuneGlobalOptimRotationAmplitude()",5)
   unsigned long initialConfig=this->CreateParamSet("Initial Configuration");
   const unsigned int nbTest=100;

   // First we store in mBaseRotationAmplitude the cumulated atomic displacement 
   // for nbTest rotations of 0.01 rad
   for(list<StretchModeBondAngle>::iterator pos=mvStretchModeBondAngle.begin();
       pos!=mvStretchModeBondAngle.end();++pos) pos->mBaseAmplitude=0;
   for(list<StretchModeTorsion>::iterator pos=mvStretchModeTorsion.begin();
       pos!=mvStretchModeTorsion.end();++pos) pos->mBaseAmplitude=0;

   REAL displacement=0;//For the global Molecule rotation

   for(unsigned int j=0;j<nbTest;j++)
   {
      this->RandomizeConfiguration();
      // Atomic positions, orthonormal coordinates
      vector<REAL> x0(this->GetNbComponent());
      vector<REAL> y0(this->GetNbComponent());
      vector<REAL> z0(this->GetNbComponent());
      // Center of molecule coords
      REAL xc=0.,yc=0.,zc=0.;
      for(long i=0;i<this->GetNbComponent();++i)
      {
         x0[i]=mvpAtom[i]->GetX(); xc += x0[i];
         y0[i]=mvpAtom[i]->GetY(); yc += y0[i];
         z0[i]=mvpAtom[i]->GetZ(); zc += z0[i];
      }
      xc /= (REAL)(this->GetNbComponent());
      yc /= (REAL)(this->GetNbComponent());
      zc /= (REAL)(this->GetNbComponent());
      // Record displacement amplitude for torsion angles
      REAL dx,dy,dz;
      for(list<StretchModeBondAngle>::iterator pos=mvStretchModeBondAngle.begin();
          pos!=mvStretchModeBondAngle.end();++pos)
      {
         const REAL dx10=pos->mpAtom0->GetX()-pos->mpAtom1->GetX();
         const REAL dy10=pos->mpAtom0->GetY()-pos->mpAtom1->GetY();
         const REAL dz10=pos->mpAtom0->GetZ()-pos->mpAtom1->GetZ();
         const REAL dx12=pos->mpAtom2->GetX()-pos->mpAtom1->GetX();
         const REAL dy12=pos->mpAtom2->GetY()-pos->mpAtom1->GetY();
         const REAL dz12=pos->mpAtom2->GetZ()-pos->mpAtom1->GetZ();

         const REAL vx=dy10*dz12-dz10*dy12;
         const REAL vy=dz10*dx12-dx10*dz12;
         const REAL vz=dx10*dy12-dy10*dx12;
         this->RotateAtomGroup(*(pos->mpAtom1),vx,vy,vz,pos->mvRotatedAtomList,0.01,false);
         for(long i=0;i<this->GetNbComponent();++i)
         {
            dx=x0[i]-mvpAtom[i]->GetX();
            dy=y0[i]-mvpAtom[i]->GetY();
            dz=z0[i]-mvpAtom[i]->GetZ();
            pos->mBaseAmplitude+=sqrt(abs(dx*dx+dy*dy+dz*dz));
         }
         this->RotateAtomGroup(*(pos->mpAtom1),vx,vy,vz,pos->mvRotatedAtomList,-0.01,false);
      }
      for(list<StretchModeTorsion>::iterator pos=mvStretchModeTorsion.begin();
          pos!=mvStretchModeTorsion.end();++pos)
      {
         this->RotateAtomGroup(*(pos->mpAtom1),*(pos->mpAtom2),pos->mvRotatedAtomList,0.01,false);
         for(long i=0;i<this->GetNbComponent();++i)
         {
            dx=x0[i]-mvpAtom[i]->GetX();
            dy=y0[i]-mvpAtom[i]->GetY();
            dz=z0[i]-mvpAtom[i]->GetZ();
            pos->mBaseAmplitude+=sqrt(abs(dx*dx+dy*dy+dz*dz));
         }
         this->RotateAtomGroup(*(pos->mpAtom1),*(pos->mpAtom2),pos->mvRotatedAtomList,-0.01,false);
      }
      // Record displacement amplitude for global rotation, for 10 random rot axis
      for(unsigned int k=0;k<10;++k)
      {
         Quaternion quat=Quaternion::RotationQuaternion
                     (mBaseRotationAmplitude,(REAL)rand(),(REAL)rand(),(REAL)rand());
         for(long i=0;i<this->GetNbComponent();++i)
         {
            REAL x=x0[i]-xc;
            REAL y=y0[i]-yc;
            REAL z=z0[i]-zc;
            quat.RotateVector(x,y,z);
            dx=(x0[i]-xc)-x;
            dy=(y0[i]-yc)-y;
            dz=(z0[i]-zc)-z;
            displacement+=sqrt(abs(dx*dx+dy*dy+dz*dz));
         }
      }
   }
   // Modify base rotation amplitudes, for an average 0.02 Angstroem displacement
   for(list<StretchModeBondAngle>::iterator pos=mvStretchModeBondAngle.begin();
       pos!=mvStretchModeBondAngle.end();++pos)
   {
      pos->mBaseAmplitude/=(REAL)(nbTest*this->GetNbComponent());//(REAL)(nbTest*(0.5*pos->mvRotatedAtomList.size()+0.5*this->GetNbComponent()));
      pos->mBaseAmplitude=0.1*0.01/pos->mBaseAmplitude;
      if(pos->mBaseAmplitude>.17) pos->mBaseAmplitude=.17;
      bool free=true;
      if((pos->mvpBrokenBond.size()+pos->mvpBrokenBondAngle.size()+pos->mvpBrokenDihedralAngle.size())>0)
      {
         pos->mBaseAmplitude/=10;
         free=false;
      }
      #if 1// def __DEBUG__
      cout<<"ANGLE :"
          <<pos->mpAtom0->GetName()<<"-"
          <<pos->mpAtom1->GetName()<<"-"
          <<pos->mpAtom2->GetName()<<":";
      for(set<MolAtom*>::const_iterator atom=pos->mvRotatedAtomList.begin();
          atom!=pos->mvRotatedAtomList.end();++atom) cout<<(*atom)->GetName()<<",";
      cout <<": base rotation="<<pos->mBaseAmplitude*RAD2DEG<<" free="<<free<<endl;
      #endif
   }
   // Modify base rotation amplitudes, for an average 0.1 Angstroem displacement
   for(list<StretchModeTorsion>::iterator pos=mvStretchModeTorsion.begin();
       pos!=mvStretchModeTorsion.end();++pos)
   {
      pos->mBaseAmplitude/=(REAL)(nbTest*this->GetNbComponent());//(REAL)(nbTest*(0.5*pos->mvRotatedAtomList.size()+0.5*this->GetNbComponent()));
      pos->mBaseAmplitude=0.1*0.01/pos->mBaseAmplitude;
      if(pos->mBaseAmplitude>.17) pos->mBaseAmplitude=.17;
      bool free=true;
      if((pos->mvpBrokenBond.size()+pos->mvpBrokenBondAngle.size()+pos->mvpBrokenDihedralAngle.size())>0)
      {
         pos->mBaseAmplitude/=10;
         free=false;
      }
      #if 1// def __DEBUG__
      cout<<"TORSION :"
          <<pos->mpAtom1->GetName()<<"-"
          <<pos->mpAtom2->GetName()<<":";
      for(set<MolAtom*>::const_iterator atom=pos->mvRotatedAtomList.begin();
          atom!=pos->mvRotatedAtomList.end();++atom) cout<<(*atom)->GetName()<<",";
      cout <<": base rotation="<<pos->mBaseAmplitude*RAD2DEG<<" free="<<free<<endl;
      #endif
   }
   // Modify base rotation amplitudes of twist modes, for an average 0.1 Angstroem displacement
   for(list<StretchModeTwist>::iterator pos=mvStretchModeTwist.begin();
       pos!=mvStretchModeTwist.end();++pos)
   {
      pos->mBaseAmplitude/=(REAL)(nbTest*this->GetNbComponent());//(REAL)(nbTest*(0.5*pos->mvRotatedAtomList.size()+0.5*this->GetNbComponent()));
      pos->mBaseAmplitude=0.1*0.01/pos->mBaseAmplitude;
      if(pos->mBaseAmplitude>.17) pos->mBaseAmplitude=.17;
      bool free=true;
      if((pos->mvpBrokenBond.size()+pos->mvpBrokenBondAngle.size()+pos->mvpBrokenDihedralAngle.size())>0)
      {
         pos->mBaseAmplitude/=10;
         free=false;
      }
      #if 1// def __DEBUG__
      cout<<"TWIST :"
          <<pos->mpAtom1->GetName()<<"-"
          <<pos->mpAtom2->GetName()<<":";
      for(set<MolAtom*>::const_iterator atom=pos->mvRotatedAtomList.begin();
          atom!=pos->mvRotatedAtomList.end();++atom) cout<<(*atom)->GetName()<<",";
      cout <<": base rotation="<<pos->mBaseAmplitude*RAD2DEG<<" free="<<free<<endl;
      #endif
   }
   // Same for global rotation
      displacement/=(REAL)(10*nbTest*this->GetNbComponent());
      //cout<<"Overall Atomic Displacement for Global Rotation:<d>="<<displacement;
      if(displacement>0) mBaseRotationAmplitude*=0.1/displacement;
      if(mBaseRotationAmplitude<(0.02*M_PI/20.))
      {
         mBaseRotationAmplitude=0.02*M_PI/20.;
         //cout <<"WARNING - too low Global BaseRotationAmplitude - setting to: "
         //     << mBaseAmplitude*RAD2DEG<< " �"<<endl;
      }
      if(mBaseRotationAmplitude>(0.02*M_PI*20.))
      {
         mBaseRotationAmplitude=0.02*M_PI*20.;
         //cout <<"WARNING - too high Global BaseRotationAmplitude - setting to: "
         //     << mBaseAmplitude*RAD2DEG<< " �"<<endl;
      }
      //cout <<" -> Base rotation="<<mBaseAmplitude*RAD2DEG<<"�"<<endl;

   // Move back atoms to initial position
   this->RestoreParamSet(initialConfig);
   VFN_DEBUG_EXIT("Molecule::TuneGlobalOptimRotationAmplitude()",5)
}

void Molecule::BuildFlipGroup()
{
   this->BuildConnectivityTable();
   if(  (mClockFlipGroup>mClockConnectivityTable)
      &&(mClockFlipGroup>mClockRigidGroup)) return;
   VFN_DEBUG_ENTRY("Molecule::BuildFlipGroup()",5)
   TAU_PROFILE("Molecule::BuildFlipGroup()","void ()",TAU_DEFAULT);
   mvFlipGroup.clear();
   
   for(vector<MolAtom*>::const_iterator atom0=this->GetAtomList().begin();
       atom0!=this->GetAtomList().end();++atom0)
   {
      const set<MolAtom*> *pConn=&(mConnectivityTable[*atom0]);
      if(pConn->size()<3) continue;
      // Build all chains
      for(set<MolAtom*>::const_iterator pos1=pConn->begin();pos1!=pConn->end();++pos1)
      {
         for(set<MolAtom*>::const_iterator pos2=pos1;pos2!=pConn->end();++pos2)
         {
            if(*pos2==*pos1) continue;
            if(mFlexModel.GetChoice()==0)
            {
               mvFlipGroup.push_back(FlipGroup(**atom0,**pos1,**pos2));
               bool foundRing=false;
               for(set<MolAtom*>::const_iterator pos=pConn->begin();pos!=pConn->end();++pos)
               {
                  if((pos==pos1)||(pos==pos2)) continue;
                  mvFlipGroup.back().mvRotatedChainList.push_back(
                     make_pair(*pos,set<MolAtom*>()));
                  mvFlipGroup.back().mvRotatedChainList.back().second.insert(*atom0);
                  ExpandAtomGroupRecursive(*pos,mConnectivityTable,
                                           mvFlipGroup.back().mvRotatedChainList.back().second);
                  mvFlipGroup.back().mvRotatedChainList.back().second.erase(*atom0);
                  set<MolAtom*>::const_iterator ringdetect1,ringdetect2;
                  ringdetect1=find(mvFlipGroup.back().mvRotatedChainList.back().second.begin(),
                                   mvFlipGroup.back().mvRotatedChainList.back().second.end(),
                                   *pos1);
                  ringdetect2=find(mvFlipGroup.back().mvRotatedChainList.back().second.begin(),
                                   mvFlipGroup.back().mvRotatedChainList.back().second.end(),
                                   *pos2);
                  if(  (ringdetect1!=mvFlipGroup.back().mvRotatedChainList.back().second.end())
                     ||(ringdetect2!=mvFlipGroup.back().mvRotatedChainList.back().second.end()))
                     foundRing=true;
               }
               unsigned long flipSize=0;
               for(list<pair<const MolAtom *,set<MolAtom*> > >::const_iterator 
                   chain=mvFlipGroup.back().mvRotatedChainList.begin();
                   chain!=mvFlipGroup.back().mvRotatedChainList.end();++chain)
                   flipSize+=chain->second.size();

               if(((flipSize*2)>mvpAtom.size())||foundRing) mvFlipGroup.pop_back();
            }
            // Add the entry which will exchange atom1 and atom2 (this entry can be a ring)
            mvFlipGroup.push_back(FlipGroup(**atom0,**pos1,**pos2));
            mvFlipGroup.back().mvRotatedChainList.push_back(
               make_pair(*atom0,set<MolAtom*>()));
               mvFlipGroup.back().mvRotatedChainList.back().second.insert(*atom0);
            ExpandAtomGroupRecursive(*pos1,mConnectivityTable,
                                     mvFlipGroup.back().mvRotatedChainList.back().second);
            ExpandAtomGroupRecursive(*pos2,mConnectivityTable,
                                     mvFlipGroup.back().mvRotatedChainList.back().second);
            mvFlipGroup.back().mvRotatedChainList.back().second.erase(*atom0);
            if((mvFlipGroup.back().mvRotatedChainList.back().second.size()*2)>mvpAtom.size())
               mvFlipGroup.pop_back();
         }
      }
   }
   // Exclude flip groups that include only a part of any rigid group
   for(list<FlipGroup>::iterator pos=mvFlipGroup.begin(); pos!=mvFlipGroup.end();)
   {
      // This should not be necessary ? Why keep one list for each chain, and not one big ?
      set<MolAtom*> fullset;
      for(list<pair<const MolAtom *,set<MolAtom *> > >::iterator chain=pos->mvRotatedChainList.begin();
            chain!=pos->mvRotatedChainList.end();++chain)
         for(set<MolAtom *>::const_iterator at=chain->second.begin();at!=chain->second.end();++at)
            fullset.insert(*at);
      
      bool keep=true;
      for(vector<RigidGroup*>::const_iterator group=mvRigidGroup.begin(); group!=mvRigidGroup.end();++group)
      {
         unsigned long ct=0;
         for(set<MolAtom *>::const_iterator at=fullset.begin();at!=fullset.end();++at)
            ct+=(*group)->count(*at);
         
         if((ct>0)&&(ct<(*group)->size())) {keep=false; break;}
      }
      
      if(!keep)
      {
         cout <<"EXCLUDING flip group (breaking a rigid group): "
              <<pos->mpAtom0->GetName()<<",exchanging bonds with "
              <<pos->mpAtom1->GetName()<<" and "
              <<pos->mpAtom2->GetName()<<", resulting in a 180deg rotation of atoms : ";
         for(set<MolAtom*>::iterator pos1=pos->mvRotatedChainList.begin()->second.begin();
             pos1!=pos->mvRotatedChainList.begin()->second.end();++pos1)
            cout<<(*pos1)->GetName()<<"  ";

         pos=mvFlipGroup.erase(pos);
      }
      else pos++;
   }
   // List them
   this->SaveParamSet(mLocalParamSet);
   #if 1//def __DEBUG__
   const REAL llk0=this->GetLogLikelihood();
   for(list<FlipGroup>::iterator pos=mvFlipGroup.begin();
       pos!=mvFlipGroup.end();++pos)
   {
      if(pos->mvRotatedChainList.begin()->first==pos->mpAtom0)
      {
         cout <<"Flip group from atom "
              <<pos->mpAtom0->GetName()<<",exchanging bonds with "
              <<pos->mpAtom1->GetName()<<" and "
              <<pos->mpAtom2->GetName()<<", resulting in a 180� rotation of atoms : ";
         for(set<MolAtom*>::iterator pos1=pos->mvRotatedChainList.begin()->second.begin();
             pos1!=pos->mvRotatedChainList.begin()->second.end();++pos1)
            cout<<(*pos1)->GetName()<<"  ";
      }
      else
      {
         cout <<"Flip group with respect to: "
              <<pos->mpAtom1->GetName()<<"-"
              <<pos->mpAtom0->GetName()<<"-"
              <<pos->mpAtom2->GetName()<<" : ";
         for(list<pair<const MolAtom *,set<MolAtom*> > >::const_iterator 
             chain=pos->mvRotatedChainList.begin();
             chain!=pos->mvRotatedChainList.end();++chain)
         {
            cout<<"    -"<<chain->first->GetName()<<":";
            for(set<MolAtom*>::const_iterator pos1=chain->second.begin();
                pos1!=chain->second.end();++pos1)
               cout<<(*pos1)->GetName()<<"  ";
         }
      }
      #if 0
      // test if they do not break something (dihedral angle restraint) ?
      // We seldom try flippping, so don't test - test is done during optimization
      this->FlipAtomGroup(*pos);
      const REAL dllk=this->GetLogLikelihood()-llk0;
      if(dllk>1000.)
      {
         pos = mvFlipGroup.erase(pos);
         --pos;
         cout <<" -> NOT a free flip, d(llk)="<<dllk;
         this->RestraintStatus(cout);
      }
      else cout <<" -> free flip, d(llk)="<<dllk;
      this->RestoreParamSet(mLocalParamSet);
      #endif
      cout<<endl;
   }
   #endif
   mClockFlipGroup.Click();
   VFN_DEBUG_EXIT("Molecule::BuildFlipGroup()",5)
}

void Molecule::BuildStretchModeGroups()
{
   // Assume Stretch modes have already been built
   map<const MolBond*,         set<const StretchMode*> > vpBond;
   for(list<StretchModeBondLength>::const_iterator mode=mvStretchModeBondLength.begin();
       mode!=mvStretchModeBondLength.end();++mode)
   {
      if(mode->mpBond!=0) vpBond[mode->mpBond].insert(&(*mode));
      for(map<const MolBond*,REAL>::const_iterator pos=mode->mvpBrokenBond.begin();
          pos!=mode->mvpBrokenBond.end();++pos)
          vpBond[pos->first].insert(&(*mode));
   }
   for(list<StretchModeBondAngle>::const_iterator mode=mvStretchModeBondAngle.begin();
       mode!=mvStretchModeBondAngle.end();++mode)
      for(map<const MolBond*,REAL>::const_iterator pos=mode->mvpBrokenBond.begin();
          pos!=mode->mvpBrokenBond.end();++pos)
          vpBond[pos->first].insert(&(*mode));

   for(list<StretchModeTorsion>::const_iterator mode=mvStretchModeTorsion.begin();
       mode!=mvStretchModeTorsion.end();++mode)
      for(map<const MolBond*,REAL>::const_iterator pos=mode->mvpBrokenBond.begin();
          pos!=mode->mvpBrokenBond.end();++pos)
          vpBond[pos->first].insert(&(*mode));
   for(list<StretchModeTwist>::const_iterator mode=mvStretchModeTwist.begin();
       mode!=mvStretchModeTwist.end();++mode)
      for(map<const MolBond*,REAL>::const_iterator pos=mode->mvpBrokenBond.begin();
          pos!=mode->mvpBrokenBond.end();++pos)
          vpBond[pos->first].insert(&(*mode));



   map<const MolBondAngle*,    set<const StretchMode*> > vpAngle;
   for(list<StretchModeBondLength>::const_iterator mode=mvStretchModeBondLength.begin();
       mode!=mvStretchModeBondLength.end();++mode)
      for(map<const MolBondAngle*,REAL>::const_iterator pos=mode->mvpBrokenBondAngle.begin();
          pos!=mode->mvpBrokenBondAngle.end();++pos)
          vpAngle[pos->first].insert(&(*mode));

   for(list<StretchModeBondAngle>::const_iterator mode=mvStretchModeBondAngle.begin();
       mode!=mvStretchModeBondAngle.end();++mode)
   {
      if(mode->mpBondAngle!=0) vpAngle[mode->mpBondAngle].insert(&(*mode));
      for(map<const MolBondAngle*,REAL>::const_iterator pos=mode->mvpBrokenBondAngle.begin();
          pos!=mode->mvpBrokenBondAngle.end();++pos)
          vpAngle[pos->first].insert(&(*mode));
   }
   for(list<StretchModeTorsion>::const_iterator mode=mvStretchModeTorsion.begin();
       mode!=mvStretchModeTorsion.end();++mode)
      for(map<const MolBondAngle*,REAL>::const_iterator pos=mode->mvpBrokenBondAngle.begin();
          pos!=mode->mvpBrokenBondAngle.end();++pos)
          vpAngle[pos->first].insert(&(*mode));
   for(list<StretchModeTwist>::const_iterator mode=mvStretchModeTwist.begin();
       mode!=mvStretchModeTwist.end();++mode)
      for(map<const MolBondAngle*,REAL>::const_iterator pos=mode->mvpBrokenBondAngle.begin();
          pos!=mode->mvpBrokenBondAngle.end();++pos)
          vpAngle[pos->first].insert(&(*mode));

   
   
   map<const MolDihedralAngle*,set<const StretchMode*> > vpDihed;
   for(list<StretchModeBondLength>::const_iterator mode=mvStretchModeBondLength.begin();
       mode!=mvStretchModeBondLength.end();++mode)
      for(map<const MolDihedralAngle*,REAL>::const_iterator pos=mode->mvpBrokenDihedralAngle.begin();
          pos!=mode->mvpBrokenDihedralAngle.end();++pos)
          vpDihed[pos->first].insert(&(*mode));

   for(list<StretchModeBondAngle>::const_iterator mode=mvStretchModeBondAngle.begin();
       mode!=mvStretchModeBondAngle.end();++mode)
      for(map<const MolDihedralAngle*,REAL>::const_iterator pos=mode->mvpBrokenDihedralAngle.begin();
          pos!=mode->mvpBrokenDihedralAngle.end();++pos)
          vpDihed[pos->first].insert(&(*mode));

   for(list<StretchModeTorsion>::const_iterator mode=mvStretchModeTorsion.begin();
       mode!=mvStretchModeTorsion.end();++mode)
   {
      if(mode->mpDihedralAngle!=0) vpDihed[mode->mpDihedralAngle].insert(&(*mode));
      for(map<const MolDihedralAngle*,REAL>::const_iterator pos=mode->mvpBrokenDihedralAngle.begin();
          pos!=mode->mvpBrokenDihedralAngle.end();++pos)
          vpDihed[pos->first].insert(&(*mode));
   }
   for(list<StretchModeTwist>::const_iterator mode=mvStretchModeTwist.begin();
       mode!=mvStretchModeTwist.end();++mode)
      for(map<const MolDihedralAngle*,REAL>::const_iterator pos=mode->mvpBrokenDihedralAngle.begin();
          pos!=mode->mvpBrokenDihedralAngle.end();++pos)
          vpDihed[pos->first].insert(&(*mode));
   #if 0
   for(map<const MolBond*,set<const StretchMode*> >::const_iterator pos=vpBond.begin();pos!=vpBond.end();++pos)
   {
      cout<<"Bond "<<pos->first->GetName()<<" is modified by the stretch modes:"<<endl;
      for(set<const StretchMode*>::const_iterator mode=pos->second.begin();mode!=pos->second.end();++mode)
      {
         cout<<"      ";
         (*mode)->Print(cout);
         cout<<endl;
      }
   }
   for(map<const MolBondAngle*,set<const StretchMode*> >::const_iterator pos=vpAngle.begin();pos!=vpAngle.end();++pos)
   {
      cout<<"Bond Angle "<<pos->first->GetName()<<" is modified by the stretch modes:"<<endl;
      for(set<const StretchMode*>::const_iterator mode=pos->second.begin();mode!=pos->second.end();++mode)
      {
         cout<<"      ";
         (*mode)->Print(cout);
         cout<<endl;
      }
   }
   for(map<const MolDihedralAngle*,set<const StretchMode*> >::const_iterator pos=vpDihed.begin();pos!=vpDihed.end();++pos)
   {
      cout<<"Dihedral Angle "<<pos->first->GetName()<<" is modified by the stretch modes:"<<endl;
      for(set<const StretchMode*>::const_iterator mode=pos->second.begin();mode!=pos->second.end();++mode)
      {
         cout<<"      ";
         (*mode)->Print(cout);
         cout<<endl;
      }
   }
   #endif
   mvpStretchModeFree.clear();
   mvpStretchModeNotFree.clear();
   
   for(list<StretchModeBondLength>::iterator mode=mvStretchModeBondLength.begin();
       mode!=mvStretchModeBondLength.end();++mode)
   {
      int nb=mode->mvpBrokenDihedralAngle.size()+mode->mvpBrokenBondAngle.size()+mode->mvpBrokenBond.size();
      if(mode->mpBond!=0) nb -= 1;
      if(nb==0) mvpStretchModeFree.push_back(&(*mode));
      else mvpStretchModeNotFree.push_back(&(*mode));
   }
   
   for(list<StretchModeBondAngle>::iterator mode=mvStretchModeBondAngle.begin();
       mode!=mvStretchModeBondAngle.end();++mode)
   {
      int nb=mode->mvpBrokenDihedralAngle.size()+mode->mvpBrokenBondAngle.size()+mode->mvpBrokenBond.size();
      if(mode->mpBondAngle!=0) nb -= 1;
      if(nb==0) mvpStretchModeFree.push_back(&(*mode));
      else mvpStretchModeNotFree.push_back(&(*mode));
   }

   for(list<StretchModeTorsion>::iterator mode=mvStretchModeTorsion.begin();
       mode!=mvStretchModeTorsion.end();++mode)
   {
      int nb=mode->mvpBrokenDihedralAngle.size()+mode->mvpBrokenBondAngle.size()+mode->mvpBrokenBond.size();
      if(mode->mpDihedralAngle!=0) nb -= 1;
      if(nb==0) mvpStretchModeFree.push_back(&(*mode));
      else mvpStretchModeNotFree.push_back(&(*mode));
   }
   #if 1
   for(list<StretchModeTwist>::iterator mode=mvStretchModeTwist.begin();
       mode!=mvStretchModeTwist.end();++mode)
      if(mode->mvpBrokenDihedralAngle.size()+mode->mvpBrokenBondAngle.size()+mode->mvpBrokenBond.size()==0)
         mvpStretchModeFree.push_back(&(*mode));
      else mvpStretchModeNotFree.push_back(&(*mode));
   #endif
}

void Molecule::BuildMDAtomGroups()
{
   // For each atom, list all atoms that are never moved relatively to it.
   // (not moved == distance cannot change)
   map<MolAtom*,set<MolAtom*> > vBoundAtoms;
   set<MolAtom*> set0;
   for(vector<MolAtom*>::iterator pos=this->GetAtomList().begin();pos!=this->GetAtomList().end();++pos)
      set0.insert(*pos);
   
   for(vector<MolAtom*>::iterator pat1=this->GetAtomList().begin();pat1!=this->GetAtomList().end();++pat1)
   {
      vBoundAtoms[*pat1]=set0;
      for(vector<MolAtom*>::iterator pat2=this->GetAtomList().begin();pat2!=this->GetAtomList().end();++pat2)
      {
         bool cont=false;
         
         for(list<StretchModeBondLength>::iterator pstretch=this->GetStretchModeBondLengthList().begin();
            pstretch!=this->GetStretchModeBondLengthList().end();++pstretch)
         {
            set<MolAtom *>::iterator pos1=pstretch->mvTranslatedAtomList.find(*pat1),
                                     pos2=pstretch->mvTranslatedAtomList.find(*pat2);
            if(  ((pos1==pstretch->mvTranslatedAtomList.end())&&(pos2!=pstretch->mvTranslatedAtomList.end()))
               ||((pos1!=pstretch->mvTranslatedAtomList.end())&&(pos2==pstretch->mvTranslatedAtomList.end())))
            {
               vBoundAtoms[*pat1].erase(*pat2);
               //cout<<(*pat1)->GetName()<<" moves (b) relatively to "<<(*pat2)->GetName()<<"  /";
               //pstretch->Print(cout);cout<<endl;
               cont=true;
               break;
            }
         }
         if(cont) continue;
         
         for(list<StretchModeBondAngle>::iterator pstretch=this->GetStretchModeBondAngleList().begin();
            pstretch!=this->GetStretchModeBondAngleList().end();++pstretch)
         {
            //pstretch->mpAtom1 does not move relatively to any atom
            if((*pat1==pstretch->mpAtom1)||(*pat2==pstretch->mpAtom1)) continue;
            
            set<MolAtom *>::iterator pos1=pstretch->mvRotatedAtomList.find(*pat1),
                                     pos2=pstretch->mvRotatedAtomList.find(*pat2);
            //:TODO: Take into account the special case of the atoms defining the bond angle
            if(  ((pos1==pstretch->mvRotatedAtomList.end())&&(pos2!=pstretch->mvRotatedAtomList.end()))
               ||((pos1!=pstretch->mvRotatedAtomList.end())&&(pos2==pstretch->mvRotatedAtomList.end())))
            {
               vBoundAtoms[*pat1].erase(*pat2);
               //cout<<(*pat1)->GetName()<<" moves (a) relatively to "<<(*pat2)->GetName()<<"  /";
               //pstretch->Print(cout);cout<<endl;
               cont=true;
               break;
            }
         }
         if(cont) continue;
         
         for(list<StretchModeTorsion>::iterator pstretch=this->GetStretchModeTorsionList().begin();
            pstretch!=this->GetStretchModeTorsionList().end();++pstretch)
         {
            set<MolAtom *>::iterator pos1=pstretch->mvRotatedAtomList.find(*pat1),
                                     pos2=pstretch->mvRotatedAtomList.find(*pat2);
            
            if(  (pos1!=pstretch->mvRotatedAtomList.end())&&(pos2==pstretch->mvRotatedAtomList.end())
               &&(*pat2!=pstretch->mpAtom1) &&(*pat2!=pstretch->mpAtom2) )
            {
               vBoundAtoms[*pat1].erase(*pat2);
               //cout<<(*pat1)->GetName()<<" moves (d1) relatively to "<<(*pat2)->GetName()<<"  /";
               //pstretch->Print(cout);cout<<endl;
               break;
            }
            if(  (pos1==pstretch->mvRotatedAtomList.end())&&(pos2!=pstretch->mvRotatedAtomList.end())
               &&(*pat1!=pstretch->mpAtom1) && (*pat1!=pstretch->mpAtom2) )
            {
               vBoundAtoms[*pat1].erase(*pat2);
               //cout<<(*pat1)->GetName()<<" moves (d2) relatively to "<<(*pat2)->GetName()<<"  /";
               //pstretch->Print(cout);cout<<endl;
               break;
            }
         }
      }
   }
   
   // List remaining group of atoms, take care of rigid groups
   set<set<MolAtom*> > vBoundGroups;
   for(map<MolAtom*,set<MolAtom*> >::iterator pos=vBoundAtoms.begin();pos!=vBoundAtoms.end();++pos)
   {
      #if 0
      cout<<"Non-flexible group from "<<pos->first->GetName()<<": ";
      for(set<MolAtom*>::const_iterator atom=pos->second.begin();atom!=pos->second.end();++atom)
         cout<<(*atom)->GetName()<<",";
      cout<<endl;
      #endif
      // Remove atoms belonging to a rigid group
      for(vector<RigidGroup *>::iterator pr=this->GetRigidGroupList().begin();pr!=this->GetRigidGroupList().end();++pr)
         for(set<MolAtom *>::iterator at=(*pr)->begin();at!=(*pr)->end();++at)
            pos->second.erase(*at);
      if(pos->second.size()>1) vBoundGroups.insert(pos->second);
   }
   #if 0
   for(set<set<MolAtom*> >::iterator pos=vBoundGroups.begin();pos!=vBoundGroups.end();++pos)
   {
      cout<<"Non-flexible group:";
      for(set<MolAtom*>::const_iterator atom=pos->begin();atom!=pos->end();++atom)
         cout<<(*atom)->GetName()<<",";
      cout<<endl;
   }
   #endif
   // Create relevant MDAtomGroup, listing associated restraints
   mvMDAtomGroup.clear();
   for(set<set<MolAtom*> >::iterator pos=vBoundGroups.begin();pos!=vBoundGroups.end();++pos)
   {
      set<MolBond*> vb;
      for(vector<MolBond*>::iterator pr=this->GetBondList().begin();pr!=this->GetBondList().end();++pr)
         if(  (pos->find(&(*pr)->GetAtom1())!=pos->end())
            ||(pos->find(&(*pr)->GetAtom2())!=pos->end())) vb.insert(*pr);
      
      set<MolBondAngle*> va;
      for(vector<MolBondAngle*>::iterator pr=this->GetBondAngleList().begin();pr!=this->GetBondAngleList().end();++pr)
         if(  (pos->find(&(*pr)->GetAtom1())!=pos->end())
            ||(pos->find(&(*pr)->GetAtom2())!=pos->end())
            ||(pos->find(&(*pr)->GetAtom3())!=pos->end())) va.insert(*pr);
      
      set<MolDihedralAngle*> vd;
      for(vector<MolDihedralAngle*>::iterator pr=this->GetDihedralAngleList().begin();pr!=this->GetDihedralAngleList().end();++pr)
         if(  (pos->find(&(*pr)->GetAtom1())!=pos->end())
            ||(pos->find(&(*pr)->GetAtom2())!=pos->end())
            ||(pos->find(&(*pr)->GetAtom3())!=pos->end())
            ||(pos->find(&(*pr)->GetAtom4())!=pos->end())) vd.insert(*pr);
      
      set<MolAtom*> tmp= *pos;// Cannot pass *pos directly ? gcc bug evaluating const-ness?
      mvMDAtomGroup.push_back(MDAtomGroup(tmp,vb,va,vd));
      mvMDAtomGroup.back().Print(cout);
   }
   
   // Create mvMDFullAtomGroup
   mvMDFullAtomGroup.clear();
   #if 1
   // All atoms except those in rigid groups
   for(vector<MolAtom*>::iterator at=this->GetAtomList().begin();at!=this->GetAtomList().end();++at)
      mvMDFullAtomGroup.insert(*at);
   for(vector<RigidGroup *>::iterator pr=this->GetRigidGroupList().begin();pr!=this->GetRigidGroupList().end();++pr)
      for(set<MolAtom *>::iterator at=(*pr)->begin();at!=(*pr)->end();++at)
         mvMDFullAtomGroup.erase(*at);
   cout<<"Full MD atom group:"<<endl<<" ";
   for(set<MolAtom*>::const_iterator pos=mvMDFullAtomGroup.begin();pos!=mvMDFullAtomGroup.end();++pos)
      cout<<(*pos)->GetName()<<" ";
   cout<<endl;
   #else
   // All atoms listed in at leat one mvMDAtomGroup
   mvMDFullAtomGroup.clear();
   for(list<MDAtomGroup>::const_iterator pos= mvMDAtomGroup.begin();pos!= mvMDAtomGroup.end();++pos)
      for(set<MolAtom*>::const_iterator at=pos->mvpAtom.begin();at!=pos->mvpAtom.end();++at)
         mvMDFullAtomGroup.insert(*at);
   cout<<"Full MD atom group:"<<endl<<" ";
   for(set<MolAtom*>::const_iterator pos=mvMDFullAtomGroup.begin();pos!=mvMDFullAtomGroup.end();++pos)
      cout<<(*pos)->GetName()<<" ";
   cout<<endl;
   #endif
   mClockMDAtomGroup.Click();
}

void Molecule::UpdateScattCompList()const
{
   if(  (mClockAtomPosition<mClockScattCompList)
      &&(mClockOrientation <mClockScattCompList)
      &&(mClockAtomScattPow<mClockScattCompList)
      &&(mClockScatterer   <mClockScattCompList))return;
   VFN_DEBUG_ENTRY("Molecule::UpdateScattCompList()",5)
   TAU_PROFILE("Molecule::UpdateScattCompList()","void ()",TAU_DEFAULT);
   const long nb=this->GetNbComponent();
   // Get internal coords
   for(long i=0;i<nb;++i)
   {
      if(mvpAtom[i]->IsDummy()) mScattCompList(i).mpScattPow=0;
      else mScattCompList(i).mpScattPow=&(mvpAtom[i]->GetScatteringPower());
      mScattCompList(i).mX=mvpAtom[i]->GetX();
      mScattCompList(i).mY=mvpAtom[i]->GetY();
      mScattCompList(i).mZ=mvpAtom[i]->GetZ();
      mScattCompList(i).mOccupancy=mvpAtom[i]->GetOccupancy()*mOccupancy;
   }
  
  #ifdef RIGID_BODY_STRICT_EXPERIMENTAL
  // During an optimization, apply the translations & rotations of the rigid group parameters
  if(true)//this->IsBeingRefined())
  {
    for(vector<RigidGroup *>::const_iterator pos=this->GetRigidGroupList().begin();pos!=this->GetRigidGroupList().end();++pos)
    {
        (*pos)->mQuat.Normalize();
        // Center of the atom group
        REAL x0=0,y0=0,z0=0;
        for(set<unsigned int>::iterator at=(*pos)->mvIdx.begin();at!=(*pos)->mvIdx.end();++at)
        {
          x0+=mvpAtom[*at]->GetX();
          y0+=mvpAtom[*at]->GetY();
          z0+=mvpAtom[*at]->GetZ();
        }
        x0/=(*pos)->size();
        y0/=(*pos)->size();
        z0/=(*pos)->size();
        
        // Apply rotation & translation to all atoms
        for(set<unsigned int>::iterator at=(*pos)->mvIdx.begin();at!=(*pos)->mvIdx.end();++at)
        {
          REAL x=mvpAtom[*at]->GetX()-x0, y=mvpAtom[*at]->GetY()-y0, z=mvpAtom[*at]->GetZ()-z0;
          (*pos)->mQuat.RotateVector(x,y,z);
          mScattCompList(*at).mX=x+x0+(*pos)->mX;
          mScattCompList(*at).mY=y+y0+(*pos)->mY;
          mScattCompList(*at).mZ=z+z0+(*pos)->mZ;
        }
    }
  }
  #endif
   // translate center to (0,0,0)
   REAL x0=0,y0=0,z0=0;
   if((mMoleculeCenter.GetChoice()==0) || (mpCenterAtom==0))
   {
      for(long i=0;i<nb;++i)
      {
         x0 += mScattCompList(i).mX;
         y0 += mScattCompList(i).mY;
         z0 += mScattCompList(i).mZ;
      }
      x0 /= nb;
      y0 /= nb;
      z0 /= nb;
   }
   else
   {
      x0=mpCenterAtom->GetX();
      y0=mpCenterAtom->GetY();
      z0=mpCenterAtom->GetZ();
   }
   for(long i=0;i<nb;++i)
   {
      mScattCompList(i).mX -= x0;
      mScattCompList(i).mY -= y0;
      mScattCompList(i).mZ -= z0;
   }
   //VFN_DEBUG_MESSAGE("Molecule::UpdateScattCompList()",10)
   // rotate
   mQuat.Normalize();
   for(long i=0;i<nb;++i)
   {
      //#error the vector must not be normalized !
      mQuat.RotateVector(mScattCompList(i).mX,mScattCompList(i).mY,mScattCompList(i).mZ);
   }
   // Convert to fractionnal coordinates
   for(long i=0;i<nb;++i)
   {
      this->GetCrystal().OrthonormalToFractionalCoords(mScattCompList(i).mX,
                                                       mScattCompList(i).mY,
                                                       mScattCompList(i).mZ);
   }
   // translate center to position in unit cell
   for(long i=0;i<nb;++i)
   {
      mScattCompList(i).mX += mXYZ(0);
      mScattCompList(i).mY += mXYZ(1);
      mScattCompList(i).mZ += mXYZ(2);
   }
   mClockScattCompList.Click();
   VFN_DEBUG_EXIT("Molecule::UpdateScattCompList()",5)
}
vector<MolAtom*>::reverse_iterator Molecule::FindAtom(const string &name)
{
   VFN_DEBUG_ENTRY("Molecule::FindAtom():"<<name,4)
   vector<MolAtom*>::reverse_iterator rpos;
   for(rpos=mvpAtom.rbegin();rpos!=mvpAtom.rend();++rpos)
      if(name==(*rpos)->GetName())
      {
         VFN_DEBUG_EXIT("Molecule::FindAtom():"<<name<<"...NOT FOUND !",4)
         return rpos;
      }
   VFN_DEBUG_EXIT("Molecule::FindAtom():"<<name<<"...NOT FOUND !",4)
   return rpos;
}
vector<MolAtom*>::const_reverse_iterator Molecule::FindAtom(const string &name)const
{
   vector<MolAtom*>::const_reverse_iterator rpos;
   rpos=mvpAtom.rbegin();
   for(rpos=mvpAtom.rbegin();rpos!=mvpAtom.rend();++rpos)
      if(name==(*rpos)->GetName()) return rpos;
   return rpos;
}
void Molecule::InitOptions()
{
   VFN_DEBUG_ENTRY("Molecule::InitOptions",7)
   static string Flexname;
   static string Flexchoices[3];

   static string autoOptimizeConformationName;
   static string autoOptimizeConformationChoices[2];

   static string optimizeOrientationName;
   static string optimizeOrientationChoices[2];
   
   static string moleculeCenterName;
   static string moleculeCenterChoices[2];
   
   static bool needInitNames=true;
   if(true==needInitNames)
   {
      Flexname="Flexibility Model";
      Flexchoices[0]="Automatic from Restraints, relaxed - RECOMMENDED";
      Flexchoices[1]="Rigid Body";
      Flexchoices[2]="Automatic from Restraints, strict";
      //Flexchoices[3]="Molecular Dynamics";
      
      autoOptimizeConformationName="Auto Optimize Starting Conformation";
      autoOptimizeConformationChoices[0]="Yes";
      autoOptimizeConformationChoices[1]="No";

      optimizeOrientationName="Optimize Orientation";
      optimizeOrientationChoices[0]="Yes";
      optimizeOrientationChoices[1]="No";
      
      moleculeCenterName="Rotation Center";
      moleculeCenterChoices[0]="Geometrical center (recommended)";
      moleculeCenterChoices[1]="User-chosen Atom";
      
      needInitNames=false;
   }
   mFlexModel.Init(3,&Flexname,Flexchoices);
   mFlexModel.SetChoice(0);
   this->AddOption(&mFlexModel);
   
   mAutoOptimizeConformation.Init(2,&autoOptimizeConformationName,
                                  autoOptimizeConformationChoices);
   this->AddOption(&mAutoOptimizeConformation);

   mOptimizeOrientation.Init(2,&optimizeOrientationName,optimizeOrientationChoices);
   this->AddOption(&mOptimizeOrientation);

   mMoleculeCenter.Init(2,&moleculeCenterName,moleculeCenterChoices);
   this->AddOption(&mMoleculeCenter);
   
   VFN_DEBUG_EXIT("Molecule::InitOptions",7)
}

Molecule::FlipGroup::FlipGroup(const MolAtom &at0,const MolAtom &at1,const MolAtom &at2):
mpAtom0(&at0),mpAtom1(&at1),mpAtom2(&at2),mNbTest(0),mNbAccept(0)
{
}

void Molecule::FlipAtomGroup(const FlipGroup& group, const bool keepCenter)
{
   TAU_PROFILE("Molecule::FlipAtomGroup(FlipGroup&)","void (...)",TAU_DEFAULT);
   if(group.mpAtom0==group.mvRotatedChainList.back().first)
   {// We are doing a 180� rotation exchanging two bonds
      const REAL vx=group.mpAtom0->X()-(group.mpAtom1->X()+group.mpAtom2->X())/2.;
      const REAL vy=group.mpAtom0->Y()-(group.mpAtom1->Y()+group.mpAtom2->Y())/2.;
      const REAL vz=group.mpAtom0->Z()-(group.mpAtom1->Z()+group.mpAtom2->Z())/2.;
      this->RotateAtomGroup(*(group.mpAtom0),vx,vy,vz,
                            group.mvRotatedChainList.back().second,M_PI,keepCenter);
   }
   else
   {// we are flipping bonds with respect to a plane defined by other bonds
      REAL v01x=group.mpAtom1->X()-group.mpAtom0->X();
      REAL v01y=group.mpAtom1->Y()-group.mpAtom0->Y();
      REAL v01z=group.mpAtom1->Z()-group.mpAtom0->Z();
      const REAL norm01=sqrt(v01x*v01x+v01y*v01y+v01z*v01z+1e-7);
      v01x /= norm01;v01y /= norm01;v01z /= norm01;

      REAL v02x=group.mpAtom2->X()-group.mpAtom0->X();
      REAL v02y=group.mpAtom2->Y()-group.mpAtom0->Y();
      REAL v02z=group.mpAtom2->Z()-group.mpAtom0->Z();
      const REAL norm02=sqrt(v02x*v02x+v02y*v02y+v02z*v02z+1e-7);
      v02x /= norm02;v02y /= norm02;v02z /= norm02;
      
      REAL v12x=group.mpAtom2->X()-group.mpAtom1->X();
      REAL v12y=group.mpAtom2->Y()-group.mpAtom1->Y();
      REAL v12z=group.mpAtom2->Z()-group.mpAtom1->Z();
      const REAL norm12=sqrt(v12x*v12x+v12y*v12y+v12z*v12z+1e-7);
      v12x /= norm12;v12y /= norm12;v12z /= norm12;
      
      REAL v0mx=group.mpAtom0->X()-(group.mpAtom1->X()+group.mpAtom2->X())/2.;
      REAL v0my=group.mpAtom0->Y()-(group.mpAtom1->Y()+group.mpAtom2->Y())/2.;
      REAL v0mz=group.mpAtom0->Z()-(group.mpAtom1->Z()+group.mpAtom2->Z())/2.;
      const REAL norm0m=sqrt(v0mx*v0mx+v0my*v0my+v0mz*v0mz+1e-7);
      v0mx /= norm0m;v0my /= norm0m;v0mz /= norm0m;

      if(fabs(v01x*v02x+v01y*v02y+v01z*v02z)
         >0.05*sqrt(abs( (v01x*v01x+v01y*v01y+v01z*v01z)
                        *(v02x*v02x+v02y*v02y+v02z*v02z))))
      {
         REAL v012x=v01y*v02z-v01z*v02y;
         REAL v012y=v01z*v02x-v01x*v02z;
         REAL v012z=v01x*v02y-v01y*v02x;
         const REAL norm012=sqrt(v012x*v012x+v012y*v012y+v012z*v012z+1e-7);
         v012x /= norm012;v012y /= norm012;v012z /= norm012;
         
      
         for(list<pair<const MolAtom *,set<MolAtom*> > >::const_iterator
             chain=group.mvRotatedChainList.begin();
             chain!=group.mvRotatedChainList.end();++chain)
         {
            REAL v03x=chain->first->X()-group.mpAtom0->X();
            REAL v03y=chain->first->Y()-group.mpAtom0->Y();
            REAL v03z=chain->first->Z()-group.mpAtom0->Z();
            const REAL norm03=sqrt( v03x*v03x + v03y*v03y + v03z*v03z +1e-7);
            v03x /= norm03;v03y /= norm03;v03z /= norm03;

            const REAL a1=v012x*v03x+v012y*v03y+v012z*v03z;
            const REAL a2= v0mx*v03x+ v0my*v03y+ v0mz*v03z;
            const REAL a3= v12x*v03x+ v12y*v03y+ v12z*v03z;
            REAL angle = -a1/sqrt(1-a3*a3+1e-7);
            if(angle>=1.)
               angle = M_PI/2.;
            else
            {
               if(angle<=-1.)
               {
                  angle = -M_PI/2.;
               }
               else angle = asin(angle);
            }
            if(a2<0) angle=M_PI-angle;
            this->RotateAtomGroup(*(group.mpAtom0),v12x,v12y,v12z,
                                  chain->second,2*angle,keepCenter);
         }
      }
   }
}

void Molecule::SetDeleteSubObjInDestructor(const bool b) {
    mDeleteSubObjInDestructor=b;
}

void Molecule::ResetRigidGroupsPar()const
{
   #ifdef RIGID_BODY_STRICT_EXPERIMENTAL
   // Apply the translations & rotations of all rigid group parameters, and
   // use this as the newly stored atomic coordinates.
   for(vector<RigidGroup *>::const_iterator pos=this->GetRigidGroupList().begin();pos!=this->GetRigidGroupList().end();++pos)
   {
      (*pos)->mQuat.Normalize();
      // Center of atom group
      REAL x0=0,y0=0,z0=0;
      for(set<MolAtom *>::iterator at=(*pos)->begin();at!=(*pos)->end();++at)
      {
        x0+=(*at)->GetX();
        y0+=(*at)->GetY();
        z0+=(*at)->GetZ();
      }
      x0/=(*pos)->size();
      y0/=(*pos)->size();
      z0/=(*pos)->size();
      
      // Apply rotation & translation to all atoms
      for(set<MolAtom *>::iterator at=(*pos)->begin();at!=(*pos)->end();++at)
      {
        REAL x=(*at)->GetX()-x0, y=(*at)->GetY()-y0, z=(*at)->GetZ()-z0;
        (*pos)->mQuat.RotateVector(x,y,z);
        (*at)->SetX(x+x0+(*pos)->mX);
        (*at)->SetY(y+y0+(*pos)->mY);
        (*at)->SetZ(z+z0+(*pos)->mZ);
      }
      
      // Reset the translation & rotation parameters, only useful during an optimization
      (*pos)->mX=0;
      (*pos)->mY=0;
      (*pos)->mZ=0;
      (*pos)->mQuat.Q0()=1;
      (*pos)->mQuat.Q1()=0;
      (*pos)->mQuat.Q2()=0;
      (*pos)->mQuat.Q3()=0;
   }
   #endif
}

#ifdef __WX__CRYST__
WXCrystObjBasic* Molecule::WXCreate(wxWindow* parent)
{
   VFN_DEBUG_ENTRY("Molecule::WXCreate()",5)
   mpWXCrystObj=new WXMolecule(parent,this);
   VFN_DEBUG_EXIT("Molecule::WXCreate()",5)
   return mpWXCrystObj;
}
#endif

}//namespace
