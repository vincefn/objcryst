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
/*   Molecule.cpp
*  Source file for the Molecule scatterer
*
*/
#include <sstream>
#include <iterator>
#include <algorithm>

#include "Quirks/VFNStreamFormat.h"
#include "ObjCryst/Molecule.h"
#include "RefinableObj/GlobalOptimObj.h"

#ifdef OBJCRYST_GL
   #ifdef __DARWIN__
      #include <OpenGL/glu.h>
   #else
      #include <GL/glu.h>
   #endif
#endif

#ifdef __WX__CRYST__
   #include "wxCryst/wxMolecule.h"
#endif

using namespace std;

namespace ObjCryst
{
REAL GetBondLength(const MolAtom&at1,const MolAtom&at2)
{
   //TAU_PROFILE("GetBondLength()","REAL (...)",TAU_DEFAULT);
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
   const REAL norm21= x21*x21+y21*y21+z21*z21;
   const REAL norm23= x23*x23+y23*y23+z23*z23;
   const REAL angle=(x21*x23+y21*y23+z21*z23)/sqrt(norm21*norm23+1e-6);
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
   const REAL norm123= x123*x123+y123*y123+z123*z123;
   
   // v32 x v34 (= -v23 x v34)
   const REAL x234= -(y23*z34-z23*y34);
   const REAL y234= -(z23*x34-x23*z34);
   const REAL z234= -(x23*y34-y23*x34);
   const REAL norm234= x234*x234+y234*y234+z234*z234;
   
   REAL angle=(x123*x234+y123*y234+z123*z234)/sqrt(norm123*norm234+1e-6);
   if(angle>= 1) angle=0;
   else 
   {
      if(angle<=-1) angle=M_PI;
      else angle=acos(angle);
   }
   if((x21*x34 + y21*y34 + z21*z34)<0) return -angle;
   return angle;
}

//######################################################################
//
//      MolAtom
//
//######################################################################

MolAtom::MolAtom(const REAL x, const REAL y, const REAL z,
                 const ScatteringPower *pPow, const string &name, Molecule &parent):
mName(name),mX(x),mY(y),mZ(z),mOccupancy(1.),mpScattPow(pPow),mpMol(&parent)
#ifdef __WX__CRYST__
,mpWXCrystObj(0)
#endif
{
   VFN_DEBUG_MESSAGE("MolAtom::MolAtom()",4)
}

MolAtom::~MolAtom(){}

void MolAtom::SetName(const string &name){mName=name;}
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

void MolAtom::SetX(const REAL a){ mX=a;}
void MolAtom::SetY(const REAL a){ mY=a;}
void MolAtom::SetZ(const REAL a){ mZ=a;}
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
      ss <<mX;
      tag.AddAttribute("x",ss.str());
   }
   {
      stringstream ss;
      ss <<mY;
      tag.AddAttribute("y",ss.str());
   }
   {
      stringstream ss;
      ss <<mZ;
      tag.AddAttribute("z",ss.str());
   }
   {
      stringstream ss;
      ss <<mOccupancy;
      tag.AddAttribute("Occup",ss.str());
   }
   os <<tag<<endl;
   VFN_DEBUG_EXIT("MolAtom::XMLOutput()",4)
}

void MolAtom::XMLInput(istream &is,const XMLCrystTag &tag)
{
   VFN_DEBUG_ENTRY("MolAtom::XMLInput()",10)
   for(unsigned int i=0;i<tag.GetNbAttribute();i++)
   {
      if("Name"==tag.GetAttributeName(i))
      {
         mName=tag.GetAttributeValue(i);
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
   VFN_DEBUG_EXIT("MolAtom::XMLInput()",10)
}

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
MolBond::MolBond(const MolAtom &atom1, const MolAtom &atom2,
                 const REAL length0, const REAL sigma, const REAL delta,
                 Molecule &parent,const REAL bondOrder):
mAtomPair(make_pair(&atom1,&atom2)),
mLength0(length0),mDelta(delta),mSigma(sigma),
mBondOrder(bondOrder),mIsFreeTorsion(false),mpMol(&parent)
{}

MolBond::~MolBond()
{}

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
      ss <<mLength0;
      tag.AddAttribute("Length",ss.str());
   }
   {
      stringstream ss;
      ss <<mDelta;
      tag.AddAttribute("Delta",ss.str());
   }
   {
      stringstream ss;
      ss <<mSigma;
      tag.AddAttribute("Sigma",ss.str());
   }
   {
      stringstream ss;
      ss <<mBondOrder;
      tag.AddAttribute("BondOrder",ss.str());
   }
   {
      stringstream ss;
      ss <<mIsFreeTorsion;
      tag.AddAttribute("FreeTorsion",ss.str());
   }
   os <<tag<<endl;
   VFN_DEBUG_EXIT("MolBond::XMLOutput()",4)
}

void MolBond::XMLInput(istream &is,const XMLCrystTag &tag)
{
   VFN_DEBUG_ENTRY("MolBond::XMLInput():",10)
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
   VFN_DEBUG_EXIT("MolBond::XMLInput():",10)
}

REAL MolBond::GetLogLikelihood()const
{
   VFN_DEBUG_ENTRY("MolBond::GetLogLikelihood():",2)
   const REAL length=this->GetLength();
   REAL tmp=length-(mLength0+mDelta);
   if(tmp>0)
   {
      tmp/=mSigma;
      //tmp=fabs(tmp)/mSigma;
      //if(tmp>30) return 2.8550185e25*(tmp-29);
      //tmp=sinh(tmp);
      VFN_DEBUG_EXIT("MolBond::GetLogLikelihood():",2)
      return tmp*tmp;
   }
   tmp=length-(mLength0-mDelta);
   if(tmp<0)
   {
      tmp/=mSigma;
      //tmp=fabs(tmp)/mSigma;
      //if(tmp>30) return 2.8550185e25*(tmp-29);
      //tmp=sinh(tmp);
      VFN_DEBUG_EXIT("MolBond::GetLogLikelihood():",2)
      return tmp*tmp;
   }
   VFN_DEBUG_EXIT("MolBond::GetLogLikelihood():",2)
   return 0;
}

const MolAtom& MolBond::GetAtom1()const{return *(mAtomPair.first);}
const MolAtom& MolBond::GetAtom2()const{return *(mAtomPair.second);}
void MolBond::SetAtom1(const MolAtom &at){mAtomPair.first =&at;}
void MolBond::SetAtom2(const MolAtom &at){mAtomPair.second=&at;}
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
void MolBond::SetLengthSigma(const REAL a){mBondOrder=a;}
void MolBond::SetBondOrder(const REAL a){mSigma=a;}

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
MolBondAngle::MolBondAngle(const MolAtom &atom1,const  MolAtom &atom2,const  MolAtom &atom3,
                           const REAL angle, const REAL sigma, const REAL delta,
                           Molecule &parent):
mAngle0(angle),mDelta(delta),mSigma(sigma),mpMol(&parent)
{
   mvpAtom.push_back(&atom1);
   mvpAtom.push_back(&atom2);
   mvpAtom.push_back(&atom3);
}

MolBondAngle::~MolBondAngle(){}

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
      ss <<mAngle0*RAD2DEG;
      tag.AddAttribute("Angle",ss.str());
   }
   {
      stringstream ss;
      ss <<mDelta*RAD2DEG;
      tag.AddAttribute("Delta",ss.str());
   }
   {
      stringstream ss;
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

REAL MolBondAngle::GetAngle()const
{
   return GetBondAngle(this->GetAtom1(),this->GetAtom2(),this->GetAtom3());
}

REAL MolBondAngle::GetLogLikelihood()const
{
   VFN_DEBUG_ENTRY("MolBondAngle::GetLogLikelihood():",2)
   const REAL angle=this->GetAngle();
   REAL tmp=angle-(mAngle0+mDelta);
   if(tmp>0)
   {
      tmp/=mSigma;
      //tmp=fabs(tmp)/mSigma;
      //if(tmp>30) return 2.8550185e25*(tmp-29);
      //tmp=sinh(tmp);
      VFN_DEBUG_EXIT("MolBondAngle::GetLogLikelihood():",2)
      return tmp*tmp;
   }
   tmp=angle-(mAngle0-mDelta);
   if(tmp<0)
   {
      tmp/=mSigma;
      //tmp=fabs(tmp)/mSigma;
      //if(tmp>30) return 2.8550185e25*(tmp-29);
      //tmp=sinh(tmp);
      VFN_DEBUG_EXIT("MolBondAngle::GetLogLikelihood():",2)
      return tmp*tmp;
   }
   VFN_DEBUG_EXIT("MolBondAngle::GetLogLikelihood():",2)
   return 0;
}
const MolAtom& MolBondAngle::GetAtom1()const{return *(mvpAtom[0]);}
const MolAtom& MolBondAngle::GetAtom2()const{return *(mvpAtom[1]);}
const MolAtom& MolBondAngle::GetAtom3()const{return *(mvpAtom[2]);}
void MolBondAngle::SetAtom1(const MolAtom& at){mvpAtom[0]=&at;}
void MolBondAngle::SetAtom2(const MolAtom& at){mvpAtom[1]=&at;}
void MolBondAngle::SetAtom3(const MolAtom& at){mvpAtom[2]=&at;}
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
MolDihedralAngle::MolDihedralAngle(const MolAtom &atom1, const MolAtom &atom2,
                                   const MolAtom &atom3, const MolAtom &atom4,
                                   const REAL angle, const REAL sigma, const REAL delta,
                                   Molecule &parent):
mAngle0(angle),mDelta(delta),mSigma(sigma),mpMol(&parent)
{
   VFN_DEBUG_ENTRY("MolDihedralAngle::MolDihedralAngle()",5)
   mvpAtom.push_back(&atom1);
   mvpAtom.push_back(&atom2);
   mvpAtom.push_back(&atom3);
   mvpAtom.push_back(&atom4);
   VFN_DEBUG_EXIT("MolDihedralAngle::MolDihedralAngle()",5)
}

MolDihedralAngle::~MolDihedralAngle(){}

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
      ss <<mAngle0*RAD2DEG;
      tag.AddAttribute("Angle",ss.str());
   }
   {
      stringstream ss;
      ss <<mDelta*RAD2DEG;
      tag.AddAttribute("Delta",ss.str());
   }
   {
      stringstream ss;
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
   return GetDihedralAngle(this->GetAtom1(),this->GetAtom2(),this->GetAtom3(),this->GetAtom4());
}

REAL& MolDihedralAngle::Angle0(){return mAngle0;}

REAL MolDihedralAngle::GetLogLikelihood()const
{
   VFN_DEBUG_ENTRY("MolDihedralAngle::GetLogLikelihood():",2)
   const REAL angle=this->GetAngle();
   REAL tmp=angle-(mAngle0+mDelta);
   if(fabs(tmp+2*M_PI)<fabs(tmp)) tmp += 2*M_PI;
   if(fabs(tmp-2*M_PI)<fabs(tmp)) tmp -= 2*M_PI;
   if(tmp>0)
   {
      tmp/=mSigma;
      //tmp=fabs(tmp)/mSigma;
      //if(tmp>30) return 2.8550185e25*(tmp-29);
      //tmp=sinh(tmp);
      VFN_DEBUG_EXIT("MolDihedralAngle::GetLogLikelihood():",2)
      return tmp*tmp;
   }
   tmp=angle-(mAngle0-mDelta);
   if(fabs(tmp+2*M_PI)<fabs(tmp)) tmp=tmp+2*M_PI;
   if(fabs(tmp-2*M_PI)<fabs(tmp)) tmp=tmp-2*M_PI;
   if(tmp<0)
   {
      tmp/=mSigma;
      //tmp=fabs(tmp)/mSigma;
      //if(tmp>30) return 2.8550185e25*(tmp-29);
      //tmp=sinh(tmp);
      VFN_DEBUG_EXIT("MolDihedralAngle::GetLogLikelihood():",2)
      return tmp*tmp;
   }
   VFN_DEBUG_EXIT("MolDihedralAngle::GetLogLikelihood():",2)
   return 0;
}

const MolAtom& MolDihedralAngle::GetAtom1()const{return *(mvpAtom[0]);}
const MolAtom& MolDihedralAngle::GetAtom2()const{return *(mvpAtom[1]);}
const MolAtom& MolDihedralAngle::GetAtom3()const{return *(mvpAtom[2]);}
const MolAtom& MolDihedralAngle::GetAtom4()const{return *(mvpAtom[3]);}
void MolDihedralAngle::SetAtom1(const MolAtom& at){mvpAtom[0]=&at;}
void MolDihedralAngle::SetAtom2(const MolAtom& at){mvpAtom[1]=&at;}
void MolDihedralAngle::SetAtom3(const MolAtom& at){mvpAtom[2]=&at;}
void MolDihedralAngle::SetAtom4(const MolAtom& at){mvpAtom[3]=&at;}
//MolAtom& MolDihedralAngle::GetAtom1();
//MolAtom& MolDihedralAngle::GetAtom2();
//MolAtom& MolDihedralAngle::GetAtom3();
//MolAtom& MolDihedralAngle::GetAtom4();
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
//      MolRing
//
//######################################################################
MolRing::MolRing(vector<MolBond*> &vRingBond):
mvpBond(vRingBond)
{}

const vector<MolBond*>& MolRing::GetBondList()const
{return mvpBond;}
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
   VFN_DEBUG_MESSAGE("Quaternion::~RotationQuaternion()",5)
}

Quaternion Quaternion::RotationQuaternion(const REAL ang,
                                          const REAL v1,
                                          const REAL v2,
                                          const REAL v3)
{
   VFN_DEBUG_MESSAGE("Quaternion::RotationQuaternion()",4)
   const REAL norm=sqrt(v1*v1+v2*v2+v3*v3);
   return Quaternion(cos(ang/2.),
                     sin(ang/2.)*v1/norm,
                     sin(ang/2.)*v2/norm,
                     sin(ang/2.)*v3/norm,
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
      ss <<mQ0;
      tag.AddAttribute("Q0",ss.str());
   }
   {
      stringstream ss;
      ss <<mQ1;
      tag.AddAttribute("Q1",ss.str());
   }
   {
      stringstream ss;
      ss <<mQ2;
      tag.AddAttribute("Q2",ss.str());
   }
   {
      stringstream ss;
      ss <<mQ3;
      tag.AddAttribute("Q3",ss.str());
   }
   {
      stringstream ss;
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

void Quaternion::Normalize()
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
//      Molecule
//
//######################################################################
Molecule::Molecule(Crystal &cryst, const string &name):
mIsSelfOptimizing(false)
{
   VFN_DEBUG_MESSAGE("Molecule::Molecule()",5)
   this->SetName(name);
   mpCryst=&cryst;
   {
      RefinablePar tmp(this->GetName()+"_x",&mXYZ(0),0.,1.,
                        gpRefParTypeScattTranslX,
                        REFPAR_DERIV_STEP_ABSOLUTE,false,false,true,true,1.,1.);
      tmp.AssignClock(mClockScatterer);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp(this->GetName()+"_y",&mXYZ(1),0,1,
                        gpRefParTypeScattTranslY,
                        REFPAR_DERIV_STEP_ABSOLUTE,false,false,true,true,1.,1.);
      tmp.AssignClock(mClockScatterer);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp(this->GetName()+"_z",&mXYZ(2),0,1,
                        gpRefParTypeScattTranslZ,
                        REFPAR_DERIV_STEP_ABSOLUTE,false,false,true,true,1.,1.);
      tmp.AssignClock(mClockScatterer);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp(this->GetName()+"_Occ",&mOccupancy,0,1,
                        gpRefParTypeScattOccup,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,1.,1.);
      tmp.AssignClock(mClockScatterer);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp(this->GetName()+"Q0",&(mQuat.Q0()),0,1,
                        gpRefParTypeScattOrient,
                        REFPAR_DERIV_STEP_ABSOLUTE,false,false,true,false,1.,1.);
      tmp.AssignClock(mClockScatterer);
      tmp.SetGlobalOptimStep(0.04);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp(this->GetName()+"Q1",&(mQuat.Q1()),0,1,
                        gpRefParTypeScattOrient,
                        REFPAR_DERIV_STEP_ABSOLUTE,false,false,true,false,1.,1.);
      tmp.AssignClock(mClockScatterer);
      tmp.SetGlobalOptimStep(0.04);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp(this->GetName()+"Q2",&(mQuat.Q2()),0,1,
                        gpRefParTypeScattOrient,
                        REFPAR_DERIV_STEP_ABSOLUTE,false,false,true,false,1.,1.);
      tmp.AssignClock(mClockScatterer);
      tmp.SetGlobalOptimStep(0.04);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp(this->GetName()+"Q3",&(mQuat.Q3()),0,1,
                        gpRefParTypeScattOrient,
                        REFPAR_DERIV_STEP_ABSOLUTE,false,false,true,false,1.,1.);
      tmp.AssignClock(mClockScatterer);
      tmp.SetGlobalOptimStep(0.04);
      this->AddPar(tmp);
   }
   mLocalParamSet=this->CreateParamSet("saved parameters for local minimization");
   this->InitOptions();
   mClockScatterer.AddChild(mClockAtomList);
   mClockScatterer.AddChild(mClockBondList);
   mClockScatterer.AddChild(mClockBondAngleList);
   mClockScatterer.AddChild(mClockDihedralAngleList);
   mClockScatterer.AddChild(mClockRingList);
   mClockScatterer.AddChild(mClockAtomPosition);
   mClockScatterer.AddChild(mClockAtomScattPow);
   mClockScatterer.AddChild(mClockOrientation);
}

Molecule::Molecule(const Molecule &old):
mIsSelfOptimizing(false)
{
   VFN_DEBUG_ENTRY("Molecule::Molecule(old&)",5)
   // a hack, but const-correct
   mpCryst=&(gCrystalRegistry.GetObj(gCrystalRegistry.Find(old.GetCrystal())));
   {
      RefinablePar tmp(this->GetName()+"_x",&mXYZ(0),0.,1.,
                        gpRefParTypeScattTranslX,
                        REFPAR_DERIV_STEP_ABSOLUTE,false,false,true,true,1.,1.);
      tmp.AssignClock(mClockScatterer);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp(this->GetName()+"_y",&mXYZ(1),0,1,
                        gpRefParTypeScattTranslY,
                        REFPAR_DERIV_STEP_ABSOLUTE,false,false,true,true,1.,1.);
      tmp.AssignClock(mClockScatterer);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp(this->GetName()+"_z",&mXYZ(2),0,1,
                        gpRefParTypeScattTranslZ,
                        REFPAR_DERIV_STEP_ABSOLUTE,false,false,true,true,1.,1.);
      tmp.AssignClock(mClockScatterer);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp(this->GetName()+"_Occ",&mOccupancy,0,1,
                        gpRefParTypeScattOccup,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,1.,1.);
      tmp.AssignClock(mClockScatterer);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp(this->GetName()+"Q0",&(mQuat.Q0()),0,1,
                        gpRefParTypeScattOrient,
                        REFPAR_DERIV_STEP_ABSOLUTE,false,false,true,false,1.,1.);
      tmp.AssignClock(mClockScatterer);
      tmp.SetGlobalOptimStep(0.04);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp(this->GetName()+"Q1",&(mQuat.Q1()),0,1,
                        gpRefParTypeScattOrient,
                        REFPAR_DERIV_STEP_ABSOLUTE,false,false,true,false,1.,1.);
      tmp.AssignClock(mClockScatterer);
      tmp.SetGlobalOptimStep(0.04);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp(this->GetName()+"Q2",&(mQuat.Q2()),0,1,
                        gpRefParTypeScattOrient,
                        REFPAR_DERIV_STEP_ABSOLUTE,false,false,true,false,1.,1.);
      tmp.AssignClock(mClockScatterer);
      tmp.SetGlobalOptimStep(0.04);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp(this->GetName()+"Q3",&(mQuat.Q3()),0,1,
                        gpRefParTypeScattOrient,
                        REFPAR_DERIV_STEP_ABSOLUTE,false,false,true,false,1.,1.);
      tmp.AssignClock(mClockScatterer);
      tmp.SetGlobalOptimStep(0.04);
      this->AddPar(tmp);
   }
   mLocalParamSet=this->CreateParamSet("saved parameters for local minimization");
   this->InitOptions();
   mClockScatterer.AddChild(mClockAtomList);
   mClockScatterer.AddChild(mClockBondList);
   mClockScatterer.AddChild(mClockBondAngleList);
   mClockScatterer.AddChild(mClockDihedralAngleList);
   mClockScatterer.AddChild(mClockRingList);
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

void Molecule::Print()const
{
   VFN_DEBUG_MESSAGE("Molecule::Print()",5)
   this->XMLOutput(cout);
}

void Molecule::XMLOutput(ostream &os,int indent)const
{
   VFN_DEBUG_ENTRY("Molecule::XMLOutput()",4)
   for(int i=0;i<indent;i++) os << "  " ;
   XMLCrystTag tag("Molecule");
   tag.AddAttribute("Name",mName);
   os <<tag<<endl;
   indent++;
   
   mQuat.XMLOutput(os,indent);
   
   this->GetPar(mXYZ.data()+0).XMLOutput(os,"x",indent);
   os <<endl;
   
   this->GetPar(mXYZ.data()+1).XMLOutput(os,"y",indent);
   os <<endl;
   
   this->GetPar(mXYZ.data()+2).XMLOutput(os,"z",indent);
   os <<endl;
   
   this->GetPar(&mOccupancy).XMLOutput(os,"Occup",indent);
   os <<endl;
   
   for(unsigned int i=0;i<this->GetNbOption();i++)
   {
      this->GetOption(i).XMLOutput(os,indent);
      os <<endl<<endl;
   }
   
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
      if("Option"==tag.GetName())
      {
         for(unsigned int i=0;i<tag.GetNbAttribute();i++)
            if("Name"==tag.GetAttributeName(i)) 
               mOptionRegistry.GetObj(tag.GetAttributeValue(i)).XMLInput(is,tag);
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

void Molecule::BeginOptimization(const bool allowApproximations,const bool enableRestraints)
{
   #if 1 // Is doing this automatically too dangerous ?
   if((!mIsSelfOptimizing) &&(this->GetLogLikelihood()>(mvpRestraint.size()*500)))
   {
      (*fpObjCrystInformUser)("Optimizing initial conformation of Molecule:"+this->GetName());
      this->OptimizeConformation(100000,(REAL)(mvpRestraint.size()));
      (*fpObjCrystInformUser)("");
   }
   #endif
   if(!mIsSelfOptimizing)
   {
      this->BuildRotorGroup();
      this->BuildFlipGroup();
   }
   this->RefinableObj::BeginOptimization(allowApproximations,enableRestraints);
   mRandomConformChangeNbTest=0;
   mRandomConformChangeNbTest=0;
   mRandomConformChangeTemp=1.;//(REAL)this->GetNbComponent();
}

void Molecule::RandomizeConfiguration()
{
   VFN_DEBUG_ENTRY("Molecule::RandomizeConfiguration()",4)
   if((!mIsSelfOptimizing) &&(this->GetLogLikelihood()>(mvpRestraint.size()*50)))
   {
      (*fpObjCrystInformUser)("Optimizing initial conformation of Molecule:"+this->GetName());
      this->OptimizeConformation(100000,(REAL)(mvpRestraint.size()));
      (*fpObjCrystInformUser)("");
   }

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
   // this will only change limited parameters i.e. translation
   this->RefinableObj::RandomizeConfiguration();
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
      VFN_DEBUG_EXIT("Molecule::GlobalOptRandomMove()",4)
      return;
   }
   TAU_PROFILE("Molecule::GlobalOptRandomMove()","void (REAL,RefParType*)",TAU_DEFAULT);
   TAU_PROFILE_TIMER(timer1,"Molecule::GlobalOptRandomMove 1","", TAU_FIELD);
   TAU_PROFILE_TIMER(timer2,"Molecule::GlobalOptRandomMove 2","", TAU_FIELD);
   TAU_PROFILE_TIMER(timer3,"Molecule::GlobalOptRandomMove 3","", TAU_FIELD);
   TAU_PROFILE_TIMER(timer4,"Molecule::GlobalOptRandomMove 4","", TAU_FIELD);
   TAU_PROFILE_START(timer1);
   VFN_DEBUG_ENTRY("Molecule::GlobalOptRandomMove()",4)
   mClockScatterer.Click();
   {//Rotate around an arbitrary vector
      static const REAL amp=2.*M_PI/100./RAND_MAX;
      mQuat *= Quaternion::RotationQuaternion
                  ((2.*(REAL)rand()-(REAL)RAND_MAX)*amp*mutationAmplitude,
                   (REAL)rand(),(REAL)rand(),(REAL)rand());
      mQuat.Normalize();
      mClockOrientation.Click();
   }
   //translation
   if(gpRefParTypeScattTransl->IsDescendantFromOrSameAs(type))
      this->RefinableObj::GlobalOptRandomMove(mutationAmplitude,gpRefParTypeScattTransl);
   TAU_PROFILE_STOP(timer1);
   if(gpRefParTypeScattConform->IsDescendantFromOrSameAs(type))
   {
      switch(mFlexModel.GetChoice())
      {
         case 0://Free atoms + restraints
         {
            TAU_PROFILE_START(timer2);
            this->SaveParamSet(mLocalParamSet);
            REAL lastll=this->GetLogLikelihood();
            //cout <<mutationAmplitude <<"oldLL="<<lastll<<" NewLL= ";
            unsigned int ct=0;
            TAU_PROFILE_STOP(timer2);
            while(true)
            {
               TAU_PROFILE_START(timer3);
               ct++;
               mRandomMoveIsDone=false;
               for(list<RotorGroup>::const_iterator pos=mvRotorGroupTorsion.begin();
                   pos!=mvRotorGroupTorsion.end();++pos)
               {//rotate around torsion bonds
                  if((rand()%2)==0) continue;
                  const REAL angle=(2.*(REAL)rand()-(REAL)RAND_MAX)
                                   *M_PI/25./(REAL)RAND_MAX*mutationAmplitude;
                  this->RotateAtomGroup(*(pos->mpAtom1),*(pos->mpAtom2),
                                        pos->mvRotatedAtomList,angle);
               }
               for(list<RotorGroup>::const_iterator pos=mvRotorGroupTorsionSingleChain.begin();
                   pos!=mvRotorGroupTorsionSingleChain.end();++pos)
               {//rotate around torsion bonds (single chain)
                  if((rand()%2)==0) continue;
                  const REAL angle=(2.*(REAL)rand()-(REAL)RAND_MAX)
                                   *M_PI/25./(REAL)RAND_MAX*mutationAmplitude;
                  this->RotateAtomGroup(*(pos->mpAtom1),*(pos->mpAtom2),
                                        pos->mvRotatedAtomList,angle);
               }
               bool doneFlip=false;
               list<FlipGroup>::const_iterator posFlip;
               if(((rand()%20)==0)&&(mvFlipGroup.size()>0))
               {// Try a flip from time to time
                  const unsigned long i=rand() % mvFlipGroup.size();
                  posFlip=mvFlipGroup.begin();
                  for(unsigned long j=0;j<i;++j)++posFlip;
                  #if 0
                  // If seems to break restraints, don't try it too often.
                  if(posFlip->mNbTest>100)
                     if(  (((REAL)(posFlip->mNbAccept)/(REAL)(posFlip->mNbTest))<0.1)
                        &&((rand()%10)!=0)) break;
                  posFlip->mNbTest++;
                  if((rand()%1000)==0)
                     for(list<FlipGroup>::const_iterator pos=mvFlipGroup.begin();
                         pos!=mvFlipGroup.end();++pos)
                     {
                        cout <<"Flip group with respect to: "
                             <<pos->mpAtom1->GetName()<<"-"
                             <<pos->mpAtom0->GetName()<<"-"
                             <<pos->mpAtom2->GetName()<<" : ";
                        for(list<pair<const MolAtom *,set<unsigned long> > >::const_iterator 
                            chain=pos->mvRotatedChainList.begin();
                            chain!=pos->mvRotatedChainList.end();++chain)
                        {
                           cout<<"    -"<<chain->first->GetName()<<":";
                           for(set<unsigned long>::iterator pos1=chain->second.begin();
                               pos1!=chain->second.end();++pos1)
                              cout<<mvpAtom[*pos1]->GetName()<<"  ";
                        }
                        cout<<"accept="<<pos->mNbAccept<<"/"<<pos->mNbTest<<endl;
                     }
                  #endif
                  this->FlipAtomGroup(*posFlip);
                  doneFlip=true;
               }
               TAU_PROFILE_STOP(timer3);
               TAU_PROFILE_START(timer4);
               this->RefinableObj::GlobalOptRandomMove(0.2,gpRefParTypeScattConform);
               //this->RefinableObj::GlobalOptRandomMove(mutationAmplitude,gpRefParTypeScattConform);
               mClockAtomPosition.Click();
               mRandomConformChangeNbTest++;
               mQuat.Normalize();
               const REAL newll=this->GetLogLikelihood();
               //cout <<newll<<" ";
               TAU_PROFILE_STOP(timer4);
               if(newll<lastll) break;
               if( log((rand()+1)/(REAL)RAND_MAX) 
                   < (-(newll-lastll)/mRandomConformChangeTemp ))
               {
                  if(doneFlip) posFlip->mNbAccept++;
                  break;
               }
               //if( log((rand()+1)/(REAL)RAND_MAX) 
               //    < (-(newll-lastll)/(mRandomConformChangeTemp*mutationAmplitude) )) break;
               this->RestoreParamSet(mLocalParamSet);
               mClockAtomPosition.Click();
               if(ct>20) break;
            }
            mRandomConformChangeNbAccept++;
            if(mRandomConformChangeNbTest>=1000)
            {
               if(((REAL)mRandomConformChangeNbAccept/(REAL)mRandomConformChangeNbTest)<0.20)
               {  
                  cout<<"mRandomConformChangeTemp="<<mRandomConformChangeTemp<<endl;
                  this->RestraintStatus(cout);
                  mRandomConformChangeTemp*=2.;
               }
               if(((REAL)mRandomConformChangeNbAccept/(REAL)mRandomConformChangeNbTest)<0.65)
               {  
                  cout<<"mRandomConformChangeTemp="<<mRandomConformChangeTemp<<endl;
                  this->RestraintStatus(cout);
                  mRandomConformChangeTemp*=1.3;
               }
               if(mRandomConformChangeTemp>0.1)
               {
                  if(((REAL)mRandomConformChangeNbAccept/(REAL)mRandomConformChangeNbTest)>0.75)
                  {
                     this->RestraintStatus(cout);
                     mRandomConformChangeTemp/=1.3;
                     cout<<"mRandomConformChangeTemp="<<mRandomConformChangeTemp<<endl;
                  }
                  if(((REAL)mRandomConformChangeNbAccept/(REAL)mRandomConformChangeNbTest)>0.90)
                  {
                     this->RestraintStatus(cout);
                     mRandomConformChangeTemp/=2.;
                     cout<<"mRandomConformChangeTemp="<<mRandomConformChangeTemp<<endl;
                  }
               }
               mRandomConformChangeNbTest=0;
               mRandomConformChangeNbAccept=0;
            }
            break;
         }
         case 1:break;//Rigid body
         case 2:// user-chosen torsion bonds
         {
            for(list<RotorGroup>::const_iterator pos=mvRotorGroupTorsion.begin();
                pos!=mvRotorGroupTorsion.end();++pos)
            {//rotate around torsion bonds
               const REAL angle=(2.*(REAL)rand()-(REAL)RAND_MAX)
                                *M_PI/50./(REAL)RAND_MAX*mutationAmplitude;
               this->RotateAtomGroup(*(pos->mpAtom1),*(pos->mpAtom2),
                                     pos->mvRotatedAtomList,angle);
            }
            break;
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
      &&(mClockLogLikelihood>mClockScatterer)) return mLogLikelihood;
   TAU_PROFILE("Molecule::GetLogLikelihood()","REAL ()",TAU_DEFAULT);
   mLogLikelihood=this->RefinableObj::GetLogLikelihood();
   mClockLogLikelihood.Click();
   return mLogLikelihood;
}
void Molecule::TagNewBestConfig()const
{
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
   if(mvpAtom[i]->IsDummy()) return "Dummy";
   return mvpAtom[i]->GetScatteringPower().GetName();
} 

ostream& Molecule::POVRayDescription(ostream &os,const bool noSymmetrics)const
{
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
   REAL en=1;
   if(displayEnantiomer==true) en=-1;
   this->UpdateScattCompList();
   
   const GLfloat colour_bondnonfree[]= { 0.3, .3, .3, 1.0 };
   const GLfloat colour_bondfree[]= { 0.7, .7, .7, 1.0 };
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
               if((r>0.8)&&(g>0.8)&&(b>0.8))
               {
                  colourChar[0] = 0.5;
                  colourChar[1] = 0.5;
                  colourChar[2] = 0.5;
               }
               glMaterialfv(GL_FRONT, GL_AMBIENT,   colour0); 
               glMaterialfv(GL_FRONT, GL_DIFFUSE,   colour0); 
               glMaterialfv(GL_FRONT, GL_SPECULAR,  colour0); 
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
               gluSphere(pQuadric,(*pos)->GetScatteringPower().GetRadius()/3.,10,10);
            }
         glPopMatrix();
      }
   }//Only independent atoms ?
   else
   {
      VFN_DEBUG_ENTRY("Molecule::GLInitDisplayList():Show all symmetrics",3)
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
         VFN_DEBUG_ENTRY("ZScatterer::GLInitDisplayList():Symmetric#"<<i,3)
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
            if(   (x(0)>xMin) && (x(0)<xMax)
                &&(y(0)>yMin) && (y(0)<yMax)
                &&(z(0)>zMin) && (z(0)<zMax))
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
                        gluSphere(pQuadric,
                           mvpAtom[k]->GetScatteringPower().GetRadius()/3.,10,10);
                     }
                  glPopMatrix();
               }
               if(displayNames==false)
               {
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
                     if(mvpBond[k]->IsFreeTorsion())
                        glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,colour_bondfree);
                     else
                        glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,colour_bondnonfree);
                     glMaterialfv(GL_FRONT, GL_SPECULAR,  colour0); 
                     glMaterialfv(GL_FRONT, GL_EMISSION,  colour0); 
                     glMaterialfv(GL_FRONT, GL_SHININESS, colour0);
                     glPolygonMode(GL_FRONT, GL_FILL);
                     glPushMatrix();
                        glTranslatef(x(n1)*en, y(n1), z(n1));
                        GLUquadricObj *quadobj = gluNewQuadric();
                        //glColor4f(1.0f,1.0f,1.0f,1.0);
                        const REAL height= sqrt(  (x(n2)-x(n1))*(x(n2)-x(n1))
                                                 +(y(n2)-y(n1))*(y(n2)-y(n1))
                                                 +(z(n2)-z(n1))*(z(n2)-z(n1)));
                        glRotatef(180,(x(n2)-x(n1))*en,y(n2)-y(n1),z(n2)-z(n1)+height);// ?!?!?!
                        gluCylinder(quadobj,.1,.1,height,10,1 );
                        gluDeleteQuadric(quadobj);
                     glPopMatrix();
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
   mvpAtom.push_back(new MolAtom(x,y,z,pPow,name,*this));
   mClockAtomPosition.Click();
   mClockAtomScattPow.Click();
   ++mScattCompList;
   {
      RefinablePar tmp(name+"_x",&(mvpAtom.back()->X()),0.,1.,
                        gpRefParTypeScattConformX,
                        REFPAR_DERIV_STEP_ABSOLUTE,false,false,true,false,1.,1.);
      tmp.AssignClock(mClockAtomPosition);
      tmp.SetGlobalOptimStep(0.05);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp(name+"_y",&(mvpAtom.back()->Y()),0.,1.,
                        gpRefParTypeScattConformY,
                        REFPAR_DERIV_STEP_ABSOLUTE,false,false,true,false,1.,1.);
      tmp.AssignClock(mClockAtomPosition);
      tmp.SetGlobalOptimStep(0.05);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp(name+"_z",&(mvpAtom.back()->Z()),0.,1.,
                        gpRefParTypeScattConformZ,
                        REFPAR_DERIV_STEP_ABSOLUTE,false,false,true,false,1.,1.);
      tmp.AssignClock(mClockAtomPosition);
      tmp.SetGlobalOptimStep(0.05);
      this->AddPar(tmp);
   }
   mClockScatterer.Click();
   if(updateDisplay) this->UpdateDisplay();
   VFN_DEBUG_EXIT("Molecule::AddAtom()",5)
}

vector<MolAtom*>::iterator Molecule::RemoveAtom(const MolAtom &atom)
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
   for(vector<MolBond*>::iterator posb=mvpBond.begin();posb!=mvpBond.end();++posb)
   {
      if( (&atom==&((*posb)->GetAtom1())) || (&atom==&((*posb)->GetAtom2())) )
      {
         posb=this->RemoveBond(**posb);
         --posb;
      }
   }
   for(vector<MolBondAngle*>::iterator posb=mvpBondAngle.begin();posb!=mvpBondAngle.end();++posb)
   {
      if(  (&atom==&((*posb)->GetAtom1())) || (&atom==&((*posb)->GetAtom2()))
         ||(&atom==&((*posb)->GetAtom3())))
      {
         posb=this->RemoveBondAngle(**posb);
         --posb;
      }
   }
   for(vector<MolDihedralAngle*>::iterator posb=mvpDihedralAngle.begin();
       posb!=mvpDihedralAngle.end();++posb)
   {
      if(  (&atom==&((*posb)->GetAtom1())) || (&atom==&((*posb)->GetAtom2()))
         ||(&atom==&((*posb)->GetAtom3())) || (&atom==&((*posb)->GetAtom4())))
      {
         posb=this->RemoveDihedralAngle(**posb);
         --posb;
      }
   }
   mClockAtomList.Click();
   mClockScatterer.Click();
   pos=mvpAtom.erase(pos);
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

vector<MolBond*>::iterator Molecule::RemoveBond(const MolBond &bond)
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

vector<MolBondAngle*>::iterator Molecule::RemoveBondAngle(const MolBondAngle &angle)
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

void Molecule::AddDihedralAngle(const MolAtom &atom1, const MolAtom &atom2,
                                const MolAtom &atom3, const MolAtom &atom4,
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

vector<MolDihedralAngle*>::iterator Molecule::RemoveDihedralAngle(const MolDihedralAngle &angle)
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

const vector<MolAtom*>& Molecule::GetAtomList()const{return mvpAtom;}
const vector<MolBond*>& Molecule::GetBondList()const{return mvpBond;}
const vector<MolBondAngle*>& Molecule::GetBondAngleList()const{return mvpBondAngle;}
const vector<MolDihedralAngle*>& Molecule::GetDihedralAngleList()const{return mvpDihedralAngle;}

vector<MolAtom*>& Molecule::GetAtomList(){return mvpAtom;}
vector<MolBond*>& Molecule::GetBondList(){return mvpBond;}
vector<MolBondAngle*>& Molecule::GetBondAngleList(){return mvpBondAngle;}
vector<MolDihedralAngle*>& Molecule::GetDihedralAngleList(){return mvpDihedralAngle;}

void Molecule::RotateAtomGroup(const MolAtom &at1,const MolAtom &at2,
                               const set<unsigned long> &atoms, const REAL angle)
{
   const REAL vx=at2.X()-at1.X();
   const REAL vy=at2.Y()-at1.Y();
   const REAL vz=at2.Z()-at1.Z();
   this->RotateAtomGroup(at1,vx,vy,vz,atoms,angle);
}
void Molecule::RotateAtomGroup(const MolAtom &at,const REAL vx,const REAL vy,const REAL vz,
                               const set<unsigned long> &atoms, const REAL angle)
{
   TAU_PROFILE("Molecule::RotateAtomGroup(MolAtom&,vx,vy,vz,...)","void (...)",TAU_DEFAULT);
   const REAL x0=at.X();
   const REAL y0=at.Y();
   const REAL z0=at.Z();
   // :KLUDGE: ? Refuse to do anything if vector is not well defined
   if((fabs(vx)+fabs(vy)+fabs(vz))<1e-6) return;
   const Quaternion quat=Quaternion::RotationQuaternion(angle,vx,vy,vz);
   for(set<unsigned long>::const_iterator pos=atoms.begin();pos!=atoms.end();++pos)
   {
      mvpAtom[*pos]->X() -= x0;
      mvpAtom[*pos]->Y() -= y0;
      mvpAtom[*pos]->Z() -= z0;
      quat.RotateVector(mvpAtom[*pos]->X(),mvpAtom[*pos]->Y(),mvpAtom[*pos]->Z());
      mvpAtom[*pos]->X() += x0;
      mvpAtom[*pos]->Y() += y0;
      mvpAtom[*pos]->Z() += z0;
   }
   mClockAtomPosition.Click();
   mClockScatterer.Click();
}

void Molecule::RestraintStatus(ostream &os)const
{
   VFN_DEBUG_ENTRY("Molecule::RestraintStatus()",5)
   for(vector<MolBond*>::const_iterator pos=mvpBond.begin();pos!=mvpBond.end();++pos)
      cout <<"Bond "<<(*pos)->GetName()
           <<"IdealLength="<<(*pos)->GetLength0()
           <<", Length="<<(*pos)->GetLength()
           <<", log(likelihood)="<<(*pos)->GetLogLikelihood()<<endl;
   for(vector<MolBondAngle*>::const_iterator pos=mvpBondAngle.begin();
       pos!=mvpBondAngle.end();++pos)
      cout <<"Bond Angle "<<(*pos)->GetName()
           <<"IdealAngle="<<(*pos)->Angle0()*180/M_PI
           <<", Angle="<<(*pos)->GetAngle()*180/M_PI
           <<", log(likelihood)="<<(*pos)->GetLogLikelihood()<<endl;
   for(vector<MolDihedralAngle*>::const_iterator pos=mvpDihedralAngle.begin();
       pos!=mvpDihedralAngle.end();++pos)
      cout <<"Dihedral Angle "<<(*pos)->GetName()
           <<"IdealAngle="<<(*pos)->Angle0()*180/M_PI
           <<", Angle="<<(*pos)->GetAngle()*180/M_PI
           <<", log(likelihood)="<<(*pos)->GetLogLikelihood()<<endl;
   VFN_DEBUG_EXIT("Molecule::RestraintStatus()",5)
}

const map<unsigned long,set<unsigned long> > &Molecule::GetConnectivityTable()const
{
   this->BuildConnectivityTable();
   return mConnectivityTable;
}

RefinableObjClock& Molecule::GetBondListClock(){return mClockBondList;}
const RefinableObjClock& Molecule::GetBondListClock()const{return mClockBondList;}

void Molecule::RigidifyWithDihedralAngles()
{
   this->BuildConnectivityTable();
   // Relationship between MolAtom adress and order
   map<const MolAtom*,unsigned long> index;
   {
      vector<MolAtom*>::const_iterator pos;
      unsigned long i=0;
      for(pos=mvpAtom.begin();pos<mvpAtom.end();++pos)
      {
         index[*pos]=i++;
      }
   }
   
   for(vector<MolBond*>::iterator bond=mvpBond.begin();bond!=mvpBond.end();++bond)
   {
      const MolAtom* at2=&((*bond)->GetAtom1());
      const MolAtom* at3=&((*bond)->GetAtom2());
      const unsigned long i2=index[at2];
      const unsigned long i3=index[at3];
      for(set<unsigned long>::const_iterator c2=mConnectivityTable[i2].begin();
          c2!=mConnectivityTable[i2].end();++c2)
      {
         MolAtom* at1=mvpAtom[*c2];
         if(at1==at3) continue;
         if(GetBondAngle(*at1,*at2,*at3)<(10 *DEG2RAD)) continue;
         if(GetBondAngle(*at1,*at2,*at3)>(180*DEG2RAD)) continue;
         for(set<unsigned long>::const_iterator c3=mConnectivityTable[i3].begin();
             c3!=mConnectivityTable[i3].end();++c3)
         {
            MolAtom* at4=mvpAtom[*c3];
            if((at4==at2)||(at4==at1)) continue;
            if(GetBondAngle(*at2,*at3,*at4)<(10 *DEG2RAD)) continue;
            if(GetBondAngle(*at2,*at3,*at4)>(180*DEG2RAD)) continue;
            if(this->FindDihedralAngle(*at1,*at2,*at3,*at4)==mvpDihedralAngle.end())
            {
               const REAL dihed=GetDihedralAngle(*at1,*at2,*at3,*at4);
               this->AddDihedralAngle(*at1,*at2,*at3,*at4,dihed,.01,.05,false);
            }
         }
      }
   }
   this->UpdateDisplay();
}

void Molecule::InitRefParList()
{
}

void Molecule::BuildRingList()const
{
   VFN_DEBUG_ENTRY("Molecule::BuildRingList()",5)
   VFN_DEBUG_EXIT("Molecule::BuildRingList()",5)
}

void Molecule::BuildConnectivityTable()const
{
   if(  (mClockConnectivityTable>mClockBondList)
      &&(mClockConnectivityTable>mClockAtomList)) return;
   VFN_DEBUG_ENTRY("Molecule::BuildConnectivityTable()",5)
   TAU_PROFILE("Molecule::BuildConnectivityTable()","void ()",TAU_DEFAULT);
   // Relationship between MolAtom adress and order
   map<const MolAtom*,unsigned long> index;
   {
      vector<MolAtom*>::const_iterator pos;
      unsigned long i=0;
      for(pos=mvpAtom.begin();pos<mvpAtom.end();++pos)
      {
         index[*pos]=i++;
      }
   }
   
   mConnectivityTable.clear();
   for(unsigned long i=0;i<mvpBond.size();++i)
   {
      mConnectivityTable[index[&(mvpBond[i]->GetAtom1())]]
                 .insert(index[&(mvpBond[i]->GetAtom2())]);
      mConnectivityTable[index[&(mvpBond[i]->GetAtom2())]]
                 .insert(index[&(mvpBond[i]->GetAtom1())]);
   }
   
   #if 1
   {
      map<unsigned long,set<unsigned long> >::const_iterator pos;
      unsigned long at=0;
      for(pos=mConnectivityTable.begin();pos!=mConnectivityTable.end();++pos)
      {
         cout<<"Atom "<<mvpAtom[at++]->GetName()<<" is connected to atoms: ";
         set<unsigned long>::const_iterator pos1;
         for(pos1=pos->second.begin();pos1!=pos->second.end();++pos1)
         {
            cout<<mvpAtom[*pos1]->GetName()<<"  ";
         }
         cout<<endl;
      }
   }
   #endif
   mClockConnectivityTable.Click();
   VFN_DEBUG_EXIT("Molecule::BuildConnectivityTable()",5)
}

/** Build recursively a list of atoms, starting from a one atom, and given
* a connectivity table.
*
* \param atom: the starting atom
* \param connect: the connectivity table
* \param atomlist: the list of atoms to which will be appended the atoms newly found.
* \param finalAtom: if specified, the list buildin will stop after finding this atom.
* This can be used to build the list of atoms between two given atoms. Otherwise,
* the list is expanded until the end of the chain(s), or until an atom already
* in the list is encountered (i.e. a ring has been found).
*/
void ExpandAtomGroupRecursive(const unsigned long atom,
                              const map<unsigned long,set<unsigned long> > &connect,
                              set<unsigned long> &atomlist,const long finalAtom=-1)
{
   const pair<set<unsigned long>::iterator,bool> status=atomlist.insert(atom);
   if(false==status.second) return;
   if(finalAtom==(long)atom) return;
   map<unsigned long,set<unsigned long> >::const_iterator c=connect.find(atom);
   set<unsigned long>::const_iterator pos;
   for(pos=c->second.begin();pos!=c->second.end();++pos)
   {
      ExpandAtomGroupRecursive(*pos,connect,atomlist,finalAtom);
   }
}

Molecule::RotorGroup::RotorGroup(const MolAtom &at1,const MolAtom &at2):
mpAtom1(&at1),mpAtom2(&at2)
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
   
   // Relationship between MolAtom adress and order
   map<const MolAtom*,unsigned long> index;
   {
      vector<MolAtom*>::const_iterator pos;
      unsigned long i=0;
      for(pos=mvpAtom.begin();pos<mvpAtom.end();++pos)
      {
         index[*pos]=i++;
      }
   }
   // Build Rotation groups around bonds
   for(unsigned long i=0;i<mvpBond.size();++i)
   {
      if((mFlexModel.GetChoice()!=0)&&(false==mvpBond[i]->IsFreeTorsion())) continue;
      const unsigned long atom1=index[&(mvpBond[i]->GetAtom1())];
      const unsigned long atom2=index[&(mvpBond[i]->GetAtom2())];
      for(unsigned int j=1;j<=2;++j)
      {
         const set<unsigned long> *pConn;
         if(j==1) pConn=&(mConnectivityTable[atom1]);
         else pConn=&(mConnectivityTable[atom2]);
         
         mvRotorGroupTorsion.push_back(RotorGroup(mvpBond[i]->GetAtom1(),
                                                  mvpBond[i]->GetAtom2()));
         mvRotorGroupTorsion.back().mvRotatedAtomList.insert(atom1);
         mvRotorGroupTorsion.back().mvRotatedAtomList.insert(atom2);
         
         for(set<unsigned long>::const_iterator pos=pConn->begin();pos!=pConn->end();++pos)
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
      for(unsigned long atom1=0;atom1<this->GetAtomList().size();++atom1)
      {
         const set<unsigned long> *pConn=&(mConnectivityTable[atom1]);
         for(unsigned long atom2=atom1+1;atom2<this->GetAtomList().size();++atom2)
         {
            for(set<unsigned long>::const_iterator pos=pConn->begin();pos!=pConn->end();++pos)
            {
               if(*pos==atom2) continue;
               mvRotorGroupInternal.push_back(RotorGroup(*(this->GetAtomList()[atom1]),
                                                         *(this->GetAtomList()[atom2])));
               mvRotorGroupInternal.back().mvRotatedAtomList.insert(atom1);
               ExpandAtomGroupRecursive(*pos,mConnectivityTable,
                                        mvRotorGroupInternal.back().mvRotatedAtomList,
                                        atom2);
               //Check if this chains leads to atom2
               set<unsigned long>::const_iterator check
                     =find(mvRotorGroupInternal.back().mvRotatedAtomList.begin(),
                           mvRotorGroupInternal.back().mvRotatedAtomList.end(),atom2);
               if(  (check==mvRotorGroupInternal.back().mvRotatedAtomList.end())
                  ||(mvRotorGroupInternal.back().mvRotatedAtomList.size()<3)
                  ||(mvRotorGroupInternal.back().mvRotatedAtomList.size()>=((mvpAtom.size()+1)/2)))
               {
                  mvRotorGroupInternal.pop_back();
               }
               else
               {
                  mvRotorGroupInternal.back().mvRotatedAtomList.erase(atom1);
                  mvRotorGroupInternal.back().mvRotatedAtomList.erase(atom2);
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
                pos2!=pRotorGroup2->end();++pos2)
            {
               if(pos2==pos1) continue;
               if((  ((pos1->mpAtom1 == pos2->mpAtom1) && (pos1->mpAtom2 == pos2->mpAtom2))
                   ||((pos1->mpAtom2 == pos2->mpAtom1) && (pos1->mpAtom1 == pos2->mpAtom2))) 
                  &&pos1->mvRotatedAtomList.size() == pos2->mvRotatedAtomList.size())
               {
                  bool ident=true;
                  for(set<unsigned long>::const_iterator pos=pos1->mvRotatedAtomList.begin();
                      pos!=pos1->mvRotatedAtomList.end();++pos)
                  {
                     set<unsigned long>::const_iterator tmp=pos2->mvRotatedAtomList.find(*pos);
                     if(tmp == pos2->mvRotatedAtomList.end())
                     {
                        ident=false;
                        break;
                     }
                  }
                  if(ident)
                  {
                     cout<<"Identical groups:"<<endl;
                     cout<<"    G1:"
                         <<pos1->mpAtom1->GetName()<<"-"
                         <<pos1->mpAtom2->GetName()<<" : ";
                     for(set<unsigned long>::iterator pos=pos1->mvRotatedAtomList.begin();
                         pos!=pos1->mvRotatedAtomList.end();++pos)
                        cout<<mvpAtom[*pos]->GetName()<<"  ";
                     cout<<endl;
                     cout<<"    G2:"
                         <<pos2->mpAtom1->GetName()<<"-"
                         <<pos2->mpAtom2->GetName()<<" : ";
                     for(set<unsigned long>::iterator pos=pos2->mvRotatedAtomList.begin();
                         pos!=pos2->mvRotatedAtomList.end();++pos)
                        cout<<mvpAtom[*pos]->GetName()<<"  ";
                     cout<<endl;
                     pos2=pRotorGroup2->erase(pos2);
                     --pos2;
                  }
               }
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
          pos!=pRotorGroup1->end();++pos)
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

         switch(i)
         {
            case 1: cout<<"Rotation Group around bond :";break;
            case 2: cout<<"Rotation Group (single chain) around bond :";break;
            case 3: cout<<"Rotation Group (internal) between :";break;
         }
         cout <<pos->mpAtom1->GetName()<<"-"
             <<pos->mpAtom2->GetName()<<" : ";
         for(set<unsigned long>::iterator pos1=pos->mvRotatedAtomList.begin();
             pos1!=pos->mvRotatedAtomList.end();++pos1)
            cout<<mvpAtom[*pos1]->GetName()<<"  ";
         cout<<"   <d(LLK)>="<< llk/36.;

         if((llk/50.)>100.)
         {
            pos = pRotorGroup1->erase(pos);
            --pos;
            cout <<" -> NOT a free torsion"<<endl;
         }
         else
            cout <<" -> free torsion"<<endl;
      }
   }
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

void Molecule::BuildFlipGroup()
{
   this->BuildConnectivityTable();
   if(mClockFlipGroup>mClockConnectivityTable) return;
   VFN_DEBUG_ENTRY("Molecule::BuildFlipGroup()",5)
   TAU_PROFILE("Molecule::BuildFlipGroup()","void ()",TAU_DEFAULT);
   mvFlipGroup.clear();
   
   // Relationship between MolAtom adress and order
   map<const MolAtom*,unsigned long> index;
   {
      vector<MolAtom*>::const_iterator pos;
      unsigned long i=0;
      for(pos=mvpAtom.begin();pos<mvpAtom.end();++pos)
      {
         index[*pos]=i++;
      }
   }
   
   for(unsigned long atom0=0;atom0<this->GetAtomList().size();++atom0)
   {
      const set<unsigned long> *pConn=&(mConnectivityTable[atom0]);
      if(pConn->size()<3) continue;
      // Build all chains
      for(set<unsigned long>::const_iterator pos1=pConn->begin();pos1!=pConn->end();++pos1)
      {
         for(set<unsigned long>::const_iterator pos2=pos1;pos2!=pConn->end();++pos2)
         {
            if(pos2==pos1) continue;
            if(mFlexModel.GetChoice()==0)
            {
               mvFlipGroup.push_back(FlipGroup(*mvpAtom[atom0],*mvpAtom[*pos1],*mvpAtom[*pos2]));
               bool foundRing=false;
               for(set<unsigned long>::const_iterator pos=pConn->begin();pos!=pConn->end();++pos)
               {
                  if((pos==pos1)||(pos==pos2)) continue;
                  mvFlipGroup.back().mvRotatedChainList.push_back(
                     make_pair(mvpAtom[*pos],set<unsigned long>()));
                  mvFlipGroup.back().mvRotatedChainList.back().second.insert(atom0);
                  ExpandAtomGroupRecursive(*pos,mConnectivityTable,
                                           mvFlipGroup.back().mvRotatedChainList.back().second);
                  mvFlipGroup.back().mvRotatedChainList.back().second.erase(atom0);
                  set<unsigned long>::const_iterator ringdetect1,ringdetect2;
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
               for(list<pair<const MolAtom *,set<unsigned long> > >::const_iterator 
                   chain=mvFlipGroup.back().mvRotatedChainList.begin();
                   chain!=mvFlipGroup.back().mvRotatedChainList.end();++chain)
                   flipSize+=chain->second.size();

               if(((flipSize*2)>mvpAtom.size())||foundRing) mvFlipGroup.pop_back();
            }
            // Add the entry which will exchange atom1 and atom2 (this entry can be a ring)
            mvFlipGroup.push_back(FlipGroup(*mvpAtom[atom0],*mvpAtom[*pos1],*mvpAtom[*pos2]));
            mvFlipGroup.back().mvRotatedChainList.push_back(
               make_pair(mvpAtom[atom0],set<unsigned long>()));
               mvFlipGroup.back().mvRotatedChainList.back().second.insert(atom0);
            ExpandAtomGroupRecursive(*pos1,mConnectivityTable,
                                     mvFlipGroup.back().mvRotatedChainList.back().second);
            ExpandAtomGroupRecursive(*pos2,mConnectivityTable,
                                     mvFlipGroup.back().mvRotatedChainList.back().second);
            mvFlipGroup.back().mvRotatedChainList.back().second.erase(atom0);
            if((mvFlipGroup.back().mvRotatedChainList.back().second.size()*2)>mvpAtom.size())
               mvFlipGroup.pop_back();
         }
      }
   }
   // List them
   this->SaveParamSet(mLocalParamSet);
   #if 0
   const REAL llk0=this->GetLogLikelihood();
   #endif
   for(list<FlipGroup>::iterator pos=mvFlipGroup.begin();
       pos!=mvFlipGroup.end();++pos)
   {
      if(pos->mvRotatedChainList.begin()->first==pos->mpAtom0)
      {
         cout <<"Flip group from atom "
              <<pos->mpAtom0->GetName()<<",exchanging bonds with "
              <<pos->mpAtom1->GetName()<<" and "
              <<pos->mpAtom2->GetName()<<", resulting in a 180° rotation of atoms : ";
         for(set<unsigned long>::iterator pos1=pos->mvRotatedChainList.begin()->second.begin();
             pos1!=pos->mvRotatedChainList.begin()->second.end();++pos1)
            cout<<mvpAtom[*pos1]->GetName()<<"  ";
      }
      else
      {
         cout <<"Flip group with respect to: "
              <<pos->mpAtom1->GetName()<<"-"
              <<pos->mpAtom0->GetName()<<"-"
              <<pos->mpAtom2->GetName()<<" : ";
         for(list<pair<const MolAtom *,set<unsigned long> > >::const_iterator 
             chain=pos->mvRotatedChainList.begin();
             chain!=pos->mvRotatedChainList.end();++chain)
         {
            cout<<"    -"<<chain->first->GetName()<<":";
            for(set<unsigned long>::const_iterator pos1=chain->second.begin();
                pos1!=chain->second.end();++pos1)
               cout<<mvpAtom[*pos1]->GetName()<<"  ";
         }
      }
      #if 0
      // test if they do not break something (dihedral angle restraint) ?
      // We seldom try flippping, so don't test
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
   mClockFlipGroup.Click();
   VFN_DEBUG_EXIT("Molecule::BuildFlipGroup()",5)
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
   // translate center to (0,0,0)
   {
      REAL x0=0,y0=0,z0=0;
      for(long i=0;i<nb;++i)
      {
         x0 += mScattCompList(i).mX;
         y0 += mScattCompList(i).mY;
         z0 += mScattCompList(i).mZ;
      }
      x0 /= nb;
      y0 /= nb;
      z0 /= nb;
      for(long i=0;i<nb;++i)
      {
         mScattCompList(i).mX -= x0;
         mScattCompList(i).mY -= y0;
         mScattCompList(i).mZ -= z0;
      }
   }
   //VFN_DEBUG_MESSAGE("Molecule::UpdateScattCompList()",10)
   // rotate
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
         VFN_DEBUG_EXIT("Molecule::FindAtom():"<<name<<"...found",4)
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
   VFN_DEBUG_ENTRY("Molecule::InitOptions",10)
   static string Flexname;
   static string Flexchoices[3];
   
   static bool needInitNames=true;
   if(true==needInitNames)
   {
      Flexname="Flexibility Model";
      Flexchoices[0]="Automatic from  Bond and Angle Restraints (recommended)";
      Flexchoices[1]="Rigid Body";
      Flexchoices[2]="User-Chosen Free Torsions";
      //Flexchoices[3]="Rigid body (relaxed)";
      //Flexchoices[4]="Torsion Angles (relaxed)";
      
      needInitNames=false;
   }
   mFlexModel.Init(3,&Flexname,Flexchoices);
   mFlexModel.SetChoice(0);
   //this->AddOption(&mFlexModel);
   
   VFN_DEBUG_EXIT("Molecule::InitOptions",10)
}

Molecule::FlipGroup::FlipGroup(const MolAtom &at0,const MolAtom &at1,const MolAtom &at2):
mpAtom0(&at0),mpAtom1(&at1),mpAtom2(&at2),mNbTest(0),mNbAccept(0)
{
}

void Molecule::FlipAtomGroup(const FlipGroup& group)
{
   TAU_PROFILE("Molecule::FlipAtomGroup(FlipGroup&)","void (...)",TAU_DEFAULT);
   if(group.mpAtom0==group.mvRotatedChainList.back().first)
   {// We are doing a 180° rotation exchanging two bonds
      const REAL vx=group.mpAtom0->X()-(group.mpAtom1->X()+group.mpAtom2->X())/2.;
      const REAL vy=group.mpAtom0->Y()-(group.mpAtom1->Y()+group.mpAtom2->Y())/2.;
      const REAL vz=group.mpAtom0->Z()-(group.mpAtom1->Z()+group.mpAtom2->Z())/2.;
      this->RotateAtomGroup(*(group.mpAtom0),vx,vy,vz,
                            group.mvRotatedChainList.back().second,M_PI);
   }
   else
   {// we are flipping bonds with respect to a plane defined by other bonds
      REAL v01x=group.mpAtom1->X()-group.mpAtom0->X();
      REAL v01y=group.mpAtom1->Y()-group.mpAtom0->Y();
      REAL v01z=group.mpAtom1->Z()-group.mpAtom0->Z();
      const REAL norm01=sqrt(v01x*v01x+v01y*v01y+v01z*v01z);
      v01x /= norm01;v01y /= norm01;v01z /= norm01;

      REAL v02x=group.mpAtom2->X()-group.mpAtom0->X();
      REAL v02y=group.mpAtom2->Y()-group.mpAtom0->Y();
      REAL v02z=group.mpAtom2->Z()-group.mpAtom0->Z();
      const REAL norm02=sqrt(v02x*v02x+v02y*v02y+v02z*v02z);
      v02x /= norm02;v02y /= norm02;v02z /= norm02;
      
      REAL v12x=group.mpAtom2->X()-group.mpAtom1->X();
      REAL v12y=group.mpAtom2->Y()-group.mpAtom1->Y();
      REAL v12z=group.mpAtom2->Z()-group.mpAtom1->Z();
      const REAL norm12=sqrt(v12x*v12x+v12y*v12y+v12z*v12z);
      v12x /= norm12;v12y /= norm12;v12z /= norm12;
      
      REAL v0mx=group.mpAtom0->X()-(group.mpAtom1->X()+group.mpAtom2->X())/2.;
      REAL v0my=group.mpAtom0->Y()-(group.mpAtom1->Y()+group.mpAtom2->Y())/2.;
      REAL v0mz=group.mpAtom0->Z()-(group.mpAtom1->Z()+group.mpAtom2->Z())/2.;
      const REAL norm0m=sqrt(v0mx*v0mx+v0my*v0my+v0mz*v0mz);
      v0mx /= norm0m;v0my /= norm0m;v0mz /= norm0m;

      if(fabs(v01x*v02x+v01y*v02y+v01z*v02z)
         >0.05*sqrt( (v01x*v01x+v01y*v01y+v01z*v01z)
                    *(v02x*v02x+v02y*v02y+v02z*v02z)))
      {
         REAL v012x=v01y*v02z-v01z*v02y;
         REAL v012y=v01z*v02x-v01x*v02z;
         REAL v012z=v01x*v02y-v01y*v02x;
         const REAL norm012=sqrt(v012x*v012x+v012y*v012y+v012z*v012z);
         v012x /= norm012;v012y /= norm012;v012z /= norm012;
         
      
         for(list<pair<const MolAtom *,set<unsigned long> > >::const_iterator
             chain=group.mvRotatedChainList.begin();
             chain!=group.mvRotatedChainList.end();++chain)
         {
            REAL v03x=chain->first->X()-group.mpAtom0->X();
            REAL v03y=chain->first->Y()-group.mpAtom0->Y();
            REAL v03z=chain->first->Z()-group.mpAtom0->Z();
            const REAL norm03=sqrt( v03x*v03x + v03y*v03y + v03z*v03z );
            v03x /= norm03;v03y /= norm03;v03z /= norm03;

            const REAL a1=v012x*v03x+v012y*v03y+v012z*v03z;
            const REAL a2= v0mx*v03x+ v0my*v03y+ v0mz*v03z;
            const REAL a3= v12x*v03x+ v12y*v03y+ v12z*v03z;
            REAL angle = -a1/sqrt(1-a3*a3);
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
                                  chain->second,2*angle);
         }
      }
   }
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
