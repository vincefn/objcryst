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

#ifdef OBJCRYST_GL
#include <GL/glut.h>
#endif

#include "ObjCryst/Molecule.h"
#include "RefinableObj/GlobalOptimObj.h"
#ifdef __WX__CRYST__
   #include "wxCryst/wxMolecule.h"
#endif

using namespace std;

namespace ObjCryst
{
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
   tag.AddAttribute("ScattPow",this->GetScatteringPower().GetName());
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
      //tmp=abs(tmp)/mSigma;
      //if(tmp>30) return 2.8550185e25*(tmp-29);
      //tmp=sinh(tmp);
      VFN_DEBUG_EXIT("MolBond::GetLogLikelihood():",2)
      return tmp*tmp;
   }
   tmp=length-(mLength0-mDelta);
   if(tmp<0)
   {
      tmp/=mSigma;
      //tmp=abs(tmp)/mSigma;
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
   return  sqrt( (this->GetAtom1().GetX()-this->GetAtom2().GetX())
                 *(this->GetAtom1().GetX()-this->GetAtom2().GetX())
                +(this->GetAtom1().GetY()-this->GetAtom2().GetY())
                 *(this->GetAtom1().GetY()-this->GetAtom2().GetY())
                +(this->GetAtom1().GetZ()-this->GetAtom2().GetZ())
                 *(this->GetAtom1().GetZ()-this->GetAtom2().GetZ()) );
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
void MolBond::SetFreeTorsion(const bool isFreeTorsion){mIsFreeTorsion=isFreeTorsion;}
#ifdef __WX__CRYST__
WXCrystObjBasic* MolBond::WXCreate(wxWindow* parent)
{
   VFN_DEBUG_ENTRY("MolBond::WXCreate()",5)
   mpWXCrystObj=new WXMolBond(parent,this);
   VFN_DEBUG_EXIT("MolBond::WXCreate()",5)
   return mpWXCrystObj;
}
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
      ss <<mAngle0;
      tag.AddAttribute("Angle",ss.str());
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
   os <<tag<<endl;
   VFN_DEBUG_EXIT("MolBondAngle::XMLOutput()",4)
}

void MolBondAngle::XMLInput(istream &is,const XMLCrystTag &tag)
{
   VFN_DEBUG_ENTRY("MolBondAngle::XMLInput():",10)
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
   }
   VFN_DEBUG_EXIT("MolBondAngle::XMLInput():",10)
}
REAL& MolBondAngle::Angle0()
{
   return mAngle0;
}

REAL MolBondAngle::GetAngle()const
{
   const REAL x21=this->GetAtom1().GetX()-this->GetAtom2().GetX();
   const REAL y21=this->GetAtom1().GetY()-this->GetAtom2().GetY();
   const REAL z21=this->GetAtom1().GetZ()-this->GetAtom2().GetZ();
   const REAL x23=this->GetAtom3().GetX()-this->GetAtom2().GetX();
   const REAL y23=this->GetAtom3().GetY()-this->GetAtom2().GetY();
   const REAL z23=this->GetAtom3().GetZ()-this->GetAtom2().GetZ();
   const REAL norm21= x21*x21+y21*y21+z21*z21;
   const REAL norm23= x23*x23+y23*y23+z23*z23;
   return acos( (x21*x23+y21*y23+z21*z23)/sqrt(norm21*norm23));
}

REAL MolBondAngle::GetLogLikelihood()const
{
   VFN_DEBUG_ENTRY("MolBondAngle::GetLogLikelihood():",2)
   const REAL angle=this->GetAngle();
   REAL tmp=angle-(mAngle0+mDelta);
   if(tmp>0)
   {
      tmp/=mSigma;
      //tmp=abs(tmp)/mSigma;
      //if(tmp>30) return 2.8550185e25*(tmp-29);
      //tmp=sinh(tmp);
      VFN_DEBUG_EXIT("MolBondAngle::GetLogLikelihood():",2)
      return tmp*tmp;
   }
   tmp=angle-(mAngle0-mDelta);
   if(tmp<0)
   {
      tmp/=mSigma;
      //tmp=abs(tmp)/mSigma;
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
   vector<const MolAtom*>::iterator pos;
   for(pos=mvpAtom.begin();pos!=mvpAtom.end();++pos) cout << (*pos)->GetName()<<endl;
   VFN_DEBUG_EXIT("MolDihedralAngle::MolDihedralAngle()",5)
}

MolDihedralAngle::~MolDihedralAngle(){}

const Molecule& MolDihedralAngle::GetMolecule()const{return *mpMol;}
      Molecule& MolDihedralAngle::GetMolecule()     {return *mpMol;}

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
      ss <<mAngle0;
      tag.AddAttribute("DihedralAngle",ss.str());
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
   }
   VFN_DEBUG_EXIT("MolDihedralAngle::XMLInput():",5)
}

REAL MolDihedralAngle::GetAngle()const
{
   const REAL x21=this->GetAtom1().GetX()-this->GetAtom2().GetX();
   const REAL y21=this->GetAtom1().GetY()-this->GetAtom2().GetY();
   const REAL z21=this->GetAtom1().GetZ()-this->GetAtom2().GetZ();
   const REAL x34=this->GetAtom4().GetX()-this->GetAtom3().GetX();
   const REAL y34=this->GetAtom4().GetY()-this->GetAtom3().GetY();
   const REAL z34=this->GetAtom4().GetZ()-this->GetAtom3().GetZ();
   const REAL norm21= x21*x21+y21*y21+z21*z21;
   const REAL norm34= x34*x34+y34*y34+z34*z34;
   const REAL angle=acos( (x21*x34+y21*y34+z21*z34)/sqrt(norm21*norm34));
   
   const REAL x23=this->GetAtom3().GetX()-this->GetAtom2().GetX();
   const REAL y23=this->GetAtom3().GetY()-this->GetAtom2().GetY();
   const REAL z23=this->GetAtom3().GetZ()-this->GetAtom2().GetZ();
   
   // v23 x v21
   const REAL x123= y23*z21-z23*y21;
   const REAL y123= z23*x21-x23*z21;
   const REAL z123= x23*y21-y23*x21;
   
   if((x123*x34 + y123*y34 + z123*z34)<0) return -angle;
   return angle;
}

REAL& MolDihedralAngle::Angle0(){return mAngle0;}

REAL MolDihedralAngle::GetLogLikelihood()const
{
   VFN_DEBUG_ENTRY("MolDihedralAngle::GetLogLikelihood():",2)
   const REAL angle=this->GetAngle();
   REAL tmp=angle-(mAngle0+mDelta);
   if(tmp>0)
   {
      tmp/=mSigma;
      //tmp=abs(tmp)/mSigma;
      //if(tmp>30) return 2.8550185e25*(tmp-29);
      //tmp=sinh(tmp);
      VFN_DEBUG_EXIT("MolDihedralAngle::GetLogLikelihood():",2)
      return tmp*tmp;
   }
   tmp=angle-(mAngle0-mDelta);
   if(tmp<0)
   {
      tmp/=mSigma;
      //tmp=abs(tmp)/mSigma;
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
   return Quaternion(cos(ang/2.),
                     sin(ang/2.)*v1,
                     sin(ang/2.)*v2,
                     sin(ang/2.)*v3,
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
   //#error P should not be a _UNIT_ quaternion...
   Quaternion P(0,v1,v2,v3,false);
   //cout<<"RotQuat:(n="<<this->GetNorm()<<")";this->XMLOutput(cout);
   //cout<<"before :(n="<<P.GetNorm()<<")";P.XMLOutput(cout);
   P= (*this)* P * this->GetConjugate();
   //cout<<"rotated:(n="<<P.GetNorm()<<")";P.XMLOutput(cout);
   v1=P.Q1();
   v2=P.Q2();
   v3=P.Q3();
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
      tmp.SetGlobalOptimStep(0.1);
      this->AddPar(tmp);
   }
   mLocalParamSet=this->CreateParamSet("saved parameters for local minimization");
   this->InitOptions();
}

Molecule::Molecule(const Molecule &old)
{
   VFN_DEBUG_MESSAGE("Molecule::Molecule(old&)",5)
   //:TODO:
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
         this->XMLOutput(cout);
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
         this->AddAtom(0.,0.,0.,(ScatteringPower *)0,"");
         mvpAtom.back()->XMLInput(is,tagg);
      }
      if("Bond"==tagg.GetName())
      {
         this->AddBond(this->GetAtom(0),this->GetAtom(1),1.5,.01,.05,1.);
         mvpBond.back()->XMLOutput(cout);
         mvpBond.back()->XMLInput(is,tagg);
         mvpBond.back()->XMLOutput(cout);
      }
      if("BondAngle"==tagg.GetName())
      {
         this->AddBondAngle(this->GetAtom(0),this->GetAtom(1),this->GetAtom(2),1.5,.01,.05);
         mvpBondAngle.back()->XMLInput(is,tagg);
      }
      if("DihedralAngle"==tagg.GetName())
      {
         this->AddDihedralAngle(this->GetAtom(0),this->GetAtom(1),
                                this->GetAtom(2),this->GetAtom(3),1.5,.01,.05);
         mvpDihedralAngle.back()->XMLInput(is,tagg);
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
   this->BuildTorsionAtomGroupTable();
   this->RefinableObj::BeginOptimization(allowApproximations,enableRestraints);
}

void Molecule::RandomizeConfiguration()
{
   VFN_DEBUG_ENTRY("Molecule::RandomizeConfiguration()",4)
   this->BuildTorsionAtomGroupTable();
   for(unsigned long torsion=0;torsion<mvpBond.size();torsion++)
   {
      if((rand()%2)==0)
      {
         list<set<unsigned long> >::const_iterator pos;
         for(pos =mTorsionAtomGroupTable[torsion].first.begin();
             pos!=mTorsionAtomGroupTable[torsion].first.end();
             ++pos)
         {
            this->RotateAtomGroup(mvpBond[torsion]->GetAtom1(),mvpBond[torsion]->GetAtom2(),
                                  *pos,(REAL)rand()*2.*M_PI/(REAL)RAND_MAX);
         }
      }
      else
      {
         list<set<unsigned long> >::const_iterator pos;
         for(pos =mTorsionAtomGroupTable[torsion].second.begin();
             pos!=mTorsionAtomGroupTable[torsion].second.end();
             ++pos)
         {
            this->RotateAtomGroup(mvpBond[torsion]->GetAtom1(),mvpBond[torsion]->GetAtom2(),
                                  *pos,(REAL)rand()*2.*M_PI/(REAL)RAND_MAX);
         }
      }
   }
   this->RefinableObj::RandomizeConfiguration();
   VFN_DEBUG_EXIT("Molecule::RandomizeConfiguration()",4)
}

void Molecule::GlobalOptRandomMove(const REAL mutationAmplitude,
                                       const RefParType *type)
{
   if(mRandomMoveIsDone) return;
   TAU_PROFILE("Molecule::GlobalOptRandomMove()","void (REAL,RefParType*)",TAU_DEFAULT);
   VFN_DEBUG_ENTRY("Molecule::GlobalOptRandomMove()",4)
   //:TODO: random moves using different models (free atoms, torsions, rigid body)
   //switch(mFlexModel.GetChoice())
   {
      //case 0://Free atoms + restraints
      {
      //#error faire tourner directement les coordonnées des atomes , sans garder comme membre un quaternion -> limiter le nombre de paramètres
         //if(rand()<(0.2*RAND_MAX))
         if(false)
         {//Rotate around an arbitrary vector
            static const REAL amp=M_PI/1000./RAND_MAX;
            mQuat *= Quaternion::RotationQuaternion((rand()-RAND_MAX/2)*amp*mutationAmplitude,
                                                        (REAL)rand(),(REAL)rand(),(REAL)rand());
            mClockOrientation.Click();
         }
         
         {//rotate around torsion bond
            const unsigned long torsion=rand() % mvpBond.size();
            const REAL angle=(rand()-RAND_MAX/2)*M_PI/50./RAND_MAX*mutationAmplitude;
            if((rand()%2)==0)
            {
               list<set<unsigned long> >::const_iterator pos;
               for(pos =mTorsionAtomGroupTable[torsion].first.begin();
                   pos!=mTorsionAtomGroupTable[torsion].first.end();
                   ++pos)
               {
                  this->RotateAtomGroup(mvpBond[torsion]->GetAtom1(),mvpBond[torsion]->GetAtom2(),
                                        *pos,angle);
               }
            }
            else
            {
               list<set<unsigned long> >::const_iterator pos;
               for(pos =mTorsionAtomGroupTable[torsion].second.begin();
                   pos!=mTorsionAtomGroupTable[torsion].second.end();
                   ++pos)
               {
                  this->RotateAtomGroup(mvpBond[torsion]->GetAtom1(),mvpBond[torsion]->GetAtom2(),
                                        *pos,angle);
               }
            }
         }
         #if 1
         if(true)
         {
            this->SaveParamSet(mLocalParamSet);
            const REAL lastll=this->GetLogLikelihood();
            //cout <<mutationAmplitude <<"oldLL="<<lastll<<" NewLL= ";
            unsigned int ct=0;
            while(true)
            {
               ct++;
               mRandomMoveIsDone=false;
               this->RefinableObj::GlobalOptRandomMove(mutationAmplitude,type);
               mQuat.Normalize();
               if(mIsSelfOptimizing) break;
               const REAL newll=this->GetLogLikelihood();
               //cout <<newll<<" ";
               if(newll<lastll) break;
               if( log((rand()+1)/(REAL)RAND_MAX) 
                   < (-(newll-lastll)/(200.*(0.+mutationAmplitude)*(1+mutationAmplitude)) )) break;
               this->RestoreParamSet(mLocalParamSet);
            }
            //cout <<endl;
            //if(ct>20) cout<<"Molecule::GlobalOptRandomMove:"<<mutationAmplitude<<", ct="<<ct<<endl;
            //break;
         }
         else this->RefinableObj::GlobalOptRandomMove(mutationAmplitude,type);
         #endif
         mClockScatterer.Click();
      }
   }
   #if 0
   /* :TODO: steepest/SA descent if too high ?*/
   REAL ll=this->GetLogLikelihood();
   const REAL maxll=this->GetNbComponent()*500.;
   this->SaveParamSet(mLocalParamSet);
   while(ll>maxll)
   {
      mLastLogLike=ll;
      mRandomMoveIsDone=false;
      this->RefinableObj::GlobalOptRandomMove(1.,type);
      ll=this->GetLogLikelihood();
      cout<<ll<<","<<mLastLogLike<<","<<maxll<<endl;
      if(ll>mLastLogLike)
      {
         if( log((rand()+1)/(REAL)RAND_MAX) < (-(ll-mLastLogLike)/100.) )
            this->SaveParamSet(mLocalParamSet);
         else
         {
            this->RestoreParamSet(mLocalParamSet);
            ll=mLastLogLike;
         }
      }
      else this->SaveParamSet(mLocalParamSet);
   }
   #endif
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
   #if 0
   const REAL nb=(REAL)this->GetNbComponent();
   mLogLikelihood /=nb;
   //:KLUDGE: avoid overflow
   if(mLogLikelihood>30) mLogLikelihood= 2.8550185e25*(mLogLikelihood-29)*nb;
   else mLogLikelihood= sinh(mLogLikelihood)*nb;
   #endif
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
                               const bool displayEnantiomer)const
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
   
   const GLfloat colour_bond[]= { 0.5, .5, .5, 1.0 };
   
   GLUquadricObj* pQuadric = gluNewQuadric();
   
   if(true==onlyIndependentAtoms)//
   {
      REAL xc=mXYZ(0),yc=mXYZ(1),zc=mXYZ(2);
      this->GetCrystal().FractionalToOrthonormalCoords(xc,yc,zc);
      vector<MolAtom*>::const_iterator pos;
      for(pos=mvpAtom.begin();pos!=mvpAtom.end();pos++)
      {
         
         if((*pos)->IsDummy())continue;
         glMaterialfv (GL_FRONT, GL_AMBIENT_AND_DIFFUSE,
                       (*pos)->GetScatteringPower().GetColourRGB());
         glPushMatrix();
            glTranslatef((*pos)->X()*en+xc, (*pos)->Y()+yc, (*pos)->Z()+zc);
            gluSphere(pQuadric,
               (*pos)->GetScatteringPower().GetRadius()/3.,10,10);
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
                  glMaterialfv (GL_FRONT,GL_AMBIENT_AND_DIFFUSE,
                     mvpAtom[k]->GetScatteringPower().GetColourRGB());
                  glPushMatrix();
                     glTranslatef(x(k)*en, y(k), z(k));
                     gluSphere(pQuadric,
                        mvpAtom[k]->GetScatteringPower().GetRadius()/3.,10,10);
                  glPopMatrix();
                  /*
                  glPushMatrix();
                     glRasterPos3f(x(k)*en+0.5, y(k), z(k));
                     glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18,*(mvpAtom[k]->GetName().c_str()+2));
                  glPopMatrix();
                  */
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
                  glMaterialfv (GL_FRONT,GL_AMBIENT_AND_DIFFUSE,colour_bond);
                  glPushMatrix();
                     glTranslatef(x(n1)*en, y(n1), z(n1));
                     GLUquadricObj *quadobj = gluNewQuadric();
                     glColor3f(1.0f,1.0f,1.0f);
                     const REAL height= sqrt(  (x(n2)-x(n1))*(x(n2)-x(n1))
                                              +(y(n2)-y(n1))*(y(n2)-y(n1))
                                              +(z(n2)-z(n1))*(z(n2)-z(n1)));
                     glRotatef(180,(x(n2)-x(n1))*en,y(n2)-y(n1),z(n2)-z(n1)+height);// ?!?!?!
                     gluCylinder(quadobj,.1,.1,height,10,1 );
                     gluDeleteQuadric(quadobj);
                  glPopMatrix();
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
             const ScatteringPower *pPow, const string &name)
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
   VFN_DEBUG_EXIT("Molecule::AddAtom()",5)
}

void Molecule::AddBond(MolAtom &atom1, MolAtom &atom2,
             const REAL length, const REAL sigma, const REAL delta,
             const REAL bondOrder)
{
   VFN_DEBUG_ENTRY("Molecule::AddBond()",5)
   mvpBond.push_back(new MolBond(atom1,atom2,length,sigma,delta,*this,bondOrder));
   this->AddRestraint(mvpBond.back());
   VFN_DEBUG_EXIT("Molecule::AddBond()",5)
}

void Molecule::AddBondAngle(MolAtom &atom1, MolAtom &atom2, MolAtom &atom3,
                  const REAL angle, const REAL sigma, const REAL delta)
{
   VFN_DEBUG_ENTRY("Molecule::AddBondAngle()",5)
   mvpBondAngle.push_back(new MolBondAngle(atom1,atom2,atom3,angle,sigma,delta,*this));
   this->AddRestraint(mvpBondAngle.back());
   VFN_DEBUG_EXIT("Molecule::AddBondAngle()",5)
}

void Molecule::AddDihedralAngle(MolAtom &atom1, MolAtom &atom2, MolAtom &atom3, MolAtom &atom4,
                      const REAL angle, const REAL sigma, const REAL delta)
{
   VFN_DEBUG_ENTRY("Molecule::AddDihedralAngle()",5)
   mvpDihedralAngle.push_back(new MolDihedralAngle(atom1,atom2,atom3,atom4,
                                                   angle,sigma,delta,*this));
   this->AddRestraint(mvpDihedralAngle.back());
   VFN_DEBUG_EXIT("Molecule::AddDihedralAngle()",5)
}

MolAtom &Molecule::GetAtom(unsigned int i){return *mvpAtom[i];}

const MolAtom &Molecule::GetAtom(unsigned int i)const{return *mvpAtom[i];}

MolAtom &Molecule::GetAtom(const string &name){return **(this->FindAtom(name));}

const MolAtom &Molecule::GetAtom(const string &name)const{return **(this->FindAtom(name));}

void Molecule::OptimizeConformation(const long nbTrial)
{
   VFN_DEBUG_ENTRY("Molecule::OptimizeConformation()",5)
   MonteCarloObj globalOptObj(true);
   globalOptObj.AddRefinableObj(*this);
   globalOptObj.SetAlgorithmParallTempering(ANNEALING_EXPONENTIAL,10000.,1.,
                                            ANNEALING_EXPONENTIAL,10,.1);      

   long nb=nbTrial;
   mIsSelfOptimizing=true;
   globalOptObj.Optimize(nb);
   mIsSelfOptimizing=false;
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
   TAU_PROFILE("Molecule::RotateAtomGroup()","void (...)",TAU_DEFAULT);
   const REAL x0=at1.X();
   const REAL y0=at1.Y();
   const REAL z0=at1.Z();
   const REAL vx=at2.X()-at1.X();
   const REAL vy=at2.Y()-at1.Y();
   const REAL vz=at2.Z()-at1.Z();
   const Quaternion quat=Quaternion::RotationQuaternion(angle,vx,vy,vz);
   //quat.XMLOutput(cout);
   set<unsigned long>::const_iterator pos;
   //cout<<"Center of rotation @"<<x0<<" "<<y0<<" "<<z0<<endl;
   for(pos=atoms.begin();pos!=atoms.end();++pos)
   {
      //cout<<"Rotating:"<<mvpAtom[*pos]->X()<<" "<<mvpAtom[*pos]->Y()<<" "<<mvpAtom[*pos]->Z()<<endl;
      mvpAtom[*pos]->X() -= x0;
      mvpAtom[*pos]->Y() -= y0;
      mvpAtom[*pos]->Z() -= z0;
      //cout<<"      ->"<<mvpAtom[*pos]->X()<<" "<<mvpAtom[*pos]->Y()<<" "<<mvpAtom[*pos]->Z()<<endl;
      quat.RotateVector(mvpAtom[*pos]->X(),mvpAtom[*pos]->Y(),mvpAtom[*pos]->Z());
      //cout<<"      ->"<<mvpAtom[*pos]->X()<<" "<<mvpAtom[*pos]->Y()<<" "<<mvpAtom[*pos]->Z()<<endl;
      mvpAtom[*pos]->X() += x0;
      mvpAtom[*pos]->Y() += y0;
      mvpAtom[*pos]->Z() += z0;
      //cout<<"      ->"<<mvpAtom[*pos]->X()<<" "<<mvpAtom[*pos]->Y()<<" "<<mvpAtom[*pos]->Z()<<endl;
   }
   mClockAtomPosition.Click();
   mClockScatterer.Click();
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

void ExpandAtomGroupRecursive(const unsigned long atom,
                              const map<unsigned long,set<unsigned long> > &connect,
                              set<unsigned long> &atomlist)
{
   const pair<set<unsigned long>::iterator,bool> status=atomlist.insert(atom);
   if(false==status.second) return;
   map<unsigned long,set<unsigned long> >::const_iterator c=connect.find(atom);
   set<unsigned long>::const_iterator pos;
   for(pos=c->second.begin();pos!=c->second.end();++pos)
   {
      ExpandAtomGroupRecursive(*pos,connect,atomlist);
   }
}

void Molecule::BuildTorsionAtomGroupTable()const
{
   if(  (mClockmTorsionAtomGroupTable>mClockBondList)
      &&(mClockmTorsionAtomGroupTable>mClockAtomList)) return;
   VFN_DEBUG_ENTRY("Molecule::BuildTorsionAtomGroupTable()",5)
   TAU_PROFILE("Molecule::BuildTorsionAtomGroupTable()","void ()",TAU_DEFAULT);
   this->BuildConnectivityTable();
   mTorsionAtomGroupTable.clear();
   
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

   for(unsigned long i=0;i<mvpBond.size();++i)
   {
      // First atom
      {
         const unsigned long atom=index[&(mvpBond[i]->GetAtom1())];
         list<set<unsigned long> > *pList=&(mTorsionAtomGroupTable[i].first);
         const set<unsigned long> *pConn=&(mConnectivityTable[atom]);
         set<unsigned long>::const_iterator pos;
         for(pos=pConn->begin();pos!=pConn->end();++pos)
         {
            if(*pos==index[&(mvpBond[i]->GetAtom2())]) continue;
            list<set<unsigned long> >::const_iterator pos1;
            bool foundring=false;
            for(pos1=pList->begin();pos1!=pList->end();++pos1)
            {
               if(pos1->find(*pos)!=pos1->end()) foundring=true;
            }
            if(foundring) continue;
            pList->push_back();
            pList->back().insert(atom);
            ExpandAtomGroupRecursive(*pos,mConnectivityTable,pList->back());
         }
      }
      // Second atom
      {
         const unsigned long atom=index[&(mvpBond[i]->GetAtom2())];
         list<set<unsigned long> > *pList=&(mTorsionAtomGroupTable[i].second);
         const set<unsigned long> *pConn=&(mConnectivityTable[atom]);
         set<unsigned long>::const_iterator pos;
         for(pos=pConn->begin();pos!=pConn->end();++pos)
         {
            if(*pos==index[&(mvpBond[i]->GetAtom1())]) continue;
            list<set<unsigned long> >::const_iterator pos1;
            bool foundring=false;
            for(pos1=pList->begin();pos1!=pList->end();++pos1)
            {
               if(pos1->find(*pos)!=pos1->end()) foundring=true;
            }
            if(foundring) continue;
            pList->push_back();
            pList->back().insert(atom);
            ExpandAtomGroupRecursive(*pos,mConnectivityTable,pList->back());
         }
      }
      #if 1
      {
         cout<<"Atom groups around bond:"
             <<mvpBond[i]->GetAtom1().GetName()<<"-"
             <<mvpBond[i]->GetAtom2().GetName()<<endl;
         list<set<unsigned long> >::const_iterator pos;
         for(pos=mTorsionAtomGroupTable[i].first.begin();
             pos!=mTorsionAtomGroupTable[i].first.end();++pos)
         {
            cout<<"   Atom1 group:";
            set<unsigned long>::const_iterator pos1;
            for(pos1=pos->begin();pos1!=pos->end();++pos1)
            {
               cout<<mvpAtom[*pos1]->GetName()<<"  ";
            }
            cout<<endl;
         }
         for(pos=mTorsionAtomGroupTable[i].second.begin();
             pos!=mTorsionAtomGroupTable[i].second.end();++pos)
         {
            cout<<"   Atom2 group:";
            set<unsigned long>::const_iterator pos1;
            for(pos1=pos->begin();pos1!=pos->end();++pos1)
            {
               cout<<mvpAtom[*pos1]->GetName()<<"  ";
            }
            cout<<endl;
         }
      }
      #endif
   }
   mClockmTorsionAtomGroupTable.Click();
   VFN_DEBUG_EXIT("Molecule::BuildTorsionAtomGroupTable()",5)
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
   static string Flexchoices[5];
   
   static bool needInitNames=true;
   if(true==needInitNames)
   {
      Flexname="Flexibility Model";
      Flexchoices[0]="Free Atoms & Restraints";
      Flexchoices[1]="Rigid Body";
      Flexchoices[2]="Rigid body (relaxed)";
      Flexchoices[3]="Torsion Angles";
      Flexchoices[4]="Torsion Angles (relaxed)";
      
      needInitNames=false;
   }
   mFlexModel.Init(5,&Flexname,Flexchoices);
   mFlexModel.SetChoice(0);
   this->AddOption(&mFlexModel);
   
   VFN_DEBUG_EXIT("Molecule::InitOptions",10)
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
