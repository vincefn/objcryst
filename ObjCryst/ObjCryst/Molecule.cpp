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
#include "ObjCryst/Molecule.h"

using namespace std;

namespace ObjCryst
{
//######################################################################
//
//      MolAtom
//
//######################################################################

MolAtom::MolAtom(const REAL x, const REAL y, const REAL z,
              const ScatteringPower *pPow):
mX(x),mY(y),mZ(z),mOccupancy(1.),mpScattPow(pPow)
{
}

MolAtom::~MolAtom(){}

REAL MolAtom::GetX()const{return mX;}
REAL MolAtom::GetY()const{return mY;}
REAL MolAtom::GetZ()const{return mZ;}
REAL MolAtom::GetOccupancy()const{return mOccupancy;}

void MolAtom::SetX(const REAL a){ mX=a;}
void MolAtom::SetY(const REAL a){ mY=a;}
void MolAtom::SetZ(const REAL a){ mZ=a;}
void MolAtom::SetOccupancy(const REAL a){ mOccupancy=a;}

const ScatteringPower& MolAtom::GetScatteringPower()const{return *mpScattPow;}
void MolAtom::SetScatteringPower(const ScatteringPower& pow){mpScattPow=&pow;}

void MolAtom::XMLOutput(ostream &os,int indent)const
{
   VFN_DEBUG_ENTRY("MolAtom::XMLOutput()",4)
   for(int i=0;i<indent;i++) os << "  " ;
   XMLCrystTag tag("Atom",false,true);
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
{}

//######################################################################
//
//      MolBond
//
//######################################################################
MolBond::MolBond(const MolAtom &atom1, const MolAtom &atom2,
     const REAL length0, const REAL sigma, const REAL delta,
     const REAL bondOrder):
mAtomPair(make_pair(&atom1,&atom2)),
mLengthIdeal(length0),mDelta(delta),mSigma(sigma),
mBondOrder(bondOrder),mIsInRing(false)
{}

MolBond::~MolBond()
{}

void MolBond::XMLOutput(ostream &os,int indent)const
{
   VFN_DEBUG_ENTRY("MolBond::XMLOutput()",4)
   for(int i=0;i<indent;i++) os << "  " ;
   XMLCrystTag tag("Bond",false,true);
   //#error "which atoms for this bond ?"
   {
      stringstream ss;
      ss <<mLengthIdeal;
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
   VFN_DEBUG_ENTRY("MolBond::XMLInput():",5)
   for(unsigned int i=0;i<tag.GetNbAttribute();i++)
   {
      if("Length"==tag.GetAttributeName(i))
      {
         stringstream ss(tag.GetAttributeValue(i));
         ss >>mLengthIdeal;
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
   VFN_DEBUG_EXIT("MolBond::XMLInput():",5)
}

REAL MolBond::GetLogLikelihood()const
{
   const REAL length=this->GetLength();
   REAL tmp=length-(mLengthIdeal+mDelta);
   if(tmp>0)
   {
      tmp/=mSigma;
      return tmp*tmp;
   }
   tmp=length-(mLengthIdeal-mDelta);
   if(tmp<0)
   {
      tmp/=mSigma;
      return tmp*tmp;
   }
   return 0;
}

const MolAtom& MolBond::GetAtom1()const{return *(mAtomPair.first);}
const MolAtom& MolBond::GetAtom2()const{return *(mAtomPair.second);}
REAL MolBond::GetLength()const
{
   return  sqrt( (this->GetAtom1().GetX()-this->GetAtom2().GetX())
                 *(this->GetAtom1().GetX()-this->GetAtom2().GetX())
                +(this->GetAtom1().GetY()-this->GetAtom2().GetY())
                 *(this->GetAtom1().GetY()-this->GetAtom2().GetY())
                +(this->GetAtom1().GetZ()-this->GetAtom2().GetZ())
                 *(this->GetAtom1().GetZ()-this->GetAtom2().GetZ()) );
}

REAL MolBond::GetIdealLength()const{return mLengthIdeal;}
REAL MolBond::GetLengthDelta()const{return mDelta;}
REAL MolBond::GetLengthSigma()const{return mSigma;}
REAL MolBond::GetBondOrder()const{return mBondOrder;}

void MolBond::SetIdealLength(const REAL a){mLengthIdeal=a;}
void MolBond::SetLengthDelta(const REAL a){mDelta=a;}
void MolBond::SetLengthSigma(const REAL a){mBondOrder=a;}
void MolBond::SetBondOrder(const REAL a){mSigma=a;}

bool MolBond::IsInRing()const{return mIsInRing;}
void MolBond::SetInRing(const bool isInRing){mIsInRing=isInRing;}
//######################################################################
//
//      MolBondAngle
//
//######################################################################
MolBondAngle::MolBondAngle(const MolAtom &atom1,const  MolAtom &atom2,const  MolAtom &atom3,
     const REAL angle, const REAL sigma, const REAL delta):
mAngleIdeal(angle),mDelta(delta),mSigma(sigma)
{
   mvpAtom.push_back(&atom1);
   mvpAtom.push_back(&atom2);
   mvpAtom.push_back(&atom3);
}

MolBondAngle::~MolBondAngle(){}

void MolBondAngle::XMLOutput(ostream &os,int indent)const
{
   VFN_DEBUG_ENTRY("MolBondAngle::XMLOutput()",4)
   for(int i=0;i<indent;i++) os << "  " ;
   XMLCrystTag tag("BondAngle",false,true);
   //#error "which atoms for this bond ?"
   {
      stringstream ss;
      ss <<mAngleIdeal;
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
   VFN_DEBUG_ENTRY("MolBondAngle::XMLInput():",5)
   for(unsigned int i=0;i<tag.GetNbAttribute();i++)
   {
      if("Angle"==tag.GetAttributeName(i))
      {
         stringstream ss(tag.GetAttributeValue(i));
         ss >>mAngleIdeal;
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
   VFN_DEBUG_EXIT("MolBondAngle::XMLInput():",5)
}

REAL MolBondAngle::GetAngle()const
{
   const REAL x21=this->GetAtom1().GetX()-this->GetAtom2().GetX();
   const REAL y21=this->GetAtom1().GetY()-this->GetAtom2().GetY();
   const REAL z21=this->GetAtom1().GetZ()-this->GetAtom2().GetZ();
   const REAL x23=this->GetAtom3().GetX()-this->GetAtom2().GetX();
   const REAL y23=this->GetAtom3().GetY()-this->GetAtom2().GetY();
   const REAL z23=this->GetAtom3().GetZ()-this->GetAtom2().GetZ();
   const REAL norm21=sqrt( (x21*x21+y21*y21+z21*z21)*(x21*x21+y21*y21+z21*z21));
   const REAL norm23=sqrt( (x23*x23+y23*y23+z23*z23)*(x23*x23+y23*y23+z23*z23));
   return acos( (x21*x23+y21*y23+z21*z23)/(norm21*norm23));
   
}

REAL MolBondAngle::GetLogLikelihood()const
{
   const REAL angle=this->GetAngle();
   REAL tmp=angle-(mAngleIdeal+mDelta);
   if(tmp>0)
   {
      tmp/=mSigma;
      return tmp*tmp;
   }
   tmp=angle-(mAngleIdeal-mDelta);
   if(tmp<0)
   {
      tmp/=mSigma;
      return tmp*tmp;
   }
   return 0;
}
const MolAtom& MolBondAngle::GetAtom1()const{return *(mvpAtom.at(0));}
const MolAtom& MolBondAngle::GetAtom2()const{return *(mvpAtom.at(1));}
const MolAtom& MolBondAngle::GetAtom3()const{return *(mvpAtom.at(2));}
//MolAtom& MolBondAngle::GetAtom1(){return *(mvpAtom.at(0));}
//MolAtom& MolBondAngle::GetAtom2(){return *(mvpAtom.at(1));}
//MolAtom& MolBondAngle::GetAtom3(){return *(mvpAtom.at(2));}
//######################################################################
//
//      MolDihedralAngle
//
//######################################################################
MolDihedralAngle::MolDihedralAngle(MolAtom &atom1, MolAtom &atom2, MolAtom &atom3, MolAtom &atom4,
     const REAL angle, const REAL sigma, const REAL delta):
mAngleIdeal(angle),mDelta(delta),mSigma(sigma)
{
   mvpAtom.push_back(&atom1);
   mvpAtom.push_back(&atom2);
   mvpAtom.push_back(&atom3);
   mvpAtom.push_back(&atom4);
}

MolDihedralAngle::~MolDihedralAngle(){}

void MolDihedralAngle::XMLOutput(ostream &os,int indent)const
{
   VFN_DEBUG_ENTRY("MolDihedralAngle::XMLOutput()",4)
   for(int i=0;i<indent;i++) os << "  " ;
   XMLCrystTag tag("DihedralAngle",false,true);
   //#error "which atoms for this bond ?"
   {
      stringstream ss;
      ss <<mAngleIdeal;
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
   for(unsigned int i=0;i<tag.GetNbAttribute();i++)
   {
      if("Angle"==tag.GetAttributeName(i))
      {
         stringstream ss(tag.GetAttributeValue(i));
         ss >>mAngleIdeal;
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
   const REAL norm21=sqrt( (x21*x21+y21*y21+z21*z21)*(x21*x21+y21*y21+z21*z21));
   const REAL norm34=sqrt( (x34*x34+y34*y34+z34*z34)*(x34*x34+y34*y34+z34*z34));
   const REAL angle=acos( (x21*x34+y21*y34+z21*z34)/(norm21*norm34));
   
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

REAL MolDihedralAngle::GetLogLikelihood()const
{
   const REAL angle=this->GetAngle();
   REAL tmp=angle-(mAngleIdeal+mDelta);
   if(tmp>0)
   {
      tmp/=mSigma;
      return tmp*tmp;
   }
   tmp=angle-(mAngleIdeal-mDelta);
   if(tmp<0)
   {
      tmp/=mSigma;
      return tmp*tmp;
   }
   return 0;
}

const MolAtom& MolDihedralAngle::GetAtom1()const{return *(mvpAtom.at(0));}
const MolAtom& MolDihedralAngle::GetAtom2()const{return *(mvpAtom.at(1));}
const MolAtom& MolDihedralAngle::GetAtom3()const{return *(mvpAtom.at(2));}
const MolAtom& MolDihedralAngle::GetAtom4()const{return *(mvpAtom.at(3));}
//MolAtom& MolDihedralAngle::GetAtom1();
//MolAtom& MolDihedralAngle::GetAtom2();
//MolAtom& MolDihedralAngle::GetAtom3();
//MolAtom& MolDihedralAngle::GetAtom4();
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
//      UnitQuaternion
//
//######################################################################
UnitQuaternion::UnitQuaternion():
mQ0(1),mQ1(0),mQ2(0),mQ3(0)
{
   VFN_DEBUG_MESSAGE("UnitQuaternion::UnitQuaternion()",5)
}

UnitQuaternion::UnitQuaternion(const REAL q0,
                               const REAL q1,
                               const REAL q2,
                               const REAL q3):
mQ0(q0),mQ1(q1),mQ2(q2),mQ3(q3)
{
   VFN_DEBUG_MESSAGE("UnitQuaternion::UnitQuaternion()",5)
   this->Normalize();
}

UnitQuaternion::~UnitQuaternion()
{
   VFN_DEBUG_MESSAGE("UnitQuaternion::~RotationQuaternion()",5)
}

UnitQuaternion UnitQuaternion::RotationQuaternion(const REAL ang,
                                                  const REAL v1,
                                                  const REAL v2,
                                                  const REAL v3)
{
   VFN_DEBUG_MESSAGE("UnitQuaternion::RotationQuaternion()",4)
   return UnitQuaternion(cos(ang),
                         sin(ang)*v1,
                         sin(ang)*v2,
                         sin(ang)*v3);
}

UnitQuaternion UnitQuaternion::GetConjugate()const
{
   return UnitQuaternion(mQ0,-mQ1,-mQ2,mQ3);
}
UnitQuaternion UnitQuaternion::operator*(const UnitQuaternion &q)const
{
   return UnitQuaternion
      (this->Q0()*q.Q0()-this->Q1()*q.Q1()-this->Q2()*q.Q2()-this->Q3()*q.Q3(),
       this->Q0()*q.Q1()+this->Q1()*q.Q0()+this->Q2()*q.Q3()-this->Q3()*q.Q2(),
       this->Q0()*q.Q2()+this->Q2()*q.Q0()-this->Q1()*q.Q3()+this->Q3()*q.Q1(),
       this->Q0()*q.Q3()+this->Q3()*q.Q0()+this->Q1()*q.Q2()-this->Q2()*q.Q1());
}

void UnitQuaternion::XMLOutput(ostream &os,int indent)const
{
   VFN_DEBUG_ENTRY("UnitQuaternion::XMLOutput()",4)
   for(int i=0;i<indent;i++) os << "  " ;
   XMLCrystTag tag("UnitQuaternion",false,true);
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
   os <<tag<<endl;
   VFN_DEBUG_EXIT("UnitQuaternion::XMLOutput()",4)
}

void UnitQuaternion::XMLInput(istream &is,const XMLCrystTag &tag)
{
   VFN_DEBUG_ENTRY("UnitQuaternion::XMLInput()",4)
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
   }
   VFN_DEBUG_EXIT("UnitQuaternion::XMLInput()",4)
}

void UnitQuaternion::RotateVector(REAL &v1,REAL &v2, REAL &v3)const
{
   UnitQuaternion P(0,v1,v2,v3);
   P= (*this)* P * this->GetConjugate();
   v1=P.Q1();
   v2=P.Q2();
   v3=P.Q3();
}

void UnitQuaternion::Normalize()
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

const REAL& UnitQuaternion::Q0()const{return mQ0;}
const REAL& UnitQuaternion::Q1()const{return mQ1;}
const REAL& UnitQuaternion::Q2()const{return mQ2;}
const REAL& UnitQuaternion::Q3()const{return mQ3;}
REAL& UnitQuaternion::Q0(){return mQ0;}
REAL& UnitQuaternion::Q1(){return mQ1;}
REAL& UnitQuaternion::Q2(){return mQ2;}
REAL& UnitQuaternion::Q3(){return mQ3;}
//######################################################################
//
//      Molecule
//
//######################################################################
Molecule::Molecule(Crystal &cryst, const string &name):
Scatterer()
{
   VFN_DEBUG_MESSAGE("Molecule::XMLOutput()",5)
   this->SetName(name);
   mpCryst=&cryst;
}

Molecule::Molecule(const Molecule &old)
{
   VFN_DEBUG_MESSAGE("Molecule::Molecule(old&)",5)
   //:TODO:
}

Molecule::~Molecule()
{
   VFN_DEBUG_MESSAGE("Molecule::~Molecule()",5)
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
   
   //this->GetPar(...).XMLOutput(os,"a",indent);
   //os <<endl;
   {
      vector<MolAtom>::const_iterator pos;
      for(pos=mvAtom.begin();pos!=mvAtom.end();++pos)
         pos->XMLOutput(os,indent);
   }
   {
      vector<MolBond>::const_iterator pos;
      for(pos=mvBond.begin();pos!=mvBond.end();++pos)
         pos->XMLOutput(os,indent);
   }
   {
      vector<MolBondAngle>::const_iterator pos;
      for(pos=mvBondAngle.begin();pos!=mvBondAngle.end();++pos)
         pos->XMLOutput(os,indent);
   }
   {
      vector<MolDihedralAngle>::const_iterator pos;
      for(pos=mvDihedralAngle.begin();pos!=mvDihedralAngle.end();++pos)
         pos->XMLOutput(os,indent);
   }
   
   indent--;
   tag.SetIsEndTag(true);
   for(int i=0;i<indent;i++) os << "  " ;
   os <<tag<<endl;
   VFN_DEBUG_EXIT("Molecule::XMLOutput()",4)
}

void Molecule::XMLInput(istream &is,const XMLCrystTag &tag)
{
   VFN_DEBUG_ENTRY("Molecule::XMLInput()",4)
   VFN_DEBUG_EXIT("Molecule::XMLInput()",4)
}

int Molecule::GetNbComponent() const { return mvAtom.size();}

const ScatteringComponentList& Molecule::GetScatteringComponentList() const
{
   VFN_DEBUG_ENTRY("Molecule::GetScatteringComponentList()",3)
   this->UpdateScattCompList();
   VFN_DEBUG_EXIT("Molecule::GetScatteringComponentList()",3)
   return mScattCompList;
}

string Molecule::GetComponentName(const int i) const
{
   return mvAtom.at(i).GetScatteringPower().GetName();
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
   VFN_DEBUG_ENTRY("Molecule::GLInitDisplayList()",3)
   VFN_DEBUG_EXIT("Molecule::GLInitDisplayList()",3)
}

void Molecule::AddAtom(const REAL x, const REAL y, const REAL z,
             const ScatteringPower *pPow)
{
   VFN_DEBUG_ENTRY("Molecule::AddAtom()",5)
   mvAtom.push_back(MolAtom(x,y,z,pPow));
   ++mScattCompList;
   VFN_DEBUG_EXIT("Molecule::AddAtom()",5)
}

void Molecule::AddBond(MolAtom &atom1, MolAtom &atom2,
             const REAL length, const REAL sigma, const REAL delta,
             const REAL bondOrder)
{
   VFN_DEBUG_ENTRY("Molecule::AddBond()",5)
   mvBond.push_back(MolBond(atom1,atom2,length,sigma,delta,bondOrder));
   VFN_DEBUG_EXIT("Molecule::AddBond()",5)
}

void Molecule::AddBondAngle(MolAtom &atom1, MolAtom &atom2, MolAtom &atom3,
                  const REAL angle, const REAL sigma, const REAL delta)
{
   VFN_DEBUG_ENTRY("Molecule::AddBondAngle()",5)
   mvBondAngle.push_back(MolBondAngle(atom1,atom2,atom3,angle,sigma,delta));
   VFN_DEBUG_EXIT("Molecule::AddBondAngle()",5)
}

void Molecule::AddDihedralAngle(MolAtom &atom1, MolAtom &atom2, MolAtom &atom3, MolAtom &atom4,
                      const REAL angle, const REAL sigma, const REAL delta)
{
   VFN_DEBUG_ENTRY("Molecule::AddDihedralAngle()",5)
   mvDihedralAngle.push_back(MolDihedralAngle(atom1,atom2,atom3,atom4,angle,sigma,delta));
   VFN_DEBUG_EXIT("Molecule::AddDihedralAngle()",5)
}

MolAtom &Molecule::GetAtom(unsigned int i){return mvAtom.at(i);}

const MolAtom &Molecule::GetAtom(unsigned int i)const{return mvAtom.at(i);}

void MinimizeConfiguration()
{
}

void Molecule::InitRefParList()
{
}

void Molecule::BuildRingList()const
{
   VFN_DEBUG_ENTRY("Molecule::BuildRingList()",5)
   VFN_DEBUG_EXIT("Molecule::BuildRingList()",5)
}

void Molecule::UpdateScattCompList()const
{
   if(  (mClockAtomPosition<mClockScattCompList)
      &&(mClockOrientation <mClockScattCompList)
      &&(mClockAtomScattPow<mClockScattCompList))return;
   VFN_DEBUG_ENTRY("Molecule::UpdateScattCompList()",5)
   // Get internal coords
   for(long i;i<this->GetNbComponent();++i)
   {
      mScattCompList(i).mpScattPow=&(mvAtom.at(i).GetScatteringPower());
      mScattCompList(i).mX=mvAtom.at(i).GetX();
      mScattCompList(i).mY=mvAtom.at(i).GetX();
      mScattCompList(i).mZ=mvAtom.at(i).GetX();
      mScattCompList(i).mOccupancy=mvAtom.at(i).GetOccupancy();
   }
   // translate center to (0,0,0)
   {
      REAL x0=0,y0=0,z0=0;
      for(long i;i<this->GetNbComponent();++i)
      {
         x0 += mScattCompList(i).mX;
         y0 += mScattCompList(i).mY;
         z0 += mScattCompList(i).mZ;
      }
      x0 /= this->GetNbComponent();
      y0 /= this->GetNbComponent();
      z0 /= this->GetNbComponent();
      for(long i;i<this->GetNbComponent();++i)
      {
         mScattCompList(i).mX -= x0;
         mScattCompList(i).mY -= y0;
         mScattCompList(i).mZ -= z0;
      }
   }
   // rotate
   for(long i;i<this->GetNbComponent();++i)
   {
      mQuat.RotateVector(mScattCompList(i).mX,mScattCompList(i).mY,mScattCompList(i).mZ);
   }
   
   // Convert to fractionnal coordinates
   mClockScattCompList.Click();
   VFN_DEBUG_EXIT("Molecule::UpdateScattCompList()",5)
}


}//namespace
