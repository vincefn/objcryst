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
/*   Atom.h
*  header file for the Atom scatterer
*
*/

#include <cmath>
#include <typeinfo>

#include <stdio.h> //for sprintf()

//#include "ObjCryst/ScatteringPower.h"
//#include "ObjCryst/Scatterer.h"
//#include "ObjCryst/Atom.h"
#include "ObjCryst/ZScatterer.h"
#include "ObjCryst/ScatteringData.h"

#include "Quirks/VFNStreamFormat.h" //simple formatting of integers, REALs..

#include "Quirks/VFNDebug.h"
#ifdef OBJCRYST_GL
#include <GL/glut.h>
#endif

#include <fstream>
#include <iomanip>

#ifdef __WX__CRYST__
   #include "wxCryst/wxZScatterer.h"
   #undef GetClassName // Conflict from wxMSW headers ? (cygwin)
#endif

namespace ObjCryst
{
//######################################################################
//
//  ZAtom   
//
//
//######################################################################
ZAtom::ZAtom(ZScatterer &scatt,const ScatteringPower *pow,
             const long atomBond, const REAL bondLength,
             const long atomAngle, const REAL bondAngle,
             const long atomDihedral, const REAL dihedralAngle,
             const REAL popu, const string &name):
mpScattPow(pow),
mAtomBond(atomBond),mAtomAngle(atomAngle),mAtomDihed(atomDihedral),
mBondLength(bondLength),mAngle(bondAngle),mDihed(dihedralAngle),
mOccupancy(popu),mName(name),mpScatt(&scatt)
{
   VFN_DEBUG_MESSAGE("ZAtom::ZAtom():("<<mName<<")",5)
}
ZAtom::~ZAtom()
{
   VFN_DEBUG_MESSAGE("ZAtom::~ZAtom():("<<mName<<")",5)
}
const string& ZAtom::GetClassName()const
{
   static string className="ZAtom";
   return className;
}
const string& ZAtom::GetName()const {return mName;}

ZScatterer& ZAtom::GetZScatterer(){return *mpScatt;}
const ZScatterer& ZAtom::GetZScatterer()const{return *mpScatt;}

void ZAtom::SetName(const string& name) {mName=name;}
long ZAtom::GetZBondAtom()const {return mAtomBond;}
long ZAtom::GetZAngleAtom()const {return mAtomAngle;}
long ZAtom::GetZDihedralAngleAtom()const {return mAtomDihed;}
REAL ZAtom::GetZBondLength()const {return mBondLength;}
REAL ZAtom::GetZAngle()const {return mAngle;}
REAL ZAtom::GetZDihedralAngle()const {return mDihed;}
REAL ZAtom::GetOccupancy()const {return mOccupancy;}
const ScatteringPower* ZAtom::GetScatteringPower()const{return mpScattPow;}
//:TODO: fix the following so that their clocks are clicked accordingly
void ZAtom::SetZBondLength(const REAL bond) {mBondLength=bond;}
void ZAtom::SetZAngle(const REAL angle) {mAngle=angle;}
void ZAtom::SetZDihedralAngle(const REAL dihed) {mDihed=dihed;}
void ZAtom::SetOccupancy(const REAL pop) {mOccupancy=pop;}
void ZAtom::SetScatteringPower(const ScatteringPower* scatt) {mpScattPow=scatt;}
#ifdef __WX__CRYST__
WXCrystObjBasic* ZAtom::WXCreate(wxWindow *parent)
{
   VFN_DEBUG_ENTRY("ZAtom::WXCreate()",7)
   mpWXCrystObj=new WXZAtom (parent,this);
   VFN_DEBUG_EXIT("ZAtom::WXCreate()",7)
   return mpWXCrystObj;
}
WXCrystObjBasic* ZAtom::WXGet()
{
   return mpWXCrystObj;
}
void ZAtom::WXDelete()
{
   if(0!=mpWXCrystObj)
   {
      VFN_DEBUG_MESSAGE("ZAtom::WXDelete()",5)
      delete mpWXCrystObj;
   }
   mpWXCrystObj=0;
}
void ZAtom::WXNotifyDelete()
{
   VFN_DEBUG_MESSAGE("ZAtom::WXNotifyDelete():"<<mName,10)
   mpWXCrystObj=0;
}
#endif
//######################################################################
//
//     ZMoveMinimizer
//
//
//######################################################################
ZMoveMinimizer::ZMoveMinimizer(ZScatterer &scatt):
RefinableObj(true),
mpZScatt(&scatt),mOptimObj(true)
{
   mOptimObj.SetAlgorithmSimulAnnealing(ANNEALING_EXPONENTIAL,.1,.001,
                                            ANNEALING_EXPONENTIAL,16,.25);
   mOptimObj.AddRefinableObj(*this);
   mOptimObj.AddCostFunction(*this,0);
}
ZMoveMinimizer::~ZMoveMinimizer(){}
unsigned int ZMoveMinimizer::GetNbCostFunction()const {return 1;}
const string& ZMoveMinimizer::GetCostFunctionName(const unsigned int)const
{ 
   static const string *str= new string("Conformation change");
   return *str;
}
const string& ZMoveMinimizer::GetCostFunctionDescription(const unsigned int)const
{ 
   static const string *str= new string("Conformation change distances squared");
   return *str;
}
REAL ZMoveMinimizer::GetCostFunctionValue(const unsigned int)
{
   TAU_PROFILE("ZMoveMinimizer::GetCostFunctionValue()","void ()",TAU_DEFAULT);
   const REAL *pX1=mpZScatt->GetXCoord().data();
   const REAL *pY1=mpZScatt->GetYCoord().data();
   const REAL *pZ1=mpZScatt->GetZCoord().data();
   const REAL *pX0=mXCoord0.data();
   const REAL *pY0=mYCoord0.data();
   const REAL *pZ0=mZCoord0.data();
   const REAL *pW=mAtomWeight.data();
   REAL dist=0;
   //for(int i=mXCoord0.numElements()-1;i>=0;i--)
   //   dist+= abs(*pX1++ - *pX0++) + abs(*pY1++ - *pY0++) + abs(*pZ1++ - *pZ0++);
   for(int i=mXCoord0.numElements()-1;i>=0;i--)
   {
      dist+= *pW++* ( (*pX1 - *pX0)*(*pX1 - *pX0) 
                     +(*pY1 - *pY0)*(*pY1 - *pY0)
                     +(*pZ1 - *pZ0)*(*pZ1 - *pZ0));
      pX1++;pY1++;pZ1++;
      pX0++;pY0++;pZ0++;
   }
   
   #if 0
   const CrystVector_REAL *pXcoord=&(mpZScatt->GetXCoord());
   const CrystVector_REAL *pYcoord=&(mpZScatt->GetYCoord());
   const CrystVector_REAL *pZcoord=&(mpZScatt->GetZCoord());
   REAL dist=0;
   for(int i=pXcoord->numElements()-1;i>=0;i--)
   {
      dist+=mAtomWeight(i)*( ((*pXcoord)(i)-mXCoord0(i))*((*pXcoord)(i)-mXCoord0(i))
                            +((*pYcoord)(i)-mYCoord0(i))*((*pYcoord)(i)-mYCoord0(i))
                            +((*pZcoord)(i)-mZCoord0(i))*((*pZcoord)(i)-mZCoord0(i)));
   }
   #endif
   return dist/mAtomWeight.sum();
}
void ZMoveMinimizer::RecordConformation()
{
   mXCoord0=mpZScatt->GetXCoord();
   mYCoord0=mpZScatt->GetYCoord();
   mZCoord0=mpZScatt->GetZCoord();
   if(mAtomWeight.numElements() != mXCoord0.numElements())
   {
      mAtomWeight.resize(mXCoord0.numElements());
      mAtomWeight=1;
   }
}
void ZMoveMinimizer::SetZAtomWeight(const CrystVector_REAL weight) {mAtomWeight=weight;}
void ZMoveMinimizer::MinimizeChange(long nbTrial)
{
   if(mAtomWeight.max()<1e-3) return;
   mOptimObj.Optimize(nbTrial,true);
}

//######################################################################
//
//  ZScatterer   
//
//
//######################################################################

ZScatterer::ZScatterer(const string &name,Crystal &cryst, 
                       const REAL x,const REAL y,const REAL z,
                       const REAL phi,const REAL chi, const REAL psi):
mScattCompList(0),mNbAtom(0),mNbDummyAtom(0),
mPhi(0),mChi(0),mPsi(0),
mZAtomRegistry("List of ZAtoms"),
mCenterAtomIndex(0),
mPhiChiPsiMatrix(3,3),
mUseGlobalScattPow(false),mpGlobalScattPow(0),
mpZMoveMinimizer(0)
{
   VFN_DEBUG_MESSAGE("ZScatterer::ZScatterer():("<<mName<<")",5)
   mName=name;
   mXYZ(0)=x;
   mXYZ(1)=y;
   mXYZ(2)=z;
   mPhi=phi;
   mChi=chi;
   mPsi=psi;
   this->SetCrystal(cryst);
   this->InitRefParList();
   VFN_DEBUG_MESSAGE("ZScatterer::ZScatterer():("<<mName<<")",5)
}

ZScatterer::ZScatterer(const ZScatterer &old):
Scatterer(old),m3DDisplayIndex(old.m3DDisplayIndex),
mScattCompList(old.mScattCompList),
mNbAtom(old.mNbAtom),mNbDummyAtom(old.mNbDummyAtom),
mPhi(old.mPhi),mChi(old.mChi),mPsi(old.mPsi),
mCenterAtomIndex(old.mCenterAtomIndex),
mPhiChiPsiMatrix(old.mPhiChiPsiMatrix),
mUseGlobalScattPow(false),mpGlobalScattPow(0),
mpZMoveMinimizer(0)
{
   VFN_DEBUG_MESSAGE("ZScatterer::ZScatterer(&old):("<<mName<<")",5)
   
   throw ObjCrystException("ZScatterer::ZScatterer(&old): the copy of all atoms has \
   not been implemented yet...");
   
   mName=old.GetName();
   mXYZ(0)=old.GetX();
   mXYZ(1)=old.GetY();
   mXYZ(2)=old.GetZ();
   mPhi=old.GetPhi();
   mChi=old.GetChi();
   mPsi=old.GetPsi();
   this->SetCrystal(*(old.mpCryst));
   
   this->InitRefParList();
   
   this->SetUseGlobalScatteringPower(old.mUseGlobalScattPow);
}

ZScatterer::~ZScatterer()
{
   VFN_DEBUG_MESSAGE("ZScatterer::~ZScatterer():("<<mName<<")",5)
   if(0 != mpGlobalScattPow) delete mpGlobalScattPow;
   mZAtomRegistry.DeleteAll();
}

ZScatterer* ZScatterer::CreateCopy() const
{
   VFN_DEBUG_MESSAGE("ZScatterer::CreateCopy()"<<mName<<")",5)
   return new ZScatterer(*this);
}
const string& ZScatterer::GetClassName() const
{
   const static string className="ZScatterer";
   return className;
}

void ZScatterer::AddAtom(const string &name,const ScatteringPower *pow,
             const long atomBond, const REAL bondLength,
             const long atomAngle, const REAL bondAngle,
             const long atomDihedral, const REAL dihedralAngle,
             const REAL popu)
{
   VFN_DEBUG_MESSAGE("ZScatterer::AddAtom():"<<name<<")",5)
   ZAtom *zatom =new ZAtom(*this,pow,
                           atomBond,bondLength,
                           atomAngle,bondAngle,
                           atomDihedral,dihedralAngle,
                           popu,name);
   
   if(true==mUseGlobalScattPow)
   {
      throw ObjCrystException("ZScatterer::AddAtom(() Cannot add an atom ! \
You are using the Global ScatteringPower approximation !!");

   }
   bool usedBond=true,usedAngle=true,usedDihed=true;
   if(mZAtomRegistry.GetNb()<1) usedBond=false;
   if(mZAtomRegistry.GetNb()<2) usedAngle=false;
   if(mZAtomRegistry.GetNb()<3) usedDihed=false;
   char buf [20];
   {
      sprintf(buf,"%d-%d",(int)mNbAtom,(int)atomBond);
      RefinablePar tmp("Length"+(string)buf,&(zatom->mBondLength),
                        1.,5.,
//                        bondLength*.9,bondLength*1.1,
                        gpRefParTypeScattConformBondLength,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,false,usedBond,false,1.);
      tmp.AssignClock(mClockScatterer);
      this->AddPar(tmp);
   }
   {
      sprintf(buf,"%d-%d-%d",(int)mNbAtom,(int)atomBond,(int)atomAngle);
      RefinablePar tmp("Angle"+(string)buf,&(zatom->mAngle),
                        0,2*M_PI,
//                        zatom->mAngle-.2,zatom->mAngle+.2,
                        gpRefParTypeScattConformBondAngle,
                        REFPAR_DERIV_STEP_ABSOLUTE,false,false,usedAngle,true,RAD2DEG,2*M_PI);
      tmp.AssignClock(mClockScatterer);
      this->AddPar(tmp);
   }
   {
      sprintf(buf,"%d-%d-%d-%d",(int)mNbAtom,(int)atomBond,(int)atomAngle,(int)atomDihedral);
      RefinablePar tmp("Dihed"+(string)buf,&(zatom->mDihed),
                        0,2*M_PI,
//                        zatom->mDihed-.2,zatom->mDihed+.2,
                        gpRefParTypeScattConformDihedAngle,
                        REFPAR_DERIV_STEP_ABSOLUTE,false,false,usedDihed,true,RAD2DEG,2*M_PI);
      tmp.AssignClock(mClockScatterer);
      this->AddPar(tmp);
   }
   {//fixed by default
      sprintf(buf,"%d",(int)mNbAtom);
      RefinablePar tmp("Occupancy"+(string)buf, 
                        &(zatom->mOccupancy),0,1,
                        gpRefParTypeScattOccup,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,1.,1.);
      tmp.AssignClock(mClockScatterer);
      this->AddPar(tmp);
   }
   
   //:NOTE: we *must* add it in the registry after declaring the parameters,
   // since it triggers the creation of the WXZAtom, which requires the parameters...
   VFN_DEBUG_MESSAGE("ZScatterer::AddAtom():Registering...",2)
   mZAtomRegistry.Register(*zatom);

   // Add an entry for this atom in the list of all components. The actual values are
   // written in Update().No entry for a dummy atom
   VFN_DEBUG_MESSAGE("ZScatterer::AddAtom():Adding to the ScattCompList...",2)
      if(0!=pow) 
      {
         ++mScattCompList;
         mScattCompList(mNbAtom-mNbDummyAtom).mpScattPow=pow;
         mComponentIndex.resizeAndPreserve(mNbAtom-mNbDummyAtom+1);
         mComponentIndex(mNbAtom-mNbDummyAtom)=mNbAtom;
      }
      else mNbDummyAtom++;
  
   //Finish
   cout << "Added atom:"<<mNbAtom
        << " : "<<this->GetZBondAtom(mNbAtom)
        << " : "<<this->GetZAngleAtom(mNbAtom)
        << " : "<<this->GetZDihedralAngleAtom(mNbAtom)
        << " : "<<this->GetZBondLength(mNbAtom)
        << " : "<<this->GetZAngle(mNbAtom)
        << " : "<<this->GetZDihedralAngle(mNbAtom)<<endl;
      mNbAtom++;
      mClockScatterer.Click();
   VFN_DEBUG_MESSAGE("ZScatterer::AddAtom():End",3)
}

int ZScatterer::GetNbComponent() const
{
   if(true==mUseGlobalScattPow) return 1;
   return mNbAtom-mNbDummyAtom;
}
const ScatteringComponentList& ZScatterer::GetScatteringComponentList() const
{
   this->UpdateScattCompList();
   
   return mScattCompList;
}

string ZScatterer::GetComponentName(const int i) const
{
   return mZAtomRegistry.GetObj(mComponentIndex(i)).GetName();
}

void ZScatterer::Print() const
{
   VFN_DEBUG_MESSAGE("ZScatterer::Print()",5)
   //:TODO: 
   //this->UpdateScattCompList();
   //for(int i=0;i<this->mNbAtom;i++) ??;
}

REAL ZScatterer::GetPhi()    const {return mPhi;}
REAL ZScatterer::GetChi()    const {return mChi;}
REAL ZScatterer::GetPsi()    const {return mPsi;}

void ZScatterer::SetPhi(const REAL x) { mClockScatterer.Click();mPhi=x;}
void ZScatterer::SetChi(const REAL y) { mClockScatterer.Click();mChi=y;}
void ZScatterer::SetPsi(const REAL z) { mClockScatterer.Click();mPsi=z;}

long ZScatterer::GetZBondAtom(const int i)const 
{return mZAtomRegistry.GetObj(i).GetZBondAtom();}

long ZScatterer::GetZAngleAtom(const int i)const 
{return mZAtomRegistry.GetObj(i).GetZAngleAtom();}

long ZScatterer::GetZDihedralAngleAtom(const int i)const
{return mZAtomRegistry.GetObj(i).GetZDihedralAngleAtom();}

REAL ZScatterer::GetZBondLength(const int i)const 
{return mZAtomRegistry.GetObj(i).GetZBondLength();}
REAL ZScatterer::GetZAngle(const int i)const
{return mZAtomRegistry.GetObj(i).GetZAngle();}
REAL ZScatterer::GetZDihedralAngle(const int i)const
{return mZAtomRegistry.GetObj(i).GetZDihedralAngle();}

void ZScatterer::SetZBondLength(const int i,const REAL a)
{mClockScatterer.Click();mZAtomRegistry.GetObj(i).SetZBondLength(a);}

void ZScatterer::SetZAngle(const int i,const REAL a)
{mClockScatterer.Click();mZAtomRegistry.GetObj(i).SetZAngle(a);}

void ZScatterer::SetZDihedralAngle(const int i,const REAL a)
   {mClockScatterer.Click();mZAtomRegistry.GetObj(i).SetZDihedralAngle(a);}

const ObjRegistry<ZAtom>& ZScatterer::GetZAtomRegistry()const
{return mZAtomRegistry;}

ostream& ZScatterer::POVRayDescription(ostream &os,
                                         const bool onlyIndependentAtoms)const
{
   VFN_DEBUG_MESSAGE("ZScatterer::POVRayDescription()",5)
   //throw ObjCrystException("ZScatterer::POVRayDescription() not implemented! "+this->GetName());
   //:TODO:
   this->UpdateScattCompList();

   //for(long i=0;i<mNbAtom;i++) mpAtom[i]->POVRayDescription(os,onlyIndependentAtoms);
   os << "// Description of ZScatterer :" << this->GetName() <<endl;
   
#if 0
   if(true==mUseGlobalScattPow) 
   {
      mpAtom[mCenterAtomIndex]->POVRayDescription(os,onlyIndependentAtoms);
      return os;
   }
   //Declare colours
   for(int i=0;i<mNbAtom;i++)
   {
      if(mpAtom[i]->IsDummy()) continue;
      os <<"   #declare colour_"<<mpAtom[i]->GetName()<<"="<<mpAtom[i]->Colour()<<";"<<endl;
   }
   
   if(true==onlyIndependentAtoms)
   {
      CrystVector_REAL x(mNbAtom),y(mNbAtom),z(mNbAtom);
      for(int i=0;i<mNbAtom;i++)
      {
         x(i)=mpAtom[i]->GetX();
         y(i)=mpAtom[i]->GetY();
         z(i)=mpAtom[i]->GetZ();
         this->GetCrystal().FractionalToOrthonormalCoords(x(i),y(i),z(i));
      }
      for(int i=0;i<mNbAtom;i++)
      {
         if(mpAtom[i]->IsDummy())
         {
            os << "   // Skipping Dummy Atom :" << mpAtom[i]->GetName() <<endl<<endl;
            continue;
         }
         os << "   // Atom :" << mpAtom[i]->GetName() <<endl;
         os << "      sphere " << endl;
         os << "      { <"<<x(i)<<","<<y(i)<<","<<z(i)<< ">,"
            << mpAtom[i]->GetRadius()/3<<endl;
         os << "          finish { ambient 0.2 diffuse 0.8 phong 1}" <<endl;
         os << "          pigment { colour colour_"<< mpAtom[i]->GetName() <<" }" << endl;
         os << "      }" <<endl;
         //Draw the bond for this Atom,if it's not linked to a dummy
         int bond=ZBondAtom(i);
         if((mpAtom[bond]->IsDummy()) || (i==0)) continue;
         os << "      cylinder"<<endl;
         os << "      { <"<<x(i)<<","<<y(i)<<","<<z(i)<< ">,"<<endl;
         os << "      <"<<x(bond)<<","<<y(bond)<<","<<z(bond)<< ">,"<<endl;
         os << "      0.1"<<endl;
         os << "      pigment { colour Gray }"<<endl;
         os << "      }"<<endl;
      }
   }
   else
   {
      CrystMatrix_REAL* xyzCoords=new CrystMatrix_REAL[mNbAtom];
      for(int i=0;i<mNbAtom;i++)
         *(xyzCoords+i)=this->GetCrystal().GetSpaceGroup().GetAllSymmetrics(mpAtom[i]->GetX(),
                                                            mpAtom[i]->GetY(),
                                                            mpAtom[i]->GetZ(),
                                                            false,false,false);
      const int nbSymmetrics=(xyzCoords+0)->rows();
      int symNum=0;
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
      CrystVector_REAL x(mNbAtom),y(mNbAtom),z(mNbAtom);
      CrystVector_REAL xSave,ySave,zSave;
      for(int i=0;i<nbSymmetrics;i++)
      {
         for(int j=0;j<mNbAtom;j++)
         {
            x(j)=(*(xyzCoords+j))(i,0);
            y(j)=(*(xyzCoords+j))(i,1);
            z(j)=(*(xyzCoords+j))(i,2);
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
            for(int j=1;j<mNbAtom;j++)
            {
               x(j) += dx;
               y(j) += dy;
               z(j) += dz;
            }
         //Generate also translated atoms near the unit cell
         const REAL limit =0.1;
         for(int j=0;j<translate.rows();j++)
         {
            xSave=x;
            ySave=y;
            zSave=z;
            x += translate(j,0);
            y += translate(j,1);
            z += translate(j,2);
            if(   (x(0)>(-limit)) && (x(0)<(1+limit))
                &&(y(0)>(-limit)) && (y(0)<(1+limit))
                &&(z(0)>(-limit)) && (z(0)<(1+limit)))
            {
               for(int k=0;k<mNbAtom;k++)
                  this->GetCrystal().FractionalToOrthonormalCoords(x(k),y(k),z(k));
               os << "  // Symmetric&Translated #" << symNum++ <<endl;
               //:NOTE: The code below is the same as without symmetrics
               for(int i=0;i<mNbAtom;i++)
               {
                  if(mpAtom[i]->IsDummy())
                  {
                     os << "   // Skipping Dummy Atom :" << mpAtom[i]->GetName() <<endl<<endl;
                     continue;
                  }
                  os << "   // Atom :" << mpAtom[i]->GetName() <<endl;
                  os << "      sphere " << endl;
                  os << "      { <"<<x(i)<<","<<y(i)<<","<<z(i)<< ">,"
                     << mpAtom[i]->GetRadius()/3<<endl;
                  os << "          finish { ambient 0.2 diffuse 0.8 phong 1}" <<endl;
                  os << "          pigment { colour colour_"<< mpAtom[i]->GetName() <<" }" << endl;
                  os << "      }" <<endl;
                  //Draw the bond for this Atom,if it's not linked to a dummy
                  int bond=ZBondAtom(i);
                  if((mpAtom[bond]->IsDummy()) || (i==0)) continue;
                  os << "      cylinder"<<endl;
                  os << "      { <"<<x(i)<<","<<y(i)<<","<<z(i)<< ">,"<<endl;
                  os << "      <"<<x(bond)<<","<<y(bond)<<","<<z(bond)<< ">,"<<endl;
                  os << "      0.1"<<endl;
                  os << "      pigment { colour Gray }"<<endl;
                  os << "      }"<<endl;
               }
            }//if in limits
            x=xSave;
            y=ySave;
            z=zSave;
         }//for translation
      }//for symmetrics
      delete[] xyzCoords;
   }//else
#endif
   return os;
}

void ZScatterer::GLInitDisplayList(const bool onlyIndependentAtoms,
                                   const REAL xMin,const REAL xMax,
                                   const REAL yMin,const REAL yMax,
                                   const REAL zMin,const REAL zMax,
                                   const bool displayEnantiomer)const
{
   #ifdef OBJCRYST_GL
   VFN_DEBUG_ENTRY("ZScatterer::GLInitDisplayList()",4)
   REAL en=1;
   if(displayEnantiomer==true) en=-1;
   this->UpdateScattCompList();
   if(true==mUseGlobalScattPow) 
   {
      //mpAtom[mCenterAtomIndex]->GLInitDisplayList(onlyIndependentAtoms);
      return;
   }
   
   GLfloat colour_bond[]= { 0.5, .5, .5, 1.0 };
   GLfloat colour_side[]= { 0.0, .0, .0, 1.0 };
   
   GLUquadricObj* pQuadric = gluNewQuadric();
   
   if(true==onlyIndependentAtoms)
   {
      //cout << m3DDisplayIndex <<endl;
      CrystVector_REAL x,y,z;
      x=mXCoord;
      y=mYCoord;
      z=mZCoord;
      if(m3DDisplayIndex.numElements()>0)
      {
         for(long k=0;k<m3DDisplayIndex.rows();k++)
         {
            switch(m3DDisplayIndex(k,0))
            {
               case 0:break;
               case 1: //Draw a sphere
               {
                  const int n1=m3DDisplayIndex(k,1);
                  if(0==mZAtomRegistry.GetObj(n1).GetScatteringPower())break;
                  glMaterialfv (GL_FRONT, GL_AMBIENT_AND_DIFFUSE,
                                mZAtomRegistry.GetObj(n1).GetScatteringPower()->GetColourRGB());
                  glPushMatrix();
                     glTranslatef(x(n1)*en, y(n1), z(n1));
                     //glutSolidSphere
                     gluSphere(pQuadric,mZAtomRegistry.GetObj(n1).GetScatteringPower()
                                       ->GetRadius()/3.,10,10);
                  glPopMatrix();
               }
               case 2: // Draw a bond
               {
                  long n1,n2;
                  n1=m3DDisplayIndex(k,1);
                  n2=m3DDisplayIndex(k,2);
                  glPushMatrix();
                     glTranslatef(x(n1)*en, y(n1), z(n1));
                     glMaterialfv (GL_FRONT, GL_AMBIENT_AND_DIFFUSE,colour_bond);
                     GLUquadricObj *quadobj = gluNewQuadric();
                     glColor3f(1.0f,1.0f,1.0f);
                     const REAL height= sqrt( (x(n2)-x(n1))*(x(n2)-x(n1))
                                              +(y(n2)-y(n1))*(y(n2)-y(n1))
                                              +(z(n2)-z(n1))*(z(n2)-z(n1)));
                     glRotatef(180,(x(n2)-x(n1))*en,y(n2)-y(n1),z(n2)-z(n1)+height);// ?!?!?!
                     gluCylinder(quadobj,.1,.1,height,10,1 );
                     gluDeleteQuadric(quadobj);
                  glPopMatrix();

               }
               case 3: // Draw a triangular face
               {
                  long n1,n2,n3;
                  REAL x1,y1,z1,x2,y2,z2,xn,yn,zn,xc,yc,zc;
                  n1=m3DDisplayIndex(k,1);
                  n2=m3DDisplayIndex(k,2);
                  n3=m3DDisplayIndex(k,3);
                  x1=x(n1)-x(n2);
                  y1=y(n1)-y(n2);
                  z1=z(n1)-z(n2);

                  x2=x(n1)-x(n3);
                  y2=y(n1)-y(n3);
                  z2=z(n1)-z(n3);

                  xn=y1*z2-z1*y2;
                  yn=z1*x2-x1*z2;
                  zn=x1*y2-y1*x2;

                  xc=(x(n1)+x(n2)+x(n3))/3.-x(mCenterAtomIndex);
                  yc=(y(n1)+y(n2)+y(n3))/3.-y(mCenterAtomIndex);
                  zc=(z(n1)+z(n2)+z(n3))/3.-z(mCenterAtomIndex);

                  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,
                                mZAtomRegistry.GetObj(0).GetScatteringPower()->GetColourRGB());
                  //glColor3f(1.0f,1.0f,1.0f);      // White
                  glBegin(GL_TRIANGLES);            // Bottom
                     if((xn*xc+yn*yc+zn*zc)>0) glNormal3f(xn*en, yn, zn);
                     else glNormal3f(-xn*en, -yn, -zn);
                     glVertex3f(x(n1)*en,y(n1),z(n1));
                     glVertex3f(x(n2)*en,y(n2),z(n2));
                     glVertex3f(x(n3)*en,y(n3),z(n3));
                  glEnd();
                  glMaterialfv (GL_FRONT, GL_AMBIENT_AND_DIFFUSE,colour_side);
                  glBegin(GL_LINE_LOOP);
                     glVertex3f(x(n1)*en,y(n1),z(n1));
                     glVertex3f(x(n2)*en,y(n2),z(n2));
                     glVertex3f(x(n3)*en,y(n3),z(n3));
                  glEnd();
               }
               case 4: // Draw a quadric face
               {
                  long n1,n2,n3,n4;
                  REAL x1,y1,z1,x2,y2,z2,xn,yn,zn,xc,yc,zc;
                  n1=m3DDisplayIndex(k,1);
                  n2=m3DDisplayIndex(k,2);
                  n3=m3DDisplayIndex(k,3);
                  n4=m3DDisplayIndex(k,3);
                  x1=x(n1)-x(n2);
                  y1=y(n1)-y(n2);
                  z1=z(n1)-z(n2);

                  x2=x(n1)-x(n3);
                  y2=y(n1)-y(n3);
                  z2=z(n1)-z(n3);

                  xn=y1*z2-z1*y2;
                  yn=z1*x2-x1*z2;
                  zn=x1*y2-y1*x2;

                  xc=(x(n1)+x(n2)+x(n3)+x(n4))/4.-x(mCenterAtomIndex);
                  yc=(y(n1)+y(n2)+y(n3)+y(n4))/4.-y(mCenterAtomIndex);
                  zc=(z(n1)+z(n2)+z(n3)+z(n4))/4.-z(mCenterAtomIndex);

                  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,
                                mZAtomRegistry.GetObj(0).GetScatteringPower()->GetColourRGB());
                  //glColor3f(1.0f,1.0f,1.0f);      // White
                  glBegin(GL_TRIANGLES);            // Bottom
                     if((xn*xc+yn*yc+zn*zc)>0) glNormal3f(xn*en, yn, zn);
                     else glNormal3f(-xn*en, -yn, -zn);
                     glVertex3f(x(n1)*en,y(n1),z(n1));
                     glVertex3f(x(n2)*en,y(n2),z(n2));
                     glVertex3f(x(n3)*en,y(n3),z(n3));
                     glVertex3f(x(n4)*en,y(n4),z(n4));
                  glEnd();
                  glMaterialfv (GL_FRONT, GL_AMBIENT_AND_DIFFUSE,colour_side);
                  glBegin(GL_LINE_LOOP);
                     glVertex3f(x(n1)*en,y(n1),z(n1));
                     glVertex3f(x(n2)*en,y(n2),z(n2));
                     glVertex3f(x(n3)*en,y(n3),z(n3));
                     glVertex3f(x(n4)*en,y(n4),z(n4));
                  glEnd();
               }
            }
         }
      }
      else
      {
         for(int k=0;k<mNbAtom;k++)
         {
            if(0==mZAtomRegistry.GetObj(k).GetScatteringPower())continue;
            glMaterialfv (GL_FRONT, GL_AMBIENT_AND_DIFFUSE,
                          mZAtomRegistry.GetObj(k).GetScatteringPower()->GetColourRGB());
            glPushMatrix();
               glTranslatef(x(k)*en, y(k), z(k));
               gluSphere(pQuadric,
                  mZAtomRegistry.GetObj(k).GetScatteringPower()->GetRadius()/3.,10,10);
               //Draw the bond for this Atom,if it's not linked to a dummy
               int bond=this->GetZBondAtom(k);
               if((0!=mZAtomRegistry.GetObj(bond).GetScatteringPower()) && (k>0))
               {
                  glMaterialfv (GL_FRONT, GL_AMBIENT_AND_DIFFUSE,colour_bond);
                  GLUquadricObj *quadobj = gluNewQuadric();
                  glColor3f(1.0f,1.0f,1.0f);
                  const REAL height= sqrt( (x(bond)-x(k))*(x(bond)-x(k))
                                           +(y(bond)-y(k))*(y(bond)-y(k))
                                           +(z(bond)-z(k))*(z(bond)-z(k)));
                  glRotatef(180,(x(bond)-x(k))*en,y(bond)-y(k),z(bond)-z(k)+height);// !!!
                  gluCylinder(quadobj,.1,.1,height,10,1 );
                  gluDeleteQuadric(quadobj);
               }
            glPopMatrix();
         }
      }//Use triangle faces ?
   }//Only independent atoms ?
   else
   {
      VFN_DEBUG_ENTRY("ZScatterer::GLInitDisplayList():Show all symmetrics",3)
      CrystMatrix_REAL xyzCoords[100]; //:TODO:
      {
         REAL x0,y0,z0;
         for(int i=0;i<mNbAtom;i++)
         {//We also generate the positions of dummy atoms.. This may be needed...
            x0=mXCoord(i);
            y0=mYCoord(i);
            z0=mZCoord(i);
            this->GetCrystal().OrthonormalToFractionalCoords(x0,y0,z0);
            xyzCoords[i]=this->GetCrystal().GetSpaceGroup().
                           GetAllSymmetrics(x0,y0,z0,false,false,false);
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
      CrystVector_REAL x(mNbAtom),y(mNbAtom),z(mNbAtom);
      CrystVector_REAL xSave,ySave,zSave;
      const int nbSymmetrics=xyzCoords[0].rows();
      for(int i=0;i<nbSymmetrics;i++)
      {
         VFN_DEBUG_ENTRY("ZScatterer::GLInitDisplayList():Symmetric#"<<i,3)
         for(int j=0;j<mNbAtom;j++)
         {
            x(j)=xyzCoords[j](i,0);
            y(j)=xyzCoords[j](i,1);
            z(j)=xyzCoords[j](i,2);
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
            for(int j=1;j<mNbAtom;j++)
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
               for(int k=0;k<mNbAtom;k++)
                  this->GetCrystal().FractionalToOrthonormalCoords(x(k),y(k),z(k));
               //:NOTE: The code below is the same as without symmetrics
               if(m3DDisplayIndex.numElements()>0)
               {
                  long n1,n2,n3;
                  REAL x1,y1,z1,x2,y2,z2,xn,yn,zn,xc,yc,zc;
                  for(long k=0;k<m3DDisplayIndex.rows();k++)
                  {
                     n1=m3DDisplayIndex(k,0);
                     n2=m3DDisplayIndex(k,1);
                     n3=m3DDisplayIndex(k,2);
                     x1=x(n1)-x(n2);
                     y1=y(n1)-y(n2);
                     z1=z(n1)-z(n2);

                     x2=x(n1)-x(n3);
                     y2=y(n1)-y(n3);
                     z2=z(n1)-z(n3);

                     xn=y1*z2-z1*y2;
                     yn=z1*x2-x1*z2;
                     zn=x1*y2-y1*x2;

                     xc=(x(n1)+x(n2)+x(n3))/3.-x(mCenterAtomIndex);
                     yc=(y(n1)+y(n2)+y(n3))/3.-y(mCenterAtomIndex);
                     zc=(z(n1)+z(n2)+z(n3))/3.-z(mCenterAtomIndex);

                     glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,
                                 mZAtomRegistry.GetObj(0).GetScatteringPower()->GetColourRGB());
                     glBegin(GL_TRIANGLES);
                        if((xn*xc+yn*yc+zn*zc)>0) glNormal3f( xn*en,  yn,  zn);
                        else                      glNormal3f(-xn*en, -yn, -zn);
                        
                        glVertex3f(x(n1)*en,y(n1),z(n1));
                        glVertex3f(x(n2)*en,y(n2),z(n2));
                        glVertex3f(x(n3)*en,y(n3),z(n3));
                     glEnd();
                     glMaterialfv (GL_FRONT, GL_AMBIENT_AND_DIFFUSE,colour_side);
                     glBegin(GL_LINE_LOOP);
                        glVertex3f(x(n1)*en,y(n1),z(n1));
                        glVertex3f(x(n2)*en,y(n2),z(n2));
                        glVertex3f(x(n3)*en,y(n3),z(n3));
                     glEnd();
                  }
               }
               else
               {
                  for(int k=0;k<mNbAtom;k++)
                  {
                     if(0==mZAtomRegistry.GetObj(k).GetScatteringPower())continue;
                     glMaterialfv (GL_FRONT,GL_AMBIENT_AND_DIFFUSE,
                        mZAtomRegistry.GetObj(k).GetScatteringPower()->GetColourRGB());
                     glPushMatrix();
                        glTranslatef(x(k)*en, y(k), z(k));
                        gluSphere(pQuadric,
                           mZAtomRegistry.GetObj(k).GetScatteringPower()->GetRadius()/3.,10,10);
                        //Draw the bond for this Atom,if it's not linked to a dummy
                        int bond=this->GetZBondAtom(k);
                        if((0!=mZAtomRegistry.GetObj(bond).GetScatteringPower()) && (k>0))
                        {
                           glMaterialfv (GL_FRONT, GL_AMBIENT_AND_DIFFUSE,colour_bond);
                           GLUquadricObj *quadobj = gluNewQuadric();
                           glColor3f(1.0f,1.0f,1.0f);
                           const REAL height= sqrt( (x(bond)-x(k))*(x(bond)-x(k))
                                                    +(y(bond)-y(k))*(y(bond)-y(k))
                                                    +(z(bond)-z(k))*(z(bond)-z(k)));
                           glRotatef(180,(x(bond)-x(k))*en,y(bond)-y(k),z(bond)-z(k)+height);// !!!
                           gluCylinder(quadobj,.1,.1,height,10,1 );
                           gluDeleteQuadric(quadobj);
                        }
                     glPopMatrix();
                  }
               }//Use triangle faces ?

               /*
               for(int k=0;k<1;k++)//mNbAtom
               {
                  if(0==mZAtomRegistry.GetObj(k).GetScatteringPower())continue;
                  glColor3f(0.8f,0.4f,0.4f);      // Red
                  glMaterialfv (GL_FRONT, GL_AMBIENT_AND_DIFFUSE,
                                mZAtomRegistry.GetObj(k).GetScatteringPower()->GetColourRGB());
                  glPushMatrix();
                     glTranslatef(x(k), y(k), z(k));
                     glutSolidSphere(
                        mZAtomRegistry.GetObj(k).GetScatteringPower()->GetRadius()/3.,10,10);
                     //Draw the bond for this Atom,if it's not linked to a dummy
                     const int bond=ZBondAtom(k);
                     if((0!=mZAtomRegistry.GetObj(bond).GetScatteringPower()->) && (k>0))
                     {
                        glMaterialfv (GL_FRONT, GL_AMBIENT_AND_DIFFUSE,colour_bond);
                        GLUquadricObj *quadobj = gluNewQuadric();
                        glColor3f(1.0f,1.0f,1.0f);      // Red
                        const REAL height= sqrt( (x(bond)-x(k))*(x(bond)-x(k))
                                                 +(y(bond)-y(k))*(y(bond)-y(k))
                                                 +(z(bond)-z(k))*(z(bond)-z(k)));
                        glRotatef(180,x(bond)-x(k),y(bond)-y(k),z(bond)-z(k)+height);// !!!
                        gluCylinder(quadobj,.1,.1,height,10,1 );
                        gluDeleteQuadric(quadobj);
                     }
                     glPopMatrix();
               }
               */
            }//if in limits
            x=xSave;
            y=ySave;
            z=zSave;
         }//for translation
      }//for symmetrics
      VFN_DEBUG_EXIT("ZScatterer::GLInitDisplayList():Show all symmetrics",3)
   }//else
   gluDeleteQuadric(pQuadric);
   VFN_DEBUG_EXIT("ZScatterer::GLInitDisplayList()",4)
   #endif //GLCryst
}

void ZScatterer::SetUseGlobalScatteringPower(const bool useIt)
{
   VFN_DEBUG_MESSAGE("ZScatterer::SetUseGlobalScatteringPower(bool):"<<this->GetName()<<":"<<useIt,5)
   if(true==useIt)
   {
      if(false==mUseGlobalScattPow)
      {
         mpGlobalScattPow=new GlobalScatteringPower(*this);
         mUseGlobalScattPow=true;
         mScattCompList.Reset();
         ++mScattCompList;
         this->InitRefParList();
         mClockScatterer.Click();
      }
      //Just change the functions which return the component list
   }
   else
   {
      if(true==mUseGlobalScattPow)
      {
         delete mpGlobalScattPow;
         mpGlobalScattPow=0;
         mUseGlobalScattPow=false;
         mScattCompList.Reset();
         for(long i=0;i<(mNbAtom-mNbDummyAtom);i++) ++mScattCompList;
         this->InitRefParList();
         mClockScatterer.Click();
      }
   }
}

void ZScatterer::GetGeneGroup(const RefinableObj &obj,
                                CrystVector_uint & groupIndex,
                                unsigned int &first) const
{
   // One group for all translation parameters, another for orientation.
   // All conformation parameters (bond, angle,torsion) are in individual groups.
   unsigned int posIndex=0;
   unsigned int orientIndex=0;
   VFN_DEBUG_MESSAGE("ZScatterer::GetGeneGroup()",4)
   for(long i=0;i<obj.GetNbPar();i++)
      for(long j=0;j<this->GetNbPar();j++)
         if(&(obj.GetPar(i)) == &(this->GetPar(j)))
         {
            if(this->GetPar(j).GetType()->IsDescendantFromOrSameAs(gpRefParTypeScattTransl))
            {
               if(posIndex==0) posIndex=first++;
               groupIndex(i)=posIndex;
            }
            else 
               if(this->GetPar(j).GetType()->IsDescendantFromOrSameAs(gpRefParTypeScattOrient))
               {
                  if(orientIndex==0) orientIndex=first++;
                  groupIndex(i)=orientIndex;
               }
               else groupIndex(i)= first++;
         }
}
const CrystVector_REAL& ZScatterer::GetXCoord() const
{
   this->UpdateCoordinates();
   return mXCoord;
}
const CrystVector_REAL& ZScatterer::GetYCoord() const
{
   this->UpdateCoordinates();
   return mYCoord;
}
const CrystVector_REAL& ZScatterer::GetZCoord() const
{
   this->UpdateCoordinates();
   return mZCoord;
}

void ZScatterer::EndOptimization()
{
   if(0!=mpZMoveMinimizer) delete mpZMoveMinimizer;
   mpZMoveMinimizer=0;
   this->RefinableObj::EndOptimization();
}
void ZScatterer::ImportFenskeHallZMatrix(istream &is)
{
   // Get read of "KEYWORD GO HERE", just in case...
   {
      const char c=is.peek();
      if ( (c < '0') || (c > '9') )
      {
         cout<<"ZScatterer::ImportFenskeHallZMatrix()"
             <<":getting rid of first line..."<<endl;
         char buf[100];
         is.getline(buf,100);
      }
   }
   // 17
   //C  1
   //N   1 1.465
   //C   2 1.366  1 119.987
   //N   3 1.321  2 120.030  1   6.0
   //C   4 1.355  3 119.982  2   6.8
   //N   5 1.136  4 180.000  3  46.3
   //N   3 1.366  2 120.022  1 186.0
   //C   7 1.466  3 119.988  2 354.9
   // ...
   int nbAtoms=0;
   is >> nbAtoms;
   string symbol;
   int bondAtom=0,angleAtom=0,dihedAtom=0,junk=0;
   float bond=0,angle=0,dihed=0;
   int scattPow;
   char buf [10];
   //first
      is >> symbol >> junk;
      cout << symbol <<" "
           << bondAtom  <<" "<< bond <<" "
           << angleAtom <<" "<< angle<<" "
           << dihedAtom <<" "<< dihed<<endl;
   {
      scattPow=mpCryst->GetScatteringPowerRegistry().Find
                  (symbol,"ScatteringPowerAtom",true);
      if(scattPow==-1)
      {
         cout<<"Scattering power"<<symbol<<"not found, creating it..."<<endl;
         mpCryst->AddScatteringPower(new ScatteringPowerAtom(symbol,symbol));
      }
      scattPow=mpCryst->GetScatteringPowerRegistry().Find
                  (symbol,"ScatteringPowerAtom");
      sprintf(buf,"%d",1);
      this->AddAtom(symbol+(string)buf,
                    &(mpCryst->GetScatteringPowerRegistry().GetObj(scattPow)),
                    0,0,
                    0,0,
                    0,0);
   }
   //second
      is >> symbol 
         >> bondAtom  >> bond;
   cout << symbol <<" "
        << bondAtom  <<" "<< bond <<" "
        << angleAtom <<" "<< angle<<" "
        << dihedAtom <<" "<< dihed<<endl;
   {
      scattPow=mpCryst->GetScatteringPowerRegistry().Find
                  (symbol,"ScatteringPowerAtom",true);
      if(scattPow==-1)
      {
         cout<<"Scattering power"<<symbol<<"not found, creating it..."<<endl;
         mpCryst->AddScatteringPower(new ScatteringPowerAtom(symbol,symbol));
      }
      scattPow=mpCryst->GetScatteringPowerRegistry().Find
                  (symbol,"ScatteringPowerAtom");
         sprintf(buf,"%d",2);
      this->AddAtom(symbol+(string)buf,
                    &(mpCryst->GetScatteringPowerRegistry().GetObj(scattPow)),
                    bondAtom-1,bond,
                    0,0,
                    0,0);
   }
   //third
      is >> symbol 
         >> bondAtom  >> bond
         >> angleAtom >> angle;
   cout << symbol <<" "
        << bondAtom  <<" "<< bond <<" "
        << angleAtom <<" "<< angle<<endl;
   {
      scattPow=mpCryst->GetScatteringPowerRegistry().Find
                  (symbol,"ScatteringPowerAtom",true);
      if(scattPow==-1)
      {
         cout<<"Scattering power"<<symbol<<"not found, creating it..."<<endl;
         mpCryst->AddScatteringPower(new ScatteringPowerAtom(symbol,symbol));
      }
      scattPow=mpCryst->GetScatteringPowerRegistry().Find
                  (symbol,"ScatteringPowerAtom");
         sprintf(buf,"%d",3);
      this->AddAtom(symbol+(string)buf,
                    &(mpCryst->GetScatteringPowerRegistry().GetObj(scattPow)),
                    bondAtom-1,bond,
                    angleAtom-1,angle*DEG2RAD,
                    0,0);
   }
   for(int i=3;i<nbAtoms;i++)
   {
      is >> symbol 
         >> bondAtom  >> bond
         >> angleAtom >> angle
         >> dihedAtom >> dihed;
      cout << symbol <<" "
           << bondAtom  <<" "<< bond <<" "
           << angleAtom <<" "<< angle
           << dihedAtom <<" "<< dihed<<endl;
      {
         scattPow=mpCryst->GetScatteringPowerRegistry().Find
                     (symbol,"ScatteringPowerAtom",true);
         if(scattPow==-1)
         {
            cout<<"Scattering power"<<symbol<<"not found, creating it..."<<endl;
            mpCryst->AddScatteringPower(new ScatteringPowerAtom(symbol,symbol));
         }
         scattPow=mpCryst->GetScatteringPowerRegistry().Find
                     (symbol,"ScatteringPowerAtom");
         cout<<
            sprintf(buf,"%d",i+1);
         this->AddAtom(symbol+(string)buf,
                       &(mpCryst->GetScatteringPowerRegistry().GetObj(scattPow)),
                       bondAtom-1,bond,
                       angleAtom-1,angle*DEG2RAD,
                       dihedAtom-1,dihed*DEG2RAD);
      }
   }
   this->SetLimitsRelative(gpRefParTypeScattConformBondLength,-.03,.03);
   this->SetLimitsRelative(gpRefParTypeScattConformBondAngle,-.01,.01);
   this->SetLimitsRelative(gpRefParTypeScattConformDihedAngle,-.01,.01);
}
void ZScatterer::ExportFenskeHallZMatrix(ostream &os)
{
   if(mNbAtom<1) return;
   os << mNbAtom<<endl;
   os << this->GetZAtomRegistry().GetObj(0).mpScattPow->GetName()
      << " 1"<<endl;
   if(mNbAtom<2) return;
   os << this->GetZAtomRegistry().GetObj(1).mpScattPow->GetName()
      << " "<<this->GetZBondAtom(1)+1<< " "<<this->GetZBondLength(1)
      <<endl;
   if(mNbAtom<3) return;
   os << this->GetZAtomRegistry().GetObj(2).mpScattPow->GetName()
      << " "<<this->GetZBondAtom(2)+1 << " "<<this->GetZBondLength(2)
      << " "<<this->GetZAngleAtom(2)+1<< " "<<this->GetZAngle(2)*RAD2DEG
      <<endl;
   if(mNbAtom<3) return;
   os << this->GetZAtomRegistry().GetObj(3).mpScattPow->GetName()
      << " "<<this->GetZBondAtom(3)+1 << " "<<this->GetZBondLength(3)
      << " "<<this->GetZAngleAtom(3)+1<< " "<<this->GetZAngle(3)*RAD2DEG
      << " "<<this->GetZDihedralAngleAtom(3)+1<< " "<<this->GetZDihedralAngle(3)*RAD2DEG
      <<endl;
   for(int i=3;i<mNbAtom;i++)
      os << this->GetZAtomRegistry().GetObj(i).mpScattPow->GetName()
         << " "<<this->GetZBondAtom(i)+1 << " "<<this->GetZBondLength(i)
         << " "<<this->GetZAngleAtom(i)+1<< " "<<this->GetZAngle(i)*RAD2DEG
         << " "<<this->GetZDihedralAngleAtom(i)+1<< " "<<this->GetZDihedralAngle(i)*RAD2DEG
         <<endl;
}

void ZScatterer::GlobalOptRandomMove(const REAL mutationAmplitude)
{
   if(mRandomMoveIsDone) return;
   VFN_DEBUG_ENTRY("ZScatterer::GlobalOptRandomMove()",3)
   TAU_PROFILE("ZScatterer::GlobalOptRandomMove()","void ()",TAU_DEFAULT);
   // give a 2% chance of either moving a single atom, or move
   // all atoms before a given torsion angle.
   // Only try this if there are more than 10 atoms (else it's not worth the speed cost)
   
   if((mNbAtom>=10) && ((rand()/(REAL)RAND_MAX)<.02))//.01
   {
      TAU_PROFILE_TIMER(timer1,\
                     "ZScatterer::GlobalOptRandomMoveSmart1(prepare ref par & mutate)"\
                     ,"", TAU_FIELD);
      TAU_PROFILE_TIMER(timer2,\
                     "ZScatterer::GlobalOptRandomMoveSmart2(optimize if necessary)"\
                     ,"", TAU_FIELD);
      TAU_PROFILE_START(timer1);
      // Do we have *any* dihedral angle to really move ?
         CrystVector_long dihed(mNbAtom);
         dihed=0;
         int nbDihed=0;
         RefinablePar *par;
         dihed(nbDihed++)=2;//This is the Psi angle, in fact. Should Chi be added, too ?
         for(int i=3;i<mNbAtom;i++)
         {
            par=&(this->GetPar(&(mZAtomRegistry.GetObj(i).mDihed)));
            if( !(par->IsFixed()) ) //&& !(par->IsLimited())
               dihed(nbDihed++)=i;
         }
         if(nbDihed<2) //Can't play :-(
         {
            this->RefinableObj::GlobalOptRandomMove(mutationAmplitude);
            VFN_DEBUG_EXIT("ZScatterer::GlobalOptRandomMove():End",3)
            return;
         }
      // Build mpZMoveMinimizer object
      if(0==mpZMoveMinimizer)
      {
         mpZMoveMinimizer=new ZMoveMinimizer(*this);
         // copy all parameters (and not reference them, we will change the fix status..
         // we could remove all popu parameters...
         for(int i=0; i<this->GetNbPar();i++) mpZMoveMinimizer->AddPar(this->GetPar(i));
      }
      // Pick one to move and get the relevant parameter
      // (maybe we should random-move also the associated bond lengths an angles,
      // but for now we'll concentrate on dihedral (torsion) angles.
         const int atom=dihed((int) (rand()/((REAL)RAND_MAX+1)*nbDihed));
         //cout<<endl;
         VFN_DEBUG_MESSAGE("ZScatterer::GlobalOptRandomMove(): Changing atom #"<<atom ,3)
         if(atom==2)
            par=&(this->GetPar(&mPsi));
         else
            par=&(this->GetPar(&(mZAtomRegistry.GetObj(atom).mDihed)));
         VFN_DEBUG_MESSAGE("ZScatterer::GlobalOptRandomMove(): initial value:"<<par->GetHumanValue() ,3)
      // Record the current conformation
         mpZMoveMinimizer->RecordConformation();
      // Set up
         const int moveType= rand()%3;
         mpZMoveMinimizer->FixAllPar();
         REAL x0,y0,z0;
         //cout << " Move Type:"<<moveType<<endl;
         CrystVector_REAL weight(mNbAtom);
         switch(moveType)
         {// :TODO: rather build a connectivity table, to include more atoms
            case 0:// (0) Try to move only one atom
            {
               VFN_DEBUG_MESSAGE("ZScatterer::GlobalOptRandomMove():Move only one atom",3)
               weight=1;
               weight(atom)=0;
               for(int i=0;i<mNbAtom;i++)
                  if(  (i==this->GetZBondAtom(atom))
                     ||(i==this->GetZAngleAtom(atom))
                     ||(i==this->GetZDihedralAngleAtom(atom))
                     ||(atom==this->GetZBondAtom(i))
                     ||(atom==this->GetZAngleAtom(i))
                     ||(atom==this->GetZDihedralAngleAtom(i)))
                  {
                     if(!(this->GetPar(&(mZAtomRegistry.GetObj(i).mBondLength)).IsFixed()))
                        mpZMoveMinimizer->GetPar(&(mZAtomRegistry.GetObj(i).mBondLength)).SetIsFixed(false);
                     if(!(this->GetPar(&(mZAtomRegistry.GetObj(i).mAngle)).IsFixed()))
                        mpZMoveMinimizer->GetPar(&(mZAtomRegistry.GetObj(i).mAngle)).SetIsFixed(false);
                     if(!(this->GetPar(&(mZAtomRegistry.GetObj(i).mDihed)).IsFixed()))
                        mpZMoveMinimizer->GetPar(&(mZAtomRegistry.GetObj(i).mDihed)).SetIsFixed(false);
                  }
               if(!(this->GetPar(&mPhi).IsFixed()))
                  mpZMoveMinimizer->GetPar(&mPhi).SetIsFixed(false);
               if(!(this->GetPar(&mChi).IsFixed()))
                  mpZMoveMinimizer->GetPar(&mChi).SetIsFixed(false);
               if( !(this->GetPar(&mPsi).IsFixed()) && (atom!=2))
                  mpZMoveMinimizer->GetPar(&mPsi).SetIsFixed(false);
               if(!(this->GetPar(&mXYZ(0)).IsFixed()))
                  mpZMoveMinimizer->GetPar(&mXYZ(0)).SetIsFixed(false);
               if(!(this->GetPar(&mXYZ(1)).IsFixed()))
                  mpZMoveMinimizer->GetPar(&mXYZ(1)).SetIsFixed(false);
               if(!(this->GetPar(&mXYZ(2)).IsFixed()))
                  mpZMoveMinimizer->GetPar(&mXYZ(2)).SetIsFixed(false);
               break;
            }
            case 1:// (1) Try to move the atoms *before* the rotated bond
            {
               VFN_DEBUG_MESSAGE("ZScatterer::GlobalOptRandomMove():Move before the rotated bond",3)
               const int atom1=this->GetZBondAtom(atom);
               const int atom2=this->GetZAngleAtom(atom);
               weight=0;
               weight(atom1)=1;
               for(int i=0;i<mNbAtom;i++) if(weight(this->GetZBondAtom(i))>.1) weight(i)=1;
               weight(atom2)=1;
               for(int i=0;i<mNbAtom;i++)
                  if(  (atom2 ==this->GetZAngleAtom(i))
                     &&(atom1 !=this->GetZBondAtom(i)) 
                     &&(i != atom) )
                  {
                     if(!(this->GetPar(&(mZAtomRegistry.GetObj(i).mBondLength)).IsFixed()))
                        mpZMoveMinimizer->GetPar(&(mZAtomRegistry.GetObj(i).mBondLength)).SetIsFixed(false);
                     if(!(this->GetPar(&(mZAtomRegistry.GetObj(i).mAngle)).IsFixed()))
                        mpZMoveMinimizer->GetPar(&(mZAtomRegistry.GetObj(i).mAngle)).SetIsFixed(false);
                     if(!(this->GetPar(&(mZAtomRegistry.GetObj(i).mDihed)).IsFixed()))
                        mpZMoveMinimizer->GetPar(&(mZAtomRegistry.GetObj(i).mDihed)).SetIsFixed(false);
                  }
               if(!(this->GetPar(&mPhi).IsFixed()))
                  mpZMoveMinimizer->GetPar(&mPhi).SetIsFixed(false);
               if(!(this->GetPar(&mChi).IsFixed()))
                  mpZMoveMinimizer->GetPar(&mChi).SetIsFixed(false);
               if( !(this->GetPar(&mPsi).IsFixed()) && (atom!=2))
                  mpZMoveMinimizer->GetPar(&mPsi).SetIsFixed(false);
                  
               if(!(this->GetPar(&mXYZ(0)).IsFixed()))
                  mpZMoveMinimizer->GetPar(&mXYZ(0)).SetIsFixed(false);
               if(!(this->GetPar(&mXYZ(1)).IsFixed()))
                  mpZMoveMinimizer->GetPar(&mXYZ(1)).SetIsFixed(false);
               if(!(this->GetPar(&mXYZ(2)).IsFixed()))
                  mpZMoveMinimizer->GetPar(&mXYZ(2)).SetIsFixed(false);
                  
               if(!(this->GetPar(&(mZAtomRegistry.GetObj(atom).mBondLength)).IsFixed()))
                  mpZMoveMinimizer->GetPar(&(mZAtomRegistry.GetObj(atom).mBondLength)).SetIsFixed(false);
               if(!(this->GetPar(&(mZAtomRegistry.GetObj(atom).mAngle)).IsFixed()))
                  mpZMoveMinimizer->GetPar(&(mZAtomRegistry.GetObj(atom).mAngle)).SetIsFixed(false);
               
               break;
            }
            case 2:// (2) Try to move the atoms *after* the rotated bond
            {
               if(mCenterAtomIndex>=atom)
               {
               VFN_DEBUG_MESSAGE("ZScatterer::GlobalOptRandomMove():Move after the bond (translate)",3)
                  x0=this->GetXCoord()(0);
                  y0=this->GetYCoord()(0);
                  z0=this->GetZCoord()(0);
               }
               else
               {
                  VFN_DEBUG_MESSAGE("ZScatterer::GlobalOptRandomMove():Move after the bond(nothing to do)",3)
               }
               break;
            }
         }
      // Move it, and with some probability use flipping to some
      // not-so-random angles., and then minimize the conformation change
         mpZMoveMinimizer->SetZAtomWeight(weight);
         REAL change;
         if( (rand()%5)==0)
         {
            switch(rand()%5)
            {
               case 0: change=-120*DEG2RAD;break;
               case 1: change= -90*DEG2RAD;break;
               case 2: change=  90*DEG2RAD;break;
               case 3: change= 120*DEG2RAD;break;
               default:change= 180*DEG2RAD;break;
            }
         }
         else
         {
            change= par->GetGlobalOptimStep()
                         *2*(rand()/(REAL)RAND_MAX-0.5)*mutationAmplitude*16;
         }
      TAU_PROFILE_STOP(timer1);
         VFN_DEBUG_MESSAGE("ZScatterer::GlobalOptRandomMove(): mutation:"<<change*RAD2DEG,3)
         par->Mutate(change);
         if(2==moveType)
         {
            if(mCenterAtomIndex>=atom)
            {
               this->UpdateCoordinates();
               x0 -= this->GetXCoord()(0);
               y0 -= this->GetYCoord()(0);
               z0 -= this->GetZCoord()(0);
               mpCryst->OrthonormalToFractionalCoords(x0,y0,z0);
               this->GetPar(mXYZ.data()).Mutate(x0);
               this->GetPar(mXYZ.data()+1).Mutate(y0);
               this->GetPar(mXYZ.data()+2).Mutate(z0);
            }
         }
         else
         {
            //cout <<"ZMoveMinimizer(flip):cost="<<mpZMoveMinimizer->GetCostFunctionValue(0);
            const REAL tmp=mpZMoveMinimizer->GetCostFunctionValue(0);
            if(tmp>.05)
            {
               TAU_PROFILE_START(timer2);
               if(tmp<1) mpZMoveMinimizer->MinimizeChange(100);
               else if(tmp<5) mpZMoveMinimizer->MinimizeChange(200);
                    else mpZMoveMinimizer->MinimizeChange(500);
               TAU_PROFILE_STOP(timer2);
            }
            //cout <<" -> "<<mpZMoveMinimizer->GetCostFunctionValue(0)<<endl;
         }
         VFN_DEBUG_MESSAGE("ZScatterer::GlobalOptRandomMove(): final value:"<<par->GetHumanValue(),3)
      //
   }
   #if 0
   {
      // find unfixed,  dihedral angles to play with //unlimited?
      CrystVector_long dihed(mNbAtom);
      dihed=0;
      int nbDihed=0;
      RefinablePar *par;
      dihed(nbDihed++)=2;
      for(int i=3;i<mNbAtom;i++)
      {
         par=&(this->GetPar(&(mZAtomRegistry.GetObj(i).mDihed)));
         if( !(par->IsFixed()) ) //&& !(par->IsLimited())
            dihed(nbDihed++)=i;
      }
      if(nbDihed<2) //Can't play :-(
         this->RefinableObj::GlobalOptRandomMove(mutationAmplitude);
      // Pick one
      const int atom=dihed((int) (rand()/((REAL)RAND_MAX+1)*nbDihed));
      VFN_DEBUG_MESSAGE("ZScatterer::GlobalOptRandomMove(): "<<FormatHorizVector<long>(dihed) ,10)
      VFN_DEBUG_MESSAGE("ZScatterer::GlobalOptRandomMove(): Changing atom #"<<atom ,10)
      if(atom==2)
         par=&(this->GetPar(&mPsi));
      else
         par=&(this->GetPar(&(mZAtomRegistry.GetObj(atom).mDihed)));
      // Get the old value
      const REAL old=par->GetValue();
      // Move it, with a max amplitude 8x greater than usual
      if( (rand()/(REAL)RAND_MAX)<.1)
      {// give some probability to use certain angles: -120,-90,90,120,180
         switch(rand()%5)
         {
            case 0: par->Mutate(-120*!DEG2RAD);break;
            case 1: par->Mutate( -90*!DEG2RAD);break;
            case 2: par->Mutate(  90*!DEG2RAD);break;
            case 3: par->Mutate( 120*!DEG2RAD);break;
            default:par->Mutate( 180*!DEG2RAD);break;
         }
      }
      else
         par->Mutate( par->GetGlobalOptimStep()
                      *2*(rand()/(REAL)RAND_MAX-0.5)*mutationAmplitude*8);
      const REAL change=mZAtomRegistry.GetObj(atom).GetZDihedralAngle()-old;
      // Now move all atoms using this changed bond as a reference
      //const int atom2=   mZAtomRegistry.GetObj(atom).GetZAngleAtom();
      for(int i=atom;i<mNbAtom;i++)
         //if(  (mZAtomRegistry.GetObj(i).GetZBondAtom()==atom)
         //   &&(mZAtomRegistry.GetObj(i).GetZAngleAtom()==atom2))
         //      this->GetPar(&(mZAtomRegistry.GetObj(i).mDihed)).Mutate(-change);
         if(mZAtomRegistry.GetObj(i).GetZAngleAtom()==atom)
               this->GetPar(&(mZAtomRegistry.GetObj(i).mDihed)).Mutate(-change);
      //cout <<"ZScatterer::GlobalOptRandomMove:"<<nbDihed
      //     <<" "<<atom
      //     <<" "<<atom2
      //     <<" "<<rand()
      //     <<endl
      //     <<" "<<FormatHorizVector<long>(dihed,4)
      //     <<endl;
   }
   #endif
   else
   {
      this->RefinableObj::GlobalOptRandomMove(mutationAmplitude);
   }
   mRandomMoveIsDone=true;
   VFN_DEBUG_EXIT("ZScatterer::GlobalOptRandomMove():End",3)
}

void ZScatterer::UpdateCoordinates() const
{
   if(mClockCoord>mClockScatterer) return;

   VFN_DEBUG_ENTRY("ZScatterer::UpdateCoordinates():"<<this->GetName(),3)
   TAU_PROFILE("ZScatterer::UpdateCoordinates()","void ()",TAU_DEFAULT);
   //if(0==mNbAtom) throw ObjCrystException("ZScatterer::Update() No Atoms in Scatterer !");
   if(0==mNbAtom) return;
   
   {
      CrystMatrix_REAL phiMatrix(3,3),chiMatrix(3,3),psiMatrix(3,3);
      phiMatrix= cos(mPhi)   , -sin(mPhi)   , 0,
                 sin(mPhi)   , cos(mPhi)    , 0,
                 0           ,0             ,1;

      chiMatrix= cos(mChi)   ,0             ,-sin(mChi),
                 0           ,1             ,0,
                 sin(mChi)   ,0             ,cos(mChi);

      psiMatrix= 1           , 0            , 0,
                 0           ,cos(mPsi)     ,-sin(mPsi),
                 0           ,sin(mPsi)     ,cos(mPsi);
                 
      mPhiChiPsiMatrix=product(chiMatrix,product(phiMatrix,psiMatrix));
      //cout << phiMatrix <<endl<< chiMatrix <<endl<< psiMatrix <<endl<<mPhiChiPsiMatrix<<endl;
   }

   mXCoord.resize(mNbAtom);
   mYCoord.resize(mNbAtom);
   mZCoord.resize(mNbAtom);
   
   // Atom 0
   mXCoord(0)=0.;
   mYCoord(0)=0.;
   mZCoord(0)=0.;
   VFN_DEBUG_MESSAGE("->Atom #0:"<<mXCoord(0)<<" : "<<mYCoord(0)<<" : "<<mZCoord(0),1)
   
   if(mNbAtom>1)
   {// Atom 1
      mXCoord(1)=GetZBondLength(1);
      mYCoord(1)=0.;
      mZCoord(1)=0.;
   }
   VFN_DEBUG_MESSAGE("->Atom #1:"<<mXCoord(1)<<" : "<<mYCoord(1)<<" : "<<mZCoord(1),1)
   if(mNbAtom>2)
   {// Atom 2
      if(0==GetZBondAtom(2)) //Linked with Atom 1
         mXCoord(2)=GetZBondLength(2)*cos(GetZAngle(2));
      else //Linked with Atom 1
         mXCoord(2)=mXCoord(1)-GetZBondLength(2)*cos(GetZAngle(2));
      mYCoord(2)=GetZBondLength(2)*sin(GetZAngle(2));
      mZCoord(2)=0.;
   }
   VFN_DEBUG_MESSAGE("->Atom #2:"<<mXCoord(2)<<" : "<<mYCoord(2)<<" : "<<mZCoord(2),1)
   for(int i=1;i<3;i++)// Global rotation of scatterer
   {
      const REAL x=mXCoord(i);
      const REAL y=mYCoord(i);
      const REAL z=mZCoord(i);
      mXCoord(i)=mPhiChiPsiMatrix(0,0)*x+mPhiChiPsiMatrix(0,1)*y+mPhiChiPsiMatrix(0,2)*z;
      mYCoord(i)=mPhiChiPsiMatrix(1,0)*x+mPhiChiPsiMatrix(1,1)*y+mPhiChiPsiMatrix(1,2)*z;
      mZCoord(i)=mPhiChiPsiMatrix(2,0)*x+mPhiChiPsiMatrix(2,1)*y+mPhiChiPsiMatrix(2,2)*z;
      if(i==mNbAtom-1) break;
   }
   if(mNbAtom>3)
   {
      REAL xa,ya,za,xb,yb,zb,xd,yd,zd,cosph,sinph,costh,sinth,coskh,sinkh,cosa,sina;
      REAL xpd,ypd,zpd,xqd,yqd,zqd;
      REAL rbc,xyb,yza,tmp,xpa,ypa,zqa;
      int na,nb,nc;
      bool flag;
      REAL dist,angle,dihed;
      for(int i=3;i<mNbAtom;i++)
      {
         na=GetZBondAtom(i);
         nb=GetZAngleAtom(i);
         nc=GetZDihedralAngleAtom(i);
         dist=GetZBondLength(i);
         angle=GetZAngle(i);
         dihed=GetZDihedralAngle(i);
         
         xb = mXCoord(nb) - mXCoord(na);
         yb = mYCoord(nb) - mYCoord(na);
         zb = mZCoord(nb) - mZCoord(na);

         rbc= 1./sqrt(xb*xb + yb*yb + zb*zb);

         cosa = cos(angle);
         sina = sin(angle);

         if( fabs(cosa) >= 0.999999 )
         {   // Colinear
            mXCoord(i)=mXCoord(na)+cosa*dist*rbc*xb;
            mYCoord(i)=mYCoord(na)+cosa*dist*rbc*yb;
            mZCoord(i)=mZCoord(na)+cosa*dist*rbc*zb;
            VFN_DEBUG_MESSAGE("->Atom #"<<i<<":"<<mXCoord(i)<<" : "<<mYCoord(i)<<" : " <<mZCoord(i)<<"(colinear)",1)
         }
         else
         {
            xa = mXCoord(nc) - mXCoord(na);
            ya = mYCoord(nc) - mYCoord(na);
            za = mZCoord(nc) - mZCoord(na);
            xd = dist*cosa;
            yd = dist*sina*cos(dihed);
            zd = -dist*sina*sin(dihed);

            xyb = sqrt(xb*xb + yb*yb);
            if( xyb < 0.1 )
            {   // Rotate about y-axis
                tmp = za; za = -xa; xa = tmp;
                tmp = zb; zb = -xb; xb = tmp;
                xyb = sqrt(xb*xb + yb*yb);
                flag = true;
            }
            else flag = false;

            costh = xb/xyb;
            sinth = yb/xyb;
            xpa = costh*xa + sinth*ya;
            ypa = costh*ya - sinth*xa;
            sinph = zb*rbc;
            cosph = sqrt(1.0 - sinph*sinph);
            zqa = cosph*za  - sinph*xpa;
            yza = sqrt(ypa*ypa + zqa*zqa);
            
            //cout <<endl<<ypa<<" "<<zqa<<" "<<yza<<" "<<endl;
            if( yza > 1e-5 )
            {
               coskh = ypa/yza;
               sinkh = zqa/yza;
               ypd = coskh*yd - sinkh*zd;
               zpd = coskh*zd + sinkh*yd;
            }
            else
            {
               ypd = yd;
               zpd = zd;
            }

            xpd = cosph*xd  - sinph*zpd;
            zqd = cosph*zpd + sinph*xd;
            xqd = costh*xpd - sinth*ypd;
            yqd = costh*ypd + sinth*xpd;

            if( true==flag )
            {   // Rotate about y-axis ?
               mXCoord(i)=mXCoord(na) - zqd;
               mYCoord(i)=mYCoord(na) + yqd;
               mZCoord(i)=mZCoord(na) + xqd;
            } else
            { 
               mXCoord(i)=mXCoord(na) + xqd;
               mYCoord(i)=mYCoord(na) + yqd;
               mZCoord(i)=mZCoord(na) + zqd;
            }
            VFN_DEBUG_MESSAGE("->Atom #"<<i<<":"<<mXCoord(i)<<" : "<<mYCoord(i)<<" : " <<mZCoord(i),1)
         }
      }
   }
   //shift atom around Central atom
   REAL x,y,z;
   x=this->GetX();
   y=this->GetY();
   z=this->GetZ();
   mpCryst->FractionalToOrthonormalCoords(x,y,z);
   const REAL x0=x-mXCoord(mCenterAtomIndex);
   const REAL y0=y-mYCoord(mCenterAtomIndex);
   const REAL z0=z-mZCoord(mCenterAtomIndex);
   for(int i=0;i<mNbAtom;i++)
   {
      mXCoord(i) += x0;
      mYCoord(i) += y0;
      mZCoord(i) += z0;
   }
   mClockCoord.Click();
   VFN_DEBUG_EXIT("ZScatterer::UpdateCoordinates()"<<this->GetName(),3)
}
   
void ZScatterer::UpdateScattCompList() const
{
   VFN_DEBUG_ENTRY("ZScatterer::UpdateScattCompList()"<<this->GetName(),3)
   this->UpdateCoordinates();
   if(  (mClockScattCompList>mClockCoord)
      &&(mClockScattCompList>mpCryst->GetClockLatticePar())) return;

   if(true==mUseGlobalScattPow)
   {
      mScattCompList(0).mX=mXYZ(0);
      mScattCompList(0).mY=mXYZ(1);
      mScattCompList(0).mZ=mXYZ(2);
      mScattCompList(0).mOccupancy=mOccupancy;
      mScattCompList(0).mpScattPow=mpGlobalScattPow;
      
      // Update central atom for Display.
      //mXCoord(mCenterAtomIndex)=this->GetX();
      //mYCoord(mCenterAtomIndex)=this->GetY();
      //mZCoord(mCenterAtomIndex)=this->GetZ();
      //mpCryst->FractionalToOrthonormalCoords(mXCoord(mCenterAtomIndex),
      //                                       mXCoord(mCenterAtomIndex),
      //                                       mXCoord(mCenterAtomIndex));
      
      mClockScattCompList.Click();
      VFN_DEBUG_MESSAGE("ZScatterer::UpdateScattCompList()->Global Scatterer:End",3)
      return;
   }

   long j=0;
   REAL x,y,z;
   VFN_DEBUG_MESSAGE("ZScatterer::UpdateScattCompList(bool):Finishing"<<mNbAtom<<","<<mNbDummyAtom,3)
   for(long i=0;i<mNbAtom;i++)
   {
      if(0!=mZAtomRegistry.GetObj(i).GetScatteringPower())
      {
         x=mXCoord(i);
         y=mYCoord(i);
         z=mZCoord(i);
         mpCryst->OrthonormalToFractionalCoords(x,y,z);
         mScattCompList(j  ).mX=x;
         mScattCompList(j  ).mY=y;
         mScattCompList(j  ).mZ=z;
         mScattCompList(j++).mOccupancy=mZAtomRegistry.GetObj(i).GetOccupancy();
      }
   }
   #ifdef __DEBUG__
   if(gVFNDebugMessageLevel<3) mScattCompList.Print();
   #endif
   mClockScattCompList.Click();
   VFN_DEBUG_EXIT("ZScatterer::UpdateScattCompList()"<<this->GetName(),3)
}

void ZScatterer::InitRefParList()
{
   VFN_DEBUG_MESSAGE("ZScatterer::InitRefParList():"<<this->GetName(),5)
   //throw ObjCrystException("ZScatterer::InitRefParList() not implemented! "+this->GetName());
   //:TODO:
   this->ResetParList();
   {
      RefinablePar tmp("x",&mXYZ(0),0.,1.,
                        gpRefParTypeScattTranslX,
                        REFPAR_DERIV_STEP_ABSOLUTE,false,false,true,true,1.,1.);
      tmp.AssignClock(mClockScatterer);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("y",&mXYZ(1),0,1,
                        gpRefParTypeScattTranslY,
                        REFPAR_DERIV_STEP_ABSOLUTE,false,false,true,true,1.,1.);
      tmp.AssignClock(mClockScatterer);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("z",&mXYZ(2),0,1,
                        gpRefParTypeScattTranslZ,
                        REFPAR_DERIV_STEP_ABSOLUTE,false,false,true,true,1.,1.);
      tmp.AssignClock(mClockScatterer);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("Occupancy",&mOccupancy,0,1,
                        gpRefParTypeScattOccup,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,1.,1.);
      tmp.AssignClock(mClockScatterer);
      this->AddPar(tmp);
   }
   if(false==mUseGlobalScattPow)
   {
      {
         RefinablePar tmp("Phi",&mPhi,0,2*M_PI,
                           gpRefParTypeScattOrient,
                           REFPAR_DERIV_STEP_ABSOLUTE,false,false,true,true,RAD2DEG,2*M_PI);
         tmp.AssignClock(mClockScatterer);
         this->AddPar(tmp);
      }
      {
         RefinablePar tmp("Chi",&mChi,0,2*M_PI,
                           gpRefParTypeScattOrient,
                           REFPAR_DERIV_STEP_ABSOLUTE,false,false,true,true,RAD2DEG,2*M_PI);
         tmp.AssignClock(mClockScatterer);
         this->AddPar(tmp);
      }
      {
         RefinablePar tmp("Psi",&mPsi,0,2*M_PI,
                           gpRefParTypeScattOrient,
                           REFPAR_DERIV_STEP_ABSOLUTE,false,false,true,true,RAD2DEG,2*M_PI);
         tmp.AssignClock(mClockScatterer);
         this->AddPar(tmp);
      }
      char buf [20];
      for(long i=0;i<mNbAtom;i++)
      {
         {
            sprintf(buf,"%d-%d",(int)i,(int)(mZAtomRegistry.GetObj(i).GetZBondAtom()));
            RefinablePar tmp((string)"Length"+(string)buf,
                              &(mZAtomRegistry.GetObj(i).mBondLength),
                              mZAtomRegistry.GetObj(i).mBondLength*.9,
                              mZAtomRegistry.GetObj(i).mBondLength*1.1,
                              gpRefParTypeScattConformBondLength,
                              REFPAR_DERIV_STEP_ABSOLUTE,true,false,true,false,1.);
            tmp.AssignClock(mClockScatterer);
            this->AddPar(tmp);
         }
         {
            sprintf(buf,"%d-%d-%d",(int)i,(int)(mZAtomRegistry.GetObj(i).GetZBondAtom()),
                                          (int)(mZAtomRegistry.GetObj(i).GetZAngleAtom()));
            RefinablePar tmp("Angle"+(string)buf,
                              &(mZAtomRegistry.GetObj(i).mAngle),0,2*M_PI,
                              gpRefParTypeScattConformBondAngle,
                              REFPAR_DERIV_STEP_ABSOLUTE,false,false,true,true,RAD2DEG,2*M_PI);
            tmp.AssignClock(mClockScatterer);
            this->AddPar(tmp);
         }
         {
            sprintf(buf,"%d-%d-%d-%d",(int)i,(int)(mZAtomRegistry.GetObj(i).GetZBondAtom()),
                                      (int)(mZAtomRegistry.GetObj(i).GetZAngleAtom()),
                                      (int)(mZAtomRegistry.GetObj(i).GetZDihedralAngleAtom()));
            RefinablePar tmp("Dihed"+(string)buf,
                              &(mZAtomRegistry.GetObj(i).mDihed),0,2*M_PI,
                              gpRefParTypeScattConformDihedAngle,
                              REFPAR_DERIV_STEP_ABSOLUTE,false,false,true,true,RAD2DEG,2*M_PI);
            tmp.AssignClock(mClockScatterer);
            this->AddPar(tmp);
         }
         if(0!=mZAtomRegistry.GetObj(i).GetScatteringPower())
         {//fixed by default
            sprintf(buf,"%d",(int)i);
            RefinablePar tmp("Occupancy"+(string)buf, 
                              &(mZAtomRegistry.GetObj(i).mOccupancy),0,1,
                              gpRefParTypeScattOccup,
                              REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,1.,1.);
            tmp.AssignClock(mClockScatterer);
            this->AddPar(tmp);
         }
      }
   }//if(false==mUseGlobalScatteringPower)
}

#ifdef __WX__CRYST__
WXCrystObjBasic* ZScatterer::WXCreate(wxWindow* parent)
{
   //:TODO: Check mpWXCrystObj==0
   mpWXCrystObj=new WXZScatterer(parent,this);
   return mpWXCrystObj;
}
#endif

//######################################################################
//
//  ZPolyhedron   
//
//
//######################################################################
ZPolyhedron::ZPolyhedron( const RegularPolyhedraType type,Crystal &cryst,
      const REAL x, const REAL y, const REAL z,
      const string &name, const ScatteringPower *centralAtomSymbol,
      const ScatteringPower *periphAtomSymbol,const REAL centralPeriphDist,
      const REAL ligandPopu,
      const REAL phi, const REAL chi, const REAL psi):
ZScatterer(name,cryst,x,y,z,phi,chi,psi),mPolyhedraType(type)
{
   VFN_DEBUG_MESSAGE("ZPolyhedron::ZPolyhedron(..)",5)
   // Additioning string and char arrays takes a huge lot of mem (gcc 2.95.2)
   const string name_=name+"_";
   const string name_central=name_+centralAtomSymbol->GetName();
   const string name_periph=name_+periphAtomSymbol->GetName();
   switch(mPolyhedraType)
   {
      case TETRAHEDRON :
      {
         REAL ang=2*asin(sqrt(2./3.));
         this ->AddAtom (name_central,   centralAtomSymbol,
                           0,0.,
                           0,0.,
                           0,0.,
                           1.);
         this ->AddAtom (name_periph+"1",periphAtomSymbol,
                           0,centralPeriphDist,
                           0,0.,
                           0,0.,
                           1.);
         this ->AddAtom (name_periph+"2",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,ang,
                           0,0.,
                           1.);
         this ->AddAtom (name_periph+"3",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,ang,
                           2, M_PI*2./3.,
                           1.);
         this ->AddAtom (name_periph+"4",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,ang,
                           2,-M_PI*2./3.,
                           1.);
         m3DDisplayIndex.resize(4,3);
         m3DDisplayIndex=  1,2,3,
                                 1,2,4,
                                 1,3,4,
                                 2,3,4;
         break;
      }
      case  OCTAHEDRON:
      {
         this ->AddAtom (name_central,   centralAtomSymbol,
                           0,0.,
                           0,0.,
                           0,0.,
                           1.);
         this ->AddAtom (name_periph+"1",periphAtomSymbol,
                           0,centralPeriphDist,
                           0,0.,
                           0,0.,
                           1.);
         this ->AddAtom (name_periph+"2",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,M_PI/2,
                           0,0.,
                           1.);
         this ->AddAtom (name_periph+"3",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,M_PI/2,
                           2,M_PI/2,
                           1.);
         this ->AddAtom (name_periph+"4",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,M_PI/2,
                           2,-M_PI/2,
                           1.);
         this ->AddAtom (name_periph+"5",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,M_PI/2,
                           2,M_PI,
                           1.);
         this ->AddAtom (name_periph+"6",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,M_PI,
                           2,0.,
                           1.);
         m3DDisplayIndex.resize(8,3);
         m3DDisplayIndex=  1,2,3,
                                 1,2,4,
                                 1,5,3,
                                 1,5,4,
                                 6,2,3,
                                 6,2,4,
                                 6,5,3,
                                 6,5,4;
         break;
      }
      case  SQUARE_PLANE:
      {
         this ->AddAtom (name_central,   centralAtomSymbol,
                           0,0.,
                           0,0.,
                           0,0.,
                           1.);
         this ->AddAtom (name+"_X",0,            //Dummy atom on top
                           0,1.,
                           0,0.,
                           0,0.,
                           1.);
         this ->AddAtom (name_periph+"1",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,M_PI/2,
                           0,0.,
                           1.);
         this ->AddAtom (name_periph+"2",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,M_PI/2,
                           2,M_PI/2,
                           1.);
         this ->AddAtom (name_periph+"3",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,M_PI/2,
                           2,M_PI,
                           1.);
         this ->AddAtom (name_periph+"4",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,M_PI/2,
                           2,-M_PI/2,
                           1.);
         //:TODO: GL Display with squares
         //m3DDisplayIndex.resize(2,3);
         //m3DDisplayIndex=  2,4,2,
         //                        2,4,5;
         break;
      }
      case  CUBE:
      {
         this ->AddAtom (name_central,   centralAtomSymbol,
                           0,0.,
                           0,0.,
                           0,0.,
                           1.);
         this ->AddAtom (name+"_X",0, //Dummy atom in middle of face
                           0,1.,
                           0,0.,
                           0,0.,
                           1.);
         this ->AddAtom (name_periph+"1",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,M_PI/4.,
                           0,0.,
                           1.);
         this ->AddAtom (name_periph+"2",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,M_PI/4.,
                           2,M_PI/2.,
                           1.);
         this ->AddAtom (name_periph+"3",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,M_PI/4.,
                           2,M_PI,
                           1.);
         this ->AddAtom (name_periph+"4",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,M_PI/4.,
                           2,-M_PI/2.,
                           1.);
         this ->AddAtom (name_periph+"5",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,M_PI*3./4.,
                           0,0.,
                           1.);
         this ->AddAtom (name_periph+"6",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,M_PI*3./4.,
                           2,M_PI/2.,
                           1.);
         this ->AddAtom (name_periph+"7",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,M_PI*3./4.,
                           2,M_PI,
                           1.);
         this ->AddAtom (name_periph+"8",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,M_PI*3./4.,
                           2,M_PI*3./2.,
                           1.);
         break;
      }
      case  ANTIPRISM_TETRAGONAL:
      {
         this ->AddAtom (name_central,   centralAtomSymbol,
                           0,0.,
                           0,0.,
                           0,0.,
                           1.);
         this ->AddAtom (name+"_X",0, //Dummy atom in middle of face
                           0,1.,
                           0,0.,
                           0,0.,
                           1.);
         this ->AddAtom (name_periph+"1",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,M_PI/4.,
                           0,0.,
                           1.);
         this ->AddAtom (name_periph+"2",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,M_PI/4.,
                           2,M_PI/2.,
                           1.);
         this ->AddAtom (name_periph+"3",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,M_PI/4.,
                           2,M_PI,
                           1.);
         this ->AddAtom (name_periph+"4",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,M_PI/4.,
                           2,-M_PI/2.,
                           1.);
         this ->AddAtom (name_periph+"5",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,M_PI*3./4.,
                           0,M_PI/4.,
                           1.);
         this ->AddAtom (name_periph+"6",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,M_PI*3./4.,
                           2,M_PI*3./4.,
                           1.);
         this ->AddAtom (name_periph+"7",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,M_PI*3./4.,
                           2,M_PI*5./4.,
                           1.);
         this ->AddAtom (name_periph+"8",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,M_PI*3./4.,
                           2,M_PI*7./4.,
                           1.);
         break;
      }
      case  PRISM_TETRAGONAL_MONOCAP:
      {
         this ->AddAtom (name_central,   centralAtomSymbol,
                           0,0.,
                           0,0.,
                           0,0.,
                           1.);
         this ->AddAtom (name_periph+"0",periphAtomSymbol,
                           0,centralPeriphDist,
                           0,0.,
                           0,0.,
                           1.);
         this ->AddAtom (name_periph+"1",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,70*DEG2RAD,
                           0,0.,
                           1.);
         this ->AddAtom (name_periph+"2",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,70*DEG2RAD,
                           2,M_PI/2.,
                           1.);
         this ->AddAtom (name_periph+"3",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,70*DEG2RAD,
                           2,M_PI,
                           1.);
         this ->AddAtom (name_periph+"4",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,70*DEG2RAD,
                           2,-M_PI/2.,
                           1.);
         this ->AddAtom (name_periph+"5",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,145*DEG2RAD,
                           0,0.,
                           1.);
         this ->AddAtom (name_periph+"6",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,145*DEG2RAD,
                           2,M_PI/2.,
                           1.);
         this ->AddAtom (name_periph+"7",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,145*DEG2RAD,
                           2,M_PI,
                           1.);
         this ->AddAtom (name_periph+"8",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,145*DEG2RAD,
                           2,M_PI*3./2.,
                           1.);
         break;
      }
      case  PRISM_TETRAGONAL_DICAP:
      {
         this ->AddAtom (name_central,   centralAtomSymbol,
                           0,0.,
                           0,0.,
                           0,0.,
                           1.);
         this ->AddAtom (name_periph+"0",periphAtomSymbol,
                           0,centralPeriphDist,
                           0,0.,
                           0,0.,
                           1.);
         this ->AddAtom (name_periph+"1",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,60*DEG2RAD,
                           0,0.,
                           1.);
         this ->AddAtom (name_periph+"2",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,60*DEG2RAD,
                           2,M_PI/2.,
                           1.);
         this ->AddAtom (name_periph+"3",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,60*DEG2RAD,
                           2,M_PI,
                           1.);
         this ->AddAtom (name_periph+"4",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,60*DEG2RAD,
                           2,-M_PI/2.,
                           1.);
         this ->AddAtom (name_periph+"5",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,120*DEG2RAD,
                           0,0.,
                           1.);
         this ->AddAtom (name_periph+"6",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,120*DEG2RAD,
                           2,M_PI/2.,
                           1.);
         this ->AddAtom (name_periph+"7",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,120*DEG2RAD,
                           2,M_PI,
                           1.);
         this ->AddAtom (name_periph+"8",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,120*DEG2RAD,
                           2,M_PI*3./2.,
                           1.);
         this ->AddAtom (name_periph+"8",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,M_PI,
                           2,0.,
                           1.);
         break;
      }
      case  PRISM_TRIGONAL:
      {
         const REAL ang=55.*DEG2RAD;
         const REAL ang2=120.*DEG2RAD;
         this ->AddAtom (name_central,   centralAtomSymbol,
                           0,0.,
                           0,0.,
                           0,0.,
                           1.);
         this ->AddAtom (name+"_X",0, //Dummy atom in middle of top face
                           0,1.,
                           0,0.,
                           0,0.,
                           0.);
         this ->AddAtom (name_periph+"0",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,ang,
                           0,0.,
                           1.);
         this ->AddAtom (name_periph+"1",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,ang,
                           2,ang2,
                           1.);
         this ->AddAtom (name_periph+"2",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,ang,
                           2,-ang2,
                           1.);
         this ->AddAtom (name_periph+"3",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,M_PI-ang,
                           0,0.,
                           1.);
         this ->AddAtom (name_periph+"4",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,M_PI-ang,
                           2,ang2,
                           1.);
         this ->AddAtom (name_periph+"5",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,M_PI-ang,
                           2,-ang2,
                           1.);
         break;
      }
      case  PRISM_TRIGONAL_TRICAPPED:
      {
         const REAL ang=55.*DEG2RAD;
         const REAL ang2=120.*DEG2RAD;
         this ->AddAtom (name_central,   centralAtomSymbol,
                           0,0.,
                           0,0.,
                           0,0.,
                           1.);
         this ->AddAtom (name+"_X",0, //Dummy atom in middle of top face
                           0,1.,
                           0,0.,
                           0,0.,
                           0.);
         this ->AddAtom (name_periph+"0",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,ang,
                           0,0.,
                           1.);
         this ->AddAtom (name_periph+"1",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,ang,
                           2,ang2,
                           1.);
         this ->AddAtom (name_periph+"2",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,ang,
                           2,-ang2,
                           1.);
         this ->AddAtom (name_periph+"3",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,M_PI-ang,
                           0,0.,
                           1.);
         this ->AddAtom (name_periph+"4",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,M_PI-ang,
                           2,ang2,
                           1.);
         this ->AddAtom (name_periph+"5",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,M_PI-ang,
                           2,-ang2,
                           1.);
         this ->AddAtom (name_periph+"6",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,M_PI/2.,
                           2,M_PI,
                           1.);
         this ->AddAtom (name_periph+"7",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,M_PI/2.,
                           2,M_PI/3.,
                           1.);
         this ->AddAtom (name_periph+"8",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,M_PI/2.,
                           2,-M_PI/3.,
                           1.);
         break;
      }
      case  ICOSAHEDRON:
      {
         const REAL ang=acos(sqrt(.2));
         const REAL ang2=M_PI*2./5.;
         this ->AddAtom (name_central,   centralAtomSymbol,
                           0,0.,
                           0,0.,
                           0,0.,
                           1.);
         this ->AddAtom (name_periph+"0",periphAtomSymbol,
                           0,centralPeriphDist,
                           0,0.,
                           0,0.,
                          1.);
         this ->AddAtom (name_periph+"1",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,ang,
                           1,0.,
                           1.);
         this ->AddAtom (name_periph+"2",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,ang,
                           2,ang2,
                           1.);
         this ->AddAtom (name_periph+"3",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,ang,
                           2,2*ang2,
                           1.);
         this ->AddAtom (name_periph+"4",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,ang,
                           2,-2*ang2,
                           1.);
         this ->AddAtom (name_periph+"5",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,ang,
                           2,-ang2,
                           1.);
         this ->AddAtom (name_periph+"6",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,M_PI,
                           0,0.,
                           1.);
         this ->AddAtom (name_periph+"7",periphAtomSymbol,
                           0,centralPeriphDist,
                           7,ang,
                           3,M_PI,
                           1.);
         this ->AddAtom (name_periph+"8",periphAtomSymbol,
                           0,centralPeriphDist,
                           7,ang,
                           8,ang2,
                           1.);
         this ->AddAtom (name_periph+"9",periphAtomSymbol,
                           0,centralPeriphDist,
                           7,ang,
                           8,2*ang2,
                           1.);
         this ->AddAtom (name_periph+"10",periphAtomSymbol,
                           0,centralPeriphDist,
                           7,ang,
                           8,-2*ang2,
                           1.);
         this ->AddAtom (name_periph+"11",periphAtomSymbol,
                           0,centralPeriphDist,
                           7,ang,
                           8,-ang2,
                           1.);
         break;
      }
      case  TRIANGLE_PLANE:
      {
         this ->AddAtom (name_central,   centralAtomSymbol,
                           0,0.,
                           0,0.,
                           0,0.,
                           1.);
         this ->AddAtom (name+"_X",0,            //Dummy atom on top
                           0,1.,
                           0,0.,
                           0,0.,
                           1.);
         this ->AddAtom (name_periph+"1",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,M_PI/2,
                           0,0.,
                           1.);
         this ->AddAtom (name_periph+"2",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,M_PI/2,
                           2,2*M_PI/3,
                           1.);
         this ->AddAtom (name_periph+"3",periphAtomSymbol,
                           0,centralPeriphDist,
                           1,M_PI/2,
                           2,-2*M_PI/3,
                           1.);
         //:TODO: GL Display with squares
         //m3DDisplayIndex.resize(2,3);
         //m3DDisplayIndex=  2,4,2,
         //                        2,4,5;
         break;
      }
      default :   throw ObjCrystException("ZPolyhedron::ZPolyhedron():Unknown Polyhedra type !");
   }
   //We want to keep an approximate geometry
   //this->RefinableObj::Print();
   this->SetLimitsRelative(gpRefParTypeScattConformBondLength,-.2,.2);
   this->SetLimitsRelative(gpRefParTypeScattConformBondAngle ,-.2,.2);
   this->SetLimitsRelative(gpRefParTypeScattConformDihedAngle,-.2,.2);
   //this->RefinableObj::Print();
   
   VFN_DEBUG_MESSAGE("ZPolyhedron::ZPolyhedron():End:"<<mName<<")",5)
}

ZPolyhedron::ZPolyhedron(const ZPolyhedron &old):
ZScatterer(old){}

ZPolyhedron* ZPolyhedron::CreateCopy() const
{
   VFN_DEBUG_MESSAGE("ZPolyhedron::CreateCopy()"<<mName<<")",5)
   return new ZPolyhedron(*this);
}

//######################################################################
//
//      GLOBAL SCATTERING POWER
//
//######################################################################

GlobalScatteringPower::GlobalScatteringPower():mpZScatterer(0)
{
   VFN_DEBUG_MESSAGE("GlobalScatteringPower::GlobalScatteringPower():"<<mName,5)
}

GlobalScatteringPower::GlobalScatteringPower(const ZScatterer &scatt):mpZScatterer(0)
{
   VFN_DEBUG_MESSAGE("GlobalScatteringPower::GlobalScatteringPower(&scatt):"<<mName,5)
   this->Init(scatt);
}

GlobalScatteringPower::GlobalScatteringPower(const GlobalScatteringPower& old):mpZScatterer(0)
{
   VFN_DEBUG_MESSAGE("GlobalScatteringPower::GlobalScatteringPower(&old):"<<mName,5)
   this->Init(*(old.mpZScatterer));
}

GlobalScatteringPower::~GlobalScatteringPower()
{
   VFN_DEBUG_MESSAGE("GlobalScatteringPower::~GlobalScatteringPower():"<<mName,5)
   if(mpZScatterer!=0) delete mpZScatterer;
}

void GlobalScatteringPower::Init(const ZScatterer &scatt)
{
   this->SetName(scatt.GetName()+(string)"_Global");
   VFN_DEBUG_MESSAGE("GlobalScatteringPower::Init(&Scatt):"<<mName,5)
   if(mpZScatterer!=0) delete mpZScatterer;
   mpZScatterer=new ZScatterer(scatt);
   mpZScatterer->SetUseGlobalScatteringPower(false);
   mpZScatterer->SetName(scatt.GetName()+"_GlobalCopy");
   
   //Set the DynPopCorrIndex to the sum of the DynPopCorrIndexs (eg the sum of atomic numbers)
   mDynPopCorrIndex=0;
   
   const ScatteringComponentList* tmp=&(mpZScatterer->GetScatteringComponentList());
   for(long i=0;i<tmp->GetNbComponent();i++)
      mDynPopCorrIndex += (*tmp)(i).mpScattPow->GetDynPopCorrIndex();
}

CrystVector_REAL GlobalScatteringPower::
   GetScatteringFactor(const ScatteringData &data,
                       const int spgSymPosIndex) const
{
   // Here comes the hard work
   VFN_DEBUG_MESSAGE("GlobalScatteringPower::GetScatteringFactor():"<<mName,10)
   TAU_PROFILE("GlobalScatteringPower::GetScatteringFactor()","void ()",TAU_DEFAULT);
   
   // copy both the scatterer and the DiffractionData object to determine
   // the average isotropic scattering power for each reflection
   ScatteringData* pData=data.CreateCopy();
   CrystVector_REAL sf(data.GetNbRefl()),rsf,isf;
   sf=0;
   
   Crystal cryst(data.GetCrystal().GetLatticePar(0),data.GetCrystal().GetLatticePar(1),
                 data.GetCrystal().GetLatticePar(2),data.GetCrystal().GetLatticePar(3),
                 data.GetCrystal().GetLatticePar(4),data.GetCrystal().GetLatticePar(5),"P1");
   cryst.AddScatterer(mpZScatterer);
   cryst.SetUseDynPopCorr(false);
   pData->SetCrystal(cryst);
   pData->SetName("GlobalScatteringPowerData!");
   VFN_DEBUG_MESSAGE("GlobalScatteringPower::GetScatteringFactor():No DEBUG Messages",5)
   VFN_DEBUG_LOCAL_LEVEL(10)
   
   const long nbStep=4;//number of steps over 90 degrees for phi,chi psi
   REAL norm=0;
   for(int i=-nbStep;i<=nbStep;i++)
   {
      mpZScatterer->SetChi(i*M_PI/2/nbStep);
      for(int j=-nbStep;j<=nbStep;j++)
      {
         mpZScatterer->SetPhi(j*M_PI/2/nbStep);
         for(int k=-nbStep;k<=nbStep;k++)
         {
            //cout <<i<<","<<j<<","<<k<<endl;
            mpZScatterer->SetPsi(k*M_PI/2/nbStep);
            
            rsf=pData->GetFhklCalcReal();
            rsf*=rsf;
            isf=pData->GetFhklCalcImag();
            isf*=isf;
            rsf+=isf;
            rsf=sqrt(rsf);
            
            //correct for the solid angle (dChi*dPhi) corresponding to this orientation
            rsf*= cos(mpZScatterer->GetPhi());//:TODO: needs checking !!!
            norm += cos(mpZScatterer->GetPhi());
            sf += rsf;
         }
      }
      //cout << FormatHorizVector<REAL>(sf) <<endl<<norm<<endl;
   }
   
   VFN_DEBUG_LOCAL_LEVEL(-1)
   sf /= norm;
   VFN_DEBUG_MESSAGE("GlobalScatteringPower::GetScatteringFactor():"<<mName<<":End.",10)
   delete pData;
   return sf;
}
REAL GlobalScatteringPower::GetForwardScatteringFactor(const RadiationType type) const
{
   REAL sf=0;
   const ScatteringComponentList *pList=&(mpZScatterer->GetScatteringComponentList());
   for(int i=0;i<pList->GetNbComponent();i++)
      sf += (*pList)(i).mpScattPow->GetForwardScatteringFactor(type);
   return sf;
}
CrystVector_REAL GlobalScatteringPower::
   GetTemperatureFactor(const ScatteringData &data,
                        const int spgSymPosIndex) const
{
   VFN_DEBUG_MESSAGE("GlobalScatteringPower::GetTemperatureFactor(data):"<<mName,5)
   CrystVector_REAL temp(data.GetNbRefl());
   temp=1.;
   return temp;
}
CrystMatrix_REAL GlobalScatteringPower::
   GetResonantScattFactReal(const ScatteringData &data,
                            const int spgSymPosIndex) const
{
   VFN_DEBUG_MESSAGE("GlobalScatteringPower::GetResonantScattFactReal(data):"<<mName,5)
   CrystMatrix_REAL res(1,1);
   res=0.;
   return res;
}
CrystMatrix_REAL GlobalScatteringPower::
   GetResonantScattFactImag(const ScatteringData &data,
                            const int spgSymPosIndex) const
{
   VFN_DEBUG_MESSAGE("GlobalScatteringPower::GetResonantScattFactImag(data):"<<mName,5)
   CrystMatrix_REAL res(1,1);
   res=0.;
   return res;
}
REAL GlobalScatteringPower::GetRadius()const
{
   //:TODO:
   return 3.;
}
void GlobalScatteringPower::InitRefParList()
{
   //nothing to do, nothing to refine 8-))
}

}//namespace
