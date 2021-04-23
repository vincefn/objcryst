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
*  source file LibCryst++ Crystal class
*
*/

#include <cmath>
#include <set>
#include <vector>
#include <typeinfo>
#include <boost/format.hpp>

#include "cctbx/sgtbx/space_group.h"

#include "ObjCryst/ObjCryst/Crystal.h"
#include "ObjCryst/ObjCryst/Molecule.h"
#include "ObjCryst/ObjCryst/Atom.h"

#include "ObjCryst/Quirks/VFNStreamFormat.h" //simple formatting of integers, REALs..
#include "ObjCryst/Quirks/VFNDebug.h"

#ifdef __WX__CRYST__
   #include "ObjCryst/wxCryst/wxCrystal.h"
#endif

#include <fstream>
#include <iomanip>

namespace ObjCryst
{
const RefParType *gpRefParTypeCrystal=0;
long NiftyStaticGlobalObjectsInitializer_Crystal::mCount=0;
////////////////////////////////////////////////////////////////////////
//
//    CRYSTAL : the crystal (Unit cell, spaceGroup, scatterers)
//
////////////////////////////////////////////////////////////////////////
ObjRegistry<Crystal> gCrystalRegistry("List of all Crystals");

Crystal::Crystal():
mScattererRegistry("List of Crystal Scatterers"),
mBumpMergeCost(0.0),mBumpMergeScale(1.0),
mDistTableMaxDistance(1.0),
mScatteringPowerRegistry("List of Crystal ScatteringPowers"),
mBondValenceCost(0.0),mBondValenceCostScale(1.0),mDeleteSubObjInDestructor(1),
mInterMolDistCostScale(1.0),mInterMolDistCost(0.0)
{
   VFN_DEBUG_MESSAGE("Crystal::Crystal()",10)
   this->InitOptions();
   this->Init(10,11,12,M_PI/2+.1,M_PI/2+.2,M_PI/2+.3,"P 1","");
   gCrystalRegistry.Register(*this);
   gTopRefinableObjRegistry.Register(*this);
   mClockMaster.AddChild(mLatticeClock);
   mClockMaster.AddChild(this->mScattererRegistry.GetRegistryClock());
   mClockMaster.AddChild(this->mScatteringPowerRegistry.GetRegistryClock());
}

Crystal::Crystal(const REAL a, const REAL b, const REAL c, const string &SpaceGroupId):
mScattererRegistry("List of Crystal Scatterers"),
mBumpMergeCost(0.0),mBumpMergeScale(1.0),
mDistTableMaxDistance(1.0),
mScatteringPowerRegistry("List of Crystal ScatteringPowers"),
mBondValenceCost(0.0),mBondValenceCostScale(1.0),mDeleteSubObjInDestructor(1),
mInterMolDistCostScale(1.0),mInterMolDistCost(0.0)
{
   VFN_DEBUG_MESSAGE("Crystal::Crystal(a,b,c,Sg)",10)
   this->Init(a,b,c,M_PI/2,M_PI/2,M_PI/2,SpaceGroupId,"");
   this->InitOptions();
   gCrystalRegistry.Register(*this);
   gTopRefinableObjRegistry.Register(*this);
   mClockMaster.AddChild(mLatticeClock);
   mClockMaster.AddChild(this->mScattererRegistry.GetRegistryClock());
   mClockMaster.AddChild(this->mScatteringPowerRegistry.GetRegistryClock());
}

Crystal::Crystal(const REAL a, const REAL b, const REAL c, const REAL alpha,
              const REAL beta, const REAL gamma,const string &SpaceGroupId):
mScattererRegistry("List of Crystal Scatterers"),
mBumpMergeCost(0.0),mBumpMergeScale(1.0),
mDistTableMaxDistance(1.0),
mScatteringPowerRegistry("List of Crystal ScatteringPowers"),
mBondValenceCost(0.0),mBondValenceCostScale(1.0),mDeleteSubObjInDestructor(1),
mInterMolDistCostScale(1.0),mInterMolDistCost(0.0)
{
   VFN_DEBUG_MESSAGE("Crystal::Crystal(a,b,c,alpha,beta,gamma,Sg)",10)
   this->Init(a,b,c,alpha,beta,gamma,SpaceGroupId,"");
   this->InitOptions();
   gCrystalRegistry.Register(*this);
   gTopRefinableObjRegistry.Register(*this);
   mClockMaster.AddChild(mLatticeClock);
   mClockMaster.AddChild(this->mScattererRegistry.GetRegistryClock());
   mClockMaster.AddChild(this->mScatteringPowerRegistry.GetRegistryClock());
}

Crystal::Crystal(const Crystal &old):
mScattererRegistry("List of Crystal Scatterers"),
mBumpMergeCost(0.0),mBumpMergeScale(1.0),
mDistTableMaxDistance(1.0),
mScatteringPowerRegistry("List of Crystal ScatteringPowers"),
mBondValenceCost(0.0),mBondValenceCostScale(1.0),mDeleteSubObjInDestructor(1),
mInterMolDistCostScale(1.0),mInterMolDistCost(0.0)
{
   VFN_DEBUG_MESSAGE("Crystal::Crystal()",10)
   // Only create a default crystal, then copy old using XML
   this->InitOptions();
   this->Init(10,11,12,M_PI/2+.1,M_PI/2+.2,M_PI/2+.3,"P 1","");
   gCrystalRegistry.Register(*this);
   gTopRefinableObjRegistry.Register(*this);
   mClockMaster.AddChild(mLatticeClock);
   mClockMaster.AddChild(this->mScattererRegistry.GetRegistryClock());
   mClockMaster.AddChild(this->mScatteringPowerRegistry.GetRegistryClock());

   stringstream sst;
   old.XMLOutput(sst);
   XMLCrystTag tag(sst);
   this->XMLInput(sst,tag);
}

Crystal::~Crystal()
{
   VFN_DEBUG_ENTRY("Crystal::~Crystal()",5)
   for(long i=0;i<mScattererRegistry.GetNb();i++)
   {
      VFN_DEBUG_MESSAGE("Crystal::~Crystal(&scatt):1:"<<i,5)
      this->RemoveSubRefObj(mScattererRegistry.GetObj(i));
      mScattererRegistry.GetObj(i).DeRegisterClient(*this);
   }
   if( mDeleteSubObjInDestructor )
   {
       mScattererRegistry.DeleteAll();
   }
   else
   {
       mScattererRegistry.DeRegisterAll();
   }
   for(long i=0;i<mScatteringPowerRegistry.GetNb();i++)
   {
      VFN_DEBUG_MESSAGE("Crystal::~Crystal(&scatt):2:"<<i,5)
      this->RemoveSubRefObj(mScatteringPowerRegistry.GetObj(i));
      mScatteringPowerRegistry.GetObj(i).DeRegisterClient(*this);
      // :TODO: check if it is not used by another Crystal (forbidden!)
   }
   if( mDeleteSubObjInDestructor )
   {
       mScatteringPowerRegistry.DeleteAll();
   }
   else
   {
       mScatteringPowerRegistry.DeRegisterAll();
   }
   gCrystalRegistry.DeRegister(*this);
   gTopRefinableObjRegistry.DeRegister(*this);
   VFN_DEBUG_EXIT("Crystal::~Crystal()",5)
}

const string& Crystal::GetClassName() const
{
   const static string className="Crystal";
   return className;
}

void Crystal::AddScatterer(Scatterer *scatt)
{
   VFN_DEBUG_ENTRY("Crystal::AddScatterer(&scatt)",5)
   mScattererRegistry.Register(*scatt);
   scatt->RegisterClient(*this);
   this->AddSubRefObj(*scatt);
   scatt->SetCrystal(*this);
   mClockScattererList.Click();
   VFN_DEBUG_EXIT("Crystal::AddScatterer(&scatt):Finished",5)
}

void Crystal::RemoveScatterer(Scatterer *scatt, const bool del)
{
   VFN_DEBUG_MESSAGE("Crystal::RemoveScatterer(&scatt)",5)
   mScattererRegistry.DeRegister(*scatt);
   scatt->DeRegisterClient(*this);
   this->RemoveSubRefObj(*scatt);
   if(del) delete scatt;
   mClockScattererList.Click();
   VFN_DEBUG_MESSAGE("Crystal::RemoveScatterer(&scatt):Finished",5)
}

long Crystal::GetNbScatterer()const {return mScattererRegistry.GetNb();}

Scatterer& Crystal::GetScatt(const string &scattName)
{
   return mScattererRegistry.GetObj(scattName);
}

const Scatterer& Crystal::GetScatt(const string &scattName) const
{
   return mScattererRegistry.GetObj(scattName);
}

Scatterer& Crystal::GetScatt(const long scattIndex)
{
   return mScattererRegistry.GetObj(scattIndex);
}

const Scatterer& Crystal::GetScatt(const long scattIndex) const
{
   return mScattererRegistry.GetObj(scattIndex);
}

ObjRegistry<Scatterer>& Crystal::GetScattererRegistry() {return mScattererRegistry;}

const ObjRegistry<Scatterer>& Crystal::GetScattererRegistry() const {return mScattererRegistry;}

ObjRegistry<ScatteringPower>& Crystal::GetScatteringPowerRegistry()
{return mScatteringPowerRegistry;}
const ObjRegistry<ScatteringPower>& Crystal::GetScatteringPowerRegistry() const
{return mScatteringPowerRegistry;}

void Crystal::AddScatteringPower(ScatteringPower *scattPow)
{
   mScatteringPowerRegistry.Register(*scattPow);
   scattPow->RegisterClient(*this);//:TODO: Should register as (unique) 'owner'.
   this->AddSubRefObj(*scattPow);
   mClockMaster.AddChild(scattPow->GetClockMaster());
   mClockMaster.AddChild(scattPow->GetMaximumLikelihoodParClock());
   mMasterClockScatteringPower.AddChild(scattPow->GetClockMaster());
}

void Crystal::RemoveScatteringPower(ScatteringPower *scattPow, const bool del)
{
   VFN_DEBUG_ENTRY("Crystal::RemoveScatteringPower()",2)
   mScatteringPowerRegistry.DeRegister(*scattPow);
   this->RemoveSubRefObj(*scattPow);
   mClockMaster.RemoveChild(scattPow->GetClockMaster());
   mClockMaster.RemoveChild(scattPow->GetMaximumLikelihoodParClock());
   mMasterClockScatteringPower.RemoveChild(scattPow->GetClockMaster());
   if(del) delete scattPow;

   for(Crystal::VBumpMergePar::iterator pos=mvBumpMergePar.begin();pos!=mvBumpMergePar.end();)
   {
      if((pos->first.first==scattPow)||(pos->first.second==scattPow))
      {
         mvBumpMergePar.erase(pos++);
         mBumpMergeParClock.Click();
      }
      else ++pos;// See Josuttis Std C++ Lib p.205 for safe method
   }

   for(map<pair<const ScatteringPower*,const ScatteringPower*>, REAL>::iterator
         pos=mvBondValenceRo.begin();pos!=mvBondValenceRo.end();)
   {
      if((pos->first.first==scattPow)||(pos->first.second==scattPow))
      {
         mvBondValenceRo.erase(pos++);
         mBondValenceParClock.Click();
      }
      else ++pos;// See Josuttis Std C++ Lib p.205 for safe method
   }
   VFN_DEBUG_EXIT("Crystal::RemoveScatteringPower()",2)
}

ScatteringPower& Crystal::GetScatteringPower(const string &name)
{
   return mScatteringPowerRegistry.GetObj(name);
}

const ScatteringPower& Crystal::GetScatteringPower(const string &name)const
{
   return mScatteringPowerRegistry.GetObj(name);
}

const RefinableObjClock& Crystal::GetMasterClockScatteringPower()const
{ return mMasterClockScatteringPower;}


const ScatteringComponentList& Crystal::GetScatteringComponentList()const
{
   if(mClockScattCompList>mClockMaster) return mScattCompList;
   bool update=false;
   if(mClockScattCompList<this->GetClockLatticePar()) update=true;
   if(update==false)
      for(long i=0;i<mScattererRegistry.GetNb();i++)
      {
         //mClockScattCompList.Print();
         //this->GetScatt(i).GetClockScatterer().Print();
         if(mClockScattCompList<this->GetScatt(i).GetClockScatterer()) {update=true;break;}
      }
   if(true==update)
   {
      VFN_DEBUG_MESSAGE("Crystal::GetScatteringComponentList()",2)
      TAU_PROFILE("Crystal::GetScatteringComponentList()","ScattCompList& ()",TAU_DEFAULT);
      mScattCompList.Reset();
      for(long i=0;i<mScattererRegistry.GetNb();i++)
      {
         //this->GetScatt(i).GetScatteringComponentList().Print();
         mScattCompList += this->GetScatt(i).GetScatteringComponentList();
      }

      //:KLUDGE: this must be *before* calling CalcDynPopCorr() to avoid an infinite loop..
      mClockScattCompList.Click();

      if(1==mUseDynPopCorr.GetChoice())
         this->CalcDynPopCorr(1.,.5); else this->ResetDynPopCorr();
      VFN_DEBUG_MESSAGE("Crystal::GetScatteringComponentList():End",2)
   }
   #ifdef __DEBUG__
   if(gVFNDebugMessageLevel<2) mScattCompList.Print();
   #endif
   return mScattCompList;
}

const RefinableObjClock& Crystal::GetClockScattCompList()const
{
   return mClockScattCompList;
}

void Crystal::Print(ostream &os)const
{
   VFN_DEBUG_MESSAGE("Crystal::Print()",5)
   this->UnitCell::Print(os);

   this->GetScatteringComponentList();
   this->CalcBondValenceSum();

   os << "List of scattering components (atoms): " << mScattCompList.GetNbComponent() << endl ;

   long k=0;
   std::map<long, REAL>::const_iterator posBV;
   for(int i=0;i<mScattererRegistry.GetNb();i++)
   {
      //mpScatterrer[i]->Print();
      const ScatteringComponentList list=this->GetScatt(i).GetScatteringComponentList();
      for(int j=0;j<list.GetNbComponent();j++)
      {
         os   << FormatString(this->GetScatt(i).GetComponentName(j),16) << " at : "
              << FormatFloat(list(j).mX,7,4)
              << FormatFloat(list(j).mY,7,4)
              << FormatFloat(list(j).mZ,7,4)
              << ", Occup=" << FormatFloat(list(j).mOccupancy,6,4)
              << " * " << FormatFloat(mScattCompList(k).mDynPopCorr,6,4);
         // Check for dummy atoms
         if( NULL != list(j).mpScattPow )
         {
           os << ", ScattPow:" << FormatString(list(j).mpScattPow->GetName(),16)
              << ", Biso=" << FormatFloat(list(j).mpScattPow->GetBiso());
         }
         else
         {
           os   << ", ScattPow: dummy";
         }
         // Check for dummy atoms
         if( NULL != this->mScattCompList(k).mpScattPow )
         {
             posBV=this->mvBondValenceCalc.find(k);
             if(posBV!=this->mvBondValenceCalc.end())
                os <<": Valence="<<posBV->second<<" (expected="
                   <<this->mScattCompList(k).mpScattPow->GetFormalCharge()<<")";
         }
         os << endl;
         k++;
      }
   }
   os <<endl
      << "Occupancy = occ * dyn, where:"<<endl
      << "        - occ is the 'real' occupancy"<< endl
      << "        - dyn is the dynamical occupancy correction, indicating  either"<< endl
      << "          an atom on a special position, or several identical atoms "<< endl
      << "          overlapping (dyn=0.5 -> atom on a symetry plane / 2fold axis.."<< endl
      << "                               -> OR 2 atoms strictly overlapping)"<< endl
      <<endl;
   REAL nbAtoms=0;
   const long genMult=this->GetSpaceGroup().GetNbSymmetrics();
   for(int i=0;i<mScattCompList.GetNbComponent();i++)
      nbAtoms += genMult * mScattCompList(i).mOccupancy * mScattCompList(i).mDynPopCorr;
   os << " Total number of components (atoms) in one unit cell : " << nbAtoms<<endl
      << " Chemical formula: "<< this->GetFormula() << endl
      << " Weight: "<< this->GetWeight()<< " g/mol" << endl;

   VFN_DEBUG_MESSAGE("Crystal::Print():End",5)
}

std::string Crystal::GetFormula() const
{
   this->GetScatteringComponentList();
   if(mScattCompList.GetNbComponent() == 0) return "";
   std::map<std::string,float> velts;
   for(unsigned int i=0; i<mScattCompList.GetNbComponent(); ++i)
   {
      const ScatteringComponent* psi = &mScattCompList(i);
      if(psi->mpScattPow == 0) continue;
      if(psi->mpScattPow->GetClassName().compare("ScatteringPowerAtom")!=0) continue;
      const ScatteringPowerAtom *pat=dynamic_cast<const ScatteringPowerAtom*>(psi->mpScattPow);
      string p=pat->GetSymbol();
      if(velts.count(p)==0)
         velts[p] = psi->mOccupancy * psi->mDynPopCorr ;
      else velts[p] += psi->mOccupancy * psi->mDynPopCorr;
   }
   stringstream s;
   s<<std::setprecision(2);
   for(std::map<std::string,float>::const_iterator pos=velts.begin();pos!=velts.end();++pos)
   {
      if(pos!=velts.begin()) s<<" ";
      float nb=pos->second;
      if(abs(round(nb)-nb)<0.005)
      {
         if(int(round(nb))==1) s<<pos->first;
         else s<<pos->first<<int(round(nb));
      }
      else s<<pos->first<<nb;
   }
   return s.str();
}


REAL Crystal::GetWeight() const
{
   this->GetScatteringComponentList();
   if(mScattCompList.GetNbComponent() == 0) return 0;
   REAL w=0;
   for(unsigned int i=0; i<mScattCompList.GetNbComponent(); ++i)
   {
      const ScatteringComponent* psi = &mScattCompList(i);
      if(psi->mpScattPow == 0) continue;
      if(psi->mpScattPow->GetClassName().compare("ScatteringPowerAtom")!=0) continue;
      const ScatteringPowerAtom *pat=dynamic_cast<const ScatteringPowerAtom*>(psi->mpScattPow);
      w += pat->GetAtomicWeight() * psi->mOccupancy * psi->mDynPopCorr ;
   }
   return w;
}


CrystMatrix_REAL Crystal::GetMinDistanceTable(const REAL minDistance) const
{
   VFN_DEBUG_MESSAGE("Crystal::MinDistanceTable()",5)
   this->CalcDistTable(true);
   const long nbComponent=mScattCompList.GetNbComponent();

   CrystMatrix_REAL minDistTable(nbComponent,nbComponent);
   REAL dist;
   REAL tmp;
   REAL min=minDistance*minDistance;
   if(minDistance<0) min = -1.;// no min distance
   minDistTable=10000.;
   for(int i=0;i<nbComponent;i++)
   {
      //VFN_DEBUG_MESSAGE("Crystal::MinDistanceTable():comp="<<i,10)
      std::vector<Crystal::Neighbour>::const_iterator pos;
      for(pos=mvDistTableSq[i].mvNeighbour.begin();
          pos<mvDistTableSq[i].mvNeighbour.end();pos++)
      {
         tmp=pos->mDist2;
         dist=minDistTable(i,pos->mNeighbourIndex);
         if(  (tmp<dist)
            && ((tmp>min) || (  (mvDistTableSq[i].mIndex !=pos->mNeighbourIndex)
                              &&(mvDistTableSq[i].mUniquePosSymmetryIndex
                                 !=pos->mNeighbourSymmetryIndex))))
            minDistTable(i,pos->mNeighbourIndex)=tmp;
      }
   }
   for(int i=0;i<nbComponent;i++)
   {
      for(int j=0;j<=i;j++)
      {
         if(9999.>minDistTable(i,j)) minDistTable(i,j)=sqrt(minDistTable(i,j));
         else minDistTable(i,j)=-1;
         minDistTable(j,i)=minDistTable(i,j);
      }
   }
   VFN_DEBUG_MESSAGE("Crystal::MinDistanceTable():End",3)
   return minDistTable;
}

void Crystal::PrintMinDistanceTable(const REAL minDistance,ostream &os) const
{
   VFN_DEBUG_MESSAGE("Crystal::PrintMinDistanceTable()",5)
   CrystMatrix_REAL minDistTable;
   minDistTable=this->GetMinDistanceTable(minDistance);
   VFN_DEBUG_MESSAGE("Crystal::PrintMinDistanceTable():0",5)
   os << "Table of minimal distances between all components (atoms)"<<endl;
   os << "               ";
   for(long i=0;i<mScattererRegistry.GetNb();i++)
   {
      VFN_DEBUG_MESSAGE("Crystal::PrintMinDistanceTable()1:Scatt:"<<i,3)
      for(long j=0;j<this->GetScatt(i).GetNbComponent();j++)
         os << FormatString(this->GetScatt(i).GetComponentName(j),7);
   }
   os << endl;
   long l=0;
   const long nbComponent=mScattCompList.GetNbComponent();
   for(long i=0;i<mScattererRegistry.GetNb();i++)
      for(long j=0;j<this->GetScatt(i).GetNbComponent();j++)
      {
         VFN_DEBUG_MESSAGE("Crystal::PrintMinDistanceTable()2:Scatt,comp:"<<i<<","<<j,3)
         os << FormatString(this->GetScatt(i).GetComponentName(j),14);
         for(long k=0;k<nbComponent;k++)
         {
            if(minDistTable(l,k)>0) os << FormatFloat(minDistTable(l,k),6,3) ;
            else os<<"  >10  ";
         }
         os << endl;
         l++;
      }
   VFN_DEBUG_MESSAGE("Crystal::PrintMinDistanceTable():End",3)
}

ostream& Crystal::POVRayDescription(ostream &os,const CrystalPOVRayOptions &options)const
{
   VFN_DEBUG_MESSAGE("Crystal::POVRayDescription(os,bool)",5)
   os <<"/////////////////////// MACROS////////////////////"<<endl;

   os << "#macro ObjCrystAtom(atomx,atomy,atomz,atomr,atomc,occ,atten)"
      << "   sphere"<<endl
      << "   { <atomx,atomy,atomz>,atomr"<<endl
      << "      finish {ambient 0.5*occ*atten diffuse 0.4*occ*atten phong occ*atten specular 0.2*occ*atten "
      << "roughness 0.02 metallic reflection 0.0}"<<endl
      << "      pigment { colour atomc transmit (1-occ*atten)}"<<endl
      << "      no_shadow"<<endl
      << "   }"<<endl
      << "#end"<<endl<<endl;

   os << "#macro ObjCrystBond(x1,y1,z1,x2,y2,z2,bondradius,bondColour,occ,atten)"<<endl
      << "   cylinder"<<endl
      << "   {  <x1,y1,z1>,"<<endl
      << "      <x2,y2,z2>,"<<endl
      << "      bondradius"<<endl
      << "      finish {ambient 0.5*occ*atten diffuse 0.4*occ*atten phong occ*atten specular 0.2*occ*atten "
      << "roughness 0.02 metallic reflection 0.0}"<<endl
      << "      pigment { colour bondColour transmit (1-occ*atten)}"<<endl
      << "      no_shadow"<<endl
      << "   }"<<endl
      << "#end"<<endl<<endl;

   os <<"//////////// Crystal Unit Cell /////////////////" <<endl;
   REAL x=1,y=1,z=1;
   this->FractionalToOrthonormalCoords(x,y,z);
   os << "   //box{ <0,0,0>, <"<< x << "," << y << "," << z << ">" <<endl;
   os << "   //      pigment {colour rgbf<1,1,1,0.9>}" << endl;
   os << "   //      hollow" << endl;
   os << "   //}" <<endl<<endl;

   #define UNITCELL_EDGE(X0,Y0,Z0,X1,Y1,Z1)\
   {\
      REAL x0=X0,y0=Y0,z0=Z0,x1=X1,y1=Y1,z1=Z1;\
      this->FractionalToOrthonormalCoords(x0,y0,z0);\
      this->FractionalToOrthonormalCoords(x1,y1,z1);\
      os << "    ObjCrystBond("\
         <<x0<<","<<y0<<","<<z0<< ","\
         <<x1<<","<<y1<<","<<z1<< ","\
         << "0.02,rgb<1.0,1.0,1.0>,1.0,1.0)"<<endl;\
   }

   UNITCELL_EDGE(0,0,0,1,0,0)
   UNITCELL_EDGE(0,0,0,0,1,0)
   UNITCELL_EDGE(0,0,0,0,0,1)

   UNITCELL_EDGE(1,1,1,0,1,1)
   UNITCELL_EDGE(1,1,1,1,0,1)
   UNITCELL_EDGE(1,1,1,1,1,0)

   UNITCELL_EDGE(1,0,0,1,1,0)
   UNITCELL_EDGE(1,0,0,1,0,1)

   UNITCELL_EDGE(0,1,0,1,1,0)
   UNITCELL_EDGE(0,1,0,0,1,1)

   UNITCELL_EDGE(0,0,1,1,0,1)
   UNITCELL_EDGE(0,0,1,0,1,1)

   os <<endl<<"/////////////// GLOBAL DECLARATIONS FOR ATOMS & BONDS ///////"<<endl;
   os << "// Atom colours"<<endl;
   for(int i=0;i<mScatteringPowerRegistry.GetNb();i++)
   {
      const float r=mScatteringPowerRegistry.GetObj(i).GetColourRGB()[0];
      const float g=mScatteringPowerRegistry.GetObj(i).GetColourRGB()[1];
      const float b=mScatteringPowerRegistry.GetObj(i).GetColourRGB()[2];
      os << "   #declare colour_"<< mScatteringPowerRegistry.GetObj(i).GetName()
         <<"= rgb <"<< r<<","<<g<<","<<b<<">;"<< endl;
   }
   os << "// Bond colours"<<endl
      << "   #declare colour_freebond   = rgb <0.7,0.7,0.7>;"<< endl
      << "   #declare colour_nonfreebond= rgb <0.3,0.3,0.3>;"<< endl;

   os<<endl;
   os << "/////////////// SCATTERERS ///////"<<endl;
   for(int i=0;i<mScattererRegistry.GetNb();i++)
      this->GetScatt(i).POVRayDescription(os,options) ;
   return os;
}

#ifdef OBJCRYST_GL
void Crystal::GLInitDisplayList(const bool onlyIndependentAtoms,
                                const REAL xMin,const REAL xMax,
                                const REAL yMin,const REAL yMax,
                                const REAL zMin,const REAL zMax,
                                const bool displayNames,
                                const bool hideHydrogens,
                                const REAL fadeDistance,
                                const bool fullMoleculeInLimits)const
{
   VFN_DEBUG_ENTRY("Crystal::GLInitDisplayList()",5)
      REAL en=1;// if -1, display enantiomeric structure
      if(mDisplayEnantiomer.GetChoice()==1) en=-1;

      //Center of displayed unit
         REAL xc=(xMin+xMax)/2.;
         REAL yc=(yMin+yMax)/2.;
         REAL zc=(zMin+zMax)/2.;
      if(false==displayNames)
      {
         //Describe Unit Cell
            REAL x111= 1.;
            REAL y111= 1.;
            REAL z111= 1.;
            this->FractionalToOrthonormalCoords(x111,y111,z111);
            REAL x110= 1.;
            REAL y110= 1.;
            REAL z110= 0.;
            this->FractionalToOrthonormalCoords(x110,y110,z110);
            REAL x101= 1.;
            REAL y101= 0.;
            REAL z101= 1.;
            this->FractionalToOrthonormalCoords(x101,y101,z101);
            REAL x100= 1.;
            REAL y100= 0.;
            REAL z100= 0.;
            this->FractionalToOrthonormalCoords(x100,y100,z100);
            REAL x011= 0.;
            REAL y011= 1.;
            REAL z011= 1.;
            this->FractionalToOrthonormalCoords(x011,y011,z011);
            REAL x010= 0.;
            REAL y010= 1.;
            REAL z010= 0.;
            this->FractionalToOrthonormalCoords(x010,y010,z010);
            REAL x001= 0.;
            REAL y001= 0.;
            REAL z001= 1.;
            this->FractionalToOrthonormalCoords(x001,y001,z001);
            REAL x000= 0.;
            REAL y000= 0.;
            REAL z000= 0.;
            this->FractionalToOrthonormalCoords(x000,y000,z000);
            REAL xM= 0.5;
            REAL yM= 0.5;
            REAL zM= 0.5;
            this->FractionalToOrthonormalCoords(xM,yM,zM);
            xM*=2;yM*=2;zM*=2;
         glPushMatrix();
         //Add Axis & axis names
            const GLfloat colour0 [] = {0.00, 0.00, 0.00, 0.00};
            const GLfloat colour1 [] = {0.50, 0.50, 0.50, 1.00};
            const GLfloat colour2 [] = {1.00, 1.00, 1.00, 1.00};
            glMaterialfv(GL_FRONT, GL_AMBIENT,   colour1);
            glMaterialfv(GL_FRONT, GL_DIFFUSE,   colour0);
            glMaterialfv(GL_FRONT, GL_SPECULAR,  colour0);
            glMaterialfv(GL_FRONT, GL_EMISSION,  colour1);
            glMaterialfv(GL_FRONT, GL_SHININESS, colour0);
            REAL x,y,z;
            x=1.2-xc;y=-yc;z=-zc;
            this->FractionalToOrthonormalCoords(x,y,z);
            glRasterPos3f(en*x,y,z);
            crystGLPrint("a");

            x=-xc;y=1.2-yc;z=-zc;
            this->FractionalToOrthonormalCoords(x,y,z);
            glRasterPos3f(en*x,y,z);
            crystGLPrint("b");

            x=-xc;y=-yc;z=1.2-zc;
            this->FractionalToOrthonormalCoords(x,y,z);
            glRasterPos3f(en*x,y,z);
            crystGLPrint("c");
         // Cell
            glMaterialfv(GL_FRONT, GL_AMBIENT,   colour1);
            glMaterialfv(GL_FRONT, GL_DIFFUSE,   colour1);
            glMaterialfv(GL_FRONT, GL_SPECULAR,  colour1);
            glMaterialfv(GL_FRONT, GL_EMISSION,  colour0);
            glMaterialfv(GL_FRONT, GL_SHININESS, colour0);
            this->FractionalToOrthonormalCoords(xc,yc,zc);
            glTranslatef(-xc*en, -yc, -zc);
            glBegin(GL_LINES);
               //top
               glNormal3f((x110+x010-xM)*en,y110+y010-yM,z110+z010-zM);
               glVertex3f(    x110*en,    y110,    z110);
               glVertex3f(    x010*en,    y010,    z010);

               glNormal3f((x011+x010-xM)*en,y011+y010-yM,z011+z010-zM);
               glVertex3f(    x010*en,    y010,    z010);
               glVertex3f(    x011*en,    y011,    z011);

               glNormal3f((x011+x111-xM)*en,y011+y111-yM,z011+z111-zM);
               glVertex3f(    x011*en,    y011,    z011);
               glVertex3f(    x111*en,    y111,    z111);

               glNormal3f((x110+x111-xM)*en,y110+y111-yM,z110+z111-zM);
               glVertex3f(    x111*en,    y111,    z111);
               glVertex3f(    x110*en,    y110,    z110);
               //bottom
               glNormal3f((x101+x001-xM)*en,y101+y001-yM,z101+z001-zM);
               glVertex3f(    x101*en,    y101,    z101);
               glVertex3f(    x001*en,    y001,    z001);

               glNormal3f((x000+x001-xM)*en,y000+y001-yM,z000+z001-zM);
               glVertex3f(    x001*en,    y001,    z001);
               glVertex3f(    x000*en,    y000,    z000);

               glNormal3f((x000+x100-xM)*en,y000+y100-yM,z000+z100-zM);
               glVertex3f(    x000*en,    y000,    z000);
               glVertex3f(    x100*en,    y100,    z100);

               glNormal3f((x101+x100-xM)*en,y101+y100-yM,z101+z100-zM);
               glVertex3f(    x100*en,    y100,    z100);
               glVertex3f(    x101*en,    y101,    z101);
               //sides
               glNormal3f((x101+x111-xM)*en,y101+y111-yM,z101+z111-zM);
               glVertex3f(    x101*en,    y101,    z101);
               glVertex3f(    x111*en,    y111,    z111);

               glNormal3f((x001+x011-xM)*en,y001+y011-yM,z001+z011-zM);
               glVertex3f(    x001*en,    y001,     z001);
               glVertex3f(    x011*en,    y011,     z011);

               glNormal3f((x000+x010-xM)*en,y000+y010-yM,z000+z010-zM);
               glVertex3f(    x000*en,    y000,    z000);
               glVertex3f(    x010*en,    y010,    z010);

               glNormal3f((x100+x110-xM)*en,y100+y110-yM,z100+z110-zM);
               glVertex3f(    x100*en,    y100,    z100);
               glVertex3f(    x110*en,    y110,    z110);
            glEnd();
         glPopMatrix();
      }

      //Describe all Scatterers
      VFN_DEBUG_MESSAGE("Crystal::GLView(bool):Scatterers...",5)
      glPushMatrix();
         if(displayNames) this->FractionalToOrthonormalCoords(xc,yc,zc);
         glTranslatef(-xc*en, -yc, -zc);
         glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
         {
            bool displayEnantiomer=false;
            if(mDisplayEnantiomer.GetChoice()==1) displayEnantiomer=true;
            for(int i=0;i<mScattererRegistry.GetNb();i++)
               this->GetScatt(i).GLInitDisplayList(onlyIndependentAtoms,
                                                   xMin,xMax,yMin,yMax,zMin,zMax,
                                                   displayEnantiomer,displayNames,hideHydrogens,fadeDistance,fullMoleculeInLimits);
         }
      glPopMatrix();
   VFN_DEBUG_EXIT("Crystal::GLInitDisplayList(bool)",5)
}
#endif  // OBJCRYST_GL

void Crystal::CalcDynPopCorr(const REAL overlapDist, const REAL mergeDist) const
{
   VFN_DEBUG_ENTRY("Crystal::CalcDynPopCorr(REAL)",4)
   TAU_PROFILE("Crystal::CalcDynPopCorr()","void (REAL)",TAU_DEFAULT);

   this->CalcDistTable(true);
   if(mClockDynPopCorr>mDistTableClock) return;

   const long nbComponent=mScattCompList.GetNbComponent();
   const int nbSymmetrics=this->GetSpaceGroup().GetNbSymmetrics();
   CrystVector_REAL neighborsDist(nbComponent*nbSymmetrics);
   long nbNeighbors=0;
   REAL corr;

   int atomicNumber;
   const REAL overlapDistSq=overlapDist*overlapDist;
   REAL dist;
   for(long i=0;i<nbComponent;i++)
   {
      VFN_DEBUG_MESSAGE("Crystal::CalcDynPopCorr(): Component:"<<i,0)
      if(0==mScattCompList(i).mpScattPow)
      {
         mScattCompList(i).mDynPopCorr=1.;
         continue;
      }
      atomicNumber=mScattCompList(i).mpScattPow->GetDynPopCorrIndex();
      nbNeighbors=0;
      std::vector<Crystal::Neighbour>::const_iterator pos;
      for(pos=mvDistTableSq[i].mvNeighbour.begin();
          pos<mvDistTableSq[i].mvNeighbour.end();pos++)
      {
         VFN_DEBUG_MESSAGE("Crystal::CalcDynPopCorr(): Component:"<<i<<"Neighbour:"<<pos->mNeighbourIndex,0)
         if(0==mScattCompList(pos->mNeighbourIndex).mpScattPow)continue;
         if(atomicNumber==mScattCompList(pos->mNeighbourIndex).mpScattPow->GetDynPopCorrIndex())
         {
            if(overlapDistSq > pos->mDist2)
            {
               //resizing can be necessary if the unit cell is small, so that an atom two unit cells away is
               //still considered a neighbor...
               if(nbNeighbors==neighborsDist.numElements()) neighborsDist.resizeAndPreserve(nbNeighbors+20);
               neighborsDist(nbNeighbors++)=sqrt(pos->mDist2);
            }
         }
      }
      corr=0.;
      for(long j=0;j<nbNeighbors;j++)
      {
         dist=neighborsDist(j)-mergeDist;
         if(dist<0.) dist=0.;
         corr += fabs((overlapDist-dist-mergeDist)/(overlapDist-mergeDist));
      }
      mScattCompList(i).mDynPopCorr=1./(1+corr);
   }
   mClockDynPopCorr.Click();
   VFN_DEBUG_EXIT("Crystal::CalcDynPopCorr(REAL):End.",4)
}

void Crystal::ResetDynPopCorr() const
{
   //:NOTE: This is useless !!
   mClockDynPopCorr.Reset();
   const long nbComponent=mScattCompList.GetNbComponent();
   for(long i=0;i<nbComponent;i++) mScattCompList(i).mDynPopCorr=1.;
}

REAL Crystal::GetDynPopCorr(const Scatterer* pscatt, unsigned int component)const
{
   if(this->GetUseDynPopCorr()==0) return 1.0;
   this->GetScatteringComponentList();
   unsigned int j=0;
   for(long i=0;i<mScattererRegistry.GetNb();i++)
   {
      if(pscatt==&(this->GetScatt(i)))
      {
         return mScattCompList(j+component).mDynPopCorr;
      }
      else j+=this->GetScatt(i).GetScatteringComponentList().GetNbComponent();
   }
   // Something is wrong if we got here !
   throw ObjCrystException("Crystal::GetDynPopCorr(): did not find this scatterer !");
   return 1.0;
}

void Crystal::SetUseDynPopCorr(const int b)
{
   VFN_DEBUG_MESSAGE("Crystal::SetUseDynPopCorr()",1)
   mUseDynPopCorr.SetChoice(b);
   mClockDynPopCorr.Reset();
}

int Crystal::GetUseDynPopCorr() const
{
   return mUseDynPopCorr.GetChoice();
}

int Crystal::FindScatterer(const string &scattName)const
{
   VFN_DEBUG_MESSAGE("Crystal::FindScatterer(name)",0)
   for(int i=0;i<this->GetNbScatterer();i++)
   {
      if( mScattererRegistry.GetObj(i).GetName() == scattName) return i;
   }
   throw ObjCrystException("Crystal::FindScatterer(string)\
      Cannot find this scatterer:"+scattName);
}
vector<int> Crystal::FindScatterersInComponentList(const string &scattName)const
{
   vector<int> res;
   for(int i=0;i<mScattCompList.GetNbComponent();i++) {
        if(mScattCompList(i).mName.compare(scattName)==0) {
            res.push_back(i);
        }
   }
   return res;
}

Crystal::BumpMergePar::BumpMergePar():
   mDist2(1.),mCanOverlap(false){}

Crystal::BumpMergePar::BumpMergePar(const REAL dist, const bool canOverlap):
   mDist2(dist*dist),mCanOverlap(canOverlap){}



REAL Crystal::GetBumpMergeCost() const
{
   if(mvBumpMergePar.size()==0) return 0;
   if(mBumpMergeScale<1e-5) return 0;
   this->CalcDistTable(true);
   VFN_DEBUG_ENTRY("Crystal::GetBumpMergeCost()",4)
   if(  (mBumpMergeCostClock>mBumpMergeParClock)
      &&(mBumpMergeCostClock>mDistTableClock)) return mBumpMergeCost*mBumpMergeScale;
   TAU_PROFILE("Crystal::GetBumpMergeCost()","REAL (REAL)",TAU_DEFAULT);

   mBumpMergeCost=0;

   std::vector<NeighbourHood>::const_iterator pos;
   std::vector<Crystal::Neighbour>::const_iterator neigh;
   REAL tmp;
   const ScatteringPower *i1,*i2;
   VBumpMergePar::const_iterator par;
   for(pos=mvDistTableSq.begin();pos<mvDistTableSq.end();pos++)
   {
      i1=mScattCompList(pos->mIndex).mpScattPow;
      for(neigh=pos->mvNeighbour.begin();neigh<pos->mvNeighbour.end();neigh++)
      {
         i2=mScattCompList(neigh->mNeighbourIndex).mpScattPow;
         if(i1<i2) par=mvBumpMergePar.find(std::make_pair(i1,i2));
         else par=mvBumpMergePar.find(std::make_pair(i2,i1));
         if(par==mvBumpMergePar.end()) continue;
         if(neigh->mDist2 > par->second.mDist2) continue;
         if(true==par->second.mCanOverlap)
            tmp = 0.5*sin(M_PI*(1.-sqrt(neigh->mDist2/par->second.mDist2)))/0.1;
         else
            tmp = tan(M_PI*0.49999*(1.-sqrt(neigh->mDist2/par->second.mDist2)))/0.1;
         mBumpMergeCost += tmp*tmp;
      }
   }
   mBumpMergeCost *= this->GetSpaceGroup().GetNbSymmetrics();
   mBumpMergeCostClock.Click();
   VFN_DEBUG_EXIT("Crystal::GetBumpMergeCost():"<<mBumpMergeCost,4)
   return mBumpMergeCost*mBumpMergeScale;
}

void Crystal::SetBumpMergeDistance(const ScatteringPower &scatt1,
                                   const ScatteringPower &scatt2,const REAL dist)
{
   VFN_DEBUG_MESSAGE("Crystal::SetBumpMergeDistance()",5)
   if(&scatt1 == &scatt2) this->SetBumpMergeDistance(scatt1,scatt2,dist,true) ;
   else this->SetBumpMergeDistance(scatt1,scatt2,dist,false);
}

void Crystal::SetBumpMergeDistance(const ScatteringPower &scatt1,
                                   const ScatteringPower &scatt2,const REAL dist,
                                   const bool allowMerge)
{
   VFN_DEBUG_MESSAGE("Crystal::SetBumpMergeDistance("<<scatt1.GetName()<<","<<scatt2.GetName()<<")="<<dist<<","<<allowMerge,3)
   if(&scatt1 < &scatt2)
      mvBumpMergePar[std::make_pair(&scatt1,&scatt2)]=BumpMergePar(dist,allowMerge);
   else
      mvBumpMergePar[std::make_pair(&scatt2,&scatt1)]=BumpMergePar(dist,allowMerge);
   mBumpMergeParClock.Click();
}
void Crystal::RemoveBumpMergeDistance(const ScatteringPower &scatt1,
                                      const ScatteringPower &scatt2)
{
   Crystal::VBumpMergePar::iterator pos;
   if(&scatt1 < &scatt2) pos=mvBumpMergePar.find(make_pair(&scatt1 , &scatt2));
   else pos=mvBumpMergePar.find(make_pair(&scatt2 , &scatt1));
   if(pos!=mvBumpMergePar.end()) mvBumpMergePar.erase(pos);
   mBumpMergeParClock.Click();
}

const Crystal::VBumpMergePar& Crystal::GetBumpMergeParList()const{return mvBumpMergePar;}
Crystal::VBumpMergePar& Crystal::GetBumpMergeParList(){return mvBumpMergePar;}

Crystal::InterMolDistPar::InterMolDistPar():
    mActDist(-1),mDist2(-1),mSig(0.01),mDelta(0.01)
{}
Crystal::InterMolDistPar::InterMolDistPar(const string At1, const string At2, const REAL actualDist, const REAL dist, const REAL sigma, const REAL delta):
    mAt1(At1),mAt2(At2),mActDist(actualDist),mDist2(dist*dist),mSig(sigma),mDelta(delta)
{}
void Crystal::SetNewInterMolDist(const string At1, const string At2, const REAL dist, const REAL sigma, const REAL delta) const
{
    mInterMolDistList.push_back(InterMolDistPar(At1, At2, 0, dist, sigma, delta));
}
int Crystal::GetIntermolDistNb() const
{
    return mInterMolDistList.size();
}
Crystal::InterMolDistPar Crystal::GetIntermolDistPar(int Index) const
{
    if((Index>=mInterMolDistList.size()) || (Index<0)) {
        return InterMolDistPar();
    }
    return mInterMolDistList[Index];
}
REAL Crystal::GetInterMolDistCost() const
{

   //REMOVE THIS PART
   //just for testing - fill something to the intermoldist list
   /*
    if(mInterMolDistList.size()==0) {
       std::cout<<"create testing data for mInterMolDistList\n";
       if((FindScatterersInComponentList("N6").size()!=0) && (FindScatterersInComponentList("N7").size()!=0)) {
           SetNewInterMolDist("N6", "N7", 2.83, 0.1, 0.5);
       }
   }
   std::cout<<"mInterMolDistList.size()=="<<mInterMolDistList.size()<<"\n";
   */
   //REMOVE THIS PART - END

   if(mInterMolDistList.size()==0) return 0;
   if(mInterMolDistCostScale<=0) return 0;

   this->CalcDistTableForInterMolDistCost();
   
   VFN_DEBUG_ENTRY("Crystal::GetInterMolDistCost()",4)
   if(mInterMolDistCostClock>mDistTableForInterMolDistClock) return mInterMolDistCost*mInterMolDistCostScale;
   TAU_PROFILE("Crystal::GetInterMolDistCost()","REAL (REAL)",TAU_DEFAULT);
   
   mInterMolDistCost=0;

   std::cout<<"list of atoms:\n";
   for(int i=0;i<mScattCompList.GetNbComponent();i++) {
       std::cout<<"["<<i<<"]"<<mScattCompList(i).mName<<"\n";
   }
   
   std::vector<NeighbourHood>::const_iterator pos;
   std::vector<Crystal::Neighbour>::const_iterator neigh;
   REAL tmp;
   string i1,i2;  
    
   cout<<"DistanceTable:\n";
   for(pos=imdTable.begin();pos<imdTable.end();pos++) {
       i1=mScattCompList(pos->mIndex).mName;
       cout<<i1<<"\n";       
       for(neigh=pos->mvNeighbour.begin();neigh<pos->mvNeighbour.end();neigh++) {
           i2=mScattCompList(neigh->mNeighbourIndex).mName;
           cout<<"   "<<i2<<" d="<<sqrt(neigh->mDist2)<<"\n";
       }
   }
 
   std::vector<InterMolDistPar>::const_iterator imd;
   //std::vector<NeighbourHood>::const_iterator pos;

   float actdiff=0;
   float mindiff=10000000;
   float bestDist=0;

   //searching all defined mInterMolDistList
   for(imd = mInterMolDistList.begin();imd<mInterMolDistList.end();imd++) {
       actdiff = 0;
       mindiff = 10000000;
       bestDist = mDistTableMaxDistance;    //TODO this is not the correct maxDist!!! It will create jumps in cost function   
       for(pos=imdTable.begin();pos<imdTable.end();pos++)
       {      
          i1=mScattCompList(pos->mIndex).mName;
          //i1=mScattCompList(pos->mIndex).
          
          if(i1 != imd->mAt1) continue;
          std::cout<<"found imd->mAt1 ("<<imd->mAt1<<")\n";
          //cout<<mScattCompList(pos->mIndex).

          for(neigh=pos->mvNeighbour.begin();neigh<pos->mvNeighbour.end();neigh++)
          {             
             //i2=mScattCompList(neigh->mNeighbourIndex).mpScattPow;
             i2=mScattCompList(neigh->mNeighbourIndex).mName;
             if(i2 != imd->mAt2) continue;
             std::cout<<"found imd->mAt2 ("<<imd->mAt2<<")\n";

             if(imd->mDist2 > neigh->mDist2) actdiff = imd->mDist2 - neigh->mDist2;
             else actdiff = neigh->mDist2 - imd->mDist2;
             std::cout<<"actualdiff = "<<actdiff<<"\n";
             if(actdiff<mindiff) {
                 mindiff = actdiff;
                 bestDist = neigh->mDist2;
                 std::cout<<"saving, bestDist="<<bestDist<<"\n";                 
             }         
          }          
       }
       if(bestDist<=0) bestDist = 0;
       else bestDist = sqrt(bestDist);
       float d = imd->mDist2;
       if(d<=0) d = 0;
       else d = sqrt(d);
       std::cout<<"Nx-Ny best dist ("<<d<<") = "<<bestDist<<"\n";

       if((bestDist < (d + imd->mDelta)) && (bestDist > (d - imd->mDelta))) {
           mInterMolDistCost += 0;
       } else if(bestDist<=(d-imd->mDelta)) {
           mInterMolDistCost += pow((bestDist-(d-imd->mDelta))/imd->mSig, 2);
       } else {
           mInterMolDistCost += pow((bestDist-(d+imd->mDelta))/imd->mSig, 2);
       }              
   }
   
   mInterMolDistCost *= this->GetSpaceGroup().GetNbSymmetrics();         
   mInterMolDistCostClock.Click();   
   std::cout<<"mInterMolDistCost="<<mInterMolDistCost<<"\n";
   std::cout<<"mInterMolDistCostScale="<<mInterMolDistCostScale<<"\n";

   return mInterMolDistCost*mInterMolDistCostScale;
}


const RefinableObjClock& Crystal::GetClockScattererList()const {return mClockScattererList;}

void Crystal::GlobalOptRandomMove(const REAL mutationAmplitude,
                                  const RefParType *type)
{
   if(mRandomMoveIsDone) return;
   VFN_DEBUG_ENTRY("Crystal::GlobalOptRandomMove()",2)
   //Either a random move or a permutation of two scatterers
   const unsigned long nb=(unsigned long)this->GetNbScatterer();
   if( ((rand()/(REAL)RAND_MAX)<.02) && (nb>1))
   {
      // This is safe even if one scatterer is partially fixed,
      // since we the SetX/SetY/SetZ actually use the MutateTo() function.
      const unsigned long n1=rand()%nb;
      const unsigned long n2=(  (rand()%(nb-1)) +n1+1) %nb;
      const float x1=this->GetScatt(n1).GetX();
      const float y1=this->GetScatt(n1).GetY();
      const float z1=this->GetScatt(n1).GetZ();
      this->GetScatt(n1).SetX(this->GetScatt(n2).GetX());
      this->GetScatt(n1).SetY(this->GetScatt(n2).GetY());
      this->GetScatt(n1).SetZ(this->GetScatt(n2).GetZ());
      this->GetScatt(n2).SetX(x1);
      this->GetScatt(n2).SetY(y1);
      this->GetScatt(n2).SetZ(z1);
   }
   else
   {
      this->RefinableObj::GlobalOptRandomMove(mutationAmplitude,type);
   }
   mRandomMoveIsDone=true;
   VFN_DEBUG_EXIT("Crystal::GlobalOptRandomMove()",2)
}

REAL Crystal::GetLogLikelihood()const
{
   return this->GetBumpMergeCost()+this->GetBondValenceCost()+this->GetInterMolDistCost();
}

void Crystal::CIFOutput(ostream &os, double mindist)const
{
   VFN_DEBUG_ENTRY("Crystal::OutputCIF()",5)

   //Data block name (must have no spaces)
   string tempname = mName;
   int where, size;
   size = tempname.size();
   if (size > 32)
   {
      tempname = tempname.substr(0,32);
      size = tempname.size();
   }
   if (size == 0)
   {
      tempname = "3D";
      size = 2;
   }
   where = tempname.find(" ",0);
   while (where != (int)string::npos)
   {
     tempname.replace(where,1,"_");
     where = tempname.find(" ",0);
     //cout << tempname << endl;
   }
   os << "data_" << tempname <<endl<<endl;

   //Program
   os <<"_computing_structure_solution     'FOX http://objcryst.sourceforge.net'"<<endl<<endl;

   //Scattering powers
   /* This is making troubles when the cif file is imported to the JANA2006
   os << "loop_"<<endl
      << "    _atom_type_symbol" <<endl
      << "    _atom_type_description"<<endl
      << "    _atom_type_scat_source" <<endl;
   for(int i=0;i<this->GetScatteringPowerRegistry().GetNb();i++)
      os << "    "
         << this->GetScatteringPowerRegistry().GetObj(i).GetName()<<" "
         <<this->GetScatteringPowerRegistry().GetObj(i).GetSymbol()<<" "
         <<"'International Tables for Crystallography (Vol. IV)'"
         <<endl;
   os <<endl;
   */

   //Symmetry
   os <<"_symmetry_space_group_name_H-M    '"
      << this->GetSpaceGroup().GetCCTbxSpg().match_tabulated_settings().hermann_mauguin();
   const char ext = this->GetSpaceGroup().GetExtension();
   if(isalnum(ext))
      os <<":"<<ext;
   os <<"'"<<endl;
   os <<"_symmetry_space_group_name_Hall   '"
      << this->GetSpaceGroup().GetCCTbxSpg().match_tabulated_settings().hall()<<"'"<<endl;
   os <<endl;

   os << "_cell_length_a    " << FormatFloat(this->GetLatticePar(0),8,5) << endl
      << "_cell_length_b    " << FormatFloat(this->GetLatticePar(1),8,5) << endl
      << "_cell_length_c    " << FormatFloat(this->GetLatticePar(2),8,5) << endl
      << "_cell_angle_alpha " << FormatFloat(this->GetLatticePar(3)*RAD2DEG,7,3) << endl
      << "_cell_angle_beta  " << FormatFloat(this->GetLatticePar(4)*RAD2DEG,7,3) << endl
      << "_cell_angle_gamma " << FormatFloat(this->GetLatticePar(5)*RAD2DEG,7,3) << endl
      << "_cell_volume      " << FormatFloat(this->GetVolume(),7,2);
   os <<endl;
   this->GetScatteringComponentList();

   /*
   TODO:
   loop_
 _symmetry_equiv_pos_site_id
 _symmetry_equiv_pos_as_xyz

   */

   os << "loop_" << endl
      << "    _atom_site_label" <<endl
      << "    _atom_site_type_symbol" <<endl
      << "    _atom_site_fract_x"<<endl
      << "    _atom_site_fract_y" <<endl
      << "    _atom_site_fract_z" <<endl
      << "    _atom_site_U_iso_or_equiv" <<endl
      << "    _atom_site_occupancy" <<endl
      << "    _atom_site_adp_type" <<endl;

   const double BtoU = 1.0 / (8 * M_PI * M_PI);
   std::vector<const ScatteringPower*> anisovec;
   std::vector<std::string> namevec;
   CrystMatrix_REAL minDistTable;
   minDistTable=this->GetMinDistanceTable(-1.);
   unsigned long k=0;
   for(int i=0;i<mScattererRegistry.GetNb();i++)
   {
      //mpScatterrer[i]->Print();
      const ScatteringComponentList list=this->GetScatt(i).GetScatteringComponentList();
      for(int j=0;j<list.GetNbComponent();j++)
      {
         if(0==list(j).mpScattPow) continue;
         bool redundant=false;
         for(unsigned long l=0;l<k;++l) if(abs(minDistTable(l,k))<mindist) redundant=true;//-1 means dist > 10A
         if(!redundant)
         {
            // We can't have spaces in atom labels
            string s=this->GetScatt(i).GetComponentName(j);
            size_t posc=s.find(' ');
            while(posc!=string::npos){s[posc]='~';posc=s.find(' ');}

            bool isiso = list(j).mpScattPow->IsIsotropic();
            if(!isiso)
            {
               anisovec.push_back(list(j).mpScattPow);
               namevec.push_back(s);
            }

            os   << "    "
                 << FormatString(s,10) << " "
                 << FormatString(list(j).mpScattPow->GetSymbol(),8) << " "
                 << FormatFloat(list(j).mX,9,6) << " "
                 << FormatFloat(list(j).mY,9,6) << " "
                 << FormatFloat(list(j).mZ,9,6) << " "
                 << FormatFloat(list(j).mpScattPow->GetBiso()*BtoU,9,6) << " "
                 << FormatFloat(list(j).mOccupancy,6,4)
                 << (isiso ? " Uiso" : " Uani")
                 << endl;
         }
         k++;
      }
   }


   // Handle anisotropic atoms
   if( anisovec.size() > 0 )
   {

      os << endl
         << "loop_" << endl
         << "    _atom_site_aniso_label" << endl
         << "    _atom_site_aniso_U_11" << endl
         << "    _atom_site_aniso_U_22" << endl
         << "    _atom_site_aniso_U_33" << endl
         << "    _atom_site_aniso_U_12" << endl
         << "    _atom_site_aniso_U_13" << endl
         << "    _atom_site_aniso_U_23" << endl;


      for(size_t i = 0; i < anisovec.size(); ++i)
      {
         os << "    " << FormatString(namevec[i],8) << " ";
         for(int j=0; j<6; ++j)
            os << FormatFloat(anisovec[i]->GetBij(j)*BtoU,9,6) << " ";
         os << endl;
      }
   }


   bool first=true;
   k=0;
   for(int i=0;i<mScattererRegistry.GetNb();i++)
   {
      const ScatteringComponentList list=this->GetScatt(i).GetScatteringComponentList();
      for(int j=0;j<list.GetNbComponent();j++)
      {
         if(0==list(j).mpScattPow) continue;
         bool redundant=false;
         for(unsigned long l=0;l<k;++l) if(abs(minDistTable(l,k))<mindist) redundant=true;//-1 means dist > 10A
         if(redundant)
         {
            if(first)
            {
               first=false;
               os <<endl
                  << "# The following atoms have been excluded by Fox because they are"<<endl
                  << "# almost fully overlapping with another atom (d<" << mindist << "A)"<< endl;
            }
            os   << "#    "
                 << FormatString(list(j).mpScattPow->GetName(),8) << " "
                 << FormatString(this->GetScatt(i).GetComponentName(j),10) << " "
                 << FormatFloat(list(j).mX,9,6) << " "
                 << FormatFloat(list(j).mY,9,6) << " "
                 << FormatFloat(list(j).mZ,9,6) << " "
                 << FormatFloat(list(j).mpScattPow->GetBiso()*BtoU,9,6) << " "
                 << FormatFloat(list(j).mOccupancy,6,4)
                 << " Uiso"
                 << endl;
         }
         k++;
      }
   }

   os <<endl;
   k=0;
   if(1==mUseDynPopCorr.GetChoice())
   {
      os << "#  Dynamical occupancy corrections found by ObjCryst++:"<<endl
         << "#  values below 1. (100%) indicate a correction,"<<endl
         << "#  which means either that the atom is on a special position,"<<endl
         << "#  or that it is overlapping with another identical atom."<<endl;
      for(int i=0;i<mScattererRegistry.GetNb();i++)
      {
         //mpScatterrer[i]->Print();
         const ScatteringComponentList list=this->GetScatt(i).GetScatteringComponentList();
         for(int j=0;j<list.GetNbComponent();j++)
         {
            os   << "#   "
                 << FormatString(this->GetScatt(i).GetComponentName(j),16)
                 << " : " << FormatFloat(mScattCompList(k).mDynPopCorr,6,4)
                 << endl;
            k++;
         }
      }
      os << "#"<<endl;
   }

   VFN_DEBUG_EXIT("Crystal::OutputCIF()",5)
}

void Crystal::GetGeneGroup(const RefinableObj &obj,
                                CrystVector_uint & groupIndex,
                                unsigned int &first) const
{
   // One group for all lattice parameters
   unsigned int latticeIndex=0;
   VFN_DEBUG_MESSAGE("Crystal::GetGeneGroup()",4)
   for(long i=0;i<obj.GetNbPar();i++)
      for(long j=0;j<this->GetNbPar();j++)
         if(&(obj.GetPar(i)) == &(this->GetPar(j)))
         {
            //if(this->GetPar(j).GetType()->IsDescendantFromOrSameAs(gpRefParTypeUnitCell))
            //{
               if(latticeIndex==0) latticeIndex=first++;
               groupIndex(i)=latticeIndex;
            //}
            //else //no parameters other than unit cell
         }
}
void Crystal::BeginOptimization(const bool allowApproximations,const bool enableRestraints)
{
   if(this->IsBeingRefined())
   {
      // RefinableObj::BeginOptimization() will increase the optimization depth
      this->RefinableObj::BeginOptimization(allowApproximations,enableRestraints);
      return;
   }
   for(int i=0;i<this->GetScattererRegistry().GetNb();i++)
   {
      this->GetScattererRegistry().GetObj(i).
         SetGlobalOptimStep(gpRefParTypeScattTranslX,0.1/this->GetLatticePar(0));
      this->GetScattererRegistry().GetObj(i).
         SetGlobalOptimStep(gpRefParTypeScattTranslY,0.1/this->GetLatticePar(1));
      this->GetScattererRegistry().GetObj(i).
         SetGlobalOptimStep(gpRefParTypeScattTranslZ,0.1/this->GetLatticePar(2));
   }
   this->RefinableObj::BeginOptimization(allowApproximations,enableRestraints);
   // Calculate mDistTableMaxDistance: Default 1 Angstroem, for dynamical occupancy correction
   mDistTableMaxDistance=10.0;
   // Up to 4 Angstroem if bond-valence is used
   if((mvBondValenceRo.size()>0) && (mBondValenceCostScale>1e-5)) mDistTableMaxDistance=4;
   // Up to whatever antibump distance the user requires (hopefully not too large !)
   for(Crystal::VBumpMergePar::iterator pos=mvBumpMergePar.begin();pos!=mvBumpMergePar.end();++pos)
      if(sqrt(pos->second.mDist2)>mDistTableMaxDistance) mDistTableMaxDistance=sqrt(pos->second.mDist2);
   VFN_DEBUG_MESSAGE("Crystal::BeginOptimization():mDistTableMaxDistance="<<mDistTableMaxDistance,10)
}

void Crystal::AddBondValenceRo(const ScatteringPower &pow1,const ScatteringPower &pow2,const REAL ro)
{
   if(&pow1 < &pow2) mvBondValenceRo[make_pair(&pow1,&pow2)]=ro;
   else mvBondValenceRo[make_pair(&pow2,&pow1)]=ro;
   mBondValenceParClock.Click();
   this->UpdateDisplay();
}

void Crystal::RemoveBondValenceRo(const ScatteringPower &pow1,const ScatteringPower &pow2)
{
   map<pair<const ScatteringPower*,const ScatteringPower*>, REAL>::iterator pos;
   if(&pow1 < &pow2) pos=mvBondValenceRo.find(make_pair(&pow1 , &pow2));
   else pos=mvBondValenceRo.find(make_pair(&pow2 , &pow1));
   if(pos!=mvBondValenceRo.end()) mvBondValenceRo.erase(pos);
   mBondValenceParClock.Click();
}

REAL Crystal::GetBondValenceCost() const
{
   VFN_DEBUG_MESSAGE("Crystal::GetBondValenceCost()?",4)
   if(mBondValenceCostScale<1e-5) return 0.0;
   if(mvBondValenceRo.size()==0)
   {
      mBondValenceCost=0.0;
      return mBondValenceCost*mBondValenceCostScale;
   }
   this->CalcBondValenceSum();
   if(  (mBondValenceCostClock>mBondValenceCalcClock)
      &&(mBondValenceCostClock>this->GetMasterClockScatteringPower())) return mBondValenceCost*mBondValenceCostScale;
   VFN_DEBUG_MESSAGE("Crystal::GetBondValenceCost():"<<mvBondValenceCalc.size()<<" valences",4)
   TAU_PROFILE("Crystal::GetBondValenceCost()","REAL ()",TAU_DEFAULT);
   mBondValenceCost=0.0;
   std::map<long, REAL>::const_iterator pos;
   for(pos=mvBondValenceCalc.begin();pos!=mvBondValenceCalc.end();++pos)
   {
      const REAL a=pos->second-mScattCompList(pos->first).mpScattPow->GetFormalCharge();
      mBondValenceCost += a*a;
      VFN_DEBUG_MESSAGE("Crystal::GetBondValenceCost():"
                        <<mScattCompList(pos->first).mpScattPow->GetName()
                        <<"="<<pos->second,4)
   }
   mBondValenceCostClock.Click();
   return mBondValenceCost*mBondValenceCostScale;
}

std::map<pair<const ScatteringPower*,const ScatteringPower*>, REAL>& Crystal::GetBondValenceRoList()
{ return mvBondValenceRo;}

const std::map<pair<const ScatteringPower*,const ScatteringPower*>, REAL>& Crystal::GetBondValenceRoList()const
{ return mvBondValenceRo;}

void Crystal::CalcBondValenceSum()const
{
   if(mvBondValenceRo.size()==0) return;
   this->CalcDistTable(true);
   VFN_DEBUG_MESSAGE("Crystal::CalcBondValenceSum()?",4)
   if(   (mBondValenceCalcClock>mDistTableClock)
       &&(mBondValenceCalcClock>mBondValenceParClock)) return;
   VFN_DEBUG_MESSAGE("Crystal::CalcBondValenceSum()",4)
   TAU_PROFILE("Crystal::CalcBondValenceSum()","void ()",TAU_DEFAULT);
   mvBondValenceCalc.clear();
   for(long i=0;i<mScattCompList.GetNbComponent();i++)
   {
      const ScatteringPower *pow1=mScattCompList(i).mpScattPow;
      int nb=0;
      REAL val=0.0;
      std::vector<Crystal::Neighbour>::const_iterator pos;
      for(pos=mvDistTableSq[i].mvNeighbour.begin();
          pos<mvDistTableSq[i].mvNeighbour.end();pos++)
      {
         const REAL dist=sqrt(pos->mDist2);
         const REAL occup= mScattCompList(pos->mNeighbourIndex).mOccupancy
                          *mScattCompList(pos->mNeighbourIndex).mDynPopCorr;
         const ScatteringPower *pow2=mScattCompList(pos->mNeighbourIndex).mpScattPow;
         map<pair<const ScatteringPower*,const ScatteringPower*>,REAL>::const_iterator pos;
         if(pow1<pow2) pos=mvBondValenceRo.find(make_pair(pow1,pow2));
         else pos=mvBondValenceRo.find(make_pair(pow2,pow1));
         if(pos!=mvBondValenceRo.end())
         {
            const REAL v=exp((pos->second-dist)/0.37);
            val += occup * v;
            nb++;
         }
      }
      if(nb!=0) mvBondValenceCalc[i]=val;
   }
   mBondValenceCalcClock.Click();
}

void Crystal::Init(const REAL a, const REAL b, const REAL c, const REAL alpha,
                   const REAL beta, const REAL gamma,const string &SpaceGroupId,
                   const string& name)
{
   VFN_DEBUG_MESSAGE("Crystal::Init(a,b,c,alpha,beta,gamma,Sg,name)",10)
   this->UnitCell::Init(a,b,c,alpha,beta,gamma,SpaceGroupId,name);
   mClockScattCompList.Reset();
   mClockNeighborTable.Reset();
   mClockDynPopCorr.Reset();

   VFN_DEBUG_MESSAGE("Crystal::Init(a,b,c,alpha,beta,gamma,Sg,name):End",10)
}

bool CompareBondDist(MolBond* b1, MolBond* b2)
{
   return b1->GetLength()<b2->GetLength();
}

bool ComparePairSecond(const std::pair<int,float> &b1, const std::pair<int,float> &b2)
{
   return b1.second < b2.second;
}

void Crystal::ConnectAtoms(const REAL min_relat_dist, const REAL max_relat_dist, const bool warnuser_fail)
{
   VFN_DEBUG_ENTRY("Crystal::ConnectAtoms(...)",10)
   // Make sure there are only atoms
   for(unsigned int i=0;i<mScattererRegistry.GetNb();i++)
   {
      if(mScattererRegistry.GetObj(i).GetClassName()!="Atom")
      {
         if(warnuser_fail) (*fpObjCrystInformUser)("Crystal::ConnectAtoms(): cannot connect atoms unless there are only Atoms in a Crystal");
         return;
      }
   }
   this->CalcDistTable(true);
   this->GetScatteringComponentList();

   // Create first Molecule
   // start from start_carbon
   Molecule *pmol=NULL;
   std::set<int> vAssignedAtoms;//atoms already assigned to a Molecule
   std::set<int> vTriedAtom0;//atoms already tried as starting point for a Molecule
   VFN_DEBUG_MESSAGE("Crystal::ConnectAtoms(...)",10)
   while(long(vAssignedAtoms.size()) != mScattererRegistry.GetNb())
   {
      VFN_DEBUG_MESSAGE("Crystal::ConnectAtoms(...): new Molecule ?",7)
      // We need at least one carbon atom to start
      int atom0=-1;
      int maxAtomicNumber = 0;
      for(unsigned int i=0;i<mScattCompList.GetNbComponent();i++)
      {
         if((vAssignedAtoms.count(i)>0) || (vTriedAtom0.find(i)!=vTriedAtom0.end()) ) continue;
         if(mScattCompList(i).mpScattPow->GetClassName()=="ScatteringPowerAtom")
         {
            const ScatteringPowerAtom *p=dynamic_cast<const ScatteringPowerAtom*>(mScattCompList(i).mpScattPow);
            if((p->GetAtomicNumber()==6)&&(maxAtomicNumber!=6))
            {// Start from the first Carbon found
               atom0=i;
               maxAtomicNumber=6;
            }
            else if((p->GetAtomicNumber()>maxAtomicNumber) && (maxAtomicNumber!=6))
            {// Else we'll try from the heaviest atom
               maxAtomicNumber=p->GetAtomicNumber();
               atom0=i;
            }
         }
         else
         {
            if(warnuser_fail)
               (*fpObjCrystInformUser)("Crystal::ConnectAtoms(): cannot connect atoms unless there are only Atoms in a Crystal");
            VFN_DEBUG_EXIT("Crystal::ConnectAtoms(...):cannot connect atoms unless there are only Atoms in the structure:"<<i<<":"<<mScattCompList(i).mpScattPow->GetClassName(),10)
            return;
         }
      }
      if(atom0<0)
      {
         if((pmol==NULL) && warnuser_fail) // We did not create any Molecule :-(
         {
            (*fpObjCrystInformUser)("Crystal::ConnectAtoms(): cannot connect atoms unless there is at least one carbon atom");
            VFN_DEBUG_EXIT("Crystal::ConnectAtoms(...):cannot connect atoms unless there is at least one carbon atom",10)
            return;
         }
         break;
      }
      vTriedAtom0.insert(atom0);
      // Atoms in Molecule but for which neighbors have not yet been searched
      // first: index in the Crystal's scatt comp list, second: index in the Molecule
      std::map<int,int>newAtoms;
      // List of all atoms in the Molecule. First is the MolAtom* in the molecule, second is the index in the Crystal
      std::map<MolAtom*,int> molAtoms;

      pmol=new Molecule(*this);

      // Add atom0 to Molecule.
      newAtoms[atom0]=0;
      vAssignedAtoms.insert(atom0);
      const ScatteringPowerAtom *p0=dynamic_cast<const ScatteringPowerAtom*>(mScattCompList(atom0).mpScattPow);
      {
         REAL x=mScattCompList(atom0).mX;
         REAL y=mScattCompList(atom0).mY;
         REAL z=mScattCompList(atom0).mZ;
         const REAL occ=mScattCompList(atom0).mOccupancy;
         this->FractionalToOrthonormalCoords(x,y,z);
         pmol->AddAtom(x,y,z,p0,mScattererRegistry.GetObj(atom0).GetName(),false);
         pmol->GetAtomList().back()->SetOccupancy(occ);
      }
      molAtoms[pmol->GetAtomList().back()]=atom0;
      // Count atoms in the Molecule per element type
      vector<unsigned int> vElementCount(140);// Should be safe pending a trans-uranian breakthrough
      for(vector<unsigned int>::iterator pos=vElementCount.begin();pos!=vElementCount.end();++pos) *pos=0;
      vElementCount[p0->GetAtomicNumber()]+=1;
      while(newAtoms.size()>0)
      {
         atom0=newAtoms.begin()->first;
         p0=dynamic_cast<const ScatteringPowerAtom*>(mScattCompList(atom0).mpScattPow);
         VFN_DEBUG_ENTRY("Crystal::ConnectAtoms(...):atom0="<<atom0<<","<<mScattererRegistry.GetObj(atom0).GetName()<<"["<<p0->GetSymbol()<<"],"<<newAtoms.size()<<" new atoms",7)
         std::vector<std::pair<int,float> > vatomsneighbours;
         // Find neigbours between min and max * sum of covalent bonds
         for(std::vector<Crystal::Neighbour>::const_iterator pos=mvDistTableSq[atom0].mvNeighbour.begin();
             pos!=mvDistTableSq[atom0].mvNeighbour.end();pos++)
         {
            const ScatteringPowerAtom *p1=dynamic_cast<const ScatteringPowerAtom*>(mScattCompList(pos->mNeighbourIndex).mpScattPow);
            const REAL dcov= p0->GetCovalentRadius()+p1->GetCovalentRadius();
            //if(  (vAssignedAtoms.count(pos->mNeighbourIndex)!=0) || (p1->GetMaxCovBonds()==0))
            if(p1->GetMaxCovBonds()==0)
               continue;
            if(  ((min_relat_dist*dcov)<sqrt(pos->mDist2))
               &&((max_relat_dist*dcov)>sqrt(pos->mDist2)))
            {
               VFN_DEBUG_MESSAGE("Crystal::ConnectAtoms(...):atom0="<<mScattererRegistry.GetObj(atom0).GetName()<<"-"<<p1->GetName()<<"("<<p1->GetMaxCovBonds()<<"):dcov="<<dcov<<",d="<<sqrt(pos->mDist2),7)
               vatomsneighbours.push_back(std::make_pair(pos->mNeighbourIndex,sqrt(pos->mDist2)));
            }
            else
            {
               VFN_DEBUG_MESSAGE("Crystal::ConnectAtoms(...):atom0="<<mScattererRegistry.GetObj(atom0).GetName()<<"-"<<p1->GetName()<<"("<<p1->GetMaxCovBonds()<<"):dcov="<<dcov<<",d="<<sqrt(pos->mDist2),5)
            }
         }
         // Remove farthest neighbours if in excess of the maximum coordination.
         // But keep excess neighbours if close (5%) of the last neighbour within the normal coordination number.
         const unsigned int maxbonds=p0->GetMaxCovBonds();
         float extra=vatomsneighbours.size()-maxbonds;
         if(extra>0)
         {
            // Check real number of bonds taking into account occupancy, and sort bonds by length
            REAL nbbond=0;
            for(std::vector<std::pair<int,float> >::const_iterator pos=vatomsneighbours.begin();pos!=vatomsneighbours.end();++pos)
            {
               nbbond+=mScattCompList(pos->first).mOccupancy;
            }
            extra = nbbond-maxbonds;
            if(extra>0.2)
            {
               std::sort(vatomsneighbours.begin(), vatomsneighbours.end(),ComparePairSecond);
               VFN_DEBUG_ENTRY("Crystal::ConnectAtoms(...): too many bonds for"<<mScattererRegistry.GetObj(atom0).GetName()
                               <<" ?(allowed="<<maxbonds<<",nb="<<vatomsneighbours.size()<<",nb_occ="<<nbbond<<",longest="<<vatomsneighbours.back().second<<"["
                               <<mScattererRegistry.GetObj(atom0).GetName()<<"-"<<mScattererRegistry.GetObj(vatomsneighbours.back().first).GetName()<<"])", 10)
               const REAL maxdist=vatomsneighbours[vatomsneighbours.size()-extra].second*1.05;
               while(vatomsneighbours.back().second>maxdist)
               {
                  VFN_DEBUG_MESSAGE("Crystal::ConnectAtoms(...): Remove bond="<<mScattererRegistry.GetObj(atom0).GetName()
                                    <<"-"<<mScattererRegistry.GetObj(vatomsneighbours.back().first).GetName()<<", d="<<vatomsneighbours.back().second,10)
                  vatomsneighbours.pop_back();
               }
               VFN_DEBUG_EXIT("Crystal::ConnectAtoms(...): too many bonds for"<<mScattererRegistry.GetObj(atom0).GetName()
                              <<" ?(allowed="<<maxbonds<<",nb="<<vatomsneighbours.size()<<",longest="<<vatomsneighbours.back().second<<"["
                              <<mScattererRegistry.GetObj(atom0).GetName()<<"-"<<mScattererRegistry.GetObj(vatomsneighbours.back().first).GetName()<<"])", 10)
            }
         }

         // Add remaining atoms to molecule and mark them as asigned
         //if(extra>0) cout<<"Crystal::ConnectAtoms(): adding neighbours around "<<mScattererRegistry.GetObj(atom0).GetName()<<" :";
         for(std::vector<std::pair<int,float> >::const_iterator pos=vatomsneighbours.begin();pos!=vatomsneighbours.end();++pos)
         {
            if(vAssignedAtoms.count(pos->first)!=0)
            {
               if(extra>0) cout<<"("<<mScattererRegistry.GetObj(pos->first).GetName()<<") ";
               continue;
            }
            //if(extra>0) cout<<mScattererRegistry.GetObj(pos->first).GetName()<<" ";
            const ScatteringPowerAtom *p1=dynamic_cast<const ScatteringPowerAtom*>(mScattCompList(pos->first).mpScattPow);
            vAssignedAtoms.insert(pos->first);
            REAL x=mScattCompList(pos->first).mX;
            REAL y=mScattCompList(pos->first).mY;
            REAL z=mScattCompList(pos->first).mZ;
            this->FractionalToOrthonormalCoords(x,y,z);
            const REAL occ=mScattCompList(pos->first).mOccupancy;
            pmol->AddAtom(x,y,z,p1,mScattererRegistry.GetObj(pos->first).GetName(),false);
            pmol->GetAtomList().back()->SetOccupancy(occ);
            newAtoms[pos->first]=pmol->GetNbComponent()-1;
            molAtoms[pmol->GetAtomList().back()]=pos->first;
            vElementCount[p1->GetAtomicNumber()]+=1;
         }
         //if(extra>0) cout<<endl;
         newAtoms.erase(atom0);
         VFN_DEBUG_EXIT("Crystal::ConnectAtoms(...):atom0="<<atom0<<","<<mScattererRegistry.GetObj(atom0).GetName()<<"["<<p0->GetSymbol()<<"],"<<newAtoms.size()<<" new atoms. Temp Molecule:"<<pmol->GetFormula(),7)
      }
      // Check this is a valid Molecule object
      bool keep=false;
      if((vElementCount[6]>0) && (pmol->GetNbComponent()>=3)) keep=true;
      else
      {// no carbon ?

         std::vector<unsigned int> vnb;
         #ifdef __DEBUG__
         cout<<"  Crystal::ConnectAtoms(..): Molecule ?";
         #endif
         for(unsigned int i=0;i<vElementCount.size();++i)
            if(vElementCount[i]!=0)
            {
               vnb.push_back(vElementCount[i]);
               #ifdef __DEBUG__
               cout<<"Z="<<i<<"("<< vElementCount[i]<<") ";
               #endif
            }
         #ifdef __DEBUG__
         cout<<endl;
         #endif
         if(vnb.size()==2)
         {
            #if 0
            if((vElementCount[8]==1) && (vElementCount[1]==2)) keep=true; //H2O
            if((vElementCount[8]==1) && (vElementCount[1]==3)) keep=true; //H3O+
            if((vElementCount[7]==1) && (vElementCount[1]==3)) keep=true; //NH3
            if((vElementCount[7]==1) && (vElementCount[1]==4)) keep=true; //NH4+
            if((vElementCount[7]==1) && (vElementCount[8]==2)) keep=true; //NO2
            if((vElementCount[7]==1) && (vElementCount[8]==3)) keep=true; //NO3-
            if((vElementCount[5]==1) && (vElementCount[1]==3)) keep=true; //BH3
            if((vElementCount[5]==1) && (vElementCount[1]==4)) keep=true; //BH4-
            if((vElementCount[14]==1) && (vElementCount[8]==4)) keep=true; //SiO4
            if((vElementCount[15]==1) && (vElementCount[8]==4)) keep=true; //PO4
            #endif
            #if 0
            // Accept any type of small molecule/polyedra with one center atom
            if( ((vnb[0]==1)||(vnb[1]==1)) && ((vnb[0]+vnb[1])>2)) keep=true;
            #endif
            // Accept any type of cluster with exactly two types of atoms
            keep=true;
         }
      }
      if(!keep)
      {
         VFN_DEBUG_MESSAGE("Crystal::ConnectAtoms(...):Rejected molecule: "<<pmol->GetFormula(),10)
         delete pmol;
         for(std::map<MolAtom*,int>::const_iterator pos=molAtoms.begin();pos!=molAtoms.end();++pos) vAssignedAtoms.erase(pos->second);
         continue;// Will start from another atom to build a molecule
      }

      // Add bonds
      for(unsigned int i=0;i<pmol->GetAtomList().size();i++)
      {
         for(unsigned int j=i+1;j<pmol->GetAtomList().size();j++)
         {
            const REAL d=GetBondLength(pmol->GetAtom(i), pmol->GetAtom(j));
            const REAL dcov= dynamic_cast<const ScatteringPowerAtom*>(&(pmol->GetAtom(i).GetScatteringPower()))->GetCovalentRadius()
                            +dynamic_cast<const ScatteringPowerAtom*>(&(pmol->GetAtom(j).GetScatteringPower()))->GetCovalentRadius();
            if(  ((min_relat_dist*dcov)<d) && ((max_relat_dist*dcov)>d))
            {
               VFN_DEBUG_MESSAGE("Crystal::ConnectAtoms(...):Add bond="<<pmol->GetAtom(i).GetName()<<"-"<<pmol->GetAtom(j).GetName()<<", d="<<d<<"(dcov="<<dcov<<")",6)
               pmol->AddBond(pmol->GetAtom(i),pmol->GetAtom(j),d,.01,.02,false);
            }
         }
      }
      // Remove longest bonds if it exceeds the expected coordination
      // :TODO: combined with the check already made, this is not fullproof, for atoms where the coordination number is not so well-defined,
      // e.g. Li and Na is defined as 1, but there could be more linked atoms...
      // If we still find a too great coordination number, remove excess ones but still keep those that are very close (5%) of the cutoff distance.
      for(vector<MolAtom*>::iterator pos=pmol->GetAtomList().begin();pos!=pmol->GetAtomList().end();)
      {
         pmol->BuildConnectivityTable();
         map<MolAtom *,set<MolAtom *> >::const_iterator p=pmol->GetConnectivityTable().find(*pos);
         if(p==pmol->GetConnectivityTable().end())
         {// While cleaning the longest bond, this atom had all his bonds removed !
            VFN_DEBUG_MESSAGE("Crystal::ConnectAtoms(...):no bond remaining for:"<<(*pos)->GetName()<<"! Removing atom from Molecule",10)
            //Remove MolAtom from Molecule and keep in Crystal.
            vAssignedAtoms.erase(molAtoms[*pos]);
            molAtoms.erase(*pos);
            pos=pmol->RemoveAtom(**pos);
            continue;
         }
         const unsigned int maxbonds=dynamic_cast<const ScatteringPowerAtom*>(&(p->first->GetScatteringPower()))->GetMaxCovBonds();
         int extra=p->second.size()-maxbonds;
         if(extra>0)
         {
            // Check real number of bonds taking into account occupancy, and sort bonds by length
            std::vector<MolBond*> vbonds;
            REAL nbbond=0;
            for(std::set<MolAtom*>::iterator p1=p->second.begin();p1!=p->second.end();++p1)
            {
               vbonds.push_back(*(pmol->FindBond(**pos,**p1)));// We can assume that exactly one bond is found
               nbbond+=(*p1)->GetOccupancy();
            }
            VFN_DEBUG_MESSAGE("Crystal::ConnectAtoms(...): too many bonds for"<<(*pos)->GetName()<<" ?(allowed="<<maxbonds<<",nb="<<p->second.size()<<",nb_occ="<<nbbond<<")", 10)
            const int extra= (int)(nbbond-maxbonds);
            if(extra>0)
            {
               std::sort(vbonds.begin(), vbonds.end(),CompareBondDist);
               if(size_t(extra) < vbonds.size()) // Am I paranoid ?
               {
                  const REAL maxdist=vbonds[vbonds.size()-extra]->GetLength()*1.05;
                  while(vbonds.back()->GetLength()>maxdist)
                  {
                     VFN_DEBUG_MESSAGE("Crystal::ConnectAtoms(...): Remove bond="<<vbonds.back()->GetAtom1().GetName()<<"-"<<vbonds.back()->GetAtom2().GetName()<<", d="<<vbonds.back()->GetLength(),10)
                     pmol->RemoveBond(*(vbonds.back()));
                     vbonds.pop_back();
                     VFN_DEBUG_MESSAGE("Crystal::ConnectAtoms(...): Next bond  ="<<vbonds.back()->GetAtom1().GetName()<<"-"<<vbonds.back()->GetAtom2().GetName()<<", d="<<vbonds.back()->GetLength(),10)
                  }
               }
            }
         }
         ++pos;
      }
      // Add all bond angles
      pmol->BuildConnectivityTable();
      for(map<MolAtom*,set<MolAtom*> >::const_iterator pos=pmol->GetConnectivityTable().begin();
          pos!=pmol->GetConnectivityTable().end();++pos)
      {
         for(set<MolAtom*>::const_iterator pos1=pos->second.begin();
             pos1!=pos->second.end();++pos1)
         {
            for(set<MolAtom*>::const_iterator pos2=pos1;
                pos2!=pos->second.end();++pos2)
            {
               if(pos2==pos1) continue;
               if(pmol->FindBondAngle(**pos1,*(pos->first),**pos2)== pmol->GetBondAngleList().end())
               pmol->AddBondAngle(**pos1,*(pos->first),**pos2,
                                 GetBondAngle(**pos1,*(pos->first),**pos2),0.01,0.02,false);
            }
         }
      }
      // Correct center of Molecule
      REAL xc=0,yc=0,zc=0;
      for(std::map<MolAtom*,int>::const_iterator pos=molAtoms.begin();pos!=molAtoms.end();++pos)
      {
         REAL x=mScattCompList(pos->second).mX;
         REAL y=mScattCompList(pos->second).mY;
         REAL z=mScattCompList(pos->second).mZ;
         this->FractionalToOrthonormalCoords(x,y,z);
         xc+=x;
         yc+=y;
         zc+=z;
      }
      xc /= pmol->GetNbComponent();
      yc /= pmol->GetNbComponent();
      zc /= pmol->GetNbComponent();
      this->OrthonormalToFractionalCoords(xc,yc,zc);
      VFN_DEBUG_MESSAGE("Crystal::ConnectAtoms(...): center?"<<pmol->GetNbComponent()<<","<<molAtoms.size()<<":"<<xc<<","<<yc<<","<<zc,10)
      pmol->SetX(xc);
      pmol->SetY(yc);
      pmol->SetZ(zc);
      this->AddScatterer(pmol);
      (*fpObjCrystInformUser)("ConnectAtoms: found Molecule: "+pmol->GetFormula());
   }
   std::set<Scatterer*> vpAtom;
   for(std::set<int>::const_iterator pos=vAssignedAtoms.begin();pos!=vAssignedAtoms.end();++pos)
   {
      Scatterer *s=&(this->GetScattererRegistry().GetObj(*pos));
      vpAtom.insert(s);
   }
   while(vpAtom.size()>0)
   {
      VFN_DEBUG_MESSAGE("Crystal::ConnectAtoms(...): remove atom:"<<(*vpAtom.begin())->GetName()<<","<<vpAtom.size()<<" remaining...",6)
      this->RemoveScatterer(*vpAtom.begin(),true);
      vpAtom.erase(vpAtom.begin());
   }
   VFN_DEBUG_EXIT("Crystal::ConnectAtoms(...)",10)
}

void Crystal::MergeEqualScatteringPowers(const bool oneScatteringPowerPerElement)
{
   VFN_DEBUG_ENTRY("Crystal::MergeEqualScatteringPowers("<<oneScatteringPowerPerElement<<")", 10)
   // Find identical scattering powers.
   std::set<ScatteringPower*> vremovedpow;
   std::map<ScatteringPower*,std::set<ScatteringPower*> > vequivpow;
   for(unsigned int i=0;i<this->GetScatteringPowerRegistry().GetNb();i++)
   {
      ScatteringPower *p1 = &(this->GetScatteringPowerRegistry().GetObj(i));
      if(vremovedpow.find(p1)!=vremovedpow.end()) continue;
      vequivpow[p1] = std::set<ScatteringPower*>();
      for(unsigned int j=i+1;j<this->GetScatteringPowerRegistry().GetNb();j++)
      {
         ScatteringPower *p2 = &(this->GetScatteringPowerRegistry().GetObj(j));
         if(oneScatteringPowerPerElement)
         {
            if(p1->GetClassName() != p2->GetClassName()) continue;
            if(p1->GetSymbol() != p2->GetSymbol()) continue;
         }
         else
         {
            if(*p1 != *p2) continue;
         }
         vequivpow[p1].insert(p2);
         vremovedpow.insert(p2);
      }
   }
   if(oneScatteringPowerPerElement)
   {
      // Average Biso and Bij
      for(std::map<ScatteringPower*,std::set<ScatteringPower*> >::iterator pos=vequivpow.begin();pos!=vequivpow.end();++pos)
      {
         if(pos->second.size()==0) continue;
         REAL b = pos->first->GetBiso();
         CrystVector_REAL bij(6);
         for(unsigned int i=0;i<6;i++) bij(i) = pos->first->GetBij(i);
         for(std::set<ScatteringPower*>::const_iterator pos2=pos->second.begin(); pos2!=pos->second.end();++pos2)
         {
            b += (*pos2)->GetBiso();
            for(unsigned int i=0;i<6;i++) bij(i) += (*pos2)->GetBij(i);
         }
         b   /= pos->second.size() + 1;
         bij /= pos->second.size() + 1;
         pos->first->SetBiso(b);
         for(unsigned int i=0;i<6;i++) if(abs(bij(i)) > 1e-6) pos->first->SetBij(i,bij(i));
      }
   }
   // Update Atoms or MolAtoms with new ScatteringPower
   for(std::map<ScatteringPower*,std::set<ScatteringPower*> >::iterator pos=vequivpow.begin();pos!=vequivpow.end();++pos)
   {
      const unsigned int nb = pos->second.size();
      if(oneScatteringPowerPerElement) pos->first->SetName(pos->first->GetSymbol());
      if(nb>0)
         (*fpObjCrystInformUser)((boost::format("Merging ScatteringPower: %s[%s] (%d identical scattering powers)") % pos->first->GetName().c_str() % pos->first->GetSymbol().c_str() % pos->second.size()).str());
      for(std::set<ScatteringPower*>::const_iterator pos2=pos->second.begin(); pos2!=pos->second.end();++pos2)
      {
         for(unsigned int i=0;i<this->GetNbScatterer();++i)
         {
            Scatterer *p = &(this->GetScatt(i));
            if(p->GetClassName()=="Atom")
            {
               Atom *pat=dynamic_cast<Atom*>(p);
               if(&(pat->GetScatteringPower()) == (*pos2))
               {
                  VFN_DEBUG_MESSAGE("Crystal:MergeEqualScatteringPowers() Atom "<<pat->GetName()<<": "<<pat->GetScatteringPower().GetName()<<"->"<<pos->first->GetName(), 10)
                  pat->SetScatteringPower(*(pos->first));
               }
            }
            else if (p->GetClassName()=="Molecule")
            {
               Molecule *pmol=dynamic_cast<Molecule*>(p);
               for(std::vector<MolAtom*>::iterator pat=pmol->GetAtomList().begin();pat!=pmol->GetAtomList().end();++pat)
               {
                  if(&((*pat)->GetScatteringPower()) ==  (*pos2))
                     (*pat)->SetScatteringPower(*(pos->first));
               }
            }
            else
            {
               // This should only happen if a new type of scatterer was derived
               cout<<__FILE__<<":"<<__LINE__<<":Crystal::MergeEqualScatteringPowers(): unidentified scatterer, cannot merge scattering power..."
                   <<(*pos2)->GetName()<<"["<<(*pos2)->GetClassName()<<"]"<<endl;
            }
         }
      }
   }
   // Delete duplicate scattering powers
   for(std::set<ScatteringPower*>::iterator pos=vremovedpow.begin();pos!=vremovedpow.end();++pos)
   {
      #ifdef __DEBUG__
      const unsigned int nb=(*pos)->GetClientRegistry().GetNb();
      if(nb>0)
      {
         VFN_DEBUG_MESSAGE("Crystal::MergeEqualScatteringPowers(): "<<nb<<" clients remaining for scattering power: "<<(*pos)->GetName()<<"["<<(*pos)->GetClassName()<<"]", 5)
         for(unsigned int i=0; i<nb;i++)
         {
            VFN_DEBUG_MESSAGE("                                       "<<&((*pos)->GetClientRegistry().GetObj(i))<<":"<<(*pos)->GetClientRegistry().GetObj(i).GetName()<<"["<<(*pos)->GetClientRegistry().GetObj(i).GetClassName()<<"]", 5)
         }
      }
      #endif
      this->RemoveScatteringPower(*pos,true);
   }
   this->UpdateDisplay();
   VFN_DEBUG_EXIT("Crystal::MergeEqualScatteringPowers()", 10)
}

void Crystal::InitOptions()
{
   VFN_DEBUG_ENTRY("Crystal::InitOptions",10)
   static string UseDynPopCorrname;
   static string UseDynPopCorrchoices[2];

   static string DisplayEnantiomername;
   static string DisplayEnantiomerchoices[2];

   static bool needInitNames=true;
   if(true==needInitNames)
   {
      UseDynPopCorrname="Use Dynamical Occupancy Correction";
      UseDynPopCorrchoices[0]="No";
      UseDynPopCorrchoices[1]="Yes";

      DisplayEnantiomername="Display Enantiomer";
      DisplayEnantiomerchoices[0]="No";
      DisplayEnantiomerchoices[1]="Yes";

      needInitNames=false;//Only once for the class
   }
   VFN_DEBUG_MESSAGE("Crystal::Init(a,b,c,alpha,beta,gamma,Sg,name):Init options",5)
   mUseDynPopCorr.Init(2,&UseDynPopCorrname,UseDynPopCorrchoices);
   mUseDynPopCorr.SetChoice(1);
   this->AddOption(&mUseDynPopCorr);

   mDisplayEnantiomer.Init(2,&DisplayEnantiomername,DisplayEnantiomerchoices);
   mDisplayEnantiomer.SetChoice(0);
   this->AddOption(&mDisplayEnantiomer);
   VFN_DEBUG_EXIT("Crystal::InitOptions",10)
}

Crystal::Neighbour::Neighbour(const unsigned long neighbourIndex,const int sym,
                              const REAL dist2):
mNeighbourIndex(neighbourIndex),mNeighbourSymmetryIndex(sym),mDist2(dist2)
{}

struct DistTableInternalPosition
{
   DistTableInternalPosition(const long atomIndex, const int sym,
                             const REAL x,const REAL y,const REAL z):
   mAtomIndex(atomIndex),mSymmetryIndex(sym),mX(x),mY(y),mZ(z)
   {}
   /// Index of the atom (order) in the component list
   long mAtomIndex;
   /// Which symmetry operation does this symmetric correspond to ?
   int mSymmetryIndex;
   /// Fractionnal coordinates
   REAL mX,mY,mZ;
};
/*
void Crystal::CalcMyDistTable() const
{
   this->GetScatteringComponentList();

   cout<<"finding the maximal distance in the user-defined list\n";
   //finding the maximal distance in the user-defined list
   float maxDist = -1;
   for(int i=0;i<mInterMolDistList.size();i++) {
       if(mInterMolDistList[i].mDist2>0) {
           float d = sqrt(mInterMolDistList[i].mDist2);
           if(d>maxDist) {
               maxDist = d;
           }
       }
   }
   cout<<"maxDist = "<<maxDist<<"\n";
   if(maxDist<0) return;
   maxDist*=2;
   cout<<"maxDist = "<<maxDist<<"\n";

   if(true) {//calculate the table (imdTable) just for user-defined atoms 
    //this code is almost the same as in this->CalcDistTable(true);

       //this->GetScatteringComponentList();

       const long nbComponent=mScattCompList.GetNbComponent();

       // Get range and origin of the (pseudo) asymmetric unit
        const REAL asux0=this->GetSpaceGroup().GetAsymUnit().Xmin();
        const REAL asuy0=this->GetSpaceGroup().GetAsymUnit().Ymin();
        const REAL asuz0=this->GetSpaceGroup().GetAsymUnit().Zmin();

        const REAL asux1=this->GetSpaceGroup().GetAsymUnit().Xmax();
        const REAL asuy1=this->GetSpaceGroup().GetAsymUnit().Ymax();
        const REAL asuz1=this->GetSpaceGroup().GetAsymUnit().Zmax();

        const REAL halfasuxrange=(asux1-asux0)*0.5+1e-5;
        const REAL halfasuyrange=(asuy1-asuy0)*0.5+1e-5;
        const REAL halfasuzrange=(asuz1-asuz0)*0.5+1e-5;

        const REAL asuxc=0.5*(asux0+asux1);
        const REAL asuyc=0.5*(asuy0+asuy1);
        const REAL asuzc=0.5*(asuz0+asuz1);

        const REAL maxdx=halfasuxrange+maxDist/GetLatticePar(0);
        const REAL maxdy=halfasuyrange+maxDist/GetLatticePar(1);
        const REAL maxdz=halfasuzrange+maxDist/GetLatticePar(2);

       // List of all positions within or near the first atom generated
       std::vector<DistTableInternalPosition> vPos;
       // index of unique atoms in vPos, which are strictly in the asymmetric unit
       std::vector<unsigned long> vUniqueIndex(nbComponent);

       const REAL asymUnitMargin2 = maxDist*maxDist;

       {
          // No need to loop on a,b,c translations if maxDist is small enough
          bool loopOnLattice=true;
          if(  ((this->GetLatticePar(0)*.5)>maxDist)
             &&((this->GetLatticePar(1)*.5)>maxDist)
             &&((this->GetLatticePar(2)*.5)>maxDist)) loopOnLattice=false;

          
          //allocate the table 
          imdTable.resize(nbComponent);
          {
              std::vector<NeighbourHood>::iterator pos;
              for(pos=imdTable.begin();pos<imdTable.end();pos++)
                 pos->mvNeighbour.clear();
          }

          CrystMatrix_REAL symmetricsCoords;
          const int nbSymmetrics=this->GetSpaceGroup().GetNbSymmetrics(false,false);

          // Get the list of all atoms within or near the asymmetric unit
          for(long i=0;i<nbComponent;i++)
          {
             VFN_DEBUG_MESSAGE("Crystal::CalcDistTable(fast):3:component "<<i,0)
             // generate all symmetrics, excluding translations
             symmetricsCoords=this->GetSpaceGroup().GetAllSymmetrics(mScattCompList(i).mX,
                                                                     mScattCompList(i).mY,
                                                                     mScattCompList(i).mZ,
                                                                     false,false,false);
             imdTable[i].mIndex=i;//USELESS ?
             bool hasUnique=false;
             for(int j=0;j<nbSymmetrics;j++)
             {
                // take the closest position (using lattice translations) to the center of the ASU
                REAL x=fmod(symmetricsCoords(j,0)-asuxc,(REAL)1.0);if(x<-.5)x+=1;else if(x>.5)x-=1;
                REAL y=fmod(symmetricsCoords(j,1)-asuyc,(REAL)1.0);if(y<-.5)y+=1;else if(y>.5)y-=1;
                REAL z=fmod(symmetricsCoords(j,2)-asuzc,(REAL)1.0);if(z<-.5)z+=1;else if(z>.5)z-=1;

                //cout<<i<<","<<j<<":"<<FormatFloat(x,8,5)<<","<<FormatFloat(y,8,5)<<","<<FormatFloat(z,8,5)<<endl;
                if( (abs(x)<maxdx) && (abs(y)<maxdy) && (abs(z)<maxdz) )
                   vPos.push_back(DistTableInternalPosition(i, j, x+asuxc, y+asuyc, z+asuzc));
                // Get one reference atom strictly within the pseudo-ASU
                if(!hasUnique)
                   if( (abs(x)<halfasuxrange) && (abs(y)<halfasuyrange) && (abs(z)<halfasuzrange) )
                   {
                      hasUnique=true;
                      vUniqueIndex[i]=vPos.size()-1;
                      imdTable[i].mUniquePosSymmetryIndex=j;
                   }
             }
             if(!hasUnique)
             {
                throw ObjCrystException("One atom did not have any symmetric in the ASU !");
             }
          }

          const CrystMatrix_REAL* pOrthMatrix=&(this->GetOrthMatrix());

          const REAL m00=(*pOrthMatrix)(0,0);
          const REAL m01=(*pOrthMatrix)(0,1);
          const REAL m02=(*pOrthMatrix)(0,2);
          const REAL m11=(*pOrthMatrix)(1,1);
          const REAL m12=(*pOrthMatrix)(1,2);
          const REAL m22=(*pOrthMatrix)(2,2);

          cout<<"Get the list of all atoms in the asymmetric unit in mInterMolDistList\n";
          // Get the list of all atoms in the asymmetric unit in mInterMolDistList
          vector<int> listOfCurrentScatterers;
          for(long i=0;i<mInterMolDistList.size();i++) {
              vector<int> p = FindScatterersInComponentList(mInterMolDistList[i].mAt1);
              listOfCurrentScatterers.insert(listOfCurrentScatterers.end(), p.begin(), p.end());
          }
          if(listOfCurrentScatterers.size()==0) return;
          cout<<"listOfCurrentScatterers="<<listOfCurrentScatterers.size()<<"\n";

          for(long i=0;i<nbComponent;i++)
          {
             bool found = false;
             for(int q=0;q<listOfCurrentScatterers.size();q++) {
                 if(i==listOfCurrentScatterers[q]) {
                     found = true;
                     break;
                 }
             }
             //calculate the table just for user-define atoms 
             if(!found) continue;

             VFN_DEBUG_MESSAGE("Crystal::CalcDistTable(fast):4:component "<<i,0)
             #if 0
             if(!this->IsBeingRefined()) cout<<endl<<"Unique pos:"<<vUniqueIndex[i]<<":"
                 <<vPos[vUniqueIndex[i]].mAtomIndex<<":"
                 <<mScattCompList(vPos[vUniqueIndex[i]].mAtomIndex).mpScattPow->GetName()<<":"
                 <<vPos[vUniqueIndex[i]].mSymmetryIndex<<":"
                 <<FormatFloat(vPos[vUniqueIndex[i]].mX,8,5)<<","
                 <<FormatFloat(vPos[vUniqueIndex[i]].mY,8,5)<<","
                 <<FormatFloat(vPos[vUniqueIndex[i]].mZ,8,5)<<endl;
             #endif
             std::vector<Crystal::Neighbour> * const vnb=&(imdTable[i].mvNeighbour);
             const REAL x0i=vPos[vUniqueIndex[i] ].mX;
             const REAL y0i=vPos[vUniqueIndex[i] ].mY;
             const REAL z0i=vPos[vUniqueIndex[i] ].mZ;
             for(unsigned long j=0;j<vPos.size();j++)
             {
                if((vUniqueIndex[i]==j) && (!loopOnLattice)) continue;// distance to self !
                // Start with the smallest absolute coordinates possible
                REAL x=fmod(vPos[j].mX - x0i,(REAL)1.0);if(x<-.5)x+=1;if(x>.5)x-=1;
                REAL y=fmod(vPos[j].mY - y0i,(REAL)1.0);if(y<-.5)y+=1;if(y>.5)y-=1;
                REAL z=fmod(vPos[j].mZ - z0i,(REAL)1.0);if(z<-.5)z+=1;if(z>.5)z-=1;

                const REAL x0=m00 * x + m01 * y + m02 * z;
                const REAL y0=          m11 * y + m12 * z;
                const REAL z0=                    m22 * z;

                if(loopOnLattice)// distance to self !
                {//Now loop over lattice translations
                   for(int sz=-1;sz<=1;sz+=2)// Sign of translation
                   {
                      for(int nz=(sz+1)/2;;++nz)
                      {
                         const REAL z=z0+sz*nz*m22;
                         if(abs(z)>maxDist) break;
                         for(int sy=-1;sy<=1;sy+=2)// Sign of translation
                         {
                            for(int ny=(sy+1)/2;;++ny)
                            {
                               const REAL y=y0 + sy*ny*m11 + sz*nz*m12;
                               if(abs(y)>maxDist) break;
                               for(int sx=-1;sx<=1;sx+=2)// Sign of translation
                               {
                                  for(int nx=(sx+1)/2;;++nx)
                                  {
                                     if((vUniqueIndex[i]==j) && (nx==0) && (ny==0) && (nz==0)) continue;// distance to self !
                                     const REAL x=x0 + sx*nx*m00 + sy*ny*m01 + sz*nz*m02;
                                     if(abs(x)>maxDist) break;
                                     const REAL d2=x*x+y*y+z*z;
                                     if(d2<=asymUnitMargin2)
                                     {
                                        Neighbour neigh(vPos[j].mAtomIndex,vPos[j].mSymmetryIndex,d2);
                                        vnb->push_back(neigh);
                                        #if 0
                                        if(!this->IsBeingRefined()) cout<<"    "<<vPos[j].mAtomIndex<<":"
                                              <<mScattCompList(vPos[j].mAtomIndex).mpScattPow->GetName()<<":"
                                              <<vPos[j].mSymmetryIndex<<":"
                                              <<FormatFloat(vPos[j].mX,8,5)<<","
                                              <<FormatFloat(vPos[j].mY,8,5)<<","<<","
                                              <<FormatFloat(vPos[j].mZ,8,5)<<","<<" vector="
                                              <<FormatFloat(x,8,5)<<","
                                              <<FormatFloat(y,8,5)<<","
                                              <<FormatFloat(z,8,5)<<":"<<sqrt(d2)<<","
                                              <<"("<<sx*nx<<","<<sy*ny<<","<<sz*nz<<")"<<endl;
                                        #endif
                                     }
                                  }
                               }
                            }
                         }
                      }
                   }
                }
                else
                {
                   const REAL d2=x0*x0+y0*y0+z0*z0;
                   if(d2<=asymUnitMargin2)
                   {
                      Neighbour neigh(vPos[j].mAtomIndex,vPos[j].mSymmetryIndex,d2);
                      vnb->push_back(neigh);
                      #if 0
                      if(!this->IsBeingRefined()) cout<<vPos[j].mAtomIndex<<":"
                            <<mScattCompList(vPos[j].mAtomIndex).mpScattPow->GetName()<<":"
                            <<vPos[j].mSymmetryIndex<<":"
                            <<vPos[j].mX<<","
                            <<vPos[j].mY<<","
                            <<vPos[j].mZ<<" vector="
                            <<x0<<","<<y0<<","<<z0<<":"<<sqrt(d2)<<","<<asymUnitMargin2<<endl;
                      #endif
                   }
                }

             }
          }
       }
    }
}
*/
void Crystal::CalcDistTableForInterMolDistCost() const
{
   this->GetScatteringComponentList();

   if(  (mDistTableForInterMolDistClock>mClockScattCompList)
      /*&&(mDistTableForInterMolDistClock>this->GetClockMetricMatrix())*/) return;
   

   cout<<"finding the maximal distance in the user-defined list\n";
   //finding the maximal distance in the user-defined list
   float maxDist = -1;
   for(int i=0;i<mInterMolDistList.size();i++) {
       if(mInterMolDistList[i].mDist2>0) {
           float d = sqrt(mInterMolDistList[i].mDist2);
           if(d>maxDist) {
               maxDist = d;
           }
       }
   }
   cout<<"maxDist = "<<maxDist<<"\n";
   if(maxDist<0) return;
   maxDist*=2;
   cout<<"maxDist = "<<maxDist<<"\n";

   if(true) {//calculate the table (imdTable) just for user-defined atoms 
    //this code is almost the same as in this->CalcDistTable(true);

       //this->GetScatteringComponentList();

       const long nbComponent=mScattCompList.GetNbComponent();

       // Get range and origin of the (pseudo) asymmetric unit
        const REAL asux0=this->GetSpaceGroup().GetAsymUnit().Xmin();
        const REAL asuy0=this->GetSpaceGroup().GetAsymUnit().Ymin();
        const REAL asuz0=this->GetSpaceGroup().GetAsymUnit().Zmin();

        const REAL asux1=this->GetSpaceGroup().GetAsymUnit().Xmax();
        const REAL asuy1=this->GetSpaceGroup().GetAsymUnit().Ymax();
        const REAL asuz1=this->GetSpaceGroup().GetAsymUnit().Zmax();

        const REAL halfasuxrange=(asux1-asux0)*0.5+1e-5;
        const REAL halfasuyrange=(asuy1-asuy0)*0.5+1e-5;
        const REAL halfasuzrange=(asuz1-asuz0)*0.5+1e-5;

        const REAL asuxc=0.5*(asux0+asux1);
        const REAL asuyc=0.5*(asuy0+asuy1);
        const REAL asuzc=0.5*(asuz0+asuz1);

        const REAL maxdx=halfasuxrange+maxDist/GetLatticePar(0);
        const REAL maxdy=halfasuyrange+maxDist/GetLatticePar(1);
        const REAL maxdz=halfasuzrange+maxDist/GetLatticePar(2);

       // List of all positions within or near the first atom generated
       std::vector<DistTableInternalPosition> vPos;
       // index of unique atoms in vPos, which are strictly in the asymmetric unit
       std::vector<unsigned long> vUniqueIndex(nbComponent);

       const REAL asymUnitMargin2 = maxDist*maxDist;

       {
          // No need to loop on a,b,c translations if maxDist is small enough
          bool loopOnLattice=true;
          if(  ((this->GetLatticePar(0)*.5)>maxDist)
             &&((this->GetLatticePar(1)*.5)>maxDist)
             &&((this->GetLatticePar(2)*.5)>maxDist)) loopOnLattice=false;

          
          //allocate the table 
          imdTable.resize(nbComponent);
          {
              std::vector<NeighbourHood>::iterator pos;
              for(pos=imdTable.begin();pos<imdTable.end();pos++)
                 pos->mvNeighbour.clear();
          }

          CrystMatrix_REAL symmetricsCoords;
          const int nbSymmetrics=this->GetSpaceGroup().GetNbSymmetrics(false,false);

          // Get the list of all atoms within or near the asymmetric unit
          for(long i=0;i<nbComponent;i++)
          {
             VFN_DEBUG_MESSAGE("Crystal::CalcDistTable(fast):3:component "<<i,0)
             // generate all symmetrics, excluding translations
             symmetricsCoords=this->GetSpaceGroup().GetAllSymmetrics(mScattCompList(i).mX,
                                                                     mScattCompList(i).mY,
                                                                     mScattCompList(i).mZ,
                                                                     false,false,false);
             imdTable[i].mIndex=i;//USELESS ?
             bool hasUnique=false;
             for(int j=0;j<nbSymmetrics;j++)
             {
                // take the closest position (using lattice translations) to the center of the ASU
                REAL x=fmod(symmetricsCoords(j,0)-asuxc,(REAL)1.0);if(x<-.5)x+=1;else if(x>.5)x-=1;
                REAL y=fmod(symmetricsCoords(j,1)-asuyc,(REAL)1.0);if(y<-.5)y+=1;else if(y>.5)y-=1;
                REAL z=fmod(symmetricsCoords(j,2)-asuzc,(REAL)1.0);if(z<-.5)z+=1;else if(z>.5)z-=1;

                //cout<<i<<","<<j<<":"<<FormatFloat(x,8,5)<<","<<FormatFloat(y,8,5)<<","<<FormatFloat(z,8,5)<<endl;
                if( (abs(x)<maxdx) && (abs(y)<maxdy) && (abs(z)<maxdz) )
                   vPos.push_back(DistTableInternalPosition(i, j, x+asuxc, y+asuyc, z+asuzc));
                // Get one reference atom strictly within the pseudo-ASU
                if(!hasUnique)
                   if( (abs(x)<halfasuxrange) && (abs(y)<halfasuyrange) && (abs(z)<halfasuzrange) )
                   {
                      hasUnique=true;
                      vUniqueIndex[i]=vPos.size()-1;
                      imdTable[i].mUniquePosSymmetryIndex=j;
                   }
             }
             if(!hasUnique)
             {
                throw ObjCrystException("One atom did not have any symmetric in the ASU !");
             }
          }

          const CrystMatrix_REAL* pOrthMatrix=&(this->GetOrthMatrix());

          const REAL m00=(*pOrthMatrix)(0,0);
          const REAL m01=(*pOrthMatrix)(0,1);
          const REAL m02=(*pOrthMatrix)(0,2);
          const REAL m11=(*pOrthMatrix)(1,1);
          const REAL m12=(*pOrthMatrix)(1,2);
          const REAL m22=(*pOrthMatrix)(2,2);

          cout<<"Get the list of all atoms in the asymmetric unit in mInterMolDistList\n";
          // Get the list of all atoms in the asymmetric unit in mInterMolDistList
          vector<int> listOfCurrentScatterers;
          for(long i=0;i<mInterMolDistList.size();i++) {
              vector<int> p = FindScatterersInComponentList(mInterMolDistList[i].mAt1);
              listOfCurrentScatterers.insert(listOfCurrentScatterers.end(), p.begin(), p.end());
          }
          if(listOfCurrentScatterers.size()==0) return;
          cout<<"listOfCurrentScatterers="<<listOfCurrentScatterers.size()<<"\n";

          // Get the list of all mAt2 atoms in mInterMolDistList
          vector<int> listOfNeigScatterers;
          for(long i=0;i<mInterMolDistList.size();i++) {
              vector<int> p = FindScatterersInComponentList(mInterMolDistList[i].mAt2);
              listOfNeigScatterers.insert(listOfNeigScatterers.end(), p.begin(), p.end());
          }
          if(listOfNeigScatterers.size()==0) return;
          cout<<"listOfNeigScatterers="<<listOfNeigScatterers.size()<<"\n";

          for(long i=0;i<nbComponent;i++)
          {
             //calculate the table just for user-define atoms mAt1
             bool found = false;
             for(int q=0;q<listOfCurrentScatterers.size();q++) {
                 if(i==listOfCurrentScatterers[q]) {
                     found = true;
                     break;
                 }
             }
             if(!found) continue;

             VFN_DEBUG_MESSAGE("Crystal::CalcDistTable(fast):4:component "<<i,0)
             #if 0
             if(!this->IsBeingRefined()) cout<<endl<<"Unique pos:"<<vUniqueIndex[i]<<":"
                 <<vPos[vUniqueIndex[i]].mAtomIndex<<":"
                 <<mScattCompList(vPos[vUniqueIndex[i]].mAtomIndex).mpScattPow->GetName()<<":"
                 <<vPos[vUniqueIndex[i]].mSymmetryIndex<<":"
                 <<FormatFloat(vPos[vUniqueIndex[i]].mX,8,5)<<","
                 <<FormatFloat(vPos[vUniqueIndex[i]].mY,8,5)<<","
                 <<FormatFloat(vPos[vUniqueIndex[i]].mZ,8,5)<<endl;
             #endif
             std::vector<Crystal::Neighbour> * const vnb=&(imdTable[i].mvNeighbour);
             const REAL x0i=vPos[vUniqueIndex[i] ].mX;
             const REAL y0i=vPos[vUniqueIndex[i] ].mY;
             const REAL z0i=vPos[vUniqueIndex[i] ].mZ;
             for(unsigned long j=0;j<vPos.size();j++)
             {
                
                 //Calculate it just for user-defied atoms mAt2
                found = false;
                for(int q=0;q<listOfNeigScatterers.size();q++) {
                     if(vPos[j].mAtomIndex==listOfNeigScatterers[q]) {
                         found = true;
                         break;
                     }
                }
                if(!found) continue;

                if((vUniqueIndex[i]==j) && (!loopOnLattice)) continue;// distance to self !
                // Start with the smallest absolute coordinates possible
                REAL x=fmod(vPos[j].mX - x0i,(REAL)1.0);if(x<-.5)x+=1;if(x>.5)x-=1;
                REAL y=fmod(vPos[j].mY - y0i,(REAL)1.0);if(y<-.5)y+=1;if(y>.5)y-=1;
                REAL z=fmod(vPos[j].mZ - z0i,(REAL)1.0);if(z<-.5)z+=1;if(z>.5)z-=1;

                const REAL x0=m00 * x + m01 * y + m02 * z;
                const REAL y0=          m11 * y + m12 * z;
                const REAL z0=                    m22 * z;

                if(loopOnLattice)// distance to self !
                {//Now loop over lattice translations
                   for(int sz=-1;sz<=1;sz+=2)// Sign of translation
                   {
                      for(int nz=(sz+1)/2;;++nz)
                      {
                         const REAL z=z0+sz*nz*m22;
                         if(abs(z)>maxDist) break;
                         for(int sy=-1;sy<=1;sy+=2)// Sign of translation
                         {
                            for(int ny=(sy+1)/2;;++ny)
                            {
                               const REAL y=y0 + sy*ny*m11 + sz*nz*m12;
                               if(abs(y)>maxDist) break;
                               for(int sx=-1;sx<=1;sx+=2)// Sign of translation
                               {
                                  for(int nx=(sx+1)/2;;++nx)
                                  {
                                     if((vUniqueIndex[i]==j) && (nx==0) && (ny==0) && (nz==0)) continue;// distance to self !
                                     const REAL x=x0 + sx*nx*m00 + sy*ny*m01 + sz*nz*m02;
                                     if(abs(x)>maxDist) break;
                                     const REAL d2=x*x+y*y+z*z;
                                     if(d2<=asymUnitMargin2)
                                     {
                                        Neighbour neigh(vPos[j].mAtomIndex,vPos[j].mSymmetryIndex,d2);
                                        vnb->push_back(neigh);
                                        #if 0
                                        if(!this->IsBeingRefined()) cout<<"    "<<vPos[j].mAtomIndex<<":"
                                              <<mScattCompList(vPos[j].mAtomIndex).mpScattPow->GetName()<<":"
                                              <<vPos[j].mSymmetryIndex<<":"
                                              <<FormatFloat(vPos[j].mX,8,5)<<","
                                              <<FormatFloat(vPos[j].mY,8,5)<<","<<","
                                              <<FormatFloat(vPos[j].mZ,8,5)<<","<<" vector="
                                              <<FormatFloat(x,8,5)<<","
                                              <<FormatFloat(y,8,5)<<","
                                              <<FormatFloat(z,8,5)<<":"<<sqrt(d2)<<","
                                              <<"("<<sx*nx<<","<<sy*ny<<","<<sz*nz<<")"<<endl;
                                        #endif
                                     }
                                  }
                               }
                            }
                         }
                      }
                   }
                }
                else
                {
                   const REAL d2=x0*x0+y0*y0+z0*z0;
                   if(d2<=asymUnitMargin2)
                   {
                      Neighbour neigh(vPos[j].mAtomIndex,vPos[j].mSymmetryIndex,d2);
                      vnb->push_back(neigh);
                      #if 0
                      if(!this->IsBeingRefined()) cout<<vPos[j].mAtomIndex<<":"
                            <<mScattCompList(vPos[j].mAtomIndex).mpScattPow->GetName()<<":"
                            <<vPos[j].mSymmetryIndex<<":"
                            <<vPos[j].mX<<","
                            <<vPos[j].mY<<","
                            <<vPos[j].mZ<<" vector="
                            <<x0<<","<<y0<<","<<z0<<":"<<sqrt(d2)<<","<<asymUnitMargin2<<endl;
                      #endif
                   }
                }

             }
          }
       }
    }
    mDistTableForInterMolDistClock.Click();
}
void Crystal::CalcDistTable(const bool fast) const
{
   this->GetScatteringComponentList();
   if(!this->IsBeingRefined())
   {
      if(mDistTableMaxDistance!=10) mDistTableClock.Reset();
      mDistTableMaxDistance=10;
   }

   if(  (mDistTableClock>mClockScattCompList)
      &&(mDistTableClock>this->GetClockMetricMatrix())) return;
   VFN_DEBUG_ENTRY("Crystal::CalcDistTable(fast="<<fast<<"),maxDist="<<mDistTableMaxDistance,4)

   const long nbComponent=mScattCompList.GetNbComponent();

   mvDistTableSq.resize(nbComponent);
   {
      std::vector<NeighbourHood>::iterator pos;
      for(pos=mvDistTableSq.begin();pos<mvDistTableSq.end();pos++)
         pos->mvNeighbour.clear();
   }
   VFN_DEBUG_MESSAGE("Crystal::CalcDistTable():1",3)

   // Get range and origin of the (pseudo) asymmetric unit
      const REAL asux0=this->GetSpaceGroup().GetAsymUnit().Xmin();
      const REAL asuy0=this->GetSpaceGroup().GetAsymUnit().Ymin();
      const REAL asuz0=this->GetSpaceGroup().GetAsymUnit().Zmin();

      const REAL asux1=this->GetSpaceGroup().GetAsymUnit().Xmax();
      const REAL asuy1=this->GetSpaceGroup().GetAsymUnit().Ymax();
      const REAL asuz1=this->GetSpaceGroup().GetAsymUnit().Zmax();

      const REAL halfasuxrange=(asux1-asux0)*0.5+1e-5;
      const REAL halfasuyrange=(asuy1-asuy0)*0.5+1e-5;
      const REAL halfasuzrange=(asuz1-asuz0)*0.5+1e-5;

      const REAL asuxc=0.5*(asux0+asux1);
      const REAL asuyc=0.5*(asuy0+asuy1);
      const REAL asuzc=0.5*(asuz0+asuz1);

      const REAL maxdx=halfasuxrange+mDistTableMaxDistance/GetLatticePar(0);
      const REAL maxdy=halfasuyrange+mDistTableMaxDistance/GetLatticePar(1);
      const REAL maxdz=halfasuzrange+mDistTableMaxDistance/GetLatticePar(2);

   // List of all positions within or near the first atom generated
   std::vector<DistTableInternalPosition> vPos;
   // index of unique atoms in vPos, which are strictly in the asymmetric unit
   std::vector<unsigned long> vUniqueIndex(nbComponent);

   const REAL asymUnitMargin2 = mDistTableMaxDistance*mDistTableMaxDistance;

   if(true)//(true==fast)
   {
      VFN_DEBUG_MESSAGE("Crystal::CalcDistTable(fast):2",3)
      TAU_PROFILE("Crystal::CalcDistTable(fast=true)","Matrix (string&)",TAU_DEFAULT);
      TAU_PROFILE_TIMER(timer1,"DiffractionData::CalcDistTable1","", TAU_FIELD);
      TAU_PROFILE_TIMER(timer2,"DiffractionData::CalcDistTable2","", TAU_FIELD);

      TAU_PROFILE_START(timer1);

      // No need to loop on a,b,c translations if mDistTableMaxDistance is small enough
      bool loopOnLattice=true;
      if(  ((this->GetLatticePar(0)*.5)>mDistTableMaxDistance)
         &&((this->GetLatticePar(1)*.5)>mDistTableMaxDistance)
         &&((this->GetLatticePar(2)*.5)>mDistTableMaxDistance)) loopOnLattice=false;

      CrystMatrix_REAL symmetricsCoords;

      const int nbSymmetrics=this->GetSpaceGroup().GetNbSymmetrics(false,false);

      // Get the list of all atoms within or near the asymmetric unit
      for(long i=0;i<nbComponent;i++)
      {
         VFN_DEBUG_MESSAGE("Crystal::CalcDistTable(fast):3:component "<<i,0)
         // generate all symmetrics, excluding translations
         symmetricsCoords=this->GetSpaceGroup().GetAllSymmetrics(mScattCompList(i).mX,
                                                                 mScattCompList(i).mY,
                                                                 mScattCompList(i).mZ,
                                                                 false,false,false);
         mvDistTableSq[i].mIndex=i;//USELESS ?
         bool hasUnique=false;
         for(int j=0;j<nbSymmetrics;j++)
         {
            // take the closest position (using lattice translations) to the center of the ASU
            REAL x=fmod(symmetricsCoords(j,0)-asuxc,(REAL)1.0);if(x<-.5)x+=1;else if(x>.5)x-=1;
            REAL y=fmod(symmetricsCoords(j,1)-asuyc,(REAL)1.0);if(y<-.5)y+=1;else if(y>.5)y-=1;
            REAL z=fmod(symmetricsCoords(j,2)-asuzc,(REAL)1.0);if(z<-.5)z+=1;else if(z>.5)z-=1;

            //cout<<i<<","<<j<<":"<<FormatFloat(x,8,5)<<","<<FormatFloat(y,8,5)<<","<<FormatFloat(z,8,5)<<endl;
            if( (abs(x)<maxdx) && (abs(y)<maxdy) && (abs(z)<maxdz) )
               vPos.push_back(DistTableInternalPosition(i, j, x+asuxc, y+asuyc, z+asuzc));
            // Get one reference atom strictly within the pseudo-ASU
            if(!hasUnique)
               if( (abs(x)<halfasuxrange) && (abs(y)<halfasuyrange) && (abs(z)<halfasuzrange) )
               {
                  hasUnique=true;
                  vUniqueIndex[i]=vPos.size()-1;
                  mvDistTableSq[i].mUniquePosSymmetryIndex=j;
               }
         }
         if(!hasUnique)
         {
            throw ObjCrystException("One atom did not have any symmetric in the ASU !");
         }
      }
      TAU_PROFILE_STOP(timer1);
      TAU_PROFILE_START(timer2);
      // Compute interatomic vectors & distance
      // between (i) unique atoms and (ii) all remaining atoms

      const CrystMatrix_REAL* pOrthMatrix=&(this->GetOrthMatrix());

      const REAL m00=(*pOrthMatrix)(0,0);
      const REAL m01=(*pOrthMatrix)(0,1);
      const REAL m02=(*pOrthMatrix)(0,2);
      const REAL m11=(*pOrthMatrix)(1,1);
      const REAL m12=(*pOrthMatrix)(1,2);
      const REAL m22=(*pOrthMatrix)(2,2);

      for(long i=0;i<nbComponent;i++)
      {
         VFN_DEBUG_MESSAGE("Crystal::CalcDistTable(fast):4:component "<<i,0)
         #if 0
         if(!this->IsBeingRefined()) cout<<endl<<"Unique pos:"<<vUniqueIndex[i]<<":"
             <<vPos[vUniqueIndex[i]].mAtomIndex<<":"
             <<mScattCompList(vPos[vUniqueIndex[i]].mAtomIndex).mpScattPow->GetName()<<":"
             <<vPos[vUniqueIndex[i]].mSymmetryIndex<<":"
             <<FormatFloat(vPos[vUniqueIndex[i]].mX,8,5)<<","
             <<FormatFloat(vPos[vUniqueIndex[i]].mY,8,5)<<","
             <<FormatFloat(vPos[vUniqueIndex[i]].mZ,8,5)<<endl;
         #endif
         std::vector<Crystal::Neighbour> * const vnb=&(mvDistTableSq[i].mvNeighbour);
         const REAL x0i=vPos[vUniqueIndex[i] ].mX;
         const REAL y0i=vPos[vUniqueIndex[i] ].mY;
         const REAL z0i=vPos[vUniqueIndex[i] ].mZ;
         for(unsigned long j=0;j<vPos.size();j++)
         {
            if((vUniqueIndex[i]==j) && (!loopOnLattice)) continue;// distance to self !
            // Start with the smallest absolute coordinates possible
            REAL x=fmod(vPos[j].mX - x0i,(REAL)1.0);if(x<-.5)x+=1;if(x>.5)x-=1;
            REAL y=fmod(vPos[j].mY - y0i,(REAL)1.0);if(y<-.5)y+=1;if(y>.5)y-=1;
            REAL z=fmod(vPos[j].mZ - z0i,(REAL)1.0);if(z<-.5)z+=1;if(z>.5)z-=1;

            const REAL x0=m00 * x + m01 * y + m02 * z;
            const REAL y0=          m11 * y + m12 * z;
            const REAL z0=                    m22 * z;

            if(loopOnLattice)// distance to self !
            {//Now loop over lattice translations
               for(int sz=-1;sz<=1;sz+=2)// Sign of translation
               {
                  for(int nz=(sz+1)/2;;++nz)
                  {
                     const REAL z=z0+sz*nz*m22;
                     if(abs(z)>mDistTableMaxDistance) break;
                     for(int sy=-1;sy<=1;sy+=2)// Sign of translation
                     {
                        for(int ny=(sy+1)/2;;++ny)
                        {
                           const REAL y=y0 + sy*ny*m11 + sz*nz*m12;
                           if(abs(y)>mDistTableMaxDistance) break;
                           for(int sx=-1;sx<=1;sx+=2)// Sign of translation
                           {
                              for(int nx=(sx+1)/2;;++nx)
                              {
                                 if((vUniqueIndex[i]==j) && (nx==0) && (ny==0) && (nz==0)) continue;// distance to self !
                                 const REAL x=x0 + sx*nx*m00 + sy*ny*m01 + sz*nz*m02;
                                 if(abs(x)>mDistTableMaxDistance) break;
                                 const REAL d2=x*x+y*y+z*z;
                                 if(d2<=asymUnitMargin2)
                                 {
                                    Neighbour neigh(vPos[j].mAtomIndex,vPos[j].mSymmetryIndex,d2);
                                    vnb->push_back(neigh);
                                    #if 0
                                    if(!this->IsBeingRefined()) cout<<"    "<<vPos[j].mAtomIndex<<":"
                                          <<mScattCompList(vPos[j].mAtomIndex).mpScattPow->GetName()<<":"
                                          <<vPos[j].mSymmetryIndex<<":"
                                          <<FormatFloat(vPos[j].mX,8,5)<<","
                                          <<FormatFloat(vPos[j].mY,8,5)<<","<<","
                                          <<FormatFloat(vPos[j].mZ,8,5)<<","<<" vector="
                                          <<FormatFloat(x,8,5)<<","
                                          <<FormatFloat(y,8,5)<<","
                                          <<FormatFloat(z,8,5)<<":"<<sqrt(d2)<<","
                                          <<"("<<sx*nx<<","<<sy*ny<<","<<sz*nz<<")"<<endl;
                                    #endif
                                 }
                              }
                           }
                        }
                     }
                  }
               }
            }
            else
            {
               const REAL d2=x0*x0+y0*y0+z0*z0;
               if(d2<=asymUnitMargin2)
               {
                  Neighbour neigh(vPos[j].mAtomIndex,vPos[j].mSymmetryIndex,d2);
                  vnb->push_back(neigh);
                  #if 0
                  if(!this->IsBeingRefined()) cout<<vPos[j].mAtomIndex<<":"
                        <<mScattCompList(vPos[j].mAtomIndex).mpScattPow->GetName()<<":"
                        <<vPos[j].mSymmetryIndex<<":"
                        <<vPos[j].mX<<","
                        <<vPos[j].mY<<","
                        <<vPos[j].mZ<<" vector="
                        <<x0<<","<<y0<<","<<z0<<":"<<sqrt(d2)<<","<<asymUnitMargin2<<endl;
                  #endif
               }
            }

         }
      }
      TAU_PROFILE_STOP(timer2);
   }
   mDistTableClock.Click();
   VFN_DEBUG_EXIT("Crystal::CalcDistTable()",4)
}

void Crystal::SetDeleteSubObjInDestructor(const bool b) {
    mDeleteSubObjInDestructor=b;
}

#ifdef __WX__CRYST__
WXCrystObjBasic* Crystal::WXCreate(wxWindow* parent)
{
   VFN_DEBUG_ENTRY("Crystal::WXCreate(wxWindow*)",6)
   //:TODO: Check mpWXCrystObj==0
   mpWXCrystObj=new WXCrystal(parent,this);
   VFN_DEBUG_EXIT("Crystal::WXCreate(wxWindow*)",6)
   return mpWXCrystObj;
}
#endif

}//namespace
