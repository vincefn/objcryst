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

#include <typeinfo>

#include "ObjCryst/Crystal.h"

#include "Quirks/VFNStreamFormat.h" //simple formatting of integers, REALs..
#include "Quirks/VFNDebug.h"

#ifdef OBJCRYST_GL
#include <GL/glut.h>
#endif
#ifdef __WX__CRYST__
   #include "wxCryst/wxCrystal.h"
#endif

#include <fstream>
#include <iomanip>

namespace ObjCryst
{
const RefParType *gpRefParTypeCrystal=
   new RefParType (gpRefParTypeObjCryst,"Crystal");

////////////////////////////////////////////////////////////////////////
//
//    CRYSTAL : the crystal (Unit cell, spaceGroup, scatterers)
//
////////////////////////////////////////////////////////////////////////
ObjRegistry<Crystal> gCrystalRegistry("List of all Crystals");

Crystal::Crystal():
mScattererRegistry("List of Crystal Scatterers"),
mScatteringPowerRegistry("List of Crystal ScatteringPowers")
{
   VFN_DEBUG_MESSAGE("Crystal::Crystal()",10)
   this->InitOptions();
   this->Init(10,11,12,M_PI/2+.1,M_PI/2+.2,M_PI/2+.3,"P1","");
   gCrystalRegistry.Register(*this);
   gTopRefinableObjRegistry.Register(*this);
}

Crystal::Crystal(const REAL a, const REAL b, const REAL c, const string &SpaceGroupId):
mScattererRegistry("List of Crystal Scatterers"),
mScatteringPowerRegistry("List of Crystal ScatteringPowers")
{
   VFN_DEBUG_MESSAGE("Crystal::Crystal(a,b,c,Sg)",10)
   this->Init(a,b,c,M_PI/2,M_PI/2,M_PI/2,SpaceGroupId,"");
   this->InitOptions();
   gCrystalRegistry.Register(*this);
   gTopRefinableObjRegistry.Register(*this);
}

Crystal::Crystal(const REAL a, const REAL b, const REAL c, const REAL alpha,
              const REAL beta, const REAL gamma,const string &SpaceGroupId):
mScattererRegistry("List of Crystal Scatterers"),
mScatteringPowerRegistry("List of Crystal ScatteringPowers")
{
   VFN_DEBUG_MESSAGE("Crystal::Crystal(a,b,c,alpha,beta,gamma,Sg)",10)
   this->Init(a,b,c,alpha,beta,gamma,SpaceGroupId,"");
   this->InitOptions();
   gCrystalRegistry.Register(*this);
   gTopRefinableObjRegistry.Register(*this);
}

Crystal::Crystal(const Crystal &old):
mScattererRegistry(old.mScattererRegistry),
mScatteringPowerRegistry(old.mScatteringPowerRegistry)
{
   VFN_DEBUG_MESSAGE("Crystal::Crystal(&oldCrystal)",10)
   for(long i=0;i<old.GetNbScatterer();i++)
   {
      mScattererRegistry.Register(*(old.GetScatt(i).CreateCopy()));
   }
   
   mUseDynPopCorr.SetChoice(old.mUseDynPopCorr.GetChoice());
   mDisplayEnantiomer.SetChoice(old.mDisplayEnantiomer.GetChoice());
   
   gCrystalRegistry.Register(*this);
   gTopRefinableObjRegistry.Register(*this);
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
   mScattererRegistry.DeleteAll();
   for(long i=0;i<mScatteringPowerRegistry.GetNb();i++)
   {
      VFN_DEBUG_MESSAGE("Crystal::~Crystal(&scatt):2:"<<i,5)
      this->RemoveSubRefObj(mScatteringPowerRegistry.GetObj(i));
      mScatteringPowerRegistry.GetObj(i).DeRegisterClient(*this);
      // :TODO: check if it is not used by another Crystal (forbidden!)
   }
   mScatteringPowerRegistry.DeleteAll();
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

void Crystal::RemoveScatterer(Scatterer *scatt)
{
   VFN_DEBUG_MESSAGE("Crystal::RemoveScatterer(&scatt)",5)
   mScattererRegistry.DeRegister(*scatt);
   scatt->DeRegisterClient(*this);
   this->RemoveSubRefObj(*scatt);
   delete scatt;
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

ObjRegistry<ScatteringPower>& Crystal::GetScatteringPowerRegistry() 
{return mScatteringPowerRegistry;}
const ObjRegistry<ScatteringPower>& Crystal::GetScatteringPowerRegistry() const
{return mScatteringPowerRegistry;}

void Crystal::AddScatteringPower(ScatteringPower *scattPow)
{
   mScatteringPowerRegistry.Register(*scattPow);
   scattPow->RegisterClient(*this);//:TODO: Should register as (unique) 'owner'.
   this->AddSubRefObj(*scattPow);
}

void Crystal::RemoveScatteringPower(ScatteringPower *scattPow)
{
   VFN_DEBUG_ENTRY("Crystal::RemoveScatteringPower()",2)
   mScatteringPowerRegistry.DeRegister(*scattPow);
   this->RemoveSubRefObj(*scattPow);
   delete scattPow;
   VFN_DEBUG_EXIT("Crystal::RemoveScatteringPower()",2)
}

ScatteringPower& Crystal::GetScatteringPower(const string &name)
{
   return mScatteringPowerRegistry.GetObj(name);
}

const ScatteringComponentList& Crystal::GetScatteringComponentList()const
{
   //:TODO: only update when necessary..
   bool update=false;
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
         this->CalcDynPopCorr(1.,.1); else this->ResetDynPopCorr();
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
   
   os << "List of scattering components (atoms): " << mScattCompList.GetNbComponent() << endl ;
   
   long k=0;
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
              << " * " << FormatFloat(mScattCompList(k).mDynPopCorr,6,4)
              << " ,ScattPow:" << FormatString(list(j).mpScattPow->GetName(),16)
              << ", Biso=" << FormatFloat(list(j).mpScattPow->GetBiso())
              << endl;
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
   os << " Total number of components (atoms) in one unit cell : " << nbAtoms<<endl<<endl;
   
   VFN_DEBUG_MESSAGE("Crystal::Print():End",5)
}

CrystMatrix_REAL Crystal::GetMinDistanceTable(const REAL minDistance) const
{
   VFN_DEBUG_MESSAGE("Crystal::MinDistanceTable()",5)
   this->CalcDistTable(true);
   const long nbComponent=mScattCompList.GetNbComponent();
   
   CrystMatrix_REAL minDistTable(nbComponent,nbComponent);
   REAL dist;
   REAL tmp;
   const REAL min=minDistance*minDistance;
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
         if(9999.<minDistTable(i,j)) minDistTable(i,j)=0;
         minDistTable(i,j)=sqrt(minDistTable(i,j));
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
            os << FormatFloat(minDistTable(l,k),6,3) ;
         os << endl;
         l++;
      }
   VFN_DEBUG_MESSAGE("Crystal::PrintMinDistanceTable():End",3)
}

ostream& Crystal::POVRayDescription(ostream &os,bool onlyIndependentAtoms)const
{
   VFN_DEBUG_MESSAGE("Crystal::POVRayDescription(os,bool)",5)
   os << "#include \"colours.inc\"  " <<endl;
   os << "// Description of Crystal :" << mName <<endl;
   os << "global_settings { assumed_gamma 2.2 ambient_light rgb <2,2,2>}"<<endl;
   os << "camera" <<endl;
   os << "{"<<endl;
   os << "    location  <200, 0, 0>"<<endl;
   os << "    rotate  <0, -30, -30-clock>" <<endl;
   os << "    sky   <0, 0, 1>"<<endl;
   REAL x=0.5;
   REAL y=0.5;
   REAL z=0.5;
   this->FractionalToOrthonormalCoords(x,y,z);
   os << "    look_at   <" << x << "," << y << "," << z <<">"<<endl;
   REAL maxDim=x;
   if(y>maxDim) maxDim=y;
   if(z>maxDim) maxDim=z;
   os << "    angle   "<< 2.5*atan(2*maxDim/200.)*RAD2DEG <<endl;
   os << "}"<<endl;
   os << "light_source"<<endl;
   os << "{" <<endl;
   os << "   <100, 0, 0> " <<endl;
   os << "   colour White" <<endl;
   os << "   shadowless" <<endl;
   os << "   rotate  <0, -60, -10-clock>" <<endl;
   os << "}" <<endl;
   os << "background { colour rgb <0.0, 0.0, 0.0> }"<<endl;
   
   os << endl <<"//Crystal Unit Cell" <<endl;
   x=1;
   y=1;
   z=1;
   this->FractionalToOrthonormalCoords(x,y,z);
   os << "   box{ <0,0,0>, <"<< x << "," << y << "," << z << ">" <<endl;
   os << "         pigment {colour rgbf<1,1,0.9,0.5>}" << endl;
   os << "         hollow" << endl;
   os << "   }" <<endl<<endl;
   
  for(int i=0;i<mScattererRegistry.GetNb();i++)
      this->GetScatt(i).POVRayDescription(os,onlyIndependentAtoms) ;
   return os;
}

void Crystal::GLInitDisplayList(const bool onlyIndependentAtoms,
                                const REAL xMin,const REAL xMax,
                                const REAL yMin,const REAL yMax,
                                const REAL zMin,const REAL zMax)const
{
   VFN_DEBUG_ENTRY("Crystal::GLInitDisplayList()",5)
   #ifdef OBJCRYST_GL
      REAL en=1;// if -1, display enantiomeric structure
      if(mDisplayEnantiomer.GetChoice()==1) en=-1;
      
      const GLfloat colorAmbient [] = {0.50, 0.50, 0.50, 1.00}; 
      const GLfloat colorDiffuse [] = {0.80, 0.80, 0.80, 1.00}; 
      const GLfloat colorSpecular [] = {1.00, 1.00, 1.00, 1.00}; 
      
      glMaterialfv(GL_FRONT, GL_AMBIENT,   colorAmbient); 
      glMaterialfv(GL_FRONT, GL_DIFFUSE,   colorDiffuse); 
      glMaterialfv(GL_FRONT, GL_SPECULAR,  colorSpecular); 
      glMaterialf( GL_FRONT, GL_SHININESS, 5.0); 
   //cout << xMin << ":"<<xMax <<endl;
   //cout << yMin << ":"<<yMax <<endl;
   //cout << zMin << ":"<<zMax <<endl;
   //Center of displayed unit
      REAL xc=(xMin+xMax)/2.;
      REAL yc=(yMin+yMax)/2.;
      REAL zc=(zMin+zMax)/2.;
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
         REAL x,y,z;

         x=1.2-xc;y=-yc;z=-zc;
         this->FractionalToOrthonormalCoords(x,y,z);
         glRasterPos3f(en*x,y,z);
         glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18,'a');

         x=-xc;y=1.2-yc;z=-zc;
         this->FractionalToOrthonormalCoords(x,y,z);
         glRasterPos3f(en*x,y,z);
         glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18,'b');

         x=-xc;y=-yc;z=1.2-zc;
         this->FractionalToOrthonormalCoords(x,y,z);
         glRasterPos3f(en*x,y,z);
         glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18,'c');
      // Cell
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


      //Describe all Scatterers
      VFN_DEBUG_MESSAGE("Crystal::GLView(bool):Scatterers...",5)
      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
      {
         bool displayEnantiomer=false;
         if(mDisplayEnantiomer.GetChoice()==1) displayEnantiomer=true;
         for(int i=0;i<mScattererRegistry.GetNb();i++) 
            this->GetScatt(i).GLInitDisplayList(onlyIndependentAtoms,
                                                xMin,xMax,yMin,yMax,zMin,zMax,
                                                displayEnantiomer);
      }
   glPopMatrix();
   #else
   cout << "Crystal::GLView(): Compiled without OpenGL support !" <<endl;
   #endif
   VFN_DEBUG_EXIT("Crystal::GLInitDisplayList(bool)",5)
}

void Crystal::CalcDynPopCorr(const REAL overlapDist, const REAL mergeDist) const
{
   VFN_DEBUG_ENTRY("Crystal::CalcDynPopCorr(REAL)",4)
   TAU_PROFILE("Crystal::CalcDynPopCorr()","void (REAL)",TAU_DEFAULT);
   
   this->CalcDistTable(true,overlapDist);
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
      atomicNumber=mScattCompList(i).mpScattPow->GetDynPopCorrIndex();
      nbNeighbors=0;
      std::vector<Crystal::Neighbour>::const_iterator pos;
      for(pos=mvDistTableSq[i].mvNeighbour.begin();
          pos<mvDistTableSq[i].mvNeighbour.end();pos++)
      {
         VFN_DEBUG_MESSAGE("Crystal::CalcDynPopCorr(): Component:"<<i<<"Neighbour:"<<pos->mNeighbourIndex,0)
         if(atomicNumber==mScattCompList(pos->mNeighbourIndex).mpScattPow->GetDynPopCorrIndex())
         {
            if(overlapDistSq > pos->mDist2)
               neighborsDist(nbNeighbors++)=sqrt(pos->mDist2);
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

void Crystal::SetUseDynPopCorr(const int b)
{
   VFN_DEBUG_MESSAGE("Crystal::SetUseDynPopCorr()",1)
   mUseDynPopCorr.SetChoice(b);
   mClockDynPopCorr.Reset();
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

Crystal::BumpMergePar::BumpMergePar():
   mDist2(1.),mCanOverlap(false){}

Crystal::BumpMergePar::BumpMergePar(const REAL dist, const bool canOverlap):
   mDist2(dist*dist),mCanOverlap(canOverlap){}

REAL Crystal::GetBumpMergeCost() const
{
   if(mvBumpMergePar.size()==0) return 0;
   this->CalcDistTable(true,3);
   VFN_DEBUG_ENTRY("Crystal::GetBumpMergeCost()",4)
   if(  (mBumpMergeCostClock>mBumpMergeParClock)
      &&(mBumpMergeCostClock>mDistTableClock)) return mBumpMergeCost;
   TAU_PROFILE("Crystal::GetBumpMergeCost()","REAL (REAL)",TAU_DEFAULT);
   
   mBumpMergeCost=0;
   
   // we don't want to search the corresponding antibump distance for each
   // couple of neighbours, so build an index table
      const unsigned int nbScattPow=this->GetScatteringPowerRegistry().GetNb();
      long maxGetDynPopCorrIndex=0;
      for(unsigned int i=0;i<nbScattPow;i++)
         if(this->GetScatteringPowerRegistry().GetObj(i).GetDynPopCorrIndex()>maxGetDynPopCorrIndex)
            maxGetDynPopCorrIndex=this->GetScatteringPowerRegistry().GetObj(i).GetDynPopCorrIndex();
      CrystVector_long index(maxGetDynPopCorrIndex+1);
      for(unsigned int i=0;i<nbScattPow;i++)
         index(this->GetScatteringPowerRegistry().GetObj(i).GetDynPopCorrIndex())=i;
      CrystMatrix_long paridx(nbScattPow,nbScattPow);
      paridx=-1;
      {
         for(unsigned int i=0;i<mvBumpMergePar.size();i++)
         {
            paridx(index(mvBumpMergePar[i].first.first ->GetDynPopCorrIndex()),
                   index(mvBumpMergePar[i].first.second->GetDynPopCorrIndex()))=i;
         }
      }
      
   std::vector<NeighbourHood>::const_iterator pos;
   std::vector<Crystal::Neighbour>::const_iterator neigh;
   REAL tmp;
   long par;
   const ScatteringPower *pow1;
   const ScatteringPower *pow2;
   for(pos=mvDistTableSq.begin();pos<mvDistTableSq.end();pos++)
   {
      pow1=mScattCompList(pos->mIndex).mpScattPow;
      for(neigh=pos->mvNeighbour.begin();neigh<pos->mvNeighbour.end();neigh++)
      {
         pow2=mScattCompList(neigh->mNeighbourIndex).mpScattPow;
         par=paridx(index( pow1->GetDynPopCorrIndex() ),
                    index( pow2->GetDynPopCorrIndex() ));
         if(-1==par) continue;
         if(neigh->mDist2 > mvBumpMergePar[par].second.mDist2) continue;
         if(true==mvBumpMergePar[par].second.mCanOverlap)
            tmp = 0.5*sin(M_PI*(1.-sqrt(neigh->mDist2/mvBumpMergePar[par].second.mDist2)))/0.1;
         else
            tmp = tan(M_PI*0.5*(1.-sqrt(neigh->mDist2/mvBumpMergePar[par].second.mDist2)))/0.1;
         mBumpMergeCost += tmp*tmp;
      }
   }
   mBumpMergeCost *= this->GetSpaceGroup().GetNbSymmetrics();
   mBumpMergeCostClock.Click();
   VFN_DEBUG_EXIT("Crystal::GetBumpMergeCost():"<<mBumpMergeCost,4)
   return mBumpMergeCost;
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
   VFN_DEBUG_MESSAGE("Crystal::SetBumpMergeDistance()",8)
   mvBumpMergePar.push_back(std::make_pair(std::make_pair(&scatt1,&scatt2),
                                           BumpMergePar(dist,allowMerge)));
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
   return this->GetBumpMergeCost();
}

void Crystal::CIFOutput(ostream &os)const
{
   VFN_DEBUG_ENTRY("Crystal::OutputCIF()",5)
   //Program
   os <<"_computing_structure_solution     'FOX http://objcryst.sourceforge.net'"<<endl<<endl;

   //Scattering powers
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
   
   //Symmetry
   os <<"_symmetry_space_group_name_H-M    "
      << this->GetSpaceGroup().GetHM_as_Hall().HM<<""<<endl;
   os <<"_symmetry_space_group_name_Hall   '"
      << this->GetSpaceGroup().GetHM_as_Hall().Hall<<"'"<<endl;
   os <<endl;
      
   os << "_cell_length_a "   << FormatFloat(this->GetLatticePar(0),8,5) << endl
      << "_cell_length_b "   << FormatFloat(this->GetLatticePar(1),8,5) << endl
      << "_cell_length_c "   << FormatFloat(this->GetLatticePar(2),8,5) << endl
      << "_cell_angle_alpha "<< FormatFloat(this->GetLatticePar(3)*RAD2DEG,8,5) << endl 
      << "_cell_angle_beta " << FormatFloat(this->GetLatticePar(4)*RAD2DEG,8,5) << endl 
      << "_cell_angle_gamma "<< FormatFloat(this->GetLatticePar(5)*RAD2DEG,8,5) << endl ;
   os <<endl;
   
   this->GetScatteringComponentList();
   
   os << "loop_" << endl
      << "    _atom_site_label" <<endl
      << "    _atom_site_fract_x"<<endl
      << "    _atom_site_fract_y" <<endl
      << "    _atom_site_fract_z" <<endl
      << "    _atom_site_U_iso_or_equiv" <<endl
      << "    _atom_site_thermal_displace_type" <<endl
      << "    _atom_site_occupancy" <<endl
      << "    _atom_site_type_symbol" <<endl;
   
   long k=0;
   for(int i=0;i<mScattererRegistry.GetNb();i++) 
   {
      //mpScatterrer[i]->Print();
      const ScatteringComponentList list=this->GetScatt(i).GetScatteringComponentList();
      for(int j=0;j<list.GetNbComponent();j++)
      {
         os   << "    "
              << FormatString(this->GetScatt(i).GetComponentName(j),16)
              << FormatFloat(list(j).mX,7,4) 
              << FormatFloat(list(j).mY,7,4) 
              << FormatFloat(list(j).mZ,7,4) 
              << FormatFloat(list(j).mpScattPow->GetBiso()/8./M_PI/M_PI)
              << " 1 "
              << FormatFloat(list(j).mOccupancy,6,4)
              << FormatString(list(j).mpScattPow->GetName(),16)
              << endl;
         k++;
      }
      
   }
   os <<endl;
   k=0;
   if(1==mUseDynPopCorr.GetChoice())
   {
      os << ";  Dynamical occupancy corrections found by ObjCryst++:"<<endl
         << "   values below 1. (100%) indicate a correction,"<<endl
         << "   which means either that the atom is on a special position,"<<endl
         << "   or that it is overlapping with another identical atom."<<endl;
      for(int i=0;i<mScattererRegistry.GetNb();i++) 
      {
         //mpScatterrer[i]->Print();
         const ScatteringComponentList list=this->GetScatt(i).GetScatteringComponentList();
         for(int j=0;j<list.GetNbComponent();j++)
         {
            os   << "    "
                 << FormatString(this->GetScatt(i).GetComponentName(j),16)
                 << " : " << FormatFloat(mScattCompList(k).mDynPopCorr,6,4)
                 << endl;
            k++;
         }
      }
      os << ";"<<endl;
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
   DistTableInternalPosition(const long atomIndex, const int sym,
                             const long x,const long y,const long z):
   mAtomIndex(atomIndex),mSymmetryIndex(sym),mXL(x),mYL(y),mZL(z)
   {}
   long mAtomIndex;
   int mSymmetryIndex;
   REAL mX,mY,mZ;
   long mXL,mYL,mZL;
};

void Crystal::CalcDistTable(const bool fast, const REAL asymUnitMargin) const
{
   this->GetScatteringComponentList();
   
   if(  (mDistTableClock>mClockScattCompList)
      &&(mDistTableClock>this->GetClockMetricMatrix())) return;
   VFN_DEBUG_ENTRY("Crystal::CalcDistTable()",4)
      
   const long nbComponent=mScattCompList.GetNbComponent();
   const int nbSymmetrics=this->GetSpaceGroup().GetNbSymmetrics();
   
   mvDistTableSq.resize(nbComponent);
   {
      std::vector<NeighbourHood>::iterator pos;
      for(pos=mvDistTableSq.begin();pos<mvDistTableSq.end();pos++)
         pos->mvNeighbour.clear();
   }
   VFN_DEBUG_MESSAGE("Crystal::CalcDistTable():1",3)
   
   // Get limits of the (pseudo) asymmetric unit
      // strict limits
      const REAL xMax0=this->GetSpaceGroup().GetAsymUnit().Xmax();
      const REAL yMax0=this->GetSpaceGroup().GetAsymUnit().Ymax();
      const REAL zMax0=this->GetSpaceGroup().GetAsymUnit().Zmax();
   
      // limits with a margin, within [0;1[
      const REAL xMax=this->GetSpaceGroup().GetAsymUnit().Xmax()+asymUnitMargin/GetLatticePar(0);
      const REAL yMax=this->GetSpaceGroup().GetAsymUnit().Ymax()+asymUnitMargin/GetLatticePar(1);
      const REAL zMax=this->GetSpaceGroup().GetAsymUnit().Zmax()+asymUnitMargin/GetLatticePar(2);
      const REAL xMin=1.-asymUnitMargin/GetLatticePar(0);
      const REAL yMin=1.-asymUnitMargin/GetLatticePar(1);
      const REAL zMin=1.-asymUnitMargin/GetLatticePar(2);
      
      bool useAsymUnit=false;
      if((xMax0*yMax0*zMax0)<0.6) useAsymUnit=true;
   
   
   // List of all positions within or near the asymmetric unit
   std::vector<DistTableInternalPosition> vPos;
   // index of unique atoms in vPos, which are strictly in the asymmetric unit
   std::vector<unsigned long> vUniqueIndex(nbComponent);
   
   if(true==fast)//
   {
      VFN_DEBUG_MESSAGE("Crystal::CalcDistTable(fast):2",3)
      TAU_PROFILE("Crystal::CalcDistTable(fast=true)","Matrix (string&)",TAU_DEFAULT);
      TAU_PROFILE_TIMER(timer1,"DiffractionData::CalcDistTable1","", TAU_FIELD);
      TAU_PROFILE_TIMER(timer2,"DiffractionData::CalcDistTable2","", TAU_FIELD);
      
      TAU_PROFILE_START(timer1);
      #define FRAC2LONG 0x4000
      #define FRAC2LONGMASK 0x3FFF
      #define HALF_FRAC2LONG 0x2000
      #define HALF_FRAC2LONGMASK 0x1FFF
      
      // Get limits of the (pseudo) asymmetric unit in long
         // strict limits
         const long xMax0l=(long)(xMax0*(REAL)FRAC2LONG);
         const long yMax0l=(long)(yMax0*(REAL)FRAC2LONG);
         const long zMax0l=(long)(zMax0*(REAL)FRAC2LONG);

         // limits with a margin, within [0;FRAC2LONG[
         const long xMaxl=(long)(xMax*(REAL)FRAC2LONG);
         const long yMaxl=(long)(yMax*(REAL)FRAC2LONG);
         const long zMaxl=(long)(zMax*(REAL)FRAC2LONG);
         const long xMinl=(long)(xMin*(REAL)FRAC2LONG);
         const long yMinl=(long)(yMin*(REAL)FRAC2LONG);
         const long zMinl=(long)(zMin*(REAL)FRAC2LONG);
      
      long xl,yl,zl;
      
      CrystMatrix_REAL symmetricsCoords;
      
      // Get the list of all atoms within or near the asymmetric unit
      for(long i=0;i<nbComponent;i++)
      {
         VFN_DEBUG_MESSAGE("Crystal::CalcDistTable(fast):3:component "<<i,0)
         symmetricsCoords=this->GetSpaceGroup().GetAllSymmetrics(mScattCompList(i).mX,
                                                                 mScattCompList(i).mY,
                                                                 mScattCompList(i).mZ);
         mvDistTableSq[i].mIndex=i;//USELESS ?
         for(int j=0;j<nbSymmetrics;j++)
         {
            // Convert to long [0;1[ -> [0 ; FRAC2LONG[ with some bit twidlling
            xl=(long)(symmetricsCoords(j,0)*(REAL)FRAC2LONG);
            yl=(long)(symmetricsCoords(j,1)*(REAL)FRAC2LONG);
            zl=(long)(symmetricsCoords(j,2)*(REAL)FRAC2LONG);
            xl= xl & FRAC2LONGMASK; //xl %= FRAC2LONG;if(xl<0) xl += FRAC2LONG;
            yl= yl & FRAC2LONGMASK; //yl %= FRAC2LONG;if(yl<0) yl += FRAC2LONG;
            zl= zl & FRAC2LONGMASK; //zl %= FRAC2LONG;if(zl<0) zl += FRAC2LONG;
            
            if(useAsymUnit)
            {
               if( (zl>zMinl) || (zl<zMaxl)) 
                  if( (xl>xMinl) || (xl<xMaxl)) 
                     if( (yl>yMinl) || (yl<yMaxl))
                     {
                        vPos.push_back(DistTableInternalPosition(i,j,xl,yl,zl));
                        if((xl<xMax0l)&&(yl<yMax0l)&&(zl<zMax0l))
                        {
                           vUniqueIndex[i]=vPos.size()-1;
                           mvDistTableSq[i].mUniquePosSymmetryIndex=j;
                        }
                     }
            }
            else
            {
               vPos.push_back(DistTableInternalPosition(i,j,xl,yl,zl));
               vUniqueIndex[i]=vPos.size()-1;
               mvDistTableSq[i].mUniquePosSymmetryIndex=j;
            }
         }
      }
      TAU_PROFILE_STOP(timer1);
      TAU_PROFILE_START(timer2);
      // Compute interatomic vectors & distance 
      // between (i) unique atoms and (ii) all remaining atoms
      
      const CrystMatrix_REAL* pOrthMatrix=&(this->GetOrthMatrix());
      const REAL M00=(*pOrthMatrix)(0,0)/(REAL)FRAC2LONG;
      const REAL M01=(*pOrthMatrix)(0,1)/(REAL)FRAC2LONG;
      const REAL M02=(*pOrthMatrix)(0,2)/(REAL)FRAC2LONG;
      //const REAL M10=(*pOrthMatrix)(1,0)/(REAL)FRAC2LONG;
      const REAL M11=(*pOrthMatrix)(1,1)/(REAL)FRAC2LONG;
      const REAL M12=(*pOrthMatrix)(1,2)/(REAL)FRAC2LONG;
      //const REAL M20=(*pOrthMatrix)(2,0)/(REAL)FRAC2LONG;
      //const REAL M21=(*pOrthMatrix)(2,1)/(REAL)FRAC2LONG;
      const REAL M22=(*pOrthMatrix)(2,2)/(REAL)FRAC2LONG;
      
      for(long i=0;i<nbComponent;i++)
      {
         VFN_DEBUG_MESSAGE("Crystal::CalcDistTable(fast):4:component "<<i,0)
         #if 0
         cout<<endl<<"Unique pos:"<<vUniqueIndex[i]<<":"
             <<vPos[vUniqueIndex[i]].mAtomIndex<<":"
             <<mScattCompList(vPos[vUniqueIndex[i]].mAtomIndex).mpScattPow->GetName()<<":"
             <<vPos[vUniqueIndex[i]].mSymmetryIndex<<":"
             <<vPos[vUniqueIndex[i]].mXL/(REAL)FRAC2LONG<<","
             <<vPos[vUniqueIndex[i]].mYL/(REAL)FRAC2LONG<<","
             <<vPos[vUniqueIndex[i]].mZL/(REAL)FRAC2LONG<<endl;
         #endif
         std::vector<Crystal::Neighbour> * const vnb=&(mvDistTableSq[i].mvNeighbour);
         for(unsigned long j=0;j<vPos.size();j++)
         {
            if(vUniqueIndex[i]==j) continue;
            xl=vPos[j].mXL - vPos[ vUniqueIndex[i] ].mXL;
            yl=vPos[j].mYL - vPos[ vUniqueIndex[i] ].mYL;
            zl=vPos[j].mZL - vPos[ vUniqueIndex[i] ].mZL;
            
            xl= xl & FRAC2LONGMASK;  //xl %= FRAC2LONG;if(xl<0) xl += FRAC2LONG;
            yl= yl & FRAC2LONGMASK;  //yl %= FRAC2LONG;if(yl<0) yl += FRAC2LONG;
            zl= zl & FRAC2LONGMASK;  //zl %= FRAC2LONG;if(zl<0) zl += FRAC2LONG;
            
            if(xl & HALF_FRAC2LONG) xl= (~xl) & HALF_FRAC2LONGMASK;//if(xl>HALF_FRAC2LONG) xl = FRAC2LONG-xl;
            if(yl & HALF_FRAC2LONG) yl= (~yl) & HALF_FRAC2LONGMASK;//if(yl>HALF_FRAC2LONG) yl = FRAC2LONG-yl;
            if(zl & HALF_FRAC2LONG) zl= (~zl) & HALF_FRAC2LONGMASK;//if(zl>HALF_FRAC2LONG) zl = FRAC2LONG-zl;
            
            const REAL x=M00 * (REAL)xl + M01 * (REAL)yl + M02 * (REAL)zl;
            const REAL y=                 M11 * (REAL)yl + M12 * (REAL)zl;
            const REAL z=                                  M22 * (REAL)zl;
            Neighbour neigh(vPos[j].mAtomIndex,vPos[j].mSymmetryIndex,x*x+y*y+z*z);
            vnb->push_back(neigh);
            #if 0
            cout<<vPos[j].mAtomIndex<<":"
                <<mScattCompList(vPos[j].mAtomIndex).mpScattPow->GetName()<<":"
                <<vPos[j].mSymmetryIndex<<":"
                <<vPos[j].mXL/(REAL)FRAC2LONG<<","
                <<vPos[j].mYL/(REAL)FRAC2LONG<<","
                <<vPos[j].mZL/(REAL)FRAC2LONG<<" vector="
                <<x<<","<<y<<","<<z<<endl;
            #endif
         }
      }
      TAU_PROFILE_STOP(timer2);
   }
   else
   {
      VFN_DEBUG_MESSAGE("Crystal::CalcDistTable(fast):2",3)
      TAU_PROFILE("Crystal::CalcDistTable(fast=false)","Matrix (string&)",TAU_DEFAULT);
      TAU_PROFILE_TIMER(timer1,"DiffractionData::CalcDistTable1","", TAU_FIELD);
      TAU_PROFILE_TIMER(timer2,"DiffractionData::CalcDistTable2","", TAU_FIELD);
      
      TAU_PROFILE_START(timer1);
      CrystMatrix_REAL symmetricsCoords;
      
      REAL x,y,z;
      double junk;
      // Get the list of all atoms within or near the asymmetric unit
      for(long i=0;i<nbComponent;i++)
      {
         VFN_DEBUG_MESSAGE("Crystal::CalcDistTable(fast):3:component "<<i,3)
         symmetricsCoords=this->GetSpaceGroup().GetAllSymmetrics(mScattCompList(i).mX,
                                                                 mScattCompList(i).mY,
                                                                 mScattCompList(i).mZ);
         mvDistTableSq[i].mIndex=i;//USELESS ?
         for(int j=0;j<nbSymmetrics;j++)
         {
            x=modf(symmetricsCoords(j,0),&junk);
            y=modf(symmetricsCoords(j,1),&junk);
            z=modf(symmetricsCoords(j,2),&junk);
            if(x<0) x +=1.;
            if(y<0) y +=1.;
            if(z<0) z +=1.;
            
            if(useAsymUnit)
            {
               if( (z>zMin) || (z<zMax)) 
                  if( (x>xMin) || (x<xMax)) 
                     if( (y>yMin) || (y<yMax))
                     {
                        vPos.push_back(DistTableInternalPosition(i,j,x,y,z));
                        if((x<xMax0)&&(y<yMax0)&&(z<zMax0))
                        {
                           vUniqueIndex[i]=vPos.size()-1;
                           mvDistTableSq[i].mUniquePosSymmetryIndex=j;
                        }
                     }
            }
            else
            {
               vPos.push_back(DistTableInternalPosition(i,j,x,y,z));
               vUniqueIndex[i]=vPos.size()-1;
               mvDistTableSq[i].mUniquePosSymmetryIndex=j;
            }
         }
      }
      TAU_PROFILE_STOP(timer1);
      TAU_PROFILE_START(timer2);
      // Compute interatomic vectors & distance 
      // between (i) unique atoms and (ii) all remaining atoms
      for(long i=0;i<nbComponent;i++)
      {
         VFN_DEBUG_MESSAGE("Crystal::CalcDistTable(fast):4:component "<<i,3)
         #if 0
         cout<<endl<<"Unique pos:"
             <<vPos[vUniqueIndex[i]].mAtomIndex<<":"
             <<mScattCompList(vPos[vUniqueIndex[i]].mAtomIndex).mpScattPow->GetName()<<":"
             <<vPos[vUniqueIndex[i]].mSymmetryIndex<<":"
             <<vPos[vUniqueIndex[i]].mX<<","
             <<vPos[vUniqueIndex[i]].mY<<","
             <<vPos[vUniqueIndex[i]].mZ<<endl;
         #endif
         std::vector<Crystal::Neighbour> * const vnb=&(mvDistTableSq[i].mvNeighbour);
         VFN_DEBUG_MESSAGE("Crystal::CalcDistTable(fast):4:vector "<<vnb->size(),3)
         for(unsigned long j=0;j<vPos.size();j++)
         {
            if(vUniqueIndex[i]==j) continue;
            x=modf(vPos[j].mX - vPos[ vUniqueIndex[i] ].mX,&junk);
            y=modf(vPos[j].mY - vPos[ vUniqueIndex[i] ].mY,&junk);
            z=modf(vPos[j].mZ - vPos[ vUniqueIndex[i] ].mZ,&junk);
            if(x<0) x+=1 ;
            if(y<0) y+=1 ;
            if(z<0) z+=1 ;
            if(x>0.5) x=1-x ;
            if(y>0.5) y=1-y ;
            if(z>0.5) z=1-z ;
            this->FractionalToOrthonormalCoords(x,y,z);
            Neighbour neigh(vPos[j].mAtomIndex,vPos[j].mSymmetryIndex,x*x+y*y+z*z);
            vnb->push_back(neigh);
            #if 0
            cout<<vPos[j].mAtomIndex<<":"
                <<mScattCompList(vPos[j].mAtomIndex).mpScattPow->GetName()<<":"
                <<vPos[j].mSymmetryIndex<<":"
                <<vPos[j].mX<<","
                <<vPos[j].mY<<","
                <<vPos[j].mZ<<" vector="
                <<x<<","<<y<<","<<z<<endl;
            #endif
         }
      }
      TAU_PROFILE_STOP(timer2);
   }
   mDistTableClock.Click();
   VFN_DEBUG_EXIT("Crystal::CalcDistTable()",4)
}

#ifdef __WX__CRYST__
WXCrystObjBasic* Crystal::WXCreate(wxWindow* parent)
{
   //:TODO: Check mpWXCrystObj==0
   mpWXCrystObj=new WXCrystal(parent,this);
   return mpWXCrystObj;
}
#endif

}//namespace
