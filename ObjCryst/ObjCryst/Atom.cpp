/* 
* ObjCryst++ : a Crystallographic computing library in C++
*			http://objcryst.sourceforge.net
*			http://www.ccp14.ac.uk/ccp/web-mirrors/objcryst/
*
*  (c) 2000-2001 Vincent FAVRE-NICOLIN vincefn@users.sourceforge.net
*
*/
/*   Atom.cpp
*  source file for the Atom scatterer
*
*/
#include <cmath>
#include <typeinfo>
#include <stdio.h> //for sprintf()

#include "ObjCryst/Atom.h"

#include "Quirks/VFNStreamFormat.h" //simple formatting of integers, REALs..
#include "Quirks/VFNDebug.h"

#ifdef OBJCRYST_GL
extern "C"
{
#include <GL/glut.h>
}
#endif

#ifdef __WX__CRYST__
   #include "wxCryst/wxAtom.h"
#endif

#include <fstream>
#include <iomanip>

namespace ObjCryst
{

////////////////////////////////////////////////////////////////////////
//
// Using external functions from 'atominfo' package
//
//-> Coefficients for the analytical approximation of scattering factors.
//AtomInfo (c) 1994-96 Ralf W. Grosse-Kunstleve 
//
////////////////////////////////////////////////////////////////////////

extern "C"
{
#include "atominfo/atominfo.h"
}

////////////////////////////////////////////////////////////////////////
//
//    ATOM : the basic atom, within the crystal
//
//includes thermic factor (betas) and x,y,z,population
//
////////////////////////////////////////////////////////////////////////
Atom::Atom()
:mScattCompList(1),mpScattPowAtom(0)
{
   VFN_DEBUG_MESSAGE("Atom::Atom()",5)
   this->InitRefParList();
   mScattCompList(0).mpScattPow=mpScattPowAtom;
}

Atom::Atom( const REAL x, const REAL y, const REAL z,
            const string &name,const ScatteringPowerAtom *pow)
:mScattCompList(1),mpScattPowAtom(pow)
{
   VFN_DEBUG_MESSAGE("Atom::Atom(x,y,z,name,ScatteringPower):"<<name,5)
   this->Init(x,y,z,name,pow,1);
}

Atom::Atom( const REAL x, const REAL y, const REAL z, const string &name,
            const ScatteringPowerAtom *pow,const REAL popu)
:mScattCompList(1),mpScattPowAtom(pow)
{
   VFN_DEBUG_MESSAGE("Atom::Atom(x,y,z,P,B,name,ScatteringPower):"<<name,5)
   this->Init(x,y,z,name,mpScattPowAtom,popu);
}

Atom::Atom(const Atom &old)
:Scatterer(old),mScattCompList(1),mpScattPowAtom(old.mpScattPowAtom)
{
   VFN_DEBUG_MESSAGE("Atom::Atom(&old):/Name="<<old.mName,5)
   //:TODO: Check if this is enough (copy constructor for Scatterer ??)
   this->Init(old.mXYZ(0),old.mXYZ(1),old.mXYZ(2),
              old.mName,mpScattPowAtom,
              old.mOccupancy);
}

Atom* Atom::CreateCopy() const
{
   VFN_DEBUG_MESSAGE("Atom::CreateCopy():/Name="<<mName,5)
   return new Atom(*this);
}


Atom::~Atom()
{
   VFN_DEBUG_MESSAGE("Atom::~Atom():("<<mName<<")",5)
}

const string& Atom::GetClassName()const
{
	const static string className="Atom";
	return className;
}

void Atom::operator=(const Atom &rhs)
{
   VFN_DEBUG_MESSAGE("Atom::operator=():/Name="<<rhs.mName,5)
   //:TODO: Check if this is enough (copy constructor for Scatterer ??)
   mScattCompList.Reset();
   ++mScattCompList;
   
   this->Init(rhs.mXYZ(0),rhs.mXYZ(1),rhs.mXYZ(2),
              rhs.mName,rhs.mpScattPowAtom,
              rhs.mOccupancy);
}

void Atom::Init(const REAL x, const REAL y, const REAL z,
            const string &name, const ScatteringPowerAtom *pow,
            const REAL popu)
{
   VFN_DEBUG_MESSAGE("Atom::Init():"<<name,3)
   mName=name;
   mpScattPowAtom=pow;
   mScattCompList(0).mpScattPow=mpScattPowAtom;
   mXYZ(0)=x;
   mXYZ(1)=y;
   mXYZ(2)=z;
   if(0==mpScattPowAtom)
   {//Dummy atom
      VFN_DEBUG_MESSAGE("Atom::Init()Dummy Atom:/Name="<<this->GetName(),5)
      mOccupancy=0;
      return;
   }
   mpScattPowAtom->RegisterClient(*this);
   mOccupancy=popu;
   
   this->InitRefParList();
   
   mClockScatterer.Click();
   VFN_DEBUG_MESSAGE("Atom::Init():End.",5)
}

int Atom::GetNbComponent() const {return 1;}
const ScatteringComponentList& Atom::GetScatteringComponentList() const
{
   mScattCompList(0).mX=mXYZ(0);
   mScattCompList(0).mY=mXYZ(1);
   mScattCompList(0).mZ=mXYZ(2);
   mScattCompList(0).mOccupancy=mOccupancy;
   mClockScattCompList.Click();
   return mScattCompList;
}

string Atom::GetComponentName(const int i) const{ return this->GetName();}

void Atom::Print() const
{
   VFN_DEBUG_MESSAGE("Atom::Print()",1)
   cout << "Atom (" 
        << FormatString(mpScattPowAtom->GetSymbol(),4) << ") :"
        << FormatString(this->GetName(),16)  << " at : " 
        << FormatFloat(this->GetX()) 
        << FormatFloat(this->GetY()) 
        << FormatFloat(this->GetZ());
   if(this->IsDummy()) cout << " DUMMY! ";
   else
   {
      cout << ", Biso=" 
           << FormatFloat(mpScattPowAtom->GetBiso())
           << ", Popu=" << FormatFloat(this->GetOccupancy());
   }
   cout << endl;
}

REAL Atom::GetMass() const
{
   if(true==this->IsDummy()) return 0;
   return mpScattPowAtom->GetBiso();
}

REAL Atom::GetRadius() const
{
   if(this->IsDummy()) return 0.5;
   return mpScattPowAtom->GetRadius();
}

ostream& Atom::POVRayDescription(ostream &os,
                                 const bool onlyIndependentAtoms)const
{
   if(this->IsDummy()) return os;
   if(true==onlyIndependentAtoms)
   {
      REAL x,y,z;
      x=mXYZ(0);
      y=mXYZ(1);
      z=mXYZ(2);
      x = fmod((REAL)x,(int)1); if(x<0) x+=1.;
      y = fmod((REAL)y,(int)1); if(y<0) y+=1.;
      z = fmod((REAL)z,(int)1); if(z<0) z+=1.;
      this->GetCrystal().FractionalToOrthonormalCoords(x,y,z);
      os << "// Description of Atom :" << this->GetName() <<endl;
      os << "   #declare colour_"<< this ->GetName() <<"="<< this ->mColourName<<";"<< endl;
      os << "   sphere " << endl;
      os << "   { <" << x << ","<< y << ","<< z << ">," << this->GetRadius()/3 <<endl;
      os << "        finish { ambient 0.2 diffuse 0.8 phong 1}" <<endl;
      os << "        pigment { colour colour_"<< this->GetName() <<" }" << endl;
      os << "   }" <<endl;
   }
   else
   {
      REAL x0,y0,z0;
      x0=mXYZ(0);
      y0=mXYZ(1);
      z0=mXYZ(2);
      CrystMatrix_REAL xyzCoords ;
      xyzCoords=this->GetCrystal().GetSpaceGroup().GetAllSymmetrics(x0,y0,z0,false,false,true);
      int nbSymmetrics=xyzCoords.rows();
      os << "// Description of Atom :" << this->GetName();
      os << "(" << nbSymmetrics << " symmetric atoms)" << endl;
      os << "   #declare colour_"<< this ->GetName() <<"="<< this ->mColourName<<";"<< endl;
      int symNum=0;
      for(int i=0;i<nbSymmetrics;i++)
      {
         x0=xyzCoords(i,0);
         y0=xyzCoords(i,1);
         z0=xyzCoords(i,2);
         x0 = fmod((REAL) x0,(int)1); if(x0<0) x0+=1.;
         y0 = fmod((REAL) y0,(int)1); if(y0<0) y0+=1.;
         z0 = fmod((REAL) z0,(int)1); if(z0<0) z0+=1.;
         //Generate also translated atoms near the unit cell
         const REAL limit =0.1;
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
         for(int j=0;j<translate.rows();j++)
         {
            REAL x=x0+translate(j,0);
            REAL y=y0+translate(j,1);
            REAL z=z0+translate(j,2);
            if(   (x>(-limit)) && (x<(1+limit))
                &&(y>(-limit)) && (y<(1+limit))
                &&(z>(-limit)) && (z<(1+limit)))
            {
               this->GetCrystal().FractionalToOrthonormalCoords(x,y,z);
               os << "  // Symmetric #" << symNum++ <<endl;
               os << "   sphere " << endl;
               os << "   { <" << x << ","<< y << ","<< z << ">," << this->GetRadius()/3 <<endl;
               os << "        finish { ambient 0.2 diffuse 0.8 phong 1}" <<endl;
               //os << "        pigment { colour red 1 green 0 blue 1 }" << endl;
               os << "        pigment { colour colour_"<< this->GetName() <<" }" << endl;
               os << "   }" <<endl;
            }
         }
      }
   }
   return os;
}

void Atom::GLInitDisplayList(const bool onlyIndependentAtoms,
                             const REAL xMin,const REAL xMax,
                             const REAL yMin,const REAL yMax,
                             const REAL zMin,const REAL zMax)const
{
   #ifdef OBJCRYST_GL
   VFN_DEBUG_MESSAGE("Atom::GLInitDisplayList():"<<this->GetName(),5)
   glMaterialfv (GL_FRONT, GL_AMBIENT_AND_DIFFUSE, mpScattPowAtom->GetColourRGB());
   
   if(this->IsDummy()) return ;
   GLUquadricObj* pQuadric = gluNewQuadric();
   if(true==onlyIndependentAtoms)
   {
      REAL x,y,z;
      x=mXYZ(0);
      y=mXYZ(1);
      z=mXYZ(2);
      x = fmod((REAL)x,(int)1); if(x<0) x+=1.;
      y = fmod((REAL)y,(int)1); if(y<0) y+=1.;
      z = fmod((REAL)z,(int)1); if(z<0) z+=1.;
      this->GetCrystal().FractionalToOrthonormalCoords(x,y,z);
      glPushMatrix();
         glTranslatef(x, y, z);
         gluSphere(pQuadric,this->GetRadius()/3,10,10);
      glPopMatrix();
   }
   else
   {
      REAL x0,y0,z0;
      x0=mXYZ(0);
      y0=mXYZ(1);
      z0=mXYZ(2);
      CrystMatrix_REAL xyzCoords ;
      xyzCoords=this->GetCrystal().GetSpaceGroup().GetAllSymmetrics(x0,y0,z0,false,false,true);
      int nbSymmetrics=xyzCoords.rows();
      for(int i=0;i<nbSymmetrics;i++)
      {
         x0=xyzCoords(i,0);
         y0=xyzCoords(i,1);
         z0=xyzCoords(i,2);
         x0 = fmod((REAL) x0,(int)1); if(x0<0) x0+=1.;
         y0 = fmod((REAL) y0,(int)1); if(y0<0) y0+=1.;
         z0 = fmod((REAL) z0,(int)1); if(z0<0) z0+=1.;
         //Generate also translated atoms near the unit cell
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
         for(int j=0;j<translate.rows();j++)
         {
            REAL x=x0+translate(j,0);
            REAL y=y0+translate(j,1);
            REAL z=z0+translate(j,2);
            if(   (x>xMin) && (x<xMax)
                &&(y>yMin) && (y<yMax)
                &&(z>zMin) && (z<zMax))
            {
               this->GetCrystal().FractionalToOrthonormalCoords(x,y,z);
               glPushMatrix();
                  glTranslatef(x, y, z);
                  gluSphere(pQuadric,this->GetRadius()/3.,10,10);
               glPopMatrix();
            }
         }
      }
   }
   gluDeleteQuadric(pQuadric);
   VFN_DEBUG_MESSAGE("Atom::GLInitDisplayList():End",5)
   #endif
}

bool Atom::IsDummy()const { if(0==mScattCompList(0).mpScattPow) return true; return false;}

const ScatteringPowerAtom& Atom::GetScatteringPower()const
{ return *mpScattPowAtom;}

void Atom::GetGeneGroup(const RefinableObj &obj,
										  CrystVector_uint & groupIndex,
										  unsigned int &first) const
{
	//One group for all translation parameters
	unsigned int posIndex=0;
   VFN_DEBUG_MESSAGE("Atom::GetGeneGroup()",4)
	for(long i=0;i<obj.GetNbPar();i++)
		for(long j=0;j<this->GetNbPar();j++)
			if(&(obj.GetPar(i)) == &(this->GetPar(j)))
			{
				if(this->GetPar(j).GetType()->IsDescendantFromOrSameAs(gpRefParTypeScattTransl))
				{
					if(posIndex==0) posIndex=first++;
					groupIndex(i)=posIndex;
				}
				else groupIndex(i)= first++;
			}
}

void Atom::InitRefParList()
{
   mClockScatterer.Click();
   VFN_DEBUG_MESSAGE("Atom::InitRefParList()",5)
   this->ResetParList();
   //:TODO: Add thermic factors
   {
      RefinablePar tmp("x",&mXYZ(0),0,1.,gpRefParTypeScattTranslX,
                        REFPAR_DERIV_STEP_ABSOLUTE,false,false,true,true,1.,1.);
      tmp.AssignClock(mClockScatterer);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("y",&mXYZ(1),0,1.,gpRefParTypeScattTranslY,
                        REFPAR_DERIV_STEP_ABSOLUTE,false,false,true,true,1.,1.);
      tmp.AssignClock(mClockScatterer);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("z",&mXYZ(2),0,1.,gpRefParTypeScattTranslZ,
                        REFPAR_DERIV_STEP_ABSOLUTE,false,false,true,true,1.,1.);
      tmp.AssignClock(mClockScatterer);
      this->AddPar(tmp);
   }
   if(false==this->IsDummy())
   {
      {
         RefinablePar tmp(this->GetName()+(string)"occup",&mOccupancy,0.01,1.,
                           gpRefParTypeScattOccup,REFPAR_DERIV_STEP_ABSOLUTE,true,true);
         tmp.AssignClock(mClockScatterer);
         this->AddPar(tmp);
      }
   }
}
#ifdef __WX__CRYST__
WXCrystObjBasic* Atom::WXCreate(wxWindow* parent)
{
   //:TODO: Check mpWXCrystObj==0
   mpWXCrystObj=new WXAtom(parent,this);
   return mpWXCrystObj;
}
#endif

}//namespace
