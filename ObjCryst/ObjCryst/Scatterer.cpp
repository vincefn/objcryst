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
//
// source file for LibCryst++ Scatterer classes

#include <cmath>
#include <typeinfo>

#include <stdio.h> //for sprintf()

#include "ObjCryst/Crystal.h"
#include "ObjCryst/Scatterer.h"

#include "Quirks/VFNStreamFormat.h" //simple formatting of integers, REALs..

#include "Quirks/VFNDebug.h"
#ifdef OBJCRYST_GL
extern "C"
{
#include <GL/glut.h>
}
#endif

#ifdef __WX__CRYST__
   #include "wxCryst/wxScatterer.h"
#endif

#include "ObjCryst/Colours.h"

#include <iomanip>
#include <fstream>

namespace ObjCryst
{

const RefParType *gpRefParTypeScatt=
   new RefParType(gpRefParTypeObjCryst,"Scatterer");
const RefParType *gpRefParTypeScattTransl=
   new RefParType(gpRefParTypeScatt,"Translation");
const RefParType *gpRefParTypeScattTranslX=
   new RefParType(gpRefParTypeScattTransl,"Translation along X");
const RefParType *gpRefParTypeScattTranslY=
   new RefParType(gpRefParTypeScattTransl,"Translation along Y");
const RefParType *gpRefParTypeScattTranslZ=
   new RefParType(gpRefParTypeScattTransl,"Translation along Z");
const RefParType *gpRefParTypeScattOrient=
   new RefParType(gpRefParTypeScatt,"Orientation");
const RefParType *gpRefParTypeScattConform=
   new RefParType(gpRefParTypeScatt,"Conformation");
const RefParType *gpRefParTypeScattConformBondLength=
   new RefParType(gpRefParTypeScattConform,"BondLengths");
const RefParType *gpRefParTypeScattConformBondAngle=
   new RefParType(gpRefParTypeScattConform,"Bond Angles");
const RefParType *gpRefParTypeScattConformDihedAngle=
   new RefParType(gpRefParTypeScattConform,"Dihedral Angles ");
const RefParType *gpRefParTypeScattConformX=
   new RefParType(gpRefParTypeScatt,"Orth. X coordinates");
const RefParType *gpRefParTypeScattConformY=
   new RefParType(gpRefParTypeScatt,"Orth. Y coordinates");
const RefParType *gpRefParTypeScattConformZ=
   new RefParType(gpRefParTypeScatt,"Orth. Z coordinates");
const RefParType *gpRefParTypeScattOccup=
   new RefParType(gpRefParTypeScatt,"Occupancy");

////////////////////////////////////////////////////////////////////////
//
//      SCATTERER : generic type of scatterer
//
//
////////////////////////////////////////////////////////////////////////
ObjRegistry<Scatterer> gScattererRegistry("Global Scatterer Registry");

Scatterer::Scatterer():mXYZ(3),mOccupancy(1.0),mColourName("White"),mpCryst(0)
{
   VFN_DEBUG_MESSAGE("Scatterer::Scatterer()",5)
   mXYZ=0;
   this->InitRGBColour();
   gScattererRegistry.Register(*this);
}

Scatterer::Scatterer(const Scatterer &old):
//RefinableObj(),
mXYZ(old.mXYZ),
mOccupancy(old.mOccupancy),
mColourName(old.mColourName),mpCryst(old.mpCryst)
{
   VFN_DEBUG_MESSAGE("Scatterer::Scatterer(&old)",5)
	mName=old.GetName();
   this->InitRGBColour();
   gScattererRegistry.Register(*this);
}

Scatterer::~Scatterer()
{
   VFN_DEBUG_MESSAGE("Scatterer::~Scatterer():("<<mName<<")",5)
   gScattererRegistry.DeRegister(*this);
}

const string& Scatterer::GetClassName() const
{
	const static string className="Scatterer";
	return className;
}

REAL Scatterer::GetX()    const {return mXYZ(0);}
REAL Scatterer::GetY()    const {return mXYZ(1);}
REAL Scatterer::GetZ()    const {return mXYZ(2);}
REAL Scatterer::GetOccupancy() const {return mOccupancy;}


void Scatterer::SetX(const REAL x) { this->GetPar(mXYZ.data()).MutateTo(x);}
void Scatterer::SetY(const REAL y) { this->GetPar(mXYZ.data()+1).MutateTo(y);}
void Scatterer::SetZ(const REAL z) { this->GetPar(mXYZ.data()+2).MutateTo(z);}
void Scatterer::SetOccupancy(const REAL occupancy)
{
   mClockScatterer.Click();
   mOccupancy=occupancy;
}

Scatterer::operator string() const {return this->GetName();}
/*
void Scatterer::Print() const
{
   cout << "Scatterer : " << FormatString(this->GetName()) << " at :" ;
   cout <<  FormatFloat(this->X());
   cout <<  FormatFloat(this->Y());
   cout <<  FormatFloat(this->Z());
   cout << endl;
}
*/

const string& Scatterer::GetColour()const{ return mColourName;}
const float* Scatterer::GetColourRGB()const{ return mColourRGB;}

const RefinableObjClock& Scatterer::GetClockScatterer()const {return mClockScatterer;}

void Scatterer::SetCrystal(Crystal &cryst){mpCryst=&cryst;}
const Crystal& Scatterer::GetCrystal()const{return *mpCryst;}

void Scatterer::InitRGBColour()
{
   //mClockScatterer.Click();
   VFN_DEBUG_MESSAGE("Scatterer::InitRGBColour()",2)
   for(long i=0;;)
   {
      if(gPOVRayColours[i].mName==mColourName)
      {
         mColourRGB[0]=gPOVRayColours[i].mRGB[0];
         mColourRGB[1]=gPOVRayColours[i].mRGB[1];
         mColourRGB[2]=gPOVRayColours[i].mRGB[2];
         break;
      }
      i++;
      if(gPOVRayColours[i].mName=="")
      {//could not find colour !
         cout << "Could not find colour:"<<mColourName<<" for scaterrer "<<mName<<endl;
         mColourRGB[0]=1;
         mColourRGB[1]=1;
         mColourRGB[2]=1;
         break;
      }
   }
   VFN_DEBUG_MESSAGE("->RGBColour:"<<mColourName<<mColourRGB[0]<<" "<<mColourRGB[1]<<" "<<mColourRGB[2],2)
}
#ifdef __WX__CRYST__
WXCrystObjBasic* Scatterer::WXCreate(wxWindow* parent)
{
   //:TODO: Check mpWXCrystObj==0
   mpWXCrystObj=new WXScatterer(parent,this);
   return mpWXCrystObj;
}
#endif

}//namespace
