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
/*
*  source file LibCryst++ Crystal class
*
*/

#include <cmath>

#include <typeinfo>

#include "ObjCryst/Crystal.h"

#include "Quirks/VFNStreamFormat.h" //simple formatting of integers, doubles..
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
const RefParType *gpRefParTypeUnitCell=
   new RefParType (gpRefParTypeCrystal,"Unit Cell");
const RefParType *gpRefParTypeUnitCellLength=
   new RefParType (gpRefParTypeUnitCell,"Length");
const RefParType *gpRefParTypeUnitCellAngle=
   new RefParType (gpRefParTypeUnitCell,"Angle");

////////////////////////////////////////////////////////////////////////
//
//    CRYSTAL : the crystal (Unit cell, spaceGroup, scatterers)
//
////////////////////////////////////////////////////////////////////////
ObjRegistry<Crystal> gCrystalRegistry("List of all Crystals");

Crystal::Crystal():
mCellDim(6),mScattererRegistry("List of Crystal Scatterers"),
mBMatrix(3,3),mOrthMatrix(3,3),mOrthMatrixInvert(3,3),
mScatteringPowerRegistry("List of Crystal ScatteringPowers")
{
   VFN_DEBUG_MESSAGE("Crystal::Crystal()",10)
   this->Init(1,1,1,M_PI/2,M_PI/2,M_PI/2,"P1","");
   gCrystalRegistry.Register(*this);
   gTopRefinableObjRegistry.Register(*this);
}

Crystal::Crystal(const double a, const double b, const double c, const string &SpaceGroupId):
mCellDim(6),mScattererRegistry("List of Crystal Scatterers"),
mBMatrix(3,3),mOrthMatrix(3,3),mOrthMatrixInvert(3,3),
mScatteringPowerRegistry("List of Crystal ScatteringPowers")
{
   VFN_DEBUG_MESSAGE("Crystal::Crystal(a,b,c,Sg)",10)
   this->Init(a,b,c,M_PI/2,M_PI/2,M_PI/2,SpaceGroupId,"");
   gCrystalRegistry.Register(*this);
   gTopRefinableObjRegistry.Register(*this);
}

Crystal::Crystal(const double a, const double b, const double c, const double alpha,
              const double beta, const double gamma,const string &SpaceGroupId):
mCellDim(6),mScattererRegistry("List of Crystal Scatterers"),
mBMatrix(3,3),mOrthMatrix(3,3),mOrthMatrixInvert(3,3),
mScatteringPowerRegistry("List of Crystal ScatteringPowers")
{
   VFN_DEBUG_MESSAGE("Crystal::Crystal(a,b,c,alpha,beta,gamma,Sg)",10)
   this->Init(a,b,c,alpha,beta,gamma,SpaceGroupId,"");
   gCrystalRegistry.Register(*this);
   gTopRefinableObjRegistry.Register(*this);
}

Crystal::Crystal(const Crystal &old):
mCellDim(old.mCellDim),mSpaceGroup(old.GetSpaceGroup()),
mScattererRegistry(old.mScattererRegistry),
mBMatrix(3,3),mOrthMatrix(3,3),mOrthMatrixInvert(3,3),
mScatteringPowerRegistry(old.mScatteringPowerRegistry)
{
   VFN_DEBUG_MESSAGE("Crystal::Crystal(&oldCrystal)",10)
   this ->InitRefParList();
   for(long i=0;i<old.GetNbScatterer();i++)
   {
      mScattererRegistry.Register(*(old.GetScatt(i).CreateCopy()));
   }
   
   mUseDynPopCorr.SetChoice(old.mUseDynPopCorr.GetChoice());
   
   this->InitMatrices();
   this->UpdateLatticePar();
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

const string Crystal::GetClassName() const {return "Crystal";}

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
   VFN_DEBUG_MESSAGE("Crystal::RemoveScatteringPower()",2)
   mScatteringPowerRegistry.DeRegister(*scattPow);
   delete &scattPow;
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
         this->CalcDynPopCorr(2.,.5); else this->ResetDynPopCorr();
      VFN_DEBUG_MESSAGE("Crystal::GetScatteringComponentList():End",2)
   }
   #ifdef __DEBUG__
   if(gVFNDebugMessageLevel<2) mScattCompList.Print();
   #endif
   return mScattCompList;
}

CrystVector_double Crystal::GetLatticePar() const
{
   VFN_DEBUG_MESSAGE("Crystal::GetLatticePar()",0)
   
   if(mClockLatticeParUpdate>mClockLatticePar) return mCellDim;
   else
   {
      //:NOTE: cannot use this->UpdateLatticePar() because it is not a const member function
      int num = mSpaceGroup.GetSpaceGroupNumber();
      CrystVector_double cellDim=mCellDim;
      if(num <=2) return cellDim;
      if((num <=15) && (0==mSpaceGroup.GetUniqueAxis()))
      {
         cellDim(4)=M_PI/2.;
         cellDim(5)=M_PI/2.;
         return cellDim;
      }
      if((num <=15) && (1==mSpaceGroup.GetUniqueAxis()))
      {
         cellDim(3)=M_PI/2.;
         cellDim(5)=M_PI/2.;
         return cellDim;
      }
      if((num <=15) && (2==mSpaceGroup.GetUniqueAxis()))
      {
         cellDim(3)=M_PI/2.;
         cellDim(4)=M_PI/2.;
         return cellDim;
      }
      if(num <=74)
      {
         cellDim(3)=M_PI/2.;
         cellDim(4)=M_PI/2.;
         cellDim(5)=M_PI/2.;
         return cellDim;
      }
      if(num <= 142) 
      {
         cellDim(3)=M_PI/2.;
         cellDim(4)=M_PI/2.;
         cellDim(5)=M_PI/2.;
         cellDim(1) = mCellDim(0) ;
         return cellDim;
      }
      if(num <= 194) //Trigonal & Hexagonal
      {
         cellDim(3) = M_PI/2.;
         cellDim(4) = M_PI/2.;
         cellDim(5) = M_PI*2./3.;
         cellDim(1) = mCellDim(0) ;
         return cellDim;
      }
      cellDim(3)=M_PI/2.;
      cellDim(4)=M_PI/2.;
      cellDim(5)=M_PI/2.;
      cellDim(1) = mCellDim(0) ;
      cellDim(2) = mCellDim(0) ;
      return cellDim;
   }
}

double Crystal::GetLatticePar(int whichPar)const
{
   VFN_DEBUG_MESSAGE("Crystal::LatticePar(i)",0)
   if( (whichPar<0) || (whichPar>5))
      throw ObjCrystException("Crystal::LatticePar(int) :trying to access parameter>5!");
      
   if(mClockLatticeParUpdate>mClockLatticePar) return mCellDim(whichPar);
   else
   {
      const int num = mSpaceGroup.GetSpaceGroupNumber();

      static CrystVector_double cellDim;
      cellDim=mCellDim;
      if(num <=2)
         return cellDim(whichPar);
      if((num <=15) && (0==mSpaceGroup.GetUniqueAxis()))
      {
         cellDim(4)=M_PI/2.;
         cellDim(5)=M_PI/2.;
         return cellDim(whichPar);
      }
      if((num <=15) && (1==mSpaceGroup.GetUniqueAxis()))
      {
         cellDim(3)=M_PI/2.;
         cellDim(5)=M_PI/2.;
         return cellDim(whichPar);
      }
      if((num <=15) && (2==mSpaceGroup.GetUniqueAxis()))
      {
         cellDim(3)=M_PI/2.;
         cellDim(4)=M_PI/2.;
         return cellDim(whichPar);
      }

      if(num <=74)
      {
         cellDim(3)=M_PI/2.;
         cellDim(4)=M_PI/2.;
         cellDim(5)=M_PI/2.;
         return cellDim(whichPar);
      }
      if(num <= 142) 
      {
         cellDim(3)=M_PI/2.;
         cellDim(4)=M_PI/2.;
         cellDim(5)=M_PI/2.;
         cellDim(1) = mCellDim(0) ;
         return cellDim(whichPar);
      }
      if(num <= 194) 
      {
         cellDim(3) = M_PI/2.;
         cellDim(4) = M_PI/2.;
         cellDim(5) = M_PI*2./3.;
         cellDim(1) = mCellDim(0) ;
         return cellDim(whichPar);
      }
      cellDim(3)=M_PI/2.;
      cellDim(4)=M_PI/2.;
      cellDim(5)=M_PI/2.;
      cellDim(1) = mCellDim(0) ;
      cellDim(2) = mCellDim(0) ;
      return cellDim(whichPar);
   }
}

const CrystMatrix_double& Crystal::GetBMatrix()const
{
   VFN_DEBUG_MESSAGE("Crystal::GetBMatrix()",0)
   this->InitMatrices();
   return mBMatrix;
}

CrystVector_double Crystal::GetOrthonormalCoords(const double x,
                                                const double y,
                                                const double z) const
{
   this->InitMatrices();
   CrystVector_double coords(3);
   coords(0)=mOrthMatrix(0,0)*x+mOrthMatrix(0,1)*y+mOrthMatrix(0,2)*z;
   coords(1)=mOrthMatrix(1,0)*x+mOrthMatrix(1,1)*y+mOrthMatrix(1,2)*z;
   coords(2)=mOrthMatrix(2,0)*x+mOrthMatrix(2,1)*y+mOrthMatrix(2,2)*z;
   return coords;
}

void Crystal::FractionalToOrthonormalCoords(double &x,double &y,double &z) const
{
   this->InitMatrices();
   const double oldx=x;
   const double oldy=y;
   x=mOrthMatrix(0,0)*oldx+mOrthMatrix(0,1)*oldy+mOrthMatrix(0,2)*z;
   y=mOrthMatrix(1,0)*oldx+mOrthMatrix(1,1)*oldy+mOrthMatrix(1,2)*z;
   z=mOrthMatrix(2,0)*oldx+mOrthMatrix(2,1)*oldy+mOrthMatrix(2,2)*z;
}

void Crystal::OrthonormalToFractionalCoords(double &x,double &y,double &z) const
{
   //cout << x << " " << y << " " << z <<endl;
   //cout << endl << mOrthMatrixInvert <<endl;
   this->InitMatrices();
   const double oldx=x;
   const double oldy=y;
   x=mOrthMatrixInvert(0,0)*oldx+mOrthMatrixInvert(0,1)*oldy+mOrthMatrixInvert(0,2)*z;
   y=mOrthMatrixInvert(1,0)*oldx+mOrthMatrixInvert(1,1)*oldy+mOrthMatrixInvert(1,2)*z;
   z=mOrthMatrixInvert(2,0)*oldx+mOrthMatrixInvert(2,1)*oldy+mOrthMatrixInvert(2,2)*z;
   //cout << x << " " << y << " " << z <<endl;
}

void Crystal::Print(ostream &os)const
{
   VFN_DEBUG_MESSAGE("Crystal::Print()",5)
   this->InitMatrices();
   int width =8 ;
   os << "crystal : " << mName <<"("<<this->GetSpaceGroup().GetName()<<")"<< endl;
   os.width(width);
   os   << "    Cell dimensions : " 
        << FormatFloat(mCellDim(0),8,5) << "  " 
        << FormatFloat(mCellDim(1),8,5) << "  " 
        << FormatFloat(mCellDim(2),8,5) << "  "
        << FormatFloat(mCellDim(3)*RAD2DEG,8,5) << "  " 
        << FormatFloat(mCellDim(4)*RAD2DEG,8,5) << "  " 
        << FormatFloat(mCellDim(5)*RAD2DEG,8,5) << endl ;
   //mSpaceGroup.Print();
   
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
   double nbAtoms=0;
   const long genMult=mSpaceGroup.GetNbSymmetrics();
   for(int i=0;i<mScattCompList.GetNbComponent();i++) 
      nbAtoms += genMult * mScattCompList(i).mOccupancy * mScattCompList(i).mDynPopCorr;
   os << " Total number of components (atoms) in one unit cell : " << nbAtoms<<endl<<endl;
   
   VFN_DEBUG_MESSAGE("Crystal::Print():End",5)
}

const SpaceGroup & Crystal::GetSpaceGroup() const {return mSpaceGroup;}
SpaceGroup & Crystal::GetSpaceGroup()  {return mSpaceGroup;}
 
CrystMatrix_double Crystal::GetMinDistanceTable(const double minDistance) const
{
   VFN_DEBUG_MESSAGE("Crystal::MinDistanceTable()",5)
   this->CalcDistTable(true);
   const long nbComponent=mScattCompList.GetNbComponent();
   CrystMatrix_double minDistTable(nbComponent,nbComponent);
   const int nbAtoms=mDistTableIndex.numElements();
   double dist;
   double tmp;
   const double min=minDistance*minDistance;
   //cout<<mDistTableSq.rows()<<" "<<mDistTableSq.cols()<<" "<<nbComponent<<" "<<nbAtoms<<endl;
   minDistTable=10000.;
   for(int i=0;i<nbComponent;i++)
   {
      //VFN_DEBUG_MESSAGE("Crystal::MinDistanceTable():comp="<<i,10)
      for(int j=0;j<nbAtoms;j++)
      {
         tmp=mDistTableSq(j,i);
         dist=minDistTable(i,mDistTableIndex(j));
         if( (tmp<dist) && ((tmp>min) || (i!=j))) minDistTable(i,mDistTableIndex(j))=tmp;
      }
   }
   //cout << mDistTableSq<<endl<<minDistTable<<endl;
   for(int i=0;i<nbComponent;i++)
   {
      for(int j=0;j<=i;j++)
      {
         if(9999.<minDistTable(i,j)) minDistTable(i,j)=0;
         minDistTable(i,j)=sqrt(minDistTable(i,j));
         minDistTable(j,i)=minDistTable(i,j);
      }
   }
   //cout << mDistTableSq<<endl<<minDistTable<<endl;
   VFN_DEBUG_MESSAGE("Crystal::MinDistanceTable():End",3)
   return minDistTable;
}

void Crystal::PrintMinDistanceTable(const double minDistance,ostream &os) const
{
   VFN_DEBUG_MESSAGE("Crystal::PrintMinDistanceTable()",5)
   CrystMatrix_double minDistTable;
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
   double x=0.5;
   double y=0.5;
   double z=0.5;
   this->FractionalToOrthonormalCoords(x,y,z);
   os << "    look_at   <" << x << "," << y << "," << z <<">"<<endl;
   double maxDim=x;
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
                                const double xMin,const double xMax,
                                const double yMin,const double yMax,
                                const double zMin,const double zMax)const
{
   VFN_DEBUG_ENTRY("Crystal::GLInitDisplayList()",5)
   #ifdef OBJCRYST_GL
   glMaterialf (GL_FRONT, GL_SHININESS, 50);
   GLfloat colour_axis[]= { 0.5, .5, .5, 1.0 };
   glMaterialfv (GL_FRONT, GL_AMBIENT_AND_DIFFUSE,colour_axis);
	//cout << xMin << ":"<<xMax <<endl;
	//cout << yMin << ":"<<yMax <<endl;
	//cout << zMin << ":"<<zMax <<endl;
	//Center of displayed unit
		double xc=(xMin+xMax)/2.;
		double yc=(yMin+yMax)/2.;
		double zc=(zMin+zMax)/2.;
   	//Describe Unit Cell
      	double x111= 1.;
      	double y111= 1.;
      	double z111= 1.;
      	this->FractionalToOrthonormalCoords(x111,y111,z111);
      	double x110= 1.;
      	double y110= 1.;
      	double z110= 0.;
      	this->FractionalToOrthonormalCoords(x110,y110,z110);
      	double x101= 1.;
      	double y101= 0.;
      	double z101= 1.;
      	this->FractionalToOrthonormalCoords(x101,y101,z101);
      	double x100= 1.;
      	double y100= 0.;
      	double z100= 0.;
      	this->FractionalToOrthonormalCoords(x100,y100,z100);
      	double x011= 0.;
      	double y011= 1.;
      	double z011= 1.;
      	this->FractionalToOrthonormalCoords(x011,y011,z011);
      	double x010= 0.;
      	double y010= 1.;
      	double z010= 0.;
      	this->FractionalToOrthonormalCoords(x010,y010,z010);
      	double x001= 0.;
      	double y001= 0.;
      	double z001= 1.;
      	this->FractionalToOrthonormalCoords(x001,y001,z001);
      	double x000= 0.;
      	double y000= 0.;
      	double z000= 0.;
      	this->FractionalToOrthonormalCoords(x000,y000,z000);
      glPushMatrix();
   	//Add Axis & axis names
      	GLfloat colour_font[]= { 1.0,1.0,1.0, 1.0 };
      	glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,colour_font);
      	//glMaterialfv (GL_FRONT, GL_EMISSION,colour_font);
      	double x,y,z;

      	x=1.2-xc;y=-yc;z=-zc;
      	this->FractionalToOrthonormalCoords(x,y,z);
      	glRasterPos3f(x,y,z);
      	glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18,'a');

      	x=-xc;y=1.2-yc;z=-zc;
      	this->FractionalToOrthonormalCoords(x,y,z);
      	glRasterPos3f(x,y,z);
      	glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18,'b');

      	x=-xc;y=-yc;z=1.2-zc;
      	this->FractionalToOrthonormalCoords(x,y,z);
      	glRasterPos3f(x,y,z);
      	glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18,'c');
		// Cell
      	this->FractionalToOrthonormalCoords(xc,yc,zc);
         glTranslatef(-xc, -yc, -zc);
      	glColor3f(1.0f,1.0f,1.0f);	   // White
      	glBegin(GL_LINE_LOOP);				// Bottom
         	glVertex3f(    x110,    y110,    z110);
         	glVertex3f(    x010,    y010,    z010);
         	glVertex3f(    x011,    y011,    z011);
         	glVertex3f(    x111,    y111,    z111);
      	glEnd();	
      	glBegin(GL_LINE_LOOP);				// Top
         	glVertex3f(    x101,    y101,    z101);
         	glVertex3f(    x001,    y001,    z001);
         	glVertex3f(    x000,    y000,    z000);
         	glVertex3f(    x100,    y100,    z100);
      	glEnd();
      	glBegin(GL_LINES);				// Top
         	glVertex3f(    x101,    y101,    z101);
         	glVertex3f(    x111,    y111,    z111);

         	glVertex3f(    x001,    y001,     z001);
         	glVertex3f(    x011,    y011,     z011);

         	glVertex3f(    x000,    y000,    z000);
         	glVertex3f(    x010,    y010,    z010);

         	glVertex3f(    x100,    y100,    z100);
         	glVertex3f(    x110,    y110,    z110);
      	glEnd();


   	//Describe all Scatterers
   	VFN_DEBUG_MESSAGE("Crystal::GLView(bool):Scatterers...",5)
   	for(int i=0;i<mScattererRegistry.GetNb();i++) 
      	this->GetScatt(i).GLInitDisplayList(onlyIndependentAtoms,
                                          	xMin,xMax,yMin,yMax,zMin,zMax);

   glPopMatrix();
   #else
   cout << "Crystal::GLView(): Compiled without OpenGL support !" <<endl;
   #endif
   VFN_DEBUG_EXIT("Crystal::GLInitDisplayList(bool)",5)
}

void Crystal::CalcDynPopCorr(const double overlapDist, const double mergeDist) const
{
   VFN_DEBUG_MESSAGE("Crystal::CalcDynPopCorr(double)",4)
   TAU_PROFILE("Crystal::CalcDynPopCorr()","void (double)",TAU_DEFAULT);
   
   this->CalcDistTable(true,overlapDist);
   if(mClockDynPopCorr>mClockNeighborTable) return;
   
   const long nbComponent=mScattCompList.GetNbComponent();
   const int nbSymmetrics=mSpaceGroup.GetNbSymmetrics();
   CrystVector_double neighborsDist(nbComponent*nbSymmetrics);
   long nbNeighbors=0;
   double corr;
   
   int atomicNumber;
   const long nbAtomSym=mDistTableIndex.numElements();
   const double overlapDistSq=overlapDist*overlapDist;
   double dist;
   for(long i=0;i<nbComponent;i++)
   {
      atomicNumber=mScattCompList(i).mpScattPow->GetDynPopCorrIndex();
      nbNeighbors=0;
      for(long j=0;j<nbAtomSym;j++)
         if(atomicNumber==mScattCompList(mDistTableIndex(j)).mpScattPow->GetDynPopCorrIndex())
         {
            if(overlapDistSq > mDistTableSq(j,i))
               neighborsDist(nbNeighbors++)=sqrt(mDistTableSq(j,i));
         }
      corr=0.;
      if(nbNeighbors==0)
      {
         cout << "Crystal::CalcDynPopCorr():no neighbors !"<<i<<":"<< sqrt(overlapDistSq)<<endl;
         #ifdef __DEBUG__
         cout <<mDistTableSq<<endl<<endl;
         for(long j=0;j<nbAtomSym;j++)
            if(atomicNumber
                ==mScattCompList(mDistTableIndex(j)).mpScattPow->GetDynPopCorrIndex())
            {
               cout << sqrt(mDistTableSq(j,i)) <<" "<<(overlapDistSq>mDistTableSq(j,i)) <<endl;
            }
         #endif
         abort();
      }
      for(long j=0;j<nbNeighbors;j++) 
      {
         dist=neighborsDist(j)-mergeDist;
         if(dist<0.) dist=0.;
         corr += fabs((overlapDist-dist-mergeDist)/(overlapDist-mergeDist));
      }
      mScattCompList(i).mDynPopCorr=1./corr;
   }
   mClockDynPopCorr.Click();
   VFN_DEBUG_MESSAGE("Crystal::CalcDynPopCorr(double):End.",3)
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

double Crystal::GetBumpMergeCostFunction() const
{
   VFN_DEBUG_MESSAGE("Crystal::GetBumpMergeCostFunction()",2)
   TAU_PROFILE("Crystal::GetBumpMergeCostFunction()","double (double)",TAU_DEFAULT);
   
   if(mBumpDistanceMatrix.numElements()==0) return 0;

   this->GetScatteringComponentList();//update list if necessary
   
   const long nbComponent=mScattCompList.GetNbComponent();
   //const int nbSymmetrics=mSpaceGroup.GetNbSymmetrics();
   long nbNeighbors=0;
   
   //this->CalcDistTable(true,overlapDist);//:KLUDGE: Assume computing has already been done
   
   //Average number of neighbors, used to normalize the cost function (max env. =1)
   //double avNbNeighbor;
   //{
   // const double a=GetLatticePar(0);
   // const double b=GetLatticePar(1);
   // const double c=GetLatticePar(2);
   // const double alpha=GetLatticePar(3);
   // const double beta=GetLatticePar(4);
   // const double gamma=GetLatticePar(5);
   // const double v=a*b*c*sqrt(1-cos(alpha)*cos(alpha)-cos(beta)*cos(beta)-cos(gamma)*cos(gamma)
   //            +2*cos(alpha)*cos(beta)*cos(gamma));
   // const double m=3;//sqrt(mBumpDistanceMatrix.max());
   // avNbNeighbor=nbComponent*nbSymmetrics*4./3.*M_PI*m*m*m/v;
   //}
   
   const long nbAtomSym=mDistTableIndex.numElements();
   //Init the cubic splines
      CrystVector_double x(3),y(3);
      x(0)=0.;x(1)=.5;x(2)=1.;
      y(0)=0.;y(1)=.5;y(2)=0.;
      static const CubicSpline splineBumpMerge(x,y,0,0);
      x.resize(2);
      y.resize(2);
      x(0)=.5;x(1)=1.;
      y(0)=1.;y(1)=0.;
      static const CubicSpline splineBump(x,y,0,0);
   double cost=0;
   
   CrystVector_int atomicNumber(nbComponent);
   for(long i=0;i<nbComponent;i++)
      atomicNumber(i)=mScattCompList(i).mpScattPow->GetDynPopCorrIndex();
   CrystVector_int atomicNumber2(nbAtomSym);
   for(long j=0;j<nbAtomSym;j++)
      atomicNumber2(j)=mScattCompList(mDistTableIndex(j)).mpScattPow->GetDynPopCorrIndex();

   //cout << nbComponent <<" "<<nbAtomSym<<endl;
   //cout << mDistTableSq.cols()<<" "<<mDistTableSq.rows()<<endl;
   
   //:TODO: faster store the splines in arrays
   double tmp;
   for(long i=0;i<nbComponent;i++)
   {
      nbNeighbors=0;
      //cout << mScattCompList(i).mpScattPow->GetName()<<":::";
      for(long j=0;j<nbAtomSym;j++)
      {
         if(mBumpDistanceMatrix(atomicNumber(i),atomicNumber2(j))<.1) continue;
         mDistTableSq(j,i);
         tmp=mDistTableSq(j,i)/mBumpDistanceMatrix(atomicNumber(i),atomicNumber2(j));
         if(tmp>=1.) continue;
         //cout <<mScattCompList(mDistTableIndex(j)).mpScattPow->GetName()<<" "
         //     <<sqrt(mDistTableSq(j,i))<<"->";
         //if(atomicNumber(i)==atomicNumber2(j))
         if(true==mAllowMerge(atomicNumber(i),atomicNumber2(j)))
         {
            cost+= splineBumpMerge(sqrt(tmp));
            //cout <<splineBumpMerge(sqrt(tmp))<<"M :: ";
         }
         else
         {
            cost+= splineBump(sqrt(tmp));
            //cout <<splineBump(sqrt(tmp))<<"B :: ";
         }
      }
      //cout <<endl;
   }
   //cout << "BumpMergeCost:"<<cost<<" "<<avNbNeighbor<<endl;
   cost /= nbComponent*nbComponent;//avNbNeighbor*nbComponent;
   //cout << "BumpMergeCost:"<<cost<<endl;
   //cout << mDistTableSq<<endl;
   //this->PrintMinDistanceTable(.001);
   //abort();
   VFN_DEBUG_MESSAGE("Crystal::GetBumpMergeCostFunction():End",2)
   return cost;
}

void Crystal::SetBumpMergeDistance(const ScatteringPower &scatt1,
                                   const ScatteringPower &scatt2,const double dist)
{
   VFN_DEBUG_MESSAGE("Crystal::SetBumpMergeDistance()",5)
   const int num1=scatt1.GetDynPopCorrIndex();
   const int num2=scatt2.GetDynPopCorrIndex();
   if(num1==num2) this->SetBumpMergeDistance(scatt1,scatt2,dist,true) ;
   else this->SetBumpMergeDistance(scatt1,scatt2,dist,false);
}

void Crystal::SetBumpMergeDistance(const ScatteringPower &scatt1,
                                   const ScatteringPower &scatt2,const double dist,
                                   const bool allowMerge)
{
   VFN_DEBUG_MESSAGE("Crystal::SetBumpMergeDistance()",5)
   const int num1=scatt1.GetDynPopCorrIndex();
   const int num2=scatt2.GetDynPopCorrIndex();
   if(mBumpDistanceMatrix.numElements()==0) 
   {
      mBumpDistanceMatrix.resize(200,200);
      mBumpDistanceMatrix=0;
      mAllowMerge.resize(200,200);
      mAllowMerge=false;
   }
   VFN_DEBUG_MESSAGE("Crystal::SetBumpMergeDistance():"<<num1<<","<<num2,8)
   mBumpDistanceMatrix(num1,num2)=dist*dist;
   mBumpDistanceMatrix(num2,num1)=dist*dist;
   mAllowMerge(num1,num2)=allowMerge;
   mAllowMerge(num2,num1)=allowMerge;
}

const RefinableObjClock& Crystal::GetClockLatticePar()const {return mClockLatticePar;}
const RefinableObjClock& Crystal::GetClockScattererList()const {return mClockScattererList;}

unsigned int Crystal::GetNbCostFunction()const {return 1;}

const string& Crystal::GetCostFunctionName(const unsigned int id)const
{
   static string costFunctionName[2];
   if(0==costFunctionName[0].length())
   {
      costFunctionName[0]="BumpMergeCost()";
   }
   switch(id)
   {
      case 0: return costFunctionName[0];
      default:
      {
         cout << "Crystal::GetCostFunctionName(): Not Found !" <<endl;
         throw 0;
      }
   }
}

const string& Crystal::GetCostFunctionDescription(const unsigned int id)const
{
   static string costFunctionDescription[2];
   if(0==costFunctionDescription[0].length())
   {
      costFunctionDescription[0]="Anti-bump, pro-merge cost function";
   }
   switch(id)
   {
      case 0: return costFunctionDescription[0];
      default:
      {
         cout << "Crystal::GetCostFunctionDescription(): Not Found !" 
              <<endl;
         throw 0;
      }
   }
}

double Crystal::GetCostFunctionValue(const unsigned int n)
{
   VFN_DEBUG_MESSAGE("Crystal::GetCostFunctionValue():"<<mName,4)
   switch(n)
   {
      case 0: return this->GetBumpMergeCostFunction();
      default:
      {
         cout << "Crystal::GetCostFunctionValue(): Not Found !" <<endl;
         throw 0;
      }
   }
}

void Crystal::GlobalOptRandomMove(const double mutationAmplitude)
{
   VFN_DEBUG_MESSAGE("Crystal::GlobalOptRandomMove()",2)
   //Either a random move or a permutation of two scatterers
   if( (rand()/(double)RAND_MAX)<-.1)//disabled :TODO: exclude fixed scatterers
   {
      const long nb=this->GetNbScatterer();
      const unsigned long n1=(unsigned long)((rand()/(double)(RAND_MAX-1))*nb);
      const long n2=( (long)((rand()/(double)(RAND_MAX-1))*(nb-1))+1+n1)%nb;
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
      this->RefinableObj::GlobalOptRandomMove(mutationAmplitude);
   }
   VFN_DEBUG_MESSAGE("Crystal::GlobalOptRandomMove():End",1)
}

void Crystal::CIFOutput(ostream &os)const
{
   VFN_DEBUG_ENTRY("Crystal::OutputCIF()",5)
   this->InitMatrices();
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
      
   os << "_cell_length_a "   << FormatFloat(mCellDim(0),8,5) << endl
      << "_cell_length_b "   << FormatFloat(mCellDim(1),8,5) << endl
      << "_cell_length_c "   << FormatFloat(mCellDim(2),8,5) << endl
      << "_cell_angle_alpha "<< FormatFloat(mCellDim(3)*RAD2DEG,8,5) << endl 
      << "_cell_angle_beta " << FormatFloat(mCellDim(4)*RAD2DEG,8,5) << endl 
      << "_cell_angle_gamma "<< FormatFloat(mCellDim(5)*RAD2DEG,8,5) << endl ;
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

void Crystal::Init(const double a, const double b, const double c, const double alpha,
                   const double beta, const double gamma,const string &SpaceGroupId,
                   const string& name)
{
   VFN_DEBUG_MESSAGE("Crystal::Init(a,b,c,alpha,beta,gamma,Sg,name)",10)
   //mSpaceGroup.Print();
   mSpaceGroup.ChangeSpaceGroup(SpaceGroupId);
   VFN_DEBUG_MESSAGE("Crystal::Init(a,b,c,alpha,beta,gamma,Sg,name)1",5)
   //mSpaceGroup.Print();
   mCellDim(0)=a;
   mCellDim(1)=b;
   mCellDim(2)=c;
   mCellDim(3)=alpha;
   mCellDim(4)=beta;
   mCellDim(5)=gamma;
   
   mClockMetricMatrix.Reset();
   mClockScattCompList.Reset();
   mClockNeighborTable.Reset();
   mClockLatticeParUpdate.Reset();
   mClockDynPopCorr.Reset();
   
   mClockLatticePar.Click();
   VFN_DEBUG_MESSAGE("Crystal::Init(a,b,c,alpha,beta,gamma,Sg,name)2",5)
   
   static string UseDynPopCorrname;
   static string UseDynPopCorrchoices[2];
   static bool needInitNames=true;
   if(true==needInitNames)
   {
      UseDynPopCorrname="Use Dynamical Occupancy Correction";
      UseDynPopCorrchoices[0]="No";
      UseDynPopCorrchoices[1]="Yes";
      needInitNames=false;//Only once for the class
   }
   VFN_DEBUG_MESSAGE("Crystal::Init(a,b,c,alpha,beta,gamma,Sg,name):Init options",5)
   mUseDynPopCorr.Init(2,&UseDynPopCorrname,UseDynPopCorrchoices);
   mUseDynPopCorr.SetChoice(1);
   this->AddOption(&mUseDynPopCorr);
   
   this->InitRefParList();
   this->InitMatrices();
   this->SetName(name);
   
   VFN_DEBUG_MESSAGE("Crystal::Init(a,b,c,alpha,beta,gamma,Sg,name):End",10)
}

void Crystal::InitMatrices() const
{
   //:NOTE: The Matrices must remain upper triangular, since this is assumed for
   //optimization purposes in some procedures.
   if(mClockMetricMatrix>mClockLatticePar) return;//no need to update
   //this->UpdateLatticePar(); we should be able to do this...
   
   VFN_DEBUG_MESSAGE("Crystal::InitMatrices() for crystal : "+this->GetName(),5)
   //mClockMetricMatrix.Print();
   //mClockLatticePar.Print();

   double a,b,c,alpha,beta,gamma;//direct space parameters
   double aa,bb,cc,alphaa,betaa,gammaa;//reciprocal space parameters
   double v;//volume of the unit cell
   a=GetLatticePar(0);
   b=GetLatticePar(1);
   c=GetLatticePar(2);
   alpha=GetLatticePar(3);
   beta=GetLatticePar(4);
   gamma=GetLatticePar(5);
   
   //cout <<a<<" "<<b<<" "<<c<<" "<<alpha<<" "<<beta<<" "<<gamma<<endl;
   
   v=sqrt(1-cos(alpha)*cos(alpha)-cos(beta)*cos(beta)-cos(gamma)*cos(gamma)
            +2*cos(alpha)*cos(beta)*cos(gamma));
            
   aa=sin(alpha)/a/v;
   bb=sin(beta )/b/v;
   cc=sin(gamma)/c/v;
   
   alphaa=acos( (cos(beta )*cos(gamma)-cos(alpha))/sin(beta )/sin(gamma) );
   betaa =acos( (cos(alpha)*cos(gamma)-cos(beta ))/sin(alpha)/sin(gamma) );
   gammaa=acos( (cos(alpha)*cos(beta )-cos(gamma))/sin(alpha)/sin(beta ) );
   
   //cout <<aa<<" "<<bb<<" "<<cc<<" "<<alphaa<<" "<<betaa<<" "<<gammaa<<endl;

   // Formula is OK
   mBMatrix = aa ,  bb*cos(gammaa) , cc*cos(betaa) ,
               0  , bb*sin(gammaa) ,-cc*sin(betaa)*cos(alpha),
               0  , 0              ,1/c;
   //cout <<"B Matrix :"<<endl<< mBMatrix <<endl;
   // Formula is OK
   mOrthMatrix = a  , b*cos(gamma) , c*cos(beta) ,
                 0  , b*sin(gamma) ,-c*sin(beta)*cos(alphaa),
                 0  , 0              ,1/cc;
   
   //Seems OK- NO !
   //mOrthMatrixInvert= 1/a, -1/tan(gamma)/a , -cos(beta)/a/v/sin(gamma)*(1+cos(alpha)*cos(gamma)),
   //						 0  , 1/b/sin(gamma)  , bb*cos(alphaa),
   //						 0  , 0  				 , cc;
	mOrthMatrixInvert=InvertMatrix(mOrthMatrix);
   //cout << "Orth Matrix :"<<endl<<mOrthMatrix <<endl;
   //cout << "InvOrth Matrix :"<<endl<<mOrthMatrixInvert <<endl;
   //cout << "Orth * InvOrth matrix :"<<endl<<product(mOrthMatrix,mOrthMatrixInvert) <<endl;
   //cout << "InvOrth * Orth Matrix :"<<endl<<product(mOrthMatrixInvert,mOrthMatrix) <<endl;
   mClockMetricMatrix.Click();
   VFN_DEBUG_MESSAGE("Crystal::InitMatrices():End.",5)
}

void Crystal::UpdateLatticePar()
{
   if(  (mClockLatticeParUpdate>mSpaceGroup.GetClockSpaceGroup())
      &&(mClockLatticeParUpdate>mClockLatticePar)) return;
   VFN_DEBUG_MESSAGE("Crystal::UpdateLatticePar().",3)

   int num = mSpaceGroup.GetSpaceGroupNumber();
   if(num <=2)
   {
      mClockLatticeParUpdate.Click();
      return;
   }
   if((num <=15) && (0==mSpaceGroup.GetUniqueAxis()))
   {
      mCellDim(4)=M_PI/2.;
      mCellDim(5)=M_PI/2.;
      mClockLatticeParUpdate.Click();
      return;
   }
   if((num <=15) && (1==mSpaceGroup.GetUniqueAxis()))
   {
      mCellDim(3)=M_PI/2.;
      mCellDim(5)=M_PI/2.;
      mClockLatticeParUpdate.Click();
      return;
   }
   if((num <=15) && (2==mSpaceGroup.GetUniqueAxis()))
   {
      mCellDim(3)=M_PI/2.;
      mCellDim(4)=M_PI/2.;
      mClockLatticeParUpdate.Click();
      return;
   }
   if(num <=74)
   {
      mCellDim(3)=M_PI/2.;
      mCellDim(4)=M_PI/2.;
      mCellDim(5)=M_PI/2.;
      mClockLatticeParUpdate.Click();
      return;
   }
   if(num <= 142) 
   {
      mCellDim(3)=M_PI/2.;
      mCellDim(4)=M_PI/2.;
      mCellDim(5)=M_PI/2.;
      mCellDim(1) = mCellDim(0) ;
      mClockLatticeParUpdate.Click();
      return;
   }
   if(num <= 194) 
   {
      mCellDim(3) = M_PI/2.;
      mCellDim(4) = M_PI/2.;
      mCellDim(5) = M_PI*2./3.;
      mCellDim(1) = mCellDim(0) ;
      mClockLatticeParUpdate.Click();
      return;
   }
   mCellDim(3)=M_PI/2.;
   mCellDim(4)=M_PI/2.;
   mCellDim(5)=M_PI/2.;
   mCellDim(1) = mCellDim(0) ;
   mCellDim(2) = mCellDim(0) ;
   mClockLatticeParUpdate.Click();
   return;
}

void Crystal::InitRefParList()
{
   VFN_DEBUG_MESSAGE("Crystal::InitRefParList()",5)
   //this->ResetParList();
   int num = mSpaceGroup.GetSpaceGroupNumber();
   bool a=true;
   bool b=true;
   bool c=true;
   bool alpha=true;
   bool beta=true;
   bool gamma=true;
   if(num <=2)
   {
   }
   else if((num <=15) && (0==mSpaceGroup.GetUniqueAxis()))
   {
      beta=false;
      gamma=false;
   }
   else if((num <=15) && (1==mSpaceGroup.GetUniqueAxis()))
   {
      alpha=false;
      gamma=false;
   }
   else if((num <=15) && (2==mSpaceGroup.GetUniqueAxis()))
   {
      alpha=false;
      beta=false;
   }
   else if(num <=74)
   {
      alpha=false;
      beta=false;
      gamma=false;
   }
   else if(num <= 142) 
   {
      b=false;
      alpha=false;
      beta=false;
      gamma=false;
   }
   else if(num <= 167)
   {//Hexagonal axes !
      b=false;
      alpha=false;
      beta=false;
      gamma=false;
   }
   else if(num <= 194) 
   {
      b=false;
      alpha=false;
      beta=false;
      gamma=false;
   }
   else
   {
      b=false;
      c=false;
      alpha=false;
      beta=false;
      gamma=false;
   }
   double *pLatPar=mCellDim.data();
   if(this->GetNbPar()==0)
   {//:KLUDGE:
      {
         RefinablePar tmp(this->GetName()+"_A",pLatPar,1.,100.,
                           gpRefParTypeUnitCellLength,REFPAR_DERIV_STEP_ABSOLUTE,
                           true,true,a,false,1.0);
         tmp.SetDerivStep(1e-4);
         tmp.AssignClock(mClockLatticePar);
         this->AddPar(tmp);
      }
      {
         RefinablePar tmp(this->GetName()+"_B",pLatPar+1,1.,100.,
                           gpRefParTypeUnitCellLength,REFPAR_DERIV_STEP_ABSOLUTE,
                           true,true,b,false,1.0);
         tmp.SetDerivStep(1e-4);
         tmp.AssignClock(mClockLatticePar);
         this->AddPar(tmp);
      }
      {
         RefinablePar tmp(this->GetName()+"_C",pLatPar+2,1.,100.,
                           gpRefParTypeUnitCellLength,REFPAR_DERIV_STEP_ABSOLUTE,
                           true,true,c,false,1.0);
         tmp.SetDerivStep(1e-4);
         tmp.AssignClock(mClockLatticePar);
         this->AddPar(tmp);
      }
      {
         RefinablePar tmp(this->GetName()+"_ALPHA",pLatPar+3,.5,3.,
                           gpRefParTypeUnitCellAngle,REFPAR_DERIV_STEP_ABSOLUTE,
                           true,true,alpha,false,RAD2DEG);
         tmp.SetDerivStep(1e-4);
         tmp.AssignClock(mClockLatticePar);
         this->AddPar(tmp);
      }
      {
         RefinablePar tmp(this->GetName()+"_BETA",pLatPar+4,.5,3.,
                           gpRefParTypeUnitCellAngle,REFPAR_DERIV_STEP_ABSOLUTE,
                           true,true,beta,false,RAD2DEG);
         tmp.SetDerivStep(1e-4);
         tmp.AssignClock(mClockLatticePar);
         this->AddPar(tmp);
      }
      {
         RefinablePar tmp(this->GetName()+"_GAMMA",pLatPar+5,.5,3.,
                           gpRefParTypeUnitCellAngle,REFPAR_DERIV_STEP_ABSOLUTE,
                           true,true,gamma,false,RAD2DEG);
         tmp.SetDerivStep(1e-4);
         tmp.AssignClock(mClockLatticePar);
         this->AddPar(tmp);
      }
   }
   else
   {//Just Fix the 'used' status
      this->GetPar(pLatPar+0).SetIsUsed(a);
      this->GetPar(pLatPar+1).SetIsUsed(b);
      this->GetPar(pLatPar+2).SetIsUsed(c);
      this->GetPar(pLatPar+3).SetIsUsed(alpha);
      this->GetPar(pLatPar+4).SetIsUsed(beta);
      this->GetPar(pLatPar+5).SetIsUsed(gamma);
   }

   VFN_DEBUG_MESSAGE("Crystal::InitRefParList():Finished",5)
}

void Crystal::CalcDistTable(const bool fast, const double asymUnitMargin) const
{
   this->GetScatteringComponentList();
   this->InitMatrices();
   
   if(  (mClockNeighborTable>mClockScattCompList)
      &&(mClockNeighborTable>mClockMetricMatrix)) return;
   
   const long nbComponent=mScattCompList.GetNbComponent();
   if(true==fast)
   {
      VFN_DEBUG_MESSAGE("Crystal::CalcDistTable(fast=true)",3)
      TAU_PROFILE("Crystal::CalcDistTable(fast=true)","Matrix (string&)",TAU_DEFAULT);
      TAU_PROFILE_TIMER(timer1,"DiffractionData::CalcDistTable1","", TAU_FIELD);
      TAU_PROFILE_TIMER(timer2,"DiffractionData::CalcDistTable2","", TAU_FIELD);
      TAU_PROFILE_TIMER(timer3,"DiffractionData::CalcDistTable3","", TAU_FIELD);
      TAU_PROFILE_TIMER(timer4,"DiffractionData::CalcDistTable4","", TAU_FIELD);

      TAU_PROFILE_START(timer1);
      const long intScale=16384;//Do not change this. Used later with AND (&)
      const long intScale2=8192;
      const int nbSymmetrics=mSpaceGroup.GetNbSymmetrics();
      long nbAtoms=0;
      CrystVector_long xCoords(nbComponent*nbSymmetrics);
      CrystVector_long yCoords(nbComponent*nbSymmetrics);
      CrystVector_long zCoords(nbComponent*nbSymmetrics);
      CrystVector_long xCoords2(nbComponent);
      CrystVector_long yCoords2(nbComponent);
      CrystVector_long zCoords2(nbComponent);
      CrystMatrix_double symmetricsCoords;
      mDistTableIndex.resize(nbComponent*nbSymmetrics);
      long *pIndex=mDistTableIndex.data();
      //get the fractionnal coordinates of all atoms
      double junk;
      for(long i=0;i<nbComponent;i++)
      {
         symmetricsCoords=mSpaceGroup.GetAllSymmetrics(mScattCompList(i).mX,
                                                       mScattCompList(i).mY,
                                                       mScattCompList(i).mZ);
         
         //:KLUDGE: Ensure that all coordinates are >0. This assumes that fractional
         //coordinates are > -8.
         symmetricsCoords += 8;
         
         // Ensure that coords are all between 0 and 1
         for(int j=0;j<nbSymmetrics;j++)
         {//'Optimizing' this loop with pointers is useless
            xCoords(nbAtoms)=(long)(modf(symmetricsCoords(j,0),&junk) * intScale);
            yCoords(nbAtoms)=(long)(modf(symmetricsCoords(j,1),&junk) * intScale);
            zCoords(nbAtoms)=(long)(modf(symmetricsCoords(j,2),&junk) * intScale);
            nbAtoms++;
            *pIndex++ = i;
         }
         xCoords2(i)=(long) (mScattCompList(i).mX*intScale);
         yCoords2(i)=(long) (mScattCompList(i).mY*intScale);
         zCoords2(i)=(long) (mScattCompList(i).mZ*intScale);
      }
      VFN_DEBUG_MESSAGE("-->Size of DistTable:"<<nbComponent<<" x "<<nbAtoms,2)
      //:TODO: will be needed when special positions will be implemented
      xCoords.resizeAndPreserve(nbAtoms);
      yCoords.resizeAndPreserve(nbAtoms);
      zCoords.resizeAndPreserve(nbAtoms);
      mDistTableIndex.resizeAndPreserve(nbAtoms);
      
      //cout << xCoords<<endl<<endl<<yCoords<<endl<<zCoords<<endl<<mDistTableIndex<<endl;
      
      TAU_PROFILE_STOP (timer1);
      TAU_PROFILE_START(timer2);
      // Only keep atoms which are in the asymmetric unit, or less than 3 Angstroem
      // from the asymmetric unit Cell
      if( mSpaceGroup.GetAsymUnit().Xmax()
         *mSpaceGroup.GetAsymUnit().Ymax()
         *mSpaceGroup.GetAsymUnit().Zmax() < 0.6) 
      {
         // Get rid of atoms outside asymmetric unit+3 Angstroems
         // Make sure that the list of independant atoms contains only atoms
         // inside the Asymmetric unit (else shit happens).
         
         // Prepare
            //CrystVector_long xCoordsOld,yCoordsOld,zCoordsOld,indexOld;
            //xCoordsOld=xCoords;
            //yCoordsOld=yCoords;// i>nbAtoms, so useless
            //zCoordsOld=zCoords;
            CrystVector_long indexOld;
            indexOld=mDistTableIndex;
            const long nbAtomsOld=nbAtoms;
            
            const long xMax0=(long)(mSpaceGroup.GetAsymUnit().Xmax()*intScale);
            const long yMax0=(long)(mSpaceGroup.GetAsymUnit().Ymax()*intScale);
            const long zMax0=(long)(mSpaceGroup.GetAsymUnit().Zmax()*intScale);
            const long xMax=(long)((mSpaceGroup.GetAsymUnit().Xmax()
                                    +asymUnitMargin/GetLatticePar(0))*intScale);
            const long yMax=(long)((mSpaceGroup.GetAsymUnit().Ymax()
                                    +asymUnitMargin/GetLatticePar(1))*intScale);
            const long zMax=(long)((mSpaceGroup.GetAsymUnit().Zmax()+
                                    asymUnitMargin/GetLatticePar(2))*intScale);
            // xMin < x < xMax ...but xMin is <0 so it's actually wrapped near 1...
            // So a valid point is either 0 < x < xMax or 1+xMin < x < 1
            const long xMin=(long)((1-asymUnitMargin/GetLatticePar(0))*intScale);
            const long yMin=(long)((1-asymUnitMargin/GetLatticePar(1))*intScale);
            const long zMin=(long)((1-asymUnitMargin/GetLatticePar(2))*intScale);
            VFN_DEBUG_MESSAGE("-> "<< xMax <<" "<< xMax0 <<" "<< xMin,0)
            VFN_DEBUG_MESSAGE("-> "<< yMax <<" "<< yMax0 <<" "<< yMin,0)
            VFN_DEBUG_MESSAGE("-> "<< zMax <<" "<< zMax0 <<" "<< zMin,0)
         // Do it
         nbAtoms=0;
         mDistTableIndex=0;
         for(long i=0;i<nbAtomsOld;i++)
         {
            if(  ((xCoords(i) < xMax) || (xCoords(i) > xMin))
               &&((yCoords(i) < yMax) || (yCoords(i) > yMin))
               &&((zCoords(i) < zMax) || (zCoords(i) > zMin)))
            {
               xCoords(nbAtoms)=xCoords(i);
               yCoords(nbAtoms)=yCoords(i);
               zCoords(nbAtoms)=zCoords(i);
               mDistTableIndex(nbAtoms)=indexOld(i);
               //Is it strictly inside the Asymmetric Unit ?
               //=> use as new origin
               if((xCoords(nbAtoms) <= xMax0) && 
                  (yCoords(nbAtoms) <= yMax0) && 
                  (zCoords(nbAtoms) <= zMax0))
               {
                  xCoords2(mDistTableIndex(nbAtoms))=xCoords(nbAtoms);
                  yCoords2(mDistTableIndex(nbAtoms))=yCoords(nbAtoms);
                  zCoords2(mDistTableIndex(nbAtoms))=zCoords(nbAtoms);
               }
               nbAtoms++;
            }
         }
         VFN_DEBUG_MESSAGE("-> Keeping "<<nbAtoms<<" out of "<< nbAtomsOld <<" atoms for distance calculations",2)
      }
      xCoords.resizeAndPreserve(nbAtoms);
      yCoords.resizeAndPreserve(nbAtoms);
      zCoords.resizeAndPreserve(nbAtoms);
      mDistTableIndex.resizeAndPreserve(nbAtoms);
      TAU_PROFILE_STOP (timer2);
      TAU_PROFILE_START(timer3);
      
      mDistTableIndex.resizeAndPreserve(nbAtoms);

      // Now get all distances, in fractionnal coordinates
      CrystMatrix_long xDist(nbAtoms,nbComponent);
      CrystMatrix_long yDist(nbAtoms,nbComponent);
      CrystMatrix_long zDist(nbAtoms,nbComponent);
      mDistTableSq.resize(nbAtoms,nbComponent);
      {
         const long *f1;
         const long *p1;
         long *p2;
         p2=xDist.data();
         f1=xCoords.data();
         for(long i=0;i<nbAtoms;i++)
         {
            p1=xCoords2.data();
            for(long j=0;j<nbComponent;j++)
            {
               *p2 = (*f1 - *p1++ ) & 16383;
               if(*p2>intScale2) *p2 -= intScale;
               p2++;
            }
            f1++;
         }
         p2=yDist.data();
         f1=yCoords.data();
         for(long i=0;i<nbAtoms;i++)
         {
            p1=yCoords2.data();
            for(long j=0;j<nbComponent;j++)
            {
               *p2 = (*f1 - *p1++ ) & 16383;
               if(*p2>intScale2) *p2 -= intScale;
               p2++;
            }
            f1++;
         }
         p2=zDist.data();
         f1=zCoords.data();
         for(long i=0;i<nbAtoms;i++)
         {
            p1=zCoords2.data();
            for(long j=0;j<nbComponent;j++)
            {
               *p2 = (*f1 - *p1++ ) & 16383;
               if(*p2>intScale2) *p2 -= intScale;
               p2++;
            }
            f1++;
         }
      }
      TAU_PROFILE_STOP (timer3);
      //mDistTableSq=0.;
      //convert these distances to euclidian coordinates
      //:NOTE: Assume mOrthMatrix is upper triangular
      TAU_PROFILE_START(timer4);
      {
         const double M00=mOrthMatrix(0,0);
         const double M01=mOrthMatrix(0,1);
         const double M02=mOrthMatrix(0,2);
         //const double M10=mOrthMatrix(1,0);
         const double M11=mOrthMatrix(1,1);
         const double M12=mOrthMatrix(1,2);
         //const double M20=mOrthMatrix(2,0);
         //const double M21=mOrthMatrix(2,1);
         const double M22=mOrthMatrix(2,2);
         const long *x=xDist.data();
         const long *y=yDist.data();
         const long *z=zDist.data();
         double *d=mDistTableSq.data();
         double xo,yo,zo;
         const double back= pow((double)intScale,2);
         for(long i=0;i<nbAtoms;i++)
         {
            for(long j=0;j<nbComponent;j++)
            {
               xo=M00* *x++ +M01* *y +M02* *z;
               yo=/* M10* *x + */ M11* *y++ +M12* *z;
               zo=/* M20* *x +M21* *y + */ M22* *z++;
               //mDistTableSq(i,j)=(xo*xo + yo*yo + zo*zo)/back;
               *d++ = (xo*xo + yo*yo + zo*zo)/back;
            }
         }
      }
      TAU_PROFILE_STOP (timer4);
   }
   else
   {
      VFN_DEBUG_MESSAGE("Crystal::CalcDistTable(fast=false)",3)
      TAU_PROFILE("Crystal::CalcDistTable(fast=false)","Matrix (string&)",TAU_DEFAULT);
      const int nbSymmetrics=mSpaceGroup.GetNbSymmetrics();
      CrystVector_double xCoords(nbComponent*nbSymmetrics);
      CrystVector_double yCoords(nbComponent*nbSymmetrics);
      CrystVector_double zCoords(nbComponent*nbSymmetrics);
      long nbAtoms=0;
      CrystMatrix_double symmetricsCoords(nbSymmetrics,3);
      mDistTableIndex.resize(nbComponent*nbSymmetrics);
      long *pIndex=mDistTableIndex.data();
      //get the reduced coordinates of all atoms
      for(long i=0;i<nbComponent;i++)
      {
         symmetricsCoords=mSpaceGroup.GetAllSymmetrics(mScattCompList(i).mX,
                                                     mScattCompList(i).mY,
                                                     mScattCompList(i).mZ);
         for(int j=0;j<nbSymmetrics;j++)
         {
            xCoords(nbAtoms)=symmetricsCoords(j,0);
            yCoords(nbAtoms)=symmetricsCoords(j,1);
            zCoords(nbAtoms)=symmetricsCoords(j,2);
            nbAtoms++;
            *pIndex++ = i;
         }
      }
      VFN_DEBUG_MESSAGE("-->Size of DistTable:"<<nbAtoms*nbAtoms,2)

      //:TODO: will be needed when special positions will be implemented
      xCoords.resizeAndPreserve(nbAtoms);
      yCoords.resizeAndPreserve(nbAtoms);
      zCoords.resizeAndPreserve(nbAtoms);
      mDistTableIndex.resizeAndPreserve(nbAtoms);

      // Now get all distances, in reduced coordinates
      CrystMatrix_double xDist(nbAtoms,nbAtoms);
      CrystMatrix_double yDist(nbAtoms,nbAtoms);
      CrystMatrix_double zDist(nbAtoms,nbAtoms);
      mDistTableSq.resize(nbAtoms,nbAtoms);
      /*
      {
         for(long i=0;i<nbAtoms;i++)
         {
            for(long j=0;j<i;j++)
            {  // +1000 ensures that the value is >0
               xDist(i,j)=xCoords(i)-xCoords(j);
               yDist(i,j)=yCoords(i)-yCoords(j);
               zDist(i,j)=zCoords(i)-zCoords(j);
            }
         }
      }
      */
      {
         double f1;
         const double *p1;
         double *p2;
         double junk;
         p2=xDist.data();
         for(long i=0;i<nbAtoms;i++)
         {
            f1=xCoords(i)+1000.;//:KLUDGE: +1000 ensures that the value is >0
            p1=xCoords.data();
            for(long j=0;j<i;j++)
            {
               *p2 = modf(f1 - *p1++,&junk);
               if(*p2>0.5) *p2 -= 1.;
               p2++;
            }
            p2 += nbAtoms-i;
         }
         p2=yDist.data();
         for(long i=0;i<nbAtoms;i++)
         {
            f1=yCoords(i)+1000.;
            p1=yCoords.data();
            for(long j=0;j<i;j++)
            {
               *p2 = modf(f1 - *p1++,&junk);
               if(*p2>0.5) *p2 -= 1.;
               p2++;
            }
            p2 += nbAtoms-i;
         }
         p2=zDist.data();
         for(long i=0;i<nbAtoms;i++)
         {
            f1=zCoords(i)+1000.;
            p1=zCoords.data();
            for(long j=0;j<i;j++)
            {
               *p2 = modf(f1 - *p1++,&junk);
               if(*p2>0.5) *p2 -= 1.;
               p2++;
            }
            p2 += nbAtoms-i;
         }
      }
      //convert these distances to euclidian coordinates
      {//:KLUDGE: ? Assume mOrthMatrix is upper triangular for optimization
         const double M00=mOrthMatrix(0,0);
         const double M01=mOrthMatrix(0,1);
         const double M02=mOrthMatrix(0,2);
         //const double M10=mOrthMatrix(1,0);
         const double M11=mOrthMatrix(1,1);
         const double M12=mOrthMatrix(1,2);
         //const double M20=mOrthMatrix(2,0);
         //const double M21=mOrthMatrix(2,1);
         const double M22=mOrthMatrix(2,2);
         const double *x=xDist.data();
         const double *y=yDist.data();
         const double *z=zDist.data();
         double *d=mDistTableSq.data();
         double xo,yo,zo;
         for(long i=0;i<nbAtoms;i++)
         {
            for(long j=0;j<i;j++)
            {
               xo=M00* *x +M01* *y +M02* *z;
               yo=/* M10* *x + */ M11* *y +M12* *z;
               zo=/* M20* *x +M21* *y + */ M22* *z;
               *d++ = (xo*xo+ yo*yo + zo*zo);
               x++ ; y++ ;z++;
            }
            x += nbAtoms-i;
            y += nbAtoms-i;
            z += nbAtoms-i;
            d += nbAtoms-i;
            mDistTableSq(i,i)=0;
         }
      }
      {
         const double *d=mDistTableSq.data();
         double *d2;
         for(long i=0;i<nbAtoms;i++)
         {
            d2=mDistTableSq.data()+i;
            for(long j=0;j<i;j++)
            {
               *d2 = *d++;
               d2 += nbAtoms;
            }
            d += nbAtoms-i;
         }
      }
   }
   mClockNeighborTable.Click();
   VFN_DEBUG_MESSAGE("Crystal::CalcDistTable():End",3)
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
