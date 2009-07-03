/*  ObjCryst++ Object-Oriented Crystallographic Library
    (c) 2006- Vincent Favre-Nicolin vincefn@users.sourceforge.net

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
*  source file for Indexing classes & functions
*
*/
#include <algorithm>

#include "ObjCryst/Indexing.h"
#include "Quirks/VFNDebug.h"
#include "Quirks/VFNStreamFormat.h"
#include "Quirks/Chronometer.h"

using namespace std;

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

#ifndef DEG2RAD
#define DEG2RAD (M_PI/180.)
#endif
#ifndef RAD2DEG
#define RAD2DEG (180./M_PI)
#endif

// Alternate parametrization of triclinic cells
#define TRI1

namespace ObjCryst
{

float EstimateCellVolume(const float dmin, const float dmax, const float nbrefl, 
                         const CrystalSystem system,const CrystalCentering centering,const float kappa)
{
   const float q1=dmin*dmin*dmin-dmax*dmax*dmax;
   const float q2=dmin*dmin     -dmax*dmax;
   float D0,C0;
   if(system==TRICLINIC)
   {
      C0=2.095;
      return nbrefl/(C0*kappa*q1);
   }
   if(system==CUBIC)
   {
      if(centering==LATTICE_P) D0=0.862;
      if(centering==LATTICE_I) D0=0.475;
      if(centering==LATTICE_F) D0=0.354;
      return pow(nbrefl/(D0*kappa*q2),(float)1.5);
   }
   // "*.85" means using D0_min rather than D0
   if(system==MONOCLINIC) {C0=1.047;D0=0.786*.85;}
   if(system==ORTHOROMBIC){C0=0.524;D0=1.36 *.85;}
   if(system==HEXAGONAL)  {C0=0.150;D0=1.04 *.85;}
   if(system==RHOMBOEDRAL){C0=0.230;D0=1.04 *.85;}
   if(system==TETRAGONAL) {C0=0.214;D0=1.25 *.85;}
   if((centering==LATTICE_I)||(centering==LATTICE_A)||(centering==LATTICE_B)||(centering==LATTICE_C)) {C0/=2;D0/=2;}
   if(centering==LATTICE_F){C0/=4;D0/=4;}
   const double alpha=D0*q2/(3*C0*q1), beta=nbrefl/(2*kappa*C0*q1);
   const double eta=beta-alpha*alpha*alpha,gamma=sqrt(beta*beta-2*beta*alpha*alpha*alpha);
   const float v=pow(pow(eta+gamma,(double)1/3.)+pow((eta-gamma),(double)1/3.)-alpha,(double)3);
   return v;
}

/** light-weight class storing the reciprocal space unitcell
*/

RecUnitCell::RecUnitCell(const float zero,const float p0,const float p1,const float p2,
                         const float p3,const float p4,const float p5,CrystalSystem lattice,
                         const CrystalCentering cent):
mlattice(lattice),mCentering(cent)
{
   par[0]=zero;
   par[1]=p0;
   par[2]=p1;
   par[3]=p2;
   par[4]=p3;
   par[5]=p4;
   par[6]=p5;
}

RecUnitCell::RecUnitCell(const RecUnitCell &old)
{
   *this=old;
}

void RecUnitCell::operator=(const RecUnitCell &rhs)
{
   for(unsigned int i=0;i<7;++i) par[i]=rhs.par[i];
   mlattice=rhs.mlattice;
   mCentering=rhs.mCentering;
}

float RecUnitCell::hkl2d(const float h,const float k,const float l,REAL *derivpar,const unsigned int derivhkl) const
{
   if((derivpar==NULL)&&(derivhkl==0))
   {
      switch(mlattice)
      {
         case TRICLINIC:
         {
            #ifdef TRI1
            return par[0]+par[1]*par[1]*h*h + par[2]*par[2]*k*k + par[3]*par[3]*l*l + par[4]*h*k + par[5]*k*l + par[6]*h*l;
            #else
            return par[0]+par[1]*par[1]*h*h + par[2]*par[2]*k*k + par[3]*par[3]*l*l
                  + 2*par[1]*par[2]*par[4]*h*k + 2*par[2]*par[3]*par[5]*k*l + 2*par[1]*par[3]*par[6]*h*l;
            #endif
            break;
         }
         case MONOCLINIC:
         {
            return par[0]+par[1]*par[1]*h*h + par[2]*par[2]*k*k + par[3]*par[3]*l*l + 2*par[1]*par[3]*par[4]*h*l;
            break;
         }
         case ORTHOROMBIC:
         {
            return par[0]+par[1]*par[1]*h*h + par[2]*par[2]*k*k + par[3]*par[3]*l*l;
            break;
         }
         case HEXAGONAL:
         {
            return par[0]+par[1]*par[1]*(h*h + k*k + h*k)+ par[2]*par[2]*l*l ;
            break;
         }
         case RHOMBOEDRAL:
         {
            return par[0]+par[1]*par[1]*(h*h + k*k + l*l + 2*par[2]*(h*k + k*l + h*l));
            break;
         }
         case TETRAGONAL:
         {
            return par[0]+par[1]*par[1]*(h*h + k*k) + par[2]*par[2]*l*l;
            break;
         }
         case CUBIC:
         {
            return par[0]+par[1]*par[1]*(h*h+k*k+l*l);
            break;
         }
      }
   }
   if(derivhkl==1)
   {
      switch(mlattice)
      {
         case TRICLINIC:
         {
            #ifdef TRI1
            return 2*par[1]*par[1]*h + par[4]*k + par[6]*l;
            #else
            return 2*par[1]*par[1]*h + 2*par[1]*par[2]*par[4]*k + 2*par[1]*par[3]*par[6]*l;
            #endif
            break;
         }
         case MONOCLINIC:
         {
            return 2*par[1]*par[1]*h + 2*par[1]*par[3]*par[4]*l;
            break;
         }
         case ORTHOROMBIC:
         {
            return 2*par[1]*par[1]*h;
            break;
         }
         case HEXAGONAL:
         {
            return par[1]*par[1]*(2*h + k);
            break;
         }
         case RHOMBOEDRAL:
         {
            return par[1]*par[1]*(2*h + 2*par[2]*(k + l));
            break;
         }
         case TETRAGONAL:
         {
            return 2*par[1]*par[1]*h;
            break;
         }
         case CUBIC:
         {
            return 2*par[1]*par[1]*h;
            break;
         }
      }
   }
   if(derivhkl==2)
   {
      switch(mlattice)
      {
         case TRICLINIC:
         {
            #ifdef TRI1
            return 2*par[2]*par[2]*k + par[4]*h + par[5]*l;
            #else
            return 2*par[2]*par[2]*k + 2*par[1]*par[2]*par[4]*h + 2*par[2]*par[3]*par[5]*l;
            #endif
            break;
         }
         case MONOCLINIC:
         {
            return 2*par[2]*par[2]*k;
            break;
         }
         case ORTHOROMBIC:
         {
            return 2*par[2]*par[2]*k;
            break;
         }
         case HEXAGONAL:
         {
            return par[1]*par[1]*(2*k + h);
            break;
         }
         case RHOMBOEDRAL:
         {
            return par[1]*par[1]*(2*k + l*l + 2*par[2]*(h + l));
            break;
         }
         case TETRAGONAL:
         {
            return 2*par[1]*par[1]*k;
            break;
         }
         case CUBIC:
         {
            return 2*par[1]*par[1]*k;
            break;
         }
      }
   }
   if(derivhkl==3)
   {
      switch(mlattice)
      {
         case TRICLINIC:
         {
            #ifdef TRI1
            return 2*par[3]*par[3]*l + par[5]*k + par[6]*h;
            #else
            return 2*par[3]*par[3]*l + 2*par[2]*par[3]*par[5]*k + 2*par[1]*par[3]*par[6]*h;
            #endif
            break;
         }
         case MONOCLINIC:
         {
            return 2*par[3]*par[3]*l + 2*par[1]*par[3]*par[4]*h;
            break;
         }
         case ORTHOROMBIC:
         {
            return 2*par[3]*par[3]*l;
            break;
         }
         case HEXAGONAL:
         {
            return 2*par[2]*par[2]*l;
            break;
         }
         case RHOMBOEDRAL:
         {
            return par[1]*par[1]*(2*l + 2*par[2]*(k + h));
            break;
         }
         case TETRAGONAL:
         {
            return 2*par[2]*par[2]*l;
            break;
         }
         case CUBIC:
         {
            return 2*par[1]*par[1]*l;
            break;
         }
      }
   }

   if(derivpar==&par[0]) return 1.0;
   if(derivpar==(par+1))
   {
      switch(mlattice)
      {
         case TRICLINIC:
         {
            #ifdef TRI1
            return 2*par[1]*h*h;
            #else
            return 2*par[1]*h*h + 2*par[2]*par[4]*h*k + 2*par[3]*par[6]*h*l;
            #endif
            break;
         }
         case MONOCLINIC:
         {
            return 2*par[1]*h*h + 2*par[3]*par[4]*h*l;
            break;
         }
         case ORTHOROMBIC:
         {
            return 2*par[1]*h*h;
            break;
         }
         case HEXAGONAL:
         {
            return 2*par[1]*(h*h + k*k + h*k);
            break;
         }
         case RHOMBOEDRAL:
         {
            return 2*par[1]*(h*h + k*k + l*l + 2*par[2]*(h*k + k*l + h*l));
            break;
         }
         case TETRAGONAL:
         {
            return 2*par[1]*(h*h + k*k);
            break;
         }
         case CUBIC:
         {
            return 2*par[1]*(h*h+k*k+l*l);
            break;
         }
      }
   }
   if(derivpar==(par+2))
   {
      switch(mlattice)
      {
         case TRICLINIC:
         {
            #ifdef TRI1
            return 2*par[2]*k*k;
            #else
            return 2*par[2]*k*k + 2*par[1]*par[4]*h*k + 2*par[2]*par[5]*k*l;
            #endif
            break;
         }
         case MONOCLINIC:
         {
            return 2*par[2]*k*k;
            break;
         }
         case ORTHOROMBIC:
         {
            return 2*par[2]*k*k;
            break;
         }
         case HEXAGONAL:
         {
            return 2*par[2]*l*l ;
            break;
         }
         case RHOMBOEDRAL:
         {
            return par[1]*par[1]*(h*h + k*k + l*l + 2*(h*k + k*l + h*l));
            break;
         }
         case TETRAGONAL:
         {
            return 2*par[2]*l*l;
            break;
         }
         case CUBIC:
         {
            throw 0;
            break;
         }
      }
   }
   if(derivpar==(par+3))
   {
      switch(mlattice)
      {
         case TRICLINIC:
         {
            #ifdef TRI1
            return 2*par[3]*l*l;
            #else
            return 2*par[3]*l*l + 2*par[2]*par[5]*k*l + 2*par[1]*par[6]*h*l;
            #endif
            break;
         }
         case MONOCLINIC:
         {
            return 2*par[3]*l*l + 2*par[1]*par[4]*h*l;
            break;
         }
         case ORTHOROMBIC:
         {
            return 2*par[3]*l*l;
            break;
         }
         case HEXAGONAL:
         {
            throw 0;
            break;
         }
         case RHOMBOEDRAL:
         {
            throw 0;
            break;
         }
         case TETRAGONAL:
         {
            throw 0;
            break;
         }
         case CUBIC:
         {
            throw 0;
            break;
         }
      }
   }
   if(derivpar==(par+4))
   {
      switch(mlattice)
      {
         case TRICLINIC:
         {
            #ifdef TRI1
            return h*k;
            #else
            return 2*par[1]*par[2]*h*k;
            #endif
            break;
         }
         case MONOCLINIC:
         {
            return 2*par[1]*par[3]*h*l;
            break;
         }
         default:
         {
            throw 0;
            break;
         }
      }
   }
   if(derivpar==(par+5))
   {
      switch(mlattice)
      {
         case TRICLINIC:
         {
            #ifdef TRI1
            return k*l;
            #else
            return 2*par[2]*par[3]*k*l;
            #endif
            break;
         }
         default:
         {
            throw 0;
            break;
         }
      }
   }
   if(derivpar==(par+6))
   {
      switch(mlattice)
      {
         case TRICLINIC:
         {
            #ifdef TRI1
            return h*l;
            #else
            return 2*par[1]*par[3]*h*l;
            #endif
            break;
         }
         default:
         {
            throw 0;
            break;
         }
      }
   }
   throw 0;
   return 0.0;
}

void RecUnitCell::hkl2d_delta(const float h,const float k,const float l,
                              const RecUnitCell &delta, float & dmin, float &dmax) const
{
   const float p0m=par[0]-delta.par[0] , p0p=par[0]+delta.par[0],
               p1m=par[1]-delta.par[1] , p1p=par[1]+delta.par[1],
               p2m=par[2]-delta.par[2] , p2p=par[2]+delta.par[2],
               p3m=par[3]-delta.par[3] , p3p=par[3]+delta.par[3],
               p4m=par[4]-delta.par[4] , p4p=par[4]+delta.par[4],
               p5m=par[5]-delta.par[5] , p5p=par[5]+delta.par[5],
               p6m=par[6]-delta.par[6] , p6p=par[6]+delta.par[6];
   switch(mlattice)
   {
      case TRICLINIC:
      {//return par[0]+par[1]*par[1]*h*h + par[2]*par[2]*k*k + par[3]*par[3]*l*l + par[4]*h*k + par[5]*k*l + par[6]*h*l;
         #ifdef TRI1
         float p4mm,p5mm,p6mm,p4pp,p5pp,p6pp;
         if((h*k)>0){p4mm=p4m;p4pp=p4p;}else{p4mm=p4p;p4pp=p4m;}
         if((k*l)>0){p5mm=p5m;p5pp=p5p;}else{p5mm=p5p;p5pp=p5m;}
         if((h*l)>0){p6mm=p6m;p6pp=p6p;}else{p6mm=p6p;p6pp=p6m;}
         dmin=p0m+p1m*p1m*h*h+p2m*p2m*k*k+p3m*p3m*l*l+p4mm*h*k+p5mm*k*l+p6mm*h*l;
         dmax=p0p+p1p*p1p*h*h+p2p*p2p*k*k+p3p*p3p*l*l+p4pp*h*k+p5pp*k*l+p6pp*h*l;
         #else
         float p1mm,p2mm,p3mm,p4mm,p5mm,p6mm,p1pp,p2pp,p3pp,p4pp,p5pp,p6pp;
         if((h*(par[1]*h+par[2]*par[4]*k+par[3]*par[6]*l))>0) {p1mm=p1m;p1pp=p1p;}else{p1mm=p1p;p1pp=p1m;}
         if((k*(par[2]*k+par[1]*par[4]*h+par[3]*par[5]*l))>0) {p2mm=p2m;p2pp=p2p;}else{p2mm=p2p;p2pp=p2m;}
         if((l*(par[3]*l+par[2]*par[5]*k+par[1]*par[6]*h))>0) {p3mm=p3m;p3pp=p3p;}else{p3mm=p3p;p3pp=p3m;}
         if((h*k)>0){p4mm=p4m;p4pp=p4p;}else{p4mm=p4p;p4pp=p4m;}
         if((k*l)>0){p5mm=p5m;p5pp=p5p;}else{p5mm=p5p;p5pp=p5m;}
         if((h*l)>0){p6mm=p6m;p6pp=p6p;}else{p6mm=p6p;p6pp=p6m;}
         dmin=p0m+p1mm*p1mm*h*h+p2mm*p2mm*k*k+p3mm*p3mm*l*l+2*p1mm*p2mm*p4mm*h*k+2*p2mm*p3mm*p5mm*k*l+2*p1mm*p3mm*p6mm*h*l;
         dmax=p0p+p1pp*p1pp*h*h+p2pp*p2pp*k*k+p3pp*p3pp*l*l+2*p1pp*p2pp*p4pp*h*k+2*p2pp*p3pp*p5pp*k*l+2*p1pp*p3pp*p6pp*h*l;
         #endif
         return;
      }
      case MONOCLINIC: //OK
      {
         if((h*l)>0)
         {
            dmin = p0m + p1m*p1m*h*h + p2m*p2m*k*k + p3m*p3m*l*l + 2*p1m*p3m*p4m*h*l;
            dmax = p0p + p1p*p1p*h*h + p2p*p2p*k*k + p3p*p3p*l*l + 2*p1p*p3p*p4p*h*l;
            return;
         }
         else
         {
            const bool b1=(h*(par[1]*h+par[3]*par[4]*l))>0;// d(d*^2)/dp1
            const bool b3=(l*(par[3]*l+par[1]*par[4]*h))>0;// d(d*^2)/dp2
            if(b1 && b3)
            {
               dmin = p0m + p1m*p1m*h*h + p2m*p2m*k*k + p3m*p3m*l*l + 2*p1m*p3m*p4p*h*l;
               dmax = p0p + p1p*p1p*h*h + p2p*p2p*k*k + p3p*p3p*l*l + 2*p1p*p3p*p4m*h*l;
               return;
            }
            else if(b1 && (!b3))
               {
                  dmin = p0m + p1m*p1m*h*h + p2m*p2m*k*k + p3p*p3p*l*l + 2*p1m*p3p*p4p*h*l;
                  dmax = p0p + p1p*p1p*h*h + p2p*p2p*k*k + p3m*p3m*l*l + 2*p1p*p3m*p4m*h*l;
                  return;
               }
               else if((!b1) && b3)
                  {
                     dmin = p0m + p1p*p1p*h*h + p2m*p2m*k*k + p3m*p3m*l*l + 2*p1p*p3m*p4p*h*l;
                     dmax = p0p + p1m*p1m*h*h + p2p*p2p*k*k + p3p*p3p*l*l + 2*p1m*p3p*p4m*h*l;
                     return;
                  }
                  else
                  {
                     dmin = p0m + p1p*p1p*h*h + p2m*p2m*k*k + p3p*p3p*l*l + 2*p1p*p3p*p4p*h*l;
                     dmax = p0p + p1m*p1m*h*h + p2p*p2p*k*k + p3m*p3m*l*l + 2*p1m*p3m*p4m*h*l;
                     return;
                  }
            }
      }
      case ORTHOROMBIC: //OK
      {
         dmin= p0m + p1m*p1m*h*h + p2m*p2m*k*k + p3m*p3m*l*l;
         dmax= p0p + p1p*p1p*h*h + p2p*p2p*k*k + p3p*p3p*l*l;
         return;
      }
      case HEXAGONAL: //OK
      {
         dmin=p0m+p1m*p1m*(h*h + k*k + h*k)+ p2m*p2m*l*l ;
         dmax=p0p+p1p*p1p*(h*h + k*k + h*k)+ p2p*p2p*l*l ;
         return;
      }
      case RHOMBOEDRAL:
      {
         if((h*k + k*l + h*l)>=0)
         {
            dmin= p0m+p1m*p1m*(h*h + k*k + l*l + 2*p2m*(h*k + k*l + h*l));
            dmax= p0p+p1p*p1p*(h*h + k*k + l*l + 2*p2p*(h*k + k*l + h*l));
         }
         else
         {
            dmin= p0m+p1m*p1m*(h*h + k*k + l*l + 2*p2p*(h*k + k*l + h*l));
            dmax= p0p+p1p*p1p*(h*h + k*k + l*l + 2*p2m*(h*k + k*l + h*l));
         }
         return;
      }
      case TETRAGONAL: //OK
      {
         dmin= p0m + p1m*p1m*(h*h + k*k) + p2m*p2m*l*l;
         dmax= p0p + p1p*p1p*(h*h + k*k) + p2p*p2p*l*l;
         return;
      }
      case CUBIC: //OK
      {
         dmin=p0m + p1m*p1m*(h*h+k*k+l*l);
         dmax=p0p + p1p*p1p*(h*h+k*k+l*l);
         return;
      }
   }
}

vector<float> RecUnitCell::DirectUnitCell(const bool equiv)const
{
   // reciprocal unit cell parameters
   float aa,bb,cc,calphaa,cbetaa,cgammaa,salphaa,sbetaa,sgammaa;
   switch(mlattice)
   {
      case TRICLINIC:
      {
         #ifdef TRI1
         aa=par[1];
         bb=par[2];
         cc=par[3];
         calphaa=par[5]/(2*bb*cc);
         cbetaa=par[6]/(2*aa*cc);
         cgammaa=par[4]/(2*aa*bb);
         salphaa=sqrt(abs(1-calphaa*calphaa));
         sbetaa =sqrt(abs(1-cbetaa *cbetaa));
         sgammaa=sqrt(abs(1-cgammaa*cgammaa));
         #else
         aa=par[1];
         bb=par[2];
         cc=par[3];
         calphaa=par[5];
         cbetaa=par[6];
         cgammaa=par[4];
         salphaa=sqrt(abs(1-calphaa*calphaa));
         sbetaa =sqrt(abs(1-cbetaa *cbetaa));
         sgammaa=sqrt(abs(1-cgammaa*cgammaa));
         #endif
         break;
      }
      case MONOCLINIC:
      {
         aa=par[1];
         bb=par[2];
         cc=par[3];
         calphaa=0;
         cbetaa=par[4];
         cgammaa=0;
         salphaa=1;
         sbetaa =sqrt(abs(1-cbetaa *cbetaa));
         sgammaa=1;
         break;
      }
      case ORTHOROMBIC:
      {
         aa=par[1];
         bb=par[2];
         cc=par[3];
         calphaa=0;
         cbetaa =0;
         cgammaa=0;
         salphaa=1;
         sbetaa =1;
         sgammaa=1;
         break;
      }
      case HEXAGONAL:
      {
         aa=par[1];
         bb=par[1];
         cc=par[2];
         calphaa=0;
         cbetaa =0;
         cgammaa=0.5;
         salphaa=1;
         sbetaa =1;
         sgammaa=0.8660254037844386;
         break;
      }
      case RHOMBOEDRAL:
      {
         aa=par[1];
         bb=par[1];
         cc=par[1];
         calphaa=par[4];
         cbetaa =par[4];
         cgammaa=par[4];
         salphaa=sqrt(abs(1-calphaa *calphaa));
         sbetaa =salphaa;
         sgammaa=salphaa;
         break;
      }
      case TETRAGONAL:
      {
         aa=par[1];
         bb=par[1];
         cc=par[2];
         calphaa=0;
         cbetaa =0;
         cgammaa=0;
         salphaa=1;
         sbetaa =1;
         sgammaa=1;
         break;
      }
      case CUBIC:
      {
         aa=par[1];
         bb=par[1];
         cc=par[1];
         calphaa=0;
         cbetaa =0;
         cgammaa=0;
         salphaa=1;
         sbetaa =1;
         sgammaa=1;
         break;
      }
   }
   // Volume of reciprocal unit cell
   const float vv=sqrt(abs(1-calphaa*calphaa-cbetaa*cbetaa-cgammaa*cgammaa+2*calphaa*cbetaa*cgammaa));

   const float a=salphaa/(aa*vv);
   const float b=sbetaa /(bb*vv);
   const float c=sgammaa/(cc*vv);

   const float calpha=(cbetaa *cgammaa-calphaa)/(sbetaa *sgammaa);
   const float cbeta =(calphaa*cgammaa-cbetaa )/(salphaa*sgammaa);
   const float cgamma=(calphaa*cbetaa -cgammaa)/(salphaa*sbetaa );
   
   const float alpha=acos( calpha );
   const float beta =acos( cbeta );
   const float gamma=acos( cgamma );

   const float v=a*b*c*sqrt(1-calpha*calpha-cbeta*cbeta-cgamma*cgamma+2*calpha*cbeta*cgamma);
   
   vector<float> uc(7);
   uc[0]=a;
   uc[1]=b;
   uc[2]=c;
   uc[3]=alpha;
   uc[4]=beta;
   uc[5]=gamma;
   uc[6]=v;
   return uc;
}
///////////////////////////////////////////////// PEAKLIST:HKL0 /////////////////////
PeakList::hkl0::hkl0(const int h0,const int k0, const int l0):
h(h0),k(k0),l(l0)
{}

///////////////////////////////////////////////// PEAKLIST:HKL /////////////////////
PeakList::hkl::hkl(const float d,const float i,const float ds,const float is,
                   const int h0,const int k0, const int l0,const float dc0):
dobs(d),dobssigma(ds),d2obs(d*d),d2obsmin((d-ds/2)*(d-ds/2)),d2obsmax((d+ds/2)*(d+ds/2)),iobs(i),iobssigma(is),
h(h0),k(k0),l(l0),isIndexed(false),isSpurious(false),stats(0),
d2calc(dc0),d2diff(0)
{}

bool compareHKL_d(const PeakList::hkl &d1, const PeakList::hkl &d2)
{
   return d1.dobs < d2.dobs;
}


///////////////////////////////////////////////// PEAKLIST /////////////////////

PeakList::PeakList()
{}

PeakList::PeakList(const PeakList &old)
{
   *this=old;
}

void PeakList::operator=(const PeakList &rhs)
{
   VFN_DEBUG_ENTRY("PeakList::operator=(PeakList &old)",10);
   mvHKL=rhs.mvHKL;
   VFN_DEBUG_EXIT("PeakList::operator=(PeakList &old)",10);
}

PeakList::~PeakList()
{}

void PeakList::ImportDhklDSigmaIntensity(istream &is,float defaultsigma)
{
   float d,sigma,iobs;
   while(true)
   {// :TODO: use readline to make sure when the end is reached
      is >>d;
      if(is.eof()) break;
      is>>sigma;
      if(is.eof()) break;
      is>>iobs;
      if(sigma<=0) sigma=d*defaultsigma;
      if(iobs<=0) iobs=1.0;
      mvHKL.push_back(hkl(1/d,iobs,1/(d-sigma/2)-1/(d+sigma/2)));
      cout<<__FILE__<<":"<<__LINE__<<"  "<<mvHKL.size()<<":d="<<d<<"+/-"<<sigma<<", I="<<iobs<<" 1/d="<<1/d<<endl;
      if(is.eof()) break;
   }
   sort(mvHKL.begin(),mvHKL.end(),compareHKL_d);
   cout<<"Imported "<<mvHKL.size()<<" observed reflection positions."<<endl;
}

void PeakList::ImportDhklIntensity(istream &is)
{
   float d,iobs;
   while(true)
   {// :TODO: use readline to make sure when the end is reached
      is >>d;
      if(is.eof()) break;
      is>>iobs;
      mvHKL.push_back(hkl(1/d,iobs));
      cout<<__FILE__<<":"<<__LINE__<<"  "<<mvHKL.size()<<":d="<<d<<", I="<<iobs<<" 1/d="<<1/d<<endl;
      if(is.eof()) break;
   }
   sort(mvHKL.begin(),mvHKL.end(),compareHKL_d);
   cout<<"Imported "<<mvHKL.size()<<" observed reflection positions."<<endl;
}

void PeakList::ImportDhkl(istream &is)
{
   std::vector<std::pair<float,float> > v;
   float d;
   while(true)
   {// :TODO: use readline to make sure when the end is reached
      is >>d;
      if(is.eof()) break;
      mvHKL.push_back(hkl(1/d));
      cout<<__FILE__<<":"<<__LINE__<<"  "<<mvHKL.size()<<":d="<<d<<" 1/d="<<1/d<<endl;
      if(is.eof()) break;
   }
   sort(mvHKL.begin(),mvHKL.end(),compareHKL_d);
   cout<<"Imported "<<v.size()<<" observed reflection positions."<<endl;
}


template<class T,class U> bool comparePairFirst(std::pair<T,U> &p1, std::pair<T,U> &p2)
{
   return p1.first < p2.first;
}

void PeakList::Import2ThetaIntensity(istream &is, const float wavelength)
{
   std::list<std::pair<float,float> > v;
   float d,iobs;
   while(true)
   {// :TODO: use readline to make sure when the end is reached
      is >>d;
      if(is.eof()) break;
      is>>iobs;
      d=2*sin(d/2*DEG2RAD)/wavelength;
      mvHKL.push_back(hkl(1/d,iobs));
      cout<<__FILE__<<":"<<__LINE__<<"  "<<mvHKL.size()<<":d="<<1/d<<", I="<<iobs<<" 1/d="<<d<<endl;
      if((is.eof())||v.size()>=20) break;
   }
   sort(mvHKL.begin(),mvHKL.end(),compareHKL_d);
   cout<<"Imported "<<v.size()<<" observed reflection positions."<<endl;
}

void PeakList::ExportDhklDSigmaIntensity(std::ostream &os)const
{
   char buf[100];
   for(vector<PeakList::hkl>::const_iterator pos=mvHKL.begin();pos!=mvHKL.end();++pos)
   {
      const float sigma=1/(pos->dobs-pos->dobssigma/2)-1/(pos->dobs+pos->dobssigma/2);
      sprintf(buf,"%6.3f %6.3f %f",1/pos->dobs,sigma,pos->iobs);
      os<<buf<<endl;
   }
}

void PeakList::AddPeak(const float d, const float iobs,const float dobssigma,const float iobssigma,
                       const int h,const int k, const int l,const float d2calc)
{
   if(dobssigma<=0)
   {// Manually added peak ? Use other reflection's sigmas to evaluate sigma for this reflection
      float s=0;
      for(vector<hkl>::const_iterator pos=mvHKL.begin();pos!=mvHKL.end();++pos)
         s+= pos->dobssigma;
      s/=mvHKL.size();
      mvHKL.push_back(hkl(d,iobs,s,iobssigma,h,k,l,d2calc));
   }
   else mvHKL.push_back(hkl(d,iobs,dobssigma,iobssigma,h,k,l,d2calc));
   sort(mvHKL.begin(),mvHKL.end(),compareHKL_d);
   //this->Print(cout);
}

void PeakList::RemovePeak(unsigned int idx)
{
   for(unsigned int i=idx;i<(mvHKL.size()-1);++i) mvHKL[i]=mvHKL[i+1];
   mvHKL.pop_back();
}

void PeakList::Print(std::ostream &os) const
{
   unsigned int i=0;
   char buf[200];
   os<<"PeakList, with "<<mvHKL.size()<<" peaks"<<endl;
   for(vector<PeakList::hkl>::const_iterator pos=mvHKL.begin();pos!=mvHKL.end();++pos)
   {
      const float sigma=1/(pos->dobs-pos->dobssigma/2)-1/(pos->dobs+pos->dobssigma/2);
      if(pos->isIndexed)
         sprintf(buf,"#%3d d=%6.3f+/-%7.4f dcalc=%6.3f, diff=%7.4f, iobs=%6.3f HKL=%2d %2d %2d Spurious=%1d stats=%6d",
                 i++,1/pos->dobs,sigma,
                 1/sqrt(abs(pos->d2calc)),1/sqrt(abs(pos->d2calc))-1/pos->dobs,
                 pos->iobs,pos->h,pos->k,pos->l,pos->isSpurious,pos->stats);
      else
         sprintf(buf,"#%3d d=%6.3f+/-%6.3f              iobs=%6.3f  UNINDEXED   Spurious=%1d stats=%6d",
                 i++,1/pos->dobs,1/(pos->dobs-pos->dobssigma/2)-1/(pos->dobs+pos->dobssigma/2),
                 pos->iobs,pos->isSpurious,pos->stats);
      os<<buf<<endl;
   }
}

vector<PeakList::hkl> & PeakList::GetPeakList(){return mvHKL;}

const vector<PeakList::hkl> & PeakList::GetPeakList()const {return mvHKL;}

/////////////////////////////////////////////////////// SCORE ///////////////////////////////////////

float Score(const PeakList &dhkl, const RecUnitCell &ruc, const unsigned int nbSpurious,
            const bool verbose,const bool storehkl,const bool storePredictedHKL)
{
   const bool autozero=false;
   vector<PeakList::hkl>::const_iterator pos,first,last;
   for(pos=dhkl.GetPeakList().begin();pos!=dhkl.GetPeakList().end();++pos)
   {
      if(storehkl) pos->isIndexed=false;
      pos->d2calc=0;
      pos->d2diff=1000;
   }
   const unsigned long nb=dhkl.GetPeakList().size();
   if(storePredictedHKL) dhkl.mvPredictedHKL.clear();
   
   unsigned long nbCalc=0;
   int h,k,l;
   float predict_coeff=1;
   if(storePredictedHKL)predict_coeff=2;
   const float dmax=dhkl.mvHKL[nb-1].d2obs*predict_coeff*1.05;
   int sk0,sl0;// do we need >0 *and* <0 indices for k,l ?
   switch(ruc.mlattice)
   {
      case TRICLINIC:
         sk0=-1;sl0=-1;
	 break;
      case MONOCLINIC:
         sk0=1;sl0=-1;
         break;
      case ORTHOROMBIC:
         sk0=1;sl0=1;
         break;
      case HEXAGONAL:
         sk0=-1;sl0=1;
         break;
      case RHOMBOEDRAL:
         sk0=-1;sl0=-1;
         break;
      case TETRAGONAL:
         sk0=1;sl0=1;
         break;
      case CUBIC:
         sk0=1;sl0=1;
         break;
   }
   int stepk,stepl;// steps in k,l to use for centered lattices
   switch(ruc.mCentering)
   {
      case LATTICE_P:stepk=1;stepl=1;break;
      case LATTICE_I:stepk=1;stepl=2;break;
      case LATTICE_A:stepk=1;stepl=2;break;
      case LATTICE_B:stepk=1;stepl=2;break;
      case LATTICE_C:stepk=2;stepl=1;break;
      case LATTICE_F:stepk=2;stepl=2;break;
   }
   first=dhkl.GetPeakList().begin();last=dhkl.GetPeakList().end();
   unsigned long nbCalcH,nbCalcK;// Number of calculated lines below dmax for each h,k
   for(h=0;;++h)
   {
      nbCalcH=0;
      for(int sk=sk0;sk<=1;sk+=2)
      {
         if(h==0) sk=1;// no need to explore 0kl with both sk -1 and 1
         if(stepk==2) k=(h%2);// For LATTICE_C,LATTICE_F: h odd => k odd
         else k=0;
         for(;;k+=stepk)
         {
            nbCalcK=0;
            for(int sl=sl0;sl<=1;sl+=2)
            {
               if((h+k)==0)
               {
                  sl=1;// No need to list 0 0 l with l<0
                  l=1;
               }
               else
               {
                  if(h==0)
                  {
                     if(ruc.mlattice==MONOCLINIC) sl=1;// 0 k l and 0 k -l are equivalent
                     if((sk<0)||(sl<0)) l=1;// Do not list 0 k 0 with k<0
                     else l=0;// h==k==0 already covered
                  }
                  else
                  {
                     if(sl<0) l=1;// Do not list h k 0 twice
                     else l=0;
                  }
               }
               if(stepl==2)
               {
                  if(ruc.mCentering==LATTICE_I) l+=(h+k+l)%2;
                  if(ruc.mCentering==LATTICE_A) l+=(k+l)%2;// Start at hk1 if k odd
                  if(  (ruc.mCentering==LATTICE_B)
                     ||(ruc.mCentering==LATTICE_F)) l+=(h+l)%2;// Start at hk1 if h odd
               }
               for(;;l+=stepl)
               {
                  const float d2=ruc.hkl2d(h,sk*k,sl*l);
                  if(d2>dmax)
                  {
                     //cout<<__FILE__<<":"<<__LINE__<<" hkl: "<<h<<" "<<sk*k<<" "<<sl*l<<":"<<sqrt(d2)<<" deriv="<<sl*ruc.hkl2d(h,sk*k,sl*l,NULL,3)<<"/"<<sqrt(dmax)<<endl;
                     // Only break if d is increasing with l 
                     if((sl*ruc.hkl2d(h,sk*k,sl*l,NULL,3))>=0) break;
                     else continue;
                  }
                  nbCalc++;nbCalcK++;nbCalcH++;
                  if(storePredictedHKL)
                  {
                     dhkl.mvPredictedHKL.push_back(PeakList::hkl(0,0,0,0,h,sk*k,sl*l,d2));
                     //continue;
                  }
                  for(pos=first;pos!=last;++pos)
                  {
                     const float tmp=d2-pos->d2obs;
                     if(tmp<.1)
                     {
                        if(tmp<-.1) break;
                        if(fabs(tmp)<fabs(pos->d2diff))
                        {  
                           pos->d2diff=tmp;
                           if(storehkl)
                           {
                              pos->h=h;
                              pos->k=sk*k;
                              pos->l=sl*l;
                              pos->isIndexed=true;
                              pos->d2calc=d2;
                           }
                        }
                        /*
                        if((verbose)&&(fabs(tmp)<.01))
                           cout<<__FILE__<<":"<<__LINE__<<"      hkl: "<<h<<" "<<k<<" "<<l
                              <<"#"<<i<<": calc="<<sqrt(d2)<<", obs="<<sqrt(*pd2)<<", min_epsilon="<<*pdq2<<", tmp="<<tmp<<endl;
                        */
                     }
                  }
               }
            }
            if(nbCalcK==0) //d(hk0)>dmax
            {
               //cout<<__FILE__<<":"<<__LINE__<<" hkl: "<<h<<" "<<sk*k<<" "<<0<<" deriv="<<sk*ruc.hkl2d(h,sk*k,0,NULL,2)<<endl;
               if((sk*ruc.hkl2d(h,sk*k,0,NULL,2))>=0) break;
            }
         }
      }
      if(nbCalcH==0) break;//h00 beyond limit
   }
   float epsilon=0.0,zero=0.0;
   if(autozero)
   {
      for(pos=dhkl.GetPeakList().begin();pos!=dhkl.GetPeakList().end();++pos) zero+=pos->d2diff;
      zero/=nb;
   }
   for(pos=dhkl.GetPeakList().begin();pos!=dhkl.GetPeakList().end();++pos)
   {
      epsilon +=fabs(pos->d2diff-zero);
   }
   if(nbSpurious>0)
   {// find worst fitting lines and remove them from epsilon calculation
      list<pair<float,unsigned int> > vdiff_idx;
      unsigned int i=0;
      for(pos=dhkl.GetPeakList().begin();pos!=dhkl.GetPeakList().end();++pos)
         vdiff_idx.push_back(make_pair(fabs(pos->d2diff),i++));
      vdiff_idx.sort(comparePairFirst<float,unsigned int>);
      i=0;
      for(list<pair<float,unsigned int> >::reverse_iterator rpos=vdiff_idx.rbegin();rpos!=vdiff_idx.rend();++rpos)
      {// :TODO: correct zero after removing spurious lines
         epsilon -= fabs(rpos->first-zero);
         if(storehkl) dhkl.GetPeakList()[rpos->second].isIndexed=false;
         dhkl.GetPeakList()[rpos->second].stats++;
         if(++i==nbSpurious) break;
      }
   }
   if(verbose)
   {
      float epstmp=0;
      unsigned long i=0;
      for(pos=dhkl.GetPeakList().begin();pos!=dhkl.GetPeakList().end();++pos)
      {
         epstmp+=fabs(pos->d2diff-zero);
         //cout<<"Line #"<<i++<<": obs="<<pos->d2obs<<", diff="<<pos->d2diff<<" -> epsilon="<<epstmp<<endl;
      }
      cout<<"epsilon="<<epstmp<<", dmax="<<dmax<<" ,nb="<<nb<<" ,nbcalc="<<nbCalc<<endl;
   }
   /*
   else
   {//Only stat+ the worst
      float max=-1;
      unsigned int worst=0;
      unsigned int i=0;
      for(pos=dhkl.GetPeakList().begin();pos!=dhkl.GetPeakList().end();++pos)
         if(abs(pos->d2diff)>max) {worst=i++;max=abs(pos->d2diff);}
         else i++;
      dhkl.GetPeakList()[worst].stats++;
   }
   */
   if(nbCalc==0) return 0;
   const float score=(dmax/predict_coeff)*nb/(2*epsilon*nbCalc);
   if(verbose)
   {
      dhkl.Print(cout);
      cout<<"Final score:"<<score<<", nbCalc="<<nbCalc<<" ,<epsilon>="<<epsilon<<" nb="<<nb<<" Qn="<<sqrt(dmax)<<endl;
   }
   return score;
}

/////////////////////////////////////////////////////// CellExplorer ///////////////////////////////////////

CellExplorer::CellExplorer(const PeakList &dhkl, const CrystalSystem lattice, const unsigned int nbSpurious):
mnpar(3),mpPeakList(&dhkl),
mLengthMin(4),mLengthMax(25),
mAngleMin(M_PI),mAngleMax(2*M_PI/3),
mVolumeMin(0),mVolumeMax(1600),
mZeroShiftMin(0),mZeroShiftMax(0),
mlattice(lattice),mNbSpurious(nbSpurious),
mObs(0),mCalc(0),mWeight(0),mDeriv(0),mBestScore(0.0),
mMinScoreReport(10),mMaxDicVolDepth(7),mDicVolDepthReport(6)
{
   this->Init();
}

void CellExplorer::Evolution(unsigned int ng,const bool randomize,const float f,const float cr,unsigned int np)
{
   this->Init();
   const bool autozero=true;
   //cout<<__FILE__<<":"<<__LINE__<<"<CellExplorer::Evolution(...): randomizing,ng="<<ng
   //    <<"random="<<randomize<<"f="<<f<<"cr="<<cr<<"np="<<np<<endl;
   vector<pair<RecUnitCell,float> > vRUC(np);
   vector<pair<RecUnitCell,float> > vTrial(np);
   float bestScore=-1e20;
   vector<pair<RecUnitCell,float> >::iterator bestpos=vRUC.begin();

   const clock_t mTime0=clock();
   
   if(randomize)
   {
      for(unsigned int i=0;i<vRUC.size();++i)
      {
         vRUC[i].first.mlattice=mlattice;
         vTrial[i].first.mlattice=mlattice;
         for(unsigned int k=0;k<mnpar;++k) vRUC[i].first.par[k]=mMin[k]+mAmp[k]*rand()/(float)RAND_MAX;
         vRUC[i].second=Score(*mpPeakList,vRUC[i].first,mNbSpurious);
      }
   }

   for(unsigned long i=ng;i>0;--i)
   {
      for(unsigned j=0;j<np;j++)
      {
         if(true)
         {// DE/rand/1/exp
            unsigned int r1=j,r2=j,r3=j;
            while(r1==j)r1=rand()%np;
            while((r2==j)||(r1==r2))r2=rand()%np;
            while((r3==j)||(r3==r1)||(r3==r2))r3=rand()%np;
            unsigned int ncr=1;//+(int)(cr*mnpar*rand()/(float)RAND_MAX);
            unsigned int ncr0=rand()%mnpar;
            RecUnitCell *t0=&(vTrial[j].first);
            RecUnitCell *c0=&(vRUC[j].first);
            RecUnitCell *c1=&(vRUC[r1].first);
            RecUnitCell *c2=&(vRUC[r2].first);
            RecUnitCell *c3=&(vRUC[r3].first);
            for(unsigned int k=0;k<mnpar;++k)t0->par[k] = c0->par[k];
            for(unsigned int k=0;k<ncr;++k)
            {
               const unsigned l=(ncr0+k)%mnpar;
               const float v1=c1->par[l]-mMin[l];
               const float v2=c2->par[l]-mMin[l];
               const float v3=c3->par[l]-mMin[l];
               t0->par[l]=mMin[l]+fmod(v1+f*(v2-v3)+3*mAmp[l],mAmp[l]);
            }
         }
         if(false)
         {// DE/rand-to-best/1/exp
            unsigned int r1=j,r2=j,r3=j;
            while(r1==j)r1=rand()%np;
            while((r2==j)||(r1==r2))r2=rand()%np;
            while((r3==j)||(r3==r1)||(r3==r2))r3=rand()%np;
            unsigned int ncr=1+(int)(cr*(mnpar-1)*rand()/(float)RAND_MAX);
            unsigned int ncr0=rand()%mnpar;
            RecUnitCell *t0=&(vTrial[j].first);
            RecUnitCell *c0=&(vRUC[j].first);
            //RecUnitCell *c1=&(vRUC[r1].first);
            RecUnitCell *c2=&(vRUC[r2].first);
            RecUnitCell *c3=&(vRUC[r3].first);
            RecUnitCell *best=&(bestpos->first);
            for(unsigned int k=0;k<6;++k)t0->par[k] = c0->par[k];//mMin[k]+mAmp[k]*rand()/(float)RAND_MAX;
            for(unsigned int k=0;k<ncr;++k)
            {
               const unsigned l=(ncr0+k)%mnpar;
               const float v0=c0->par[l]-mMin[l];
               //const float v1=c1->par[l]-mMin[l];
               const float v2=c2->par[l]-mMin[l];
               const float v3=c3->par[l]-mMin[l];
               const float vb=best->par[l]-mMin[l];
               t0->par[l]=mMin[l]+fmod(vb+f*(vb-v0)+f*(v2-v3)+5*mAmp[l],mAmp[l]);
            }
         }
         if(false)
         {// MC
            const float amp=.05/(1+i*.01);
            RecUnitCell *t0=&(vTrial[j].first);
            for(unsigned int k=0;k<6;++k)
            {
               
               t0->par[k] = mMin[k]+ fmod((float)(amp*mAmp[k]*(rand()/(float)RAND_MAX-0.5)+5*mAmp[k]),(float)mAmp[k]);
            }
         }
      }
      // Compute cost for all trials and select best
      vector<pair<RecUnitCell,float> >::iterator posTrial,pos;
      posTrial=vTrial.begin();
      pos=vRUC.begin();
      for(;posTrial!=vTrial.end();)
      {
         // If using auto-zero, fix zero parameter
         if(autozero) posTrial->first.par[0]=0;
         // Did we go beyond allowed volume ?
         switch(mlattice)
         {
            case TRICLINIC:
               break;
            case MONOCLINIC:
            {
               float v0=posTrial->first.par[1]*posTrial->first.par[2]*posTrial->first.par[3];
               while(v0<1/mVolumeMax)
               {
                  const unsigned int i=rand()%3+1;
                  posTrial->first.par[i]*=1/(mVolumeMax*v0)+1e-4;
                  if(posTrial->first.par[i]>(mMin[i]+mAmp[i])) posTrial->first.par[i]=mMin[i]+mAmp[i];
                  v0=posTrial->first.par[1]*posTrial->first.par[2]*posTrial->first.par[3];
               }
               break;
            }
            case ORTHOROMBIC:
               break;
            case HEXAGONAL:
               break;
            case RHOMBOEDRAL:
               break;
            case TETRAGONAL:
               break;
            case CUBIC:
               break;
         }

         const float score=Score(*mpPeakList,posTrial->first,mNbSpurious);
         if(score > pos->second)
         {
            pos->second=score;
            const REAL *p0=posTrial->first.par;
            REAL *p1=pos->first.par;
            for(unsigned int k=0;k<mnpar;++k) *p1++ = *p0++;
            if(score>bestScore)
            {
               bestScore=score;
               bestpos=pos;
            }
         }
         /*
         else
         {
            if(log(rand()/(float)RAND_MAX)>(-(score-pos->second)))
            {
               pos->second=score;
               const float *p0=posTrial->first.par;
               float *p1=pos->first.par;
               for(unsigned int k=0;k<mnpar;++k) *p1++ = *p0++;
            }
         }
         */
         ++pos;++posTrial;
      }
      if((i%100000)==0)
      {
         vector<float> uc=bestpos->first.DirectUnitCell();
         cout<<"Generation #"<<ng-i<<", Best score="<<bestScore
             <<" Trial: a="<<uc[0]<<", b="<<uc[1]<<", c="<<uc[2]<<", alpha="
             <<uc[3]*RAD2DEG<<", beta="<<uc[4]*RAD2DEG<<", gamma="<<uc[5]*RAD2DEG<<", V="<<uc[6]
             <<"   "<<(ng-i)*np/((clock()-mTime0)/(float)CLOCKS_PER_SEC)<<" trials/s"
             <<endl;
      }
      if(false)//((i%10000)==0)
      {// Randomize periodically
         for(vector<pair<RecUnitCell,float> >::iterator pos=vRUC.begin();pos!=vRUC.end();++pos)
         {
            if(pos==bestpos) continue;
            for(unsigned int k=0;k<mnpar;++k) pos->first.par[k]=mMin[k]+mAmp[k]*rand()/(float)RAND_MAX;
         }
      }
   }
   /*
   for(vector<pair<RecUnitCell,float> >::iterator pos=vRUC.begin();pos!=vRUC.end();++pos)
   {
      // Final cost
      vector<float> uc=pos->first.DirectUnitCell();
      cout<<__FILE__<<":"<<__LINE__<<" Trial: a="<<uc[0]<<", b="<<uc[1]<<", c="<<uc[2]<<", alpha="
       <<uc[3]*RAD2DEG<<", beta="<<uc[4]*RAD2DEG<<", gamma="<<uc[5]*RAD2DEG<<", V="<<uc[6]
          <<", score="<<pos->second<<endl;
   }
   Score(*mpPeakList,bestpos->first,mNbSpurious,true);
   */
   
   //this->ReduceSolutions(true);

   mRecUnitCell=bestpos->first;
   float score=Score(*mpPeakList,mRecUnitCell,mNbSpurious,false,true);
   vector<float> uc=mRecUnitCell.DirectUnitCell();
   cout<<__FILE__<<":"<<__LINE__<<" Best-DE : a="<<uc[0]<<", b="<<uc[1]<<", c="<<uc[2]<<", alpha="
       <<uc[3]*RAD2DEG<<", beta="<<uc[4]*RAD2DEG<<", gamma="<<uc[5]*RAD2DEG<<", V="<<uc[6]
       <<", score="<<bestpos->second
       <<"     ("<<ng*np/((clock()-mTime0)/(float)CLOCKS_PER_SEC)<<" trials/s)"<<endl;
   if(score>mMinScoreReport*.5)
   {
      // Now, do a least-squares refinement on best
      mRecUnitCell=bestpos->first;
      this->LSQRefine(10,true,true);
      uc=mRecUnitCell.DirectUnitCell();
      score=Score(*mpPeakList,mRecUnitCell,mNbSpurious,false,true);
      cout<<__FILE__<<":"<<__LINE__<<" Best-LSQ: a="<<uc[0]<<", b="<<uc[1]<<", c="<<uc[2]<<", alpha="
         <<uc[3]*RAD2DEG<<", beta="<<uc[4]*RAD2DEG<<", gamma="<<uc[5]*RAD2DEG<<", V="<<uc[6]
         <<", score="<<score<<endl;
      if((score>mMinScoreReport)&&(score>(mBestScore/3)))
      {
         if(score>mBestScore) mBestScore=score;
         mvSolution.push_back(make_pair(mRecUnitCell,score));
         this->ReduceSolutions(true);// We may have solutions from previous runs
      }
   }
}

void CellExplorer::SetLengthMinMax(const float min,const float max)
{
   mLengthMin=min;
   mLengthMax=max;
}
void CellExplorer::SetAngleMinMax(const float min,const float max)
{
   mAngleMin=min;
   mAngleMax=max;
}
void CellExplorer::SetVolumeMinMax(const float min,const float max)
{
   mVolumeMin=min;
   mVolumeMax=max;
}
void CellExplorer::SetNbSpurious(const unsigned int nb)
{
   mNbSpurious=nb;
}
void CellExplorer::SetMinMaxZeroShift(const float min,const float max)
{
   mZeroShiftMin=min;
   mZeroShiftMax=max;
}

void CellExplorer::SetCrystalSystem(const CrystalSystem system)
{
   mlattice=system;
}

void CellExplorer::SetCrystalCentering(const CrystalCentering cent)
{
   mCentering=cent;
}

void CellExplorer::SetD2Error(const float err){mD2Error=err;}

const string& CellExplorer::GetClassName() const
{
   const static string className="CellExplorer";
   return className;
}
const string& CellExplorer::GetName() const
{
   const static string name="Some CellExplorer Object";
   return name;
}
void CellExplorer::Print() const
{
   this->RefinableObj::Print();
}
unsigned int CellExplorer::GetNbLSQFunction() const
{return 1;}

const CrystVector_REAL& CellExplorer::GetLSQCalc(const unsigned int) const
{
   VFN_DEBUG_ENTRY("CellExplorer::GetLSQCalc()",2)
   unsigned int j=0;
   for(vector<PeakList::hkl>::const_iterator pos=mpPeakList->GetPeakList().begin();pos!=mpPeakList->GetPeakList().end();++pos)
   {
      if(pos->isIndexed)
         mCalc(j++)=mRecUnitCell.hkl2d(pos->h,pos->k,pos->l);
   }
   //cout<<__FILE__<<":"<<__LINE__<<"LSQCalc : Score:"<<Score(*mpPeakList,mRecUnitCell,mNbSpurious,false,true,false)<<endl;
   VFN_DEBUG_EXIT("CellExplorer::GetLSQCalc()",2)
   return mCalc;
}
const CrystVector_REAL& CellExplorer::GetLSQObs(const unsigned int) const
{
   VFN_DEBUG_MESSAGE("CellExplorer::GetLSQObs()",2)
   return mObs;
}
const CrystVector_REAL& CellExplorer::GetLSQWeight(const unsigned int) const
{
   VFN_DEBUG_MESSAGE("CellExplorer::GetLSQWeight()",2)
   //:TODO: exclude the worst points (user-chosen number)
   return mWeight;
}
const CrystVector_REAL& CellExplorer::GetLSQDeriv(const unsigned int, RefinablePar &refpar)
{
   VFN_DEBUG_ENTRY("CellExplorer::GetLSQDeriv()",2)
   REAL *par=NULL;
   if(refpar.GetName()=="Reciprocal unit cell par #0") par=mRecUnitCell.par+1;
   else
      if(refpar.GetName()=="Reciprocal unit cell par #1") par=mRecUnitCell.par+2;
      else
         if(refpar.GetName()=="Reciprocal unit cell par #2") par=mRecUnitCell.par+3;
         else
            if(refpar.GetName()=="Reciprocal unit cell par #3") par=mRecUnitCell.par+4;
            else
               if(refpar.GetName()=="Reciprocal unit cell par #4") par=mRecUnitCell.par+5;
               else
                  if(refpar.GetName()=="Reciprocal unit cell par #5") par=mRecUnitCell.par+6;
                  else 
                     if(refpar.GetName()=="Zero") par=mRecUnitCell.par+0;
                     else cout<<__FILE__<<":"<<__LINE__<<":Parameter not found:"<<refpar.GetName()<<endl;
   unsigned int j=0;
   for(vector<PeakList::hkl>::const_iterator pos=mpPeakList->GetPeakList().begin();pos!=mpPeakList->GetPeakList().end();++pos)
   {
      VFN_DEBUG_MESSAGE("CellExplorer::GetLSQDeriv():"<<j<<"/"<<mpPeakList->GetPeakList().size(),2)
      VFN_DEBUG_MESSAGE("CellExplorer::GetLSQDeriv():"<<pos->h<<","<<pos->k<<","<<pos->l,2)
      if(pos->isIndexed) 
         mDeriv(j++)=mRecUnitCell.hkl2d(pos->h,pos->k,pos->l,par);
   }
   VFN_DEBUG_EXIT("CellExplorer::GetLSQDeriv()",2)
   return mDeriv;
}

void CellExplorer::BeginOptimization(const bool allowApproximations, const bool enableRestraints)
{
   VFN_DEBUG_ENTRY("CellExplorer::BeginOptimization()",10)
   Score(*mpPeakList,mRecUnitCell,mNbSpurious,false,true,false);
   const unsigned int nb=mpPeakList->GetPeakList().size();
   mCalc.resize(nb-mNbSpurious);
   mObs.resize(nb-mNbSpurious);
   mWeight.resize(nb-mNbSpurious);
   mDeriv.resize(nb-mNbSpurious);
   int j=0;
   float thres=0.0;
   for(vector<PeakList::hkl>::const_iterator pos=mpPeakList->GetPeakList().begin();pos!=mpPeakList->GetPeakList().end();++pos)
      if(thres<pos->iobs) thres=pos->iobs;
   thres/=10;// weight=1 for intensities up to Imax/10

   //cout <<"Beginning optimization with reflexions:"<<endl;
   //char buf[100];
   for(vector<PeakList::hkl>::const_iterator pos=mpPeakList->GetPeakList().begin();pos!=mpPeakList->GetPeakList().end();++pos)
   {
      if(pos->isIndexed)
      {
         mObs(j)=pos->d2obs;
         if(mObs(j)>thres) mWeight(j)=1;
         else mWeight(j)=mObs(j)/thres;
         /*
         sprintf(buf,"#%2d  (%3d %3d %3d) dobs=%6.3f dcalc=%6.3f iobs=%6.3f weight=%6.4f",
                 i,mpPeakList->mvHKL[i].h,mpPeakList->mvHKL[i].k,mpPeakList->mvHKL[i].l,
                 1/mpPeakList->mvdobs[i],1/sqrt(mRecUnitCell.hkl2d(mpPeakList->mvHKL[i].h,mpPeakList->mvHKL[i].k,mpPeakList->mvHKL[i].l)),
                 mpPeakList->mviobs[i],mWeight(j));
         */
         j++;
      }
      /*
      else
      {
         sprintf(buf,"#%2d  (%3d %3d %3d) dobs=%6.3f dcalc=%6.3f iobs=%6.3f               SPURIOUS",
                 i,mpPeakList->mvHKL[i].h,mpPeakList->mvHKL[i].k,mpPeakList->mvHKL[i].l,
                 1/mpPeakList->mvdobs[i],1/sqrt(mRecUnitCell.hkl2d(mpPeakList->mvHKL[i].h,mpPeakList->mvHKL[i].k,mpPeakList->mvHKL[i].l)),
                 mpPeakList->mviobs[i]);
      }
      cout<<buf<<endl;
      */
   }
   this->RefinableObj::BeginOptimization(allowApproximations,enableRestraints);
   VFN_DEBUG_EXIT("CellExplorer::BeginOptimization()",10)
}

void CellExplorer::LSQRefine(int nbCycle, bool useLevenbergMarquardt, const bool silent)
{
   VFN_DEBUG_ENTRY("CellExplorer::LSQRefine()",5)
   mLSQObj.SetRefinedObj(*this);
   //this->BeginOptimization();
   //cout<<FormatVertVector<REAL>(this->GetLSQObs(0),this->GetLSQCalc(0),this->GetLSQWeight(0),this->GetLSQDeriv(0,this->GetPar((long)0)))<<endl;
   mLSQObj.Refine(nbCycle,useLevenbergMarquardt,silent);
   if(!silent) mpPeakList->Print(cout);
   VFN_DEBUG_EXIT("CellExplorer::LSQRefine()",5)
}

/** Number of reflexions found in the intervals calculated between uc+duc and uc-duc
*
* \param useStoredHKL:
*    - if equal to 0, explore all possible hkl values to find possible Miller indices.
*    - if useStoredHKL=1, use the Miller indices already stored in hkl.vDicVolHKL
*    for each observed line as the only possible indices.
*    - if useStoredHKL=2, search all the possible Miller indices for all reflections
*    and store them in hkl.vDicVolHKL for each observed line.
* \param maxNbMissingBelow5: the maximum number of lines that have been calculated before 
* the d-value of the 5th observed line, but which have not been observed. 
* If nbMissingBelow5>0 and more than nbMissingBelow5 lines have not been observed,
* return false. Recommended to speed up triclinic, with nbMissingBelow5=5
*/
bool DichoIndexed(const PeakList &dhkl, const RecUnitCell &uc,const RecUnitCell &duc,
                  const unsigned int nbUnindexed=0,const bool verbose=false,unsigned int useStoredHKL=0,
                  const unsigned int maxNbMissingBelow5=0)
{
   const unsigned int nb=dhkl.GetPeakList().size();
   int nbIndexed=nb-nbUnindexed;// Number of reflections we require to be indexed
   float d5=0;
   if(maxNbMissingBelow5>0) d5=dhkl.GetPeakList()[4].d2obs;
   // number of missing reflections calculated below 5th observed line
   unsigned int nbMissingBelow5=0;
   // List of indexed reflections
   vector<PeakList::hkl>::const_iterator pos,first,last,end;
   if(useStoredHKL==1)
   {// We already now possible Miller indices for all reflections
      int nbUnIx=0;
      for(pos=dhkl.GetPeakList().begin();pos!=dhkl.GetPeakList().end();++pos)
      {
         pos->isIndexed=false;
         for(list<PeakList::hkl0>::const_iterator phkl0=pos->vDicVolHKL.begin();phkl0!=pos->vDicVolHKL.end();++phkl0)
         {
            float d0,d1;
            uc.hkl2d_delta(phkl0->h,phkl0->k,phkl0->l,duc,d0,d1);
            if((pos->d2obsmax>=d0) && (d1>=pos->d2obsmin))
            {
               pos->d2calc=(d0+d1)/2;
               pos->isIndexed=true;
               if(--nbIndexed==0) return true;
               break;
            }
         }
         if(!(pos->isIndexed)) if(++nbUnIx>nbUnindexed) return false;
      }
      return false;
   }
   const bool storePossibleHKL=(useStoredHKL==2);
   
   if(storePossibleHKL)
      for(pos=dhkl.GetPeakList().begin();pos!=dhkl.GetPeakList().end();++pos)
      {
         pos->isIndexed=false;
         pos->vDicVolHKL.clear();
      }
   else
      for(pos=dhkl.GetPeakList().begin();pos!=dhkl.GetPeakList().end();++pos) pos->isIndexed=false;
   
   int h,k,l;
   float dmax=dhkl.GetPeakList()[nb-1].d2obs;
   float dmin=dhkl.GetPeakList()[0   ].d2obs;
   
   
   int sk0,sl0;// do we need >0 *and* <0 indices for k,l ?
   switch(uc.mlattice)
   {
      case TRICLINIC:
         sk0=-1;sl0=-1;
	 break;
      case MONOCLINIC:
      {
         sk0=1;sl0=-1;
         break;
      }
      case ORTHOROMBIC:
         sk0=1;sl0=1;
         break;
      case HEXAGONAL:
         sk0=-1;sl0=1;
         break;
      case RHOMBOEDRAL:
         sk0=-1;sl0=-1;
         break;
      case TETRAGONAL:
         sk0=1;sl0=1;
         break;
      case CUBIC:
         sk0=1;sl0=1;
         break;
   }
   int stepk,stepl;// steps in k,l to use for centered lattices
   switch(uc.mCentering)
   {
      case LATTICE_P:stepk=1;stepl=1;break;
      case LATTICE_I:stepk=1;stepl=2;break;
      case LATTICE_A:stepk=1;stepl=2;break;
      case LATTICE_B:stepk=1;stepl=2;break;
      case LATTICE_C:stepk=2;stepl=1;break;
      case LATTICE_F:stepk=2;stepl=2;break;
   }
   RecUnitCell uc0(uc),uc1(uc);
   for(unsigned int i=0;i<7;++i) {uc0.par[i]-=duc.par[i];uc1.par[i]+=duc.par[i];}
   
   //currently first & last unindexed dhkl
   first=dhkl.GetPeakList().begin(),last=dhkl.GetPeakList().end(),end=dhkl.GetPeakList().end();
   
   bool break_h,break_k;//0: continue, 1:break if only >0 to be explored, 2: break
   for(h=0;;++h)
   {
      break_h=false;
      for(int sk=sk0;sk<=1;sk+=2)
      {
         if(h==0) sk=1;
         if(stepk==2) k=(h%2);// For LATTICE_C,LATTICE_F: h odd => k odd
         else k=0;
         for(;;k+=stepk)
         {
            break_k=false;
            for(int sl=sl0;sl<=1;sl+=2)
            {
               int l0;
               if((h+k)==0)
               {
                  sl=1;// No need to list 0 0 l with l<0
                  l0=1;
               }
               else
               {
                  if(h==0)
                  {
                     if(uc.mlattice==MONOCLINIC) sl=1;// 0 k l and 0 k -l are equivalent
                     if((sk<0)||(sl<0)) l0=1;// Do not list 0 k 0 with k<0
                     else l0=0;// h==k==0 already covered
                  }
                  else
                  {
                     if(sl<0) l0=1;// Do not list h k 0 twice
                     else l0=0;
                  }
               }
               if(stepl==2)
               {
                  if(uc.mCentering==LATTICE_I) l0+=(h+k+l0)%2;
                  if(uc.mCentering==LATTICE_A) l0+=(k+l0)%2;// Start at k+l even
                  if(  (uc.mCentering==LATTICE_B)
                     ||(uc.mCentering==LATTICE_F)) l0+=(h+l0)%2;// Start at h+l even
               }
               for(l=l0;;l+=stepl)
               {
                  float d0,d1;
                  uc.hkl2d_delta(h,sk*k,sl*l,duc,d0,d1);
                  if((d1<dmin)&&(maxNbMissingBelow5==0)) continue;
                  if(d0>dmax)
                  {
                     if(uc.mlattice==TRICLINIC)
                     {
                        // Must check that d is increasing with l, otherwise we still need to increase it
                        if((uc.hkl2d(h,sk*k,sl*l,NULL,3)*sl)>0)
                        {
                           if(l==l0) break_k=true; // d(hk0) is always increasing with k
                           break;
                        }
                     }
                     else
                     {
                        if(l==l0) break_k=true;
                        break;
                     }
                  }
                  bool missing=(d0<d5)&&(maxNbMissingBelow5>0);
                  for(pos=first;pos!=end;++pos)
                  {
                     if(pos==last) break;
                     if((!storePossibleHKL)&&(pos->isIndexed)&&missing) continue;
                     const float d2obs=pos->d2obs,d2obsmin=pos->d2obsmin, d2obsmax=pos->d2obsmax;
                     if((d2obsmax>=d0) && (d1>=d2obsmin))
                     {
                        missing=false;
                        if(!(pos->isIndexed))
                        {
                           pos->d2calc=(d0+d1)/2;
                           --nbIndexed;
                           pos->isIndexed=true;
                        }
                        if(verbose) cout<<d1<<" < ? <"<<d0<<"("<<h<<","<<sk*k<<","<<sl*l<<"): "<<d2obs<<" (remaining to index:"<<nbIndexed<<")"<<endl;
                        if(storePossibleHKL)
                           pos->vDicVolHKL.push_back(PeakList::hkl0(h,sk*k,sl*l));
                        else
                        {
                           if((!storePossibleHKL)&&(nbIndexed==0)) return true;
                           if(pos==first){first++;dmin=first->d2obsmin;}
                           if(pos==last){last--;dmax=last->d2obsmax;}
                        }
                     }
                  }
                  if(missing) if(++nbMissingBelow5>=maxNbMissingBelow5)return false;
               }
            }
            if(break_k)
            {
               if(k==0) break_h=true; // We hit the limit for k==0, h too large
               break; // hk0 beyond limit
            }
         }
      }
      if(break_h) break;//h00 beyond limit
   }
   if(verbose)
   {
      dhkl.Print(cout);
   }
   return nbIndexed<=0;
}

float CellExplorer::GetBestScore()const{return mBestScore;}
const std::list<std::pair<RecUnitCell,float> >& CellExplorer::GetSolutions()const {return mvSolution;}
std::list<std::pair<RecUnitCell,float> >& CellExplorer::GetSolutions() {return mvSolution;}

unsigned int CellExplorer::RDicVol(RecUnitCell uc0,RecUnitCell duc, unsigned int depth,unsigned long &nbCalc,const float minV,const float maxV,vector<unsigned int> vdepth)
{
   static bool localverbose=false;
   if(mlattice==TRICLINIC)
   {
      const float p1m=uc0.par[1]-duc.par[1], p2m=uc0.par[2]-duc.par[2], p3m=uc0.par[3]-duc.par[3], p4m=uc0.par[4]-duc.par[4], p5m=uc0.par[5]-duc.par[5], p6m=uc0.par[6]-duc.par[6];
      const float p1p=uc0.par[1]+duc.par[1], p2p=uc0.par[2]+duc.par[2], p3p=uc0.par[3]+duc.par[3], p4p=uc0.par[4]+duc.par[4], p5p=uc0.par[5]+duc.par[5], p6p=uc0.par[6]+duc.par[6];
      
      #if 0
      // XCell-07-triclinic-00375-07.98-09.30-05.85-105.14-099.43-063.67.dat
      if((  ((p1m-.1229)*(p1p-.1229)<0) && ((p2m-.14024)*(p2p-.14024)<0) && ((p3m-.1776)*(p3p-.1776)<0)
          &&((p4m-.4223)*(p4p-.4223)<0) && ((p5m+.05395)*(p5p+.05395)<0) && ((p6m-.2166)*(p6p-.2166)<0))
          ||(localverbose && (depth==3)))
      {
         bool indexed=DichoIndexed(*mpPeakList,uc0,duc,mNbSpurious,false,2,5);
         
         RecUnitCell ucm=uc0,ucp=uc0;
         for(unsigned int i=0;i<4;++i) {ucm.par[i]-=duc.par[i];ucp.par[i]+=duc.par[i];}
         for(unsigned int i=4;i<7;++i) {ucm.par[i]+=duc.par[i];ucp.par[i]-=duc.par[i];}
         char buf[200];
         sprintf(buf,"a=%7.5f-%7.5f b=%7.5f-%7.5f c=%7.5f-%7.5f alpha=%7.5f-%7.5f beta=%7.5f-%7.5f gamma=%7.5f-%7.5f",
                  ucp.par[1],ucm.par[1],ucp.par[2],ucm.par[2],ucp.par[3],ucm.par[3],ucp.par[4],ucm.par[4],
                  ucp.par[5],ucm.par[5],ucp.par[6],ucm.par[6]);
         for(int i=0;i<depth;++i) cout<<" ";
         cout<<"AA"<<buf<<"level="<<depth<<",indexed="<<indexed<<"("<<mvSolution.size()<<" sol.)"<<endl;
         //if(indexed && depth==3) mpPeakList->Print(cout);
         //if(!indexed) mpPeakList->Print(cout);
         if(depth>=2)localverbose=true;
      }
      else if(depth<2)localverbose=false;
      #endif
      // a*<b*<c*
      if((p1m>p2p)||(p2m>p3p)) return 0;
      
      // max/min absolute values for p4,p5,p6
      float ap5m,ap5p;
      if(abs(p5m)<abs(p5p)) {ap5p=abs(p5p);ap5m=abs(p5m);} else {ap5p=abs(p5m);ap5m=abs(p5p);}
      float ap6m,ap6p;
      if(abs(p6m)<abs(p6p)) {ap6p=abs(p6p);ap6m=abs(p6m);} else {ap6p=abs(p6m);ap6m=abs(p6p);}
      float ap4m,ap4p;
      if(abs(p4m)<abs(p4p)) {ap4p=abs(p4p);ap4m=abs(p4m);} else {ap4p=abs(p4m);ap4m=abs(p4p);}
      
      #ifdef TRI1
      if(ap4m>(p1p*p1p)) return 0;// |b| <= |b +/- a|
      if(ap6m>(p1p*p1p)) return 0;// |c| <= |c +/- a|
      if(ap5m>(p2p*p2p)) return 0;// |c| <= |c +/- b|
      
      if((-p6p-p5p-p4p)>(p1p*p1p+p2p*p2p)) return 0;// |c| <= |c + a + b|
      if((+p6p-p5p+p4p)>(p1p*p1p+p2p*p2p)) return 0;// |c| <= |c - a + b|
      if((-p6p+p5p+p4p)>(p1p*p1p+p2p*p2p)) return 0;// |c| <= |c + a - b|
      if((+p6p+p5p-p4p)>(p1p*p1p+p2p*p2p)) return 0;// |c| <= |c - a - b|
      
      const float cam=p5m/(2*p2p*p3p),cbm=p6m/(2*p1p*p3p),ccm=p4m/(2*p1p*p2p),
                  cap=p5p/(2*p2m*p3m),cbp=p6p/(2*p1m*p3m),ccp=p4m/(2*p1m*p2m);
      // Min square volume of reciprocal UC
      const float rvminsq=p1m*p1m*p2m*p2m*p3m*p3m*(1-cap*cap-cbp*cbp-ccp*ccp+2*cap*cbp*ccp);
      // Max square volume of reciprocal UC
      const float rvmaxsq=p1p*p1p*p2p*p2p*p3p*p3p*(1-cam*cam-cbm*cbm-ccm*ccm+2*cam*cbm*ccm);
      if(( (1/rvmaxsq)>(maxV*maxV) ) || ( (1/rvminsq)<(minV*minV) ))
      {
         //cout<<"Volume out of range:"<<rvminsq<<", "<<rvmaxsq<<endl;
         return 0;
      }
      #if 0
      // Max/min value for direct cell a=sin(alpha*) * b* * c* / v*
      const float aminsq=(p2m*p2m*p3m*p3m-0.25*p5p)/abs(rvmaxsq);
      if(aminsq>mLengthMax*mLengthMax)
      {
         //cout<<"Maxed a:"<<sqrt(aminsq)<<":"<<ucmd[0]<<"-"<<ucpd[0]<< " :"<<p4p<<","<<p5p<<","<<p6p<<","<<rvmaxsq<<","<<rvminsq<<endl;
         return 0;
      }
      const float amaxsq=(p2p*p2p*p3p*p3p-0.25*p5m)/abs(rvminsq);
      if(amaxsq<mLengthMin*mLengthMin)
      {
         //cout<<"Maxed a:"<<sqrt(aminsq)<<":"<<ucmd[0]<<"-"<<ucpd[0]<< " :"<<p4p<<","<<p5p<<","<<p6p<<","<<rvmaxsq<<","<<rvminsq<<endl;
         return 0;
      }
      
      // Max/min value for direct cell b=sin(beta*) * a* * c* / v*
      const float bminsq=(p1m*p1m*p3m*p3m-0.25*p6p)/abs(rvmaxsq);
      if(bminsq>mLengthMax*mLengthMax)
      {
         return 0;
      }
      const float bmaxsq=(p1p*p1p*p3p*p3p-0.25*p6m)/abs(rvminsq);
      if(bmaxsq<mLengthMin*mLengthMin)
      {
         return 0;
      }
      
      // Max/min value for direct cell a=sin(alpha*) * b* * c* / v*
      const float cminsq=(p1m*p1m*p2m*p2m-0.25*p4p)/abs(rvmaxsq);
      if(cminsq>mLengthMax*mLengthMax)
      {
         return 0;
      }
      const float cmaxsq=(p1p*p1p*p2p*p2p-0.25*p4m)/abs(rvminsq);
      if(cmaxsq<mLengthMin*mLengthMin)
      {
         return 0;
      }
      #endif
      #else
      // Min square volume of reciprocal UC
      const float rvminsq=p1m*p1m*p2m*p2m*p3m*p3m*(1-p4m*p4m-p5m*p5m-p6m*p6m+2*p4m*p5m*p6m);
      // Max square volume of reciprocal UC
      const float rvmaxsq=p1p*p1p*p2p*p2p*p3p*p3p*(1-p4p*p4p-p5p*p5p-p6p*p6p+2*p4p*p5p*p6p);
      if(( (1/rvmaxsq)>(maxV*maxV) ) || ( (1/rvminsq)<(minV*minV) ))
      {
         //cout<<"Volume out of range:"<<rvminsq<<", "<<rvmaxsq<<endl;
         return 0;
      }
      
      //RecUnitCell ucm=uc0,ucp=uc0;
      //for(unsigned int i=0;i<6;++i) {ucm.par[i]-=duc.par[i];ucp.par[i]+=duc.par[i];}
      //vector<float> ucmd=ucm.DirectUnitCell();
      //vector<float> ucpd=ucp.DirectUnitCell();
      if(ap4m>(p1p/(2*p2m))) return 0;// |b| <= |b +/- a|
      if(ap6m>(p1p/(2*p3m))) return 0;// |c| <= |c +/- a|
      if(ap5m>(p2p/(2*p3m))) return 0;// |c| <= |c +/- b|
      
      // 0.52 factor to keep a margin of error ?
      float abm,abp,acm,acp,bcm,bcp;// min, max of scalar products
      if(uc0.par[4]>0){abm=p1m*p2m*p4m;abp=p1p*p2p*p4p;}
      else {abm=p1p*p2p*p4m;abp=p1m*p2m*p4p;}
      if(uc0.par[6]>0){acm=p1m*p3m*p6m;acp=p1p*p3p*p6p;}
      else {acm=p1p*p3p*p6m;acp=p1m*p3m*p6p;}
      if(uc0.par[5]>0){bcm=p2m*p3m*p5m;bcp=p2p*p3p*p5p;}
      else {bcm=p2p*p3p*p5m;bcp=p2m*p3m*p5p;}
      
      if((-abp-acp-bcp)>(0.52*(p1p*p1p+p2p*p2p))) return 0;// |c| <= |c + a + b|
      if((+abm+acm-bcp)>(0.52*(p1p*p1p+p2p*p2p))) return 0;// |c| <= |c - a + b|
      if((+abm-acp+bcm)>(0.52*(p1p*p1p+p2p*p2p))) return 0;// |c| <= |c + a - b|
      if((-abp+acm+bcm)>(0.52*(p1p*p1p+p2p*p2p))) return 0;// |c| <= |c - a - b|
      
      // Max/min value for direct cell a=sin(alpha*) * b* * c* / v*
      const float aminsq=(1-ap5p*ap5p)*p2m*p2m*p3m*p3m/abs(rvmaxsq);
      if(aminsq>mLengthMax*mLengthMax)
      {
         //cout<<"Maxed a:"<<sqrt(aminsq)<<":"<<ucmd[0]<<"-"<<ucpd[0]<< " :"<<p4p<<","<<p5p<<","<<p6p<<","<<rvmaxsq<<","<<rvminsq<<endl;
         return 0;
      }
      const float amaxsq=(1-ap5m*ap5m)*p2p*p2p*p3p*p3p/abs(rvminsq);
      if(amaxsq<mLengthMin*mLengthMin)
      {
         //cout<<"Maxed a:"<<sqrt(aminsq)<<":"<<ucmd[0]<<"-"<<ucpd[0]<< " :"<<p4p<<","<<p5p<<","<<p6p<<","<<rvmaxsq<<","<<rvminsq<<endl;
         return 0;
      }
      
      // Max/min value for direct cell b=sin(beta*)/(b* * v*)
      const float bminsq=(1-ap6p*ap6p)*p1m*p1m*p3m*p3m/abs(rvmaxsq);
      if(bminsq>mLengthMax*mLengthMax)
      {
         //cout<<"Maxed b:"<<sqrt(bminsq)<<":"<<ucmd[1]<<"-"<<ucpd[1]<<endl;
         return 0;
      }
      const float bmaxsq=(1-ap6m*ap6m)*p1p*p1p*p3p*p3p/abs(rvminsq);
      if(bmaxsq<mLengthMin*mLengthMin)
      {
         //cout<<"Maxed b:"<<sqrt(bminsq)<<":"<<ucmd[1]<<"-"<<ucpd[1]<<endl;
         return 0;
      }
      
      // Max/min value for direct cell c=sin(gamma*)/(c* * v*)
      const float cminsq=(1-ap4p*ap4p)*p1m*p1m*p2m*p2m/abs(rvmaxsq);
      if(cminsq>mLengthMax*mLengthMax)
      {
         //cout<<"Maxed c:"<<sqrt(cminsq)<<":"<<ucmd[2]<<"-"<<ucpd[2]<<endl;
         return 0;
      }
      const float cmaxsq=(1-ap4m*ap4m)*p1p*p1p*p2p*p2p/abs(rvminsq);
      if(cmaxsq<mLengthMin*mLengthMin)
      {
         //cout<<"Maxed c:"<<sqrt(cminsq)<<":"<<ucmd[2]<<"-"<<ucpd[2]<<endl;
         return 0;
      }
      #endif
   }
   else
      if((depth>0)&&(depth<=2))// test if volume is within range
      {
         RecUnitCell ucm=uc0,ucp=uc0;
         for(unsigned int i=0;i<6;++i) {ucm.par[i]-=duc.par[i];ucp.par[i]+=duc.par[i];}
         vector<float> ucmd=ucm.DirectUnitCell();
         vector<float> ucpd=ucp.DirectUnitCell();
         if((ucpd[6]>maxV)||(ucmd[6]<minV))return 0;
      }
   unsigned int useStoredHKL=1;//Use already stored hkl
   if(depth==0) useStoredHKL=2; //Store possible hkl for all observed lines
   
   unsigned int maxMissingBelow5=0;
   // In the triclinic case, accept a maximum of 5 missing reflections below the 5th observed line
   if(mlattice==TRICLINIC) maxMissingBelow5=5;
   
   bool indexed=DichoIndexed(*mpPeakList,uc0,duc,mNbSpurious,false,useStoredHKL,maxMissingBelow5);
   
   #if 0
   // If indexation failed but depth>=4, try adding a zero ?
   if( (!indexed) && (depth>=4))
   {//:TODO: Check if this is OK ! Vary value with depth 
      duc.par[0]=.0001;
      indexed=DichoIndexed(*mpPeakList,uc0,duc,mNbSpurious,false,useStoredHKL,maxMissingBelow5);
      //if(indexed) cout<<"Added zero - SUCCESS !"<<endl;
   }
   #endif
   
   if((indexed)&&(useStoredHKL==2))
   {
      // Test if two successive lines have been indexed exclusively with the same hkl
      unsigned int nbident=0;
      for(vector<PeakList::hkl>::const_iterator pos=mpPeakList->GetPeakList().begin();pos!=mpPeakList->GetPeakList().end();)
      {
         if(pos->vDicVolHKL.size()==1)
         {
            const PeakList::hkl0 h0=pos->vDicVolHKL.front();
            if(++pos==mpPeakList->GetPeakList().end()) break;
            if(pos->vDicVolHKL.size()==1)
            {
               const PeakList::hkl0 h1=pos->vDicVolHKL.front();
               if((h0.h==h1.h)&&(h0.k==h1.k)&&(h0.l==h1.l)) nbident++;
               if(nbident>mNbSpurious) {indexed=false;break;}
            }
         }
         else ++pos;
      }
   }
   
   // if we can zoom in for one parameter directly, we need per-parameter depth
   if(vdepth.size()==0)
   {
      vdepth.resize(mnpar-1);
      for(vector<unsigned int>::iterator pos=vdepth.begin();pos!=vdepth.end();) *pos++=depth;
   }
   else
      for(vector<unsigned int>::iterator pos=vdepth.begin();pos!=vdepth.end();++pos) if(*pos<depth)*pos=depth;
   #if 1
   if((useStoredHKL==2)&&(mNbSpurious==0)&&indexed)
   {  // If high-d lines have been associated to a single reflection which is either h00, 0k0 or 00l,
      // jump the corresponding parameter to higher depth (mDicVolDepthReport, lowest depth report) immediately
      vector<pair<unsigned int,float> > vq0(3);
      for(unsigned int i=0;i<3;++i) {vq0[i].first=0;vq0[i].second=0.0;}
      RecUnitCell uc0orig=uc0,ducorig=duc;
      for(vector<PeakList::hkl>::const_iterator pos=mpPeakList->GetPeakList().begin();pos!=mpPeakList->GetPeakList().end();++pos)
      {
         if(pos->vDicVolHKL.size()==1)
         {
            const PeakList::hkl0 h0=pos->vDicVolHKL.front();
            if((h0.k==0)&&(h0.l==0)) {vq0[0].first+=1 ; vq0[0].second+=pos->dobs/h0.h;}
            else
               if((h0.h==0)&&(h0.l==0)) {vq0[1].first+=1 ; vq0[1].second+=pos->dobs/h0.k;}
               else
                  if((h0.h==0)&&(h0.k==0)) {vq0[2].first+=1 ; vq0[2].second+=pos->dobs/h0.l;if(localverbose) cout<<h0.h<<" "<<h0.k<<" "<<h0.l<<": d="<<pos->dobs<<endl;}
         }
      }
      switch(mlattice)
      {
         case TRICLINIC:
         {// In the triclinic case we may start with p1 and p2 already at depth>0 (triclinic quick tests)
            if(vq0[0].first>0) {uc0.par[1]=vq0[0].second/vq0[0].first ; duc.par[1]*=pow((float)0.5,float(mDicVolDepthReport-vdepth[0]));vdepth[0]=mDicVolDepthReport;}
            if(vq0[1].first>0) {uc0.par[2]=vq0[1].second/vq0[1].first ; duc.par[2]*=pow((float)0.5,float(mDicVolDepthReport-vdepth[1]));vdepth[1]=mDicVolDepthReport;}
            if(vq0[2].first>0) {uc0.par[3]=vq0[2].second/vq0[2].first ; duc.par[3]*=pow((float)0.5,float(mDicVolDepthReport-vdepth[2]));vdepth[2]=mDicVolDepthReport;}
            break;
         }
         case MONOCLINIC:
         {
            if(vq0[0].first>0) {uc0.par[1]=vq0[0].second/vq0[0].first ; vdepth[0]=mDicVolDepthReport;duc.par[1]*=.0625;}
            if(vq0[1].first>0) {uc0.par[2]=vq0[1].second/vq0[1].first ; vdepth[1]=mDicVolDepthReport;duc.par[2]*=.0625;}
            if(vq0[2].first>0) {uc0.par[3]=vq0[2].second/vq0[2].first ; vdepth[2]=mDicVolDepthReport;duc.par[3]*=.0625;}
            break;                                           
         }                                                   
         case ORTHOROMBIC:
         {
            if(vq0[0].first>0) {uc0.par[1]=vq0[0].second/vq0[0].first ; vdepth[0]=mDicVolDepthReport;duc.par[1]*=.0625;}//pow((float)0.5,(int)(mDicVolDepthReport-depth));}
            if(vq0[1].first>0) {uc0.par[2]=vq0[1].second/vq0[1].first ; vdepth[1]=mDicVolDepthReport;duc.par[2]*=.0625;}//pow((float)0.5,(int)(mDicVolDepthReport-depth));}
            if(vq0[2].first>0) {uc0.par[3]=vq0[2].second/vq0[2].first ; vdepth[2]=mDicVolDepthReport;duc.par[3]*=.0625;}//pow((float)0.5,(int)(mDicVolDepthReport-depth));}
            break;
         }
         case HEXAGONAL:break;
         case RHOMBOEDRAL:break;
         case TETRAGONAL:break;
         case CUBIC:break;
      }
      // If all parameters are at a higher depth, jump the global depth immediately
      unsigned int newdepth=40;
      for(vector<unsigned int>::iterator pos=vdepth.begin();pos!=vdepth.end();++pos) if(*pos<newdepth) newdepth=*pos;
      if(newdepth>depth) depth=newdepth;
      if((vq0[0].first>0)||(vq0[1].first>0)||(vq0[2].first>0))
      {
         indexed=DichoIndexed(*mpPeakList,uc0,duc,mNbSpurious,false,1,maxMissingBelow5);
         if(false)
         {
            {
               RecUnitCell ucm=uc0orig,ucp=uc0;
               for(unsigned int i=0;i<6;++i) {ucm.par[i]-=ducorig.par[i];ucp.par[i]+=ducorig.par[i];}
               vector<float> ucmd=ucm.DirectUnitCell();
               vector<float> ucpd=ucp.DirectUnitCell();
               char buf[200];
               sprintf(buf,"orig:   a=%5.2f-%5.2f b=%5.2f-%5.2f c=%5.2f-%5.2f alpha=%5.2f-%5.2f beta=%5.2f-%5.2f gamma=%5.2f-%5.2f V=%5.2f-%5.2f",
                        ucpd[0],ucmd[0],ucpd[1],ucmd[1],ucpd[2],ucmd[2],ucpd[3]*RAD2DEG,ucmd[3]*RAD2DEG,
                        ucpd[4]*RAD2DEG,ucmd[4]*RAD2DEG,ucpd[5]*RAD2DEG,ucmd[5]*RAD2DEG,ucpd[6],ucmd[6]);
               for(int i=0;i<depth;++i) cout<<" ";
               cout<<buf<<"level="<<depth<<", indexed="<<indexed<<endl;
            }
            {
               RecUnitCell ucm=uc0,ucp=uc0;
               for(unsigned int i=0;i<6;++i) {ucm.par[i]-=duc.par[i];ucp.par[i]+=duc.par[i];}
               vector<float> ucmd=ucm.DirectUnitCell();
               vector<float> ucpd=ucp.DirectUnitCell();
               char buf[200];
               sprintf(buf,"bypass: a=%5.2f-%5.2f b=%5.2f-%5.2f c=%5.2f-%5.2f alpha=%5.2f-%5.2f beta=%5.2f-%5.2f gamma=%5.2f-%5.2f V=%5.2f-%5.2f",
                        ucpd[0],ucmd[0],ucpd[1],ucmd[1],ucpd[2],ucmd[2],ucpd[3]*RAD2DEG,ucmd[3]*RAD2DEG,
                        ucpd[4]*RAD2DEG,ucmd[4]*RAD2DEG,ucpd[5]*RAD2DEG,ucmd[5]*RAD2DEG,ucpd[6],ucmd[6]);
               for(int i=0;i<depth;++i) cout<<" ";
               cout<<buf<<"level="<<depth<<", indexed="<<indexed<<endl;
            }
         }
      }
   }
   #endif
   if(false)//(depth==1)&&(rand()%10==0))
   {
      RecUnitCell ucm=uc0,ucp=uc0;
      for(unsigned int i=0;i<4;++i) {ucm.par[i]-=duc.par[i];ucp.par[i]+=duc.par[i];}
      for(unsigned int i=4;i<7;++i) {ucm.par[i]+=duc.par[i];ucp.par[i]-=duc.par[i];}
      vector<float> ucmd=ucm.DirectUnitCell();
      vector<float> ucpd=ucp.DirectUnitCell();
      char buf[200];
      sprintf(buf,"a=%5.2f-%5.2f b=%5.2f-%5.2f c=%5.2f-%5.2f alpha=%6.2f-%6.2f beta=%6.2f-%6.2f gamma=%6.2f-%6.2f V=%6.2f-%6.2f",
               ucpd[0],ucmd[0],ucpd[1],ucmd[1],ucpd[2],ucmd[2],ucpd[3]*RAD2DEG,ucmd[3]*RAD2DEG,
               ucpd[4]*RAD2DEG,ucmd[4]*RAD2DEG,ucpd[5]*RAD2DEG,ucmd[5]*RAD2DEG,ucpd[6],ucmd[6]);
      for(int i=0;i<depth;++i) cout<<" ";
      cout<<buf<<"level="<<depth<<", indexed="<<indexed<<"("<<mvSolution.size()<<" sol.)"<<endl;
   }
   if(localverbose)//(depth==1)&&(rand()%10==0))
   {
      RecUnitCell ucm=uc0,ucp=uc0;
      for(unsigned int i=0;i<4;++i) {ucm.par[i]-=duc.par[i];ucp.par[i]+=duc.par[i];}
      for(unsigned int i=4;i<7;++i) {ucm.par[i]+=duc.par[i];ucp.par[i]-=duc.par[i];}
      char buf[200];
      sprintf(buf,"a=%7.5f-%7.5f b=%7.5f-%7.5f c=%7.5f-%7.5f alpha=%7.5f-%7.5f beta=%7.5f-%7.5f gamma=%7.5f-%7.5f",
               ucp.par[1],ucm.par[1],ucp.par[2],ucm.par[2],ucp.par[3],ucm.par[3],ucp.par[4],ucm.par[4],
               ucp.par[5],ucm.par[5],ucp.par[6],ucm.par[6]);
      for(int i=0;i<depth;++i) cout<<" ";
      cout<<buf<<"level="<<depth<<", indexed="<<indexed<<"("<<mvSolution.size()<<" sol.)"<<endl;
   }
   nbCalc++;
   // :TODO: if we failed the dichotomy and reached some depth, try guessing a zero shift from the indexed reflections
   /*
   if((!indexed)&&(depth>=2))
   {
      vector<float> shifts(mpPeakList->GetPeakList().size());
      vector<PeakList::hkl>::const_iterator peakpos=mpPeakList->GetPeakList().begin();
      for(vector<float>::iterator spos=shifts.begin();spos!=shifts.end();)
      {   *spos++ = peakpos->d2diff * (float)(peakpos->isIndexed&&(!peakpos->isSpurious));peakpos++;}
      sort(shifts.begin(),shifts.end());
      uc0.par[0]=shifts[mpPeakList->GetPeakList().size()/2];//use median value
      indexed=DichoIndexed(*mpPeakList,uc0,duc,mNbSpurious);
      if(indexed) cout<<"Failed Dicho ? Trying auto-zero shifting :Worked !"<<endl;
   }
   */
   if(indexed)
   {
      unsigned int deeperSolutions=0;
      if(depth<mMaxDicVolDepth)
      {
         if(false)//depth>=3)
         {
            RecUnitCell ucm=uc0,ucp=uc0;
            for(unsigned int i=0;i<6;++i) {ucm.par[i]-=duc.par[i];ucp.par[i]+=duc.par[i];}
            vector<float> ucmd=ucm.DirectUnitCell();
            vector<float> ucpd=ucp.DirectUnitCell();
            char buf[200];
            sprintf(buf,"a=%5.2f-%5.2f b=%5.2f-%5.2f c=%5.2f-%5.2f alpha=%5.2f-%5.2f beta=%5.2f-%5.2f gamma=%5.2f-%5.2f V=%5.2f-%5.2f",
                     ucpd[0],ucmd[0],ucpd[1],ucmd[1],ucpd[2],ucmd[2],ucpd[3]*RAD2DEG,ucmd[3]*RAD2DEG,
                     ucpd[4]*RAD2DEG,ucmd[4]*RAD2DEG,ucpd[5]*RAD2DEG,ucmd[5]*RAD2DEG,ucpd[6],ucmd[6]);
            cout<<"Depth="<<depth<<buf<<endl;
         }
         RecUnitCell uc=uc0;
         // zero (if used...)
         duc.par[0]=0.5*duc.par[0];
         // Divide interval by 2, except if this parameter is already at a higher depth 
         // because a main axis has been indexed already.
         for(unsigned int i=1;i<mnpar;++i) duc.par[i]*=(0.5+0.5*(vdepth[i-1]>depth));
         
         for(int i0=-1;i0<=1;i0+=2)
         {
            //:TODO: dichotomy on zero shift ?
            if(localverbose) cout<<__FILE__<<":"<<__LINE__<<":"<<uc.par[3]<<" +/- "<<duc.par[3]<<" ("<<vdepth[2]<<")"<<endl;
            // Don't change parameter if it is already determined at a higher depth
            if(vdepth[0]==depth) {uc.par[1]=uc0.par[1]+i0*duc.par[1];}
            else {i0=2;}// no need to dicho this parameter which is already at higher depth
            if(mnpar==2)
               deeperSolutions+=RDicVol(uc,duc, depth+1,nbCalc,minV,maxV,vdepth);
            else
               for(int i1=-1;i1<=1;i1+=2)
               {
                  if(vdepth[1]==depth) {uc.par[2]=uc0.par[2]+i1*duc.par[2];}
                  else {i1=2;}// no need to dicho this parameter which is already at higher depth
                  if(mnpar==3)
                     deeperSolutions+=RDicVol(uc,duc, depth+1,nbCalc,minV,maxV,vdepth);
                  else
                     for(int i2=-1;i2<=1;i2+=2)
                     {
                        if(vdepth[2]==depth) {uc.par[3]=uc0.par[3]+i2*duc.par[3];}
                        else {i2=2;}// no need to dicho this parameter which is already at higher depth
                        if(mnpar==4)
                           deeperSolutions+=RDicVol(uc,duc, depth+1,nbCalc,minV,maxV,vdepth);
                        else
                           for(int i3=-1;i3<=1;i3+=2)
                           {
                              if(vdepth[3]==depth)uc.par[4]=uc0.par[4]+i3*duc.par[4];
                              else i3=2;
                              if(mnpar==5)
                                 deeperSolutions+=RDicVol(uc,duc, depth+1,nbCalc,minV,maxV,vdepth);
                              else
                                 for(int i4=-1;i4<=1;i4+=2)
                                 {
                                    uc.par[5]=uc0.par[5]+i4*duc.par[5];
                                    //if(mnpar==7)
                                    //   deeperSolutions+=RDicVol(uc,duc, depth+1,nbCalc,minV,maxV,vdepth);
                                    //else
                                       for(int i5=-1;i5<=1;i5+=2)
                                       {
                                          uc.par[6]=uc0.par[6]+i5*duc.par[6];
                                          //if(localverbose) cout<<__FILE__<<":"<<__LINE__<<":"<<uc.par[3]<<" +/- "<<duc.par[3]<<" ("<<vdepth[2]<<")"<<endl;
                                          deeperSolutions+=RDicVol(uc,duc, depth+1,nbCalc,minV,maxV,vdepth);
                                       }
                                 }
                           }
                     }
               }
          }
      }
      if((deeperSolutions==0) &&(depth>=mDicVolDepthReport))
      {
         mRecUnitCell=uc0;
         vector<float> uc=mRecUnitCell.DirectUnitCell();
         float score=Score(*mpPeakList,mRecUnitCell,mNbSpurious,false,true,false);
         // If we already have enough reports at higher depths (depth+2), don't bother record this one
         bool report=true;
         if(depth<(mMaxDicVolDepth-1))
            if(mvNbSolutionDepth[depth+2]>100)report=false;
         if(report && (((score>(mMinScoreReport*.5))&&(depth>=mDicVolDepthReport)) || (depth>=mMaxDicVolDepth)))
         {
            if(false)//score>) mBestScore//((score>mMinScoreReport)||(depth>=mDicVolDepthReport))
               cout<<__FILE__<<":"<<__LINE__<<" Depth="<<depth<<" (DIC) ! a="<<uc[0]<<", b="<<uc[1]<<", c="<<uc[2]<<", alpha="
                  <<uc[3]*RAD2DEG<<", beta="<<uc[4]*RAD2DEG<<", gamma="<<uc[5]*RAD2DEG<<", V="<<uc[6]
                  <<", score="<<score<<endl;
            this->LSQRefine(5,true,true);
            
            // Re-score (may change to a better hkl indexing), and refine again
            score=Score(*mpPeakList,mRecUnitCell,mNbSpurious,false,true,false);
            this->LSQRefine(5,true,true);
            
            uc=mRecUnitCell.DirectUnitCell();
            score=Score(*mpPeakList,mRecUnitCell,mNbSpurious,false,true,false);
            if(  ((score>mMinScoreReport)||(depth>=mDicVolDepthReport))
               &&((mvSolution.size()<50)||(score>(mBestScore/3)))
               &&((mvSolution.size()<50)||(score>mMinScoreReport)))
            {
               if((score>(mBestScore))||((score>(mBestScore*0.8))&&(mvSolution.size()<50)))//||(rand()%100==0))
               {
                  char buf[200];
                  {
                     RecUnitCell ucm=uc0,ucp=uc0;
                     for(unsigned int i=0;i<4;++i) {ucm.par[i]-=duc.par[i];ucp.par[i]+=duc.par[i];}
                     for(unsigned int i=4;i<7;++i) {ucm.par[i]+=duc.par[i];ucp.par[i]-=duc.par[i];}
                     vector<float> ucmd=ucm.DirectUnitCell();
                     vector<float> ucpd=ucp.DirectUnitCell();
                     sprintf(buf,"a=%5.2f-%5.2f b=%5.2f-%5.2f c=%5.2f-%5.2f alpha=%6.2f-%6.2f beta=%6.2f-%6.2f gamma=%6.2f-%6.2f V=%6.2f-%6.2f",
                              ucpd[0],ucmd[0],ucpd[1],ucmd[1],ucpd[2],ucmd[2],ucpd[3]*RAD2DEG,ucmd[3]*RAD2DEG,
                              ucpd[4]*RAD2DEG,ucmd[4]*RAD2DEG,ucpd[5]*RAD2DEG,ucmd[5]*RAD2DEG,ucpd[6],ucmd[6]);
                     for(int i=0;i<depth;++i) cout<<" ";
                     
                     cout<<buf<<"level="<<depth<<", indexed="<<indexed<<"("<<mvSolution.size()<<" sol.)"<<endl;
                     sprintf(buf,"a=%7.5f-%7.5f b=%7.5f-%7.5f c=%7.5f-%7.5f alpha=%7.5f-%7.5f beta=%7.5f-%7.5f gamma=%7.5f-%7.5f",
                              ucp.par[1],ucm.par[1],ucp.par[2],ucm.par[2],ucp.par[3],ucm.par[3],ucp.par[4],ucm.par[4],
                              ucp.par[5],ucm.par[5],ucp.par[6],ucm.par[6]);
                     for(int i=0;i<depth;++i) cout<<" ";
                     cout<<buf<<"level="<<depth<<", indexed="<<indexed<<"("<<mvSolution.size()<<" sol.)"<<endl;
                  }
                  sprintf(buf," Solution ? a=%7.3f b=%7.3f c=%7.3f alpha=%7.3f beta=%7.3f gamma=%7.3f V=%8.2f score=%6.2f #%4i",
                          uc[0],uc[1],uc[2],uc[3]*RAD2DEG,uc[4]*RAD2DEG,uc[5]*RAD2DEG,uc[6],score,mvSolution.size());
                  cout<<buf<<endl;
                  mBestScore=score;
               }
               mvSolution.push_back(make_pair(mRecUnitCell,score));
               mvNbSolutionDepth[depth]+=1;
               if((mvSolution.size()>1100)&&(rand()%1000==0))
               {
                  cout<<mvSolution.size()<<" solutions ! Reducing..."<<endl;
                  this->ReduceSolutions(true);// This will update the min report score
                  cout<<"-> "<<mvSolution.size()<<" remaining"<<endl;
               }
            }
         }
      }
      return deeperSolutions+1;
   }
   return 0;
}

void CellExplorer::DicVol(const float minScore,const unsigned int minDepth,const float stopOnScore,const unsigned int stopOnDepth)
{
   mDicVolDepthReport=minDepth;
   mMinScoreReport=minScore;
   this->Init();
   if(minDepth>mMaxDicVolDepth) mMaxDicVolDepth=minDepth;
   mvNbSolutionDepth.resize(mMaxDicVolDepth+1);
   for(unsigned int i=0;i<=mMaxDicVolDepth;++i) mvNbSolutionDepth[i]=0;
   
   float latstep=0.5,
         vstep=(mVolumeMax-mVolumeMin)/(ceil((mVolumeMax-mVolumeMin)/500)-0.0001);
   mCosAngMax=abs(cos(mAngleMax));
   const float cosangstep=mCosAngMax/(ceil(mCosAngMax/.08)-.0001);
   if(((mVolumeMax-mVolumeMin)/vstep)>10) vstep=(mVolumeMax-mVolumeMin)/9.999;
   if(((mLengthMax-mLengthMin)/latstep)>25) latstep=(mLengthMax-mLengthMin)/24.9999;
   
   cout<<mLengthMin<<"->"<<mLengthMax<<","<<latstep<<","<<(mLengthMax-mLengthMin)/latstep<<endl;
   cout<<mAngleMin<<"->"<<mAngleMax<<","<<cosangstep<<","<<mCosAngMax<<","<<(mAngleMax-mAngleMin)/cosangstep<<endl;
   cout<<mVolumeMin<<"->"<<mVolumeMax<<","<<vstep<<","<<(mVolumeMax-mVolumeMin)/vstep<<endl;
   RecUnitCell uc0,duc;
   uc0.mlattice=mlattice;
   duc.mlattice=mlattice;
   uc0.mCentering=mCentering;
   duc.mCentering=mCentering;
   //Zero shift parameter - not used for dicvol right now ? :TODO:
   uc0.par[0]=0.0;
   duc.par[0]=0.0;
   unsigned long nbCalc=0;
   Chronometer chrono;
   float bestscore=0;
   list<pair<RecUnitCell,float> >::iterator bestpos;
   bool breakDepth=false;
   // In the triclinic case, first try assigning a* and b* from the first reflections
   if(mlattice==TRICLINIC)
      for(float minv=mVolumeMin;minv<mVolumeMax;minv+=vstep)
      {
         float maxv=minv+vstep;
         if(maxv>mVolumeMax)maxv=mVolumeMax;
         cout<<"Starting: V="<<minv<<"->"<<maxv<<endl;
         const float minr=1/mLengthMax;
         const float maxr=1/mLengthMin;
         const float stepr=(maxr-minr)/24.999;
         float p1,p2;
         for(unsigned int i=0;i<=5;i++)
         {
            switch(i)
            {// Try to find a and b from the first observed reflections
               case 0: p1=mpPeakList->GetPeakList()[0].dobs  ;p2=mpPeakList->GetPeakList()[1].dobs  ; break;
               case 1: p1=mpPeakList->GetPeakList()[0].dobs  ;p2=mpPeakList->GetPeakList()[2].dobs  ; break;
               case 2: p1=mpPeakList->GetPeakList()[1].dobs/2;p2=mpPeakList->GetPeakList()[0].dobs  ; break;
               case 3: p1=mpPeakList->GetPeakList()[1].dobs/2;p2=mpPeakList->GetPeakList()[2].dobs  ; break;
               case 4: p1=mpPeakList->GetPeakList()[2].dobs/2;p2=mpPeakList->GetPeakList()[0].dobs  ; break;
               case 5: p1=mpPeakList->GetPeakList()[2].dobs/2;p2=mpPeakList->GetPeakList()[1].dobs  ; break;
            }
            //if(i>0) exit(0);
            if(p1>p2) continue;
            cout<<"Trying #"<<i<<": a*="<<p1<<", b*="<<p2<<endl;
            float min3r=p2,
                  max3r=maxr*2;// Due to the angles, use 2/mLengthMin
            const float step3r=(max3r-min3r)/(ceil((max3r-min3r)/stepr)-.001);
            vector<unsigned int> vdepth(mnpar-1);
            for(vector<unsigned int>::iterator pos=vdepth.begin();pos!=vdepth.end();) *pos++=0;
            vdepth[0]=3;
            vdepth[1]=3;
            for(float p3=min3r;p3<max3r;p3+=step3r)
            {//p3=c*
               //cout<<"    p3="<<p3<<endl;
               #ifdef TRI1
               float max4r=(p1+stepr)*(p1+stepr);
               const float step4r=max4r/(14.99);
               #else
               float max4r=(p1+stepr)/(2*p2);
               if(max4r>0.5) max4r=0.5;
               const float step4r=max4r/(ceil(max4r/.05)-.001);
               #endif
               for(float p4=-max4r;p4<max4r;p4+=step4r)
               {//p4=cos(gamma*)
                  //cout<<"      p4="<<p4<<endl;
                  #ifdef TRI1
                  float max5r=(p2+stepr)*(p2+stepr);
                  const float step5r=max5r/(14.99);
                  #else
                  float max5r=(p2+stepr)/(2*p3);
                  if(max5r>0.5) max5r=0.5;
                  const float step5r=max5r/(ceil(max5r/.05)-.001);
                  #endif
                  for(float p5=0;p5<max5r;p5+=step5r)
                  {//p5=cos(alpha*)
                     //cout<<"        p5="<<p5<<endl;
                     #ifdef TRI1
                     float max6r=(p1+stepr)*(p1+stepr);
                     const float step6r=max6r/(14.99);
                     #else
                     float max6r=(p1+stepr)/(2*p3);
                     if(max6r>0.5) max6r=0.5;
                     const float step6r=max6r/(ceil(max6r/.05)-.001);
                     #endif
                     for(float p6=0;p6<max6r;p6+=step6r)
                     {//p6=cos(beta*)
                        //cout<<"          p6="<<p6<<"/"<<p1<<"/"<<p3<<endl;
                        duc.par[1]=stepr*pow(float(0.51),int(vdepth[0]));
                        duc.par[2]=stepr*pow(float(0.51),int(vdepth[1]));
                        duc.par[3]=step3r*0.51;
                        duc.par[4]=step4r*0.51;
                        duc.par[5]=step5r*0.51;
                        duc.par[6]=step6r*0.51;
                        
                        uc0.par[0]=0;
                        uc0.par[1]=p1;
                        uc0.par[2]=p2;
                        uc0.par[3]=p3+step3r/2;
                        uc0.par[4]=p4+step4r/2;
                        uc0.par[5]=p5+step5r/2;
                        uc0.par[6]=p6+step6r/2;
                        //for(int i=0;i<=6;++i)cout<<uc0.par[i]<<",";
                        //cout<<endl;
                        //for(int i=0;i<=6;++i)cout<<duc.par[i]<<",";
                        //cout<<endl;
                        RDicVol(uc0,duc,0,nbCalc,minv,maxv,vdepth);
                     }
                  }
               }
            }
            cout<<"Finished trying: a*="<<p1<<" A, b*="<<p2<<" A, "<<nbCalc
               <<" unit cells tested, "<<nbCalc/chrono.seconds()<<" tests/s,   Elapsed time="
               <<chrono.seconds()<<"s, Best score="<<mBestScore<<", "<<stopOnScore<<", "<<breakDepth<<endl;
            breakDepth=false;
            if(stopOnDepth>0)
               for(unsigned int i=stopOnDepth; i<mvNbSolutionDepth.size();++i)
                  if(mvNbSolutionDepth[i]>1) {breakDepth=true;break;}
            if((mBestScore>stopOnScore)&&(breakDepth)) break;
         }//cases
         cout<<"Finished triclinic QUICK tests for: V="<<minv<<"->"<<maxv<<" A^3, "<<nbCalc
            <<" unit cells tested, "<<nbCalc/chrono.seconds()<<" tests/s,   Elapsed time="
            <<chrono.seconds()<<"s, Best score="<<mBestScore<<endl;
         if((mBestScore>stopOnScore)&&(breakDepth)) break;
      }//volume
   if((mBestScore<stopOnScore)||(!breakDepth))
   for(float minv=mVolumeMin;minv<mVolumeMax;minv+=vstep)
   {
      float maxv=minv+vstep;
      if(maxv>mVolumeMax)maxv=mVolumeMax;
      cout<<"Starting: V="<<minv<<"->"<<maxv<<endl;
      switch(mlattice)
      {
         case TRICLINIC:
         {
            const float minr=1/mLengthMax;
            const float maxr=1/mLengthMin*2;// *2 due to angles which can be up to cos(angle)=0.5
            const float stepr=(maxr-minr)/24.999;
            RecUnitCell uclarge,//ucm: smallest reciprocal, largest direct cell
                        ucsmall;//ucp: largest reciprocal, smallest direct cell
            vector<float> uclarged,ucsmalld;
            float min1r=minr;
            // The first peak cannot be at lower d* than 100
            float max1r=mpPeakList->GetPeakList()[0].dobs+stepr/2;
            const float step1r=(max1r-min1r)/(ceil((max1r-min1r)/stepr)-.001);
            for(float p1=min1r;p1<=max1r;p1+=step1r)
            {
               cout<<"p1="<<p1<<endl;
               float min2r=p1,max2r=maxr;
               const float step2r=(max2r-min2r)/(ceil((max2r-min2r)/stepr)-.001);
               for(float p2 =min2r;p2<max2r;p2+=step2r)
               {
                  // The second peak cannot be at lower d* than 200 or 010
                  if((mpPeakList->GetPeakList()[1].dobs<p2)&&((mpPeakList->GetPeakList()[1].dobs<(2*p1)))) break;
                  // The third peak cannot be at lower d* than 300 or 010
                  if((mpPeakList->GetPeakList()[2].dobs<p2)&&((mpPeakList->GetPeakList()[1].dobs<(3*p1)))) break;
                  cout<<"  p2="<<p2<<endl;
                  float min3r=p2,max3r=maxr;
                  const float step3r=(max3r-min3r)/(ceil((max3r-min3r)/stepr)-.001);
                  for(float p3=min3r;p3<max3r;p3+=step3r)
                  {
                     //cout<<"    p3="<<p3<<endl;
                     #ifdef TRI1
                     float max4r=(p1+step1r)*(p1+step1r);
                     const float step4r=max4r/(14.99);
                     #else
                     float max4r=(p1+stepr)/(2*p2);
                     const float step4r=max4r/(ceil(max4r/.05)-.001);
                     #endif
                     for(float p4=-max4r;p4<max4r;p4+=step4r)
                     {
                        //cout<<"      p4="<<p4<<endl;
                        #ifdef TRI1
                        float max5r=(p2+stepr)*(p2+step2r);
                        const float step5r=max5r/(14.99);
                        #else
                        float max5r=(p2+step2r)/(2*p3);
                        const float step5r=max5r/(ceil(max5r/.05)-.001);
                        #endif
                        for(float p5=0;p5<max5r;p5+=step5r)
                        {
                           //cout<<"        p5="<<p5<<endl;
                           #ifdef TRI1
                           float max6r=(p1+step1r)*(p1+step1r);
                           const float step6r=max6r/(14.99);
                           #else
                           float max6r=(p1+stepr)/(2*p3);
                           const float step6r=max6r/(ceil(max6r/.05)-.001);
                           #endif
                           for(float p6=0;p6<max6r;p6+=step6r)
                           {
                              //cout<<"          p6="<<p6<<endl;
                              duc.par[1]=step1r*0.51;
                              duc.par[2]=step2r*0.51;
                              duc.par[3]=step3r*0.51;
                              duc.par[4]=step4r*0.51;
                              duc.par[5]=step5r*0.51;
                              duc.par[6]=step6r*0.51;
                              
                              uc0.par[0]=0;
                              uc0.par[1]=p1+step1r/2;
                              uc0.par[2]=p2+step2r/2;
                              uc0.par[3]=p3+step3r/2;
                              uc0.par[4]=p4+step4r/2;
                              uc0.par[5]=p5+step5r/2;
                              uc0.par[6]=p6+step6r/2;
                              //for(int i=0;i<=6;++i)cout<<uc0.par[i]<<",";
                              //cout<<endl;
                              //for(int i=0;i<=6;++i)cout<<duc.par[i]<<",";
                              //cout<<endl;
                              RDicVol(uc0,duc,0,nbCalc,minv,maxv);
                           }
                        }
                     }
                  }
               }
            }
            break;
         }
         case MONOCLINIC:
         {
            RecUnitCell uclarge,//ucm: smallest reciprocal, largest direct cell
                        ucsmall;//ucp: largest reciprocal, smallest direct cell
            vector<float> uclarged,ucsmalld;
            latstep=(mLengthMax-mLengthMin)/24.999;
            for(float x4=0;x4<mCosAngMax+cosangstep;x4+=cosangstep)
            {
               const float sinbeta=sqrt(abs(1-x4*x4));
               float x1=mLengthMin;
               for(;x1<mLengthMax;x1+=latstep)
               {
                  float x2=mLengthMin;
                  for(;x2<mLengthMax;x2+=latstep)
                  {
                     float x3=x1;
                     const float x3step=(mLengthMax-x1)/(ceil((mLengthMax-x1)/latstep)-0.001);
                     for(;x3<mLengthMax;x3+=x3step) //x3+=(latstep+x3*sin4)
                     {
                        if((x3*x4)>x1) break;// | c * cos(beta) | <a
                        duc.par[1]=(1/(x1)-1/(x1+latstep))*0.5/sinbeta;
                        duc.par[2]=(1/(x2)-1/(x2+latstep))*0.5/sinbeta;
                        duc.par[3]=(1/(x3)-1/(x3+x3step ))*0.5/sinbeta;
                        duc.par[4]=cosangstep*0.5;
                        
                        uc0.par[0]=0;
                        uc0.par[1]=(1/(x1)+1/(x1+latstep))*0.5/sinbeta;
                        uc0.par[2]=(1/(x2)+1/(x2+latstep))*0.5/sinbeta;
                        uc0.par[3]=(1/(x3)+1/(x3+x3step ))*0.5/sinbeta;
                        uc0.par[4]=x4+cosangstep*.5;
                        
                        const float smallv=x1*x2*x3*sinbeta;
                        if(smallv>maxv) break;
                        const float largev=(x1+latstep)*(x2+latstep)*(x3+latstep)*(sinbeta+cosangstep);
                        if(largev<minv) continue;
                        /*
                        char buf[200];
                        sprintf(buf,"a=%5.2f-%5.2f b=%5.2f-%5.2f c=%5.2f-%5.2f alpha=%5.2f-%5.2f beta=%5.2f-%5.2f gamma=%5.2f-%5.2f V=%5.2f-%5.2f",
                                 ucsmalld[0],uclarged[0],ucsmalld[1],uclarged[1],ucsmalld[2],uclarged[2],ucsmalld[3]*RAD2DEG,uclarged[3]*RAD2DEG,
                                 ucsmalld[4]*RAD2DEG,uclarged[4]*RAD2DEG,ucsmalld[5]*RAD2DEG,uclarged[5]*RAD2DEG,ucsmalld[6],uclarged[6]);
                        cout<<buf<<"   VM="<<maxv<<", x3="<<x3<<endl;
                        */
                        RDicVol(uc0,duc,0,nbCalc,minv,maxv);
                     }//x3
                     //if(((ucsmalld[6]>maxv)&&(x3==x1))||(uclarged[1]>mLengthMax)) break;
                  }//x2
               }//x1
               // Test if we have one solution before going to the next angle range
               for(list<pair<RecUnitCell,float> >::iterator pos=mvSolution.begin();pos!=mvSolution.end();++pos)
               {
                  const float score=pos->second;//Score(*mpPeakList,pos->first,mNbSpurious);
                  if(score>bestscore) {bestscore=score;bestpos=pos;}
               }
               bool breakDepth=false;
               if(stopOnDepth>0)
                  for(unsigned int i=stopOnDepth; i<mvNbSolutionDepth.size();++i)
                     if(mvNbSolutionDepth[i]>1) {breakDepth=true;break;}
               if((bestscore>stopOnScore)&&(breakDepth)) break;
            }//x4
            break;
         }
         case ORTHOROMBIC:
         {
            if(false)
            {
               // Test 7.677350  5.803770  10.313160   V=480
               //const float a=7.75,b=5.75,c=10.25;
               // Test 6.062000 16.779400 16.8881 v=1750
               const float a=6.25,b=16.75,c=16.75;
               duc.par[1]=(1/(a-.25)-1/(a+.25))*0.5;
               duc.par[2]=(1/(b-.25)-1/(b+.25))*0.5;
               duc.par[3]=(1/(c-.25)-1/(c+.25))*0.5;
               uc0.par[0]=0;
               uc0.par[1]=1/a;
               uc0.par[2]=1/b;
               uc0.par[3]=1/c;
               RDicVol(uc0,duc,0,nbCalc,minv,maxv);
               break;
            }
            latstep=(mLengthMax-mLengthMin)/24.999;
            for(float x1=mLengthMin;x1<mLengthMax;x1+=latstep)
            {
               for(float x2=x1;x2<mLengthMax;x2+=latstep)
               {
                  for(float x3=x2;x3<mLengthMax;x3+=latstep)
                  {
                     duc.par[1]=(1/(x1)-1/(x1+latstep))*0.5;
                     duc.par[2]=(1/(x2)-1/(x2+latstep))*0.5;
                     duc.par[3]=(1/(x3)-1/(x3+latstep))*0.5;
                     
                     uc0.par[0]=0;
                     uc0.par[1]=(1/(x1)+1/(x1+latstep))*0.5;
                     uc0.par[2]=(1/(x2)+1/(x2+latstep))*0.5;
                     uc0.par[3]=(1/(x3)+1/(x3+latstep))*0.5;
                     
                     const float vmin=x1*x2*x3,vmax=(x1+latstep)*(x2+latstep)*(x3+latstep);
                     if(vmin>maxv) break;
                     if(vmax>=minv) RDicVol(uc0,duc,0,nbCalc,minv,maxv);
                  }
                  if((x1*x2*x2)>maxv) break;
               }
               if((x1*x1*x1)>maxv) break;
            }
            break;
         }
         case HEXAGONAL:
         {
            vector<float> uclarged,ucsmalld;// Small & large UC in direct space
            latstep=(mLengthMax-mLengthMin)/24.999;
            for(float x1=mLengthMin;;x1+=latstep)
            {
               for(float x2=mLengthMin;x2<(mLengthMax+latstep);x2+=latstep)
               {
                  duc.par[1]=(1/(x1)-1/(x1+latstep))*0.5;
                  duc.par[2]=(1/(x2)-1/(x2+latstep))*0.5;
                     
                  uc0.par[0]=0;
                  uc0.par[1]=(1/(x1)+1/(x1+latstep))*0.5;
                  uc0.par[2]=(1/(x2)+1/(x2+latstep))*0.5;
                  
                  RecUnitCell uclarge=uc0,ucsmall=uc0;
                  for(unsigned int i=0;i<6;++i) {uclarge.par[i]-=duc.par[i];ucsmall.par[i]+=duc.par[i];}
                  uclarged=uclarge.DirectUnitCell();
                  ucsmalld=ucsmall.DirectUnitCell();
                  /*
                  char buf[200];
                  sprintf(buf,"a=%5.2f-%5.2f b=%5.2f-%5.2f c=%5.2f-%5.2f alpha=%5.2f-%5.2f beta=%5.2f-%5.2f gamma=%5.2f-%5.2f V=%5.2f-%5.2f",
                           ucsmalld[0],uclarged[0],ucsmalld[1],uclarged[1],ucsmalld[2],uclarged[2],ucsmalld[3]*RAD2DEG,uclarged[3]*RAD2DEG,
                           ucsmalld[4]*RAD2DEG,uclarged[4]*RAD2DEG,ucsmalld[5]*RAD2DEG,uclarged[5]*RAD2DEG,ucsmalld[6],uclarged[6]);
                  */
                  if((ucsmalld[6]<maxv)&&(uclarged[6]>minv))
                  {
                     //cout<<buf<<endl;
                     RDicVol(uc0,duc,0,nbCalc,minv,maxv);
                  }
                  //else cout<<buf<<" BREAK"<<endl;
               }
               if(uclarged[0]>mLengthMax) break;
            }
            break;
         }
         case RHOMBOEDRAL:   //:TODO:
         {
            latstep=(mLengthMax-mLengthMin)/24.999;
            for(float x1=mLengthMin;x1<(mLengthMax+latstep);x1+=latstep)
            {
               for(float x2=0;x2<mCosAngMax+cosangstep;x2+=cosangstep)
               {
                  duc.par[1]=latstep/2*1.1;
                  duc.par[2]=cosangstep/2*1.1;
                  
                  uc0.par[0]=0;
                  uc0.par[1]=x1-latstep/2*1.1;
                  uc0.par[2]=x2-cosangstep/2*1.1;
                  vector<float> uc=uc0.DirectUnitCell();
                  if((uc[6]<maxv)&&(uc[6]>minv))
                  {
                     RDicVol(uc0,duc,0,nbCalc,minv,maxv);
                  }
               }
            }
            break;
         }
         case TETRAGONAL:
         {
            vector<float> uclarged,ucsmalld;// Small & large UC in direct space
            latstep=(mLengthMax-mLengthMin)/24.999;
            for(float x1=mLengthMin;x1<mLengthMax;x1+=latstep)
            {
               for(float x2=mLengthMin;x2<mLengthMax;x2+=latstep)
               {
                  duc.par[1]=(1/(x1)-1/(x1+latstep))*0.5;
                  duc.par[2]=(1/(x2)-1/(x2+latstep))*0.5;
                     
                  uc0.par[0]=0;
                  uc0.par[1]=(1/(x1)+1/(x1+latstep))*0.5;
                  uc0.par[2]=(1/(x2)+1/(x2+latstep))*0.5;
                  
                  RecUnitCell uclarge=uc0,ucsmall=uc0;
                  for(unsigned int i=0;i<6;++i) {uclarge.par[i]-=duc.par[i];ucsmall.par[i]+=duc.par[i];}
                  uclarged=uclarge.DirectUnitCell();
                  ucsmalld=ucsmall.DirectUnitCell();
                  /*
                  char buf[200];
                  sprintf(buf,"a=%5.2f-%5.2f b=%5.2f-%5.2f c=%5.2f-%5.2f alpha=%5.2f-%5.2f beta=%5.2f-%5.2f gamma=%5.2f-%5.2f V=%5.2f-%5.2f",
                           ucsmalld[0],uclarged[0],ucsmalld[1],uclarged[1],ucsmalld[2],uclarged[2],ucsmalld[3]*RAD2DEG,uclarged[3]*RAD2DEG,
                           ucsmalld[4]*RAD2DEG,uclarged[4]*RAD2DEG,ucsmalld[5]*RAD2DEG,uclarged[5]*RAD2DEG,ucsmalld[6],uclarged[6]);
                  */
                  if((ucsmalld[6]<maxv)&&(uclarged[6]>minv))
                  {
                     RDicVol(uc0,duc,0,nbCalc,minv,maxv);
                  }
                  if(ucsmalld[6]>maxv) break;
               }
               if((x1*mLengthMin*mLengthMin)>maxv) break;
            }
            break;
         }
         case CUBIC:
         {
            latstep=(mLengthMax-mLengthMin)/24.999;
            cout<<mLengthMax<<","<<mLengthMin<<","<<latstep<<endl;
            for(float x1=mLengthMin;x1<(mLengthMax+latstep);x1+=latstep)
            {
               duc.par[1]=(1/(x1)-1/(x1+latstep))*0.5;
                  
               uc0.par[0]=0;
               uc0.par[1]=(1/(x1)+1/(x1+latstep))*0.5;
               
               const float vmin=x1*x1*x1,vmax=(x1+latstep)*(x1+latstep)*(x1+latstep);
               if(vmin>maxv)break;
               if(vmax>minv) RDicVol(uc0,duc,0,nbCalc,minv,maxv);
            }
            break;
         }
      }
      cout<<"Finished: V="<<minv<<"->"<<maxv<<" A^3, "<<nbCalc
          <<" unit cells tested, "<<nbCalc/chrono.seconds()<<" tests/s,   Elapsed time="
          <<chrono.seconds()<<"s"<<endl;
      for(list<pair<RecUnitCell,float> >::iterator pos=mvSolution.begin();pos!=mvSolution.end();++pos)
      {
         const float score=pos->second;//Score(*mpPeakList,pos->first,mNbSpurious);
         if(score>bestscore) {bestscore=score;bestpos=pos;}
      }
      bool breakDepth=false;
      if(stopOnDepth>0)
         for(unsigned int i=stopOnDepth; i<mvNbSolutionDepth.size();++i)
            if(mvNbSolutionDepth[i]>1) {breakDepth=true;break;}
      if((bestscore>stopOnScore)&&(breakDepth)) break;
   }
   /*
   {// Tag spurious lines
      vector<int> vSpuriousScore;
      for(vector<PeakList::hkl>::const_iterator pos=mpPeakList->GetPeakList().begin();pos!=mpPeakList->GetPeakList().end();++pos)
         vSpuriousScore.push_back(pos->stats);
      sort(vSpuriousScore.begin(),vSpuriousScore.end());
      const int threshold=vSpuriousScore[vSpuriousScore.size()/2]*5;
      for(vector<PeakList::hkl>::iterator pos=mpPeakList->mvHKL.begin();pos!=mpPeakList->mvHKL.end();++pos)
         if(pos->stats > threshold) pos->isSpurious=true;
         else pos->isSpurious=false;
      mpPeakList->Print(cout);
   }
   */
   this->ReduceSolutions(true);
   bestscore=0;bestpos=mvSolution.end();
   for(list<pair<RecUnitCell,float> >::iterator pos=mvSolution.begin();pos!=mvSolution.end();++pos)
   {
      const float score=Score(*mpPeakList,pos->first,mNbSpurious);
      if(score>bestscore) {bestpos=pos;bestscore=score;}
      vector<float> uc=pos->first.DirectUnitCell();
      cout<<__FILE__<<":"<<__LINE__<<" Solution ? a="<<uc[0]<<", b="<<uc[1]<<", c="<<uc[2]
          <<", alpha="<<uc[3]*RAD2DEG<<", beta="<<uc[4]*RAD2DEG<<", gamma="<<uc[5]*RAD2DEG
          <<", V="<<uc[6]<<", score="<<score<<endl;
   }
   if(bestpos!=mvSolution.end())
   {
      vector<float> uc=bestpos->first.DirectUnitCell();
      cout<<__FILE__<<":"<<__LINE__<<" BEST ? a="<<uc[0]<<", b="<<uc[1]<<", c="<<uc[2]
            <<", alpha="<<uc[3]*RAD2DEG<<", beta="<<uc[4]*RAD2DEG<<", gamma="<<uc[5]*RAD2DEG
            <<", V="<<uc[6]<<", score="<<bestscore<<endl;
      cout<<nbCalc<<"unit cells tested, "<<nbCalc/chrono.seconds()<<" tests/s,   Elapsed time="
          <<chrono.seconds()<<"s"<<endl;
   }
}

bool SimilarRUC(const RecUnitCell &c0,const RecUnitCell &c1, const float delta=0.005)
{
   vector<float> uc0=c0.DirectUnitCell();
   vector<float> uc1=c1.DirectUnitCell();
   float diff=0;
   for(unsigned int i=0;i<6;++i) diff += abs(uc0[i]-uc1[i]);
   return (diff/6)<delta;
}

bool compareRUCScore(std::pair<RecUnitCell,float> &p1, std::pair<RecUnitCell,float> &p2)
{
   return p1.second > p2.second;
}

void CellExplorer::ReduceSolutions(const bool updateReportThreshold)
{
   const bool verbose=false;
   std::list<std::pair<RecUnitCell,float> > vSolution2;
   
   // keep only solutions above mBestScore/5
   for(list<pair<RecUnitCell,float> >::iterator pos=mvSolution.begin();pos!=mvSolution.end();)
   {
      if(pos->second<(mBestScore/5)) pos=mvSolution.erase(pos);
      else ++pos;
   }
   if(updateReportThreshold&& ((mBestScore/5)>mMinScoreReport)) 
   {
      cout<<"CellExplorer::ReduceSolutions(): update threshold for report from "
          <<mMinScoreReport<<" to "<<mBestScore/5<<endl;
      mMinScoreReport=mBestScore/5;
   }
   
   while(mvSolution.size()>0)
   {
      vSolution2.push_back(mvSolution.front());
      mvSolution.pop_front();
      vector<float> uc=vSolution2.back().first.DirectUnitCell();
      if(verbose)
         cout<<__FILE__<<":"<<__LINE__<<" SOLUTION: a="<<uc[0]<<", b="<<uc[1]<<", c="<<uc[2]
               <<", alpha="<<uc[3]*RAD2DEG<<", beta="<<uc[4]*RAD2DEG<<", gamma="<<uc[5]*RAD2DEG
               <<", V="<<uc[6]<<", score="<<vSolution2.back().second<<",   SIMILAR TO:"<<endl;
      for(list<pair<RecUnitCell,float> >::iterator pos=mvSolution.begin();pos!=mvSolution.end();)
      {
         if(SimilarRUC(pos->first,vSolution2.back().first))
         {
            uc=pos->first.DirectUnitCell();
            if(verbose)
               cout<<__FILE__<<":"<<__LINE__<<"        1: a="<<uc[0]<<", b="<<uc[1]<<", c="<<uc[2]
                     <<", alpha="<<uc[3]*RAD2DEG<<", beta="<<uc[4]*RAD2DEG<<", gamma="<<uc[5]*RAD2DEG
                     <<", V="<<uc[6]<<", score="<<pos->second<<"       ("<<mvSolution.size()<<")"<<endl;
            if(vSolution2.back().first.mlattice==pos->first.mlattice)
            {
               if(pos->second>vSolution2.back().second) vSolution2.back()=*pos;
            }
            else if(vSolution2.back().first.mlattice<pos->first.mlattice) vSolution2.back()=*pos;
            pos=mvSolution.erase(pos);
         }
         else
         {
            uc=pos->first.DirectUnitCell();
            if(verbose)
               cout<<__FILE__<<":"<<__LINE__<<"        0: a="<<uc[0]<<", b="<<uc[1]<<", c="<<uc[2]
                     <<", alpha="<<uc[3]*RAD2DEG<<", beta="<<uc[4]*RAD2DEG<<", gamma="<<uc[5]*RAD2DEG
                     <<", V="<<uc[6]<<", score="<<pos->second<<"       ("<<mvSolution.size()<<")"<<endl;
            ++pos;
         }
      }
   }
   mvSolution=vSolution2;
   mvSolution.sort(compareRUCScore);
   
   // keep at most 100 solutions, update mDicVolDepthReport and mMinScoreReport if necessary
   if(mvSolution.size()>100)
   {
      mvSolution.resize(100);
      if(updateReportThreshold && (mvSolution.back().second>mMinScoreReport))
      {
         cout<<"CellExplorer::ReduceSolutions(): update threshold for report from "
             <<mMinScoreReport<<" to "<<mvSolution.back().second<<endl;
         mMinScoreReport=mvSolution.back().second;
      }
   }
}


void CellExplorer::Init()
{
   // Prepare global optimisation
   //for(unsigned int i=0;i<mpPeakList->nb;++i)
   //   cout<<__FILE__<<":"<<__LINE__<<":d*="<<mpPeakList->mvdobs[i]<<", d*^2="<<mpPeakList->mvd2obs[i]<<endl;
   srand(time(NULL));
   vector<pair<RecUnitCell,float> >::iterator pos;
   const float min_latt=1./mLengthMax;
   const float max_latt=1./mLengthMin;
   const float amp_crossp=abs(cos(mAngleMax));
   //mMin[0]=-.002;mAmp[0]=.004;
   mMin[0]=.00;mAmp[0]=.00;
   switch(mlattice)
   {
      case TRICLINIC:
         mMin[1]=min_latt;mAmp[1]=max_latt-min_latt;
         mMin[2]=min_latt;mAmp[2]=max_latt-min_latt;
         mMin[3]=min_latt;mAmp[3]=max_latt-min_latt;
         mMin[4]=0;mAmp[4]=amp_crossp;
         mMin[5]=0;mAmp[5]=amp_crossp;
         mMin[6]=0;mAmp[6]=amp_crossp;
         mnpar=7;
	 break;
      case MONOCLINIC:
         mMin[1]=min_latt;mAmp[1]=max_latt-min_latt;
         mMin[2]=min_latt;mAmp[2]=max_latt-min_latt;
         mMin[3]=min_latt;mAmp[3]=max_latt-min_latt;
         mMin[4]=0;mAmp[4]=amp_crossp;
	 mnpar=5;
         break;
      case ORTHOROMBIC:
         mMin[1]=min_latt;mAmp[1]=max_latt-min_latt;
         mMin[2]=min_latt;mAmp[2]=max_latt-min_latt;
         mMin[3]=min_latt;mAmp[3]=max_latt-min_latt;
	 mnpar=4;
         break;
      case HEXAGONAL:
         mMin[1]=min_latt;mAmp[1]=max_latt-min_latt;
         mMin[2]=min_latt;mAmp[2]=max_latt-min_latt;
         mnpar=3;
         break;
      case RHOMBOEDRAL:
         mMin[1]=min_latt;mAmp[1]=max_latt-min_latt;
         mMin[2]=-amp_crossp;mAmp[2]=2*amp_crossp;
         mnpar=3;
         break;
      case TETRAGONAL:
         mMin[1]=min_latt;mAmp[1]=max_latt-min_latt;
         mMin[2]=min_latt;mAmp[2]=max_latt-min_latt;
         mnpar=3;
         break;
      case CUBIC:
         mMin[1]=min_latt;mAmp[1]=max_latt-min_latt;
         mnpar=2;
         break;
   }
   //for(unsigned int k=0;k<mnpar;++k) cout<<"par["<<k<<"]: "<<mMin[k]<<"->"<<mMin[k]+mAmp[k]<<endl;
   
   unsigned int nb1=0,nb2=0;
   switch(mlattice)
   {
      case TRICLINIC:
      {
         nb1=3;
         nb2=3;
         break;
      }
      case MONOCLINIC:
      {
         nb1=3;
         nb2=1;
         break;
      }
      case ORTHOROMBIC:
      {
         nb1=3;
         break;
      }
      case HEXAGONAL:
      {
         nb1=2;
         break;
      }
      case RHOMBOEDRAL:
      {
         nb1=2;
         break;
      }
      case TETRAGONAL:
      {
         nb1=2;
         break;
      }
      case CUBIC:
      {
         nb1=1;
         break;
      }
   }
   this->ResetParList();
   {
      RefinablePar tmp("Zero",mRecUnitCell.par+0,-0.01,0.01,gpRefParTypeObjCryst,REFPAR_DERIV_STEP_ABSOLUTE,
                       true,false,true,false);
      tmp.SetDerivStep(1e-4);
      this->AddPar(tmp);
   }
   char buf [50];
   string str="Reciprocal unit cell par #";
   for(unsigned int i=0;i<nb1;++i)
   {
      sprintf(buf,"%i",i);
      RefinablePar tmp(str+(string)buf,mRecUnitCell.par+i+1,
                       0.01,1.,gpRefParTypeObjCryst,REFPAR_DERIV_STEP_ABSOLUTE,
                       false,false,true,false);
      tmp.SetDerivStep(1e-4);
      this->AddPar(tmp);
   }
   for(unsigned int i=nb1;i<(nb1+nb2);++i)
   {
      sprintf(buf,"%i",i);
      RefinablePar tmp(str+(string)buf,mRecUnitCell.par+i+1,
                       0.0,0.5,gpRefParTypeObjCryst,REFPAR_DERIV_STEP_ABSOLUTE,
                       false,false,true,false);
      tmp.SetDerivStep(1e-4);
      this->AddPar(tmp);
   }
}

}//namespace
