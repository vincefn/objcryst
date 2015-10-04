/*  ObjCryst++ Object-Oriented Crystallographic Library
    (c) 2007-2008 Vincent Favre-Nicolin vincefn@users.sourceforge.net

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
/*   PDF.cpp - Pair Distribution calculations
*
*/
#include "ObjCryst/ObjCryst/PDF.h"
#include "ObjCryst/Quirks/VFNStreamFormat.h"
#ifdef __WX__CRYST__
   #include "ObjCryst/wxCryst/wxPDF.h"
#endif

#include <string>

namespace ObjCryst
{
const RefParType *gpRefParTypePDF=0;

PDF::PDF():
mRadiationType(RAD_XRAY)
{
}

PDF::PDF(const PDF &old):
mRadiationType(old.mRadiationType),
mPDFR(old.mPDFR),mPDFCalc(old.mPDFCalc),mPDFObs(old.mPDFObs),
mvPDFPhase(old.mvPDFPhase)
{}

PDF::~PDF(){}

const string& PDF::GetClassName() const
{
    static const string className = "PDF";
    return className;
}

void PDF::SetPDFObs(const CrystVector_REAL &r,const CrystVector_REAL &obs)
{
   mPDFR=r;
   mPDFObs=obs,
   mPDFCalc.resize(r.numElements());
   mPDFCalc=0;
}

RadiationType PDF::GetRadiationType()const{return mRadiationType;}

void PDF::SetRadiationType(RadiationType type){mRadiationType=type;}

REAL PDF::GetRMax()const
{
   const unsigned long nb=mPDFR.numElements();
   if(nb>0) return mPDFR(nb-1);
   return 0;
}

const CrystVector_REAL &PDF::GetPDFCalc()const
{
   //:TODO: re-calculate only when needed !
   mPDFCalc.resize(mPDFR.numElements());
   mPDFCalc=0;
   CrystVector_REAL tmp;
   for(list<pair<PDFPhase*,REAL> >::const_iterator pos=mvPDFPhase.begin();pos!=mvPDFPhase.end();++pos)
   {
      tmp=pos->first->GetPDFCalc();
      tmp*=pos->second;
      mPDFCalc+=tmp;
   }
   return mPDFCalc;
}

const CrystVector_REAL &PDF::GetPDFObs()const
{
   return mPDFObs;
}

const CrystVector_REAL &PDF::GetPDFR()const
{
   return mPDFR;
}

void PDF::AddPDFPhase(PDFPhase &phase)
{
   mvPDFPhase.push_back(make_pair(&phase,1.0));
}

#ifdef __WX__CRYST__
WXCrystObjBasic* PDF::WXCreate(wxWindow* parent)
{
   if(mpWXCrystObj==0) mpWXCrystObj=new WXPDF(parent,this);
   return mpWXCrystObj;
}
#endif


////////////////////////// PDFPhase /////////////////////////////
PDFPhase::PDFPhase(const PDF &pdf):
mpPDF(&pdf)
{
}

const CrystVector_REAL &PDFPhase::GetPDFCalc()const
{
   this->CalcPDF();
   return mPDFCalc;
}

////////////////////////// PDFCrystal /////////////////////////////

PDFCrystal::PDFCrystal(const PDF &pdf, const Crystal &cryst):
PDFPhase(pdf),mpCrystal(&cryst),mDelta1(0.0),mDelta2(0.0),mQbroad(0.0),mQdamp(0.0)
{
}

PDFCrystal::pdfAtom::pdfAtom():
fx0(0),fy0(0),fz0(0),x0(0),y0(0),z0(0),occupBi(1),hasChanged(true)
{}

void PDFCrystal::Init(const PDF &pdf, const Crystal &cryst)
{
   mpCrystal=&cryst;
   mpPDF=&pdf;
   this->ResetParList();
   {
      RefinablePar tmp("Delta1",&mDelta1,0,1.0,gpRefParTypePDF,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,1);
      tmp.AssignClock(mClockMaster);
      tmp.SetDerivStep(1e-3);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("Delta2",&mDelta2,0,1.0,gpRefParTypePDF,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,1);
      tmp.AssignClock(mClockMaster);
      tmp.SetDerivStep(1e-3);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("Qbroad",&mQbroad,0,1.0,gpRefParTypePDF,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,1);
      tmp.AssignClock(mClockMaster);
      tmp.SetDerivStep(1e-3);
      this->AddPar(tmp);
   }
   {
      RefinablePar tmp("Qdamp",&mQdamp,0,1.0,gpRefParTypePDF,
                        REFPAR_DERIV_STEP_ABSOLUTE,true,true,true,false,1);
      tmp.AssignClock(mClockMaster);
      tmp.SetDerivStep(1e-3);
      this->AddPar(tmp);
   }
}

void PDFCrystal::CalcPDF()const
{
   const unsigned long nbr=mpPDF->GetPDFR().numElements();
   mPDFCalc.resize(nbr);
   if(mpCrystal==0)
   {
      mPDFCalc=0;
      return;
   }
   // Get current fractionnal coordinates from Crystal, check if any has changed
   const ScatteringComponentList *pScatt=&(mpCrystal->GetScatteringComponentList());
   unsigned long nb=pScatt->GetNbComponent();
   pScatt->Print();

   // Calc <b> and rho0
   REAL rho0=0,b_av=0;

   // Note: We cannot use the dynamical occupancy as computed in a Crystal Object,
   //as we need the real occupancy for each *unique* poisition.
   // :TODO: So we need to compute a new dynamical occupancy that only corrects the overlap
   //between one unique atom and different atoms, excluding symetrics of the unique atom.

   if(nb!=mvPDFAtom.size())
   {
      mvPDFAtom.resize(nb);
      unsigned int i=0;
      for(vector<pdfAtom>::iterator pos=mvPDFAtom.begin();pos!=mvPDFAtom.end();++pos)
      {
         pos->fx0=(*pScatt)(i).mX;
         pos->fy0=(*pScatt)(i).mY;
         pos->fz0=(*pScatt)(i).mZ;
         if((*pScatt)(i).mpScattPow==0) pos->occupBi=0;
         else
         {
            const REAL occ=(*pScatt)(i).mOccupancy;
            const REAL b=(*pScatt)(i).mpScattPow->GetForwardScatteringFactor(mpPDF->GetRadiationType());
            rho0+=occ;
            b_av+=occ*b;
            pos->occupBi= occ*b;
            pos->pScattPow=(*pScatt)(i).mpScattPow;
            pos->hasChanged=true;
         }
         ++i;
      }
   }
   else
   {
      unsigned int i=0;
      for(vector<pdfAtom>::iterator pos=mvPDFAtom.begin();pos!=mvPDFAtom.end();++pos)
      {
         if(pos->fx0!=(*pScatt)(i).mX){pos->fx0=(*pScatt)(i).mX;pos->hasChanged=true;}
         if(pos->fy0!=(*pScatt)(i).mY){pos->fy0=(*pScatt)(i).mY;pos->hasChanged=true;}
         if(pos->fz0!=(*pScatt)(i).mZ){pos->fz0=(*pScatt)(i).mZ;pos->hasChanged=true;}
         if(pos->pScattPow!=(*pScatt)(i).mpScattPow){pos->pScattPow=(*pScatt)(i).mpScattPow;pos->hasChanged=true;}
         REAL occupBi;
         if((*pScatt)(i).mpScattPow==0) occupBi=0;
         else
         {
            const REAL occ=(*pScatt)(i).mOccupancy;
            const REAL b=(*pScatt)(i).mpScattPow->GetForwardScatteringFactor(mpPDF->GetRadiationType());
            rho0+=occ;
            b_av+=occ*b;
            occupBi= occ*b;
            //:TODO: check if B-factor has changed
         }
         if(pos->occupBi!=occupBi){pos->occupBi=occupBi;pos->hasChanged=true;}
      }
   }

   // Determine number of translations needed
   //to get all interatomic distances up to Rmax
   // :TODO: faster by smarter limits on translations ?
   const int nx=int(ceil(mpPDF->GetRMax()/mpCrystal->GetLatticePar(0))+.001);
   const int ny=int(ceil(mpPDF->GetRMax()/mpCrystal->GetLatticePar(1))+.001);
   const int nz=int(ceil(mpPDF->GetRMax()/mpCrystal->GetLatticePar(2))+.001);
   const unsigned int nbSymmetrics=mpCrystal->GetSpaceGroup().GetNbSymmetrics();
   const unsigned int neq=(2*nx+1)*(2*ny+1)*(2*nz+1)*nbSymmetrics;

   b_av/=rho0;

   cout<<"rho0="<<rho0<<"*"<<nbSymmetrics<<"="<<rho0*nbSymmetrics
       <<"/"<<mpCrystal->GetVolume()<<"="<<rho0*nbSymmetrics/mpCrystal->GetVolume()<<endl;

   rho0*=nbSymmetrics/mpCrystal->GetVolume();

   // Calc all equivalent positions & translations
   // TODO: Use knowledge of special positions, rather than use dynamical occupancy ?
   {
      const CrystMatrix_REAL *pOrth=&(mpCrystal->GetOrthMatrix());
      const REAL m00=(*pOrth)(0,0);
      const REAL m01=(*pOrth)(0,1);
      const REAL m02=(*pOrth)(0,2);
      const REAL m11=(*pOrth)(1,1);
      const REAL m12=(*pOrth)(1,2);
      const REAL m22=(*pOrth)(2,2);
      cout<<m00<<" "<<m01<<" "<<m02<<endl<<m11<<" "<<m12<<endl<<m22<<endl<<neq<<endl;
      CrystMatrix_REAL symmetricsCoords;
      for(vector<pdfAtom>::iterator pos=mvPDFAtom.begin();pos!=mvPDFAtom.end();++pos)
      {
         if(pos->hasChanged)
         {
            symmetricsCoords=mpCrystal->GetSpaceGroup().GetAllSymmetrics(pos->fx0,pos->fy0,pos->fz0);
            pos->x.resize(neq);
            pos->y.resize(neq);
            pos->z.resize(neq);
            REAL *px=pos->x.data();
            REAL *py=pos->y.data();
            REAL *pz=pos->z.data();
            for(unsigned int j=0;j<nbSymmetrics;++j)
            {
               REAL x=fmod(symmetricsCoords(j,0),(REAL)1.0);
               REAL y=fmod(symmetricsCoords(j,1),(REAL)1.0);
               REAL z=fmod(symmetricsCoords(j,2),(REAL)1.0);
               mpCrystal->FractionalToOrthonormalCoords(x,y,z);
               if(pos==mvPDFAtom.begin()) cout<<j<<": "<<x<<","<<y<<","<<z<<endl;
               if(j==0){pos->x0=x;pos->y0=y;pos->z0=z;}
               for(int ix=-nx;ix<=nx;++ix)
                  for(int iy=-ny;iy<=ny;++iy)
                     for(int iz=-nz;iz<=nz;++iz)
                     {
                        //cout<<ix<<","<<iy<<","<<iz<<endl;
                        *px++ = x + ix*m00 + iy*m01 + iz*m02;
                        *py++ = y          + iy*m11 + iz*m12;
                        *pz++ = z                   + iz*m22;
                     }
            }
         }
      }
   }
   // Calculate pdf
   // :TODO: only recalculate the contributions that have changed !
   mPDFCalc=0;
   const REAL r2max=(mpPDF->GetRMax()+0.2)*(mpPDF->GetRMax()+0.2);
   const REAL nsigcut=5;// Cut gaussian at abs(r_ij-r)<3*sigma
   const REAL norm=1/sqrt(2*M_PI);
   for(vector<pdfAtom>::iterator pos=mvPDFAtom.begin();pos!=mvPDFAtom.end();++pos)
   {
      if(pos->occupBi==0) continue;
      const REAL x0=pos->x0;
      const REAL y0=pos->y0;
      const REAL z0=pos->z0;
      for(vector<pdfAtom>::iterator pos1=pos;pos1!=mvPDFAtom.end();++pos1)
      {
         if(pos1->occupBi==0) continue;
         REAL normij=norm*pos->occupBi*pos1->occupBi/(b_av*b_av)/nb;
         if(pos!=pos1)normij*=2;// i!j should be counted twice in the loop
         const REAL *px=pos1->x.data();
         const REAL *py=pos1->y.data();
         const REAL *pz=pos1->z.data();
         const REAL sigma2=(pos->pScattPow->GetBiso()+pos1->pScattPow->GetBiso())/(8*M_PI*M_PI);
         for(unsigned long i=0;i<neq;++i)
         {
            const REAL dx=*px++-x0, dy=*py++-y0, dz=*pz++-z0;
            const REAL d2=dx*dx+dy*dy+dz*dz;
            if((d2<r2max)&&(d2>1))
            {
               const REAL rij=sqrt(d2);
               REAL s2=sigma2*(1-mDelta1/rij-mDelta2/d2+mQbroad*d2);
               if(s2<.01) s2=0.01;
               const REAL sig=sqrt(s2);
               REAL *p=mPDFCalc.data();
               const REAL *pr=mpPDF->GetPDFR().data();
               const REAL rmin=rij-nsigcut*sig,rmax=rij+nsigcut*sig;
               //cout<<"    "<<rij<<":"<<rmin<<"->"<<rmax<<endl;
               const REAL n=normij/sig*exp(-0.5*rij*mQdamp*mQdamp);
               for(unsigned long i=0;i<nbr;++i)
               {
                  if(*pr<rmin){pr++;p++;continue;}
                  if(*pr>rmax) break;
                  const REAL dr=rij-*pr;
                  *p += n*exp(-dr*dr/(2*s2));
                  p++;pr++;
               }
            }
         }
      }
   }
   cout<<rho0<<","<<b_av<<endl;
   mPDFCalc/=mpPDF->GetPDFR();

   CrystVector_REAL tmp;
   tmp=mpPDF->GetPDFR();
   tmp*=4*M_PI*rho0;
   mPDFCalc-=tmp;
}

#ifdef __WX__CRYST__
WXCrystObjBasic* PDFCrystal::WXCreate(wxWindow* parent)
{
   if(mpWXCrystObj==0) mpWXCrystObj=new WXPDFCrystal(parent,this);
   return mpWXCrystObj;
}
#endif

}//namespace
