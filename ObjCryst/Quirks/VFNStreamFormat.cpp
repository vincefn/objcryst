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
// This file includes some "effectors" used to format
//strings (and especially numbers)

#include <iomanip>

#include "Quirks/VFNStreamFormat.h"
#include "Quirks/VFNDebug.h"

////////////////////////////////////////////////////////////////////////
//
//    FormatInt
//
////////////////////////////////////////////////////////////////////////

FormatInt::FormatInt(const long num,const int width): mValue(num),mWidth(width) {}
FormatInt::~FormatInt(){}

ostream& operator<< (ostream& os, const FormatInt& fInt)
{
   return os << setiosflags(ios::right) << setw(fInt.mWidth) <<
         fInt.mValue << " ";
}

////////////////////////////////////////////////////////////////////////
//
//    FormatFloat
//
////////////////////////////////////////////////////////////////////////

FormatFloat::FormatFloat(const REAL num,const int width,const int precision)
         :mValue(num),mWidth(width),mPrecision(precision) 
         {
         }
FormatFloat::~FormatFloat(){}

ostream& operator<< (ostream& os, const FormatFloat &fFloat)
{
   //:TODO: check the formatting ->DOESN'T WORK
   //#ifdef __MWERKS__
   //return os << fFloat.mValue << " ";
   //#endif
   return os << setiosflags(ios::right) <<
                setiosflags(ios::fixed)<<
                setprecision(fFloat.mPrecision) <<
                setiosflags(ios::showpoint)<<
                setw(fFloat.mWidth) <<
                fFloat.mValue << " ";
}

////////////////////////////////////////////////////////////////////////
//
//    FormatString
//
////////////////////////////////////////////////////////////////////////
FormatString::FormatString(const string &str,const unsigned int width): 
            /*mString(str,1,width),*/mWidth(width) 
{
   mString=str;
   //:KLUDGE:
   if(mString.size() > mWidth) mString.resize(mWidth);
   else
      for(unsigned int i=mString.size();i<width;i++) mString += " ";
   //cout << mString.size()<<endl;
   //cout << str << " " << str.length()<< " " << str.size()<<endl;
}
FormatString::~FormatString(){}
int FormatString::length() const {return mString.length();}

ostream& operator<< (ostream& os, const FormatString& fStr)
{
   return os << fStr.mString ;
}

////////////////////////////////////////////////////////////////////////
//
//    FormatVertVector
//
////////////////////////////////////////////////////////////////////////

template<class T> FormatVertVector<T>::FormatVertVector( const CrystVector<T> &fVect,
                                    const int width,
                                    const int precision):
mNbVectors(1),mWidth(width),mPrecision(precision)
{
   mpVectors=new const CrystVector<T>*[mNbVectors];
   mpVectors[0]=&fVect;
}

template<class T> FormatVertVector<T>::FormatVertVector( const CrystVector<T> &fVect1,
                                                         const CrystVector<T> &fVect2,
                                    const int width,
                                    const int precision):
mNbVectors(2),mWidth(width),mPrecision(precision)
{
   mpVectors=new const CrystVector<T>*[mNbVectors];
   mpVectors[0]=&fVect1;
   mpVectors[1]=&fVect2;
}

template<class T> FormatVertVector<T>::FormatVertVector( const CrystVector<T> &fVect1,
                                                         const CrystVector<T> &fVect2,
                                                         const CrystVector<T> &fVect3,
                                    const int width,
                                    const int precision):
mNbVectors(3),mWidth(width),mPrecision(precision)
{
   mpVectors=new const CrystVector<T>*[mNbVectors];
   mpVectors[0]=&fVect1;
   mpVectors[1]=&fVect2;
   mpVectors[2]=&fVect3;
}

template<class T> FormatVertVector<T>::FormatVertVector( const CrystVector<T> &fVect1,
                                                         const CrystVector<T> &fVect2,
                                                         const CrystVector<T> &fVect3,
                                                         const CrystVector<T> &fVect4,
                                    const int width,
                                    const int precision):
mNbVectors(4),mWidth(width),mPrecision(precision)
{
   mpVectors=new const CrystVector<T>*[mNbVectors];
   mpVectors[0]=&fVect1;
   mpVectors[1]=&fVect2;
   mpVectors[2]=&fVect3;
   mpVectors[3]=&fVect4;
}
template<class T> FormatVertVector<T>::FormatVertVector( const CrystVector<T> &fVect1,
                                                         const CrystVector<T> &fVect2,
                                                         const CrystVector<T> &fVect3,
                                                         const CrystVector<T> &fVect4,
                                                         const CrystVector<T> &fVect5,
                                    const int width,
                                    const int precision):
mNbVectors(5),mWidth(width),mPrecision(precision)
{
   mpVectors=new const CrystVector<T>*[mNbVectors];
   mpVectors[0]=&fVect1;
   mpVectors[1]=&fVect2;
   mpVectors[2]=&fVect3;
   mpVectors[3]=&fVect4;
   mpVectors[4]=&fVect5;
}
template<class T> FormatVertVector<T>::FormatVertVector( const CrystVector<T> &fVect1,
                                                         const CrystVector<T> &fVect2,
                                                         const CrystVector<T> &fVect3,
                                                         const CrystVector<T> &fVect4,
                                                         const CrystVector<T> &fVect5,
                                                         const CrystVector<T> &fVect6,
                                    const int width,
                                    const int precision):
mNbVectors(6),mWidth(width),mPrecision(precision)
{
   mpVectors=new const CrystVector<T>*[mNbVectors];
   mpVectors[0]=&fVect1;
   mpVectors[1]=&fVect2;
   mpVectors[2]=&fVect3;
   mpVectors[3]=&fVect4;
   mpVectors[4]=&fVect5;
   mpVectors[5]=&fVect6;
}

template<class T> FormatVertVector<T>::FormatVertVector( const CrystVector<T> *pVect,
                                   const int nbVect,
                                   const int width,
                                   const int precision):
mNbVectors(nbVect),mWidth(width),mPrecision(precision)
{
   mpVectors=new const CrystVector<T>*[mNbVectors];
   for(int i=0;i<mNbVectors;i++) mpVectors[i]=&(pVect[i]);
}

template<class T> FormatVertVector<T>::FormatVertVector( const CrystVector<T> &fVect1,
                                   const CrystVector<T> *pVect,
                                   const int nbVect,
                                   const int width,
                                   const int precision):
mNbVectors(nbVect+1),mWidth(width),mPrecision(precision)
{
   mpVectors=new const CrystVector<T>*[mNbVectors];
   mpVectors[0]=&fVect1;
   for(int i=1;i<mNbVectors;i++) mpVectors[i]=&(pVect[i-1]);
}

template<class T> FormatVertVector<T>::~FormatVertVector()
{
   delete[] mpVectors;
}

template<class T> ostream& operator<< (ostream &os, const FormatVertVector<T> &fVect)
{
   VFN_DEBUG_MESSAGE("ostream& operator<<(os,FormatVertVector<T>)",3)
   long i;
   int j,width,precision;
   width=fVect.mWidth;
   precision=fVect.mPrecision;
   for(i=0;i<(fVect.mpVectors[0])->numElements();i++)
   {
      for(j=0;j<fVect.mNbVectors;j++)
      {
         os << setiosflags(ios::right) <<
               setiosflags(ios::fixed)<<
               setprecision(precision) <<
               setiosflags(ios::showpoint)<<
               setw(width) <<
               (*fVect.mpVectors[j])(i) << " ";
      }
      os << endl;
   }
   VFN_DEBUG_MESSAGE("ostream& operator<<(os,FormatVertVector<T>):End",2)
   return os;
}

////////////////////////////////////////////////////////////////////////
//
//    FormatHorizVector
//
////////////////////////////////////////////////////////////////////////

template<class T> FormatHorizVector<T>::FormatHorizVector( const CrystVector<T> &fVect,
                                    const int width,
                                    const int precision):
mWidth(width),mPrecision(precision)
{
   mpVectors=&fVect;
}
template<class T> FormatHorizVector<T>::~FormatHorizVector()
{}

template<class T> ostream& operator<< (ostream &os, const FormatHorizVector<T> &fVect)
{
   long i;
   int width,precision;
   width=fVect.mWidth;
   precision=fVect.mPrecision;
   for(i=0;i<(*(fVect.mpVectors)).numElements();i++)
   {
      os << setiosflags(ios::right) <<
            setiosflags(ios::fixed)<<
            setprecision(precision) <<
            setiosflags(ios::showpoint)<<
            setw(width) <<
            (*(fVect.mpVectors))(i) << " ";
   }
   os<<endl;
   return os;
}

////////////////////////////////////////////////////////////////////////
//
//    FormatVertVectorHKLFloats
//
////////////////////////////////////////////////////////////////////////
template<class T> FormatVertVectorHKLFloats<T>::FormatVertVectorHKLFloats( 
                                                      const CrystVector<T> &h,
                                                      const CrystVector<T> &k,
                                                      const CrystVector<T> &l,
                                                      const int width,
                                                      const int precision):
mNbVectors(3),mWidth(width),mPrecision(precision)
{
   mpVectors[0]=&h;
   mpVectors[1]=&k;
   mpVectors[2]=&l;
}
template<class T> FormatVertVectorHKLFloats<T>::FormatVertVectorHKLFloats( 
                                                      const CrystVector<T> &h,
                                                      const CrystVector<T> &k,
                                                      const CrystVector<T> &l,
                                                      const CrystVector<T> &m,
                                                      const int width,
                                                      const int precision):
mNbVectors(4),mWidth(width),mPrecision(precision)
{
   mpVectors[0]=&h;
   mpVectors[1]=&k;
   mpVectors[2]=&l;
   mpVectors[3]=&m;
}
template<class T> FormatVertVectorHKLFloats<T>::FormatVertVectorHKLFloats(
                                                      const CrystVector<T> &h,
                                                      const CrystVector<T> &k,
                                                      const CrystVector<T> &l,
                                                      const CrystVector<T> &m,
                                                      const CrystVector<T> &n,
                                                      const int width,
                                                      const int precision):
mNbVectors(5),mWidth(width),mPrecision(precision)
{
   mpVectors[0]=&h;
   mpVectors[1]=&k;
   mpVectors[2]=&l;
   mpVectors[3]=&m;
   mpVectors[4]=&n;
}
template<class T> FormatVertVectorHKLFloats<T>::FormatVertVectorHKLFloats(
                                                      const CrystVector<T> &h,
                                                      const CrystVector<T> &k,
                                                      const CrystVector<T> &l,
                                                      const CrystVector<T> &m,
                                                      const CrystVector<T> &n,
                                                      const CrystVector<T> &o,
                                                      const int width,
                                                      const int precision):
mNbVectors(6),mWidth(width),mPrecision(precision)
{
   mpVectors[0]=&h;
   mpVectors[1]=&k;
   mpVectors[2]=&l;
   mpVectors[3]=&m;
   mpVectors[4]=&n;
   mpVectors[5]=&o;
}
template<class T> FormatVertVectorHKLFloats<T>::FormatVertVectorHKLFloats(
                                                      const CrystVector<T> &h,
                                                      const CrystVector<T> &k,
                                                      const CrystVector<T> &l,
                                                      const CrystVector<T> &m,
                                                      const CrystVector<T> &n,
                                                      const CrystVector<T> &o,
                                                      const CrystVector<T> &p,
                                                      const int width,
                                                      const int precision):
mNbVectors(7),mWidth(width),mPrecision(precision)
{
   mpVectors[0]=&h;
   mpVectors[1]=&k;
   mpVectors[2]=&l;
   mpVectors[3]=&m;
   mpVectors[4]=&n;
   mpVectors[5]=&o;
   mpVectors[6]=&p;
}
template<class T> FormatVertVectorHKLFloats<T>::FormatVertVectorHKLFloats(
                                                      const CrystVector<T> &h,
                                                      const CrystVector<T> &k,
                                                      const CrystVector<T> &l,
                                                      const CrystVector<T> &m,
                                                      const CrystVector<T> &n,
                                                      const CrystVector<T> &o,
                                                      const CrystVector<T> &p,
                                                      const CrystVector<T> &q,
                                                      const int width,
                                                      const int precision):
mNbVectors(8),mWidth(width),mPrecision(precision)
{
   mpVectors[0]=&h;
   mpVectors[1]=&k;
   mpVectors[2]=&l;
   mpVectors[3]=&m;
   mpVectors[4]=&n;
   mpVectors[5]=&o;
   mpVectors[6]=&p;
   mpVectors[7]=&q;
}
template<class T> FormatVertVectorHKLFloats<T>::FormatVertVectorHKLFloats(
                                                      const CrystVector<T> &h,
                                                      const CrystVector<T> &k,
                                                      const CrystVector<T> &l,
                                                      const CrystVector<T> &m,
                                                      const CrystVector<T> &n,
                                                      const CrystVector<T> &o,
                                                      const CrystVector<T> &p,
                                                      const CrystVector<T> &q,
                                                      const CrystVector<T> &r,
                                                      const int width,
                                                      const int precision):
mNbVectors(9),mWidth(width),mPrecision(precision)
{
   mpVectors[0]=&h;
   mpVectors[1]=&k;
   mpVectors[2]=&l;
   mpVectors[3]=&m;
   mpVectors[4]=&n;
   mpVectors[5]=&o;
   mpVectors[6]=&p;
   mpVectors[7]=&q;
   mpVectors[8]=&r;
}
template<class T> FormatVertVectorHKLFloats<T>::FormatVertVectorHKLFloats(
                                                      const CrystVector<T> &h,
                                                      const CrystVector<T> &k,
                                                      const CrystVector<T> &l,
                                                      const CrystVector<T> &m,
                                                      const CrystVector<T> &n,
                                                      const CrystVector<T> &o,
                                                      const CrystVector<T> &p,
                                                      const CrystVector<T> &q,
                                                      const CrystVector<T> &r,
                                                      const CrystVector<T> &s,
                                                      const int width,
                                                      const int precision):
mNbVectors(10),mWidth(width),mPrecision(precision)
{
   mpVectors[0]=&h;
   mpVectors[1]=&k;
   mpVectors[2]=&l;
   mpVectors[3]=&m;
   mpVectors[4]=&n;
   mpVectors[5]=&o;
   mpVectors[6]=&p;
   mpVectors[7]=&q;
   mpVectors[8]=&r;
   mpVectors[9]=&s;
}
template<class T> FormatVertVectorHKLFloats<T>::FormatVertVectorHKLFloats(
                                                      const CrystVector<T> &h,
                                                      const CrystVector<T> &k,
                                                      const CrystVector<T> &l,
                                                      const CrystVector<T> &m,
                                                      const CrystVector<T> &n,
                                                      const CrystVector<T> &o,
                                                      const CrystVector<T> &p,
                                                      const CrystVector<T> &q,
                                                      const CrystVector<T> &r,
                                                      const CrystVector<T> &s,
                                                      const CrystVector<T> &t,
                                                      const int width,
                                                      const int precision):
mNbVectors(11),mWidth(width),mPrecision(precision)
{
   mpVectors[0]=&h;
   mpVectors[1]=&k;
   mpVectors[2]=&l;
   mpVectors[3]=&m;
   mpVectors[4]=&n;
   mpVectors[5]=&o;
   mpVectors[6]=&p;
   mpVectors[7]=&q;
   mpVectors[8]=&r;
   mpVectors[9]=&s;
   mpVectors[10]=&t;
}
template<class T> FormatVertVectorHKLFloats<T>::FormatVertVectorHKLFloats(
                                                      const CrystVector<T> &h,
                                                      const CrystVector<T> &k,
                                                      const CrystVector<T> &l,
                                                      const CrystVector<T> &m,
                                                      const CrystVector<T> &n,
                                                      const CrystVector<T> &o,
                                                      const CrystVector<T> &p,
                                                      const CrystVector<T> &q,
                                                      const CrystVector<T> &r,
                                                      const CrystVector<T> &s,
                                                      const CrystVector<T> &t,
                                                      const CrystVector<T> &u,
                                                      const int width,
                                                      const int precision):
mNbVectors(12),mWidth(width),mPrecision(precision)
{
   mpVectors[0]=&h;
   mpVectors[1]=&k;
   mpVectors[2]=&l;
   mpVectors[3]=&m;
   mpVectors[4]=&n;
   mpVectors[5]=&o;
   mpVectors[6]=&p;
   mpVectors[7]=&q;
   mpVectors[8]=&r;
   mpVectors[9]=&s;
   mpVectors[10]=&t;
   mpVectors[11]=&u;
}
template<class T> FormatVertVectorHKLFloats<T>::~FormatVertVectorHKLFloats()
{}

template<class T> ostream& operator<< (ostream& os, const FormatVertVectorHKLFloats<T> &fVect)
{
   long i;
   int j,width,precision;
   width=fVect.mWidth;
   precision=fVect.mPrecision;
   for(i=0 ; i<(fVect.mpVectors[0])->numElements() ; i++)
   {
      for(j=0;j<3;j++)
      {
         os << setiosflags(ios::right) << setw(width-precision) <<
            (int) ((*fVect.mpVectors[j])(i)) << " ";
      }
      for(j=3;j<fVect.mNbVectors;j++)
      {
         os << setiosflags(ios::right) <<
               setiosflags(ios::fixed)<<
               setprecision(precision) <<
               setiosflags(ios::showpoint)<<
               setw(width) <<
               (*fVect.mpVectors[j])(i) << " ";
      }
      os << endl;
   }
   return os;
}


//Explicit instantiations
template class FormatVertVector<REAL>;
template ostream& operator<< (ostream&,const FormatVertVector<REAL>&);
template class FormatHorizVector<REAL>;
template ostream& operator<< (ostream&,const FormatHorizVector<REAL>&);
template class FormatVertVectorHKLFloats<REAL>;
template ostream& operator<< (ostream&,const FormatVertVectorHKLFloats<REAL>&);
template class FormatVertVector<long>;
template ostream& operator<< (ostream&,const FormatVertVector<long>&);
template class FormatHorizVector<long>;
template ostream& operator<< (ostream&,const FormatHorizVector<long>&);
template class FormatVertVectorHKLFloats<long>;
template ostream& operator<< (ostream&,const FormatVertVectorHKLFloats<long>&);
