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

#include "ObjCryst/Quirks/VFNStreamFormat.h"
#include "ObjCryst/Quirks/VFNDebug.h"

using namespace std;

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
   std::istream::fmtflags old_flags=os.flags();
   os.setf( std::istream::fixed | std::istream::right | std::istream::showpoint,std::istream::floatfield );
   std::streamsize old_prec = os.precision(fFloat.mPrecision);
   std::streamsize old_width = os.width(fFloat.mWidth);
   os << fFloat.mValue;
   os.flags( old_flags );
   os.precision( old_prec );
   os.width(old_width);
   return os;
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
                                    const int precision, const int nb):
mWidth(width),mPrecision(precision),mNb(nb)
{
   mvpVectors.push_back(&fVect);
}

template<class T> FormatVertVector<T>::FormatVertVector( const CrystVector<T> &fVect1,
                                                         const CrystVector<T> &fVect2,
                                    const int width,
                                    const int precision, const int nb):
mWidth(width),mPrecision(precision),mNb(nb)
{
   mvpVectors.push_back(&fVect1);
   mvpVectors.push_back(&fVect2);
}

template<class T> FormatVertVector<T>::FormatVertVector( const CrystVector<T> &fVect1,
                                                         const CrystVector<T> &fVect2,
                                                         const CrystVector<T> &fVect3,
                                    const int width,
                                    const int precision, const int nb):
mWidth(width),mPrecision(precision),mNb(nb)
{
   mvpVectors.push_back(&fVect1);
   mvpVectors.push_back(&fVect2);
   mvpVectors.push_back(&fVect3);
}

template<class T> FormatVertVector<T>::FormatVertVector( const CrystVector<T> &fVect1,
                                                         const CrystVector<T> &fVect2,
                                                         const CrystVector<T> &fVect3,
                                                         const CrystVector<T> &fVect4,
                                    const int width,
                                    const int precision, const int nb):
mWidth(width),mPrecision(precision),mNb(nb)
{
   mvpVectors.push_back(&fVect1);
   mvpVectors.push_back(&fVect2);
   mvpVectors.push_back(&fVect3);
   mvpVectors.push_back(&fVect4);
}
template<class T> FormatVertVector<T>::FormatVertVector( const CrystVector<T> &fVect1,
                                                         const CrystVector<T> &fVect2,
                                                         const CrystVector<T> &fVect3,
                                                         const CrystVector<T> &fVect4,
                                                         const CrystVector<T> &fVect5,
                                    const int width,
                                    const int precision, const int nb):
mWidth(width),mPrecision(precision),mNb(nb)
{
   mvpVectors.push_back(&fVect1);
   mvpVectors.push_back(&fVect2);
   mvpVectors.push_back(&fVect3);
   mvpVectors.push_back(&fVect4);
   mvpVectors.push_back(&fVect5);
}
template<class T> FormatVertVector<T>::FormatVertVector( const CrystVector<T> &fVect1,
                                                         const CrystVector<T> &fVect2,
                                                         const CrystVector<T> &fVect3,
                                                         const CrystVector<T> &fVect4,
                                                         const CrystVector<T> &fVect5,
                                                         const CrystVector<T> &fVect6,
                                    const int width,
                                    const int precision, const int nb):
mWidth(width),mPrecision(precision),mNb(nb)
{
   mvpVectors.push_back(&fVect1);
   mvpVectors.push_back(&fVect2);
   mvpVectors.push_back(&fVect3);
   mvpVectors.push_back(&fVect4);
   mvpVectors.push_back(&fVect5);
   mvpVectors.push_back(&fVect6);
}

template<class T> FormatVertVector<T>::FormatVertVector( const CrystVector<T> *pVect,
                                   const int nbVect,
                                   const int width,
                                   const int precision, const int nb):
mWidth(width),mPrecision(precision),mNb(nb)
{
   for(int i=0;i<nbVect;i++) mvpVectors.push_back(&(pVect[i]));
}

template<class T> FormatVertVector<T>::FormatVertVector( const CrystVector<T> &fVect1,
                                   const CrystVector<T> *pVect,
                                   const int nbVect,
                                   const int width,
                                   const int precision, const int nb):
mWidth(width),mPrecision(precision),mNb(nb)
{
   mvpVectors.push_back(&fVect1);
   for(int i=1;i<nbVect;i++) mvpVectors.push_back(&(pVect[i-1]));
}

template<class T> FormatVertVector<T>::FormatVertVector(
                              vector<const CrystVector<T> *>& v,
                              const int width,
                              const int precision, const int nb):
mWidth(width),mPrecision(precision),mNb(nb)
{
   mvpVectors=v;
}

template<class T> FormatVertVector<T>::~FormatVertVector()
{
}

template<class T> ostream& operator<< (ostream &os, const FormatVertVector<T> &fVect)
{
   VFN_DEBUG_MESSAGE("ostream& operator<<(os,FormatVertVector<T>)",3)
   long i;
   int j;
   std::istream::fmtflags old_flags=os.flags();
   os.setf( std::istream::fixed | std::istream::right | std::istream::showpoint);
   std::streamsize old_prec = os.precision(fVect.mPrecision);
   std::streamsize old_width = os.width();
   long nb=fVect.mNb;
   if(nb==0)nb=(fVect.mvpVectors[0])->numElements();
   for(i=0;i<nb;i++)
   {
      for(j=0;j<fVect.mvpVectors.size();j++)
      {
         os <<setw(fVect.mWidth)<< (*fVect.mvpVectors[j])(i) << " ";
      }
      os << endl;
   }
   os.flags( old_flags );
   os.precision( old_prec );
   os<<setw(old_width);
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
   std::istream::fmtflags old_flags=os.flags();
   os.setf( std::istream::fixed | std::istream::right | std::istream::showpoint);
   std::streamsize old_prec = os.precision(fVect.mPrecision);
   std::streamsize old_width = os.width();
   for(i=0;i<(*(fVect.mpVectors)).numElements();i++)
   {
      os <<setw(fVect.mWidth)<< (*(fVect.mpVectors))(i) << " ";
   }
   os<<endl;
   os.flags( old_flags );
   os.precision( old_prec );
   os<<setw(old_width);
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
                                                      const int precision,
                                                      const int nb):
mWidth(width),mPrecision(precision),mNb(nb)
{
   mvpVectors.push_back(&h);
   mvpVectors.push_back(&k);
   mvpVectors.push_back(&l);
}
template<class T> FormatVertVectorHKLFloats<T>::FormatVertVectorHKLFloats(
                                                      const CrystVector<T> &h,
                                                      const CrystVector<T> &k,
                                                      const CrystVector<T> &l,
                                                      const CrystVector<T> &m,
                                                      const int width,
                                                      const int precision,
                                                      const int nb):
mWidth(width),mPrecision(precision),mNb(nb)
{
   mvpVectors.push_back(&h);
   mvpVectors.push_back(&k);
   mvpVectors.push_back(&l);
   mvpVectors.push_back(&m);
}
template<class T> FormatVertVectorHKLFloats<T>::FormatVertVectorHKLFloats(
                                                      const CrystVector<T> &h,
                                                      const CrystVector<T> &k,
                                                      const CrystVector<T> &l,
                                                      const CrystVector<T> &m,
                                                      const CrystVector<T> &n,
                                                      const int width,
                                                      const int precision,
                                                      const int nb):
mWidth(width),mPrecision(precision),mNb(nb)
{
   mvpVectors.push_back(&h);
   mvpVectors.push_back(&k);
   mvpVectors.push_back(&l);
   mvpVectors.push_back(&m);
   mvpVectors.push_back(&n);
}
template<class T> FormatVertVectorHKLFloats<T>::FormatVertVectorHKLFloats(
                                                      const CrystVector<T> &h,
                                                      const CrystVector<T> &k,
                                                      const CrystVector<T> &l,
                                                      const CrystVector<T> &m,
                                                      const CrystVector<T> &n,
                                                      const CrystVector<T> &o,
                                                      const int width,
                                                      const int precision,
                                                      const int nb):
mWidth(width),mPrecision(precision),mNb(nb)
{
   mvpVectors.push_back(&h);
   mvpVectors.push_back(&k);
   mvpVectors.push_back(&l);
   mvpVectors.push_back(&m);
   mvpVectors.push_back(&n);
   mvpVectors.push_back(&o);
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
                                                      const int precision,
                                                      const int nb):
mWidth(width),mPrecision(precision),mNb(nb)
{
   mvpVectors.push_back(&h);
   mvpVectors.push_back(&k);
   mvpVectors.push_back(&l);
   mvpVectors.push_back(&m);
   mvpVectors.push_back(&n);
   mvpVectors.push_back(&o);
   mvpVectors.push_back(&p);
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
                                                      const int precision,
                                                      const int nb):
mWidth(width),mPrecision(precision),mNb(nb)
{
   mvpVectors.push_back(&h);
   mvpVectors.push_back(&k);
   mvpVectors.push_back(&l);
   mvpVectors.push_back(&m);
   mvpVectors.push_back(&n);
   mvpVectors.push_back(&o);
   mvpVectors.push_back(&p);
   mvpVectors.push_back(&q);
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
                                                      const int precision,
                                                      const int nb):
mWidth(width),mPrecision(precision),mNb(nb)
{
   mvpVectors.push_back(&h);
   mvpVectors.push_back(&k);
   mvpVectors.push_back(&l);
   mvpVectors.push_back(&m);
   mvpVectors.push_back(&n);
   mvpVectors.push_back(&o);
   mvpVectors.push_back(&p);
   mvpVectors.push_back(&q);
   mvpVectors.push_back(&r);
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
                                                      const int precision,
                                                      const int nb):
mWidth(width),mPrecision(precision),mNb(nb)
{
   mvpVectors.push_back(&h);
   mvpVectors.push_back(&k);
   mvpVectors.push_back(&l);
   mvpVectors.push_back(&m);
   mvpVectors.push_back(&n);
   mvpVectors.push_back(&o);
   mvpVectors.push_back(&p);
   mvpVectors.push_back(&q);
   mvpVectors.push_back(&r);
   mvpVectors.push_back(&s);
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
                                                      const int precision,
                                                      const int nb):
mWidth(width),mPrecision(precision),mNb(nb)
{
   mvpVectors.push_back(&h);
   mvpVectors.push_back(&k);
   mvpVectors.push_back(&l);
   mvpVectors.push_back(&m);
   mvpVectors.push_back(&n);
   mvpVectors.push_back(&o);
   mvpVectors.push_back(&p);
   mvpVectors.push_back(&q);
   mvpVectors.push_back(&r);
   mvpVectors.push_back(&s);
   mvpVectors.push_back(&t);
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
                                                      const int precision,
                                                      const int nb):
mWidth(width),mPrecision(precision),mNb(nb)
{
   mvpVectors.push_back(&h);
   mvpVectors.push_back(&k);
   mvpVectors.push_back(&l);
   mvpVectors.push_back(&m);
   mvpVectors.push_back(&n);
   mvpVectors.push_back(&o);
   mvpVectors.push_back(&p);
   mvpVectors.push_back(&q);
   mvpVectors.push_back(&r);
   mvpVectors.push_back(&s);
   mvpVectors.push_back(&t);
   mvpVectors.push_back(&u);
}
template<class T> FormatVertVectorHKLFloats<T>::FormatVertVectorHKLFloats(
                              vector<const CrystVector<T> *>& v,
                              const int width,
                              const int precision, const int nb):
mWidth(width),mPrecision(precision),mNb(nb)
{
   mvpVectors=v;
}

template<class T> FormatVertVectorHKLFloats<T>::~FormatVertVectorHKLFloats()
{}

template<class T> ostream& operator<< (ostream& os, const FormatVertVectorHKLFloats<T> &fVect)
{
   long i;
   unsigned int j;
   std::istream::fmtflags old_flags=os.flags();
   os.setf( std::istream::fixed | std::istream::right | std::istream::showpoint);
   std::streamsize old_prec = os.precision(fVect.mPrecision);
   std::streamsize old_width = os.width();
   long nb=fVect.mNb;
   if(nb==0) nb=(fVect.mvpVectors[0])->numElements();
   for(i=0;i<nb;i++)
   {
      for(j=0;j<3;j++)
      {
         os <<setw(fVect.mWidth-fVect.mPrecision)<< (int) ((*fVect.mvpVectors[j])(i)) << " ";
      }
      for(j=3;j<fVect.mvpVectors.size();j++)
      {
         os <<setw(fVect.mWidth)<<(*fVect.mvpVectors[j])(i) << " ";
      }
      os << endl;
   }
   os.flags( old_flags );
   os.precision( old_prec );
   os<<setw(old_width);
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
