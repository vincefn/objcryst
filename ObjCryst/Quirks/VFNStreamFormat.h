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
//strings, numbers, some arrays...
#ifndef _VFN_STREAM_FORMAT_H_
#define _VFN_STREAM_FORMAT_H_

#include <string>
#include <iostream>
#include <vector>

using namespace std;

#include "CrystVector/CrystVector.h"

/** output a number as a formatted integer: 
*
* \code os << FormatInt(mynumber,5);\endcode
*
*
*/
class FormatInt
{
   public:
      FormatInt(const long num,const int width=5);
      ~FormatInt();
   private:
      const long mValue;
      const int mWidth;
   
   friend ostream& operator<< (ostream& os, const FormatInt &fInt);
};

ostream& operator<< (ostream& os, const FormatInt& fInt);


/** output a number as a formatted float:
*
*\code os << FormatFloat(mynumber,10,4); \endcode
*
*
*/
class FormatFloat
{
   public:
      FormatFloat(const REAL num,const int width=10,const int precision=4);
      ~FormatFloat();
      
   private:
      const REAL mValue;
      const int mWidth;
      const int mPrecision;
  friend ostream& operator<< (ostream&, const FormatFloat&);
};

ostream& operator<< (ostream& os, const FormatFloat &fFloat);

/** output a string with a fixed length (adding necessary space or removing
* excess characters) : 
*
* \code os << FormatString(myString,15);\endcode
*/
class FormatString
{
   public:
      FormatString(const string &str,const unsigned int width=5);
      ~FormatString();
      int length() const;
   private:
      string mString;
      const unsigned int mWidth;
   
   friend ostream& operator<< (ostream&, const FormatString&);
};

ostream& operator<< (ostream& os, const FormatString& fStr);

/** output one or several vectors as (a) column(s): 
*
* \code 
*  os << FormatVertVector<REAL>(vect,8,3);
*  os << FormatVertVector<REAL>(vect1,vect2,vetc3,12,6);
*  // For 7 vectors with width 12 and precision 4,
*  // pVect being a pointer to an array of 7 vectors:
*  os << FormatVertVector<REAL>(pVect,7,12,4);
* \endcode
*/
template<class T> class FormatVertVector
{
   public:
      FormatVertVector( const CrystVector<T> &fVect,
                        const int width=10,
                        const int precision=4);
      FormatVertVector( const CrystVector<T> &fVect1,
                        const CrystVector<T> &fVect2,
                        const int width=10,
                        const int precision=4);
      FormatVertVector( const CrystVector<T> &fVect1,
                        const CrystVector<T> &fVect2,
                        const CrystVector<T> &fVect3,
                        const int width=10,
                        const int precision=4);
      FormatVertVector( const CrystVector<T> &fVect1,
                        const CrystVector<T> &fVect2,
                        const CrystVector<T> &fVect3,
                        const CrystVector<T> &fVect4,
                        const int width=10,
                        const int precision=4);
      FormatVertVector( const CrystVector<T> &fVect1,
                        const CrystVector<T> &fVect2,
                        const CrystVector<T> &fVect3,
                        const CrystVector<T> &fVect4,
                        const CrystVector<T> &fVect5,
                        const int width=10,
                        const int precision=4);
      FormatVertVector( const CrystVector<T> &fVect1,
                        const CrystVector<T> &fVect2,
                        const CrystVector<T> &fVect3,
                        const CrystVector<T> &fVect4,
                        const CrystVector<T> &fVect5,
                        const CrystVector<T> &fVect6,
                        const int width=10,
                        const int precision=4);
      FormatVertVector( const CrystVector<T> *pVect,
                        const int nbVect,
                        const int width=10,
                        const int precision=4);
      FormatVertVector( const CrystVector<T> &fVect1,
                        const CrystVector<T> *pVect,
                        const int nbVect,
                        const int width=10,
                        const int precision=4);
      ~FormatVertVector();
      //int length() const;
   //private:
      const CrystVector<T> **mpVectors;
      const int mNbVectors;
      const int mWidth;
      const int mPrecision;
   
   friend ostream& operator<< <T>(ostream&, const FormatVertVector<T>&);
};

template<class T> ostream& operator<< (ostream &os, const FormatVertVector<T> &fVect);

/** Format vector as horiz array:
*
* \code os << FormatHorizVector<REAL>(vect,8,3);\endcode
*/
template<class T> class FormatHorizVector
{
   public:
      FormatHorizVector(const CrystVector<T> &fVect,
                        const int width=10,
                        const int precision=4);
      ~FormatHorizVector();
      //int length() const;
   //private:
      const CrystVector<T> *mpVectors;
      const int mWidth;
      const int mPrecision;
   
   friend ostream& operator<< <T>(ostream&, const FormatHorizVector<T>&);
};

template<class T> ostream& operator<< (ostream &os, const FormatHorizVector<T> &fVect);

/** Output vectors as column arrays, with the first 3 columns printed as integers.
*
* \code cout << FormatVertVectorHKLFloats<REAL>(vH,vK,vL,vIobs,vIcalc,vSigma,12,4);\endcode
*/
template<class T>class FormatVertVectorHKLFloats
{
   public:
      FormatVertVectorHKLFloats( const CrystVector<T> &h,
                                 const CrystVector<T> &k,
                                 const CrystVector<T> &l,
                                 const int width=10,
                                 const int precision=4);
      FormatVertVectorHKLFloats( const CrystVector<T> &h,
                                 const CrystVector<T> &k,
                                 const CrystVector<T> &l,
                                 const CrystVector<T> &m,
                                 const int width=10,
                                 const int precision=4);
      FormatVertVectorHKLFloats( const CrystVector<T> &h,
                                 const CrystVector<T> &k,
                                 const CrystVector<T> &l,
                                 const CrystVector<T> &m,
                                 const CrystVector<T> &n,
                                 const int width=10,
                                 const int precision=4);
      FormatVertVectorHKLFloats( const CrystVector<T> &h,
                                 const CrystVector<T> &k,
                                 const CrystVector<T> &l,
                                 const CrystVector<T> &m,
                                 const CrystVector<T> &n,
                                 const CrystVector<T> &o,
                                 const int width=10,
                                 const int precision=4);
      FormatVertVectorHKLFloats( const CrystVector<T> &h,
                                 const CrystVector<T> &k,
                                 const CrystVector<T> &l,
                                 const CrystVector<T> &m,
                                 const CrystVector<T> &n,
                                 const CrystVector<T> &o,
                                 const CrystVector<T> &p,
                                 const int width=10,
                                 const int precision=4);
      FormatVertVectorHKLFloats( const CrystVector<T> &h,
                                 const CrystVector<T> &k,
                                 const CrystVector<T> &l,
                                 const CrystVector<T> &m,
                                 const CrystVector<T> &n,
                                 const CrystVector<T> &o,
                                 const CrystVector<T> &p,
                                 const CrystVector<T> &q,
                                 const int width=10,
                                 const int precision=4);
      FormatVertVectorHKLFloats( const CrystVector<T> &h,
                                 const CrystVector<T> &k,
                                 const CrystVector<T> &l,
                                 const CrystVector<T> &m,
                                 const CrystVector<T> &n,
                                 const CrystVector<T> &o,
                                 const CrystVector<T> &p,
                                 const CrystVector<T> &q,
                                 const CrystVector<T> &r,
                                 const int width=10,
                                 const int precision=4);
      FormatVertVectorHKLFloats( const CrystVector<T> &h,
                                 const CrystVector<T> &k,
                                 const CrystVector<T> &l,
                                 const CrystVector<T> &m,
                                 const CrystVector<T> &n,
                                 const CrystVector<T> &o,
                                 const CrystVector<T> &p,
                                 const CrystVector<T> &q,
                                 const CrystVector<T> &r,
                                 const CrystVector<T> &s,
                                 const int width=10,
                                 const int precision=4);
      FormatVertVectorHKLFloats( const CrystVector<T> &h,
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
                                 const int width=10,
                                 const int precision=4);
      FormatVertVectorHKLFloats( const CrystVector<T> &h,
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
                                 const int width=10,
                                 const int precision=4);
      FormatVertVectorHKLFloats( vector<const CrystVector<T> *>& v,
                                 const int width=10,
                                 const int precision=4);
      ~FormatVertVectorHKLFloats();
      //int length() const;
   //private:
      vector<const CrystVector<T> *>mvpVectors;
      const int mWidth;
      const int mPrecision;
   
   friend ostream& operator<< <T>(ostream&, const FormatVertVectorHKLFloats<T>&);
};

template<class T> ostream& operator<< (ostream& os, const FormatVertVectorHKLFloats<T> &fStr);


#endif

