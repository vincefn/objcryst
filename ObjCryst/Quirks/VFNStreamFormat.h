// This file includes some "effectors" used to format
//strings, numbers, some arrays...
#ifndef _VFN_STREAM_FORMAT_H_
#define _VFN_STREAM_FORMAT_H_

#include <string>
#include <iostream>

using namespace std;

#include "CrystVector/CrystVector.h"


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


class FormatFloat
{
   public:
      FormatFloat(const double num,const int width=10,const int precision=4);
      ~FormatFloat();
      
   private:
      const double mValue;
      const int mWidth;
      const int mPrecision;
  friend ostream& operator<< (ostream&, const FormatFloat&);
};

ostream& operator<< (ostream& os, const FormatFloat &fFloat);

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

/**
*  Format vectors as column arrays
*
* \example: cout << FormatVertVector<double>(vect,8,3);
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

/**
*  Format vector as horiz array
*
* \example: cout << FormatHorizVector<double>(vect,8,3);
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

/**
*  Format vectors as column arrays
*Here the first 3 columns are printed as integers
*
* \example: cout << FormatVertVectorHKLFloats<double>(vH,vK,vL,vIobs,vIcalc,vSigma,12,4);
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
      ~FormatVertVectorHKLFloats();
      //int length() const;
   //private:
      const CrystVector<T> *mpVectors[20];
      const int mNbVectors;
      const int mWidth;
      const int mPrecision;
   
   friend ostream& operator<< <T>(ostream&, const FormatVertVectorHKLFloats<T>&);
};

template<class T> ostream& operator<< (ostream& os, const FormatVertVectorHKLFloats<T> &fStr);


#endif

