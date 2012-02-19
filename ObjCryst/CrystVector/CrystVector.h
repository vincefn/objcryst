/*  ObjCryst++ Object-Oriented Crystallographic Library
    (c) 2000-2008 Vincent Favre-Nicolin vincefn@users.sourceforge.net
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
// This header defines which description should be used
// for vectors.

#ifndef __LIBCRYST_VECTOR_H
#define __LIBCRYST_VECTOR_H

#include "ObjCryst/ObjCryst/General.h"

#undef __LIBCRYST_VECTOR_USE_BLITZ__
//#define __LIBCRYST_VECTOR_USE_BLITZ__

#ifdef __LIBCRYST_VECTOR_USE_BLITZ__

#include <blitz/array.h>
using namespace blitz;

//Still use pointers for the Geometrical Structure computation ?
//(due to the complexity of the formulas, a blitz coding requires
//to much memory when compiling)
#define __VFN_GEOM_STRUCT_FACTOR_USE_POINTERS

#define CrystVector_REAL Array<REAL,1>
#define CrystVector_float  Array<float,1>
#define CrystVector_long   Array<long,1>
#define CrystVector_int    Array<int,1>
#define CrystVector_uint   Array<unsigned int,1>
#define CrystVector_bool   Array<bool,1>

#define CrystMatrix_REAL Array<REAL,2>
#define CrystMatrix_float  Array<float,2>
#define CrystMatrix_long   Array<long,2>
#define CrystMatrix_int    Array<int,2>
#define CrystMatrix_uint   Array<unsigned int,2>
#define CrystMatrix_bool   Array<bool,2>

#define CrystVector_T    Array<T,1> 
#define CrystMatrix_T    Array<T,2> 
#define CrystArray3D_T    Array<T,3> 

template<class T> T MaxDifference(const Array<T,1> &a,const Array<T,1> &b);

template<class T> T MaxDifference(const Array<T,2> &a,const Array<T,2> &b);

#else  // __VFN_VECTOR_USE_BLITZ__


#define CrystVector_REAL CrystVector<REAL>
#define CrystVector_double  CrystVector<double>
#define CrystVector_float  CrystVector<float>
#define CrystVector_long   CrystVector<long>
#define CrystVector_int    CrystVector<int>
#define CrystVector_uint   CrystVector<unsigned int>
#define CrystVector_bool   CrystVector<bool>

#define CrystMatrix_REAL CrystMatrix<REAL>
#define CrystMatrix_double CrystMatrix<double>
#define CrystMatrix_float  CrystMatrix<float>
#define CrystMatrix_long   CrystMatrix<long>
#define CrystMatrix_int    CrystMatrix<int>
#define CrystMatrix_uint   CrystMatrix<unsigned int>
#define CrystMatrix_bool   CrystMatrix<bool>

#define CrystArray3D_REAL CrystArray3D<REAL>

#define CrystVector_T    CrystVector<T> 
#define CrystMatrix_T    CrystMatrix<T> 
#define CrystArray3D_T    CrystMatrix<T> 

#define __VFN_GEOM_STRUCT_FACTOR_USE_POINTERS

#include <iostream>
#include <cmath>
using namespace std;

#ifdef _MSC_VER // MS VC++ predefined macros....
#undef min
#undef max
#endif

//######################################################################
//  CrystVector
//######################################################################
/** \brief Vector library (Blitz++ mimic) for ObjCryst++
*  
*  The CrystVector library is not a new array computation library, despite
* the appearances. ObjCryst++ should used the 
*  <a href="http://www.oonumerics.org/blitz/"> Blitz++ array library</a> , which yields
* excellent performance \e and simple array expressions. Unfortunately, the memory
* required to \e compile the library using gcc is far too high to be reasonable
* when using complex expressions and optimizing code. So until this has changed,
* The CrystVector and CrystMatrix library have been created, and these emulate
* (supposedly exactly) the Blitz++ interface (but not the smart handling of
* mathematical expressions, so pointers must be used). For documentation about these
* two libraries you should read the <a href="http://www.oonumerics.org/blitz/manual/"> 
* Blitz++ documentation</a>. CrystVector and CrystMatrix use the same kind of storage
* in memory.
*
* You can use CrystVector_REAL, CrystVector_long,etc... to declare 1D vectors. Macros
* ensure (well, should) ensure compatibility with Blitz++. (as of april 2001 support of
* blitz++ is broken).
*
* \todo check again Blitz++ support in newer gcc versions.
*/
template<class T> class CrystVector
{
   public:
   CrystVector();
   
   CrystVector(const long nbElements);
   
   CrystVector(const CrystVector &old);
   
   ~CrystVector();

   void operator=(const CrystVector &old);
   
   template<class U> void operator=(const CrystVector<U> &old)
   {
      if(mNumElements != old.numElements())
      {
         mNumElements = old.numElements();
         delete[] mpData;
         mpData=new T[mNumElements];
      };
      mIsAreference=false;
      T *p1=mpData;
      const U *p2=old.data();
      for(long i=0;i<mNumElements;i++) *p1++ = (T) *p2++;
   }
   
   #ifdef __MWERKS__
   operator CrystVector<bool>() const
   {
      CrystVector<bool> vect;
      vect=*this;
      return vect;
   }
   operator CrystVector<int>() const
   {
      CrystVector<int> vect;
      vect=*this;
      return vect;
   }
   operator CrystVector<long>() const
   {
      CrystVector<long> vect;
      vect=*this;
      return vect;
   }
   operator CrystVector<float>() const
   {
      CrystVector<float> vect;
      vect=*this;
      return vect;
   }
   operator CrystVector<double>() const
   {
      CrystVector<double> vect;
      vect=*this;
      return vect;
   }
   #else
   template<class U> operator CrystVector<U>() const
   {
      CrystVector<U> vect;
      vect=*this;
      return vect;
   }
   #endif
   /** Define this vector as a reference to another vector, or part of another vector
   *
   * \param old: the vector this one should point to
   * \param imin,imax: imax>imin, this vector will point to a part of the other
   * vector, from old(i) to old(imax-1)
   */
   void reference(CrystVector &old, const long imin=0, const long imax=0);
   
   long numElements()const;
   long size()const;
   T sum()const;
   T min()const;
   T max()const;
   /// Find index of minimum, between start and end (if start==end, use full vector)
   unsigned long imin(const unsigned long start=0,const unsigned long finish=0)const;
   /// Find index of maximum, between start and end (if start==end, use full vector)
   unsigned long imax(const unsigned long start=0,const unsigned long finish=0)const;
   
   T * data();
   const T * data() const;
   
   void resize(const long newNbElements);
   
   void resizeAndPreserve(const long newNbElements);
   
   void operator*=(const T num);
   
   void operator*=(const CrystVector &vect);
   
   void operator/=(const T num);
   
   void operator/=(const CrystVector &vect);
   
   void operator+=(const T num);
   
   void operator+=(const CrystVector &vect);
   
   void operator-=(const CrystVector &vect);
/* Buggy ? (my version, not Blitz's!)
   // ListInitializer & ListInitializerSwitch are a simplified
   //version borrowed from the blitz++ library (see the blitz/listinit.h file)
   //
   // It allows a very convenient way of initializing arrays like this:
   // CrystVector arr(5); arr = 1,3,6,8,9;
   class ListInitializer
   {
      public:
         ListInitializer(T *pData):mpData(pData){};
         ~ListInitializer(){cout << "toto";};
         ListInitializer operator,(T x)
         {
            *mpData=x;
            return ListInitializer(mpData+1);
         }
      private:
         ListInitializer(){};
         T *mpData;
   };
   class ListInitializerSwitch
   {
      public:
         ListInitializerSwitch(CrystVector &vect, const T value):
            mVectData(vect.data()),mValue(value),
            mNumElements(vect.numElements()),wipeOnDestruct(true)
         {};

         ~ListInitializerSwitch()
         {
            if(wipeOnDestruct)
            {  //only one value given -> set all elements
               for(int i=0;i<mNumElements;i++)*mVectData++ = mValue;
            };
         };

         ListInitializer operator,(T x)
         {
            wipeOnDestruct=false;//operator, is used
            *mVectData=mValue;
            mVectData++;
            *mVectData=x;
            return ListInitializer(mVectData+1);
         }
      private:
         T *mVectData;
         T mValue;
         long mNumElements;
         bool wipeOnDestruct;
   };
   
   
   ListInitializerSwitch operator=(const T num)
   {
      return ListInitializerSwitch(*this,num);
   }
*/
   void operator=(const T num);

   T operator()(const long i) const;

   T& operator()(const long i);

   
   protected:
   private:
   T *mpData;
   long mNumElements;
   bool mIsAreference;//is a reference to another vector ?
   //friend ostream& operator<<(ostream &os,const CrystVector &vect);

};

template<class T> ostream& operator<<(ostream &os, CrystVector<T> &vect);

//Return the sorted subscripts of the array
template<class T> CrystVector<long> SortSubs(const CrystVector<T> &vect);
//Sort the array in place and also return the sorted subscripts
template<class T> long QuickSortSubs(CrystVector<T> &vect,
                                     CrystVector<long> &subscript,
                                     long last,long first=0, int depth=0);
                                     
///Cosinus (slow routine, not memory-savy..)
template<class T> CrystVector<T> cos(const CrystVector<T> &vect);
///Sinus (slow routine, not memory-savy...)
template<class T> CrystVector<T> sin(const CrystVector<T> &vect);
///Tangent (slow routine, not memory-savy...)
template<class T> CrystVector<T> tan(const CrystVector<T> &vect);
///Square root (slow routine, not memory-savy...)
template<class T> CrystVector<T> sqrt(const CrystVector<T> &vect);

//######################################################################
//  CrystMatrix
//######################################################################
/** \brief 2D Vector library (Blitz++ mimic) for ObjCryst++
*  
*  The CrystVector library is not a new array computation library, despite
* the appearances. ObjCryst++ should used the 
*  <a href="http://www.oonumerics.org/blitz/"> Blitz++ array library</a> , which yields
* excellent performance \e and simple array expressions. Unfortunately, the memory
* required to \e compile the library using gcc is far too high to be reasonable
* when using complex expressions and optimizing code. So until this has changed,
* The CrystVector and CrystMatrix library have been created, and these emulate
* (supposedly exactly) the Blitz++ interface (but not the smart handling of
* mathematical expressions, so pointers must be used). For documentation about these
* two libraries you should read the <a href="http://www.oonumerics.org/blitz/manual/"> 
* Blitz++ documentation</a>. CrystVector and CrystMatrix use the same kind of storage
* in memory.
*
* You can use CrystMatrix_REAL, CrystMatrix_long,etc... to declare 2D vectors. Macros
* ensure (well, should) ensure compatibility with Blitz++. (as of april 2001 support of
* blitz++ is broken).
*/
template<class T> class CrystMatrix
{
   public:
   CrystMatrix();
   //ySize : number of rows, xSize : number of columns
   CrystMatrix(const long ySize,const long xSize);
   
   CrystMatrix(const CrystMatrix &old);
   
   ~CrystMatrix();
   
   void operator=(const CrystMatrix &old);

   void reference(CrystMatrix &old);
   long numElements()const;
   long size()const;
   T sum()const;
   T min()const;
   T max()const;
   long rows()const;
   long cols()const;
   
   T * data();
   const T * data() const;
   
   void resize(const long ySize,const long xSize);

   void resizeAndPreserve(const long ySize,const long xSize);
   
   //void operator=(const T num);
   /// Element-by element multiplication (array-like)
   void operator*=(const T num);
   /// Element-by element multiplication (array-like)
   void operator*=(const CrystMatrix &vect);
   /// Element-by element division (array-like)
   void operator/=(const T num);
   /// Element-by element addition (array-like)
   void operator+=(const T num);
   /// Element-by element substraction (array-like)
   void operator-=(const T num);
   
   /// matrix multiplication (linear algebra)
   CrystMatrix Mult(const CrystMatrix &rhs);
   //:TODO: Check the following...
   
   // ListInitializer & ListInitializerSwitch are a simplified
   //version borrowed from the blitz++ library (see the blitz/listinit.h file)
   //
   // It allows a very convenient way of initializing arrays like this:
   // CrystVector arr(5); arr = 1,3,6,8,9;
   class ListInitializer
   {
      public:
         ListInitializer(T *pData):mpData(pData){};
         ~ListInitializer(){};
         ListInitializer operator,(T x)
         {
            *mpData=x;
            return ListInitializer(mpData+1);
         }
      private:
         ListInitializer(){};
         T *mpData;
   };
   class ListInitializerSwitch
   {
      public:
         ListInitializerSwitch(CrystMatrix &vect, const T value):
            mVectData(vect.data()),mValue(value),
            mNumElements(vect.numElements()),wipeOnDestruct(true)
         {};

         ~ListInitializerSwitch()
         {
            if(wipeOnDestruct)
            {  //only one value given -> set all elements
               for(int i=0;i<mNumElements;i++)*mVectData++ = mValue;
            };
         };

         ListInitializer operator,(T x)
         {
            wipeOnDestruct=false;//operator, is used
            *mVectData=mValue;
            mVectData++;
            *mVectData=x;
            return ListInitializer(mVectData+1);
         }
      private:
         T *mVectData;
         T mValue;
         long mNumElements;
         bool wipeOnDestruct;
   };
   
   
   ListInitializerSwitch operator=(const T num)
   {
      return ListInitializerSwitch(*this,num);
   }
   
   
   T operator()(const long i) const;

   T operator()(const long row,const long col) const;
   
   T& operator()(const long i);

   T& operator()(const long i,const long j);

   CrystMatrix transpose()const;
   
   protected:
   private:
   T *mpData;
   long mNumElements;
   long mXSize,mYSize;
   bool mIsAreference;//is a reference to another vector ?
   
   //friend ostream& operator<<(ostream &os,const CrystMatrix &vect);

};

template<class T> ostream& operator<<(ostream &os, const CrystMatrix<T> &vect);

template<class T> T MaxDifference(const CrystVector<T> &a,const CrystVector<T> &b);

template<class T> T MaxDifference(const CrystMatrix<T> &a,const CrystMatrix<T> &b);

template<class T> CrystMatrix<T> product(const CrystMatrix<T> &a,const CrystMatrix<T> &b);

//######################################################################
//  CrystArray3D
//######################################################################
/** \brief 3D Vector (Blitz++ mimic) for ObjCryst++
*  
*  The CrystVector library is not a new array computation library, despite
* the appearances. ObjCryst++ should used the 
*  <a href="http://www.oonumerics.org/blitz/"> Blitz++ array library</a> , which yields
* excellent performance \e and simple array expressions. Unfortunately, the memory
* required to \e compile the library using gcc is far too high to be reasonable
* when using complex expressions and optimizing code. So until this has changed,
* The CrystVector and CrystMatrix library have been created, and these emulate
* (supposedly exactly) the Blitz++ interface (but not the smart handling of
* mathematical expressions, so pointers must be used). For documentation about these
* two libraries you should read the <a href="http://www.oonumerics.org/blitz/manual/"> 
* Blitz++ documentation</a>. CrystVector and CrystMatrix use the same kind of storage
* in memory.
*
* You can use CrystArray3D_REAL, CrystArray3D_long,etc... to declare 3D vectors. Macros
* ensure (well, should) ensure compatibility with Blitz++. (as of april 2001 support of
* blitz++ is broken).
*/
template<class T> class CrystArray3D
{
   public:
   CrystArray3D();
   //ySize : number of rows, xSize : number of columns
   CrystArray3D(const long zSize,const long ySize,const long xSize);
   
   CrystArray3D(const CrystArray3D &old);
   
   ~CrystArray3D();
   
   void operator=(const CrystArray3D &old);

   void reference(CrystArray3D &old);
   long numElements()const;
   long size()const;
   T sum()const;
   T min()const;
   T max()const;
   long rows()const;
   long cols()const;
   long depth()const;
   
   T * data();
   const T * data() const;
   
   void resize(const long zSize,const long ySize,const long xSize);

   void resizeAndPreserve(const long zSize,const long ySize,const long xSize);
   
   void operator=(const T num);
   void operator*=(const T num);
   void operator*=(const CrystArray3D &vect);
   void operator/=(const T num);
   void operator+=(const T num);
   void operator-=(const T num);
      
   //void operator=(const T num);
   
   T operator()(const long i) const;

   T operator()(const long depth,const long row,const long col) const;
   
   T& operator()(const long i);

   T& operator()(const long depth,const long row,const long col);

   protected:
   private:
   T *mpData;
   long mNumElements;
   long mXSize,mYSize,mZSize;
   bool mIsAreference;//is a reference to another vector ?
};
template<class T> ostream& operator<<(ostream &os, const CrystArray3D<T> &vect);

#endif // __LIBCRYST_VECTOR_USE_BLITZ__

//Basic Gauss-Jordan elimination with partial pivot (rows only, using max pivot)
//Definitly *not* optimized !
CrystMatrix_REAL InvertMatrix(const CrystMatrix_REAL &m);
template<class T> void MatrixExchangeRows(CrystMatrix_T &m, const long row1, const long row2);

///Maximum absolute value of vector
template<class T> T MaxAbs(const CrystVector_T &vector);
///Minimum absolute value of vector
template<class T> T MinAbs(const CrystVector_T &vector);

//######################################################################
//  CubicSpline
//######################################################################
/// Cubic spline interpolation.
class CubicSpline
{
   public:
      /// Default constructor - CubicSpline::Init() should be called afterwards
      CubicSpline();
      /// Spline with given extremum derivatives
      CubicSpline(const CrystVector_REAL &x, const CrystVector_REAL &y, 
                  const REAL yp1, const REAL ypn);
      /// Spline with given extremum derivatives
      CubicSpline(const REAL *px, const REAL *py, const unsigned long nbPoints,
                  const REAL yp1, const REAL ypn);
      /// Natural cubic spline
      CubicSpline(const CrystVector_REAL &x, const CrystVector_REAL &y);
      /// Natural cubic spline
      CubicSpline(const REAL *px, const REAL *py, const unsigned long nbPoints);
      ~CubicSpline();
      
      /// Spline with given extremum derivatives
      void Init(const CrystVector_REAL &x, const CrystVector_REAL &y, 
                const REAL yp1, const REAL ypn);
      /// Spline with given extremum derivatives
      void Init(const REAL *px, const REAL *py, const unsigned long nbPoints,
                const REAL yp1, const REAL ypn);
      /// Natural cubic spline
      void Init(const CrystVector_REAL &x, const CrystVector_REAL &y);
      /// Natural cubic spline
      void Init(const REAL *px, const REAL *py, const unsigned long nbPoints);
      
      /// Get spline value at a series of point - x is assumed to be sorted by increasing values
      CrystVector_REAL operator()(const CrystVector_REAL &x) const;
      /// Get spline value on a range of values with a fixed step
      CrystVector_REAL operator()(const REAL min,const REAL step, const long nbpoint) const;
      /// Get spline value at one point
      REAL operator()(const REAL x) const;
   private:
      void InitSpline(const REAL yp1, const REAL ypn);
      void InitNaturalSpline();
      CrystVector_REAL mX;
      CrystVector_REAL mY;
      CrystVector_REAL mYsecond;
};
//######################################################################
//  Savitzky-Golay interpolation
//######################################################################
/** Savitzky-Golay computing of smoothed data, using 2nd order polynomial coefficients
*
*/
CrystVector_REAL SavitzkyGolay(const CrystVector_REAL &v, const unsigned int m, const unsigned int deriv);

#endif   // __LIBCRYST_VECTOR_H
