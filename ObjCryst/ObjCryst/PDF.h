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
/*   PDF.h header file for Pair Distribution calculations
*
*/
#ifndef _OBJCRYST_PDF_H_
#define _OBJCRYST_PDF_H_

#include "CrystVector/CrystVector.h"

#include "ObjCryst/Crystal.h"
#ifdef __WX__CRYST__
   #include "wxCryst/wxCryst.h"
#endif

#include <string>
#include <list>
#include <vector>

namespace ObjCryst
{
extern const RefParType *gpRefParTypePDF;

// Forward declaration
class PDFPhase;


/** Main class for Pair distribution function calculations and comparison to observed one
*
*/
class PDF:public RefinableObj
{
   public:
      PDF();
      /// Copy constructor
      PDF(const PDF &old);
      /// Crystal destructor
      ~PDF();
      virtual const string& GetClassName() const;
      /// Experimental data
      void SetPDFObs(const CrystVector_REAL &r,const CrystVector_REAL &obs);
      /// GetRadiationType
      RadiationType GetRadiationType()const;
      /// SetRadiationType
      void SetRadiationType(RadiationType type);
      /// Cutoff r value
      REAL GetRMax()const;
      /// Get the calculated PDF
      const CrystVector_REAL &GetPDFCalc()const;
      /// Get the observed PDF
      const CrystVector_REAL &GetPDFObs()const;
      /// Get the r coordinates (in Angstroem) for the PDF
      const CrystVector_REAL &GetPDFR()const;
      /// Add PDF phase
      void AddPDFPhase(PDFPhase &phase);
   private:
      /// Radiation type
      RadiationType mRadiationType;
      /// r coordinates (in Angstroem) of the PDF
      CrystVector_REAL mPDFR;
      /// The PDF itself,calculated
      mutable CrystVector_REAL mPDFCalc;
      /// The PDF itself,observed
      CrystVector_REAL mPDFObs;
      /// Contributions to the PDF, with their scale factor
      std::list<std::pair<PDFPhase*,REAL> > mvPDFPhase;
   #ifdef __WX__CRYST__
   public:
      virtual WXCrystObjBasic* WXCreate(wxWindow*);
   #endif
};

/// Global registry for all PDF objects
extern ObjRegistry<PDF> gPDFRegistry;

/** Contribution to a PDF
*
* Right now this can only be a PDFCrystal.
*/
class PDFPhase:public RefinableObj
{
   public:
      /// Constructor
      PDFPhase(const PDF &pdf);
      /// Get the calculated PDF
      const CrystVector_REAL &GetPDFCalc()const;
   protected:
      /// Calculate the pdf
      virtual void CalcPDF()const=0;
      /// Parent PDF
      const PDF *mpPDF;
      /// The calculated PDF
      mutable CrystVector_REAL mPDFCalc;
   #ifdef __WX__CRYST__
   public:
      virtual WXCrystObjBasic* WXCreate(wxWindow*)=0;
   #endif
};

/** Class for Pair Distribution Function calculations for a single Crystal object
*
*/
class PDFCrystal:public PDFPhase
{
   public:
      /// Constructor
      PDFCrystal(const PDF &pdf, const Crystal &cryst);
   private:
      /// Initialize all parameters
      void Init(const PDF &pdf, const Crystal &cryst);
      /// Calculate the pdf
      virtual void CalcPDF()const;
      /// The Crystal
      const Crystal *mpCrystal;
      // Parameters to describe the PDF
         /// Parameters affecting the width & intensity of peaks
         REAL mDelta1,mDelta2,mQbroad,mQdamp;
      // Cached data to speed up computations
         /// Container of temp data for each atom
         struct pdfAtom
         {
            /// default constructor
            pdfAtom();
            /// Fractionnal & cartesian unique coordinates
            REAL fx0,fy0,fz0,x0,y0,z0;
            /// The scattering power
            const ScatteringPower *pScattPow;
            /// Scattering amplitudes, multiplied by occupancy
            REAL occupBi;
            /// Has this atom changed since last time ?
            bool hasChanged;
            /// List of all equivalent & translated positions,
            /// in cartesian coordinates
            CrystVector_REAL x,y,z;
         };
         /// List of all temp data
         mutable std::vector<pdfAtom> mvPDFAtom;
   #ifdef __WX__CRYST__
   public:
      virtual WXCrystObjBasic* WXCreate(wxWindow*);
   #endif
};
}//namespace
#endif
