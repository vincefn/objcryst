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
#ifndef _OBJCRYST_SCATTERING_CORR_H_
#define _OBJCRYST_SCATTERING_CORR_H_

#include "ObjCryst/ScatteringData.h"

namespace ObjCryst
{

/** Base class to compute all kind of corrections to intensities:  Lorentz, Polar,
* absorption, texcture, extinction, etc...
*
* The computed intensities are to be multiplied by all the ScatteringCorr calculated.
* 
*
* This is an abstract base class.
*/
class ScatteringCorr
{
   public:
      /// Constructor, with the associated ScatteringData object.
      ScatteringCorr(const ScatteringData & data);
      virtual ~ScatteringCorr();
      /// Get the name of this object
      virtual const string & GetName() const=0;
      /// Get the name of the class
      virtual const string & GetClassName() const=0;
      /// Get the vector of corrections for all reflections. Calculated values must
      /// be multiplied by these values.
      const CrystVector_REAL& GetCorr() const;
      /// Get the value of the clock corresponding to the last time the correction
      /// was actually computed
      const RefinableObjClock& GetClockCorr()const;
   protected:
      /// Do the computation of corrected intensities
      virtual void CalcCorr() const=0;
      /// The associated ScatteringData object
      const ScatteringData *mpData;
      /// The vector of correction to intensities.
      mutable CrystVector_REAL mCorr;
      /// The clock marking the last time the correction was calculated
      mutable RefinableObjClock mClockCorrCalc;
};

/** Lorentz Correction
*
* So far, it only considers the correction for equatorial diffraction:
* \f$ L = \frac{1}{\sin(2\theta)} \f$
*/
class LorentzCorr:public ScatteringCorr
{
   public:
      LorentzCorr(const ScatteringData & data);
      virtual ~LorentzCorr();
      virtual const string & GetName() const;
      virtual const string & GetClassName() const;
   protected:
      virtual void CalcCorr() const;
};

/** Polarization Correction
*
* So far, it only considers the correction for equatorial diffraction:
*\f$ P = \frac{1}{1+A}\left(1+A\cos^2(2\theta)\right) \f$ (Polarization factor), with
* \f$ A = \frac{1-f}{1+f} \f$, where f is the polarization rate of the incident
*beam in the plane which (i) includes the incident beam, and (ii) is perpendicular to  
*the diffracting plane. For an X-Ray Tube without monochromator, A=1, and 
*if there is a monochromator : \f$ A = \cos^2(2\theta_{mono}) \f$
*
* Currently, the linear polarization factor is directly read from the radiation object,
* and the linear polarization (if any) is assumed to be perpendicular to the diffracting
* plane (standard synchrotron geometry).
*
* \todo: extend this to take into account other diffracting & monochromatic geometries.
*/
class PolarizationCorr:public ScatteringCorr
{
   public:
      PolarizationCorr(const ScatteringData & data);
      virtual ~PolarizationCorr();
      virtual const string & GetName() const;
      virtual const string & GetClassName() const;
   protected:
      virtual void CalcCorr() const;
      mutable REAL mPolarAfactor;
};

/** Slit aperture correction (for powder)
*
* This correction takes into account the fact that diffraction
* rings (cones) have a portion of the ring proportionnal to 
* \f$ SlitAp = \frac{1}{\sin(\theta)} \f$ which falls into the detector
* (due to slits in the direction perpendicular to the incident beam/ detector plane).
*/
class PowderSlitApertureCorr:public ScatteringCorr
{
   public:
      PowderSlitApertureCorr(const ScatteringData & data);
      virtual ~PowderSlitApertureCorr();
      virtual const string & GetName() const;
      virtual const string & GetClassName() const;
   protected:
      virtual void CalcCorr() const;
};

/** Texture correction using the March-Dollase model.
*
* This can include several phases.
*
* \todo the actual implementation of the texture classes...
*/
class TextureMarchDollaseCorr:public ScatteringCorr
{
   public:
      TextureMarchDollaseCorr(const ScatteringData & data);
      virtual ~TextureMarchDollaseCorr();
      virtual const string & GetName() const;
      virtual const string & GetClassName() const;
   protected:
      virtual void CalcCorr() const;
};
}//namespace
#endif //_OBJCRYST_SCATTERING_CORR_H_
