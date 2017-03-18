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
#ifndef _OBJCRYST_SCATTERINGDATA_H_
#define _OBJCRYST_SCATTERINGDATA_H_

#include "ObjCryst/CrystVector/CrystVector.h"

#include "ObjCryst/ObjCryst/General.h"

#include "ObjCryst/ObjCryst/SpaceGroup.h"
#include "ObjCryst/ObjCryst/ScatteringPower.h"
#include "ObjCryst/ObjCryst/Scatterer.h"
#include "ObjCryst/ObjCryst/Crystal.h"

//#include <stdlib.h>
#include <string>
//#include <iomanip>
//#include <cmath>
//#include <typeinfo>
//#include <fstream>
//#include <ctime>

namespace ObjCryst
{
#ifndef HAVE_SSE_MATHFUN
//initialize tabulated values of cosine
void InitLibCrystTabulCosine();
void DeleteLibCrystTabulCosine();
void InitLibCrystTabulExp();
void DeleteLibCrystTabulExp();
#endif
/// Generic type for scattering data
extern const RefParType *gpRefParTypeScattData;
/// Type for scattering data scale factors
extern const RefParType *gpRefParTypeScattDataScale;
/// Type for reflection profile
extern const RefParType *gpRefParTypeScattDataProfile;
/// Type for reflection profiles type (e.g. gaussian/lorentzian mix)
extern const RefParType *gpRefParTypeScattDataProfileType;
/// Type for reflection profile width
extern const RefParType *gpRefParTypeScattDataProfileWidth;
/// Type for reflection profile asymmetry
extern const RefParType *gpRefParTypeScattDataProfileAsym;
/// Generic type for scattering data correction parameter
extern const RefParType *gpRefParTypeScattDataCorr;
/// Generic type for correction to calculated intensities
extern const RefParType *gpRefParTypeScattDataCorrInt;
/// Parameter type for preferred orientation direction
extern const RefParType *gpRefParTypeScattDataCorrIntPO_Direction;
/// Parameter type for fraction of preferred orientation
extern const RefParType *gpRefParTypeScattDataCorrIntPO_Fraction;
/// Parameter type for the amplitude of preferred orientation
extern const RefParType *gpRefParTypeScattDataCorrIntPO_Amplitude;
/// Parameter type for the ellipsoid coefficient
extern const RefParType *gpRefParTypeScattDataCorrInt_Ellipsoid;
/// Parameter type for absorption correction
extern const RefParType *gpRefParTypeScattDataCorrIntAbsorp;
/// Parameter type for polarization correction
extern const RefParType *gpRefParTypeScattDataCorrIntPolar;
/// Parameter type for extinction correction
extern const RefParType *gpRefParTypeScattDataCorrIntExtinc;
/// Parameter type for correction to peak positions
extern const RefParType *gpRefParTypeScattDataCorrPos;
/// Parameter type for background intensity
extern const RefParType *gpRefParTypeScattDataBackground;

extern const RefParType *gpRefParTypeRadiation;
extern const RefParType *gpRefParTypeRadiationWavelength;

class NiftyStaticGlobalObjectsInitializer_ScatteringData
{
   public:
      NiftyStaticGlobalObjectsInitializer_ScatteringData()
      {
         if (mCount++ == 0)
         {
            #ifndef HAVE_SSE_MATHFUN
            InitLibCrystTabulCosine();
            InitLibCrystTabulExp();
            #endif
            gpRefParTypeScattData= new RefParType(gpRefParTypeObjCryst,"Scattering Data");
            gpRefParTypeScattDataScale= new RefParType(gpRefParTypeObjCryst,"Scale Factor");
            gpRefParTypeScattDataProfile= new RefParType(gpRefParTypeScattData,"Profile");
            gpRefParTypeScattDataProfileType= new RefParType(gpRefParTypeScattDataProfile,"Type");
            gpRefParTypeScattDataProfileWidth= new RefParType(gpRefParTypeScattDataProfile,"Width");
            gpRefParTypeScattDataProfileAsym= new RefParType(gpRefParTypeScattDataProfile,"Asymmetry");
            gpRefParTypeScattDataCorr= new RefParType(gpRefParTypeScattData,"Correction");
            gpRefParTypeScattDataCorrInt= new RefParType(gpRefParTypeScattDataCorr,"Intensities");
            gpRefParTypeScattDataCorrIntPO_Direction= new RefParType(gpRefParTypeScattDataCorrIntPO_Direction,"Preferred orientation direction");
            gpRefParTypeScattDataCorrIntPO_Fraction= new RefParType(gpRefParTypeScattDataCorrIntPO_Fraction,"Preferred orientation fraction");
            gpRefParTypeScattDataCorrIntPO_Amplitude= new RefParType(gpRefParTypeScattDataCorrIntPO_Amplitude,"Preferred orientation amplitude");
            gpRefParTypeScattDataCorrInt_Ellipsoid= new RefParType(gpRefParTypeScattDataCorrInt_Ellipsoid,"Preferred orientation ellipsoid");
            gpRefParTypeScattDataCorrIntAbsorp= new RefParType(gpRefParTypeScattDataCorrInt,"Absorption");
            gpRefParTypeScattDataCorrIntPolar= new RefParType(gpRefParTypeScattDataCorrInt,"Polarization");
            gpRefParTypeScattDataCorrIntExtinc= new RefParType(gpRefParTypeScattDataCorrInt,"Extinction");
            gpRefParTypeScattDataCorrPos= new RefParType(gpRefParTypeScattDataCorr,"Reflections Positions");
            gpRefParTypeScattDataBackground= new RefParType(gpRefParTypeScattData,"Background");
            gpRefParTypeRadiation= new RefParType(gpRefParTypeObjCryst,"Radiation");
            gpRefParTypeRadiationWavelength= new RefParType(gpRefParTypeRadiation,"Wavelength");
         }
      }
      ~NiftyStaticGlobalObjectsInitializer_ScatteringData()
      {
         if (--mCount == 0)
         {
            #ifndef HAVE_SSE_MATHFUN
            DeleteLibCrystTabulCosine();
            DeleteLibCrystTabulExp();
            #endif
            delete gpRefParTypeScattData;
            delete gpRefParTypeScattDataScale;
            delete gpRefParTypeScattDataProfile;
            delete gpRefParTypeScattDataProfileType;
            delete gpRefParTypeScattDataProfileWidth;
            delete gpRefParTypeScattDataProfileAsym;
            delete gpRefParTypeScattDataCorr;
            delete gpRefParTypeScattDataCorrInt;
            delete gpRefParTypeScattDataCorrIntAbsorp;
            delete gpRefParTypeScattDataCorrIntPolar;
            delete gpRefParTypeScattDataCorrIntExtinc;
            delete gpRefParTypeScattDataCorrPos;
            delete gpRefParTypeScattDataBackground;
            delete gpRefParTypeRadiation;
            delete gpRefParTypeRadiationWavelength;
            gpRefParTypeScattData=0;
            gpRefParTypeScattDataScale=0;
            gpRefParTypeScattDataProfile=0;
            gpRefParTypeScattDataProfileType=0;
            gpRefParTypeScattDataProfileWidth=0;
            gpRefParTypeScattDataProfileAsym=0;
            gpRefParTypeScattDataCorr=0;
            gpRefParTypeScattDataCorrInt=0;
            gpRefParTypeScattDataCorrIntAbsorp=0;
            gpRefParTypeScattDataCorrIntPolar=0;
            gpRefParTypeScattDataCorrIntExtinc=0;
            gpRefParTypeScattDataCorrPos=0;
            gpRefParTypeScattDataBackground=0;
            gpRefParTypeRadiation=0;
            gpRefParTypeRadiationWavelength=0;
         }
      }
   private:
      static long mCount;
};
static NiftyStaticGlobalObjectsInitializer_ScatteringData NiftyStaticGlobalObjectsInitializer_ScatteringData_counter;
//######################################################################
/** \brief Class to define the radiation (type, monochromaticity, wavelength(s)) of an experiment
*
* This can be developped for more complex experiments, hence the \e vector of
* wavelengths (so far it is not possible to use several wavelengths, though).
*
* X-Rays and Neutrons are available. Electrons are not available yet in
* ScatteringData classes.
*
* \todo also add here information about the polarization of the beam.
*/
//######################################################################
class Radiation: public RefinableObj
{
   public:
      /// Default constructor
      Radiation();
      /** \ brief Constructor
      *
      * \param rad the RadiationType used (X-Rays, neutrons)
      * \param wavelength the wavelength (in Angstroems) of the monochromatic
      * radiation.
      */
      Radiation(const RadiationType rad,const REAL wavelength);
      /** \ brief Constructor for X-Ray tube radiation
      *
      *\param XRayTubeElementName : name of the anticathode element name. Known
      *ones are Cr, Fe, Cu, Mo, Ag.
      *\param alpha2Alpha2ratio: Kalpha2/Kalpha1 ratio (0.5 by default)
      *
      *the average wavelength is calculated
      *using the alpha2/alpha1 weight. All structure factors computation are made
      *using the average wavelength, and for powder diffraction, profiles are output
      *at the alpha1 and alpha2 ratio for the calculated pattern.
      *
      *NOTE : if the name of the wavelength is generic (eg"Cu"),
      *then the program considers that
      *there are both Alpha1 and Alpha2, and thus automatically changes the WavelengthType
      *to WAVELENGTH_ALPHA12. If instead either alpha1 or alpha2 (eg "CuA1") is asked for,
      *the WavelengthType is set to WAVELENGTH_MONOCHROMATIC. In both cases,
      * the radiation type is set to X-Ray.
      */
      Radiation(const string &XRayTubeElementName,const REAL alpha2Alpha2ratio=0.5);
      /// Copy constructor
      Radiation(const Radiation&);
      ~Radiation();
      virtual const string& GetClassName() const;

      void operator=(const Radiation&);

      /// Get the radiation type (X-Rays, Neutron)
      RadiationType GetRadiationType()const;
      /// Set the radiation type (X-Rays, Neutron)
      void SetRadiationType(const RadiationType);
      /// Set the Wavelength type (monochromatic, Alpha1+Alpha2, Time Of Flight...)
      void SetWavelengthType(const WavelengthType &type);
      /// Get the Wavelength type (monochromatic, Alpha1+Alpha2, Time Of Flight...)
      WavelengthType GetWavelengthType()const;
      /// Get the wavelength(s) in Angstroems. Currently only
      /// monochromatic is used, so the vector should only return
      /// only one wavelength.
      const CrystVector_REAL& GetWavelength()const;
      /// Set the (monochromatic) wavelength of the beam.
      void SetWavelength(const REAL );
      /** \ brief Set X-Ray tube radiation.
      *
      *\param XRayTubeElementName : name of the anticathode element name. Known
      *ones are Cr, Fe, Cu, Mo, Ag.
      *\param alpha2Alpha2ratio: Kalpha2/Kalpha1 ratio (0.5 by default)
      *
      *the average wavelength is calculated
      *using the alpha2/alpha1 weight. All structure factors computation are made
      *using the average wavelength, and for powder diffraction, profiles are output
      *at the alpha1 and alpha2 ratio for the calculated pattern.
      *
      *NOTE : if the name of the wavelength is generic (eg"Cu"),
      *then the program considers that
      *there are both Alpha1 and Alpha2, and thus automatically changes the WavelengthType
      *to WAVELENGTH_ALPHA12. If instead either alpha1 or alpha2 (eg "CuA1") is asked for,
      *the WavelengthType is set to WAVELENGTH_MONOCHROMATIC. In both cases,
      * the radiation type is set to X-Ray.
      */
      void SetWavelength(const string &XRayTubeElementName,const REAL alpha2Alpha2ratio=0.5);

      /// Get the wavelength difference for Alpha1 and Alpha2
      REAL GetXRayTubeDeltaLambda()const;
      /// Get the Kalpha2/Kalpha1 ratio
      REAL GetXRayTubeAlpha2Alpha1Ratio()const;

      /// Last time the wavelength has been changed
      const RefinableObjClock& GetClockWavelength()const ;
      /// Last time the nature (X-Rays/Neutron, number of wavelengths)radiation has been changed
      const RefinableObjClock& GetClockRadiation()const ;
      virtual void XMLOutput(ostream &os,int indent=0)const;
      virtual void XMLInput(istream &is,const XMLCrystTag &tag);
      //virtual void XMLInputOld(istream &is,const IOCrystTag &tag);
      /// Print to screen/console the charcteristics of the radiation.
      void Print()const;
      REAL GetLinearPolarRate()const;
      void SetLinearPolarRate(const REAL f);
   private:
      void InitOptions();
      /// Neutron ? X-Ray ? (Electron: unimplemented)
      RefObjOpt mRadiationType;
      /// monochromatic ? Alpha1 & Alpha2 ? Multi-Wavelength ?
      RefObjOpt mWavelengthType;
      ///Wavelength of the Experiment, in Angstroems.
      CrystVector_REAL mWavelength;
      ///Name of the X-Ray tube used, if relevant. ie "Cu", "Fe",etc...
      /// "CuA1" for Cu-alpha1, etc...
      string mXRayTubeName;
      ///Absolute difference between alpha1 and alpha2, in angstroems
      REAL mXRayTubeDeltaLambda;
      ///Ratio alpha2/alpha1 (should be 0.5)
      REAL mXRayTubeAlpha2Alpha1Ratio;
      //Clocks
         RefinableObjClock mClockWavelength;
         RefinableObjClock mClockRadiation;
      /// Linear Polarization Rate (default:0, X-Ray tube unmonochromatized)
      REAL mLinearPolarRate;
   #ifdef __WX__CRYST__
   public:
      virtual WXCrystObjBasic* WXCreate(wxWindow*);
      friend class WXRadiation;
   #endif
};

//######################################################################
/** \brief Class to compute structure factors for a set of reflections and a Crystal.
*
* This class only computes structure factor, but no intensity. i.e. it does
* not include any correction such as absorption, Lorentz or Polarization.
*
* Does this really need to be a RefinableObj ?
* \todo Optimize computation for Bijvoet/Friedel mates. To do this, generate
* an internal list of 'true independent reflections', with two entries for each,
* for both mates, and make the 'real' reflections only a reference to these reflections.
*
* \todo a \b lot of cleaning is necessary in the computing of structure
* factors, for (1) the 'preparation' part (deciding what needs to be recomputed)
* and (2) to allow anisotropic temperature factors (or other anisotropic parts)
*/
//######################################################################
class ScatteringData: virtual public RefinableObj
{
   public:
      ScatteringData();
      ScatteringData(const ScatteringData &old);
      ~ScatteringData();
      /// So-called virtual copy constructor
      virtual ScatteringData* CreateCopy()const=0;

      /** \brief input H,K,L
      *
      * \param h,k,l: REAL arrays (vectors with NbRefl elements -same size),
      *with the h, k and l coordinates of all reflections.
      */
      virtual void SetHKL( const CrystVector_REAL &h,
                           const CrystVector_REAL &k,
                           const CrystVector_REAL &l);
      /** \brief Generate a list of h,k,l to describe a full reciprocal space,
      * up to a given maximum theta value
      *
      * \param maxTheta:maximum theta value
      * \param unique: if set to true, only unique reflections will be listed.
      * Bijvoet (Friedel) pairs
      * are NOT merged, for 'anomalous' reasons, unless you have chosen to ignore the
      * imaginary part of the scattering factor.
      *
      * The multiplicity is always stored in ScatteringData::mMultiplicity.
      *
      * \warning The ScatteringData object must already have been assigned
      * a crystal object using SetCrystal(), and the experimental wavelength
      * must also have been set before calling this function.
      */
      virtual void GenHKLFullSpace2(const REAL maxsithsl,
                                   const bool unique=false);
      /** \brief Generate a list of h,k,l to describe a full reciprocal space,
      * up to a given maximum theta value
      *
      * \param maxsithsl:maximum sin(theta)/lambda=1/2d value
      * \param unique: if set to true, only unique reflections will be listed.
      * Bijvoet (Friedel) pairs
      * are NOT merged, for 'anomalous' reasons, unless you have chosen to ignore the
      * imaginary part of the scattering factor.
      *
      * The multiplicity is always stored in ScatteringData::mMultiplicity.
      *
      * \warning The ScatteringData object must already have been assigned
      * a crystal object using SetCrystal(), and the experimental wavelength
      * must also have been set before calling this function.
      *
      * \deprecated Rather use PowderPattern::GenHKLFullSpace2,
      * with a maximum sin(theta)/lambda value, which also works for dispersive experiments.
      */
      virtual void GenHKLFullSpace(const REAL maxTheta,
                                   const bool unique=false);

      ///Neutron or x-ray experiment ? Wavelength ?
      RadiationType GetRadiationType()const;
      /// Get the radiation object for this data
      virtual const Radiation& GetRadiation()const=0;

      /**Set the crystal for this experiment
      *
      */
      virtual void SetCrystal(Crystal &crystal);
      /// Const access to the data's crystal
      const Crystal& GetCrystal()const ;
      /// Access to the data's crystal
      Crystal& GetCrystal() ;
      /// Has a Crystal structure associated yet ?
      bool HasCrystal()const;

      ///Return the number of reflections in this experiment.
      long GetNbRefl() const;
      ///Return the 1D array of H coordinates for all reflections
      const CrystVector_REAL& GetH() const;
      ///Return the 1D array of K coordinates for all reflections
      const CrystVector_REAL& GetK() const;
      ///Return the 1D array of L coordinates for all reflections
      const CrystVector_REAL& GetL() const;
      /// Return the 1D array of H coordinates for all reflections, multiplied by 2*pi
      /// \internal  Should be private
      const CrystVector_REAL& GetH2Pi() const;
      ///Return the 1D array of K coordinates for all reflections, multiplied by 2*pi
      /// \internal  Should be private
      const CrystVector_REAL& GetK2Pi() const;
      ///Return the 1D array of L coordinates for all reflections, multiplied by 2*pi
      /// \internal  Should be private
      const CrystVector_REAL& GetL2Pi() const;
      ///Return the 1D array of orthonormal x coordinates for all reflections (recipr. space)
      const CrystVector_REAL& GetReflX() const;
      ///Return the 1D array of orthonormal y coordinates for all reflections (recipr. space)
      const CrystVector_REAL& GetReflY() const;
      ///Return the 1D array of orthonormal z coordinates for all reflections (recipr. space)
      const CrystVector_REAL& GetReflZ() const;

      /// Return an array with \f$ \frac{sin(\theta)}{\lambda} = \frac{1}{2d_{hkl}}\f$
      ///for all reflections
      const CrystVector_REAL& GetSinThetaOverLambda()const;
      /// Return an array with theta values for all reflections
      const CrystVector_REAL& GetTheta()const;
      /// Clock the last time the sin(theta)/lambda and theta arrays were re-computed
      const RefinableObjClock& GetClockTheta()const;

      ///  Returns the Array of calculated |F(hkl)|^2 for all reflections.
      const CrystVector_REAL& GetFhklCalcSq() const;
      std::map<RefinablePar*, CrystVector_REAL> & GetFhklCalcSq_FullDeriv(std::set<RefinablePar *> &vPar);
      /// Access to real part of F(hkl)calc
      const CrystVector_REAL& GetFhklCalcReal() const;
      /// Access to imaginary part of F(hkl)calc
      const CrystVector_REAL& GetFhklCalcImag() const;

      ///  Returns the vector of observed |F(hkl)|^2 for all reflections.
      const CrystVector_REAL& GetFhklObsSq() const;
      ///  Set the vector of observed |F(hkl)|^2 for all reflections. The supplied vector must have the same size
      /// as the mH, mK, mL vectors.
      void SetFhklObsSq(const CrystVector_REAL &obs);

      /// Scattering factors for each ScatteringPower, as vectors with NbRefl elements
      const map<const ScatteringPower*,CrystVector_REAL> &GetScatteringFactor() const;

      ///wavelength of the experiment (in Angstroems)
      CrystVector_REAL GetWavelength()const;

      /// If true, then the imaginary part of the scattering factor is ignored during
      /// Structure factor computation. (default value=false)
      ///
      /// \todo this should become useless once we take fully advantage of coupled
      /// computation of Structure Factors for Fridel/Bijvoet mates using an internal
      /// list of 'fully unique' reflections. Then only one of the mates need to be computed..
      void SetIsIgnoringImagScattFact(const bool b);
      /// If true, then the imaginary part of the scattering factor is ignored during
      /// Structure factor computation.
      bool IsIgnoringImagScattFact() const;
      // Set an option so that only low-amgle reflections (theta < angle)
      // are used. See DiffractionData::mUseOnlyLowAngleData
      //virtual void SetUseOnlyLowAngleData(const bool useOnlyLowAngle,const REAL angle)=0;
      /** \brief Print H, K, L F^2 Re(F) Im(F) theta sin(theta)/lambda for all reflections
      *
      */
      virtual void PrintFhklCalc(ostream &os=cout)const;
      /** \brief Print H, K, L sin(theta)/lambda theta F^2 Re(F) Im(F) [Re(F) Im(F)]_i,
      * where [Re(F) Im(F)]_i are the real and imaginary contribution of the different
      * scattering powers to the overall structure factor
      */
      virtual void PrintFhklCalcDetail(ostream &os=cout)const;

      virtual void BeginOptimization(const bool allowApproximations=false,
                                     const bool enableRestraints=false);
      virtual void EndOptimization();
      virtual void SetApproximationFlag(const bool allow);
      /// Set the maximum value for sin(theta)/lambda. All data (reflections,..) still
      /// exist but are ignored for all calculations.
      virtual void SetMaxSinThetaOvLambda(const REAL max);
      /// Get the maximum value for sin(theta)/lambda.
      REAL GetMaxSinThetaOvLambda()const;
      /// Recalc, and get the number of reflections which should be actually used,
      /// due to the maximuml sin(theta)/lambda value set.
      virtual long GetNbReflBelowMaxSinThetaOvLambda()const;
      /// Clock the last time the number of reflections used was changed
      const RefinableObjClock& GetClockNbReflBelowMaxSinThetaOvLambda()const;
   protected:
      /// \internal This function is called after H,K and L arrays have
      /// been initialized or modified.
      virtual void PrepareHKLarrays() ;
      /// \internal sort reflections by theta values (also get rid of [0,0,0] if present)
      /// If maxSTOL >0, then only reflections where sin(theta)/lambda<maxSTOL are kept
      /// \return an array with the subscript of the kept reflections (for inherited classes)
      virtual CrystVector_long SortReflectionBySinThetaOverLambda(const REAL maxSTOL=-1.);
      /// \internal Get rid of extinct reflections. Useful after GenHKLFullSpace().
      /// Do not use this if you have a list of observed reflections !
      ///
      /// Currently done using (brute-force) numerical evaluation. Should rather use
      /// SpaceGroup info... To do !
      ///
      /// \return an array with the subscript of the kept reflections (for inherited classes)
      CrystVector_long EliminateExtinctReflections();

      //The following functions are used during the calculation of structure factors,
         /// \internal Compute sin(theta)/lambda as well a orthonormal coordinates
         /// for all reflections. theta and tan(theta),
         /// are also re-computed, provided a wavelength has been supplied.
         virtual void CalcSinThetaLambda()const;
         /// Calculate sin(theta)/lambda for a single reflection
         REAL CalcSinThetaLambda(REAL h, REAL k, REAL l)const;
         /// Get access to the B matrix used to compute reflection positions.
         /// May be overridden by derived classes to use different lattice parameters
         /// than the Crystal's unitcell (used for multiple datasets).
         virtual const CrystMatrix_REAL& GetBMatrix()const;
         /// \internal Get scattering factors for all ScatteringPower & reflections
         void CalcScattFactor()const;
         /// \internal Compute thermic factors for all ScatteringPower & reflections
         void CalcTemperatureFactor()const;
         /// \internal get f' and f" for ScatteringPower of the crystal, at the exp. wavelength
         ///
         /// This \e could be specialized for multi-wavelength experiments...
         virtual void CalcResonantScattFactor()const;
         /**\brief Compute the overall temperature factor affecting all reflections
         */
         void CalcGlobalTemperatureFactor() const;

      /**\brief Compute the overall structure factor (real \b and imaginary part).
      *This function is \e optimized \e for \e speed (geometrical structure factors are
      *computed for all atoms and all reflections in two loops, avoiding re-calculation).
      *So use this function for repetitive calculations.
      *
      *This function recognizes the type of radiation (XRay or neutron) and
      *uses the corresponding scattering factor/length.
      *
      *  \return the result (real and imaginary part of the structure factor)
      * (mRealFhklCalc, mImagFhklCalc) are stored in ScatteringData.
      */
      void CalcStructFactor() const;
      void CalcStructFactor_FullDeriv(std::set<RefinablePar *> &vPar);
      /** \brief Compute the 'Geometrical Structure Factor' for each ScatteringPower
      * of the Crystal
      *
      */
      void CalcGeomStructFactor() const;
      void CalcGeomStructFactor_FullDeriv(std::set<RefinablePar*> &vPar);
      /** Calculate the Luzzati factor associated to each ScatteringPower and
      * each reflection, for maximum likelihood optimization.
      *
      */
      void CalcLuzzatiFactor()const;
      /** Calculate the variance associated to the calculated structure factor
      *
      */
      void CalcStructFactVariance()const;

      /// Number of H,K,L reflections
      long mNbRefl;
      /// H,K,L coordinates
      CrystVector_REAL mH, mK, mL ;
      /// H,K,L integer coordinates
      mutable CrystVector_long mIntH, mIntK, mIntL ;
      /// H,K,L coordinates, multiplied by 2PI
      mutable CrystVector_REAL mH2Pi, mK2Pi, mL2Pi ;
      /// reflection coordinates in an orthonormal base
      mutable CrystVector_REAL mX, mY, mZ ;

      ///Multiplicity for each reflections (mostly for powder diffraction)
      CrystVector_int mMultiplicity ;

      /** Expected intensity factor for all reflections.
      *
      * See SpaceGroup::GetExpectedIntensityFactor()
      */
      CrystVector_int mExpectedIntensityFactor;

      /// real &imaginary parts of F(HKL)calc
      mutable CrystVector_REAL mFhklCalcReal, mFhklCalcImag ;
      mutable std::map<RefinablePar*, CrystVector_REAL> mFhklCalcReal_FullDeriv, mFhklCalcImag_FullDeriv ;
      /// F(HKL)^2 calc for each reflection
      mutable CrystVector_REAL mFhklCalcSq ;
      mutable std::map<RefinablePar*, CrystVector_REAL> mFhklCalcSq_FullDeriv;

      /** Pointer to the crystal corresponding to this experiment.
      *
      *  This gives an access to the UB matrix for the crystal,
      * as well as to the list of Scatterer.
      */
      Crystal *mpCrystal;

      /** Global Biso, affecting the overall structure factor for all
      * reflections (but not the structure factors of individual atoms or
      * type of atomes).
      *
      */
      REAL mGlobalBiso;
      /// Global Biso factor
      mutable CrystVector_REAL  mGlobalTemperatureFactor;

      ///Use faster, but less precise, approximations for functions? (integer
      ///approximations to compute sin and cos in structure factors, and also
      ///to compute interatomic distances).
      /// This is activated by global optimization algortithms, only during the
      /// optimization.
      bool mUseFastLessPreciseFunc;

      //The Following members are only kept to avoid useless re-computation
      //during global refinements. They are used \b only by CalcStructFactor()

         ///  \f$ \frac{sin(\theta)}{\lambda} = \frac{1}{2d_{hkl}}\f$
         ///for the crystal and the reflections in ReciprSpace
         mutable CrystVector_REAL mSinThetaLambda;

         /// theta for the crystal and the HKL in ReciprSpace (in radians)
         mutable CrystVector_REAL mTheta;

         /// Anomalous X-Ray scattering term f' and f" are stored here for each ScatteringPower
         /// We store here only a value. For multi-wavelength support this should be changed
         /// to a vector... or to a matrix to take into account anisotropy of anomalous
         /// scattering...
         mutable map<const ScatteringPower*,REAL> mvFprime,mvFsecond;

         /// Thermic factors for each ScatteringPower, as vectors with NbRefl elements
         mutable map<const ScatteringPower*,CrystVector_REAL> mvTemperatureFactor;

         /// Scattering factors for each ScatteringPower, as vectors with NbRefl elements
         mutable map<const ScatteringPower*,CrystVector_REAL> mvScatteringFactor;

         /// Geometrical Structure factor for each ScatteringPower, as vectors with NbRefl elements
         mutable map<const ScatteringPower*,CrystVector_REAL> mvRealGeomSF,mvImagGeomSF;
         mutable map<RefinablePar*,map<const ScatteringPower*,CrystVector_REAL> > mvRealGeomSF_FullDeriv,mvImagGeomSF_FullDeriv;

      //Public Clocks
         /// Clock for the list of hkl
         RefinableObjClock mClockHKL;
         /// Clock for the structure factor
         mutable RefinableObjClock mClockStructFactor;
         /// Clock for the square modulus of the structure factor
         mutable RefinableObjClock mClockStructFactorSq;
      //Internal Clocks
         /// Clock the last time theta was computed
         mutable RefinableObjClock mClockTheta;
         /// Clock the last time scattering factors were computed
         mutable RefinableObjClock mClockScattFactor;
         /// Clock the last time resonant scattering factors were computed
         mutable RefinableObjClock mClockScattFactorResonant;
         /// Clock the last time the geometrical structure factors were computed
         mutable RefinableObjClock mClockGeomStructFact;
         /// Clock the last time temperature factors were computed
         mutable RefinableObjClock mClockThermicFact;

         /// last time the global Biso factor was modified
         RefinableObjClock mClockGlobalBiso;
         /// last time the global temperature factor was computed
         mutable RefinableObjClock mClockGlobalTemperatureFact;

      /** \brief Ignore imaginary part of scattering factor.
      *
      * This can be used either to speed up computation, or when f"
      * has a small effect on calculated intensities, mostly for powder
      * diffraction (GenHKLFullSpace will not generate Friedel pairs, reducing
      * the number of reflections by a factor up to 2 for some structures).
      *
      * Practically this makes f"=0 during computation. The real resonant contribution (f')
      * is not affected.
      *
      * This may be removed later on...
      */
      bool mIgnoreImagScattFact;

      // Maximum sin(theta)/lambda
         /** Maximum sin(theta)/lambda for all calculations (10 by default).
         *
         * This keeps all data in memory, but only the part which is below
         * the max is calculated.
         *
         * This affects the computing of structure factors, intensities (for single
         * crystal and powder patterns), R and Rw.
         *
         * The reflections \b must be sorted by increasing sin(theta)/lambda for
         * this to work correctly.
         */
         REAL mMaxSinThetaOvLambda;
         /// Number of reflections which are below the max. This is updated automatically
         /// from ScatteringData::mMaxSinThetaOvLambda
         mutable long mNbReflUsed;
         /// Clock recording the last time the number of reflections used has increased.
         mutable RefinableObjClock mClockNbReflUsed;

      // Maximum Likelihood
         /// The Luzzati 'D' factor for each scattering power and each reflection
         mutable map<const ScatteringPower*,CrystVector_REAL> mvLuzzatiFactor;
         /** The variance on all calculated structure factors, taking into account
         * the positionnal errors and the expected intensity factor.
         *
         * Actually this is the variance on both real and imaginary parts.
         */
         mutable CrystVector_REAL mFhklCalcVariance;
         mutable RefinableObjClock mClockLuzzatiFactor;
         mutable RefinableObjClock mClockFhklCalcVariance;
      /// Observed squared structure factors (zero-sized if none)
      CrystVector_REAL mFhklObsSq;
      /// Last time observed squared structure factors were altered
      RefinableObjClock mClockFhklObsSq;
   #ifdef __WX__CRYST__
      //to access mMaxSinThetaOvLambda
      friend class WXDiffractionSingleCrystal;
      friend class WXPowderPattern;
   #endif
};

}//namespace ObjCryst
#endif // _OBJCRYST_SCATTERINGDATA_H_
