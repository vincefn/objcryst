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

#include "CrystVector/CrystVector.h"

#include "ObjCryst/General.h"

#include "ObjCryst/SpaceGroup.h"
#include "ObjCryst/ScatteringPower.h"
#include "ObjCryst/Scatterer.h"
#include "ObjCryst/Crystal.h"

//#include <stdlib.h>
#include <string>
//#include <iomanip>
//#include <cmath>
//#include <typeinfo>
//#include <fstream>
//#include <ctime>

namespace ObjCryst
{
extern const RefParType *gpRefParTypeScattData;
extern const RefParType *gpRefParTypeScattDataScale;
extern const RefParType *gpRefParTypeScattDataProfile;
extern const RefParType *gpRefParTypeScattDataProfileType;
extern const RefParType *gpRefParTypeScattDataProfileWidth;
extern const RefParType *gpRefParTypeScattDataProfileAsym;
extern const RefParType *gpRefParTypeScattDataCorr;
extern const RefParType *gpRefParTypeScattDataCorrInt;
extern const RefParType *gpRefParTypeScattDataCorrIntAbsorp;
extern const RefParType *gpRefParTypeScattDataCorrIntPolar;
extern const RefParType *gpRefParTypeScattDataCorrIntExtinc;
extern const RefParType *gpRefParTypeScattDataCorrPos;
extern const RefParType *gpRefParTypeScattDataBackground;

extern const RefParType *gpRefParTypeRadiation;
extern const RefParType *gpRefParTypeRadiationWavelength;

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
      //Get the Wavelength type (monochromatic, Alpha1+Alpha2, ...)
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
      * \param useMultiplicity: if set to true, equivalent reflections will be removed.
      * Bijvoet (Friedel) pairs
      * are NOT merged, for 'anomalous' reasons, unless you have chosen to ignore the
      * imaginary part of the scattering factor. If true, then multiplicity is stored
      * in the mMultiplicity data member.
      *
      * \warning The ScatteringData object must already have been assigned 
      * a crystal object using SetCrystal(), and the experimental wavelength 
      * must also have been set before calling this function.
      *
      * \todo smarter generation, using spacegroup information to remove extinct reflection
      * rather than brute-force computation.
      */
      virtual void GenHKLFullSpace(const REAL maxTheta,
                                   const bool useMultiplicity=false);
      
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
      
      ///Return the number of reflections in this experiment.
      long GetNbRefl() const;
      ///Return the 1D array of H coordinates for all reflections
      const CrystVector_REAL& GetH() const;
      ///Return the 1D array of K coordinates
      const CrystVector_REAL& GetK() const;
      ///Return the 1D array of L coordinates
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
      
      /// Return an array with \f$ \frac{sin(\theta)}{\lambda} = \frac{1}{2d_{hkl}}\f$ 
      ///for all reflections
      const CrystVector_REAL& GetSinThetaOverLambda()const;
   
      ///  Returns the Array of calculated |F(hkl)|^2 for all reflections.
      const CrystVector_REAL& GetFhklCalcSq() const;
      /// Access to real part of F(hkl)calc
      const CrystVector_REAL& GetFhklCalcReal() const;
      /// Access to imaginary part of F(hkl)calc
      const CrystVector_REAL& GetFhklCalcImag() const;
      
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

      virtual void BeginOptimization(const bool allowApproximations=false,
                                     const bool enableRestraints=false);
      virtual void EndOptimization();
      /// Set the maximum value for sin(theta)/lambda. All data (reflections,..) still
      /// exist but are ignored for all calculations.
      virtual void SetMaxSinThetaOvLambda(const REAL max);
      /// Get the maximum value for sin(theta)/lambda.
      REAL GetMaxSinThetaOvLambda()const;
   protected:
      /// \internal This function is called after H,K and L arrays have 
      /// been initialized or modified.
      virtual void PrepareHKLarrays() ;
      /// \internal sort reflections by theta values (also get rid of [0,0,0] if present)
      /// If maxTheta >0, then only reflections where theta<maxTheta are kept
      /// \return an array with the subscript of the kept reflections (for inherited classes)
      virtual CrystVector_long SortReflectionByTheta(const REAL maxTheta=-1.);
      /// \internal Get rid of extinct reflections. Useful after GenHKLFullSpace().
      /// Do not use this if you have a list of observed reflections !
      ///
      /// Currently done using (brute-force) numerical evaluation. Should rather use
      /// SpaceGroup info... To do !
      ///
      /// \return an array with the subscript of the kept reflections (for inherited classes)
      CrystVector_long EliminateExtinctReflections();
      
      //The following functions are used during the calculation of structure factors,
         /// \internal Get the list of scattering components, and check what needs to
         /// be recomputed to get the new structure factors. No calculation is made in
         /// this function. Just getting prepared...
         /// \todo Clean up the code, which is a really unbelievable mess (but working!)
         ///
         /// Currently using flags to decide what should be recomputed, whereas
         /// Clocks should be used. a LOT of cleaning is necessary
         virtual void PrepareCalcStructFactor()const;
         /// \internal Compute sin(theta)/lambda. 
         /// theta and tan(theta) values are also re-computed, provided a wavelength has
         /// been supplied.
         virtual void CalcSinThetaLambda()const;
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
      * (mRealFhklCalc, mImagFhklCalc) are stored within ScatteringData.
      */
      void CalcStructFactor() const;

      /** \brief Compute the 'Geometrical Structure Factor' for each ScatteringPower
      * of the Crystal
      *
      */
      void CalcGeomStructFactor(const ScatteringComponentList &scattCompList,
                                const SpaceGroup &spg,
                                const CrystVector_long &structFactorIndex,
                                CrystVector_REAL* rsf2,
                                CrystVector_REAL* isf2,
                                bool useFastTabulatedTrigFunctions=false) const;
      
      
      /// Number of H,K,L reflections
      long mNbRefl;
      ///H,K,L coordinates
      CrystVector_REAL mH, mK, mL ;
      ///H,K,L integer coordinates
      mutable CrystVector_long mIntH, mIntK, mIntL ;
      ///H,K,L coordinates, multiplied by 2PI
      mutable CrystVector_REAL mH2Pi, mK2Pi, mL2Pi ;

      ///Multiplicity for each reflections (mostly for powder diffraction)
      CrystVector_int mMultiplicity ;
      
      /// real &imaginary parts of F(HKL)calc
      mutable CrystVector_REAL mFhklCalcReal, mFhklCalcImag ;
      ///F(HKL)^2 calc for each reflection
      mutable CrystVector_REAL mFhklCalcSq ;
      
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

         /// tan(theta) for the crystal and the HKL in ReciprSpace (for Caglioti's law)
         /// \note this should be moved to DiffractionDataPowder
         mutable CrystVector_REAL mTanTheta;

         /// Anomalous X-Ray scattering term f' and f" are stored here for each ScatteringPower 
         /// We assume yet that data is monochromatic, but this could be specialized.
         mutable CrystVector_REAL mFprime,mFsecond;

         ///Thermic factors as mNbScatteringPower vectors with NbRefl elements
         mutable CrystVector_REAL* mpTemperatureFactor;

         ///Scattering factors as mNbScatteringPower vectors with NbRefl elements
         mutable CrystVector_REAL* mpScatteringFactor;
      
         ///Geometrical Structure factor for all reflection and ScatteringPower
         mutable CrystVector_REAL* mpRealGeomSF,*mpImagGeomSF;
      
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
      
      // Info about the Scattering components and powers
         /// Pointer to the ScatteringComponentList of the crystal.
         mutable const ScatteringComponentList* mpScattCompList;
         /// Index of the storage for scattering information. This array as the same size as
         /// the total number of ScatteringPower, and for each we give the index of where
         /// their scattering, temperature and resonant factors are stored.
         /// THIS IS AWFULLY KLUDGE-ESQUE !!!
         mutable CrystVector_long mScatteringPowerIndex;
         /// This is the reverse index. KLUDGEEEEEEEE !!!
         mutable CrystVector_long mScatteringPowerIndex2;
         /// Total number os ScatteringPower used for this DiffractionData
         mutable long mNbScatteringPower;
      // :TODO: This must be replaced by Clocks...
         mutable ScatteringComponentList mLastScattCompList;
         mutable bool mAnomalousNeedRecalc;
         mutable bool mThermicNeedRecalc;
         mutable bool mScattFactNeedRecalc;
         mutable bool mGeomFhklCalcNeedRecalc;
         mutable bool mFhklCalcNeedRecalc;

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
   #ifdef __WX__CRYST__
      //to access mMaxSinThetaOvLambda
      friend class WXDiffractionSingleCrystal;
      friend class WXPowderPattern;
   #endif
};

}//namespace ObjCryst
#endif // _OBJCRYST_SCATTERINGDATA_H_
