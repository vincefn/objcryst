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
#ifndef _OBJCRYST_POWDERPATTERN_H_
#define _OBJCRYST_POWDERPATTERN_H_

#include <utility>
#include <list>
#include <string>

#include "ObjCryst/CrystVector/CrystVector.h"
#include "ObjCryst/ObjCryst/General.h"
#include "ObjCryst/ObjCryst/Crystal.h"
#include "ObjCryst/ObjCryst/ScatteringCorr.h"
#include "ObjCryst/ObjCryst/ReflectionProfile.h"
#include "ObjCryst/ObjCryst/CIF.h"
#include "ObjCryst/ObjCryst/Indexing.h"
#include "ObjCryst/ObjCryst/DiffractionDataSingleCrystal.h"

namespace ObjCryst
{

class PowderPattern;

//######################################################################
/** \brief Generic class to compute components (eg the contribution of
* a given phase, or background) of a powder pattern. This is an abstract base class.
*
* Most functions are protected, only to be accessed,
* internally or from the friend class PowderPattern.
*/
//######################################################################
class PowderPatternComponent : virtual public RefinableObj
{
   public:
      PowderPatternComponent();
      PowderPatternComponent(const PowderPatternComponent&);
      virtual ~PowderPatternComponent();
      virtual const string& GetClassName() const;

      /// Get the PowderPattern object which uses this component.
      /// This allows to know the observed powder pattern to evaluate
      /// the background.
      const PowderPattern& GetParentPowderPattern()const;
      /// Get the PowderPattern object which uses this component.
      /// This allows to know the observed powder pattern to evaluate
      /// the background.
      PowderPattern& GetParentPowderPattern();
      /// Set the PowderPattern object which uses this component.
      /// This sets all necessary pattern parameters (2theta/tof range,
      /// wavelength, radiation type...) accordingly.
      ///
      virtual void SetParentPowderPattern(PowderPattern&)=0;
      /// Get the calculated powder pattern for this component.
      /// Note that the pattern is \e not scaled.
      ///
      virtual const CrystVector_REAL& GetPowderPatternCalc()const=0;
      virtual std::map<RefinablePar*,CrystVector_REAL>& GetPowderPattern_FullDeriv(std::set<RefinablePar *> &vPar);
      /** Get the integrated values of the powder pattern
      *
      * \note: the integration intervals are those given by the parent
      *   PowderPattern, so that all PowderPatternComponent's intervals
      *   are taken into account
      *
      *   This avoids explicitely calculating the full profile powder pattern.
      */
      virtual pair<const CrystVector_REAL*,const RefinableObjClock*> GetPowderPatternIntegratedCalc()const=0;
      virtual std::map<RefinablePar*,CrystVector_REAL>& GetPowderPatternIntegrated_FullDeriv(std::set<RefinablePar *> &vPar);
      /** \brief Is this component scalable ?
      *
      * This is used by the PowderPattern class, which fits all
      * pattern components using scale factors. Some components may not
      * need to be scaled: background components, which are assumed
      * to be absolute.
      */
      bool IsScalable()const;
      /** Get the variance associated to each point of the
      * calculated powder pattern, for this component.
      *
      * \warning: this is experimental, with the aim of using Maximum Likelihood
      * to improve structure determination.
      */
      virtual const CrystVector_REAL& GetPowderPatternCalcVariance()const=0;
      /** Get the variance associated to each point of the
      * calculated powder pattern, for this component (integrated version).
      *
      * \warning: this is experimental, with the aim of using Maximum Likelihood
      * to improve structure determination.
      */
      virtual pair<const CrystVector_REAL*,const RefinableObjClock*>
         GetPowderPatternIntegratedCalcVariance() const=0;
      /// Does this component have a variance associated with each calculated
      /// point ? i.e., do we use maximum likelihood to take into account
      /// incomplete models ?
      virtual bool HasPowderPatternCalcVariance()const=0;
      /// Last time the powder pattern was calculated.
      const RefinableObjClock& GetClockPowderPatternCalc()const;
      /// Get a list of labels for the pattern (usually reflection indexes). This
      /// returns the list generated during the last computation of the powder pattern.
      const list<pair<const REAL ,const string > >& GetPatternLabelList() const;
   protected:
      /// Last time the variance on the pattern was actually calculated.
      const RefinableObjClock& GetClockPowderPatternCalcVariance()const;

      /// Calc the powder pattern. As always, recomputation is only
      /// done if necessary (ie if a parameter has changed since the last
      /// computation)
      virtual void CalcPowderPattern() const=0;
      virtual void CalcPowderPattern_FullDeriv(std::set<RefinablePar *> &vPar);
      /// Calc the integrated powder pattern. This should be optimized so that
      /// the full powder pattern is not explicitely computed.
      virtual void CalcPowderPatternIntegrated_FullDeriv(std::set<RefinablePar *> &vPar);

      /** Get the pixel positions separating the integration intervals around reflections.
      *
      * \returns: an array with the pixel positions, empty if this component
      * has no peaks. The positions should be in increasing order, but
      * could go beyond the pattern limits.
      */
      virtual const CrystVector_long& GetBraggLimits()const=0;
      /// Get last time the Bragg Limits were changed
      const RefinableObjClock& GetClockBraggLimits()const;

       /// Set the maximum value for sin(theta)/lambda. All data above still
      /// exist but are ignored for all calculations.
      virtual void SetMaxSinThetaOvLambda(const REAL max)=0;

      /// The calculated component of a powder pattern. It is mutable since it is
      /// completely defined by other parameters (eg it is not an 'independent parameter')
      mutable CrystVector_REAL mPowderPatternCalc;
      /// The calculated powder pattern, integrated.
      mutable CrystVector_REAL mPowderPatternIntegratedCalc;

      /// The variance associated to each point of the calculated powder pattern.
      mutable CrystVector_REAL mPowderPatternCalcVariance;
      /// The variance associated to each point of the calculated powder pattern, integrated
      mutable CrystVector_REAL mPowderPatternIntegratedCalcVariance;

      /// Interval limits around each reflection, for integrated R-factors
      mutable CrystVector_long mIntegratedReflLimits;

      /// \internal
      /// This will be called by the parent PowderPattern object, before
      /// calculating the first powder pattern. Or maybe it should be called
      /// automatically by the object itself...
      virtual void Prepare()=0;

      /// Scalable ? (crystal phase = scalable, background= not scalable)
      bool mIsScalable;

      //Clocks
         /// When was the powder pattern last computed ?
         mutable RefinableObjClock mClockPowderPatternCalc;
         /// When was the 'integrated' powder pattern last computed ?
         mutable RefinableObjClock mClockPowderPatternIntegratedCalc;
         /// When was the powder pattern variance last computed ?
         mutable RefinableObjClock mClockPowderPatternVarianceCalc;
         /// When was the 'integrated' powder pattern variance last computed ?
         mutable RefinableObjClock mClockPowderPatternIntegratedVarianceCalc;

      /// The PowderPattern object in which this component is included
      PowderPattern *mpParentPowderPattern;
      /// Get last time the Bragg Limits were changed
      mutable RefinableObjClock mClockBraggLimits;

      /// The labels associated to different points of the pattern
      mutable list<pair<const REAL ,const string > > mvLabel;

      mutable std::map<RefinablePar*,CrystVector_REAL> mPowderPattern_FullDeriv;
      mutable std::map<RefinablePar*,CrystVector_REAL> mPowderPatternIntegrated_FullDeriv;
      //Eventually this should be removed (?)
      friend class PowderPattern;
};
/// Global registry for all PowderPatternComponent objects
extern ObjRegistry<PowderPatternComponent> gPowderPatternComponentRegistry;

//######################################################################
/** \brief Phase to compute a background contribution to a powder
* pattern using an interpolation. Currently only linear interpolation is
* available. (in the works: cubic spline interpolation background)
*/
//######################################################################
class PowderPatternBackground : public PowderPatternComponent
{
   public:
      PowderPatternBackground();
      PowderPatternBackground(const PowderPatternBackground&);
      virtual ~PowderPatternBackground();
      virtual const string& GetClassName() const;

      virtual void SetParentPowderPattern(PowderPattern&);
      virtual const CrystVector_REAL& GetPowderPatternCalc()const;
      virtual pair<const CrystVector_REAL*,const RefinableObjClock*>
         GetPowderPatternIntegratedCalc()const;
      /// Import background points from a file (with two columns 2theta (or tof), intensity)
      void ImportUserBackground(const string &filename);
      void SetInterpPoints(const CrystVector_REAL tth, const CrystVector_REAL backgd);
      const std::pair<const CrystVector_REAL*,const CrystVector_REAL*> GetInterpPoints()const;
      virtual void XMLOutput(ostream &os,int indent=0)const;
      virtual void XMLInput(istream &is,const XMLCrystTag &tag);
      //virtual void XMLInputOld(istream &is,const IOCrystTag &tag);
      virtual void GetGeneGroup(const RefinableObj &obj,
                                CrystVector_uint & groupIndex,
                                unsigned int &firstGroup) const;
      virtual void BeginOptimization(const bool allowApproximations=false,
                                     const bool enableRestraints=false);
      virtual const CrystVector_REAL& GetPowderPatternCalcVariance()const;
      virtual pair<const CrystVector_REAL*,const RefinableObjClock*>
         GetPowderPatternIntegratedCalcVariance() const;
      virtual bool HasPowderPatternCalcVariance()const;
      virtual void TagNewBestConfig()const;
      /** Optimize the background using a Bayesian approach. The background parameters
      * must be un-fixed before.
      *
      * The minimization will a maximum of 50 Simplex runs (see the SimplexObj documentation),
      * each with 200 cycles.
      *
      * See the class documentation for PowderPatternBackgroundBayesianMinimiser.
      */
      void OptimizeBayesianBackground();
      /** Fix parameters corresponding to points of the pattern that are not actually calculated.
      * This is necessary for modelling using splines, to avoid divergence of interpolation
      * points during least squares optimization.
      *
      * \param obj: the object in which are parameters to be fixed. Normally this will be
      * the PowderPatternBackground object itself, but it can also be the parameter list
      * copied such as in a LSQNumObj.
      */
      void FixParametersBeyondMaxresolution(RefinableObj &obj);
   protected:
      virtual void CalcPowderPattern() const;
      virtual void CalcPowderPattern_FullDeriv(std::set<RefinablePar *> &vPar);
      virtual void CalcPowderPatternIntegrated() const;
      virtual void CalcPowderPatternIntegrated_FullDeriv(std::set<RefinablePar *> &vPar);
      virtual void Prepare();
      virtual const CrystVector_long& GetBraggLimits()const;
      virtual void SetMaxSinThetaOvLambda(const REAL max);
      void InitRefParList();
      void InitOptions();
      void InitSpline()const;
      /// Number of fitting points for background
      int mBackgroundNbPoint;
      /// Vector of 2theta values for the fitting points of the background
      CrystVector_REAL mBackgroundInterpPointX;
      /// Values of background at interpolating points
      CrystVector_REAL mBackgroundInterpPointIntensity;
      /// Subscript of the points, sorted the correct order,
      ///taking into account the type of radiation (monochromatic/TOF).
      mutable CrystVector_long mPointOrder;
      /// Vector of pixel values between each interval, for faster CubicSpline calculations.
      /// Mutable since it copies information from mBackgroundInterpPointX.
      mutable CrystVector_REAL mvSplinePixel;
      /// Spline used for interpolation.
      /// Mutable since it copies information from mBackgroundInterpPointX
      ///and mBackgroundInterpPointIntensity.
      mutable CubicSpline mvSpline;
      // Clocks
         /// Modification of the interpolated points
         RefinableObjClock mClockBackgroundPoint;
         /// Initialization of the spline
         mutable RefinableObjClock mClockSpline;

      /** Maximum sin(theta)/lambda for all calculations (10 by default).
      *
      * This keeps all data in memory, but only the part which is below
      * the max is calculated.
      */
      REAL mMaxSinThetaOvLambda;

      /// Constant error (sigma) on the calculated pattern, due to an incomplete
      /// model
      REAL mModelVariance;

      /// Type of interpolation performed: linear or cubic spline
      RefObjOpt mInterpolationModel;
      //To be removed
      friend class PowderPattern;
   #ifdef __WX__CRYST__
   public:
      virtual WXCrystObjBasic* WXCreate(wxWindow*);
      friend class WXPowderPatternBackground;
   #endif
};

//######################################################################
/** \brief Class to compute the contribution to a powder pattern from
* a crystalline phase.
*
*/
//######################################################################
class PowderPatternDiffraction : virtual public PowderPatternComponent,public ScatteringData
{
   public:
      PowderPatternDiffraction();
      PowderPatternDiffraction(const PowderPatternDiffraction&);
      virtual ~PowderPatternDiffraction();
      virtual PowderPatternDiffraction* CreateCopy()const;
      virtual const string& GetClassName() const;

      virtual void SetParentPowderPattern(PowderPattern&);
      virtual const CrystVector_REAL& GetPowderPatternCalc()const;
      virtual pair<const CrystVector_REAL*,const RefinableObjClock*>
         GetPowderPatternIntegratedCalc()const;

      /** Set reflection profile parameters
      *
      * :TODO: assymmetric profiles
      * \param fwhmCagliotiW,fwhmCagliotiU,fwhmCagliotiV : these are the U,V and W
      * parameters in the Caglioti's law :
      * \f$ fwhm^2= U \tan^2(\theta) + V \tan(\theta) +W \f$
      * if only W is given, the width is constant
      * \param eta0,eta1: these are the mixing parameters in the case of a
      * pseudo-Voigt function.
      */
      void SetReflectionProfilePar(const ReflectionProfileType prof,
                                   const REAL fwhmCagliotiW,
                                   const REAL fwhmCagliotiU=0,
                                   const REAL fwhmCagliotiV=0,
                                   const REAL eta0=0.5,
                                   const REAL eta1=0.);
      /** Assign a new profile
      *
      */
      void SetProfile(ReflectionProfile *prof);
      /// Get reflection profile
      const ReflectionProfile& GetProfile()const;
      /// Get reflection profile
      ReflectionProfile& GetProfile();
      virtual void GenHKLFullSpace();
      virtual void XMLOutput(ostream &os,int indent=0)const;
      virtual void XMLInput(istream &is,const XMLCrystTag &tag);
      //virtual void XMLInputOld(istream &is,const IOCrystTag &tag);
      virtual void GetGeneGroup(const RefinableObj &obj,
                                CrystVector_uint & groupIndex,
                                unsigned int &firstGroup) const;
      virtual void BeginOptimization(const bool allowApproximations=false,
                                     const bool enableRestraints=false);
      virtual void EndOptimization();
      virtual void SetApproximationFlag(const bool allow);
      virtual const Radiation& GetRadiation()const;
      virtual const CrystVector_REAL& GetPowderPatternCalcVariance()const;
      virtual pair<const CrystVector_REAL*,const RefinableObjClock*>
         GetPowderPatternIntegratedCalcVariance() const;
      virtual bool HasPowderPatternCalcVariance()const;
      virtual void SetCrystal(Crystal &crystal);
      /** Prepare intensity extraction (Le Bail or Pawley)
      *
      *
      *\param extract: if true, begin extraction mode, else enable structure factor calculations
      *\param init: if true and extract=true, intensities are set to 100. Otherwise if extract==true
      * and init=false, the program will try to re-use existing extracted data (in mpLeBailData),
      * but only if the list of HKL is unchanged. Otherwise initilization to 100 will be forced.
      */
      void SetExtractionMode(const bool extract=true,const bool init=false);
      /// Return true if in extraction mode, i.e. using extracted intensities instead of computed structure factors.
      bool GetExtractionMode()const;
      /** Extract intensities using Le Bail method
      *
      *\param nbcycle: number of cycles
      */
      void ExtractLeBail(unsigned int nbcycle=1);
      /// Recalc, and get the number of reflections which should be actually used,
      /// due to the maximuml sin(theta)/lambda value set.
      virtual long GetNbReflBelowMaxSinThetaOvLambda()const;
      /// Change one parameter in mFrozenLatticePar. This triggers a call to CalcLocalBMatrix() if the parameter has changed
      void SetFrozenLatticePar(const unsigned int i, REAL v);
      /// Access to one parameter in mFrozenLatticePar
      REAL GetFrozenLatticePar(const unsigned int i) const;
      /// Set the use local cell parameters ? (see mFrozenLatticePar)
      /// If this changes mUseLocalLatticePar from false to true, this triggers a copy
      /// of the Crystal's lattice parameters into mFrozenLatticePar
      void FreezeLatticePar(const bool use);
      /// Do we use local cell parameters ? (see mFrozenLatticePar)
      bool FreezeLatticePar() const;
   protected:
      virtual void CalcPowderPattern() const;
      virtual void CalcPowderPattern_FullDeriv(std::set<RefinablePar *> &vPar);
      virtual void CalcPowderPatternIntegrated() const;
      virtual void CalcPowderPatternIntegrated_FullDeriv(std::set<RefinablePar *> &vPar);

      /// \internal Calc reflection profiles for ALL reflections (powder diffraction)
      void CalcPowderReflProfile()const;
      /// \internal Calc derivatives of reflection profiles for all used reflections,
      /// for a given list of refinable parameters
      void CalcPowderReflProfile_FullDeriv(std::set<RefinablePar *> &vPar);
      /// \internal Calc Lorentz-Polarisation-Aperture correction
      void CalcIntensityCorr()const;
      /// \internal Compute the intensity for all reflections (taking into account
      /// corrections, but not the multiplicity)
      virtual void CalcIhkl() const;
      virtual void CalcIhkl_FullDeriv(std::set<RefinablePar*> &vPar);
      virtual void Prepare();
      virtual void InitOptions();
      virtual const CrystVector_long& GetBraggLimits()const;
      virtual void SetMaxSinThetaOvLambda(const REAL max);
      /// This can use either locally stored lattice parameters from mLocalLatticePar,
      /// or the Crystal's, depending on mUseLocalLatticePar.
      virtual const CrystMatrix_REAL& GetBMatrix()const;
      /// Calculate the local BMatrix, used if mFreezeLatticePar is true.
      void CalcFrozenBMatrix()const;
      void PrepareIntegratedProfile()const;
      //Clocks
         /// Last time the reflection parameters were changed
         RefinableObjClock mClockProfilePar;
         /// Last time the Lorentz-Polarization and slit parameters were changed
         RefinableObjClock mClockLorentzPolarSlitCorrPar;
      //Clocks (internal, mutable)
         /// Last time the Lorentz-Polar-Slit correction was computed
         mutable RefinableObjClock mClockIntensityCorr;
         /// Last time the reflection profiles were computed
         mutable RefinableObjClock mClockProfileCalc;
         /// Last time intensities were computed
         mutable RefinableObjClock mClockIhklCalc;
      /// Profile
         ReflectionProfile *mpReflectionProfile;
      // Corrections
         /** \brief Calculated corrections for all reflections. Calc F^2 must be multiplied
         *by this factor to yield intensities.
         *
         * Thus we have : \f$ I_{hkl} = L \times P \times SlitAp \times F_{hkl}^2 \f$
         *
         *with \f$ L = \frac{1}{\sin(2\theta)} \f$ (Lorentz factor).
         *\f$ P = \frac{1}{1+A}\left(1+A\cos^2(2\theta)\right) \f$ (Polarization factor), with
         * \f$ A = \frac{1-f}{1+f} \f$, where f is the polarization rate of the incident
         *beam in the plane which (i) includes the incident beam, and (ii) is perpendicular to
         *the diffracting plane. For an X-Ray Tube without monochromator, A=1, and
         *if there is a monochromator : \f$ A = \cos^2(2\theta_{mono}) \f$
         *The factor \f$ SlitAp = \frac{1}{\sin(\theta)} \f$ takes into account the
         *fraction of the diffracted cone which falls in the detector slit.
         *
         * If there is prefereed orientation, this also holds the associated correction.
         *
         * \todo: store all corrections in a registry, so that other corrections
         * can more easily be added (? Maybe not that useful, especially since these
         * correction do not need to be displayed to the user ?).
         */
         mutable CrystVector_REAL mIntensityCorr;
         /// Lorentz correction
         LorentzCorr mCorrLorentz;
         /// Polarization correction
         PolarizationCorr mCorrPolar;
         /// Slit aperture correction
         PowderSlitApertureCorr mCorrSlitAperture;
         /// Preferred orientation (texture) correction following the March-Dollase model
         TextureMarchDollase mCorrTextureMarchDollase;
         /// Preferred orientation (texture) correction following the Ellipsoidal function
         TextureEllipsoid mCorrTextureEllipsoid;
         /// Time-Of-Flight intensity correction
         TOFCorr mCorrTOF;

      /// Computed intensities for all reflections
         mutable CrystVector_REAL mIhklCalc;
         mutable std::map<RefinablePar*,CrystVector_REAL> mIhkl_FullDeriv;
      /// Variance on computed intensities for all reflections
         mutable CrystVector_REAL mIhklCalcVariance;

      // Saved arrays to speed-up computations
         /// Profile of a single reflection
         struct ReflProfile
         {
            /// First point of the pattern for which the profile is calculated
            long first;
            /// Last point of the pattern for which the profile is calculated
            long last;
            /// The profile
            CrystVector_REAL profile;
         };
         ///Reflection profiles for ALL reflections during the last powder pattern generation
         mutable vector<ReflProfile> mvReflProfile;
         /// Derivatives of reflection profiles versus a list of parameters. This will be limited
         /// to the reflections actually used. First and last point of each profile
         /// are the same as in mvReflProfile.
         mutable std::map<RefinablePar*,vector<CrystVector_REAL> > mvReflProfile_FullDeriv;

      // When using integrated profiles
         /** For each reflection, store the integrated value of the normalized
         * profile over all integration intervals.
         *
         * The first field is the first integration interval to which the reflection
         * contributes, and the second field is a vector with all the integrated
         * values for the intervals, listed in ascending 2theta(tof) order.
         */
         mutable vector< pair<unsigned long, CrystVector_REAL> > mIntegratedProfileFactor;
         /// Last time the integrated values of normalized profiles was calculated.
         mutable RefinableObjClock mClockIntegratedProfileFactor;
      /// Extraction mode (Le Bail, Pawley)
      bool mExtractionMode;
      /// Single crystal data extracted from the powder pattern.
      DiffractionDataSingleCrystal *mpLeBailData;
      /// a,b and c in Angstroems, angles (stored) in radians
      /// This is used to override lattice parameter from the Crystal structure,
      /// e.g. for multiple datasets collected at different temperatures
      /// Ignored unless mFreezeLatticePar is true
      mutable CrystVector_REAL mFrozenLatticePar;
      /// If true, use local cell parameters from mFrozenLatticePar rather than the Crystal
      bool mFreezeLatticePar;
      /// Local B Matrix, used if mFreezeLatticePar is true
      mutable CrystMatrix_REAL mFrozenBMatrix;
   #ifdef __WX__CRYST__
   public:
      virtual WXCrystObjBasic* WXCreate(wxWindow*);
      friend class WXPowderPatternDiffraction;
   #endif
};


//######################################################################
/** \brief Powder pattern class, with an observed pattern and several
* calculated components to modelize the pattern.
*
* This can also be used for simulation, using a fake Iobs. Supports
* multiple phases.
*
*/
//######################################################################
class PowderPattern : public RefinableObj
{
    public:
      PowderPattern();
      PowderPattern(const PowderPattern&);
      ~PowderPattern();
      virtual const string& GetClassName() const;
      /** Add a component (phase, backround) to this pattern.
      *
      * It must have been allocated in the heap. The pattern parameters (2theta min,
      * step, nbpoints, wavelength, radiation type) of the component
      * are automatically changed to that of the PowderPattern object.
      */
      void AddPowderPatternComponent(PowderPatternComponent &);
      /// Number of components
      unsigned int GetNbPowderPatternComponent()const;
      /// Access to a component of the powder pattern
      const PowderPatternComponent& GetPowderPatternComponent(const string &name)const;
      /// Access to a component of the powder pattern
      const PowderPatternComponent& GetPowderPatternComponent(const int)const;
      /// Access to a component of the powder pattern
      PowderPatternComponent& GetPowderPatternComponent(const string &name);
      /// Access to a component of the powder pattern
      PowderPatternComponent& GetPowderPatternComponent(const int);
      /// Access to the scale factor of components (will be 1 for background components)
      REAL GetScaleFactor(const int i)const;
      /// Access to the scale factor of components (will be 1 for background components)
      REAL GetScaleFactor(const PowderPatternComponent &comp)const;
      /// Access to the scale factor of components (will be 1 for background components)
      void SetScaleFactor(const int i, REAL s);
      /// Access to the scale factor of components (will be 1 for background components)
      void SetScaleFactor(const PowderPatternComponent &comp, REAL s);

      // Pattern parameters (2theta range, wavelength, radiation)
         /** \briefSet the powder pattern angular range & resolution parameter.
         * this will affect all components (phases) of the pattern.
         *
         *   Use this with caution, as the number of points must be correct with
         * respect to the observed data (Iobs).
         *
         * \param min: min 2theta (in radians) or time-of-flight
         *(in microseconds) value,
         * \param step: step (assumed constant) in 2theta or time-of-flight
         * (in microseconds).
         * \param nbPoints: number of points in the pattern.
         *
         * \warning : use only this for constant-step patterns. Otherwise, use
         * PowderPattern::SetPowderPatternX()
         */
         void SetPowderPatternPar(const REAL min,
                                  const REAL step,
                                  unsigned long nbPoint);
         /** Set the x coordinate of the powder pattern : either the
         * 2theta or time-of-flight values for each recorded point. The
         * step need not be constant, but the variation must be strictly
         * monotonous.
         *
         * 2theta must be in radians and time-of-flight in microseconds
         */
         void SetPowderPatternX(const CrystVector_REAL &x);
         ///Number of points ?
         unsigned long GetNbPoint()const;
         ///Number of points actually calculated (below the chosen max(sin(theta)/lambda)) ?
         unsigned long GetNbPointUsed()const;
         /// Clock corresponding to the last time the number of points used was changed
         const RefinableObjClock& GetClockNbPointUsed()const;

         /// Set the radiation
         void SetRadiation(const Radiation &radiation);
         ///Neutron or x-ray experiment ?
         const Radiation& GetRadiation()const;
         ///Neutron or x-ray experiment ?
         Radiation& GetRadiation();

         /// Set the radiation type
         void SetRadiationType(const RadiationType radiation);
         ///Neutron or x-ray experiment ?
         RadiationType GetRadiationType()const;
         /** Set the wavelength of the experiment (in Angstroems).
         *
         * \note: this is only useful for a monochromatic (X-Ray or neutron)
         * powder pattern.
         */
         void SetWavelength(const REAL lambda);

         /** \brief Set the wavelength of the experiment to that of an X-Ray tube.
         *
         *\param XRayTubeElementName : name of the anticathode element name. Known
         *ones are Cr, Fe, Cu, Mo, Ag.
         *\param alpha2Alpha2ratio: Kalpha2/Kalpha1 ratio (0.5 by default)
         *
         *Alpha1 and alpha2 wavelength are taken
         *from R. Grosse-Kunstleve package, and the average wavelength is calculated
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
         void SetWavelength(const string &XRayTubeElementName,const REAL alpha12ratio=0.5);

         /// Set the energy of the experiment [in keV, lambda(A)=12398/E(keV)].
         void SetEnergy(const REAL energy);
         /// wavelength of the experiment (in Angstroems)
         REAL GetWavelength()const;

      //Access to pattern data
         /// Get the calculated powder pattern
         const CrystVector_REAL& GetPowderPatternCalc()const;
         std::map<RefinablePar*,CrystVector_REAL>& GetPowderPattern_FullDeriv(std::set<RefinablePar *> &vPar);
         /// Get the observed powder pattern
         const CrystVector_REAL& GetPowderPatternObs()const;
         /// Get the sigma for each point of the observed powder pattern
         const CrystVector_REAL& GetPowderPatternObsSigma()const;
         /// Get the variance (obs+model) for each point of the powder pattern
         const CrystVector_REAL& GetPowderPatternVariance()const;
         /// Get the weight for each point of the powder pattern
         const CrystVector_REAL& GetPowderPatternWeight()const;
         /// Get the Minimum 2theta
         REAL GetPowderPatternXMin()const;
         /** Get the average step in 2theta
         *
         *  \warning : this will only return (2ThetaMax-2ThetaMin)/(nbPoints-1),
         * so this is the 2theta step only if the step is fixed.
         *
         * \deprecated
         */
         REAL GetPowderPatternXStep()const;
         /// Get the maximum 2theta
         REAL GetPowderPatternXMax()const;
         /// Get the vector of X (2theta or time-of-flight) coordinates
         const CrystVector_REAL& GetPowderPatternX()const;
         /** Get the powder pattern cumulative Chi^2. Depending on the chosen option,
         *it will be calculated in an integrated manner or not.
         *
         * The vector is recomputed on every call, so this is \e slow.
         */
         const CrystVector_REAL& GetChi2Cumul()const;

      // Clocks
         /// Last time the pattern was calculated
         const RefinableObjClock& GetClockPowderPatternCalc()const;
         /// When were the pattern parameters (2theta range, step) changed ?
         const RefinableObjClock& GetClockPowderPatternPar()const;
         /// When were the radiation parameter (radiation type, wavelength) changed ?
         const RefinableObjClock& GetClockPowderPatternRadiation()const;
         /// When were the parameters for 2theta/TOF correction (zero, transparency,
         /// displacement) last changed ?
         const RefinableObjClock& GetClockPowderPatternXCorr()const;

      // Corrections to the x (2theta, tof) coordinate
         ///Change Zero in x (2theta,tof)
         void SetXZero(const REAL newZero);
         /// Change displacement correction
         /// \f$ (2\theta)_{obs} = (2\theta)_{real} + a\cos(\theta) \f$
         void Set2ThetaDisplacement(const REAL displacement);
         ///Change transparency correction
         /// \f$ (2\theta)_{obs} = (2\theta)_{real} + b\sin(2\theta) \f$
         void Set2ThetaTransparency(const REAL transparency);

      // Import & export powder pattern
         /** \brief Import fullprof-style diffraction data.
         *\param fullprofFileName: filename
         */
         void ImportPowderPatternFullprof(const string &fullprofFileName);
         /** \brief Import powder pattern, format DMC from PSI
         */
         void ImportPowderPatternPSI_DMC(const string &filename);
         /** \brief Import powder pattern, format from ILL D1A/D2B
         * (format without counter info)
         */
         void ImportPowderPatternILL_D1A5(const string &filename);
         /** \brief Import *.xdd diffraction data (Topas,...).
         *\param fileName: filename
         */
         void ImportPowderPatternXdd(const string &fileName);
         /** \brief Import *.cpi Sietronics diffraction data
         *\param fileName: filename
         */
         void ImportPowderPatternSietronicsCPI(const string &fileName);
         /** \brief Import file with 3 columns 2Theta Iobs Sigma.
         *\param fileName: the filename (surprise!)
         *\param nbSkip: the number of lines to skip at the beginning of the file (default=0)
         */
         void ImportPowderPattern2ThetaObsSigma(const string &fileName,const int nbSkip=0);
         /** \brief Import diffraction data from a file, with the first line
         * has 2ThetaMin, step, 2thetaMax, and the following lines alternate
         * 10 Iobs and 10 sigma. Ends with null entries (to fill last Iobs line
         * to reach last sigme line).

         * That's fullprof format #4.
         *\param fileName: filename
         */
         void ImportPowderPatternFullprof4(const string &fileName);
         /** \brief diffraction data in a multi-detector format (fullprof format #6).
         *
         * First line is text. Third entry of second line is the 2theta step. Third
         *line has the 2thetamin, fourth line has monitors and temperatures.
         * Then each line has ten pairs (I2,I8)of NbCounters,intensity. Ends with
         * negative entries.
         *\param fileName: filename
         */
         void ImportPowderPatternMultiDetectorLLBG42(const string &fileName);
         /** \brief Import file with 2 columns 2Theta Iobs.
         *
         *\param fileName: the filename (surprise!)
         *\param nbSkip: the number of lines to skip at the beginning of the file (default=0)
         */
         void ImportPowderPattern2ThetaObs(const string &fileName,const int nbSkip=0);
         /** \brief Import TOF file (ISIS type, 3 columns t, Iobs, sigma(Iobs))
         *\param fileName: the filename
         */
         void ImportPowderPatternTOF_ISIS_XYSigma(const string &fileName);
         /** Import GSAS standard powder pattern data (see GSAS manual).
         * \warning : partial support (only CONST-constant wavelength- data so far)
         */
         void ImportPowderPatternGSAS(const string &fileName);
         /** Import CIF powder pattern data.
         */
         void ImportPowderPatternCIF(const CIF &cif);
         /** \brief Set observed powder pattern from vector array.
         *
         * Note: powder pattern parameters must have been set before calling this function,
         * for example by calling DiffractionDataPowder::InitPowderPatternPar().
         */
         void SetPowderPatternObs(const CrystVector_REAL& obs);

         ///Save powder pattern to one file, text format, 3 columns theta Iobs Icalc.
         ///If Iobs is missing, the column is omitted.
         ///
         /// \todo export in other formats (.prf,...), with a list of reflection
         /// position for all phases...
         void SavePowderPattern(const string &filename="powderPattern.out") const;
         /// Print to thee screen/console the observed and calculated pattern (long,
         /// mostly useful for debugging)
         void PrintObsCalcData(ostream&os=cout)const;

      // Statistics..
         /** \brief  Unweighted R-factor
         *
         * \return \f$ R= \sqrt {\frac{\sum_i \left( I_i^{obs}-I_i^{calc} \right)^2}
         * {\sum_i (I_i^{obs})^2} }\f$
         */
         REAL GetR()const ;
         REAL GetIntegratedR()const ;
         /** Get the weighted R-factor
         * \return \f$ R_{w}= \sqrt {\frac{\sum_i w_i\left( I_i^{obs}-I_i^{calc} \right)^2}
         * {\sum_i w_i (I_i^{obs})^2} }\f$
         */
         REAL GetRw()const;
         REAL GetIntegratedRw()const;
         /** \brief  Return conventionnal Chi^2
         *
         * \return \f$ \chi^2 = \sum_i w_i \left(I_i^{obs}-I_i^{calc} \right)^2
         *  \f$
         */
         REAL GetChi2()const;
         /** \brief  Return integrated Chi^2
         *
         */
         REAL GetIntegratedChi2()const;
         /** Return the conventionnal or integrated Chi^2, depending on the option.
         *
         */
         REAL GetChi2_Option()const;
         /// Fit the scale(s) factor of each component to minimize R
         void FitScaleFactorForR()const;
         void FitScaleFactorForIntegratedR()const;
         /// Fit the scale(s) factor of each component to minimize Rw
         void FitScaleFactorForRw()const;
         void FitScaleFactorForIntegratedRw()const;
         /// Set sigma=sqrt(Iobs)
         void SetSigmaToSqrtIobs();
         /// Set w = 1/sigma^2.
         ///
         /// To filter too small or null intensities :If sigma< minRelatSigma* max(sigma),
         /// then w=1/(minRelatSigma* max(sigma))^2
         void SetWeightToInvSigmaSq(const REAL minRelatSigma=1e-3);
         /// Set w = 1
         void SetWeightToUnit();
         /// Set w = 1/(a+ Iobs + b*Iobs^2+c*Iobs^3)
         ///
         /// To filter too small or null intensities:
         /// if Iobs < [minRelatIobs * max(Iobs)], then use Iobs=minRelatIobs * max(Iobs)
         /// to compute the weight.
         ///
         /// Typical values: a=2*min(Iobs) b=2/max(Iobs) c=0
         void SetWeightPolynomial(const REAL a, const REAL b, const REAL c,
                                     const REAL minRelatIobs=1e-3);

         /// Add an Exclusion region, in 2theta, which will be ignored when computing R's
         /// XMLInput values must be, as always, in radians. Does not work yet with
         /// integrated R factors.
         /// Note that the pattern is still computed in these regions. They are only ignored
         /// by statistics functions (R, Rws).
         void AddExcludedRegion(const REAL min2Theta,const REAL max2theta);

      virtual void BeginOptimization(const bool allowApproximations=false,
                                     const bool enableRestraints=false);
      //virtual void SetApproximationFlag(const bool allow);
      virtual void GlobalOptRandomMove(const REAL mutationAmplitude,
                                       const RefParType *type=gpRefParTypeObjCryst);
      virtual REAL GetLogLikelihood()const;
      //LSQ functions
         virtual unsigned int GetNbLSQFunction()const;
         virtual const CrystVector_REAL& GetLSQCalc(const unsigned int) const;
         virtual const CrystVector_REAL& GetLSQObs(const unsigned int) const;
         virtual const CrystVector_REAL& GetLSQWeight(const unsigned int) const;
         virtual std::map<RefinablePar*, CrystVector_REAL>& GetLSQ_FullDeriv(const unsigned int,std::set<RefinablePar *> &vPar);
      // I/O
         virtual void XMLOutput(ostream &os,int indent=0)const;
         virtual void XMLInput(istream &is,const XMLCrystTag &tag);
         //virtual void XMLInputOld(istream &is,const IOCrystTag &tag);
         void Prepare();
      virtual void GetGeneGroup(const RefinableObj &obj,
                                CrystVector_uint & groupIndex,
                                unsigned int &firstGroup) const;
      /// Set the maximum value for sin(theta)/lambda. All data (reflections,..) still
      /// exist but are ignored for all calculations.
      virtual void SetMaxSinThetaOvLambda(const REAL max);
      /// Get the maximum value for sin(theta)/lambda.
      REAL GetMaxSinThetaOvLambda()const;

      // For integrated pattern calculations
         /// Get the list of first pixels for the integration intervals
         const CrystVector_long& GetIntegratedProfileMin()const;
         /// Get the list of last pixels for the integration intervals
         const CrystVector_long& GetIntegratedProfileMax()const;
         /// When were the integration intervals last changed ?
         const RefinableObjClock& GetIntegratedProfileLimitsClock()const;
      /// Get the experimental x (2theta, tof) from the theoretical value, taking
      /// into account all corrections (zero, transparency,..).
      /// \internal
      /// \param ttheta: the theoretical x (2theta, tof) value.
      /// \return the x (2theta, tof) value as it appears on the pattern.
      REAL X2XCorr(const REAL x)const;
      /// Get the pixel number on the experimental pattern, from the
      /// theoretical (uncorrected) x coordinate, taking into account all corrections.
      /// (zero, transparency,..).
      /// \internal
      /// \param x: the theoretical x (2theta, tof) value.
      /// \return the x (2theta, tof) value as it appears on the pattern.
      ///
      /// \warning: this can be real slow, especially for non-fixed steps.
      ///
      /// \warning: this returns the exact pixel coordinate, as a floating-point
      /// value, and \e not the closest pixel coordinate.
      REAL X2PixelCorr(const REAL x)const;
      /// Get the pixel number on the experimental pattern, corresponding
      /// to a given (experimental) x coordinate
      /// \param x: the x (2theta, tof) value.
      /// \return the x (2theta, tof) value as it appears on the pattern.
      ///
      /// \warning: this can be real slow, especially for non-fixed steps.
      ///
      /// \warning: this returns the exact pixel coordinate, as a floating-point
      /// value, and \e not the closest pixel coordinate.
      REAL X2Pixel(const REAL x)const;

      /// Convert sin(theta)/lambda to X (i.e. either to 2theta or to TOF),
      /// depending on the type of radiation.
      ///
      /// This does not take into account any zero/transparency, etc... correction
         REAL STOL2X(const REAL stol)const;
      /// Convert X (either 2theta or TOF) to sin(theta)/lambda,
      /// depending on the type of radiation.
      ///
      /// This does not take into account any zero/transparency, etc... correction
         REAL X2STOL(const REAL x)const;
      /// Convert sin(theta)/lambda to pixel,
      /// depending on the type of radiation.
      ///
      /// This does not take into account any zero/transparency, etc... correction
         REAL STOL2Pixel(const REAL stol)const;
      /// Find peaks in the pattern
      PeakList FindPeaks(const float dmin=2.0,const float maxratio=0.01,const unsigned int maxpeak=100);
      /// Access the scale factors (see PowderPattern::mScaleFactor)
      const CrystVector_REAL &GetScaleFactor() const;
      /// Access the scale factors (see PowderPattern::mScaleFactor)
      CrystVector_REAL &GetScaleFactor();
      /** Export powder pattern & crystal structure in Fullprof format.
      *
      * This will create two files - the .pcr file (including the crystal structure
      * and all pattern parameters), and the .dat file with the powder pattern,
      * written using the "Ins=10" file format.
      * \param prefix: the prefix used to output the two files, 'prefix'.pcr and 'prefix'.dat
      *
      * \note: in development. Only supports constant wavelength neutron & X-ray patterns.
      */
      void ExportFullprof(const std::string &prefix)const;
   protected:
      /// Calc the powder pattern
      void CalcPowderPattern() const;
      void CalcPowderPattern_FullDeriv(std::set<RefinablePar *> &vPar);
      /// Calc the integrated powder pattern
      void CalcPowderPatternIntegrated() const;
      void CalcPowderPatternIntegrated_FullDeriv(std::set<RefinablePar *> &vPar);
      /// Init parameters and options
      virtual void Init();
      /// Prepare  the calculation of the integrated R-factors
      void PrepareIntegratedRfactor()const;
      /// Calculate the number of points of the pattern actually used, from the maximum
      /// value of sin(theta)/lambda
      void CalcNbPointUsed()const;
      /// Initialize options
      virtual void InitOptions();

      /// The calculated powder pattern. It is mutable since it is
      /// completely defined by other parameters (eg it is not an 'independent parameter')
      mutable CrystVector_REAL mPowderPatternCalc;
      mutable std::map<RefinablePar*,CrystVector_REAL> mPowderPattern_FullDeriv;
      /// The calculated powder pattern, integrated
      mutable CrystVector_REAL mPowderPatternIntegratedCalc;
      mutable std::map<RefinablePar*,CrystVector_REAL> mPowderPatternIntegrated_FullDeriv;
      /// The calculated powder pattern part which corresponds to 'background'
      /// (eg non-scalable components). It is already included in mPowderPatternCalc
      mutable CrystVector_REAL mPowderPatternBackgroundCalc;
      /// The calculated powder pattern part which corresponds to 'background'
      /// (eg non-scalable components), integrated
      mutable CrystVector_REAL mPowderPatternBackgroundIntegratedCalc;
      /// The observed powder pattern.
      CrystVector_REAL mPowderPatternObs;
      /// The sigma of the observed pattern.
      CrystVector_REAL mPowderPatternObsSigma;
      /// The weight for each point of the pattern.
      mutable CrystVector_REAL mPowderPatternWeight;
      /// The complete variance associated to each point of the powder pattern,
      /// taking into account observation and model errors.
      mutable CrystVector_REAL mPowderPatternVariance;
      /// The complete variance associated to each point of the powder pattern,
      /// taking into account observation and model errors. Integrated.
      mutable CrystVector_REAL mPowderPatternVarianceIntegrated;
      /// The cumulative Chi^2 (integrated or not, depending on the option)
      mutable CrystVector_REAL mChi2Cumul;


      /// The calculated powder pattern. Cropped to the maximum sin(theta)/lambda for LSQ
      mutable CrystVector_REAL mPowderPatternUsedCalc;
      mutable std::map<RefinablePar*,CrystVector_REAL> mPowderPatternUsed_FullDeriv;
      /// The calculated powder pattern. Cropped to the maximum sin(theta)/lambda for LSQ
      mutable CrystVector_REAL mPowderPatternUsedObs;
      /// The weight for each point of the pattern. Cropped to the maximum sin(theta)/lambda for LSQ
      mutable CrystVector_REAL mPowderPatternUsedWeight;

      /** Vector of x coordinates (either 2theta or time-of-flight) for the pattern
      *
      * Stored in ascending order for 2theta, and descending for TOF, i.e. always
      * in ascending order for the corresponding sin(theta)/lambda.
      */
      CrystVector_REAL mX;
      /// Is the mX vector sorted in ascending order ? (true for 2theta, false for TOF)
      bool mIsXAscending;
      /// Number of points in the pattern
      unsigned long mNbPoint;

      /// The Radiation corresponding to this experiment
      Radiation mRadiation;

      // Clocks
         /// When were the pattern parameters (2theta or time-of-flight range) changed ?
         RefinableObjClock mClockPowderPatternPar;
         /// When were the radiation parameter (radiation type, wavelength) changed ?
         RefinableObjClock mClockPowderPatternRadiation;
         /// When was the powder pattern last computed ?
         mutable RefinableObjClock mClockPowderPatternCalc;
         /// When was the powder pattern (integrated) last computed ?
         mutable RefinableObjClock mClockPowderPatternIntegratedCalc;
         /// Corrections to 2Theta
         RefinableObjClock mClockPowderPatternXCorr;
         /// Last modification of the scale factor
         mutable RefinableObjClock mClockScaleFactor;

      //Excluded regions in the powder pattern, for statistics.
         /// Min coordinate for for all excluded regions
         CrystVector_REAL mExcludedRegionMinX;
         /// Max coordinate for 2theta for all excluded regions
         CrystVector_REAL mExcludedRegionMaxX;

      //Various corrections to 2theta-to be used by the components
         /// Zero correction :
         /// \f$ (2\theta)_{obs} = (2\theta)_{real} +(2\theta)_{0}\f$
         ///Thus mPowderPattern2ThetaMin=(mPowderPattern2ThetaMin-m2ThetaZero)
         REAL mXZero;
         /// Displacement correction :
         ///\f$ (2\theta)_{obs} = (2\theta)_{real} + \frac{a}{\cos(\theta)} \f$
         REAL m2ThetaDisplacement;
         /// Transparency correction :
         ///\f$ (2\theta)_{obs} = (2\theta)_{real} + b\sin(2\theta) \f$
         REAL m2ThetaTransparency;
         /// Time Of Flight (TOF) parameters :
         ///\f$ t = DIFC*\frac{\sin(\theta)}{\lambda} + DIFA*\left(\frac{\sin(\theta)}{\lambda}\right)^2 + mXZero\f$
         REAL mDIFC,mDIFA;
      // Components of the powder pattern
         /// The components (crystalline phases, background,...) of the powder pattern
         ObjRegistry<PowderPatternComponent> mPowderPatternComponentRegistry;
         /// The scale factors for each component. For unscalable phases,
         /// this is set to 1 (constant).
         ///
         /// This is mutable because generally we use the 'best' scale factor, but
         /// it should not be...
         mutable CrystVector_REAL mScaleFactor;

      /// Use faster, less precise functions ?
      bool mUseFastLessPreciseFunc;

      // For statistics
         /// Should Statistics (R, Rw,..) exclude the background ?
         bool mStatisticsExcludeBackground;
         /// \internal To compute scale factors, which are the components (phases) that
         /// can be scaled ?
         mutable CrystVector_int mScalableComponentIndex;
         /// \internal Used to fit the components' scale factors
         mutable CrystMatrix_REAL mFitScaleFactorM,mFitScaleFactorB,mFitScaleFactorX;

      /// Use Integrated profiles for Chi^2, R, Rwp...
         RefObjOpt mOptProfileIntegration;

      // Integrated R-factors
         mutable CrystVector_long mIntegratedPatternMin,mIntegratedPatternMax;
         mutable CrystVector_REAL mIntegratedObs;
         mutable CrystVector_REAL mIntegratedWeight;
         mutable CrystVector_REAL mIntegratedWeightObs;
         mutable CrystVector_REAL mIntegratedVarianceObs;
         mutable RefinableObjClock mClockIntegratedFactorsPrep;
      // Statistical indicators
         mutable REAL mChi2,mIntegratedChi2;
         /// This is the logarithm of the part of log(Likelihood) which corresponds
         /// to the normalization terms of gaussian distribution for each obs/calc
         /// point. In practice, this is the sum of 1/2*log(2pi*sig(i)^2), although
         /// we discard the 2pi terms.
         mutable REAL mChi2LikeNorm,mIntegratedChi2LikeNorm;
         mutable REAL mR;
         mutable REAL mRw;
         ///Clock the last time Chi^2 was computed
         mutable RefinableObjClock mClockChi2,mClockIntegratedChi2;
      /** Maximum sin(theta)/lambda for all calculations (10 by default).
      *
      * This keeps all data in memory, but only the part which is below
      * the max is calculated.
      */
      REAL mMaxSinThetaOvLambda;
      /// Number of points actually used, due to the maximum value of
      /// sin(theta)/lambda.
      mutable unsigned long mNbPointUsed;
      /// Number of integration intervals actually used, due to the maximum value of
      /// sin(theta)/lambda.
      mutable unsigned long mNbIntegrationUsed;
      /// Clock recording the last time the number of points used (PowderPattern::mNbPointUsed)
      /// was changed.
      mutable RefinableObjClock mClockNbPointUsed;
   #ifdef __WX__CRYST__
   public:
      virtual WXCrystObjBasic* WXCreate(wxWindow*);
      friend class WXPowderPattern;
      // This should be removed
      friend class WXPowderPatternGraph;
   #endif
};
/// Global registry for all PowderPattern objects
extern ObjRegistry<PowderPattern> gPowderPatternRegistry;
//######################################################################
//    PROFILE FUNCTIONS (for powder diffraction)
//######################################################################

///Gaussian, normalized (ie integral is equal to 1), as a function of theta
/// and of the FWHM. The input is an array of the theta values. The maximum of the
///function is in theta=0. If asymmetry is used, negative tth values must be first.
CrystVector_REAL PowderProfileGauss  (const CrystVector_REAL theta,
                                      const REAL fwhm,
                                      const REAL asymmetryPar=1.);
///Lorentzian, normalized (ie integral is equal to 1), as a function of theta
/// and of the FWHM. The input is an array of the theta values. The maximum of the
///function is in theta=0. If asymmetry is used, negative tth values must be first.
CrystVector_REAL PowderProfileLorentz(const CrystVector_REAL theta,
                                      const REAL fwhm,
                                      const REAL asymmetryPar=1.);


}//namespace ObjCryst
#endif // _OBJCRYST_POWDERPATTERN_H_
