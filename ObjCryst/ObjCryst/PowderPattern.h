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

#include "CrystVector/CrystVector.h"
#include "ObjCryst/General.h"
//#include "ObjCryst/SpaceGroup.h"
//#include "ObjCryst/ScatteringPower.h"
//#include "ObjCryst/Scatterer.h"
#include "ObjCryst/Crystal.h"
#include "ObjCryst/ScatteringCorr.h"

//#include <stdlib.h>
#include <string>
//#include <iomanip>
//#include <cmath>
//#include <typeinfo>
//#include <fstream>
//#include <ctime>

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
      
      /// \internal
      /// Set the PowderPattern object which uses this component.
      /// This sets all necessary spectrum parameters (2theta range,
      /// wavelength, radiation type...) accordingly.
      /// 
      virtual void SetParentPowderPattern(const PowderPattern&)=0;
      /// Get the calculated powder spectrum for this component.
      /// Note that the spectrum is \e not scaled.
      /// 
      virtual const CrystVector_REAL& GetPowderPatternCalc()const=0;
      /** \brief Is this component scalable ?
      *
      * This is used by the PowderPattern class, which fits all
      * spectrum components using scale factors. Some components may not
      * need to be scaled: background components, which are assumed
      * to be absolute.
      */
      bool IsScalable()const;
   protected:
      /// Last time the spectrum was actually calculated.
      const RefinableObjClock& GetClockPowderPatternCalc()const;
         
      /// Calc the powder spectrum. As always, recomputation is only
      /// done if necessary (ie if a parameter has changed since the last
      /// computation)
      virtual void CalcPowderPattern() const=0;
      
      /// Get the integration limits (first and last pixels) around each reflection,
      /// if this component has Bragg reflections. Used for integrated R(w) factors.
      /// The limits currently go from -2*FWHM to +2*FWHM.
      /// returns a pointer to the min and max pixels arrays (null pointers if
      /// no reflection for this phase).
      virtual void GetBraggLimits(CrystVector_long *&min,CrystVector_long *&max)const=0;
      /// Get last time the Bragg Limits were changed
      const RefinableObjClock& GetClockBraggLimits()const;
      
       /// Set the maximum value for sin(theta)/lambda. All data above still
      /// exist but are ignored for all calculations.
      virtual void SetMaxSinThetaOvLambda(const REAL max)=0;
      
      /// The calculated component of a powder spectrum. It is mutable since it is
      /// completely defined by other parameters (eg it is not an 'independent parameter')
      mutable CrystVector_REAL mPowderPatternCalc;
      
      /// \internal
      /// This will be called by the parent PowderPattern object, before
      /// calculating the first powder spectrum. Or maybe it should be called
      /// automatically by the object itself...
      virtual void Prepare()=0;
      
      /// Scalable ? (crystal phase = scalable, background= not scalable)
      bool mIsScalable;
      
      //Clocks
         /// When was the powder spectrum last computed ?
         mutable RefinableObjClock mClockPowderPatternCalc;
      
      /// The PowderPattern object in which this component is included
      const PowderPattern *mpParentPowderPattern;
      /// Get last time the Bragg Limits were changed
      mutable RefinableObjClock mClockBraggLimits;
      
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
      
      virtual void SetParentPowderPattern(const PowderPattern&);
      virtual const CrystVector_REAL& GetPowderPatternCalc()const;
      /// Import background points from a file (with two columns 2theta, intensity)
      void ImportUserBackground(const string &filename);
      virtual void XMLOutput(ostream &os,int indent=0)const;
      virtual void XMLInput(istream &is,const XMLCrystTag &tag);
      //virtual void XMLInputOld(istream &is,const IOCrystTag &tag);
      virtual void GetGeneGroup(const RefinableObj &obj, 
                                CrystVector_uint & groupIndex,
                                unsigned int &firstGroup) const;
   protected:
      virtual void CalcPowderPattern() const;
      virtual void Prepare();
      virtual void GetBraggLimits(CrystVector_long *&min,CrystVector_long *&max)const;
      virtual void SetMaxSinThetaOvLambda(const REAL max);
      
      /// The kind of interpolation used
      PowderBackgroundInterpType mBackgroundType;
      /// Number of fitting points for background
      int mBackgroundNbPoint;
      /// Vector of 2theta values for the fitting points of the background
      CrystVector_REAL mBackgroundInterpPoint2Theta;
      /// Values of background at interpolating points
      CrystVector_REAL mBackgroundInterpPointIntensity;
      
      // Clocks
         /// Modification of the interpolated points
         RefinableObjClock mClockBackgroundPoint;
      // mutable clocks
         /// Last time Splines were generated
         mutable RefinableObjClock mClockBackgroundSpline;
      
      /** Maximum sin(theta)/lambda for all calculations (10 by default).
      *
      * This keeps all data in memory, but only the part which is below
      * the max is calculated.
      */
      REAL mMaxSinThetaOvLambda;
      
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
      
      virtual void SetParentPowderPattern(const PowderPattern&);
      virtual const CrystVector_REAL& GetPowderPatternCalc()const;
      
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
      virtual const Radiation& GetRadiation()const;
   protected:
      virtual void CalcPowderPattern() const;
      
      /// \internal Calc reflection profiles for ALL reflections (powder diffraction)
      void CalcPowderReflProfile()const;
      /// \internal Calc Lorentz-Polarisation-APerture correction
      void CalcIntensityCorr()const;
      /// \internal Compute the intensity for all reflections (taking into account
      /// corrections, but not the multiplicity)
      virtual void CalcIhkl() const;
      virtual void Prepare();
      virtual void InitOptions();
      virtual void GetBraggLimits(CrystVector_long *&min,CrystVector_long *&max)const;
      virtual void SetMaxSinThetaOvLambda(const REAL max);
      //Clocks
         /// Last time the reflection parameters were changed
         RefinableObjClock mClockProfilePar;
         /// Last time the 
         RefinableObjClock mClockLorentzPolarSlitCorrPar;
      //Clocks (internal, mutable)
         /// Last time the Lorentz-Polar-Slit correction was computed
         mutable RefinableObjClock mClockIntensityCorr;
         /// Last time the reflection profiles were computed
         mutable RefinableObjClock mClockProfileCalc;
         /// Last time intensities were computed
         mutable RefinableObjClock mClockIhklCalc;
      //Profile parameters
         ///Gaussian ? Lorentzian ? Pseudo-Voigt ?
         RefObjOpt mReflectionProfileType;
         /// The 'full' profile of each reflection os taken equal
         /// to the FWHM, multiplied by this factor.
         REAL mFullProfileWidthFactor;
         ///FWHM parameters, following Caglioti's law
         REAL mCagliotiU,mCagliotiV,mCagliotiW;
         ///Pseudo-Voigt mixing parameter : eta=eta0 +2*theta*eta1
         /// eta=1 -> pure Lorentzian ; eta=0 -> pure Gaussian
         REAL mPseudoVoigtEta0,mPseudoVoigtEta1;
         /// Use asymmetric profiles ? (unused yet)
         bool mUseAsymmetricProfile;
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
         
      /// Computed intensities for all reflections
         mutable CrystVector_REAL mIhklCalc;
      
      // Saved arrays to speed-up computations
         ///Reflection profiles for ALL reflections during the last powder spectrum generation
         mutable CrystMatrix_REAL mSavedPowderReflProfile;
         /// \internal Number of points used to describe each individual profile
         mutable long mSavedPowderReflProfileNbPoint;
         /// \internal The 1st pixel for each reflection
         mutable CrystVector_long mReflectionProfileFirstPixel;
         /// First and last pixel for integrated R-factors around each reflection
         mutable CrystVector_long mIntegratedReflMin,mIntegratedReflMax;
   #ifdef __WX__CRYST__
   public:
      virtual WXCrystObjBasic* WXCreate(wxWindow*);
      friend class WXPowderPatternDiffraction;
   #endif
};


//######################################################################
/** \brief Powder pattern class, with an observed pattern and several
* calculated components to modelize the spectrum.
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
      /** Add a component (phase, backround) to this spectrum. 
      *
      * It must have been allocated in the heap. The spectrum parameters (2theta min,
      * step, nbpoints, wavelength, radiation type) of the component 
      * are automatically changed to that of the PowderPattern object.
      */
      void AddPowderPatternComponent(PowderPatternComponent &);
      /// Number of components
      unsigned int GetNbPowderPatternComponent()const;
      /// Access to a component of the powder spectrum
      const PowderPatternComponent& GetPowderPatternComponent(const string &name)const;
      /// Access to a component of the powder spectrum
      const PowderPatternComponent& GetPowderPatternComponent(const int)const;
      /// Access to a component of the powder spectrum
      PowderPatternComponent& GetPowderPatternComponent(const string &name);
      /// Access to a component of the powder spectrum
      PowderPatternComponent& GetPowderPatternComponent(const int);
      
      // Pattern parameters (2theta range, wavelength, radiation)
         /** \briefSet the powder spectrum angular range & resolution parameter.
         * this will affect all components (phases) of the spectrum.
         *
         *   Use this with caution, as the number of points must be correct with
         * respect to the observed data (Iobs).
         *
         * \param tthetaMin: min 2theta value, in radians
         * \param tthetaStep: step (assumed constant) in 2theta.
         * \param nbPoints: number of points in the spectrum.
         */
         void SetPowderPatternPar(const REAL tthetaMin,
                                           const REAL tthetaStep,
                                           unsigned long nbPoint);
         ///Number of points ?
         unsigned long GetNbPoint()const;
         
         /// Set the radiation
         void SetRadiation(const Radiation &radiation);
         ///Neutron or x-ray experiment ?
         const Radiation& GetRadiation()const;
         
         /// Set the radiation type
         void SetRadiationType(const RadiationType radiation);
         ///Neutron or x-ray experiment ?
         RadiationType GetRadiationType()const;
         ///Set the wavelength of the experiment (in Angstroems).
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
      
      //Access to spectrum data
         /// Get the calculated powder spectrum
         const CrystVector_REAL& GetPowderPatternCalc()const;
         /// Get the observed powder spectrum
         const CrystVector_REAL& GetPowderPatternObs()const;
         /// Get the sigma for each point of the observed powder spectrum 
         const CrystVector_REAL& GetPowderPatternObsSigma()const;
         /// Get the weight for each point of the powder spectrum 
         const CrystVector_REAL& GetPowderPatternWeight()const;
         // Get the 2theta values corresponding to the powder spectrum.
         //Values are returned in radian, and are experimental 2theta values.
         //CrystVector_REAL Get2Theta()const;
         /// Get the Minimum 2theta
         REAL Get2ThetaMin()const;
         /// Get the step in 2theta
         REAL Get2ThetaStep()const;
         /// Get the maximum 2theta
         REAL Get2ThetaMax()const;
      
      // Clocks
         /// Last time the spectrum was calculated
         const RefinableObjClock& GetClockPowderPatternCalc()const;
         /// When were the spectrum parameters (2theta range, step) changed ?
         const RefinableObjClock& GetClockPowderPatternPar()const;
         /// When were the radiation parameter (radiation type, wavelength) changed ?
         const RefinableObjClock& GetClockPowderPatternRadiation()const;
         /// When were the parameters for 2theta correction (zero, transparency,
         /// displacement) last changed ?
         const RefinableObjClock& GetClockPowderPattern2ThetaCorr()const;
      
      // Corrections to 2theta
         ///Change Zero in 2Theta
         void Set2ThetaZero(const REAL newZero);
         /// Change displacement correction
         /// \f$ (2\theta)_{obs} = (2\theta)_{real} + \frac{a}{\cos(\theta)} \f$
         void Set2ThetaDisplacement(const REAL displacement);
         ///Change transparency correction
         /// \f$ (2\theta)_{obs} = (2\theta)_{real} + b\sin(2\theta) \f$
         void Set2ThetaTransparency(const REAL transparency);
         /// Get the experimental 2theta from the theoretical value, taking
         /// into account all corrections (zero, transparency,..).
         /// \internal
         /// \param ttheta: the theoretical 2theta value.
         /// \return the 2theta value as it appears on the spectrum.
         REAL Get2ThetaCorr(const REAL ttheta)const;
         /// Get the pixel number on the experimental spectrum, from the
         /// theoretical (uncorrected) value of 2theta, taking into account all corrections.
         /// (zero, transparency,..).
         /// \internal
         /// \param ttheta: the theoretical 2theta value.
         /// \return the 2theta value as it appears on the spectrum.
         long Get2ThetaCorrPixel(const REAL ttheta)const;

      // Import & export powder spectrum
         /** \brief Import fullprof-style diffraction data.
         *\param fullprofFileName: filename
         */
         void ImportPowderPatternFullprof(const string &fullprofFileName);
         /** \brief Import powder spectrum, format DMC from PSI
         */
         void ImportPowderPatternPSI_DMC(const string &filename);
         /** \brief Import powder spectrum, format from ILL D1A/D2B 
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
         *\warning the 2theta step is assumed to be constant !!!
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
         *\warning the 2theta step is assumed to be constant !!!
         *\param fileName: the filename (surprise!)
         *\param nbSkip: the number of lines to skip at the beginning of the file (default=0)
         */
         void ImportPowderPattern2ThetaObs(const string &fileName,const int nbSkip=0);
         /** \brief Set observed powder spectrum from vector array.
         *
         * Note: powder spectrum parameters must have been set before calling this function,
         * for example by calling DiffractionDataPowder::InitPowderPatternPar().
         */
         void SetPowderPatternObs(const CrystVector_REAL& obs);
         
         ///Save powder spectrum to one file, text format, 3 columns theta Iobs Icalc.
         ///If Iobs is missing, the column is omitted.
         ///
         /// \todo export in other formats (.prf,...), with a list of reflection
         /// position for all phases...
         void SavePowderPattern(const string &filename="powderPattern.out") const;
         /// Print to thee screen/console the observed and calculated spectrum (long,
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
         REAL GetChiSq()const;
         /// Fit the scale(s) factor of each component to minimize R
         void FitScaleFactorForR();
         void FitScaleFactorForIntegratedR();
         /// Fit the scale(s) factor of each component to minimize Rw
         void FitScaleFactorForRw();
         void FitScaleFactorForIntegratedRw();
         /// Set sigma=sqrt(Iobs)
         void SetSigmaToSqrtIobs();
         /// Set w = 1/sigma^2.
         ///
         /// To filter too small or null intensities :If sigma< minRelatSigma* max(sigma),
         /// then w=1/(minRelatSigma* max(sigma))^2
         void SetWeightToInvSigmaSq(const REAL minRelatSigma=1e-3);
         /// Set w = sin(theta). Not really usful, huh ?
         void SetWeightToSinTheta(const REAL power=1.);
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
         void Add2ThetaExcludedRegion(const REAL min2Theta,const REAL max2theta);
         
      //Cost functions
         /// Number of Cost functions
         unsigned int GetNbCostFunction()const;
         const string& GetCostFunctionName(const unsigned int)const;
         const string& GetCostFunctionDescription(const unsigned int)const;
         virtual REAL GetCostFunctionValue(const unsigned int);
      //LSQ functions
         virtual unsigned int GetNbLSQFunction()const;
         virtual const CrystVector_REAL& GetLSQCalc(const unsigned int) const;
         virtual const CrystVector_REAL& GetLSQObs(const unsigned int) const;
         virtual const CrystVector_REAL& GetLSQWeight(const unsigned int) const;
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
   protected:
      /// Calc the powder spectrum
      void CalcPowderPattern() const;
      /// Init parameters and options
      virtual void Init();
      /// Prepare  the calculation of the integrated R-factors
      void PrepareIntegratedRfactor()const;
      /// Calculate the number of points of the pattern actually used, from the maximum
      /// value of sin(theta)/lambda
      void CalcNbPointUsed()const;
      
      /// The calculated powder spectrum. It is mutable since it is
      /// completely defined by other parameters (eg it is not an 'independent parameter')
      mutable CrystVector_REAL mPowderPatternCalc;
      /// The calculated powder spectrum part which corresponds to 'background' 
      /// (eg non-scalable components). It is already included in mPowderPatternCalc
      mutable CrystVector_REAL mPowderPatternBackgroundCalc;
      /// The observed powder spectrum.
      CrystVector_REAL mPowderPatternObs;
      /// The sigma of the observed spectrum.
      CrystVector_REAL mPowderPatternObsSigma;
      /// The weight for each point of the spectrum.
      CrystVector_REAL mPowderPatternWeight;
      
      /// 2theta min and step for the spectrum
      REAL m2ThetaMin,m2ThetaStep;
      /// Number of points in the spectrum
      unsigned long mNbPoint;
      
      /// The wavelength of the experiment, in Angstroems.
      /// \warning This should be removed, as it is also available in mRadiation.
      REAL mWavelength;
      /// The wavelength of the experiment, in Angstroems.
      Radiation mRadiation;
      
      // Clocks
         /// When were the spectrum parameters (2theta range, step) changed ?
         RefinableObjClock mClockPowderPatternPar;
         /// When were the radiation parameter (radiation type, wavelength) changed ?
         RefinableObjClock mClockPowderPatternRadiation;
         /// When was the powder spectrum last computed ?
         mutable RefinableObjClock mClockPowderPatternCalc;
         /// Corrections to 2Theta
         RefinableObjClock mClockPowderPattern2ThetaCorr;
         /// Last modification of the scale factor
         RefinableObjClock mClockScaleFactor;
      
      //Excluded 2Theta regions in the powder spectrum, for statistics.
         /// Min value for 2theta for all excluded regions
         CrystVector_REAL mExcludedRegionMin2Theta;
         /// Max value for 2theta for all excluded regions
         CrystVector_REAL mExcludedRegionMax2Theta;
  
      //Various corrections to 2theta-to be used by the components
         /// Zero correction :
         /// \f$ (2\theta)_{obs} = (2\theta)_{real} +(2\theta)_{0}\f$
         ///Thus mPowderPattern2ThetaMin=(mPowderPattern2ThetaMin-m2ThetaZero)
         REAL m2ThetaZero;
         /// Displacement correction :
         ///\f$ (2\theta)_{obs} = (2\theta)_{real} + \frac{a}{\cos(\theta)} \f$
         REAL m2ThetaDisplacement;
         /// Transparency correction : 
         ///\f$ (2\theta)_{obs} = (2\theta)_{real} + b\sin(2\theta) \f$
         REAL m2ThetaTransparency;
      // Components of the powder spectrum
         /// The components (crystalline phases, background,...) of the powder spectrum
         ObjRegistry<PowderPatternComponent> mPowderPatternComponentRegistry;
         /// The scale factors for each component. For unscalable phases,
         /// this is set to 1 (constant).
         CrystVector_REAL mScaleFactor;
         
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
         
      // Integrated R-factors
         mutable CrystVector_long mIntegratedPatternMin,mIntegratedPatternMax;
         mutable CrystVector_REAL mIntegratedObs;
         mutable CrystVector_REAL mIntegratedWeight;
         mutable RefinableObjClock mClockIntegratedFactorsPrep;
      /** Maximum sin(theta)/lambda for all calculations (10 by default).
      *
      * This keeps all data in memory, but only the part which is below
      * the max is calculated.
      */
      REAL mMaxSinThetaOvLambda;
      /// Number of points actually used, due to the maximum value of
      /// sin(theta)/lambda.
      mutable unsigned long mNbPointUsed;
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
