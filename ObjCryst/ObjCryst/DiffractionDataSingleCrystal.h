/* 
* ObjCryst++ : a Crystallographic computing library in C++
*			http://objcryst.sourceforge.net
*			http://www.ccp14.ac.uk/ccp/web-mirrors/objcryst/
*
*  (c) 2000-2001 Vincent FAVRE-NICOLIN vincefn@users.sourceforge.net
*
*/
/*   DiffractionDataSingleCrystal.h header file for single-crystal
* diffraction objects.
*
*/
#ifndef _OBJCRYST_DIFFDATA_SINGLECRYSTAL_H_
#define _OBJCRYST_DIFFDATA_SINGLECRYSTAL_H_

#include "CrystVector/CrystVector.h"

#include "ObjCryst/General.h"

#include "ObjCryst/ScatteringPower.h"
#include "ObjCryst/Crystal.h"
#include "ObjCryst/ScatteringData.h"
//#include "ObjCryst/SpaceGroup.h"
//#include "ObjCryst/Scatterer.h"

//#include <stdlib.h>
#include <string>
//#include <iomanip>
//#include <cmath>
//#include <typeinfo>
//#include <fstream>
//#include <ctime>

namespace ObjCryst
{
//######################################################################
//                   DiffractionDataSingleCrystal
/**
* \brief DiffractionData object for Single Crystal analysis.
*
* Currently this handles only in the simplest way single crystal dat: ie
* only data which has been completely corrected for Lorentz/Polarization
* and absorption.
*
* What needs to be developped: define the geometry of the experiment
* (incident and emerging angles), the polarization of the beam, etc...
*/
//######################################################################

class DiffractionDataSingleCrystal:public ScatteringData
{
   public:
      DiffractionDataSingleCrystal();
		/** Constructor, with an assigned crystal structure.
		*
		*/
      DiffractionDataSingleCrystal(Crystal &cryst);
		/// Copy constructor
      DiffractionDataSingleCrystal(const DiffractionDataSingleCrystal &old);
      ~DiffractionDataSingleCrystal();
      virtual DiffractionDataSingleCrystal* CreateCopy()const;
      virtual const string GetClassName() const;
      
      /**  \brief returns the calculated diffracted intensity.
      *
      * This is an array of calculated intensities for each reflections in the
      * single crystal case, and the array with the full powder powder profile
      * for powder diffraction.
      */
      const CrystVector_REAL& GetIcalc() const;

      /// Return the array of observed intensities for all peaks
      const CrystVector_REAL& GetIobs() const;
      /// Return the array of observed intensities for all peaks
      void SetIobs(const CrystVector_REAL&);
      
      /// Return the array of sigmas for observed intensities, for all peaks.
      const CrystVector_REAL& GetSigma() const;
      /// Return the array of sigmas for observed intensities, for all peaks.
      void SetSigma(const CrystVector_REAL&);
      
      /// Set Iobs to current values of Icalc. Mostly used for tests.
      void SetIobsToIcalc();
      
      /// Return the weights (for each reflection) used for computing Rw.
      const CrystVector_REAL& GetWeight() const;
      /// Change the weights (for each reflection) used for computing Rw.
      void SetWeight(const CrystVector_REAL&);
      
      /** \brief input H,K,L, Iobs and Sigma
      *
      * \param h,k,l: REAL arrays (vectors with NbRefl elements -same size)
      *with the h, k and l coordinates of all reflections.
      * \param iobs,sigma: REAL arrays (vectors with NbRefl elements -same size)
      *with the Observed intensity and sigma for all reflections.
      *
      */
      void SetHklIobs(  CrystVector_long const &h,
                        CrystVector_long const &k,
                        CrystVector_long const &l,
                        CrystVector_REAL const &iObs,
                        CrystVector_REAL const &sigma);
      
      /** \brief Import h,k,l,I form a file
      *
      *The file is assumed to correspond to a single crystal diffraction file.
      *  \param fileName The name of the data file. This file should be formatted
      *with H,k,l, Iobs separated by spaces.
      *  \param nbRefl The number of reflections to extract.
      *  \param skipLines The number of lines to skip at the beginning of the file.
      */
      void ImportHklIobs(const string &fileName,const long nbRefl,const int skipLines=0);
      /** \brief Import h,k,l,I,Sigma form a file
      *
      *The file is assumed to correspond to a single crystal diffraction file.
      *  \param fileName The name of the data file. This file should be formatted
      *with H,k,l, Iobs and Sigma separated by spaces.
      *  \param nbRefl The number of reflections to extract.
      *  \param skipLines The number of lines to skip at the beginning of the file.
      */
      void ImportHklIobsSigma(const string &fileName,const long nbRefl,const int skipLines=0);
      /** \brief Import h,k,l,I,Sigma form a Jana98 '*.m91' file
      *
      *The file is assumed to correspond to a single crystal diffraction file.
      *  \param fileName The name of the data file.
      */
      void ImportHklIobsSigmaJanaM91(const string &fileName);
      
      /** \brief  Return the Crystal R-factor (weighted)
      *
      * \return \f$ R_{w}= \sqrt {\frac{\sum_i w_i\left( I_i^{obs}-I_i^{calc} \right)^2}
      * {\sum_i w_i (I_i^{obs})^2} }\f$
      */
      virtual REAL GetRw()const;
      /** \brief  Return the Crystal R-factor (unweighted)
      *
      * \return \f$ R= \sqrt {\frac{\sum_i \left( I_i^{obs}-I_i^{calc} \right)^2}
      * {\sum_i (I_i^{obs})^2} }\f$
      */
      virtual REAL GetR()const;
      
      /** \brief  Return  conventionnal Chi^2
      *
       \return \f$ \chi^2 = \sum_i w_i \left(I_i^{obs}-I_i^{calc} \right)^2
      *  \f$
      */
      // \return \f$ \chi^2 = \sum_i \left( \frac{ y_i^{obs}-y_i^{calc}}{\sigma_i} \right)^2
      virtual REAL GetChi2()const;
      
      /* \brief  Return  Goodness of Fit
      *
      * \return \f$ GoF =  \frac{\sum_i w_i\left( y_i^{obs}-y_i^{calc} \right)^2}
      * {N_{obs}-N_{indep.par}} \f$
      virtual REAL GoF()const;
      */
      
      /** Compute the best scale factor minimising Rw.
      *
      * The computed scale factor is \e immediatly applied to Icalc
      */
      virtual void FitScaleFactorForRw();
      /** Compute the best scale factor minimising R.
      *
      * The computed scale factor is \e immediatly applied to Icalc
      */
      virtual void FitScaleFactorForR();
      /// Compute the best scale factor to minimize R, apply this scale factor and return
      /// the R value obtained.
      virtual REAL GetBestRFactor();

      /** \brief Set sigma for all observed intensities to sqrt(obs)
      *
      */
      virtual void SetSigmaToSqrtIobs();
      /** \brief Set the weight for all observed intensities to 1/sigma^2
      *
      *For sigmas which are smaller than minRelatSigma times the max value of sigma,
      *the output weight is set to 0.
      */
      virtual void SetWeightToInvSigma2(const REAL minRelatSigma=1e-4);
      
      /// Scale factor (applied to Icalc to match Iobs)
      REAL GetScaleFactor()const;
      // Set the Scale factor (applied to Icalc to match Iobs)
      //void SetScaleFactor(const REAL);

      /** \brief Print H, K, L Iobs sigma for all reflections
      *
      */
      virtual void PrintObsData()const;
      /** \brief Print H, K, L Iobs sigma Icalc for all reflections
      *Iobs and sigma (if given) are scaled to Icalc (if available).
      *
      */
      virtual void PrintObsCalcData()const;
      
      virtual void SetUseOnlyLowAngleData(const bool useOnlyLowAngle,const REAL angle=0.);
      
      ///Save H,K,L Iobs Icalc to a file, text format, 3 columns theta Iobs Icalc.
      ///If Iobs is missing, the column is omitted.
      void SaveHKLIobsIcalc(const string &filename="hklIobsIcalc.out");
      //Cost functions
         unsigned int GetNbCostFunction()const;
         const string& GetCostFunctionName(const unsigned int)const;
         const string& GetCostFunctionDescription(const unsigned int)const;
         virtual REAL GetCostFunctionValue(const unsigned int);
      //LSQ functions
         virtual unsigned int GetNbLSQFunction()const;
         virtual const CrystVector_REAL& GetLSQCalc(const unsigned int) const;
         virtual const CrystVector_REAL& GetLSQObs(const unsigned int) const;
         virtual const CrystVector_REAL& GetLSQWeight(const unsigned int) const;
      virtual void XMLOutput(ostream &os,int indent=0)const;
      virtual void XMLInput(istream &is,const XMLCrystTag &tag);
      //virtual void XMLInputOld(istream &is,const IOCrystTag &tag);
   protected:
   private:
      virtual void InitRefParList();
      /// Calc intensities
      void CalcIcalc() const;
      /// Are there observed intensities ?
      bool mHasObservedData;

      /** \brief Observed intensity (after ABS and LP corrections)
      *
      * In the single crystal case, this is a list of intensity corresponding to (h,k,l).
      * For a powder sample, this is a list of all peaks intensities.
      */
      CrystVector_REAL mObsIntensity ;  
      /// Sigma for observed intensities (either individual reflections or spectrum)
      CrystVector_REAL mObsSigma ;
      /// weight for computing R-Factor, for each observed value.
      CrystVector_REAL mWeight ;
      /// Calculated intensities
      mutable CrystVector_REAL mCalcIntensity ;  
      /// Scale factor. It is applied when computing intensities. The scale
      ///applies to intensities
      REAL mScaleFactor;
   #ifdef __WX__CRYST__
   public:
      virtual WXCrystObjBasic* WXCreate(wxWindow*);
      //friend class WXDiffractionSinglCrystal;

   #endif
};
/// Global registry for all PowderPattern objects
extern ObjRegistry<DiffractionDataSingleCrystal> gDiffractionDataSingleCrystalRegistry;

} //namespace ObjCryst

#endif //_OBJCRYST_DIFFDATA_SINGLECRYSTAL_H_
