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
*/
//######################################################################

class DiffractionDataSingleCrystal:public ScatteringData
{
   public:
      DiffractionDataSingleCrystal();
      ~DiffractionDataSingleCrystal();
      virtual DiffractionDataSingleCrystal* CreateCopy()const;
      virtual const string GetClassName() const;
      
      /**  \brief returns the calculated diffracted intensity.
      *
      * This is an array of calculated intensities for each reflections in the
      * single crystal case, and the array with the full powder powder profile
      * for powder diffraction.
      */
      const CrystVector_double& GetIcalc() const;

      /// Return the array of observed intensities for all peaks
      const CrystVector_double& GetIobs() const;
      /// Return the array of observed intensities for all peaks
      void SetIobs(const CrystVector_double&);
      
      /// Return the array of sigmas for observed intensities, for all peaks.
      const CrystVector_double& GetSigma() const;
      /// Return the array of sigmas for observed intensities, for all peaks.
      void SetSigma(const CrystVector_double&);
      
      /// Set Iobs to current values of Icalc. Mostly used for tests.
      void SetIobsToIcalc();
      
      /// Return the weights (for each reflection) used for computing Rw.
      const CrystVector_double& GetWeight() const;
      /// Change the weights (for each reflection) used for computing Rw.
      void SetWeight(const CrystVector_double&);
      
      /** \brief input H,K,L, Iobs and Sigma
      *
      * \param h,k,l: double arrays (vectors with NbRefl elements -same size)
      *with the h, k and l coordinates of all reflections.
      * \param iobs,sigma: double arrays (vectors with NbRefl elements -same size)
      *with the Observed intensity and sigma for all reflections.
      *
      */
      void SetHklIobs(  CrystVector_long const &h,
                        CrystVector_long const &k,
                        CrystVector_long const &l,
                        CrystVector_double const &iObs,
                        CrystVector_double const &sigma);
      
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
      virtual double GetRw()const;
      /** \brief  Return the Crystal R-factor (unweighted)
      *
      * \return \f$ R= \sqrt {\frac{\sum_i \left( I_i^{obs}-I_i^{calc} \right)^2}
      * {\sum_i (I_i^{obs})^2} }\f$
      */
      virtual double GetR()const;
      
      /** \brief  Return  conventionnal Chi^2
      *
       \return \f$ \chi^2 = \sum_i w_i \left(I_i^{obs}-I_i^{calc} \right)^2
      *  \f$
      */
      // \return \f$ \chi^2 = \sum_i \left( \frac{ y_i^{obs}-y_i^{calc}}{\sigma_i} \right)^2
      virtual double GetChi2()const;
      
      /* \brief  Return  Goodness of Fit
      *
      * \return \f$ GoF =  \frac{\sum_i w_i\left( y_i^{obs}-y_i^{calc} \right)^2}
      * {N_{obs}-N_{indep.par}} \f$
      virtual double GoF()const;
      */
      
      /** Compute the best scale factor minimising Rw.
      *
      * The computed scale factor is \e immediatly applied to Icalc
      */
      virtual void CalcBestScaleFactorForRw();
      /** Compute the best scale factor minimising R.
      *
      * The computed scale factor is \e immediatly applied to Icalc
      */
      virtual void CalcBestScaleFactorForR();
      /// Compute the best scale factor to minimize R, apply this scale factor and return
      /// the R value obtained.
      virtual double GetBestRFactor();

      /** \brief Set sigma for all observed intensities to sqrt(obs)
      *
      */
      virtual void SetSigmaToSqrtIobs();
      /** \brief Set the weight for all observed intensities to 1/sigma^2
      *
      *For sigmas which are smaller than minRelatSigma times the max value of sigma,
      *the output weight is set to 0.
      */
      virtual void SetWeightToInvSigma2(const double minRelatSigma=1e-4);
      
      /// Scale factor (applied to Icalc to match Iobs)
      double GetScaleFactor()const;
      // Set the Scale factor (applied to Icalc to match Iobs)
      //void SetScaleFactor(const double);

      /** \brief Print H, K, L Iobs sigma for all reflections
      *
      */
      virtual void PrintObsData()const;
      /** \brief Print H, K, L Iobs sigma Icalc for all reflections
      *Iobs and sigma (if given) are scaled to Icalc (if available).
      *
      */
      virtual void PrintObsCalcData()const;
      
      virtual void SetUseOnlyLowAngleData(const bool useOnlyLowAngle,const double angle=0.);
      
      ///Save H,K,L Iobs Icalc to a file, text format, 3 columns theta Iobs Icalc.
      ///If Iobs is missing, the column is omitted.
      void SaveHKLIobsIcalc(const string &filename="hklIobsIcalc.out");
      //Cost functions
         unsigned int GetNbCostFunction()const;
         const string& GetCostFunctionName(const unsigned int)const;
         const string& GetCostFunctionDescription(const unsigned int)const;
         virtual double GetCostFunctionValue(const unsigned int);
      virtual void Output(ostream &os,int indent=0)const;
      virtual void Input(istream &is,const XMLCrystTag &tag);
      virtual void InputOld(istream &is,const IOCrystTag &tag);
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
      CrystVector_double mObsIntensity ;  
      /// Sigma for observed intensities (either individual reflections or spectrum)
      CrystVector_double mObsSigma ;
      /// weight for computing R-Factor, for each observed value.
      CrystVector_double mWeight ;
      /// Calculated intensities
      mutable CrystVector_double mCalcIntensity ;  
      /// Scale factor. It is applied when computing intensities. The scale
      ///applies to intensities
      double mScaleFactor;
};
/// Global registry for all PowderPattern objects
extern ObjRegistry<DiffractionDataSingleCrystal> gDiffractionDataSingleCrystalRegistry;

} //namespace ObjCryst

#endif //_OBJCRYST_DIFFDATA_SINGLECRYSTAL_H_
