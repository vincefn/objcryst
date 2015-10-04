/*  ObjCryst++ Object-Oriented Crystallographic Library
    (c) 2000-2005 Vincent Favre-Nicolin vincefn@users.sourceforge.net
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
#ifndef _OBJCRYST_REFLECTIONPROFILE_H_
#define _OBJCRYST_REFLECTIONPROFILE_H_

#include <complex>
#include "ObjCryst/CrystVector/CrystVector.h"
#include "ObjCryst/ObjCryst/General.h"
#include "ObjCryst/ObjCryst/UnitCell.h"
#include "ObjCryst/RefinableObj/RefinableObj.h"

namespace ObjCryst
{
//######################################################################
//    basic functions used for profiles
//######################################################################

///Gaussian, normalized (ie integral is equal to 1), as a function of theta
/// and of the FWHM. The input is an array of the theta values. The maximum of the
///function is in theta=center. If asymmetry is used, negative tth values must be first.
CrystVector_REAL PowderProfileGauss  (const CrystVector_REAL theta,
                                      const REAL fwhm, const REAL center, const REAL asym=1.0);
///Lorentzian, normalized (ie integral is equal to 1), as a function of theta
/// and of the FWHM. The input is an array of the theta values. The maximum of the
///function is in theta=center. If asymmetry is used, negative tth values must be first.
CrystVector_REAL PowderProfileLorentz(const CrystVector_REAL theta,
                                      const REAL fwhm, const REAL center, const REAL asym=1.0);
/// Asymmetry function [Ref J. Appl. Cryst 26 (1993), 128-129
CrystVector_REAL AsymmetryBerarBaldinozzi(const CrystVector_REAL theta,
                                          const REAL fwhm, const REAL center,
                                          const REAL A0, const REAL A1,
                                          const REAL B0, const REAL B1);

/** Complex exponential integral E1(z) (Abramowitz & Stegun,  chap. 5)
*
* Using A&S 5.1.11 (series) and 5.1.51 (asymptotic) expansions
*/
template <class T> std::complex<T>ExponentialIntegral1(const complex<T> z);
/** E1(z)*exp(z)
*
* This can be computed for large z values to avoid floating-point exceptions.
*/
template <class T> std::complex<T>ExponentialIntegral1_ExpZ(const complex<T> z);

/** Abstract base class for reflection profiles.
*
*/
class ReflectionProfile:public RefinableObj
{
   public:
      ReflectionProfile();
      ReflectionProfile(const ReflectionProfile &old);
      virtual ~ReflectionProfile();
      virtual ReflectionProfile* CreateCopy()const=0;
      /** Get the reflection profile
      *
      *\param x: the vector of x  coordinates (i.e. either 2theta or time-of-flight)
      *\param xcenter: coordinate (2theta, tof) of the center of the peak
      *\param h,k,l: reflection Miller indices
      *\note: derived classes who need either d_hkl or the orthonormal coordinates
      * of the scattering vector should be passed a ObjCryst::UnitCell object in
      * the constructor so that they can use ObjCryst::UnitCell::MillerToOrthonormalCoords()
      */
      virtual CrystVector_REAL GetProfile(const CrystVector_REAL &x, const REAL xcenter,
                                  const REAL h, const REAL k, const REAL l)const=0;
      /// Get the (approximate) full profile width at a given percentage
      /// of the profile maximum (e.g. FWHM=GetFullProfileWidth(0.5)).
      virtual REAL GetFullProfileWidth(const REAL relativeIntensity, const REAL xcenter,
                                       const REAL h, const REAL k, const REAL l)=0;
      /// Is the profile anisotropic ?
      virtual bool IsAnisotropic()const;
      virtual void XMLOutput(ostream &os,int indent=0)const=0;
      virtual void XMLInput(istream &is,const XMLCrystTag &tag)=0;
   private:
#ifdef __WX__CRYST__
   public:
      virtual WXCrystObjBasic* WXCreate(wxWindow* parent)=0;
#endif
};

/** Pseudo-Voigt reflection profile.
*
*/
class ReflectionProfilePseudoVoigt:public ReflectionProfile
{
   public:
      ReflectionProfilePseudoVoigt();
      ReflectionProfilePseudoVoigt(const ReflectionProfilePseudoVoigt &old);
      virtual ~ReflectionProfilePseudoVoigt();
      virtual ReflectionProfilePseudoVoigt* CreateCopy()const;
      virtual const string& GetClassName()const;
      CrystVector_REAL GetProfile(const CrystVector_REAL &x, const REAL xcenter,
                                  const REAL h, const REAL k, const REAL l)const;
      /** Set reflection profile parameters
      *
      * \param fwhmCagliotiW,fwhmCagliotiU,fwhmCagliotiV : these are the U,V and W
      * parameters in the Caglioti's law :
      * \f$ fwhm^2= U \tan^2(\theta) + V \tan(\theta) +W \f$
      * if only W is given, the width is constant
      * \param eta0,eta1: these are the mixing parameters.
      */
      void SetProfilePar(const REAL fwhmCagliotiW,
                         const REAL fwhmCagliotiU=0,
                         const REAL fwhmCagliotiV=0,
                         const REAL eta0=0.5,
                         const REAL eta1=0.);
      virtual REAL GetFullProfileWidth(const REAL relativeIntensity, const REAL xcenter,
                                       const REAL h, const REAL k, const REAL l);
      bool IsAnisotropic()const;
      virtual void XMLOutput(ostream &os,int indent=0)const;
      virtual void XMLInput(istream &is,const XMLCrystTag &tag);
   private:
      /// Initialize parameters
      void InitParameters();
      ///FWHM parameters, following Caglioti's law
      REAL mCagliotiU,mCagliotiV,mCagliotiW;
      ///Pseudo-Voigt mixing parameter : eta=eta0 +2*theta*eta1
      /// eta=1 -> pure Lorentzian ; eta=0 -> pure Gaussian
      REAL mPseudoVoigtEta0,mPseudoVoigtEta1;
      /** Asymmetry parameters, following the Bérar \& Baldinozzi approach
      * ( Bérar \& baldinozzi, J. Appl. Cryst 26 (1993), 128-129)
      *
      * \note: these are not used right now.
      */
      REAL mAsymBerarBaldinozziA0,mAsymBerarBaldinozziA1,
           mAsymBerarBaldinozziB0,mAsymBerarBaldinozziB1;
      /** Asymmetry parameters, following the analytical function for asymmetric pseudo-voigt
      * given by Toraya in J. Appl. Cryst 23(1990),485-491
      *
      * Asymmetric coefficient:
      * \f[ A=A_0+A_1/\sin(2\vartheta)+A_2/\sin^2(2\vartheta) \f]
      *
      * Asymmetric profile:
      * \f[ Prof(\vartheta-\vartheta_{max}<=0)=Prof_0(\frac{1+A}{A})(\vartheta-\vartheta_{max})) \f]
      * \f[ Prof(\vartheta-\vartheta_{max}>=0)=Prof_0((1+A)         (\vartheta-\vartheta_{max})) \f]
      */
      REAL mAsym0,mAsym1,mAsym2;
#ifdef __WX__CRYST__
   public:
      virtual WXCrystObjBasic* WXCreate(wxWindow* parent);
#endif
};

/** Pseudo-Voigt reflection profile, with 6-parameters anisotropic Lorentzian broadening and Toraya asymmetric modelling.
 *
 */
class ReflectionProfilePseudoVoigtAnisotropic:public ReflectionProfile
{
   public:
      ReflectionProfilePseudoVoigtAnisotropic();
      ReflectionProfilePseudoVoigtAnisotropic(const ReflectionProfilePseudoVoigtAnisotropic &old);
      virtual ~ReflectionProfilePseudoVoigtAnisotropic();
      virtual ReflectionProfilePseudoVoigtAnisotropic* CreateCopy()const;
      virtual const string& GetClassName()const;
      CrystVector_REAL GetProfile(const CrystVector_REAL &x, const REAL xcenter,
                                  const REAL h, const REAL k, const REAL l)const;
      /** Set reflection profile parameters
       *
       * if only W is given, the width is constant
       *
       * See documentation for mCagliotiU,mCagliotiV,mCagliotiW,mScherrerP,mLorentzX,mLorentzY,
       * mLorentzGammaHH,mLorentzGammaKK,mLorentzGammaLL,mLorentzGammaHK,mLorentzGammaHL,mLorentzGammaKL
       * mAsym0,mAsym1,mAsym2.
       */
      void SetProfilePar(const REAL fwhmCagliotiW,
                         const REAL fwhmCagliotiU=0,
                         const REAL fwhmCagliotiV=0,
                         const REAL fwhmGaussP=0,
                         const REAL fwhmLorentzX=0,
                         const REAL fwhmLorentzY=0,
                         const REAL fwhmLorentzGammaHH=0,
                         const REAL fwhmLorentzGammaKK=0,
                         const REAL fwhmLorentzGammaLL=0,
                         const REAL fwhmLorentzGammaHK=0,
                         const REAL fwhmLorentzGammaHL=0,
                         const REAL fwhmLorentzGammaKL=0,
                         const REAL pseudoVoigtEta0=0,
                         const REAL pseudoVoigtEta1=0,
                         const REAL asymA0=0,
                         const REAL asymA1=0,
                         const REAL asymA2=0
                         );
      virtual REAL GetFullProfileWidth(const REAL relativeIntensity, const REAL xcenter,
                                       const REAL h, const REAL k, const REAL l);
      bool IsAnisotropic()const;
      virtual void XMLOutput(ostream &os,int indent=0)const;
      virtual void XMLInput(istream &is,const XMLCrystTag &tag);
   private:
      /// Initialize parameters
      void InitParameters();
      /** FWHM parameters:
       *
       * Profile:\f$ Prof(\vartheta) = \eta Lorentz(\vartheta) + (1-\eta) Gauss(\vartheta) \f$
       *
       * This pseudo-Voigt profile is the linear combination of a Gaussian and a Lorentzian:
       * \f[ Profile(\vartheta)= Lorentz(\vartheta) + Gauss(\vartheta)\f]
       *
       * Gaussian FWHM parameters, following Caglioti's law and Scherrer broadening coefficient:
       * \f[ fwhm_g^2= U \tan^2(\theta) + V \tan(\theta) +W +\frac{P}{\cos^2\vartheta} \f]
       * See Cagliotti, Pauletti & Ricci, Nucl. Instrum. 3 (1958), 223.
       *
       * Lorentzian FWHM parameters
       *
       * \f[ fwhm_l= \frac{X}{\cos\vartheta} +(Y + \frac{\gamma_L}{\sin^2\vartheta})\tan\vartheta \f]
       * with anisotropic broadening factor:
       * \f[ \gamma_L = \gamma_{HH}h^2 + \gamma_{KK}k^2 + \gamma_{LL}l^2 + 2\gamma_{HK}hk + 2\gamma_{hl} + 2\gamma_{kl}kl\f]
       *
       * \note: in the above formula we use \f$ \frac{\gamma_L}{\sin^2\vartheta} \f$ rather than the more usual\f$ \gamma_L d^2 \f$
       * for computing purposes (the model being purely phenomenological).
       */
      REAL mCagliotiU,mCagliotiV,mCagliotiW,mScherrerP,mLorentzX,mLorentzY,
           mLorentzGammaHH,mLorentzGammaKK,mLorentzGammaLL,mLorentzGammaHK,mLorentzGammaHL,mLorentzGammaKL;
      ///Pseudo-Voigt mixing parameter : \f$ \eta=\eta_0 +2*\vartheta*\eta_1\f$ 
      /// eta=1 -> pure Lorentzian ; eta=0 -> pure Gaussian
      REAL mPseudoVoigtEta0,mPseudoVoigtEta1;
      /** Asymmetry parameters, following the analytical function for asymmetric pseudo-voigt
       * given by Toraya in J. Appl. Cryst 23(1990),485-491
       *
       * Asymmetric coefficient:
       * \f[ A=A_0+A_1/\sin(2\vartheta)+A_2/\sin^2(2\vartheta) \f]
       *
       * Asymmetric profile:
       * \f[ Prof(\vartheta-\vartheta_{max}<=0)=Prof_0(\frac{1+A}{A})(\vartheta-\vartheta_{max})) \f]
       * \f[ Prof(\vartheta-\vartheta_{max}>=0)=Prof_0((1+A)         (\vartheta-\vartheta_{max})) \f]
       */
      REAL mAsym0,mAsym1,mAsym2;
   #ifdef __WX__CRYST__
   public:
      virtual WXCrystObjBasic* WXCreate(wxWindow* parent);
   #endif
};

/** Double-Exponential Pseudo-Voigt profile for TOF.
*
* Ref Mark Pitt
*/
class ReflectionProfileDoubleExponentialPseudoVoigt:public ReflectionProfile
{
   public:
      /** Constructor, without unit cell
      */
      ReflectionProfileDoubleExponentialPseudoVoigt();
      /** Constructor, with unit cell
      */
      ReflectionProfileDoubleExponentialPseudoVoigt(const UnitCell &cell);
      ReflectionProfileDoubleExponentialPseudoVoigt
         (const ReflectionProfileDoubleExponentialPseudoVoigt &old);
      virtual ~ReflectionProfileDoubleExponentialPseudoVoigt();
      virtual ReflectionProfileDoubleExponentialPseudoVoigt* CreateCopy()const;
      virtual const string& GetClassName()const;
      CrystVector_REAL GetProfile(const CrystVector_REAL &x, const REAL xcenter,
                                  const REAL h, const REAL k, const REAL l)const;
      /** Set reflection profile parameters
      *
      */
      void SetProfilePar(const REAL instrumentAlpha0,
                         const REAL instrumentAlpha1,
                         const REAL instrumentBeta0,
                         const REAL instrumentBeta1,
                         const REAL gaussianSigma0,
                         const REAL gaussianSigma1,
                         const REAL gaussianSigma2,
                         const REAL lorentzianGamma0,
                         const REAL lorentzianGamma1,
                         const REAL lorentzianGamma2);
      virtual REAL GetFullProfileWidth(const REAL relativeIntensity, const REAL xcenter,
                                       const REAL h, const REAL k, const REAL l);
      bool IsAnisotropic()const;
      virtual void XMLOutput(ostream &os,int indent=0)const;
      virtual void XMLInput(istream &is,const XMLCrystTag &tag);
      /** Set unit cell. This is used to compute d_hkl for reflections.
      *
      */
      void SetUnitCell(const UnitCell &cell);
   private:
      /// Initialize parameters
      void InitParameters();
      REAL mInstrumentAlpha0;
      REAL mInstrumentAlpha1;
      REAL mInstrumentBeta0;
      REAL mInstrumentBeta1;
      REAL mGaussianSigma0;
      REAL mGaussianSigma1;
      REAL mGaussianSigma2;
      REAL mLorentzianGamma0;
      REAL mLorentzianGamma1;
      REAL mLorentzianGamma2;
      const UnitCell *mpCell;
#ifdef __WX__CRYST__
   public:
      virtual WXCrystObjBasic* WXCreate(wxWindow* parent);
#endif
};
/// Global registry for all ReflectionProfile objects
extern ObjRegistry<ReflectionProfile> gReflectionProfileRegistry;
}//namespace
#endif
