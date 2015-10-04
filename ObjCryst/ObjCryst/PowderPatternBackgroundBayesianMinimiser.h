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
/*   PowderPatternBackgroundBayesianMinimiser.h
*  Bayesian estimation of powder pattern background
*
*/
#ifndef __POWDERPATTERNBACKGROUNDBAYESIANMINIMISER_H
#define __POWDERPATTERNBACKGROUNDBAYESIANMINIMISER_H

#include "ObjCryst/CrystVector/CrystVector.h"
#include "ObjCryst/RefinableObj/RefinableObj.h"
#include "ObjCryst/ObjCryst/PowderPattern.h"

namespace ObjCryst
{
/** This object is used to estimate the background in a powder pattern,
* using a Bayesian approach (David & Sivia, Acta Cryst A50 (1994), 703)
*
*/
class PowderPatternBackgroundBayesianMinimiser:public RefinableObj
{
   public:
      PowderPatternBackgroundBayesianMinimiser(PowderPatternBackground &backgd);
      ~PowderPatternBackgroundBayesianMinimiser();
      virtual const string& GetClassName()const;
      virtual REAL GetLogLikelihood()const;
      virtual unsigned int  GetNbLSQFunction () const ;
      virtual const CrystVector_REAL &  GetLSQCalc (const unsigned int) const;
      virtual const CrystVector_REAL &  GetLSQObs (const unsigned int) const;
      virtual const CrystVector_REAL &  GetLSQWeight (const unsigned int) const;
   //private:
   /** Returns the log(likelihood) of a background by marginalizing the effect of Bragg peaks,
   following the method described by David and Sivia (\e J.Appl.Cryst. \b 34(2001), 318).
   *
   * \returns
   * - \f$ L(t) = \frac{\left(y^{calc}_{Background}-y^{obs}\right)^2}{2\sigma^2}\f$
   * for \t<0
   * - \f$ L(t) =  A-\log{\int_\epsilon^{\infty} {\frac{1}{u} e^{-(t-u)^2} du}}\f$ for t>0,
   * with: \f$\epsilon = 1e-6\f$, \f$ t = \frac{y^{calc}_{Background}-y^{obs}}{\sigma\sqrt{2}}\f$
   * and A a normalizing constant so that the function is continuous for t=0 (i.e. \e L(0)=0).
   *
   * \param \f$ t = \frac{y^{calc}_{Background}-y^{obs}}{\sigma\sqrt{2}}\f$
   *
   * For small \e t>0 values, \e L(t) behaves like a quadratic function, and for large
   * positive values it is equivalent to \f$ \log{\frac{\sqrt{\pi}}{t}}\f$.
   *
   * As the integral diverges for \f$\epsilon=0\f$, it is necessary to use
   * a small, non-null \f$\epsilon\f$. The use of a smaller \f$\epsilon\f$ changes
   * the range over which the function behaves quadratically, as well as the maximum
   * value (at \e t=0), but does not affect the asymptotic value.
   *
   * See tabulated values in the source code. The function is approximated with a cubic
   * spline until \e t=8, and then uses the asymptotic\f$A-\log{\frac{\sqrt{\pi}}{t}}\f$ value.
   *
   * \note For a more strict calculation, we should include a normalizing constant (?)
   */
   static REAL BayesianBackgroundLogLikelihood(const REAL t);
   PowderPatternBackground *mpBackground;
   /// Bayesian cost (-log(likelihood)) for each point
   mutable CrystVector_REAL mBayesianCalc;
   /// Obs==0 (desired -log(likelihood))
   mutable CrystVector_REAL mBayesianObs;
   /// Weight=1
   mutable CrystVector_REAL mBayesianWeight;
};

}//namespace
#endif //__POWDERPATTERNBACKGROUNDBAYESIANMINIMISER_H
