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
/*   ScatteringPowerSphere.h
*  header file for the Spherical scattering power (fullerenes, etc..)
*
*/

#ifndef _OBJCRYST_SCATTPOWERSPHERE_H_
#define _OBJCRYST_SCATTPOWERSPHERE_H_

#include "ObjCryst/CrystVector/CrystVector.h"
#include "ObjCryst/ObjCryst/ScatteringPower.h"

namespace ObjCryst
{
/** \ brief ScatteringPower for a spherical particule
*
* This can be used to modelize the form factor of disordered (or low-resolution)
* of fullerene-type compounds, where all atoms are located on a
* sphere.
*
* This actually modelizes a spherical distribution of a \e single electron, so to modelize
* C60 the occupancy must be set to 60*6.
*/
class ScatteringPowerSphere: public ScatteringPower
{
   public:
      /// Default constructor
      ScatteringPowerSphere();
      /**   \brief constructor
      *  \param name : name of the ScatteringPower ('C60','France 98'...).
      *The name can have \e any format 
      * \param nbAtom: the number of atoms 
      *  \param biso : Isotropic thermic coefficient
      * \param AxisLengthX,AxisLengthY,AxisLengthZ: length of the different
      * main axis of the ellipsoid
      * \param symbol: the symbol of the element associated to this fullerene. By
      * default it is assumed to be carbon
      */
      ScatteringPowerSphere(const string &name,const REAL radius,const REAL bIso=1.0);
      ScatteringPowerSphere(const ScatteringPowerSphere& old);
      ~ScatteringPowerSphere();
      void Init(const string &name,const REAL radius,const REAL bIso=1.0);
      virtual const string& GetClassName() const;
      virtual CrystVector_REAL GetScatteringFactor(const ScatteringData &data,
                                                     const int spgSymPosIndex=0) const;
      virtual REAL GetForwardScatteringFactor(const RadiationType) const;
      virtual CrystVector_REAL GetTemperatureFactor(const ScatteringData &data,
                                                     const int spgSymPosIndex=0) const;
      virtual CrystMatrix_REAL GetResonantScattFactReal(const ScatteringData &data,
                                                     const int spgSymPosIndex=0) const;
      virtual CrystMatrix_REAL GetResonantScattFactImag(const ScatteringData &data,
                                                     const int spgSymPosIndex=0) const;
      REAL GetRadius()const;
      virtual void Print()const;
      virtual void XMLOutput(ostream &os,int indent=0)const;
      virtual void XMLInput(istream &is,const XMLCrystTag &tag);
   private:
      virtual void InitRefParList();
      /// Isotropic temperature B-factor.
      REAL mBiso;
      /** Radius of the sphere.
      */
      REAL mRadius;
   #ifdef __WX__CRYST__
   public:
      virtual WXCrystObjBasic* WXCreate(wxWindow*);
   friend class WXScatteringPowerSphere;
   #endif
};

}//namespace
#endif
