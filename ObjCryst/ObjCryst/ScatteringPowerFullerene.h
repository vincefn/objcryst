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
/*   ScatteringPowerFullerene.h
*  header file for the Fullerene-type scattering power
*
*/

#ifndef _OBJCRYST_SCATTPOWERFULLERENE_H_
#define _OBJCRYST_SCATTPOWERFULLERENE_H_

#include "CrystVector/CrystVector.h"
#include "ObjCryst/ScatteringPower.h"

namespace ObjCryst
{
/** \ brief ScatteringPower for a fullerene-type particule (spherical or
* ellipsoidal.
*
* This can be used to modelize the form factor of disordered (or low-resolution)
* of fullerene-type compounds, where all atoms are located on an ellipsoid
* or sphere.
*
*/
class ScatteringPowerFullerene: public ScatteringPower
{
   public:
      /// Default constructor
      ScatteringPowerFullerene();
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
      ScatteringPowerFullerene(const string &name,
                               const unsigned int nbAtom,
                               const REAL axisLengthX,
                               const REAL axisLengthY,
                               const REAL axisLengthZ,
                               const REAL bIso=1.0,
                               const string symbol="C");
      ScatteringPowerFullerene(const ScatteringPowerFullerene& old);
      ~ScatteringPowerFullerene();
      void Init(const string &name,
                const unsigned int nbAtom,
                const REAL axisLengthX,
                const REAL axisLengthY,
                const REAL axisLengthZ,
                const REAL bIso=1.0,
                const string symbol="C");
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
      REAL GetAxisLengthX() const;
      REAL GetAxisLengthY() const;
      REAL GetAxisLengthZ() const;
      virtual REAL GetRadius()const;
      virtual void Print()const;
      virtual void XMLOutput(ostream &os,int indent=0)const;
      virtual void XMLInput(istream &is,const XMLCrystTag &tag);
   private:
      virtual void InitRefParList();
      /// Number of atoms on the ellipsoid
      REAL mNbAtom;
      /** Lengths of the ellipsoid 3 axis
      * The X,Y,Z axis correspond to the XYZ directions when the orientation
      * angles are equal to 0.
      */
      REAL mAxisLengthX,mAxisLengthY,mAxisLengthZ;
      /** \brief Angles giving the orientation of the ellipsoid (stored in radian)
      *
      * \f[     \left[ \begin{array}{ccc} \cos(\chi) & 0 & -\sin(\chi) \\
                                           0 & 1 & 0 \\
                                           \sin(\chi) & 0 & \cos(\chi) \end{array} \right]
         \times \left[ \begin{array}{ccc} \cos(\phi) & -\sin(\phi) & 0 \\
                                           \sin(\phi) & \cos(\phi) & 0 \\
                                           0 & 0 & 1 \end{array} \right]
         \times \left[ \begin{array}{ccc} 1 & 0 & 0 \\
                                           0 & \cos(\psi) & -\sin(\psi) \\
                                           0 & \sin(\psi) & \cos(\psi) \end{array} \right]
      * \f]
      *
      *
      */
      REAL mPhi,mChi,mPsi;
      /// The scattering power of a single element
      ScatteringPower *mpScatteringPower;
   #ifdef __WX__CRYST__
   public:
      virtual WXCrystObjBasic* WXCreate(wxWindow*);
   friend class WXScatteringPowerFullerene;
   #endif
};

}//namespace
#endif
