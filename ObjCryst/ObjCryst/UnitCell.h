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
/*   UnitCell.h header file for the UnitCell object
*
*/
#ifndef _OBJCRYST_UNITCELL_H_
#define _OBJCRYST_UNITCELL_H_

#include "CrystVector/CrystVector.h"

#include "ObjCryst/General.h"
#include "RefinableObj/RefinableObj.h"
#include "ObjCryst/SpaceGroup.h"


namespace ObjCryst
{
extern const RefParType *gpRefParTypeCrystal;// Defined in Crystal.cpp
extern const RefParType *gpRefParTypeUnitCell;
extern const RefParType *gpRefParTypeUnitCellLength;
extern const RefParType *gpRefParTypeUnitCellAngle;
class NiftyStaticGlobalObjectsInitializer_UnitCell
{
   public:
      NiftyStaticGlobalObjectsInitializer_UnitCell()
      {
         if (mCount++ == 0)
         {
            gpRefParTypeUnitCell=new RefParType (gpRefParTypeCrystal,"Unit Cell");
            gpRefParTypeUnitCellLength=new RefParType (gpRefParTypeUnitCell,"Unit Cell Length");
            gpRefParTypeUnitCellAngle=new RefParType (gpRefParTypeUnitCell,"Unit Cell Angle");
         }
      }
      ~NiftyStaticGlobalObjectsInitializer_UnitCell()
      {
         if (--mCount == 0)
         {
            delete gpRefParTypeUnitCell;
            delete gpRefParTypeUnitCellLength;
            delete gpRefParTypeUnitCellAngle;
            gpRefParTypeUnitCell=0;
            gpRefParTypeUnitCellLength=0;
            gpRefParTypeUnitCellAngle=0;
         }
      }
   private:
      static long mCount;
};
static NiftyStaticGlobalObjectsInitializer_UnitCell NiftyStaticGlobalObjectsInitializer_UnitCell_counter;

//######################################################################
/** \brief Unit Cell class: Unit cell with  spacegroup information.
*
*/
//######################################################################
class UnitCell:public RefinableObj
{
   public:
      /// Default Constructor
      UnitCell();
      /** \brief UnitCell Constructor (orthorombic)
      *  \param a,b,c : unit cell dimension, in angstroems
      *  \param SpaceGroupId: space group symbol or number
      */
      UnitCell(const REAL a, const REAL b, const REAL c,
              const string &SpaceGroupId);
      /** \brief UnitCell Constructor (triclinic)
      *  \param a,b,c : unit cell dimension, in angstroems
      *  \param alpha,beta,gamma : unit cell angles, in radians.
      *  \param SpaceGroupId: space group symbol or number
      */
      UnitCell(const REAL a, const REAL b, const REAL c, const REAL alpha,
              const REAL beta, const REAL gamma,const string &SpaceGroupId);
              
      /// UnitCell copy constructor
      UnitCell(const UnitCell &oldCryst);
      /// Destructor
      ~UnitCell();
      virtual const string& GetClassName() const;
      /// Lattice parameters (a,b,c,alpha,beta,gamma) as a 6-element vector in Angstroems 
      /// and radians.
      CrystVector_REAL GetLatticePar() const;
      /// Return one of the 6 Lattice parameters, 0<= whichPar <6 (a,b,c,alpha,beta,gamma),
      /// returned in Angstroems and radians. 
      REAL GetLatticePar(const int whichPar)const;
      /// last time the Lattice parameters were changed
      const RefinableObjClock& GetClockLatticePar()const;
      /** \brief Get the 'B' matrix (UnitCell::mBMatrix)for the UnitCell (orthogonalization 
      * matrix for the given lattice, in the reciprocal space)
      *
      * The convention is taken following Giacovazzo, "Fundamentals of Crystallography", p.69
      * "e1 is chosen along a*, e2 in the (a*,b*) plane, then e3 is along c".
      */
      const CrystMatrix_REAL& GetBMatrix() const;
      /** \brief Get the orthogonalization matrix (UnitCell::mOrthMatrix)for the UnitCell 
      * in real space
      *
      */
      const CrystMatrix_REAL& GetOrthMatrix() const;
      /// last time the metric matrices were changed
      const RefinableObjClock& GetClockMetricMatrix()const;
      /** \brief Get orthonormal cartesian coordinates for a set of (x,y,z)
      * fractional coordinates.
      *
      * Results are given in Angstroems.
      * The convention is taken following :
      * e1 is chosen along a, e2 in the (a,b) plane, then e3 is along c*
      */
      CrystVector_REAL GetOrthonormalCoords(const REAL x,const REAL y,const REAL z) const;
      /** \brief Get orthonormal cartesian coordinates for a set of (x,y,z)
      * fractional coordinates.
      *
      * X,y and z input are changed to Amgstroems values
      * The convention is taken following :
      * e1 is chosen along a, e2 in the (a,b) plane, then e3 is along c*
      */
      void FractionalToOrthonormalCoords(REAL &x,REAL &y,REAL &z) const;
      /** \brief Get fractional cartesian coordinates for a set of (x,y,z)
      * orthonormal coordinates.
      *
      * Result is stored into x,y and z
      * The convention is taken following :
      * e1 is chosen along a, e2 in the (a,b) plane, then e3 is along c*
      */
      void OrthonormalToFractionalCoords(REAL &x,REAL &y,REAL &z) const;
      /** \brief Get Miller H,K, L indices from orthonormal coordinates 
      * in reciprocal space.
      *
      * Result is stored into x,y and z
      */
      void MillerToOrthonormalCoords(REAL &x,REAL &y,REAL &z) const;
      /** \brief Get orthonormal coordinates given a set of H,K, L indices
      * in reciprocal space.
      *
      * Result is stored into x,y and z
      */
      void OrthonormalToMillerCoords(REAL &x,REAL &y,REAL &z) const;
      /** Prints some info about the UnitCell
      *
      * \param os the stream to which the information is outputed (default=cout)
      */
      virtual void Print(ostream &os=cout) const;
            
      /// Access to the SpaceGroup object
      const SpaceGroup & GetSpaceGroup()const;
      /// Access to the SpaceGroup object.
      SpaceGroup & GetSpaceGroup();
      
      // :TODO: ?
      //virtual void XMLOutput(ostream &os,int indent=0)const;
      //virtual void XMLInput(istream &is,const XMLCrystTag &tag);
      
      /// Volume of Unit Cell (in Angstroems)
      REAL GetVolume()const;
   protected:
      /** \brief Init all UnitCell parameters
      *  \param a,b,c : unit cell dimension, in angstroems
      *  \param alpha,beta,gamma : unit cell angles
      *  \param SpcGroup: space group number (1..230)
      *  \param name: name for the UnitCell, : '(TaSe4)2I'
      */
      virtual void Init(const REAL a, const REAL b, const REAL c, const REAL alpha,
                const REAL beta, const REAL gamma,const string &SpaceGroupId,
                const string& name);
      /** \brief Prepare the refinable parameters list
      *
      * This is called once when creating the UnitCell.
      */
      void InitRefParList();
   private:
      /** Init options.
      *
      * Need only be done once per UnitCell.
      */
      virtual void InitOptions();
      /// \internal.Init the (de)orthogonalization matrices. 
      /// They are re-computed only if parameters have changed since last call.
      void InitMatrices() const;
      /// \internal
      /// Update cell parameters for tetragonal, trigonal, hexagonal, cubic lattices.
      /// Also set angular parameters for those group which need it. This is needed during
      /// Refinement, since for example in a quadratic spg, only a is refined and
      /// we need to have b=a...
      void UpdateLatticePar();

      /// a,b and c in Angstroems, angles (stored) in radians
      /// For cubic, rhomboedric UnitCells, only the 'a' parameter is relevant.
      /// For quadratic and hexagonal UnitCells, only a and c parameters are relevant.
      /// The MUTABLE is temporary ! It should not be !
      CrystVector_REAL mCellDim;
      /// The space group of the UnitCell
      SpaceGroup mSpaceGroup ;
      
      /** \brief B Matrix (Orthogonalization matrix for reciprocal space)
      * \f[ B= \left[ \begin {array}{ccc} a^* & b^*\cos(\gamma^*) & c^*\cos(\beta^*) \\
      *                            0 & b^*\sin(\gamma^*) & -c^*\sin(\beta^*)\cos(\alpha) \\
      *                            0 & 0 & \frac{1}{c}\end{array} \right]\f]
      *\f[ \left[ \begin{array}{c} k_x \\ k_y \\ k_z \end{array} \right]_{orthonormal}  
      *  = B \times \left[ \begin{array}{c} h \\ k \\ l \end{array}\right]_{integer} \f]
      * \note this matrix is and must remain upper triangular. this is assumed for
      * some optimizations.
      */
      mutable CrystMatrix_REAL mBMatrix;
      /// inverse of B Matrix (i.e. inverse of orthogonalization matrix for direct space)
      mutable CrystMatrix_REAL mBMatrixInvert;
      /** \brief Eucl Matrix (Orthogonalization matrix for direct space)
      * \f[ M_{orth}= \left[ \begin {array}{ccc} a & b\cos(\gamma) & c\cos(\beta) \\
      *                            0 & b\sin(\gamma) & -c\sin(\beta)\cos(\alpha^*) \\
      *                            0 & 0 & \frac{1}{c^*}\end{array} \right]\f]
      * \f[ \left[ \begin{array}{c} x \\ y \\ z \end{array} \right]_{orthonormal}  
      *  = M_{orth} \times \left[ \begin{array}{c} x \\ y \\ z \end{array}
      *                    \right]_{fractional coords}\f]
      * \note this matrix is and must remain upper triangular. this is assumed for
      * some optimizations.
      */
      mutable CrystMatrix_REAL mOrthMatrix;
      /// inverse of Eucl Matrix (i.e. inverse of de-orthogonalization matrix for direct space)
      mutable CrystMatrix_REAL mOrthMatrixInvert;

      /// Last time lattice parameters were changed
      RefinableObjClock mClockLatticePar;
      /// \internal Last time metric matrices were computed
      mutable RefinableObjClock mClockMetricMatrix;
      /// \internal Last time the lattice parameters whre updated
      mutable RefinableObjClock mClockLatticeParUpdate;
      
      /** Option to override lattice parameters constraints from spacegroup choice.
      *
      * \warning EXPERIMENTAL
      *
      * Normally lattice parameters are constrained by the space group choice
      * (e.g. a=b=c and angles =90° for cubic spacegroups). Using this option
      * allows you to override this, and choose any lattice parameter. THis works
      * as long as symmetry operations are applied to fractionnal coordinates.
      *
      * This is useful duting global optimization when searching the structure in a UnitCell
      * which has (or is expected to have) a known pseudo-crystallographic
      * symmetry, to reduce dramatically the number of parameters. Of course
      * for final refinement the 'real' symmetry should be imposed.
      */
      RefObjOpt mConstrainLatticeToSpaceGroup;
      
};


}// namespace


#endif //_OBJCRYST_UNITCELL_H_
