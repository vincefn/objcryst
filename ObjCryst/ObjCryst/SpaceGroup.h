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
/*   Spacegroup.h header for Spacegroup and AsymmetricUnit classes
*
*/
#ifndef _OBJCRYST_SPACEGROUP_H_
#define _OBJCRYST_SPACEGROUP_H_

#include "CrystVector/CrystVector.h"

#include "ObjCryst/General.h"
#include "RefinableObj/RefinableObj.h"


namespace ObjCryst
{

//R. Gross-Kunstleve "sglite/sglite.h"
//   struct T_SgInfo;
extern "C"
{
#include "sglite/sglite.h"
}
class SpaceGroup;
//######################################################################
//  AsymmetricUnit.
/**
*
* \brief The basic description of spacegroup asymmetric unit.
*
* Only xmin,xmax,ymin,ymax and zmin,zmax are recorded, thus resulting
* in a parallelepipedic unit with one (0,0,0) corner.
* It is not really 'asymmetric' since more than the crystallographic asymmetric
* unit can be included in it.
*
* \todo Currently the initialization of the asymmetric unit is done
* numerically, slowly. A faster algorithm should be used (using
* dichotomy), or we could switch to using a table of asymmetric units.
*/
//######################################################################

class AsymmetricUnit
{
   public:
      /// Default Constructor.
      AsymmetricUnit();
      /// Constructor, for a given spacegroup.
      AsymmetricUnit(const SpaceGroup &spg);
      ~AsymmetricUnit();
      /// Assign a SpaceGroup and generate the corrsponding Xmax, Ymax, ZMax.
      void SetSpaceGroup(const SpaceGroup &spg);
      /// Test if (x,y,z) is in the asymmetric unit
      bool IsInAsymmetricUnit(const REAL x, const REAL y, const REAL z)const;
      REAL Xmin() const;
      REAL Xmax() const;
      REAL Ymin() const;
      REAL Ymax() const;
      REAL Zmin() const;
      REAL Zmax() const;
   protected:
   private:
      REAL mXmin,mXmax,mYmin,mYmax,mZmin,mZmax;
};

//######################################################################
//                SpaceGroup
/**
* \brief The crystallographic space group, and the cell choice.
*
* This class includes functions to get basic information about
* the symmetries, as well as getting all symmetrics for a given
* position in a unit cell.
*
* This class included a pointer to a function calculating the "geometrical
* structure factor" (ie the sum of sin() and cos() for all symetrics, as
* could be found in the old version of the (red) International Tables),
* which was used to speed up computation of structure factors 
* by using pre-factorised formulas.
* This is not used anymore, since methods can be used to speed up computations.
*
* This class uses R. Grosse-Kunstleve 'SgLite' package,
* which is part of the Pymol package : http://pymol.sourceforge.net/
*
*\warning: the interface of the class will somewhat change when switching 
* from sgLite to cctbx (http://cctbx.sourceforge.net). Particularly
* functions Spacegroup::GetSgOps() and Spacegroup::GetHM_as_Hall() will
* be removed.
*/
//######################################################################

   
class SpaceGroup
{
   public:
		/// Default Constructor (initializes in P1)
		///
		/// You can use later SpaceGroup::ChangeSpaceGroup() to set the spacegroup.
      SpaceGroup();
      /** \brief Constructor with a specified spacegroup symbol or number
      *
      *  \param spgId The space group identifier, either an Hermann-Maugin,
      * or Hall, or Schonflies symbol.
      */
      SpaceGroup(const string &spgId);
      /// Destructor
      ~SpaceGroup();
      /// Change the Spacegroup
      void ChangeSpaceGroup(const string &spgId);
      /// Get the name of this spacegroup (its name, as supplied initially by
		/// the calling program or user)
      const string& GetName()const;
      /// Test if a given scatterer at (x,y,z) is in the asymmetric unit.
      bool IsInAsymmetricUnit(const REAL x, const REAL y, const REAL z) const;
      /// Move (x,y,z) coordinates to their equivalent in the asym unit
      /// \warning Not implemented yet.
		/// \todo SpaceGroup::IsInAsymmetricUnit()
      void ChangeToAsymmetricUnit(REAL x, REAL y, REAL z) const;//:TODO:
      /// Get the AsymmetricUnit for this spacegroup
      const AsymmetricUnit& GetAsymUnit() const;
      
      /// Id number of the spacegroup
      int GetSpaceGroupNumber()const;
      
      /// Is the crystal centrosymmetric ?
      bool IsCentrosymmetric()const;
      
      /** \brief Number of translation vectors
      * (1 for 'P' cells, 2 for 'I', 4 for 'F',etc..)
      *
      */
      int GetNbTranslationVectors()const;
      
      /** Return all Translation Vectors, as a 3 columns-array
      *
      * The first vector is always [0,0,0]
      *  \return 
      *    \f$ \left[ \begin {array}{ccc}  0 & 0 & 0 \end{array} \right] \f$
      * for a 'P' Cell,
      *    \f$ \left[ \begin {array}{ccc}  0 & 0 & 0 \\
      *                \frac{1}{2} & \frac{1}{2} & \frac{1}{2} \\ \end{array} \right] \f$
      * for a 'I' cell, and 
      *    \f$ \left[ \begin {array}{ccc}  0 & 0 & 0 \\
      *                       \frac{1}{2} & \frac{1}{2} & 0 \\
      *                       \frac{1}{2} & 0 & \frac{1}{2} \\
      *                       0 & \frac{1}{2} & \frac{1}{2} \\ \end{array} \right] \f$
      *for a 'F' cell,etc...
      */
      CrystMatrix_REAL GetTranslationVectors()const;
      
      /** \brief Get all equivalent positions of a (xyz) position
      *
      *  \param x,y,z fractional coordinates of the position
      *  \param  noCenter if set to 'false' (the default), then the center of
      * symmetry (if any) is used to generate ALL positions. If 'true', then
      * only one half of equivalent positions are generated. This has 
      * no influence if the group is not centrosymmetric. (\b note Not generating
      * symmetrical positions from center of symmetry is useful to speed up computation
      * of structure factor, but is a bit tricky if the inversion is not at the origin.
      * This is taken into account)
      *  \param  noTransl if set to 'false' (the default), then translation are
      * taken into account to generate all atom positions. This affect
      * only body or face(s)-centered spacegroups.
      *  \param  noIdentical if set to true, then atom in special positions
      * will only return the distinct atomic positions. Currently two atoms are considered
      * distinct if the difference for all of their fractionnal coordinates is less than 1e-5
      *  \return a 3-column (x,y,z) matrix with as many rows as symmetric atoms
      *  \warning 'special' positions are not taken into account. (ie an
      * atom in special position will return duplicate entries. This may be
      * corrected automatically later.) You can use the 'noIdentical' option for that,
      */
      CrystMatrix_REAL GetAllSymmetrics(const REAL x, const REAL y, const REAL z,
                                const bool noCenter=false,const bool noTransl=false,
                                const bool noIdentical=false) const;
      
      /** \brief Return the number of equivalent positions in the spacegroup,
      *ie the multilicity of the general position.
      *
      *  \param noCenter if 'true', do not take into account the center of symmetry
      *  \param noTransl if 'true', do not take into account translations
      */
      int GetNbSymmetrics(const bool noCenter=false,const bool noTransl=false)const;
      
      /// Prints a description of the spacegroup (symbol, properties).
		///
		/// \todo 
      void Print()const;
      /// Is centrosymmetric ?
      bool HasInversionCenter()const;
      /// Is the center of symmetry at the origin ?
      bool IsInversionCenterAtOrigin()const;
      /// Get the SgOps structure. This will be removed when switching to cctbx.
      const T_SgOps& GetSgOps()const;
      /// Get the SpaceGroup Clock (corresponding to the time of the
		/// initialization of the SpaceGroup)
      const RefinableObjClock& GetClockSpaceGroup() const;
      /// Access to the HM_As_Hall structure. This will be removed when switching to cctbx.
      const T_HM_as_Hall& GetHM_as_Hall()const;
      /// Which is the unique axis (for monoclinic space groups )
      unsigned int GetUniqueAxis()const;
   protected:
   private:
      /** \brief Init the spaceGroup object from its name
      *
      *Initialize the SgOps & HM_as_Hall structures (SgLite),
      *and set the function pointers to the functions used to
      *compute the geometrical structure factors.
      */
      void InitSpaceGroup(const string &spgId);
      
      /// Spacegroup's name ( 'I422', 'D2^8','230')
		/// Maybe we should only store the Hermann-Mauguin symbol, rather than storing
		/// the string which was initially given by the user/program for the initialization.
      string mId;
      
      /** \brief  SgOps structure for this spacegroup. (Symmetry operations)
      *
      * See sglite subdirectory for more information.
      * This is (c) R. Gross-Kunstleve, part of PyMol software
      * http://pymol.sourceforge.net/
      */
      T_SgOps mSgOps;
      
      /** \brief Is spacegroup centrosymmetric ?
      *
      */
      bool mHasInversionCenter;
      /** \brief Is center of symmetry at the origin ?
      *
      */
      bool mIsInversionCenterAtOrigin;
      
      /** \brief  SgOps structure for this spacegroup. (Symmetry operations)
      *
      * See sglite subdirectory for more information.
      * This is (c) R. Gross-Kunstleve, part of PyMol software
      * http://pymol.sourceforge.net/
      */
      T_HM_as_Hall mHM_as_Hall;


      /// The spacegroup asymmetric unit
      AsymmetricUnit mAsymmetricUnit;
      
      /** Use geometrical structure factor ?
      *
		* \deprecated
      * when all atoms have an isotropic 
      * thermic factor, use a sophisticated formula to compute the structure 
      * factor for a given independent atom, instead of generating all symmetric 
      * positions and summing all contributions. This NOT USED CURRENTLY.
      * Given the speed improvement 
      * obtained by using tabulated sine and cosines, the use of geometrical 
      * structure factors has become useless and will be removed from the class.
      */
      const static bool mUseGeomStructFactor=false;
      
      /// The Spacegroup clock
      RefinableObjClock mClock;
      /// Unique axis number (0=a,1=b,2=c)
      unsigned int mUniqueAxisId;
};

}//namespace
#endif //_OBJCRYST_SPACEGROUP_H_
