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
//######################################################################
//  AsymmetricUnit.
/**
*
* \brief The basic description of spacegroup asymmetric unit.
*
* Only xmin,xmax,ymin,ymax and zmin,zmax are recorded, thus resulting
* in a parallelepipedic unit with one (0,0,0) corner.
* It is not really 'asymmetric' since more than the crystallographi asymmetric
* unit can be included in it.
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
      bool IsInAsymmetricUnit(const double x, const double y, const double z)const;
      double Xmin() const;
      double Xmax() const;
      double Ymin() const;
      double Ymax() const;
      double Zmin() const;
      double Zmax() const;
   protected:
   private:
      double mXmin,mXmax,mYmin,mYmax,mZmin,mZmax;
};

//######################################################################
//                SpaceGroup
/**
* \brief The crystallographic space group, and the cell choice.
*
* This class includes functions to compute the geometrical structure
*factor (real and imaginary part),using the formulas given in
*the International Tables for X-Ray Crystallography (1969).
*
* This class is based on R. Grosse-Kunstleve 'SgLite' package,
* which is part of Pymol : http://pymol.sourceforge.net/
*/
//######################################################################

   
class SpaceGroup
{
   public:
      SpaceGroup();
      /** \brief Spacegroup constructor
      *
      *  \param spgId The space group identifier, either an Hermann-Maugin,
      * or Hall, or Schonflies symbol.
      */
      SpaceGroup(const string &spgId);
      ///Spacegroup destructor
      ~SpaceGroup();
      ///Change the Spacegroup
      void ChangeSpaceGroup(const string &spgId);
      /// Get the name of this spacegroup (well... its conventionnal name, as supplied)
      const string& GetName()const;
      ///Test if a given scatterer at (x,y,z) is in the asymmetric unit.
      ///This is not really implemented yet
      bool IsInAsymmetricUnit(const double x, const double y, const double z) const;
      ///Move (x,y,z) coordinates to their equivalent in the asym unit
      ///Not implemented yet.
      void ChangeToAsymmetricUnit(double x, double y, double z) const;//:TODO:
      //Returns the AsymmetricUnit of this spacegroup
      const AsymmetricUnit& GetAsymUnit() const;
      
      /// Id number of the spacegroup
      int GetSpaceGroupNumber()const;
      
      /// Is the crystal centrosymmetric ?
      bool IsCentrosymmetric()const;
      
      /** \brief Number of translation vectors
      *(1 for 'P' cells, 2 for 'I', 4 for 'F',etc..)
      *
      *The first vector is always [0,0,0]
      */
      int GetNbTranslationVectors()const;
      
      /** Return all Translation Vactors, as a 3 columns-array
      *
      *  \return 
      *    \f$ \left[ \begin {array}{ccc}  0 & 0 & 0 \end{array} \right] \f$
      *for a 'P' Cell,
      *    \f$ \left[ \begin {array}{ccc}  0 & 0 & 0 \\
      *                \frac{1}{2} & \frac{1}{2} & \frac{1}{2} \\ \end{array} \right] \f$
      *for a 'I' cell, and 
      *    \f$ \left[ \begin {array}{ccc}  0 & 0 & 0 \\
      *                       \frac{1}{2} & \frac{1}{2} & 0 \\
      *                       \frac{1}{2} & 0 & \frac{1}{2} \\
      *                       0 & \frac{1}{2} & \frac{1}{2} \\ \end{array} \right] \f$
      *for a 'F' cell,etc...
      */
      CrystMatrix_double GetTranslationVectors()const;
      
      /** \brief Get all equivalent position of a scatterer
      *
      *  \param x,y,z doubleing-point fractional coordinates of the scatterer
      *  \param  noCenter if set to 'false' (the default), then the center of
      * symmetry (if any) is used to generate ALL positions. IF 'true', then
      * only one half of equivalent positions are generated. This has 
      * no influence if the group is not centrosymmetric. (\b note Not generating
      * symmetrical positions from center of symmetry is useful to speed up computation
      * of structure factor, but is a bit tricky if the inversion is not at the origin.
      * This is taken into account in SpaceGroup::GeomStructFactor)
      *  \param  noTransl if set to 'false' (the default), then translation are
      *taken into account to generate all atom positions. This affect
      *only body or face(s)-centered spacegroups.
      *  \param  noIdentical if set to true, then atom in special positions
      * will only return the distinct atomic positions. Currently two atoms are considered
      * distinct if the difference for all of their fractionnal coordinates is less than 1e-5
      *  \return a 3-column (x,y,z) matrix with as many rows as symetric atoms
      *  \warning 'special' positions are not taken into account. (ie an
      *atom in special position will return duplicate entries. This may be
      *corrected automatically later.) Use the 'noIdentical' option for that.
      */
      CrystMatrix_double GetAllSymetrics(const double x, const double y, const double z,
                                const bool noCenter=false,const bool noTransl=false,
                                const bool noIdentical=false) const;
      
      /** \brief Return the number of equivalent positions in the spacegroup,
      *ie the multilicity of the general position.
      *
      *  \param noCenter if 'true', do not take into account the center of symmetry
      *  \param noTransl if 'true', do not take into account translations
      */
      int GetNbSymetrics(const bool noCenter=false,const bool noTransl=false)const;
      
      ///Prints a short description of the spacegroup (one line). 
      void Print()const;
      /// Is centrosymmetric ?
      bool HasInversionCenter()const;
      /// Is the center of symmetry at the origin
      bool IsInversionCenterAtOrigin()const;
      /// Get the SgOps structure. This will be removed to something nicer.
      const T_SgOps& GetSgOps()const;
      /// Get the SpaceGroup Clock (corresponding to the initialization of the SpaceGroup)
      const RefinableObjClock& GetClockSpaceGroup() const;
      /// Access to the HM_As_Hall structure
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
      
      /// Spacegroup's name ( 'I422', 'D2^8')
      string mId;
      
      /** \brief  SgOps structure for this spacegroup. (Symmetry operations)
      *
      *See sglite subdirectory for more information.
      *This is (c) R. Gross-Kunstleve, part of PyMol software
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
      *See sglite subdirectory for more information.
      *This is (c) R. Gross-Kunstleve, part of PyMol software
      * http://pymol.sourceforge.net/
      */
      T_HM_as_Hall mHM_as_Hall;


      /// The spacegroup asymmetric unit
      AsymmetricUnit mAsymmetricUnit;
      
      /**Use geometrical structure factor ? OBSOLETE ?
      *
      *when all atoms have an isotropic 
      *thermic factor, use a sophisticated formula to compute the structure 
      *factor for a given independent atom, instead of generating all symetric 
      *positions and summing all contributions.
      *
      * \b warning about Geometrical structure factors : given the speed improvement 
      *obtained by using tabulated sine and cosines, the use of geometrical 
      *structure factors may become useless and be removed from the class.
      */
      const static bool mUseGeomStructFactor=false;
      
      ///The Spacegroup clock
      RefinableObjClock mClock;
      /// Unique axis number (0=a,1=b,2=c)
      unsigned int mUniqueAxisId;
};

}//namespace
#endif //_OBJCRYST_SPACEGROUP_H_
