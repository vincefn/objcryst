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
#ifndef _OBJCRYST_ZSCATTERER_H_
#define _OBJCRYST_ZSCATTERER_H_

#include "CrystVector/CrystVector.h"

#include "ObjCryst/General.h"

#include "ObjCryst/ScatteringPower.h"
#include "ObjCryst/Scatterer.h"
#include "RefinableObj/GlobalOptimObj.h"

#include <string>

namespace ObjCryst
{
class ZScatterer;
//######################################################################
//
//      GLOBAL SCATTERING POWER
/**
* \brief Global Scattering Power. Used to approximate the scattering
* power of a multi-atom ZScatterer (polyhedron,...) to an isotropic
* scattering power.
*
* The scattering power will be (\b slowly) approximated to the average
* scattering power of this ZScatterer take in all orientations.
* The approximated scattering factor will include the temperature factor
* and resonant scattering factors contributions.
*
* This is only used in ZScatterer, and maybe obsolete (?).
*/
//######################################################################

class GlobalScatteringPower:virtual public ScatteringPower
{
   public:
      GlobalScatteringPower();
      GlobalScatteringPower(const ZScatterer &scatt);
      GlobalScatteringPower(const GlobalScatteringPower& old);
      ~GlobalScatteringPower();
      /// Re-initialize parameters (after using the default constructor).
      void Init(const ZScatterer &scatt);
      virtual CrystVector_REAL GetScatteringFactor(const ScatteringData &data,
                                                     const int spgSymPosIndex=0) const;
      virtual REAL GetForwardScatteringFactor(const RadiationType) const;
      virtual CrystVector_REAL GetTemperatureFactor(const ScatteringData &data,
                                                          const int spgSymPosIndex=0) const;
      virtual CrystMatrix_REAL GetResonantScattFactReal(const ScatteringData &data,
                                                          const int spgSymPosIndex=0) const;
      virtual CrystMatrix_REAL GetResonantScattFactImag(const ScatteringData &data,
                                                          const int spgSymPosIndex=0) const;
      virtual REAL GetRadius()const;
   protected:
      virtual void InitRefParList();
      /// a copy of the ZScatterer associated to this object
      ZScatterer *mpZScatterer;
   private:
};

//######################################################################
///  Class for individual atoms in a ZScatterer Object. This class
///  is \e purely \e internal to ZScatterer, so should not be used
///  for any other purpose...
//######################################################################

class ZAtom
{
   public:
      ZAtom(ZScatterer &scatt,const ScatteringPower *pow,
            const long atomBond=0, const REAL bondLength=1,
            const long atomAngle=0, const REAL bondAngle=M_PI,
            const long atomDihedral=0, const REAL dihedralAngle=M_PI,
            const REAL popu=1., const string &name="");
      ~ZAtom();
      const string& GetClassName()const;
      const string& GetName()const;
      void SetName(const string&);
      /// Get the ZScatterer associated to this ZAtom
      const ZScatterer& GetZScatterer()const;
      /// Get the ZScatterer associated to this ZAtom
      ZScatterer& GetZScatterer();
      
      /// Index of the 1st atom used to define the atom in the Z-Matrix (the one from
      /// which the bondlength is calculated)
      long GetZBondAtom()const;
      /// Index of the 2nd atom used to define the atom in the Z-Matrix (the one from
      /// which the angle is calculated)
      long GetZAngleAtom()const;
      /// Index of the 3rd atom used to define the atom in the Z-Matrix (the one from
      /// which the dihedral angle is calculated)
      long GetZDihedralAngleAtom()const;
      
      ///Const access to bondlength parameter.
      REAL GetZBondLength()const;
      ///Const access to the angle parameter.
      REAL GetZAngle()const;
      ///Const access to the dihedral angle parameter.
      REAL GetZDihedralAngle()const;
      ///Const access to the ocupancy parameter.
      REAL GetOccupancy()const;
      ///ScatteringPower for this atom.
      const ScatteringPower* GetScatteringPower()const;
      
      ///Access to bondlength parameter.
      void SetZBondLength(const REAL);
      ///Access to the angle parameter.
      void SetZAngle(const REAL);
      ///Access to the dihedral angle parameter.
      void SetZDihedralAngle(const REAL);
      ///Access to the dihedral angle parameter.
      void SetOccupancy(const REAL);
      ///Set the ScatteringPower.
      void SetScatteringPower(const ScatteringPower*);
      void XMLOutput(ostream &os,int indent=0)const;
      void XMLInput(istream &is,const XMLCrystTag &tag);
   private:
      /// The ScatteringPower corresponding to this atom.
      const ScatteringPower *mpScattPow;
      /// The index (in the ZScatterer) of the atoms which are used to define the
      /// position of this atom.
      long mAtomBond,mAtomAngle,mAtomDihed;
      /// Bond length, angle and dihedral angle.
      REAL mBondLength,mAngle,mDihed,mOccupancy;
      /// Name for this atom
      string mName;
      /// the ZScatterer in which this atom is included.
      ZScatterer *mpScatt;
      
      friend class ZScatterer; //So that RefinablePar can be declared in ZScatterer
      
   #ifdef __WX__CRYST__
   public:
      WXCrystObjBasic* WXCreate(wxWindow *parent);
      WXCrystObjBasic* WXGet();
      void WXDelete();
      void WXNotifyDelete();
   private:
      WXCrystObjBasic *mpWXCrystObj;
      friend class WXZAtom;
   #endif
};
//######################################################################
/**  \brief Class to minimize conformation changes for random moves. Very
* experimental !!!
* \internal Only used within ZScatterer.
*/
//######################################################################
class ZMoveMinimizer:public RefinableObj
{
   public:
      ZMoveMinimizer(ZScatterer &scatt);
      ~ZMoveMinimizer();
      virtual unsigned int GetNbCostFunction()const;
      virtual const string& GetCostFunctionName(const unsigned int)const;
      virtual const string& GetCostFunctionDescription(const unsigned int)const;
      virtual REAL GetCostFunctionValue(const unsigned int);
      void RecordConformation();
      void SetZAtomWeight(const CrystVector_REAL weight);
      void MinimizeChange(long nbTrial=10000);
   private:
      ZScatterer *mpZScatt;
      MonteCarloObj mOptimObj;
      CrystVector_REAL mXCoord0,mYCoord0,mZCoord0;
      CrystVector_REAL mAtomWeight;
};

//######################################################################
///  ZScatterer: the basic type of complex scatterers, where atom positions
///  are defined using a standard "Z-Matrix" description. This is used
/// to describe inorganic polyhedras, as well as molecules.
//######################################################################

class ZScatterer: public Scatterer
{
   public:
      /**   \brief ZScatterer constructor
      *
      *  \param name: the name of the scatterer
      *  \param cryst: the crystal in which the scatterer is (needed to convert
      *  from cartesian to fractionnal coordinates).
      *  \param x,y,z: fractionnal coordinates of the scatterer
      *  \param phi,chi: angles defining the orientation of the scatterer
      */
      ZScatterer(const string &name,Crystal &cryst, 
                 const REAL x=0.,const REAL y=0.,const REAL z=0.,
                 const REAL phi=0.,const REAL chi=0., const REAL psi=0.);
      /** \brief Copy constructor
      *
      */
      ZScatterer(const ZScatterer &old);
     ~ZScatterer();
      /// \internal so-called Virtual copy constructor, needed to make copies
      /// of arrays of Scatterers
      virtual ZScatterer* CreateCopy() const;
      virtual const string& GetClassName() const;
      /// Add an atom to the Zscatterer. If &ScatteringPower=0, then it is a 'dummy'
      /// atom and will be ignored for any scattering analysis. The 'name' supplied may
      /// not be respected, and can be replaced by 'ZScatterer_name'+'AtomNum'+'ScattPowName'
      void AddAtom(const string &name,const ScatteringPower *pow,
                   const long atomBond, const REAL bondLength,
                   const long atomAngle, const REAL bondAngle,
                   const long atomDihedral, const REAL dihedralAngle,
                   const REAL popu=1.);
      
      virtual int GetNbComponent() const;
      virtual const ScatteringComponentList& GetScatteringComponentList() const;
      virtual string GetComponentName(const int i) const;
      
      ///Print a single line of information about this scatterer
      void Print() const;
            
      ///Access to phi parameter (overall orientation of the scatterer)
      REAL GetPhi()const;
      ///Access to chi parameter (overall orientation of the scatterer)
      REAL GetChi()const;
      ///Access to psi parameter (overall orientation of the scatterer)
      REAL GetPsi()const;
      ///Access to phi parameter (overall orientation of the scatterer)
      void SetPhi(const REAL);
      ///Access to chi parameter (overall orientation of the scatterer)
      void SetChi(const REAL);
      ///Access to psi parameter (overall orientation of the scatterer)
      void SetPsi(const REAL);

      /// Index of the 1st atom used to define the i-th atom in the Z-Matrix (the one from
      /// which the bondlength is calculated)
      long GetZBondAtom(const int i)const;
      /// Index of the 2nd atom used to define the i-th atom in the Z-Matrix (the one from
      /// which the angle is calculated)
      long GetZAngleAtom(const int i)const;
      /// Index of the 3rd atom used to define the i-th atom in the Z-Matrix (the one from
      /// which the dihedral angle is calculated)
      long GetZDihedralAngleAtom(const int i)const;
      
      ///Const access to bondlength parameter, for the i-th row in the Z-Matrix.
      REAL GetZBondLength(const int i)const;
      ///Const access to the angle parameter, for the i-th row in the Z-Matrix.
      REAL GetZAngle(const int i)const;
      ///Const access to the dihedral angle parameter, for the i-th row in the Z-Matrix.
      REAL GetZDihedralAngle(const int i)const;
      
      ///Access to bondlength parameter, for the i-th row in the Z-Matrix.
      void SetZBondLength(const int i,const REAL);
      ///Access to the angle parameter, for the i-th row in the Z-Matrix.
      void SetZAngle(const int i,const REAL);
      ///Access to the dihedral angle parameter, for the i-th row in the Z-Matrix.
      void SetZDihedralAngle(const int i,const REAL);

      ///Access to the registry of ZAtoms
      const ObjRegistry<ZAtom>& GetZAtomRegistry()const;
      /// \warning Not implemented for ZScatterer
      virtual ostream& POVRayDescription(ostream &os,
                                         const bool onlyIndependentAtoms=false)const;

      virtual void GLInitDisplayList(const bool onlyIndependentAtoms=false,
                                     const REAL xMin=-.1,const REAL xMax=1.1,
                                     const REAL yMin=-.1,const REAL yMax=1.1,
                                     const REAL zMin=-.1,const REAL zMax=1.1,
                                     const bool displayEnantiomer=false)const;
      /** \brief use a Global scattering power for this scatterer ?
      *
      * If true, then the overall scattering power of this ZScatterer will be 
      * approximated to an isotropic scattering power computed for this scatterer.
      * Of course, only use this if the "isotropic" approximation is reasonable for
      * this scatterer (typically true for 'large' polyhedra). See GlobalScatteringPower.
      *
      * \warning EXPERIMENTAL
      */
      virtual void SetUseGlobalScatteringPower(const bool useIt);
      virtual void XMLOutput(ostream &os,int indent=0)const;
      virtual void XMLInput(istream &is,const XMLCrystTag &tag);
      //virtual void XMLInputOld(istream &is,const IOCrystTag &tag);
      virtual void GetGeneGroup(const RefinableObj &obj, 
                                CrystVector_uint & groupIndex,
                                unsigned int &firstGroup) const;
      virtual void GlobalOptRandomMove(const REAL mutationAmplitude);
      /// Get the list of all ZAtom cartesian x coordinates.
      const CrystVector_REAL& GetXCoord() const;
      /// Get the list of all ZAtom cartesian x coordinates.
      const CrystVector_REAL& GetYCoord() const;
      /// Get the list of all ZAtom cartesian x coordinates.
      const CrystVector_REAL& GetZCoord() const;
      virtual void EndOptimization();
      /** Import "Fenske-Hall" ZMatrix file (fhz in the babel
      * program http://www.eyesopen.com/babel.html\
      * example: use "./babel -ipdb foo.pdb -ofhz foo.fhz -d",
      * to convert a pdb file to a Z-Matrix file (the -d removes
      * hydrogen atoms)
      *
      * \warning: this should be called before any atom has been
      * added (if there are already atoms, they should be removed
      * but this has not been tested...)
      *
      * \note: this will search in the Crystal associated with this
      * ZScatterer the relevant ScatteringPowerAtom, which should have
      * the name of the corresponding symbol (eg 'H', 'C',...)
      * if these are not found then they will be added to the Crystal
      * with a default isotropic B-factor equal to 1.
      * \note: this also sets relative limits of +/-.03 Angstroems
      * for all bond distances, and +/-3.6 degress for bond and
      * dihedral angles.
      * \todo: identify which dihedral angles should \e not be limited,
      * by analysing a coordination table.
      */
      void ImportFenskeHallZMatrix(istream &is);
      /** Export to Fenske-Hall ZMatrix file
      * 
      * \todo USe more strict formatting than space-delimited.
      */
      void ExportFenskeHallZMatrix(ostream &os);
   protected:
      /** Update the atom coordinates (in real units, in Angstroems).
      *
      * This takes into account the translation and global 
      * rotation of the scatterer (ie this does not generate 'internal
      * coordinates).
      */
      void UpdateCoordinates() const;
      /** Update the scattering component list, ie compute all atom
      * positions from the bonds/angles/dihedral angles, and convert
      * the coordinates to fractionnal coordinates of the Crystal.
      *
      * 
      */
      void UpdateScattCompList() const;
      /** For 3D display of the structure, bonds, triangular and quadric
      * faces can be displayed. This matrix determines what is drawn.
      * This is a 5-column matrix. The first column indicates the type of
      * drawing (0: : nothing, 1: display the atom (a sphere),  2: bond, 
      * 3: triangular face, 4: quadric face)
      * the other columns indicate the index of the atoms involved in the 
      * drawing (2 atoms for a bond, 3 for....)
      *
      * If the matrix is empty only the individual atoms are displayed.
      *
      * \todo This is still experimental. This is only used for the display
      * of ZPolyhedron, and should be more developped (and it should also
      * be saved in XML files !)
      */
      CrystMatrix_long m3DDisplayIndex;
      /// The list of scattering components.
      mutable ScatteringComponentList mScattCompList;
      ///Total number of atoms in the structure
      long mNbAtom;
   private:
      ///Prepare refinable parameters for the scatterer object
      virtual void InitRefParList();
            
      ///Number of "dummy" atoms in the structure
      long mNbDummyAtom;
      
      /// Index of atoms in the ScatteringComponentList. We need this list because
      /// Dummy atoms are not included in it. So this is an integer array with
      /// [ 0 1 2 4 5 6].. where the missing '3' marks a Dummy atom. Thus
      /// to get the name of component 'i' in the component list, you must
      /// take the mComponentIndex(i) atom in the ZAtom Registry.
      CrystVector_int mComponentIndex;
      
      /** \brief Angles giving the orientation of the ZScatterer (stored in radian)
      *
      * The position of any atom can be transformed from internal coordinates (orthonormal
      * coordinates derived from the ZMatrix, with first atom at (0,0,0), second
      * atom at (x,0,0), third atom at (x,y,0),...) to orthonormal coordinates in
      * the crystal reference frame (ie with orientation of the ZScatterer) using :
      * \f[ \left[ \begin{array}{c} x(i) \\ y(i) \\ z(i) \end{array} \right]_{orthonormal}
              = \left[ \begin{array}{ccc} \cos(\chi) & 0 & -\sin(\chi) \\
                                           0 & 1 & 0 \\
                                           \sin(\chi) & 0 & \cos(\chi) \end{array} \right]
         \times \left[ \begin{array}{ccc} \cos(\phi) & -\sin(\phi) & 0 \\
                                           \sin(\phi) & \cos(\phi) & 0 \\
                                           0 & 0 & 1 \end{array} \right]
         \times \left[ \begin{array}{ccc} 1 & 0 & 0 \\
                                           0 & \cos(\psi) & -\sin(\psi) \\
                                           0 & \sin(\psi) & \cos(\psi) \end{array} \right]
         \times \left[ \begin{array}{c} x_0(i) \\ y_0(i) \\ z_0(i) \end{array} \right]
      * \f]
      *, where x0(i), y0(i) and z0(i) describe the position for atom (i) in
      * internal coordinates, and x(i), y(i), z(i) are coordinates of the rotated ZScatterer.
      *
      *
      * The rotation is performed around a 'pivot' atom (see ZScatterer::mPivotAtom)
      */
      REAL mPhi,mChi,mPsi;
      
      /// Registry for ZAtoms in this Scatterer.
      ObjRegistry<ZAtom> mZAtomRegistry;
      /// Index of the atom used as a pivot (the scatterer is rotated around this atom).
      /// This should more or less be at the center of the Scatterer.
      long mCenterAtomIndex;
      
      /// Rotation matrix for the orientation of the scatterer
      mutable CrystMatrix_REAL mPhiChiPsiMatrix;
      
      /// Does the ZScatterer use a global scattering power ?
      /// \warning EXPERIMENTAL.
      bool mUseGlobalScattPow;
      
      /// the global scattering power used, if mUseGlobalScattPow=true
      /// \warning EXPERIMENTAL.
      GlobalScatteringPower* mpGlobalScattPow;
      
      /// Storage for Cartesian coordinates. The (0,0,0) is on the central atom. This
      /// includes Dummy atoms.
      mutable CrystVector_REAL mXCoord,mYCoord,mZCoord;
      /// Last time the cartesian coordinates were computed
      mutable RefinableObjClock mClockCoord;
      ZMoveMinimizer *mpZMoveMinimizer;
   #ifdef __WX__CRYST__
   public:
      virtual WXCrystObjBasic* WXCreate(wxWindow*);
   friend class WXZScatterer;
   #endif
};

//######################################################################
//
//    Different types of regular polyhedra
//
//######################################################################

enum RegularPolyhedraType { TETRAHEDRON, OCTAHEDRON, SQUARE_PLANE, CUBE, ANTIPRISM_TETRAGONAL,
                            PRISM_TETRAGONAL_MONOCAP, PRISM_TETRAGONAL_DICAP,
                            PRISM_TRIGONAL,PRISM_TRIGONAL_TRICAPPED,
                            ICOSAHEDRON, TRIANGLE_PLANE};

//######################################################################
///  \class ZPolyhedron include.h ObjCryst/ZScatterer.h
///
///  ZPolyhedron: a Scatterer to describe polyhedras such as octahedron,
/// tetrahedron, square plane, etc... These are ZScatterer objects, so that
/// even if they are initialized with constraints, these can be removed
/// to make any configuration.
//######################################################################
class ZPolyhedron: public ZScatterer
{
   public:
      /**   \brief ZPolyhedron constructor
      *  \param type : OCTAHEDRON, TETRAHEDRON,...
      *  \param cryst : a crystal is necessary to transform the orthornormal
      * coordinates (in which the relative atom positions are computed) to fractional ones.
      *  \param x,y,z : \e fractional coordinates of the center of the polyhedra
      *  \param name : name of the Polyhedra ('WO6','TaSe4_a', 'Tetra01'...).
      * The name can have \e any format but spaces should be avoided, since it
      * will generate problems when reading the names from a file. And there must \b not
      * be identical name for two scatterers in a given crystal object.
      *  \param  periphAtomPow,centralAtomPow: the ScatteringPower corresponding to the
      * central and peripheral atoms, respectively.
      *  \param centralPeriphDist: the distance, in angstroems, from the central
      * to the peripheral atoms
      * \param ligandPopu : the relative population of ligand atoms, equal to 1/n if 
      * ligand atoms are shared by n polyhedra. If you are using the 'Dynamical Polpulation
      * Correction', then keep this parameter to 1.
      * \param phi,chi,psi: initial angles for this polyhedron
      */
      ZPolyhedron( const RegularPolyhedraType type, Crystal &cryst,
            const REAL x, const REAL y, const REAL z,
            const string &name, const ScatteringPower *centralAtomPow,
            const ScatteringPower *periphAtomPow,const REAL centralPeriphDist,
            const REAL ligandPopu=1,
            const REAL phi=0., const REAL chi=0., const REAL psi=0.);
      /// Copy Constructor
      ZPolyhedron(const ZPolyhedron&);
      /// \internal so-called Virtual copy constructor, needed to make copies
      ///of arrays of Scatterers
      virtual ZPolyhedron* CreateCopy() const;
      
      /* \brief Copy constructor
      *
      ZPolyhedron(const ZPolyhedron &old);
      */
   protected:
   private:
      //Prepare refinable parameters for the scatterer object
      //virtual void InitRefParList();
      
      ///
      RegularPolyhedraType mPolyhedraType;
};

}//namespace
#include "ObjCryst/Crystal.h"

#endif //_OBJCRYST_ZSCATTERER_H_
