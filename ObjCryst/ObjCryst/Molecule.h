/*  ObjCryst++ Object-Oriented Crystallographic Library
    (c) 2000-2002 Vincent Favre-Nicolin vincefn@users.sourceforge.net
        2000-2001 University of Geneva (Switzerland)

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation;  version 2 of the License.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/*   Molecule.h
*  header file for the Molecule scatterer
*
*/
#ifndef _OBJCRYST_MOLECULE_H_
#define _OBJCRYST_MOLECULE_H_

#include <string>
#include <vector>
#include <map>
#include <set>
#include <list>

#include "ObjCryst/General.h"
#include "ObjCryst/ScatteringPower.h"
#include "ObjCryst/Scatterer.h"


namespace ObjCryst
{
class MolAtom;
class MolBond;
class Molecule;

/** MolAtom : atom inside a Molecule
*
* This keeps coordinates, recorded in a cartesian frame (in Angstroem),
* the associated scattering power 
* and it also keeps in a list of all bonds in which this atom is involved.
*
* \note maybe it's not a great idea to keep a reference of bonds for this
* atom in here
*/
class MolAtom
{
   public:
      /** Constructor for a MolAtom
      *
      */
      MolAtom(const REAL x, const REAL y, const REAL z,
              const ScatteringPower *pPow, const string &name,
              Molecule &parent);
      /** Destructor
      *
      * Tells the parent Molecule and all Bond that it is being destroyed.
      */
      virtual ~MolAtom();
      void SetName(const string &name);
      const string& GetName()const;
      string& GetName();
      const Molecule& GetMolecule()const;
      Molecule& GetMolecule();
      const REAL& X()const;
      const REAL& Y()const;
      const REAL& Z()const;
      REAL& X();
      REAL& Y();
      REAL& Z();
      REAL GetX()const;
      REAL GetY()const;
      REAL GetZ()const;
      REAL GetOccupancy()const;
      void SetX(const REAL);
      void SetY(const REAL);
      void SetZ(const REAL);
      void SetOccupancy(const REAL);
      /** Returns true if this is a dummy atom, i.e. without an associated scattering power.
      *
      * Dummy atoms can be used to mark positions, or for restraints.
      */
      bool IsDummy()const;
      const ScatteringPower& GetScatteringPower()const;
      void SetScatteringPower(const ScatteringPower&);
      virtual void XMLOutput(ostream &os,int indent=0)const;
      virtual void XMLInput(istream &is,const XMLCrystTag &tag);
   private:
      /// Name for this atom
      string mName;
      /* Get the atom at the other end of bond #i
      MolAtom & GetBondedAtom(unsigned int i);
      */
      /// Cartesian oordinates in the Molecule reference frame.
      REAL mX,mY,mZ;
      /// Occupancy
      REAL mOccupancy;
      /// ScatteringPower
      const ScatteringPower* mpScattPow;
      /// Parent Molecule
      Molecule *mpMol;
   #ifdef __WX__CRYST__
   public:
      WXCrystObjBasic *mpWXCrystObj;
      virtual WXCrystObjBasic* WXCreate(wxWindow*);
   #endif
};

/** Bond between two atoms, also a restraint on the associated bond length.
*
*/
class MolBond:public Restraint
{
   public:
      /** Constructor 
      *
      * Both atoms of the bond are told of the creation of the bond, so that
      * they can keep a list of bonds they are involved in.
      *
      * \param atom1, atom2: the atoms of the bond
      * \param length: the expected bond length
      * \param sigma,delta: depending on the calculated length, the log(likelihood) is equal to:
      * - within \f$[length_{0}-\delta;length_{0}+\delta]\f$:
      * \f$ -\log(likelihood)= \log\left(2\delta+\sqrt{2\pi\sigma^2}\right)\f$
      * - if \f$length > length_{0}+\delta\f$:
      * \f$ -\log(likelihood)= \log\left(2\delta+\sqrt{2\pi\sigma^2}\right)
      * + \left(\frac{(length-\delta)-length_{0}}{\sigma} \right)^2\f$
      * - if \f$length < length_{0}-\delta\f$:
      * \f$ -\log(likelihood)= \log\left(2\delta+\sqrt{2\pi\sigma^2}\right)
      * + \left(\frac{(length+\delta)-length_{0}}{\sigma} \right)^2\f$
      * 
      */
      MolBond(const MolAtom &atom1, const MolAtom &atom2,
              const REAL length, const REAL sigma, const REAL delta,
              Molecule &parent,const REAL bondOrder=1.);
      /** Destructor
      *
      * Notifies the atoms that the bond has disappeared.
      */
      virtual ~MolBond();
      const Molecule& GetMolecule()const;
      Molecule& GetMolecule();
      virtual void XMLOutput(ostream &os,int indent=0)const;
      virtual void XMLInput(istream &is,const XMLCrystTag &tag);
      virtual REAL GetLogLikelihood()const;
      const MolAtom& GetAtom1()const;
      const MolAtom& GetAtom2()const;
      void SetAtom1(const MolAtom &at1);
      void SetAtom2(const MolAtom &at2);
      REAL GetLength()const;
      REAL GetLength0()const;
      REAL GetLengthDelta()const;
      REAL GetLengthSigma()const;
      REAL GetBondOrder()const;
      REAL& Length0();
      REAL& LengthDelta();
      REAL& LengthSigma();
      REAL& BondOrder();
      void SetLength0(const REAL length);
      void SetLengthDelta(const REAL length);
      void SetLengthSigma(const REAL length);
      void SetBondOrder(const REAL length);
      bool IsFreeTorsion()const;
      void SetFreeTorsion(const bool isInRing);
  private:
      pair<const MolAtom*,const MolAtom*> mAtomPair;
      REAL mLength0,mDelta,mSigma;
      REAL mBondOrder;
      bool mIsFreeTorsion;
      /// Parent Molecule
      Molecule *mpMol;
   #ifdef __WX__CRYST__
   public:
      WXCrystObjBasic *mpWXCrystObj;
      virtual WXCrystObjBasic* WXCreate(wxWindow*);
   #endif
};

/** Bond angle restraint between 3 atoms.
*
* The atoms involved need not be bonded.
*
*
*/
class MolBondAngle:public Restraint
{
   public:
      /** Constructor 
      *
      */
      MolBondAngle(const MolAtom &atom1,const  MolAtom &atom2,const  MolAtom &atom3,
                   const REAL angle, const REAL sigma, const REAL delta,
                   Molecule &parent);
      /** Destructor
      *
      */
      virtual ~MolBondAngle();
      const Molecule& GetMolecule()const;
      Molecule& GetMolecule();
      virtual void XMLOutput(ostream &os,int indent=0)const;
      virtual void XMLInput(istream &is,const XMLCrystTag &tag);
      virtual REAL GetLogLikelihood()const;
      REAL GetAngle()const;
      REAL& Angle0();
      const MolAtom& GetAtom1()const;
      const MolAtom& GetAtom2()const;
      const MolAtom& GetAtom3()const;
      void SetAtom1(const MolAtom &at);
      void SetAtom2(const MolAtom &at);
      void SetAtom3(const MolAtom &at);
      //MolAtom& GetAtom1();
      //MolAtom& GetAtom2();
      //MolAtom& GetAtom3();
   private:
      /// The vector of the 3 atoms involved in the bond angle.
      vector<const MolAtom*> mvpAtom;
      REAL mAngle0,mDelta,mSigma;
      /// Parent Molecule
      Molecule *mpMol;
   #ifdef __WX__CRYST__
   public:
      WXCrystObjBasic *mpWXCrystObj;
      virtual WXCrystObjBasic* WXCreate(wxWindow*);
   #endif
};

/** Dihedral angle restraint between 4 atoms.
*
* The atoms involved need not be bonded.
*
*
*/
class MolDihedralAngle:public Restraint
{
   public:
      /** Constructor 
      *
      */
      MolDihedralAngle(const MolAtom &atom1, const MolAtom &atom2,
                       const MolAtom &atom3, const MolAtom &atom4,
                       const REAL angle, const REAL sigma, const REAL delta,
                       Molecule &parent);
      /** Destructor
      *
      */
      virtual ~MolDihedralAngle();
      const Molecule& GetMolecule()const;
      Molecule& GetMolecule();
      virtual void XMLOutput(ostream &os,int indent=0)const;
      virtual void XMLInput(istream &is,const XMLCrystTag &tag);
      virtual REAL GetLogLikelihood()const;
      REAL GetAngle()const;
      REAL& Angle0();
      const MolAtom& GetAtom1()const;
      const MolAtom& GetAtom2()const;
      const MolAtom& GetAtom3()const;
      const MolAtom& GetAtom4()const;
      void SetAtom1(const MolAtom& at);
      void SetAtom2(const MolAtom& at);
      void SetAtom3(const MolAtom& at);
      void SetAtom4(const MolAtom& at);
      //MolAtom& GetAtom1();
      //MolAtom& GetAtom2();
      //MolAtom& GetAtom3();
      //MolAtom& GetAtom4();
   private:
      /// The vector of the 4 atoms involved in the bond angle.
      vector<const MolAtom*> mvpAtom;
      REAL mAngle,mAngle0,mDelta,mSigma;
      /// Parent Molecule
      Molecule *mpMol;
   #ifdef __WX__CRYST__
   public:
      WXCrystObjBasic *mpWXCrystObj;
      virtual WXCrystObjBasic* WXCreate(wxWindow*);
   #endif
};

/** Ring class
*
* \note This class could be used to restrain atom positions for
* aromatic rings (planar restriction).
* \note We could also store the aromatic status of the ring.
*/
class MolRing
{
   public:
      MolRing(vector<MolBond*> &vRingBond);
      const vector<MolBond*>& GetBondList()const;
   private:
      vector<MolBond*> mvpBond;
};

/** A quaternion class, used to represent the orientation of the molecule.
* It may or may not be a unit quaternion.
* 
*/
class Quaternion
{
   public:
      /// Default constructor, yields q=(1,0,0,0)
      Quaternion();
      /// Creates a unit quaternion from its components (normalized automatically)
      Quaternion(const REAL q0,const REAL q1,const REAL q2,const REAL q3,bool unit=true);
      ~Quaternion();
      /// Create a rotation quaternion around a given vector for a given angle 
      static Quaternion  RotationQuaternion(const REAL ang,const REAL v1,const REAL v2,const REAL v3);
      /// Get the conjugate of this quaternion (== the inverse if unit quaternion)
      Quaternion  GetConjugate()const;
      /// Quaternion multiplication
      Quaternion operator*(const Quaternion &q)const;
      void operator*=(const Quaternion &q);
      void XMLOutput(ostream &os,int indent=0)const;
      void XMLInput(istream &is,const XMLCrystTag &tag);
      /// Rotate vector v=(v1,v2,v3). The rotated components are directly written
      void RotateVector(REAL &v1,REAL &v2, REAL &v3)const;
      /// Re-normalize the quaternion to unity. This should not be useful, except
      /// on individual component input, or after long calculations. And even
      /// if wrong, the rotation is independent of the norm of the quaternion.
      void Normalize();
      REAL GetNorm()const;
      const REAL& Q0()const;
      const REAL& Q1()const;
      const REAL& Q2()const;
      const REAL& Q3()const;
      REAL& Q0();
      REAL& Q1();
      REAL& Q2();
      REAL& Q3();
   private:
      /// The components of the quaternion z=(q0,v) with v=(q1,q2,q3)
      REAL mQ0,mQ1,mQ2,mQ3;
      bool mIsUniQuaternion;
};


/** Molecule : class for complex scatterer descriptions using
* cartesian coordinates with bond length/angle restraints, and
* moves either of individual atoms or using torsion bonds.
*
* This can also be used for non-organic compounds (polyhedras etc...)
* \note the parametrization is very different from ZScatterer: we keep 
* a list of x,y,z which do not use limits (they must not), but the 
* coordinates must be restrained or constrained from the expected
* bond lengths, angles and dihedral angles. The list of parameters is
* re-created in BeginOptimization() (except for the global x y z parameters
* for the global position of the Molecule, in fractionnal coordinates).
*
* \note : all atoms must be somehow connected
*/
class Molecule: public Scatterer
{
   public:
      /** Constructor
      *
      */
      Molecule(Crystal &cryst, const string &name="");
      /** Copy constructor
      *
      */
      Molecule(const Molecule &old);
      /** Destructor
      *
      */
      ~Molecule();
      virtual Molecule* CreateCopy() const;
      virtual const string& GetClassName() const;
      virtual void Print()const;
      virtual void XMLOutput(ostream &os,int indent=0)const;
      virtual void XMLInput(istream &is,const XMLCrystTag &tag);
      virtual void BeginOptimization(const bool allowApproximations=false,const bool enableRestraints=false);
      virtual void RandomizeConfiguration();
      virtual void GlobalOptRandomMove(const REAL mutationAmplitude,
                                       const RefParType *type);
      virtual REAL GetLogLikelihood()const;
      virtual void TagNewBestConfig()const;
      virtual int GetNbComponent() const;
      virtual const ScatteringComponentList& GetScatteringComponentList() const;
      virtual string GetComponentName(const int i) const;
      virtual ostream& POVRayDescription(ostream &os,
                                         const bool noSymmetrics=false)const;
      virtual void GLInitDisplayList(const bool onlyIndependentAtoms=false,
                                     const REAL xMin=-.1,const REAL xMax=1.1,
                                     const REAL yMin=-.1,const REAL yMax=1.1,
                                     const REAL zMin=-.1,const REAL zMax=1.1,
                                     const bool displayEnantiomer=false)const;
      /** Add an atom
      *
      *
      */
      void AddAtom(const REAL x, const REAL y, const REAL z,
                   const ScatteringPower *pPow,const string &name);
      /** Remove an atom
      *
      * This also removes all corresponding bonds, bond angles, etc...
      */
      void RemoveAtom(const MolAtom&);
      /** Add a bond
      *
      *
      */
      void AddBond(MolAtom &atom1, MolAtom &atom2,
                   const REAL length, const REAL sigma, const REAL delta,
                   const REAL bondOrder=1.);
      /** Remove a bond
      *
      */
      void RemoveBond(const MolBond&);
      /** Add a bond angle restraint
      *
      *
      */
      void AddBondAngle(MolAtom &atom1, MolAtom &atom2, MolAtom &atom3,
                        const REAL angle, const REAL sigma, const REAL delta);
      /** Remove a BondAngle
      *
      */
      void RemoveBondAngle(const MolBondAngle&);
      /** Add a dihedral angle restraint
      *
      *
      */
      void AddDihedralAngle(MolAtom &atom1, MolAtom &atom2, MolAtom &atom3, MolAtom &atom4,
                            const REAL angle, const REAL sigma, const REAL delta);
      /** Remove a dihedral angle
      *
      */
      void RemoveDihedralAngle(const MolDihedralAngle&);
      MolAtom &GetAtom(unsigned int i);
      const MolAtom &GetAtom(unsigned int i)const;
      MolAtom &GetAtom(const string &name);
      const MolAtom &GetAtom(const string &name)const;
      /** Minimize configuration from internal restraints (bond lengths, angles
      * and dihedral angles). Useful when adding manually atoms to get an initial
      * reasonable configuration. 
      */
      void OptimizeConformation(const long nbTrial=10000);
      const vector<MolAtom*>& GetAtomList()const;
      const vector<MolBond*>& GetBondList()const;
      const vector<MolBondAngle*>& GetBondAngleList()const;
      const vector<MolDihedralAngle*>& GetDihedralAngleList()const;
      vector<MolAtom*>& GetAtomList();
      vector<MolBond*>& GetBondList();
      vector<MolBondAngle*>& GetBondAngleList();
      vector<MolDihedralAngle*>& GetDihedralAngleList();
      /// Rotate a group of atoms around an axis defined by two atoms
      void RotateAtomGroup(const MolAtom &at1,const MolAtom &at2,
                           const set<unsigned long> &atoms, const REAL angle);
   private:
      virtual void InitRefParList();
      /** Build the list of rings in the molecule.
      *
      * The list is \e only rebuilt if the bond or atom list has changed,so
      * it should be safe to call again this function.
      *
      * \note So far this is a const method as the ring list just reflects the bond list 
      * and therefore is mutable (see Molecule::mvRing)... but maybe this could
      * change...
      */
      void BuildRingList()const;
      /** Build the Connectivity table
      *
      */
      void BuildConnectivityTable()const;
      /** Build the groups of atoms around each bond
      *
      */
      void BuildTorsionAtomGroupTable()const;
      /** Update the Molecule::mScattCompList from the cartesian coordinates
      * of all atoms, and the orientation parameters.
      */
      void UpdateScattCompList()const;
      /// Search a MolAtom from its name. Search begins at the end, and the
      /// first match is returned. returns mvAtom.rend() if no atom matches
      vector<MolAtom*>::reverse_iterator FindAtom(const string &name);
      /// Search a MolAtom from its name. Search begins at the end, and the
      /// first match is returned. returns mvAtom.rend() if no atom matches
      vector<MolAtom*>::const_reverse_iterator FindAtom(const string &name)const;
      /// Build options for this object
      void InitOptions();
      /** The list of scattering components
      *
      * this is mutable since it only reflects the list of atoms.
      */
      mutable ScatteringComponentList mScattCompList;
      /** The list of atoms
      *
      */
      vector<MolAtom*> mvpAtom;
      /** The list of bonds
      *
      */
      vector<MolBond*> mvpBond;
      /** The list of bond angles
      *
      */
      vector<MolBondAngle*> mvpBondAngle;
      /** The list of dihedral angles
      *
      */
      vector<MolDihedralAngle*> mvpDihedralAngle;
      /** List of Bonds for each atom.
      *
      * This duplicates the information in Molecule::mvBond
      */
      map<MolAtom* , std::vector<MolBond*> > mvAtomBond;

      /** The list of rings
      *
      * \note this only reflects the bond list, so it is mutable.
      */
      mutable vector<MolRing> mvRing;
      /** The unit quaternion defining the orientation
      *
      */
      Quaternion mQuat;
      // Clocks
         RefinableObjClock mClockAtomList;
         RefinableObjClock mClockBondList;
         RefinableObjClock mClockBondAngleList;
         RefinableObjClock mClockDihedralAngleList;
         RefinableObjClock mClockRingList;
         RefinableObjClock mClockAtomPosition;
         RefinableObjClock mClockAtomScattPow;
         RefinableObjClock mClockOrientation;
         mutable RefinableObjClock mClockLogLikelihood;
         mutable RefinableObjClock mClockConnectivityTable;
         mutable RefinableObjClock mClockmTorsionAtomGroupTable;
         
      // For local minimization (EXPERIMENTAL)
         unsigned long mLocalParamSet;
         mutable unsigned long mRandomConformChangeNbTest;
         mutable unsigned long mRandomConformChangeNbAccept;
         mutable REAL mRandomConformChangeTemp;
         REAL mLastLogLike;
         bool mIsSelfOptimizing;
      /// OPtion for the different types of flexibility possible for this
      /// molecule: rigid body, free atoms + restraints, torsion angles...
      /// \warning still EXPERIMENTAL !
      RefObjOpt mFlexModel;
      
      /// Connectivity table: for each atom, keep the list of atoms
      /// bonded to it. All atoms are referenced from their index.
      mutable map<unsigned long,set<unsigned long> > mConnectivityTable;
      /** Groups of atoms around each bond: for each bond, keep a list of group of atoms
      * on each side of the bond.
      *
      * The map runs on all bonds, referenced from their index in the bond list. For each bond,
      * there are two lists of group of atoms, corresponding to the two ends of the bond.
      * Each list references all the group of atoms starting from one end of the bond.
      * e.g. if A-B make one bond, and A is directly connected to 3 atoms A1 A2 A3, then
      * there are 3 atom groups listed for A, (A,A1,...), (A,A2,...) and (A,A3,...).
      *
      * This is useful to rotate atoms around a torsion bond, with the ability to rotate
      * all atoms on one side of the bond, or only one group of atoms starting from one
      * end of the bond.
      */
      mutable map<unsigned long,pair<list<set<unsigned long> >,list<set<unsigned long> > > > mTorsionAtomGroupTable;
      /** The index of bonds with free rotation
      *
      */
      mutable vector<unsigned long> mTorsionBondIndex;
   /// The current log(likelihood)
   mutable REAL mLogLikelihood;
   #ifdef __WX__CRYST__
   public:
      virtual WXCrystObjBasic* WXCreate(wxWindow*);
   #endif
};

}//namespace
#endif
