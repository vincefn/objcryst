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

/// Structure holding 3 coordinates, or deriviatives with respect to each of these coordinates 
struct XYZ
{
   XYZ(REAL x=0,REAL y=0,REAL z=0);
   REAL x,y,z;
};

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
      /// Flag this atom as being in a ring (or not). This is a const method
      /// because the existence of a ring is only a consequence of the connectivity
      /// of the Molecule.
      void SetIsInRing(const bool r)const;
      bool IsInRing()const;
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
      /// Is the atom in a ring ?
      mutable bool mIsInRing;
   #ifdef __WX__CRYST__
   public:
      WXCrystObjBasic *mpWXCrystObj;
      virtual WXCrystObjBasic* WXCreate(wxWindow*);
      WXCrystObjBasic* WXGet();
      void WXDelete();
      void WXNotifyDelete();
   #endif
};

/// Get The Bond Length between two atoms
REAL GetBondLength(const MolAtom&,const MolAtom&);
/// Get The Bond Angle of 3 atoms
REAL GetBondAngle(const MolAtom&,const MolAtom&,const MolAtom&);
/// Get The dihedral angle defined by 4 atoms
REAL GetDihedralAngle(const MolAtom&,const MolAtom&,const MolAtom&,const MolAtom&);

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
      MolBond(MolAtom &atom1, MolAtom &atom2,
              const REAL length, const REAL sigma, const REAL delta,
              Molecule &parent,const REAL bondOrder=1.);
      /** Destructor
      *
      * Notifies the atoms that the bond has disappeared.
      */
      virtual ~MolBond();
      const Molecule& GetMolecule()const;
      Molecule& GetMolecule();
      /// Name of the bond, e.g. "C3-O4"
      string GetName()const;
      virtual void XMLOutput(ostream &os,int indent=0)const;
      virtual void XMLInput(istream &is,const XMLCrystTag &tag);
      virtual REAL GetLogLikelihood(const bool calcDeriv=false, const bool recalc=true)const;
      /// Get the derivative of the bond length, given the derivatives of the atom positions
      /// This requires that GetLogLikelihood(calcDeriv=true) be called first.
      /// If llk=true, this will return the derivative of the llk rather than the derivative of the length or angle
      REAL GetDeriv(const std::map<const MolAtom*,XYZ> &m, const bool llk=false)const;
      const MolAtom& GetAtom1()const;
      const MolAtom& GetAtom2()const;
      MolAtom& GetAtom1();
      MolAtom& GetAtom2();
      void SetAtom1(MolAtom &at1);
      void SetAtom2(MolAtom &at2);
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
      pair<MolAtom*,MolAtom*> mAtomPair;
      REAL mLength0,mDelta,mSigma;
      REAL mBondOrder;
      bool mIsFreeTorsion;
      /// Parent Molecule
      Molecule *mpMol;
      /// Stored log(likelihood)
      mutable REAL mLLK;
      /** Derivatives of the bond length with respect to the coordinates of the atoms
      *
      * The derivatives are calculated in MolBond::GetLogLikelihood(true)
      */
      mutable XYZ mDerivAtom1,mDerivAtom2;
      /** The factor used to change the derivative of the length/angle, to the derivative
      * of the log(likelihood). e.g. (for mDelta=0) \f$ mDerivLLKCoeff = \frac{L-L_0}{\sigma^2} \f$
      */
      mutable REAL mDerivLLKCoeff;
   #ifdef __WX__CRYST__
   public:
      WXCrystObjBasic *mpWXCrystObj;
      virtual WXCrystObjBasic* WXCreate(wxWindow*);
      WXCrystObjBasic* WXGet();
      void WXDelete();
      void WXNotifyDelete();
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
      MolBondAngle(MolAtom &atom1,MolAtom &atom2,MolAtom &atom3,
                   const REAL angle, const REAL sigma, const REAL delta,
                   Molecule &parent);
      /** Destructor
      *
      */
      virtual ~MolBondAngle();
      const Molecule& GetMolecule()const;
      Molecule& GetMolecule();
      string GetName()const;
      virtual void XMLOutput(ostream &os,int indent=0)const;
      virtual void XMLInput(istream &is,const XMLCrystTag &tag);
      virtual REAL GetLogLikelihood(const bool calcDeriv=false, const bool recalc=true)const;
      /// Get the derivative of the angle, given the derivatives of the atom positions
      /// This requires that GetLogLikelihood(calcDeriv=true) be called first
      /// If llk=true, this will return the derivative of the llk rather than the derivative of the length or angle
      REAL GetDeriv(const std::map<const MolAtom*,XYZ> &m, const bool llk=false)const;
      REAL GetAngle()const;
      REAL& Angle0();
      REAL& AngleDelta();
      REAL& AngleSigma();
      REAL GetAngle0()const;
      REAL GetAngleDelta()const;
      REAL GetAngleSigma()const;
      void SetAngle0(const REAL angle);
      void SetAngleDelta(const REAL delta);
      void SetAngleSigma(const REAL sigma);
      const MolAtom& GetAtom1()const;
      const MolAtom& GetAtom2()const;
      const MolAtom& GetAtom3()const;
      void SetAtom1(MolAtom &at);
      void SetAtom2(MolAtom &at);
      void SetAtom3(MolAtom &at);
      MolAtom& GetAtom1();
      MolAtom& GetAtom2();
      MolAtom& GetAtom3();
      bool IsFlexible()const;
      void SetFlexible(const bool isInRing);
   private:
      /// The vector of the 3 atoms involved in the bond angle.
      vector<MolAtom*> mvpAtom;
      REAL mAngle0,mDelta,mSigma;
      /// Parent Molecule
      Molecule *mpMol;
      /// When using the user-chosen flexibility model, this allows
      /// some flexibility for this bond angle, i.e. the bond angle
      /// does not remain strictly rigid, and is still restrained.
      bool mIsFlexible;
      /// Stored log(likelihood)
      mutable REAL mLLK;
      /** Partial derivatives of the angle with respect to the coordinates of the atoms
      *
      * The derivatives are calculated in MolBondAngle::GetLogLikelihood(true)
      */
      mutable XYZ mDerivAtom1,mDerivAtom2,mDerivAtom3;
      /** The factor used to change the derivative of the length/angle, to the derivative
      * of the log(likelihood). e.g. (for mDelta=0) \f$ mDerivLLKCoeff = \frac{L-L_0}{\sigma^2} \f$
      */
      mutable REAL mDerivLLKCoeff;
   #ifdef __WX__CRYST__
   public:
      WXCrystObjBasic *mpWXCrystObj;
      virtual WXCrystObjBasic* WXCreate(wxWindow*);
      WXCrystObjBasic* WXGet();
      void WXDelete();
      void WXNotifyDelete();
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
      MolDihedralAngle(MolAtom &atom1, MolAtom &atom2,
                       MolAtom &atom3, MolAtom &atom4,
                       const REAL angle, const REAL sigma, const REAL delta,
                       Molecule &parent);
      /** Destructor
      *
      */
      virtual ~MolDihedralAngle();
      const Molecule& GetMolecule()const;
      Molecule& GetMolecule();
      string GetName()const;
      virtual void XMLOutput(ostream &os,int indent=0)const;
      virtual void XMLInput(istream &is,const XMLCrystTag &tag);
      virtual REAL GetLogLikelihood(const bool calcDeriv=false, const bool recalc=true)const;
      /// Get the derivative of the Angle, given the derivatives of the atom positions
      /// This requires that GetLogLikelihood(calcDeriv=true) be called first
      /// If llk=true, this will return the derivative of the llk rather than the derivative of the length or angle
      REAL GetDeriv(const std::map<const MolAtom*,XYZ> &m, const bool llk=false)const;
      REAL GetAngle()const;
      REAL& Angle0();
      REAL& AngleDelta();
      REAL& AngleSigma();
      REAL GetAngle0()const;
      REAL GetAngleDelta()const;
      REAL GetAngleSigma()const;
      void SetAngle0(const REAL angle);
      void SetAngleDelta(const REAL delta);
      void SetAngleSigma(const REAL sigma);
      const MolAtom& GetAtom1()const;
      const MolAtom& GetAtom2()const;
      const MolAtom& GetAtom3()const;
      const MolAtom& GetAtom4()const;
      void SetAtom1(MolAtom& at);
      void SetAtom2(MolAtom& at);
      void SetAtom3(MolAtom& at);
      void SetAtom4(MolAtom& at);
      MolAtom& GetAtom1();
      MolAtom& GetAtom2();
      MolAtom& GetAtom3();
      MolAtom& GetAtom4();
   private:
      /// The vector of the 4 atoms involved in the bond angle.
      vector<MolAtom*> mvpAtom;
      REAL mAngle,mAngle0,mDelta,mSigma;
      /// Parent Molecule
      Molecule *mpMol;
      /// Stored log(likelihood)
      mutable REAL mLLK;
      /** Partial derivatives of the angle with respect to the coordinates of the atoms
      *
      * The derivatives are calculated in MolBondAngle::GetLogLikelihood(true)
      */
      mutable XYZ mDerivAtom1,mDerivAtom2,mDerivAtom3,mDerivAtom4;
      /** The factor used to change the derivative of the length/angle, to the derivative
      * of the log(likelihood). e.g. (for mDelta=0) \f$ mDerivLLKCoeff = \frac{L-L_0}{\sigma^2} \f$
      */
      mutable REAL mDerivLLKCoeff;
   #ifdef __WX__CRYST__
   public:
      WXCrystObjBasic *mpWXCrystObj;
      virtual WXCrystObjBasic* WXCreate(wxWindow*);
      WXCrystObjBasic* WXGet();
      void WXDelete();
      void WXNotifyDelete();
   #endif
};

//typedef std::set<MolAtom *> RigidGroup;
class RigidGroup:public std::set<MolAtom *>
{
   public:
      std::string GetName()const;
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
      MolRing();
      const std::list<MolAtom*>& GetAtomList()const;
      std::list<MolAtom*>& GetAtomList();
   private:
      std::list<MolAtom*> mvpAtom;
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
/** Abstract base Stretch Mode for Molecule objects
*
*/
struct StretchMode
{
   virtual ~StretchMode();
   /// Calculate the derivative of the Molecule's Log(likelihood) and atomic psoitions
   /// versus a change of the bond length. The result is stored in mLLKDeriv and mLLKDerivXYZ,
   /// as well as in the various lists of restraints broken by this mode.
   virtual void CalcDeriv()const=0;
   /// Print one-line list of atoms moved
   virtual void Print(ostream &os,bool full=true)const=0;
   /// Move the atoms according to this mode
   virtual void Stretch(const REAL change)=0;
   /// Move the atoms according to this mode, randomly
   virtual void RandomStretch(const REAL amplitude)=0;
   /// List of bond restraints affected by this mode
   /// The key is the restraint, the value is the derivative of the LLK associated
   std::map<const MolBond*,REAL> mvpBrokenBond;
   /// List of bond angle restraints modified by this mode
   /// The key is the restraint, the value is the derivative of the LLK associated
   std::map<const MolBondAngle*,REAL> mvpBrokenBondAngle;
   /// List of dihedral angle restraints modified by this mode
   /// The key is the restraint, the value is the derivative of the LLK associated
   std::map<const MolDihedralAngle*,REAL> mvpBrokenDihedralAngle;
   /// Derivative of the Molecule's Log(likelihood) versus a change of the bond length.
   mutable REAL mLLKDeriv;
   /// Derivative of the atomic positions versus a change of the bond length.
   mutable std::map<const MolAtom*,XYZ> mDerivXYZ;
   /// The Molecule corresponding to this stretch mode
   Molecule *mpMol;
   /** The recommended change amplitude, for a base global optimization
   * displacement, to obtain an average 0.1 Angstroem displacement.
   *
   * This is learnt at the beginning of an optimization.
   *
   * This can be superseeded to respect any restraint.
   */
   REAL mBaseAmplitude;
};

/** Group of atoms for random moves changing a bond length. 
*
* This should be merged (or have an inheritance relation) with MolBond.
*/
struct StretchModeBondLength:public StretchMode
{
   /** Constructor
   * If pBond!=0, the bond length restraint is respected.
   */
   StretchModeBondLength(MolAtom &at0,MolAtom &at1,const MolBond *pBond);
   virtual ~StretchModeBondLength();
   /// Calculate the derivative of the Molecule's Log(likelihood) and atomic psoitions
   /// versus a change of the bond length. The result is stored in mLLKDeriv and mLLKDerivXYZ,
   /// as well as in the various lists of restraints broken by this mode.
   virtual void CalcDeriv()const;
   /// Print one-line list of atoms moved
   virtual void Print(ostream &os,bool full=true)const;
   /// Move the atoms according to this mode
   virtual void Stretch(const REAL change);
   /// Move the atoms according to this mode, randomly
   virtual void RandomStretch(const REAL amplitude);
   /// The first atom (fixed).
   MolAtom * mpAtom0;
   /// The second atom  (first atom moved)
   MolAtom * mpAtom1;
   /// The (optional) bond length which this stretch mode should respect.
   const MolBond *mpBond;
   /// The set of atoms that are to be translated, including at1.
   set<MolAtom *> mvTranslatedAtomList;
};

/** Atoms moved when changing a bond angle.
*
* This should be merged (or have an inheritance relation) with MolBondAngle.
*/
struct StretchModeBondAngle:public StretchMode
{
   /** Constructor
   * If pBondAngle!=0, the bond angle length restraint is respected.
   */
   StretchModeBondAngle(MolAtom &at0,MolAtom &at1,MolAtom &at2,
                        const MolBondAngle *pBondAngle);
   virtual ~StretchModeBondAngle();
   /// Calculate the derivative of the Molecule's Log(likelihood) and atomic psoitions
   /// versus a change of the angle. The result is stored in mLLKDeriv and mLLKDerivXYZ,
   /// as well as in the various lists of restraints broken by this mode.
   virtual void CalcDeriv()const;
   /// Print one-line list of atoms moved
   virtual void Print(ostream &os,bool full=true)const;
   /// Move the atoms according to this mode
   virtual void Stretch(const REAL change);
   /// Move the atoms according to this mode, randomly
   virtual void RandomStretch(const REAL amplitude);
   /// The first atom
   MolAtom * mpAtom0;
   /// The second atom 
   MolAtom * mpAtom1;
   /// The third atom
   MolAtom * mpAtom2;
   /// The (optional) bond angle restraint which this stretch mode should respect.
   const MolBondAngle *mpBondAngle;
   /// The set of atoms that are to be rotated around the direction going
   /// through at1 and perpendicular to the at0-at1-at2 plane.
   set<MolAtom *> mvRotatedAtomList;
};
/** Atoms moved when rotated around a bond at0-at1-at2-at3
*
* This should be merged (or have an inheritance relation) with MolDihedralAngle
*/
struct StretchModeTorsion:public StretchMode
{
   /** Constructor
   * If pDihedralAngle!=0, the dihedral angle length restraint is respected.
   */
   StretchModeTorsion(MolAtom &at1,MolAtom &at2,
                      const MolDihedralAngle *pDihedralAngle);
   virtual ~StretchModeTorsion();
   /// Calculate the derivative of the Molecule's Log(likelihood) and atomic psoitions
   /// versus a change of the angle. The result is stored in mLLKDeriv and mLLKDerivXYZ,
   /// as well as in the various lists of restraints broken by this mode.
   virtual void CalcDeriv()const;
   /// Print one-line list of atoms moved
   virtual void Print(ostream &os,bool full=true)const;
   /// Move the atoms according to this mode
   virtual void Stretch(const REAL change);
   /// Move the atoms according to this mode, randomly
   virtual void RandomStretch(const REAL amplitude);
   /// The first atom
   MolAtom * mpAtom1;
   /// The second atom 
   MolAtom * mpAtom2;
   /// The (optional) bond angle restraint which this stretch mode should respect.
   /// The mpAtom1 and mpAtom2 must be the central atoms of this restraint.
   const MolDihedralAngle *mpDihedralAngle;
   /// The set of atoms that are to be rotated around at1-at2
   set<MolAtom *> mvRotatedAtomList;
};

/** Atoms moved *between* two other atoms, using a "twist"
*of their positions - only small twists of their positions
* are allowed to avoid breaking restraints too much. Also,
* the atoms in the middle between the two reference atoms are more
* displaced than the one closer, to help avoid breaking restraints.
*/
struct StretchModeTwist:public StretchMode
{
   /** Constructor
   * If pDihedralAngle!=0, the dihedral angle length restraint is respected.
   */
   StretchModeTwist(MolAtom &at1,MolAtom &at2);
   virtual ~StretchModeTwist();
   /// Calculate the derivative of the Molecule's Log(likelihood) and atomic psoitions
   /// versus a change of the angle. The result is stored in mLLKDeriv and mLLKDerivXYZ,
   /// as well as in the various lists of restraints broken by this mode.
   virtual void CalcDeriv()const;
   /// Print one-line list of atoms moved
   virtual void Print(ostream &os,bool full=true)const;
   /// Move the atoms according to this mode
   virtual void Stretch(const REAL change);
   /// Move the atoms according to this mode, randomly
   virtual void RandomStretch(const REAL amplitude);
   /// The first atom
   MolAtom * mpAtom1;
   /// The second atom 
   MolAtom * mpAtom2;
   /// The set of atoms that are to be rotated around at1-at2
   set<MolAtom *> mvRotatedAtomList;
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
                                         const CrystalPOVRayOptions &options)const;
      virtual void GLInitDisplayList(const bool onlyIndependentAtoms=false,
                                     const REAL xMin=-.1,const REAL xMax=1.1,
                                     const REAL yMin=-.1,const REAL yMax=1.1,
                                     const REAL zMin=-.1,const REAL zMax=1.1,
                                     const bool displayEnantiomer=false,
                                     const bool displayNames=false)const;
      /** Add an atom
      *
      *
      */
      void AddAtom(const REAL x, const REAL y, const REAL z,
                   const ScatteringPower *pPow,const string &name,
                   const bool updateDisplay=true);
      /** Remove an atom. Returns the iterator to the next atom in the list.
      *
      * This also removes all corresponding bonds, bond angles, etc...
      */
      vector<MolAtom*>::iterator RemoveAtom(MolAtom&);
      /** Add a bond
      *
      *
      */
      void AddBond(MolAtom &atom1, MolAtom &atom2,
                   const REAL length, const REAL sigma, const REAL delta,
                   const REAL bondOrder=1.,
                   const bool updateDisplay=true);
      /** Remove a bond. Returns the iterator to the next bond in the list.
      *
      */
      vector<MolBond*>::iterator RemoveBond(const MolBond&);
      /** Searches whether a bond between two atoms already exists.
      *
      * If no bond is found, returns Molecule::mvpAtom.end().
      */
      vector<MolBond*>::const_iterator FindBond(const MolAtom&,const MolAtom&)const;
      /** Searches whether a bond between two atoms already exists.
      *
      * If no bond is found, returns Molecule::mvpAtom.end().
      */
      vector<MolBond*>::iterator FindBond(const MolAtom&,const MolAtom&);
      /** Add a bond angle restraint
      *
      *
      */
      void AddBondAngle(MolAtom &atom1, MolAtom &atom2, MolAtom &atom3,
                        const REAL angle, const REAL sigma, const REAL delta,
                        const bool updateDisplay=true);
      /** Remove a BondAngle
      *
      */
      vector<MolBondAngle*>::iterator RemoveBondAngle(const MolBondAngle&);
      /** Searches whether a bond between three atoms already exists,
      * searching for either (at1,at2,at3) and (at3,at2,at1), as these are equivalent.
      *
      * If no bond angle is found, returns Molecule::mvpBondAngle.end().
      */
      vector<MolBondAngle*>::const_iterator FindBondAngle(const MolAtom& at1,const MolAtom&at0,const MolAtom&at2)const;
      /** Add a dihedral angle restraint
      *
      *
      */
      void AddDihedralAngle(MolAtom &atom1, MolAtom &atom2,
                            MolAtom &atom3, MolAtom &atom4,
                            const REAL angle, const REAL sigma, const REAL delta,
                            const bool updateDisplay=true);
      /** Remove a dihedral angle
      *
      */
      vector<MolDihedralAngle*>::iterator RemoveDihedralAngle(const MolDihedralAngle&);
      /** Searches whether a dihedral between four atoms already exists,
      * searching for either (at1,at2,at3,at4) and (at4,at3,at2,at1), as these are equivalent.
      *
      * If no dihedral angle is found, returns Molecule::mvpDihedralAngle.end().
      */
      vector<MolDihedralAngle*>::const_iterator FindDihedralAngle(const MolAtom &at1,
                                                                  const MolAtom &at2,
                                                                  const MolAtom &at3,
                                                                  const MolAtom &at4)const;
      /** Add a rigid group of atoms. See Molecule::mvRigidGroup
      */
      void AddRigidGroup(const RigidGroup&,const bool updateDisplay=true);
      /** Remove a rigid group of atoms. See Molecule::mvRigidGroup
      */
      std::vector<RigidGroup*>::iterator RemoveRigidGroup(const RigidGroup &group);

      MolAtom &GetAtom(unsigned int i);
      const MolAtom &GetAtom(unsigned int i)const;
      MolAtom &GetAtom(const string &name);
      const MolAtom &GetAtom(const string &name)const;

      /// Search a MolAtom from its name. Search begins at the end, and the
      /// first match is returned. returns mvAtom.rend() if no atom matches
      vector<MolAtom*>::reverse_iterator FindAtom(const string &name);
      /// Search a MolAtom from its name. Search begins at the end, and the
      /// first match is returned. returns mvAtom.rend() if no atom matches
      vector<MolAtom*>::const_reverse_iterator FindAtom(const string &name)const;

      /** Minimize configuration from internal restraints (bond lengths, angles
      * and dihedral angles). Useful when adding manually atoms to get an initial
      * reasonable configuration. 
      */
      void OptimizeConformation(const long nbTrial=10000,const REAL stopCost=0.);
      const vector<MolAtom*>& GetAtomList()const;
      const vector<MolBond*>& GetBondList()const;
      const vector<MolBondAngle*>& GetBondAngleList()const;
      const vector<MolDihedralAngle*>& GetDihedralAngleList()const;
      vector<MolAtom*>& GetAtomList();
      vector<MolBond*>& GetBondList();
      vector<MolBondAngle*>& GetBondAngleList();
      vector<MolDihedralAngle*>& GetDihedralAngleList();
      
      list<StretchModeBondLength>& GetStretchModeBondLengthList();
      list<StretchModeBondAngle>& GetStretchModeBondAngleList();
      list<StretchModeTorsion>& GetStretchModeTorsionList();

      const list<StretchModeBondLength>& GetStretchModeBondLengthList()const;
      const list<StretchModeBondAngle>& GetStretchModeBondAngleList()const;
      const list<StretchModeTorsion>& GetStretchModeTorsionList()const;

      /** List of rigid group of atoms. See Molecule::mvRigidGroup
      */
      const std::vector<RigidGroup *>& GetRigidGroupList()const;
      
      /** Rotate a group of atoms around an axis defined by two atoms
      *
      * \param keepCenter: if true, the coordinates of the molecule are modified
      * so that only the rotated atoms are moved.
      */
      void RotateAtomGroup(const MolAtom &at1,const MolAtom &at2,
                           const set<MolAtom *> &atoms, const REAL angle,
                           const bool keepCenter=true);
      /** Rotate a group of atoms around an axis defined by one atom and a vector
      *
      * \param keepCenter: if true, the coordinates of the molecule are modified
      * so that only the rotated atoms are moved.
      */
      void RotateAtomGroup(const MolAtom &at,const REAL vx,const REAL vy,const REAL vz,
                           const set<MolAtom *> &atoms, const REAL angle,
                           const bool keepCenter=true);
      /** Translate a group of atoms in a given direction
      *
      * \param keepCenter: if true, the coordinates of the molecule are modified
      * so that only the translated atoms are moved.
      */
      void TranslateAtomGroup(const set<MolAtom *> &atoms, 
                              const REAL dx,const REAL dy,const REAL dz,
                              const bool keepCenter=true);
      /// Print the status of all restraints (bond length, angles...)
      void RestraintStatus(ostream &os)const;
      /// Get the connectivity table
      const map<MolAtom *,set<MolAtom *> > & GetConnectivityTable();
      /// get the clock associated to the list of bonds
      RefinableObjClock& GetBondListClock();
      /// get the clock associated to the list of bonds
      const RefinableObjClock& GetBondListClock()const;
      /** Add dihedral angles so as to rigidify the Molecule.
      *
      * In practice, for every sequence of atoms A-B-C-D, add the dihedral angle
      * defined by these 4 atoms, unless either ABC or BCD are aligned (angle below 10°).
      *
      * No duplicate dihedral angle is generated.
      */
      void RigidifyWithDihedralAngles();
      /** Stretch a bond, while respecting the Restraint (if any).
      *
      *\return the \e actual change in bond length.
      *\param: the desired change in bond length. This will be the actual change \e if
      * there is no restraint \e or if the restraint is constant in this range.
      */
      REAL BondLengthRandomChange(const StretchModeBondLength& mode, const REAL amplitude,
                                  const bool respectRestraint=true);
      /** change a bond angle, while respecting the Restraint (if any).
      *
      *\return the \e actual change in bond angle.
      *\param: the desired angular change. This will be the actual change \e if
      * there is no restraint \e or if the restraint is constant in this range.
      */
      REAL BondAngleRandomChange(const StretchModeBondAngle& mode, const REAL amplitude,
                                 const bool respectRestraint=true);
      /** Change a dihedral angle, while respecting the Restraint (if any).
      *
      *\return the \e actual change in bond angle.
      *\param: the desired angular change. This will be the actual change \e if
      * there is no restraint \e or if the restraint is constant in this range.
      */
      REAL DihedralAngleRandomChange(const StretchModeTorsion& mode, const REAL amplitude,
                                     const bool respectRestraint=true);
      /// Get the atom defining the origin of the Molecule
      /// Equal to 0 if no atom as been set.
      const MolAtom* GetCenterAtom()const;
      /// Get the atom defining the origin of the Molecule
      /// Equal to 0 if no atom as been set.
      void SetCenterAtom(const MolAtom &at);
   public:
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
      void BuildRingList();
      /** Build the Connectivity table
      *
      */
      void BuildConnectivityTable();
      /** Build the groups of atoms that will be rotated during global optimization.
      *
      * This is not const because we temporarily modify the molecule conformation
      * to test which RotorGroups are forbidden by restraints (but it should be const).
      */
      void BuildRotorGroup();
      /** Tune the rotation amplitude for free torsions and for the overall Molecule
      * Rotation.
      *
      * This should be done after Molecule::BuildRotorGroup();
      */
      void TuneGlobalOptimRotationAmplitude();
      /** Build the groups of atoms that can be flipped.
      *
      * This is not const because we temporarily modify the molecule conformation
      * to test which FlipGroups are forbidden by restraints (but it should be const).
      */
      void BuildFlipGroup();
      /** Build the groups of atoms moved when stretching a bond length, while
      * respecting the Molecule restraints.
      */
      void BuildStretchModeBondLength();
      /** Build the groups of atoms moved when changing a bond angle, while
      * respecting the Molecule restraints.
      */
      void BuildStretchModeBondAngle();
      /** Build the groups of atoms moved when changing a dihedral angle, while
      * respecting the Molecule restraints.
      */
      void BuildStretchModeTorsion();
      /** Build the groups of atoms used to twist internally the Molecule, e.g. by
      * rotating one chain of atoms between 2 given atoms.
      */
      void BuildStretchModeTwist();
      /** Separate StretchMode that break more than their assigned restraint from others.
      * See Molecule::mvpStretchModeFree and Molecule::mvpStretchModeNotFree
      */
      void BuildStretchModeGroups();
      /** Update the Molecule::mScattCompList from the cartesian coordinates
      * of all atoms, and the orientation parameters.
      */
      void UpdateScattCompList()const;
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
      /** Rigid groups of atoms. This group will be kept \e strictly rigid,
      * preventing the use of any stretch mode altering their relative position.
      * The entire group of atoms can however be rotated or translated.
      */
      std::vector<RigidGroup*> mvRigidGroup;
      /** The list of rings
      *
      * \note this only reflects the bond list, so it is mutable.
      */
      mutable list<MolRing> mvRing;
      /** The unit quaternion defining the orientation
      *
      */
      Quaternion mQuat;
      /** Base Rotation amplitude (in radians) for the Molecule, so that
      * the average atomic displacement is equal to 0.1 A
      *
      * Default=0.02*pi
      */
      REAL mBaseRotationAmplitude;
      // Clocks
         RefinableObjClock mClockAtomList;
         RefinableObjClock mClockBondList;
         RefinableObjClock mClockBondAngleList;
         RefinableObjClock mClockDihedralAngleList;
         RefinableObjClock mClockAtomPosition;
         RefinableObjClock mClockAtomScattPow;
         RefinableObjClock mClockOrientation;
         mutable RefinableObjClock mClockLogLikelihood;
         mutable RefinableObjClock mClockConnectivityTable;
         mutable RefinableObjClock mClockRingList;
         mutable RefinableObjClock mClockRotorGroup;
         mutable RefinableObjClock mClockFlipGroup;
         mutable RefinableObjClock mClockStretchModeBondLength;
         mutable RefinableObjClock mClockStretchModeBondAngle;
         mutable RefinableObjClock mClockStretchModeTorsion;
         mutable RefinableObjClock mClockStretchModeTwist;
         mutable RefinableObjClock mClockRigidGroup;
         
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

      /** Option to automatically optimize the starting conformation, if the
      * total restraint cost is too high. This is done in BeginOptimization().
      *
      * This is enabled by default, and should be disabled by people who already
      * supply a good starting conformation for their molecule.
      */
      RefObjOpt mAutoOptimizeConformation;

      /** Option to optimize the Molecule's orientation. Useful to completely
      * fix the Molecule.
      */
      RefObjOpt mOptimizeOrientation;

      /** Option to choose the center of rotation of the Molecule for the global orientation
      * either as the geometrical center, or as a given atom.
      */
      RefObjOpt mMoleculeCenter;
      
      /** Atom chosen as center of rotation, if mRotationCenter is set to use
      * an atom rather than the geometrical center.
      */
      const MolAtom* mpCenterAtom;
      
      /// Connectivity table: for each atom, keep the list of atoms
      /// bonded to it. All atoms are referenced from their index.
      mutable map<MolAtom *,set<MolAtom *> > mConnectivityTable;
      /** Defines a group of atoms which can be rotated around an axis defined
      * by two other atoms.
      */
      struct RotorGroup
      {
         /** Constructor, with the two atoms around which the rotation
         * shall be made. The list of atoms to be rotated is initially empty.
         */
         RotorGroup(const MolAtom &at1,const MolAtom &at2);
         /// The first atom defining the rotation axis
         const MolAtom * mpAtom1;
         /// The second atom defining the rotation axis
         const MolAtom * mpAtom2;
         /// The set of atoms that are to be rotated
         set<MolAtom *> mvRotatedAtomList;
         /** The recommended rotation amplitude, for a base global optimization
         * displacement, to obtain an average 0.1 Angstroem displacement 
         * per atom (pi*0.04 by default)
         *
         * This is learnt at the beginning of an optimization, i.e. in
         * Molecule::BuildRotorGroup()
         */
         REAL mBaseRotationAmplitude;
      };
      /** List of RotorGroups corresponding to free torsion bonds.
      *
      * In this list are list of atoms on one side of a bond, that can be rotated
      * freely around this bond. Each bond is listed only once, with
      * the side which has the smallest number of atoms.
      */
      mutable list<RotorGroup> mvRotorGroupTorsion;
      /** List of RotorGroups corresponding to free torsion bonds, but with only
      * one chain of atoms listed.
      *
      * The difference with Molecule::mRotorGroupTorsion is that if the bond is
      * A-B, with atom A linked with atoms A1,A2,A3, in this list only one chain
      * (starting either from A1, A2 or A3) will be rotated, instead of the 3 chains.
      * This is useful when searching for the absolute configuration of atoms.
      */
      mutable list<RotorGroup> mvRotorGroupTorsionSingleChain;
      /** List of RotorGroups for internal rotations. This lists groups of atoms
      * that can be rotated \e between two given atoms. This is useful to alter
      * the conformation of large rings, where no free torsion bonds exists, and
      * also for long flexible chains.
      */
      mutable list<RotorGroup> mvRotorGroupInternal;
      /** When 3(A1..1n) or more atoms are connected to a same atom A, it defines
      * a 'flip' group, where it is possible to rotate bonds to their symmetric
      * with respect to one plane defined by atoms Ai-A-Aj. This is useful to
      * flip the absolute configuration for asymmetric centers. Note that the bond
      * is only rotated, so that the entire group is not mirrored (no absolute configuration 
      * is broken in the group).
      *
      * Also, a FlipGroup can correspond to a 180° rotation exchanging Ai and Aj
      * (rotating the two chains around the bissecting angle of bonds A-Ai and A-Aj)
      */
      struct FlipGroup
      {
         /** Constructor, with the central atom.
         */
         FlipGroup(const MolAtom &at0,const MolAtom &at1,const MolAtom &at2);
         /// The atom which is an asymmetric center
         const MolAtom * mpAtom0;
         /// The first atom defining the rotation axis
         const MolAtom * mpAtom1;
         /// The second atom defining the rotation axis
         const MolAtom * mpAtom2;
         /// The set of atoms that are to be rotated during the flip. The first
         /// atom is the one bonded to the central atom, whose bond will be flipped
         /// with respect to the plane defined by (at1,at0,at2).
         ///
         /// However, if this atom is identical to mpAtom0, then this indicates that
         /// a 180° rotation exchanging atom1 and atom2 is to be performed.
         list<pair<const MolAtom *,set<MolAtom *> > > mvRotatedChainList;
         /// Number of times this flip has been tried, and the number of times
         /// it has been accepted. Used in Molecule::GlobalOptRandomMove,
         /// to avoid flips that break some restraint (and deciding which flips
         /// break some restraint is difficult before having a real conformation).
         mutable unsigned long mNbTest,mNbAccept;
      };
      /// Flip a group of atom. See Molecule::FlipGroup.
      void FlipAtomGroup(const FlipGroup&);
      /** The list of FlipGroups.
      *
      */
      mutable list<FlipGroup> mvFlipGroup;
      
      // Group of atoms for random moves naturally respecting restraints
         /// List of StretchModeBondLength
         mutable list<StretchModeBondLength> mvStretchModeBondLength;
         /// List of StretchModeBondLength
         mutable list<StretchModeBondAngle> mvStretchModeBondAngle;
         /// List of StretchModeBondLength
         mutable list<StretchModeTorsion> mvStretchModeTorsion;
         /// List of StretchModeTwist
         mutable list<StretchModeTwist> mvStretchModeTwist;
      
      /// Groups of StretchMode not breaking any restraint (unless the one they are associated to)
      mutable std::list<StretchMode*> mvpStretchModeFree;
      /// Groups of StretchMode breaking restraints (beyond the one they are associated to)
      mutable std::list<StretchMode*> mvpStretchModeNotFree;
      /// Group of concurrent StretchModes (affecting common restraints)
      /// A given stretch mode can only belong to one group.
      struct StretchModeGroup
      {
         std::set<StretchMode*> mvpStretchMode;
         std::set<const MolBond*> mvpBrokenBond;
         std::set<const MolBondAngle*> mvpBrokenBondAngle;
         std::set<const MolDihedralAngle*> mvpBrokenDihedralAngle;
      };

   /// The current log(likelihood)
   mutable REAL mLogLikelihood;
   #ifdef __WX__CRYST__
   public:
      virtual WXCrystObjBasic* WXCreate(wxWindow*);
   #endif
};

}//namespace
#endif
