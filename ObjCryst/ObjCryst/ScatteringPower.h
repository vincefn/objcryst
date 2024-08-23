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
#ifndef _OBJCRYST_SCATTPOWER_H_
#define _OBJCRYST_SCATTPOWER_H_

#include "ObjCryst/CrystVector/CrystVector.h"
#include "ObjCryst/ObjCryst/General.h"
#include "ObjCryst/RefinableObj/RefinableObj.h"


//#include <stdlib.h>
//#include <string>
//#include <iomanip>
//#include <cmath>
//#include <typeinfo>
//#include <fstream>
//#include <ctime>

// forward declaration to avoid including boost headers
namespace cctbx { namespace eltbx { namespace xray_scattering {class gaussian;}}}

namespace ObjCryst
{
extern const RefParType *gpRefParTypeScattPow;
extern const RefParType *gpRefParTypeScattPowResonant;
extern const RefParType *gpRefParTypeScattPowTemperature;
extern const RefParType *gpRefParTypeScattPowTemperatureIso;
extern const RefParType *gpRefParTypeScattPowTemperatureAniso;
class NiftyStaticGlobalObjectsInitializer_ScatteringPower
{
   public:
      NiftyStaticGlobalObjectsInitializer_ScatteringPower()
      {
         if (mCount++ == 0)
         {
            gpRefParTypeScattPow=new RefParType(gpRefParTypeObjCryst,"Scattering Power");
            gpRefParTypeScattPowResonant=new RefParType(gpRefParTypeScattPow,"Resonant Scatt.");
            gpRefParTypeScattPowTemperature=new RefParType(gpRefParTypeScattPow,"Temperature");
            gpRefParTypeScattPowTemperatureIso=new RefParType(gpRefParTypeScattPowTemperature,"Isotropic");
            gpRefParTypeScattPowTemperatureAniso=new RefParType(gpRefParTypeScattPowTemperatureIso,"Anisotropic");
         }
      }
      ~NiftyStaticGlobalObjectsInitializer_ScatteringPower()
      {
         if (--mCount == 0)
         {
            delete gpRefParTypeScattPow;
            delete gpRefParTypeScattPowResonant;
            delete gpRefParTypeScattPowTemperature;
            delete gpRefParTypeScattPowTemperatureIso;
            delete gpRefParTypeScattPowTemperatureAniso;
            gpRefParTypeScattPow=0;
            gpRefParTypeScattPowResonant=0;
            gpRefParTypeScattPowTemperature=0;
            gpRefParTypeScattPowTemperatureIso=0;
            gpRefParTypeScattPowTemperatureAniso=0;
         }
      }
   private:
      static long mCount;
};
static NiftyStaticGlobalObjectsInitializer_ScatteringPower NiftyStaticGlobalObjectsInitializer_ScatteringPower_counter;

   class ScatteringData;//forward declaration :KLUDGE: ?
//######################################################################
//
//      SCATTERING POWER
/** \brief Abstract Base Class to describe the scattering power of any
* Scatterer component in a crystal.
*
* This includes:
* - the scattering factor,
* - the temperature factor
* - real and imaginary parts of the resonant
* scattering factor.
*
* The interface is independent of the radiation type.
*
* This base class is designed to handle both isotropic and anisotropic
* versions of scattering, temperature and anomalous factors.
*
* \todo Anisotropic scattering (temperature factor especially) code, using derived
* classes
* \todo Clarify organization by removing any 'real' data from the top, abstract
* base class (eg remove Biso and Betaij), and by creating derived classes.
* Optionnaly 3 classes (used as members of ScatteringPower) could be created,
* TemperatureFactor, ScatteringFactor, and ResonantScatteringFactor. In any way
* the design of this class should not evolve, so that code using the ScatteringPower
* interface will remain compatible whatever modifications are made.
* \warning: there is currently a storage for Anisotropic Displacement Parameters,
* but Debye-Waller calculation is \e only isotropic.
*/
//######################################################################
class ScatteringPower:virtual public RefinableObj
{
   public:
      ScatteringPower();
      ScatteringPower(const ScatteringPower& old);
      virtual ~ScatteringPower();
      virtual const string& GetClassName() const;
      virtual void operator=(const ScatteringPower& rhs);
      /** Comparison operator. Two scattering powers are equal if they have the same
      * displacement parameter, correspond to same element, and are of the same class.
      */
      virtual bool operator==(const ScatteringPower& rhs) const;
      /** Comparison operator. Two scattering powers are equal if they have the same
      * displacement parameter, correspond to same element, and are of the same class.
      */
      virtual bool operator!=(const ScatteringPower& rhs) const;
      /** \brief Get the Scattering factor for all reflections of a given
      * ScatteringData object.
      *
      * \return a vector with the scattering factor for all reflections, in the same
      * order as in the ScatteringData object. This format is independent of the radiation
      * type (X-Ray, neutron..).
      * \param data: the ScatteringData object, giving access to all the reflections.
      * \param spgSymPosIndex: if the ScatteringPower is anisotropic, then the
      * different symmetrics will not have the same scattering power for all reflections.
      * This parameter is the index of the symmetric position in the Spacegroup.
      * If spgSymPosIndex=-1, the isotropic values are returned.
      * \warning There is no anisotropic code yet, so spgSymPosIndex is simply ignored so far
      * , but the design of this function is general for any anisotropic scattering.
      */
      virtual CrystVector_REAL GetScatteringFactor(const ScatteringData &data,
                                                     const int spgSymPosIndex=-1) const=0;
      /// Get the scattering factor at (0,0,0). Used for scatterer (electron, nucleus)
      /// density generation.
      virtual REAL GetForwardScatteringFactor(const RadiationType) const=0;
      /** \brief Get the temperature factor for all reflections of a given
      * ScatteringData object.
      *
      * \return a vector with the temperature factor for all reflections, in the same
      * order as in the ScatteringData object.
      * \param data: the ScatteringData object, giving access to all the reflections.
      * \param spgSymPosIndex: if the ScatteringPower is anisotropic, then the
      * different symmetrics will not have the same scattering power for all reflections.
      * This parameter is the index of the symmetric position in the Spacegroup.
      * If spgSymPosIndex=-1, the isotropic values are returned.
      * \warning There is no anisotropic code yet, so spgSymPosIndex is simply ignored so far
      * , but the design of this function is general for any anisotropic scattering.
      */
      virtual CrystVector_REAL GetTemperatureFactor(const ScatteringData &data,
                                                      const int spgSymPosIndex=-1) const=0;
      /** \brief Get the real part of the resonant scattering factor.
      *
      * \return a matrix where each row corresponds to each wavelength (currently only
      * monochromatic experiments are made so there is only one row), and each column
      * corresponds to each reflection \e only if the scattering term is anisotropic, which
      * is not the case so far...
      * \param data: the ScatteringData object, giving access to all the reflections and
      * a list of wavelengths).
      * \param spgSymPosIndex: if the ScatteringPower is anisotropic, then the
      * different symmetrics will not have the same scattering power for all reflections.
      * This parameter is the index of the symmetric position in the Spacegroup.
      * If spgSymPosIndex=-1, the isotropic values are returned.
      * \warning There is no anisotropic code yet, so spgSymPosIndex is simply ignored so far
      * , but the design of this function is general for any anisotropic scattering.
      */
      virtual CrystMatrix_REAL GetResonantScattFactReal(const ScatteringData &data,
                                                          const int spgSymPosIndex=-1) const=0;
      /** \brief Get the imaginary part of the resonant scattering factor.
      *
      * \return a matrix where each row corresponds to each wavelength (currently only
      * monochromatic experiments are made so there is only one row), and each column
      * corresponds to each reflection \e only if the scattering term is anisotropic, which
      * is not the case so far...
      * \param data: the ScatteringData object, giving access to all the reflections,
      * and a list of wavelengths.
      * \param spgSymPosIndex: if the ScatteringPower is anisotropic, then the
      * different symmetrics will not have the same scattering power for all reflections.
      * This parameter is the index of the symmetric position in the Spacegroup.
      * If spgSymPosIndex=-1, the isotropic values are returned.
      * \warning There is no anisotropic code yet, so spgSymPosIndex is simply ignored so far
      * , but the design of this function is general for any anisotropic scattering.
      */
      virtual CrystMatrix_REAL GetResonantScattFactImag(const ScatteringData &data,
                                                          const int spgSymPosIndex=-1) const=0;
      /// Is the scattering factor anisotropic ?
      virtual bool IsScatteringFactorAnisotropic()const;
      /// Is the thermic factor anisotropic ?
      virtual bool IsTemperatureFactorAnisotropic()const;
      /// Are the resonant scattering terms anisotropic ?
      virtual bool IsResonantScatteringAnisotropic()const;
      /** \brief Symbol for this Scattering power (the atom name for atoms)
      *
      */
      virtual const string& GetSymbol() const;
      /** \brief Returns the isotropic temperature B factor
      *
      */
      REAL GetBiso() const;
      /** \brief Returns the isotropic temperature B factor
      *
      */
      REAL& GetBiso();
      /** \brief Sets the isotropic temperature B factor
      *
      */
      virtual void SetBiso(const REAL newB);
      /** \brief Returns the anisotropic temperature B factor for (i, j) pair.
      *
      * \warning: this is ambiguous, as it is Beta_ij which are stored, and not Bij...
      */
      REAL GetBij(const size_t &i, const size_t &j) const;
      /** \brief Returns the anisotropic temperature B factor for given index.
      *
      * 0 -> (1, 1)
      * 1 -> (2, 2)
      * 2 -> (3, 3)
      * 3 -> (1, 2)
      * 4 -> (1, 3)
      * 5 -> (2, 3)
      *
      * \warning: this is ambiguous, as it is Beta_ij which are stored, and not Bij...
      */
      REAL GetBij(const size_t &idx) const;
      /** \brief Sets the anisotropic temperature B factor for (i, j) pair.
      *
      * \warning: this is ambiguous, as it is Beta_ij which are stored, and not Bij...
      */
      virtual void SetBij(const size_t &i, const size_t &j, const REAL newB);
      /** \brief Sets the anisotropic temperature B factor for given index.
      *
      * 0 -> (1, 1)
      * 1 -> (2, 2)
      * 2 -> (3, 3)
      * 3 -> (1, 2)
      * 4 -> (1, 3)
      * 5 -> (2, 3)
      *
      * \warning: this is ambiguous, as it is Beta_ij which are stored, and not Bij...
      */
      virtual void SetBij(const size_t &idx, const REAL newB);
      /** \brief Returns true if the scattering power is isotropic, else false.
      *
      *
      */
      bool IsIsotropic()const ;
      /// Get the number identifying this kind of scatterer, used to decide whether two
      /// scatterers are equivalent, for the dynamical occupancy correction.
      long GetDynPopCorrIndex()const;
      /// Total number of ScatteringPower object
      long GetNbScatteringPower()const;
      /// ObjCrystClock time when the last modification was made to the object
      const RefinableObjClock& GetLastChangeClock()const;
      ///Get the (POV-Ray) name associated to the color (if any)
      const string& GetColourName()const;
      /// Get the float[3] array of RGB components defining the colour of this scattering power
      const float* GetColourRGB()const;
      /// Set the colour from the associated POV-Ray name
      void SetColour(const string& colorName);
      /// Set the colour from RGB components (all between 0 and 1.)
      void SetColour(const float r,const float g,const float b);
      /// Return the physical radius of this type of scatterer (for 3D display purposes).
      /// \warning this may be removed later.
      virtual REAL GetRadius()const=0;
      virtual void GetGeneGroup(const RefinableObj &obj,
                                CrystVector_uint & groupIndex,
                                unsigned int &firstGroup) const;
      /// Maximum Likelihood: get the estimated error (sigma) on the positions
      /// for this kind of element.
      REAL GetMaximumLikelihoodPositionError()const;
      /// Maximum Likelihood: set the estimated error (sigma) on the positions
      /// for this kind of element.
      void SetMaximumLikelihoodPositionError(const REAL mle);
      /// Maximum Likelihood: get the number of ghost elements per asymmetric unit.
      REAL GetMaximumLikelihoodNbGhostAtom()const;
      /// Maximum Likelihood: set the number of ghost elements per asymmetric unit.
      void SetMaximumLikelihoodNbGhostAtom(const REAL nb);
      /// Get the clock value for the last change on the maximum likelihood parameters
      /// (positionnal error, number of ghost atoms).
      const RefinableObjClock& GetMaximumLikelihoodParClock()const;
      virtual REAL GetFormalCharge()const;
      virtual void SetFormalCharge(const REAL charge);
   protected:
      virtual void InitRefParList()=0;
      /// Initialization of the object, used by all constructors, and operator=.
      virtual void Init();
      ///Get RGB Colour coordinates from Colour Name
      virtual void InitRGBColour();
      /// number identifying this kind of scatterer, for the dynamical occupancy correction.
      /// Right now it is the atomic number.
      long mDynPopCorrIndex;
      /// Temperature isotropic B factor.
      REAL mBiso;
      /// Is the scattering isotropic ?
      bool mIsIsotropic;
      /** Anisotropic Beta(ij)
      *
      * \internal
      * These are stored temporarily, and derived from the Bij
      */
      mutable CrystVector_REAL mBeta;
      /// Anisotropic B(ij)
      CrystVector_REAL mB;
      /// Clock.
      RefinableObjClock mClock;
      /// Colour for this ScatteringPower (from POVRay)
      string mColourName;
      /// Colour for this ScatteringPower using RGB
      float mColourRGB[3];
      // Maximum Likelihood
         /// estimated error (sigma) on the positions for this type of element.
         REAL mMaximumLikelihoodPositionError;
         ///
         RefinableObjClock mMaximumLikelihoodParClock;
         /// Number of ghost atoms in the asymmetric unit.
         /// These contribute to the variance of the structure factor, but not to the structure
         /// factor as the uncertainty on their position is infinite.
         REAL mMaximumLikelihoodNbGhost;
      /** Formal Charge. This can be used for bond valence analysis,
      * or energy calculations.
      *
      * Default value is 0.
      */
      REAL mFormalCharge;
};

/// Global registry for all ScatteringPower objects
extern ObjRegistry<ScatteringPower> gScatteringPowerRegistry;

//######################################################################
//
//      SCATTERING POWER ATOM
/** \brief The Scattering Power for an Atom
*
*/
//######################################################################

class ScatteringPowerAtom:virtual public ScatteringPower
{
   public:
      ScatteringPowerAtom();
      /**   \brief Atom constructor
      *  \param symbol : 'Ti' , 'Ti4+', 'Cl1-'
      *These symbols \e must correspond to one of the entries of the international tables
      *for crystallography (1995) giving the analytical approximation for scattering factors.
      *  \param name : name of the atom ('Ta1','Sm2', 'Tungsten_1'...).
      *The name can have \e any format but spaces should be avoided, since it
      *will generate problems when reading the names from a file...
      *  \param biso : Isotropic thermic coefficient
      */
      ScatteringPowerAtom(const string &name,const string &symbol,const REAL bIso=1.0);
      ScatteringPowerAtom(const ScatteringPowerAtom& old);
      ~ScatteringPowerAtom();
      virtual const string& GetClassName() const;
      /// Re-initialize parameters (after using the default constructor).
      void Init(const string &name,const string &symbol,const REAL bIso=1.0);
      virtual CrystVector_REAL GetScatteringFactor(const ScatteringData &data,
                                                     const int spgSymPosIndex=0) const;
      virtual REAL GetForwardScatteringFactor(const RadiationType) const;
      virtual CrystVector_REAL GetTemperatureFactor(const ScatteringData &data,
                                                     const int spgSymPosIndex=0) const;
      virtual CrystMatrix_REAL GetResonantScattFactReal(const ScatteringData &data,
                                                     const int spgSymPosIndex=0) const;
      virtual CrystMatrix_REAL GetResonantScattFactImag(const ScatteringData &data,
                                                     const int spgSymPosIndex=0) const;
      /// Set the symbol for this atom
      void SetSymbol(const string &symbol) ;
      /// Returns the symbol ('Ta', 'O2-',...) of the atom
      virtual const string& GetSymbol() const;
      /** \brief Returns the standard name of the element (ie "hydrogen", "tantalum",..).
      *
      *  Names are extracted form Grosse-Kunstleve 'atominfo' package,
      *which uses data from the CRC Handbook of Chemistry & Physics, 63rd & 70th editions
      */
      string GetElementName() const;
      /// Atomic number for this atom
      int GetAtomicNumber() const;
      /// Atomic weight (g/mol) for this atom
      REAL GetAtomicWeight() const;
      /// Atomic radius for this atom or ion, in Angstroems (ICSD table from cctbx)
      REAL GetRadius() const;
      /// Covalent Radius for this atom, in Angstroems (from cctbx)
      REAL GetCovalentRadius() const;
      /// Maximum number of covalent bonds (from openbabel element.txt)
      unsigned int GetMaxCovBonds() const;
      virtual void Print()const;
      virtual void XMLOutput(ostream &os,int indent=0)const;
      virtual void XMLInput(istream &is,const XMLCrystTag &tag);
      //virtual void XMLInputOld(istream &is,const IOCrystTag &tag);
   protected:
      /** \internal
      * Fetch the coefficients for analytical approximation of the
      *atomic scattering factor.
      */
      void InitAtScattCoeffsWK95();
      /**   \internal
      * Fetch the coefficients neutron scattering.
      */
      void InitAtNeutronScattCoeffs();
      virtual void InitRefParList();
      /** \brief Symbol of this atom
      *
      *This symbol \e must correspond to one of the entries of the international tables
      *for crystallography (1995) giving the analytical approximation for scattering factors.
      */
      string mSymbol;
      /// atomic number (Z) for the atom
      int mAtomicNumber;
      /// atomic weight (g/mol) for the atom
      REAL mAtomicWeight;

      /** Pointer to cctbx's gaussian describing the thomson x-ray
      * scattering factor.
      */
      cctbx::eltbx::xray_scattering::gaussian *mpGaussian;

      /** \brief Neutron Bond Coherent Scattering lengths
      *
      *Real and imaginary (for atoms who have an imaginary part)
      *
      *Reference : Neutron News, Vol. 3, No. 3, 1992, pp. 29-37.
      */
      REAL mNeutronScattLengthReal,mNeutronScattLengthImag;

      /// Radius of the atom or ion, in Angstroems (ICSD table from cctbx)
      REAL mRadius;
      /// Covalent Radius for this atom, in Angstroems (from cctbx)
      REAL mCovalentRadius;
      /// Maximum number of covalent bonds
      unsigned int mMaxCovBonds;

      /** \brief Neutron Absorption cross section (barn)
      *
      *For 2200 m/s neutrons.
      *
      *Reference : Neutron News, Vol. 3, No. 3, 1992, pp. 29-37.
      */
      REAL mNeutronAbsCrossSection;
   private:
      // Avoid compiler warnings.  Explicitly hide the base-class method.
      void Init();
   #ifdef __WX__CRYST__
   public:
      virtual WXCrystObjBasic* WXCreate(wxWindow*);
   friend class WXScatteringPowerAtom;
   #endif
};
/// Global registry for all ScatteringPowerAtom objects
extern ObjRegistry<ScatteringPowerAtom> gScatteringPowerAtomRegistry;

//######################################################################
//
//      SCATTERING COMPONENT
/** \brief A scattering position in a crystal, associated with the corresponding
* occupancy and a pointer to the ScatteringPower. Also given is the
*
*/
//######################################################################
struct ScatteringComponent
{
   // Default Constructor
   ScatteringComponent();
   bool operator==(const ScatteringComponent& rhs)const;
   bool operator!=(const ScatteringComponent& rhs)const;
   ///Print  one line oabout this component
   void Print() const;
   /// Coordinates of scattering positions i the crystal with the corresponding occupancy
   REAL mX,mY,mZ,mOccupancy;
   /// The ScatteringPower associated with this position
   const ScatteringPower *mpScattPow;

   /// Dynamical Population Correction.
   ///
   /// The population of any atom is given by mOccupancy*mDynPopCorr.
   /// mPopu is the \e real mOccupancy (0<.<1), and should be the only one
   /// used during a refinement. However during a \e model \e search for the structure,
   /// atoms may fall unexpectedly in a special position or with an overlap of
   /// two atoms (the shared oxygen between two polyhedras, for example). In these
   /// cases it is necessary to dynamically correct the population during the
   /// generation of structural models.
   /// See also Crystal::CalcDynPopCorr
   ///
   /// \note this parameter is mutable, and is computed by the Crystal object
   mutable REAL mDynPopCorr;
   // The scatterer to which this component is associated
   // Scatterer *mScatterer;
};

//######################################################################
//
//      SCATTERING COMPONENT LIST
/** \brief list of scattering positions in a crystal, associated with the corresponding
* occupancy and a pointer to the ScatteringPower.
*
*/
//######################################################################
class ScatteringComponentList
{
   public:
      ScatteringComponentList();
      ScatteringComponentList(const long nbComponent);
      ScatteringComponentList(const ScatteringComponentList &old);
      ~ScatteringComponentList();
      /// Reset the list. This does \b not free the memory, but simply forgets
      /// that there already are some entries.
      void Reset();
      /// Access to a component
      const ScatteringComponent& operator()(const long i) const;
      ScatteringComponent& operator()(const long i);
      /// Number of components
      long GetNbComponent() const;
      /// Assignement operator
      void operator=(const ScatteringComponentList &rhs);
      /// Compare two lists.
      bool operator==(const ScatteringComponentList &rhs)const;
      /// Add another list of components
      void operator+=(const ScatteringComponentList &rhs);
      /// Add component
      void operator+=(const ScatteringComponent &rhs);
      /// Add component (the whole list should be updated after that)
      void operator++();
      /// Remove component (the whole list should be updated after that)
      void operator--();
      /// Print the list of Scattering components. For debugging
      void Print() const;
   protected:
      /// The vector of components
      vector<ScatteringComponent>  mvScattComp;
};

}//namespace
#include "ObjCryst/ObjCryst/ScatteringData.h"

#endif //_OBJCRYST_SCATTPOWER_H_
