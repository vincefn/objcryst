#ifndef _OBJCRYST_CRYSTAL_H_
#define _OBJCRYST_CRYSTAL_H_

#include "CrystVector/CrystVector.h"

#include "ObjCryst/General.h"
#include "RefinableObj/RefinableObj.h"
#include "ObjCryst/SpaceGroup.h"
#include "ObjCryst/ScatteringPower.h"
#include "ObjCryst/Scatterer.h"

//#include <stdlib.h>
#include <string>
//#include <iomanip>
//#include <cmath>
//#include <typeinfo>
//#include <fstream>
//#include <ctime>

namespace ObjCryst
{
extern const RefParType *gpRefParTypeCrystal;
extern const RefParType *gpRefParTypeUnitCell;
extern const RefParType *gpRefParTypeUnitCellLength;
extern const RefParType *gpRefParTypeUnitCellAngle;

//######################################################################
//
//                CRYSTAL
///The crystal (Unit cell, spaceGroup, scatterers)
//######################################################################
class Crystal:public RefinableObj
{
   public:
      Crystal();
      
      /** \brief Crystal Constructor (orthorombic)
      *  \param a,b,c : unit cell dimension, in angstroems
      *  \param SpcGroup: space group number (1..230)
      *  \param name: name for the crystal, : '(TaSe4)2I'
      */
      Crystal(const double a, const double b, const double c,
              const string &SpaceGroupId);
      /** \brief Crystal Constructor (triclinic)
      *  \param a,b,c : unit cell dimension, in angstroems
      *  \param alpha,beta,gamma : unit cell angles
      *  \param SpcGroup: space group number (1..230)
      *  \param name: name for the crystal, : '(TaSe4)2I'
      */
      Crystal(const double a, const double b, const double c, const double alpha,
              const double beta, const double gamma,const string &SpaceGroupId);
              
      ///Crystal copy constructor
      Crystal(const Crystal &oldCryst);
      ///Crystal destructor
      ~Crystal();
      virtual const string GetClassName() const;
      
      /** \brief Add a scatterer to the crystal.
      *
      * \warning the scatterer \e must be allocated in the heap, since the scatterer
      * will not be copied but used directly. A Scatterer can only belong to one Crystal. It
      * will be detroyed when removed or when the Crystal is destroyed.
      * \param scatt : the address of the scatterer to be included in the crystal
      * scatterer names \b must be unique in a given crystal.
      *
      */
      void AddScatterer(Scatterer *scatt);
      /// Remove a Scatterer. This also deletes the scatterer.
      void RemoveScatterer(Scatterer *scatt);
      
      /// Number of scatterers in the crystal
      long GetNbScatterer()const;
      /** \brief Provides an access to the scatterers
      *
      * \param scattName the name of the scatterer to access
      */
      Scatterer & GetScatt(const string &scattName);
      /** \brief Provides a const access to the scatterers
      *
      * \param scattName the name of the scatterer to access
      */
      const Scatterer & GetScatt(const string &scattName) const;
      /** \brief Provides an access to the scatterers
      *
      * \param scattIndex the number of the scatterer to access
      */
      Scatterer & GetScatt(const long scattIndex);
      /** \brief Provides a const access to the scatterers
      *
      * \param scattIndex the number of the scatterer to access
      */
      const Scatterer & GetScatt(const long scattIndex) const;
      
      /// Get the registry of scatterers
      ObjRegistry<Scatterer>& GetScattererRegistry();
      /// Get the registry of ScatteringPower included in this Crystal.
      ObjRegistry<ScatteringPower>& GetScatteringPowerRegistry();
      /// Get the registry of ScatteringPower included in this Crystal.
      const ObjRegistry<ScatteringPower>& GetScatteringPowerRegistry()const;

      /// Add a ScatteringPower for this Crystal. It must be allocated in the heap,
      /// and not used by any other Crystal.
      void AddScatteringPower(ScatteringPower *scattPow);
      /// Remove a ScatteringPower for this Crystal. (the Scattering power is deleted).
      /// This function should check that it is not used any more before removing it.
      void RemoveScatteringPower(ScatteringPower *scattPow);
      /// Find a ScatteringPower from its name. Names must be unique in a given Crystal.
      ScatteringPower& GetScatteringPower(const string &name);
      
      /** \brief Get the list of all scattering components
      * \param useDynPopCorr: if true, then the dynamical polpulation correction
      * will be computed and stored in the list for each component.
      */
      virtual const ScatteringComponentList& GetScatteringComponentList()const;
      ///Lattice parameters (a,b,c,alpha,beta,gamma) as a 6-element vector in Angstroems 
      ///and radians.
      CrystVector_double GetLatticePar() const;
      ///Return one of the 6 Lattice parameters, 0<= whichPar <6 (a,b,c,alpha,beta,gamma),
      ///returned in Angstroems. 
      double GetLatticePar(const int whichPar)const;
      /** \brief Get the 'B' matrix (Crystal::mBMatrix)for the crystal (orthogonalization 
      *matrix for the given lattice, in the reciprocal space)
      *
      *The convention is taken following Giacovazzo, "Fundamentals of Crystallography", p.69
      *"e1 is chosen along a*, e2 in the (a*,b*) plane, then e3 is along c".
      */
      const CrystMatrix_double& GetBMatrix() const;
      /** \brief Get orthonormal cartesian coordinates for a set of (x,y,z)
      *fractional coordinates.
      *
      * Results are given in Amgstroems.
      *The convention is taken following :
      *e1 is chosen along a, e2 in the (a,b) plane, then e3 is along c*
      */
      CrystVector_double GetOrthonormalCoords(const double x,const double y,const double z) const;
      /** \brief Get orthonormal cartesian coordinates for a set of (x,y,z)
      *fractional coordinates.
      *
      * X,y and z input are changed to Amgstroems values
      *The convention is taken following :
      *e1 is chosen along a, e2 in the (a,b) plane, then e3 is along c*
      */
      void FractionalToOrthonormalCoords(double &x,double &y,double &z) const;
      /** \brief Get fractional cartesian coordinates for a set of (x,y,z)
      *orthonormal coordinates.
      *
      *Result is stored into x,y and z
      *The convention is taken following :
      *e1 is chosen along a, e2 in the (a,b) plane, then e3 is along c*
      */
      void OrthonormalToFractionalCoords(double &x,double &y,double &z) const;
      /// Prints some info about the crystal
      /// \todo one function to print on one line and a PrintLong() function
      void Print(ostream &os=cout) const;
            
      /// Access to the crystal SpaceGroup object
      const SpaceGroup & GetSpaceGroup()const;
      /// Access to the crystal SpaceGroup object
      SpaceGroup & GetSpaceGroup();
      
      /** \brief Minimum interatomic distance between all scattering components (atoms) in
      * the crystal.
      *
      * This will return a symetrical matrix with NbComp rows and cols, where 
      * NbComp is the number of independent scattering components in the unit cell. 
      * All distances are given in Angstroems.
      *
      * Note that the distance of a given atom with 'itself' is not generally equal
      * to 0 (except full special position), but equal to the min distance with its 
      * symetrics.
      * \param minDistance : atoms who are less distant than (minDistance,in Angstroems) 
      * are considered equivalent. So the smallest distance between any atoms will
      * be at least minDistance.
      */
      CrystMatrix_double GetMinDistanceTable(const double minDistance=0.1) const;
      /** \brief Print the minimum distance table between all scattering centers
      * (atoms) in the crystal.
      */
      void PrintMinDistanceTable(const double minDistance=0.1,ostream &os=cout) const;

      /** \brief Output POV-Ray Description for this Crystal
      *
      *
      */
      ostream& POVRayDescription(ostream &os,bool onlyIndependentAtoms=false)const;
      
      /** Create an OpenGL DisplayList of the crystal.
      *
      */
      virtual void GLInitDisplayList(const bool onlyIndependentAtoms=false,
                                     const double xMin=-.1,const double xMax=1.1,
                                     const double yMin=-.1,const double yMax=1.1,
                                     const double zMin=-.1,const double zMax=1.1)const;
      
      /** \brief Compute the 'Dynamical population correction for all atoms.
      *
      *\param overlapDist : distance below which atoms are considered overlapping and
      *should be corrected. 
      *\param mergeDist : distance below which atoms are considered fully overlapping.
      * If 3 atoms are 'fully' overlapping, then all have a dynamical population 
      * correction equal to 1/3
      *
      *\note This is const as long as ScatteringComponent::mDynPopCorr is mutable.
      */
      void CalcDynPopCorr(const double overlapDist=1., const double mergeDist=.0)const ;
      ///Reset Dynamical Population Correction factors
      void ResetDynPopCorr()const ;
      /// Set the use of dynamical population correction (Crystal::mUseDynPopCorr).
      /// This \e seriously affects the speed of the calculation, since computing
      /// distances is lenghty.
      void SetUseDynPopCorr(const int use);
      /// Get the Anti-bumping/pro-Merging cost function
      double GetBumpMergeCostFunction() const;
      /// Set the Anti-bumping distance between two scattering types
      void SetBumpMergeDistance(const ScatteringPower &scatt1,
                                const ScatteringPower &scatt2, const double dist=1.5);
      /// Set the Anti-bumping distance between two scattering types.
      void SetBumpMergeDistance(const ScatteringPower &scatt1,
                                const ScatteringPower &scatt2, const double dist,
                                const bool allowMerge);
      /// When were lattice parameters changed ?
      const RefinableObjClock& GetClockLatticePar()const;
      /// When was the list of scatterers last changed ?
      const RefinableObjClock& GetClockScattererList()const;
      //Cost functions
         unsigned int GetNbCostFunction()const;
         const string& GetCostFunctionName(const unsigned int)const;
         const string& GetCostFunctionDescription(const unsigned int)const;
         virtual double GetCostFunctionValue(const unsigned int);
         
      virtual void Output(ostream &os,int indent=0)const;
      virtual void Input(istream &is,const XMLCrystTag &tag);
      virtual void InputOld(istream &is,const IOCrystTag &tag);
      
      virtual void GlobalOptRandomMove(const double mutationAmplitude);
      /** \brief output Crystal structure as a cif file (crude !)
      *
      * This is very crude so far: only isotropic scattering power
      * are supported, and there is not much information beside atom
      * positions... Still a lot to do !!!
      */
      virtual void OutputCIF(ostream &os)const;
      
   protected:
   
   private:
      /** \brief Init all Crystal parameters
      *  \param a,b,c : unit cell dimension, in angstroems
      *  \param alpha,beta,gamma : unit cell angles
      *  \param SpcGroup: space group number (1..230)
      *  \param name: name for the crystal, : '(TaSe4)2I'
      */
      void Init(const double a, const double b, const double c, const double alpha,
                const double beta, const double gamma,const string &SpaceGroupId,
                const string& name);
      
      /// Find a scatterer (its index # in mpScatterrer[]) with a given name
      /// \warning There should be no duplicate names !!! :TODO: test in AddScatterer()
      int FindScatterer(const string &scattName)const;
      
      /// Init the UBMatrix and EuclMatrix. They are re-computed only if parameters
      ///have changed since last call.
      
      void InitMatrices() const;
      /// \internal
      /// Update cell parameters for tetragonal, trigonal, hexagonal, cubic lattices.
      ///Also set angular parameters for those group which need it. This is needed during
      ///Refinement, since for example in a quadratic spg, only a is refined and
      ///we need to have b=a updated very often...
      void UpdateLatticePar();

      /** \brief Prepare the refinable parameters list
      *
      *This is called once when creating the crystal. Scatterers parameters
      *are added by  Crystal::AddScatterer()
      */
      void InitRefParList();
      
      /** \brief Compute the distance Table (mDistTable) for all scattering components
      * \param fast : if true, the distance calculations will be made using
      * integers, thus with a lower precision but faster. Less atoms will also
      * be involved (using the AsymmetricUnit) to make it even faster.
      * \param asymUnitMargin (in Angstroem). This is used only if fast=true.
      * In that case, the distance is calculated between (i) independent atoms in the
      * asymmetric unit cell and (ii) all atoms which are inside the asymmetric unit cell
      * or less than 'asymUnitMargin' distant from the asymmetric unit borders.
      * This parameter should be used when only the shortest distances need to be calculated
      * (typically for dynamical population correction). Using a too short  margin will
      * result in having some distances calculated wrongly (ie one atom1 in the unit cell
      * could have an atom2 neighbor just outside the sym unit: if margin=0, then the distance
      * is calculated between atom1 and the atom2 symetric inside the asym unit).
      * \warning Crystal::GetScatteringComponentList() \b must be called beforehand.
      */
      void CalcDistTable(const bool fast,const double asymUnitMargin=4)const;
            
      ///a,b and c in Angstroems, angles (stored) in radians
      ///For cubic, rhomboedric crystals, only the 'a' parameter is relevant.
      ///For quadratic and hexagonal crystals, only a and c parameters are relevant.
      /// The MUTABLE is temporary ! It should not be !
      CrystVector_double mCellDim;
      /// The space group of the crystal
      SpaceGroup mSpaceGroup ;
      /// The registry of scatterers for this crystal
      ObjRegistry<Scatterer> mScattererRegistry ;
      
      /** \brief B Matrix (Orthogonalization matrix for reciprocal space)
      * \f[ B= \left[ \begin {array}{ccc} a^* & b^*\cos(\gamma^*) & c^*\cos(\beta^*) \\
      *                            0 & b^*\sin(\gamma^*) & -c^*\sin(\beta^*)\cos(\alpha) \\
      *                            0 & 0 & \frac{1}{c}\end{array} \right]\f]
      *\f[ \left[ \begin{array}{c} k_x \\ k_y \\ k_z \end{array} \right]_{orthonormal}  
      *  = B \times \left[ \begin{array}{c} h \\ k \\ l \end{array}\right]_{integer} \f]
      *\note this matrix is and must remain upper triangular. this is assumed for
      *some optimizations.
      */
      mutable CrystMatrix_double mBMatrix;
      /** \brief Eucl Matrix (Orthogonalization matrix for direct space)
      * \f[ M_{orth}= \left[ \begin {array}{ccc} a & b\cos(\gamma) & c\cos(\beta) \\
      *                            0 & b\sin(\gamma) & -c\sin(\beta)\cos(\alpha^*) \\
      *                            0 & 0 & \frac{1}{c^*}\end{array} \right]\f]
      *\f[ \left[ \begin{array}{c} x \\ y \\ z \end{array} \right]_{orthonormal}  
      *  = M_{orth} \times \left[ \begin{array}{c} x \\ y \\ z \end{array}
      *                    \right]_{fractional coords}\f]
      *\note this matrix is and must remain upper triangular. this is assumed for
      *some optimizations.
      */
      mutable CrystMatrix_double mOrthMatrix;
      /// inverse of Eucl Matrix (de-orthogonalization matrix for direct space)
      mutable CrystMatrix_double mOrthMatrixInvert;
            
      /** \brief Distance table (squared) between all scattering components in the crystal
      *
      *Symetrical matrix, in Angstroems^2 (the square root is not computed
      *for optimization purposes).
      */
      mutable CrystMatrix_double mDistTableSq;
      /** \brief Index of scattering components for the Distance table
      *
      *These are the index of the scattering components corresponding to each row/column in
      *the distance table. Each component has as many entries as its number of symetrics.
      *\note (kludge) this will only be valid if the order of components does not change...
      */
      mutable CrystVector_long  mDistTableIndex;
      /// The list of all scattering components in the crystal
      mutable ScatteringComponentList mScattCompList;
         
      /// Clock for lattice paramaters.
      RefinableObjClock mLatticeClock;
      ///Use Dynamical population correction (ScatteringComponent::mDynPopCorr) during Structure
      ///factor calculation ?
      RefObjOpt mUseDynPopCorr;
      /// Matrix of "bumping" (squared) distances
      CrystMatrix_double mBumpDistanceMatrix;
      /// Allow merging of atoms in the bump/merge function(should be true for identical atoms)?
      CrystMatrix_bool mAllowMerge;
      
      /// The registry of ScatteringPower for this Crystal.
      ObjRegistry<ScatteringPower> mScatteringPowerRegistry;
      
      //Clocks
      /// Last time the list of Scatterers was changed
      RefinableObjClock mClockScattererList;
      /// Last time lattice parameters were changed
      RefinableObjClock mClockLatticePar;
      /// \internal Last time metric matrices were computed
      mutable RefinableObjClock mClockMetricMatrix;
      /// \internal Last time the ScatteringComponentList was generated
      mutable RefinableObjClock mClockScattCompList;
      /// \internal Last time the Neighbor Table was generated
      mutable RefinableObjClock mClockNeighborTable;
      /// \internal Last time the lattice parameters whre updated
      mutable RefinableObjClock mClockLatticeParUpdate;
      /// \internal Last time the dynamical population correction was computed
      mutable RefinableObjClock mClockDynPopCorr;
   #ifdef __WX__CRYST__
   public:
      virtual WXCrystObjBasic* WXCreate(wxWindow*);
      friend class WXCrystal;
   #endif
};

/// Global registry for all Crystal objects
extern ObjRegistry<Crystal> gCrystalRegistry;


}//namespace


#endif //_OBJCRYST_CRYSTAL_H_
