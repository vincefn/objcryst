/* 
* ObjCryst++ : a Crystallographic computing library in C++
*			http://objcryst.sourceforge.net
*			http://www.ccp14.ac.uk/ccp/web-mirrors/objcryst/
*
*  (c) 2000-2001 Vincent FAVRE-NICOLIN vincefn@users.sourceforge.net
*
*/
/*   Atom.h
*  header file for the Atom scatterer
*
*/
#ifndef _OBJCRYST_ATOM_H_
#define _OBJCRYST_ATOM_H_

#include "CrystVector/CrystVector.h"

#include "ObjCryst/General.h"

#include "ObjCryst/ScatteringPower.h"
#include "ObjCryst/Scatterer.h"

//#include <stdlib.h>
#include <string>
//#include <iomanip>
//#include <cmath>
//#include <typeinfo>
//#include <fstream>
//#include <ctime>
//using namespace std;

namespace ObjCryst
{

//######################################################################
//
///    ATOM : the basic atom, within the crystal.
///
/// This class records the position of the atom, and has a pointer to its
/// ScatteringPowerAtom.
///
/// Note that there can be 'Dummy' atoms, for which the used symbol is "X",
/// and which have no scattering power (use with caution: dummy atoms
/// are only supposed to be used within ZScatterer)
//######################################################################

class Atom: public Scatterer
{
   public:
      /// Default constructor. The Atom \e must be initialized thereafter
		/// using Atom::Init()
      Atom();
      /**   \brief Atom constructor
      *  \param x,y,z : \e fractional coordinates of the atom
      *  \param pow : the ScatteringPower associated to this atom. Must be allocated separately.
      *  \param name : name of the atom ('Ta1','Sm2', 'Tungsten_1'...).
      * The name can have \e any format but spaces should be avoided, since it
      * will generate problems when reading the names from a file...
      */
      Atom( const double x, const double y, const double z,
            const string &name, const ScatteringPowerAtom *pow);
      /**   \brief Atom constructor
      *  \param x,y,z : \e fractional coordinates of the atom
      *  \param popu : the population of the atom (0.0->1.0)
      * This should take into account the multiplicity of the atom. For
      * an atom in group P2 and on the 2 axis, this should be set to 0.5,
		* \b unless you are using the dynamical occupancy correction (recommended
		* for global optimizations). See Crystal::CalcDynPopCorr() and
		* Crystal::mUseDynPopCorr
		*
      *  \param pow : the ScatteringPower associated to this atom. Must be allocated separatly.
      *  \param name : name of the atom ('Ta1','Sm2', 'Tungsten_1'...).
      * The name can have \e any format but spaces should be avoided (just a precaution)
      */
      Atom( const double x, const double y, const double z,const string &name,
             const ScatteringPowerAtom *pow, const double popu);
      /// Copy constructor
      Atom(const Atom &old);
      /// \internal so-called Virtual copy constructor, needed to make copies
      /// of arrays of Scatterers
      virtual Atom* CreateCopy() const;
      /// Atom desintegrator...
     ~Atom();
      virtual const string GetClassName() const;
      ///
      virtual void operator=(const Atom & rhs);

      /** initialize the atom (used for arrays of atoms).
      *  \param x,y,z : \e fractional coordinates of the atom
      *  \param pow : the ScatteringPower associated to this atom. Must be allocated separately.
      *  \param name : name of the atom ('Ta1','Sm2', 'Tungsten_1'...).
		*/
      void Init(const double x, const double y, const double z,
            const string &name, const ScatteringPowerAtom *pow,
            const double popu=1);
      
      virtual int GetNbComponent() const;
      virtual const ScatteringComponentList& GetScatteringComponentList() const;
      virtual string GetComponentName(const int i) const;

      virtual void Print() const;
      /// \internal This should be (and soon will be) private.
      virtual void Update() const;

      /** \brief Returns the molar mass of the atom.
      *
      *  Values are extracted from the 'atominfo' package,
      * which uses data from the CRC Handbook of Chemistry & Physics, 63rd & 70th editions
      * The Mass is actually extracted from the ScatteringPowerAtom.
      */
      double GetMass() const;
      /** \brief Returns the radius (in Angstroems) of the atom.
      *
      *  Values are extracted from the 'atominfo' package,
      * which uses data from the ICSD/CRYSTIN Manual
      * The Radius is extracted from the ScatteringPowerAtom.
      */
      double GetRadius() const;
      /** \brief Output a description of the scatterer for POVRay
      *
      */
      virtual ostream& POVRayDescription(ostream &os,const Crystal &cryst,
                                         bool onlyIndependentAtoms=false)const;
      /** Create the OpenGL Display List for this Atom.
      *
      */
      virtual void GLInitDisplayList(const Crystal &cryst,
                                     const bool onlyIndependentAtoms=false,
                                     const double xMin=-.1,const double xMax=1.1,
                                     const double yMin=-.1,const double yMax=1.1,
                                     const double zMin=-.1,const double zMax=1.1)const;

      /// Is this a dummy atom ? (ie no ScatteringPower)
      /// Dummy atoms should not exist !
      bool IsDummy()const;

      virtual void Output(ostream &os,int indent=0)const;
      virtual void Input(istream &is,const XMLCrystTag &tag);
      virtual void InputOld(istream &is,const IOCrystTag &tag);
		/// Get the ScatteringPowerAtom corresponding to this atom.
      const ScatteringPowerAtom& GetScatteringPower()const;
   protected:
   private:
      /// Prepare refinable parameters for the scatterer object
      virtual void InitRefParList();
   
      /// The list of scattering components.
      mutable ScatteringComponentList mScattCompList;
      /// The ScatteringPowerAtom associated to that atom
      const ScatteringPowerAtom *mpScattPowAtom;
   #ifdef __WX__CRYST__
   public:
      virtual WXCrystObjBasic* WXCreate(wxWindow*);
   friend class WXAtom;
   #endif
};

}//namespace Objcryst

// do we need this ?
#include "ObjCryst/Crystal.h"

#endif //_OBJCRYST_ATOM_H_
