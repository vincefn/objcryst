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
#ifndef _OBJCRYST_SCATTERER_H_
#define _OBJCRYST_SCATTERER_H_

#include "ObjCryst/CrystVector/CrystVector.h"

#include "ObjCryst/ObjCryst/General.h"

#include "ObjCryst/ObjCryst/ScatteringPower.h"

#include <string>

namespace ObjCryst
{
class Crystal ; //forward declaration.:KLUDGE: ?
extern const RefParType *gpRefParTypeScatt;
extern const RefParType *gpRefParTypeScattTransl;
extern const RefParType *gpRefParTypeScattTranslX;
extern const RefParType *gpRefParTypeScattTranslY;
extern const RefParType *gpRefParTypeScattTranslZ;
extern const RefParType *gpRefParTypeScattOrient;
extern const RefParType *gpRefParTypeScattConform;
extern const RefParType *gpRefParTypeScattConformBondLength;
extern const RefParType *gpRefParTypeScattConformBondAngle;
extern const RefParType *gpRefParTypeScattConformDihedAngle;
extern const RefParType *gpRefParTypeScattConformX;
extern const RefParType *gpRefParTypeScattConformY;
extern const RefParType *gpRefParTypeScattConformZ;
extern const RefParType *gpRefParTypeScattOccup;
class NiftyStaticGlobalObjectsInitializer_Scatterer
{
   public:
      NiftyStaticGlobalObjectsInitializer_Scatterer()
      {
         if (mCount++ == 0)
         {
            gpRefParTypeScatt=new RefParType(gpRefParTypeObjCryst,"Scatterer");
            gpRefParTypeScattTransl=new RefParType(gpRefParTypeScatt,"Translation");
            gpRefParTypeScattTranslX=new RefParType(gpRefParTypeScattTransl,"Translation along X");
            gpRefParTypeScattTranslY=new RefParType(gpRefParTypeScattTransl,"Translation along Y");
            gpRefParTypeScattTranslZ=new RefParType(gpRefParTypeScattTransl,"Translation along Z");
            gpRefParTypeScattOrient=new RefParType(gpRefParTypeScatt,"Orientation");
            gpRefParTypeScattConform=new RefParType(gpRefParTypeScatt,"Conformation");
            gpRefParTypeScattConformBondLength=new RefParType(gpRefParTypeScattConform,"BondLengths");
            gpRefParTypeScattConformBondAngle=new RefParType(gpRefParTypeScattConform,"Bond Angles");
            gpRefParTypeScattConformDihedAngle=new RefParType(gpRefParTypeScattConform,"Dihedral Angles ");
            gpRefParTypeScattConformX=new RefParType(gpRefParTypeScattConform,"Orth. X coordinates");
            gpRefParTypeScattConformY=new RefParType(gpRefParTypeScattConform,"Orth. Y coordinates");
            gpRefParTypeScattConformZ=new RefParType(gpRefParTypeScattConform,"Orth. Z coordinates");
            gpRefParTypeScattOccup=new RefParType(gpRefParTypeScatt,"Occupancy");
         }
      }
      ~NiftyStaticGlobalObjectsInitializer_Scatterer()
      {
         if (--mCount == 0)
         {
            delete gpRefParTypeScatt;
            delete gpRefParTypeScattTransl;
            delete gpRefParTypeScattTranslX;
            delete gpRefParTypeScattTranslY;
            delete gpRefParTypeScattTranslZ;
            delete gpRefParTypeScattOrient;
            delete gpRefParTypeScattConform;
            delete gpRefParTypeScattConformBondLength;
            delete gpRefParTypeScattConformBondAngle;
            delete gpRefParTypeScattConformDihedAngle;
            delete gpRefParTypeScattConformX;
            delete gpRefParTypeScattConformY;
            delete gpRefParTypeScattConformZ;
            delete gpRefParTypeScattOccup;
            gpRefParTypeScatt=0;
            gpRefParTypeScattTransl=0;
            gpRefParTypeScattTranslX=0;
            gpRefParTypeScattTranslY=0;
            gpRefParTypeScattTranslZ=0;
            gpRefParTypeScattOrient=0;
            gpRefParTypeScattConform=0;
            gpRefParTypeScattConformBondLength=0;
            gpRefParTypeScattConformBondAngle=0;
            gpRefParTypeScattConformDihedAngle=0;
            gpRefParTypeScattConformX=0;
            gpRefParTypeScattConformY=0;
            gpRefParTypeScattConformZ=0;
            gpRefParTypeScattOccup=0;
         }
      }
   private:
      static long mCount;
};
static NiftyStaticGlobalObjectsInitializer_Scatterer NiftyStaticGlobalObjectsInitializer_Scatterer_counter;


//######################################################################
//
//      SCATTERER
/** \brief Generic type of scatterer: can be an atom, or a more complex assembly of atoms.
*
* A Scatterer is able to give its position (in fractionnal coordinates)
* in the unit cell, and more generally the position of all point
* scattering centers (ScatteringComponent), along with the ScatteringPower
* associated with each position.
*
* For simple atoms, there is only one scattering position with the associated
* scattering power (scattering factor, anomalous, thermic). For complex
* scatterers (molecules: ZScatterer) there are as many positions as atoms.
*
* A scatterer also has a few functions to display itself in 3D
*
* This is an abstract base class.
*/
//######################################################################

class Scatterer:virtual public RefinableObj
{
   public:
      /// Constructor
      Scatterer();
      /// Copy Constructor
      Scatterer(const Scatterer &old);
      /// Destructor
      virtual ~Scatterer();
      /// \internal so-called Virtual copy constructor, needed to make copies
      /// of arrays of Scatterers
      virtual Scatterer* CreateCopy() const=0;
      virtual const string& GetClassName() const;

      /// Number of components in the scatterer (eg number of point scatterers)
      virtual int GetNbComponent() const=0;

      /** \brief Get the list of all scattering components for this scatterer.
      *
      * This is the most important function of this class, giving the list of
      * scattering positions along with the associated ScatteringPower.
      *
      */
      virtual const ScatteringComponentList& GetScatteringComponentList() const=0;
      /// Name for the i-th component of this scatterer. If the component is an Atom,
      /// Then the name is that of the atom. Else, it is the name of the scatterer
      /// plus the component number in the scatterer plus the name of the ScatteringPower.
      /// \note It would be better to return a reference, but we don't want to keep
      /// a name for all components... Weeelll, needs some more thinking... see what
      /// performance hit results (if any).
      ///
      /// \bug does not take into account dummy atoms !!
      virtual string GetComponentName(const int i) const=0;

      /// X coordinate (fractionnal) of the scatterer (for complex scatterers,
      /// this corresponds to the position of one atom of the Scatterer, ideally
      /// it should be near the center of the Scatterer.
      REAL GetX() const;
      /// Y coordinate (fractionnal) of the scatterer (for complex scatterers,
      /// this corresponds to the position of one atom of the Scatterer, ideally
      /// it should be near the center of the Scatterer.
      REAL GetY() const;
      /// Z coordinate (fractionnal) of the scatterer (for complex scatterers,
      /// this corresponds to the position of one atom of the Scatterer, ideally
      /// it should be near the center of the Scatterer.
      REAL GetZ() const;

      /** \brief Get the occupancy of the scatterer (0. -> 1.0)
      *
      * The occupancy is given in %, and must take into account the
      * 'special position' character of atoms (ie an atom on a 2fold axis
      * should have at most a .5 population, etc...).
      * For a multi-atom scatterer (polyhedra), this is the \b overall occupancy
      * of the scatterer, affecting all atoms.
      */
      REAL GetOccupancy() const ;

      /// X coordinate (fractionnal) of the scatterer (for complex scatterers,
      /// this corresponds to the position of one atom of the Scatterer, ideally
      /// it should be near the center of the Scatterer.
      virtual void SetX(const REAL x);
      /// Y coordinate (fractionnal) of the scatterer (for complex scatterers,
      /// this corresponds to the position of one atom of the Scatterer, ideally
      /// it should be near the center of the Scatterer.
      virtual void SetY(const REAL y);
      /// Z coordinate (fractionnal) of the scatterer (for complex scatterers,
      /// this corresponds to the position of one atom of the Scatterer, ideally
      /// it should be near the center of the Scatterer.
      virtual void SetZ(const REAL z);
      /** \brief Change the occupancy of the scatterer (0. -> 1.0)
      *
      * The occupancy is given in %, and must take into account the
      * 'special position' character of atoms (ie an atom on a 2fold axis
      * should have at most a .5 population, etc...).
      * For a multi-atom scatterer (polyhedra), this is the \b overall occupancy
      * of the scatterer, affecting all atoms.
      */
      virtual void SetOccupancy(const REAL occupancy) ;

      /// Conversion function.-> returns the scatt name
      ///
      /// \warning EXPERIMENTAL. DO NOT USE, as this may be removed.
      operator string()const;
      /// Print some info about the scatterer (ideally this should be one line...).
      virtual void Print() const=0;

      /** \brief Colour associated to this scatterer (using POVRay names)
      *
      */
      virtual const string& GetColour()const;
      /** \brief Colour associated to this scatterer, 3 RGB Coordinates.
      *
      */
      virtual const float* GetColourRGB()const;

      /** \brief \internal Output a description of the scatterer for POVRay.
      * This should only be called by the Crystal Object to which belongs this
      * scatterer.
      *
      */
      virtual ostream& POVRayDescription(ostream &os,
                                         const CrystalPOVRayOptions &options)const=0;

#ifdef OBJCRYST_GL
      /** \internal Create an OpenGL Display List of the scatterer. This should only
      * be called by a Crystal object.
      *
      * \param noSymmetrics: if false (the default), then all symmetrics are shown in the
      * 3D display, within the limits defined by the min/max parameters
      * \ param xMin,xMax,yMin,yMax,zMin,zMax: in fractionnal coordinates, the region
      * in which we want scaterrer to be displayed. The test is made on the center
      * of the scatterer (eg a ZScatterer (molecule) will not be 'cut' on the border).
      * \param displayNames: if true, only the names of the scatterers will be displayed,
      * at the position of the scatterers (to actually see them, they will have to
      * be translated with respect to the drawing of the scatterers).
      * \param hideHydrogens: if true, do not display hydrogens/deuterium and their bonds
      * \param fadeDistance: atoms which are beyond the display limits are still showm, but
      * with transparency which is progressively fading up to a certain distance.
      * \param fullMoleculeInLimits: if true, molecules which are centered inside the
      * display limits will be shown completely without fading. (default=false)
      */
      virtual void GLInitDisplayList(const bool noSymmetrics=false,
                                     const REAL xMin=-.1,const REAL xMax=1.1,
                                     const REAL yMin=-.1,const REAL yMax=1.1,
                                     const REAL zMin=-.1,const REAL zMax=1.1,
                                     const bool displayEnantiomer=false,
                                     const bool displayNames=false,
                                     const bool hideHydrogens=false,
                                     const REAL fadeDistance=0,
                                     const bool fullMoleculeInLimits=false)const=0;
#endif  // OBJCRYST_GL

      /// Last time anything in the scatterer was changed (atoms, positions, scattering power)
      const RefinableObjClock& GetClockScatterer()const;
      /// Last time anything in the scatterer was changed (atoms, positions, scattering power)
      RefinableObjClock& GetClockScatterer();
      /// Set the crystal in which is included this Scatterer
      void SetCrystal(Crystal&);
      /// In which crystal is this Scatterer included ?
      const Crystal& GetCrystal()const;
      /// In which crystal is this Scatterer included ?
      Crystal& GetCrystal();
   protected:
      /// \internal Prepare refinable parameters for the scatterer object
      virtual void InitRefParList()=0;
      /** Get RGB Colour coordinates from Colour Name. Note that the colour
      * used for display is usually that of the ScatteringPower, not that of the Scatterer
      *
      */
      virtual void InitRGBColour();
      /// Last time the ScatteringComponentList was generated
      const RefinableObjClock& GetClockScattCompList()const;
      ///coordinates of the scatterer (or of its center..)
      CrystVector_REAL mXYZ;

      /** \brief  Occupancy : 0 <= occ <= 1
      * For a multi-atom scatterer (polyhedron,..), this is the \b overall occupancy
      * of the scatterer (affects all components of the scatterer).
      */
      REAL mOccupancy;
      /// Colour for this scatterer (from POVRay)
      string mColourName;
      /// Colour for this scatterer using RGB
      float mColourRGB[3];
      /// Last time anything (number of atoms, positions, scattering power) was changed
      RefinableObjClock mClockScatterer;
      /// \internal Last time the ScatteringComponentList was generated
      mutable RefinableObjClock mClockScattCompList;
      /// The crystal in which the Scatterer is
      /// This is needed so that we can know which scattering powers are available
      /// in the crystal, and also to convert fractionnal to orthonormal
      /// coordinates (for some scatterers only).
      /// It cannot be const since some scatterers may want to modify
      /// the list of scattering powers of the crystal...
      Crystal *mpCryst;

   private:
   #ifdef __WX__CRYST__
   public:
      virtual WXCrystObjBasic* WXCreate(wxWindow*);
      friend class WXScatterer;
   #endif
};

/// Global registry for all Scatterer objects
extern ObjRegistry<Scatterer> gScattererRegistry;


}//namespace

#include "ObjCryst/ObjCryst/Crystal.h"

#endif //_OBJCRYST_SCATTERER_H_
