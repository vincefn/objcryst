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
/*
*  header file for the RefinablePar and RefinableObj classes
*
* This is still in early development stages !! Not secure !
*
*/

#ifndef _VFN_WX_CRYSTAL_H_
#define _VFN_WX_CRYSTAL_H_

#include "wxCryst/wxRefinableObj.h"

#include "wx/glcanvas.h"

#include "ObjCryst/Crystal.h"

#include "wxCryst/MC.h"

namespace ObjCryst
{
typedef struct
  {
    float       xMin;
    float       xMax;
    float       yMin;
    float       yMax;
    float       zMin;
    float       zMax;
  }
  BBox;

typedef struct
  {
    float       x;
    float       y;
    float       z;
  }
  Triple;

class WXGLCrystalCanvas;

/// wxCryst class for Crystals
class WXCrystal: public WXRefinableObj
{
   public:
      WXCrystal(wxWindow *parent, Crystal*);
      virtual void CrystUpdate();
      #ifdef OBJCRYST_GL
      /// Update the OpenGL Display List
      void UpdateGL(const bool onlyIndependentAtoms=false,
                    const REAL xMin=-.1,const REAL xMax=1.1,
                    const REAL yMin=-.1,const REAL yMax=1.1,
                    const REAL zMin=-.1,const REAL zMax=1.1);
      /// Gets the integer index of the OpenGL display list. Wait, if necessary, for the list
      /// not to be used any more. When finished, ReleaseCrystalGLDisplayList() must be called.
      /// if AtomName=true, then the display list returned is the one with the atom names.
      int GrabCrystalGLDisplayList(const bool atomName=false)const;
      void ReleaseCrystalGLDisplayList()const;
      bool GLDisplayListIsLocked()const;
      /// Create OpenGL Display of the Crystal Structure
      void OnMenuCrystalGL(wxCommandEvent & WXUNUSED(event));
      /// Tell this object that its 3D OpenGL display has been destroyed
      void NotifyCrystalGLDelete();
      /// get a pointer to the 3D OpenGL display object
      WXGLCrystalCanvas * GetCrystalGL();
      #endif
      void OnMenuSaveCIF(wxCommandEvent & WXUNUSED(event));
      void OnMenuSaveText(wxCommandEvent & WXUNUSED(event));
      void OnMenuAddScattPowAtom(wxCommandEvent & WXUNUSED(event));
      void OnMenuAddScattPowSphere(wxCommandEvent & WXUNUSED(event));
      void OnMenuRemoveScattPow(wxCommandEvent & WXUNUSED(event));
      void OnMenuAddScatterer(wxCommandEvent & event);
      void OnMenuRemoveScatterer(wxCommandEvent & WXUNUSED(event));
      void OnMenuDuplicateScatterer(wxCommandEvent & WXUNUSED(event));
      void OnMenuImportMoleculeFromFenskeHallZMatrix(wxCommandEvent &event);
      void OnMenuAddAntiBumpDist(wxCommandEvent & WXUNUSED(event));
      void OnMenuSetRelativeXYZLimits(wxCommandEvent & WXUNUSED(event));
      bool OnChangeName(const int id);
      void UpdateUI();
      Crystal& GetCrystal();
      const Crystal& GetCrystal()const;
   private:
      Crystal* mpCrystal;
      /// SpaceGroup
         WXFieldName* mpFieldSpacegroup;
      /// Scatterers
         WXRegistry<Scatterer>* mpWXScattererRegistry;
      /// Scattering Powers
         WXRegistry<ScatteringPower>* mpWXScatteringPowerRegistry;
      /// Storage for the Anti-bump/merge cost
         REAL mBumpmergeCost;

      #ifdef OBJCRYST_GL
      //OpenGl
         /// OpenGL Display of the Crystal-Display List. Updated each time CrystUpdate() is called.
         unsigned int mCrystalGLDisplayList;
         /// OpenGL Display of the Crystal-Display List. Updated each time CrystUpdate() is called.
         /// This is the list with all the scatterer (atoms) names
         unsigned int mCrystalGLNameDisplayList;
         /// This is true when the display list is being used
         mutable bool mCrystalGLDisplayListIsLocked;
         /// the frame in which the crystal is displayed. There can only be one...
         WXGLCrystalCanvas* mpCrystalGL;
      #endif
   DECLARE_EVENT_TABLE()
};

#ifdef OBJCRYST_GL
/// Class to store a Fourier map, imported from another package.
/// This can also generate a 3d view (OpenGL display list) of the map.
class UnitCellMapImport
{
   public:
      /** Creator
      *
      * \param crystal: the crystal correponding to this map 
      */
      UnitCellMapImport(const Crystal&crystal);
      ~UnitCellMapImport();
      /// Perform the OpenGL drawing, to be stored in an OpenGL Display List.
      /// Assumes the color and type of drawing (GL_LINE or GL_FILL) 
      /// is chosen before calling this display list.
      void GLInitDisplayList(const float contourValue,
			     WXGLCrystalCanvas * parentCrystal) const;
      /** Import map from a '.grd' GSAS/EXPGUI map.
      * Returns 0 on error, 1 on success
      * \param filename: the file with the fourier map
      */ 
      int ImportGRD(const string&filename);
      /// Name associated to this map (the filename)
      const string & GetName()const;
      /// Get the value of the map at a given set of fractionnal coordinates
      REAL GetValue(const REAL x,const REAL y,const REAL z)const;
   private:
      /// The crystal corresponding to this map
      const Crystal *mpCrystal;
      /// The map data points
      CrystArray3D_REAL mPoints;
      /// Name associated to this map (the filename)
      string mName;
};

/// Class to store and execute OpenGL Display Lists of fourier maps.
/// Only display information is kept here.
class UnitCellMapGLList
{
   public:
      UnitCellMapGLList(const bool showWire=true,
                        const float r=1.0,const float g=0.0,const float b=0.0,
                        const float t=0.5);
      ~UnitCellMapGLList();
      /// Generates, or changes the display list.
      void GenList(const UnitCellMapImport &ucmap,
		   WXGLCrystalCanvas * parent,
                   const float contourValue=1.);
      /// Change name for this map
      void SetName(const string &name);
      /// Name for this map
      const string &GetName()const;
      /// Change the color.
      void SetColour(const float r=1.0,const float g=0.0,const float b=0.0,
                    const float t=1.0);
      /// Get the color, as a float[4] array
      const float* GetColour()const;
      /// Toggle show Wire/Filled
      void ToggleShowWire();
      /// Perform the OpenGL drawing
      void Draw()const;
   private:
      /// The index of the OpenGL display list
      unsigned int mGLDisplayList;
      /// The color to display the map
      float mColour[4];
      /// Show as wireframe (if true) or fill (if false).
      bool mShowWire;
      /// The name associated with this map display
      string mName;
};
/// Class for 3D OpenGL display of Crystal structures
class WXGLCrystalCanvas : public wxGLCanvas
{
   public:
      WXGLCrystalCanvas(WXCrystal *wxcryst,
                        wxFrame *parent, wxWindowID id=-1,
                        const wxPoint &pos=wxDefaultPosition,
                        const wxSize &size=wxDefaultSize);
      ~WXGLCrystalCanvas();
      void OnExit(wxCommandEvent &event);
      void OnPaint(wxPaintEvent &event);
      void OnSize(wxSizeEvent& event);
      void OnEraseBackground(wxEraseEvent& event);
      void OnKeyDown(wxKeyEvent& event);
      void OnKeyUp(wxKeyEvent& event);
      void OnEnterWindow( wxMouseEvent& event );
      void OnMouse( wxMouseEvent& event );
      /// This forces a new Display List (user-asked)
      void OnUpdate(wxCommandEvent & WXUNUSED(event));
      /// This is called by the Crystal in WXCrystal::UpdateGL().
      void CrystUpdate();
      void OnChangeLimits(wxCommandEvent & WXUNUSED(event));
      /// Redraw the structure (special function to ensure complete redrawing under windows...)
      void OnUpdateUI(wxUpdateUIEvent& WXUNUSED(event));
      void OnShowCrystal();
      void OnShowAtomLabel();
      void OnShowCursor();
      void OnSetCursor();
      void OnLoadFourier();
      void LoadFourier(const string&filename);
      void OnAddContour();
      void OnChangeContour();
      void OnShowFourier();
      void OnFourierChangeColor();
      void OnUnloadFourier();
      void OnShowWire();
      void OnFourierChangeBbox();
      // get bounding box for atoms display
      BBox GetCellBBox();
      // get bounding box for display of Fourier map
      BBox GetMapBBox();
   private:
      void InitGL();
      /// Shows a dialog to choose a displayed fourier map from one of those
      /// available.
      int UserSelectUnitCellMapGLList()const;
      /// Shows a dialog to choose a fourier map from one of those
      /// available.
      int UserSelectUnitCellMapImport()const;
      /** Transform (x,y) window coordinates to (x,y,z) coordinates in the
      * Crystal's orthonormal frame (with the origin at the center of the Unit Cell).
      *
      * \param x,y: on input, these contain the mouse screen coordinates,
      * and in return these are the x and y coordinates in the Crystal orthonormal
      * frame.
      * \param z: ignored on input (the z coordinate taken is that of the 
      * center of rotation).
      *
      * \note: this should be const, but we use SetCurrent() and build_rotmatrix( ,mQuat)...
      * which are not const-correct...
      */
      void UnProject(REAL &x, REAL &y, REAL &z);
      /// The parent wxFrame
      wxFrame* mpParentFrame;
      /// The owner WXCrystal
      WXCrystal* mpWXCrystal;
      /// \internal
      bool mIsGLInit;
      /// quaternion for the orientation of the display
      float mQuat [4];
      /// \internal
      float mTrackBallLastX,mTrackBallLastY;
      /// Distance from viewer to crystal (Z)
      float mDist;
      float mX0, mY0,mZ0;
      /// View Angle, in degrees
      float mViewAngle;
      /// Pop-up menu
      wxMenu* mpPopUpMenu;
      
      /// To display Fourier map
      bool mShowFourier, mShowCrystal, mShowAtomName, mShowCursor;
      // bounding box for atoms to be included
      BBox mcellbbox;
      // bounding box for display of Fourier map
      BBox mmapbbox;
      // position to use as center
      Triple mViewCntr;

      /** Imported fourier maps
      *
      * :TODO: use an auto_ptr<UnitCellMapGLList>
      */
      vector<UnitCellMapImport* > mvpUnitCellMapImport;
      /** OpenGL display lists corresponding to Fourier maps.
      *
      * We store 3 things in this list, a const pointer to the imported map,
      * the contour level used, and a pointer to the display list used.
      *
      * :TODO: use an auto_ptr<UnitCellMapGLList>
      */
      vector<pair<pair<const UnitCellMapImport*,float>,UnitCellMapGLList* > > mvpUnitCellMapGLList;
      
   DECLARE_EVENT_TABLE()
};
#endif


} //namespace

#endif //_VFN_WX_CRYSTAL_H_
