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

#include <list>

#include "wxCryst/wxRefinableObj.h"

#include "wx/glcanvas.h"
#include "wx/grid.h"

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
class WXCrystal;

class WXCrystalScrolledGridWindow:public wxGrid
{
   public:
      WXCrystalScrolledGridWindow(wxWindow* parent, WXCrystal* pWXMol,long id=-1);
      virtual ~WXCrystalScrolledGridWindow();
   private:
      /// The WXCrystal window which created this window, and who should be told
      /// if it is destroyed.
      WXCrystal* mpWXCrystal;
};


/// wxCryst class for Crystals
class WXCrystal: public WXRefinableObj
{
   public:
      WXCrystal(wxWindow *parent, Crystal*);
      ~WXCrystal();
      virtual void CrystUpdate(const bool updateUI=false,const bool mutexlock=false);
      #ifdef OBJCRYST_GL
      /// Update the OpenGL Display List
      void UpdateGL(const bool onlyIndependentAtoms=false,
                    const REAL xMin=-.1,const REAL xMax=1.1,
                    const REAL yMin=-.1,const REAL yMax=1.1,
                    const REAL zMin=-.1,const REAL zMax=1.1);
      /// Gets the integer index of the OpenGL display list. Wait, if necessary, for the list
      /// not to be used any more. When finished, ReleaseCrystalGLDisplayList() must be called.
      /// if AtomName=true, then the display list returned is the one with the atom names.
      int GetCrystalGLDisplayList(const bool atomName=false)const;
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
      void OnMenuSetRelativeXYZLimits(wxCommandEvent & WXUNUSED(event));
      bool OnChangeName(const int id);
      void UpdateUI(const bool mutexlock=false);
      Crystal& GetCrystal();
      const Crystal& GetCrystal()const;
      void OnMenuShowScattPowWindow(wxCommandEvent &event);
      void OnEditGridScattPow(wxGridEvent &e);
      void OnEditGridScattPowAntiBump(wxGridEvent &e);
      void OnEditGridScattPowBondValence(wxGridEvent &e);
      void NotifyDeleteListWin(WXCrystalScrolledGridWindow *win);
      virtual bool Enable(bool enable=true);
   private:
      Crystal* mpCrystal;
      /// SpaceGroup
         WXFieldName* mpFieldSpacegroup;
      /// Scatterers
         WXRegistry<Scatterer>* mpWXScattererRegistry;
      /// Scattering Powers
         WXRegistry<ScatteringPower>* mpWXScatteringPowerRegistry;

      /// Structure to store the scattering power parameters
      struct RowScattPow
      {
         RowScattPow();
         std::string mName;
         /// Last displayed values
         REAL mBiso,mFormalCharge,mR,mG,mB,mMaximumLikelihoodError,mNbGhostAtoms;
         /// Last displayed values for antibump
         std::vector<REAL> mvAntiBumpDistance;
         /// Last displayed values for bond valence
         std::vector<REAL> mvBondValenceRo;
         /// True if we need to update the displayed values
         bool mNeedUpdateUI;
         /// Current position in the list of rows or columns
         int mIdx;
      };
      
      WXCrystalScrolledGridWindow* mpScattPowWin;
      WXCrystalScrolledGridWindow* mpAntiBumpWin;
      WXCrystalScrolledGridWindow* mpBondValenceWin;
      
      std::map<ScatteringPowerAtom*,RowScattPow> mvpRowScattPow;
      
      /// Flag to indicate that we are updating values in the wxGrid data,
      /// and that it is not the user inputing data.
      bool mIsSelfUpdating;
      
      #ifdef OBJCRYST_GL
      //OpenGl
         /// OpenGL Display of the Crystal-Display List. Updated each time CrystUpdate() is called.
         unsigned int mCrystalGLDisplayList;
         /// OpenGL Display of the Crystal-Display List. Updated each time CrystUpdate() is called.
         /// This is the list with all the scatterer (atoms) names
         unsigned int mCrystalGLNameDisplayList;
         /// the frame in which the crystal is displayed. There can only be one...
         WXGLCrystalCanvas* mpCrystalGL;
      #endif
      /// Mutex used when updating the OpenGL display List, between background and main thread
      wxMutex mMutexGLUpdate;
      /// wxCondition used when updating the OpenGL display List, between background and main thread
      wxCondition *mpConditionGLUpdate;
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
      /// Generate POVRay script to draw
      void POVRayDescription(ostream &os,const float contourValue,
                             const CrystalPOVRayOptions &options)const;
      /** Import map from a '.grd' GSAS/EXPGUI map.
      * Returns 0 on error, 1 on success
      * \param filename: the file with the fourier map
      */ 
      int ImportGRD(const string&filename);
      /** Import map with DSN6 format.
      * Returns 0 on error, 1 on success
      * \param filename: the file with the fourier map
      */ 
      int ImportDSN6(const string&filename);
      /// Name associated to this map (the filename)
      const string & GetName()const;
      /// Get the value of the map at a given set of fractionnal coordinates
      REAL GetValue(const REAL x,const REAL y,const REAL z)const;
      /// Max value of the map
      REAL Max()const;
      /// Min value of the map
      REAL Min()const;
      /// Mean value of the map
      REAL Mean()const;
      /// Standard Deviation of the map
      REAL StandardDeviation()const;
   private:
      /// The crystal corresponding to this map
      const Crystal *mpCrystal;
      /// The map data points
      CrystArray3D_REAL mPoints;
      /// Name associated to this map (the filename)
      string mName;
      /// Min and max value of the map
      REAL mMin,mMax;
      /// Mean value of the map
      REAL mMean;
      /// Standard Deviation of the map
      REAL mStandardDeviation;
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
      /// Show Wire/Filled ?
      bool ShowWire()const;
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
      void OnChangeLimits(wxCommandEvent &event);
      /// Redraw the structure (special function to ensure complete redrawing under windows...)
      void OnUpdateUI(wxUpdateUIEvent& WXUNUSED(event));
      void OnShowCrystal(wxCommandEvent & WXUNUSED(event));
      void OnShowAtomLabel(wxCommandEvent & WXUNUSED(event));
      void OnShowCursor(wxCommandEvent & WXUNUSED(event));
      void OnSetCursor(wxCommandEvent & WXUNUSED(event));
      void OnLoadFourierGRD(wxCommandEvent & WXUNUSED(event));
      void OnLoadFourierDSN6(wxCommandEvent & WXUNUSED(event));
      void AddFourier(UnitCellMapImport*);
      void OnAddContour(wxCommandEvent & WXUNUSED(event));
      void OnChangeContour(wxCommandEvent & WXUNUSED(event));
      void OnShowFourier(wxCommandEvent & WXUNUSED(event));
      void OnFourierChangeColor(wxCommandEvent & WXUNUSED(event));
      void OnUnloadFourier(wxCommandEvent & WXUNUSED(event));
      void OnShowWire(wxCommandEvent & WXUNUSED(event));
      void OnFourierChangeBbox(wxCommandEvent & WXUNUSED(event));
      /// Save view as povray file
      void OnPOVRay(wxCommandEvent & WXUNUSED(event));
      // get bounding box for atoms display
      BBox GetCellBBox();
      // get bounding box for display of Fourier map
      BBox GetMapBBox();
      virtual void SetCurrent();
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
      /// Build the 96 display lists for the font
      void BuildGLFont();
      /// Delete the 96 display lists for the font
      void DeleteGLFont()const;
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
      
      // Display lists for the font used to display strings
      mutable bool mIsGLFontBuilt;
      mutable int mGLFontDisplayListBase;
   DECLARE_EVENT_TABLE()
};
#endif


} //namespace

#endif //_VFN_WX_CRYSTAL_H_
