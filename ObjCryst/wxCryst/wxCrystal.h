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

#include "wx/glcanvas.h"

#include "wxCryst/wxRefinableObj.h"
#include "ObjCryst/Crystal.h"

namespace ObjCryst
{
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
      int GrabCrystalGLDisplayList()const;
      void ReleaseCrystalGLDisplayList()const;
      bool GLDisplayListIsLocked()const;
      /// Create OpenGL Display of the Crystal Structure
      void OnMenuCrystalGL(wxCommandEvent & WXUNUSED(event));
      /// Tell this object that its 3D OpenGL display has been destroyed
      void NotifyCrystalGLDelete();
      #endif
      void OnMenuSaveCIF(wxCommandEvent & WXUNUSED(event));
      void OnMenuSaveText(wxCommandEvent & WXUNUSED(event));
      void OnMenuAddScattPowAtom(wxCommandEvent & WXUNUSED(event));
      void OnMenuRemoveScattPow(wxCommandEvent & WXUNUSED(event));
      void OnMenuAddScatterer(wxCommandEvent & event);
      void OnMenuRemoveScatterer(wxCommandEvent & WXUNUSED(event));
      void OnMenuDuplicateScatterer(wxCommandEvent & WXUNUSED(event));
      void OnMenuAddAntiBumpDist(wxCommandEvent & WXUNUSED(event));
      bool OnChangeName(const int id);
      void UpdateUI();
   private:
      Crystal* mpCrystal;
      /// SpaceGroup
         WXFieldName* mpFieldSpacegroup;
      /// Scatterers
         WXRegistry<Scatterer>* mpWXScattererRegistry;
      /// Scattering Powers
         WXRegistry<ScatteringPower>* mpWXScatteringPowerRegistry;

      #ifdef OBJCRYST_GL
      //OpenGl
         /// OpenGL Display of the Crystal-Display List. Updated each time CrystUpdate() is called.
         unsigned int mCrystalGLDisplayList;
         /// This is true when the display list is being used
         mutable bool mCrystalGLDisplayListIsLocked;
         /// the frame in which the crystal is displayed. There can only be one...
         WXGLCrystalCanvas* mpCrystalGL;
      #endif
   DECLARE_EVENT_TABLE()
};

#ifdef OBJCRYST_GL
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
   private:
      void InitGL();
      /// The owner WXCrystal
      WXCrystal* mpWXCrystal;
      /// \internal
      bool mIsGLInit;
      /// quaternion for the orientation of the display
      float mQuat [4];
      /// \internal
      float mTrackBallLastX,mTrackBallLastY;
      /// Distance from viewer to crystal
      float mDist;
      /// View Angle, in degrees
      float mViewAngle;
      /// Pop-up menu
      wxMenu* mpPopUpMenu;
      float mXmin,mXmax,mYmin,mYmax,mZmin,mZmax;
   DECLARE_EVENT_TABLE()
};
#endif


} //namespace

#endif //_VFN_WX_CRYSTAL_H_
