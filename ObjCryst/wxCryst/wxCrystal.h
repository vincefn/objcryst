/*
* LibCryst++ : a Crystallographic computing library in C++
*
*  (c) 2000 Vincent FAVRE-NICOLIN
*           Laboratoire de Cristallographie
*           24, quai Ernest-Ansermet, CH-1211 Geneva 4, Switzerland
*  Contact: Vincent.Favre-Nicolin@cryst.unige.ch
*           Radovan.Cerny@cryst.unige.ch
*
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
      /// Update the OpenGL Display List
      void UpdateGL();
      /// Gets the integer index of the OpenGL display list. Wait, if necessary, for the list
      /// not to be used any more. When finished, ReleaseCrystalGLDisplayList() must be called.
      int GrabCrystalGLDisplayList()const;
      void ReleaseCrystalGLDisplayList()const;
      /// Create OpenGL Display of the Crystal Structure
      void OnMenuCrystalGL(wxCommandEvent & WXUNUSED(event));
      /// Tell this object that its 3D OpenGL display has been destroyed
      void NotifyCrystalGLDelete();
      void OnMenuSaveCIF(wxCommandEvent & WXUNUSED(event));
      void OnMenuSaveText(wxCommandEvent & WXUNUSED(event));
      void OnMenuAddScattPowAtom(wxCommandEvent & WXUNUSED(event));
      void OnMenuRemoveScattPow(wxCommandEvent & WXUNUSED(event));
      void OnMenuAddScatterer(wxCommandEvent & event);
      void OnMenuRemoveScatterer(wxCommandEvent & WXUNUSED(event));
      void OnMenuAddAntiBumpDist(wxCommandEvent & WXUNUSED(event));
      bool OnChangeName(const int id);
   private:
      Crystal* mpCrystal;
      /// Lattice
         WXFieldRefPar* mpFieldLatticeA;
         WXFieldRefPar* mpFieldLatticeB;
         WXFieldRefPar* mpFieldLatticeC;
         WXFieldRefPar* mpFieldLatticeAlpha;
         WXFieldRefPar* mpFieldLatticeBeta;
         WXFieldRefPar* mpFieldLatticeGamma;
      /// SpaceGroup
         WXFieldName* mpFieldSpacegroup;
      /// Scatterers
         WXRegistry<Scatterer>* mpWXScattererRegistry;
      /// Scattering Powers
         WXRegistry<ScatteringPower>* mpWXScatteringPowerRegistry;
         
      //OpenGl
         /// OpenGL Display of the Crystal-Display List. Updated each time CrystUpdate() is called.
         unsigned int mCrystalGLDisplayList;
         /// This is true when the display list is being used
         mutable bool mCrystalGLDisplayListIsLocked;
         /// the frame in which the crystal is displayed. There can only be one...
         WXGLCrystalCanvas* mpCrystalGL;
   DECLARE_EVENT_TABLE()
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
      void CrystUpdate();
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
   DECLARE_EVENT_TABLE()
};


} //namespace

#endif //_VFN_WX_CRYSTAL_H_
