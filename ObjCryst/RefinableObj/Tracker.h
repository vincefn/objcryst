/*  ObjCryst++ Object-Oriented Crystallographic Library
    (c) 2005- Vincent Favre-Nicolin vincefn@users.sourceforge.net

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
*  header file ObjCryst++ Tracker class
*
*/

#ifndef _REFINABLEOBJ_TRACKER_H_
#define _REFINABLEOBJ_TRACKER_H_

#include <list>
#include <vector>
#include <utility>
#include "RefinableObj/RefinableObj.h"
#include "wx/wxprec.h"

#ifdef __BORLANDC__
    #pragma hdrstop
#endif
#ifndef WX_PRECOMP
    #include "wx/wx.h"
#endif

namespace ObjCryst
{

/** A class to track the variation of parameters as a function
* of a number of cycles/trials.
*
* This is an abstract base class.
*/
class Tracker
{
   public:
      Tracker(const std::string &name);
      virtual ~Tracker();
      const std::string& GetName()const;
      void AppendValue(const long trial);
      /// Removes all stored values
      void Clear();
      const std::map<long,REAL>& GetValues() const;
      std::map<long,REAL>& GetValues();
   protected:
      virtual REAL ReadValue()=0;
      std::map<long,REAL> mvValues;
      std::string mName;
};

/** A class to hold all trackers
*
* All trackers need not have the same 
*/
class MainTracker
{
   public:
      MainTracker();
      ~MainTracker();
      void AddTracker(Tracker *t);
      void AppendValues(const long trial);
      /// Removes all Trackers
      void ClearTrackers();
      /// Removes all stored values
      void ClearValues();
      /// Will save to a single file if all recorded trial numbers are the same
      /// Otherwise ?
      void SaveAll(std::ostream &out)const;
   private:
      std::list<Tracker*> mvpTracker;
};

/** Tracker for objects (RefinableObj, Crystal, PowderPattern, RefPar,...) 
*
*/
template <class T> class TrackerObject:public Tracker
{
   public:
      TrackerObject(const std::string &name, const T&obj, REAL (T::*f)() const):
         Tracker(name),mpObj(&obj),mfp(f)
         {}
   private:
      const T *mpObj;
      REAL (T::*mfp)() const;
      REAL ReadValue(){return (mpObj->*mfp)();}
};

/* WX

Select colour using wxColourDatabase (wxTheColourDatabase

All display is normalized to max-min, i.e. each tracked value has
its own scale. Zooming is therefore relative to [0.0-1.0]

Choose displayed Y values using a wxChoice, for both Y axis

The list of values to be displayed can be chosen from a wxListBox

*/

class WXMainTracker
{
};

template<class T,class U> class WXMultiGraph:public wxWindow
{
   public:
      void OnPaint(wxPaintEvent &event);
      void OnMouse(wxMouseEvent &event);
      void OnMouseWheel(wxMouseEvent &event);
      void SetData(const std::vector<const string> &vName,
                   const std::vector<std::list<pair<U,T> > > &data);
   private:
      /// Convert (screen) pixel coordinates to reduced coordinates between 0 and 1
      void Screen2Reduced(float &x,float &y);
      /// Convert (screen) pixel coordinates to reduced coordinates between 0 and 1
      void Reduced2Screen(float &x,float &y);
      /// Convert data #i to reduced coordinates (i.e. coordinates between 0 and 1.0)
      void Data2Reduced(unsigned long i,float &x,float &y);
      /// Convert reduced to data #i coordinates
      void Reduced2Data(unsigned long i,float &x,float &y);
      /// Convert data #i to screen (pixel) coordinates
      void Data2Screen(unsigned long i,float &x,float &y);
      /// Convert screen (pixel) to data #i coordinates
      void Screen2Data(unsigned long i,float &x,float &y);
      std::vector<const std::string*> mvName;
      std::vector<std::list<pair<T,U> > > mvData;
      mutable std::vector<T> mvMinX,mvMaxX;
      mutable std::vector<U> mvMinY,mvMaxY;
      /// The \e current min & max values along x and y, for reduced coordinates.
      /// This means that graph #i will have xmin=mvMinX[i]+mMinX*(mvMaxX[i]-mvMinX[i]),
      /// etc...
      float mMinX,mMaxX,mMinY,mMaxY;
      /// Pop-up menu
      wxMenu* mpPopUpMenu;
      /// Are we within a dragging event ?
      bool mIsDragging;
      /// dragging origin (in reduced coordinates)
      float mDragX0,mDragY0;
      DECLARE_EVENT_TABLE()
};

}//namespace
#endif
