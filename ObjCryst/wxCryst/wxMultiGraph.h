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

#ifndef _WX_WXMULTIGRAPH_H_
#define _WX_WXMULTIGRAPH_H_

#include <valarray>
#include <map>
#include <string>

#include "wx/wxprec.h"
#ifdef __BORLANDC__
    #pragma hdrstop
#endif
#ifndef WX_PRECOMP
    #include "wx/wx.h"
#endif

/* WX

Select colour using wxColourDatabase (wxTheColourDatabase

All display is normalized to max-min, i.e. each tracked value has
its own scale. Zooming is therefore relative to [0.0-1.0]

Choose displayed Y values using a wxChoice, for both Y axis

The list of values to be displayed can be chosen from a wxListBox

*/
namespace ObjCryst
{
class WXMultiGraph:public wxWindow
{
   public:
      WXMultiGraph(wxFrame *frame);
      virtual ~WXMultiGraph();
      void OnPaint(wxPaintEvent &event);
      void OnMouse(wxMouseEvent &event);
      void OnMouseWheel(wxMouseEvent &event);
      void OnKeyDown(wxKeyEvent &event);
      /** Add a graph. This returns an ID identifying the graph
      *
      */
      unsigned long AddGraph(const std::string &name);
      /** Set data for a given graph.
      *
      * The two arrays \b must have the same number of elements
      */
      void SetGraphData(const unsigned long id,const std::valarray<float> &vx,
                        const std::valarray<float> &vy);
      /** Remove graph.
      *
      */
      virtual void DeleteGraph(const unsigned long id);
      /** Auto-scale graph, i.e. bring min& max along both axes to
      * the min&max of a given graph. However
      *
      */
      void AutoScale(const long id=-1,const bool xmin=true,const bool xmax=true,
                                      const bool ymin=true,const bool ymax=true);
      void OnUpdateUI(wxUpdateUIEvent &event);
      virtual void UpdateDisplay();
   private:
      /// Convert data to screen (pixel) coordinates
      void Data2Screen(float &x,float &y);
      /// Convert screen (pixel) to data coordinates
      void Screen2Data(float &x,float &y);
      struct GraphData
      {
         std::string name;
         std::valarray<float> vx;
         std::valarray<float> vy;
         float xmin,xmax,ymin,ymax;
      };
      std::map<unsigned long, GraphData> mvData;
      /// The \e current min & max values along x and y.
      float mMinX,mMaxX,mMinY,mMaxY;
      /// The margins in pixels around the graph
      long mLeft,mRight,mTop,mBottom;
      /// Pop-up menu
      wxMenu* mpPopUpMenu;
      /// Are we within a dragging event ?
      bool mIsDragging;
      /// dragging origin (in reduced coordinates)
      float mDragX0,mDragY0;
      /// Mutex for the data
      wxMutex mMutexData;
      /// parent frame
      wxFrame *mpParentFrame;
   DECLARE_EVENT_TABLE()
};


}//namespace
#endif
