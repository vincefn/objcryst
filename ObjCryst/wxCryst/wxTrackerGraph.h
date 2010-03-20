/*  ObjCryst++ Object-Oriented Crystallographic Library
    (c) 2005-     Vincent Favre-Nicolin vincefn@users.sourceforge.net

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

#ifndef _WX_WXTRACKER_H_
#define _WX_WXTRACKER_H_

#include <map>

#include "wx/wxprec.h"

#ifdef __BORLANDC__
    #pragma hdrstop
#endif
#ifndef WX_PRECOMP
    #include "wx/wx.h"
#endif

namespace ObjCryst
{
class WXTrackerGraph;
} //namespace

#include "ObjCryst/wxCryst/wxMultiGraph.h"
#include "ObjCryst/RefinableObj/Tracker.h"

namespace ObjCryst
{

class WXTrackerGraph:public WXMultiGraph
{
   public:
      WXTrackerGraph(wxFrame *frame, MainTracker *tracker);
      virtual ~WXTrackerGraph();
      /// reads new values from the MainTracker, and asks for a repaint.
      virtual void UpdateDisplay();
      virtual void DeleteGraph(const unsigned long id);
   private:
      MainTracker *mpMainTracker;
      std::map<Tracker*,long> mvId;
      /// Last time a tracker was added to the graph
      RefinableObjClock mClockGraphList;
      /// Last time values were added to the graph
      RefinableObjClock mClockGraphValues;
};


}//namespace
#endif
