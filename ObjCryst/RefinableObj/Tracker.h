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

#include <set>
#include <vector>
#include <utility>
#include "ObjCryst/RefinableObj/RefinableObj.h"

#ifdef __WX__CRYST__
namespace ObjCryst
{
class MainTracker;
class Tracker;
}
#include "ObjCryst/wxCryst/wxTrackerGraph.h"
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
      const std::set<Tracker*> &GetTrackerList()const;
      /// Update display, if any
      void UpdateDisplay()const;
      /// Get last time a tracker was added
      const RefinableObjClock& GetClockTrackerList()const;
      /// Get last time values were whanged
      const RefinableObjClock& GetClockValues()const;
   private:
      std::set<Tracker*> mvpTracker;
      /// Last time a tracker was added
      RefinableObjClock mClockTrackerList;
      /// Last time values were whanged
      RefinableObjClock mClockValues;
   #ifdef __WX__CRYST__
   public:
      WXTrackerGraph* WXCreate(wxFrame*);
      WXTrackerGraph* WXGet();
      void WXDelete();
      void WXNotifyDelete();
   protected:
      WXTrackerGraph *mpWXTrackerGraph;
   #endif
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

}//namespace
#endif
