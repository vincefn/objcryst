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
*  source file ObjCryst++ Tracker class
*
*/
#include "RefinableObj/Tracker.h"

using namespace std;

namespace ObjCryst
{
////////////////////////////////////////////////////////////////////////
//
//    Tracker
//
////////////////////////////////////////////////////////////////////////
Tracker::Tracker(const string &name)
:mName(name)
{}

Tracker::~Tracker()
{}

const string& Tracker::GetName()const{return mName;}
void Tracker::AppendValue(const long n)
{mvValues[n]=this->ReadValue();}

void Tracker::Clear()
{
   mvValues.clear();
}

const std::map<long,REAL>& Tracker::GetValues()const{return mvValues;}
std::map<long,REAL>& Tracker::GetValues(){return mvValues;}



////////////////////////////////////////////////////////////////////////
//
//    MainTracker
//
////////////////////////////////////////////////////////////////////////

MainTracker::MainTracker()
{
   #ifdef __WX__CRYST__
   mpWXTrackerGraph=0;
   #endif
}

MainTracker::~MainTracker()
{
   #ifdef __WX__CRYST__
   this->WXDelete();
   #endif
   this->ClearTrackers();
}
void MainTracker::AddTracker(Tracker *t)
{
   mvpTracker.insert(t);
   mClockTrackerList.Click();
   this->UpdateDisplay();
}

void MainTracker::AppendValues(const long nb)
{
   for(std::set<Tracker*>::iterator pos=mvpTracker.begin(); pos!=mvpTracker.end();++pos)
      (*pos)->AppendValue(nb);
   mClockValues.Click();
}

void MainTracker::ClearTrackers()
{
   std::set<Tracker*>::iterator pos;
   for(pos=mvpTracker.begin();pos!=mvpTracker.end();++pos) delete *pos;
   mvpTracker.clear();
   mClockTrackerList.Click();
   this->UpdateDisplay();
}

void MainTracker::ClearValues()
{
   std::set<Tracker*>::iterator pos;
   for(pos=mvpTracker.begin();pos!=mvpTracker.end();++pos) (*pos)->Clear();
   mClockValues.Click();
   this->UpdateDisplay();
}

void MainTracker::SaveAll(std::ostream &os)const
{
   std::set<Tracker*>::const_iterator posT,posT0;
   os<<"#Trial ";
   for(posT=mvpTracker.begin();posT!=mvpTracker.end();++posT) os<<(*posT)->GetName()<<" ";
   os<<endl;
   
   posT0=mvpTracker.begin();
   std::map<long,REAL>::const_iterator pos0,pos;
   for(pos0=(*posT0)->GetValues().begin();pos0!=(*posT0)->GetValues().end();++pos0)
   {
      const long k=pos0->first;
      os<<k<<" ";
      for(posT=mvpTracker.begin();posT!=mvpTracker.end();posT++)
      {
         pos=(*posT)->GetValues().find(k);
         if(pos==(*posT)->GetValues().end()) os << -1.0 <<" ";
         else os << pos->second <<" ";
      }
      os<<endl;
   }
}

const std::set<Tracker*> &MainTracker::GetTrackerList()const
{
   return mvpTracker;
}

void MainTracker::UpdateDisplay()const
{
   #ifdef __WX__CRYST__
   if(0!=mpWXTrackerGraph)mpWXTrackerGraph->UpdateDisplay();
   #endif
}
const RefinableObjClock& MainTracker::GetClockTrackerList()const{return mClockTrackerList;}
const RefinableObjClock& MainTracker::GetClockValues()const{ return mClockValues;}
#ifdef __WX__CRYST__
WXTrackerGraph* MainTracker::WXCreate(wxFrame *frame)
{
   if(0==mpWXTrackerGraph) mpWXTrackerGraph=new WXTrackerGraph(frame,this);
   return mpWXTrackerGraph;
}
WXTrackerGraph* MainTracker::WXGet(){return mpWXTrackerGraph;}
void MainTracker::WXDelete()
{
   if(0!=mpWXTrackerGraph)
   {
      delete mpWXTrackerGraph;
      mpWXTrackerGraph=0;
   }
}
void MainTracker::WXNotifyDelete()
{
   mpWXTrackerGraph=0;
}
#endif


}//namespace
