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
#include "wxCryst/wxTrackerGraph.h"

namespace ObjCryst
{

WXTrackerGraph::WXTrackerGraph(wxFrame *frame, MainTracker *tracker):
WXMultiGraph(frame),mpMainTracker(tracker)
{
   this->UpdateDisplay();
}

WXTrackerGraph::~WXTrackerGraph()
{
   if(mpMainTracker!=0) mpMainTracker->WXNotifyDelete();
}

void WXTrackerGraph::UpdateDisplay()
{
   VFN_DEBUG_ENTRY("WXTrackerGraph::UpdateDisplay()",4)
   const std::set<Tracker*> *pList=&(mpMainTracker->GetTrackerList());
   std::set<Tracker*>::const_iterator pos;
   
   // Remove orphan graphs
   if(mClockGraphList<mpMainTracker->GetClockTrackerList())
   {
      std::map<Tracker*,long>::iterator pos1;
      for(pos1=mvId.begin();pos1!=mvId.end();++pos1)
         if(pList->find(pos1->first)==pList->end())
         {
            this->DeleteGraph(pos1->second);
            mClockGraphList.Click();
         }
   }
   unsigned long nbxMax=0;
   // Add new graphs and update data for all
   if(  (mClockGraphList  <mpMainTracker->GetClockTrackerList())
      ||(mClockGraphValues<mpMainTracker->GetClockValues()))
   {
      for(pos=pList->begin();pos!=pList->end();++pos)
      {
         if(mvId.find(*pos)==mvId.end())
         {
            mvId[*pos]=this->AddGraph((*pos)->GetName());
            mClockGraphList.Click();
         }
         const unsigned long nb=(*pos)->GetValues().size();
         std::map<long,REAL>::const_iterator pos2;
         valarray<float> vx(nb),vy(nb);
         unsigned long i=0;
         for(pos2=(*pos)->GetValues().begin();pos2!=(*pos)->GetValues().end();++pos2)
         {
            vx[i]=pos2->first;
            vy[i]=pos2->second;
            i++;
         }
         this->SetGraphData(mvId[*pos],vx,vy);
         mClockGraphValues.Click();
         if(nbxMax<vx.size()) nbxMax=vx.size();
      }
      if(nbxMax==1) this->AutoScale(-1);// first drawing, so rescale everything
      else this->AutoScale(-1,false,true,false,false);//just rescale xmax
   }
   this->WXMultiGraph::UpdateDisplay();
   VFN_DEBUG_EXIT("WXTrackerGraph::UpdateDisplay()",4)
}

void WXTrackerGraph::DeleteGraph(const unsigned long id)
{
   // remove tracker. Maybe should keep a reverse map.
   std::map<Tracker*,long>::iterator pos;
   for(pos=mvId.begin();pos!=mvId.end();++pos)
   {
      if((long)id==pos->second)
      {
         mvId.erase(pos);
         break;
      }
   }
   this->WXMultiGraph::DeleteGraph(id);
}

   
}//namespace
