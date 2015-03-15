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

#include <iostream>

// wx headers, with or without precompilation
#include "wx/wxprec.h"
#ifdef __BORLANDC__
    #pragma hdrstop
#endif
#ifndef WX_PRECOMP
    #include "wx/wx.h"
#endif
#include "wx/dcbuffer.h"
#include "wx/gdicmn.h"
#include "ObjCryst/ObjCryst/General.h"
#include "ObjCryst/wxCryst/wxMultiGraph.h"
#include "ObjCryst/wxCryst/wxCryst.h"

using namespace std;

namespace ObjCryst
{
static const char* swxColourNameList[]={
"BLACK",
"BLUE",
"RED",
"GREEN",
"BROWN",
"CYAN",
"ORANGE",
"YELLOW",
"MAGENTA",
"MAROON",
"GOLD",
"GREY",
"PINK",
"SALMON",
"PURPLE",
"DARK GREY",
"CORAL",
"AQUAMARINE",
"VIOLET",

"BLUE VIOLET",
"CADET BLUE",
"CORNFLOWER BLUE",
"DARK GREEN",
"DARK OLIVE GREEN",
"DARK ORCHID",
"DARK SLATE BLUE",
"DARK SLATE GREY DARK TURQUOISE",
"DIM GREY",
"FIREBRICK",
"FOREST GREEN",
"GOLDENROD",
"GREEN YELLOW",
"INDIAN RED",
"KHAKI",
"LIGHT BLUE",
"LIGHT GREY",
"LIGHT STEEL BLUE",
"LIME GREEN",
"MEDIUM AQUAMARINE",
"MEDIUM BLUE",
"MEDIUM FOREST GREEN",
"MEDIUM GOLDENROD",
"MEDIUM ORCHID",
"MEDIUM SEA GREEN",
"MEDIUM SLATE BLUE",
"MEDIUM SPRING GREEN",
"MEDIUM TURQUOISE",
"MEDIUM VIOLET RED",
"MIDNIGHT BLUE",
"NAVY",
"ORANGE RED",
"ORCHID",
"PALE GREEN",
"PLUM",
"SEA GREEN",
"SIENNA",
"SKY BLUE",
"SLATE BLUE",
"SPRING GREEN",
"STEEL BLUE",
"TAN",
"THISTLE",
"TURQUOISE",
"VIOLET RED",
"WHEAT",
//"WHITE",
"YELLOW GREEN."};

////////////////////////////////////////////////////////////////////////
//
//    WXMultiGraph
//
////////////////////////////////////////////////////////////////////////
static const long ID_UPDATEUI=              WXCRYST_ID();
static const long ID_MENU_AUTOSCALE=        WXCRYST_ID();

BEGIN_EVENT_TABLE(WXMultiGraph, wxWindow)
   EVT_PAINT(                                   WXMultiGraph::OnPaint)
   EVT_MOUSE_EVENTS(                            WXMultiGraph::OnMouse)
   EVT_CHAR(                                    WXMultiGraph::OnKeyDown)
   EVT_MOUSEWHEEL(                              WXMultiGraph::OnMouseWheel)
   EVT_UPDATE_UI(ID_UPDATEUI,                   WXMultiGraph::OnUpdateUI)
   EVT_SIZE(                                    WXMultiGraph::OnSize)
END_EVENT_TABLE()


WXMultiGraph::WXMultiGraph(wxFrame *frame):
wxWindow(frame,-1,wxPoint(-1,-1),wxSize(-1,-1)),
mMinX(0.0),mMaxX(1.0),mMinY(0.0),mMaxY(1.0),
mLeft(40),mRight(10),mTop(10),mBottom(25),
mIsDragging(false),mpParentFrame(frame)
{
   #ifdef VFN_CRYST_MUTEX
   cout <<"new CrystMutex("<<&mMutexData<<")for WXMultiGraph:"<<this<<endl;
   #endif
   mpPopUpMenu=new wxMenu(_T("Graph"));
   mpPopUpMenu->Append(ID_MENU_AUTOSCALE, _T("&AutoScale"));
   //mpPopUpMenu->Append(ID_POWDERGRAPH_MENU_TOGGLELABEL, "&Hide Labels");
}

WXMultiGraph::~WXMultiGraph()
{
   VFN_DEBUG_MESSAGE("WXMultiGraph::~WXMultiGraph():",4)
   delete mpPopUpMenu;
   #ifdef VFN_CRYST_MUTEX
   cout <<"Deleting CrystMutex("<<&mMutexData<<")for WXMultiGraph:"<<this<<endl;
   #endif
}

unsigned long WXMultiGraph::AddGraph(const string &name)
{
   wxMutexLocker mlock(mMutexData);
   long id;
   for(id=mvData.size();id>=0;id--)
      if(mvData.end()==mvData.find(id)) break;
   mvData[id].name=name;
   mvData[id].xmin=0.0;
   mvData[id].xmax=-1.0;
   mvData[id].ymin=0.0;
   mvData[id].ymax=-1.0;
   return id;
}

void WXMultiGraph::SetGraphData
   (const unsigned long id,const std::valarray<float> &vx,
    const std::valarray<float> &vy)
{
   bool rescale=false;
   mMutexData.Lock();
   if(mvData[id].xmin>mvData[id].xmax)rescale=true;// A new graph has been added
   mvData[id].vx.resize(vx.size());
   mvData[id].vy.resize(vx.size());
   mvData[id].vx=vx;
   mvData[id].vy=vy;
   if(vx.size()!=0)
   {
      mvData[id].xmin=mvData[id].vx.min();
      mvData[id].xmax=mvData[id].vx.max();
      mvData[id].ymin=mvData[id].vy.min();
      mvData[id].ymax=mvData[id].vy.max();
      VFN_DEBUG_MESSAGE("WXMultiGraph::SetGraphData():"<<mvData[id].vx.size()<<":"
                        << mvData[id].xmin<<","<< mvData[id].xmax<<","
                        << mvData[id].ymin<<","<< mvData[id].ymax<<",",3)
      mMutexData.Unlock();
      if(rescale)
         this->AutoScale(-1);
   }
   else mMutexData.Unlock();
}


void WXMultiGraph::DeleteGraph(const unsigned long id)
{
   mvData.erase(id);
}

void WXMultiGraph::OnPaint(wxPaintEvent &event)
{
   wxMutexLocker mlock(mMutexData);
   if(mvData.size()<1) return;
   VFN_DEBUG_ENTRY("WXMultiGraph::OnPaint()",4)
   wxBufferedPaintDC dc(this);
   this->PrepareDC(dc);
   mpParentFrame->PrepareDC(dc);
   
   dc.DestroyClippingRegion();
   dc.SetBackground(wxBrush(_T("white"), wxSOLID));
   dc.Clear();

   wxString fontInfo;
   #ifdef __WIN32__
   dc.SetFont(*wxNORMAL_FONT);
   #else
   dc.SetFont(*wxSMALL_FONT);
   #endif
   // Get Window Size
   wxCoord width,height;
   this->GetSize(&width, &height);
   
   // Draw Axis
   VFN_DEBUG_MESSAGE("WXMultiGraph::OnPaint():Axis",3)
   {
      dc.SetPen(*wxBLACK_PEN);
      dc.SetTextForeground(*wxBLACK);
      int nbTick=5;//approx.
      wxCoord tmpW,tmpH;
      float xs,ys;
      // X & Y margins.
         float yStep,xStep,dx,dy;
         dx=mMaxX-mMinX;
         dy=mMaxY-mMinY;
         if(dx<1e-6)dx=1e-6;
         if(dy<1e-6)dy=1e-6;
         yStep=pow((float)10,(float)floor(log10(dy/nbTick)));
         yStep *= floor((dy/yStep+0.1)/nbTick);
         
         xStep=pow((float)10,(float)floor(log10(dx/nbTick)));
         xStep *= floor((dx/xStep+0.1)/nbTick);

         mLeft=0;
         for(float y=yStep*ceil(mMinY/yStep);y<mMaxY;y+=yStep)
         {//get left margin from tick labels
            fontInfo.Printf(_T("%g"),y);
            dc.GetTextExtent(fontInfo, &tmpW, &tmpH);
            if((tmpW+3)>mLeft) mLeft=tmpW+3;
         }
         
         fontInfo.Printf(_T("%g"),xStep);
         dc.GetTextExtent(fontInfo, &tmpW, &tmpH);
         mBottom=tmpH*3/2+3;
         
      //Y axis
         dc.DrawLine(mLeft,height-mBottom,mLeft,mTop);
         VFN_DEBUG_MESSAGE("WXMultiGraph::OnPaint():AxisStep="<<yStep<<","<<mMinY<<","<<mMaxY,3)
         nbTick=int(dy/yStep);
         for(int i=0;i<nbTick;i++)
         {
            float y=yStep*ceil(mMinY/yStep)+i*yStep;
            xs=mMinX;
            ys=y;
            this->Data2Screen(xs,ys);
            VFN_DEBUG_MESSAGE("WXMultiGraph::OnPaint():Axis:"<<xs<<","<<ys,3)
            dc.DrawLine(wxCoord(xs-3),wxCoord(ys),wxCoord(xs+3),wxCoord(ys));
            fontInfo.Printf(_T("%g"),y);
            dc.GetTextExtent(fontInfo, &tmpW, &tmpH);
            dc.DrawText(fontInfo,wxCoord(xs-tmpW-3),wxCoord(ys-tmpH/2));
         }
      //X axis
         dc.DrawLine(mLeft,height-mBottom,width-mRight,height-mBottom);
         nbTick=int(dx/xStep);
         for(int i=0;i<nbTick;i++)
         {
            float x=xStep*ceil(mMinX/xStep)+i*xStep;
            xs=x;
            ys=mMinY;
            this->Data2Screen(xs,ys);
            dc.DrawLine(wxCoord(xs),wxCoord(ys-3),wxCoord(xs),wxCoord(ys+3));
            fontInfo.Printf(_T("%g"),x);
            dc.GetTextExtent(fontInfo, &tmpW, &tmpH);
            dc.DrawText(fontInfo,wxCoord(xs-tmpW/2),wxCoord(ys+tmpH/2));
         }
      // Axis labels;
      dc.GetTextExtent(mYLabel, &tmpW, &tmpH);
      dc.DrawText(mYLabel,wxCoord(0),wxCoord(0));
      xs=xStep*ceil(mMinX/xStep)+(nbTick-1.5)*xStep;
      ys=mMinY;
      this->Data2Screen(xs,ys);
      dc.GetTextExtent(mXLabel, &tmpW, &tmpH);
      dc.DrawText(mXLabel,wxCoord(xs-tmpW/2),wxCoord(ys+tmpH/2+3));
   }
   // Draw data
   map<unsigned long, GraphData >::const_iterator pos;
   long ix=-1,ixdrawn=-1;
   for(pos=mvData.begin();pos!=mvData.end();++pos)
   {
      ix++;
      VFN_DEBUG_MESSAGE("WXMultiGraph::OnPaint():Data#"<<ix,3)
      if((pos->second.vx.size()<1)||(pos->second.vx.size()!=pos->second.vy.size())) continue;
      ixdrawn++;
      dc.SetPen(wxPen(wxTheColourDatabase->Find(wxString::FromAscii(swxColourNameList[ix])),1,wxSOLID));
      float x1,y1,x2,y2;
      x2=pos->second.vx[0];
      y2=pos->second.vy[0];
      this->Data2Screen(x2,y2);
      for(unsigned long i=0;i<pos->second.vx.size();i++)
      {
         x1=x2;
         y1=y2;
         x2=pos->second.vx[i];
         y2=pos->second.vy[i];
         this->Data2Screen(x2,y2);
         if(  ((x1>=mLeft)&&(x1<=(width-mRight))&&(y1>=mBottom)&&(y1<=(height-mTop)))
            ||((x2>=mLeft)&&(x2<=(width-mRight))&&(y2>=mBottom)&&(y2<=(height-mTop))))
            dc.DrawLine(wxCoord(x1),wxCoord(y1),wxCoord(x2),wxCoord(y2));
      }
      // Print Name
      dc.SetTextForeground(wxPen(wxTheColourDatabase->Find(wxString::FromAscii(swxColourNameList[ix])),1,wxSOLID).GetColour());
      wxCoord tmpW,tmpH;
      fontInfo.Printf(wxString::FromAscii(pos->second.name.c_str()));
      dc.GetTextExtent(fontInfo, &tmpW, &tmpH);
      dc.DrawText(fontInfo,wxCoord(width-tmpW-2),wxCoord(tmpH*(ixdrawn)+2));
   }

   VFN_DEBUG_EXIT("WXMultiGraph::OnPaint()",4)
}

void WXMultiGraph::OnMouse(wxMouseEvent &event)
{
   event.Skip();// Make sure the default handler gets the event too
   
   if(event.Leaving()) return;// ?
   wxCoord width,height;
   this->GetSize(&width, &height);
   // Write mouse pointer coordinates
      wxClientDC dc(this);
      PrepareDC(dc);
      mpParentFrame->PrepareDC(dc);

      wxPoint pos=event.GetPosition();
      float x= float(dc.DeviceToLogicalX(pos.x));
      float y= float(dc.DeviceToLogicalY(pos.y));

      if((x>width)||(y>height))
      {
         return;
      }
      this->Screen2Data(x,y);
      wxString str;
      str.Printf(_T("x=%f    ,y=%f"),x,y);
      mpParentFrame->SetStatusText(str);

   if(event.RightIsDown())
   {
      this->PopupMenu(mpPopUpMenu, event.GetX(), event.GetY() );
      return;
   }
   if (event.Dragging() && event.LeftIsDown() && (!mIsDragging))
   {//Begin zooming
      mIsDragging=true;
      mDragX0=x;
      mDragY0=y;
      return;
   }
   if(event.LeftUp() && mIsDragging)
   {//Finished zooming !
      mMutexData.Lock();
      mIsDragging=false;
      if(x>mDragX0)
      {
         mMinX=mDragX0;
         mMaxX=x;
      }
      else
      {
         mMinX=x;
         mMaxX=mDragX0;
      }
      if(y>mDragY0)
      {
         mMinY=mDragY0;
         mMaxY=y;
      }
      else
      {
         mMinY=y;
         mMaxY=mDragY0;
      }
      if(mMaxX<=mMinX) mMaxX=mMinX+1e-6;
      if(mMaxY<=mMinY) mMaxY=mMinY+1e-6;
      mMutexData.Unlock();
      this->UpdateDisplay();
      return;
   }
   if(false==event.Dragging()) mIsDragging=false;

   if(event.LeftDClick())
   {//Reset axis range
      this->AutoScale();
      this->UpdateDisplay();
      return;
   }
}

void WXMultiGraph::OnMouseWheel(wxMouseEvent &event)
{
   if(event.GetWheelRotation()>=event.GetWheelDelta())
   {
      const REAL range=mMaxX-mMinX;
      mMaxX += range/8;
      mMinX += range/8;
   }
   if(event.GetWheelRotation()<=(-event.GetWheelDelta()))
   {
      const REAL range=mMaxX-mMinX;
      mMaxX -= range/8;
      mMinX -= range/8;
   }
   this->UpdateDisplay();
}

void WXMultiGraph::AutoScale(const long id,const bool xmin,const bool xmax,
                                           const bool ymin,const bool ymax)
{
   wxMutexLocker mlock(mMutexData);
   if(mvData.size()==0) return;
   std::map<unsigned long, GraphData>::const_iterator pos;
   if(id<0)pos=mvData.end();
   else pos=mvData.find((unsigned long)id);
   if(pos==mvData.end())
   {
      pos=mvData.begin();
      if(xmax) mMaxX=pos->second.xmax;
      if(xmin) mMinX=pos->second.xmin;
      if(ymax) mMaxY=pos->second.ymax;
      if(ymin) mMinY=pos->second.ymin;
      ++pos;
      for(;pos!=mvData.end();++pos)
      {
         if(xmax) if(mMaxX<pos->second.xmax) mMaxX=pos->second.xmax;
         if(xmin) if(mMinX>pos->second.xmin) mMinX=pos->second.xmin;
         if(ymax) if(mMaxY<pos->second.ymax) mMaxY=pos->second.ymax;
         if(ymin) if(mMinY>pos->second.ymin) mMinY=pos->second.ymin;
      }
   }
   else
   {
      if(xmax) mMaxX=pos->second.xmax;
      if(xmin) mMinX=pos->second.xmin;
      if(ymax) mMaxY=pos->second.ymax;
      if(ymin) mMinY=pos->second.ymin;
   }
   if(mMaxX<=mMinX) mMaxX=mMinX+1e-6;
   if(mMaxY<=mMinY) mMaxY=mMinY+1e-6;
   //cout<<"Autoscale to:"<<mMinX<<"->"<<mMaxX<<" , "<<mMinY<<"->"<<mMaxY<<endl;
}

void WXMultiGraph::OnKeyDown(wxKeyEvent& event)
{
   switch(event.GetKeyCode())
   {
      case(WXK_LEFT):
      {
         const REAL range=mMaxX-mMinX;
         mMinX-=range/8;
         mMaxX-=range/8;
         break;
      }
      case(WXK_RIGHT):
      {
         const REAL range=mMaxX-mMinX;
         mMinX+=range/8;
         mMaxX+=range/8;
         break;
      }
      case(WXK_UP):
      {
         const REAL range=mMaxY-mMinY;
         mMinY+=range/8;
         mMaxY+=range/8;
         break;
      }
      case(WXK_DOWN):
      {
         const REAL range=mMaxY-mMinY;
         mMinY-=range/8;
         mMaxY-=range/8;
         break;
      }
      case(43):// WXK_ADD ?
      {
         const REAL range=mMaxX-mMinX;
         const REAL center=(mMaxX+mMinX)/2;
         mMinX=center-range/3.0;
         mMaxX=center+range/3.0;
         break;
      }
      case(45):// WXK_SUBTRACT ?
      {
         const REAL range=mMaxX-mMinX;
         const REAL center=(mMaxX+mMinX)/2;
         mMinX=center-range*2.0/3.0;
         mMaxX=center+range*2.0/3.0;
         break;
      }
      case(42):// WXK_MULTIPLY
      {
         const REAL range=mMaxY-mMinY;
         mMaxY=mMinY+range*2.0/3.0;
         break;
      }
      case(47):// WXK_DIVIDE
      {
         const REAL range=mMaxY-mMinY;
         mMaxY=mMinY+range*4.0/3.0;
         break;
      }
      default: 
      {
         VFN_DEBUG_MESSAGE("WXMultiGraph::OnKeyDown(): no command for key #"<<event.GetKeyCode(),5);
      }
   }
   this->UpdateDisplay();
   event.Skip();// Make sure the default handler gets the event too
}

void WXMultiGraph::OnUpdateUI(wxUpdateUIEvent &event)
{
   VFN_DEBUG_MESSAGE("WXMultiGraph::OnUpdateUI()",4)
   this->Refresh(false);
   event.Skip();
}

void WXMultiGraph::OnSize(wxSizeEvent &event)
{
   this->Refresh(false);
}
void WXMultiGraph::SetXLabel(const wxString &xlabel)
{
   mXLabel=xlabel;
}

void WXMultiGraph::SetYLabel(const wxString &ylabel)
{
   mYLabel=ylabel;
}

void WXMultiGraph::UpdateDisplay()
{
   VFN_DEBUG_ENTRY("WXMultiGraph::UpdateDisplay()",4)
   if(wxThread::IsMain()) this->Refresh(false);
   else
   {
      wxUpdateUIEvent event(ID_UPDATEUI);
      wxPostEvent(this,event);
   }
   VFN_DEBUG_EXIT("WXMultiGraph::UpdateDisplay()",4)
}


void WXMultiGraph::Screen2Data(float &x,float &y)
{
   wxCoord width,height;
   this->GetSize(&width, &height);
   
   float range=float(width-(mLeft+mRight));
   if(range<=0) range=1.0;
   x=(x-mLeft)/range;
   x=mMinX+x*(mMaxX-mMinX);
   
   range=float(height-(mTop+mBottom));
   if(range<=0) range=1.0;
   y=(height-mBottom-y)/range;
   y=mMinY+y*(mMaxY-mMinY);
}

void WXMultiGraph::Data2Screen(float &x,float &y)
{
   VFN_DEBUG_ENTRY("WXMultiGraph::Data2Screen()"<<x<<","<<y,2)
   wxCoord width,height;
   this->GetSize(&width, &height);

   float range=float(width-(mLeft+mRight));
   x=(x-mMinX)/(mMaxX-mMinX);
   x=mLeft+x*range;

   range=float(height-(mTop+mBottom));
   y=(y-mMinY)/(mMaxY-mMinY);
   y=height-mBottom-y*range;
   VFN_DEBUG_EXIT("WXMultiGraph::Data2Screen()->"<<x<<","<<y,2)
}

}//namespace
