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
//#include <sstream> //for stringstream
#include <fstream>

#include "wx/wx.h"
#include "wxCryst/wxScatterer.h"

//Fixes for Cygwin; where do those stupid macros come from ? Somewhere in wxMSW headers
#ifdef max
#undef max
#endif
#ifdef min
#undef min
#endif
#ifdef DrawText
#undef DrawText
#endif
 
namespace ObjCryst
{
////////////////////////////////////////////////////////////////////////
//
//    WXScatterer
//
////////////////////////////////////////////////////////////////////////
WXScatterer::WXScatterer(wxWindow* parent, Scatterer *obj):
WXRefinableObj(parent,(RefinableObj*)obj),mpScatterer(obj)
{
   VFN_DEBUG_MESSAGE("WXScatterer::WXScatterer()",6)
   mpWXTitle->SetForegroundColour(wxColour(0,100,0));
   //Lattice
      wxBoxSizer* sizer=new wxBoxSizer(wxHORIZONTAL);
      mpScatterer->RefinableObj::Print();
      WXCrystObjBasic* mpFieldX
         =mpScatterer->GetPar(mpScatterer->mXYZ.data()+0).WXCreate(this);

      WXCrystObjBasic* mpFieldY
         =mpScatterer->GetPar(mpScatterer->mXYZ.data()+1).WXCreate(this);

      WXCrystObjBasic* mpFieldZ
         =mpScatterer->GetPar(mpScatterer->mXYZ.data()+2).WXCreate(this);
          
      WXCrystObjBasic* mpFieldPopu 
         =mpScatterer->GetPar(&(mpScatterer->mOccupancy)).WXCreate(this);

      sizer->Add(mpFieldX    ,0,wxALIGN_CENTER);
      sizer->Add(mpFieldY    ,0,wxALIGN_CENTER);
      sizer->Add(mpFieldZ    ,0,wxALIGN_CENTER);
      sizer->Add(mpFieldPopu ,0,wxALIGN_CENTER);
      
      mpSizer->Add(sizer,0,wxALIGN_LEFT);
      mList.Add(mpFieldX);
      mList.Add(mpFieldY);
      mList.Add(mpFieldZ);
      mList.Add(mpFieldPopu);
   
   //mWXParent->Layout();
   mpSizer->Layout();
   mpTopSizer->Layout();
   mpTopSizer->Fit(this);
   this->CrystUpdate();
}

}// namespace 

