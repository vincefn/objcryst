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

// wx headers, with or without precompilation
#include "wx/wxprec.h"
#ifdef __BORLANDC__
    #pragma hdrstop
#endif
#ifndef WX_PRECOMP
    #include "wx/wx.h"
#endif

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
#if 1
      WXFieldRefPar* pFieldX    =new WXFieldRefPar(this,"x:",
                                     &(mpScatterer->GetPar(mpScatterer->mXYZ.data()+0)) );
      WXFieldRefPar* pFieldY    =new WXFieldRefPar(this,"y:",
                                     &(mpScatterer->GetPar(mpScatterer->mXYZ.data()+1)) );
      WXFieldRefPar* pFieldZ    =new WXFieldRefPar(this,"z:",
                                     &(mpScatterer->GetPar(mpScatterer->mXYZ.data()+2)) );
      WXFieldRefPar* pFieldPopu    =new WXFieldRefPar(this,"Occup:",
                                     &(mpScatterer->GetPar(&(mpScatterer->mOccupancy))) );
#else
      WXCrystObjBasic* pFieldX
         =mpScatterer->GetPar(mpScatterer->mXYZ.data()+0).WXCreate(this);
      WXCrystObjBasic* pFieldY
         =mpScatterer->GetPar(mpScatterer->mXYZ.data()+1).WXCreate(this);
      WXCrystObjBasic* pFieldZ
         =mpScatterer->GetPar(mpScatterer->mXYZ.data()+2).WXCreate(this);
      WXCrystObjBasic* pFieldPopu 
         =mpScatterer->GetPar(&(mpScatterer->mOccupancy)).WXCreate(this);
#endif
      sizer->Add(pFieldX    ,0,wxALIGN_CENTER);
      sizer->Add(pFieldY    ,0,wxALIGN_CENTER);
      sizer->Add(pFieldZ    ,0,wxALIGN_CENTER);
      sizer->Add(pFieldPopu ,0,wxALIGN_CENTER);
      
      mpSizer->Add(sizer,0,wxALIGN_LEFT);
      mList.Add(pFieldX);
      mList.Add(pFieldY);
      mList.Add(pFieldZ);
      mList.Add(pFieldPopu);
   //layout & update is done in derived objects
}

}// namespace 

