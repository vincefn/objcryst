/*  ObjCryst++ Object-Oriented Crystallographic Library
    (c) 2008 Vincent Favre-Nicolin vincefn@users.sourceforge.net

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
*  code file for the wxLSQ class
*
*/

#include "ObjCryst/wxCryst/wxLSQ.h"

namespace ObjCryst
{

BEGIN_EVENT_TABLE(WXLSQ, wxWindow)
END_EVENT_TABLE()

WXLSQ::WXLSQ(wxWindow *parent, LSQNumObj*p):
WXCrystObj(parent, wxVERTICAL, false),mpLSQ(p)
{
   // Add all individual refinable parameters
   for(unsigned int i=0;i<mpLSQ->GetCompiledRefinedObj().GetNbPar();++i)
   {
      if(mpLSQ->GetCompiledRefinedObj().GetPar(i).IsUsed())
      {
         WXCrystObjBasic* pRefPar=mpLSQ->GetCompiledRefinedObj().GetPar(i).WXCreate(this);
         mpSizer->Add(pRefPar, 1, wxALIGN_RIGHT);
         mList.Add(pRefPar);
      }
   }
   this->CrystUpdate(true,true);
}

WXLSQ::~WXLSQ(){mpLSQ->WXNotifyDelete();}

bool WXLSQ::OnChangeName(const int id)
{
   return false;
}

} //namespace
