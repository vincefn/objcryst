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
*  header file for the wxLSQ class
*
*/

#ifndef _VFN_WX_LSQ_H_
#define _VFN_WX_LSQ_H_

#include "RefinableObj/LSQNumObj.h"
#include "wxCryst/wxCryst.h"

namespace ObjCryst
{
class WXLSQ: public WXCrystObj
{
   public:
      WXLSQ(wxWindow *parent, LSQNumObj*);
      virtual ~WXLSQ();
      virtual bool OnChangeName(const int id);
   private:
      LSQNumObj* mpLSQ;
   DECLARE_EVENT_TABLE()
};

} //namespace

#endif //_VFN_WX_ATOM_H_
