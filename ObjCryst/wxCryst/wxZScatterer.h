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

#ifndef _VFN_WX_ZSCATTERER_H_
#define _VFN_WX_ZSCATTERER_H_

#include "wxCryst/wxScatterer.h"
#include "ObjCryst/ZScatterer.h"

namespace ObjCryst
{

class WXZAtom:public WXCrystObjBasic
{
   public:
      WXZAtom(wxWindow *parent, ZAtom*);
      virtual void CrystUpdate();
      virtual void UpdateUI();
      virtual bool Layout();
      void OnChangeScattPow(wxCommandEvent & WXUNUSED(event));
   private:
      ZAtom *mpZAtom;
      wxBoxSizer *mpSizer;
      WXCrystObjBasicList mList;
      WXFieldString *mpFieldName;
      WXFieldChoice* mpFieldScattPower;
   DECLARE_EVENT_TABLE()
};

/// wxCryst class for ZScatterer objects
class WXZScatterer: public WXScatterer
{
   public:
      WXZScatterer(wxWindow *parent, ZScatterer*);
      void OnMenuAddZAtom(wxCommandEvent & WXUNUSED(event));
      void OnMenuSetLimits(wxCommandEvent &event);
      void OnMenuChangePivotAtom(wxCommandEvent &WXUNUSED(event));
      void OnMenuImportZMatrix(wxCommandEvent &WXUNUSED(event));
      void OnMenuExportZMatrix(wxCommandEvent &WXUNUSED(event));
   private:
      ZScatterer* mpZScatterer;
      WXRegistry<ZAtom> *mpWXZAtomRegistry;
   DECLARE_EVENT_TABLE()
};

} //namespace

#endif //_VFN_WX_ZSCATTERER_H_
