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

#ifndef _VFN_WX_MOLECULE_H_
#define _VFN_WX_MOLECULE_H_

#include "wxCryst/wxScatterer.h"
#include "ObjCryst/Molecule.h"

namespace ObjCryst
{

/// wx class for MolAtom objects
class WXMolAtom:public WXCrystObjBasic
{
   public:
      WXMolAtom(wxWindow *parent, MolAtom*);
      virtual void CrystUpdate();
      virtual void UpdateUI();
      virtual bool Layout();
      void OnChangeScattPow(wxCommandEvent &);
   private:
      MolAtom *mpMolAtom;
      wxBoxSizer *mpSizer;
      WXCrystObjBasicList mList;
      WXFieldString *mpFieldName;
      WXFieldChoice* mpFieldScattPower;
   DECLARE_EVENT_TABLE()
};

/// wx class for MolBond objects
class WXMolBond:public WXCrystObjBasic
{
   public:
      WXMolBond(wxWindow *parent, MolBond*);
      virtual void CrystUpdate();
      virtual void UpdateUI();
      virtual bool Layout();
      void OnChangeAtom(wxCommandEvent &);
   private:
      MolBond *mpMolBond;
      wxBoxSizer *mpSizer;
      WXCrystObjBasicList mList;
      WXFieldChoice* mpFieldAtom1;
      WXFieldChoice* mpFieldAtom2;
   DECLARE_EVENT_TABLE()
};

/// wx class for MolBondAngle objects
class WXMolBondAngle:public WXCrystObjBasic
{
   public:
      WXMolBondAngle(wxWindow *parent, MolBondAngle*);
      virtual void CrystUpdate();
      virtual void UpdateUI();
      virtual bool Layout();
      void OnChangeAtom(wxCommandEvent &);
   private:
      MolBondAngle *mpMolBondAngle;
      wxBoxSizer *mpSizer;
      WXCrystObjBasicList mList;
      WXFieldChoice* mpFieldAtom1;
      WXFieldChoice* mpFieldAtom2;
      WXFieldChoice* mpFieldAtom3;
   DECLARE_EVENT_TABLE()
};

/// wx class for MolDihedralAngle objects
class WXMolDihedralAngle:public WXCrystObjBasic
{
   public:
      WXMolDihedralAngle(wxWindow *parent, MolDihedralAngle*);
      virtual void CrystUpdate();
      virtual void UpdateUI();
      virtual bool Layout();
      void OnChangeAtom(wxCommandEvent &);
   private:
      MolDihedralAngle *mpMolDihedralAngle;
      wxBoxSizer *mpSizer;
      WXCrystObjBasicList mList;
      WXFieldChoice* mpFieldAtom1;
      WXFieldChoice* mpFieldAtom2;
      WXFieldChoice* mpFieldAtom3;
      WXFieldChoice* mpFieldAtom4;
   DECLARE_EVENT_TABLE()
};

/// wxCryst class for Molecule objects
class WXMolecule: public WXScatterer
{
   public:
      WXMolecule(wxWindow *parent, Molecule*);
      virtual ~WXMolecule();
      void OnMenuOptimizeConformation(wxCommandEvent & WXUNUSED(event));
      void OnMenuAddAtom(wxCommandEvent & WXUNUSED(event));
      void OnMenuAddBond(wxCommandEvent & WXUNUSED(event));
      void OnMenuAddAngle(wxCommandEvent & WXUNUSED(event));
      void OnMenuAddDihedralAngle(wxCommandEvent & WXUNUSED(event));
      void OnMenuRemoveAtom(wxCommandEvent & WXUNUSED(event));
      void OnMenuRemoveBond(wxCommandEvent & WXUNUSED(event));
      void OnMenuRemoveAngle(wxCommandEvent & WXUNUSED(event));
      void OnMenuRemoveDihedralAngle(wxCommandEvent & WXUNUSED(event));
      void OnMenuSetLimits(wxCommandEvent &event);
      void OnMenuTest(wxCommandEvent &event);
      virtual void CrystUpdate();
      virtual void UpdateUI();
   private:
      Molecule* mpMolecule;
      wxBoxSizer* mpSizerAtomList;
      wxBoxSizer* mpSizerBondList;
      wxBoxSizer* mpSizerAngleList;
      wxBoxSizer* mpSizerDihedralAngleList;
      /** Displayed list of atoms
      */
      vector<MolAtom*> mvpAtom;
      /** Displayed list of bonds
      */
      vector<MolBond*> mvpBond;
      /** Displayed list of bond angle
      */
      vector<MolBondAngle*> mvpBondAngle;
      /** Displayed list of Dihedral angles
      */
      vector<MolDihedralAngle*> mvpDihedralAngle;
   DECLARE_EVENT_TABLE()
};

} //namespace

#endif //_VFN_WX_MOLECULE_H_
