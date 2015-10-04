/*  ObjCryst++ Object-Oriented Crystallographic Library
    (c) 2000-2007 Vincent Favre-Nicolin vincefn@users.sourceforge.net
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

#include "wx/grid.h"
#include "ObjCryst/wxCryst/wxScatterer.h"
#include "ObjCryst/ObjCryst/Molecule.h"

namespace ObjCryst
{

/** Storage of a Molecule description
*/
struct MoleculeFragment
{
    /// Constructor
    MoleculeFragment(Molecule * pmol,std::string formula);
    /// Points to the real Molecule object as long as it resides in memory and is displayed
    Molecule *mpMolecule;
    /// The formula
    std::string mFormula;
    /// The XML description as a string
    std::string mXML;
    #ifdef __DEBUG__
    /// Count how many times we updated this (for testing)
    unsigned long ct;
    #endif
};
/** Global storage for molecule descriptions. Each molecule used is stored as a
* stringstream of its XML description. It can then be used in another Crystal.
* The molecule descriptions are stored as long as the GUI application runs,
* but unless a large number of Molecules are displayed this should not use
* too much memory.
*/
extern std::list<MoleculeFragment> gvMoleculeFragment;

class WXMolecule;
// Scrolled window for bonds, bond angles and dihedral angles
class WXMolScrolledWindow:public wxGrid
{
   public:
      WXMolScrolledWindow(wxWindow* parent, WXMolecule* pWXMol,long id=-1);
      virtual ~WXMolScrolledWindow();
   private:
      /// The WXMolecule window which created this window, and who should be told
      /// if it is destroyed.
      WXMolecule* mpWXMolecule;
};

/// wx class for MolAtom objects
class WXMolAtom:public WXCrystObjBasic
{
   public:
      WXMolAtom(wxWindow *parent, MolAtom*);
      virtual ~WXMolAtom();
      virtual void CrystUpdate(const bool updateUI=false,const bool mutexlock=false);
      virtual void UpdateUI(const bool mutexlock=false);
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
      virtual ~WXMolBond();
      virtual void CrystUpdate(const bool updateUI=false,const bool mutexlock=false);
      virtual void UpdateUI(const bool mutexlock=false);
      void OnChangeAtom(wxCommandEvent &);
      /// Toggle the 'free' status of the bond.
      void OnToggleFree(wxCommandEvent & WXUNUSED(event));
   private:
      MolBond *mpMolBond;
      wxBoxSizer *mpSizer;
      WXCrystObjBasicList mList;
      WXFieldChoice* mpFieldAtom1;
      WXFieldChoice* mpFieldAtom2;
      wxCheckBox *mpButtonFree;
      /// The current value
      REAL mValue;
   DECLARE_EVENT_TABLE()
};

/// wx class for MolBondAngle objects
class WXMolBondAngle:public WXCrystObjBasic
{
   public:
      WXMolBondAngle(wxWindow *parent, MolBondAngle*);
      virtual ~WXMolBondAngle();
      virtual void CrystUpdate(const bool updateUI=false,const bool mutexlock=false);
      virtual void UpdateUI(const bool mutexlock=false);
      void OnChangeAtom(wxCommandEvent &);
   private:
      MolBondAngle *mpMolBondAngle;
      wxBoxSizer *mpSizer;
      WXCrystObjBasicList mList;
      WXFieldChoice* mpFieldAtom1;
      WXFieldChoice* mpFieldAtom2;
      WXFieldChoice* mpFieldAtom3;
      /// The current value
      REAL mValue;
   DECLARE_EVENT_TABLE()
};

/// wx class for MolDihedralAngle objects
class WXMolDihedralAngle:public WXCrystObjBasic
{
   public:
      WXMolDihedralAngle(wxWindow *parent, MolDihedralAngle*);
      virtual ~WXMolDihedralAngle();
      virtual void CrystUpdate(const bool updateUI=false,const bool mutexlock=false);
      virtual void UpdateUI(const bool mutexlock=false);
      void OnChangeAtom(wxCommandEvent &);
   private:
      MolDihedralAngle *mpMolDihedralAngle;
      wxBoxSizer *mpSizer;
      WXCrystObjBasicList mList;
      WXFieldChoice* mpFieldAtom1;
      WXFieldChoice* mpFieldAtom2;
      WXFieldChoice* mpFieldAtom3;
      WXFieldChoice* mpFieldAtom4;
      /// The current value
      REAL mValue;
   DECLARE_EVENT_TABLE()
};

/// wxCryst class for Molecule objects
class WXMolecule: public WXScatterer
{
   public:
      WXMolecule(wxWindow *parent, Molecule*);
      virtual ~WXMolecule();
      void OnMenuOptimizeConformation(wxCommandEvent & WXUNUSED(event));
      void OnMenuPrintRestraintStatus(wxCommandEvent & WXUNUSED(event));
      void OnMenuExportRestraints(wxCommandEvent & WXUNUSED(event));
      void OnMenuAddAtom(wxCommandEvent & WXUNUSED(event));
      void OnMenuAddBond(wxCommandEvent & WXUNUSED(event));
      void OnMenuAddAngle(wxCommandEvent & WXUNUSED(event));
      void OnMenuAddDihedralAngle(wxCommandEvent & WXUNUSED(event));
      void OnMenuAddRigidGroup(wxCommandEvent & WXUNUSED(event));
      void OnMenuAddNonFlipAtom(wxCommandEvent & WXUNUSED(event));
      void OnMenuRigidfyWithDihedralAngles(wxCommandEvent & WXUNUSED(event));
      void OnMenuRemoveAtom(wxCommandEvent & WXUNUSED(event));
      void OnMenuRemoveBond(wxCommandEvent & WXUNUSED(event));
      void OnMenuRemoveAngle(wxCommandEvent & WXUNUSED(event));
      void OnMenuRemoveDihedralAngle(wxCommandEvent & WXUNUSED(event));
      void OnMenuRemoveNonFlipAtom(wxCommandEvent & WXUNUSED(event));
      void OnMenuRemoveRigidGroup(wxCommandEvent & WXUNUSED(event));   
      void OnMenuSetLimits(wxCommandEvent &event);
      void OnMenuShowRestraintWindow(wxCommandEvent &event);
      void OnMenuSetDeltaSigma(wxCommandEvent &event);
      void OnChangeCenterAtom(wxCommandEvent &event);
      void OnEditGridAtom(wxGridEvent &e);
      void OnEditGridBondLength(wxGridEvent &e);
      void OnEditGridBondAngle(wxGridEvent &e);
      void OnEditGridDihedralAngle(wxGridEvent &e);
      void OnEditGridRigidGroup(wxGridEvent &e);
      void OnMenuExport2ZMatrix(wxCommandEvent &event);
      void OnMenuTest(wxCommandEvent &event);
      void OnMenuMDTest(wxCommandEvent &event);
      void OnMenuRotate(wxCommandEvent &event);
      /// Notify that either the bond, bond angle or dihedral angle list window has
      /// been destroyed
      void NotifyDeleteListWin(WXMolScrolledWindow *win);
      virtual void CrystUpdate(const bool updateUI=false,const bool mutexlock=false);
      virtual void UpdateUI(const bool mutexlock=false);
      virtual bool Enable(bool enable=true);
   private:
      Molecule* mpMolecule;
      WXMolScrolledWindow* mpAtomWin;
      WXMolScrolledWindow* mpBondWin;
      WXMolScrolledWindow* mpAngleWin;
      WXMolScrolledWindow* mpDihedralAngleWin;
      WXMolScrolledWindow* mpRigidGroupWin;
      WXMolScrolledWindow* mpNonFlipAtomWin;
      /// Structure to store the Atom parameters
      struct CellAtom
      {
         CellAtom();
         MolAtom* mpAtom;
         std::string mName;
         const ScatteringPower* mpScatteringPower;
         REAL mX,mY,mZ,mOcc;
         /// True if we need to update the displayed values
         bool mNeedUpdateUI;
      };
      /** Displayed list of atoms
      */
      std::list<CellAtom> mvpAtom;
      
      /// Structure to store the bond current values
      struct CellBond
      {
         CellBond();
         MolBond* mpBond;
         std::string mAtom1;
         std::string mAtom2;
         REAL mLength;
         REAL mLength0;
         REAL mSigma;
         REAL mDelta;
         /// True if we need to update the displayed values
         bool mNeedUpdateUI;
      };
      /** Displayed list of bonds, in the order they appear
      */
      std::list<CellBond> mvpBond;
      /// Structure to store the bond angles current values
      struct CellBondAngle
      {
         CellBondAngle();
         MolBondAngle* mpBondAngle;
         std::string mAtom1;
         std::string mAtom2;
         std::string mAtom3;
         REAL mAngle;
         REAL mAngle0;
         REAL mSigma;
         REAL mDelta;
         /// True if we need to update the displayed values
         bool mNeedUpdateUI;
      };
      /** Displayed list of bond angle
      */
      std::list<CellBondAngle> mvpBondAngle;
      /// Structure to store the dihedral angles current values
      struct CellDihedralAngle
      {
         CellDihedralAngle();
         MolDihedralAngle* mpDihedralAngle;
         std::string mAtom1;
         std::string mAtom2;
         std::string mAtom3;
         std::string mAtom4;
         REAL mAngle;
         REAL mAngle0;
         REAL mSigma;
         REAL mDelta;
         /// True if we need to update the displayed values
         bool mNeedUpdateUI;
      };
      /** Displayed list of Dihedral angles
      */
      std::list<CellDihedralAngle> mvpDihedralAngle;
      struct CellRigidGroup
      {
         CellRigidGroup();
         /// Rigid group in the Molecule
         RigidGroup *mpGroup;
         /// Copy of the set of atoms, as it was last displayed
         RigidGroup mGroupCopy;
         /// True if we need to update the displayed values
         bool mNeedUpdateUI;
      };
      /** Displayed list of Dihedral angles
      */
      std::list<CellRigidGroup> mvpRigidGroup;
      /// Flag to indicate whether we are updating values in the wxGrid data.
      /// (enabled in wxMolecule::UpdateUI()).
      bool mIsSelfUpdating;
      /// Center atom
      WXFieldChoice* mpFieldCenterAtom;
   DECLARE_EVENT_TABLE()
};

} //namespace

#endif //_VFN_WX_MOLECULE_H_
