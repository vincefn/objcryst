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
/*
*  header file for the RefinablePar and RefinableObj classes
*
* This is still in early development stages !! Not secure !
*
*/

#ifndef _VFN_WX_REFINABLE_OBJ_H_
#define _VFN_WX_REFINABLE_OBJ_H_

namespace ObjCryst
{
template<class T> class WXRegistry;
class WXFieldOption;
class WXFieldRefPar;
} //namespace

#include "wxCryst/wxCryst.h"

// We should not have to put this here. !! :TODO:
#include "RefinableObj/RefinableObj.h"

namespace ObjCryst
{
class RefinablePar;
class RefinableObj;
/// This function allows to pick up one object in a registry. Shows a modal
/// dialog with a list of objects
template<class T> T* WXDialogChooseFromRegistry(ObjRegistry<T> &reg,wxWindow*parent,
                                                const string &message,int &);
/// This function allows to pick up one object in a registry. Shows a modal
/// dialog with a list of objects. This is a const access.
template<class T> const T* WXDialogChooseFromRegistry(const ObjRegistry<T> &reg,
                                                      wxWindow*parent, const string &message,
                                                      int &);

/// A field for a RefinablePar. This shows the 'human' value of the parameter, and allows
/// the modification of the parameter. A button allows to fix/unfix the parameter.
/// \todo: allow acces to the parameters limits
class WXFieldRefPar:public WXField
{
   public:
      WXFieldRefPar(wxWindow *parent,const string& label, 
                    RefinablePar *refpar,const int hsize=50,
                    const bool enableFixButton=true, const bool enableLimitedButton=true);
      /// When a new value is entered (must type it and then hit the 'enter' key).
      /// The Field reads the new value, 
      /// and directly changes the RefinablePar value (contrary to what happens
      /// for WXFieldName)by using RefinablePar::SetHumanValue().
      ~WXFieldRefPar();
      void OnEnter(wxCommandEvent & WXUNUSED(event));
      /// Records when text is entered (either from self-updating or user input)
      void OnText(wxCommandEvent & WXUNUSED(event));
      /// Toggle the 'fixed' status of the parameter.
      void OnToggleFix(wxCommandEvent & WXUNUSED(event));
      /// Toggle the 'limited' status of the parameter.
      void OnToggleLimited(wxCommandEvent & WXUNUSED(event));
      /// Opens the popu menu, to allow changing limits
      void OnPopupMenu(wxMouseEvent & event);
      /// Opens the popu menu, to allow changing limits
      void OnPopupMenuChoice(wxCommandEvent& event);
      virtual void CrystUpdate();
      virtual void UpdateUI();
      /// Get the RefinablePar associated to this field
      RefinablePar& GetRefPar();
      void Revert();
      virtual void ValidateUserInput();
   protected:
      REAL mValue;
      wxCheckBox *mpButtonFix;
      wxCheckBox *mpButtonLimited;
      wxTextCtrl *mpField;
      RefinablePar *mpRefPar;
      REAL mValueOld;
      bool mIsSelfUpdating;
   DECLARE_EVENT_TABLE()
};

class RefObjOpt;// Declared in RefinableObj.h
/// WX representation of a RefObj option. This displays the names of the different choices.
class WXFieldOption:public WXField
{
   public:
      WXFieldOption(wxWindow *parent,const int field_id,
                    RefObjOpt* option);
      virtual ~WXFieldOption();
      void OnUpdateUI(wxUpdateUIEvent & WXUNUSED(event));
      /// When a new value is entered. The Field reads the new value, then
      /// forwards the event to its owner, who will take care of anything
      /// that must be done.
      void OnChoice(wxCommandEvent & WXUNUSED(event));
      virtual void CrystUpdate();
      virtual void UpdateUI();
      void Revert();
      /// Does nothing. Any user input is directly validated (OnChoice).
      virtual void ValidateUserInput();
   protected:
      int mChoice;
      int mChoiceOld;
      RefObjOpt* mpOption;
      wxChoice *mpList;
   DECLARE_EVENT_TABLE()
};

/// This displays all components of a ObjCryst++ Registry.
template<class T> class WXRegistry:public WXCrystObj
{
   public:
      WXRegistry(wxWindow *parent,ObjRegistry<T>* reg);
      ~WXRegistry();
      void Add(WXCrystObjBasic *obj);
      void Remove(WXCrystObjBasic *obj);
      virtual bool OnChangeName(const int id);
   private:
      ObjRegistry<T> *mpRegistry;
};
                                                      
/// The base wxCryst class for all RefinableObj objects. This shows the title,
/// a menu for XMLInput/XMLOutput, and all RefObjOpt.
class WXRefinableObj: public WXCrystObj
{
   public:
      WXRefinableObj(wxWindow *parent, RefinableObj*);
      ~WXRefinableObj();
      virtual void CrystUpdate();
      virtual void UpdateUI();
      virtual bool OnChangeName(const int id);
      void OnMenuSave(wxCommandEvent & WXUNUSED(event));
      void OnMenuLoad(wxCommandEvent & WXUNUSED(event));
      void OnMenuFixAllPar(wxCommandEvent & WXUNUSED(event));
      void OnMenuUnFixAllPar(wxCommandEvent & WXUNUSED(event));
      void OnMenuParRandomize(wxCommandEvent & WXUNUSED(event));
      virtual void OnUpdateUI(wxUpdateUIEvent& event);
   protected:
      WXCrystMenuBar* mpMenuBar;
   private:
      RefinableObj* mpRefinableObj;
   DECLARE_EVENT_TABLE()
};
} //namespace

#endif //_VFN_WX_REFINABLE_OBJ_H_
