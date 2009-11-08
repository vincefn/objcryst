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
*  header file for the general wxCryst functions.
*
*
*/

#ifndef _VFN_WX_CRYST_H_
#define _VFN_WX_CRYST_H_

#include <string>
#include <iostream>
#include <set>
#include <map>

// wx headers, with or without precompilation
#include "wx/wxprec.h"
#ifdef __BORLANDC__
    #pragma hdrstop
#endif
#ifndef WX_PRECOMP
    #include "wx/wx.h"
#endif

#include "CrystVector/CrystVector.h"

namespace ObjCryst
{
#if 1
/** Provides the same functionnality as wxMultiChoiceDialog, but always using
* a wxListBox, which is much easier when selecting a large number of successive
* choices (using shift-click).
*/
class wxMultiChoiceDialog_ListBox:public wxDialog
{
  public:
    wxMultiChoiceDialog_ListBox(wxWindow* parent, const wxString& message, const wxString& caption, 
                                int n, const wxString* choices);
    wxArrayInt GetSelections() const;
  private:
    wxListBox mListBox;
};
#endif
////////////////////////////////////////////////////////////////////////
//
// Unique ID for menus incrementer
//
////////////////////////////////////////////////////////////////////////
/** Class to automatically assign a unique wxID to each window
*
*/
class WXCRYST_ID
{
   public:
      WXCRYST_ID();
      operator long();
   private:
      long mIndex;
      static long mCounter;
};
extern const long ID_WXOBJ_ENABLE; //These are used in ObjCryst/RefinableObj.cpp
extern const long ID_WXOBJ_DISABLE;// and defined in wxCryst/wxCryst.cpp

#undef VFN_CRYST_MUTEX
#ifdef VFN_CRYST_MUTEX
/** Derived Mutex class, to keep track of Mutexes
*
*/
class CrystMutex:public wxMutex
{
   public:
      CrystMutex();
      ~CrystMutex();
      wxMutexError Lock();
      wxMutexError Unlock();
   private:
      unsigned long mNbLock;
};
#else
#define CrystMutex wxMutex
#endif
#ifndef DOXYGEN_SKIP_THIS
// constants that are shared between several classes must
// be defined here. Constants declared as "extern..." are
// not correctly used in EVT_MENU() macros etc..
enum
{
   ID_REFOBJ_MENU_OBJ=wxID_HIGHEST+1,
   ID_REFOBJ_MENU_OBJ_SAVE,
   ID_REFOBJ_MENU_OBJ_LOAD,
   ID_REFOBJ_MENU_PAR,
   ID_REFOBJ_MENU_PAR_FIXALL,
   ID_REFOBJ_MENU_PAR_UNFIXALL,
   ID_REFOBJ_MENU_PAR_RANDOMIZE,
   ID_CRYST_UPDATEUI,
   ID_WXOBJ_COLLAPSE,
   ID_WXOBJ_NAME,
   ID_WXFIELD
};
#endif
// Forward declaration
class WXCrystObjBasicList;

/// Abstract base class for all objects in wxCryst
class WXCrystObjBasic: public wxWindow
{
   public:
      /// Constructor
      WXCrystObjBasic(wxWindow* parent);
      /// Destructor
      virtual ~WXCrystObjBasic();
      /** Get new values to be displayed from the underlying object,
      * and raise flag if an UI update is necessary.
      * The actual GUI update is not made here. UpdateUI() should be
      * called separately, from the main thread.
      *
      * \param updateUI: if true, this will call UpdateUI, either directly
      *(if in the main thread), or by sending a message.
      * \param mutexlock: if true, a Mutex will be used to lock the data shared
      * between main and background thread. The idea is to only use a few Mutexes
      * to lock data from the top objects (wxRefinableObj,...), when
      * calling CrystUpdate() and UpdateUI(). As sub-objects (WXField,...)
      * are only updated from within a top object, the mutex lock in the top object
      * will also lock the data in the sub-objects.
      */
      virtual void CrystUpdate(const bool updateUI=false,const bool mutexlock=false)=0;
      /** Update the User Interface, if necessary
      *
      * \param mutexlock: if true, a Mutex will be used to lock the data shared
      * between main and background thread.
      *
      * The idea is to only use a few Mutexes
      * to lock data from the top objects (wxRefinableObj,...), when
      * calling CrystUpdate() and UpdateUI(). As sub-objects (WXField,...)
      * are only updated from within a top object, the mutex lock in the top object
      * will also lock the data in the sub-objects.
      */
      virtual void UpdateUI(const bool mutexlock=false)=0;
      /** Ask for a new Layout with recalculated size hints, because
      * a child has been changed or added
      *
      * \param pChild: the modified child. If null, this will mean we are
      * asking for a new layout of the window itself (useful when a sub-window
      * is deleted).
      */
      virtual void BottomLayout(WXCrystObjBasic *pChild);
      /// Notify that a new children has been added, also adding it to the
      /// correct sizer (which can be the top sizer or not).
      ///
      /// \param doBottomLayout ask for a new Layout of the window and of its parents.
      virtual void AddChild(WXCrystObjBasic *pChild, bool doBottomLayout=true);
      /// \internal Tell the object it has been added to a list
      void AddedToList(WXCrystObjBasicList* list);
      /// \internal Tell the object it has been removed from a list
      void RemovedFromList(WXCrystObjBasicList* list);
      virtual bool Layout();
      /// Set tooltip for this window. The WXCrystObjBasicList::SetToolTip() will just
      /// call wxWindow::SetToolTip(), but is virtual and will correctly call specialized
      /// derived function. This is needed to avoid settin a tooltip on widgets
      /// not receiving motion events (e.g. static text under windows).
      virtual void SetToolTip(const wxString& tip);
   protected:
      /// Parent, if a WXCrystObjBasic itself
      WXCrystObjBasic *mWXCrystParent;
      /// Is the the window currently shown ?
      bool mIsShown;
      /// Do we need to update the display ?
      bool mNeedUpdateUI;
      /// Mutex used to lock data when preparing to update the UI in non-main thread
      CrystMutex mMutex;
      /// WXCrystObjBasicList which are aware of this object, 
      /// and which should be told on destruction.
      std::set<WXCrystObjBasicList*> mvpList;
};

/// A List of WXCrystObjBasic.
class WXCrystObjBasicList
{
   public:
      /// Constructor
      WXCrystObjBasicList();
      /// Destructor
      ~WXCrystObjBasicList();
      /// Number of objects.
      unsigned int GetNb()const;
      /// Add an object to the list. The object is just referenced,
      /// not copied.
      void Add(WXCrystObjBasic *);
      /// remove an object from the list
      void Remove(WXCrystObjBasic *);
      /// Show or hide all of the windows
      bool Show(bool);
      /** Forces all objects in the list to update. See WXCrystObjBasic::CrystUpdate()
      *
      * see WXCrystObjBasic::CrystUpdate() on the use of updateUI and mutexlock
      *
      * Normally WXCrystObjBasicList::CrystUpdate() should never be used with mutexlock=true,
      * as the mutex locking should be done from the calling object.
      */
      void CrystUpdate(const bool updateUI=false,const bool mutexlock=false);
      /** Forces all objects in the list to update. See WXCrystObjBasic::CrystUpdate()
      *
      * see WXCrystObjBasic::UpdateUI() on the use of updateUI and mutexlock.
      *
      * Normally WXCrystObjBasicList::UpdateUI() should never be used with mutexlock=true,
      * as the mutex locking should be done from the calling object.
      */
      void UpdateUI(const bool mutexlock=false);
      void Enable(bool enable);
   private:
      /// List of pointers to the objects.
      std::set<WXCrystObjBasic*> mvpWXCrystObj;
};

class WXFieldName;

/** Base class for all displayed ObjCryst objects (with a title,
* and a sizer to stack objects).
*
* A button (which should be used to collapse the object) is used
* to create an indentation for the sub-objects.
* 
* \todo Allow the objects to be collabsable. The difficulty is that
* even if the object is not shown, it is not removed by the Sizer as long
* as it is not deleted... Needs some testing ! Otherwise it would also
* be possible to delete and re-create sub-objects when collapsing, but
* that would be more difficult.
*/
class WXCrystObj: public WXCrystObjBasic
{
   public:
      /// Constructor, with a 
      WXCrystObj(wxWindow* parent,int orient=wxHORIZONTAL,bool showName=true);
      virtual ~WXCrystObj();
      /// Only display the title, and collapse everything else.
      /// \bug : the windows do collapse, but the size of the window is not
      /// changed, so it is pretty useless so far...
      void OnToggleCollapse(wxCommandEvent & WXUNUSED(event));
      /// When a WXFieldName has been changed by the user, it is handled
      /// here. This returns true if the value has been handled (for inheritance
      /// purposes).
      virtual bool OnChangeName(const int id)=0;
      virtual void CrystUpdate(const bool updateUI=false,const bool mutexlock=false);
      virtual void UpdateUI(const bool mutexlock=false);
      virtual void OnEnable(wxUpdateUIEvent &event);
      virtual bool Enable(bool enable);
      virtual void BottomLayout(WXCrystObjBasic *pChild);
      virtual void AddChild(WXCrystObjBasic *pChild, bool doBottomLayout=true);
   protected:
      /// Top sizer including the title and WXCrystObj::mpSizer
      wxBoxSizer *mpTopSizer;
      /// Sizer including all sub-objects
      wxBoxSizer *mpSizer;
      /// The title
      WXFieldName *mpWXTitle;
      /// To be used for collapsing the sub-objects.
      bool mIsExpanded;
      ///All windows but the title and collapse button are in this list
      WXCrystObjBasicList mList;
      /// The collapse button
      wxButton * mpCollapseButton;
   DECLARE_EVENT_TABLE()
};

/** This is the abstract base class for all fields, wether they contain
* a floating-point parameter, or a string,... All WXField have a title
* and an entry field.
*
* Note that WXField::CrystUpdate() and WXField::UpdateUI() should be done 
* from the parent object. Notably using WXField::CrystUpdate(updateui=true)
* will not trigger an update of the UI.
*/
class WXField: public WXCrystObjBasic
{
   public:
      /** Constructor, specifying the label of the field.
      *
      */
      WXField(wxWindow *parent,const string& label,const int field_id);
      /// Change the field's label.
      void SetLabel(const string&);
      /// After a user entry, this allows to go back to the last value, if for some reason
      /// the entry was rejected (because the object is currently busy, ...)
      virtual void Revert()=0;
      /// Change the colour of the field's title. Can be used (with parcimony)
      /// to clarify the interface.
      virtual bool SetForegroundColour(const wxColour& colour);
      /// This function shall be called when a new value has been entered.
      virtual void ValidateUserInput()=0;
      /// Change the size of the field (excluding the title)
      virtual void SetSize(int width, int height);
   protected:
      /// The horizontal sizer in which the title, button, fields, are put. 
      wxBoxSizer *mpSizer;
      /// The label
      wxStaticText *mpLabel;
      /// The Id of this field
      const int mId;
};

/// This function validates all user input (in a WXField) 
/// not yet taken into account, if needs be. This should be called by \b ALL
/// functions using data stored in fields (basically all functions !)
void WXCrystValidateAllUserInput();

/** A field which directly links to a string.
*
* 
*/
class WXFieldString:public WXField
{
   public:
      WXFieldString(wxWindow *parent,string& st,const int field_id,
                  const int hsize=50, bool isEditable=true);
      /// When a new value is entered (must type it and then hit the 'enter' key).
      /// The Field reads the new value, then
      /// forwards the event to its owner, who will take care of anything
      /// that must be done.
      void OnEnter(wxCommandEvent & event);
      /// Records when text is entered (either from self-updating or user input)
      void OnText(wxCommandEvent & WXUNUSED(event));
      /// This actually posts an UpdateUI event, so that it is safe to call it
      /// from a non-graphic thread.
      void SetValue(const string&);
      /// Get the current name.
      const string GetValue() const;
      virtual void CrystUpdate(const bool updateUI=false,const bool mutexlock=false);
      virtual void UpdateUI(const bool mutexlock=false);
      void Revert();
      virtual void ValidateUserInput();
      virtual void SetSize(int width, int height);
      /// Set tooltip for this window. It will be activated when going over the entry field.
      virtual void SetToolTip(const wxString& tip);
   protected:
      /// The WXCrystObj whose name is shown here
      string* mpString;
      /// Last name displayed.
      string mValue;
      /// The text window
      wxTextCtrl *mpField;
      /// Last name displayed, before the value was changed by the user. Not used yet,
      /// could be useful for undo.
      string mValueOld;
      /// Set to true if the Field is being updated, so that no 
      /// 'EVT_TEXT' is understood as user input.
      bool mIsSelfUpdating;
   DECLARE_EVENT_TABLE()
};
/** A field with the name of a WXCrystObj. Updating must be done by the WXCrystObj owner.
* For a simple string field linked directly to a string, use ObjCryst::WXFieldString
* 
*/
class WXFieldName:public WXField
{
   public:
      WXFieldName(wxWindow *parent,const string& label, WXCrystObj* owner,const int field_id,
                  const int hsize=50, bool isEditable=true);
      /// When a new value is entered (must type it and then hit the 'enter' key).
      /// The Field reads the new value, then
      /// forwards the event to its owner, who will take care of anything
      /// that must be done.
      void OnEnter(wxCommandEvent & event);
      /// Records when text is entered (either from self-updating or user input)
      void OnText(wxCommandEvent & WXUNUSED(event));
      /// This actually posts an UpdateUI event, so that it is safe to call it
      /// from a non-graphic thread.
      void SetValue(const string&);
      /// Get the current name.
      const string GetValue() const;
      /// This does nothing. Updates should be done by the owner in the particular
      /// case of names.
      virtual void CrystUpdate(const bool updateUI=false,const bool mutexlock=false);
      virtual void UpdateUI(const bool mutexlock=false);
      void Revert();
      virtual void ValidateUserInput();
      virtual void SetSize(int width, int height);
      /// Set tooltip for this window. It will be activated when going over the entry field.
      virtual void SetToolTip(const wxString& tip);
   protected:
      /// The WXCrystObj whose name is shown here
      WXCrystObj* mpWXObj;
      /// Last name displayed.
      string mValue;
      /// The text window
      wxTextCtrl *mpField;
      /// Last name displayed, before the value was changed by the user. Not used yet,
      /// could be useful for undo.
      string mValueOld;
      /// Set to true if the Field is being updated, so that no 
      /// 'EVT_TEXT' is understood as user input.
      bool mIsSelfUpdating;
   DECLARE_EVENT_TABLE()
};

/// A field for a parameter. This is a an abstract bas class, which can
/// handle events (the real classes to use is the templated WXFieldPar class).
/// If the parameter is a RefinablePar, use WXFieldRefPar.
class WXFieldParBase:public WXField
{
   public:
      /// Constructor
      WXFieldParBase(wxWindow *parent,const string& label, const int field_id,
                     const int hsize=65);
      /// When a new value is entered (must type it and then hit the 'enter' key).
      /// The Field reads the new value, 
      /// and directly changes the RefinablePar value (contrary to what happens
      /// for WXFieldName)by using RefinablePar::SetHumanValue().
      void OnEnter(wxCommandEvent & WXUNUSED(event));
      /// Records when text is entered (either from self-updating or user input)
      void OnText(wxCommandEvent & WXUNUSED(event));
      /// This gets a new value from the parameter.
      virtual void CrystUpdate(const bool updateUI=false,const bool mutexlock=false)=0;
      virtual void Revert()=0;
      virtual void ValidateUserInput();
      /// Set tooltip for this window. It will be activated when going over the entry field.
      virtual void SetToolTip(const wxString& tip);
      /// Set Format
      void SetFormat(const wxString &format);
   protected:
      /// Reads the new value when the Enter key is hit
      virtual void ReadNewValue()=0;
      /// The field in which the value is written.
      wxTextCtrl *mpField;
      /// Set to true if the Field is being updated, so that no 
      /// 'EVT_TEXT' is understood as user input.
      bool mIsSelfUpdating;
      /// Format to be used, default = _T("%8f")
      wxString mFormat;
   DECLARE_EVENT_TABLE()
};

/// A field for a parameter. Template version.
/// If the parameter is a RefinablePar, use WXFieldRefPar instead.
template<class T>class WXFieldPar:public WXFieldParBase
{
   public:
      /// Constructor
      WXFieldPar(wxWindow *parent,const string& label, const int field_id,
                    T *par,const int hsize=65);
      /// This gets a new value from the parameter.
      virtual void CrystUpdate(const bool updateUI=false,const bool mutexlock=false);
      virtual void UpdateUI(const bool mutexlock=false);
      virtual void Revert();
      /// Set Coefficient between the value used by ObjCryst++ and the one
      /// to be displayed to the user. Typically, 180/pi
      void SetHumanValueScale(const T s);
   protected:
      /// Reads the new value when the Enter key is hit
      virtual void ReadNewValue();
      /// A pointer to the value displayed
      T* mpValue;
      /// The value displayed
      T mValue;
      /// Last value
      T mValueOld;
      /// Coefficient between the value used by ObjCryst++ and the one
      /// to be displayed to the user. Typically, 180/pi
      T mHumanScale;
};

/// Class to pick one choice... Choice change/update is handled by the WXCrystObj owner, who
/// should grab the incoming event.
/// Useful, for example, to change the scattering power associated to an atom.
class WXFieldChoice:public WXField
{
   public:
      /// Constructor
      WXFieldChoice(wxWindow *parent,const int field_id,
                            const string &name,const int hsize=80);
      /// Does nothing
      virtual void CrystUpdate(const bool updateUI=false,const bool mutexlock=false);
      /// Does nothing
      virtual void UpdateUI(const bool mutexlock=false);
      void Revert();
      /// Used by the owner to change the name of the choice
      void SetValue(const string&);
      /// Unnecessary here. Any change is immediately taken into account.
      virtual void ValidateUserInput();
   protected:
      /// The button to be clicked to change the value.
      wxButton *mpButton;
};

/// Our own local menu bar, using buttons and Popup menus
class WXCrystMenuBar: public WXCrystObjBasic
{
   public:
      /// Ctor
      WXCrystMenuBar(wxWindow *parent, WXCrystObj* owner);
      /// Add a menu
      void AddMenu(const string &name,const int menuId, const string& help="");
      /// Get access to a menu
      wxMenu& GetMenu(const int menuId);
      /// Add an entry to a menu
      void AddMenuItem(const int menuId, int id, const string& item, const string& help="",
                       const bool checkable= false);
      /// Add a sub-menu to a menu
      void AddMenuItem(const int menuId,int id, const wxString&  item,
                       wxMenu *subMenu, const wxString& helpString = _T(""));
      virtual void CrystUpdate(const bool updateUI=false,const bool mutexlock=false);
      virtual void UpdateUI(const bool mutexlock=false);
      /// Event handler to popu the menu when the button is clicked.
      void OnPopupMenu(wxCommandEvent & event);
      /// Set tooltip for each menu.
      virtual void SetToolTip(const wxString& tip,long menu=0);
   protected:
      /// The sizer of the menu
      wxBoxSizer* mpSizer;
      /// List of menus, first is the menu Id and second is a pair of
      /// <pointer to the menu, pointer to the button of the menu>
      std::map<long,pair<wxMenu *,wxButton*> > mvpMenu;
   DECLARE_EVENT_TABLE()
};

} //namespace

#endif //_VFN_WX_CRYST_H_
