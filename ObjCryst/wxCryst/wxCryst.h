/*
* LibCryst++ : a Crystallographic computing library in C++
*
*  (c) 2000 Vincent FAVRE-NICOLIN
*           Laboratoire de Cristallographie
*           24, quai Ernest-Ansermet, CH-1211 Geneva 4, Switzerland
*  Contact: Vincent.Favre-Nicolin@cryst.unige.ch
*           Radovan.Cerny@cryst.unige.ch
*
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

#include "wx/window.h"
#include "wx/button.h"
#include "wx/sizer.h"
#include "wx/stattext.h"
#include "wx/textctrl.h"
#include "wx/checkbox.h"
#include "wx/choice.h"
#include "CrystVector/CrystVector.h"

namespace ObjCryst
{
#ifndef DOXYGEN_SKIP_THIS
enum
{
   ID_FIX_SCROLLBARS,//for FOX
   ID_WXOBJ_COLLAPSE,
   ID_WXOBJ_NAME,
   ID_WXFIELD,
   ID_WXFIELD_REFPAR,//id of all refpar window
   ID_WXFIELD_REFPAR_FIXBUTTON,
   ID_REFPAR_POPUP_SET_LIMITS,
   ID_REFPAR_POPUP_REMOVE_LIMITS,
   ID_WXFIELD_OPTION,//Actually not the id of a window
   ID_WXFIELD_COSTFUNC,//Actually not the id of a window
   ID_CRYST_MENU1,
   ID_CRYST_MENU2,   // This is done to 'reserve' some id for wxCrystMenu
   ID_CRYST_MENU3,
   ID_CRYST_MENU4,
   ID_CRYST_MENU5,
   ID_CRYST_MENU6,
   ID_CRYST_MENU7,
   ID_CRYST_MENU8,
   ID_CRYST_MENU9,
   ID_CRYST_MENU10,
   ID_CRYST_MENU11,
   ID_CRYST_MENU12,
   ID_CRYST_MENU13,
   ID_CRYST_MENU14,
   ID_CRYST_MENU15,
   ID_CRYST_MENU16,
   ID_REFOBJ_MENU_OBJ,
   ID_REFOBJ_MENU_OBJ_SAVE,
   ID_REFOBJ_MENU_OBJ_LOAD,
   ID_REFOBJ_MENU_PAR,
   ID_REFOBJ_MENU_PAR_FIXALL,
   ID_REFOBJ_MENU_PAR_UNFIXALL,
   ID_REFOBJ_MENU_PAR_RANDOMIZE,
   ID_WXSCATTPOWATOM_SYMBOL,
   ID_SCATTPOWATOM_MENU_COLOUR,
   ID_SCATTPOWATOM_MENU_COLOUR_SETRGB,
   ID_ATOM_SCATTPOW,
   ID_CRYSTAL_MENU_SAVECIF,
   ID_CRYSTAL_MENU_SAVETEXT,
   ID_CRYSTAL_MENU_DISPLAY,
   ID_CRYSTAL_MENU_DISPLAY_3DVIEW,
   ID_CRYSTAL_MENU_SCATT,
   ID_CRYSTAL_MENU_PAR_ADDANTIBUMP,
   ID_CRYSTAL_MENU_SCATT_REMOVESCATTPOW,
   ID_CRYSTAL_MENU_SCATT_ADDSCATTPOWATOM,
   ID_CRYSTAL_MENU_SCATT_REMOVESCATTERER,
   ID_CRYSTAL_MENU_SCATT_ADDATOM,
   ID_CRYSTAL_MENU_SCATT_ADDZSCATTERER,
   ID_CRYSTAL_MENU_SCATT_ADDTETRAHEDRON,
   ID_CRYSTAL_MENU_SCATT_ADDOCTAHEDRON,
   ID_CRYSTAL_MENU_SCATT_ADDTRIANGLE,
   ID_CRYSTAL_MENU_SCATT_ADDSQUAREPLANE,
   ID_CRYSTAL_MENU_SCATT_ADDCUBE,
   ID_CRYSTAL_MENU_SCATT_ADDANTIPRISMTETRAGONAL,
   ID_CRYSTAL_MENU_SCATT_ADDPRISMTRIGONAL,
   ID_CRYSTAL_SPACEGROUP,
   ID_GLCRYSTAL_MENU_UPDATE,
	ID_GLCRYSTAL_MENU_CHANGELIMITS,
   ID_ZATOM_NAME,
   ID_ZATOM_SCATTPOW,
   ID_ZATOM_BOND,
   ID_ZATOM_ANGLE,
   ID_ZATOM_DIHED,
   ID_ZSCATTERER_MENU_ATOM,
   ID_ZSCATTERER_MENU_ATOM_ADD,
   ID_ZSCATTERER_MENU_ATOM_CHANGE_PIVOT,
   ID_ZSCATTERER_MENU_PAR_LIMITS_RELAT_BOND,
   ID_ZSCATTERER_MENU_PAR_LIMITS_RELAT_ANGLE,
   ID_ZSCATTERER_MENU_PAR_LIMITS_RELAT_DIHED,
   ID_POWDERSPECTRUM_MENU_SCATT_ADDCOMPBACKGD,
   ID_POWDERSPECTRUM_MENU_SCATT_ADDCOMPCRYST,
   ID_POWDERSPECTRUM_MENU_GRAPH,
   ID_POWDERSPECTRUM_MENU_SAVETEXT,
   ID_POWDERSPECTRUM_MENU_SIMULATE,
   ID_POWDERSPECTRUM_MENU_IMPORT_FULLPROF,
   ID_POWDERSPECTRUM_MENU_IMPORT_PSI_DMC,
   ID_POWDERSPECTRUM_MENU_IMPORT_ILL_D1A5,
   ID_POWDERSPECTRUM_MENU_IMPORT_XDD,
   ID_POWDERSPECTRUM_MENU_IMPORT_CPI,
   ID_POWDERSPECTRUM_MENU_IMPORT_2THETAOBSSIGMA,
   ID_POWDERSPECTRUM_MENU_IMPORT_2THETAOBS,
   ID_POWDERSPECTRUM_MENU_FITSCALE_R,
   ID_POWDERSPECTRUM_MENU_FITSCALE_RW,
   ID_POWDERSPECTRUM_MENU_WAVELENGTH,
   ID_POWDERSPECTRUM_MENU_WAVELENGTH_XRAY,
   ID_POWDERSPECTRUM_MENU_WAVELENGTH_NEUTRON,
   ID_POWDERSPECTRUM_MENU_WAVELENGTH_SET,
   ID_POWDERSPECTRUM_MENU_WAVELENGTH_SET_AG,
   ID_POWDERSPECTRUM_MENU_WAVELENGTH_SET_MO,
   ID_POWDERSPECTRUM_MENU_WAVELENGTH_SET_CU,
   ID_POWDERSPECTRUM_MENU_WAVELENGTH_SET_FE,
   ID_POWDERSPECTRUM_MENU_WAVELENGTH_SET_CR,
   ID_POWDERSPECTRUM_MENU_WAVELENGTH_SET_AGA1,
   ID_POWDERSPECTRUM_MENU_WAVELENGTH_SET_MOA1,
   ID_POWDERSPECTRUM_MENU_WAVELENGTH_SET_CUA1,
   ID_POWDERSPECTRUM_MENU_WAVELENGTH_SET_FEA1,
   ID_POWDERSPECTRUM_MENU_WAVELENGTH_SET_CRA1,
   ID_POWDERSPECTRUM_MENU_ADD_2THETA_EXCLUDE,
   ID_POWDERSPECTRUMGRAPH_MENU_UPDATE,
   ID_POWDERSPECTRUMBACKGROUND_IMPORT,
   ID_POWDERSPECTRUMDIFF_CRYSTAL,
   ID_POWDERSPECTRUM_GRAPH_NEW_PATTERN,   
   ID_GLOBALOPT_MENU_GLOBAlOPT,
   ID_GLOBALOPT_MENU_GLOBAlOPT_ADDOBJ,
   ID_GLOBALOPT_MENU_GLOBAlOPT_ADDCOSTFUNC,
   ID_GLOBALOPT_MENU_GLOBAlOPT_RUN,
   ID_GLOBALOPT_MENU_GLOBAlOPT_STOP,
	ID_DIFFSINGLECRYST_MENU_SAVEHKLIOBSICALC,   
	ID_DIFFSINGLECRYST_MENU_SIMULATE,           
	ID_DIFFSINGLECRYST_MENU_IMPORT_HKLIOBS,     
	ID_DIFFSINGLECRYST_MENU_IMPORT_HKLIOBSSIGMA,
	ID_DIFFSINGLECRYST_MENU_IMPORT_JANAM91,     
	ID_DIFFSINGLECRYST_MENU_FITSCALE_R,         
	ID_DIFFSINGLECRYST_MENU_FITSCALE_RW,
   ID_DIFFSINGLECRYST_MENU_WAVELENGTH,
   ID_DIFFSINGLECRYST_MENU_WAVELENGTH_XRAY,
   ID_DIFFSINGLECRYST_MENU_WAVELENGTH_NEUTRON,
   ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET,
   ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_AG,
   ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_MO,
   ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_CU,
   ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_FE,
   ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_CR,
   ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_AGA1,
   ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_MOA1,
   ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_CUA1,
   ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_FEA1,
   ID_DIFFSINGLECRYST_MENU_WAVELENGTH_SET_CRA1,
	ID_DIFFSINGLECRYST_CRYSTAL
};
#endif
/// Abstract base class for all objects in wxCryst
class WXCrystObjBasic: public wxWindow
{
   public:
		/// Constructor
      WXCrystObjBasic(wxWindow* parent);
		/// Destructor
      virtual ~WXCrystObjBasic();
      /// Update the display, by getting new values from the object. 
      /// New values should be grabbed from the object, and then a
      /// wxUpdateUI event should be posted (this for multi-thread, since
		/// only one thread should do graphical changes)
      virtual void CrystUpdate()=0;
   protected:
		/// Parent 
      wxWindow *mWXParent;
		/// Is the the window currently shown ?
      bool mIsShown;
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
      void Remove(const WXCrystObjBasic *);
      /// Get an object from the list
      WXCrystObjBasic* Get(const unsigned int i);
      /// Show or hide all of the windows
      bool Show(bool);
      /// Forces all objects in the list to update. See WXCrystObjBasic::CrystUpdate()
      void CrystUpdate();
   private:
      /// Number of objects.
      unsigned int mNbWXCrystObj;
		/// \internal Maximum number of objects
      unsigned int mMaxNbWXCrystObj;
		/// Array of pointers to the objects.
      WXCrystObjBasic** mpWXCrystObj;
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
      /// Fix the Layout of the window, resize if necessary.
      bool Layout();
      /// Only display the title, and collapse everything else.
      /// \bug : the windows do collapse, but the size of the window is not
      /// changed, so it is pretty useless so far...
      void OnToggleCollapse(wxCommandEvent & WXUNUSED(event));
      /// When a WXFieldName has been changed by the user, it is handled
      /// here. This returns true if the value has been handled (for inheritance
      /// purposes).
      virtual bool OnChangeName(const int id)=0;
      virtual void CrystUpdate();
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
* \todo Currently for all fields, the user must hit the 'enter' key to validate
* the input, which is annoying. So it would be nice to use another way. For 
* example by recording when some typing has been made in a field, and 
* doing the actual reading afterwards (eg a global function should be 
* created to read the input for all fields, and that could be 
* complicated to handle).
*/
class WXField: public WXCrystObjBasic
{
   public:
		/** Constructor, specifying the label of the field.
		*
		*/
      WXField(wxWindow *parent,const string& label,const int field_id);
		/// Redo the layout of the field.
      bool Layout();
		/// Change the field's label.
      void SetLabel(const string&);
      /// After a user entry, this allows to go back to the last value, if for some reason
      /// the entry was rejected (because the object is currently busy, ...)
      virtual void Revert()=0;
		/// Change the colour of the field's title. Can be used (with parcimony)
		/// to clarify the interface.
      virtual bool SetForegroundColour(const wxColour& colour);
   protected:
		/// The horizontal sizer in which the title, button, fields, are put. 
      wxBoxSizer *mpSizer;
		/// The label
      wxStaticText *mpLabel;
		/// The Id of this field
      const int mId;
};

/** A field with the name of a WXCrystObj. Updating must be done by the WXCrystObj owner.
*
* 
*/
class WXFieldName:public WXField
{
   public:
      WXFieldName(wxWindow *parent,const string& label, WXCrystObj* owner,const int field_id,
                  const int hsize=50, bool isEditable=true);
      /// UpdateUI does not grab new values in the underlying object,
      /// but only updates the values which have been supplied.
      void OnUpdateUI(wxUpdateUIEvent & WXUNUSED(event));
      /// When a new value is entered (must type it and then hit the 'enter' key).
      /// The Field reads the new value, then
      /// forwards the event to its owner, who will take care of anything
      /// that must be done.
      void OnEnter(wxCommandEvent & event);
      /// This posts an UpdateUI event.
      void SetValue(const string&);
		/// Get the current name.
      const string GetValue() const;
      /// This does nothing. Updates should be done by the owner in the particular
      /// case of names.
      virtual void CrystUpdate();
      void Revert();
   protected:
		/// The WXCrystObj whose name is shown here
      WXCrystObj* mpWXObj;
		/// Last name displayed.
      wxString mValue;
		/// The text window
      wxTextCtrl *mpField;
		/// Last name displayed, before the value was changed by the user. Not used yet,
		/// could be useful for undo.
      wxString mValueOld;
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
                     const int hsize=50);
      /// UpdateUI does not grab new values in the underlying object,
      /// but only updates the values which have been supplied.
      void OnUpdateUI(wxUpdateUIEvent & WXUNUSED(event));
      /// When a new value is entered (must type it and then hit the 'enter' key).
      /// The Field reads the new value, 
      /// and directly changes the RefinablePar value (contrary to what happens
      /// for WXFieldName)by using RefinablePar::SetHumanValue().
      void OnEnter(wxCommandEvent & WXUNUSED(event));
      /// This gets a new value from the parameter.
      virtual void CrystUpdate()=0;
      /// Revert to the previous value
      virtual void Revert()=0;
   protected:
      /// Reads the new value when the Enter key is hit
      virtual void ReadNewValue()=0;
		/// Apply a new value
      virtual void ApplyNewValue()=0;
		/// The field in which the value is written.
      wxTextCtrl *mpField;
   DECLARE_EVENT_TABLE()
};

/// A field for a parameter. Template version.
/// If the parameter is a RefinablePar, use WXFieldRefPar instead.
template<class T>class WXFieldPar:public WXFieldParBase
{
   public:
		/// Constructor
      WXFieldPar(wxWindow *parent,const string& label, const int field_id,
                    T *par,const int hsize=50);
      /// This gets a new value from the parameter.
      virtual void CrystUpdate();
      /// Revert to the previous value
      virtual void Revert();
   protected:
      /// Reads the new value when the Enter key is hit
      virtual void ReadNewValue();
      /// Applies a new value and shows it.
      virtual void ApplyNewValue();
		/// A pointer to the value displayed
      T* mpValue;
		/// The last value displayed, before it was changed by the user.
      T mValueOld;
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
      bool Layout();
      /// UpdateUI does not grab new values in the underlying object,
      /// but only updates the values which have been supplied.
      void OnUpdateUI(wxUpdateUIEvent & WXUNUSED(event));
      /// Does nothing
      virtual void CrystUpdate();
      /// Does nothing
      void Revert();
      /// Used by the owner to change the name of the choice
      void SetValue(const string&);
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
		/// Redo the Layout
      bool Layout();
		/// Add a menu
      void AddMenu(const string &name,const int menuId, const string& help="");
		/// Add an entry to a menu
      void AddMenuItem(const int menuId, int id, const string& item, const string& help="",
                       const bool checkable= false);
		/// Add a sub-menu to a menu
      void AddMenuItem(const int menuId,int id, const wxString&  item,
                       wxMenu *subMenu, const wxString& helpString = "");
      virtual void CrystUpdate();
		/// Event handler to popu the menu when the button is clicked.
      void OnPopupMenu(wxCommandEvent & event);
   protected:
		/// The sizer of the menu
      wxBoxSizer* mpSizer;
		/// The list of menu IDs
      CrystVector_int mMenuId;
		/// Number of menus
      unsigned int mNbMenu;
		/// Max number of menus
      unsigned int mMaxNbMenu;
		/// Array of menus
      wxMenu **mpMenu;
		/// The buttons corresponding to each menu
      wxButton **mpButton;
   DECLARE_EVENT_TABLE()
};

} //namespace

#endif //_VFN_WX_CRYST_H_
