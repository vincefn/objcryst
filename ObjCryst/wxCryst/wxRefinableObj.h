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
} //namespace

#include "wxCryst/wxCryst.h"

// We should not have to put this here. !! :TODO:
#include "RefinableObj/RefinableObj.h"

namespace ObjCryst
{

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
                    RefinablePar *refpar,const int hsize=50);
      /// UpdateUI does not grab new values in the underlying object,
      /// but only updates the values which have been supplied.
      void OnUpdateUI(wxUpdateUIEvent & WXUNUSED(event));
      /// When a new value is entered (must type it and then hit the 'enter' key).
      /// The Field reads the new value, 
      /// and directly changes the RefinablePar value (contrary to what happens
      /// for WXFieldName)by using RefinablePar::SetHumanValue().
      void OnEnter(wxCommandEvent & WXUNUSED(event));
      /// Toggle the 'fixed' status of the parameter.
      void OnToggleFix(wxCommandEvent & WXUNUSED(event));
      /// Opens the popu menu, to allow changing limits
      void OnPopupMenu(wxCommandEvent & event);
      /// Opens the popu menu, to allow changing limits
      void OnPopupMenuChoice(wxMenuEvent& event);
      /// This gets a new value from the RefinablePar, and then posts an
      /// OnUpdateUI event. (never immediately update the GUI)
      void CrystUpdate();
      /// Get the RefinablePar associated to this field
      RefinablePar& GetRefPar();
      void Revert();
   protected:
      double mValue;
      wxCheckBox *mpButtonFix;
      wxTextCtrl *mpField;
      RefinablePar *mpRefPar;
      double mValueOld;
      wxMenu *mpPopUpMenu;
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
      void Revert();
   protected:
      int mChoice;
      int mChoiceOld;
      RefObjOpt* mpOption;
      wxChoice *mpList;
   DECLARE_EVENT_TABLE()
};

/// Field for a RefinableObj cost function
class WXCostFunction:public WXField
{
   public:
      WXCostFunction(wxWindow *parent,RefinableObj *obj, const int field_id,
                     const int funcNum,double * weight);
      void OnUpdateUI(wxUpdateUIEvent & WXUNUSED(event));
      void OnEnter(wxCommandEvent & WXUNUSED(event));
      virtual void CrystUpdate();
      virtual void Revert();
   protected:
      wxTextCtrl *mpValue;
      double mValue;
      RefinableObj *mpObj;
      const int mFuncNum;
      WXFieldPar<double> *mpWeight;
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
/// a menu for Input/Output, and all RefObjOpt.
class WXRefinableObj: public WXCrystObj
{
   public:
      WXRefinableObj(wxWindow *parent, RefinableObj*);
      ~WXRefinableObj();
      bool Layout();
      virtual void CrystUpdate();
      virtual bool OnChangeName(const int id);
      void OnMenuSave(wxCommandEvent & WXUNUSED(event));
      void OnMenuLoad(wxCommandEvent & WXUNUSED(event));
      void OnMenuFixAllPar(wxCommandEvent & WXUNUSED(event));
      void OnMenuUnFixAllPar(wxCommandEvent & WXUNUSED(event));
      void OnMenuParRandomize(wxCommandEvent & WXUNUSED(event));
   protected:
      WXCrystMenuBar* mpMenuBar;
   private:
      RefinableObj* mpRefinableObj;
   DECLARE_EVENT_TABLE()
};
} //namespace

#endif //_VFN_WX_REFINABLE_OBJ_H_
