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

#include "ObjCryst/wxCryst/wxRefinableObj.h"

#include "ObjCryst/wxCryst/wxRefinableObj.h"
#include "ObjCryst/Quirks/VFNStreamFormat.h"

#include "ObjCryst/RefinableObj/GlobalOptimObj.h"

//These are only for explicit instantiation
#include "ObjCryst/ObjCryst/Atom.h"
#include "ObjCryst/ObjCryst/Crystal.h"
#include "ObjCryst/ObjCryst/DiffractionDataSingleCrystal.h"
#include "ObjCryst/ObjCryst/PowderPattern.h"
#include "ObjCryst/ObjCryst/ScatteringPower.h"
#include "ObjCryst/ObjCryst/ZScatterer.h"
#include "ObjCryst/ObjCryst/ScatteringCorr.h"
#include "ObjCryst/ObjCryst/ReflectionProfile.h"

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
/// This pointer records the last wxField in which something was enetered,
/// so that it can be validated when inpu is finished (either when another
/// input has begun in another field, or when an action requires to purge
/// all input
extern WXField *spLastWXFieldInputNotValidated;
////////////////////////////////////////////////////////////////////////
//
//    WXFieldRefPar
//
////////////////////////////////////////////////////////////////////////
static const long ID_WXFIELD_REFPAR              =WXCRYST_ID();
static const long ID_WXFIELD_REFPAR_FIXBUTTON    =WXCRYST_ID();
static const long ID_WXFIELD_REFPAR_LIMITEDBUTTON=WXCRYST_ID();
static const long ID_REFPAR_POPUP_SET_LIMITS     =WXCRYST_ID();

BEGIN_EVENT_TABLE(WXFieldRefPar,wxWindow)
   EVT_TEXT_ENTER(ID_WXFIELD,                   WXFieldRefPar::OnEnter)
   EVT_TEXT(ID_WXFIELD,                         WXFieldRefPar::OnText)
   EVT_CHECKBOX(ID_WXFIELD_REFPAR_FIXBUTTON,    WXFieldRefPar::OnToggleFix)
   EVT_CHECKBOX(ID_WXFIELD_REFPAR_LIMITEDBUTTON,WXFieldRefPar::OnToggleLimited)
   EVT_RIGHT_DOWN(                              WXFieldRefPar::OnPopupMenu)
   EVT_MENU(ID_REFPAR_POPUP_SET_LIMITS,         WXFieldRefPar::OnPopupMenuChoice)
END_EVENT_TABLE()

WXFieldRefPar::WXFieldRefPar(wxWindow *parent,const string& label,
                             RefinablePar *par, const int hsize,
                             const bool enableFixButton, const bool enableLimitedButton):
WXField(parent,label,ID_WXFIELD_REFPAR),mValue(0.),mpRefPar(par),mIsSelfUpdating(false),mFormat(_T("%8f"))
{
   VFN_DEBUG_MESSAGE("WXFieldRefPar::WXFieldName():End",6)
   if(enableFixButton)
   {
      this->SetLabel(label+"R");
      mpButtonFix=new wxCheckBox(this,ID_WXFIELD_REFPAR_FIXBUTTON,_T("L"),wxDefaultPosition, wxDefaultSize);
      mpButtonFix->Fit();
      mpButtonFix->SetToolTip(_T("Check this box to enable optimizing this parameter.\n") 
                              _T("(some parameters may be automatically fixed for global optimization)"));
      mpSizer->Add(mpButtonFix,0,wxALIGN_CENTER);
   }else mpButtonFix=0;
   if(enableLimitedButton)
   {
      mpButtonLimited=new wxCheckBox(this,ID_WXFIELD_REFPAR_LIMITEDBUTTON,_T(""),
                                     wxDefaultPosition, wxSize(16,20));
      mpButtonLimited->SetToolTip(_T("Check this box to use limits for this parameter"));
      mpSizer->Add(mpButtonLimited,0,wxALIGN_CENTER);
   }else mpButtonLimited=0;
   
   mpField=new wxTextCtrl(this,ID_WXFIELD,_T(""),
                            wxDefaultPosition,wxSize(hsize,-1),wxTE_PROCESS_ENTER,
                            wxTextValidator(wxFILTER_NUMERIC));
   mpSizer->Add(mpField,0,wxALIGN_CENTER);
   this->BottomLayout(0);
}
WXFieldRefPar::~WXFieldRefPar()
{
   mpRefPar->WXNotifyDelete();
}

void WXFieldRefPar::OnEnter(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXFieldRefPar::OnEnter()",6)
   WXCrystValidateAllUserInput();
}
void WXFieldRefPar::OnText(wxCommandEvent & WXUNUSED(event))
{   
   if(true==mIsSelfUpdating) return;
   VFN_DEBUG_MESSAGE("WXFieldRefPar::OnEnter()",6)
   if(spLastWXFieldInputNotValidated!=this)
   {
      WXCrystValidateAllUserInput();
      spLastWXFieldInputNotValidated=this;
   }
}

void WXFieldRefPar::OnToggleFix(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXFieldRefPar::OnToggleFix()",6)
   if(0!=mpButtonFix) mpRefPar->SetIsFixed(!(mpButtonFix->GetValue()));
   mpRefPar->Print();
}

void WXFieldRefPar::OnToggleLimited(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXFieldRefPar::OnToggleLimited()",6)
   if(0!=mpButtonLimited) mpRefPar->SetIsLimited(mpButtonLimited->GetValue());
   mpRefPar->Print();
}

void WXFieldRefPar::OnPopupMenu(wxMouseEvent & WXUNUSED(event))
{
   static wxMenu sWXFieldRefParPopupMenu;//("Refinable Parameter");
   static bool needInitMenu=true;
   if(needInitMenu)
   {
      needInitMenu=false;
      sWXFieldRefParPopupMenu.Append(ID_REFPAR_POPUP_SET_LIMITS, _T("Change Limits"));
   }
   this->PopupMenu(&sWXFieldRefParPopupMenu,0,0);
}

void WXFieldRefPar::OnPopupMenuChoice(wxCommandEvent& event)
{
   VFN_DEBUG_MESSAGE("WXFieldRefPar::OnPopupMenuChoice()",7)
   if(event.GetId()== ID_REFPAR_POPUP_SET_LIMITS)
   {
      double min,max;
      {
         wxString str;
         str << mpRefPar->GetHumanMin();
         wxTextEntryDialog limitDialog(this,_T("Enter the minimum value"),
                                 _T("Minimum"),str,wxOK | wxCANCEL);
         if(wxID_OK!=limitDialog.ShowModal())
         {
            VFN_DEBUG_EXIT("WXZScatterer::OnMenuSetLimits():Cancelled",6)
            return;
         }
         limitDialog.GetValue().ToDouble(&min);
      }
      {
         wxString str;
         str << mpRefPar->GetHumanMax();
         wxTextEntryDialog limitDialog(this,_T("Enter the maximum value"),
                                 _T("Maximum"),str,wxOK | wxCANCEL);
         if(wxID_OK!=limitDialog.ShowModal())
         {
            VFN_DEBUG_EXIT("WXZScatterer::OnMenuSetLimits():Cancelled",6)
            return;
         }
         limitDialog.GetValue().ToDouble(&max);
      }
      if(max<=min)
      {
         wxMessageDialog dumbUser(this,_T("max <= min !!!"),
                                  _T("Whooops"),wxOK|wxICON_EXCLAMATION);
         dumbUser.ShowModal();
         return;
      }
      mpRefPar->SetHumanMin(min);
      mpRefPar->SetHumanMax(max);
      mpRefPar->SetIsLimited(true);
      mpRefPar->Print();
   }
}

void WXFieldRefPar::CrystUpdate(const bool uui,const bool lock)
{
   VFN_DEBUG_MESSAGE("WXFieldRefPar::CrystUpdate()",6)
   if(lock) mMutex.Lock();
   bool needUpdate=false;
   if(wxThread::IsMain())
   {
      if(mpRefPar->IsUsed()!=this->IsShown()) needUpdate=true;
      if(0!=mpButtonFix) if(mpButtonFix->GetValue()==mpRefPar->IsFixed()) needUpdate=true;
      if(0!=mpButtonLimited) if(mpButtonLimited->GetValue()==mpRefPar->IsLimited()) needUpdate=true;
   }
   if(mValue!=mpRefPar->GetHumanValue()) needUpdate=true;
   if(!needUpdate)
   {
      if(lock) mMutex.Unlock();
      return;
   }
   mValueOld=mValue;
   mValue=mpRefPar->GetHumanValue();
   mNeedUpdateUI=true;
   if(lock) mMutex.Unlock();
   if(mNeedUpdateUI && wxThread::IsMain()) this->UpdateUI(lock);
}

void WXFieldRefPar::UpdateUI(const bool lock)
{
   VFN_DEBUG_MESSAGE("WXFieldRefPar::UpdateUI("<<lock<<")"<<"MainThread="<<wxThread::IsMain(),3)
   if(lock) mMutex.Lock();
   if(mNeedUpdateUI==false)
   {
      if(lock) mMutex.Unlock();
      return;
   }
   if(false==mpRefPar->IsUsed()) this->Show(false);
   else this->Show(true);
   
   if(mpField==0)
   {
      if(lock) mMutex.Unlock();
      return;
   }
   
   //mpField->SetValue(wxString::Printf("%f",mValue));
   wxString tmp;
   //if((abs(mValue)<100)&&(abs(mValue)>0.01)) tmp.Printf(_T("%6.4f"),mValue);
   tmp.Printf(mFormat,mValue);
   mIsSelfUpdating=true;
   mpField->SetValue(tmp);
   mIsSelfUpdating=false;
   if(0!=mpButtonFix) mpButtonFix->SetValue(!(mpRefPar->IsFixed()));
   if(0!=mpButtonLimited) mpButtonLimited->SetValue(mpRefPar->IsLimited());
   mNeedUpdateUI=false;
   if(lock) mMutex.Unlock();
}

void WXFieldRefPar::Revert()
{
   VFN_DEBUG_MESSAGE("WXFieldRefPar::Revert()",6)
   wxMutexLocker mlock(mMutex);
   mValue=mValueOld;
   mNeedUpdateUI=true;
}
void WXFieldRefPar::ValidateUserInput()
{
   VFN_DEBUG_MESSAGE("WXFieldRefPar::ValidateUserInput()",6)
   wxMutexLocker mlock(mMutex);
   mValueOld=mValue;
   wxString s=mpField->GetValue();
   double tmp;
   s.ToDouble(&tmp);
   mValue=tmp;
   mpRefPar->SetHumanValue(mValue);
}

void WXFieldRefPar::SetToolTip(const wxString& tip){mpField->SetToolTip(tip);}

void WXFieldRefPar::SetFormat(const wxString &format)
{
   mFormat=format;
}

////////////////////////////////////////////////////////////////////////
//
//    WXFieldOption
//
////////////////////////////////////////////////////////////////////////
BEGIN_EVENT_TABLE(WXFieldOption,wxWindow)
   EVT_CHOICE(ID_WXFIELD,WXFieldOption::OnChoice)
END_EVENT_TABLE()

WXFieldOption::WXFieldOption(wxWindow *parent,
                             const int field_id,RefObjOpt* option):
WXField(parent,option->GetName(),field_id),
mChoice(-1),mChoiceOld(-1),mpOption(option),mpList(0)
{
   wxString choices[20];//:TODO: dynamically choose correct number
   for(int i=0;i<mpOption->GetNbChoice();i++)
      choices[i]=wxString::FromAscii(mpOption->GetChoiceName(i).c_str());
   
   mpList= new wxChoice(this,ID_WXFIELD,wxDefaultPosition,wxDefaultSize,
                        mpOption->GetNbChoice(),choices);
   mpSizer->Add(mpList,0,wxALIGN_CENTER);
   this->BottomLayout(0);
}
WXFieldOption::~WXFieldOption()
{
   mpOption->WXNotifyDelete();
}
void WXFieldOption::OnChoice(wxCommandEvent & WXUNUSED(event))
{
   if(mChoice==mpList->GetSelection()) return;
   mChoiceOld=mChoice;
   mChoice=mpList->GetSelection();
   mpOption->SetChoice(mChoice);
}

void WXFieldOption::CrystUpdate(const bool uui,const bool lock)
{
   VFN_DEBUG_MESSAGE("WXFieldOption::CrystUpdate("<<uui<<lock<<")",6)
   if(lock) mMutex.Lock();
   if(mChoice==mpOption->GetChoice())
   {
      if(lock) mMutex.Unlock();
      return;
   }
   mChoice=mpOption->GetChoice();
   mNeedUpdateUI=true;
   if(lock) mMutex.Unlock();
}

void WXFieldOption::UpdateUI(const bool lock)
{
   VFN_DEBUG_MESSAGE("WXFieldOption::UpdateUI("<<lock<<")"<<"MainThread="<<wxThread::IsMain(),6)
   if(lock) mMutex.Lock();
   if(mNeedUpdateUI==false)
   {
      if(lock) mMutex.Unlock();
      return;
   }
   mpList->SetSelection(mChoice);
   mNeedUpdateUI=false;
   if(lock) mMutex.Unlock();
}

void WXFieldOption::Revert()
{
   wxMutexLocker mlock(mMutex);
   mChoice=mChoiceOld;
   mNeedUpdateUI=true;
}
void WXFieldOption::ValidateUserInput()
{
}

void WXFieldOption::SetToolTip(const wxString& tip){mpList->SetToolTip(tip);}

////////////////////////////////////////////////////////////////////////
//
//    WXRegistry
//
////////////////////////////////////////////////////////////////////////

template<class T> WXRegistry<T>::WXRegistry(wxWindow *parent,ObjRegistry<T>* reg):
WXCrystObj(parent,wxHORIZONTAL,false),mpRegistry(reg)
{
   VFN_DEBUG_ENTRY("WXCrystRegistry::WXCrystRegistry(wxWindow*)",6)
   #ifdef VFN_CRYST_MUTEX
   cout <<"new CrystMutex("<<&mMutex<<")for WXCrystRegistry:"<<reg->GetName()<<endl;
   #endif
   wxStaticText* mpLabel=new wxStaticText(this,-1,wxString::FromAscii(reg->GetName().c_str()));
   mpSizer->Add(mpLabel,0,wxALIGN_LEFT);
   mpLabel->SetForegroundColour(wxColour(0,0,255));
   this->BottomLayout(0);
   VFN_DEBUG_EXIT("WXCrystRegistry::WXCrystRegistry(wxWindow*):End",6)
}
template<class T> WXRegistry<T>::~WXRegistry()
{
   mpRegistry->WXNotifyDelete();
   #ifdef VFN_CRYST_MUTEX
   cout <<"Deleting CrystMutex("<<&mMutex<<")for WXCrystRegistry:"<<mpRegistry->GetName()<<endl;
   #endif
}
template<class T> void WXRegistry<T>::Add(WXCrystObjBasic *obj)
{
   VFN_DEBUG_ENTRY("WXCrystRegistry::AddWXCrystObj(WXCrystObj*)",6)
   mList.Add(obj);
   obj->Show(mIsExpanded);
   this->AddChild(obj);
   VFN_DEBUG_EXIT("WXCrystRegistry::AddWXCrystObj(WXCrystObj*)",6)
}
template<class T> void WXRegistry<T>::Remove(WXCrystObjBasic *obj)
{
   if(obj==0) return;
   VFN_DEBUG_ENTRY("WXCrystRegistry::RemoveWXCrystObj(WXCrystObj*)",6)
   mList.Remove(obj);
   mpSizer->Detach((wxWindow*)obj);
   obj->Destroy();
   this->BottomLayout(0);
   VFN_DEBUG_EXIT("WXCrystRegistry::RemoveWXCrystObj(WXCrystObj*):End",6)
}

template<class T> bool WXRegistry<T>::OnChangeName(const int id)
{
   VFN_DEBUG_MESSAGE("WXRegistry<T>::OnChangeName()",6)
   if(id==ID_WXOBJ_NAME)
   {
      mpRegistry->SetName(mpWXTitle->GetValue());
      return true;
   }
   return false;
}


//Explicit instantiation
template class WXRegistry<RefinableObj>;
template class WXRegistry<RefObjOpt>;
template class WXRegistry<Crystal>;
template class WXRegistry<Scatterer>;
template class WXRegistry<ScatteringPower>;
template class WXRegistry<ScatteringPowerAtom>;
template class WXRegistry<PowderPattern>;
template class WXRegistry<PowderPatternComponent>;
template class WXRegistry<DiffractionDataSingleCrystal>;
template class WXRegistry<OptimizationObj>;
template class WXRegistry<XMLCrystTag>;
//template class WXRegistry<IOCrystTag>;//to be removed
template class WXRegistry<ZAtom>;
template class WXRegistry<TexturePhaseMarchDollase>;
template class WXRegistry<TextureEllipsoid>;
template class WXRegistry<ReflectionProfile>;
////////////////////////////////////////////////////////////////////////
//
//    WXDialogChooseFromRegistry
//
////////////////////////////////////////////////////////////////////////
template<class T> T* WXDialogChooseFromRegistry(ObjRegistry<T> &reg,wxWindow*parent,
                                                const string &message,int &choice)
{
   wxString* choices=new wxString[reg.GetNb()];
   for(int i=0;i<reg.GetNb();i++) 
      *(choices+i)=wxString::FromAscii((reg.GetObj(i).GetClassName()+":"+reg.GetObj(i).GetName()).c_str());
   #if wxCHECK_VERSION( 2, 9, 4 )
   wxSingleChoiceDialog dialog
         (parent,wxString::FromAscii(message.c_str()),_T("Choose"),reg.GetNb(),choices,(void**)NULL,wxOK | wxCANCEL);
   #else
   wxSingleChoiceDialog dialog
         (parent,wxString::FromAscii(message.c_str()),_T("Choose"),reg.GetNb(),choices,NULL,wxOK | wxCANCEL);
   #endif
   dialog.SetSize(300,300);
   if(wxID_OK!=dialog.ShowModal())
   {
      delete[] choices;
      return 0;
   }
   delete[] choices;
   choice=dialog.GetSelection();
   return &(reg.GetObj(choice));
}

template RefinableObj* 
   WXDialogChooseFromRegistry(ObjRegistry<RefinableObj> &,wxWindow*,const string &,int &);
template Crystal* 
   WXDialogChooseFromRegistry(ObjRegistry<Crystal> &,wxWindow*,const string &,int &);
template Scatterer* 
   WXDialogChooseFromRegistry(ObjRegistry<Scatterer> &,wxWindow*,const string &,int &);
template ScatteringPower* 
   WXDialogChooseFromRegistry(ObjRegistry<ScatteringPower> &,wxWindow*,const string &,int &);
template ScatteringPowerAtom* 
   WXDialogChooseFromRegistry(ObjRegistry<ScatteringPowerAtom> &,wxWindow*,
                              const string &,int &);
template ZAtom* 
   WXDialogChooseFromRegistry(ObjRegistry<ZAtom> &,wxWindow*,const string &,int &);
template PowderPattern* 
   WXDialogChooseFromRegistry(ObjRegistry<PowderPattern> &,wxWindow*,const string &,int &);
template PowderPatternComponent* 
   WXDialogChooseFromRegistry(ObjRegistry<PowderPatternComponent>&,wxWindow*,
                              const string&,int &);
template DiffractionDataSingleCrystal* 
   WXDialogChooseFromRegistry(ObjRegistry<DiffractionDataSingleCrystal>&,wxWindow*,
                              const string &,int &);
template OptimizationObj* 
   WXDialogChooseFromRegistry(ObjRegistry<OptimizationObj> &,wxWindow*,const string &,int &);
template XMLCrystTag* 
   WXDialogChooseFromRegistry(ObjRegistry<XMLCrystTag> &,wxWindow*,const string &,int &);


template<class T> const T* WXDialogChooseFromRegistry(const ObjRegistry<T> &reg,
                                                      wxWindow*parent,const string &message
                                                      ,int &choice)
{
   wxString* choices=new wxString[reg.GetNb()];
   for(int i=0;i<reg.GetNb();i++) 
      *(choices+i)=wxString::FromAscii((reg.GetObj(i).GetClassName()+":"+reg.GetObj(i).GetName()).c_str());
   #if wxCHECK_VERSION( 2, 9, 4 )
   wxSingleChoiceDialog dialog
         (parent,wxString::FromAscii(message.c_str()),_T("Choose"),reg.GetNb(),choices,(void**)NULL,wxOK | wxCANCEL);
   #else
   wxSingleChoiceDialog dialog
         (parent,wxString::FromAscii(message.c_str()),_T("Choose"),reg.GetNb(),choices,NULL,wxOK | wxCANCEL);
   #endif
   dialog.SetSize(300,300);
   if(wxID_OK!=dialog.ShowModal())
   {
      delete[] choices;
      return 0;
   }
   delete[] choices;
   choice=dialog.GetSelection();
   return &(reg.GetObj(choice));
}

template const RefinableObj* 
   WXDialogChooseFromRegistry(const ObjRegistry<RefinableObj> &,wxWindow*,const string &,int &);
   
template const Crystal* 
   WXDialogChooseFromRegistry(const ObjRegistry<Crystal> &,wxWindow*,const string &,int &);
   
template const Scatterer* 
   WXDialogChooseFromRegistry(const ObjRegistry<Scatterer> &,wxWindow*,const string &,int &);
   
template const ScatteringPower* 
   WXDialogChooseFromRegistry(const ObjRegistry<ScatteringPower> &,wxWindow*,
                              const string &,int &);
                              
template const ScatteringPowerAtom* 
   WXDialogChooseFromRegistry(const ObjRegistry<ScatteringPowerAtom> &,wxWindow*,
                              const string &,int &);
                              
template const ZAtom* 
   WXDialogChooseFromRegistry(const ObjRegistry<ZAtom> &,wxWindow*,const string &,int &);
   
template const PowderPattern* 
   WXDialogChooseFromRegistry(const ObjRegistry<PowderPattern> &,wxWindow*,
                              const string &,int &);
                              
template const PowderPatternComponent* 
   WXDialogChooseFromRegistry(const ObjRegistry<PowderPatternComponent>&,wxWindow*,
                              const string&,int &);
                              
template const DiffractionDataSingleCrystal* 
   WXDialogChooseFromRegistry(const ObjRegistry<DiffractionDataSingleCrystal>&,wxWindow*,
                              const string &,int &);
                              
template const OptimizationObj* 
   WXDialogChooseFromRegistry(const ObjRegistry<OptimizationObj> &,wxWindow*,
                              const string &,int &);
                              
template const XMLCrystTag* 
   WXDialogChooseFromRegistry(const ObjRegistry<XMLCrystTag> &,wxWindow*,
                              const string &,int &);

////////////////////////////////////////////////////////////////////////
//
//    WXRefinableObj
//
////////////////////////////////////////////////////////////////////////
BEGIN_EVENT_TABLE(WXRefinableObj,wxWindow)
   EVT_MENU(ID_REFOBJ_MENU_OBJ_SAVE,        WXRefinableObj::OnMenuSave)
   EVT_MENU(ID_REFOBJ_MENU_OBJ_LOAD,        WXRefinableObj::OnMenuLoad)
   EVT_MENU(ID_REFOBJ_MENU_PAR_FIXALL,      WXRefinableObj::OnMenuFixAllPar)
   EVT_MENU(ID_REFOBJ_MENU_PAR_UNFIXALL,    WXRefinableObj::OnMenuUnFixAllPar)
   EVT_MENU(ID_REFOBJ_MENU_PAR_RANDOMIZE,   WXRefinableObj::OnMenuParRandomize)
   EVT_UPDATE_UI(ID_CRYST_UPDATEUI,         WXRefinableObj::OnUpdateUI)
END_EVENT_TABLE()

WXRefinableObj::WXRefinableObj(wxWindow* parent, RefinableObj*obj):
WXCrystObj(parent,wxHORIZONTAL),mpRefinableObj(obj)
{
   VFN_DEBUG_ENTRY("WXRefinableObj::WXRefinableObj():"<<obj->GetName(),6)
   #ifdef VFN_CRYST_MUTEX
   cout <<"new CrystMutex("<<&mMutex<<")for :"<<obj->GetClassName()<<":"<<obj->GetName()<<endl;
   #endif
   mpWXTitle->SetLabel(mpRefinableObj->GetClassName());
   
   // Menu
      mpMenuBar=new WXCrystMenuBar(this,this);
   
      mpSizer->Add(mpMenuBar);
      mList.Add(mpMenuBar);
   //:TODO: Rather use a WXRegistry for the options ?
   for(unsigned int i=0;i<mpRefinableObj->GetNbOption();i++)
   {
   VFN_DEBUG_MESSAGE("WXRefinableObj::WXRefinableObj():Adding option "<<i,6)
      WXCrystObjBasic *opt=mpRefinableObj->GetOption(i).WXCreate(this);
      mpSizer->Add(opt,0,wxALIGN_LEFT);
      mList.Add(opt);
   }
   //this->BottomLayout(0);
   //this->CrystUpdate(true,true);
   VFN_DEBUG_EXIT("WXRefinableObj::WXRefinableObj()",6)
}

WXRefinableObj::~WXRefinableObj()
{
   VFN_DEBUG_MESSAGE("WXRefinableObj::~WXRefinableObj():"<<mpRefinableObj->GetName(),6)
   mpRefinableObj->WXNotifyDelete();
   #ifdef VFN_CRYST_MUTEX
   cout <<"Deleting CrystMutex("<<&mMutex<<")for :"<<mpRefinableObj->GetClassName()
        <<":"<<mpRefinableObj->GetName()<<endl;
   #endif
}

void WXRefinableObj::CrystUpdate(const bool uui,const bool lock)
{
   VFN_DEBUG_ENTRY("WXRefinableObj::CrystUpdate("<<uui<<lock<<"):"<<mpRefinableObj->GetName(),6)
   //if(lock) mMutex.Lock();
   this->WXCrystObj::CrystUpdate(false,lock);
   //if(lock) mMutex.Unlock();
   if(uui)
   {
      if(true==wxThread::IsMain()) this->UpdateUI(lock);
      else
      {
         wxUpdateUIEvent event(ID_CRYST_UPDATEUI);
         wxPostEvent(this,event);
      }
   }
   VFN_DEBUG_EXIT("WXRefinableObj::CrystUpdate():"<<mpRefinableObj->GetName(),6)
}

bool WXRefinableObj::OnChangeName(const int id)
{
   VFN_DEBUG_MESSAGE("WXRefinableObj::OnChangeName()",6)
   if(id==ID_WXOBJ_NAME)
   {
   VFN_DEBUG_MESSAGE("WXRefinableObj::OnChangeName():Changing RefinableObj Name",6)
      mpRefinableObj->SetName(mpWXTitle->GetValue());
      return true;
   }
   return false;
}

void WXRefinableObj::OnMenuSave(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXRefinableObj::OnButtonSave()",6)
   wxFileDialog save(this,_T("Choose a file"),_T(""),_T(""),_T("*.xml"),wxFD_SAVE | wxFD_OVERWRITE_PROMPT);
   if(save.ShowModal() != wxID_OK) return;
   
   ofstream out(save.GetPath().ToAscii());
   if(!out) return;//:TODO:
   {
      mpRefinableObj->XMLOutput(out);
   }
   out.close();
}

void WXRefinableObj::OnMenuLoad(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXRefinableObj::OnButtonLoad()",6)
   wxFileDialog *open= new wxFileDialog(this,_T("Choose a file"),_T(""),_T(""),_T("*.xml"),
                                        wxFD_OPEN | wxFD_FILE_MUST_EXIST);
   if(open->ShowModal() != wxID_OK) return;
   
   ifstream fin(open->GetPath().ToAscii());
   if(!fin) return;//:TODO:
   {
      XMLCrystTag tag(fin);//:TODO: load all tags and find the right ones for this class
      mpRefinableObj->XMLInput(fin,tag);
   }
   fin.close();
   open->Destroy();
}

void WXRefinableObj::OnMenuFixAllPar(wxCommandEvent & WXUNUSED(event))
{
   mpRefinableObj->FixAllPar();
   this->CrystUpdate(true);
}

void WXRefinableObj::OnMenuUnFixAllPar(wxCommandEvent & WXUNUSED(event))
{
   mpRefinableObj->UnFixAllPar();
   this->CrystUpdate(true);
}

void WXRefinableObj::OnMenuParRandomize(wxCommandEvent & WXUNUSED(event))
{
   mpRefinableObj->RandomizeConfiguration();
   this->CrystUpdate(true);
}

void WXRefinableObj::OnUpdateUI(wxUpdateUIEvent& event)
{
   this->UpdateUI(true);
   event.Skip();
}
void WXRefinableObj::UpdateUI(const bool lock)
{
   VFN_DEBUG_ENTRY("WXRefinableObj::UpdateUI("<<lock<<")"<<"MainThread="<<wxThread::IsMain()<<":"<<mpRefinableObj->GetName(),6)
   if(!mpRefinableObj->IsBeingRefined())
   {
      if(lock) mMutex.Lock();
      mpWXTitle->SetValue(mpRefinableObj->GetName());
      mpWXTitle->UpdateUI(false);
      //mpSizer->SetItemMinSize
      //      (mpWXTitle, mpWXTitle->GetSize().GetWidth(),mpWXTitle->GetSize().GetHeight());
      //mpWXTitle->Layout();
      if(lock) mMutex.Unlock();
   }
   this->WXCrystObj::UpdateUI(lock);
   VFN_DEBUG_EXIT("WXRefinableObj::UpdateUI()"<<mpRefinableObj->GetName(),6)
}

}// namespace 
