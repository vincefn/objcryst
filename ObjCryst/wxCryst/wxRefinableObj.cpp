//#include <sstream> //for stringstream
#include <fstream>

#include "wx/wx.h"

#include "wxCryst/wxRefinableObj.h"

#include "wxCryst/wxRefinableObj.h"
#include "Quirks/VFNStreamFormat.h"

#include "RefinableObj/GlobalOptimObj.h"

//These are only for explicit instantiation
#include "ObjCryst/Atom.h"
#include "ObjCryst/Crystal.h"
#include "ObjCryst/DiffractionDataSingleCrystal.h"
#include "ObjCryst/PowderPattern.h"
#include "ObjCryst/ScatteringPower.h"
#include "ObjCryst/ZScatterer.h"

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
BEGIN_EVENT_TABLE(WXFieldRefPar,wxEvtHandler)
   EVT_TEXT_ENTER(ID_WXFIELD,                   WXFieldRefPar::OnEnter)
   EVT_TEXT(ID_WXFIELD,                   		WXFieldRefPar::OnText)
   EVT_CHECKBOX(ID_WXFIELD_REFPAR_FIXBUTTON,    WXFieldRefPar::OnToggleFix)
   EVT_RIGHT_DOWN(                              WXFieldRefPar::OnPopupMenu)
   EVT_UPDATE_UI(ID_WXFIELD_REFPAR,             WXFieldRefPar::OnUpdateUI)
   EVT_MENU(ID_REFPAR_POPUP_SET_LIMITS,         WXFieldRefPar::OnPopupMenuChoice)
   EVT_MENU(ID_REFPAR_POPUP_REMOVE_LIMITS,      WXFieldRefPar::OnPopupMenuChoice)
END_EVENT_TABLE()

WXFieldRefPar::WXFieldRefPar(wxWindow *parent,const string& label,
                         RefinablePar *par, const int hsize):
WXField(parent,label,ID_WXFIELD_REFPAR),mValue(0.),mpRefPar(par),mIsSelfUpdating(false)
{
   VFN_DEBUG_MESSAGE("WXFieldRefPar::WXFieldName():End",6)

   mpButtonFix=new wxCheckBox(this,ID_WXFIELD_REFPAR_FIXBUTTON,"");
   //mpButtonFix->PushEventHandler(this);
   mpSizer->Add(mpButtonFix,0,wxALIGN_CENTER);
   
   mpField=new wxTextCtrl(this,ID_WXFIELD,"",
                            wxDefaultPosition,wxSize(hsize,-1),wxTE_PROCESS_ENTER,
                            wxTextValidator(wxFILTER_NUMERIC));
   //mpField->PushEventHandler(this);
   mpSizer->Add(mpField,0,wxALIGN_CENTER);

   mpPopUpMenu=new wxMenu("Refinable Parameter");
   mpPopUpMenu->Append(ID_REFPAR_POPUP_SET_LIMITS, "Set Limits");
   mpPopUpMenu->Append(ID_REFPAR_POPUP_REMOVE_LIMITS, "Remove Limits");
   
   this->Layout();
}

void WXFieldRefPar::OnUpdateUI(wxUpdateUIEvent & WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXFieldRefPar::OnUpdateUI()"<<mValue,3)
   //The 'fix' state should not change during refinement so this should be safe
   if(false==mpRefPar->IsUsed()) this->Show(false);
   else this->Show(true);
   
   if(mpField==0) return;
   
   //mpField->SetValue(wxString::Printf("%f",mValue));
   wxString tmp;
   tmp.Printf("%f",mValue);
	mIsSelfUpdating=true;
   mpField->SetValue(tmp);
	mIsSelfUpdating=false;
   mpButtonFix->SetValue(!(mpRefPar->IsFixed()));
   VFN_DEBUG_MESSAGE("WXFieldRefPar::OnUpdateUI():End",2)
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
   mpRefPar->SetIsFixed(!(mpButtonFix->GetValue()));
   mpRefPar->Print();
}

void WXFieldRefPar::OnPopupMenu(wxCommandEvent & WXUNUSED(event))
{
   this->PopupMenu(mpPopUpMenu,0,0);
}

void WXFieldRefPar::OnPopupMenuChoice(wxMenuEvent& event)
{
   VFN_DEBUG_MESSAGE("WXFieldRefPar::OnPopupMenuChoice()",7)
   switch(event.GetId())
   {
      case ID_REFPAR_POPUP_SET_LIMITS:
      {
         double min,max;
         {
            wxString str;
            str << mpRefPar->GetHumanMin();
            wxTextEntryDialog limitDialog(this,"Enter the minimum value",
                                    "Minimum",str,wxOK | wxCANCEL);
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
            wxTextEntryDialog limitDialog(this,"Enter the maximum value",
                                    "Maximum",str,wxOK | wxCANCEL);
            if(wxID_OK!=limitDialog.ShowModal())
            {
               VFN_DEBUG_EXIT("WXZScatterer::OnMenuSetLimits():Cancelled",6)
               return;
            }
            limitDialog.GetValue().ToDouble(&max);
         }
         if(max<=min)
         {
            wxMessageDialog dumbUser(this,"max <= min !!!",
                                     "Whooops",wxOK|wxICON_EXCLAMATION);
            dumbUser.ShowModal();
            return;
         }
         mpRefPar->SetHumanMin(min);
         mpRefPar->SetHumanMax(max);
         mpRefPar->SetIsLimited(true);
         mpRefPar->Print();
         break;
      }
      case ID_REFPAR_POPUP_REMOVE_LIMITS:mpRefPar->SetIsLimited(false);mpRefPar->Print();break;
   }
 
}

void WXFieldRefPar::CrystUpdate()
{
   VFN_DEBUG_MESSAGE("WXFieldRefPar::CrystUpdate()",6)
   //cout << mpField <<endl;
   mValue=mpRefPar->GetHumanValue();
   wxUpdateUIEvent event(ID_WXFIELD_REFPAR);
   wxPostEvent(this,event);
}
void WXFieldRefPar::Revert()
{
   VFN_DEBUG_MESSAGE("WXFieldRefPar::Revert()",6)
   mValue=mValueOld;
   wxUpdateUIEvent event(ID_WXFIELD_REFPAR);
   wxPostEvent(this,event);
}
void WXFieldRefPar::ValidateUserInput()
{
   VFN_DEBUG_MESSAGE("WXFieldRefPar::ValidateUserInput()",6)
   mValueOld=mValue;
   wxString s=mpField->GetValue();
	double tmp;
   s.ToDouble(&tmp);
	mValue=tmp;
   mpRefPar->SetHumanValue(mValue);
}
////////////////////////////////////////////////////////////////////////
//
//    WXFieldOption
//
////////////////////////////////////////////////////////////////////////
BEGIN_EVENT_TABLE(WXFieldOption,wxEvtHandler)
   EVT_CHOICE(ID_WXFIELD,WXFieldOption::OnChoice)
   EVT_UPDATE_UI(ID_WXFIELD_OPTION,WXFieldOption::OnUpdateUI)
END_EVENT_TABLE()

WXFieldOption::WXFieldOption(wxWindow *parent,
                             const int field_id,RefObjOpt* option):
WXField(parent,option->GetName(),field_id),mpOption(option)
{
   wxString choices[20];//:TODO: dynamically choose correct number
   for(int i=0;i<mpOption->GetNbChoice();i++)
      choices[i]=mpOption->GetChoiceName(i).c_str();
   
   mpList= new wxChoice(this,ID_WXFIELD,wxDefaultPosition,wxDefaultSize,
                        mpOption->GetNbChoice(),choices);
   mpSizer->Add(mpList,0,wxALIGN_CENTER);
   this->Layout();
}
WXFieldOption::~WXFieldOption()
{
   mpOption->WXNotifyDelete();
}
void WXFieldOption::OnUpdateUI(wxUpdateUIEvent & WXUNUSED(event))
{
   mpList->SetSelection(mChoice);
}

void WXFieldOption::OnChoice(wxCommandEvent & WXUNUSED(event))
{
   if(mChoice==mpList->GetSelection()) return;
   mChoiceOld=mChoice;
   mChoice=mpList->GetSelection();
   mpOption->SetChoice(mChoice);
}

void WXFieldOption::CrystUpdate()
{
   VFN_DEBUG_MESSAGE("WXFieldOption::CrystUpdate()",6)
   // The list of choices cannot change once initialized
   mChoice=mpOption->GetChoice();
   wxUpdateUIEvent event(ID_WXFIELD_OPTION);
   wxPostEvent(this,event);
}

void WXFieldOption::Revert()
{
   mChoice=mChoiceOld;
   wxUpdateUIEvent event(ID_WXFIELD_OPTION);
   wxPostEvent(this,event);
}
void WXFieldOption::ValidateUserInput()
{
}
////////////////////////////////////////////////////////////////////////
//
//    WXCostFunction
//
////////////////////////////////////////////////////////////////////////
BEGIN_EVENT_TABLE(WXCostFunction,wxWindow)
   EVT_UPDATE_UI(ID_WXFIELD_COSTFUNC,WXCostFunction::OnUpdateUI)
END_EVENT_TABLE()

WXCostFunction::WXCostFunction(wxWindow *parent,RefinableObj *obj, const int field_id,
               const int funcNum,REAL * weight):
WXField(parent,obj->GetName()+":"+obj->GetCostFunctionName(funcNum)+"=",-1),
mpObj(obj),mFuncNum(funcNum)
{
   mpValue=new wxTextCtrl(this,ID_WXFIELD,"",wxDefaultPosition,wxDefaultSize,wxTE_READONLY);
   mpSizer->Add(mpValue,0,wxALIGN_CENTER);
   
   mpWeight=new WXFieldPar<REAL>(this,",weight=",-1,weight);
   mpSizer->Add(mpWeight,0,wxALIGN_CENTER);
   
   this->CrystUpdate();
   this->Layout();
}

void WXCostFunction::OnUpdateUI(wxUpdateUIEvent & WXUNUSED(event))
{
   wxString tmp;
   tmp.Printf("%f",mValue);
   mpValue->SetValue(tmp);
}

void WXCostFunction::OnEnter(wxCommandEvent & WXUNUSED(event))
{
   //Nothing to do here (will never happen, value is not editable)
}

void WXCostFunction::CrystUpdate()
{
   mValue=mpObj->GetCostFunctionValue(mFuncNum);
   wxUpdateUIEvent event(ID_WXFIELD_COSTFUNC);
   wxPostEvent(this,event);
}

void WXCostFunction::Revert()
{
   //Nothing to do here
}
void WXCostFunction::ValidateUserInput()
{
}

////////////////////////////////////////////////////////////////////////
//
//    WXRegistry
//
////////////////////////////////////////////////////////////////////////

template<class T> WXRegistry<T>::WXRegistry(wxWindow *parent,ObjRegistry<T>* reg):
WXCrystObj(parent,wxHORIZONTAL,false),mpRegistry(reg)
{
   VFN_DEBUG_MESSAGE("WXCrystRegistry::WXCrystRegistry(wxWindow*)",6)
   wxStaticText* mpLabel=new wxStaticText(this,-1,reg->GetName().c_str());
   mpSizer->Add(mpLabel,0,wxALIGN_LEFT);
   mpLabel->SetBackgroundColour(wxColour(250,250,200));
   mpLabel->SetForegroundColour(wxColour(0,0,255));

   //mpWXTitle->SetValue(mpRegistry->GetName());
   //mpWXTitle->SetBackgroundColour(wxColour(250,250,200));
   this->Layout();
   VFN_DEBUG_MESSAGE("WXCrystRegistry::WXCrystRegistry(wxWindow*):End",6)
}
template<class T> WXRegistry<T>::~WXRegistry(){mpRegistry->WXNotifyDelete();}

template<class T> void WXRegistry<T>::Add(WXCrystObjBasic *obj)
{
   VFN_DEBUG_MESSAGE("WXCrystRegistry::AddWXCrystObj(WXCrystObj*)",6)
   mpSizer->Add(obj,0,wxALIGN_LEFT);
   mList.Add(obj);
   obj->Show(mIsExpanded);
   this->Layout();
}
template<class T> void WXRegistry<T>::Remove(WXCrystObjBasic *obj)
{
   VFN_DEBUG_MESSAGE("WXCrystRegistry::RemoveWXCrystObj(WXCrystObj*)",6)
   mList.Remove(obj);
   mpSizer->Remove(obj);
   obj->Destroy();
   this->Layout();
   VFN_DEBUG_MESSAGE("WXCrystRegistry::RemoveWXCrystObj(WXCrystObj*):End",6)
}

template<class T> bool WXRegistry<T>::OnChangeName(const int id)
{
   VFN_DEBUG_MESSAGE("WXRegistry<T>::OnChangeName()",6)
   switch(id)
   {
      case ID_WXOBJ_NAME:
      {
         mpRegistry->SetName(mpWXTitle->GetValue());
         return true;
      }
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
////////////////////////////////////////////////////////////////////////
//
//    WXDialogChooseFromRegistry
//
////////////////////////////////////////////////////////////////////////
template<class T> T* WXDialogChooseFromRegistry(ObjRegistry<T> &reg,wxWindow*parent,
                                                const string &message,int &choice)
{
   wxString choices[50];//:TODO: give correctt number...
   for(int i=0;i<reg.GetNb();i++) 
      choices[i]=(reg.GetObj(i).GetClassName()+":"+reg.GetObj(i).GetName()).c_str();
   wxSingleChoiceDialog dialog
         (parent,message.c_str(),"Choose",reg.GetNb(),choices,0,wxOK | wxCANCEL);
   dialog.SetSize(300,300);
   if(wxID_OK!=dialog.ShowModal()) return 0;
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
   wxString choices[50]; //:TODO: 
   for(int i=0;i<reg.GetNb();i++) 
      choices[i]=(reg.GetObj(i).GetClassName()+":"+reg.GetObj(i).GetName()).c_str();
   wxSingleChoiceDialog dialog
         (parent,message.c_str(),"Choose",reg.GetNb(),choices,0,wxOK | wxCANCEL);
   if(wxID_OK!=dialog.ShowModal()) return 0;
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
   EVT_MENU(ID_REFOBJ_MENU_PAR_RANDOMIZE,   WXRefinableObj::OnMenuUnFixAllPar)
   EVT_UPDATE_UI(ID_CRYST_UPDATEUI, 		  WXRefinableObj::OnUpdateUI)
END_EVENT_TABLE()

WXRefinableObj::WXRefinableObj(wxWindow* parent, RefinableObj*obj):
WXCrystObj(parent,wxHORIZONTAL),mpRefinableObj(obj)
{
   VFN_DEBUG_MESSAGE("WXRefinableObj::WXRefinableObj():"<<obj->GetName(),6)
   mpWXTitle->SetLabel(mpRefinableObj->GetClassName());
   mpWXTitle->Layout();
   //cout<<"0:"<<mpWXTitle->GetSize().GetWidth()<<":"<<mpWXTitle->GetSize().GetHeight()<<endl;
   mpSizer->SetItemMinSize
               (mpWXTitle, mpWXTitle->GetSize().GetWidth(),mpWXTitle->GetSize().GetHeight());
   mpWXTitle->Layout();
   // Menu
      mpMenuBar=new WXCrystMenuBar(this,this);
   
      mpSizer->Add(mpMenuBar);
      mList.Add(mpMenuBar);
   //:TODO: Rather use a WXRegistry for the options ?
   for(unsigned int i=0;i<mpRefinableObj->GetNbOption();i++)
   {
   VFN_DEBUG_MESSAGE("WXRefinableObj::WXRefinableObj():Adding option "<<i,6)
      WXFieldOption *opt=new WXFieldOption(this,-1,&(mpRefinableObj->GetOption(i)));
      mpSizer->Add(opt,0,wxALIGN_LEFT);
      mList.Add(opt);
   }
   this->CrystUpdate();
   this->Layout();
   VFN_DEBUG_MESSAGE("WXRefinableObj::WXRefinableObj():End",6)
}

WXRefinableObj::~WXRefinableObj()
{
   VFN_DEBUG_MESSAGE("WXRefinableObj::~WXRefinableObj():"<<mpRefinableObj->GetName(),6)
   mpRefinableObj->WXNotifyDelete();
}

bool WXRefinableObj::Layout()
{
   VFN_DEBUG_MESSAGE("WXRefinableObj::Layout():"<<mpRefinableObj->GetName(),6)
   mpMenuBar->Layout();
   mpSizer->SetItemMinSize(mpMenuBar,
                           mpMenuBar->GetSize().GetWidth(),
                           mpMenuBar->GetSize().GetHeight());
   return this->WXCrystObj::Layout();
}

void WXRefinableObj::CrystUpdate()
{
   VFN_DEBUG_MESSAGE("WXRefinableObj::CrystUpdate():"<<mpRefinableObj->GetName(),6)
   wxUpdateUIEvent event(ID_CRYST_UPDATEUI);
   wxPostEvent(this,event);
   this->WXCrystObj::CrystUpdate();
}

bool WXRefinableObj::OnChangeName(const int id)
{
   VFN_DEBUG_MESSAGE("WXRefinableObj::OnChangeName()",6)
   switch(id)
   {
      case ID_WXOBJ_NAME:
      {
      VFN_DEBUG_MESSAGE("WXRefinableObj::OnChangeName():Changing RefinableObj Name",6)
         mpRefinableObj->SetName(mpWXTitle->GetValue());
         return true;
      }
   }
   return false;
}

void WXRefinableObj::OnMenuSave(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXRefinableObj::OnButtonSave()",6)
   wxFileDialog save(this,"Choose a file","","","*.xml",wxSAVE | wxOVERWRITE_PROMPT);
   if(save.ShowModal() != wxID_OK) return;
   
   ofstream out(save.GetPath().c_str());
   if(!out) return;//:TODO:
   {
      mpRefinableObj->XMLOutput(out);
   }
   out.close();
}

void WXRefinableObj::OnMenuLoad(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXRefinableObj::OnButtonLoad()",6)
   wxFileDialog *open= new wxFileDialog(this,"Choose a file","","","*.xml",
                                        wxOPEN | wxFILE_MUST_EXIST);
   if(open->ShowModal() != wxID_OK) return;
   
   ifstream fin(open->GetPath().c_str());
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
   this->CrystUpdate();
}

void WXRefinableObj::OnMenuUnFixAllPar(wxCommandEvent & WXUNUSED(event))
{
   mpRefinableObj->UnFixAllPar();
   this->CrystUpdate();
}

void WXRefinableObj::OnMenuParRandomize(wxCommandEvent & WXUNUSED(event))
{
   mpRefinableObj->RandomizeConfiguration();
   mpRefinableObj->RefinableObj::Print();
   this->CrystUpdate();
}

void WXRefinableObj::OnUpdateUI(wxUpdateUIEvent& event)
{
	mpWXTitle->SetValue(mpRefinableObj->GetName());
}

}// namespace 
