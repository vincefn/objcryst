//#include <sstream> //for stringstream
#include <fstream>

#include "wxCryst/wxCryst.h"
#include "wx/wx.h"

//#include "Quirks/VFNStreamFormat.h"
#include "Quirks/VFNDebug.h"

//#include "RefinableObj/GlobalOptimObj.h"

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
////////////////////////////////////////////////////////////////////////
//
//    WXCrystObjBasic
//
////////////////////////////////////////////////////////////////////////
WXCrystObjBasic::WXCrystObjBasic(wxWindow* parent):
wxWindow(parent,-1),mWXParent(parent),mIsShown(true),mNeedUpdateUI(true)
{
   VFN_DEBUG_MESSAGE("WXCrystObjBasic::WXCrystObjBasic() at "<<this,6)
   //mWXParent->Layout();
   VFN_DEBUG_MESSAGE("WXCrystObjBasic::WXCrystObjBasic():End",6)
}

WXCrystObjBasic::~WXCrystObjBasic(){}
//void WXCrystObjBasic::CrystUpdate()
//{
//   cout <<"Just called  WXCrystObjBasic::CrystUpdate(), which is pure virtual !!"<<endl;
//   abort();
//}

//wxWindow* WXCrystObjBasic::GetParent(){return mWXParent;}

////////////////////////////////////////////////////////////////////////
//
//    WXCrystObjBasicList
//
////////////////////////////////////////////////////////////////////////
WXCrystObjBasicList::WXCrystObjBasicList():
mNbWXCrystObj(0),mMaxNbWXCrystObj(32),mpWXCrystObj(0)
{
   mpWXCrystObj= new WXCrystObjBasic*[mMaxNbWXCrystObj];
}

WXCrystObjBasicList::~WXCrystObjBasicList()
{
   delete[] mpWXCrystObj;
}

unsigned int WXCrystObjBasicList::GetNb()const {return mNbWXCrystObj;}

void WXCrystObjBasicList::Add(WXCrystObjBasic *win)
{
   VFN_DEBUG_MESSAGE("WXCrystObjBasicList::Add()",6)
   if(mNbWXCrystObj==mMaxNbWXCrystObj)
   {
      WXCrystObjBasic** tmp= new WXCrystObjBasic*[mMaxNbWXCrystObj+16];
      for(unsigned int i=0;i<mNbWXCrystObj;i++) tmp[i]=mpWXCrystObj[i];
      delete[] mpWXCrystObj;
      mpWXCrystObj= tmp;
		mMaxNbWXCrystObj+=16;
   }
   mpWXCrystObj[mNbWXCrystObj++]=win;
}

void WXCrystObjBasicList::Remove(const WXCrystObjBasic *win)
{
   VFN_DEBUG_MESSAGE("WXCrystObjBasicList::Remove():"<<win,6)
   unsigned int i;
   for(i=0;i<mNbWXCrystObj;i++) if(mpWXCrystObj[i]==win) break;
   if(i<(mNbWXCrystObj-1)) mpWXCrystObj[i]=mpWXCrystObj[--mNbWXCrystObj];
   else
   {
      if(i==(mNbWXCrystObj-1)) mNbWXCrystObj--;
      else
      {
         if(i==mNbWXCrystObj)
         {
            cout << "Trying to remove a non-existing window...aborting"<<endl;
            for(i=0;i<mNbWXCrystObj;i++) cout <<i<<":"<<mpWXCrystObj[i]<<endl;
            abort();
         }
      }
   }
}

bool WXCrystObjBasicList::Show(bool show)
{
   VFN_DEBUG_MESSAGE("WXCrystObjBasicList::Show(bool)",3)
   for(unsigned int i=0;i<mNbWXCrystObj;i++)
      mpWXCrystObj[i]->Show(show);
   //this->CrystUpdate();
   VFN_DEBUG_MESSAGE("WXCrystObjBasicList::Show(bool):End",3)
   return true;
}

WXCrystObjBasic* WXCrystObjBasicList::Get(const unsigned int i)
{
   if(i>= mNbWXCrystObj)
   {
      cout << "WXCrystObjBasicList::Get(i): i out of range !" <<endl;
      throw 0;
   }
   return mpWXCrystObj[i];
}

void WXCrystObjBasicList::CrystUpdate()
{
   VFN_DEBUG_MESSAGE("WXCrystObjBasicList::CrystUpdate()",3)
   for(unsigned int i=0;i<mNbWXCrystObj;i++)
      mpWXCrystObj[i]->CrystUpdate();
}
void WXCrystObjBasicList::UpdateUI()
{
   VFN_DEBUG_MESSAGE("WXCrystObjBasicList::UpdateUI()",3)
   for(unsigned int i=0;i<mNbWXCrystObj;i++)
      mpWXCrystObj[i]->UpdateUI();
}

////////////////////////////////////////////////////////////////////////
//
//    For the automatic validation of user input
//
////////////////////////////////////////////////////////////////////////
/// This pointer records the last wxField in which something was enetered,
/// so that it can be validated when inpu is finished (either when another
/// input has begun in another field, or when an action requires to purge
/// all input
WXField *spLastWXFieldInputNotValidated=0;

void WXCrystValidateAllUserInput()
{
	if(0==spLastWXFieldInputNotValidated) return;
   VFN_DEBUG_ENTRY("WXCrystValidateAllUserInput()...",6)
	static WXField *pField;
	pField=spLastWXFieldInputNotValidated;
	spLastWXFieldInputNotValidated=0;//avoid loops
	pField->ValidateUserInput();
   VFN_DEBUG_EXIT("WXCrystValidateAllUserInput()...",6)
}

////////////////////////////////////////////////////////////////////////
//
//    WXField
//
////////////////////////////////////////////////////////////////////////
WXField::WXField(wxWindow *parent,const string& label,const int id):
WXCrystObjBasic(parent),mId(id)
{
   VFN_DEBUG_MESSAGE("WXField::WXField()",6)
   mpSizer = new wxBoxSizer(wxHORIZONTAL);
   mpLabel=new wxStaticText(this,-1,label.c_str());
   mpLabel->SetEventHandler(this);
   mpSizer->Add(mpLabel,0,wxALIGN_CENTER);
   this->SetSizer(mpSizer);
   this->Layout();
}
bool WXField::Layout()
{
   VFN_DEBUG_MESSAGE("WXField::Layout()",3)
   mpSizer->Layout();
   mpSizer->Fit(this);
   return this->wxWindow::Layout();
}
void WXField::SetLabel(const string& s)
{
   VFN_DEBUG_MESSAGE("WXField::SetLabel()",3)
   mpLabel->SetLabel(s.c_str());
   mpSizer->SetItemMinSize
               (mpLabel, mpLabel->GetSize().GetWidth(),mpLabel->GetSize().GetHeight());
   this->Layout();
}
bool WXField::SetForegroundColour(const wxColour& colour)
{
   VFN_DEBUG_MESSAGE("WXField::SetLabel()",3)
   mpLabel->SetForegroundColour(colour);
   return this->wxWindow::SetForegroundColour(colour);
}

////////////////////////////////////////////////////////////////////////
//
//    WXFieldName
//
////////////////////////////////////////////////////////////////////////
BEGIN_EVENT_TABLE(WXFieldName,wxEvtHandler)
   EVT_TEXT_ENTER(ID_WXFIELD, 		WXFieldName::OnEnter)
   EVT_TEXT(		ID_WXFIELD, 		WXFieldName::OnText)
END_EVENT_TABLE()

WXFieldName::WXFieldName(wxWindow *parent,const string& label, WXCrystObj* owner,
                         const int id,const int hsize, bool isEditable):
WXField(parent,label,id),mpWXObj(owner),mValue(""),mIsSelfUpdating(false)
{
   VFN_DEBUG_MESSAGE("WXFieldName::WXFieldName():End",6)

   if(true==isEditable)
      mpField=new wxTextCtrl(this,ID_WXFIELD,mValue,
                             wxDefaultPosition,wxSize(hsize,-1),wxTE_PROCESS_ENTER,
                             wxTextValidator(wxFILTER_ASCII));
   else
      mpField=new wxTextCtrl(this,ID_WXFIELD,mValue,
                             wxDefaultPosition,wxSize(hsize,-1),wxTE_READONLY,
                             wxTextValidator(wxFILTER_ASCII));

   mpSizer->Add(mpField,0,wxALIGN_CENTER);
   this->Layout();
}

void WXFieldName::OnEnter(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXFieldName::OnEnter()",6)
	WXCrystValidateAllUserInput();
}
void WXFieldName::OnText(wxCommandEvent & WXUNUSED(event))
{
	if(true==mIsSelfUpdating) return;
   VFN_DEBUG_MESSAGE("WXFieldName::OnText():",6)
	if(spLastWXFieldInputNotValidated!=this)
	{
		WXCrystValidateAllUserInput();
		spLastWXFieldInputNotValidated=this;
	}
}

void WXFieldName::SetValue(const string&s)
{
	if(mValue==(wxString)(s.c_str()))return;
   VFN_DEBUG_MESSAGE("WXFieldName::SetValue()",3)
   mValue=s.c_str();
	mNeedUpdateUI=true;
}

const string WXFieldName::GetValue() const
{
   VFN_DEBUG_MESSAGE("WXFieldName::GetValue()"<<mValue<<":"<<mpField->GetValue(),6)
   return mValue.c_str();
}
void WXFieldName::CrystUpdate()
{
   VFN_DEBUG_MESSAGE("WXFieldName::CrystUpdate()",3)
   // The name must be updated by the owner
}
void WXFieldName::UpdateUI()
{
	if(mNeedUpdateUI==false) return;
   VFN_DEBUG_MESSAGE("WXFieldName::UpdateUI()",10)
	mIsSelfUpdating=true;
   mpField->SetValue(mValue);
	mIsSelfUpdating=false;
	mNeedUpdateUI=false;
}
void WXFieldName::Revert()
{
   VFN_DEBUG_MESSAGE("WXFieldName::Revert()",3)
   mValue=mValueOld;
	mNeedUpdateUI=true;
}
void WXFieldName::ValidateUserInput()
{
   VFN_DEBUG_MESSAGE("WXFieldName::ValidateUserInput()",6)
   //:TODO: Check that the object is not busy (input should be frozen)
   mValueOld=mValue;
   wxString s=mpField->GetValue();
   mValue=s.c_str();
   mpWXObj->OnChangeName(mId);
}

////////////////////////////////////////////////////////////////////////
//
//    WXFieldParBase
//
////////////////////////////////////////////////////////////////////////
BEGIN_EVENT_TABLE(WXFieldParBase,wxWindow)
   EVT_TEXT_ENTER(ID_WXFIELD, 		WXFieldParBase::OnEnter)
   EVT_TEXT(		ID_WXFIELD, 		WXFieldParBase::OnText)
END_EVENT_TABLE()

WXFieldParBase::WXFieldParBase(wxWindow *parent,const string& label,
                               const int id, const int hsize):
WXField(parent,label,id),mIsSelfUpdating(false)
{
   VFN_DEBUG_MESSAGE("WXFieldParBase::WXFieldName():End",6)

   mpField=new wxTextCtrl(this,ID_WXFIELD,"",
                            wxDefaultPosition,wxSize(hsize,-1),wxTE_PROCESS_ENTER,
                            wxTextValidator(wxFILTER_NUMERIC));
   mpSizer->Add(mpField,0,wxALIGN_CENTER);
   
   this->Layout();
}


void WXFieldParBase::OnEnter(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXFieldRefPar::OnEnter()",6)
   WXCrystValidateAllUserInput();
}
void WXFieldParBase::OnText(wxCommandEvent & WXUNUSED(event))
{	
	if(true==mIsSelfUpdating) return;
   VFN_DEBUG_MESSAGE("WXFieldRefPar::OnText()",6)
	if(spLastWXFieldInputNotValidated!=this)
	{
		WXCrystValidateAllUserInput();
		spLastWXFieldInputNotValidated=this;
	}
}
void WXFieldParBase::ValidateUserInput()
{
   VFN_DEBUG_MESSAGE("WXFieldRefPar::ValidateUserInput()",6)
	this->ReadNewValue();
}

////////////////////////////////////////////////////////////////////////
//
//    WXFieldPar<T>
//
////////////////////////////////////////////////////////////////////////
template<class T> WXFieldPar<T>::WXFieldPar(wxWindow *parent,const string& label, 
                                            const int id,T *par,const int hsize):
WXFieldParBase(parent,label,id,hsize),mpValue(par),mValue(*par),mValueOld(*par)
{
   this->CrystUpdate();
}

template<class T> void WXFieldPar<T>::CrystUpdate()
{
   VFN_DEBUG_MESSAGE("WXFieldPar<T>::CrystUpdate()",6)
	mValueOld=mValue;
	mValue=*mpValue;
	mNeedUpdateUI=true;
}

template<> void WXFieldPar<REAL>::UpdateUI()
{
	if(mNeedUpdateUI==false) return;
   VFN_DEBUG_MESSAGE("WXFieldPar<REAL>::UpdateUI()",6)
   wxString tmp;
   tmp.Printf("%f",mValue);
	mIsSelfUpdating=true;
   mpField->SetValue(tmp);
	mIsSelfUpdating=false;
	mNeedUpdateUI=false;
}

template<> void WXFieldPar<long>::UpdateUI()
{
	if(mNeedUpdateUI==false) return;
   VFN_DEBUG_MESSAGE("WXFieldPar<long>::UpdateUI()",6)
   wxString tmp;
   tmp.Printf("%d",mValue);
	mIsSelfUpdating=true;
   mpField->SetValue(tmp);
	mIsSelfUpdating=false;
	mNeedUpdateUI=false;
}
/*
template<class T> void WXFieldPar<T>::UpdateUI()
{
	if(mNeedUpdateUI==false) return;
   stringstream s;
   s <<*mpValue;
	mIsSelfUpdating=true;
   mpField->SetValue(s.str().c_str());;
   mpField->SetValue(wxString::Printf("%f",mValue));
	mIsSelfUpdating=false;
	mNeedUpdateUI=false;
}
*/

template<class T> void WXFieldPar<T>::Revert()
{
   VFN_DEBUG_MESSAGE("WXFieldPar<T>::Revert()",6)
   *mpValue=mValueOld;
	mValue=mValueOld;
	mNeedUpdateUI=true;
}

template<> void WXFieldPar<REAL>::ReadNewValue()
{
   VFN_DEBUG_MESSAGE("WXFieldPar<REAL>::ReadNewValue()",6)
   mValueOld=*mpValue;
   wxString s=mpField->GetValue();
	double tmp;
   s.ToDouble(&tmp);
	*mpValue=tmp;
}
template<> void WXFieldPar<long>::ReadNewValue()
{
   VFN_DEBUG_MESSAGE("WXFieldPar<long>::ReadNewValue()",6)
   mValueOld=*mpValue;
   wxString s=mpField->GetValue();
   s.ToLong(mpValue);
}


template class WXFieldPar<REAL>;
template class WXFieldPar<long>;


////////////////////////////////////////////////////////////////////////
//
//    WXFieldChoice
//
////////////////////////////////////////////////////////////////////////
WXFieldChoice::WXFieldChoice
   (wxWindow *parent,const int field_id,const string &name,const int hsize):
WXField(parent,name,field_id)
{
   mpButton=new wxButton(this,field_id,name.c_str(),wxDefaultPosition,wxSize(hsize,-1));
   mpSizer->Add(mpButton,0,wxALIGN_CENTER);
   this->Layout();
}


bool WXFieldChoice::Layout()
{
   mpSizer->SetItemMinSize
               (mpButton,mpButton->GetSize().GetWidth(),mpButton->GetSize().GetHeight());
   mpSizer->Layout();
   mpSizer->Fit(this);
   return this->wxWindow::Layout();
}

void WXFieldChoice::CrystUpdate()
{
   //Nothing to do. This should be done by the owner
}
void WXFieldChoice::UpdateUI()
{
   //Nothing to do. This should be done by the owner
}

void WXFieldChoice::Revert()
{
   //Nothing to do. This should be done by the owner
}
void WXFieldChoice::SetValue(const string&name)
{
   mpButton->SetLabel(name.c_str());
}
void WXFieldChoice::ValidateUserInput(){}

////////////////////////////////////////////////////////////////////////
//
//    WXCrystObj
//
////////////////////////////////////////////////////////////////////////
BEGIN_EVENT_TABLE(WXCrystObj,wxEvtHandler)
   EVT_BUTTON(ID_WXOBJ_COLLAPSE,WXCrystObj::OnToggleCollapse)
END_EVENT_TABLE()

WXCrystObj::WXCrystObj(wxWindow* parent,int orient,bool showName):
WXCrystObjBasic(parent),mIsExpanded(true)
{
   VFN_DEBUG_MESSAGE("WXCrystObj::WXCrystObj()",6)
   mpTopSizer= new wxBoxSizer(orient);
   this->SetSizer(mpTopSizer);
   
   mpCollapseButton=new wxButton(this,ID_WXOBJ_COLLAPSE,"-",
                                 wxDefaultPosition,wxSize(14,14));
   mpTopSizer->Add(mpCollapseButton,0, wxALIGN_TOP);//wxRIGHT | wxTOP | wxALIGN_TOP,4
   //mpCollapseButton->PushEventHandler(this);
   
   mpSizer=new wxBoxSizer(wxVERTICAL);
   mpTopSizer->Add(mpSizer,0, wxALIGN_TOP);
   
   //wxStaticText* test=new wxStaticText(mWXParent,-1,"Test!!!");
   //mpSizer->Add(test,0,wxALIGN_LEFT);
   if(showName)
   {
      mpWXTitle = new WXFieldName(this,"name:",this,ID_WXOBJ_NAME,200);
      mpSizer->Add(mpWXTitle,0,wxALIGN_LEFT);
   }else mpWXTitle=0;
   
   this->Layout();
   VFN_DEBUG_MESSAGE("WXCrystObj::WXCrystObj():End",6)
}

WXCrystObj::~WXCrystObj(){}

bool WXCrystObj::Layout()
{
   VFN_DEBUG_MESSAGE("WXCrystObj::Layout()",3)
   //There are probably easier ways to do this
   if(0!=mpWXTitle)
   {
      mpWXTitle->Layout();
      mpSizer->SetItemMinSize
               (mpWXTitle, mpWXTitle->GetSize().GetWidth(),mpWXTitle->GetSize().GetHeight());
   }
   for(unsigned int i=0;i<mList.GetNb();i++)
      mpSizer->SetItemMinSize(mList.Get(i),
                              mList.Get(i)->GetSize().GetWidth(),
                              mList.Get(i)->GetSize().GetHeight());
      
   mpSizer->Layout();
   mpTopSizer->Layout();
   mpTopSizer->Fit(this);
   wxSizer* s=mWXParent->GetSizer();
   if(s != 0)
   {// Need to do it that way, in case  the parent is not a WXCrystObj
    // with an adequate Layout() function
      s->SetItemMinSize(this,this->GetSize().GetWidth(),this->GetSize().GetHeight());
      s->Fit(mWXParent);
   }
   //wxCommandEvent event(1758,-1);
   //wxPostEvent(this->GetParent(),event);
	this->UpdateUI();
   mWXParent->Layout();
   return this->wxWindow::Layout();
}

void WXCrystObj::OnToggleCollapse(wxCommandEvent & WXUNUSED(event))
{
   #if 0
   VFN_DEBUG_MESSAGE("WXCrystObj::OnToggleCollapse()",6)
   mIsExpanded = !mIsExpanded;
   mList.Show(mIsExpanded);
   if(true==mIsExpanded) mpCollapseButton->SetLabel("-");
   else mpCollapseButton->SetLabel("+");
   this->Layout();
   VFN_DEBUG_MESSAGE("WXCrystObj::OnToggleCollapse():End",6)
   #endif
}

void WXCrystObj::CrystUpdate()
{
   VFN_DEBUG_MESSAGE("WXCrystObj::CrystUpdate()",6)
   mList.CrystUpdate();
}
void WXCrystObj::UpdateUI()
{
   VFN_DEBUG_MESSAGE("WXCrystObj::UpdateUI()",6)
	if(mpWXTitle!=0) mpWXTitle->UpdateUI();
   mList.UpdateUI();
}
////////////////////////////////////////////////////////////////////////
//
//    WXCrystMenuBar
//
////////////////////////////////////////////////////////////////////////
BEGIN_EVENT_TABLE(WXCrystMenuBar,wxWindow)
   EVT_BUTTON(ID_CRYST_MENU1,WXCrystMenuBar::OnPopupMenu)
   EVT_BUTTON(ID_CRYST_MENU1+1 ,WXCrystMenuBar::OnPopupMenu)
   EVT_BUTTON(ID_CRYST_MENU1+2 ,WXCrystMenuBar::OnPopupMenu)
   EVT_BUTTON(ID_CRYST_MENU1+3 ,WXCrystMenuBar::OnPopupMenu)
   EVT_BUTTON(ID_CRYST_MENU1+4 ,WXCrystMenuBar::OnPopupMenu)
   EVT_BUTTON(ID_CRYST_MENU1+5 ,WXCrystMenuBar::OnPopupMenu)
   EVT_BUTTON(ID_CRYST_MENU1+6 ,WXCrystMenuBar::OnPopupMenu)
   EVT_BUTTON(ID_CRYST_MENU1+7 ,WXCrystMenuBar::OnPopupMenu)
   EVT_BUTTON(ID_CRYST_MENU1+8 ,WXCrystMenuBar::OnPopupMenu)
END_EVENT_TABLE()

WXCrystMenuBar::WXCrystMenuBar(wxWindow *parent, WXCrystObj* owner):
WXCrystObjBasic(parent),mNbMenu(0),mMaxNbMenu(16),mpMenu(0),mpButton(0)
{
   VFN_DEBUG_MESSAGE("WXCrystMenuBar::WXCrystMenuBar():",6)
   mMenuId.resize(mMaxNbMenu);
   mpSizer=new wxBoxSizer(wxHORIZONTAL);
   mpMenu=new wxMenu*[16];//:TODO:
   mpButton=new wxButton*[16];
   this->SetSizer(mpSizer);
   this->Layout();
}

bool WXCrystMenuBar::Layout()
{
   VFN_DEBUG_MESSAGE("WXCrystMenuBar::Layout()",3)
   mpSizer->Layout();
   mpSizer->Fit(this);
   return this->wxWindow::Layout();
}

void WXCrystMenuBar::AddMenu(const string &name,const int menuId, const string& help)
{
   VFN_DEBUG_MESSAGE("WXCrystMenuBar::AddMenu()",6)
   mpMenu[mNbMenu]= new wxMenu(name.c_str());
   VFN_DEBUG_MESSAGE("WXCrystMenuBar::AddMenu():1",6)
   mpButton[mNbMenu]= new wxButton(this,ID_CRYST_MENU1+mNbMenu,name.c_str());
   VFN_DEBUG_MESSAGE("WXCrystMenuBar::AddMenu():2",6)
   mpSizer->Add(mpButton[mNbMenu],0);
   VFN_DEBUG_MESSAGE("WXCrystMenuBar::AddMenu():3",6)
   mMenuId(mNbMenu++)=menuId;
   this->Layout();
   VFN_DEBUG_MESSAGE("WXCrystMenuBar::AddMenu():End",6)
}

void WXCrystMenuBar::AddMenuItem(const int menuId, int id, const string& item, const string& help,const bool checkable)
{
   VFN_DEBUG_MESSAGE("WXCrystMenuBar::AddMenuItem():",6)
   //:TODO: Check we found the correct menu
   unsigned int i;
   for(i=0;i<mNbMenu;i++) if(menuId==mMenuId(i)) break;
   mpMenu[i]->Append(id,item.c_str(),help.c_str(),checkable);
}

void WXCrystMenuBar::AddMenuItem(const int menuId,int id, const wxString&  item,
                 wxMenu *subMenu, const wxString& helpString)
{
   VFN_DEBUG_MESSAGE("WXCrystMenuBar::AddMenuItem():",6)
   //:TODO: Check we found the correct menu
   unsigned int i;
   for(i=0;i<mNbMenu;i++) if(menuId==mMenuId(i)) break;
   mpMenu[i]->Append(id,item.c_str(),subMenu,helpString.c_str());
}

void WXCrystMenuBar::CrystUpdate()
{
}
void WXCrystMenuBar::UpdateUI()
{
}

void WXCrystMenuBar::OnPopupMenu(wxCommandEvent & event)
{
   VFN_DEBUG_MESSAGE("WXCrystMenuBar::OnPopupMenu():",6)
   const int i=event.GetId()-ID_CRYST_MENU1;
   this->PopupMenu(mpMenu[i],mpButton[i]->GetPosition());
}

}//namespace
