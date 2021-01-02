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

#include "ObjCryst/wxCryst/wxCryst.h"
#include "ObjCryst/ObjCryst/Undo.h"

//#include "ObjCryst/Quirks/VFNStreamFormat.h"
#include "ObjCryst/Quirks/VFNDebug.h"

//#include "ObjCryst/RefinableObj/GlobalOptimObj.h"

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
std::map<wxWindowID,std::pair<wxPoint,wxSize> > gvWindowPosition;

wxMultiChoiceDialog_ListBox::wxMultiChoiceDialog_ListBox(wxWindow* parent, const wxString& message, const wxString& caption,
                            int n, const wxString* choices):
wxDialog(parent,-1,_T("Choose the molecule's atoms"),wxDefaultPosition,wxSize(300,400),wxCAPTION|wxSTAY_ON_TOP),
mListBox(this,-1,wxDefaultPosition,wxSize(250,350),n,choices,wxLB_MULTIPLE)
{
   wxBoxSizer *sizer=new wxBoxSizer(wxVERTICAL);
   sizer->Add(&mListBox);
   sizer->Add(this->CreateSeparatedButtonSizer(wxOK | wxCANCEL));
   this->SetSizer(sizer);
   sizer->SetSizeHints(this);
   sizer->Fit(this);
}

wxArrayInt wxMultiChoiceDialog_ListBox::GetSelections() const
{
   wxArrayInt choice;
   mListBox.GetSelections(choice);
   return choice;
}

////////////////////////////////////////////////////////////////////////
//
// Unique ID for menus incrementer
//
////////////////////////////////////////////////////////////////////////
WXCRYST_ID::WXCRYST_ID(){mIndex=mCounter++;}
WXCRYST_ID::operator long(){return mIndex;}
long WXCRYST_ID::mCounter=wxID_HIGHEST+100;

#ifdef VFN_CRYST_MUTEX
////////////////////////////////////////////////////////////////////////
//
// Mutex
//
////////////////////////////////////////////////////////////////////////
CrystMutex::CrystMutex():mNbLock(0){}
CrystMutex::~CrystMutex()
{
   cout <<"~CrystMutex("<<this<<"), Total number of Locks="<<mNbLock<<endl;
}
wxMutexError CrystMutex::Lock()
{
   #if 1
   return this->wxMutex::Lock();
   #else
   cout <<"Thread:"<<wxThread::IsMain()<<":CrystMutex("<<this<<")::Lock() ?"<<endl;
   wxMutexError res=this->wxMutex::TryLock();
   if(res==wxMUTEX_NO_ERROR)
      cout <<"Thread:"<<wxThread::IsMain()<<":CrystMutex("<<this<<")::Lock()-OK"<<endl;
   if(res==wxMUTEX_DEAD_LOCK)
   {
      cout <<"Thread:"<<wxThread::IsMain()<<":CrystMutex("<<this<<")::Lock()-DEADLOCK!"<<endl;
      exit(0);
   }
   if(res==wxMUTEX_BUSY)
   {
      cout <<"Thread:"<<wxThread::IsMain()<<":CrystMutex("<<this<<")::Lock()-Busy?"<<endl;
      res=this->wxMutex::Lock();
      cout <<"Thread:"<<wxThread::IsMain()<<":CrystMutex("<<this<<")::Lock()-OK2"<<endl;
   }
   return res;
   #endif
}
wxMutexError CrystMutex::Unlock()
{
   //cout <<"Thread:"<<wxThread::IsMain()<<":CrystMutex("<<this<<"::Unlock()"<<endl;
   mNbLock++;
   return this->wxMutex::Unlock();
}
#endif

////////////////////////////////////////////////////////////////////////
//
//    WXCrystObjBasic
//
////////////////////////////////////////////////////////////////////////
WXCrystObjBasic::WXCrystObjBasic(wxWindow* parent):
wxWindow(parent,-1),mWXCrystParent(0),mIsShown(true),mNeedUpdateUI(true)
{
   VFN_DEBUG_MESSAGE("WXCrystObjBasic::WXCrystObjBasic() at "<<this,6)
   if(parent !=0) mWXCrystParent=dynamic_cast<WXCrystObjBasic*>(parent);
#ifdef __WXGTK__
#ifndef __WXGTK20__
   // Why is this necessary for GTK 1.2 ? wxWidgets 2.6.3
   this->SetBackgroundColour(wxColour(240,240,240));// wxLIGHT_GREY
#endif
#endif
   VFN_DEBUG_MESSAGE("WXCrystObjBasic::WXCrystObjBasic():End",6)
}

WXCrystObjBasic::~WXCrystObjBasic()
{
   // Every time we destroy a widget, validate all input to make sure the destroyed
   // widget does not have some unread info.
   WXCrystValidateAllUserInput();
   set<WXCrystObjBasicList*> vpList=mvpList;//use a copy
   for(set<WXCrystObjBasicList*>::iterator pos=vpList.begin();pos!=vpList.end();++pos)
      (*pos)->Remove(this);
}

void WXCrystObjBasic::BottomLayout(WXCrystObjBasic *pChild)
{
   VFN_DEBUG_ENTRY("WXCrystObjBasic::BottomLayout(...)"<<this->GetSize().GetWidth()<<","<<this->GetSize().GetHeight(),5);
   wxTheApp->GetTopWindow()->SendSizeEvent();
   VFN_DEBUG_EXIT("WXCrystObjBasic::BottomLayout(...)"<<this->GetSize().GetWidth()<<","<<this->GetSize().GetHeight(),5);
}
void WXCrystObjBasic::AddChild(WXCrystObjBasic *pChild, bool doBottomLayout)
{
   VFN_DEBUG_ENTRY("WXCrystObjBasic::AddChild(...)"<<this->GetSize().GetWidth()<<","<<this->GetSize().GetHeight(),5);
   wxSizer *pSizer=this->GetSizer();
   if(pSizer!=0)
   {
      pSizer->Add(pChild);
   }
   if(doBottomLayout) wxTheApp->GetTopWindow()->SendSizeEvent();
   VFN_DEBUG_EXIT("WXCrystObjBasic::AddChild(...)"<<this->GetSize().GetWidth()<<","<<this->GetSize().GetHeight(),5);
}

void WXCrystObjBasic::AddedToList(WXCrystObjBasicList* list){mvpList.insert(list);}
void WXCrystObjBasic::RemovedFromList(WXCrystObjBasicList* list){mvpList.erase(list);}
bool WXCrystObjBasic::Layout()
{
   //static unsigned long ct=0;
   //cout<<"WXCrystObjBasic::Layout() #"<<++ct<<endl;
   return wxWindow::Layout();
}
void WXCrystObjBasic::SetToolTip(const wxString& tip){this->wxWindow::SetToolTip(tip);}

////////////////////////////////////////////////////////////////////////
//
//    WXCrystObjBasicList
//
////////////////////////////////////////////////////////////////////////
WXCrystObjBasicList::WXCrystObjBasicList()
{}

WXCrystObjBasicList::~WXCrystObjBasicList()
{
   set<WXCrystObjBasic*> vpWXCrystObj=mvpWXCrystObj;
   for(set<WXCrystObjBasic*>::iterator pos=vpWXCrystObj.begin();pos!=vpWXCrystObj.end();pos++)
      (*pos)->RemovedFromList(this);
   mvpWXCrystObj.clear();
}

unsigned int WXCrystObjBasicList::GetNb()const {return mvpWXCrystObj.size();}

void WXCrystObjBasicList::Add(WXCrystObjBasic *win)
{
   VFN_DEBUG_MESSAGE("WXCrystObjBasicList::Add()",6)
   win->AddedToList(this);
   mvpWXCrystObj.insert(win);
}

void WXCrystObjBasicList::Remove(WXCrystObjBasic *win)
{
   VFN_DEBUG_MESSAGE("WXCrystObjBasicList::Remove():"<<win,6)
   mvpWXCrystObj.erase(win);
   win->RemovedFromList(this);
}

bool WXCrystObjBasicList::Show(bool show)
{
   VFN_DEBUG_MESSAGE("WXCrystObjBasicList::Show(bool)",3)
   for(set<WXCrystObjBasic*>::iterator pos=mvpWXCrystObj.begin();pos!=mvpWXCrystObj.end();pos++)
      (*pos)->Show(show);
   //this->CrystUpdate();
   VFN_DEBUG_MESSAGE("WXCrystObjBasicList::Show(bool):End",3)
   return true;
}

void WXCrystObjBasicList::CrystUpdate(const bool updateUI,const bool mutexlock)
{
   VFN_DEBUG_ENTRY("WXCrystObjBasicList::CrystUpdate()",5)
   for(set<WXCrystObjBasic*>::iterator pos=mvpWXCrystObj.begin();pos!=mvpWXCrystObj.end();pos++)
   {
      VFN_DEBUG_MESSAGE("WXCrystObjBasicList::CrystUpdate("<<updateUI<<mutexlock<<")@"<<*pos,5)
      (*pos)->CrystUpdate(updateUI,mutexlock);
   }
   VFN_DEBUG_EXIT("WXCrystObjBasicList::CrystUpdate()",5)
}
void WXCrystObjBasicList::UpdateUI(const bool mutexlock)
{
   VFN_DEBUG_ENTRY("WXCrystObjBasicList::UpdateUI("<<mutexlock<<")"<<"MainThread="<<wxThread::IsMain(),5)
   for(set<WXCrystObjBasic*>::iterator pos=mvpWXCrystObj.begin();pos!=mvpWXCrystObj.end();pos++)
      (*pos)->UpdateUI(mutexlock);
   VFN_DEBUG_EXIT("WXCrystObjBasicList::UpdateUI()",5)
}
void WXCrystObjBasicList::Enable(bool enable)
{
   for(set<WXCrystObjBasic*>::iterator pos=mvpWXCrystObj.begin();pos!=mvpWXCrystObj.end();pos++)
      (*pos)->Enable(enable);
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
   /*
   // Hijack this function to record changes in the configuration
   gConfigHistory.Store();
    */
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
   VFN_DEBUG_ENTRY("WXField::WXField()",6)
   mpSizer = new wxBoxSizer(wxHORIZONTAL);
   mpLabel=new wxStaticText(this,-1,wxString::FromAscii(label.c_str()));
   mpSizer->Add(mpLabel,0,wxALIGN_CENTER);
   this->SetSizer(mpSizer);
   //mpSizer->SetSizeHints(this);
   //this->Layout();
   VFN_DEBUG_EXIT("WXField::WXField()",6)
}
void WXField::SetLabel(const string& s)
{
   VFN_DEBUG_MESSAGE("WXField::SetLabel()",3)
   mpLabel->SetLabel(wxString::FromAscii(s.c_str()));
   mpLabel->Layout();
}
bool WXField::SetForegroundColour(const wxColour& colour)
{
   VFN_DEBUG_MESSAGE("WXField::SetLabel()",3)
   mpLabel->SetForegroundColour(colour);
   return this->wxWindow::SetForegroundColour(colour);
}

void WXField::SetSize(int width, int height)
{}
////////////////////////////////////////////////////////////////////////
//
//    WXFieldString
//
////////////////////////////////////////////////////////////////////////
BEGIN_EVENT_TABLE(WXFieldString,wxWindow)
   EVT_TEXT_ENTER(ID_WXFIELD,       WXFieldString::OnEnter)
   EVT_TEXT(      ID_WXFIELD,       WXFieldString::OnText)
END_EVENT_TABLE()

WXFieldString::WXFieldString(wxWindow *parent,string& st,
                         const int id,const int hsize, bool isEditable):
WXField(parent,"",id),mpString(&st),mValue(st),mIsSelfUpdating(false)
{
   VFN_DEBUG_MESSAGE("WXFieldString::WXFieldName():End",6)

   if(true==isEditable)
      mpField=new wxTextCtrl(this,ID_WXFIELD,wxString::FromAscii(mValue.c_str()),
                             wxDefaultPosition,wxSize(hsize,-1),wxTE_PROCESS_ENTER,
                             wxTextValidator(wxFILTER_ASCII));
   else
      mpField=new wxTextCtrl(this,ID_WXFIELD,wxString::FromAscii(mValue.c_str()),
                             wxDefaultPosition,wxSize(hsize,-1),wxTE_READONLY,
                             wxTextValidator(wxFILTER_ASCII));

   mpSizer->Add(mpField,0,wxALIGN_CENTER);
   mpSizer->SetSizeHints(this);
   this->Layout();
}

void WXFieldString::OnEnter(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXFieldString::OnEnter()",6)
   WXCrystValidateAllUserInput();
}
void WXFieldString::OnText(wxCommandEvent & WXUNUSED(event))
{
   if(true==mIsSelfUpdating) return;
   VFN_DEBUG_MESSAGE("WXFieldString::OnText():",6)
   if(spLastWXFieldInputNotValidated!=this)
   {
      WXCrystValidateAllUserInput();
      spLastWXFieldInputNotValidated=this;
   }
}

void WXFieldString::SetValue(const string&s)
{
   VFN_DEBUG_ENTRY("WXFieldString::SetValue()",3)
   wxMutexLocker mlock(mMutex);
   if(mValue==s)
   {
      VFN_DEBUG_EXIT("WXFieldString::SetValue(): string unchanged",3)
      return;
   }
   mValue=s;
   mNeedUpdateUI=true;
   VFN_DEBUG_EXIT("WXFieldString::SetValue()",3)
}

const string WXFieldString::GetValue() const
{
   VFN_DEBUG_MESSAGE("WXFieldString::GetValue()"<<mValue<<":"<<mpField->GetValue(),6)
   return mValue;
}
void WXFieldString::CrystUpdate(const bool uui,const bool lock)
{
   VFN_DEBUG_ENTRY("WXFieldString::CrystUpdate()",3)
   if(lock) mMutex.Lock();
   if(mValue==*mpString)
   {
      if(lock)mMutex.Unlock();
      return;
   }
   mValueOld=mValue;
   mValue=*mpString;
   mNeedUpdateUI=true;
   if(lock) mMutex.Unlock();
   if(uui) if(true==wxThread::IsMain()) this->UpdateUI(lock);
   VFN_DEBUG_EXIT("WXFieldString::CrystUpdate()",3)
}
void WXFieldString::UpdateUI(const bool lock)
{
   if(lock) mMutex.Lock();
   if(mNeedUpdateUI==false)
   {
      if(lock) mMutex.Unlock();
      return;
   }
   VFN_DEBUG_ENTRY("WXFieldString::UpdateUI("<<lock<<")"<<"MainThread="<<wxThread::IsMain(),4)
   mIsSelfUpdating=true;
   mpField->SetValue(wxString::FromAscii(mValue.c_str()));
   mIsSelfUpdating=false;
   mNeedUpdateUI=false;
   if(lock) mMutex.Unlock();
   VFN_DEBUG_EXIT("WXFieldString::UpdateUI()",4)
}
void WXFieldString::Revert()
{
   VFN_DEBUG_MESSAGE("WXFieldString::Revert()",3)
   wxMutexLocker mlock(mMutex);
   mValue=mValueOld;
   mNeedUpdateUI=true;
}
void WXFieldString::ValidateUserInput()
{
   VFN_DEBUG_MESSAGE("WXFieldString::ValidateUserInput()",6)
   //:TODO: Check that the object is not busy (input should be frozen)
   wxMutexLocker mlock(mMutex);
   mValueOld=mValue;
   mValue=mpField->GetValue().ToAscii();
   *mpString=mValue;
}
void WXFieldString::SetSize(int width, int height)
{
   mpField->SetSize(width,height);
   //wxTheApp->GetTopWindow()->SendSizeEvent();
}

void WXFieldString::SetToolTip(const wxString& tip){mpField->SetToolTip(tip);}

////////////////////////////////////////////////////////////////////////
//
//    WXFieldName
//
////////////////////////////////////////////////////////////////////////
BEGIN_EVENT_TABLE(WXFieldName,wxWindow)
   EVT_TEXT_ENTER(ID_WXFIELD,       WXFieldName::OnEnter)
   EVT_TEXT(      ID_WXFIELD,       WXFieldName::OnText)
END_EVENT_TABLE()

WXFieldName::WXFieldName(wxWindow *parent,const string& label, WXCrystObj* owner,
                         const int id,const int hsize, bool isEditable):
WXField(parent,label,id),mpWXObj(owner),mValue(""),mIsSelfUpdating(false)
{
   VFN_DEBUG_ENTRY("WXFieldName::WXFieldName()",6)
   if(true==isEditable)
      mpField=new wxTextCtrl(this,ID_WXFIELD,wxString::FromAscii(mValue.c_str()),
                             wxDefaultPosition,wxSize(hsize,-1),wxTE_PROCESS_ENTER,
                             wxTextValidator(wxFILTER_ASCII));
   else
      mpField=new wxTextCtrl(this,ID_WXFIELD,wxString::FromAscii(mValue.c_str()),
                             wxDefaultPosition,wxSize(hsize,-1),wxTE_READONLY,
                             wxTextValidator(wxFILTER_ASCII));

   mpSizer->Add(mpField,0,wxALIGN_CENTER);
   //mpSizer->SetSizeHints(this);
   this->Layout();
   VFN_DEBUG_EXIT("WXFieldName::WXFieldName()",6)
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
   VFN_DEBUG_ENTRY("WXFieldName::SetValue()",3)
   wxMutexLocker mlock(mMutex);
   if(mValue==s)
   {
      VFN_DEBUG_EXIT("WXFieldName::SetValue():name unchanged",3)
      return;
   }
   mValue=s;
   mNeedUpdateUI=true;
   VFN_DEBUG_EXIT("WXFieldName::SetValue()",3)
}

const string WXFieldName::GetValue() const
{
   VFN_DEBUG_MESSAGE("WXFieldName::GetValue()"<<mValue<<":"<<mpField->GetValue(),6)
   return mValue;
}
void WXFieldName::CrystUpdate(const bool updateUI,const bool lock)
{
   VFN_DEBUG_MESSAGE("WXFieldName::CrystUpdate()",3)
   // The name must be updated by the owner
}
void WXFieldName::UpdateUI(const bool lock)
{
   if(lock) mMutex.Lock();
   if(mNeedUpdateUI==false)
   {
      if(lock) mMutex.Unlock();
      return;
   }
   VFN_DEBUG_ENTRY("WXFieldName::UpdateUI("<<lock<<")"<<"MainThread="<<wxThread::IsMain(),4)
   mIsSelfUpdating=true;
   mpField->SetValue(wxString::FromAscii(mValue.c_str()));
   //:TODO: Some way of resizing the field which is less of a kludge...
   int w=mpField->GetTextExtent(wxString::FromAscii(mValue.c_str())).GetWidth();
   const int wmax=wxTheApp->GetTopWindow()->GetSize().GetWidth();
   if(w>wmax)w=wmax;
   if(w>mpField->GetSize().GetWidth())
      this->GetSizer()->SetItemMinSize(mpField,w+30,-1);
   this->GetSizer()->Fit(this);

   mIsSelfUpdating=false;
   mNeedUpdateUI=false;
   if(lock) mMutex.Unlock();
   VFN_DEBUG_EXIT("WXFieldName::UpdateUI()",4)
}
void WXFieldName::Revert()
{
   VFN_DEBUG_MESSAGE("WXFieldName::Revert()",3)
   wxMutexLocker mlock(mMutex);
   mValue=mValueOld;
   mNeedUpdateUI=true;
}
void WXFieldName::ValidateUserInput()
{
   VFN_DEBUG_MESSAGE("WXFieldName::ValidateUserInput()",6)
   mMutex.Lock();
   mValueOld=mValue;
   mValue=mpField->GetValue().ToAscii();
   mMutex.Unlock();
   mpWXObj->OnChangeName(mId);
}
void WXFieldName::SetSize(int width, int height)
{//:TODO: deprecate this function ?
   mpField->SetSize(width,height);
}

void WXFieldName::SetToolTip(const wxString& tip){mpField->SetToolTip(tip);}

////////////////////////////////////////////////////////////////////////
//
//    WXFieldParBase
//
////////////////////////////////////////////////////////////////////////
BEGIN_EVENT_TABLE(WXFieldParBase,wxWindow)
   EVT_TEXT_ENTER(ID_WXFIELD,       WXFieldParBase::OnEnter)
   EVT_TEXT(      ID_WXFIELD,       WXFieldParBase::OnText)
END_EVENT_TABLE()

WXFieldParBase::WXFieldParBase(wxWindow *parent,const string& label,
                               const int id, const int hsize):
WXField(parent,label,id),mIsSelfUpdating(false),mFormat(_T("%8f"))
{
   VFN_DEBUG_MESSAGE("WXFieldParBase::WXFieldName():End",6)

   mpField=new wxTextCtrl(this,ID_WXFIELD,_T(""),
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

void WXFieldParBase::SetToolTip(const wxString& tip){mpField->SetToolTip(tip);}

void WXFieldParBase::SetFormat(const wxString &format)
{
   mFormat=format;
}

////////////////////////////////////////////////////////////////////////
//
//    WXFieldPar<T>
//
////////////////////////////////////////////////////////////////////////
template<class T> WXFieldPar<T>::WXFieldPar(wxWindow *parent,const string& label,
                                            const int id,T *par,const int hsize):
WXFieldParBase(parent,label,id,hsize),mpValue(par),mValue(*par),mValueOld(*par),mHumanScale(1)
{
   this->CrystUpdate(true,true);
}

template<> WXFieldPar<long>::WXFieldPar(wxWindow *parent,const string& label,
                                            const int id,long *par,const int hsize):
WXFieldParBase(parent,label,id,hsize),mpValue(par),mValue(*par),mValueOld(*par),mHumanScale(1)
{
   mFormat=_T("%ld");
   this->CrystUpdate(true,true);
}

template<class T> void WXFieldPar<T>::CrystUpdate(const bool uui,const bool lock)
{
   if(lock) mMutex.Lock();
   if(mValue==*mpValue)
   {
      if(lock) mMutex.Unlock();
      return;
   }
   VFN_DEBUG_ENTRY("WXFieldPar<T>::CrystUpdate()",6)
   mValueOld=mValue;
   mValue=*mpValue;
   mNeedUpdateUI=true;
   if(lock) mMutex.Unlock();
   if(uui) if(true==wxThread::IsMain()) this->UpdateUI(lock);
   VFN_DEBUG_EXIT("WXFieldPar<T>::CrystUpdate()",6)
}

template<> void WXFieldPar<REAL>::UpdateUI(const bool lock)
{
   if(lock)mMutex.Lock();
   if(mNeedUpdateUI==false)
   {
      if(lock)mMutex.Unlock();
      return;
   }
   VFN_DEBUG_ENTRY("WXFieldPar<REAL>::UpdateUI("<<lock<<")"<<"MainThread="<<wxThread::IsMain(),4)
   wxString tmp;
   if((abs(mValue*mHumanScale)<1000)&&(abs(mValue*mHumanScale)>0.01)) tmp.Printf(_T("%6.4f"),mValue*mHumanScale);
   else tmp.Printf(mFormat,mValue*mHumanScale);
   mIsSelfUpdating=true;
   mpField->SetValue(tmp);
   mIsSelfUpdating=false;
   mNeedUpdateUI=false;
   if(lock)mMutex.Unlock();
   VFN_DEBUG_EXIT("WXFieldPar<REAL>::UpdateUI()",4)
}

template<> void WXFieldPar<long>::UpdateUI(const bool lock)
{
   if(lock)mMutex.Lock();
   if(mNeedUpdateUI==false)
   {
      if(lock)mMutex.Unlock();
      return;
   }
   VFN_DEBUG_ENTRY("WXFieldPar<long>::UpdateUI("<<lock<<")"<<"MainThread="<<wxThread::IsMain(),4)
   wxString tmp;
   tmp.Printf(mFormat,mValue*mHumanScale);
   mIsSelfUpdating=true;
   mpField->SetValue(tmp);
   mIsSelfUpdating=false;
   mNeedUpdateUI=false;
   if(lock)mMutex.Unlock();
   VFN_DEBUG_EXIT("WXFieldPar<long>::UpdateUI()",4)
}
/*
template<class T> void WXFieldPar<T>::UpdateUI(const bool lock)
{
   if(lock)mMutex.Lock();
   if(mNeedUpdateUI==false)
   {
      if(lock)mMutex.Unlock();
      return;
   }
   stringstream s;
   s <<*mpValue;
   mIsSelfUpdating=true;
   mpField->SetValue(s.str().c_str());;
   mpField->SetValue(wxString::Printf("%f",mValue));
   mIsSelfUpdating=false;
   mNeedUpdateUI=false;
   if(lock)mMutex.Unlock();
}
*/

template<class T> void WXFieldPar<T>::Revert()
{
   VFN_DEBUG_MESSAGE("WXFieldPar<T>::Revert()",6)
   mMutex.Lock();
   *mpValue=mValueOld;
   mValue=mValueOld;
   mNeedUpdateUI=true;
   mMutex.Unlock();
   if(true==wxThread::IsMain()) this->UpdateUI(true);
}

template<class T> void WXFieldPar<T>::SetHumanValueScale(const T s)
{
   mHumanScale=s;
}

template<> void WXFieldPar<REAL>::ReadNewValue()
{
   VFN_DEBUG_MESSAGE("WXFieldPar<REAL>::ReadNewValue()",6)
   wxMutexLocker mlock(mMutex);
   mValueOld=*mpValue;
   wxString s=mpField->GetValue();
   double tmp;
   s.ToDouble(&tmp);
   *mpValue=tmp;
   *mpValue /= mHumanScale;
}
template<> void WXFieldPar<long>::ReadNewValue()
{
   VFN_DEBUG_MESSAGE("WXFieldPar<long>::ReadNewValue()",6)
   wxMutexLocker mlock(mMutex);
   mValueOld=*mpValue;
   wxString s=mpField->GetValue();
   s.ToLong(mpValue);
   *mpValue /= mHumanScale;
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
   mpButton=new wxButton(this,field_id,wxString::FromAscii(name.c_str()),wxDefaultPosition,wxSize(hsize,-1));
   mpSizer->Add(mpButton,0,wxALIGN_CENTER);
   //mpSizer->SetItemMinSize(mpButton,
   //                        mpButton->GetSize().GetWidth(),
   //                        mpButton->GetSize().GetHeight());
   //mpSizer->SetSizeHints(this);
   this->Layout();
}

void WXFieldChoice::CrystUpdate(const bool uui,const bool lock)
{
   //Nothing to do. This should be done by the owner
}
void WXFieldChoice::UpdateUI(const bool lock)
{
   //Nothing to do. This should be done by the owner
}

void WXFieldChoice::Revert()
{
   //Nothing to do. This should be done by the owner
}
void WXFieldChoice::SetValue(const string&name)
{
   mpButton->SetLabel(wxString::FromAscii(name.c_str()));
}
void WXFieldChoice::ValidateUserInput(){}

////////////////////////////////////////////////////////////////////////
//
//    WXCrystObj
//
////////////////////////////////////////////////////////////////////////
const long ID_WXOBJ_ENABLE=WXCRYST_ID(); //These are used in ObjCryst/RefinableObj.cpp
const long ID_WXOBJ_DISABLE=WXCRYST_ID();
BEGIN_EVENT_TABLE(WXCrystObj,wxWindow)
   EVT_BUTTON(ID_WXOBJ_COLLAPSE,WXCrystObj::OnToggleCollapse)
   EVT_UPDATE_UI(ID_WXOBJ_ENABLE,WXCrystObj::OnEnable)
   EVT_UPDATE_UI(ID_WXOBJ_DISABLE,WXCrystObj::OnEnable)
END_EVENT_TABLE()

WXCrystObj::WXCrystObj(wxWindow* parent,int orient,bool showName):
WXCrystObjBasic(parent),mIsExpanded(true)
{
   VFN_DEBUG_ENTRY("WXCrystObj::WXCrystObj()",6)
   mpTopSizer= new wxBoxSizer(orient);
   this->SetSizer(mpTopSizer);

   mpCollapseButton=new wxButton(this,ID_WXOBJ_COLLAPSE,_T("-"),
                                 wxDefaultPosition,wxSize(14,14));
   mpTopSizer->Add(mpCollapseButton,0, wxALIGN_TOP);//wxRIGHT | wxTOP | wxALIGN_TOP,4
   //mpCollapseButton->PushEventHandler(this);

   mpSizer=new wxBoxSizer(wxVERTICAL);
   mpTopSizer->Add(mpSizer,0, wxALIGN_TOP);

   if(showName)
   {
      mpWXTitle = new WXFieldName(this,"name:",this,ID_WXOBJ_NAME,100);
      mpSizer->Add(mpWXTitle,0,wxALIGN_LEFT);
   }else mpWXTitle=0;

   //this->Layout();
   VFN_DEBUG_EXIT("WXCrystObj::WXCrystObj():End",6)
}

WXCrystObj::~WXCrystObj()
{
   //if(0!=mpWXTitle)
   //{
   //   delete mpWXTitle;
   //}
   //delete mpCollapseButton;
   //sizers are automatically deleted
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

void WXCrystObj::CrystUpdate(const bool uui,const bool lock)
{
   VFN_DEBUG_ENTRY("WXCrystObj::CrystUpdate("<<uui<<lock<<")",6)
   //if(lock) mMutex.Lock();
   mList.CrystUpdate(uui,lock);
   //if(lock) mMutex.Unlock();
   VFN_DEBUG_EXIT("WXCrystObj::CrystUpdate()",6)
}
void WXCrystObj::UpdateUI(const bool lock)
{
   VFN_DEBUG_ENTRY("WXCrystObj::UpdateUI("<<lock<<")"<<"MainThread="<<wxThread::IsMain(),6)
   if(lock) mMutex.Lock();
   if(mpWXTitle!=0) mpWXTitle->UpdateUI(false);
   mList.UpdateUI(false);
   if(lock) mMutex.Unlock();
   VFN_DEBUG_EXIT("WXCrystObj::UpdateUI()",6)
}
void WXCrystObj::OnEnable(wxUpdateUIEvent &event)
{
   if(ID_WXOBJ_ENABLE==event.GetId()) this->Enable(true);
   else this->Enable(false);
   event.Skip();
}
bool WXCrystObj::Enable(bool enable)
{
   mList.Enable(enable);
   return this->wxWindow::Enable(enable);
}

void WXCrystObj::AddChild(WXCrystObjBasic *pChild, bool doBottomLayout)
{
   VFN_DEBUG_ENTRY("WXCrystObj::AddChild(..)"<<this->GetSize().GetWidth()<<","<<this->GetSize().GetHeight(),5);
   if(mpSizer!=0)
   {
      mpSizer->Add(pChild);
   }
   if(doBottomLayout) wxTheApp->GetTopWindow()->SendSizeEvent();
   VFN_DEBUG_EXIT("WXCrystObj::AddChild(..)"<<this->GetSize().GetWidth()<<","<<this->GetSize().GetHeight(),5);
}

////////////////////////////////////////////////////////////////////////
//
//    WXCrystMenuBar
//
////////////////////////////////////////////////////////////////////////
WXCRYST_ID ID_CRYST_MENU1;
WXCRYST_ID ID_CRYST_MENU2;
WXCRYST_ID ID_CRYST_MENU3;
WXCRYST_ID ID_CRYST_MENU4;
WXCRYST_ID ID_CRYST_MENU5;
WXCRYST_ID ID_CRYST_MENU6;
WXCRYST_ID ID_CRYST_MENU7;
WXCRYST_ID ID_CRYST_MENU8;
WXCRYST_ID ID_CRYST_MENU9;
WXCRYST_ID ID_CRYST_MENU10;
WXCRYST_ID ID_CRYST_MENU11;
WXCRYST_ID ID_CRYST_MENU12;
WXCRYST_ID ID_CRYST_MENU13;
WXCRYST_ID ID_CRYST_MENU14;
WXCRYST_ID ID_CRYST_MENU15;
WXCRYST_ID ID_CRYST_MENU16;


BEGIN_EVENT_TABLE(WXCrystMenuBar,wxWindow)
   EVT_BUTTON(ID_CRYST_MENU1 ,WXCrystMenuBar::OnPopupMenu)
   EVT_BUTTON(ID_CRYST_MENU2 ,WXCrystMenuBar::OnPopupMenu)
   EVT_BUTTON(ID_CRYST_MENU3 ,WXCrystMenuBar::OnPopupMenu)
   EVT_BUTTON(ID_CRYST_MENU4 ,WXCrystMenuBar::OnPopupMenu)
   EVT_BUTTON(ID_CRYST_MENU5 ,WXCrystMenuBar::OnPopupMenu)
   EVT_BUTTON(ID_CRYST_MENU6 ,WXCrystMenuBar::OnPopupMenu)
   EVT_BUTTON(ID_CRYST_MENU7 ,WXCrystMenuBar::OnPopupMenu)
   EVT_BUTTON(ID_CRYST_MENU8 ,WXCrystMenuBar::OnPopupMenu)
   EVT_BUTTON(ID_CRYST_MENU9 ,WXCrystMenuBar::OnPopupMenu)
END_EVENT_TABLE()

WXCrystMenuBar::WXCrystMenuBar(wxWindow *parent, WXCrystObj* owner):
WXCrystObjBasic(parent)
{
   VFN_DEBUG_MESSAGE("WXCrystMenuBar::WXCrystMenuBar():",6)
   mpSizer=new wxBoxSizer(wxHORIZONTAL);
   this->SetSizer(mpSizer);
}

void WXCrystMenuBar::AddMenu(const string &name,const int menuId, const string& help)
{
   VFN_DEBUG_MESSAGE("WXCrystMenuBar::AddMenu()",6)
   const long id=ID_CRYST_MENU1+mvpMenu.size();
   mvpMenu[menuId]=make_pair(new wxMenu(wxString::FromAscii(name.c_str())),
                             new wxButton(this,id,wxString::FromAscii(name.c_str())));
   VFN_DEBUG_MESSAGE("WXCrystMenuBar::AddMenu():2",6)
   mvpMenu[menuId].second->Layout();
   mpSizer->Add(mvpMenu[menuId].second,0,wxALIGN_CENTER);
   //wxTheApp->GetTopWindow()->SendSizeEvent();
   VFN_DEBUG_MESSAGE("WXCrystMenuBar::AddMenu():End",6)
}

wxMenu& WXCrystMenuBar::GetMenu(const int menuId)
{
   //:TODO: Check we found the correct menu
   return *(mvpMenu[menuId].first);
}

void WXCrystMenuBar::AddMenuItem(const int menuId, int id, const string& item, const string& help,const bool checkable)
{
   VFN_DEBUG_MESSAGE("WXCrystMenuBar::AddMenuItem():",6)
   //:TODO: Check we found the correct menu
   this->GetMenu(menuId).Append(id,wxString::FromAscii(item.c_str()),wxString::FromAscii(help.c_str()),checkable);
}

void WXCrystMenuBar::AddMenuItem(const int menuId,int id, const wxString&  item,
                 wxMenu *subMenu, const wxString& helpString)
{
   VFN_DEBUG_MESSAGE("WXCrystMenuBar::AddMenuItem():",6)
   //:TODO: Check we found the correct menu
   this->GetMenu(menuId).Append(id,item,subMenu,helpString);
}

void WXCrystMenuBar::CrystUpdate(const bool updateUI,const bool mutexlock)
{
}
void WXCrystMenuBar::UpdateUI(const bool mutexlock)
{
}

void WXCrystMenuBar::OnPopupMenu(wxCommandEvent & event)
{
   VFN_DEBUG_MESSAGE("WXCrystMenuBar::OnPopupMenu():",6)
   // The ID of the button is not the same as the ID of the menu,
   // so look it up. (to have the same ID we'd need to delegate event handling,
   // which would be a pain.
   std::map<long,pair<wxMenu *,wxButton*> >::iterator pos;
   const int i=event.GetId();
   for(pos=mvpMenu.begin();pos!=mvpMenu.end();++pos) if(i==pos->second.second->GetId()) break;
   if(pos!=mvpMenu.end())
      this->PopupMenu(pos->second.first,pos->second.second->GetPosition());
}
void WXCrystMenuBar::SetToolTip(const wxString& tip,long menuId){mvpMenu[menuId].second->SetToolTip(tip);}

}//namespace
