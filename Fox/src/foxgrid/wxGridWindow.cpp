
#ifdef __WX__CRYST__
   #include "ObjCryst/wxCryst/wxCrystal.h"

   //#if defined(__WXGTK__) || defined(__WXMOTIF__) || defined(__WXMAC__) || defined(__WXMGL__) || defined(__WXX11__)
   //   #include "Fox.xpm"
   //#endif
#endif

using namespace ObjCryst;
using namespace std;

#include "wxGridWindow.h"

WXGrigWindow::WXGrigWindow(wxWindow *parent):
wxWindow(parent,-1)
{
   dataLoaded = false;
   m_WXFoxMaster = NULL;
   m_WXFoxSlave = NULL;
}
WXGrigWindow::~WXGrigWindow(void)
{
    if (m_WXFoxMaster != NULL) {
        m_WXFoxMaster->Clear();
        delete m_WXFoxMaster;
    }
    if (m_WXFoxSlave != NULL) {
        m_WXFoxSlave->Clear();
        delete m_WXFoxSlave;
    }
}
WXFoxMaster *WXGrigWindow::StartServer()
{
   //if client or server exist return NULL
   if(m_WXFoxSlave!=NULL) return NULL;
   if(m_WXFoxMaster!=NULL) return NULL;
   m_WXFoxMaster = new WXFoxMaster(this, m_working_dir);
   m_WXFoxMaster->m_dataLoaded = dataLoaded;
   if(this->GetSizer()==NULL) this->SetSizer(new wxBoxSizer(wxVERTICAL));
   this->GetSizer()->Add(m_WXFoxMaster);
   this->Layout();
   wxTopLevelWindow *ptopwin = dynamic_cast<wxTopLevelWindow*>( wxTheApp->GetTopWindow());
   if(ptopwin!=NULL)
   {
      wxString title=ptopwin->GetTitle();
      if(wxNOT_FOUND==title.Find("GRID"))
      {
         title +=" [GRID SERVER]";
         ptopwin->SetTitle(title);
      }
   }
   wxTheApp->GetTopWindow()->Layout();
   wxTheApp->GetTopWindow()->SendSizeEvent();
   return m_WXFoxMaster;
}

WXFoxSlave *WXGrigWindow::StartClientWindow()
{
   //if client or server exist return NULL
   if(m_WXFoxSlave!=NULL) return NULL;
   if(m_WXFoxMaster!=NULL) return NULL;
   m_WXFoxSlave = new WXFoxSlave(this, m_working_dir);
   if(this->GetSizer()==NULL) this->SetSizer(new wxBoxSizer(wxVERTICAL));
   this->GetSizer()->Add(m_WXFoxSlave);
   this->Layout();
   wxTopLevelWindow *ptopwin = dynamic_cast<wxTopLevelWindow*>( wxTheApp->GetTopWindow());
   if(ptopwin!=NULL)
   {
      wxString title=ptopwin->GetTitle();
      if(wxNOT_FOUND==title.Find("GRID"))
      {
         title +=" [GRID CLIENT]";
         ptopwin->SetTitle(title);
      }
   }
   wxTheApp->GetTopWindow()->Layout();
   wxTheApp->GetTopWindow()->SendSizeEvent();
   return m_WXFoxSlave;
}
void WXGrigWindow::DataLoaded()
{
   dataLoaded = true;

   if(m_WXFoxMaster!=NULL) {
      m_WXFoxMaster->m_dataLoaded = true;
   }
}
void WXGrigWindow::Clear()
{
   if(m_WXFoxMaster!=NULL) {
      m_WXFoxMaster->Clear();
      delete m_WXFoxMaster;
      m_WXFoxMaster=NULL;
   }

   if(m_WXFoxSlave!=NULL) {
      m_WXFoxSlave->Clear();
      delete m_WXFoxSlave;
      m_WXFoxSlave=NULL;
   }
}
