














#ifdef __WX__CRYST__
   #include "ObjCryst/wxCryst/wxCrystal.h"

   //#if defined(__WXGTK__) || defined(__WXMOTIF__) || defined(__WXMAC__) || defined(__WXMGL__) || defined(__WXX11__)
   //   #include "Fox.xpm"
   //#endif
#endif

using namespace ObjCryst;
using namespace std;

#include "WXGridWindow.h"

WXGrigWindow::WXGrigWindow(wxWindow *parent):
wxWindow(parent,-1)
{
   dataLoaded = false;
   m_WXFoxServer = NULL;
   m_WXFoxClient = NULL;
}
WXGrigWindow::~WXGrigWindow(void)
{
   
}
WXFoxServer *WXGrigWindow::StartServer()
{
   //if client or server exist return NULL
   if(m_WXFoxClient!=NULL) return NULL;
   if(m_WXFoxServer!=NULL) return NULL;
   m_WXFoxServer = new WXFoxServer(this);
   m_WXFoxServer->m_dataLoaded = dataLoaded;
   return m_WXFoxServer;
}

WXFoxClient *WXGrigWindow::StartClientWindow()
{
   //if client or server exist return NULL
   if(m_WXFoxClient!=NULL) return NULL;
   if(m_WXFoxServer!=NULL) return NULL;
   m_WXFoxClient = new WXFoxClient(this);
   return m_WXFoxClient;
}
void WXGrigWindow::DataLoaded()
{
   dataLoaded = true;
   
   if(m_WXFoxServer!=NULL) {
      m_WXFoxServer->m_dataLoaded = true;
   }   
   
}
void WXGrigWindow::Clear()
{
   if(m_WXFoxServer!=NULL) {
      m_WXFoxServer->Clear();
      m_WXFoxServer=NULL;
   }
   
   if(m_WXFoxClient!=NULL) {
      m_WXFoxClient->Clear();
      m_WXFoxClient=NULL;
   }
}
