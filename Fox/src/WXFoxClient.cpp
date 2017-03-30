#ifdef __WX__CRYST__
   #include "ObjCryst/wxCryst/wxCrystal.h"
#endif

// Fox client/server grid
/*
#include "wx/socket.h"
#include <wx/file.h>
#include <wx/filefn.h>
*/

using namespace ObjCryst;
using namespace std;

#include "WXFoxClient.h"

static long CONNECT_CLIENT_BUTTON=                        WXCRYST_ID();
static long CONNECT_TIMER=                          WXCRYST_ID();

BEGIN_EVENT_TABLE(WXFoxClient, wxWindow)
   EVT_BUTTON(CONNECT_CLIENT_BUTTON,                  WXFoxClient::OnConnectClient)
   EVT_TIMER(CONNECT_TIMER,                         WXFoxClient::OnConnectTimer)
END_EVENT_TABLE()

WXFoxClient::WXFoxClient(wxWindow* parent, wxString working_dir):
wxWindow(parent,-1)
{
   m_parent = parent;
   m_working_dir = working_dir;
   m_FoxClient = new FoxClient(m_working_dir);
   InitClient();
   m_ConnectTimer = new wxTimer(this, CONNECT_TIMER);
   m_connecting = false;
}
WXFoxClient::~WXFoxClient(void)
{
    Clear();
}
void WXFoxClient::Clear()
{
   if(m_connecting) {
      m_connecting=false;
      if(m_ConnectTimer->IsRunning()) m_ConnectTimer->Stop();
   }
   if (m_FoxClient != 0) {
       delete m_FoxClient;
       m_FoxClient = 0;
   }
}
void WXFoxClient::InitClient()
{
   unsigned int xsize=400;

   wxBoxSizer *topSizer = new wxBoxSizer( wxVERTICAL);
   wxBoxSizer *IPSizer = new wxBoxSizer( wxVERTICAL);

   //Editbox for IP
   wxStaticText *labelIP = new wxStaticText(this, NULL, _T("Server IP:"), wxDefaultPosition, wxDefaultSize, 0 , _T("label"));
   IPSizer->Add(labelIP, 0, wxALL|wxALIGN_LEFT|wxALIGN_TOP ,3);

   m_IPWindow = new wxComboBox(this, NULL, _T("localhost"), wxDefaultPosition,
                       wxDefaultSize, 0,0,
                       wxCB_DROPDOWN, wxDefaultValidator, _T("TextBox"));

   IPSizer->Add(m_IPWindow, 0, wxALL|wxALIGN_LEFT|wxALIGN_TOP ,3);

   //nbCPUs window
   wxString nbCPUs;
   nbCPUs << wxThread::GetCPUCount();
   wxStaticText *label2 = new wxStaticText(this, NULL, _T("Set number of available CPUs or cores: "), wxDefaultPosition, wxDefaultSize, 0 , _T("label"));
   IPSizer->Add(label2, 0, wxALL|wxALIGN_LEFT|wxALIGN_TOP ,3);
   m_nbCPUs = new wxTextCtrl(this, NULL, nbCPUs, wxDefaultPosition,
                                                wxDefaultSize, 0,
                                                wxDefaultValidator, _T("TextBox"));

   IPSizer->Add(m_nbCPUs, 0, wxALL|wxALIGN_LEFT|wxALIGN_TOP ,3);

   //Connect Button
   wxBoxSizer *IPButtonSizer = new wxBoxSizer( wxVERTICAL);
   m_ConnectButton = new wxButton(this, CONNECT_CLIENT_BUTTON, _T("Connect"), wxDefaultPosition, wxDefaultSize, 0, wxDefaultValidator, _T("Button1"));
   m_ConnectButton->Show();
   IPButtonSizer->Add(m_ConnectButton,0, wxALL|wxALIGN_LEFT|wxALIGN_TOP ,3);


   wxBoxSizer *ConnectSizer = new wxBoxSizer( wxVERTICAL);
   ConnectSizer->Add(IPSizer, 0, wxALL|wxALIGN_TOP);
   ConnectSizer->Add(IPButtonSizer, 0, wxALL|wxALIGN_BOTTOM);

   wxBoxSizer *eventSizer = new wxBoxSizer( wxVERTICAL );

   //events window
   wxStaticText *label1 = new wxStaticText(this, NULL, _T("Client events: "), wxDefaultPosition, wxDefaultSize, 0 , _T("label"));
   eventSizer->Add(label1,0, wxALL|wxALIGN_LEFT|wxALIGN_TOP ,3);

   m_EventsWindow = new wxTextCtrl(this, NULL, wxEmptyString, wxDefaultPosition,
                                                wxSize(xsize,50), wxTE_MULTILINE|wxTE_READONLY,
                                                wxDefaultValidator, _T("TextBox"));

   eventSizer->Add(m_EventsWindow,0, wxALL|wxALIGN_LEFT|wxALIGN_TOP ,3);



   topSizer->Add(eventSizer, 0, wxALL|wxALIGN_TOP);
   topSizer->Add(ConnectSizer, 0, wxALL|wxALIGN_TOP);

   SetSizer(topSizer);

   //this->LoadUsedIPs();
   wxTheApp->GetTopWindow()->Layout();
   wxTheApp->GetTopWindow()->SendSizeEvent();
}
void WXFoxClient::OnConnectClient(wxCommandEvent& event)
{
   if(m_connecting) {
      m_connecting=false;
      if(m_ConnectTimer->IsRunning()) m_ConnectTimer->Stop();
      m_ConnectButton->SetLabel(_T("Connect"));
      m_EventsWindow->AppendText(_T("connecting canceled."));
   }else if(m_FoxClient->IsClientConnected()) {
      wxMessageDialog d(this,_T("Are you sure you want to disconnect?\n It aborts conmuting of this client."), _T(""), wxYES | wxNO);
      if(wxID_YES!=d.ShowModal()) return;

      m_FoxClient->Disconnect();
      m_FoxClient->KillProcesses();
      m_ConnectButton->SetLabel(_T("Connect"));
      m_EventsWindow->SetInsertionPointEnd();
      m_EventsWindow->AppendText(_T("Client was disconnected.\n"));
   } else {
      long nbCPUs;
      m_nbCPUs->GetValue().ToLong((long *) &nbCPUs);
      this->setNbCPU((int) nbCPUs);
      m_connecting = true;
      m_ConnectTimer->Start(2000, false);
   }
}
void WXFoxClient::setNbCPU(int nb)
{
    m_FoxClient->setNbOfAvailCPUs(nb);
    wxString tmp;
    tmp << nb;
    m_nbCPUs->SetValue(tmp);
}
void WXFoxClient::CloseClient()
{
    //kill processes before exit
    m_FoxClient->KillProcesses();

    if(m_connecting) {
        m_connecting=false;
        if(m_ConnectTimer->IsRunning()) m_ConnectTimer->Stop();
    } else if(m_FoxClient->IsClientConnected()) {
        m_FoxClient->Disconnect();
    }
}
void WXFoxClient::OnConnectTimer(wxTimerEvent& event)
{
   if(!m_connecting) return;
   wxString host = m_IPWindow->GetValue();
   ConnectClient(host);
}
void WXFoxClient::ConnectClient(wxString IP)
{
   m_EventsWindow->Clear();
   m_EventsWindow->SetValue(_T("Connecting to: \"") + IP + _T("\"... "));
   m_ConnectButton->SetLabel(_T("Stop"));
   m_FoxClient->ConnectClient(1, IP);

   if(m_FoxClient->IsClientConnected()){
      m_ConnectButton->SetLabel(_T("Disconnect"));
      m_EventsWindow->Clear();
      m_EventsWindow->SetValue(_T(".:Client is connected to the server: \"") + IP + _T("\":.\n"));
      wxString CPUs;
      CPUs.Printf(_T("Number of available CPUs is: %d\n"), m_FoxClient->getNbOfAvailCPUs());
      m_EventsWindow->AppendText(CPUs);
      m_connecting=false;
   }else{
      m_ConnectTimer->Start(10000, true);
   }
}
