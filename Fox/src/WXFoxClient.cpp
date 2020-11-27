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
static long UPDATE_TIMER=                           WXCRYST_ID();

BEGIN_EVENT_TABLE(WXFoxClient, wxWindow)
   EVT_BUTTON(CONNECT_CLIENT_BUTTON,                  WXFoxClient::OnConnectClient)
   EVT_TIMER(CONNECT_TIMER,                         WXFoxClient::OnConnectTimer)
   EVT_TIMER(UPDATE_TIMER,                         WXFoxClient::OnUpdateProcessTimer)
END_EVENT_TABLE()

WXFoxClient::WXFoxClient(wxWindow* parent, wxString working_dir):
wxWindow(parent,-1)
{
   m_parent = parent;
   m_working_dir = working_dir;
   m_FoxClient = new FoxClient(m_working_dir);
   InitClient();
   m_ConnectTimer = new wxTimer(this, CONNECT_TIMER);
   m_UpdateProcessTimer = new wxTimer(this, UPDATE_TIMER);
   m_UpdateProcessTimer->Start(10000, false);
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
   unsigned int xsize=600;

   wxBoxSizer *topSizer = new wxStaticBoxSizer( wxVERTICAL, this, "");
   wxBoxSizer *IPSizer = new wxBoxSizer( wxHORIZONTAL);
   wxGridSizer *connect_sizer = new wxGridSizer(2, 2, 10, 10);

   //Editbox for IP
   wxStaticText *labelIP = new wxStaticText(this, NULL, _T("Server IP:"), wxDefaultPosition, wxDefaultSize, 0 , _T("label"));
   //IPSizer->Add(labelIP, 0, wxALL|wxALIGN_LEFT|wxALIGN_TOP ,3);
   connect_sizer->Add(labelIP, 0, wxALL|wxALIGN_LEFT|wxALIGN_TOP ,3);

   m_IPWindow = new wxComboBox(this, NULL, _T("localhost"), wxDefaultPosition,
                       wxDefaultSize, 0,0,
                       wxCB_DROPDOWN, wxDefaultValidator, _T("TextBox"));

   //IPSizer->Add(m_IPWindow, 0, wxALL|wxALIGN_LEFT|wxALIGN_TOP ,3);
   connect_sizer->Add(m_IPWindow, 0, wxALL|wxALIGN_LEFT|wxALIGN_TOP ,3);

   //nbCPUs window
   wxString nbCPUs;
   nbCPUs << wxThread::GetCPUCount();
   wxStaticText *label2 = new wxStaticText(this, NULL, wxString::Format(_T("CPU cores allowed for computing (max=%d): "), wxThread::GetCPUCount()), wxDefaultPosition, wxDefaultSize, 0 , _T("label"));
   //IPSizer->Add(label2, 0, wxALL|wxALIGN_LEFT|wxALIGN_TOP ,3);
   connect_sizer->Add(label2, 0, wxALL|wxALIGN_LEFT|wxALIGN_TOP ,3);
   m_nbCPUs = new wxTextCtrl(this, NULL, nbCPUs, wxDefaultPosition,
                                                wxDefaultSize, 0,
                                                wxDefaultValidator, _T("TextBox"));

   //IPSizer->Add(m_nbCPUs, 0, wxALL|wxALIGN_LEFT|wxALIGN_TOP ,3);
   connect_sizer->Add(m_nbCPUs, 0, wxALL|wxALIGN_LEFT|wxALIGN_TOP ,3);

   //Connect Button
   wxBoxSizer *IPButtonSizer = new wxBoxSizer( wxVERTICAL);
   m_ConnectButton = new wxButton(this, CONNECT_CLIENT_BUTTON, _T("Connect"), wxDefaultPosition, wxDefaultSize, 0, wxDefaultValidator, _T("Button1"));
   m_ConnectButton->Show();
   IPButtonSizer->Add(m_ConnectButton,0, wxALL|wxALIGN_RIGHT|wxALIGN_TOP ,3);


   wxStaticBoxSizer *ConnectSizer = new wxStaticBoxSizer( wxVERTICAL, this, "Connection");
   //ConnectSizer->Add(IPSizer, 0, wxALL|wxALIGN_TOP);
   ConnectSizer->Add(connect_sizer, 0, wxALL|wxALIGN_TOP);
   ConnectSizer->Add(IPButtonSizer, 0, wxALL|wxALIGN_BOTTOM|wxALIGN_RIGHT);
   /*
   wxBoxSizer *eventSizer = new wxBoxSizer( wxVERTICAL );

   //events window

   wxStaticText *label1 = new wxStaticText(this, NULL, _T("Client events: "), wxDefaultPosition, wxDefaultSize, 0 , _T("label"));
   eventSizer->Add(label1,0, wxALL|wxALIGN_LEFT|wxALIGN_TOP ,3);

   m_EventsWindow = new wxTextCtrl(this, NULL, wxEmptyString, wxDefaultPosition,
                                                wxSize(xsize,50), wxTE_MULTILINE|wxTE_READONLY,
                                                wxDefaultValidator, _T("TextBox"));

   eventSizer->Add(m_EventsWindow,0, wxALL|wxALIGN_LEFT|wxALIGN_TOP ,3);


   topSizer->Add(eventSizer, 0, wxALL|wxALIGN_TOP);
   */
   topSizer->Add(ConnectSizer, 0, wxALL|wxALIGN_TOP);
   topSizer->Add(new wxStaticText(this, NULL, "List of Processes", wxDefaultPosition, wxDefaultSize, 0 , _T("label")), 0, wxALL|wxALIGN_TOP, 3);
   m_process_table = new wxGrid(this, NULL, wxDefaultPosition, wxSize(xsize,200), wxWANTS_CHARS, _T("List of Processes"));
   m_process_table->CreateGrid(1,4,wxGrid::wxGridSelectRows);
   m_process_table->SetColLabelValue(0, _T("No."));
   m_process_table->SetColLabelValue(1, _T("Job ID"));
   m_process_table->SetColLabelValue(2, _T("Started"));
   m_process_table->SetColLabelValue(3, _T("Status"));
   m_process_table->SetColLabelSize(20);
   m_process_table->SetRowLabelSize(0);
   m_process_table->SetColumnWidth(0, 1.0 * xsize/10.0);
   m_process_table->SetColumnWidth(1, 3.0 * xsize/10.0);
   m_process_table->SetColumnWidth(2, 3.0 * xsize/10.0);
   m_process_table->SetColumnWidth(3, 3.0 * xsize/10.0);
   m_process_table->DeleteRows(0, 1, false);
   topSizer->Add(m_process_table, 0, wxALL|wxALIGN_BOTTOM);

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
      //m_EventsWindow->AppendText(_T("connecting canceled."));
   }else if(m_FoxClient->IsClientConnected()) {
      wxMessageDialog d(this,_T("Are you sure you want to disconnect?\n It aborts conmuting of this client."), _T(""), wxYES | wxNO);
      if(wxID_YES!=d.ShowModal()) return;

      m_FoxClient->Disconnect();
      m_FoxClient->KillProcesses();
      m_ConnectButton->SetLabel(_T("Connect"));
      //m_EventsWindow->SetInsertionPointEnd();
      //m_EventsWindow->AppendText(_T("Client was disconnected.\n"));
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
void WXFoxClient::OnUpdateProcessTimer(wxTimerEvent& event)
{
    if(m_FoxClient==0) return;

    vector<FoxProcess> p = m_FoxClient->get_copy_of_processes();
    
    int nbRow = m_process_table->GetRows();
    
    //delete it just in the case of different nb of processes then before
    if(nbRow != p.size()) {
        if(nbRow>0) {
            m_process_table->DeleteRows(0, nbRow, true);
        }
        for(int i=0;i<p.size();i++) {
            m_process_table->InsertRows(i, 1, false);   
        }
    }

    m_process_table->ClearGrid();
    for(int i=0;i<p.size();i++) {
        if(p[i].isRunning()) {
            m_process_table->SetCellValue(i,0,wxString::Format("%d",i));
            m_process_table->SetReadOnly(i,0);
            m_process_table->SetCellBackgroundColour(wxColor(255, 200, 200), i, 0);

            m_process_table->SetCellValue(i,1,wxString::Format("%d",p[i].getJobID()));
            m_process_table->SetReadOnly(i,1);           
            m_process_table->SetCellBackgroundColour(wxColor(255, 200, 200), i, 1);

            m_process_table->SetCellValue(i,2,wxString::Format("%s",p[i].getStartingTime().FormatTime()));
            m_process_table->SetReadOnly(i,2);           
            m_process_table->SetCellBackgroundColour(wxColor(255, 200, 200), i, 2);

            m_process_table->SetCellValue(i,3,"running");
            m_process_table->SetReadOnly(i,3);
            m_process_table->SetCellBackgroundColour(wxColor(255, 200, 200), i, 3);
        } else {
            m_process_table->SetCellValue(i,0,wxString::Format("%d",i));
            m_process_table->SetReadOnly(i,0);
            m_process_table->SetCellBackgroundColour(wxColor(200, 255, 200), i, 0);

            m_process_table->SetCellValue(i,1,"");
            m_process_table->SetReadOnly(i,1);           
            m_process_table->SetCellBackgroundColour(wxColor(200, 255, 200), i, 1);

            m_process_table->SetCellValue(i,2,"");
            m_process_table->SetReadOnly(i,2);
            m_process_table->SetCellBackgroundColour(wxColor(200, 255, 200), i, 2);

            m_process_table->SetCellValue(i,3,"waiting for job");
            m_process_table->SetReadOnly(i,3);
            m_process_table->SetCellBackgroundColour(wxColor(200, 255, 200), i, 3);
        }
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
   //m_EventsWindow->Clear();
   //m_EventsWindow->SetValue(_T("Connecting to: \"") + IP + _T("\"... "));
   m_ConnectButton->SetLabel(_T("Stop"));
   m_FoxClient->ConnectClient(1, IP);

   if(m_FoxClient->IsClientConnected()){
      m_ConnectButton->SetLabel(_T("Disconnect"));
      //m_EventsWindow->Clear();
      //m_EventsWindow->SetValue(_T(".:Client is connected to the server: \"") + IP + _T("\":.\n"));
      wxString CPUs;
      CPUs.Printf(_T("Number of available CPUs is: %d\n"), m_FoxClient->getNbOfAvailCPUs());
      //m_EventsWindow->AppendText(CPUs);
      m_connecting=false;
   }else{
      m_ConnectTimer->Start(10000, true);
   }
}
