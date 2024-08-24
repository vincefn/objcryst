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

#include "WXFoxSlave.h"

static long CONNECT_CLIENT_BUTTON=                        WXCRYST_ID();
static long CONNECT_TIMER=                          WXCRYST_ID();
static long UPDATE_TIMER=                           WXCRYST_ID();

BEGIN_EVENT_TABLE(WXFoxSlave, wxWindow)
   EVT_BUTTON(CONNECT_CLIENT_BUTTON,                  WXFoxSlave::OnConnectClient)
   EVT_TIMER(CONNECT_TIMER,                         WXFoxSlave::OnConnectTimer)
   EVT_TIMER(UPDATE_TIMER,                         WXFoxSlave::OnUpdateProcessTimer)
END_EVENT_TABLE()

WXFoxSlave::WXFoxSlave(wxWindow* parent, wxString working_dir):
wxWindow(parent,-1)
{
   m_parent = parent;
   m_working_dir = working_dir;
   //m_FoxClient = new FoxClient(m_working_dir);
   InitClient();
   m_ConnectTimer = new wxTimer(this, CONNECT_TIMER);
   m_UpdateProcessTimer = new wxTimer(this, UPDATE_TIMER);
   m_UpdateProcessTimer->Start(3*1000, true);
   m_connecting = false;

   m_grid_slave = new FoxGridSlave(m_working_dir);
}
WXFoxSlave::~WXFoxSlave(void)
{
    Clear();
}
void WXFoxSlave::Clear()
{
   if(m_connecting) {
      m_connecting=false;
      if(m_ConnectTimer->IsRunning()) m_ConnectTimer->Stop();
   }
   /*
   if (m_FoxClient != 0) {
       delete m_FoxClient;
       m_FoxClient = 0;
   }
   */
}
void WXFoxSlave::InitClient()
{
   unsigned int xsize=600;

   wxBoxSizer *topSizer = new wxStaticBoxSizer( wxVERTICAL, this, "");
   wxBoxSizer *IPSizer = new wxBoxSizer( wxHORIZONTAL);
   wxGridSizer *connect_sizer = new wxGridSizer(2, 4, 2, 2);

   //Editbox for IP
   wxStaticText *labelIP = new wxStaticText(this, NULL, _T("Server IP:"), wxDefaultPosition, wxDefaultSize, 0 , _T("label"));
   //IPSizer->Add(labelIP, 0, wxALL|wxALIGN_LEFT|wxALIGN_TOP ,3);
   connect_sizer->Add(labelIP, 0, wxALL|wxALIGN_RIGHT|wxALIGN_TOP ,3);
   

   m_IPWindow = new wxComboBox(this, NULL, _T("localhost"), wxDefaultPosition,
                       wxDefaultSize, 0,0,
                       wxCB_DROPDOWN, wxDefaultValidator, _T("TextBox"));
   connect_sizer->Add(m_IPWindow, 0, wxALL|wxALIGN_LEFT|wxALIGN_TOP ,3);
   
   wxStaticText *labelPort = new wxStaticText(this, NULL, _T("Port:"), wxDefaultPosition, wxDefaultSize, 0 , _T("label"));
   connect_sizer->Add(labelPort, 0, wxALL|wxALIGN_RIGHT|wxALIGN_TOP ,3);

   m_portWindow = new wxTextCtrl(this, NULL, "2854", wxDefaultPosition,
                                                wxDefaultSize, 0,
                                                wxDefaultValidator, _T("TextBox"));

   connect_sizer->Add(m_portWindow, 0, wxALL|wxALIGN_LEFT|wxALIGN_TOP ,3);

   //nbCPUs window
   wxString nbCPUs;
   nbCPUs << wxThread::GetCPUCount();
   wxStaticText *label2 = new wxStaticText(this, NULL, wxString::Format(_T("CPUs (max=%d): "), wxThread::GetCPUCount()), wxDefaultPosition, wxDefaultSize, 0 , _T("label"));
   //IPSizer->Add(label2, 0, wxALL|wxALIGN_LEFT|wxALIGN_TOP ,3);
   connect_sizer->Add(label2, 0, wxALL|wxALIGN_RIGHT|wxALIGN_TOP ,3);
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
   ConnectSizer->Add(IPButtonSizer, 0, wxALL|wxALIGN_RIGHT);
   
   topSizer->Add(ConnectSizer, 0, wxALL|wxALIGN_TOP);
   topSizer->Add(new wxStaticText(this, NULL, "List of Processes", wxDefaultPosition, wxDefaultSize, 0 , _T("label")), 0, wxALL|wxALIGN_TOP, 3);
   m_process_table = new wxGrid(this, NULL, wxDefaultPosition, wxSize(xsize,200), wxWANTS_CHARS, _T("List of Processes"));
   m_process_table->CreateGrid(1,5,wxGrid::wxGridSelectRows);
   m_process_table->SetColLabelValue(0, _T("No."));
   m_process_table->SetColLabelValue(1, _T("Job ID"));
   m_process_table->SetColLabelValue(2, _T("Started"));
   m_process_table->SetColLabelValue(3, _T("Estimated progress"));
   m_process_table->SetColLabelValue(4, _T("Status"));
   m_process_table->SetColLabelSize(20);
   m_process_table->SetRowLabelSize(0);
   m_process_table->SetColSize(0, 1.0 * xsize/13.0);
   m_process_table->SetColSize(1, 3.0 * xsize/13.0);
   m_process_table->SetColSize(2, 3.0 * xsize/13.0);
   m_process_table->SetColSize(3, 3.0 * xsize/13.0);
   m_process_table->SetColSize(4, 3.0 * xsize/13.0);
   m_process_table->DeleteRows(0, 1, false);
   topSizer->Add(m_process_table, 0, wxALL);

   topSizer->Add(new wxStaticText(this, NULL, "List of processed jobs", wxDefaultPosition, wxDefaultSize, 0 , _T("label")), 0, wxALL|wxALIGN_TOP, 3);
   m_job_table = new wxGrid(this, NULL, wxDefaultPosition, wxSize(xsize,100), wxWANTS_CHARS, _T("List of Processed Jobs"));
   m_job_table->CreateGrid(1,5,wxGrid::wxGridSelectRows);
   m_job_table->SetColLabelValue(0, _T("No."));
   m_job_table->SetColLabelValue(1, _T("Job ID"));
   m_job_table->SetColLabelValue(2, _T("no. runs"));
   m_job_table->SetColLabelValue(3, _T("Times processed"));
   m_job_table->SetColLabelValue(4, _T("Average calc. time"));
   m_job_table->SetColLabelSize(20);
   m_job_table->SetRowLabelSize(0);
   m_job_table->SetColSize(0, 1.0 * xsize/10.0);
   m_job_table->SetColSize(1, 2.0 * xsize/10.0);
   m_job_table->SetColSize(2, 2.0 * xsize/10.0);
   m_job_table->SetColSize(3, 2.0 * xsize/10.0);
   m_job_table->SetColSize(4, 2.0 * xsize/10.0);
   m_job_table->DeleteRows(0, 1, false);
   topSizer->Add(m_job_table, 0, wxALL);

   SetSizer(topSizer);

   //this->LoadUsedIPs();
   wxTheApp->GetTopWindow()->Layout();
   wxTheApp->GetTopWindow()->SendSizeEvent();
}
void WXFoxSlave::OnConnectClient(wxCommandEvent& event)
{
   if(m_connecting) {
      m_connecting=false;
      if(m_ConnectTimer->IsRunning()) m_ConnectTimer->Stop();
      m_ConnectButton->SetLabel(_T("Connect"));
      m_ConnectButton->Show();
   } else {
      long nbCPUs;
      long port;
      m_nbCPUs->GetValue().ToLong(&nbCPUs);
      m_portWindow->GetValue().ToLong(&port);
      this->setNbCPU((int) nbCPUs);
      m_connecting = true;      
      ConnectClient(m_IPWindow->GetValue(), port);
   }   
}
void WXFoxSlave::setNbCPU(int nb)
{
    m_grid_slave->ResetNbCPUsAll(nb);
    wxString tmp;
    tmp << nb;
    m_nbCPUs->SetValue(tmp);
}
void WXFoxSlave::CloseClient()
{
    /*
    //kill processes before exit
    m_FoxClient->KillProcesses();

    if(m_connecting) {
        m_connecting=false;
        if(m_ConnectTimer->IsRunning()) m_ConnectTimer->Stop();
    } else if(m_FoxClient->IsClientConnected()) {
        m_FoxClient->Disconnect();
    }
    */
}
void WXFoxSlave::OnUpdateProcessTimer(wxTimerEvent& event)
{    
    if(m_grid_slave==0) {
        m_ConnectButton->Show(true);
        m_ConnectButton->SetLabel(_T("Connect"));
        Layout();
        m_UpdateProcessTimer->Start(3*1000, true);
        return;
    }

    if(!m_grid_slave->isConnectedToMaster()) {
        m_ConnectButton->Show(true);
        m_ConnectButton->SetLabel(_T("Connect"));
        Layout();
    }

    vector<FOX_PROCESS> p = m_grid_slave->getProcesses();            
    vector<MasterJob> jobs = m_grid_slave->getJobs();
    
    int nbRow = m_process_table->GetNumberRows();
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
        if(p[i].running) {
            m_process_table->SetCellValue(i,0,wxString::Format("%d",i));
            m_process_table->SetReadOnly(i,0);
            m_process_table->SetCellBackgroundColour(i, 0, wxColour(255, 200, 200));

            m_process_table->SetCellValue(i,1,wxString::Format("%d",p[i].jobID));
            m_process_table->SetReadOnly(i,1);           
            m_process_table->SetCellBackgroundColour(i, 1, wxColour(255, 200, 200));

            m_process_table->SetCellValue(i,2,wxString::Format("%s",p[i].startingtime.FormatTime()));
            m_process_table->SetReadOnly(i,2);           
            m_process_table->SetCellBackgroundColour(i, 2, wxColour(255, 200, 200));

            int cprogress = -1;            
            for(int j=0;j<jobs.size();j++) {
                if(p[i].jobID==jobs[j].ID) {
                    cprogress = p[i].getProgressInPercents(jobs[j].average_calc_time);
                }
            }
            
            if(cprogress==-1) {
                m_process_table->SetCellValue(i,3,"not available yet");
            } else {
                m_process_table->SetCellValue(i,3,wxString::Format("%d %%", cprogress));
            }
            m_process_table->SetReadOnly(i,3);
            m_process_table->SetCellBackgroundColour(i, 3, wxColour(255, 200, 200));


            m_process_table->SetCellValue(i,4,"running");
            m_process_table->SetReadOnly(i,4);
            m_process_table->SetCellBackgroundColour(i, 4, wxColour(255, 200, 200));
        } else {
            m_process_table->SetCellValue(i,0,wxString::Format("%d",i));
            m_process_table->SetReadOnly(i,0);
            m_process_table->SetCellBackgroundColour(i, 0, wxColour(200, 255, 200));

            m_process_table->SetCellValue(i,1,"");
            m_process_table->SetReadOnly(i,1);           
            m_process_table->SetCellBackgroundColour(i, 1, wxColour(200, 255, 200));

            m_process_table->SetCellValue(i,2,"");
            m_process_table->SetReadOnly(i,2);
            m_process_table->SetCellBackgroundColour(i, 2, wxColour(200, 255, 200));

            m_process_table->SetCellValue(i,3,"");
            m_process_table->SetReadOnly(i,3);
            m_process_table->SetCellBackgroundColour(i, 3, wxColour(200, 255, 200));

            m_process_table->SetCellValue(i,4,"waiting for job");
            m_process_table->SetReadOnly(i,4);
            m_process_table->SetCellBackgroundColour(i, 4, wxColour(200, 255, 200));
        }
    }

    //m_job_table
    nbRow = m_job_table->GetNumberRows();
    //delete it just in the case of different nb of processes then before
    if(nbRow != jobs.size()) {
        if(nbRow>0) {
            m_job_table->DeleteRows(0, nbRow, true);
        }
        for(int i=0;i<jobs.size();i++) {
            m_job_table->InsertRows(i, 1, false);   
        }
    }
    m_job_table->ClearGrid();
    for(int i=0;i<jobs.size();i++) {
        m_job_table->SetCellValue(i,0,wxString::Format("%d",i));
        m_job_table->SetReadOnly(i,0);

        m_job_table->SetCellValue(i,1,wxString::Format("%d",jobs[i].ID));
        m_job_table->SetReadOnly(i,1);  

         m_job_table->SetCellValue(i,2,wxString::Format("%d", jobs[i].nb_runs));
        m_job_table->SetReadOnly(i,2);

        m_job_table->SetCellValue(i,3,wxString::Format("%d", jobs[i].results.size()));
        m_job_table->SetReadOnly(i,3);           

        m_job_table->SetCellValue(i,4,jobs[i].average_calc_time.Format());
        m_job_table->SetReadOnly(i,4);
    }      

    m_UpdateProcessTimer->Start(3*1000, true);
}
void WXFoxSlave::OnConnectTimer(wxTimerEvent& event)
{
   if(!m_connecting) return;
   wxString host = m_IPWindow->GetValue();
   long port;
   m_portWindow->GetValue().ToLong(&port);
   if(!m_grid_slave->isConnectedToMaster()) {
        ConnectClient(host, port);
   }
}
void WXFoxSlave::ConnectClient(wxString IP, int port)
{   
   if(!m_grid_slave->isConnectedToMaster()) {       
       m_ConnectButton->SetLabel(_T("Cancel Connecting"));   
       if(m_grid_slave->ConnectToMaster(1, IP, port)) {           
          m_ConnectButton->Hide();
          Layout();
          wxString CPUs;
          m_connecting=false;
       }else{
          m_ConnectTimer->Start(10000, true);       
       }
   }   
}
