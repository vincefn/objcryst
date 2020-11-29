
#include <wx/textdlg.h>
#include "wx/filename.h"
#include "wx/sstream.h"

#ifdef __WX__CRYST__
   #include "ObjCryst/wxCryst/wxCrystal.h"
#endif

// Fox client/server grid


using namespace ObjCryst;
using namespace std;

#include "WXFoxServer.h"
#include <iostream>


static long STOP_TO_CLIENT=                        WXCRYST_ID();
static long RUN_LOCAL_CLIENT=                     WXCRYST_ID();
static long RUN_TO_CLIENT=                        WXCRYST_ID();
static long NEW_JOB=                            WXCRYST_ID();
static long LOAD_JOB=                           WXCRYST_ID();
static long EDIT_JOB=                           WXCRYST_ID();
static long DELETE_JOB=                           WXCRYST_ID();
static long ID_UPDATE_TIMER=                     WXCRYST_ID();
static long GRID_RESULT_LIST=                     WXCRYST_ID();
static long GRID_JOB_LIST=                        WXCRYST_ID();
static long SHOW_RESULTS=                        WXCRYST_ID();
static long SHOW_RESULTS_SERVER=                        WXCRYST_ID();


BEGIN_EVENT_TABLE(WXFoxServer, wxWindow)
  //EVT_BUTTON(STOP_TO_CLIENT,                  WXFoxServer::StopAllClients)
  EVT_BUTTON(RUN_LOCAL_CLIENT,                 WXFoxServer::RunLocalClient)
  //EVT_BUTTON(RUN_TO_CLIENT,                 WXFoxServer::RunALLClient)
  EVT_BUTTON(NEW_JOB,                        WXFoxServer::OnNewJob)
  EVT_BUTTON(EDIT_JOB,                        WXFoxServer::OnEditJob)
  EVT_BUTTON(DELETE_JOB,                     WXFoxServer::OnDeleteJob)
  EVT_BUTTON(LOAD_JOB,                        WXFoxServer::OnLoadJob)
  EVT_TIMER(ID_UPDATE_TIMER,                  WXFoxServer::UpdateLists)
  EVT_GRID_CMD_CELL_LEFT_CLICK(GRID_RESULT_LIST,   WXFoxServer::OnGridResultClick)
  EVT_GRID_CMD_CELL_LEFT_CLICK(GRID_JOB_LIST,      WXFoxServer::OnGridJobClick)
  EVT_BUTTON(SHOW_RESULTS,                     WXFoxServer::OnShowResults)
  EVT_BUTTON(SHOW_RESULTS_SERVER,                     WXFoxServer::OnShowResultsServer)

END_EVENT_TABLE()

WXFoxServer::WXFoxServer(wxWindow* parent, wxString workingDir):
wxWindow(parent,-1),m_parent(parent)
{
   m_dataLoaded = false;
   //m_WXFoxServerDataMutex = new wxMutex();
   m_working_dir = workingDir;
   InitServer();
}
WXFoxServer::~WXFoxServer(void)
{

}
void WXFoxServer::Clear()
{
   m_UpdateTimer->Stop();
   delete m_UpdateTimer;
   //todo: clear m_jobs, m_results
   //delete m_WXFoxServerDataMutex;
   delete m_FoxServer;
}
void WXFoxServer::InitServer()
{
   //starting server
   m_FoxServer = new FoxServer();
   m_FoxServer->SetWorkingDir(m_working_dir);
   m_FoxServer->StartGridServer();

   //Start update timer
   m_UpdateTimer = new wxTimer(this, ID_UPDATE_TIMER);
   m_UpdateTimer->Start(10*1000, false);

   //display all necessary controls
   wxBoxSizer *topSizer = new wxBoxSizer( wxVERTICAL );

   unsigned int xsize=650;
   //Job List table
   wxBoxSizer *JobSizer = new wxBoxSizer( wxVERTICAL );
   wxStaticText *JobLabel = new wxStaticText(m_parent, NULL, _T("Job List: "), wxDefaultPosition, wxDefaultSize, 0 , _T("label"));
   JobSizer->Add(JobLabel,0, wxALL|wxALIGN_LEFT|wxALIGN_TOP ,3);


   m_JobListTable = new wxGrid(this,GRID_JOB_LIST, wxDefaultPosition, wxSize(xsize,150), wxWANTS_CHARS, _T("m_JobListTable"));
   m_JobListTable->CreateGrid(1,5,wxGrid::wxGridSelectRows);
   m_JobListTable->SetColLabelValue(0, _T("Job ID"));
   m_JobListTable->SetColLabelValue(1, _T("Job name"));
   m_JobListTable->SetColLabelValue(2, _T("Nb trials"));
   m_JobListTable->SetColLabelValue(3, _T("Randomize ?"));
   m_JobListTable->SetColLabelValue(4, _T("runs/done/active"));

   m_JobListTable->SetColLabelSize(20);
   m_JobListTable->SetRowLabelSize(0);
   m_JobListTable->SetColumnWidth(0, (xsize-10)/5);
   m_JobListTable->SetColumnWidth(1, (xsize-10)/5+20);
   m_JobListTable->SetColumnWidth(2, (xsize-10)/5-20);
   m_JobListTable->SetColumnWidth(3, (xsize-10)/5-10);
   m_JobListTable->SetColumnWidth(4, (xsize-10)/5+10);
   m_JobListTable->DeleteRows(0, 1, false);
   JobSizer->Add(m_JobListTable,0, wxALL|wxALIGN_LEFT|wxALIGN_TOP ,3);

   wxBoxSizer *JobButtonSizer = new wxBoxSizer( wxHORIZONTAL );
   wxButton *LoadJobButton = new wxButton(this, LOAD_JOB, _T("Create Job (from xml)"), wxDefaultPosition, wxDefaultSize, 0, wxDefaultValidator, _T(""));
   JobButtonSizer->Add(LoadJobButton,0, wxALL|wxALIGN_LEFT|wxALIGN_TOP ,3);
   LoadJobButton->SetToolTip(_T("Create a job from an existing xml file"));

   wxButton *AddJobButton = new wxButton(this, NEW_JOB, _T("Create Job (from mem)"), wxDefaultPosition, wxDefaultSize, 0, wxDefaultValidator, _T(""));
   JobButtonSizer->Add(AddJobButton,0, wxALL|wxALIGN_LEFT|wxALIGN_TOP ,3);
   AddJobButton->SetToolTip(_T("Create job from objects currently in memory\n\n")
                            _T("This requires that you have one global \n")
                            _T("optimization obect in memory."));

   wxButton *EditJobButton = new wxButton(this, EDIT_JOB, _T("Edit Job"), wxDefaultPosition, wxDefaultSize, 0, wxDefaultValidator, _T(""));
   JobButtonSizer->Add(EditJobButton,0, wxALL|wxALIGN_LEFT|wxALIGN_TOP ,3);

   wxButton *DeleteJobButton = new wxButton(this, DELETE_JOB, _T("Delete Job"), wxDefaultPosition, wxDefaultSize, 0, wxDefaultValidator, _T(""));
   JobButtonSizer->Add(DeleteJobButton,0, wxALL|wxALIGN_LEFT|wxALIGN_TOP ,3);
   JobSizer->Add(JobButtonSizer,0, wxALL|wxALIGN_LEFT|wxALIGN_TOP ,3);

   //run button
   //wxButton *RunButton = new wxButton(this, RUN_TO_CLIENT, _T("Ping clients"), wxDefaultPosition, wxDefaultSize, 0, wxDefaultValidator, _T("Button2"));
   //JobButtonSizer->Add(RunButton,0, wxALL|wxALIGN_LEFT|wxALIGN_TOP ,3);
   //RunButton->SetToolTip(_T("It asks connected clients whether they compute.\n If not, it sends them job\n\n")
   //                      _T("This is useful when some client does not compute.\n Use it carefully, it increases traffic on the server."));

   //stop button
   /*
   wxButton *StopButton = new wxButton(this, STOP_TO_CLIENT, _T("Stop"), wxDefaultPosition, wxDefaultSize, 0, wxDefaultValidator, _T("Button1"));
   JobButtonSizer->Add(StopButton,0, wxALL|wxALIGN_LEFT|wxALIGN_TOP ,3);
   */

   //client list label
   wxBoxSizer *listClientSizer = new wxBoxSizer( wxVERTICAL );
   wxStaticText *label2 = new wxStaticText(this, NULL, _T("Client list: "), wxDefaultPosition, wxDefaultSize, 0 , _T("label"));
   listClientSizer->Add(label2,0, wxALL|wxALIGN_LEFT|wxALIGN_TOP ,3);

   xsize = xsize/2-5;
   //client list table
   m_ClientTable = new wxGrid(this, NULL, wxDefaultPosition, wxSize(xsize,150), wxWANTS_CHARS, _T("m_ClientTable"));
   m_ClientTable->CreateGrid(1,3,wxGrid::wxGridSelectRows);
   m_ClientTable->SetColLabelValue(0, _T("Name"));
   m_ClientTable->SetColLabelValue(1, _T("ID"));
   m_ClientTable->SetColLabelValue(2, _T("CPUs/using"));
   //m_ClientTable->SetColLabelValue(3, _T("Status"));
   m_ClientTable->SetColLabelSize(20);
   m_ClientTable->SetRowLabelSize(0);
   m_ClientTable->SetColumnWidth(0, xsize/3);
   m_ClientTable->SetColumnWidth(1, xsize/3);
   m_ClientTable->SetColumnWidth(2, xsize/3);
   //m_ClientTable->SetColumnWidth(3, xsize/8);
   m_ClientTable->DeleteRows(0, 1, false);
   listClientSizer->Add(m_ClientTable,0, wxALL|wxALIGN_LEFT|wxALIGN_TOP ,3);

   //run Run local client button
   wxButton *RunClientButton = new wxButton(this, RUN_LOCAL_CLIENT, _T("Run local client"), wxDefaultPosition, wxDefaultSize, 0, wxDefaultValidator, _T("Button2"));
   listClientSizer->Add(RunClientButton,0, wxALL|wxALIGN_LEFT|wxALIGN_TOP ,3);
   RunClientButton->SetToolTip(_T("Gives you the option to run the job \n")
                               _T("locally, using multiple processors or cores"));



   //result list label
   wxBoxSizer *listResultSizer = new wxBoxSizer( wxVERTICAL );
   wxStaticText *label3 = new wxStaticText(this, NULL, _T("Result list: "), wxDefaultPosition, wxDefaultSize, 0 , _T("label"));
   listResultSizer->Add(label3,0, wxALL|wxALIGN_LEFT|wxALIGN_TOP ,3);

   //result list table
   m_ResultTable = new wxGrid(this, GRID_RESULT_LIST, wxDefaultPosition, wxSize(xsize,150), wxWANTS_CHARS, _T("m_ResultTable"));
   m_ResultTable->CreateGrid(1,3,wxGrid::wxGridSelectCells);
   m_ResultTable->SetColLabelValue(0, _T("Nb"));
   m_ResultTable->SetColLabelValue(1, _T("Job ID"));
   m_ResultTable->SetColLabelValue(2, _T("Cost"));
   //m_ResultTable->SetColLabelValue(3, _T("Show"));
   m_ResultTable->SetColLabelSize(20);
   m_ResultTable->SetRowLabelSize(0);
   m_ResultTable->SetColumnWidth(0, 50);
   m_ResultTable->SetColumnWidth(1, (xsize-50)/2);
   m_ResultTable->SetColumnWidth(2, (xsize-50)/2);
   //m_ResultTable->SetColumnWidth(3, 50);
   m_ResultTable->DeleteRows(0, 1, false);
   listResultSizer->Add(m_ResultTable,0, wxALL|wxALIGN_LEFT|wxALIGN_TOP ,3);

   wxBoxSizer *ResultListButtonsSizer = new wxBoxSizer( wxHORIZONTAL );
   wxButton *ShowButtom = new wxButton(this, SHOW_RESULTS, _T("Open by new FOX"), wxDefaultPosition, wxDefaultSize, 0, wxDefaultValidator, _T("Button3"));
   ResultListButtonsSizer->Add(ShowButtom,0, wxALL|wxALIGN_LEFT|wxALIGN_TOP ,3);

   wxButton *ShowButton = new wxButton(this, SHOW_RESULTS_SERVER, _T("Open by server"), wxDefaultPosition, wxDefaultSize, 0, wxDefaultValidator, _T("Button3"));
   ResultListButtonsSizer->Add(ShowButton,0, wxALL|wxALIGN_LEFT|wxALIGN_TOP ,3);

   listResultSizer->Add(ResultListButtonsSizer, 0, wxALL|wxALIGN_LEFT|wxALIGN_TOP ,3);

   wxBoxSizer *listsSizer = new wxBoxSizer( wxHORIZONTAL );

   listsSizer->Add(listClientSizer, 0, wxALL|wxALIGN_TOP);
   listsSizer->Add(listResultSizer, 0, wxALL|wxALIGN_TOP);

   //topSizer->Add(eventSizer, 0, wxALL|wxALIGN_TOP);
   topSizer->Add(JobSizer, 0, wxALL|wxALIGN_TOP);
   topSizer->Add(listsSizer, 0, wxALL|wxALIGN_TOP);
   //topSizer->Add(ButtonSizer, 0, wxALL|wxALIGN_TOP);

   this->SetSizer( topSizer );
   wxTheApp->GetTopWindow()->Layout();
   wxTheApp->GetTopWindow()->SendSizeEvent();
}
int WXFoxServer::GenerateJobID()
{
   wxString filename;
   int ID;
   bool again;
   do{
      again=false;
      long long date = (long long) time(0);
      while(date>1000000000){date-=1000000000;}//int is only to 2,147,483,647...
      ID = (int) date;
      for(int i=0;i<m_jobs.size();i++)
      {
         if(ID==m_jobs[i].getM_ID()) again=true;
      }
   }while(again);

   return ID;
}
void WXFoxServer::saveJobHeader(wxString filename, int ID, wxString name, long nbOfTrial, int nbRun, bool rand) {
    wxString data;
    wxString file;
    if(this->LoadFile(filename, file)) {
        int r = (int) rand;
        data.Printf(_T("<FoxJob Name=\"%s\" ID=\"%d\" nbOfTrial=\"%ld\" nbRun=\"%d\" rand=\"%d\"></FoxJob>\n"), name.c_str(), ID, nbOfTrial, nbRun, r);
        data +=file;
        SaveDataAsFile(data, filename);
    }
}
void WXFoxServer::ChangeJobHeader(wxString filename, int ID, wxString name, long nbOfTrial, int nbRun, bool rand) {
    wxString file;
    wxString data;
    if(LoadFile(filename, file)) {
        int r = (int) rand;
        int pos = file.First(_T("</FoxJob>"));
        if(pos!=-1) {
            pos +=9;
            file = file.Mid(pos);
        }
        data.Printf(_T("<FoxJob Name=\"%s\" ID=\"%d\" nbOfTrial=\"%d\" nbRun=\"%d\" rand=\"%ld\"></FoxJob>\n"), name.c_str(), ID, nbOfTrial, nbRun, r);
        data +=file;
        SaveDataAsFile(data, filename);
    }
}
bool WXFoxServer::LoadFile(wxString filename, wxString &in)
{
   wxFile infile(filename, wxFile::read);
   in.Clear();
   if(!infile.IsOpened()) return false;
   long len = infile.Length();
   char *buffer;
   buffer = (char*) calloc(len+1, sizeof(char));
   infile.Read(buffer, len);
   infile.Close();
   in = wxString::FromAscii(buffer);
   free(buffer);
   return true;
}
void WXFoxServer::SaveDataAsFile(wxString out, wxString filename)
{
   wxFile outFile(filename, wxFile::write);
   if(outFile.IsOpened())
   {
      outFile.Write(out);
      outFile.Close();
   }
}
void WXFoxServer::OnNewJob(wxCommandEvent& event)
{
   //if(!m_dataLoaded)
   if(gOptimizationObjRegistry.GetNb()==0)
   {
      wxMessageBox(_T("You need at least one optimization object !"), _T("No data to optimize"), wxOK | wxICON_INFORMATION, this);
      return;
   }

   WXCrystValidateAllUserInput();

   int newID = GenerateJobID();
   wxString Name = _T("jobname");
   long nbOfTrial = 1000000;
   long nbRun = 10;
   bool randomize = true;

   do{//show setup window
      if(!ShowSetJobWindow(newID, Name, nbOfTrial, nbRun, randomize)) return;
      cout<<nbRun<<","<<nbOfTrial<<","<<Name<<endl;
   }while((nbRun<=0)||(nbOfTrial<=0)||(Name==_T("")));

   /*
   wxString filename;
   filename.Printf(_T("JOB_%d.xml"), newID);
   #ifdef WIN32
       filename = m_working_dir + _T("\\") + filename;
   #else
       filename = m_working_dir + _T("/") + filename;
   #endif
   */
   //for unicode only
   std::string tmp("");
   std::stringstream str;
   #ifdef WIN32
   str<<m_working_dir<<"\\"<<"JOB_"<<newID<<".xml";
   #else
   str<<m_working_dir<<"/"<<"JOB_"<<newID<<".xml";
   #endif
   tmp = str.str();
   //unicode
   XMLCrystFileSaveGlobal(tmp);
   saveJobHeader(wxString::FromAscii(tmp.c_str()), newID, Name, nbOfTrial, nbRun, randomize);
   AddServerJob(wxString::FromAscii(tmp.c_str()), Name, newID, nbOfTrial, nbRun, randomize);
}
bool WXFoxServer::isFileFoxJob(wxString path, wxString &name, int &id, long &nbOfTrial, long &nbRun, bool &rand) {
    wxString ID, Tr, Run, Rand;

    ifstream is(path.ToAscii());
    if(!is){
        return false;
    }
    XMLCrystTag tag;
    do {
      if(true==is.eof()) return false;
      is>>tag;
    } while("FoxJob"!=tag.GetName());
    for(int i=tag.GetNbAttribute()-1;i>=0;i--) {
        if(tag.GetAttributeName(i)=="Name"){
            name = wxString::FromAscii(tag.GetAttributeValue(i).c_str());
        }
        if(tag.GetAttributeName(i)=="ID"){
           ID = wxString::FromAscii(tag.GetAttributeValue(i).c_str());
           ID.ToLong((long*) &id);
        }
        if(tag.GetAttributeName(i)=="nbOfTrial"){
           Tr = wxString::FromAscii(tag.GetAttributeValue(i).c_str());
           Tr.ToLong((long*) &nbOfTrial);
        }
        if(tag.GetAttributeName(i)=="nbRun"){
           Run = wxString::FromAscii(tag.GetAttributeValue(i).c_str());
           Run.ToLong((long*) &nbRun);
        }
        if(tag.GetAttributeName(i)=="rand"){
           long r;
           Rand = wxString::FromAscii(tag.GetAttributeValue(i).c_str());
           Rand.ToLong((long*) &r);
           rand = (bool) r;
        }
    }
    is.close();
    return true;
}
bool WXFoxServer::isJobLoaded(long ID) {
    //if(m_WXFoxServerDataMutex->Lock()!=wxMUTEX_NO_ERROR) return true;
    for(int i=0;i<m_jobs.size();i++) {
        if(m_jobs[i].getM_ID()==ID) return true;
    }
    //m_WXFoxServerDataMutex->Unlock();
    return false;
}
void WXFoxServer::AddServerJob(wxString filename, wxString name, int id, long nbOfTrial, long nbRun, bool rand) {
    //save it in joblist
    FoxJob newJob;
    newJob.setFileName(filename);
    newJob.setM_ID(id);
    newJob.setNbDone(0);
    newJob.setNbRuns((int)nbRun);
    newJob.setNbTrial(nbOfTrial);
    newJob.setStatus(0);
    newJob.setName(name);
    newJob.setRand(rand);

    //if(m_WXFoxServerDataMutex->Lock()!=wxMUTEX_NO_ERROR) return;
    //save it and save it in FoxServer
    //m_jobs.push_back(newJob);
    m_FoxServer->AddJobToList(newJob);
    //m_WXFoxServerDataMutex->Unlock();

    //m_UpdateTimer->Stop();
    //Update GUI
    wxTimerEvent evt;
    evt.SetId(ID_UPDATE_TIMER);
    this->UpdateLists(evt);

    //int interval = m_UpdateTimer->GetInterval();
    //m_UpdateTimer->Start(interval);
}
void WXFoxServer::OnLoadJob(wxCommandEvent& event)
{
    wxFileDialog *dlg;
    dlg = new wxFileDialog(this,_T("Choose File :"),
                             _T(""),_T(""),_T("FOX files (*.xml,*.xmlgz)|*.xml;*.xmlgz;*.gz"),wxFD_OPEN | wxFD_FILE_MUST_EXIST);
    if(dlg->ShowModal() != wxID_OK) return;
    wxString path=dlg->GetPath();

    int newID;
    wxString Name=wxFileName::FileName(path).GetName();
    long nbOfTrial;
    long nbRun;
    wxString filename;
    bool randomize=true;

    if(path.Mid(path.size()-4)==wxString(_T(".xml")) && isFileFoxJob(path, Name, newID, nbOfTrial, nbRun, randomize)) {
        if(!isJobLoaded(newID)) {
            do{//show setup window
                if(!ShowSetJobWindow(newID, Name, nbOfTrial, nbRun, randomize)) return;
                cout<<nbRun<<","<<nbOfTrial<<","<<Name<<endl;
            }while((nbRun<=0)||(nbOfTrial<=0)||(Name==_T("")));
            filename.Printf(_T("JOB_%d.xml"), newID);
            #ifdef WIN32
            filename = m_working_dir + _T("\\") + filename;
            if(!wxFileExists(filename)) wxCopyFile(path, filename);
            #else
             filename = m_working_dir + _T("/") + filename;
            if(!wxFileExists(filename)) wxCopyFile(path, filename);
            #endif
            ChangeJobHeader(path, newID, Name, nbOfTrial, nbRun, randomize);
            AddServerJob(filename, Name, newID, nbOfTrial, nbRun, randomize);
        } else {
            wxMessageBox(_T("Job was probably loaded. Change job ID to load this job"), _T("Notice"), wxOK, this);
        }
    } else if(path.Mid(path.size()-4)==wxString(_T(".xml"))) {

        newID = GenerateJobID();
        nbOfTrial = 1000000;
        nbRun = 10;
        do{//show setup window
            if(!ShowSetJobWindow(newID, Name, nbOfTrial, nbRun, randomize)) return;
            cout<<nbRun<<","<<nbOfTrial<<","<<Name<<endl;
        }while((nbRun<=0)||(nbOfTrial<=0)||(Name==_T("")));

        filename.Printf(_T("JOB_%d.xml"), newID);
        #ifdef WIN32
        filename = m_working_dir + _T("\\") + filename;
        wxCopyFile(path, filename);
        #else
        filename = m_working_dir + _T("/") + filename;
        wxCopyFile(path, filename);
        #endif
        saveJobHeader(filename, newID, Name, nbOfTrial, nbRun, randomize);
        AddServerJob(filename, Name, newID, nbOfTrial, nbRun, randomize);
    }
    else
    {
        bool gz = false;
        if(path.size()>7)
           if(path.Mid(path.size()-6)==wxString(_T(".xml.gz"))) gz=true;
        if(path.size()>6)
           if(path.Mid(path.size()-6)==wxString(_T(".xmlgz"))) gz=true;
        if(gz)
        {
           wxFileInputStream is(path);
           wxZlibInputStream zstream(is);
           wxStringOutputStream wxos;
           zstream.Read(wxos);

           newID = GenerateJobID();
           nbOfTrial = 1000000;
           nbRun = 10;

           do{//show setup window
               if(!ShowSetJobWindow(newID, Name, nbOfTrial, nbRun, randomize)) return;
               cout<<nbRun<<","<<nbOfTrial<<","<<Name<<endl;
           }while((nbRun<=0)||(nbOfTrial<=0)||(Name==_T("")));

           filename.Printf(_T("JOB_%d.xml"), newID);
           #ifdef WIN32
           filename = m_working_dir + _T("\\") + filename;
           #else
           filename = m_working_dir + _T("/") + filename;
           #endif
           SaveDataAsFile(wxos.GetString(), filename);
           saveJobHeader(filename, newID, Name, nbOfTrial, nbRun, randomize);
           AddServerJob(filename, Name, newID, nbOfTrial, nbRun, randomize);
        }
    }
    dlg->Destroy();
}
bool WXFoxServer::ShowSetJobWindow(int ID, wxString &name, long &trials, long &runs, bool &rand)
{
   wxString tmp;
   tmp.Printf(_T("%d"), ID);
   wxDialog *InfoWindow = new wxDialog(NULL, -1,_T("Set job parameters"),wxDefaultPosition, wxDefaultSize, wxCAPTION);
   wxBoxSizer *sizer = new wxBoxSizer( wxVERTICAL );

   wxStaticText *label1 = new wxStaticText(InfoWindow, NULL, _T("Job ID:"));
   sizer->Add(label1, 0, wxALL|wxALIGN_LEFT|wxALIGN_TOP ,3);
   wxTextCtrl *IDText = new wxTextCtrl(InfoWindow, NULL, tmp, wxDefaultPosition, wxDefaultSize, wxTE_READONLY);
   sizer->Add(IDText, 0, wxALL|wxALIGN_LEFT|wxALIGN_TOP ,3);

   wxStaticText *label1a = new wxStaticText(InfoWindow, NULL, _T("Job Name"));
   sizer->Add(label1a, 0, wxALL|wxALIGN_LEFT|wxALIGN_TOP ,3);
   wxTextCtrl *NameText = new wxTextCtrl(InfoWindow, NULL, name, wxDefaultPosition, wxDefaultSize);
   sizer->Add(NameText, 0, wxALL|wxALIGN_LEFT|wxALIGN_TOP ,3);

   wxStaticText *label2 = new wxStaticText(InfoWindow, NULL, _T("Number of trials:"));
   sizer->Add(label2, 0, wxALL|wxALIGN_LEFT|wxALIGN_TOP ,3);
   tmp.Printf(_T("%ld"), trials);
   wxTextCtrl *TrialText = new wxTextCtrl(InfoWindow, NULL, tmp, wxDefaultPosition, wxDefaultSize);
   sizer->Add(TrialText, 0, wxALL|wxALIGN_LEFT|wxALIGN_TOP ,3);

   wxStaticText *label3 = new wxStaticText(InfoWindow, NULL, _T("Number of runs:"));
   sizer->Add(label3, 0, wxALL|wxALIGN_LEFT|wxALIGN_TOP ,3);
   tmp.Printf(_T("%ld"), runs);
   wxTextCtrl *RunText = new wxTextCtrl(InfoWindow, NULL, tmp, wxDefaultPosition, wxDefaultSize);
   sizer->Add(RunText, 0, wxALL|wxALIGN_LEFT|wxALIGN_TOP ,3);

   wxCheckBox *xrand = new wxCheckBox(InfoWindow, NULL, _T( "randomize"), wxDefaultPosition, wxDefaultSize);
   xrand->SetValue(rand);
   sizer->Add(xrand, 0, wxALL|wxALIGN_LEFT|wxALIGN_TOP ,3);

   wxButton *OKButton = new wxButton(InfoWindow, wxID_OK, _T("Save"));
   wxButton *NOButton = new wxButton(InfoWindow, wxID_CANCEL, _T("Cancel"));
   sizer->Add(OKButton, 0, wxALL|wxALIGN_LEFT|wxALIGN_TOP ,3);
   sizer->Add(NOButton, 0, wxALL|wxALIGN_LEFT|wxALIGN_TOP ,3);

   InfoWindow->SetAutoLayout( TRUE );
   InfoWindow->SetSizer( sizer );
   sizer->SetSizeHints( InfoWindow );
   sizer->Fit( InfoWindow );
   if(InfoWindow->ShowModal() == wxID_OK)
   {
      wxString Trials = TrialText->GetValue();
      wxString Runs = RunText->GetValue();
      name = NameText->GetValue();
      rand = xrand->GetValue();
      Trials.ToLong((long *) &trials);
      Runs.ToLong((long *) &runs);
      InfoWindow->Destroy();
      return true;
   }
   InfoWindow->Destroy();
   return false;
}
void WXFoxServer::EditJob()
{
   int nb = m_JobListTable->GetSelectedRows().Count();
   if(nb!=1) return;

   int r = m_JobListTable->GetSelectedRows().Item(0);

   if(r>= m_jobs.size()) {
      return;
   }
   
   int ID        = m_jobs[r].getM_ID();
   long trials   = m_jobs[r].getNbTrial();
   long runs     = m_jobs[r].getNbRuns();
   wxString name = m_jobs[r].getName();
   wxString filename = m_jobs[r].getFileName();
   bool randomize = m_jobs[r].randomize();
   do{//show setup window
      if(!ShowSetJobWindow(ID, name, trials, runs, randomize)) {
         return;
      }
   }while((runs<=0)||(trials<=0)||(name==_T("")));

   ChangeJobHeader(filename, ID, name, trials, runs, randomize);
   FoxJob newjob;
   newjob.setNbTrial(trials);
   newjob.setNbRuns(runs);
   newjob.setName(name);
   newjob.setRand(randomize);

   m_FoxServer->UpdateJob(r,&newjob);

}
void WXFoxServer::OnEditJob(wxCommandEvent& event)
{
   m_UpdateTimer->Stop();

   this->EditJob();

   //Update GUI List
   wxTimerEvent evt;
   evt.SetId(ID_UPDATE_TIMER);
   this->UpdateLists(evt);

   //start timer for updating gui...
   int interval = m_UpdateTimer->GetInterval();
   m_UpdateTimer->Start(interval);
}
void WXFoxServer::OnDeleteJob(wxCommandEvent& event)
{
   int nb = m_JobListTable->GetSelectedRows().Count();
   if(nb!=1) return;

   //not sure what happens woth the timer events if the modal window is open for some time...
   //better to stop the timer
   m_UpdateTimer->Stop();

   wxMessageDialog d(this,_T("Do you really want to delete selected job?"), _T("Alert"), wxYES | wxNO);
   if(wxID_YES!=d.ShowModal()) {       
       int interval = m_UpdateTimer->GetInterval();
       m_UpdateTimer->Start(interval);
       return;
   }

   int r = m_JobListTable->GetSelectedRows().Item(0);

   int er = m_FoxServer->DeleteJob(r);
   if(er==1)  wxMessageBox(_T("Job was deleted"), _T("Notice"), wxOK, this);
   if(er==0)  wxMessageBox(_T("Can't delete this job.\n Job was sent to clients or job finished. Only 'nbRuns' was changed."), _T("Notice"), wxOK, this);
   if(er==-1) wxMessageBox(_T("Can't delete this job! Job not found or server is busy."), _T("Notice"), wxOK, this);

   
   //Update Lists
   wxTimerEvent evt;
   evt.SetId(ID_UPDATE_TIMER);
   this->UpdateLists(evt);
   int interval = m_UpdateTimer->GetInterval();
   m_UpdateTimer->Start(interval);
}
void WXFoxServer::RunLocalClient(wxCommandEvent& event)
{
   m_UpdateTimer->Stop();

   int nCPU = wxThread::GetCPUCount();
   wxString nbCPUs;
   nbCPUs.Printf(_T("%d"), nCPU);

   wxString message;
   wxStandardPaths sp=wxStandardPaths::Get();

   message.Printf(_T("Would you like to run a client on this PC?\nSet the number of CPUs available for client or cancel this operation.\n%d CPUs has been detected on this PC.") , nCPU);
   wxTextEntryDialog dlg(m_parent, message, _T("Set a number of available CPUs"), nbCPUs, wxCANCEL | wxOK );
   if(wxID_OK==dlg.ShowModal()){
       nbCPUs = dlg.GetValue();
	   //wxString appname = wxApp::GetInstance()->argv[0];
       wxString appname = sp.GetExecutablePath();
       #ifdef WIN32
       wxString ClientDir = m_working_dir + _T("\\client");
       if(!wxDirExists(ClientDir)) wxMkdir(ClientDir);
       wxString cmd = appname + _T(" --runclient localhost --CPUs ") + nbCPUs + _T(" --working_dir ") + ClientDir;
       wxExecute(cmd);
       #else
       //if(appname(0,1)!=_T("/")) appname=wxGetCwd()+_T("/")+appname;
       //wxExecute(appname+_T(" --runclient localhost --CPUs ") + nbCPUs);
       wxString ClientDir = m_working_dir + _T("/client");
       if(!wxDirExists(ClientDir)) wxMkdir(ClientDir);
       long result= wxExecute(appname+_T(" --runclient localhost --CPUs ") + nbCPUs + _T(" --working_dir ") + ClientDir);
       //if(result==0) result=wxExecute(wxGetCwd()+_T("/")+appname+_T(" --runclient localhost --CPUs ") + nbCPUs);
       //if(result==0) result=wxExecute(_T("/usr/bin/")+appname+_T(" --runclient localhost --CPUs ") + nbCPUs);
       //if(result==0) result=wxExecute(_T("/usr/local/bin/")+appname+_T(" --runclient localhost --CPUs ") + nbCPUs);
       //if(appname(0,1)!=_T("/")) appname=wxGetCwd()+_T("/")+appname;
       //wxExecute(appname+_T(" --runclient localhost --CPUs ") + nbCPUs);
       #endif
   }

   int interval = m_UpdateTimer->GetInterval();
   m_UpdateTimer->Start(interval);
}
void WXFoxServer::UpdateJobList()
{
   int nbRow = m_JobListTable->GetRows();
   if(nbRow>0) m_JobListTable->DeleteRows(0, nbRow, true);
   wxString tmp;

   for(int i=0;i<m_jobs.size();i++){
      m_JobListTable->InsertRows(i,1,false);

      tmp.Printf(_T("%d"),m_jobs[i].getM_ID());
      m_JobListTable->SetCellValue(i,0,tmp);//ID
      m_JobListTable->SetReadOnly(i,0);
      m_JobListTable->SetCellValue(i,1,m_jobs[i].getName());//name
      m_JobListTable->SetReadOnly(i,1);
      tmp.Printf(_T("%ld"), m_jobs[i].getNbTrial());
      m_JobListTable->SetCellValue(i,2,tmp);//nbtrial
      m_JobListTable->SetReadOnly(i,1);
      //rand
      if(m_jobs[i].randomize()) {
          tmp.Printf(_T("YES"), m_jobs[i].getNbTrial());//???? :TODO:
      } else {
          tmp.Printf(_T("NO"), m_jobs[i].getNbTrial());//???? :TODO:
      }
      m_JobListTable->SetCellValue(i,3,tmp);//nbtrial
      m_JobListTable->SetReadOnly(i,1);

      tmp.Printf(_T("%d/%d/%d"), m_jobs[i].getNbRuns(), m_jobs[i].getNbDone(), m_jobs[i].GetSolvingNb());
      m_JobListTable->SetCellValue(i,4,tmp);//runs
      m_JobListTable->SetReadOnly(i,2);
      wxColor ccolor(255, 255, 255);        
      if(m_jobs[i].getNbDone() >= m_jobs[i].getNbRuns()) {
          ccolor.Set(200, 255, 200);
      } else if(m_jobs[i].GetSolvingNb()>0) {
          ccolor.Set(255, 200, 200);
      } else if(m_jobs[i].getNbDone()>0) {
          ccolor.Set(255, 255, 200);
      } 
      m_JobListTable->SetCellBackgroundColour(ccolor, i, 0);
      m_JobListTable->SetCellBackgroundColour(ccolor, i, 1);
      m_JobListTable->SetCellBackgroundColour(ccolor, i, 2);
      m_JobListTable->SetCellBackgroundColour(ccolor, i, 3);
      m_JobListTable->SetCellBackgroundColour(ccolor, i, 4);

   }
}
void WXFoxServer::UpdateResultList()
{
   int nb = m_ResultTable->GetRows(); 
   std::sort(m_results.begin(), m_results.end());
   
   wxString tmp;
   //rewrite already present rows (because we sorted/changed results before) and insert new ones.
   for(int i=0;i<m_results.size();i++) {
      if(i>=nb) {
          m_ResultTable->InsertRows(i, 1, false);
      }
      tmp.Printf(_T("%d"), m_results[i].order);
      m_ResultTable->SetCellValue(i, 0, tmp);//Nb
      tmp.Printf(_T("%d"), m_results[i].JobID);
      m_ResultTable->SetCellValue(i, 1, tmp);//Job ID
      tmp.Printf(_T("%0.2f"), m_results[i].Cost);
      m_ResultTable->SetCellValue(i, 2, tmp);//Cost
      //if(m_results[i].Show)   m_ResultTable->SetCellValue(i, 3, _T("TRUE"));//Show
      //else m_ResultTable->SetCellValue(i, 3, _T("FALSE"));//Show
   }
}
void WXFoxServer::UpdateClientList()
{
       //update client list
   int nbRow = m_ClientTable->GetRows();
   if(nbRow>0) m_ClientTable->DeleteRows(0, nbRow, true);
   for(int i=0;i<m_clients.size();i++){
        m_ClientTable->InsertRows(i,1,false);
        wxString tmp;

        //Name
        m_ClientTable->SetCellValue(i,0,m_clients[i].name);
        m_ClientTable->SetReadOnly(i,0);

        //ID
        tmp.Printf(_T("%d"), (int) m_clients[i].id);
        m_ClientTable->SetCellValue(i,1,tmp);
        m_ClientTable->SetReadOnly(i,1);

        //CPUs
        tmp.Printf(_T("%d/%d"), (int) m_clients[i].allCPUs, (int) (m_clients[i].allCPUs-m_clients[i].availCPUs));
        m_ClientTable->SetCellValue(i,2,tmp);
        m_ClientTable->SetReadOnly(i,2);        
        //status
        //m_ClientTable->SetCellValue(i,3,clients[i].status);
        //m_ClientTable->SetReadOnly(i,3);
        //set colors
        if(m_clients[i].availCPUs < m_clients[i].allCPUs) {
            m_ClientTable->SetCellBackgroundColour(wxColor(255, 200, 200), i, 0);
            m_ClientTable->SetCellBackgroundColour(wxColor(255, 200, 200), i, 1);
            m_ClientTable->SetCellBackgroundColour(wxColor(255, 200, 200), i, 2);
        } else {
            m_ClientTable->SetCellBackgroundColour(wxColor(200, 255, 200), i, 0);
            m_ClientTable->SetCellBackgroundColour(wxColor(200, 255, 200), i, 1);
            m_ClientTable->SetCellBackgroundColour(wxColor(200, 255, 200), i, 2);
        }
   }

}
void WXFoxServer::UpdateLists(wxTimerEvent& event)
{
   std::vector<GridClient>      newclients;
   std::vector<FoxJob >         newjobs;
   if(m_FoxServer->IsServerRunning()){
      m_FoxServer->GetData(newclients, m_results, newjobs);
   }


   //check if there is something new in client list
   bool update = false;
   if(m_clients.size() == newclients.size()) {
       for(int i=0;i<m_clients.size();i++) {
           if(m_clients[i].name.compare(newclients[i].name)!=0) update = true;
           if(m_clients[i].status.compare(newclients[i].status)!=0) update = true;
           if(m_clients[i].id != newclients[i].id) update = true;
           if(m_clients[i].availCPUs != newclients[i].availCPUs) update = true;
           if(m_clients[i].allCPUs != newclients[i].allCPUs) update = true;
       }
   } else {
       update = true;
   }
   m_clients = newclients;
   //update client list...
   if(update) {
       UpdateClientList();
   }
   //update job list only when something is new there...
   update = false;
   if(m_jobs.size() == newjobs.size()) {
       for(int i=0;i<m_jobs.size();i++) {
           if(m_jobs[i].getName().compare(newjobs[i].getName())!=0) update = true;
           if(m_jobs[i].getNbDone() != newjobs[i].getNbDone()) update = true;
           if(m_jobs[i].getNbRuns() != newjobs[i].getNbRuns()) update = true;
           if(m_jobs[i].GetSolvingNb() != newjobs[i].GetSolvingNb()) update = true;
           if(m_jobs[i].getNbTrial() != newjobs[i].getNbTrial()) update = true;
           if(m_jobs[i].randomize() != newjobs[i].randomize()) update = true;
       }
   } else {
       update = true;
   }
   m_jobs = newjobs;
   if(update) {
       UpdateJobList();
   }

   //there is a check whethear to update or not in the function...
   UpdateResultList();
}
void WXFoxServer::OnGridJobClick(wxGridEvent &event)
{
   int r = event.GetRow();
   int c = event.GetCol();

   if(m_JobListTable->GetSelectedRows().Count()==1){
      m_JobListTable->DeselectRow(m_JobListTable->GetSelectedRows().Item(0));
   }

   m_JobListTable->SelectRow(r, true);
}
void WXFoxServer::OnGridResultClick(wxGridEvent &event)
{
   int r = event.GetRow();
   int c = event.GetCol();

   //if(m_WXFoxServerDataMutex->Lock()!=wxMUTEX_NO_ERROR) return;

   if(m_ResultTable->GetSelectedRows().Count()==1){
      m_ResultTable->DeselectRow(m_ResultTable->GetSelectedRows().Item(0));
   }

   m_ResultTable->SelectRow(r, true);
   //m_WXFoxServerDataMutex->Unlock();

}
void WXFoxServer::OnShowResultsServer(wxCommandEvent& event)
{
   int nb = m_ResultTable->GetSelectedRows().Count();
   #ifdef __DEBUG__
   (*fpObjCrystInformUser)(wxString::Format("Show Results: %d rows selected",nb).ToStdString());
   #endif
   if(nb!=1)
   {
      wxMessageBox(_T("Select one result!"), _T("Error"), wxOK, this);
      return;
   }

   int r = m_ResultTable->GetSelectedRows().Item(0);
   #ifdef __DEBUG__
   (*fpObjCrystInformUser)(wxString::Format("Show Results: result %d/%zd selected",r,m_results.size()).ToStdString());
   #endif
   if(r>=m_results.size()) return;

   wxString file = m_results[r].filename;
   (*fpObjCrystInformUser)(wxString::Format("Show Results: opening file: "+file).ToStdString());

/*
   if(safe)
   {
      bool saved=true;
      for(int i=0;i<gRefinableObjRegistry.GetNb();i++)
         if(gRefinableObjRegistry.GetObj(i).GetClockMaster()>mClockLastSave)
         {
            saved=false;
            break;
         }
      if(!saved)
      {
         wxString msg;
         msg.Printf( _T("Some objects have not been saved\n")
                  _T("Do you really want to close all ?"));

         wxMessageDialog d(this,msg, _T("Really Close ?"), wxYES | wxNO);
         if(wxID_YES!=d.ShowModal()) return;
      }
   }

   */
   if(m_dataLoaded) {
       gOptimizationObjRegistry.DeleteAll();
       gDiffractionDataSingleCrystalRegistry.DeleteAll();
       gPowderPatternRegistry.DeleteAll();
       gCrystalRegistry.DeleteAll();
       (*fpObjCrystInformUser)("");
       (*fpObjCrystInformUser)("Closed all objects");
       (*fpObjCrystInformUser)("");
       m_dataLoaded = false;
   }

   if(file.size()>4 && file.Mid(file.size()-4)==wxString(_T(".xml"))) {    
        wxFileInputStream is(file);
        stringstream in;
        if(is.GetSize()>0)
        {
          char * tmpbuf=new char[is.GetSize()];
          is.Read(tmpbuf,is.GetSize());
          in.write(tmpbuf,is.GetSize());
          delete[] tmpbuf;
        }
        else while (!is.Eof()) in<<(char)is.GetC();
        try{XMLCrystFileLoadAllObject(in);}
        catch(const ObjCrystException &except)
        {
          wxMessageDialog d(this,_T("Failed loading file1:\n")+file,_T("Error loading file"),wxOK|wxICON_ERROR);
          d.ShowModal();
          return;
        };
        m_dataLoaded = true;
   }

}
void WXFoxServer::OnShowResults(wxCommandEvent& event)
{
   //if(m_WXFoxServerDataMutex->Lock()!=wxMUTEX_NO_ERROR) return;
   int nb = m_ResultTable->GetSelectedRows().Count();
   #ifdef __DEBUG__
   (*fpObjCrystInformUser)(wxString::Format("Show Results: %d rows selected",nb).ToStdString());
   #endif
   if(nb!=1)
   {
      wxMessageBox(_T("Select one result!"), _T("Error"), wxOK, this);
      return;
   }

   int r = m_ResultTable->GetSelectedRows().Item(0);
   #ifdef __DEBUG__
   (*fpObjCrystInformUser)(wxString::Format("Show Results: result %d/%zd selected",r,m_results.size()).ToStdString());
   #endif
   if(r>=m_results.size()) return;

   wxString file = m_results[r].filename;
   (*fpObjCrystInformUser)(wxString::Format("Show Results: opening file: "+file).ToStdString());

   //m_WXFoxServerDataMutex->Unlock();


   wxString cmd;
   #ifdef WIN32
   cmd = wxApp::GetInstance()->argv[0];
   cmd +=_T(" ");
   cmd += file;
   wxExecute(cmd);
   #else
   wxString appname = wxStandardPaths::Get().GetExecutablePath();
   wxString com=appname+_T(" ")+file;
   long result= wxExecute(com);
   (*fpObjCrystInformUser)(com.ToStdString());
   //if(result==0) result=wxExecute(wxGetCwd()+_T("/")+appname+_T(" ")+file);
   //if(result==0) result=wxExecute(_T("/usr/bin/")+appname+_T(" ")+file);
   //if(result==0) result=wxExecute(_T("/usr/local/bin/")+appname+_T(" ")+file);
   #endif
}
