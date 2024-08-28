
#include <wx/textdlg.h>
#include "wx/filename.h"
#include "wx/sstream.h"

#ifdef __WX__CRYST__
   #include "ObjCryst/wxCryst/wxCrystal.h"
#endif

// Fox client/server grid


using namespace ObjCryst;
using namespace std;

#include "wxFoxMaster.h"
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


BEGIN_EVENT_TABLE(WXFoxMaster, wxWindow)
  //EVT_BUTTON(STOP_TO_CLIENT,                  WXFoxMaster::StopAllClients)
  EVT_BUTTON(RUN_LOCAL_CLIENT,                 WXFoxMaster::RunLocalClient)
  //EVT_BUTTON(RUN_TO_CLIENT,                 WXFoxMaster::RunALLClient)
  EVT_BUTTON(NEW_JOB,                        WXFoxMaster::OnNewJob)
  EVT_BUTTON(EDIT_JOB,                        WXFoxMaster::OnEditJob)
  EVT_BUTTON(DELETE_JOB,                     WXFoxMaster::OnDeleteJob)
  EVT_BUTTON(LOAD_JOB,                        WXFoxMaster::OnLoadJob)
  EVT_TIMER(ID_UPDATE_TIMER,                  WXFoxMaster::UpdateLists)
  EVT_GRID_CMD_CELL_LEFT_CLICK(GRID_RESULT_LIST,   WXFoxMaster::OnGridResultClick)
  EVT_GRID_CMD_CELL_LEFT_CLICK(GRID_JOB_LIST,      WXFoxMaster::OnGridJobClick)
  EVT_BUTTON(SHOW_RESULTS,                     WXFoxMaster::OnShowResults)
  EVT_BUTTON(SHOW_RESULTS_SERVER,                     WXFoxMaster::OnShowResultsServer)

END_EVENT_TABLE()

WXFoxMaster::WXFoxMaster(wxWindow* parent, wxString workingDir):
wxWindow(parent,-1),m_parent(parent)
{
   m_dataLoaded = false;
   m_working_dir = workingDir;
   m_grid_master = NULL;
   m_ServerStatus = NULL;
   InitServer();
}
WXFoxMaster::~WXFoxMaster(void)
{

}
void WXFoxMaster::Clear()
{
   m_UpdateTimer->Stop();
   delete m_UpdateTimer;
   //todo: clear m_jobs, m_results
   //delete m_WXFoxMasterDataMutex;
   //delete m_FoxServer;
}
void WXFoxMaster::InitServer()
{
   
   //starting server   
   m_grid_master = new FoxGridMaster(m_working_dir);
   m_grid_master->InitializeCommunication();   

   //Start update timer
   m_UpdateTimer = new wxTimer(this, ID_UPDATE_TIMER);
   m_UpdateTimer->Start(3*1000, true);

   //display all necessary controls
   wxBoxSizer *topSizer = new wxBoxSizer( wxVERTICAL );

   unsigned int xsize=650;
   //Job List table
   wxBoxSizer *JobSizer = new wxBoxSizer( wxVERTICAL );

   m_ServerStatus = new wxStaticText(this, NULL, _T("Server status: "), wxDefaultPosition, wxDefaultSize, 0 , _T("server status"));
   JobSizer->Add(m_ServerStatus,0, wxALL|wxALIGN_LEFT|wxALIGN_TOP ,3);
   wxStaticText *JobLabel = new wxStaticText(this, NULL, _T("Job List: "), wxDefaultPosition, wxDefaultSize, 0 , _T("label"));
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
   m_JobListTable->SetColSize(0, (xsize-10)/5);
   m_JobListTable->SetColSize(1, (xsize-10)/5+20);
   m_JobListTable->SetColSize(2, (xsize-10)/5-20);
   m_JobListTable->SetColSize(3, (xsize-10)/5-10);
   m_JobListTable->SetColSize(4, (xsize-10)/5+10);
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


   //client list label
   wxBoxSizer *listClientSizer = new wxBoxSizer( wxVERTICAL );
   wxStaticText *label2 = new wxStaticText(this, NULL, _T("Client list: "), wxDefaultPosition, wxDefaultSize, 0 , _T("label"));
   listClientSizer->Add(label2,0, wxALL|wxALIGN_LEFT|wxALIGN_TOP ,3);

   xsize = xsize/2-5;
   //client list table
   m_ClientTable = new wxGrid(this, NULL, wxDefaultPosition, wxSize(xsize,150), wxWANTS_CHARS, _T("m_ClientTable"));
   m_ClientTable->CreateGrid(1,2,wxGrid::wxGridSelectRows);
   m_ClientTable->SetColLabelValue(0, _T("Client IP"));   
   m_ClientTable->SetColLabelValue(1, _T("CPUs/using"));
   //m_ClientTable->SetColLabelValue(3, _T("Status"));
   m_ClientTable->SetColLabelSize(20);
   m_ClientTable->SetRowLabelSize(0);
   m_ClientTable->SetColSize(0, xsize/2);
   m_ClientTable->SetColSize(1, xsize/2);     
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
   m_ResultTable->SetColLabelValue(0, _T("Job ID"));
   m_ResultTable->SetColLabelValue(1, _T("Cost"));
   m_ResultTable->SetColLabelValue(2, _T("Rwp"));
   //m_ResultTable->SetColLabelValue(3, _T("Show"));
   m_ResultTable->SetColLabelSize(20);
   m_ResultTable->SetRowLabelSize(0);
   m_ResultTable->SetColSize(0, 50);
   m_ResultTable->SetColSize(1, (xsize-50)/2);
   m_ResultTable->SetColSize(2, (xsize-50)/2);
   //m_ResultTable->SetColSize(3, 50);
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

   
   topSizer->Add(JobSizer, 0, wxALL|wxALIGN_TOP);
   topSizer->Add(listsSizer, 0, wxALL|wxALIGN_TOP);
   
   

   this->SetSizer( topSizer );
   wxTheApp->GetTopWindow()->Layout();
   wxTheApp->GetTopWindow()->SendSizeEvent();   
}
void WXFoxMaster::saveJobHeader(wxString filename, long ID, wxString name, int nbOfTrial, int nbRun, bool rand, wxString &data) {    
    wxString file;
    if(this->LoadFile(filename, data)) {
        int r = (int) rand;
        file.Printf(_T("<FoxJob Name=\"%s\" ID=\"%d\" nbOfTrial=\"%ld\" nbRun=\"%d\" rand=\"%d\"></FoxJob>\n"), name.c_str(), ID, nbOfTrial, nbRun, r);
        file +=data;
        SaveDataAsFile(file, filename);
    }   
}
void WXFoxMaster::ChangeJobHeader(wxString filename, long ID, wxString name, int nbOfTrial, int nbRun, bool rand, wxString &data) 
{
    wxString file;    
    if(LoadFile(filename, data)) {
        int r = (int) rand;
        int pos = data.First(_T("</FoxJob>"));
        if(pos!=-1) {
            pos +=9;
            data = data.Mid(pos);
        }
        file.Printf(_T("<FoxJob Name=\"%s\" ID=\"%d\" nbOfTrial=\"%ld\" nbRun=\"%d\" rand=\"%d\"></FoxJob>\n"), name.c_str(), ID, nbOfTrial, nbRun, r);
        file += data;
        SaveDataAsFile(file, filename);
    }
}
bool WXFoxMaster::LoadFile(wxString filename, wxString &in)
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
void WXFoxMaster::SaveDataAsFile(wxString out, wxString filename)
{
   wxFile outFile(filename, wxFile::write);
   if(outFile.IsOpened())
   {
      outFile.Write(out);
      outFile.Close();
   }
}
void WXFoxMaster::OnNewJob(wxCommandEvent& event)
{
    //if(!m_dataLoaded)
    if(gOptimizationObjRegistry.GetNb()==0)
    {
        wxMessageBox(_T("You need at least one optimization object !"), _T("No data to optimize"), wxOK | wxICON_INFORMATION, this);
        return;
    }

    WXCrystValidateAllUserInput();

    long newID = m_grid_master->generateJobID();
    wxString Name = _T("jobname");
    int nbOfTrial = 1000000;
    int nbRun = 10;
    bool randomize = true;

    do {//show setup window
        if(!ShowEditJobWindow(newID, Name, nbOfTrial, nbRun, randomize, false)) return;
    } while((nbRun<=0)||(nbOfTrial<=0)||(Name==_T("")));
   
    string tmp;
    stringstream str;
    #ifdef WIN32
    str<<m_working_dir<<"\\"<<"JOB_"<<newID<<".xml";
    #else
    str<<m_working_dir<<"/"<<"JOB_"<<newID<<".xml";
    #endif
    tmp = str.str();

    XMLCrystFileSaveGlobal(tmp);
    wxString jobdata;
    saveJobHeader(wxString::FromAscii(tmp.c_str()), newID, Name, nbOfTrial, nbRun, randomize, jobdata);
    AddServerJob(wxString::FromAscii(tmp.c_str()), Name, jobdata, newID, nbOfTrial, nbRun, randomize);

    UpdateJobList();
}
bool WXFoxMaster::isFileFoxJob(wxString path, wxString &name, long &id, int &nbOfTrial, int &nbRun, bool &rand) {
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
           ID.ToLong(&id);
        }
        if(tag.GetAttributeName(i)=="nbOfTrial"){
           Tr = wxString::FromAscii(tag.GetAttributeValue(i).c_str());
           Tr.ToInt(&nbOfTrial);
        }
        if(tag.GetAttributeName(i)=="nbRun"){
           Run = wxString::FromAscii(tag.GetAttributeValue(i).c_str());
           Run.ToInt(&nbRun);
        }
        if(tag.GetAttributeName(i)=="rand"){
           long r;
           Rand = wxString::FromAscii(tag.GetAttributeValue(i).c_str());
           Rand.ToLong(&r);
           rand = (bool) r;
        }
    }
    is.close();
    return true;
}
bool WXFoxMaster::isJobLoaded(long ID) 
{
    return m_grid_master->existsJob(ID);    
}
void WXFoxMaster::AddServerJob(wxString filename, wxString name, wxString data, long id, int nbOfTrial, int nbRun, bool rand) 
{
    m_grid_master->addJob(filename, name, data, id, nbOfTrial, nbRun, rand);   
}
void WXFoxMaster::loadMultipleJobs(wxArrayString jobfiles)
{
    for (size_t i = 0; i < jobfiles.GetCount(); i++) {
        long newID;
        wxString Name=wxFileName::FileName(jobfiles[i]).GetName();
        int nbOfTrial;
        int nbRun;
        wxString filename;
        bool randomize=true;
        wxString jobdata;
        if(jobfiles[i].Mid(jobfiles[i].size()-4)==wxString(_T(".xml")) && isFileFoxJob(jobfiles[i], Name, newID, nbOfTrial, nbRun, randomize)) {
            if(!isJobLoaded(newID)) {                
                filename.Printf(_T("JOB_%d.xml"), newID);
                #ifdef WIN32
                filename = m_working_dir + _T("\\") + filename;
                if(!wxFileExists(filename)) wxCopyFile(jobfiles[i], filename);
                #else
                 filename = m_working_dir + _T("/") + filename;
                if(!wxFileExists(filename)) wxCopyFile(jobfiles[i], filename);
                #endif          
                ChangeJobHeader(filename, newID, Name, nbOfTrial, nbRun, randomize, jobdata);
                AddServerJob(filename, Name, jobdata, newID, nbOfTrial, nbRun, randomize);
            } else {
                wxMessageBox(_T("Job " + Name + "was probably loaded. Change job ID to load this job"), _T("Notice"), wxOK, this);
            }
        }
    }
}
void WXFoxMaster::OnLoadJob(wxCommandEvent& event)
{
    wxFileDialog *dlg;
    dlg = new wxFileDialog(this,_T("Choose File :"),
                             _T(""),_T(""),_T("FOX files (*.xml,*.xmlgz)|*.xml;*.xmlgz;*.gz"),wxFD_OPEN | wxFD_FILE_MUST_EXIST| wxFD_MULTIPLE);
    if(dlg->ShowModal() != wxID_OK) return;
    wxArrayString filepaths;
    dlg->GetPaths(filepaths);

    wxString path;
    if(filepaths.GetCount()==1) {
        path = filepaths[0];
    } else {
        loadMultipleJobs(filepaths);
        dlg->Destroy();
        UpdateJobList();
        return;
    }
    
    long newID;
    wxString Name=wxFileName::FileName(path).GetName();
    int nbOfTrial;
    int nbRun;
    wxString filename;
    bool randomize=true;
    wxString jobdata;

    if(path.Mid(path.size()-4)==wxString(_T(".xml")) && isFileFoxJob(path, Name, newID, nbOfTrial, nbRun, randomize)) {
        if(!isJobLoaded(newID)) {
            do {//show setup window
                if(!ShowEditJobWindow(newID, Name, nbOfTrial, nbRun, randomize, false)) return;
                VFN_DEBUG_MESSAGE("WXFoxMaster::OnLoadJob() nbrun="<<nbRun<<", nbtrial="<<nbOfTrial<<", name="<<Name,10);
            } while((nbRun<=0)||(nbOfTrial<=0)||(Name==_T("")));
            filename.Printf(_T("JOB_%d.xml"), newID);
            #ifdef WIN32
            filename = m_working_dir + _T("\\") + filename;
            if(!wxFileExists(filename)) wxCopyFile(path, filename);
            #else
             filename = m_working_dir + _T("/") + filename;
            if(!wxFileExists(filename)) wxCopyFile(path, filename);
            #endif
            ChangeJobHeader(filename, newID, Name, nbOfTrial, nbRun, randomize, jobdata);
            AddServerJob(filename, Name, jobdata, newID, nbOfTrial, nbRun, randomize);
        } else {
            wxMessageBox(_T("Job was probably loaded. Change job ID to load this job"), _T("Notice"), wxOK, this);
        }
    } else if(path.Mid(path.size()-4)==wxString(_T(".xml"))) {

        newID = m_grid_master->generateJobID();
        nbOfTrial = 1000000;
        nbRun = 10;
        do{//show setup window
            if(!ShowEditJobWindow(newID, Name, nbOfTrial, nbRun, randomize, false)) return;
            VFN_DEBUG_MESSAGE("WXFoxMaster::OnLoadJob() nbrun="<<nbRun<<", nbtrial="<<nbOfTrial<<", name="<<Name,10);
        }while((nbRun<=0)||(nbOfTrial<=0)||(Name==_T("")));

        filename.Printf(_T("JOB_%d.xml"), newID);
        #ifdef WIN32
        filename = m_working_dir + _T("\\") + filename;
        wxCopyFile(path, filename);
        #else
        filename = m_working_dir + _T("/") + filename;
        wxCopyFile(path, filename);
        #endif
        saveJobHeader(filename, newID, Name, nbOfTrial, nbRun, randomize, jobdata);
        AddServerJob(filename, Name, jobdata, newID, nbOfTrial, nbRun, randomize);
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

           newID = m_grid_master->generateJobID();
           nbOfTrial = 1000000;
           nbRun = 10;

           do{//show setup window
               if(!ShowEditJobWindow(newID, Name, nbOfTrial, nbRun, randomize, false)) return;
               VFN_DEBUG_MESSAGE("WXFoxMaster::OnLoadJob() nbrun="<<nbRun<<", nbtrial="<<nbOfTrial<<", name="<<Name,10);
           }while((nbRun<=0)||(nbOfTrial<=0)||(Name==_T("")));

           filename.Printf(_T("JOB_%d.xml"), newID);
           #ifdef WIN32
           filename = m_working_dir + _T("\\") + filename;
           #else
           filename = m_working_dir + _T("/") + filename;
           #endif
           SaveDataAsFile(wxos.GetString(), filename);
           saveJobHeader(filename, newID, Name, nbOfTrial, nbRun, randomize, jobdata);
           AddServerJob(filename, Name, jobdata, newID, nbOfTrial, nbRun, randomize);
        }
    }
    dlg->Destroy();

    UpdateJobList();
}
bool WXFoxMaster::ShowEditJobWindow(long ID, wxString &name, int &trials, int &runs, bool &rand, bool onlyRuns)
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

   if(onlyRuns) {
       NameText->SetEditable(false);
       TrialText->SetEditable(false);
       xrand->Enable(false);
   }

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
      Trials.ToInt(&trials);
      Runs.ToInt(&runs);
      InfoWindow->Destroy();
      VFN_DEBUG_MESSAGE("WXFoxMaster::ShowEditJobWindow:runs="<<runs<<", trials="<<trials<<", name="<<name, 10);
      return true;
   }
   InfoWindow->Destroy();
   return false;
}
void WXFoxMaster::EditJob()
{
    vector<MasterJob> jobs = m_grid_master->getJobs();


    int nb = m_JobListTable->GetSelectedRows().Count();
    if(nb!=1) return;
    int r = m_JobListTable->GetSelectedRows().Item(0);
    if(r>= jobs.size()) {
        return;
    }

    if(jobs[r].getNbFreeJobs() == jobs[r].nb_runs) {
        //the job was not processed yet, we can change it completely
        do {//show setup window
            if(!ShowEditJobWindow(jobs[r].ID, jobs[r].name, jobs[r].nb_trial, jobs[r].nb_runs, jobs[r].randomize, false)) {
                return;
            }
        } while((jobs[r].nb_runs<=0)||(jobs[r].nb_trial<=0)||(jobs[r].name==_T("")));
    } else {
        //the job is in process or done and we can change only the nb of runs
        do {//show setup window
            if(!ShowEditJobWindow(jobs[r].ID, jobs[r].name, jobs[r].nb_trial, jobs[r].nb_runs, jobs[r].randomize, true)) {
                return;
            }
        } while(jobs[r].nb_runs<=0);
    }

    ChangeJobHeader(jobs[r].filename, jobs[r].ID, jobs[r].name, jobs[r].nb_trial, jobs[r].nb_runs, jobs[r].randomize, jobs[r].data);

    m_grid_master->UpdateJobHeader(jobs[r]);

}
void WXFoxMaster::OnEditJob(wxCommandEvent& event)
{
   m_UpdateTimer->Stop();

   this->EditJob();

   //Update GUI List
   wxTimerEvent evt;
   evt.SetId(ID_UPDATE_TIMER);
   this->UpdateLists(evt);

   //start timer for updating gui...
   m_UpdateTimer->Start(1, true);
}
void WXFoxMaster::OnDeleteJob(wxCommandEvent& event)
{
   
   int nb = m_JobListTable->GetSelectedRows().Count();
   if(nb!=1) return;
   int r = m_JobListTable->GetSelectedRows().Item(0);

   vector<MasterJob> jobs = m_grid_master->getJobs();

   if(r>=jobs.size()) return;

   wxMessageDialog d(this,_T("Do you really want to delete selected job?"), _T("Alert"), wxYES | wxNO);
   if(wxID_YES!=d.ShowModal()) {       
       return;
   }

   if(jobs[r].getNbFreeJobs() == jobs[r].nb_runs) {
       m_grid_master->DeleteJob(jobs[r].ID);
       wxMessageBox(_T("Job was deleted"), _T("Notice"), wxOK, this);
   } else {
       jobs[r].nb_runs = jobs[r].nb_runs-jobs[r].getNbFreeJobs();
       ChangeJobHeader(jobs[r].filename, jobs[r].ID, jobs[r].name, jobs[r].nb_trial, jobs[r].nb_runs, jobs[r].randomize, jobs[r].data);
       m_grid_master->UpdateJobHeader(jobs[r]);
       wxMessageBox(_T("Can't delete this job.\n Job was sent to clients or job had already finished. Only 'nbRuns' was changed."), _T("Notice"), wxOK, this);
   }
   m_UpdateTimer->Start(1, true);
}
void WXFoxMaster::RunLocalClient(wxCommandEvent& event)
{
   m_UpdateTimer->Stop();

   int nCPU = wxThread::GetCPUCount();
   wxString nbCPUs;
   nbCPUs.Printf(_T("%d"), nCPU);

   wxString message;
   wxStandardPaths sp=wxStandardPaths::Get();

   message.Printf(_T("Would you also like to run client on this computer?\nSet the number of CPUs available for client or cancel this operation.\n%d CPUs have been detected.") , nCPU);
   wxTextEntryDialog dlg(m_parent, message, _T("Set a number of available CPUs"), nbCPUs, wxCANCEL | wxOK );
   if(wxID_OK==dlg.ShowModal()){
       nbCPUs = dlg.GetValue();
       int port = m_grid_master->getServersPort();

	   //wxString appname = wxApp::GetInstance()->argv[0];
       wxString appname = sp.GetExecutablePath();
       #ifdef WIN32
       wxString ClientDir = m_working_dir + _T("\\client");
       if(!wxDirExists(ClientDir)) wxMkdir(ClientDir);
       wxString cmd = "\"" + appname +"\"" + _T(" --runclient localhost:") + to_string(port) + _T(" --CPUs ") + nbCPUs + _T(" --working_dir ") + "\"" + ClientDir + "\"";
       wxExecute(cmd);
       #else
       //if(appname(0,1)!=_T("/")) appname=wxGetCwd()+_T("/")+appname;
       //wxExecute(appname+_T(" --runclient localhost --CPUs ") + nbCPUs);
       wxString ClientDir = m_working_dir + _T("/client");
       if(!wxDirExists(ClientDir)) wxMkdir(ClientDir);
       wxString cmd = appname+ _T(" --runclient localhost:") + to_string(port) + _T(" --CPUs ") + nbCPUs + _T(" --working_dir ") + ClientDir;
       VFN_DEBUG_MESSAGE("WXFoxMaster::RunLocalClient() command="<<cmd,10);
       long result= wxExecute(cmd);
       VFN_DEBUG_MESSAGE("WXFoxMaster::RunLocalClient() result="<<result,10);
       //if(result==0) result=wxExecute(wxGetCwd()+_T("/")+appname+_T(" --runclient localhost --CPUs ") + nbCPUs);
       //if(result==0) result=wxExecute(_T("/usr/bin/")+appname+_T(" --runclient localhost --CPUs ") + nbCPUs);
       //if(result==0) result=wxExecute(_T("/usr/local/bin/")+appname+_T(" --runclient localhost --CPUs ") + nbCPUs);
       //if(appname(0,1)!=_T("/")) appname=wxGetCwd()+_T("/")+appname;
       //wxExecute(appname+_T(" --runclient localhost --CPUs ") + nbCPUs);
       #endif
   }

   m_UpdateTimer->Start(1000, true);
}
void WXFoxMaster::UpdateJobList()
{    
    vector<MasterJob> jobs = m_grid_master->getJobs();

    int nbRow = m_JobListTable->GetNumberRows();
    if(nbRow < jobs.size()) {
        m_JobListTable->InsertRows(0,jobs.size()-nbRow,false);
    } else if(nbRow > jobs.size()) {
        m_JobListTable->DeleteRows(0, nbRow - jobs.size(), true);
    }

    //if(nbRow>0) m_JobListTable->DeleteRows(0, nbRow, true);
    wxString tmp;

    for(int i=0;i<jobs.size();i++){
        //m_JobListTable->InsertRows(i,1,false);

        tmp.Printf(_T("%d"),jobs[i].ID);
        m_JobListTable->SetCellValue(i,0,tmp);//ID
        m_JobListTable->SetReadOnly(i,0);
        m_JobListTable->SetCellValue(i,1,jobs[i].name);//name
        m_JobListTable->SetReadOnly(i,1);
        tmp.Printf(_T("%ld"), jobs[i].nb_trial);
        m_JobListTable->SetCellValue(i,2,tmp);//nbtrial
        m_JobListTable->SetReadOnly(i,1);
        //rand
        if(jobs[i].randomize) {
            tmp = _T("YES");
        } else {
            tmp = _T("NO");
        }
        m_JobListTable->SetCellValue(i,3,tmp);//nbtrial
        m_JobListTable->SetReadOnly(i,1);

        tmp.Printf(_T("%d/%d/%d"), jobs[i].nb_runs, jobs[i].results.size(), jobs[i].runningIPs.size());
        m_JobListTable->SetCellValue(i,4,tmp);//runs
        m_JobListTable->SetReadOnly(i,2);
        wxColor ccolor(255, 255, 255);    
        
        if(jobs[i].results.size() >= jobs[i].nb_runs) {
            ccolor.Set(200, 255, 200);
        } else if(jobs[i].runningIPs.size()>0) {
            ccolor.Set(255, 200, 200);
        } else if(jobs[i].results.size()>0) {
            ccolor.Set(255, 255, 200);
        } 
        
        m_JobListTable->SetCellBackgroundColour(i, 0, ccolor);
        m_JobListTable->SetCellBackgroundColour(i, 1, ccolor);
        m_JobListTable->SetCellBackgroundColour(i, 2, ccolor);
        m_JobListTable->SetCellBackgroundColour(i, 3, ccolor);
        m_JobListTable->SetCellBackgroundColour(i, 4, ccolor);

    }   
}
void WXFoxMaster::UpdateResultList()
{
    int nb = m_ResultTable->GetNumberRows();
    vector<MasterJob> jobs = m_grid_master->getJobs();
    int q = 0;
    for(int i=0;i<jobs.size();i++) {
        q+=jobs[i].results.size();
    }
    if(nb == q) return;
    
    m_Results.clear();
    for(int i=0;i<jobs.size();i++) {
        m_Results.insert(m_Results.begin(), jobs[i].results.begin(), jobs[i].results.end());
    }

    std::sort(m_Results.begin(), m_Results.end());
    wxString tmp;
    //rewrite already present rows (because we sorted/changed results before) and insert new ones.
    for(int i=0;i<m_Results.size();i++) {        
        if(i>=nb) {
            m_ResultTable->InsertRows(i, 1, false);
        }
        tmp.Printf(_T("%d"), m_Results[i].JobID);
        m_ResultTable->SetCellValue(i, 0, tmp);//Nb
        tmp.Printf(_T("%0.2f"), m_Results[i].Cost);
        m_ResultTable->SetCellValue(i, 1, tmp);//Job ID
        tmp.Printf(_T("%0.4f"), m_Results[i].Rwp);
        m_ResultTable->SetCellValue(i, 2, tmp);//Cost
    }
}
void WXFoxMaster::UpdateClientList(vector<FoxGridMaster::SLAVE_FOX_INFO> SI)
{
    
   //update client list
   int nbRow = m_ClientTable->GetNumberRows();
   if(nbRow>0) m_ClientTable->DeleteRows(0, nbRow, true);
   for(int i=0;i<SI.size();i++){
        m_ClientTable->InsertRows(i,1,false);
        wxString tmp;

        //Name
        m_ClientTable->SetCellValue(i,0,SI[i].sinfo.ip+":"+to_string(SI[i].sinfo.port));
        m_ClientTable->SetReadOnly(i,0);
        

        //CPUs
        tmp.Printf(_T("%d/%d"), (int) SI[i].nb_CPU_all, (int) (SI[i].nb_CPU_all-SI[i].nb_CPU_idle));
        m_ClientTable->SetCellValue(i,1,tmp);
        m_ClientTable->SetReadOnly(i,1);        

        //set colors
        wxColour bg = wxColour(255, 255, 255);
        if(SI[i].connected) {
            if((SI[i].nb_CPU_idle == 0) && (SI[i].nb_CPU_all == 0)) {
                bg = wxColour(255, 255, 200);                
            } else if(SI[i].nb_CPU_idle == 0) {
                bg = wxColour(255, 200, 200);                
            } else {
                bg = wxColour(200, 255, 200);                
            }
        } else {
            bg = wxColour(200, 200, 200);            
        }
        m_ClientTable->SetCellBackgroundColour(i, 0, bg);
        m_ClientTable->SetCellBackgroundColour(i, 1, bg);
        
   }   
}
void WXFoxMaster::UpdateLists(wxTimerEvent& event)
{
   if(m_grid_master==NULL) return;

   int port = m_grid_master->getServersPort();

   if(m_grid_master->isServerListening()) {
       if(m_ServerStatus!=NULL) {
            m_ServerStatus->SetForegroundColour(wxColour(0, 128, 0));
            m_ServerStatus->SetLabelText("Master status: Server is listening on port "+to_string(port));  

       }
   } else {
       if(m_ServerStatus!=NULL) {
            m_ServerStatus->SetForegroundColour(wxColour(128, 0, 0));
            m_ServerStatus->SetLabelText("Master status: ERROR - server is not listenning");            
       }
   }
   
   vector<FoxGridMaster::SLAVE_FOX_INFO> SI = m_grid_master->getFoxSlaveInfo();
   UpdateClientList(SI);

   UpdateJobList();
   UpdateResultList();
   
   m_UpdateTimer->Start(5*1000, true);
}

void WXFoxMaster::OnGridJobClick(wxGridEvent &event)
{
   int r = event.GetRow();
   int c = event.GetCol();

   if(m_JobListTable->GetSelectedRows().Count()==1){
      m_JobListTable->DeselectRow(m_JobListTable->GetSelectedRows().Item(0));
   }

   m_JobListTable->SelectRow(r, true);
}
void WXFoxMaster::OnGridResultClick(wxGridEvent &event)
{
   int r = event.GetRow();
   int c = event.GetCol();

   if(m_ResultTable->GetSelectedRows().Count()==1){
      m_ResultTable->DeselectRow(m_ResultTable->GetSelectedRows().Item(0));
   }

   m_ResultTable->SelectRow(r, true);
}
void WXFoxMaster::OnShowResultsServer(wxCommandEvent& event)
{
    int nb = m_ResultTable->GetSelectedRows().Count();
    if(nb!=1)
    {
        wxMessageBox(_T("Select one result!"), _T("Error"), wxOK, this);
        return;
    }    

    int r = m_ResultTable->GetSelectedRows().Item(0);
    long JobID;
    double Cost, Rwp;
    
    if(r<0 || r>=m_Results.size()) {
        wxMessageBox(_T("The selection os out of range!"), _T("Error"), wxOK, this);
        return;
    }   
    
    (*fpObjCrystInformUser)(wxString::Format("Show Results: opening file: "+m_Results[r].filename).ToStdString());

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

    if(m_Results[r].filename.size()>4 && m_Results[r].filename.Mid(m_Results[r].filename.size()-4)==wxString(_T(".xml"))) {    
        wxFileInputStream is(this->m_working_dir + "\\GridRslt\\" + m_Results[r].filename);
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
            wxMessageDialog d(this,_T("Failed loading file1:\n")+m_Results[r].filename,_T("Error loading file"),wxOK|wxICON_ERROR);
            d.ShowModal();
            return;
        };
        m_dataLoaded = true;
    }
}
void WXFoxMaster::OnShowResults(wxCommandEvent& event)
{    
    int nb = m_ResultTable->GetSelectedRows().Count();
    if(nb!=1)
    {
        wxMessageBox(_T("Select one result!"), _T("Error"), wxOK, this);
        return;
    }    

    int r = m_ResultTable->GetSelectedRows().Item(0);
    if(r<0 || r>=m_Results.size()) {
        wxMessageBox(_T("The selection os out of range!"), _T("Error"), wxOK, this);
        return;
    }                       

    wxString cmd;
    #ifdef WIN32
    cmd = "\"" + wxApp::GetInstance()->argv[0] + "\"";
    cmd +=_T(" ");
    cmd += "\"" + this->m_working_dir + "\\GridRslt\\" + m_Results[r].filename + "\"";
    wxExecute(cmd);
    (*fpObjCrystInformUser)(cmd.ToStdString());
    #else
    wxString appname = wxStandardPaths::Get().GetExecutablePath();
    wxString com=appname+_T(" ")+ this->m_working_dir + "/GridRslt/" +m_Results[r].filename;
   VFN_DEBUG_MESSAGE("WXFoxMaster::OnShowResults() command: "<<com, 10)
    long result= wxExecute(com);
    (*fpObjCrystInformUser)(com.ToStdString());    
    #endif
}
