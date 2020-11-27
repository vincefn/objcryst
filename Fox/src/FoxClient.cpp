#ifdef __WX__CRYST__
   #include "ObjCryst/wxCryst/wxCrystal.h"
#endif

using namespace ObjCryst;
using namespace std;

#include "FoxClient.h"
#include "wx/stdpaths.h"

static const long GRID_CLIENT_SOCKET_ID=                  WXCRYST_ID();
static const long ID_UPDATE_TIMER_CLIENT=                   WXCRYST_ID();
static const long ID_SEND_TIMER=                         WXCRYST_ID();

#define EVT_PROCESS_MY(func) DECLARE_EVENT_TABLE_ENTRY( wxEVT_PROCESS_MY, -1, -1, (wxObjectEventFunction) (wxEventFunction) (ProcessEventFunction) & func, (wxObject *) NULL ),

typedef void (wxEvtHandler::*ProcessEventFunction)(ProcessMyEvent&);



BEGIN_EVENT_TABLE(FoxClient, wxEvtHandler)
   EVT_SOCKET(GRID_CLIENT_SOCKET_ID,                FoxClient::OnSocketEvent)
   //EVT_TIMER(ID_SEND_TIMER,                  FoxClient::OnSendResults)
    EVT_TIMER(ID_SEND_TIMER,                  FoxClient::OnTimerEvent)
    EVT_PROCESS_MY(FoxClient::OnProcessEvent)
    //EVT_UPDATE_UI(ID_CRYST_UPDATEUI,                FoxClient::OnUpdateUI)
END_EVENT_TABLE()

///////////////////////////////////////////////
MyProcess::MyProcess(FoxClient *parent, const wxString& cmd, wxString dir)
    : m_cmd(cmd)
{
    m_parent = parent;
   m_dir = dir;
}
void MyProcess::OnTerminate(int pid, int status)
{
    ProcessMyEvent *e = new ProcessMyEvent(this);
    e->m_exit = true;
    e->m_pid = pid;
    e->m_dir = m_dir;
    e->m_status = status;
    //m_parent->AddPendingEvent( e );
    m_parent->QueueEvent(e);
    //m_parent->onProcessTerminate(pid, status, m_dir);
    delete this;
}
///////////////////////////////////////////////
FoxProcess::FoxProcess(wxString tmpDir)
{
    tmpDIR = tmpDir;
    running = false;
    pid = -1;
    jobID = -1;
    startingtime = wxDateTime::Now();//just to initiate it
}
FoxProcess::~FoxProcess()
{
}
void FoxProcess::setPid(int pid)
{
    this->pid = pid;
}
int FoxProcess::getPid()
{
    return pid;
}
void FoxProcess::setTmpDir(wxString dir)
{
    tmpDIR = dir;
}
wxString FoxProcess::getTmpDir()
{
    return tmpDIR;
}
void FoxProcess::setRunning(bool run)
{
    running = run;
}
bool FoxProcess::isRunning()
{
    return running;
}
void FoxProcess::setJobID(int id)
{
    jobID = id;
}
int FoxProcess::getJobID()
{
    return jobID;
}
void FoxProcess::setStarted(wxDateTime t)
{
    startingtime = t;
}
wxDateTime FoxProcess::getStartingTime()
{
    return startingtime;
}
///////////////////////////////////////////////
GrdRslt::GrdRslt(int ID, wxString cost, wxString content)
{
    this->ID = ID;
    this->cost = cost;
    this->content = content;
    sent = false;
	pending = false;
}
GrdRslt::~GrdRslt()
{
}
///////////////////////////////////////////////


///////////////////////////////////////////////
FoxClient::FoxClient(wxString working_dir):
wxEvtHandler()
{
   m_working_dir = working_dir;
   wxString dirName = addToPath(getWorkingDir(), _T("processes"));
   if(!wxDirExists(dirName)) wxMkdir(dirName);
   //m_DataMutex = new wxMutex();
   //m_ResultsMutex = new wxMutex();
   //m_ProcessMutex = new wxMutex();
   mpClient = new wxSocketClient();
   //m_Connecting = false;
   m_sendingTimer = new wxTimer(this, ID_SEND_TIMER);
   m_exit = false;
   m_nbOfAvailCPUs = wxThread::GetCPUCount();
   resetProcesses(m_nbOfAvailCPUs);
}
FoxClient::~FoxClient()
{
   if(m_sendingTimer!=NULL) {
      m_sendingTimer->Stop();
      delete m_sendingTimer;
   }

   if (mpClient != 0) {
      mpClient->Destroy();
	  mpClient = 0;
   }
}
void FoxClient::setNbOfAvailCPUs(int nb)
{
    m_nbOfAvailCPUs = nb;
    resetProcesses(m_nbOfAvailCPUs);
}
int FoxClient::getNbOfAvailCPUs()
{
    return m_nbOfAvailCPUs;
}
void FoxClient::resetProcesses(int nbProcesses)
{
    wxCriticalSectionLocker locker(m_ProcessCriticalSection);
    
    m_processes.clear();
    for(int i=0;i<nbProcesses;i++) {
      wxString dir;
      dir = addToPath(getWorkingDir(), _T("processes"));
      wxString tmpDir;
      tmpDir.Printf(_T("process_%d"), i);
      dir = addToPath(dir, tmpDir);
      if(!wxDirExists(dir.c_str())) {
         wxMkdir(dir);
      }
/*
      #ifdef __WIN32__
      dir.Printf(_T("\\ObjCryst-FoxGrid-process_%d"), i);
      #else
      dir.Printf(_T("/ObjCryst-FoxGrid-process_%d"), i);
      #endif
      dir=wxStandardPaths::Get().GetTempDir()+dir;
      if(!wxDirExists(dir.c_str())) {
         wxMkdir(dir);
      }
*/
      FoxProcess proc(dir);
      m_processes.push_back(proc);
   }
}
void FoxClient::WriteMessageLog(wxString msg)
{
#if __CLIENT_LOGS
   wxString filename = addToPath(this->getWorkingDir(), _T("client.log"));
   wxFile logfile(filename, wxFile::write_append);
   wxDateTime datetime = wxDateTime::Now();
   logfile.Write(datetime.Format(_T("%X ")) + msg + _T("\n"));
   logfile.Close();
#endif
   (*fpObjCrystInformUser)(msg.ToStdString());
}
void FoxClient::OnProcessEvent(ProcessMyEvent& pEvent)
{
    onProcessTerminate(pEvent.m_pid, pEvent.m_status, pEvent.m_dir);
}
void FoxClient::onProcessTerminate(int pid, int status, wxString dir)
{
   int nb=0;
   //WriteMessageLog(_T("Going to lock"));
   wxCriticalSectionLocker locker(m_ProcessCriticalSection);
   //WriteMessageLog(_T("Going to locked"));

   wxString st = _T("");
   st.Printf(_T("pid=%d, status=%d, dir="), pid, status);
   WriteMessageLog(_T("Process terminated: ") + st + dir);
   int ID;
   //identify process ID
   
   for(int i=0;i<m_processes.size();i++) {
       if(m_processes[i].getPid() == pid) {
           m_processes[i].setRunning(false);
           m_processes[i].setPid(-1);
           ID = m_processes[i].getJobID();
           break;
       }
   }
   
   
   //WriteMessageLog(_T("onProcessTerminate: Mutex Unlocked"));
   //if some error return
   if(status!=0) {
       WriteMessageLog(_T("process terminated with an error => empty result will be sent to the server"));
       SaveResult("", "-1", ID, true);
       return;
   }

   wxString filename = getOutputFile(dir);
   wxString cost = getCost(filename);
   #ifdef WIN32
   //save result to the memory
   SaveResult(dir + _T("\\") + filename, cost, ID, false);
   //delete output and input files
   wxRemoveFile(dir + _T("\\") + filename);
   wxRemoveFile(dir + _T("\\input.xml"));

   #else
   //save result to the memory
   SaveResult(dir + _T("/") + filename, cost, ID, false);
   //delete output and input files
   wxRemoveFile(dir + _T("/") + filename);
   wxRemoveFile(dir + _T("/input.xml"));
   wxRemoveFile(dir + _T("/out.txt"));
   #endif
   //WriteMessageLog(_T("Leaving the locked part..."));

   //wxMutexError qre = mmh.Unlock();
   //WriteMessageLog(_T("Mutex unlocked"+wxString::Format("%d", qre)));
}
wxString FoxClient::getOutputFile(wxString dir)
{
   wxDir tmpDir(dir);
   wxString filename;

   if(!tmpDir.IsOpened()) return _T("-1");
   //there should be only one file out-Cost-*.xml
   if(!tmpDir.GetFirst(&filename, _T("out-Cost*.xml"), wxDIR_FILES)) return _T("-1");
   return filename;
}
wxString FoxClient::getCost(wxString filename)
{
    if(filename.size()<=13) return _T("-1");
    return filename.SubString(9, filename.size()-5);
}
void FoxClient::WriteProtocol()
{

}
bool FoxClient::IsClientConnected()
{
   if(mpClient!=0)   return mpClient->IsConnected();
   else return false;
}
/*
void FoxClient::Reconnect()
{
    
    //stop timer
    if(m_sendingTimer!=NULL) {
       m_sendingTimer->Stop();
    }
    WriteMessageLog(_T("Reconnecting..."));
    Disconnect();

    //clear results
    //m_results.clear();

    //wait a few seconds...
    wxSleep(10);
    ConnectClient(3, m_hostname);

    //start timer
    if(m_sendingTimer!=NULL) {
        m_sendingTimer->Start(30*1000, false);
    }
    
}
*/
void FoxClient::KillProcesses()
{
    wxCriticalSectionLocker locker(m_ProcessCriticalSection);
    WriteMessageLog(_T("Killing processes"));
    for(int i=0;i<m_processes.size();i++) {
        if(m_processes[i].isRunning()) {
            wxKillError er;
            wxKill(m_processes[i].getPid(), wxSIGKILL, &er, wxKILL_CHILDREN);
            wxString tmp;
            tmp.Printf(_T("Process terminated with error code: %d"), er);
            WriteMessageLog(tmp);
        }
    }
}
void FoxClient::Disconnect()
{
    WriteMessageLog(_T("Disconnecting"));
    //m_sendingTimer->Stop();
    if(mpClient!=0) {
       mpClient->Destroy();
       mpClient = 0;
    }
}
bool FoxClient::ConnectClient(int nbOfTrial, wxString hostname)
{
   //m_Connecting = true;
   WriteMessageLog(_T("Client try to connect"));
   wxIPV4address ip;
   int i=0;

   m_hostname = hostname;
   ip.Service(2854);
   if(!ip.Hostname(m_hostname)) {
       return false;
   }

   if(mpClient==0) mpClient = new wxSocketClient();
   mpClient->SetEventHandler(*this, GRID_CLIENT_SOCKET_ID);
   mpClient->SetNotify( wxSOCKET_CONNECTION_FLAG |
                        wxSOCKET_INPUT_FLAG |
                        wxSOCKET_LOST_FLAG | wxSOCKET_OUTPUT_FLAG);
   mpClient->Notify(true);

   do
   {
      i++;
      mpClient->Connect(ip, false);
      mpClient->WaitOnConnect(3,0);
      if(i==nbOfTrial) break;
      //if(!m_Connecting) break;
      if(m_exit) break;
   } while(!mpClient->IsConnected());
   
   //m_Connecting = false;
   
   if (mpClient->IsConnected()){
      mpClient->SaveState();
      //better to run it just once and then start it again after all precedures are done
      m_sendingTimer->Start(10*1000, true);
      WriteMessageLog(_T("Client is connected to the server"));
      return true;
   }
   else {
     mpClient->Close();
     return false;
   }
}

void FoxClient::OnSocketEvent(wxSocketEvent &event)
{
   //WriteMessageLog(_T("OnSocketEvent Begin"));
   wxSocketBase *tmpSock = event.GetSocket();

   tmpSock->SetNotify(wxSOCKET_LOST_FLAG);

   switch(event.GetSocketEvent())
   {
     case wxSOCKET_INPUT:
        {
           //WriteMessageLog(_T("INPUT - This should never happen!"));

           //WriteMessageLog(_T("INPUT_END"));
           break;
        }
     case wxSOCKET_LOST:
        {
          //ClientEvent("Connection failed.\n");
           //WriteMessageLog(_T("LOST - This can sometimes happen!"));
           //m_sendingTimer->Stop();
           if(mpClient!=0) {
               mpClient->Destroy();
               mpClient = NULL;
           }
           break;
        }
     case wxSOCKET_OUTPUT:
         {
             //WriteMessageLog(_T("OUTPUT - This should never happen!"));
             break;
         }
     default:
       break;
   }
   tmpSock->SetNotify(wxSOCKET_LOST_FLAG | wxSOCKET_INPUT_FLAG | wxSOCKET_OUTPUT_FLAG);
   //WriteMessageLog(_T("OnSocketEvent End"));
}
 
bool FoxClient::AnalyzeMessage(wxString msg)
{
   wxString ID, nbTrial, nbRuns, Rand;
   stringstream in_string;

   bool newJob = false;
   vector<wxString> asks;
   vector<wxString> jobs;
   vector<long> jobRuns;
   vector<long> trials;
   vector<long> ids;
   vector<bool> rands;

   SaveDataAsFile(msg, addToPath(getWorkingDir(), _T("client_msg_in.txt")));
   //jobs = getJobs(inmsg);
   in_string<<msg;
   //WriteMessageLog(_T("Start parsing file"));
   while(true) {
          XMLCrystTag tag;
          in_string>>tag;
          if(true==in_string.eof()) break;
          if( ("ClientJob"==tag.GetName()) && (!tag.IsEndTag()) ){

             //WriteMessageLog(_T("New job found"));
             newJob = true;
             long runs=0, trial=0, id=0, rand;
             for(int i=tag.GetNbAttribute()-1;i>=0;i--){
                if(tag.GetAttributeName(i)=="ID"){
                   //WriteMessageLog(_T("ID found"));
                   ID = wxString::FromAscii(tag.GetAttributeValue(i).c_str());
                   if(!ID.ToLong((long *) &id)) WriteMessageLog(_T("Can't convert ID attribute to long"));
                }
                if(tag.GetAttributeName(i)=="nbTrials"){
                   //WriteMessageLog(_T("nbTrials found"));
                   nbTrial = wxString::FromAscii(tag.GetAttributeValue(i).c_str());
                   if(!nbTrial.ToLong((long *) &trial)) WriteMessageLog(_T("Can't convert nbTrials attribute to long"));
                }
                if(tag.GetAttributeName(i)=="nbOfRuns") {
                    //WriteMessageLog(_T("nbOfRuns found"));
                    wxString nbRuns = wxString::FromAscii(tag.GetAttributeValue(i).c_str());
                    if(!nbRuns.ToLong((long *) &runs)) WriteMessageLog(_T("Can't convert nbOfRuns attribute to long"));
                }
                if(tag.GetAttributeName(i)=="rand") {
                    //WriteMessageLog(_T("rand found"));
                    wxString Rand = wxString::FromAscii(tag.GetAttributeValue(i).c_str());
                    if(!Rand.ToLong((long *) &rand)) WriteMessageLog(_T("Can't convert rand attribute to long"));
                }

             }
             long pos = in_string.tellg();
             wxString tmp = getJob(msg, pos);
             if(tmp.Cmp(_T(""))!=0) {
                 jobs.push_back(tmp);
                 jobRuns.push_back(runs);
                 trials.push_back(trial);
                 ids.push_back(id);
                 rands.push_back(rand);
             } else {
                 WriteMessageLog(_T("ERROR: job was not load"));
             }
          }
   }

   if(newJob){
	WriteMessageLog(_T("New job was received from server"));
        std::vector<int> jobsForRejecting;
        for(int i=0;i<jobs.size();i++) {
            for(int run=0;run<jobRuns[i];run++) {
                //if job does not run, reject it...
                if(runNewJob(jobs[i], ids[i], (int) trials[i], rands[i])!=0) {
                    jobsForRejecting.push_back(ids[i]);
                }
            }
        }
        //TODO
        //rejectJobs(jobsForRejecting);
   }
   return true;
}
wxString FoxClient::getMyHostname()
{
    wxIPV4address addr;
    addr.Hostname(wxGetFullHostName());
    return addr.IPAddress();
}
void FoxClient::SendCurrentState()
{
    //WriteMessageLog(_T("Sending a current state info to server"));

    wxString out;
    //starting tag
    out = _T("<FoxGrid>\n <currentstate ");
    out << _T("name=\"") << getMyHostname() << _T("\" ");
    out << _T("freeCPUs=\"") << getNbOfUnusedProcesses() << _T("\" ");
    out << _T("availableCPUs=\"") << getNbOfAvailCPUs() << _T("\" ");
    out += _T(" />\n");

    bool newResults = false;
    //send all result with value sent=false
    for(int i=0;i<m_results.size();i++) {
        if(m_results[i].sent==false && m_results[i].pending==false) {
            //WriteMessageLog(_T("Result found - will be sent as well"));
            out << m_results[i].content+_T("\n");
            m_results[i].pending=true;
            newResults = true;
        }
    }
    out << _T("</FoxGrid>\n");

    if(newResults) {
        WriteMessageLog(_T("Sending result(s) to server"));
    }

    //WriteMessageLog(_T("Saving answer to the file"));
    //SaveDataAsFile(out, addToPath(getWorkingDir(), _T("client_out.txt")));

    //WriteMessageLog(_T("sending message..."));
    if(!m_IOSocket.WriteStringToSocket(mpClient, out)) {
        WriteMessageLog(m_IOSocket.getError());
        WriteMessageLog(_T("ERROR: answer was not send"));
        //change back status of results and send it later...
        for(int i=0;i<m_results.size();i++) {
            if(m_results[i].sent==false && m_results[i].pending==true) {
                m_results[i].pending=false;
            }
        }
        //Reconnect();
    } else {
        //if sent, change the status of results
        for(int i=0;i<m_results.size();i++) {
            if(m_results[i].sent==false && m_results[i].pending==true) {
                m_results[i].sent=true;
                m_results[i].pending==false;
            }
        }
    }

    //WriteMessageLog(_T("Message sent"));
}
wxString FoxClient::getJob(wxString inmsg, long pos)
{
    wxString job;
    //copy of inmsg
    wxString in = inmsg;
    in.Remove(0,pos);
    //find end of the job
    int p = in.First(_T("</ClientJob>"));
    if(p==-1) return _T("");
    //save it
    job = in.Left(p);
    return job;
}
vector<FoxProcess> FoxClient::get_copy_of_processes()
{
    WriteMessageLog("get_copy_of_processes: Locking");
    vector<FoxProcess> res;
    wxCriticalSectionLocker locker(m_ProcessCriticalSection);  
    WriteMessageLog("get_copy_of_processes: Locked");
    for(int i=0;i<m_processes.size();i++) {
        res.push_back(m_processes[i]);
    }
    WriteMessageLog("get_copy_of_processes: Locking");
    return res;
}
int FoxClient::getNbOfUnusedProcesses()
{
    WriteMessageLog("getNbOfUnusedProcesses: Locking");
    wxCriticalSectionLocker locker(m_ProcessCriticalSection);
    WriteMessageLog("getNbOfUnusedProcesses: Locked");
    int nb=0;
    for(int i=0;i<m_processes.size();i++) {
        if (!m_processes[i].isRunning()) nb++;
    }
    WriteMessageLog("getNbOfUnusedProcesses: Unlocked");
    return nb;
}
FoxProcess * FoxClient::getUnusedProcess()
{//this function itself can't lock m_ProcessMutex, because it returns pointer
 //the mutex has to be used by the caller
    
    for(int i=0;i<m_processes.size();i++) {
        if(!m_processes[i].isRunning()) {
            return &m_processes[i];
        }
    }
    return 0;
}
int FoxClient::runNewJob(wxString job, int id, int nbTrial, bool rand)
{
    WriteMessageLog(_T("running new job"));
    //we have to use protection for m_processes here!
    //if(m_ProcessMutex->Lock()!=wxMUTEX_NO_ERROR) return -1;
    wxCriticalSectionLocker locker(m_ProcessCriticalSection);

    FoxProcess *proc = getUnusedProcess();
    if(proc!=0) {
        //WriteMessageLog(_T("process found"));
        #ifdef WIN32
        //saving input file to the temporary directory
        SaveDataAsFile(job,proc->getTmpDir() + _T("\\input.xml"));
        //creating bat file
        //wxFile *batFile = new wxFile(proc->getTmpDir() + _T("\\!run.bat"), wxFile::write);
        //creating bat file content
        wxString cmd_content;
		wxString appname = wxApp::GetInstance()->argv[0];
        if(rand) {
			cmd_content.Printf(_T("%3$s %1$s\\input.xml --nogui -n %2$d --randomize -o %1$s\\out#cost.xml"), proc->getTmpDir(), nbTrial, appname);
        } else {
            cmd_content.Printf(_T("%3$s %1$s\\input.xml --nogui -n %2$d -o %1$s\\out#cost.xml"), proc->getTmpDir(), nbTrial, appname);
        }
        #else
        //WriteMessageLog(_T("saving data as file - input.xml"));
        SaveDataAsFile(job, proc->getTmpDir() + _T("/input.xml"));
        wxFile *batFile = new wxFile();
        batFile->Create(proc->getTmpDir() + _T("/run.sh"), true, wxS_IRUSR|wxS_IWUSR|wxS_IXUSR|wxS_IRGRP|wxS_IWGRP|wxS_IXGRP|wxS_IROTH|wxS_IWOTH|wxS_IXOTH);
        if(!batFile->IsOpened()) {
            batFile->Open(proc->getTmpDir() + _T("/run.sh"),  wxFile::write);
        }
        wxString content;
        wxString appname = wxStandardPaths::Get().GetExecutablePath();
        WriteMessageLog(appname);
        //if(appname(0,1)!=_T("/")) appname=wxGetCwd()+_T("/")+appname;
        //WriteMessageLog(appname);
        wxString tr;
        tr.Printf(_T("%d"), nbTrial);
        if(rand) {
            content = appname+_T(" ")+proc->getTmpDir()+_T("/input.xml --nogui -n ")+tr+_T(" --randomize --silent -o ")+proc->getTmpDir()+_T("/out#cost.xml > ")+proc->getTmpDir()+_T("/out.txt");
        } else {
            content = appname+_T(" ")+proc->getTmpDir()+_T("/input.xml --nogui -n ")+tr+_T(" --silent -o ")            +proc->getTmpDir()+_T("/out#cost.xml > ")+proc->getTmpDir()+_T("/out.txt");
        }
        //wxString path_to_input;
        //path_to_input = proc->getTmpDir();
        //WriteMessageLog(path_to_input);
/*
        if(rand) {
            content = appname+_T(" ")+path_to_input+_T("/input.xml --nogui -n ")+tr+_T(" --randomize -o ")+path_to_input+_T("/out#cost.xml");
        } else {
            content = _T(".") + appname+_T(" ")+path_to_input+_T("/input.xml --nogui -n ")+tr+_T(" -o ")+path_to_input+_T("/out#cost.xml");
        }
*/
        WriteMessageLog(content);
        //cout<<content.ToAscii()<<endl;
	    batFile->Write(content);
        batFile->Close();
        delete batFile;
        #endif


        #ifdef WIN32
        //wxString cmd = proc->getTmpDir() + _T("\\!run.bat");
        #else
        wxString cmd_content = _T("sh ")+ proc->getTmpDir() + _T("/run.sh");
        #endif
        //WriteMessageLog(_T("creating process"));
        wxProcess *process = new MyProcess(this, cmd_content, proc->getTmpDir());
        //WriteMessageLog(_T("executing process"));
        int pid = wxExecute(cmd_content, wxEXEC_ASYNC|wxEXEC_MAKE_GROUP_LEADER, process);
        if ( !pid ) {
            delete process;
            WriteMessageLog(_T("Process not started: ") + cmd_content);
            //m_ProcessMutex->Unlock();
            return -1;
        } else {
            proc->setPid(pid);
            proc->setRunning(true);
            proc->setJobID(id);
            proc->setStarted(wxDateTime::Now());
            wxString tmp;
            tmp.Printf(_T("%d"), pid);
            WriteMessageLog(_T("Process started: pid=") + tmp + _T(", dir=") + cmd_content);
            //m_ProcessMutex->Unlock();
            return 0;
        }
    } else {
        WriteMessageLog(_T("No unsused process is available, job will be rejected to server"));
        //m_ProcessMutex->Unlock();
        return -1;
    }
    //m_ProcessMutex->Unlock();
    return 0;
}
/*
void FoxClient::rejectJobs(std::vector<int> ids)
{
    
    if(ids.size()==0) return;

    WriteMessageLog(_T("Rejecting jobs..."));
    wxString out;
    out.Printf(_T("<FoxGrid>\n"), ids.size());

    for(int i=0;i<ids.size();i++) {
        wxString jobTag;
        jobTag.Printf(_T("  <rejectedJob id=\"%d\" />\n"), ids[i]);
        out+=jobTag;
    }
    out+=_T("</FoxGrid>\n");
    WriteMessageLog(_T("Saving output message to the file"));
    SaveDataAsFile(out, addToPath(getWorkingDir(), _T("client_out.txt")));

    WriteMessageLog(_T("sending message..."));
    wxSleep(2);
    if(!m_IOSocket.WriteStringToSocket(mpClient, string(out.ToAscii()))) {
        WriteMessageLog(m_IOSocket.getError());
        Reconnect();
    }
    WriteMessageLog(_T("Rejecting jobs...end"));
}
*/
void FoxClient::SaveDataAsFile(wxString out, wxString filename)
{
   wxFile outFile(filename, wxFile::write);
   if(outFile.IsOpened())
   {
      outFile.Write(out);
      outFile.Close();
   }
}
wxString FoxClient::getWorkingDir()
{
    return m_working_dir;
}
wxString FoxClient::addToPath(wxString str1, wxString str2)
{
    wxString res;
    #ifdef WIN32
       res = str1 + _T("\\") + str2;
    #else
       res = str1 + _T("/") + str2;
    #endif
    return res;
}
bool FoxClient::LoadFile(wxString filename, wxString &in)
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
void FoxClient::SaveResult(wxString fileName, wxString Cost, int ID, bool error)
{
    WriteMessageLog(_T("Saving new result"));
    wxString in;
    if(!error) {
        if(!LoadFile(fileName, in)) {
            WriteMessageLog(_T("can't load the file"));
            return;
        }
        wxString out;
        out.Printf(_T("<result ID=\"%d\" Cost=\"%s\">\n"), ID, Cost.c_str());
        out += in;
        out += _T("</result>\n");
        
        GrdRslt res(ID, Cost, out);
        m_results.push_back(res);
    } else {
        wxString out;
        out.Printf(_T("<result ID=\"%d\" Cost=\"-1\">\n"), ID);
        out += in;
        out += _T("</result>\n");

        GrdRslt res(ID, Cost, out);
        m_results.push_back(res);
    }
    
    //WriteMessageLog(_T("Result saved"));
}
void FoxClient::OnTimerEvent(wxTimerEvent& event)
{
    DoManyThingsOnTimer();
    //timer should be called  only for one shot. 
    //Imagine, that the communication with server is not done and another timer event occur...
    m_sendingTimer->Start(10*1000, true);
}
void FoxClient::DoManyThingsOnTimer()
{
    //stop timer until this function finishes - then we have to run it again
    
    if(mpClient == NULL) {  
        
        return;
    }
    if (!mpClient->IsConnected()) {
        //TODO: try to connect...

        return;
    }
    //1. send info to server (regularly)
    SendCurrentState();

    //2. wait for answer
    wxString msg;
    if(!m_IOSocket.ReadStringFromSocket(mpClient, msg)) {
        return;
    }

    //3. analyze answer, todo
    if(!AnalyzeMessage(msg)) {
        return;
    }
}
