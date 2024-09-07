#include <wx/stdpaths.h>
#include "FoxGridSlave.h"

#define EVT_PROCESS_MY(func) DECLARE_EVENT_TABLE_ENTRY( wxEVT_PROCESS_MY, -1, -1, (wxObjectEventFunction) (wxEventFunction) (ProcessEventFunction) & func, (wxObject *) NULL ),
typedef void (wxEvtHandler::*ProcessEventFunction)(ProcessMyEvent&);

BEGIN_EVENT_TABLE(FoxGridSlave, wxEvtHandler)
    EVT_TIMER(ID_GRID_CHECK_SLAVE,                         FoxGridSlave::OnCheckSlaveTimerEvent)
    EVT_PROCESS_MY(FoxGridSlave::OnProcessEvent)
END_EVENT_TABLE()

MyProcess::MyProcess(FoxGridSlave *parent, const wxString& cmd, wxString dir)
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


FoxGridSlave::FoxGridSlave(wxString m_working_dir): GridSlaveBase(m_working_dir)
{
    wxFileName dirname = wxFileName::DirName(getWorkingDir());
    dirname.AppendDir("processes");
    wxString dir =  dirname.GetPath();
    if(!wxDirExists(dir)) wxMkdir(dir);

    m_checkSlaveTimer = new wxTimer(this, ID_GRID_CHECK_SLAVE);
    m_checkSlaveTimer->Start(3000, true);
    m_last_job_msg_processed = 0;
}
FoxGridSlave::~FoxGridSlave()
{
    wxMutexLocker locker(m_processes_mutex);
    //if there are running processes, kill them all
    for(int i=0;i<m_processes.size();i++) {
        if(m_processes[i].running) {
            wxKill(m_processes[i].pid, wxSIGKILL, NULL, wxKILL_CHILDREN);
        }
    }
}
void FoxGridSlave::resetProcesses(int nbProcesses)
{
    wxMutexLocker locker(m_processes_mutex);

    if(nbProcesses == m_processes.size()) return;

    //todo: make it smart
    // - if nbProcesses > m_processes.size(), then just add and do not kill all previous.
    // - else try to delete not running processes. if not possible kill only necessary number of them and not all of them

    //if there are running processes, kill them all
    for(int i=0;i<m_processes.size();i++) {
        if(m_processes[i].running) {
            wxKill(m_processes[i].pid, wxSIGKILL, NULL, wxKILL_CHILDREN);
        }
    }
    m_processes.clear();

    for(int i=0;i<nbProcesses;i++) {
        wxString dir;
        wxFileName dirname = wxFileName::DirName(getWorkingDir());
        dirname.AppendDir("processes");
        wxString tmpDir;
        tmpDir.Printf(_T("process_%d"), i);
        dirname.AppendDir(tmpDir);
        dir =  dirname.GetPath();

        if(!wxDirExists(dir.c_str())) {
            wxMkdir(dir);
        }
        FOX_PROCESS proc;
        proc.dir = dir;
        m_processes.push_back(proc);
    }
}
FOX_PROCESS *FoxGridSlave::getUnusedProcess()
{//this function itself can't lock m_ProcessMutex, because it returns pointer
 //the mutex has to be used by the caller

    for(int i=0;i<m_processes.size();i++) {
        if(!m_processes[i].running) {
            return &m_processes[i];
        }
    }
    return NULL;
}
void FoxGridSlave::setProcessUnused(int pid)
{
    wxMutexLocker locker(m_processes_mutex);
    for(int i=0;i<m_processes.size();i++) {
       if(m_processes[i].pid == pid) {
           m_processes[i].running = false;
           m_processes[i].pid = -1;
           break;
       }
   }
}
vector<FOX_PROCESS> FoxGridSlave::getProcesses()
{
    wxMutexLocker locker(m_processes_mutex);
    return m_processes;
}
bool FoxGridSlave::findResultedFile(wxString dir, wxString &Rwp, wxString &Cost, wxString &filename)
{
    wxDir tmpDir(dir);

    if(!tmpDir.IsOpened()) return false;
    //there should be only one file out-Cost-*.xml

    if(!tmpDir.GetFirst(&filename, _T("out-Cost*.xml"), wxDIR_FILES)) return false;

    //filename = tmpfile;

    int pos = filename.find("Cost");
    if(pos==-1) return false;
    int pos2 = filename.find("Rwp");
    if(pos2==-1) return false;

    Cost = filename.substr(pos+5,pos2-1-(pos+5));
    Cost.Replace(",", ".", true);

    pos = filename.find("Rwp");
    if(pos==-1) return false;
    pos2 = filename.find(".xml");
    if(pos2==-1) return false;

    Rwp = filename.substr(pos+4,pos2-(pos+4));
    Rwp.Replace(",", ".", true);

    return true;
}
void FoxGridSlave::OnProcessEvent(ProcessMyEvent& pEvent)
{
   wxString st = _T("");
   wxString filename, Cost, Rwp;
   st.Printf(_T("pid=%d, status=%d, dir="), pEvent.m_pid, pEvent.m_status);
   WriteLogMessage(_T("Process terminated: ") + st + pEvent.m_dir);

   int jobID=-1;
   wxDateTime startingtime = wxDateTime();
   //identify process ID
   {
       wxMutexLocker locker(m_processes_mutex);
       for(int i=0;i<m_processes.size();i++) {
           if(m_processes[i].pid == pEvent.m_pid) {
               jobID = m_processes[i].jobID;
               startingtime = m_processes[i].startingtime;
               break;
           }
       }
   }

   if(jobID==-1) {

       WriteLogMessage(_T("ERROR: Process JobID not found!"));
       setProcessUnused(pEvent.m_pid);

   } else if(pEvent.m_status!=0) {

       WriteLogMessage(_T("ERROR: process terminated with an error => empty result will be sent to the server"));
       SaveResult("", "-1", "-1", jobID, startingtime, true);

   } else if(!findResultedFile(pEvent.m_dir, Rwp, Cost, filename)) {

       WriteLogMessage("ERROR: parsing filename error: filename=" + filename + ", cost=" + Cost + ", Rwp=" + Rwp);
       SaveResult("", "-1", "-1", jobID, startingtime, true);

   } else {

       #ifdef WIN32
       SaveResult(pEvent.m_dir + _T("\\") + filename, Cost, Rwp, jobID, startingtime, false);
       #else
       SaveResult(pEvent.m_dir + _T("/") + filename, Cost, Rwp, jobID, startingtime, false);
       #endif
   }

    #ifdef WIN32
    if(filename.length()!=0) {
        if(wxFileExists(pEvent.m_dir + _T("\\") + filename)) wxRemoveFile(pEvent.m_dir + _T("\\") + filename);
    }
    if(wxFileExists(pEvent.m_dir + _T("\\input.xml"))) wxRemoveFile(pEvent.m_dir + _T("\\input.xml"));
    #else
    if(filename.length()!=0) {
        wxRemoveFile(pEvent.m_dir + _T("/") + filename);
    }
    wxRemoveFile(pEvent.m_dir + _T("/input.xml"));
    wxRemoveFile(pEvent.m_dir + _T("/out.txt"));
    #endif

    setProcessUnused(pEvent.m_pid);

    return;

}
bool FoxGridSlave::LoadFile(wxString filename, wxString &in)
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
bool FoxGridSlave::SaveResult(wxString fileName, wxString Cost, wxString Rwp, int ID, wxDateTime startingtime, bool error)
{
    WriteLogMessage(_T("Saving new result"));
    wxString in;
    wxString out;

    if(!error) {
        if(!LoadFile(fileName, in)) {
            WriteLogMessage("ERROR: can't load the file with result: " + fileName);
            return false;
        }

        out.Printf(_T("<result ID=\"%d\" Cost=\"%s\" Rwp=\"%s\">\n"), ID, Cost.c_str(), Rwp.c_str());
        out += in;
        out += _T("</result>\n");
    } else {
        wxString out;
        out.Printf(_T("<result ID=\"%d\" Cost=\"-1\" Rwp=\"-1\" >\n"), ID);
        out += in;
        out += _T("</result>\n");
    }

    MasterResult mr;
    mr.Cost = atof(Cost);
    mr.Rwp = atof(Rwp);
    mr.data = out;

    bool found = false;
    wxMutexLocker l(m_job_mutex);
    for(int i=0;i<m_jobs.size();i++) {
        if(m_jobs[i].ID == ID) {
            m_jobs[i].results.push_back(mr);
            m_jobs[i].removeIPFromRunningIPs("local", 1);
            wxLongLong Avtime = m_jobs[i].average_calc_time.GetSeconds();
            wxTimeSpan newt = wxDateTime::Now() - startingtime;
            double p = (Avtime.ToDouble()*m_jobs[i].results.size() + newt.GetSeconds().ToDouble()) / ((double) m_jobs[i].results.size()+1.0);
            m_jobs[i].average_calc_time = wxTimeSpan::Seconds(p);
            found = true;
            break;
        }
    }
    return found;
}
/*
void FoxGridSlave::SendResult(wxString result, long duration, long jobID, long resultID)
{
    wxString outmsg = "duration="+to_string(duration)+" jobID="+to_string(jobID)+" resultID="+to_string(resultID)+"\n"+result;
    WriteLogMessage("GridSlave::SendResult");

    if(!sendMessage(outmsg)) {
        WriteLogMessage("ERROR: sendMessage returns false!");

        //TODO process result later
        return;
    }

    ResultInfo ri;
    ri.data = " ";
    ri.duration_seconds = duration;
    ri.jobID = jobID;
    ri.resultID = resultID;
    wxMutexLocker l(m_results_mutex);
    m_results.push_back(ri);
}
*/
bool FoxGridSlave::SaveDataAsFile(wxString out, wxString filename)
{
   wxFile outFile(filename, wxFile::write);
   if(outFile.IsOpened()) {
      outFile.Write(out);
      outFile.Close();
      return true;
   }
   return false;
}
bool FoxGridSlave::DeleteXMLFilesInDirectory(const wxString& directoryPath) {
    wxDir dir(directoryPath);
    if (!dir.IsOpened()) {
        // Directory couldn't be opened
        return false;
    }

    wxString filename;
    bool hasFiles = dir.GetFirst(&filename, "*.xml", wxDIR_FILES);
    while (hasFiles) {
        wxFileName fullPath(directoryPath, filename);
        if (!wxRemoveFile(fullPath.GetFullPath())) {
            // Failed to remove a file, might want to handle this case
        }
        hasFiles = dir.GetNext(&filename);
    }

    return true;
}
void FoxGridSlave::CheckResultsAndSendOne()
{
    wxMutexLocker l(m_job_mutex);
    for(int i=0;i<m_jobs.size();i++) {
        for(int j=0;j<m_jobs[i].results.size();j++) {
            if(m_jobs[i].results[j].data.length()!=0) {
                if(sendMessage(m_jobs[i].results[j].data)>0) {
                    m_jobs[i].results[j].data = "";
                }
                return;
            }
        }
    }
}
void FoxGridSlave::CheckJobsAndStartCalculation()
{
    wxMutexLocker l(m_job_mutex);
    for(int i=0;i<m_jobs.size();i++) {
        while(m_jobs[i].getNbFreeJobs()>0) {
            wxMutexLocker ll(m_processes_mutex);
            FOX_PROCESS *proc = getUnusedProcess();
            if(proc==NULL) return;

            #ifdef WIN32
            //saving input file to the temporary directory
            DeleteXMLFilesInDirectory(proc->dir);
            if(!SaveDataAsFile(m_jobs[i].data,proc->dir + _T("\\input.xml"))) {
                WriteLogMessage("ERROR: Can't write a file " + proc->dir + _T("\\input.xml"));
                return;
            }
            wxString cmd_content;
            wxString appname = wxApp::GetInstance()->argv[0];
            if(rand) {
			    cmd_content.Printf(_T("%3$s \"%1$s\\input.xml\" --nogui -n %2$d --randomize -o \"%1$s\\out#cost#Rwp.xml\""), proc->dir, m_jobs[i].nb_trial, appname);
            } else {
                cmd_content.Printf(_T("%3$s \"%1$s\\input.xml\" --nogui -n %2$d -o \"%1$s\\out#cost#Rwp.xml\""), proc->dir, m_jobs[i].nb_trial, appname);
            }
            #else
            //WriteLogMessage(_T("saving data as file - input.xml"));
            DeleteXMLFilesInDirectory(proc->dir);
            if(!SaveDataAsFile(m_jobs[i].data,proc->dir + _T("/input.xml"))) {
               WriteLogMessage("ERROR: Can't write a file " + proc->dir + _T("/input.xml"));
               return;
            }
            wxFile *batFile = new wxFile();
            batFile->Create(proc->dir + _T("/run.sh"), true, wxS_IRUSR|wxS_IWUSR|wxS_IXUSR|wxS_IRGRP|wxS_IWGRP|wxS_IXGRP|wxS_IROTH|wxS_IWOTH|wxS_IXOTH);
            if(!batFile->IsOpened()) {
                batFile->Open(proc->dir + _T("/run.sh"),  wxFile::write);
            }
            wxString content;
            wxString appname = wxStandardPaths::Get().GetExecutablePath();
            WriteLogMessage(appname);
            //if(appname(0,1)!=_T("/")) appname=wxGetCwd()+_T("/")+appname;
            //WriteLogMessage(appname);
            wxString tr;
            tr.Printf(_T("%d"), m_jobs[i].nb_trial);
            if(rand) {
                content = appname+_T(" ")+proc->dir+_T("/input.xml --nogui -n ")+tr+_T(" --randomize --silent -o ")+proc->dir+_T("/out#cost#Rwp.xml > ")+proc->dir+_T("/out.txt");
            } else {
                content = appname+_T(" ")+proc->dir+_T("/input.xml --nogui -n ")+tr+_T(" --silent -o ")            +proc->dir+_T("/out#cost#Rwp.xml > ")+proc->dir+_T("/out.txt");
            }
            //wxString path_to_input;
            //path_to_input = proc->getTmpDir();
            //WriteLogMessage(path_to_input);
    /*
            if(rand) {
                content = appname+_T(" ")+path_to_input+_T("/input.xml --nogui -n ")+tr+_T(" --randomize -o ")+path_to_input+_T("/out#cost.xml");
            } else {
                content = _T(".") + appname+_T(" ")+path_to_input+_T("/input.xml --nogui -n ")+tr+_T(" -o ")+path_to_input+_T("/out#cost.xml");
            }
    */
            WriteLogMessage(content);
            cout<<content.ToAscii()<<endl;
	         batFile->Write(content);
            batFile->Close();
            delete batFile;
            #endif


            #ifdef WIN32
            //wxString cmd = proc->getTmpDir() + _T("\\!run.bat");
            #else
            wxString cmd_content = _T("sh ")+ proc->dir + _T("/run.sh");
            #endif
            //WriteLogMessage(_T("creating process"));
            wxProcess *process = new MyProcess(this, cmd_content, proc->dir);
            //WriteLogMessage(_T("executing process"));
            int pid = wxExecute(cmd_content, wxEXEC_ASYNC|wxEXEC_MAKE_GROUP_LEADER, process);
            if ( !pid ) {
                delete process;
                WriteLogMessage(_T("Process not started: ") + cmd_content);
                //m_ProcessMutex->Unlock();
                return;
            } else {
                m_jobs[i].runningIPs.push_back("local");
                proc->pid = pid;
                proc->running = true;
                proc->jobID = m_jobs[i].ID;
                proc->startingtime =  wxDateTime::Now();
                wxString tmp;
                tmp.Printf(_T("%d"), pid);
                WriteLogMessage(_T("Process started: pid=") + tmp + _T(", dir=") + cmd_content);
                //m_ProcessMutex->Unlock();
                return;
            }
        }
    }
}
void FoxGridSlave::OnCheckSlaveTimerEvent(wxTimerEvent& event)
{
    WriteLogMessage("GridSlave::OnCheckSlaveTimerEvent");
    vector<GridCommunication::MSGINFO_REC> msgs = getReceivedMessages();
    WriteLogMessage("GridSlave: There are new "+to_string(msgs.size())+" received messages");

    if(msgs.size()!=0) {
        ProcessMsgs(msgs);
    }

    CheckJobsAndStartCalculation();
    CheckResultsAndSendOne();

    m_checkSlaveTimer->Start(1000, true);
}
void FoxGridSlave::addJob(MasterJob job)
{
    wxMutexLocker l(m_job_mutex);
    for(int i=0;i<m_jobs.size();i++) {
        if(m_jobs[i].ID == job.ID) {
            m_jobs[i].nb_runs += job.nb_runs;
            return;
        }
    }
    m_jobs.push_back(job);
}
wxString FoxGridSlave::getMyHostname()
{
    wxIPV4address addr;
    addr.Hostname(wxGetFullHostName());
    return addr.IPAddress();
}

wxString FoxGridSlave::getMyStateMsg()
{
    wxString out;
    out = _T("<FoxGrid>\n <currentstate ");
    out << _T("IP=\"") << getMyHostname()<< _T("\" ");
    out << _T("freeCPUs=\"") << getCPUsFree() << _T("\" ");
    out << _T("AllCPUs=\"") << getCPUsAll() << _T("\" ");
    out << _T("LastProcessedJobMsgID=\"") << to_string(m_last_job_msg_processed) << _T("\" ");
    out << _T(" />\n");
    out << _T("</FoxGrid>\n");
    return out;
}
wxString FoxGridSlave::getJobData(wxString inmsg, long pos)
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

void FoxGridSlave::ProcessMsgs(vector<GridCommunication::MSGINFO_REC> &msgs)
{
    //just for testing...
    for(int i=0;i<msgs.size();i++) {
        if(msgs[i].msg.compare("state")==0) {
            WriteLogMessage("ProcessMsgs: asking for my state, sending answer");
            sendMessage(getMyStateMsg());
        } else {
            WriteLogMessage("Job received");
            //WriteLogMessage(msgs[i].msg);
            stringstream in_string;
            in_string<<msgs[i].msg;
            int id=0, trial=0, runs=0, rand=0;
            wxString name;

            while(true) {
                XMLCrystTag tag;
                in_string>>tag;
                if(true==in_string.eof()) break;
                if( ("ClientJob"==tag.GetName()) && (!tag.IsEndTag()) ) {

                    //WriteLogMessage(_T("New job found"));
                    for(int j=tag.GetNbAttribute()-1;j>=0;j--) {
                        if(tag.GetAttributeName(j)=="Name") {
                            //WriteLogMessage(_T("ID found"));
                            name = wxString::FromAscii(tag.GetAttributeValue(j).c_str());
                        }
                        if(tag.GetAttributeName(j)=="ID") {
                            //WriteLogMessage(_T("ID found"));
                            wxString ID = wxString::FromAscii(tag.GetAttributeValue(j).c_str());
                            if(!ID.ToInt(&id)) WriteLogMessage(_T("Can't convert ID attribute to int"));
                        }
                        if(tag.GetAttributeName(j)=="nbOfTrial"){
                            //WriteLogMessage(_T("nbTrials found"));
                            wxString nbTrial = wxString::FromAscii(tag.GetAttributeValue(j).c_str());
                            if(!nbTrial.ToInt(&trial)) WriteLogMessage(_T("Can't convert nbTrials attribute to int"));
                        }
                        if(tag.GetAttributeName(j)=="nbRun") {
                            //WriteLogMessage(_T("nbOfRuns found"));
                            wxString nbRuns = wxString::FromAscii(tag.GetAttributeValue(j).c_str());
                            if(!nbRuns.ToInt(&runs)) WriteLogMessage(_T("Can't convert nbOfRuns attribute to int"));
                        }
                        if(tag.GetAttributeName(j)=="rand") {
                            //WriteLogMessage(_T("rand found"));
                            wxString Rand = wxString::FromAscii(tag.GetAttributeValue(j).c_str());
                            if(!Rand.ToInt(&rand)) WriteLogMessage(_T("Can't convert rand attribute to int"));
                        }
                    }

                    long pos = in_string.tellg();
                    wxString tmp = getJobData(msgs[i].msg, pos);
                    tmp.Trim();
                    if(tmp.length()!=0) {

                        MasterJob mj;
                        mj.data = tmp;
                        mj.ID = id;
                        mj.name = name;
                        mj.nb_runs = runs;
                        mj.nb_trial = trial;
                        mj.randomize = rand;
                        mj.msgID=msgs[i].ID;
                        addJob(mj);

                        WriteLogMessage("New Job Saved");
                    } else {
                        WriteLogMessage(_T("ERROR: job was not load"));
                    }
                    m_last_job_msg_processed = msgs[i].ID;
                    break;
                }
            }
        }
    }
}
vector<MasterJob> FoxGridSlave::getJobs()
{
    wxMutexLocker l(m_job_mutex);
    return m_jobs;
}
void FoxGridSlave::ResetNbCPUsAll(int CPUs)
{
    resetProcesses(CPUs);
}
int FoxGridSlave::getCPUsAll()
{
    wxMutexLocker l(m_processes_mutex);
    return m_processes.size();
}
int FoxGridSlave::getCPUsFree()
{
    int nb = 0;
    {
        wxMutexLocker l(m_processes_mutex);
        for(int i=0;i<m_processes.size();i++) {
            if(!m_processes[i].running) {
                nb++;
            }
        }
    }
    {
        wxMutexLocker l(m_job_mutex);
        for(int i=0;i<m_jobs.size();i++) {
            nb -= m_jobs[i].getNbFreeJobs();
        }
    }
    if (nb<0) nb = 0;

    return nb;
}
