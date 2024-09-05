#include "FoxGridMaster.h"

BEGIN_EVENT_TABLE(FoxGridMaster, wxEvtHandler)
    EVT_TIMER(ID_GRID_MASTER_CHECK_SLAVES,                         FoxGridMaster::OnCheckSlavesTimerEvent)
END_EVENT_TABLE()

FoxGridMaster::FoxGridMaster(wxString working_dir): GridMasterBase(working_dir)
{
    m_checkSlaveStateTimer = new wxTimer(this, ID_GRID_MASTER_CHECK_SLAVES);
    m_timer_iter = 1;
    m_checkSlaveStateTimer->Start(1000, true);
}
FoxGridMaster::~FoxGridMaster()
{

}
int FoxGridMaster::generateJobID()
{
    return generateUniqueID();
}
bool FoxGridMaster::addJob(wxString filename, wxString name, wxString data, int id, long nbOfTrial, long nbRun, bool rand)
{
    wxMutexLocker l(m_jobs_mutex);
    MasterJob mj;
    mj.filename = filename;
    mj.name = name;
    mj.data = data;
    mj.ID = id;
    mj.nb_trial = nbOfTrial;
    mj.nb_runs = nbRun;
    mj.randomize = rand;
    m_jobs.push_back(mj);
    WriteLogMessage("new job added");

    return true;
}
void FoxGridMaster::CheckIfJobsReceived()
{
    wxMutexLocker l(m_jobs_mutex);
    for(int i=0;i<m_jobs.size();i++) {
        for(int j=0;j<m_jobs[i].sentToClientsMsgInfo.size();j++) {
            if(isMsgReceived(m_jobs[i].sentToClientsMsgInfo[j].msgID)) {
                //m_jobs[i].nb_running += m_jobs[i].sentToClientsMsgInfo[j].nb_runs;
                m_jobs[i].runningIPs.insert(m_jobs[i].runningIPs.end(), m_jobs[i].sentToClientsMsgInfo[j].nb_runs, m_jobs[i].sentToClientsMsgInfo[j].SlaveIP);
                m_jobs[i].sentToClientsMsgInfo.erase(m_jobs[i].sentToClientsMsgInfo.begin()+j);
                j--;
            }
        }
    }

}
void FoxGridMaster::OnCheckSlavesTimerEvent(wxTimerEvent& event)
{
    WriteLogMessage("OnCheckSlavesTimerEvent");

    CheckIfJobsReceived();
    ProcessMessagesFromSlaves();

    if((m_timer_iter % 10)==0) {
        AskSlavesState();
    }

    m_checkSlaveStateTimer->Start(1000, true);
    m_timer_iter++;
}

void FoxGridMaster::SlaveDisconnected(wxString SlaveIP)
{
    wxMutexLocker l1(m_jobs_mutex);
    for(int j=0;j<m_jobs.size();j++) {
        for(int k=0;k<m_jobs[j].runningIPs.size();k++) {
            m_jobs[j].removeIPFromRunningIPs(SlaveIP, -1);
        }
    }
}
vector<FoxGridMaster::SLAVE_FOX_INFO> FoxGridMaster::getFoxSlaveInfo()
{
    vector<SLAVE_INFO_PUBLIC> pi = getSlavesPublicInfo();

    wxMutexLocker l(m_fox_slaves_info_mutex);
    //refresh the old m_fox_slaves_info
    for(int i=0;i<m_fox_slaves_info.size();i++) {
        bool found = false;
        int j=0;
        for(j=0;j<pi.size();j++) {
            if(pi[j].ip.compare(m_fox_slaves_info[i].sinfo.ip)==0) {
                found = true;
                break;
            }
        }
        if(found) {
            m_fox_slaves_info[i].sinfo = pi[j];
            m_fox_slaves_info[i].connected = true;
        } else if(m_fox_slaves_info[i].connected) {
            m_fox_slaves_info[i].connected = false;
            m_fox_slaves_info[i].nb_CPU_all = 0;
            m_fox_slaves_info[i].nb_CPU_idle = 0;
            m_fox_slaves_info[i].LasSentJobMsgID = 0;
            m_fox_slaves_info[i].LastProcessedJobMsgID = 0;

            SlaveDisconnected(m_fox_slaves_info[i].sinfo.ip);
        }
    }
    //finding new
    for(int j=0;j<pi.size();j++) {
        bool found = false;
        for(int i=0;i<m_fox_slaves_info.size();i++) {
            if(pi[j].ip.compare(m_fox_slaves_info[i].sinfo.ip)==0) {
                found = true;
                break;
            }
        }
        if(!found) {
            SLAVE_FOX_INFO sfi;
            sfi.connected = true;
            sfi.sinfo = pi[j];
            m_fox_slaves_info.push_back(sfi);
        }
    }

    return m_fox_slaves_info;
}
wxString FoxGridMaster::getResult(wxString message, long pos)
{
    wxString result;
    //copy of inmsg
    wxString in = message;
    in.Remove(0,pos);
    //find end of the job
    int p = in.First(_T("</result>"));
    if(p==-1) return _T("");
    //save it
    result = in.Left(p);
   //VFN_DEBUG_MESSAGE(__FUNCTION__<<":"<<message<<","<<in,10)
    return result;
}
bool FoxGridMaster::SaveResult(wxString filename, wxString data)
{
   #ifdef WIN32
      filename = m_working_dir + "\\GridRslt\\" + filename;
   #else
      filename = m_working_dir + "/GridRslt/" + filename;
   #endif
   WriteLogMessage(_T("Saving result as file " + filename));

   wxFile outFile(filename, wxFile::write);
   if(outFile.IsOpened()) {
      outFile.Write(data);
      outFile.Close();
      return true;
   }
   return false;
}
void FoxGridMaster::ProcessMessagesFromSlaves()
{
    vector<GridMasterBase::REC_MSG> msgs = getReceivedMsgs();

    for(int i=0;i<msgs.size();i++) {
        stringstream in_string;
        in_string<<msgs[i].msg;

        bool newResult = false;

        while(true) {
            XMLCrystTag tag;
            in_string>>tag;
            if(true==in_string.eof()) break;

            if("currentstate"==tag.GetName()) {
                int freeCPUs = 0;
                int AllCPUs = 0;
                long long LastProcessedJobMsgID = 0;
                wxString IP;

                for(int j=tag.GetNbAttribute()-1;j>=0;j--){
                    if(tag.GetAttributeName(j)=="freeCPUs") {
                        freeCPUs = atoi(tag.GetAttributeValue(j).c_str());
                    }
                    if(tag.GetAttributeName(j)=="IP") {
                        IP = wxString(tag.GetAttributeValue(j));
                    }
                    if(tag.GetAttributeName(j)=="AllCPUs") {
                        AllCPUs = atoi(tag.GetAttributeValue(j).c_str());
                    }
                    if(tag.GetAttributeName(j)=="LastProcessedJobMsgID") {
                        LastProcessedJobMsgID = stoll(tag.GetAttributeValue(j).c_str());
                    }
                }
                long long LasSentJobMsgID = 0;
                updateSlaveCPUsInfo(msgs[i].IP, freeCPUs, AllCPUs, LastProcessedJobMsgID, LasSentJobMsgID);

                if((freeCPUs!=0) && (LasSentJobMsgID==LastProcessedJobMsgID)) {
                    SendSomeJob(msgs[i].IP, freeCPUs);
                }
                break;
            }

            //<result ID=\"%d\" Cost=\"%s\" Rwp=\"%s\">
            if( ("result"==tag.GetName()) && (!tag.IsEndTag()) ) {
                newResult = true;
                long id=-1;
                double cost=-1, Rwp=-1;
                int count=0;
                for(int j=tag.GetNbAttribute()-1;j>=0;j--){
                    if(tag.GetAttributeName(j)=="ID") {
                        count++;
                        id = atol(tag.GetAttributeValue(j).c_str());
                    }
                    if(tag.GetAttributeName(j)=="Cost") {
                        count++;
                        cost = atof(tag.GetAttributeValue(j).c_str());
                    }
                    if(tag.GetAttributeName(j)=="Rwp") {
                        count++;
                        Rwp = atof(tag.GetAttributeValue(j).c_str());
                    }
                }
                long pos = in_string.tellg();
                if(cost<0) {
                    //error occured during the calculation and empty result was obtained from the slave, just delete from the job list one run
                    wxMutexLocker l(m_jobs_mutex);
                    for(int j=0;j<m_jobs.size();j++) {
                        if(m_jobs[j].ID == id) {
                            //just decrease nb_running...
                            m_jobs[j].removeIPFromRunningIPs(msgs[i].IP);
                            break;
                        }
                    }
                } else if(count!=3) {
                    WriteLogMessage("ERROR: Result received, but not all attributes read properly!");
                    //try to atleast delete it from the job list...
                    if(id!=-1) {
                        wxMutexLocker l(m_jobs_mutex);
                        for(int j=0;j<m_jobs.size();j++) {
                            if(m_jobs[j].ID == id) {
                                //just decrease nb_running...
                                m_jobs[j].removeIPFromRunningIPs(msgs[i].IP);
                                break;
                            }
                        }
                    }
                } else {
                    //everything seems to be OK -> save the result
                    wxString res = getResult(msgs[i].msg, pos);
                    MasterResult mr;
                    wxMutexLocker l(m_jobs_mutex);
                    for(int j=0;j<m_jobs.size();j++) {
                        if(m_jobs[j].ID == id) {
                            //just decrease nb_running...
                            m_jobs[j].removeIPFromRunningIPs(msgs[i].IP);
                            mr.Cost = cost;
                            mr.data = "";
                            mr.Rwp = Rwp;
                            mr.JobID = m_jobs[j].ID;
                            mr.filename = mr.generateFileName(m_jobs[j].ID, cost, Rwp);
                            m_jobs[j].results.push_back(mr);
                            break;
                        }
                    }
                    if(mr.Cost!=-1) {
                        if(SaveResult(mr.filename, res)) {
                            mr.data = "";
                        } else {
                            //TODO: save it later!
                            mr.data = "";
                        }
                    }
                }
                break;
            }
        }
    }
}
void FoxGridMaster::updateSlaveLastJobMsgsInfo(wxString SlaveIP, long long lastSentJobMsgID)
{
    wxMutexLocker l(m_fox_slaves_info_mutex);
    for(int i=0;i<m_fox_slaves_info.size();i++) {
        if(m_fox_slaves_info[i].sinfo.ip==SlaveIP) {
            m_fox_slaves_info[i].LasSentJobMsgID = lastSentJobMsgID;
        }
    }
}
void FoxGridMaster::updateSlaveCPUsInfo(wxString SlaveIP, int CPUfree, int CPUall, long long lastProcessedJobMsgID, long long &lastSentJobMsgID)
{
    lastSentJobMsgID = 0;
    wxMutexLocker l(m_fox_slaves_info_mutex);
    for(int i=0;i<m_fox_slaves_info.size();i++) {
        if(m_fox_slaves_info[i].sinfo.ip==SlaveIP) {
            m_fox_slaves_info[i].nb_CPU_all = CPUall;
            m_fox_slaves_info[i].nb_CPU_idle = CPUfree;
            m_fox_slaves_info[i].LastProcessedJobMsgID = lastProcessedJobMsgID;
            lastSentJobMsgID=m_fox_slaves_info[i].LasSentJobMsgID;
        }
    }
}
bool FoxGridMaster::existsJob(int ID)
{
    wxMutexLocker l(m_jobs_mutex);
    for(int i=0;i<m_jobs.size();i++) {
        if(m_jobs[i].ID == ID) return true;
    }
    return false;
}
vector<MasterJob> FoxGridMaster::getJobs()
{
    wxMutexLocker l(m_jobs_mutex);
    return m_jobs;
}
bool FoxGridMaster::DeleteJob(long JobID)
{
    wxMutexLocker l(m_jobs_mutex);
    for(int i=0;i<m_jobs.size();i++) {
        if(m_jobs[i].ID == JobID) {
            m_jobs.erase(m_jobs.begin()+i);
            return true;
        }
    }
    return false;
}
bool FoxGridMaster::UpdateJobHeader(MasterJob mj)
{
    bool found = false;
    wxMutexLocker l(m_jobs_mutex);
    for(int i=0;i<m_jobs.size();i++) {
        if(m_jobs[i].ID == mj.ID) {
            found = true;
            m_jobs[i].name = mj.name;
            m_jobs[i].nb_runs = mj.nb_runs;
            m_jobs[i].nb_trial = mj.nb_trial;
            m_jobs[i].randomize = mj.randomize;
        }
    }
    return found;
}
wxString FoxGridMaster::createJobHeader(wxString name, int ID, int nb_trial, int nb_runs, int randomize)
{
    wxString header;
    header.Printf(_T("<ClientJob Name=\"%s\" ID=\"%d\" nbOfTrial=\"%ld\" nbRun=\"%d\" rand=\"%d\">\n"), name, ID, nb_trial, nb_runs, randomize);
    return header;
}
bool FoxGridMaster::SendSomeJob(wxString SlaveIP, int freeCPUs)
{ //cannot use GetSlaveIndex here, because we have to be under mutex

    MasterJob mj;
    bool found = false;

    wxMutexLocker l(m_jobs_mutex);
    for(int i=0;i<m_jobs.size();i++) {
        if(m_jobs[i].getNbFreeJobs() > 0) {
            mj = m_jobs[i];
            found = true;

            int min = m_jobs[i].getNbFreeJobs() < freeCPUs ? m_jobs[i].getNbFreeJobs() : freeCPUs;
            wxString header = createJobHeader(mj.name, mj.ID, mj.nb_trial, min, mj.randomize);
            long long msgID = SendMsgToSlave(header + mj.data + "\n</ClientJob>\n", SlaveIP);
            if (msgID==-1) {
                return false;
            } else {
                m_jobs[i].sentToClientsMsgInfo.push_back(SentMsgInfo(msgID, min, SlaveIP));
                updateSlaveLastJobMsgsInfo(SlaveIP, msgID);
                freeCPUs = freeCPUs - min;
            }
            if(freeCPUs<=0) break;
        }
    }

    if(!found) {
        WriteLogMessage("No job found");
        return false;
    }
    return true;
}
void FoxGridMaster::AskSlavesState()
{
    SendMsgtoAll("state");
}
