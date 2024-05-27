#include "GridMasterBase.h"



GridMasterBase::GridMasterBase(wxString working_dir)
{
    m_working_dir = working_dir;
    m_server = NULL;
    m_port = 2854;
}
GridMasterBase::~GridMasterBase()
{
    if(m_server!=NULL) {
        delete(m_server);
    }
    wxMutexLocker l(m_slaves_mutex);
    for(int i=0;i<m_slaves.size();i++) {
        delete(m_slaves[i].grid_client);
    }
}
int GridMasterBase::GetSlaveIndex(wxString IP)
{  
    wxMutexLocker l(m_slaves_mutex);
    for(int i=0;i<m_slaves.size();i++) {
        if(m_slaves[i].ip.compare(IP)==0) return i;
    }
    return -1;
}
int GridMasterBase::getServersPort()
{
    return m_port;
}
bool GridMasterBase::InitializeCommunication() //connections for sending
{
    //search for the free port
    vector<int> used_ports = GridCommunication::getUsedPorts();
    
    while (find(used_ports.begin(), used_ports.end(), m_port) != used_ports.end()) {
        m_port++;
    }

    WriteLogMessage("Starting the Master server for incomming msgs");
    if(m_server==NULL) {
        m_server = new GridServer(m_working_dir);
        m_server->setExpectPortSendingAfterConnection(true);

        if(!m_server->RunGridServer(m_port)) {
            delete(m_server);
            m_server = NULL;
            return false;
        }
    }
    
    return true;
}
bool GridMasterBase::isServerListening()
{
    if(m_server!=NULL) {
        return m_server->isListening();
    }
    return false;
}
long long GridMasterBase::SendMsgtoAll(wxString msg, long long msgID)
{   
    refreshSlaveList();
    if(msgID==-1) {
        msgID = GridCommunication::getTimeStampNanoSeconds();
    }
    wxMutexLocker l(m_slaves_mutex);
    for(int i=0;i<m_slaves.size();i++) {
        m_slaves[i].grid_client->addMessageToSend(msg, msgID);
    }
    return msgID;
}
long long GridMasterBase::SendMsgToSlave(wxString msg, wxString slaveIP, long long msgID, bool refresh_slave_list)
{    
    if(refresh_slave_list) {
        refreshSlaveList();
    }

    wxMutexLocker l(m_slaves_mutex);

    int slaveindex=-1;
    for(int i=0;i<m_slaves.size();i++) {
        if(m_slaves[i].ip.compare(slaveIP)==0) slaveindex = i;
    }

    if(slaveindex<0) {
        WriteLogMessage("ERROR: Unknown slave (ID="+slaveIP+"). Job will not be sent...");
        return -1;
    }
                            
    if(msgID==-1) {
        msgID = GridCommunication::getTimeStampNanoSeconds();
    }
    if(m_slaves[slaveindex].grid_client!=NULL) {
        m_slaves[slaveindex].grid_client->addMessageToSend(msg, msgID);
    } else {
        WriteLogMessage("ERROR: grid_client=0! (ID="+slaveIP+"). Job will not be sent...");
        return -1;
    }
    WriteLogMessage("Job sent");
    
    return msgID;
}
vector<GridMasterBase::REC_MSG> GridMasterBase::getReceivedMsgs()
{
    vector<GridMasterBase::REC_MSG> res;

    if(m_server==NULL) {
        return res;
    }

    vector<GridCommunication::MSGINFO_REC> msgs = m_server->getReceivedMsgs();
    for(int i=0;i<msgs.size();i++) {
        REC_MSG rm;
        rm.msgID = msgs[i].ID;
        rm.IP = msgs[i].IP;
        //rm.slaveID = msgs[i].socketID;
        rm.msg = msgs[i].msg;
        res.push_back(rm);
    }
    return res;
}
vector<wxString> GridMasterBase::getSlavesIPs()
{
    refreshSlaveList();

    vector<wxString> res;
    wxMutexLocker l(m_slaves_mutex);
    for(int i=0;i<m_slaves.size();i++) {
        res.push_back(m_slaves[i].ip);
    }
    return res;
}
vector<GridMasterBase::SLAVE_INFO_PUBLIC> GridMasterBase::getSlavesPublicInfo()
{
    refreshSlaveList();

    vector<SLAVE_INFO_PUBLIC> res;
    wxMutexLocker l(m_slaves_mutex);
    for(int i=0;i<m_slaves.size();i++) {
        SLAVE_INFO_PUBLIC pi;
        //pi.ID = m_slaves[i].ID;
        pi.ip = m_slaves[i].ip;
        pi.port = m_slaves[i].port;
        res.push_back(pi);
    }
    return res;
}
void GridMasterBase::refreshSlaveList()
{
    //to avoid refreshing every time
    long long now = GridCommunication::getTimeStampNanoSeconds();    
    if((now - m_slaves_refresh_time) < 3) return;

    if(m_server!=NULL) {        
        vector<SocketThreadInfo> stis = m_server->getSocketThreadsInfo();
        //checking new connections
        for(int i=0;i<stis.size();i++) {
            if(GetSlaveIndex(stis[i].address.IPAddress())==-1) {
                WriteLogMessage("New server thread detected (ip="+stis[i].address.IPAddress()+")");
                //there is a new connection!
                SLAVE_INFO_BASE si;

                si.grid_client = new GridClient(m_working_dir, stis[i].address.IPAddress());
                si.grid_client->setAutoReconnectWhenConnectionLost(true);

                if(!si.grid_client->ConnectClient(1, stis[i].address.IPAddress(), stis[i].port)) {
                    WriteLogMessage("ERROR: Starting outcomming port failed");
                    si.grid_client->Delete();
                    si.grid_client = NULL;
                } else {
                    //si.ID = stis[i].ID;
                    si.ip = stis[i].address.IPAddress();
                    si.port = stis[i].port;
                    wxMutexLocker l(m_slaves_mutex);
                    m_slaves.push_back(si);
                    WriteLogMessage("New outcomming port with (IP="+stis[i].address.IPAddress()+" established");
                }
            }            
        }
        
        //checking disconnections
        wxMutexLocker l(m_slaves_mutex);
        for (auto it_slaves = m_slaves.begin(); it_slaves != m_slaves.end(); ) {
            bool found = false;
            for (auto it_stis = stis.begin(); it_stis != stis.end(); ++it_stis) {
                if(it_stis->address.IPAddress().compare(it_slaves->ip)==0) {
                    found = true;
                    break;
                }
            }
            if(!found) {
                //just delete                 
                //Clients usually try to reconnect, but here we are deleting it, so, just let him know and it will delete itself..
                it_slaves->grid_client->Delete();                                
                it_slaves = m_slaves.erase(it_slaves);
            } else {
                ++it_slaves;
            }
        }        
    }
}
bool GridMasterBase::isMsgReceived(long long msgID)
{
    refreshSlaveList();

    wxMutexLocker l(m_slaves_mutex);
    GridCommunication::MSGINFO_SENT msg;
    for (int i=0;i<m_slaves.size();i++) {
        if(m_slaves[i].grid_client!=NULL) {
            if(m_slaves[i].grid_client->getCopyMsgToBeSent(msgID, msg)) {
                if(msg.delivery_confirmation_obtained>0) {
                    return true;
                } else {
                    return false;
                }
            }
        }        
    }        
    return false;
}
void GridMasterBase::WriteLogMessage(wxString msg)
{
#ifdef ALLOW_GRID_LOGS
   wxString name;
#ifdef WIN32
   name = m_working_dir + _T("\\") + "Master_log.txt";
#else
   name = m_working_dir + _T("/") + "Master_log.txt";
#endif
   wxFile logfile(name, wxFile::write_append);
   if(logfile.IsOpened())
   {
      wxDateTime datetime = wxDateTime::Now();
      logfile.Write(datetime.Format(_T("%X ")) + msg + _T("\n"));
      logfile.Close();
   } 
#endif
}

