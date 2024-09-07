#include "GridSlaveBase.h"



GridSlaveBase::GridSlaveBase(wxString working_dir)
{
    m_working_dir = working_dir;
    m_client = NULL;
    m_server = NULL;
    m_port_host = 2854;
    m_port_server = 2853;
}
GridSlaveBase::~GridSlaveBase()
{
    delete(m_client);
    delete(m_server);
}
wxString GridSlaveBase::getWorkingDir()
{
    return m_working_dir;
}
wxString GridSlaveBase::getHostName()
{
    return m_hostname;
}
bool GridSlaveBase::ConnectToMaster(int nbOfTrial, wxString hostname, int port)
{
    m_port_host = port;
    m_hostname = hostname;
    m_connection_attempts = nbOfTrial;
    if(m_server==NULL) {
        m_port_server = 2853;
        //search for the free port
        vector<int> used_ports = GridCommunication::getUsedPorts();

        while (find(used_ports.begin(), used_ports.end(), m_port_server) != used_ports.end()) {
            m_port_server++;
        }

        m_server = new GridServer(m_working_dir);
        m_server->setExpectPortSendingAfterConnection(false);
        if(!m_server->RunGridServer(m_port_server)) {
            delete(m_server);
            m_server = NULL;
            return false;
        }
    }
    if(m_client==NULL) {
        m_client = new GridClient(m_working_dir, "");
        m_client->setAutoReconnectWhenConnectionLost(true);

        if(!m_client->ConnectClient(nbOfTrial, hostname, m_port_host, m_port_server)) {
            delete(m_client);
            m_client = NULL;
            return false;
        }
        //todo
        //sendMessage("newConnection port="+to_string(m_port_server));
    }
    return true;
}
long long GridSlaveBase::sendMessage(wxString msg, long long msgID)
{
    if(msgID==-1) {
        msgID = GridCommunication::getTimeStampNanoSeconds();
    }
    if(m_client!=0) {
        m_client->addMessageToSend(msg, msgID);
        return msgID;
    }
    return -1;
}
bool GridSlaveBase::isMsgDelivered(long long msgID)
{
    if(m_client!=0) {
        vector<GridCommunication::MSGINFO_SENT> msgs = m_client->getCopyMsgsToBeSent();
        for (auto it = msgs.rbegin(); it != msgs.rend(); ++it) {
            if (it->ID == msgID) {
                if(it->delivery_confirmation_obtained>0) return true;
                else return false;
            }
        }
    }
    //if not foun returns false;
    return false;
}
vector<GridCommunication::MSGINFO_REC> GridSlaveBase::getReceivedMessages()
{
    if(m_server!=NULL) {
         return m_server->getReceivedMsgs();
    }
    vector<GridCommunication::MSGINFO_REC> dummy;
    return dummy;
}
bool GridSlaveBase::isOutcommingConnectionRunning()
{
    if(m_client==0) return false;
    return m_client->isConnected();
}
bool GridSlaveBase::isReceivingConnectionRunning()
{
    if(m_server==NULL) return false;
    return m_server->isConnected();
}
bool GridSlaveBase::isConnectedToMaster()
{
    if(isOutcommingConnectionRunning() && isReceivingConnectionRunning()) {
        return true;
    }
    return false;
}
void GridSlaveBase::WriteLogMessage(wxString msg)
{
#ifdef ALLOW_GRID_LOGS
   wxString name;
#ifdef WIN32
   name = m_working_dir + _T("\\") + "Slave_log.txt";
#else
   name = m_working_dir + _T("/") + "Slave_log.txt";
#endif

   wxMutexLocker l(m_log_msg_mutex);

   wxFile logfile(name, wxFile::write_append);
   if(logfile.IsOpened())
   {
      wxDateTime datetime = wxDateTime::Now();
      logfile.Write(datetime.Format(_T("%X ")) + msg + _T("\n"));
      logfile.Close();
   }
#endif
}
