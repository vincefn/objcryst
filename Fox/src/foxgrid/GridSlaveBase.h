#ifndef __GRIDSLAVE_H // Make sure to only declare these classes once
#define __GRIDSLAVE_H

#include "GridServer.h"
#include "GridClient.h"
#include <vector>

using namespace std;

class GridSlaveBase
{
    public:        
        GridSlaveBase(wxString working_dir);
        ~GridSlaveBase();
        bool ConnectToMaster(int nbOfTrial, wxString hostname, int port=2854);                
        bool isConnectedToMaster();
        bool isOutcommingConnectionRunning();
        bool isReceivingConnectionRunning();
        wxString getHostName();
        wxString getWorkingDir();
        long long sendMessage(wxString msg, long long msgID=-1);
        bool isMsgDelivered(long long msgID);
        vector<GridCommunication::MSGINFO_REC> getReceivedMessages();

    protected:
        void WriteLogMessage(wxString msg); 
        
        wxMutex   m_log_msg_mutex;

    private:
        wxString     m_working_dir;
        wxString     m_hostname;
        int          m_port_host;
        int          m_connection_attempts;
        GridClient  *m_client;//for outcomming msgs
        GridServer  *m_server;//for incomming msgs
        int          m_port_server;
        
        
};
#endif

