#ifndef __GRIDMASTER_H // Make sure to only declare these classes once
#define __GRIDMASTER_H

#include "GridServer.h"
#include "GridClient.h"
#include <vector>
#include "wx/event.h"

using namespace std;

class GridMasterBase
{
    public:

        struct REC_MSG {
            //IP of the sender
            wxString    IP;
            long long msgID;
            //long slaveID;
            wxString msg;

            REC_MSG() {
                msgID = 0;
                //slaveID = 0;
            }
        };

         struct SLAVE_INFO_PUBLIC {
                //long          ID;
                wxString      ip;
                int port;

                SLAVE_INFO_PUBLIC() {
                    port = -1;
                }
            };


        GridMasterBase(wxString working_dir);
        ~GridMasterBase();

        bool isServerListening();
        bool InitializeCommunication();
        int getServersPort();
        //call think about calling refreshSlaveList();
        long long SendMsgtoAll(wxString msg, long long msgID=-1);
        long long SendMsgToSlave(wxString msg, wxString slaveIP, long long msgID=-1, bool refresh_slave_list=false);
        vector<GridMasterBase::REC_MSG> getReceivedMsgs();
        vector<wxString> getSlavesIPs();
        vector<SLAVE_INFO_PUBLIC> getSlavesPublicInfo();
        bool isMsgReceived(long long msgID);

    protected:
        void refreshSlaveList();
        void WriteLogMessage(wxString msg);
        int GetSlaveIndex(wxString IP);

        wxString                 m_working_dir;
        int                      m_port;

    private:
        struct SLAVE_INFO_BASE {
                //long          ID;
                wxString      ip;
                int           port;
                GridClient    *grid_client;

                SLAVE_INFO_BASE() {
                    //ID = 0;
                    grid_client = 0;
                    port = -1;
                }
            };

        vector<SLAVE_INFO_BASE>  m_slaves; //for outcomming msgs
        long                     m_slaves_refresh_time;
        wxMutex                  m_slaves_mutex;

        GridServer              *m_server; //for incomming msgs


};
#endif
