#ifndef __GRIDSERVER_H // Make sure to only declare these classes once
#define __GRIDSERVER_H

#include <wx/frame.h>
#include <wx/socket.h>
#include <wx/file.h>
#include "GridCommunication.h"
#include "wx/timer.h"
#include "SocketThreadServer.h"
#include <random>



class GridServer: public wxEvtHandler, GridCommunication
{
    public:
        GridServer(wxString working_dir);
        ~GridServer();

        bool RunGridServer(unsigned short port);
        void OnServerEvent(wxSocketEvent &event);
        void HandleConnection(wxSocketBase* sock);
        void OnTimerEvent(wxTimerEvent& event);
        vector<SocketThreadInfo> getSocketThreadsInfo();
        void SetWorkingDir(wxString path);
        bool isConnected();
        bool isListening();

        //void OnWorkerEvent(WorkerEvent& pEvent);

        vector<GridCommunication::MSGINFO_REC> getReceivedMsgs();

        void setExpectPortSendingAfterConnection(bool exp);


    protected:
        void getReceivedMsgsFromAllThreads();
        void refreshServerThreadList();
        int readPortToConnectClient(wxSocketBase* sock);

        wxString GenerateRandomText(size_t length);
        wxSocketServer            * m_server;
        vector<SocketThreadServer*> m_server_threads;
        wxMutex                     m_server_threads_mutex;
        //bool                        m_isRunning;
        vector<MSGINFO_REC>         m_messages_received;
        wxMutex                     m_messages_received_mutex;
        wxTimer                   * m_checkingTimer;
        bool                        m_expect_port_sending_after_connection;


        DECLARE_EVENT_TABLE()

};
#endif
