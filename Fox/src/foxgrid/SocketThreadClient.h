#ifndef __SOCKETTHREADCLIENT_H // Make sure to only declare these classes once
#define __SOCKETTHREADCLIENT_H

#include <wx/thread.h>
#include "GridCommunication.h"

class SocketThreadClient : public wxThread, GridCommunication
{
    public:
        SocketThreadClient(wxSocketBase         *pSocket,
                     wxString                    workingDir,
                     wxEvtHandler               *parent,
                     wxIPV4address               address,
                     vector<MSGINFO_SENT>       *messages_out,
                     wxMutex                    *messages_out_mutex);

        virtual ExitCode Entry();

        long GetId();
        wxIPV4address GetAddress();


    private:
        void WriteLogMessage(wxString msg);
        bool isMessageToBeSent();
        bool SndMsg();

        vector<MSGINFO_SENT> *m_messages_out;
        wxMutex              *m_messages_out_mutex;

        wxSocketBase        *m_socket;
        wxEvtHandler        *m_parent;
        wxIPV4address        m_address;
        bool                 m_socket_error;
};
#endif
