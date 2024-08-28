#include "GridServer.h"

BEGIN_EVENT_TABLE(GridServer, wxEvtHandler)
   EVT_SOCKET(GRID_SERVER_ID,                       GridServer::OnServerEvent)
   EVT_TIMER(ID_SEND_TIMER,                         GridServer::OnTimerEvent)
END_EVENT_TABLE()

GridServer::GridServer(wxString working_dir)
{
    m_server = NULL;
    m_working_dir = working_dir;
    m_checkingTimer = new wxTimer(this, ID_SEND_TIMER);
    m_expect_port_sending_after_connection = false;
}
GridServer::~GridServer()
{
    delete(m_checkingTimer);
    m_checkingTimer = NULL;

    { 
        //releasing all these threads may take some time...
        wxMutexLocker l(m_server_threads_mutex);
        for(int i=0;i<m_server_threads.size();i++) {
            m_server_threads[i]->Delete();
        }
    }

    if(m_server!=NULL) {
        m_server->Destroy();
    }
}
void GridServer::SetWorkingDir(wxString path)
{
    m_working_dir = path;
}
void GridServer::setExpectPortSendingAfterConnection(bool exp)
{
    m_expect_port_sending_after_connection = exp;
}
bool GridServer::RunGridServer(unsigned short port)
{
   if(m_server!=NULL) return false;

   wxIPV4address ip;
   ip.Service(port);

   WriteLogMessage(_T("Starting server..."));
   m_server = new wxSocketServer(ip, wxSOCKET_REUSEADDR);
   m_server->SetEventHandler(*this, GRID_SERVER_ID);   
   m_server->SetNotify(wxSOCKET_CONNECTION_FLAG);
   m_server->Notify(true);
   m_server->SetTimeout(60);

   if (! m_server->Ok()) {
      WriteLogMessage(_T("ERROR: Server starting failure"));
      m_server->Destroy();
      m_server = NULL;
      return false;
   } else {
        WriteLogMessage(_T("Server started"));
        m_checkingTimer->Start(10*1000, true);
        return true;
   }
}
void GridServer::OnServerEvent(wxSocketEvent &event)
{
   WriteLogMessage(_T("On Server Event - start"));
   switch(event.GetSocketEvent())
   {
        case wxSOCKET_CONNECTION :
            {
                wxSocketBase* sock = m_server->Accept();
                HandleConnection(sock);
            }
            break;
        case wxSOCKET_INPUT :
            //This should never happend!
            WriteLogMessage(_T("WARNING: On Server Event: wxSOCKET_INPUT! Why?"));
            break;
        case wxSOCKET_LOST:
            WriteLogMessage(_T("wxSOCKET_LOST"));
            //HandleDisconnection(event.GetSocket());
            break;
        default:
            break;
   }
   WriteLogMessage(_T("On Server Event - end"));
}
int GridServer::readPortToConnectClient(wxSocketBase* sock)
{//reading 10 bytes with port number
    int port = -1;

    if(sock->WaitForRead(3, 0)) {       
        WriteLogMessage("Reading the port number to connect client");
        sock->SetFlags(wxSOCKET_WAITALL);
        char msg[10];
        sock->Read(msg, 10);
        if(sock->LastReadCount()!=10) {
            WriteLogMessage("ERROR: last Read counts: "+to_string(sock->LastReadCount()) +" != " + to_string(10));
            return -1;
        }
        if(sock->Error()==true) {
            WriteLogMessage("ERROR: last Read error: "+ to_string(sock->LastError()));
            return -1;
        }
        msg[9] = '\0';
        wxString tmpport(msg);
        tmpport.Trim();
        if(!tmpport.ToInt(&port)) return -1;
    }    
    return port;
}
void GridServer::HandleConnection(wxSocketBase* sock)
{
    wxIPV4address addr;
    
    if (!sock->GetPeer(addr))
    {
        WriteLogMessage("ERROR: Server cannot get peer info of a new connection!");
        sock->Close();
        return;
    } else {
        WriteLogMessage(wxString::Format("Got connection from %s:%d", addr.IPAddress().c_str(), addr.Service()));
    }
    sock->SetEventHandler(*this, GRID_SERVER_ID);
    sock->SetNotify(wxSOCKET_LOST_FLAG);
    sock->Notify(true);

    {
        wxMutexLocker l(m_server_threads_mutex);
        for(int i=0;i<m_server_threads.size();i++) {
            if(addr.IPAddress().compare(m_server_threads[0]->GetAddress().IPAddress())==0) {
                WriteLogMessage("ERROR: Server cannot accept this connection, there is something already connected to the server with the same IP!");
                sock->Close();
                return;
            }
        }
    }


    int port = -1;
    if(m_expect_port_sending_after_connection) {
        port = readPortToConnectClient(sock);
        if(port<0) {
            WriteLogMessage("ERROR: Server cannot accept this connection. Port number for connection of the client was not sent or error occured!");
            return;
        }
    }   

    {
        WriteLogMessage("Server: Creating Thread");
        SocketThreadServer *st = new SocketThreadServer(sock,
                                            m_working_dir, 
                                            this,
                                            addr,
                                            port);

        if (st->Create() == wxTHREAD_NO_ERROR) {
            WriteLogMessage("Server: Running Thread");
            st->Run();
            wxMutexLocker l(m_server_threads_mutex);
            m_server_threads.push_back(st);
        } else {
            //TODO release *st?
            WriteLogMessage("Server: cannot create next thread ");
            sock->Close();
            st->Delete();
        };
    }    
}
vector<SocketThreadInfo> GridServer::getSocketThreadsInfo()
{       
    vector<SocketThreadInfo> res;

    //refreshServerThreadList();

    wxMutexLocker l(m_server_threads_mutex);
    for(int i=0;i<m_server_threads.size();i++) {
        SocketThreadInfo si;
        //si.ID = m_server_threads[i]->GetId();
        si.address = m_server_threads[i]->GetAddress();
        si.port = m_server_threads[i]->getPortToConnectClient();
        res.push_back(si);
    }

    return res;
}
vector<GridCommunication::MSGINFO_REC> GridServer::getReceivedMsgs()
{
    vector<GridCommunication::MSGINFO_REC> msgs;
    wxMutexLocker l(m_messages_received_mutex);
    if (!l.IsOk()) {
        return msgs;
    }
    
    for(int i=0;i<m_messages_received.size();i++) {
        if(m_messages_received[i].processed==0) {
            msgs.push_back(m_messages_received[i]);
            m_messages_received[i].processed = getTimeStampMinutes();
        }
    }

    //Delete all processed messages
    m_messages_received.erase(std::remove_if(m_messages_received.begin(), m_messages_received.end(),
                              [](const MSGINFO_REC& msg) { return msg.processed != 0; }),
                              m_messages_received.end());

    return msgs;
}
bool GridServer::isConnected()
{//returns true if at least one connection (running thread) detected 

    if(m_server==NULL) return false;
    if(!m_server->IsConnected()) return false;
    
    wxMutexLocker l(m_server_threads_mutex);
    for(int i=0;i<m_server_threads.size();i++) {        
        if(m_server_threads[i]->IsRunning()) {
            return true;
        }
    }

    return false;
}
bool GridServer::isListening()
{
    if(m_server!=NULL) {
        return m_server->IsOk();
    }
    return false;
}
void GridServer::refreshServerThreadList()
{
    WriteLogMessage("refreshServerThreadList() - start");
    vector<MSGINFO_REC> msgs;
    {
        wxMutexLocker l(m_server_threads_mutex); 
        //if(m_server_threads.size()==0) return;
        for (auto it = m_server_threads.begin(); it != m_server_threads.end(); ) {
            if (!(*it)->IsRunning()) { 
                WriteLogMessage("Thread " + to_string((*it)->GetId()) + "is not running -> Save messages and delete it");
                //first save messages
                vector<MSGINFO_REC> ms;
                ms = (*it)->getReceivedMessages();
                msgs.insert(msgs.end(), ms.begin(), ms.end());
                //then delete the thread
                (*it)->Delete();
                it = m_server_threads.erase(it);
            } else {
                ++it;
            }
        }
    }   

    {//save it to received
        wxMutexLocker l(m_messages_received_mutex);
        m_messages_received.insert(m_messages_received.end(), msgs.begin(), msgs.end());
    }
     WriteLogMessage("refreshServerThreadList() - end");
}
void GridServer::getReceivedMsgsFromAllThreads()
{
    WriteLogMessage("getReceivedMsgsFromAllThreads() - start");
    vector<MSGINFO_REC> msgs;
    {//get all msgs from threads
        wxMutexLocker l(m_server_threads_mutex);
        vector<MSGINFO_REC> ms;
        for(int i=0;i<m_server_threads.size();i++) {
            ms = m_server_threads[i]->getReceivedMessages();
            msgs.insert(msgs.end(), ms.begin(), ms.end());
        }
    }
    WriteLogMessage("There are "+to_string(msgs.size())+" new messages");

    {//save it to received
        wxMutexLocker l(m_messages_received_mutex);
        m_messages_received.insert(m_messages_received.end(), msgs.begin(), msgs.end());
    }
    WriteLogMessage("getReceivedMsgsFromAllThreads() - end");
}
void GridServer::OnTimerEvent(wxTimerEvent& event)
{
    WriteLogMessage("OnTimerEvent - start");
    
    refreshServerThreadList();
    getReceivedMsgsFromAllThreads();

    m_checkingTimer->Start(5*1000, true);
    WriteLogMessage("OnTimerEvent - end");
}
wxString GridServer::GenerateRandomText(size_t length) {
    const char characters[] =
        "0123456789 ABCDEFGHI JKLM NOPQRSTUVWXYZ abcdefghijkl mnopqr stuvw xyz";
    std::random_device random_device;
    std::mt19937 generator(random_device());
    std::uniform_int_distribution<> distribution(0, sizeof(characters) - 2);

    wxString result;
    result.reserve(length);

    for(size_t i = 0; i < length; ++i) {
        result += characters[distribution(generator)];
    }

    return result;
}
