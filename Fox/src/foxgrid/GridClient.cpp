#include "GridClient.h"

BEGIN_EVENT_TABLE(GridClient, wxEvtHandler)
   EVT_SOCKET(GRID_CLIENT_SOCKET_ID,                GridClient::OnSocketEvent)
    EVT_TIMER(ID_SEND_TIMER,                  GridClient::OnTimerEvent)
END_EVENT_TABLE()


GridClient::GridClient(wxString working_dir, wxString IP)
{
    m_working_dir = working_dir;
    m_sendingTimer = new wxTimer(this, ID_SEND_TIMER);
    m_reconnect_when_connection_lost = true;
    m_need_socket_connect = false;
    m_socket_thread_client = NULL;
    m_IP = IP;
    m_keep_at_least_last_m_messages_to_be_send = 50;
    m_keep_minutes_m_messages_to_be_send = 10;
    m_delete = false; 
}
GridClient::~GridClient()
{
    delete(m_sendingTimer);
    m_sendingTimer = NULL;

    {
        wxMutexLocker l(m_socket_thread_client_mutex);
        if(m_socket_thread_client!=NULL) {
            m_socket_thread_client->Delete();
            m_socket_thread_client = NULL;
        }
    }
}
void GridClient::Delete()
{
    m_delete = true;
}
void GridClient::DeleteMyself()
{
    delete this;
}
void GridClient::setAutoReconnectWhenConnectionLost(bool reconnect)
{
    m_reconnect_when_connection_lost = reconnect;
}
wxString GridClient::GetIP()
{
    return m_IP;
}
bool GridClient::sendPortNumber(wxSocketBase *sock, int port)
{
    sock->SetFlags(wxSOCKET_WAITALL);
    char msg[10];
    stringstream ss;
    ss << std::left << std::setw(10) << std::setfill(' ') << port; 
    string formattedString = ss.str();    
    formattedString.copy(msg, 10);
    sock->Write(msg, 10);
    if(sock->LastWriteCount()!=10) {
        WriteLogMessage("ERROR: last Write counts: "+to_string(sock->LastWriteCount()) +" != " + to_string(10));
        return false;
    }
    if(sock->Error()==true) {
        WriteLogMessage("ERROR: last Write error: "+ to_string(sock->LastError()));
        return false;
    }

    return true;
}
bool GridClient::ConnectClient(int nbOfTrial, wxString hostname, int hostport, int my_server_port)
{
   WriteLogMessage(_T("Client try to connect"));
   wxIPV4address ip;
   
   m_port = hostport;
   m_hostname = hostname;
   if(!ip.Service(m_port)) return false;
   if(!ip.Hostname(m_hostname)) return false;

   wxMutexLocker l(m_socket_thread_client_mutex);
   if(m_socket_thread_client!=NULL) { 
        WriteLogMessage("ERROR: The thread is still working, m_socket_thread!=0, cannot start client again"); 
        return false;                        
   }

   wxSocketClient *client = new wxSocketClient();
   client->SetEventHandler(*this, GRID_CLIENT_SOCKET_ID);   
   client->SetNotify(wxSOCKET_LOST_FLAG);
      
   client->Notify(true);
   client->SetTimeout(60);

   int i=0;
   do
   {
      i++;
      client->Connect(ip, false);
      client->WaitOnConnect(3,0);
      if(i==nbOfTrial) break;
   } while(!client->IsConnected());
   
    if (client->IsConnected()) {

        if(my_server_port>0) {
            if(!sendPortNumber(client, my_server_port)) {
                client->Destroy();
                return false;
            }
        }

        client->SaveState();
        //better to run it just once and then start it again after all precedures are done
        m_sendingTimer->Start(3*1000, true);
        WriteLogMessage(_T("Client is connected to the server"));                

        wxIPV4address tmp;
        WriteLogMessage("Creating Thread");
        
        m_socket_thread_client = new SocketThreadClient(client,
                                           m_working_dir, 
                                           this,
                                           tmp,
                                           &m_messages_to_be_send,
                                           &m_messages_to_be_send_mutex);

        if (m_socket_thread_client->Create() == wxTHREAD_NO_ERROR) {           
            WriteLogMessage("Thread is running");
            m_socket_thread_client->Run();
        } else {
            WriteLogMessage("ERROR: Cannot create next thread");
            m_socket_thread_client->Delete();
            m_socket_thread_client=NULL;
        }        
      return true;
   }
   else {
      client->Destroy();
      return false;
   }
}
void GridClient::OnSocketEvent(wxSocketEvent &event)
{
   WriteLogMessage("OnSocketEvent Begin");
   wxSocketBase *tmpSock = event.GetSocket();

   switch(event.GetSocketEvent())
   {
       //this should never happen. All this is realized in thread
     case wxSOCKET_INPUT:
        {
           WriteLogMessage("INPUT");
           break;
        }
     case wxSOCKET_LOST:
        {
           WriteLogMessage("Connection lost");
           break;
        }
     case wxSOCKET_OUTPUT:
         {
             WriteLogMessage(_T("OUTPUT"));
             break;
         }
     default:
       break;
   }
   
   WriteLogMessage("OnSocketEvent end");
}

bool GridClient::isConnected()
{
    wxMutexLocker l(m_socket_thread_client_mutex);
    if(m_socket_thread_client!=NULL) {
        return m_socket_thread_client->IsRunning();
    }
    return false;
}
void GridClient::refreshClientThreadState()
{    
    {
        wxMutexLocker l(m_socket_thread_client_mutex);
        if(m_socket_thread_client==NULL) return;
    
        //when thread is not running, delete it and reconnect
        if(!m_socket_thread_client->IsRunning()) {
            m_socket_thread_client->Delete();
            m_socket_thread_client = NULL;
            m_need_socket_connect = true;
        }
    }
}
void GridClient::CleanMessagesToBeSent()
{//just keep a few last messages to allow confirmation of delivery but not to consume a lot of memory
    wxMutexLocker locker(m_messages_to_be_send_mutex);
    for(int i=0;i<m_messages_to_be_send.size();i++) {
        if(m_messages_to_be_send.size() < m_keep_at_least_last_m_messages_to_be_send) return;
        if((m_messages_to_be_send[i].delivery_confirmation_obtained > 0) && (m_messages_to_be_send[i].delivery_confirmation_obtained > m_keep_minutes_m_messages_to_be_send)) {
            m_messages_to_be_send.erase(m_messages_to_be_send.begin()+i);
            i--;
        }
    }
}
void GridClient::OnTimerEvent(wxTimerEvent& event)
{
    WriteLogMessage("OnTimerEvent");        
    
    if(m_reconnect_when_connection_lost && m_need_socket_connect) {
        WriteLogMessage("Socket reconnection detected");    
        if(ConnectClient(1, m_hostname, m_port)) {
            m_need_socket_connect = false;
        }        
    } else {
        refreshClientThreadState();
    }

    CleanMessagesToBeSent();

    if(m_delete) {
        DeleteMyself();
        return;
    } 

    m_sendingTimer->Start(1*1000, true);
    
}
long long GridClient::addMessageToSend(wxString msg, long long msgID)
{
    wxMutexLocker locker(m_messages_to_be_send_mutex);
    
    MSGINFO_SENT mi;
    mi.ID = msgID;
    mi.msg = msg;
    m_messages_to_be_send.push_back(mi);
    return mi.ID;
}
void GridClient::printMsgs()
{
    wxMutexLocker locker(m_messages_to_be_send_mutex);

    wxString outmsg ="\n----------------\nm_messages_to_be_send:\n";
    for(int i=0;i<m_messages_to_be_send.size();i++) {
        outmsg+= "["+to_string(i)+"]:\n";
        outmsg+= "ID: "+to_string(m_messages_to_be_send[i].ID)+"\n";
        outmsg+= "sent: "+to_string(m_messages_to_be_send[i].sent)+"\n";
        outmsg+= "receiving confirmed: "+to_string(m_messages_to_be_send[i].delivery_confirmation_obtained)+"\n";        
        outmsg+= "msg len: "+to_string(m_messages_to_be_send[i].msg.length())+"\n";
    }    
    outmsg+="------------\n";
    WriteLogMessage(outmsg);
}
vector<GridCommunication::MSGINFO_SENT> GridClient::getCopyMsgsToBeSent()
{
    wxMutexLocker locker(m_messages_to_be_send_mutex);    
    return m_messages_to_be_send;
}
bool GridClient::getCopyMsgToBeSent(long long msgID, GridCommunication::MSGINFO_SENT &msg)
{
    wxMutexLocker locker(m_messages_to_be_send_mutex);
    for(int i=m_messages_to_be_send.size()-1;i>=0;i--) {
        if(m_messages_to_be_send[i].ID == msgID) {
            msg = m_messages_to_be_send[i];
            return true;
        }
    }
    return false;
}
