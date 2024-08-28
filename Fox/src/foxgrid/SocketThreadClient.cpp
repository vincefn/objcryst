#include "SocketThreadClient.h"

SocketThreadClient::SocketThreadClient(  wxSocketBase         *pSocket, 
                             wxString              workingDir, 
                             wxEvtHandler         *parent,
                             wxIPV4address         address,
                             vector<MSGINFO_SENT>       *messages_out,
                             wxMutex                    *messages_out_mutex)
:wxThread(wxTHREAD_JOINABLE )
{
    m_socket = pSocket;
    //Notify() cannot be called in thread context. We have to detach from main loop
    //before switching thread contexts.
    m_socket->Notify(false);
    m_socket->SetFlags(wxSOCKET_WAITALL|wxSOCKET_BLOCK);
    m_working_dir = workingDir;
    m_parent = parent;
    m_address = address;
    m_socket_error = false;
    m_messages_out = messages_out;
    m_messages_out_mutex = messages_out_mutex;
}
long SocketThreadClient::GetId()
{
   return wxThread::GetId();
}
wxIPV4address SocketThreadClient::GetAddress()
{
    return m_address;
}
wxThread::ExitCode SocketThreadClient::Entry()
{
    WriteLogMessage("Thread Entry");    
    while (!TestDestroy()) { //Connection lost event is handled by owner!
        if(!m_socket->IsConnected()) {
            break;
        }
     
        wxThread::Yield(); // this is important to call it before WaitForRead()
        wxThread::Sleep(100);
        if(isMessageToBeSent()) {
            WriteLogMessage("Something to send...");
            SndMsg();
        }
        if(m_socket_error) {
            break;
        }
    }

    m_socket->Destroy();
    m_socket = 0;
    WriteLogMessage("Leaving thread");
    return 0;
}
void SocketThreadClient::WriteLogMessage(wxString msg)
{
#ifdef ALLOW_GRID_LOGS
   wxString name;
#ifdef WIN32
   name = GetWorkingDir() + _T("\\") + "Thread_client_" + to_string(GetId())+".txt";
#else
   name = GetWorkingDir() + _T("/") + filename;
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
bool SocketThreadClient::SndMsg()
{//send data may take some time, so lock mutex only for finding the msg and changing the msg status. 
 //it does not lock the mutex during sending the data

    MSGINFO_SENT ms;
    bool found = false;

    WriteLogMessage("SndMsg() - start");
    { //searching the msg
        wxMutexLocker locker(*m_messages_out_mutex);
        if (!locker.IsOk()) {
            WriteLogMessage("ERROR: SndMsg() - !locker.IsOk()");
            return false;
        }
        for(int i=0;i<m_messages_out->size();i++) {
            if((*m_messages_out)[i].sent==0) { 
                ms = (*m_messages_out)[i];
                found = true;
                break;
            }
        }
    }
    //this is not under mutex - this is what we want
    short er;
    if(found) {
        //if found send it
        er = SendData2(m_socket, ms.ID, ms.msg.c_str(), ms.msg.length(), true);
   
        //change the state of the msg based on the resulting er
        wxMutexLocker locker(*m_messages_out_mutex);
        if (!locker.IsOk()) {
            WriteLogMessage("ERROR: SndMsg() - !locker.IsOk() 2");
            return false;
        }
        for(int i=0;i<m_messages_out->size();i++) {
            if((*m_messages_out)[i].ID==ms.ID) { 
                if(er==0) {
                    (*m_messages_out)[i].sent = getTimeStampMinutes();
                    (*m_messages_out)[i].delivery_confirmation_obtained = (*m_messages_out)[i].sent;
                    (*m_messages_out)[i].msg = ""; // empty the memory
                    WriteLogMessage("Send");
                } else {
                    (*m_messages_out)[i].sent = 0;
                    (*m_messages_out)[i].delivery_confirmation_obtained = 0;
                    WriteLogMessage("ERROR: "+to_string(er));
                    m_socket_error = true;
                    return false;
                }
                break; 
            }
        }
    }
                        

    WriteLogMessage("SndMsg() - end");
    return true;
}

bool SocketThreadClient::isMessageToBeSent()
{
    wxMutexLocker locker(*m_messages_out_mutex);
    if (!locker.IsOk()) {
        WriteLogMessage("ERROR: isMessageToBeSent() - !locker.IsOk()");
        return false;
    }

    for(int i=0;i<m_messages_out->size();i++) {
        if((*m_messages_out)[i].sent==0) {
            return true;
        }
    }
    return false;
}


