#ifdef __WX__CRYST__
   #include "ObjCryst/wxCryst/wxCrystal.h"
#endif

using namespace ObjCryst;
using namespace std;

#include "IOSocket.h"
#include <wx/stopwatch.h>

MyMutexHelper::MyMutexHelper(wxMutex *m, long milliseconds, int &nb)
{
    m_mutex = m;
    wxStopWatch sw;
    nb = 0;
    do {
        nb++;
        m_error = m_mutex->Lock();
        if(m_error == wxMUTEX_NO_ERROR) return;
        wxSleep(1);
        if(nb>10) return;
    } while(sw.Time()!=milliseconds);
}
MyMutexHelper::~MyMutexHelper()
{
    m_mutex->Unlock();
}
wxMutexError MyMutexHelper::getError()
{
    return m_error;
}
wxMutexError MyMutexHelper::Unlock()
{
    return m_mutex->Unlock();
}


IOSocket::IOSocket(void)
{
}

IOSocket::~IOSocket(void)
{
}
bool IOSocket::ReadStringFromSocket(wxSocketBase *pSocket, wxString &message)
{
    //m_log = "";
    //m_log = "ThreadWorker: Reading the socket ...";
    pSocket->SetFlags(wxSOCKET_WAITALL);
    unsigned int len;
    pSocket->Read(&len, sizeof(int));
    if (pSocket->Error()) {
        m_error = "ThreadWorker: Read error";
        m_log += m_error;
        return false;
    }
    int processed = pSocket->LastCount();
    //m_log += wxString::Format("ThreadWorker: %d bytes read in the header (len = %d)", processed, len);           

    if (len == 0) {
        //m_log += "ThreadWorker: 0 bytes in socket, nothing to read...";
        return true;
    }

    wxCharBuffer buf(len);
    //m_log += wxString::Format("Message header was: size = %d (bytes)",len);
    //WriteLogMessage("ThreadWorker: Reading message ...");
    pSocket->Read(buf.data(), len);

    if (pSocket->Error())
    {
        //WriteLogMessage("ThreadWorker: Read error");
        //wxGetApp().AddPendingEvent(e);
        return false;
    }
    processed = pSocket->LastCount();
    //m_log += wxString::Format("ThreadWorker: %d bytes readed", processed);
    wxString tmp(buf);
    message = tmp;
    return true;
 }
bool IOSocket::WriteStringToSocket(wxSocketBase *pSocket, wxString msg)
{
      pSocket->SetFlags(wxSOCKET_WAITALL);
     
      wxCharBuffer buffer=msg.ToAscii();
      // Note that len is in bytes here!
      unsigned int len = strlen(buffer.data()) * sizeof(char);
      //WriteLogMessage(wxString::Format("ThreadWorker: Sending header of the message of %d kilobytes", len));
      pSocket->Write(&len, sizeof(int));
      if (pSocket->Error()) {
         //WriteLogMessage("ThreadWorker: Write error");
         return false;
      }
      //WriteLogMessage("ThreadWorker: Sending the message ...");
      pSocket->Write(buffer.data(), len);
      if (pSocket->Error()) {
         //WriteLogMessage("ThreadWorker: Write error");
         return false;
      }
      //WriteLogMessage(m_socket->Error() ? _("ThreadWorker: failed !\n") : _("ThreadWorker: done\n"));
      return true;
}
unsigned int IOSocket::getMessageLen(wxSocketBase *pSocket)
{
    wxUint32 len, sig;

    struct
    {
        unsigned char sig[4];
        unsigned char len[4];
    } msg;

    //peek socket
    pSocket->Peek(&msg, sizeof(msg));
    if (pSocket->Error()) {
       m_error << _T("getMessageLen error: ") << (int) pSocket->LastError() << _T(", ") << pSocket->GetLastIOSize() <<_T(", ") << pSocket->LastCount();
       m_error << _T(" (Peeking)\n");
       return -1;
    }

    sig = (wxUint32)msg.sig[0];
    sig |= (wxUint32)(msg.sig[1] << 8);
    sig |= (wxUint32)(msg.sig[2] << 16);
    sig |= (wxUint32)(msg.sig[3] << 24);

    if (sig != 0xfeeddead)
    {
        m_error << _T("Invalid signature\n");
        return -1;
    }

    len = (wxUint32)msg.len[0];
    len |= (wxUint32)(msg.len[1] << 8);
    len |= (wxUint32)(msg.len[2] << 16);
    len |= (wxUint32)(msg.len[3] << 24);

    return len;
}
wxString IOSocket::getError()
{
    return this->m_error;
}
