#ifdef __WX__CRYST__
   #include "ObjCryst/wxCryst/wxCrystal.h"
#endif

using namespace ObjCryst;
using namespace std;

#include "IOSocket.h"

IOSocket::IOSocket(void)
{
}

IOSocket::~IOSocket(void)
{
}
bool IOSocket::ReadStringFromSocket(wxSocketBase *pSocket, std::string &message, bool receipt)
{
    //clearing error message
    m_error.Clear();

    VFN_DEBUG_MESSAGE(__FUNCTION__<<":"<<message,10)
    //set one minute to wait
    pSocket->SetTimeout(20);

    //the same as in writestringtosocket(...)
    if (pSocket->Error()) {
       m_error << _T("ReadStringFromSocket: Previous error detected: ") << (int) pSocket->LastError() << _T(", ") << pSocket->GetLastIOSize() <<_T(", ") << pSocket->LastCount();
       m_error << _T(" try to continue\n");
    }

    //get length
    unsigned int len = getMessageLen(pSocket);

    //alloc memory for message
#if 1 //def __WIN32__
    char *buf;
    buf = (char *) calloc(len+1, sizeof(char) );
    if(buf==NULL)  {
        m_error << _T("Can't alloc memory\n");
        return false;
    }
#else
    char buf[len+1];
#endif

    //try to read message
    pSocket->ReadMsg(buf,len);
    if (pSocket->Error()) {
        pSocket->Discard();
        m_error << _T("ReadStringFromSocket error: ") << (int) pSocket->LastError() << _T(", ") << pSocket->GetLastIOSize() <<_T(", ") << pSocket->LastCount();
        m_error << _T(" (sending message)\n");
        #if 1//def __WIN32__
            free(buf);
        #endif
        return false;
    }
   /*
    if(receipt) {
        //send receipt
        wxString sreceipt;
        sreceipt << (int) len;
        //wait for write
        if( !pSocket->WaitForWrite(10)) {
            m_error << _T("ReadStringFromSocket error: WaitForWrite return false (sending receipt)\n");
            #ifdef __WIN32__
                free(buf);
            #endif
            return false;
        }
        if( !WriteStringToSocket(pSocket, string(sreceipt.ToAscii()), false) ) {
            m_error << _T("ReadStringFromSocket error: sending receipt return false (sending receipt)\n");
            #ifdef __WIN32__
                free(buf);
            #endif
            return false;
        }
    }
*/
    //saving message
    buf[len]='\0';
    message=string(buf);
    VFN_DEBUG_MESSAGE(__FUNCTION__<<":"<<message,10)

    //free memory
#if 1 //def __WIN32__
    free(buf);
#endif

    return true;
}
bool IOSocket::WriteStringToSocket(wxSocketBase *pSocket, std::string s, bool receipt)
{
    //clearing error message
    m_error.Clear();

    pSocket->SetTimeout(20);

    //inform about previous error...
    if (pSocket->Error()) {
       m_error << _T("WriteStringToSocket: Previous error detected: ") << (int) pSocket->LastError() << _T(", ") << pSocket->GetLastIOSize() <<_T(", ") << pSocket->LastCount();
       m_error << _T(" try to continue\n");
    }
    VFN_DEBUG_MESSAGE(__FUNCTION__<<":"<<s,10)

    //send message
    pSocket->WriteMsg((void*) s.c_str(),s.size());
    if (pSocket->Error()) {
       m_error << _T("WriteStringToSocket error: ") << (int) pSocket->LastError() << _T(", ") << pSocket->GetLastIOSize() <<_T(", ") << pSocket->LastCount();
       m_error << _T(" (sending message)\n");
       return false;
    }

    if(!receipt) return true;
/*
    //wait for receipt
    pSocket->WaitForRead(10);
    if(!pSocket->IsData()) {
        m_error << _T("WriteStringToSocket error: waiting for receipt - timeout\n");
        return false;
    }
    //read receipt
    string sreceipt;
    long lreceipt;
    if(!ReadStringFromSocket(pSocket, sreceipt, false)) {
        m_error << _T("WriteStringToSocket error: reading receipt return false\n");
        return false;
    }
    //convert to a number and compare
    wxString tmp = wxString::FromAscii(sreceipt.c_str());
    tmp.ToLong((long *) &lreceipt);
    if(lreceipt!=s.size()) {
        m_error << _T("WriteStringToSocket error: returned size does not match ") << (int) s.size() <<_T(" != ") << (int) lreceipt << _T("\n");
        return false;
    }
*/
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
