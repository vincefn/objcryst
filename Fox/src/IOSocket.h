#ifdef __WX__CRYST__
   // For compilers that support precompilation, includes "wx/wx.h".
   #ifndef __DARWIN__ // work around MacOSX type_info bug (??)
      #include "wx/wxprec.h"
   #endif

   #ifdef __BORLANDC__
       #pragma hdrstop
   #endif

   // for all others, include the necessary headers (this file is usually all you
   // need because it includes almost all "standard" wxWindows headers)
   #ifndef WX_PRECOMP
       #include "wx/wx.h"
   #endif

   //#include "wx/tooltip.h"
   //#include "wx/notebook.h"
   #include "wx/wfstream.h"
   #include "wx/zstream.h"
   #include "wx/fileconf.h"
   #include "wx/socket.h"
   #include "wx/process.h"
   #include "wx/sckstrm.h"
   #include "wx/sstream.h"
   #include "wx/zstream.h"
   #include "wx/wfstream.h"
   #include "wx/thread.h"
   #include "wx/stream.h"
   #include "wx/dir.h"
   //#include "wx/dynarray.h"
#endif

#ifndef __IO_SOCKET__
#define __IO_SOCKET__

class IOSocket
{
public:
    IOSocket(void);
    ~IOSocket(void);

    //reads message from the socket
    //returns true if successful, otherwise returns false
    //use getError() message to get error details;
    //receipt - send receipt?
    bool ReadStringFromSocket(wxSocketBase *pSocket, std::string &message, bool receipt=true);

    //writes message to the socket
    //returns true if successful, otherwise returns false
    //use getError() message to get error details;
    //receipt - send receipt?
    bool WriteStringToSocket(wxSocketBase *pSocket, std::string s, bool receipt=true);

    //returns error message
    //Please note that this function merely returns the last error message,
    //but it should not be used to determine if an error has occurred
    wxString getError();
private:

    unsigned int getMessageLen(wxSocketBase *pSocket);

    wxString m_error;

};

#endif
