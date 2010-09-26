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

   #include "wx/tooltip.h"
   #include "wx/notebook.h"
   #include "wx/wfstream.h"
   #include "wx/zstream.h"
   #include "wx/fileconf.h"
   #include "wx/socket.h"
   #include "wx/dynarray.h"
   #include "wx/thread.h"
   #include "wx/sckstrm.h"
   #include "wx/sstream.h"
   #include "wx/zstream.h"
   #include "wx/wfstream.h"

    
#endif

#ifndef __FOX_THREAD__
#define __FOX_THREAD__ 
#include "GridResult.h"
#include "FoxJob.h"
#include "IOSocket.h"

#define __SERVER_LOGS 1

enum FoxServerEvents {INPUT_MSG, LOST_CONNECTION, NEW_CONNECTION, SEND_JOB, EMPTY_MSG};
enum ServerThreadStatus {FG_N_A, FG_CONNECTED, FG_EXPECTING_ANSWER, FG_EXPECTING_RESULT};

class FoxServerThread : public wxThread
{
////////////////////////////////////////////////////////////////////////////////////////
//Set the Socket->SetNotify(wxSOCKET_LOST_FLAG) and wxMutex before you run this thread//  
////////////////////////////////////////////////////////////////////////////////////////

public:

   FoxServerThread(   wxSocketBase* pSocket,
                  FoxServerEvents evt,
                  wxMutex *pMutex,
                  std::vector<GridResult > *pResults,
                  std::vector<FoxJob > *pJobs);              

   ~FoxServerThread();
   //bool GetJobID(int &ID);
   bool NewEvent(FoxServerEvents evt, wxMutex *MutexToUnlock);

   wxSocketBase* GetSocket();

   virtual void *Entry();
   virtual void OnExit(); 
   int GetId();
   wxString getName();
   long     getAvailCPUs();
   long     getAllCPUs();
   ServerThreadStatus getStatus();
    
private:
   void   CloseConnection();
   void   OnInput();
   bool   AnalyzeMessage(std::string message);
   void   SendJob(int nbOfJobs);
   bool   SendAsk(bool getClientInfo=false);   
   void   SaveDataAsFile(wxString out, wxString filename);
   bool   LoadFile(wxString filename, wxString &in);
   void   WriteLogMessage(wxString msg);
   void   SaveResult(wxString result, int JobID, float ResultCost);
   void   rejectJobs(vector<long> ids);
   wxString getResult(wxString message, long pos);

   wxSocketBase       * m_pSocket;
   wxMutex            * m_tMutexObj;
   wxMutex            * m_tThreadMutex;
   int                 m_No;
   FoxServerEvents     m_sckt_ntf;
   bool                m_exit;
   bool                m_newEvt;
   vector<GridResult > * m_results;
   vector<FoxJob >     * m_jobs;
   IOSocket             m_IOSocket;
   wxString             m_name;
   long                 m_allCPUs;
   long                 m_availableCPUs;   
   ServerThreadStatus   m_status;

};

#endif