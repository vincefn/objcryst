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
   #include "ObjCryst/RefinableObj/GlobalOptimObj.h"
   #include "wx/socket.h"
   #include "wx/dynarray.h"
   #include "wx/arrstr.h"
   #include "wx/list.h"
#endif

//#include <locale.h>
//#include <fstream>
//#include <sstream>
//#include <list>
#include <vector>

#define __FOX_SERVER__

#include "FoxJob.h"
#include "GridResult.h"
#include "IOSocket.h"

//#include "FoxServerThread.h"
//#include "wxCryst/wxCryst.h"
//DECLARE_EVENT_TYPE(SERVER_THREAD_EVENT, -1)
const wxEventType wxEVT_WORKER = wxNewEventType();

struct MY_STHREAD
{
    int             m_availableCPUs;
    int             m_freeCPUs;
    wxString        m_name;
    int             m_id;
};

class WorkerEvent : public wxEvent
{
public:
    WorkerEvent(void* pSender)
    {
        SetId(-1);
        SetEventType(wxEVT_WORKER);
        m_sender = pSender;
        m_exit = false;
        m_workerFailed = false;
    }

    virtual wxEvent* Clone() const
    {
        return new WorkerEvent(*this);
    }

    void*           m_sender;
    bool            m_exit;
    bool            m_workerFailed;
    wxString        m_msg;
    MY_STHREAD      m_thread_info;
};


//#include <wx/frame.h>
class ThreadWorker : public wxThread
{
public:
    ThreadWorker(wxSocketBase         *pSocket, 
                 wxString              workingDir, 
                 wxFrame              *parent,
                 vector<GridResult >  *results,
                 vector<FoxJob >      *jobs,
                 wxMutex              *mutexProtecting_Jobs_Results);
    virtual ExitCode Entry();
    int GetId();

private:
    bool analyze_message_and_get_answer(wxString msg, wxString &answer);
    //bool read_socket(wxSocketBase* m_socket, wxCharBuffer &message);
    //bool write_socket(wxSocketBase* m_socket, wxCharBuffer message);
    void WriteLogMessage(wxString msg);
    bool LoadFile(wxString filename, wxString &in);
    void SaveDataAsFile(wxString out, wxString filename);
    wxString getResult(wxString message, long pos);
    void SaveResult(wxString result, int JobID, float ResultCost);
    //this will return N of available jobs (or lower no.)
    vector<FoxJob > getJobsToCalculate(int n);
    void AddResultToJobList(int id);
    void RemoveActiveFromJobList(int id);
    void RemoveThreadFromJobList();
    void cannot_calculate_this_job(FoxJob fj);
    
    wxString             m_working_directory;
    wxSocketBase        *m_socket;
    wxIPV4address        m_peer;
    wxFrame             *m_parent;
    vector<GridResult > *m_global_results;
    vector<FoxJob >     *m_global_jobs;
    IOSocket             m_IO;
    wxMutex             *m_mutexProtecting_Jobs_Results;
    wxString             m_name;
    int                  m_availableCPUs;
    int                  m_freeCPUs;
};

//WX_DECLARE_LIST(ThreadWorker, TList);

class FoxServer: public wxFrame
{
   public:
     FoxServer();
     ~FoxServer();
     void StartGridServer();
     void SetWorkingDir(wxString path);
     wxString GetWorkingDir();
     void OnServerEvent(wxSocketEvent &event);
     //void OnThreadEvent(wxCommandEvent &event);
     //void OnSocketEvent(wxSocketEvent &event);
     //void WriteProtocol();
     void GetData(std::vector<GridClient> &Clients, std::vector<GridResult > &results,  std::vector<FoxJob > &Joblist);
     bool IsServerRunning();
     //void RunAllClients();
     void AddJobToList(FoxJob newjob);
     void UpdateJob(int index, FoxJob *cjob);
     int DeleteJob(int index);
     void OnWorkerEvent(WorkerEvent& pEvent);

   protected:

     void WriteLogMessage(wxString msg);

     //TList m_threadWorkers;

     //protecting jobs and results
     wxMutex            * m_mutexProtecting_Jobs_Results;

     //wxMutex            * mutexMessageLog;
     //wxMutex            * m_threadMutex;

     wxSocketServer     * mpServer;

     //here we save pointers to threads
     std::vector<ThreadWorker *>   m_threads;

     //all results protected by mutex
     std::vector<GridResult >      m_results;

     //all jobs proctected by mutex
     std::vector<FoxJob > m_jobs;

     //basic info for GUI, can access without mutex
     vector<MY_STHREAD>   m_threads_info;

     bool                 m_isRunning;
     //bool                 m_needUpdate;
     wxString             m_working_dir;
     DECLARE_EVENT_TABLE()
};
