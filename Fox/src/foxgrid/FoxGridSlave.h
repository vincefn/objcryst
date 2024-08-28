#ifndef _FOXGRIDSLAVE_H 
#define _FOXGRIDSLAVE_H
#include "GridSlaveBase.h"

#ifdef __WX__CRYST__
   #include "ObjCryst/wxCryst/wxCrystal.h"
   #include "wx/process.h"
   #include "wx/filename.h"
   #include "wx/dir.h"
#endif


using namespace ObjCryst;
using namespace std;

const wxEventType wxEVT_PROCESS_MY = wxNewEventType();

class ProcessMyEvent : public wxEvent
{
public:
    ProcessMyEvent(void* pSender)
    {
        SetId(-1);
        SetEventType(wxEVT_PROCESS_MY);
        m_sender = pSender;
        m_exit = false;
        m_pid = -1;
        m_status = -1;
    }

    virtual wxEvent* Clone() const
    {
        return new ProcessMyEvent(*this);
    }

    void*           m_sender;
    int             m_pid;
    int             m_status;
    wxString        m_dir;
    bool            m_exit;
};

struct FOX_PROCESS {
    int      pid;
    wxString dir;
    bool     running;
    int      jobID;
    wxDateTime startingtime;

    FOX_PROCESS() {
        pid = 0;
        running = false;
        jobID = 0;
    }
    int getProgressInPercents(wxTimeSpan avCalcTime) {
        if(avCalcTime == wxTimeSpan()) return -1;
        wxDateTime ct = wxDateTime::Now();
        wxTimeSpan duration = ct - startingtime;
        int p = (int) (100*(duration.GetSeconds().ToDouble() / avCalcTime.GetSeconds().ToDouble()));
        if(p>100) p=100;
        return p;
    }

};

class FoxGridSlave: public GridSlaveBase, public wxEvtHandler
{
    public:
        FoxGridSlave(wxString m_working_dir);
        ~FoxGridSlave();

        void OnCheckSlaveTimerEvent(wxTimerEvent& event);
        void OnProcessEvent(ProcessMyEvent& pEvent);               

        vector<MasterJob> getJobs();        
        
        //Sets nb of available CPUs
        //Use it carefully, it kills all current processes and allocate them again!
        void ResetNbCPUsAll(int CPUs);      

        //void SendResult(wxString result, long duration, long jobID, long resultID);
        
        int getCPUsAll();
        int getCPUsFree();
        wxString getMyStateMsg();
        wxString getMyHostname();

        vector<FOX_PROCESS> getProcesses();
        
    protected:

        void addJob(MasterJob job);
        void CheckJobsAndStartCalculation();
        void CheckResultsAndSendOne();       
        FOX_PROCESS *getUnusedProcess();
        bool SaveDataAsFile(wxString out, wxString filename);
        void setProcessUnused(int pid);
        bool SaveResult(wxString fileName, wxString Cost, wxString Rwp, int ID, wxDateTime startingtime, bool error);                
        bool findResultedFile(wxString dir, wxString &Rwp, wxString &Cost, wxString &filename);
        bool LoadFile(wxString filename, wxString &in);
        bool DeleteXMLFilesInDirectory(const wxString& directoryPath);

        //Sets nb of CPUs
        //Use it carefully, it kills all current processes and allocate them again!
        void resetProcesses(int nbProcesses);
        void ProcessMsgs(vector<GridCommunication::MSGINFO_REC> &msgs);
        wxString getJobData(wxString inmsg, long pos);

        wxTimer     *m_checkSlaveTimer;
        
        //vector<ResultInfo> m_results;  
        //wxMutex            m_results_mutex;
        
        vector<MasterJob>  m_jobs;
        wxMutex            m_job_mutex;
        long long          m_last_job_msg_processed;
        


        vector<FOX_PROCESS> m_processes;
        wxMutex             m_processes_mutex;



        DECLARE_EVENT_TABLE()
};

class MyProcess : public wxProcess
{
public:
    MyProcess(FoxGridSlave *parent, const wxString& cmd, wxString dir);
    virtual void OnTerminate(int pid, int status);

protected:
   wxString     m_cmd;
   FoxGridSlave  * m_parent;
   wxString     m_dir;
};

#endif