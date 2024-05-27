#ifndef __FOXGRIDMASTER_H // Make sure to only declare these classes once
#define __FOXGRIDMASTER_H
#include "GridMasterBase.h"

#ifdef __WX__CRYST__
   #include "ObjCryst/wxCryst/wxCrystal.h"
#endif


using namespace ObjCryst;
using namespace std;

class FoxGridMaster: public GridMasterBase, wxEvtHandler
{
    public:
        struct SLAVE_FOX_INFO {
                
                GridMasterBase::SLAVE_INFO_PUBLIC sinfo;                
                unsigned short nb_CPU_all;
                unsigned short nb_CPU_idle;
                bool connected;
                long long LastProcessedJobMsgID;
                long long LasSentJobMsgID;
                int port;
            
                SLAVE_FOX_INFO() {
                    nb_CPU_all = 0;
                    nb_CPU_idle = 0;
                    connected = false;
                    LastProcessedJobMsgID = 0;
                    LasSentJobMsgID = 0;
                    port = -1;
                }
            };

        FoxGridMaster(wxString working_dir);
        ~FoxGridMaster();

        void OnCheckSlavesTimerEvent(wxTimerEvent& event);
        bool addJob(wxString filename, wxString name, wxString data, int id, long nbOfTrial, long nbRun, bool rand);
        int  generateJobID();
        void AskSlavesState();
        void ProcessMessagesFromSlaves();
        //vector<wxString> getResults();
        bool SendSomeJob(wxString SlaveIP, int freeCPUs);
        void updateSlaveCPUsInfo(wxString SlaveIP, int CPUfree, int CPUall, long long lastProcessedJobMsgID, long long &lastSentJobMsgID);
        void updateSlaveLastJobMsgsInfo(wxString SlaveIP, long long lastSentJobMsgID);
        vector<MasterJob> getJobs();
        bool DeleteJob(long JobID);
        bool UpdateJobHeader(MasterJob mj);
        bool existsJob(int ID);
        wxString createJobHeader(wxString name, int ID, int nb_trial, int nb_runs, int randomize);
        void CheckIfJobsReceived();

        vector<SLAVE_FOX_INFO> FoxGridMaster::getFoxSlaveInfo();

    protected:

        wxString getResult(wxString message, long pos);
        bool SaveResult(wxString filename, wxString data);
        void SlaveDisconnected(wxString SlaveIP);

        wxTimer                 *m_checkSlaveStateTimer;
        vector<MasterJob>        m_jobs;
        wxMutex                  m_jobs_mutex;
        unsigned long            m_timer_iter;
        //vector<wxString>         m_results;
        //wxMutex                  m_results_mutex;
        vector<SLAVE_FOX_INFO>   m_fox_slaves_info;
        wxMutex                  m_fox_slaves_info_mutex;

        DECLARE_EVENT_TABLE()
};
#endif

