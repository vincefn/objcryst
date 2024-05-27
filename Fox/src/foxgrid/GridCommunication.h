#ifndef __GRIDCOMM_H // Make sure to only declare these classes once
#define __GRIDCOMM_H

#include <wx/utils.h>
#include <chrono>
#include <sstream>
#include <iomanip>
#include <wx/socket.h>
#include <wx/file.h>
#include <wx/sckstrm.h>
#include <wx/sstream.h>
#include <wx/wfstream.h>
#include <wx/zstream.h>
#include <wx/txtstrm.h>
#include <array>

#include "ObjCryst/wxCryst/wxCryst.h" //because WXCRYST_ID


using namespace std;
using namespace ObjCryst; //becaouse WXCRYST_ID

// Determine the platform and set appropriate pipe functions
#if defined(_WIN32) || defined(_WIN64)
#define POPEN _popen
#define PCLOSE _pclose
#else
#define POPEN popen
#define PCLOSE pclose
#endif

#define _ALLOW_GRID_LOGS



static const long GRID_SERVER_ID=                  WXCRYST_ID();
static const long GRID_CLIENT_SOCKET_ID=                   WXCRYST_ID();
static const long ID_SEND_TIMER=                         WXCRYST_ID();
static const long ID_GRID_MASTER_CHECK_SERVER=                         WXCRYST_ID();
static const long ID_GRID_MASTER_CHECK_SLAVES=                         WXCRYST_ID();
static const long ID_GRID_MASTER_CHECK_SLAVE=                         WXCRYST_ID();
static const long ID_GRID_CHECK_SLAVE=                         WXCRYST_ID();

struct SocketThreadInfo {
    //long ID;
    wxIPV4address address;
    int port;

    SocketThreadInfo() {
        port = -1;
    }

};
struct SentMsgInfo {
    long long msgID;
    int nb_runs;
    wxString SlaveIP;

    SentMsgInfo() {
        msgID = 0;
        nb_runs = 0;
    }
    SentMsgInfo(long long ID, int runs, wxString IP) {
        msgID = ID;
        nb_runs = runs;
        SlaveIP = IP;
    }
};
struct MasterResult {    
    wxString data;
    wxString filename;
    float Cost;
    float Rwp;
    long JobID;
    MasterResult () {
        Cost = -1.0;
        Rwp = -1.0;
        JobID = -1;
    }
    wxString generateFileName(long ID, float Cost, float Rwp)
    {
        int t = time(0);
        filename = wxString::Format("ID-%ld_Cost-%.1f_Rwp-%7.5f_Time-%d.xml", ID, Cost, Rwp, t);
        return filename;
    }
    bool operator<(const MasterResult &m) const
   {
	    return (Cost < m.Cost);
   }
};
/*
struct SlaveJob {
    long ID;
    wxString data;
    int nb_runs;
    int nb_trial;
    int nb_running;
    int nb_done;
    bool randomize;
    wxString name;
    wxTimeSpan average_calc_time;

    SlaveJob() {
        ID = 0;
        nb_runs = 0;
        nb_trial = 0;
        nb_running = 0;
        nb_done = 0;
        randomize = true;        
        average_calc_time = wxTimeSpan();
    }
    int getNbFreeJobs() {
        return (nb_runs - nb_done - nb_running);
    }
};
*/
struct MasterJob {
    bool       deleted;
    long       ID; 
    wxString   data;
    int        nb_runs;
    vector<SentMsgInfo> sentToClientsMsgInfo; //msgIDs that were sent to clients with this job and receiving was not confirmed yet and nb of running expected    
    vector<MasterResult> results;
    int        nb_trial;
    vector<wxString> runningIPs;
    bool       randomize;
    wxString   name;
    wxString   filename;
    wxTimeSpan average_calc_time;
    long long  msgID;

    MasterJob () {
        deleted = false;
        ID = 0;
        nb_runs = 0;
        //nb_done = 0;
        nb_trial = 0;
        //nb_running = 0;
        randomize = true;
        average_calc_time = wxTimeSpan();
        msgID = 0;
    }   
    int getNbSentToClients() {
        int sum = 0;        
        for (const auto& p : sentToClientsMsgInfo) {
            sum += p.nb_runs;
        }
        return sum;
    }
    int getNbFreeJobs() {
        return (nb_runs - results.size() - runningIPs.size() - getNbSentToClients());
    }

    //removes just one ore several of them from runningIPs or from sentToClientsMsgInfo
    //use nb=-1 to delete all from runningIPs and also from sentToClientsMsgInfo
    bool removeIPFromRunningIPs(wxString IP, int nb=1) {         
        for(int i=0;i<runningIPs.size();i++) {            
            if(runningIPs[i].compare(IP)==0) {
                runningIPs.erase(runningIPs.begin()+i);
                i--;
                nb--;
                if(nb==0) return true;
            }
        }
        //in the case that the IP was not found in runninIPs or there are still some leftovers, remove IP from sentToClientsMsgInfo
        for(int i=0;i<sentToClientsMsgInfo.size();i++) {
            if(sentToClientsMsgInfo[i].SlaveIP.compare(IP)==0) {
                sentToClientsMsgInfo[i].nb_runs--;
                if(sentToClientsMsgInfo[i].nb_runs==0) {
                    sentToClientsMsgInfo.erase(sentToClientsMsgInfo.begin()+i);
                    i--;
                }
                nb--;
                if(nb==0) return true;
            }
        }

        //not sure, what to do if this happen...:(
        return false;
    }


};
struct ResultInfo {
    long          resultID;
    long          jobID;
    wxString      data;
    long          duration_seconds;
    bool          sent;
	bool          pending;
    long long     msgID;
            
    ResultInfo() {
        resultID = -1; 
        jobID = -1;
        duration_seconds = -1;
        sent = false;
        pending = false;
        msgID = -1;
    }
};

int generateUniqueID();

class GridCommunication
{
    public:
        GridCommunication();
        ~GridCommunication();

        char CalculateXORChecksum(const string& data);
        string GenerateUniqueIdentifier(const string& userSpecificString);

        static long long getTimeStampNanoSeconds();
        unsigned int getTimeStampMinutes();
        static long long getTimeStampSeconds();

        static vector<int> getUsedPorts();

        short SendData(wxSocketBase *socket, long long msgID, const char* data, wxUint32 const dataLen);
        short ReadData(wxSocketBase *socket, long long &msgID, vector<char> &data);

        short SendDataThread(wxSocketBase *socket, long long msgID, wxString text);
        short ReceiveDataThread(wxSocketBase *socket, long long &msgID, wxString &text);

        short SendData2(wxSocketBase *socket, long long msgID, const char* data, wxUint32 const dataLen, bool thread);
        short ReceiveData2(wxSocketBase *socket, long long &msgID, vector<char> &data, bool thread);

        void WriteLogMessage(wxString msg, wxString filename="log.txt");
        wxString GetWorkingDir();
        wxString m_working_dir;

    protected:

        bool lastReadOK(wxSocketBase *socket, wxUint32 len);
        bool lastWriteOK(wxSocketBase *socket, wxUint32 len);

    public:
        struct MessageHeader {
            wxUint32 chunkLength; //length of the message after the header, if 0 it is the information about the msgID and msgLen
            char chunkChecksum;
            wxUint32 msgLen;
            long long msgID;
        };

        struct MSGINFO_REC {
            long long ID;// Message ID
            //long socketID;

            //IP of the sender 
            wxString IP;
            wxString msg;            
            unsigned int recieved; //time stamp in minutes, when it was delivered
            unsigned int delivery_confirmation_sent; //time stamp in minutes, when the delivery confirmation was sent.
            unsigned int processed; //time stamp in minutes, when the message was proccessed or copied somewhere else.
            wxSocketBase* m_socket;

            MSGINFO_REC() {
                ID = 0;                
                recieved = 0;
                delivery_confirmation_sent = 0;
                processed = 0;
                m_socket = 0;
                //socketID = 0;
            }
        };
        struct MSGINFO_SENT {
            long long ID; // Message ID
            wxString msg;
            unsigned int sent; //time stamp in minutes, when the msg was sent            
            unsigned int delivery_confirmation_obtained; //time stamp in minutes, when delivery confirmation was received                        

            MSGINFO_SENT() {
                ID = 0;
                sent = 0;
                delivery_confirmation_obtained = 0;                
            }
        };
        
        
};
#endif

