#include "FoxServer.h"

#ifdef __WX__CRYST__
   #include "ObjCryst/wxCryst/wxCrystal.h"
#endif


using namespace ObjCryst;
using namespace std;

// Fox client/server grid
//#include "wx/socket.h"
//#include <wx/file.h>
//#include <wx/filefn.h>

static const long GRID_SERVER_ID=                       WXCRYST_ID();
static const long GRID_SERVER_SOCKET_ID=                WXCRYST_ID();
static const long SERVER_THREAD_EXIT=                   WXCRYST_ID();

BEGIN_EVENT_TABLE(FoxServer, wxFrame)
   EVT_SOCKET(GRID_SERVER_ID,                       FoxServer::OnServerEvent)
   EVT_SOCKET(GRID_SERVER_SOCKET_ID,                FoxServer::OnSocketEvent)
END_EVENT_TABLE()

FoxServer::FoxServer():
mpServer(0)
{
   s_mutexProtectingTheGlobalData = new wxMutex();
   mutexMessageLog = new wxMutex();
   m_threadMutex = new wxMutex();
   m_needUpdate = false;
   m_isRunning = false;
   srand( (unsigned)time( NULL ) );
   m_working_dir = _T("");
}
FoxServer::~FoxServer()
{
   //todo: clear m_results
   delete mpServer;
   delete s_mutexProtectingTheGlobalData;
   delete mutexMessageLog;
}
void FoxServer::WriteLogMessage(wxString msg)
{
#if __SERVER_LOGS
   if(mutexMessageLog->Lock()!=wxMUTEX_NO_ERROR) return;
   wxString filename;
#ifdef WIN32
   filename = GetWorkingDir() + _T("\\server.log");
#else
   filename = GetWorkingDir() + _T("/server.log");
#endif
   wxFile *logfile = new wxFile(filename, wxFile::write_append);
   if(logfile != 0)
   {
      wxDateTime datetime = wxDateTime::Now();
      logfile->Write(datetime.Format(_T("%X ")) + msg + _T("\n"));
      logfile->Close();
   }
   delete logfile;
   mutexMessageLog->Unlock();
#endif
   (*fpObjCrystInformUser)(msg.ToStdString());
}
void FoxServer::SetWorkingDir(wxString path)
{
    m_working_dir = path;
}
wxString FoxServer::GetWorkingDir()
{
    return m_working_dir;
}
void FoxServer::StartGridServer()
{
   VFN_DEBUG_MESSAGE(__FUNCTION__,10)
   wxIPV4address ip;
   ip.Service(2854);

   WriteLogMessage(_T("Starting server..."));
   mpServer = new wxSocketServer(ip,wxSOCKET_REUSEADDR);
   mpServer->SetEventHandler(*this, GRID_SERVER_ID);
   mpServer->SetNotify(wxSOCKET_CONNECTION_FLAG);
   mpServer->Notify(true);

   if (! mpServer->Ok())
   {
      WriteLogMessage(_T("Server starting failure"));
      m_isRunning = false;
   }
   else
   {
        WriteLogMessage(_T("Server started"));
        m_isRunning = true;
   }

}
void FoxServer::WriteProtocol()
{
   if(s_mutexProtectingTheGlobalData->Lock()!=wxMUTEX_NO_ERROR) return;

      wxFile *logfile = new wxFile(_T("srvr_prtcl.txt"), wxFile::write);
      if(logfile != 0){
         logfile->Write(_T("List of threads - clients\n"));

         if(m_threadMutex->Lock()!=wxMUTEX_NO_ERROR) {
            logfile->Close();
            delete logfile;
            s_mutexProtectingTheGlobalData->Unlock();
            return;
         }
         for(int i=0;i<m_threads.size();i++){
            wxString tmp;
            tmp.Printf(_T("\nClient nb: %d\n Thread ID: %d\n"),  i, m_threads[i]->GetId());
            logfile->Write(tmp);
         }
         m_threadMutex->Unlock();

         logfile->Write(_T("End List of Clients\n\n"));
         logfile->Write(_T("List of Results\n"));
         for(int i=0;i<m_results.size();i++){
            wxString tmp;

            tmp.Printf(_T("\nJob nb: %d\n Job ID: %d\n Thread ID: %d\n Filename: %s\n"), i,
                         m_results[i].JobID,
                         m_results[i].threadID,
                         m_results[i].filename.c_str());
            logfile->Write(tmp);
         }
         logfile->Write(_T("End List of Results\n"));

         logfile->Close();
      }
      delete logfile;
    s_mutexProtectingTheGlobalData->Unlock();
}
void FoxServer::GetData( std::vector<GridClient> &clients, std::vector<GridResult > &results, std::vector<FoxJob > &Joblist)
{
   VFN_DEBUG_MESSAGE(__FUNCTION__,10)
   if(s_mutexProtectingTheGlobalData->Lock()!=wxMUTEX_NO_ERROR) return;

   if(m_threadMutex->Lock()!=wxMUTEX_NO_ERROR) {
      s_mutexProtectingTheGlobalData->Unlock();
      WriteLogMessage(_T("lock thread error (GetData)"));
      return;
   }
    //Clients info
    for(int i=0;i<m_threads.size();i++){
        GridClient client;
        client.name = m_threads[i]->getName();
        client.id = m_threads[i]->GetId();
        client.allCPUs = m_threads[i]->getAllCPUs();
        client.availCPUs = m_threads[i]->getAvailCPUs();
        switch(m_threads[i]->getStatus()) {
            case FG_CONNECTED:
                client.status = _T("connected");
                break;
            case FG_EXPECTING_ANSWER:
                client.status = _T("asking");
                break;
            case FG_EXPECTING_RESULT:
                client.status = _T("computing");
                break;
            default:
                client.status = _T("n/a");
                break;
        }
        clients.push_back(client);
    }
    m_threadMutex->Unlock();

    //Results. It adds only new results.
    for(int i=results.size();i<m_results.size();i++){
       results.push_back(m_results[i]);
    }
    //Jobs
    Joblist.clear();
    //save data to the outgoing list.
    for(int i=0;i<m_jobs.size();i++) {
       Joblist.push_back(m_jobs[i]);
    }

   s_mutexProtectingTheGlobalData->Unlock();
}
void FoxServer::ChangeJob(int index, FoxJob *cjob)
{
   VFN_DEBUG_MESSAGE(__FUNCTION__,10)
//change only nbRuns and nbTrial
//You can't change jobID

   if(s_mutexProtectingTheGlobalData->Lock()!=wxMUTEX_NO_ERROR) return;

   if(index >= m_jobs.size()) {
      s_mutexProtectingTheGlobalData->Unlock();
      return;
   }
   //nbsolve + nbDone <= nbRuns
   if(cjob->getNbRuns() >= (m_jobs[index].getNbDone() + m_jobs[index].getNbThread()))
      m_jobs[index].setNbRuns(cjob->getNbRuns());
   else m_jobs[index].setNbRuns(m_jobs[index].getNbDone() + m_jobs[index].getNbThread());

   m_jobs[index].setNbTrial(cjob->getNbTrial());
   m_jobs[index].setName(cjob->getName());
   m_jobs[index].setRand(cjob->randomize());

   s_mutexProtectingTheGlobalData->Unlock();
}
int FoxServer::DeleteJob(int index)
{
   VFN_DEBUG_MESSAGE(__FUNCTION__,10)
//you should call GetData(...) after this function...
/* return: -1 = can't delete - not found,
 *          0 = can't delete - found, but sent to client. Number of runs changed.
 *         1 = Job deleted
 */
   if(s_mutexProtectingTheGlobalData->Lock()!=wxMUTEX_NO_ERROR) return -1;

   if(index >= m_jobs.size()) {
      s_mutexProtectingTheGlobalData->Unlock();
      return -1;
   }
   //if job was sent to client, we can't erase it. We can only change the numer of runs...
   if((m_jobs[index].getNbDone() + m_jobs[index].getNbThread())>0) {
      m_jobs[index].setNbRuns(m_jobs[index].getNbDone() + m_jobs[index].getNbThread());
      s_mutexProtectingTheGlobalData->Unlock();
      return 0;
   }

   //todo: test this removing from vector!!!
   std::vector<FoxJob >::iterator it = m_jobs.begin()+index;
   m_jobs.erase(it);

   s_mutexProtectingTheGlobalData->Unlock();
   return 1;
}
void FoxServer::AddJobToList(FoxJob newjob)
{
   VFN_DEBUG_MESSAGE(__FUNCTION__,10)
   if(s_mutexProtectingTheGlobalData->Lock()!=wxMUTEX_NO_ERROR) return;

   m_jobs.push_back(newjob);

   s_mutexProtectingTheGlobalData->Unlock();
}
bool FoxServer::IsServerRunning()
{
   return m_isRunning;
}
void FoxServer::OnServerEvent(wxSocketEvent &event)
{
   WriteLogMessage(_T("On Server Event"));
   VFN_DEBUG_MESSAGE(__FUNCTION__,10)
   switch(event.GetSocketEvent())
   {
      case wxSOCKET_CONNECTION :
     {
         WriteLogMessage(_T("New Client"));
         if(s_mutexProtectingTheGlobalData->Lock()!=wxMUTEX_NO_ERROR) return;

         wxIPV4address tmp_addr;
         wxSocketBase *pSocket =  mpServer->Accept(false);
         pSocket->SetEventHandler(*this, GRID_SERVER_SOCKET_ID);
         pSocket->SetNotify(wxSOCKET_INPUT_FLAG | wxSOCKET_LOST_FLAG);
         pSocket->Notify(true);

         //run thread
         FoxServerThread* pThread = new FoxServerThread(   pSocket,
                                             EMPTY_MSG,
                                             s_mutexProtectingTheGlobalData,
                                             &m_results,
                                             &m_jobs,
                                             m_working_dir);

         if(pThread->Create()!=wxTHREAD_NO_ERROR) {
            s_mutexProtectingTheGlobalData->Unlock();
            return;
         }
         else pThread->Run();

         if(m_threadMutex->Lock()!=wxMUTEX_NO_ERROR) {
             s_mutexProtectingTheGlobalData->Unlock();
            return;
         }
         //add new thread to the list.
         m_threads.push_back(pThread);
         pThread->NewEvent(NEW_CONNECTION, m_threadMutex);

         m_threadMutex->Unlock();
         s_mutexProtectingTheGlobalData->Unlock();

         break;
     }
      default:
         break;
   }
}
void FoxServer::OnSocketEvent(wxSocketEvent &event)
{
   VFN_DEBUG_MESSAGE(__FUNCTION__,10)
   WriteLogMessage(_T("OnSocketEvent"));
   wxSocketBase* tmpSock = event.GetSocket();

   //do not setNotify back to LOST | INPUT in this function, it does the thread itself...
   tmpSock->SetNotify(wxSOCKET_LOST_FLAG);

   int No = -1;
   int i;

   //locking mutex for the thread list.
   if(m_threadMutex->Lock()!=wxMUTEX_NO_ERROR) {
       WriteLogMessage(_T("error: m_threadMutex locked => return"));
       return;
   }
   //who sent this socket?
   for(i=0; i<m_threads.size(); i++) {
      if(tmpSock==m_threads[i]->GetSocket()) {
         No = i;
         break;
      }
   }
   //bad client identification
   if(No==-1) {
      WriteLogMessage(_T("error: Bad client identification"));
      m_threadMutex->Unlock();
      return;
   }
   switch(event.GetSocketEvent())
   {
     case wxSOCKET_INPUT:
     {
        WriteLogMessage(_T("wxSOCKET_INPUT"));
        //call the thread.
        m_threads[i]->NewEvent(INPUT_MSG, m_threadMutex);
        m_threadMutex->Unlock();
        break;
     }
     case wxSOCKET_LOST:
     {
         wxString tmp;
         int id;
         id = m_threads[i]->GetId();
         tmp.Printf(_T("Client disconnected, thread id=%d"), id);
         WriteLogMessage(tmp);
         //call the thread "socket lost".
         m_threads[i]->NewEvent(LOST_CONNECTION, m_threadMutex);

         std::vector<FoxServerThread *>::iterator it;
         it = m_threads.begin()+i;
         m_threads.erase(it);

         //ulocking before locking global data!
         m_threadMutex->Unlock();

         if(s_mutexProtectingTheGlobalData->Lock()!=wxMUTEX_NO_ERROR) return;

         tmp.Printf(_T("removing thread %d from jobs..."), id);
         WriteLogMessage(tmp);
         //reduce nbSolving number
         for(int q=0;q<m_jobs.size();q++){
            WriteLogMessage(m_jobs[q].getListOfThreads());
            m_jobs[q].RemoveThread(id, -1);
            WriteLogMessage(m_jobs[q].getListOfThreads());
         }

         s_mutexProtectingTheGlobalData->Unlock();

         break;
     }
     default:
     {
        //just unlock the mutex
        m_threadMutex->Unlock();
        WriteLogMessage(_T("unrecognized event"));
        break;
     }
   }

   m_needUpdate=true;
   WriteLogMessage(_T("OnSocketEvent_end"));
}
void FoxServer::RunAllClients()
{
   VFN_DEBUG_MESSAGE(__FUNCTION__,10)
   if(m_threadMutex->Lock()!=wxMUTEX_NO_ERROR) return;

   for(int i=0; i<m_threads.size(); i++)
   {
      m_threads[i]->NewEvent(SEND_JOB, m_threadMutex);
   }
   m_threadMutex->Unlock();
}
