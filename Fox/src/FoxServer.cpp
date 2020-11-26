#include "FoxServer.h"
#include <wx/listimpl.cpp>

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
//static const long SERVER_THREAD_EVENT=                   WXCRYST_ID();


#define EVT_WORKER(func) DECLARE_EVENT_TABLE_ENTRY( wxEVT_WORKER, -1, -1, (wxObjectEventFunction) (wxEventFunction) (WorkerEventFunction) & func, (wxObject *) NULL ),

typedef void (wxEvtHandler::*WorkerEventFunction)(WorkerEvent&);


//DEFINE_EVENT_TYPE(SERVER_THREAD_EVENT)


BEGIN_EVENT_TABLE(FoxServer, wxFrame)
   EVT_SOCKET(GRID_SERVER_ID,                       FoxServer::OnServerEvent)
   //EVT_COMMAND(wxID_ANY, SERVER_THREAD_EVENT,                 FoxServer::OnThreadEvent)
   EVT_WORKER(FoxServer::OnWorkerEvent)
   //EVT_SOCKET(GRID_SERVER_SOCKET_ID,                FoxServer::OnSocketEvent)
END_EVENT_TABLE()

//WX_DEFINE_LIST(TList);

ThreadWorker::ThreadWorker(wxSocketBase         *pSocket, 
                           wxString              workingDir, 
                           wxFrame              *parent,
                           vector<GridResult >  *results,
                           vector<FoxJob >      *jobs,
                           wxMutex              *mutexProtecting_Jobs_Results)
{
    m_socket = pSocket;
    //Notify() cannot be called in thread context. We have to detach from main loop
    //before switching thread contexts.
    m_socket->Notify(false);
    m_socket->SetFlags(wxSOCKET_WAITALL|wxSOCKET_BLOCK);
    pSocket->GetPeer(m_peer);
    m_working_directory = workingDir;
    m_parent = parent;
    m_global_results = results;
    m_global_jobs = jobs;
    m_mutexProtecting_Jobs_Results = mutexProtecting_Jobs_Results;
}
bool ThreadWorker::LoadFile(wxString filename, wxString &in)
{
   wxFile infile(filename, wxFile::read);
   in.Clear();
   if(!infile.IsOpened()) return false;
   long len = infile.Length();
   char *buffer;
   buffer = (char*) calloc(len+1, sizeof(char));
   infile.Read(buffer, len);
   infile.Close();
   in = wxString::FromAscii(buffer);
   free(buffer);
   return true;
}
/*
bool ThreadWorker::write_socket(wxSocketBase* m_socket, wxCharBuffer message)
{
      //deprechated

      m_socket->SetFlags(wxSOCKET_WAITALL);

      // Note that len is in kbytes here!
      unsigned int len  = message.length();
      WriteLogMessage(wxString::Format("ThreadWorker: Sending header of the message of %d kilobytes", len));
      m_socket->Write(&len, sizeof(int));
      if (m_socket->Error()) {
         WriteLogMessage("ThreadWorker: Write error");
         return false;
      }
      WriteLogMessage("ThreadWorker: Sending the message ...");
      m_socket->Write(message, len);
      WriteLogMessage(m_socket->Error() ? _("ThreadWorker: failed !\n") : _("ThreadWorker: done\n"));
      return true;
}
*/
int ThreadWorker::GetId()
{
   int i= (int) wxThread::GetId();
   return i;
}
/*
bool ThreadWorker::read_socket(wxSocketBase* m_socket, wxCharBuffer &message)
{
    //deprechated

    WriteLogMessage("ThreadWorker: Reading the socket ...");
    m_socket->SetFlags(wxSOCKET_WAITALL);
    unsigned int len;
    m_socket->Read(&len, sizeof(int));
    if (m_socket->Error()) {
        WriteLogMessage("ThreadWorker: Read error");
        //wxGetApp().AddPendingEvent(e);
        return false;
    }
    int processed = m_socket->LastCount();
    WriteLogMessage(wxString::Format("ThreadWorker: %d bytes read in the header", processed));           

    if (len == 0)
    {
        //e.m_exit = true;
        //return 0;
        WriteLogMessage("ThreadWorker: 0 bytes in socket, nothing to read...");
        return true;
    }

    
    message.extend(len);
    WriteLogMessage(wxString::Format("Message header was: size = %d (bytes)",len));
    WriteLogMessage("ThreadWorker: Reading message ...");
    m_socket->Read(message.data(), len);

    if (m_socket->Error())
    {
        WriteLogMessage("ThreadWorker: Read error");
        //wxGetApp().AddPendingEvent(e);
        return 0;
    }
    processed = m_socket->LastCount();
    WriteLogMessage(wxString::Format("ThreadWorker: %d bytes readed", processed));

    return true;
           
}
*/
void ThreadWorker::WriteLogMessage(wxString msg)
{
#if __SERVER_LOGS
   wxString filename;
   filename.Printf(_T("thread_log_%d.txt"), GetId());
#ifdef WIN32
   filename = m_working_directory + _T("\\") + filename;
#else
   filename = m_working_directory + _T("/") + filename;
#endif

   wxFile logfile(filename, wxFile::write_append);
   if(logfile.IsOpened())
   {
      wxDateTime datetime = wxDateTime::Now();
      logfile.Write(datetime.Format(_T("%X ")) + msg + _T("\n"));
      logfile.Close();
   }
#endif
}
void ThreadWorker::cannot_calculate_this_job(FoxJob fj)
{
    if(m_mutexProtecting_Jobs_Results->Lock()!=wxMUTEX_NO_ERROR) return;
    for(int i=0;i<(*m_global_jobs).size();i++) {
        if((*m_global_jobs)[i].getM_ID() == fj.getM_ID()) {
            (*m_global_jobs)[i].RemoveThread(GetId(), fj.getNbRuns());
            break;
        }
    }
    m_mutexProtecting_Jobs_Results->Unlock();
}
void ThreadWorker::RemoveThreadFromJobList()
{//run this before closing this thread
    if(m_mutexProtecting_Jobs_Results->Lock()!=wxMUTEX_NO_ERROR) return;
    for(int j=0;j<(*m_global_jobs).size();j++) {
       (*m_global_jobs)[j].RemoveThread(GetId(), -1);
    }
    m_mutexProtecting_Jobs_Results->Unlock();
}
void ThreadWorker::AddResultToJobList(int id)
{
    if(m_mutexProtecting_Jobs_Results->Lock()!=wxMUTEX_NO_ERROR) return;
    for(int j=0;j<(*m_global_jobs).size();j++) {
        if((*m_global_jobs)[j].getM_ID() == id) {
            (*m_global_jobs)[j].RemoveThread(GetId());
            (*m_global_jobs)[j].setNbDone((*m_global_jobs)[j].getNbDone()+1);
        }
    }
    m_mutexProtecting_Jobs_Results->Unlock();
}
vector<FoxJob > ThreadWorker::getJobsToCalculate(int n)
{//returns copies of found available jobs. 
 //Use their changed getNbRuns() to see how many times you have to run it

    vector<FoxJob > res;
    //m_global_jobs is shared variable!
    if(m_mutexProtecting_Jobs_Results->Lock()!=wxMUTEX_NO_ERROR) return res;
    int all=0;
    //get available jobs and save it to the res
    for(int i=0;i<(*m_global_jobs).size();i++) {
        int nb = (*m_global_jobs)[i].getAvailableRuns();
        if(nb<=0) continue;
        //do not save more than 'n'
        if(nb>(n-all)) nb = n-all;
        all += nb;
        FoxJob tmp = (*m_global_jobs)[i];

        //just change the numbers of done and runs in the copy...
        tmp.setNbDone(0);
        tmp.setNbRuns(nb);       
        res.push_back(tmp);

        //we have to say, that this job is already reserved
        //if sending is not successfull, we have to return it back!
        (*m_global_jobs)[i].AddThread(GetId(), nb);

        if(all>=n) break;
    }
    m_mutexProtecting_Jobs_Results->Unlock();
    return res;
}
bool ThreadWorker::analyze_message_and_get_answer(wxString msg, wxString &answer)
{
    /*
    out = _T("<FoxGrid>\n <currentstate ");
    out << _T("name=\"") << getMyHostname() << _T("\" ");
    out << _T("freeCPUs=\"") << getNbOfUnusedProcesses() << _T("\" ");
    out << _T("availableCPUs=\"") << getNbOfAvailCPUs() << _T("\" ");
    out += _T(" />\n</FoxGrid>\n");
    */
    stringstream in_string;

   wxString ID, Cost, name;
   int availableCPUs=-1, freeCPUs=-1;
   vector<long> ids;
   vector<double> costs;
   vector<wxString> results;
   bool newResult = false;
   //vector<long> rejectedJobs;

   in_string<<msg;
   while(true)
   {
      XMLCrystTag tag;
      in_string>>tag;
      if(true==in_string.eof()) break;
      
      if( ("result"==tag.GetName()) && (!tag.IsEndTag()) ){
         newResult = true;
         long id;
         double cost;
         for(int i=tag.GetNbAttribute()-1;i>=0;i--){
            if(tag.GetAttributeName(i)=="ID"){
               ID = wxString::FromAscii(tag.GetAttributeValue(i).c_str());
               ID.ToLong((long *) &id);
            }
            if(tag.GetAttributeName(i)=="Cost"){
               Cost = wxString::FromAscii(tag.GetAttributeValue(i).c_str());
               Cost.ToDouble(&cost);
            }
         }
         long pos = in_string.tellg();
         wxString res = getResult(msg, pos);
         if(res.Cmp(_T(""))!=0) {
             ids.push_back(id);
             costs.push_back(cost);
             results.push_back(res);
         }
      }
      if("currentstate"==tag.GetName()){
          for(int i=tag.GetNbAttribute()-1;i>=0;i--){
              if(tag.GetAttributeName(i)=="availableCPUs") {
                    wxString nb = wxString::FromAscii(tag.GetAttributeValue(i).c_str());
                    nb.ToLong((long *) &availableCPUs);
              }
              if(tag.GetAttributeName(i)=="name") {
                    name = wxString::FromAscii(tag.GetAttributeValue(i).c_str());
              }
              if(tag.GetAttributeName(i)=="freeCPUs") {
                    wxString nb = wxString::FromAscii(tag.GetAttributeValue(i).c_str());
                    nb.ToLong((long *) &freeCPUs);
              }
          }
      }
   }

   if(availableCPUs!=-1) {
       m_availableCPUs = availableCPUs;
   }
   if(name.Length()!=0) {
       m_name = name;
   }
   if(freeCPUs!=-1) {
       m_freeCPUs = freeCPUs;
   }

   if(newResult) {
      for(int i=0;i<results.size();i++) {
          SaveResult(results[i], ids[i], (float) costs[i]);
          //Updating JobList...
          AddResultToJobList(ids[i]);
      }
   }

   if(m_freeCPUs>0) {
       WriteLogMessage("Looking for a new job to send ... ");
       //send some job ...
       //getting copies of global jobs - > we will not need mutex for them.
       vector<FoxJob> jobs = getJobsToCalculate(m_freeCPUs);
       if(jobs.size()==0) {
           WriteLogMessage("Nothing found!");
           return true;
       }

       answer += _T("<FoxGrid>\n");
       for(int i=0;i<jobs.size();i++) {
           wxString msg;
           if(!LoadFile(jobs[i].getFileName(), msg)) {
               //this should not happen...
               cannot_calculate_this_job(jobs[i]);
               continue;
           }
           m_freeCPUs -= jobs[i].getNbRuns();
           wxString header;
           header.Printf(_T("<ClientJob ID=\"%d\" nbTrials=\"%d\" nbOfRuns=\"%d\" rand=\"%d\">\n"), jobs[i].getM_ID(), jobs[i].getNbTrial(), jobs[i].getNbRuns(), (int) jobs[i].randomize());
           answer += header;
           answer += msg;
           answer += _T("\n</ClientJob>\n");
           WriteLogMessage("Job found: " + header);
       }
       answer += _T("</FoxGrid>\n");
   }
   return true;
}
void ThreadWorker::SaveResult(wxString result, int JobID, float ResultCost)
{
   int t = time(0);
   wxString name;
   int r = (int) ResultCost;
   #ifdef WIN32
      name.Printf(_T("GridRslt\\ID-%d_Cost-%d_Thread-%d_Time-%d.xml"), JobID, r, this->GetId(), t);
      name = m_working_directory + _T("\\") + name;
   #else
      name.Printf(_T("GridRslt/ID-%d_Cost-%d_Thread-%d_Time-%d.xml"), JobID, r, this->GetId(), t);
      name = m_working_directory + _T("/") + name;
   #endif
   WriteLogMessage(_T("Saving result as file"));
   VFN_DEBUG_MESSAGE(__FUNCTION__<<name.ToAscii(),10)
   SaveDataAsFile(result, name);

   if(m_mutexProtecting_Jobs_Results->Lock()!=wxMUTEX_NO_ERROR) return;

   WriteLogMessage(_T("Saving result in glodal data"));
   GridResult newres;
   newres.ClientName = _T("--");
   newres.filename = name;
   newres.JobID = JobID;
   newres.threadID = this->GetId();
   newres.Cost = ResultCost;
   m_global_results->push_back(newres);
   WriteLogMessage(_T("Result saved in global data"));

   m_mutexProtecting_Jobs_Results->Unlock();
}
wxString ThreadWorker::getResult(wxString message, long pos)
{
    wxString result;
    //copy of inmsg
    wxString in = message;
    in.Remove(0,pos);
    //find end of the job
    int p = in.First(_T("</result>"));
    if(p==-1) return _T("");
    //save it
    result = in.Left(p);
   VFN_DEBUG_MESSAGE(__FUNCTION__<<":"<<message<<","<<in,10)
    return result;
}
void ThreadWorker::SaveDataAsFile(wxString out, wxString filename)
{
   wxFile outFile(filename, wxFile::write);
   if(outFile.IsOpened())
   {
      outFile.Write(out);
      outFile.Close();
   }
}
wxThread::ExitCode ThreadWorker::Entry()
{
    //WorkerEvent e(this);
    if (!m_socket->IsConnected())
    {
        //LogWorker("ThreadWorker: not connected",wxLOG_Error);
        return 0;
    }
    //int to_process = -1;
    while (m_socket->IsConnected())
    {
        WriteLogMessage("Waiting for reading 5s");
        if(m_socket->WaitForRead(5, 0)) {
            if(!m_socket->IsConnected()) {
                WriteLogMessage("Connection lost");
                break;
            }
            wxCharBuffer ch_buf;
            wxString msg;
            WriteLogMessage("Reading message from the socket");
            if(!m_IO.ReadStringFromSocket(m_socket, msg)) {
                WriteLogMessage("ERROR: reading from the socket");
                continue;
            }

            
            wxString tmp_path;
            #ifdef WIN32
            tmp_path = m_working_directory + _T("\\server_msg_in.txt");
            #else
            tmp_path = m_working_directory + _T("/server_msg_in.txt");
            #endif
            SaveDataAsFile(msg, tmp_path);

            wxString answ;
            if(!analyze_message_and_get_answer(msg, answ)) {
                WriteLogMessage("ERROR: something wrong during analyzing the message...");
            }

            if(answ.length()==0) {
                //send just dummy answer, because something has to be send (client is waiting for it...)
                answ = "Nothing to do...";
            }
            WriteLogMessage("Sending answer ...");
            if(!m_IO.WriteStringToSocket(m_socket, answ)) {
                WriteLogMessage("ERROR: sending answer returned false");
                continue;
            }

            //inform parent thread about news. Sending event is a save way how to do it.
            WorkerEvent e(this);
            e.m_exit = false;
            e.m_thread_info.m_name = m_name;
            e.m_thread_info.m_freeCPUs = m_freeCPUs;
            e.m_thread_info.m_availableCPUs = m_availableCPUs;
            e.m_thread_info.m_id = GetId();
            m_parent->GetEventHandler()->AddPendingEvent( e );
        }
    }

    WriteLogMessage("Connection lost ...");

    //Remove this thread ID from all jobs before leaving
    RemoveThreadFromJobList();

    //inform parent thread about leaving
    WorkerEvent e(this);
    e.m_exit = true;
    e.m_thread_info.m_name = m_name;
    e.m_thread_info.m_freeCPUs = m_freeCPUs;
    e.m_thread_info.m_availableCPUs = m_availableCPUs;
    e.m_thread_info.m_id = GetId();
    m_parent->GetEventHandler()->AddPendingEvent( e );

    m_socket->Destroy();
    return 0;
}


FoxServer::FoxServer():
mpServer(0)
{
   s_mutexProtecting_Jobs_Results = new wxMutex();
   //m_needUpdate = false;
   m_isRunning = false;
   srand( (unsigned)time( NULL ) );
   m_working_dir = _T("");
}
FoxServer::~FoxServer()
{
   //todo: clear m_results
   delete s_mutexProtecting_Jobs_Results;
}
void FoxServer::WriteLogMessage(wxString msg)
{
#if __SERVER_LOGS
   wxString filename;
#ifdef WIN32
   filename = GetWorkingDir() + _T("\\server.log");
#else
   filename = GetWorkingDir() + _T("/server.log");
#endif
   wxFile logfile(filename, wxFile::write_append);
   if(logfile.IsOpened())
   {
      wxDateTime datetime = wxDateTime::Now();
      logfile.Write(datetime.Format(_T("%X ")) + msg + _T("\n"));
      logfile.Close();
   }
#endif
   (*fpObjCrystInformUser)(msg.ToStdString());
}
void FoxServer::OnWorkerEvent(WorkerEvent& pEvent)
{
     if(pEvent.m_exit == true) {
         WriteLogMessage(_T("Thread closed [") + wxString::Format("%d", pEvent.m_thread_info.m_id) + _T(", ") + pEvent.m_thread_info.m_name + _T("]"));
         std::vector<MY_STHREAD>::iterator it;
         for(int i=0;i<m_threads_info.size();i++) {
             it = m_threads_info.begin() + i;
             if(pEvent.m_thread_info.m_id == m_threads_info[i].m_id) {
                 m_threads_info.erase(it);
                 break;
             }
         }
     } else {
         //updating basic thread info
         for(int i=0;i<m_threads_info.size();i++) {
             if(pEvent.m_thread_info.m_id == m_threads_info[i].m_id) {
                 m_threads_info[i] = pEvent.m_thread_info;
             }
         }
     }
}
/*
void FoxServer::OnThreadEvent(wxCommandEvent &event)
{
    WriteLogMessage("Nova udalost z Threadu");
    WriteLogMessage(event.GetString());
}
*/
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
/*
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
 */
void FoxServer::GetData( std::vector<GridClient> &clients, std::vector<GridResult > &results, std::vector<FoxJob > &Joblist)
{
    VFN_DEBUG_MESSAGE(__FUNCTION__,10)
    //Clients info
    for(int i=0;i<m_threads_info.size();i++){
        GridClient client;
        client.name = m_threads_info[i].m_name;
        client.id = m_threads_info[i].m_id;
        client.allCPUs = m_threads_info[i].m_availableCPUs;
        client.availCPUs = m_threads_info[i].m_freeCPUs;

        /*
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
        */
        clients.push_back(client);
    }

    if(s_mutexProtecting_Jobs_Results->Lock()!=wxMUTEX_NO_ERROR) return;

    //Results. It adds only new results.
    for(int i=results.size();i<m_results.size();i++){
       results.push_back(m_results[i]);
    }
    
    Joblist.clear();
    //save data to the outgoing list.
    for(int i=0;i<m_jobs.size();i++) {
       Joblist.push_back(m_jobs[i]);
    }
    s_mutexProtecting_Jobs_Results->Unlock();
  
}
void FoxServer::UpdateJob(int index, FoxJob *cjob)
{
    
   VFN_DEBUG_MESSAGE(__FUNCTION__,10)
//change only nbRuns and nbTrial
//You can't change jobID

   if(s_mutexProtecting_Jobs_Results->Lock()!=wxMUTEX_NO_ERROR) return;

   if(index >= m_jobs.size()) {
      s_mutexProtecting_Jobs_Results->Unlock();
      return;
   }
   //nbsolve + nbDone <= nbRuns
   if(cjob->getNbRuns() >= (m_jobs[index].getNbDone() + m_jobs[index].getNbThread()))
      m_jobs[index].setNbRuns(cjob->getNbRuns());
   else m_jobs[index].setNbRuns(m_jobs[index].getNbDone() + m_jobs[index].getNbThread());

   m_jobs[index].setNbTrial(cjob->getNbTrial());
   m_jobs[index].setName(cjob->getName());
   m_jobs[index].setRand(cjob->randomize());

   s_mutexProtecting_Jobs_Results->Unlock();
   
}
int FoxServer::DeleteJob(int index)
{
   VFN_DEBUG_MESSAGE(__FUNCTION__,10)
//you should call GetData(...) after this function...
/* return: -1 = can't delete - not found,
 *          0 = can't delete - found, but number of runs changed.
 *          1 = Job deleted
 */
    
   if(s_mutexProtecting_Jobs_Results->Lock()!=wxMUTEX_NO_ERROR) return -1;

   if(index >= m_jobs.size()) {
      s_mutexProtecting_Jobs_Results->Unlock();
      return -1;
   }
   //if job was sent to client, we can't erase it. We can only change the numer of runs...
   if((m_jobs[index].getNbDone() + m_jobs[index].getNbThread())>0) {
      m_jobs[index].setNbRuns(m_jobs[index].getNbDone() + m_jobs[index].getNbThread());
      s_mutexProtecting_Jobs_Results->Unlock();
      return 0;
   }

   //remove it
   std::vector<FoxJob >::iterator it = m_jobs.begin()+index;
   m_jobs.erase(it);

   s_mutexProtecting_Jobs_Results->Unlock();
   return 1;
}
void FoxServer::AddJobToList(FoxJob newjob)
{
   VFN_DEBUG_MESSAGE(__FUNCTION__,10)
   if(s_mutexProtecting_Jobs_Results->Lock()!=wxMUTEX_NO_ERROR) return;
   m_jobs.push_back(newjob);
   s_mutexProtecting_Jobs_Results->Unlock();
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
            wxSocketBase* sock = mpServer->Accept();
            wxIPV4address addr;
            if (!sock->GetPeer(addr))
            {
                WriteLogMessage("Server: cannot get peer info");
            } else {
                WriteLogMessage(wxString::Format("Got connection from %s:%d", addr.IPAddress().c_str(), addr.Service()));
            }
            bool createThread = true;

            if (createThread)
            {
                WriteLogMessage("Server: Creating Thread");
                ThreadWorker* c = new ThreadWorker(sock,
                                                   m_working_dir, 
                                                   this,
                                                   &m_results,
                                                   &m_jobs,
                                                   s_mutexProtecting_Jobs_Results);
                if (c->Create() == wxTHREAD_NO_ERROR)
                {
                    MY_STHREAD tmp;
                    tmp.m_availableCPUs = 0;
                    tmp.m_freeCPUs = 0;
                    tmp.m_id = c->GetId();
                    //saving basic info, that will be used for communication with GUI (without any need of mutex)
                    m_threads_info.push_back(tmp);
                    //here we save pointers to threads
                    m_threads.push_back(c);
                    //m_threadWorkers.Append(c);
                    /*
                    if (m_threadWorkers.GetCount() > m_maxThreadWorkers)
                    m_maxThreadWorkers++;
                    m_threadWorkersCreated++;
                    */
                    WriteLogMessage("Server: Running Thread");
                    c->Run();
                }
                else
                {
                    WriteLogMessage("Server: cannot create next thread ");
                };
            }
     }
      default:
         break;
   }
}
/*
void FoxServer::OnSocketEvent(wxSocketEvent &event)
{
   VFN_DEBUG_MESSAGE(__FUNCTION__,10)
   WriteLogMessage(_T("OnSocketEvent - this should not happen!!!"));
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
*/

