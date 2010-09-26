
#ifdef __WX__CRYST__
   #include "ObjCryst/wxCryst/wxCrystal.h"
#endif

// Fox client/server grid
#ifndef __FOX_THREAD__
#include "FoxServerThread.h"
#endif

using namespace ObjCryst;
using namespace std;

FoxServerThread::FoxServerThread(   wxSocketBase* pSocket,
                           FoxServerEvents evt,
                           wxMutex *pMutex,
                           std::vector<GridResult > *pResults,
                           std::vector<FoxJob > *pJobs) :
m_pSocket(pSocket),
m_sckt_ntf(evt),
m_tMutexObj(pMutex),
m_results(pResults),
m_jobs(pJobs)
//m_parent(parent),
//m_parentEvnt(parentEvnt)
{
   m_newEvt =         false;
   m_tThreadMutex =   new wxMutex();
   m_allCPUs = 0;
   m_availableCPUs = 0;
   m_name = _T("n/a");
   m_status = FG_N_A;
} 

FoxServerThread::~FoxServerThread()
{
   delete m_tThreadMutex;
}
bool FoxServerThread::NewEvent(FoxServerEvents evt, wxMutex *MutexToUnlock )
{
   m_sckt_ntf = evt;
   m_newEvt = true;
   return true;
}
int FoxServerThread::GetId()
{
   int i= (int) wxThread::GetId();
   return i;
}
wxSocketBase* FoxServerThread::GetSocket()
{
   return m_pSocket;
}
void *FoxServerThread::Entry()
{ 
   m_exit = false;
   do{

      if(TestDestroy()) 
         break;

      //is still connected?
      //TODO: this destroy does not sent the CLOSE event to the server
      //It has no effect in special cases...
      if(!m_pSocket->IsConnected()) {
          CloseConnection();
      }

      if(m_tThreadMutex->Lock()!=wxMUTEX_NO_ERROR) {
          WriteLogMessage(_T("m_tThreadMutex locking error (Entry)"));
          return false;
      }
      
      WriteLogMessage(_T("Looping"));

      if(m_newEvt){
          m_newEvt=false;
         switch(m_sckt_ntf)
         {
           case INPUT_MSG:
           {
              WriteLogMessage(_T("INPUT"));
              OnInput();
              WriteLogMessage(_T("INPUT end"));
              m_pSocket->SetNotify(wxSOCKET_LOST_FLAG | wxSOCKET_INPUT_FLAG);
              break;
           }
           case LOST_CONNECTION:
           {
              WriteLogMessage(_T("LOST"));
              m_pSocket->Destroy();
              m_exit=true;
              break;
           }
           case NEW_CONNECTION:
           {
              m_status = FG_CONNECTED;
              wxThread::Sleep(3000);
              SendAsk(true);
              break;
           }
           case SEND_JOB:
           {
              wxThread::Sleep(3000);
              SendAsk();
              break;
           }
           default:
             break;
         }
      }
      m_tThreadMutex->Unlock();

      wxThread::Sleep(1000);

      if(m_exit) break;
   }while(true);

   return 0;
} 
void FoxServerThread::OnInput()
{
   VFN_DEBUG_MESSAGE(__FUNCTION__,10)
   std::string in;

   //try to read input
   WriteLogMessage(_T("Reading socket"));
   if(!m_IOSocket.ReadStringFromSocket(m_pSocket, in)) {
       WriteLogMessage(m_IOSocket.getError());
       VFN_DEBUG_MESSAGE(__FUNCTION__<<":"<<in,10)
       return;
   }

   WriteLogMessage(_T("Locking mutex"));
   if(m_tMutexObj->Lock()!=wxMUTEX_NO_ERROR) {
      WriteLogMessage(_T("m_tMutexObj error"));
      return;
   }

   AnalyzeMessage(in);
   m_tMutexObj->Unlock();

   VFN_DEBUG_MESSAGE(__FUNCTION__,10)
   return; 
}
void FoxServerThread::CloseConnection()
{   
    WriteLogMessage(_T("Closing Connection"));
    m_pSocket->Destroy();
}
void FoxServerThread::OnExit() 
{
   VFN_DEBUG_MESSAGE(__FUNCTION__,10)
   WriteLogMessage(_T("This thread terminates..."));
   wxThread::Sleep(1000);
}
bool FoxServerThread::AnalyzeMessage(std::string message)
{
   VFN_DEBUG_MESSAGE(__FUNCTION__<<":"<<message,10)
   WriteLogMessage(_T("Saving message from client"));
   //wxString xtmp;
   //xtmp.Printf(_T("Server_msg_in%d.txt"), (long long) time(0));
   SaveDataAsFile(wxString::FromAscii(message.c_str()), _T("Server_msg_in.txt"));

   stringstream in_string;
   
   wxString ID, Cost;
   long nbCPUs=0;
   vector<long> ids;
   vector<double> costs;
   vector<wxString> results;
   bool newResult = false, answer=false;
   vector<long> rejectedJobs;


   in_string<<message;
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
         wxString res = getResult(wxString::FromAscii(message.c_str()), pos);
         if(res.Cmp(_T(""))!=0) {
             ids.push_back(id);
             costs.push_back(cost);
             results.push_back(res);
         }
      }
      if("answer"==tag.GetName()){
          answer = true;
          for(int i=tag.GetNbAttribute()-1;i>=0;i--){
              if(tag.GetAttributeName(i)=="nbOfAvailCPUs") {
                    wxString nb = wxString::FromAscii(tag.GetAttributeValue(i).c_str());
                    nb.ToLong((long *) &nbCPUs);
                    m_availableCPUs = nbCPUs;
              }
              if(tag.GetAttributeName(i)=="name") {
                    m_name = wxString::FromAscii(tag.GetAttributeValue(i).c_str());
              }
              if(tag.GetAttributeName(i)=="AllCPUs") {
                    wxString nb = wxString::FromAscii(tag.GetAttributeValue(i).c_str());
                    nb.ToLong((long *) &m_allCPUs);
              }
          }
      }
      if("rejectedJob"==tag.GetName()){
           long rejctd=0;
           for(int i=tag.GetNbAttribute()-1;i>=0;i--){
              if(tag.GetAttributeName(i)=="id") {
                    wxString nb = wxString::FromAscii(tag.GetAttributeValue(i).c_str());
                    nb.ToLong((long *) &rejctd);
                    rejectedJobs.push_back(rejctd);
              }
           }
      }
   }

   if(newResult){
      for(int i=0;i<results.size();i++) {
          WriteLogMessage(_T("Saving result"));
          SaveResult(results[i], ids[i], (float) costs[i]);
          //Updating JobList...
          WriteLogMessage(_T("Updating JobList"));
          for(int j=0;j<(*m_jobs).size();j++) {
             if((*m_jobs)[j].getM_ID() == (int)ids[i]) {
                wxString tmp;
                tmp.Printf(_T("result=%d (%d), m_jobs->Item(%d).getM_ID()=%d, m_jobs->Item(j)->getNbThread()=%d, done=%d"), i, results.size(), j, (*m_jobs)[j].getM_ID(), (*m_jobs)[j].getNbThread(), (*m_jobs)[j].getNbDone()); 
                WriteLogMessage(tmp);
                WriteLogMessage((*m_jobs)[j].getListOfThreads());
                tmp.Printf(_T("removing thread: %d"), GetId());
                WriteLogMessage(tmp);
                (*m_jobs)[j].RemoveThread(GetId());
                WriteLogMessage((*m_jobs)[j].getListOfThreads());
                (*m_jobs)[j].setNbDone((*m_jobs)[j].getNbDone()+1);
                tmp.Printf(_T("m_jobs->Item(j)->getNbThread()=%d, done=%d"), (*m_jobs)[j].getNbThread(), (*m_jobs)[j].getNbDone()); 
                WriteLogMessage(tmp);
             }
          } 
          //WriteLogMessage(_T("Joblist Updated"));
      }
      wxThread::Sleep(2000);
      SendAsk();
   }
   if(rejectedJobs.size()!=0) {
      WriteLogMessage(_T("rejecting jobs"));
      rejectJobs(rejectedJobs);
   } 
   if(answer) {      
      wxThread::Sleep(2000);
      SendJob(nbCPUs);
   }
   VFN_DEBUG_MESSAGE(__FUNCTION__,10)
   return true;
}
void FoxServerThread::rejectJobs(vector<long> ids)
{
    for(int i=0;i<ids.size();i++) {
        for(int q=0;q<(*m_jobs).size();q++){
            if((*m_jobs)[q].getM_ID() == (int) ids[i]) {
                (*m_jobs)[q].RemoveThread(GetId());
            }
        }
    }
}
wxString FoxServerThread::getResult(wxString message, long pos)
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
void FoxServerThread::SaveResult(wxString result, int JobID, float ResultCost)
{//this function must be under m_tMutexObj->Lock()!!
   int t = time(0);
   wxString name;
   int r = (int) ResultCost;
   #ifdef WIN32
      name.Printf(_T("GridRslt\\ID-%d_Cost-%d_Thread-%d_Time-%d.xml"), JobID, r, this->GetId(), t);
   #else
      name.Printf(_T("GridRslt/ID-%d_Cost-%d_Thread-%d_Time-%d.xml"), JobID, r, this->GetId(), t);
   #endif
   WriteLogMessage(_T("Saving result as file"));
   VFN_DEBUG_MESSAGE(__FUNCTION__<<name.ToAscii(),10)
   SaveDataAsFile(result, name);

   WriteLogMessage(_T("Saving result in glodal data"));
   GridResult newres;
   newres.ClientName = _T("--");
   newres.filename = name;
   newres.JobID = JobID;
   newres.threadID = this->GetId();
   newres.Cost = ResultCost;
   m_results->push_back(newres);
   WriteLogMessage(_T("Result saved in global data"));
}
void FoxServerThread::SendJob(int nbOfJobs)
{//this function must be under m_tMutexObj->Lock()!!
   VFN_DEBUG_MESSAGE(__FUNCTION__,10)
   vector<int> jobsToSend;

   if(nbOfJobs==0) {
       WriteLogMessage(_T("SendJob(): nbOfJobs is 0 -> no job was sent"));
       m_status = FG_CONNECTED;
       return;
   }

   WriteLogMessage(_T("Finding available jobs"));
   //find available jobs
   for(int i=0;i<(*m_jobs).size();i++){
      int available = (*m_jobs)[i].getNbRuns() - ((*m_jobs)[i].GetSolvingNb() + (*m_jobs)[i].getNbDone());
      if(available<=0) continue;
      if(available>(nbOfJobs-jobsToSend.size())) available = nbOfJobs-jobsToSend.size();
      for(int q=0;q<available;q++) {
        jobsToSend.push_back(i);
      }
      if(jobsToSend.size()==nbOfJobs) break;
   }
   if(jobsToSend.size()==0) {
      WriteLogMessage(_T("SendJob(): No job available"));
      m_status = FG_CONNECTED;
      return; //No job available
   }
   
   wxString out = _T("<FoxGrid>\n");
   int current = -1;
    for(int i=0;i<jobsToSend.size();i++) {
        int count = 0;
        current = jobsToSend[i];
        //get the number of the identical jobs => nbOfRuns
        do {
            count++;
            if(i+count>=jobsToSend.size()) break;
        }while(current == jobsToSend[i+count]);
        i = i+count-1;
        WriteLogMessage(_T("Loading job from file"));
        wxString in;
        if(!LoadFile((*m_jobs)[jobsToSend[i]].getFileName(), in)) {
            WriteLogMessage(_T("Can't load file") + (*m_jobs)[jobsToSend[i]].getFileName() );
            m_status = FG_CONNECTED;
            return;
        }
        wxString header;
        header.Printf(_T("<ClientJob ID=\"%d\" nbTrials=\"%d\" nbOfRuns=\"%d\" rand=\"%d\">\n"), (*m_jobs)[jobsToSend[i]].getM_ID(), (*m_jobs)[jobsToSend[i]].getNbTrial(), count, (int) (*m_jobs)[jobsToSend[i]].randomize());
        out += header;
        out += in;
        out += _T("\n</ClientJob>\n");
        wxString tmp;
        tmp.Printf(_T("AddThread info: m_jobs[%d], GetID()=%d, count=%d"), jobsToSend[i], GetId(), count);
        WriteLogMessage(tmp);
        (*m_jobs)[jobsToSend[i]].AddThread(GetId(), count);
    }
    out += _T("</FoxGrid>\n");
    WriteLogMessage(_T("Sending Job"));
    if(!m_IOSocket.WriteStringToSocket(m_pSocket, string(out.ToAscii()))) {
        WriteLogMessage(m_IOSocket.getError());
        m_status = FG_CONNECTED;
        //remove thread from jobs
        for(int i=0;i<jobsToSend.size();i++) {
            (*m_jobs)[jobsToSend[i]].RemoveThread(GetId());
        }
        return;
    }  
    WriteLogMessage(_T("Job sent"));
    m_status = FG_EXPECTING_RESULT;
}
bool FoxServerThread::SendAsk(bool getClientInfo)
{
    m_pSocket->SetNotify(wxSOCKET_LOST_FLAG);
    VFN_DEBUG_MESSAGE(__FUNCTION__,10)
    wxString out;
    out = _T("<FoxGrid>\n <ask nbOfAvailCPUs=\"?\" ");
    if(getClientInfo) {
        out+= _T("info=\"?\" ");
    }
    out+= _T("/>\n</FoxGrid>\n");
    VFN_DEBUG_MESSAGE(__FUNCTION__<<":"<<out,10)
    
    WriteLogMessage(_T("Sending ask..."));
    if(!m_IOSocket.WriteStringToSocket(m_pSocket, string(out.ToAscii()))) {
        WriteLogMessage(m_IOSocket.getError());
        return false;
    }
    m_status = FG_EXPECTING_ANSWER;
    WriteLogMessage(_T("ask sent"));
    m_pSocket->SetNotify(wxSOCKET_LOST_FLAG | wxSOCKET_INPUT_FLAG);
    return true;
}
void FoxServerThread::SaveDataAsFile(wxString out, wxString filename)
{
   wxFile outFile(filename, wxFile::write);
   if(outFile.IsOpened())
   {
      outFile.Write(out);
      outFile.Close();
   }
}   
bool FoxServerThread::LoadFile(wxString filename, wxString &in)
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
void FoxServerThread::WriteLogMessage(wxString msg)
{
#if __SERVER_LOGS
   wxString filename;
   filename.Printf(_T("thread_log_%d.txt"), GetId());
   wxFile logfile(filename, wxFile::write_append);
   if(logfile.IsOpened())
   {
      wxDateTime datetime = wxDateTime::Now();
      logfile.Write(datetime.Format(_T("%X ")) + msg + _T("\n"));
      logfile.Close();
   }
#endif
}
wxString FoxServerThread::getName()
{
    return m_name;
}
long FoxServerThread::getAvailCPUs()
{
    return m_availableCPUs;
}
long FoxServerThread::getAllCPUs()
{
    return m_allCPUs;
}
ServerThreadStatus FoxServerThread::getStatus()
{
    return m_status;
}
