#include "FoxJob.h"

FoxJob::FoxJob(void)
{
}
FoxJob::~FoxJob(void)
{
}
void FoxJob::AddThread(int threadID, int nbCPUs)
{  
    for(int i=0;i<nbCPUs;i++){
       m_ThreadID.push_back(threadID);
    }
}
void FoxJob::RemoveThread(int threadID, int nbThreads)
{
    if(nbThreads==0) return;
    std::vector<int>::iterator it;
    for(int i=0;i<m_ThreadID.size();i++){
        if(m_ThreadID[i]==threadID) {
            it=m_ThreadID.begin()+i;
            m_ThreadID.erase(it);
            i--;
            nbThreads--;
            if(nbThreads==0) break;
        }
    }
}
wxString FoxJob::getListOfThreads()
{
    wxString out=_T("Thread List:\n");
    for(int i=0;i<m_ThreadID.size();i++){
        wxString tmp;
        tmp.Printf(_T("[%d]: %d\n"), i, m_ThreadID[i]);
        out+=tmp;
    }
    return out;
}
void FoxJob::replaceThreadID(std::vector<int> threadIDs) {
    m_ThreadID.clear();
    for(int i=0;i<threadIDs.size();i++){
        m_ThreadID.push_back(threadIDs[i]);
    }
}
std::vector<int> FoxJob::getThreadID() {
    return m_ThreadID;
}
int FoxJob::GetSolvingNb() {
   return this->getNbThread();
}
void FoxJob::setName(wxString name) {
    this->m_name = name;
}
void FoxJob::setM_ID(int id){
    this->m_ID = id;
}
void FoxJob::setNbTrial(long nbTrial){
    this->m_nbTrial = nbTrial;
}
void FoxJob::setNbRuns(int nbRuns){
    this->m_nbRuns = nbRuns;
}
void FoxJob::setNbDone(int nbDone){
    this->m_nbDone = nbDone;
}
void FoxJob::setStatus(short status){
    this->m_status = status;
}
void FoxJob::setFileName(wxString name){
    this->m_fileName = name;
}

int FoxJob::getM_ID(){
    return m_ID;
}
wxString FoxJob::getName(){
    return m_name;
}
long FoxJob::getNbTrial(){
    return m_nbTrial;
}
int FoxJob::getNbRuns(){
    return m_nbRuns;
}
int FoxJob::getNbDone(){
    return m_nbDone;
}
short FoxJob::getStatus(){
    return m_status;
}
int FoxJob::getNbThread(){
    return m_ThreadID.size();
}
wxString FoxJob::getFileName(){
    return m_fileName;
}
bool FoxJob::randomize() {
    return m_randomize;
}
void FoxJob::setRand(bool rand) {
    m_randomize = rand;
}