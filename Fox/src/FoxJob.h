//#include <wx/string.h>
//#include <wx/thread.h>

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
#endif

#include <vector>

class FoxJob
{
public:
   FoxJob(void);
   ~FoxJob(void);

   void AddThread(int threadID, int nbCPUs);

   //remove thread's id from thread list (default: remove only one)
   //nbThreads=-1 remove all thread's ids in the list (called when lost_socket occurs)
   void RemoveThread(int threadID, int nbThreads = 1);

   wxString getListOfThreads();
   void setName(wxString name);
   void setM_ID(int id);
   void setNbTrial(long nbTrial);
   void setNbRuns(int nbRuns);
   void setNbDone(int nbDone);
   void setStatus(short status);
//   void setNbThread(int nbThread);
   void replaceThreadID(std::vector<int> threadID);
   void setFileName(wxString name);
   bool randomize();

   void       setRand(bool rand);
   int        getM_ID();
   wxString   getName();
   long       getNbTrial();
   int        getNbRuns();
   int        getNbDone();
   short      getStatus();
   int        getNbThread();
   int        GetSolvingNb();
   std::vector<int> getThreadID();
   wxString getFileName();

private:
   wxString    m_name;
   wxString    m_fileName;
   int         m_ID;
   long        m_nbTrial;
   int         m_nbRuns;
   int         m_nbDone;
   short       m_status;
   bool        m_randomize;
   std::vector<int> m_ThreadID;
};
