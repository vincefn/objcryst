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
#endif

//#include <locale.h>
//#include <fstream>
//#include <sstream>
//#include <list>
#include <vector>

#define __FOX_SERVER__ 

#include "FoxServerThread.h"
//#include "wxCryst/wxCryst.h"

//#include <wx/frame.h>


class FoxServer: public wxFrame
{
   public:
     FoxServer();
     ~FoxServer();
     void StartGridServer();
     void OnServerEvent(wxSocketEvent &event);
     void OnSocketEvent(wxSocketEvent &event);
     void WriteProtocol();
     void GetData(std::vector<GridClient> &Clients, std::vector<GridResult > &results,  std::vector<FoxJob > &Joblist);
     bool IsServerRunning();
     void RunAllClients();
     void AddJobToList(FoxJob newjob);
     void ChangeJob(int index, FoxJob *cjob);
     int DeleteJob(int index);

   protected:

     void WriteLogMessage(wxString msg);

     wxMutex            * s_mutexProtectingTheGlobalData;
     wxMutex            * mutexMessageLog;
     wxMutex            * m_threadMutex;
     wxSocketServer     * mpServer;
     std::vector<FoxServerThread * > m_threads;
     std::vector<GridResult >      m_results;
     bool                 m_isRunning;
     bool                 m_needUpdate;
     std::vector<FoxJob > m_jobs;
     DECLARE_EVENT_TABLE()
};
