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
   #include "wx/socket.h"
   #include "wx/grid.h"
   #include "wx/dynarray.h"
   #include "wx/stdpaths.h"
#endif

#include "ObjCryst/ObjCryst/IO.h"
#include "ObjCryst/ObjCryst/Crystal.h"
#include "ObjCryst/ObjCryst/PowderPattern.h"
#include "ObjCryst/ObjCryst/DiffractionDataSingleCrystal.h"
#include "ObjCryst/RefinableObj/GlobalOptimObj.h"
#include "FoxServer.h"


#define __FOX_SERVER__

class WXFoxServer : public wxWindow
{
public:
   WXFoxServer(wxWindow* parent, wxString workingDir);
   ~WXFoxServer(void);
   void Clear();
   
   bool            m_dataLoaded;

private:
   void InitServer();
   void OnGridResultClick(wxGridEvent &event);
   void OnGridJobClick(wxGridEvent &event);
   void OnShowResults(wxCommandEvent& event);
   void RunLocalClient(wxCommandEvent& event); 
   void RunALLClient(wxCommandEvent& event); 
   void UpdateLists(wxTimerEvent& event);
   void OnNewJob(wxCommandEvent& event);
   void OnEditJob(wxCommandEvent& event);
   void OnDeleteJob(wxCommandEvent& event);
   void OnLoadJob(wxCommandEvent& event);
   int  GenerateJobID();
   void UpdateJobList();
   void UpdateResultList();
   bool ShowSetJobWindow(int ID, wxString &name, long &trials, long &runs, bool &rand);
   bool isFileFoxJob(wxString path, wxString &name, int &id, long &nbOfTrial, long &nbRun, bool &rand);
   void saveJobHeader(wxString filename, int ID, wxString name, long nbOfTrial, int nbRun, bool rand);
   void ChangeJobHeader(wxString filename, int ID, wxString name, long nbOfTrial, int nbRun, bool rand);
   void AddJob(wxString filename, wxString name, int id, long nbOfTrial, long nbRun, bool rand);
   void SaveDataAsFile(wxString out, wxString filename);
   bool LoadFile(wxString filename, wxString &in);
   bool isJobLoaded(long ID);

   wxWindow         * m_parent;
   wxGrid            * m_ResultTable;
   wxGrid            * m_ClientTable;
   wxGrid            * m_JobListTable;
   //wxTextCtrl         * m_EventsWindow;
   //wxTextCtrl         * m_ClientWindow;
   //wxTextCtrl         * m_JobListWindow;
   wxTimer            * m_UpdateTimer;
   FoxServer         * m_FoxServer;
   std::vector<FoxJob >         m_jobs;
   std::vector<GridResult >     m_results;
   wxMutex            * m_WXFoxServerDataMutex;
   wxString             m_working_dir;

   DECLARE_EVENT_TABLE()
};
