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
#include "FoxGridMaster.h"


#define __FOX_SERVER__

class WXFoxMaster : public wxWindow
{
public:
   WXFoxMaster(wxWindow* parent, wxString workingDir);
   ~WXFoxMaster(void);
   void Clear();

   bool            m_dataLoaded;

private:
   void InitServer();
   void OnGridResultClick(wxGridEvent &event);
   void OnGridJobClick(wxGridEvent &event);
   void OnShowResults(wxCommandEvent& event);
   void OnShowResultsServer(wxCommandEvent& event);
   void RunLocalClient(wxCommandEvent& event);
   //void RunALLClient(wxCommandEvent& event);
   void UpdateLists(wxTimerEvent& event);
   void OnNewJob(wxCommandEvent& event);
   void OnEditJob(wxCommandEvent& event);
   void EditJob();
   void OnDeleteJob(wxCommandEvent& event);
   void OnLoadJob(wxCommandEvent& event);
   void UpdateJobList();
   void UpdateResultList();
   void UpdateClientList(vector<FoxGridMaster::SLAVE_FOX_INFO> SI);
   bool ShowEditJobWindow(long ID, wxString &name, int &trials, int &runs, bool &rand, bool onlyRuns);
   bool isFileFoxJob(wxString path, wxString &name, long &id, int &nbOfTrial, int &nbRun, bool &rand);
   void saveJobHeader(wxString filename, long ID, wxString name, int nbOfTrial, int nbRun, bool rand, wxString &data);
   void ChangeJobHeader(wxString filename, long ID, wxString name, int nbOfTrial, int nbRun, bool rand, wxString &data);
   void AddServerJob(wxString filename, wxString name, wxString data, long id, int nbOfTrial, int nbRun, bool rand);
   void SaveDataAsFile(wxString out, wxString filename);
   bool LoadFile(wxString filename, wxString &in);
   bool isJobLoaded(long ID);
   void loadMultipleJobs(wxArrayString jobfiles);

   wxWindow          * m_parent;
   wxGrid            * m_ResultTable;
   vector<MasterResult> m_Results;
   wxGrid            * m_ClientTable;
   wxGrid            * m_JobListTable;
   wxStaticText      * m_ServerStatus;
   //wxTextCtrl         * m_EventsWindow;
   //wxTextCtrl         * m_ClientWindow;
   //wxTextCtrl         * m_JobListWindow;
   wxTimer            * m_UpdateTimer;
   //FoxServer         * m_FoxServer;
   //std::vector<FoxJob >         m_jobs;
   //std::vector<GridResult >     m_results;
   //std::vector<GridClient>      m_clients;
   //wxMutex            * m_WXFoxMasterDataMutex;
   wxString             m_working_dir;
   FoxGridMaster      * m_grid_master;

   DECLARE_EVENT_TABLE()
};
