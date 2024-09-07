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
   #include "wx/timer.h"
#endif
/*
#include <locale.h>
#include <fstream>
#include <sstream>
#include <list>
*/

#include <wx/datetime.h>
#include "ObjCryst/ObjCryst/IO.h"
#include "ObjCryst/ObjCryst/Crystal.h"
#include "ObjCryst/ObjCryst/PowderPattern.h"
#include "ObjCryst/ObjCryst/DiffractionDataSingleCrystal.h"
#include "ObjCryst/RefinableObj/GlobalOptimObj.h"
#include "foxgrid/FoxGridSlave.h"


class WXFoxSlave : public wxWindow
{
public:
   WXFoxSlave(wxWindow* parent, wxString working_dir);
   ~WXFoxSlave(void);
   void Clear();
   void ConnectClient(wxString IP);
   void OnConnectClient(wxCommandEvent& event);
   void setNbCPU(int nb);
   void CloseClient();

   wxComboBox  * m_IPWindow;
private:

   void OnConnectTimer(wxTimerEvent& event);
   void OnUpdateProcessTimer(wxTimerEvent& event);
   //void UpdateListOfProcesses(vector<FoxProcess> p);
   void InitClient();


   wxWindow     * m_parent;
   //wxTextCtrl   * m_EventsWindow;
   wxTextCtrl   * m_nbCPUs;
   //wxTextCtrl   * m_TryConnectWindow;
   wxButton     * m_ConnectButton;
   wxTimer      * m_ConnectTimer;
   wxTimer      * m_UpdateProcessTimer;
   bool           m_connecting;
   wxString       m_working_dir;
   wxGrid       * m_process_table;
   wxGrid       * m_job_table;
   FoxGridSlave * m_grid_slave;

      DECLARE_EVENT_TABLE()
};
