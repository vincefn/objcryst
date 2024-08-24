#define __WXGRIDWINDOW__

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

   #include "wx/frame.h"
#endif


#include "wxFoxMaster.h"
#include "wxFoxSlave.h"


class WXGrigWindow: public wxWindow
{
public:
   WXGrigWindow(wxWindow *parent);
   ~WXGrigWindow(void);
   WXFoxMaster *StartServer();
   WXFoxSlave *StartClientWindow();
   //WXFoxSlave *StartClientWindow(wxString IP, int nbCPU=-1);
   void DataLoaded();
   void Clear();

   WXFoxMaster *m_WXFoxMaster;
   WXFoxSlave *m_WXFoxSlave;
   bool dataLoaded;
   wxString m_working_dir;


};
