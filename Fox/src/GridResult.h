//#include <wx/string.h>

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

class GridResult
{
public:
   GridResult(void);
   ~GridResult(void);

   wxString   filename;
   wxString   ClientName;
   int         threadID;
   int         JobID;
   float      Cost;
   bool      Show;
};

class GridClient
{
public:
   GridClient(void);
   ~GridClient(void);

   wxString   name;
   long       id;
   long       allCPUs;
   long       availCPUs;
   wxString   status;
};
