
#include "GridResult.h"

#ifdef __WX__CRYST__
   #include "ObjCryst/wxCryst/wxCrystal.h"
#endif

GridResult::GridResult(void)
{
   filename = _T("");
   ClientName = _T("");
   threadID = -1;
   JobID = -1;
   Cost = -1;
   Show = false;
}
GridResult::~GridResult(void)
{
}
GridClient::GridClient(void) {
   id = -1;
   allCPUs = -1;
   availCPUs = -1;
}
GridClient::~GridClient(void){
}

