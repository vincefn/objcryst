//#include <sstream> //for stringstream
#include <fstream>

#include "wx/wx.h"
#include "wxCryst/wxScatterer.h"

//Fixes for Cygwin; where do those stupid macros come from ? Somewhere in wxMSW headers
#ifdef max
#undef max
#endif
#ifdef min
#undef min
#endif
#ifdef DrawText
#undef DrawText
#endif
 
namespace ObjCryst
{
////////////////////////////////////////////////////////////////////////
//
//    WXScatterer
//
////////////////////////////////////////////////////////////////////////
WXScatterer::WXScatterer(wxWindow* parent, Scatterer *obj):
WXRefinableObj(parent,(RefinableObj*)obj),mpScatterer(obj)
{
   VFN_DEBUG_MESSAGE("WXScatterer::WXScatterer()",6)
   mpWXTitle->SetForegroundColour(wxColour(0,100,0));
   //Lattice
      wxBoxSizer* sizer=new wxBoxSizer(wxHORIZONTAL);
      mpScatterer->RefinableObj::Print();
      mpFieldX    =new WXFieldRefPar(this,"x:",
                                     &(mpScatterer->GetPar(mpScatterer->mXYZ.data()+0)) );

      mpFieldY    =new WXFieldRefPar(this,"y:",
                                     &(mpScatterer->GetPar(mpScatterer->mXYZ.data()+1)) );

      mpFieldZ    =new WXFieldRefPar(this,"z:",
                                     &(mpScatterer->GetPar(mpScatterer->mXYZ.data()+2)) );
          
      mpFieldPopu =new WXFieldRefPar(this,"Popu:",
                                     &(mpScatterer->GetPar(&(mpScatterer->mOccupancy))) );

      sizer->Add(mpFieldX    ,0,wxALIGN_CENTER);
      sizer->Add(mpFieldY    ,0,wxALIGN_CENTER);
      sizer->Add(mpFieldZ    ,0,wxALIGN_CENTER);
      sizer->Add(mpFieldPopu ,0,wxALIGN_CENTER);
      
      mpSizer->Add(sizer,0,wxALIGN_LEFT);
      mList.Add(mpFieldX);
      mList.Add(mpFieldY);
      mList.Add(mpFieldZ);
      mList.Add(mpFieldPopu);
   
   //mWXParent->Layout();
   mpSizer->Layout();
   mpTopSizer->Layout();
   mpTopSizer->Fit(this);
   this->CrystUpdate();
}

}// namespace 

