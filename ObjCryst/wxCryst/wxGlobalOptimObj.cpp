//#include <sstream> //for stringstream
#include <fstream>

#include "wx/wx.h"
#include "wxCryst/wxGlobalOptimObj.h"

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
//    WXGlobalOptimObj
//
////////////////////////////////////////////////////////////////////////
BEGIN_EVENT_TABLE(WXGlobalOptimObj, wxWindow)
   EVT_BUTTON(ID_WXOBJ_COLLAPSE,                       WXCrystObj::OnToggleCollapse)
   EVT_MENU(ID_REFOBJ_MENU_OBJ_SAVE,                   WXGlobalOptimObj::OnSave)
   EVT_MENU(ID_REFOBJ_MENU_OBJ_LOAD,                   WXGlobalOptimObj::OnLoad)
   EVT_MENU(ID_GLOBALOPT_MENU_GLOBAlOPT_ADDOBJ,        WXGlobalOptimObj::OnAddRefinedObject)
   EVT_MENU(ID_GLOBALOPT_MENU_GLOBAlOPT_ADDCOSTFUNC,   WXGlobalOptimObj::OnAddCostFunction)
   EVT_MENU(ID_GLOBALOPT_MENU_GLOBAlOPT_RUN,           WXGlobalOptimObj::OnRunOptimization)
   EVT_MENU(ID_GLOBALOPT_MENU_GLOBAlOPT_STOP,           WXGlobalOptimObj::OnStopOptimization)
END_EVENT_TABLE()

WXGlobalOptimObj::WXGlobalOptimObj(wxWindow* parent, GlobalOptimObj*obj):
WXCrystObj(parent),mpGlobalOptimObj(obj),mpGlobalOptimRunThread(0),mNbTrial(100000000)
{
   VFN_DEBUG_MESSAGE("WXGlobalOptimObj::WXGlobalOptimObj(wxWindow*,GlobalOptimObj*,)",6)
   mpWXTitle->SetForegroundColour(wxColour(255,0,0));
   mpWXTitle->SetLabel("Global Optimization:");
   // Menu
      mpMenuBar=new WXCrystMenuBar(this,this);
      mpSizer->Add(mpMenuBar);
      mList.Add(mpMenuBar);
      
      mpMenuBar->AddMenu("Object",ID_REFOBJ_MENU_OBJ);
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_REFOBJ_MENU_OBJ_SAVE,"Save");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_REFOBJ_MENU_OBJ_LOAD,"Load");
      mpMenuBar->AddMenu("Optimize",ID_GLOBALOPT_MENU_GLOBAlOPT);
         mpMenuBar->AddMenuItem(ID_GLOBALOPT_MENU_GLOBAlOPT,ID_GLOBALOPT_MENU_GLOBAlOPT_ADDOBJ,
                                "Add object to optimize");
         mpMenuBar->AddMenuItem(ID_GLOBALOPT_MENU_GLOBAlOPT,
                                ID_GLOBALOPT_MENU_GLOBAlOPT_ADDCOSTFUNC,"Add Cost Function");
         mpMenuBar->AddMenuItem(ID_GLOBALOPT_MENU_GLOBAlOPT,
                                ID_GLOBALOPT_MENU_GLOBAlOPT_RUN,"Run Optimization");
         mpMenuBar->AddMenuItem(ID_GLOBALOPT_MENU_GLOBAlOPT,
                                ID_GLOBALOPT_MENU_GLOBAlOPT_STOP,"Stop Optimization");
   //options
      WXFieldOption *opt;
      
      opt=new WXFieldOption(this,-1,&(mpGlobalOptimObj->mGlobalOptimType));
      mpSizer->Add(opt,0,wxALIGN_LEFT);
      mList.Add(opt);
      
      opt=new WXFieldOption(this,-1,&(mpGlobalOptimObj->mAnnealingScheduleTemp));
      mpSizer->Add(opt,0,wxALIGN_LEFT);
      mList.Add(opt);
      
      wxBoxSizer *sizerTemp=new wxBoxSizer(wxHORIZONTAL);
      WXFieldPar<double> *tempMax=
         new WXFieldPar<double>(this,"Temperature Max:",-1,&(mpGlobalOptimObj->mTemperatureMax));
      WXFieldPar<double> *tempMin=
         new WXFieldPar<double>(this,"Temperature Min:",-1,&(mpGlobalOptimObj->mTemperatureMin));
      sizerTemp->Add(tempMax);
      sizerTemp->Add(tempMin);
      mpSizer->Add(sizerTemp);
      mList.Add(tempMax);
      mList.Add(tempMin);
      
      opt=new WXFieldOption(this,-1,&(mpGlobalOptimObj->mAnnealingScheduleMutation));
      mpSizer->Add(opt,0,wxALIGN_LEFT);
      mList.Add(opt);
      
      wxBoxSizer *sizerAmp=new wxBoxSizer(wxHORIZONTAL);
      WXFieldPar<double> *ampMax=
         new WXFieldPar<double>(this,"Amplitude Max:",-1,&(mpGlobalOptimObj->mMutationAmplitudeMax));
      WXFieldPar<double> *ampMin=
         new WXFieldPar<double>(this,"Amplitude Min:",-1,&(mpGlobalOptimObj->mMutationAmplitudeMin));
      sizerAmp->Add(ampMax);
      sizerAmp->Add(ampMin);
      mpSizer->Add(sizerAmp);
      mList.Add(ampMax);
      mList.Add(ampMin);

      mpWXFieldNbTrial=new WXFieldPar<long>(this,"Number of trials to go:",-1,&mNbTrial,70);
      mpSizer->Add(mpWXFieldNbTrial);
      mList.Add(mpWXFieldNbTrial);
      
    //Refined Objects
   for(int i=0;i<mpGlobalOptimObj->mRefinedObjList.GetNb();i++)
   {
      WXFieldName *refobj=new WXFieldName(this,"Optimized object:",this,-1,300,false);
      refobj->SetValue(mpGlobalOptimObj->mRefinedObjList.GetObj(i).GetClassName()
                       +":"+mpGlobalOptimObj->mRefinedObjList.GetObj(i).GetName());
      mpSizer->Add(refobj);
      mList.Add(refobj);
   }
   // Cost Functions
   for(unsigned int i=0;i<mpGlobalOptimObj->mNbCostFunction;i++)
   {
      WXCostFunction *cf=new
         WXCostFunction(this,mpGlobalOptimObj->mpCostFunctionRefinableObj[i],
                        -1,mpGlobalOptimObj->mpCostFunctionId(i),
                        mpGlobalOptimObj->mCostFunctionWeight.data()+i);
      mpSizer->Add(cf);
      mList.Add(cf);
   }
   
   this->CrystUpdate();
   this->Layout();
}

void WXGlobalOptimObj::UpdateDisplayNbTrial()
{
   VFN_DEBUG_ENTRY("WXGlobalOptimObj::OnChangeName()",5)
   mpWXFieldNbTrial->CrystUpdate();
   VFN_DEBUG_EXIT("WXGlobalOptimObj::OnChangeName()",5)
}

void WXGlobalOptimObj::CrystUpdate()
{
   mpWXTitle->SetValue(mpGlobalOptimObj->GetName());
   this->WXCrystObj::CrystUpdate();
}

bool WXGlobalOptimObj::OnChangeName(const int id)
{
   VFN_DEBUG_MESSAGE("WXGlobalOptimObj::OnChangeName()",6)
   switch(id)
   {
      case ID_WXOBJ_NAME:
      {
      VFN_DEBUG_MESSAGE("WXGlobalOptimObj::OnChangeName():Changing GlobalOptimObj Name",6)
         mpGlobalOptimObj->SetName(mpWXTitle->GetValue());
         return true;
      }
   }
   return false;
}
void WXGlobalOptimObj::OnSave()
{
}

void WXGlobalOptimObj::OnLoad()
{
}

void WXGlobalOptimObj::OnAddRefinedObject()
{
   int choice;
   RefinableObj *obj=
      WXDialogChooseFromRegistry(gTopRefinableObjRegistry,this,
                                 "Choose object to optimize:",choice);
   if(0==obj) return;
   mpGlobalOptimObj->AddRefinableObj(*obj);
}

void WXGlobalOptimObj::AddRefinedObject(RefinableObj &obj)
{
   WXFieldName *refobj=new WXFieldName(this,"Optimized object:",this,-1,300,false);
   refobj->SetValue(obj.GetClassName()+":"+obj.GetName());
   mpSizer->Add(refobj);
   mList.Add(refobj);
   this->Layout();
}

void WXGlobalOptimObj::OnAddCostFunction()
{
   VFN_DEBUG_MESSAGE("WXGlobalOptimObj::OnAddCostFunction()",6)
   const int nbObj=mpGlobalOptimObj->mRefinedObjList.GetNb();
   int nbCostFunc=0;
   for(int i=0;i<nbObj;i++)
   {
      VFN_DEBUG_MESSAGE(mpGlobalOptimObj->mRefinedObjList.GetObj(i).GetName()<<":"<<mpGlobalOptimObj->mRefinedObjList.GetObj(i).GetNbCostFunction()<< " cost functions",6)
      nbCostFunc += mpGlobalOptimObj->mRefinedObjList.GetObj(i).GetNbCostFunction();
   }
   if(0==nbCostFunc) 
   {
      cout << nbObj <<":"<<nbCostFunc<<endl;
      return;
   }
   CrystVector_int refobj(nbCostFunc);
   CrystVector_int costFunc(nbCostFunc);
   wxString choices[10];//:TODO: use nbCostFunc
   int k=0;
   for(int i=0;i<nbObj;i++) 
    for(int j=0;j<(int) mpGlobalOptimObj->mRefinedObjList.GetObj(i).GetNbCostFunction();j++)
    {
      refobj(k)=i;
      costFunc(k)=j;
      choices[k++]=(mpGlobalOptimObj->mRefinedObjList.GetObj(i).GetName()+":"+
         mpGlobalOptimObj->mRefinedObjList.GetObj(i).GetCostFunctionName(j)
         ).c_str();
    }
   
   wxSingleChoiceDialog dialog
         (this,"Choose a new cost function","Choose",nbCostFunc,choices,0,wxOK | wxCANCEL);
   if(wxID_OK!=dialog.ShowModal()) return;
   int choice=dialog.GetSelection();
   mpGlobalOptimObj->AddCostFunction(
         mpGlobalOptimObj->mRefinedObjList.GetObj(refobj(choice)),
         costFunc(choice),1.0);
}
void WXGlobalOptimObj::AddCostFunction(RefinableObj &obj, const int costFunc)
{
   VFN_DEBUG_MESSAGE("WXGlobalOptimObj::AddCostFunction()",6)
   //:KLUDGE: Here we are assuming that this is the *last* function added
   WXCostFunction *cf=new
         WXCostFunction(this,&obj,-1,costFunc,
                        mpGlobalOptimObj->mCostFunctionWeight.data()
                        +mpGlobalOptimObj->mNbCostFunction-1);
   mpSizer->Add(cf);
   mList.Add(cf);
   this->Layout();
}

void WXGlobalOptimObj::OnRunOptimization()
{
   //VFN_DEBUG_GLOBAL_LEVEL(5);
   if(true==mpGlobalOptimObj->mIsOptimizing)
   {
      wxMessageDialog dumbUser(this,"The optimization is already running !","Huh ?",
                               wxOK|wxICON_EXCLAMATION);
      dumbUser.ShowModal();
      return;
   }
   mpGlobalOptimRunThread = new WXGlobalOptimRunThread(mpGlobalOptimObj,mNbTrial);
   if(mpGlobalOptimRunThread->Create() != wxTHREAD_NO_ERROR) 
      wxLogError("Can't create optimization thread");
   else mpGlobalOptimRunThread->Run();
   
}
void WXGlobalOptimObj::OnStopOptimization()
{
   mpGlobalOptimObj->StopAfterCycle();
}

////////////////////////////////////////////////////////////////////////
//
//    WXGlobalOptimRunThread
//
////////////////////////////////////////////////////////////////////////
WXGlobalOptimRunThread::WXGlobalOptimRunThread(GlobalOptimObj *globalOptObj,long &nbTrial):
wxThread(wxTHREAD_DETACHED),mpGlobalOptObj(globalOptObj),mpNbTrial(&nbTrial)
{
}

void *WXGlobalOptimRunThread::Entry()
{
   cout <<endl<<"Entering refinement thread "<<endl<<endl;
   mpGlobalOptObj->Optimize(*mpNbTrial);
   return NULL;
}
void WXGlobalOptimRunThread::OnExit()
{
   cout <<endl<<"Exiting refinement thread "<<endl<<endl;
}

}// namespace 

