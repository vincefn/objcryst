//#include <sstream> //for stringstream
#include <fstream>

#include "wx/wx.h"

#include "wxCryst/wxCrystal.h"

#include "ObjCryst/Atom.h"
#include "ObjCryst/ZScatterer.h"

extern "C" {
#include "GL/glu.h"
#include "wxCryst/trackball.h"
}

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

//#include "ObjCryst/Map.cpp"

namespace ObjCryst
{

////////////////////////////////////////////////////////////////////////
//
//    WXCrystal
//
////////////////////////////////////////////////////////////////////////
BEGIN_EVENT_TABLE(WXCrystal,wxEvtHandler)
   EVT_BUTTON(ID_WXOBJ_COLLAPSE,                      WXCrystObj::OnToggleCollapse)
   EVT_MENU(ID_REFOBJ_MENU_OBJ_SAVE,                  WXRefinableObj::OnMenuSave)
   EVT_MENU(ID_REFOBJ_MENU_OBJ_LOAD,                  WXRefinableObj::OnMenuLoad)
   EVT_MENU(ID_CRYSTAL_MENU_SAVECIF,                WXCrystal::OnMenuSaveCIF)
   EVT_MENU(ID_CRYSTAL_MENU_SAVETEXT,                WXCrystal::OnMenuSaveText)
   EVT_MENU(ID_REFOBJ_MENU_PAR_FIXALL,                WXRefinableObj::OnMenuFixAllPar)
   EVT_MENU(ID_REFOBJ_MENU_PAR_UNFIXALL,              WXRefinableObj::OnMenuUnFixAllPar)
   EVT_MENU(ID_REFOBJ_MENU_PAR_RANDOMIZE,             WXRefinableObj::OnMenuParRandomize)
   EVT_MENU(ID_CRYSTAL_MENU_PAR_ADDANTIBUMP,          WXCrystal::OnMenuAddAntiBumpDist)
   EVT_MENU(ID_CRYSTAL_MENU_DISPLAY_3DVIEW,           WXCrystal::OnMenuCrystalGL)
   EVT_MENU(ID_CRYSTAL_MENU_SCATT_ADDSCATTPOWATOM,    WXCrystal::OnMenuAddScattPowAtom)
   EVT_MENU(ID_CRYSTAL_MENU_SCATT_REMOVESCATTPOW,     WXCrystal::OnMenuRemoveScattPow)
   EVT_MENU(ID_CRYSTAL_MENU_SCATT_ADDATOM,            WXCrystal::OnMenuAddScatterer)
   EVT_MENU(ID_CRYSTAL_MENU_SCATT_ADDZSCATTERER,      WXCrystal::OnMenuAddScatterer)
   EVT_MENU(ID_CRYSTAL_MENU_SCATT_ADDTETRAHEDRON,     WXCrystal::OnMenuAddScatterer)
   EVT_MENU(ID_CRYSTAL_MENU_SCATT_ADDOCTAHEDRON,      WXCrystal::OnMenuAddScatterer)
   EVT_MENU(ID_CRYSTAL_MENU_SCATT_ADDTRIANGLE,        WXCrystal::OnMenuAddScatterer)
   EVT_MENU(ID_CRYSTAL_MENU_SCATT_ADDSQUAREPLANE,     WXCrystal::OnMenuAddScatterer)
   EVT_MENU(ID_CRYSTAL_MENU_SCATT_ADDCUBE,            WXCrystal::OnMenuAddScatterer)
   EVT_MENU(ID_CRYSTAL_MENU_SCATT_ADDANTIPRISMTETRAGONAL,WXCrystal::OnMenuAddScatterer)
   EVT_MENU(ID_CRYSTAL_MENU_SCATT_ADDPRISMTRIGONAL,   WXCrystal::OnMenuAddScatterer)
   EVT_MENU(ID_CRYSTAL_MENU_SCATT_ADDICOSAHEDRON,     WXCrystal::OnMenuAddScatterer)
   EVT_MENU(ID_CRYSTAL_MENU_SCATT_REMOVESCATTERER,    WXCrystal::OnMenuRemoveScatterer)
   EVT_UPDATE_UI(ID_CRYST_UPDATEUI, 						WXRefinableObj::OnUpdateUI)
END_EVENT_TABLE()

WXCrystal::WXCrystal(wxWindow* parent, Crystal *obj):
WXRefinableObj(parent,(RefinableObj*)obj),mpCrystal(obj),
mCrystalGLDisplayList(0),
/*mCrystalGLDisplayList(gGLDisplayListNb++),
mCrystalGLDisplayList(glGenLists(1)),*/
mCrystalGLDisplayListIsLocked(false),mpCrystalGL(0)
{
   VFN_DEBUG_MESSAGE("WXCrystal::WXCrystal()",6)
   //this->SetBackgroundColour("Red");
   //mpWXTitle->SetBackgroundColour(wxColour(255,200,200));
   mpWXTitle->SetForegroundColour(wxColour(255,0,0));
   // Menu
      mpMenuBar->AddMenu("Object",ID_REFOBJ_MENU_OBJ);
         //mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_REFOBJ_MENU_OBJ_LOAD,"Load");
         //mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_REFOBJ_MENU_OBJ_SAVE,"Save");
         //mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_CRYSTAL_MENU_SAVECIF,"Save (CIF)");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_OBJ,ID_CRYSTAL_MENU_SAVETEXT,"Save (Text)");
      mpMenuBar->AddMenu("Parameters",ID_REFOBJ_MENU_PAR);
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_PAR,ID_CRYSTAL_MENU_PAR_ADDANTIBUMP,
                                "Add Antibump distance");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_PAR,ID_REFOBJ_MENU_PAR_FIXALL,"Fix all");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_PAR,ID_REFOBJ_MENU_PAR_UNFIXALL,"Unfix all");
         mpMenuBar->AddMenuItem(ID_REFOBJ_MENU_PAR,ID_REFOBJ_MENU_PAR_RANDOMIZE,
                                "Randomize Configuration");
      mpMenuBar->AddMenu("Scatterers",ID_CRYSTAL_MENU_SCATT);
         mpMenuBar->AddMenuItem(ID_CRYSTAL_MENU_SCATT,ID_CRYSTAL_MENU_SCATT_ADDSCATTPOWATOM,
                                "Add Atomic Scattering Power");
         mpMenuBar->AddMenuItem(ID_CRYSTAL_MENU_SCATT,ID_CRYSTAL_MENU_SCATT_REMOVESCATTPOW,
                                "Remove Scattering Power");
         mpMenuBar->AddMenuItem(ID_CRYSTAL_MENU_SCATT,ID_CRYSTAL_MENU_SCATT_ADDATOM,
                                "Add Atom");
         mpMenuBar->AddMenuItem(ID_CRYSTAL_MENU_SCATT,ID_CRYSTAL_MENU_SCATT_ADDZSCATTERER,
                                "Add Z-Matrix Scatterer");
         mpMenuBar->AddMenuItem(ID_CRYSTAL_MENU_SCATT,ID_CRYSTAL_MENU_SCATT_ADDTETRAHEDRON,
                                "Add Tetrahedron");
         mpMenuBar->AddMenuItem(ID_CRYSTAL_MENU_SCATT,ID_CRYSTAL_MENU_SCATT_ADDOCTAHEDRON,
                                "Add Octahedron");
         mpMenuBar->AddMenuItem(ID_CRYSTAL_MENU_SCATT,ID_CRYSTAL_MENU_SCATT_ADDTRIANGLE,
                                "Add Triangle Plane");
         mpMenuBar->AddMenuItem(ID_CRYSTAL_MENU_SCATT,ID_CRYSTAL_MENU_SCATT_ADDSQUAREPLANE,
                                "Add Square Plane");
         mpMenuBar->AddMenuItem(ID_CRYSTAL_MENU_SCATT,ID_CRYSTAL_MENU_SCATT_ADDCUBE,
                                "Add Cube");
         mpMenuBar->AddMenuItem(ID_CRYSTAL_MENU_SCATT,
                                ID_CRYSTAL_MENU_SCATT_ADDANTIPRISMTETRAGONAL,
                                "Add Antiprism Tetragonal");
         mpMenuBar->AddMenuItem(ID_CRYSTAL_MENU_SCATT,ID_CRYSTAL_MENU_SCATT_ADDPRISMTRIGONAL,
                                "Add Prism Trigonal");
         mpMenuBar->AddMenuItem(ID_CRYSTAL_MENU_SCATT,ID_CRYSTAL_MENU_SCATT_ADDICOSAHEDRON,
                                "Add Icosahedron");
         mpMenuBar->AddMenuItem(ID_CRYSTAL_MENU_SCATT,ID_CRYSTAL_MENU_SCATT_REMOVESCATTERER,
                                "Remove Scatterer");
      mpMenuBar->AddMenu("Display",ID_CRYSTAL_MENU_DISPLAY);
         mpMenuBar->AddMenuItem(ID_CRYSTAL_MENU_DISPLAY,ID_CRYSTAL_MENU_DISPLAY_3DVIEW,
                                "3D Display");
   // Lattice
      wxBoxSizer* lattice=new wxBoxSizer(wxHORIZONTAL);
      mpFieldLatticeA    =new WXFieldRefPar(this,"a:",
                                     &(mpCrystal->GetPar(mpCrystal->mCellDim.data()+0)) );
          
      mpFieldLatticeB    =new WXFieldRefPar(this,"b:",
                                     &(mpCrystal->GetPar(mpCrystal->mCellDim.data()+1)) );
          
      mpFieldLatticeC    =new WXFieldRefPar(this,"c:",
                                     &(mpCrystal->GetPar(mpCrystal->mCellDim.data()+2)) );
          
      mpFieldLatticeAlpha=new WXFieldRefPar(this,"alpha:",
                                     &(mpCrystal->GetPar(mpCrystal->mCellDim.data()+3)) );
          
      mpFieldLatticeBeta =new WXFieldRefPar(this,"beta:",
                                     &(mpCrystal->GetPar(mpCrystal->mCellDim.data()+4)) );
          
      mpFieldLatticeGamma=new WXFieldRefPar(this,"gamma:",
                                     &(mpCrystal->GetPar(mpCrystal->mCellDim.data()+5)) );

      lattice->Add(mpFieldLatticeA    ,0,wxALIGN_CENTER);
      lattice->Add(mpFieldLatticeB    ,0,wxALIGN_CENTER);
      lattice->Add(mpFieldLatticeC    ,0,wxALIGN_CENTER);
      lattice->Add(mpFieldLatticeAlpha,0,wxALIGN_CENTER);
      lattice->Add(mpFieldLatticeBeta ,0,wxALIGN_CENTER);
      lattice->Add(mpFieldLatticeGamma,0,wxALIGN_CENTER);
      lattice->Layout();
      
      mpSizer->Add(lattice,0,wxALIGN_LEFT);
      mList.Add(mpFieldLatticeA);
      mList.Add(mpFieldLatticeB);
      mList.Add(mpFieldLatticeC);
      mList.Add(mpFieldLatticeAlpha);
      mList.Add(mpFieldLatticeBeta);
      mList.Add(mpFieldLatticeGamma);
      
   // SpaceGroup
      mpFieldSpacegroup=new WXFieldName(this,"SpaceGroup:",this,ID_CRYSTAL_SPACEGROUP);
      mpSizer->Add(mpFieldSpacegroup,0,wxALIGN_LEFT);
      mList.Add(mpFieldSpacegroup);
      
   // Scattering Powers
      mpWXScatteringPowerRegistry=mpCrystal
                    ->GetScatteringPowerRegistry().WXCreate(this);
      mpSizer->Add(mpWXScatteringPowerRegistry,0,wxALIGN_LEFT);
      mList.Add(mpWXScatteringPowerRegistry);
      
   // Scatterers
      mpWXScattererRegistry=mpCrystal
                    ->GetScattererRegistry().WXCreate(this);
      mpSizer->Add(mpWXScattererRegistry,0,wxALIGN_LEFT);
      mList.Add(mpWXScattererRegistry);
      
   this->CrystUpdate();
   this->Layout();
   VFN_DEBUG_MESSAGE("WXCrystal::WXCrystal():End",6)
}

void WXCrystal::CrystUpdate()
{
   VFN_DEBUG_ENTRY("WXCrystal::CrystUpdate()",7)
   this->WXRefinableObj::CrystUpdate();
   //mWXParent->Layout();
   this->UpdateGL();
   VFN_DEBUG_EXIT("WXCrystal::CrystUpdate():End",7)
}

void WXCrystal::UpdateGL(const bool onlyIndependentAtoms,
                         const REAL xMin,const REAL xMax,
                         const REAL yMin,const REAL yMax,
                         const REAL zMin,const REAL zMax)
{
   VFN_DEBUG_ENTRY("WXCrystal::UpdateGL()",8)
	WXCrystValidateAllUserInput();
   #ifdef OBJCRYST_GL
   if(mpCrystalGL!=0)
   {
   	VFN_DEBUG_MESSAGE("WXCrystal::UpdateGL():mpCrystalGL",7)
   	this->GrabCrystalGLDisplayList();
   	if(mCrystalGLDisplayList==0) mCrystalGLDisplayList=glGenLists(1);
      
		// During a refinement (multi-threaded)
		// Wait until the display list has been updated by the main thread...
		static bool cont;//:TODO: not static, but mutable member function (if >1 crystal,...)
		if(false==wxThread::IsMain())
		{
			this->ReleaseCrystalGLDisplayList();
			cont=false;
			wxCommandEvent event(wxEVT_COMMAND_MENU_SELECTED,ID_GLCRYSTAL_MENU_UPDATE);
   		wxPostEvent(mpCrystalGL,event);
			while(!cont) wxUsleep(10);
   		VFN_DEBUG_EXIT("WXCrystal::UpdateGL()-Not in main thread :End",8)
			return;
		}
   	mpCrystalGL->SetCurrent();
      glNewList(mCrystalGLDisplayList,GL_COMPILE);
         glPushMatrix();
            mpCrystal->GLInitDisplayList(onlyIndependentAtoms,xMin,xMax,yMin,yMax,zMin,zMax);
				//ScatteringPowerMap map1(mpCrystal->GetScatteringPowerRegistry().GetObj(0),
				//									 *mpCrystal,.02,.05,.05,RAD_XRAY);
				//map1.GLInitDisplayList(xMin,xMax,yMin,yMax,zMin,zMax);
				//UnitCellScattererDensityMap map2(*mpCrystal,21,21,21);
				//cout << map2.GetMap3D()<<endl;
				//map2.GLInitDisplayList(xMin,xMax,yMin,yMax,zMin,zMax);
         glPopMatrix();
      glEndList();
      //#ifdef __WINDOWS__
      cont=true;
		//#endif
      this->ReleaseCrystalGLDisplayList();
      mpCrystalGL->CrystUpdate();
   }
   else
   {
   	VFN_DEBUG_MESSAGE("WXCrystal::UpdateGL():No mpCrystalGL",7)
   }
   #endif
   VFN_DEBUG_EXIT("WXCrystal::UpdateGL():End",8)
}

int WXCrystal::GrabCrystalGLDisplayList()const
{
   VFN_DEBUG_MESSAGE("WXCrystal::GrabCrystalGLDisplayList()",7)
   //:KLUDGE: ? or OK ?
   while(mCrystalGLDisplayListIsLocked) wxUsleep(5);
   mCrystalGLDisplayListIsLocked=true;
   VFN_DEBUG_MESSAGE("WXCrystal::GrabCrystalGLDisplayList():"<<mCrystalGLDisplayList,7)
   return mCrystalGLDisplayList;
}
void WXCrystal::ReleaseCrystalGLDisplayList()const
{
   VFN_DEBUG_MESSAGE("WXCrystal::ReleaseCrystalGLDisplayList()",7)
   mCrystalGLDisplayListIsLocked=false;
}
bool WXCrystal::GLDisplayListIsLocked()const
{
   return mCrystalGLDisplayListIsLocked;
}

void WXCrystal::OnMenuCrystalGL(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXCrystal::OnMenuCrystalGL()",6)
   if(mpCrystalGL!=0) return;
   wxFrame* frame= new wxFrame(this,-1,mpCrystal->GetName().c_str());
   mpCrystalGL=new WXGLCrystalCanvas(this,frame,-1);
   frame->Show(true);
   this->UpdateGL();
}

void WXCrystal::NotifyCrystalGLDelete()
{
   VFN_DEBUG_MESSAGE("WXCrystal::NotifyCrystalGLDelete()",7)
   mpCrystalGL=0;
}

void WXCrystal::OnMenuSaveCIF(wxCommandEvent & WXUNUSED(event))
{
	WXCrystValidateAllUserInput();
   wxFileDialog save(this,"Choose a file","","","*.cif",wxSAVE | wxOVERWRITE_PROMPT);
   if(save.ShowModal() != wxID_OK) return;
   
   ofstream out(save.GetPath().c_str());
   if(!out) return;//:TODO:
   mpCrystal->CIFOutput(out);
   out.close();
}

void WXCrystal::OnMenuSaveText(wxCommandEvent & WXUNUSED(event))
{
	WXCrystValidateAllUserInput();
   wxFileDialog save(this,"Choose a file","","","*.txt",wxSAVE | wxOVERWRITE_PROMPT);
   if(save.ShowModal() != wxID_OK) return;
   
   ofstream out(save.GetPath().c_str());
   if(!out) return;//:TODO:
   mpCrystal->Print(out);
   mpCrystal->PrintMinDistanceTable(.05,out);
   out.close();
}

void WXCrystal::OnMenuAddScattPowAtom(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXCrystal::OnMenuAddScattPowAtom()",6)
   ScatteringPowerAtom *scatt=new ScatteringPowerAtom("Change me","H");
   mpCrystal->AddScatteringPower(scatt);
   VFN_DEBUG_MESSAGE("WXCrystal::OnMenuAddScattPowAtom():End",6)
   this->Layout();
}
void WXCrystal::OnMenuRemoveScattPow(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXCrystal::OnButtonRemoveScattPow()",6)
	WXCrystValidateAllUserInput();
   int choice;
   ScatteringPower *scatt=
      WXDialogChooseFromRegistry(mpCrystal->GetScatteringPowerRegistry(),this,
                                 "Choose Scattering Power to remove:",choice);
   if(0==scatt) return;
   mpCrystal->RemoveScatteringPower(scatt);
   mpCrystal->XMLOutput(cout);
   VFN_DEBUG_MESSAGE("WXCrystal::OnButtonRemoveScattPow():End",6)
   this->Layout();
}

void WXCrystal::OnMenuAddScatterer(wxCommandEvent &event)
{
   VFN_DEBUG_ENTRY("WXCrystal::OnMenuAddScatterer()",6)
	WXCrystValidateAllUserInput();
   Scatterer *scatt;
   switch(event.GetId())
   {
      case ID_CRYSTAL_MENU_SCATT_ADDATOM:
      {
         int choice;
         ScatteringPowerAtom *scattPow=dynamic_cast<ScatteringPowerAtom*>(
            WXDialogChooseFromRegistry(mpCrystal->GetScatteringPowerRegistry(),this,
                                       "Choose an atom type (ScatteringPower):",choice));
         if(0==scattPow)
         {
            VFN_DEBUG_EXIT("WXCrystal::OnMenuAddScatterer():Canceled",6)
            return;
         }
         scatt=new Atom(0,0,0,"Change Me!",scattPow);
         break;
      }
      case ID_CRYSTAL_MENU_SCATT_ADDZSCATTERER:
      {
         scatt=new ZScatterer("Change Me!",*mpCrystal);
         break;
      }
      case ID_CRYSTAL_MENU_SCATT_ADDTETRAHEDRON:
      {
         int choice;
         //Scattering power 1
            const ScatteringPower *scattPow1=WXDialogChooseFromRegistry(
                                       mpCrystal->GetScatteringPowerRegistry(),
                                       this,"Central atom type (ScatteringPower):",choice);
            if(0==scattPow1)
            {
               VFN_DEBUG_EXIT("WXCrystal::OnMenuAddScatterer():Canceled",6)
               return;
            }
         //Scattering power 2
            const ScatteringPower *scattPow2=WXDialogChooseFromRegistry(
                                       mpCrystal->GetScatteringPowerRegistry(),
                                       this,"Corner atom type (ScatteringPower):",choice);
            if(0==scattPow2)
            {
               VFN_DEBUG_EXIT("WXCrystal::OnMenuAddScatterer():Canceled",6)
               return;
            }
         //Bond length
            wxTextEntryDialog bondLengthDialog(this,"Bond length",
                                    "Enter bond length (Angstroems)","1",wxOK | wxCANCEL);
            if(wxID_OK!=bondLengthDialog.ShowModal())
            {
               VFN_DEBUG_EXIT("WXZScatterer::OnMenuAddZAtom():Cancelled",6)
               return;
            }
            double bondLength;
            bondLengthDialog.GetValue().ToDouble(&bondLength);
            
         scatt=new ZPolyhedron(TETRAHEDRON,*mpCrystal,0,0,0,"Change Me!",
                               scattPow1,scattPow2,bondLength);
         break;
      }
      case ID_CRYSTAL_MENU_SCATT_ADDOCTAHEDRON:
      {
         int choice;
         //Scattering power 1
            const ScatteringPower *scattPow1=WXDialogChooseFromRegistry(
                                       mpCrystal->GetScatteringPowerRegistry(),
                                       this,"Central atom type (ScatteringPower):",choice);
            if(0==scattPow1)
            {
               VFN_DEBUG_EXIT("WXCrystal::OnMenuAddScatterer():Canceled",6)
               return;
            }
         //Scattering power 2
            const ScatteringPower *scattPow2=WXDialogChooseFromRegistry(
                                       mpCrystal->GetScatteringPowerRegistry(),
                                       this,"Corner atom type (ScatteringPower):",choice);
            if(0==scattPow2)
            {
               VFN_DEBUG_EXIT("WXCrystal::OnMenuAddScatterer():Canceled",6)
               return;
            }
         //Bond length
            wxTextEntryDialog bondLengthDialog(this,"Bond length",
                                    "Enter bond length (Angstroems)","1",wxOK | wxCANCEL);
            if(wxID_OK!=bondLengthDialog.ShowModal())
            {
               VFN_DEBUG_EXIT("WXZScatterer::OnMenuAddZAtom():Cancelled",6)
               return;
            }
            double bondLength;
            bondLengthDialog.GetValue().ToDouble(&bondLength);
            
         scatt=new ZPolyhedron(OCTAHEDRON,*mpCrystal,0,0,0,"Change Me!",
                               scattPow1,scattPow2,bondLength);
         break;
      }
      case ID_CRYSTAL_MENU_SCATT_ADDTRIANGLE:
      {
         VFN_DEBUG_MESSAGE("WXCrystal::OnMenuAddScatterer():Add triangle plane",6)
         int choice;
         //Scattering power 1
            const ScatteringPower *scattPow1=WXDialogChooseFromRegistry(
                                       mpCrystal->GetScatteringPowerRegistry(),
                                       this,"Central atom type (ScatteringPower):",choice);
            if(0==scattPow1)
            {
               VFN_DEBUG_EXIT("WXCrystal::OnMenuAddScatterer():Canceled",6)
               return;
            }
         //Scattering power 2
            const ScatteringPower *scattPow2=WXDialogChooseFromRegistry(
                                       mpCrystal->GetScatteringPowerRegistry(),
                                       this,"Corner atom type (ScatteringPower):",choice);
            if(0==scattPow2)
            {
               VFN_DEBUG_EXIT("WXCrystal::OnMenuAddScatterer():Canceled",6)
               return;
            }
         //Bond length
            wxTextEntryDialog bondLengthDialog(this,"Bond length",
                                    "Enter bond length (Angstroems)","1",wxOK | wxCANCEL);
            if(wxID_OK!=bondLengthDialog.ShowModal())
            {
               VFN_DEBUG_EXIT("WXZScatterer::OnMenuAddZAtom():Cancelled",6)
               return;
            }
            double bondLength;
            bondLengthDialog.GetValue().ToDouble(&bondLength);
            
         scatt=new ZPolyhedron(TRIANGLE_PLANE,*mpCrystal,0,0,0,"Change Me!",
                               scattPow1,scattPow2,bondLength);
         break;
      }
      case ID_CRYSTAL_MENU_SCATT_ADDSQUAREPLANE:
      {
         VFN_DEBUG_MESSAGE("WXCrystal::OnMenuAddScatterer():Add square plane",6)
         int choice;
         //Scattering power 1
            const ScatteringPower *scattPow1=WXDialogChooseFromRegistry(
                                       mpCrystal->GetScatteringPowerRegistry(),
                                       this,"Central atom type (ScatteringPower):",choice);
            if(0==scattPow1)
            {
               VFN_DEBUG_EXIT("WXCrystal::OnMenuAddScatterer():Canceled",6)
               return;
            }
         //Scattering power 2
            const ScatteringPower *scattPow2=WXDialogChooseFromRegistry(
                                       mpCrystal->GetScatteringPowerRegistry(),
                                       this,"Corner atom type (ScatteringPower):",choice);
            if(0==scattPow2)
            {
               VFN_DEBUG_EXIT("WXCrystal::OnMenuAddScatterer():Canceled",6)
               return;
            }
         //Bond length
            wxTextEntryDialog bondLengthDialog(this,"Bond length",
                                    "Enter bond length (Angstroems)","1",wxOK | wxCANCEL);
            if(wxID_OK!=bondLengthDialog.ShowModal())
            {
               VFN_DEBUG_EXIT("WXZScatterer::OnMenuAddZAtom():Cancelled",6)
               return;
            }
            double bondLength;
            bondLengthDialog.GetValue().ToDouble(&bondLength);
            
         scatt=new ZPolyhedron(SQUARE_PLANE,*mpCrystal,0,0,0,"Change Me!",
                               scattPow1,scattPow2,bondLength);
         break;
      }
      case ID_CRYSTAL_MENU_SCATT_ADDCUBE:
      {
         int choice;
         //Scattering power 1
            const ScatteringPower *scattPow1=WXDialogChooseFromRegistry(
                                       mpCrystal->GetScatteringPowerRegistry(),
                                       this,"Central atom type (ScatteringPower):",choice);
            if(0==scattPow1)
            {
               VFN_DEBUG_EXIT("WXCrystal::OnMenuAddScatterer():Canceled",6)
               return;
            }
         //Scattering power 2
            const ScatteringPower *scattPow2=WXDialogChooseFromRegistry(
                                       mpCrystal->GetScatteringPowerRegistry(),
                                       this,"Corner atom type (ScatteringPower):",choice);
            if(0==scattPow2)
            {
               VFN_DEBUG_EXIT("WXCrystal::OnMenuAddScatterer():Canceled",6)
               return;
            }
         //Bond length
            wxTextEntryDialog bondLengthDialog(this,"Bond length",
                                    "Enter bond length (Angstroems)","1",wxOK | wxCANCEL);
            if(wxID_OK!=bondLengthDialog.ShowModal())
            {
               VFN_DEBUG_EXIT("WXZScatterer::OnMenuAddZAtom():Cancelled",6)
               return;
            }
            double bondLength;
            bondLengthDialog.GetValue().ToDouble(&bondLength);
            
         scatt=new ZPolyhedron(CUBE,*mpCrystal,0,0,0,"Change Me!",
                               scattPow1,scattPow2,bondLength);
         break;
      }
      case ID_CRYSTAL_MENU_SCATT_ADDANTIPRISMTETRAGONAL:
      {
         int choice;
         //Scattering power 1
            const ScatteringPower *scattPow1=WXDialogChooseFromRegistry(
                                       mpCrystal->GetScatteringPowerRegistry(),
                                       this,"Central atom type (ScatteringPower):",choice);
            if(0==scattPow1)
            {
               VFN_DEBUG_EXIT("WXCrystal::OnMenuAddScatterer():Canceled",6)
               return;
            }
         //Scattering power 2
            const ScatteringPower *scattPow2=WXDialogChooseFromRegistry(
                                       mpCrystal->GetScatteringPowerRegistry(),
                                       this,"Corner atom type (ScatteringPower):",choice);
            if(0==scattPow2)
            {
               VFN_DEBUG_EXIT("WXCrystal::OnMenuAddScatterer():Canceled",6)
               return;
            }
         //Bond length
            wxTextEntryDialog bondLengthDialog(this,"Bond length",
                                    "Enter bond length (Angstroems)","1",wxOK | wxCANCEL);
            if(wxID_OK!=bondLengthDialog.ShowModal())
            {
               VFN_DEBUG_EXIT("WXZScatterer::OnMenuAddZAtom():Cancelled",6)
               return;
            }
            double bondLength;
            bondLengthDialog.GetValue().ToDouble(&bondLength);
            
         scatt=new ZPolyhedron(ANTIPRISM_TETRAGONAL,*mpCrystal,0,0,0,"Change Me!",
                               scattPow1,scattPow2,bondLength);
         break;
      }
      case ID_CRYSTAL_MENU_SCATT_ADDPRISMTRIGONAL:
      {
         int choice;
         //Scattering power 1
            const ScatteringPower *scattPow1=WXDialogChooseFromRegistry(
                                       mpCrystal->GetScatteringPowerRegistry(),
                                       this,"Central atom type (ScatteringPower):",choice);
            if(0==scattPow1)
            {
               VFN_DEBUG_EXIT("WXCrystal::OnMenuAddScatterer():Canceled",6)
               return;
            }
         //Scattering power 2
            const ScatteringPower *scattPow2=WXDialogChooseFromRegistry(
                                       mpCrystal->GetScatteringPowerRegistry(),
                                       this,"Corner atom type (ScatteringPower):",choice);
            if(0==scattPow2)
            {
               VFN_DEBUG_EXIT("WXCrystal::OnMenuAddScatterer():Canceled",6)
               return;
            }
         //Bond length
            wxTextEntryDialog bondLengthDialog(this,"Bond length",
                                    "Enter bond length (Angstroems)","1",wxOK | wxCANCEL);
            if(wxID_OK!=bondLengthDialog.ShowModal())
            {
               VFN_DEBUG_EXIT("WXZScatterer::OnMenuAddZAtom():Cancelled",6)
               return;
            }
            double bondLength;
            bondLengthDialog.GetValue().ToDouble(&bondLength);
            
         scatt=new ZPolyhedron(PRISM_TRIGONAL,*mpCrystal,0,0,0,"Change Me!",
                               scattPow1,scattPow2,bondLength);
         break;
      }
      case ID_CRYSTAL_MENU_SCATT_ADDICOSAHEDRON:
      {
         int choice;
         //Scattering power 1
            const ScatteringPower *scattPow1=WXDialogChooseFromRegistry(
                                       mpCrystal->GetScatteringPowerRegistry(),
                                       this,"Central atom type (ScatteringPower):",choice);
            if(0==scattPow1)
            {
               VFN_DEBUG_EXIT("WXCrystal::OnMenuAddScatterer():Canceled",6)
               return;
            }
         //Scattering power 2
            const ScatteringPower *scattPow2=WXDialogChooseFromRegistry(
                                       mpCrystal->GetScatteringPowerRegistry(),
                                       this,"Corner atom type (ScatteringPower):",choice);
            if(0==scattPow2)
            {
               VFN_DEBUG_EXIT("WXCrystal::OnMenuAddScatterer():Canceled",6)
               return;
            }
         //Bond length
            wxTextEntryDialog bondLengthDialog(this,"Bond length",
                                    "Enter bond length (Angstroems)","1",wxOK | wxCANCEL);
            if(wxID_OK!=bondLengthDialog.ShowModal())
            {
               VFN_DEBUG_EXIT("WXZScatterer::OnMenuAddZAtom():Cancelled",6)
               return;
            }
            double bondLength;
            bondLengthDialog.GetValue().ToDouble(&bondLength);
            
         scatt=new ZPolyhedron(ICOSAHEDRON,*mpCrystal,0,0,0,"Change Me!",
                               scattPow1,scattPow2,bondLength);
         break;
      }
   }
   mpCrystal->AddScatterer(scatt);
   //mpCrystal->XMLOutput(cout);
   this->CrystUpdate();
   //this->Layout();
   VFN_DEBUG_EXIT("WXCrystal::OnMenuAddScatterer():End",6)
}

void WXCrystal::OnMenuRemoveScatterer(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXCrystal::OnButtonRemoveScatterer()",6)
	WXCrystValidateAllUserInput();
   int choice;
   Scatterer *scatt=WXDialogChooseFromRegistry(mpCrystal->GetScattererRegistry(),this,
                                             "Select the Scatterer to remove:",choice);
   if(0==scatt) return;
   mpCrystal->RemoveScatterer(scatt);
   VFN_DEBUG_MESSAGE("WXCrystal::OnButtonRemoveScatterer():End",6)
   mpCrystal->XMLOutput(cout);
   this->Layout();
}

void WXCrystal::OnMenuAddAntiBumpDist(wxCommandEvent & WXUNUSED(event))
{
	WXCrystValidateAllUserInput();
   int choice;
   //Scattering power 1
      const ScatteringPower *scattPow1=WXDialogChooseFromRegistry(
                                 mpCrystal->GetScatteringPowerRegistry(),
                                 this,"Atom type (ScatteringPower) #1",choice);
      if(0==scattPow1)
      {
         VFN_DEBUG_EXIT("WXCrystal::OnMenuAddAntiBumpDist():Canceled",6)
         return;
      }
   //Scattering power 2
      const ScatteringPower *scattPow2=WXDialogChooseFromRegistry(
                                 mpCrystal->GetScatteringPowerRegistry(),
                                 this,"Atom type (ScatteringPower) #2",choice);
      if(0==scattPow2)
      {
         VFN_DEBUG_EXIT("WXCrystal::OnMenuAddAntiBumpDist():Canceled",6)
         return;
      }
   //Distance
      wxTextEntryDialog bondLengthDialog(this,"Antibump Distance",
                              "Enter antibmup distance (<bond length) (Angstroems)",
                              "1.2",wxOK | wxCANCEL);
      if(wxID_OK!=bondLengthDialog.ShowModal())
      {
         VFN_DEBUG_EXIT("WXCrystal::OnMenuAddAntiBumpDist():Cancelled",6)
         return;
      }
      double bondLength;
      bondLengthDialog.GetValue().ToDouble(&bondLength);
      
   mpCrystal->SetBumpMergeDistance(*scattPow1,*scattPow2,bondLength);
}

bool WXCrystal::OnChangeName(const int id)
{
   VFN_DEBUG_MESSAGE("WXCrystal::OnChangeName()",6)
   if(this->WXRefinableObj::OnChangeName(id)==true) return true;
   switch(id)
   {
      case ID_CRYSTAL_SPACEGROUP:
      {
         VFN_DEBUG_MESSAGE("WXCrystal::OnChangeName():Changing SpaceGroup",6)
         mpCrystal->Init(mpCrystal->GetLatticePar(0),
                         mpCrystal->GetLatticePar(1),
                         mpCrystal->GetLatticePar(2),
                         mpCrystal->GetLatticePar(3),
                         mpCrystal->GetLatticePar(4),
                         mpCrystal->GetLatticePar(5),
                         mpFieldSpacegroup->GetValue(),
                         mpCrystal->GetName());
         mpCrystal->UpdateLatticePar();
         //Update the lattice parameters
         this->CrystUpdate();
         this->Layout();
         return true;
      }
   }
   return false;
}

void WXCrystal::UpdateUI()
{
	mpFieldSpacegroup->SetValue(mpCrystal->GetSpaceGroup().GetName());
	if(0!=mpCrystalGL) mpCrystalGL->GetParent()->SetTitle(mpCrystal->GetName().c_str());
	this->WXRefinableObj::UpdateUI();
}

////////////////////////////////////////////////////////////////////////
//
//    WXGLCrystalCanvas
//
////////////////////////////////////////////////////////////////////////
BEGIN_EVENT_TABLE(WXGLCrystalCanvas, wxGLCanvas)
   EVT_SIZE             (WXGLCrystalCanvas::OnSize)
   EVT_PAINT            (WXGLCrystalCanvas::OnPaint)
   EVT_ERASE_BACKGROUND (WXGLCrystalCanvas::OnEraseBackground)
   EVT_MOUSE_EVENTS     (WXGLCrystalCanvas::OnMouse)
   EVT_MENU(ID_GLCRYSTAL_MENU_UPDATE,WXGLCrystalCanvas::OnUpdate)
   EVT_MENU(ID_GLCRYSTAL_MENU_CHANGELIMITS,WXGLCrystalCanvas::OnChangeLimits)
   EVT_CHAR             (WXGLCrystalCanvas::OnKeyDown)
   EVT_KEY_DOWN         (WXGLCrystalCanvas::OnKeyDown)
   EVT_KEY_UP           (WXGLCrystalCanvas::OnKeyUp)
END_EVENT_TABLE()

WXGLCrystalCanvas::WXGLCrystalCanvas(WXCrystal *wxcryst,
                                     wxFrame *parent, wxWindowID id,
                                     const wxPoint &pos,
                                     const wxSize &size):
wxGLCanvas(parent,id,pos,size,wxDEFAULT_FRAME_STYLE),//
mpWXCrystal(wxcryst),mIsGLInit(false),mDist(60),mViewAngle(15),
mXmin(-.1),mXmax(1.1),mYmin(-.1),mYmax(1.1),mZmin(-.1),mZmax(1.1)
{
   VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::WXGLCrystalCanvas()",3)
   mpPopUpMenu=new wxMenu("Crystal");
   mpPopUpMenu->Append(ID_GLCRYSTAL_MENU_UPDATE, "&Update");
   mpPopUpMenu->Append(ID_GLCRYSTAL_MENU_CHANGELIMITS, "Change display &Limits");
}

WXGLCrystalCanvas::~WXGLCrystalCanvas()
{
   mpWXCrystal->NotifyCrystalGLDelete();
}

void WXGLCrystalCanvas::OnExit(wxCommandEvent &event)
{
   
}

void WXGLCrystalCanvas::OnPaint(wxPaintEvent &event)
{
   VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnPaint()",7)
	// This means that another update of the display list is being done, so...
	if(true==mpWXCrystal->GLDisplayListIsLocked()) return;
	
   wxPaintDC dc(this);
   PrepareDC(dc);
   this->GetParent()->PrepareDC(dc);

   #ifndef __WXMOTIF__
   if (!GetContext()) return;
   #endif

   this->SetCurrent();
   if(false==mIsGLInit)
   {
      mIsGLInit=true;
      this->InitGL();
   }

   //glMatrixMode( GL_PROJECTION );
   //glLoadIdentity();
   //gluPerspective( 15, .5, 1, 100 );
   glMatrixMode( GL_MODELVIEW );

   //clear
   glClearColor(0.0f, 0.0f, 0.0f, 1.0f);//Black background
   glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

   //Orientation using the trackball
   GLfloat m[4][4];
   glLoadIdentity();
   glTranslatef( 0, 0, -mDist );
   build_rotmatrix( m,mQuat);
   glMultMatrixf( &m[0][0] );
   
   //Draw
		// another update of the display list is being done, so...
		if(true==mpWXCrystal->GLDisplayListIsLocked()) return;
	
   	glCallList(mpWXCrystal->GrabCrystalGLDisplayList());  //Draw Crystal
   	mpWXCrystal->ReleaseCrystalGLDisplayList();
   
   glFlush();
   SwapBuffers();
   VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnPaint():End",7)
}

void WXGLCrystalCanvas::OnSize(wxSizeEvent& event)
{
   VFN_DEBUG_ENTRY("WXGLCrystalCanvas::OnSize()",7)
   int width, height;
   GetClientSize(& width, & height);

   #ifndef __WXMOTIF__
   if (GetContext())
   #endif
   {
      SetCurrent();
      glViewport(0, 0, width, height);
      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();
   	VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnSize():"<<mViewAngle<<","<<width<<","<<height<<","<<mDist,2)
      if( (width>0)&&(height>0)) //in case the window is docked...
      	gluPerspective(mViewAngle,(float)width/(float)height,1.f,2.*mDist);	
   }
   VFN_DEBUG_EXIT("WXGLCrystalCanvas::OnSize():End",7)
}

void WXGLCrystalCanvas::OnEraseBackground(wxEraseEvent& event)
{
}

void WXGLCrystalCanvas::OnKeyDown(wxKeyEvent& event)
{
   VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnKeyDown()",2)
   switch(event.GetKeyCode())
   {
      case(45):// +
      {
         VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnKeyDown():Bigger",2)
         int width, height;
         GetClientSize(& width, & height);
         mDist *= 1.05;
         SetCurrent();
         glMatrixMode(GL_PROJECTION);
         glLoadIdentity();
         if( (width>0)&&(height>0)) //in case size is null...
         	gluPerspective(mViewAngle,(float)width/(float)height,
                        (mDist>100)?(mDist-100):1.,mDist+100);	
         Refresh(FALSE);
         break;
      }
      case(43):// -
      {
         VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnKeyDown():Smaller",2)
         int width, height;
         GetClientSize(& width, & height);
         mDist *= .95;
         SetCurrent();
         glMatrixMode(GL_PROJECTION);
         glLoadIdentity();
         if( (width>0)&&(height>0)) //in case size is null...
         	gluPerspective(mViewAngle,(float)width/(float)height,
                        (mDist>100)?(mDist-100):1.,mDist+100);	
         Refresh(FALSE);
         break;
      }
      case(52):// 4
      {
         VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnKeyDown():rotate left",2)
         float spin_quat[4];
         trackball(spin_quat,0,0,-.05,0);
         add_quats( spin_quat, mQuat, mQuat );
         Refresh(FALSE);
         break;
      }
      case(54):// 6
      {
         VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnKeyDown():rotate right",2)
         float spin_quat[4];
         trackball(spin_quat,0,0,.05,0);
         add_quats( spin_quat, mQuat, mQuat );
         Refresh(FALSE);
         break;
      }
      case(50):// 2
      {
         VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnKeyDown():rotate down",2)
         float spin_quat[4];
         trackball(spin_quat,0,0,0,-.05);
         add_quats( spin_quat, mQuat, mQuat );
         Refresh(FALSE);
         break;
      }
      case(56):// 8
      {
         VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnKeyDown():rotate up",2)
         float spin_quat[4];
         trackball(spin_quat,0,0,0,.05);
         add_quats( spin_quat, mQuat, mQuat );
         Refresh(FALSE);
         break;
      }
      case(68):// D
      {
         VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnKeyDown():rotate left",2)
         float spin_quat[4];
         trackball(spin_quat,0,0,-.05,0);
         add_quats( spin_quat, mQuat, mQuat );
         Refresh(FALSE);
         break;
      }
      case(70):// F
      {
         VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnKeyDown():rotate right",2)
         float spin_quat[4];
         trackball(spin_quat,0,0,.05,0);
         add_quats( spin_quat, mQuat, mQuat );
         Refresh(FALSE);
         break;
      }
      case(67):// C
      {
         VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnKeyDown():rotate down",2)
         float spin_quat[4];
         trackball(spin_quat,0,0,0,-.05);
         add_quats( spin_quat, mQuat, mQuat );
         Refresh(FALSE);
         break;
      }
      case(82):// R
      {
         VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnKeyDown():rotate up",2)
         float spin_quat[4];
         trackball(spin_quat,0,0,0,.05);
         add_quats( spin_quat, mQuat, mQuat );
         Refresh(FALSE);
         break;
      }
   }
}

void WXGLCrystalCanvas::OnKeyUp(wxKeyEvent& event)
{
   VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnKeyUp():"<<event.GetKeyCode(),2)
}

void WXGLCrystalCanvas::OnEnterWindow( wxMouseEvent& event )
{
   VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnEnterWindow()",10)
}

void WXGLCrystalCanvas::OnMouse( wxMouseEvent& event )
{
   VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnMouse()",7)
   if (event.Dragging())
   {
      int width, height;
      GetClientSize(& width, & height);
      if(event.LeftIsDown())
      {
   		VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnMouse():Dragging Left Button",2)
         // drag in progress, simulate trackball
         float spin_quat[4];
         trackball(spin_quat,
         (2.0*mTrackBallLastX -       width) / (width+.001),
         (     height - 2.0*mTrackBallLastY) / (height+.001),
         (     2.0*event.GetX() - width) / (width+.001),
         (    height - 2.0*event.GetY()) / (height+.001));

         add_quats( spin_quat, mQuat, mQuat );
         Refresh(FALSE);
      }
      if(event.MiddleIsDown())
      {
   		VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnMouse():Dragging Middle Button",2)
         const float v= (mTrackBallLastY-event.GetY())/(float)height;
         const float h= (mTrackBallLastX-event.GetX())/(float)width;
         
         mDist *= (1.+v)/(1.+h);
         mViewAngle *=(1.+h);
         SetCurrent();
         glMatrixMode(GL_PROJECTION);
         glLoadIdentity();
         if( (width>0)&&(height>0)) //in case size is null...
         	gluPerspective(mViewAngle,(float)width/(float)height,
                        (mDist>100)?(mDist-100):1.,mDist+100);	
         Refresh(FALSE);
         VFN_DEBUG_MESSAGE(mViewAngle <<" "<<mDist,2)
      }
   }
   if(event.RightIsDown())
   {
   	VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnMouse():Right Button",2)
      this->PopupMenu(mpPopUpMenu, event.GetX(), event.GetY() );
   }

   mTrackBallLastX = event.GetX();
   mTrackBallLastY = event.GetY();
}

void WXGLCrystalCanvas::OnUpdate(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnUpdate()",4)
   mpWXCrystal->UpdateGL(false,mXmin,mXmax,mYmin,mYmax,mZmin,mZmax);
}

void WXGLCrystalCanvas::CrystUpdate()
{
   VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::CrystUpdate()",7)
   wxPaintEvent event;
   wxPostEvent(this,event);
}

void WXGLCrystalCanvas::InitGL()
{
   VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::InitGL()",8)
   this->SetCurrent();
    
   glEnable(GL_DEPTH_TEST);
      glEnable(GL_LIGHTING);
		
      const GLfloat color_Ambient [] = {0.3, 0.3, 0.3, 1.00}; 
      const GLfloat color_Diffuse [] = {0.6, 0.6, 0.6, 1.00}; 
      const GLfloat color_Specular[] = {0.8, 0.8, 0.8, 1.00}; 
		
		glLightfv( GL_LIGHT0, GL_AMBIENT,  color_Ambient); 
		glLightfv( GL_LIGHT0, GL_DIFFUSE,  color_Diffuse); 
		glLightfv( GL_LIGHT0, GL_SPECULAR, color_Specular); 
		glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, 1.0); 

		glEnable( GL_LIGHT0 ); 

   	glEnable(GL_NORMALIZE);
      glHint(GL_PERSPECTIVE_CORRECTION_HINT,GL_FASTEST);
      glHint(GL_POLYGON_SMOOTH_HINT,GL_FASTEST);
   
   //Initialize Trackball
   trackball(mQuat,0.,0.,0.,0.);
   
   wxSizeEvent event;
   this->OnSize(event);
   
   //First display
   this->CrystUpdate();
   VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::InitGL():End",10)
}
void WXGLCrystalCanvas::OnChangeLimits(wxCommandEvent & WXUNUSED(event))
{
   VFN_DEBUG_MESSAGE("WXGLCrystalCanvas::OnChangeLimits():End",10)
   double xMin,xMax,yMin,yMax,zMin,zMax;
   {
      wxString str;
      str << mXmin;
      wxTextEntryDialog limitDialog(this,"Enter the minimum x value (reduced coords)",
                              "x min",str,wxOK | wxCANCEL);
      if(wxID_OK!=limitDialog.ShowModal())
      {
         VFN_DEBUG_EXIT("WXGLCrystalCanvas::OnChangeLimits():Cancelled",6)
         return;
      }
      limitDialog.GetValue().ToDouble(&xMin);
   }
   {
      wxString str;
      str << mXmax;
      wxTextEntryDialog limitDialog(this,"Enter the maximum x value (reduced coords)",
                              "x max",str,wxOK | wxCANCEL);
      if(wxID_OK!=limitDialog.ShowModal())
      {
         VFN_DEBUG_EXIT("WXGLCrystalCanvas::OnChangeLimits():Cancelled",6)
         return;
      }
      limitDialog.GetValue().ToDouble(&xMax);
   }
   if(xMax<=xMin)
   {
      wxMessageDialog dumbUser(this,"max <= min !!!",
                               "Whooops",wxOK|wxICON_EXCLAMATION);
      dumbUser.ShowModal();
      return;
   }
   {
      wxString str;
      str << mYmin;
      wxTextEntryDialog limitDialog(this,"Enter the minimum y value (reduced coords)",
                              "y min",str,wxOK | wxCANCEL);
      if(wxID_OK!=limitDialog.ShowModal())
      {
         VFN_DEBUG_EXIT("WXGLCrystalCanvas::OnChangeLimits():Cancelled",6)
         return;
      }
      limitDialog.GetValue().ToDouble(&yMin);
   }
   {
      wxString str;
      str << mYmax;
      wxTextEntryDialog limitDialog(this,"Enter the maximum y value (reduced coords)",
                              "y max",str,wxOK | wxCANCEL);
      if(wxID_OK!=limitDialog.ShowModal())
      {
         VFN_DEBUG_EXIT("WXGLCrystalCanvas::OnChangeLimits():Cancelled",6)
         return;
      }
      limitDialog.GetValue().ToDouble(&yMax);
   }
   if(yMax<=yMin)
   {
      wxMessageDialog dumbUser(this,"max <= min !!!",
                               "Whooops",wxOK|wxICON_EXCLAMATION);
      dumbUser.ShowModal();
      return;
   }
   {
      wxString str;
      str << mZmin;
      wxTextEntryDialog limitDialog(this,"Enter the minimum z value (reduced coords)",
                              "z min",str,wxOK | wxCANCEL);
      if(wxID_OK!=limitDialog.ShowModal())
      {
         VFN_DEBUG_EXIT("WXGLCrystalCanvas::OnChangeLimits():Cancelled",6)
         return;
      }
      limitDialog.GetValue().ToDouble(&zMin);
   }
   {
      wxString str;
      str << mZmax;
      wxTextEntryDialog limitDialog(this,"Enter the maximum z value (reduced coords)",
                              "z max",str,wxOK | wxCANCEL);
      if(wxID_OK!=limitDialog.ShowModal())
      {
         VFN_DEBUG_EXIT("WXGLCrystalCanvas::OnChangeLimits():Cancelled",6)
         return;
      }
      limitDialog.GetValue().ToDouble(&zMax);
   }
   if(zMax<=zMin)
   {
      wxMessageDialog dumbUser(this,"max <= min !!!",
                               "Whooops",wxOK|wxICON_EXCLAMATION);
      dumbUser.ShowModal();
      return;
   }
	mXmin=xMin;mXmax=xMax;
	mYmin=yMin;mYmax=yMax;
	mZmin=zMin;mZmax=zMax;
   mpWXCrystal->UpdateGL(false,mXmin,mXmax,mYmin,mYmax,mZmin,mZmax);
	this->CrystUpdate();
}

}// namespace 

