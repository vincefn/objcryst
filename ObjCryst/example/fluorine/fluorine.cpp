/*
*  Minimal example on CaF2.
*
*/

#include <stdlib.h>
#include <fstream>
#include "Quirks/VFNStreamFormat.h"
#include "Quirks/Chronometer.h"
#include "Quirks/VFNDebug.h"
#include "RefinableObj/GlobalOptimObj.h"

#include "ObjCryst/General.h"

#include "ObjCryst/IO.h"
#include "ObjCryst/ScatteringPower.h"
#include "ObjCryst/Atom.h"
#include "ObjCryst/Crystal.h"
#include "ObjCryst/PowderPattern.h"


using namespace ObjCryst;

void fluorine()
{
	//Create crystal structure
      Crystal caf2(5.46,5.46,5.46,"Fm3m");
      caf2.SetName("CaF2");
		// Create atom types
		ScatteringPowerAtom *ScattPowCa=new ScatteringPowerAtom("Ca","Pb",0.6);
		ScatteringPowerAtom *ScattPowF =new ScatteringPowerAtom("F" ,"F" ,0.8);
		//add atom types to the crystal
      caf2.AddScatteringPower(ScattPowCa);
      caf2.AddScatteringPower(ScattPowF);
		//create the atoms
	   Atom *ca=new Atom(.0 ,.0 ,.0 ,"Ca",ScattPowCa ,1.);
	   Atom *f=new Atom( .25,.25,.25,"F" ,ScattPowF  ,1.);
      caf2.AddScatterer(ca);
      caf2.AddScatterer(f);
		
   //Create Diffraction data object, for Cu-Alpha1
      PowderPattern data;
      data.SetWavelength("CuA1");
		//add CaF2 as a Crystalline phase
      PowderPatternDiffraction * diffData=new PowderPatternDiffraction;
      diffData->SetCrystal(caf2);
      data.AddPowderPatternComponent(*diffData);
		//we don't have data, so just simulate (0->Pi/2)..
		//give a constant 'obs pattern of unit intensity
   	data.SetPowderPatternPar(0,M_PI/10000.,5000);
   	CrystVector_REAL obs(5000);
   	obs=1;
   	data.SetPowderPatternObs(obs);
		data.Prepare();
		
	// Save the powder pattern in text format
		data.SavePowderPattern("caf2.dat");
	// Save everything in xml so that we can reload it later
   	XMLCrystFileSaveGlobal("caf2.xml");
}

int main (int argc, char *argv[])
{
   TAU_PROFILE_SET_NODE(0); // sequential code 
   TAU_PROFILE("main()","int()",TAU_DEFAULT);

	cout << " Beginning CaF2 example...." << endl ;
   
   int level =10;
   if(argc==2)//debug level hase been supplied
   {
      level=atoi(argv[1]);
   }
   VFN_DEBUG_GLOBAL_LEVEL(level);
   
   fluorine();
   
	cout << " End of CaF2 example." << endl ;
   TAU_REPORT_STATISTICS();
	return 0;
}
