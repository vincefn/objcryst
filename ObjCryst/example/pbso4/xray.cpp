/*
* ObjCryst++ : a Crystallographic computing library in C++
*
*  (c) 2000 Vincent FAVRE-NICOLIN
*           Laboratoire de Cristallographie
*           24, quai Ernest-Ansermet, CH-1211 Geneva 4, Switzerland
*  Contact: Vincent.Favre-Nicolin@cryst.unige.ch
*           Radovan.Cerny@cryst.unige.ch
*
*/
/*
*  Global refinement example on PbSO4.
*
*/

#include <stdlib.h>
#include <fstream>
#include "ObjCryst/ZScatterer.h"
#include "ObjCryst/Atom.h"
#include "ObjCryst/PowderPattern.h"
#include "Quirks/VFNStreamFormat.h"
#include "Quirks/Chronometer.h"
#include "Quirks/VFNDebug.h"
#include "RefinableObj/GlobalOptimObj.h"

//#ifdef __MWERKS__
//#include <Profiler.h>
//#else
//#include "Profile/Profiler.h" //TAU profiling package
//#endif

using namespace ObjCryst;

void testPbSO4()
{

   //Create 'real' PBSO4 structure (for reference only-no cheating !)
      Crystal cryst_orig(8.482,5.398,6.959,"Pnma");
      cryst_orig.SetName("PbSO4 (reference structure)");
      cryst_orig.Print();
      {//Create 'real' PBSO4 structure (for reference only)
	      ScatteringPowerAtom *ScattPowPb=new ScatteringPowerAtom("Pb","Pb",1.48);
	      ScatteringPowerAtom *ScattPowS =new ScatteringPowerAtom("S" ,"S",0.74);
	      ScatteringPowerAtom *ScattPowO1=new ScatteringPowerAtom("O1","O",1.87);
	      ScatteringPowerAtom *ScattPowO2=new ScatteringPowerAtom("O2","O",1.76);
	      ScatteringPowerAtom *ScattPowO3=new ScatteringPowerAtom("O3","O",1.34);
         cryst_orig.AddScatteringPower(ScattPowPb);
         cryst_orig.AddScatteringPower(ScattPowS);
         cryst_orig.AddScatteringPower(ScattPowO1);
         cryst_orig.AddScatteringPower(ScattPowO2);
         cryst_orig.AddScatteringPower(ScattPowO3);
	      Atom *Pb=new Atom(.188,.250,.167,"Pb",ScattPowPb   ,.5);
	      Atom *S=new Atom (.437,.750,.186,"S" ,ScattPowS    ,.5);
	      Atom *O1=new Atom(.595,.750,.100,"O1",ScattPowO1   ,.5);
	      Atom *O2=new Atom(.319,.750,.043,"O2",ScattPowO2   ,.5);
	      Atom *O3=new Atom(.415,.974,.306,"O3",ScattPowO3   ,1.);
         cryst_orig.AddScatterer(Pb);
         cryst_orig.AddScatterer(S);
         cryst_orig.AddScatterer(O1);
         cryst_orig.AddScatterer(O2);
         cryst_orig.AddScatterer(O3);
      }
   cryst_orig.PrintMinDistanceTable(.05);
   //cryst_orig.GLDisplayCrystal();
   
   //Creating PbSo4 refined crystal, either with independent atoms or a SO4 tetrahedron
   cout << "Creating PBSO4 crystal..."<<endl;
      // PbSo4 crystal, version with SO4 tetrahedron
      Crystal cryst(8.482,5.398,6.959,"Pnma");
      cryst.SetName("PbSO4");
      {//create the refined crystal
	      ScatteringPowerAtom *ScattPowPb=new ScatteringPowerAtom("Pb","Pb",1.5);
	      ScatteringPowerAtom *ScattPowS =new ScatteringPowerAtom("S" ,"S" ,1);
	      ScatteringPowerAtom *ScattPowO =new ScatteringPowerAtom("O" ,"O" ,1.5);
         cryst.AddScatteringPower(ScattPowPb);
         cryst.AddScatteringPower(ScattPowS);
         cryst.AddScatteringPower(ScattPowO);
	      Atom *Pb=new Atom(.0,.0,.0,"Pb",ScattPowPb   ,1.);
	      ZPolyhedron *SO4=new ZPolyhedron(TETRAHEDRON,cryst,.1,.2,.3,"SO4",
                                          ScattPowS,ScattPowO,1.5,1);
         cryst.AddScatterer(Pb);
         //SO4.SetUseGlobalScatteringPower(true);
         cryst.AddScatterer(SO4);
      }
   cryst.RefinableObj::Print();
   
   
   //Create Diffraction data object
      PowderPattern data;
      data.SetRadiationType(RAD_XRAY);
      data.SetWavelength("CuA1");
      data.SetName("PbSO4-RR");
      data.ImportPowderPatternFullprof("xray-pattern.dat");
   
   //Components
      //diffraction
      PowderPatternDiffraction * diffData=new PowderPatternDiffraction;
      diffData->SetCrystal(cryst);
      data.AddPowderPatternComponent(*diffData);
      diffData->SetName("PbSo4-diffraction");
      //approximate (hand-determined) background
      PowderPatternBackground *backgdData= new PowderPatternBackground;
      backgdData->ImportUserBackground("xray-background.dat");
      backgdData->SetName("PbSo4-background");
      data.AddPowderPatternComponent(*backgdData);
      
   //Powder data can reasonably (?) ignore f" terms
   diffData->SetIsIgnoringImagScattFact(true);
   
   //Set sigma and weight to be used
   data.SetSigmaToSqrtIobs();
   data.SetWeightToInvSigmaSq();
   
   //Exclude the background for any statistics calculation (R, Rw)
//   data.SetStatisticsExcludeBackground(true);
   //Compare 'integrated' (sort of) intensities. Be less sensitive to profile
//   data.SetUseIntegratedRfactor(true);
   
   //Profile (approximate parameters,again)
   diffData->SetReflectionProfilePar(PROFILE_PSEUDO_VOIGT,
                                     .03*DEG2RAD*DEG2RAD,
                                     0*DEG2RAD*DEG2RAD,
                                     0*DEG2RAD*DEG2RAD,
                                     0.3,0);

   //Options for faster calculations
      //Use Dynamical population correction for special positions / shared atoms
      diffData->GetCrystal().SetUseDynPopCorr(true);
      //faster for global optimization
      data.SetUseFastLessPreciseFunc(true);

   //Create the global optimization object
      MonteCarloObj globalOptObj;
      globalOptObj.AddRefinableObj(data);
      globalOptObj.AddCostFunction(data,0,1.0);
   
   //Refine only positionnal parameters
      globalOptObj.FixAllPar();
      globalOptObj.SetParIsFixed(gpRefParTypeScattTransl,false);
      globalOptObj.SetParIsFixed(gpRefParTypeScattOrient,false);
      globalOptObj.SetParIsFixed(gpRefParTypeScattConformBondLength,false);
      globalOptObj.SetParIsFixed(gpRefParTypeScattConformBondAngle,false);
      globalOptObj.SetParIsFixed(gpRefParTypeScattConformDihedAngle,false);
      globalOptObj.SetLimitsRelative(gpRefParTypeScattConformBondLength,-.1,.1);
      globalOptObj.SetLimitsRelative(gpRefParTypeScattConformBondAngle,-10.*DEG2RAD,10.*DEG2RAD);
      globalOptObj.SetLimitsRelative(gpRefParTypeScattConformDihedAngle,-10.*DEG2RAD,10.*DEG2RAD);

   //Don't cheat ;-)
      globalOptObj.RandomizeStartingConfig();
   //Print Crystal structure
   cout << "Random starting configuration"<<endl;
   cryst.Print();
   
   //Only use 2theta<2x35 to begin
      data.SetUseOnlyLowAngleData(true,35.*DEG2RAD);
   
   //Annealing parameters (schedule, Tmax, Tmin, displacement schedule, 
      globalOptObj.SetAlgorithmSimulAnnealing(ANNEALING_EXPONENTIAL,.05,.005,
                                              ANNEALING_EXPONENTIAL,10,1,
                                              100000,.25,40000);      
   //Global Optimization
      //The real job-first test
      long nbTrial=50000;
      globalOptObj.Optimize(nbTrial);
   /*   
      //Use full powder spectrum
      data.SetUseOnlyLowAngleData(false);
      //Annealing parameters
      globalOptObj.SetAlgorithmSimulAnnealing(ANNEALING_EXPONENTIAL,.01,.001,
                                              ANNEALING_EXPONENTIAL,1,1,
                                              .30,40000,20000);
      //Finer optimization
      globalOptObj.Optimize(10000);
   */
   
   //Save obtained spectrum
   //data.SavePowderPattern("calc.out");
   
   //Print calculated reflections
   diffData->PrintFhklCalc();
   
   //Print Crystal structure
   cryst.Print();
   //Print minimum distance between different atoms 
   // (<.05 are considered identical, if same element)
   cryst.PrintMinDistanceTable(.05);
   //Also print real structure
   cryst_orig.Print();

   //Save crystal POVRay file
      //ofstream out1("crystal.pov");
      //cryst.POVRayDescription(out1);
      //out1.close();
      //ofstream out2("pbso4/crystal-real.pov");
      //cryst_orig.POVRayDescription(out2);
      //out2.close();
}


int main (int argc, char *argv[])
{
   TAU_PROFILE_SET_NODE(0); // sequential code 
   TAU_PROFILE("main()","int()",TAU_DEFAULT);

	cout << " Beginning PbSO4 example...." << endl ;
   
   int level =10;
   if(argc==2)//debug level hase been supplied
   {
      level=atoi(argv[1]);
   }
   VFN_DEBUG_GLOBAL_LEVEL(level);
   
   testPbSO4();
   
	cout << " End of PbSO4 example." << endl ;
   TAU_REPORT_STATISTICS();
	return 0;
}
