/*  ObjCryst++ Object-Oriented Crystallographic Library
    (c) 2000-2002 Vincent Favre-Nicolin vincefn@users.sourceforge.net
        2000-2001 University of Geneva (Switzerland)

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/*
*  Global refinement example on PbSO4.
*
*/

#include <stdlib.h>
#include <fstream>
#include "Quirks/VFNStreamFormat.h"
#include "Quirks/Chronometer.h"
#include "Quirks/VFNDebug.h"
#include "RefinableObj/GlobalOptimObj.h"

#include "ObjCryst/General.h"

#include "ObjCryst/ScatteringPower.h"
#include "ObjCryst/Atom.h"
#include "ObjCryst/ZScatterer.h"
#include "ObjCryst/Crystal.h"
#include "ObjCryst/DiffractionDataSingleCrystal.h"


using namespace ObjCryst;

void testPbSO4()
{

   //Create 'real' PBSO4 structure (for reference only-no cheating !)
      Crystal cryst_orig(8.482,5.398,6.959,"Pnma");
      cryst_orig.SetName("PbSO4 (reference structure)");
      cryst_orig.GetSpaceGroup().Print();
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
         Atom *Pb=new Atom(.188,.250,.167,"Pb",ScattPowPb   ,1.);
         Atom *S=new Atom (.437,.750,.186,"S" ,ScattPowS    ,1.);
         Atom *O1=new Atom(.595,.750,.100,"O1",ScattPowO1   ,1.);
         Atom *O2=new Atom(.319,.750,.043,"O2",ScattPowO2   ,1.);
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
   //cryst.GLDisplayCrystal();
   cryst.RefinableObj::Print();
   
   
   //Create Diffraction data object
   DiffractionDataSingleCrystal data;
      data.SetWavelength("CuA1");
      cout << "Creating PBSO4 crystal..."<<endl;
      data.SetCrystal(cryst);
   
   //Experimental data
      //Integrated F(hkl)
      data.ImportHklIobsSigma("xray-single.hkl",80);//83,374
      

   //Options for faster or better calculations
      //Use Dynamical population correction for special positions / shared atoms
      data.GetCrystal().SetUseDynPopCorr(1);

   //Create the global optimization object
      MonteCarloObj globalOptObj;
      globalOptObj.AddRefinableObj(data);
      globalOptObj.AddRefinableObj(cryst);

   //Refine only positionnal parameters
      globalOptObj.FixAllPar();
      globalOptObj.SetParIsFixed(gpRefParTypeScattTransl,false);
      globalOptObj.SetParIsFixed(gpRefParTypeScattOrient,false);

   //Don't cheat ;-)
      globalOptObj.RandomizeStartingConfig();
   //Print Crystal structure
   cout << "Random starting configuration"<<endl;
   cryst.Print();
   
   
   
   //Calc intensity before doing anything
      data.GetIcalc();
      data.FitScaleFactorForRw();
      data.PrintObsCalcData();
   
   //Annealing parameters (schedule, Tmax, Tmin, displacement schedule, 
      globalOptObj.SetAlgorithmParallTempering(ANNEALING_SMART,1,.00001,
                                              ANNEALING_EXPONENTIAL,8,.125);      
      
   //Global Optimization
      //The real job-first test
      long nbTrial=50000;
      globalOptObj.Optimize(nbTrial);
      
   //Print calculated reflections
   data.PrintObsCalcData();
   
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
