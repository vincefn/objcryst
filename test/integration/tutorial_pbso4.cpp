// Integration test that reproduces, via direct (non-GUI) API calls, the
// PbSO4 tutorial described in Fox/doc/source/tutorial-pbso4.rst: a joint
// optimization using both an x-ray and a neutron powder pattern. That
// document is written for a GUI (wxWidgets Fox) user; this test implements
// the equivalent workflow purely with the ObjCryst++ API, by following the
// exact sequence of API calls made internally by the GUI handlers (see
// ObjCryst/wxCryst/wxCrystal.cpp and ObjCryst/wxCryst/wxPowderPattern.cpp,
// which are themselves thin wrappers around the non-GUI API exercised here).
//
// Steps performed:
//   1. Create the orthorhombic Crystal (unit cell + Pnma spacegroup).
//   2. Create the Pb, S, O atom types (ScatteringPowerAtom), add the Pb
//      atom, and add an SO4 tetrahedron (MakeTetrahedron()), matching the
//      GUI's Scatterers->Add Atomic Scattering Power / Add Atom / Add
//      Tetrahedron menus.
//   3. Create two PowderPattern objects (x-ray, neutron), each importing its
//      Fullprof data file, with an automatic Bayesian background, the
//      crystalline phase, and the profile/wavelength/sin(theta)/lambda
//      parameters suggested by the tutorial.
//   4. Run a single, very short (1000-trial) Parallel Tempering global
//      optimization -- the MonteCarloObj default algorithm -- with the
//      crystal and both powder patterns as refined objects, and automatic
//      LSQ refinement enabled at the end of the run.
//
// Unlike the GUI tutorial (which suggests running until convergence, e.g.
// 50000-100000 trials), this test only runs 1000 trials: the goal is to
// exercise the joint x-ray/neutron optimization pipeline end-to-end
// quickly and check that it produces finite results, not to actually solve
// the structure.
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

#include "ObjCryst/ObjCryst/Crystal.h"
#include "ObjCryst/ObjCryst/Molecule.h"
#include "ObjCryst/ObjCryst/Polyhedron.h"
#include "ObjCryst/ObjCryst/PowderPattern.h"
#include "ObjCryst/ObjCryst/ReflectionProfile.h"
#include "ObjCryst/RefinableObj/GlobalOptimObj.h"

#include "../unit/test_common.h"

using namespace ObjCryst;
using objcryst_test::Check;

namespace
{

// Step 1-2: create the Crystal (unit cell + spacegroup), the Pb/S/O atom
// types, the Pb atom, and the SO4 tetrahedron, mirroring
// WXCrystal::OnMenuAddScatterer() (cases ID_CRYSTAL_MENU_SCATT_ADDATOM and
// ID_CRYSTAL_MENU_SCATT_ADDTETRAHEDRON) in wxCrystal.cpp.
Crystal* CreatePbSO4Crystal()
{
   Crystal* cryst = new Crystal(8.482, 5.398, 6.959, "Pnma");
   cryst->SetName("PbSO4");

   cryst->AddScatteringPower(new ScatteringPowerAtom("Pb", "Pb"));
   cryst->AddScatteringPower(new ScatteringPowerAtom("S", "S"));
   cryst->AddScatteringPower(new ScatteringPowerAtom("O", "O"));

   // Lead atom, at the default (0,0,0) position mentioned in the tutorial.
   // Dynamical Occupancy Correction (the Crystal default) takes care of the
   // special-position multiplicity automatically, so Occup=1 is used as-is.
   cryst->AddScatterer(new Atom(0., 0., 0., "Pb1", &(cryst->GetScatteringPower("Pb")), 1.));

   // SO4 tetrahedron (~1.5 Angstroem S-O bond length), with its geometry
   // held together by bond-length/bond-angle restraints.
   Molecule* so4 = MakeTetrahedron(*cryst, "SO4", &(cryst->GetScatteringPower("S")),
                                    &(cryst->GetScatteringPower("O")), 1.5);
   cryst->AddScatterer(so4);

   return cryst;
}

// Automatic Bayesian (David-Sivia) background, reproducing
// WXPowderPattern::OnMenuAddCompBackgdBayesian() in wxPowderPattern.cpp: a
// linear model is optimized first, then a spline with the same
// interpolation points.
void AddBayesianBackground(PowderPattern& pattern)
{
   PowderPatternBackground* pBckgd = new PowderPatternBackground;
   pattern.AddPowderPatternComponent(*pBckgd);
   {
      const long nbPointSpline = 20;
      CrystVector_REAL x(nbPointSpline), backgd(nbPointSpline);
      const CrystVector_REAL& obs = pattern.GetPowderPatternObs();
      const unsigned long nbPoint = pattern.GetNbPoint();
      const REAL xmin = pattern.GetPowderPatternX()(0);
      const REAL xmax = pattern.GetPowderPatternX()(nbPoint - 1);
      for(long i = 0; i < nbPointSpline; ++i)
      {
         x(i) = xmin + (xmax - xmin) / (REAL)(nbPointSpline - 1) * (REAL)i;
         const REAL x1 = xmin + (xmax - xmin) / (REAL)(nbPointSpline - 1) * (REAL)(i - .2);
         const REAL x2 = xmin + (xmax - xmin) / (REAL)(nbPointSpline - 1) * (REAL)(i + .2);
         long n1 = (long)pattern.X2Pixel(x1);
         long n2 = (long)pattern.X2Pixel(x2);
         if(n1 < 0) n1 = 0;
         if(n2 > (long)nbPoint) n2 = nbPoint;
         backgd(i) = obs(n1);
         for(long j = n1; j < n2; ++j)
            if(obs(j) < backgd(i)) backgd(i) = obs(j);
      }
      pBckgd->SetInterpPoints(x, backgd);
   }
   pBckgd->UnFixAllPar();
   pBckgd->GetOption(0).SetChoice(0); // linear
   pBckgd->OptimizeBayesianBackground();
   pBckgd->GetOption(0).SetChoice(1); // spline
   pBckgd->OptimizeBayesianBackground();
   pBckgd->FixAllPar();
}

// Step 3 (x-ray): PowderPattern with the x-ray data, Cu Ka12 tube
// wavelength, the profile parameters and the 2theta-zero shift and
// sin(theta)/lambda cutoff recommended by the tutorial.
PowderPattern* CreateXRayPattern(Crystal& cryst, const std::string& dataPath)
{
   PowderPattern* pattern = new PowderPattern;
   pattern->SetName("PbSO4 X-ray");
   pattern->SetRadiationType(RAD_XRAY);
   pattern->ImportPowderPatternFullprof(dataPath);
   pattern->SetXZero(-0.02 * DEG2RAD);

   AddBayesianBackground(*pattern);

   PowderPatternDiffraction* diffData = new PowderPatternDiffraction;
   diffData->SetCrystal(cryst);
   pattern->AddPowderPatternComponent(*diffData);

   // Cu tube, Ka1+Ka2 doublet ("Ka12"), as selected by the tutorial's
   // Radiation->X-Ray Cu Tube Ka12 menu.
   pattern->SetWavelength("Cu");

   // Pseudo-Voigt profile: W=.01, U=V=0, Eta0=0.5, Eta1=0. The tutorial's
   // Caglioti U/V/W values are entered in the GUI in degrees^2 (the GUI
   // applies a RAD2DEG^2 "human" scale on top of the internal radians^2
   // parameters -- see ReflectionProfilePseudoVoigt::InitParameters()), so
   // they must be converted to radians^2 here since SetProfilePar() sets
   // the internal (unscaled) values directly.
   dynamic_cast<ReflectionProfilePseudoVoigt&>(diffData->GetProfile())
      .SetProfilePar(0.01 * DEG2RAD * DEG2RAD, 0, 0, 0.5, 0.);

   pattern->Prepare();
   pattern->FitScaleFactorForRw();

   // Restrict to the low-angle data for the global optimization.
   pattern->SetMaxSinThetaOvLambda(0.25);
   pattern->Prepare();

   return pattern;
}

// Step 3-bis (neutron): PowderPattern with the neutron data, matching the
// wavelength and profile parameters recommended by the tutorial.
PowderPattern* CreateNeutronPattern(Crystal& cryst, const std::string& dataPath)
{
   PowderPattern* pattern = new PowderPattern;
   pattern->SetName("PbSO4 Neutron");
   pattern->SetRadiationType(RAD_NEUTRON);
   pattern->ImportPowderPatternFullprof(dataPath);

   AddBayesianBackground(*pattern);

   PowderPatternDiffraction* diffData = new PowderPatternDiffraction;
   diffData->SetCrystal(cryst);
   pattern->AddPowderPatternComponent(*diffData);

   pattern->SetWavelength(1.909);

   // Pseudo-Voigt profile: W=.25, U=V=0, Eta0=0.15, Eta1=0 (degrees^2 -> radians^2,
   // see the comment in CreateXRayPattern()).
   dynamic_cast<ReflectionProfilePseudoVoigt&>(diffData->GetProfile())
      .SetProfilePar(0.25 * DEG2RAD * DEG2RAD, 0, 0, 0.15, 0.);

   pattern->Prepare();
   pattern->FitScaleFactorForRw();

   // Restrict to the low-angle data for the global optimization.
   pattern->SetMaxSinThetaOvLambda(0.25);
   pattern->Prepare();

   return pattern;
}

void RunTutorial(const std::string& rootDir)
{
   const std::string xrayPath = rootDir + "/Fox/example/tutorial-pbso4/xray.dat";
   const std::string neutronPath = rootDir + "/Fox/example/tutorial-pbso4/neutron.dat";

   Crystal* cryst = CreatePbSO4Crystal();
   Check(cryst->GetNbScatterer() == 2, "PbSO4 crystal should have 2 scatterers: the Pb atom and the SO4 tetrahedron");

   PowderPattern* xrayPattern = CreateXRayPattern(*cryst, xrayPath);
   PowderPattern* neutronPattern = CreateNeutronPattern(*cryst, neutronPath);
   Check(xrayPattern->GetNbPoint() > 0, "X-ray powder pattern should have been imported with a non-zero number of points");
   Check(neutronPattern->GetNbPoint() > 0,
         "Neutron powder pattern should have been imported with a non-zero number of points");

   // Step 4: Global Optimization -- Parallel Tempering (the MonteCarloObj
   // default algorithm, left untouched as recommended by the tutorial), with
   // the crystal and both powder patterns as refined objects, and automatic
   // LSQ refinement enabled at the end of the run.
   MonteCarloObj mc("PbSO4 Joint Powder Optimization");
   mc.AddRefinableObj(*cryst);
   mc.AddRefinableObj(*xrayPattern);
   mc.AddRefinableObj(*neutronPattern);
   mc.GetOption("Automatic Least Squares Refinement").SetChoice(1); // At the end of each run

   mc.RandomizeStartingConfig();

   long nbTrial = 1000;
   mc.Optimize(nbTrial, true, 0);

   const REAL finalLogLikelihood = mc.GetLogLikelihood();
   Check(std::isfinite(finalLogLikelihood), "Final Log(Likelihood) after optimization should be finite");
   Check(finalLogLikelihood > 0, "Final Log(Likelihood) after optimization should be strictly positive "
                                 "(a null value would mean the cost function was never actually evaluated)");

   std::cout << "tutorial_pbso4: final Log(Likelihood)=" << finalLogLikelihood << std::endl;
}

} // namespace

int main(int argc, char* argv[])
{
   if(argc < 2)
   {
      std::cerr << "Usage: " << argv[0] << " <ROOT_DIR>" << std::endl;
      return 1;
   }
   RunTutorial(argv[1]);
   std::cout << "tutorial_pbso4: PASS" << std::endl;
   return 0;
}
