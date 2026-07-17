// Integration test that reproduces, via direct (non-GUI) API calls, the
// Cimetidine tutorial described in Fox/doc/source/tutorial-cimetidine.rst.
// That document is written for a GUI (wxWidgets Fox) user; this test
// implements the equivalent workflow purely with the ObjCryst++ API, by
// following the exact sequence of API calls made internally by the GUI
// handlers (see ObjCryst/wxCryst/wxCrystal.cpp and
// ObjCryst/wxCryst/wxPowderPattern.cpp, which are themselves thin wrappers
// around the non-GUI API exercised here).
//
// Steps performed:
//   1. Create the triclinic Crystal (unit cell + spacegroup P121/a1).
//   2. Import the Cimetidine molecule from a Fenske-Hall Z-matrix file.
//   3. Create a PowderPattern and import the x-ray data (Fullprof format).
//   4. Add an automatic Bayesian (David-Sivia) background.
//   5. Add the crystalline phase, set profile/wavelength/polarization
//      parameters and the max sin(theta)/lambda cutoff used in the tutorial.
//   6. Run a single, very short (1000-trial) Parallel Tempering global
//      optimization -- the MonteCarloObj default algorithm -- with automatic
//      LSQ refinement enabled at the end of the run.
//
// Unlike the GUI tutorial (which suggests running for a long time, e.g.
// several million trials), this test only runs 1000 trials: the goal is to
// exercise the full pipeline end-to-end quickly and check that it produces
// finite results, not to actually solve the structure.
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <list>
#include <string>

#include "ObjCryst/ObjCryst/Crystal.h"
#include "ObjCryst/ObjCryst/Molecule.h"
#include "ObjCryst/ObjCryst/PowderPattern.h"
#include "ObjCryst/ObjCryst/ReflectionProfile.h"
#include "ObjCryst/ObjCryst/ZScatterer.h"
#include "ObjCryst/RefinableObj/GlobalOptimObj.h"
#include "ObjCryst/RefinableObj/LSQNumObj.h"

#include "../unit/test_common.h"

using namespace ObjCryst;
using objcryst_test::Check;

namespace
{

// Step 1-2: create the Crystal (unit cell + spacegroup) and import the
// Cimetidine molecule from the Fenske-Hall Z-matrix file, mirroring
// WXCrystal::OnMenuImportMoleculeFromFenskeHallZMatrix() in wxCrystal.cpp.
Crystal* CreateCimetidineCrystal(const std::string& fhzPath)
{
   Crystal* cryst = new Crystal(10.394, 18.819, 6.825, 90 * DEG2RAD, 106.44 * DEG2RAD, 90 * DEG2RAD, "P121/a1");
   cryst->SetName("Cimetidine");

   std::ifstream fin(fhzPath.c_str());
   Check(static_cast<bool>(fin), "Could not open the cime.fhz Fenske-Hall Z-matrix file");

   ZScatterer scatt("Cimetidine", *cryst);
   scatt.ImportFenskeHallZMatrix(fin, false);
   fin.close();
   cryst->AddScatterer(ZScatterer2Molecule(&scatt));

   return cryst;
}

// Full-profile Le Bail fitting, reproducing the "Le Bail + Fit Profile !"
// automatic button of WXProfileFitting::OnFit() (see
// ObjCryst/wxCryst/wxPowderPattern.cpp): Le Bail extraction mode is
// activated (reflection intensities are then extracted from the data
// instead of being calculated from the -- still unrefined -- crystal
// structure), 20 initial Le Bail cycles bring the extracted intensities
// close to convergence, and then progressively more parameters are
// unfixed and least-squares-refined -- each stage followed by a couple
// more Le Bail cycles -- exactly mirroring OnFit()'s sequence: zero
// shift + constant width (W), then variable width (U,V) + fixed
// gaussian/lorentzian mix (Eta0), then variable mix (Eta1), then peak
// asymmetry (Asym0-2) + sample displacement/transparency, then
// background, then unit cell. Le Bail mode is deactivated again at the
// end, so that later steps (e.g. the global structure optimization)
// calculate intensities from the structure as usual.
void FitProfileLeBail(PowderPattern& pattern, PowderPatternDiffraction& diffData)
{
   diffData.SetExtractionMode(true, true);

   // 20 initial Le Bail cycles (10 calls x 2 cycles each), refitting the
   // overall scale factor after each call, as done in OnFit() before the
   // progressive profile fit starts.
   for(int i = 0; i < 10; ++i)
   {
      diffData.ExtractLeBail(2);
      pattern.FitScaleFactorForRw();
   }

   LSQNumObj lsq("Cimetidine Le Bail profile fitting");
   lsq.SetRefinedObj(pattern, 0, true, true);
   lsq.PrepareRefParList(true);
   // Only fit pattern-level parameters (scale, profile, background, cell),
   // not the scatterer (atom/molecule) parameters: those are irrelevant in
   // Le Bail mode, where intensities are extracted rather than calculated
   // from the structure.
   lsq.SetParIsUsed(gpRefParTypeScatt, false);
   lsq.SetParIsUsed(gpRefParTypeScattPow, false);

   std::list<RefinablePar*> vnewpar;
   std::list<const RefParType*> vnewpartype;

   // Scale factor first.
   lsq.SetParIsFixed(gpRefParTypeScattDataScale, false);
   lsq.SafeRefine(vnewpar, vnewpartype, 2, 1, true, true);

   // Stage 1: zero shift + constant Caglioti width (W).
   vnewpar.push_back(&lsq.GetCompiledRefinedObj().GetPar("Zero"));
   vnewpar.push_back(&lsq.GetCompiledRefinedObj().GetPar("W"));
   lsq.SafeRefine(vnewpar, vnewpartype, 1.01, 5, true, true);
   vnewpar.clear();
   diffData.ExtractLeBail(2);

   // Stage 2: variable Caglioti width (U,V) + fixed gaussian/lorentzian
   // mix (Eta0).
   vnewpar.push_back(&lsq.GetCompiledRefinedObj().GetPar("U"));
   vnewpar.push_back(&lsq.GetCompiledRefinedObj().GetPar("V"));
   vnewpar.push_back(&lsq.GetCompiledRefinedObj().GetPar("Eta0"));
   lsq.SafeRefine(vnewpar, vnewpartype, 1.01, 5, true, true);
   vnewpar.clear();
   diffData.ExtractLeBail(2);

   // Stage 3: variable gaussian/lorentzian mix (Eta1).
   lsq.SetParIsFixed(diffData.GetProfile().GetPar("Eta1"), false);
   lsq.SafeRefine(vnewpar, vnewpartype, 1.01, 5, true, true);
   diffData.ExtractLeBail(2);

   // Stage 4: peak asymmetry + sample displacement/transparency.
   vnewpar.push_back(&lsq.GetCompiledRefinedObj().GetPar("Asym0"));
   vnewpar.push_back(&lsq.GetCompiledRefinedObj().GetPar("Asym1"));
   vnewpar.push_back(&lsq.GetCompiledRefinedObj().GetPar("Asym2"));
   vnewpar.push_back(&lsq.GetCompiledRefinedObj().GetPar("2ThetaDispl"));
   vnewpar.push_back(&lsq.GetCompiledRefinedObj().GetPar("2ThetaTransp"));
   lsq.SafeRefine(vnewpar, vnewpartype, 1.01, 5, true, true);
   vnewpar.clear();
   diffData.ExtractLeBail(2);

   // Stage 5: background (skipping points beyond the max resolution, as
   // OnFit() does).
   lsq.SetParIsFixed(gpRefParTypeScattDataBackground, false);
   for(unsigned int i = 0; i < pattern.GetNbPowderPatternComponent(); ++i)
      if(pattern.GetPowderPatternComponent(i).GetClassName() == "PowderPatternBackground")
      {
         PowderPatternBackground* pBack
            = dynamic_cast<PowderPatternBackground*>(&(pattern.GetPowderPatternComponent(i)));
         pBack->FixParametersBeyondMaxresolution(lsq.GetCompiledRefinedObj());
      }
   for(unsigned int i = 0; i < lsq.GetCompiledRefinedObj().GetNbPar(); ++i)
      if(!lsq.GetCompiledRefinedObj().GetPar(i).IsFixed()
         && lsq.GetCompiledRefinedObj().GetPar(i).GetType() == gpRefParTypeScattDataBackground)
         vnewpar.push_back(&lsq.GetCompiledRefinedObj().GetPar(i));
   lsq.SafeRefine(vnewpar, vnewpartype, 1.01, 5, true, true);
   vnewpar.clear();
   diffData.ExtractLeBail(2);

   // Stage 6: unit cell parameters.
   vnewpartype.push_back(gpRefParTypeUnitCell);
   lsq.SafeRefine(vnewpar, vnewpartype, 1.01, 5, true, true);
   vnewpartype.clear();
   diffData.ExtractLeBail(2);

   pattern.FitScaleFactorForRw();
   std::cout << "tutorial_cimetidine: Le Bail + profile fit Rwp=" << pattern.GetRw() * 100 << "%" << std::endl;

   // Deactivate Le Bail mode: reflection intensities are once again
   // calculated from the (still unrefined) crystal structure, as needed
   // for the global structure optimization performed afterwards.
   diffData.SetExtractionMode(false);
}

// Step 3-5: create the PowderPattern, import the x-ray data, add the
// automatic Bayesian background and the crystalline phase, and set the
// profile/wavelength/polarization/sin(theta)/lambda parameters described in
// the tutorial. Mirrors WXPowderPattern::OnMenuAddCompBackgdBayesian() and
// WXPowderPattern::OnMenuAddCompCryst() in wxPowderPattern.cpp.
PowderPattern* CreateCimetidinePowderPattern(Crystal& cryst, const std::string& dataPath)
{
   PowderPattern* pattern = new PowderPattern;
   pattern->SetName("Cimetidine X-ray");
   pattern->SetRadiationType(RAD_XRAY);
   pattern->ImportPowderPatternFullprof(dataPath);

   // Automatic Bayesian (David-Sivia) background: a linear model is fitted
   // first, then a spline using the same interpolation points, matching the
   // exact sequence used by the GUI's "Add Bayesian background" menu.
   PowderPatternBackground* pBckgd = new PowderPatternBackground;
   pattern->AddPowderPatternComponent(*pBckgd);
   {
      const long nbPointSpline = 20;
      CrystVector_REAL x(nbPointSpline), backgd(nbPointSpline);
      const CrystVector_REAL& obs = pattern->GetPowderPatternObs();
      const unsigned long nbPoint = pattern->GetNbPoint();
      const REAL xmin = pattern->GetPowderPatternX()(0);
      const REAL xmax = pattern->GetPowderPatternX()(nbPoint - 1);
      for(long i = 0; i < nbPointSpline; ++i)
      {
         x(i) = xmin + (xmax - xmin) / (REAL)(nbPointSpline - 1) * (REAL)i;
         const REAL x1 = xmin + (xmax - xmin) / (REAL)(nbPointSpline - 1) * (REAL)(i - .2);
         const REAL x2 = xmin + (xmax - xmin) / (REAL)(nbPointSpline - 1) * (REAL)(i + .2);
         long n1 = (long)pattern->X2Pixel(x1);
         long n2 = (long)pattern->X2Pixel(x2);
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

   // Crystalline phase.
   PowderPatternDiffraction* diffData = new PowderPatternDiffraction;
   diffData->SetCrystal(cryst);
   pattern->AddPowderPatternComponent(*diffData);

   // Synchrotron wavelength & linear polarization rate, as in the tutorial.
   pattern->SetWavelength(1.52904);
   pattern->GetRadiation().SetLinearPolarRate(0.98);

   // Pseudo-Voigt profile: W=.001, U=V=0, Eta0=0.5. The tutorial's Caglioti
   // U/V/W value is entered in the GUI in degrees^2 (the GUI applies a
   // RAD2DEG^2 "human" scale on top of the internal radians^2 parameters --
   // see ReflectionProfilePseudoVoigt::InitParameters()), so it must be
   // converted to radians^2 here since SetProfilePar() sets the internal
   // (unscaled) value directly.
   dynamic_cast<ReflectionProfilePseudoVoigt&>(diffData->GetProfile())
      .SetProfilePar(0.001 * DEG2RAD * DEG2RAD, 0, 0, 0.5, 0.);

   pattern->Prepare();
   pattern->FitScaleFactorForRw();

   // Full-profile Le Bail fitting: since the crystal structure (molecule
   // position/orientation) is not known yet at this stage, reflection
   // intensities are extracted from the data (Le Bail mode) rather than
   // calculated from the structure, allowing the profile shape, background
   // and unit cell parameters to be refined beforehand -- exactly what the
   // GUI's "Le Bail + Fit Profile !" button (WXProfileFitting::OnFit()) does.
   FitProfileLeBail(*pattern, *diffData);

   // Restrict refinement to the low-angle data, as recommended in the tutorial.
   pattern->SetMaxSinThetaOvLambda(0.25);
   pattern->Prepare();

   return pattern;
}

void RunTutorial(const std::string& rootDir)
{
   const std::string fhzPath = rootDir + "/Fox/example/tutorial-cimetidine/cime.fhz";
   const std::string dataPath = rootDir + "/Fox/example/tutorial-cimetidine/cime.dat";

   Crystal* cryst = CreateCimetidineCrystal(fhzPath);
   Check(cryst->GetNbScatterer() > 0, "Cimetidine crystal should have a scatterer after Z-matrix import");

   PowderPattern* pattern = CreateCimetidinePowderPattern(*cryst, dataPath);
   Check(pattern->GetNbPoint() > 0, "Powder pattern should have been imported with a non-zero number of points");

   // Step 6: Global Optimization -- Parallel Tempering (the MonteCarloObj
   // default algorithm, left untouched as recommended by the tutorial), with
   // the crystal and the powder pattern as refined objects, and automatic
   // LSQ refinement enabled at the end of the run.
   MonteCarloObj mc("Cimetidine Powder Optimization");
   mc.AddRefinableObj(*cryst);
   mc.AddRefinableObj(*pattern);
   mc.GetOption("Automatic Least Squares Refinement").SetChoice(1); // At the end of each run

   mc.RandomizeStartingConfig();

   long nbTrial = 1000;
   mc.Optimize(nbTrial, true, 0);

   const REAL finalLogLikelihood = mc.GetLogLikelihood();
   Check(std::isfinite(finalLogLikelihood), "Final Log(Likelihood) after optimization should be finite");

   std::cout << "tutorial_cimetidine: final Log(Likelihood)=" << finalLogLikelihood << std::endl;
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
   std::cout << "tutorial_cimetidine: PASS" << std::endl;
   return 0;
}
