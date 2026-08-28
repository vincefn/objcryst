// Extended unit tests for PowderPattern, PowderPatternBackground,
// PowderPatternDiffraction, the ReflectionProfile* classes and the
// ScatteringCorr subclasses.
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

#include "ObjCryst/ObjCryst/Crystal.h"
#include "ObjCryst/ObjCryst/DiffractionDataSingleCrystal.h"
#include "ObjCryst/ObjCryst/PowderPattern.h"
#include "ObjCryst/ObjCryst/ReflectionProfile.h"
#include "ObjCryst/ObjCryst/ScatteringCorr.h"

#include "test_common.h"

using namespace objcryst_test;

namespace
{

constexpr REAL kPowderGroundTruthAbsTol = 1e-2f;
// Regression fixtures are generated in single precision by default, while
// libobjcryst/pyobjcryst uses double precision; keep relTol loose enough so
// the same ground-truth checks remain stable in both build modes.
constexpr REAL kPowderGroundTruthRelTol = 5e-4f;

void TestPowderPatternBackground()
{
   using namespace ObjCryst;
   PowderPattern p;
   p.SetRadiationType(RAD_XRAY);
   p.SetWavelength("CuA1");
   p.SetPowderPatternPar(0, 0.01 * DEG2RAD, 800);
   CrystVector_REAL obs(800);
   obs = 10;
   p.SetPowderPatternObs(obs);

   auto* bg = new PowderPatternBackground;
   CrystVector_REAL x(2), y(2);
   x(0) = 0;
   x(1) = 8.0 * DEG2RAD;
   y(0) = 5;
   y(1) = 7;
   bg->SetInterpPoints(x, y);
   p.AddPowderPatternComponent(*bg);
   p.Prepare();

   Check(p.GetNbPowderPatternComponent() == 1, "PowderPattern background component not added");
   Check(p.GetPowderPatternCalc().numElements() == p.GetNbPoint(), "PowderPattern calc size mismatch");
}

void TestPowderPatternBackgroundHist()
{
   using namespace ObjCryst;
   PowderPattern p;
   p.SetPowderPatternPar(0, 0.01 * DEG2RAD, 5);
   CrystVector_REAL obs(5);
   obs = 0;
   p.SetPowderPatternObs(obs);

   auto* bg = new PowderPatternBackgroundHist;
   CrystVector_REAL hist(5);
   for(long i=0;i<hist.numElements();++i) hist(i)=i+1;
   bg->SetHistogram(hist);
   p.AddPowderPatternComponent(*bg);

   Check(p.GetNbPowderPatternComponent() == 1,
         "Histogram background component not added");
   Check(!bg->IsScalable(),
         "Histogram background should use its own Scale parameter");
   Check(bg->GetNbPar() == 1 && bg->GetPar(0L).GetName() == "Scale",
         "Histogram background should have one Scale parameter");

   const CrystVector_REAL stored = bg->GetHistogram();
   for(long i=0;i<hist.numElements();++i)
      CheckNearAbs(stored(i),hist(i),1e-12,
                   "Histogram background data mismatch");

   bg->GetPar("Scale").SetValue(2.5);
   const CrystVector_REAL componentCalc = bg->GetPowderPatternCalc();
   const CrystVector_REAL patternCalc = p.GetPowderPatternCalc();
   for(long i=0;i<hist.numElements();++i)
   {
      const REAL expected=2.5*hist(i);
      CheckNearAbs(componentCalc(i),expected,1e-12,
                   "Histogram component scaling mismatch");
      CheckNearAbs(patternCalc(i),expected,1e-12,
                   "Histogram contribution to powder pattern mismatch");
   }

   hist=3;
   bg->SetHistogram(hist);
   const CrystVector_REAL updatedCalc = bg->GetPowderPatternCalc();
   for(long i=0;i<hist.numElements();++i)
      CheckNearAbs(updatedCalc(i),7.5,1e-12,
                   "Updated histogram did not invalidate calculation");

   CrystVector_REAL wrongSize(4);
   wrongSize=1;
   bg->SetHistogram(wrongSize);
   bool mismatchThrown=false;
   try
   {
      bg->GetPowderPatternCalc();
   }
   catch(const ObjCrystException&)
   {
      mismatchThrown=true;
   }
   Check(mismatchThrown,
         "Histogram size mismatch should raise ObjCrystException");
}

void TestPowderPatternBackgroundHistXML()
{
   using namespace ObjCryst;
   PowderPatternBackgroundHist source;
   source.SetName("dense background");
   CrystVector_REAL hist(4);
   for(long i=0;i<hist.numElements();++i) hist(i)=0.25*(i+1);
   source.SetHistogram(hist);
   source.GetPar("Scale").SetValue(3.5);
   source.GetPar("Scale").SetIsFixed(false);

   std::stringstream xml;
   source.XMLOutput(xml,0);
   XMLCrystTag tag(xml);
   PowderPatternBackgroundHist restored;
   restored.XMLInput(xml,tag);

   Check(restored.GetName() == source.GetName(),
         "Histogram background XML did not restore the name");
   Check(restored.GetNbPar() == 1,
         "Histogram background XML did not restore the Scale parameter");
   CheckNearAbs(restored.GetPar("Scale").GetValue(),3.5,1e-12,
                "Histogram background XML did not restore Scale");
   Check(!restored.GetPar("Scale").IsFixed(),
         "Histogram background XML did not restore Scale refinement state");
   const CrystVector_REAL restoredHist = restored.GetHistogram();
   Check(restoredHist.numElements() == hist.numElements(),
         "Histogram background XML restored the wrong histogram size");
   for(long i=0;i<hist.numElements();++i)
      CheckNearAbs(restoredHist(i),hist(i),1e-12,
                   "Histogram background XML data mismatch");
}

void TestPowderPatternDiffraction()
{
   using namespace ObjCryst;
   Crystal c = MakePbso4Crystal();

   PowderPattern p;
   p.SetRadiationType(RAD_XRAY);
   p.SetWavelength("CuA1");
   p.SetPowderPatternPar(10 * DEG2RAD, 0.02 * DEG2RAD, 4000);
   CrystVector_REAL obs(4000);
   obs = 1;
   p.SetPowderPatternObs(obs);

   auto* phase = new PowderPatternDiffraction;
   phase->SetCrystal(c);
   phase->SetReflectionProfilePar(PROFILE_PSEUDO_VOIGT, .03 * DEG2RAD * DEG2RAD, 0, 0, 0.3, 0);
   p.AddPowderPatternComponent(*phase);
   p.Prepare();

   Check(p.GetNbPowderPatternComponent() == 1, "PowderPatternDiffraction component not added");
   Check(phase->GetNbRefl() > 0, "PowderPatternDiffraction reflection list is empty");
   Check(p.GetPowderPatternCalc().max() > 0, "PowderPattern diffraction calc is empty");
}

void TestPowderPatternDiffractionMuRAbsorption()
{
   using namespace ObjCryst;
   Crystal c = MakePbso4Crystal();

   PowderPattern p;
   p.SetRadiationType(RAD_XRAY);
   p.SetWavelength("CuA1");
   p.SetPowderPatternPar(10 * DEG2RAD, 0.05 * DEG2RAD, 1000);
   CrystVector_REAL obs(1000);
   obs = 1;
   p.SetPowderPatternObs(obs);

   auto* phase = new PowderPatternDiffraction;
   phase->SetCrystal(c);
   phase->SetReflectionProfilePar(PROFILE_PSEUDO_VOIGT, .03 * DEG2RAD * DEG2RAD, 0, 0, 0.3, 0);
   p.AddPowderPatternComponent(*phase);
   Check(p.GetMuR() == 0, "Default cylinder-absorption muR should be 0");
   p.SetMuR(0.5);
   CheckNearAbs(p.GetMuR(), 0.5, 1e-6, "PowderPattern::SetMuR/GetMuR mismatch");
   p.Prepare();
   Check(p.GetPowderPatternCalc().max() > 0, "PowderPattern calc with muR absorption should be non-empty");
}

void TestPowderPatternDiffractionLeBailFhklObsSq()
{
   using namespace ObjCryst;
   Crystal c = MakePbso4Crystal();

   PowderPattern p;
   p.SetRadiationType(RAD_XRAY);
   p.SetWavelength("CuA1");
   p.SetPowderPatternPar(10 * DEG2RAD, 0.02 * DEG2RAD, 4000);
   CrystVector_REAL obs(4000);
   obs = 1;
   p.SetPowderPatternObs(obs);

   auto* phase = new PowderPatternDiffraction;
   phase->SetCrystal(c);
   phase->SetReflectionProfilePar(PROFILE_PSEUDO_VOIGT, .03 * DEG2RAD * DEG2RAD, 0, 0, 0.3, 0);
   p.AddPowderPatternComponent(*phase);
   p.Prepare();

   Check(!phase->GetExtractionMode(), "Extraction mode should be disabled by default");
   phase->SetExtractionMode(true, true);
   Check(phase->GetExtractionMode(), "Extraction mode should be enabled");
   phase->ExtractLeBail(1);
   Check(phase->HasFhklObsSq(), "Le Bail extraction should populate FhklObsSq");

   const unsigned long nbrefl = phase->GetNbReflBelowMaxSinThetaOvLambda();
   Check(nbrefl > 0, "Le Bail test requires at least one reflection");
   const CrystVector_REAL extracted = phase->GetFhklObsSq();
   Check(static_cast<unsigned long>(extracted.numElements()) == nbrefl,
         "Extracted FhklObsSq size mismatch");

   CrystVector_REAL assigned(nbrefl);
   for(unsigned long i = 0; i < nbrefl; ++i) assigned(i) = extracted(i) + 1.0;
   phase->SetFhklObsSq(assigned);
   const CrystVector_REAL assignedOut = phase->GetFhklObsSq();
   Check(static_cast<unsigned long>(assignedOut.numElements()) == nbrefl,
         "Assigned FhklObsSq size mismatch");
   const unsigned long nbcheck = (nbrefl < 5) ? nbrefl : 5;
   for(unsigned long i = 0; i < nbcheck; ++i)
      CheckNearAbs(assignedOut(i), assigned(i), 1e-10, "Assigned FhklObsSq value mismatch");

   CrystVector_REAL negative = assigned;
   negative(0) = -1.0;
   bool negativeThrown = false;
   try
   {
      phase->SetFhklObsSq(negative);
   }
   catch(const ObjCrystException&)
   {
      negativeThrown = true;
   }
   Check(negativeThrown, "SetFhklObsSq should reject negative intensities");

   phase->SetExtractionMode(false);
   Check(!phase->GetExtractionMode(), "Extraction mode should be disabled after request");
}

void TestPowderPatternImport()
{
   using namespace ObjCryst;
   PowderPattern p;
   p.SetRadiationType(RAD_XRAY);
   p.SetWavelength("CuA1");
   p.ImportPowderPatternFullprof("../../ObjCryst/example/pbso4/xray-pattern.dat");
   Check(p.GetNbPoint() > 1000, "Imported powder pattern should contain many points");
   Check(p.GetPowderPatternObs().numElements() == p.GetNbPoint(), "Imported powder pattern size mismatch");
   p.SetSigmaToSqrtIobs();
   p.SetWeightToInvSigmaSq();
   Check(p.GetPowderPatternObsSigma().numElements() == p.GetNbPoint(), "Sigma vector size mismatch");
}

void TestScatteringCorrSubclasses()
{
   using namespace ObjCryst;
   Crystal c = MakePbso4Crystal();
   DiffractionDataSingleCrystal sc(c, false);
   sc.SetRadiationType(RAD_XRAY);
   sc.SetWavelength(1.54056);
   sc.GenHKLFullSpace2(0.35, true);

   LorentzCorr lorentz(sc);
   Check(lorentz.GetCorr().numElements() == sc.GetNbRefl(), "LorentzCorr size mismatch");
   Check(lorentz.GetCorr().min() > 0, "LorentzCorr values should be positive");

   PolarizationCorr polar(sc);
   Check(polar.GetCorr().numElements() == sc.GetNbRefl(), "PolarizationCorr size mismatch");
   Check(polar.GetCorr().min() > 0, "PolarizationCorr values should be positive");

   PowderSlitApertureCorr slit(sc);
   Check(slit.GetCorr().numElements() == sc.GetNbRefl(), "PowderSlitApertureCorr size mismatch");
   Check(slit.GetCorr().min() > 0, "PowderSlitApertureCorr values should be positive");
}

void TestReflectionProfilePseudoVoigt()
{
   using namespace ObjCryst;
   ReflectionProfilePseudoVoigt profile;
   profile.SetProfilePar(.03f * DEG2RAD * DEG2RAD, 0, 0, 0.5f, 0);

   const REAL center = 30 * DEG2RAD;
   const REAL width = profile.GetFullProfileWidth(0.01f, center, 1, 0, 0);
   Check(width > 0, "PseudoVoigt full profile width should be positive");

   const CrystVector_REAL x = MakePowderPatternXGrid(center - width, width / 20, 41);
   const CrystVector_REAL prof = profile.GetProfile(x, center, 1, 0, 0);
   Check(prof.numElements() == x.numElements(), "PseudoVoigt profile size mismatch");
   Check(prof.max() > 0, "PseudoVoigt profile should be positive somewhere");
   for(long i = 0; i < prof.numElements(); ++i)
      Check(std::isfinite(prof(i)), "PseudoVoigt profile should be finite everywhere");
}

void TestReflectionProfilePseudoVoigtTCH()
{
   using namespace ObjCryst;
   ReflectionProfilePseudoVoigtTCH profile;
   Check(!profile.IsAnisotropic(), "TCH pseudo-Voigt should be isotropic");

   const REAL center=30*DEG2RAD;
   const REAL hg=.08f*DEG2RAD;
   const REAL hl=.04f*DEG2RAD;
   profile.SetProfilePar(hg*hg,0,0,0,0,hl);
   const REAL hg2=hg*hg;
   const REAL hl2=hl*hl;
   const REAL h5=hg2*hg2*hg+2.69269f*hg2*hg2*hl
                 +2.42843f*hg2*hg*hl2+4.47163f*hg2*hl2*hl
                 +0.07842f*hg*hl2*hl2+hl2*hl2*hl;
   const REAL expectedFwhm=pow(h5,(REAL).2);
   CheckNearAbs(profile.GetFullProfileWidth(.5,center,1,0,0),expectedFwhm,1e-9,
                "TCH profile should have the calculated common FWHM");

   CrystVector_REAL x(3);
   x(0)=center-expectedFwhm/2;
   x(1)=center;
   x(2)=center+expectedFwhm/2;
   const CrystVector_REAL prof=profile.GetProfile(x,center,1,0,0);
   CheckNearAbs(prof(0)/prof(1),.5,2e-5,
                "TCH profile should be half-height at -FWHM/2");
   CheckNearAbs(prof(2)/prof(1),.5,2e-5,
                "TCH profile should be half-height at +FWHM/2");

   profile.SetProfilePar(hg*hg);
   CheckNearAbs(profile.GetFullProfileWidth(.5,center,1,0,0),hg,1e-9,
                "TCH pure-Gaussian limit should retain H_G");
   profile.SetProfilePar(0,0,0,0,0,hl);
   CheckNearAbs(profile.GetFullProfileWidth(.5,center,1,0,0),hl,1e-9,
                "TCH pure-Lorentzian limit should retain H_L");
}

void TestReflectionProfileDoubleExponentialPseudoVoigt()
{
   using namespace ObjCryst;
   // This profile needs a unit cell to compute d_hkl for each reflection;
   // without one, the instrument parameters divide by a zero d-spacing.
   Crystal c = MakePbso4Crystal();
   ReflectionProfileDoubleExponentialPseudoVoigt profile(c);
   Check(!profile.IsAnisotropic(), "DoubleExponentialPseudoVoigt should not be anisotropic by default");

   const REAL center = 0;
   const REAL h = 1, k = 1, l = 1;
   // Use a modest, fixed-size grid rather than GetFullProfileWidth(): the
   // latter's search loop assumes a d-spacing scale inconsistent with
   // GetProfile()'s and is not reliable for arbitrary instrument parameters.
   const CrystVector_REAL x = MakePowderPatternXGrid(center - 1.0, 0.01, 201);
   const CrystVector_REAL prof = profile.GetProfile(x, center, h, k, l);
   Check(prof.numElements() == x.numElements(), "DoubleExponentialPseudoVoigt profile size mismatch");
   for(long i = 0; i < prof.numElements(); ++i)
      Check(std::isfinite(prof(i)), "DoubleExponentialPseudoVoigt profile should be finite everywhere");
}

void TestReflectionProfilePseudoVoigtAnisotropicDirect()
{
   using namespace ObjCryst;
   ReflectionProfilePseudoVoigtAnisotropic parameterProfile;
   const REAL w = 1e-6f;
   const REAL p = 2e-6f;
   parameterProfile.SetProfilePar(w, 0, 0, p);
   CheckNearAbs(parameterProfile.GetPar("P").GetValue(), p, 1e-12,
                "Anisotropic SetProfilePar should assign Gaussian P");
   CheckNearAbs(parameterProfile.GetPar("Asym0").GetValue(), 1, 1e-12,
                "Anisotropic SetProfilePar should default to symmetric Asym0");

   const REAL wavelength = 1.54056f;
   const REAL xMin = 10 * DEG2RAD;
   const REAL step = 0.05f * DEG2RAD;
   const unsigned long nbPoint = 240;
   const REAL xMax = xMin + step * static_cast<REAL>(nbPoint - 1);
   const CrystVector_REAL patternX = MakePowderPatternXGrid(xMin, step, nbPoint);

   Crystal c = MakePbso4Crystal();
   DiffractionDataSingleCrystal sc(c, false);
   sc.SetRadiationType(RAD_XRAY);
   sc.SetWavelength(wavelength);
   sc.GenHKLFullSpace2(0.32f, true);

   ReflectionProfilePseudoVoigtAnisotropic profile = MakePbso4AnisotropicProfile();

   const auto& h = sc.GetH();
   const auto& k = sc.GetK();
   const auto& l = sc.GetL();
   const auto& stol = sc.GetSinThetaOverLambda();

   unsigned long computedProfiles = 0;
   for(long i = 0; i < sc.GetNbRefl(); ++i)
   {
      const REAL center = StolToTwoTheta(stol(i), wavelength);
      if(center < xMin || center > xMax) continue;

      const REAL halfwidth = profile.GetFullProfileWidth(0.04f, center, h(i), k(i), l(i)) * 5.0f;
      Check(std::isfinite(halfwidth), "Direct anisotropic profile halfwidth is not finite");
      Check(halfwidth > 0, "Direct anisotropic profile halfwidth should be positive");

      long first = static_cast<long>((center - halfwidth - xMin) / step) - 1;
      long last = static_cast<long>((center + halfwidth - xMin) / step) + 1;
      if(last < 0 || first >= static_cast<long>(nbPoint)) continue;
      if(first < 0) first = 0;
      if(last >= static_cast<long>(nbPoint)) last = static_cast<long>(nbPoint) - 1;

      CrystVector_REAL vx(last - first + 1);
      for(long j = first; j <= last; ++j) vx(j - first) = patternX(j);

      const CrystVector_REAL reflProfile = profile.GetProfile(vx, center, h(i), k(i), l(i));
      Check(reflProfile.numElements() == vx.numElements(), "Direct anisotropic profile size mismatch");
      Check(reflProfile.max() > 0, "Direct anisotropic profile should have positive intensity");
      for(long j = 0; j < reflProfile.numElements(); ++j)
      {
         Check(std::isfinite(reflProfile(j)), "Direct anisotropic profile contains non-finite values");
      }
      ++computedProfiles;
   }
   Check(computedProfiles > 0, "Direct anisotropic profile test did not compute any in-range reflections");
}

void TestPowderGroundTruthXrayPseudoVoigtGaussian()
{
   using namespace ObjCryst;
   auto* profile = new ReflectionProfilePseudoVoigt();
   profile->SetProfilePar(.03f * DEG2RAD * DEG2RAD, 0, 0, 0.0f, 0);
   ComparePowderSimulationToGroundTruth(RAD_XRAY, 1.54056f, profile,
                                        "powder-groundtruth-xray-pv-gaussian",
                                        "../../test/data/ground_truth/powder_xray_pv_gaussian_pbso4.txt",
                                        kPowderGroundTruthAbsTol, kPowderGroundTruthRelTol);
}

void TestPowderGroundTruthXrayPseudoVoigtLorentzian()
{
   using namespace ObjCryst;
   auto* profile = new ReflectionProfilePseudoVoigt();
   profile->SetProfilePar(.03f * DEG2RAD * DEG2RAD, 0, 0, 1.0f, 0);
   ComparePowderSimulationToGroundTruth(RAD_XRAY, 1.54056f, profile,
                                        "powder-groundtruth-xray-pv-lorentzian",
                                        "../../test/data/ground_truth/powder_xray_pv_lorentzian_pbso4.txt",
                                        kPowderGroundTruthAbsTol, kPowderGroundTruthRelTol);
}

void TestPowderGroundTruthXrayAnisotropic()
{
   using namespace ObjCryst;
   auto* profile = new ReflectionProfilePseudoVoigtAnisotropic(MakePbso4AnisotropicProfile());
   ComparePowderSimulationToGroundTruth(RAD_XRAY, 1.54056f, profile,
                                        "powder-groundtruth-xray-anisotropic",
                                        "../../test/data/ground_truth/powder_xray_anisotropic_pbso4.txt",
                                        kPowderGroundTruthAbsTol, kPowderGroundTruthRelTol);
}

void TestPowderGroundTruthNeutronPseudoVoigtGaussian()
{
   using namespace ObjCryst;
   auto* profile = new ReflectionProfilePseudoVoigt();
   profile->SetProfilePar(.03f * DEG2RAD * DEG2RAD, 0, 0, 0.0f, 0);
   ComparePowderSimulationToGroundTruth(RAD_NEUTRON, 1.54056f, profile,
                                        "powder-groundtruth-neutron-pv-gaussian",
                                        "../../test/data/ground_truth/powder_neutron_pv_gaussian_pbso4.txt",
                                        kPowderGroundTruthAbsTol, kPowderGroundTruthRelTol);
}

void TestPowderGroundTruthNeutronPseudoVoigtLorentzian()
{
   using namespace ObjCryst;
   auto* profile = new ReflectionProfilePseudoVoigt();
   profile->SetProfilePar(.03f * DEG2RAD * DEG2RAD, 0, 0, 1.0f, 0);
   ComparePowderSimulationToGroundTruth(RAD_NEUTRON, 1.54056f, profile,
                                        "powder-groundtruth-neutron-pv-lorentzian",
                                        "../../test/data/ground_truth/powder_neutron_pv_lorentzian_pbso4.txt",
                                        kPowderGroundTruthAbsTol, kPowderGroundTruthRelTol);
}

} // namespace

int main(int argc, char* argv[])
{
   if(argc != 2)
   {
      std::cerr << "Usage: api_powderpattern <test-case>" << std::endl;
      return 2;
   }

   const std::string testName = argv[1];

   if(testName == "powderpattern-background") TestPowderPatternBackground();
   else if(testName == "powderpattern-background-hist") TestPowderPatternBackgroundHist();
   else if(testName == "powderpattern-background-hist-xml") TestPowderPatternBackgroundHistXML();
   else if(testName == "powderpattern-diffraction") TestPowderPatternDiffraction();
   else if(testName == "powderpattern-diffraction-mur") TestPowderPatternDiffractionMuRAbsorption();
   else if(testName == "powderpattern-diffraction-lebail-fhklobs") TestPowderPatternDiffractionLeBailFhklObsSq();
   else if(testName == "powderpattern-import") TestPowderPatternImport();
   else if(testName == "scatteringcorr-subclasses") TestScatteringCorrSubclasses();
   else if(testName == "reflectionprofile-pseudo-voigt") TestReflectionProfilePseudoVoigt();
   else if(testName == "reflectionprofile-pseudo-voigt-tch") TestReflectionProfilePseudoVoigtTCH();
   else if(testName == "reflectionprofile-double-exponential-pv") TestReflectionProfileDoubleExponentialPseudoVoigt();
   else if(testName == "reflectionprofile-pv-anisotropic-direct") TestReflectionProfilePseudoVoigtAnisotropicDirect();
   else if(testName == "powder-groundtruth-xray-pv-gaussian") TestPowderGroundTruthXrayPseudoVoigtGaussian();
   else if(testName == "powder-groundtruth-xray-pv-lorentzian") TestPowderGroundTruthXrayPseudoVoigtLorentzian();
   else if(testName == "powder-groundtruth-xray-anisotropic") TestPowderGroundTruthXrayAnisotropic();
   else if(testName == "powder-groundtruth-neutron-pv-gaussian") TestPowderGroundTruthNeutronPseudoVoigtGaussian();
   else if(testName == "powder-groundtruth-neutron-pv-lorentzian") TestPowderGroundTruthNeutronPseudoVoigtLorentzian();
   else
   {
      std::cerr << "Unknown test case: " << testName << std::endl;
      return 2;
   }

   return 0;
}
