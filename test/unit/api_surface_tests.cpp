#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "ObjCryst/ObjCryst/Atom.h"
#include "ObjCryst/ObjCryst/CIF.h"
#include "ObjCryst/ObjCryst/Crystal.h"
#include "ObjCryst/ObjCryst/DiffractionDataSingleCrystal.h"
#include "ObjCryst/ObjCryst/Indexing.h"
#include "ObjCryst/ObjCryst/Molecule.h"
#include "ObjCryst/ObjCryst/PowderPattern.h"
#include "ObjCryst/ObjCryst/ScatteringPower.h"
#include "ObjCryst/ObjCryst/SpaceGroup.h"
#include "ObjCryst/RefinableObj/GlobalOptimObj.h"
#include "ObjCryst/RefinableObj/LSQNumObj.h"
#include "ObjCryst/RefinableObj/RefinableObj.h"

namespace
{
void Check(const bool ok, const char* const msg)
{
   if(!ok)
   {
      std::cerr << msg << std::endl;
      std::exit(1);
   }
}

void CheckNearAbs(const REAL actual, const REAL expected, const REAL tol, const std::string& msg)
{
   if(std::fabs(actual - expected) > tol)
   {
      std::ostringstream os;
      os << msg << " (actual=" << actual << ", expected=" << expected << ", tol=" << tol << ")";
      std::cerr << os.str() << std::endl;
      std::exit(1);
   }
}

void CheckNearAbsRel(const REAL actual, const REAL expected, const REAL absTol, const REAL relTol, const std::string& msg)
{
   const REAL diff = std::fabs(actual - expected);
   const REAL scale = std::max<REAL>(1.0f, std::fabs(expected));
   const REAL tol = std::max(absTol, relTol * scale);
   if(diff > tol)
   {
      std::ostringstream os;
      os << msg << " (actual=" << actual << ", expected=" << expected
         << ", absTol=" << absTol << ", relTol=" << relTol << ")";
      std::cerr << os.str() << std::endl;
      std::exit(1);
   }
}

struct SingleCrystalPoint
{
   long h;
   long k;
   long l;
   REAL f2;
};

struct PowderPoint
{
   unsigned long idx;
   REAL intensity;
};

ObjCryst::Crystal MakePbso4Crystal();

std::vector<SingleCrystalPoint> LoadSingleCrystalGroundTruth(const std::string& path)
{
   std::ifstream in(path.c_str());
   Check(in.good(), ("Unable to open single-crystal ground truth file: " + path).c_str());

   std::vector<SingleCrystalPoint> out;
   std::string line;
   while(std::getline(in, line))
   {
      if(line.empty() || line[0] == '#') continue;
      std::istringstream is(line);
      SingleCrystalPoint p{};
      is >> p.h >> p.k >> p.l >> p.f2;
      Check(!is.fail(), ("Invalid line in single-crystal ground truth file: " + line).c_str());
      out.push_back(p);
   }
   Check(!out.empty(), ("No data in single-crystal ground truth file: " + path).c_str());
   return out;
}

std::vector<PowderPoint> LoadPowderGroundTruth(const std::string& path)
{
   std::ifstream in(path.c_str());
   Check(in.good(), ("Unable to open powder ground truth file: " + path).c_str());

   std::vector<PowderPoint> out;
   std::string line;
   while(std::getline(in, line))
   {
      if(line.empty() || line[0] == '#') continue;
      std::istringstream is(line);
      PowderPoint p{};
      is >> p.idx >> p.intensity;
      Check(!is.fail(), ("Invalid line in powder ground truth file: " + line).c_str());
      out.push_back(p);
   }
   Check(!out.empty(), ("No data in powder ground truth file: " + path).c_str());
   return out;
}

void CompareSingleCrystalToGroundTruth(const ObjCryst::DiffractionDataSingleCrystal& sc,
                                       const std::string& filePath,
                                       const REAL absTol)
{
   const std::vector<SingleCrystalPoint> gt = LoadSingleCrystalGroundTruth(filePath);
   const auto& h = sc.GetH();
   const auto& k = sc.GetK();
   const auto& l = sc.GetL();
   const auto& f2 = sc.GetFhklCalcSq();

   Check(static_cast<long>(gt.size()) <= sc.GetNbRefl(), "Ground-truth has more entries than generated reflections");

   for(size_t i = 0; i < gt.size(); ++i)
   {
      Check(static_cast<long>(h(i)) == gt[i].h, "Single-crystal H mismatch against ground truth");
      Check(static_cast<long>(k(i)) == gt[i].k, "Single-crystal K mismatch against ground truth");
      Check(static_cast<long>(l(i)) == gt[i].l, "Single-crystal L mismatch against ground truth");
      CheckNearAbs(f2(i), gt[i].f2, absTol, "Single-crystal F^2 mismatch against ground truth");
   }
}

void CompareSingleCrystalSimulationToGroundTruth(const ObjCryst::RadiationType radiation,
                                                 const REAL wavelength,
                                                 const std::string& filePath,
                                                 const REAL absTol)
{
   using namespace ObjCryst;
   Crystal c = MakePbso4Crystal();
   DiffractionDataSingleCrystal sc(c, false);
   sc.SetRadiationType(radiation);
   sc.SetWavelength(wavelength);
   sc.GenHKLFullSpace2(0.32, true);
   CompareSingleCrystalToGroundTruth(sc, filePath, absTol);
}

void ComparePowderSimulationToGroundTruth(const ObjCryst::RadiationType radiation,
                                          const REAL wavelength,
                                          ObjCryst::ReflectionProfile* profile,
                                          const std::string& filePath,
                                          const REAL absTol,
                                          const REAL relTol)
{
   using namespace ObjCryst;
   Crystal c = MakePbso4Crystal();

   PowderPattern p;
   p.SetRadiationType(radiation);
   p.SetWavelength(wavelength);
   p.SetPowderPatternPar(10 * DEG2RAD, 0.05 * DEG2RAD, 240);
   CrystVector_REAL obs(240);
   obs = 1;
   p.SetPowderPatternObs(obs);

   auto* phase = new PowderPatternDiffraction;
   phase->SetCrystal(c);
   phase->SetProfile(profile);
   p.AddPowderPatternComponent(*phase);
   p.Prepare();

   const std::vector<PowderPoint> gt = LoadPowderGroundTruth(filePath);
   const auto& calc = p.GetPowderPatternCalc();
   const unsigned long nbPoint = p.GetNbPoint();
   for(const auto& pt : gt)
   {
      Check(pt.idx < nbPoint, "Powder ground-truth index out of range");
      CheckNearAbsRel(calc(pt.idx), pt.intensity, absTol, relTol, "Powder intensity mismatch against ground truth");
   }
}

void DumpGroundTruthData()
{
   using namespace ObjCryst;
   std::cout << std::setprecision(9) << std::fixed;

   Crystal cx = MakePbso4Crystal();
   DiffractionDataSingleCrystal scX(cx, false);
   scX.SetRadiationType(RAD_XRAY);
   scX.SetWavelength(1.54056f);
   scX.GenHKLFullSpace2(0.32, true);
   std::cout << "singlecrystal-xray\n";
   for(int i = 0; i < 20; ++i)
   {
      std::cout << static_cast<long>(scX.GetH()(i)) << " "
                << static_cast<long>(scX.GetK()(i)) << " "
                << static_cast<long>(scX.GetL()(i)) << " "
                << scX.GetFhklCalcSq()(i) << "\n";
   }

   Crystal cn = MakePbso4Crystal();
   DiffractionDataSingleCrystal scN(cn, false);
   scN.SetRadiationType(RAD_NEUTRON);
   scN.SetWavelength(1.54056f);
   scN.GenHKLFullSpace2(0.32, true);
   std::cout << "singlecrystal-neutron\n";
   for(int i = 0; i < 20; ++i)
   {
      std::cout << static_cast<long>(scN.GetH()(i)) << " "
                << static_cast<long>(scN.GetK()(i)) << " "
                << static_cast<long>(scN.GetL()(i)) << " "
                << scN.GetFhklCalcSq()(i) << "\n";
   }

   Crystal cXg = MakePbso4Crystal();
   PowderPattern pXg;
   pXg.SetRadiationType(RAD_XRAY);
   pXg.SetWavelength(1.54056f);
   pXg.SetPowderPatternPar(10 * DEG2RAD, 0.05 * DEG2RAD, 240);
   CrystVector_REAL obsXg(240); obsXg = 1; pXg.SetPowderPatternObs(obsXg);
   auto* phaseXg = new PowderPatternDiffraction;
   phaseXg->SetCrystal(cXg);
   auto* pvGx = new ReflectionProfilePseudoVoigt();
   pvGx->SetProfilePar(.03f * DEG2RAD * DEG2RAD, 0, 0, 0.0f, 0);
   phaseXg->SetProfile(pvGx);
   pXg.AddPowderPatternComponent(*phaseXg);
   pXg.Prepare();
   std::cout << "powder-xray-pseudo-voigt-gaussian\n";
   for(unsigned long i = 0; i < pXg.GetNbPoint(); i += 12) std::cout << i << " " << pXg.GetPowderPatternCalc()(i) << "\n";

   Crystal cXl = MakePbso4Crystal();
   PowderPattern pXl;
   pXl.SetRadiationType(RAD_XRAY);
   pXl.SetWavelength(1.54056f);
   pXl.SetPowderPatternPar(10 * DEG2RAD, 0.05 * DEG2RAD, 240);
   CrystVector_REAL obsXl(240); obsXl = 1; pXl.SetPowderPatternObs(obsXl);
   auto* phaseXl = new PowderPatternDiffraction;
   phaseXl->SetCrystal(cXl);
   auto* pvLx = new ReflectionProfilePseudoVoigt();
   pvLx->SetProfilePar(.03f * DEG2RAD * DEG2RAD, 0, 0, 1.0f, 0);
   phaseXl->SetProfile(pvLx);
   pXl.AddPowderPatternComponent(*phaseXl);
   pXl.Prepare();
   std::cout << "powder-xray-pseudo-voigt-lorentzian\n";
   for(unsigned long i = 0; i < pXl.GetNbPoint(); i += 12) std::cout << i << " " << pXl.GetPowderPatternCalc()(i) << "\n";

   Crystal cXa = MakePbso4Crystal();
   PowderPattern pXa;
   pXa.SetRadiationType(RAD_XRAY);
   pXa.SetWavelength(1.54056f);
   pXa.SetPowderPatternPar(10 * DEG2RAD, 0.05 * DEG2RAD, 240);
   CrystVector_REAL obsXa(240); obsXa = 1; pXa.SetPowderPatternObs(obsXa);
   auto* phaseXa = new PowderPatternDiffraction;
   phaseXa->SetCrystal(cXa);
   auto* anisoX = new ReflectionProfilePseudoVoigtAnisotropic();
   anisoX->SetProfilePar(.02f * DEG2RAD * DEG2RAD, .003f * DEG2RAD * DEG2RAD, 0, 0,
                         0.002f * DEG2RAD, 0.002f * DEG2RAD, 0, 0, 0, 0, 0, 0, 0.4f, 0, 0, 0, 0);
   phaseXa->SetProfile(anisoX);
   pXa.AddPowderPatternComponent(*phaseXa);
   pXa.Prepare();
   std::cout << "powder-xray-anisotropic\n";
   for(unsigned long i = 0; i < pXa.GetNbPoint(); i += 12) std::cout << i << " " << pXa.GetPowderPatternCalc()(i) << "\n";

   Crystal cNg = MakePbso4Crystal();
   PowderPattern pNg;
   pNg.SetRadiationType(RAD_NEUTRON);
   pNg.SetWavelength(1.54056f);
   pNg.SetPowderPatternPar(10 * DEG2RAD, 0.05 * DEG2RAD, 240);
   CrystVector_REAL obsNg(240); obsNg = 1; pNg.SetPowderPatternObs(obsNg);
   auto* phaseNg = new PowderPatternDiffraction;
   phaseNg->SetCrystal(cNg);
   auto* pvGn = new ReflectionProfilePseudoVoigt();
   pvGn->SetProfilePar(.03f * DEG2RAD * DEG2RAD, 0, 0, 0.0f, 0);
   phaseNg->SetProfile(pvGn);
   pNg.AddPowderPatternComponent(*phaseNg);
   pNg.Prepare();
   std::cout << "powder-neutron-pseudo-voigt-gaussian\n";
   for(unsigned long i = 0; i < pNg.GetNbPoint(); i += 12) std::cout << i << " " << pNg.GetPowderPatternCalc()(i) << "\n";

   Crystal cNl = MakePbso4Crystal();
   PowderPattern pNl;
   pNl.SetRadiationType(RAD_NEUTRON);
   pNl.SetWavelength(1.54056f);
   pNl.SetPowderPatternPar(10 * DEG2RAD, 0.05 * DEG2RAD, 240);
   CrystVector_REAL obsNl(240); obsNl = 1; pNl.SetPowderPatternObs(obsNl);
   auto* phaseNl = new PowderPatternDiffraction;
   phaseNl->SetCrystal(cNl);
   auto* pvLn = new ReflectionProfilePseudoVoigt();
   pvLn->SetProfilePar(.03f * DEG2RAD * DEG2RAD, 0, 0, 1.0f, 0);
   phaseNl->SetProfile(pvLn);
   pNl.AddPowderPatternComponent(*phaseNl);
   pNl.Prepare();
   std::cout << "powder-neutron-pseudo-voigt-lorentzian\n";
   for(unsigned long i = 0; i < pNl.GetNbPoint(); i += 12) std::cout << i << " " << pNl.GetPowderPatternCalc()(i) << "\n";
}

ObjCryst::Crystal MakePbso4Crystal()
{
   using namespace ObjCryst;
   Crystal c(8.516, 5.399, 6.989, "Pnma");
   auto* pb = new ScatteringPowerAtom("Pb", "Pb", 1.48);
   auto* s = new ScatteringPowerAtom("S", "S", 0.74);
   auto* o = new ScatteringPowerAtom("O", "O", 1.2);
   c.AddScatteringPower(pb);
   c.AddScatteringPower(s);
   c.AddScatteringPower(o);
   c.AddScatterer(new Atom(.1882, .250, .167, "Pb1", pb, 1.));
   c.AddScatterer(new Atom(.4370, .750, .186, "S1", s, 1.));
   c.AddScatterer(new Atom(.5950, .750, .100, "O1", o, 1.));
   c.AddScatterer(new Atom(.3190, .750, .043, "O2", o, 1.));
   c.AddScatterer(new Atom(.4150, .974, .306, "O3", o, 1.));
   return c;
}

void TestSpaceGroup()
{
   using namespace ObjCryst;
   SpaceGroup sg("Pnma");
   Check(sg.GetSpaceGroupNumber() == 62, "SpaceGroup number mismatch for Pnma");
   Check(sg.IsCentrosymmetric(), "Pnma should be centrosymmetric");
   Check(sg.GetNbTranslationVectors() == 1, "Pnma should have one translation vector");
   Check(sg.GetNbSymmetrics(false, false) >= 4, "Pnma should expose at least 4 symmetrics");
}

void TestScatteringPowerAtom()
{
   using namespace ObjCryst;
   ScatteringPowerAtom o("O", "O", 1.2);
   Check(o.GetAtomicNumber() == 8, "Atomic number mismatch for oxygen");
   Check(o.GetAtomicWeight() > 10.0, "Atomic weight for oxygen invalid");
   Check(o.GetRadius() > 0.1, "Atomic radius should be positive");
   Check(o.GetCovalentRadius() > 0.1, "Covalent radius should be positive");
   Check(o.GetForwardScatteringFactor(RAD_XRAY) > 0, "Forward scattering factor should be positive");
}

void TestCrystalAndAtom()
{
   using namespace ObjCryst;
   Crystal c = MakePbso4Crystal();
   Check(c.GetNbScatterer() == 5, "Crystal should contain 5 scatterers");
   Check(c.GetScatteringPowerRegistry().GetNb() == 3, "Crystal should contain 3 scattering powers");

   Atom& pb = dynamic_cast<Atom&>(c.GetScatt(0));
   pb.SetX(0.2);
   pb.SetY(0.3);
   pb.SetZ(0.4);
   Check(std::fabs(pb.GetX() - 0.2) < 1e-6, "Atom X coordinate update failed");
   Check(pb.GetMass() > 1., "Atom mass should be positive");
   Check(!pb.IsDummy(), "Pb atom should not be dummy");
}

void TestMolecule()
{
   using namespace ObjCryst;
   Crystal c(10, 10, 10, "P1");
   auto* spO = new ScatteringPowerAtom("O", "O", 1.0);
   c.AddScatteringPower(spO);
   Molecule m(c, "water-like");
   m.AddAtom(0.0, 0.0, 0.0, spO, "O1", false);
   m.AddAtom(0.1, 0.0, 0.0, spO, "O2", false);
   Check(m.GetAtomList().size() == 2, "Molecule atom insertion failed");
   m.AddBond(m.GetAtom(0), m.GetAtom(1), 1.0, 0.01, 0.02, 1.0, false);
   Check(m.GetBondList().size() == 1, "Molecule bond insertion failed");
   Check(m.FindBond(m.GetAtom(0), m.GetAtom(1)) != m.GetBondList().end(), "Molecule bond lookup failed");
}

void TestScatteringDataAndSingleCrystal()
{
   using namespace ObjCryst;
   Crystal c = MakePbso4Crystal();
   DiffractionDataSingleCrystal sc(c, false);
   sc.SetRadiationType(RAD_XRAY);
   sc.SetWavelength(1.54056);
   sc.GenHKLFullSpace2(0.35, true);
   Check(sc.GetNbRefl() > 50, "Single-crystal reflection generation failed");
   Check(sc.GetH().numElements() == sc.GetNbRefl(), "H array size mismatch");
   Check(sc.GetK().numElements() == sc.GetNbRefl(), "K array size mismatch");
   Check(sc.GetL().numElements() == sc.GetNbRefl(), "L array size mismatch");
}

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

void TestCifImport()
{
   using namespace ObjCryst;
   std::ifstream cifIn("../../test/data/cif/PbSO4-COD-1528837.cif");
   Check(cifIn.good(), "Unable to open COD CIF test file");
   CIF cif(cifIn, true, false);
   Crystal c;
   Crystal* result = CreateCrystalFromCIF(cif, false, true, false, false, &c);
   Check(result == &c, "CreateCrystalFromCIF returned wrong pointer");
   Check(c.GetSpaceGroup().GetSpaceGroupNumber() == 62, "CIF spacegroup mismatch");
   Check(c.GetNbScatterer() >= 5, "CIF import created too few scatterers");
}

void TestRefinablePar()
{
   using namespace ObjCryst;
   REAL v = 1.0;
   RefinablePar p("test-par", &v, 0.0, 2.0, gpRefParTypeObjCryst);
   p.SetValue(1.5);
   Check(std::fabs(p.GetValue() - 1.5) < 1e-12, "RefinablePar SetValue failed");
   p.Mutate(1.0);
   Check(p.GetValue() <= 2.0, "RefinablePar should respect upper limit");
   p.SetHumanScale(2.0);
   p.SetHumanValue(2.0);
   Check(std::fabs(p.GetValue() - 1.0) < 1e-12, "RefinablePar human scaling failed");
   p.SetIsFixed(true);
   Check(p.IsFixed(), "RefinablePar fixed flag not applied");
}

void TestRefinableObj()
{
   using namespace ObjCryst;
   RefinableObj obj;
   REAL a = 0.1, b = 0.2;
   obj.AddPar(new RefinablePar("a", &a, -1, 1, gpRefParTypeObjCryst));
   obj.AddPar(new RefinablePar("b", &b, -1, 1, gpRefParTypeObjCryst));
   Check(obj.GetNbPar() == 2, "RefinableObj parameter count mismatch");
   obj.SetParIsFixed("a", true);
   obj.PrepareForRefinement();
   Check(obj.GetNbParNotFixed() == 1, "RefinableObj fixed parameter handling failed");
   const unsigned long setId = obj.CreateParamSet("baseline");
   Check(obj.GetParamSetName(setId) == "baseline", "RefinableObj parameter-set naming failed");
   Check(obj.GetParamSet(setId).numElements() > 0, "RefinableObj parameter-set should not be empty");
}

void TestOptimizationObj()
{
   using namespace ObjCryst;
   RefinableObj obj;
   REAL a = 0.1;
   obj.AddPar(new RefinablePar("a", &a, -1, 1, gpRefParTypeObjCryst));
   MonteCarloObj mc;
   mc.AddRefinableObj(obj);
   RefinableObj& full = mc.GetFullRefinableObj(true);
   Check(full.GetNbPar() >= 1, "OptimizationObj failed to aggregate refined parameters");
   mc.FixAllPar();
   mc.UnFixAllPar();
}

void TestLsqNumObj()
{
   using namespace ObjCryst;
   PowderPattern p;
   p.SetRadiationType(RAD_XRAY);
   p.SetWavelength("CuA1");
   p.SetPowderPatternPar(5 * DEG2RAD, 0.05 * DEG2RAD, 1000);
   CrystVector_REAL obs(1000);
   obs = 10;
   p.SetPowderPatternObs(obs);
   auto* bg = new PowderPatternBackground;
   CrystVector_REAL x(2), y(2);
   x(0) = 0;
   x(1) = 3.0 * DEG2RAD;
   y(0) = 8;
   y(1) = 9;
   bg->SetInterpPoints(x, y);
   p.AddPowderPatternComponent(*bg);
   p.Prepare();

   LSQNumObj lsq("lsq-test");
   lsq.SetRefinedObj(p, 0, true, false);
   lsq.PrepareRefParList(true);
   Check(lsq.GetCompiledRefinedObj().GetNbPar() >= 1, "LSQNumObj compiled parameter list is empty");
}

void TestCellExplorer()
{
   using namespace ObjCryst;
   PeakList peaks;
   peaks.Simulate(0, 8.516f, 5.399f, 6.989f, 90.f, 90.f, 90.f, true, 15, 0, 0, 0, false);
   Check(peaks.GetPeakList().size() == 15, "PeakList simulation failed");
   const float estimatedVolume = EstimateCellVolume(peaks.GetPeakList().back().dobs,
                                                    peaks.GetPeakList().front().dobs / 10.f,
                                                    static_cast<float>(peaks.GetPeakList().size()),
                                                    ORTHOROMBIC, LATTICE_P, 1.f);
   Check(estimatedVolume > 0, "EstimateCellVolume failed");
   RecUnitCell cell(0.f, 1.f / 8.516f, 1.f / 5.399f, 1.f / 6.989f, 0.f, 0.f, 0.f,
                    ORTHOROMBIC, LATTICE_P, 0);
   const float score = Score(peaks, cell, 0, false, false, false);
   Check(score > 0, "CellExplorer score path failed");
}

void TestSingleCrystalGroundTruthXray()
{
   CompareSingleCrystalSimulationToGroundTruth(ObjCryst::RAD_XRAY, 1.54056f,
                                               "../../test/data/ground_truth/singlecrystal_xray_pbso4.txt", 1e-4f);
}

void TestSingleCrystalGroundTruthNeutron()
{
   CompareSingleCrystalSimulationToGroundTruth(ObjCryst::RAD_NEUTRON, 1.54056f,
                                               "../../test/data/ground_truth/singlecrystal_neutron_pbso4.txt", 1e-4f);
}

void TestPowderGroundTruthXrayPseudoVoigtGaussian()
{
   using namespace ObjCryst;
   auto* profile = new ReflectionProfilePseudoVoigt();
   profile->SetProfilePar(.03f * DEG2RAD * DEG2RAD, 0, 0, 0.0f, 0);
   ComparePowderSimulationToGroundTruth(RAD_XRAY, 1.54056f, profile,
                                        "../../test/data/ground_truth/powder_xray_pv_gaussian_pbso4.txt", 1e-2f, 1e-5f);
}

void TestPowderGroundTruthXrayPseudoVoigtLorentzian()
{
   using namespace ObjCryst;
   auto* profile = new ReflectionProfilePseudoVoigt();
   profile->SetProfilePar(.03f * DEG2RAD * DEG2RAD, 0, 0, 1.0f, 0);
   ComparePowderSimulationToGroundTruth(RAD_XRAY, 1.54056f, profile,
                                        "../../test/data/ground_truth/powder_xray_pv_lorentzian_pbso4.txt", 1e-2f, 1e-5f);
}

void TestPowderGroundTruthXrayAnisotropic()
{
   using namespace ObjCryst;
   auto* profile = new ReflectionProfilePseudoVoigtAnisotropic();
   profile->SetProfilePar(.02f * DEG2RAD * DEG2RAD, .003f * DEG2RAD * DEG2RAD, 0, 0,
                          0.002f * DEG2RAD, 0.002f * DEG2RAD, 0, 0, 0, 0, 0, 0, 0.4f, 0, 0, 0, 0);
   ComparePowderSimulationToGroundTruth(RAD_XRAY, 1.54056f, profile,
                                        "../../test/data/ground_truth/powder_xray_anisotropic_pbso4.txt", 1e-2f, 1e-5f);
}

void TestPowderGroundTruthNeutronPseudoVoigtGaussian()
{
   using namespace ObjCryst;
   auto* profile = new ReflectionProfilePseudoVoigt();
   profile->SetProfilePar(.03f * DEG2RAD * DEG2RAD, 0, 0, 0.0f, 0);
   ComparePowderSimulationToGroundTruth(RAD_NEUTRON, 1.54056f, profile,
                                        "../../test/data/ground_truth/powder_neutron_pv_gaussian_pbso4.txt", 1e-2f, 1e-5f);
}

void TestPowderGroundTruthNeutronPseudoVoigtLorentzian()
{
   using namespace ObjCryst;
   auto* profile = new ReflectionProfilePseudoVoigt();
   profile->SetProfilePar(.03f * DEG2RAD * DEG2RAD, 0, 0, 1.0f, 0);
   ComparePowderSimulationToGroundTruth(RAD_NEUTRON, 1.54056f, profile,
                                        "../../test/data/ground_truth/powder_neutron_pv_lorentzian_pbso4.txt", 1e-2f, 1e-5f);
}

} // namespace

int main(int argc, char* argv[])
{
   if(argc != 2)
   {
      std::cerr << "Usage: api_surface_tests <test-case>" << std::endl;
      return 2;
   }

   const std::string testName = argv[1];

   if(testName == "spacegroup") TestSpaceGroup();
   else if(testName == "scattering-power-atom") TestScatteringPowerAtom();
   else if(testName == "crystal-atom") TestCrystalAndAtom();
   else if(testName == "molecule") TestMolecule();
   else if(testName == "scatteringdata-singlecrystal") TestScatteringDataAndSingleCrystal();
   else if(testName == "powderpattern-background") TestPowderPatternBackground();
   else if(testName == "powderpattern-diffraction") TestPowderPatternDiffraction();
   else if(testName == "powderpattern-import") TestPowderPatternImport();
   else if(testName == "cif-import") TestCifImport();
   else if(testName == "refinablepar") TestRefinablePar();
   else if(testName == "refinableobj") TestRefinableObj();
   else if(testName == "optimizationobj") TestOptimizationObj();
   else if(testName == "lsqnumobj") TestLsqNumObj();
   else if(testName == "cellexplorer") TestCellExplorer();
   else if(testName == "singlecrystal-groundtruth-xray") TestSingleCrystalGroundTruthXray();
   else if(testName == "singlecrystal-groundtruth-neutron") TestSingleCrystalGroundTruthNeutron();
   else if(testName == "powder-groundtruth-xray-pv-gaussian") TestPowderGroundTruthXrayPseudoVoigtGaussian();
   else if(testName == "powder-groundtruth-xray-pv-lorentzian") TestPowderGroundTruthXrayPseudoVoigtLorentzian();
   else if(testName == "powder-groundtruth-xray-anisotropic") TestPowderGroundTruthXrayAnisotropic();
   else if(testName == "powder-groundtruth-neutron-pv-gaussian") TestPowderGroundTruthNeutronPseudoVoigtGaussian();
   else if(testName == "powder-groundtruth-neutron-pv-lorentzian") TestPowderGroundTruthNeutronPseudoVoigtLorentzian();
   else if(testName == "dump-ground-truth") DumpGroundTruthData();
   else
   {
      std::cerr << "Unknown test case: " << testName << std::endl;
      return 2;
   }

   return 0;
}
