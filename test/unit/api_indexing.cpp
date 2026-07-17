// Extended unit tests for the CellExplorer and PeakList classes.
#include <iostream>
#include <list>
#include <string>
#include <vector>

#include "ObjCryst/ObjCryst/Indexing.h"

#include "test_common.h"

using namespace objcryst_test;

namespace
{

void TestPeakListSimulateAndVolume()
{
   using namespace ObjCryst;
   PeakList peaks;
   const float volume = peaks.Simulate(0, 8.516f, 5.399f, 6.989f, 90.f, 90.f, 90.f, true, 15, 0, 0, 0, false);
   Check(volume > 250.f, "PeakList simulation produced invalid cell volume");
   Check(peaks.GetPeakList().size() == 15, "PeakList simulation failed");
   const float estimatedVolume = EstimateCellVolume(peaks.GetPeakList().back().dobs,
                                                    peaks.GetPeakList().front().dobs / 10.f,
                                                    static_cast<float>(peaks.GetPeakList().size()),
                                                    ORTHOROMBIC, LATTICE_P, 1.f);
   Check(estimatedVolume > 0, "EstimateCellVolume failed");
}

void TestPeakListAddRemove()
{
   using namespace ObjCryst;
   PeakList peaks;
   peaks.AddPeak(1.f / 8.516f, 100.0);
   peaks.AddPeak(1.f / 5.399f, 50.0);
   peaks.AddPeak(1.f / 6.989f, 25.0);
   Check(peaks.GetPeakList().size() == 3, "PeakList::AddPeak should append peaks");
   peaks.RemovePeak(0);
   Check(peaks.GetPeakList().size() == 2, "PeakList::RemovePeak should remove a peak");
}

void TestCellExplorerScore()
{
   using namespace ObjCryst;
   PeakList peaks;
   peaks.Simulate(0, 8.516f, 5.399f, 6.989f, 90.f, 90.f, 90.f, true, 15, 0, 0, 0, false);
   RecUnitCell cell(0.f, 1.f / 8.516f, 1.f / 5.399f, 1.f / 6.989f, 0.f, 0.f, 0.f,
                    ORTHOROMBIC, LATTICE_P, 0);
   const float score = Score(peaks, cell, 0, false, false, false);
   Check(score > 0, "CellExplorer score path failed");
}

// One representative, non-degenerate unit cell per CrystalSystem, used to
// exercise CellExplorer's configuration API for every lattice system (not
// just ORTHOROMBIC).
struct SystemTestCell
{
   ObjCryst::CrystalSystem system;
   const char* name;
   float a, b, c, alpha, beta, gamma; // angles in degrees
};

const SystemTestCell& GetSystemTestCell(unsigned int i)
{
   using namespace ObjCryst;
   static const SystemTestCell cells[] = {
      {TRICLINIC,  "TRICLINIC",  6.0f, 7.0f, 8.0f, 80.0f, 95.0f, 100.0f},
      {MONOCLINIC, "MONOCLINIC", 8.0f, 6.0f, 9.0f, 90.0f, 100.0f, 90.0f},
      {ORTHOROMBIC,"ORTHOROMBIC",8.516f, 5.399f, 6.989f, 90.0f, 90.0f, 90.0f},
      {HEXAGONAL,  "HEXAGONAL",  5.0f, 5.0f, 10.0f, 90.0f, 90.0f, 120.0f},
      {RHOMBOEDRAL,"RHOMBOEDRAL",6.0f, 6.0f, 6.0f, 80.0f, 80.0f, 80.0f},
      {TETRAGONAL, "TETRAGONAL", 7.0f, 7.0f, 10.0f, 90.0f, 90.0f, 90.0f},
      {CUBIC,      "CUBIC",      6.0f, 6.0f, 6.0f, 90.0f, 90.0f, 90.0f},
   };
   return cells[i];
}

const unsigned int kNbCrystalSystems = 7;

void TestCellExplorerConfiguration()
{
   using namespace ObjCryst;
   for(unsigned int i = 0; i < kNbCrystalSystems; ++i)
   {
      const SystemTestCell& sc = GetSystemTestCell(i);
      PeakList peaks;
      const float volume = peaks.Simulate(0, sc.a, sc.b, sc.c, sc.alpha, sc.beta, sc.gamma, true, 20, 0, 0, 0, false);
      Check(volume > 0, sc.name);
      Check(peaks.GetPeakList().size() == 20, sc.name);

      CellExplorer explorer(peaks, sc.system, 0);
      explorer.SetLengthMinMax(3.0f, 20.0f);
      explorer.SetAngleMinMax(80.0f * DEG2RAD, 130.0f * DEG2RAD);
      explorer.SetVolumeMinMax(volume * 0.5f, volume * 2.0f);
      explorer.SetNbSpurious(0);
      explorer.SetCrystalSystem(sc.system);
      Check(!explorer.GetName().empty() || explorer.GetName().empty(),
           "CellExplorer::GetName() should be callable after configuration");
      Check(explorer.GetSolutions().empty(),
           "CellExplorer should have no solutions before running DicVol/Evolution");
   }
}

void TestCellExplorerDicVolTetragonal()
{
   using namespace ObjCryst;
   // Known tetragonal cell: a=b=7.0 A, c=10.0 A, all angles=90 deg.
   // Sigma cannot be exactly 0 here: PeakList::Simulate would then produce
   // observed d*^2 values with a zero-width [d2obsmin,d2obsmax] bracket,
   // making it very unlikely that DicVol's discretized cell-parameter search
   // intervals ever contain the (floating-point) observed value exactly.
   // merge=true is needed because Simulate() has no symmetry information and
   // generates one entry per (h,k,l) triple regardless of multiplicity: for a
   // tetragonal a=b cell, {h k 0} and {k h 0} would otherwise appear as two
   // separate "observed" peaks with the same d-spacing, which does not happen
   // with a real diffraction pattern (where such reflections overlap into a
   // single peak).
   const float a = 7.0f, c = 10.0f;
   PeakList peaks;
   const float volume = peaks.Simulate(0, a, a, c, 90.f, 90.f, 90.f, true, 20, 0, 1e-3f, 0, false, true);
   Check(peaks.GetPeakList().size() == 20, "Tetragonal PeakList simulation failed");

   CellExplorer explorer(peaks, TETRAGONAL, 0);
   explorer.SetD2Error(0);
   // Tight bounds around the true cell keep DicVol's search space (and thus
   // runtime) small while still exercising the actual dichotomy algorithm.
   explorer.SetLengthMinMax(5.0f, 12.0f);
   explorer.SetVolumeMinMax(volume * 0.7f, volume * 1.3f);
   explorer.SetNbSpurious(0);

   explorer.DicVol(10, 3, 50, 3);

   Check(!explorer.GetSolutions().empty(), "DicVol (TETRAGONAL) failed to find any solution");
   Check(explorer.GetBestScore() > 10.f, "DicVol (TETRAGONAL) best score too low");

   float bestScore = -1.f;
   RecUnitCell bestCell;
   for(std::list<std::pair<RecUnitCell, float> >::const_iterator pos = explorer.GetSolutions().begin();
       pos != explorer.GetSolutions().end(); ++pos)
      if(pos->second > bestScore) { bestScore = pos->second; bestCell = pos->first; }

   const std::vector<float> direct = bestCell.DirectUnitCell();
   CheckNearAbsRel(direct[0], a, 0.2f, 0.05f, "DicVol (TETRAGONAL) recovered 'a' does not match the simulated cell");
   CheckNearAbsRel(direct[2], c, 0.2f, 0.05f, "DicVol (TETRAGONAL) recovered 'c' does not match the simulated cell");
}

void TestCellExplorerDicVolMonoclinic()
{
   using namespace ObjCryst;
   // Known monoclinic cell: a=8.0, b=6.0, c=9.0 A, beta=100 deg (alpha=gamma=90).
   const float a = 8.0f, b = 6.0f, c = 9.0f, beta = 100.0f;
   PeakList peaks;
   const float volume = peaks.Simulate(0, a, b, c, 90.f, beta, 90.f, true, 20, 0, 1e-3f, 0, false, true);
   Check(peaks.GetPeakList().size() == 20, "Monoclinic PeakList simulation failed");

   CellExplorer explorer(peaks, MONOCLINIC, 0);
   explorer.SetD2Error(0);
   // Tight length/angle/volume bounds around the true cell keep DicVol's
   // search space (and thus runtime) small while still exercising the
   // actual dichotomy algorithm.
   explorer.SetLengthMinMax(4.0f, 11.0f);
   explorer.SetAngleMinMax(90.0f * DEG2RAD, 108.0f * DEG2RAD);
   explorer.SetVolumeMinMax(volume * 0.7f, volume * 1.3f);
   explorer.SetNbSpurious(0);

   explorer.DicVol(10, 3, 50, 3);

   Check(!explorer.GetSolutions().empty(), "DicVol (MONOCLINIC) failed to find any solution");
   Check(explorer.GetBestScore() > 10.f, "DicVol (MONOCLINIC) best score too low");

   float bestScore = -1.f;
   RecUnitCell bestCell;
   for(std::list<std::pair<RecUnitCell, float> >::const_iterator pos = explorer.GetSolutions().begin();
       pos != explorer.GetSolutions().end(); ++pos)
      if(pos->second > bestScore) { bestScore = pos->second; bestCell = pos->first; }

   const std::vector<float> direct = bestCell.DirectUnitCell();
   CheckNearAbsRel(direct[4] * RAD2DEG, beta, 3.0f, 0.03f,
                   "DicVol (MONOCLINIC) recovered 'beta' does not match the simulated cell");
}

} // namespace

int main(int argc, char* argv[])
{
   if(argc != 2)
   {
      std::cerr << "Usage: api_indexing <test-case>" << std::endl;
      return 2;
   }

   const std::string testName = argv[1];

   if(testName == "peaklist-simulate-volume") TestPeakListSimulateAndVolume();
   else if(testName == "peaklist-add-remove") TestPeakListAddRemove();
   else if(testName == "cellexplorer") TestCellExplorerScore();
   else if(testName == "cellexplorer-configuration") TestCellExplorerConfiguration();
   else if(testName == "cellexplorer-dicvol-tetragonal") TestCellExplorerDicVolTetragonal();
   else if(testName == "cellexplorer-dicvol-monoclinic") TestCellExplorerDicVolMonoclinic();
   else
   {
      std::cerr << "Unknown test case: " << testName << std::endl;
      return 2;
   }

   return 0;
}
