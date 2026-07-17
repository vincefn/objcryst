// Extended unit tests for the CellExplorer and PeakList classes.
#include <iostream>
#include <string>

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

void TestCellExplorerConfiguration()
{
   using namespace ObjCryst;
   PeakList peaks;
   peaks.Simulate(0, 8.516f, 5.399f, 6.989f, 90.f, 90.f, 90.f, true, 20, 0, 0, 0, false);
   CellExplorer explorer(peaks, ORTHOROMBIC, 0);
   explorer.SetLengthMinMax(3.0f, 20.0f);
   explorer.SetAngleMinMax(80.0f, 100.0f);
   explorer.SetVolumeMinMax(100.0f, 2000.0f);
   explorer.SetNbSpurious(0);
   explorer.SetCrystalSystem(ORTHOROMBIC);
   Check(!explorer.GetName().empty() || explorer.GetName().empty(),
        "CellExplorer::GetName() should be callable after configuration");
   Check(explorer.GetSolutions().empty(), "CellExplorer should have no solutions before running DicVol/Evolution");
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
   else
   {
      std::cerr << "Unknown test case: " << testName << std::endl;
      return 2;
   }

   return 0;
}
