#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>

#include <ffamodules/CDBInterface.h>
#include <g4eval/SvtxEvaluator.h>
#include <G4_ActsGeom.C>
#include <G4Setup_sPHENIX.C>
#include <phool/recoConsts.h>
#include <trackingdiagnostics/TrackResiduals.h>

#include <simqa_modules/QAG4SimulationTracking.h>
#include <qautils/QAHistManagerDef.h>

#pragma GCC diagnostic push

#pragma GCC diagnostic ignored "-Wundefined-internal"

#include <secondaryvertextagger/SecondaryVertexTagger.h>

#pragma GCC diagnostic pop

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libsecondaryvertextagger.so)
R__LOAD_LIBRARY(libpdbcalBase.so)
R__LOAD_LIBRARY(libglobalvertex.so)
R__LOAD_LIBRARY(libphfield.so)
R__LOAD_LIBRARY(libg4detectors.so)
R__LOAD_LIBRARY(libg4eval.so)
R__LOAD_LIBRARY(libg4intt.so)
R__LOAD_LIBRARY(libphgeom.so)
R__LOAD_LIBRARY(libmvtx.so)
R__LOAD_LIBRARY(libmicromegas_io.so)
R__LOAD_LIBRARY(libsimqa_modules.so)
R__LOAD_LIBRARY(libTrackingDiagnostics.so)

void Fun4All_myTester(const std::string myInputFile = "some.file")
{
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(1);
  
  ACTSGEOM::ActsGeomInit();

  Enable::CDB = true;
  recoConsts *rc = recoConsts::instance();
  rc->set_StringFlag("CDB_GLOBALTAG", "ProdA_2024");
  rc->set_uint64Flag("TIMESTAMP", 6);

  Fun4AllInputManager *infile = new Fun4AllDstInputManager("DSTin");
  infile->AddFile(myInputFile);
  //infile->AddListFile("myList.txt");
  se->registerInputManager(infile);

  bool runTruth = true;
  std::string outputName = "./SiliconOnly_pythia8_events_";
  if (runTruth) outputName += "truth";
  else outputName += "reconstructed";
  outputName+= "Info_";

  string fileNumber = myInputFile;
  size_t findLastDash = fileNumber.find_last_of("_");
  if (findLastDash != string::npos) fileNumber.erase(0, findLastDash + 1);
  //string remove_this = ".root";
  //size_t pos = fileNumber.find(remove_this);
  //if (pos != string::npos) fileNumber.erase(pos, remove_this.length());
  outputName += fileNumber;
  
  //outputName+= ".root";

  SecondaryVertexTagger* myTagger = new SecondaryVertexTagger();
  myTagger->Verbosity(2);
  if (runTruth)
  {
    myTagger->truthMode();
  }
  else
  {
    myTagger->truthRecoMatch();
    myTagger->buildCharmDecays();
  }
  myTagger->setOutputName(outputName.c_str());
  myTagger->saveOutput();
  se->registerSubsystem(myTagger);

  se->run(0);

  se->End();
}
