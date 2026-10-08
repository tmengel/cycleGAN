// Pass 1: GEANT4 simulation of JEWEL HepMC events -> G4Hits DST.
//
// Reads <nEvents> events, after skipping <skip>, from one JEWEL HepMC2 file, runs them through
// the sPHENIX calorimeters (EMCal, iHCal, oHCal) with magnet, beam pipe, MBD/EPD and plug door,
// and writes a DST with the G4Hits and G4TruthInfo. Nothing is reconstructed here: towers,
// vertex and jets are all built in pass 2, so pass 2 can be rerun with different calorimeter
// settings without repeating GEANT4.
//
// Adapted from cycleGAN/code/pass1/Fun4All_G4_JEWEL_pass1.C (itself derived from the MDC2
// Fun4All_G4_Pass1.C). Differences, all deliberate:
//   * Beam/vertex: that macro took the vertex settings from RunSettings(runnumber), which maps
//     every run number >= 100 to Input::mRad_00 -- the 2024 pp zero-crossing-angle beam with a
//     sigma_z = 65 cm vertex (its trees show |vz| up to 125 cm). Here the configuration is an
//     explicit argument. Default "fixed": every event at vertex (0, 0, 0), head-on beams, no
//     smearing. Also available: "AuAu" (sigma_z 13.5 cm, 1 mrad crossing, 2024 beam spot),
//     "pp" (16 cm, 1.5 mrad) and "pp_zeroangle" (the old 65 cm behaviour).
//   * No truth-jet finding. The old macro clustered truth jets with embedding flag 1, but HepMC
//     input particles carry embedding id 0, so those jet containers were always empty. Truth
//     jets are made in pass 2 from G4TruthInfo with the right flag.
//   * No GlobalVertexReco and no Production_* output handling: pass 2 builds the vertex, and
//     the DST is written straight to the working directory (the run script moves it).
//   * A fixed random seed (argument), so a job can be rerun exactly.
//   * No tracking detectors (as before): the calorimeter images never use them, and they are
//     most of the GEANT4 CPU time.

#ifndef MACRO_FUN4ALL_G4_PASS1_C
#define MACRO_FUN4ALL_G4_PASS1_C

#include <GlobalVariables.C>

#include <G4Setup_sPHENIX.C>
#include <G4_Input.C>
#include <G4_RunSettings.C>

#include <ffamodules/CDBInterface.h>
#include <ffamodules/FlagHandler.h>
#include <ffamodules/HeadReco.h>
#include <ffamodules/SyncReco.h>

#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllSyncManager.h>

#include <phhepmc/PHHepMCGenHelper.h>

#include <phool/PHRandomSeed.h>
#include <phool/recoConsts.h>

#include <string>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)

int Fun4All_G4_pass1(
    const int nEvents = 2,
    const int skip = 0,
    const std::string &inputFile = "/sphenix/tg/tg01/jets/tmengel/JEWEL_hepmc/pthat16/vacuum/raw/vacuum_000001.hepmc",
    const std::string &outputFile = "G4Hits_vacuum-0000001000-000000.root",
    const int runnumber = 1000,
    const int segment = 0,
    const unsigned int randomSeed = 12345,
    const std::string &beamConfig = "fixed",
    const std::string &cdbtag = "MDC2_ana.435")
{
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);

  recoConsts *rc = recoConsts::instance();
  // Every random number in the job (vertex smearing, GEANT4) derives from this seed.
  rc->set_IntFlag("RANDOMSEED", static_cast<int>(randomSeed));
  PHRandomSeed::Verbosity(1);

  rc->set_IntFlag("RUNNUMBER", runnumber);
  se->getSyncManager()->SegmentNumber(segment);

  // RunSettings() for a generic run number (>= 100) only sets the beam configuration, to
  // mRad_00; it is overridden right below ("fixed" does not use it at all).
  RunSettings(runnumber);
  const bool fixedVertex = (beamConfig == "fixed");
  if (fixedVertex)
  {
    // vertex fixed at (0, 0, 0), t = 0, head-on beams: set after InputInit() below
  }
  else if (beamConfig == "AuAu")
  {
    Input::BEAM_CONFIGURATION = Input::AuAu_COLLISION;  // sigma_z 13.5 cm, 1 mrad crossing
  }
  else if (beamConfig == "pp")
  {
    Input::BEAM_CONFIGURATION = Input::pp_COLLISION;  // sigma_z 16 cm, 1.5 mrad crossing
  }
  else if (beamConfig == "pp_zeroangle")
  {
    Input::BEAM_CONFIGURATION = Input::pp_ZEROANGLE;  // sigma_z 65 cm (old pass1 behaviour)
  }
  else
  {
    std::cout << "Fun4All_G4_pass1: unknown beamConfig " << beamConfig
              << " (use fixed, AuAu, pp or pp_zeroangle)" << std::endl;
    gSystem->Exit(1);
  }

  // conditions DB (field map etc.)
  Enable::CDB = true;
  rc->set_StringFlag("CDB_GLOBALTAG", cdbtag);
  rc->set_uint64Flag("TIMESTAMP", runnumber);

  //===============
  // Input
  //===============
  Input::VERBOSITY = 0;
  Input::HEPMC = true;
  INPUTHEPMC::filename = inputFile;
  // JEWEL events are single nucleon-nucleon collisions: no flow afterburner, no HIJING
  // vertex flip, no Fermi motion.
  InputInit();
  if (fixedVertex)
  {
    // JEWEL writes every vertex at the origin; no smearing, no shift, no crossing angle
    PHHepMCGenHelper *gen = INPUTMANAGER::HepMCInputManager;
    gen->set_vertex_distribution_function(PHHepMCGenHelper::Gaus, PHHepMCGenHelper::Gaus,
                                          PHHepMCGenHelper::Gaus, PHHepMCGenHelper::Gaus);
    gen->set_vertex_distribution_mean(0, 0, 0, 0);
    gen->set_vertex_distribution_width(0, 0, 0, 0);
    gen->set_beam_direction_theta_phi(0, 0, M_PI, 0);
  }
  else
  {
    Input::ApplysPHENIXBeamParameter(INPUTMANAGER::HepMCInputManager);
  }
  InputRegister();

  se->registerSubsystem(new SyncReco());
  se->registerSubsystem(new FlagHandler());
  se->registerSubsystem(new HeadReco());

  //===============
  // Detectors
  //===============
  Enable::MBD = true;  // keeps MBD hits on the DST for a later MBD vertex/centrality
  Enable::EPD = true;
  Enable::PIPE = true;
  Enable::CEMC = true;
  Enable::HCALIN = true;
  Enable::MAGNET = true;
  Enable::HCALOUT = true;
  Enable::PLUGDOOR = true;
  Enable::BLACKHOLE = true;
  Enable::BLACKHOLE_FORWARD_SAVEHITS = false;

  G4Init();
  G4Setup();

  //===============
  // Output
  //===============
  InputManagers();
  Fun4AllDstOutputManager *out = new Fun4AllDstOutputManager("DSTOUT", outputFile);
  se->registerOutputManager(out);

  se->skip(skip);
  se->run(nEvents);

  CDBInterface::instance()->Print();
  se->End();
  se->PrintTimer();
  std::cout << "All done" << std::endl;
  delete se;
  gSystem->Exit(0);
  return 0;
}
#endif
