// Pass 2: G4Hits DST -> flat trees (towertree / jettree / globaltree).
//
// Builds calorimeter cells and towers from the pass 1 G4Hits (digitization + calibration),
// retowers the EMCal onto the HCal grid, reconstructs the event vertex and anti-kT R=0.4 tower
// and truth jets, and writes everything with JewelTreeWriter (src/jeweltreewriter; its header
// documents the branches). The per-event TH2 images are made from these trees by
// make_images.C, run right after this macro by scripts/run_pass2.sh.
//
// Adapted from cycleGAN/code/trees/Fun4All_JEWEL_MakeTrees.C. Differences:
//   * GlobalVertexReco runs here (pass 1 no longer does). With no MBD/tracking vertex on the
//     DST it falls back to the primary truth vertex from G4TruthInfo, so the stored vertex and
//     the vertex the tower jets are corrected with are the true one -- the same behaviour the
//     old pass1+trees chain had.
//   * Fixed random seed (argument): the tower digitizers draw photon statistics.
//   * The calorimeter macros are the local snapshot in macros/calo/ (see below), not included
//     from another user's directory.
//
// Calorimeter macros: yeonjugo's G4_CEmc_Spacal.C / G4_HcalIn_ref.C / G4_HcalOut_ref.C, copied
// unchanged (md5 22b241c5..., fc35a5b5..., 681e5702..., identical to the versions chosen on
// 2026-09-15). No noise (pedestal width 0 in all three calorimeters) and a fixed EMCal
// calibration from the TowerCalibCombinedParams_2020 XML (Enable::CDB is false). The build's
// central <G4_*.C> macros give +43% tower energy plus noise -- do not swap them in without a
// decision to change the images.

#pragma once
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 00, 0)

#include <calo/G4_CEmc_Spacal.C>
#include <calo/G4_HcalIn_ref.C>
#include <calo/G4_HcalOut_ref.C>

#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <globalvertex/GlobalVertexReco.h>

#include <jetbase/FastJetAlgo.h>
#include <jetbase/Jet.h>
#include <jetbase/JetReco.h>
#include <jetbase/TowerJetInput.h>
#include <g4jets/TruthJetInput.h>

#include <jetbackground/RetowerCEMC.h>

#include <phool/PHRandomSeed.h>
#include <phool/recoConsts.h>

#include <jeweltreewriter/JewelTreeWriter.h>

#include <string>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4jets.so)
R__LOAD_LIBRARY(libjetbase.so)
R__LOAD_LIBRARY(libjetbackground.so)
R__LOAD_LIBRARY(libg4dst.so)
R__LOAD_LIBRARY(libglobalvertex.so)
R__LOAD_LIBRARY(libglobalvertex_io.so)
R__LOAD_LIBRARY(libjeweltreewriter.so)

#endif

void Fun4All_pass2(
    const std::string &inputDST = "G4Hits_vacuum-0000001000-000000.root",
    const std::string &outputTrees = "trees_vacuum-0000001000-000000.root",
    const unsigned int randomSeed = 54321)
{
  Fun4AllServer *se = Fun4AllServer::instance();
  const int verbosity = 0;
  se->Verbosity(verbosity);

  recoConsts *rc = recoConsts::instance();
  rc->set_IntFlag("RANDOMSEED", static_cast<int>(randomSeed));
  PHRandomSeed::Verbosity(1);

  // vertex first: TowerJetInput requires GlobalVertexMap
  GlobalVertexReco *gvertex = new GlobalVertexReco();
  gvertex->Verbosity(verbosity);
  se->registerSubsystem(gvertex);

  // cells -> towers -> digitization -> calibration (yeonjugo's macros)
  CEMC_Cells();
  HCALInner_Cells();
  HCALOuter_Cells();
  CEMC_Towers();
  HCALInner_Towers();
  HCALOuter_Towers();

  RetowerCEMC *rcemc = new RetowerCEMC();
  rcemc->Verbosity(verbosity);
  rcemc->set_towerinfo(true);
  se->registerSubsystem(rcemc);

  JetReco *towerjetreco = new JetReco("TOWERJETRECO");
  towerjetreco->add_input(new TowerJetInput(Jet::HCALIN_TOWERINFO));
  towerjetreco->add_input(new TowerJetInput(Jet::HCALOUT_TOWERINFO));
  towerjetreco->add_input(new TowerJetInput(Jet::CEMC_TOWERINFO_RETOWER));
  towerjetreco->add_algo(new FastJetAlgo(Jet::ANTIKT, 0.4, verbosity), "AntiKt_Tower_r04");
  towerjetreco->set_algo_node("ANTIKT");
  towerjetreco->set_input_node("TOWER");
  towerjetreco->Verbosity(verbosity);
  se->registerSubsystem(towerjetreco);

  // Truth jets from all final-state primaries (HepMC input particles carry embedding id 0).
  // The node name is the one JewelTreeWriter reads by default; it also keeps this macro
  // runnable on the old pass1 DSTs, which carry an (empty) AntiKt_Truth_r04 of their own.
  JetReco *truthjetreco = new JetReco("TRUTHJETRECO");
  TruthJetInput *tji = new TruthJetInput(Jet::PARTICLE);
  tji->add_embedding_flag(0);
  truthjetreco->add_input(tji);
  truthjetreco->add_algo(new FastJetAlgo(Jet::ANTIKT, 0.4, verbosity), "AntiKt_TruthFixed_r04");
  truthjetreco->set_algo_node("ANTIKTFIXED");
  truthjetreco->set_input_node("TRUTHFIXED");
  truthjetreco->Verbosity(verbosity);
  se->registerSubsystem(truthjetreco);

  // The writer's defaults match the inputs registered above (IHCal + OHCal + retowered EMCal
  // towers, truth embedding flag 0, AntiKt_Tower_r04 / AntiKt_TruthFixed_r04, jets >= 10 GeV).
  JewelTreeWriter *writer = new JewelTreeWriter("JewelTreeWriter", outputTrees);
  writer->Verbosity(verbosity);
  se->registerSubsystem(writer);

  Fun4AllInputManager *in = new Fun4AllDstInputManager("DSTin");
  in->AddFile(inputDST);
  se->registerInputManager(in);

  se->run();
  se->End();
  delete se;
  gSystem->Exit(0);
}
