# Fun4All GEANT4 simulation of the JEWEL samples

Takes the HepMC files from `../jewel` (JEWEL vacuum/medium) through the sPHENIX detector
simulation in two condor passes:

| pass | input | output | cost (condor) |
|---|---|---|---|
| **1** GEANT4 | HepMC events | `G4Hits/G4Hits_<sample>-<run>-<seg>.root` (G4Hits + G4TruthInfo DST) | ~27 s/event, ~2.2 MB/event, 5.6 GB peak memory |
| **2** reco | one G4Hits DST | `trees/trees_*.root` (towertree, jettree, globaltree) and `images/images_*.root` (one TH2D per event) | ~5 min/job, mostly startup |

```
fun4all/
├── env.sh                      source for the runtime env (pinned: sPHENIX ana.572)
├── install.sh                  builds src/jeweltreewriter into install/ (done; rerun after editing it)
├── macros/
│   ├── Fun4All_G4_pass1.C      pass 1: HepMC -> GEANT4 -> G4Hits DST
│   ├── Fun4All_pass2.C         pass 2: towers, vertex, tower + truth jets -> trees
│   ├── make_images.C           pass 2: trees -> per-event 24x64 TH2D calorimeter images
│   └── calo/                   yeonjugo's no-noise calorimeter macros (snapshot, see below)
├── scripts/
│   ├── submit_pass1.sh         build the job list from HepMC files, submit pass 1
│   ├── submit_pass2.sh         one pass 2 job per pass 1 DST
│   └── run_pass1.sh / run_pass2.sh   what each condor job runs (fine interactively too)
└── src/jeweltreewriter/        the tree-writer Fun4All module
```

## Running a production

```bash
cd /sphenix/user/tmengel/cycleGAN/fun4all
H=/sphenix/tg/tg01/jets/tmengel/JEWEL_hepmc/pthat16

# pass 1: <sample> <tag> <runnumber> <events/job> <max events> <hepmc files...>
scripts/submit_pass1.sh vacuum pthat16 1000 10 100000 $H/vacuum/raw/*.hepmc
scripts/submit_pass1.sh medium pthat16 2000 10 100000 $H/medium/raw/*.hepmc

# resubmit only the pass 1 jobs whose DST is missing (same arguments + --only-missing)
scripts/submit_pass1.sh --only-missing vacuum pthat16 1000 10 100000 $H/vacuum/raw/*.hepmc

# pass 2, once pass 1 is done (or partly done: it takes whatever DSTs exist)
scripts/submit_pass2.sh vacuum pthat16
scripts/submit_pass2.sh medium pthat16
scripts/submit_pass2.sh --only-missing vacuum pthat16
```

Output: `/sphenix/tg/tg01/jets/tmengel/JEWEL_sim/<tag>/<sample>/{G4Hits,trees,images,logs,condor}`
plus `production.txt` (what was submitted, with which release and seeds). 100k events are
~220 GB of DSTs per sample, which is why output goes to tg01. Each job checks its own output
(DST event count = requested, every tree has one entry per input event, one image per event)
and exits non-zero otherwise, so failures show up as missing files plus a non-empty
`logs/*.err`.

Pass 1 needs ~5.4 GB no matter how many events a job has (measured 2026-10-02: flat
5.36-5.37 GB over 15 events). Request 6000 MB and use 10 events/job, as genesis did: a
4000 MB request gets every job killed by the 4096 MB cgroup limit, and condor counted 7.3 GB
for a 50-event job (cgroup accounting includes file cache). Only part of the farm offers
6 GB slots, so pass 1 jobs can sit idle before matching.

Run numbers: one per sample, 100-9999 (convention vacuum 1000, medium 2000). With the segment
they fix the random seeds.

## Reproducibility

* **Events:** `condor/pass1_jobs.list` maps segment k to a fixed block of events
  (`hepmc skip nevents`). It is written once per production and the submit script refuses to
  change it, so segment k always means the same events.
* **Seeds:** pass 1 = `run*100000 + segment`, pass 2 = `1000000000 + run*100000 + segment`, set
  through `RANDOMSEED`, so all of a job's random numbers (vertex smearing, GEANT4, tower
  digitization) are fixed. Any job can be rerun with `run_pass1.sh` / `run_pass2.sh` and gives
  the same result with the same release.
* **Release:** everything runs on the frozen `ana.572`, never the floating `new` nightly (a
  nightly once broke RetowerCEMC with no change of ours). CDB global tag `MDC2_ana.435` (as the
  genesis production), override with `F4A_CDBTAG`.

## What the passes do, and what changed from `../code/` (the genesis-era macros)

**Pass 1** (`Fun4All_G4_pass1.C`, from `code/pass1/Fun4All_G4_JEWEL_pass1.C`): HepMC input with
the sPHENIX beam parameters, then GEANT4 with beam pipe, MBD, EPD, EMCal, iHCal, magnet,
oHCal, plug door and black hole; no tracking detectors (unused by the images, most of the CPU).
Fixed relative to the old macro:

* **Vertex distribution.** The old macro took the beam settings from `RunSettings(1000)`,
  which maps every run >= 100 to `mRad_00`: 2024 pp zero-crossing-angle beam, vertex
  sigma_z = **65 cm** (the genesis trees have |vz| up to 125 cm). Now an explicit choice,
  default **`fixed`: every event at vertex (0, 0, 0)**, head-on beams, no smearing (JEWEL
  writes all vertices at the origin, and the sPHENIX beam smearing is simply not applied).
  Alternatives via `--beam`: `AuAu` (sigma_z 13.5 cm, 1 mrad crossing, 2024 beam spot), `pp`
  (16 cm, 1.5 mrad), `pp_zeroangle` (the old 65 cm behaviour).
* **Truth jets removed from pass 1.** The old macro clustered them with embedding flag 1, but
  HepMC particles carry embedding id 0, so every pass 1 truth-jet container was empty (the
  trees stage worked around it). Pass 2 now makes all jets.
* GlobalVertexReco moved to pass 2, the `Production_*` output shuffling removed (the DST is
  written locally and moved by the run script), fixed seeds added.

**Pass 2** (`Fun4All_pass2.C` + `make_images.C`, from `code/trees/Fun4All_JEWEL_MakeTrees.C`
and yeonjugo's `draw_event_v2.C`):

* Vertex: `GlobalVertexReco` finds no MBD/tracking vertex in a simulation-only DST and falls
  back to the **true** primary vertex. Tower jets are corrected for it, and it is stored in
  `globaltree` (vx, vy, vz).
* Towers: yeonjugo's `G4_CEmc_Spacal.C` / `G4_HcalIn_ref.C` / `G4_HcalOut_ref.C` (snapshot in
  `macros/calo/`, md5 identical to the versions chosen on 2026-09-15): **no noise** (pedestal
  width 0) and a **fixed EMCal calibration** from the TowerCalibCombinedParams_2020 XML.
  The build's central `<G4_*.C>` macros give +43% tower energy plus noise; don't swap them in
  without deciding to change the images.
* EMCal retowered onto the 24 x 64 HCal grid; anti-kT R=0.4 jets from towers
  (`AntiKt_Tower_r04`) and from truth final-state particles (`AntiKt_TruthFixed_r04`,
  embedding flag 0).
* `JewelTreeWriter` writes `towertree` (EMCal/iHCal/oHCal energy per tower), `jettree` (reco
  and truth jets >= 10 GeV with constituents; reco constituents carry their calorimeter,
  truth ones energy and PDG id) and `globaltree` (vertex). The branch layout is documented in
  `src/jeweltreewriter/JewelTreeWriter.h` and is unchanged, so the existing QA macros in
  `../code/qa` read these trees.
* `make_images.C`: identical images to `draw_event_v2.C` (24 x 64, eta [-1.1, 1.1], phi
  [0, 2 pi], EMCal + iHCal + oHCal energy, negatives clipped, histograms named
  `h_eta_phi_cent0_file<seg>_evt<i>`), with explicit paths. Image i = entry i of the trees.
  (The phi axis now uses the exact pi; draw_event_v2.C's `3.14156` only shifted axis labels.)

## Checked (2026-10-01)

Interactive test, 3 vacuum events (with the then-default AuAu beam): pass 1 -> pass 2 -> images
all ran; the vertex was at the AuAu beam spot (x, y = -0.06, 0.13 cm; z within +-16 cm); truth
jets present and matching the reco jets in eta; image sums equal the tower sums; 100% of jet
constituents resolved.

## Before a big submission

Run a few jobs first (e.g. `submit_pass1.sh vacuum test 1000 3 6 <one hepmc>`, then
`submit_pass2.sh vacuum test`) and check `logs/` and the outputs, especially after any change
of release.
