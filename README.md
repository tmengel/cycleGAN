# cycleGAN JEWEL simulation

Code used to produce the v3 JEWEL pp200 signal dataset
(`/sphenix/tg/tg01/jets/tmengel/JEWEL_pp200_signal/v3_2026-10-03/`: vacuum, medium, vacuum_hijing, medium_hijing).

- [`jewel/`](jewel/README.md) — JEWEL 2.6.0 HepMC generation (vacuum and medium/norecoil configs, pthat16).
- [`fun4all/`](fun4all/README.md) — Fun4All two-pass sPHENIX simulation: pass 1 HepMC → GEANT4 hits DST,
  pass 2 towers + jets → JewelTreeWriter trees + calorimeter images, plus the HIJING image overlay.

Both pin sPHENIX release `ana.572` and expect to live at `/sphenix/user/tmengel/cycleGAN` (see each `env.sh`).
Run each `install.sh` once after cloning; JEWEL source, PDF sets and JEWEL tables are downloaded/generated, not stored here.
