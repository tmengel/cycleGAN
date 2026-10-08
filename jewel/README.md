# JEWEL 2.6.0 HepMC generation

Generates JEWEL vacuum and medium-modified (AuAu 200 GeV, 0-10%) HepMC2 files.
The defaults reproduce the setup of Luis's "genesis" production (with one bug
fixed, see below).

```
jewel/
├── install.sh          one-time: download + build JEWEL 2.6.0, fetch PDF set, build jetFilter
├── env.sh              source for a runtime env (pinned: sPHENIX ana.572 + cvmfs LHAPDF 6.5.4)
├── configs/
│   ├── vacuum/         jewel-2.6.0-vac
│   └── medium/         jewel-2.6.0-simple, no recoils (+ medium.params.dat)
│       each: config.sh, params.tmpl.dat, tables/ (built by make_tables.sh)
├── scripts/
│   ├── make_tables.sh  build a config's integration tables (once per config, ~30 s)
│   ├── run_jewel.sh    run one job (what condor runs; also fine interactively)
│   ├── submit.sh       submit a production to condor
│   └── merge_hepmc.sh  concatenate per-job files into one HepMC file
├── filter/jetFilter.cc optional anti-kT event filter
├── seeds.ledger        every NJOB (seed) range ever submitted
└── src/ bin/ pdfsets/  build products
```

## Quick start

```bash
cd /sphenix/user/tmengel/cycleGAN/jewel
./install.sh                         # already done; safe to re-run
scripts/make_tables.sh vacuum        # already done for vacuum and medium
scripts/make_tables.sh medium

# interactive test: 200 events
scripts/run_jewel.sh medium 899999 200 /tmp/$USER/jewel_test

# production: <config> <tag> <njobs> <nevent/job>; seed = job number 1..njobs
# (prepared for pthat16 with --dry-run; run without it to submit)
scripts/submit.sh vacuum pthat16 100 10000
scripts/submit.sh medium pthat16 100 10000
# more statistics later: jobs 101-200
scripts/submit.sh medium pthat16 100 10000 --first-job 101

# afterwards: one file with exactly 100k events
scripts/merge_hepmc.sh vacuum_100k.hepmc --max 100000 \
    /sphenix/tg/tg01/jets/tmengel/JEWEL_hepmc/pthat16/vacuum/raw/*.hepmc

# Optional: add "--filter R ptmin etamax" to submit.sh (or the same three
# args to run_jewel.sh) to also write filtered/ files keeping only events
# with an anti-kT jet above ptmin inside |eta| < etamax, e.g. the genesis
# selection "--filter 0.4 20 0.7".
```

Output goes to `$JEWEL_OUTPUT_ROOT/<tag>/<config>/` (default root
`/sphenix/tg/tg01/jets/tmengel/JEWEL_hepmc`), with `raw/` (every generated
event), `filtered/` (with `--filter`, plus `.stats` files of event counts and
summed weights), `logs/`, `config/` (snapshot of the exact params used), and
`condor/`. Every per-job file is a complete HepMC2 file.

Measured cost (one core, 2026-10-01): with PTMIN 5-16 and ETAMAX 2.0 the
medium takes ~0.1-0.2 s/event (pT-hat 5-10: 0.12; PTMIN 5 open-ended: 0.21),
vacuum ~0.01 s/event, so a 10k-event medium job takes ~20-35 min. Raw output
is ~6.5 kB/event.

## Seeds and reproducibility

JEWEL seeds PYTHIA with `NJOB*1000` (NJOB < 900000). `submit.sh` sets
**NJOB = job number** (1, 2, ...; `--first-job` to extend), so job k of a
production can always be regenerated:

    scripts/run_jewel.sh <config> <k> <nevent> <outdir>

This reproduces `raw/<config>_<k>.hepmc` byte for byte, provided the config
(fingerprint), the tables, `<nevent>` and the binary are unchanged. Checked
2026-10-01: same seed run twice locally and on a condor node gives identical
files. A different `<nevent>` gives different events, because the pp/pn/np/nn
event counts are drawn first. Each production directory has a
`production.txt` recording all of this (fingerprint, tables md5, binary md5,
environment), plus `config/` with the exact params.

`seeds.ledger` records every submitted range with its fingerprint. Reusing
job numbers is refused only for the same config with identical physics, since
only that would duplicate events (`--reproduce` overrides it for a
deliberate rerun). Vacuum and medium may use the same job numbers: with the
same NJOB they produce different events, because the medium setup consumes
random numbers. Smoke tests used 899100+.

## Physics settings (and what changed relative to genesis)

Both configs: `SQRTS 200`, `PDFSET 901200` (EPPS16nlo_CT14nlo_Au197, downloaded
to `pdfsets/`), `MASS 197 NPROTON 79`, `PTMIN 16 PTMAX -1` (parton pT-hat, no upper limit; as genesis), `ETAMAX 2.0` (medium extent; genesis used 1.1),
`WEIGHTED T WEXPO 8`. Medium: `TI 0.375 TAUI 0.4 TC 0.17`, 0-10%,
`SIGMANN 4.2` (42 mb), `MDSCALEFAC 1 MDFACTOR 0.45`, `KEEPRECOILS F`.

* **Fixed: integration tables.** JEWEL caches `xsecs.dat`/`pdfs.dat`/`splitint.dat`.
  `xsecs.dat` is tabulated over the Debye-mass range `[3*T_c, 3*T_max]`, which
  depends on the medium. Genesis reused one table set, built for a TI = 0.26
  medium (m_D grid 0.51-0.78 GeV), for both its vacuum and its TI = 0.375
  run, whose m_D reaches ~1.12 GeV. Here each config builds its own
  tables. In practice the impact was negligible: Luis's medium runs logged
  9418 out-of-range lookups over ~270k events, and a rerun with correct tables
  but genesis ETAMAX reproduces his jet observables (see Validation). The tables are fingerprinted against the params, so editing a config
  without rerunning `make_tables.sh` is an error, not a silent mismatch.
* **Explicit medium parameters.** 2.6.0 changed the built-in medium defaults to
  5 TeV values (TI 0.59, MDSCALEFAC 1). All keys are now set explicitly.
  Luis's `A 197` line was dropped because medium-simple.f has no such key (A
  comes from `MASS`).
* **Events are weighted.** With `WEXPO 8` the pT-hat spectrum is flattened and
  each event carries a weight (first weight on the HepMC `E` line). Any
  physical distribution (jet pT spectrum, R_AA, ...) needs these weights. If
  you'd rather train on unweighted samples, set `WEIGHTED F` (much slower to
  populate high pT) or resample by weight.
* **TI is the main knob.** 0.375 GeV at TAUI 0.4 fm is kept for continuity with
  genesis. For reference, the pre-2.6 JEWEL RHIC setting (TAUI 0.6, TI 0.36)
  Bjorken-scales to TI ≈ 0.41 at TAUI 0.4. The stock sPHENIX JEWEL 2.2.0 file
  uses 0.26, which is barely above T_c.
* NEVENT is approximate: with a nuclear PDF (ISOCHANNEL XX) JEWEL draws the
  pp/pn/np/nn event counts from Poisson distributions around NEVENT·Z²/A² etc.,
  so a job gives NEVENT ± ~√NEVENT events.
* ETAMAX does **not** restrict the hard scattering, which is generated at all
  rapidities. It sets the medium's rapidity
  extent, and partons beyond |eta| > ETAMAX+1 are not showered by JEWEL.

## Validation against Luis's genesis sample (2026-10-01)

`../qa/jewel/run_compare.sh <new_vacuum.hepmc> <new_medium.hepmc> <outdir>`
compares the first 100k events of a production with Luis's genesis files. It
applies the genesis selection (leading R=0.4 jet, |eta|<0.7, pT>20) to both and
writes `compare.pdf` plus a summary table. Results for `pthat16` (plots in
`/sphenix/tg/tg01/jets/tmengel/JEWEL_pp200_signal/v3_2026-10-03/qa/generator_vs_genesis/`):

* vacuum: all observables agree (pass rate 51.5% vs 51.5%; leading-jet pT,
  eta, mass, girth, z, radial profile, multiplicities, weights: chi2/ndf
  0.5-1.9, identical means).
* medium: leading-jet observables and the medium/vacuum yield ratio agree (0.425
  vs 0.436). Event-wide multiplicity is lower (57.1 vs 59.1) and additional
  jets slightly softer. That comes from ETAMAX 2.0, which quenches partons at
  1.1<|y|<2 whose lost energy (no recoils) leaves the event. A 20k-event
  check with ETAMAX 1.1 and new tables gives 59.9 +- 0.5, matching Luis.

## Adding a config

Copy a config directory, edit `params.tmpl.dat` / `medium.params.dat`, and run
`scripts/make_tables.sh <new>`. It prints the parameter block JEWEL actually
read. Check it: JEWEL silently keeps the default for any value it fails to
parse. For example, a medium with recoils needs `KEEPRECOILS T` (and usually
`WRITESCATCEN T`) plus a background subtraction downstream.

## Notes

* Build: stock hepforge tarball (byte-identical `jewel-2.6.0.f` to Luis's
  build), gfortran 14.2 from `ana.572`, linked against
  `/cvmfs/.../opt/sphenix/core/LHAPDF-6.5.4` with an rpath, so the binaries
  don't depend on the floating `new` nightly.
* JEWEL 2.6.0 uses LHAPDF **6** (`pythia6425mod-lhapdf6.f`), not LHAPDF 5.
