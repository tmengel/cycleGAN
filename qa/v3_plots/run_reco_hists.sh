#!/bin/bash
# one condor job of reco_hists.C: run_reco_hists.sh <sample> <run> <first> <last> <out>
source /sphenix/user/tmengel/cycleGAN/fun4all/env.sh
cd /sphenix/user/tmengel/cycleGAN/qa/v3_plots
root.exe -b -q -l "reco_hists.C(\"$1\",$2,$3,$4,\"$5\")"
