#!/bin/bash
source /sphenix/user/tmengel/cycleGAN/fun4all/env.sh
cd /sphenix/user/tmengel/cycleGAN/qa/fun4all
root.exe -b -q 'compare_sim.C+("/sphenix/tg/tg01/jets/tmengel/JEWEL_pp200_signal/v3_2026-10-03/qa/compare_v2", "/sphenix/tg/tg01/jets/tmengel/JEWEL_pp200_signal/v3_2026-10-03/trees/vacuum/trees_*.root", "/sphenix/tg/tg01/jets/tmengel/JEWEL_pp200_signal/v3_2026-10-03/trees/medium/trees_*.root", "/sphenix/tg/tg01/jets/tmengel/JEWEL_pp200_signal/v2_2026-09-14_genesis100k/trees/merged/vacuum_trees_merged.root", "/sphenix/tg/tg01/jets/tmengel/JEWEL_pp200_signal/v2_2026-09-14_genesis100k/trees/merged/norecoil_trees_merged.root")'
echo "EXIT $?"
