#!/bin/bash
# Compares a new JEWEL production against Luis's genesis sample.
#
#   qa/jewel/run_compare.sh <new_vacuum.hepmc> <new_medium.hepmc> <outdir> [nevents=100000]
#
# Builds jetqa if needed, fills jet-QA histograms for the first <nevents>
# events of each of the four samples (selection: leading R=0.4 jet in
# |eta|<0.7 with pT>20 GeV = the genesis filter), and writes
# <outdir>/compare.pdf + page*.png, with a summary table on stdout.
#
# The new files are assumed UNFILTERED (N_gen = nevents). Luis's files are
# filtered, so his N_gen is reconstructed from the filter .stats files of
# the blocks that make up his first <nevents> events.

set -euo pipefail
here=$(dirname "$(readlink -f "$0")")
set +u; source /sphenix/user/tmengel/cycleGAN/jewel/env.sh; set -u

[[ $# -ge 3 ]] || { echo "usage: $0 <new_vacuum.hepmc> <new_medium.hepmc> <outdir> [nevents]" >&2; exit 1; }
NEWVAC=$1; NEWMED=$2; OUT=$3; N=${4:-100000}
GEN=/sphenix/tg/tg01/jets/luisvale/genesis
mkdir -p "$OUT"

if [[ ! -x "$here/jetqa" || "$here/jetqa.cc" -nt "$here/jetqa" ]]; then
  g++ -O2 -std=c++17 "$here/jetqa.cc" $(fastjet-config --cxxflags --libs) $(root-config --cflags --libs) \
    -Wl,-rpath,"$(fastjet-config --prefix)/lib" -Wl,-rpath,"$(root-config --libdir)" -o "$here/jetqa"
fi

# N_gen behind Luis's first N filtered events (mergeGenesis.sh takes blocks in
# version order; the last contributing block is used fractionally).
luis_ngen() {
  ls "$GEN/$1/blocks/"*.stats | sort -V | xargs cat | paste - - | \
    awk -v n="$N" '{tot += $2; pas += $4; if (pas >= n) { f = ($4 - (pas - n)) / $4; printf "%.0f\n", tot - $2 + $2 * f; exit }}'
}

"$here/jetqa" "$NEWVAC" "$OUT/new_vacuum.root" "$N" &
"$here/jetqa" "$NEWMED" "$OUT/new_medium.root" "$N" &
"$here/jetqa" "$GEN/vacuum/eventFile/vacuum.hepmc" "$OUT/luis_vacuum.root" "$N" &
"$here/jetqa" "$GEN/norecoil/eventFile/norecoil.hepmc" "$OUT/luis_norecoil.root" "$N" &
NGV=$(luis_ngen vacuum); NGM=$(luis_ngen norecoil)
wait
echo "Luis N_gen: vacuum $NGV, norecoil $NGM"

cd "$here"
root.exe -b -q "compare.C+(\"$OUT\", $N, $N, $NGV, $NGM)"
