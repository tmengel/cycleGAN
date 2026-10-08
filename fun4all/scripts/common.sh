# Shared helpers for the Fun4All scripts (sourced, not executed).

F4A_SCRIPTS=$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")
set +u; source "$F4A_SCRIPTS/../env.sh"; set -u

die() { echo "ERROR: $*" >&2; exit 1; }

# Fixed random seeds, so every job can be rerun exactly: pass1 run*100000 + segment, pass2
# 1000000000 + run*100000 + segment (both fit in an int for run < 10000, segment < 1000000).
# Seeds of two samples cannot collide as long as their run numbers differ by at least 10
# (convention: vacuum 1000, medium 2000).
check_run_segment() {
  [[ "$1" =~ ^[0-9]+$ ]] && (( $1 >= 100 && $1 < 10000 )) || die "run number $1 must be in [100, 9999]"
  [[ "$2" =~ ^[0-9]+$ ]] && (( $2 < 1000000 )) || die "segment $2 must be in [0, 999999]"
}
pass1_seed() { echo $(( $1 * 100000 + $2 )); }
pass2_seed() { echo $(( 1000000000 + $1 * 100000 + $2 )); }

# file names: <kind>_<sample>-<run, 10 digits>-<segment, 6 digits>.root
f4a_name() { printf '%s_%s-%010d-%06d.root' "$1" "$2" "$3" "$4"; }

# work in condor's scratch dir (local disk), or a temp dir interactively
enter_workdir() {
  if [[ -n "${_CONDOR_SCRATCH_DIR:-}" && -d "${_CONDOR_SCRATCH_DIR}" ]]; then
    cd "$_CONDOR_SCRATCH_DIR"
  else
    F4A_TMP=$(mktemp -d)
    trap 'rm -rf "$F4A_TMP"' EXIT
    cd "$F4A_TMP"
  fi
}

# copy to a temp name, then rename, so a half-copied file never looks finished
ship() {
  mkdir -p "$2"
  cp "$1" "$2/.$(basename "$1").part"
  mv "$2/.$(basename "$1").part" "$2/$(basename "$1")"
}

# number of entries of a tree in a ROOT file (-1 if missing)
tree_entries() {
  root.exe -q -b -l -e 'TFile f("'"$1"'"); TTree *t = (TTree*) f.Get("'"$2"'"); printf("%lld\n", t ? t->GetEntries() : -1LL); exit(0);' 2>/dev/null | tail -1
}

# names of the files in a directory, one per line, from a single unsorted readdir (no per-file
# stat: much faster than [[ -f ]] per file on a busy tg01)
dir_names() { ls -f "$1" 2>/dev/null | grep -v '^\.\.\?$' || true; }
