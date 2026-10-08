# Shared helpers for the JEWEL scripts (sourced, not executed).

JEWEL_SCRIPTS=$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")
set +u; source "$JEWEL_SCRIPTS/../env.sh"; set -u

# Default location for produced datasets (big files: keep them off /sphenix/user).
: "${JEWEL_OUTPUT_ROOT:=/sphenix/tg/tg01/jets/tmengel/JEWEL_hepmc}"

die() { echo "ERROR: $*" >&2; exit 1; }

# load_config <name>: sets CONFIG_DIR, JEWEL_BINARY, JEWEL_EXE
load_config() {
  CONFIG_NAME=$1
  CONFIG_DIR="$JEWEL_ROOT/configs/$CONFIG_NAME"
  [[ -f "$CONFIG_DIR/config.sh" ]] || die "no such config '$CONFIG_NAME' (looked in $CONFIG_DIR)"
  [[ -f "$CONFIG_DIR/params.tmpl.dat" ]] || die "$CONFIG_DIR has no params.tmpl.dat"
  source "$CONFIG_DIR/config.sh"
  JEWEL_EXE="$JEWEL_ROOT/bin/jewel-${JEWEL_VERSION}-${JEWEL_BINARY}"
  [[ -x "$JEWEL_EXE" ]] || die "$JEWEL_EXE not found -- run install.sh first"
}

# Fingerprint of everything the cached tables depend on: the physics lines of
# the params template (per-job @...@ lines and comments stripped) plus the
# medium params file, if any.
config_fingerprint() {
  {
    grep -v -e '^#' -e '@' -e '^[[:space:]]*$' "$CONFIG_DIR/params.tmpl.dat"
    if [[ -f "$CONFIG_DIR/medium.params.dat" ]]; then
      grep -v -e '^#' -e '^[[:space:]]*$' "$CONFIG_DIR/medium.params.dat"
    fi
    echo "binary=$JEWEL_BINARY version=$JEWEL_VERSION"
  } | md5sum | cut -d' ' -f1
}

# check_tables: fail unless the config's tables exist and match the params.
check_tables() {
  local t="$CONFIG_DIR/tables"
  for f in xsecs.dat pdfs.dat splitint.dat fingerprint; do
    [[ -s "$t/$f" ]] || die "config '$CONFIG_NAME' has no tables ($t/$f missing) -- run scripts/make_tables.sh $CONFIG_NAME"
  done
  [[ "$(cat "$t/fingerprint")" == "$(config_fingerprint)" ]] || \
    die "tables in $t were built for different parameters -- rerun scripts/make_tables.sh $CONFIG_NAME"
}

# write_params <out> <njob> <nevent> <hepmcfile> <logfile>
write_params() {
  sed -e "s|@NJOB@|$2|" -e "s|@NEVENT@|$3|" \
      -e "s|@HEPMCFILE@|$4|" -e "s|@LOGFILE@|$5|" \
      "$CONFIG_DIR/params.tmpl.dat" > "$1"
  ! grep -v '^#' "$1" | grep -q '@' || die "unfilled placeholder in $1"
}

# JEWEL seeds PYTHIA's RANMAR with MRPY(1) = NJOB*1000, which must stay below
# 9e8, and two jobs with the same NJOB produce identical event streams.
check_njob() {
  [[ "$1" =~ ^[0-9]+$ ]] && (( $1 >= 1 && $1 < 900000 )) || die "NJOB=$1 must be an integer in [1, 899999]"
}
