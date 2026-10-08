#!/bin/bash
# One-time install of JEWEL 2.6.0 (stock hepforge release) + the PDF sets the
# configs need + the optional jet filter. Safe to re-run: each step is
# skipped if its output already exists.
#
#   ./install.sh
#
# Produces:
#   src/jewel-2.6.0/          unpacked + compiled source
#   bin/jewel-2.6.0-vac       vacuum (pp-like) binary
#   bin/jewel-2.6.0-simple    medium binary (simple Bjorken/Glauber medium)
#   bin/jetFilter             optional anti-kT event filter (filter/jetFilter.cc)
#   pdfsets/<set>/            LHAPDF6 sets not available in cvmfs

set -euo pipefail

this_dir=$(dirname "$(readlink -f "$BASH_SOURCE")")
set +u; source "$this_dir/env.sh"; set -u

TARBALL_URL="https://jewel.hepforge.org/downloads/?f=jewel-${JEWEL_VERSION}.tgz"
PDFSET_URL="https://lhapdfsets.web.cern.ch/current"
PDFSETS=(EPPS16nlo_CT14nlo_Au197)   # LHAPDF id 901200

SRC="$JEWEL_ROOT/src/jewel-${JEWEL_VERSION}"
mkdir -p "$JEWEL_ROOT/src" "$JEWEL_ROOT/bin" "$JEWEL_ROOT/pdfsets"

# --- 1. JEWEL source --------------------------------------------------------
if [[ ! -f "$SRC/jewel-${JEWEL_VERSION}.f" ]]; then
  echo ">> downloading JEWEL ${JEWEL_VERSION}"
  curl -fsSL -o "$JEWEL_ROOT/src/jewel-${JEWEL_VERSION}.tgz" "$TARBALL_URL"
  tar -xzf "$JEWEL_ROOT/src/jewel-${JEWEL_VERSION}.tgz" -C "$JEWEL_ROOT/src"
fi

# --- 2. build ---------------------------------------------------------------
# Point the stock Makefile at the cvmfs LHAPDF 6 and bake in an rpath for it
# and for the gcc-14 runtime (libgfortran/libstdc++), so the binaries run on
# a condor node without needing LD_LIBRARY_PATH.
if [[ ! -x "$SRC/jewel-${JEWEL_VERSION}-vac" || ! -x "$SRC/jewel-${JEWEL_VERSION}-simple" ]]; then
  echo ">> building JEWEL ${JEWEL_VERSION}"
  # (the stock link rule is "$(FC) -o $@ -L$(LHAPDF_PATH) ...", so the rpath
  # rides along in FC; -Wl options are ignored by the -c compile steps)
  GCC_LIB=$(dirname "$(readlink -f "$(gfortran -print-file-name=libgfortran.so)")")
  make -C "$SRC" \
    FC="gfortran -Wl,-rpath,${LHAPDF6_PREFIX}/lib -Wl,-rpath,${GCC_LIB}" \
    LHAPDF_PATH="${LHAPDF6_PREFIX}/lib" \
    2>&1 | tee "$SRC/build.log"
fi
for v in vac simple; do
  ln -sfn "$SRC/jewel-${JEWEL_VERSION}-$v" "$JEWEL_ROOT/bin/jewel-${JEWEL_VERSION}-$v"
done

# --- 3. PDF sets ------------------------------------------------------------
for set in "${PDFSETS[@]}"; do
  if [[ ! -f "$JEWEL_ROOT/pdfsets/$set/$set.info" ]]; then
    echo ">> downloading PDF set $set"
    curl -fsSL "$PDFSET_URL/$set.tar.gz" | tar -xz -C "$JEWEL_ROOT/pdfsets"
  fi
done

# --- 4. optional jet filter -------------------------------------------------
if [[ ! -x "$JEWEL_ROOT/bin/jetFilter" || "$JEWEL_ROOT/filter/jetFilter.cc" -nt "$JEWEL_ROOT/bin/jetFilter" ]]; then
  echo ">> building jetFilter"
  g++ -std=c++17 -O2 "$JEWEL_ROOT/filter/jetFilter.cc" \
    $(fastjet-config --cxxflags --libs) \
    -Wl,-rpath,"$(fastjet-config --prefix)/lib" \
    -o "$JEWEL_ROOT/bin/jetFilter"
fi

echo ">> checking binaries resolve all libraries"
for b in "$JEWEL_ROOT"/bin/*; do
  if ldd "$b" | grep -q "not found"; then
    echo "ERROR: $b has unresolved libraries:"; ldd "$b" | grep "not found"; exit 1
  fi
done
echo ">> install OK"
ls -l "$JEWEL_ROOT/bin"
