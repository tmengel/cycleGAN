#!/bin/bash
# Builds the JewelTreeWriter Fun4All module (src/jeweltreewriter) against the pinned release
# into install/. Re-run after changing the module's source; the build directory is recreated.
#
#   ./install.sh

set -euo pipefail
here=$(dirname "$(readlink -f "$0")")
set +u; source "$here/env.sh"; set -u

build="$here/build/jeweltreewriter"
rm -rf "$build"
mkdir -p "$build"
cd "$build"
"$here/src/jeweltreewriter/autogen.sh" --prefix="$F4A_INSTALL" > autogen.log 2>&1 || { tail -20 autogen.log; exit 1; }
make -j4 install > make.log 2>&1 || { tail -30 make.log; exit 1; }

echo ">> installed into $F4A_INSTALL (built against $F4A_RELEASE: $OFFLINE_MAIN)"
ls -l "$F4A_INSTALL/lib"
# the prefix must hold only this module (see env.sh)
extra=$(ls "$F4A_INSTALL/lib" | grep -v '^libjeweltreewriter' || true)
[[ -z "$extra" ]] || { echo "ERROR: unexpected libraries in $F4A_INSTALL/lib: $extra" >&2; exit 1; }
ldd "$F4A_INSTALL/lib/libjeweltreewriter.so" | grep "not found" && { echo "ERROR: unresolved libraries" >&2; exit 1; }
echo ">> install OK"
