# Source this (do not execute) to get the Fun4All runtime environment for both passes.
#   source /sphenix/user/tmengel/cycleGAN/fun4all/env.sh
#
# Pinned to a frozen sPHENIX release, never the floating "new" nightly: a nightly once broke
# RetowerCEMC with no change on our side.

export F4A_ROOT=/sphenix/user/tmengel/cycleGAN/fun4all
export F4A_RELEASE=ana.572

source /opt/sphenix/core/bin/sphenix_setup.sh -n ${F4A_RELEASE} > /dev/null 2>&1

# The JewelTreeWriter module (src/jeweltreewriter), built against $F4A_RELEASE by install.sh into
# a prefix that holds nothing else. APPENDED to the search paths, never prepended (as
# setup_local.sh would): a prefix that also held e.g. an old libjetbackground.so would otherwise
# shadow the release's copy.
export F4A_INSTALL=${F4A_ROOT}/install
export ROOT_INCLUDE_PATH=${ROOT_INCLUDE_PATH}:${F4A_INSTALL}/include:${F4A_ROOT}/macros
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${F4A_INSTALL}/lib

# Default location for simulation output (big files: keep them off /sphenix/user).
export F4A_OUTPUT_ROOT=${F4A_OUTPUT_ROOT:-/sphenix/tg/tg01/jets/tmengel/JEWEL_sim}
