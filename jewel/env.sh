# Source this (do not execute) to get a JEWEL 2.6.0 runtime environment.
#   source /sphenix/user/tmengel/cycleGAN/jewel/env.sh
#
# The compiler/FastJet come from a frozen sPHENIX release (not the floating
# "new" nightly), and LHAPDF from its versioned cvmfs path, so a build made
# today still runs the same way next month.

export JEWEL_ROOT=/sphenix/user/tmengel/cycleGAN/jewel
export JEWEL_VERSION=2.6.0
export JEWEL_SPHENIX_RELEASE=ana.572

source /opt/sphenix/core/bin/sphenix_setup.sh -n ${JEWEL_SPHENIX_RELEASE} > /dev/null 2>&1

export LHAPDF6_PREFIX=/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/opt/sphenix/core/LHAPDF-6.5.4
export LD_LIBRARY_PATH=${LHAPDF6_PREFIX}/lib:${LD_LIBRARY_PATH}

# PDF sets: locally downloaded sets (e.g. the EPPS16 Au nPDF, not shipped in
# cvmfs) first, then the cvmfs sets (CT14nlo, cteq6l1, ...).
export LHAPDF_DATA_PATH=${JEWEL_ROOT}/pdfsets:${LHAPDF6_PREFIX}/share/LHAPDF
export LHAPATH=${LHAPDF_DATA_PATH}

export PATH=${JEWEL_ROOT}/bin:${PATH}
