#!/bin/bash
# Runs one QA/plotting macro in its own condor scratch dir (ACLiC writes its build products next
# to the macro, so each job compiles a private copy).
# $1: macro file name, found in qa/genesis
# $2..: macro arguments in order. Each is quoted as a string unless it is an integer.
source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.572
hostname
macro=$1; shift
work=${_CONDOR_SCRATCH_DIR:-$(mktemp -d)}
for d in genesis; do
  [ -f /sphenix/user/tmengel/cycleGAN/qa/$d/$macro ] && cp /sphenix/user/tmengel/cycleGAN/qa/$d/$macro $work/
done
cd $work
args=""
for a in "$@"; do
  if [[ "$a" =~ ^-?[0-9]+$ ]]; then args="${args:+$args,}$a"; else args="${args:+$args,}\"$a\""; fi
done
echo "running $macro($args)"
root.exe -l -b -q "$macro+($args)"
echo "script done"
