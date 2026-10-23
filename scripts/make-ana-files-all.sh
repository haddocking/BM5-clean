#!/bin/bash
# Regenerate the analysis files of every complex: ana_scripts/ from XXXX_target.pdb, and the
# symmetry-related references ana_scripts_K/ from their target.pdb (contacts5, zones).
# Usage: make-ana-files-all.sh [HADDOCK-ready directory]
H=$(cd "$(dirname "$0")" && pwd)
D=${1:-${BM5_HADDOCK_READY:-$H/../HADDOCK-ready}}
TMP=$(mktemp -d)
one() {
  c=$1
  bash $H/make-ana-files.sh $D/$c/${c}_target.pdb $D/$c/ana_scripts || echo "FAILED $c"
  for k in $D/$c/ana_scripts_*; do
    [ -d $k ] || continue
    bash $H/make-ana-files.sh $k/target.pdb $TMP/$c-$(basename $k) || { echo "FAILED $k"; continue; }
    cp $TMP/$c-$(basename $k)/{target.contacts5,target.izone,target.izoneA,target.izoneB,target.lzone} $k/
  done
}
export -f one; export H D TMP
ls -d $D/*/ana_scripts | xargs -n1 dirname | xargs -n1 basename | xargs -P ${NPROC:-8} -I{} bash -c 'one {}'
rm -rf $TMP
