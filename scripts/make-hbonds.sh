#!/bin/bash
# XXXX_hbonds.tbl of every complex: haddock-restraints restraint on the unbound receptor and ligand
# (XXXX_r_u.pdb, XXXX_l_u.pdb), concatenated. The restraints keep the disconnected parts of each
# partner together (chain breaks, merged chains and co-factors). No file is written (an existing
# one is removed) when there is nothing to restrain. Columns 77-80 are rebuilt in the copy read
# by the tool (some files hold other data there).
# Usage: make-hbonds.sh [HADDOCK-ready directory]
D=${1:-${BM5_HADDOCK_READY:-$(cd "$(dirname "$0")/../HADDOCK-ready" && pwd)}}
TMP=$(mktemp -d)
one() {
  c=$1; out=$TMP/${c}_hbonds.tbl
  : > $out
  for s in r l; do
    cut -c1-76 $D/$c/${c}_${s}_u.pdb | pdb_element > $TMP/${c}_${s}.pdb
    haddock-restraints restraint $TMP/${c}_${s}.pdb >> $out || echo "FAILED $c $s"
  done
  if grep -q assign $out; then cp $out "$D/$c/${c}_hbonds.tbl"; else rm -f "$D/$c/${c}_hbonds.tbl"; fi
}
export -f one; export D TMP
ls -d $D/*/ana_scripts | xargs -n1 dirname | xargs -n1 basename | xargs -P ${NPROC:-8} -I{} bash -c 'one {}'
rm -rf $TMP
