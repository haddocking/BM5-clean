#!/bin/bash
# True-interface ambiguous restraints of every complex, from the bound reference XXXX_target.pdb:
# XXXX_ambig.tbl (3.9 A cutoff) and XXXX_ambig5.tbl (5 A), with haddock-restraints ti.
# HETATM records are removed first (no accessibility for co-factors), and columns 77-80 rebuilt
# (some files hold other data there). Usage: make-ambig.sh [HADDOCK-ready directory]
D=${1:-${BM5_HADDOCK_READY:-$(cd "$(dirname "$0")/../HADDOCK-ready" && pwd)}}
TMP=$(mktemp -d)
one() {
  c=$1
  pdb_delhetatm $D/$c/${c}_target.pdb | cut -c1-76 | pdb_element > $TMP/$c.pdb
  haddock-restraints ti $TMP/$c.pdb 3.9 > $D/$c/${c}_ambig.tbl || echo "FAILED $c 3.9"
  haddock-restraints ti $TMP/$c.pdb 5 > $D/$c/${c}_ambig5.tbl || echo "FAILED $c 5"
}
export -f one; export D TMP
ls -d $D/*/ana_scripts | xargs -n1 dirname | xargs -n1 basename | xargs -P ${NPROC:-8} -I{} bash -c 'one {}'
rm -rf $TMP
