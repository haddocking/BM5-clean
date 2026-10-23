#!/bin/bash
# Recompute the i-RMSD tables of HADDOCK-ready with ProFit (irmsd.sh):
#   i-RMSD.dat   : XXXX_target.pdb vs ana_scripts/target-unbound.pdb over target.izone
#   i-RMSD_r.dat : XXXX_r_b-matched.pdb vs XXXX_r_u.pdb over target.izoneA
#   i-RMSD_l.dat : XXXX_l_b-matched.pdb vs XXXX_l_u.pdb over target.izoneB
# and the -sorted versions (smallest first). Usage: make-irmsd-tables.sh [HADDOCK-ready directory]
H=$(cd "$(dirname "$0")" && pwd)
D=${1:-${BM5_HADDOCK_READY:-$H/../HADDOCK-ready}}
TMP=$(mktemp -d)
one() {
  c=$1; cd $D
  a=$(bash $H/irmsd.sh $c/${c}_target.pdb $c/ana_scripts/target-unbound.pdb $c/ana_scripts/target.izone </dev/null)
  r=$(bash $H/irmsd.sh $c/${c}_r_b-matched.pdb $c/${c}_r_u.pdb $c/ana_scripts/target.izoneA </dev/null)
  l=$(bash $H/irmsd.sh $c/${c}_l_b-matched.pdb $c/${c}_l_u.pdb $c/ana_scripts/target.izoneB </dev/null)
  echo "$c $a $r $l" > $TMP/$c
}
export -f one; export H D TMP
ls -d $D/*/ana_scripts | xargs -n1 dirname | xargs -n1 basename | xargs -P ${NPROC:-8} -I{} bash -c 'one {}'
cd $D
{ echo "#struc i-RMSD"; for f in $(ls $TMP); do read c a r l < $TMP/$f; echo "$c/ana_script/target-unbound.pdb $a"; done; } > i-RMSD.dat
{ echo "#struc i-RMSD"; for f in $(ls $TMP); do read c a r l < $TMP/$f; echo "$c/${c}_r_u.pdb $r"; done; } > i-RMSD_r.dat
{ echo "#struc i-RMSD"; for f in $(ls $TMP); do read c a r l < $TMP/$f; echo "$c/${c}_l_u.pdb $l"; done; } > i-RMSD_l.dat
for f in i-RMSD i-RMSD_r i-RMSD_l; do sort -s -k2,2g $f.dat > $f-sorted.dat; done
rm -rf $TMP
