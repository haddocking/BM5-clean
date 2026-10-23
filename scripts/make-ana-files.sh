#!/bin/bash
# Interface files for the analysis, from a reference complex (receptor chain A first, ligand chain B):
# target.contacts3.9, target.contacts5 (intermolecular contacts), target.izone, target.izoneA,
# target.izoneB (ProFit interface zones, from the contacts at 10 A) and target.lzone (ProFit
# zones for the ligand RMSD, from the contacts at 200 A).
#
# Usage: make-ana-files.sh TARGET.pdb OUTDIR
#
# Needs contact-chainID (HADDOCK 2.4 tools; cutoffs up to 200 A) and make_izone.csh
# (haddock-CSB-tools/CAPRI_analysis), from the PATH or the CONTACT_CHAINID and MAKE_IZONE
# variables, and join_zones.py (in this directory).
set -e
CONTACT=${CONTACT_CHAINID:-contact-chainID}
IZONE=${MAKE_IZONE:-make_izone.csh}
JOIN="python3 $(cd "$(dirname "$0")" && pwd)/join_zones.py"
target=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
mkdir -p "$2"; cd "$2"
cp "$target" target.tmp.pdb
$CONTACT target.tmp.pdb 3.9 >target.contacts3.9
$CONTACT target.tmp.pdb 5 >target.contacts5
$CONTACT target.tmp.pdb 10 >target.contacts10
$CONTACT target.tmp.pdb 200 >target.contacts200
csh $IZONE target.contacts10 >tmp
$JOIN tmp >target.izone 2>/dev/null
grep "zone A" target.izone >target.izoneA || true
grep "zone B" target.izone >target.izoneB || true
csh $IZONE target.contacts200 >tmp
$JOIN tmp >tt 2>/dev/null
grep "zone A" tt >target.lzone; echo fit >>target.lzone
grep "zone B" tt | sed -e 's/zone/rzone/g' >>target.lzone
rm -f tmp tt target.contacts10 target.contacts200 target.tmp.pdb
