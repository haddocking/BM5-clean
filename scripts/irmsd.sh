#!/bin/bash
# Interface RMSD with ProFit (backbone N, CA, C, O over a zone file), as in the i-RMSD*.dat files.
# Usage: irmsd.sh REFERENCE.pdb MOBILE.pdb ZONEFILE
{ echo "refe $1"; echo "mobi $2"; echo "atoms CA,C,N,O"; cat $3; echo "fit"; echo "quit"; } |
  profit 2>/dev/null | awk '/RMS:/{r=$2} END{printf "%.3f\n", r}'
