#!/bin/bash
# .info statistics files (pdb_wc, plus whether the residue numbering has gaps) for PDB files.
# Usage: make-info.sh FILE.pdb...   (writes FILE.info next to each file)
for f in "$@"; do
  gaps=$(awk '/^(ATOM|HETATM)/{k=substr($0,22,5); if(k!=p){n=substr($0,23,4)+0; c=substr($0,22,1); if(p!="" && c==pc && n>pn+1) g=1; pn=n; pc=c; p=k}} END{print g?"True":"False"}' "$f")
  { pdb_wc "$f"; printf "Has seq. gaps:\t%s\n" $gaps; } > "${f%.pdb}.info"
done
