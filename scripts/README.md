This directory contains a few basic scripts used to create single chainID PDBs 
with non-overlapping numbering in the matching of unbound and bound structure
and to match the bound and unbound structures

* `correct-chains.csh` : for entries with a single chain change the chainID to A and B, for receptor and ligand, respectively
* `make-docking-files.csh` : based on the matched entries, generate HADDOCK-ready directories with PDBs and input files for RMSD calculations
* `match-chains.csh` : for entries with multiple chain change the chainID to A and B, for receptor and ligand, respectively, and remove any overlap in numbering

Note that most scripts make use of the `pdb-tools` scripts available from [https://github.com/haddocking/pdb-tools](https://github.com/haddocking/pdb-tools)
