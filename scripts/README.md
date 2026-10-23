# scripts

## Regenerating the HADDOCK-ready restraint and analysis files

These scripts regenerate, for all complexes, the files derived from the PDB files of
`HADDOCK-ready`. Run them after any change to the structures. They work on `../HADDOCK-ready`
by default (or the directory given as argument, or in `BM5_HADDOCK_READY`) and run 8 complexes in
parallel (`NPROC` to change).

Requirements:

* [pdb-tools](https://github.com/haddocking/pdb-tools) (`pdb_delhetatm`, `pdb_element`)
* `haddock-restraints` ([haddocking/haddock-restraints](https://github.com/haddocking/haddock-restraints))
* `contact-chainID` from the HADDOCK 2.4 tools (cutoffs up to 200 A are needed) and `make_izone.csh`
  from haddock-CSB-tools (`CAPRI_analysis`), in the `PATH` or given by the `CONTACT_CHAINID` and
  `MAKE_IZONE` variables
* [ProFit](http://www.bioinf.org.uk/software/profit/) (`profit`) for the i-RMSD values
* Python 3 with numpy (`rebuild-symmetry-refs.py`)

In order:

0. `superimpose-unbound.py [COMPLEX...]` : superimposes each unbound partner onto its matched bound partner
   over its interface (backbone of `target.izoneA` or `target.izoneB`)
1. `make-target-unbound.py` : rebuild `ana_scripts/target-unbound.pdb` from `XXX_r_u.pdb` and `XXX_l_u.pdb`
   (`--write`; without it, only lists the complexes whose file differs)
2. `make-ana-files-all.sh` : intermolecular contacts at 3.9 and 5A (`target.contacts3.9`, `target.contacts5`),
   ProFit interface zones (`target.izone`, `target.izoneA`, `target.izoneB`, from contacts at 10A) and ligand
   RMSD zones (`target.lzone`) for every complex (`ana_scripts`, from `XXX_target.pdb`) and every
   symmetry-related reference (`ana_scripts_K`, from their `target.pdb`). It calls:
   * `make-ana-files.sh TARGET.pdb OUTDIR` : the same files for one reference complex
     (receptor chain A first, ligand chain B)
   * `join_zones.py` : merges consecutive residues into zones (HADDOCK analysis tools, Python 3 version)
3. `make-ambig.sh` : true-interface ambiguous restraints `XXX_ambig.tbl` (3.9A) and `XXX_ambig5.tbl` (5A),
   with `haddock-restraints ti` on `XXX_target.pdb` without `HETATM`
4. `make-hbonds.sh` : `XXX_hbonds.tbl`, `haddock-restraints restraint` on `XXX_r_u.pdb` and `XXX_l_u.pdb`
   concatenated (restraints keeping chain breaks, merged chains and co-factors together); no file when
   there is nothing to restrain
5. `make-info.sh FILE.pdb...` : the `.info` statistics files (pdb_wc and whether the numbering has gaps)
6. `make-irmsd-tables.sh` : the i-RMSD tables `i-RMSD.dat`, `i-RMSD_r.dat`, `i-RMSD_l.dat` and their
   `-sorted` versions, with
   * `irmsd.sh REFERENCE.pdb MOBILE.pdb ZONEFILE` : ProFit backbone (N, CA, C, O) RMSD over a zone file

haddock3 benchmark lists:

* `make-haddock3-lists.py` : writes `haddock3-scripts/input.list` (all complexes except the AB ones) and one
  `haddock3-scripts/input-<NAME>.list` per `HADDOCK-ready/categories/<NAME>.txt` (haddock-runner format)

Symmetry-related references (`ana_scripts_K`):

* `rebuild-symmetry-refs.py OUT COMPLEX...` : rebuilds `ana_scripts_K/target.pdb` of the given complexes
  from the current `XXX_target.pdb`. One partner is placed on the K-th equivalent binding site by a
  symmetry operation of the other partner, computed from the reference itself (superposition of its
  copies, on their rigid core, or on the binding site when the copies differ by domain motions). The
  operation reproducing the existing reference best is chosen (so that `ana_scripts_K` keeps its site).
  Copy the new `target.pdb` into place, then run `make-ana-files-all.sh`.

* `make-symmetry-refs.py OUT COMPLEX...` : creates the symmetry-related references of complexes that have
  none yet: each symmetry operation of a partner made of symmetric copies is applied to the other partner,
  and kept when it gives a distinct binding site without clashes and with an interface of comparable size.

Some PDB files hold other data than the element and charge in columns 77-80: `make-ambig.sh` and
`make-hbonds.sh` rebuild these columns in the copies read by `haddock-restraints`.


## curation

Scripts used to check the `HADDOCK-ready` structures against the PDB entries and correct them
(missing residues, selenomethionines, modified residues, numbering, `TER` records). They document
how the current files were produced; they were run once and their results checked by hand.

Requirements: Python 3 with numpy and biopython. Run them from a working directory, which holds the
downloaded entries and the intermediate results. The BM5 table entries are in `unbound_ids.json`
(unbound receptor and ligand per complex, `PDB_chains`) and `bound_ids.json`
(`PDB_receptorchains:ligandchains`), next to the scripts.

Missing residues:

1. `fetch_entries.sh` : downloads the unbound (`cif/`) and bound (`cifb/`) entries as mmCIF from RCSB
2. `gaps.py` : maps each unbound file onto its PDB entry (superposition seeded by sequence alignment,
   CA within 0.1A; symmetry copies and NMR models handled) and lists the internal gaps (missing termini
   ignored) &rarr; `gaps.json`
3. `bound_check.py` : classifies the residues of each gap: not modelled in the unbound entry, modelled in
   the unbound but not the bound entry, or modelled in both &rarr; `gap_report.json`
4. `lost_why.py` : for residues modelled in both entries, why they were lost (selenomethionine, modified
   residue, alternate locations, insertion code, type differs between bound and unbound...) &rarr; `lost.json`
5. `make_tables.py` : `lost_residues.tsv` (residues modelled in both entries but missing from both files)
   and `crystal_gaps.tsv` (all gaps, per class)

Corrections (each writes the corrected files to an output directory, to be checked and copied into
`HADDOCK-ready`; the derived files are then regenerated with the scripts above):

* `renumber.py` : renumbers bound partners whose numbering differed from the unbound one (1OPH ligand,
  1FAK and 1AZS receptors)
* `fix_mse.py OUT` : adds the selenomethionines dropped from both files, as MET
* `fix_mod.py OUT` : adds dropped modified residues as their parent residue (methyl-His, Cys
  derivatives, seleno-Ser)
* `fix_alt.py OUT PARTNER... [--all-same-type] [--unbound-only]` : adds residues dropped because of
  alternate locations (or, with `--all-same-type`, any residue of the same type in both entries, e.g.
  insertion codes), renumbering the later residues and `XXX_hbonds.tbl`; writes `renumbering.json`
* `fix_bound.py NEWDIR PARTNER` : completes the bound file and the reference complex to match a new
  unbound file (used with `fix_alt.py --unbound-only` when the bound entry numbers residues differently)
* `swap_protomers.py FILE CHAIN N` : swaps the numbering of protomers 2 and 3 (4FQI ligand trimer)
* `fix_se.py FILE...` : remaining selenium in MET (atom `SE` renamed `SD`, element `S`)

Adding the Docking Benchmark 5.5 entries:

* `add_bm55.py STRUCTURES TABLE.xlsx OUTDIR COMPLEX...` : converts BM5.5 entries (`benchmark5.5/structures`
  and `Table_BM5.5.xlsx` from [https://zlab.wenglab.org/benchmark/](https://zlab.wenglab.org/benchmark/)) to the
  HADDOCK-ready layout: bound and unbound chains paired (sequence identity and CA positions), only residues
  present in both forms with the same type kept, numbering of the unbound chains without insertion codes,
  chains merged into one chain per partner (A receptor, B ligand; later chains shifted by multiples of 500),
  and the reference complex `XXX_target.pdb`. The derived files are then made with the scripts above
  (`make-info.sh`, `make-target-unbound.py`, `make-ana-files.sh`, `make-ambig.sh`, `make-hbonds.sh`,
  `make-irmsd-tables.sh`).

`TER` records:

* `ter_chains.py` : origin (PDB entry chain and symmetry copy) of every residue of the unbound and bound
  files &rarr; `origins.json`
* `apply_ters.py [--write]` : removes `TER` records at gaps within a chain and keeps those between
  different chains (dry run: `ter_decisions.tsv`)
