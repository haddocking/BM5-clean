"""Write the summary tables of the missing-residue analysis (run in the working directory):

- lost_residues.tsv: residues modelled in both the unbound and the bound PDB entries but missing
  from both files (from lost.json, lost_why.py), with whether the file numbering has room for them;
- crystal_gaps.tsv: the internal gaps of the unbound files (residues really missing from the file) with the number of residues per class
  (from gap_report.json, bound_check.py).
"""
import json
from collections import Counter

lost = json.load(open("lost.json"))
per_gap = Counter((x["partner"], x["file_after"], x["file_before"]) for x in lost)
with open("lost_residues.tsv", "w") as f:
    f.write("partner\tunbound_entry\tunbound_residue\tbound_residue\tfile_after\tfile_before\tsame_type\t"
            "numbering_has_room\treasons\n")
    for x in lost:
        entry, res = x["unbound"].split(" ", 1)
        a, b = x["file_after"], x["file_before"]
        room = int(b) - int(a) - 1 >= per_gap[(x["partner"], a, b)]
        f.write("\t".join([x["partner"], entry, res, x["bound"], a, b, str(x["same_type"]), str(room),
                           "; ".join(x["reasons"])]) + "\n")

report = json.load(open("gap_report.json"))
with open("crystal_gaps.tsv", "w") as f:
    f.write("partner\tunbound_entry\tfile_after\tfile_before\tresidues_in_gap\tnot_modelled_in_unbound_entry\t"
            "modelled_in_unbound_not_in_bound\tmodelled_in_both\n")
    for g in report:
        if g["missing"] <= 0:  # the residues in the gap are in the file, only not mapped (e.g. altlocs)
            continue
        c = g["classes"]
        f.write("\t".join(map(str, [g["partner"], g["unbound_pdb"], g["file_after"], g["file_before"], g["n"],
                                    c["unmodeled"], c["not_in_bound"], c["in_both"]])) + "\n")
