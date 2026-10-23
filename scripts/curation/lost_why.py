"""Why were residues modelled in both PDB entries lost from the unbound files?"""
import json, collections
from pathlib import Path
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
r = json.load(open("gap_report.json"))
cache = {}
def atoms(pid, folder="cif"):
    if (folder, pid) not in cache:
        d = MMCIF2Dict(f"{folder}/{pid}.cif")
        a = collections.defaultdict(list)
        for grp, name, alt, comp, ch, num, ins, model in zip(d["_atom_site.group_PDB"], d["_atom_site.label_atom_id"],
                d["_atom_site.label_alt_id"], d["_atom_site.auth_comp_id"], d["_atom_site.auth_asym_id"],
                d["_atom_site.auth_seq_id"], d["_atom_site.pdbx_PDB_ins_code"], d["_atom_site.pdbx_PDB_model_num"]):
            if model != d["_atom_site.pdbx_PDB_model_num"][0]: continue
            a[(ch, num + ("" if ins in ".?" else ins))].append((grp, name, alt, comp))
        cache[(folder, pid)] = a
    return cache[(folder, pid)]
STD = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split())
rows = []
for g in r:
    if g["missing"] <= 0 or not g["classes"]["in_both"]: continue
    pid = g["unbound_pdb"][:4]
    for x in g["residues"]:
        if x["class"] != "in_both": continue
        ch, rest = x["unbound"].split(":"); comp, rid = rest[:3], rest[3:]
        at = atoms(pid).get((ch, rid), [])
        names = {a[1] for a in at}; alts = {a[2] for a in at} - {"."}
        het = any(a[0] == "HETATM" for a in at)
        reasons = []
        if rid[-1].isalpha(): reasons.append("insertion code")
        if comp == "MSE": reasons.append("selenomethionine (HETATM MSE)")
        elif comp not in STD: reasons.append(f"modified residue {comp}" + (" (HETATM)" if het else ""))
        if alts and not (alts & {"A"}): reasons.append(f"only altloc {''.join(sorted(alts))}")
        elif alts: reasons.append("altlocs")
        if not {"N", "CA", "C"} <= names: reasons.append("incomplete backbone " + ",".join(sorted(names & {"N","CA","C","O"})))
        # bound side
        bch, brest = x["bound"].split(":"); bcomp, bid = brest[:3], brest[3:]
        bat = atoms(g["bound_pdb"][:4], "cifb").get((bch, bid), [])
        bnames = {a[1] for a in bat}; balts = {a[2] for a in bat} - {"."}
        same = comp == bcomp or {comp, bcomp} <= {"MET", "MSE"}
        if not same: reasons.append(f"type differs ({comp} unbound, {bcomp} bound)")
        if bid[-1].isalpha() and not rid[-1].isalpha(): reasons.append("bound insertion code")
        if bcomp == "MSE" and comp != "MSE": reasons.append("bound selenomethionine")
        elif bcomp not in STD and bcomp != "MSE" and bcomp != comp: reasons.append(f"bound modified {bcomp}")
        if balts: reasons.append("bound altlocs")
        heavy = lambda ns: {n for n in ns if not n.startswith("H")}
        if len(heavy(bnames)) < len(heavy(names)): reasons.append(f"bound incomplete ({len(heavy(bnames))}/{len(heavy(names))} atoms)")
        elif len(heavy(names)) < len(heavy(bnames)): reasons.append(f"unbound incomplete ({len(heavy(names))}/{len(heavy(bnames))} atoms)")
        if not reasons: reasons.append("unexplained")
        x["same_type"] = same
        rows.append({"partner": g["partner"], "unbound": f"{pid} {x['unbound']}", "bound": x["bound"],
                     "file_after": g["file_after"], "file_before": g["file_before"], "reasons": reasons, "same_type": same,
                     "bmatched_in_gap": g["bmatched_in_gap"]})
json.dump(rows, open("lost.json", "w"), indent=1)
c = collections.Counter("; ".join(x["reasons"]) for x in rows)
for k, v in c.most_common(): print(v, k)
print(sorted({x["partner"] for x in rows}))
