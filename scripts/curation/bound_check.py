"""Classify the internal gaps of the BM5-clean unbound files (gaps.json) against the
unbound PDB entry, the bound PDB entry and the bound b-matched file.

For each residue missing between two residues of the same chain in an unbound file:
- unmodeled: not modelled in the unbound PDB entry (a gap in the deposited structure);
- not_in_bound: modelled in the unbound entry, but not in the bound entry (removed so
  that bound and unbound match);
- in_both: modelled in both entries (lost when preparing the files: candidate to add);
- no_bound_match: the unbound chain has no counterpart among the bound chains.
"""
import os
import json
import re
import warnings
from pathlib import Path

from Bio.Align import PairwiseAligner
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.SeqUtils import seq1

warnings.filterwarnings("ignore")
D = Path(os.environ.get("BM5_HADDOCK_READY", Path(__file__).resolve().parents[2] / "HADDOCK-ready"))
HERE = Path(__file__).resolve().parent  # BM5 table IDs (unbound_ids.json, bound_ids.json)
gaps = json.load(open("gaps.json"))
unbound_ids = json.load(open(HERE / "unbound_ids.json"))
bound_ids = json.load(open(HERE / "bound_ids.json"))
aligner = PairwiseAligner(mode="global", match_score=2, mismatch_score=-1,
                          open_gap_score=-4, extend_gap_score=-0.5)
aligner.target_end_gap_score = aligner.query_end_gap_score = 0
_cache = {}


def scheme(cif):
    """Per author chain: SEQRES residues with author number and whether modelled."""
    if cif in _cache:
        return _cache[cif]
    d = MMCIF2Dict(str(cif))
    out = {}
    keys = ["seq_id", "mon_id", "pdb_seq_num", "pdb_ins_code", "pdb_strand_id", "pdb_mon_id", "auth_seq_num"]
    for seq, mon, num, ins, strand, pmon, auth in zip(*[d[f"_pdbx_poly_seq_scheme.{k}"] for k in keys]):
        lst = out.setdefault(strand, [])
        if lst and lst[-1]["seq"] == int(seq):
            continue
        lst.append({"seq": int(seq), "mon": mon, "id": num + ("" if ins in ".?" else ins),
                    "modeled": pmon != "?" and auth != "?"})
    _cache[cif] = out
    return out


def chains_of(spec_str):
    m = re.match(r"(\w{4})_?([A-Za-z0-9]*)", spec_str)
    return m.group(1), list(m.group(2))


def residues(pdb):
    out = []
    for line in Path(pdb).read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            k = line[22:27]
            if not out or out[-1] != k:
                out.append(k)
    return out


def num(key):
    return int(re.match(r"-?\d+", key.strip()).group())


report = []
for key, v in sorted(gaps.items()):
    if "gaps" not in v:
        continue
    name, side = key.rsplit("_", 1)
    entry_side = ("l" if side == "r" else "r") if v.get("swapped") else side
    upid = chains_of(v["pdb"])[0]
    usch = scheme(Path("cif") / f"{upid}.cif")
    bpid, rest = bound_ids[name].split("_", 1)
    rchains, lchains = rest.split(":")
    bchains = list(rchains if entry_side == "r" else lchains)
    bsch = scheme(Path("cifb") / f"{bpid}.cif")
    bmatched = [num(k) for k in residues(D / name / f"{name}_{side}_b-matched.pdb")]
    ufile = [k.strip() for k in residues(D / name / f"{name}_{side}_u.pdb")]
    for g in v["gaps"]:
        if g.get("order"):
            continue
        ch = g["chain"]
        lst = usch[ch]
        ids = [x["id"] for x in lst]
        p1, p2 = ids.index(g["pdb_from"]), ids.index(g["pdb_to"])
        between = lst[p1 + 1:p2]
        # counterpart bound chain: best SEQRES alignment
        useq = "".join(seq1(x["mon"]) or "X" for x in lst)
        best = None
        for bc in bchains:
            if bc not in bsch:
                continue
            bseq = "".join(seq1(x["mon"]) or "X" for x in bsch[bc])
            aln = aligner.align(useq, bseq)[0]
            if best is None or aln.score > best[0]:
                best = (aln.score, bc, aln)
        umap = {}
        if best:
            for (a0, a1), (b0, b1) in zip(*best[2].aligned):
                for a, b in zip(range(a0, a1), range(b0, b1)):
                    umap[a] = b
        residues_info = []
        for k, x in enumerate(between, start=p1 + 1):
            if not x["modeled"]:
                cls, bres = "unmodeled", None
            elif k not in umap:
                cls, bres = "no_bound_match", None
            else:
                bx = bsch[best[1]][umap[k]]
                bres = f"{best[1]}:{bx['mon']}{bx['id']}"
                cls = "in_both" if bx["modeled"] else "not_in_bound"
            residues_info.append({"unbound": f"{ch}:{x['mon']}{x['id']}", "bound": bres, "class": cls})
        lo, hi = num(g["after"]), num(g["before"])
        i1, i2 = ufile.index(g["after"]), ufile.index(g["before"])
        in_file = i2 - i1 - 1  # residues in the file between the two (present but not mapped)
        report.append({
            "partner": key, "unbound_pdb": v["pdb"], "bound_pdb": f"{bpid}_{''.join(bchains)}",
            "bound_chain": best[1] if best else None, "swapped": v.get("swapped", False),
            "file_after": g["after"], "file_before": g["before"],
            "bmatched_in_gap": [n for n in bmatched if lo < n < hi],
            "n": len(between), "file_residues_in_gap": in_file, "missing": len(between) - in_file,
            "residues": residues_info,
            "classes": {c: sum(r["class"] == c for r in residues_info)
                        for c in ("unmodeled", "not_in_bound", "in_both", "no_bound_match")},
        })
json.dump(report, open("gap_report.json", "w"), indent=1)
tot = {c: sum(r["classes"][c] for r in report) for c in ("unmodeled", "not_in_bound", "in_both", "no_bound_match")}
print(len(report), "gaps;", tot)
print("gaps with residues in both entries:", sum(1 for r in report if r["classes"]["in_both"]))
