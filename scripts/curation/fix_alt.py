"""Add back the residues dropped from BM5-clean because of alternate locations (report section C).

python fix_alt.py OUT PARTNER... [--all-same-type] [--unbound-only]   (e.g. 1JMO_l 1OPH_l 2VXT_r)

With --all-same-type, every residue of the partner present in both entries with the same
type is added, whatever the reason it was lost (altlocs, insertion codes, unexplained).
With --unbound-only, only the unbound file (and hbonds) is written; the bound side is then
completed with fix_bound.py (needed when the bound entry numbers residues differently, e.g.
repeated residue types, or the bound file is built from symmetry copies).
LOST_TSV=file.tsv replaces lost_residues.tsv (e.g. after an earlier renumbering of the partner).

For each residue in lost_residues.tsv lost because of altlocs (same type in both entries):
- unbound file: the residue of the unbound entry (altloc A or the only conformation);
- b-matched and _target.pdb: the residue of the bound entry;
placed by a fit on the neighbouring residues (matched by order and residue name).
The files number residues sequentially, so the partner's later residues are renumbered
(+1 per inserted residue before them), in the unbound file, b-matched, _target.pdb and
XXXX_hbonds.tbl. Writes the new files to OUT/<complex>/ and the renumbering
(old -> new, per partner) to OUT/renumbering.json.
"""
import csv
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
src = (HERE / "fix_mse.py").read_text()
exec(src[:src.index("def transform(")])  # kabsch, atoms, pick, file_residues, G, D
exec(src[src.index("def fmt_name("):src.index("def new_lines(")])
exec(src[src.index("def renumber_serials("):src.index("def main(")])

OUT = Path(sys.argv[1])
ALL = "--all-same-type" in sys.argv
UNBOUND_ONLY = "--unbound-only" in sys.argv
PARTNERS = [a for a in sys.argv[2:] if not a.startswith("--")]
WINDOW = 4


def parse_res(s):
    ch, name, num = re.match(r"(\w):([A-Z0-9]{3})(-?\d+[A-Z]?)$", s).groups()
    return ch, name, num


def entry_order(ent, chain):
    """Residue keys of one chain of an entry, in file order."""
    return [k for k in ent if k[0] == chain]


def blocks(lines, chain):
    """[(residue number, [line indices])] for ATOM/HETATM lines of chain, in order."""
    out = []
    for i, l in enumerate(lines):
        if l.startswith(("ATOM", "HETATM")) and l[21] == chain:
            n = int(l[22:26])
            if not out or out[-1][0] != n:
                out.append((n, []))
            out[-1][1].append(i)
    return out


def fit_group(lines, chain, after, before, ent, ech, first, last):
    """Superimpose entry -> file using up to WINDOW residues on each side of the gap,
    walking outwards in both and stopping at a residue-name mismatch."""
    blk = blocks(lines, chain)
    nums = [n for n, _ in blk]
    order = entry_order(ent, ech)
    ia, ib = nums.index(after), nums.index(before)
    ea, eb = order.index(first) - 1, order.index(last) + 1
    pairs = []
    for step, (fi, ei) in [(-1, (ia, ea)), (1, (ib, eb))]:
        for k in range(WINDOW):
            f, e = fi + step * k, ei + step * k
            if not (0 <= f < len(blk) and 0 <= e < len(order)):
                break
            fl = [lines[i] for i in blk[f][1]]
            if fl[0][17:20] != ent[order[e]]["_comp"]:
                break
            fat = {l[12:16].strip(): np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]) for l in fl}
            for a in ("N", "CA", "C", "O"):
                if a in fat and a in ent[order[e]]:
                    pairs.append((pick(ent[order[e]][a])[3], fat[a]))
    if len(pairs) < 3:
        return np.eye(3), np.zeros(3), float("inf"), len(pairs)
    X, Y = np.array([p[0] for p in pairs]), np.array([p[1] for p in pairs])
    R, t = kabsch(X, Y)
    rmsd = float(np.sqrt(((X @ R.T + t - Y) ** 2).sum(1).mean()))
    return R, t, rmsd, len(pairs)


def residue_lines(template, res, R, t):
    out = []
    for name, confs in res.items():
        if name == "_comp" or name.startswith("H") or name in ("D",) or name.startswith("D") and len(name) > 1 and name[1].isdigit():
            continue
        alt, occ, b, xyz = pick(confs)
        x, y, z = xyz @ R.T + t
        l = template.rstrip("\n").ljust(80)
        el = name[0]
        l = ("ATOM  " + NEW + l[11:12] + fmt_name(name) + " " + res["_comp"] + l[20:22] + "{NUM}" + " " + l[27:30]
             + f"{x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{b:6.2f}" + l[66:76] + f"{el:>2}" + "  ")
        out.append(l.rstrip() + "\n")
    return out


def main():
    gaps = json.load(open(G / "gaps.json"))
    bids = json.load(open(HERE / "bound_ids.json"))
    import os
    rows = [r for r in csv.DictReader(open(os.environ.get("LOST_TSV", G / "lost_residues.tsv")), delimiter="\t")
            if r["partner"] in PARTNERS and r["same_type"] == "True"
            and not any(x in r["reasons"] for x in ("selenomethionine", "modified", "type differs"))
            and (ALL or ("altloc" in r["reasons"] and "insertion" not in r["reasons"]))]
    renumbering, log = {}, []
    for partner in PARTNERS:
        name, side = partner.rsplit("_", 1)
        tch = "A" if side == "r" else "B"
        upid = gaps[partner]["pdb"][:4]
        bpid = bids[name].split("_")[0]
        uent, bent = atoms(G / "cif" / f"{upid}.cif"), atoms(G / "cifb" / f"{bpid}.cif")
        prow = [r for r in rows if r["partner"] == partner]
        groups = defaultdict(list)
        for r in prow:
            groups[(int(r["file_after"]), int(r["file_before"]))].append(r)
        paths = {"u": D / name / f"{name}_{side}_u.pdb", "b": D / name / f"{name}_{side}_b-matched.pdb",
                 "t": D / name / f"{name}_target.pdb"}
        new_files = {}
        for kind, path in paths.items():
            if UNBOUND_ONLY and kind != "u":
                continue
            lines = path.read_text().splitlines(keepends=True)
            chain = tch if kind == "t" else next(l[21] for l in lines if l.startswith("ATOM"))
            ent = uent if kind == "u" else bent
            inserts = {}  # after -> list of residue line lists
            for (after, before), grp in sorted(groups.items()):
                col = "unbound_residue" if kind == "u" else "bound_residue"
                keys = []
                for r in grp:
                    ech, rname, num = parse_res(r[col])
                    keys.append((ech, num))
                    assert ent[(ech, num)]["_comp"] == rname, (r, ent[(ech, num)]["_comp"])
                order = entry_order(ent, keys[0][0])
                keys.sort(key=order.index)
                R, t, rmsd, n = fit_group(lines, chain, after, before, ent, keys[0][0], keys[0], keys[-1])
                shifted = 0
                if rmsd >= 0.05 and kind != "u":
                    # repeated residue types (e.g. Asn-Asn) can make the table pick the neighbouring
                    # copy, which then sits on the other side of a file residue: try the same residue
                    # types one position along the entry, and neighbouring insertion points
                    nums = [m for m, _ in blocks(lines, chain)]
                    ia = nums.index(after)
                    found = None
                    for sh in (1, -1):
                        idx = [order.index(k) + sh for k in keys]
                        if min(idx) < 0 or max(idx) >= len(order):
                            continue
                        alt_keys = [order[i] for i in idx]
                        if [ent[k]["_comp"] for k in alt_keys] != [ent[k]["_comp"] for k in keys]:
                            continue
                        for da in (0, 1, -1, 2, -2):
                            j = ia + da
                            if not (0 <= j < len(nums) - 1):
                                continue
                            R2, t2, r2, n2 = fit_group(lines, chain, nums[j], nums[j + 1], ent,
                                                       alt_keys[0][0], alt_keys[0], alt_keys[-1])
                            if r2 < 0.05:
                                found = (alt_keys, R2, t2, r2, n2, sh, nums[j], nums[j + 1])
                                break
                        if found:
                            break
                    if found:
                        keys, R, t, rmsd, n, shifted, after, before = found
                assert rmsd < 0.05, f"{path.name} {after}-{before}: rmsd {rmsd:.3f} on {n} atoms"
                blk = dict(blocks(lines, chain))
                template = lines[blk[after][-1]]
                inserts.setdefault(after, []).extend((order.index(k), residue_lines(template, ent[k], R, t)) for k in keys)
                log.append({"partner": partner, "file": kind, "after": after, "before": before,
                            "residues": [f"{k[0]}:{ent[k]['_comp']}{k[1]}" for k in keys],
                            "rmsd": round(rmsd, 3), "fit_atoms": n, "shifted": shifted})
            # rebuild the chain with insertions and new numbers
            out, shift, mapping = [], 0, {}
            blk = blocks(lines, chain)
            last_of = {idx[-1]: num for num, idx in blk}
            inserts = {a: [r for _, r in sorted(v, key=lambda x: x[0])] for a, v in inserts.items()}
            for i, l in enumerate(lines):
                if l.startswith(("ATOM", "HETATM")) and l[21] == chain:
                    old = int(l[22:26])
                    mapping[old] = old + shift
                    l = l[:22] + f"{old + shift:4d}" + l[26:]
                out.append(l)
                if i in last_of and last_of[i] in inserts:
                    base = last_of[i] + shift
                    for k, res in enumerate(inserts[last_of[i]], start=1):
                        out += [x.replace("{NUM}", f"{base + k:4d}") for x in res]
                    # free numbers between the neighbours are used first; only the excess shifts
                    nums_ = [m for m, _ in blk]
                    j_ = nums_.index(last_of[i])
                    room = nums_[j_ + 1] - last_of[i] - 1 if j_ + 1 < len(nums_) else len(inserts[last_of[i]])
                    shift += max(0, len(inserts[last_of[i]]) - room)
            if kind == "u":
                renumbering[partner] = mapping
                useq = [(int(l[22:26]), l[17:20]) for l in out if l.startswith(("ATOM", "HETATM")) and l[21] == chain]
                useq = [x for i, x in enumerate(useq) if i == 0 or x != useq[i - 1]]
            else:
                # bound residues take the unbound numbers, position by position (same residue types)
                bres = []
                for i, l in enumerate(out):
                    if l.startswith(("ATOM", "HETATM")) and l[21] == chain:
                        k = (l[22:26], l[17:20])
                        if not bres or bres[-1][0] != k:
                            bres.append((k, []))
                        bres[-1][1].append(i)
                assert len(bres) == len(useq), f"{path.name}: {len(bres)} residues, unbound {len(useq)}"
                for (k, idx), (num, name3) in zip(bres, useq):
                    assert k[1] == name3, f"{path.name}: {k} vs unbound {num}{name3}"
                    for i in idx:
                        out[i] = out[i][:22] + f"{num:4d}" + out[i][26:]
            out, _ = renumber_serials(out)
            new_files[path] = out
        for path, out in new_files.items():
            dst = OUT / name / path.name
            dst.parent.mkdir(parents=True, exist_ok=True)
            dst.write_text("".join(out))
        # hbonds.tbl: residue numbers of this partner's segid
        hb = D / name / f"{name}_hbonds.tbl"
        dst = OUT / name / hb.name
        text = (dst if dst.exists() else hb).read_text() if hb.exists() else None
        if text is not None:
            m = renumbering[partner]
            text = re.sub(rf"(segid {tch} and resi )(-?\d+)",
                          lambda x: x.group(1) + str(m.get(int(x.group(2)), int(x.group(2)))), text)
            dst.write_text(text)
    json.dump({p: {str(k): v for k, v in m.items() if k != v} for p, m in renumbering.items()},
              open(OUT / "renumbering.json", "w"), indent=1)
    json.dump(log, open(OUT / "added.json", "w"), indent=1)
    for e in log:
        print(e["partner"], e["file"], f"{e['after']}-{e['before']}", e["residues"], e["rmsd"], e["fit_atoms"])


if __name__ == "__main__":
    main()
