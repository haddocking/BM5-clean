"""Complete a bound partner (b-matched and its chain in _target.pdb) to match a new unbound file.

python fix_bound.py NEWDIR PARTNER

NEWDIR/<complex>/<partner>_u.pdb is the unbound file with residues added (fix_alt.py, unbound
side). The bound entry is superimposed on the bound file and each file residue mapped to its entry
residue (CA within 0.1 A); copies built from symmetry mates are fitted one at a time. Entry residues
lying between two consecutive mapped file residues of the same entry chain, but absent from the
file, are candidates. The new unbound sequence is aligned to the bound residues + candidates:
every unbound residue must be matched, every bound file residue must be used, residue types must
be identical; candidates may be skipped. Among equivalent alignments, bound file residues keep
their old unbound partner where possible. Matched candidates are added (altloc A), and every
bound residue takes the number of its unbound partner.
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
src = (HERE / "fix_mse.py").read_text()
exec(src[:src.index("def transform(")])
exec(src[src.index("def fmt_name("):src.index("def new_lines(")])
exec(src[src.index("def renumber_serials("):src.index("def main(")])
gsrc = (Path(__file__).resolve().parent / "gaps.py").read_text()
exec(gsrc[gsrc.index("def fit("):gsrc.index("def analyse(")])
from Bio.Align import PairwiseAligner  # noqa: E402
from Bio.SeqUtils import seq1  # noqa: E402
aligner = PairwiseAligner(mode="global", match_score=2, mismatch_score=-1, open_gap_score=-4, extend_gap_score=-0.5)
aligner.target_end_gap_score = aligner.query_end_gap_score = 0

NEWDIR, PARTNER = Path(sys.argv[1]), sys.argv[2]
name, side = PARTNER.rsplit("_", 1)
tch = "A" if side == "r" else "B"
bids = json.load(open(HERE / "bound_ids.json"))
bpid, rest = bids[name].split("_", 1)
ent = atoms(G / "cifb" / f"{bpid}.cif")


def residues_of(lines, chain):
    """[(number, name, [line indices])] in file order."""
    out = []
    for i, l in enumerate(lines):
        if l.startswith(("ATOM", "HETATM")) and l[21] == chain:
            k = int(l[22:26])
            if not out or out[-1][0] != k:
                out.append((k, l[17:20], []))
            out[-1][2].append(i)
    return out


def map_to_entry(lines, chain):
    """Per copy rigid transforms entry -> file, {file residue index: entry key} and
    {file residue index: copy}. Copies of the entry (symmetry mates built into the file) are
    fitted one at a time, as in gaps.py."""
    res = residues_of(lines, chain)
    fres = [{"name": n, "atoms": {lines[i][12:16].strip(): np.array([float(lines[i][30:38]), float(lines[i][38:46]),
                                                                  float(lines[i][46:54])]) for i in idx}}
            for _, n, idx in res]
    fca_all = [(i, r["atoms"]["CA"]) for i, r in enumerate(fres) if "CA" in r["atoms"]]
    mapping, copy_of, transforms = {}, {}, []
    for copy in range(8):
        fca = [(i, c) for i, c in fca_all if i not in mapping]
        if len(fca) < 10:
            break
        best = None
        for chains in rest.split(":"):
            keys = [k for k in ent if k[0] in chains and "CA" in ent[k]]
            pdb_ca = [(k[0], k[1], ent[k]["_comp"], pick(ent[k]["CA"])[3]) for k in keys]
            coords = np.array([x[3] for x in pdb_ca])
            kk, ok = fit(fca, fres, pdb_ca, coords, list(chains))
            if kk is None:
                continue
            if best is None or ok.sum() > best[0]:
                best = (int(ok.sum()), {fca[j][0]: (pdb_ca[kk[j]][0], pdb_ca[kk[j]][1]) for j in range(len(fca)) if ok[j]})
        if best is None or best[0] < 10:
            break
        m = best[1]
        X = np.array([pick(ent[m[i]]["CA"])[3] for i in m]); Y = np.array([fres[i]["atoms"]["CA"] for i in m])
        R, t = kabsch(X, Y)
        rmsd = float(np.sqrt(((X @ R.T + t - Y) ** 2).sum(1).mean()))
        transforms.append((R, t, rmsd, len(m)))
        for i, k in m.items():
            mapping[i] = k; copy_of[i] = copy
    return res, mapping, copy_of, transforms


def candidates(res, mapping, copy_of):
    """Sequence of ('file', residue index) and ('cand', entry key) items."""
    items = []
    order = {c: [k for k in ent if k[0] == c] for c in {v[0] for v in mapping.values()}}
    prev = None
    for i in range(len(res)):
        if i in mapping and prev is not None and mapping[prev][0] == mapping[i][0] and copy_of[prev] == copy_of[i]:
            o = order[mapping[i][0]]
            a, b = o.index(mapping[prev]), o.index(mapping[i])
            items += [("cand", k, copy_of[i]) for k in o[a + 1:b] if "CA" in ent[k]]
        items.append(("file", i))
        if i in mapping:
            prev = i
    return items


def align(useq, items, res, old_u_of):
    """DP: every useq position matched; every 'file' item matched; 'cand' items optional."""
    INF = 10 ** 9
    n, m = len(useq), len(items)
    name_of = lambda it: res[it[1]][1] if it[0] == "file" else ent[it[1]]["_comp"]
    cost = np.full((n + 1, m + 1), INF, dtype=np.int64)
    cost[0, 0] = 0
    back = {}
    for j in range(m + 1):
        for i in range(n + 1):
            if cost[i, j] >= INF:
                continue
            if j < m and items[j][0] == "cand" and cost[i, j] < cost[i, j + 1]:  # skip candidate
                cost[i, j + 1] = cost[i, j]; back[(i, j + 1)] = (i, j, False)
            if i < n and j < m and useq[i][1] == name_of(items[j]):
                c = cost[i, j] + (0 if items[j][0] == "cand" or old_u_of.get(items[j][1]) == useq[i][2] else 1)
                if c < cost[i + 1, j + 1]:
                    cost[i + 1, j + 1] = c; back[(i + 1, j + 1)] = (i, j, True)
    assert cost[n, m] < INF, "no complete alignment"
    pairs, ij = [], (n, m)
    while ij != (0, 0):
        i, j, matched = back[ij]
        if matched:
            pairs.append((i, j))
        ij = (i, j)
    return pairs[::-1], int(cost[n, m])


def main():
    uf = NEWDIR / name / f"{name}_{side}_u.pdb"
    ulines = uf.read_text().splitlines(keepends=True)
    uchain = next(l[21] for l in ulines if l.startswith("ATOM"))
    ren = {int(k): v for k, v in json.load(open(NEWDIR / "renumbering.json")).get(PARTNER, {}).items()}
    # new unbound sequence: (number, name, old number; None for added residues)
    orig_u = residues_of((D / name / uf.name).read_text().splitlines(keepends=True), uchain)
    pre = {ren.get(num, num): num for num, _, _ in orig_u}
    useq = [(num, nm, pre.get(num)) for num, nm, _ in residues_of(ulines, uchain)]
    assert sum(x[2] is not None for x in useq) == len(orig_u)
    report = {}
    for kind, path, chain in (("b", D / name / f"{name}_{side}_b-matched.pdb", None),
                              ("t", D / name / f"{name}_target.pdb", tch)):
        lines = path.read_text().splitlines(keepends=True)
        chain = chain or next(l[21] for l in lines if l.startswith("ATOM"))
        res, mapping, copy_of, transforms = map_to_entry(lines, chain)
        rmsd = max(x[2] for x in transforms)
        assert rmsd < 0.05, f"{path.name}: entry -> file rmsd {rmsd:.3f}"
        items = candidates(res, mapping, copy_of)
        old_u_of = {i: res[i][0] for i in range(len(res))}  # bound and unbound shared numbers before
        pairs, cost = align(useq, items, res, old_u_of)
        # build new chain lines
        newnum = {j: useq[i][0] for i, j in pairs}
        used = {j for _, j in pairs}
        assert all(j in used for j, it in enumerate(items) if it[0] == "file"), "bound residue left unmatched"
        out, added = [], []
        first_chain_line = min(idx[0] for _, _, idx in res)
        last_chain_line = max(idx[-1] for _, _, idx in res)
        chain_block = []
        for j, it in enumerate(items):
            if it[0] == "file":
                for i in res[it[1]][2]:
                    l = lines[i]
                    chain_block.append(l[:22] + f"{newnum[j]:4d}" + l[26:])
            elif j in used:
                template = lines[res[0][2][0]]
                k = it[1]
                R, t = transforms[it[2]][:2]
                for a, confs in ent[k].items():
                    if a == "_comp" or a.startswith("H"):
                        continue
                    alt, occ, b, xyz = pick(confs)
                    x, y, z = xyz @ R.T + t
                    l = template.rstrip("\n").ljust(80)
                    l = ("ATOM  " + "*****" + l[11:12] + fmt_name(a) + " " + ent[k]["_comp"] + l[20:22] + f"{newnum[j]:4d}"
                         + " " + l[27:30] + f"{x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{b:6.2f}" + l[66:76] + f"{a[0]:>2}" + "  ")
                    chain_block.append(l.rstrip() + "\n")
                added.append(f"{k[0]}:{ent[k]['_comp']}{k[1]}->{newnum[j]}")
        # chain lines are contiguous in these files except TERs; keep non-chain lines in place
        between = [l for l in lines[first_chain_line:last_chain_line + 1]
                   if not (l.startswith(("ATOM", "HETATM")) and l[21] == chain)]
        assert all(l.startswith("TER") for l in between), "chain lines not contiguous"
        out = lines[:first_chain_line] + chain_block + between + lines[last_chain_line + 1:]
        out, _ = renumber_serials(out)
        dst = NEWDIR / name / path.name
        dst.write_text("".join(out))
        report[kind] = {"entry_fit_rmsd": round(rmsd, 3), "copies": len(transforms), "mapped": len(mapping), "residues": len(res),
                        "added": len(added), "kept_old_partner_mismatches": cost, "added_list": added}
        print(path.name, {k: v for k, v in report[kind].items() if k != "added_list"})
    json.dump(report, open(NEWDIR / f"{PARTNER}_bound.json", "w"), indent=1)


if __name__ == "__main__":
    main()
