"""Convert Docking Benchmark 5.5 entries (zlab.wenglab.org/benchmark, benchmark5.5/structures) to the
HADDOCK-ready layout of BM5-clean.

python add_bm55.py STRUCTURES TABLE.xlsx OUTDIR COMPLEX...

For each complex and partner (r: receptor, l: ligand), from XXXX_p_u.pdb and XXXX_p_b.pdb (the unbound
structures are superimposed on the bound complex):
- bound chains are paired one-to-one with unbound chains (sequence identity >= 0.8, then the smallest
  CA distances between aligned residues); unbound chains without bound counterpart are dropped;
- each chain pair is aligned and only residues present in both forms with the same type are kept;
- numbering follows the unbound chains (in the order of the BM5.5 table): insertion codes are removed by
  shifting the following residues (gaps keep their size), and each further chain is shifted by a multiple
  of 500 so that numbers do not overlap; all chains are merged into chain A (receptor) or B (ligand),
  with a TER record between the original chains;
- XXXX_target.pdb is the matched bound receptor + ligand.
Writes OUTDIR/XXXX/{XXXX_r_u,XXXX_l_u,XXXX_r_b-matched,XXXX_l_b-matched,XXXX_target}.pdb and
OUTDIR/XXXX/match.json (chain pairing, residues kept and dropped, per chain fit).
"""
import itertools
import json
import re
import sys
from pathlib import Path

import numpy as np
from Bio.Align import PairwiseAligner
from Bio.SeqUtils import seq1

ST, TABLE, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
COMPLEXES = sys.argv[4:]
al = PairwiseAligner(mode="global", match_score=2, mismatch_score=-1, open_gap_score=-4, extend_gap_score=-0.5)
al.target_end_gap_score = al.query_end_gap_score = 0


def table_chains():
    """{complex: {"r": unbound chain order, "l": ...}} from the BM5.5 table."""
    import openpyxl
    out = {}
    for row in openpyxl.load_workbook(TABLE, read_only=True).worksheets[0].iter_rows(values_only=True):
        cells = [c for c in row if c is not None]
        if cells and isinstance(cells[0], str) and re.match(r"^[0-9A-Z]{4}_", cells[0]):
            # unbound receptor and ligand: the next two cells looking like PDBID_chains
            ub = [c for c in cells[1:] if isinstance(c, str) and re.match(r"^[0-9A-Z]{4}_[A-Za-z0-9]*$", c.strip())]
            if len(ub) < 2:
                continue
            out[cells[0][:4]] = {"r": ub[0].split("_")[1], "l": ub[1].split("_")[1], "bound": cells[0],
                                 "unbound": ub[:2]}
    return out


def read(f):
    """{chain: [residue]} with residue = {"key", "name", "lines", "ca"} in file order."""
    chains = {}
    for l in open(f):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        res = chains.setdefault(l[21], [])
        key = l[22:27]
        if not res or res[-1]["key"] != key:
            res.append({"key": key, "name": l[17:20], "lines": [], "ca": None})
        res[-1]["lines"].append(l.rstrip("\n"))
        if l[12:16] == " CA ":
            res[-1]["ca"] = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return chains


def align(b, u):
    """Aligned residue index pairs (bound, unbound)."""
    bs = "".join(seq1(r["name"]) or "X" for r in b); us = "".join(seq1(r["name"]) or "X" for r in u)
    a = al.align(bs, us)[0]
    return [(i, j) for (b0, b1), (u0, u1) in zip(*a.aligned) for i, j in zip(range(b0, b1), range(u0, u1))]


def pair_score(b, u):
    pairs = align(b, u)
    ident = sum(b[i]["name"] == u[j]["name"] for i, j in pairs) / max(1, len(b))
    d = [np.linalg.norm(b[i]["ca"] - u[j]["ca"]) for i, j in pairs if b[i]["ca"] is not None and u[j]["ca"] is not None]
    return ident, float(np.median(d)) if d else 1e9


def pair_chains(bch, uch):
    """One-to-one bound -> unbound chain assignment: identity >= 0.8, minimal summed median CA distance."""
    sc = {(bc, uc): pair_score(bch[bc], uch[uc]) for bc in bch for uc in uch}
    best = None
    for perm in itertools.permutations(uch, len(bch)):
        if any(sc[(bc, uc)][0] < 0.8 for bc, uc in zip(bch, perm)):
            continue
        cost = sum(sc[(bc, uc)][1] for bc, uc in zip(bch, perm))
        if best is None or cost < best[0]:
            best = (cost, dict(zip(bch, perm)))
    assert best, f"no chain pairing: {sc}"
    return best[1], {f"{b}->{u}": [round(sc[(b, u)][0], 3), round(sc[(b, u)][1], 2)] for b, u in best[1].items()}


def numbering(res):
    """New numbers for a chain: insertion codes removed by shifting the following residues."""
    out, shift, prev = [], 0, None
    for r in res:
        num, ins = int(r["key"][:4]), r["key"][4]
        if prev is not None:
            pnum, pins = prev
            if num == pnum:  # insertion code (or repeated number): next number
                shift += 1
            elif num < pnum:  # numbering goes back (should not happen): continue sequentially
                shift += pnum - num + 1
        out.append(num + shift)
        prev = (num, ins)
    return out


def fmt(line, chain, num, serial):
    return f"{line[:6]}{serial:5d}{line[11:21]}{chain}{num:4d} {line[27:72]}{chain:<4}{line[76:]}".rstrip()


def write(path, chain, blocks):
    """blocks: list of [(new number, residue)] per original chain; TER between chains."""
    out, serial = [], 1
    for blk in blocks:
        for num, r in blk:
            for l in r["lines"]:
                out.append(fmt(l.ljust(80), chain, num, serial)); serial += 1
        last = out[-1]
        out.append(f"TER   {serial:5d}      {last[17:27]}".rstrip()); serial += 1
    out.append("END")
    path.write_text("\n".join(out) + "\n")


def kabsch_rmsd(X, Y):
    cx, cy = X.mean(0), Y.mean(0)
    U, S, Vt = np.linalg.svd((X - cx).T @ (Y - cy))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return float(np.sqrt((((X - cx) @ R.T - (Y - cy)) ** 2).sum(1).mean()))


def main():
    tab = table_chains()
    for c in COMPLEXES:
        o = OUT / c; o.mkdir(parents=True, exist_ok=True)
        info = {"table": tab[c]["bound"]}
        for p, cid in (("r", "A"), ("l", "B")):
            u = read(ST / f"{c}_{p}_u.pdb"); b = read(ST / f"{c}_{p}_b.pdb")
            pairing, scores = pair_chains(b, u)
            order = [ch for ch in tab[c][p] if ch in pairing.values()] or list(pairing.values())
            order += [ch for ch in pairing.values() if ch not in order]
            inv = {uc: bc for bc, uc in pairing.items()}
            ublocks, bblocks, maxnum, chain_info = [], [], None, {}
            for uc in order:
                bc = inv[uc]; ures, bres = u[uc], b[bc]
                pairs = [(i, j) for i, j in align(bres, ures) if bres[i]["name"] == ures[j]["name"]]
                keep_u = [j for _, j in pairs]
                nums = numbering([ures[j] for j in keep_u])
                offset = 0
                if maxnum is not None and nums[0] <= maxnum:
                    offset = 500 * ((maxnum - nums[0]) // 500 + 1)
                nums = [n + offset for n in nums]
                maxnum = max(nums) if maxnum is None else max(maxnum, max(nums))
                ublocks.append([(n, ures[j]) for n, (_, j) in zip(nums, pairs)])
                bblocks.append([(n, bres[i]) for n, (i, _) in zip(nums, pairs)])
                ca = [(bres[i]["ca"], ures[j]["ca"]) for i, j in pairs if bres[i]["ca"] is not None and ures[j]["ca"] is not None]
                chain_info[f"{bc}->{uc}"] = {
                    "kept": len(pairs), "bound": len(bres), "unbound": len(ures),
                    "type_differs": sum(1 for i, j in align(bres, ures) if bres[i]["name"] != ures[j]["name"]),
                    "offset": offset, "numbers": [min(nums), max(nums)],
                    "ca_rmsd_after_fit": round(kabsch_rmsd(np.array([x for x, _ in ca]), np.array([y for _, y in ca])), 2)}
            assert maxnum <= 9999, f"{c}_{p}: numbering above 9999"
            write(o / f"{c}_{p}_u.pdb", cid, ublocks)
            write(o / f"{c}_{p}_b-matched.pdb", cid, bblocks)
            info[p] = {"pairing": scores, "unbound_dropped": [ch for ch in u if ch not in pairing.values()],
                       "chains": chain_info}
        # reference complex: matched bound receptor + ligand
        lines = [l for l in (o / f"{c}_r_b-matched.pdb").read_text().splitlines() if l != "END"]
        lines += [l for l in (o / f"{c}_l_b-matched.pdb").read_text().splitlines() if l != "END"]
        n = 1
        for i, l in enumerate(lines):
            lines[i] = l[:6] + f"{n:5d}" + l[11:]; n += 1
        (o / f"{c}_target.pdb").write_text("\n".join(lines) + "\nEND\n")
        json.dump(info, open(o / "match.json", "w"), indent=1)
        print(c, {p: {k: (v["kept"], v["bound"], v["unbound"], v["type_differs"], v["offset"], v["ca_rmsd_after_fit"])
                      for k, v in info[p]["chains"].items()} for p in "rl"},
              "dropped:", {p: info[p]["unbound_dropped"] for p in "rl" if info[p]["unbound_dropped"]})


if __name__ == "__main__":
    main()
