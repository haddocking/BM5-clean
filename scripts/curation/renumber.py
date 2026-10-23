"""Renumber bound partners of BM5-clean HADDOCK-ready to the unbound numbering.

- 1OPH ligand: the bound file skips 49; mapped by sequence alignment (identical sequences).
- 1FAK receptor: insertion codes in the bound file; mapped by sequence alignment.
- 1AZS receptor: the bound domain superimposes on the second copy of the unbound
  homodimer (numbered +1000), not on the first copy whose numbers it used.
The b-matched file and the same chain of XXXX_target.pdb are rewritten (columns 23-27).
"""
from pathlib import Path
from Bio.Align import PairwiseAligner
from Bio.SeqUtils import seq1

D = Path(os.environ.get("BM5_HADDOCK_READY", Path(__file__).resolve().parents[2] / "HADDOCK-ready"))
HERE = Path(__file__).resolve().parent  # BM5 table IDs (unbound_ids.json, bound_ids.json)
FIXES = {("1OPH", "l", "B"): "align", ("1FAK", "r", "A"): "align", ("1AZS", "r", "A"): 1000}
al = PairwiseAligner(mode="global", match_score=2, mismatch_score=-1, open_gap_score=-4, extend_gap_score=-0.5)
al.target_end_gap_score = al.query_end_gap_score = 0


def residues(pdb):
    out = []
    for line in Path(pdb).read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            k = (line[22:27], line[17:20].strip())
            if not out or out[-1] != k:
                out.append(k)
    return out


def mapping(bound, unbound, how):
    if how == "align":
        s = lambda rs: "".join(seq1(r[1]) or "X" for r in rs)
        a = al.align(s(bound), s(unbound))[0]
        m = {}
        for (b0, b1), (u0, u1) in zip(*a.aligned):
            for i, j in zip(range(b0, b1), range(u0, u1)):
                m[bound[i][0]] = int(unbound[j][0][:4])
        assert len(m) == len(bound), "unmatched bound residues"
    else:
        m = {k: int(k[:4]) + how for k, _ in bound}
        assert all(k[4] == " " for k, _ in bound)
    new = [m[k] for k, _ in bound]
    assert new == sorted(new) and len(set(new)) == len(new), "mapping not increasing"
    utype = {int(k[:4]): t for k, t in unbound if k[4] == " "}
    diff = [(k.strip(), t, m[k], utype[m[k]]) for k, t in bound if m[k] in utype and utype[m[k]] != t]
    return m, diff


def rewrite(pdb, m, chain=None):
    lines = []
    for line in Path(pdb).read_text().splitlines(keepends=True):
        if line.startswith(("ATOM", "HETATM", "ANISOU", "TER")) and len(line) > 27 \
                and (chain is None or line[21] == chain) and line[22:27] in m:
            line = f"{line[:22]}{m[line[22:27]]:>4} {line[27:]}"
        lines.append(line)
    Path(pdb).write_text("".join(lines))


for (n, side, chain), how in FIXES.items():
    b, u, t = D / n / f"{n}_{side}_b-matched.pdb", D / n / f"{n}_{side}_u.pdb", D / n / f"{n}_target.pdb"
    bound = residues(b)
    m, diff = mapping(bound, residues(u), how)
    changed = sum(1 for k, _ in bound if k.strip() != str(m[k]))
    rewrite(b, m)
    rewrite(t, m, chain)
    print(f"{n}_{side}: {changed} of {len(bound)} residues renumbered; type differences {diff}")
