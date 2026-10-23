"""Create the symmetry-related references (ana_scripts_K/target.pdb) of complexes that have none yet.

python make-symmetry-refs.py OUT COMPLEX...   writes OUT/<complex>/ana_scripts_K/target.pdb (K = 2, 3, ...)

For each partner Y that consists of symmetric copies (repeated residue numbering, as in the merged
chains: copies shifted by multiples of 500), the symmetry operations of Y are taken from Y itself (the
same as rebuild-symmetry-refs.py: superposition of its copies on their rigid core, or on the binding site
when the copies differ by domain motions). Each operation is applied to the other partner X; a placement
is kept as a new reference when it is a distinct binding site (centroid more than 5 A from the original
site and from the sites already kept), X does not clash with Y (no CA pair closer than 3 A) and the
interface keeps a comparable size (Y residues within 8 A of X: 50 to 200 % of the original).
Copy the new target.pdb files into place, then run make-ana-files.sh on them.
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
src = (HERE / "rebuild-symmetry-refs.py").read_text()
argv = sys.argv
sys.argv = [argv[0], "/dev/null"]
exec(src[:src.index("\nfor c in sys.argv[2:]")])  # D, kabsch, cas, segments, operations, local_operations
sys.argv = argv
OUT = Path(sys.argv[1])


def interface(Xc, Yc, cutoff=8.0):
    return int((np.sqrt(((Yc[:, None] - Xc[None]) ** 2).sum(2)).min(1) < cutoff).sum())


for c in sys.argv[2:]:
    target = (D / c / f"{c}_target.pdb").read_text().splitlines(keepends=True)
    chains = {"A": cas(target, "A"), "B": cas(target, "B")}
    k = 1
    for xch, ych in (("A", "B"), ("B", "A")):
        X, Y = chains[xch], chains[ych]
        ops = operations(Y) or local_operations(Y, X)
        if not ops:
            continue
        Xc = np.array([v[1] for v in X.values()]); Yc = np.array([v[1] for v in Y.values()])
        ref_if = interface(Xc, Yc)
        sites = [Xc.mean(0)]
        for name, R, t, rmsd, n in sorted(ops, key=lambda o: o[3]):
            Xm = Xc @ R.T + t
            cen = Xm.mean(0)
            if min(np.linalg.norm(cen - s) for s in sites) < 5.0:
                continue
            clash = int((np.sqrt(((Yc[:, None] - Xm[None]) ** 2).sum(2)) < 3.0).sum())
            iface = interface(Xm, Yc)
            if clash or not (0.5 * ref_if <= iface <= 2.0 * ref_if):
                print(f"{c}: operation {name} of chain {ych} rejected (clashes {clash}, interface {iface}/{ref_if})")
                continue
            sites.append(cen)
            k += 1
            out = []
            for l in target:
                if l.startswith(("ATOM", "HETATM")) and l[21] == xch:
                    xyz = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]) @ R.T + t
                    l = l[:30] + f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}" + l[54:]
                out.append(l)
            dst = OUT / c / f"ana_scripts_{k}" / "target.pdb"
            dst.parent.mkdir(parents=True, exist_ok=True)
            dst.write_text("".join(out))
            print(f"{c} ana_scripts_{k}: chain {xch} moved by operation {name} of chain {ych} "
                  f"(self-fit {rmsd:.2f} A on {n} CA; interface {iface}/{ref_if} residues)")
    if k == 1:
        print(f"{c}: no symmetric partner")
