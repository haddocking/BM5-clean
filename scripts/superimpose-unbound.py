"""Superimpose each unbound partner (XXXX_r_u.pdb, XXXX_l_u.pdb) onto its matched bound partner
(XXXX_r_b-matched.pdb, XXXX_l_b-matched.pdb) over its interface: backbone atoms (N, CA, C, O) of the
residues in ana_scripts/target.izoneA (receptor) or target.izoneB (ligand). All atoms of the unbound
file (including HETATM) are moved. Partners whose interface zone has fewer than 4 residues present in
both files are fitted on all their backbone atoms instead.

Afterwards, rebuild target-unbound.pdb (make-target-unbound.py --write) and the i-RMSD tables
(make-irmsd-tables.sh).

python superimpose-unbound.py [COMPLEX...]     (all complexes by default)
"""
import os
import re
import sys
from pathlib import Path

import numpy as np

D = Path(os.environ.get("BM5_HADDOCK_READY", Path(__file__).resolve().parent.parent / "HADDOCK-ready"))
BB = ("N", "CA", "C", "O")


def backbone(path):
    out = {}
    for l in path.read_text().splitlines():
        if l.startswith("ATOM") and l[12:16].strip() in BB:
            out[(int(l[22:26]), l[12:16].strip())] = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return out


def zone_residues(path):
    res = set()
    for l in path.read_text().splitlines():
        m = re.findall(r"[AB]\.?(-?\d+)", l)
        if len(m) == 2:
            res.update(range(int(m[0]), int(m[1]) + 1))
    return res


def kabsch(X, Y):
    cx, cy = X.mean(0), Y.mean(0)
    U, S, Vt = np.linalg.svd((X - cx).T @ (Y - cy))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, cy - cx @ R.T


def superimpose(c, p):
    u_path, b_path = D / c / f"{c}_{p}_u.pdb", D / c / f"{c}_{p}_b-matched.pdb"
    u, b = backbone(u_path), backbone(b_path)
    zone = zone_residues(D / c / "ana_scripts" / f"target.izone{'A' if p == 'r' else 'B'}")
    keys = [k for k in u if k in b and k[0] in zone]
    used = "interface"
    if len({k[0] for k in keys}) < 4:
        keys, used = [k for k in u if k in b], "all backbone"
    X = np.array([u[k] for k in keys]); Y = np.array([b[k] for k in keys])
    R, t = kabsch(X, Y)
    rmsd = float(np.sqrt(((X @ R.T + t - Y) ** 2).sum(1).mean()))
    moved = float(np.sqrt(((X @ R.T + t - X) ** 2).sum(1).mean()))
    out = []
    for l in u_path.read_text().splitlines():
        if l.startswith(("ATOM", "HETATM")):
            x = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]) @ R.T + t
            l = l[:30] + f"{x[0]:8.3f}{x[1]:8.3f}{x[2]:8.3f}" + l[54:]
        out.append(l)
    u_path.write_text("\n".join(out) + "\n")
    return used, len({k[0] for k in keys}), rmsd, moved


def main():
    complexes = sys.argv[1:] or sorted(p.name for p in D.iterdir() if (p / "ana_scripts").is_dir())
    for c in complexes:
        for p in ("r", "l"):
            used, n, rmsd, moved = superimpose(c, p)
            print(f"{c}_{p}\t{used}\t{n} residues\tfit rmsd {rmsd:.3f}\tinterface moved {moved:.3f}")


if __name__ == "__main__":
    main()
