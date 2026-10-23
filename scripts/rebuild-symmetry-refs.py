"""Rebuild out-of-date symmetry references (ana_scripts_K/target.pdb) from the current XXXX_target.pdb.

A symmetry reference is the complex with one partner (X) on another binding site of the
other, symmetric partner (Y): X moved by a symmetry operation S of Y, Y unchanged.
The operations of Y are taken from Y itself: residue-number offsets d at which Y repeats
(same residue name at n and n + d) give the superposition Y[n] -> Y[n + d]. For each old
reference, the operation (and the choice of X/Y) that best reproduces the old placement of X
relative to Y is kept: the old references keep Y in place, so the score is the median
distance between the moved X and X in the old reference (same residue numbers).
Copies are superimposed on their rigid core; if they differ too much (domain motions), on the
binding site only (residues of Y within 12 A of X).

python rebuild-symmetry-refs.py OUT COMPLEX...   writes OUT/<complex>/ana_scripts_K/target.pdb
(then copy it into place and run make-ana-files.sh on it)
"""
import os
import sys
from pathlib import Path

import numpy as np

D = Path(os.environ.get("BM5_HADDOCK_READY", Path(__file__).resolve().parent.parent / "HADDOCK-ready"))
OUT = Path(sys.argv[1])


def kabsch(X, Y):
    cx, cy = X.mean(0), Y.mean(0)
    U, S, Vt = np.linalg.svd((X - cx).T @ (Y - cy))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, cy - cx @ R.T


def cas(lines, chain):
    return {int(l[22:26]): (l[17:20], np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]))
            for l in lines if l.startswith("ATOM") and l[12:16] == " CA " and l[21] == chain}


def segments(nums, gap=30):
    """Contiguous numbering segments (split at jumps > gap)."""
    segs, cur = [], [nums[0]]
    for a, b in zip(nums, nums[1:]):
        if b - a > gap:
            segs.append(cur); cur = []
        cur.append(b)
    segs.append(cur)
    return segs


def operations(Y):
    """Symmetry operations of Y: for each numbering offset d and each segment of Y, the
    superposition of the segment's residues n onto residues n + d (same residue name).
    Per segment, so that non-cyclic arrangements (e.g. D2 tetramers) give separate operations."""
    nums = sorted(Y)
    segs = segments(nums)
    ops = []
    for seg in segs:
        lo, hi = seg[0], seg[-1]
        best = {}
        for d in range(nums[0] - hi, nums[-1] - lo + 1):
            if d == 0:
                continue
            pairs = [n for n in seg if n + d in Y and Y[n + d][0] == Y[n][0]]
            if len(pairs) < max(20, 0.6 * len(seg)):
                continue
            key = round(d / 25)
            if key not in best or len(pairs) > best[key][1]:
                best[key] = (d, len(pairs), pairs)
        for d, n, pairs in best.values():
            X1 = np.array([Y[m][1] for m in pairs]); X2 = np.array([Y[m + d][1] for m in pairs])
            # rigid core: iterative outlier rejection (copies may differ by domain motions)
            keep = np.ones(len(pairs), bool)
            for _ in range(10):
                R, t = kabsch(X1[keep], X2[keep])
                dist = np.linalg.norm(X1 @ R.T + t - X2, axis=1)
                new = dist < max(1.0, 2 * np.median(dist))
                if (new == keep).all() or new.sum() < 0.5 * len(pairs):
                    break
                keep = new
            rmsd = float(np.sqrt((dist[keep] ** 2).mean()))
            if rmsd < 1.5:
                ops.append((f"{lo}-{hi}{d:+d}", R, t, rmsd, int(keep.sum())))
    return ops


def local_operations(Y, X, cutoff=12.0):
    """Fallback when the copies of Y differ too much for a global fit: superposition of the
    residues of Y within cutoff of X (the binding site) onto their copies at offset d."""
    Xc = np.array([v[1] for v in X.values()])
    site = [n for n, v in Y.items() if np.linalg.norm(Xc - v[1], axis=1).min() < cutoff]
    nums = sorted(Y)
    ops = []
    for d in range(nums[0] - max(site), nums[-1] - min(site) + 1):
        pairs = [n for n in site if d and n + d in Y and Y[n + d][0] == Y[n][0]]
        if len(pairs) < max(10, 0.8 * len(site)):
            continue
        X1 = np.array([Y[m][1] for m in pairs]); X2 = np.array([Y[m + d][1] for m in pairs])
        R, t = kabsch(X1, X2)
        rmsd = float(np.sqrt(((X1 @ R.T + t - X2) ** 2).sum(1).mean()))
        if rmsd < 1.5:
            ops.append((f"site{d:+d}", R, t, rmsd, len(pairs)))
    return ops


def score(candX, old, oldXchain, X):
    """Median distance between the moved X and X in the old reference (same frame: Y is not moved)."""
    oldX = cas(old, oldXchain)
    common = [n for n in X if n in oldX and oldX[n][0] == X[n][0]]
    return float(np.median([np.linalg.norm(candX[n] - oldX[n][1]) for n in common]))


for c in sys.argv[2:]:
    target = (D / c / f"{c}_target.pdb").read_text().splitlines(keepends=True)
    chains = {"A": cas(target, "A"), "B": cas(target, "B")}
    for k in range(2, 7):
        oldf = D / c / f"ana_scripts_{k}" / "target.pdb"
        if not oldf.exists():
            continue
        old = oldf.read_text().splitlines()
        best = None
        for xch, ych in (("A", "B"), ("B", "A")):
            X, Y = chains[xch], chains[ych]
            # old chain holding X: the one sharing most numbered residues with X
            oldXchain = max("AB", key=lambda ch: sum(1 for n, v in cas(old, ch).items() if n in X and X[n][0] == v[0]))
            for name, R, t, rmsd, n in (operations(Y) or local_operations(Y, X)):
                candX = {m: v[1] @ R.T + t for m, v in X.items()}
                s = score(candX, old, oldXchain, X)
                if best is None or s < best[0]:
                    best = (s, xch, name, R, t, rmsd, n)
        if best is None:
            print(f"{c} ana_scripts_{k}: no symmetry operation found"); continue
        s, xch, name, R, t, rmsd, n = best
        ych = "B" if xch == "A" else "A"
        Ycur = chains[ych]
        oldYchain = "B" if max("AB", key=lambda ch: sum(1 for m, v in cas(old, ch).items() if m in chains[xch] and chains[xch][m][0] == v[0])) == "A" else "A"
        oY = cas(old, oldYchain)
        P = np.array([v[1] for v in Ycur.values()]); Q = np.array([v[1] for v in oY.values()])
        yres = float(np.median(np.sqrt(((P[:, None] - Q[None]) ** 2).sum(2)).min(1)))
        out = []
        for l in target:
            if l.startswith(("ATOM", "HETATM")) and l[21] == xch:
                xyz = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]) @ R.T + t
                l = l[:30] + f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}" + l[54:]
            out.append(l)
        dst = OUT / c / f"ana_scripts_{k}" / "target.pdb"
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.write_text("".join(out))
        print(f"{c} ana_scripts_{k}: moved chain {xch} by operation {name} of chain {'B' if xch == 'A' else 'A'} "
              f"(self-fit {rmsd:.2f} A on {n} CA); moved partner vs old: median {s:.2f} A; fixed partner vs old frame {yres:.2f} A")
