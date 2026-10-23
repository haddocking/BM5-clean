"""Add back the selenomethionines dropped from BM5-clean (report section A2), as MET.

For each residue in lost_residues.tsv that is a HETATM MSE in the unbound entry:
- unbound file: the MSE atoms of the unbound entry, SE renamed SD (element S),
  superimposed into the file frame by a fit on the neighbouring residues;
- b-matched file and _target.pdb: the MET of the bound entry, superimposed the same way.
The residue takes the free number between its two neighbours in the file.
Writes the new files to OUT/<complex>/ (BM5-clean is not touched).
"""
import os
import csv
import json
import re
import shutil
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

G = Path.cwd()  # working directory: cif/, cifb/, gaps.json, lost_residues.tsv
D = Path(os.environ.get("BM5_HADDOCK_READY", Path(__file__).resolve().parents[2] / "HADDOCK-ready"))
HERE = Path(__file__).resolve().parent  # BM5 table IDs (unbound_ids.json, bound_ids.json)
OUT = Path(sys.argv[1])
WINDOW = 4  # residues on each side used for the fit
MET_ATOMS = ["N", "CA", "C", "O", "CB", "CG", "SD", "CE"]
RENAME = {"SE": "SD"}
NEW = "*****"  # serial of inserted atoms until assign_serials
ELEMENT = {"N": "N", "CA": "C", "C": "C", "O": "O", "CB": "C", "CG": "C", "SD": "S", "CE": "C"}


def kabsch(X, Y):
    cx, cy = X.mean(0), Y.mean(0)
    U, S, Vt = np.linalg.svd((X - cx).T @ (Y - cy))
    dd = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, dd]) @ U.T
    return R, cy - cx @ R.T


_cif = {}


def atoms(cif):
    """{(auth chain, auth number + ins): {atom name: [(altloc, occ, b, xyz)]}} for model 1."""
    if cif in _cif:
        return _cif[cif]
    d = MMCIF2Dict(str(cif))
    out = defaultdict(lambda: defaultdict(list))
    first = d["_atom_site.pdbx_PDB_model_num"][0]
    cols = ["pdbx_PDB_model_num", "auth_asym_id", "auth_seq_id", "pdbx_PDB_ins_code", "auth_comp_id",
            "auth_atom_id", "label_alt_id", "occupancy", "B_iso_or_equiv", "Cartn_x", "Cartn_y", "Cartn_z"]
    for m, ch, num, ins, comp, name, alt, occ, b, x, y, z in zip(*[d[f"_atom_site.{c}"] for c in cols]):
        if m != first:
            continue
        key = (ch, num + ("" if ins in ".?" else ins))
        out[key]["_comp"] = comp
        out[key][name].append((alt if alt not in ".?" else "", float(occ), float(b),
                               np.array([float(x), float(y), float(z)])))
    _cif[cif] = out
    return out


def pick(confs):
    """Altloc A (or the only conformation)."""
    for c in confs:
        if c[0] in ("", "A"):
            return c
    return confs[0]


def file_residues(lines):
    """{residue number: {atom name: (line index, xyz)}} for ATOM/HETATM lines."""
    res = defaultdict(dict)
    for i, l in enumerate(lines):
        if l.startswith(("ATOM", "HETATM")):
            res[int(l[22:26])][l[12:16].strip()] = (i, np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]))
    return res


def transform(lines, entry, chain, offset, slot):
    """Fit entry -> file on backbone atoms of the residues around slot (file number n <-> entry n+offset)."""
    fres = file_residues(lines)
    X, Y = [], []
    for n in range(slot - WINDOW, slot + WINDOW + 1):
        if n == slot or n not in fres:
            continue
        e = entry.get((chain, str(n + offset)))
        if not e:
            continue
        for a in ("N", "CA", "C", "O"):
            if a in fres[n] and a in e:
                X.append(pick(e[a])[3]); Y.append(fres[n][a][1])
    X, Y = np.array(X), np.array(Y)
    R, t = kabsch(X, Y)
    rmsd = float(np.sqrt(((X @ R.T + t - Y) ** 2).sum(1).mean()))
    return R, t, rmsd, len(X)


def fmt_name(name):
    return f" {name:<3}" if len(name) < 4 else name


def new_lines(template, residue, R, t, slot):
    """MET lines in the layout of template (a line of the neighbouring residue)."""
    out = []
    for name in MET_ATOMS:
        src = next((residue[s] for s, d in [(name, name)] + [(k, v) for k, v in RENAME.items() if v == name]
                    if s in residue), None)
        if src is None:
            continue
        alt, occ, b, xyz = pick(src)
        x, y, z = xyz @ R.T + t
        l = template.rstrip("\n").ljust(80)
        l = ("ATOM  " + NEW + l[11:12] + fmt_name(name) + " " + "MET" + l[20:22] + f"{slot:4d}" + " " + l[27:30]
             + f"{x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{b:6.2f}" + l[66:76] + f"{ELEMENT[name]:>2}" + "  ")
        out.append(l.rstrip() + "\n" if not template.rstrip("\n")[76:].strip() and len(template.rstrip("\n")) <= 76
                   else l + "\n")
    return out


def insert(lines, slot, new, chain=None):
    """Insert new after the preceding residue of the same chain.

    The files have a TER at each chain break: when the residue closes the break the
    TER is removed; when a break remains before it (unmodelled residues) the residue
    goes after the TER."""
    keep = [i for i, l in enumerate(lines) if l.startswith(("ATOM", "HETATM"))
            and (chain is None or l[21] == chain)]
    num = {i: int(lines[i][22:26]) for i in keep}
    after = max(n for n in num.values() if n < slot)
    before = min(n for n in num.values() if n > slot)
    i_after = max(i for i in keep if num[i] == after)
    i_before = min(i for i in keep if num[i] == before)
    assert i_after < i_before and lines[i_after][21] == lines[i_before][21], "neighbours not adjacent"
    between = lines[i_after + 1:i_before]
    assert all(l.startswith("TER") for l in between), between
    if slot == after + 1:
        ter = [] if slot + 1 == before else between
        return lines[:i_after + 1] + new + ter + lines[i_before:]
    assert slot + 1 == before, "residue in the middle of a gap"
    return lines[:i_before] + new + lines[i_before:]


def renumber_serials(lines):
    """Sequential atom serials (from the first serial of the file)."""
    out, n = [], next(int(l[6:11]) for l in lines if l.startswith(("ATOM", "HETATM")) and l[6:11] != NEW)
    for l in lines:
        if l.startswith(("ATOM", "HETATM")):
            l = l[:6] + f"{n:5d}" + l[11:]; n += 1
        elif l.startswith("TER") and len(l) > 11 and l[6:11].strip().isdigit():
            l = l[:6] + f"{n:5d}" + l[11:]; n += 1
        out.append(l)
    return out, True


def assign_serials(orig, lines):
    """Serials for the inserted atoms (marked NEW).

    A file numbered sequentially is renumbered. Otherwise the new atoms take free
    serials between their neighbours; the file is renumbered only if there is no room."""
    serials = [int(l[6:11]) for l in orig if l.startswith(("ATOM", "HETATM"))]
    sequential = serials == list(range(serials[0], serials[0] + len(serials)))
    if not sequential:
        out, ok = [], True
        atom_idx = [i for i, l in enumerate(lines) if l.startswith(("ATOM", "HETATM"))]
        new_idx = [i for i in atom_idx if lines[i][6:11] == NEW]
        lines = list(lines)
        # group consecutive new atoms
        groups, cur = [], []
        for i in new_idx:
            if cur and i != cur[-1] + 1:
                groups.append(cur); cur = []
            cur.append(i)
        if cur:
            groups.append(cur)
        for g in groups:
            prev = max(j for j in atom_idx if j < g[0] and lines[j][6:11] != NEW)
            nxt = [j for j in atom_idx if j > g[-1] and lines[j][6:11] != NEW]
            lo = int(lines[prev][6:11])
            hi = int(lines[nxt[0]][6:11]) if nxt else lo + len(g) + 1
            if hi - lo - 1 < len(g):
                ok = False; break
            for k, i in enumerate(g, start=1):
                lines[i] = lines[i][:6] + f"{lo + k:5d}" + lines[i][11:]
        if ok:
            return lines, False
    return renumber_serials(lines, force=True)


def main():
    gaps = json.load(open(G / "gaps.json"))
    ids = json.load(open(HERE / "bound_ids.json"))
    rows = [r for r in csv.DictReader(open(G / "lost_residues.tsv"), delimiter="\t")
            if r["reasons"].startswith("selenomethionine")]
    log = []
    files = {}  # path -> lines (edited in memory)

    def get(path):
        if path not in files:
            files[path] = path.read_text().splitlines(keepends=True)
        return files[path]

    for r in rows:
        name, side = r["partner"].rsplit("_", 1)
        uchain, unum = re.match(r"(\w+):MSE(\S+)", r["unbound_residue"]).groups()
        bchain, bnum = re.match(r"(\w+):MET(\S+)", r["bound_residue"]).groups()
        after, before = int(r["file_after"]), int(r["file_before"])
        # file number of the residue: same offset as the preceding file residue in the entry
        g = next(x for x in gaps[r["partner"]]["gaps"] if unum in x.get("lost", []))
        uoff = int(re.match(r"-?\d+", g["pdb_from"]).group()) - after
        slot = int(unum) - uoff
        assert after < slot < before, (r, slot)
        upid = gaps[r["partner"]]["pdb"][:4]
        bpid = ids[name].split("_")[0]
        uent, bent = atoms(G / "cif" / f"{upid}.cif"), atoms(G / "cifb" / f"{bpid}.cif")
        ures, bres = uent[(uchain, unum)], bent[(bchain, bnum)]
        assert ures["_comp"] == "MSE" and bres["_comp"] == "MET", (ures["_comp"], bres["_comp"])
        boff = int(bnum) - slot
        entry = {"partner": r["partner"], "slot": slot, "unbound": f"{upid} {uchain}:MSE{unum}",
                 "bound": f"{bpid} {bchain}:MET{bnum}"}
        for kind, path, ent, ch, off, res in [
                ("u", D / name / f"{name}_{side}_u.pdb", uent, uchain, uoff, ures),
                ("b", D / name / f"{name}_{side}_b-matched.pdb", bent, bchain, boff, bres),
                ("t", D / name / f"{name}_target.pdb", bent, bchain, boff, bres)]:
            lines = get(path)
            if kind == "t":  # only the chain of this partner
                tch = "A" if side == "r" else "B"
                sub = [l if l[21:22] == tch else "REMARK\n" for l in lines]
                R, t, rmsd, n = transform(sub, ent, ch, off, slot)
                fres = file_residues(sub)
            else:
                R, t, rmsd, n = transform(lines, ent, ch, off, slot)
                fres = file_residues(lines)
            assert slot not in fres, f"{path.name}: {slot} already present"
            assert rmsd < 0.05, f"{path.name}: fit rmsd {rmsd:.3f} on {n} atoms"
            template = lines[max(i for i, _ in fres[max(k for k in fres if k < slot)].values())]
            new = new_lines(template, res, R, t, slot)
            lines[:] = insert(lines, slot, new, tch if kind == "t" else None)
            alts = sorted({c[0] for v in res.values() if isinstance(v, list) for c in v if c[0]})
            entry[kind] = {"rmsd": round(rmsd, 3), "fit_atoms": n, "atoms": len(new), "altlocs": alts}
        log.append(entry)

    for path, lines in files.items():
        lines, ren = renumber_serials(lines)
        dst = OUT / path.parent.name / path.name
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.write_text("".join(lines))
        print(f"{path.parent.name}/{path.name}: serials {'renumbered' if ren else 'kept, new atoms in free serials'}")
    json.dump(log, open(OUT / "added.json", "w"), indent=1)
    for e in log:
        print(e["partner"], e["slot"], e["unbound"], e["bound"],
              {k: (e[k]["rmsd"], e[k]["atoms"], e[k]["altlocs"]) for k in "ubt"})


if __name__ == "__main__":
    main()
