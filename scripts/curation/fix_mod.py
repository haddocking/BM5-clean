"""Add back the modified residues dropped from BM5-clean (report section B), as their parent residue.

Unbound file: the modified residue of the unbound entry, reduced to the parent's atoms.
b-matched and _target.pdb: the residue of the bound entry (also reduced when modified).
Placement, numbering, TER handling and serials as in fix_mse.py. D-Gln35 (1S1Q, 1XD3,
2OOB) is left out on purpose (bound is Gly35).
Writes the new files to OUT/<complex>/ (BM5-clean is not touched).
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
src = (HERE / "fix_mse.py").read_text()
exec(src[:src.index("def new_lines(")])
exec(src[src.index("def insert("):src.index("def main(")])

OUT = Path(sys.argv[1])
PARENT_ATOMS = {
    "HIS": ["N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2"],
    "CYS": ["N", "CA", "C", "O", "CB", "SG"],
    "SER": ["N", "CA", "C", "O", "CB", "OG"],
}
# partner, file slot, unbound (pdb, chain, number), bound (pdb, chain, number), parent
RESIDUES = [
    ("1Y64_l", 73, ("2FXU", "A", "73"), ("1Y64", "A", "73"), "HIS"),
    ("1R8S_l", 174, ("1R8M", "E", "174"), ("1R8S", "E", "174"), "CYS"),
    ("3RVW_l", 34, ("3F5V", "A", "34"), ("3RVW", "A", "34"), "CYS"),
    ("2SNI_r", 221, ("1UBN", "A", "221"), ("2SNI", "E", "221"), "SER"),
]


def parent_coords(res, parent):
    """{atom: (occ, b, xyz)} with the parent's atom names; SE7 SEG -> SER OG at 1.43 A from CB."""
    out = {}
    for name in PARENT_ATOMS[parent]:
        if name in res:
            _, occ, b, xyz = pick(res[name])
            out[name] = (occ, b, xyz)
        elif parent == "SER" and name == "OG" and "SEG" in res:
            _, occ, b, se = pick(res["SEG"])
            cb = pick(res["CB"])[3]
            v = se - cb
            out[name] = (occ, b, cb + 1.43 * v / np.linalg.norm(v))
    missing = [a for a in PARENT_ATOMS[parent] if a not in out]
    assert not missing, (res["_comp"], missing)
    return out


def lines_for(template, coords, parent, R, t, slot):
    out = []
    for name in PARENT_ATOMS[parent]:
        occ, b, xyz = coords[name]
        x, y, z = xyz @ R.T + t
        l = template.rstrip("\n").ljust(80)
        l = ("ATOM  " + NEW + l[11:12] + fmt_name(name) + " " + parent + l[20:22] + f"{slot:4d}" + " " + l[27:30]
             + f"{x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{b:6.2f}" + l[66:76] + f"{name[0]:>2}" + "  ")
        out.append(l.rstrip() + "\n")
    return out


def main():
    files, log = {}, []

    def get(path):
        if path not in files:
            files[path] = path.read_text().splitlines(keepends=True)
        return files[path]

    for partner, slot, (upid, uch, unum), (bpid, bch, bnum), parent in RESIDUES:
        name, side = partner.rsplit("_", 1)
        ures = atoms(G / "cif" / f"{upid}.cif")[(uch, unum)]
        bres = atoms(G / "cifb" / f"{bpid}.cif")[(bch, bnum)]
        entry = {"partner": partner, "slot": slot, "parent": parent,
                 "unbound": f"{upid} {uch}:{ures['_comp']}{unum}", "bound": f"{bpid} {bch}:{bres['_comp']}{bnum}"}
        tch = "A" if side == "r" else "B"
        for kind, path, ent, ch, num, res in [
                ("u", D / name / f"{name}_{side}_u.pdb", upid, uch, unum, ures),
                ("b", D / name / f"{name}_{side}_b-matched.pdb", bpid, bch, bnum, bres),
                ("t", D / name / f"{name}_target.pdb", bpid, bch, bnum, bres)]:
            lines = get(path)
            sub = [l if (kind != "t" or l[21:22] == tch) else "REMARK\n" for l in lines]
            entd = atoms(G / ("cif" if kind == "u" else "cifb") / f"{ent}.cif")
            R, t, rmsd, n = transform(sub, entd, ch, int(num) - slot, slot)
            fres = file_residues(sub)
            assert slot not in fres, f"{path.name}: {slot} already present"
            assert rmsd < 0.05, f"{path.name}: fit rmsd {rmsd:.3f} on {n} atoms"
            template = lines[max(i for i, _ in fres[max(k for k in fres if k < slot)].values())]
            new = lines_for(template, parent_coords(res, parent), parent, R, t, slot)
            lines[:] = insert(lines, slot, new, tch if kind == "t" else None)
            entry[kind] = {"rmsd": round(rmsd, 3), "fit_atoms": n, "atoms": len(new)}
        log.append(entry)

    for path, lines in files.items():
        lines, _ = renumber_serials(lines)
        dst = OUT / path.parent.name / path.name
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.write_text("".join(lines))
    json.dump(log, open(OUT / "added.json", "w"), indent=1)
    for e in log:
        print(e["partner"], e["slot"], e["parent"], e["unbound"], e["bound"], {k: (e[k]["rmsd"], e[k]["atoms"]) for k in "ubt"})


if __name__ == "__main__":
    main()
