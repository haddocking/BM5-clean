"""Write the haddock3 benchmark input lists (haddock-runner format) in haddock3-scripts/:

- input.list : all complexes of HADDOCK-ready except those of categories/AB.txt (antibody-antigen
  complexes without unbound antibody);
- input-<NAME>.list : the complexes of each HADDOCK-ready/categories/<NAME>.txt list.

Each list contains one block per complex (unbound receptor and ligand, reference complex, ambig5
restraints), then the ligand topology and parameter files, the hbonds restraint files and the ambig
restraint files of these complexes. Paths are relative to haddock3-scripts/.

python make-haddock3-lists.py
"""
import os
from pathlib import Path

ROOT = Path(os.environ.get("BM5_ROOT", Path(__file__).resolve().parent.parent))
D = ROOT / "HADDOCK-ready"
OUT = ROOT / "haddock3-scripts"


def write(path, complexes):
    rel = lambda c, f: f"../HADDOCK-ready/{c}/{f}"
    lines = []
    for c in complexes:
        lines += [f"# {c}"] + [rel(c, f"{c}_{s}") for s in ("r_u.pdb", "l_u.pdb", "target.pdb", "ambig5.tbl")] + [""]
    sections = (("# ligand top list", "ligand.top"), ("# ligand param list", "ligand.param"))
    for header, suffix in sections:
        lines.append(header)
        lines += [rel(c, f"{c}_{suffix}") for c in complexes if (D / c / f"{c}_{suffix}").exists()]
        lines.append("")
    lines.append("# hbonds files containing the restraints to keep bodies together")
    lines += [rel(c, f"{c}_hbonds.tbl") for c in complexes if (D / c / f"{c}_hbonds.tbl").exists()]
    lines += [rel(c, f"{c}_ambig.tbl") for c in complexes]
    path.write_text("\n".join(lines) + "\n")


def main():
    every = sorted(p.name for p in D.iterdir() if (p / "ana_scripts").is_dir())
    ab = set((D / "categories" / "AB.txt").read_text().split())
    write(OUT / "input.list", [c for c in every if c not in ab])
    print(f"input.list: {len(every) - len(ab)} complexes (all but AB)")
    for f in sorted((D / "categories").glob("*.txt")):
        names = [l.strip() for l in f.read_text().splitlines() if l.strip()]
        missing = [c for c in names if c not in every]
        assert not missing, f"{f.name}: unknown complexes {missing}"
        write(OUT / f"input-{f.stem}.list", sorted(names))
        print(f"input-{f.stem}.list: {len(names)} complexes")


if __name__ == "__main__":
    main()
