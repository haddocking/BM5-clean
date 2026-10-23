"""Rebuild ana_scripts/target-unbound.pdb = XXXX_r_u.pdb + XXXX_l_u.pdb (chain IDs A/B kept,
segid = chain), with the TER records of the inputs, a TER after each partner, sequential serials
and END.

python make-target-unbound.py           report the complexes whose target-unbound.pdb differs
python make-target-unbound.py --write   rewrite those files
"""
import os
import sys
from pathlib import Path
D = Path(os.environ.get("BM5_HADDOCK_READY", Path(__file__).resolve().parent.parent / "HADDOCK-ready"))


def build(c):
    out = []
    for s in ("r", "l"):
        part = [l for l in (D / c / f"{c}_{s}_u.pdb").read_text().splitlines()
                if l.startswith(("ATOM", "HETATM", "TER"))]
        while part and part[-1].startswith("TER"):
            part.pop()
        out += part + ["TER"]
    n, res = 1, []
    for l in out:
        if l.startswith(("ATOM", "HETATM")):
            l = l.ljust(76)
            l = l[:6] + f"{n:5d}" + l[11:72] + f"{l[21]:<4}" + l[76:]
            last = l
        else:  # TER: serial and residue of the preceding atom
            l = "TER   " + f"{n:5d}" + "      " + last[17:27]
        res.append(l.rstrip() if l.startswith("TER") else l)
        n += 1
    return "\n".join(res) + "\nEND\n"


def atoms(text):
    return [l[12:66] for l in text.splitlines() if l.startswith(("ATOM", "HETATM"))]


write = "--write" in sys.argv
for cdir in sorted(p for p in D.iterdir() if (p / "ana_scripts" / "target-unbound.pdb").exists()):
    c = cdir.name
    new = build(c)
    tu = cdir / "ana_scripts" / "target-unbound.pdb"
    old = tu.read_text()
    if atoms(new) != atoms(old) or sum(l.startswith("TER") for l in new.splitlines()) != sum(l.startswith("TER") for l in old.splitlines()):
        a, b = atoms(new), atoms(old)
        print(c, f"atoms new {len(a)} old {len(b)}, content differs: {a != b}")
        if write:
            tu.write_text(new)
