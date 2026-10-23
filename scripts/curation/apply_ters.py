"""Remove TER records at gaps inside a chain of the BM5-clean files, keep those between chains.

A TER between two ATOM residues with the same chain ID is removed when both residues
come from the same chain (and symmetry copy) of the PDB entry (origins.json). It is kept
when they come from different chains, or when either residue could not be mapped
(listed as undecided). Atom serials are renumbered when a TER is removed.

python apply_ters.py           dry run: writes ter_decisions.tsv
python apply_ters.py --write   edits the files
"""
import os
import json
import sys
from pathlib import Path

D = Path(os.environ.get("BM5_HADDOCK_READY", Path(__file__).resolve().parents[2] / "HADDOCK-ready"))
HERE = Path(__file__).resolve().parent  # BM5 table IDs (unbound_ids.json, bound_ids.json)
write = "--write" in sys.argv
orig = json.load(open("origins.json"))


def labels(name, kind):
    """Residue key (chain + number + ins) -> origin, for one of the four file kinds."""
    if kind in ("r_u", "l_u", "r_b-matched", "l_b-matched"):
        side, t = kind.split("_", 1)
        return orig.get(f"{name}_{side}_{'u' if t == 'u' else 'b'}", {}).get("res", {})
    # target (bound) / target-unbound: chain A = receptor, B = ligand, same numbering
    t = "b" if kind == "target" else "u"
    out = {}
    for side, ch in (("r", "A"), ("l", "B")):
        for k, v in orig.get(f"{name}_{side}_{t}", {}).get("res", {}).items():
            out[ch + k[1:]] = f"{side}:{v}"
    return out


rows = []
nfiles = 0
for cdir in sorted(p for p in D.iterdir() if (p / "ana_scripts").is_dir()):
    name = cdir.name
    files = [(cdir / f"{name}_{k}.pdb", k) for k in ("r_u", "l_u", "r_b-matched", "l_b-matched", "target")]
    files.append((cdir / "ana_scripts" / "target-unbound.pdb", "target-unbound"))
    for path, kind in files:
        if not path.exists():
            continue
        lab = labels(name, kind)
        lines = path.read_text().splitlines(keepends=True)
        out, removed = [], 0
        for i, l in enumerate(lines):
            if l.startswith("TER"):
                prev = next((x for x in reversed(lines[:i]) if x.startswith(("ATOM", "HETATM", "ENDMDL", "MODEL"))), "")
                nxt = next((x for x in lines[i + 1:] if x.startswith(("ATOM", "HETATM", "ENDMDL", "MODEL", "END"))), "")
                if prev.startswith("ATOM") and nxt.startswith("ATOM") and prev[21] == nxt[21]:
                    a, b = lab.get(prev[21:27]), lab.get(nxt[21:27])
                    if a and b and a == b:
                        decision = "remove"
                    elif a and b:
                        decision = "keep: different chains"
                    else:
                        decision = "keep: undecided (unmapped residue)"
                    rows.append([name, path.name, prev[17:20] + prev[21:27].strip(), nxt[17:20] + nxt[21:27].strip(),
                                 str(a), str(b), decision])
                    if decision == "remove":
                        removed += 1
                        continue
            out.append(l)
        if removed and write:
            first = next(l for l in out if l.startswith(("ATOM", "HETATM")))
            n = int(first[6:11])
            for j, l in enumerate(out):
                if l.startswith("MODEL"):
                    n = int(first[6:11])  # each model numbered from the start
                if l.startswith(("ATOM", "HETATM")) or (l.startswith("TER") and l[6:11].strip().isdigit()):
                    out[j] = l[:6] + f"{n:5d}" + l[11:]
                    n += 1
            path.write_text("".join(out))
            nfiles += 1

Path("ter_decisions.tsv").write_text("complex\tfile\tbefore\tafter\torigin_before\torigin_after\tdecision\n"
                                     + "".join("\t".join(r) + "\n" for r in rows))
from collections import Counter
print(Counter(r[6] for r in rows))
print("files edited:", nfiles if write else "(dry run)")
