"""Swap the numbering of protomers 2 and 3 of the 4FQI ligand in the bound files, so that the bound
trimer is numbered like the unbound one (unbound protomer k sits on bound protomer k).
python swap_protomers.py FILE CHAIN N   (N residues per protomer, numbered 1..3N)"""
import sys
from pathlib import Path
f, chain, n = Path(sys.argv[1]), sys.argv[2], int(sys.argv[3])
lines = f.read_text().splitlines(keepends=True)
idx = [i for i, l in enumerate(lines) if l.startswith(("ATOM", "HETATM")) and l[21] == chain]
assert idx == list(range(idx[0], idx[-1] + 1)) or all(lines[i].startswith(("ATOM", "HETATM", "TER")) for i in range(idx[0], idx[-1] + 1))
block = {1: [], 2: [], 3: []}
for i in idx:
    k = int(lines[i][22:26]); p = (k - 1) // n + 1
    block[p].append(lines[i])
names = {p: [(int(l[22:26]) - (p - 1) * n, l[17:20]) for l in block[p] if l[12:16] == " CA "] for p in block}
assert names[2] == names[3], "protomers 2 and 3 differ in residues"
renum = lambda l, d: l[:22] + f"{int(l[22:26]) + d:4d}" + l[26:]
new_chain = block[1] + [renum(l, -n) for l in block[3]] + [renum(l, n) for l in block[2]]
inner = [lines[i] for i in range(idx[0], idx[-1] + 1) if i not in set(idx)]
out = lines[:idx[0]] + new_chain + inner + lines[idx[-1] + 1:]
# sequential serials
s = next(int(l[6:11]) for l in out if l.startswith(("ATOM", "HETATM")))
for j, l in enumerate(out):
    if l.startswith(("ATOM", "HETATM")) or (l.startswith("TER") and l[6:11].strip().isdigit()):
        out[j] = l[:6] + f"{s:5d}" + l[11:]; s += 1
f.write_text("".join(out))
print(f.name, "swapped protomers 2 and 3:", len(block[2]), len(block[3]), "atoms")
