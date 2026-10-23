"""MET selenium leftovers: atom SE -> SD, and element S for MET SD (columns 77-80)."""
import sys
from collections import Counter
total = Counter()
for f in sys.argv[1:]:
    lines = open(f).read().splitlines(keepends=True)
    changed = Counter()
    for i, l in enumerate(lines):
        if not l.startswith(("ATOM", "HETATM")) or l[17:20] != "MET":
            continue
        name = l[12:16].strip()
        if name not in ("SE", "SD"):
            continue
        body = l.rstrip("\n")
        new = body
        if name == "SE":
            new = new[:12] + " SD " + new[16:]
            changed["SE renamed SD"] += 1
        if len(new) > 76 and new[76:80].strip() != "S":
            changed[f"element {new[76:80].strip()!r} -> S"] += 1
            new = new[:76].ljust(76) + " S"
        if new != body:
            lines[i] = new + "\n"
    if changed:
        open(f, "w").write("".join(lines))
        print(f, dict(changed))
        total.update(changed)
print("total", dict(total))
