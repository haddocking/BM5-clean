#!/bin/bash
# Download the unbound (cif/) and bound (cifb/) PDB entries of the BM5 table (unbound_ids.json,
# bound_ids.json in this directory) as mmCIF from RCSB, into the current working directory.
H=$(cd "$(dirname "$0")" && pwd)
mkdir -p cif cifb
python3 - "$H" <<'PY' | while read dir id; do
import json, sys
h = sys.argv[1]
u = json.load(open(f"{h}/unbound_ids.json")); b = json.load(open(f"{h}/bound_ids.json"))
for v in u.values():
    for spec in v.values():
        print("cif", spec[:4])
for spec in b.values():
    print("cifb", spec[:4])
PY
  [ -s $dir/$id.cif ] || curl -sf https://files.rcsb.org/download/$id.cif -o $dir/$id.cif || echo "failed $id"
done
