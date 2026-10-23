"""Origin (PDB entry chain and symmetry copy) of every residue of the BM5-clean
unbound (_u) and bound (_b-matched) files, by superposition on the PDB entries
(same method as gaps.py). Used to tell TER records at gaps (same chain on both
sides) from TER records between merged chains.

Writes origins.json: {"XXXX_r_u": {"entry": ..., "res": {"A  12 ": "B/0", ...}}, ...}
"""
import os
import json
import re
import sys
import warnings
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser

sys.argv = sys.argv[:1]
warnings.filterwarnings("ignore")
src = Path(__file__).with_name("gaps.py").read_text()
# reuse kabsch, spec, file_residues, fit (everything before analyse) from gaps.py
exec(src[:src.index("def analyse(")])

parser = MMCIFParser(QUIET=True)
D = Path(os.environ.get("BM5_HADDOCK_READY", Path(__file__).resolve().parents[2] / "HADDOCK-ready"))
HERE = Path(__file__).resolve().parent  # BM5 table IDs (unbound_ids.json, bound_ids.json)
uids = json.load(open(HERE / "unbound_ids.json"))
bids = json.load(open(HERE / "bound_ids.json"))
_st = {}


def structure(pid, cifdir):
    if pid not in _st:
        _st[pid] = list(parser.get_structure(pid, Path(cifdir) / f"{pid}.cif"))
    return _st[pid]


def origins(pdb, pid, chains, cifdir, model=0):
    """{file residue key (chain + number + ins): "entrychain/copy"} and number of CA residues mapped."""
    models = structure(pid, cifdir)
    mdl = models[min(model, len(models) - 1)]
    chains = list(chains) if chains else [c.id for c in mdl]
    pdb_ca = [(ch, f"{r.id[1]}{r.id[2].strip()}", r.get_resname(), r["CA"].coord)
              for ch in chains if ch in mdl for r in mdl[ch] if "CA" in r]
    if not pdb_ca:
        return {}, 0, 0
    coords = np.array([x[3] for x in pdb_ca])
    fres = []
    for line in Path(pdb).read_text().splitlines():
        if line.startswith("ENDMDL"):
            break
        if line.startswith(("ATOM", "HETATM")):
            key = line[21:27]
            if not fres or fres[-1]["key"] != key:
                fres.append({"key": key, "name": line[17:20], "atoms": {}, "het": line.startswith("HETATM")})
            fres[-1]["atoms"][line[12:16].strip()] = np.array(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])])
    fca_all = [(i, r["atoms"]["CA"]) for i, r in enumerate(fres) if "CA" in r["atoms"]]
    out = {}
    for copy in range(8):
        fca = [(i, c) for i, c in fca_all if fres[i]["key"] not in out]
        if len(fca) < 10:
            break
        k, ok = fit(fca, fres, pdb_ca, coords, chains)
        if ok is None or ok.sum() < 10:
            break
        for j, (i, _) in enumerate(fca):
            if ok[j]:
                out[fres[i]["key"]] = f"{pdb_ca[k[j]][0]}/{copy}"
    return out, len(out), len(fca_all)


def best_origin(pdb, candidates, cifdir):
    """Try each (pid, chains, model) candidate; for NMR entries also the other models."""
    best = ({}, -1, 0, None)
    for pid, chains, model in candidates:
        try:
            o, n, tot = origins(pdb, pid, chains, cifdir, (model or 1) - 1)
        except Exception:
            continue
        if n > best[1]:
            best = (o, n, tot, f"{pid}_{chains}")
    if best[1] < best[2]:
        for pid, chains, _ in candidates:
            for m in range(min(len(structure(pid, cifdir)), 60)):
                o, n, tot = origins(pdb, pid, chains, cifdir, m)
                if n > best[1]:
                    best = (o, n, tot, f"{pid}_{chains}({m + 1})")
    return best


result = {}
for name in sorted(p.name for p in D.iterdir() if p.is_dir() and (p / "ana_scripts").is_dir()):
    for side in ("r", "l"):
        u = D / name / f"{name}_{side}_u.pdb"
        if u.exists() and name in uids:
            cands = [spec(uids[name][s]) for s in ("r", "l")]
            o, n, tot, entry = best_origin(u, cands, "cif")
            result[f"{name}_{side}_u"] = {"entry": entry, "mapped": n, "ca": tot, "res": o}
        b = D / name / f"{name}_{side}_b-matched.pdb"
        if b.exists() and name in bids:
            pid, rest = bids[name].split("_", 1)
            cands = [(pid, c, None) for c in rest.split(":")]
            o, n, tot, entry = best_origin(b, cands, "cifb")
            result[f"{name}_{side}_b"] = {"entry": entry, "mapped": n, "ca": tot, "res": o}
    print(name, {k: f"{v['mapped']}/{v['ca']}" for k, v in result.items() if k.startswith(name)}, flush=True)
json.dump(result, open("origins.json", "w"))
