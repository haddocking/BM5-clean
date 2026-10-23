"""Internal missing residues of the BM5-clean unbound structures, against the PDB entries."""
import os
import json, re, sys, warnings
import numpy as np
from pathlib import Path
from Bio.PDB import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.Align import PairwiseAligner
from Bio.SeqUtils import seq1
aligner = PairwiseAligner(mode="global", match_score=2, mismatch_score=-1, open_gap_score=-4, extend_gap_score=-0.5)
aligner.target_end_gap_score = aligner.query_end_gap_score = 0


def kabsch(X, Y):
    """Rotation R and translation t minimising |X R^T + t - Y|."""
    cx, cy = X.mean(0), Y.mean(0)
    U, S, Vt = np.linalg.svd((X - cx).T @ (Y - cy))
    dd = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, dd]) @ U.T
    return R, cy - cx @ R.T
warnings.filterwarnings("ignore")
D = Path(os.environ.get("BM5_HADDOCK_READY", Path(__file__).resolve().parents[2] / "HADDOCK-ready"))
HERE = Path(__file__).resolve().parent  # BM5 table IDs (unbound_ids.json, bound_ids.json)
ids = json.load(open(HERE / "unbound_ids.json"))
parser = MMCIFParser(QUIET=True)


def spec(s):
    m = re.match(r"(\w{4})_?([A-Za-z0-9]*)(?:\((\d+)\))?$", s.strip())
    return m.group(1), m.group(2), int(m.group(3)) if m.group(3) else None


def file_residues(pdb):
    res = []
    for line in Path(pdb).read_text().splitlines():
        if line.startswith("ENDMDL"): break
        if line.startswith(("ATOM", "HETATM")):
            key = line[22:27]
            if not res or res[-1]["key"] != key:
                res.append({"key": key, "name": line[17:20], "atoms": {}, "het": line.startswith("HETATM")})
            res[-1]["atoms"][line[12:16].strip()] = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return res


def scheme(cif):
    d = MMCIF2Dict(str(cif))
    out = {}
    keys = ["asym_id", "seq_id", "mon_id", "pdb_seq_num", "pdb_ins_code", "pdb_strand_id", "pdb_mon_id", "auth_seq_num"]
    cols = [d[f"_pdbx_poly_seq_scheme.{k}"] for k in keys]
    for asym, seq, mon, num, ins, strand, pmon, auth in zip(*cols):
        lst = out.setdefault(strand, [])
        if lst and lst[-1]["seq"] == int(seq): continue  # microheterogeneity: keep first
        lst.append({"seq": int(seq), "mon": mon, "num": num, "ins": "" if ins in ".?" else ins,
                    "modeled": pmon != "?" and auth != "?"})
    return out


def fit(fca, fres, pdb_ca, coords, chains):
    """Superimpose file CAs on the PDB entry (seeded from each chain and from all
    chains, rejecting outliers); return the nearest PDB CA per file CA and whether
    it lies within 0.1 A, for the fit matching the most residues."""
    fseq = "".join(seq1(fres[i]["name"]) or "X" for i, _ in fca)
    F = np.array([c for _, c in fca])
    best = (-1, None, None)
    for seed in [None] + chains:
        sel = [k for k, x in enumerate(pdb_ca) if seed is None or x[0] == seed]
        if len(sel) < 3: continue
        pseq = "".join(seq1(pdb_ca[k][2]) or "X" for k in sel)
        aln = aligner.align(fseq, pseq)[0]
        pairs = [(a, sel[b]) for (a0, a1), (b0, b1) in zip(*aln.aligned) for a, b in zip(range(a0, a1), range(b0, b1))]
        if len(pairs) < 3: continue
        X = F[[a for a, _ in pairs]]; Y = coords[[b for _, b in pairs]]
        keep = np.ones(len(pairs), bool)
        for _ in range(8):
            if keep.sum() < 3: break
            R, t = kabsch(X[keep], Y[keep])
            d = np.linalg.norm(X @ R.T + t - Y, axis=1)
            keep = d < max(0.3, np.percentile(d, 50))
        R, t = kabsch(X[keep], Y[keep]) if keep.sum() >= 3 else kabsch(X, Y)
        dist = np.linalg.norm((F @ R.T + t)[:, None, :] - coords[None, :, :], axis=2)
        k = dist.argmin(1)
        ok = dist[np.arange(len(F)), k] < 0.1
        if ok.sum() > best[0]:
            best = (int(ok.sum()), k, ok)
    return best[1], best[2]


def analyse(name, side, entry_side, model_override=None):
    pid, chains, model = spec(ids[name][entry_side])
    cif = Path("cif") / f"{pid}.cif"
    st = parser.get_structure(pid, cif)
    models = list(st)
    if model_override is not None:
        mdl = models[model_override]
    else:
        mdl = models[model - 1] if model and len(models) >= model else models[0]
    sch = scheme(cif)
    chains = list(chains) if chains else [c.id for c in mdl]
    # CA coordinates of the PDB entry
    pdb_ca = []
    for ch in chains:
        if ch not in mdl: continue
        for r in mdl[ch]:
            if "CA" in r:
                pdb_ca.append((ch, f"{r.id[1]}{r.id[2].strip()}", r.get_resname(), r["CA"].coord))
    if not pdb_ca: return {"error": f"no chains {chains} in {pid}"}
    coords = np.array([x[3] for x in pdb_ca])
    fres = file_residues(D / name / f"{name}_{side}_u.pdb")
    fca_all = [(i, r["atoms"]["CA"]) for i, r in enumerate(fres) if "CA" in r["atoms"]]
    mapped = [None] * len(fres)
    copy_of = [None] * len(fres)
    for copy in range(8):  # symmetry copies are fitted one rigid body at a time
        fca = [(i, c) for i, c in fca_all if mapped[i] is None]
        if len(fca) < 10: break
        k, ok = fit(fca, fres, pdb_ca, coords, chains)
        if ok is None or ok.sum() < 10: break
        for j, (i, _) in enumerate(fca):
            if ok[j]:
                mapped[i] = pdb_ca[k[j]][:3]; copy_of[i] = copy
    nmap = sum(m is not None for m in mapped)
    nca = sum("CA" in r["atoms"] for r in fres)
    gaps = []
    for ch in chains:
        if ch not in sch: continue
        index = {f"{x['num']}{x['ins']}": i for i, x in enumerate(sch[ch]) if x["modeled"]}
        # positions in the SEQRES of the file residues from this chain, in file order
        pos = [(i, index[m[1]]) for i, m in enumerate(mapped) if m and m[0] == ch and m[1] in index]
        pos = [(i, p) for i, p in pos]
        for (i1, p1), (i2, p2) in zip(pos, pos[1:]):
            if copy_of[i1] != copy_of[i2]:
                continue
            if p2 - p1 > 1:
                between = sch[ch][p1 + 1:p2]
                lost = [x for x in between if x["modeled"]]
                gaps.append({"chain": ch, "after": fres[i1]["key"].strip(), "before": fres[i2]["key"].strip(),
                             "pdb_from": sch[ch][p1]["num"] + sch[ch][p1]["ins"],
                             "pdb_to": sch[ch][p2]["num"] + sch[ch][p2]["ins"],
                             "n": len(between), "n_unmodeled": len(between) - len(lost),
                             "n_lost": len(lost),
                             "seq": "".join(x["mon"][:3] + " " for x in between).strip(),
                             "lost": [x["num"] + x["ins"] for x in lost]})
            elif p2 < p1:
                gaps.append({"chain": ch, "after": fres[i1]["key"].strip(), "before": fres[i2]["key"].strip(),
                             "order": "reversed"})
    # residues without CA (incomplete backbone) in the file
    noca = [r["key"].strip() + r["name"] for r in fres if "CA" not in r["atoms"] and not r["het"]]
    return {"pdb": ids[name][entry_side], "swapped": entry_side != side, "model_used": model_override, "copies": len({c for c in copy_of if c is not None}), "model": model, "file_res": len(fres), "file_ca": nca, "mapped": nmap,
            "gaps": gaps, "no_ca": noca}


out = {}
for name in sorted(ids):
    if not (D / name).is_dir(): continue
    for side in ("r", "l"):
        results = []
        for entry_side in ("r", "l"):
            try:
                results.append(analyse(name, side, entry_side))
            except Exception as e:
                results.append({"error": repr(e)[:200], "mapped": -1})
        best = max(results, key=lambda v: v.get("mapped", -1))
        if best.get("mapped", -1) < best.get("file_ca", 0):
            # NMR entries: try the other models
            entry_side = "l" if best.get("swapped") != (side == "r") else "r"
            pid = spec(ids[name][entry_side])[0]
            n_models = len(list(parser.get_structure(pid, Path("cif") / f"{pid}.cif")))
            for m in range(min(n_models, 60)):
                try:
                    v = analyse(name, side, entry_side, m)
                except Exception:
                    continue
                if v["mapped"] > best["mapped"]:
                    best = v
        out[f"{name}_{side}"] = best
json.dump(out, open("gaps.json", "w"), indent=1)
bad = {k: v for k, v in out.items() if "error" in v or v["mapped"] < v["file_ca"]}
print(len(out), "partners;", len(bad), "not fully mapped by coordinates")
for k, v in bad.items(): print(" ", k, v.get("error") or f"{v['mapped']}/{v['file_ca']} {v['pdb']}")
g = [(k, x) for k, v in out.items() if "gaps" in v for x in v["gaps"]]
print(len(g), "internal gaps;", sum(1 for _, x in g if x.get("n_lost")), "with residues present in the PDB entry")
