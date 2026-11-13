#!/usr/bin/env python3
# =====================================================================
# 07 药物开发状态注释 (ChEMBL WebResource)
# 匹配层级: 1) standard_inchi_key 逐个查询  2) canonical_smiles 精确
#           3) 归一化 pref_name
# 输出: metadata/compound_development_status.csv
# =====================================================================
import json, urllib.request, time, csv, os, datetime, re
from concurrent.futures import ThreadPoolExecutor, as_completed

BASE = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"
OUT = "metadata/compound_development_status.csv"
RAW = "metadata/compound_perturbagens_raw.csv"
CACHE = "metadata/chembl_phase_gte1_cache.json"

def fetch(url, tries=4):
    for t in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "ais-repurposing/1.0", "Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.loads(r.read().decode())
        except Exception as e:
            print(f"retry {t+1}: {e}", flush=True)
            time.sleep(8)
    return None

def norm(s):
    s = (s or "").lower().strip()
    s = re.sub(r"[^a-z0-9]+", "-", s).strip("-")
    return s

# ---- 1) 缓存 max_phase>=1 分子 (smiles + pref_name + 状态) ----
if os.path.exists(CACHE):
    rows = json.load(open(CACHE))
    print(f"cache loaded: {len(rows)}")
else:
    mols, url, page = [], f"{BASE}?max_phase__gte=1&limit=1000", 1
    while True:
        d = fetch(url)
        if d is None:
            raise RuntimeError("fetch failed")
        mols.extend(d["molecules"])
        nxt = d.get("page_meta", {}).get("next")
        print(f"page {page}: {len(mols)}", flush=True)
        if not nxt:
            break
        url = "https://www.ebi.ac.uk" + nxt if nxt.startswith("/") else nxt
        page += 1
        time.sleep(0.5)
    rows = []
    for m in mols:
        st = m.get("molecule_structures") or {}
        rows.append({"chembl_id": m.get("molecule_chembl_id"),
                     "pref_name": (m.get("pref_name") or "").strip(),
                     "max_phase": m.get("max_phase"),
                     "first_approval": m.get("first_approval"),
                     "canonical_smiles": st.get("canonical_smiles") or ""})
    json.dump(rows, open(CACHE, "w"))
    print(f"cache saved: {len(rows)}")

smiles2mol = {}
name2mol = {}
for r in rows:
    if r["canonical_smiles"]:
        smiles2mol.setdefault(r["canonical_smiles"], r)
    if r["pref_name"]:
        name2mol.setdefault(norm(r["pref_name"]), r)

# ---- 2) trt_cp 化合物 ----
cp = list(csv.DictReader(open(RAW)))
cp = [r for r in cp if r["pert_type"] == "trt_cp"]
today = datetime.date.today().isoformat()

def lookup_by_key(ik):
    d = fetch(f"{BASE}?molecule_structures__standard_inchi_key={urllib.parse.quote(ik)}"
              f"&limit=5&only=molecule_chembl_id,pref_name,max_phase,first_approval")
    if d and d["molecules"]:
        m = d["molecules"][0]
        return {"chembl_id": m["molecule_chembl_id"], "pref_name": m.get("pref_name") or "",
                "max_phase": m.get("max_phase"), "first_approval": m.get("first_approval")}
    return None

import urllib.parse
KEYCACHE = "metadata/chembl_key_lookup_cache.json"
keycache = json.load(open(KEYCACHE)) if os.path.exists(KEYCACHE) else {}

results = [None] * len(cp)

def work(i_r):
    i, r = i_r
    ik = r["inchi_key"] if r["inchi_key"] != "-666" else ""
    mol, how = None, "none"
    if ik:
        if ik in keycache:
            mol = keycache[ik]
        else:
            mol = lookup_by_key(ik)
            keycache[ik] = mol
            if len(keycache) % 50 == 0:
                json.dump(keycache, open(KEYCACHE, "w"))
        if mol:
            how = "inchi_key_full"
    if mol is None and r["canonical_smiles"] not in ("-666", "") and r["canonical_smiles"] in smiles2mol:
        mol, how = smiles2mol[r["canonical_smiles"]], "canonical_smiles"
    if mol is None and norm(r["pert_iname"]) in name2mol:
        mol, how = name2mol[norm(r["pert_iname"])], "pref_name_normalized"
    return i, mol, how

with ThreadPoolExecutor(max_workers=4) as ex:
    futs = [ex.submit(work, x) for x in enumerate(cp)]
    done = 0
    for fu in as_completed(futs):
        i, mol, how = fu.result()
        results[i] = (mol, how)
        done += 1
        if done % 200 == 0:
            json.dump(keycache, open(KEYCACHE, "w"))
            print(f"{done}/{len(cp)}", flush=True)
json.dump(keycache, open(KEYCACHE, "w"))

n_full = n_smiles = n_name = n_research = 0
out = []
for r, (mol, how) in zip(cp, results):
    if how == "inchi_key_full":
        n_full += 1
    elif how == "canonical_smiles":
        n_smiles += 1
    elif how == "pref_name_normalized":
        n_name += 1
    else:
        n_research += 1
    out.append({"pert_id": r["pert_id"], "pert_iname": r["pert_iname"], "inchi_key": r["inchi_key"],
                "matched_by": how,
                "chembl_id": mol["chembl_id"] if mol else "",
                "pref_name": mol["pref_name"] if mol else "",
                "max_phase": mol["max_phase"] if mol else "0_or_research",
                "first_approval": mol["first_approval"] if mol else "",
                "status_source": f"ChEMBL WebResource molecule.json (max_phase>=1 cache + per-key lookup), checked {today}"})

with open(OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(out[0].keys()))
    w.writeheader(); w.writerows(out)

n_cp = len(out)
def ph(o):
    try:
        return int(float(o["max_phase"]))
    except Exception:
        return 0
n_marketed = sum(1 for o in out if ph(o) >= 4)
n_clinical = sum(1 for o in out if 1 <= ph(o) <= 3)
print(f"trt_cp {n_cp}: InChIKey {n_full}, smiles {n_smiles}, pref_name {n_name}, 未匹配 {n_research}")
print(f"上市(max_phase>=4) {n_marketed}, 临床(1-3期) {n_clinical}, 研究/未核实 {n_cp - n_marketed - n_clinical}")
print(f"写入 {OUT}")
