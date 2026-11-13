#!/usr/bin/env python3
# =====================================================================
# 18 重点候选的药物—靶点/MOA 可追溯证据 (ChEMBL mechanism 端点)
# 仅收录 ChEMBL curated mechanism_of_action 记录 (含文献来源),
# 不做任何推断; 无记录者如实记"无可追溯记录"。
# 输出 results/tables/featured_target_evidence.csv
# =====================================================================
import json, subprocess, urllib.parse
import csv, os

# 切到仓库根目录（脚本位于 scripts/ 下），也可用 AIS_ROOT 覆盖
ROOT = os.environ.get("AIS_ROOT") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)
status = {}
with open("metadata/compound_development_status.csv") as f:
    for r in csv.DictReader(f):
        status[r["pert_iname"].lower()] = r

feat = []
with open("results/tables/featured_candidates.csv") as f:
    for r in csv.DictReader(f):
        feat.append(r)

rows = []
for fc in feat:
    iname = fc["iname"].lower()
    st = status.get(iname)
    chembl_id = st["chembl_id"] if st and st.get("chembl_id") else None
    if not chembl_id:
        rows.append(dict(iname=fc["iname"], chembl_id="", moa="", action_type="",
                         target_chembl_id="", source="ChEMBL mechanism.json", checked="2026-11-10",
                         note="开发状态表中无 chembl_id"))
        continue
    url = f"https://www.ebi.ac.uk/chembl/api/data/mechanism.json?molecule_chembl_id={chembl_id}&format=json"
    try:
        out = subprocess.run(["curl", "-sSL", "--max-time", "60", url],
                             capture_output=True, text=True, timeout=90)
        d = json.loads(out.stdout)
        mechs = d.get("mechanisms", [])
        if not mechs:
            rows.append(dict(iname=fc["iname"], chembl_id=chembl_id, moa="", action_type="",
                             target_chembl_id="", source="ChEMBL mechanism.json", checked="2026-11-10",
                             note="ChEMBL 无 curated MOA 记录"))
        else:
            for m in mechs:
                rows.append(dict(iname=fc["iname"], chembl_id=chembl_id,
                                 moa=m.get("mechanism_of_action", ""),
                                 action_type=m.get("action_type", ""),
                                 target_chembl_id=m.get("target_chembl_id", ""),
                                 source="ChEMBL mechanism.json (curated, 含文献)", checked="2026-11-10",
                                 note=""))
    except Exception as e:
        rows.append(dict(iname=fc["iname"], chembl_id=chembl_id, moa="", action_type="",
                         target_chembl_id="", source="ChEMBL mechanism.json", checked="2026-11-10",
                         note=f"查询失败: {e}"))

with open("results/tables/featured_target_evidence.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["iname", "chembl_id", "moa", "action_type",
                                      "target_chembl_id", "source", "checked", "note"])
    w.writeheader()
    w.writerows(rows)
for r in rows:
    print(r["iname"], "|", r["moa"][:60], "|", r["action_type"], "|", r["target_chembl_id"], "|", r["note"])
