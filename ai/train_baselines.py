# train_baselines.py — 两个基线: 训练集条件均值 / Morgan指纹+条件one-hot + 岭回归
# 只用 signature_train 训练, signature_dev 调参(alpha), 测试标签不读取。
import os, sys, pickle
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from data_common import load_signature, onehot_conds, load_drug_smiles, morgan_fp
from sklearn.linear_model import Ridge
from scipy.stats import spearmanr

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")
os.makedirs(OUT, exist_ok=True)
np.random.seed(20261109)

def per_cond_rho(true, pred):
    out = np.full(len(true), np.nan)
    for i in range(len(true)):
        if pred[i].std() < 1e-6 or true[i].std() < 1e-6:
            continue
        out[i] = spearmanr(true[i], pred[i]).statistic
    return out

print("loading...", flush=True)
m_tr, lab_tr, _ = load_signature("train")
m_de, lab_de, _ = load_signature("dev")
key_tr = list(zip(m_tr["pert_id"], m_tr["cell_id"], m_tr["pert_idose"]))
key_de = list(zip(m_de["pert_id"], m_de["cell_id"], m_de["pert_idose"]))

# ---------- 基线1: 训练集条件均值 (cell,dose) -> cell 均值 -> 全局均值 ----------
group_mean, cell_sum, cell_n = {}, {}, {}
for i, k in enumerate(key_tr):
    g = group_mean.setdefault(k, [np.zeros(lab_tr.shape[1], np.float64), 0])
    g[0] += lab_tr[i]; g[1] += 1
    c = cell_sum.setdefault(k[1], [np.zeros(lab_tr.shape[1], np.float64), 0])
    c[0] += lab_tr[i]; c[1] += 1
group_mean = {k: (v[0] / v[1]).astype(np.float32) for k, v in group_mean.items()}
cell_mean = {k: (v[0] / v[1]).astype(np.float32) for k, v in cell_sum.items()}
glob = lab_tr.mean(axis=0)

def predict_mean(keys):
    P = np.empty((len(keys), lab_tr.shape[1]), dtype=np.float32)
    for i, k in enumerate(keys):
        P[i] = group_mean.get(k, cell_mean.get(k[1], glob))
    return P

Pm_de = predict_mean(key_de)
rho_m = per_cond_rho(lab_de, Pm_de)
print(f"[cond-mean] dev rho_mean={np.nanmean(rho_m):.4f} undef={np.mean(np.isnan(rho_m)):.3f} "
      f"rmse={np.sqrt(np.mean((lab_de - Pm_de) ** 2)):.4f}", flush=True)

# ---------- 基线2: Morgan指纹 + 条件 one-hot -> 岭回归 ----------
smiles = load_drug_smiles()
uids = sorted(set(m_tr["pert_id"]) | set(m_de["pert_id"]))
fp_map = {u: morgan_fp(smiles[u]) for u in uids}
print("fingerprints done", flush=True)

def design(meta):
    cell, dose = onehot_conds(meta)
    X = np.stack([fp_map[p] for p in meta["pert_id"]])
    return np.hstack([X, cell, dose]).astype(np.float64)

Xtr, Xde = design(m_tr), design(m_de)
print("design:", Xtr.shape, Xde.shape, flush=True)

best = {"alpha": None, "rho": -9}
for a in [0.1, 1.0, 10.0, 100.0]:
    r = Ridge(alpha=a, solver="cholesky", fit_intercept=True)
    r.fit(Xtr, lab_tr)
    P = r.predict(Xde)
    rho = per_cond_rho(lab_de, P.astype(np.float32))
    m = float(np.nanmean(rho))
    print(f"[ridge a={a}] dev rho_mean={m:.4f} undef={np.mean(np.isnan(rho)):.3f}", flush=True)
    if m > best["rho"]:
        best = {"alpha": a, "rho": m}
print("best alpha:", best, flush=True)

ridge = Ridge(alpha=best["alpha"], solver="cholesky", fit_intercept=True).fit(Xtr, lab_tr)

with open(f"{OUT}/baselines.pkl", "wb") as f:
    pickle.dump({
        "ridge": ridge, "ridge_alpha": best["alpha"],
        "group_mean": group_mean, "cell_mean": cell_mean, "glob": glob,
        "fp_map": fp_map,
    }, f)
print("saved baselines.pkl", flush=True)
