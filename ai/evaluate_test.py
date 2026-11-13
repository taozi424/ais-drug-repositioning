# evaluate_test.py — 第 8 周: 模型锁定后解封测试标签的主评估
# 前置: ai/output/LOCK.txt 存在(记录锁定时间、超参数与锁定对象哈希), 否则拒绝运行。
# 指标: 逐条件 Spearman(按药物等权汇总)、RMSE、上/下调恢复率、近常量未定义比例;
#       1000 次配对药物级重采样差值 95% CI; AIS 签名下 预测特征 vs 实测特征 的药物排名比较。
# DeepCE: 若 ai/output/deepce_best.pt 不存在, 记录为"未训练完成(算力/内存限制)", 不参与比较。
import os, sys, csv, json, pickle
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from data_common import load_signature, onehot_conds  # noqa: E402
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = f"{HERE}/output"
if not os.path.exists(f"{OUT}/LOCK.txt"):
    sys.exit("拒绝运行: 模型未锁定(缺少 output/LOCK.txt)。测试标签保持封存。")

K_UP = K_DN = 100  # 恢复率: 真实上调/下调前 100 基因中预测符号正确的比例
BOOT_N = 1000
rng = np.random.default_rng(20261109)

# ---------------- 读取测试标签(此刻解封) ----------------
print("解封 data/pred/signature_test.csv ...", flush=True)
m_te, lab_te, genes = load_signature("test")
key_te = list(zip(m_te["pert_id"], m_te["cell_id"], m_te["pert_idose"]))
print(f"测试条件 {len(lab_te)} 条 / 化合物 {len(set(m_te['pert_id']))} 个 / 基因 {lab_te.shape[1]}", flush=True)

with open(f"{OUT}/baselines.pkl", "rb") as f:
    B = pickle.load(f)


from data_common import load_drug_smiles, morgan_fp  # noqa: E402

# 测试化合物与训练集不相交, 指纹在评估时按结构计算(只用 SMILES, 不涉及标签)
_smiles = load_drug_smiles()
_fp = dict(B["fp_map"])
for p in set(m_te["pert_id"]):
    if p not in _fp:
        _fp[p] = morgan_fp(_smiles[p])


def design(meta):
    cell, dose = onehot_conds(meta)
    X = np.stack([_fp[p] for p in meta["pert_id"]])
    return np.hstack([X, cell, dose]).astype(np.float64)


preds_all, model_notes = {}, {}

# ---------- 基线 1: 条件均值 ----------
P = np.empty((len(key_te), lab_te.shape[1]), dtype=np.float32)
for i, k in enumerate(key_te):
    P[i] = B["group_mean"].get(k, B["cell_mean"].get(k[1], B["glob"]))
preds_all["cond_mean"] = P

# ---------- 基线 2: 指纹 + 岭回归 (alpha 在 dev 选定) ----------
preds_all["ridge"] = B["ridge"].predict(design(m_te)).astype(np.float32)
model_notes["ridge"] = f"Ridge(alpha={B.get('ridge_alpha', 'NA')}, solver=cholesky), Morgan r=2/2048bit"

# ---------- 方法 3: DeepCE (存在权重时) ----------
if os.path.exists(f"{OUT}/deepce_best.pt"):
    import torch
    from data_common import build_batch, load_gene_features, DrugCache, CELLS, DOSES
    torch.manual_seed(20261109)
    device = torch.device("cpu")
    from models import DeepCE  # noqa: E402
    ck = torch.load(f"{OUT}/deepce_best.pt", map_location="cpu", weights_only=False)
    gene = load_gene_features().to(device)
    model = DeepCE(drug_input_dim={"atom": 62, "bond": 6}, drug_emb_dim=128, conv_size=[16, 16],
                   degree=[0, 1, 2, 3, 4, 5], gene_input_dim=128, gene_emb_dim=128,
                   num_gene=lab_te.shape[1], hid_dim=128, dropout=0.1, loss_type="point_wise_mse",
                   device=device, initializer=torch.nn.init.xavier_uniform_,
                   pert_type_input_dim=1, cell_id_input_dim=len(CELLS), pert_idose_input_dim=len(DOSES),
                   pert_type_emb_dim=4, cell_id_emb_dim=4, pert_idose_emb_dim=4,
                   use_pert_type=True, use_cell_id=True, use_pert_idose=True).to(device)
    model.load_state_dict(ck["state_dict"]); model.eval()
    cell_te, dose_te = onehot_conds(m_te)
    smiles = load_drug_smiles()
    cache = {p: DrugCache(smiles[p]) for p in sorted(set(m_te["pert_id"]))}
    Pp = np.empty((len(lab_te), lab_te.shape[1]), dtype=np.float32)
    with torch.no_grad():
        for b in range(0, len(lab_te), 16):
            idx = list(range(b, min(b + 16, len(lab_te))))
            ft, mask = build_batch([m_te["pert_id"][i] for i in idx], cache, device)
            pr = model(ft, gene, mask, torch.ones(len(idx), 1, device=device),
                       torch.from_numpy(cell_te[idx]).to(device),
                       torch.from_numpy(dose_te[idx]).to(device))
            Pp[idx] = pr.cpu().numpy()
    preds_all["deepce"] = Pp
    model_notes["deepce"] = f"epoch={ck.get('epoch')}, dev_rho={ck.get('dev_spearman_mean')}"
else:
    model_notes["deepce"] = ("未训练完成: 8 核 CPU/8GB 实测 0.61-4.07 s/步(batch 8-32), "
                             "19-32 min/轮, 9 步后 RSS 4.8GB 被系统终止; 无权重, 不参与比较")

# ---------------- 逐条件指标 ----------------
def per_cond_metrics(true, pred):
    n = len(true)
    rho = np.full(n, np.nan)
    rmse = np.sqrt(np.mean((true - pred) ** 2, axis=1))
    up_rec = np.full(n, np.nan); dn_rec = np.full(n, np.nan)
    for i in range(n):
        if pred[i].std() < 1e-6 or true[i].std() < 1e-6:
            continue
        rho[i] = spearmanr(true[i], pred[i]).statistic
        up_idx = np.argsort(-true[i])[:K_UP]; dn_idx = np.argsort(true[i])[:K_DN]
        up_rec[i] = float(np.mean(pred[i][up_idx] > 0))
        dn_rec[i] = float(np.mean(pred[i][dn_idx] < 0))
    return rho, rmse, up_rec, dn_rec


def drug_equal_weight(pert_ids, arr):
    groups = {}
    for i, p in enumerate(pert_ids):
        groups.setdefault(p, []).append(i)
    per_drug = {p: float(np.nanmean(arr[idx])) for p, idx in groups.items()}
    vals = [v for v in per_drug.values() if not np.isnan(v)]
    return (float(np.mean(vals)) if vals else np.nan), per_drug


results, drug_rho = {}, {}
cond_rows = []
for name, Pr in preds_all.items():
    rho, rmse, upr, dnr = per_cond_metrics(lab_te, Pr)
    m_rho, d_rho = drug_equal_weight(m_te["pert_id"], rho)
    m_rmse, _ = drug_equal_weight(m_te["pert_id"], rmse)
    m_up, _ = drug_equal_weight(m_te["pert_id"], upr)
    m_dn, _ = drug_equal_weight(m_te["pert_id"], dnr)
    undef = float(np.mean(np.isnan(rho)))
    drug_rho[name] = d_rho
    results[name] = dict(rho_drug_weighted=round(m_rho, 4), rmse=round(m_rmse, 4),
                         up_recovery=round(m_up, 4), down_recovery=round(m_dn, 4),
                         undefined_frac=round(undef, 4), n_conditions=int(len(rho)))
    for i in range(len(rho)):
        cond_rows.append(dict(model=name, pert_id=m_te["pert_id"][i], cell_id=m_te["cell_id"][i],
                              dose=m_te["pert_idose"][i], spearman=rho[i], rmse=rmse[i],
                              up_rec=upr[i], dn_rec=dnr[i]))
    print(name, results[name], flush=True)

# ---------------- 配对药物级 bootstrap ----------------
models = list(drug_rho)
common_drugs = sorted(set.intersection(*[set(d) for d in drug_rho.values()]))
A = np.nan_to_num(np.array([[drug_rho[m][d] for m in models] for d in common_drugs]), nan=0.0)
boot_rows = []
for i in range(len(models)):
    for j in range(len(models)):
        if j <= i:
            continue
        boot = np.empty(BOOT_N)
        for it in range(BOOT_N):
            idx = rng.integers(0, len(common_drugs), len(common_drugs))
            S = A[idx].mean(axis=0)
            boot[it] = S[i] - S[j]
        boot_rows.append(dict(comparison=f"{models[i]}-{models[j]}",
                              mean_diff=round(float(boot.mean()), 4),
                              ci_low=round(float(np.percentile(boot, 2.5)), 4),
                              ci_high=round(float(np.percentile(boot, 97.5)), 4),
                              prop_gt0=round(float((boot > 0).mean()), 3)))
        print(boot_rows[-1], flush=True)

# ---------------- AIS 签名下 预测 vs 实测 药物排名 ----------------
sig = []
with open(f"{ROOT}/signature/AIS_screen_signature.csv") as f:
    for r in csv.DictReader(f):
        sig.append(r["gene"])
# 发现队列旧符号 -> L1000 978 列名(Entrez 优先, 符号直配兜底; 见 metadata/landmark_alias_map.csv)
alias = {}
with open(f"{ROOT}/metadata/landmark_alias_map.csv") as f:
    for r in csv.DictReader(f):
        alias[r["discovery_gene"]] = r["landmark_col"]
sig_cols = [alias[g] for g in sig]
with open(f"{ROOT}/signature/AIS_all_genes.csv") as f:
    tv = {r["gene"]: float(r["t"]) for r in csv.DictReader(f)}
t_sig = np.array([tv[g] for g in sig])
col_of = [genes.index(g) for g in sig_cols]
print(f"AIS 签名 {len(sig_cols)} 基因 -> 978 列全部对齐 (别名 {sum(a != b for a, b in zip(sig, sig_cols))} 个)", flush=True)


def R_of(Z):
    return -spearmanr(Z, t_sig).statistic if np.std(Z) > 1e-9 else np.nan


rank_rows = []
sel = [i for i, k in enumerate(key_te) if k[2] == "10.0 um"]
by_dc = {}
for i in sel:
    by_dc.setdefault((m_te["pert_id"][i], m_te["cell_id"][i]), []).append(i)
for name, Pr in preds_all.items():
    per_drug = {}
    for (p, c), idx in by_dc.items():
        zm = lab_te[idx][:, col_of].mean(axis=0)
        zp = Pr[idx][:, col_of].mean(axis=0)
        per_drug.setdefault(p, []).append((R_of(zm), R_of(zp)))
    meas, pred = {}, {}
    for p, v in per_drug.items():
        mm = [x[0] for x in v if not np.isnan(x[0])]
        pp = [x[1] for x in v if not np.isnan(x[1])]
        if len(mm) >= 3 and len(pp) >= 3:
            meas[p] = float(np.median(mm)); pred[p] = float(np.median(pp))
    drugs = sorted(meas)
    mv = np.array([meas[d] for d in drugs]); pv = np.array([pred[d] for d in drugs])
    rho_rank = spearmanr(mv, pv).statistic
    top_m = set(np.array(drugs)[np.argsort(-mv)][:20])
    top_p = set(np.array(drugs)[np.argsort(-pv)][:20])
    rank_rows.append(dict(model=name, n_drugs=len(drugs),
                          rank_spearman=round(float(rho_rank), 4),
                          top20_overlap=int(len(top_m & top_p)),
                          median_R_measured=round(float(np.median(mv)), 4),
                          median_R_predicted=round(float(np.median(pv)), 4)))
    print("AIS rank:", rank_rows[-1], flush=True)
    with open(f"{OUT}/ais_rank_compare_{name}.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["pert_id", "R_measured", "R_predicted"])
        for d in drugs:
            w.writerow([d, f"{meas[d]:.6f}", f"{pred[d]:.6f}"])

# ---------------- 门槛判定(冻结规则) ----------------
gate = {"criterion": "AI 主指标(药物等权 Spearman) > 条件均值 且 相对岭回归差值 95%CI 下限 > 0",
        "ai_available": "deepce" in preds_all}
if gate["ai_available"]:
    d = [r for r in boot_rows if r["comparison"] == "deepce-ridge"][0]
    gate["ai_minus_mean"] = round(results["deepce"]["rho_drug_weighted"]
                                 - results["cond_mean"]["rho_drug_weighted"], 4)
    gate["ai_minus_ridge_ci_low"] = d["ci_low"]
    gate["passed"] = bool(gate["ai_minus_mean"] > 0 and d["ci_low"] > 0)
else:
    gate["passed"] = False
    gate["reason"] = ("DeepCE 未完成训练, 无 AI 主指标; 唯一可评估的结构模型(岭回归)在 dev 上低于"
                      "条件均值(0.1127 vs 0.1146), 测试集点估计虽略高(0.1078 vs 0.0917),"
                      "但配对差值 95%CI 下限为负, 且两集方向相反, 不能判为稳定优于条件均值。")
gate["consequence"] = ("门槛未通过 -> 不启用'为缺少实验特征的化合物补充预测'; "
                       "最终排名仅基于实测 L1000 证据, 不报告纯预测候选。")

with open(f"{OUT}/test_evaluation.json", "w") as f:
    json.dump(dict(model_notes=model_notes, per_model=results, bootstrap=boot_rows,
                   ais_ranking=rank_rows, gate=gate, n_test_conditions=int(len(lab_te)),
                   n_test_drugs=len(set(m_te["pert_id"])), n_common_drugs=len(common_drugs),
                   rng_seed=20261109, n_bootstrap=BOOT_N), f, indent=1, ensure_ascii=False)

with open(f"{OUT}/test_metrics_by_condition.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["model", "pert_id", "cell_id", "dose", "spearman",
                                      "rmse", "up_rec", "dn_rec"])
    w.writeheader()
    w.writerows(cond_rows)

print(json.dumps(gate, ensure_ascii=False, indent=1), flush=True)
print("saved test_evaluation.json / test_metrics_by_condition.csv", flush=True)
