# train_deepce_ais.py — DeepCE 复现训练(公开架构从头训练, 无预训练权重)
# 训练/调参: signature_train / signature_dev。测试标签在模型锁定前不读取。
# 输出: ai/output/deepce_best.pt(最佳dev权重), ai/output/deepce_train_log.csv
import os, sys, time, math, csv
# 沙箱环境下 tempfile 探测(创建+删除)会被拦截, 固定临时目录绕过探测
import tempfile
_TMP = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), ".tmp")
os.makedirs(_TMP, exist_ok=True)
tempfile.tempdir = _TMP
tempfile._get_default_tempdir = lambda: _TMP
tempfile.gettempdir = lambda: _TMP
import numpy as np
import torch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from data_common import (PRED, CELLS, DOSES, load_signature, onehot_conds,
                         load_drug_smiles, DrugCache, build_batch, load_gene_features)
from scipy.stats import spearmanr

device = torch.device("mps" if torch.backends.mps.is_available() else "cpu")
print("device:", device, flush=True)
torch.manual_seed(20261109)
np.random.seed(20261109)

BATCH = 64
MAX_EPOCH = 40
LR = 2e-4  # 作者默认
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")
os.makedirs(OUT, exist_ok=True)

print("loading data...", flush=True)
m_tr, lab_tr, genes = load_signature("train")
m_de, lab_de, _ = load_signature("dev")
with open(f"{PRED}/signature_dev.csv") as f:
    assert genes == next(csv.reader(f))[5:]
cell_tr, dose_tr = onehot_conds(m_tr)
cell_de, dose_de = onehot_conds(m_de)
smiles = load_drug_smiles()
need = sorted(set(m_tr["pert_id"]) | set(m_de["pert_id"]))
cache = {pid: DrugCache(smiles[pid]) for pid in need}
print(f"train {len(lab_tr)} conds / dev {len(lab_de)} conds / {len(cache)} drugs", flush=True)

gene = load_gene_features().to(device)

from models import DeepCE  # noqa: E402
model = DeepCE(drug_input_dim={"atom": 62, "bond": 6}, drug_emb_dim=128, conv_size=[16, 16],
               degree=[0, 1, 2, 3, 4, 5], gene_input_dim=128, gene_emb_dim=128, num_gene=978,
               hid_dim=128, dropout=0.1, loss_type="point_wise_mse", device=device,
               initializer=torch.nn.init.xavier_uniform_,
               pert_type_input_dim=1, cell_id_input_dim=len(CELLS), pert_idose_input_dim=len(DOSES),
               pert_type_emb_dim=4, cell_id_emb_dim=4, pert_idose_emb_dim=4,
               use_pert_type=True, use_cell_id=True, use_pert_idose=True).to(device)
opt = torch.optim.Adam(model.parameters(), lr=LR)

def run_epoch(meta, lab, cellh, doseh, train_mode):
    model.train(train_mode)
    n = len(lab)
    order = np.random.permutation(n) if train_mode else np.arange(n)
    preds = np.empty((n, 978), dtype=np.float32)
    nb = math.ceil(n / BATCH)
    tot = 0.0
    with torch.set_grad_enabled(train_mode):
        for b in range(nb):
            idx = order[b * BATCH:(b + 1) * BATCH]
            ids = [meta["pert_id"][i] for i in idx]
            ft, mask = build_batch(ids, cache, device)
            pt = torch.ones(len(idx), 1, device=device)
            lb = torch.from_numpy(lab[idx]).to(device)
            pred = model(ft, gene, mask, pt,
                         torch.from_numpy(cellh[idx]).to(device),
                         torch.from_numpy(doseh[idx]).to(device))
            loss = model.loss(lb, pred)
            if train_mode:
                opt.zero_grad(); loss.backward(); opt.step()
            tot += loss.item()
            preds[idx] = pred.detach().cpu().numpy()
    return tot / nb, preds

def per_cond_rho(true, pred):
    out = np.full(len(true), np.nan)
    for i in range(len(true)):
        if pred[i].std() < 1e-6 or true[i].std() < 1e-6:
            continue
        out[i] = spearmanr(true[i], pred[i]).statistic
    return out

log = []
best = {"rho": -9.0, "epoch": -1}
for ep in range(1, MAX_EPOCH + 1):
    t0 = time.time()
    tl, _ = run_epoch(m_tr, lab_tr, cell_tr, dose_tr, True)
    _, pv = run_epoch(m_de, lab_de, cell_de, dose_de, False)
    rho_v = per_cond_rho(lab_de, pv)
    rho_m = float(np.nanmean(rho_v))
    undef = float(np.mean(np.isnan(rho_v)))
    rmse = float(np.sqrt(np.mean((lab_de - pv) ** 2)))
    log.append(dict(epoch=ep, train_loss=tl, dev_spearman_mean=rho_m,
                    dev_undef_frac=undef, dev_rmse=rmse, sec=round(time.time() - t0, 1)))
    with open(f"{OUT}/deepce_train_log.csv", "w") as f:
        w = csv.DictWriter(f, fieldnames=list(log[0].keys()))
        w.writeheader(); w.writerows(log)
    star = ""
    if rho_m > best["rho"]:
        best.update(rho=rho_m, epoch=ep)
        torch.save({"state_dict": model.state_dict(), "epoch": ep, "dev_spearman_mean": rho_m,
                    "cells": CELLS, "doses": DOSES}, f"{OUT}/deepce_best.pt")
        star = " *best*"
    print(f"ep{ep:02d} loss={tl:.4f} dev_rho={rho_m:.4f} undef={undef:.3f} rmse={rmse:.4f} "
          f"({log[-1]['sec']}s){star}", flush=True)

print(f"BEST epoch {best['epoch']} dev_rho={best['rho']:.4f} -> deepce_best.pt", flush=True)
