# data_common.py — 三方法共享的数据加载与特征模块(无 pandas 依赖)
# 规则: 与 DeepCE 作者一致的重复条件代表选择(choose_mean_example)、one-hot 条件编码、
#       分子图特征预缓存(等价于作者 convert_smile_to_feature, 但每药物只算一次)。
import csv
import os
import sys
import numpy as np
import torch
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]          # 仓库根目录，可用 AIS_ROOT 覆盖
DEEPCEDIR = str(Path(os.environ.get("AIS_ROOT", ROOT)) / "ai" / "DeepCE")
if DEEPCEDIR not in sys.path:
    sys.path.insert(0, DEEPCEDIR)

PRED = str(Path(os.environ.get("AIS_ROOT", ROOT)) / "data" / "pred")
META_COLS = ("sig_id", "pert_id", "pert_type", "cell_id", "pert_idose")

CELLS = sorted(["A375", "HA1E", "HELA", "HT29", "MCF7", "PC3", "YAPC"])
DOSES = sorted(["0.04 um", "0.12 um", "0.37 um", "1.11 um", "3.33 um", "10.0 um"])  # 与作者 DataReader 相同的 sorted 顺序


def load_signature(split):
    """split: 'train' | 'dev'。读取并按作者规则对重复条件去重。
    返回 meta(dict of lists), label [n,978] float32, genes(list)。"""
    fp = f"{PRED}/signature_{split}.csv"
    meta = {k: [] for k in META_COLS}
    rows = []
    with open(fp) as f:
        r = csv.reader(f)
        header = next(r)
        genes = header[5:]
        for line in r:
            if not line:
                continue
            for i, k in enumerate(META_COLS):
                meta[k].append(line[i])
            rows.append(line[5:])
    lab = np.asarray(rows, dtype=np.float32)
    del rows
    key = [f"{p}|{c}|{d}" for p, c, d in zip(meta["pert_id"], meta["cell_id"], meta["pert_idose"])]
    groups = {}
    for i, k in enumerate(key):
        groups.setdefault(k, []).append(i)
    keep = []
    for k, idx in groups.items():
        if len(idx) == 1:
            keep.append(idx[0])
        else:
            ex = lab[idx]
            rk = np.argsort(np.argsort(ex, axis=0), axis=0).mean(axis=1)  # 作者 choose_mean_example 规则
            med = (len(idx) - 1) / 2
            keep.append(idx[int(np.argmin((rk - med) ** 2))])
    keep = sorted(keep)
    for k in META_COLS:
        meta[k] = [meta[k][i] for i in keep]
    lab = lab[keep]
    return meta, lab, genes


def onehot_conds(meta):
    """与作者 DataReader 相同: one-hot(sorted cell set), one-hot(sorted dose set)。"""
    n = len(meta["cell_id"])
    cell = np.zeros((n, len(CELLS)), dtype=np.float32)
    cell[np.arange(n), [CELLS.index(c) for c in meta["cell_id"]]] = 1
    dose = np.zeros((n, len(DOSES)), dtype=np.float32)
    dose[np.arange(n), [DOSES.index(d) for d in meta["pert_idose"]]] = 1
    return cell, dose


def load_drug_smiles():
    out = {}
    with open(f"{PRED}/drugs_smiles.csv") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            pid, smi = line.split(",", 1)
            out[pid] = smi
    return out


# ---------------- 分子图缓存(与作者 atom_features/bond_features 完全一致) ----------------
from utils.molecule_utils import atom_features, bond_features  # noqa: E402  (DeepCE utils)


class DrugCache:
    """每药物一次: 度排序的原子/键特征数组 + 邻接索引(全局索引由 batch 组装时加偏移)。"""

    def __init__(self, smiles):
        from utils.molecules import Molecule, degrees
        mol = Molecule(smiles)  # 复用作者的解析与度排序
        self.atom_feat = np.asarray([n.data for n in mol.get_node_list("atom")], dtype=np.float32)
        self.bond_feat = np.asarray([n.data for n in mol.get_node_list("bond")], dtype=np.float32)
        atom_idx = {n.ext_id: i for i, n in enumerate(mol.get_node_list("atom"))}
        bond_idx = {n.ext_id: i for i, n in enumerate(mol.get_node_list("bond"))}
        atoms = mol.get_node_list("atom")  # 度排序后的原子列表
        self.atom_nb_atom = [sorted(atom_idx[n.ext_id] for n in a.get_neighbors() if n.ntype == "atom")
                             for a in atoms]
        self.atom_nb_bond = [sorted(bond_idx[n.ext_id] for n in a.get_neighbors() if n.ntype == "bond")
                             for a in atoms]
        self.degree_groups = {}
        for d in degrees:
            self.degree_groups[d] = [i for i, n in enumerate(atoms)
                                     if len([x for x in n.get_neighbors() if x.ntype == "atom"]) == d]
        self.n_atom, self.n_bond = len(atoms), len(mol.get_node_list("bond"))


class FastBatch:
    """与作者 Molecules(batch) 等价的批组装, 但只做索引偏移, 不重新解析 SMILES。
    提供 DeepCE/NeuralFingerprint 用到的 3 个接口。"""

    def __init__(self, keys, cache, device):
        from utils.molecules import degrees
        self.batch_size = len(keys)
        atom_parts, bond_parts = [], []
        batch_atom_idx = []
        nb = {d: {"node": [], "edge": []} for d in degrees}
        aoff = boff = 0
        for k in keys:
            dc = cache[k]
            aoff0, boff0 = aoff, boff
            atom_parts.append(dc.atom_feat)
            bond_parts.append(dc.bond_feat)
            batch_atom_idx.append(list(range(aoff0, aoff0 + dc.n_atom)))
            for d in degrees:
                nb[d]["node"] += [[i + aoff0 for i in dc.atom_nb_atom[i]] for i in dc.degree_groups[d]]
                nb[d]["edge"] += [[j + boff0 for j in dc.atom_nb_bond[i]] for i in dc.degree_groups[d]]
            aoff += dc.n_atom
            boff += dc.n_bond
        self.atom = torch.from_numpy(np.concatenate(atom_parts)).to(device)
        self.bond = torch.from_numpy(np.concatenate(bond_parts)).to(device)
        self._nb = nb
        self._batch_atom_idx = batch_atom_idx

    def get_neighbor_idx_by_batch(self, ntype):
        assert ntype == "atom"
        return self._batch_atom_idx

    def get_neighbor_idx_by_degree(self, ntype, degree):
        return self._nb[degree]["node" if ntype == "atom" else "edge"]


def build_batch(pert_ids, cache, device):
    mb = FastBatch(pert_ids, cache, device)
    mask = torch.zeros(len(pert_ids), max(len(x) for x in mb._batch_atom_idx), device=device)
    for i, idxs in enumerate(mb._batch_atom_idx):
        mask[i, :len(idxs)] = 1
    return {"molecules": mb, "atom": mb.atom, "bond": mb.bond}, mask


# ---------------- Morgan 指纹(岭回归用) ----------------
_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def morgan_fp(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"bad smiles: {smiles}")
    return np.asarray(_gen.GetFingerprint(mol), dtype=np.float32)


def load_gene_features():
    rows = []
    with open(f"{PRED}/gene_vector.csv") as f:
        for line in f:
            line = line.strip()
            if line:
                rows.append([float(x) for x in line.split(",")[1:]])
    return torch.from_numpy(np.asarray(rows, dtype=np.float32))  # [978,128]
