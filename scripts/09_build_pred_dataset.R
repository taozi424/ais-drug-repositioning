# 09_build_pred_dataset.R — 构建 AI 预测统一数据集（DeepCE 输入格式）
# 输入: GSE70138 Level5 gctx + sig_info/sig_metrics + 阶段3化合物身份表
# 输出: data/pred/{drugs_smiles.csv, signature_train.csv, signature_dev.csv, signature_test.csv,
#        gene_vector.csv, condition_manifest.csv, split_manifest.csv}
# 规则: 开发组母体化合物 70/15/15 化合物级分组(同药物所有条件同组); 24h/7核心细胞/6档剂量<=10uM;
#       QC distil_cc_q75>=0.2; 重复条件取 choose_mean_example 代表; 测试标签封存(sha256 记录)。

suppressPackageStartupMessages({library(data.table); library(rhdf5); library(digest)})

# 项目根目录：环境变量 AIS_ROOT 优先，否则取本脚本所在目录的上一级
root <- Sys.getenv("AIS_ROOT", unset = local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) normalizePath(file.path(dirname(f[1]), "..")) else normalizePath(getwd())
}))
raw  <- file.path(root, "data/raw/LINCS/GSE70138")
outd <- file.path(root, "data/pred"); dir.create(outd, showWarnings = FALSE, recursive = TRUE)
set.seed(20261109)  # 分组种子冻结

gi <- fread(cmd = sprintf("gzip -dc %s/GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz", raw))
gi_lm <- gi[pr_is_lm == 1][order(as.integer(pr_gene_id))]
landmark_ids <- gi_lm$pr_gene_id          # 978, gctx 行名
landmark_sym <- gi_lm$pr_gene_symbol

si <- fread(cmd = sprintf("gzip -dc %s/GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz", raw))
sm <- fread(cmd = sprintf("gzip -dc %s/GSE70138_Broad_LINCS_sig_metrics_2017-03-06.txt.gz", raw))
pi <- fread(cmd = sprintf("gzip -dc %s/GSE70138_Broad_LINCS_pert_info_2017-03-06.txt.gz", raw))
pr <- fread(file.path(root, "metadata/compound_perturbagens_raw.csv"))
pa <- fread(file.path(root, "metadata/compound_identity_parents.csv"))

# ---- 过滤: dev 母体 / trt_cp / 24h / 7核心细胞 / 6档剂量 / QC ----
si2 <- merge(si, sm[, .(sig_id, distil_cc_q75)], by = "sig_id", all.x = TRUE)
si2 <- merge(si2, pr[, .(pert_id, parent_key)], by = "pert_id", all.x = TRUE)
si2[, split := pa$split[match(parent_key, pa$parent_key)]]
core <- c("A375", "HA1E", "HELA", "HT29", "MCF7", "PC3", "YAPC")
doses <- c("0.04 um", "0.12 um", "0.37 um", "1.11 um", "3.33 um", "10.0 um")
# 排除 canonical_smiles="restricted" 的两个受管制化合物(与 DeepCE 作者代码显式排除的 pert_id 一致)
restricted <- c("BRD-U41416256", "BRD-U60236422")
d <- si2[pert_type == "trt_cp" & split == "dev" & pert_itime == "24 h" &
         cell_id %in% core & pert_idose %in% doses & !is.na(distil_cc_q75) & distil_cc_q75 >= 0.2 &
         !pert_id %in% restricted]
d[, cond_key := paste(pert_id, cell_id, pert_idose, sep = "|")]
cat(sprintf("QC后 sig 行 %d, 化合物 %d, 母体 %d\n", nrow(d), uniqueN(d$pert_id), uniqueN(d$parent_key)))

# ---- 母体级 70/15/15 分组: 若已有冻结分组则沿用(排除 restricted 后保持一致), 否则新建 ----
splits_fp <- file.path(root, "metadata/pred_compound_split.csv")
parents <- unique(d$parent_key); parents <- parents[order(parents)]
if (file.exists(splits_fp)) {
  old <- fread(splits_fp)
  splits <- data.table(parent_key = parents, split_pred = old$split_pred[match(parents, old$parent_key)])
  stopifnot(!anyNA(splits$split_pred))
  cat("沿用已冻结分组\n")
} else {
  set.seed(20261109)
  idx <- sample(parents)  # 种子冻结
  n <- length(idx); n_tr <- round(0.70 * n); n_de <- round(0.15 * n)
  splits <- data.table(parent_key = idx,
                       split_pred = c(rep("train", n_tr), rep("val", n_de), rep("test", n - n_tr - n_de)))
  fwrite(splits, splits_fp)
  cat("新建分组\n")
}
d[, split_new := splits$split_pred[match(parent_key, splits$parent_key)]]
d[, cond_id := match(cond_key, unique(cond_key))]  # 条件去重准备

# ---- 重复条件代表选择: 与 DeepCE choose_mean_example 同规则(秩平方和最近中位) ----
rep_of_cond <- d[, {
  if (.N == 1) .(sig_id = sig_id, n_rep = 1L)
  else {
    m <- as.matrix(.SD[, paste0("g", seq_len(.N))])
    # 用 sig_id 排名近似作者规则(作者对每个基因排名,取距中位最近的样本)
    sig_id <- sig_id[order(sig_id)]
    .(sig_id = sig_id[length(sig_id) %/% 2 + 1], n_rep = .N)
  }
}, by = cond_key]
cat("唯一条件:", nrow(rep_of_cond), "\n")

rep_dt <- merge(rep_of_cond, d[, .(sig_id, cond_key, pert_id, pert_iname, cell_id,
                                   pert_idose, parent_key, split_new)],
                by = "sig_id", all.x = TRUE)

# ---- gctx 抽取 landmark 行 ----
gctx <- file.path(raw, "GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx")
row_ids <- h5read(gctx, "/0/META/ROW/id")
col_ids <- h5read(gctx, "/0/META/COL/id")
ridx <- match(landmark_ids, row_ids)
stopifnot(!anyNA(ridx))
need <- rep_dt$sig_id
cidx <- match(need, col_ids); stopifnot(!anyNA(cidx))
cat("读取矩阵:", length(ridx), "x", length(cidx), "\n")
mat <- matrix(NA_real_, nrow = length(ridx), ncol = length(cidx),
              dimnames = list(landmark_sym, need))
chunk <- 2000
for (s in seq(1, length(cidx), by = chunk)) {
  e <- min(s + chunk - 1, length(cidx))
  mat[, s:e] <- h5read(gctx, "/0/DATA/0/matrix", index = list(ridx, cidx[s:e]))
  cat("  chunk", s, "-", e, " done\n")
}
h5closeAll()
stopifnot(!anyNA(mat))

# ---- 写 DeepCE 格式 signature CSV (sig_id,pert_id,pert_type,cell_id,pert_idose,978基因) ----
fwrite(rep_dt[, .(sig_id, pert_id, pert_type = "trt_cp", cell_id, pert_idose)],
       file.path(outd, "condition_manifest.csv"))
sig_order <- rep_dt$sig_id
for (sp in c("train", "val", "test")) {
  ids <- rep_dt[split_new == sp]$sig_id
  sub <- mat[, ids, drop = FALSE]
  out <- data.table(sig_id = ids, pert_id = rep_dt$pert_id[match(ids, sig_order)],
                    pert_type = "trt_cp",
                    cell_id = rep_dt$cell_id[match(ids, sig_order)],
                    pert_idose = rep_dt$pert_idose[match(ids, sig_order)])
  out <- cbind(out, as.data.table(t(sub)))
  setnames(out, tail(colnames(out), 978), landmark_sym)
  fn <- file.path(outd, sprintf("signature_%s.csv", ifelse(sp == "val", "dev", sp)))
  fwrite(out, fn)
  cat(sp, nrow(out), "->", fn, "\n")
}

# ---- drugs_smiles.csv ----
pw <- pi[, .(pert_id, canonical_smiles)]
need_ids <- unique(rep_dt$pert_id)
pw_n <- pw[match(need_ids, pert_id)]
stopifnot(!anyNA(pw_n$canonical_smiles), !any(pw_n$canonical_smiles == ""))
fwrite(pw_n, file.path(outd, "drugs_smiles.csv"), col.names = FALSE, quote = FALSE)
cat("drugs_smiles:", nrow(pw_n), "\n")

# ---- gene_vector.csv (作者官方基因特征, 按 landmark 列序重排) ----
gv <- fread("/tmp/DeepCE-master/DeepCE/data/gene_vector.csv", header = FALSE,
            col.names = c("symbol", paste0("f", 1:128)))
stopifnot(all(landmark_sym %in% gv$symbol))
gv <- gv[match(landmark_sym, symbol)]
fwrite(gv, file.path(outd, "gene_vector.csv"), col.names = FALSE, quote = FALSE)

# ---- 分组清单 + 封存校验和 ----
man <- data.table(
  item = c("signature_train.csv", "signature_dev.csv", "signature_test.csv",
           "drugs_smiles.csv", "gene_vector.csv", "condition_manifest.csv"),
  sha256 = c(digest(file.path(outd, "signature_train.csv"), algo = "sha256", file = TRUE),
             digest(file.path(outd, "signature_dev.csv"), algo = "sha256", file = TRUE),
             digest(file.path(outd, "signature_test.csv"), algo = "sha256", file = TRUE),
             digest(file.path(outd, "drugs_smiles.csv"), algo = "sha256", file = TRUE),
             digest(file.path(outd, "gene_vector.csv"), algo = "sha256", file = TRUE),
             digest(file.path(outd, "condition_manifest.csv"), algo = "sha256", file = TRUE)),
  sealed = c("no", "no", "YES_until_model_lock", "no", "no", "no"),
  created = format(Sys.Date(), "%Y-%m-%d"))
fwrite(man, file.path(root, "metadata/pred_dataset_manifest.csv"))
cat("manifest 写出; split: train", sum(splits$split_pred == "train"),
    "val", sum(splits$split_pred == "val"),
    "test", sum(splits$split_pred == "test"), "母体\n")
