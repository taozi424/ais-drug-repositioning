# =====================================================================
# 10 构建药物条件矩阵 (dev + test 全部化合物; 第 8-10 周共用底表)
# 与 08_drug_ranking.R 完全相同的: QC 阈值、landmark->发现队列基因映射、
# 重复条件聚合(基因级均值)、普遍扰动方向 consensus。
# 输出 data/processed/drug_zmatrix.rds
#   Zc      952 基因 x 条件 (核心 7 细胞, 24h, 每个 parent x cell x dose 一条)
#   meta    条件元数据 (parent_key, cell_id, dose, n_sig_merged, iname, split, ...)
#   Zc_3h   3h 层 (核心细胞) 同结构, 仅用于时间层敏感性
#   consensus / genes / qc 统计
# =====================================================================
suppressPackageStartupMessages({
  library(data.table); library(yaml); library(rhdf5)
})
cfg <- yaml::read_yaml("signature/drug_ranking_config.yaml")
raw_dir <- "data/raw/LINCS/GSE70138"
t0 <- Sys.time()

gctx <- file.path(raw_dir, "GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx")
stopifnot(file.exists(gctx), file.size(gctx) > 1e9)
cells_core <- unlist(cfg$cells$core_panel)

# ---------- 1) 基因映射 (与 08 一致) ----------
lm_map <- fread("results/tables/L1000_landmark_to_discovery_genes.csv")
gi <- fread(cmd = sprintf("gzip -dc %s/GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz", raw_dir),
            sep = "\t")
gi_lm <- gi[pr_is_lm == 1]
gi_lm[, m_entrez := lm_map$gene[match(pr_gene_id, lm_map$entrez)]]
gi_lm[, m_symbol := lm_map$gene[match(pr_gene_symbol, lm_map$landmark_symbol)]]
gi_lm[, gene := fifelse(!is.na(m_entrez), m_entrez, m_symbol)]
gi_use <- gi_lm[!is.na(gene)]
gi_use <- gi_use[order(-as.integer(!is.na(m_entrez)))][!duplicated(gene)]

# ---------- 2) 样本元数据 / 化合物身份 / QC ----------
si <- fread(cmd = sprintf("gzip -dc %s/GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz", raw_dir),
            sep = "\t")
sm <- fread(cmd = sprintf("gzip -dc %s/GSE70138_Broad_LINCS_sig_metrics_2017-03-06.txt.gz", raw_dir),
            sep = "\t")
cp <- si[pert_type == "trt_cp"][sm[, .(sig_id, distil_cc_q75, distil_nsample)],
                                 on = "sig_id", nomatch = 0L,
                                 .(sig_id, pert_id, pert_iname, cell_id, pert_idose, pert_itime,
                                   distil_cc_q75, distil_nsample)]
cp <- cp[distil_cc_q75 >= cfg$qc$distil_cc_q75_min]
pr <- fread("metadata/compound_perturbagens_raw.csv")[pert_type == "trt_cp"]
parents <- fread("metadata/compound_identity_parents.csv")
cp[, parent_key := pr$parent_key[match(pert_id, pr$pert_id)]]
cp[, split := parents$split[match(parent_key, parents$parent_key)]]
cp <- cp[!is.na(parent_key)]
cp[, dose_um := suppressWarnings(as.numeric(sub(" um", "", pert_idose)))]
cat(sprintf("QC 合格 trt_cp %d 条 / 母体 %d / 未分组 %d\n", nrow(cp), uniqueN(cp$parent_key),
            sum(is.na(cp$split))))

# 条件键: parent x cell x dose x time (24h 核心细胞; 3h 核心细胞)
mk_cond <- function(d) {
  d <- d[cell_id %in% cells_core]
  d[, cond := paste(parent_key, cell_id, dose_um, pert_itime, sep = "___")]
  d
}
cp24 <- mk_cond(cp[pert_itime == cfg$qc$time_primary])
cp3 <- mk_cond(cp[pert_itime == "3 h"])
cat(sprintf("24h 核心细胞特征 %d 条 / 条件 %d; 3h %d 条 / 条件 %d\n",
            nrow(cp24), uniqueN(cp24$cond), nrow(cp3), uniqueN(cp3$cond)))

# ---------- 3) 读取 gctx 行子集(一次性), 再在内存中按条件聚合 ----------
# 注: 分块重复调用 h5read 会累积内存(rhdf5 未及时释放), 在 8GB 机器上被系统终止;
#     改为与 08 相同的一次性行子集读取, 随后用列索引聚合。
h5i <- h5ls(gctx)
drow <- h5i[h5i$name == "matrix" & h5i$group == "/0/DATA/0", ]
dims <- as.numeric(trimws(unlist(strsplit(drow$dim, "x"))))
stopifnot(dims[1] == 12328)
row_ids <- h5read(gctx, "0/META/ROW/id")
col_ids <- h5read(gctx, "0/META/COL/id")
gene_row_idx <- match(gi_use$pr_gene_id, row_ids)
stopifnot(!anyNA(gene_row_idx))
Zm_all <- h5read(gctx, "/0/DATA/0/matrix", index = list(gene_row_idx, seq_len(dims[2])))
h5closeAll()
rownames(Zm_all) <- gi_use$gene; colnames(Zm_all) <- col_ids
cat(sprintf("gctx 子集 %d x %d (%.0f MB)\n", nrow(Zm_all), ncol(Zm_all),
            as.numeric(object.size(Zm_all)) / 1e6))

aggregate_conditions <- function(d, Zall) {
  conds <- unique(d$cond)
  Zc <- matrix(NA_real_, nrow = nrow(Zall), ncol = length(conds),
               dimnames = list(rownames(Zall), conds))
  n_merged <- integer(length(conds)); names(n_merged) <- conds
  idx_of <- setNames(match(d$cond, conds), d$sig_id)
  pos <- match(colnames(Zall), names(idx_of))
  keep <- !is.na(pos)
  sub <- idx_of[pos[keep]]; cols <- which(keep)
  for (cid in unique(sub)) {
    jj <- cols[sub == cid]
    Zc[, cid] <- if (length(jj) == 1) Zall[, jj] else rowMeans(Zall[, jj, drop = FALSE])
    n_merged[cid] <- length(jj)
  }
  list(Zc = Zc, n_merged = n_merged)
}
cat("聚合 24h ...\n"); a24 <- aggregate_conditions(cp24, Zm_all)
cat("聚合 3h ...\n");  a3  <- aggregate_conditions(cp3, Zm_all)

# consensus: 全部 QC 合格 trt_cp 的基因级均值 (普遍扰动方向)
pos_all <- match(unique(cp$sig_id), col_ids); pos_all <- pos_all[!is.na(pos_all)]
consensus <- rowMeans(Zm_all[, pos_all, drop = FALSE]); names(consensus) <- gi_use$gene
rm(Zm_all); gc()

# ---------- 4) 条件元数据 ----------
meta_of <- function(d, a) {
  m <- d[, .(parent_key = parent_key[1], cell_id = cell_id[1], dose_um = dose_um[1],
             pert_itime = pert_itime[1], n_sig_merged = a$n_merged[cond[1]]),
         by = cond]
  m[, iname := parents$iname_std[match(parent_key, parents$parent_key)]]
  m[, pert_inames := parents$pert_inames[match(parent_key, parents$parent_key)]]
  m[, split := parents$split[match(parent_key, parents$parent_key)]]
  # ChEMBL 开发状态: 以大小写不敏感的名称键匹配(pert_iname 与 iname_std 写法不一致)
  st <- fread("metadata/compound_development_status.csv")
  st[, iname_key := tolower(trimws(pert_iname))]
  st <- st[!duplicated(iname_key)]
  m[, iname_key := tolower(trimws(iname))]
  m[, max_phase := as.numeric(st$max_phase[match(iname_key, st$iname_key)])]
  m[, first_approval := st$first_approval[match(iname_key, st$iname_key)]]
  m[, matched_by := st$matched_by[match(iname_key, st$iname_key)]]
  m[, iname_key := NULL]
  m[cond %in% colnames(a$Zc)]
}
meta24 <- meta_of(cp24, a24)
meta3 <- meta_of(cp3, a3)

# 剂量选择 (冻结规则: 10uM 优先, 否则 24h 内 log 距离最近): 每个 parent x cell 一条
sel_dose <- meta24[, {
  d <- dose_um
  pick <- if (10 %in% d) 10 else d[which.min(abs(log10(d) - 1))]
  .(cond = cond[match(pick, dose_um)], dose_used = pick)
}, by = .(parent_key, cell_id)]
meta24[, primary_dose := cond %in% sel_dose$cond]

saveRDS(list(Zc = a24$Zc, meta = meta24, Zc_3h = a3$Zc, meta_3h = meta3,
             consensus = consensus, genes = gi_use$gene,
             n_landmark_mapped = nrow(gi_use), cells_core = cells_core,
             config_frozen = "signature/drug_ranking_config.yaml",
             built_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
        "data/processed/drug_zmatrix.rds")
cat(sprintf("完成: %d 基因 x %d 条件 (24h) + %d 条件 (3h); 用时 %.1f min\n",
            nrow(a24$Zc), ncol(a24$Zc), ncol(a3$Zc),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
print(meta24[, .(n_cond = .N, n_parent = uniqueN(parent_key)), by = cell_id][order(-n_cond)])
