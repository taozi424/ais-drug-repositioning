# =====================================================================
# 08 开发集药物排名 (GSE70138 Level5, 规则见 drug_ranking_config.yaml)
# 测试化合物(split=="test")不计算不输出, 至 AI 最终评估解封
# =====================================================================
suppressPackageStartupMessages({
  library(data.table); library(yaml); library(rhdf5); library(matrixStats)
})
cfg <- yaml::read_yaml("signature/drug_ranking_config.yaml")
raw_dir <- "data/raw/LINCS/GSE70138"
tab_dir <- "results/tables"
set.seed(20261109)

gctx_gz <- file.path(raw_dir, "GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz")
gctx <- sub("\\.gz$", "", gctx_gz)
if (!file.exists(gctx)) {
  cat("解压 gctx ...\n")
  stopifnot(system(sprintf("gunzip -k %s", gctx_gz)) == 0)
}
stopifnot(file.exists(gctx), file.size(gctx) > 1e9)

# ---------- 1) AIS 侧输入 ----------
allg <- fread("signature/AIS_all_genes.csv")
lm_map <- fread("results/tables/L1000_landmark_to_discovery_genes.csv")   # landmark_symbol -> gene
sig_genes <- allg[in_signature == TRUE]
cat(sprintf("签名基因 %d (上 %d / 下 %d)\n", nrow(sig_genes),
            sum(sig_genes$direction == "Up"), sum(sig_genes$direction == "Down")))

# 随机签名匹配源: GSE58294 (全血 log 尺度, 提供均值代理) + GSE16561 sd (变异)
# 注: GSE16561 为基线变换矩阵, 均值无信息, 故均值取自 GSE58294
e58 <- readRDS("data/processed/GSE58294_expr_gene.rds")
e61 <- readRDS("data/processed/GSE16561_expr_gene.rds")
pool_genes <- lm_map$gene                       # 952 个发现队列映射基因
stopifnot(all(pool_genes %in% rownames(e61)))
# 表达均值代理: GSE58294 覆盖 862/952; 缺失基因(平台无该探针)退化为仅按 sd 分箱
has_mean <- pool_genes %in% rownames(e58)
g_mean <- setNames(rep(NA_real_, length(pool_genes)), pool_genes)
g_mean[has_mean] <- rowMeans(e58[pool_genes[has_mean], ])
g_sd   <- apply(e61[pool_genes, ], 1, sd)
tvals  <- allg[match(pool_genes, gene), t]; names(tvals) <- pool_genes
tvals_sig <- sig_genes$t; names(tvals_sig) <- sig_genes$gene

qcut <- function(x) quantile(x, probs = seq(.1, .9, .1), na.rm = TRUE)
sd_bin  <- setNames(findInterval(g_sd, qcut(g_sd), rightmost.closed = TRUE) + 1, pool_genes)
mean_bin <- setNames(rep(NA_integer_, length(pool_genes)), pool_genes)
mean_bin[has_mean] <- findInterval(g_mean[has_mean], qcut(g_mean[has_mean]), rightmost.closed = TRUE) + 1
bin_key <- ifelse(is.na(mean_bin), NA_character_,
                  paste0(mean_bin, "_", sd_bin))
bin_members <- split(pool_genes[!is.na(bin_key)], bin_key[!is.na(bin_key)])
sd_members  <- split(pool_genes, sd_bin)
random_t_vec <- function() {
  vapply(names(tvals_sig), function(g) {
    cand <- if (!is.na(bin_key[g])) bin_members[[bin_key[g]]] else sd_members[[sd_bin[g]]]
    if (is.null(cand) || length(cand) == 0) cand <- pool_genes
    unname(tvals[sample(cand, 1)])
  }, numeric(1))
}

# ---------- 2) 样本元数据 + 化合物身份 + QC ----------
si <- fread(cmd = sprintf("gzip -dc %s/GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz", raw_dir))
sm <- fread(cmd = sprintf("gzip -dc %s/GSE70138_Broad_LINCS_sig_metrics_2017-03-06.txt.gz", raw_dir))
cp_meta <- si[pert_type == "trt_cp"]
cp_meta <- cp_meta[sm[, .(sig_id, distil_cc_q75, distil_nsample)],
                   on = "sig_id", nomatch = 0L,
                   .(sig_id, pert_id, pert_iname, cell_id, pert_idose, pert_itime,
                     distil_cc_q75, distil_nsample)]
cp_meta <- cp_meta[distil_cc_q75 >= cfg$qc$distil_cc_q75_min]
pr <- fread("metadata/compound_perturbagens_raw.csv")[pert_type == "trt_cp"]
parents <- fread("metadata/compound_identity_parents.csv")
cp_meta[, parent_key := pr$parent_key[match(pert_id, pr$pert_id)]]
cp_meta[, split := parents$split[match(parent_key, parents$parent_key)]]
cat(sprintf("QC 合格 trt_cp %d 条\n", nrow(cp_meta)))

dev <- cp_meta[split == "dev"]
dev[, dose_um := as.numeric(sub(" um", "", pert_idose))]
dev24 <- dev[pert_itime == cfg$qc$time_primary]
dose_opts <- sort(unique(dev24$dose_um))
dev24[, dose_used := { d <- unique(dose_um)
                       if (10 %in% d) 10 else d[which.min(abs(log10(d) - 1))] },
      by = .(parent_key, cell_id)]

# ---------- 3) 读取 gctx landmark 行 ----------
h5i <- h5ls(gctx)
drow <- h5i[h5i$name == "matrix" & h5i$group == "/0/DATA/0", ]
stopifnot(nrow(drow) == 1)
dpath <- "/0/DATA/0/matrix"
dims <- as.numeric(trimws(unlist(strsplit(drow$dim, "x"))))
stopifnot(length(dims) == 2, dims[1] == 12328)
cat("gctx 矩阵维度:", dims, "\n")
row_ids <- h5read(gctx, "0/META/ROW/id")
col_ids <- h5read(gctx, "0/META/COL/id")
gi <- fread(cmd = sprintf("gzip -dc %s/GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz", raw_dir))
gi_lm <- gi[pr_is_lm == 1]                                   # 978 landmark
# 映射 GSE70138 landmark -> 发现队列 gene:
# 优先 Entrez(pr_gene_id), 辅以官方符号; lm_map 含阶段二的 Entrez(ALIAS) 解析
lm_map <- fread("results/tables/L1000_landmark_to_discovery_genes.csv")
gi_lm[, m_entrez := lm_map$gene[match(pr_gene_id, lm_map$entrez)]]
gi_lm[, m_symbol := lm_map$gene[match(pr_gene_symbol, lm_map$landmark_symbol)]]
gi_lm[, gene := fifelse(!is.na(m_entrez), m_entrez, m_symbol)]
gi_use <- gi_lm[!is.na(gene)]
# 同一基因符号对应多个 L1000 探针 ID (RPN1/LAP3): 优先保留 Entrez 匹配行
gi_use <- gi_use[order(-as.integer(!is.na(m_entrez)))][!duplicated(gene)]
cat(sprintf("landmark 978 -> 与阶段二 952 映射空间并集命中 %d (Entrez %d / 符号 %d)\n",
            nrow(gi_use), sum(!is.na(gi_lm$m_entrez)), sum(!is.na(gi_lm$m_symbol))))
gene_row_idx <- match(gi_use$pr_gene_id, row_ids)
Z <- h5read(gctx, dpath, index = list(gene_row_idx, seq_len(dims[2])))
h5closeAll()
colnames(Z) <- col_ids; rownames(Z) <- gi_use$gene
cat(sprintf("提取 Z: %d 基因 x %d 特征\n", nrow(Z), ncol(Z)))

# ---------- 4) 普遍扰动方向 (全部合格 trt_cp, 不分开发/测试) ----------
pos_qc <- match(cp_meta$sig_id, col_ids)
stopifnot(!anyNA(pos_qc))
consensus <- rowMeans(Z[, pos_qc, drop = FALSE], na.rm = TRUE)

# ---------- 5) 条件级聚合与指标 (仅开发组) ----------
pos_dev <- match(dev24$sig_id, col_ids)
stopifnot(!anyNA(pos_dev))
cond_key <- paste(dev24$parent_key, dev24$cell_id, dev24$dose_used, sep = "___")
uiq <- unique(cond_key)
grp <- split(seq_len(nrow(dev24)), cond_key)
Zc <- matrix(NA_real_, nrow = nrow(Z), ncol = length(uiq), dimnames = list(rownames(Z), uiq))
for (k in seq_along(grp)) Zc[, k] <- rowMeans(Z[, pos_dev[grp[[k]]], drop = FALSE], na.rm = TRUE)
Zc_sig <- Zc[sig_genes$gene, , drop = FALSE]                 # 185 x 条件
t_sig <- tvals_sig[sig_genes$gene]

Rs     <- apply(Zc_sig, 2, function(x) cor(rank(x), rank(t_sig)))
R_sens <- apply(Zc,     2, function(x) cor(rank(x), rank(t_all <- tvals[rownames(Zc)])))
med_up <- colMedians(Zc_sig[sig_genes$direction == "Up",   , drop = FALSE])
med_dn <- colMedians(Zc_sig[sig_genes$direction == "Down", , drop = FALSE])
rev_up <- colMeans(Zc_sig[sig_genes$direction == "Up",   , drop = FALSE] < 0)
rev_dn <- colMeans(Zc_sig[sig_genes$direction == "Down", , drop = FALSE] > 0)
r_cons  <- apply(Zc, 2, function(x) cor(x, consensus))
medabsz <- colMedians(abs(Zc))

parts <- tstrsplit(uiq, "___", fixed = TRUE)
res <- data.table(parent_key = parts[[1]], cell_id = parts[[2]],
                  dose_used = as.numeric(parts[[3]]))
res[, `:=`(n_sig_merged = lengths(grp)[uiq],
           R_primary = -Rs[uiq],
           R_sens_all_lm = -R_sens[uiq],
           median_z_up = med_up[uiq],
           median_z_down = med_dn[uiq],
           reversal_up_pct = round(100 * rev_up[uiq], 1),
           reversal_dn_pct = round(100 * rev_dn[uiq], 1),
           r_consensus = r_cons[uiq],
           median_absz = medabsz[uiq])]
res[, flag_stress := r_cons > 0.5]
res[, flag_strong_effect := median_absz > 3]
res[, candidate := R_primary > 0 & median_z_up < 0 & median_z_down > 0]
# 高暴露标记: 该条件实际剂量为 20uM
res[dev24[, .(dose_used = unique(dose_used)), by = .(parent_key, cell_id)],
    high_exposure := i.dose_used >= 20, on = c("parent_key", "cell_id")]

res[, iname := parents$iname_std[match(parent_key, parents$parent_key)]]
res[, pert_inames := parents$pert_inames[match(parent_key, parents$parent_key)]]
st <- fread("metadata/compound_development_status.csv")
res[, max_phase := st$max_phase[match(iname, st$pert_iname)]]
res[, first_approval := st$first_approval[match(iname, st$pert_iname)]]
res[, matched_by := st$matched_by[match(iname, st$pert_iname)]]
setorder(res, parent_key, cell_id, dose_used)
cat(sprintf("开发集条件明细 %d 条; 候选规则命中 %d 条\n", nrow(res), sum(res$candidate)))

# ---------- 6) 随机签名 null (各核心细胞; R 尺度 = -Spearman) ----------
cells <- unlist(cfg$cells$core_panel)
null_rows <- list()
for (cl in cells) {
  cl_of <- parts[[2]]
  cond_cl <- uiq[cl_of == cl]
  if (!length(cond_cl)) next
  Zs <- scale(Zc_sig[, cond_cl, drop = FALSE], center = TRUE, scale = TRUE)
  Zs <- pmin(pmax(Zs, -10), 10)
  nullR <- NULL
  for (b in 1:100) {
    Tm <- replicate(10, random_t_vec())
    Tr <- scale(apply(Tm, 2, rank), center = TRUE, scale = TRUE)
    nullR <- c(nullR, -as.numeric(crossprod(Zs, Tr)) / nrow(Zs))
  }
  null_rows[[cl]] <- data.table(cell_id = cl, n_cond = length(cond_cl),
                                null_R_p50 = quantile(nullR, .50, names = FALSE),
                                null_R_p95 = quantile(nullR, .95, names = FALSE),
                                null_R_p99 = quantile(nullR, .99, names = FALSE))
}
null_tab <- rbindlist(null_rows)
fwrite(null_tab, file.path(tab_dir, "drug_ranking_null_reference.csv"))
print(null_tab)
res[, null_R_p95 := null_tab$null_R_p95[match(cell_id, null_tab$cell_id)]]
res[, exceeds_null95 := R_primary > null_R_p95]
fwrite(res, file.path(tab_dir, "drug_condition_detail_dev.csv"))
cat("明细表(含 null 标记)已写出\n")

# ---------- 7) 固定条件排名 (分析 A: 各核心细胞) ----------
for (cl in cells) {
  d <- res[cell_id == cl]
  if (!nrow(d)) next
  setorder(d, -R_primary)
  fwrite(d, file.path(tab_dir, sprintf("drug_ranking_dev_%s_24h.csv", cl)))
}

# ---------- 8) 共同面板汇总 (分析 B: 7 核心细胞齐备) ----------
cov <- res[, .(n_cells = uniqueN(cell_id)), by = parent_key]
full_parents <- cov[n_cells == length(cells)]$parent_key
common <- res[parent_key %in% full_parents]
summ <- common[, .(R_median = median(R_primary),
                   R_min = min(R_primary), R_max = max(R_primary),
                   n_cells = uniqueN(cell_id),
                   cell_consistent_pct = round(100 * mean(candidate)),
                   median_z_up = median(median_z_up),
                   median_z_down = median(median_z_down),
                   r_consensus_med = median(r_consensus),
                   max_phase = max_phase[1], first_approval = first_approval[1]),
               by = .(parent_key, iname)]
summ[, candidate_common := R_median > 0 & median_z_up < 0 & median_z_down > 0]
setorder(summ, -R_median)
fwrite(summ, file.path(tab_dir, "drug_ranking_dev_common_panel.csv"))
cat(sprintf("共同面板化合物 %d 个; 候选 %d 个\n", nrow(summ), sum(summ$candidate_common)))

# ---------- 9) Top20 展示表 (不含疗效解读) ----------
top20_A <- rbindlist(lapply(cells, function(cl) {
  d <- res[cell_id == cl][order(-R_primary)][1:min(20, .N)]
  d[, analysis := "A_fixed_cell"]
}), use.names = TRUE)
top20_B <- summ[1:min(20, .N), .(parent_key, iname, R_median, R_min, R_max,
                                 median_z_up, median_z_down, max_phase,
                                 candidate_common)]
top20_B[, analysis := "B_common_panel_median"]
fwrite(top20_A, file.path(tab_dir, "top20_display_analysis_A.csv"))
fwrite(top20_B, file.path(tab_dir, "top20_display_analysis_B.csv"))
saveRDS(list(consensus = consensus, null_tab = null_tab), "data/processed/drug_ranking_aux.rds")
cat("完成: 08_drug_ranking.R\n")
