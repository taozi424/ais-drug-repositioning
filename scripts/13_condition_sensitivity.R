# =====================================================================
# 13 第 9 周: 条件与组成敏感性
# A 签名基因数 (50/100/150)  B 细胞背景  C 剂量/时间层  D 单条件反向  E 组成调整
# 所有变体都只改一个维度, 其余保持冻结; 结果只解释为敏感性, 不用于挑选最有利结果
# =====================================================================
suppressPackageStartupMessages({
  library(data.table); library(matrixStats); library(yaml)
})
tab <- "results/tables"
cfg <- yaml::read_yaml("signature/drug_ranking_config.yaml")
allg <- fread("signature/AIS_all_genes.csv")
sig <- fread("signature/AIS_screen_signature.csv")
Z <- readRDS("data/processed/drug_zmatrix.rds")
meta <- copy(Z$meta)
m <- meta[primary_dose == TRUE]
Zc <- Z$Zc[, m$cond, drop = FALSE]
cells <- unlist(cfg$cells$core_panel)

cor_rank <- function(Zr, rows, tvec) {
  A <- Zr[rows, , drop = FALSE]
  A <- sweep(A, 2, colMeans(A), "-")
  b <- rank(tvec); b <- b - mean(b)
  as.numeric(colSums(A * b)) / sqrt(colSums(A^2) * sum(b^2))
}
score_set <- function(genes, tvals) {
  g <- genes[genes %in% rownames(Zc)]
  Zs <- Zc[g, , drop = FALSE]
  Zr <- apply(Zs, 2, rank)
  -cor_rank(Zr, g, tvals[g])
}
agg_parent <- function(R) {
  dt <- data.table(parent_key = m$parent_key, R = R)
  dt[, .(R_med = median(R), n_cells = .N), by = parent_key]
}
final <- fread(file.path(tab, "final_drug_ranking_all.csv"))
base <- final[, .(parent_key, R_base = R_median, rank_base = rank_overall, candidate_base = candidate_common)]

# ---------- A 签名基因数 ----------
fdr_ok <- allg[FDR < 0.05]
resA <- list()
for (k in c(50, 100, 150)) {
  up <- head(fdr_ok[t > 0][order(-abs(t))], k)
  dn <- head(fdr_ok[t < 0][order(-abs(t))], k)
  g <- c(up$gene, dn$gene); tv <- setNames(c(up$t, dn$t), g)
  R <- score_set(g, tv); a <- agg_parent(R)
  x <- merge(base, a, by = "parent_key")
  x[, cand := R_med > 0]
  resA[[as.character(k)]] <- data.table(
    variant = sprintf("top%d_per_direction", k), n_up = nrow(up), n_dn = nrow(dn),
    rho_vs_base = cor(x$R_base, x$R_med, method = "spearman"),
    top20_overlap = length(intersect(x[order(-R_base)][1:20]$parent_key, x[order(-R_med)][1:20]$parent_key)),
    overlap_Rpositive_with_184 = length(intersect(base[candidate_base == TRUE]$parent_key,
                                                  x[R_med > 0]$parent_key)),
    n_R_positive = sum(x$R_med > 0))
  fwrite(x[, .(parent_key, R_base, R_variant = R_med)], file.path(tab, sprintf("sens_genes_%d.csv", k)))
}
resA <- rbindlist(resA)
fwrite(resA, file.path(tab, "sensitivity_signature_size.csv")); print(resA)

# ---------- B 细胞背景 ----------
cond <- fread(file.path(tab, "final_drug_condition_detail.csv"))
wide <- dcast(cond, parent_key ~ cell_id, value.var = "R_primary")
pw <- list()
for (i in 1:(length(cells) - 1)) for (j in (i + 1):length(cells)) {
  a <- cells[i]; b <- cells[j]
  ok <- !is.na(wide[[a]]) & !is.na(wide[[b]])
  pw[[length(pw) + 1]] <- data.table(cell_a = a, cell_b = b, n = sum(ok),
                                     rho = cor(wide[[a]][ok], wide[[b]][ok], method = "spearman"))
}
pw <- rbindlist(pw)
fwrite(pw, file.path(tab, "sensitivity_cell_pairwise.csv"))
cat(sprintf("细胞间 R 相关: 中位 %.3f (范围 %.3f - %.3f)\n", median(pw$rho), min(pw$rho), max(pw$rho)))

ov <- list()
for (cl in cells) {
  d <- cond[cell_id == cl][order(-R_primary)]
  top <- d[1:20]$parent_key
  ov[[cl]] <- data.table(cell_id = cl, n_compounds = nrow(d),
                         top20_in_final_candidates = sum(top %in% base[candidate_base == TRUE]$parent_key),
                         top20_in_final_top20 = sum(top %in% base[order(-R_base)][1:20]$parent_key))
}
ov <- rbindlist(ov)
fwrite(ov, file.path(tab, "sensitivity_cell_top20.csv")); print(ov)

# ---------- C 剂量层 ----------
alldose <- meta[cell_id %in% cells]
Rall <- NULL
{
  Zl <- Z$Zc[, alldose$cond, drop = FALSE]
  g <- sig$gene; tv <- setNames(allg$t[match(g, allg$gene)], g)
  Rl <- score_set(g, tv)   # 注意 score_set 使用 Zc, 这里手工算
  Zs <- Zl[g, , drop = FALSE]; Zr <- apply(Zs, 2, rank)
  Rl <- -cor_rank(Zr, g, tv)
  dtl <- data.table(parent_key = alldose$parent_key, R = Rl)
  per <- dtl[, .(R_dose_med = median(R)), by = parent_key]
  x <- merge(final[, .(parent_key, R_base = R_median, rank_base = rank_overall)], per, by = "parent_key")
  cat(sprintf("剂量层: 主剂量 vs 全剂量中位数 R 相关 %.3f\n", cor(x$R_base, x$R_dose_med, method = "spearman")))
  fwrite(x, file.path(tab, "sensitivity_dose_all_doses.csv"))
  fwrite(data.table(variant = "all_doses_median_24h", n_conditions = ncol(Zl),
                    rho_vs_primary = cor(x$R_base, x$R_dose_med, method = "spearman")),
         file.path(tab, "sensitivity_dose_summary.csv"))
}
# 3h 覆盖
cat(sprintf("3h 条件 %d 条 / 母体 %d (覆盖不足, 仅报告, 不做排名)\n", ncol(Z$Zc_3h), uniqueN(Z$meta_3h$parent_key)))
fwrite(data.table(n_conditions_3h = ncol(Z$Zc_3h), n_parents_3h = uniqueN(Z$meta_3h$parent_key),
                  note = "3h 层条件过少, 无法支撑剂量-时间层比较; 主分析仅 24h"),
       file.path(tab, "sensitivity_time_layer.csv"))

# ---------- D 单条件反向 ----------
sc <- cond[, .(n_cells = uniqueN(cell_id), n_cand = sum(candidate),
               n_pos = sum(R_primary > 0), R_min = min(R_primary), R_max = max(R_primary)),
           by = parent_key]
sc <- merge(sc, final[, .(parent_key, iname, R_median, candidate_common)], by = "parent_key")
sc[, single_condition_only := n_cells >= 4 & n_cand == 1]
sc[, opposite_in_one_cell := n_cells >= 4 & n_pos == 1 & n_cand == 0]
sc[, flag := fifelse(single_condition_only, "candidate_in_one_cell_only",
              fifelse(opposite_in_one_cell, "positive_R_in_one_cell_only", "none"))]
fwrite(sc[flag != "none"][order(-R_median)], file.path(tab, "single_condition_reversal.csv"))
cat(sprintf("单条件候选 %d 个; 仅一个细胞 R>0 的 %d 个\n", sum(sc$single_condition_only),
            sum(sc$opposite_in_one_cell)))

# ---------- E 组成调整 ----------
adj <- fread(file.path(tab, "AIS_all_genes_composition_adjusted.csv"))
adj_ok <- adj[FDR < 0.05]
up <- head(adj_ok[t > 0][order(-abs(t))], 100); dn <- head(adj_ok[t < 0][order(-abs(t))], 100)
g <- c(up$gene, dn$gene); tv <- setNames(c(up$t, dn$t), g)
R <- score_set(g, tv); a <- agg_parent(R)
x <- merge(base, a, by = "parent_key")
x[, cand_adj := R_med > 0]
cat(sprintf("组成调整: 签名 %d 基因(上 %d/下 %d); 与主签名重叠 %d; R 相关 %.3f\n", length(g), nrow(up), nrow(dn),
            length(intersect(g, sig$gene)), cor(x$R_base, x$R_med, method = "spearman")))
summE <- data.table(variant = "composition_adjusted", n_genes = length(g),
                    n_up = nrow(up), n_dn = nrow(dn),
                    gene_overlap_with_primary = length(intersect(g, sig$gene)),
                    rho_vs_base = cor(x$R_base, x$R_med, method = "spearman"),
                    top20_overlap = length(intersect(x[order(-R_base)][1:20]$parent_key, x[order(-R_med)][1:20]$parent_key)),
                    n_candidate_adj = sum(x$cand_adj),
                    candidates_surviving = sum(x$parent_key %in% base[candidate_base == TRUE]$parent_key & x$cand_adj))
fwrite(summE, file.path(tab, "sensitivity_composition_adjusted.csv"))
fwrite(x[, .(parent_key, R_base, R_adjusted = R_med)], file.path(tab, "sens_composition_all.csv"))
print(summE)
cat("完成: 13_condition_sensitivity.R\n")
