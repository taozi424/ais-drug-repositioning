# =====================================================================
# 11 第 8 周: 最终药物排名 (dev + test 全部有实测特征的化合物)
# 规则: signature/drug_ranking_config.yaml (冻结于 2026-11-09 15:20)
# 与第 4 周脚本的差异(版本修正, 不改规则只改实现):
#   08 的 dose_used 在 parent x cell 组内为常数, 导致同一 parent-cell 的不同剂量
#   被合并进同一条件; 冻结配置要求"同一 (parent, cell, dose, time) 内聚合"。
#   本脚本按 (parent, cell, dose, time) 聚合后再取主剂量, 与配置一致。
# =====================================================================
suppressPackageStartupMessages({
  library(data.table); library(matrixStats); library(yaml)
})
cfg <- yaml::read_yaml("signature/drug_ranking_config.yaml")
tab <- "results/tables"
Z <- readRDS("data/processed/drug_zmatrix.rds")
Zc_all <- Z$Zc; meta <- copy(Z$meta); consensus <- Z$consensus
sig <- fread("signature/AIS_screen_signature.csv")
allg <- fread("signature/AIS_all_genes.csv")
cells <- unlist(cfg$cells$core_panel)
set.seed(20261109)

# ---------- 1) 主剂量条件 ----------
m <- meta[primary_dose == TRUE]
Zc <- Zc_all[, m$cond, drop = FALSE]
m[, condition_id := cond][, cond := NULL][, dose_used := dose_um]
cat(sprintf("主剂量条件 %d 条 / 母体 %d (dev %d / test %d)\n", nrow(m), uniqueN(m$parent_key),
            uniqueN(m[split == "dev"]$parent_key), uniqueN(m[split == "test"]$parent_key)))

t_all <- allg$t[match(rownames(Zc), allg$gene)]; names(t_all) <- rownames(Zc)
t_sig <- allg$t[match(sig$gene, allg$gene)]; names(t_sig) <- sig$gene
up <- sig$direction == "Up"; dn <- sig$direction == "Down"

# ---------- 2) 秩相关: 中心必须用被比较子集自身的均值, 不能用全局秩均值 ----------
cor_rank <- function(Zr, rows, tvec) {
  A <- Zr[rows, , drop = FALSE]
  A <- sweep(A, 2, colMeans(A), "-")
  b <- rank(tvec); b <- b - mean(b)
  as.numeric(colSums(A * b)) / sqrt(colSums(A^2) * sum(b^2))
}
Zr <- apply(Zc, 2, rank)
R_primary <- -cor_rank(Zr, names(t_sig), t_sig)
R_sens <- -cor_rank(Zr, rownames(Zc), t_all)

m[, `:=`(R_primary = R_primary, R_sens_all_lm = R_sens,
         median_z_up = colMedians(Zc[sig$gene[up], , drop = FALSE]),
         median_z_down = colMedians(Zc[sig$gene[dn], , drop = FALSE]),
         reversal_up_pct = round(100 * colMeans(Zc[sig$gene[up], , drop = FALSE] < 0), 1),
         reversal_dn_pct = round(100 * colMeans(Zc[sig$gene[dn], , drop = FALSE] > 0), 1),
         r_consensus = as.numeric(cor(Zc, consensus)),
         median_absz = colMedians(abs(Zc)))]
m[, `:=`(flag_stress = r_consensus > 0.5, flag_strong_effect = median_absz > 3,
         high_exposure = dose_used >= 20)]
m[, candidate := R_primary > 0 & median_z_up < 0 & median_z_down > 0]

# ---------- 3) 随机签名 null (配置: 952 映射 landmark 内按表达均值/变异分箱有放回抽) ----------
e58 <- readRDS("data/processed/GSE58294_expr_gene.rds")
e61 <- readRDS("data/processed/GSE16561_expr_gene.rds")
pool <- rownames(Zc)
has_mean <- pool %in% rownames(e58)
g_mean <- setNames(rep(NA_real_, length(pool)), pool)
g_mean[has_mean] <- rowMeans(e58[pool[has_mean], , drop = FALSE])
g_sd <- setNames(apply(e61[pool, , drop = FALSE], 1, sd), pool)
qcut <- function(x) quantile(x, probs = seq(.1, .9, .1), na.rm = TRUE)
sd_bin <- setNames(findInterval(g_sd, qcut(g_sd), rightmost.closed = TRUE) + 1, pool)
mean_bin <- setNames(rep(NA_integer_, length(pool)), pool)
mean_bin[has_mean] <- findInterval(g_mean[has_mean], qcut(g_mean[has_mean]), rightmost.closed = TRUE) + 1
bin_key <- ifelse(is.na(mean_bin), NA_character_, paste0(mean_bin, "_", sd_bin))
bin_members <- split(pool[!is.na(bin_key)], bin_key[!is.na(bin_key)])
sd_members <- split(pool, sd_bin)
random_t <- function() {
  vapply(names(t_sig), function(g) {
    cand <- if (!is.na(bin_key[g])) bin_members[[bin_key[g]]] else sd_members[[sd_bin[g]]]
    if (is.null(cand) || !length(cand)) cand <- pool
    unname(t_all[sample(cand, 1)])
  }, numeric(1))
}
null_rows <- list()
for (cl in cells) {
  cols <- which(m$cell_id == cl)
  if (!length(cols)) next
  Zs <- Zr[names(t_sig), cols, drop = FALSE]
  Zs <- sweep(Zs, 2, colMeans(Zs), "-")
  Zs <- sweep(Zs, 2, sqrt(colSums(Zs^2) / nrow(Zs)), "/")
  Zs <- pmin(pmax(Zs, -10), 10)
  nullR <- numeric(0)
  for (b in 1:100) {
    Tm <- replicate(10, random_t())
    for (k in 1:10) {
      bb <- rank(Tm[, k]); bb <- bb - mean(bb); bb <- bb / sqrt(sum(bb^2))
      nullR <- c(nullR, -as.numeric(crossprod(Zs, bb)) / sqrt(nrow(Zs) - 1))
    }
  }
  null_rows[[cl]] <- data.table(cell_id = cl, n_cond = length(cols),
                                null_R_p50 = quantile(nullR, .50, names = FALSE),
                                null_R_p95 = quantile(nullR, .95, names = FALSE),
                                null_R_p99 = quantile(nullR, .99, names = FALSE))
}
null_tab <- rbindlist(null_rows)
fwrite(null_tab, file.path(tab, "final_null_reference.csv"))
print(null_tab)
m[, null_R_p95 := null_tab$null_R_p95[match(cell_id, null_tab$cell_id)]]
m[, exceeds_null95 := R_primary > null_R_p95]
setorder(m, parent_key, cell_id)
fwrite(m, file.path(tab, "final_drug_condition_detail.csv"))
cat(sprintf("条件 %d; 候选条件 %d; 超 null95 %d (%.1f%%)\n", nrow(m), sum(m$candidate),
            sum(m$exceeds_null95), 100 * mean(m$exceeds_null95)))

# ---------- 4) 固定细胞排名 ----------
for (cl in cells) {
  d <- m[cell_id == cl]
  if (!nrow(d)) next
  setorder(d, -R_primary)
  d[, rank_in_cell := .I]
  fwrite(d, file.path(tab, sprintf("final_ranking_%s_24h.csv", cl)))
}

# ---------- 5) 汇总 + 证据分级 ----------
summ <- m[, .(R_median = median(R_primary), R_min = min(R_primary), R_max = max(R_primary),
              n_cells = uniqueN(cell_id), cell_consistent_pct = round(100 * mean(candidate)),
              n_candidate_cells = sum(candidate),
              median_z_up = median(median_z_up), median_z_down = median(median_z_down),
              r_consensus_med = median(r_consensus),
              any_stress = any(flag_stress), any_strong_effect = any(flag_strong_effect),
              max_absz = max(median_absz), doses = paste(sort(unique(dose_used)), collapse = "/"),
              iname = iname[1], max_phase = max_phase[1], first_approval = first_approval[1],
              split = split[1]),
          by = parent_key]
summ[, candidate_common := R_median > 0 & median_z_up < 0 & median_z_down > 0]
summ[, evidence_tier := fifelse(n_cells == length(cells), "A_common_panel",
                        fifelse(n_cells >= 4, "B_partial_panel", "C_single_or_few"))]
summ[, exceeds_null := R_median > max(null_tab$null_R_p95)]
setorder(summ, -R_median)
summ[, rank_overall := .I]
summ[, rank_A_only := NA_integer_]
summ[evidence_tier == "A_common_panel", rank_A_only := seq_len(.N)]
fwrite(summ, file.path(tab, "final_drug_ranking_all.csv"))
cat(sprintf("化合物 %d (A %d / B %d / C %d); 共同面板候选 %d\n", nrow(summ),
            sum(summ$evidence_tier == "A_common_panel"), sum(summ$evidence_tier == "B_partial_panel"),
            sum(summ$evidence_tier == "C_single_or_few"), sum(summ$candidate_common)))

# ---------- 6) 实测候选 ----------
cand <- summ[candidate_common == TRUE]
fwrite(cand, file.path(tab, "final_candidates_measured.csv"))
cat(sprintf("实测候选 %d; A级 %d; 上市(phase4) %d; 临床期(1-3) %d; 研究阶段 %d\n", nrow(cand),
            sum(cand$evidence_tier == "A_common_panel"), sum(cand$max_phase == 4, na.rm = TRUE),
            sum(cand$max_phase %in% 1:3, na.rm = TRUE),
            sum(is.na(cand$max_phase) | cand$max_phase < 1, na.rm = TRUE)))

# ---------- 7) 纯预测分支 ----------
gate <- jsonlite::fromJSON("ai/output/test_evaluation.json")
status <- data.table(branch = "predicted_only_supplement",
                     gate_rule = gate$gate$criterion, gate_passed = gate$gate$passed,
                     n_candidates = 0L, reason = gate$gate$consequence,
                     note = "未对缺少实验特征的化合物补充预测; 实测候选与纯预测候选分开报告, 不因分数相近视为证据等价")
fwrite(status, file.path(tab, "prediction_branch_status.csv"))
fwrite(data.table(parent_key = character(0), iname = character(0), R_predicted = numeric(0),
                  applicability_domain = character(0), max_tanimoto_train = numeric(0),
                  uncertainty = character(0)),
       file.path(tab, "final_candidates_predicted.csv"))

st <- fread("metadata/compound_development_status.csv")
st[, iname_key := tolower(trimws(pert_iname))]
miss <- st[!iname_key %in% tolower(trimws(summ$iname))]
fwrite(miss[, .(pert_id, pert_iname, max_phase, first_approval, matched_by)],
       file.path(tab, "no_lincs_feature_drugs.csv"))
cat(sprintf("无 L1000 特征的 ChEMBL 记录(仅统计, 不产生候选) %d\n", nrow(miss)))

# ---------- 8) 与第 4 周开发集排名的对照 ----------
old <- fread(file.path(tab, "drug_ranking_dev_common_panel.csv"))
cmp <- merge(old[, .(parent_key, R_median_wk4 = R_median, iname_wk4 = iname)],
             summ[, .(parent_key, R_median_final = R_median, rank_final = rank_overall,
                      evidence_tier)], by = "parent_key")
cmp[, rank_wk4 := frank(-R_median_wk4, ties.method = "first")]
fwrite(cmp[order(rank_wk4)], file.path(tab, "ranking_wk4_vs_final.csv"))
cat(sprintf("共同面板 %d 个: R 相关 %.3f; 第4周 Top20 在最终排名的中位位次 %.0f\n",
            nrow(cmp), cor(cmp$R_median_wk4, cmp$R_median_final),
            median(cmp[rank_wk4 <= 20]$rank_final)))
cat("完成: 11_final_drug_ranking.R\n")
