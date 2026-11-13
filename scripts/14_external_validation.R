# =====================================================================
# 14 第 10 周: GSE58294 外部队列重新计算
# 冻结要素: 基因列表(AIS 签名 185)、药物集合(1714 全部有实测特征化合物)、
#           药物响应特征(24h 主剂量 z)、计分(-Spearman)、汇总(7 细胞中位数)
# 缺失基因处理: 事先固定为"两平台共同可测基因 ∩ 949 映射 landmark"
# 主要评估背景: <=3h (23 CES / 23 对照, 采血前); 5h/24h 单独展示, 患者重复样本保持成组
# 注意: 外部队列只用于重复性检验, 不用于修改发现阶段规则
# =====================================================================
suppressPackageStartupMessages({
  library(data.table); library(limma); library(yaml)
})
tab <- "results/tables"
cfg <- yaml::read_yaml("signature/AIS_signature_config.yaml")
sig <- fread("signature/AIS_screen_signature.csv")
allg <- fread("signature/AIS_all_genes.csv")
Z <- readRDS("data/processed/drug_zmatrix.rds")
meta <- copy(Z$meta)[primary_dose == TRUE]
Zc <- Z$Zc[, meta$cond, drop = FALSE]
e58 <- readRDS("data/processed/GSE58294_expr_gene.rds")
e61 <- readRDS("data/processed/GSE16561_expr_gene.rds")
man <- fread("metadata/sample_manifest.csv")[cohort == "GSE58294"]
man <- man[match(colnames(e58), sample_id)]
stopifnot(!anyNA(man$sample_id))

# ---------- 1) 预先固定的共同基因空间 ----------
common <- Reduce(intersect, list(rownames(Zc), rownames(e58), rownames(e61)))
cat(sprintf("共同可测基因: %d (landmark 映射 %d / GSE58294 %d / GSE16561 %d)\n",
            length(common), nrow(Zc), nrow(e58), nrow(e61)))
sig_in <- sum(sig$gene %in% common)
cat(sprintf("签名 185 基因中共同空间可测 %d 个 (%.1f%%)\n", sig_in, 100 * sig_in / nrow(sig)))

cor_rank <- function(Zr, rows, tvec) {
  A <- Zr[rows, , drop = FALSE]
  A <- sweep(A, 2, colMeans(A), "-")
  b <- rank(tvec); b <- b - mean(b)
  as.numeric(colSums(A * b)) / sqrt(colSums(A^2) * sum(b^2))
}
score_external <- function(tvec, genes) {
  Zs <- Zc[genes, , drop = FALSE]
  Zr <- apply(Zs, 2, rank)
  R <- -cor_rank(Zr, genes, tvec[genes])
  dt <- data.table(parent_key = meta$parent_key, cell_id = meta$cell_id, R = R)
  dt[, .(R_med = median(R)), by = parent_key]
}

# ---------- 2) 各时间点的疾病对比 ----------
run_contrast <- function(labels, tag) {
  idx <- which(man$timepoint_label %in% labels)
  g <- factor(ifelse(man$disease_group[idx] == "AIS", "AIS", "Control"), levels = c("Control", "AIS"))
  des <- model.matrix(~ g)
  fit <- eBayes(lmFit(e58[common, idx], des), trend = TRUE)
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  tvec <- setNames(tt$t, rownames(tt))
  s <- score_external(tvec, common)
  setnames(s, "R_med", paste0("R_", tag))
  attr(s, "t") <- tvec
  attr(s, "n_up") <- sum(tt$t > 0 & p.adjust(tt$P.Value, "BH") < 0.05)
  attr(s, "n_dn") <- sum(tt$t < 0 & p.adjust(tt$P.Value, "BH") < 0.05)
  s
}
res_le3 <- run_contrast(c("control", "<=3h"), "le3h")
res_5h <- run_contrast(c("control", "5h"), "t5h")
res_24 <- run_contrast(c("control", "24h"), "t24h")
cat(sprintf("<=3h: FDR<0.05 上 %d / 下 %d; 5h: 上 %d / 下 %d; 24h: 上 %d / 下 %d\n",
            attr(res_le3, "n_up"), attr(res_le3, "n_dn"), attr(res_5h, "n_up"), attr(res_5h, "n_dn"),
            attr(res_24, "n_up"), attr(res_24, "n_dn")))

# 时间进程(仅 AIS, 患者成组): 24h - 5h
idx_tc <- which(man$disease_group == "AIS" & man$timepoint_label %in% c("5h", "24h"))
pat <- factor(man$subject_id[idx_tc])
tp <- factor(man$timepoint_label[idx_tc], levels = c("5h", "24h"))
des_tc <- model.matrix(~ pat + tp)
fit_tc <- eBayes(lmFit(e58[common, idx_tc], des_tc), trend = TRUE)
tt_tc <- topTable(fit_tc, coef = ncol(des_tc), number = Inf, sort.by = "none")
t_tc <- setNames(tt_tc$t, rownames(tt_tc))
res_tc <- score_external(t_tc, common); setnames(res_tc, "R_med", "R_tc_24h_minus_5h")

# ---------- 3) 发现阶段 R (同一共同基因空间重算, 保证可比) ----------
fin <- fread(file.path(tab, "final_drug_ranking_all.csv"))
t_disc <- allg$t[match(common, allg$gene)]; names(t_disc) <- common
disc_common <- score_external(t_disc, common); setnames(disc_common, "R_med", "R_disc_common")

x <- Reduce(function(a, b) merge(a, b, by = "parent_key", all = TRUE),
            list(fin[, .(parent_key, iname, R_primary = R_median, rank_primary = rank_overall,
                         evidence_tier, candidate_final = candidate_common, max_phase)],
                 disc_common, res_le3, res_5h, res_24, res_tc))
fwrite(x, file.path(tab, "external_validation_scores.csv"))

# ---------- 4) 指标 ----------
metrics <- list()
for (col in c("R_le3h", "R_t5h", "R_t24h", "R_tc_24h_minus_5h")) {
  ok <- !is.na(x$R_disc_common) & !is.na(x[[col]]) & x$evidence_tier == "A_common_panel"
  rho <- cor(x$R_disc_common[ok], x[[col]][ok], method = "spearman")
  sign_conc <- mean(sign(x$R_disc_common[ok]) == sign(x[[col]][ok]))
  top_d <- x[ok][order(-R_disc_common)][1:20]$parent_key
  top_e <- x[ok][order(-get(col))][1:20]$parent_key
  cand <- x[ok][candidate_final == TRUE]
  metrics[[col]] <- data.table(
    timepoint = col, n_compounds = sum(ok), rho_rank = round(rho, 3),
    sign_concordance = round(sign_conc, 3),
    top20_overlap = length(intersect(top_d, top_e)),
    candidate_median_percentile = round(median(rank(-x[[col]][ok])[x$candidate_final[ok] == TRUE]) /
                                          sum(ok), 3),
    candidate_top_quartile_pct = round(100 * mean(rank(-x[[col]][ok])[x$candidate_final[ok] == TRUE] <=
                                                    sum(ok) / 4), 1))
}
metrics <- rbindlist(metrics)
fwrite(metrics, file.path(tab, "external_validation_summary.csv")); print(metrics)

# 候选名次变化
chg <- x[evidence_tier == "A_common_panel" & candidate_final == TRUE,
         .(parent_key, iname, rank_primary, pct_le3h = round(100 * rank(-R_le3h) / .N, 1),
           pct_t5h = round(100 * rank(-R_t5h) / .N, 1), pct_t24h = round(100 * rank(-R_t24h) / .N, 1))]
setorder(chg, rank_primary)
fwrite(chg, file.path(tab, "external_candidate_rank_change.csv"))
cat(sprintf("候选(A级 %d 个): <=3h 名次百分位中位 %.1f; 进入前 25%% 的 %d 个\n",
            nrow(chg), median(chg$pct_le3h), sum(chg$pct_le3h <= 25)))

# 时间进程: 5h 与 24h 的一致性 (治疗后, 不能归因为药物作用)
ok2 <- !is.na(x$R_t5h) & !is.na(x$R_t24h) & x$evidence_tier == "A_common_panel"
cat(sprintf("5h vs 24h 药物排序相关 %.3f (n=%d); <=3h vs 24h %.3f\n",
            cor(x$R_t5h[ok2], x$R_t24h[ok2], method = "spearman"), sum(ok2),
            cor(x$R_le3h[ok2], x$R_t24h[ok2], method = "spearman")))

fwrite(data.table(common_genes = length(common), signature_genes_measurable = sig_in,
                  note = "缺失基因按事先固定的跨平台可测交集处理; 未重训、未调阈值"),
       file.path(tab, "external_validation_config.csv"))
cat("完成: 14_external_validation.R\n")
