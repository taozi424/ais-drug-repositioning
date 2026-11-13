# =====================================================================
# 17 Reactome 通路层面: 疾病 vs 药物方向比较 + 两张热图
# 规则 (第 11 周指令冻结, 先于看结果):
#   * 预注册集合 = Reactome 顶级 "Metabolism" (R-HSA-1430728) 与
#     "Immune System" (R-HSA-168256) 的全部人类后代通路, 每条 ≥5 个
#     可测(949 landmark)基因; 另允许对全通路库做校正后探索性分析。
#   * 检验: 竞争性 Wilcoxon 平均秩 (背景 = 949 个可测基因)。
#     疾病侧统计量 = AIS limma t; 药物侧 = L1000 z (主剂量/核心细胞/24h,
#     候选取细胞中位)。z>0 = 通路整体上调。
#   * FDR: 预注册集合内按类别 BH; 探索性 = 全库(全部通路) BH,
#     药物侧跨候选合并后 BH。
#   * landmark 覆盖不足的通路只如实报覆盖数, 不据此判"整条通路被逆转"。
# 输出:
#   results/tables/pathway_disease_reactome.csv
#   results/tables/pathway_drug_candidates.csv
#   results/tables/pathway_agreement_summary.csv
#   results/tables/featured_candidates.csv
#   results/figures/fig_five_gene_heatmap.png
#   results/figures/fig_pathway_heatmap.png
# =====================================================================
suppressPackageStartupMessages({ library(data.table); library(matrixStats); library(pheatmap) })
Z <- readRDS("data/processed/drug_zmatrix.rds")
Zc <- Z$Zc; meta <- Z$meta; GENES <- rownames(Zc); N <- length(GENES)
cat("universe:", N, "\n")
prim <- meta[primary_dose == TRUE]
Zp <- Zc[, prim$cond, drop = FALSE]

# ---------- 0) 证据卡重点候选 (先于通路结果选择, 规则冻结) ----------
cand <- fread("results/tables/final_candidates_measured.csv")
stab <- fread("results/tables/stability_rank_intervals.csv")
ext  <- fread("results/tables/external_candidate_rank_change.csv")
feat <- merge(cand, stab[, .(parent_key, R_gt0_frac, rank_med = rank_med,
                             rank_p2.5, rank_p97.5)], by = "parent_key")
feat <- merge(feat, ext[, .(parent_key, pct_le3h)], by = "parent_key", all.x = TRUE)
sel <- feat[evidence_tier %in% c("A_common_panel", "B_partial_panel") &
              (is.na(any_stress) | any_stress == FALSE) & n_cells >= 4]
sel <- sel[order(-max_phase, -R_median)]   # 先上市/临床, 后按 R
feat5 <- head(sel[!duplicated(iname)], 5)
fwrite(feat5[, .(parent_key, iname, R_median, n_cells, cell_consistent_pct, max_phase,
                 evidence_tier, rank_overall, R_gt0_frac, rank_med, rank_p2.5, rank_p97.5,
                 pct_le3h)],
       "results/tables/featured_candidates.csv")
cat("featured:", paste(feat5$iname, collapse = ", "), "\n")

# ---------- 1) Reactome 通路集合 ----------
gmt <- fread("reference/ReactomePathways.gmt", sep = "\t", header = FALSE,
             fill = Inf)
gmt <- gmt[, lapply(.SD, as.character)]
pw <- melt(gmt, id.vars = c("V1", "V2"), variable.name = "k", value.factor = FALSE)
pw <- pw[!is.na(value) & value != ""]
setnames(pw, c("V1", "V2", "value"), c("name", "id", "gene"))
pw <- pw[, .(name, id, gene)]
pw <- pw[grepl("^R-HSA-", id)]
# 注意: 该文件为 "parent\tchild" 方向 (经人工核对: R-HSA-1430728 的行
# 列出 TCA/脂质代谢等子通路)。descendants = 沿第二列向下扩展。
rel <- fread("reference/ReactomePathwaysRelation.txt", header = FALSE,
             col.names = c("parent", "child"))
rel <- rel[grepl("^R-HSA-", child) & grepl("^R-HSA-", parent)]
desc_of <- function(root) {
  keep <- root; repeat {
    ch <- rel[parent %in% keep]$child
    ch <- setdiff(ch, keep)
    if (!length(ch)) break
    keep <- c(keep, ch)
  }
  keep
}
met_ids  <- desc_of("R-HSA-1430728")
imm_ids  <- desc_of("R-HSA-168256")
pw[, category := fifelse(id %in% met_ids, "Metabolism",
                  fifelse(id %in% imm_ids, "Immune System", "Other"))]
pw <- pw[gene %in% GENES]
pw_n <- pw[, .(n_universe = uniqueN(gene)), by = .(id, name, category)]
pw_use <- pw_n[n_universe >= 5]
cat("通路数: 全库", nrow(pw_n), "| 代谢", pw_use[category=="Metabolism", .N],
    "| 免疫", pw_use[category=="Immune System", .N], "\n")

# 稀疏成员矩阵 (949 x nP)
pw_use[, pid := .GRP, by = id]
P <- matrix(0, N, nrow(pw_use), dimnames = list(GENES, pw_use$id))
for (i in seq_len(nrow(pw_use))) {
  P[pw[id == pw_use$id[i], unique(gene)], pw_use$id[i]] <- 1
}
ng <- colSums(P)

# 竞争性 Wilcoxon 平均秩 z (向量化)
mu_of  <- function(ng, N) ng * (N + 1) / 2
sd_of  <- function(ng, N) sqrt(ng * (N - ng) * (N + 1) / 12)
z_of   <- function(W, ng, N) (W - mu_of(ng, N)) / sd_of(ng, N)
# ---------- 2) 疾病侧 ----------
allg <- fread("signature/AIS_all_genes.csv")
t_u <- allg$t[match(GENES, allg$gene)]
stopifnot(!anyNA(t_u))
zd <- z_of(as.numeric(crossprod(P, matrix(rank(t_u), ncol = 1))), ng, N)
pw_res <- pw_use[, .(id, name, category, n_universe = ng)]
pw_res[, disease_z := round(zd, 3)]
pw_res[, pv := 2 * pnorm(-abs(zd), lower.tail = FALSE)]
pw_res[, disease_fdr_category := p.adjust(pv, "BH"), by = category]
pw_res[, disease_fdr_full := p.adjust(pv, "BH")]
setorder(pw_res, disease_fdr_full)
fwrite(pw_res, "results/tables/pathway_disease_reactome.csv")
cat("\n疾病侧显著通路 (全库 FDR<0.05):", pw_res[which(disease_fdr_full < 0.05), .N], "\n")
print(pw_res[category != "Other"][order(-abs(disease_z))][1:10,
      .(name, category, n_universe, disease_z, disease_fdr_full)])

# ---------- 3) 药物侧 (候选级基因统计 + 按通路规模的置换经验 p) ----------
# 候选级基因统计 = 其主剂量核心细胞条件的秩百分位按基因取中位 (949 维),
# 再排秩做竞争性检验。p 值用置换: 通路 n 个基因在 1..949 秩下随机秩和,
# 每个通路规模 2000 次 (秩可交换 -> 对 L1000 z 的重尾/并列稳健)。
Rm <- apply(Zp, 2, rank) / N                        # 基因 x 条件, [0,1]
Gm <- matrix(NA_real_, N, length(unique(meta[primary_dose == TRUE]$parent_key)),
             dimnames = list(GENES, sort(unique(meta[primary_dose == TRUE]$parent_key))))
for (pk in colnames(Gm)) {
  Gm[, pk] <- matrixStats::rowMedians(Rm[, meta[primary_dose == TRUE &
                                                 parent_key == pk]$cond, drop = FALSE])
}
rm(Rm); invisible(gc())

# 每个通路规模的置换零分布
sizes <- unique(ng); names(sizes) <- as.character(sizes)
set.seed(20261110)
null_w <- lapply(sizes, function(n) replicate(2000, sum(sample.int(N, n))))
emp_p <- function(W, n, mu, sd) {
  nw <- null_w[[as.character(n)]]
  (1 + sum(abs(nw - mu) >= abs(W - mu))) / (length(nw) + 1)
}
mu <- mu_of(ng, N); sgv <- sd_of(ng, N)

cand_pw <- rbindlist(lapply(colnames(Gm), function(pk) {
  r <- rank(Gm[, pk])
  W <- as.numeric(crossprod(P, r))
  z <- z_of(W, ng, N)
  pv <- mapply(emp_p, W, ng, mu, sgv)
  data.table(parent_key = pk, id = names(ng), drug_pw_z = round(z, 3), pv = pv)
}))
names(ng) <- pw_use$id
cand_pw <- merge(cand_pw, pw_res[, .(id, name, category, n_universe, disease_z)], by = "id")
cand_pw[, fdr_category := p.adjust(pv, "BH"), by = .(parent_key, category)]
cand_pw[, fdr_focus_pool := p.adjust(pv, "BH")]   # 预注册集合(代谢+免疫)跨候选合并 BH
expl <- cand_pw[id %in% pw_res$id]                 # 探索性: 加全库通路
expl_all <- rbindlist(lapply(colnames(Gm), function(pk) {
  r <- rank(Gm[, pk]); W <- as.numeric(crossprod(P, r))
  data.table(parent_key = pk, id = names(ng),
             pv = mapply(emp_p, W, ng, mu, sgv))
}))
expl_all[, fdr_exploratory := p.adjust(pv, "BH")]
cand_pw <- merge(cand_pw, expl_all[, .(id, parent_key, fdr_exploratory)],
                 by = c("id", "parent_key"))
cand_pw <- merge(cand_pw, cand[, .(parent_key, iname, evidence_tier, max_phase, R_median)],
                 by = "parent_key")
cand_pw[, direction_vs_disease := fifelse(sign(drug_pw_z) != sign(disease_z), "反向", "同向")]
setorder(cand_pw, iname, category, fdr_focus_pool)
fwrite(cand_pw, "results/tables/pathway_drug_candidates.csv")
# ---------- 4) 通路一致性汇总 ----------
agr <- cand_pw[, .(n_cand = .N,
                   pct_reverse = round(100 * mean(direction_vs_disease == "反向"), 1),
                   med_z = round(median(drug_pw_z), 2)), by = .(id, name, category, disease_z)]
agr <- merge(agr, pw_res[, .(id, disease_fdr_full, n_universe)], by = "id")
setorder(agr, pct_reverse)
fwrite(agr, "results/tables/pathway_agreement_summary.csv")
cat("\n候选反向比例最高的 10 条预注册通路:\n")
print(agr[category != "Other"][order(-pct_reverse)][1:10,
      .(name, category, n_universe, disease_z, med_z, pct_reverse, disease_fdr_full)])

# ---------- 5) 热图 1: 五基因方向 (行=候选, 列=基因) ----------
fz <- readRDS("data/processed/five_gene_z_conditions.rds")
summ5 <- fread("results/tables/five_gene_direction_summary.csv")
m <- summ5[, .(v = ifelse(ais_direction == "Up", 1, -1) * median_drug_z_consensus), by = .(iname, rank_overall, gene)]
hm <- dcast(m, iname + rank_overall ~ gene, value.var = "v")
setorder(hm, rank_overall)
genes_ord <- c("SLC22A4", "ACSL1", "PDK4", "IRS2", "CCR7")
mat <- as.matrix(hm[, genes_ord, with = FALSE])
rownames(mat) <- hm$iname
ais_dir <- summ5[!duplicated(gene), .(gene, d = paste0(ais_direction, " (FDR<0.05)"))]
ais_dir[, d := ifelse(grepl("^Up", d), "Up (FDR<0.05)", "Down (FDR<0.05)")]
ann_col <- data.frame(AIS_direction = ais_dir$d[match(colnames(mat), ais_dir$gene)])
rownames(ann_col) <- colnames(mat)
ann_row <- data.frame(Development = ifelse(is.na(cand$max_phase[match(rownames(mat), cand$iname)]), "unknown",
                          ifelse(cand$max_phase[match(rownames(mat), cand$iname)] >= 4, "marketed (phase>=4)",
                          ifelse(cand$max_phase[match(rownames(mat), cand$iname)] >= 1, "clinical", "research"))))
rownames(ann_row) <- rownames(mat)
grDevices::png("results/figures/fig_five_gene_heatmap.png", width = 1600, height = 2600, res = 150, family = "Hiragino Sans GB")
pheatmap(mat, color = c("#2166AC", "#F7F7F7", "#B2182B"), breaks = c(-1.5, -0.5, 0.5, 1.5),
         cluster_rows = FALSE, cluster_cols = FALSE, border_color = "grey90",
         annotation_col = ann_col, annotation_row = ann_row, fontsize = 9, fontsize_row = 6,
         legend_breaks = c(-1, 1), legend_labels = c("same as AIS", "opposite to AIS"),
         main = "Five-gene direction, 184 candidates (v = sign(AIS t) x median drug z, consensus-corrected)\nEvidence: L1000 bing inferred z (same source as LINCS, not independent validation)")
dev.off()
cat("\nfig_five_gene_heatmap.png 完成\n")

# ---------- 6) 热图 2: 通路方向 (featured + top A/B 候选 x 预注册通路) ----------
show_pw <- pw_res[category != "Other" &
                    (disease_fdr_full < 0.05 | abs(disease_z) > 2) &
                    n_universe >= 8][order(category, -abs(disease_z))]
show_pw <- show_pw[, head(.SD, 12), by = category]   # 每类别 |疾病z| 前 12
cp2 <- cand_pw[id %in% show_pw$id & parent_key %in% feat5$parent_key]
pm <- dcast(cp2, parent_key ~ id, value.var = "drug_pw_z")
pm <- pm[, intersect(c("parent_key", show_pw$id), names(pm)), with = FALSE]
mat2 <- as.matrix(pm[, -1])
rownames(mat2) <- feat5$iname[match(pm$parent_key, feat5$parent_key)]  # as.matrix 无行名, 须用 pm 原列
colnames(mat2) <- show_pw$name[match(colnames(mat2), show_pw$id)]
mat2 <- t(mat2)  # 行=通路, 列=候选
ann2 <- data.frame(Disease_direction = ifelse(show_pw$disease_z[match(rownames(mat2), show_pw$name)] > 0,
                            "AIS up", "AIS down"))
rownames(ann2) <- rownames(mat2)
brk <- seq(-3, 3, length.out = 101)
pal <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
grDevices::png("results/figures/fig_pathway_heatmap.png", width = 1500, height = 1800, res = 150, family = "Hiragino Sans GB")
pheatmap(mat2, color = pal, breaks = brk, cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_col = data.frame(row.names = colnames(mat2), Evidence_tier =
           feat5$evidence_tier[match(colnames(mat2), feat5$iname)]),
         annotation_row = ann2, fontsize = 9, fontsize_row = 8, fontsize_col = 9,
         main = "Pre-registered Reactome pathways: drug pathway z (median across cells) x 5 featured candidates\nRow annotation = disease-side direction (competitive rank test, permutation p)")
dev.off()
cat("fig_pathway_heatmap.png 完成\n")
cat("\n全部完成.\n")
