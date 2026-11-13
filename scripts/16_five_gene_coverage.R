# =====================================================================
# 16 五基因覆盖表与方向表 (SLC22A4 / ACSL1 / PDK4 / IRS2 / CCR7)
# 规则 (第 11 周指令冻结):
#   * 逐基因记录证据类型: 平台直接测量(landmark) / 平台推断(bing,
#     与 L1000 同源, 不算独立实验验证) / AI 预测(第8周门槛未通过, 无) /
#     缺失(未知)。不补零; 不因药物未覆盖五基因而淘汰任何化合物。
#   * 方向判断只对有证据的格子: AIS 方向 = 发现队列 limma t 符号;
#     药物方向 = L1000 bing 推断 z 符号 (parent x cell x dose x 24h,
#     主剂量, 核心细胞), 每格标注细胞/剂量/时间与证据类型。
# 输出:
#   results/tables/five_gene_coverage.csv            覆盖表
#   results/tables/five_gene_direction_detail.csv    条件级方向表(184 候选)
#   results/tables/five_gene_direction_summary.csv   候选 x 基因汇总
#   data/processed/five_gene_z_conditions.rds        条件矩阵(供热图)
# =====================================================================
suppressPackageStartupMessages({ library(data.table); library(yaml); library(rhdf5) })
cfg <- yaml::read_yaml("signature/drug_ranking_config.yaml")
raw_dir <- "data/raw/LINCS/GSE70138"
gctx <- file.path(raw_dir, "GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx")
cells_core <- unlist(cfg$cells$core_panel)

FIVE <- data.table(
  gene = c("SLC22A4", "ACSL1", "PDK4", "IRS2", "CCR7"),
  entrez = c("6583", "2180", "5166", "8660", "1236"))
stopifnot(!duplicated(FIVE$gene), !duplicated(FIVE$entrez))

# ---------- 1) AIS 方向 (发现队列, 全转录组) ----------
allg <- fread("signature/AIS_all_genes.csv")
FIVE[, ais_t := allg$t[match(gene, allg$gene)]]
FIVE[, ais_fdr := allg$FDR[match(gene, allg$gene)]]
FIVE[, ais_direction := fifelse(is.na(ais_t), "未知", fifelse(ais_t > 0, "Up", "Down"))]

# ---------- 2) 平台证据类型 ----------
gi <- fread(cmd = sprintf("gzip -dc %s/GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz", raw_dir),
            sep = "\t")
gi5 <- gi[match(FIVE$entrez, pr_gene_id)]
FIVE[, platform_status := fifelse(is.na(gi5$pr_is_lm), "缺失(未知)",
        fifelse(gi5$pr_is_lm == 1, "平台直接测量(landmark)",
        fifelse(gi5$pr_is_bing == 1, "平台推断(bing, 与L1000同源, 非独立实验)", "缺失(未知)")))]
FIVE[, ai_predicted := "无(第8周门槛未通过, 预测分支未启用)"]
FIVE[, z_available := fifelse(platform_status == "缺失(未知)", FALSE, TRUE)]

# ---------- 3) 条件映射 (与脚本 10 同一 QC 与条件定义) ----------
si <- fread(cmd = sprintf("gzip -dc %s/GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz", raw_dir),
            sep = "\t")
sm <- fread(cmd = sprintf("gzip -dc %s/GSE70138_Broad_LINCS_sig_metrics_2017-03-06.txt.gz", raw_dir),
            sep = "\t")
cp <- si[pert_type == "trt_cp"][sm[, .(sig_id, distil_cc_q75)], on = "sig_id", nomatch = 0L,
        .(sig_id, pert_id, cell_id, pert_idose, pert_itime, distil_cc_q75)]
cp <- cp[distil_cc_q75 >= cfg$qc$distil_cc_q75_min]
pr <- fread("metadata/compound_perturbagens_raw.csv")[pert_type == "trt_cp"]
cp[, parent_key := pr$parent_key[match(pert_id, pr$pert_id)]]
cp <- cp[!is.na(parent_key)]
cp[, dose_um := suppressWarnings(as.numeric(sub(" um", "", pert_idose)))]
cp24 <- cp[pert_itime == cfg$qc$time_primary & cell_id %in% cells_core]
cp24[, cond := paste(parent_key, cell_id, dose_um, pert_itime, sep = "___")]

# 主剂量 (与 11 脚本一致: 每个最细剂量组条数最多者为该 parent 的主剂量)
ndose <- cp24[, .N, by = .(parent_key, cell_id, dose_um)]
ndose <- ndose[order(-N)][!duplicated(parent_key)]
cp24[, primary_dose := dose_um == ndose$dose_um[match(parent_key, ndose$parent_key)]]
cp24 <- cp24[primary_dose == TRUE]
cat(sprintf("主剂量 24h 核心细胞特征 %d 条 / 条件 %d\n", nrow(cp24), uniqueN(cp24$cond)))

# ---------- 4) 读取 gctx 中五基因行 (bing 推断 z) ----------
row_ids <- h5read(gctx, "0/META/ROW/id")
col_ids <- h5read(gctx, "0/META/COL/id")
frow <- match(FIVE$entrez, row_ids)
stopifnot(!anyNA(frow))
sids <- cp24$sig_id
fcol <- match(sids, col_ids)
stopifnot(!anyNA(fcol))
Zm <- h5read(gctx, "/0/DATA/0/matrix", index = list(frow, fcol))
dimnames(Zm) <- list(FIVE$gene, sids)
# consensus(普遍扰动方向): 全部 QC 合格 trt_cp 的基因级均值, 用于扣除
# L1000 全局下调偏置后再判方向 (与脚本 10 的 consensus 同一定义)
cp_all <- cp  # 全部 QC trt_cp (不限细胞/时间)
fcol_all <- match(cp_all$sig_id, col_ids)
fcol_all <- fcol_all[!is.na(fcol_all)]
Zall5 <- h5read(gctx, "/0/DATA/0/matrix", index = list(frow, fcol_all))
h5closeAll()
consensus5 <- setNames(rowMeans(Zall5), FIVE$gene)
cat(sprintf("读取五基因 z: %d x %d; consensus(均值): %s\n", nrow(Zm), ncol(Zm),
    paste(round(consensus5, 2), collapse = ", ")))

# ---------- 5) 聚合到条件 (重复 sig 取均值, 与 10 一致) ----------
conds <- unique(cp24$cond)
cond_of <- setNames(cp24$cond, cp24$sig_id)
Zc5 <- matrix(NA_real_, nrow = 5, ncol = length(conds), dimnames = list(FIVE$gene, conds))
idx <- match(colnames(Zm), names(cond_of))
grp <- cond_of[idx]
for (cc in unique(grp)) {
  j <- which(grp == cc)
  Zc5[, cc] <- if (length(j) == 1) Zm[, j] else rowMeans(Zm[, j, drop = FALSE])
}

# ---------- 6) 条件级方向表 (184 候选 x 主剂量核心细胞条件 x 5 基因) ----------
cand <- fread("results/tables/final_candidates_measured.csv")[, .(parent_key, iname, evidence_tier,
           max_phase, R_median, rank_overall)]
# 同一条件的多条 sig 已在 Zc5 聚合; 这里按条件去重 (每 cond 一行)
condtab <- unique(cp24[, .(cond, parent_key, cell_id, dose_um, pert_itime)])
detail <- as.data.table(cbind(condtab, t(Zc5[, condtab$cond, drop = FALSE])))
detail <- merge(detail, cand, by = "parent_key")
stopifnot(nrow(detail) == sum(condtab$parent_key %in% cand$parent_key))  # 内连接只保留 184 候选
long <- melt(detail,
             id.vars = c("cond", "parent_key", "cell_id", "dose_um", "pert_itime",
                         "iname", "evidence_tier", "max_phase", "R_median", "rank_overall"),
             variable.name = "gene", value.name = "drug_z")
long <- long[!is.na(drug_z)]
long <- merge(long, FIVE[, .(gene, ais_direction, platform_status, ai_predicted)], by = "gene")
long[, drug_z_consensus := drug_z - consensus5[gene]]
long[, drug_direction := fifelse(drug_z > 0, "Up", "Down")]
long[, drug_direction_consensus := fifelse(drug_z_consensus > 0, "Up", "Down")]
long[, direction_relation := fifelse(ais_direction == "未知" | is.na(drug_z), "证据不足",
        fifelse(drug_direction == ais_direction, "同向", "反向"))]
long[, direction_relation_consensus := fifelse(ais_direction == "未知" | is.na(drug_z), "证据不足",
        fifelse(drug_direction_consensus == ais_direction, "同向", "反向"))]
long[, evidence_note := paste0("L1000 bing 推断 z | ", cell_id, " | ", dose_um, " um | ", pert_itime)]
setorder(long, rank_overall, gene, cell_id)
fwrite(long, "results/tables/five_gene_direction_detail.csv")

# ---------- 7) 候选 x 基因汇总 (中位 z, 按条件计数) ----------
summ <- long[, .(n_cells = uniqueN(cond),
                 median_drug_z = median(drug_z),
                 median_drug_z_consensus = median(drug_z_consensus),
                 n_reverse = sum(direction_relation == "反向"),
                 n_reverse_consensus = sum(direction_relation_consensus == "反向"),
                 n_same = sum(direction_relation == "同向"),
                 cells_reverse = paste(sort(unique(cell_id[direction_relation == "反向"])),
                                       collapse = ";"),
                 relation = fifelse(all(direction_relation == "同向"), "全部同向",
                 fifelse(all(direction_relation == "反向"), "全部反向", "混合/部分"))),
             by = .(parent_key, iname, evidence_tier, max_phase, R_median, rank_overall,
                    gene, ais_direction, platform_status)]
summ[, evidence_type := fifelse(platform_status == "缺失(未知)", "未知", platform_status)]
setorder(summ, rank_overall, gene)
fwrite(summ, "results/tables/five_gene_direction_summary.csv")

# ---------- 8) 覆盖表 ----------
FIVE[, in_signature := allg$gene[match(gene, allg$gene)] %in%
       fread("signature/AIS_screen_signature.csv")$gene]
cov_out <- FIVE[, .(gene, entrez, ais_t = round(ais_t, 3),
                    ais_fdr = signif(ais_fdr, 3), ais_direction,
                    platform_status, ai_predicted, in_signature,
                    note = fifelse(platform_status == "平台推断(bing, 与L1000同源, 非独立实验)",
                    "z 值来自同一 L1000 数据库的 bing 推断, 不算独立实验验证; AIS 方向为全转录组实测",
                    ""))]
fwrite(cov_out, "results/tables/five_gene_coverage.csv")

saveRDS(list(Zc5 = Zc5, conds_meta = cp24[, .(cond, parent_key, cell_id, dose_um, pert_itime)]),
        "data/processed/five_gene_z_conditions.rds")

cat("\n===== 覆盖表 =====\n"); print(cov_out)
cat("\n===== 候选反向比例汇总 (按基因, 原始符号) =====\n")
print(summ[, .(gene, n_cand = .N,
               pct_reverse = round(100 * sum(n_reverse > 0) / .N, 1),
               pct_reverse_consensus = round(100 * sum(n_reverse_consensus > 0) / .N, 1)), by = gene])
cat("\n五基因全部反向(原始)的候选数:", summ[, sum(n_reverse == n_cells)], "/", nrow(summ)/5, "\n")
cat("扣除 consensus 后仍全部反向的候选数:", summ[, sum(n_reverse_consensus == n_cells)], "\n")
cat("\n反向(原始)最多的 8 个候选:\n")
print(summ[, .(inv = sum(n_reverse), inv_c = sum(n_reverse_consensus)),
           by = .(parent_key, iname, rank_overall)][order(-inv)][1:8])
cat("\n完成.\n")
