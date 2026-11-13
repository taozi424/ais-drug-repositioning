#!/usr/bin/env Rscript
# 03_preprocess_qc.R
# 目的：每个队列独立完成 探针→基因映射 → log2 判定 → 重复探针合并 → 质控指标与图。
# 铁律：不跨队列合并（急性全血 vs 慢性 PBMC 永不共矩阵）。
# 合并规则在看药物结果之前固定：每个基因保留"队列内平均表达最高"的探针（maxMean）。

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(matrixStats)
})

# 中文绘图：ragg + 系统中文字体（避免默认设备丢字形）
FONT <- "Hiragino Sans GB"
save_png <- function(p, file, width, height, dpi = 160) {
  p <- p + theme(text = element_text(family = FONT))
  ggsave(file, plot = p, width = width, height = height, dpi = dpi,
         device = ragg::agg_png)
}

# 项目根目录：环境变量 AIS_ROOT 优先，否则取本脚本所在目录的上一级
root <- Sys.getenv("AIS_ROOT", unset = local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) normalizePath(file.path(dirname(f[1]), "..")) else normalizePath(getwd())
}))
parsed_dir <- file.path(root, "data", "parsed")
proc_dir   <- file.path(root, "data", "processed")
tab_dir    <- file.path(root, "results", "tables")
fig_dir    <- file.path(root, "results", "figures")
qc_dir     <- file.path(root, "results", "qc")
for (d in c(proc_dir, tab_dir, fig_dir, qc_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

manifest <- fread(file.path(root, "metadata", "sample_manifest.csv"))

# ---------------- 平台注释 ----------------
read_gpl6883 <- function(path) {
  # Illumina 注释文本：前 27 行为 "#字段说明"，第 28 行为表头
  dt <- fread(path, sep = "\t", header = TRUE, skip = 27,
              quote = "", data.table = FALSE, fill = TRUE)
  # 兼容两种表头：series 内导出的平台表用 Symbol，GEO AnnotGPL 用 "Gene symbol"
  sym_col <- intersect(c("Symbol", "Gene symbol"), colnames(dt))[1]
  if (is.na(sym_col)) stop("GPL6883 注释缺少 Symbol / Gene symbol 列: ", path)
  data.table(probe_id = dt$ID, symbol = dt[[sym_col]])
}

read_gpl570 <- function(path) {
  # GEO AnnotGPL 格式：前 27 行为注释说明，第 28 行为表头 (ID / Gene title / Gene symbol ...)
  dt <- fread(path, sep = "\t", header = TRUE, skip = 27,
              quote = "", data.table = FALSE, fill = TRUE)
  stopifnot("Gene symbol" %in% colnames(dt), "ID" %in% colnames(dt))
  data.table(probe_id = dt$ID, symbol = dt[["Gene symbol"]])
}

anno_list <- list(
  GPL6883 = read_gpl6883(file.path(root, "data/raw/GSE16561/GPL6883-11606.txt")),
  GPL570  = read_gpl570(file.path(root, "data/raw/GPL570/GPL570.annot.gz"))
)

# ---------------- 各队列的处理参数（写死在脚本里，可审计） ----------------
cohort_cfg <- list(
  GSE16561 = list(platform = "GPL6883", already_log2 = TRUE,
                  log2_evidence = "GEO Sample_data_processing: RMA in GeneSpring, log base 2 scale; 实测 median=0, max=6.1"),
  GSE58294 = list(platform = "GPL570",  already_log2 = TRUE,
                  log2_evidence = "GEO Sample_data_processing: Partek 内部归一化, log2 尺度; 实测 median=-2.3, max=7.9"),
  GSE22255 = list(platform = "GPL570",  already_log2 = TRUE,
                  log2_evidence = "GEO Sample_data_processing: RMA; 实测 median=4.9, max=14.6")
)

qc_rows <- list(); map_rows <- list(); proc_rows <- list()

for (gse in names(cohort_cfg)) {
  cfg <- cohort_cfg[[gse]]
  message("\n===== ", gse, " (", cfg$platform, ") =====")

  expr  <- readRDS(file.path(parsed_dir, paste0(gse, "_expr_probe.rds")))
  annot <- manifest[cohort == gse]
  stopifnot(setequal(annot$sample_id, colnames(expr)))
  annot <- annot[match(colnames(expr), sample_id)]

  # --- 1) log2 判定：以"是否已是 log 尺度"为准，绝不重复转换 ---
  rng <- range(expr, na.rm = TRUE); med <- median(expr, na.rm = TRUE)
  is_log <- cfg$already_log2
  # 硬性校验：若声称已 log2，则 99 分位不应超过 30（否则几乎必然还是线性尺度）
  q99 <- quantile(expr, .99, na.rm = TRUE)
  if (is_log && q99 > 30) stop(gse, ": 声称已 log2 但 99 分位=", round(q99, 2), "，需人工复核")

  # --- 2) 缺失处理：先按探针缺失率过滤，再以探针中位数填补 ---
  na_rate <- rowMeans(is.na(expr))
  keep_na <- na_rate <= 0.20
  expr <- expr[keep_na, , drop = FALSE]
  n_probe_raw <- nrow(expr)
  if (anyNA(expr)) {
    for (i in which(rowSums(is.na(expr)) > 0)) {
      expr[i, is.na(expr[i, ])] <- median(expr[i, ], na.rm = TRUE)
    }
  }

  # --- 3) 探针 → 基因 ---
  anno <- copy(anno_list[[cfg$platform]])
  anno <- anno[!is.na(symbol) & symbol != "" & symbol != "---"]
  anno[, symbol := sub(" ?///.*$", "", symbol)]
  anno[, symbol := toupper(trimws(symbol))]
  anno <- unique(anno, by = "probe_id")

  common <- intersect(rownames(expr), anno$probe_id)
  expr <- expr[common, , drop = FALSE]
  anno <- anno[match(common, probe_id)]

  # --- 4) 重复探针合并：保留平均表达最高的探针（规则在看药物结果前冻结） ---
  mean_expr <- rowMeans(expr)
  map <- data.table(probe_id = common, symbol = anno$symbol,
                    mean_expr = round(mean_expr, 4),
                    n_probes_for_gene = NA_integer_)
  map[, n_probes_for_gene := .N, by = symbol]
  setorder(map, symbol, -mean_expr, probe_id)
  map[, kept := seq_len(.N) == 1L, by = symbol]

  kept <- map[kept == TRUE]
  expr_gene <- expr[kept$probe_id, , drop = FALSE]
  rownames(expr_gene) <- kept$symbol
  expr_gene <- expr_gene[order(rownames(expr_gene)), , drop = FALSE]

  map_rows[[gse]] <- map
  proc_rows[[gse]] <- data.table(
    cohort = gse, platform = cfg$platform,
    n_probes_input = nrow(readRDS(file.path(parsed_dir, paste0(gse, "_expr_probe.rds")))),
    n_probes_after_na_filter = n_probe_raw,
    n_probes_with_symbol = length(common),
    n_genes = nrow(expr_gene),
    n_genes_multi_probe = map[kept == TRUE & n_probes_for_gene > 1, .N],
    already_log2 = is_log,
    log2_evidence = cfg$log2_evidence,
    range_min = round(rng[1], 2), range_max = round(rng[2], 2), median = round(med, 2),
    collapse_rule = "每基因保留队列内平均表达最高的探针; 并列取 probe_id 字典序最小"
  )

  saveRDS(expr_gene, file.path(proc_dir, paste0(gse, "_expr_gene.rds")))
  fwrite(as.data.table(expr_gene, keep.rownames = "gene"),
         file.path(tab_dir, paste0(gse, "_expr_gene.csv.gz")))

  # --- 5) 样本级质控指标 ---
  # 先判断矩阵是否已按基因做基线变换（GeneSpring "normalize to median" 的产物）：
  # 若每个基因跨样本均值都接近 0，则"样本 vs 队列中位谱"的相关性指标退化，不能用。
  gene_means <- rowMeans(expr_gene)
  baseline_centered <- median(abs(gene_means)) < 0.1 && max(abs(gene_means)) < 2

  gene_med <- rowMedians(expr_gene, na.rm = TRUE)
  cors <- apply(expr_gene, 2, function(x) cor(x, gene_med, method = "spearman"))
  if (baseline_centered) cors[] <- NA_real_

  qc <- data.table(
    sample_id = colnames(expr_gene),
    cohort = gse,
    disease_group = annot$disease_group,
    timepoint = annot$timepoint_label,
    batch = annot$batch,
    baseline_centered = baseline_centered,
    median = round(apply(expr_gene, 2, median), 3),
    iqr = round(apply(expr_gene, 2, IQR), 3),
    sd = round(apply(expr_gene, 2, sd), 3),
    n_detected = apply(expr_gene, 2, function(x) sum(x > median(expr_gene))),
    cor_to_cohort_median = round(cors, 4)
  )

  # PCA（top 5000 变量基因）
  v <- rowVars(expr_gene)
  top <- order(v, decreasing = TRUE)[seq_len(min(5000, length(v)))]
  pca <- prcomp(t(expr_gene[top, ]), scale. = TRUE, center = TRUE)
  pv <- summary(pca)$importance[2, 1:3] * 100
  qc[, `:=`(PC1 = round(pca$x[, 1], 2), PC2 = round(pca$x[, 2], 2),
            PC3 = round(pca$x[, 3], 2))]

  # 异常判定：分布离群（median/IQR 超出 median±3MAD）+ PCA 稳健距离
  robust_flag <- function(x) {
    m <- median(x); s <- mad(x)
    if (s == 0) return(rep(FALSE, length(x)))
    abs(x - m) > 5 * s
  }
  qc[, flag_dist := robust_flag(median) | robust_flag(iqr)]  # 阈值在下方统一收紧为 5*MAD
  # PCA 离群：PC1-PC3 标准化后的平方和，超过卡方 97.5% 分位（df=3, 阈 9.35）判为离群
  z <- scale(pca$x[, 1:3])
  d2 <- rowSums(z^2)
  qc[, pc_d2 := round(as.numeric(d2), 2)]
  qc[, pc_distance := round(sqrt(as.numeric(d2)), 2)]
  # 采用 99.9% 分位（df=3 -> 16.27），只抓"明确的技术性离群"，不把生物学离散当异常
  qc[, flag_pca := d2 > qchisq(0.999, df = 3)]
  qc[, outlier := flag_dist | flag_pca]
  # 观察名单（不排除）：动态范围最低 10% 的阵列，留待差异分析做敏感性检验
  qc[, low_range := iqr <= quantile(iqr, 0.10)]

  # 样本间相关性（重复测量的一致性检查）
  cm <- cor(expr_gene, method = "spearman")
  qc_rows[[gse]] <- qc

  # --- 6) 图 ---
  long <- melt(as.data.table(expr_gene, keep.rownames = "gene"), id.vars = "gene",
               variable.name = "sample_id", value.name = "expr")
  long <- merge(long, qc[, .(sample_id, disease_group, timepoint)], by = "sample_id")

  p1 <- ggplot(long, aes(sample_id, expr, fill = disease_group)) +
    geom_boxplot(outlier.size = .3, linewidth = .25) +
    scale_fill_manual(values = c(AIS = "#C0392B", Control = "#2E86C1")) +
    labs(title = paste0(gse, " — 样本表达分布（", nrow(expr_gene), " 基因）"),
         x = NULL, y = "log2 表达") +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  save_png(p1, file.path(fig_dir, paste0(gse, "_boxplot.png")), 9, 4, 150)

  p2 <- ggplot(long, aes(expr, colour = disease_group, group = sample_id)) +
    geom_density(linewidth = .3, alpha = .7) +
    scale_colour_manual(values = c(AIS = "#C0392B", Control = "#2E86C1")) +
    labs(title = paste0(gse, " — 密度曲线"), x = "log2 表达", y = "密度") +
    theme_bw(base_size = 9)
  save_png(p2, file.path(fig_dir, paste0(gse, "_density.png")), 7, 4, 150)

  pd <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                   group = qc$disease_group, tp = qc$timepoint,
                   out = qc$outlier)
  p3 <- ggplot(pd, aes(PC1, PC2, colour = group, shape = tp)) +
    geom_point(size = 2.2, alpha = .85) +
    scale_colour_manual(values = c(AIS = "#C0392B", Control = "#2E86C1")) +
    labs(title = paste0(gse, " — PCA"),
         subtitle = sprintf("PC1 %.1f%% | PC2 %.1f%%", pv[1], pv[2]),
         x = sprintf("PC1 (%.1f%%)", pv[1]), y = sprintf("PC2 (%.1f%%)", pv[2])) +
    theme_bw(base_size = 9)
  save_png(p3, file.path(fig_dir, paste0(gse, "_pca.png")), 6.5, 4.6, 150)

  # 样本相关性热图
  ord <- order(qc$disease_group, qc$timepoint, qc$PC1)
  cml <- melt(as.data.table(cm[ord, ord], keep.rownames = "s1"), id.vars = "s1",
              variable.name = "s2", value.name = "rho")
  p4 <- ggplot(cml, aes(s1, s2, fill = rho)) + geom_raster() +
    scale_fill_gradient2(low = "#2E86C1", mid = "white", high = "#C0392B", midpoint = median(cml$rho)) +
    labs(title = paste0(gse, " — 样本间 Spearman 相关"), x = NULL, y = NULL) +
    theme_minimal(base_size = 7) +
    theme(axis.text = element_blank(), panel.grid = element_blank())
  save_png(p4, file.path(fig_dir, paste0(gse, "_sample_cor.png")), 5.5, 5, 150)

  message(sprintf("  %d 探针 -> %d 基因 | 基线变换=%s | 异常样本 %d | PC1 %.1f%%",
                  length(common), nrow(expr_gene), baseline_centered, sum(qc$outlier), pv[1]))
}

qc_all <- rbindlist(qc_rows)
fwrite(qc_all, file.path(qc_dir, "qc_per_sample.csv"))
fwrite(rbindlist(map_rows), file.path(tab_dir, "probe_to_gene_map.csv.gz"))
fwrite(rbindlist(proc_rows), file.path(qc_dir, "preprocessing_summary.csv"))

# 队列内重复测量一致性：GSE58294 同一患者三个时点应高度相关
g58 <- qc_all[cohort == "GSE58294"]
if (nrow(g58)) {
  e <- readRDS(file.path(proc_dir, "GSE58294_expr_gene.rds"))
  cm <- cor(e, method = "spearman")
  man <- manifest[cohort == "GSE58294"]
  subj <- setNames(man$subject_id, man$sample_id)
  pairs <- which(upper.tri(cm), arr.ind = TRUE)
  same <- subj[rownames(cm)[pairs[, 1]]] == subj[colnames(cm)[pairs[, 2]]]
  rep_stats <- data.table(
    comparison = c("same_subject_repeat", "different_subject"),
    n_pairs = c(sum(same), sum(!same)),
    median_rho = round(c(median(cm[pairs[same, ]]), median(cm[pairs[!same, ]])), 4)
  )
  fwrite(rep_stats, file.path(qc_dir, "GSE58294_repeat_concordance.csv"))
  cat("\nGSE58294 重复测量一致性:\n"); print(rep_stats)
}

cat("\n--- 预处理汇总 ---\n"); print(rbindlist(proc_rows)[, .(cohort, n_probes_input, n_probes_with_symbol, n_genes, n_genes_multi_probe)])
cat("\n--- 各队列异常样本 ---\n")
print(qc_all[, .(n = .N, n_outlier = sum(outlier),
                 baseline_centered = unique(baseline_centered),
                 min_cor = if (all(is.na(cor_to_cohort_median))) NA_real_ else round(min(cor_to_cohort_median, na.rm = TRUE), 3)),
             by = cohort])
print(qc_all[outlier == TRUE, .(sample_id, cohort, disease_group, timepoint, median, iqr, pc_distance)])
