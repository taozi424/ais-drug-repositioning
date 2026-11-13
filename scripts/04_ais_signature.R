#!/usr/bin/env Rscript
# 04_ais_signature.R
# 阶段二：在发现队列 GSE16561 内建立 AIS 整体表达特征。
#   1) limma 全转录组差异分析，校正年龄与性别（批次不可得，记为限制）
#   2) 先落盘全基因结果 AIS_all_genes.csv（logFC / 带方向 t / FDR 全保留）
#   3) 再映射到 L1000 978 landmark 基因空间，按冻结规则选签名 AIS_screen_signature.csv
#   4) 火山图 + 签名得分 + 敏感性变体（50/150、logFC 排序、不校正模型）

suppressPackageStartupMessages({
  library(data.table); library(limma); library(yaml)
  library(AnnotationDbi); library(org.Hs.eg.db)
  library(ggplot2); library(ggrepel)
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
sig_dir <- file.path(root, "signature")
fig_dir <- file.path(root, "results", "figures")
tab_dir <- file.path(root, "results", "tables")
dir.create(sig_dir, recursive = TRUE, showWarnings = FALSE)

cfg <- yaml::read_yaml(file.path(sig_dir, "AIS_signature_config.yaml"))

expr <- readRDS(file.path(root, "data/processed/GSE16561_expr_gene.rds"))
man  <- fread(file.path(root, "metadata/sample_manifest.csv"))[cohort == cfg$discovery$cohort]
man  <- man[match(colnames(expr), sample_id)]
stopifnot(!anyNA(man$age), !any(man$sex == "unknown"))

# ---------------- 1) 设计矩阵与可估计性 ----------------
man[, group := factor(disease_group, levels = c("Control", "AIS"))]
man[, sex := factor(sex, levels = c("F", "M"))]

design <- model.matrix(~ group + age + sex, data = man)
cat("设计矩阵列:", paste(colnames(design), collapse = " | "), "\n")
cat("秩 =", qr(design)$rank, " 列数 =", ncol(design),
    ifelse(qr(design)$rank == ncol(design), " -> 可估计\n", " -> 不可估计！\n"))

# 混杂程度量化（写进报告，不用于模型）
conf <- man[, .(n = .N, age_mean = round(mean(age), 1), age_sd = round(sd(age), 1),
                pct_F = round(100 * mean(sex == "F"))), by = group]
cat("\n协变量分布:\n"); print(conf)
cat("年龄组间 p =", signif(t.test(age ~ group, data = man)$p.value, 3),
    "；性别组间 p =", signif(chisq.test(table(man$sex, man$group))$p.value, 3), "\n")

# ---------------- 2) limma ----------------
fit  <- lmFit(expr, design)
fit2 <- eBayes(fit, trend = TRUE)
tt <- topTable(fit2, coef = "groupAIS", number = Inf, adjust.method = cfg$model$multiple_testing,
               sort.by = "none")
tt$gene <- rownames(tt)
setDT(tt)
setnames(tt, "adj.P.Val", "FDR")
tt[, t := round(t, 4)][, logFC := round(logFC, 4)]
tt[, P.Value := signif(P.Value, 6)][, FDR := signif(FDR, 6)]
tt[, B := round(B, 3)]

# 敏感性：不校正模型
design0 <- model.matrix(~ group, data = man)
fit0 <- eBayes(lmFit(expr, design0), trend = TRUE)
tt0 <- topTable(fit0, coef = "groupAIS", number = Inf, sort.by = "none")
tt <- merge(tt, data.table(gene = rownames(tt0), t_unadjusted = round(tt0$t, 4),
                           FDR_unadjusted = signif(tt0$adj.P.Val, 6)),
            by = "gene", all.x = TRUE, sort = FALSE)

# ---------------- 3) L1000 landmark 映射（含别名回解） ----------------
lmk <- fread(file.path(root, "reference", "GPL20573_L1000_landmark_table.csv"))
lmk <- unique(lmk[, .(landmark_symbol = pr_gene_symbol)])
universe <- rownames(expr)

direct <- intersect(lmk$landmark_symbol, universe)
missing <- setdiff(lmk$landmark_symbol, universe)

# 别名回解：用 org.Hs.eg.db 把 landmark 符号的别名对回矩阵中的符号
alias_hits <- rbindlist(lapply(missing, function(s) {
  al <- tryCatch(suppressMessages(
    AnnotationDbi::select(org.Hs.eg.db, keys = s, keytype = "SYMBOL",
                          columns = c("ALIAS", "SYMBOL"))), error = function(e) NULL)
  if (is.null(al) || !nrow(al)) return(NULL)
  hit <- intersect(unique(al$ALIAS), universe)
  if (!length(hit)) return(NULL)
  data.table(landmark_symbol = s, gene = hit[1])
}), fill = TRUE)

landmark_map <- rbind(
  data.table(landmark_symbol = direct, gene = direct, mapped_by = "exact"),
  if (nrow(alias_hits)) alias_hits[, .(landmark_symbol, gene, mapped_by = "alias")]
)
landmark_map <- unique(landmark_map, by = "gene")
fwrite(landmark_map, file.path(tab_dir, "L1000_landmark_to_discovery_genes.csv"))

cat(sprintf("\nlandmark 基因: %d | 直接匹配 %d | 别名回解 %d | 合计覆盖 %d (%.1f%%)\n",
            nrow(lmk), length(direct), nrow(alias_hits), nrow(landmark_map),
            100 * nrow(landmark_map) / nrow(lmk)))

tt[, is_landmark := gene %in% landmark_map$gene]
tt[, landmark_symbol := landmark_map$landmark_symbol[match(gene, landmark_map$gene)]]

# ---------------- 4) 全基因结果落盘（先于签名） ----------------
tt[, direction := fifelse(FDR < cfg$signature_rule$fdr_max,
                          fifelse(logFC > 0, "Up", "Down"), "NS")]
tt <- tt[order(-abs(t))]
out_cols <- c("gene", "landmark_symbol", "is_landmark", "logFC", "AveExpr", "t",
              "P.Value", "FDR", "B", "direction", "t_unadjusted", "FDR_unadjusted")
fwrite(tt[, ..out_cols], file.path(sig_dir, "AIS_all_genes.csv"))

sig_counts <- tt[, .(n_total = .N,
                     n_sig = sum(FDR < cfg$signature_rule$fdr_max),
                     n_sig_up = sum(FDR < cfg$signature_rule$fdr_max & logFC > 0),
                     n_sig_down = sum(FDR < cfg$signature_rule$fdr_max & logFC < 0),
                     n_landmark_tested = sum(is_landmark),
                     n_landmark_sig = sum(is_landmark & FDR < cfg$signature_rule$fdr_max))]
cat("\n全基因差异结果:\n"); print(sig_counts)

# ---------------- 5) 主签名（冻结规则） ----------------
r <- cfg$signature_rule
cand <- tt[is_landmark == TRUE & FDR < r$fdr_max]
up   <- cand[logFC > 0][order(-abs(t))][seq_len(min(r$max_per_direction, .N))]
down <- cand[logFC < 0][order(-abs(t))][seq_len(min(r$max_per_direction, .N))]

n_up <- nrow(up); n_down <- nrow(down)
insufficient <- n_up < r$min_per_direction || n_down < r$min_per_direction
if (insufficient) warning("签名信息不足：上调 ", n_up, " / 下调 ", n_down,
                          "，低于最小要求 ", r$min_per_direction)

signature <- rbind(
  up[, .(gene, landmark_symbol, direction = "Up", logFC, t, FDR, rank = seq_len(.N))],
  down[, .(gene, landmark_symbol, direction = "Down", logFC, t, FDR, rank = seq_len(.N))]
)
signature[, `:=`(cohort = cfg$discovery$cohort, gene_space = "L1000_landmark",
                 signature_version = "v1-frozen-2026-11-09")]
fwrite(signature, file.path(sig_dir, "AIS_screen_signature.csv"))
tt[, in_signature := gene %in% signature$gene]
fwrite(tt[, c(out_cols, "in_signature"), with = FALSE],
       file.path(sig_dir, "AIS_all_genes.csv"))  # 回写 in_signature 标记

cat(sprintf("\n签名: 上调 %d / 下调 %d | 信息充足 = %s\n", n_up, n_down, !insufficient))
cat("上调前 10:", paste(head(signature[direction == "Up"]$gene, 10), collapse = ", "), "\n")
cat("下调前 10:", paste(head(signature[direction == "Down"]$gene, 10), collapse = ", "), "\n")

# 敏感性变体：50 / 150 / logFC 排序 / 不校正模型
sens <- list()
for (k in cfg$sensitivity$max_per_direction) {
  for (rk in cfg$sensitivity$rank_by_variants) {
    ord_up <- cand[logFC > 0][order(-abs(get(rk)))][seq_len(min(k, .N))]
    ord_dn <- cand[logFC < 0][order(-abs(get(rk)))][seq_len(min(k, .N))]
    g <- c(ord_up$gene, ord_dn$gene)
    sens[[paste0("top", k, "_by_", rk)]] <- data.table(
      variant = paste0("top", k, "_by_", rk), n_up = nrow(ord_up), n_down = nrow(ord_dn),
      overlap = length(intersect(g, signature$gene)),
      jaccard = round(length(intersect(g, signature$gene)) / length(union(g, signature$gene)), 3))
  }
}
# 不校正模型下按主规则重选
cand0 <- tt[is_landmark == TRUE & FDR_unadjusted < r$fdr_max]
u0 <- cand0[logFC > 0][order(-abs(t_unadjusted))][seq_len(min(r$max_per_direction, .N))]
d0 <- cand0[logFC < 0][order(-abs(t_unadjusted))][seq_len(min(r$max_per_direction, .N))]
g0 <- c(u0$gene, d0$gene)
sens[["unadjusted_model"]] <- data.table(
  variant = "unadjusted_model", n_up = nrow(u0), n_down = nrow(d0),
  overlap = length(intersect(g0, signature$gene)),
  jaccard = round(length(intersect(g0, signature$gene)) / length(union(g0, signature$gene)), 3))
# 全共同基因空间排序（不设截断）作为对照
sens[["full_common_space"]] <- data.table(
  variant = "full_common_space", n_up = sum(cand$logFC > 0), n_down = sum(cand$logFC < 0),
  overlap = nrow(signature), jaccard = round(nrow(signature) / nrow(cand), 3))
sens_dt <- rbindlist(sens)
fwrite(sens_dt, file.path(root, "results/qc/signature_sensitivity.csv"))
cat("\n敏感性变体:\n"); print(sens_dt)

# 模型校正前后 t 的一致性
rho_models <- cor(tt$t, tt$t_unadjusted, method = "spearman")
cat(sprintf("\n校正 vs 未校正 t 的 Spearman 相关: %.3f\n", rho_models))

# ---------------- 6) 签名得分（每样本） ----------------
z <- t(scale(t(expr)))
up_g   <- intersect(signature[direction == "Up"]$gene, rownames(z))
down_g <- intersect(signature[direction == "Down"]$gene, rownames(z))
score <- data.table(
  sample_id = colnames(z),
  score = colMeans(z[up_g, , drop = FALSE]) - colMeans(z[down_g, , drop = FALSE])
)
score <- merge(score, man[, .(sample_id, disease_group, age, sex)], by = "sample_id")
fwrite(score, file.path(tab_dir, "AIS_signature_score_GSE16561.csv"))
cat("\n签名得分组间差 (AIS - Control):",
    round(mean(score[disease_group == "AIS"]$score) - mean(score[disease_group == "Control"]$score), 3),
    "；Wilcoxon p =", signif(wilcox.test(score ~ disease_group, data = score)$p.value, 3), "\n")

# ---------------- 7) 火山图 ----------------
tt[, class := fifelse(!is_landmark, "非 landmark",
                      fifelse(FDR < r$fdr_max & logFC > 0, "landmark·上调",
                              fifelse(FDR < r$fdr_max & logFC < 0, "landmark·下调", "landmark·不显著")))]
pal <- c("非 landmark" = "grey85", "landmark·不显著" = "grey60",
         "landmark·上调" = "#C0392B", "landmark·下调" = "#2E86C1")
lab <- rbind(signature[order(-abs(t))][seq_len(12)],
             signature[direction == "Down"][order(-abs(t))][seq_len(8)])

p <- ggplot(tt, aes(logFC, -log10(FDR), colour = class)) +
  geom_point(size = .7, alpha = .75) +
  scale_colour_manual(values = pal, name = NULL) +
  geom_vline(xintercept = 0, linetype = 3, colour = "grey50") +
  geom_hline(yintercept = -log10(r$fdr_max), linetype = 3, colour = "grey50") +
  geom_text_repel(data = lab, aes(logFC, -log10(FDR), label = gene),
                  inherit.aes = FALSE, size = 2.4, max.overlaps = 25,
                  min.segment.length = 0, segment.size = .2) +
  labs(title = "GSE16561 — AIS vs Control（limma，校正年龄/性别）",
       subtitle = sprintf("landmark 基因 %d 个；签名 上调 %d / 下调 %d（FDR<%.2f，|t| 排序）",
                          sum(tt$is_landmark), n_up, n_down, r$fdr_max),
       x = "log2 fold-change (AIS - Control)", y = "-log10 FDR") +
  theme_bw(base_size = 9) + theme(legend.position = "bottom")
save_png(p, file.path(fig_dir, "AIS_volcano.png"), 7.2, 5.2, 160)

# 签名得分分布
p2 <- ggplot(score, aes(disease_group, score, fill = disease_group)) +
  geom_boxplot(width = .5, outlier.size = .8) +
  geom_jitter(width = .12, size = .9, alpha = .55) +
  scale_fill_manual(values = c(AIS = "#C0392B", Control = "#2E86C1"), guide = "none") +
  labs(title = "AIS 签名得分（landmark 上调均值 − 下调均值）", x = NULL, y = "签名得分") +
  theme_bw(base_size = 9)
save_png(p2, file.path(fig_dir, "AIS_signature_score.png"), 4.2, 3.6, 160)

saveRDS(list(tt = tt, signature = signature, model = fit2, design = design,
             man = man, landmark_map = landmark_map, score = score),
        file.path(root, "data/processed/AIS_signature_objects.rds"))

cat("\n完成：AIS_all_genes.csv / AIS_screen_signature.csv 已写入 signature/\n")
