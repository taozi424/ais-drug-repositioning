#!/usr/bin/env Rscript
# 19_redraw_figures.R
# 按统一规范重绘全部 18 张图（输出 results/figures_v2/），不改变任何统计结果，
# 只改视觉：画布宽统一 5.8 in(=14.7 cm)、300 dpi、全英文标注、无图内标题、
# 低饱和配色（暖 #A8534E / 冷 #4C7597）、浅网格。原始数据与计算结果均复用既有文件。

# 项目根目录：环境变量 AIS_ROOT 优先，否则取本脚本所在目录的上一级
root <- Sys.getenv("AIS_ROOT", unset = local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) normalizePath(file.path(dirname(f[1]), "..")) else normalizePath(getwd())
}))
source(file.path(root, "scripts/fig_style.R"))
tab <- file.path(root, "results/tables")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
cat("== A. Cohort QC ==\n")

qc <- fread(file.path(root, "results/qc/qc_per_sample.csv"))
qc[, disease_group := factor(disease_group, levels = c("Control", "AIS"))]

for (gse in c("GSE16561", "GSE58294", "GSE22255")) {
  e <- readRDS(file.path(root, "data/processed", paste0(gse, "_expr_gene.rds")))
  q <- qc[cohort == gse]
  q <- q[match(colnames(e), sample_id)]
  stopifnot(!anyNA(q$sample_id))

  # ---- A.x PCA ----
  v <- rowVars(e)
  top <- order(v, decreasing = TRUE)[seq_len(min(5000, length(v)))]
  pca <- prcomp(t(e[top, ]), scale. = TRUE, center. = TRUE)
  pv <- summary(pca)$importance[2, 1:2] * 100
  tp_lev <- unique(as.character(q$timepoint))
  shape_vals <- c(16, 17, 15, 18, 3)[seq_along(tp_lev)]
  pd <- data.table(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                   group = q$disease_group, tp = factor(q$timepoint, levels = tp_lev))
  p <- ggplot(pd, aes(PC1, PC2, colour = group, shape = tp)) +
    geom_hline(yintercept = 0, linewidth = 0.18, colour = "#E3E3E3") +
    geom_vline(xintercept = 0, linewidth = 0.18, colour = "#E3E3E3") +
    geom_point(size = 1.7, alpha = 0.9, stroke = 0.3) +
    scale_colour_manual(values = PAL_GROUP, name = "Group") +
    scale_shape_manual(values = shape_vals, name = "Timepoint") +
    labs(x = sprintf("PC1 (%.1f%%)", pv[1]), y = sprintf("PC2 (%.1f%%)", pv[2])) +
    guides(colour = guide_legend(override.aes = list(size = 1.9)),
           shape = guide_legend(override.aes = list(size = 1.9))) +
    theme_pub()
  save_fig(p, paste0(gse, "_pca.png"), 3.9, 5.2)

  # ---- A.x per-sample distribution summary (boxplot 替代表现更紧凑) ----
  qs <- t(apply(e, 2, quantile, probs = c(.05, .25, .5, .75, .95), na.rm = TRUE))
  ds <- data.table(sample_id = colnames(e), q05 = qs[, 1], q25 = qs[, 2],
                   med = qs[, 3], q75 = qs[, 4], q95 = qs[, 5])
  ds <- merge(ds, q[, .(sample_id, disease_group)], by = "sample_id")
  setorder(ds, disease_group, med)
  ds[, idx := .I]
  p <- ggplot(ds, aes(idx, colour = disease_group)) +
    geom_linerange(aes(ymin = q05, ymax = q95), linewidth = 0.35) +
    geom_linerange(aes(ymin = q25, ymax = q75), linewidth = 1.5) +
    geom_point(aes(y = med), size = 0.5, colour = "white") +
    scale_colour_manual(values = PAL_GROUP, name = "Group") +
    labs(x = "Samples (ordered by group, then median)", y = expression(log[2]~expression)) +
    theme_pub() +
    theme(panel.grid.major.x = element_blank(), legend.position = "top")
  save_fig(p, paste0(gse, "_boxplot.png"), 2.7)

  # ---- A.x density ----
  rng <- range(e, na.rm = TRUE)
  gx <- seq(rng[1], rng[2], length.out = 400)
  dmat <- sapply(seq_len(ncol(e)), function(j) {
    d <- density(e[, j], n = 400, from = rng[1], to = rng[2])
    approx(d$x, d$y, xout = gx)$y
  })
  dl <- data.table(
    x = rep(gx, times = ncol(dmat)),
    y = as.numeric(dmat),
    sample_id = rep(colnames(e), each = nrow(dmat))
  )
  dl <- merge(dl, q[, .(sample_id, disease_group)], by = "sample_id")
  pdens <- dl[, .(y = mean(y)), by = .(x, disease_group)]
  p <- ggplot(dl, aes(x, y, group = sample_id, colour = disease_group)) +
    geom_line(linewidth = 0.12, alpha = 0.30) +
    geom_line(data = pdens, aes(group = disease_group), linewidth = 0.7) +
    scale_colour_manual(values = PAL_GROUP, name = "Group") +
    labs(x = expression(log[2]~expression), y = "Density") +
    theme_pub() +
    theme(panel.grid.major = element_blank(), legend.position = "top")
  save_fig(p, paste0(gse, "_density.png"), 3.2, 4.9)

  # ---- A.x sample correlation ----
  cm <- cor(e, method = "spearman")
  ord <- order(as.integer(q$disease_group), as.character(q$timepoint), q$PC1)
  ids <- colnames(e)[ord]
  cl <- as.data.table(as.table(cm[ids, ids]))
  setnames(cl, c("s1", "s2", "rho"))
  cl[, s1 := factor(s1, levels = ids)][, s2 := factor(s2, levels = ids)]
  mids <- median(cl$rho)
  rng2 <- max(abs(cl$rho - mids))
  p <- ggplot(cl, aes(s1, s2, fill = rho)) +
    geom_raster() +
    scale_fill_gradient2(low = PAL$cool, mid = "#F6F4F1", high = PAL$warm,
                         midpoint = mids, limits = c(mids - rng2, mids + rng2),
                         name = "Spearman\nrho") +
    coord_fixed() +
    labs(x = NULL, y = NULL) +
    theme_pub(7) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank(), axis.line = element_blank(),
          legend.key.width = unit(0.22, "cm"), legend.key.height = unit(0.85, "cm"))
  save_fig(p, paste0(gse, "_sample_cor.png"), 4.4, 4.4)
}

cat("== B. Disease signature ==\n")
obj <- readRDS(file.path(root, "data/processed/AIS_signature_objects.rds"))
tt <- as.data.table(obj$tt); sig <- as.data.table(obj$signature); score <- as.data.table(obj$score)

# ---- B1 volcano ----
tt[, cls := fifelse(!is_landmark, "Other genes",
            fifelse(FDR < 0.05 & logFC > 0, "AIS up (landmark)",
            fifelse(FDR < 0.05 & logFC < 0, "AIS down (landmark)", "Landmark, n.s.")))]
lv <- c("Other genes", "Landmark, n.s.", "AIS up (landmark)", "AIS down (landmark)")
tt[, cls := factor(cls, levels = lv)]
pal_v <- setNames(c(PAL$faint, "#BFBFBF", PAL$warm, PAL$cool), lv)
lab <- rbind(sig[direction == "Up"][order(-abs(t))][1:10],
             sig[direction == "Down"][order(-abs(t))][1:8], fill = TRUE)
# 只取基因名再从 tt 取坐标，避免 sig 与 tt 都有 logFC/FDR 时 merge 出 .x/.y 后缀
lab <- merge(lab[, .(gene)], tt[, .(gene, logFC, FDR)], by = "gene")
p <- ggplot(tt, aes(logFC, -log10(FDR), colour = cls)) +
  geom_point(size = 0.55, alpha = 0.65, stroke = 0) +
  scale_colour_manual(values = pal_v, name = NULL) +
  geom_vline(xintercept = 0, linewidth = 0.25, linetype = 2, colour = "#9A9A9A") +
  geom_hline(yintercept = -log10(0.05), linewidth = 0.25, linetype = 2, colour = "#9A9A9A") +
  ggrepel::geom_text_repel(data = lab, aes(label = gene), size = 1.9, colour = "#333333",
                           max.overlaps = 30, min.segment.length = 0,
                           segment.size = 0.15, segment.colour = "#AAAAAA",
                           box.padding = 0.15, point.padding = 0.06, show.legend = FALSE) +
  labs(x = expression(log[2]~fold~change~(AIS - Control)), y = expression(-log[10]~FDR)) +
  guides(colour = guide_legend(override.aes = list(size = 1.6, alpha = 1))) +
  theme_pub() +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  theme(legend.text = element_text(size = 6.6))
save_fig(p, "AIS_volcano.png", 4.3)

# ---- B2 signature score ----
ws <- wilcox.test(score ~ disease_group, data = score)$p.value
ymax <- max(score$score) * 1.02
p <- ggplot(score, aes(disease_group, score, fill = disease_group)) +
  geom_boxplot(width = 0.42, outlier.shape = NA, linewidth = 0.3, alpha = 0.55,
               colour = "#4A4A4A") +
  geom_jitter(width = 0.10, size = 1.0, alpha = 0.75, stroke = 0,
              aes(colour = disease_group)) +
  scale_fill_manual(values = PAL_GROUP, guide = "none") +
  scale_colour_manual(values = PAL_GROUP, guide = "none") +
  annotate("segment", x = 1.18, xend = 1.82, y = ymax, yend = ymax, linewidth = 0.3) +
  annotate("text", x = 1.5, y = ymax + 0.22, size = 2.4,
           label = sprintf("Wilcoxon p = %.1e", ws)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.14))) +
  labs(x = NULL, y = "Signature score") +
  theme_pub() +
  theme(panel.grid.major.x = element_blank())
save_fig(p, "AIS_signature_score.png", 3.1, 3.3)

# ---- B3 GSEA ----
go <- fread(file.path(tab, "GSEA_GO_BP_all.csv"))[, ontology := "GO:BP"]
kg <- fread(file.path(tab, "GSEA_KEGG_all.csv"))[, ontology := "KEGG"]
fg <- rbind(go[, .(pathway, padj, NES, size, ontology)],
            kg[, .(pathway, padj, NES, size, ontology)], fill = TRUE)
fg <- fg[!is.na(NES)]
top_plot <- rbind(fg[NES > 0][order(padj)][1:10], fg[NES < 0][order(padj)][1:10])
top_plot[, dir := fifelse(NES > 0, "Higher in AIS", "Lower in AIS")]
top_plot[, lab := ifelse(nchar(pathway) > 46, paste0(substr(pathway, 1, 44), "..."), pathway)]
top_plot <- top_plot[order(NES)]
top_plot[, lab := factor(lab, levels = lab)]
p <- ggplot(top_plot, aes(NES, lab, size = -log10(padj), colour = dir)) +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "#9A9A9A", linetype = 2) +
  geom_point(alpha = 0.9, stroke = 0) +
  scale_colour_manual(values = c("Higher in AIS" = PAL$warm, "Lower in AIS" = PAL$cool),
                      name = NULL) +
  scale_size_continuous(range = c(0.7, 3.2), name = expression(-log[10]~FDR)) +
  labs(x = "Normalised enrichment score (NES)", y = NULL) +
  theme_pub(8) +
  theme(axis.text.y = element_text(size = 6.4),
        panel.grid.major.y = element_blank(),
        legend.position = "bottom", legend.direction = "horizontal",
        legend.box = "horizontal")
save_fig(p, "AIS_GSEA_top20.png", 4.9)

# ---- B4 cell composition ----
comp <- fread(file.path(tab, "cell_composition_diagnostic.csv"))
comp[, module := gsub("_", " ", module)]
comp[, dir := fifelse(rho_with_signature > 0, "Positive", "Negative")]
setorder(comp, rho_with_signature)
comp[, module := factor(module, levels = module)]
p <- ggplot(comp, aes(rho_with_signature, module, fill = dir)) +
  geom_vline(xintercept = 0, linewidth = 0.28, colour = "#7A7A7A") +
  geom_col(width = 0.58) +
  geom_text(aes(label = sprintf("%.2f", rho_with_signature),
                hjust = ifelse(rho_with_signature > 0, -0.25, 1.25)),
            size = 2.2, colour = "#3C3C3C") +
  scale_fill_manual(values = c(Positive = PAL$warm, Negative = PAL$cool), guide = "none") +
  scale_x_continuous(limits = c(-1.05, 1.05), breaks = c(-1, -0.5, 0, 0.5, 1)) +
  labs(x = "Spearman rho with AIS signature score", y = NULL) +
  theme_pub() +
  theme(panel.grid.major.y = element_blank())
save_fig(p, "AIS_cell_composition.png", 2.9, 4.6)

cat("== C. Interpretation layer ==\n")
# ---- C1 five-gene heatmap ----
s5 <- fread(file.path(tab, "five_gene_direction_summary.csv"))
m5 <- s5[, .(v = ifelse(ais_direction == "Up", 1, -1) * median_drug_z_consensus),
         by = .(iname, rank_overall, gene)]
hm <- dcast(m5, iname + rank_overall ~ gene, value.var = "v")
setorder(hm, rank_overall)
genes_ord <- c("SLC22A4", "ACSL1", "PDK4", "IRS2", "CCR7")
adir <- unique(s5[, .(gene, ais_direction)])
hm_long <- melt(hm, id.vars = c("iname", "rank_overall"), measure.vars = genes_ord,
                variable.name = "gene", value.name = "v")
hm_long[, gene := factor(gene, levels = genes_ord)]
hm_long[, rank_ord := frank(rank_overall, ties.method = "first")]
hm_long[, gene_lab := factor(gene, levels = genes_ord,
        labels = sprintf("%s (%s)", genes_ord,
                         ifelse(adir$ais_direction[match(genes_ord, adir$gene)] == "Up", "up", "down")))]
p <- ggplot(hm_long, aes(gene_lab, -rank_ord, fill = v)) +
  geom_raster() +
  scale_fill_gradient2(low = PAL$cool, mid = DIVERGING["mid"], high = PAL$warm,
                       midpoint = 0, limits = c(-1.2, 1.2), oob = squish,
                       name = "Reversal\ndirection") +
  scale_y_continuous(breaks = -c(1, 50, 100, 150, 184), labels = c(1, 50, 100, 150, 184)) +
  labs(x = NULL, y = "Candidates (ranked by R)") +
  theme_pub(7.5) +
  theme(panel.grid = element_blank(), axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 35, hjust = 0, size = 6.6),
        legend.key.height = unit(0.9, "cm"), legend.key.width = unit(0.25, "cm"))
save_fig(p, "fig_five_gene_heatmap.png", 6.4, 4.6)

# ---- C2 pathway heatmap ----
cpw <- fread(file.path(tab, "pathway_drug_candidates.csv"))
pwr <- fread(file.path(tab, "pathway_disease_reactome.csv"))
feat <- fread(file.path(tab, "featured_candidates.csv"))
show_pw <- pwr[category != "Other" & (disease_fdr_full < 0.05 | abs(disease_z) > 2) & n_universe >= 8
               ][order(category, -abs(disease_z))]
show_pw <- show_pw[, head(.SD, 10), by = category]
sub <- cpw[id %in% show_pw$id & parent_key %in% feat$parent_key]
sub <- merge(sub, show_pw[, .(id, name, category, disease_z)], by = "id", suffixes = c("", ".pw"))
sub[, drug := feat$iname[match(parent_key, feat$parent_key)]]
sub[, pw_lab := sprintf("%s %s", ifelse(disease_z > 0, "\u25b2", "\u25bc"),
                        ifelse(nchar(name) > 40, paste0(substr(name, 1, 38), "..."), name))]
lvl <- sub[!duplicated(id)][order(category, -abs(disease_z))]
sub[, pw_lab := factor(pw_lab, levels = rev(unique(lvl$pw_lab)))]
sub[, drug := factor(drug, levels = feat$iname)]
p <- ggplot(sub, aes(drug, pw_lab, fill = drug_pw_z)) +
  geom_raster() +
  scale_fill_gradient2(low = PAL$cool, mid = DIVERGING["mid"], high = PAL$warm,
                       midpoint = 0, limits = c(-2, 2), oob = squish,
                       name = "Drug pathway\nz (median)") +
  facet_grid(category ~ ., scales = "free_y", space = "free_y") +
  labs(x = NULL, y = NULL) +
  theme_pub(7.5) +
  theme(panel.grid = element_blank(), axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 40, hjust = 0, size = 6.6),
        axis.text.y = element_text(size = 5.8),
        strip.text = element_text(size = 7, hjust = 0, face = "bold", colour = "#333333"),
        strip.background = element_rect(fill = "#F4F4F4", colour = NA),
        panel.spacing = unit(0.12, "cm"),
        legend.key.height = unit(0.9, "cm"), legend.key.width = unit(0.25, "cm"))
save_fig(p, "fig_pathway_heatmap.png", 5.8)

cat("\n完成: ", FIG_DIR, "\n", sep = "")
cat("文件数: ", length(list.files(FIG_DIR, "\\.png$")), "\n", sep = "")
