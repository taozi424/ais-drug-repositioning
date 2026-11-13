#!/usr/bin/env Rscript
# 05_pathway_composition.R
# 阶段二补：通路解释与细胞组成诊断。
#   1) 全基因 t 排序做 GSEA（GO BP，fgsea），KEGG 尝试下载失败则跳过并如实报告
#   2) 签名上/下调基因分别做 ORA（enrichGO/enrichKEGG，universe = 全部受检基因）
#   3) 代谢与免疫通路专项摘要（关键词过滤）
#   4) 细胞组成诊断：白细胞标志模块得分与签名得分的相关——判断特征是否主要是细胞比例

suppressPackageStartupMessages({
  library(data.table); library(fgsea); library(clusterProfiler); library(limma)
  library(org.Hs.eg.db); library(AnnotationDbi); library(ggplot2)
})

# 项目根目录：环境变量 AIS_ROOT 优先，否则取本脚本所在目录的上一级
root <- Sys.getenv("AIS_ROOT", unset = local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) normalizePath(file.path(dirname(f[1]), "..")) else normalizePath(getwd())
}))
sig_dir <- file.path(root, "signature"); fig_dir <- file.path(root, "results/figures")
tab_dir <- file.path(root, "results/tables")
FONT <- "Hiragino Sans GB"
save_png <- function(p, file, width, height, dpi = 160) {
  p <- p + theme(text = element_text(family = FONT))
  ggsave(file, plot = p, width = width, height = height, dpi = dpi, device = ragg::agg_png)
}

obj <- readRDS(file.path(root, "data/processed/AIS_signature_objects.rds"))
tt  <- obj$tt
signature <- obj$signature

# 排序列表：全基因带方向 t（用同一批患者的全基因排序，不截断）
rnk <- setNames(tt$t, tt$gene)
rnk <- sort(rnk, decreasing = TRUE)
rnk <- rnk[!duplicated(names(rnk))]

# ---------------- GO BP 基因集（org.Hs.eg.db） ----------------
go_map <- suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db, keys(org.Hs.eg.db, "GO"),
  columns = c("SYMBOL", "ONTOLOGY"), keytype = "GO"))
bp <- as.data.table(go_map[go_map$ONTOLOGY == "BP" & !is.na(go_map$SYMBOL), ])
# 本版 org.Hs.eg.db 无 TERM 列，改从 GO.db 取术语名
gt <- as.list(GO.db::GOTERM)          # GOTerm 对象列表
term_of <- function(id) { x <- gt[[id]]; if (is.null(x)) NA_character_ else AnnotationDbi::Term(x) }
bp_ids <- unique(bp$GO)
term_tbl <- data.table(GO = bp_ids, TERM = vapply(bp_ids, term_of, character(1)))
bp <- merge(bp, term_tbl, by = "GO", all.x = TRUE)
pathways <- split(bp$SYMBOL, bp$TERM)
pathways <- lapply(pathways, unique)
cat("GO BP 基因集数:", length(pathways), "\n")

fg <- fgseaMultilevel(pathways, stats = rnk, minSize = 15, maxSize = 500, eps = 1e-50)
fg <- as.data.table(fg)
fg[, leadingEdge := vapply(leadingEdge, function(x) paste(x, collapse = "; "), character(1))]
fg <- fg[order(padj, -abs(NES))]
fg[, ontology := "GO_BP"]
fwrite(fg, file.path(tab_dir, "GSEA_GO_BP_all.csv"))

# ---------------- KEGG（联网尝试，失败则如实跳过） ----------------
kegg_ok <- TRUE
kegg <- tryCatch({
  kk <- suppressMessages(download_KEGG("hsa"))
  t2g <- as.data.table(kk$KEGGPATHID2EXTID)   # 列: from(通路ID), to(Entrez ID)
  t2n <- as.data.table(kk$KEGGPATHID2NAME)    # 列: from, name
  ez <- suppressMessages(bitr(unique(t2g$to), fromType = "ENTREZID", toType = "SYMBOL",
                              OrgDb = org.Hs.eg.db))
  t2g[, gene_symbol := ez$SYMBOL[match(to, ez$ENTREZID)]]
  t2g <- t2g[!is.na(gene_symbol)]
  pl <- split(t2g$gene_symbol, t2g$from)
  pl <- lapply(pl, unique)
  f2 <- fgseaMultilevel(pl, stats = rnk, minSize = 10, maxSize = 500, eps = 1e-50)
  f2 <- as.data.table(f2)
  f2[, pathway := t2n$to[match(pathway, t2n$from)]]   # NAME 表第 2 列名为 to，内容是通路名
  f2[, leadingEdge := vapply(leadingEdge, function(x) paste(x, collapse = "; "), character(1))]
  f2 <- f2[order(padj, -abs(NES))]
  f2[, ontology := "KEGG"]
  fwrite(f2, file.path(tab_dir, "GSEA_KEGG_all.csv"))
  f2
}, error = function(e) { kegg_ok <<- FALSE; message("KEGG 不可用: ", conditionMessage(e)); NULL })
cat("KEGG 通路库:", ifelse(kegg_ok, "已下载并完成 GSEA", "不可用（联网失败），仅 GO BP"), "\n")

# ---------------- 代谢与免疫专项 ----------------
kw <- c("metaboli", "lipid", "cholesterol", "glucose", "insulin", "fatty acid",
        "mitochond", "oxidative", "immune", "inflamm", "leukocyte", "neutrophil",
        "mononuclear", "T cell", "interleukin", "cytokine", "chemokine",
        "complement", "platelet", "coagulation")
focus <- fg[grepl(paste(kw, collapse = "|"), pathway, ignore.case = TRUE)]
if (!is.null(kegg)) focus <- rbind(focus, kegg[grepl(paste(kw, collapse = "|"), pathway, ignore.case = TRUE)], fill = TRUE)
fwrite(focus, file.path(tab_dir, "GSEA_metabolism_immune_focus.csv"))
cat("代谢/免疫相关显著通路:", nrow(focus[padj < 0.05]), "\n")

# ---------------- ORA：签名上/下调分别富集 ----------------
universe <- tt$gene
ora_res <- list()
for (d in c("Up", "Down")) {
  genes <- signature[direction == d]$gene
  eg <- suppressMessages(enrichGO(gene = genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                                  universe = universe, ont = "BP", pAdjustMethod = "BH",
                                  pvalueCutoff = 0.05, qvalueCutoff = 0.1, readable = TRUE))
  dt_eg <- as.data.table(eg@result)
  if (nrow(dt_eg)) {
    dt_eg[, direction := d]
    ora_res[[d]] <- dt_eg
    fwrite(dt_eg, file.path(tab_dir, sprintf("ORA_%s_GO_BP.csv", d)))
  }
  cat(sprintf("ORA %s: %d 个显著 GO BP 通路 (n=%d)\n", d, nrow(dt_eg), length(genes)))
}
ek <- tryCatch({
  ez <- suppressMessages(bitr(signature$gene, fromType = "SYMBOL", toType = "ENTREZID",
                              OrgDb = org.Hs.eg.db))
  res <- list()
  for (d in c("Up", "Down")) {
    ids <- ez$ENTREZID[ez$SYMBOL %in% signature[direction == d]$gene]
    u <- suppressMessages(bitr(universe, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db))$ENTREZID
    er <- suppressMessages(enrichKEGG(gene = ids, universe = u, organism = "hsa",
                                      pAdjustMethod = "BH", pvalueCutoff = 0.05))
    dt_er <- as.data.table(er@result)
    if (nrow(dt_er)) { dt_er[, direction := d]; res[[d]] <- dt_er
      fwrite(dt_er, file.path(sprintf("%s/ORA_%s_KEGG.csv", tab_dir, d))) }
    cat(sprintf("ORA KEGG %s: %d 个显著通路\n", d, nrow(dt_er)))
  }
  rbindlist(res, fill = TRUE)
}, error = function(e) { message("KEGG ORA 不可用: ", conditionMessage(e)); NULL })

# ---------------- 细胞组成诊断 ----------------
modules <- list(
  Neutrophil = c("FCGR3B", "CSF3R", "S100A8", "S100A9", "S100A12", "MMP9", "ELANE", "CEACAM8", "LTF", "CAMP"),
  T_cell     = c("CD3D", "CD3E", "CD2", "IL7R", "LTB", "CCR7", "CD27", "CD5", "LCK"),
  B_cell     = c("CD19", "MS4A1", "CD79A", "CD79B", "TCL1A", "BLK", "BANK1"),
  Monocyte   = c("CD14", "LYZ", "FCGR1A", "CD163", "MNDA", "CD300E"),
  NK_cell    = c("NKG7", "GNLY", "KLRD1", "PRF1", "NCAM1"),
  Platelet   = c("PPBP", "PF4", "ITGA2B", "GP9", "TUBB1", "CXCL4")
)
z <- t(scale(t(readRDS(file.path(root, "data/processed/GSE16561_expr_gene.rds")))))

mod_scores <- sapply(modules, function(g) {
  g <- intersect(g, rownames(z))
  if (length(g) < 3) return(rep(NA_real_, ncol(z)))
  colMeans(z[g, , drop = FALSE])
})
mod_dt <- as.data.table(mod_scores)
mod_dt[, sample_id := colnames(z)]

sc <- fread(file.path(tab_dir, "AIS_signature_score_GSE16561.csv"))
cd <- merge(sc, mod_dt, by = "sample_id")

comp <- rbindlist(lapply(names(modules), function(m) {
  x <- cd[[m]]
  if (all(is.na(x))) return(NULL)
  wt <- suppressWarnings(wilcox.test(x ~ disease_group, data = cd)$p.value)
  data.table(module = m,
             n_markers = length(intersect(modules[[m]], rownames(z))),
             median_AIS = round(median(x[cd$disease_group == "AIS"]), 3),
             median_Control = round(median(x[cd$disease_group == "Control"]), 3),
             p_group = signif(wt, 3),
             rho_with_signature = round(cor(x, cd$score, method = "spearman", use = "complete.obs"), 3))
}))
fwrite(comp, file.path(tab_dir, "cell_composition_diagnostic.csv"))
cat("\n细胞组成诊断:\n"); print(comp)

# ---------------- 图：GSEA 点图 + 组成相关 ----------------
top_plot <- rbind(fg[NES > 0][seq_len(10)], fg[NES < 0][seq_len(10)])
top_plot[, label := paste0(pathway, " (", ontology, ")")]
top_plot <- top_plot[order(NES)]
top_plot[, label := factor(label, levels = label)]
p <- ggplot(top_plot, aes(NES, label, size = -log10(padj), colour = NES > 0)) +
  geom_point(alpha = .85) +
  scale_colour_manual(values = c(`TRUE` = "#C0392B", `FALSE` = "#2E86C1"),
                      labels = c("下调 (AIS 中降低)", "上调 (AIS 中升高)"), name = NULL) +
  labs(title = "GSEA（全基因 t 排序，GO BP + KEGG）Top20",
       x = "NES", y = NULL, size = "-log10 padj") +
  theme_bw(base_size = 8) + theme(legend.position = "bottom")
save_png(p, file.path(fig_dir, "AIS_GSEA_top20.png"), 7, 5)

p2 <- ggplot(comp[!is.na(rho_with_signature)],
             aes(rho_with_signature, reorder(module, rho_with_signature))) +
  geom_col(aes(fill = rho_with_signature > 0), width = .6) +
  scale_fill_manual(values = c("TRUE" = "#C0392B", "FALSE" = "#2E86C1"), guide = "none") +
  labs(title = "细胞标志模块得分与 AIS 签名得分的 Spearman 相关",
       subtitle = "|rho| 越大越提示签名主要由细胞比例差异驱动",
       x = "Spearman rho", y = NULL) +
  theme_bw(base_size = 9)
save_png(p2, file.path(fig_dir, "AIS_cell_composition.png"), 6, 3.5)

cat("\n完成：GSEA/ORA/组成诊断已写入 results/tables/\n")

# ---------------- 组成调整敏感性分析（协议要求） ----------------
# 模型加入最强的两个细胞模块得分（中性粒 / T 细胞）作协变量，
# 观察 landmark 显著基因与主签名的保留情况。
cd2 <- cd[, .(sample_id, Neutrophil, T_cell)]
man2 <- obj$man
man2 <- merge(man2, cd2, by = "sample_id", all.x = TRUE, sort = FALSE)
expr_z <- t(scale(t(readRDS(file.path(root, "data/processed/GSE16561_expr_gene.rds")))))
man2 <- man2[match(colnames(expr_z), sample_id)]
stopifnot(!anyNA(man2$Neutrophil), !anyNA(man2$T_cell))

man2[, group := factor(disease_group, levels = c("Control", "AIS"))]
man2[, sex := factor(sex, levels = c("F", "M"))]
design_adj <- model.matrix(~ group + age + sex + Neutrophil + T_cell, data = man2)
cat("\n组成调整模型列:", paste(colnames(design_adj), collapse = " | "),
    " 秩 =", qr(design_adj)$rank, "/", ncol(design_adj), "\n")
fit_adj <- eBayes(lmFit(expr_z, design_adj), trend = TRUE)
tt_adj <- as.data.table(topTable(fit_adj, coef = "groupAIS", number = Inf, adjust.method = "BH", sort.by = "none"))
tt_adj$gene <- rownames(expr_z)
setnames(tt_adj, "adj.P.Val", "FDR")

lmk_genes <- landmark_map <- obj$landmark_map
cand_adj <- tt_adj[gene %in% lmk_genes$gene & FDR < 0.05]
n_up_adj <- sum(cand_adj$logFC > 0); n_dn_adj <- sum(cand_adj$logFC < 0)
sig_old <- signature$gene
new_sel <- rbind(
  cand_adj[logFC > 0][order(-abs(t))][seq_len(min(100, .N))],
  cand_adj[logFC < 0][order(-abs(t))][seq_len(min(100, .N))]
)
adj_summary <- data.table(
  metric = c("landmark_FDR_lt_0.05", "signature_overlap", "signature_retained_pct"),
  value = c(nrow(cand_adj), length(intersect(new_sel$gene, sig_old)),
            round(100 * length(intersect(new_sel$gene, sig_old)) / length(sig_old), 1)))
fwrite(adj_summary, file.path(tab_dir, "composition_adjusted_sensitivity.csv"))
cat("\n组成调整敏感性:\n"); print(adj_summary)
cat(sprintf("  组成调整后 landmark 显著: %d (上调 %d / 下调 %d)\n", nrow(cand_adj), n_up_adj, n_dn_adj))

fwrite(tt_adj[order(-abs(t))][, .(gene, logFC = round(logFC,4), t = round(t,4),
       P.Value = signif(P.Value,6), FDR = signif(FDR,6), in_signature = gene %in% sig_old)],
       file.path(tab_dir, "AIS_all_genes_composition_adjusted.csv"))
