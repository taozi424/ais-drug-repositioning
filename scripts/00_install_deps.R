#!/usr/bin/env Rscript
# 00_install_deps.R — 安装 R 依赖（Bioconductor + CRAN）
# 用法： Rscript scripts/00_install_deps.R
# 版本为项目实际运行环境（2026-11，R 4.5.3）。

cran <- c("data.table", "ggplot2", "matrixStats", "scales", "ggrepel",
          "pheatmap", "yaml", "digest", "ragg")
bioc <- c("limma", "rhdf5", "clusterProfiler", "fgsea",
          "org.Hs.eg.db", "AnnotationDbi", "BiocManager")

need <- function(p) !requireNamespace(p, quietly = TRUE)
todo_cran <- cran[vapply(cran, need, logical(1))]
todo_bioc <- bioc[vapply(bioc, need, logical(1))]

if (length(todo_cran)) install.packages(todo_cran, repos = "https://cloud.r-project.org")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (length(todo_bioc)) BiocManager::install(todo_bioc, ask = FALSE, update = FALSE)

pkgs <- c(cran, bioc)
cat("\n== 已安装版本 ==\n")
for (p in pkgs) {
  v <- tryCatch(as.character(packageVersion(p)), error = function(e) "缺失")
  cat(sprintf("  %-16s %s\n", p, v))
}
