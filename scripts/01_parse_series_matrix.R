#!/usr/bin/env Rscript
# 01_parse_series_matrix.R
# 目的：把 GEO series matrix 解析成 (1) 探针级表达矩阵 RDS (2) 样本注释长表。
# 不做任何跨队列合并、不做标准化、不做 log 变换——只做"原样读取 + 结构解析"。
# 输出：
#   data/parsed/<GSE>_expr_probe.rds   探针 x 样本，数值型矩阵
#   data/parsed/<GSE>_sample_annot.csv 样本 x 字段（含 characteristics 拆解）
#   data/parsed/<GSE>_series_meta.csv  系列级元信息

suppressPackageStartupMessages({
  library(data.table)
})

# 项目根目录：环境变量 AIS_ROOT 优先，否则取本脚本所在目录的上一级
root <- Sys.getenv("AIS_ROOT", unset = local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) normalizePath(file.path(dirname(f[1]), "..")) else normalizePath(getwd())
}))
raw_dir <- file.path(root, "data", "raw")
out_dir <- file.path(root, "data", "parsed")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cohorts <- list(
  GSE16561 = file.path(raw_dir, "GSE16561", "GSE16561_series_matrix.txt"),
  GSE58294 = file.path(raw_dir, "GSE58294", "GSE58294_series_matrix.txt"),
  GSE22255 = file.path(raw_dir, "GSE22255", "GSE22255_series_matrix.txt")
)

# ---- 1. 读 header 行，拆解 "key: value" 形式的 characteristics ----
parse_series <- function(path, gse) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")

  begin <- grep("^!series_matrix_table_begin", lines)
  end   <- grep("^!series_matrix_table_end", lines)
  head_lines <- lines[seq_len(begin - 1)]

  # 系列级元信息：!Series_xxx
  series_meta <- head_lines[grepl("^!Series_", head_lines)]
  sm <- data.table(
    field = sub("^!Series_", "", sub("\t.*$", "", series_meta)),
    value = sub("^[^\t]*\t", "", series_meta)
  )
  sm[, value := gsub('^"|"$', "", value)]

  # 样本级字段（按列对齐，每行一个字段）
  sample_lines <- head_lines[grepl("^!Sample_", head_lines)]
  split_row <- function(l) {
    v <- strsplit(l, "\t", fixed = TRUE)[[1]]
    vals <- gsub('^"|"$', "", v[-1])
    list(field = sub("^!Sample_", "", v[1]), values = vals)
  }
  parsed <- lapply(sample_lines, split_row)
  fields <- vapply(parsed, `[[`, character(1), "field")
  n_col  <- length(parsed[[1]]$values)
  stopifnot(all(vapply(parsed, function(p) length(p$values), integer(1)) == n_col))

  # 基础字段（同名字段只出现一次）
  base_fields <- c("geo_accession", "title", "source_name_ch1", "description",
                   "platform_id", "data_processing", "extract_protocol_ch1",
                   "treatment_protocol_ch1", "growth_protocol_ch1")
  base <- list()
  for (f in base_fields) {
    idx <- which(fields == f)
    if (length(idx) >= 1) base[[f]] <- parsed[[idx[1]]]$values
  }

  # characteristics_ch1 可能有多个同名行，逐行拆 "key: value"
  ch_idx <- which(fields == "characteristics_ch1")
  ch_list <- list()
  ch_raw <- character(n_col)
  for (i in ch_idx) {
    vals <- parsed[[i]]$values
    keys <- sub(":.*$", "", vals)
    keys[!grepl(":", vals)] <- ""
    for (k in unique(keys[keys != ""])) {
      pos <- which(keys == k)
      if (is.null(ch_list[[k]])) ch_list[[k]] <- rep(NA_character_, n_col)
      ch_list[[k]][pos] <- sub("^[^:]*:\\s*", "", vals[pos])
    }
    ch_raw <- ifelse(vals == "", ch_raw, ifelse(ch_raw == "", vals,
                                                paste(ch_raw, vals, sep = " | ")))
  }

  annot <- data.table(sample_id = base$geo_accession,
                      title = base$title,
                      source_name = base$source_name_ch1)
  if (!is.null(base$description)) annot[, description := base$description]
  if (!is.null(base$treatment_protocol_ch1))
    annot[, treatment_protocol := base$treatment_protocol_ch1]
  annot[, platform := base$platform_id]
  annot[, characteristics_raw := ch_raw]
  for (k in names(ch_list)) {
    col <- make.names(tolower(gsub("[^A-Za-z0-9]+", "_", k)))
    annot[[col]] <- ch_list[[k]]
  }

  # 表达矩阵
  tbl_lines <- lines[(begin + 1):(end - 1)]
  dt <- fread(text = paste(tbl_lines, collapse = "\n"),
              sep = "\t", header = TRUE, data.table = FALSE,
              na.strings = c("null", "NULL", "NA", ""),
              check.names = FALSE, quote = "\"")
  id_col <- colnames(dt)[1]
  probes <- dt[[1]]
  mat <- as.matrix(dt[, -1, drop = FALSE])
  rownames(mat) <- probes
  storage.mode(mat) <- "double"

  stopifnot(identical(colnames(mat), annot$sample_id))

  list(annot = annot, expr = mat, series = sm)
}

summary_rows <- list()

for (gse in names(cohorts)) {
  message("== ", gse)
  res <- parse_series(cohorts[[gse]], gse)

  saveRDS(res$expr, file.path(out_dir, paste0(gse, "_expr_probe.rds")))
  fwrite(res$annot, file.path(out_dir, paste0(gse, "_sample_annot.csv")))
  fwrite(res$series, file.path(out_dir, paste0(gse, "_series_meta.csv")))

  vals <- res$expr[is.finite(res$expr)]
  summary_rows[[gse]] <- data.table(
    cohort = gse,
    n_samples = ncol(res$expr),
    n_probes = nrow(res$expr),
    n_missing = sum(!is.finite(res$expr)),
    pct_missing = round(100 * mean(!is.finite(res$expr)), 3),
    min = round(min(vals), 3), q01 = round(quantile(vals, .01), 3),
    median = round(median(vals), 3), q99 = round(quantile(vals, .99), 3),
    max = round(max(vals), 3),
    platform = unique(res$annot$platform)[1]
  )
  message(sprintf("  %d probes x %d samples | missing %.3f%% | range [%.2f, %.2f] | median %.2f",
                  nrow(res$expr), ncol(res$expr), 100 * mean(!is.finite(res$expr)),
                  min(vals), max(vals), median(vals)))
}

summ <- rbindlist(summary_rows)
fwrite(summ, file.path(out_dir, "matrix_summary.csv"))
cat("\n--- 矩阵摘要 ---\n"); print(summ)
