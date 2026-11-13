# =====================================================================
# 12 第 9 周: 患者级重采样稳定性
# - 病例/对照分别有放回抽取(患者为单位), 按冻结规则重建疾病签名
# - 先 200 次试运行, 再 1000 次正式; 同一脚本按轮次增量写盘
# - 每轮: limma ~group+age+sex, BH FDR<0.05, |t| 排序, 上下各至多 100,
#         任一方 <10 记为"信息不足"(无效轮次, 不删除)
# - 打分: 冻结的药物 z 条件矩阵上重算 R = -Spearman(t_sig, z)
# =====================================================================
suppressPackageStartupMessages({
  library(data.table); library(limma); library(matrixStats); library(yaml)
})
tab <- "results/tables"
set.seed(20261109)
cfg <- yaml::read_yaml("signature/AIS_signature_config.yaml")

Z <- readRDS("data/processed/drug_zmatrix.rds")
meta <- copy(Z$meta)[primary_dose == TRUE]
Zc <- Z$Zc[, meta$cond, drop = FALSE]
cpat <- paste(meta$parent_key, meta$cell_id, meta$dose_um, "24 h", sep = "___")

e61 <- readRDS("data/processed/GSE16561_expr_gene.rds")
man <- fread("metadata/sample_manifest.csv")[cohort == "GSE16561"]
man <- man[match(colnames(e61), sample_id)]
stopifnot(!anyNA(man$sample_id))
grp <- factor(man$disease_group, levels = c("Control", "AIS"))
age <- as.numeric(man$age); sex <- factor(man$sex)
case_idx <- which(grp == "AIS"); ctrl_idx <- which(grp == "Control")
cat(sprintf("队列: %d 例 / %d 对照; 基因 %d\n", length(case_idx), length(ctrl_idx), nrow(e61)))

build_signature <- function(idx) {
  E <- e61[, idx, drop = FALSE]
  g <- droplevels(grp[idx]); a <- age[idx]; s <- sex[idx]
  des <- model.matrix(~ g + a + s)
  if (qr(des)$rank < ncol(des)) return(NULL)
  fit <- eBayes(lmFit(E, des), trend = TRUE)
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  tt$FDR <- p.adjust(tt$P.Value, "BH")
  sig <- tt[tt$FDR < cfg$signature_rule$fdr_max, ]
  up <- sig[sig$t > 0, , drop = FALSE]; dn <- sig[sig$t < 0, , drop = FALSE]
  up <- head(up[order(-abs(up$t)), ], cfg$signature_rule$max_per_direction)
  dn <- head(dn[order(-abs(dn$t)), ], cfg$signature_rule$max_per_direction)
  list(up = rownames(up), dn = rownames(dn), t_up = setNames(up$t, rownames(up)),
       t_dn = setNames(dn$t, rownames(dn)),
       n_sig_total = nrow(sig), n_sig_valid = nrow(sig) >= 1)
}

score_round <- function(sg) {
  genes <- c(sg$up, sg$dn)
  ts <- c(sg$t_up, sg$t_dn)
  genes <- genes[genes %in% rownames(Zc)]
  ts <- ts[genes]
  A <- Zc[genes, , drop = FALSE]
  Ar <- apply(A, 2, rank)
  Ar <- sweep(Ar, 2, colMeans(Ar), "-")
  b <- rank(ts); b <- b - mean(b)
  R <- -as.numeric(colSums(Ar * b)) / sqrt(colSums(Ar^2) * sum(b^2))
  # 条件 -> 化合物 (7 细胞中位数)
  dt <- data.table(parent_key = meta$parent_key, R = R)
  dt[, .(R_med = median(R), n_cells = .N), by = parent_key]
}

N_TOTAL <- as.integer(Sys.getenv("N_ITER", "1000"))
rows <- list(); rank_mat <- matrix(NA_real_, nrow = 0, ncol = 0)
parents_ref <- NULL
t0 <- Sys.time()
for (it in seq_len(N_TOTAL)) {
  ci <- sample(case_idx, length(case_idx), replace = TRUE)
  co <- sample(ctrl_idx, length(ctrl_idx), replace = TRUE)
  idx <- c(ci, co)
  sg <- build_signature(idx)
  if (is.null(sg)) {
    rows[[it]] <- data.table(iter = it, valid = FALSE, n_up = NA_integer_, n_dn = NA_integer_,
                             n_sig_total = NA_integer_, reason = "design_not_estimable")
    next
  }
  ok <- length(sg$up) >= cfg$signature_rule$min_per_direction &&
        length(sg$dn) >= cfg$signature_rule$min_per_direction
  if (!ok) {
    rows[[it]] <- data.table(iter = it, valid = FALSE, n_up = length(sg$up), n_dn = length(sg$dn),
                             n_sig_total = sg$n_sig_total, reason = "insufficient_direction")
    next
  }
  s <- score_round(sg)
  if (is.null(parents_ref)) {
    parents_ref <- s$parent_key
    rank_mat <- matrix(NA_real_, nrow = length(parents_ref), ncol = N_TOTAL,
                       dimnames = list(parents_ref, paste0("it", seq_len(N_TOTAL))))
  }
  setkey(s, parent_key)
  rank_mat[, it] <- s[.(parents_ref), R_med]
  rows[[it]] <- data.table(iter = it, valid = TRUE, n_up = length(sg$up), n_dn = length(sg$dn),
                           n_sig_total = sg$n_sig_total, reason = "ok")
  if (it %% 20 == 0) {
    cat(sprintf("iter %d/%d (%.1f min) valid=%d\n", it, N_TOTAL,
                as.numeric(difftime(Sys.time(), t0, units = "mins")),
                sum(sapply(rows, function(x) x$valid[1]))))
    fwrite(rbindlist(rows), file.path(tab, "stability_round_log.csv"))
    saveRDS(list(R = rank_mat, parents = parents_ref), "data/processed/stability_ranks.rds")
  }
  if (it == 200) cat("200 次试运行完成, 继续正式 1000 次\n")
}
fwrite(rbindlist(rows), file.path(tab, "stability_round_log.csv"))
saveRDS(list(R = rank_mat, parents = parents_ref), "data/processed/stability_ranks.rds")
cat(sprintf("完成 %d 轮, 用时 %.1f min\n", N_TOTAL,
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
