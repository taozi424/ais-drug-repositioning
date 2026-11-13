# =====================================================================
# 06 化合物身份清洗与开发/测试分组 (GSE70138 / CMAP02 Level5)
# 规则: signature/drug_ranking_config.yaml (冻结于任何排名计算之前)
# =====================================================================
suppressPackageStartupMessages({library(data.table); library(yaml)})
set.seed(20261109)

cfg <- yaml::read_yaml("signature/drug_ranking_config.yaml")
raw_dir <- "data/raw/LINCS/GSE70138"
out_dir <- "metadata"
tab_dir <- "results/tables"
dir.create(out_dir, showWarnings = FALSE)

# ---------- 1) 读取 pert_info, 仅保留小分子 trt_cp ----------
pert <- fread(cmd = sprintf("gzip -dc %s/GSE70138_Broad_LINCS_pert_info_2017-03-06.txt.gz", raw_dir))
stopifnot(all(c("pert_id","canonical_smiles","inchi_key","pert_iname","pert_type") %in% colnames(pert)))
cp <- pert[pert_type == "trt_cp"]
cat(sprintf("pert_info 共 %d 条, trt_cp 小分子 %d 条\n", nrow(pert), nrow(cp)))

# ---------- 2) 身份标准化 ----------
cp[, iname_std := tolower(trimws(gsub("[[:space:]-]+", "-", pert_iname)))]
cp[, iname_std := gsub("^-|-$", "", iname_std)]
# InChIKey 第一段 = 分子骨架 (盐形式共享同一骨架)
ik <- cp$inchi_key
skeleton <- fifelse(!is.na(ik) & ik != "-666" & nzchar(ik),
                    sub("-.*$", "", ik), NA_character_)
cp[, parent_key := fifelse(!is.na(skeleton) & nzchar(skeleton), skeleton,
                           paste0("NOIK-", pert_id))]
n_noik <- cp[grepl("^NOIK-", parent_key), .N]
cat(sprintf("无有效 InChIKey 的记录 %d 条 (不参与骨架合并, 单列 review)\n", n_noik))

# 同一 parent_key 下多个 pert_id = 盐/重复条目, 合并为一个母体化合物
cp_by_parent <- cp[, .(n_pert_ids = uniqueN(pert_id),
                       n_inchi = uniqueN(inchi_key),
                       iname_std = min(iname_std),
                       pert_inames = paste(sort(unique(pert_iname)), collapse = ";"),
                       has_smiles = any(canonical_smiles != "-666" & nzchar(canonical_smiles))),
                   by = parent_key]
cp_by_parent[, identity_review := (n_inchi > 1)]
cat(sprintf("母体化合物 %d 个; 其中骨架内 InChIKey 不一致需复核 %d 个\n",
            nrow(cp_by_parent), sum(cp_by_parent$identity_review)))

# ---------- 3) 特征覆盖 (通过 QC 阈值前先计数, 排名后才决定纳入) ----------
si <- fread(cmd = sprintf("gzip -dc %s/GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz", raw_dir))
sm <- fread(cmd = sprintf("gzip -dc %s/GSE70138_Broad_LINCS_sig_metrics_2017-03-06.txt.gz", raw_dir))
cp_sig <- si[pert_type == "trt_cp"][sm[pert_type == "trt_cp", .(sig_id, distil_cc_q75, distil_nsample)],
                                    on = "sig_id", nomatch = 0L]
cp_sig[, qc_pass := distil_cc_q75 >= cfg$qc$distil_cc_q75_min]
cp_sig[, parent_key := cp$parent_key[match(pert_id, cp$pert_id)]]
cov <- cp_sig[qc_pass == TRUE,
              .(n_sig_qc = .N), by = .(parent_key, cell_id, pert_itime)]
cov_parent <- cov[, .(n_sig_qc_total = sum(n_sig_qc),
                      n_cells_any = uniqueN(cell_id)), by = parent_key]
core <- unlist(cfg$cells$core_panel)
cov_core <- cov[cell_id %in% core & pert_itime == cfg$qc$time_primary,
                .(n_core_cells_24h = uniqueN(cell_id)), by = parent_key]
cp_by_parent <- merge(cp_by_parent, cov_parent, by = "parent_key", all.x = TRUE)
cp_by_parent <- merge(cp_by_parent, cov_core, by = "parent_key", all.x = TRUE)
cp_by_parent[is.na(n_sig_qc_total), `:=`(n_sig_qc_total = 0L, n_cells_any = 0L)]
cp_by_parent[is.na(n_core_cells_24h), n_core_cells_24h := 0L]
cat(sprintf("QC 后特征 >0 的母体化合物 %d 个; 24h 核心细胞覆盖 >=1 个的 %d 个\n",
            sum(cp_by_parent$n_sig_qc_total > 0), sum(cp_by_parent$n_core_cells_24h > 0)))

# ---------- 4) 开发/测试分组 (母体化合物为单位, 冻结) ----------
setorder(cp_by_parent, parent_key)
n <- nrow(cp_by_parent)
idx <- sample.int(n, size = floor(0.7 * n))
cp_by_parent[, split := fifelse(seq_len(n) %in% idx, "dev", "test")]
cat(sprintf("分组: 开发 %d / 测试 %d\n",
            sum(cp_by_parent$split == "dev"), sum(cp_by_parent$split == "test")))

# ---------- 5) 输出 ----------
fwrite(cp, file.path(out_dir, "compound_perturbagens_raw.csv"))
fwrite(cp_by_parent, file.path(out_dir, "compound_identity_parents.csv"))
test_sealed <- cp_by_parent[split == "test", .(parent_key, iname_std, pert_inames)]
fwrite(test_sealed, file.path(out_dir, "test_compounds_sealed.csv"))
cat(sprintf("测试化合物清单已封存: %s (%d 个)\n",
            file.path(out_dir, "test_compounds_sealed.csv"), nrow(test_sealed)))

# 汇总表
smry <- data.table(
  item = c("pert_info_rows", "trt_cp_rows", "parent_compounds",
           "compounds_no_inchikey", "identity_review_conflicts",
           "compounds_with_qc_sig", "compounds_core_cells_24h_ge1",
           "dev_compounds", "test_compounds"),
  value = c(nrow(pert), nrow(cp), nrow(cp_by_parent), n_noik,
            sum(cp_by_parent$identity_review),
            sum(cp_by_parent$n_sig_qc_total > 0),
            sum(cp_by_parent$n_core_cells_24h > 0),
            sum(cp_by_parent$split == "dev"), sum(cp_by_parent$split == "test")))
fwrite(smry, file.path(tab_dir, "compound_identity_summary.csv"))
print(smry)
