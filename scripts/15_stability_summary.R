# =====================================================================
# 15 第 9 周: 重采样稳定性汇总
# 输入 data/processed/stability_ranks.rds (1000 轮 x 1714 化合物 R_median)
# 输出: 排名区间、Top20 双分母频率、候选出现频率
# =====================================================================
suppressPackageStartupMessages({ library(data.table) })
tab <- "results/tables"
S <- readRDS("data/processed/stability_ranks.rds")
R <- S$R; parents <- S$parents
log <- fread(file.path(tab, "stability_round_log.csv"))
valid <- log[valid == TRUE]$iter
R <- R[, paste0("it", valid), drop = FALSE]
cat(sprintf("轮次 %d (无效 %d, %.1f%%)\n", length(valid), nrow(log) - length(valid),
            100 * (nrow(log) - length(valid)) / nrow(log)))

final <- fread(file.path(tab, "final_drug_ranking_all.csv"))
final <- final[match(parents, parent_key)]
stopifnot(nrow(final) == nrow(R))

# 每轮在全体化合物内的排名(1 = 最强反向)
ranks <- apply(-R, 2, function(x) frank(x, ties.method = "first"))
rk <- ranks   # 化合物 x 轮次 (apply 按列, 行=化合物)
# Top20 频次 (双分母)
top20 <- rowSums(rk <= 20)
cand_round <- rowSums(R > 0)   # 该轮 R>0 的轮数(候选必要条件之一)

summ <- data.table(parent_key = parents, iname = final$iname,
                   R_median_final = final$R_median, rank_final = final$rank_overall,
                   evidence_tier = final$evidence_tier, candidate_final = final$candidate_common,
                   max_phase = final$max_phase,
                   rank_med = apply(rk, 1, median),
                   rank_p2.5 = apply(rk, 1, quantile, .025),
                   rank_p97.5 = apply(rk, 1, quantile, .975),
                   rank_iqr = apply(rk, 1, IQR),
                   top20_n = top20,
                   top20_freq_all = round(top20 / ncol(rk), 3),
                   top20_freq_valid = round(top20 / length(valid), 3),
                   R_gt0_frac = round(cand_round / ncol(R), 3))
setorder(summ, rank_final)
fwrite(summ, file.path(tab, "stability_rank_intervals.csv"))

# 汇总指标
q <- summ[, .(median_rank_iqr = median(rank_iqr),
              median_rank_range = median(rank_p97.5 - rank_p2.5),
              frac_rank_range_le_100 = mean((rank_p97.5 - rank_p2.5) <= 100))]
print(q)
top_final <- summ[rank_final <= 20]
cat(sprintf("主排名 Top20: 中位 IQR %.0f, 排名区间中位宽度 %.0f; 在 >=50%% 轮次仍居 Top20 的 %d 个\n",
            median(top_final$rank_iqr), median(top_final$rank_p97.5 - top_final$rank_p2.5),
            sum(top_final$top20_freq_all >= 0.5)))
cand <- summ[candidate_final == TRUE]
cat(sprintf("184 候选: R>0 的轮次比例中位 %.2f; 排名区间中位宽度 %.0f; Top20 出现频率中位 %.2f\n",
            median(cand$R_gt0_frac), median(cand$rank_p97.5 - cand$rank_p2.5),
            median(cand$top20_freq_all)))

fwrite(summ[rank_final <= 20], file.path(tab, "stability_top20_detail.csv"))
cat("完成: 15_stability_summary.R\n")
