#!/usr/bin/env Rscript
# fig_style.R —— 全项目统一出图规范（供 19_redraw_figures.R 调用）
#
# 规范要点（与 docs/08_figure_spec.md 一致）：
#   1. 画布宽固定 5.8 in = 14.7 cm，等于 A4 正文宽度；插入 Word 时 1:1 显示，
#      图内字号即为实际印刷字号（不再因缩放而变小）。
#   2. dpi = 300（1740 px 宽），印刷级；不做二次重采样。
#   3. 全部英文标注、无图内标题（标题由图注 caption 承担），避免中文字形缺失并统一观感。
#   4. 低饱和配色：warm #A8534E / cool #4C7597 / 中性灰阶；色盲友好，红蓝可分。
#   5. 图形要素：仅保留浅色主网格与左/下轴线，图例无边框、字号小于正文。

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(matrixStats); library(scales)
})

FIG_W <- 5.8
FIG_DPI <- 300
FIG_DIR <- file.path(get0("root", ifnotfound = Sys.getenv("AIS_ROOT", unset = ".")),
                     "results/figures_v2")

PAL <- list(
  warm    = "#A8534E",   # AIS / up / 正向
  cool    = "#4C7597",   # Control / down / 负向
  mid     = "#9E9E9E",
  pale    = "#D9D9D9",
  faint   = "#EFEFEF",
  ink     = "#333333"
)
PAL_GROUP <- c(AIS = PAL$warm, Control = PAL$cool)
PAL_DIR   <- c(Up = PAL$warm, Down = PAL$cool)
PAL_DEV   <- c("marketed" = "#6E7F6A", "clinical" = "#8C7B62", "research" = "#B0A79E")
DIVERGING <- c(low = PAL$cool, mid = "#F6F4F1", high = PAL$warm)

theme_pub <- function(base_size = 8.5) {
  theme_minimal(base_size = base_size, base_family = "Helvetica") +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major   = element_line(linewidth = 0.18, colour = "#EDEDED"),
      panel.background   = element_rect(fill = "white", colour = NA),
      axis.line          = element_line(linewidth = 0.28, colour = "#555555"),
      axis.ticks         = element_line(linewidth = 0.28, colour = "#555555"),
      axis.text          = element_text(colour = "#3C3C3C", size = base_size - 1),
      axis.title         = element_text(colour = "#222222", size = base_size),
      plot.title         = element_blank(),
      plot.subtitle      = element_blank(),
      plot.caption       = element_blank(),
      legend.title       = element_text(size = base_size - 0.5, colour = "#222222"),
      legend.text        = element_text(size = base_size - 1, colour = "#3C3C3C"),
      legend.key.size    = unit(0.30, "cm"),
      legend.background  = element_blank(),
      legend.key         = element_blank(),
      plot.margin        = margin(3, 6, 3, 3)
    )
}

# 输出画布统一 FIG_W 宽；data_width 为图形实际设计宽，差值用左右留白补齐，
# 保证插入 Word(强制撑满页宽)后所有图的字号比例一致。
save_fig <- function(p, name, height = 4.0, data_width = FIG_W) {
  dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
  if (data_width < FIG_W) {
    pad <- (FIG_W - data_width) / 2 * 72   # inch -> pt
    p <- p + theme(plot.margin = margin(3, pad, 3, pad))
  }
  ggsave(file.path(FIG_DIR, name), p, width = FIG_W, height = height,
         dpi = FIG_DPI, device = ragg::agg_png, bg = "white")
  cat(sprintf("  %-34s %.1f in canvas, plot %.1f x %.1f in\n", name, FIG_W, data_width, height))
}
