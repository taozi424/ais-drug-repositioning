# 08 统一出图规范（v2，2026-11-11）

所有最终图表由 `scripts/fig_style.R` + `scripts/19_redraw_figures.R` 统一生成到
`results/figures_v2/`，原始数据与统计结果不变（复用 results/tables、data/processed），
仅重设视觉。规范要点：

| 项目 | 规定 |
|---|---|
| 画布宽度 | 统一 5.8 in = 14.7 cm（A4 正文宽），插入 Word 后 1:1 显示，字号即印刷字号 |
| 分辨率 | 300 dpi（1740 px 宽），印刷级 PNG |
| 图内文字 | 全英文；不设图内标题（标题由文档图注承担），避免中文字形缺失 |
| 基准字号 | 8.5 pt（刻度 7.5 pt），图例 7–7.5 pt |
| 字体 | Helvetica（无衬线） |
| 分组配色 | AIS/上调 = 暖砖红 `#A8534E`；对照/下调 = 冷灰蓝 `#4C7597`（低饱和、色盲友好） |
| 发散色带 | `#4C7597 → #F6F4F1 → #A8534E`，中点固定为统计零点或相关中位数 |
| 开发状态 | marketed `#6E7F6A` / clinical `#8C7B62` / research `#B0A79E` |
| 背景/网格 | 白底、仅浅灰主网格 `#EDEDED`、左右下轴线 `#555`；热图无网格 |
| 居中留白 | 设计宽 < 5.8 in 的图以左右留白补齐画布，保证 Word 强制撑满页宽后字号比例一致 |
| 内容映射 | QC 图改为分位数带（5–95% 线 + 25–75% 带 + 白点中位数）；箱线图不再逐样本画须 |
| 图注 | 每图在图集/报告中配英文 Figure caption：内容、关键数字、证据来源与注意事项 |

命名与编号（图集文档）：
- Group A（A1–A12）：三队列 QC（pca / boxplot / density / sample_cor × GSE16561/GSE58294/GSE22255）
- Group B（B1–B4）：疾病签名（volcano / score / GSEA / composition）
- Group C（C1–C2）：解释层（five-gene heatmap / pathway heatmap）

与 v1 的差异说明：v1 各图尺寸、字号、标题语言不统一（中英混排、dpi 150–160），
且五基因/通路热图用 pheatmap（图内中文标题在导出设备上渲染失败被迫改英文）。
v2 全部改为 ggplot 绘制，坐标/字体/配色完全可控，并去掉图内标题。
