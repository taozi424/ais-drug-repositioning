# 阶段二报告：AIS 整体表达特征与疾病签名

日期：2026-11-09 ｜ 发现队列：GSE16561（39 AIS / 24 对照，全血，GPL6883）｜ 签名版本：`v1-frozen-2026-11-09`

**规则冻结时序**：阈值与规则先写入 `signature/AIS_signature_config.yaml`（2026-11-09 13:05），**早于**任何 LINCS 药物排名计算——满足"看药物排名前锁定"的要求。100/10/FDR 阈值是项目拟定规则，非经验证的最优参数；预设敏感性见 §4。

---

## 1. 模型与协变量

- 主模型：`limma`，`~ group + age + sex`，对比方向 **AIS − Control**（logFC>0 = AIS 中上调），eBayes(trend=TRUE)，BH 校正。设计矩阵满秩（4/4），可估计。
- 协变量判定（实测）：
  - **年龄强混杂**（AIS 73.1±14.0 vs 对照 59.9±9.7 岁，p=4.8e-05）→ 必须校正；
  - **性别平衡**（AIS 17M/22F vs 对照 10M/14F，p=1.0）→ 纳入；
  - **批次不可校正**：GSE16561 为单系列单平台，GEO 无芯片/批次字段——如实记为限制，不假设统计模型能消除。
- 校正 vs 未校正模型的全基因 t 秩相关 0.827；未校正模型作为敏感性变体一并保存。

## 2. 全转录组差异结果（AIS_all_genes.csv，先落盘）

18,575 个基因全部保留 logFC、带方向 t、P、FDR、B、方向标签，以及 `is_landmark` / `in_signature` / 未校正模型对照列：

| 指标 | 数值 |
|---|---|
| 受检基因 | 18,575 |
| FDR<0.05 | 2,464（上调 1,487 / 下调 977） |
| L1000 landmark 受检 | 952/978（直接匹配 891 + org.Hs.eg.db 别名回解 61，覆盖 97.3%；映射表 `results/tables/L1000_landmark_to_discovery_genes.csv`） |
| landmark 中 FDR<0.05 | 238 |

**L1000 证据分级**：978 landmark 为直接测量（进入主筛选）；其余约 1.1 万 L1000 基因为推断值——`AIS_all_genes.csv` 用 `is_landmark` 列区分，推断基因不进入主签名，只在全基因层保留证据。

## 3. 筛选签名（AIS_screen_signature.csv，按冻结规则）

FDR<0.05 的共同 landmark 基因中按 |t| 各取至多 100 个：

- **上调 100 个**（满足上限）；**下调 85 个**（该方向 FDR<0.05 的 landmark 共 85 个，实际数量如实记录，未凑数、未放宽阈值）。
- 两方向均 ≥10 → 签名信息充足（`insufficient = FALSE`）。
- 上调前 10：RGS2, RAB31, **HMGCR**, PAK1, RAI17, RAB27A, HIST2H2BE, **TSC22D3**, PRCP, BIRC2；下调前 10：DFFB, TIMM22, DXS9879E, PH-4, PPIE, MIF, CDK6, PRKX, HEATR1, TARBP1（个别为旧版注释符号，已在映射表中保留对应关系）。
- 签名得分（landmark 上调均值 − 下调均值）组间差 1.93（log2 比值 z 尺度），Wilcoxon p=3.5e-14（`results/tables/AIS_signature_score_GSE16561.csv`）。

## 4. 敏感性变体（预设，全部锁定在配置中）

| 变体 | 上/下调数 | 与主签名重叠 | Jaccard |
|---|---|---|---|
| top50_by_t / by_logFC | 50/50 | 100 / 93 | 0.541 / 0.484 |
| top150_by_t / by_logFC | 150/85 | 185 | 0.787 |
| 不校正模型（~group） | 100/100 | 144 | 0.598 |
| 全共同基因空间排序 | 153/85 | 185 | 0.777 |

主规则与变体方向高度一致；top150 与全空间排序几乎就是主签名的超集，说明排序稳健。

## 5. 通路解释（全基因 t 排序，GSEA + ORA）

- **GSEA GO BP**：检验 1,761 个基因集（15–500 基因），显著 321（上调 279 / 下调 42）；**GSEA KEGG**：352 条通路，显著 184。全表 `results/tables/GSEA_GO_BP_all.csv`、`GSEA_KEGG_all.csv`。
- **免疫方向（上调）**：innate immune response（NES 2.44）、positive regulation of inflammatory response（2.47）、IL-1β/IL-6/IL-8 产生调节、趋化因子调节、中性粒细胞趋化（KEGG：中性粒细胞胞外诱捕网 NET 形成 NES 2.76、吞噬 2.76、溶酶体 2.75、FcγR 吞噬体）——即"先天免疫/吞噬细胞激活"主导。
- **代谢与翻译方向（下调）**：线粒体翻译（NES −3.97）、翻译/核糖体（KEGG Ribosome −4.25）、核糖体生物合成、剪接体、DNA 复制；代谢相关（脂质/胆固醇）未见显著上调——代谢变化更多呈**下调/重组**而非激活（焦点表 `GSEA_metabolism_immune_focus.csv`，显著 66 条）。
- **ORA**（签名基因 vs 全受检基因）：上调 2,776 个 GO BP 显著项、下调 1,957 项；KEGG ORA 上调 207 / 下调 154。

## 6. 细胞组成诊断（重要限制）

标志模块得分与签名得分的 Spearman 相关：**中性粒 +0.786、单核 +0.768、血小板 +0.613、T 细胞 −0.858、B 细胞 −0.665**；AIS 组中性粒/单核模块显著升高、T/B 显著降低（均 p<1e-4）。

**组成调整敏感性**（模型加入中性粒 + T 细胞模块得分）：landmark 显著基因从 238 → **32**（28 上 / 4 下），主签名仅保留 **19/185（10.3%）**（`results/tables/composition_adjusted_sensitivity.csv`、全基因表 `AIS_all_genes_composition_adjusted.csv`）。

**结论限定**：发现队列特征**相当大部分反映急性期白细胞比例变化**（中性粒升高、淋巴降低）。主签名仍按冻结规则用于药物反向匹配，但机制解读必须限定为"外周血关联"；后续药物解释需同时给出组成调整后的对照签名（其下调方向仅 4 个显著基因 <10，按规则标记为信息不足，不放宽阈值）。

## 7. 五基因面板（单独展示，不并入签名）

旧项目代谢相关五基因（SLC22A4/IRS2/CCR7/PDK4/ACSL1）在发现队列中全部 FDR<0.05：

| 基因 | logFC | t | FDR | 方向 | L1000 证据级别 |
|---|---|---|---|---|---|
| PDK4 | +0.920 | 6.26 | 1.0e-05 | 上调 | 非 landmark（推断/不测） |
| SLC22A4 | +0.664 | 5.17 | 1.6e-04 | 上调 | 非 landmark |
| IRS2 | +0.725 | 4.53 | 8.3e-04 | 上调 | 非 landmark |
| ACSL1 | +0.621 | 3.37 | 1.5e-02 | 上调 | 非 landmark |
| CCR7 | −0.897 | −4.20 | 2.0e-03 | 下调 | 非 landmark |

注：这五个基因均不在 L1000 landmark 内——后续与药物库比对时只能用推断基因证据（Level 3/4 INF 数据），需与 landmark 主签名分开解读。统计见 `results/tables/five_gene_panel_stats.csv`。

## 8. 本阶段交付

| 文件 | 内容 |
|---|---|
| `signature/AIS_all_genes.csv` | 18,575 基因全量差异结果（logFC/带方向 t/FDR + landmark 与签名标记） |
| `signature/AIS_screen_signature.csv` | 185 基因筛选签名（100 上 / 85 下，landmark 空间） |
| `signature/AIS_signature_config.yaml` | 冻结的队列、模型、注释版本、阈值与敏感性预设 |
| `results/tables/GSEA_*.csv`、`ORA_*.csv` | 通路结果（GO BP + KEGG，含代谢/免疫焦点表） |
| `results/tables/cell_composition_diagnostic.csv` | 细胞组成诊断 |
| `results/figures_v2/AIS_volcano.png` 等 | 火山图、签名得分、GSEA Top20、组成相关图 |
| `scripts/04_ais_signature.R`、`05_pathway_composition.R` | 全部代码（可复跑） |

## 9. 下一步

1. 外部重复性评估：GSE58294 ≤3h（23 CES vs 23 对照）检验特征方向一致性——按患者去重，不重训阈值。
2. LINCS/L1000 取数：化合物身份冻结训练/调参/测试分组后，先用主签名做反向排序，开发集试跑；测试药物的实验特征与排名封存。
3. 药物解释时并行报告两套签名（主签名 + 组成调整对照）的排名差异。
