# 阶段一报告：数据清单、样本清单与质量控制

日期：2026-11-09 ｜ 项目：基于疾病转录特征与 AI 药物扰动预测的 AIS 候选药物重定位 ｜ 数据根目录：仓库根（脚本自动定位，可用 `AIS_ROOT` 覆盖）

**通过标准达成情况**：每个样本、每个基因均可追溯（样本→`metadata/sample_manifest.csv` 的 evidence 列；基因→`results/tables/probe_to_gene_map.csv.gz` 逐一记录保留探针）；发现队列与外部评估队列用途已固定（见数据用途表）。

---

## 1. 数据用途表（冻结）

| 编号 | 数据 | 用途 | 关键约束 |
|---|---|---|---|
| D1 | GSE16561 | 发现队列：AIS 整体疾病表达特征（39 卒中 / 24 对照） | 全血；已按基因基线变换，勿再 log2 |
| D2 | GSE58294 ≤3h | 主要外部重复性评估（23 CES / 23 血管危险因素对照） | 治疗前首采；按 subject 去重 |
| D3 | GSE58294 5h/24h | 时间与治疗背景敏感性分析（同一 23 例重复采血） | 治疗后；**不算独立队列** |
| D4 | GSE22255 | 可选慢性期补充（≥6 个月 PBMC） | **不并入急性全血主分析** |
| D5 | LINCS/L1000 | 药物特征库与 AI 训练标签 | 第 2 周取数；先冻结化合物分组 |

明细见 `metadata/data_inventory.csv`。GSE58294 按要求称为"外部队列重复性评估"（该数据已用于前期研究，非未接触过的独立验证）。

## 2. 样本清单摘要

详见 `metadata/sample_manifest.csv`（195 样本、149 受试者）与 `metadata/sample_manifest_dictionary.csv`。每样本记录：样本编号、患者编号、疾病分组、卒中亚型、发病至采血时间、采血前治疗、性别、年龄、平台、批次、时点、标本类型、用途、字段来源（evidence）。缺失一律 `unknown`：

| 队列 | 样本 | 受试者 | AIS/对照 | 亚型 | 发病-采血 | 治疗信息 | 性别/年龄 | 批次 |
|---|---|---|---|---|---|---|---|---|
| GSE16561 | 63 | 63 | 39/24 | unknown | unknown | unknown | 有（36F/27M，40–96 岁） | 系列级（无芯片级信息） |
| GSE58294 | 92 | 46 | 69/23 | 心源性栓塞（CES） | <3h / 5h / 24h | 按设计推定：≤3h 治疗前，5h/24h 治疗后 | unknown（GEO 未提交） | 系列级 |
| GSE22255 | 40 | 40 | 20/20 | unknown | unknown（文献：≥6 个月） | unknown | 有（性别年龄齐全） | 杂交批次 A–F |

要点：
- GSE58294 的 23 例患者各 3 个时点共享 `subject_id`，`analysis_role` 区分 `external_replication_primary`（≤3h+对照）与 `sensitivity_repeat_measure`（5h/24h），从结构上防止把重复测量当独立队列。
- GSE16561 的亚型、发病-采血时间、治疗信息均未进入 GEO 字段 → 按规则记 `unknown`，不得回填。
- GSE58294 的"采血前治疗"为**设计推定**（原始研究首采在治疗前），GEO 无对应字段，evidence 列已注明推定依据。

## 3. 表达矩阵产出与预处理（每队列独立）

| 队列 | 平台 | 输入探针 | 有基因注释 | 基因级矩阵 | 多探针基因 | log 尺度证据 | 缺失 |
|---|---|---|---|---|---|---|---|
| GSE16561 | GPL6883 | 24,424 | 24,424 | **18,575** | 3,724 | GEO 注明 RMA/log2（GeneSpring）；实测 median=0，max=6.1 | 1.59%（null），已过滤+填补 |
| GSE58294 | GPL570 | 54,675 | 45,118 | **21,753** | 11,517 | Partek log2；实测 median=-2.33，max=7.85 | 0 |
| GSE22255 | GPL570 | 54,675 | 45,118 | **21,753** | 11,517 | RMA log2；实测 median=4.92，max=14.63 | 0 |

- **log2 判定**：三个队列均已是 log 尺度，均未做二次转换（脚本内含 99 分位硬校验）。与"先检查是否已转换"的要求一致。
- **重复探针合并规则（在看任何药物结果之前冻结）**：每个基因保留队列内平均表达最高的探针（maxMean），并列时取 probe_id 字典序最小。全部保留决策逐条记录于 `results/tables/probe_to_gene_map.csv.gz`（含 mean_expr、kept 标记）。
- **未跨队列合并**：急性全血（GSE16561/GSE58294）与慢性 PBMC（GSE22255）分别落在独立矩阵；旧流程用 ComBat 把 GSE16561+GSE22255 合并的做法在本项目中已废止。
- 跨队列基因交集（仅供后续特征匹配参考）：GSE16561∩GSE58294 = 13,182；三队列交集 = 13,182。

## 4. 关键发现：GSE16561 是基线变换后的 log2 比值

实测每个基因跨样本均值 ≈ 0（中位 |gene mean| < 0.1），即 GEO 分发矩阵已经过 GeneSpring 的"按基因/阵列中位数基线变换"。影响：
1. 样本与队列中位谱的相关性指标在该队列**不具信息量**（实测约 0.04，等于噪声），质控改用分布统计 + PCA 稳健距离；
2. 该矩阵只能表达**相对**水平，后续疾病特征应使用秩/相对量方法（如 rank-based），不能解读绝对强度；
3. 这也解释了旧流程 ComBat 合并的高风险——两队列数值语义不同。

## 5. 质量控制结果

判定规则（已写入脚本、可审计）：分布离群 = 样本 median 或 IQR 超出队列 median±5×MAD；PCA 离群 = top5000 变异基因、PC1–3 标准化平方和 > χ²(3) 99.9% 分位（16.27）。低动态范围观察名单 = 队列内 IQR 最低 10%（**仅标注、不排除**，留待差异分析做敏感性检验）。

| 队列 | 样本 | 硬性离群 | 观察名单 | 样本间 Spearman 相关（中位/最小） |
|---|---|---|---|---|
| GSE16561 | 63 | 0 | 7 | 不适用（基线变换） |
| GSE58294 | 92 | 0 | 10 | 0.97 / 0.933 |
| GSE22255 | 40 | 0 | 4 | 0.983 / 0.960 |

- **GSE58294 重复测量一致性**：同一患者跨 3 时点 69 对样本相关中位 **0.981**，不同受试者 0.970 —— subject_id 映射正确，重复测量结构成立（`results/qc/GSE58294_repeat_concordance.csv`）。
- **PCA**：GSE16561 病例/对照在前两主成分上混合（差异信号弱于技术分散，属正常，差异分析在低维之外进行）；GSE58294 对照组在 PC2 上与卒中患者清晰分离，23 例 CES 三个时点呈轨迹分布；GSE22255 分组混合。无样本孤立于主群。
- 图：`results/figures_v2/` 下每队列 boxplot / density / PCA / 样本相关热图；逐样本指标见 `results/qc/qc_per_sample.csv`。

## 6. 既有资产审计（本阶段复核结论）

- 旧差异基因结果（桌面 `ais-biomarker-ml/results/tables/limma_DEGs.csv`，48 基因）产生自 **GSE16561+GSE22255 ComBat 合并矩阵**（`scripts/04_limma_deg.R` 输入 `training_combat.rds`，`03_batch_correction_pca.R` 将两队列合并）。因 GSE22255 为慢性期 PBMC，该结果**不可复用于本项目**；疾病特征将在 GSE16561 单队列内重新计算。
- 旧配置文件（`config.yaml`：training=GSE16561+GSE22255）只能证明旧流程的写法，不能证明最终分析按新方案运行——与预判一致，不再单独采信。
- 桌面"原始数据"仅含 GSE16561 矩阵与 GPL6883 注释；GSE58294、GSE22255 本轮从 NCBI 补齐（含 GPL570 官方注释），MD5 与下载 URL 已记录于 raw 目录。

## 7. 风险与下一步

1. GSE16561 无亚型/时间/治疗信息 → 发现特征只能代表"急性全血整体"，亚型解释依赖 GSE58294（CES）。
2. GSE58294 无性别/年龄字段 → 外部重复性评估无法做人口学匹配，敏感性分析需在方法层面说明。
3. LINCS/L1000 第 2 周取数：先读元数据与小规模示例；按化合物身份冻结训练/调参/测试分组，封存测试药物的实验特征与排名；技术试跑与初筛只用开发集。
4. 观察名单 21 个样本（7+10+4）在差异分析阶段做"含/不含"敏感性对比。

## 8. 文件索引

```
metadata/data_inventory.csv              数据用途表（冻结）
metadata/sample_manifest.csv             195 样本统一样本表
metadata/sample_manifest_dictionary.csv  字段定义与 unknown 规则
metadata/cohort_composition.csv          队列构成汇总
data/raw/                                原始 series matrix + 平台注释
data/parsed/                             解析后的探针级矩阵与样本注释
data/processed/                          基因级表达矩阵 RDS（每队列独立）
results/tables/GSE*_expr_gene.csv.gz     交付用基因级矩阵
results/tables/probe_to_gene_map.csv.gz  探针→基因溯源（含保留标记）
results/qc/qc_per_sample.csv             逐样本质控指标与标记
results/qc/preprocessing_summary.csv     预处理决策汇总
results/qc/GSE58294_repeat_concordance.csv 重复测量一致性
results/figures_v2/*.png                    boxplot/密度/PCA/相关热图
scripts/01–03                            解析、样本表、预处理与质控脚本
```
