# 03 药物特征整理与开发集排名（阶段三报告）

日期：2026-11-09 ｜ 状态：**已被 `04_model_eval_and_final_ranking.md` 取代（superseded）**

> 作废原因：本报告使用的 `08_drug_ranking.R` 中，`dose_used` 在 (parent, cell) 组内为常数，
> 导致同一化合物-细胞的不同剂量被合并进一个条件，与冻结配置“同一 (parent, cell, dose, time) 内聚合”不一致。
> 第 8 周已按配置重新聚合（脚本 `10`/`11`），本报告的共同面板排名与 37 个候选**不再有效**，
> 与最终排名的 R 相关约 0.00。文件保留以保持可追溯，结论请以 `04` 为准。

## 1. 数据层级与选型

- **首选检查**：GSE92742（L1000 Phase 1）Level5 `COMPZ.MODZ_n473647x12328`，21.3 GB gz / 解压约 25–30 GB，本机可用磁盘 31 GB，**无法容纳，记录为不可行**（写入配置）。
- **实际使用**：GSE70138（Broad LINCS 官方 CMAP02 发布，2017-03-06 构建）`Level5_COMPZ_n118050x12328`（5.4 GB gz / 解压约 6 GB）。仅化学扰动（trt_cp），1,796 个小分子。
- **Level 5 含义（必须遵守）**：同一 (化合物, 细胞, 剂量, 时间) 下重复实验合并后的**标准化修饰分数（consensus MODZ/z-score），不是原始表达量，也不是处理−对照差值**。
- 基因空间：12,328 = 978 landmark（直接测量）+ 11,350 推断。**仅使用 978 landmark 中与发现队列映射的 952 个**；推断基因不参与主指标。

## 2. 化合物身份清洗（scripts/06）

| 步骤 | 结果 |
|---|---|
| pert_info 全表 | 2,170 条 |
| trt_cp 小分子 | 1,796（本数据集无生物制剂；如出现 trt_cp 以外类型一律排除） |
| InChIKey 骨架合并（盐形式） | 1,760 个母体化合物 |
| 骨架内 InChIKey 不一致（标记复核） | 31 |
| QC 后有特征 | 1,714 |
| 开发/测试分组（母体为单位, seed=20261109） | **开发 1,232 / 测试 528（封存 `metadata/test_compounds_sealed.csv`）** |

## 3. 药物开发状态注释（scripts/07）

- 来源：ChEMBL WebResource（max_phase≥1 全量缓存 + 按 `molecule_structures__standard_inchi_key` 逐个查询），检索日期 2026-11-09。
- 匹配：InChIKey 全匹配 1,602、pref_name 归一匹配 140、未匹配 54。
- 结果：**上市（max_phase≥4）763，临床 1–3 期 582，研究/未核实 451**（记录于 `metadata/compound_development_status.csv`，含来源与日期）。

## 4. 冻结规则（signature/drug_ranking_config.yaml，在任何排名之前锁定）

- QC：trt_cp 且 `distil_cc_q75 ≥ 0.2`；主条件 **24 h**（覆盖 94%），主剂量 **10 µM**（按覆盖选，不按表现选）；化合物无 10 µM 时回退 24 h 内对数剂量最近的可用剂量，逐条记录 `dose_used`。
- 细胞：核心面板 7 个（MCF7/A375/PC3/HT29/HA1E/YAPC/HELA，各约 1.0–1.3 万条特征）。
  - 分析 A：7 个细胞各自固定条件排名，分别发布。
  - 分析 B：仅 7 细胞全部齐备的化合物，按 R 中位数汇总，报范围与方向一致比例。
  - 跨不同细胞覆盖集合的中位数仅探索性展示。
- 指标：**R = −Spearman(AIS 185 签名基因 t 值, 药物 z 分数)**；敏感性 R 用全部 952 landmark；同时报告上/下调组符号逆转比例与组内 z 中位数。
- 候选规则：R>0 且 上调组中位 z<0 且 下调组中位 z>0。
- 对照：1,000 组随机签名（同数量 100+85，按 GSE58294 表达均值十位 × GSE16561 sd 十位分箱有放回抽取 t 值；90/952 基因在 GSE58294 无探针，退化为仅 sd 分箱）→ null 分位。
- 应激/毒性核查：与"全部合格 trt_cp 共识方向"的 r_consensus>0.5 标记普遍应激驱动；median|z|>3 标记强效应；剂量 20 µM 标 high_exposure。
- **细胞背景声明**：7 个细胞均为癌细胞系/永生化细胞系，无患者免疫细胞；不得解释为"患者免疫细胞的响应"。
- **R 解读边界**：R 不是官方 tau，不是疗效概率；临床可达暴露/毒性/给药背景仅为候选解释字段。

## 5. 开发集排名结果（scripts/08，2026-11-09）

**数据规模**：QC 后合格 trt_cp 特征 57,091 / 107,404（`distil_cc_q75≥0.2`）；landmark 与阶段二映射空间并集命中 951/978 探针 → 949 个唯一基因（RPN1、LAP3 各对应 2 个 L1000 探针 ID，取 Entrez 匹配行；Entrez 命中 931 + 官方符号命中 938 取并集；**185 个签名基因全部可用**）。14 个旧别名（ERO1L→ERO1A、CTSL1→CTSL 等）经 Entrez/官方符号双路解析。

**开发组条件**（药物×细胞×剂量，24h）：9,349 条，其中核心 7 细胞 6,970 条；剂量回退（非 10µM）占 25.7%（逐条记录 `dose_used`）。

**冻结候选规则命中**：全部条件 1,189 条；核心细胞内 908 条（占 13%）。

**随机签名 null（各核心细胞，1,000 次）**：R 的 95% 分位 0.118–0.122、99% 分位 0.170–0.174——各细胞高度一致，说明阈值稳定；核心细胞条件中 **7.8% 超过各自 null 95% 分位**。

**稳健性**：R（185 签名基因）与 R_sens（949 landmark 全空间）跨条件 Spearman **0.832**——基因空间选择不改变总体排序结构。

**应激/毒性核查**：r_consensus>0.5（普遍应激驱动）3.9%；median|z|>3（强效应）0.1%；20µM 高暴露 0.07%——绝大多数条件非细胞毒性假象，但仍需逐个候选复核。

**分析 B（共同面板，7 细胞齐备 451 个化合物）**：候选 37 个——**上市 15、临床 1–3 期 6、研究/未核实 16**。上市候选按 R 中位数前 10：

| 化合物 | R 中位 | R 范围 | 7 细胞一致比例 | 备注 |
|---|---|---|---|---|
| palovarotene | 0.074 | −0.07–0.18 | 43% | RARγ 激动剂 |
| teriflunomide | 0.065 | −0.07–0.19 | 43% | 免疫抑制（DHODH） |
| lenvatinib | 0.052 | −0.22–0.16 | 57% | 多激酶 |
| aripiprazole | 0.045 | −0.20–0.15 | 29% | |
| indoprofen | 0.041 | −0.04–0.06 | 14% | NSAID |
| quinidine | 0.036 | −0.22–0.18 | 29% | |
| erythromycin-ethylsuccinate | 0.029 | −0.32–0.25 | 43% | |
| amodiaquine | 0.026 | −0.06–0.23 | 29% | |
| clofazimine | 0.023 | −0.07–0.12 | 43% | 抗麻风/免疫调节 |
| fingolimod | 0.022 | −0.09–0.23 | 14% | S1P 免疫调节 |

分析 A 各细胞 Top20 中超 null95 的上市药举例：pomalidomide（A375, 0.219）、nintedanib（MCF7, 0.307）、fostamatinib（HT29, 0.227，SYK 抑制剂）、teriflunomide（PC3, 0.185）、spironolactone/eberolimus 等（完整列表见 `top20_display_analysis_A.csv`）。共同面板 Top15 中还出现糖皮质激素 triamcinolone-acetonide（0.072）。

**解读边界（重申）**：以上为表达方向反向关联的筛选基准，不是疗效结论；fingolimod/teriflunomide 等免疫调节药靠前与"签名主要由白细胞比例驱动"的阶段二发现一致，恰恰提示可能部分反映免疫组成逆转而非神经元保护机制——将在五基因/通路解释阶段与组成调整对照签名并行判断。

## 6. 输出文件

| 文件 | 内容 |
|---|---|
| `results/tables/drug_condition_detail_dev.csv` | 开发组 药物×细胞×剂量 全部条件明细 |
| `results/tables/drug_ranking_dev_{CELL}_24h.csv` | 分析 A：7 细胞固定条件排名 |
| `results/tables/drug_ranking_dev_common_panel.csv` | 分析 B：共同面板汇总 |
| `results/tables/top20_display_analysis_{A,B}.csv` | Top20 展示（不代表有效） |
| `results/tables/drug_ranking_null_reference.csv` | 随机签名 null 分位 |
| `metadata/compound_identity_parents.csv` 等 | 身份/分组/状态注释 |

## 7. 限制

- CMAP02 为早期芯片重处理数据，扰动谱覆盖与 L1000 Phase 1 不完全等同；GSE92742 因磁盘限制未用（已记录）。
- 细胞模型非免疫细胞，且多为肿瘤系；外推到患者外周血仅是"表达方向反向关联"假设。
- 随机签名 null 只说明分数是否超出一般扰动噪声，不构成疗效检验。
- 测试组 528 个化合物的排名已封存，第 4 周不计算，至 AI 最终评估解封。
