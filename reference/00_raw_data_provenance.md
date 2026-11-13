# 原始数据溯源记录

本仓库不分发原始数据。`scripts/00_fetch_raw_data.sh` 按下列地址重新下载；
MD5 为 2026-11-09 首次取数时的校验值，可用于确认拿到的是同一版本。

## GEO 表达矩阵与平台注释

| 文件 | 来源 | 大小 | MD5 |
|---|---|---|---|
| GSE16561/GSE16561_series_matrix.txt | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE16nnn/GSE16561/ | 17 MB | c9f08d3b9b6a60bab17de08db94fb3e9 |
| GSE16561/GPL6883-11606.txt | GPL6883 AnnotGPL（https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6883/annot/） | 16 MB | 0ae280140e265fa6066ba257b7e16d8e |
| GSE58294/GSE58294_series_matrix.txt.gz | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE58nnn/GSE58294/matrix/ | 16 MB | f4fc63965431342b668d231487afdfa0 |
| GSE22255/GSE22255_series_matrix.txt.gz | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE22nnn/GSE22255/matrix/ | 8.1 MB | 8c96c36f84f9e5c257ce1ea73b1d71d7 |
| GPL570/GPL570.annot.gz | https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/ | 9.0 MB | be9fdbdb62e45257f2d3177a32f9a518 |

## L1000 / LINCS（GSE70138，Level5 COMPZ）

目录：https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/

| 文件 | 大小 | 说明 |
|---|---|---|
| GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz | 5.0 GB | 主特征库，解压后约 6 GB |
| GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz | 1.9 MB | 特征元数据（细胞/剂量/时间） |
| GSE70138_Broad_LINCS_sig_metrics_2017-03-06.txt.gz | 3.0 MB | 共识分数质量控制（distil_cc_q75） |
| GSE70138_Broad_LINCS_pert_info_2017-03-06.txt.gz | 80 KB | 化合物信息（含 InChIKey） |
| GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz | 212 KB | 基因空间（landmark / bing 标记） |
| GSE70138_Broad_LINCS_cell_info_2017-04-28.txt.gz | 2.5 KB | 细胞系信息 |

> 未使用 GSE92742：其 Level5 为 473,647 × 12,328（gz 21.3 GB，解压约 25–30 GB），
> 超出当时可用磁盘；GSE70138 为同一 Broad 官方发布的 CMAP02 构建，仅含化学扰动。

## 参考数据（已随仓库固定版本）

| 文件 | 来源 | 版本 |
|---|---|---|
| reference/ReactomePathways.gmt | https://reactome.org/download/current/ReactomePathways.gmt | 2026-10-19 下载，2,874 条通路 |
| reference/ReactomePathwaysRelation.txt | https://reactome.org/download/current/ReactomePathwaysRelation.txt | 同上，首列父通路、次列子通路 |
| reference/GPL20573_L1000_landmark_table.csv | 由 GPL20573 平台记录提取的 978 个 landmark 基因 | 978 行 |

Reactome 通路集合随版本变化会改变通路检验结果，因此固定副本入库，不随 `current` 更新。

## 化合物开发状态

ChEMBL WebResource（max_phase ≥ 1 分子列表），检索日期 2026-11-09。
查询结果缓存于 `metadata/chembl_phase_gte1_cache.json` 与
`metadata/chembl_key_lookup_cache.json`，无需联网即可复算。
