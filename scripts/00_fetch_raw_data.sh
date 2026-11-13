#!/usr/bin/env bash
# ==============================================================================
# 00_fetch_raw_data.sh — 重新获取全部原始数据到 data/raw/
#
# 用法： bash scripts/00_fetch_raw_data.sh [--with-lincs-level5]
#
# 默认跳过 5 GB 的 L1000 Level5 矩阵（只有第 4、8-11 周的药物排名需要它）。
# 全部下载约 6 GB / 11 GB（含 Level5），网络状况好时 20-60 分钟。
# 已存在的文件会跳过，可反复运行。
# ==============================================================================
set -euo pipefail

cd "$(dirname "$0")/.."
RAW="data/raw"
GEO="https://ftp.ncbi.nlm.nih.gov/geo"
LINCS_DIR="${RAW}/LINCS/GSE70138"
mkdir -p "${RAW}/GSE16561" "${RAW}/GSE58294" "${RAW}/GSE22255" "${RAW}/GPL570" \
         "${RAW}/GPL6883" "${LINCS_DIR}"

fetch() {  # fetch <url> <out>  —— 断点续传，已存在则跳过
  local url="$1" out="$2"
  if [[ -s "$out" ]]; then echo "  = 已存在 $(basename "$out")"; return 0; fi
  echo "  ↓ $(basename "$out")"
  curl -fSL --retry 3 --retry-delay 2 -C - -o "$out" "$url"
}

gunzip_to() {  # gunzip_to <gz> <out> —— 解压后删源，保留明文供脚本直接读取
  local gz="$1" out="$2"
  [[ -s "$out" ]] && { echo "  = 已存在 $(basename "$out")"; return 0; }
  echo "  解压 -> $(basename "$out")"
  gunzip -c "$gz" > "$out"
}

echo "== 1. GEO 表达矩阵（三个队列） =="
fetch "${GEO}/series/GSE16nnn/GSE16561/matrix/GSE16561_series_matrix.txt.gz" \
      "${RAW}/GSE16561/GSE16561_series_matrix.txt.gz"
fetch "${GEO}/series/GSE58nnn/GSE58294/matrix/GSE58294_series_matrix.txt.gz" \
      "${RAW}/GSE58294/GSE58294_series_matrix.txt.gz"
fetch "${GEO}/series/GSE22nnn/GSE22255/matrix/GSE22255_series_matrix.txt.gz" \
      "${RAW}/GSE22255/GSE22255_series_matrix.txt.gz"
for g in GSE16561 GSE58294 GSE22255; do
  gunzip_to "${RAW}/${g}/${g}_series_matrix.txt.gz" "${RAW}/${g}/${g}_series_matrix.txt"
done

echo "== 2. 平台注释 =="
fetch "${GEO}/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz" "${RAW}/GPL570/GPL570.annot.gz"
# 03_preprocess_qc.R 读的是 series 内命名（GPL6883-11606.txt）；此处用 GEO AnnotGPL，
# 两者同为 27 行注释 + "ID/Gene symbol" 表头，脚本已做列名兼容。
fetch "${GEO}/platforms/GPL6nnn/GPL6883/annot/GPL6883.annot.gz" "${RAW}/GPL6883/GPL6883.annot.gz"
gunzip_to "${RAW}/GPL6883/GPL6883.annot.gz" "${RAW}/GSE16561/GPL6883-11606.txt"

echo "== 3. L1000 / LINCS（GSE70138，CMAP02 化学扰动） =="
for f in GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz \
         GSE70138_Broad_LINCS_sig_metrics_2017-03-06.txt.gz \
         GSE70138_Broad_LINCS_pert_info_2017-03-06.txt.gz \
         GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz \
         GSE70138_Broad_LINCS_cell_info_2017-04-28.txt.gz; do
  fetch "${GEO}/series/GSE70nnn/GSE70138/suppl/${f}" "${LINCS_DIR}/${f}"
done

if [[ "${1:-}" == "--with-lincs-level5" ]]; then
  echo "  （Level5 约 5 GB，解压后约 6 GB）"
  f="GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz"
  fetch "${GEO}/series/GSE70nnn/GSE70138/suppl/${f}" "${LINCS_DIR}/${f}"
  gunzip_to "${LINCS_DIR}/${f}" "${LINCS_DIR}/${f%.gz}"
else
  echo "  - 跳过 Level5 矩阵（需要时加 --with-lincs-level5）"
fi

cat <<'EOF'

== 完成 ==
下一步：
  Rscript scripts/00_install_deps.R     # 安装 R 依赖
  Rscript scripts/01_parse_series_matrix.R

说明：
  * 978 个 L1000 landmark 清单、Reactome 通路集合已随仓库放在 reference/，
    不需要额外下载（Reactome 版本固定，换版本会改变通路检验结果）。
  * 药物开发状态由 ChEMBL 查询得到，结果缓存在 metadata/chembl_*.json，已入库。
EOF
