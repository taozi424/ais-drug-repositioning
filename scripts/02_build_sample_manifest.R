#!/usr/bin/env Rscript
# 02_build_sample_manifest.R
# 目的：构建跨队列统一样本元数据表。每个样本必录 11 个字段，缺失一律写 unknown，
#       并在 evidence 列留下"这个字段的值从哪来"的痕迹。
# 关键约束：GSE58294 同一患者 3 个时点必须共享 subject_id，且标注为重复测量，
#           不得当作三个独立队列。

suppressPackageStartupMessages({ library(data.table) })

# 项目根目录：环境变量 AIS_ROOT 优先，否则取本脚本所在目录的上一级
root <- Sys.getenv("AIS_ROOT", unset = local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) normalizePath(file.path(dirname(f[1]), "..")) else normalizePath(getwd())
}))
parsed <- file.path(root, "data", "parsed")
meta_dir <- file.path(root, "metadata")
dir.create(meta_dir, recursive = TRUE, showWarnings = FALSE)

UNK <- "unknown"

# ---------------- GSE16561：发现队列 ----------------
a <- fread(file.path(parsed, "GSE16561_sample_annot.csv"))
g1 <- data.table(
  sample_id        = a$sample_id,
  subject_id       = a$sample_id,          # 每例仅一次采样，样本即受试者；无重复测量
  cohort           = "GSE16561",
  disease_group    = fifelse(a$description == "Stroke", "AIS", "Control"),
  stroke_subtype   = UNK,
  onset_to_draw_h  = UNK,
  pretreatment     = UNK,
  sex              = fifelse(tolower(a$gender) == "female", "F",
                      fifelse(tolower(a$gender) == "male", "M", UNK)),
  age              = suppressWarnings(as.numeric(a$age)),
  platform         = a$platform,
  batch            = "GPL6883_single_series",
  timepoint_label  = "single",
  tissue           = "whole blood (PAXgene)",
  specimen         = "whole blood",
  analysis_role    = "discovery",
  evidence = paste(
    "group=Sample_description; sex/age=characteristics_ch1;",
    "subtype/onset/treatment 未在GEO字段中出现,标unknown")
)

# ---------------- GSE58294：外部重复性评估 ----------------
b <- fread(file.path(parsed, "GSE58294_sample_annot.csv"))
tp <- b$time_after_stroke_h_
b_is_case <- b$group != "Control"
g2 <- data.table(
  sample_id       = b$sample_id,
  subject_id      = b$subject_id,
  cohort          = "GSE58294",
  disease_group   = fifelse(b_is_case, "AIS", "Control"),
  stroke_subtype  = fifelse(b_is_case, "cardioembolic (CES)", "not applicable"),
  onset_to_draw_h = fifelse(!b_is_case, "not applicable",
                     fifelse(tp == "3", "<3", tp)),
  # 治疗信息不在 GEO 字段中：按原始研究设计（首采在治疗前）赋值，并显式标注为设计假设
  pretreatment    = fifelse(!b_is_case, "not applicable",
                     fifelse(tp == "3",
                             "no (pre-treatment, design assumption)",
                             "yes (post-treatment, design assumption)")),
  sex             = UNK,   # 该系列未提供性别/年龄字段
  age             = NA_real_,
  platform        = b$platform,
  batch           = "GPL570_single_series",
  timepoint_label = fifelse(!b_is_case, "control",
                     fifelse(tp == "3", "<=3h", paste0(tp, "h"))),
  tissue          = "whole blood",
  specimen        = "whole blood",
  analysis_role   = fifelse(!b_is_case, "external_replication_primary",
                     fifelse(tp == "3", "external_replication_primary",
                             "sensitivity_repeat_measure")),
  evidence = paste0(
    "group/time/subject=characteristics_ch1; ",
    "sex/age 该系列未提交,标unknown; treatment 由设计推定(GEO无字段)")
)

# ---------------- GSE22255：慢性期补充（不入主分析） ----------------
cc <- fread(file.path(parsed, "GSE22255_sample_annot.csv"))
is_case <- grepl("IS patient", cc$affected_status_disease_state_, ignore.case = TRUE)
g3 <- data.table(
  sample_id       = cc$sample_id,
  subject_id      = cc$sample_id,
  cohort          = "GSE22255",
  disease_group   = fifelse(is_case, "AIS", "Control"),
  stroke_subtype  = UNK,
  # 原始研究为卒中后 >=6 个月采样；GEO 字段无时间，仅保留定性说明
  onset_to_draw_h = UNK,
  pretreatment    = UNK,
  sex             = fifelse(tolower(cc$gender) == "female", "F",
                      fifelse(tolower(cc$gender) == "male", "M", UNK)),
  age             = suppressWarnings(as.numeric(cc$age_at_examination)),
  platform        = cc$platform,
  batch           = sub("^group (.) of hybridized samples$", "hyb_group_\\1", cc$description),
  timepoint_label = "chronic_single",
  tissue          = "PBMC",
  specimen        = "PBMC",
  analysis_role   = "chronic_supplementary",
  evidence = paste(
    "group/sex/age=characteristics_ch1; batch=Sample_description(杂交批次);",
    "采样时点GEO未给,文献为>=6个月,标unknown")
)

manifest <- rbindlist(list(g1, g2, g3), fill = TRUE)
setcolorder(manifest, c(
  "sample_id", "subject_id", "cohort", "disease_group", "stroke_subtype",
  "onset_to_draw_h", "pretreatment", "sex", "age", "platform", "batch",
  "timepoint_label", "tissue", "specimen", "analysis_role", "evidence"))

# 完整性自检
stopifnot(!anyDuplicated(manifest$sample_id))
stopifnot(all(!is.na(manifest$disease_group)))
# 同一 subject 的样本必须落在同一队列同一分组
chk <- manifest[, .(n_cohort = uniqueN(cohort), n_group = uniqueN(disease_group)),
                by = subject_id]
stopifnot(all(chk$n_cohort == 1), all(chk$n_group == 1))

fwrite(manifest, file.path(meta_dir, "sample_manifest.csv"))

# 字段字典（写清 unknown 的含义与判定规则）
dict <- data.table(
  field = c("sample_id", "subject_id", "cohort", "disease_group", "stroke_subtype",
            "onset_to_draw_h", "pretreatment", "sex", "age", "platform", "batch",
            "timepoint_label", "tissue", "specimen", "analysis_role", "evidence"),
  definition = c(
    "GEO 样本登录号 (GSM*)",
    "受试者编号；同一受试者多次采样共享同一编号",
    "GEO 系列号",
    "AIS = 急性缺血性卒中；Control = 对照",
    "卒中亚型；未报告则为 unknown",
    "发病至首次采血时间（小时）；<3 表示不足 3 小时；unknown 表示未报告",
    "采血时是否已接受再灌注/抗栓治疗",
    "F/M/unknown",
    "采血时年龄（岁）；未报告为 NA",
    "表达谱平台",
    "批次：系列级为测序/芯片批次；GSE22255 为杂交批次 A-F",
    "采样时点标签",
    "标本类型",
    "标本类别（全血 / PBMC）",
    "本样本在本研究中的用途",
    "每个字段取值的来源与推定说明")
)
fwrite(dict, file.path(meta_dir, "sample_manifest_dictionary.csv"))

# 队列角色表
roles <- manifest[, .(
  n_samples = .N,
  n_subjects = uniqueN(subject_id),
  n_AIS = sum(disease_group == "AIS"),
  n_ctrl = sum(disease_group == "Control"),
  specimen = unique(specimen),
  roles = paste(sort(unique(analysis_role)), collapse = "; ")
), by = cohort][order(cohort)]
fwrite(roles, file.path(meta_dir, "cohort_composition.csv"))

cat("样本表行数:", nrow(manifest), "\n")
print(manifest[, .N, by = .(cohort, disease_group, timepoint_label, analysis_role)][order(cohort, timepoint_label)])
cat("\n队列构成:\n"); print(roles)
