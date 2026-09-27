# ==============================================================================
# 境外输入序列补充完整分析.R
#
# 分析目的：
#   在原始广州市基孔肯雅热患者数据的基础上，合并境外输入病例的流行病学信息，
#   并补充3例境外输入病例（张璐、JOEL GA ELUDO、PEREZ GONZALEA）的基因组序列，
#   利用 outbreaker2 贝叶斯框架对所有有序列病例进行传播链推断，
#   最终输出传播关系、代际间隔分布、传播网络图及时间序列图等结果。
#
# 数据说明：
#   - 原始患者数据：广州市本地病例 + 既有境外输入病例
#   - 境外序列补充：从 imported.fas 中提取3例境外输入病例的CHIKV基因组序列
#   - 仅纳入有基因序列的病例参与传播链推断
# ==============================================================================

# ============ 0. 加载依赖包 ============
library(outbreaker2)
library(ape)
library(readxl)
library(dplyr)
library(openxlsx)
library(ggplot2)
library(igraph)


# ==============================================================================
# Part 1：数据合并
# 读取本地患者数据与境外输入病例数据，合并为统一的患者总表；
# 同时读取原始序列采样时间表，为后续序列-患者匹配做准备。
# ==============================================================================
cat("\n========================================\n")
cat("Part 1: 数据合并\n")
cat("========================================\n")

# 定义数据根目录及各子目录路径
data_path      <- "D:/Privacy/广州市基孔肯雅热数据/Data"
data_1204_path <- file.path(data_path, "1204")
data_1230_path <- file.path(data_path, "1230")

# 1.1 读取广州市原始患者数据（本地病例 + 既有输入病例）
patient_original <- read_excel(
  file.path(data_1204_path, "GZCHIKF数据-2025年截止11月30日V1.xlsx")
)
cat("原始患者数据:", nrow(patient_original), "例\n")

# 1.2 读取境外输入病例流行病学数据
patient_import <- read_excel(
  file.path(data_1204_path, "广州市CHIKF境外输入数据-2025年.xlsx")
)
cat("境外输入数据:", nrow(patient_import), "例\n")

# 1.3 统一提取两个数据集的核心字段并纵向合并
#     保留字段：患者姓名、发病日期、病例类型
common_cols <- c("患者姓名", "发病日期", "病例类型")

patient_original_subset <- patient_original %>% select(all_of(common_cols))
patient_import_subset   <- patient_import   %>% select(all_of(common_cols))

patient_combined <- bind_rows(patient_original_subset, patient_import_subset)
cat("合并后患者总数:", nrow(patient_combined), "例\n")

# 1.4 读取已有序列的采样时间信息表（seq、name、date三列）
sampling <- read_excel(file.path(data_1204_path, "采样时间.xlsx"))
cat("原始序列采样数据:", nrow(sampling), "条\n")


# ==============================================================================
# 1.5：合并基因组序列
# 从 imported.fas 中提取3例境外输入病例的CHIKV基因组序列，
# 与原有 GZCDC_CHIKV.fasta 合并，构建完整序列集合；
# 同时为新增病例创建对应的采样记录，补充入采样时间表。
# ==============================================================================
cat("\n========================================\n")
cat("Part 1.5: 合并基因组序列\n")
cat("========================================\n")

# 1.5.1 读取原有CHIKV基因组序列（FASTA格式）
sequences_original <- read.dna(
  file.path(data_path, "GZCDC_CHIKV.fasta"),
  format = "fasta"
)
cat("原始FASTA序列数:", nrow(sequences_original), "\n")

# 1.5.2 读取境外输入补充序列文件（包含多条境外来源序列）
sequences_imported <- read.dna(
  file.path(data_1230_path, "imported.fas"),
  format = "fasta"
)
cat("imported.fas 序列数:", nrow(sequences_imported), "\n")
cat("imported.fas 序列ID列表:\n")
print(labels(sequences_imported))

# 1.5.3 定义本次需要纳入分析的3例境外输入病例及其对应序列ID
#        这3例病例在既往分析中没有基因组数据，本次从 imported.fas 中补充
import_cases <- data.frame(
  name   = c("张璐", "JOEL GA ELUDO", "PEREZ GONZALEA"),
  seq_id = c("25GZ15250", "25GZ20965", "25GZ22074"),
  stringsAsFactors = FALSE
)

cat("\n需要补充序列的3例境外输入病例:\n")
print(import_cases)

# 1.5.4 验证目标序列在 imported.fas 中是否存在
imported_labels <- labels(sequences_imported)
for (i in 1:nrow(import_cases)) {
  found <- import_cases$seq_id[i] %in% imported_labels
  cat(sprintf("  %s (%s): %s\n",
              import_cases$name[i], import_cases$seq_id[i],
              ifelse(found, "✓ 找到", "✗ 未找到")))
}

# 1.5.5 从 imported.fas 中定位并提取3条目标序列
target_seq_ids    <- import_cases$seq_id
target_idx        <- match(target_seq_ids, imported_labels)
target_idx_valid  <- target_idx[!is.na(target_idx)]

if (length(target_idx_valid) != length(target_seq_ids)) {
  warning("部分序列在 imported.fas 中未找到！")
  cat("未找到的序列:", target_seq_ids[is.na(target_idx)], "\n")
}

sequences_to_add <- sequences_imported[target_idx_valid, ]
cat("\n成功提取的目标序列数:", nrow(sequences_to_add), "\n")

# 1.5.6 检查新旧序列长度是否一致，若不一致则裁剪末尾多余列
#        依据：25GZ15130_CHIKV 同时存在于两个文件，对比后发现
#        imported.fas 序列末尾比 GZCDC_CHIKV.fasta 多出 77 个空缺字符 (-)，
#        通过截取前 ncol(sequences_original) 列使两者长度对齐。
cat("原始序列长度:", ncol(sequences_original), "bp\n")
cat("新增序列长度:", ncol(sequences_to_add), "bp\n")

if (ncol(sequences_to_add) != ncol(sequences_original)) {
  len_diff <- ncol(sequences_to_add) - ncol(sequences_original)
  if (len_diff > 0) {
    cat(sprintf(
      "⚠ 新增序列比原始序列长 %d bp（末尾空缺字符），自动裁剪至 %d bp\n",
      len_diff, ncol(sequences_original)
    ))
    sequences_to_add <- sequences_to_add[, 1:ncol(sequences_original)]
    cat("裁剪后新增序列长度:", ncol(sequences_to_add), "bp\n")
  } else {
    warning(sprintf(
      "新增序列比原始序列短 %d bp，请检查序列文件是否正确！", abs(len_diff)
    ))
  }
}

# 1.5.7 合并序列：去除已存在于原始FASTA中的重复序列后再追加
original_labels <- labels(sequences_original)
already_exist   <- target_seq_ids[target_seq_ids %in% original_labels]

if (length(already_exist) > 0) {
  # 若部分目标序列已存在原始FASTA中，则只追加尚未存在的新序列
  cat("\n以下序列已存在原始FASTA（跳过）:", paste(already_exist, collapse = ", "), "\n")
  new_idx <- which(!labels(sequences_to_add) %in% original_labels)
  if (length(new_idx) > 0) {
    sequences_combined <- rbind(sequences_original, sequences_to_add[new_idx, ])
  } else {
    sequences_combined <- sequences_original
  }
} else {
  # 全部3条序列均为新增，直接追加
  sequences_combined <- rbind(sequences_original, sequences_to_add)
}

cat("合并后总序列数:", nrow(sequences_combined), "\n")
cat("合并后序列长度:", ncol(sequences_combined), "bp\n")

# 1.5.8 为3例新增病例构建采样记录
#        由于这3例病例没有独立的采样时间，使用发病日期作为采样日期的近似值
patient_combined_temp <- patient_combined %>%
  mutate(发病日期 = as.Date(发病日期))

new_sampling_rows <- list()
for (i in 1:nrow(import_cases)) {
  patient_match <- patient_combined_temp %>%
    filter(trimws(患者姓名) == trimws(import_cases$name[i]))

  if (nrow(patient_match) > 0) {
    new_row <- data.frame(
      seq  = import_cases$seq_id[i],
      name = import_cases$name[i],
      date = patient_match$发病日期[1],
      stringsAsFactors = FALSE
    )
    new_sampling_rows[[i]] <- new_row
    cat(sprintf("  %s：发病日期 %s，序列 %s\n",
                import_cases$name[i],
                as.character(patient_match$发病日期[1]),
                import_cases$seq_id[i]))
  } else {
    cat(sprintf("  ⚠ 患者 %s 在合并患者数据中未找到\n", import_cases$name[i]))
  }
}

# 将新增采样记录并入采样时间表
new_sampling_df  <- bind_rows(new_sampling_rows)
sampling_combined <- bind_rows(sampling, new_sampling_df)

cat("\n新增采样记录:\n")
print(new_sampling_df)
cat("\n合并后采样数据总数:", nrow(sampling_combined), "条\n")


# ==============================================================================
# Part 1.6：序列-患者匹配
# 由于存在一个患者多次采样或者同名患者的情况，因此需要进行匹配
# 将采样时间表与患者总表关联，按照"发病日期 ≤ 采样日期"的时序约束，
# 对每条序列找到时间差最小的对应患者，构建序列-患者一一对应关系。
# ==============================================================================
cat("\n========================================\n")
cat("Part 1.6: 序列-患者匹配\n")
cat("========================================\n")

# 为合并后的患者数据分配唯一数字ID（后续传播链推断的基础标识）
patient_combined <- patient_combined %>%
  mutate(patient_uid = row_number())

# 执行序列-患者时序匹配：
#   - 允许一对多关联（同名患者/同一患者多条序列）
#   - 滤除"采样日期早于发病日期"的非合理记录
#   - 每条序列保留时间差最小的匹配患者
#   - 最终保证最后分析不存在重复的患者或序列
result <- sampling_combined %>%
  left_join(
    patient_combined %>% select(patient_uid, 患者姓名, 发病日期, 病例类型),
    by = c("name" = "患者姓名"),
    relationship = "many-to-many"
  ) %>%
  filter(发病日期 <= date | is.na(发病日期)) %>%
  mutate(time_diff = as.integer(date - 发病日期)) %>%
  group_by(seq) %>%
  slice_min(time_diff, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  filter(!is.na(发病日期)) %>%
  select(seq, name, date, 发病日期, 病例类型, patient_uid, time_diff)

cat("\n匹配结果:\n")
cat("匹配到的序列数:", n_distinct(result$patient_uid), "\n")

# 输出匹配患者的病例类型构成
cat("\n匹配患者的病例类型分布:\n")
print(table(result$病例类型))

# 逐一确认3例新增病例是否匹配成功
cat("\n3例新增境外输入病例匹配情况:\n")
for (i in 1:nrow(import_cases)) {
  matched <- result %>% filter(seq == import_cases$seq_id[i])
  if (nrow(matched) > 0) {
    cat(sprintf("  ✓ %s (%s) → patient_uid=%d，发病日期=%s\n",
                import_cases$name[i], import_cases$seq_id[i],
                matched$patient_uid[1], as.character(matched$发病日期[1])))
  } else {
    cat(sprintf("  ✗ %s (%s) → 未匹配\n",
                import_cases$name[i], import_cases$seq_id[i]))
  }
}

# ==============================================================================
# Part 2：传播链推断数据准备
# 基于合并后的患者总表和序列匹配结果，构建 outbreaker2 所需的：
#   - 流行病学数据框（发病天数、是否输入病例等）
#   - DNA序列矩阵（行=病例，列=碱基位点）
# 仅纳入有基因序列的病例参与推断。
# ==============================================================================
cat("\n========================================\n")
cat("Part 2: 传播链推断数据准备\n")
cat("========================================\n")

# 赋值：使用合并后的患者总表及匹配结果、合并后的序列集合
all_patients <- patient_combined
genome_epi   <- result
sequences    <- sequences_combined

# 2.1 整理流行病学数据框
#     - 统一日期格式
#     - 标记输入病例（输入/境外输入均视为输入来源）
#     - 按发病日期升序排列，计算相对发病天数（以最早发病日为第0天）
epi_data_all <- all_patients %>%
  mutate(
    case_id     = row_number(),
    发病日期    = as.Date(发病日期),
    is_imported = (病例类型 == "输入" | 病例类型 == "境外输入")
  ) %>%
  arrange(发病日期) %>%
  mutate(
    case_id    = row_number(),
    onset_days = as.integer(发病日期 - min(发病日期, na.rm = TRUE))
  )

# 2.2 分配患者唯一数字ID（与采样匹配时一致）
epi_data_all <- epi_data_all %>%
  mutate(patient_uid = row_number())

# 2.3 检查同名患者情况（同名可能导致错误匹配，需在此确认规模）
cat("\n=== 检查同名患者情况 ===\n")

dup_names_epi <- epi_data_all %>%
  group_by(患者姓名) %>% filter(n() > 1) %>% ungroup()
cat("患者总表中同名患者数:", nrow(dup_names_epi), "\n")

dup_names_genome <- genome_epi %>%
  group_by(name) %>% filter(n() > 1) %>% ungroup()
cat("基因组匹配表中多序列患者数:", nrow(dup_names_genome), "\n")

# 2.4 将序列与患者匹配结果回写至流行病学数据框
genome_for_match <- genome_epi %>%
  mutate(sample_date = as.Date(date)) %>%
  select(seq, name, sample_date)

epi_for_match <- epi_data_all %>%
  select(patient_uid, name = 患者姓名, onset_date = 发病日期)

# 二次匹配：在全量患者数据上再次执行序列-患者时序关联
# 若同一患者有多条序列，保留采样日期最早的一条
match_result <- genome_for_match %>%
  inner_join(epi_for_match, by = "name", relationship = "many-to-many") %>%
  filter(onset_date <= sample_date) %>%
  mutate(time_diff = as.integer(sample_date - onset_date)) %>%
  group_by(seq) %>%
  slice_min(time_diff, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  group_by(patient_uid) %>%
  slice_min(sample_date, n = 1, with_ties = FALSE) %>%
  ungroup()

cat("\n最终序列-患者匹配数:", nrow(match_result), "个患者有序列\n")

# 将序列名（seq）和采样日期（sample_date）合并进流行病学数据框
epi_data_all <- epi_data_all %>%
  left_join(
    match_result %>% select(patient_uid, seq, sample_date),
    by = "patient_uid"
  )

# 2.5 筛选：仅保留有基因序列的病例纳入传播链推断
epi_data <- epi_data_all %>%
  filter(!is.na(seq)) %>%
  arrange(发病日期) %>%
  mutate(
    case_id    = row_number(),   # 在有序列病例子集中重新编号（1..n_cases）
    onset_days = as.integer(发病日期 - min(发病日期, na.rm = TRUE))
  )

n_cases <- nrow(epi_data)
cat("\n========================================\n")
cat("纳入传播链推断的病例（仅有基因序列）\n")
cat("========================================\n")
cat("纳入病例总数:", n_cases, "\n")
cat("  输入病例:", sum(epi_data$is_imported), "\n")
cat("  本地病例:", sum(!epi_data$is_imported), "\n")

# 确认3例新增境外输入病例已成功纳入
cat("\n3例新增境外输入病例纳入情况:\n")
for (i in 1:nrow(import_cases)) {
  in_data <- epi_data %>% filter(seq == import_cases$seq_id[i])
  if (nrow(in_data) > 0) {
    cat(sprintf("  ✓ %s (case_id=%d, onset_days=%d, imported=%s)\n",
                import_cases$name[i], in_data$case_id[1],
                in_data$onset_days[1], in_data$is_imported[1]))
  } else {
    cat(sprintf("  ✗ %s 未纳入\n", import_cases$name[i]))
  }
}

# 2.6 构建DNA序列矩阵（n_cases 行 × n_sites 列，元素为 raw 类型碱基编码）
#     对无法匹配序列的病例行，以 0xf0（N/缺失）填充
seq_labels <- labels(sequences)
n_sites    <- ncol(sequences)

dna_matrix <- matrix(as.raw(0xf0), nrow = n_cases, ncol = n_sites)

matched_count <- 0
for (i in 1:n_cases) {
  seq_idx <- match(epi_data$seq[i], seq_labels)
  if (!is.na(seq_idx)) {
    dna_matrix[i, ] <- as.raw(sequences[seq_idx, ])
    matched_count   <- matched_count + 1
  } else {
    cat(sprintf("  ⚠ 序列 %s 在合并FASTA中未找到\n", epi_data$seq[i]))
  }
}

class(dna_matrix)       <- "DNAbin"
rownames(dna_matrix)    <- as.character(1:n_cases)

cat("\nDNA矩阵维度:", nrow(dna_matrix), "×", ncol(dna_matrix), "\n")
cat("成功匹配序列数:", matched_count, "/", n_cases, "\n")


# ==============================================================================
# Part 3：分布参数计算（感染自然史参数）
# 基于文献的基孔肯雅热自然史参数，计算 outbreaker2 所需的：
#   - Serial Interval（代际间隔）离散密度向量 w
#   - 潜伏期（Incubation Period）离散密度向量 f
# 两者均假设服从 Gamma 分布，由均值和范围反推形状参数与尺度参数。
# ==============================================================================
# Serial Interval：均值 23天，参考范围 17–38天
si_mean        <- 23
si_range_width <- 38 - 17
si_sd          <- si_range_width / 4   # 95%区间 ≈ ±2SD
si_var         <- si_sd^2
si_shape       <- si_mean^2 / si_var
si_scale       <- si_mean / si_shape

w  <- dgamma(1:50, shape = si_shape, scale = si_scale)
w  <- w / sum(w)

# 潜伏期：均值 3天，参考范围 0.5-3.1天
inc_mean        <- 3
inc_range_width <- 3.1 - 0.5
inc_sd          <- inc_range_width / 4
inc_var         <- inc_sd^2
inc_shape       <- inc_mean^2 / inc_var
inc_scale       <- inc_mean / inc_shape

f  <- dgamma(1:30, shape = inc_shape, scale = inc_scale)
f  <- f / sum(f)

cat("\n=== 感染自然史分布参数（基孔肯雅热） ===\n")
cat("Serial Interval : shape =", round(si_shape, 3), "，scale =", round(si_scale, 3),
    "，均值 =", round(si_shape * si_scale, 2), "天\n")
cat("潜伏期          : shape =", round(inc_shape, 3), "，scale =", round(inc_scale, 3),
    "，均值 =", round(inc_shape * inc_scale, 2), "天\n")


# ==============================================================================
# Part 4：初始化传播链祖先（init_alpha）
# 对每个本地病例，以发病时间最接近且早于自身的病例作为初始祖先；
# 输入病例（境外来源）的祖先设置为 NA，表示其传播来源在观测队列之外。
# ==============================================================================
init_alpha <- rep(NA_integer_, n_cases)

for (i in 1:n_cases) {
  if (!epi_data$is_imported[i]) {
    # 候选祖先：发病时间早于当前病例的所有病例
    candidates <- which(epi_data$onset_days < epi_data$onset_days[i])
    if (length(candidates) > 0) {
      time_diffs     <- epi_data$onset_days[i] - epi_data$onset_days[candidates]
      init_alpha[i]  <- candidates[which.min(time_diffs)]
    }
  }
}

cat("\n初始化祖先设置:\n")
cat("  输入病例（祖先=NA）:", sum(is.na(init_alpha) & epi_data$is_imported), "\n")
cat("  本地病例（有初始祖先）:", sum(!is.na(init_alpha)), "\n")


# ==============================================================================
# Part 5：配置并运行 outbreaker2
# outbreaker2 采用贝叶斯 MCMC 方法推断传播链。
# pi（观测到序列的病例比例）固定为"有序列病例数/患者总数"，不允许 MCMC 更新，
# 以避免因样本量不足导致的 pi 估计偏差。
# ==============================================================================

# 构建 outbreaker2 输入数据对象
data_ob2 <- outbreaker_data(
  dates  = epi_data$onset_days,
  dna    = dna_matrix,
  w_dens = w,
  f_dens = f
)

# 计算 pi 先验值：有基因序列的病例数 / 患者总数
n_seq_cases   <- nrow(genome_epi)
n_total_cases <- nrow(patient_combined)
pi_prior      <- n_cases / n_total_cases

cat("\n=== pi 先验设置（序列覆盖率） ===\n")
cat("pi 初始值:", round(pi_prior, 4),
    "(", n_cases, "有序列 /", n_total_cases, "总病例)\n")
cat("pi move: 禁用（固定不更新）\n")

# 配置 MCMC 运行参数
config <- create_config(
  n_iter      = 10000,   # MCMC 总迭代次数
  sample_every = 50,     # 每隔50步采样一次（共200个后验样本）
  init_alpha  = init_alpha,
  init_pi     = pi_prior,
  move_pi     = FALSE,   # 固定 pi，不允许 MCMC 更新
  pb          = TRUE     # 显示进度条
)

# 运行 outbreaker2 传播链推断
cat("\n开始运行 outbreaker2...\n")
results <- outbreaker(data = data_ob2, config = config)


# ==============================================================================
# Part 6：模型拟合诊断
# 通过链迹图（trace plot）检查 MCMC 是否收敛；
# 输出后验分布和关键参数（mu：突变率；pi：序列覆盖率）的统计摘要。
# ==============================================================================
plot(results)           # 后验参数链迹图
plot(results, "prior")  # 先验分布
plot(results, "mu")     # 突变率后验
plot(results, "pi")     # 序列覆盖率后验
summary(results)

# 提取 MCMC 后验样本为数据框，用于后续统计
res_df <- as.data.frame(results)

# 突变率（mu）后验统计
mu_samples <- res_df$mu
cat("\n=== mu（每位点每代突变率）后验统计 ===\n")
cat("均值   :", round(mean(mu_samples), 10), "\n")
cat("中位数 :", round(median(mu_samples), 10), "\n")
cat("95% CI : [",
    round(quantile(mu_samples, 0.025), 10), ",",
    round(quantile(mu_samples, 0.975), 6), "]\n")


# ==============================================================================
# Part 7：传播链 Support 统计
# 对每个病例的传播关系，计算后验样本中最高频祖先的出现概率（support），
# 反映传播链推断的可信度。
# ==============================================================================
alpha_cols    <- grep("^alpha_", names(res_df), value = TRUE)
n_samples     <- nrow(res_df)

support_values <- sapply(alpha_cols, function(col) {
  alphas <- res_df[[col]]
  if (all(is.na(alphas))) return(NA)
  tab <- table(alphas, useNA = "no")
  if (length(tab) == 0) return(NA)
  max(tab) / n_samples
})

support_clean <- support_values[!is.na(support_values)]
cat("\n传播链 support（最高频祖先后验概率）:\n")
cat("  均值   =", round(mean(support_clean), 3), "\n")
cat("  95%CI  = [",
    round(quantile(support_clean, 0.025), 3), ",",
    round(quantile(support_clean, 0.975), 3), "]\n")


# ==============================================================================
# Part 8：提取最可能传播关系
# 对每个病例，取后验样本中出现频率最高的祖先作为"最可能传播来源"，
# 并记录对应的后验概率；同时计算传播对之间的代际间隔（天）。
# ==============================================================================
transmission_df <- data.frame(
  case_id      = 1:n_cases,
  onset_date   = epi_data$发病日期,
  is_imported  = epi_data$is_imported,
  name         = epi_data$患者姓名,
  ancestor     = NA_integer_,
  ancestor_prob = NA_real_
)

for (i in 1:n_cases) {
  col_name <- paste0("alpha_", i)
  if (col_name %in% names(res_df)) {
    alphas <- res_df[[col_name]]
    if (!all(is.na(alphas))) {
      tab <- table(alphas, useNA = "no")
      if (length(tab) > 0) {
        transmission_df$ancestor[i]      <- as.integer(names(tab)[which.max(tab)])
        transmission_df$ancestor_prob[i] <- max(tab) / n_samples
      }
    }
  }
}

# 将祖先的发病日期合并进来，计算代际间隔
transmission_df <- transmission_df %>%
  left_join(
    transmission_df %>% select(case_id, ancestor_onset = onset_date),
    by = c("ancestor" = "case_id")
  ) %>%
  mutate(serial_interval = as.integer(onset_date - ancestor_onset))


# ==============================================================================
# Part 9：地址匹配与传播链筛选
# 读取患者居住街道/镇信息（包括原始及境外补充两份地址数据），
# 筛选满足以下条件的传播对作为"高可信度本地传播"：
#   (1) 传播者与被传播者居住于同一街镇
#   (2) 传播源后验概率 ≥ 0.95
# ==============================================================================

# 9.1 读取并合并两份地址信息文件
address_info_original <- read_excel(
  file.path(data_path, "患者地址信息提取结果.xlsx")
) %>% mutate(姓名_clean = trimws(姓名))

address_info_import <- read_excel(
  file.path(data_path, "患者信息提取结果境外补充.xlsx")
) %>% mutate(姓名_clean = trimws(姓名))

address_info <- bind_rows(address_info_original, address_info_import)
cat("\n地址信息合并:\n")
cat("  原始地址信息:", nrow(address_info_original), "条\n")
cat("  境外补充地址信息:", nrow(address_info_import), "条\n")
cat("  合并后总数:", nrow(address_info), "条\n")

# 9.2 将被传播者和传播者的街道/镇信息分别关联进传播关系表
transmission_with_address <- transmission_df %>%
  mutate(name_clean = trimws(name)) %>%
  left_join(
    address_info %>% select(姓名_clean, case_street = `街道/镇`),
    by = c("name_clean" = "姓名_clean")
  ) %>%
  left_join(
    transmission_df %>%
      mutate(name_clean = trimws(name)) %>%
      left_join(address_info %>% select(姓名_clean, source_street = `街道/镇`),
                by = c("name_clean" = "姓名_clean")) %>%
      select(case_id, source_street),
    by = c("ancestor" = "case_id")
  )

# 9.3 双重筛选：街镇相同 + 传播源后验概率 ≥ 0.95
same_address_transmission <- transmission_with_address %>%
  filter(
    !is.na(ancestor),
    !is.na(case_street),
    !is.na(source_street),
    case_street == source_street,
    !is.na(ancestor_prob),
    ancestor_prob >= 0.95
  ) %>%
  left_join(
    transmission_df %>% select(case_id, source_name = name),
    by = c("ancestor" = "case_id")
  ) %>%
  select(
    被传播者case_id   = case_id,
    被传播者姓名      = name,
    被传播者街道镇    = case_street,
    传播者case_id     = ancestor,
    传播者姓名        = source_name,
    传播者街道镇      = source_street,
    被传播者发病日期  = onset_date,
    传播者发病日期    = ancestor_onset,
    代际间隔_天       = serial_interval,
    传播源概率        = ancestor_prob
  ) %>%
  arrange(被传播者发病日期)

cat("\n========================================\n")
cat("街镇相同 + 高可信度(≥0.95) 传播对统计\n")
cat("========================================\n")
cat("满足条件的传播对数:", nrow(same_address_transmission), "对\n")
print(same_address_transmission)


# ==============================================================================
# Part 10：确认输出目录
# ==============================================================================
output_dir <- file.path(data_path, "../Output")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)


# ==============================================================================
# Part 11：Serial Interval（代际间隔）描述统计与 Gamma 分布拟合
# 提取筛选后传播对的代际间隔数据，计算描述统计量，
# 并以最大似然法拟合 Gamma 分布，输出形状/尺度参数及95%分位区间。
# ==============================================================================
si_data <- same_address_transmission %>%
  filter(!is.na(代际间隔_天) & 代际间隔_天 > 0) %>%
  pull(代际间隔_天)

cat("\n=== Serial Interval 描述统计 ===\n")
cat("样本量:", length(si_data), "\n")

if (length(si_data) > 0) {
  cat("均值   :", round(mean(si_data), 2), "天\n")
  cat("中位数 :", median(si_data), "天\n")
  cat("标准差 :", round(sd(si_data), 2), "天\n")
  cat("范围   :", min(si_data), "–", max(si_data), "天\n")

  # 以矩估计值为初值，用 L-BFGS-B 最优化最大对数似然，拟合 Gamma 分布
  si_mean_obs  <- mean(si_data)
  si_var_obs   <- var(si_data)
  shape_init   <- si_mean_obs^2 / si_var_obs
  scale_init   <- si_var_obs / si_mean_obs

  gamma_fit <- optim(
    c(shape_init, scale_init),
    function(p) -sum(dgamma(si_data, shape = p[1], scale = p[2], log = TRUE)),
    method = "L-BFGS-B",
    lower  = c(0.1, 0.1)
  )

  fitted_shape <- gamma_fit$par[1]
  fitted_scale <- gamma_fit$par[2]

  cat("\n=== Gamma 分布拟合结果 ===\n")
  cat("Shape 参数 (k) :", round(fitted_shape, 4), "\n")
  cat("Scale 参数 (θ) :", round(fitted_scale, 4), "\n")
  cat("均值   E[X]    :", round(fitted_shape * fitted_scale, 4), "天\n")
  cat("标准差 SD      :", round(sqrt(fitted_shape) * fitted_scale, 4), "天\n")

  ci_lower <- qgamma(0.025, fitted_shape, scale = fitted_scale)
  ci_upper <- qgamma(0.975, fitted_shape, scale = fitted_scale)
  cat("95% CI         : [", round(ci_lower, 2), ",", round(ci_upper, 2), "] 天\n")


  # ============================================================================
  # Part 12：Bootstrap 采样绘制 Serial Interval 分布图
  # 重复1000次 Bootstrap 重采样并拟合 Gamma 分布，绘制不确定性区间；
  # 以黑色粗线叠加原始拟合曲线，生成发表级别的 SI 分布图。
  # ============================================================================
  set.seed(123)
  n_boot  <- 1000
  x_range <- seq(0, 42, length.out = 500)   # x 轴范围：0–42天

  # 存储每次 Bootstrap 拟合的密度曲线（行=bootstrap次，列=x轴点）
  boot_densities <- matrix(NA, nrow = n_boot, ncol = length(x_range))
  # 同时保存每次 Bootstrap 拟合的 Gamma 参数（shape, scale），供后续 CDF 绘图使用
  boot_params    <- matrix(NA, nrow = n_boot, ncol = 2)

  cat("\n正在进行 Bootstrap 拟合（n=1000）...\n")
  for (b in 1:n_boot) {
    boot_sample <- sample(si_data, replace = TRUE)
    boot_mean   <- mean(boot_sample)
    boot_var    <- var(boot_sample)
    shape_init  <- boot_mean^2 / boot_var
    scale_init  <- boot_var / boot_mean

    tryCatch({
      boot_fit <- optim(
        c(shape_init, scale_init),
        function(p) -sum(dgamma(boot_sample, shape = p[1], scale = p[2], log = TRUE)),
        method = "L-BFGS-B", lower = c(0.1, 0.1)
      )
      boot_densities[b, ] <- dgamma(
        x_range, shape = boot_fit$par[1], scale = boot_fit$par[2]
      )
      boot_params[b, ]    <- boot_fit$par
    }, error = function(e) {})
  }
  cat("Bootstrap 完成！\n")

  # 随机抽取50条 Bootstrap 曲线绘制不确定性背景
  n_show    <- min(50, n_boot)
  boot_idx  <- sample(which(!is.na(boot_densities[, 1])), n_show)
  boot_lines_df <- do.call(rbind, lapply(boot_idx, function(i) {
    data.frame(si = x_range, density = boot_densities[i, ], boot_id = i)
  }))

  # 原始最大似然拟合的密度曲线
  pdf_data <- data.frame(
    si      = x_range,
    density = dgamma(x_range, shape = fitted_shape, scale = fitted_scale)
  )

  # 绘图：灰色 Bootstrap 曲线背景 + 黑色主曲线
  p_si <- ggplot() +
    geom_line(data = boot_lines_df,
              aes(x = si, y = density, group = boot_id),
              color = "grey50", alpha = 0.30, linewidth = 0.4) +
    geom_line(data = pdf_data,
              aes(x = si, y = density),
              color = "black", linewidth = 1.2) +
    labs(x = "Serial Interval (days)", y = "Relative Frequency") +
    scale_x_continuous(limits = c(0, 42),  breaks = seq(0, 42, by = 7), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 0.1), breaks = seq(0, 0.1, by = 0.02), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 14) +
    theme(
      panel.border       = element_rect(fill = NA, color = "black", linewidth = 0.8),
      panel.background   = element_rect(fill = "white"),
      panel.grid.major   = element_line(color = "grey85", linewidth = 0.4, linetype = "dashed"),
      panel.grid.minor   = element_blank(),
      plot.background    = element_rect(fill = "white", color = NA),
      axis.line          = element_blank(),
      axis.ticks         = element_line(color = "black", linewidth = 0.8),
      axis.ticks.length  = unit(0.3, "cm"),
      axis.text          = element_text(color = "black", size = 13, family = "sans"),
      axis.title.x       = element_text(color = "black", size = 15, face = "bold",
                                        margin = margin(t = 12)),
      axis.title.y       = element_text(color = "black", size = 15, face = "bold",
                                        margin = margin(r = 12)),
      plot.margin        = margin(25, 25, 20, 20)
    )

  print(p_si)

  # ============================================================================
  # Part 13：保存 Serial Interval 分布图
  # ============================================================================
  ggsave(file.path(output_dir, "si_import_seq_supp_3cases.pdf"),
         p_si, width = 8, height = 5.5)
  ggsave(file.path(output_dir, "si_import_seq_supp_3cases.png"),
         p_si, width = 8, height = 5.5, dpi = 300)
  cat("\n✓ Serial Interval 分布图已保存\n")


  # ============================================================================
  # Part 13.1：Bootstrap 采样绘制 Serial Interval 累积分布函数（CDF）图
  # 利用 Part 12 中保存的 Bootstrap Gamma 参数，计算每次重采样的 CDF 曲线，
  # 叠加经验 CDF 散点（红色）和最大似然拟合 CDF 曲线（黑色），全面评估拟合效果。
  # ============================================================================

  # 从已保存的 boot_params 中抽取50条有效 Bootstrap CDF 曲线
  valid_boots     <- which(!is.na(boot_params[, 1]))
  boot_cdf_idx    <- sample(valid_boots, min(50, length(valid_boots)))

  boot_cdf_lines_df <- do.call(rbind, lapply(boot_cdf_idx, function(i) {
    data.frame(
      si      = x_range,
      cdf     = pgamma(x_range, shape = boot_params[i, 1], scale = boot_params[i, 2]),
      boot_id = i
    )
  }))

  # 最大似然拟合的 CDF 曲线
  cdf_data <- data.frame(
    si  = x_range,
    cdf = pgamma(x_range, shape = fitted_shape, scale = fitted_scale)
  )

  # 经验 CDF 数据（基于观测 SI 值的阶梯函数）
  ecdf_fn   <- ecdf(si_data)
  ecdf_data <- data.frame(
    si  = sort(unique(si_data)),
    cdf = ecdf_fn(sort(unique(si_data)))
  )

  # 绘制 CDF 图：灰色 Bootstrap 曲线 + 黑色拟合 CDF + 红色经验 CDF 散点
  p_cdf <- ggplot() +
    geom_line(data = boot_cdf_lines_df,
              aes(x = si, y = cdf, group = boot_id),
              color = "grey50", alpha = 0.30, linewidth = 0.4) +
    geom_line(data = cdf_data,
              aes(x = si, y = cdf),
              color = "black", linewidth = 1.2) +
    geom_point(data = ecdf_data,
               aes(x = si, y = cdf),
               color = "#D9534F", size = 2.5) +
    labs(x = "Serial Interval (days)", y = "Cumulative Probability") +
    scale_x_continuous(limits = c(0, 42), breaks = seq(0, 42, by = 7), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1),  breaks = seq(0, 1, by = 0.2), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 14) +
    theme(
      panel.border       = element_rect(fill = NA, color = "black", linewidth = 0.8),
      panel.background   = element_rect(fill = "white"),
      panel.grid.major   = element_line(color = "grey85", linewidth = 0.4, linetype = "dashed"),
      panel.grid.minor   = element_blank(),
      plot.background    = element_rect(fill = "white", color = NA),
      axis.line          = element_blank(),
      axis.ticks         = element_line(color = "black", linewidth = 0.8),
      axis.ticks.length  = unit(0.3, "cm"),
      axis.text          = element_text(color = "black", size = 13, family = "sans"),
      axis.title.x       = element_text(color = "black", size = 15, face = "bold",
                                        margin = margin(t = 12)),
      axis.title.y       = element_text(color = "black", size = 15, face = "bold",
                                        margin = margin(r = 12)),
      plot.margin        = margin(25, 25, 20, 20)
    )

  print(p_cdf)

  # 保存 CDF 图
  ggsave(file.path(output_dir, "si_cdf_import_seq_supp_3cases.pdf"),
         p_cdf, width = 8, height = 5.5)
  ggsave(file.path(output_dir, "si_cdf_import_seq_supp_3cases.png"),
         p_cdf, width = 8, height = 5.5, dpi = 300)
  cat("\n✓ Serial Interval CDF 图已保存\n")

} else {
  cat("\n警告：没有有效的 SI 数据，跳过分布拟合与绘图\n")
}


# ==============================================================================
# Part 14：保存传播链结果至 Excel
# 输出两个版本：
#   (1) 全量传播关系表（含所有病例的最可能祖先及后验概率）
#   (2) 筛选后的高可信度同街镇传播对（供后续分析使用）
# ==============================================================================
write.xlsx(transmission_df,
           file.path(output_dir, "transmission_import_seq_supp_3cases_all.xlsx"))
write.xlsx(same_address_transmission,
           file.path(output_dir, "transmission_import_seq_supp_3cases_same_street_high_conf.xlsx"))
cat("\n✓ 传播链结果已保存\n")


# ==============================================================================
# Part 15：绘制传播链时间序列图
# 以水平线段展示每个传播对的时间跨度（传播者发病日 → 被传播者发病日），
# 节点颜色区分输入来源与本地来源，直观呈现传播的时间动态。
# 额外输出一张排除"传播者与被传播者发病日相同"传播对的版本，
# 以避免零代际间隔在图形上的叠压。
# ==============================================================================
if (nrow(same_address_transmission) > 0) {

  # 15.1 全量时间序列图（含零代际间隔传播对）
  transmission_timeline <- same_address_transmission %>%
    arrange(被传播者发病日期) %>%
    mutate(
      pair_rank  = row_number(),
      pair_label = paste0("P", row_number())
    ) %>%
    left_join(transmission_df %>% select(case_id, is_imported),
              by = c("传播者case_id" = "case_id")) %>%
    mutate(传播者类型 = ifelse(is_imported, "Imported case", "Local case"))

  p_timeline <- ggplot(transmission_timeline) +
    geom_segment(aes(x = 传播者发病日期, xend = 被传播者发病日期,
                     y = pair_rank, yend = pair_rank),
                 color = "#2E7DB8", linewidth = 0.7) +
    geom_point(aes(x = 传播者发病日期, y = pair_rank, color = 传播者类型),
               size = 2.5) +
    geom_point(aes(x = 被传播者发病日期, y = pair_rank),
               color = "#2E7DB8", size = 2.5) +
    scale_y_reverse(breaks  = transmission_timeline$pair_rank,
                    labels  = transmission_timeline$pair_label) +
    scale_x_date(date_labels = "%m-%d", date_breaks = "2 days") +
    scale_color_manual(
      values = c("Imported case" = "#D9534F", "Local case" = "#2E7DB8"),
      name   = "Source case type"
    ) +
    labs(x = "Onset Date", y = "Transmission Pair") +
    theme_classic() +
    theme(
      axis.text.x        = element_text(angle = 45, hjust = 1),
      legend.position    = c(0.99, 0.99),
      legend.justification = c(1, 1),
      legend.background  = element_rect(fill = "white", color = "black")
    )

  ggsave(file.path(output_dir, "timeline_import_seq_supp_3cases.pdf"),
         p_timeline, width = 10, height = 8)
  ggsave(file.path(output_dir, "timeline_import_seq_supp_3cases.png"),
         p_timeline, width = 10, height = 8, dpi = 300)
  cat("\n✓ 传播链时间序列图已保存\n")

  # 16.2 排除发病日相同的传播对（代际间隔 > 0 版本）
  transmission_timeline_diff <- same_address_transmission %>%
    filter(传播者发病日期 != 被传播者发病日期) %>%
    arrange(被传播者发病日期) %>%
    mutate(
      pair_rank  = row_number(),
      pair_label = paste0("P", row_number())
    ) %>%
    left_join(transmission_df %>% select(case_id, is_imported),
              by = c("传播者case_id" = "case_id")) %>%
    mutate(传播者类型 = ifelse(is_imported, "Imported case", "Local case"))

  cat("\n排除发病日相同传播对后剩余:", nrow(transmission_timeline_diff), "对\n")

  if (nrow(transmission_timeline_diff) > 0) {
    p_timeline_diff <- ggplot(transmission_timeline_diff) +
      geom_segment(aes(x = 传播者发病日期, xend = 被传播者发病日期,
                       y = pair_rank, yend = pair_rank),
                   color = "#2E7DB8", linewidth = 0.7) +
      geom_point(aes(x = 传播者发病日期, y = pair_rank, color = 传播者类型),
                 size = 2.5) +
      geom_point(aes(x = 被传播者发病日期, y = pair_rank),
                 color = "#2E7DB8", size = 2.5) +
      scale_y_reverse(breaks  = transmission_timeline_diff$pair_rank,
                      labels  = transmission_timeline_diff$pair_label) +
      scale_x_date(date_labels = "%m-%d", date_breaks = "2 days") +
      scale_color_manual(
        values = c("Imported case" = "#D9534F", "Local case" = "#2E7DB8"),
        name   = "Source case type"
      ) +
      labs(x = "Onset Date", y = "Transmission Pair") +
      theme_classic() +
      theme(
        axis.text.x          = element_text(angle = 45, hjust = 1),
        legend.position      = c(0.99, 0.99),
        legend.justification = c(1, 1),
        legend.background    = element_rect(fill = "white", color = "black")
      )

    ggsave(file.path(output_dir, "timeline_import_seq_supp_3cases_diff_onset.pdf"),
           p_timeline_diff, width = 10, height = 8)
    ggsave(file.path(output_dir, "timeline_import_seq_supp_3cases_diff_onset.png"),
           p_timeline_diff, width = 10, height = 8, dpi = 300)
    cat("\n✓ 传播链时间序列图（排除发病日相同）已保存\n")
  }

} else {
  cat("\n⚠ 没有满足条件的传播对，跳过时间序列图绘制\n")
}


# ==============================================================================
# Part 16：保存分析环境
# 将当前 R 工作环境（所有对象）保存为 .RData 文件，
# 方便后续加载复现结果，无需重跑耗时的 outbreaker2 推断。
# ==============================================================================
save.image(file = file.path(output_dir, "chikv_import_seq_supp_3cases_environment.RData"))
cat("\n✓ 分析环境已保存\n")
