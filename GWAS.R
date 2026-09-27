# =========================
# Block 2: Adjusted GLM per SNP (from consensus dna_sequence)
# Reference: NC_004162 (genome FASTA)
# =========================

library(dplyr)
library(purrr)
library(broom)
library(splines)
library(stringr)
library(Biostrings)
library(ggplot2)

# ---------- 用户需要改的路径 ----------
ref_fasta <- "/Users/junwenzhou/Desktop/1204/reference_NC_004162.fasta"

# ---------- 0) 检查参考序列 ----------
stopifnot(file.exists(ref_fasta))
ref <- readDNAStringSet(ref_fasta)
if (length(ref) != 1) stop("Reference FASTA should contain exactly 1 sequence (whole genome).")
names(ref) <- "REF_NC_004162"

# ---------- 1) 临床数据清洗（仅保留分析需要的字段与非缺失） ----------
dat0 <- clean_data %>%
  mutate(
    gender = as.factor(gender),
    case_type = as.factor(case_type),
    age = as.numeric(age),
    days_onset_to_test = as.numeric(days_onset_to_test),
    ct_value = as.numeric(ct_value)
  ) %>%
  filter(
    !is.na(seq),
    !is.na(dna_sequence),
    !is.na(ct_value),
    !is.na(age),
    !is.na(days_onset_to_test),
    !is.na(gender),
    !is.na(case_type)
  )

# 可选：如果你决定主分析排除 days < 0（你现在看起来没有负值了）
dat0 <- dat0 %>% filter(days_onset_to_test >= 0)

# 基础检查：是否都是 DNA 字符（ACGTN）
if (any(!grepl("^[ACGTN]+$", dat0$dna_sequence))) {
  stop("Some dna_sequence contain characters outside A/C/G/T/N.")
}

# ---------- 2) MAFFT 比对（参考 + 所有样本） ----------
if (Sys.which("mafft") == "") {
  stop("Cannot find `mafft` in PATH. Please install MAFFT (e.g., brew install mafft).")
}

Sys.which("mafft")

Sys.getenv("PATH")
Sys.setenv(PATH = paste("/opt/anaconda3/bin", Sys.getenv("PATH"), sep = ":"))
Sys.which("mafft")
system("mafft --version")

qry <- DNAStringSet(dat0$dna_sequence)
names(qry) <- dat0$seq

in_fa  <- tempfile(fileext = ".fasta")
out_fa <- tempfile(fileext = ".fasta")
writeXStringSet(c(ref, qry), filepath = in_fa, format = "fasta")

cmd <- paste("mafft --auto --thread -1", shQuote(in_fa), ">", shQuote(out_fa))
status <- system(cmd)
if (status != 0) stop("MAFFT failed. Check MAFFT installation and input FASTA.")

aln <- readDNAStringSet(out_fa)
if (!("REF_NC_004162" %in% names(aln))) stop("Reference not found in alignment output.")

# ---------- 3) 从比对提取 callable 位点、SNP 0/1 矩阵，并计算 L_callable ----------
aln_mat <- as.matrix(aln)  # rows: sequences, cols: alignment columns
ref_row <- which(names(aln) == "REF_NC_004162")
ref_aln <- aln_mat[ref_row, ]
sample_rows <- setdiff(seq_len(nrow(aln_mat)), ref_row)
sample_mat <- aln_mat[sample_rows, , drop = FALSE]

is_acgt <- function(x) x %in% c("A","C","G","T")

# 3.1 建立 alignment column -> reference position 的映射（参考去 gap 后的位置）
ref_pos <- rep(NA_integer_, length(ref_aln))
pos_counter <- 0L
for (k in seq_along(ref_aln)) {
  if (ref_aln[k] != "-") {
    pos_counter <- pos_counter + 1L
    ref_pos[k] <- pos_counter
  }
}

# 3.2 定义 callable 位点：
# - 参考该列是 A/C/G/T 且不是 gap
# - 样本中该列 A/C/G/T 的比例 >= 0.95（避免 N/gap 太多）
ref_callable <- is_acgt(ref_aln) & !is.na(ref_pos)
sample_acgt_prop <- colMeans(apply(sample_mat, 2, is_acgt), na.rm = TRUE)
callable_col <- ref_callable & (sample_acgt_prop >= 0.95)

callable_positions <- ref_pos[callable_col]
L_callable <- length(unique(callable_positions))

# 3.3 构建 0/1/NA 基因型矩阵（0=与参考相同；1=与参考不同；NA=非 A/C/G/T）
geno_cols <- which(callable_col)

geno_mat <- matrix(NA_integer_, nrow = nrow(sample_mat), ncol = length(geno_cols))
colnames(geno_mat) <- paste0("pos_", ref_pos[geno_cols])
rownames(geno_mat) <- rownames(sample_mat)

for (j in seq_along(geno_cols)) {
  k <- geno_cols[j]
  ref_base <- ref_aln[k]
  b <- sample_mat[, k]
  
  g <- rep(NA_integer_, length(b))
  ok <- is_acgt(b)
  g[ok & b == ref_base] <- 0L
  g[ok & b != ref_base] <- 1L
  geno_mat[, j] <- g
}

geno_df <- as.data.frame(geno_mat) %>% mutate(seq = rownames(geno_mat))

# 3.4 计算每个位点 MAF（对 0/1 来说，MAF=mean(1)）
maf_tbl <- geno_df %>%
  summarise(across(starts_with("pos_"), ~ mean(.x == 1, na.rm = TRUE))) %>%
  tidyr::pivot_longer(everything(), names_to = "snp", values_to = "maf") %>%
  mutate(
    n_mut = map_int(snp, ~ sum(geno_df[[.x]] == 1, na.rm = TRUE)),
    n_wt  = map_int(snp, ~ sum(geno_df[[.x]] == 0, na.rm = TRUE))
  )

# 按 Methods：仅保留 MAF > 1% 的 SNP
snp_keep <- maf_tbl %>% filter(maf > 0.01) %>% pull(snp)

# 3.4 修改版
maf_tbl <- geno_df %>%
  summarise(across(starts_with("pos_"), ~ mean(.x == 1, na.rm = TRUE))) %>%
  tidyr::pivot_longer(everything(), names_to = "snp", values_to = "p_alt") %>%
  mutate(
    n_alt = map_int(snp, ~ sum(geno_df[[.x]] == 1, na.rm = TRUE)),
    n_ref = map_int(snp, ~ sum(geno_df[[.x]] == 0, na.rm = TRUE)),
    n_called = n_alt + n_ref,
    maf = pmin(p_alt, 1 - p_alt),
    mac = pmin(n_alt, n_ref)
  )
snp_keep <- maf_tbl %>%
  filter(maf > 0.01, mac >= 5) %>%
  pull(snp)


# 合并回临床数据
dat1 <- dat0 %>%
  select(seq, ct_value, age, gender, days_onset_to_test, case_type) %>%
  left_join(geno_df %>% select(seq, all_of(snp_keep)), by = "seq")

# ---------- 4) Adjusted GLM：逐 SNP 拟合 ----------
fit_one_snp <- function(snp_name, data, df_age = 3, df_days = 3) {
  
  # 若该位点全是 0 或全是 1（或有效值太少），直接跳过
  x <- data[[snp_name]]
  ux <- sort(unique(x[!is.na(x)]))
  if (length(ux) < 2) {
    return(tibble(
      snp = snp_name, beta = NA_real_, conf.low = NA_real_, conf.high = NA_real_,
      p = NA_real_, n = nrow(data),
      n_mut = sum(x == 1, na.rm = TRUE), n_wt = sum(x == 0, na.rm = TRUE)
    ))
  }
  
  fml <- as.formula(paste0(
    "ct_value ~ ", snp_name,
    " + ns(age, df = ", df_age, ")",
    " + gender",
    " + ns(days_onset_to_test, df = ", df_days, ")",
    " + case_type"
  ))
  
  m <- lm(fml, data = data)
  
  broom::tidy(m, conf.int = TRUE) %>%
    filter(term == snp_name) %>%
    transmute(
      snp = snp_name,
      beta = estimate,
      conf.low = conf.low,
      conf.high = conf.high,
      p = p.value,
      n = nobs(m),
      n_mut = sum(data[[snp_name]] == 1, na.rm = TRUE),
      n_wt  = sum(data[[snp_name]] == 0, na.rm = TRUE)
    )
}

res_adj <- map_dfr(snp_keep, fit_one_snp, data = dat1, df_age = 3, df_days = 3)

# ---------- 5) Bonferroni（导师要求：按 callable 链条长度） ----------
alpha <- 0.05
res_adj <- res_adj %>%
  mutate(
    L_callable = L_callable,
    bonf_thr_callable = alpha / L_callable,
    bonf_sig_callable = ifelse(is.na(p), NA, p < bonf_thr_callable)
  ) %>%
  arrange(p)

# ---------- 6) 输出：主结果表 + 关键汇总 ----------
cat("Callable length (L_callable) =", L_callable, "\n")
cat("Tested SNPs (MAF>1%) =", length(snp_keep), "\n")
cat("Bonferroni threshold = ", alpha, "/", L_callable, " = ", alpha / L_callable, "\n\n")

res_adj