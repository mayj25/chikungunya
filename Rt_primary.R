
#20260301Rt计算（已按最新SI更新）
# ============================================================
# Rt estimation (alpha=0.5 only): Guangzhou city + districts
# Input:
#   City     : D:/基孔肯亚热/大文章数据分析/蚊媒可视化/Rt和蚊媒/df_city.csv
#   District : D:/基孔肯亚热/大文章数据分析/蚊媒可视化/Rt和蚊媒/District_Cases_Level.csv
# Output:
#   D:/基孔肯亚热/大文章数据分析/蚊媒可视化/Rt和蚊媒/output
#
# SI (UPDATED, final):
#   mean = 14.89 days, SD = 9.95
#   Gamma params: shape = 2.44, scale = 6.65
# ============================================================

# -------------------------
# 0) Packages
# -------------------------
pkgs <- c("dplyr","readr","stringr","writexl","purrr")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(dplyr)
library(readr)
library(stringr)
library(writexl)
library(purrr)

# -------------------------
# 1) Paths
# -------------------------
city_path <- "D:/基孔肯雅热/大文章数据分析/蚊媒可视化/Rt和蚊媒/df_city.csv"
dist_path <- "D:/基孔肯雅热/大文章数据分析/蚊媒可视化/Rt和蚊媒/District_Cases_Level.csv"
out_dir   <- "D:/基孔肯雅热/大文章数据分析/蚊媒可视化/Rt和蚊媒/output"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# -------------------------
# 2) Core functions
# -------------------------

# 2.1 discrete SI distribution (from mean & sd)
get_discrete_si <- function(mean_si, std_si, max_si) {
  shape <- (mean_si / std_si)^2
  scale <- (std_si^2) / mean_si
  si_distr <- sapply(0:max_si, function(k) {
    pgamma(k + 1, shape = shape, scale = scale) -
      pgamma(k, shape = shape, scale = scale)
  })
  return(si_distr / sum(si_distr))
}

# (Optional) discrete SI distribution (directly from shape & scale)
get_discrete_si_shape_scale <- function(shape, scale, max_si) {
  si_distr <- sapply(0:max_si, function(k) {
    pgamma(k + 1, shape = shape, scale = scale) -
      pgamma(k, shape = shape, scale = scale)
  })
  return(si_distr / sum(si_distr))
}

# 2.2 lambda (infectiousness), with alpha
calculate_lambda <- function(df, si_distr, alpha) {
  n <- nrow(df)
  lambda <- numeric(n)
  si_len <- length(si_distr)
  
  for (t in 1:n) {
    lambda_sum <- 0
    max_s <- min(t - 1, si_len)
    if (max_s >= 1) {
      for (s in 1:max_s) {
        idx <- t - s
        loc <- ifelse(is.na(df$local[idx]), 0, df$local[idx])
        imp <- ifelse(is.na(df$imported[idx]), 0, df$imported[idx])
        effective_cases <- loc + (imp * alpha)
        lambda_sum <- lambda_sum + effective_cases * si_distr[s]
      }
    }
    lambda[t] <- lambda_sum
  }
  return(lambda)
}

# 2.3 Bayesian Rt with window
calculate_rt_bayes <- function(df, lambda, window = 10, si_mean) {
  n <- nrow(df)
  t_start <- seq(2, n - window + 1)
  t_end <- t_start + window - 1
  
  res <- data.frame(
    t_end    = t_end,
    date     = df$date[t_end],
    Mean_R   = NA_real_,
    Median_R = NA_real_,
    Q0.025   = NA_real_,
    Q0.975   = NA_real_
  )
  
  # prior
  prior_mean <- 5; prior_sd <- 5
  a_prior <- (prior_mean / prior_sd)^2
  b_prior <- prior_sd^2 / prior_mean
  
  for (i in seq_along(t_start)) {
    ts <- t_start[i]; te <- t_end[i]
    if (te > si_mean/2) {
      total_local  <- sum(df$local[ts:te],  na.rm = TRUE)
      total_lambda <- sum(lambda[ts:te],    na.rm = TRUE)
      
      a_post <- a_prior + total_local
      b_post <- 1 / (1/b_prior + total_lambda)
      
      if (total_lambda > 0) {
        res$Mean_R[i]   <- a_post * b_post
        res$Median_R[i] <- qgamma(0.5,   shape = a_post, scale = b_post)
        res$Q0.025[i]   <- qgamma(0.025, shape = a_post, scale = b_post)
        res$Q0.975[i]   <- qgamma(0.975, shape = a_post, scale = b_post)
      }
    }
  }
  return(res)
}

# 2.4 batch regions Rt
estimate_rt_for_regions <- function(data, group_col_name, si_distr, alpha, window = 14, si_mean) {
  
  regions <- unique(data[[group_col_name]])
  all_results <- list()
  
  cat("开始批量计算，共", length(regions), "个区域...\n")
  
  for (region in regions) {
    sub_data <- data %>%
      filter(.data[[group_col_name]] == region) %>%
      arrange(date)
    
    sub_lambda <- calculate_lambda(sub_data, si_distr, alpha)
    sub_rt <- calculate_rt_bayes(sub_data, sub_lambda, window = window, si_mean = si_mean)
    
    if (nrow(sub_rt) > 0) {
      sub_rt[[group_col_name]] <- region
      all_results[[as.character(region)]] <- sub_rt
    }
  }
  
  final_df <- bind_rows(all_results) %>%
    select(all_of(group_col_name), everything())
  
  cat("计算完成。\n")
  return(final_df)
}

# -------------------------
# 3) Read data (handle date like 2025/7/9)
# -------------------------
df_city <- read_csv(city_path, show_col_types = FALSE)
df_district <- read_csv(dist_path, show_col_types = FALSE)

parse_date_slash <- function(x) as.Date(str_replace_all(as.character(x), "/", "-"))

df_city <- df_city %>%
  mutate(
    date = parse_date_slash(date),
    local = as.numeric(local),
    imported = as.numeric(imported)
  )

df_district <- df_district %>%
  mutate(
    date = parse_date_slash(date),
    district = as.character(district),
    local = as.numeric(local),
    imported = as.numeric(imported)
  )

# checks
need_city <- c("date","local","imported")
miss_city <- setdiff(need_city, names(df_city))
if (length(miss_city) > 0) stop("df_city 缺少列：", paste(miss_city, collapse = ", "))

need_dist <- c("district","date","local","imported")
miss_dist <- setdiff(need_dist, names(df_district))
if (length(miss_dist) > 0) stop("District_Cases_Level 缺少列：", paste(miss_dist, collapse = ", "))

df_city <- df_city %>% arrange(date)
df_district <- df_district %>% arrange(district, date)

# -------------------------
# 4) Fixed parameters (UPDATED SI)
# -------------------------
alpha_val   <- 0.5
window_val  <- 14

# ✅ Updated SI (final)
mu_si       <- 14.89
sigma_si    <- 9.95
max_si_days <- 50

# Provided Gamma parameters (for cross-check)
shape_si <- 2.44
scale_si <- 6.65

# Main: use mean/sd to generate SI (keeps your original pipeline unchanged)
si_distr <- get_discrete_si(mean_si = mu_si, std_si = sigma_si, max_si = max_si_days)

# Optional check (not used downstream unless you switch)
si_distr_shape <- get_discrete_si_shape_scale(shape = shape_si, scale = scale_si, max_si = max_si_days)

cat("\n[SI check]\n")
cat("From mean/sd -> shape =", round((mu_si / sigma_si)^2, 4),
    ", scale =", round((sigma_si^2) / mu_si, 4), "\n")
cat("Given shape/scale -> shape =", shape_si, ", scale =", scale_si, "\n")
cat("L1 difference between two discrete SI distr =",
    signif(sum(abs(si_distr - si_distr_shape)), 6), "\n\n")

# -------------------------
# 5) City Rt (alpha=0.5, window=14)
# -------------------------
cat("开始计算：广州市 Rt（alpha=0.5, window=14）...\n")
lambda_city <- calculate_lambda(df_city, si_distr, alpha = alpha_val)
Rt_city <- calculate_rt_bayes(df_city, lambda_city, window = window_val, si_mean = mu_si) %>%
  mutate(alpha = alpha_val, window = window_val) %>%
  arrange(date)

# -------------------------
# 6) District Rt (alpha=0.5, window=14)
# -------------------------
cat("开始计算：各区 Rt（alpha=0.5, window=14）...\n")
Rt_district <- estimate_rt_for_regions(
  data = df_district,
  group_col_name = "district",
  si_distr = si_distr,
  alpha = alpha_val,
  window = window_val,
  si_mean = mu_si
) %>%
  mutate(alpha = alpha_val, window = window_val) %>%
  arrange(district, date)

# -------------------------
# 7) Output (CSV + Excel)
# -------------------------
city_csv <- file.path(out_dir, paste0("Rt_city_alpha", alpha_val, "_window", window_val, "_SImean", mu_si, "_SIsd", sigma_si, ".csv"))
dist_csv <- file.path(out_dir, paste0("Rt_district_alpha", alpha_val, "_window", window_val, "_SImean", mu_si, "_SIsd", sigma_si, ".csv"))
xlsx_out <- file.path(out_dir, paste0("Rt_alpha", alpha_val, "_window", window_val, "_SImean", mu_si, "_SIsd", sigma_si, ".xlsx"))

write_csv(Rt_city, city_csv)
write_csv(Rt_district, dist_csv)

write_xlsx(
  list(
    Rt_city = Rt_city,
    Rt_district = Rt_district
  ),
  path = xlsx_out
)

cat("\n✅ 计算完成！输出目录：", out_dir, "\n")
cat("✅ City CSV :", city_csv, "\n")
cat("✅ Dist CSV :", dist_csv, "\n")
cat("✅ Excel    :", xlsx_out, "\n")

