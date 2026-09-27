# Rt sensitivity analysis for imported-case relative infectiousness (alpha)
#
# Reproduces the calculation in Rt_computation.R while varying alpha from
# 0.0 to 1.0 by 0.1.  Dates are used exactly as supplied in the input files;
# this script does not add missing calendar dates.

required_pkgs <- c("dplyr", "readr", "stringr")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace,
                                      logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required R package(s): ", paste(missing_pkgs, collapse = ", "))
}

library(dplyr)
library(readr)
library(stringr)

# Locate the Rt data folder from the working directory. This avoids R's
# --file path-encoding issue on Windows when parent directories contain Chinese
# characters. Run either from the project root or from the Rt folder itself.
candidate_dirs <- c(file.path(getwd(), "Rt"), getwd())
candidate_dirs <- candidate_dirs[file.exists(file.path(candidate_dirs, "df_city.csv"))]
if (length(candidate_dirs) == 0) {
  stop("Could not locate df_city.csv. Run this script from the project root or the Rt folder.")
}
script_dir <- normalizePath(candidate_dirs[1])

city_path <- file.path(script_dir, "df_city.csv")
dist_path <- file.path(script_dir, "District_Cases_Level.csv")
out_dir <- file.path(script_dir, "output_rt_alpha_0.0_to_1.0")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Discrete Gamma serial-interval distribution from mean and standard deviation.
get_discrete_si <- function(mean_si, std_si, max_si) {
  shape <- (mean_si / std_si)^2
  scale <- (std_si^2) / mean_si
  si_distr <- vapply(0:max_si, function(k) {
    pgamma(k + 1, shape = shape, scale = scale) -
      pgamma(k, shape = shape, scale = scale)
  }, numeric(1))
  si_distr / sum(si_distr)
}

# Total infectiousness at t. Imported cases contribute alpha times as much as
# local cases, exactly as in the supplied Rt_computation.R.
calculate_lambda <- function(df, si_distr, alpha) {
  n <- nrow(df)
  lambda <- numeric(n)
  si_len <- length(si_distr)

  for (t in seq_len(n)) {
    max_s <- min(t - 1, si_len)
    if (max_s >= 1) {
      for (s in seq_len(max_s)) {
        idx <- t - s
        local_cases <- ifelse(is.na(df$local[idx]), 0, df$local[idx])
        imported_cases <- ifelse(is.na(df$imported[idx]), 0, df$imported[idx])
        lambda[t] <- lambda[t] +
          (local_cases + alpha * imported_cases) * si_distr[s]
      }
    }
  }
  lambda
}

# Bayesian renewal-model Rt estimate with the same Gamma(1, scale = 5) prior
# used in Rt_computation.R (prior mean = 5; prior SD = 5).
calculate_rt_bayes <- function(df, lambda, window, si_mean) {
  n <- nrow(df)
  if (n <= window) stop("Input has too few rows for the requested Rt window.")

  t_start <- seq(2, n - window + 1)
  t_end <- t_start + window - 1
  res <- data.frame(
    t_end = t_end,
    date = df$date[t_end],
    Mean_R = NA_real_,
    Median_R = NA_real_,
    Q0.025 = NA_real_,
    Q0.975 = NA_real_
  )

  prior_mean <- 5
  prior_sd <- 5
  a_prior <- (prior_mean / prior_sd)^2
  b_prior <- prior_sd^2 / prior_mean

  for (i in seq_along(t_start)) {
    ts <- t_start[i]
    te <- t_end[i]
    if (te > si_mean / 2) {
      total_local <- sum(df$local[ts:te], na.rm = TRUE)
      total_lambda <- sum(lambda[ts:te], na.rm = TRUE)
      a_post <- a_prior + total_local
      b_post <- 1 / (1 / b_prior + total_lambda)

      if (total_lambda > 0) {
        res$Mean_R[i] <- a_post * b_post
        res$Median_R[i] <- qgamma(0.5, shape = a_post, scale = b_post)
        res$Q0.025[i] <- qgamma(0.025, shape = a_post, scale = b_post)
        res$Q0.975[i] <- qgamma(0.975, shape = a_post, scale = b_post)
      }
    }
  }
  res
}

estimate_rt_for_regions <- function(data, si_distr, alpha, window, si_mean) {
  regions <- unique(data$district)
  results <- lapply(regions, function(region) {
    sub_data <- data %>% filter(district == region) %>% arrange(date)
    sub_lambda <- calculate_lambda(sub_data, si_distr, alpha)
    calculate_rt_bayes(sub_data, sub_lambda, window, si_mean) %>%
      mutate(district = region, .before = 1)
  })
  bind_rows(results)
}

parse_date_slash <- function(x) as.Date(str_replace_all(as.character(x), "/", "-"))

df_city <- read_csv(city_path, show_col_types = FALSE) %>%
  mutate(date = parse_date_slash(date), local = as.numeric(local), imported = as.numeric(imported)) %>%
  arrange(date)
df_district <- read_csv(dist_path, show_col_types = FALSE) %>%
  mutate(date = parse_date_slash(date), district = as.character(district),
         local = as.numeric(local), imported = as.numeric(imported)) %>%
  arrange(district, date)

stopifnot(all(c("date", "local", "imported") %in% names(df_city)))
stopifnot(all(c("district", "date", "local", "imported") %in% names(df_district)))
stopifnot(!anyDuplicated(df_city$date))
stopifnot(!anyDuplicated(df_district[c("district", "date")]))
stopifnot(sum(df_city$local, na.rm = TRUE) == sum(df_district$local, na.rm = TRUE))
stopifnot(sum(df_city$imported, na.rm = TRUE) == sum(df_district$imported, na.rm = TRUE))

# Fixed analysis settings retained from Rt_computation.R.
alpha_values <- seq(0, 1, by = 0.1)
window_val <- 14
mu_si <- 14.89
sigma_si <- 9.95
max_si_days <- 50
si_distr <- get_discrete_si(mu_si, sigma_si, max_si_days)

city_results <- list()
district_results <- list()

for (alpha_val in alpha_values) {
  alpha_label <- sprintf("%.1f", alpha_val)
  message("Calculating alpha = ", alpha_label)

  city_rt <- calculate_rt_bayes(
    df_city,
    calculate_lambda(df_city, si_distr, alpha_val),
    window = window_val,
    si_mean = mu_si
  ) %>%
    mutate(alpha = alpha_val, window = window_val) %>%
    arrange(date)

  district_rt <- estimate_rt_for_regions(
    df_district, si_distr, alpha_val, window_val, mu_si
  ) %>%
    mutate(alpha = alpha_val, window = window_val) %>%
    arrange(district, date)

  city_results[[alpha_label]] <- city_rt
  district_results[[alpha_label]] <- district_rt

  write_csv(city_rt, file.path(out_dir, paste0("Rt_city_alpha", alpha_label, ".csv")))
  write_csv(district_rt, file.path(out_dir, paste0("Rt_district_alpha", alpha_label, ".csv")))
}

city_all <- bind_rows(city_results)
district_all <- bind_rows(district_results)
run_parameters <- data.frame(
  alpha = alpha_values,
  imported_case_relative_infectiousness = alpha_values,
  window_days = window_val,
  si_mean_days = mu_si,
  si_sd_days = sigma_si,
  si_max_days = max_si_days,
  prior_mean_R = 5,
  prior_sd_R = 5,
  city_input_rows = nrow(df_city),
  district_input_rows = nrow(df_district),
  stringsAsFactors = FALSE
)

write_csv(city_all, file.path(out_dir, "Rt_city_all_alpha.csv"))
write_csv(district_all, file.path(out_dir, "Rt_district_all_alpha.csv"))
write_csv(run_parameters, file.path(out_dir, "run_parameters.csv"))

# Reproducibility check against the supplied alpha=0.5 city result, if present.
reference_path <- file.path(script_dir, "Rt_city_alpha0.5_window14_SImean14.89_SIsd9.95_xueyi.csv")
if (file.exists(reference_path)) {
  reference <- read_csv(reference_path, show_col_types = FALSE) %>% arrange(date)
  new_alpha_05 <- city_results[["0.5"]] %>% arrange(date)
  value_cols <- c("Mean_R", "Median_R", "Q0.025", "Q0.975")
  stopifnot(nrow(reference) == nrow(new_alpha_05))
  stopifnot(all(reference$date == as.character(new_alpha_05$date)))
  max_abs_difference <- max(abs(as.matrix(reference[value_cols]) -
                                as.matrix(new_alpha_05[value_cols])), na.rm = TRUE)
  if (!is.finite(max_abs_difference) || max_abs_difference > 1e-10) {
    stop("Alpha=0.5 reproducibility check failed; maximum absolute difference = ",
         max_abs_difference)
  }
  write_csv(data.frame(alpha = 0.5, max_abs_difference = max_abs_difference),
            file.path(out_dir, "alpha0.5_reproducibility_check.csv"))
}

message("Completed. Results saved to: ", normalizePath(out_dir))
