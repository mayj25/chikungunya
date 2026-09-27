#!/usr/bin/env Rscript

# District-level mosquito indices and Rt uncertainty-propagated Spearman analysis
#
# Optional command-line arguments:
#   1. district daily case CSV
#   2. point-level mosquito Excel workbook
#   3. output directory
#   4. number of Monte Carlo draws (default: 10000)

options(stringsAsFactors = FALSE, scipen = 999)
if (.Platform$OS.type == "windows") {
  invisible(suppressWarnings(Sys.setlocale("LC_ALL", "Chinese")))
}

args <- commandArgs(trailingOnly = TRUE)
case_file <- if (length(args) >= 1) args[[1]] else
  "Rt_uncertainty_spearman_district/input_snapshot/District_Cases_Level.csv"
mosquito_file <- if (length(args) >= 2) args[[2]] else
  "Rt_uncertainty_spearman_district/input_snapshot/mosquito_clean.xlsx"
output_dir <- if (length(args) >= 3) args[[3]] else
  file.path("Rt_uncertainty_spearman_district", "results")
n_draws <- if (length(args) >= 4) as.integer(args[[4]]) else 10000L

seed <- 20250901L
window_size <- 14L
mean_si <- 14.89
sd_si <- 9.95
max_si <- 50L
alpha_imported <- 0.5
prior_mean <- 5
prior_sd <- 5
lags <- 0:9
minimum_pairs <- 5L

required_packages <- c("data.table", "readxl", "ggplot2", "patchwork", "writexl")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(ggplot2)
  library(patchwork)
})

if (!file.exists(case_file)) stop("Case file not found: ", case_file)
if (!file.exists(mosquito_file)) stop("Mosquito file not found: ", mosquito_file)
if (!is.finite(n_draws) || n_draws < 1000L) stop("n_draws must be at least 1000")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(seed)

trim_text <- function(x) {
  gsub("[[:space:]]+", "", trimws(as.character(x)))
}

monday_of_week <- function(x) {
  x <- as.Date(x)
  x - ((as.POSIXlt(x)$wday + 6L) %% 7L)
}

iso_week_start <- function(year, week) {
  year <- as.integer(year)
  week <- as.integer(week)
  jan4 <- as.Date(sprintf("%04d-01-04", year))
  week1_monday <- monday_of_week(jan4)
  week1_monday + 7L * (week - 1L)
}

geometric_mean <- function(x) {
  if (length(x) == 0L || any(!is.finite(x)) || any(x <= 0)) return(NA_real_)
  exp(mean(log(x)))
}

make_discrete_si <- function(mean_si, sd_si, max_si) {
  shape <- (mean_si / sd_si)^2
  rate <- mean_si / (sd_si^2)
  out <- vapply(
    0:max_si,
    function(k) pgamma(k + 1, shape = shape, rate = rate) -
      pgamma(k, shape = shape, rate = rate),
    numeric(1)
  )
  out / sum(out)
}

calculate_lambda <- function(cases, si_distribution) {
  n <- length(cases)
  lambda <- numeric(n)
  max_lag <- length(si_distribution)
  if (n <= 1L) return(lambda)
  for (t in 2:n) {
    available_lags <- min(t - 1L, max_lag)
    lag_index <- seq_len(available_lags)
    lambda[t] <- sum(cases[t - lag_index] * si_distribution[lag_index])
  }
  lambda
}

estimate_daily_rt <- function(district_cases, si_distribution) {
  setorder(district_cases, date)
  local_cases <- district_cases$local
  imported_cases <- district_cases$imported
  effective_cases <- local_cases + alpha_imported * imported_cases
  lambda <- calculate_lambda(effective_cases, si_distribution)

  prior_shape <- (prior_mean / prior_sd)^2
  prior_scale <- prior_sd^2 / prior_mean
  n <- nrow(district_cases)
  if (n <= window_size) return(data.table())

  t_end <- seq.int(window_size + 1L, n)
  t_start <- t_end - window_size + 1L
  total_local <- vapply(
    seq_along(t_end),
    function(i) sum(local_cases[t_start[i]:t_end[i]]),
    numeric(1)
  )
  total_lambda <- vapply(
    seq_along(t_end),
    function(i) sum(lambda[t_start[i]:t_end[i]]),
    numeric(1)
  )

  posterior_shape <- prior_shape + total_local
  posterior_scale <- 1 / (1 / prior_scale + total_lambda)

  data.table(
    district = district_cases$district[t_end],
    date = district_cases$date[t_end],
    window_start = district_cases$date[t_start],
    window_end = district_cases$date[t_end],
    local_cases_window = total_local,
    infectiousness_window = total_lambda,
    posterior_shape = posterior_shape,
    posterior_scale = posterior_scale,
    rt_mean = posterior_shape * posterior_scale,
    rt_median = qgamma(0.5, shape = posterior_shape, scale = posterior_scale),
    rt_q025 = qgamma(0.025, shape = posterior_shape, scale = posterior_scale),
    rt_q975 = qgamma(0.975, shape = posterior_shape, scale = posterior_scale),
    week_start = monday_of_week(district_cases$date[t_end])
  )
}

rank_columns <- function(x) {
  if (ncol(x) == 1L) return(matrix(rank(x[, 1L], ties.method = "average"), ncol = 1L))
  vapply(
    seq_len(ncol(x)),
    function(j) rank(x[, j], ties.method = "average"),
    numeric(nrow(x))
  )
}

spearman_for_draws <- function(exposure, rt_draws) {
  if (length(exposure) < minimum_pairs || nrow(rt_draws) != length(exposure)) {
    return(rep(NA_real_, ncol(rt_draws)))
  }
  exposure_rank <- rank(exposure, ties.method = "average")
  exposure_centered <- exposure_rank - mean(exposure_rank)
  exposure_ss <- sum(exposure_centered^2)
  if (!is.finite(exposure_ss) || exposure_ss <= 0) {
    return(rep(NA_real_, ncol(rt_draws)))
  }

  rt_rank <- rank_columns(rt_draws)
  rt_centered <- sweep(rt_rank, 2L, colMeans(rt_rank), FUN = "-")
  rt_ss <- colSums(rt_centered^2)
  numerator <- as.vector(crossprod(exposure_centered, rt_centered))
  rho <- numerator / sqrt(exposure_ss * rt_ss)
  rho[!is.finite(rho)] <- NA_real_
  rho
}

summarize_correlation <- function(
    district, indicator, lag, exposure, rt_point, rt_draws, exposure_weeks) {
  n_pairs <- length(exposure)
  if (n_pairs < minimum_pairs || length(unique(exposure)) < 2L ||
      length(unique(rt_point)) < 2L) {
    return(data.table(
      district = district, indicator = indicator, lag_weeks = lag,
      n_pairs = n_pairs, first_exposure_week = min(exposure_weeks),
      last_exposure_week = max(exposure_weeks), rho_point = NA_real_,
      rho_median = NA_real_, rho_q025 = NA_real_, rho_q975 = NA_real_,
      probability_positive = NA_real_, probability_negative = NA_real_,
      directional_p = NA_real_, probability_same_direction = NA_real_,
      interval_excludes_zero = FALSE, significance = ""
    ))
  }

  rho_point <- suppressWarnings(cor(exposure, rt_point, method = "spearman"))
  rho_draw <- spearman_for_draws(exposure, rt_draws)
  rho_draw <- rho_draw[is.finite(rho_draw)]
  if (length(rho_draw) == 0L) {
    return(data.table(
      district = district, indicator = indicator, lag_weeks = lag,
      n_pairs = n_pairs, first_exposure_week = min(exposure_weeks),
      last_exposure_week = max(exposure_weeks), rho_point = rho_point,
      rho_median = NA_real_, rho_q025 = NA_real_, rho_q975 = NA_real_,
      probability_positive = NA_real_, probability_negative = NA_real_,
      directional_p = NA_real_, probability_same_direction = NA_real_,
      interval_excludes_zero = FALSE, significance = ""
    ))
  }

  interval <- unname(quantile(rho_draw, c(0.025, 0.5, 0.975), na.rm = TRUE))
  prob_positive <- mean(rho_draw > 0)
  prob_negative <- mean(rho_draw < 0)
  directional_p <- min(1, 2 * min(prob_positive, prob_negative))
  probability_same_direction <- if (is.na(rho_point) || rho_point == 0) {
    max(prob_positive, prob_negative)
  } else if (rho_point > 0) {
    prob_positive
  } else {
    prob_negative
  }
  significance <- if (directional_p <= 0.001) {
    "***"
  } else if (directional_p <= 0.01) {
    "**"
  } else if (directional_p <= 0.05) {
    "*"
  } else {
    ""
  }

  data.table(
    district = district,
    indicator = indicator,
    lag_weeks = lag,
    n_pairs = n_pairs,
    first_exposure_week = min(exposure_weeks),
    last_exposure_week = max(exposure_weeks),
    rho_point = rho_point,
    rho_median = interval[2L],
    rho_q025 = interval[1L],
    rho_q975 = interval[3L],
    probability_positive = prob_positive,
    probability_negative = prob_negative,
    directional_p = directional_p,
    probability_same_direction = probability_same_direction,
    interval_excludes_zero = interval[1L] > 0 || interval[3L] < 0,
    significance = significance
  )
}

# -----------------------------------------------------------------------------
# 1. Input and quality control
# -----------------------------------------------------------------------------

cases_raw <- as.data.table(read.csv(
  case_file, fileEncoding = "UTF-8-BOM", check.names = FALSE
))
required_case_columns <- c("district", "date", "local", "imported")
if (!all(required_case_columns %in% names(cases_raw))) {
  stop("Case file must contain: ", paste(required_case_columns, collapse = ", "))
}
cases_raw[, district := trim_text(district)]
cases_raw[, date := as.IDate(date)]
cases_raw[, `:=`(local = as.numeric(local), imported = as.numeric(imported))]
if (anyNA(cases_raw[, ..required_case_columns])) stop("Missing values in required case columns")
if (any(cases_raw$local < 0 | cases_raw$imported < 0)) stop("Negative case counts found")
if (anyDuplicated(cases_raw, by = c("district", "date"))) {
  stop("Duplicate district-date rows found in case file")
}

all_dates <- seq(min(cases_raw$date), max(cases_raw$date), by = "day")
case_grid <- CJ(district = sort(unique(cases_raw$district)), date = as.IDate(all_dates))
cases <- merge(
  case_grid,
  cases_raw[, .(district, date, local, imported)],
  by = c("district", "date"), all.x = TRUE, sort = TRUE
)
cases[, imputed_zero_day := is.na(local) & is.na(imported)]
cases[imputed_zero_day == TRUE, `:=`(local = 0, imported = 0)]
if (anyNA(cases[, .(local, imported)])) {
  stop("Only one of local/imported was missing on at least one observed row")
}
cases[, day_index_complete := as.integer(date - min(date)) + 1L]

case_completion_audit <- cases[imputed_zero_day == TRUE]

mosquito_raw <- as.data.table(read_excel(mosquito_file, sheet = 1))
required_mosquito_columns <- c("week", "district", "BI", "MOI")
if (!all(required_mosquito_columns %in% names(mosquito_raw))) {
  stop("Mosquito file must contain: ", paste(required_mosquito_columns, collapse = ", "))
}
mosquito_raw[, district := trim_text(district)]
mosquito_raw[, `:=`(
  week = as.integer(week),
  BI = suppressWarnings(as.numeric(BI)),
  MOI = suppressWarnings(as.numeric(MOI))
)]
if (anyNA(mosquito_raw$week)) stop("Missing/non-numeric week in mosquito file")
if (any(mosquito_raw$week < 1L | mosquito_raw$week > 53L)) stop("Invalid week number")

case_districts <- sort(unique(cases$district))
mosquito_districts <- sort(unique(mosquito_raw$district))
unmatched_case_districts <- setdiff(case_districts, mosquito_districts)
unmatched_mosquito_districts <- setdiff(mosquito_districts, case_districts)
if (length(unmatched_case_districts) > 0L || length(unmatched_mosquito_districts) > 0L) {
  stop(
    "District mismatch. Case-only: ", paste(unmatched_case_districts, collapse = ", "),
    "; mosquito-only: ", paste(unmatched_mosquito_districts, collapse = ", ")
  )
}

mosquito_raw[, year := 2025L]
mosquito_raw[, week_start := as.IDate(iso_week_start(year, week))]
mosquito_weekly <- mosquito_raw[, .(
  BI_mean = if (all(is.na(BI))) NA_real_ else mean(BI, na.rm = TRUE),
  MOI_mean = if (all(is.na(MOI))) NA_real_ else mean(MOI, na.rm = TRUE),
  n_points = .N,
  n_BI_nonmissing = sum(!is.na(BI)),
  n_MOI_nonmissing = sum(!is.na(MOI))
), by = .(district, year, week, week_start)]
setorder(mosquito_weekly, district, week_start)

# -----------------------------------------------------------------------------
# 2. Exact district Rt model from the supplied code, retaining Gamma posteriors
# -----------------------------------------------------------------------------

si_distribution <- make_discrete_si(mean_si, sd_si, max_si)
daily_rt <- rbindlist(
  lapply(split(cases, by = "district", keep.by = TRUE), estimate_daily_rt,
         si_distribution = si_distribution),
  use.names = TRUE
)
setorder(daily_rt, district, date)

# -----------------------------------------------------------------------------
# 3. Posterior trajectories, weekly geometric means, and lagged correlations
# -----------------------------------------------------------------------------

weekly_rt_list <- vector("list", length(case_districts))
correlation_list <- vector("list", length(case_districts) * 2L * length(lags))
cor_index <- 0L

for (district_name in case_districts) {
  district_rt <- daily_rt[district == district_name]
  n_days <- nrow(district_rt)
  daily_draws <- matrix(
    rgamma(
      n_days * n_draws,
      shape = rep(district_rt$posterior_shape, times = n_draws),
      scale = rep(district_rt$posterior_scale, times = n_draws)
    ),
    nrow = n_days,
    ncol = n_draws
  )

  week_values <- sort(unique(district_rt$week_start))
  weekly_draws <- matrix(
    NA_real_, nrow = length(week_values), ncol = n_draws,
    dimnames = list(as.character(week_values), NULL)
  )
  weekly_point <- numeric(length(week_values))
  weekly_n_days <- integer(length(week_values))

  for (i in seq_along(week_values)) {
    rows <- which(district_rt$week_start == week_values[i])
    weekly_draws[i, ] <- exp(colMeans(log(daily_draws[rows, , drop = FALSE])))
    weekly_point[i] <- geometric_mean(district_rt$rt_median[rows])
    weekly_n_days[i] <- length(rows)
  }

  weekly_quantiles <- t(apply(
    weekly_draws, 1L, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE
  ))
  weekly_rt <- data.table(
    district = district_name,
    week_start = as.IDate(week_values),
    week = as.integer(format(as.Date(week_values), "%V")),
    n_daily_rt_values = weekly_n_days,
    rt_geometric_point = weekly_point,
    rt_geometric_median = weekly_quantiles[, 2L],
    rt_geometric_q025 = weekly_quantiles[, 1L],
    rt_geometric_q975 = weekly_quantiles[, 3L]
  )
  weekly_rt_list[[match(district_name, case_districts)]] <- weekly_rt

  district_mosquito <- mosquito_weekly[district == district_name]
  for (indicator_name in c("BI", "MOI")) {
    indicator_column <- paste0(indicator_name, "_mean")
    for (lag_value in lags) {
      exposure_data <- district_mosquito[
        !is.na(get(indicator_column)),
        .(
          exposure_week = week_start,
          target_rt_week = week_start + 7L * lag_value,
          exposure = get(indicator_column)
        )
      ]
      exposure_data[, rt_row := match(target_rt_week, weekly_rt$week_start)]
      exposure_data <- exposure_data[!is.na(rt_row)]

      cor_index <- cor_index + 1L
      if (nrow(exposure_data) == 0L) {
        correlation_list[[cor_index]] <- data.table(
          district = district_name, indicator = indicator_name,
          lag_weeks = lag_value, n_pairs = 0L,
          first_exposure_week = as.IDate(NA), last_exposure_week = as.IDate(NA),
          rho_point = NA_real_, rho_median = NA_real_,
          rho_q025 = NA_real_, rho_q975 = NA_real_,
          probability_positive = NA_real_, probability_negative = NA_real_,
          directional_p = NA_real_, probability_same_direction = NA_real_,
          interval_excludes_zero = FALSE, significance = ""
        )
      } else {
        correlation_list[[cor_index]] <- summarize_correlation(
          district = district_name,
          indicator = indicator_name,
          lag = lag_value,
          exposure = exposure_data$exposure,
          rt_point = weekly_rt$rt_geometric_point[exposure_data$rt_row],
          rt_draws = weekly_draws[exposure_data$rt_row, , drop = FALSE],
          exposure_weeks = exposure_data$exposure_week
        )
      }
    }
  }
}

weekly_rt_all <- rbindlist(weekly_rt_list, use.names = TRUE)
correlation_results <- rbindlist(correlation_list, use.names = TRUE, fill = TRUE)
setorder(correlation_results, indicator, district, lag_weeks)
correlation_results[, positive_significance := fcase(
  !is.na(probability_positive) & probability_positive >= 0.999, "***",
  !is.na(probability_positive) & probability_positive >= 0.99, "**",
  !is.na(probability_positive) & probability_positive >= 0.95, "*",
  default = ""
)]

# -----------------------------------------------------------------------------
# 4. Publication-style heat maps
# -----------------------------------------------------------------------------

district_en_map <- c(
  "\u767d\u4e91\u533a" = "Baiyun", "\u4ece\u5316\u533a" = "Conghua",
  "\u6d77\u73e0\u533a" = "Haizhu", "\u82b1\u90fd\u533a" = "Huadu",
  "\u9ec4\u57d4\u533a" = "Huangpu", "\u8354\u6e7e\u533a" = "Liwan",
  "\u5357\u6c99\u533a" = "Nansha", "\u756a\u79ba\u533a" = "Panyu",
  "\u5929\u6cb3\u533a" = "Tianhe", "\u8d8a\u79c0\u533a" = "Yuexiu",
  "\u589e\u57ce\u533a" = "Zengcheng"
)
district_order <- c(
  "Baiyun", "Conghua", "Haizhu", "Huadu", "Huangpu", "Liwan",
  "Nansha", "Panyu", "Tianhe", "Yuexiu", "Zengcheng"
)

plot_data <- copy(correlation_results)
plot_data[, district_en := unname(district_en_map[district])]
if (anyNA(plot_data$district_en)) {
  stop("English district label missing for: ",
       paste(unique(plot_data[is.na(district_en), district]), collapse = ", "))
}
plot_data[, district_en := factor(district_en, levels = rev(district_order))]
plot_data[, lag_label := factor(
  paste("Lag", lag_weeks), levels = paste("Lag", lags)
)]

make_heatmap <- function(indicator_name, show_legend = TRUE) {
  p <- ggplot(plot_data[indicator == indicator_name],
              aes(x = lag_label, y = district_en, fill = rho_median)) +
    geom_tile(color = "white", linewidth = 0.45) +
    geom_text(aes(label = positive_significance), size = 3.1,
              fontface = "bold", na.rm = TRUE) +
    scale_fill_gradient2(
      low = "#3B73B9", mid = "#F7F7F7", high = "#E65F51",
      midpoint = 0, limits = c(-1, 1), na.value = "#ECECEC",
      name = "Spearman\nposterior median"
    ) +
    labs(
      title = NULL,
      x = "Lag (weeks)", y = "District"
    ) +
    theme_minimal(base_size = 10.5) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.title = element_text(face = "plain"),
      legend.position = if (show_legend) "right" else "none",
      plot.margin = margin(6, 8, 6, 6)
    )
  p
}

combined_plot <- make_heatmap("BI", show_legend = FALSE) +
  make_heatmap("MOI", show_legend = TRUE) +
  plot_layout(widths = c(1, 1.13))

original_working_dir <- getwd()
setwd(output_dir)
ggsave(
  "figure_district_spearman_rt_uncertainty_clean.png",
  combined_plot, width = 13.5, height = 5.7, units = "in", dpi = 400,
  device = grDevices::png, type = "cairo", bg = "white"
)
ggsave(
  "figure_district_spearman_rt_uncertainty_clean.pdf",
  combined_plot, width = 13.5, height = 5.7, units = "in", device = cairo_pdf,
  bg = "white"
)
ggsave(
  "figure_district_spearman_rt_uncertainty_clean.tiff",
  combined_plot, width = 13.5, height = 5.7, units = "in", dpi = 600,
  device = grDevices::tiff, type = "cairo", compression = "lzw", bg = "white"
)
setwd(original_working_dir)

# -----------------------------------------------------------------------------
# 5. Reproducible tables and audit record
# -----------------------------------------------------------------------------

write_utf8_csv <- function(x, path) {
  out <- copy(x)
  character_columns <- names(out)[vapply(out, is.character, logical(1))]
  if (length(character_columns) > 0L) {
    out[, (character_columns) := lapply(.SD, enc2utf8), .SDcols = character_columns]
  }
  fwrite(out, path, bom = TRUE)
}

write_utf8_csv(mosquito_weekly,
               file.path(output_dir, "district_weekly_mosquito_means.csv"))
write_utf8_csv(cases, file.path(output_dir, "district_daily_cases_completed.csv"))
write_utf8_csv(case_completion_audit,
               file.path(output_dir, "case_date_completion_audit.csv"))
write_utf8_csv(daily_rt,
               file.path(output_dir, "district_daily_rt_posterior_parameters.csv"))
write_utf8_csv(weekly_rt_all,
               file.path(output_dir, "district_weekly_rt_uncertainty.csv"))
write_utf8_csv(correlation_results,
               file.path(output_dir, "district_spearman_rt_uncertainty_results.csv"))

writexl::write_xlsx(
  list(
    mosquito_weekly_means = as.data.frame(mosquito_weekly),
    completed_daily_cases = as.data.frame(cases),
    imputed_zero_dates = as.data.frame(case_completion_audit),
    daily_rt_posteriors = as.data.frame(daily_rt),
    weekly_rt_uncertainty = as.data.frame(weekly_rt_all),
    spearman_uncertainty = as.data.frame(correlation_results)
  ),
  path = file.path(output_dir, "district_rt_uncertainty_analysis.xlsx")
)

partial_rt_weeks <- weekly_rt_all[n_daily_rt_values < 7L]
audit_lines <- c(
  "District Rt uncertainty-propagated Spearman analysis",
  paste("Run time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste("R version:", R.version.string),
  paste("Case input:", normalizePath(case_file, winslash = "/")),
  paste("Mosquito input:", normalizePath(mosquito_file, winslash = "/")),
  paste("Output directory:", normalizePath(output_dir, winslash = "/")),
  "",
  paste("Districts:", length(case_districts)),
  paste("Raw case rows:", nrow(cases_raw)),
  paste("Completed case rows:", nrow(cases)),
  paste("Calendar dates added as zero-case days:", nrow(case_completion_audit)),
  paste("Distinct added dates per district:",
        length(unique(case_completion_audit$date))),
  paste("Case date range:", min(cases$date), "to", max(cases$date)),
  paste("Raw mosquito point rows:", nrow(mosquito_raw)),
  paste("District-week mosquito rows:", nrow(mosquito_weekly)),
  paste("Mosquito week range:", min(mosquito_weekly$week), "to",
        max(mosquito_weekly$week)),
  paste("Raw missing BI:", sum(is.na(mosquito_raw$BI))),
  paste("Raw missing MOI:", sum(is.na(mosquito_raw$MOI))),
  "",
  paste("Serial interval mean/sd/max:", mean_si, "/", sd_si, "/", max_si),
  paste("Imported-case relative infectiousness:", alpha_imported),
  paste("Rt window (days):", window_size),
  paste("Gamma prior mean/sd:", prior_mean, "/", prior_sd),
  paste("Monte Carlo draws:", n_draws),
  paste("Random seed:", seed),
  paste("Lag range:", paste(range(lags), collapse = " to "), "weeks"),
  paste("Minimum paired weeks:", minimum_pairs),
  "Lag definition: mosquito index in week w is correlated with Rt in week w + lag.",
  "Weekly Rt: geometric mean of available daily Rt values within each ISO week.",
  "Uncertainty propagation: independent draws from each daily Gamma Rt posterior,",
  "followed by weekly geometric aggregation and Spearman correlation for each draw.",
  "The correlation interval therefore represents Rt-estimation uncertainty only;",
  "it does not include uncertainty in BI/MOI, serial interval, case reporting, or",
  "cross-day posterior dependence induced by overlapping Rt windows.",
  "",
  paste("Partial district-weeks with <7 daily Rt values:", nrow(partial_rt_weeks)),
  if (nrow(partial_rt_weeks) > 0L) paste(capture.output(print(
    partial_rt_weeks[, .(district, week_start, week, n_daily_rt_values)]
  )), collapse = "\n") else "None",
  "",
  "Posterior directional probability is a Monte Carlo posterior tail-area summary,",
  "not a classical frequentist p-value. The 95% interval is the primary uncertainty result."
)
audit_connection <- file(
  file.path(output_dir, "analysis_audit.txt"), open = "wt", encoding = "UTF-8"
)
writeLines(audit_lines, audit_connection)
close(audit_connection)

cat("Analysis completed successfully.\n")
cat("Output directory:", normalizePath(output_dir, winslash = "/"), "\n")
cat("Districts:", length(case_districts), "\n")
cat("Monte Carlo draws:", n_draws, "\n")
cat("Correlation rows:", nrow(correlation_results), "\n")
