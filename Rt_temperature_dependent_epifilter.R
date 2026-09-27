# Temperature-dependent EpiFilter Rt analysis for the 2025 Guangzhou
# chikungunya outbreak.  This script deliberately does not overwrite the
# fixed-serial-interval EpiFilter outputs in this folder.
#
# Method:
#   G(T) = I_H + T_HM + I_M(T) + T_MH
# and Lambda_t = sum_{u>=1} I*_{t-u} w_u(T_{t-u}), where I* gives imported
# cases a user-specified relative infectiousness.  Thus weather on the source
# case's date, rather than weather on the reporting date, determines its
# generation-interval kernel.
#
# Run (from this folder):
#   "C:/Program Files/R/R-4.4.1/bin/Rscript.exe" run_temperature_dependent_epifilter.R

options(stringsAsFactors = FALSE)

# ---- User-visible settings -------------------------------------------------
city_case_file <- "df_city.csv"
district_case_file <- "District_Cases_Level.csv"
city_weather_file <- "广州市全市_2025年7-11月_气象数据.csv"
district_weather_file <- "广州市各区_2025年7-11月_气象数据_详细.csv"
output_dir <- "output_temperature_dependent_epifilter"
start_date_arg <- sub("^--start-date=", "", commandArgs(trailingOnly = TRUE)[grepl("^--start-date=", commandArgs(trailingOnly = TRUE))])
analysis_start_date <- if (length(start_date_arg) == 0L) as.Date(NA) else as.Date(start_date_arg[1L])
if (length(start_date_arg) > 1L || (!is.na(analysis_start_date) && is.na(as.numeric(analysis_start_date)))) stop("Use --start-date=YYYY-MM-DD.", call. = FALSE)
if (!is.na(analysis_start_date)) output_dir <- paste0(output_dir, "_start_", format(analysis_start_date, "%Y%m%d"))

# Observation model and EpiFilter state model.
imported_infectivity_weight <- 0.50
Rmin <- 0.01
Rmax <- 20.00
m <- 400L
eta <- 0.20
credibility <- 0.95

# Temperature-dependent human--mosquito--human generation interval.
max_generation_lag <- 60L
temperature_lower_anchor <- 18
temperature_upper_anchor <- 30
eip50_at_lower_anchor <- 8.7  # days, CHIKV/Ae. albopictus literature anchor
eip50_at_upper_anchor <- 3.7  # days, CHIKV/Ae. albopictus literature anchor
eip_gamma_shape <- 4          # spread around the temperature-specific median
tmh_mean_days <- 1            # mosquito-to-human waiting time, primary analysis

# Relative-onset-day human-to-mosquito infectiousness kernel.  The values are
# normalized; it is a transparent proxy in the absence of Guangzhou viremia
# measurements and is varied in future work when such data become available.
thm_days <- -1:6
thm_weights <- c(0.10, 0.20, 0.18, 0.16, 0.13, 0.10, 0.08, 0.05)
district_plot_labels <- c(
  "白云区" = "Baiyun", "从化区" = "Conghua", "番禺区" = "Panyu",
  "海珠区" = "Haizhu", "花都区" = "Huadu", "黄埔区" = "Huangpu",
  "荔湾区" = "Liwan", "南沙区" = "Nansha", "天河区" = "Tianhe",
  "越秀区" = "Yuexiu", "增城区" = "Zengcheng"
)

# ---- Paths and checks ------------------------------------------------------
script_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", script_args[grepl("^--file=", script_args)])
base_dir <- if (length(script_path) == 1L) dirname(normalizePath(script_path)) else getwd()
path_here <- function(path) file.path(base_dir, path)
stop_if <- function(condition, message) if (isTRUE(condition)) stop(message, call. = FALSE)
require_columns <- function(data, columns, label) {
  missing <- setdiff(columns, names(data))
  stop_if(length(missing) > 0L, paste0(label, " is missing column(s): ", paste(missing, collapse = ", ")))
}

# ---- Incidence and weather preparation ------------------------------------
complete_daily_series <- function(data, group_name = NULL) {
  require_columns(data, c("date", "local", "imported"), "Case input")
  data$date <- as.Date(data$date)
  data$local <- as.numeric(data$local)
  data$imported <- as.numeric(data$imported)
  stop_if(anyNA(data$date) || anyNA(data$local) || anyNA(data$imported), "Case dates and counts must be non-missing.")
  stop_if(any(data$local < 0 | data$imported < 0), "Case counts cannot be negative.")
  form <- if (is.null(group_name)) cbind(local, imported) ~ date else as.formula(paste("cbind(local, imported) ~", group_name, "+ date"))
  collapsed <- aggregate(form, data = data, FUN = sum)
  dates <- seq(min(collapsed$date), max(collapsed$date), by = "day")
  if (is.null(group_name)) {
    out <- merge(data.frame(date = dates), collapsed, by = "date", all.x = TRUE, sort = TRUE)
  } else {
    groups <- sort(unique(as.character(collapsed[[group_name]])))
    calendar <- expand.grid(group = groups, date = dates, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    names(calendar)[1L] <- group_name
    out <- merge(calendar, collapsed, by = c(group_name, "date"), all.x = TRUE, sort = TRUE)
  }
  out$local[is.na(out$local)] <- 0
  out$imported[is.na(out$imported)] <- 0
  out
}

prepare_weather <- function(weather, temperature_column, group_name = NULL) {
  required <- c("date", temperature_column)
  if (!is.null(group_name)) required <- c(required, group_name)
  require_columns(weather, required, "Weather input")
  out <- weather[, required, drop = FALSE]
  names(out)[names(out) == temperature_column] <- "temperature_raw"
  out$date <- as.Date(out$date)
  out$temperature_raw <- as.numeric(out$temperature_raw)
  stop_if(anyNA(out$date) || anyNA(out$temperature_raw), "Weather date and temperature must be non-missing.")
  key <- if (is.null(group_name)) "date" else c(group_name, "date")
  stop_if(anyDuplicated(out[, key, drop = FALSE]) > 0L, "Weather has duplicate area-date rows.")
  out$temperature_used <- pmin(pmax(out$temperature_raw, temperature_lower_anchor), temperature_upper_anchor)
  out
}

# ---- Temperature-dependent generation interval ---------------------------
discretise_distribution <- function(cdf_fun, max_day) {
  day <- 0:max_day
  p <- cdf_fun(day + 1) - cdf_fun(day)
  p[p < 0] <- 0
  p
}

convolve_pmf <- function(x, y) {
  out <- numeric(length(x) + length(y) - 1L)
  for (i in seq_along(x)) out[i - 1L + seq_along(y)] <- out[i - 1L + seq_along(y)] + x[i] * y
  out
}

eip50_from_temperature <- function(temp_c) {
  # Log-linear interpolation through the two published temperature anchors.
  beta <- log(eip50_at_upper_anchor / eip50_at_lower_anchor) /
    (temperature_upper_anchor - temperature_lower_anchor)
  eip50_at_lower_anchor * exp(beta * (temp_c - temperature_lower_anchor))
}

temperature_kernel <- function(temp_c, tmh_mean = tmh_mean_days, max_lag = max_generation_lag) {
  # I_H has median 5.4 days.  A 60-day internal horizon makes residual mass
  # negligible before the final max_lag truncation check.
  ih <- discretise_distribution(function(x) plnorm(x, meanlog = log(5.4), sdlog = 0.387), 60L)
  eip50 <- eip50_from_temperature(temp_c)
  eip_scale <- eip50 / qgamma(0.5, shape = eip_gamma_shape, scale = 1)
  eip <- discretise_distribution(function(x) pgamma(x, shape = eip_gamma_shape, scale = eip_scale), 60L)
  tmh <- discretise_distribution(function(x) pexp(x, rate = 1 / tmh_mean), 60L)

  # Support starts: I_H=0, T_HM=-1, EIP=0, T_MH=0.
  combined <- convolve_pmf(convolve_pmf(convolve_pmf(ih, thm_weights), eip), tmh)
  support <- -1:(-1L + length(combined) - 1L)
  keep <- support >= 1L & support <= max_lag
  kernel <- combined[keep]
  names(kernel) <- support[keep]
  captured_mass <- sum(kernel)
  stop_if(captured_mass <= 0, "Generation-interval kernel has zero mass at positive lags.")
  kernel <- kernel / captured_mass
  list(kernel = as.numeric(kernel), captured_mass = captured_mass, eip50_days = eip50)
}

build_kernel_matrix <- function(temperature, tmh_mean = tmh_mean_days, max_lag = max_generation_lag) {
  built <- lapply(temperature, temperature_kernel, tmh_mean = tmh_mean, max_lag = max_lag)
  list(
    weights = do.call(rbind, lapply(built, `[[`, "kernel")),
    captured_mass = vapply(built, `[[`, numeric(1), "captured_mass"),
    eip50_days = vapply(built, `[[`, numeric(1), "eip50_days")
  )
}

calculate_dynamic_infectiousness <- function(source_incidence, kernel_weights) {
  n <- length(source_incidence)
  max_lag <- ncol(kernel_weights)
  lambda <- numeric(n)
  for (t in seq_len(n)) {
    usable_lags <- seq_len(min(t - 1L, max_lag))
    if (length(usable_lags) > 0L) {
      source_rows <- t - usable_lags
      lambda[t] <- sum(source_incidence[source_rows] * kernel_weights[cbind(source_rows, usable_lags)])
    }
  }
  lambda
}

# ---- EpiFilter filter and smoother ----------------------------------------
weighted_quantile_grid <- function(grid, probabilities, probs = c(0.025, 0.5, 0.975)) {
  probabilities <- pmax(probabilities, 0)
  probabilities <- probabilities / sum(probabilities)
  cumulative <- cumsum(probabilities)
  vapply(probs, function(p) grid[which(cumulative >= p)[1L]], numeric(1))
}

make_transition_matrix <- function(grid, eta_value) {
  transition <- vapply(grid, function(previous) dnorm(grid, mean = previous, sd = sqrt(previous) * eta_value), numeric(length(grid)))
  transition <- t(transition)
  transition / rowSums(transition)
}

epi_filter <- function(incidence, infectiousness, grid, transition, prior) {
  n <- length(incidence); n_grid <- length(grid)
  posterior <- matrix(NA_real_, nrow = n, ncol = n_grid)
  predictive <- matrix(NA_real_, nrow = n, ncol = n_grid)
  posterior[1L, ] <- prior / sum(prior)
  predictive[1L, ] <- posterior[1L, ]
  if (n >= 2L) for (t in 2:n) {
    pred <- drop(posterior[t - 1L, ] %*% transition)
    loglik <- dpois(incidence[t], lambda = pmax(infectiousness[t] * grid, .Machine$double.xmin), log = TRUE)
    lik <- exp(loglik - max(loglik))
    unnormalised <- pred * lik
    posterior[t, ] <- if (sum(unnormalised) > 0 && is.finite(sum(unnormalised))) unnormalised / sum(unnormalised) else pred
    predictive[t, ] <- pred
  }
  list(posterior = posterior, predictive = predictive)
}

epi_smoother <- function(posterior, predictive, transition) {
  n <- nrow(posterior)
  smoothed <- posterior
  if (n >= 2L) for (t in (n - 1L):1L) {
    back <- smoothed[t + 1L, ] / pmax(predictive[t + 1L, ], .Machine$double.xmin)
    smoothed[t, ] <- posterior[t, ] * drop(transition %*% back)
    smoothed[t, ] <- smoothed[t, ] / sum(smoothed[t, ])
  }
  smoothed
}

summarise_posterior <- function(data, infectiousness, posterior, grid, estimate_type, area, kernel_info) {
  limits <- t(vapply(seq_len(nrow(posterior)), function(i) weighted_quantile_grid(grid, posterior[i, ], c((1 - credibility) / 2, 0.5, 1 - (1 - credibility) / 2)), numeric(3)))
  data.frame(
    date = data$date, area = area, estimate_type = estimate_type,
    local_cases = data$local, imported_cases = data$imported,
    temperature_raw_c = data$temperature_raw, temperature_used_c = data$temperature_used,
    eip50_days = kernel_info$eip50_days, generation_mass_within_truncation_days = kernel_info$captured_mass,
    total_infectiousness = infectiousness,
    rt_mean = drop(posterior %*% grid), rt_median = limits[, 2L],
    rt_lower = limits[, 1L], rt_upper = limits[, 3L],
    pr_rt_gt_1 = rowSums(posterior[, grid > 1, drop = FALSE]),
    stringsAsFactors = FALSE
  )
}

estimate_area <- function(data, area, grid, transition, prior, tmh_mean = tmh_mean_days, max_lag = max_generation_lag) {
  kernel_info <- build_kernel_matrix(data$temperature_used, tmh_mean = tmh_mean, max_lag = max_lag)
  source <- data$local + imported_infectivity_weight * data$imported
  infectiousness <- calculate_dynamic_infectiousness(source, kernel_info$weights)
  fit <- epi_filter(data$local, infectiousness, grid, transition, prior)
  smoother <- epi_smoother(fit$posterior, fit$predictive, transition)
  rbind(
    summarise_posterior(data, infectiousness, fit$posterior, grid, "filtered_realtime", area, kernel_info),
    summarise_posterior(data, infectiousness, smoother, grid, "smoothed_retrospective", area, kernel_info)
  )
}

plot_city_rt <- function(result, filename) {
  real_time <- result[result$estimate_type == "filtered_realtime", ]
  smooth <- result[result$estimate_type == "smoothed_retrospective", ]
  grDevices::pdf(filename, width = 10, height = 5.5)
  on.exit(grDevices::dev.off(), add = TRUE)
  ymax <- max(c(real_time$rt_upper, smooth$rt_upper, 1), na.rm = TRUE)
  plot(real_time$date, real_time$rt_median, type = "n", xaxt = "n", ylim = c(0, ymax * 1.05), xlab = "Date", ylab = expression(R[t]), main = "Guangzhou temperature-dependent EpiFilter Rt")
  axis.Date(1, at = seq(min(real_time$date), max(real_time$date), length.out = 6), format = "%Y-%m-%d")
  polygon(c(real_time$date, rev(real_time$date)), c(real_time$rt_lower, rev(real_time$rt_upper)), col = grDevices::adjustcolor("#0072B2", alpha.f = .18), border = NA)
  lines(real_time$date, real_time$rt_median, col = "#0072B2", lwd = 2)
  lines(smooth$date, smooth$rt_median, col = "#D55E00", lwd = 2, lty = 2)
  abline(h = 1, lty = 3, col = "grey35")
  legend("topright", c("Filtered real-time (95% CrI)", "Smoothed retrospective"), col = c("#0072B2", "#D55E00"), lty = c(1, 2), lwd = 2, bty = "n")
}

plot_district_panels <- function(result, filename) {
  result <- result[result$estimate_type == "smoothed_retrospective", ]
  areas <- sort(unique(result$area))
  grDevices::pdf(filename, width = 12, height = 12)
  on.exit(grDevices::dev.off(), add = TRUE)
  old_par <- par(mfrow = c(4, 3), mar = c(3.2, 3.3, 2.4, 0.7), oma = c(0, 0, 1.2, 0))
  on.exit(par(old_par), add = TRUE)
  for (area in areas) {
    x <- result[result$area == area, ]
    estimable <- x[x$total_infectiousness >= 1, ]
    ymax <- max(c(estimable$rt_upper, 1), na.rm = TRUE)
    plot(x$date, x$rt_median, type = "n", xaxt = "n", ylim = c(0, ymax * 1.08), xlab = "", ylab = expression(R[t]), main = unname(district_plot_labels[area]))
    axis.Date(1, at = seq(min(x$date), max(x$date), length.out = 4), format = "%m-%d", cex.axis = .75)
    polygon(c(x$date, rev(x$date)), c(x$rt_lower, rev(x$rt_upper)), col = grDevices::adjustcolor("#D55E00", alpha.f = .18), border = NA)
    lines(x$date, x$rt_median, col = "#D55E00", lwd = 1.5)
    abline(h = 1, lty = 3, col = "grey35")
  }
  mtext("District temperature-dependent EpiFilter Rt (smoothed retrospective; 95% CrI)", outer = TRUE, cex = 1.05)
}

summarise_areas <- function(result) {
  result <- result[result$estimate_type == "smoothed_retrospective", ]
  split_result <- split(result, result$area)
  out <- lapply(split_result, function(x) {
    # When Lambda_t is zero the posterior is driven by the diffuse initial
    # prior, so those dates must never be labelled as an epidemic peak.
    estimable <- x[x$total_infectiousness >= 1, ]
    if (nrow(estimable) == 0L) estimable <- x[x$total_infectiousness > 0, ]
    peak <- estimable[which.max(estimable$rt_median), ]
    data.frame(
      area = peak$area,
      local_cases = sum(x$local_cases), imported_cases = sum(x$imported_cases),
      peak_date = peak$date, peak_rt_median = peak$rt_median,
      peak_rt_lower = peak$rt_lower, peak_rt_upper = peak$rt_upper,
      peak_selection_min_infectiousness = 1,
      first_estimable_date = min(estimable$date),
      days_pr_rt_gt_1_ge_0_5 = sum(estimable$pr_rt_gt_1 >= 0.5),
      first_date_pr_rt_gt_1_ge_0_5 = if (any(estimable$pr_rt_gt_1 >= 0.5)) min(estimable$date[estimable$pr_rt_gt_1 >= 0.5]) else as.Date(NA),
      last_date_pr_rt_gt_1_ge_0_5 = if (any(estimable$pr_rt_gt_1 >= 0.5)) max(estimable$date[estimable$pr_rt_gt_1 >= 0.5]) else as.Date(NA),
      mean_temperature_c = mean(x$temperature_raw_c), mean_eip50_days = mean(x$eip50_days),
      min_generation_mass_within_truncation_days = min(x$generation_mass_within_truncation_days),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

# ---- Run -------------------------------------------------------------------
dir.create(path_here(output_dir), showWarnings = FALSE, recursive = TRUE)
city_cases <- read.csv(path_here(city_case_file), check.names = FALSE, fileEncoding = "UTF-8-BOM")
district_cases <- read.csv(path_here(district_case_file), check.names = FALSE, fileEncoding = "UTF-8-BOM")
city_weather <- read.csv(path_here(city_weather_file), check.names = FALSE, fileEncoding = "UTF-8-BOM")
district_weather <- read.csv(path_here(district_weather_file), check.names = FALSE, fileEncoding = "UTF-8-BOM")
require_columns(district_cases, "district", "District case input")
if (!is.na(analysis_start_date)) {
  city_cases$date <- as.Date(city_cases$date)
  district_cases$date <- as.Date(district_cases$date)
  city_cases <- city_cases[city_cases$date >= analysis_start_date, , drop = FALSE]
  district_cases <- district_cases[district_cases$date >= analysis_start_date, , drop = FALSE]
  stop_if(nrow(city_cases) == 0L || nrow(district_cases) == 0L, "No case records remain after start-date restriction.")
}

city_daily <- complete_daily_series(city_cases)
district_daily <- complete_daily_series(district_cases, "district")
city_weather <- prepare_weather(city_weather, "temperature")
district_weather <- prepare_weather(district_weather, "temp_mean", "district")
city_daily <- merge(city_daily, city_weather, by = "date", all.x = TRUE, sort = TRUE)
district_daily <- merge(district_daily, district_weather, by = c("district", "date"), all.x = TRUE, sort = TRUE)
stop_if(anyNA(city_daily$temperature_raw), "City case dates are not fully covered by city weather data.")
stop_if(anyNA(district_daily$temperature_raw), "District case dates are not fully covered by district weather data.")

grid <- seq(Rmin, Rmax, length.out = m)
transition <- make_transition_matrix(grid, eta)
prior <- rep(1 / m, m)
city_result <- estimate_area(city_daily, "City", grid, transition, prior)
district_result <- do.call(rbind, lapply(sort(unique(district_daily$district)), function(area) {
  estimate_area(district_daily[district_daily$district == area, c("date", "local", "imported", "temperature_raw", "temperature_used")], area, grid, transition, prior)
}))
all_result <- rbind(city_result, district_result)
area_summary <- summarise_areas(all_result)

# One-factor city-level sensitivity analyses requested in the analysis plan.
# The district-level primary analysis retains the full 60-day kernel to avoid
# unnecessary instability where counts are sparse.
sensitivity_specs <- rbind(
  data.frame(scenario = "T_MH mean = 1 day (primary)", tmh_mean_days = 1, max_generation_lag = 60),
  data.frame(scenario = "T_MH mean = 2 days", tmh_mean_days = 2, max_generation_lag = 60),
  data.frame(scenario = "T_MH mean = 3 days", tmh_mean_days = 3, max_generation_lag = 60),
  data.frame(scenario = "Truncation U = 28 days", tmh_mean_days = 1, max_generation_lag = 28),
  data.frame(scenario = "Truncation U = 35 days", tmh_mean_days = 1, max_generation_lag = 35),
  data.frame(scenario = "Truncation U = 42 days", tmh_mean_days = 1, max_generation_lag = 42)
)
sensitivity_fits <- lapply(seq_len(nrow(sensitivity_specs)), function(i) {
  spec <- sensitivity_specs[i, ]
  fit <- estimate_area(city_daily, "City", grid, transition, prior,
                       tmh_mean = spec$tmh_mean_days, max_lag = spec$max_generation_lag)
  fit$scenario <- spec$scenario
  fit$scenario_tmh_mean_days <- spec$tmh_mean_days
  fit$scenario_max_generation_lag <- spec$max_generation_lag
  fit
})
city_sensitivity_result <- do.call(rbind, sensitivity_fits)
city_sensitivity_summary <- do.call(rbind, lapply(seq_along(sensitivity_fits), function(i) {
  summary <- summarise_areas(sensitivity_fits[[i]])
  summary$scenario <- sensitivity_specs$scenario[i]
  summary$scenario_tmh_mean_days <- sensitivity_specs$tmh_mean_days[i]
  summary$scenario_max_generation_lag <- sensitivity_specs$max_generation_lag[i]
  summary
}))

write.csv(city_daily, path_here(file.path(output_dir, "city_daily_cases_weather.csv")), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(district_daily, path_here(file.path(output_dir, "district_daily_cases_weather.csv")), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(city_result, path_here(file.path(output_dir, "rt_city_temperature_dependent_epifilter.csv")), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(district_result, path_here(file.path(output_dir, "rt_district_temperature_dependent_epifilter.csv")), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(area_summary, path_here(file.path(output_dir, "rt_temperature_dependent_area_summary.csv")), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(city_sensitivity_result, path_here(file.path(output_dir, "rt_city_temperature_dependent_sensitivity.csv")), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(city_sensitivity_summary, path_here(file.path(output_dir, "rt_city_temperature_dependent_sensitivity_summary.csv")), row.names = FALSE, fileEncoding = "UTF-8")
plot_city_rt(city_result, path_here(file.path(output_dir, "rt_city_temperature_dependent_epifilter.pdf")))
plot_district_panels(district_result, path_here(file.path(output_dir, "rt_district_temperature_dependent_epifilter_panels.pdf")))

metadata <- c(
  "Method: temperature-dependent generation interval within discrete-grid EpiFilter and backward smoother.",
  "Generation interval: I_H + T_HM + I_M(T) + T_MH.",
  "I_H: Lognormal(meanlog=log(5.4), sdlog=0.387).",
  paste0("T_HM: discrete onset-relative kernel on days ", min(thm_days), " to ", max(thm_days), "; weights=", paste(thm_weights, collapse = ","), "."),
  paste0("I_M(T): Gamma(shape=", eip_gamma_shape, ") with median EIP interpolated log-linearly from ", eip50_at_lower_anchor, " d at ", temperature_lower_anchor, " C to ", eip50_at_upper_anchor, " d at ", temperature_upper_anchor, " C."),
  paste0("T_MH: Exponential(mean=", tmh_mean_days, " day)."),
  paste0("Temperature clipped only outside ", temperature_lower_anchor, "-", temperature_upper_anchor, " C; raw and used temperatures are both retained."),
  paste0("Generation-interval kernel truncated at ", max_generation_lag, " days and renormalised; retained mass reported per date."),
  paste0("Imported-case infectivity weight: ", imported_infectivity_weight, "."),
  paste0("Rt grid: ", Rmin, "-", Rmax, " (", m, " points); eta=", eta, "."),
  paste0("Date range: ", min(city_daily$date), " to ", max(city_daily$date), "."),
  paste0("Analysis start-date restriction: ", ifelse(is.na(analysis_start_date), "none", as.character(analysis_start_date)), "."),
  "filtered_realtime uses data through each date only; smoothed_retrospective uses the full epidemic series."
)
writeLines(metadata, path_here(file.path(output_dir, "run_metadata.txt")), useBytes = TRUE)
message("Temperature-dependent EpiFilter analysis completed: ", normalizePath(path_here(output_dir)))
