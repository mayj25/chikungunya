#!/usr/bin/env Rscript

# Serial-interval uncertainty sensitivity analysis for Guangzhou CHIKV Rt.
# This script is additive: it never overwrites files under Rt/.
#
# Primary estimand:
#   propagation of sampling uncertainty in the locally inferred serial
#   interval, conditional on the nine selected transmission pairs.
#
# The supplied workbook already contains the 1,000 fitted Gamma bootstrap
# distributions.  It is therefore used as the bootstrap input; the nine raw
# pair values are not reconstructed from the parameter workbook.

options(stringsAsFactors = FALSE, scipen = 999)
if (.Platform$OS.type == "windows") {
  invisible(suppressWarnings(Sys.setlocale("LC_ALL", "Chinese")))
}

required <- c("readxl", "EpiEstim", "ggplot2", "writexl")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))

suppressPackageStartupMessages({
  library(readxl)
  library(EpiEstim)
  library(ggplot2)
  library(writexl)
})

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
final_out_dir <- file.path(root, "Rt_SI_uncertainty_analysis_20260829")
# R's Windows native-encoding layer may reject writes below a Chinese path.
# All artifacts are therefore staged under an ASCII temp directory and copied
# to the requested workspace directory after the analysis completes.
out_dir <- file.path(tempdir(), "Rt_SI_uncertainty_analysis_20260829")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

bootstrap_file <- file.path(root, "bootstrap_gamma_parameters.xlsx")
city_file <- file.path(root, "Rt", "df_city.csv")
district_file <- file.path(root, "Rt", "District_Cases_Level.csv")
for (f in c(bootstrap_file, city_file, district_file)) {
  if (!file.exists(f)) stop("Required input not found: ", f)
}

seed <- 20260829L
alpha_imported <- 0.5
window_size <- 14L
prior_mean <- 5
prior_sd <- 5
n_boot_expected <- 1000L
n2_epiestim <- 50L

# EpiEstim convention: day 0 has zero serial-interval mass.  For k>=1,
# discretise with Gamma CDF interval probabilities P(k-1 <= SI < k), then
# renormalise after an upper tail chosen from the fitted distributions.
gamma_to_epiestim_pmf <- function(shape, scale, max_lag) {
  if (!is.finite(shape) || !is.finite(scale) || shape <= 0 || scale <= 0)
    stop("Invalid Gamma parameters")
  p <- numeric(max_lag + 1L)
  p[1L] <- 0
  p[-1L] <- pgamma(seq_len(max_lag), shape = shape, scale = scale) -
    pgamma(seq_len(max_lag) - 1, shape = shape, scale = scale)
  if (!all(is.finite(p)) || sum(p) <= 0) stop("Invalid SI PMF")
  p / sum(p)
}

read_bootstrap <- function(path) {
  # readxl's unzip helper on this Windows/R locale cannot open an xlsx whose
  # parent path contains Chinese characters; copy to the ASCII temp path while
  # retaining the original as the declared input.
  tmp_path <- file.path(tempdir(), "bootstrap_gamma_parameters_input.xlsx")
  if (!file.copy(path, tmp_path, overwrite = TRUE)) stop("Could not stage bootstrap workbook")
  x <- as.data.frame(read_excel(tmp_path, sheet = 1), check.names = FALSE)
  required <- c("bootstrap_id", "shape", "scale", "gamma_mean", "gamma_sd")
  if (!all(required %in% names(x))) stop("Bootstrap workbook must contain: ", paste(required, collapse = ", "))
  x$valid <- is.finite(x$shape) & is.finite(x$scale) & x$shape > 0 & x$scale > 0 &
    is.finite(x$gamma_mean) & is.finite(x$gamma_sd) & x$gamma_sd > 0
  x$failure_reason <- ifelse(x$valid, "", "non-finite/non-positive Gamma parameter")
  x
}

boot <- read_bootstrap(bootstrap_file)
if (nrow(boot) != n_boot_expected) warning("Workbook rows: ", nrow(boot), "; expected ", n_boot_expected)
valid_boot <- boot[boot$valid, , drop = FALSE]
if (!nrow(valid_boot)) stop("No valid bootstrap Gamma distribution")

# Select one common lag limit that retains >=99.9% of every valid fitted Gamma.
max_lag <- max(1L, ceiling(qgamma(0.999, valid_boot$shape, scale = valid_boot$scale)))
tail_before_norm <- pgamma(max_lag, valid_boot$shape, scale = valid_boot$scale)
if (any(tail_before_norm < 0.999)) stop("SI lag limit failed 0.999 tail rule")

pmf_boot <- sapply(seq_len(nrow(valid_boot)), function(i)
  gamma_to_epiestim_pmf(valid_boot$shape[i], valid_boot$scale[i], max_lag))
colnames(pmf_boot) <- sprintf("bootstrap_%04d", valid_boot$bootstrap_id)
if (any(abs(colSums(pmf_boot) - 1) > 1e-10)) stop("PMF columns do not sum to 1")

boot$pmf_included <- FALSE
boot$pmf_max_lag <- NA_integer_
boot$tail_probability_at_max_lag <- NA_real_
boot$pmf_included[boot$valid] <- TRUE
boot$pmf_max_lag[boot$valid] <- max_lag
boot$tail_probability_at_max_lag[boot$valid] <- tail_before_norm
write.csv(boot, file.path(out_dir, "si_bootstrap_parameters.csv"), row.names = FALSE, fileEncoding = "UTF-8")
pmf_df <- data.frame(day = 0:max_lag, pmf_boot, check.names = FALSE)
write.csv(pmf_df, file.path(out_dir, "si_bootstrap_pmf.csv"), row.names = FALSE, fileEncoding = "UTF-8")
saveRDS(pmf_boot, file.path(out_dir, "si_bootstrap_pmf.rds"))

fixed_local <- list(name = "fixed_local", label = "Fixed local SI", shape = (14.89 / 9.95)^2,
                    scale = 9.95^2 / 14.89, mean = 14.89, sd = 9.95)
riou <- list(name = "riou_external", label = "Riou et al. external SI", shape = (13.3 / 5.4)^2,
             scale = 5.4^2 / 13.3, mean = 13.3, sd = 5.4)
short_user <- list(name = "short_si_user_specified", label = "User-specified short SI", shape = (11.2 / 4.2)^2,
                   scale = 4.2^2 / 11.2, mean = 11.2, sd = 4.2)
scenario_info <- data.frame(
  scenario = c(fixed_local$name, "bootstrap_integrated", short_user$name, riou$name),
  label = c(fixed_local$label, "Bootstrap-integrated local SI", short_user$label, riou$label),
  mean_si = c(fixed_local$mean, mean(valid_boot$gamma_mean), short_user$mean, riou$mean),
  sd_si = c(fixed_local$sd, mean(valid_boot$gamma_sd), short_user$sd, riou$sd),
  stringsAsFactors = FALSE
)
write.csv(scenario_info, file.path(out_dir, "si_scenario_parameters.csv"), row.names = FALSE, fileEncoding = "UTF-8")
si_summary <- data.frame(
  quantity = c("raw_n", "fixed_gamma_shape", "fixed_gamma_scale", "fixed_gamma_mean", "fixed_gamma_sd",
               "fixed_gamma_central_95_lower", "fixed_gamma_central_95_upper",
               "bootstrap_valid_n", "bootstrap_mean_ci_lower", "bootstrap_mean_ci_upper",
               "bootstrap_sd_ci_lower", "bootstrap_sd_ci_upper"),
  estimate = c(NA, fixed_local$shape, fixed_local$scale, fixed_local$mean, fixed_local$sd,
               qgamma(0.025, fixed_local$shape, scale = fixed_local$scale),
               qgamma(0.975, fixed_local$shape, scale = fixed_local$scale), nrow(valid_boot),
               quantile(valid_boot$gamma_mean, 0.025), quantile(valid_boot$gamma_mean, 0.975),
               quantile(valid_boot$gamma_sd, 0.025), quantile(valid_boot$gamma_sd, 0.975)),
  note = c("Raw nine SI values are not present in the supplied workbook", "", "", "", "",
           "central 95% range of fitted Gamma distribution", "central 95% range of fitted Gamma distribution", "", "bootstrap percentile interval", "bootstrap percentile interval", "bootstrap percentile interval", "bootstrap percentile interval"),
  stringsAsFactors = FALSE
)
write.csv(si_summary, file.path(out_dir, "si_summary_statistics.csv"), row.names = FALSE, fileEncoding = "UTF-8")

read_cases <- function(path, district = FALSE) {
  x <- read.csv(path, check.names = FALSE, fileEncoding = "UTF-8-BOM")
  need <- if (district) c("district", "date", "local", "imported") else c("date", "local", "imported")
  if (!all(need %in% names(x))) stop("Case file missing: ", paste(setdiff(need, names(x)), collapse = ", "))
  x$date <- as.Date(gsub("/", "-", as.character(x$date)))
  x$local <- as.numeric(x$local)
  x$imported <- as.numeric(x$imported)
  if (anyNA(x[, need, drop = FALSE])) stop("Missing required case values in ", path)
  if (any(x$local < 0 | x$imported < 0)) stop("Negative case counts in ", path)
  if (district) {
    x$district <- as.character(x$district)
    x <- x[order(x$district, seq_len(nrow(x))), , drop = FALSE]
  } else x <- x[order(x$date), , drop = FALSE]
  x
}

city <- read_cases(city_file, FALSE)
district <- read_cases(district_file, TRUE)

# Preserve the existing analysis' row sequence.  The source district file has
# 125 rows over a 142-day calendar span; dates are retained only as labels.
make_incid <- function(x) {
  z <- data.frame(local = x$local, imported = x$imported * alpha_imported)
  if (z$local[1] > 0) { z$imported[1] <- z$imported[1] + z$local[1]; z$local[1] <- 0 }
  z
}

posterior_probability_gt1 <- function(incid, pmf, t_start, t_end) {
  lambda <- EpiEstim:::overall_infectivity(incid, pmf)
  a0 <- (prior_mean / prior_sd)^2
  b0 <- prior_sd^2 / prior_mean
  vapply(seq_along(t_end), function(j) {
    if (t_end[j] <= sum(pmf * (0:max_lag))) return(NA_real_)
    a <- a0 + sum(incid$local[t_start[j]:t_end[j]])
    b <- 1 / (1 / b0 + sum(lambda[t_start[j]:t_end[j]], na.rm = TRUE))
    1 - pgamma(1, shape = a, scale = b)
  }, numeric(1))
}

posterior_probability_integrated <- function(incid, pmf_matrix, t_start, t_end) {
  all_p <- vapply(seq_len(ncol(pmf_matrix)), function(k)
    posterior_probability_gt1(incid, pmf_matrix[, k], t_start, t_end), numeric(length(t_end)))
  vapply(seq_along(t_end), function(j) mean(all_p[j, is.finite(all_p[j, ])]), numeric(1))
}

run_one <- function(x, area, scenario, pmf, integrated = FALSE, seed_value = 1L) {
  incid <- make_incid(x)
  n <- nrow(incid)
  t_start <- 2:(n - window_size + 1L)
  t_end <- t_start + window_size - 1L
  cfg_args <- list(t_start = t_start, t_end = t_end, mean_prior = prior_mean,
                   std_prior = prior_sd, seed = seed_value)
  if (integrated) {
    cfg_args$n2 <- n2_epiestim
    fit <- estimate_R(incid, method = "si_from_sample", si_sample = pmf,
                      config = make_config(cfg_args))
    pgt <- posterior_probability_integrated(incid, pmf, t_start, t_end)
  } else {
    fit <- estimate_R(incid, method = "non_parametric_si", config =
      make_config(c(cfg_args, list(si_distr = pmf))))
    pgt <- posterior_probability_gt1(incid, pmf, t_start, t_end)
  }
  r <- fit$R
  out <- data.frame(
    area = area, scenario = scenario, time_index = r$t_end,
    date = x$date[r$t_end], window_start = x$date[r$t_start], window_end = x$date[r$t_end],
    rt_mean = r$`Mean(R)`, rt_median = r$`Median(R)`, rt_sd = r$`Std(R)`,
    rt_q025 = r$`Quantile.0.025(R)`, rt_q975 = r$`Quantile.0.975(R)`,
    cri_width = r$`Quantile.0.975(R)` - r$`Quantile.0.025(R)`,
    posterior_probability_rt_gt1 = pgt,
    rt_class = ifelse(r$`Quantile.0.025(R)` > 1, "clearly_above_1",
      ifelse(r$`Quantile.0.975(R)` < 1, "clearly_below_1", "indeterminate")),
    stringsAsFactors = FALSE
  )
  out
}

run_area <- function(x, area) {
  pmf_fixed <- gamma_to_epiestim_pmf(fixed_local$shape, fixed_local$scale, max_lag)
  pmf_short <- gamma_to_epiestim_pmf(short_user$shape, short_user$scale, max_lag)
  pmf_riou <- gamma_to_epiestim_pmf(riou$shape, riou$scale, max_lag)
  do.call(rbind, list(
    run_one(x, area, fixed_local$name, pmf_fixed, FALSE, seed + 1L),
    run_one(x, area, "bootstrap_integrated", pmf_boot, TRUE, seed + 2L),
    run_one(x, area, short_user$name, pmf_short, FALSE, seed + 3L),
    run_one(x, area, riou$name, pmf_riou, FALSE, seed + 4L)
  ))
}

rt_city <- run_area(city, "全市")
rt_district <- do.call(rbind, lapply(unique(district$district), function(d)
  run_area(district[district$district == d, , drop = FALSE], d)))
write.csv(rt_city, file.path(out_dir, "rt_citywide_all_SI_scenarios.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(rt_district, file.path(out_dir, "rt_district_all_SI_scenarios.csv"), row.names = FALSE, fileEncoding = "UTF-8")

safe_cor <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 3 || length(unique(a[ok])) < 2 || length(unique(b[ok])) < 2) NA_real_
  else suppressWarnings(cor(a[ok], b[ok], method = "spearman"))
}
first_sustained_below1 <- function(z, k = 7L) {
  z <- z[order(z$date), ]
  ok <- is.finite(z$rt_median) & z$rt_median < 1
  for (i in which(ok)) {
    idx <- i:min(nrow(z), i + k - 1L)
    if (length(idx) == k && all(ok[idx]) && all(diff(z$date[idx]) == 1)) return(z$date[i])
  }
  as.Date(NA)
}
first_clearly_below1 <- function(z) {
  d <- z$date[is.finite(z$rt_q975) & z$rt_q975 < 1]
  if (length(d)) min(d) else as.Date(NA)
}
two_peaks <- function(z, min_gap_days = 21L) {
  z <- z[is.finite(z$rt_median), , drop = FALSE]
  if (!nrow(z)) return(data.frame(peak_rank = 1:2, peak_date = as.Date(NA), peak_rt = NA_real_))
  z <- z[order(-z$rt_median), ]
  chosen <- integer()
  for (i in seq_len(nrow(z))) {
    if (!length(chosen) || all(abs(as.numeric(z$date[i] - z$date[chosen])) >= min_gap_days)) chosen <- c(chosen, i)
    if (length(chosen) == 2) break
  }
  data.frame(peak_rank = 1:2,
             peak_date = as.Date(c(z$date[chosen], rep(NA, 2 - length(chosen)))),
             peak_rt = c(z$rt_median[chosen], rep(NA_real_, 2 - length(chosen))))
}

compare_one <- function(dat, area, scenario) {
  f <- dat[dat$area == area & dat$scenario == fixed_local$name, ]
  s <- dat[dat$area == area & dat$scenario == scenario, ]
  m <- merge(f, s, by = "time_index", suffixes = c("_fixed", "_scenario"))
  m <- m[is.finite(m$rt_median_fixed) & is.finite(m$rt_median_scenario), ]
  width_pct <- 100 * (m$cri_width_scenario - m$cri_width_fixed) / m$cri_width_fixed
  cat_fixed <- m$rt_class_fixed; cat_scenario <- m$rt_class_scenario
  peak_f <- two_peaks(f); peak_s <- two_peaks(s)
  data.frame(
    area = area, scenario = scenario, n_paired_windows = nrow(m),
    spearman_median_rt = safe_cor(m$rt_median_fixed, m$rt_median_scenario),
    median_absolute_difference_rt = if (nrow(m)) median(abs(m$rt_median_scenario - m$rt_median_fixed), na.rm = TRUE) else NA,
    median_cri_width_scenario = if (nrow(m)) median(m$cri_width_scenario, na.rm = TRUE) else NA,
    median_percentage_increase_cri_width = if (nrow(m)) median(width_pct, na.rm = TRUE) else NA,
    rt_class_concordance = if (nrow(m)) mean(cat_fixed == cat_scenario) else NA,
    fixed_peak1_date = peak_f$peak_date[1], fixed_peak1_rt = peak_f$peak_rt[1],
    scenario_peak1_date = peak_s$peak_date[1], scenario_peak1_rt = peak_s$peak_rt[1],
    fixed_peak2_date = peak_f$peak_date[2], fixed_peak2_rt = peak_f$peak_rt[2],
    scenario_peak2_date = peak_s$peak_date[2], scenario_peak2_rt = peak_s$peak_rt[2],
    fixed_first_sustained_below1 = first_sustained_below1(f),
    scenario_first_sustained_below1 = first_sustained_below1(s),
    fixed_first_clearly_below1 = first_clearly_below1(f),
    scenario_first_clearly_below1 = first_clearly_below1(s),
    stringsAsFactors = FALSE
  )
}

scenarios_to_compare <- c("bootstrap_integrated", short_user$name, riou$name)
city_comp <- do.call(rbind, lapply(scenarios_to_compare, function(s) compare_one(rt_city, "全市", s)))
district_comp <- do.call(rbind, lapply(unique(rt_district$area), function(a)
  do.call(rbind, lapply(scenarios_to_compare, function(s) compare_one(rt_district, a, s)))))
write.csv(city_comp, file.path(out_dir, "rt_SI_comparison_summary_citywide.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(district_comp, file.path(out_dir, "rt_SI_comparison_summary_district.csv"), row.names = FALSE, fileEncoding = "UTF-8")

district_aggregate <- do.call(rbind, lapply(scenarios_to_compare, function(s) {
  q <- district_comp[district_comp$scenario == s, ]
  data.frame(scenario = s, n_districts = nrow(q),
    median_spearman = median(q$spearman_median_rt, na.rm = TRUE),
    range_spearman_min = min(q$spearman_median_rt, na.rm = TRUE), range_spearman_max = max(q$spearman_median_rt, na.rm = TRUE),
    median_cri_width_increase_pct = median(q$median_percentage_increase_cri_width, na.rm = TRUE),
    range_cri_width_increase_pct_min = min(q$median_percentage_increase_cri_width, na.rm = TRUE),
    range_cri_width_increase_pct_max = max(q$median_percentage_increase_cri_width, na.rm = TRUE),
    median_class_concordance = median(q$rt_class_concordance, na.rm = TRUE),
    range_class_concordance_min = min(q$rt_class_concordance, na.rm = TRUE), range_class_concordance_max = max(q$rt_class_concordance, na.rm = TRUE),
    stringsAsFactors = FALSE)
}))
write.csv(district_aggregate, file.path(out_dir, "rt_SI_comparison_summary_district_aggregate.csv"), row.names = FALSE, fileEncoding = "UTF-8")

pair_status <- data.frame(
  status = "not_available",
  reason = "The supplied bootstrap_gamma_parameters.xlsx contains Gamma parameters only; it has no source/recipient IDs, onset dates, posterior transmission probabilities, same-street/town fields, sequence IDs, or genetic distances. No pair-level table was fabricated.",
  stringsAsFactors = FALSE
)
write_xlsx(list(status = pair_status), file.path(out_dir, "Table_SI_pair_evidence.xlsx"))

# Figures: fitted Gamma curves and Rt medians/95% CrI.
grid <- seq(0, min(max_lag, 160), by = 0.1)
curve_list <- lapply(seq_len(nrow(valid_boot)), function(i)
  data.frame(si = grid, density = dgamma(grid, valid_boot$shape[i], scale = valid_boot$scale[i]),
             distribution = "Bootstrap local SI", id = valid_boot$bootstrap_id[i]))
curve_list <- c(curve_list,
  list(data.frame(si = grid, density = dgamma(grid, fixed_local$shape, scale = fixed_local$scale), distribution = fixed_local$label, id = NA),
       data.frame(si = grid, density = dgamma(grid, short_user$shape, scale = short_user$scale), distribution = short_user$label, id = NA),
       data.frame(si = grid, density = dgamma(grid, riou$shape, scale = riou$scale), distribution = riou$label, id = NA)))
curves <- do.call(rbind, curve_list)
p_si <- ggplot() + geom_line(data = curves[curves$distribution == "Bootstrap local SI", ], aes(si, density, group = id), colour = "grey60", alpha = 0.12) +
  geom_line(data = curves[curves$distribution != "Bootstrap local SI", ], aes(si, density, colour = distribution), linewidth = 0.9) +
  theme_classic(base_size = 12) + labs(x = "Serial interval (days)", y = "Density", colour = NULL) +
  scale_x_continuous(limits = c(0, min(max_lag, 100)))
ggsave(file.path(out_dir, "Figure_SI_observations_bootstrap_and_external_distributions.pdf"), p_si, width = 9, height = 6)
ggsave(file.path(out_dir, "Figure_SI_observations_bootstrap_and_external_distributions.png"), p_si, width = 9, height = 6, dpi = 300)

plot_rt <- function(dat, areas, path_stem) {
  z <- dat[dat$area %in% areas, ]
  p <- ggplot(z, aes(date, rt_median, colour = scenario, fill = scenario)) +
    geom_ribbon(aes(ymin = rt_q025, ymax = rt_q975), alpha = 0.10, colour = NA) +
    geom_line(linewidth = 0.45) + geom_hline(yintercept = 1, linetype = 2) +
    facet_wrap(~area, scales = "free_y") + theme_classic(base_size = 10) +
    labs(x = "Date label (source observation rows)", y = "Rt median (95% CrI)", colour = "SI scenario", fill = "SI scenario")
  ggsave(paste0(path_stem, ".pdf"), p, width = 12, height = if (length(areas) > 1) 9 else 6)
  ggsave(paste0(path_stem, ".png"), p, width = 12, height = if (length(areas) > 1) 9 else 6, dpi = 300)
}
plot_rt(rt_city, "全市", file.path(out_dir, "Figure_Rt_citywide_SI_sensitivity"))
plot_rt(rt_district, unique(rt_district$area), file.path(out_dir, "Figure_Rt_district_SI_sensitivity"))

# The source workbook contains no pair-level identifiers/onset dates, and the
# source case files contain no BI/MOI/intervention join keys sufficient to
# rerun downstream models without importing additional, separately specified
# analysis data.  Record this explicitly rather than fabricating a workbook.
downstream_status <- data.frame(
  analysis = c("BI/MOI--Rt lagged correlation", "intervention--Rt mixed model", "intervention--BI--Rt mediation", "pair-level SI evidence table"),
  status = "not_run",
  reason = c("No unambiguous BI/MOI join specification supplied in these inputs", "Existing model script uses separate derived panels; not altered here", "Existing mediation workflow requires its prespecified panel and model objects", "bootstrap parameter workbook has no source/recipient IDs, onset dates, posterior probabilities, or sequence IDs"),
  stringsAsFactors = FALSE
)
write_xlsx(list(status = downstream_status), file.path(out_dir, "downstream_analysis_SI_sensitivity.xlsx"))

boot_city <- city_comp[city_comp$scenario == "bootstrap_integrated", ]
short_city <- city_comp[city_comp$scenario == short_user$name, ]
summary_lines <- c(
  "Results summary (Chinese translation is provided in the final hand-off)",
  sprintf("The workbook supplied %d valid Gamma bootstrap distributions; each was discretised by Gamma CDF intervals into a daily PMF with a common %d-day upper limit and column sums of 1.", nrow(valid_boot), max_lag),
  sprintf("Citywide bootstrap-integrated versus fixed-local Rt median Spearman correlation was %.4f; median 95%% CrI width increased by %.1f%%; Rt class concordance was %.1f%%.", boot_city$spearman_median_rt, boot_city$median_percentage_increase_cri_width, 100 * boot_city$rt_class_concordance),
  sprintf("The user-specified short-SI scenario (mean 11.2 days; SD 4.2 days) gave correlation %.4f, median CrI-width change %.1f%%, and class concordance %.1f%%.", short_city$spearman_median_rt, short_city$median_percentage_increase_cri_width, 100 * short_city$rt_class_concordance),
  "These results quantify SI sampling uncertainty conditional on the nine selected transmission pairs and do not fully propagate transmission-pair reconstruction, under-ascertainment, or onset-date uncertainty.",
  "",
  "Methods (English)",
  "We propagated sampling uncertainty in the locally informed serial-interval distribution conditional on the nine selected transmission pairs. For each of 1,000 supplied bootstrap Gamma fits, we calculated daily serial-interval probabilities from Gamma CDF interval probabilities, set day 0 to zero according to the EpiEstim convention, renormalised after retaining at least 99.9% of the fitted distribution, and supplied the resulting PMF matrix to EpiEstim (method=si_from_sample). Incidence, the 14-observation sliding window, the Gamma prior with mean 5 and SD 5, imported-case relative infectiousness (0.5), and the source row sequence were held constant.",
  "",
  "Results (English)",
  "Across serial-interval assumptions, the absolute magnitude and posterior uncertainty of Rt varied. The bootstrap-integrated local-SI analysis retained the principal temporal pattern relative to the fixed local-SI analysis while widening the citywide credible intervals; the fixed shorter and external scenarios produced additional changes in Rt magnitude and threshold timing. These are descriptive sensitivity comparisons rather than evidence that any serial-interval scenario is the true distribution.",
  "",
  "Limitations (English)",
  "The local serial-interval distribution was informed by only nine selected high-confidence transmission pairs. The bootstrap analysis quantifies sampling uncertainty conditional on those pairs and does not fully propagate uncertainty in outbreaker2 reconstruction, case ascertainment, onset-date error, spatial transmission, or temporal variation in the serial interval. The external scenarios were derived under locations, temperatures, vector ecology, and modelling assumptions that may not represent Guangzhou.",
  "",
  "Response to reviewer (English)",
  "We thank the reviewer for noting that estimating the serial interval from nine inferred transmission pairs could understate uncertainty in Rt. We therefore added a conditional sampling-uncertainty analysis: all 1,000 bootstrap Gamma serial-interval distributions were discretised as EpiEstim-compatible PMFs and propagated through citywide and district-level Rt estimation using si_from_sample, with the incidence series and all other Rt settings unchanged. We report posterior medians, 95% credible intervals, interval widths, Rt>1 classifications, peak summaries, and comparisons with the fixed local, user-specified shorter, and external SI scenarios."
)
writeLines(summary_lines, file.path(out_dir, "results_summary_cn_and_writing_text.txt"), useBytes = TRUE)

readme <- c(
  "Serial-interval uncertainty Rt sensitivity analysis (additive output)",
  paste("Run date:", as.character(Sys.Date())), paste("R:", R.version.string), paste("EpiEstim:", as.character(packageVersion("EpiEstim"))),
  paste("Random seed:", seed), paste("Bootstrap rows read:", nrow(boot)), paste("Valid Gamma rows:", nrow(valid_boot)),
  paste("Common SI PMF maximum lag:", max_lag, "days; minimum retained Gamma CDF:", signif(min(tail_before_norm), 8)),
  paste("Window:", window_size, "observation rows; prior mean/SD:", prior_mean, "/", prior_sd), paste("Imported relative infectiousness:", alpha_imported),
  "The EpiEstim si_from_sample analysis uses all valid bootstrap PMF columns and n2=50 posterior draws per PMF.",
  "Case rows are not calendar-completed: both source case files contain 125 observations, while the date span is 142 days. Dates are labels; time_index is the sequence used by the existing Rt script.",
  "Fixed local SI: mean=14.89, SD=9.95 days. Riou et al.: mean=13.3, SD=5.4 days; literature model is temperature-dependent and the 25 C interpretation is not Guangzhou-specific.",
  "User-specified short-SI scenario: mean=11.2 days, SD=4.2 days; this is a sensitivity scenario without a supplied external citation.",
  "Interpretation: bootstrap-integrated results quantify sampling uncertainty conditional on the nine selected transmission pairs; they do not propagate outbreaker2 reconstruction, onset-date, reporting, or pair-selection uncertainty.",
  "The pair-evidence table and downstream analyses are marked not_run because the supplied parameter workbook lacks pair-level evidence and no unambiguous downstream join specification was provided.",
  paste("Output directory:", final_out_dir)
)
writeLines(readme, file.path(out_dir, "README_methods_and_reproducibility.txt"), useBytes = TRUE)

dir.create(final_out_dir, recursive = TRUE, showWarnings = FALSE)
staged_files <- list.files(out_dir, full.names = TRUE, recursive = FALSE)
if (length(staged_files)) file.copy(staged_files, final_out_dir, overwrite = TRUE)

cat("Completed additive SI sensitivity analysis. Outputs: ", final_out_dir, "\n", sep = "")
