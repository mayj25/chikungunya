#!/usr/bin/env Rscript

# Reproducible intervention, mosquito-index, Rt, and mediation analyses
# Guangzhou, 2025
#
# Analysis order
#   1. Primary analysis: current-week temperature + relative humidity
#   2. Mediation analysis: current-week temperature + relative humidity
#   3. Stratified analysis: central urban versus other districts
#   4. Sensitivity analysis: current-week temperature + total rainfall
#
# The mediation analysis uses district random-intercept models only. District
# fixed-effect mediation models are intentionally excluded from this submission.

options(stringsAsFactors = FALSE, digits = 10, scipen = 999)
if (.Platform$OS.type == "windows") {
  suppressWarnings(try(Sys.setlocale("LC_ALL", "Chinese"), silent = TRUE))
}

required_packages <- c("data.table", "lme4", "lmerTest", "mediation", "writexl")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Missing required R packages: ", paste(missing_packages, collapse = ", "),
    ". Install them before running this script."
  )
}
suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
  library(lmerTest)
  library(mediation)
})

# -----------------------------------------------------------------------------
# 0. Prespecified settings, paths, and reusable functions
# -----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
N_RT_DRAWS <- if (length(args) >= 1L) as.integer(args[[1L]]) else 10000L
N_BETA_DRAWS <- if (length(args) >= 2L) as.integer(args[[2L]]) else 10L
N_WORKERS_REQUESTED <- if (length(args) >= 3L) as.integer(args[[3L]]) else 4L
N_MEDIATION_SIMS <- if (length(args) >= 4L) as.integer(args[[4L]]) else 10000L

if (!is.finite(N_RT_DRAWS) || N_RT_DRAWS < 100L) {
  stop("N_RT_DRAWS must be at least 100")
}
if (!is.finite(N_BETA_DRAWS) || N_BETA_DRAWS < 1L) {
  stop("N_BETA_DRAWS must be positive")
}
if (!is.finite(N_WORKERS_REQUESTED) || N_WORKERS_REQUESTED < 1L) {
  stop("N_WORKERS_REQUESTED must be positive")
}
if (!is.finite(N_MEDIATION_SIMS) || N_MEDIATION_SIMS < 1000L) {
  stop("N_MEDIATION_SIMS must be at least 1000")
}

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
SCRIPT_DIR <- if (length(script_arg) == 1L) {
  dirname(sub("^--file=", "", script_arg))
} else {
  "."
}
DATA_DIR <- file.path(SCRIPT_DIR, "data")
OUTPUT_DIR <- file.path(SCRIPT_DIR, "output")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

RT_OFFSET <- 0.1
BI_OFFSET <- 1
OUTER_RT_SEED <- 20250901L
SECONDARY_BETA_SEED <- 20260829L
MEDIATION_SEED_BASE <- 20260829L
SAVE_INTERMEDIATE_DRAWS <- identical(
  tolower(Sys.getenv("SAVE_INTERMEDIATE_DRAWS", unset = "false")), "true"
)

central_districts <- c(
  "\u8d8a\u79c0\u533a", "\u8354\u6e7e\u533a",
  "\u6d77\u73e0\u533a", "\u5929\u6cb3\u533a"
)
interventions <- data.frame(
  Intervention_index = 1:3,
  Intervention = c(
    "Routine intervention", "Enhanced intervention", "Dynamic intervention"
  ),
  Rt_exposure = c("intervention_2", "enhance_2", "dynamic_2"),
  BI_exposure = c(
    "intervention_last_week", "enhance_last_week", "dynamic_last_week"
  ),
  stringsAsFactors = FALSE
)

input_files <- c(
  panel = file.path(DATA_DIR, "data33.csv"),
  daily_rt = file.path(DATA_DIR, "daily_rt_posterior_parameters_reconstructed.csv"),
  weekly_rt_reference = file.path(DATA_DIR, "weekly_rt_uncertainty_10000_draws.csv"),
  weather = file.path(DATA_DIR, "guangzhou_district_weekly_weather_2025_07_11.csv")
)
if (any(!file.exists(input_files))) {
  stop(
    "Missing input file(s): ",
    paste(input_files[!file.exists(input_files)], collapse = "; ")
  )
}

message("[Setup] Reading and validating the 2025 district-week panel")
panel <- read.csv(
  input_files[["panel"]], check.names = FALSE, fileEncoding = "UTF-8"
)
panel <- panel[, nzchar(names(panel)), drop = FALSE]
required_columns <- c(
  "district", "week", "weekly_Rt", "BI_this_week", "BI_last_week",
  "temp_weekly_avg_lag0", "rh_weekly_avg_lag0",
  interventions$Rt_exposure, interventions$BI_exposure
)
missing_columns <- setdiff(required_columns, names(panel))
if (length(missing_columns) > 0L) {
  stop("Missing panel columns: ", paste(missing_columns, collapse = ", "))
}
if (nrow(panel) != 182L || anyDuplicated(panel[c("district", "week")])) {
  stop("Unexpected district-week panel structure")
}
if (
  anyNA(panel[required_columns]) || any(panel$weekly_Rt <= 0) ||
  any(panel$BI_this_week < 0) || any(panel$BI_last_week < 0)
) {
  stop("Invalid or missing required panel values")
}
for (exposure in c(interventions$Rt_exposure, interventions$BI_exposure)) {
  values <- sort(unique(as.numeric(panel[[exposure]])))
  if (!identical(values, c(0, 1))) {
    stop(exposure, " must contain exactly 0 and 1")
  }
  panel[[exposure]] <- as.integer(panel[[exposure]])
}

weather <- data.table::fread(
  input_files[["weather"]], encoding = "UTF-8", data.table = FALSE
)
weather <- weather[weather$iso_year == 2025, c(
  "district_name", "iso_week", "total_precipitation_mm"
)]
weather <- unique(weather)
if (anyDuplicated(weather[c("district_name", "iso_week")])) {
  stop("Duplicate district-week keys in weather data")
}
panel_key <- paste(panel$district, panel$week, sep = "||")
weather_key <- paste(weather$district_name, weather$iso_week, sep = "||")
weather_index <- match(panel_key, weather_key)
if (anyNA(weather_index)) stop("Rainfall matching failed")
panel$total_precipitation_mm <- as.numeric(
  weather$total_precipitation_mm[weather_index]
)
if (
  any(!is.finite(panel$total_precipitation_mm)) ||
  any(panel$total_precipitation_mm < 0)
) {
  stop("Invalid total-rainfall values")
}

panel$district <- factor(panel$district)
panel$week <- as.integer(panel$week)
panel$week_centered <- panel$week - mean(panel$week)
panel$Stratum_code <- ifelse(
  as.character(panel$district) %in% central_districts, "central", "other"
)
panel$log_BI_current <- log(panel$BI_this_week + BI_OFFSET)
panel$log_BI_previous <- log(panel$BI_last_week + BI_OFFSET)
panel$log_Rt_plus0p1 <- log(panel$weekly_Rt + RT_OFFSET)

subset_panel <- function(group_code) {
  d <- if (group_code == "citywide") {
    panel
  } else {
    panel[panel$Stratum_code == group_code, , drop = FALSE]
  }
  droplevels(d)
}

group_labels <- c(
  citywide = "All Guangzhou districts",
  central = "Central urban districts",
  other = "Other districts"
)

lmer_control <- lme4::lmerControl(
  optimizer = "bobyqa",
  optCtrl = list(maxfun = 100000),
  check.conv.singular = lme4::.makeCC(action = "ignore", tol = 1e-4)
)

write_csv <- function(x, filename) {
  data.table::fwrite(
    data.table::as.data.table(x), file.path(OUTPUT_DIR, filename), bom = TRUE
  )
}

reconstruct_weekly_rt <- function(n_draws) {
  message("[Setup] Reconstructing ", n_draws, " weekly Rt posterior trajectories")
  daily_posterior <- data.table::fread(
    input_files[["daily_rt"]], encoding = "UTF-8", data.table = TRUE
  )
  daily_posterior[, valid := as.logical(valid)]
  valid_daily <- daily_posterior[
    valid == TRUE & is.finite(posterior_shape) & is.finite(posterior_scale)
  ]
  data.table::setorder(valid_daily, district, date)
  district_names <- unique(valid_daily$district)

  set.seed(OUTER_RT_SEED)
  weekly_key_list <- vector("list", length(district_names))
  weekly_matrix_list <- vector("list", length(district_names))
  for (i in seq_along(district_names)) {
    district_name <- district_names[i]
    district_rt <- valid_daily[district == district_name]
    n_days <- nrow(district_rt)
    daily_draws <- matrix(
      stats::rgamma(
        n_days * n_draws,
        shape = rep(district_rt$posterior_shape, times = n_draws),
        scale = rep(district_rt$posterior_scale, times = n_draws)
      ),
      nrow = n_days, ncol = n_draws
    )
    district_weeks <- sort(unique(district_rt$week))
    weekly_draws <- matrix(
      NA_real_, nrow = length(district_weeks), ncol = n_draws
    )
    for (j in seq_along(district_weeks)) {
      rows <- which(district_rt$week == district_weeks[j])
      weekly_draws[j, ] <- exp(colMeans(log(
        daily_draws[rows, , drop = FALSE]
      )))
    }
    weekly_key_list[[i]] <- data.frame(
      district = district_name, week = district_weeks,
      stringsAsFactors = FALSE
    )
    weekly_matrix_list[[i]] <- weekly_draws
  }

  weekly_key <- do.call(rbind, weekly_key_list)
  weekly_matrix <- do.call(rbind, weekly_matrix_list)
  weekly_id <- paste(weekly_key$district, weekly_key$week, sep = "||")
  panel_id <- paste(as.character(panel$district), panel$week, sep = "||")
  weekly_index <- match(panel_id, weekly_id)
  if (anyNA(weekly_index)) stop("Weekly Rt posterior matching failed")
  panel_weekly_draws <- weekly_matrix[weekly_index, , drop = FALSE]
  if (!identical(dim(panel_weekly_draws), c(nrow(panel), n_draws))) {
    stop("Unexpected weekly Rt posterior matrix dimensions")
  }

  validation <- data.frame()
  if (n_draws == 10000L) {
    q_matrix <- t(apply(
      panel_weekly_draws, 1L, stats::quantile,
      probs = c(0.025, 0.5, 0.975), names = FALSE
    ))
    current <- data.frame(
      district = as.character(panel$district), week = panel$week,
      current_q025 = q_matrix[, 1L], current_median = q_matrix[, 2L],
      current_q975 = q_matrix[, 3L]
    )
    reference <- data.table::fread(
      input_files[["weekly_rt_reference"]],
      encoding = "UTF-8", data.table = FALSE
    )
    validation <- merge(
      current,
      reference[, c(
        "district", "week", "rt_mc_q025", "rt_mc_median", "rt_mc_q975"
      )],
      by = c("district", "week"), all.x = TRUE, sort = FALSE
    )
    validation$max_row_difference <- apply(cbind(
      abs(validation$current_q025 - validation$rt_mc_q025),
      abs(validation$current_median - validation$rt_mc_median),
      abs(validation$current_q975 - validation$rt_mc_q975)
    ), 1L, max)
    if (
      anyNA(validation$max_row_difference) ||
      max(validation$max_row_difference) > 1e-10
    ) {
      stop("Weekly Rt posterior reconstruction validation failed")
    }
    write_csv(validation, "00_weekly_rt_reconstruction_validation.csv")
  }
  list(draws = panel_weekly_draws, validation = validation)
}

fit_rt_uncertainty_models <- function(analysis_sets, section_label) {
  model_rows <- vector(
    "list", nrow(analysis_sets) * nrow(interventions)
  )
  counter <- 0L
  for (a in seq_len(nrow(analysis_sets))) {
    for (i in seq_len(nrow(interventions))) {
      counter <- counter + 1L
      weather_terms <- if (analysis_sets$Weather_code[a] == "temp_rh") {
        "temp_weekly_avg_lag0 + rh_weekly_avg_lag0"
      } else {
        "temp_weekly_avg_lag0 + total_precipitation_mm"
      }
      model_rows[[counter]] <- data.frame(
        Model_index = counter,
        Seed_model_index = analysis_sets$Seed_model_offset[a] + i,
        Analysis_code = analysis_sets$Analysis_code[a],
        Analysis_scope = group_labels[[analysis_sets$Group_code[a]]],
        Group_code = analysis_sets$Group_code[a],
        Weather_code = analysis_sets$Weather_code[a],
        Weather_adjustment = analysis_sets$Weather_adjustment[a],
        Intervention_index = i,
        Intervention = interventions$Intervention[i],
        Exposure_variable = interventions$Rt_exposure[i],
        Outcome = "log(Rt + 0.1)",
        Rt_offset = RT_OFFSET,
        Week_adjustment = "Centered linear week",
        Formula = paste(
          "log_Rt_plus0p1_work ~", interventions$Rt_exposure[i], "+",
          weather_terms, "+ week_centered + (1 | district)"
        ),
        stringsAsFactors = FALSE
      )
    }
  }
  model_metadata <- do.call(rbind, model_rows)
  model_row_indices <- lapply(seq_len(nrow(model_metadata)), function(j) {
    group_code <- model_metadata$Group_code[j]
    if (group_code == "citywide") {
      seq_len(nrow(panel))
    } else {
      which(panel$Stratum_code == group_code)
    }
  })

  model_templates <- lapply(seq_len(nrow(model_metadata)), function(j) {
    rows <- model_row_indices[[j]]
    d <- droplevels(panel[rows, , drop = FALSE])
    d$log_Rt_plus0p1_work <- log(d$weekly_Rt + RT_OFFSET)
    suppressMessages(suppressWarnings(lme4::lmer(
      stats::as.formula(model_metadata$Formula[j]), data = d,
      REML = TRUE, na.action = na.fail, control = lmer_control
    )))
  })

  point_results <- do.call(rbind, lapply(seq_len(nrow(model_metadata)), function(j) {
    fit <- model_templates[[j]]
    term <- model_metadata$Exposure_variable[j]
    beta <- unname(lme4::fixef(fit)[term])
    se <- sqrt(unname(stats::vcov(fit)[term, term]))
    data.frame(
      model_metadata[j, ], N = stats::nobs(fit),
      Districts = nlevels(fit@frame$district),
      Point_Beta = beta, Point_Std_error = se,
      Point_Rt_plus0p1_reduction_percent = 100 * (1 - exp(beta)),
      Point_singular = lme4::isSingular(fit, tol = 1e-4),
      stringsAsFactors = FALSE
    )
  }))

  message(
    "[", section_label, "] Fitting ",
    N_RT_DRAWS * nrow(model_metadata), " mixed-model refits"
  )
  process_draw_chunk <- function(draw_indices) {
    n_models <- nrow(model_metadata)
    n_rows <- length(draw_indices) * n_models
    core <- data.frame(
      Draw = integer(n_rows), Model_index = integer(n_rows),
      Status = character(n_rows), Beta = numeric(n_rows),
      Std_error = numeric(n_rows), Singular = logical(n_rows),
      Error = character(n_rows), stringsAsFactors = FALSE
    )
    secondary <- matrix(NA_real_, nrow = n_rows, ncol = N_BETA_DRAWS)
    row_counter <- 0L
    for (draw_id in draw_indices) {
      for (model_index in seq_len(n_models)) {
        row_counter <- row_counter + 1L
        core$Draw[row_counter] <- draw_id
        core$Model_index[row_counter] <- model_index
        rows <- model_row_indices[[model_index]]
        new_response <- log(
          weekly_rt$draws[rows, draw_id] + RT_OFFSET
        )
        result <- tryCatch({
          fit <- suppressMessages(suppressWarnings(lme4::refit(
            model_templates[[model_index]], newresp = new_response
          )))
          term <- model_metadata$Exposure_variable[model_index]
          beta <- unname(lme4::fixef(fit)[term])
          se <- sqrt(unname(stats::vcov(fit)[term, term]))
          if (!is.finite(beta) || !is.finite(se) || se < 0) {
            stop("Non-finite intervention estimate")
          }
          seed_index <- model_metadata$Seed_model_index[model_index]
          set.seed(SECONDARY_BETA_SEED + draw_id * 100L + seed_index)
          list(
            beta = beta, se = se,
            singular = lme4::isSingular(fit, tol = 1e-4),
            secondary = beta + se * stats::rnorm(N_BETA_DRAWS),
            error = ""
          )
        }, error = function(e) list(
          beta = NA_real_, se = NA_real_, singular = NA,
          secondary = rep(NA_real_, N_BETA_DRAWS),
          error = conditionMessage(e)
        ))
        core$Status[row_counter] <- if (is.finite(result$beta)) "OK" else "ERROR"
        core$Beta[row_counter] <- result$beta
        core$Std_error[row_counter] <- result$se
        core$Singular[row_counter] <- result$singular
        core$Error[row_counter] <- result$error
        secondary[row_counter, ] <- result$secondary
      }
    }
    list(core = core, secondary = secondary)
  }

  n_workers <- min(
    N_WORKERS_REQUESTED,
    max(1L, parallel::detectCores(logical = FALSE)),
    N_RT_DRAWS
  )
  draw_chunks <- split(seq_len(N_RT_DRAWS), cut(
    seq_len(N_RT_DRAWS),
    breaks = min(N_RT_DRAWS, n_workers * 8L), labels = FALSE
  ))
  if (.Platform$OS.type == "windows" && n_workers > 1L) {
    cluster <- parallel::makeCluster(n_workers, outfile = "")
    tryCatch({
      parallel::clusterEvalQ(
        cluster, suppressPackageStartupMessages(library(lme4))
      )
      parallel::clusterExport(
        cluster,
        c(
          "weekly_rt", "RT_OFFSET", "model_metadata", "model_templates",
          "model_row_indices", "N_BETA_DRAWS", "SECONDARY_BETA_SEED",
          "process_draw_chunk"
        ),
        envir = environment()
      )
      chunk_results <- parallel::parLapply(
        cluster, draw_chunks, process_draw_chunk
      )
    }, finally = parallel::stopCluster(cluster))
  } else {
    chunk_results <- lapply(draw_chunks, process_draw_chunk)
  }

  outer_results <- do.call(rbind, lapply(chunk_results, `[[`, "core"))
  secondary_matrix <- do.call(rbind, lapply(chunk_results, `[[`, "secondary"))
  order_index <- order(outer_results$Draw, outer_results$Model_index)
  outer_results <- outer_results[order_index, , drop = FALSE]
  secondary_matrix <- secondary_matrix[order_index, , drop = FALSE]
  if (nrow(outer_results) != N_RT_DRAWS * nrow(model_metadata)) {
    stop("Unexpected number of fitted Rt-model rows")
  }

  summarise_model <- function(j) {
    rows <- outer_results$Model_index == j
    ok <- rows & outer_results$Status == "OK" &
      is.finite(outer_results$Beta) & is.finite(outer_results$Std_error)
    beta_outer <- outer_results$Beta[ok]
    beta_nested <- as.vector(secondary_matrix[ok, , drop = FALSE])
    beta_nested <- beta_nested[is.finite(beta_nested)]
    if (length(beta_outer) == 0L || length(beta_nested) == 0L) {
      stop("No successful Rt-model draws for model ", j)
    }
    reduction <- 100 * (1 - exp(beta_nested))
    p_negative <- mean(beta_nested < 0)
    p_positive <- mean(beta_nested > 0)
    q_outer <- stats::quantile(
      beta_outer, c(0.025, 0.5, 0.975), names = FALSE
    )
    q_nested <- stats::quantile(
      beta_nested, c(0.025, 0.5, 0.975), names = FALSE
    )
    q_reduction <- stats::quantile(
      reduction, c(0.025, 0.5, 0.975), names = FALSE
    )
    data.frame(
      point_results[j, ],
      Outer_Rt_draws_requested = N_RT_DRAWS,
      Outer_Rt_draws_successful = length(beta_outer),
      Secondary_draws_per_fit = N_BETA_DRAWS,
      Secondary_draws_successful = length(beta_nested),
      Outer_Beta_q025 = q_outer[1], Outer_Beta_median = q_outer[2],
      Outer_Beta_q975 = q_outer[3],
      Nested_Beta_q025 = q_nested[1], Nested_Beta_median = q_nested[2],
      Nested_Beta_q975 = q_nested[3],
      Rt_plus0p1_reduction_percent_q025 = q_reduction[1],
      Rt_plus0p1_reduction_percent_median = q_reduction[2],
      Rt_plus0p1_reduction_percent_q975 = q_reduction[3],
      Probability_reduction = p_negative,
      Directional_p = min(1, 2 * min(p_negative, p_positive)),
      Interval_excludes_zero = q_nested[1] > 0 || q_nested[3] < 0,
      Fit_failure_rate = mean(outer_results$Status[rows] != "OK"),
      Singular_fit_rate = mean(outer_results$Singular[ok], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
  summary_results <- do.call(
    rbind, lapply(seq_len(nrow(model_metadata)), summarise_model)
  )
  if (any(summary_results$Fit_failure_rate > 0)) {
    warning("One or more Rt models had failed posterior-draw fits")
  }
  list(
    summary = summary_results,
    formulas = model_metadata,
    outer = outer_results,
    secondary = secondary_matrix
  )
}

fit_bi_models <- function(analysis_sets) {
  result_rows <- list()
  formula_rows <- list()
  model_objects <- list()
  counter <- 0L
  for (a in seq_len(nrow(analysis_sets))) {
    d <- subset_panel(analysis_sets$Group_code[a])
    weather_terms <- if (analysis_sets$Weather_code[a] == "temp_rh") {
      "temp_weekly_avg_lag0 + rh_weekly_avg_lag0"
    } else {
      "temp_weekly_avg_lag0 + total_precipitation_mm"
    }
    for (i in seq_len(nrow(interventions))) {
      counter <- counter + 1L
      exposure <- interventions$BI_exposure[i]
      formula_text <- paste(
        "log_BI_current ~", exposure, "+", weather_terms,
        "+ week_centered + (1 | district)"
      )
      fit <- suppressMessages(suppressWarnings(lmerTest::lmer(
        stats::as.formula(formula_text), data = d, REML = TRUE,
        na.action = na.fail, control = lmer_control
      )))
      coefficient_table <- summary(fit)$coefficients
      beta <- unname(coefficient_table[exposure, "Estimate"])
      se <- unname(coefficient_table[exposure, "Std. Error"])
      df <- unname(coefficient_table[exposure, "df"])
      p_value <- unname(coefficient_table[exposure, "Pr(>|t|)"])
      critical <- stats::qt(0.975, df = df)
      ci_lower <- beta - critical * se
      ci_upper <- beta + critical * se
      result_rows[[counter]] <- data.frame(
        Analysis_code = analysis_sets$Analysis_code[a],
        Analysis_scope = group_labels[[analysis_sets$Group_code[a]]],
        Group_code = analysis_sets$Group_code[a],
        Weather_code = analysis_sets$Weather_code[a],
        Weather_adjustment = analysis_sets$Weather_adjustment[a],
        Outcome = "log(BI + 1)",
        Intervention_index = i,
        Intervention = interventions$Intervention[i],
        Exposure_variable = exposure,
        Week_adjustment = "Centered linear week",
        Formula = formula_text,
        N = stats::nobs(fit), Districts = nlevels(fit@frame$district),
        Beta = beta, Standard_error = se, Degrees_freedom = df,
        CI_lower = ci_lower, CI_upper = ci_upper, P_value = p_value,
        Reduction_percent = 100 * (1 - exp(beta)),
        Reduction_CI_lower = 100 * (1 - exp(ci_upper)),
        Reduction_CI_upper = 100 * (1 - exp(ci_lower)),
        Singular_fit = lme4::isSingular(fit, tol = 1e-4),
        stringsAsFactors = FALSE
      )
      formula_rows[[counter]] <- data.frame(
        Analysis_scope = group_labels[[analysis_sets$Group_code[a]]],
        Weather_adjustment = analysis_sets$Weather_adjustment[a],
        Intervention = interventions$Intervention[i],
        Formula = formula_text,
        stringsAsFactors = FALSE
      )
      model_objects[[counter]] <- fit
    }
  }
  list(
    summary = do.call(rbind, result_rows),
    formulas = do.call(rbind, formula_rows),
    models = model_objects
  )
}

fit_random_intercept_mediation <- function() {
  groups <- data.frame(
    Group_index = 1:3,
    Group_code = c("citywide", "central", "other"),
    Group = unname(group_labels[c("citywide", "central", "other")]),
    stringsAsFactors = FALSE
  )
  draw_rows <- list()
  diagnostic_rows <- list()
  formula_rows <- list()
  model_objects <- list()
  counter <- 0L

  message(
    "[2. Mediation] Fitting nine random-intercept model pairs with ",
    N_MEDIATION_SIMS, " quasi-Bayesian simulations each"
  )
  for (g in seq_len(nrow(groups))) {
    d <- subset_panel(groups$Group_code[g])
    for (i in seq_len(nrow(interventions))) {
      counter <- counter + 1L
      exposure <- interventions$Rt_exposure[i]
      mediator_formula <- stats::as.formula(paste(
        "log_BI_previous ~", exposure,
        "+ temp_weekly_avg_lag0 + rh_weekly_avg_lag0 +",
        "week_centered + (1 | district)"
      ))
      outcome_formula <- stats::as.formula(paste(
        "log_Rt_plus0p1 ~", exposure,
        "+ log_BI_previous + temp_weekly_avg_lag0 + rh_weekly_avg_lag0 +",
        "week_centered + (1 | district)"
      ))
      mediator_model <- suppressMessages(suppressWarnings(lme4::lmer(
        mediator_formula, data = d, REML = FALSE,
        na.action = na.fail, control = lmer_control
      )))
      outcome_model <- suppressMessages(suppressWarnings(lme4::lmer(
        outcome_formula, data = d, REML = FALSE,
        na.action = na.fail, control = lmer_control
      )))

      simulation_seed <- MEDIATION_SEED_BASE + 1000L + g * 100L + i
      set.seed(simulation_seed)
      mediation_fit <- suppressMessages(suppressWarnings(mediation::mediate(
        model.m = mediator_model,
        model.y = outcome_model,
        treat = exposure,
        mediator = "log_BI_previous",
        sims = N_MEDIATION_SIMS,
        boot = FALSE,
        group.out = "district"
      )))
      acme <- (
        as.numeric(mediation_fit$d0.sims) + as.numeric(mediation_fit$d1.sims)
      ) / 2
      ade <- (
        as.numeric(mediation_fit$z0.sims) + as.numeric(mediation_fit$z1.sims)
      ) / 2
      total <- as.numeric(mediation_fit$tau.sims)
      proportion <- (
        as.numeric(mediation_fit$n0.sims) + as.numeric(mediation_fit$n1.sims)
      ) / 2
      if (!all(lengths(list(acme, ade, total, proportion)) == N_MEDIATION_SIMS)) {
        stop("Unexpected mediation simulation count")
      }

      path_a <- unname(lme4::fixef(mediator_model)[exposure])
      path_a_se <- sqrt(unname(stats::vcov(mediator_model)[exposure, exposure]))
      path_b <- unname(lme4::fixef(outcome_model)["log_BI_previous"])
      path_b_se <- sqrt(unname(
        stats::vcov(outcome_model)["log_BI_previous", "log_BI_previous"]
      ))

      draw_rows[[counter]] <- data.table::data.table(
        Simulation = seq_len(N_MEDIATION_SIMS), Model_index = counter,
        Specification_code = "random_intercept",
        Specification = "District random intercept",
        Group_index = g, Group_code = groups$Group_code[g],
        Group = groups$Group[g],
        Intervention_index = i, Intervention = interventions$Intervention[i],
        Exposure_variable = exposure,
        ACME_log_Rt_plus0p1 = acme,
        ADE_log_Rt_plus0p1 = ade,
        Total_log_Rt_plus0p1 = total,
        Proportion_mediated = proportion
      )
      diagnostic_rows[[counter]] <- data.frame(
        Model_index = counter,
        Specification_code = "random_intercept",
        Specification = "District random intercept",
        Group_index = g, Group_code = groups$Group_code[g],
        Group = groups$Group[g],
        Intervention_index = i, Intervention = interventions$Intervention[i],
        Exposure_variable = exposure,
        N = stats::nobs(outcome_model), Districts = nlevels(d$district),
        Simulation_seed = simulation_seed,
        Path_a = path_a, Path_a_SE = path_a_se,
        Path_a_normal_P = 2 * stats::pnorm(-abs(path_a / path_a_se)),
        Path_b = path_b, Path_b_SE = path_b_se,
        Path_b_normal_P = 2 * stats::pnorm(-abs(path_b / path_b_se)),
        Mediator_singular = lme4::isSingular(mediator_model, tol = 1e-4),
        Outcome_singular = lme4::isSingular(outcome_model, tol = 1e-4),
        stringsAsFactors = FALSE
      )
      formula_rows[[counter]] <- data.frame(
        Model_index = counter,
        Specification = "District random intercept",
        Group = groups$Group[g],
        Intervention = interventions$Intervention[i],
        Mediator_formula = paste(deparse(mediator_formula), collapse = " "),
        Outcome_formula = paste(deparse(outcome_formula), collapse = " "),
        stringsAsFactors = FALSE
      )
      model_objects[[counter]] <- list(
        mediator_model = mediator_model, outcome_model = outcome_model
      )
    }
  }

  mediation_draws <- data.table::rbindlist(draw_rows)
  diagnostics <- do.call(rbind, diagnostic_rows)
  formulas <- do.call(rbind, formula_rows)
  summarise_reduction <- function(x) {
    stats::quantile(100 * (1 - exp(x)), c(0.025, 0.5, 0.975))
  }
  summary_results <- mediation_draws[, {
    acme_log <- stats::quantile(
      ACME_log_Rt_plus0p1, c(0.025, 0.5, 0.975)
    )
    acme_reduction <- summarise_reduction(ACME_log_Rt_plus0p1)
    ade_reduction <- summarise_reduction(ADE_log_Rt_plus0p1)
    total_reduction <- summarise_reduction(Total_log_Rt_plus0p1)
    p_indirect <- mean(ACME_log_Rt_plus0p1 < 0)
    p_direct <- mean(ADE_log_Rt_plus0p1 < 0)
    p_total <- mean(Total_log_Rt_plus0p1 < 0)
    list(
      Rt_offset = RT_OFFSET,
      Outcome = "log(Rt + 0.1)",
      Rt_uncertainty_propagated = FALSE,
      Coefficient_simulations = .N,
      ACME_log_q025 = acme_log[1],
      ACME_log_median = acme_log[2],
      ACME_log_q975 = acme_log[3],
      Indirect_Rt_plus0p1_reduction_q025 = acme_reduction[1],
      Indirect_Rt_plus0p1_reduction_median = acme_reduction[2],
      Indirect_Rt_plus0p1_reduction_q975 = acme_reduction[3],
      Probability_indirect_reduction = p_indirect,
      Indirect_directional_p = min(1, 2 * min(p_indirect, 1 - p_indirect)),
      Direct_Rt_plus0p1_reduction_q025 = ade_reduction[1],
      Direct_Rt_plus0p1_reduction_median = ade_reduction[2],
      Direct_Rt_plus0p1_reduction_q975 = ade_reduction[3],
      Probability_direct_reduction = p_direct,
      Total_Rt_plus0p1_reduction_q025 = total_reduction[1],
      Total_Rt_plus0p1_reduction_median = total_reduction[2],
      Total_Rt_plus0p1_reduction_q975 = total_reduction[3],
      Probability_total_reduction = p_total,
      Proportion_mediated_q025 = stats::quantile(
        Proportion_mediated, 0.025
      ),
      Proportion_mediated_median = stats::median(Proportion_mediated),
      Proportion_mediated_q975 = stats::quantile(
        Proportion_mediated, 0.975
      )
    )
  }, by = .(
    Model_index, Specification_code, Specification,
    Group_index, Group_code, Group,
    Intervention_index, Intervention, Exposure_variable
  )]
  summary_results <- merge(
    summary_results, diagnostics,
    by = c(
      "Model_index", "Specification_code", "Specification",
      "Group_index", "Group_code", "Group",
      "Intervention_index", "Intervention", "Exposure_variable"
    ),
    all.x = TRUE, sort = FALSE
  )
  data.table::setorder(summary_results, Group_index, Intervention_index)
  list(
    summary = summary_results,
    formulas = formulas,
    diagnostics = diagnostics,
    draws = mediation_draws,
    models = model_objects
  )
}

weekly_rt <- reconstruct_weekly_rt(N_RT_DRAWS)

# -----------------------------------------------------------------------------
# 1. Primary analysis: temperature + relative humidity (both lag 0)
# -----------------------------------------------------------------------------

message("[1. Primary] Intervention associations adjusted for temperature and RH")
main_sets <- data.frame(
  Analysis_code = "city_temp_rh",
  Group_code = "citywide",
  Weather_code = "temp_rh",
  Weather_adjustment = "Temperature + relative humidity (lag 0)",
  Seed_model_offset = 0L,
  stringsAsFactors = FALSE
)
main_rt <- fit_rt_uncertainty_models(main_sets, "1. Primary Rt")
main_bi <- fit_bi_models(main_sets)
write_csv(main_rt$summary, "01_main_intervention_logRt_plus0p1_temp_rh.csv")
write_csv(main_bi$summary, "01_main_intervention_logBI1_temp_rh.csv")

# -----------------------------------------------------------------------------
# 2. Mediation analysis: temperature + relative humidity (both lag 0)
# -----------------------------------------------------------------------------

# Temporal sequence:
#   intervention(d,t-2) -> log[BI(d,t-1)+1] -> log[Rt(d,t)+0.1]
# The weekly Rt posterior median is used as the outcome. mediate() propagates
# coefficient uncertainty (sims=10000, boot=FALSE by default) but does not
# propagate Rt posterior uncertainty through the mediation analysis.
mediation_results <- fit_random_intercept_mediation()
write_csv(
  mediation_results$summary,
  "02_mediation_temp_rh_random_intercept.csv"
)
write_csv(
  mediation_results$diagnostics,
  "02_mediation_temp_rh_random_intercept_diagnostics.csv"
)

# -----------------------------------------------------------------------------
# 3. Stratified analysis: central urban versus other districts
# -----------------------------------------------------------------------------

message("[3. Stratified] Central urban versus other districts")
stratified_sets <- data.frame(
  Analysis_code = c("central_temp_rh", "other_temp_rh"),
  Group_code = c("central", "other"),
  Weather_code = c("temp_rh", "temp_rh"),
  Weather_adjustment = rep(
    "Temperature + relative humidity (lag 0)", 2
  ),
  # Offsets preserve the model-specific secondary-draw seeds from the original
  # combined 12-model analysis: central models 7-9; other models 10-12.
  Seed_model_offset = c(6L, 9L),
  stringsAsFactors = FALSE
)
stratified_rt <- fit_rt_uncertainty_models(
  stratified_sets, "3. Stratified Rt"
)
stratified_bi <- fit_bi_models(stratified_sets)
write_csv(
  stratified_rt$summary,
  "03_stratified_intervention_logRt_plus0p1_temp_rh.csv"
)
write_csv(
  stratified_bi$summary,
  "03_stratified_intervention_logBI1_temp_rh.csv"
)

# -----------------------------------------------------------------------------
# 4. Sensitivity analysis: temperature + total rainfall (both lag 0)
# -----------------------------------------------------------------------------

message("[4. Sensitivity] Substituting total rainfall for relative humidity")
sensitivity_sets <- data.frame(
  Analysis_code = "city_temp_rain",
  Group_code = "citywide",
  Weather_code = "temp_rain",
  Weather_adjustment = "Temperature + total rainfall (lag 0)",
  # Offset preserves original combined-model secondary-draw seeds 4-6.
  Seed_model_offset = 3L,
  stringsAsFactors = FALSE
)
sensitivity_rt <- fit_rt_uncertainty_models(
  sensitivity_sets, "4. Sensitivity Rt"
)
sensitivity_bi <- fit_bi_models(sensitivity_sets)
write_csv(
  sensitivity_rt$summary,
  "04_sensitivity_intervention_logRt_plus0p1_temp_rain.csv"
)
write_csv(
  sensitivity_bi$summary,
  "04_sensitivity_intervention_logBI1_temp_rain.csv"
)

# -----------------------------------------------------------------------------
# 5. Formula audit, validation, combined workbook, and software record
# -----------------------------------------------------------------------------

formula_audit <- rbind(
  data.frame(
    Analysis_order = 1L, Analysis_family = "Primary intervention-log(Rt+0.1)",
    Scope = main_rt$formulas$Analysis_scope,
    Weather = main_rt$formulas$Weather_adjustment,
    Intervention = main_rt$formulas$Intervention,
    Mediator_formula = NA_character_,
    Outcome_formula = main_rt$formulas$Formula,
    stringsAsFactors = FALSE
  ),
  data.frame(
    Analysis_order = 1L, Analysis_family = "Primary intervention-log(BI+1)",
    Scope = main_bi$formulas$Analysis_scope,
    Weather = main_bi$formulas$Weather_adjustment,
    Intervention = main_bi$formulas$Intervention,
    Mediator_formula = NA_character_,
    Outcome_formula = main_bi$formulas$Formula,
    stringsAsFactors = FALSE
  ),
  data.frame(
    Analysis_order = 2L, Analysis_family = "Random-intercept mediation",
    Scope = mediation_results$formulas$Group,
    Weather = "Temperature + relative humidity (lag 0)",
    Intervention = mediation_results$formulas$Intervention,
    Mediator_formula = mediation_results$formulas$Mediator_formula,
    Outcome_formula = mediation_results$formulas$Outcome_formula,
    stringsAsFactors = FALSE
  ),
  data.frame(
    Analysis_order = 3L, Analysis_family = "Stratified intervention-log(Rt+0.1)",
    Scope = stratified_rt$formulas$Analysis_scope,
    Weather = stratified_rt$formulas$Weather_adjustment,
    Intervention = stratified_rt$formulas$Intervention,
    Mediator_formula = NA_character_,
    Outcome_formula = stratified_rt$formulas$Formula,
    stringsAsFactors = FALSE
  ),
  data.frame(
    Analysis_order = 3L, Analysis_family = "Stratified intervention-log(BI+1)",
    Scope = stratified_bi$formulas$Analysis_scope,
    Weather = stratified_bi$formulas$Weather_adjustment,
    Intervention = stratified_bi$formulas$Intervention,
    Mediator_formula = NA_character_,
    Outcome_formula = stratified_bi$formulas$Formula,
    stringsAsFactors = FALSE
  ),
  data.frame(
    Analysis_order = 4L, Analysis_family = "Sensitivity intervention-log(Rt+0.1)",
    Scope = sensitivity_rt$formulas$Analysis_scope,
    Weather = sensitivity_rt$formulas$Weather_adjustment,
    Intervention = sensitivity_rt$formulas$Intervention,
    Mediator_formula = NA_character_,
    Outcome_formula = sensitivity_rt$formulas$Formula,
    stringsAsFactors = FALSE
  ),
  data.frame(
    Analysis_order = 4L, Analysis_family = "Sensitivity intervention-log(BI+1)",
    Scope = sensitivity_bi$formulas$Analysis_scope,
    Weather = sensitivity_bi$formulas$Weather_adjustment,
    Intervention = sensitivity_bi$formulas$Intervention,
    Mediator_formula = NA_character_,
    Outcome_formula = sensitivity_bi$formulas$Formula,
    stringsAsFactors = FALSE
  )
)
write_csv(formula_audit, "05_formula_audit.csv")

validation <- data.frame(
  Check = c(
    "Primary Rt analysis has three intervention results",
    "Primary BI analysis has three intervention results",
    "Mediation has citywide and stratified random-intercept results only",
    "Stratified Rt analysis has six intervention results",
    "Stratified BI analysis has six intervention results",
    "Sensitivity Rt analysis has three intervention results",
    "Sensitivity BI analysis has three intervention results",
    "Every Rt model completed all requested posterior draws",
    "Every Rt model completed all requested secondary coefficient draws",
    "All mediation models used the requested quasi-Bayesian simulations",
    "No district fixed-effect mediation model is present",
    "Mediation did not propagate Rt posterior uncertainty",
    "All models adjust centered linear week",
    "Rainfall appears only in sensitivity models",
    "All current-week weather variables use lag 0"
  ),
  Passed = c(
    nrow(main_rt$summary) == 3L,
    nrow(main_bi$summary) == 3L,
    nrow(mediation_results$summary) == 9L &&
      all(mediation_results$summary$Specification_code == "random_intercept") &&
      all(c("citywide", "central", "other") %in%
        mediation_results$summary$Group_code),
    nrow(stratified_rt$summary) == 6L,
    nrow(stratified_bi$summary) == 6L,
    nrow(sensitivity_rt$summary) == 3L,
    nrow(sensitivity_bi$summary) == 3L,
    all(c(
      main_rt$summary$Outer_Rt_draws_successful,
      stratified_rt$summary$Outer_Rt_draws_successful,
      sensitivity_rt$summary$Outer_Rt_draws_successful
    ) == N_RT_DRAWS),
    all(c(
      main_rt$summary$Secondary_draws_successful,
      stratified_rt$summary$Secondary_draws_successful,
      sensitivity_rt$summary$Secondary_draws_successful
    ) == N_RT_DRAWS * N_BETA_DRAWS),
    all(mediation_results$summary$Coefficient_simulations == N_MEDIATION_SIMS),
    !any(grepl(
      "fixed", c(
        mediation_results$summary$Specification,
        mediation_results$formulas$Mediator_formula,
        mediation_results$formulas$Outcome_formula
      ), ignore.case = TRUE
    )),
    all(!mediation_results$summary$Rt_uncertainty_propagated),
    all(grepl("week_centered", formula_audit$Outcome_formula, fixed = TRUE)) &&
      all(grepl(
        "week_centered",
        formula_audit$Mediator_formula[!is.na(formula_audit$Mediator_formula)],
        fixed = TRUE
      )),
    all(grepl(
      "total_precipitation_mm",
      formula_audit$Outcome_formula[formula_audit$Analysis_order == 4L],
      fixed = TRUE
    )) && !any(grepl(
      "total_precipitation_mm",
      formula_audit$Outcome_formula[formula_audit$Analysis_order != 4L],
      fixed = TRUE
    )),
    all(grepl(
      "temp_weekly_avg_lag0", formula_audit$Outcome_formula, fixed = TRUE
    )) && all(grepl(
      "rh_weekly_avg_lag0",
      formula_audit$Outcome_formula[formula_audit$Analysis_order %in% 1:3],
      fixed = TRUE
    ))
  ),
  stringsAsFactors = FALSE
)
write_csv(validation, "05_analysis_validation.csv")
if (!all(validation$Passed)) {
  stop(
    "Analysis validation failed: ",
    paste(validation$Check[!validation$Passed], collapse = "; ")
  )
}

writexl::write_xlsx(
  list(
    `01_main_Rt_temp_RH` = main_rt$summary,
    `01_main_BI_temp_RH` = main_bi$summary,
    `02_mediation_random` = as.data.frame(mediation_results$summary),
    `02_mediation_diagnostics` = mediation_results$diagnostics,
    `03_stratified_Rt` = stratified_rt$summary,
    `03_stratified_BI` = stratified_bi$summary,
    `04_sensitivity_Rt_rain` = sensitivity_rt$summary,
    `04_sensitivity_BI_rain` = sensitivity_bi$summary,
    `05_formula_audit` = formula_audit,
    `05_validation` = validation
  ),
  file.path(OUTPUT_DIR, "submission_analysis_results.xlsx")
)

if (SAVE_INTERMEDIATE_DRAWS) {
  saveRDS(
    list(
      main = main_rt[c("outer", "secondary")],
      stratified = stratified_rt[c("outer", "secondary")],
      sensitivity = sensitivity_rt[c("outer", "secondary")],
      mediation_draws = mediation_results$draws
    ),
    file.path(OUTPUT_DIR, "intermediate_simulation_draws.rds"),
    compress = "xz"
  )
}

metadata <- c(
  paste0("Completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("R version: ", R.version.string),
  paste0("Rt offset: ", RT_OFFSET),
  paste0("BI offset: ", BI_OFFSET),
  paste0("Rt posterior trajectories: ", N_RT_DRAWS),
  paste0("Secondary coefficient draws per Rt fit: ", N_BETA_DRAWS),
  paste0("Mediation quasi-Bayesian simulations: ", N_MEDIATION_SIMS),
  paste0("Rt posterior seed: ", OUTER_RT_SEED),
  paste0("Secondary coefficient seed base: ", SECONDARY_BETA_SEED),
  paste0("Mediation seed base: ", MEDIATION_SEED_BASE),
  "Primary covariates: current-week temperature and relative humidity",
  "Sensitivity covariates: current-week temperature and total rainfall",
  "Week adjustment: centered linear epidemiological week",
  "District term: random intercept in every model",
  "Mediation district fixed effects: excluded",
  "Mediation Rt posterior uncertainty propagated: FALSE",
  "Percent transformation for Rt models applies to Rt+0.1, not directly to Rt"
)
writeLines(
  metadata, file.path(OUTPUT_DIR, "analysis_run_metadata.txt"), useBytes = TRUE
)
writeLines(
  capture.output(sessionInfo()),
  file.path(OUTPUT_DIR, "sessionInfo.txt"), useBytes = TRUE
)

message("All analyses completed and validated successfully")
print(validation, row.names = FALSE)
