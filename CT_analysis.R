# Multivariable linear regression of Ct values and interaction analyses
#
# best_skill: epidemiology, cross-sectional multivariable linear regression,
#             publication-ready coefficient and interaction tables
# train_signal: estimate adjusted coefficients and P values for interactions
# selection_split: 1,080-record primary Ct analysis dataset
# heldout_gate: continuous-age and categorical-day sensitivity specifications
# accepted_patterns: complete-case audit, explicit references, HC3 inference
# rejected_patterns: silent row deletion, unlabelled factor references
# patch_scope: one outcome, one additive model, two interaction models
# reject_if: unexpected levels, missing model variables, or N != 1,080
#
# Primary specification:
#   Ct = beta0 + beta1(days from onset to testing)
#        + beta2(age <18) + beta3(age >=60)
#        + beta4(female) + beta5(imported case) + error
#
# Reference groups:
#   age 18-59 years, male, local case.
#
# Inference:
#   OLS point estimates with HC3 heteroskedasticity-robust standard errors.
#   The age-by-days P for interaction is a joint 2-df robust Wald F test.
#   The sex-by-days P for interaction is a 1-df robust Wald F test.

if (.Platform$OS.type == "windows") {
  suppressWarnings(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.utf8"))
}

required_packages <- c("dplyr", "openxlsx", "sandwich", "lmtest", "car")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Missing R packages: ", paste(missing_packages, collapse = ", "),
    ". Install them before running this script."
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(openxlsx)
  library(sandwich)
  library(lmtest)
  library(car)
})

all_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", all_args, value = TRUE)
script_dir <- if (length(file_arg) > 0) {
  script_path_arg <- sub("^--file=", "", file_arg[1])
  script_parent <- dirname(script_path_arg)
  if (script_parent %in% c(".", "")) "." else script_parent
} else {
  "."
}

user_args <- commandArgs(trailingOnly = TRUE)
input_path <- if (length(user_args) >= 1) {
  user_args[1]
} else {
  file.path(script_dir, "ct_new_outputs_v2", "ct_analysis_data_primary.xlsx")
}
output_dir <- if (length(user_args) >= 2) {
  user_args[2]
} else {
  file.path(script_dir, "ct_new_outputs_v2", "multivariable_regression")
}

if (!file.exists(input_path)) stop("Input file not found: ", input_path)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Stage files in an ASCII-only temporary directory because some Windows R
# locale configurations cannot unzip or create workbooks under Chinese paths.
runtime_dir <- file.path(tempdir(), "ct_multivariable_runtime")
dir.create(runtime_dir, recursive = TRUE, showWarnings = FALSE)
runtime_input <- file.path(runtime_dir, "ct_analysis_data_primary.xlsx")
if (!file.copy(input_path, runtime_input, overwrite = TRUE)) {
  stop("Failed to stage the input workbook.")
}

copy_from_runtime <- function(runtime_file, final_name) {
  final_file <- file.path(output_dir, final_name)
  if (!file.copy(runtime_file, final_file, overwrite = TRUE)) {
    stop("Failed to copy output file to: ", final_file)
  }
  invisible(final_file)
}

write_xlsx_safely <- function(x, final_name) {
  runtime_file <- file.path(runtime_dir, final_name)
  openxlsx::write.xlsx(x, runtime_file, overwrite = TRUE)
  copy_from_runtime(runtime_file, final_name)
}

write_lines_safely <- function(x, final_name) {
  runtime_file <- file.path(runtime_dir, final_name)
  writeLines(x, runtime_file, useBytes = TRUE)
  copy_from_runtime(runtime_file, final_name)
}

raw_data <- openxlsx::read.xlsx(
  runtime_input,
  sheet = 1,
  detectDates = FALSE,
  check.names = FALSE
)

required_columns <- c("CT", "days", "age", "age_group", "gender", "case_type")
if (!all(required_columns %in% names(raw_data))) {
  stop(
    "Primary analysis workbook is missing columns: ",
    paste(setdiff(required_columns, names(raw_data)), collapse = ", ")
  )
}
if (nrow(raw_data) != 1080L) {
  stop("Expected 1,080 primary-analysis records, but found ", nrow(raw_data))
}

model_data <- raw_data %>%
  transmute(
    CT = as.numeric(CT),
    days = as.numeric(days),
    age_years = as.numeric(age),
    age_group_raw = as.character(age_group),
    gender_raw = as.character(gender),
    case_type_raw = as.character(case_type)
  ) %>%
  mutate(
    age_group = factor(
      case_when(
        age_group_raw == "0-17" ~ "Age_lt18",
        age_group_raw == "18-59" ~ "Age_18_59",
        age_group_raw == "60~" ~ "Age_ge60",
        TRUE ~ NA_character_
      ),
      levels = c("Age_18_59", "Age_lt18", "Age_ge60")
    ),
    sex = factor(
      case_when(
        gender_raw == "male" ~ "Male",
        gender_raw == "female" ~ "Female",
        TRUE ~ NA_character_
      ),
      levels = c("Male", "Female")
    ),
    case_class = factor(
      case_when(
        case_type_raw == "\u672c\u5730" ~ "Local",
        case_type_raw %in% c("\u8f93\u5165", "\u5883\u5916\u8f93\u5165") ~ "Imported",
        TRUE ~ NA_character_
      ),
      levels = c("Local", "Imported")
    ),
    age_per_10y = age_years / 10,
    days_centered = days - mean(days),
    day_group = factor(
      if_else(days >= 6, "Day_ge6", paste0("Day_", as.integer(days))),
      levels = c(paste0("Day_", 0:5), "Day_ge6")
    )
  )

model_variables <- c(
  "CT", "days", "age_years", "age_group", "sex", "case_class", "day_group"
)
missing_by_variable <- vapply(
  model_data[model_variables],
  function(x) sum(is.na(x)),
  integer(1)
)

sample_audit <- tibble(
  item = c(
    "Records in primary Ct dataset",
    "Records with complete additive-model variables",
    "Records excluded from additive model",
    "Minimum days from onset to testing",
    "Maximum days from onset to testing",
    "Local cases",
    "Imported cases",
    "Male",
    "Female",
    "Age <18 years",
    "Age 18-59 years",
    "Age >=60 years"
  ),
  value = c(
    nrow(model_data),
    sum(complete.cases(model_data[c("CT", "days", "age_group", "sex", "case_class")])),
    sum(!complete.cases(model_data[c("CT", "days", "age_group", "sex", "case_class")])),
    min(model_data$days, na.rm = TRUE),
    max(model_data$days, na.rm = TRUE),
    sum(model_data$case_class == "Local", na.rm = TRUE),
    sum(model_data$case_class == "Imported", na.rm = TRUE),
    sum(model_data$sex == "Male", na.rm = TRUE),
    sum(model_data$sex == "Female", na.rm = TRUE),
    sum(model_data$age_group == "Age_lt18", na.rm = TRUE),
    sum(model_data$age_group == "Age_18_59", na.rm = TRUE),
    sum(model_data$age_group == "Age_ge60", na.rm = TRUE)
  )
)

if (any(missing_by_variable > 0)) {
  stop(
    "Missing values remain in model variables: ",
    paste(
      names(missing_by_variable)[missing_by_variable > 0],
      missing_by_variable[missing_by_variable > 0],
      sep = "=",
      collapse = ", "
    )
  )
}
if (nrow(model_data) != 1080L) stop("Model-data row count changed unexpectedly.")
if (!identical(sort(unique(model_data$case_type_raw)),
               sort(c("\u5883\u5916\u8f93\u5165", "\u672c\u5730", "\u8f93\u5165")))) {
  stop("Unexpected source case-type levels; review case-classification mapping.")
}
if (!all(model_data$days >= 0 & model_data$days <= 11)) {
  stop("Days from onset to testing fall outside the expected 0-11 range.")
}

robust_tidy <- function(model, model_name, term_labels = NULL) {
  V <- sandwich::vcovHC(model, type = "HC3")
  test <- lmtest::coeftest(model, vcov. = V)
  result <- tibble(
    model = model_name,
    term = rownames(test),
    coefficient = as.numeric(test[, 1]),
    robust_SE_HC3 = as.numeric(test[, 2]),
    statistic = as.numeric(test[, 3]),
    p_value = as.numeric(test[, 4]),
    CI95_lower = coefficient - qnorm(0.975) * robust_SE_HC3,
    CI95_upper = coefficient + qnorm(0.975) * robust_SE_HC3
  )
  if (!is.null(term_labels)) {
    result <- result %>%
      mutate(
        term_label = unname(term_labels[term]),
        term_label = if_else(is.na(term_label), term, term_label)
      ) %>%
      select(model, term, term_label, everything())
  }
  result
}

wald_test <- function(model, hypotheses, label) {
  V <- sandwich::vcovHC(model, type = "HC3")
  test <- car::linearHypothesis(
    model,
    hypotheses,
    vcov. = V,
    test = "F",
    singular.ok = FALSE
  )
  last_row <- nrow(test)
  tibble(
    interaction = label,
    numerator_df = as.numeric(test$Df[last_row]),
    denominator_df = as.numeric(test$Res.Df[last_row]),
    F_statistic = as.numeric(test$F[last_row]),
    p_for_interaction = as.numeric(test$`Pr(>F)`[last_row]),
    covariance = "HC3 heteroskedasticity-robust",
    status = "Estimated"
  )
}

linear_combination <- function(model, weights, label) {
  V <- sandwich::vcovHC(model, type = "HC3")
  b <- coef(model)
  L <- rep(0, length(b))
  names(L) <- names(b)
  unknown <- setdiff(names(weights), names(b))
  if (length(unknown) > 0) {
    stop("Unknown coefficient(s) in linear combination: ", paste(unknown, collapse = ", "))
  }
  L[names(weights)] <- weights
  estimate <- sum(L * b)
  se_value <- sqrt(as.numeric(t(L) %*% V %*% L))
  tibble(
    subgroup = label,
    coefficient_per_1_day = estimate,
    robust_SE_HC3 = se_value,
    CI95_lower = estimate - qnorm(0.975) * se_value,
    CI95_upper = estimate + qnorm(0.975) * se_value,
    p_value = 2 * pnorm(abs(estimate / se_value), lower.tail = FALSE)
  )
}

main_model <- lm(
  CT ~ days + age_group + sex + case_class,
  data = model_data,
  na.action = na.fail
)

main_labels <- c(
  "(Intercept)" = "Intercept: day 0, age 18-59, male, local case",
  "days" = "Days from symptom onset to testing (per 1 day)",
  "age_groupAge_lt18" = "Age <18 vs 18-59 years",
  "age_groupAge_ge60" = "Age >=60 vs 18-59 years",
  "sexFemale" = "Female vs male",
  "case_classImported" = "Imported vs local case"
)
main_results <- robust_tidy(main_model, "Primary additive model", main_labels)

# Sensitivity model treating age as a linear continuous variable.
continuous_age_model <- lm(
  CT ~ days + age_per_10y + sex + case_class,
  data = model_data,
  na.action = na.fail
)
continuous_age_labels <- c(
  "(Intercept)" = "Intercept: day 0, age 0, male, local case",
  "days" = "Days from symptom onset to testing (per 1 day)",
  "age_per_10y" = "Age (per 10-year increase)",
  "sexFemale" = "Female vs male",
  "case_classImported" = "Imported vs local case"
)
continuous_age_results <- robust_tidy(
  continuous_age_model,
  "Continuous-age sensitivity model",
  continuous_age_labels
)

# Age group by days interaction; the P for interaction is a joint 2-df test.
age_interaction_model <- lm(
  CT ~ days * age_group + sex + case_class,
  data = model_data,
  na.action = na.fail
)
age_interaction_results <- robust_tidy(
  age_interaction_model,
  "Age-group by days interaction model"
)
age_interaction_test <- wald_test(
  age_interaction_model,
  c(
    "days:age_groupAge_lt18 = 0",
    "days:age_groupAge_ge60 = 0"
  ),
  "Age group x days (overall 2-df test)"
)

age_specific_slopes <- bind_rows(
  linear_combination(
    age_interaction_model,
    c("days" = 1),
    "Age 18-59 years"
  ),
  linear_combination(
    age_interaction_model,
    c("days" = 1, "days:age_groupAge_lt18" = 1),
    "Age <18 years"
  ),
  linear_combination(
    age_interaction_model,
    c("days" = 1, "days:age_groupAge_ge60" = 1),
    "Age >=60 years"
  )
)

# Continuous age by days interaction as a sensitivity analysis.
continuous_age_interaction_model <- lm(
  CT ~ days * age_per_10y + sex + case_class,
  data = model_data,
  na.action = na.fail
)
continuous_age_interaction_test <- wald_test(
  continuous_age_interaction_model,
  "days:age_per_10y = 0",
  "Continuous age (per 10 years) x days"
)

# Sex by days interaction; one interaction coefficient, 1-df robust Wald test.
sex_interaction_model <- lm(
  CT ~ days * sex + age_group + case_class,
  data = model_data,
  na.action = na.fail
)
sex_interaction_results <- robust_tidy(
  sex_interaction_model,
  "Sex by days interaction model"
)
sex_interaction_test <- wald_test(
  sex_interaction_model,
  "days:sexFemale = 0",
  "Sex x days (overall 1-df test)"
)

sex_specific_slopes <- bind_rows(
  linear_combination(
    sex_interaction_model,
    c("days" = 1),
    "Male"
  ),
  linear_combination(
    sex_interaction_model,
    c("days" = 1, "days:sexFemale" = 1),
    "Female"
  )
)

# Sensitivity specification using the publication day categories
# (Day 0-5 and >=6 days), adjusted for the same covariates.
categorical_day_model <- lm(
  CT ~ day_group + age_group + sex + case_class,
  data = model_data,
  na.action = na.fail
)
categorical_day_results <- robust_tidy(
  categorical_day_model,
  "Categorical-day sensitivity model"
)
categorical_day_test <- wald_test(
  categorical_day_model,
  paste0("day_group", c(paste0("Day_", 1:5), "Day_ge6"), " = 0"),
  "Day category (overall 6-df test; reference Day 0)"
)

# Sensitivity interactions allowing non-linear day-category patterns.
age_by_day_group_model <- lm(
  CT ~ day_group * age_group + sex + case_class,
  data = model_data,
  na.action = na.fail
)
age_day_cell_counts <- model_data %>%
  count(day_group, age_group, name = "n", .drop = FALSE)
age_by_day_group_terms <- as.vector(outer(
  paste0("day_group", c(paste0("Day_", 1:5), "Day_ge6")),
  c("age_groupAge_lt18", "age_groupAge_ge60"),
  paste,
  sep = ":"
))
age_by_day_group_test <- if (min(age_day_cell_counts$n) < 5) {
  tibble(
    interaction = "Age group x 7-category day group (sensitivity; overall 12-df test)",
    numerator_df = 12,
    denominator_df = df.residual(age_by_day_group_model),
    F_statistic = NA_real_,
    p_for_interaction = NA_real_,
    covariance = "Not calculated",
    status = paste0(
      "Not reported because sparse age-by-day cells made HC3 inference unstable; ",
      "minimum cell n=", min(age_day_cell_counts$n)
    )
  )
} else {
  wald_test(
    age_by_day_group_model,
    paste0(age_by_day_group_terms, " = 0"),
    "Age group x 7-category day group (sensitivity; overall 12-df test)"
  )
}

sex_by_day_group_model <- lm(
  CT ~ day_group * sex + age_group + case_class,
  data = model_data,
  na.action = na.fail
)
sex_by_day_group_terms <- paste0(
  "day_group",
  c(paste0("Day_", 1:5), "Day_ge6"),
  ":sexFemale"
)
sex_by_day_group_test <- wald_test(
  sex_by_day_group_model,
  paste0(sex_by_day_group_terms, " = 0"),
  "Sex x 7-category day group (sensitivity; overall 6-df test)"
)

# Flexible cubic time trend: a joint test of quadratic and cubic components
# assesses departure from the primary linear day trend.
cubic_day_model <- lm(
  CT ~ days_centered + I(days_centered^2) + I(days_centered^3) +
    age_group + sex + case_class,
  data = model_data,
  na.action = na.fail
)
cubic_nonlinearity_test <- wald_test(
  cubic_day_model,
  c("I(days_centered^2) = 0", "I(days_centered^3) = 0"),
  "Departure from linear day trend (quadratic and cubic terms; 2-df test)"
)

# Flexible interaction sensitivity models use a cubic day curve without
# saturating sparse age-by-day cells.
age_flexible_interaction_model <- lm(
  CT ~ (days_centered + I(days_centered^2) + I(days_centered^3)) *
    age_group + sex + case_class,
  data = model_data,
  na.action = na.fail
)
age_flexible_interaction_terms <- c(
  "days_centered:age_groupAge_lt18",
  "days_centered:age_groupAge_ge60",
  "I(days_centered^2):age_groupAge_lt18",
  "I(days_centered^2):age_groupAge_ge60",
  "I(days_centered^3):age_groupAge_lt18",
  "I(days_centered^3):age_groupAge_ge60"
)
age_flexible_interaction_test <- wald_test(
  age_flexible_interaction_model,
  paste0(age_flexible_interaction_terms, " = 0"),
  "Age group x flexible cubic day trend (sensitivity; overall 6-df test)"
)

sex_flexible_interaction_model <- lm(
  CT ~ (days_centered + I(days_centered^2) + I(days_centered^3)) *
    sex + age_group + case_class,
  data = model_data,
  na.action = na.fail
)
sex_flexible_interaction_terms <- c(
  "days_centered:sexFemale",
  "I(days_centered^2):sexFemale",
  "I(days_centered^3):sexFemale"
)
sex_flexible_interaction_test <- wald_test(
  sex_flexible_interaction_model,
  paste0(sex_flexible_interaction_terms, " = 0"),
  "Sex x flexible cubic day trend (sensitivity; overall 3-df test)"
)

interaction_tests <- bind_rows(
  age_interaction_test,
  sex_interaction_test,
  continuous_age_interaction_test,
  age_by_day_group_test,
  sex_by_day_group_test,
  age_flexible_interaction_test,
  sex_flexible_interaction_test,
  cubic_nonlinearity_test
)

# Model diagnostics and influence audit.
V_main <- sandwich::vcovHC(main_model, type = "HC3")
bp <- lmtest::bptest(main_model)
influence_cutoff <- 4 / nobs(main_model)
cooks <- cooks.distance(main_model)
vif_values <- car::vif(main_model)
if (is.matrix(vif_values)) {
  vif_table <- tibble(
    term = rownames(vif_values),
    GVIF = vif_values[, "GVIF"],
    df = vif_values[, "Df"],
    adjusted_GVIF = vif_values[, "GVIF"]^(1 / (2 * vif_values[, "Df"]))
  )
} else {
  vif_table <- tibble(
    term = names(vif_values),
    GVIF = as.numeric(vif_values),
    df = 1,
    adjusted_GVIF = sqrt(as.numeric(vif_values))
  )
}

diagnostics <- tibble(
  diagnostic = c(
    "N",
    "R-squared",
    "Adjusted R-squared",
    "Residual degrees of freedom",
    "Breusch-Pagan statistic",
    "Breusch-Pagan P value",
    "Maximum Cook's distance",
    "Cook's distance cutoff (4/N)",
    "Records above Cook's distance cutoff"
  ),
  value = c(
    nobs(main_model),
    summary(main_model)$r.squared,
    summary(main_model)$adj.r.squared,
    df.residual(main_model),
    unname(bp$statistic),
    bp$p.value,
    max(cooks),
    influence_cutoff,
    sum(cooks > influence_cutoff)
  )
)

# Outcome-based influence deletion is a sensitivity analysis only, not a new
# primary sample restriction.
influence_sensitivity_model <- lm(
  CT ~ days + age_group + sex + case_class,
  data = model_data[cooks <= influence_cutoff, ],
  na.action = na.fail
)
influence_sensitivity_results <- robust_tidy(
  influence_sensitivity_model,
  "Sensitivity excluding Cook's distance >4/N",
  main_labels
)

model_specification <- tibble(
  item = c(
    "Analysis unit",
    "Outcome",
    "Primary exposure/time variable",
    "Age coding",
    "Sex reference",
    "Case-classification reference",
    "Missing-data handling",
    "Estimator",
    "Inference",
    "Age P for interaction",
    "Sex P for interaction"
  ),
  definition = c(
    "One patient record",
    "Arithmetic-mean Ct value reconstructed from the raw Ct field",
    "Continuous days from symptom onset (report date fallback when needed) to testing",
    "Three categories; 18-59 years is the reference",
    "Male",
    "Local case; imported combines imported and overseas-imported source categories",
    "Complete-case analysis on variables in each model",
    "Ordinary least squares linear regression",
    "HC3 heteroskedasticity-robust standard errors, two-sided tests, 95% CIs",
    "Joint robust Wald F test of both age-group-by-days product terms",
    "Robust Wald F test of the sex-by-days product term"
  )
)

write_xlsx_safely(
  list(
    Model_specification = model_specification,
    Sample_audit = sample_audit,
    Main_model_HC3 = main_results,
    Continuous_age_model = continuous_age_results,
    Interaction_tests = interaction_tests,
    Age_interaction_coefficients = age_interaction_results,
    Age_specific_day_slopes = age_specific_slopes,
    Sex_interaction_coefficients = sex_interaction_results,
    Sex_specific_day_slopes = sex_specific_slopes,
    Categorical_day_model = categorical_day_results,
    Categorical_day_overall_test = categorical_day_test,
    Age_by_day_cell_counts = age_day_cell_counts,
    Categorical_sex_interaction = robust_tidy(
      sex_by_day_group_model,
      "Sex by categorical day interaction"
    ),
    Flexible_age_interaction = robust_tidy(
      age_flexible_interaction_model,
      "Age group by flexible cubic day interaction"
    ),
    Flexible_sex_interaction = robust_tidy(
      sex_flexible_interaction_model,
      "Sex by flexible cubic day interaction"
    ),
    Cubic_day_model = robust_tidy(
      cubic_day_model,
      "Cubic day-trend sensitivity model"
    ),
    Influence_sensitivity = influence_sensitivity_results,
    Diagnostics = diagnostics,
    VIF = vif_table
  ),
  "ct_multivariable_regression_results.xlsx"
)

write.csv(
  main_results,
  file.path(runtime_dir, "ct_main_model_coefficients.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
copy_from_runtime(
  file.path(runtime_dir, "ct_main_model_coefficients.csv"),
  "ct_main_model_coefficients.csv"
)

write.csv(
  interaction_tests,
  file.path(runtime_dir, "ct_interaction_tests.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
copy_from_runtime(
  file.path(runtime_dir, "ct_interaction_tests.csv"),
  "ct_interaction_tests.csv"
)

plot_slopes_data <- bind_rows(
  age_specific_slopes %>% mutate(stratification = "Age group"),
  sex_specific_slopes %>% mutate(stratification = "Sex")
) %>%
  select(stratification, everything())
write.csv(
  plot_slopes_data,
  file.path(runtime_dir, "ct_interaction_slopes_plot_data.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
copy_from_runtime(
  file.path(runtime_dir, "ct_interaction_slopes_plot_data.csv"),
  "ct_interaction_slopes_plot_data.csv"
)

write.csv(
  categorical_day_results,
  file.path(runtime_dir, "ct_categorical_day_model_coefficients.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
copy_from_runtime(
  file.path(runtime_dir, "ct_categorical_day_model_coefficients.csv"),
  "ct_categorical_day_model_coefficients.csv"
)

format_p <- function(x) {
  ifelse(x < 0.001, "<0.001", sprintf("%.3f", x))
}

main_display <- main_results %>%
  filter(term != "(Intercept)") %>%
  transmute(
    Variable = term_label,
    Coefficient = sprintf("%.3f", coefficient),
    `Robust SE` = sprintf("%.3f", robust_SE_HC3),
    `95% CI` = sprintf("%.3f to %.3f", CI95_lower, CI95_upper),
    `P value` = format_p(p_value)
  )

interaction_display <- interaction_tests %>%
  filter(
    interaction %in% c(
      "Age group x days (overall 2-df test)",
      "Sex x days (overall 1-df test)",
      "Age group x flexible cubic day trend (sensitivity; overall 6-df test)",
      "Sex x flexible cubic day trend (sensitivity; overall 3-df test)",
      "Age group x 7-category day group (sensitivity; overall 12-df test)",
      "Sex x 7-category day group (sensitivity; overall 6-df test)"
    )
  ) %>%
  transmute(
    Interaction = interaction,
    `F statistic` = sprintf("%.3f", F_statistic),
    `P for interaction` = format_p(p_for_interaction),
    Status = status
  )

md_table <- function(data) {
  display <- as.data.frame(data, stringsAsFactors = FALSE)
  display[] <- lapply(display, function(column) {
    result <- ifelse(is.na(column), "", as.character(column))
    gsub("\\|", "\\\\|", result)
  })
  header <- paste0("| ", paste(names(display), collapse = " | "), " |")
  separator <- paste0(
    "| ", paste(rep("---", ncol(display)), collapse = " | "), " |"
  )
  rows <- apply(
    display,
    1,
    function(row) paste0("| ", paste(row, collapse = " | "), " |")
  )
  paste(c(header, separator, rows), collapse = "\n")
}

day_beta <- main_results %>% filter(term == "days")
age_int_p <- age_interaction_test$p_for_interaction
sex_int_p <- sex_interaction_test$p_for_interaction
categorical_day_p <- categorical_day_test$p_for_interaction
nonlinear_p <- cubic_nonlinearity_test$p_for_interaction

results_sentence <- paste0(
  "In the multivariable linear regression model including 1,080 patients, ",
  "each additional day from symptom onset to testing was associated with a ",
  sprintf("%.2f", day_beta$coefficient),
  "-cycle change in the mean Ct value (95% CI ",
  sprintf("%.2f", day_beta$CI95_lower), " to ",
  sprintf("%.2f", day_beta$CI95_upper), "; P",
  ifelse(day_beta$p_value < 0.001, "<0.001", paste0("=", sprintf("%.3f", day_beta$p_value))),
  "), after adjustment for age group, sex, and case classification. ",
  "There was ",
  ifelse(age_int_p < 0.05, "evidence", "no evidence"),
  " of effect modification by age group (P for interaction ",
  ifelse(age_int_p < 0.001, "<0.001", paste0("=", sprintf("%.3f", age_int_p))),
  ") and ",
  ifelse(sex_int_p < 0.05, "evidence", "no evidence"),
  " of effect modification by sex (P for interaction ",
  ifelse(sex_int_p < 0.001, "<0.001", paste0("=", sprintf("%.3f", sex_int_p))),
  ")."
)

report_lines <- c(
  "# Multivariable regression of Ct values",
  "",
  "## Analysis specification",
  "",
  "Ct value was modelled as a continuous outcome using ordinary least-squares linear regression. The primary time variable was the continuous number of days from symptom onset to testing. Age was entered using the prespecified categories of <18, 18-59, and >=60 years; 18-59 years was the reference. Male sex and local cases were the reference categories. Imported cases combined the imported and overseas-imported categories in the source data. All 1,080 primary-analysis records had complete model variables. HC3 heteroskedasticity-robust standard errors, two-sided P values, and 95% confidence intervals were reported.",
  "",
  "Separate models added age-group-by-days and sex-by-days product terms. The age P for interaction was obtained from a joint 2-df robust Wald F test of both age interaction terms. The sex P for interaction was obtained from a 1-df robust Wald F test. A continuous-age model, a seven-category day model (Day 0-5 and >=6 days), categorical-day interaction models, and a cubic time-trend model were fitted as sensitivity analyses. Deletion of observations with Cook's distance >4/N was examined only as an influence sensitivity analysis and did not redefine the primary sample.",
  "",
  "## Main adjusted model",
  "",
  md_table(main_display),
  "",
  "## Interaction tests",
  "",
  md_table(interaction_display),
  "",
  "## Suggested Results text",
  "",
  results_sentence,
  "",
  "## Diagnostics",
  "",
  paste0(
    "The primary model had R-squared=",
    sprintf("%.3f", summary(main_model)$r.squared),
    " and adjusted R-squared=",
    sprintf("%.3f", summary(main_model)$adj.r.squared),
    ". The Breusch-Pagan test P value was ",
    format_p(bp$p.value),
    "; therefore, HC3 robust inference was retained. The overall adjusted ",
    "categorical-day association had P",
    ifelse(categorical_day_p < 0.001, "<0.001", paste0("=", sprintf("%.3f", categorical_day_p))),
    ". Departure from a strictly linear time trend had P",
    ifelse(nonlinear_p < 0.001, "<0.001", paste0("=", sprintf("%.3f", nonlinear_p))),
    " in the cubic sensitivity model."
  ),
  "",
  "## Interpretation note",
  "",
  "Coefficients are adjusted mean differences in Ct cycles. The coefficient for days is the adjusted mean change in Ct for a one-day increase in the onset-to-testing interval. Interaction tests assess whether this per-day Ct slope differs across age groups or between sexes; they do not test whether the groups have different Ct values at a single day.",
  ""
)

write_lines_safely(
  report_lines,
  "CT_multivariable_regression_report_English.md"
)

message("Ct multivariable regression completed. Output directory: ", output_dir)
print(main_results)
print(interaction_tests)
print(age_specific_slopes)
print(sex_specific_slopes)
