# ==============================================================================
# Reviewer entry point: complete household-transmission analysis
#
# This script deliberately excludes all plotting. It runs:
#   1) analysis-data construction from gz_keep_first_distinct.csv
#   2) descriptive household-transmission analysis
#   3) complete-case adjusted regression
#   4) age-by-household-size interaction analysis (numeric results only)
#
# Figures are generated separately by:
#   reviewer_code/02_make_household_transmission_figures.R
# ==============================================================================

options(stringsAsFactors = FALSE, encoding = "UTF-8")
Sys.unsetenv(c("LC_ALL", "LC_CTYPE"))

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source_csv <- file.path(project_dir, "gz_keep_first_distinct.csv")

if (!file.exists(source_csv)) {
  stop("The required source file was not found: ", source_csv)
}

# The current revision uses only gz_keep_first_distinct.csv.
# Ct is derived directly from the CSV column named 'ct':
#   missing/non-numeric/<=0 -> missing; valid <30 or >=30.
Sys.setenv(
  HOUSEHOLD_ANALYSIS_PROFILE = "direct_csv",
  HOUSEHOLD_DATA_PATH = source_csv,
  HOUSEHOLD_INITIAL_DATA_PATH = source_csv,
  HOUSEHOLD_CT_DEFINITION = "direct"
)

analysis_modules <- c(
  file.path("analysis", "01_data_audit.R"),
  file.path("analysis", "02_household_descriptive.R"),
  file.path("analysis", "10_complete_case_regression.R"),
  file.path("analysis", "14_reviewer_age_household_interaction.R")
)

for (module in analysis_modules) {
  if (!file.exists(module)) {
    stop("Required analysis module was not found: ", module)
  }
  message("Running reviewer module: ", module)
  source(module, encoding = "UTF-8", chdir = FALSE)
}

# ------------------------------------------------------------------------------
# Consolidated reviewer workbook and reproducibility checks
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(openxlsx)
})

output_dir <- file.path(project_dir, "outputs_direct_csv")
interaction_dir <- file.path(
  output_dir,
  "age_household_interaction_review"
)

read_result <- function(directory, filename) {
  path <- file.path(directory, filename)
  if (!file.exists(path)) {
    stop("Expected output was not found: ", path)
  }
  utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    fileEncoding = "UTF-8"
  )
}

reviewer_tables <- list(
  "Data audit" = read_result(output_dir, "data_audit_checks.csv"),
  "Missingness" = read_result(output_dir, "missingness_summary.csv"),
  "Descriptive results" = read_result(
    output_dir,
    "household_descriptive_by_group.csv"
  ),
  "Regression flow" = read_result(
    output_dir,
    "complete_case_sample_flow.csv"
  ),
  "Regression missing" = read_result(
    output_dir,
    "complete_case_missingness.csv"
  ),
  "Adjusted RR" = read_result(
    output_dir,
    "complete_case_regression_rr.csv"
  ),
  "Regression diagnostics" = read_result(
    output_dir,
    "complete_case_model_diagnostics.csv"
  ),
  "Logbin selection gate" = read_result(
    output_dir,
    "logbinomial_selection_gate.csv"
  ),
  "Interaction cells" = read_result(
    interaction_dir,
    "01_joint_cell_counts.csv"
  ),
  "Interaction diagnostics" = read_result(
    interaction_dir,
    "02_model_diagnostics.csv"
  ),
  "Multiplicative interaction" = read_result(
    interaction_dir,
    "04_multiplicative_interaction.csv"
  ),
  "Global interaction tests" = read_result(
    interaction_dir,
    "05_global_interaction_tests.csv"
  ),
  "Standardized risks" = read_result(
    interaction_dir,
    "06_standardized_risks.csv"
  ),
  "Additive interaction" = read_result(
    interaction_dir,
    "07_additive_interaction.csv"
  ),
  "Interaction support" = read_result(
    interaction_dir,
    "09_analysis_support.csv"
  )
)

descriptive <- reviewer_tables[["Descriptive results"]]
regression_flow <- reviewer_tables[["Regression flow"]]
regression_diagnostics <- reviewer_tables[["Regression diagnostics"]]
regression_gate <- reviewer_tables[["Logbin selection gate"]]
global_tests <- reviewer_tables[["Global interaction tests"]]

overall <- descriptive[descriptive$variable_key == "overall", , drop = FALSE][1, ]
complete_row <- regression_flow[
  regression_flow$step == "Complete-case regression sample",
  ,
  drop = FALSE
]

qa_checks <- tibble::tibble(
  check = c(
    "Source data contain 932 households",
    "Descriptive analysis contains 2523 co-residents",
    "Descriptive analysis contains 28 infected co-residents",
    "Descriptive analysis contains 25 event households",
    "Complete-case regression contains 811 households",
    "Complete-case regression contains 22 event households",
    "Main regression selected log-binomial",
    "All main-regression log-binomial gates passed",
    "Interaction multiplicative global test is finite",
    "Interaction additive global test is finite",
    "No interaction figure is required by this analysis entry point"
  ),
  passed = c(
    overall$households == 932L,
    overall$co_residents == 2523L,
    overall$infected_co_residents == 28L,
    overall$event_households == 25L,
    complete_row$households == 811L,
    complete_row$event_households == 22L,
    regression_diagnostics$model[1] == "Log-binomial regression",
    all(regression_gate$passed),
    is.finite(global_tests$p_value[
      grepl("multiplicative.*joint Wald", global_tests$test)
    ][1]),
    is.finite(global_tests$p_value[
      grepl("Global additive", global_tests$test)
    ][1]),
    TRUE
  )
)

qa_path <- file.path(output_dir, "reviewer_household_analysis_QA.csv")
utils::write.csv(
  qa_checks,
  qa_path,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

if (!all(qa_checks$passed)) {
  stop(
    "Reviewer analysis QA failed: ",
    paste(qa_checks$check[!qa_checks$passed], collapse = "; ")
  )
}

workbook <- openxlsx::createWorkbook()
for (sheet_name in names(reviewer_tables)) {
  # Excel worksheet names cannot exceed 31 characters.
  safe_sheet_name <- substr(sheet_name, 1L, 31L)
  openxlsx::addWorksheet(workbook, safe_sheet_name)
  openxlsx::writeData(
    workbook,
    safe_sheet_name,
    reviewer_tables[[sheet_name]]
  )
  openxlsx::freezePane(workbook, safe_sheet_name, firstRow = TRUE)
  openxlsx::setColWidths(
    workbook,
    safe_sheet_name,
    cols = seq_len(ncol(reviewer_tables[[sheet_name]])),
    widths = "auto"
  )
}
openxlsx::addWorksheet(workbook, "QA")
openxlsx::writeData(workbook, "QA", qa_checks)
openxlsx::setColWidths(workbook, "QA", cols = 1:2, widths = "auto")

combined_workbook <- file.path(
  output_dir,
  "reviewer_household_transmission_analysis.xlsx"
)
openxlsx::saveWorkbook(workbook, combined_workbook, overwrite = TRUE)

formula_lines <- c(
  "Reviewer household-transmission analysis formulas",
  "",
  "Main adjusted model:",
  "log{P(Y=1)} = beta0 + age + sex + work environment + household size + Ct category + detection period + region.",
  "",
  "Multiplicative interaction model:",
  "log{P(Y=1)} = beta0 + age + household size + age*household size + covariates.",
  "",
  "Additive interaction:",
  "IC_a = R(a,large) - R(a,small) - R(<18,large) + R(<18,small).",
  "RERI = RR11 - RR10 - RR01 + 1.",
  "",
  "Log-binomial regression is preferred. Modified Poisson regression with HC0 robust variance is used when the prespecified log-binomial fitting gates fail."
)
writeLines(
  formula_lines,
  file.path(output_dir, "reviewer_household_analysis_formulas.txt"),
  useBytes = TRUE
)

session_lines <- c(
  paste0("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("Project directory: ", project_dir),
  paste0("Source data: ", source_csv),
  paste0("Consolidated workbook: ", combined_workbook),
  "",
  capture.output(utils::sessionInfo())
)
writeLines(
  session_lines,
  file.path(output_dir, "reviewer_household_analysis_session_info.txt"),
  useBytes = TRUE
)

message("Complete reviewer household-transmission analysis finished.")
message("Consolidated workbook: ", combined_workbook)
message("QA file: ", qa_path)
message("Interaction analysis was numeric only; no interaction figures were made.")

