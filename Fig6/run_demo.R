


args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
demo_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/")) else getwd()
output_dir <- file.path(demo_dir, "demo_output")
if (length(args)) {
  if (length(args) != 2L || args[1] != "--output") stop("Usage: Rscript run_demo.R [--output DIR]")
  output_dir <- args[2]
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
required <- c("randomForest", "rfPermute", "lme4", "lmerTest", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required R packages first: ", paste(missing, collapse = ", "))
source(file.path(demo_dir, "source", "panel_helpers.R"), encoding = "UTF-8")
started <- proc.time()[["elapsed"]]
d <- read.delim(file.path(demo_dir, "input", "samples.tsv"), check.names = FALSE, stringsAsFactors = FALSE, na.strings = "NA")
stopifnot(nrow(d) >= 40L, !anyDuplicated(d$sample_id), length(unique(d$city)) >= 3L,
          length(unique(d$site_nm)) >= 6L, all(grepl("^City[0-9]+$", d$city)),
          all(grepl("^Hospital[0-9]+$", d$site_nm)), all(grepl("^Sample[0-9]+$", d$sample_id)))
config <- fig6_network_config()
config$ntree <- 200L
config$num.rep <- 99L
config$num.cores <- 1L


config$labels["prec_sum"] <- "Precipitation proxy"
predictors <- unique(unlist(config$classes, use.names = FALSE))
for (j in seq_along(predictors)) {
  lag30_field <- paste0(predictors[j], "_lag1")
  values <- d[[lag30_field]]
  spread <- stats::sd(values, na.rm = TRUE)
  if (!is.finite(spread) || spread == 0) spread <- 1
  for (k in c(0L, 2L, 3L)) {
    perturbation <- spread * (0.08 + 0.02 * k) *
      sin(seq_len(nrow(d)) * (0.37 + j / 17) + k * 1.13 + j)
    d[[paste0(predictors[j], "_lag", k)]] <- values + perturbation
  }
}
prepared <- fig6_network_prepare_lmm(d, config)





lag_columns <- paste0(predictors, "_lag1")
rf_data <- data.frame(city = factor(prepared$city),
                      setting = factor(prepared$setting),
                      time = factor(format(as.Date(prepared$sample_date), "%Y-%m")),
                      arg_abundance = prepared$arg_abundance,
                      check.names = FALSE)
rf_data[predictors] <- prepared[lag_columns]
rf_complete <- complete.cases(rf_data)
rf_data <- rf_data[rf_complete, , drop = FALSE]
stopifnot(!anyNA(rf_data), nrow(rf_data) >= 40L)
constant <- names(rf_data)[vapply(rf_data, function(x) length(unique(x)) < 2L, logical(1))]
stopifnot(!config$outcome %in% constant)
rf_data <- rf_data[setdiff(names(rf_data), constant)]
rf_data <- droplevels(rf_data)





set.seed(config$seed)
rf_model <- fig6_rf_permute(stats::reformulate(setdiff(names(rf_data), config$outcome), response = config$outcome),
                           data = rf_data, importance = TRUE, ntree = config$ntree,
                           num.rep = config$num.rep, num.cores = config$num.cores)
importance <- as.data.frame(randomForest::importance(rf_model, scale = TRUE), check.names = FALSE)
importance$variable <- rownames(importance)
rf_table <- fig6_network_rf_table(importance, "Total", config)
lmm <- fig6_network_fit_lmm(prepared, "Total", config)
stopifnot(nrow(rf_table) == length(predictors) + 2L,
          nrow(lmm$results) == length(predictors) * length(config$lags),
          nrow(lmm$diagnostics) == length(config$classes) * length(config$lags),
          all(lmm$diagnostics$Status == "ok"),
          all(is.finite(lmm$results$coef)), all(is.finite(lmm$results$p)), all(is.finite(lmm$results$fdr)))



run_fig6(output_dir, network_rf = list(Total = rf_table), network_lmm = list(Total = lmm),
         panels = c("a", "e"), config = config, combine = FALSE, family = "sans")

session_file <- file.path(output_dir, "sessionInfo.txt")
session_text <- readLines(session_file, warn = FALSE)
session_text <- session_text[!grepl("^time zone:", session_text)]
writeLines(session_text, session_file)
write.csv(importance, file.path(output_dir, "rf_importance.csv"), row.names = FALSE, na = "")
write.csv(lmm$results, file.path(output_dir, "lmm_coefficients.csv"), row.names = FALSE, na = "")
write.csv(lmm$diagnostics, file.path(output_dir, "lmm_diagnostics.csv"), row.names = FALSE, na = "")
write.csv(lmm$preprocessing, file.path(output_dir, "preprocessing.csv"), row.names = FALSE, na = "")
write.csv(data.frame(Panel = letters[1:8],
                     Status = ifelse(letters[1:8] %in% c("a", "e"), "ACTUALLY_RUN_ADAPTED_DEMO", "NOT_TESTED")),
          file.path(output_dir, "panel_status.csv"), row.names = FALSE)
metadata <- data.frame(
  Key = c("sample_rows", "cities", "hospitals", "rf_predictors", "rf_trees", "rf_permutations",
          "rf_seed", "rf_input_window_days", "lmm_window_origin", "rf_constant_predictors_removed", "rf_imputation", "rf_complete_rows", "rf_incomplete_rows_excluded",
          "lmm_models_successful", "lmm_coefficients", "lmm_singular_models", "lmm_warning_models",
          "elapsed_seconds", "rf_model_class", "randomForest_version", "rfPermute_version", "lmerTest_version"),
  Value = c(nrow(d), length(unique(d$city)), length(unique(d$site_nm)), nrow(importance),
            config$ntree, config$num.rep, config$seed, 30, "lag1 selected; lag0/2/3 deterministic synthetic",
            paste(constant, collapse = ";"), "complete_case", sum(rf_complete), sum(!rf_complete),
            sum(lmm$diagnostics$Status == "ok"), nrow(lmm$results), sum(lmm$diagnostics$Singular),
            sum(nzchar(lmm$diagnostics$Warnings)), round(proc.time()[["elapsed"]] - started, 3),
            paste(class(rf_model), collapse = ";"), as.character(packageVersion("randomForest")),
            as.character(packageVersion("rfPermute")), as.character(packageVersion("lmerTest"))))
write.csv(metadata, file.path(output_dir, "run_metadata.csv"), row.names = FALSE)
source(file.path(demo_dir, "tests", "verify_outputs.R"), local = TRUE)
verify_outputs(output_dir)
message("PASS: Fig6 a/e; NOT_TESTED: b/c/d/f/g/h")
