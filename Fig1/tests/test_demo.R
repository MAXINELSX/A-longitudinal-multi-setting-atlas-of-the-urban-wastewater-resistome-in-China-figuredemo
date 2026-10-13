
args <- commandArgs(trailingOnly = FALSE)
test_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])
root <- dirname(dirname(normalizePath(test_file, mustWork = TRUE)))
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2); library(cowplot); library(scales); library(grid)
})
source(file.path(root, "source", "plot_functions.R"))
meta <- load_metadata(file.path(root, "inputs", "metadata.csv"))
stopifnot(nrow(meta) == 96L, length(unique(meta$Sample)) == 96L,
          length(unique(meta$City_Code)) == 3L, all(table(meta$Sample_Type) == 24L),
          all(table(meta$Sample_Type, meta$City_Code) == 8L))
expected <- read.csv(file.path(root, "expected", "summary.csv"))
stopifnot(abs(sum(meta$Clean_Reads_GBp) - expected$value[expected$metric == "total_clean_gbp"]) < 1e-10)


raw <- read.csv(file.path(root, "inputs", "metadata.csv"), stringsAsFactors = FALSE)
expect_error <- function(data, message_part) {
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path), add = TRUE)
  write.csv(data, path, row.names = FALSE, na = "")
  err <- tryCatch({ load_metadata(path); NULL }, error = identity)
  stopifnot(inherits(err, "error"), grepl(message_part, conditionMessage(err), fixed = TRUE))
}
duplicate <- raw; duplicate$Sample[2] <- duplicate$Sample[1]
expect_error(duplicate, "unique")
negative <- raw; negative$Clean_Reads_GBp[1] <- -1
expect_error(negative, "negative")
bad_date <- raw; bad_date$Sample_Date[1] <- "not-a-date"
expect_error(bad_date, "dates")
bad_city <- raw; bad_city$City[1] <- "Unknown"
expect_error(bad_city, "unknown city")
bad_setting <- raw; bad_setting$Sample_Type[1] <- "Unknown"
expect_error(bad_setting, "unknown setting")
expect_error(raw[, setdiff(names(raw), "Clean_Reads_GBp")], "Missing metadata columns")
stopifnot(identical(normalize_setting(c("wetmarket", "Wet market", " WWTP ")),
                    c("Wet Market", "Wet Market", "WWTP")))
curves <- list(setting = read.csv(file.path(root, "inputs", "precomputed_setting_curves.csv")),
               city = read.csv(file.path(root, "inputs", "precomputed_city_curves.csv")))
stopifnot(nrow(curves$setting) == 132L, nrow(curves$city) == 718L,
          setequal(curves$setting$Group, site_order), setequal(curves$city$Group, names(city_palette)))
for (curve in c(split(curves$setting, curves$setting$Group), split(curves$city, curves$city$Group))) {
  stopifnot(all(diff(curve$N) > 0), all(curve$Subtypes >= 0),
            all(is.finite(curve$Subtypes)))
}
cat("PASS: Fig1 metadata\n")
