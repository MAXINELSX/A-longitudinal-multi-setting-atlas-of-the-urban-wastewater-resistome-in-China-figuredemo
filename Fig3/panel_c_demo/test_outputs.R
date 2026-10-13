script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
panel_demo_root <- dirname(normalizePath(script, winslash = "/", mustWork = TRUE))
source(file.path(panel_demo_root, "tests", "test_outputs.R"))
