
argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
root <- file.path(dirname(dirname(normalizePath(sub("^--file=", "", argument)))), "e_technical_demo")
source(file.path(root, "run_demo.R"))
source(file.path(root, "source", "analysis_functions.R"))
input <- load_technical_inputs(root)
meta <- align_metadata(input$matrix, input$metadata)
stopifnot(nrow(input$matrix) == 72L, ncol(input$matrix) == 2646L,
          all(table(meta$Sample_Type) == 18L),
          setequal(meta$Sample_Type, c("Hospital", "WWTP", "Community", "Wet Market")),
          identical(as.character(meta$Sample), rownames(input$matrix)))

x <- input$matrix[1, ]; y <- input$matrix[2, ]
manual_bray <- sum(abs(x - y)) / sum(x + y)
vegan_bray <- as.matrix(vegan::vegdist(input$matrix[1:2, , drop = FALSE], method = "bray"))[1, 2]
stopifnot(abs(manual_bray - vegan_bray) < 1e-12)
reordered <- align_metadata(input$matrix, input$metadata[72:1, ])
stopifnot(identical(reordered$Sample, meta$Sample))
bad <- input$metadata[-1, ]
err <- tryCatch({ align_metadata(input$matrix, bad); NULL }, error = identity)
stopifnot(inherits(err, "error"), grepl("Sample IDs differ", conditionMessage(err), fixed = TRUE))
cat("PASS: Fig1 technical\n")
