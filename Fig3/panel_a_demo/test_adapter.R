args <- commandArgs(trailingOnly = FALSE)
script <- sub('^--file=', '', args[grepl('^--file=', args)][1])
root <- dirname(normalizePath(script, winslash = '/', mustWork = TRUE))
adapter <- file.path(root, 'adapter.R')
if (!file.exists(adapter)) stop('FAIL: the relative-abundance adapter is not implemented.')
source(adapter)
x <- read.delim(file.path(root, 'input', 'relative_abundance.tsv'), check.names = FALSE)
m <- read.delim(file.path(root, 'input', 'metadata.tsv'), check.names = FALSE)
d <- prepare_panel_a(x, m)
stopifnot(nrow(x) == 6L, ncol(x) == 17L, nrow(m) == 16L, nrow(d) == 48L)
stopifnot(all(grepl('^Sample[0-9]{3}$', m$sample)), all(grepl('^City[0-9]{2}$', m$city)))
stopifnot(all(d$n_samples == 2L), !'measurement' %in% names(d))
for (i in seq_len(nrow(d))) {
  ids <- m$sample[m$city == d$city[i] & m$setting == d$setting[i]]
  expected <- mean(as.numeric(x[x$species == d$species[i], ids, drop = TRUE]))
  stopifnot(abs(expected - d$mean_relative_abundance[i]) < 1e-14)
}
expect_error <- function(expr) stopifnot(inherits(tryCatch({force(expr); NULL}, error = identity), 'error'))
expect_error(prepare_panel_a(x, rbind(m, m[1, ])))
bad <- x; bad[1, 2] <- -1
expect_error(prepare_panel_a(bad, m))
expect_error(prepare_panel_a(x, m[-1, ]))
cat('PASS_FIG3A_ADAPTER\n')
