# =============================================================================
# check_reproduction.R -- compares a fresh run with the output used in the paper
#
# Author : Bruno Duarte Gomes (ORCID 0000-0002-3646-7847)
# Contact: brunodgomes@ufpa.br
#
# Usage  : Rscript check_reproduction.R <new_dir> <stored_dir>
#          e.g. copy results/ to paper_results/, re-run run_simulation_study.R
#          (which overwrites results/), then
#          Rscript check_reproduction.R results paper_results
# Reports, for every series x estimator, whether the classification and the
# fitted exponent agree with the stored output.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
new_dir <- if (length(args) >= 1) args[1] else "results"
old_dir <- if (length(args) >= 2) args[2] else "paper_results"
a <- read.csv(file.path(old_dir, "rows.csv"))
b <- read.csv(file.path(new_dir, "rows.csv"))
stopifnot(identical(a$label, b$label), identical(a$method, b$method))
same <- a$cls == b$cls & (is.na(a$alpha) & is.na(b$alpha) |
                          abs(a$alpha - b$alpha) < 1e-8 & !is.na(a$alpha) & !is.na(b$alpha))
tab <- aggregate(same, list(series = sub("_[^_]*$", "", a$label)), function(z) sprintf("%d/%d", sum(z), length(z)))
names(tab)[2] <- "identical"
print(tab, row.names = FALSE)
cat(sprintf("\n%d of %d rows identical to the stored output.\n", sum(same), length(same)))
