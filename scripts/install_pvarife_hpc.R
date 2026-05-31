#!/usr/bin/env Rscript
# =============================================================================
# Install pvarife on HPC — run ONCE before submitting any jobs
# Usage: Rscript scripts/install_pvarife_hpc.R
# =============================================================================

# Install to user library (no root needed)
user_lib <- Sys.getenv("R_LIBS_USER")
if (!nchar(user_lib)) user_lib <- "~/R/library"
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))

cat("Installing dependencies...\n")
pkgs_needed <- c("future", "furrr", "mvtnorm", "ggplot2", "rlang")
new_pkgs <- pkgs_needed[!sapply(pkgs_needed, requireNamespace, quietly = TRUE)]
if (length(new_pkgs) > 0L) {
  install.packages(new_pkgs, lib = user_lib,
                   repos = "https://cloud.r-project.org", quiet = FALSE)
}

cat("\nBuilding pvarife package tarball...\n")
# Build the package tarball from source (run from repo root)
pkg_tar <- devtools::build(quiet = FALSE)   # produces pvarife_0.1.0.tar.gz

cat("\nInstalling pvarife...\n")
install.packages(pkg_tar, repos = NULL, type = "source",
                 lib = user_lib, quiet = FALSE)

cat("\nVerification:\n")
library(pvarife, lib.loc = user_lib)
sim <- sim_pvarife(n_units = 5L, n_time = 8L, seed = 1L)
fit <- pvarife(sim$y, n_lags = 1L, n_factors = 1L, n_out = 2L, n_in = 2L)
cat("pvarife installed and working. beta length:", length(fit$beta), "\n")
