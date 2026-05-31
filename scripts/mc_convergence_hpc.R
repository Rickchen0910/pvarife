#!/usr/bin/env Rscript
# =============================================================================
# Monte Carlo Convergence Check — HPC version (no package install needed)
# Tugan (2021) DGP: CommonFactors / ShortRun, K=2, L=1, r=1
#
# Upload to HPC folder:
#   utils.R, objective_fn.R, initial_step.R, estimate_pvarife.R,
#   asymptotic.R, simulate.R, mc_convergence_hpc.R, mc_array.sh
#
# Run (SLURM array, one task per grid point):
#   sbatch mc_array.sh
#
# Run (single node, all grid points):
#   Rscript mc_convergence_hpc.R --ncores 16 --task 0 --nsds 500
#
# SLURM array (each task = one (I,T) combination):
#   Rscript mc_convergence_hpc.R --ncores 16 --task $SLURM_ARRAY_TASK_ID --nsds 500
# =============================================================================

# ---- Source all pvarife code (replaces library(pvarife)) --------------------
# These files must be in the same folder as this script.
pvarife_files <- c(
  "utils.R",
  "objective_fn.R",
  "initial_step.R",
  "estimate_pvarife.R",
  "asymptotic.R",
  "simulate.R"
)

# Get the directory containing this script so source() works from any working dir
script_dir <- tryCatch(
  dirname(sys.frame(1)$ofile),   # when run via Rscript
  error = function(e) "."        # fallback: current working dir
)

for (f in pvarife_files) {
  fp <- file.path(script_dir, f)
  if (!file.exists(fp)) stop("Required file not found: ", fp)
  source(fp, local = FALSE)       # load into global env
}

# ---- External packages (must be installed on the HPC) -----------------------
library(MonteCarlo)   # parallel MC framework
library(ggplot2)      # for plots

# ---- Parse command-line arguments -------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, env_var = "", default) {
  idx <- which(args == flag)
  if (length(idx) > 0L && idx[1L] + 1L <= length(args))
    return(as.numeric(args[idx[1L] + 1L]))
  env_val <- Sys.getenv(env_var, unset = "")
  if (nchar(env_val) > 0L) return(as.numeric(env_val))
  default
}

task_id <- as.integer(get_arg("--task",   "SLURM_ARRAY_TASK_ID", 0))
n_cores <- as.integer(get_arg("--ncores", "SLURM_CPUS_PER_TASK", 1))
n_sds   <- as.integer(get_arg("--nsds",   "",                    500))
out_dir <- {
  idx <- which(args == "--outdir")
  if (length(idx) > 0L && idx[1L] + 1L <= length(args))
    args[idx[1L] + 1L]
  else
    file.path(script_dir, "mc_hpc_results")
}

# ---- Settings ---------------------------------------------------------------
n_lags   <- 1L
n_factors <- 1L
n_vars   <- 2L
n_out    <- 50L     # Paper: 50 outer iterations
n_in     <- 10L     # Paper: 10 inner iterations

# True parameters matching sim_pvarife() DGP (beta = [intercepts, Theta_1])
true_beta <- c(1.0, 1.0,          # intercepts c = [1; 1]
               0.65, 0.20,         # Theta_1 row 1
               0.30, 0.60)         # Theta_1 row 2

# (I, T) grid — 9 combinations matching paper Table S10.1
full_grid <- expand.grid(
  n_units = c(25L, 50L, 100L),
  n_time  = c(25L, 50L, 100L)
)

# If task_id > 0: run only that single grid row (SLURM array mode)
if (task_id > 0L) {
  if (task_id > nrow(full_grid))
    stop("task_id ", task_id, " exceeds grid size ", nrow(full_grid))
  run_grid <- full_grid[task_id, , drop = FALSE]
  cat(sprintf("=== SLURM Array Task %d: I=%d, T=%d ===\n",
              task_id, run_grid$n_units, run_grid$n_time))
} else {
  run_grid <- full_grid
  cat("=== Running all", nrow(run_grid), "grid points ===\n")
}

cat(sprintf("n_sds=%d, n_out=%d, n_in=%d, n_cores=%d\n\n",
            n_sds, n_out, n_in, n_cores))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Monte Carlo function ---------------------------------------------------
# IMPORTANT: snow workers (ncpus > 1) are fresh R sessions.
# All pvarife functions must be sourced INSIDE mc_func so each worker gets them.
# We hard-code the file list to avoid closure-capture issues.

mc_func <- function(n_units, n_time) {

  # Source pvarife code into this worker's global environment
  # (script_dir is captured from the enclosing environment below)
  for (f in c("utils.R","objective_fn.R","initial_step.R",
               "estimate_pvarife.R","asymptotic.R","simulate.R")) {
    fp <- file.path(script_dir, f)
    source(fp, local = FALSE)
  }

  # Simulation (no fixed seed — MonteCarlo handles RNG)
  sim <- sim_pvarife(
    n_units   = n_units,
    n_time    = n_time,
    n_vars    = 2L,
    n_lags    = 1L,
    n_factors = 1L,
    seed      = NULL     # must be NULL so MonteCarlo's RNG controls draws
  )

  # Estimation
  fit <- tryCatch(
    pvarife(sim$y,
            n_lags        = 1L,
            n_factors     = 1L,
            n_out         = 50L,
            n_in          = 10L,
            balanced_init = FALSE),   # balanced data: skip subsample search
    error = function(e) NULL
  )

  # Return NA list on failure (MonteCarlo will mark as missing)
  na_list <- as.list(rep(NA_real_, 12L))
  names(na_list) <- c(paste0("bias_b",  1:6),
                      paste0("cover_b", 1:6))
  if (is.null(fit)) return(na_list)

  beta_est <- as.numeric(fit$beta)
  true_b   <- c(1.0, 1.0, 0.65, 0.20, 0.30, 0.60)
  bias_vec <- beta_est - true_b

  # Asymptotic 95% CI (Theorem 2.3)
  avar <- tryCatch(asymptotic_var(fit), error = function(e) NULL)
  cover_vec <- rep(NA_real_, 6L)
  if (!is.null(avar)) {
    se     <- sqrt(pmax(0.0, diag(avar$variance)))
    b_bc   <- beta_est - avar$bias       # bias-corrected centre
    ci_lo  <- b_bc - 1.96 * se
    ci_hi  <- b_bc + 1.96 * se
    cover_vec <- as.numeric((true_b >= ci_lo) & (true_b <= ci_hi))
  }

  result <- as.list(c(bias_vec, cover_vec))
  names(result) <- c(paste0("bias_b",  1:6),
                     paste0("cover_b", 1:6))
  result
}

# Make script_dir available inside mc_func's closure
environment(mc_func)$script_dir <- script_dir

# ---- Run Monte Carlo for each grid point ------------------------------------
all_results <- vector("list", nrow(run_grid))

for (gg in seq_len(nrow(run_grid))) {
  nn_units <- run_grid$n_units[gg]
  nn_time  <- run_grid$n_time[gg]

  cat(sprintf("[%d/%d] I=%d, T=%d — running %d sims on %d core(s) ...\n",
              gg, nrow(run_grid), nn_units, nn_time, n_sds, n_cores))
  t0 <- proc.time()

  mc_res <- MonteCarlo(
    func       = mc_func,
    nrep       = n_sds,
    param_list = list(n_units = nn_units, n_time = nn_time),
    ncpus      = n_cores
  )

  elapsed <- (proc.time() - t0)["elapsed"]
  cat(sprintf("  done in %.0fs\n", elapsed))

  # ---- Extract results from MonteCarlo array --------------------------------
  # mc_res$results$bias_b3 is a 3D array [1, 1, nrep] for one param combo
  extract_vec <- function(name) {
    arr <- mc_res$results[[name]]
    as.numeric(arr)          # flatten: length = nrep
  }

  bias_mat  <- sapply(paste0("bias_b",  1:6), extract_vec)  # nrep x 6
  cover_mat <- sapply(paste0("cover_b", 1:6), extract_vec)  # nrep x 6

  valid     <- !is.na(bias_mat[, 1L])
  n_valid   <- sum(valid)

  bias_mean <- colMeans(bias_mat[valid, , drop = FALSE])
  rmse      <- sqrt(colMeans(bias_mat[valid, , drop = FALSE]^2L))
  cov_rate  <- colMeans(cover_mat[valid, , drop = FALSE], na.rm = TRUE)

  cat(sprintf("  n_valid=%d | Theta_1 bias: %.4f %.4f %.4f %.4f\n",
              n_valid, bias_mean[3L], bias_mean[4L],
              bias_mean[5L], bias_mean[6L]))
  cat(sprintf("  Theta_1 RMSE: %.4f %.4f %.4f %.4f\n",
              rmse[3L], rmse[4L], rmse[5L], rmse[6L]))
  cat(sprintf("  Theta_1 coverage (95%% CI): %.3f %.3f %.3f %.3f\n",
              cov_rate[3L], cov_rate[4L], cov_rate[5L], cov_rate[6L]))

  all_results[[gg]] <- data.frame(
    n_units  = nn_units,
    n_time   = nn_time,
    n_valid  = n_valid,
    param    = paste0("beta[", 1:6, "]"),
    true_val = c(1.0, 1.0, 0.65, 0.20, 0.30, 0.60),
    bias     = bias_mean,
    rmse     = rmse,
    coverage = cov_rate
  )

  # Save results after each grid point (so partial results survive if job killed)
  partial <- do.call(rbind, all_results[!sapply(all_results, is.null)])
  out_rds <- if (task_id > 0L)
    file.path(out_dir, sprintf("mc_task%02d.rds", task_id))
  else
    file.path(out_dir, "mc_results.rds")
  saveRDS(partial, out_rds)
}

res <- do.call(rbind, all_results)
cat("\n\n=== FINAL SUMMARY ===\n")

# ---- Summary tables ---------------------------------------------------------
theta_params <- paste0("beta[", 3:6, "]")
theta_names  <- c("Theta[1,1]=0.65","Theta[2,1]=0.20",
                  "Theta[1,2]=0.30","Theta[2,2]=0.60")

res_theta <- res[res$param %in% theta_params, ]

cat("\nBias (should shrink toward 0 as I and T increase):\n")
for (ii in seq_along(theta_params)) {
  pp  <- theta_params[ii]
  sub <- res_theta[res_theta$param == pp, ]
  cat(sprintf("\n  %s:\n", theta_names[ii]))
  tbl <- reshape(sub[, c("n_units","n_time","bias")],
                 idvar = "n_units", timevar = "n_time", direction = "wide")
  print(round(tbl, 4))
}

cat("\nRMSE (should decrease with I and T):\n")
for (ii in seq_along(theta_params)) {
  pp  <- theta_params[ii]
  sub <- res_theta[res_theta$param == pp, ]
  cat(sprintf("\n  %s:\n", theta_names[ii]))
  tbl <- reshape(sub[, c("n_units","n_time","rmse")],
                 idvar = "n_units", timevar = "n_time", direction = "wide")
  print(round(tbl, 4))
}

cat("\nCoverage rate of 95% CI (target: 0.950):\n")
for (ii in seq_along(theta_params)) {
  pp  <- theta_params[ii]
  sub <- res_theta[res_theta$param == pp, ]
  cat(sprintf("\n  %s:\n", theta_names[ii]))
  tbl <- reshape(sub[, c("n_units","n_time","coverage")],
                 idvar = "n_units", timevar = "n_time", direction = "wide")
  print(round(tbl, 3))
}

# ---- Convergence diagnosis --------------------------------------------------
cat("\nConvergence check (RMSE should decrease as I*T grows):\n")
for (ii in seq_along(theta_params)) {
  pp   <- theta_params[ii]
  sub  <- res_theta[res_theta$param == pp, ]
  sub  <- sub[order(sub$n_units * sub$n_time), ]
  mono <- all(diff(sub$rmse) <= 0.02)   # allow small MC noise
  cat(sprintf("  %s: [%s] — %s\n",
              theta_names[ii],
              paste(round(sub$rmse, 4), collapse = ", "),
              if (mono) "CONVERGING OK" else "check manually"))
}

# ---- Plot -------------------------------------------------------------------
res_plot <- res_theta
res_plot$T_label <- factor(res_plot$n_time,
                            levels = c(25, 50, 100),
                            labels = c("T=25", "T=50", "T=100"))
res_plot$I_label <- factor(res_plot$n_units)
res_plot$param_label <- factor(
  res_plot$param,
  levels = theta_params,
  labels = theta_names
)

p_rmse <- ggplot(res_plot,
                 aes(x = I_label, y = rmse,
                     colour = T_label, group = T_label)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  facet_wrap(~ param_label, scales = "free_y", ncol = 2L) +
  labs(
    title    = "RMSE of VAR coefficient estimates",
    subtitle = sprintf("pvarife R package — Tugan (2021) DGP, n_sds=%d, n_out=%d",
                       n_sds, n_out),
    x        = "Units (I)",
    y        = "RMSE",
    colour   = "Time periods"
  ) +
  theme_minimal(base_size = 11L) +
  theme(panel.grid.minor  = element_blank(),
        strip.text        = element_text(face = "bold"),
        legend.position   = "bottom")

pdf_path <- file.path(out_dir, sprintf("mc_rmse%s.pdf",
                                        if (task_id > 0L)
                                          sprintf("_task%02d", task_id)
                                        else ""))
ggsave(pdf_path, p_rmse, width = 8L, height = 5L)
cat(sprintf("\nPlot saved: %s\n", pdf_path))

rds_path <- if (task_id > 0L)
  file.path(out_dir, sprintf("mc_task%02d.rds", task_id))
else
  file.path(out_dir, "mc_results.rds")
saveRDS(res, rds_path)
cat(sprintf("Results saved: %s\n", rds_path))
cat("\nDone.\n")
