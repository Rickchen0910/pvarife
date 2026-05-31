#!/usr/bin/Rscript
#______________________________________________________________________________#
# Title:  MC convergence check — Panel VAR with Interactive Fixed Effects
#         Tugan (2021), Econometrics Journal, doi:10.1093/ectj/utaa021
# DGP:    CommonFactors / ShortRun, K=2, L=1, r=1
#         Theta_1 = [0.65, 0.3; 0.2, 0.6], Sigma_e = [1, 0.5; 0.5, 1]
# Output: csv/mc_pvarife_I{n}_T{t}.csv
#
# FIX: MonteCarlo uses ls(parent.frame()) for clusterExport, which silently
#      drops ALL dot-prefixed functions (.build_yz, .initial_step, etc.).
#      Solution: source() inside mc_func so each snow worker loads everything.
#______________________________________________________________________________#
rm(list = ls())

# ---- Working directory (adjust to your HPC path) ----------------------------
wd <- "/home/bc25911/MC_PVARIFE"
setwd(wd)

# ---- Number of cores (from SGE -pe smp) -------------------------------------
n_cores <- as.numeric(Sys.getenv("NSLOTS", unset = 1))
message("Using ", n_cores, " cores")

# ---- External packages -------------------------------------------------------
library(MonteCarlo)

# ---- Source pvarife code into main session ----------------------------------
# (also sourced inside mc_func for snow workers — see below)
pvarife_files <- c("utils.R", "objective_fn.R", "initial_step.R",
                   "estimate_pvarife.R", "asymptotic.R", "simulate.R")
for (f in pvarife_files) source(f, local = FALSE)

sessionInfo()

#______________________________________________________________________________#
# Command-line arguments:  $1=n_units  $2=n_time  $3=n_sds
#______________________________________________________________________________#
cmd_args <- commandArgs(TRUE)
n_units  <- as.integer(cmd_args[1])
n_time   <- as.integer(cmd_args[2])
n_sds    <- as.integer(cmd_args[3])

message(sprintf("I=%d  T=%d  nrep=%d", n_units, n_time, n_sds))

# True parameter values (DGP matches sim_pvarife())
# Z matrix ordering: beta = [c1, c2, Theta11, Theta12, Theta21, Theta22]
#   beta[3] = Theta11 = 0.65  (coeff of y1_{t-1} in equation for y1)
#   beta[4] = Theta12 = 0.30  (coeff of y2_{t-1} in equation for y1)  <-- NOT 0.20
#   beta[5] = Theta21 = 0.20  (coeff of y1_{t-1} in equation for y2)  <-- NOT 0.30
#   beta[6] = Theta22 = 0.60  (coeff of y2_{t-1} in equation for y2)
# Verified by: sim_pvarife(seed=1)$beta_true = c(1, 1, 0.65, 0.30, 0.20, 0.60)
true_beta <- c(1.0, 1.0, 0.65, 0.30, 0.20, 0.60)

#______________________________________________________________________________#
# Simulation function
#______________________________________________________________________________#
mc_func <- function(n_units, n_time) {

  # ---- Load pvarife code in this worker's session ---------------------------
  for (f in c("utils.R", "objective_fn.R", "initial_step.R",
               "estimate_pvarife.R", "asymptotic.R", "simulate.R")) {
    source(file.path(wd, f), local = FALSE)
  }

  message(sprintf("[PID %s] MC rep | I=%d | T=%d", Sys.getpid(), n_units, n_time))

  # Return-value template (all NA on failure)
  na_out <- list(
    bias_theta11 = NA_real_, bias_theta21 = NA_real_,
    bias_theta12 = NA_real_, bias_theta22 = NA_real_,
    sq_theta11   = NA_real_, sq_theta21   = NA_real_,
    sq_theta12   = NA_real_, sq_theta22   = NA_real_,
    cover_theta11 = NA_real_, cover_theta21 = NA_real_,
    cover_theta12 = NA_real_, cover_theta22 = NA_real_
  )

  # ---- Simulate (seed=NULL: let MonteCarlo control RNG) ----------------------
  sim <- tryCatch(
    sim_pvarife(n_units = n_units, n_time = n_time,
                n_vars = 2L, n_lags = 1L, n_factors = 1L,
                seed = NULL),
    error = function(e) {
      message("sim_pvarife() failed: ", conditionMessage(e)); NULL
    }
  )
  if (is.null(sim)) return(na_out)

  # ---- Estimate --------------------------------------------------------------
  fit <- tryCatch(
    pvarife(sim$y,
            n_lags        = 1L,
            n_factors     = 1L,
            n_out         = 50L,
            n_in          = 10L,
            balanced_init = FALSE),
    error = function(e) {
      message("pvarife() failed: ", conditionMessage(e)); NULL
    }
  )
  if (is.null(fit)) return(na_out)

  # ---- Bias (Theta_1 elements = beta[3:6]) -----------------------------------
  # beta[4]=Theta12=0.30, beta[5]=Theta21=0.20  (see true_beta comment above)
  true_b      <- c(1.0, 1.0, 0.65, 0.30, 0.20, 0.60)
  beta_est    <- as.numeric(fit$beta)
  bias_all    <- beta_est - true_b
  bias_theta  <- bias_all[3:6]     # Theta11, Theta21, Theta12, Theta22

  # ---- Asymptotic 95% CI (Theorem 2.3) --------------------------------------
  avar <- tryCatch(asymptotic_var(fit), error = function(e) NULL)
  cover_theta <- rep(NA_real_, 4L)
  if (!is.null(avar)) {
    se      <- sqrt(pmax(0.0, diag(avar$variance)))
    b_bc    <- beta_est - avar$bias      # bias-corrected centre
    ci_lo   <- b_bc - 1.96 * se
    ci_hi   <- b_bc + 1.96 * se
    # cover_theta[1]=Theta11, [2]=Theta12, [3]=Theta21, [4]=Theta22
    cover_theta <- as.numeric(
      (true_b[3:6] >= ci_lo[3:6]) & (true_b[3:6] <= ci_hi[3:6])
    )
  }

  # Naming: beta[4]=Theta12, beta[5]=Theta21 (not the other way around)
  list(
    bias_theta11  = bias_theta[1],
    bias_theta12  = bias_theta[2],   # beta[4]: coeff of y2 in eq1
    bias_theta21  = bias_theta[3],   # beta[5]: coeff of y1 in eq2
    bias_theta22  = bias_theta[4],
    sq_theta11    = bias_theta[1]^2,
    sq_theta12    = bias_theta[2]^2,
    sq_theta21    = bias_theta[3]^2,
    sq_theta22    = bias_theta[4]^2,
    cover_theta11 = cover_theta[1],
    cover_theta12 = cover_theta[2],
    cover_theta21 = cover_theta[3],
    cover_theta22 = cover_theta[4]
  )
}

# Make wd available in mc_func's closure so workers can source files
environment(mc_func)$wd <- wd

#______________________________________________________________________________#
# Run Monte Carlo
#______________________________________________________________________________#
param_list <- list("n_units" = n_units, "n_time" = n_time)

mc_result <- MonteCarlo(
  func       = mc_func,
  nrep       = n_sds,
  param_list = param_list,
  ncpus      = n_cores
)

summary(mc_result)

#______________________________________________________________________________#
# Save raw results
#______________________________________________________________________________#
dir.create("csv", showWarnings = FALSE, recursive = TRUE)

res_df   <- MakeFrame(mc_result)
outfile  <- sprintf("csv/mc_pvarife_I%d_T%d.csv", n_units, n_time)
write.csv(res_df, file = outfile, row.names = FALSE)
message("Saved: ", outfile)

#______________________________________________________________________________#
# Print summary to log file
#______________________________________________________________________________#
get_mean <- function(nm) mean(as.numeric(mc_result$results[[nm]]), na.rm = TRUE)
n_valid  <- sum(!is.na(as.numeric(mc_result$results$bias_theta11)))

# Order: Theta11(0.65), Theta12(0.30), Theta21(0.20), Theta22(0.60)
cat(sprintf("\n=== Summary: I=%d, T=%d, nrep=%d, valid=%d ===\n",
            n_units, n_time, n_sds, n_valid))
cat(sprintf("True   Theta11= 0.6500  Theta12= 0.3000  Theta21= 0.2000  Theta22= 0.6000\n"))
cat(sprintf("Bias   Theta11=%7.4f  Theta12=%7.4f  Theta21=%7.4f  Theta22=%7.4f\n",
            get_mean("bias_theta11"), get_mean("bias_theta12"),
            get_mean("bias_theta21"), get_mean("bias_theta22")))
cat(sprintf("RMSE   Theta11=%7.4f  Theta12=%7.4f  Theta21=%7.4f  Theta22=%7.4f\n",
            sqrt(get_mean("sq_theta11")), sqrt(get_mean("sq_theta12")),
            sqrt(get_mean("sq_theta21")), sqrt(get_mean("sq_theta22"))))
cat(sprintf("Cov95  Theta11=%7.3f  Theta12=%7.3f  Theta21=%7.3f  Theta22=%7.3f\n",
            get_mean("cover_theta11"), get_mean("cover_theta12"),
            get_mean("cover_theta21"), get_mean("cover_theta22")))
