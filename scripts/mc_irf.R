#!/usr/bin/Rscript
#______________________________________________________________________________#
# Title:  MC IRF replication — Figure S.1 (ShortRun) / Figure S.2 (LongRun)
#         Tugan (2021), Econometrics Journal, doi:10.1093/ectj/utaa021
# DGP:    CommonFactors, K=2, L=1, r=1
#         Theta_1 = [0.65, 0.30; 0.20, 0.60], Sigma_e = [1, 0.5; 0.5, 1]
#
# Replicates ShortRunImpulseResponses.m / LongRunImpulseResponses.m:
#   For each replication:
#     1. pvarife estimate (n_out=1, n_in=1 to match MATLAB MC speed)
#     2. irf_bands (200 draws) -> bias-corrected median IRF + CI, normalised
#     3. OLS beta and Sigma stored
#   After loop:
#     - avg pvarife IRF = mean of per-rep median IRFs
#     - avg CI           = mean of per-rep CI bounds
#     - OLS IRF          = IRF from AVERAGE OLS beta and AVERAGE OLS Sigma
#     - True IRF         = IRF from known DGP parameters
#
# Output: csv/mc_irf_{id}_{I}_{T}.csv  and  figures/mc_irf_{id}_{I}_{T}.pdf
#
# Usage (same pattern as mc_pvarife.R):
#   Rscript mc_irf.R  <n_units> <n_time> <n_sds> <identification>
#   identification: "sr" = short-run (Fig S.1), "lr" = long-run (Fig S.2)
#   e.g.:  Rscript mc_irf.R 50 50 500 sr
#______________________________________________________________________________#
rm(list = ls())

wd      <- "/home/bc25911/MC_PVARIFE"
setwd(wd)

n_cores <- as.numeric(Sys.getenv("NSLOTS", unset = 1))
message("Using ", n_cores, " cores")

library(MonteCarlo)
library(ggplot2)

for (f in c("utils.R", "objective_fn.R", "initial_step.R",
             "estimate_pvarife.R", "asymptotic.R", "simulate.R",
             "irf.R", "confidence_bands.R")) {
  source(f, local = FALSE)
}

sessionInfo()

#______________________________________________________________________________#
# Command-line arguments
#   $1 = n_units  (I)
#   $2 = n_time   (T)
#   $3 = n_sds    (replications, paper uses 500)
#   $4 = identification: "sr" or "lr"  (default: "sr")
#______________________________________________________________________________#
cmd_args <- commandArgs(TRUE)
n_units <- as.integer(cmd_args[1])
n_time  <- as.integer(cmd_args[2])
n_sds   <- as.integer(cmd_args[3])
id_str  <- if (length(cmd_args) >= 4L) cmd_args[4] else "sr"

identification <- if (id_str == "lr") "long_run" else "short_run"
id_label       <- if (id_str == "lr") "LongRun"  else "ShortRun"
diff_vars      <- if (id_str == "lr") 1L else integer(0L)

message(sprintf("I=%d  T=%d  nrep=%d  identification=%s",
                n_units, n_time, n_sds, id_label))

# IRF settings (matching MATLAB: nir=21, n_random=200)
n_periods <- 21L
n_draw    <- 200L
shock_idx <- 1L

# True DGP parameters
# beta ordering in our package: [intercepts, Theta_1_elements]
# beta[3]=Theta11=0.65, beta[4]=Theta12=0.30, beta[5]=Theta21=0.20, beta[6]=Theta22=0.60
true_beta_r <- c(1.0, 1.0, 0.65, 0.30, 0.20, 0.60)  # our ordering
true_sigma  <- matrix(c(1.0, 0.5, 0.5, 1.0), nrow = 2L)

# Helper: OLS residual covariance (no factor correction)
# Faithful to Sigma_e_given_Beta_WITHOUT_F.m: Sigma = e*e' / (T*C)
compute_sigma_ols <- function(beta_ols, y_stack, z_stack, n_vars) {
  resid   <- as.numeric(y_stack) - as.numeric(z_stack %*% beta_ols)
  n_obs   <- length(resid) / n_vars            # total "observation-periods"
  e_mat   <- matrix(resid, nrow = n_vars, ncol = n_obs)
  tcrossprod(e_mat) / n_obs
}

# Helper: build minimal fake fit for compute_irf given any (beta, sigma)
make_irf_fit <- function(beta, sigma, n_vars, n_lags) {
  fit <- list(beta = beta, sigma = sigma,
              n_vars = n_vars, n_lags = n_lags)
  class(fit) <- "pvarife_result"
  fit
}

# True IRF (computed once, outside MC loop)
true_fit    <- make_irf_fit(true_beta_r, true_sigma, 2L, 1L)
true_ir_raw <- compute_irf(true_fit, n_periods = n_periods,
                            shock = shock_idx, diff_vars = diff_vars,
                            identification = identification)
true_ir     <- true_ir_raw / true_ir_raw[shock_idx, 1L]  # normalise

#______________________________________________________________________________#
# Monte Carlo function — one replication
# Returns:
#   irf_{v}_h{hh}: pvarife bias-corrected median IRF, normalised (v=1,2; h=00..20)
#   lo_{v}_h{hh}:  lower 95% CI, normalised
#   hi_{v}_h{hh}:  upper 95% CI, normalised
#   ols_b_{j}:     OLS beta elements (j=1..6) for averaging
#   ols_s_{ij}:    OLS Sigma unique elements (11, 12, 22) for averaging
#______________________________________________________________________________#
mc_func <- function(n_units, n_time) {

  # Source pvarife code for snow workers
  for (f in c("utils.R", "objective_fn.R", "initial_step.R",
               "estimate_pvarife.R", "asymptotic.R", "simulate.R",
               "irf.R", "confidence_bands.R"))
    source(file.path(wd, f), local = FALSE)

  n_periods_  <- 21L
  n_draw_     <- 200L
  shock_idx_  <- 1L
  diff_vars_  <- diff_vars
  ident_      <- identification

  # Build NA template (returned on failure)
  h_fmt_    <- sprintf("%02d", 0:20)
  irf_names <- c(
    paste0("irf_1_h", h_fmt_), paste0("irf_2_h", h_fmt_),
    paste0("lo_1_h",  h_fmt_), paste0("lo_2_h",  h_fmt_),
    paste0("hi_1_h",  h_fmt_), paste0("hi_2_h",  h_fmt_),
    paste0("ols_b_",  1:6),
    c("ols_s_11", "ols_s_12", "ols_s_22")
  )
  na_out <- as.list(rep(NA_real_, length(irf_names)))
  names(na_out) <- irf_names

  # ---- Simulate ---------------------------------------------------------------
  sim <- tryCatch(
    sim_pvarife(n_units=n_units, n_time=n_time,
                n_vars=2L, n_lags=1L, n_factors=1L,
                identification=ident_,    # DGP matches identification
                seed=NULL),
    error = function(e) NULL)
  if (is.null(sim)) return(na_out)

  yz <- tryCatch(.build_yz(sim$y, 1L), error=function(e) NULL)
  if (is.null(yz)) return(na_out)

  # ---- OLS (no factor correction) -------------------------------------------
  ols_res <- tryCatch(ols_nan(yz$y_stack, yz$z_stack), error=function(e) NULL)
  if (is.null(ols_res)) return(na_out)
  beta_ols  <- ols_res$beta
  sigma_ols <- compute_sigma_ols(beta_ols, yz$y_stack, yz$z_stack, 2L)

  # ---- pvarife (n_out=1, n_in=1 matching MATLAB MC) -------------------------
  fit <- tryCatch(
    pvarife(sim$y, n_lags=1L, n_factors=1L,
            n_out=1L, n_in=1L,            # MATLAB uses out_number=1, in_number=1
            balanced_init=FALSE),
    error=function(e) NULL)
  if (is.null(fit)) return(na_out)

  # ---- Confidence bands (200 draws from joint asymptotic distribution) -------
  bands <- tryCatch(
    irf_bands(fit, n_periods=n_periods_, shock=shock_idx_,
              diff_vars=diff_vars_, identification=ident_,
              n_draw=n_draw_, seed=NULL),
    error=function(e) NULL)
  if (is.null(bands)) return(na_out)

  # Normalise by median IRF of shock variable at h=0
  norm_f <- bands$irf[shock_idx_, 1L]
  if (is.na(norm_f) || abs(norm_f) < 1e-10) return(na_out)

  ir_norm  <- bands$irf   / norm_f
  lo_norm  <- bands$lower / norm_f
  hi_norm  <- bands$upper / norm_f

  # ---- Pack output -----------------------------------------------------------
  out <- c(
    as.list(ir_norm[1L, ]),  names = paste0("irf_1_h", formatC(0:20, width=2, flag="0")),
    as.list(ir_norm[2L, ]),  names = paste0("irf_2_h", formatC(0:20, width=2, flag="0")),
    as.list(lo_norm[1L, ]),  names = paste0("lo_1_h",  formatC(0:20, width=2, flag="0")),
    as.list(lo_norm[2L, ]),  names = paste0("lo_2_h",  formatC(0:20, width=2, flag="0")),
    as.list(hi_norm[1L, ]),  names = paste0("hi_1_h",  formatC(0:20, width=2, flag="0")),
    as.list(hi_norm[2L, ]),  names = paste0("hi_2_h",  formatC(0:20, width=2, flag="0")),
    as.list(as.numeric(beta_ols)),  names = paste0("ols_b_", 1:6),
    list(ols_s_11 = sigma_ols[1L, 1L],
         ols_s_12 = sigma_ols[1L, 2L],
         ols_s_22 = sigma_ols[2L, 2L])
  )

  # Consistent zero-padded column names: h00..h20
  h_fmt <- sprintf("%02d", 0:20)
  result <- as.list(c(
    setNames(as.numeric(ir_norm[1L,]), paste0("irf_1_h", h_fmt)),
    setNames(as.numeric(ir_norm[2L,]), paste0("irf_2_h", h_fmt)),
    setNames(as.numeric(lo_norm[1L,]), paste0("lo_1_h",  h_fmt)),
    setNames(as.numeric(lo_norm[2L,]), paste0("lo_2_h",  h_fmt)),
    setNames(as.numeric(hi_norm[1L,]), paste0("hi_1_h",  h_fmt)),
    setNames(as.numeric(hi_norm[2L,]), paste0("hi_2_h",  h_fmt)),
    setNames(as.numeric(beta_ols),     paste0("ols_b_",  1:6)),
    c(ols_s_11 = sigma_ols[1L,1L], ols_s_12 = sigma_ols[1L,2L],
      ols_s_22 = sigma_ols[2L,2L])
  ))
  result
}

environment(mc_func)$wd           <- wd
environment(mc_func)$diff_vars    <- diff_vars
environment(mc_func)$identification <- identification

#______________________________________________________________________________#
# Run Monte Carlo
#______________________________________________________________________________#
message("Starting MonteCarlo...")
mc_result <- MonteCarlo(
  func       = mc_func,
  nrep       = n_sds,
  param_list = list("n_units" = n_units, "n_time" = n_time),
  ncpus      = n_cores
)
summary(mc_result)

#______________________________________________________________________________#
# Save raw MonteCarlo output
#______________________________________________________________________________#
dir.create("csv",     showWarnings=FALSE, recursive=TRUE)
dir.create("figures", showWarnings=FALSE, recursive=TRUE)

res_df  <- MakeFrame(mc_result)
csv_out <- sprintf("csv/mc_irf_%s_I%d_T%d.csv", id_str, n_units, n_time)
write.csv(res_df, file=csv_out, row.names=FALSE)
message("Raw results saved: ", csv_out)

#______________________________________________________________________________#
# Compute averages across replications
#______________________________________________________________________________#
get_avg <- function(nm) mean(as.numeric(mc_result$results[[nm]]), na.rm=TRUE)

horizons <- 0:20

# Average pvarife IRF and bands (avg of per-rep normalised values)
avg_irf <- matrix(NA_real_, nrow=2L, ncol=21L)
avg_lo  <- matrix(NA_real_, nrow=2L, ncol=21L)
avg_hi  <- matrix(NA_real_, nrow=2L, ncol=21L)
for (vv in 1:2) {
  for (hh in horizons) {
    h_str           <- formatC(hh, width=2L, flag="0")
    avg_irf[vv, hh+1L] <- get_avg(paste0("irf_",vv,"_h",h_str))
    avg_lo [vv, hh+1L] <- get_avg(paste0("lo_", vv,"_h",h_str))
    avg_hi [vv, hh+1L] <- get_avg(paste0("hi_", vv,"_h",h_str))
  }
}

# Average OLS beta and Sigma -> compute OLS IRF (matching MATLAB approach)
avg_beta_ols  <- sapply(paste0("ols_b_", 1:6), get_avg)
avg_sigma_ols <- matrix(c(get_avg("ols_s_11"), get_avg("ols_s_12"),
                           get_avg("ols_s_12"), get_avg("ols_s_22")), nrow=2L)

ols_fit   <- make_irf_fit(avg_beta_ols, avg_sigma_ols, 2L, 1L)
ols_ir_raw <- tryCatch(
  compute_irf(ols_fit, n_periods=n_periods, shock=shock_idx,
              diff_vars=diff_vars, identification=identification),
  error=function(e) NULL)
if (!is.null(ols_ir_raw)) {
  ols_ir <- ols_ir_raw / ols_ir_raw[shock_idx, 1L]
} else {
  ols_ir <- matrix(NA_real_, nrow=2L, ncol=n_periods)
}

n_valid <- sum(!is.na(as.numeric(mc_result$results$irf_1_h00)))
message(sprintf("Valid replications: %d / %d", n_valid, n_sds))

#______________________________________________________________________________#
# Print summary
#______________________________________________________________________________#
cat(sprintf("\n=== IRF Summary: %s  I=%d  T=%d  valid=%d ===\n",
            id_label, n_units, n_time, n_valid))
cat(sprintf("%-8s %6s %6s %6s %6s %6s\n",
            "Horizon","TrIRF1","AvIRF1","OLSIRF1","TrIRF2","AvIRF2"))
for (hh in c(0,1,2,4,8,12,20)) {
  h_idx <- hh + 1L
  cat(sprintf("h=%-4d  %6.3f %6.3f %6.3f  %6.3f %6.3f\n",
              hh, true_ir[1,h_idx], avg_irf[1,h_idx],
              ifelse(is.null(ols_ir), NA, ols_ir[1,h_idx]),
              true_ir[2,h_idx], avg_irf[2,h_idx]))
}

#______________________________________________________________________________#
# Plot (replicates Figure S.1 / S.2 style)
#______________________________________________________________________________#
var_labels <- c("Variable 1", "Variable 2")

build_df <- function(irf_mat, label, type) {
  data.frame(
    variable = rep(var_labels, each=21L),
    horizon  = rep(horizons,   times=2L),
    value    = as.numeric(t(irf_mat)),
    type     = label
  )
}

df_plot <- rbind(
  build_df(true_ir, "True IRF",       "true"),
  build_df(avg_irf, "pvarife (avg)",  "est"),
  build_df(ols_ir,  "OLS (avg)",      "ols")
)
df_plot$variable <- factor(df_plot$variable, levels=var_labels)
df_plot$type     <- factor(df_plot$type,
                            levels=c("true","est","ols"),
                            labels=c("True IRF","pvarife (avg)","OLS (avg)"))

# CI band (avg of per-rep bands)
df_band <- data.frame(
  variable = rep(rep(var_labels, each=21L), times=1L),
  horizon  = rep(horizons, times=2L),
  lo       = as.numeric(t(avg_lo)),
  hi       = as.numeric(t(avg_hi))
)
df_band$variable <- factor(df_band$variable, levels=var_labels)

p <- ggplot(df_band, aes(x=horizon)) +
  geom_ribbon(aes(ymin=lo, ymax=hi), fill="grey75", alpha=0.5) +
  geom_hline(yintercept=0, linewidth=0.3, colour="black") +
  geom_line(data=df_plot,
            aes(x=horizon, y=value, linetype=type, shape=type),
            linewidth=0.9) +
  geom_point(data=df_plot[df_plot$type %in% c("True IRF","pvarife (avg)"),],
             aes(x=horizon, y=value, shape=type), size=2L) +
  scale_linetype_manual(values=c("True IRF"="solid","pvarife (avg)"="dashed","OLS (avg)"="dashed")) +
  scale_shape_manual(   values=c("True IRF"=8L,    "pvarife (avg)"=3L,      "OLS (avg)"=NA)) +
  facet_wrap(~variable, scales="free_y", ncol=2L) +
  scale_x_continuous(breaks=seq(0,20,4)) +
  labs(
    title    = sprintf("Average IRFs — %s identification (%s)",
                       id_label, sprintf("I=%d, T=%d, n=%d", n_units, n_time, n_sds)),
    subtitle = "Shaded: avg 95% CI (pvarife) | Lines: True / pvarife / OLS",
    x="Period", y="Normalised response",
    linetype=NULL, shape=NULL
  ) +
  theme_minimal(base_size=11L) +
  theme(panel.grid.minor=element_blank(),
        strip.text=element_text(face="bold"),
        legend.position="bottom")

fig_out <- sprintf("figures/mc_irf_%s_I%d_T%d.pdf", id_str, n_units, n_time)
ggsave(fig_out, p, width=9L, height=4.5)
message("Figure saved: ", fig_out)
message("Done.")
