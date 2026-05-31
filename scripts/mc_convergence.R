# =============================================================================
# Monte Carlo Convergence Check: Tugan (2021) DGP
# Checks that pvarife() recovers the true beta as (I, T) grows.
#
# Faithfully matches the paper's DGP (Section S10, CommonFactors / ShortRun):
#   Theta_1 = [0.65, 0.3; 0.2, 0.6], Sigma_e = [1, 0.5; 0.5, 1]
#   f_t ~ AR(1) with coeff 0.5, variance 0.5; lambda_{k,i} ~ N(1,1)
#
# Outputs:
#   - Bias and RMSE of each beta element at each (I, T) combination
#   - Coverage rates of 95% CI (using asymptotic_var)
#   - Saved to scripts/mc_results.rds and scripts/mc_convergence.pdf
#
# Runtime: ~5-30 min depending on n_sds and n_out
# Adjust n_sds (simulations), n_out, n_in for a quick test run.
# =============================================================================

devtools::load_all(quiet = TRUE)
library(ggplot2)

# ---- Settings ---------------------------------------------------------------
n_sds    <- 500L   # simulations per (I,T). Paper uses 500; use 100 for quick run.
n_lags   <- 1L
n_factors <- 1L
n_vars   <- 2L
n_out    <- 50L    # Paper uses 50; reduce for speed
n_in     <- 10L     # Paper uses 10
level    <- 0.95
seed0    <- 2026L

# (I, T) grid matching the paper's Table S10.1: I in {25,50,100}, T in {25,50,100}
grid <- expand.grid(
  n_units = c(25L, 50L, 100L),
  n_time  = c(25L, 50L, 100L)
)

# True parameters
true_beta <- c(1.0, 1.0,   # intercepts (c = [1;1])
               0.65, 0.20, 0.30, 0.60)  # Theta_1 row-major (matches sim_pvarife)

cat("=== Monte Carlo Convergence Check ===\n")
cat(sprintf("n_sds=%d, n_out=%d, n_in=%d\n\n", n_sds, n_out, n_in))

# ---- Run MC -----------------------------------------------------------------
results <- vector("list", nrow(grid))

for (gg in seq_len(nrow(grid))) {
  nn_units <- grid$n_units[gg]
  nn_time  <- grid$n_time[gg]
  set.seed(seed0 + gg * 100L)

  cat(sprintf("[%d/%d] I=%d, T=%d ... ", gg, nrow(grid), nn_units, nn_time))
  t0 <- proc.time()

  beta_mat   <- matrix(NA_real_, nrow = n_sds, ncol = length(true_beta))
  covered    <- matrix(NA_real_, nrow = n_sds, ncol = length(true_beta))

  for (ss in seq_len(n_sds)) {
    sim <- sim_pvarife(
      n_units   = nn_units,
      n_time    = nn_time,
      n_vars    = n_vars,
      n_lags    = n_lags,
      n_factors = n_factors,
      seed      = seed0 + gg * 100L + ss
    )

    fit <- tryCatch(
      pvarife(sim$y, n_lags = n_lags, n_factors = n_factors,
              n_out = n_out, n_in = n_in, balanced_init = FALSE),
      error = function(e) NULL
    )
    if (is.null(fit)) next

    beta_mat[ss, ] <- as.numeric(fit$beta)

    # Asymptotic CI (Theorem 2.3)
    avar <- tryCatch(asymptotic_var(fit), error = function(e) NULL)
    if (!is.null(avar)) {
      se     <- sqrt(pmax(0, diag(avar$variance)))
      b_bc   <- as.numeric(fit$beta) - avar$bias
      ci_lo  <- b_bc - qnorm(1 - (1 - level) / 2) * se
      ci_hi  <- b_bc + qnorm(1 - (1 - level) / 2) * se
      covered[ss, ] <- (true_beta >= ci_lo) & (true_beta <= ci_hi)
    }
  }

  elapsed <- (proc.time() - t0)["elapsed"]
  n_valid <- sum(!is.na(beta_mat[, 1L]))
  bias    <- colMeans(beta_mat - rep(1, n_valid), na.rm = TRUE)
  # recompute properly:
  bias    <- colMeans(sweep(beta_mat, 2L, true_beta, "-"), na.rm = TRUE)
  rmse    <- sqrt(colMeans(sweep(beta_mat, 2L, true_beta, "-")^2, na.rm = TRUE))
  cov_rate <- colMeans(covered, na.rm = TRUE)

  results[[gg]] <- data.frame(
    n_units   = nn_units,
    n_time    = nn_time,
    n_valid   = n_valid,
    param     = paste0("beta[", seq_along(true_beta), "]"),
    true_val  = true_beta,
    mean_est  = colMeans(beta_mat, na.rm = TRUE),
    bias      = bias,
    rmse      = rmse,
    coverage  = cov_rate
  )
  cat(sprintf("done (%ds, n_valid=%d)\n", round(elapsed), n_valid))
  cat(sprintf("  Bias (Theta_1 elements): %.4f %.4f %.4f %.4f\n",
              bias[3], bias[4], bias[5], bias[6]))
}

# ---- Compile results --------------------------------------------------------
res <- do.call(rbind, results)
saveRDS(res, "scripts/mc_results.rds")

# ---- Print summary table ----------------------------------------------------
cat("\n=== Bias of Theta_1 elements (should shrink with I and T) ===\n")
theta_params <- c("beta[3]","beta[4]","beta[5]","beta[6]")
res_theta <- res[res$param %in% theta_params, ]

for (pp in theta_params) {
  sub <- res_theta[res_theta$param == pp, ]
  cat(sprintf("\n%s (true=%.2f):\n", pp, unique(sub$true_val)))
  tbl <- reshape(sub[, c("n_units","n_time","bias","rmse")],
                 idvar = "n_units", timevar = "n_time", direction = "wide")
  print(round(tbl, 4))
}

cat("\n=== Coverage rates of 95% CI ===\n")
for (pp in theta_params) {
  sub <- res_theta[res_theta$param == pp, ]
  cat(sprintf("\n%s:\n", pp))
  tbl <- reshape(sub[, c("n_units","n_time","coverage")],
                 idvar = "n_units", timevar = "n_time", direction = "wide")
  print(round(tbl, 3))
}

# ---- Plot: RMSE grid --------------------------------------------------------
cat("\nSaving convergence plot to scripts/mc_convergence.pdf ...\n")

res_plot <- res[res$param %in% theta_params, ]
res_plot$IT <- paste0("I=",res_plot$n_units,", T=",res_plot$n_time)
res_plot$label <- paste0("I=", res_plot$n_units)
res_plot$n_time_f <- factor(res_plot$n_time,
                             levels = c(25, 50, 100),
                             labels = c("T=25","T=50","T=100"))
res_plot$n_units_f <- factor(res_plot$n_units)

p <- ggplot(res_plot, aes(x = n_units_f, y = rmse, colour = n_time_f, group = n_time_f)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  facet_wrap(~ param, scales = "free_y") +
  labs(
    title    = "RMSE of beta estimates (pvarife R package)",
    subtitle = sprintf("MC DGP: Tugan (2021), n_sds=%d, n_out=%d, n_in=%d", n_sds, n_out, n_in),
    x        = "Number of units (I)",
    y        = "RMSE",
    colour   = "Time periods"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

ggsave("scripts/mc_convergence.pdf", p, width = 8, height = 5)
cat("Plot saved to scripts/mc_convergence.pdf\n")
cat("Results saved to scripts/mc_results.rds\n\n")

# ---- Quick convergence diagnosis --------------------------------------------
cat("=== Convergence diagnosis ===\n")
cat("For each Theta_1 element, RMSE should decrease as I*T increases.\n\n")

for (pp in theta_params) {
  sub <- res_theta[res_theta$param == pp, ]
  sub <- sub[order(sub$n_units * sub$n_time), ]
  monotone <- all(diff(sub$rmse) <= 0.01)  # allow 0.01 tolerance for MC noise
  cat(sprintf("  %s: RMSE = [%s] — %s\n",
              pp,
              paste(round(sub$rmse, 4), collapse=", "),
              if (monotone) "CONVERGING (PASS)" else "check manually"))
}
cat("\nDone.\n")
