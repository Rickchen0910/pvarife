# =============================================================================
# Replication: Tugan (2021) Economic Application
# Replicates EmpiricalApplication.m using the pvarife R package.
#
# MATLAB settings: kmax=1, rmax=2, out_number=50, in_number=10, nir=9
# Data: Liaqat (2019) panel, 1990-2017, 4 variables, all economies
#
# Comparison target: beta and IRFs from the full MATLAB run.
# Since we don't have MATLAB's exact beta output, we check:
#   (a) sigma positive definite
#   (b) IRF signs match Figure 1 qualitatively
#   (c) beta magnitudes are economically plausible
# =============================================================================

devtools::load_all(quiet = TRUE)
library(ggplot2)

cat("=== pvarife Empirical Replication ===\n\n")

# ---- 1. Load data -----------------------------------------------------------
csv_path <- file.path(
  "documentforrepl", "replication package", "Economic Application",
  "Liaqat (2019) Data.csv"
)
raw <- read.csv(csv_path, stringsAsFactors = FALSE)
cat("Raw data dimensions:", nrow(raw), "rows x", ncol(raw), "cols\n")
cat("Years:", range(raw$year), "\n")
cat("Countries (unique):", length(unique(raw$country)), "\n\n")

# ---- 2. Construct 4 variables (MATLAB EmpiricalApplication.m lines 103-118) -
# Columns mapping:
#   debt_cent_perc_imf = DebtAsShareofGDP
#   gdp_cons           = GDPinConstantDollars
#   rir                = RealInterestRate
#   fixcapfor_per      = GrossFixedFormationAsShareofGDP
#   population         = Population

# Filter to sample period 1990-2017
samp <- raw[raw$year >= 1990 & raw$year <= 2017, ]
years  <- 1990:2017
n_time_raw <- length(years)   # 28 time periods

countries <- sort(unique(samp$country))
n_units   <- length(countries)
cat("Countries in sample:", n_units, "\n")
cat("Raw time periods:   ", n_time_raw, "\n\n")

# Build T x I arrays for each variable
make_tmat <- function(df, years, countries, col) {
  mat <- matrix(NA_real_, nrow = length(years), ncol = length(countries))
  for (ii in seq_along(countries)) {
    sub <- df[df$country == countries[ii], ]
    idx <- match(years, sub$year)
    mat[, ii] <- ifelse(is.na(idx), NA_real_, sub[[col]][idx])
  }
  mat
}

DebtAsShareGDP  <- make_tmat(samp, years, countries, "debt_cent_perc_imf")
GDPConst        <- make_tmat(samp, years, countries, "gdp_cons")
RIR             <- make_tmat(samp, years, countries, "rir")
FixCapPerGDP    <- make_tmat(samp, years, countries, "fixcapfor_per")
Pop             <- make_tmat(samp, years, countries, "population")

# Construct 4 analysis variables (MATLAB lines 103-106)
DebtPC  <- DebtAsShareGDP * GDPConst / Pop    # DebtPerCapita
RIR_raw <- RIR                                  # RealInterestRate (levels)
FixCapPC <- FixCapPerGDP * GDPConst / Pop       # GrossFixedCapitalPerCapita
GDPPC    <- GDPConst / Pop                      # GDPPerCapita

# Set country-time to NA if ANY of the 4 variables is NA for that country-time
# (MATLAB lines 120-123)
any_na <- is.na(DebtPC) | is.na(RIR_raw) | is.na(FixCapPC) | is.na(GDPPC)
DebtPC[any_na]   <- NA
RIR_raw[any_na]  <- NA
FixCapPC[any_na] <- NA
GDPPC[any_na]    <- NA

# ---- 3. First differences (MATLAB lines 137-148) ----------------------------
# lag_lead_matrix(X, 1, 1) returns one lag
lag1 <- function(x) {
  rbind(NA_real_, x[-nrow(x), , drop = FALSE])
}

dDebt   <- 100 * (log(DebtPC)   - log(lag1(DebtPC)))      # % log change
dRIR    <- RIR_raw - lag1(RIR_raw)                          # level change
dFixCap <- 100 * (log(FixCapPC) - log(lag1(FixCapPC)))     # % log change
dGDP    <- 100 * (log(GDPPC)    - log(lag1(GDPPC)))        # % log change

# After first-differencing: T=28 rows, first row is NA for all countries
# Effective: 27 periods (1991-2017)

# ---- 4. Build y[i, t, k] array (I x T x K) ---------------------------------
# Variables: k=1 Debt, k=2 RIR, k=3 FixedCapital, k=4 GDP
# Dimensions: I=n_units, T=28, K=4
n_vars <- 4L
y_arr  <- array(NA_real_, dim = c(n_units, n_time_raw, n_vars))
y_arr[, , 1L] <- t(dDebt)
y_arr[, , 2L] <- t(dRIR)
y_arr[, , 3L] <- t(dFixCap)
y_arr[, , 4L] <- t(dGDP)

# Report how many countries have any valid data
n_obs_per_unit <- apply(y_arr[, , 1L], 1L, function(x) sum(!is.na(x)))
cat("Obs per country (summary):\n")
print(summary(n_obs_per_unit))
cat("\nCountries with >= 10 obs:", sum(n_obs_per_unit >= 10L), "\n\n")

# ---- 5. Estimate (MATLAB: kmax=1, r=2, out_number=50, in_number=10) --------
cat("Running pvarife estimation...\n")
cat("  n_lags=1, n_factors=2, n_out=50, n_in=10\n")
cat("  This matches the MATLAB rmax=2 (Panel B of Figure 1)\n\n")

t_start <- proc.time()
fit <- pvarife(
  y         = y_arr,
  n_lags    = 1L,
  n_factors = 2L,
  n_out     = 50L,
  n_in      = 10L
)
t_elapsed <- proc.time() - t_start
cat(sprintf("Estimation completed in %.1f seconds.\n\n", t_elapsed["elapsed"]))

# ---- 6. Report beta ---------------------------------------------------------
cat("=== Estimated beta ===\n")
beta <- coef(fit)
n_vars2 <- fit$n_vars

# Reshape into intercepts + Theta_1 matrix
intercepts <- beta[seq_len(n_vars2)]
theta1_vec <- beta[(n_vars2 + 1L):length(beta)]
theta1     <- t(matrix(theta1_vec, nrow = n_vars2 * fit$n_lags, ncol = n_vars2))

cat("Intercepts:\n")
cat("  Debt =", round(intercepts[1], 4),
    "  RIR =",  round(intercepts[2], 4),
    "  FixCap =", round(intercepts[3], 4),
    "  GDP =",  round(intercepts[4], 4), "\n\n")

cat("Theta_1 (VAR coefficient matrix):\n")
rownames(theta1) <- c("Debt", "RIR", "FixCap", "GDP")
colnames(theta1) <- paste0(c("L.Debt", "L.RIR", "L.FixCap", "L.GDP"))
print(round(theta1, 4))
cat("\n")

cat("Reduced-form covariance (Sigma):\n")
rownames(fit$sigma) <- colnames(fit$sigma) <- c("Debt","RIR","FixCap","GDP")
print(round(fit$sigma, 4))
cat("\nSigma eigenvalues (all > 0 required):",
    round(eigen(fit$sigma, only.values = TRUE)$values, 6), "\n\n")

# ---- 7. Compute IRFs (MATLAB: nir=9, DifferencedVariables=[]) --------------
# IMPORTANT: MATLAB's Figure 1 "point estimate" is the MEDIAN of 500
# simulation draws from N(beta-bias, V) — NOT the direct IRF at beta_hat.
# Use irf_bands()$irf to match MATLAB's approach (bias-corrected median).
# compute_irf(fit) gives the uncorrected direct estimate for reference.

cat("Computing IRFs — two versions:\n")
cat("  (a) Direct (compute_irf): uncorrected, matches nothing directly in Figure 1\n")
cat("  (b) Bias-corrected median (irf_bands): matches MATLAB ConfidenceBandforIRs\n\n")
cat("Running irf_bands (n_draw=500, seed=42) — this may take a few minutes...\n")

set.seed(42)
bands <- irf_bands(fit, n_periods = 9L, shock = 1L, n_draw = 500L, seed = 42L)

# MATLAB normalises by IRs(1,1) from the SIMULATION MEDIAN
norm_factor <- bands$irf[1L, 1L]
ir_norm <- bands$irf / norm_factor
lower_norm <- bands$lower / norm_factor
upper_norm <- bands$upper / norm_factor

cat("\n=== Bias-corrected Median IRF (normalised, matches MATLAB Figure 1 approach) ===\n")
ir_df <- as.data.frame(t(ir_norm))
names(ir_df) <- c("Debt", "RIR", "FixCap", "GDP")
ir_df$horizon <- 0L:(nrow(ir_df) - 1L)
print(round(ir_df, 4))

# ---- 8. Qualitative check against Figure 1 (Panel B, rmax=2) ---------------
cat("\n=== Qualitative Check vs. Paper Figure 1 ===\n")
cat("Expected from Figure 1 (shock = ΔDebt, h=0 to 8, normalised):\n")
cat("  Debt:   starts at 1.0, then declines (mean-reverting)\n")
cat("  GDP:    should show a negative or near-zero response (crowding out)\n")
cat("  FixCap: expected negative response (investment crowded out)\n\n")

cat("Observed:\n")
for (vv in c("Debt", "RIR", "FixCap", "GDP")) {
  vals <- round(ir_df[[vv]][1:5], 3)
  sign_pat <- ifelse(vals >= 0, "+", "-")
  cat(sprintf("  %-8s h=0..4: %s  (%.3f %.3f %.3f %.3f %.3f)\n",
              vv, paste(sign_pat, collapse=""), vals[1], vals[2], vals[3], vals[4], vals[5]))
}

# ---- 9. Plot ----------------------------------------------------------------
cat("\nSaving IRF plot to scripts/replication_irf.pdf ...\n")
ir_long <- tidyr::pivot_longer(ir_df, -horizon,
                                names_to = "Variable", values_to = "IRF")
ir_long$Variable <- factor(ir_long$Variable, levels = c("Debt","RIR","FixCap","GDP"))

p <- ggplot(ir_long, aes(x = horizon, y = IRF)) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.3) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2, shape = 18) +
  facet_wrap(~Variable, scales = "free_y", ncol = 2) +
  labs(
    title    = "IRF to Debt shock — pvarife R package",
    subtitle = "n_lags=1, n_factors=2 (Panel B of Figure 1, Tugan 2021)",
    x        = "Horizon (years)", y = "Normalised response"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))

ggsave("scripts/replication_irf.pdf", p, width = 8, height = 5)
cat("Plot saved.\n\n")

cat("=== Summary ===\n")
cat("Units (I):", fit$n_units, "\n")
cat("Time (T):", fit$n_time, "\n")
cat("Variables (K):", fit$n_vars, "\n")
cat("Factors (r):", fit$n_factors, "\n")
cat("n_time_i range:", range(fit$n_time_i[fit$n_time_i > 0]), "\n")
cat("Effective obs (sum TC):", sum(fit$n_time_i), "\n\n")
cat("Done. Compare with MATLAB EmpiricalApplication.m output manually.\n")
