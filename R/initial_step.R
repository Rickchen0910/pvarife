# Obtain initial beta estimate for the pvarife iteration.
#
# Faithful translation of Initial_Step_in_Iteration.m from Tugan (2021).
# Uses OLS on the balanced subsample to get a starting value, then minimises
# .objective_fn() from three starting points and selects the best.
#
# @param y_stack Pooled (T*I) x 1 stacked outcome vector (NAs present).
# @param z_stack Pooled (T*I) x K+K^2*n_lags stacked regressors (NAs present).
# @param y_c     NT x 1 x n_units array.
# @param z_c     NT x n_cols_z x n_units array.
# @param y_arr   n_units x n_time x n_vars raw array (used for dimensions).
# @param n_factors Number of interactive fixed effects r.
#
# @return A list with:
#   beta_initial: initial coefficient vector (n_cols_z x 1)
#   w_c:          NT x 1 x n_units array of initial residuals (imputed at missing)
#   i_obs:        NT x n_units 0/1 indicator (1 = observed, 0 = missing/imputed)
#
# @noRd
.initial_step <- function(y_stack, z_stack, y_c, z_c, y_arr, n_factors,
                           balanced_init = TRUE) {
  n_units <- dim(y_arr)[1L]
  n_time  <- dim(y_arr)[2L]
  n_vars  <- dim(y_arr)[3L]

  if (balanced_init) {
    # Select units that have 10 balanced recent observations (last 11 periods).
    # MATLAB: find(sum(~any(isnan(NonStackedDVC(end-10:end-1,:,:)),2))==10)
    # Needed for unbalanced real data; for balanced MC data all units pass.
    window_end   <- n_time
    window_start <- max(1L, n_time - 10L)

    all_obs_in_window <- vapply(seq_len(n_units), function(ii) {
      slice <- y_arr[ii, window_start:(window_end - 1L), , drop = FALSE]
      all(!is.na(slice))
    }, logical(1L))

    units_bal <- which(all_obs_in_window)
    if (length(units_bal) == 0L) units_bal <- seq_len(n_units)  # fallback
  } else {
    # Skip balanced-subsample search: use all units for the initial OLS.
    # Appropriate for balanced panels (MC simulations) — faster and unambiguous.
    units_bal <- seq_len(n_units)
  }

  n_rows <- dim(y_c)[1L]

  # Extract balanced subsample y_c/z_c
  y_c_bal <- y_c[, , units_bal, drop = FALSE]
  z_c_bal <- z_c[, , units_bal, drop = FALSE]

  # Further trim to rows where no unit has NA (balanced panel across selected units)
  # MATLAB: data1(any(isnan(data1(:,1,:)),3),:,:) = []
  has_na_y <- apply(y_c_bal[, 1L, , drop = FALSE], 1L, function(x) any(is.na(x)))
  has_na_z <- apply(z_c_bal[, 2L, , drop = FALSE], 1L, function(x) any(is.na(x)))
  keep_rows <- !has_na_y & !has_na_z

  # Drop whole time periods, not individual rows: if any of the K rows of a
  # period is dropped, drop all K rows of that period. This keeps the row
  # count a multiple of K so the T x K reshape inside .objective_fn is
  # well-formed. (MATLAB trims row-wise and its reshape would error on a
  # partially-observed period; whole-period trimming handles that case.)
  period_of_row <- ceiling(seq_along(keep_rows) / n_vars)
  bad_periods   <- unique(period_of_row[!keep_rows])
  keep_rows     <- keep_rows & !(period_of_row %in% bad_periods)

  y_c_bal_trim <- y_c_bal[keep_rows, , , drop = FALSE]
  z_c_bal_trim <- z_c_bal[keep_rows, , , drop = FALSE]

  # NOTE on MATLAB bug: In the Economic Application's ObjectiveFunction.m,
  # Country_number = size(Z_c, 2) = n_cols_z (e.g., 20 for K=4, L=1) instead
  # of the correct size(Z_c, 3) = n_bal (true number of balanced countries).
  # This means MATLAB only evaluates the objective over the first n_cols_z
  # balanced countries — an unintentional underuse of the balanced sample.
  # This R implementation uses n_bal correctly via dim(z_c_bal_trim)[3].

  # OLS on full pooled (unbalanced) sample
  ols_res <- ols_nan(y_stack, z_stack)
  beta_ols <- ols_res$beta
  std_ols  <- if (!is.null(ols_res$std_err)) ols_res$std_err else
    rep(1.0, nrow(beta_ols))

  # Try 3 starting values: beta_ols + k * std_ols, k in {-2, 0, 2}
  k_vals <- c(-2.0, 0.0, 2.0)
  best_val  <- Inf
  best_beta <- beta_ols

  for (kk in k_vals) {
    beta0 <- beta_ols + kk * std_ols
    opt <- stats::optim(
      par     = as.numeric(beta0),
      fn      = function(b) {
        .objective_fn(matrix(b, ncol = 1L), n_factors,
                      y_c_bal_trim, z_c_bal_trim, y_arr)
      },
      method  = "BFGS",
      control = list(maxit = 1000L, reltol = 1e-8)
    )
    if (opt$value < best_val) {
      best_val  <- opt$value
      best_beta <- matrix(opt$par, ncol = 1L)
    }
  }

  beta_initial <- best_beta

  # Build initial W_c and observation indicator
  n_rows_yc <- dim(y_c)[1L]
  w_c <- array(0.0, dim = c(n_rows_yc, 1L, n_units))
  i_obs <- matrix(0L, nrow = n_rows_yc, ncol = n_units)

  for (ii in seq_len(n_units)) {
    yy <- y_c[, 1L, ii]
    zz <- z_c[, , ii]
    obs <- !is.na(yy)
    i_obs[obs, ii] <- 1L
    w_c[obs, 1L, ii] <- yy[obs] - zz[obs, , drop = FALSE] %*% beta_initial
    # Missing rows remain zero (as in MATLAB: W_c_initial(find(...),:,i) = zeros)
  }

  list(
    beta_initial = beta_initial,
    w_c          = w_c,
    i_obs        = i_obs
  )
}
