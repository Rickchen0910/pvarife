# Objective function for the initial step of the pvarife estimator.
#
# Faithful translation of ObjectiveFuncton.m (typo preserved in comments) from
# Tugan (2021). Computes the sum of the (n_obs - n_factors) smallest eigenvalues
# of (1 / (T * K * I)) * sum_i v_i v_i', where v_i is the T x K matrix of
# residuals for unit i.
#
# Paper reference: Eq. (2.44).
#
# @noRd
.objective_fn <- function(beta, n_factors, y_c, z_c, y_arr) {
  n_vars <- dim(y_arr)[3L]

  # y_c and z_c here are the *balanced-trimmed subsample* arrays passed from
  # .initial_step(). Their row count (n_rows_bal) equals n_vars * T_eff where
  # T_eff = T - n_lags is the effective time after lag trimming.
  # We infer T_eff from the data to avoid confusion with the full panel n_time.
  n_bal     <- dim(z_c)[3L]
  n_rows_bal <- nrow(y_c)           # n_vars * T_eff
  n_time_eff <- n_rows_bal %/% n_vars  # effective periods in balanced subsample

  vvprime <- matrix(0.0, nrow = n_time_eff, ncol = n_time_eff)

  for (ii in seq_len(n_bal)) {
    zz    <- matrix(z_c[, , ii], nrow = n_rows_bal, ncol = ncol(z_c[, , ii]))
    w_vec <- as.numeric(y_c[, 1L, ii]) -
             as.numeric(zz %*% matrix(beta, ncol = 1L))
    # Reshape (n_vars * T_eff) vector -> T_eff x n_vars matrix
    # MATLAB: reshape(W_c', n_vars, T_eff)'
    v_mat <- t(matrix(w_vec, nrow = n_vars, ncol = n_time_eff))
    vvprime <- vvprime + tcrossprod(v_mat)
  }

  scale <- as.numeric(n_time_eff * n_vars * n_bal)
  vvprime_scaled <- vvprime / scale
  ev <- sort(eigen(vvprime_scaled, symmetric = TRUE, only.values = TRUE)$values,
             decreasing = TRUE)
  sum(ev[seq(n_factors + 1L, length(ev))])
}
