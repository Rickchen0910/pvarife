# Generalization invariants across dimensions (K, r, L).
#
# Rationale: layout bugs can be invisible at K = 2, r = 1, L = 1 (e.g., the
# loading-matrix reshape is only wrong for r >= 2; the beta_true layout was
# only wrong for L >= 2). These tests pin the algebraic identities that must
# hold for ANY dimension combination, on the non-trivial corners of the grid.

.check_invariants <- function(n_vars, n_factors, n_lags) {
  sim <- sim_pvarife(n_units = 15L, n_time = 12L, n_vars = n_vars,
                     n_lags = n_lags, n_factors = n_factors,
                     seed = 100L * n_vars + 10L * n_factors + n_lags)
  fit <- pvarife(sim$y, n_lags = n_lags, n_factors = n_factors,
                 n_out = 5L, n_in = 3L, balanced_init = FALSE)
  kk <- fit$n_vars; rr <- fit$n_factors

  # (a) factors_mat row (t,n) == e_n' (x) f_t'
  for (tt in c(1L, fit$n_time)) {
    for (nn in seq_len(kk)) {
      row <- fit$factors_mat[(tt - 1L) * kk + nn, ]
      expe <- numeric(kk * rr)
      expe[((nn - 1L) * rr + 1L):(nn * rr)] <- fit$ff[tt, ]
      expect_equal(row, expe, tolerance = 1e-12)
    }
  }

  # (b) loading layout: factors_mat %*% v reproduced by f_t' Lam[n, ] with
  #     Lam = t(matrix(v, r, K)) — the layout used inside asymptotic_var
  v   <- fit$loadings[, 1L]
  lam <- t(matrix(v, nrow = rr, ncol = kk))
  cc1 <- as.numeric(fit$factors_mat %*% v)
  cc2 <- numeric(length(cc1))
  for (tt in seq_len(fit$n_time)) {
    for (nn in seq_len(kk)) {
      cc2[(tt - 1L) * kk + nn] <- sum(fit$ff[tt, ] * lam[nn, ])
    }
  }
  expect_equal(cc1, cc2, tolerance = 1e-10)

  # (c) residual identity on observed rows: y = Z beta + F lambda + u
  obs <- which(fit$i_obs[, 1L] == 1L)
  lhs <- fit$y_c[obs, 1L, 1L]
  rhs <- as.numeric(matrix(fit$z_c[obs, , 1L], nrow = length(obs)) %*% fit$beta) +
         as.numeric(fit$factors_mat %*% fit$loadings[, 1L])[obs] +
         fit$u_c[obs, 1L, 1L]
  expect_equal(lhs, rhs, tolerance = 1e-8)

  # (d) asymptotic_var: dimensions, symmetry, PSD
  av <- asymptotic_var(fit)
  p  <- kk + kk^2 * fit$n_lags
  expect_equal(dim(av$variance), c(p, p))
  expect_lt(max(abs(av$variance - t(av$variance))), 1e-8)
  expect_gt(min(eigen(av$variance, symmetric = TRUE,
                      only.values = TRUE)$values), -1e-8)

  # (e) IRF identities: h=0 equals chol(Sigma)' column; long-run C(1) is
  #     lower-triangular
  ir <- compute_irf(fit, n_periods = 4L, shock = 1L)
  expect_equal(ir[, 1L], t(chol(fit$sigma))[, 1L], tolerance = 1e-10)

  beta_v <- as.numeric(fit$beta)
  alpha  <- t(matrix(beta_v[-seq_len(kk)], nrow = kk * fit$n_lags, ncol = kk))
  theta_sum <- matrix(0, kk, kk)
  for (ll in seq_len(fit$n_lags)) {
    theta_sum <- theta_sum + alpha[, ((ll - 1L) * kk + 1L):(ll * kk)]
  }
  lr_inv <- solve(diag(kk) - theta_sum)
  d_mat  <- lr_inv %*% fit$sigma %*% t(lr_inv)
  c1     <- lr_inv %*% ((diag(kk) - theta_sum) %*% t(chol(d_mat)))
  expect_lt(max(abs(c1[upper.tri(c1)])), 1e-8)

  invisible(NULL)
}

test_that("invariants hold for K=3, r=2, L=1", {
  .check_invariants(3L, 2L, 1L)
})

test_that("invariants hold for K=2, r=3, L=1", {
  .check_invariants(2L, 3L, 1L)
})

test_that("invariants hold for K=2, r=2, L=2", {
  .check_invariants(2L, 2L, 2L)
})

test_that("invariants hold for K=3, r=3, L=1", {
  .check_invariants(3L, 3L, 1L)
})
