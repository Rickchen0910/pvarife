test_that("asymptotic_var: returns list with bias and variance", {
  sim <- sim_pvarife(n_units = 20L, n_time = 15L, n_vars = 2L,
                     n_lags = 1L, n_factors = 1L, seed = 7L)
  fit <- pvarife(sim$y, n_lags = 1L, n_factors = 1L, n_out = 5L, n_in = 5L)
  avar <- asymptotic_var(fit)

  expect_true(all(c("bias", "variance") %in% names(avar)))
  expect_equal(length(avar$bias), length(fit$beta))
  expect_equal(dim(avar$variance), c(length(fit$beta), length(fit$beta)))
})

test_that("asymptotic_var: variance is symmetric and PSD", {
  sim <- sim_pvarife(n_units = 20L, n_time = 15L, n_vars = 2L,
                     n_lags = 1L, n_factors = 1L, seed = 8L)
  fit <- pvarife(sim$y, n_lags = 1L, n_factors = 1L, n_out = 5L, n_in = 5L)
  avar <- asymptotic_var(fit)

  # Symmetric
  expect_equal(avar$variance, t(avar$variance), tolerance = 1e-10)
  # Positive semi-definite (eigenvalues >= 0)
  ev <- eigen(avar$variance, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(ev >= -1e-10))
})

test_that("extract_factors: returns correct dimensions", {
  sim <- sim_pvarife(n_units = 15L, n_time = 12L, n_vars = 2L,
                     n_lags = 1L, n_factors = 1L, seed = 9L)
  fit <- pvarife(sim$y, n_lags = 1L, n_factors = 1L, n_out = 3L, n_in = 3L)
  ef  <- extract_factors(fit$beta, fit, n_in = 3L)

  expect_equal(dim(ef$ff), c(12L, 1L))   # T x r
  expect_equal(ncol(ef$loadings), 15L)   # r*K x n_units -> ncol = n_units
})

test_that("loading reshape uses variable-major layout (critical for r >= 2)", {
  # The loading vector follows factors_mat's column layout: variable-major
  # blocks of r. The K x r loading matrix must be t(matrix(v, r, K)).
  # Regression test for a bug where matrix(v, K, r) scrambled r >= 2 layouts.
  sim <- sim_pvarife(n_units = 20L, n_time = 15L, n_vars = 2L,
                     n_lags = 1L, n_factors = 2L, seed = 7L)
  fit <- pvarife(sim$y, n_lags = 1L, n_factors = 2L, n_out = 3L, n_in = 3L,
                 balanced_init = FALSE)
  v       <- fit$loadings[, 1L]
  cc_true <- as.numeric(fit$factors_mat %*% v)
  lam     <- t(matrix(v, nrow = fit$n_factors, ncol = fit$n_vars))  # K x r
  cc_alt  <- numeric(length(cc_true))
  for (tt in seq_len(fit$n_time)) {
    for (nn in seq_len(fit$n_vars)) {
      cc_alt[(tt - 1L) * fit$n_vars + nn] <- sum(fit$ff[tt, ] * lam[nn, ])
    }
  }
  expect_equal(cc_alt, cc_true, tolerance = 1e-10)

  # asymptotic_var must run and give a PSD variance with r = 2
  av <- asymptotic_var(fit)
  ev <- eigen(av$variance, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(is.finite(av$bias)))
  expect_true(all(ev > -1e-8))
})
