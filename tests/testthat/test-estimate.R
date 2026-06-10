test_that("pvarife: smoke test on small simulated panel", {
  sim <- sim_pvarife(n_units = 15L, n_time = 12L, n_vars = 2L,
                     n_lags = 1L, n_factors = 1L, seed = 3L)
  # Very few iterations for speed; just check it runs without error
  fit <- pvarife(sim$y, n_lags = 1L, n_factors = 1L, n_out = 3L, n_in = 3L)
  expect_s3_class(fit, "pvarife_result")
})

test_that("pvarife: output structure is complete", {
  sim <- sim_pvarife(n_units = 15L, n_time = 12L, n_vars = 2L,
                     n_lags = 1L, n_factors = 1L, seed = 4L)
  fit <- pvarife(sim$y, n_lags = 1L, n_factors = 1L, n_out = 3L, n_in = 3L)

  expect_true(all(c("beta", "ff", "factors_mat", "loadings", "sigma",
                    "u_c", "y_c", "z_c", "i_obs",
                    "n_lags", "n_factors", "n_vars", "n_units", "n_time")
                  %in% names(fit)))
  # beta has K + K^2*n_lags elements
  expect_equal(length(fit$beta), 2L + 2L^2 * 1L)
  # ff: T x r
  expect_equal(dim(fit$ff), c(12L, 1L))
  # sigma: K x K, symmetric
  expect_equal(nrow(fit$sigma), 2L)
  expect_equal(fit$sigma, t(fit$sigma), tolerance = 1e-12)
})

test_that("pvarife: sigma is positive definite", {
  sim <- sim_pvarife(n_units = 20L, n_time = 15L, n_vars = 2L,
                     n_lags = 1L, n_factors = 1L, seed = 5L)
  fit <- pvarife(sim$y, n_lags = 1L, n_factors = 1L, n_out = 5L, n_in = 5L)
  ev  <- eigen(fit$sigma, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(ev > 0))
})

test_that("pvarife: coef() returns named vector of correct length", {
  sim <- sim_pvarife(n_units = 10L, n_time = 10L, n_vars = 2L,
                     n_lags = 1L, n_factors = 1L, seed = 6L)
  fit <- pvarife(sim$y, n_lags = 1L, n_factors = 1L, n_out = 3L, n_in = 3L)
  b   <- coef(fit)
  expect_equal(length(b), 2L + 4L)
  expect_false(is.null(names(b)))
})

test_that("pvarife: balanced_init=FALSE gives same result on balanced data", {
  sim  <- sim_pvarife(n_units = 15L, n_time = 12L, n_vars = 2L,
                      n_lags = 1L, n_factors = 1L, seed = 20L)
  # With balanced data, both options should select all units and give same result
  fit_t <- pvarife(sim$y, n_lags = 1L, n_factors = 1L, n_out = 5L, n_in = 5L,
                   balanced_init = TRUE)
  fit_f <- pvarife(sim$y, n_lags = 1L, n_factors = 1L, n_out = 5L, n_in = 5L,
                   balanced_init = FALSE)
  # Point estimates should be very close (same algorithm, same data)
  expect_equal(as.numeric(fit_t$beta), as.numeric(fit_f$beta), tolerance = 1e-4)
})

test_that("partial-variable missingness: Sigma uses complete periods, no warnings", {
  # Unit 3 missing var 1 only at t = 7 (with n_lags = 1 this invalidates the
  # (7, var1) row and all rows at t = 8). Observed row count becomes odd;
  # a naive reshape would silently recycle. Regression test for that bug.
  sim <- sim_pvarife(n_units = 15L, n_time = 20L, seed = 3L)
  y <- sim$y
  y[3L, 7L, 1L] <- NA
  expect_no_warning(
    fit <- pvarife(y, n_lags = 1L, n_factors = 1L, n_out = 5L, n_in = 3L)
  )
  # Complete periods for unit 3: t = 2..20 minus t = 7 (partial) and t = 8
  expect_identical(fit$n_time_i[3L], 17L)
  expect_true(all(is.finite(fit$sigma)))
  # Full inference chain runs on partially-missing data
  av <- asymptotic_var(fit)
  expect_true(all(is.finite(av$bias)))
})
