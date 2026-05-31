test_that("sim_pvarife: output dimensions are correct", {
  sim <- sim_pvarife(n_units = 20L, n_time = 15L, n_vars = 2L,
                     n_lags = 1L, n_factors = 1L, seed = 1L)
  expect_equal(dim(sim$y), c(20L, 15L, 2L))
  expect_equal(dim(sim$factors_true), c(15L, 1L))
  expect_equal(dim(sim$loadings_true), c(2L, 20L))
  expect_equal(length(sim$beta_true), 2L + 2L * 2L * 1L)  # K + K^2*n_lags
})

test_that("sim_pvarife: no NA in output", {
  sim <- sim_pvarife(n_units = 10L, n_time = 10L, n_vars = 2L,
                     n_lags = 1L, n_factors = 1L, seed = 2L)
  expect_false(anyNA(sim$y))
})

test_that("sim_pvarife: seed reproduces identical results", {
  s1 <- sim_pvarife(seed = 99L)
  s2 <- sim_pvarife(seed = 99L)
  expect_equal(s1$y, s2$y)
})

test_that("sim_pvarife: true theta_1 matches MATLAB DGP for 2-var case", {
  sim <- sim_pvarife(n_vars = 2L, n_lags = 1L, seed = 1L)
  expected <- matrix(c(0.65, 0.20, 0.30, 0.60), nrow = 2L)
  expect_equal(sim$theta_true[, , 1L], expected, tolerance = 1e-12)
})

test_that("sim_pvarife: long-run DGP returns correct dimensions and new slots", {
  sim_lr <- sim_pvarife(n_units = 20L, n_time = 15L, n_vars = 2L,
                        n_lags = 1L, n_factors = 1L,
                        identification = "long_run", seed = 30L)
  expect_equal(dim(sim_lr$y), c(20L, 15L, 2L))
  expect_equal(dim(sim_lr$a0_true), c(2L, 2L))
  expect_equal(sim_lr$identification, "long_run")
  expect_equal(sim_lr$diff_vars_suggested, 1L)
})

test_that("sim_pvarife: long-run a0_true produces lower-triangular long-run multiplier", {
  sim_lr <- sim_pvarife(n_vars = 2L, n_lags = 1L,
                        identification = "long_run", seed = 31L)
  theta1 <- sim_lr$theta_true[, , 1L]
  lr_inv <- solve(diag(2L) - theta1)
  c1     <- lr_inv %*% sim_lr$a0_true   # = chol(D)' — should be lower-triangular
  expect_equal(c1[1L, 2L], 0.0, tolerance = 1e-10)
})

test_that("sim_pvarife: short-run returns integer(0) for diff_vars_suggested", {
  sim_sr <- sim_pvarife(identification = "short_run", seed = 32L)
  expect_equal(sim_sr$diff_vars_suggested, integer(0))
})
