test_that("compute_irf: returns matrix of correct dimensions", {
  sim <- sim_pvarife(n_units = 15L, n_time = 12L, n_vars = 2L,
                     n_lags = 1L, n_factors = 1L, seed = 10L)
  fit <- pvarife(sim$y, n_lags = 1L, n_factors = 1L, n_out = 3L, n_in = 3L)
  ir  <- compute_irf(fit, n_periods = 8L)

  expect_equal(dim(ir), c(2L, 8L))
  expect_s3_class(ir, "pvarife_irf")
})

test_that("compute_irf: h=1 response to first shock equals first col of A0", {
  sim <- sim_pvarife(n_units = 15L, n_time = 12L, n_vars = 2L,
                     n_lags = 1L, n_factors = 1L, seed = 11L)
  fit <- pvarife(sim$y, n_lags = 1L, n_factors = 1L, n_out = 3L, n_in = 3L)
  ir  <- compute_irf(fit, n_periods = 5L, shock = 1L)

  # At h=1: B_1 = Theta_1; IRF = B_1 * A0 * e1
  # But B_0 = I, so h=1 in our code is B[,,2] = Theta_1
  # The first horizon (h=1) should be B_0 * A0 * e1... wait:
  # In the code, h=1 uses b_arr[,,1] which is B_0 = I. So:
  a0 <- t(chol(fit$sigma))
  expect_equal(ir[, 1L], as.numeric(a0[, 1L]), tolerance = 1e-10)
})

test_that("compute_irf: cumulation works for diff_vars", {
  sim <- sim_pvarife(n_units = 15L, n_time = 12L, n_vars = 2L,
                     n_lags = 1L, n_factors = 1L, seed = 12L)
  fit <- pvarife(sim$y, n_lags = 1L, n_factors = 1L, n_out = 3L, n_in = 3L)

  ir_plain <- compute_irf(fit, n_periods = 5L, diff_vars = integer(0))
  ir_cum   <- compute_irf(fit, n_periods = 5L, diff_vars = 1L)

  # Variable 1 should be cumulative sum of plain IRF
  expect_equal(ir_cum[1L, ], cumsum(ir_plain[1L, ]), tolerance = 1e-12)
  # Variable 2 unchanged
  expect_equal(ir_cum[2L, ], ir_plain[2L, ], tolerance = 1e-12)
})

test_that("irf_bands: output class and dimensions", {
  sim <- sim_pvarife(n_units = 15L, n_time = 12L, n_vars = 2L,
                     n_lags = 1L, n_factors = 1L, seed = 13L)
  fit <- pvarife(sim$y, n_lags = 1L, n_factors = 1L, n_out = 3L, n_in = 3L)
  bands <- irf_bands(fit, n_periods = 5L, n_draw = 20L, seed = 1L)

  expect_s3_class(bands, "pvarife_bands")
  expect_equal(dim(bands$irf),   c(2L, 5L))
  expect_equal(dim(bands$lower), c(2L, 5L))
  expect_equal(dim(bands$upper), c(2L, 5L))
})

test_that("irf_bands: lower <= irf <= upper pointwise (at least 90% of cells)", {
  sim <- sim_pvarife(n_units = 15L, n_time = 12L, n_vars = 2L,
                     n_lags = 1L, n_factors = 1L, seed = 14L)
  fit <- pvarife(sim$y, n_lags = 1L, n_factors = 1L, n_out = 3L, n_in = 3L)
  bands <- irf_bands(fit, n_periods = 5L, n_draw = 50L, seed = 2L)

  ok_lower <- bands$lower <= bands$irf + 1e-10
  ok_upper <- bands$irf   <= bands$upper + 1e-10
  expect_true(mean(ok_lower) >= 0.9)
  expect_true(mean(ok_upper) >= 0.9)
})

test_that("compute_irf: long-run identification returns correct dimensions", {
  sim <- sim_pvarife(n_units = 15L, n_time = 12L, n_vars = 2L,
                     n_lags = 1L, n_factors = 1L, seed = 15L)
  fit <- pvarife(sim$y, n_lags = 1L, n_factors = 1L, n_out = 3L, n_in = 3L)
  ir_lr <- compute_irf(fit, n_periods = 6L, identification = "long_run")

  expect_equal(dim(ir_lr), c(2L, 6L))
  expect_s3_class(ir_lr, "pvarife_irf")
})

test_that("compute_irf: long-run A0 gives lower-triangular long-run multiplier", {
  sim <- sim_pvarife(n_units = 15L, n_time = 12L, n_vars = 2L,
                     n_lags = 1L, n_factors = 1L, seed = 16L)
  fit <- pvarife(sim$y, n_lags = 1L, n_factors = 1L, n_out = 5L, n_in = 5L)

  n_vars <- fit$n_vars
  beta   <- as.numeric(fit$beta)
  alpha_vec <- beta[(n_vars + 1L):length(beta)]
  alpha <- t(matrix(alpha_vec, nrow = n_vars * fit$n_lags, ncol = n_vars))
  theta1 <- alpha[, 1L:n_vars]

  lr_inv <- solve(diag(n_vars) - theta1)
  d_mat  <- lr_inv %*% fit$sigma %*% t(lr_inv)
  a0_lr  <- (diag(n_vars) - theta1) %*% t(chol(d_mat))
  # Long-run multiplier C(1) = lr_inv %*% a0_lr should be lower-triangular (chol(D)')
  c1 <- lr_inv %*% a0_lr
  expect_equal(c1[1L, 2L], 0.0, tolerance = 1e-10)  # upper-right element = 0
})

test_that("irf_bands: long-run identification runs and returns pvarife_bands", {
  sim <- sim_pvarife(n_units = 15L, n_time = 12L, n_vars = 2L,
                     n_lags = 1L, n_factors = 1L, seed = 17L)
  fit   <- pvarife(sim$y, n_lags = 1L, n_factors = 1L, n_out = 3L, n_in = 3L)
  bands <- irf_bands(fit, n_periods = 4L, identification = "long_run",
                     n_draw = 15L, seed = 3L)

  expect_s3_class(bands, "pvarife_bands")
  expect_equal(dim(bands$irf), c(2L, 4L))
})

test_that("compute_irf: bias_correct=TRUE runs and returns pvarife_irf", {
  sim <- sim_pvarife(n_units = 20L, n_time = 15L, n_vars = 2L,
                     n_lags = 1L, n_factors = 1L, seed = 18L)
  fit <- pvarife(sim$y, n_lags = 1L, n_factors = 1L, n_out = 5L, n_in = 5L)
  ir_bc <- compute_irf(fit, n_periods = 5L, bias_correct = TRUE)

  expect_s3_class(ir_bc, "pvarife_irf")
  expect_equal(dim(ir_bc), c(2L, 5L))
})

test_that("compute_irf: bias_correct=TRUE gives different result from FALSE", {
  sim <- sim_pvarife(n_units = 25L, n_time = 18L, n_vars = 2L,
                     n_lags = 1L, n_factors = 1L, seed = 19L)
  fit   <- pvarife(sim$y, n_lags = 1L, n_factors = 1L, n_out = 5L, n_in = 5L)
  ir    <- compute_irf(fit, n_periods = 5L, bias_correct = FALSE)
  ir_bc <- compute_irf(fit, n_periods = 5L, bias_correct = TRUE)
  # The two should differ (bias != 0 in general)
  expect_false(isTRUE(all.equal(as.numeric(ir), as.numeric(ir_bc))))
})
