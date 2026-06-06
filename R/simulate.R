#' Simulate panel VAR data with interactive fixed effects
#'
#' Generates synthetic panel data from the DGP of Tugan (2021), Section S10,
#' with a single common factor following an AR(1) process and factor loadings
#' drawn from \eqn{N(1, 1)}.
#'
#' @details
#' The data generating process is:
#' \deqn{y_{i,t} = c + \Theta_1 y_{i,t-1} + f_t \lambda_i + A_0 \varepsilon_{i,t},}
#' where
#' \deqn{\Theta_1 = \begin{pmatrix} 0.65 & 0.30 \\ 0.20 & 0.60 \end{pmatrix},
#'       \quad \Sigma_e = \begin{pmatrix} 1 & 0.5 \\ 0.5 & 1 \end{pmatrix},}
#' \deqn{f_t = 0.5 f_{t-1} + \eta_t, \; \eta_t \sim N(0, 0.5),}
#' \deqn{\lambda_{k,i} \sim N(1, 1) \;\text{(independently)},}
#' and \eqn{A_0 = \mathrm{chol}(\Sigma_e)^\top} (Cholesky identification). A
#' burn-in of 1000 periods is discarded.
#'
#' For \code{n_vars = 2} and \code{n_lags = 1}, the parameters match the MATLAB
#' \code{GeneratingSynteticDataSets.m} (Short-Run identification variant).
#'
#' @param n_units Positive integer. Number of cross-sectional units \eqn{I}
#'   (default 50).
#' @param n_time Positive integer. Number of time periods \eqn{T} after
#'   burn-in (default 30).
#' @param n_vars Positive integer. Number of VAR variables \eqn{K} (default 2).
#' @param n_lags Positive integer. VAR lag order (default 1; higher lags use
#'   zero coefficient matrices).
#' @param n_factors Positive integer. Number of common factors (default 1;
#'   each factor has its own AR(1) process and independent loadings).
#' @param identification Character. \code{"short_run"} (default) uses
#'   \eqn{A_0 = \mathrm{chol}(\Sigma_e)^\top}; \code{"long_run"} uses the
#'   Blanchard-Quah impact matrix \eqn{A_0 = (I - \Theta_1)\,\mathrm{chol}(D)^\top}
#'   where \eqn{D = (I-\Theta_1)^{-1}\Sigma_e(I-\Theta_1)^{-\top}}. The
#'   long-run multiplier \eqn{C(1) = (I-\Theta_1)^{-1}A_0 = \mathrm{chol}(D)^\top}
#'   is lower-triangular. Matches the MATLAB
#'   \code{GeneratingSynteticDataSets.m} with
#'   \code{Identification='WithLongRunRestrictions'}.
#' @param seed Optional integer seed for reproducibility (default 42).
#'
#' @return A list with:
#'   \describe{
#'     \item{y}{Array of dimension \eqn{I \times T \times K} (unit \eqn{\times}
#'       time \eqn{\times} variable). No missing values.}
#'     \item{beta_true}{True coefficient vector, same format as \code{fit$beta}.}
#'     \item{theta_true}{True VAR coefficient array \eqn{K \times K \times n\_lags}.}
#'     \item{sigma_true}{True reduced-form covariance matrix \eqn{K \times K}.}
#'     \item{a0_true}{True structural impact matrix \eqn{K \times K}.}
#'     \item{factors_true}{\eqn{T \times n\_factors} matrix of true factors.}
#'     \item{loadings_true}{\eqn{K \times n\_units} matrix of true loadings.}
#'     \item{identification}{The identification scheme used (\code{"short_run"}
#'       or \code{"long_run"}).}
#'     \item{diff_vars_suggested}{Integer vector: variables that should be
#'       cumulated in \code{compute_irf()}. \code{integer(0)} for short-run;
#'       \code{1L} for long-run (to display cumulative responses, matching the
#'       MATLAB MC code's \code{DifferencedVariables = [1]}).}
#'   }
#'
#' @examples
#' sim <- sim_pvarife(n_units = 30, n_time = 20, n_vars = 2,
#'                    n_lags = 1, n_factors = 1, seed = 1)
#' dim(sim$y)        # 30 x 20 x 2
#' sim$beta_true
#'
#' @export
sim_pvarife <- function(n_units = 50L, n_time = 30L, n_vars = 2L,
                         n_lags = 1L, n_factors = 1L,
                         identification = c("short_run", "long_run"),
                         seed = 42L) {
  identification <- match.arg(identification)
  n_units   <- as.integer(n_units)
  n_time    <- as.integer(n_time)
  n_vars    <- as.integer(n_vars)
  n_lags    <- as.integer(n_lags)
  n_factors <- as.integer(n_factors)
  if (!is.null(seed)) set.seed(seed)

  burn_in <- 1000L
  n_total <- n_time + burn_in

  # --- True parameters (matching MATLAB GeneratingSynteticDataSets.m) ---
  if (n_vars == 2L && n_lags == 1L) {
    theta_arr <- array(0.0, dim = c(n_vars, n_vars, n_lags))
    theta_arr[, , 1L] <- matrix(c(0.65, 0.20, 0.30, 0.60), nrow = 2L)
    sigma_true <- matrix(c(1.0, 0.5, 0.5, 1.0), nrow = 2L)
  } else {
    # Generic stable VAR: diagonal AR with coefficient 0.5, unit variance
    theta_arr <- array(0.0, dim = c(n_vars, n_vars, n_lags))
    theta_arr[, , 1L] <- diag(0.5, n_vars)
    sigma_true <- diag(n_vars)
  }
  # Impact matrix: depends on identification scheme
  if (identification == "short_run") {
    # Short-run: Cholesky of Sigma — impact matrix itself is lower-triangular
    # MATLAB: A_0 = chol(Sigma_e)'  (WithShortRunRestrictions)
    a0 <- t(chol(sigma_true))
  } else {
    # Long-run: Blanchard-Quah type — long-run multiplier is lower-triangular
    # MATLAB: D = inv(I-Theta_1)*Sigma_e*inv(I-Theta_1)'; A_0 = (I-Theta_1)*chol(D)'
    theta1  <- theta_arr[, , 1L]
    lr_inv  <- solve(diag(n_vars) - theta1)
    d_mat   <- lr_inv %*% sigma_true %*% t(lr_inv)
    a0      <- (diag(n_vars) - theta1) %*% t(chol(d_mat))
  }

  # Stationary mean: (I - Theta_1)^{-1} * c_vec; use c_vec = rep(1, K)
  c_vec <- rep(1.0, n_vars)
  if (n_lags == 1L) {
    y_mean <- solve(diag(n_vars) - theta_arr[, , 1L]) %*% c_vec
  } else {
    y_mean <- rep(0.0, n_vars)
  }

  # --- Common factors: n_factors independent AR(1) processes ---
  # f[t] = 0.5 * f[t-1] + eta[t], eta ~ N(0, 0.5)
  f_all <- matrix(0.0, nrow = n_total + 1L, ncol = n_factors)
  for (rr in seq_len(n_factors)) {
    eta <- stats::rnorm(n_total, mean = 0.0, sd = sqrt(0.5))
    for (tt in seq(2L, n_total + 1L)) {
      f_all[tt, rr] <- 0.5 * f_all[tt - 1L, rr] + eta[tt - 1L]
    }
  }
  f_sim <- f_all[2L:(n_total + 1L), , drop = FALSE]  # n_total x n_factors

  # --- Loadings: lambda[k, i, r] ~ N(1, 1), one per variable-unit-factor ---
  lam_true <- array(stats::rnorm(n_vars * n_units * n_factors, mean = 1.0, sd = 1.0),
                    dim = c(n_vars, n_units, n_factors))

  # --- Simulate VAR ---
  # y[k, t, i] in MATLAB notation; we produce y_big[k, n_total, i]
  y_big <- array(0.0, dim = c(n_vars, n_total, n_units))
  for (ii in seq_len(n_units)) {
    y_big[, 1L, ii] <- y_mean
  }

  for (tt in seq(2L, n_total)) {
    # Common component: K x 1, same for all units except loading
    f_t <- matrix(f_sim[tt, ], ncol = 1L)  # n_factors x 1

    for (ii in seq_len(n_units)) {
      # loading for unit ii: K x n_factors
      lam_i <- lam_true[, ii, , drop = FALSE]
      lam_i <- matrix(lam_i, nrow = n_vars, ncol = n_factors)

      # VAR contribution
      y_lag <- matrix(0.0, nrow = n_vars, ncol = 1L)
      for (ll in seq_len(n_lags)) {
        if (tt - ll >= 1L)
          y_lag <- y_lag + theta_arr[, , ll] %*% y_big[, tt - ll, ii]
      }

      eps_i <- matrix(stats::rnorm(n_vars), ncol = 1L)
      y_big[, tt, ii] <- c_vec + y_lag + lam_i %*% f_t + a0 %*% eps_i
    }
  }

  # Discard burn-in; permute to [i, t, k] = [n_units, n_time, n_vars]
  y_stable <- y_big[, (burn_in + 1L):n_total, , drop = FALSE]
  # y_stable is [k, t, i]; permute to [i, t, k]
  y_out <- aperm(y_stable, c(3L, 2L, 1L))  # I x T x K

  # True factors (post burn-in): T x n_factors
  f_true <- f_sim[(burn_in + 1L):n_total, , drop = FALSE]

  # True loadings summary: K x I (sum over factors for display)
  lam_summary <- apply(lam_true, c(1L, 2L), sum)  # K x I

  # Build beta_true in the SAME layout as fit$beta (see .build_yz):
  # intercepts first, then for each equation k, the lag coefficients nested as
  # (lag 1 vars 1..K, lag 2 vars 1..K, ...). i.e. equation-major, lag-nested.
  #   beta[K + (k-1)*K*L + (l-1)*K + j] = Theta_l[k, j]
  # For n_lags = 1 this reduces to the row-major vec of Theta_1.
  theta_vec <- numeric(0L)
  for (kk in seq_len(n_vars)) {        # equation
    for (ll in seq_len(n_lags)) {      # lag
      for (jj in seq_len(n_vars)) {    # regressor variable
        theta_vec <- c(theta_vec, theta_arr[kk, jj, ll])
      }
    }
  }
  beta_true <- c(c_vec, theta_vec)

  diff_vars_suggested <- if (identification == "long_run") 1L else integer(0L)

  list(
    y                  = y_out,
    beta_true          = beta_true,
    theta_true         = theta_arr,
    sigma_true         = sigma_true,
    a0_true            = a0,
    factors_true       = f_true,
    loadings_true      = lam_summary,
    identification     = identification,
    diff_vars_suggested = diff_vars_suggested
  )
}
