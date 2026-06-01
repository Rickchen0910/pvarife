#' Parametric confidence bands for impulse response functions
#'
#' Constructs confidence bands for structural impulse responses by drawing from
#' the joint asymptotic distribution of the estimator and the common component
#' (Theorem 2.3 of Tugan 2021). This is a parametric simulation, not a
#' residual bootstrap.
#'
#' @details
#' For each of the \code{n_draw} repetitions, the function:
#' \enumerate{
#'   \item Draws \eqn{\beta^{(b)} \sim N(\hat\beta - \hat b, \hat V)} where
#'     \eqn{\hat b} and \eqn{\hat V} are the estimated bias and variance from
#'     \code{\link{asymptotic_var}}.
#'   \item Draws the common component \eqn{\tilde C_{i,t,k}} from a normal with
#'     mean \eqn{\hat F_t \hat\lambda_i} and standard deviation
#'     \deqn{\sqrt{\hat\Xi^*_{i,t,k} / I + \hat\Xi^{**}_{i,t,k} / T},}
#'     capturing cross-sectional and time-series uncertainty in factor estimation.
#'   \item Computes the IRF at \eqn{(\beta^{(b)}, \tilde C^{(b)})}.
#' }
#' The median and the \eqn{(1-\mathrm{level})/2} and
#' \eqn{(1+\mathrm{level})/2} quantiles across draws give the point estimate
#' and confidence bands.
#'
#' @param fit An object of class \code{"pvarife_result"}.
#' @param n_periods Positive integer. Number of IRF horizons.
#' @param shock Positive integer. Index of the structural shock (default 1).
#' @param diff_vars Integer vector. Variables to cumulate (default none).
#' @param identification Character. \code{"short_run"} (default) or
#'   \code{"long_run"}. Passed to \code{\link{compute_irf}} for each draw.
#'   The median IRF returned is implicitly bias-corrected regardless of this
#'   choice (draws are centred on \eqn{\hat\beta - \hat b}).
#' @param n_draw Positive integer. Number of simulation draws (default 500).
#' @param level Numeric in (0, 1). Confidence level (default 0.95).
#' @param seed Optional integer seed for reproducibility.
#'
#' @return An object of class \code{"pvarife_bands"} with components:
#'   \describe{
#'     \item{irf}{Point estimate (median across draws, bias-corrected),
#'       \eqn{K \times n\_periods}.}
#'     \item{lower}{Lower confidence bound, \eqn{K \times n\_periods}.}
#'     \item{upper}{Upper confidence bound, \eqn{K \times n\_periods}.}
#'     \item{level}{The confidence level used.}
#'     \item{method}{\code{"asymptotic"}.}
#'   }
#'
#' @references
#' Tugan, M. (2021). Panel VAR models with interactive fixed effects.
#' \emph{Econometrics Journal}, 24, 225--246.
#' \doi{10.1093/ectj/utaa021}
#'
#' @examples
#' \dontrun{
#' sim   <- sim_pvarife(n_units = 30, n_time = 20, n_vars = 2,
#'                      n_lags = 1, n_factors = 1, seed = 1)
#' fit   <- pvarife(sim$y, n_lags = 1, n_factors = 1, n_out = 10, n_in = 5)
#' bands <- irf_bands(fit, n_periods = 8, n_draw = 200, seed = 42)
#' plot(bands)
#' }
#'
#' @seealso \code{\link{bootstrap_irf_bands}}, \code{\link{compute_irf}}
#'
#' @export
irf_bands <- function(fit, n_periods, shock = 1L, diff_vars = integer(0),
                      identification = c("short_run", "long_run"),
                      n_draw = 500L, level = 0.95, seed = NULL) {
  stopifnot(inherits(fit, "pvarife_result"))
  identification <- match.arg(identification)
  n_periods <- as.integer(n_periods)
  n_draw    <- as.integer(n_draw)
  if (!is.null(seed)) set.seed(seed)

  n_units   <- fit$n_units
  n_time    <- fit$n_time
  n_vars    <- fit$n_vars
  n_factors <- fit$n_factors
  factors_mat <- fit$factors_mat  # TK x Kr
  ff          <- fit$ff           # T x r
  loadings    <- fit$loadings     # Kr x n_units
  sigma       <- fit$sigma        # K x K
  y_c         <- fit$y_c
  z_c         <- fit$z_c
  i_obs       <- fit$i_obs
  n_time_i    <- fit$n_time_i

  n_rows <- nrow(factors_mat)  # TK

  # Compute asymptotic bias and variance
  avar <- asymptotic_var(fit)
  beta_bias <- avar$bias
  beta_var  <- avar$variance

  # --- Estimated common component C_hat[NT, n_units] ---
  c_hat <- matrix(NA_real_, nrow = n_rows, ncol = n_units)
  for (ii in seq_len(n_units)) {
    c_hat[, ii] <- as.numeric(factors_mat %*% loadings[, ii])
  }

  # --- XiStar: (1/I) * uncertainty from loading estimation ---
  # Faithful to ConfidenceBandforIRs.m (Tugan 2021):
  #   XiStar[N(t-1)+n, c] = ( lambda_c' ((1/I) lambda lambda')^{-1} lambda_c )
  #                         * Sigma_e[n, n]
  # where `loadings` is the (Kr x I) raw loading matrix and lambda_c its c-th
  # column. The quadratic form is a single scalar per unit (mixing all Kr
  # loadings); it is constant across time t and is scaled by the n-th diagonal
  # element of Sigma for variable n.
  m_lam     <- (loadings %*% t(loadings)) / n_units   # (Kr) x (Kr)
  m_lam_inv <- solve(m_lam)
  q_unit <- vapply(
    seq_len(n_units),
    function(ii) as.numeric(crossprod(loadings[, ii], m_lam_inv %*% loadings[, ii])),
    numeric(1L)
  )                                                    # length I
  sig_diag <- diag(sigma)                              # length K

  xi_star <- matrix(NA_real_, nrow = n_rows, ncol = n_units)
  for (tt in seq_len(n_time)) {
    for (kk in seq_len(n_vars)) {
      row_idx <- (tt - 1L) * n_vars + kk
      xi_star[row_idx, ] <- q_unit * sig_diag[kk]
    }
  }

  # --- XiStarStar: (1/T) * uncertainty from factor estimation ---
  # A = (1/T) * sum_{t,n} f_{tn}' (Sigma_{n.}' F_t)  [r x r... no, Nr x Nr]
  # Simplified: A_mat = (1/T) sum_t f_t' Sigma f_t (r x r)
  a_mat <- matrix(0.0, nrow = n_factors * n_vars, ncol = n_factors * n_vars)
  for (tt in seq_len(n_time)) {
    for (kk in seq_len(n_vars)) {
      row_idx <- (tt - 1L) * n_vars + kk
      f_tnk <- factors_mat[row_idx, , drop = FALSE]  # 1 x Kr
      sig_krow <- sigma[kk, , drop = FALSE]           # 1 x K
      f_t_block <- factors_mat[((tt - 1L) * n_vars + 1L):(tt * n_vars), , drop = FALSE]
      a_mat <- a_mat + t(f_tnk) %*% (sig_krow %*% f_t_block)
    }
  }
  a_mat <- a_mat / n_time

  xi_starstar <- matrix(NA_real_, nrow = n_rows, ncol = n_units)
  for (tt in seq_len(n_time)) {
    for (kk in seq_len(n_vars)) {
      row_idx <- (tt - 1L) * n_vars + kk
      f_obs <- factors_mat[row_idx, , drop = FALSE]  # 1 x Kr
      val   <- as.numeric(f_obs %*% a_mat %*% t(f_obs))
      xi_starstar[row_idx, ] <- val  # same for all units
    }
  }

  # --- Simulation draws ---
  beta_mean <- fit$beta - matrix(beta_bias, ncol = 1L)
  irf_draws <- array(NA_real_, dim = c(n_vars, n_periods, n_draw))

  for (bb in seq_len(n_draw)) {
    # Draw common component
    sd_cc <- sqrt(xi_star / n_units + xi_starstar / n_time)
    cc_b  <- matrix(stats::rnorm(n_rows * n_units,
                                  mean = as.numeric(c_hat),
                                  sd   = as.numeric(sd_cc)),
                    nrow = n_rows, ncol = n_units)

    # Draw beta
    beta_b <- matrix(mvtnorm::rmvnorm(1L, mean = as.numeric(beta_mean),
                                      sigma = beta_var),
                     ncol = 1L)

    ir_b <- .irf_worker(beta_b, y_c, z_c, cc_b,
                         n_vars, n_lags = fit$n_lags,
                         n_time_i, i_obs, n_units,
                         n_periods, shock, diff_vars,
                         identification = identification)
    if (!is.null(ir_b)) irf_draws[, , bb] <- ir_b
  }

  # Quantiles: sort along draw dimension
  alpha_lo <- (1.0 - level) / 2.0
  alpha_hi <- 1.0 - alpha_lo

  ir_sorted  <- apply(irf_draws, c(1L, 2L), function(x) sort(x[!is.na(x)]))
  n_valid    <- apply(irf_draws, c(1L, 2L), function(x) sum(!is.na(x)))

  irf_med   <- apply(irf_draws, c(1L, 2L), stats::median, na.rm = TRUE)
  irf_lower <- apply(irf_draws, c(1L, 2L), stats::quantile,
                     probs = alpha_lo, na.rm = TRUE)
  irf_upper <- apply(irf_draws, c(1L, 2L), stats::quantile,
                     probs = alpha_hi, na.rm = TRUE)

  structure(
    list(
      irf    = irf_med,
      lower  = irf_lower,
      upper  = irf_upper,
      level  = level,
      method = "asymptotic"
    ),
    class = "pvarife_bands"
  )
}


#' Classical residual bootstrap confidence bands for IRFs
#'
#' Constructs confidence bands for structural impulse responses via a classical
#' residual bootstrap. For each bootstrap replicate, residuals are resampled
#' (by time index, preserving cross-sectional dependence), a new panel is
#' constructed, the full \code{pvarife} model is re-estimated, and IRFs are
#' computed. This is computationally expensive but does not rely on asymptotic
#' approximations.
#'
#' @param fit An object of class \code{"pvarife_result"}.
#' @param n_periods Positive integer. Number of IRF horizons.
#' @param shock Positive integer. Index of the structural shock (default 1).
#' @param diff_vars Integer vector. Variables to cumulate (default none).
#' @param identification Character. \code{"short_run"} (default) or
#'   \code{"long_run"}. Passed to \code{\link{compute_irf}} for each bootstrap
#'   replicate.
#' @param n_boot Positive integer. Number of bootstrap replicates (default 200).
#' @param level Numeric in (0, 1). Confidence level (default 0.95).
#' @param seed Optional integer seed for reproducibility.
#'
#' @return An object of class \code{"pvarife_bands"} with components
#'   \code{irf}, \code{lower}, \code{upper}, \code{level}, and
#'   \code{method = "bootstrap"}.
#'
#' @examples
#' \dontrun{
#' sim   <- sim_pvarife(n_units = 20, n_time = 15, n_vars = 2,
#'                      n_lags = 1, n_factors = 1, seed = 1)
#' fit   <- pvarife(sim$y, n_lags = 1, n_factors = 1, n_out = 5, n_in = 3)
#' bands <- bootstrap_irf_bands(fit, n_periods = 6, n_boot = 50, seed = 42)
#' plot(bands)
#' }
#'
#' @seealso \code{\link{irf_bands}}, \code{\link{compute_irf}}
#'
#' @export
bootstrap_irf_bands <- function(fit, n_periods, shock = 1L,
                                 diff_vars = integer(0),
                                 identification = c("short_run", "long_run"),
                                 n_boot = 200L, level = 0.95, seed = NULL) {
  stopifnot(inherits(fit, "pvarife_result"))
  identification <- match.arg(identification)
  n_periods <- as.integer(n_periods)
  n_boot    <- as.integer(n_boot)
  if (!is.null(seed)) set.seed(seed)

  n_units   <- fit$n_units
  n_time    <- fit$n_time
  n_vars    <- fit$n_vars
  n_lags    <- fit$n_lags
  n_factors <- fit$n_factors
  n_out     <- 50L   # use same defaults for re-estimation
  n_in      <- 10L

  y_c         <- fit$y_c
  z_c         <- fit$z_c
  i_obs       <- fit$i_obs
  factors_mat <- fit$factors_mat
  loadings    <- fit$loadings
  beta        <- fit$beta
  n_time_i    <- fit$n_time_i
  n_rows      <- nrow(y_c)

  # Compute residuals (T_i x K per unit)
  u_list <- vector("list", n_units)
  for (ii in seq_len(n_units)) {
    obs <- which(i_obs[, ii] == 1L)
    if (length(obs) == 0L) { u_list[[ii]] <- NULL; next }
    yy   <- matrix(y_c[obs, 1L, ii], ncol = 1L)
    zz   <- matrix(z_c[obs, , ii], nrow = length(obs))
    ci   <- matrix(factors_mat[obs, ] %*% loadings[, ii], ncol = 1L)
    uu   <- yy - zz %*% beta - ci
    tt_i <- n_time_i[ii]
    u_list[[ii]] <- matrix(uu, nrow = n_vars, ncol = tt_i)  # K x T_i
    u_list[[ii]] <- t(u_list[[ii]])                          # T_i x K
  }

  irf_draws <- array(NA_real_, dim = c(n_vars, n_periods, n_boot))

  for (bb in seq_len(n_boot)) {
    # Build bootstrapped y_arr: I x T x K
    y_boot <- fit$y_c * NA_real_  # initialise as NA, same structure as y_c
    y_arr_boot <- array(NA_real_, dim = c(n_units, n_time, n_vars))

    for (ii in seq_len(n_units)) {
      obs <- which(i_obs[, ii] == 1L)
      if (length(obs) == 0L || is.null(u_list[[ii]])) next

      tt_i   <- n_time_i[ii]
      resamp <- sample.int(tt_i, tt_i, replace = TRUE)
      u_boot <- u_list[[ii]][resamp, , drop = FALSE]  # T_i x K (resampled)

      # Reconstruct y = Z*beta + F*lambda + u_boot (observed positions)
      zz   <- matrix(z_c[obs, , ii], nrow = length(obs))
      ci   <- matrix(factors_mat[obs, ] %*% loadings[, ii], ncol = 1L)
      yy_b <- zz %*% beta + ci +
        matrix(t(u_boot), ncol = 1L)  # NT_i x 1

      y_arr_boot_ii <- array(NA_real_, dim = c(n_time, n_vars))
      # Map back: obs rows -> time indices
      obs_t_idx <- ceiling(obs / n_vars)
      for (jj in seq_along(obs_t_idx)) {
        tt_idx <- obs_t_idx[jj]
        kk_idx <- obs[jj] - (tt_idx - 1L) * n_vars
        y_arr_boot[ii, tt_idx, kk_idx] <- yy_b[jj]
      }
    }

    fit_b <- tryCatch(
      pvarife(y_arr_boot, n_lags = n_lags, n_factors = n_factors,
              n_out = n_out, n_in = n_in),
      error = function(e) NULL
    )
    if (is.null(fit_b)) next

    ir_b <- tryCatch(
      compute_irf(fit_b, n_periods = n_periods, shock = shock,
                  diff_vars = diff_vars, identification = identification),
      error = function(e) NULL
    )
    if (!is.null(ir_b)) irf_draws[, , bb] <- ir_b
  }

  alpha_lo <- (1.0 - level) / 2.0
  alpha_hi <- 1.0 - alpha_lo

  irf_med   <- apply(irf_draws, c(1L, 2L), stats::median, na.rm = TRUE)
  irf_lower <- apply(irf_draws, c(1L, 2L), stats::quantile,
                     probs = alpha_lo, na.rm = TRUE)
  irf_upper <- apply(irf_draws, c(1L, 2L), stats::quantile,
                     probs = alpha_hi, na.rm = TRUE)

  structure(
    list(
      irf    = irf_med,
      lower  = irf_lower,
      upper  = irf_upper,
      level  = level,
      method = "bootstrap"
    ),
    class = "pvarife_bands"
  )
}
