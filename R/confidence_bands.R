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
#' \strong{Scope of the bands.} This is a faithful implementation of the band
#' construction in Tugan (2021) (\code{ConfidenceBandforIRs.m}): the bands
#' propagate uncertainty in the VAR coefficients \eqn{\beta} and in the
#' estimated common component, but the reduced-form covariance \eqn{\hat\Sigma}
#' is recomputed deterministically for each draw and therefore contributes
#' little draw-to-draw variation. As a result the bands mainly reflect
#' \emph{coefficient} (dynamic) uncertainty and \emph{under-state} uncertainty
#' in the contemporaneous impact at horizon 0 (which is a function of
#' \eqn{\Sigma} alone). Impulse responses are conventionally reported
#' normalised by the shock's own horizon-0 response (see \code{plot} and the
#' \code{normalise_by_h1} argument), which fixes that element to one. For formal
#' inference on the coefficients themselves use \code{\link{asymptotic_var}} or
#' \code{\link{summary.pvarife_result}}, whose Wald intervals attain nominal
#' coverage in simulations.
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
#' sim   <- sim_pvarife(n_units = 20, n_time = 15, n_vars = 2,
#'                      n_lags = 1, n_factors = 1, seed = 1)
#' fit   <- pvarife(sim$y, n_lags = 1, n_factors = 1, n_out = 5, n_in = 3)
#' bands <- irf_bands(fit, n_periods = 6, n_draw = 100, seed = 42)
#' plot(bands)
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


#' Recursive residual bootstrap confidence bands for IRFs
#'
#' Constructs confidence bands for structural impulse responses via a recursive
#' residual bootstrap. For each bootstrap replicate, the idiosyncratic residuals
#' are resampled (by time index, preserving the contemporaneous correlation
#' across variables), a new panel is generated \emph{recursively} from the
#' estimated VAR dynamics and the (fixed) estimated common component, the full
#' \code{pvarife} model is re-estimated, and IRFs are computed. This is
#' computationally expensive but does not rely on asymptotic approximations.
#'
#' @details
#' Each bootstrap panel is generated as
#' \deqn{y^*_{i,t} = c + \sum_{l=1}^{\ell} \Theta_l y^*_{i,t-l}
#'       + \hat F_t \hat\lambda_i + e^*_{i,t},}
#' where \eqn{(c, \Theta_l)} are the estimated coefficients, the common
#' component \eqn{\hat F_t \hat\lambda_i} is held fixed at its estimate, the
#' first \eqn{\ell} periods are initialised at the observed data
#' \eqn{y^*_{i,t} = y_{i,t}}, and \eqn{e^*_{i,t}} are residuals resampled with
#' replacement (whole time rows, so the cross-variable correlation is kept).
#' Because the path is generated recursively, the bootstrap correctly
#' propagates the VAR dynamics — unlike a fixed-design scheme that reuses the
#' original lags.
#'
#' \strong{Scope of the bands.} Only the idiosyncratic errors are resampled;
#' the common component \eqn{\hat F_t \hat\lambda_i} and the reduced-form
#' covariance structure are held at their estimates. The resulting bands
#' therefore capture idiosyncratic-error uncertainty only and \emph{under-cover}
#' when the common factors account for a large share of the variation (as they
#' do in the simulation design of Tugan 2021). This routine is intended as a
#' robustness check; for coefficient inference use \code{\link{asymptotic_var}}
#' or \code{\link{summary.pvarife_result}}, and for the paper's IRF bands use
#' \code{\link{irf_bands}}.
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
#' \donttest{
#' sim   <- sim_pvarife(n_units = 20, n_time = 15, n_vars = 2,
#'                      n_lags = 1, n_factors = 1, seed = 1)
#' fit   <- pvarife(sim$y, n_lags = 1, n_factors = 1, n_out = 5, n_in = 3)
#' bands <- bootstrap_irf_bands(fit, n_periods = 6, n_boot = 20, seed = 42)
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
  beta        <- as.numeric(fit$beta)
  n_time_i    <- fit$n_time_i
  y_arr       <- fit$y_arr   # original I x T x K input

  if (is.null(y_arr))
    stop("fit$y_arr is missing; re-estimate with pvarife() (>= 0.1.1).")

  # --- VAR coefficient matrices Theta_1, ..., Theta_l and intercept c ---
  c_vec <- beta[seq_len(n_vars)]
  alpha <- t(matrix(beta[(n_vars + 1L):length(beta)],
                    nrow = n_vars * n_lags, ncol = n_vars))   # K x (K*L)
  theta_l <- array(0.0, dim = c(n_vars, n_vars, n_lags))
  for (ll in seq_len(n_lags)) {
    theta_l[, , ll] <- alpha[, ((ll - 1L) * n_vars + 1L):(ll * n_vars), drop = FALSE]
  }

  # --- Estimated common component C_i (T x K per unit) and residuals e_i ---
  # C_i stacked = factors_mat %*% loadings[, i]  (TK x 1), reshaped to T x K.
  # Residuals are taken on observed rows only (T_i x K).
  cc_list <- vector("list", n_units)   # T x K (all periods)
  u_list  <- vector("list", n_units)   # complete periods only (T_i x K)
  n_rows_full <- dim(y_c)[1L]
  for (ii in seq_len(n_units)) {
    cc_vec <- as.numeric(factors_mat %*% loadings[, ii])      # TK x 1
    cc_list[[ii]] <- t(matrix(cc_vec, nrow = n_vars, ncol = n_time))  # T x K

    obs <- which(i_obs[, ii] == 1L)
    if (length(obs) == 0L) { u_list[[ii]] <- NULL; next }
    yy <- y_c[obs, 1L, ii]
    zz <- matrix(z_c[obs, , ii], nrow = length(obs))
    # Full-vector reshape, then keep complete periods (robust to partial
    # missingness, where length(obs) is not a multiple of K)
    u_full <- rep(NA_real_, n_rows_full)
    u_full[obs] <- yy - as.numeric(zz %*% beta) - cc_vec[obs]
    u_mat  <- t(matrix(u_full, nrow = n_vars))                # T x K
    keep_t <- stats::complete.cases(u_mat)
    if (!any(keep_t)) { u_list[[ii]] <- NULL; next }
    u_list[[ii]] <- u_mat[keep_t, , drop = FALSE]
  }

  irf_draws <- array(NA_real_, dim = c(n_vars, n_periods, n_boot))

  for (bb in seq_len(n_boot)) {
    # Generate a bootstrap panel recursively: I x T x K
    y_arr_boot <- array(NA_real_, dim = c(n_units, n_time, n_vars))

    for (ii in seq_len(n_units)) {
      if (is.null(u_list[[ii]])) next
      cc_i   <- cc_list[[ii]]                                   # T x K
      e_pool <- u_list[[ii]]                                    # T_i x K
      tt_i   <- nrow(e_pool)                                    # pool size

      # Initialise the first n_lags periods at the observed data
      for (tt in seq_len(n_lags)) y_arr_boot[ii, tt, ] <- y_arr[ii, tt, ]

      # Resample residual time-rows for periods (n_lags+1) .. n_time
      n_gen  <- n_time - n_lags
      resamp <- sample.int(tt_i, n_gen, replace = TRUE)

      for (tt in seq(n_lags + 1L, n_time)) {
        y_t <- c_vec + cc_i[tt, ] + e_pool[resamp[tt - n_lags], ]
        for (ll in seq_len(n_lags)) {
          y_t <- y_t + theta_l[, , ll] %*% y_arr_boot[ii, tt - ll, ]
        }
        y_arr_boot[ii, tt, ] <- as.numeric(y_t)
      }
    }

    fit_b <- tryCatch(
      pvarife(y_arr_boot, n_lags = n_lags, n_factors = n_factors,
              n_out = n_out, n_in = n_in, balanced_init = FALSE),
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
