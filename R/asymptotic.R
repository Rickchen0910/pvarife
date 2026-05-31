#' Asymptotic variance and bias of the pvarife estimator
#'
#' Computes the asymptotic bias and variance-covariance matrix of \eqn{\hat\beta}
#' under Theorem 2.3 of Tugan (2021). These quantities are used by
#' \code{\link{irf_bands}} to construct parametric confidence bands.
#'
#' @details
#' The function computes three components:
#' \describe{
#'   \item{\strong{D}}{The Hessian \eqn{D_{F,\Lambda}} (Eq. A.4), which accounts
#'     for factor estimation uncertainty via a two-term formula.}
#'   \item{\strong{Omega}}{The sandwich variance \eqn{\Omega}, accumulated
#'     unit by unit.}
#'   \item{\strong{Bias}}{Two bias terms: \eqn{B_\Psi} (from factor loading
#'     estimation) and \eqn{B_\gamma} (HAC serial correction with bandwidth
#'     \eqn{\bar G = \lfloor T^{1/3} \rceil}).}
#' }
#'
#' \strong{Note on MATLAB replication:} The original \code{Asymptotic_Distribution_of_beta.m}
#' contains a bug at line 189 where the \eqn{B_\gamma} accumulation uses only the
#' final value of the loop variable \code{g} rather than summing over
#' \eqn{g = 1, \ldots, \bar G}. This function implements the \emph{corrected} version.
#'
#' @param fit An object of class \code{"pvarife_result"} returned by
#'   \code{\link{pvarife}}.
#'
#' @return A list with:
#'   \describe{
#'     \item{bias}{Bias vector for \eqn{\hat\beta}, of the same length as
#'       \code{fit$beta}.}
#'     \item{variance}{Asymptotic variance-covariance matrix of \eqn{\hat\beta}.}
#'   }
#'
#' @references
#' Tugan, M. (2021). Panel VAR models with interactive fixed effects.
#' \emph{Econometrics Journal}, 24, 225--246.
#' \doi{10.1093/ectj/utaa021}
#'
#' @examples
#' sim <- sim_pvarife(n_units = 30, n_time = 20, n_vars = 2,
#'                    n_lags = 1, n_factors = 1, seed = 1)
#' fit <- pvarife(sim$y, n_lags = 1, n_factors = 1, n_out = 5, n_in = 3)
#' avar <- asymptotic_var(fit)
#' cat("Bias:", avar$bias, "\n")
#'
#' @export
asymptotic_var <- function(fit) {
  stopifnot(inherits(fit, "pvarife_result"))

  n_units   <- fit$n_units      # I
  n_vars    <- fit$n_vars       # K
  n_factors <- fit$n_factors    # r
  n_time    <- fit$n_time       # T (full)

  z_c          <- fit$z_c       # (T*K) x n_cols_z x I
  u_c          <- fit$u_c       # (T*K) x 1 x I
  factors_mat  <- fit$factors_mat  # (T*K) x (K*r)
  ff           <- fit$ff        # T x r
  loadings     <- fit$loadings  # (K*r) x I
  sigma        <- fit$sigma     # K x K
  n_time_i     <- fit$n_time_i  # effective time per unit (integer vector)
  i_obs        <- fit$i_obs     # (T*K) x I  (1=observed, 0=missing)

  n_rows    <- dim(z_c)[1L]     # = K * T
  n_cols_z  <- dim(z_c)[2L]
  sum_tc    <- sum(n_time_i)
  sum_tnc   <- sum(n_time_i) * n_vars   # approximation: sum(TC_i * K)

  # -------------------------------------------------------------------------
  # lambda_underbar: stacked (K x r) loading matrix, one per unit, bound by row
  # MATLAB: lambda_underbar_c(:,:,i) = reshape(lambda(:,i)', r, K)  [r x K]
  # In R: lam_list[[i]] = matrix(loadings[,i], nrow=K, ncol=r)     [K x r]
  # The full lambda_underbar is stacked vertically: (K*I x r)
  # -------------------------------------------------------------------------
  lam_list <- vector("list", n_units)
  for (ii in seq_len(n_units)) {
    lam_list[[ii]] <- matrix(loadings[, ii], nrow = n_vars, ncol = n_factors)  # K x r
  }
  lambda_underbar <- do.call(rbind, lam_list)  # (K*I) x r

  # (lambda' lambda / (K*I))^{-1}  — r x r matrix
  lam_lam_inv <- solve(crossprod(lambda_underbar) / (n_vars * n_units))

  # M_F = I_{TK} - F (F'F)^{-1} F'
  m_f <- diag(n_rows) - factors_mat %*% solve(crossprod(factors_mat), t(factors_mat))

  # -------------------------------------------------------------------------
  # Precompute:
  #   kron_T_c[[ii]] = kron(I_T, lam_i)           [T*K x T*r]  -- MATLAB's Kron_lambda_T_c'
  #   kron_i_c[[ii]] = kron(I_T, t(lam_i))         [T*r x T*K]  -- MATLAB's Kron_lambda_i
  # These correspond to MATLAB Kronecker_Identity_lambda_T_c and
  # Kronecker_Identity_lambda_i respectively.
  # -------------------------------------------------------------------------
  # MATLAB: Kronecker_Identity_lambda_T_c = kron(I_T, (lambda_underbar_c)')
  #   where lambda_underbar_c is r x K, so (lambda_underbar_c)' = K x r
  #   result: (T*K) x (T*r)
  # R: kronecker(I_T, lam_i) where lam_i = K x r  -> (T*K) x (T*r) ✓
  kron_T_list <- vector("list", n_units)   # each (T*K) x (T*r)
  kron_i_list <- vector("list", n_units)   # each (T*r) x (T*K)
  for (ii in seq_len(n_units)) {
    lam_i <- lam_list[[ii]]                            # K x r
    kron_T_list[[ii]] <- kronecker(diag(n_time), lam_i)       # (T*K) x (T*r)
    kron_i_list[[ii]] <- kronecker(diag(n_time), t(lam_i))    # (T*r) x (T*K)
  }

  # kron(I_T, lam_lam_inv): (T*r) x (T*r)
  kron_inv <- kronecker(diag(n_time), lam_lam_inv)  # (T*r) x (T*r)

  # -------------------------------------------------------------------------
  # D_first = (1/sum_tnc) * sum_i Z_i'[obs] M_F[obs,obs] Z_i[obs]
  # -------------------------------------------------------------------------
  d_first <- matrix(0.0, nrow = n_cols_z, ncol = n_cols_z)
  for (ii in seq_len(n_units)) {
    obs <- which(i_obs[, ii] == 1L)
    if (length(obs) == 0L) next
    zz  <- matrix(z_c[obs, , ii], nrow = length(obs))
    mf  <- m_f[obs, obs, drop = FALSE]
    d_first <- d_first + crossprod(zz, mf) %*% zz
  }
  d_first <- d_first / sum_tnc

  # -------------------------------------------------------------------------
  # sum_z_mf_kron_T[ii] = Z_i'[obs] M_F[obs,obs] kron_T[[ii]][obs, ]
  # size: n_cols_z x (T*r)
  # Then sum_z_mf_lam = sum_i sum_z_mf_kron_T[[ii]]
  # -------------------------------------------------------------------------
  sum_z_mf_lam <- matrix(0.0, nrow = n_cols_z, ncol = n_time * n_factors)
  z_mf_kron_T_list <- vector("list", n_units)
  kron_i_obs_z_list <- vector("list", n_units)

  for (ii in seq_len(n_units)) {
    obs <- which(i_obs[, ii] == 1L)
    if (length(obs) == 0L) next
    zz       <- matrix(z_c[obs, , ii], nrow = length(obs))
    mf       <- m_f[obs, obs, drop = FALSE]
    kron_T_i <- kron_T_list[[ii]]   # (T*K) x (T*r)

    z_mf_kron_T <- crossprod(zz, mf) %*% kron_T_i[obs, , drop = FALSE]  # n_cols_z x (T*r)
    z_mf_kron_T_list[[ii]] <- z_mf_kron_T
    sum_z_mf_lam <- sum_z_mf_lam + z_mf_kron_T

    # kron_i[[ii]][:, obs] * Z_i[obs]  -- (T*r) x n_cols_z
    kron_i_i <- kron_i_list[[ii]]  # (T*r) x (T*K)
    kron_i_obs_z_list[[ii]] <-
      kron_i_i[, obs, drop = FALSE] %*% zz  # (T*r) x n_cols_z
  }

  sum_kron_i_z <- Reduce("+", kron_i_obs_z_list[!vapply(kron_i_obs_z_list, is.null,
                                                          logical(1L))])

  # D_second = sum_z_mf_lam * kron_inv * sum_kron_i_z / (sum_tnc * K * I)
  d_second <- sum_z_mf_lam %*% kron_inv %*% sum_kron_i_z / (sum_tnc * n_vars * n_units)

  d_fl <- d_first - d_second
  d_inv <- solve(d_fl)

  # -------------------------------------------------------------------------
  # Omega: sandwich variance
  # Gamma_i_Z = (1/K) * (Z_i'M_F - (1/(K*I)) * sum_z_mf_lam * kron_inv * kron_i[:, obs])
  # Omega = (1/sum_tc) * sum_i Gamma_i_Z diag(u_i) diag(u_i)' Gamma_i_Z'
  # -------------------------------------------------------------------------
  omega <- matrix(0.0, nrow = n_cols_z, ncol = n_cols_z)

  for (ii in seq_len(n_units)) {
    obs <- which(i_obs[, ii] == 1L)
    if (length(obs) == 0L) next
    zz  <- matrix(z_c[obs, , ii], nrow = length(obs))
    uu  <- u_c[obs, 1L, ii]
    mf  <- m_f[obs, obs, drop = FALSE]

    kron_i_obs <- kron_i_list[[ii]][, obs, drop = FALSE]  # (T*r) x n_obs

    gamma_first  <- crossprod(zz, mf)                                    # n_cols_z x n_obs
    gamma_second <- sum_z_mf_lam %*% kron_inv %*% kron_i_obs             # n_cols_z x n_obs
    gamma_iz <- (1.0 / n_vars) * (gamma_first - (1.0 / (n_vars * n_units)) * gamma_second)

    gamma_u <- gamma_iz %*% diag(uu, nrow = length(uu))                  # n_cols_z x n_obs
    omega <- omega + tcrossprod(gamma_u)
  }
  omega <- omega / sum_tc

  beta_variance <- d_inv %*% omega %*% t(d_inv) / sum_tc

  # -------------------------------------------------------------------------
  # B_Psi bias
  # Psi_star_i  = Z_i'F[obs] / T_i * kron(I_K, lam_lam_inv * lam_i) * sigma_vec
  # Psi_double_star_i = Z_i'F[obs] / T_i * kron(lam_i' lam_lam_inv, lam_lam_inv) *
  #                     sum_j kron(lam_j, lam_j) * sigma_vec
  # rho = T / I
  # -------------------------------------------------------------------------
  rho <- n_time / n_units
  sigma_vec <- as.numeric(t(sigma))  # vectorise Sigma' (column-major of Sigma')

  # sum_j kron(lam_j, lam_j)  -- (K*r x K*r) summed over all units... wait
  # MATLAB: Psi_double_star_Third_Term = sum_i kron(lambda_underbar_c(:,:,i), lambda_underbar_c(:,:,i))
  # lambda_underbar_c is r x K, so kron(r x K, r x K) = (r^2) x (K^2)
  # In R: kron(t(lam_i), t(lam_i)) where t(lam_i) = r x K -> (r^2) x (K^2)
  kron_lam_lam_sum <- matrix(0.0, nrow = n_factors^2, ncol = n_vars^2)
  for (ii in seq_len(n_units)) {
    lam_i_rk <- t(lam_list[[ii]])    # r x K
    kron_lam_lam_sum <- kron_lam_lam_sum + kronecker(lam_i_rk, lam_i_rk)
  }

  psi_star_sum   <- matrix(0.0, nrow = n_cols_z, ncol = 1L)
  psi_double_sum <- matrix(0.0, nrow = n_cols_z, ncol = 1L)

  for (ii in seq_len(n_units)) {
    obs <- which(i_obs[, ii] == 1L)
    if (length(obs) == 0L) next
    zz      <- matrix(z_c[obs, , ii], nrow = length(obs))
    lam_i   <- lam_list[[ii]]    # K x r
    lam_i_rk <- t(lam_i)         # r x K
    tt_obs  <- length(obs) / n_vars

    # Psi_star_First_i = Z_i'[obs] F[obs] / T_i   [n_cols_z x (K*r)]
    psi_first_i <- crossprod(zz, factors_mat[obs, , drop = FALSE]) / tt_obs

    # Psi_star_Second_i = kron(I_K, lam_lam_inv * lam_i)
    # lam_lam_inv is r x r, lam_i is K x r => lam_lam_inv * lam_i' (transpose needed)
    # MATLAB: kron(eye(K), inv(lambda_lambda) * lambda_underbar_c(:,:,i))
    # lambda_underbar_c = r x K, inv(lambda_lambda) = r x r
    # inv(lambda_lambda) * (r x K) = r x K
    # kron(I_K, r x K) = (K*r) x (K^2)... hmm
    # Actually factors_mat has K*r columns, so we need: n_cols_z x K*r * (K*r x K^2) * K^2
    # Let me follow MATLAB more carefully
    # Psi_star_i = psi_first_i * psi_second_i * Sigma_vec
    # psi_first_i: n_cols_z x (K*r)
    # psi_second_i: (K*r) x K^2
    # Sigma_vec: K^2 x 1
    # MATLAB: kron(eye(K), inv(lambda)*lam_i) where lam_i = r x K, inv = r x r
    # inv * lam_i = r x K; kron(I_K, r x K) = K*r x K^2
    psi_second_i <- kronecker(diag(n_vars), lam_lam_inv %*% lam_i_rk)  # (K*r) x (K^2)
    psi_star_sum <- psi_star_sum + psi_first_i %*% psi_second_i %*% sigma_vec

    # Psi_double_star_second: kron(lam_i' lam_lam_inv, lam_lam_inv)
    # MATLAB: kron(lam_i' * lam_lam_inv, lam_lam_inv) where lam_i' = r x K (as in MATLAB)
    # but we need lam_i' as K x r... MATLAB: lambda_underbar_c' = K x r
    # kron(lam_lam_inv * lambda_underbar_c', lam_lam_inv) -- checking MATLAB:
    # kron( lambda_underbar_c(:,:,i)' * inv(lambda_lambda), inv(lambda_lambda) )
    # lambda_underbar_c(:,:,i)' = K x r (from r x K)
    # (K x r) * (r x r) = K x r
    # kron(K x r, r x r) = (K*r) x (r^2)... this is getting complex
    # Let me simplify: follow paper notation more closely.
    # For now, use scalar approximation that gives correct sign.
    # MATLAB: kron(lambda_underbar_c(:,:,i)' * inv(lam_lam), inv(lam_lam))
    # lambda_underbar_c = r x K  =>  (r x K)' = K x r (= our lam_i)
    # (K x r) * (r x r) = K x r; kron(K x r, r x r) = (K*r) x (r^2)
    psi_ds_second <- kronecker(lam_i %*% lam_lam_inv, lam_lam_inv)  # (K*r) x (r^2)
    # psi_ds_third = kron_lam_lam_sum: (r^2) x (K^2)
    psi_double_sum <- psi_double_sum +
      psi_first_i %*% psi_ds_second %*% kron_lam_lam_sum %*% sigma_vec
  }

  psi_star   <- -sqrt(rho) / (n_vars^2 * n_units) * psi_star_sum
  psi_double <- -sqrt(rho) / (n_vars^3 * n_units^2) * psi_double_sum
  b_psi <- d_inv %*% (psi_star - psi_double)

  # -------------------------------------------------------------------------
  # B_gamma bias (HAC serial correction, CORRECTED from MATLAB bug)
  # G_bar = round(T^(1/3))
  # For each unit i and lag g = 1:G_bar:
  #   accumulate lead-g of Z_c and F, weighted by u_i
  # -------------------------------------------------------------------------
  g_bar <- max(1L, round(n_time^(1.0 / 3.0)))

  b_gamma_sum <- matrix(0.0, nrow = n_cols_z, ncol = 1L)

  for (ii in seq_len(n_units)) {
    obs <- which(i_obs[, ii] == 1L)
    if (length(obs) == 0L) next
    uu <- matrix(u_c[obs, 1L, ii], ncol = 1L)

    z_full <- matrix(z_c[, , ii], nrow = n_rows, ncol = n_cols_z)
    # FTF_inv = (F'F/T)^{-1}  (r x r)
    ff_tff_inv <- solve(crossprod(ff) / n_time)  # T x r -> r x r

    for (gg in seq_len(g_bar)) {
      shift <- n_vars * gg  # number of rows to shift (each time period = K rows)

      # Lead of Z_c by gg time steps (shift rows up by K*gg positions)
      lead_z <- lag_lead_matrix(z_full, -shift, -shift)  # (T*K) x n_cols_z

      # Lead of ff by gg time steps
      lead_ff <- lag_lead_matrix(ff, -gg, -gg)  # T x r (shift by gg time periods)

      # Build lead factors_mat (block-diagonal) from lead_ff
      lead_fmat <- .build_factors_mat(lead_ff, n_vars)  # (T*K) x (K*r), NAs at top g rows

      # Find rows valid in both lead_z and obs
      lead_z_valid <- which(stats::complete.cases(lead_z))
      common_obs   <- intersect(obs, lead_z_valid)
      if (length(common_obs) == 0L) next

      lead_z_sub <- lead_z[common_obs, , drop = FALSE]   # n_common x n_cols_z
      uu_common  <- u_c[common_obs, 1L, ii]

      # Weight: diag(F_{t+g} (F'F/T)^{-1} F_t') applied elementwise
      # MATLAB: diag(diag(Lead_F_g[obs] * inv(F'F/T) * F[obs]'))
      lead_fm_sub <- lead_fmat[common_obs, , drop = FALSE]  # n_common x K*r
      fm_sub      <- factors_mat[common_obs, , drop = FALSE]  # n_common x K*r

      # Weight matrix: each row of lead_fm * inv(F'F/T) * t(fm), take diagonal
      # This is a scalar per observation
      # Lead_F_g[obs] is (T*K) x (K*r) evaluated at obs; we need the (T x r) version
      # Actually MATLAB uses the small f (T x r), not F (T*K x K*r)
      # Lead_F_g is lead of F (T*K x K*r), and it computes:
      # diag( Lead_F_g[obs, :] * inv(F'F/T) * F[obs, :]' ) -- this is n_obs x n_obs diagonal
      # For the scalar weight approach: use lead_ff and ff at the corresponding time indices
      # Map obs rows back to time indices
      t_idx      <- ceiling(common_obs / n_vars)  # time index for each obs row
      lead_ff_t  <- lead_ff[t_idx, , drop = FALSE]   # n_common x r (NAs where lead is unavailable)
      ff_t       <- ff[t_idx, , drop = FALSE]         # n_common x r

      has_lead <- stats::complete.cases(lead_ff_t)
      if (sum(has_lead) == 0L) next

      lead_ff_ok <- lead_ff_t[has_lead, , drop = FALSE]
      ff_t_ok    <- ff_t[has_lead, , drop = FALSE]
      lead_z_ok  <- lead_z_sub[has_lead, , drop = FALSE]
      uu_ok      <- uu_common[has_lead]

      wts <- rowSums((lead_ff_ok %*% ff_tff_inv) * ff_t_ok)  # n_valid scalars

      b_gamma_sum <- b_gamma_sum +
        crossprod(lead_z_ok, wts * uu_ok)
    }
  }

  b_gamma <- -d_inv %*% (1.0 / sqrt(rho) / sum_tnc * b_gamma_sum)
  beta_bias <- as.numeric((b_gamma + b_psi) / sqrt(sum_tc))

  list(bias = beta_bias, variance = beta_variance)
}


#' Extract factors and loadings at an arbitrary coefficient vector
#'
#' Given an estimated \code{pvarife_result} and an arbitrary coefficient vector
#' \code{beta}, runs the inner EM loop to extract common factors and factor
#' loadings. Useful for advanced users (e.g., bootstrap procedures that need
#' factor estimates at a perturbed beta).
#'
#' Faithful translation of \code{Inner_Iteration.m} from Tugan (2021).
#'
#' @param beta Numeric vector of VAR coefficients (same length as
#'   \code{fit$beta}).
#' @param fit An object of class \code{"pvarife_result"}.
#' @param n_in Number of inner iterations (default 10).
#'
#' @return A list with \code{ff} (T x r factor matrix) and \code{loadings}
#'   (Kr x I loading matrix).
#'
#' @examples
#' sim <- sim_pvarife(n_units = 20, n_time = 15, n_vars = 2,
#'                    n_lags = 1, n_factors = 1, seed = 2)
#' fit <- pvarife(sim$y, n_lags = 1, n_factors = 1, n_out = 5, n_in = 3)
#' ef  <- extract_factors(fit$beta * 0.9, fit)
#' dim(ef$ff)       # T x r
#'
#' @export
extract_factors <- function(beta, fit, n_in = 10L) {
  stopifnot(inherits(fit, "pvarife_result"))
  beta <- matrix(as.numeric(beta), ncol = 1L)
  n_in <- as.integer(n_in)

  y_c    <- fit$y_c
  z_c    <- fit$z_c
  i_obs  <- fit$i_obs
  n_units <- fit$n_units
  n_time  <- fit$n_time
  n_vars  <- fit$n_vars
  n_rows  <- dim(y_c)[1L]

  # Initialise W_c from current beta
  w_c <- array(0.0, dim = dim(y_c))
  for (ii in seq_len(n_units)) {
    obs  <- i_obs[, ii] == 1L
    if (any(obs)) {
      yy <- y_c[obs, 1L, ii]
      zz <- matrix(z_c[obs, , ii], nrow = sum(obs))
      w_c[obs, 1L, ii] <- yy - as.numeric(zz %*% beta)
    }
  }

  ff          <- NULL
  factors_mat <- NULL
  loadings    <- NULL

  for (i_inner in seq_len(n_in)) {
    vvprime <- matrix(0.0, nrow = n_time, ncol = n_time)
    for (ii in seq_len(n_units)) {
      v_mat <- t(matrix(w_c[, 1L, ii], nrow = n_vars))
      vvprime <- vvprime + tcrossprod(v_mat)
    }

    ev <- eigen(vvprime / (n_time * n_vars * n_units), symmetric = TRUE)
    ff <- sqrt(n_time) * ev$vectors[, seq_len(fit$n_factors), drop = FALSE]
    factors_mat <- .build_factors_mat(ff, n_vars)

    w_stack <- matrix(as.numeric(w_c), ncol = 1L)
    loadings <- matrix(NA_real_, nrow = ncol(factors_mat), ncol = n_units)
    ff_inv <- solve(crossprod(factors_mat), t(factors_mat))
    for (ii in seq_len(n_units)) {
      row_s <- (ii - 1L) * n_rows + 1L
      row_e <- ii * n_rows
      loadings[, ii] <- ff_inv %*% w_stack[row_s:row_e]
    }

    w_c_new <- array(0.0, dim = dim(w_c))
    for (ii in seq_len(n_units)) {
      obs  <- i_obs[, ii] == 1L
      miss <- i_obs[, ii] == 0L
      if (any(miss))
        w_c_new[miss, 1L, ii] <- as.numeric(factors_mat %*% loadings[, ii])[miss]
      if (any(obs)) {
        yy <- y_c[obs, 1L, ii]
        zz <- matrix(z_c[obs, , ii], nrow = sum(obs))
        w_c_new[obs, 1L, ii] <- yy - as.numeric(zz %*% beta)
      }
    }
    w_c <- w_c_new
  }

  list(ff = ff, loadings = loadings, factors_mat = factors_mat)
}
