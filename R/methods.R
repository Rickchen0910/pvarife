# S3 methods for pvarife objects.
# All plotting uses ggplot2 with theme_minimal() for publication quality.

#' Print method for pvarife_result
#'
#' @param x An object of class \code{"pvarife_result"}.
#' @param ... Ignored.
#' @return Invisibly returns \code{x}.
#' @export
print.pvarife_result <- function(x, ...) {
  cat("Panel VAR with Interactive Fixed Effects\n")
  cat(sprintf("  Units (I): %d  |  Time periods (T): %d  |  Variables (K): %d\n",
              x$n_units, x$n_time, x$n_vars))
  cat(sprintf("  Lag order: %d  |  Factors: %d\n", x$n_lags, x$n_factors))
  cat("\nCoefficients:\n")
  n_vars <- x$n_vars
  n_lags <- x$n_lags
  beta   <- as.numeric(x$beta)

  intercepts <- beta[seq_len(n_vars)]
  names(intercepts) <- paste0("intercept[", seq_len(n_vars), "]")
  cat("  Intercepts:\n")
  print(round(intercepts, 6))

  lag_coefs <- beta[(n_vars + 1L):length(beta)]
  alpha_mat <- t(matrix(lag_coefs, nrow = n_vars * n_lags, ncol = n_vars))
  for (ll in seq_len(n_lags)) {
    col_s <- (ll - 1L) * n_vars + 1L
    col_e <- ll * n_vars
    cat(sprintf("\n  Theta_%d (VAR lag %d):\n", ll, ll))
    theta_l <- alpha_mat[, col_s:col_e, drop = FALSE]
    rownames(theta_l) <- paste0("y[", seq_len(n_vars), "]")
    colnames(theta_l) <- paste0("y[", seq_len(n_vars), ",t-", ll, "]")
    print(round(theta_l, 6))
  }

  cat("\nReduced-form covariance (Sigma):\n")
  print(round(x$sigma, 6))
  invisible(x)
}


#' Summary method for pvarife_result
#'
#' @param object An object of class \code{"pvarife_result"}.
#' @param ... Ignored.
#' @return Invisibly returns a list with \code{fit}, \code{avar}.
#' @export
summary.pvarife_result <- function(object, ...) {
  cat("=== Summary: pvarife ===\n\n")
  print(object)

  cat("\nComputing asymptotic variance ...\n")
  avar <- tryCatch(asymptotic_var(object),
                   error = function(e) {
                     cat("  (asymptotic_var failed:", conditionMessage(e), ")\n")
                     NULL
                   })

  if (!is.null(avar)) {
    se    <- sqrt(diag(avar$variance))
    beta  <- as.numeric(object$beta)
    z_val <- (beta - avar$bias) / se
    p_val <- 2.0 * stats::pnorm(-abs(z_val))

    coef_table <- data.frame(
      Estimate   = round(beta, 6),
      Std.Error  = round(se, 6),
      z.value    = round(z_val, 3),
      Pr.z       = round(p_val, 4),
      stringsAsFactors = FALSE
    )
    cat("\nInference (bias-corrected, Theorem 2.3):\n")
    print(coef_table)
  }

  invisible(list(fit = object, avar = avar))
}


#' Extract coefficients from a pvarife_result
#'
#' @param object An object of class \code{"pvarife_result"}.
#' @param ... Ignored.
#' @return Named numeric vector of coefficients.
#' @export
coef.pvarife_result <- function(object, ...) {
  beta <- as.numeric(object$beta)
  n_vars <- object$n_vars
  n_lags <- object$n_lags
  nms <- c(
    paste0("intercept[", seq_len(n_vars), "]"),
    unlist(lapply(seq_len(n_lags), function(ll) {
      paste0("theta", ll, "[", rep(seq_len(n_vars), each = n_vars), ",",
             rep(seq_len(n_vars), n_vars), "]")
    }))
  )
  stats::setNames(beta, nms)
}


#' Plot impulse response functions
#'
#' @param x An object of class \code{"pvarife_irf"}.
#' @param var_names Character vector of variable names (optional).
#' @param ... Ignored.
#' @return A \code{ggplot2} plot object.
#' @importFrom rlang .data
#' @export
plot.pvarife_irf <- function(x, var_names = NULL, ...) {
  ir     <- as.matrix(x)
  n_vars <- nrow(ir)
  n_periods <- ncol(ir)

  if (is.null(var_names)) var_names <- paste0("y[", seq_len(n_vars), "]")

  df <- data.frame(
    variable = rep(var_names, each = n_periods),
    horizon  = rep(seq_len(n_periods) - 1L, n_vars),
    irf      = as.numeric(t(ir)),
    stringsAsFactors = FALSE
  )
  df$variable <- factor(df$variable, levels = var_names)

  ggplot2::ggplot(df, ggplot2::aes(x = .data$horizon, y = .data$irf)) +
    ggplot2::geom_hline(yintercept = 0, colour = "black", linewidth = 0.4) +
    ggplot2::geom_line(linewidth = 0.9, colour = "black") +
    ggplot2::geom_point(size = 2, shape = 18, colour = "black") +
    ggplot2::facet_wrap(ggplot2::vars(.data$variable), scales = "free_y") +
    ggplot2::labs(x = "Horizon", y = "Response") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor  = ggplot2::element_blank(),
      strip.text        = ggplot2::element_text(face = "bold"),
      axis.title        = ggplot2::element_text(size = 10)
    )
}


#' Plot impulse response bands
#'
#' @param x An object of class \code{"pvarife_bands"}.
#' @param var_names Character vector of variable names (optional).
#' @param normalise_by_h1 Logical. If \code{TRUE}, divide all responses by the
#'   first-horizon response of the shock variable (replicates Figure 1 of
#'   Tugan 2021). Default \code{FALSE}.
#' @param ... Ignored.
#' @return A \code{ggplot2} plot object.
#' @export
plot.pvarife_bands <- function(x, var_names = NULL, normalise_by_h1 = FALSE,
                                ...) {
  n_vars    <- nrow(x$irf)
  n_periods <- ncol(x$irf)

  if (is.null(var_names)) var_names <- paste0("y[", seq_len(n_vars), "]")

  scale_factor <- if (normalise_by_h1) x$irf[1L, 1L] else 1.0

  df <- data.frame(
    variable = rep(var_names, each = n_periods),
    horizon  = rep(seq_len(n_periods) - 1L, n_vars),
    irf      = as.numeric(t(x$irf))   / scale_factor,
    lower    = as.numeric(t(x$lower)) / scale_factor,
    upper    = as.numeric(t(x$upper)) / scale_factor,
    stringsAsFactors = FALSE
  )
  df$variable <- factor(df$variable, levels = var_names)

  ci_label <- sprintf("%d%% CI (%s)", round(x$level * 100), x$method)

  ggplot2::ggplot(df, ggplot2::aes(x = .data$horizon)) +
    ggplot2::geom_hline(yintercept = 0, colour = "black", linewidth = 0.4) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
      fill = "grey75", alpha = 0.6
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = .data$irf),
      linewidth = 1.0, colour = "black"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(y = .data$irf),
      size = 2.5, shape = 18, colour = "black"
    ) +
    ggplot2::facet_wrap(ggplot2::vars(.data$variable), scales = "free_y") +
    ggplot2::labs(
      x       = "Horizon",
      y       = "Response",
      caption = ci_label
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor  = ggplot2::element_blank(),
      strip.text        = ggplot2::element_text(face = "bold"),
      axis.title        = ggplot2::element_text(size = 10),
      plot.caption      = ggplot2::element_text(size = 8, colour = "grey50")
    )
}
