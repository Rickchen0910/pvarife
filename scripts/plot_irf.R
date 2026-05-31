#!/usr/bin/env Rscript
# =============================================================================
# Plot IRF figures from mc_irf.R CSV output  (run LOCALLY, not on HPC)
#
# Run from repo root:
#   Rscript scripts/plot_irf.R sr    # Figure S.1 — short-run
#   Rscript scripts/plot_irf.R lr    # Figure S.2 — long-run
#
# Input:  scripts/csv/mc_irf_{sr|lr}_I*_T*.csv
# Output: scripts/figures/fig_s1_sr.pdf  (or fig_s2_lr.pdf)
# =============================================================================

# ---- load pvarife locally ---------------------------------------------------
devtools::load_all(quiet = TRUE)

library(ggplot2)

# ---- command-line argument --------------------------------------------------
args   <- commandArgs(trailingOnly = TRUE)
id_str <- if (length(args) >= 1L && args[1L] == "lr") "lr" else "sr"

identification <- if (id_str == "sr") "short_run" else "long_run"
id_label       <- if (id_str == "sr") "Short-Run"  else "Long-Run"
diff_vars      <- if (id_str == "lr") 1L else integer(0L)

cat(sprintf("Making Figure for %s identification (%s)\n\n", id_label, id_str))

# ---- constants --------------------------------------------------------------
horizons  <- 0L:20L
n_periods <- 21L
shock_idx <- 1L

# True DGP parameters
# Our beta ordering: [intercepts, Theta11, Theta12, Theta21, Theta22]
true_beta  <- c(1.0, 1.0, 0.65, 0.30, 0.20, 0.60)
true_sigma <- matrix(c(1.0, 0.5, 0.5, 1.0), nrow = 2L)

# Compute true IRF (known parameters, deterministic)
make_fit <- function(beta, sigma) {
  structure(list(beta = beta, sigma = sigma, n_vars = 2L, n_lags = 1L),
            class = "pvarife_result")
}
true_raw <- compute_irf(make_fit(true_beta, true_sigma),
                         n_periods      = n_periods,
                         shock          = shock_idx,
                         diff_vars      = diff_vars,
                         identification = identification)
true_ir <- true_raw / true_raw[shock_idx, 1L]   # normalise

# ---- find CSV files ---------------------------------------------------------
# Auto-find the csv/ folder — works regardless of working directory
find_csv_dir <- function() {
  # Try several locations in priority order
  candidates <- c(
    file.path("scripts", "csv"),                      # from repo root
    "csv",                                            # from scripts/
    file.path(dirname(sys.frame(1)$ofile), "csv")    # from script location
  )
  for (d in candidates) {
    if (dir.exists(d)) return(d)
  }
  # Fall back: search upward for DESCRIPTION (repo root marker)
  wd <- getwd()
  for (i in 1:5) {
    if (file.exists(file.path(wd, "DESCRIPTION"))) {
      d <- file.path(wd, "scripts", "csv")
      if (dir.exists(d)) return(d)
    }
    wd <- dirname(wd)
  }
  file.path("scripts", "csv")   # last resort
}

csv_dir   <- tryCatch(find_csv_dir(), error = function(e) file.path("scripts","csv"))
cat("Looking for CSVs in:", normalizePath(csv_dir, mustWork=FALSE), "\n")

csv_files <- sort(list.files(
  csv_dir,
  pattern    = paste0("^mc_irf_", id_str, "_I[0-9]+_T[0-9]+\\.csv$"),
  full.names = TRUE
))

if (length(csv_files) == 0L) {
  # Show what IS in the directory to help diagnose
  all_files <- list.files(csv_dir)
  cat("Files found in csv/:\n")
  if (length(all_files) == 0L) cat("  (empty)\n") else cat(paste(" ", all_files, "\n"))
  stop("\nNo CSV files matching pattern 'mc_irf_", id_str,
       "_I<n>_T<t>.csv' found in:\n  ", normalizePath(csv_dir, mustWork=FALSE),
       "\nDownload the CSV files from the HPC first.")
}

cat("CSV files found:\n")
cat(paste(" -", basename(csv_files), "\n"))
cat("\n")

# ---- process one CSV file ---------------------------------------------------
process_csv <- function(path) {
  df <- read.csv(path, check.names = FALSE)

  # Parse I and T from filename
  fname   <- basename(path)
  n_units <- as.integer(regmatches(fname, regexpr("(?<=I)\\d+", fname, perl=TRUE)))
  n_time  <- as.integer(regmatches(fname, regexpr("(?<=T)\\d+(?=\\.)", fname, perl=TRUE)))
  # n_valid: count rows where h=0 IRF is not NA (detection happens below)
  # Use the first irf_1_h* column to count valid rows
  h0_col_guess <- grep("^irf_1_h", names(df), value=TRUE)
  h0_col_guess <- if (length(h0_col_guess) > 0L)
    h0_col_guess[which.min(nchar(h0_col_guess))] else "irf_1_h_0"
  n_valid <- sum(!is.na(df[[h0_col_guess]]))
  n_total <- nrow(df)

  cat(sprintf("I=%d, T=%d: %d/%d valid replications\n", n_units, n_time, n_valid, n_total))

  # Auto-detect column naming by building a direct horizon->column-suffix map.
  # Works with any format mc_irf.R may have produced:
  #   "irf_1_h00"  (sprintf "%02d")
  #   "irf_1_h_0"  (paste0 "_", h)
  #   "irf_1_h.0"  (space-padded, read by default read.csv)
  all_cols  <- names(df)
  irf1_cols <- grep("^irf_1_h", all_cols, value = TRUE)
  if (length(irf1_cols) == 0L)
    stop("No 'irf_1_h*' columns found in ", path,
         "\nActual columns: ", paste(head(all_cols, 10L), collapse=", "))

  # Build a map: integer horizon -> actual column suffix
  # Strategy: for each candidate h (0..20), find which irf1_col ends with that h
  suffix_map <- setNames(vector("character", 21L), as.character(0:20))
  for (col in irf1_cols) {
    sfx <- sub("^irf_1_h", "", col)          # e.g. "_0", "00", ".0", "0"
    # numeric value of this suffix
    h_val <- suppressWarnings(as.integer(gsub("[^0-9]", "", sfx)))
    if (!is.na(h_val) && h_val >= 0L && h_val <= 20L)
      suffix_map[as.character(h_val)] <- sfx
  }
  if (any(suffix_map == ""))
    warning("Could not detect suffix for all horizons 0-20 in ", path)

  h_str_fn <- function(h) suffix_map[as.character(h)]

  # Verify
  h1_col <- paste0("irf_1_h", h_str_fn(1L))
  if (!h1_col %in% all_cols)
    warning("Inferred column '", h1_col, "' not found in CSV.")

  # Helper: column mean across replications
  col_avg <- function(nm) mean(df[[nm]], na.rm = TRUE)

  # Average pvarife IRF and CI bands (already normalised per-replication)
  irf <- matrix(NA_real_, 2L, 21L)
  lo  <- matrix(NA_real_, 2L, 21L)
  hi  <- matrix(NA_real_, 2L, 21L)
  for (vv in 1:2) {
    for (hh in horizons) {
      h_str          <- h_str_fn(hh)
      irf[vv, hh+1L] <- col_avg(paste0("irf_", vv, "_h", h_str))
      lo [vv, hh+1L] <- col_avg(paste0("lo_",  vv, "_h", h_str))
      hi [vv, hh+1L] <- col_avg(paste0("hi_",  vv, "_h", h_str))
    }
  }

  # Average OLS beta and Sigma → OLS IRF (from averages, matching MATLAB)
  avg_b <- sapply(paste0("ols_b_", 1:6), col_avg)
  avg_s <- matrix(c(col_avg("ols_s_11"), col_avg("ols_s_12"),
                    col_avg("ols_s_12"), col_avg("ols_s_22")), nrow = 2L)
  ols_raw <- tryCatch(
    compute_irf(make_fit(avg_b, avg_s),
                n_periods      = n_periods,
                shock          = shock_idx,
                diff_vars      = diff_vars,
                identification = identification),
    error = function(e) { warning("OLS IRF failed: ", e$message); NULL })
  ols_ir <- if (!is.null(ols_raw))
    ols_raw / ols_raw[shock_idx, 1L]
  else
    matrix(NA_real_, 2L, 21L)

  list(n_units = n_units, n_time  = n_time, n_valid = n_valid,
       irf = irf, lo = lo, hi = hi, ols_ir = ols_ir)
}

results <- lapply(csv_files, process_csv)

# Sort by (n_units ASC, n_time ASC)
ord     <- order(sapply(results, `[[`, "n_units"),
                 sapply(results, `[[`, "n_time"))
results <- results[ord]

# ---- build tidy data for ggplot2 --------------------------------------------
panel_levels <- vapply(results, function(r)
  sprintf("I=%d, T=%d", r$n_units, r$n_time), character(1L))

var_labels <- c("1" = "Variable 1", "2" = "Variable 2")

make_line_df <- function(res, tag) {
  do.call(rbind, lapply(1:2, function(vv) {
    rbind(
      data.frame(panel=tag, var=vv, horizon=horizons,
                 value=true_ir[vv, ],     type="True IRF"),
      data.frame(panel=tag, var=vv, horizon=horizons,
                 value=res$irf[vv, ],    type="pvarife (avg)"),
      data.frame(panel=tag, var=vv, horizon=horizons,
                 value=res$ols_ir[vv, ], type="OLS (avg)")
    )
  }))
}
make_band_df <- function(res, tag) {
  do.call(rbind, lapply(1:2, function(vv)
    data.frame(panel=tag, var=vv, horizon=horizons,
               lo=res$lo[vv, ], hi=res$hi[vv, ])))
}

df_lines <- do.call(rbind, Map(make_line_df, results, panel_levels))
df_bands <- do.call(rbind, Map(make_band_df, results, panel_levels))

df_lines$panel <- factor(df_lines$panel, levels = panel_levels)
df_bands$panel <- factor(df_bands$panel, levels = panel_levels)
df_lines$type  <- factor(df_lines$type,
                          levels = c("True IRF","pvarife (avg)","OLS (avg)"))
df_lines$var_f <- factor(df_lines$var, levels=1:2, labels=var_labels)
df_bands$var_f <- factor(df_bands$var, levels=1:2, labels=var_labels)

# ---- plot -------------------------------------------------------------------
type_colours   <- c("True IRF"="black", "pvarife (avg)"="black",   "OLS (avg)"="grey40")
type_linetypes <- c("True IRF"="solid", "pvarife (avg)"="dashed",  "OLS (avg)"="dashed")
type_linewidths<- c("True IRF"=0.85,    "pvarife (avg)"=0.85,      "OLS (avg)"=0.65)
type_shapes    <- c("True IRF"=8L,      "pvarife (avg)"=3L,        "OLS (avg)"=NA)

p <- ggplot(df_bands, aes(x = horizon)) +
  # --- confidence band (avg of per-rep 95% CI) ---
  geom_ribbon(aes(ymin = lo, ymax = hi),
              fill = "grey80", alpha = 0.6) +
  # --- zero line ---
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "black") +
  # --- IRF lines ---
  geom_line(data  = df_lines,
            aes(x = horizon, y = value,
                colour = type, linetype = type, linewidth = type)) +
  # --- markers (True IRF: asterisk; pvarife: plus) ---
  geom_point(
    data  = df_lines[df_lines$type != "OLS (avg)", ],
    aes(x = horizon, y = value, colour = type, shape = type),
    size  = 1.6
  ) +
  # --- scales ---
  scale_colour_manual(  values = type_colours,    name = NULL) +
  scale_linetype_manual(values = type_linetypes,  name = NULL) +
  scale_linewidth_manual(values= type_linewidths, name = NULL) +
  scale_shape_manual(   values = type_shapes,     name = NULL) +
  scale_x_continuous(breaks = seq(0L, 20L, 4L)) +
  # --- facet: rows = variable, cols = (I, T) panel ---
  facet_grid(var_f ~ panel, scales = "free_y") +
  labs(
    title    = paste0("Figure S.", if(id_str=="sr") "1" else "2",
                      " -- ", id_label, " Identification (CommonFactors DGP)"),
    subtitle = paste0(
      "Shaded: avg 95% CI (pvarife)  |  ",
      "* True IRF  |  --+ pvarife (avg)  |  -- OLS (avg)  |  ",
      "Normalised by shock-1 contemporaneous response"),
    x = "Period",
    y = "Normalised impulse response"
  ) +
  theme_minimal(base_size = 10L) +
  theme(
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(colour = "grey85", fill = NA, linewidth = 0.4),
    strip.text        = element_text(face = "bold", size = 9L),
    strip.background  = element_rect(fill = "grey95"),
    legend.position   = "bottom",
    legend.key.width  = unit(1.8, "cm"),
    plot.title        = element_text(size = 11L),
    plot.subtitle     = element_text(size = 7.5, colour = "grey45")
  ) +
  guides(linewidth = "none")

# Dimensions scale with number of panels
n_panels <- length(results)
fig_w    <- max(6.5, 2.8 * n_panels + 1.2)
fig_h    <- 5.2

out_dir <- file.path("scripts", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_pdf <- file.path(out_dir, sprintf("fig_s%s_%s.pdf",
                                       if (id_str=="sr") "1" else "2", id_str))
ggsave(out_pdf, p, width = fig_w, height = fig_h, device = "pdf")
cat(sprintf("\nFigure saved: %s  (%.1f x %.1f in)\n", out_pdf, fig_w, fig_h))

# ---- quick sanity check table -----------------------------------------------
cat("\nSanity check (h=0 should be 1.00 for var1, ~0.50 for var2):\n")
hdr <- sprintf("%-14s %8s %8s %8s %8s %8s %8s",
               "Panel","Tr_v1","pv_v1","ols_v1","Tr_v2","pv_v2","ols_v2")
cat(hdr, "\n", strrep("-", nchar(hdr)), "\n", sep="")
for (ii in seq_along(results)) {
  r <- results[[ii]]
  cat(sprintf("%-14s %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
              panel_levels[ii],
              true_ir[1,1], r$irf[1,1], r$ols_ir[1,1],
              true_ir[2,1], r$irf[2,1], r$ols_ir[2,1]))
}
