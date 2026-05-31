#!/usr/bin/env Rscript
# =============================================================================
# Merge results from SLURM array tasks
# Run AFTER all 9 tasks complete:
#   Rscript mc_merge_results.R
# =============================================================================

library(ggplot2)

out_dir <- "mc_hpc_results"
task_files <- list.files(out_dir, pattern = "^mc_task[0-9]+\\.rds$", full.names = TRUE)
cat("Found", length(task_files), "task files\n")
if (length(task_files) == 0L) stop("No task files in ", out_dir)

res <- do.call(rbind, lapply(task_files, readRDS))
saveRDS(res, file.path(out_dir, "mc_results_merged.rds"))
cat("Merged →", file.path(out_dir, "mc_results_merged.rds"), "\n\n")

theta_params <- paste0("beta[", 3:6, "]")
theta_names  <- c("Theta11=0.65","Theta21=0.20","Theta12=0.30","Theta22=0.60")
res_theta    <- res[res$param %in% theta_params, ]

cat("=== Bias ===\n")
for (pp in theta_params) {
  sub <- res_theta[res_theta$param == pp, ]
  tbl <- reshape(sub[, c("n_units","n_time","bias")],
                 idvar="n_units", timevar="n_time", direction="wide")
  cat("\n", pp, "\n"); print(round(tbl, 4))
}

cat("\n=== RMSE ===\n")
for (pp in theta_params) {
  sub <- res_theta[res_theta$param == pp, ]
  tbl <- reshape(sub[, c("n_units","n_time","rmse")],
                 idvar="n_units", timevar="n_time", direction="wide")
  cat("\n", pp, "\n"); print(round(tbl, 4))
}

cat("\n=== Coverage (target 0.950) ===\n")
for (pp in theta_params) {
  sub <- res_theta[res_theta$param == pp, ]
  tbl <- reshape(sub[, c("n_units","n_time","coverage")],
                 idvar="n_units", timevar="n_time", direction="wide")
  cat("\n", pp, "\n"); print(round(tbl, 3))
}

# Plot
res_plot <- res_theta
res_plot$T_label <- factor(res_plot$n_time, levels=c(25,50,100),
                            labels=c("T=25","T=50","T=100"))
res_plot$I_label <- factor(res_plot$n_units)

p <- ggplot(res_plot, aes(x=I_label, y=rmse, colour=T_label, group=T_label)) +
  geom_line(linewidth=0.8) + geom_point(size=2.5) +
  facet_wrap(~ param, scales="free_y", ncol=2) +
  labs(title="RMSE convergence — pvarife (Tugan 2021)",
       x="Units (I)", y="RMSE", colour="Time periods") +
  theme_minimal(base_size=11) +
  theme(panel.grid.minor=element_blank(),
        strip.text=element_text(face="bold"), legend.position="bottom")
ggsave(file.path(out_dir, "mc_convergence.pdf"), p, width=8, height=5)
cat("\nPlot:", file.path(out_dir, "mc_convergence.pdf"), "\n")
