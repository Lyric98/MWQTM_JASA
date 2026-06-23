# ============================================================================
# plot_coverage.R — the decisive coverage-vs-tau figure + decision table.
# Reads results/pilot.csv (from run_pilot.R) and renders one panel per (n,scen).
# ============================================================================
args <- commandArgs(trailingOnly = TRUE)
infile <- if (length(args) >= 1) args[1] else "results/pilot.csv"
res <- read.csv(infile)

mc_band <- function(R = 150, p = 0.95) 1.96 * sqrt(p * (1 - p) / R)

panels <- unique(res[, c("n", "scenario")])
panels <- panels[order(panels$scenario, panels$n), ]
pdf("results/coverage_vs_tau.pdf", width = 11, height = 9)
op <- par(mfrow = c(2, length(unique(res$n))), mar = c(4, 4, 2.5, 1))
for (i in seq_len(nrow(panels))) {
  s <- res[res$n == panels$n[i] & res$scenario == panels$scenario[i], ]
  s <- s[order(s$tau), ]
  plot(NA, xlim = range(res$tau), ylim = c(0.7, 1.0),
       xlab = expression(tau), ylab = "95% coverage of beta1",
       main = sprintf("n=%d, Scenario %d", panels$n[i], panels$scenario[i]))
  abline(h = 0.95, col = "grey50", lty = 2)
  R <- 200
  rect(par("usr")[1], 0.95 - mc_band(R), par("usr")[2], 0.95 + mc_band(R),
       col = adjustcolor("grey80", 0.4), border = NA)
  cols <- c(o3 = "black", bos_sw = "firebrick", bos_b = "red",
            acv_b = "blue", aus_b = "darkgreen")
  lines(s$tau, s$cov_o3_sw, col = cols["o3"], lwd = 2, type = "b", pch = 16)
  lines(s$tau, s$cov_bos_sw, col = cols["bos_sw"], lwd = 2, type = "b", pch = 17)
  if (any(!is.na(s$cov_bos_boot)))
    lines(s$tau, s$cov_bos_boot, col = cols["bos_b"], lwd = 2, type = "b", pch = 1, lty = 3)
  if (any(!is.na(s$cov_acv_boot)))
    lines(s$tau, s$cov_acv_boot, col = cols["acv_b"], lwd = 2, type = "b", pch = 1, lty = 3)
  if (any(!is.na(s$cov_aus_boot)))
    lines(s$tau, s$cov_aus_boot, col = cols["aus_b"], lwd = 2, type = "b", pch = 1, lty = 3)
  if ("cov_bos_froz" %in% names(s) && any(!is.na(s$cov_bos_froz)))
    lines(s$tau, s$cov_bos_froz, col = "orange", lwd = 1.5, type = "b", pch = 4, lty = 4)
  # condition number annotation at extreme tau
  text(s$tau, rep(0.72, nrow(s)), sprintf("k=%.0f", s$Scond_med), cex = 0.7, col = "grey30")
  if (i == 1)
    legend("bottomright", bty = "n", cex = 0.8,
           legend = c("IPCW oracle-nuis. sandwich",
                      "IPCW-B-OS-NCF sandwich (no chi^C)", "IPCW-B-OS-NCF boot",
                      "IPCW-A-CV boot", "IPCW-A-Rich boot",
                      "IPCW-B-OS-NCF frozen-G_C boot"),
           col = c(cols, "orange"), lwd = 2,
           pch = c(16, 17, 1, 1, 1, 4), lty = c(1, 1, 3, 3, 3, 4))
}
par(op); dev.off()

# ---- bias / SER / coverage summary table ----------------------------------
cat("\n==== PILOT SUMMARY ====\n")
fmt <- res[, c("n", "scenario", "tau", "true_b1",
               "bias_acv", "bias_aus", "bias_bos", "bias_o3",
               "ser_bos", "cov_bos_sw", "cov_o3_sw",
               "cov_acv_boot", "cov_aus_boot", "cov_bos_boot",
               "Scond_med", "df_med")]
print(round(fmt, 3), row.names = FALSE)

# ---- decision-table verdict (reply3 section 13) ---------------------------
verdict <- function(r) {
  # "ok" = nominal 0.95 lies within the estimator's Monte-Carlo band, i.e.
  # coverage is not significantly below nominal (over-coverage is acceptable).
  ok <- function(c) !is.na(c) && c >= 0.95 - mc_band(40)
  acv <- ok(r$cov_acv_boot); aus <- ok(r$cov_aus_boot)
  bos <- ok(r$cov_bos_boot); o3 <- ok(r$cov_o3_sw)
  if (isTRUE(acv)) "A-CV ok: implicit orthogonalization suffices"
  else if (isTRUE(aus)) "A-CV fails, A-US ok: penalty/sieve bias -> undersmoothing"
  else if (isTRUE(bos)) "A-US fails, B-OS ok: explicit orthogonalization needed"
  else if (isTRUE(o3)) "B-OS fails, O3 ok: f/g nuisance estimation is the bottleneck"
  else "O3 fails too: revisit estimand / sieve / implementation"
}
cat("\n==== PER-CELL VERDICT ====\n")
for (i in seq_len(nrow(res)))
  cat(sprintf("n=%d scen=%d tau=%.2f : %s\n",
              res$n[i], res$scenario[i], res$tau[i], verdict(res[i, ])))
