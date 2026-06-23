# Visual confirmation that chi^C = 0: the estimated-vs-oracle weight gap
# sd(beta_ipcw - beta_orc) decays like n^{-1} (second order), faster than the
# estimator SD esd ~ n^{-1/2}. Log-log with reference slopes.
r <- read.csv("results/compare.csv")
pdf("results/compare_rate.pdf", width = 10, height = 4.6)
op <- par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1))
for (sc in c(1, 2)) {
  s <- r[r$scenario == sc & r$tau == 0.50, ]; s <- s[order(s$n), ]
  plot(s$n, s$esd_ipcw, log = "xy", type = "b", pch = 16, col = "navy",
       ylim = range(c(s$esd_ipcw, s$sd_gap)), lwd = 2,
       xlab = "n", ylab = "dispersion (log scale)",
       main = sprintf("Scenario %d, tau=0.5", sc))
  lines(s$n, s$sd_gap, type = "b", pch = 17, col = "firebrick", lwd = 2)
  # reference slopes anchored at first point
  n0 <- s$n[1]
  lines(s$n, s$esd_ipcw[1] * (s$n / n0)^(-0.5), lty = 2, col = "navy")
  lines(s$n, s$sd_gap[1]  * (s$n / n0)^(-1.0), lty = 2, col = "firebrick")
  legend("bottomleft", bty = "n", cex = 0.85,
         legend = c("ESD(beta_ipcw)   ~ n^-1/2", "sd(beta_ipcw - beta_orc)  ~ n^-1",
                    "ref n^-1/2", "ref n^-1"),
         col = c("navy", "firebrick", "navy", "firebrick"),
         pch = c(16, 17, NA, NA), lty = c(1, 1, 2, 2), lwd = 2)
}
par(op); dev.off()

# fitted log-log slopes
cat("=== log-log slopes (sd_gap should be ~ -1, esd ~ -0.5) ===\n")
for (sc in c(1, 2)) for (tt in c(0.1, 0.5, 0.9)) {
  s <- r[r$scenario == sc & r$tau == tt, ]
  sl_gap <- coef(lm(log(sd_gap) ~ log(n), s))[2]
  sl_esd <- coef(lm(log(esd_ipcw) ~ log(n), s))[2]
  cat(sprintf("scen %d tau %.1f : slope(sd_gap)=%.2f  slope(esd)=%.2f\n",
              sc, tt, sl_gap, sl_esd))
}
cat("=== est-IPCW vs oracle-IPCW ESD ratio (should be ~1) ===\n")
cat(sprintf("mean esd_ipcw/esd_orc = %.3f ; mean esd_cc/esd_ipcw = %.3f\n",
            mean(r$esd_ipcw / r$esd_orc), mean(r$esd_cc / r$esd_ipcw)))
