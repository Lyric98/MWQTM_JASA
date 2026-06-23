# plot_surface.R — non-parallel quantile-trajectory recovery figure.
# True vs mean-estimated m_tau(u,t) slices at fixed t, for tau in {0.1,0.5,0.9},
# het DGP, n=2000. Shows the surfaces are non-parallel and are recovered.
g <- read.csv("results/surface_grid.csv")
g <- g[g$paramset == "het" & g$n == 2000, ]
tslices <- c(3, 5)
cols <- c("0.25"="#1b9e77","0.5"="#7570b3","0.75"="#d95f02")
pdf("results/surface_trajectory.pdf", width = 9, height = 4.2)
op <- par(mfrow = c(1, length(tslices)), mar = c(4, 4, 2.5, 1))
for (tt in tslices) {
  gg <- g[abs(g$t - tt) == min(abs(g$t - tt)), ]
  plot(NA, xlim = range(gg$u), ylim = range(c(gg$m_true, gg$m_hat)),
       xlab = "age u", ylab = expression(m[tau](u, t)),
       main = sprintf("t = %.1f", gg$t[1]))
  for (q in c(0.25, 0.5, 0.75)) {
    z <- gg[abs(gg$tau - q) < 1e-6, ]; z <- z[order(z$u), ]
    lines(z$u, z$m_true, col = cols[as.character(q)], lwd = 2)              # truth
    points(z$u, z$m_hat, col = cols[as.character(q)], pch = 1, cex = 0.8)   # estimate
  }
  if (tt == tslices[1])
    legend("topleft", bty = "n", cex = 0.85,
           legend = c(expression(tau == 0.25), expression(tau == 0.5), expression(tau == 0.75),
                      "true (line)", "est. (points)"),
           col = c(cols, "black", "black"), lwd = c(2,2,2,NA,NA), pch = c(NA,NA,NA,NA,1))
}
par(op); dev.off()
cat("wrote results/surface_trajectory.pdf\n")
