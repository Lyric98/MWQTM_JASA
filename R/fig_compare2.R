# ============================================================================
# fig_compare2.R — comparison figures from results/compare2_fig_*.rds (5 DGPs).
# Cells looked up by (dgp,n) from the saved objects (order-independent).
#   Fig 1  beta_tau(tau): DGP2 (dispersion), DGP3 (non-proportional), DGP4 (tail).
#   Fig 2  mean log-contrast Delta(u,t): true vs Q vs Sun's flat (DGP3).
#   Fig 3  conditional-mean trajectories M(u,t): true vs Sun vs Q (DGP1).
#   Fig 4  DGP5 robustness: sampling spread of Sun's mean slope vs our median slope.
# Base graphics only.
# ============================================================================
.FIGS <- lapply(list.files("results", pattern = "^compare2_fig_[0-9]+\\.rds$",
                            full.names = TRUE), readRDS)
getfig <- function(dgp, n = 1000) {
  for (o in .FIGS) if (o$dgp == dgp && o$n == n) return(o)
  stop(sprintf("no fig cell for dgp=%d n=%d", dgp, n))
}
COL <- list(true = "black", q = "#1f6feb", sun = "#d1242f")
dir.create("paper/figs", showWarnings = FALSE, recursive = TRUE)

# ---- Fig 1: beta_tau curves (DGP2, DGP3, DGP4) -----------------------------
draw_bt <- function(f, ttl) {
  tau <- f$taus; lo <- f$qb_curve_mean - f$qb_curve_sd; hi <- f$qb_curve_mean + f$qb_curve_sd
  yl <- range(c(f$true_betatau, lo, hi, f$beta_sun_mean), na.rm = TRUE) + c(-.05, .05)
  plot(tau, f$true_betatau, type = "n", ylim = yl, xlab = expression(tau),
       ylab = expression(beta[list(tau, 1)]), main = ttl)
  polygon(c(tau, rev(tau)), c(lo, rev(hi)), col = adjustcolor(COL$q, .18), border = NA)
  lines(tau, f$true_betatau, col = COL$true, lwd = 2.4)
  lines(tau, f$qb_curve_mean, col = COL$q, lwd = 2.2, lty = 2)
  points(tau, f$qb_curve_mean, col = COL$q, pch = 19, cex = .5)
  abline(h = f$beta_sun_mean, col = COL$sun, lwd = 2.2, lty = 4)
  legend("topleft", bty = "n", cex = .8, lwd = 2.2, lty = c(1, 2, 4),
         col = c(COL$true, COL$q, COL$sun),
         legend = c("truth", "quantile (Q)", expression(paste("Sun ", hat(beta)[1]))))
}
pdf("paper/figs/cmp_betatau.pdf", width = 11.4, height = 3.9)
par(mfrow = c(1, 3), mar = c(4.0, 4.0, 2.4, 0.8), cex.lab = 1.1)
draw_bt(getfig(2), "DGP2: dispersion varies")
draw_bt(getfig(3), "DGP3: non-proportional mean")
draw_bt(getfig(4), "DGP4: tail-specific effect")
dev.off(); cat("wrote paper/figs/cmp_betatau.pdf\n")

# ---- Fig 2: mean log-contrast Delta(u,t) (DGP3) ----------------------------
f <- getfig(3); tg <- sort(unique(f$pts$t)); sel <- tg[c(2, 4, 6)]; cols <- c("#444444", "#1f6feb", "#d1242f")
pdf("paper/figs/cmp_contrast.pdf", width = 6.4, height = 5.0)
par(mar = c(4.0, 4.2, 2.4, 0.8), cex.lab = 1.05)
rng <- range(c(f$trueD, f$qD, f$beta_sun_mean), na.rm = TRUE)
plot(NA, xlim = range(f$pts$u), ylim = rng + c(-.05, .05),
     xlab = "visit age u", ylab = expression(paste("mean log-contrast  ", Delta(u, t))),
     main = "DGP3: induced mean contrast varies with (u,t)")
for (j in seq_along(sel)) {
  ii <- which(abs(f$pts$t - sel[j]) < 1e-6); ii <- ii[order(f$pts$u[ii])]
  lines(f$pts$u[ii], f$trueD[ii], col = cols[j], lwd = 2.4)
  lines(f$pts$u[ii], f$qD[ii], col = cols[j], lwd = 2.0, lty = 2); points(f$pts$u[ii], f$qD[ii], col = cols[j], pch = 19, cex = .6)
  if (!is.null(f$ssD)) lines(f$pts$u[ii], f$ssD[ii], col = cols[j], lwd = 1.4, lty = 3)  # Sun-stratified (noisy)
}
abline(h = f$beta_sun_mean, col = "black", lwd = 2.0, lty = 4)
legend("topright", bty = "n", cex = .78,
       legend = c(sprintf("t=%.2f", sel), "truth (solid)", "Q (dashed)", "Sun-strat (dotted)", expression(paste("Sun-basic ", hat(beta)[1]))),
       col = c(cols, "black", "black", "black", "black"), lwd = c(2.4, 2.4, 2.4, 2.4, 2.0, 1.4, 2.0),
       lty = c(1, 1, 1, 1, 2, 3, 4), pch = c(NA, NA, NA, NA, 19, NA, NA))
dev.off(); cat("wrote paper/figs/cmp_contrast.pdf\n")

# ---- Fig 3: conditional-mean trajectories (DGP1) ---------------------------
f <- getfig(1); tg <- sort(unique(f$pts$t)); sel <- tg[c(2, 4, 6)]; cols <- c("#444444", "#1f6feb", "#d1242f")
pdf("paper/figs/cmp_meansurface.pdf", width = 6.4, height = 5.0)
par(mar = c(4.0, 4.2, 2.4, 0.8), cex.lab = 1.05)
logT <- log(f$trueM0); rng <- range(c(logT, f$qLM0, f$sLM0), na.rm = TRUE)
plot(NA, xlim = range(f$pts$u), ylim = rng + c(-.05, .05),
     xlab = "visit age u", ylab = expression(paste("log conditional mean  ", log, " E[Y|u,t,x]")),
     main = "DGP1: conditional-mean recovery (proportional mean)")
for (j in seq_along(sel)) {
  ii <- which(abs(f$pts$t - sel[j]) < 1e-6); ii <- ii[order(f$pts$u[ii])]
  lines(f$pts$u[ii], logT[ii], col = cols[j], lwd = 2.4)
  lines(f$pts$u[ii], f$qLM0[ii], col = cols[j], lwd = 2.0, lty = 2); points(f$pts$u[ii], f$qLM0[ii], col = cols[j], pch = 19, cex = .55)
  lines(f$pts$u[ii], f$sLM0[ii], col = cols[j], lwd = 1.8, lty = 3); points(f$pts$u[ii], f$sLM0[ii], col = cols[j], pch = 4, cex = .6)
}
legend("topleft", bty = "n", cex = .82,
       legend = c(sprintf("t=%.2f", sel), "truth (solid)", "Q (dashed)", "Sun (dotted)"),
       col = c(cols, "black", "black", "black"), lwd = c(2.4, 2.4, 2.4, 2.4, 2.0, 1.8),
       lty = c(1, 1, 1, 1, 2, 3), pch = c(NA, NA, NA, NA, 19, 4))
dev.off(); cat("wrote paper/figs/cmp_meansurface.pdf\n")

# ---- Fig 4: DGP5 robustness (Sun's mean slope vs our median slope) ----------
ns <- c(500, 1000, 2000)
box_data <- list(); labs <- character(0); ats <- numeric(0); bcol <- character(0); k <- 1
for (j in seq_along(ns)) {
  f <- getfig(5, ns[j])
  box_data[[k]] <- f$beta_sun1_reps; ats[k] <- j * 3 - 1; bcol[k] <- COL$sun; labs[k] <- ""; k <- k + 1
  box_data[[k]] <- f$qb50_reps;      ats[k] <- j * 3;     bcol[k] <- COL$q;   labs[k] <- ""; k <- k + 1
}
pdf("paper/figs/cmp_heavytail.pdf", width = 6.6, height = 4.6)
par(mar = c(4.2, 4.2, 2.6, 0.8), cex.lab = 1.05)
boxplot(box_data, at = ats, col = adjustcolor(bcol, .35), border = bcol, outline = TRUE,
        xaxt = "n", ylab = expression(paste("estimated slope (truth ", beta[1], "=0.5)")),
        main = expression(paste("DGP5 (heavy ", t[5], " tail): Sun mean slope vs. our ", beta[0.5])),
        pch = 20, cex = .4)
abline(h = 0.5, lty = 2, col = "grey40")
axis(1, at = (seq_along(ns)) * 3 - 0.5, labels = paste0("n=", ns))
legend("topright", bty = "n", fill = adjustcolor(c(COL$sun, COL$q), .35),
       border = c(COL$sun, COL$q), legend = c(expression(paste("Sun ", hat(beta)[1])),
       expression(paste("Q ", hat(beta)[0.5]))), cex = .9)
dev.off(); cat("wrote paper/figs/cmp_heavytail.pdf\n")
