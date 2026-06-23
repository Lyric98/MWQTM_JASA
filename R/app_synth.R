# ============================================================================
# app_synth.R -- ILLUSTRATIVE application figure on a SYNTHETIC ADNI-matched
# cohort (NOT real ADNI; demonstrates the deliverables only). Reuses the
# Sun-style survival/truncation/visit design and our quantile estimator.
# Produces paper/figs/app_panels.pdf with three panels:
#   (L) median trajectory surface  m_.5(u,t)  at three onset ages
#   (C) APOE4 carrier effect  beta_{tau,1}  across tau, with bootstrap band
#   (R) induced mean log-contrast (quantile) vs Sun's scalar mean coefficient
# ============================================================================
suppressWarnings(suppressMessages({ source("R/compare2.R"); try_compile_sun_cpp() }))
set.seed(2026)

p <- cmp2_params()
n <- 1500
d <- gen_sun(n, 3, p)            # DGP3-style: distributional (non-proportional) effect; X1 := APOE4 carrier
mk <- cc_rows2(d)

taus  <- c(.05,.1,.15,.2,.25,.3,.4,.5,.6,.7,.75,.8,.85,.9,.95)
qf    <- qfit2(mk, taus, df = 4)
qb    <- qbeta2(qf)[, 1]          # APOE4 effect across tau
sb    <- sun_beta2(d)            # Sun scalar mean slope

# evaluation grid (interior)
ug <- seq(0.35, 0.85, length.out = 30)
tsel <- c(0.55, 0.70, 0.90)
cols <- c("#444444", "#1f6feb", "#d1242f")

# --- subject bootstrap band for beta_tau1 (B resamples) ---
B <- 120; ids <- unique(mk$id)
bootb <- matrix(NA_real_, B, length(taus))
for (b in seq_len(B)) {
  samp <- sample(ids, length(ids), replace = TRUE)
  rows <- do.call(rbind, lapply(seq_along(samp), function(j) {
    r <- mk[mk$id == samp[j], ]; if (nrow(r)) { r$id <- j; r } }))
  qfb <- tryCatch(qfit2(rows, taus, df = 4), error = function(e) NULL)
  if (!is.null(qfb)) bootb[b, ] <- qbeta2(qfb)[, 1]
}
blo <- apply(bootb, 2, quantile, .025, na.rm = TRUE)
bhi <- apply(bootb, 2, quantile, .975, na.rm = TRUE)

# --- induced mean contrast (Q-LN) and Sun scalar, over u at mid onset t ---
tmid <- 0.70
qD <- qmean_ln(qf, ug, rep(tmid, length(ug)), 1, 0) - qmean_ln(qf, ug, rep(tmid, length(ug)), 0, 0)

pdf("paper/figs/app_panels.pdf", width = 12.4, height = 3.9)
par(mfrow = c(1, 3), mar = c(4.1, 4.3, 2.6, 0.8), cex.lab = 1.05, cex.main = 1.05)

## (L) median trajectory surface m_.5(u,t) at three onset ages
plot(NA, xlim = range(ug), ylim = c(-0.7, 0.7), xlab = "visit age u (centered)",
     ylab = expression(paste("median ", widehat(m)[0.5](u, t))),
     main = "Median thickness trajectory")
for (j in seq_along(tsel)) {
  m <- predict_A(qf$fits[[which(taus == .5)]], ug, rep(tsel[j], length(ug)),
                 rep(0, length(ug)), rep(0, length(ug)))
  lines(ug, m - mean(m), col = cols[j], lwd = 2.4)
}
legend("topright", bty = "n", cex = .85, lwd = 2.4, col = cols,
       legend = sprintf("onset t=%.2f", tsel))

## (C) APOE4 effect beta_tau1 across tau, bootstrap band
plot(taus, qb, type = "n", ylim = range(c(blo, bhi, qb, sb$beta[1]), na.rm = TRUE),
     xlab = expression(tau), ylab = expression(beta[list(tau, "APOE4")]),
     main = "APOE4 effect across quantiles")
polygon(c(taus, rev(taus)), c(blo, rev(bhi)), col = adjustcolor("#1f6feb", .16), border = NA)
lines(taus, qb, col = "#1f6feb", lwd = 2.4); points(taus, qb, col = "#1f6feb", pch = 19, cex = .55)
abline(h = sb$beta[1], col = "#d1242f", lwd = 2.2, lty = 4)
abline(h = 0, col = "grey60", lty = 3)
legend("topleft", bty = "n", cex = .82, lwd = c(2.4, 2.2), lty = c(1, 4),
       col = c("#1f6feb", "#d1242f"),
       legend = c(expression(paste("quantile ", beta[tau])),
                  expression(paste("Sun mean ", hat(beta)))))

## (R) induced mean contrast (Q) vs Sun scalar mean coefficient
plot(ug, qD, type = "l", col = "#1f6feb", lwd = 2.4,
     ylim = range(c(qD, sb$beta[1])) + c(-.05, .05),
     xlab = "visit age u (centered)",
     ylab = expression(paste("mean log-contrast  ", Delta(u, t))),
     main = "Induced APOE4 mean contrast")
points(ug, qD, col = "#1f6feb", pch = 19, cex = .4)
abline(h = sb$beta[1], col = "#d1242f", lwd = 2.2, lty = 4)
legend("topright", bty = "n", cex = .82, lwd = c(2.4, 2.2), lty = c(1, 4),
       col = c("#1f6feb", "#d1242f"),
       legend = c("quantile-induced", "Sun scalar mean"))
dev.off()
cat("wrote paper/figs/app_panels.pdf\n")
