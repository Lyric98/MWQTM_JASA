# ============================================================================
# compare_sun.R — fair comparison of Sun et al.'s conditional-MEAN estimator vs
# our conditional-QUANTILE estimator, exploiting that the conditional mean can be
# recovered by integrating the fitted conditional quantile process (Q-to-Mean).
# Plan: compare.md. Three Phase-1 DGPs, all on the SAME Monte-Carlo datasets:
#   DGP1  mean-correct, Sun-favorable (proportional mean holds; beta_tau flat)
#   DGP2  mean true + quantile heterogeneity (mean proportional; beta_{tau,1}(tau))
#   DGP3  quantile true, proportional mean FALSE (mean contrast varies in (u,t))
#
# Marker is on the log scale L=log Y; Y=exp(L). Both estimators use the SAME
# complete-event rows (observed-onset, pre-onset visits; omega=1), the fairness
# principle of the plan (same datasets, same rows, common evaluation domain).
#
#   Sun:   E[Y|u,t,x] = mu0(u,t) exp(beta' x)  -- profile-kernel (Nadaraya-Watson
#          smoothing of Y exp(-beta'x); 2-D estimating equation for beta).
#   Ours:  Q_tau{log Y|u,t,x} = m_tau(u,t)+beta_tau'x  (CC sieve QR per tau);
#          M_Q(u,t,x)=\int exp{Q_tau} dtau after pointwise rearrangement.
# ============================================================================
suppressWarnings(suppressMessages({source("R/dgp.R"); source("R/estimators.R")}))

# ----------------------------------------------------------------------------
# 1. analytic truths per DGP (all on the SAME mu/m0 baseline + (b01,b02))
# ----------------------------------------------------------------------------
# DGP knobs carried in p: $dgp in {1,2,3}; reuse default_params() for the rest.
# a(u,t) for DGP3 baseline dispersion (the het-surface form); a0=sig0.
a_disp <- function(u, t, p) p$sig0 + p$a1 * (t - u) + p$a2 * (u / t)

# true conditional quantile of log Y: m_tau(u,t) + beta_tau1 x1 + beta_tau2 x2
true_qsurface <- function(u, t, tau, p) {              # m_tau(u,t) at x=0
  z <- qnorm(tau); mu <- mu_surface(u, t, p$surface)
  if (p$dgp == 1) mu + p$sig0 * z
  else if (p$dgp == 2) mu - 0.5 * p$sig0^2 + p$sig0 * z
  else mu + a_disp(u, t, p) * z                        # DGP3
}
true_beta_tau <- function(tau, p) {                    # (beta_tau1, beta_tau2)
  z <- qnorm(tau)
  if (p$dgp == 1) c(p$b01, p$b02)
  else if (p$dgp == 2) c(p$b01 - 0.5 * ((p$sig0 + p$sig1)^2 - p$sig0^2) + p$sig1 * z, p$b02)
  else c(p$b01 + p$sig1 * z, p$b02)                    # DGP3
}
# true conditional MEAN of Y (the Sun target)
true_mean <- function(u, t, x1, x2, p) {
  mu <- mu_surface(u, t, p$surface)
  if (p$dgp == 1) exp(mu + 0.5 * p$sig0^2 + p$b01 * x1 + p$b02 * x2)
  else if (p$dgp == 2) exp(mu + p$b01 * x1 + p$b02 * x2)         # -1/2 sig^2 cancels
  else {                                                          # DGP3 (non-proportional in X1)
    a <- a_disp(u, t, p); s <- a + p$sig1 * x1
    exp(mu + p$b01 * x1 + p$b02 * x2 + 0.5 * s^2)
  }
}
# true mean log-contrast Delta(u,t) = log E[Y|x1=1] - log E[Y|x1=0]  (x2=0)
true_contrast <- function(u, t, p) {
  if (p$dgp == 3) { a <- a_disp(u, t, p); p$b01 + a * p$sig1 + 0.5 * p$sig1^2 }
  else rep(p$b01, length(u))                                      # constant
}

# ----------------------------------------------------------------------------
# 2. data generation: reuse the survival/truncation/visit design, plug in the
#    DGP-specific log-marker, return Y=exp(L). Same machinery as dgp.R.
# ----------------------------------------------------------------------------
gen_compare <- function(n, dgp, params = NULL, regime = "A", oversample = 6) {
  p <- if (is.null(params)) default_params() else params
  p$dgp <- dgp
  # DGP-specific marker constants
  p$sig1 <- if (dgp == 1) 0 else p$sig1_cmp           # X1 scale-interaction
  if (regime == "B") { p$A_lo <- -1; p$A_hi <- 9; p$rate_C <- 0.10; p$visit_gap <- 1.2 }

  Npool <- ceiling(n * oversample)
  X1 <- rbinom(Npool, 1, 0.5); X2 <- rnorm(Npool)
  lpT <- p$aT1 * X1 + p$aT2 * X2
  T0 <- p$scaleT * (-log(runif(Npool)) * exp(-lpT))^(1 / p$kappa)
  D0 <- p$scaleD * (-log(runif(Npool)) * exp(-(p$aD1 * X1 + p$aD2 * X2)))^(1 / p$kappaD)
  A0 <- runif(Npool, p$A_lo, p$A_hi)
  keep <- pmax(A0, 0) <= pmin(T0, D0); idx <- which(keep)
  if (length(idx) < n) stop("pool too small; raise oversample")
  idx <- idx[seq_len(n)]
  X1 <- X1[idx]; X2 <- X2[idx]; T0 <- T0[idx]; D0 <- D0[idx]; A0 <- A0[idx]
  C0 <- rexp(n, rate = p$rate_C * exp(p$cC1 * X1)); endFU <- A0 + C0
  minTD <- pmin(T0, D0); Z <- pmin(T0, D0, endFU)
  dT <- as.integer(T0 <= D0 & T0 <= endFU); dD <- as.integer(D0 < T0 & D0 <= endFU)
  dC <- as.integer(endFU < minTD); W <- T0 - A0
  subj <- data.frame(id = seq_len(n), X1 = X1, X2 = X2, A = A0, T = T0, D = D0,
                     C = C0, endFU = endFU, Z = Z, dT = dT, dD = dD, dC = dC, W = W)

  Bi <- rnorm(n); mk_list <- vector("list", n)
  for (i in seq_len(n)) {
    upper <- min(T0[i], endFU[i]); if (upper <= 0) next
    start <- max(A0[i], 0); times <- numeric(0); cur <- start
    repeat { cur <- cur + rexp(1, rate = 1 / p$visit_gap); if (cur >= upper) break
             if (cur > 0) times <- c(times, cur) }
    if (!length(times)) next
    m <- length(times)
    eps <- sqrt(p$rho) * Bi[i] + sqrt(1 - p$rho) * rnorm(m)        # N(0,1), correlated
    mu <- mu_surface(times, T0[i], p$surface)
    if (dgp == 1) {
      L <- mu + p$b01 * X1[i] + p$b02 * X2[i] + p$sig0 * eps
    } else if (dgp == 2) {
      sX <- p$sig0 + p$sig1 * X1[i]
      L <- mu + p$b01 * X1[i] + p$b02 * X2[i] - 0.5 * sX^2 + sX * eps
    } else {                                                        # DGP3
      a <- a_disp(times, T0[i], p)
      L <- mu + p$b01 * X1[i] + p$b02 * X2[i] + (a + p$sig1 * X1[i]) * eps
    }
    mk_list[[i]] <- data.frame(id = i, U = times, L = L, Y = exp(L),
                               X1 = X1[i], X2 = X2[i], T = T0[i], dT = dT[i])
  }
  mk <- do.call(rbind, mk_list); rownames(mk) <- NULL
  list(subj = subj, mk = mk, params = p, dgp = dgp)
}

# complete-event rows (observed onset, pre-onset visits) -- the CC sample
cc_rows <- function(d) { mk <- d$mk; mk[mk$dT == 1 & mk$U > 0 & mk$U < mk$T, ] }

# ----------------------------------------------------------------------------
# 3. Sun profile-kernel mean estimator: E[Y|u,t,x]=mu0(u,t) exp(beta'x).
#    Product-Gaussian NW smoothing of Y exp(-beta'x); 2-D EE for beta.
# ----------------------------------------------------------------------------
gkern <- function(d, h) exp(-0.5 * (d / h)^2)

sun_fit <- function(mk, h_mult = 1.0) {
  u <- mk$U; t <- mk$T; X <- cbind(mk$X1, mk$X2); y <- mk$Y; N <- length(y)
  hu <- h_mult * 1.06 * sd(u) * N^(-1/6); ht <- h_mult * 1.06 * sd(t) * N^(-1/6)
  hu <- max(hu, 1e-3); ht <- max(ht, 1e-3)
  # observation-to-observation kernel (leave-one-out for the EE)
  Ku <- outer(u, u, function(a, b) gkern(a - b, hu))
  Kt <- outer(t, t, function(a, b) gkern(a - b, ht))
  W <- Ku * Kt; diag(W) <- 0                     # LOO
  rs <- rowSums(W); rs[rs < 1e-12] <- 1e-12; Wn <- W / rs
  ee <- function(beta) {                         # || sum_i x_i (y_i - muhat_i exp(beta'x_i)) ||^2
    r <- y * exp(-as.numeric(X %*% beta)); muhat <- as.numeric(Wn %*% r)
    resid <- y - muhat * exp(as.numeric(X %*% beta))
    g <- colSums(X * resid); sum(g^2)
  }
  opt <- tryCatch(optim(c(0, 0), ee, method = "BFGS",
                        control = list(reltol = 1e-9, maxit = 200)), error = function(e) NULL)
  beta <- if (is.null(opt)) c(NA, NA) else opt$par
  list(beta = beta, u = u, t = t, X = X, y = y, hu = hu, ht = ht)
}
# predict mu0 and the mean surface at grid points (G x .) for a fitted Sun model
sun_predict <- function(fit, ug, tg, x1, x2) {
  if (any(!is.finite(fit$beta))) return(rep(NA_real_, length(ug)))
  Ku <- outer(ug, fit$u, function(a, b) gkern(a - b, fit$hu))
  Kt <- outer(tg, fit$t, function(a, b) gkern(a - b, fit$ht))
  W <- Ku * Kt; rs <- rowSums(W); rs[rs < 1e-12] <- 1e-12
  r <- fit$y * exp(-as.numeric(fit$X %*% fit$beta))
  mu0 <- as.numeric((W %*% r) / rs)
  mu0 * exp(fit$beta[1] * x1 + fit$beta[2] * x2)
}

# ----------------------------------------------------------------------------
# 4. Q-to-Mean: CC sieve QR per tau, exp, rearrange across tau, trapezoid.
#    Returns beta_tau matrix and a mean-surface predictor over a tau range.
# ----------------------------------------------------------------------------
qmean_fit <- function(mk, tau_grid, df = 4) {
  fits <- lapply(tau_grid, function(tau)
    tryCatch(fit_A(mk, rep(1, nrow(mk)), df = df, tau = tau), error = function(e) NULL))
  ok <- !vapply(fits, is.null, logical(1))
  list(fits = fits, ok = ok, tau = tau_grid)
}
# integrated mean on grid (ug,tg) at (x1,x2); range=c(lo,hi) over tau; trapezoid.
# rearrange = sort exp(Q_tau) across tau at each grid point (monotone).
qmean_predict <- function(qf, ug, tg, x1, x2, range = c(0, 1)) {
  tau <- qf$tau; use <- qf$ok & tau >= range[1] - 1e-9 & tau <= range[2] + 1e-9
  tt <- tau[use]; G <- length(ug)
  Qm <- matrix(NA_real_, G, length(tt)); k <- 1
  for (j in which(use)) {
    Qm[, k] <- exp(predict_A(qf$fits[[j]], ug, tg, rep(x1, G), rep(x2, G))); k <- k + 1
  }
  Qm <- t(apply(Qm, 1, sort))                              # pointwise rearrangement
  # trapezoid over tt, then normalize by (hi-lo) for a central mean; for full
  # mean (range 0,1) the integral over [min tt,max tt] approximates E[Y].
  w <- numeric(length(tt))
  if (length(tt) >= 2) {
    w[1] <- (tt[2] - tt[1]) / 2; w[length(tt)] <- (tt[length(tt)] - tt[length(tt) - 1]) / 2
    if (length(tt) > 2) for (k in 2:(length(tt) - 1)) w[k] <- (tt[k + 1] - tt[k - 1]) / 2
  }
  span <- if (identical(range, c(0, 1))) 1 else (max(tt) - min(tt))
  as.numeric(Qm %*% w) / span                              # central: /span; full: /1
}
# beta_tau point estimates (X1,X2) over the grid of tau
qmean_betatau <- function(qf) {
  t(vapply(seq_along(qf$tau), function(j) {
    if (!qf$ok[j]) return(c(NA, NA)); as.numeric(qf$fits[[j]]$beta) }, numeric(2)))
}

# ----------------------------------------------------------------------------
# 5. evaluation grid (common interior domain) + masked IMSE
# ----------------------------------------------------------------------------
eval_grid <- function(tgrid = seq(2.0, 6.0, length.out = 11), nu = 11, edge = 0.5) {
  pts <- do.call(rbind, lapply(tgrid, function(tt) {
    ug <- seq(edge, tt - edge, length.out = nu); data.frame(u = ug, t = tt) }))
  pts[pts$u > 0, ]
}
# support mask: grid cell estimable iff >= kmin marker rows within +/-hw in (u,t)
supp_mask <- function(pts, mk, hw = 0.6, kmin = 10) {
  vapply(seq_len(nrow(pts)), function(j)
    sum(abs(mk$U - pts$u[j]) <= hw & abs(mk$T - pts$t[j]) <= hw) >= kmin, logical(1))
}
imse <- function(est, tru, mask) { e <- (est - tru)[mask]; mean(e^2, na.rm = TRUE) }
