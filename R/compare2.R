# ============================================================================
# compare2.R — Sun (conditional-MEAN) vs ours (conditional-QUANTILE + Q-to-Mean),
# built on Sun et al.'s REAL code and REAL data design.
#
#   * Data design is Sun's Simulation_scen123.R, verbatim in structure:
#       X2 ~ Unif(0,1); onset T0 ~ Weibull(shape=2, scale=exp(-(X1+X2)/2))  [a0=(1,1)]
#       entry A0 ~ Exp(2.52); censoring C0 ~ Exp(0.57); left-truncate A0<=T0;
#       z=min(A0+C0,T0), delta=1{event}; FIXED visits A0+0.1k, k=0..6, kept if <=z;
#       baseline mean surface mu2(s,t)=pweibull(s,3,0.7t)+0.2.
#     One deliberate change: X1 is a BINARY group indicator (Bernoulli .5), not
#     Unif(-.5,.5). Reason: the estimand that separates the methods is a GROUP
#     contrast log E[Y|X1=1]-log E[Y|X1=0]; a binary X1 makes that contrast a
#     clean (u,t)-surface AND keeps our linear-in-x quantile model correctly
#     specified, so the comparison isolates mean-vs-quantile (not spline error).
#     Sun's estimator is agnostic to the covariate law; everything else is theirs.
#
#   * Sun arm is their GENUINE estimator: beta_hat = optim(PPL_fast) - coxph(onset)
#     using their compiled marker3.cpp (PPL_fast/vectorI/matrixR), and mu_hat()
#     from their Utilities.R for the mean surface. Pure-R port (sun_port.R) is the
#     validated fallback (matches the C++ to 9e-4).
#
#   * DGP1 mean-correct, Sun-favorable (Sun's exact Gamma-multiplicative marker;
#         proportional mean holds, beta_tau flat). HONEST baseline: Sun efficient.
#     DGP2 mean still proportional, but covariate drives DISPERSION:
#         log Y = log mu2 + b1 X1 + b2 X2 - 1/2 sx^2 + sx*eps, sx=sig0+sig1 X1.
#         E[Y] proportional (Sun correct for the mean) but beta_{tau,1}=b1
#         -1/2[(sig0+sig1)^2-sig0^2]+sig1 z_tau varies with tau -> only Q sees it.
#     DGP3 proportional mean is FALSE (the headline):
#         log Y = log mu2 + b1 X1 + b2 X2 + (a(u,t)+sig1 X1) eps, a(u,t)>0 varying.
#         beta_{tau,1}=b1+sig1 z_tau (our model correct); induced mean log-contrast
#         Delta(u,t)=b1+sig1 a(u,t)+1/2 sig1^2 VARIES in (u,t), so no proportional
#         slope exists. Sun's scalar beta1 is a pseudo-average; Q tracks the surface.
# ============================================================================
suppressWarnings(suppressMessages({
  library(survival)
  source("R/estimators.R")     # fit_A, predict_A (our sieve QR)
  source("R/sun_port.R")       # mu_hat/lambda_hat/K1/Kq/f + sun_PPL/sun_vecI (fallback)
}))

# ---- Sun arm: compile their REAL C++ once (forks inherit it) ---------------
.SUN_CPP <- FALSE
try_compile_sun_cpp <- function(
    path = "Sun_code_supp/code/marker3.cpp") {
  if (.SUN_CPP) return(TRUE)
  ok <- tryCatch({ Rcpp::sourceCpp(path); TRUE }, error = function(e) {
    message("marker3.cpp compile failed (", conditionMessage(e),
            "); using validated R port."); FALSE })
  assign(".SUN_CPP", ok, envir = .GlobalEnv); ok
}

# ----------------------------------------------------------------------------
# baseline mean surface + DGP knobs (Sun scenario 2)
# ----------------------------------------------------------------------------
mu2 <- function(s, t) pweibull(s, shape = 3, scale = 0.7 * t) + 0.2

cmp2_params <- function() list(
  b01 = 0.5, b02 = 1.0,            # Sun's b0 = c(0.5, 1)
  sig0 = 0.55, sig1 = 0.45,        # log-normal dispersion; X1 scale-interaction (DGP2/3)
  a0 = 0.50, a1 = 0.65, a2 = -0.45,# a(u,t)=a0+a1(t-u)+a2(u/t)>0; makes the DGP3 mean
                                   # log-contrast b1+sig1*a(u,t)+.5 sig1^2 vary in (u,t)
  rho = 0.5,                       # within-subject correlation of the bracket
  # ---- DGP4: TAIL-SPECIFIC effect (X1=1 has a p4 fraction of high-responders) ----
  s4 = 0.5, p4 = 0.15, m4 = 2.0,   # control N(0,s4^2); treated mixes in N(m4,s4^2) w.p. p4
  # ---- DGP5: HEAVY-TAILED error (Student-t on the log scale; pure location shift) ----
  nu5 = 5, s5 = 0.5,               # standardized t_nu obs noise; E[Y] is undefined (t tail)
  delta = 0.08, tau_sun = 1.0)     # Sun's vecI interior trimming

# baseline dispersion a(u,t) for DGP3 (varies in (u,t))
a_disp2 <- function(u, t, p) p$a0 + p$a1 * (t - u) + p$a2 * (u / t)

# DGP4 mixture quantile of W|X1=1: (1-p4)N(0,s4^2)+p4 N(m4,s4^2)
qmix4 <- function(tau, p) {
  Fm <- function(w) (1 - p$p4) * pnorm(w / p$s4) + p$p4 * pnorm((w - p$m4) / p$s4)
  vapply(tau, function(tt) uniroot(function(w) Fm(w) - tt,
    lower = -6 * p$s4, upper = p$m4 + 6 * p$s4, tol = 1e-9)$root, numeric(1))
}

# ----------------------------------------------------------------------------
# analytic truths (binary X1; contrasts are X1=1 vs X1=0 at fixed x2)
# ----------------------------------------------------------------------------
true_betatau2 <- function(tau, p, dgp) {
  z <- qnorm(tau)
  if (dgp == 1) c(p$b01, p$b02)                                   # log-Gamma: flat slope
  else if (dgp == 2) c(p$b01 - 0.5 * ((p$sig0 + p$sig1)^2 - p$sig0^2) + p$sig1 * z, p$b02)
  else if (dgp == 3) c(p$b01 + p$sig1 * z, p$b02)                 # DGP3
  else if (dgp == 4) c(qmix4(tau, p) - p$s4 * z, p$b02)           # DGP4 tail kink (mixture Q - control Q)
  else c(p$b01, p$b02)                                            # DGP5 pure location shift: flat
}
# DGP4 constant mean log-contrast log{(1-p4)+p4 e^{m4}}
ctr4 <- function(p) log((1 - p$p4) + p$p4 * exp(p$m4))
true_mean2 <- function(u, t, x1, x2, p, dgp) {
  lm <- log(mu2(u, t))
  if (dgp == 1) mu2(u, t) * exp(p$b01 * x1 + p$b02 * x2)          # E[Gamma]=1
  else if (dgp == 2) mu2(u, t) * exp(p$b01 * x1 + p$b02 * x2)     # -1/2 sx^2 cancels
  else if (dgp == 3) { a <- a_disp2(u, t, p); s <- a + p$sig1 * x1
         exp(lm + p$b01 * x1 + p$b02 * x2 + 0.5 * s^2) }          # DGP3 non-proportional
  else if (dgp == 4) exp(lm + p$b02 * x2 + 0.5 * p$s4^2 + x1 * ctr4(p))  # proportional (tail-driven avg)
  else rep(NA_real_, length(u))                                   # DGP5: E[Y] undefined (t tail)
}
true_contrast2 <- function(u, t, p, dgp) {                        # logE[Y|1]-logE[Y|0]
  if (dgp == 3) { a <- a_disp2(u, t, p); p$b01 + a * p$sig1 + 0.5 * p$sig1^2 }
  else if (dgp == 4) rep(ctr4(p), length(u))                      # constant
  else if (dgp == 5) rep(NA_real_, length(u))                     # undefined mean
  else rep(p$b01, length(u))                                      # constant = b1
}

# ----------------------------------------------------------------------------
# data generator: Sun's design exactly, binary X1, DGP-specific marker.
# Returns dl (long, sorted by z then A; cols X1,X2,delta,z,A,Y,obs) + aw (wide:
# id,X1,X2,A_0,z,delta) for the Cox onset fit -- the two objects Sun's code uses.
# ----------------------------------------------------------------------------
gen_sun <- function(n, dgp, p, oversample = 4) {
  pool <- ceiling(n * oversample)
  repeat {
    X1 <- rbinom(pool, 1, 0.5); X2 <- runif(pool)
    T0 <- rweibull(pool, shape = 2, scale = exp(-(X1 + X2) / 2))   # Sun a0=c(1,1)
    A0 <- rexp(pool, rate = 2.52); C0 <- rexp(pool, rate = 0.57)
    keep <- which(A0 <= T0)                                        # left-truncation
    if (length(keep) >= n) break
    pool <- pool * 2
  }
  idx <- keep[seq_len(n)]
  X1 <- X1[idx]; X2 <- X2[idx]; T0 <- T0[idx]; A0 <- A0[idx]; C0 <- C0[idx]
  censor <- T0 > A0 + C0
  z <- ifelse(censor, A0 + C0, T0); delta <- as.numeric(!censor)

  Bi <- rnorm(n); rows <- vector("list", n)
  for (i in seq_len(n)) {
    aks <- A0[i] + 0.1 * (0:6); aks <- aks[aks <= z[i]]            # fixed visits
    m <- length(aks); if (!m) next
    s <- aks; t <- T0[i]; lm <- log(mu2(s, t))
    if (dgp == 1) {
      Y <- mu2(s, t) * exp(p$b01 * X1[i] + p$b02 * X2[i]) * rgamma(m, shape = 10, scale = 0.1)
    } else if (dgp == 4) {                                        # tail-specific: X1=1 mixes a shifted upper tail
      eps <- sqrt(p$rho) * Bi[i] + sqrt(1 - p$rho) * rnorm(m)
      shift <- if (X1[i] == 1 && rbinom(1, 1, p$p4) == 1) p$m4 else 0
      Y <- exp(lm + p$b02 * X2[i] + p$s4 * eps + shift)           # NO base X1 location effect
    } else if (dgp == 5) {                                        # heavy-tailed: standardized t obs noise
      et <- rt(m, df = p$nu5) / sqrt(p$nu5 / (p$nu5 - 2))
      eps <- sqrt(p$rho) * Bi[i] + sqrt(1 - p$rho) * et            # normal frailty + heavy obs noise
      Y <- exp(lm + p$b01 * X1[i] + p$b02 * X2[i] + p$s5 * eps)
    } else {
      eps <- sqrt(p$rho) * Bi[i] + sqrt(1 - p$rho) * rnorm(m)      # N(0,1), correlated
      if (dgp == 2) { sx <- p$sig0 + p$sig1 * X1[i]
        L <- lm + p$b01 * X1[i] + p$b02 * X2[i] - 0.5 * sx^2 + sx * eps
      } else { a <- a_disp2(s, t, p)
        L <- lm + p$b01 * X1[i] + p$b02 * X2[i] + (a + p$sig1 * X1[i]) * eps }
      Y <- exp(L)
    }
    rows[[i]] <- data.frame(id = i, X1 = X1[i], X2 = X2[i], delta = delta[i],
                            z = z[i], A = s, Y = Y, obs = seq_len(m) - 1)
  }
  dl <- do.call(rbind, rows); dl <- dl[order(dl$z, dl$A), ]; rownames(dl) <- NULL
  aw <- data.frame(id = seq_len(n), X1 = X1, X2 = X2, A_0 = A0, z = z, delta = delta)
  list(dl = dl, aw = aw, p = p, dgp = dgp,
       frac = c(event = mean(delta), trunc_keep = length(keep) / pool))
}

# our CC sample (complete-event, observed-onset rows): mk with U,L,X1,X2,T
cc_rows2 <- function(d) {
  e <- d$dl[d$dl$delta == 1, ]
  data.frame(id = e$id, U = e$A, L = log(e$Y), X1 = e$X1, X2 = e$X2, T = e$z)
}

# ----------------------------------------------------------------------------
# Sun beta: their genuine PPL (C++ if compiled, else validated R port) minus the
# Cox onset alpha. Mirrors Simulation_scen123.R lines 102-114 exactly.
# ----------------------------------------------------------------------------
sun_beta2 <- function(d, h = NULL, init = c(0, 0)) {
  dl <- d$dl; aw <- d$aw
  fit <- coxph(Surv(A_0, z, event = delta) ~ X1 + X2, data = aw,
               control = coxph.control(timefix = FALSE))
  a_hat <- fit$coefficients
  if (is.null(h)) h <- 2.34 * sd(dl$A) * nrow(aw)^(-1/3)
  X <- as.matrix(dl[, c("X1", "X2")])
  if (.SUN_CPP) {
    vecI <- vectorI(tau = 1, delta = 0.08, Z = dl$z, A = dl$A, dlt = dl$delta)
    matR <- matrixR(h = h, Z = dl$z, A = dl$A)
    th <- optim(init, PPL_fast, X = X, Y_A = dl$Y, vecI = vecI, matR = matR)$par
  } else {
    vI <- sun_vecI(dl$A, dl$delta, dl$z, 0.08, 1)
    th <- optim(init, sun_PPL, X = X, Y = dl$Y, A = dl$A, vI = vI, h = h,
                method = "Nelder-Mead", control = list(reltol = 1e-9, maxit = 500))$par
  }
  list(beta = as.numeric(th - a_hat), theta = th, alpha = a_hat, fit = fit, h = h)
}

# ---- Sun's Sec. 3.5 STRATIFIED-BASELINE extension ---------------------------
# Sun et al.\ allow a baseline stratified by a covariate. We fit their GENUINE
# estimator (Cox onset + PPL beta + kernel mu_hat) SEPARATELY within each X1
# group, so the baseline mu_0^{(g)}(u,t) differs freely by group. The induced mean
# log-contrast log mu_1(u,t) - log mu_0(u,t) is then nonparametric in (u,t) -- the
# most flexible a mean method can be about the covariate effect. This is the
# strongest, fully-faithful Sun comparator (no proportional-X1 restriction at all).
# Fits Sun on one group (covariate = X2 only; X1 is constant within group).
sun_mean_group <- function(dl, aw, ug, tg, x2, h_mu_mult = 0.5) {
  fit <- coxph(Surv(A_0, z, event = delta) ~ X2, data = aw,
               control = coxph.control(timefix = FALSE))
  a2 <- unname(fit$coefficients)
  hb <- 2.34 * sd(dl$A) * nrow(aw)^(-1/3)
  X <- as.matrix(dl[, "X2", drop = FALSE])
  if (.SUN_CPP) {
    vecI <- vectorI(tau = 1, delta = 0.08, Z = dl$z, A = dl$A, dlt = dl$delta)
    matR <- matrixR(h = hb, Z = dl$z, A = dl$A)
    th <- optim(0, PPL_fast, X = X, Y_A = dl$Y, vecI = vecI, matR = matR,
                method = "Brent", lower = -6, upper = 6)$par
  } else {
    vI <- sun_vecI(dl$A, dl$delta, dl$z, 0.08, 1)
    th <- optim(0, sun_PPL, X = X, Y = dl$Y, A = dl$A, vI = vI, h = hb,
                method = "Brent", lower = -6, upper = 6)$par
  }
  b2 <- th - a2
  h_mu <- h_mu_mult * 2.34 * sd(dl$A) * nrow(aw)^(-1/6)
  failure <- sort(unique(aw$z[aw$delta == 1]))
  H0 <- survival::basehaz(fit, centered = FALSE)
  H0 <- H0[H0$time %in% failure, ]; H0$h0 <- diff(c(0, H0$hazard))
  mu0 <- vapply(seq_along(ug), function(j)
    tryCatch(mu_hat(s = ug[j], t = tg[j], h1 = h_mu, h2 = h_mu, h3 = h_mu,
                    a = a2, b = b2, X = X, Y_A = dl$Y, dlt = dl$delta,
                    A = dl$A, Z = dl$z, H0 = H0), error = function(e) NA_real_), numeric(1))
  mu0[!is.finite(mu0) | mu0 <= 0] <- NA                       # sparse per-group cells -> drop
  mu0 * exp(b2 * x2)
}
# Sun-stratified mean log-contrast (X1=1 vs X1=0) at grid points; non-finite -> NA.
sun_strat_contrast <- function(d, ug, tg, x2 = 0.5, h_mu_mult = 0.5) {
  d1 <- list(dl = d$dl[d$dl$X1 == 1, ], aw = d$aw[d$aw$X1 == 1, ])
  d0 <- list(dl = d$dl[d$dl$X1 == 0, ], aw = d$aw[d$aw$X1 == 0, ])
  m1 <- sun_mean_group(d1$dl, d1$aw, ug, tg, x2, h_mu_mult)
  m0 <- sun_mean_group(d0$dl, d0$aw, ug, tg, x2, h_mu_mult)
  ct <- log(m1) - log(m0); ct[!is.finite(ct)] <- NA; ct
}

# Sun mean surface mu0(u,t) exp(beta'x) at grid points, using their mu_hat.
sun_mean2 <- function(d, sb, ug, tg, x1, x2, h_mu = NULL) {
  dl <- d$dl
  if (is.null(h_mu)) h_mu <- 2.34 * sd(dl$A) * nrow(d$aw)^(-1/6)
  failure <- sort(unique(d$aw$z[d$aw$delta == 1]))
  H0 <- survival::basehaz(sb$fit, centered = FALSE)
  H0 <- H0[H0$time %in% failure, ]; H0$h0 <- diff(c(0, H0$hazard))
  X <- as.matrix(dl[, c("X1", "X2")])
  mu0 <- vapply(seq_along(ug), function(j)
    tryCatch(mu_hat(s = ug[j], t = tg[j], h1 = h_mu, h2 = h_mu, h3 = h_mu,
                    a = sb$alpha, b = sb$beta, X = X, Y_A = dl$Y, dlt = dl$delta,
                    A = dl$A, Z = dl$z, H0 = H0),
             error = function(e) NA_real_), numeric(1))
  mu0 * exp(sb$beta[1] * x1 + sb$beta[2] * x2)
}

# ----------------------------------------------------------------------------
# our quantile arm: CC sieve QR per tau, integrate exp(Q_tau) -> mean (reused).
# ----------------------------------------------------------------------------
qfit2 <- function(mk, tau_grid, df = 3) {
  fits <- lapply(tau_grid, function(tau)
    tryCatch(fit_A(mk, rep(1, nrow(mk)), df = df, tau = tau), error = function(e) NULL))
  list(fits = fits, ok = !vapply(fits, is.null, logical(1)), tau = tau_grid)
}
qbeta2 <- function(qf) t(vapply(seq_along(qf$tau), function(j)
  if (!qf$ok[j]) c(NA, NA) else as.numeric(qf$fits[[j]]$beta), numeric(2)))

# --- conditional mean from the fitted quantile process ----------------------
# Both estimators live in the multiplicative log-link regime (Sun: E[Y]=mu0 e^{b'x};
# ours: log Y has conditional quantiles Q_tau). There the conditional mean is
#   log E[Y|u,t,x] = Q_{.5}(log Y|u,t,x) + 1/2 * sigma(u,t,x)^2,
# i.e. the median LEVEL plus a Duan-style retransformation correction. We read
# BOTH the level and the local scale sigma directly off the fitted quantile
# process: regress Q_tau on z_tau=Phi^{-1}(tau) across tau -> (intercept at z=0 =
# median level, slope = local SD). This recovers E[Y] WITHOUT exponentiating the
# extreme tail quantiles (the tail-unstable step), so it is robust under the
# heavy-dispersion cells that arise when the covariate inflates the variance.
# Returns log E[Y] on the (ug,tg) grid (exp() it for the mean).
qmean_ln <- function(qf, ug, tg, x1, x2) {
  use <- which(qf$ok); tt <- qf$tau[use]; z <- qnorm(tt); G <- length(ug)
  Qm <- vapply(use, function(j)
    predict_A(qf$fits[[j]], ug, tg, rep(x1, G), rep(x2, G)), numeric(G))  # G x ntau
  if (is.null(dim(Qm))) Qm <- matrix(Qm, nrow = G)
  zc <- z - mean(z); Szz <- sum(zc^2)
  slope  <- as.numeric(Qm %*% zc) / Szz                     # local SD sigma(u,t,x)
  inter0 <- rowMeans(Qm) - slope * mean(z)                  # median level at z=0
  inter0 + 0.5 * slope^2                                     # log E[Y]
}
# Q-INT: the genuine quantile-to-mean integral  E[Y]=\int_0^1 exp(Q_tau) dtau,
# evaluated by trapezoid over the (rearranged) tau-grid PLUS a constant-tail rule
# (rectangles [0,tau_min] and [tau_max,1] at the boundary quantiles) so it targets
# the full [0,1] integral. This is the model-free reconstruction; it is tail-
# sensitive by construction (the [tau_max,1] rectangle cannot capture a heavy
# upper tail), which is exactly the instability we document for DGP4/DGP5.
# Returns log E[Y] on the (ug,tg) grid.
qmean_int <- function(qf, ug, tg, x1, x2) {
  use <- which(qf$ok); tt <- qf$tau[use]; n <- length(tt); G <- length(ug)
  Qm <- matrix(NA_real_, G, n)
  for (k in seq_len(n))
    Qm[, k] <- exp(predict_A(qf$fits[[use[k]]], ug, tg, rep(x1, G), rep(x2, G)))
  Qm <- t(apply(Qm, 1, sort))                                     # pointwise rearrangement
  w <- numeric(n)
  w[1] <- (tt[2] - tt[1]) / 2 + tt[1]                             # + rectangle [0, tau_min]
  w[n] <- (tt[n] - tt[n - 1]) / 2 + (1 - tt[n])                   # + rectangle [tau_max, 1]
  if (n > 2) for (k in 2:(n - 1)) w[k] <- (tt[k + 1] - tt[k - 1]) / 2
  log(as.numeric(Qm %*% w))                                       # log E[Y]
}

# ----------------------------------------------------------------------------
# evaluation grid on Sun's interior domain + support mask + IMSE
# ----------------------------------------------------------------------------
# common INTERIOR evaluation domain. Both estimators are nonparametric in (u,t);
# evaluating in the interior (away from the steep u->0 boundary where the marker
# surface mu2 rises sharply and ANY local kernel is boundary-biased) is the plan's
# fairness principle -- it does not flatter or penalise either method.
eval_grid2 <- function(tg = seq(0.50, 0.95, length.out = 6), nu = 6, lo = 0.35, edge = 0.10) {
  pts <- do.call(rbind, lapply(tg, function(tt) {
    ug <- seq(lo, tt - edge, length.out = nu); data.frame(u = ug, t = tt) }))
  pts[pts$u > 0 & pts$u < pts$t, ]
}
supp_mask2 <- function(pts, mk, hw = 0.12, kmin = 8)
  vapply(seq_len(nrow(pts)), function(j)
    sum(abs(mk$U - pts$u[j]) <= hw & abs(mk$T - pts$t[j]) <= hw) >= kmin, logical(1))
imse2 <- function(est, tru, mask) { e <- (est - tru)[mask]; mean(e^2, na.rm = TRUE) }
