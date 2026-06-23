# ============================================================================
# dgp.R  —  MWQTM base-case data-generating process
#
# Target population : individuals free of disease at time-zero (age 65) who
#                     develop disease before death  {D^0 >= T^0}.
# Time scale        : u, t, A, D are ages MINUS 65 (so time-zero = 0).
#                     C is duration since study entry; follow-up ends at A + C.
#
# Marker (log scale), for 0 < u < T:
#   L_i(u) = mu(u, T_i) + b01 X1 + b02 X2 + sigma(X) { sqrt(rho) B_i
#                                                      + sqrt(1-rho) E_il }
#   B_i, E_il ~ iid N(0,1)   =>  bracket ~ N(0,1) marginally.
#   sigma(X)  = sig0 + sig1 * X1     (Scenario 1: sig1 = 0 ; Scenario 2: sig1 > 0)
#
# Analytic truth of the conditional tau-quantile given (T=t, X):
#   Q_tau{L(u) | T=t, X} = m_{tau,0}(u,t) + beta_{tau,1} X1 + beta_{tau,2} X2
#     m_{tau,0}(u,t) = mu(u,t) + sig0 * z_tau
#     beta_{tau,1}   = b01 + sig1 * z_tau          ( <- tau-heterogeneity )
#     beta_{tau,2}   = b02
#   where z_tau = qnorm(tau).
#
# Event time T : Cox-Weibull, lambda_T(t|X) = (kappa/scaleT^kappa) t^{kappa-1}
#                exp(aT1 X1 + aT2 X2).  Larger linear predictor => earlier onset.
# Death  D     : Weibull, independent of marker noise given (T,X) (depends on X).
# Entry  A     : delayed entry; sample retained iff max(A,0) <= min(T,D).
# Censor C     : duration since entry, C _||_ {Y,T,D} | (A,X), Exp(rate_C(X)).
# Visits       : noninformative Poisson-gap process; markers used iff 0<U<T.
# ============================================================================

# ---- baseline mean surface mu(u,t).
#   "smooth"     : 1+0.5 log(1+u)+0.15 t-0.03 u t  (NOT in any finite spline space;
#                  approximation genuinely tested -- default & growing-sieve runs)
#   "bilinear"   : 1+0.4 u+0.15 t-0.03 u t  (a degree-1 polynomial in (u,t), hence
#                  EXACTLY in the tensor cubic B-spline span -> zero approximation
#                  bias -- the exact-sieve check)
mu_surface <- function(u, t, form = "smooth") {
  if (identical(form, "bilinear")) 1.0 + 0.4 * u + 0.15 * t - 0.03 * u * t
  else                             1.0 + 0.5 * log1p(u) + 0.15 * t - 0.03 * u * t
}

# ---- default parameters --------------------------------------------------
default_params <- function() {
  list(
    # marker
    b01 = 0.5, b02 = -0.3,
    sig0 = 0.8,
    rho  = 0.5,          # within-subject correlation of the bracket
    # event time T: onset = "cox" (Cox-Weibull) or "aft" (log-normal AFT, non-PH)
    onset = "cox",
    kappa = 2.0, scaleT = 5.0, aT1 = 0.5, aT2 = 0.3,
    # log-normal AFT onset: log T = muT_aft - (aT1 X1 + aT2 X2) + sigT_aft * N(0,1)
    muT_aft = 1.6, sigT_aft = 0.5,
    # death D (Weibull), tuned for ~10-15% death-before-onset
    kappaD = 1.5, scaleD = 9.0, aD1 = 0.2, aD2 = 0.1,
    # entry A ~ Uniform(A_lo, A_hi); truncation by max(A,0) <= min(T,D)
    A_lo = -2.0, A_hi = 6.0,
    # censoring C ~ Exp(rate_C * exp(cC1 X1)); tuned for ~20% events censored
    rate_C = 0.06, cC1 = 0.10,
    # visit process: gap ~ Exp(mean = visit_gap)
    visit_gap = 0.7,
    # baseline surface form: "smooth" (default) or "bilinear" (exact-sieve)
    surface = "smooth",
    # ---- surface-heterogeneity option (sim-review #13) --------------------
    # when het_surface=TRUE the scale depends on (u,t) as well as X1:
    #   sigma(u,t,X) = sig0 + sig1 X1 + sig2 (t-u) + sig3 (u/t)
    # so the conditional-quantile SURFACE m_{tau,0}(u,t) is NON-parallel in tau:
    #   m_{tau,0}(u,t) = mu(u,t) + [sig0 + sig2 (t-u) + sig3 (u/t)] z_tau
    # (beta_{tau,1} = b01 + sig1 z_tau is UNCHANGED). Tests dispersion that
    # grows as onset nears (u/t -> 1) -- the quantile-trajectory selling point.
    het_surface = FALSE, sig2 = 0, sig3 = 0,
  # ---- error-distribution option (sim-review #14) ----------------------
  # standardized (mean 0, var 1) error for the marker bracket; quantile
  # regression should be agnostic to its shape. Within-subject correlation is
  # preserved via a Gaussian copula (transform the correlated normal bracket).
  #   "normal" (default), "t5" (heavy symmetric tail), "skew" (right-skew lognormal)
  err_dist = "normal"
  )
}

sigma_X <- function(X1, sig0, sig1) sig0 + sig1 * X1

# standardized (mean 0, var 1) error quantile function F_G^{-1}(p)
err_q <- function(p, dist = "normal") {
  if (dist == "normal") return(qnorm(p))
  if (dist == "t5")     return(qt(p, df = 5) / sqrt(5 / 3))          # var t5 = 5/3
  if (dist == "skew") {                                              # std lognormal, s=0.5
    s <- 0.5; m <- exp(s^2 / 2); v <- (exp(s^2) - 1) * exp(s^2)
    return((exp(s * qnorm(p)) - m) / sqrt(v))
  }
  stop("unknown err_dist")
}

# full scale sigma(u,t,X1); reduces to sigma_X when het_surface is off.
sigma_uts <- function(u, t, X1, p) {
  s <- p$sig0 + p$sig1 * X1
  if (isTRUE(p$het_surface)) s <- s + p$sig2 * (t - u) + p$sig3 * (u / t)
  pmax(s, 1e-3)                     # guard positivity
}

# true conditional-quantile surface m_{tau,0}(u,t) at X=0 (for IMSE / trajectory)
true_surface <- function(u, t, tau, p) {
  z <- err_q(tau, if (is.null(p$err_dist)) "normal" else p$err_dist)
  s0 <- p$sig0 + if (isTRUE(p$het_surface)) p$sig2 * (t - u) + p$sig3 * (u / t) else 0
  mu_surface(u, t, p$surface) + s0 * z
}

# analytic truth for a given tau ------------------------------------------
true_beta <- function(tau, p) {
  z <- err_q(tau, p$err_dist %||% "normal")
  c(beta1 = p$b01 + p$sig1 * z, beta2 = p$b02)
}
`%||%` <- function(a, b) if (is.null(a)) b else a
true_intercept_shift <- function(tau, p) p$sig0 * err_q(tau, p$err_dist %||% "normal")

# ---- generate one dataset -------------------------------------------------
# scenario: 1 (homoscedastic, sig1=0) or 2 (heteroscedastic, sig1=sig1_val)
# Returns list(subj = subject-level df, mk = long marker df, params).
generate_data <- function(n, scenario = 1, params = default_params(),
                          sig1_val = 0.4, oversample = 6) {
  p <- params
  p$sig1 <- if (scenario == 1) 0.0 else sig1_val

  # oversample a pool, then keep first n that satisfy left-truncation
  Npool <- ceiling(n * oversample)
  X1 <- rbinom(Npool, 1, 0.5)
  X2 <- rnorm(Npool)

  # event time T: Cox-Weibull (PH) or log-normal AFT (non-PH, for validating
  # that the onset mechanism does not enter beta_tau estimation)
  lpT <- p$aT1 * X1 + p$aT2 * X2
  if (identical(p$onset, "aft")) {
    T0 <- exp(p$muT_aft - lpT + p$sigT_aft * rnorm(Npool))
  } else {
    U  <- runif(Npool)
    T0 <- p$scaleT * (-log(U) * exp(-lpT))^(1 / p$kappa)
  }

  # death D, independent of marker noise given (T,X)
  V  <- runif(Npool)
  lpD <- p$aD1 * X1 + p$aD2 * X2
  D0 <- p$scaleD * (-log(V) * exp(-lpD))^(1 / p$kappaD)

  # entry A (delayed entry)
  A0 <- runif(Npool, p$A_lo, p$A_hi)

  # left truncation: retain iff max(A,0) <= min(T,D)
  keep <- pmax(A0, 0) <= pmin(T0, D0)
  trunc_frac <- 1 - mean(keep)
  idx <- which(keep)
  if (length(idx) < n)
    stop(sprintf("pool too small: kept %d < n=%d (raise oversample)",
                 length(idx), n))
  idx <- idx[seq_len(n)]

  X1 <- X1[idx]; X2 <- X2[idx]; T0 <- T0[idx]; D0 <- D0[idx]; A0 <- A0[idx]

  # censoring C (duration since entry), C _||_ {Y,T,D} | (A,X)
  rateC <- p$rate_C * exp(p$cC1 * X1)
  C0 <- rexp(n, rate = rateC)
  endFU <- A0 + C0                     # age (minus 65) at loss to follow-up

  # observed event quantities
  minTD <- pmin(T0, D0)
  Z  <- pmin(T0, D0, endFU)
  dT <- as.integer(T0 <= D0 & T0 <= endFU)     # Delta^T : onset observed
  dD <- as.integer(D0 <  T0 & D0 <= endFU)     # Delta^D : death before onset
  dC <- as.integer(endFU < minTD)              # Delta^C : loss to follow-up
  W  <- T0 - A0                                # entry -> onset time (for G_C)

  subj <- data.frame(
    id = seq_len(n), X1 = X1, X2 = X2,
    A = A0, T = T0, D = D0, C = C0, endFU = endFU,
    Z = Z, dT = dT, dD = dD, dC = dC, W = W,
    rateC = rateC, sigma = sigma_X(X1, p$sig0, p$sig1)
  )

  # marker: shared random effect B_i; per-visit E_il
  Bi <- rnorm(n)
  # noninformative visit process: Poisson gaps from max(A,0), keep 0<U<min(T,endFU)
  mk_list <- vector("list", n)
  for (i in seq_len(n)) {
    upper <- min(T0[i], endFU[i])      # markers must be pre-onset AND observed
    if (upper <= 0) { mk_list[[i]] <- NULL; next }
    start <- max(A0[i], 0)
    times <- numeric(0); cur <- start
    repeat {
      cur <- cur + rexp(1, rate = 1 / p$visit_gap)
      if (cur >= upper) break
      if (cur > 0) times <- c(times, cur)
    }
    if (!length(times)) { mk_list[[i]] <- NULL; next }
    # per-visit scale: constant sigma(X) unless het_surface adds (u,t) terms
    sg  <- sigma_uts(times, T0[i], X1[i], p)
    brk <- sqrt(p$rho) * Bi[i] + sqrt(1 - p$rho) * rnorm(length(times))   # N(0,1), correlated
    # Gaussian-copula transform to a standardized non-normal marginal (#14):
    # marginal per-visit error is exactly err_q's law, correlation preserved.
    e   <- if ((p$err_dist %||% "normal") == "normal") brk else err_q(pnorm(brk), p$err_dist)
    L   <- mu_surface(times, T0[i], p$surface) + p$b01 * X1[i] + p$b02 * X2[i] + sg * e
    mk_list[[i]] <- data.frame(id = i, U = times, L = L,
                               X1 = X1[i], X2 = X2[i], T = T0[i], Z = Z[i],
                               dT = dT[i], W = W[i], A = A0[i])
  }
  mk <- do.call(rbind, mk_list)
  rownames(mk) <- NULL

  attr(subj, "trunc_frac") <- trunc_frac
  list(subj = subj, mk = mk, params = p,
       summary = c(n = n,
                   trunc_frac = round(trunc_frac, 3),
                   p_onset_obs = round(mean(dT), 3),
                   p_death = round(mean(dD), 3),
                   p_cens = round(mean(dC), 3),
                   mean_visits_event = round(
                     mean(table(factor(mk$id[mk$dT == 1],
                                       levels = which(dT == 1)))), 2)))
}
