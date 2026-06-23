# ============================================================================
# sun_port.R — faithful pure-R port of Sun et al.'s estimator (their marker3.cpp
# won't compile here: Rcpp/RcppArmadillo libstdc++ mismatch). The C++ is only for
# speed; we reproduce the EXACT objective. mu_hat/lambda_hat/K1/Kq/f are already
# pure R in their Utilities.R (no C++ dependency) and are sourced verbatim.
#
# Their estimator (supp Thm 3.1): profile partial-likelihood Phi_2 estimates
# theta = alpha + beta; the marker mean coefficient is beta_hat = theta_hat -
# alpha_hat, with alpha_hat from coxph(Surv(A0,z,delta)~X) on the ONSET process.
#   PPL_fast: logPPL = sum_{i in vecI} Y_i [ theta'X_i - log( sum_{k>=i}
#             K_h(A_i-A_k) exp(theta'X_k) ) ]    (data sorted by z then A;
#             k>=i in sorted order encodes the risk set {Z_k >= Z_i}).
#   matrixR is lower-triangular (k>=i); the kernel constant and 1/h cancel in the
#   argmax over theta, so we use the cpp kernel form 1-u^2.
#   vecI = { dlt==1, A>=delta, A<=Z-delta, Z>=2 delta, Z<=tau }.
# ============================================================================
# --- Sun's Utilities.R functions, inlined VERBATIM (their library(rmutil) header
#     won't load here and is unused by these functions). ---
f <- function(u, a) {
  (2*u/3 + u^2/2 - a*u - u^4/12 + u^3*a/3 - u^2*a^2/2 + a^3*u/3 -
     2*u^3/9 - u^4/4 + a*u^3/3 + u^6/18 - 3*u^5*a/15 + 3*u^4*a^2/12 - a^3*u^3/9) * 9/16
}
K1 <- function(u) 0.75 * (1 - u^2) * (abs(u) < 1)
Kq <- function(x, q) {
  sigk1 <- sqrt(0.2)
  2 / (q + 1) * K1(2 / (q + 1) * (x - (q - 1) / 2)) *
    (1 + ((q - 1) / (q + 1) / sigk1)^2 + 2 / sigk1^2 * (1 - q) / (1 + q)^2 * x)
}
lambda_hat <- function(t, b, haz_data) {
  haz_t <- 0
  if (t < b) {
    for (i in 1:nrow(haz_data)) haz_t <- haz_t + Kq((t - haz_data[i,'time'])/b, t/b) * haz_data[i,'h0']/b
  } else {
    for (i in 1:nrow(haz_data)) haz_t <- haz_t + K1((t - haz_data[i,'time'])/b) * haz_data[i,'h0']/b
  }
  haz_t
}
mu_hat <- function(s, t, h1, h2, h3, a, b, X, Y_A, dlt, A, Z, H0) {
  X <- as.matrix(X); xab <- X %*% (a + b); n <- length(Y_A)
  lambda <- lambda_hat(t = t, b = h3, haz_data = H0); numerator <- 0
  for (i in 1:n) {
    if (dlt[i] == 1) {
      Kh_s <- K1((s - A[i])/h1); Kh_t <- K1((t - Z[i])/h1)
      if (t - s < h1 + h1) {
        aa <- (t - s)/h1; intK <- 1 - f(1, aa) + f(aa - 1, aa)
        numerator <- numerator + Kh_s*Kh_t*Y_A[i]/h1/h1/intK
      } else numerator <- numerator + Kh_s*Kh_t*Y_A[i]/h1/h1
    }
  }
  denominator <- 0
  for (k in 1:n) {
    Kh_a <- K1((s - A[k])/h2); q <- (t - s)/h2
    Kh_a <- ifelse(q >= 1, Kh_a, Kh_a/(0.75*(q - q^3/3 + 2/3)))
    denominator <- denominator + Kh_a/h2 * exp(xab[k]) * lambda * (A[k] <= t & t <= Z[k])
  }
  numerator/denominator
}

# valid event-rows (their vectorI), 1-based indices into the sorted long data
sun_vecI <- function(A, dlt, Z, delta, tau)
  which(dlt == 1 & A >= delta & A <= Z - delta & Z >= 2 * delta & Z <= tau)

# negative profile partial log-likelihood in theta=ab (their PPL_fast, ported).
# X: matrix (rows = marker observations, sorted by z then A); Y: marker values.
sun_PPL <- function(ab, X, Y, A, vI, h) {
  n <- length(Y); xab <- as.numeric(X %*% ab); exab <- exp(xab)
  lp <- 0
  for (i in vI) {
    kk <- i:n                                   # k >= i  == risk set {Z_k >= Z_i}
    d <- (A[i] - A[kk]) / h; Kw <- 1 - d * d; Kw[Kw < 0] <- 0   # Epanechnikov (cpp form)
    den <- sum(Kw * exab[kk])
    lp <- lp + Y[i] * (xab[i] - log(den))
  }
  -lp
}

# full Sun beta: theta_hat (PPL) minus alpha_hat (Cox on onset). dl = long marker
# data with columns X1,X2,A,Y,z,delta (sorted by z,A); aw = wide subject data with
# A_0,z,delta,X1,X2 for the Cox fit.
sun_beta <- function(dl, aw, h = NULL, delta = 0.08, tau = 1, init = c(0, 0)) {
  fit <- survival::coxph(survival::Surv(A_0, z, event = delta) ~ X1 + X2, data = aw,
                         control = survival::coxph.control(timefix = FALSE))
  a_hat <- fit$coefficients
  if (is.null(h)) h <- 2.34 * sd(dl$A) * nrow(aw)^(-1/3)
  X <- as.matrix(dl[, c("X1", "X2")])
  vI <- sun_vecI(dl$A, dl$delta, dl$z, delta, tau)
  th <- optim(init, sun_PPL, X = X, Y = dl$Y, A = dl$A, vI = vI, h = h,
              method = "Nelder-Mead", control = list(reltol = 1e-9, maxit = 500))$par
  list(beta = th - a_hat, theta = th, alpha = a_hat, fit = fit, h = h)
}

# Sun baseline-mean surface mu_hat(u,t) at given (u,t) points, using their
# Utilities.R mu_hat (Breslow H0 + kernel). Returns mu0(u,t); mean = mu0 exp(beta'x).
sun_mu_surface <- function(dl, aw, fit, beta, ug, tg, h_mu = NULL) {
  if (is.null(h_mu)) h_mu <- 2.34 * sd(dl$A) * nrow(aw)^(-1/6)
  event <- aw[aw$delta == 1, ]; failure <- sort(unique(event$z))
  H0 <- survival::basehaz(fit, centered = FALSE)
  H0 <- H0[H0$time %in% failure, ]; H0$h0 <- diff(c(0, H0$hazard))
  X <- as.matrix(dl[, c("X1", "X2")])
  vapply(seq_along(ug), function(j)
    tryCatch(mu_hat(s = ug[j], t = tg[j], h1 = h_mu, h2 = h_mu, h3 = h_mu,
                    a = fit$coefficients, b = beta, X = X, Y_A = dl$Y,
                    dlt = dl$delta, A = dl$A, Z = dl$z, H0 = H0),
             error = function(e) NA_real_), numeric(1))
}
