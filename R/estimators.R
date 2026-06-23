# ============================================================================
# estimators.R — censoring model, IPCW, and the A-CV / A-US / B-OS estimators.
#
# Data objects (from dgp.R):
#   subj : one row per retained subject (id, X1, X2, A, T, Z, dT, dD, dC, W, ...)
#   mk   : long marker rows for event subjects' pre-onset visits
#          (id, U, L, X1, X2, T, Z, dT, W, A)
#
# Base case: omega^V = 1 (noninformative visits); IPCW only for loss-to-follow-up.
# Subject weight a_i = dT_i / G_C(W_i | H_i),  H_i = (A_i, X_i),  W_i = T_i - A_i.
# ============================================================================
suppressMessages({library(quantreg); library(survival)})
source("R/basis.R")

# ---- censoring distribution G_C(w | A,X) via Cox on time-since-entry --------
# Censoring "time" is the observed duration to LTFU, right-censored by onset/death:
#   dur = Z - A = min(C, min(T,D)-A) ;  event = dC (LTFU observed).
fit_GC <- function(subj) {
  dur <- pmax(subj$Z - subj$A, 1e-8)
  fit <- coxph(Surv(dur, subj$dC) ~ X1 + X2 + A, data = subj,
               ties = "breslow")
  bh  <- basehaz(fit, centered = FALSE)          # cumulative baseline hazard H0
  H0_fun <- approxfun(bh$time, bh$hazard, method = "constant",
                      yleft = 0, rule = 2)
  lp  <- as.numeric(predict(fit, newdata = subj, type = "lp", reference = "zero"))
  list(fit = fit, H0_fun = H0_fun, subj = subj)
}

# G_C(w | H_i) for each subject i, evaluated at its own w (default w = W_i).
GC_at <- function(gc, w = NULL) {
  sj <- gc$subj
  if (is.null(w)) w <- sj$W
  lp <- as.numeric(predict(gc$fit, newdata = sj, type = "lp", reference = "zero"))
  G  <- exp(-gc$H0_fun(w) * exp(lp))
  pmax(G, 1e-4)                                   # floor to stabilise weights
}

# subject IPCW weights a_i (length = nrow(subj)); 0 for non-event subjects.
ipcw_weights <- function(subj, GC_vec) {
  ifelse(subj$dT == 1, 1 / GC_vec, 0)
}

# attach per-row weight (= subject weight) to marker rows
attach_weights <- function(mk, subj, a_subj) {
  a_subj[match(mk$id, subj$id)]
}

# ---- A-estimator: IPCW-weighted joint sieve QR -----------------------------
# Returns coef beta_tau (X1,X2), the surface spec, and the full coef vector.
fit_A <- function(mk, w_row, df = 5, tau = 0.5) {
  Z <- make_tensor_basis(mk$U, mk$T, df_u = df, df_t = df)
  spec <- attr(Z, "spec")
  bn <- colnames(Z)
  dat <- data.frame(L = mk$L, X1 = mk$X1, X2 = mk$X2, Z)
  # NO separate intercept: tensor B-splines are a partition of unity (columns
  # sum to 1), so they already span the constant. An extra intercept would make
  # the design rank-deficient. The surface m_tau(u,t) carries the level.
  fm <- as.formula(paste("L ~ 0 + X1 + X2 +", paste(bn, collapse = "+")))
  fit <- rq(fm, tau = tau, data = dat, weights = w_row, method = "fn")
  cf <- coef(fit)
  list(beta = cf[c("X1", "X2")], coef = cf, spec = spec, df = df, tau = tau)
}

# predicted quantile q_tau(u,t,x) from a fitted A model
predict_A <- function(fitA, u, t, X1, X2) {
  Z <- eval_tensor_basis(u, t, fitA$spec)
  cf <- fitA$coef
  bcoef <- cf[grep("^B", names(cf))]
  as.numeric(cf["X1"] * X1 + cf["X2"] * X2 + Z %*% bcoef)
}

# ---- subject-level CV to pick df (A-CV) ------------------------------------
# minimise IPCW-weighted check loss on held-out subjects.
rho_tau <- function(r, tau) r * (tau - (r < 0))

select_df_cv <- function(mk, subj, w_row, tau, df_grid = c(3, 4, 5, 6),
                         K = 5, seed = 1) {
  set.seed(seed)
  ids <- unique(mk$id)
  fold <- sample(rep_len(seq_len(K), length(ids)))
  fold_of <- fold[match(mk$id, ids)]
  cvloss <- sapply(df_grid, function(df) {
    tot <- 0
    for (k in seq_len(K)) {
      tr <- mk[fold_of != k, ]; te <- mk[fold_of == k, ]
      wtr <- w_row[fold_of != k]
      ft <- tryCatch(fit_A(tr, wtr, df = df, tau = tau), error = function(e) NULL)
      if (is.null(ft)) { tot <- tot + 1e12; next }
      qhat <- predict_A(ft, te$U, te$T, te$X1, te$X2)
      tot <- tot + sum(w_row[fold_of == k] * rho_tau(te$L - qhat, tau))
    }
    tot
  })
  df_grid[which.min(cvloss)]
}

# ---- B-OS: quantile-spacing sparsity + density-weighted projection + one-step
# f_hat(0 | u,t,x) = 1 / [ (q_{tau+a} - q_{tau-a}) / (2a) ]
sparsity_fhat <- function(fit_lo, fit_hi, u, t, X1, X2, a_n,
                          f_min = 1e-2, f_max = 50) {
  q_hi <- predict_A(fit_hi, u, t, X1, X2)
  q_lo <- predict_A(fit_lo, u, t, X1, X2)
  s <- (q_hi - q_lo) / (2 * a_n)
  s <- pmax(s, 1e-4)                       # rearrange/guard positivity
  fhat <- 1 / s
  pmin(f_max, pmax(f_min, fhat))
}

# density-weighted spline projection g_j(u,t) = E[a f X_j|u,t]/E[a f|u,t],
# fit by weighted least squares of X_j on a (u,t) spline basis with weights a*f.
fit_projection <- function(u, t, X, w_af, df_g = 4) {
  Zg <- make_tensor_basis(u, t, df_u = df_g, df_t = df_g)
  spec <- attr(Zg, "spec")
  G <- cbind(1, Zg)
  WZ <- G * w_af
  XtX <- crossprod(WZ, G)
  diag(XtX) <- diag(XtX) + 1e-8 * mean(diag(XtX))   # tiny ridge for stability
  coefs <- solve(XtX, crossprod(WZ, X))             # ncol(G) x ncol(X)
  list(coefs = coefs, spec = spec)
}
predict_projection <- function(proj, u, t) {
  Zg <- eval_tensor_basis(u, t, proj$spec)
  cbind(1, Zg) %*% proj$coefs                        # n x ncol(X)
}
