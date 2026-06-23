# ============================================================================
# run_one.R — full estimator + inference pipeline on a single dataset.
#
# Estimators (beta_tau on X1,X2):
#   A-CV  : IPCW joint sieve QR, df by subject-level CV
#   A-US  : same, df bumped up (undersmoothed)
#   B-OS  : one-step orthogonalised correction on top of A-CV
#   O3    : projection oracle (true f(0), true g, true G_C)
# Inference: no-correction clustered sandwich (all reps, fast) + optional
#            full-refit / frozen-G_C subject bootstrap (subset of reps).
# ============================================================================
source("R/dgp.R"); source("R/basis.R"); source("R/estimators.R")

prep_mk <- function(d) {
  mk <- d$mk
  mk[mk$dT == 1 & mk$U > 0 & mk$U < mk$T, ]
}

# ---- one-step orthogonal correction ---------------------------------------
one_step <- function(mk, w_row, fitA, tau, a_n,
                     fhat = NULL, ghat = NULL, df_g = 4) {
  df <- fitA$df
  if (is.null(fhat)) {
    fit_lo <- fit_A(mk, w_row, df = df, tau = tau - a_n)
    fit_hi <- fit_A(mk, w_row, df = df, tau = tau + a_n)
    fhat <- sparsity_fhat(fit_lo, fit_hi, mk$U, mk$T, mk$X1, mk$X2, a_n)
  }
  X <- cbind(mk$X1, mk$X2)
  if (is.null(ghat)) {
    proj  <- fit_projection(mk$U, mk$T, X, w_row * fhat, df_g = df_g)
    ghat  <- predict_projection(proj, mk$U, mk$T)
  }
  Xc  <- X - ghat
  qhat <- predict_A(fitA, mk$U, mk$T, mk$X1, mk$X2)
  eps  <- mk$L - qhat
  psi  <- tau - (eps <= 0)
  nsub <- length(unique(mk$id))
  Uvec <- colSums(w_row * Xc * psi) / nsub
  S <- matrix(0, 2, 2)
  for (k in seq_len(nrow(mk))) S <- S + w_row[k] * fhat[k] * tcrossprod(Xc[k, ])
  S <- S / nsub
  beta <- as.numeric(fitA$beta + solve(S, Uvec))
  list(beta = beta, S = S, Xc = Xc, psi = psi, fhat = fhat,
       Scond = kappa(S), w_row = w_row)
}

# clustered sandwich SE for a one-step (no nuisance-estimation correction term)
sandwich_se <- function(mk, os) {
  ids <- mk$id; uid <- unique(ids)
  nsub <- length(uid)
  contrib <- os$w_row * os$Xc * os$psi          # rows x 2
  Z <- rowsum(contrib, group = ids)             # subject-level sums
  V <- crossprod(as.matrix(Z)) / nsub
  Sinv <- solve(os$S)
  Sig <- Sinv %*% V %*% t(Sinv) / nsub
  sqrt(diag(Sig))
}

# ---- feasible estimators on one dataset (point estimates + sandwich SE) ----
estimate_all <- function(d, tau, a_n = 0.05, df_us_bump = 2,
                         df_grid = c(3, 4, 5, 6), oracle = NULL) {
  subj <- d$subj; mk <- prep_mk(d)
  # feasible IPCW
  gc   <- fit_GC(subj)
  GCv  <- GC_at(gc)
  a_s  <- ipcw_weights(subj, GCv)
  w    <- attach_weights(mk, subj, a_s)
  # A-CV
  df_cv  <- select_df_cv(mk, subj, w, tau, df_grid = df_grid)
  fitcv  <- fit_A(mk, w, df = df_cv, tau = tau)
  # A-US
  df_us  <- min(df_cv + df_us_bump, max(df_grid) + df_us_bump)
  fitus  <- fit_A(mk, w, df = df_us, tau = tau)
  # B-OS (feasible f,g; feasible GC weights)
  bos    <- one_step(mk, w, fitcv, tau, a_n)
  bos_se <- sandwich_se(mk, bos)

  out <- list(
    acv = list(beta = as.numeric(fitcv$beta), df = df_cv),
    aus = list(beta = as.numeric(fitus$beta), df = df_us),
    bos = list(beta = bos$beta, se = bos_se, Scond = bos$Scond),
    imse = surface_imse(fitcv, d, tau)
  )

  # O3 projection oracle: true f(0)=phi(z_tau)/sigma(X), true g (frozen pool),
  # true G_C weights. Isolates how much nuisance estimation (f,g,G_C) costs.
  if (!is.null(oracle)) {
    p <- d$params; z <- qnorm(tau)
    a_true_s <- ifelse(subj$dT == 1, 1 / true_GC(subj, p), 0)
    w_o <- a_true_s[match(mk$id, subj$id)]
    f_true <- dnorm(z) / sigma_X(mk$X1, p$sig0, p$sig1)
    ghat_o <- predict_projection(oracle$proj, mk$U, mk$T)
    os_o <- one_step(mk, w_o, fitcv, tau, a_n, fhat = f_true, ghat = ghat_o)
    out$o3 <- list(beta = os_o$beta, se = sandwich_se(mk, os_o),
                   Scond = os_o$Scond)
  }
  out
}

# surface IMSE on the estimable interior region: mean squared error of
# m_hat_tau(u,t) (= surface part, X=0) vs truth mu(u,t)+sig0*z_tau.
surface_imse <- function(fitA, d, tau, delta = 0.3) {
  p <- d$params
  tgrid <- seq(1.0, 7.0, length.out = 16)
  pts <- do.call(rbind, lapply(tgrid, function(tt) {
    ug <- seq(delta, tt - delta, length.out = 12)
    cbind(u = ug, t = tt)
  }))
  mhat  <- predict_A(fitA, pts[, "u"], pts[, "t"], 0, 0)
  mtrue <- true_surface(pts[, "u"], pts[, "t"], tau, p)  # handles surface form + het
  mean((mhat - mtrue)^2)
}

# true censoring survival G_C(W_i | H_i) = exp(-rate_C exp(cC1 X1) * W_i).
true_GC <- function(subj, p) {
  rate <- p$rate_C * exp(p$cC1 * subj$X1)
  exp(-rate * subj$W)
}

# Build projection oracle g_{tau,0}(u,t) from a large pool with TRUE f(0) and
# TRUE G_C weights; df_g matches the feasible projection.
# NOTE: g_{tau,0}(u,t) = E[a (1/sigma(X)) X | u,t] / E[a (1/sigma(X)) | u,t]
# is TAU-FREE: f(0) = phi(z_tau)/sigma(X) and the constant phi(z_tau) cancels in
# the ratio. So build ONE oracle per scenario and reuse across tau.
build_oracle <- function(scenario, n_pool = 8000, df_g = 4,
                         params = default_params(), sig1_val = 0.4) {
  d <- generate_data(n_pool, scenario = scenario, params = params,
                     sig1_val = sig1_val)
  p <- d$params; mk <- prep_mk(d); subj <- d$subj
  a_true_s <- ifelse(subj$dT == 1, 1 / true_GC(subj, p), 0)
  w_o <- a_true_s[match(mk$id, subj$id)]
  w_af <- w_o / sigma_X(mk$X1, p$sig0, p$sig1)      # tau-free density weight
  proj <- fit_projection(mk$U, mk$T, cbind(mk$X1, mk$X2), w_af, df_g = df_g)
  list(proj = proj)
}
