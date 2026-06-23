# ============================================================================
# cf.R — fold-specific CROSS-FITTED one-step (B-OS-CF), the theory-aligned
# estimator of supplement I. None existed before (run_one.R::one_step is NCF).
#
# beta_tilde = sum_k (n_k/n) [ beta_init^(-k) + S_k^{-1} P_{n,k} phi^(-k) ],
# every nuisance (m,f,g,G_C) and the initial beta trained OFF fold k, score and
# bread evaluated ON fold k. CC (omega=1) or IPCW (omega=1/hatG_C).
# ============================================================================
source("R/run_one.R")

# G_C(W_i | H_i) predicted for NEW subjects from a model fit on a training fold
# (GC_at only predicts on the stored training subjects; we need out-of-fold).
GC_predict <- function(gc, newsubj) {
  lp <- as.numeric(predict(gc$fit, newdata = newsubj, type = "lp",
                           reference = "zero"))
  G  <- exp(-gc$H0_fun(newsubj$W) * exp(lp))
  pmax(G, 1e-4)
}

# Cross-fitted one-step. type in {"cc","ipcw"}. Returns beta (X1,X2), the
# one-step correction vs a supplied A-CV beta, and #folds actually used.
one_step_cf <- function(d, tau, a_n = 0.05, K = 5, df_grid = c(3, 4, 5),
                        type = "ipcw", df_g = 4, seed = 1, beta_acv = NULL,
                        df_fixed = NULL, fold_map = NULL, leak = FALSE) {
  subj <- d$subj; mk <- prep_mk(d)
  ev_ids <- unique(mk$id)                 # event subjects contribute to score
  n_ev <- length(ev_ids)
  if (n_ev < 6 * K) return(NULL)
  # reproducible fold assignment WITHOUT clobbering the caller's RNG stream
  # (else repeated calls inside a bootstrap loop freeze resampling)
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", .GlobalEnv)) get(".Random.seed", .GlobalEnv) else NULL
    set.seed(seed)
    on.exit(if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv))
  }
  # Fold by the ORIGINAL (cluster) id so that, in a bootstrap resample, all
  # duplicate copies of one original subject land in the SAME fold -- otherwise
  # one copy can train a nuisance while another is in the held-out score, leaking
  # information across the train/validation split (review #3). Outside bootstrap,
  # orig_id == id and this reduces to plain subject folds.
  # leak=TRUE reproduces the OLD BUG: fold by the resampled id, so duplicate
  # copies of one original subject land in DIFFERENT folds (train/validation
  # leak). leak=FALSE (default) folds by orig_id (the fix).
  cid <- if (!leak && "orig_id" %in% names(subj)) subj$orig_id else subj$id
  ucl <- unique(cid)
  if (!is.null(fold_map)) {
    # FIXED folds: each (original) cluster keeps the fold it was assigned once on
    # the original sample, so bootstrap resamples do NOT redraw the partition
    # (review #7: isolates fold-randomization variance from the leakage fix).
    cl_fold <- fold_map[as.character(ucl)]; names(cl_fold) <- as.character(ucl)
    if (anyNA(cl_fold)) cl_fold[is.na(cl_fold)] <- sample(seq_len(K), sum(is.na(cl_fold)), replace = TRUE)
  } else {
    cl_fold <- sample(rep_len(seq_len(K), length(ucl))); names(cl_fold) <- as.character(ucl)
  }
  subj_fold <- cl_fold[as.character(cid)]; names(subj_fold) <- as.character(subj$id)
  fold_ev  <- subj_fold[as.character(ev_ids)]
  fold_row <- subj_fold[as.character(mk$id)]
  fold_sj  <- subj_fold[as.character(subj$id)]

  acc_beta <- c(0, 0); acc_init <- c(0, 0); acc_os <- c(0, 0); tot_n <- 0L; okK <- 0L
  for (k in seq_len(K)) {
    te_ev <- ev_ids[fold_ev == k]; nk <- length(te_ev)
    if (nk < 2) next
    tr_subj <- subj[fold_sj != k, ]
    mk_tr <- mk[fold_row != k, ]; mk_te <- mk[fold_row == k, ]
    if (nrow(mk_tr) < 30 || nrow(mk_te) < 5) next

    # ---- training-fold subject weights (G_C fit on training subjects) ----
    gc <- NULL
    if (type == "ipcw") {
      gc <- tryCatch(fit_GC(tr_subj), error = function(e) NULL); if (is.null(gc)) next
      a_tr <- ipcw_weights(tr_subj, GC_predict(gc, tr_subj))
    } else a_tr <- as.numeric(tr_subj$dT == 1)
    w_tr <- a_tr[match(mk_tr$id, tr_subj$id)]

    # ---- A-CV (df by CV), sparsity f-hat, projection g-hat: all on training --
    df   <- if (!is.null(df_fixed)) df_fixed else
            tryCatch(select_df_cv(mk_tr, tr_subj, w_tr, tau, df_grid = df_grid),
                     error = function(e) 4)
    fitA <- tryCatch(fit_A(mk_tr, w_tr, df = df, tau = tau), error = function(e) NULL)
    if (is.null(fitA)) next
    binit  <- as.numeric(fitA$beta)
    fit_lo <- tryCatch(fit_A(mk_tr, w_tr, df = df, tau = tau - a_n), error = function(e) NULL)
    fit_hi <- tryCatch(fit_A(mk_tr, w_tr, df = df, tau = tau + a_n), error = function(e) NULL)
    if (is.null(fit_lo) || is.null(fit_hi)) next
    f_tr <- sparsity_fhat(fit_lo, fit_hi, mk_tr$U, mk_tr$T, mk_tr$X1, mk_tr$X2, a_n)
    proj <- tryCatch(fit_projection(mk_tr$U, mk_tr$T, cbind(mk_tr$X1, mk_tr$X2),
                                    w_tr * f_tr, df_g = df_g), error = function(e) NULL)
    if (is.null(proj)) next

    # ---- test-fold weights (out-of-fold G_C), score and bread ----
    subj_te <- subj[subj$id %in% te_ev, ]
    if (type == "ipcw") a_te <- ipcw_weights(subj_te, GC_predict(gc, subj_te))
    else                a_te <- as.numeric(subj_te$dT == 1)
    w_te <- a_te[match(mk_te$id, subj_te$id)]
    Xte  <- cbind(mk_te$X1, mk_te$X2)
    gte  <- predict_projection(proj, mk_te$U, mk_te$T)
    fte  <- sparsity_fhat(fit_lo, fit_hi, mk_te$U, mk_te$T, mk_te$X1, mk_te$X2, a_n)
    qte  <- predict_A(fitA, mk_te$U, mk_te$T, mk_te$X1, mk_te$X2)
    psi  <- tau - (mk_te$L - qte <= 0)
    Xc   <- Xte - gte
    Uk   <- colSums(w_te * Xc * psi) / nk
    Sk   <- matrix(0, 2, 2)
    for (r in seq_len(nrow(mk_te))) Sk <- Sk + w_te[r] * fte[r] * tcrossprod(Xc[r, ])
    Sk <- Sk / nk
    os_k <- tryCatch(as.numeric(solve(Sk, Uk)), error = function(e) NULL)
    if (is.null(os_k) || any(!is.finite(os_k))) next
    bk <- binit + os_k
    acc_beta <- acc_beta + nk * bk; acc_init <- acc_init + nk * binit
    acc_os <- acc_os + nk * os_k; tot_n <- tot_n + nk; okK <- okK + 1L
  }
  if (tot_n == 0L || okK < K - 1L) return(NULL)
  beta <- acc_beta / tot_n
  # decompose the gap to A-CV into sample-SPLIT and pure ONE-STEP parts (#2):
  # beta_tilde - beta_ACV = (sum n_k/n beta_init^{(-k)} - beta_ACV) + sum n_k/n S_k^{-1} U_k
  c_split <- if (!is.null(beta_acv)) acc_init / tot_n - beta_acv else c(NA, NA)
  c_os    <- acc_os / tot_n
  list(beta = beta, K_used = okK,
       correction = if (!is.null(beta_acv)) beta - beta_acv else c(NA, NA),
       c_split = c_split, c_os = c_os)
}
