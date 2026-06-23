# ============================================================================
# bootstrap.R — subject-level bootstrap for A-CV / A-US / B-OS.
#   full-refit : refit G_C, surface, f, g on each resample (correct).
#   frozen-GC  : reuse the ORIGINAL-sample G_C weights (deliberately wrong);
#                isolates the missing IPCW nuisance-estimation variance.
# Returns SEs (bootstrap SD) of the X1 coefficient for each estimator.
# ============================================================================
source("R/run_one.R")

# resample subjects with replacement; carry the original id so frozen weights
# can be matched. New consecutive ids are assigned for clustering.
resample_subjects <- function(d, mk_by = NULL) {
  subj <- d$subj
  bid  <- sample(subj$id, replace = TRUE)
  new_subj <- subj[match(bid, subj$id), ]
  new_subj$orig_id <- bid
  new_subj$id <- seq_len(nrow(new_subj))
  if (is.null(mk_by)) mk_by <- split(d$mk, d$mk$id)
  parts <- vector("list", length(bid))
  for (b in seq_along(bid)) {
    rows <- mk_by[[as.character(bid[b])]]
    if (!is.null(rows) && nrow(rows)) { rows$id <- b; parts[[b]] <- rows }
  }
  new_mk <- do.call(rbind, parts)
  list(subj = new_subj, mk = new_mk, params = d$params)
}

# bootstrap SEs (X1 coef): returns list(full = c(acv,aus,bos), frozen = c(...))
boot_se <- function(d, tau, B = 150, a_n = 0.05, df_fixed = NULL,
                    do_frozen = TRUE) {
  mk_by <- split(d$mk, d$mk$id)
  # original-sample G_C weights, keyed by original subject id (for frozen path)
  gc0 <- fit_GC(d$subj); GC0 <- GC_at(gc0)
  a0  <- ipcw_weights(d$subj, GC0); names(a0) <- d$subj$id

  bmat_full <- matrix(NA, B, 3); bmat_froz <- matrix(NA, B, 3)
  fit_block <- function(mk, subj, a_s) {
    w <- attach_weights(mk, subj, a_s)
    df <- if (is.null(df_fixed)) select_df_cv(mk, subj, w, tau) else df_fixed
    fitcv <- tryCatch(fit_A(mk, w, df = df, tau = tau), error = function(e) NULL)
    if (is.null(fitcv)) return(rep(NA, 3))
    fitus <- tryCatch(fit_A(mk, w, df = df + 2, tau = tau), error = function(e) NULL)
    bos   <- tryCatch(one_step(mk, w, fitcv, tau, a_n), error = function(e) NULL)
    c(fitcv$beta[1],
      if (is.null(fitus)) NA else fitus$beta[1],
      if (is.null(bos))   NA else bos$beta[1])
  }
  for (b in seq_len(B)) {
    db <- resample_subjects(d, mk_by); mk <- prep_mk(db); subj <- db$subj
    if (nrow(mk) < 30) next
    # full-refit: G_C refit on the resampled subjects
    gc <- fit_GC(subj); a_full <- ipcw_weights(subj, GC_at(gc))
    bmat_full[b, ] <- fit_block(mk, subj, a_full)
    # frozen-GC: reuse original weights via carried orig_id
    if (do_frozen) {
      a_froz <- a0[as.character(subj$orig_id)]
      bmat_froz[b, ] <- fit_block(mk, subj, a_froz)
    }
  }
  list(full   = apply(bmat_full, 2, sd, na.rm = TRUE),
       frozen = if (do_frozen) apply(bmat_froz, 2, sd, na.rm = TRUE) else rep(NA, 3))
}
