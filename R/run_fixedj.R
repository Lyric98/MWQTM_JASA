# ============================================================================
# run_fixedj.R (sim-review round-11 #5) — is the data-adaptive A-CV materially
# different from a FIXED small sieve? CV selects J=3 in 81-86% of replicates
# (Table tab:cvfreq), so we compare, on the SAME replicates (Scenario 2, CC
# estimator, df selected exactly as in the headline = IPCW-weighted subject CV):
#   * A-CV (df = CV-selected)   vs   fixed J=3 (df=3)
#     reporting bias, ESD, bootstrap coverage, and surface IMSE for each, plus the
#     mean |A-CV - J=3| point-estimate difference.
#   * the CV-loss GAPS  CV(4)-CV(3), CV(5)-CV(3): if near zero, the J=3 choice is
#     essentially a tie (adaptivity is noise); if clearly >0, CV genuinely prefers
#     the small sieve. We report the mean gaps and P(gap>0).
# Writes results/fixedj.csv.
# ============================================================================
suppressWarnings(suppressMessages({library(parallel); source("R/run_one.R"); source("R/bootstrap.R")}))

CFG <- list(ns = c(500, 1000, 2000), taus = c(0.25, 0.50, 0.75), scenario = 2,
  R = 400, R_cov = 400, B = 200, df_grid = c(3, 4, 5),
  cores = { sc <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "")); if (is.na(sc) || sc < 1) 4L else sc },
  seed = 20260616)        # SAME base seed family as the headline

cc_w <- function(mk, subj) attach_weights(mk, subj, as.numeric(subj$dT == 1))

# CV losses over the grid (mirrors select_df_cv but returns the whole vector)
cv_losses <- function(mk, subj, w_row, tau, df_grid, K = 5, seed = 1) {
  set.seed(seed); ids <- unique(mk$id)
  fold <- sample(rep_len(seq_len(K), length(ids))); fold_of <- fold[match(mk$id, ids)]
  sapply(df_grid, function(df) {
    tot <- 0
    for (k in seq_len(K)) {
      tr <- mk[fold_of != k, ]; te <- mk[fold_of == k, ]; wtr <- w_row[fold_of != k]
      ft <- tryCatch(fit_A(tr, wtr, df = df, tau = tau), error = function(e) NULL)
      if (is.null(ft)) { tot <- tot + 1e12; next }
      qhat <- predict_A(ft, te$U, te$T, te$X1, te$X2)
      tot <- tot + sum(w_row[fold_of == k] * (te$L - qhat) * (tau - ((te$L - qhat) < 0)))
    }
    tot
  })
}

one_rep <- function(r, n, tau, cfg) {
  set.seed(cfg$seed + r * 7919 + cfg$scenario * 101 + round(tau * 1000))
  p <- default_params(); p$sig1 <- 0.4
  d <- tryCatch(generate_data(n, scenario = cfg$scenario, params = p), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  mk <- prep_mk(d); subj <- d$subj; if (nrow(mk) < 30) return(NULL)
  w <- cc_w(mk, subj)                                   # CC (omega=1): no hat-G_C anywhere
  cvl <- tryCatch(cv_losses(mk, subj, w, tau, cfg$df_grid), error = function(e) rep(NA, 3))
  dfcv <- cfg$df_grid[which.min(cvl)]
  # normalized CV-loss gaps relative to df=3 (relative increase in held-out CC check loss)
  base <- cvl[1]; gap4 <- (cvl[2] - base) / abs(base); gap5 <- (cvl[3] - base) / abs(base)
  fit_acv <- tryCatch(fit_A(mk, w, df = dfcv, tau = tau), error = function(e) NULL)
  fit_j3  <- tryCatch(fit_A(mk, w, df = 3,    tau = tau), error = function(e) NULL)
  if (is.null(fit_acv) || is.null(fit_j3)) return(NULL)
  b_acv <- unname(fit_acv$beta[1]); b_j3 <- unname(fit_j3$beta[1])
  imse_acv <- tryCatch(surface_imse(fit_acv, d, tau), error = function(e) NA)
  imse_j3  <- tryCatch(surface_imse(fit_j3,  d, tau), error = function(e) NA)
  se_acv <- NA; se_j3 <- NA
  if (r <= cfg$R_cov) {
    mk_by <- split(d$mk, d$mk$id); va <- rep(NA, cfg$B); vj <- rep(NA, cfg$B)
    for (b in seq_len(cfg$B)) {
      db <- resample_subjects(d, mk_by); mkb <- prep_mk(db); sjb <- db$subj
      if (is.null(mkb) || nrow(mkb) < 30) next
      wb <- cc_w(mkb, sjb)
      va[b] <- tryCatch(fit_A(mkb, wb, df = dfcv, tau = tau)$beta[1], error = function(e) NA)
      vj[b] <- tryCatch(fit_A(mkb, wb, df = 3,    tau = tau)$beta[1], error = function(e) NA)
    }
    se_acv <- sd(va, na.rm = TRUE); se_j3 <- sd(vj, na.rm = TRUE)
  }
  c(dfcv = dfcv, gap4 = gap4, gap5 = gap5, b_acv = b_acv, b_j3 = b_j3,
    imse_acv = imse_acv, imse_j3 = imse_j3, se_acv = se_acv, se_j3 = se_j3)
}

run_cell <- function(n, tau, cfg) {
  tb <- as.numeric(true_beta(tau, { p <- default_params(); p$sig1 <- 0.4; p })["beta1"])
  M <- do.call(rbind, Filter(Negate(is.null),
        mclapply(seq_len(cfg$R), one_rep, n = n, tau = tau, cfg = cfg,
                 mc.cores = cfg$cores, mc.preschedule = FALSE)))
  ba <- M[, "b_acv"]; bj <- M[, "b_j3"]; sa <- M[, "se_acv"]; sj <- M[, "se_j3"]
  cov <- function(b, s) mean(abs(b - tb) <= 1.96 * s, na.rm = TRUE)
  dimse <- M[, "imse_acv"] - M[, "imse_j3"]; dimse <- dimse[is.finite(dimse)]  # paired A-CV - J=3 IMSE
  data.frame(scenario = cfg$scenario, n = n, tau = tau, R = nrow(M),
    p_df3 = mean(M[, "dfcv"] == 3), gap4_med = median(M[, "gap4"], na.rm = TRUE),
    gap5_med = median(M[, "gap5"], na.rm = TRUE),
    p_gap4_pos = mean(M[, "gap4"] > 0, na.rm = TRUE), p_gap5_pos = mean(M[, "gap5"] > 0, na.rm = TRUE),
    bias_acv = mean(ba, na.rm = TRUE) - tb, bias_j3 = mean(bj, na.rm = TRUE) - tb,
    esd_acv = sd(ba, na.rm = TRUE), esd_j3 = sd(bj, na.rm = TRUE),
    mean_abs_diff = mean(abs(ba - bj), na.rm = TRUE),
    max_abs_diff = max(abs(ba - bj), na.rm = TRUE),
    cov_acv = cov(ba, sa), cov_j3 = cov(bj, sj),
    imse_acv = mean(M[, "imse_acv"], na.rm = TRUE), imse_j3 = mean(M[, "imse_j3"], na.rm = TRUE),
    dimse_mean = mean(dimse), dimse_mcse = sd(dimse) / sqrt(length(dimse)),  # paired DeltaIMSE +/- MCSE
    row.names = NULL)
}

main <- function(cfg = CFG, out = "results/fixedj.csv") {
  dir.create("results", showWarnings = FALSE); rows <- list(); k <- 1
  for (n in cfg$ns) for (tau in cfg$taus) {
    cat(sprintf("[%s] fixedj n=%d tau=%.2f\n", format(Sys.time(), "%H:%M:%S"), n, tau)); flush.console()
    rows[[k]] <- run_cell(n, tau, cfg); k <- k + 1
    write.csv(do.call(rbind, rows), out, row.names = FALSE)
  }
  cat("done ->", out, "\n")
}
if (sys.nframe() == 0) main()
