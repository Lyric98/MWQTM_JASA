# ============================================================================
# run_crossing.R (sim-review #18) — DENSE quantile-crossing diagnostics.
# Fit A-CV at tau in {0.1,0.25,0.5,0.75,0.9}, predict on a dense interior grid
# (30x30) x covariate profiles, and report: raw crossing rate split INTERIOR vs
# BOUNDARY, mean/max magnitude OVER INCIDENTS, position of incidents (mean u/t,
# t-u -> are they near the low-information boundary?), and the effect of monotone
# rearrangement on surface IMSE and quantile CALIBRATION (empirical tau-coverage).
# Writes results/crossing.csv.
# ============================================================================
suppressWarnings(suppressMessages({library(parallel); source("R/run_one.R")}))

CFG <- list(ns = c(500, 1000, 2000), scenarios = c(1, 2),
  taus = c(0.10, 0.25, 0.50, 0.75, 0.90), R = 300, df_grid = c(3, 4, 5),
  cores = { sc <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK","")); if (is.na(sc)||sc<1) 4L else sc },
  seed = 20260620, edge = 0.5)

eval_grid <- function(cfg) {
  tg <- seq(1.5, 6.0, length.out = 30)
  pts <- do.call(rbind, lapply(tg, function(tt) {
    ug <- seq(0.3, tt - 0.3, length.out = 30); cbind(u = ug, t = tt) }))
  pts <- as.data.frame(pts)
  e <- cfg$edge
  pts$interior <- with(pts, u > e & u < t - e & t > min(tg) + e & t < max(tg) - e)
  # profiles: X1 in {0,1} x X2 at 5 standard-normal quantiles
  Xs <- expand.grid(X1 = c(0, 1), X2 = qnorm(c(.1, .3, .5, .7, .9)))
  list(pts = pts, Xs = Xs)
}

one_rep <- function(r, n, scenario, cfg, grd) {
  set.seed(cfg$seed + r * 7919 + scenario * 101)
  p <- default_params(); p$sig1 <- if (scenario == 1) 0 else 0.4
  d <- tryCatch(generate_data(n, scenario = scenario, params = p), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  mk <- prep_mk(d); subj <- d$subj; if (nrow(mk) < 30) return(NULL)
  w <- attach_weights(mk, subj, ipcw_weights(subj, GC_at(fit_GC(subj))))
  df <- tryCatch(select_df_cv(mk, subj, w, 0.5, df_grid = cfg$df_grid), error = function(e) 4)
  taus <- cfg$taus
  fits <- lapply(taus, function(tt) tryCatch(fit_A(mk, w, df = df, tau = tt), error = function(e) NULL))
  if (any(sapply(fits, is.null))) return(NULL)
  # crossing on the grid x profiles, split interior/boundary, with positions
  cc <- list(int = c(cr = 0, pr = 0), bnd = c(cr = 0, pr = 0))
  mags <- numeric(0); pos_ut <- numeric(0); pos_tmu <- numeric(0)
  for (xi in seq_len(nrow(grd$Xs))) {
    Q <- sapply(fits, function(ft) predict_A(ft, grd$pts$u, grd$pts$t, grd$Xs$X1[xi], grd$Xs$X2[xi]))
    dQ <- Q[, -1] - Q[, -ncol(Q)]                       # adjacent diffs (>=0 if monotone)
    cr_row <- rowSums(dQ < 0)                            # #crossings per grid point
    for (tag in c("int", "bnd")) {
      sel <- if (tag == "int") grd$pts$interior else !grd$pts$interior
      cc[[tag]]["cr"] <- cc[[tag]]["cr"] + sum(dQ[sel, ] < 0)
      cc[[tag]]["pr"] <- cc[[tag]]["pr"] + sum(is.finite(dQ[sel, ]))
    }
    neg <- which(dQ < 0, arr.ind = TRUE)
    if (nrow(neg)) { mags <- c(mags, -dQ[dQ < 0])
      pos_ut <- c(pos_ut, grd$pts$u[neg[, 1]] / grd$pts$t[neg[, 1]])
      pos_tmu <- c(pos_tmu, grd$pts$t[neg[, 1]] - grd$pts$u[neg[, 1]]) }
  }
  # calibration on the fitting sample (in-sample; partly mechanical)
  Qrow <- sapply(fits, function(ft) predict_A(ft, mk$U, mk$T, mk$X1, mk$X2))   # nrow x ntau
  emp_before <- colMeans(mk$L <= Qrow)                       # empirical tau-coverage
  Qsort <- t(apply(Qrow, 1, sort))                           # monotone rearrangement
  emp_after  <- colMeans(mk$L <= Qsort)
  # OUT-OF-SAMPLE calibration on an INDEPENDENT test draw (#12,#15): a fresh data
  # set from the SAME DGP (so the test marker rows experience the same left
  # truncation, competing death, censoring, and irregular visit process as
  # training). Calibration is the mean over test MARKER ROWS (multi-visit subjects
  # weigh more, matching how the surface is used). Reported BOTH before AND after
  # monotone rearrangement -- in-sample invariance does not imply OOS invariance.
  dte <- tryCatch(generate_data(n, scenario = scenario, params = p), error = function(e) NULL)
  cal_oos_b <- NA; cal_oos_a <- NA; cal_oos_subj <- NA; n_test <- NA
  if (!is.null(dte)) { mkte <- prep_mk(dte)
    if (nrow(mkte) >= 30) {
      n_test <- length(unique(mkte$id))
      Qte <- sapply(fits, function(ft) predict_A(ft, mkte$U, mkte$T, mkte$X1, mkte$X2))
      cal_oos_b <- mean(abs(colMeans(mkte$L <= Qte) - taus))            # before, visit-wtd
      Qte_s <- t(apply(Qte, 1, sort))                                   # rearranged
      cal_oos_a <- mean(abs(colMeans(mkte$L <= Qte_s) - taus))          # after, visit-wtd
      # subject-EQUAL-weighted (#17): average the per-row indicators within each
      # test subject first, then over subjects, so multi-visit subjects do not
      # dominate. Reported alongside the visit-weighted figure.
      ind <- (mkte$L <= Qte) + 0; cnt <- as.numeric(table(mkte$id))
      subj_cov <- rowsum(ind, mkte$id) / cnt
      cal_oos_subj <- mean(abs(colMeans(subj_cov) - taus)) } }
  imse_before <- mean((predict_A(fits[[3]], grd$pts$u, grd$pts$t, 0, 0) - true_surface(grd$pts$u, grd$pts$t, 0.5, p))^2)
  c(cr_int = cc$int["cr"] / max(cc$int["pr"], 1), cr_bnd = cc$bnd["cr"] / max(cc$bnd["pr"], 1),
    any_cross = as.numeric(length(mags) > 0),
    max_mag = if (length(mags)) max(mags) else 0, mean_mag = if (length(mags)) mean(mags) else 0,
    pos_ut = if (length(pos_ut)) mean(pos_ut) else NA, pos_tmu = if (length(pos_tmu)) mean(pos_tmu) else NA,
    cal_before = mean(abs(emp_before - taus)), cal_after = mean(abs(emp_after - taus)),
    cal_oos_b = cal_oos_b, cal_oos_a = cal_oos_a, cal_oos_subj = cal_oos_subj,
    n_test = n_test, imse_med = imse_before)
}

run_cell <- function(n, scenario, cfg, grd) {
  reps <- Filter(Negate(is.null), mclapply(seq_len(cfg$R), one_rep, n = n, scenario = scenario,
            cfg = cfg, grd = grd, mc.cores = cfg$cores, mc.preschedule = FALSE))
  M <- do.call(rbind, reps); m <- function(c) mean(M[, c], na.rm = TRUE)
  data.frame(n = n, scenario = scenario, R = nrow(M),
    cross_rate_interior = m("cr_int.cr"), cross_rate_boundary = m("cr_bnd.cr"),
    p_reps_any_cross = m("any_cross"),
    mean_mag_per_incident = m("mean_mag"), max_mag = max(M[, "max_mag"], na.rm = TRUE),
    incident_mean_u_over_t = m("pos_ut"), incident_mean_t_minus_u = m("pos_tmu"),
    cal_err_before = m("cal_before"), cal_err_after_rearr = m("cal_after"),
    cal_err_oos_before = m("cal_oos_b"), cal_err_oos_after = m("cal_oos_a"),
    cal_err_oos_subj = m("cal_oos_subj"), n_test = round(m("n_test")),
    cross_rate_after_rearr = 0, row.names = NULL)
}

main <- function(cfg = CFG, out = "results/crossing.csv") {
  dir.create("results", showWarnings = FALSE); grd <- eval_grid(cfg); rows <- list(); k <- 1
  for (sc in cfg$scenarios) for (n in cfg$ns) {
    cat(sprintf("[%s] crossing sc=%d n=%d (dense)\n", format(Sys.time(), "%H:%M:%S"), sc, n)); flush.console()
    rows[[k]] <- run_cell(n, sc, cfg, grd); k <- k + 1
    write.csv(do.call(rbind, rows), out, row.names = FALSE)
  }
  cat("done ->", out, "\n")
}
if (sys.nframe() == 0) main()
