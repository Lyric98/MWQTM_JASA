# ============================================================================
# run_surface.R (sim-review #13,#15,#16,#17,#21) — the NEW-SCIENCE batch.
# Tests the quantile-TRAJECTORY-SURFACE selling point with a (u,t)-dependent
# scale sigma(u,t,X) (het_surface) so m_{tau,0}(u,t) is non-parallel across tau.
# Reports, for the default CC--A-CV estimator:
#   * beta_{tau,1} AND beta_{tau,2} bias / ESD / bootstrap coverage  (#17)
#   * surface IMSE split into interior vs boundary, + integrated bias (#16)
#   * mean estimated vs true surface grid for a figure                (#16)
#   * design diagnostics: n_T, ESS, IPCW weight quantiles, non-estimable
#     fraction, truncation/censoring/visit summaries                  (#21)
#   * a HARD stress paramset: high truncation + censoring + sparse visits (#15)
# Writes results/surface_*.csv and results/surface_grid.csv.
# ============================================================================
suppressWarnings(suppressMessages({library(parallel); source("R/run_one.R"); source("R/bootstrap.R")}))

# lean A-CV bootstrap SE for BOTH beta1 and beta2 (df fixed, full G_C refit)
boot_acv_b12 <- function(d, tau, df, B) {
  mk_by <- split(d$mk, d$mk$id); bm <- matrix(NA, B, 2)
  for (b in seq_len(B)) {
    db <- resample_subjects(d, mk_by); mkb <- prep_mk(db); sjb <- db$subj
    if (is.null(mkb) || nrow(mkb) < 30) next
    w <- attach_weights(mkb, sjb, ipcw_weights(sjb, GC_at(fit_GC(sjb))))
    bb <- tryCatch(fit_A(mkb, w, df = df, tau = tau)$beta, error = function(e) c(NA, NA))
    bm[b, ] <- as.numeric(bb)
  }
  apply(bm, 2, sd, na.rm = TRUE)
}

# parameter sets ------------------------------------------------------------
params_het <- function() {
  p <- default_params(); p$het_surface <- TRUE; p$sig1 <- 0.4
  p$sig2 <- -0.05; p$sig3 <- 0.35; p          # sigma in [~0.48, ~1.12]
}
params_stress <- function() {                  # hard truncation/censoring/sparse
  p <- params_het()
  p$A_lo <- -1.0; p$A_hi <- 9.0                # later/ wider entry -> more truncation
  p$rate_C <- 0.16                             # heavier loss to follow-up
  p$visit_gap <- 1.6                           # sparser visits (~2 per subject)
  p
}

CFG <- list(R = 500, R_cov = 300, B = 300, a_n = 0.05, df_grid = c(3, 4, 5),
  ns = c(500, 1000, 2000), taus = c(0.25, 0.50, 0.75),
  cores = { sc <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "")); if (is.na(sc) || sc < 1) 4L else sc },
  seed = 20260617,
  # evaluation grid for surface IMSE / figure (interior vs boundary split)
  tgrid = seq(1.5, 6.5, length.out = 14), nu = 14, edge = 0.4,
  # non-estimable masking (#8): a grid cell is estimable iff >= mask_min observed
  # pre-onset visits fall within +/- mask_h (in both u and t) of it.
  mask_h = 0.5, mask_min = 10, mask_min_subj = 5)

# surface evaluation points with interior/boundary tag
surf_points <- function(cfg) {
  pts <- do.call(rbind, lapply(cfg$tgrid, function(tt) {
    ug <- seq(0.25, tt - 0.25, length.out = cfg$nu)
    data.frame(u = ug, t = tt)
  }))
  e <- cfg$edge
  pts$interior <- with(pts, u > e & u < t - e & t > min(cfg$tgrid) + e & t < max(cfg$tgrid) - e)
  pts
}

# effective sample size of the IPCW weights over event subjects
ess_weights <- function(subj) {
  w <- ipcw_weights(subj, GC_at(fit_GC(subj))); w <- w[w > 0]
  c(ess = (sum(w))^2 / sum(w^2), n_T = length(w),
    wq50 = median(w), wq95 = quantile(w, .95, names = FALSE),
    wq99 = quantile(w, .99, names = FALSE), wmax = max(w))
}

one_rep <- function(r, n, tau, pfun, cfg, pts) {
  set.seed(cfg$seed + r * 7919 + round(tau * 1000))
  p <- pfun()
  d <- tryCatch(generate_data(n, scenario = 2, params = p), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  mk <- prep_mk(d); subj <- d$subj
  if (nrow(mk) < 30) return(NULL)
  w <- attach_weights(mk, subj, ipcw_weights(subj, GC_at(fit_GC(subj))))
  df <- tryCatch(select_df_cv(mk, subj, w, tau, df_grid = cfg$df_grid), error = function(e) 4)
  fit <- tryCatch(fit_A(mk, w, df = df, tau = tau), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  # surface error on the evaluation grid
  mhat  <- predict_A(fit, pts$u, pts$t, 0, 0)
  mtrue <- true_surface(pts$u, pts$t, tau, p)
  err <- mhat - mtrue
  # local information per grid cell + non-estimable mask (#8)
  hw <- cfg$mask_h
  loc <- vapply(seq_len(nrow(pts)), function(j)
    sum(abs(mk$U - pts$u[j]) <= hw & abs(mk$T - pts$t[j]) <= hw), numeric(1))
  est <- loc >= cfg$mask_min                      # estimable by VISIT count
  imse_masked <- if (any(est)) mean(err[est]^2) else NA
  # #13: distinct-SUBJECT local support (visits cluster within subject, so a
  # 10-visit window may come from only 1-2 onset subjects). Count distinct ids
  # and mask on that -- the surface-relevant notion of local estimability.
  loc_subj <- vapply(seq_len(nrow(pts)), function(j) {
    sel <- abs(mk$U - pts$u[j]) <= hw & abs(mk$T - pts$t[j]) <= hw
    length(unique(mk$id[sel])) }, numeric(1))
  est_subj <- loc_subj >= cfg$mask_min_subj
  imse_masked_subj <- if (any(est_subj)) mean(err[est_subj]^2) else NA
  # bootstrap coverage subset (A-CV full refit, B resamples): SE of beta1 & beta2
  bse <- c(NA, NA)
  if (r <= cfg$R_cov)
    bse <- tryCatch(boot_acv_b12(d, tau, df, cfg$B), error = function(e) c(NA, NA))
  list(beta1 = fit$beta[1], beta2 = fit$beta[2], df = df,
       imse_int = mean(err[pts$interior]^2), imse_bnd = mean(err[!pts$interior]^2),
       imse_unmasked = mean(err^2), imse_masked = imse_masked, frac_nonest = mean(!est),
       imse_masked_subj = imse_masked_subj, frac_nonest_subj = mean(!est_subj),
       mean_nsubj = mean(loc_subj),
       ibias_int = mean(err[pts$interior]),
       se1 = bse[1], se2 = bse[2], mhat = mhat,
       diag = if (r == 1) ess_weights(subj) else NULL,
       summ = if (r == 1) d$summary else NULL)
}

run_cell <- function(n, tau, pfun, label, cfg, pts) {
  p0 <- pfun(); tb <- true_beta(tau, p0)
  reps <- Filter(Negate(is.null),
                 mclapply(seq_len(cfg$R), one_rep, n = n, tau = tau, pfun = pfun,
                          cfg = cfg, pts = pts, mc.cores = cfg$cores, mc.preschedule = FALSE))
  g <- function(f) sapply(reps, function(x) { v <- x[[f]]; if (is.null(v)) NA else v })
  b1 <- g("beta1"); b2 <- g("beta2"); se1 <- g("se1"); se2 <- g("se2")
  cov1 <- mean(abs(b1 - tb["beta1"]) <= 1.96 * se1, na.rm = TRUE)
  cov2 <- mean(abs(b2 - tb["beta2"]) <= 1.96 * se2, na.rm = TRUE)
  # mean estimated surface across reps (for the figure)
  M <- do.call(rbind, lapply(reps, function(x) x$mhat))
  mbar <- colMeans(M, na.rm = TRUE)
  diag <- reps[[1]]$diag; summ <- reps[[1]]$summ
  list(row = data.frame(
    paramset = label, n = n, tau = tau, R = length(reps),
    true_b1 = tb["beta1"], true_b2 = tb["beta2"],
    bias_b1 = mean(b1, na.rm = TRUE) - tb["beta1"], esd_b1 = sd(b1, na.rm = TRUE),
    bias_b2 = mean(b2, na.rm = TRUE) - tb["beta2"], esd_b2 = sd(b2, na.rm = TRUE),
    cov_b1_boot = cov1, cov_b2_boot = cov2, R_cov = sum(is.finite(se1)),
    imse_int = mean(g("imse_int"), na.rm = TRUE), imse_bnd = mean(g("imse_bnd"), na.rm = TRUE),
    imse_unmasked = mean(g("imse_unmasked"), na.rm = TRUE),
    imse_masked = mean(g("imse_masked"), na.rm = TRUE),
    frac_nonest = mean(g("frac_nonest"), na.rm = TRUE),
    imse_masked_subj = mean(g("imse_masked_subj"), na.rm = TRUE),
    frac_nonest_subj = mean(g("frac_nonest_subj"), na.rm = TRUE),
    mean_nsubj = mean(g("mean_nsubj"), na.rm = TRUE),
    ibias_int = mean(g("ibias_int"), na.rm = TRUE),
    n_T = if (!is.null(diag)) diag["n_T"] else NA,
    ess = if (!is.null(diag)) diag["ess"] else NA,
    w_med = if (!is.null(diag)) diag["wq50"] else NA,
    w_p95 = if (!is.null(diag)) diag["wq95"] else NA,
    w_p99 = if (!is.null(diag)) diag["wq99"] else NA,
    w_max = if (!is.null(diag)) diag["wmax"] else NA,
    trunc_frac = if (!is.null(summ)) summ["trunc_frac"] else NA,
    p_cens = if (!is.null(summ)) summ["p_cens"] else NA,
    mean_visits = if (!is.null(summ)) summ["mean_visits_event"] else NA,
    row.names = NULL),
    grid = data.frame(paramset = label, n = n, tau = tau,
                      u = pts$u, t = pts$t, interior = pts$interior,
                      m_true = true_surface(pts$u, pts$t, tau, p0), m_hat = mbar))
}

main <- function(cfg = CFG) {
  dir.create("results", showWarnings = FALSE)
  pts <- surf_points(cfg)
  sets <- list(het = params_het, stress = params_stress)
  rows <- list(); grids <- list(); k <- 1
  for (label in names(sets)) for (n in cfg$ns) for (tau in cfg$taus) {
    cat(sprintf("[%s] surface %s n=%d tau=%.2f\n", format(Sys.time(), "%H:%M:%S"), label, n, tau)); flush.console()
    rc <- run_cell(n, tau, sets[[label]], label, cfg, pts)
    rows[[k]] <- rc$row; grids[[k]] <- rc$grid; k <- k + 1
    write.csv(do.call(rbind, rows), "results/surface.csv", row.names = FALSE)
    write.csv(do.call(rbind, grids), "results/surface_grid.csv", row.names = FALSE)
  }
  cat("done -> results/surface.csv, results/surface_grid.csv\n")
}
if (sys.nframe() == 0) main()
