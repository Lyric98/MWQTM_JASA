# ============================================================================
# run_compare2.R — Sun (real C++) vs ours, on Sun's design. SLURM array: one
# (dgp, n) cell per task (9 cells: dgp 1:3 x n in {500,1000,2000}). NO bootstrap
# (point estimates + Monte-Carlo ESD/IMSE only). Per cell writes:
#   results/compare2_cell_XX.csv   scalar summary row
#   results/compare2_fig_XX.rds    averaged beta_tau curve + mean/contrast
#                                  surfaces over the common grid (for the figures)
# ============================================================================
suppressWarnings(suppressMessages({library(parallel); source("R/compare2.R")}))

CFG <- list(
  R = 500,
  taus = c(.05,.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.75,.8,.85,.9,.95),
  report = c(.1,.25,.5,.75,.9), df = 4,
  # Sun mu_0 surface bandwidth: HALF the density-reference rule (2.34 sd(A) N^{-1/6}).
  # The 2.34 constant is Silverman's DENSITY rule; on the curved regression surface
  # mu_2 it over-smooths (over- at low u, under- at high u), unfairly biasing Sun's
  # baseline. Undersmoothing to ~the regression-reference scale makes his kernel
  # approximately unbiased. beta/contrast are kernel-free and unaffected.
  h_mu_mult = 0.5,
  cores = { sc <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK","")); if (is.na(sc)||sc<1) 4L else sc },
  seed = 20260621)
CELLS <- expand.grid(dgp = 1:5, n = c(500, 1000, 2000))

logimse <- function(logest, tru, mask) mean(((logest) - log(tru))[mask]^2, na.rm = TRUE)

one_rep2 <- function(r, n, dgp, cfg, pts) {
  set.seed(cfg$seed + r * 7919 + dgp * 101)
  p <- cmp2_params()
  d  <- tryCatch(gen_sun(n, dgp, p), error = function(e) NULL); if (is.null(d)) return(NULL)
  mk <- cc_rows2(d); if (nrow(mk) < 60) return(NULL)
  sb <- tryCatch(sun_beta2(d), error = function(e) NULL)
  qf <- tryCatch(qfit2(mk, cfg$taus, df = cfg$df), error = function(e) NULL)
  if (is.null(sb) || is.null(qf) || any(!is.finite(sb$beta)) || sum(qf$ok) < 6) return(NULL)
  mask <- supp_mask2(pts, mk); if (sum(mask) < 8) return(NULL)

  qb <- qbeta2(qf)[, 1]                                       # Qhat beta_tau1 over full grid
  ri <- match(cfg$report, cfg$taus)
  # ---- mean & contrast surfaces on the common mask ----
  trueM0 <- true_mean2(pts$u, pts$t, 0, 0.5, p, dgp); trueM1 <- true_mean2(pts$u, pts$t, 1, 0.5, p, dgp)
  trueD  <- true_contrast2(pts$u, pts$t, p, dgp)
  # two mean reconstructions from the fitted quantile process:
  qLM0 <- qmean_ln(qf, pts$u, pts$t, 0, 0.5);  qLM1 <- qmean_ln(qf, pts$u, pts$t, 1, 0.5)  # Q-LN (log-loc-scale)
  qI0  <- qmean_int(qf, pts$u, pts$t, 0, 0.5); qI1  <- qmean_int(qf, pts$u, pts$t, 1, 0.5) # Q-INT (true integral)
  # Sun surface at his density-rule bandwidth (h1) AND the regression bandwidth (h.5):
  h1 <- 2.34 * sd(d$dl$A) * nrow(d$aw)^(-1/6); hr <- cfg$h_mu_mult * h1
  sM0 <- sun_mean2(d, sb, pts$u, pts$t, 0, 0.5, h_mu = hr); sM1 <- sun_mean2(d, sb, pts$u, pts$t, 1, 0.5, h_mu = hr)
  sM0h1 <- sun_mean2(d, sb, pts$u, pts$t, 0, 0.5, h_mu = h1); sM1h1 <- sun_mean2(d, sb, pts$u, pts$t, 1, 0.5, h_mu = h1)
  qD <- qLM1 - qLM0; qDi <- qI1 - qI0; sD <- rep(sb$beta[1], nrow(pts))   # contrasts (Sun = beta1, const)
  # Sun Sec.3.5 stratified-baseline contrast (nonparametric per X1 group):
  ssD <- tryCatch(sun_strat_contrast(d, pts$u, pts$t, 0.5, cfg$h_mu_mult),
                  error = function(e) rep(NA_real_, nrow(pts)))
  list(
    beta_sun = sb$beta, qb_report = qb[ri], qb_curve = qb,
    spread_q = qb[match(.75, cfg$taus)] - qb[match(.25, cfg$taus)],
    ctr_sun = imse2(sD, trueD, mask), ctr_sstrat = imse2(ssD, trueD, mask),
    ctr_q = imse2(qD, trueD, mask), ctr_qint = imse2(qDi, trueD, mask),
    lim_sun0 = logimse(log(sM0), trueM0, mask), lim_q0 = logimse(qLM0, trueM0, mask), lim_qint0 = logimse(qI0, trueM0, mask),
    lim_sun1 = logimse(log(sM1), trueM1, mask), lim_q1 = logimse(qLM1, trueM1, mask), lim_qint1 = logimse(qI1, trueM1, mask),
    lim_sun0_h1 = logimse(log(sM0h1), trueM0, mask), lim_sun1_h1 = logimse(log(sM1h1), trueM1, mask),  # Sun at density bandwidth
    maxY = max(d$dl$Y),                                          # heavy-tail diagnostic (DGP5)
    # surfaces (masked NA outside) for figure averaging
    sfc = list(qD = ifelse(mask, qD, NA), sD = ifelse(mask, sD, NA), ssD = ifelse(mask, ssD, NA),
               qLM0 = ifelse(mask, qLM0, NA), sLM0 = ifelse(mask, log(sM0), NA),
               qLM1 = ifelse(mask, qLM1, NA), sLM1 = ifelse(mask, log(sM1), NA)),
    nmask = sum(mask), ncc = nrow(mk), event = unname(d$frac["event"]))
}

run_cell2 <- function(n, dgp, cfg) {
  p <- cmp2_params(); pts <- eval_grid2()
  tb <- sapply(cfg$report, function(tt) true_betatau2(tt, p, dgp)[1])
  res <- Filter(Negate(is.null),
    mclapply(seq_len(cfg$R), one_rep2, n = n, dgp = dgp, cfg = cfg, pts = pts,
             mc.cores = cfg$cores, mc.preschedule = FALSE))
  R <- length(res)
  col  <- function(f) sapply(res, f)
  mat  <- function(f) do.call(rbind, lapply(res, f))                  # R x .
  bsun <- mat(function(x) x$beta_sun); qbr <- mat(function(x) x$qb_report)
  dctr <- col(function(x) x$ctr_q) - col(function(x) x$ctr_sun)       # paired Q-LN - Sun-basic
  dctr_i <- col(function(x) x$ctr_qint) - col(function(x) x$ctr_sun)  # paired Q-INT - Sun-basic
  am   <- function(M) colMeans(M, na.rm = TRUE)                       # avg surface over reps
  # averaged figure objects
  fig <- list(dgp = dgp, n = n, R = R, pts = pts, taus = cfg$taus, report = cfg$report,
    true_betatau = sapply(cfg$taus, function(tt) true_betatau2(tt, p, dgp)[1]),
    qb_curve_mean = am(mat(function(x) x$qb_curve)),
    qb_curve_sd   = apply(mat(function(x) x$qb_curve), 2, sd, na.rm = TRUE),
    beta_sun_mean = mean(bsun[, 1]),
    trueD  = true_contrast2(pts$u, pts$t, p, dgp),
    trueM0 = true_mean2(pts$u, pts$t, 0, 0.5, p, dgp),
    trueM1 = true_mean2(pts$u, pts$t, 1, 0.5, p, dgp),
    qD = am(mat(function(x) x$sfc$qD)), sD = am(mat(function(x) x$sfc$sD)), ssD = am(mat(function(x) x$sfc$ssD)),
    qLM0 = am(mat(function(x) x$sfc$qLM0)), sLM0 = am(mat(function(x) x$sfc$sLM0)),
    qLM1 = am(mat(function(x) x$sfc$qLM1)), sLM1 = am(mat(function(x) x$sfc$sLM1)),
    beta_sun1_reps = bsun[, 1], qb50_reps = qbr[, match(0.5, cfg$report)])  # per-rep (DGP5 robustness fig)
  row <- data.frame(dgp = dgp, n = n, R = R,
    event = mean(col(function(x) x$event)), ncc = mean(col(function(x) x$ncc)),
    nmask = mean(col(function(x) x$nmask)),
    beta_sun1 = mean(bsun[, 1]), esd_sun1 = sd(bsun[, 1]),
    beta_sun2 = mean(bsun[, 2]),
    # beta_tau1 bias at reported taus + recovered spread
    setNames(as.list(colMeans(qbr) - tb), paste0("bias_bq", sub("0\\.", "", cfg$report))),
    true_b10 = tb[1], true_b25 = tb[2], true_b50 = tb[3], true_b75 = tb[4], true_b90 = tb[5],
    spread_q = mean(col(function(x) x$spread_q)), spread_true = tb[4] - tb[2],
    spread90_q = mean(qbr[, 5]) - mean(qbr[, 1]), spread90_true = tb[5] - tb[1],
    esd_bq50 = sd(qbr[, 3]),
    # robustness diagnostics for Sun's slope (DGP5): median/IQR/upper quantiles/trimmed ESD + marker tail
    med_sun1 = median(bsun[, 1]), iqr_sun1 = IQR(bsun[, 1]),
    q95_sun1 = quantile(bsun[, 1], .95, names = FALSE), q99_sun1 = quantile(bsun[, 1], .99, names = FALSE),
    tesd_sun1 = sd(bsun[abs(bsun[, 1] - median(bsun[, 1])) <= 3 * IQR(bsun[, 1]), 1]),  # 3-IQR trimmed ESD
    maxY = mean(col(function(x) x$maxY)), maxY_p99 = quantile(col(function(x) x$maxY), .99, names = FALSE),
    # HEADLINE: mean-contrast IMSE (paired) -- Sun-basic, Sun-stratified, both Q reconstructions
    ctr_sun = mean(col(function(x) x$ctr_sun)), ctr_sstrat = mean(col(function(x) x$ctr_sstrat)),
    ctr_q = mean(col(function(x) x$ctr_q)), ctr_qint = mean(col(function(x) x$ctr_qint)),
    dctr_mean = mean(dctr), dctr_mcse = sd(dctr) / sqrt(R),
    dctr_int_mean = mean(dctr_i), dctr_int_mcse = sd(dctr_i) / sqrt(R),
    # log-mean SURFACE IMSE: Sun (regression bw), Sun (density bw h1), Q-LN, Q-INT
    lim_sun0 = mean(col(function(x) x$lim_sun0)), lim_q0 = mean(col(function(x) x$lim_q0)),
    lim_qint0 = mean(col(function(x) x$lim_qint0)), lim_sun0_h1 = mean(col(function(x) x$lim_sun0_h1)),
    lim_sun1 = mean(col(function(x) x$lim_sun1)), lim_q1 = mean(col(function(x) x$lim_q1)),
    lim_qint1 = mean(col(function(x) x$lim_qint1)), lim_sun1_h1 = mean(col(function(x) x$lim_sun1_h1)),
    row.names = NULL)
  list(row = row, fig = fig)
}

main <- function(cfg = CFG) {
  dir.create("results", showWarnings = FALSE)
  cat("Sun real C++ compiled:", try_compile_sun_cpp(), "\n"); flush.console()  # master compiles; forks inherit
  task <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", ""))
  cells <- if (!is.na(task)) CELLS[task, , drop = FALSE] else CELLS
  for (i in seq_len(nrow(cells))) {
    k <- if (!is.na(task)) task else i
    cat(sprintf("[%s] cell %d: dgp=%d n=%d\n", format(Sys.time(),"%H:%M:%S"),
                k, cells$dgp[i], cells$n[i])); flush.console()
    out <- run_cell2(cells$n[i], cells$dgp[i], cfg)
    write.csv(out$row, sprintf("results/compare2_cell_%02d.csv", k), row.names = FALSE)
    saveRDS(out$fig, sprintf("results/compare2_fig_%02d.rds", k))
    cat("  done cell", k, " R=", out$row$R, "\n"); flush.console()
  }
}
if (sys.nframe() == 0) main()
