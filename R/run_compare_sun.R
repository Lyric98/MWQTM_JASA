# ============================================================================
# run_compare_sun.R (compare.md, Phase 1) — Sun conditional-MEAN vs our
# conditional-QUANTILE (+ Q-to-Mean integration), three DGPs, same Monte-Carlo
# datasets. Point estimates + Monte-Carlo ESD + masked IMSE on a common interior
# domain; PAIRED Sun-vs-Q mean-contrast IMSE (the DGP3 headline). SLURM array:
# one (dgp, n) cell per task. Writes results/compare_sun_cell_XX.csv.
# ============================================================================
suppressWarnings(suppressMessages({library(parallel); source("R/compare_sun.R")}))

CFG <- list(
  dgps = 1:3, ns = c(500, 1000), R = 500, regime = "A", df = 4,
  tau = c(.02,.05,.1,.2,.25,.3,.4,.5,.6,.7,.75,.8,.9,.95,.98),
  taus_report = c(.25, .5, .75),
  h_mult = 1.0,                         # Sun bandwidth (sensitivity 0.8/1.2 separate)
  # DGP knobs (smoke-tuned so DGP3 mean-contrast clearly varies)
  sig1_cmp = 0.5, a1 = -0.10, a2 = 0.70,
  cores = { sc <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "")); if (is.na(sc) || sc < 1) 4L else sc },
  seed = 20260621)

CELLS <- expand.grid(dgp = 1:3, n = c(500, 1000))

logimse <- function(est, tru, mask) { e <- (log(est) - log(tru))[mask]; mean(e^2, na.rm = TRUE) }

one_rep <- function(r, n, dgp, cfg, pts) {
  set.seed(cfg$seed + r * 7919 + dgp * 101)
  p <- default_params(); p$sig1_cmp <- cfg$sig1_cmp; p$a1 <- cfg$a1; p$a2 <- cfg$a2
  d <- tryCatch(gen_compare(n, dgp, p, regime = cfg$regime), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  mk <- cc_rows(d); pp <- d$params
  if (nrow(mk) < 50) return(NULL)
  sf <- tryCatch(sun_fit(mk, cfg$h_mult), error = function(e) NULL)
  qf <- tryCatch(qmean_fit(mk, cfg$tau, df = cfg$df), error = function(e) NULL)
  if (is.null(sf) || is.null(qf) || any(!is.finite(sf$beta))) return(NULL)
  mask <- supp_mask(pts, mk); if (sum(mask) < 5) return(NULL)
  # quantile slopes at the reported taus
  qb <- qmean_betatau(qf); ti <- match(cfg$taus_report, cfg$tau)
  b_q <- qb[ti, 1]                                            # Qhat beta_tau1
  # mean surfaces (x=0,1) + contrast on the common mask
  trueM0 <- true_mean(pts$u, pts$t, 0, 0, pp); trueM1 <- true_mean(pts$u, pts$t, 1, 0, pp)
  trueD  <- true_contrast(pts$u, pts$t, pp)
  sunM0 <- sun_predict(sf, pts$u, pts$t, 0, 0); sunM1 <- sun_predict(sf, pts$u, pts$t, 1, 0)
  qM0 <- qmean_predict(qf, pts$u, pts$t, 0, 0, c(0, 1)); qM1 <- qmean_predict(qf, pts$u, pts$t, 1, 0, c(0, 1))
  c(beta_sun1 = sf$beta[1], beta_sun2 = sf$beta[2],
    bq_25 = b_q[1], bq_50 = b_q[2], bq_75 = b_q[3],
    logimse_sun0 = logimse(sunM0, trueM0, mask), logimse_q0 = logimse(qM0, trueM0, mask),
    logimse_sun1 = logimse(sunM1, trueM1, mask), logimse_q1 = logimse(qM1, trueM1, mask),
    ctr_sun = imse(log(sunM1) - log(sunM0), trueD, mask),
    ctr_q   = imse(log(qM1) - log(qM0), trueD, mask),
    nmask = sum(mask), ncc = nrow(mk))
}

run_cell <- function(n, dgp, cfg) {
  p <- default_params(); p$sig1_cmp <- cfg$sig1_cmp; p$a1 <- cfg$a1; p$a2 <- cfg$a2; p$dgp <- dgp
  pts <- eval_grid()
  tb <- sapply(cfg$taus_report, function(tt) true_beta_tau(tt, p)[1])   # true beta_tau1
  M <- do.call(rbind, Filter(Negate(is.null),
        mclapply(seq_len(cfg$R), one_rep, n = n, dgp = dgp, cfg = cfg, pts = pts,
                 mc.cores = cfg$cores, mc.preschedule = FALSE)))
  g <- function(col) M[, col]
  mm <- function(col) mean(M[, col], na.rm = TRUE); ss <- function(col) sd(M[, col], na.rm = TRUE)
  dctr <- M[, "ctr_q"] - M[, "ctr_sun"]                                 # paired (Q - Sun) contrast IMSE
  data.frame(dgp = dgp, n = n, R = nrow(M), regime = cfg$regime,
    # quantile-slope recovery (bias vs true beta_tau1)
    true_b25 = tb[1], true_b50 = tb[2], true_b75 = tb[3],
    bias_bq25 = mm("bq_25") - tb[1], bias_bq50 = mm("bq_50") - tb[2], bias_bq75 = mm("bq_75") - tb[3],
    esd_bq50 = ss("bq_50"), spread_q = mm("bq_75") - mm("bq_25"), spread_true = tb[3] - tb[1],
    # Sun scalar slope (pseudo-effect in DGP3)
    beta_sun1 = mm("beta_sun1"), esd_sun1 = ss("beta_sun1"),
    # log-mean surface IMSE (x=0,1): Sun vs Q-to-full-mean
    logimse_sun0 = mm("logimse_sun0"), logimse_q0 = mm("logimse_q0"),
    logimse_sun1 = mm("logimse_sun1"), logimse_q1 = mm("logimse_q1"),
    # HEADLINE: mean-contrast IMSE, paired
    ctr_sun = mm("ctr_sun"), ctr_q = mm("ctr_q"),
    dctr_mean = mean(dctr, na.rm = TRUE), dctr_mcse = sd(dctr, na.rm = TRUE) / sqrt(sum(is.finite(dctr))),
    nmask = mm("nmask"), ncc = mm("ncc"), row.names = NULL)
}

main <- function(cfg = CFG) {
  dir.create("results", showWarnings = FALSE)
  task <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", ""))
  if (!is.na(task)) {
    cell <- CELLS[task, ]
    cat(sprintf("[%s] compare_sun cell %d: dgp=%d n=%d\n", format(Sys.time(), "%H:%M:%S"), task, cell$dgp, cell$n)); flush.console()
    row <- run_cell(cell$n, cell$dgp, cfg)
    write.csv(row, sprintf("results/compare_sun_cell_%02d.csv", task), row.names = FALSE)
    cat("done cell", task, "\n")
  } else {
    rows <- list()
    for (i in seq_len(nrow(CELLS))) {
      cat(sprintf("[%s] cell %d dgp=%d n=%d\n", format(Sys.time(), "%H:%M:%S"), i, CELLS$dgp[i], CELLS$n[i])); flush.console()
      rows[[i]] <- run_cell(CELLS$n[i], CELLS$dgp[i], cfg)
      write.csv(do.call(rbind, rows), "results/compare_sun.csv", row.names = FALSE)
    }
    cat("done -> results/compare_sun.csv\n")
  }
}
if (sys.nframe() == 0) main()
