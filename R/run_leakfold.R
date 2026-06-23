# ============================================================================
# run_leakfold.R (sim-review round-12 #2) — establish CAUSALLY that the
# duplicate-ID fold-leakage FIX (not fold randomization) is what widened the
# B-OS-CF bootstrap SE. On the SAME datasets and SAME bootstrap resamples we
# compute the cross-fitted one-step two ways per resample:
#   * LEAKED : fold by the resampled id, so duplicate copies of an original
#              subject split across folds (the old bug; train/validation leak).
#   * CLEAN  : fold by orig_id, duplicates kept together (the fix).
# The point estimate is the clean B-OS-CF on the original data (no duplicates, so
# leak is irrelevant there). Reports, per cell, the bootstrap SER and coverage
# under each folding. If SER_leaked < SER_clean, the leak artificially tightened
# the SE and the fix is the cause of the (honest) widening. CC estimator, Sc 2.
# Writes results/leakfold.csv.
# ============================================================================
suppressWarnings(suppressMessages({library(parallel); source("R/cf.R"); source("R/bootstrap.R")}))

CFG <- list(ns = c(500, 1000, 2000), taus = c(0.25, 0.50, 0.75), scenario = 2,
  R = 250, B = 150, K = 5, a_n = 0.05, df_grid = c(3, 4, 5), df_g = 4,
  cores = { sc <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "")); if (is.na(sc) || sc < 1) 4L else sc },
  seed = 20260616)

cf_b1 <- function(d, tau, df, cfg, leak) {
  z <- one_step_cf(d, tau, df_fixed = df, type = "cc", K = cfg$K, a_n = cfg$a_n,
                   df_g = cfg$df_g, leak = leak)
  if (is.null(z)) NA_real_ else z$beta[1]
}

one_rep <- function(r, n, tau, cfg) {
  set.seed(cfg$seed + r * 7919 + cfg$scenario * 101 + round(tau * 1000))
  p <- default_params(); p$sig1 <- 0.4
  d <- tryCatch(generate_data(n, scenario = cfg$scenario, params = p), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  mk <- prep_mk(d); subj <- d$subj; if (nrow(mk) < 30) return(NULL)
  w_ip <- attach_weights(mk, subj, ipcw_weights(subj, GC_at(fit_GC(subj))))
  df <- tryCatch(select_df_cv(mk, subj, w_ip, tau, df_grid = cfg$df_grid), error = function(e) 4)
  b_pt <- cf_b1(d, tau, df, cfg, leak = FALSE)            # clean point (original data)
  if (!is.finite(b_pt)) return(NULL)
  mk_by <- split(d$mk, d$mk$id); vl <- rep(NA_real_, cfg$B); vc <- rep(NA_real_, cfg$B)
  for (b in seq_len(cfg$B)) {
    db <- resample_subjects(d, mk_by); mkb <- prep_mk(db); sjb <- db$subj
    if (is.null(mkb) || nrow(mkb) < 30) next
    vl[b] <- cf_b1(db, tau, df, cfg, leak = TRUE)          # leaked folding (old bug)
    vc[b] <- cf_b1(db, tau, df, cfg, leak = FALSE)         # clean folding (fix)
  }
  c(b_pt = b_pt, se_leaked = sd(vl, na.rm = TRUE), se_clean = sd(vc, na.rm = TRUE),
    nbl = sum(is.finite(vl)), nbc = sum(is.finite(vc)))
}

run_cell <- function(n, tau, cfg) {
  tb <- as.numeric(true_beta(tau, { p <- default_params(); p$sig1 <- 0.4; p })["beta1"])
  M <- do.call(rbind, Filter(Negate(is.null),
        mclapply(seq_len(cfg$R), one_rep, n = n, tau = tau, cfg = cfg,
                 mc.cores = cfg$cores, mc.preschedule = FALSE)))
  bp <- M[, "b_pt"]; sl <- M[, "se_leaked"]; sc <- M[, "se_clean"]
  esd <- sd(bp, na.rm = TRUE)
  cov <- function(s) mean(abs(bp - tb) <= 1.96 * s, na.rm = TRUE)
  rr <- (sc / sl); rr <- rr[is.finite(rr)]                          # paired per-rep clean/leaked SE ratio
  data.frame(scenario = cfg$scenario, n = n, tau = tau, R = nrow(M), esd = esd,
    se_leaked = mean(sl, na.rm = TRUE), se_clean = mean(sc, na.rm = TRUE),
    ser_leaked = mean(sl, na.rm = TRUE) / esd, ser_clean = mean(sc, na.rm = TRUE) / esd,
    se_ratio = mean(rr), se_ratio_mcse = sd(rr) / sqrt(length(rr)),  # paired clean/leaked ratio +/- MCSE
    cov_leaked = cov(sl), cov_clean = cov(sc), row.names = NULL)
}

main <- function(cfg = CFG, out = "results/leakfold.csv") {
  dir.create("results", showWarnings = FALSE); rows <- list(); k <- 1
  for (n in cfg$ns) for (tau in cfg$taus) {
    cat(sprintf("[%s] leakfold n=%d tau=%.2f\n", format(Sys.time(), "%H:%M:%S"), n, tau)); flush.console()
    rows[[k]] <- run_cell(n, tau, cfg); k <- k + 1
    write.csv(do.call(rbind, rows), out, row.names = FALSE)
  }
  cat("done ->", out, "\n")
}
if (sys.nframe() == 0) main()
