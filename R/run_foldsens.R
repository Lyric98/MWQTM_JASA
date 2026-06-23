# ============================================================================
# run_foldsens.R (sim-review round-10 #4/#7) — does the B-OS-CF bootstrap's
# fold REDRAWING inflate the SE relative to the single realized split used for
# the point estimate? The headline B-OS-CF is mildly conservative (SER 1.04-1.13,
# coverage <=0.986) after the leakage fix; this isolates how much of that is
# fold-randomization variance vs the leakage fix itself.
#
# For CC B-OS-CF (Scenario 2, central/quartile tau):
#   POINT estimators
#     b_real : one realized fold draw (seed=1), exactly the headline point
#     b_avg  : average of M independent fold draws (the fold-marginal estimand)
#   BOOTSTRAP SE
#     se_redrawn : folds redrawn each resample (headline default)
#     se_fixed   : each ORIGINAL subject keeps its realized fold across resamples
#                  (fold_map), so the SE conditions on the realized split
#   COVERAGE: (b_real, se_fixed), (b_real, se_redrawn), (b_avg, se_redrawn).
# Hypothesis: se_redrawn > se_fixed, and (b_avg, se_redrawn) is best calibrated
# -> the conservatism is fold-randomization variance, now made explicit and the
# inference target (conditional-on-split vs fold-marginal) named.
# Writes results/foldsens.csv.
# ============================================================================
suppressWarnings(suppressMessages({library(parallel); source("R/cf.R"); source("R/bootstrap.R")}))

CFG <- list(R = 250, B = 120, M = 8, K = 5, a_n = 0.05, df_grid = c(3, 4, 5), df_g = 4,
  cells = list(c(2, 500, 0.25), c(2, 500, 0.50), c(2, 500, 0.75),
               c(2, 1000, 0.50), c(2, 2000, 0.50)),
  cores = { sc <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "")); if (is.na(sc) || sc < 1) 4L else sc },
  seed = 20260616)

cc_w <- function(mk, subj) attach_weights(mk, subj, as.numeric(subj$dT == 1))

# canonical fold map for the ORIGINAL sample, reproducing cf.R's seed=1 fold draw
make_fold_map <- function(subj, K, seed = 1) {
  cid <- if ("orig_id" %in% names(subj)) subj$orig_id else subj$id
  ucl <- unique(cid)
  old <- if (exists(".Random.seed", .GlobalEnv)) get(".Random.seed", .GlobalEnv) else NULL
  set.seed(seed); fm <- sample(rep_len(seq_len(K), length(ucl)))
  if (!is.null(old)) assign(".Random.seed", old, .GlobalEnv)
  names(fm) <- as.character(ucl); fm
}

cf_b1 <- function(d, tau, df, cfg, seed = 1, fold_map = NULL) {
  z <- one_step_cf(d, tau, df_fixed = df, type = "cc", K = cfg$K, a_n = cfg$a_n,
                   df_g = cfg$df_g, seed = seed, fold_map = fold_map)
  if (is.null(z)) NA_real_ else z$beta[1]
}

one_rep <- function(r, scenario, n, tau, cfg) {
  set.seed(cfg$seed + r * 7919 + scenario * 101 + round(tau * 1000))
  p <- default_params(); p$sig1 <- if (scenario == 1) 0 else 0.4
  d <- tryCatch(generate_data(n, scenario = scenario, params = p), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  mk <- prep_mk(d); subj <- d$subj; if (nrow(mk) < 30) return(NULL)
  w_ip <- attach_weights(mk, subj, ipcw_weights(subj, GC_at(fit_GC(subj))))
  df <- tryCatch(select_df_cv(mk, subj, w_ip, tau, df_grid = cfg$df_grid), error = function(e) 4)
  b_real <- cf_b1(d, tau, df, cfg, seed = 1)                  # realized split (headline)
  if (!is.finite(b_real)) return(NULL)
  b_avg <- mean(vapply(seq_len(cfg$M), function(s) cf_b1(d, tau, df, cfg, seed = 100 + s),
                       numeric(1)), na.rm = TRUE)
  fm0 <- make_fold_map(subj, cfg$K, seed = 1)                # the realized fold partition
  mk_by <- split(d$mk, d$mk$id)
  vr <- rep(NA_real_, cfg$B); vf <- rep(NA_real_, cfg$B)
  for (b in seq_len(cfg$B)) {
    db <- resample_subjects(d, mk_by); mkb <- prep_mk(db); sjb <- db$subj
    if (is.null(mkb) || nrow(mkb) < 30) next
    vr[b] <- cf_b1(db, tau, df, cfg, seed = 1)               # redrawn (default)
    vf[b] <- cf_b1(db, tau, df, cfg, fold_map = fm0)         # fixed by original fold
  }
  c(b_real = b_real, b_avg = b_avg,
    se_redrawn = sd(vr, na.rm = TRUE), se_fixed = sd(vf, na.rm = TRUE),
    nbr = sum(is.finite(vr)), nbf = sum(is.finite(vf)))
}

run_cell <- function(scenario, n, tau, cfg) {
  tb <- as.numeric(true_beta(tau, { p <- default_params(); p$sig1 <- if (scenario == 1) 0 else 0.4; p })["beta1"])
  M <- do.call(rbind, mclapply(seq_len(cfg$R), one_rep, scenario = scenario, n = n, tau = tau,
                               cfg = cfg, mc.cores = cfg$cores, mc.preschedule = FALSE))
  br <- M[, "b_real"]; ba <- M[, "b_avg"]; sr <- M[, "se_redrawn"]; sf <- M[, "se_fixed"]
  esd_real <- sd(br, na.rm = TRUE); esd_avg <- sd(ba, na.rm = TRUE)
  cov <- function(b, s) mean(abs(b - tb) <= 1.96 * s, na.rm = TRUE)
  dser <- (sr - sf) / esd_real; dser <- dser[is.finite(dser)]      # paired per-rep SER difference
  data.frame(scenario = scenario, n = n, tau = tau, R = sum(is.finite(br)), true_b1 = tb,
    esd_real = esd_real, esd_avg = esd_avg,
    se_fixed = mean(sf, na.rm = TRUE), se_redrawn = mean(sr, na.rm = TRUE),
    ser_fixed = mean(sf, na.rm = TRUE) / esd_real,
    ser_redrawn = mean(sr, na.rm = TRUE) / esd_real,
    ser_avg = mean(sr, na.rm = TRUE) / esd_avg,
    ser_diff = mean(dser), ser_diff_mcse = sd(dser) / sqrt(length(dser)),   # paired DeltaSER +/- MCSE
    cov_real_fixed = cov(br, sf), cov_real_redrawn = cov(br, sr),
    cov_avg_redrawn = cov(ba, sr), row.names = NULL)
}

main <- function(cfg = CFG, out = "results/foldsens.csv") {
  dir.create("results", showWarnings = FALSE); rows <- list(); k <- 1
  for (cell in cfg$cells) {
    sc <- cell[1]; n <- cell[2]; tau <- cell[3]
    cat(sprintf("[%s] foldsens sc=%d n=%d tau=%.2f\n", format(Sys.time(), "%H:%M:%S"), sc, n, tau)); flush.console()
    rows[[k]] <- run_cell(sc, n, tau, cfg); k <- k + 1
    write.csv(do.call(rbind, rows), out, row.names = FALSE)
  }
  cat("done ->", out, "\n")
}
if (sys.nframe() == 0) main()
