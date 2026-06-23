# ============================================================================
# run_consolidated.R  (sim-review priorities 1,2,10) — the INTERNAL-CONSISTENCY
# batch. On the SAME Monte-Carlo replicates, fit all four estimators
#   {CC, IPCW} x {A-CV, B-OS-CF}
# so Table "pilot" (IPCW) and the CC-vs-IPCW comparison are one experiment.
# Stores PER-REPLICATE coverage indicators for a PAIRED test (McNemar /
# equivalence). Full-refit subject bootstrap (df fixed) for all four, incl. the
# cross-fitted one-step. SLURM array: one task per (scenario,n,tau) cell.
#
# Full submission scale: R=1000 point, R_cov=500 coverage, B=500 bootstrap.
# ============================================================================
suppressWarnings(suppressMessages({library(parallel); source("R/cf.R"); source("R/bootstrap.R")}))

CFG <- list(
  R = 1000, R_cov = 500, B = 500, K = 5, a_n = 0.05,
  df_grid = c(3, 4, 5), df_g = 4,
  cores = { sc <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "")); if (is.na(sc) || sc < 1) 4L else sc },
  seed = 20260616)

CELLS <- expand.grid(scenario = c(1, 2), n = c(500, 1000, 2000),
                     tau = c(0.25, 0.50, 0.75))

# ---- point estimators on one dataset (CV-selected df) ----------------------
point_all <- function(d, tau, cfg) {
  mk <- prep_mk(d); subj <- d$subj
  if (nrow(mk) < 30) return(NULL)
  # Each pipeline tunes its df by CV under ITS OWN weighting, so the CC member is
  # censoring-nuisance-free end-to-end (no hat-G_C in the CC df selection).
  w_ip <- attach_weights(mk, subj, ipcw_weights(subj, GC_at(fit_GC(subj))))
  w_cc <- attach_weights(mk, subj, as.numeric(subj$dT == 1))
  df_ip <- tryCatch(select_df_cv(mk, subj, w_ip, tau, df_grid = cfg$df_grid), error = function(e) 4)
  df_cc <- tryCatch(select_df_cv(mk, subj, w_cc, tau, df_grid = cfg$df_grid), error = function(e) 4)
  acv_ip <- tryCatch(as.numeric(fit_A(mk, w_ip, df = df_ip, tau = tau)$beta), error = function(e) c(NA, NA))
  acv_cc <- tryCatch(as.numeric(fit_A(mk, w_cc, df = df_cc, tau = tau)$beta), error = function(e) c(NA, NA))
  cf_ip <- one_step_cf(d, tau, df_fixed = df_ip, type = "ipcw", K = cfg$K, a_n = cfg$a_n,
                       df_g = cfg$df_g, beta_acv = acv_ip)
  cf_cc <- one_step_cf(d, tau, df_fixed = df_cc, type = "cc",   K = cfg$K, a_n = cfg$a_n,
                       df_g = cfg$df_g, beta_acv = acv_cc)
  list(df_ip = df_ip, df_cc = df_cc,
       acv_cc = acv_cc[1], acv_ip = acv_ip[1],
       cf_cc = if (is.null(cf_cc)) NA else cf_cc$beta[1],
       cf_ip = if (is.null(cf_ip)) NA else cf_ip$beta[1],
       corr_cc = if (is.null(cf_cc)) NA else cf_cc$correction[1],
       corr_ip = if (is.null(cf_ip)) NA else cf_ip$correction[1])
}

# ---- full-refit bootstrap of all four (df fixed) ---------------------------
# returns bootstrap SD c(acv_cc, acv_ip, cf_cc, cf_ip)
boot_all <- function(d, tau, df_cc, df_ip, cfg) {
  mk_by <- split(d$mk, d$mk$id); B <- cfg$B; bm <- matrix(NA, B, 4)
  for (b in seq_len(B)) {
    db <- resample_subjects(d, mk_by); mkb <- prep_mk(db); sjb <- db$subj
    if (is.null(mkb) || nrow(mkb) < 30) next
    w_ip <- attach_weights(mkb, sjb, ipcw_weights(sjb, GC_at(fit_GC(sjb))))
    w_cc <- attach_weights(mkb, sjb, as.numeric(sjb$dT == 1))
    bm[b, 1] <- tryCatch(fit_A(mkb, w_cc, df = df_cc, tau = tau)$beta[1], error = function(e) NA)
    bm[b, 2] <- tryCatch(fit_A(mkb, w_ip, df = df_ip, tau = tau)$beta[1], error = function(e) NA)
    cfc <- one_step_cf(db, tau, df_fixed = df_cc, type = "cc",   K = cfg$K, a_n = cfg$a_n, df_g = cfg$df_g)
    cfi <- one_step_cf(db, tau, df_fixed = df_ip, type = "ipcw", K = cfg$K, a_n = cfg$a_n, df_g = cfg$df_g)
    bm[b, 3] <- if (is.null(cfc)) NA else cfc$beta[1]
    bm[b, 4] <- if (is.null(cfi)) NA else cfi$beta[1]
  }
  apply(bm, 2, sd, na.rm = TRUE)
}

one_rep <- function(r, n, scenario, tau, cfg) {
  set.seed(cfg$seed + r * 7919 + scenario * 101 + round(tau * 1000))
  p <- default_params(); p$sig1 <- if (scenario == 1) 0 else 0.4
  d <- tryCatch(generate_data(n, scenario = scenario, params = p), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  pt <- point_all(d, tau, cfg)
  if (is.null(pt)) return(NULL)
  out <- pt
  if (r <= cfg$R_cov) {
    se <- tryCatch(boot_all(d, tau, pt$df_cc, pt$df_ip, cfg), error = function(e) rep(NA, 4))
    out$se_acv_cc <- se[1]; out$se_acv_ip <- se[2]
    out$se_cf_cc  <- se[3]; out$se_cf_ip  <- se[4]
  }
  out
}

run_cell <- function(n, scenario, tau, cfg) {
  tb <- true_beta(tau, { p <- default_params(); p$sig1 <- if (scenario == 1) 0 else 0.4; p })["beta1"]
  tb <- as.numeric(tb)
  reps <- mclapply(seq_len(cfg$R), one_rep, n = n, scenario = scenario, tau = tau,
                   cfg = cfg, mc.cores = cfg$cores, mc.preschedule = FALSE)
  reps <- Filter(Negate(is.null), reps)
  g <- function(f) sapply(reps, function(x) { v <- x[[f]]; if (is.null(v)) NA else v })
  est <- list(acv_cc = g("acv_cc"), acv_ip = g("acv_ip"),
              cf_cc = g("cf_cc"), cf_ip = g("cf_ip"))
  se  <- list(acv_cc = g("se_acv_cc"), acv_ip = g("se_acv_ip"),
              cf_cc = g("se_cf_cc"), cf_ip = g("se_cf_ip"))
  covered <- function(b, s) abs(b - tb) <= 1.96 * s          # per-rep indicator
  # paired (same-rep) coverage indicators on the bootstrap subset
  I <- lapply(names(est), function(nm) covered(est[[nm]], se[[nm]])); names(I) <- names(est)
  bias <- function(b) mean(b, na.rm = TRUE) - tb
  esd  <- function(b) sd(b, na.rm = TRUE)
  cov  <- function(ind) mean(ind, na.rm = TRUE)
  # McNemar discordance for CC vs IPCW (same rep): b = CC covers & IPCW not, etc.
  mcnemar <- function(icc, iip) {
    ok <- is.finite(icc) & is.finite(iip)
    c(b = sum(icc[ok] & !iip[ok]), c = sum(!icc[ok] & iip[ok]), nd = sum(ok))
  }
  mn_acv <- mcnemar(I$acv_cc, I$acv_ip); mn_cf <- mcnemar(I$cf_cc, I$cf_ip)
  corr_cc <- g("corr_cc"); corr_ip <- g("corr_ip")
  std_cc <- corr_cc / se$cf_cc; std_ip <- corr_ip / se$cf_ip
  data.frame(
    scenario = scenario, n = n, tau = tau, true_b1 = tb, R = length(reps),
    R_cov = sum(is.finite(se$acv_ip)),
    bias_acv_cc = bias(est$acv_cc), bias_acv_ip = bias(est$acv_ip),
    bias_cf_cc = bias(est$cf_cc),  bias_cf_ip = bias(est$cf_ip),
    esd_acv_cc = esd(est$acv_cc), esd_acv_ip = esd(est$acv_ip),
    esd_cf_cc = esd(est$cf_cc),   esd_cf_ip = esd(est$cf_ip),
    se_acv_cc = mean(se$acv_cc, na.rm = TRUE), se_acv_ip = mean(se$acv_ip, na.rm = TRUE),
    se_cf_cc = mean(se$cf_cc, na.rm = TRUE),   se_cf_ip = mean(se$cf_ip, na.rm = TRUE),
    cov_acv_cc = cov(I$acv_cc), cov_acv_ip = cov(I$acv_ip),
    cov_cf_cc = cov(I$cf_cc),   cov_cf_ip = cov(I$cf_ip),
    esd_ratio_acv = esd(est$acv_cc) / esd(est$acv_ip),
    esd_ratio_cf  = esd(est$cf_cc)  / esd(est$cf_ip),
    cov_diff_acv = cov(I$acv_cc) - cov(I$acv_ip),
    cov_diff_cf  = cov(I$cf_cc)  - cov(I$cf_ip),
    mcnemar_b_acv = mn_acv["b"], mcnemar_c_acv = mn_acv["c"], nd_acv = mn_acv["nd"],
    mcnemar_b_cf = mn_cf["b"],   mcnemar_c_cf = mn_cf["c"],   nd_cf = mn_cf["nd"],
    corr_cf_cc = mean(corr_cc, na.rm = TRUE), corr_cf_ip = mean(corr_ip, na.rm = TRUE),
    sd_corr_cc = sd(corr_cc, na.rm = TRUE),   sd_corr_ip = sd(corr_ip, na.rm = TRUE),
    std_corr_cc = mean(std_cc, na.rm = TRUE), std_corr_ip = mean(std_ip, na.rm = TRUE),
    row.names = NULL)
}

main <- function(cfg = CFG) {
  dir.create("results", showWarnings = FALSE)
  task <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", ""))
  if (!is.na(task)) {                       # one cell per array task
    cell <- CELLS[task, ]
    cat(sprintf("[%s] consolidated cell %d: sc=%d n=%d tau=%.2f\n",
                format(Sys.time(), "%H:%M:%S"), task, cell$scenario, cell$n, cell$tau)); flush.console()
    row <- run_cell(cell$n, cell$scenario, cell$tau, cfg)
    write.csv(row, sprintf("results/consolidated_cell_%02d.csv", task), row.names = FALSE)
    cat("done cell", task, "\n")
  } else {                                  # local: loop all cells
    rows <- list()
    for (i in seq_len(nrow(CELLS))) {
      cell <- CELLS[i, ]
      cat(sprintf("[%s] cell %d sc=%d n=%d tau=%.2f\n", format(Sys.time(), "%H:%M:%S"),
                  i, cell$scenario, cell$n, cell$tau)); flush.console()
      rows[[i]] <- run_cell(cell$n, cell$scenario, cell$tau, cfg)
      write.csv(do.call(rbind, rows), "results/consolidated.csv", row.names = FALSE)
    }
    cat("done -> results/consolidated.csv\n")
  }
}
if (sys.nframe() == 0) main()
