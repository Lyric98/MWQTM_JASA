# ============================================================================
# run_corrdist.R (sim-review #8) — distribution of the ABSOLUTE B-OS-CF one-step
# correction |beta_CF - beta_ACV|, to replace the signed-mean reporting (a signed
# mean near 0 can mask per-replicate cancellation). Default CC member, Scenario 2.
# No bootstrap (point estimates only) => cheap. Writes results/corrdist.csv.
# ============================================================================
suppressWarnings(suppressMessages({library(parallel); source("R/cf.R")}))

CFG <- list(R = 500, K = 5, a_n = 0.05, df_grid = c(3, 4, 5),
  ns = c(500, 1000, 2000), taus = c(0.25, 0.50, 0.75),
  cores = { sc <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK","")); if (is.na(sc)||sc<1) 4L else sc },
  seed = 20260622)

cc_w <- function(mk, subj) attach_weights(mk, subj, as.numeric(subj$dT == 1))

one_rep <- function(r, n, tau, cfg) {
  set.seed(cfg$seed + r * 7919 + round(tau * 1000))
  p <- default_params(); p$sig1 <- 0.4
  na4 <- c(corr = NA, c_os = NA, c_split = NA, acv = NA)
  d <- tryCatch(generate_data(n, scenario = 2, params = p), error = function(e) NULL)
  if (is.null(d)) return(na4)
  mk <- prep_mk(d); subj <- d$subj; if (nrow(mk) < 30) return(na4)
  w <- cc_w(mk, subj)
  df <- tryCatch(select_df_cv(mk, subj, w, tau, df_grid = cfg$df_grid), error = function(e) 4)
  acv <- tryCatch(unname(fit_A(mk, w, df = df, tau = tau)$beta[1]), error = function(e) NA)
  cf <- one_step_cf(d, tau, df_fixed = df, type = "cc", K = cfg$K, a_n = cfg$a_n, beta_acv = c(acv, NA))
  c(corr = if (is.null(cf)) NA else cf$correction[1],
    c_os = if (is.null(cf)) NA else cf$c_os[1],
    c_split = if (is.null(cf)) NA else cf$c_split[1], acv = acv)
}

run_cell <- function(n, tau, cfg) {
  M <- do.call(rbind, mclapply(seq_len(cfg$R), one_rep, n = n, tau = tau, cfg = cfg,
                               mc.cores = cfg$cores, mc.preschedule = FALSE))
  cc <- M[, "corr"]; ac <- abs(cc[is.finite(cc)]); esd <- sd(M[, "acv"], na.rm = TRUE)
  aos <- abs(M[, "c_os"][is.finite(M[, "c_os"])]); asp <- abs(M[, "c_split"][is.finite(M[, "c_split"])])
  q <- quantile(ac, c(.5, .9, .95), names = FALSE)
  data.frame(n = n, tau = tau, R = length(ac),
    mean_signed = mean(cc, na.rm = TRUE), mean_abs = mean(ac),
    median_abs = q[1], q90_abs = q[2], q95_abs = q[3],
    esd_acv = esd, mean_abs_over_esd = mean(ac) / esd,
    mean_abs_os = mean(aos), mean_abs_split = mean(asp),     # decomposition (#2)
    os_over_esd = mean(aos) / esd, split_over_esd = mean(asp) / esd, row.names = NULL)
}

main <- function(cfg = CFG, out = "results/corrdist.csv") {
  dir.create("results", showWarnings = FALSE); rows <- list(); k <- 1
  for (n in cfg$ns) for (tau in cfg$taus) {
    cat(sprintf("[%s] corrdist n=%d tau=%.2f\n", format(Sys.time(), "%H:%M:%S"), n, tau)); flush.console()
    rows[[k]] <- run_cell(n, tau, cfg); k <- k + 1
    write.csv(do.call(rbind, rows), out, row.names = FALSE)
  }
  cat("done ->", out, "\n")
}
if (sys.nframe() == 0) main()
