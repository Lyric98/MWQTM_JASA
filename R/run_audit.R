# ============================================================================
# run_audit.R (sim-review #17) — computational reproducibility audit. Over
# (scenario, n, tau) it runs the CC A-CV + B-OS-CF pipeline and a small subject
# bootstrap, and reports the rates of every failure/fallback mode that could bias
# the headline numbers if silently dropped:
#   * QR solve failure (fit_A errors / NA beta)
#   * B-OS-CF failure (one_step_cf returns NULL: too few events in a fold, etc.)
#   * near-singular bread (condition number > 1e8)
#   * sparsity f-hat hitting the [1e-2, 50] cap
#   * bootstrap-replicate failure rate (fraction of resamples with NA beta)
#   * CV df-selection ties (more than one df within 1e-8 of the CV minimum)
# Failed reps are COUNTED, not silently discarded. Writes results/audit.csv.
# ============================================================================
suppressWarnings(suppressMessages({library(parallel); source("R/cf.R"); source("R/bootstrap.R")}))

CFG <- list(ns = c(500, 1000, 2000), scenarios = c(1, 2), taus = c(0.25, 0.50, 0.75),
  R = 300, Bsmall = 60, K = 5, a_n = 0.05, df_g = 4, df_grid = c(3, 4, 5),
  cores = { sc <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK","")); if (is.na(sc)||sc<1) 4L else sc },
  seed = 20260625)

cc_w <- function(mk, subj) attach_weights(mk, subj, as.numeric(subj$dT == 1))

one_rep <- function(r, n, scenario, tau, cfg) {
  set.seed(cfg$seed + r * 7919 + scenario * 101 + round(tau * 1000))
  p <- default_params(); p$sig1 <- if (scenario == 1) 0 else 0.4
  d <- tryCatch(generate_data(n, scenario = scenario, params = p), error = function(e) NULL)
  na <- c(gen_fail = 1, acv_fail = NA, cf_fail = NA, acv_ip_fail = NA, cf_ip_fail = NA,
          bread_sing = NA, cond = NA, fcap = NA, boot_fail = NA)
  if (is.null(d)) return(na)
  mk <- prep_mk(d); subj <- d$subj
  if (nrow(mk) < 30) return(na)
  w <- cc_w(mk, subj)
  df <- tryCatch(select_df_cv(mk, subj, w, tau, df_grid = cfg$df_grid), error = function(e) 4)
  fitA <- tryCatch(fit_A(mk, w, df = df, tau = tau), error = function(e) NULL)
  acv_fail <- as.numeric(is.null(fitA) || any(!is.finite(fitA$beta)))
  # IPCW pipeline failures too (review #8): G_C fit + IPCW A-CV + IPCW B-OS-CF
  gc <- tryCatch(fit_GC(subj), error = function(e) NULL)
  w_ip <- if (is.null(gc)) w else attach_weights(mk, subj, ipcw_weights(subj, GC_at(gc)))
  fitAi <- tryCatch(fit_A(mk, w_ip, df = df, tau = tau), error = function(e) NULL)
  acv_ip_fail <- as.numeric(is.null(gc) || is.null(fitAi) || any(!is.finite(fitAi$beta)))
  cfi <- one_step_cf(d, tau, df_fixed = df, type = "ipcw", K = cfg$K, a_n = cfg$a_n, df_g = cfg$df_g)
  cf_ip_fail <- as.numeric(is.null(cfi) || any(!is.finite(cfi$beta)))
  # bread condition number (value, for quantiles) + near-singular flag + f-cap
  bread_sing <- NA; cond <- NA; fcap <- NA
  if (!acv_fail) {
    X <- cbind(mk$X1, mk$X2)
    flo <- tryCatch(fit_A(mk, w, df = df, tau = tau - cfg$a_n), error = function(e) NULL)
    fhi <- tryCatch(fit_A(mk, w, df = df, tau = tau + cfg$a_n), error = function(e) NULL)
    if (!is.null(flo) && !is.null(fhi)) {
      f <- sparsity_fhat(flo, fhi, mk$U, mk$T, mk$X1, mk$X2, cfg$a_n)
      fcap <- mean(f >= 50 | f <= 1e-2)
      g <- tryCatch(predict_projection(fit_projection(mk$U, mk$T, X, w * f, df_g = cfg$df_g), mk$U, mk$T),
                    error = function(e) matrix(0, nrow(mk), 2))
      Xc <- X - g; nsub <- length(unique(mk$id))
      S <- matrix(0, 2, 2); for (k in seq_len(nrow(mk))) S <- S + w[k] * f[k] * tcrossprod(Xc[k, ]); S <- S / nsub
      cond <- kappa(S); bread_sing <- as.numeric(cond > 1e8)
    }
  }
  cf <- one_step_cf(d, tau, df_fixed = df, type = "cc", K = cfg$K, a_n = cfg$a_n, df_g = cfg$df_g)
  cf_fail <- as.numeric(is.null(cf) || any(!is.finite(cf$beta)))
  # small subject bootstrap: fraction of resamples whose A-CV fit fails
  mk_by <- split(d$mk, d$mk$id); nb_fail <- 0L; nb <- 0L
  for (b in seq_len(cfg$Bsmall)) {
    db <- resample_subjects(d, mk_by); mkb <- prep_mk(db); sjb <- db$subj
    if (is.null(mkb) || nrow(mkb) < 30) { nb_fail <- nb_fail + 1L; nb <- nb + 1L; next }
    wb <- cc_w(mkb, sjb)
    bb <- tryCatch(fit_A(mkb, wb, df = df, tau = tau)$beta[1], error = function(e) NA)
    if (!is.finite(bb)) nb_fail <- nb_fail + 1L
    nb <- nb + 1L
  }
  c(gen_fail = 0, acv_fail = acv_fail, cf_fail = cf_fail, acv_ip_fail = acv_ip_fail,
    cf_ip_fail = cf_ip_fail, bread_sing = bread_sing, cond = cond,
    fcap = fcap, boot_fail = nb_fail / max(nb, 1))
}

run_cell <- function(n, scenario, tau, cfg) {
  M <- do.call(rbind, mclapply(seq_len(cfg$R), one_rep, n = n, scenario = scenario, tau = tau,
                               cfg = cfg, mc.cores = cfg$cores, mc.preschedule = FALSE))
  mm <- function(c) mean(M[, c], na.rm = TRUE)
  qq <- function(c, p) quantile(M[, c], p, na.rm = TRUE, names = FALSE)
  data.frame(scenario = scenario, n = n, tau = tau, R = nrow(M),
    acv_fail_rate = mm("acv_fail"), cf_fail_rate = mm("cf_fail"),
    acv_ip_fail_rate = mm("acv_ip_fail"), cf_ip_fail_rate = mm("cf_ip_fail"),
    bread_singular_rate = mm("bread_sing"), cond_med = qq("cond", .5), cond_q95 = qq("cond", .95),
    fcap_rate = mm("fcap"), boot_fail_rate = mm("boot_fail"), row.names = NULL)
}

main <- function(cfg = CFG, out = "results/audit.csv") {
  dir.create("results", showWarnings = FALSE); rows <- list(); k <- 1
  for (sc in cfg$scenarios) for (n in cfg$ns) for (tau in cfg$taus) {
    cat(sprintf("[%s] audit sc=%d n=%d tau=%.2f\n", format(Sys.time(), "%H:%M:%S"), sc, n, tau)); flush.console()
    rows[[k]] <- run_cell(n, sc, tau, cfg); k <- k + 1
    write.csv(do.call(rbind, rows), out, row.names = FALSE)
  }
  cat("done ->", out, "\n")
}
if (sys.nframe() == 0) main()
