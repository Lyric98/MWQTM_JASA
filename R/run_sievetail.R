# ============================================================================
# run_sievetail.R (sim-review #2,#9) — does the proposed REMEDY fix the
# growing-sieve tail failure? In the failure regime (deterministic growing
# sieve, extreme tau), compare three estimators on the SAME replicates:
#   A-js   : A-joint-sieve at deterministic J     (shows the tail bias)
#   B-OS-CF: cross-fitted one-step on the det-J base (the proposed debiasing)
#   A-rich : undersmoothed joint sieve, J+2        (alternative remedy)
# Report sqrt(n)*bias (+ MCSE), ESD, and bootstrap coverage for each, so the
# table shows whether one-step / undersmoothing removes the tail sieve bias.
# CC default; Scenario 2; smooth surface (growing sieve). Writes sievetail.csv.
# ============================================================================
suppressWarnings(suppressMessages({library(parallel); source("R/cf.R"); source("R/bootstrap.R")}))

CFG <- list(ns = c(1000, 2000, 4000), R = 400, R_cov = 150, B = 120, K = 5,
  taus = c(0.90, 0.10, 0.50), kappa = 0.18, cJ = 2.0, a_n = 0.05,
  cores = { sc <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK","")); if (is.na(sc)||sc<1) 4L else sc },
  seed = 20260621)

cc_w <- function(mk, subj) attach_weights(mk, subj, as.numeric(subj$dT == 1))
Jof <- function(n, cfg) max(3, ceiling(cfg$cJ * n^cfg$kappa))

# CC det-J A-CV bootstrap SE; and CC B-OS-CF full CF bootstrap SE
boot_three <- function(d, tau, J, cfg) {
  mk_by <- split(d$mk, d$mk$id); bm <- matrix(NA, cfg$B, 3)
  for (b in seq_len(cfg$B)) {
    db <- resample_subjects(d, mk_by); mkb <- prep_mk(db); sjb <- db$subj
    if (is.null(mkb) || nrow(mkb) < 30) next
    w <- cc_w(mkb, sjb)
    bm[b, 1] <- tryCatch(unname(fit_A(mkb, w, df = J, tau = tau)$beta[1]), error = function(e) NA)
    cf <- one_step_cf(db, tau, df_fixed = J, type = "cc", K = cfg$K, a_n = cfg$a_n)
    bm[b, 2] <- if (is.null(cf)) NA else cf$beta[1]
    bm[b, 3] <- tryCatch(unname(fit_A(mkb, w, df = J + 2, tau = tau)$beta[1]), error = function(e) NA)
  }
  apply(bm, 2, sd, na.rm = TRUE)
}

one_rep <- function(r, n, tau, cfg) {
  set.seed(cfg$seed + r * 7919 + round(tau * 1000))
  p <- default_params(); p$sig1 <- 0.4; p$surface <- "smooth"
  J <- Jof(n, cfg)
  d <- tryCatch(generate_data(n, scenario = 2, params = p), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  mk <- prep_mk(d); subj <- d$subj; if (nrow(mk) < 30) return(NULL)
  w <- cc_w(mk, subj)
  a_js   <- tryCatch(unname(fit_A(mk, w, df = J, tau = tau)$beta[1]), error = function(e) NA)
  a_rich <- tryCatch(unname(fit_A(mk, w, df = J + 2, tau = tau)$beta[1]), error = function(e) NA)
  cf <- one_step_cf(d, tau, df_fixed = J, type = "cc", K = cfg$K, a_n = cfg$a_n)
  b_cf <- if (is.null(cf)) NA else cf$beta[1]
  se <- if (r <= cfg$R_cov) tryCatch(boot_three(d, tau, J, cfg), error = function(e) rep(NA, 3)) else rep(NA, 3)
  c(J = J, a_js = a_js, b_cf = b_cf, a_rich = a_rich, se_js = se[1], se_cf = se[2], se_rich = se[3])
}

run_cell <- function(n, tau, cfg) {
  M <- do.call(rbind, mclapply(seq_len(cfg$R), one_rep, n = n, tau = tau, cfg = cfg,
                               mc.cores = cfg$cores, mc.preschedule = FALSE))
  tb <- as.numeric(true_beta(tau, { p <- default_params(); p$sig1 <- 0.4; p })["beta1"])
  R <- function(col) sum(is.finite(M[, col]))
  stat <- function(bcol, scol) {
    b <- M[, bcol]; bias <- mean(b, na.rm = TRUE) - tb; esd <- sd(b, na.rm = TRUE); Rb <- sum(is.finite(b))
    mse <- mean(M[, scol], na.rm = TRUE)
    cov <- mean(abs(b - tb) <= 1.96 * M[, scol], na.rm = TRUE)
    c(bias = bias, sqrtn_bias = sqrt(n) * bias, mcse = sqrt(n) * esd / sqrt(Rb),
      sqrtn_esd = sqrt(n) * esd, esd = esd, se = mse, ser = mse / esd,
      rmse = sqrt(bias^2 + esd^2), cov = cov)
  }
  s_js <- stat("a_js", "se_js"); s_cf <- stat("b_cf", "se_cf"); s_ri <- stat("a_rich", "se_rich")
  data.frame(n = n, tau = tau, J = median(M[, "J"]), R = R("a_js"),
    js_sqrtn_bias = s_js["sqrtn_bias"], js_mcse = s_js["mcse"], js_sqrtn_esd = s_js["sqrtn_esd"],
    js_esd = s_js["esd"], js_se = s_js["se"], js_ser = s_js["ser"], js_rmse = s_js["rmse"], js_cov = s_js["cov"],
    cf_sqrtn_bias = s_cf["sqrtn_bias"], cf_mcse = s_cf["mcse"], cf_sqrtn_esd = s_cf["sqrtn_esd"],
    cf_esd = s_cf["esd"], cf_se = s_cf["se"], cf_ser = s_cf["ser"], cf_rmse = s_cf["rmse"], cf_cov = s_cf["cov"],
    rich_sqrtn_bias = s_ri["sqrtn_bias"], rich_mcse = s_ri["mcse"],
    rich_esd = s_ri["esd"], rich_ser = s_ri["ser"], rich_rmse = s_ri["rmse"], rich_cov = s_ri["cov"],
    row.names = NULL)
}

main <- function(cfg = CFG, out = "results/sievetail.csv") {
  dir.create("results", showWarnings = FALSE); rows <- list(); k <- 1
  for (tau in cfg$taus) for (n in cfg$ns) {
    cat(sprintf("[%s] sievetail tau=%.2f n=%d\n", format(Sys.time(), "%H:%M:%S"), tau, n)); flush.console()
    rows[[k]] <- run_cell(n, tau, cfg); k <- k + 1
    write.csv(do.call(rbind, rows), out, row.names = FALSE)
  }
  cat("done ->", out, "\n")
}
if (sys.nframe() == 0) main()
