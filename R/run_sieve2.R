# ============================================================================
# run_sieve2.R (sim-review #3,#5) — sieve asymptotic check WITH inference.
# Extends run_sieve.R: the CC joint-sieve estimator with DETERMINISTIC J (NOT
# the data-adaptive A-CV) at tau in {0.1,0.5,0.9}. Reports scaled bias with its
# Monte-Carlo SE (= sqrt(n) ESD / sqrt(R)) and a subject-bootstrap coverage, so
# the table verifies inference, not just the n^{-1/2} ESD scaling.
#   EXACT  : m_0 bilinear (in cubic-spline span), fixed df=4.
#   GROWING: m_0 smooth, J_n = ceil(2 n^0.18).
# Writes results/sieve_cov.csv.
# ============================================================================
suppressWarnings(suppressMessages({library(parallel); source("R/run_one.R"); source("R/bootstrap.R")}))

CFG <- list(ns = c(500, 1000, 2000, 4000), R = 400, R_cov = 250, B = 250,
  taus = c(0.10, 0.50, 0.90), scenario = 2, kappa = 0.18, cJ = 2.0, df_exact = 4,
  cores = { sc <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK","")); if (is.na(sc)||sc<1) 4L else sc },
  seed = 20260618)

cc_w <- function(mk, subj) attach_weights(mk, subj, as.numeric(subj$dT == 1))

boot_cc_fixedJ <- function(d, tau, df, B) {        # CC det-J bootstrap SE of beta1
  mk_by <- split(d$mk, d$mk$id); v <- rep(NA, B)
  for (b in seq_len(B)) {
    db <- resample_subjects(d, mk_by); mkb <- prep_mk(db); sjb <- db$subj
    if (is.null(mkb) || nrow(mkb) < 30) next
    v[b] <- tryCatch(fit_A(mkb, cc_w(mkb, sjb), df = df, tau = tau)$beta[1], error = function(e) NA)
  }
  sd(v, na.rm = TRUE)
}

one_rep <- function(r, n, tau, mode, cfg) {
  set.seed(cfg$seed + r * 7919 + (mode == "grow") * 5 + round(tau * 1000))
  p <- default_params(); p$sig1 <- 0.4
  p$surface <- if (mode == "exact") "bilinear" else "smooth"
  df <- if (mode == "exact") cfg$df_exact else max(3, ceiling(cfg$cJ * n^cfg$kappa))
  d <- tryCatch(generate_data(n, scenario = cfg$scenario, params = p), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  mk <- prep_mk(d); subj <- d$subj; if (nrow(mk) < 30) return(NULL)
  b1 <- tryCatch(unname(fit_A(mk, cc_w(mk, subj), df = df, tau = tau)$beta[1]), error = function(e) NA)
  se <- if (r <= cfg$R_cov) tryCatch(boot_cc_fixedJ(d, tau, df, cfg$B), error = function(e) NA) else NA
  c(beta = b1, df = df, se = se)
}

run_cell <- function(n, tau, mode, cfg) {
  M <- do.call(rbind, mclapply(seq_len(cfg$R), one_rep, n = n, tau = tau, mode = mode,
                               cfg = cfg, mc.cores = cfg$cores, mc.preschedule = FALSE))
  tb <- true_beta(tau, { p <- default_params(); p$sig1 <- 0.4; p })["beta1"]
  b <- M[, "beta"]; se <- M[, "se"]
  bias <- mean(b, na.rm = TRUE) - tb; esd <- sd(b, na.rm = TRUE); R <- sum(is.finite(b))
  cov <- mean(abs(b - tb) <= 1.96 * se, na.rm = TRUE)
  data.frame(mode = mode, n = n, tau = tau, J = median(M[, "df"]), R = R,
    bias = bias, sqrtn_bias = sqrt(n) * bias,
    mcse_sqrtn_bias = sqrt(n) * esd / sqrt(R),     # MCSE of the scaled bias (#3)
    esd = esd, sqrtn_esd = sqrt(n) * esd,
    mean_se = mean(se, na.rm = TRUE), cov_boot = cov, R_cov = sum(is.finite(se)))
}

main <- function(cfg = CFG, out = "results/sieve_cov.csv") {
  dir.create("results", showWarnings = FALSE); rows <- list(); k <- 1
  for (mode in c("exact", "grow")) for (tau in cfg$taus) for (n in cfg$ns) {
    cat(sprintf("[%s] sieve2 %s tau=%.2f n=%d\n", format(Sys.time(), "%H:%M:%S"), mode, tau, n)); flush.console()
    rows[[k]] <- run_cell(n, tau, mode, cfg); k <- k + 1
    write.csv(do.call(rbind, rows), out, row.names = FALSE)
  }
  cat("done ->", out, "\n")
}
if (sys.nframe() == 0) main()
