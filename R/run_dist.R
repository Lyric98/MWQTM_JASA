# ============================================================================
# run_dist.R (sim-review #14) — distributional stress test. The marker error is
# a standardized non-normal law (heavy-tail t5 or right-skew lognormal) via a
# Gaussian copula (within-subject correlation preserved). QR should be agnostic
# to the error shape; we verify bias / ESD / bootstrap coverage of the default
# CC--A-CV beta_{tau,1} AND beta_{tau,2} at the tails, Scenario 2.
# Writes results/dist.csv.
# ============================================================================
suppressWarnings(suppressMessages({library(parallel); source("R/run_one.R"); source("R/bootstrap.R")}))

CFG <- list(R = 500, R_cov = 300, B = 300, df_grid = c(3, 4, 5),
  dists = c("t5", "skew"), ns = c(500, 1000, 2000), taus = c(0.25, 0.50, 0.75),
  cores = { sc <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK","")); if (is.na(sc)||sc<1) 4L else sc },
  seed = 20260619)

cc_w <- function(mk, subj) attach_weights(mk, subj, as.numeric(subj$dT == 1))

boot_cc_b12 <- function(d, tau, df, B) {       # CC A-CV bootstrap SE of (beta1,beta2)
  mk_by <- split(d$mk, d$mk$id); bm <- matrix(NA, B, 2)
  for (b in seq_len(B)) {
    db <- resample_subjects(d, mk_by); mkb <- prep_mk(db); sjb <- db$subj
    if (is.null(mkb) || nrow(mkb) < 30) next
    bm[b, ] <- tryCatch(as.numeric(fit_A(mkb, cc_w(mkb, sjb), df = df, tau = tau)$beta),
                        error = function(e) c(NA, NA))
  }
  apply(bm, 2, sd, na.rm = TRUE)
}

one_rep <- function(r, n, tau, dist, cfg) {
  set.seed(cfg$seed + r * 7919 + round(tau * 1000) + nchar(dist))
  p <- default_params(); p$sig1 <- 0.4; p$err_dist <- dist
  d <- tryCatch(generate_data(n, scenario = 2, params = p), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  mk <- prep_mk(d); subj <- d$subj; if (nrow(mk) < 30) return(NULL)
  df <- tryCatch(select_df_cv(mk, subj, cc_w(mk, subj), tau, df_grid = cfg$df_grid), error = function(e) 4)
  b <- tryCatch(as.numeric(fit_A(mk, cc_w(mk, subj), df = df, tau = tau)$beta), error = function(e) c(NA, NA))
  se <- if (r <= cfg$R_cov) tryCatch(boot_cc_b12(d, tau, df, cfg$B), error = function(e) c(NA, NA)) else c(NA, NA)
  c(b1 = b[1], b2 = b[2], se1 = se[1], se2 = se[2])
}

run_cell <- function(n, tau, dist, cfg) {
  p <- default_params(); p$sig1 <- 0.4; p$err_dist <- dist; tb <- true_beta(tau, p)
  M <- do.call(rbind, mclapply(seq_len(cfg$R), one_rep, n = n, tau = tau, dist = dist,
                               cfg = cfg, mc.cores = cfg$cores, mc.preschedule = FALSE))
  cov <- function(b, s, t0) mean(abs(b - t0) <= 1.96 * s, na.rm = TRUE)
  data.frame(dist = dist, n = n, tau = tau, R = sum(is.finite(M[, "b1"])),
    true_b1 = tb["beta1"], true_b2 = tb["beta2"],
    bias_b1 = mean(M[, "b1"], na.rm = TRUE) - tb["beta1"], esd_b1 = sd(M[, "b1"], na.rm = TRUE),
    bias_b2 = mean(M[, "b2"], na.rm = TRUE) - tb["beta2"], esd_b2 = sd(M[, "b2"], na.rm = TRUE),
    cov_b1 = cov(M[, "b1"], M[, "se1"], tb["beta1"]),
    cov_b2 = cov(M[, "b2"], M[, "se2"], tb["beta2"]),
    R_cov = sum(is.finite(M[, "se1"])), row.names = NULL)
}

main <- function(cfg = CFG, out = "results/dist.csv") {
  dir.create("results", showWarnings = FALSE); rows <- list(); k <- 1
  for (dist in cfg$dists) for (n in cfg$ns) for (tau in cfg$taus) {
    cat(sprintf("[%s] dist=%s n=%d tau=%.2f\n", format(Sys.time(), "%H:%M:%S"), dist, n, tau)); flush.console()
    rows[[k]] <- run_cell(n, tau, dist, cfg); k <- k + 1
    write.csv(do.call(rbind, rows), out, row.names = FALSE)
  }
  cat("done ->", out, "\n")
}
if (sys.nframe() == 0) main()
