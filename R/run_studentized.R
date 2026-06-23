# ============================================================================
# run_studentized.R (sim-review round-10 #3) — distributional / studentization
# diagnostic for the headline CC--A-CV estimator at the central/quartile
# quantiles tau in {0.25,0.5,0.75}, BOTH scenarios, n in {500,1000,2000}.
# Targets the lowest-coverage headline cell (Sc 2, n=500, tau=0.25: cov 0.918,
# SER~1.00, bias~0.006), where the shortfall is neither mean-SE nor bias.
# Per replicate: CC A-CV beta1_hat, a subject bootstrap (df fixed at the rep's
# CV df, exactly as the headline boot) giving SE and the bootstrap 2.5/97.5
# quantiles; then per-rep Z = (bhat - b0)/SE. Reports, per cell:
#   * Z moments: mean, sd, skewness, excess kurtosis, 2.5%/97.5% quantiles
#   * lower-side / upper-side normal-CI miss rates (asymmetry)
#   * Corr(|bhat-b0|, SE) across reps (SE-error dependence)
#   * coverage of NORMAL-SE vs PERCENTILE vs BASIC bootstrap CIs
# Uses the SAME DGP/seed scheme as run_consolidated so the cells line up with
# Table tab:estse / tab:estsim. Writes results/studentized.csv.
# ============================================================================
suppressWarnings(suppressMessages({library(parallel); source("R/run_one.R"); source("R/bootstrap.R")}))

CFG <- list(R = 400, B = 250, df_grid = c(3, 4, 5),
  scenarios = c(1, 2), ns = c(500, 1000, 2000), taus = c(0.25, 0.50, 0.75),
  cores = { sc <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "")); if (is.na(sc) || sc < 1) 4L else sc },
  seed = 20260616)            # SAME base seed as run_consolidated.R

cc_w <- function(mk, subj) attach_weights(mk, subj, as.numeric(subj$dT == 1))

# CC A-CV bootstrap (df fixed): vector of beta1* over B resamples (= headline CC path)
boot_cc_dist <- function(d, tau, df, B) {
  mk_by <- split(d$mk, d$mk$id); v <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    db <- resample_subjects(d, mk_by); mkb <- prep_mk(db); sjb <- db$subj
    if (is.null(mkb) || nrow(mkb) < 30) next
    v[b] <- tryCatch(fit_A(mkb, cc_w(mkb, sjb), df = df, tau = tau)$beta[1], error = function(e) NA)
  }
  v
}

one_rep <- function(r, n, scenario, tau, cfg) {
  set.seed(cfg$seed + r * 7919 + scenario * 101 + round(tau * 1000))
  p <- default_params(); p$sig1 <- if (scenario == 1) 0 else 0.4
  d <- tryCatch(generate_data(n, scenario = scenario, params = p), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  mk <- prep_mk(d); subj <- d$subj; if (nrow(mk) < 30) return(NULL)
  w_ip <- attach_weights(mk, subj, ipcw_weights(subj, GC_at(fit_GC(subj))))
  df <- tryCatch(select_df_cv(mk, subj, w_ip, tau, df_grid = cfg$df_grid), error = function(e) 4)
  bhat <- tryCatch(fit_A(mk, cc_w(mk, subj), df = df, tau = tau)$beta[1], error = function(e) NA)
  if (!is.finite(bhat)) return(NULL)
  bstar <- boot_cc_dist(d, tau, df, cfg$B)
  q <- quantile(bstar, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  c(bhat = unname(bhat), se = sd(bstar, na.rm = TRUE), q025 = q[1], q975 = q[2],
    nb = sum(is.finite(bstar)))
}

skew   <- function(z) { z <- z[is.finite(z)]; m <- mean(z); mean((z - m)^3) / mean((z - m)^2)^1.5 }
exkurt <- function(z) { z <- z[is.finite(z)]; m <- mean(z); mean((z - m)^4) / mean((z - m)^2)^2 - 3 }

run_cell <- function(n, scenario, tau, cfg) {
  tb <- as.numeric(true_beta(tau, { p <- default_params(); p$sig1 <- if (scenario == 1) 0 else 0.4; p })["beta1"])
  M <- do.call(rbind, mclapply(seq_len(cfg$R), one_rep, n = n, scenario = scenario, tau = tau,
                               cfg = cfg, mc.cores = cfg$cores, mc.preschedule = FALSE))
  bhat <- M[, "bhat"]; se <- M[, "se"]; q025 <- M[, "q025"]; q975 <- M[, "q975"]
  ok <- is.finite(bhat) & is.finite(se) & se > 0
  bhat <- bhat[ok]; se <- se[ok]; q025 <- q025[ok]; q975 <- q975[ok]
  Z   <- (bhat - tb) / se
  err <- abs(bhat - tb)
  data.frame(scenario = scenario, n = n, tau = tau, R = length(Z), true_b1 = tb,
    z_mean = mean(Z), z_sd = sd(Z), z_skew = skew(Z), z_exkurt = exkurt(Z),
    z_q025 = quantile(Z, 0.025, names = FALSE), z_q975 = quantile(Z, 0.975, names = FALSE),
    miss_lo = mean(Z >  1.96),     # b0 below the normal CI (Z = (bhat-b0)/SE large +)
    miss_hi = mean(Z < -1.96),     # b0 above the normal CI
    corr_err_se = cor(err, se),    # large errors paired with small SE => masks undercov
    cov_normal = mean(err <= 1.96 * se),
    cov_pct    = mean(q025 <= tb & tb <= q975),
    cov_basic  = mean((2 * bhat - q975) <= tb & tb <= (2 * bhat - q025)),
    row.names = NULL)
}

main <- function(cfg = CFG, out = "results/studentized.csv") {
  dir.create("results", showWarnings = FALSE); rows <- list(); k <- 1
  for (sc in cfg$scenarios) for (n in cfg$ns) for (tau in cfg$taus) {
    cat(sprintf("[%s] studentized sc=%d n=%d tau=%.2f\n", format(Sys.time(), "%H:%M:%S"), sc, n, tau)); flush.console()
    rows[[k]] <- run_cell(n, sc, tau, cfg); k <- k + 1
    write.csv(do.call(rbind, rows), out, row.names = FALSE)
  }
  cat("done ->", out, "\n")
}
if (sys.nframe() == 0) main()
