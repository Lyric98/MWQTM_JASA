# ============================================================================
# run_pilot.R — first-round pilot. Answers the ONE question:
#   does the feasible beta_tau have small bias and correct coverage,
#   and how does it degrade at extreme tau?
#
# Per replicate: point estimates (A-CV, A-US, B-OS, O3) + clustered-sandwich SE.
# Subset of replicates: full-refit & frozen-G_C bootstrap (IPCW-variance check).
# Aggregates to bias / SB / ESD / SER / coverage; saves CSV + figure data.
# ============================================================================
suppressMessages(library(parallel))
source("R/run_one.R")
source("R/bootstrap.R")

CFG <- list(
  ns       = c(500, 1000, 2000),
  scenarios= c(1, 2),
  taus     = c(0.10, 0.25, 0.50, 0.75, 0.90),
  R        = 200,                 # replicates
  a_n      = 0.05,                # quantile-spacing bandwidth
  df_grid  = c(3, 4, 5),          # CV grid (bias flat in df per gate0b2)
  df_us_bump = 2,
  boot_subset = 150,              # replicates that also get bootstrap (v3: larger)
  B        = 200,
  # respect the SLURM allocation; fall back to a safe small default
  cores    = {
    sc <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", ""))
    if (is.na(sc) || sc < 1) 4L else sc
  },
  seed     = 20260613
)

one_rep <- function(r, n, scenario, tau, oracle, cfg) {
  set.seed(cfg$seed + r * 7919 + scenario * 101 + round(tau * 1000))
  d <- generate_data(n, scenario = scenario)
  est <- tryCatch(estimate_all(d, tau, a_n = cfg$a_n,
                               df_us_bump = cfg$df_us_bump,
                               df_grid = cfg$df_grid, oracle = oracle),
                  error = function(e) NULL)
  if (is.null(est)) return(NULL)
  do_boot <- (r <= cfg$boot_subset)
  bse <- NULL
  if (do_boot)
    bse <- tryCatch(boot_se(d, tau, B = cfg$B, a_n = cfg$a_n,
                            df_fixed = est$acv$df, do_frozen = TRUE),
                    error = function(e) NULL)
  list(
    b_acv = est$acv$beta[1], df = est$acv$df,
    b_aus = est$aus$beta[1],
    b_bos = est$bos$beta[1], se_bos = est$bos$se[1], Scond = est$bos$Scond,
    b_o3  = est$o3$beta[1],  se_o3  = est$o3$se[1],
    imse  = est$imse,
    boot_full = if (is.null(bse)) NULL else bse$full,
    boot_froz = if (is.null(bse)) NULL else bse$frozen
  )
}

run_cell <- function(n, scenario, tau, cfg, oracle) {
  reps <- mclapply(seq_len(cfg$R), one_rep, n = n, scenario = scenario,
                   tau = tau, oracle = oracle, cfg = cfg,
                   mc.cores = cfg$cores)
  reps <- Filter(Negate(is.null), reps)
  tb <- true_beta(tau, default_params2(scenario))["beta1"]
  pick <- function(f) sapply(reps, f)
  b_acv <- pick(function(x) x$b_acv); b_aus <- pick(function(x) x$b_aus)
  b_bos <- pick(function(x) x$b_bos); b_o3 <- pick(function(x) x$b_o3)
  se_bos <- pick(function(x) x$se_bos); se_o3 <- pick(function(x) x$se_o3)
  Scond <- pick(function(x) x$Scond); dfsel <- pick(function(x) x$df)
  imse  <- pick(function(x) x$imse)
  # bootstrap subset: full-refit and frozen-G_C SEs (each B x 3)
  bfull <- lapply(reps, function(x) x$boot_full)
  bfroz <- lapply(reps, function(x) x$boot_froz)
  has   <- !sapply(bfull, is.null)
  Bf <- if (any(has)) do.call(rbind, bfull[has]) else matrix(NA, 0, 3)
  Bz <- if (any(has)) do.call(rbind, bfroz[has]) else matrix(NA, 0, 3)

  cover <- function(bhat, se) mean(abs(bhat - tb) <= 1.96 * se, na.rm = TRUE)
  esd <- function(b) sd(b, na.rm = TRUE)
  row <- data.frame(
    n = n, scenario = scenario, tau = tau, true_b1 = as.numeric(tb),
    bias_acv = mean(b_acv) - tb, bias_aus = mean(b_aus) - tb,
    bias_bos = mean(b_bos) - tb, bias_o3 = mean(b_o3) - tb,
    esd_acv = esd(b_acv), esd_aus = esd(b_aus),
    esd_bos = esd(b_bos), esd_o3 = esd(b_o3),
    cov_bos_sw = cover(b_bos, se_bos), cov_o3_sw = cover(b_o3, se_o3),
    ser_bos = mean(se_bos, na.rm = TRUE) / esd(b_bos),
    ser_o3  = mean(se_o3,  na.rm = TRUE) / esd(b_o3),
    Scond_med = median(Scond, na.rm = TRUE), df_med = median(dfsel),
    imse_med = median(imse, na.rm = TRUE),
    cov_acv_boot = NA, cov_aus_boot = NA, cov_bos_boot = NA,
    cov_acv_froz = NA, cov_bos_froz = NA
  )
  if (nrow(Bf) > 0) {
    idx <- which(has)
    row$cov_acv_boot <- mean(abs(b_acv[idx] - tb) <= 1.96 * Bf[, 1], na.rm = TRUE)
    row$cov_aus_boot <- mean(abs(b_aus[idx] - tb) <= 1.96 * Bf[, 2], na.rm = TRUE)
    row$cov_bos_boot <- mean(abs(b_bos[idx] - tb) <= 1.96 * Bf[, 3], na.rm = TRUE)
    row$cov_acv_froz <- mean(abs(b_acv[idx] - tb) <= 1.96 * Bz[, 1], na.rm = TRUE)
    row$cov_bos_froz <- mean(abs(b_bos[idx] - tb) <= 1.96 * Bz[, 3], na.rm = TRUE)
  }
  row
}

# true beta for a scenario (sig1=0 vs 0.4)
default_params2 <- function(scenario) {
  p <- default_params(); p$sig1 <- if (scenario == 1) 0 else 0.4; p
}

main <- function(cfg = CFG, out = "results/pilot_v3.csv") {
  dir.create("results", showWarnings = FALSE)
  rows <- list(); k <- 1
  for (sc in cfg$scenarios) {
    oracle <- build_oracle(sc, n_pool = 8000, df_g = 4)   # one per scenario
    for (n in cfg$ns) for (tau in cfg$taus) {
      cat(sprintf("[%s] scenario=%d n=%d tau=%.2f\n",
                  format(Sys.time(), "%H:%M:%S"), sc, n, tau)); flush.console()
      rows[[k]] <- run_cell(n, sc, tau, cfg, oracle); k <- k + 1
      write.csv(do.call(rbind, rows), out, row.names = FALSE)
    }
  }
  cat("done ->", out, "\n")
}

if (sys.nframe() == 0) main()
