# ============================================================================
# run_sievedist.R (sim-review round-10 #2) — distributional audit of the
# growing-sieve subject bootstrap at the MEDIAN tau=0.5 (in scope), with the
# exact sieve as control, n in {500,1000,2000}. Explains WHY the n=500 growing
# sieve bootstrap is degenerate (Table tab:sieve now marks it n.i.). Uses the
# SAME DGP/seed as run_sieve2.R so this DECOMPOSES that table cell's mean SER.
#
# Per replicate: CC A-joint-sieve beta1_hat (deterministic J), the weighted
# basis-Gram condition number, and the # distinct onset subjects; then a subject
# bootstrap (B resamples) recording beta1*, per-resample Gram condition number,
# and the distinct onset count. Aggregates:
#   * SE distribution: median / IQR / q90 / q95 / q99 / max, and q99/median ratio
#   * fraction of reps with SE > 3x median (concentration in few MC datasets)
#   * 95% interval-length median / q95
#   * bootstrap-estimate extreme spread (|b*-bhat| q99, max)
#   * basis-Gram condition number (original fit + worst resample)
#   * fraction of near-singular ("degenerate") bootstrap resamples
#   * distinct onset-subject count (original + worst resample)
# Writes results/sievedist.csv.
# ============================================================================
suppressWarnings(suppressMessages({library(parallel); source("R/run_one.R"); source("R/bootstrap.R")}))

CFG <- list(ns = c(500, 1000, 2000), R = 250, B = 150, tau = 0.50, scenario = 2,
  kappa = 0.18, cJ = 2.0, df_exact = 4, cond_bad = 1e6,
  cores = { sc <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "")); if (is.na(sc) || sc < 1) 4L else sc },
  seed = 20260618)            # SAME base seed as run_sieve2.R (decomposes that cell)

cc_w <- function(mk, subj) attach_weights(mk, subj, as.numeric(subj$dT == 1))

# CC joint-sieve beta1 + the weighted basis-Gram condition number (over-rich /
# locally-unsupported basis => large kappa => unstable bootstrap fit).
fit_sieve_diag <- function(mk, subj, df, tau) {
  w <- cc_w(mk, subj)
  Z <- make_tensor_basis(mk$U, mk$T, df_u = df, df_t = df)
  D <- cbind(mk$X1, mk$X2, Z)
  G <- crossprod(D * sqrt(w))
  cond <- tryCatch(kappa(G, exact = FALSE), error = function(e) NA_real_)
  fit <- tryCatch(fit_A(mk, w, df = df, tau = tau), error = function(e) NULL)
  list(b1 = if (is.null(fit)) NA_real_ else unname(fit$beta[1]), cond = cond,
       nsubj = length(unique(mk$id)))
}

one_rep <- function(r, n, mode, cfg) {
  set.seed(cfg$seed + r * 7919 + (mode == "grow") * 5 + round(cfg$tau * 1000))
  p <- default_params(); p$sig1 <- 0.4
  p$surface <- if (mode == "exact") "bilinear" else "smooth"
  df <- if (mode == "exact") cfg$df_exact else max(3, ceiling(cfg$cJ * n^cfg$kappa))
  d <- tryCatch(generate_data(n, scenario = cfg$scenario, params = p), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  mk <- prep_mk(d); subj <- d$subj; if (nrow(mk) < 30) return(NULL)
  fd <- fit_sieve_diag(mk, subj, df, cfg$tau)
  if (!is.finite(fd$b1)) return(NULL)
  mk_by <- split(d$mk, d$mk$id)
  bs <- rep(NA_real_, cfg$B); cs <- rep(NA_real_, cfg$B); ns <- rep(NA_real_, cfg$B)
  for (b in seq_len(cfg$B)) {
    db <- resample_subjects(d, mk_by); mkb <- prep_mk(db); sjb <- db$subj
    if (is.null(mkb) || nrow(mkb) < 30) next
    z <- fit_sieve_diag(mkb, sjb, df, cfg$tau)
    bs[b] <- z$b1; cs[b] <- z$cond; ns[b] <- z$nsubj
  }
  se <- sd(bs, na.rm = TRUE)
  list(df = df, b1 = fd$b1, cond0 = fd$cond, nsubj0 = fd$nsubj,
       se = se, ilen = 2 * 1.96 * se,
       bmax = max(abs(bs - fd$b1), na.rm = TRUE),
       bq99 = quantile(abs(bs - fd$b1), 0.99, na.rm = TRUE, names = FALSE),
       frac_bad = mean(cs > cfg$cond_bad, na.rm = TRUE),
       cond_bmax = max(cs, na.rm = TRUE),
       nsubj_min = min(ns, na.rm = TRUE))
}

run_cell <- function(n, mode, cfg) {
  reps <- Filter(Negate(is.null),
                 mclapply(seq_len(cfg$R), one_rep, n = n, mode = mode, cfg = cfg,
                          mc.cores = cfg$cores, mc.preschedule = FALSE))
  g  <- function(f) sapply(reps, function(x) x[[f]])
  qs <- function(v, p) quantile(v, p, na.rm = TRUE, names = FALSE)
  se <- g("se"); ilen <- g("ilen"); med <- median(se, na.rm = TRUE)
  data.frame(mode = mode, n = n, J = median(g("df")), R = length(reps),
    se_med = med, se_iqr = IQR(se, na.rm = TRUE),
    se_q90 = qs(se, .90), se_q95 = qs(se, .95), se_q99 = qs(se, .99), se_max = max(se, na.rm = TRUE),
    se_q99_over_med = qs(se, .99) / med,
    frac_se_gt3med = mean(se > 3 * med, na.rm = TRUE),
    ilen_med = median(ilen, na.rm = TRUE), ilen_q95 = qs(ilen, .95),
    bstar_q99 = qs(g("bq99"), .95), bstar_max = max(g("bmax"), na.rm = TRUE),
    cond0_med = median(g("cond0"), na.rm = TRUE), cond0_q95 = qs(g("cond0"), .95),
    frac_bad_resamp = mean(g("frac_bad"), na.rm = TRUE),
    cond_resamp_max = max(g("cond_bmax"), na.rm = TRUE),
    nsubj0_med = median(g("nsubj0"), na.rm = TRUE),
    nsubj_min_resamp = min(g("nsubj_min"), na.rm = TRUE), row.names = NULL)
}

main <- function(cfg = CFG, out = "results/sievedist.csv") {
  dir.create("results", showWarnings = FALSE); rows <- list(); k <- 1
  for (mode in c("exact", "grow")) for (n in cfg$ns) {
    cat(sprintf("[%s] sievedist %s n=%d\n", format(Sys.time(), "%H:%M:%S"), mode, n)); flush.console()
    rows[[k]] <- run_cell(n, mode, cfg); k <- k + 1
    write.csv(do.call(rbind, rows), out, row.names = FALSE)
  }
  cat("done ->", out, "\n")
}
if (sys.nframe() == 0) main()
