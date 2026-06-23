# ============================================================================
# run_tuning.R (sim-review #14) — sensitivity of the default CC estimators to the
# tuning that theory leaves to rates: the sparsity bandwidth a_n and the
# projection dimension df_g. The consolidated runs fix a_n=0.05, df_g=4; this
# checks the B-OS-CF (and A-CV) beta_{tau,1} bias/ESD are stable across a grid,
# so the headline results are not an artifact of one tuning choice. Sc 2.
# Writes results/tuning.csv.
# ============================================================================
suppressWarnings(suppressMessages({library(parallel); source("R/cf.R")}))

CFG <- list(R = 400, K = 5, df_grid = c(3, 4, 5),
  a_ns = c(0.03, 0.05, 0.08), df_gs = c(3, 4, 5),
  # expanded to n=500 (least stable) and the tails tau=0.9/0.1 (review #10/#11)
  cells = { g <- expand.grid(n = c(500, 1000, 2000), tau = c(0.10, 0.50, 0.90))
            lapply(seq_len(nrow(g)), function(i) c(n = g$n[i], tau = g$tau[i])) },
  cores = { sc <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK","")); if (is.na(sc)||sc<1) 4L else sc },
  seed = 20260623)

cc_w <- function(mk, subj) attach_weights(mk, subj, as.numeric(subj$dT == 1))

# plug-in clustered sandwich SE (X1 component), bread condition number, and f-cap
# rate at the SELECTED A-CV fit but with the tuning-specific f (bandwidth a_n) and
# g (dim df_g), so SER/coverage sensitivity to (a_n,df_g) is exposed, not just bias.
sandwich_se <- function(mk, w, fit, df, tau, a_n, df_g) {
  X <- cbind(mk$X1, mk$X2)
  qhat <- predict_A(fit, mk$U, mk$T, mk$X1, mk$X2); psi <- tau - ((mk$L - qhat) <= 0)
  flo <- fit_A(mk, w, df = df, tau = tau - a_n); fhi <- fit_A(mk, w, df = df, tau = tau + a_n)
  f <- sparsity_fhat(flo, fhi, mk$U, mk$T, mk$X1, mk$X2, a_n)
  g <- predict_projection(fit_projection(mk$U, mk$T, X, w * f, df_g = df_g), mk$U, mk$T)
  Xc <- X - g; nsub <- length(unique(mk$id))
  S <- matrix(0, 2, 2); for (k in seq_len(nrow(mk))) S <- S + w[k] * f[k] * tcrossprod(Xc[k, ]); S <- S / nsub
  Z <- rowsum(w * Xc * psi, group = mk$id); V <- crossprod(as.matrix(Z)) / nsub
  Si <- tryCatch(solve(S), error = function(e) NULL); if (is.null(Si)) return(c(se = NA, cond = NA, fcap = NA))
  se <- sqrt(diag(Si %*% V %*% t(Si) / nsub))[1]
  c(se = unname(se), cond = unname(kappa(S)), fcap = mean(f >= 50 | f <= 1e-2))
}

one_rep <- function(r, n, tau, a_n, df_g, cfg) {
  set.seed(cfg$seed + r * 7919 + round(tau * 1000) + round(a_n * 1000) + df_g)
  p <- default_params(); p$sig1 <- 0.4
  d <- tryCatch(generate_data(n, scenario = 2, params = p), error = function(e) NULL)
  na <- c(acv = NA, cf = NA, se = NA, cond = NA, fcap = NA)
  if (is.null(d)) return(na)
  mk <- prep_mk(d); subj <- d$subj; if (nrow(mk) < 30) return(na)
  w <- cc_w(mk, subj)
  df <- tryCatch(select_df_cv(mk, subj, w, tau, df_grid = cfg$df_grid), error = function(e) 4)
  fitA <- tryCatch(fit_A(mk, w, df = df, tau = tau), error = function(e) NULL)
  if (is.null(fitA)) return(na)
  acv <- as.numeric(fitA$beta[1])
  cf <- one_step_cf(d, tau, df_fixed = df, type = "cc", K = cfg$K, a_n = a_n, df_g = df_g)
  sw <- tryCatch(sandwich_se(mk, w, fitA, df, tau, a_n, df_g), error = function(e) c(se = NA, cond = NA, fcap = NA))
  c(acv = acv, cf = if (is.null(cf)) NA else cf$beta[1], se = unname(sw["se"]),
    cond = unname(sw["cond"]), fcap = unname(sw["fcap"]))
}

run_cell <- function(n, tau, a_n, df_g, cfg) {
  n <- unname(n); tau <- unname(tau)        # cell[] passes named scalars; unname
  M <- do.call(rbind, mclapply(seq_len(cfg$R), one_rep, n = n, tau = tau, a_n = a_n, df_g = df_g,
                               cfg = cfg, mc.cores = cfg$cores, mc.preschedule = FALSE))
  tb <- as.numeric(true_beta(tau, { p <- default_params(); p$sig1 <- 0.4; p })["beta1"])
  acv <- M[, "acv"]; esd_acv <- sd(acv, na.rm = TRUE); se <- M[, "se"]
  cov_acv <- mean(abs(acv - tb) <= 1.96 * se, na.rm = TRUE)         # sandwich coverage
  data.frame(n = n, tau = tau, a_n = a_n, df_g = df_g, R = sum(is.finite(M[, "cf"])),
    bias_acv = mean(acv, na.rm = TRUE) - tb, esd_acv = esd_acv,
    bias_cf = mean(M[, "cf"], na.rm = TRUE) - tb, esd_cf = sd(M[, "cf"], na.rm = TRUE),
    ser_acv = mean(se, na.rm = TRUE) / esd_acv, cov_acv = cov_acv,
    fcap = mean(M[, "fcap"], na.rm = TRUE), cond_med = median(M[, "cond"], na.rm = TRUE),
    row.names = NULL)
}

main <- function(cfg = CFG, out = "results/tuning.csv") {
  dir.create("results", showWarnings = FALSE); rows <- list(); k <- 1
  for (cell in cfg$cells) for (a_n in cfg$a_ns) for (df_g in cfg$df_gs) {
    cat(sprintf("[%s] tuning n=%d tau=%.2f a_n=%.2f df_g=%d\n", format(Sys.time(), "%H:%M:%S"),
                cell[["n"]], cell[["tau"]], a_n, df_g)); flush.console()
    rows[[k]] <- run_cell(cell[["n"]], cell[["tau"]], a_n, df_g, cfg); k <- k + 1
    write.csv(do.call(rbind, rows), out, row.names = FALSE)
  }
  cat("done ->", out, "\n")
}
if (sys.nframe() == 0) main()
