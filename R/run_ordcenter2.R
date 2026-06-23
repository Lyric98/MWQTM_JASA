# ============================================================================
# run_ordcenter2.R (sim-review round-11 #2/#3) — DW vs ORD centering, redone so
# the global significance test is LEGITIMATE and the bread-calibration claim is
# quantified.
#   * All three tau in {0.25,0.5,0.75} are fit on the SAME dataset per replicate,
#     and replicates (within each n, seeded independently across n) are the
#     INDEPENDENT clusters. The per-replicate net statistic
#         D_r = sum_tau ( I{DW covers} - I{ORD covers} )
#     is then the unit of a replicate-CLUSTERED sign test / cluster bootstrap,
#     which respects the cross-tau dependence the pooled binomial ignored.
#   * Interval geometry: mean interval length DW vs ORD, mean SE ratio, and the
#     fraction of (rep,tau) where the DW interval NESTS the ORD interval
#     (se_dw >= se_ord, since the point estimates coincide). This pins the effect
#     to "density-weighted bread => larger, better-calibrated SE", not a different
#     estimator.
# Scenario 2 only (Sc 1 is the DW==ORD sanity, omitted). Oracle f isolates the
# centering. Writes results/ordcenter2.csv (per cell) + results/ordcenter_global.csv.
# ============================================================================
suppressMessages(library(parallel)); source("R/run_one.R")

CFG <- list(ns = c(500, 1000, 2000), tau = c(0.25, 0.50, 0.75), scenario = 2,
  R = 600, a_n = 0.05, df_grid = c(3, 4, 5), nperm = 200000,
  cores = { sc <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "")); if (is.na(sc) || sc < 1) 4L else sc },
  seed = 20260619)

# proper 95% interval score (Gneiting-Raftery): width + miss penalties, lower better
intscore <- function(b, se, y, alpha = 0.05) {
  L <- b - 1.96 * se; U <- b + 1.96 * se
  (U - L) + (2 / alpha) * ((L - y) * (y < L) + (y - U) * (y > U))
}

# one replicate: fit DW & ORD at ALL taus on the SAME dataset; return, per tau,
# the coverage indicators and SEs (oracle f, IPCW weight, NCF one-step + bread).
one_rep <- function(r, n, cfg, tbv) {
  set.seed(cfg$seed + n * 7L + r * 7919L)
  p <- default_params(); p$sig1 <- 0.4
  d <- tryCatch(generate_data(n, scenario = cfg$scenario, params = p), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  mk <- prep_mk(d); subj <- d$subj; if (nrow(mk) < 30) return(NULL)
  gc <- tryCatch(fit_GC(subj), error = function(e) NULL); if (is.null(gc)) return(NULL)
  w <- attach_weights(mk, subj, ipcw_weights(subj, GC_at(gc)))
  X <- cbind(mk$X1, mk$X2); nsub <- length(unique(mk$id))
  out <- matrix(NA_real_, length(cfg$tau), 4,
               dimnames = list(NULL, c("b_dw", "se_dw", "b_ord", "se_ord")))
  for (j in seq_along(cfg$tau)) {
    tau <- cfg$tau[j]
    df <- tryCatch(select_df_cv(mk, subj, w, tau, df_grid = cfg$df_grid), error = function(e) 4)
    fitcv <- tryCatch(fit_A(mk, w, df = df, tau = tau), error = function(e) NULL)
    if (is.null(fitcv)) next
    qhat <- predict_A(fitcv, mk$U, mk$T, mk$X1, mk$X2); psi <- tau - ((mk$L - qhat) <= 0)
    fhat <- dnorm(qnorm(tau)) / sigma_X(mk$X1, p$sig0, p$sig1)   # oracle density
    onestep <- function(ghat) {
      Xc <- X - ghat; Uvec <- colSums(w * Xc * psi) / nsub
      S <- matrix(0, 2, 2); for (k in seq_len(nrow(mk))) S <- S + w[k] * fhat[k] * tcrossprod(Xc[k, ]); S <- S / nsub
      beta <- (fitcv$beta + solve(S, Uvec))[1]
      se <- sandwich_se(mk, list(w_row = w, Xc = Xc, psi = psi, S = S))[1]
      c(unname(beta), se)
    }
    g_dw  <- predict_projection(fit_projection(mk$U, mk$T, X, w * fhat, df_g = 4), mk$U, mk$T)
    g_ord <- predict_projection(fit_projection(mk$U, mk$T, X, w,        df_g = 4), mk$U, mk$T)
    dw <- onestep(g_dw); ord <- onestep(g_ord)
    out[j, ] <- c(dw[1], dw[2], ord[1], ord[2])
  }
  # per-tau coverage indicators (+ SE geometry) on this one dataset
  Idw <- as.integer(abs(out[, "b_dw"]  - tbv) <= 1.96 * out[, "se_dw"])
  Ior <- as.integer(abs(out[, "b_ord"] - tbv) <= 1.96 * out[, "se_ord"])
  dis <- intscore(out[, "b_dw"], out[, "se_dw"], tbv) -      # per-tau paired interval-score diff (DW-ORD)
         intscore(out[, "b_ord"], out[, "se_ord"], tbv)
  list(out = out, Idw = Idw, Ior = Ior, Dr = sum(Idw - Ior, na.rm = TRUE),
       dis = dis, disDr = mean(dis, na.rm = TRUE),
       se_dw = out[, "se_dw"], se_ord = out[, "se_ord"])
}

run_all <- function(cfg) {
  tbv <- sapply(cfg$tau, function(tau) as.numeric(true_beta(tau, { p <- default_params(); p$sig1 <- 0.4; p })["beta1"]))
  names(tbv) <- as.character(cfg$tau)
  percell <- list(); kc <- 1; Dr_all <- c(); disDr_all <- c()   # one entry per (n,replicate) cluster
  for (n in cfg$ns) {
    reps <- Filter(Negate(is.null),
                   mclapply(seq_len(cfg$R), one_rep, n = n, cfg = cfg, tbv = tbv,
                            mc.cores = cfg$cores, mc.preschedule = FALSE))
    Dr_all <- c(Dr_all, sapply(reps, function(z) z$Dr))
    disDr_all <- c(disDr_all, sapply(reps, function(z) z$disDr))
    for (j in seq_along(cfg$tau)) {
      tau <- cfg$tau[j]; tb <- tbv[j]
      bd  <- sapply(reps, function(z) z$out[j, "b_dw"]);  sd_ <- sapply(reps, function(z) z$out[j, "se_dw"])
      bo  <- sapply(reps, function(z) z$out[j, "b_ord"]); so  <- sapply(reps, function(z) z$out[j, "se_ord"])
      Id  <- sapply(reps, function(z) z$Idw[j]); Io <- sapply(reps, function(z) z$Ior[j])
      ok <- is.finite(bd) & is.finite(bo) & is.finite(sd_) & is.finite(so)
      dch <- Id[ok] - Io[ok]
      disv <- intscore(bd[ok], sd_[ok], tb) - intscore(bo[ok], so[ok], tb)   # paired DeltaIS per rep
      bdisc <- sum(Id[ok] == 1 & Io[ok] == 0); cdisc <- sum(Id[ok] == 0 & Io[ok] == 1)
      mcp <- function(b, c) { m <- b + c; if (m == 0) NA else min(1, 2 * pbinom(min(b, c), m, 0.5)) }
      percell[[kc]] <- data.frame(scenario = cfg$scenario, n = n, tau = tau, R = sum(ok),
        true_b1 = as.numeric(tb),
        bias_dw = mean(bd[ok]) - tb, esd_dw = sd(bd[ok]), cov_dw = mean(Id[ok]),
        bias_ord = mean(bo[ok]) - tb, esd_ord = sd(bo[ok]), cov_ord = mean(Io[ok]),
        cov_diff = mean(dch), cov_diff_mcse = sd(dch) / sqrt(length(dch)),
        mcnemar_b = bdisc, mcnemar_c = cdisc, mcnemar_p = mcp(bdisc, cdisc),
        ilen_dw = mean(2 * 1.96 * sd_[ok]), ilen_ord = mean(2 * 1.96 * so[ok]),
        len_ratio = mean((sd_[ok]) / so[ok]),    # = SE ratio = length ratio
        se_ratio = mean(sd_[ok] / so[ok]), nest_frac = mean(sd_[ok] >= so[ok]),
        # proper 95% interval score (lower = better; penalizes width AND misses)
        is_dw = mean(intscore(bd[ok], sd_[ok], tb)), is_ord = mean(intscore(bo[ok], so[ok], tb)),
        dis_mean = mean(disv), dis_mcse = sd(disv) / sqrt(length(disv)),   # paired DeltaIS +/- MCSE
        row.names = NULL); kc <- kc + 1
    }
  }
  # replicate-CLUSTERED global test on D_r (each (n,replicate) an independent cluster)
  Dr_all <- Dr_all[is.finite(Dr_all)]
  Knz <- sum(Dr_all != 0); Kpos <- sum(Dr_all > 0); Kneg <- sum(Dr_all < 0)
  # exact two-sided cluster sign test on net-direction (signs of nonzero D_r)
  signp <- if (Knz == 0) NA else min(1, 2 * pbinom(min(Kpos, Kneg), Knz, 0.5))
  # sign-flip permutation p for mean(D_r) (cluster-level), Monte-Carlo
  obs <- mean(Dr_all); set.seed(1)
  perm <- replicate(cfg$nperm, abs(mean(Dr_all * sample(c(-1, 1), length(Dr_all), replace = TRUE))))
  perm_p <- (1 + sum(perm >= abs(obs))) / (cfg$nperm + 1)
  # cluster bootstrap CI for mean(D_r)
  set.seed(2); bs <- replicate(20000, mean(sample(Dr_all, replace = TRUE)))
  # replicate-clustered paired interval-score difference (DeltaIS averaged over tau per cluster)
  disDr_all <- disDr_all[is.finite(disDr_all)]
  set.seed(3); bsd <- replicate(20000, mean(sample(disDr_all, replace = TRUE)))
  glob <- data.frame(n_clusters = length(Dr_all), K_nonzero = Knz, K_pos = Kpos, K_neg = Kneg,
    mean_Dr = obs, ci_lo = quantile(bs, .025, names = FALSE), ci_hi = quantile(bs, .975, names = FALSE),
    sign_test_p = signp, signflip_perm_p = perm_p,
    dis_clustered_mean = mean(disDr_all),
    dis_ci_lo = quantile(bsd, .025, names = FALSE), dis_ci_hi = quantile(bsd, .975, names = FALSE),
    row.names = NULL)
  list(percell = do.call(rbind, percell), glob = glob)
}

main <- function(cfg = CFG) {
  dir.create("results", showWarnings = FALSE)
  cat(sprintf("[%s] ordcenter2 starting\n", format(Sys.time(), "%H:%M:%S"))); flush.console()
  res <- run_all(cfg)
  write.csv(res$percell, "results/ordcenter2.csv", row.names = FALSE)
  write.csv(res$glob,    "results/ordcenter_global.csv", row.names = FALSE)
  cat("done -> results/ordcenter2.csv, results/ordcenter_global.csv\n")
  print(res$glob)
}
if (sys.nframe() == 0) main()
