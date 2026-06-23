# ============================================================================
# run_sieveresample.R (sim-review round-12 #3) — isolate the numerical trigger of
# the n=500 growing-sieve SE blow-up at the BOOTSTRAP-RESAMPLE level (Table 22's
# dataset-level metrics could not). For each dataset (grow n=500 primary; exact
# n=500 and grow n=1000 as controls; tau=0.5) and each subject bootstrap resample
# we record the resampled estimate beta1* and resample-level degeneracy:
#   * n_zerosupp : # basis columns with < 3 resample rows carrying mass (lost support)
#   * lam_min*   : min eigenvalue of the standardized resample design Gram
#   * n_onset*   : # DISTINCT original onset subjects represented in the resample
# Then, pooling all resamples, we compare the EXTREME-|beta*-betahat| tail (top
# 0.5%) against the bulk, and correlate |beta*-betahat| with the degeneracy, to
# test whether the blow-up resamples are exactly the support-losing ones.
# Writes results/sieveresample.csv.
# ============================================================================
suppressWarnings(suppressMessages({library(parallel); source("R/run_one.R"); source("R/bootstrap.R")}))

CFG <- list(cells = list(c("grow", 500), c("exact", 500), c("grow", 1000)),
  R = 250, B = 150, tau = 0.50, scenario = 2, kappa = 0.18, cJ = 2.0, df_exact = 4,
  supp_tol = 0.01, supp_min = 3, extreme_q = 0.995,
  cores = { sc <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "")); if (is.na(sc) || sc < 1) 4L else sc },
  seed = 20260618)

cc_w <- function(mk, subj) attach_weights(mk, subj, as.numeric(subj$dT == 1))

# beta1 + resample-level degeneracy on one (resampled) dataset
resamp_diag <- function(mk, subj, df, tau, cfg) {
  w <- cc_w(mk, subj)
  b1 <- tryCatch(unname(fit_A(mk, w, df = df, tau = tau)$beta[1]), error = function(e) NA_real_)
  Z <- make_tensor_basis(mk$U, mk$T, df_u = df, df_t = df)
  n_zerosupp <- sum(colSums(Z > cfg$supp_tol) < cfg$supp_min)
  D <- cbind(mk$X1, mk$X2, Z); sw <- sqrt(w); Dw <- D * sw
  cn <- sqrt(colSums(Dw^2)); cn[cn < 1e-300] <- 1e-300
  Ds <- sweep(Dw, 2, cn, "/")
  lam_min <- tryCatch(min(svd(crossprod(Ds), nu = 0, nv = 0)$d), error = function(e) NA_real_)
  oid <- if ("orig_id" %in% names(subj)) subj$orig_id else subj$id
  n_onset <- length(unique(oid[subj$dT == 1]))
  c(b1 = b1, n_zerosupp = n_zerosupp, lam_min = lam_min, n_onset = n_onset)
}

one_rep <- function(r, mode, n, cfg) {
  set.seed(cfg$seed + r * 7919 + (mode == "grow") * 5 + round(cfg$tau * 1000))
  p <- default_params(); p$sig1 <- 0.4
  p$surface <- if (mode == "exact") "bilinear" else "smooth"
  df <- if (mode == "exact") cfg$df_exact else max(3, ceiling(cfg$cJ * n^cfg$kappa))
  d <- tryCatch(generate_data(n, scenario = cfg$scenario, params = p), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  mk <- prep_mk(d); subj <- d$subj; if (nrow(mk) < 30) return(NULL)
  bhat <- tryCatch(unname(fit_A(mk, cc_w(mk, subj), df = df, tau = cfg$tau)$beta[1]), error = function(e) NA)
  if (!is.finite(bhat)) return(NULL)
  mk_by <- split(d$mk, d$mk$id)
  out <- matrix(NA_real_, cfg$B, 4, dimnames = list(NULL, c("absdev", "n_zerosupp", "lam_min", "n_onset")))
  for (b in seq_len(cfg$B)) {
    db <- resample_subjects(d, mk_by); mkb <- prep_mk(db); sjb <- db$subj
    if (is.null(mkb) || nrow(mkb) < 30) next
    z <- resamp_diag(mkb, sjb, df, cfg$tau, cfg)
    out[b, ] <- c(abs(z["b1"] - bhat), z["n_zerosupp"], z["lam_min"], z["n_onset"])
  }
  out[is.finite(out[, "absdev"]), , drop = FALSE]
}

run_cell <- function(mode, n, cfg) {
  M <- do.call(rbind, Filter(Negate(is.null),
        mclapply(seq_len(cfg$R), one_rep, mode = mode, n = n, cfg = cfg,
                 mc.cores = cfg$cores, mc.preschedule = FALSE)))
  ad <- M[, "absdev"]; zs <- M[, "n_zerosupp"]; lm <- M[, "lam_min"]; no <- M[, "n_onset"]
  thr <- quantile(ad, cfg$extreme_q, names = FALSE, na.rm = TRUE)
  ex <- ad >= thr; bk <- ad < thr
  data.frame(mode = mode, n = n, n_resamples = nrow(M),
    absdev_med = median(ad), absdev_max = max(ad), extreme_thr = thr, n_extreme = sum(ex),
    # extreme tail vs bulk degeneracy
    zerosupp_extreme = mean(zs[ex]), zerosupp_bulk = mean(zs[bk]),
    lammin_extreme = median(lm[ex], na.rm = TRUE), lammin_bulk = median(lm[bk], na.rm = TRUE),
    onset_extreme = mean(no[ex]), onset_bulk = mean(no[bk]),
    # correlations of blow-up magnitude with resample degeneracy
    corr_absdev_zerosupp = suppressWarnings(cor(ad, zs, use = "complete.obs")),
    corr_absdev_neglogLam = suppressWarnings(cor(log(ad + 1e-12), -log(pmax(lm, 1e-300)), use = "complete.obs")),
    corr_absdev_onset = suppressWarnings(cor(ad, no, use = "complete.obs")),
    row.names = NULL)
}

main <- function(cfg = CFG, out = "results/sieveresample.csv") {
  dir.create("results", showWarnings = FALSE); rows <- list(); k <- 1
  for (cell in cfg$cells) {
    mode <- cell[1]; n <- as.integer(cell[2])
    cat(sprintf("[%s] sieveresample %s n=%d\n", format(Sys.time(), "%H:%M:%S"), mode, n)); flush.console()
    rows[[k]] <- run_cell(mode, n, cfg); k <- k + 1
    write.csv(do.call(rbind, rows), out, row.names = FALSE)
  }
  cat("done ->", out, "\n")
}
if (sys.nframe() == 0) main()
