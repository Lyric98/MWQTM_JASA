# ============================================================================
# run_sievedist2.R (sim-review round-11, their Table 21) — establish (not just
# assert) WHY the n=500 growing-sieve bootstrap SE blows up. The round-10 audit
# used the RAW basis-Gram kappa, which is huge even where SE is stable (the basis
# columns have wildly different scales), so it cannot prove causation. Here we use
# scale-invariant degeneracy diagnostics and test the SE<->degeneracy link
# per Monte-Carlo dataset.
#
# Per dataset (CC A-joint-sieve, deterministic J, tau=0.5; same DGP/seed family as
# run_sieve2.R): bootstrap SE of beta1, and on the ORIGINAL design [X1,X2,basis]
#   * standardized Gram (columns scaled to unit norm) min eigenvalue lambda_min
#   * effective rank = (sum sv)^2 / sum(sv^2) and numerical rank (sv > 1e-8 max)
#   * min nonzero singular value
#   * # basis columns with NO local support (<5 rows with appreciable mass)
# Then per cell we test the causal link across datasets:
#   * corr( log SE , -log lambda_min )  and corr( log SE , #unsupported )
#   * 2x2 coincidence: P(large SE | rank-deficient), P(rank-deficient | large SE),
#     and whether the extreme-SE and most-degenerate datasets are the SAME ones.
# Writes results/sievedist2.csv.
# ============================================================================
suppressWarnings(suppressMessages({library(parallel); source("R/run_one.R"); source("R/bootstrap.R")}))

CFG <- list(ns = c(500, 1000, 2000), R = 250, B = 150, tau = 0.50, scenario = 2,
  kappa = 0.18, cJ = 2.0, df_exact = 4, supp_tol = 0.01, supp_min = 5,
  cores = { sc <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "")); if (is.na(sc) || sc < 1) 4L else sc },
  seed = 20260618)

cc_w <- function(mk, subj) attach_weights(mk, subj, as.numeric(subj$dT == 1))

# scale-invariant degeneracy diagnostics of the weighted design [X1,X2,basis]
gram_diag <- function(mk, subj, df, w, cfg) {
  Z <- make_tensor_basis(mk$U, mk$T, df_u = df, df_t = df)
  D <- cbind(mk$X1, mk$X2, Z)
  sw <- sqrt(w); Dw <- D * sw
  cn <- sqrt(colSums(Dw^2)); cn[cn < 1e-300] <- 1e-300
  Ds <- sweep(Dw, 2, cn, "/")                 # unit-norm columns => standardized Gram
  G <- crossprod(Ds)
  sv <- svd(G, nu = 0, nv = 0)$d
  lam_min <- min(sv); mx <- max(sv)
  eff_rank <- (sum(sv))^2 / sum(sv^2)
  num_rank <- sum(sv > 1e-8 * mx)
  min_sv <- min(sv[sv > 1e-12 * mx])
  # # basis columns with no local support: < supp_min rows carrying mass > supp_tol
  Bsupp <- colSums(Z > cfg$supp_tol)
  n_unsupp <- sum(Bsupp < cfg$supp_min)
  c(lam_min = lam_min, eff_rank = eff_rank, num_rank = num_rank, ncol = ncol(D),
    min_sv = min_sv, n_unsupp = n_unsupp)
}

one_rep <- function(r, n, mode, cfg) {
  set.seed(cfg$seed + r * 7919 + (mode == "grow") * 5 + round(cfg$tau * 1000))
  p <- default_params(); p$sig1 <- 0.4
  p$surface <- if (mode == "exact") "bilinear" else "smooth"
  df <- if (mode == "exact") cfg$df_exact else max(3, ceiling(cfg$cJ * n^cfg$kappa))
  d <- tryCatch(generate_data(n, scenario = cfg$scenario, params = p), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  mk <- prep_mk(d); subj <- d$subj; if (nrow(mk) < 30) return(NULL)
  w <- cc_w(mk, subj)
  b1 <- tryCatch(unname(fit_A(mk, w, df = df, tau = cfg$tau)$beta[1]), error = function(e) NA)
  if (!is.finite(b1)) return(NULL)
  gd <- tryCatch(gram_diag(mk, subj, df, w, cfg), error = function(e) rep(NA, 6))
  mk_by <- split(d$mk, d$mk$id); bs <- rep(NA_real_, cfg$B)
  for (b in seq_len(cfg$B)) {
    db <- resample_subjects(d, mk_by); mkb <- prep_mk(db); sjb <- db$subj
    if (is.null(mkb) || nrow(mkb) < 30) next
    bs[b] <- tryCatch(fit_A(mkb, cc_w(mkb, sjb), df = df, tau = cfg$tau)$beta[1], error = function(e) NA)
  }
  c(df = df, se = sd(bs, na.rm = TRUE), gd)
}

run_cell <- function(n, mode, cfg) {
  M <- do.call(rbind, Filter(Negate(is.null),
        mclapply(seq_len(cfg$R), one_rep, n = n, mode = mode, cfg = cfg,
                 mc.cores = cfg$cores, mc.preschedule = FALSE)))
  se <- M[, "se"]; lam <- M[, "lam_min"]; nu <- M[, "n_unsupp"]
  ok <- is.finite(se) & is.finite(lam) & se > 0
  se <- se[ok]; lam <- lam[ok]; nu <- nu[ok]; effr <- M[ok, "eff_rank"]; msv <- M[ok, "min_sv"]
  med <- median(se)
  big <- se > 3 * med                                  # "large SE" datasets
  rd  <- lam < quantile(lam, 0.10)                     # "most rank-deficient" 10%
  # 2x2 coincidence of large-SE and rank-deficient
  P_big_given_rd <- if (sum(rd) > 0) mean(big[rd]) else NA
  P_rd_given_big <- if (sum(big) > 0) mean(rd[big]) else NA
  corr_se_lam <- suppressWarnings(cor(log(se), -log(pmax(lam, 1e-300))))
  corr_se_unsupp <- suppressWarnings(cor(log(se), nu))
  data.frame(mode = mode, n = n, J = median(M[, "df"]), R = length(se),
    se_med = med, se_max = max(se), frac_se_gt3med = mean(big),
    lam_min_med = median(lam), lam_min_min = min(lam),
    eff_rank_med = median(effr), ncol = M[1, "ncol"],
    min_sv_med = median(msv), n_unsupp_med = median(nu), n_unsupp_max = max(nu),
    corr_logSE_neglogLam = corr_se_lam, corr_logSE_unsupp = corr_se_unsupp,
    P_bigSE_given_rankdef = P_big_given_rd, P_rankdef_given_bigSE = P_rd_given_big,
    n_bigSE = sum(big), n_bigSE_and_rankdef = sum(big & rd), row.names = NULL)
}

main <- function(cfg = CFG, out = "results/sievedist2.csv") {
  dir.create("results", showWarnings = FALSE); rows <- list(); k <- 1
  for (mode in c("exact", "grow")) for (n in cfg$ns) {
    cat(sprintf("[%s] sievedist2 %s n=%d\n", format(Sys.time(), "%H:%M:%S"), mode, n)); flush.console()
    rows[[k]] <- run_cell(n, mode, cfg); k <- k + 1
    write.csv(do.call(rbind, rows), out, row.names = FALSE)
  }
  cat("done ->", out, "\n")
}
if (sys.nframe() == 0) main()
