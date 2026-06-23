# ============================================================================
# run_compare.R — confirm chi^C = 0 (first order) and the role of IPCW.
#
# Three estimators of beta_{tau,1}, identical sieve/df, differing only in weight:
#   CC      : complete-event, unweighted   a_i = Delta^T
#   IPCW    : estimated censoring weights   a_i = Delta^T / hat G_C(W|H)
#   ORACLE  : true censoring weights        a_i = Delta^T / G_{C,0}(W|H)
#
# Outputs per (n,scenario,tau):
#   bias/ESD for CC, IPCW, ORACLE  -> (a) all consistent? (b) efficiency differ?
#   sd(beta_IPCW - beta_ORACLE)    -> the chi^C diagnostic: if G_C estimation
#                                     entered at first order this would be
#                                     Theta(n^{-1/2}); Neyman-orthogonality
#                                     predicts it shrinks FASTER than n^{-1/2}
#                                     (a nuisance-product / second-order term).
# ============================================================================
suppressMessages(library(parallel))
source("R/run_one.R")

CFG <- list(
  ns = c(500, 1000, 2000, 4000),     # extra n to see the rate of the gap
  scenarios = c(1, 2),
  taus = c(0.10, 0.50, 0.90),
  regimes = c("base", "hard"),        # #16: also a HARD-censoring stress regime
  R = 400,
  df_grid = c(3, 4, 5),
  cores = { sc <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "")); if (is.na(sc)||sc<1) 4L else sc },
  seed = 20260614
)

default_params2 <- function(scenario){ p<-default_params(); p$sig1<-if(scenario==1)0 else 0.4; p }
# #16 hard-censoring regime: heavier, more covariate-dependent loss to follow-up
# + wider/later entry (more truncation) -> heavier IPCW weight tails, lower ESS,
# approaching the positivity boundary.
params_regime <- function(scenario, regime){
  p <- default_params2(scenario)
  if (regime == "hard"){ p$rate_C <- 0.24; p$A_lo <- -1.0; p$A_hi <- 9.0 }
  p
}

one_rep <- function(r, n, scenario, tau, regime, cfg) {
  set.seed(cfg$seed + r*7919 + scenario*101 + round(tau*1000) + (regime=="hard")*53)
  p <- params_regime(scenario, regime)
  d <- generate_data(n, scenario = scenario, params = p)
  subj <- d$subj; mk <- prep_mk(d)
  if (nrow(mk) < 30) return(NULL)
  # weights; track Cox G_C fit failure and the IPCW weight distribution
  gc <- tryCatch(fit_GC(subj), error = function(e) NULL)
  cox_fail <- is.null(gc)
  a_est <- if (cox_fail) as.numeric(subj$dT == 1) else ipcw_weights(subj, GC_at(gc))
  a_orc <- ifelse(subj$dT == 1, 1/true_GC(subj, p), 0)
  a_cc  <- as.numeric(subj$dT == 1)
  wpos  <- a_est[a_est > 0]
  ess   <- if (length(wpos)) (sum(wpos))^2 / sum(wpos^2) else NA
  w_est <- attach_weights(mk, subj, a_est)
  w_orc <- attach_weights(mk, subj, a_orc)
  w_cc  <- attach_weights(mk, subj, a_cc)
  df <- tryCatch(select_df_cv(mk, subj, w_est, tau, df_grid = cfg$df_grid),
                 error = function(e) 4)
  bx <- function(w) tryCatch(unname(fit_A(mk, w, df = df, tau = tau)$beta[1]),
                             error = function(e) NA_real_)
  c(cc = bx(w_cc), ipcw = bx(w_est), orc = bx(w_orc),
    wq50 = median(wpos), wq95 = quantile(wpos, .95, names=FALSE),
    wq99 = quantile(wpos, .99, names=FALSE), wmax = max(wpos),
    ess_ratio = ess / length(wpos), cox_fail = as.numeric(cox_fail),
    p_cens = unname(d$summary["p_cens"]), trunc_frac = unname(d$summary["trunc_frac"]))
}

run_cell <- function(n, scenario, tau, regime, cfg) {
  reps <- mclapply(seq_len(cfg$R), one_rep, n=n, scenario=scenario, tau=tau,
                   regime=regime, cfg=cfg, mc.cores=cfg$cores)
  M <- do.call(rbind, Filter(Negate(is.null), reps))
  tb <- true_beta(tau, params_regime(scenario, regime))["beta1"]
  mm <- function(c) mean(M[,c], na.rm=T)
  data.frame(
    regime=regime, n=n, scenario=scenario, tau=tau, true_b1=as.numeric(tb), R=nrow(M),
    bias_cc=mm("cc")-tb, esd_cc=sd(M[,"cc"],na.rm=T),
    bias_ipcw=mm("ipcw")-tb, esd_ipcw=sd(M[,"ipcw"],na.rm=T),
    bias_orc=mm("orc")-tb, esd_orc=sd(M[,"orc"],na.rm=T),
    # chi^C diagnostic: dispersion of the estimated-vs-oracle-weight gap
    sd_gap=sd(M[,"ipcw"]-M[,"orc"],na.rm=T), mean_gap=mm("ipcw")-mm("orc"),
    # #16 weight-distribution / positivity diagnostics
    w_med=mm("wq50"), w_p95=mm("wq95"), w_p99=mm("wq99"), w_max=mm("wmax"),
    ess_over_nT=mm("ess_ratio"), cox_fail_rate=mm("cox_fail"),
    p_cens=mm("p_cens"), trunc_frac=mm("trunc_frac")
  )
}

main <- function(cfg=CFG, out="results/compare.csv"){
  dir.create("results", showWarnings=FALSE)
  rows <- list(); k <- 1
  for (rg in cfg$regimes) for (sc in cfg$scenarios) for (n in cfg$ns) for (tau in cfg$taus){
    cat(sprintf("[%s] %s sc=%d n=%d tau=%.2f\n", format(Sys.time(),"%H:%M:%S"), rg,sc,n,tau)); flush.console()
    rows[[k]] <- run_cell(n,sc,tau,rg,cfg); k <- k+1
    write.csv(do.call(rbind,rows), out, row.names=FALSE)
  }
  cat("done ->", out, "\n")
}
if (sys.nframe()==0) main()
