# ============================================================================
# run_diag.R — variance ladder for the plug-in sandwich SER ~ 0.7.
# For each replicate (B-OS point fixed at the A-CV fit), recompute the clustered
# sandwich SE under combinations of oracle vs plug-in (f, g), plus a row-level
# (non-clustered) meat and a Jacobian-form bread, and compare mean SE to the
# Monte-Carlo ESD. Also report R_f = f_hat/f_true. Localizes the undercoverage.
# Writes results/diag.csv.
# ============================================================================
suppressMessages(library(parallel))
source("R/run_one.R")

CFG <- list(ns = c(500, 1000, 2000), scenarios = c(1, 2), taus = c(0.25, 0.50, 0.75),
  R = 500, a_n = 0.05, df_grid = c(3,4,5), df_fixed = 4,
  cores = { sc<-as.integer(Sys.getenv("SLURM_CPUS_PER_TASK","")); if(is.na(sc)||sc<1)4L else sc },
  seed = 20260614)

# sandwich SE from centered covariate Xc, density f (vectors over rows), psi,
# weights w, subject ids. bread_jac: use Jacobian form a f (X-g) X' instead of
# symmetric; row_meat: skip subject aggregation.
se_from <- function(Xc, f, psi, w, ids, Xfull, bread_jac = FALSE, row_meat = FALSE) {
  nsub <- length(unique(ids))
  # bread
  S <- matrix(0,2,2)
  for (k in seq_along(w)) {
    A <- if (bread_jac) tcrossprod(Xc[k,], Xfull[k,]) else tcrossprod(Xc[k,])
    S <- S + w[k]*f[k]*A
  }
  S <- S/nsub
  contrib <- w * Xc * psi
  if (row_meat) {
    V <- crossprod(contrib)/nsub
  } else {
    Z <- rowsum(contrib, group=ids); V <- crossprod(as.matrix(Z))/nsub
  }
  Si <- solve(S); Sig <- Si %*% V %*% t(Si) / nsub
  sqrt(diag(Sig))[1]   # X1 component
}

one_rep <- function(r, n, scenario, tau, oracles, cfg) {
  set.seed(cfg$seed + r*7919 + scenario*101 + round(tau*1000))
  p <- default_params(); p$sig1 <- if (scenario==1) 0 else 0.4
  d <- generate_data(n, scenario=scenario, params=p)
  mk <- prep_mk(d); subj <- d$subj
  if (nrow(mk) < 30) return(NULL)
  gc <- fit_GC(subj); w <- attach_weights(mk, subj, ipcw_weights(subj, GC_at(gc)))
  df <- tryCatch(select_df_cv(mk, subj, w, tau, df_grid=cfg$df_grid), error=function(e) 4)
  fitcv <- tryCatch(fit_A(mk, w, df=df, tau=tau), error=function(e) NULL)
  if (is.null(fitcv)) return(NULL)
  X <- cbind(mk$X1, mk$X2)
  # residuals/psi from A-CV surface (independent of f,g)
  qhat <- predict_A(fitcv, mk$U, mk$T, mk$X1, mk$X2); psi <- tau - ((mk$L-qhat)<=0)
  # f: est (quantile spacing) vs true
  fit_lo <- fit_A(mk, w, df=df, tau=tau-cfg$a_n); fit_hi <- fit_A(mk, w, df=df, tau=tau+cfg$a_n)
  f_est <- sparsity_fhat(fit_lo, fit_hi, mk$U, mk$T, mk$X1, mk$X2, cfg$a_n)
  f_tru <- dnorm(qnorm(tau)) / sigma_X(mk$X1, p$sig0, p$sig1)
  # g: est (feasible projection, df matched to the SELECTED point-estimator df)
  #    vs FINITE-SIEVE oracle g_J at the SAME df (review #4): the point
  #    estimator's implicit centering lives in the df-dim spline space, so the
  #    oracle rung must project in that same space (true f, pop weights), NOT
  #    the infinite-dim pool projection -- else bread and pseudo-truth mismatch.
  proj_e <- fit_projection(mk$U, mk$T, X, w*f_est, df_g=df)
  g_est <- predict_projection(proj_e, mk$U, mk$T)
  orc <- oracles[[as.character(df)]]; if (is.null(orc)) orc <- oracles[[1]]
  g_tru <- predict_projection(orc$proj, mk$U, mk$T)
  Xc_e <- X - g_est; Xc_t <- X - g_tru
  # one-step point (est f, est g)
  Uvec <- colSums(w*Xc_e*psi)/length(unique(mk$id))
  S_e <- matrix(0,2,2); for(k in seq_len(nrow(mk))) S_e<-S_e+w[k]*f_est[k]*tcrossprod(Xc_e[k,]); S_e<-S_e/length(unique(mk$id))
  b_bos <- (fitcv$beta + solve(S_e, Uvec))[1]
  ids <- mk$id
  # fixed-J rung: refit A at deterministic df, recompute psi/f/g/sandwich, to
  # separate CV-tuning randomness from the plug-in sandwich gap (#6).
  se_fixedJ <- NA; b_fixedJ <- NA
  if (cfg$df_fixed != df) {
    fitf <- tryCatch(fit_A(mk, w, df=cfg$df_fixed, tau=tau), error=function(e) NULL)
    if (!is.null(fitf)) {
      qf <- predict_A(fitf, mk$U, mk$T, mk$X1, mk$X2); psif <- tau-((mk$L-qf)<=0)
      flo<-fit_A(mk,w,df=cfg$df_fixed,tau=tau-cfg$a_n); fhi<-fit_A(mk,w,df=cfg$df_fixed,tau=tau+cfg$a_n)
      ff <- sparsity_fhat(flo,fhi,mk$U,mk$T,mk$X1,mk$X2,cfg$a_n)
      gf <- predict_projection(fit_projection(mk$U,mk$T,X,w*ff,df_g=4), mk$U,mk$T)
      se_fixedJ <- se_from(X-gf, ff, psif, w, ids, X); b_fixedJ <- fitf$beta[1]
    }
  } else { se_fixedJ <- se_from(Xc_e,f_est,psi,w,ids,X); b_fixedJ <- fitcv$beta[1] }
  list(
    b_acv = fitcv$beta[1], b_bos = unname(b_bos), b_fixedJ = unname(b_fixedJ),
    df_sel = df, nT = length(unique(ids)), nrow = nrow(mk),
    Rf_med = median(f_est/f_tru),
    se_plug   = se_from(Xc_e, f_est, psi, w, ids, X),
    se_ftrue  = se_from(Xc_e, f_tru, psi, w, ids, X),
    se_gtrue  = se_from(Xc_t, f_est, psi, w, ids, X),
    se_oracle = se_from(Xc_t, f_tru, psi, w, ids, X),
    se_jac    = se_from(Xc_e, f_est, psi, w, ids, X, bread_jac=TRUE),
    se_row    = se_from(Xc_e, f_est, psi, w, ids, X, row_meat=TRUE),
    se_fixedJ = se_fixedJ
  )
}

run_cell <- function(n, scenario, tau, cfg) {
  # finite-sieve oracle projections at each candidate df (true f, pop weights)
  oracles <- lapply(cfg$df_grid, function(dfg) build_oracle(scenario, n_pool=8000, df_g=dfg))
  names(oracles) <- as.character(cfg$df_grid)
  reps <- Filter(Negate(is.null), mclapply(seq_len(cfg$R), one_rep, n=n,
            scenario=scenario, tau=tau, oracles=oracles, cfg=cfg, mc.cores=cfg$cores))
  g <- function(f) sapply(reps, f)
  esd <- sd(g(function(x) x$b_bos), na.rm=TRUE)
  esd_fj <- sd(g(function(x) x$b_fixedJ), na.rm=TRUE)
  m <- function(f) mean(g(f), na.rm=TRUE)
  dfsel <- g(function(x) x$df_sel)
  data.frame(n=n, scenario=scenario, tau=tau, R=length(reps),
    nT_over_n = m(function(x) x$nT)/n, Rf_med = m(function(x) x$Rf_med),
    df_mode = as.integer(names(which.max(table(dfsel)))),
    p_df3 = mean(dfsel==3), p_df4 = mean(dfsel==4), p_df5 = mean(dfsel==5),
    ESD = esd, ESD_fixedJ = esd_fj,
    SER_plug   = m(function(x) x$se_plug)/esd,
    SER_ftrue  = m(function(x) x$se_ftrue)/esd,
    SER_gtrue  = m(function(x) x$se_gtrue)/esd,
    SER_oracle = m(function(x) x$se_oracle)/esd,
    SER_jac    = m(function(x) x$se_jac)/esd,
    SER_fixedJ = m(function(x) x$se_fixedJ)/esd_fj,
    ratio_row_over_cluster = m(function(x) x$se_row)/m(function(x) x$se_plug))
}

main <- function(cfg=CFG, out="results/diag.csv"){
  dir.create("results", showWarnings=FALSE); rows<-list(); k<-1
  for (sc in cfg$scenarios) for (n in cfg$ns) for (tau in cfg$taus){
    cat(sprintf("[%s] diag sc=%d n=%d tau=%.2f\n", format(Sys.time(),"%H:%M:%S"),sc,n,tau)); flush.console()
    rows[[k]] <- run_cell(n,sc,tau,cfg); k<-k+1
    write.csv(do.call(rbind,rows), out, row.names=FALSE)
  }
  cat("done ->", out, "\n")
}
if (sys.nframe()==0) main()
