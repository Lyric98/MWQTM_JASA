# ============================================================================
# run_ordcenter.R — demonstrate that density weighting is needed (reply #12).
# Compare two one-step estimators that differ ONLY in the centering g:
#   DW  : density-weighted projection  g_0 = E[a f X|u,t]/E[a f|u,t]  (correct)
#   ORD : ordinary (unweighted)        g_ord = E[a X|u,t]/E[a|u,t]    (f set const)
# In Scenario 1 (homoscedastic) f does not depend on X so DW=ORD and both are
# fine; in Scenario 2 (f depends on X1) ORD is NOT Neyman-orthogonal to the
# surface, so surface-estimation error leaks into beta -> bias/undercoverage.
# Writes results/ordcenter.csv.
# ============================================================================
suppressMessages(library(parallel))
source("R/run_one.R")

CFG <- list(ns=c(500,1000,2000), scenarios=c(1,2), taus=c(0.25,0.50,0.75),
  R=400, a_n=0.05, df_grid=c(3,4,5),
  cores={sc<-as.integer(Sys.getenv("SLURM_CPUS_PER_TASK",""));if(is.na(sc)||sc<1)4L else sc},
  seed=20260614)

one_rep <- function(r, n, scenario, tau, cfg) {
  set.seed(cfg$seed + r*7919 + scenario*101 + round(tau*1000))
  p <- default_params(); p$sig1 <- if (scenario==1) 0 else 0.4
  d <- generate_data(n, scenario=scenario, params=p); mk <- prep_mk(d); subj <- d$subj
  if (nrow(mk) < 30) return(NULL)
  gc <- fit_GC(subj); w <- attach_weights(mk, subj, ipcw_weights(subj, GC_at(gc)))
  df <- tryCatch(select_df_cv(mk, subj, w, tau, df_grid=cfg$df_grid), error=function(e) 4)
  fitcv <- tryCatch(fit_A(mk, w, df=df, tau=tau), error=function(e) NULL)
  if (is.null(fitcv)) return(NULL)
  X <- cbind(mk$X1, mk$X2)
  qhat <- predict_A(fitcv, mk$U, mk$T, mk$X1, mk$X2); psi <- tau - ((mk$L-qhat)<=0)
  # ORACLE density: isolates the CENTERING effect from quantile-spacing f-hat
  # instability. In Scenario 1 (sig1=0) f is constant => g_dw should EXACTLY
  # equal g_ord (sanity); in Scenario 2 (sig1>0) f depends on X1 => g_dw != g_ord.
  fhat <- dnorm(qnorm(tau)) / sigma_X(mk$X1, p$sig0, p$sig1)
  nsub <- length(unique(mk$id))
  # one-step with a given centering ghat and bread density fb
  onestep <- function(ghat, fb){
    Xc <- X - ghat
    Uvec <- colSums(w*Xc*psi)/nsub
    S <- matrix(0,2,2); for(k in seq_len(nrow(mk))) S<-S+w[k]*fb[k]*tcrossprod(Xc[k,]); S<-S/nsub
    beta <- (fitcv$beta + solve(S,Uvec))[1]
    os <- list(w_row=w, Xc=Xc, psi=psi, S=S)
    se <- sandwich_se(mk, os)[1]
    c(unname(beta), se)
  }
  g_dw  <- predict_projection(fit_projection(mk$U,mk$T,X,w*fhat,df_g=4), mk$U,mk$T)
  g_ord <- predict_projection(fit_projection(mk$U,mk$T,X,w,     df_g=4), mk$U,mk$T)
  dw <- onestep(g_dw, fhat); ord <- onestep(g_ord, fhat)
  c(b_dw=dw[1], se_dw=dw[2], b_ord=ord[1], se_ord=ord[2])
}

run_cell <- function(n, scenario, tau, cfg){
  reps <- do.call(rbind, Filter(Negate(is.null), mclapply(seq_len(cfg$R), one_rep,
            n=n, scenario=scenario, tau=tau, cfg=cfg, mc.cores=cfg$cores)))
  tb <- true_beta(tau, {p<-default_params();p$sig1<-if(scenario==1)0 else 0.4;p})["beta1"]
  cov <- function(b,se) mean(abs(b-tb)<=1.96*se, na.rm=TRUE)
  # PAIRED (same-rep) coverage indicators for DW vs ORD (review #16): the paired
  # difference has much smaller MC error than two independent coverages.
  I_dw  <- abs(reps[,"b_dw"]-tb)  <= 1.96*reps[,"se_dw"]
  I_ord <- abs(reps[,"b_ord"]-tb) <= 1.96*reps[,"se_ord"]
  ok <- is.finite(I_dw) & is.finite(I_ord); d <- I_dw[ok]-I_ord[ok]
  bdisc <- sum(I_dw[ok] & !I_ord[ok]); cdisc <- sum(!I_dw[ok] & I_ord[ok])
  mcp <- function(b,c){n<-b+c; if(n==0) NA else min(1,2*pbinom(min(b,c),n,0.5))}
  data.frame(n=n, scenario=scenario, tau=tau, true_b1=as.numeric(tb), R=nrow(reps),
    bias_dw =mean(reps[,"b_dw"], na.rm=T)-tb, esd_dw =sd(reps[,"b_dw"], na.rm=T),
    cov_dw =cov(reps[,"b_dw"], reps[,"se_dw"]),
    bias_ord=mean(reps[,"b_ord"],na.rm=T)-tb, esd_ord=sd(reps[,"b_ord"],na.rm=T),
    cov_ord=cov(reps[,"b_ord"],reps[,"se_ord"]),
    cov_diff=mean(d), cov_diff_mcse=sd(d)/sqrt(length(d)),
    mcnemar_b=bdisc, mcnemar_c=cdisc, mcnemar_p=mcp(bdisc,cdisc))
}

main <- function(cfg=CFG, out="results/ordcenter.csv"){
  dir.create("results",showWarnings=FALSE); rows<-list(); k<-1
  for (sc in cfg$scenarios) for (n in cfg$ns) for (tau in cfg$taus){
    cat(sprintf("[%s] ord sc=%d n=%d tau=%.2f\n", format(Sys.time(),"%H:%M:%S"),sc,n,tau)); flush.console()
    rows[[k]] <- run_cell(n,sc,tau,cfg); k<-k+1; write.csv(do.call(rbind,rows),out,row.names=FALSE)
  }
  cat("done ->",out,"\n")
}
if (sys.nframe()==0) main()
