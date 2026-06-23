# ============================================================================
# run_sieve.R (review step 7) — two asymptotic checks for the DEFAULT
# complete-event (omega=1) A-CV estimator of beta_{tau,1}, Scenario 2, tau=0.5.
#   EXACT-sieve  : m_0 bilinear (in the cubic-spline span) => zero approximation
#                  bias; fixed df. Tests score scaling / root-n stochastic part.
#                  Pass: sqrt(n)*bias bounded, sqrt(n)*ESD ~ constant.
#   GROWING-sieve: m_0 smooth (not in any fixed span); J_n grows with n,
#                  J_n = ceil(c n^kappa). Tests that the approximation bias b_K
#                  is controlled. Pass: sqrt(n)*|bias| does not blow up (ideally
#                  decreases).
# Writes results/sieve.csv.
# ============================================================================
suppressMessages(library(parallel)); source("R/run_one.R")

CFG <- list(ns=c(500,1000,2000,4000,8000), R=300, tau=0.50, scenario=2,
  kappa=0.18, cJ=2.0,   # growing-sieve: J_n = ceil(cJ * n^kappa)
  df_exact=4,           # cubic spline df per margin (contains bilinear exactly)
  cores={sc<-as.integer(Sys.getenv("SLURM_CPUS_PER_TASK",""));if(is.na(sc)||sc<1)4L else sc},
  seed=20260615)

cc_beta1 <- function(d, tau, df){
  mk<-prep_mk(d); subj<-d$subj; if(nrow(mk)<30) return(NA_real_)
  w<-attach_weights(mk,subj,as.numeric(subj$dT==1))   # complete-event omega=1
  tryCatch(unname(fit_A(mk,w,df=df,tau=tau)$beta[1]), error=function(e) NA_real_)
}

one_rep <- function(r,n,mode,cfg){
  set.seed(cfg$seed+r*7919+(mode=="grow")*5)
  p<-default_params(); p$sig1<-0.4
  p$surface<-if(mode=="exact")"bilinear" else "smooth"
  df<-if(mode=="exact") cfg$df_exact else max(3, ceiling(cfg$cJ*n^cfg$kappa))
  d<-generate_data(n,scenario=cfg$scenario,params=p)
  c(beta=cc_beta1(d,cfg$tau,df), df=df)
}

run_cell <- function(n,mode,cfg){
  M<-do.call(rbind,mclapply(seq_len(cfg$R),one_rep,n=n,mode=mode,cfg=cfg,mc.cores=cfg$cores))
  tb<-true_beta(cfg$tau,{p<-default_params();p$sig1<-0.4;p})["beta1"]
  b<-M[,"beta"]; bias<-mean(b,na.rm=T)-tb; esd<-sd(b,na.rm=T)
  data.frame(mode=mode, n=n, df=median(M[,"df"]), R=sum(is.finite(b)),
    bias=bias, sqrtn_bias=sqrt(n)*bias, esd=esd, sqrtn_esd=sqrt(n)*esd)
}

main<-function(cfg=CFG,out="results/sieve.csv"){
  dir.create("results",showWarnings=FALSE); rows<-list(); k<-1
  for(mode in c("exact","grow")) for(n in cfg$ns){
    cat(sprintf("[%s] sieve %s n=%d\n",format(Sys.time(),"%H:%M:%S"),mode,n)); flush.console()
    rows[[k]]<-run_cell(n,mode,cfg); k<-k+1; write.csv(do.call(rbind,rows),out,row.names=FALSE)
  }
  cat("done ->",out,"\n")
}
if(sys.nframe()==0) main()
