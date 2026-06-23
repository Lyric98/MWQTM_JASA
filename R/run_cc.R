# ============================================================================
# run_cc.R (review step 6/#8) — bootstrap coverage of the DEFAULT complete-event
# (CC) A-CV estimator, alongside the IPCW A-CV estimator, by (n, tau).
#   CC   : omega = 1            (weight a_i = Delta^T)
#   IPCW : omega = 1/hat G_C    (weight a_i = Delta^T / hat G_C)
# Reports per cell: bias, ESD, mean bootstrap SE, bootstrap coverage, for both.
# Writes results/cc_coverage.csv.
# ============================================================================
suppressMessages(library(parallel)); source("R/run_one.R")

CFG <- list(ns=c(500,1000,2000), scenarios=c(1,2), taus=c(0.10,0.50,0.90),
  R=200, B=200, boot_subset=150, df_grid=c(3,4,5),
  cores={sc<-as.integer(Sys.getenv("SLURM_CPUS_PER_TASK",""));if(is.na(sc)||sc<1)4L else sc},
  seed=20260615)

# weight vector for a fitted subject table; type in {"cc","ipcw"}
wvec <- function(subj, type){
  if (type=="cc") as.numeric(subj$dT==1)
  else ipcw_weights(subj, GC_at(fit_GC(subj)))
}

acv_beta1 <- function(mk, subj, type, df){
  w <- attach_weights(mk, subj, wvec(subj,type))
  tryCatch(unname(fit_A(mk,w,df=df,tau=attr(mk,"tau"))$beta[1]), error=function(e) NA_real_)
}

# subject bootstrap SE for CC and IPCW A-CV (full refit each resample)
boot_cc <- function(d, tau, df, B){
  mk_by <- split(d$mk, d$mk$id); bm <- matrix(NA,B,2)
  for(b in 1:B){
    bid <- sample(d$subj$id, replace=TRUE)
    ns  <- d$subj[match(bid,d$subj$id),]; ns$id <- seq_len(nrow(ns))
    parts <- vector("list",length(bid))
    for(j in seq_along(bid)){r<-mk_by[[as.character(bid[j])]]; if(!is.null(r)&&nrow(r)){r$id<-j; parts[[j]]<-r}}
    mkb <- do.call(rbind,parts); mkb<-mkb[mkb$dT==1 & mkb$U>0 & mkb$U<mkb$T,]
    if(is.null(mkb)||nrow(mkb)<30) next
    attr(mkb,"tau")<-tau
    bm[b,1]<-acv_beta1(mkb,ns,"cc",df); bm[b,2]<-acv_beta1(mkb,ns,"ipcw",df)
  }
  apply(bm,2,sd,na.rm=TRUE)
}

one_rep <- function(r,n,scenario,tau,cfg){
  set.seed(cfg$seed+r*7919+scenario*101+round(tau*1000))
  p<-default_params(); p$sig1<-if(scenario==1)0 else 0.4
  d<-generate_data(n,scenario=scenario,params=p); mk<-prep_mk(d); subj<-d$subj
  if(nrow(mk)<30) return(NULL); attr(mk,"tau")<-tau
  w_ip<-attach_weights(mk,subj,wvec(subj,"ipcw"))
  df<-tryCatch(select_df_cv(mk,subj,w_ip,tau,df_grid=cfg$df_grid),error=function(e)4)
  out<-list(b_cc=acv_beta1(mk,subj,"cc",df), b_ip=acv_beta1(mk,subj,"ipcw",df))
  if(r<=cfg$boot_subset) out$bse<-tryCatch(boot_cc(d,tau,df,cfg$B),error=function(e)c(NA,NA))
  out
}

run_cell <- function(n,scenario,tau,cfg){
  reps<-Filter(Negate(is.null),mclapply(seq_len(cfg$R),one_rep,n=n,scenario=scenario,tau=tau,cfg=cfg,mc.cores=cfg$cores))
  tb<-true_beta(tau,{p<-default_params();p$sig1<-if(scenario==1)0 else 0.4;p})["beta1"]
  bcc<-sapply(reps,function(x)x$b_cc); bip<-sapply(reps,function(x)x$b_ip)
  bl<-lapply(reps,function(x)x$bse); has<-!sapply(bl,is.null); B<-if(any(has))do.call(rbind,bl[has])else matrix(NA,0,2); idx<-which(has)
  covf<-function(bh,se)if(length(se))mean(abs(bh-tb)<=1.96*se,na.rm=TRUE)else NA
  data.frame(n=n,scenario=scenario,tau=tau,true_b1=as.numeric(tb),R=length(reps),
    bias_cc=mean(bcc,na.rm=T)-tb, esd_cc=sd(bcc,na.rm=T),
    bias_ip=mean(bip,na.rm=T)-tb, esd_ip=sd(bip,na.rm=T),
    esd_ratio_cc_ip=sd(bcc,na.rm=T)/sd(bip,na.rm=T),
    cov_cc=if(nrow(B))covf(bcc[idx],B[,1])else NA,
    cov_ip=if(nrow(B))covf(bip[idx],B[,2])else NA)
}

main<-function(cfg=CFG,out="results/cc_coverage.csv"){
  dir.create("results",showWarnings=FALSE); rows<-list(); k<-1
  for(sc in cfg$scenarios) for(n in cfg$ns) for(tau in cfg$taus){
    cat(sprintf("[%s] cc sc=%d n=%d tau=%.2f\n",format(Sys.time(),"%H:%M:%S"),sc,n,tau)); flush.console()
    rows[[k]]<-run_cell(n,sc,tau,cfg); k<-k+1; write.csv(do.call(rbind,rows),out,row.names=FALSE)
  }
  cat("done ->",out,"\n")
}
if(sys.nframe()==0) main()
