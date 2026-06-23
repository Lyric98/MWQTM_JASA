# ============================================================================
# run_nonph.R — validate that the onset mechanism does not enter beta_tau:
# rerun the estimator under a NON-proportional-hazards (log-normal AFT) onset.
# If the claim holds, bias stays ~0 and bootstrap coverage stays near nominal,
# just as under Cox-Weibull. Writes results/nonph.csv.
# ============================================================================
suppressMessages(library(parallel))
source("R/run_one.R"); source("R/bootstrap.R")

CFG <- list(
  ns = c(500, 1000, 2000), scenarios = c(1, 2),
  taus = c(0.25, 0.50, 0.75), R = 200, a_n = 0.05, df_grid = c(3, 4, 5),
  boot_subset = 60, B = 150,
  cores = { sc <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK","")); if (is.na(sc)||sc<1) 4L else sc },
  seed = 20260614
)
aft_params <- function(){ p <- default_params(); p$onset <- "aft"; p }

one_rep <- function(r, n, scenario, tau, cfg) {
  set.seed(cfg$seed + r*7919 + scenario*101 + round(tau*1000))
  p <- aft_params(); p$sig1 <- if (scenario==1) 0 else 0.4
  d <- generate_data(n, scenario = scenario, params = p)
  mk <- prep_mk(d); subj <- d$subj
  if (nrow(mk) < 30) return(NULL)
  gc <- fit_GC(subj); w <- attach_weights(mk, subj, ipcw_weights(subj, GC_at(gc)))
  df <- tryCatch(select_df_cv(mk, subj, w, tau, df_grid=cfg$df_grid), error=function(e) 4)
  fitcv <- tryCatch(fit_A(mk, w, df=df, tau=tau), error=function(e) NULL)
  if (is.null(fitcv)) return(NULL)
  bos <- tryCatch(one_step(mk, w, fitcv, tau, cfg$a_n), error=function(e) NULL)
  out <- list(b_acv = fitcv$beta[1], b_bos = if (is.null(bos)) NA else bos$beta[1])
  if (r <= cfg$boot_subset)
    out$boot <- tryCatch(boot_se(d, tau, B=cfg$B, a_n=cfg$a_n, df_fixed=df,
                                 do_frozen=FALSE)$full, error=function(e) NULL)
  out
}

run_cell <- function(n, scenario, tau, cfg) {
  reps <- Filter(Negate(is.null), mclapply(seq_len(cfg$R), one_rep, n=n,
            scenario=scenario, tau=tau, cfg=cfg, mc.cores=cfg$cores))
  tb <- true_beta(tau, {p<-default_params(); p$sig1<-if(scenario==1)0 else 0.4; p})["beta1"]
  b_acv <- sapply(reps, function(x) x$b_acv); b_bos <- sapply(reps, function(x) x$b_bos)
  bl <- lapply(reps, function(x) x$boot); has <- !sapply(bl, is.null)
  Bf <- if (any(has)) do.call(rbind, bl[has]) else matrix(NA,0,3)
  idx <- which(has)
  cov <- function(bh, se) if (length(se)) mean(abs(bh-tb) <= 1.96*se, na.rm=TRUE) else NA
  data.frame(n=n, scenario=scenario, tau=tau, true_b1=as.numeric(tb), R=length(reps),
    bias_acv=mean(b_acv,na.rm=T)-tb, esd_acv=sd(b_acv,na.rm=T),
    bias_bos=mean(b_bos,na.rm=T)-tb,
    cov_acv_boot=if(nrow(Bf))cov(b_acv[idx],Bf[,1]) else NA,
    cov_bos_boot=if(nrow(Bf))cov(b_bos[idx],Bf[,3]) else NA)
}

main <- function(cfg=CFG, out="results/nonph.csv"){
  dir.create("results", showWarnings=FALSE); rows<-list(); k<-1
  for (sc in cfg$scenarios) for (n in cfg$ns) for (tau in cfg$taus){
    cat(sprintf("[%s] AFT sc=%d n=%d tau=%.2f\n", format(Sys.time(),"%H:%M:%S"),sc,n,tau)); flush.console()
    rows[[k]] <- run_cell(n,sc,tau,cfg); k<-k+1
    write.csv(do.call(rbind,rows), out, row.names=FALSE)
  }
  cat("done ->", out, "\n")
}
if (sys.nframe()==0) main()
