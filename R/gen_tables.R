# ============================================================================
# gen_tables.R — regenerate all sim-revision LaTeX tables from results/.
# Produces: table_estimators.tex (4-estimator headline), table_cc.tex (paired
# CC vs IPCW A-CV + McNemar), table_bosc.tex (B-OS-CF one-step correction),
# table_surface.tex (het+stress IMSE/beta/ESS), table_sieve.tex (sieve+coverage),
# table_ladder.tex (full variance ladder).
# ============================================================================
f3 <- function(x) ifelse(is.na(x), "--", formatC(x, format="f", digits=3))
f2 <- function(x) ifelse(is.na(x), "--", formatC(x, format="f", digits=2))
f4 <- function(x) ifelse(is.na(x), "--", formatC(x, format="f", digits=4))
f1 <- function(x) ifelse(is.na(x), "--", formatC(x, format="f", digits=1))
fexp <- function(x) {                       # LaTeX a x 10^b for big condition numbers
  if (length(x) != 1 || is.na(x) || !is.finite(x)) return("--")
  e <- formatC(x, format = "e", digits = 1)
  paste0(sub("e.*", "", e), "$\\times10^{", as.integer(sub(".*e", "", e)), "}$")
}
blk <- function(lines, prev, cur) if (!is.null(prev) && prev != cur) c(lines, "\\addlinespace") else lines

# ---------- 1. consolidated 4-estimator table ------------------------------
cc <- do.call(rbind, lapply(sprintf("results/consolidated_cell_%02d.csv", 1:18), read.csv))
cc <- cc[order(cc$scenario, cc$n, cc$tau), ]
L <- c("\\begin{table}[t]\\centering\\small",
 "\\caption{Headline estimator comparison on a \\emph{single} paired experiment ($R=1000$; bootstrap coverage over a $500$-replicate subset, $B=500$ full-refit resamples; MC SE $\\approx0.010$). Four estimators on the same replicates: complete-event (CC, $\\omega=1$) and IPCW ($\\omega=1/\\widehat G_C$), each as the joint-sieve A-CV point estimator and the cross-fitted one-step \\textbf{B-OS-CF} ($K=5$ folds). All coverages are full-refit subject bootstrap.}",
 "\\label{tab:estsim}\\begin{tabular}{rrr r rrrr rrrr}\\toprule",
 "&&&& \\multicolumn{4}{c}{bias} & \\multicolumn{4}{c}{bootstrap coverage}\\\\",
 "\\cmidrule(lr){5-8}\\cmidrule(lr){9-12}",
 "Sc.&$n$&$\\tau$&$\\beta_{\\tau,1}$& \\tiny CC-ACV & \\tiny IP-ACV & \\tiny CC-CF & \\tiny IP-CF & \\tiny CC-ACV & \\tiny IP-ACV & \\tiny CC-CF & \\tiny IP-CF\\\\\\midrule")
prev <- NULL
for (i in seq_len(nrow(cc))) { x <- cc[i,]; cur <- paste(x$scenario,x$n); L <- blk(L,prev,cur); prev <- cur
  L <- c(L, paste(paste(c(x$scenario,x$n,f2(x$tau),f3(x$true_b1),
    f3(x$bias_acv_cc),f3(x$bias_acv_ip),f3(x$bias_cf_cc),f3(x$bias_cf_ip),
    f3(x$cov_acv_cc),f3(x$cov_acv_ip),f3(x$cov_cf_cc),f3(x$cov_cf_ip)), collapse=" & "), "\\\\")) }
L <- c(L, "\\bottomrule\\end{tabular}\\end{table}")
writeLines(L, "paper/table_estimators.tex")

# ---------- 2. paired CC vs IPCW A-CV with McNemar -------------------------
# two-sided EXACT McNemar (binomial on discordant pairs); b=c => p=1.
# (The earlier continuity-corrected normal approx wrongly gave p<1 when b=c,
#  because |b-c|-1 = -1 there; the exact binomial is correct at small counts.)
mcnemar_p <- function(b, c) { n <- b + c; if (n == 0) return(NA); min(1, 2 * pbinom(min(b, c), n, 0.5)) }
L <- c("\\begin{table}[t]\\centering\\small",
 "\\caption{\\textbf{Paired} CC--A-CV vs.\\ IPCW--A-CV on the same replicates and bootstrap draws ($R=1000$; coverage over $500$; $B=500$). $\\Delta$cov $=\\mathrm{cov}_{\\rm CC}-\\mathrm{cov}_{\\rm IPCW}$; McNemar $(b,c)$ are the discordant counts (CC covers \\& IPCW not; IPCW covers \\& CC not) with the two-sided \\emph{exact} (binomial) McNemar $p$; ESD ratio is CC/IPCW. Of the $18$ cells, \\emph{one} (Sc.~1, $n{=}1000$, $\\tau{=}0.75$: $b{=}6,c{=}0$) gives an unadjusted exact McNemar $p=0.031$; \\textbf{no cell remains significant after a Holm/Bonferroni multiplicity correction} (the next smallest is $p=0.07$). The absolute coverage differences are $|\\Delta\\mathrm{cov}|\\le0.012$ with no systematic direction. At $R_{\\rm cov}=500$ this does not formally establish equivalence at a $0.02$ margin (a paired TOST would, and is left to a larger run) but is consistent with it. Separately, the CC ESD is smaller in every cell---an \\emph{observed} efficiency improvement under these DGPs (not a claim of provably greater accuracy, which a non-significant test cannot support).}",
 "\\label{tab:cc}\\begin{tabular}{rrr rr r rr r r}\\toprule",
 "Sc.&$n$&$\\tau$&cov$_{\\rm CC}$&cov$_{\\rm IPCW}$&$\\Delta$cov&$b$&$c$&$p$&$\\tfrac{\\rm ESD_{CC}}{\\rm ESD_{IPCW}}$\\\\\\midrule")
prev <- NULL
for (i in seq_len(nrow(cc))) { x <- cc[i,]; cur <- paste(x$scenario,x$n); L <- blk(L,prev,cur); prev <- cur
  p <- mcnemar_p(x$mcnemar_b_acv, x$mcnemar_c_acv)
  L <- c(L, paste(paste(c(x$scenario,x$n,f2(x$tau),f3(x$cov_acv_cc),f3(x$cov_acv_ip),
    f3(x$cov_acv_cc-x$cov_acv_ip),x$mcnemar_b_acv,x$mcnemar_c_acv,f2(p),f3(x$esd_ratio_acv)),collapse=" & "),"\\\\")) }
L <- c(L, "\\bottomrule\\end{tabular}\\end{table}")
writeLines(L, "paper/table_cc.tex")

# ---------- 3. B-OS-CF one-step correction ---------------------------------
L <- c("\\begin{table}[t]\\centering\\small",
 "\\caption{Cross-fitted one-step \\textbf{B-OS-CF} ($K=5$), CC and IPCW members. corr is the per-cell \\emph{mean} signed one-step correction $\\widetilde\\beta_{\\rm CF}-\\widehat\\beta_{\\rm ACV}$; std is the mean of the replicate-wise standardized corrections $\\frac1R\\sum_r\\mathrm{corr}_r/\\widehat{\\mathrm{SE}}_r$ (so $\\mathrm{corr}$ and $\\mathrm{std}$ can differ in sign). Both signed summaries are small ($|\\mathrm{corr}|\\le0.007$, $|\\mathrm{std}|\\le0.05$ in every cell), consistent with the A-CV first-order conditions already nearly annihilating the orthogonal score at $\\lambda=0$. \\emph{Caveat}: a signed mean near zero can mask per-replicate cancellation; the per-replicate \\emph{absolute} correction and its split/one-step decomposition are in Table~\\ref{tab:corr}.}",
 "\\label{tab:bosc}\\begin{tabular}{rrr rr rr}\\toprule",
 "Sc.&$n$&$\\tau$&corr$_{\\rm CC}$&corr$_{\\rm IPCW}$&std$_{\\rm CC}$&std$_{\\rm IPCW}$\\\\\\midrule")
prev <- NULL
for (i in seq_len(nrow(cc))) { x <- cc[i,]; cur <- paste(x$scenario,x$n); L <- blk(L,prev,cur); prev <- cur
  L <- c(L, paste(paste(c(x$scenario,x$n,f2(x$tau),f3(x$corr_cf_cc),f3(x$corr_cf_ip),
    f3(x$std_corr_cc),f3(x$std_corr_ip)),collapse=" & "),"\\\\")) }
L <- c(L, "\\bottomrule\\end{tabular}\\end{table}")
writeLines(L, "paper/table_bosc.tex")

# ---------- 4. surface (het + stress) --------------------------------------
s <- read.csv("results/surface.csv"); s <- s[order(s$paramset, s$n, s$tau), ]
L <- c("\\begin{table}[t]\\centering\\small",
 "\\caption{Surface-heterogeneity DGP $\\sigma(u,t,X)=\\sigma_0+\\sigma_1 X_1+\\sigma_2(t-u)+\\sigma_3(u/t)$, so $m_{\\tau,0}(u,t)$ is \\emph{non-parallel} across $\\tau$ (\\textbf{fitted estimator: CC--A-CV}; $R=500$, coverage over $300$, $B=300$). ``het'' is the base design; ``stress'' adds high truncation, $32\\%$ censoring, $\\sim1.3$ visits/subject. \\textbf{Both} slopes are reported: bias and bootstrap coverage of $\\beta_{\\tau,1}$ \\emph{and} the continuous-covariate $\\beta_{\\tau,2}$. IMSE columns are Monte-Carlo \\emph{means}: IMSE$_{\\rm unm}$ over the whole grid, and---now shown side by side---IMSE$_{\\rm msk_v}$ over the \\emph{visit}-estimable cells and IMSE$_{\\rm msk_s}$ over the \\emph{distinct-subject}-estimable cells. The two local-support masks: ``\\%n.e.\\ visit'' (fewer than $10$ pre-onset \\emph{visits} within $\\pm0.5$ in $(u,t)$) and the surface-relevant ``\\%n.e.\\ subj'' (fewer than $5$ \\emph{distinct subjects}); both are replicate-specific and design-only (not strictly comparable across $n$, and the two thresholds are not information-equivalent---the subject mask is looser by construction). Because visits cluster within subject, the subject mask flags about \\emph{half} as many cells as the visit mask (stress: $\\approx21\\%$ vs.\\ $\\approx41\\%$), yet the two masked IMSEs are nearly identical (columns msk$_{\\rm v}$ vs.\\ msk$_{\\rm s}$, both $\\approx0.10$ under stress), so the two operational masks give \\emph{qualitatively similar} conclusions---not a formal mask-invariance claim, since the thresholds are not information-matched (the subject mask is looser) and each mask is replicate-specific (a fixed common-domain IMSE would be needed to compare convergence across $n$). ESS$^{\\rm ip}/n_T$ is a \\emph{design diagnostic only} (the corresponding IPCW weights; the fitted CC estimator has $w\\equiv1$). $\\beta$ stays broadly covered (finite-sample bias rises somewhat in the smallest stress sample); the low-information boundary surface should be read only on the estimable domain.}",
 "\\label{tab:surface}\\begin{tabular}{lrr rr rr rrr rr r}\\toprule",
 "&&&\\multicolumn{2}{c}{$\\beta_1$}&\\multicolumn{2}{c}{$\\beta_2$}&\\multicolumn{3}{c}{IMSE}&\\multicolumn{2}{c}{\\%n.e.}&\\\\",
 "\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}\\cmidrule(lr){8-10}\\cmidrule(lr){11-12}",
 "set&$n$&$\\tau$&bias&cov&bias&cov&unm&msk$_{\\rm v}$&msk$_{\\rm s}$&visit&subj&$\\tfrac{\\rm ESS^{ip}}{n_T}$\\\\\\midrule")
prev <- NULL
for (i in seq_len(nrow(s))) { x <- s[i,]; cur <- paste(x$paramset,x$n); L <- blk(L,prev,cur); prev <- cur
  L <- c(L, paste(paste(c(x$paramset,x$n,f2(x$tau),f3(x$bias_b1),f3(x$cov_b1_boot),
    f3(x$bias_b2),f3(x$cov_b2_boot),
    f4(x$imse_unmasked),f4(x$imse_masked),f4(x$imse_masked_subj),f2(x$frac_nonest),f2(x$frac_nonest_subj),
    f2(x$ess/x$n_T)),collapse=" & "),"\\\\")) }
L <- c(L, "\\bottomrule\\end{tabular}\\end{table}")
writeLines(L, "paper/table_surface.tex")

# ---------- 5. sieve with coverage + MCSE ----------------------------------
sv <- read.csv("results/sieve_cov.csv"); sv <- sv[order(sv$mode, sv$tau, sv$n), ]
sv$sig <- abs(sv$sqrtn_bias) > 1.96 * sv$mcse_sqrtn_bias    # scaled bias significant?
L <- c("\\begin{table}[t]\\centering\\small",
 "\\caption{Sieve asymptotic check WITH inference, CC joint-sieve with \\emph{deterministic} $J$ (Scenario~2; $R=400$, coverage over $250$, $B=250$; $\\tau\\in\\{0.1,0.5,0.9\\}$). $\\sqrt n\\,$bias is shown with its Monte-Carlo SE; $\\dagger$ marks cells where $|\\sqrt n\\,\\mathrm{bias}|>1.96\\,$MCSE (a \\emph{detectable} scaled bias). SER$=\\overline{\\mathrm{SE}}_{\\rm boot}/\\mathrm{ESD}$ and the mean $95\\%$ interval length disentangle the coverage behaviour. \\textbf{At $n{=}500$ the \\emph{growing}-sieve bootstrap is degenerate} (mean SER $5$--$8$, mean interval length $3$--$5$ vs.\\ an ESD-implied length $\\approx0.6$, \\emph{at all three $\\tau$ including the median}): a few resamples produce enormous $\\widehat{\\mathrm{SE}}$ that dominate the \\emph{mean}, so the near-$0.95$ coverage is an artifact of rare very-long intervals, not valid inference. We therefore mark these cells $\\ddagger$ \\textbf{not interpretable} (n.i.). The distributional audit (Table~\\ref{tab:sievedist}) shows the \\emph{median} bootstrap SE is sound (equal to the exact sieve's) and the blow-up (SE up to $96$) is confined to $2/250$ datasets ($0.8\\%$, with substantial binomial uncertainty); that audit does \\emph{not} attribute these events to the measured Gram-rank, local-support, or onset-diversity diagnostics, and the precise numerical trigger remains \\emph{unidentified} (it audits $\\tau{=}0.5$; the same $n{=}500$ growing-sieve regime shows similar blow-ups at $\\tau{=}0.1,0.9$, but we do not claim an identical mechanism across $\\tau$). The growing sieve is a diagnostic only---the \\emph{recommended} procedure uses the small fixed sieve, where this does not occur. For $n\\ge2000$ the SER is well-behaved, and there the low-bias/low-coverage cells (e.g.\\ growing $\\tau{=}0.1$, $n{=}4000$: $\\sqrt n\\,$bias $0.52$, SER $0.96$, cov $0.896$) cannot be explained by the estimated bias ($b_K$) alone---the scaled bias is only $\\approx0.13$ ESD---and appear to reflect studentization/SE distortion as well.}",
 "\\label{tab:sieve}\\begin{tabular}{lrr r r@{$\\,\\pm\\,$}l r rr r}\\toprule",
 "sieve&$\\tau$&$n$&$J$&\\multicolumn{2}{c}{$\\sqrt n\\,$bias}&$\\sqrt n\\,$ESD&SER&len&cov$^{\\rm boot}$\\\\\\midrule")
prev <- NULL
for (i in seq_len(nrow(sv))) { x <- sv[i,]; cur <- paste(x$mode,x$tau); L <- blk(L,prev,cur); prev <- cur
  bias_str <- paste0(f2(x$sqrtn_bias), if (isTRUE(x$sig)) "$^\\dagger$" else "")
  esd <- x$sqrtn_esd / sqrt(x$n); ser <- x$mean_se / esd; ilen <- 2 * 1.96 * x$mean_se
  ni <- x$mode == "grow" && x$n == 500
  cov_str <- if (ni) "n.i.$^\\ddagger$" else f3(x$cov_boot)
  L <- c(L, paste(paste(c(x$mode,f2(x$tau),x$n,x$J,bias_str,f2(x$mcse_sqrtn_bias),
    f2(x$sqrtn_esd),f2(ser),f3(ilen),cov_str),collapse=" & "),"\\\\")) }
L <- c(L, "\\bottomrule\\end{tabular}\\end{table}")
writeLines(L, "paper/table_sieve.tex")

# ---------- 6. full variance ladder ----------------------------------------
d <- read.csv("results/diag.csv"); d <- d[order(d$scenario, d$n, d$tau), ]
L <- c("\\begin{table}[t]\\centering\\small",
 "\\caption{Variance-ladder for the plug-in clustered sandwich SER ($=\\overline{\\mathrm{SE}}/$Monte-Carlo ESD; $R=500$). $\\widehat f/f_0$ is the per-cell \\emph{median over replicates and rows} of the sparsity ratio. Rungs replace plug-in nuisances by truth one at a time, holding the A-CV point estimator fixed; the $g$-true rung uses the \\emph{finite-sieve} oracle $g_J$ at the point-estimator's df; fix-$J$ is a \\emph{separate} rung that re-selects the point estimator at a deterministic df (normalized by its own ESD). Both the $f$-true ($0.88$--$1.07$) and oracle-$(f,g)$ ($0.89$--$1.06$) rungs are near nominal, while the $g$-true rung \\emph{alone} stays low ($0.54$--$0.95$): replacing the sparsity $\\widehat f$ nearly closes the gap, replacing $g$ alone does not. Sparsity estimation is thus the \\emph{largest identified contributor} in the small- and moderate-$n$ cells, though residual distortion remains and the $f$-true improvement is not monotone in every large-$n$ cell; the bread nuisances must be replaced \\emph{jointly} to be consistent. Row/cluster $<1$ confirms clustering is already included (and matters).}",
 "\\label{tab:ladder}\\begin{tabular}{rrr r rrrrr r}\\toprule",
 "Sc.&$n$&$\\tau$&$\\widehat f/f_0$&plug&$f$-true&$g$-true&oracle&fix-$J$&$\\tfrac{\\rm row}{\\rm clus}$\\\\\\midrule")
prev <- NULL
for (i in seq_len(nrow(d))) { x <- d[i,]; cur <- paste(x$scenario,x$n); L <- blk(L,prev,cur); prev <- cur
  L <- c(L, paste(paste(c(x$scenario,x$n,f2(x$tau),f2(x$Rf_med),f2(x$SER_plug),f2(x$SER_ftrue),
    f2(x$SER_gtrue),f2(x$SER_oracle),f2(x$SER_fixedJ),f2(x$ratio_row_over_cluster)),collapse=" & "),"\\\\")) }
L <- c(L, "\\bottomrule\\end{tabular}\\end{table}")
writeLines(L, "paper/table_ladder.tex")

# ---------- 6b. CV-selected df frequency ----------------------
{
  L <- c("\\begin{table}[t]\\centering\\small",
   "\\caption{Frequency of the CV-selected spline complexity $\\widehat J\\in\\{3,4,5\\}$ (the A-CV tuning), by Scenario, $n$, and $\\tau$ ($R=500$ per cell). Here $\\widehat J$ is the per-margin \\texttt{df}~$=$~degree~$+$~(interior knots), so $\\widehat J{=}3$ is a genuine \\emph{cubic} with \\emph{zero} interior knots ($4$ B-spline functions/margin, $16$ tensor functions), not a low-degree fit; the adaptive choice is the number of interior knots (defined in the reproducibility specification). \\textbf{The selection concentrates on the grid minimum $\\widehat J=3$} ($81$--$86\\%$), with $P(\\widehat J{=}4)\\approx0.10$--$0.15$ and $P(\\widehat J{=}5)\\approx0.02$--$0.04$, and is nearly invariant in $n$ and $\\tau$. Two consequences: (i) A-CV is here close to a \\emph{fixed} $J{=}3$ estimator, so its ``adaptive'' component is small and the A-Rich ($J{+}2$) contrast is effectively $J{=}5$ vs.\\ $J{=}3$ rather than a perturbation around an interior optimum; (ii) the selection saturates at the \\emph{minimum} complexity available within the cubic-spline family ($\\widehat J{=}3$ is zero interior knots, which a wider $J$ grid cannot undercut), so whether a lower-degree or parametric surface would be preferred is not determined by this experiment. This tempers the adaptive-tuning rung of Table~\\ref{tab:ladder}.}",
   "\\label{tab:cvfreq}\\begin{tabular}{rrr rrr r}\\toprule",
   "Sc.&$n$&$\\tau$&$P(\\widehat J{=}3)$&$P(\\widehat J{=}4)$&$P(\\widehat J{=}5)$&mode\\\\\\midrule")
  prev <- NULL
  for (i in seq_len(nrow(d))) { x <- d[i,]; cur <- paste(x$scenario,x$n); L <- blk(L,prev,cur); prev <- cur
    L <- c(L, paste(paste(c(x$scenario,x$n,f2(x$tau),f2(x$p_df3),f2(x$p_df4),f2(x$p_df5),x$df_mode),collapse=" & "),"\\\\")) }
  L <- c(L, "\\bottomrule\\end{tabular}\\end{table}"); writeLines(L, "paper/table_cvfreq.tex")
}

# ---------- 7. distributional robustness (non-normal errors) ---------------
if (file.exists("results/dist.csv")) {
  dd <- read.csv("results/dist.csv"); dd <- dd[order(dd$dist, dd$n, dd$tau), ]
  L <- c("\\begin{table}[t]\\centering\\small",
   "\\caption{Distributional robustness on the \\emph{heteroscedastic} Scenario~2 DGP ($\\sigma_1=0.4$), so the truth is $\\beta_{\\tau,1}^{\\rm true}=\\beta_{01}+\\sigma_1 q_\\tau(F_G)$ and shifts with the error law $F_G$ (reported as $\\beta_1^{\\rm tr}$ for independent verification). The marker error is a standardized non-normal law via a Gaussian copula (latent-Gaussian dependence preserved; raw within-subject correlation not held fixed): ``t5'' is standardized Student-$t_5$, ``skew'' a standardized right-skew lognormal. $R=500$, coverage over $300$ ($B=300$; MC SE $\\approx0.013$). Point estimation stays approximately unbiased for both slopes; bootstrap coverage is generally satisfactory ($0.91$--$0.98$), with moderate under-coverage (min $\\approx0.91$) in a few cells.}",
   "\\label{tab:dist}\\begin{tabular}{lrr r rr r rr}\\toprule",
   "err&$n$&$\\tau$&$\\beta_1^{\\rm tr}$&bias$_{\\beta1}$&ESD$_{\\beta1}$&cov$_{\\beta1}$&bias$_{\\beta2}$&cov$_{\\beta2}$\\\\\\midrule")
  prev <- NULL
  for (i in seq_len(nrow(dd))) { x <- dd[i,]; cur <- paste(x$dist,x$n); L <- blk(L,prev,cur); prev <- cur
    L <- c(L, paste(paste(c(x$dist,x$n,f2(x$tau),f3(x$true_b1),f3(x$bias_b1),f3(x$esd_b1),f3(x$cov_b1),
      f3(x$bias_b2),f3(x$cov_b2)),collapse=" & "),"\\\\")) }
  L <- c(L, "\\bottomrule\\end{tabular}\\end{table}"); writeLines(L, "paper/table_dist.tex")
}

# ---------- 8. dense quantile-crossing diagnostics -------------------------
if (file.exists("results/crossing.csv")) {
  cr <- read.csv("results/crossing.csv"); cr <- cr[order(cr$scenario, cr$n), ]
  L <- c("\\begin{table}[t]\\centering\\small",
   "\\caption{Quantile-crossing on a \\emph{dense} $30\\times30$ interior $(u,t)$ grid $\\times$ $10$ covariate profiles ($\\tau\\in\\{.1,.25,.5,.75,.9\\}$, adjacent pairs; $R=300$). Crossing rates split interior/boundary; ``$\\bar u/t$'' is the mean $u/t$ of crossing \\emph{incidents} (small $\\Rightarrow$ near the low-age boundary). The interior crossing rate is $0$ in every cell; boundary crossing rates round to $\\approx0$ by $n=1000$, but at $n{=}500$ a non-trivial \\emph{fraction of replicates} ($\\le15\\%$, Sc.~2) show \\emph{some} boundary crossing with per-incident magnitude up to $0.65$---so we rearrange by default. cal$_{\\rm in}$ is the in-sample $\\tau$-calibration error; cal$^{\\rm OOS}_{\\rm vis}$/cal$^{\\rm OOS}_{\\rm sub}$ are the \\emph{out-of-sample} errors on an \\emph{independent test draw} (same DGP), \\emph{visit}-weighted (mean over $\\approx620$ test rows) vs.\\ \\emph{subject}-equal-weighted (so multi-visit subjects do not dominate). The two agree ($0.011$--$0.022$); cal$^{\\rm OOS}_{\\rm vis}$ is the value \\emph{before} rearrangement, and rearrangement shifts it by at most $0.001$ in every cell (largest, Sc.~2 $n{=}500$: $0.0194\\!\\to\\!0.0184$; the full before/after pair is in the released \\texttt{crossing.csv}), so rearrangement does not distort out-of-sample calibration under either weighting.}",
   "\\label{tab:crossing}\\begin{tabular}{rr rr r r r r rr}\\toprule",
   "Sc.&$n$&cr$_{\\rm int}$&cr$_{\\rm bnd}$&frac.\\ reps&max mag&$\\bar u/t$&cal$_{\\rm in}$&cal$^{\\rm OOS}_{\\rm vis}$&cal$^{\\rm OOS}_{\\rm sub}$\\\\\\midrule")
  for (i in seq_len(nrow(cr))) { x <- cr[i,]
    L <- c(L, paste(paste(c(x$scenario,x$n,f4(x$cross_rate_interior),f4(x$cross_rate_boundary),
      f2(x$p_reps_any_cross),f3(x$max_mag),ifelse(is.na(x$incident_mean_u_over_t),"--",f2(x$incident_mean_u_over_t)),
      f4(x$cal_err_before),f4(x$cal_err_oos_before),f4(x$cal_err_oos_subj)),collapse=" & "),"\\\\")) }
  L <- c(L, "\\bottomrule\\end{tabular}\\end{table}"); writeLines(L, "paper/table_crossing.tex")
}

# ---------- 9. absolute one-step correction distribution -------------------
if (file.exists("results/corrdist.csv")) {
  co <- read.csv("results/corrdist.csv"); co <- co[order(co$n, co$tau), ]
  L <- c("\\begin{table}[t]\\centering\\small",
   "\\caption{B-OS-CF one-step correction $\\widetilde\\beta_{\\rm CF}-\\widehat\\beta_{\\rm ACV}$ and its \\emph{decomposition} (CC, Scenario~2, $R=500$). The signed mean is near zero (cancellation), but the per-replicate \\emph{absolute} correction is non-trivial (mean$|c|$ is $0.28$--$0.40$ of the A-CV ESD). Decomposing $\\mathrm{corr}=c_{\\rm split}+c_{\\rm os}$ into the sample-splitting part $c_{\\rm split}=\\sum_k\\tfrac{n_k}{n}\\widehat\\beta^{(-k)}_{\\rm init}-\\widehat\\beta_{\\rm ACV}$ and the pure one-step part $c_{\\rm os}=\\sum_k\\tfrac{n_k}{n}\\widehat S_k^{-1}\\widehat U_k$: the orthogonal part \\textbf{dominates} (mean$|c_{\\rm os}|/$ESD $\\approx0.28$--$0.40$ vs mean$|c_{\\rm split}|/$ESD $\\approx0.07$--$0.13$), so the adjustment is genuinely the orthogonal correction, not cross-fitting/splitting noise. Whether it reduces \\emph{bias} is separate, addressed in the tail regime of Table~\\ref{tab:sievetail}.}",
   "\\label{tab:corr}\\begin{tabular}{rr r rr rr}\\toprule",
   "$n$&$\\tau$&mean$|c|$&$\\tfrac{|c|}{\\rm ESD}$&$\\tfrac{|c_{\\rm os}|}{\\rm ESD}$&$\\tfrac{|c_{\\rm split}|}{\\rm ESD}$&signed mean\\\\\\midrule")
  prev <- NULL
  for (i in seq_len(nrow(co))) { x <- co[i,]; cur <- x$n; L <- blk(L,prev,cur); prev <- cur
    L <- c(L, paste(paste(c(x$n,f2(x$tau),f3(x$mean_abs),f2(x$mean_abs_over_esd),
      f2(x$os_over_esd),f2(x$split_over_esd),f4(x$mean_signed)),collapse=" & "),"\\\\")) }
  L <- c(L, "\\bottomrule\\end{tabular}\\end{table}"); writeLines(L, "paper/table_corr.tex")
}

# ---------- 10. sieve-tail remedy: does the one-step fix the tail bias? -----
if (file.exists("results/sievetail.csv")) {
  st <- read.csv("results/sievetail.csv"); st <- st[order(st$tau, st$n), ]
  L <- c("\\begin{table}[t]\\centering\\small",
   "\\caption{Does the one-step reduce the growing-sieve tail \\emph{bias}? Deterministic growing sieve, same replicates through A-joint-sieve, cross-fitted \\textbf{B-OS-CF}, and a richer-sieve A-Rich ($J{+}2$) sensitivity (CC, Sc.\\ 2; $R=400$, coverage over $150$, $B=120$---a \\emph{pilot-grade} inference subset, coverage MCSE $\\approx0.018$; $\\sqrt n\\,$bias MCSE $\\approx0.15$--$0.20$). At \\emph{both} tails B-OS-CF consistently moves the scaled bias toward zero, increasingly with $n$ (e.g.\\ $\\tau{=}0.1$: $0.66{\\to}0.28$, $0.53{\\to}0.19$; $\\tau{=}0.9$: $-0.49{\\to}-0.12$ at $n{=}4000$), whereas the fixed-$(J{+}2)$ richer sieve \\emph{worsens} it; paired Monte-Carlo intervals would be needed to resolve which individual moves are significant. \\textbf{However}, SER and RMSE show a \\emph{bias}-targeted correction only: B-OS-CF's SER $\\approx0.93$--$0.97$, coverage is not improved over the deterministic A-joint-sieve (\\emph{not} the CV-selected A-CV; the comparator here is the fixed growing sieve), and RMSE is essentially unchanged (the bias drop is offset by a small variance rise)---the scaled bias is already small relative to the SE here, so debiasing is a safeguard, not a coverage fix at these $n$.}",
   "\\label{tab:sievetail}\\begin{tabular}{rr r rrrr rrrr rr}\\toprule",
   "&&& \\multicolumn{4}{c}{A-joint-sieve} & \\multicolumn{4}{c}{B-OS-CF} & \\multicolumn{2}{c}{A-Rich}\\\\",
   "\\cmidrule(lr){4-7}\\cmidrule(lr){8-11}\\cmidrule(lr){12-13}",
   "$\\tau$&$n$&$J$&$\\sqrt n\\,$bias&SER&RMSE&cov&$\\sqrt n\\,$bias&SER&RMSE&cov&$\\sqrt n\\,$bias&cov\\\\\\midrule")
  prev <- NULL
  for (i in seq_len(nrow(st))) { x <- st[i,]; cur <- x$tau; L <- blk(L,prev,cur); prev <- cur
    L <- c(L, paste(paste(c(f2(x$tau),x$n,x$J,f2(x$js_sqrtn_bias),f2(x$js_ser),f3(x$js_rmse),f3(x$js_cov),
      f2(x$cf_sqrtn_bias),f2(x$cf_ser),f3(x$cf_rmse),f3(x$cf_cov),f2(x$rich_sqrtn_bias),f3(x$rich_cov)),collapse=" & "),"\\\\")) }
  L <- c(L, "\\bottomrule\\end{tabular}\\end{table}"); writeLines(L, "paper/table_sievetail.tex")
}

# ---------- 11. consolidated ESD / mean boot SE / SER / RMSE (supplement) ---
{
  rmse <- function(b, e) sqrt(b^2 + e^2)
  L <- c("\\begin{table}[t]\\centering\\small",
   "\\caption{Supplement to Table~\\ref{tab:estsim}: SER$=\\overline{\\mathrm{SE}}_{\\rm boot}/\\mathrm{ESD}$ and RMSE for \\emph{all four} estimators ($R=1000$, $R_{\\rm cov}=500$). Since bias is negligible in these cells, ESD$\\approx$RMSE and $\\overline{\\mathrm{SE}}=\\mathrm{SER}\\times\\mathrm{ESD}$; the underlying ESD and $\\overline{\\mathrm{SE}}$ are reported directly in Table~\\ref{tab:estfull}. \\emph{Most} coverage dips coincide with a SER slightly below $1$ (a mildly small bootstrap SE), \\emph{not} bias; the exception is the single lowest-coverage cell (Sc.~2, $n{=}500$, $\\tau{=}0.25$: cov $0.918$), which has SER $\\approx1.00$ \\emph{and} negligible bias ($0.006$), so its shortfall is neither mean SE nor bias. The studentization diagnostic (Table~\\ref{tab:studentized}) shows $Z=(\\widehat\\beta-\\beta_0)/\\widehat{\\mathrm{SE}}$ is near-Gaussian in its central moments (skew $-0.01$, excess kurtosis $0.02$) with normal $=$ percentile coverage, so the shortfall is not explained by mean SE, central skew/kurtosis, or a simple SE--error association; it is \\emph{consistent with} a mild small-sample upper-tail deviation and shrinks with $n$. B-OS-CF carries a modestly larger RMSE, generally around $5$--$15\\%$ (the variance cost of the explicit orthogonal adjustment), and---after the fold-leakage fix---a SER $\\ge1$, mildly to \\emph{moderately} conservative in several cells (coverage up to $\\approx0.986$).}",
   "\\label{tab:estse}\\begin{tabular}{rrr rr rr rr rr}\\toprule",
   "&&& \\multicolumn{2}{c}{A-CV cc} & \\multicolumn{2}{c}{A-CV ip} & \\multicolumn{2}{c}{B-OS-CF cc} & \\multicolumn{2}{c}{B-OS-CF ip}\\\\",
   "\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}\\cmidrule(lr){8-9}\\cmidrule(lr){10-11}",
   "Sc.&$n$&$\\tau$&SER&RMSE&SER&RMSE&SER&RMSE&SER&RMSE\\\\\\midrule")
  prev <- NULL
  for (i in seq_len(nrow(cc))) { x <- cc[i,]; cur <- paste(x$scenario,x$n); L <- blk(L,prev,cur); prev <- cur
    L <- c(L, paste(paste(c(x$scenario,x$n,f2(x$tau),
      f2(x$se_acv_cc/x$esd_acv_cc),f3(rmse(x$bias_acv_cc,x$esd_acv_cc)),
      f2(x$se_acv_ip/x$esd_acv_ip),f3(rmse(x$bias_acv_ip,x$esd_acv_ip)),
      f2(x$se_cf_cc/x$esd_cf_cc),f3(rmse(x$bias_cf_cc,x$esd_cf_cc)),
      f2(x$se_cf_ip/x$esd_cf_ip),f3(rmse(x$bias_cf_ip,x$esd_cf_ip))),collapse=" & "),"\\\\")) }
  L <- c(L, "\\bottomrule\\end{tabular}\\end{table}"); writeLines(L, "paper/table_estse.tex")
}

# ---------- 11b. full ESD + mean bootstrap SE supplement (round-11 #9) ------
{
  L <- c("\\begin{table}[t]\\centering\\small",
   "\\caption{Reporting-completeness supplement to Tables~\\ref{tab:estsim}/\\ref{tab:estse}: the Monte-Carlo ESD and mean bootstrap SE ($\\overline{\\mathrm{SE}}$) for all four estimators, reported \\emph{directly} so neither need be backed out of SER/RMSE ($R=1000$, $R_{\\rm cov}=500$). Within each pipeline B-OS-CF has a larger ESD than A-CV. The B-OS-CF$-$A-CV \\emph{difference} is dominated in \\emph{magnitude} by the orthogonal correction, with a smaller sample-splitting contribution (Table~\\ref{tab:corr}: $|c_{\\rm os}|/$ESD $\\approx0.28$--$0.40$ vs $|c_{\\rm split}|/$ESD $\\approx0.07$--$0.13$); we do \\emph{not} decompose the variance of the ESD increase itself (that would require the covariance terms with $\\widehat\\beta_{\\rm ACV}$). that $\\overline{\\mathrm{SE}}\\gtrsim$ESD for B-OS-CF is consistent with bootstrap-SE overestimation being an important contributor to its mild over-coverage. Bias and coverage are in Table~\\ref{tab:estsim}; SER and RMSE in Table~\\ref{tab:estse}.}",
   "\\label{tab:estfull}\\begin{tabular}{rrr rr rr rr rr}\\toprule",
   "&&&\\multicolumn{2}{c}{CC A-CV}&\\multicolumn{2}{c}{CC B-OS-CF}&\\multicolumn{2}{c}{IP A-CV}&\\multicolumn{2}{c}{IP B-OS-CF}\\\\",
   "\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}\\cmidrule(lr){8-9}\\cmidrule(lr){10-11}",
   "Sc.&$n$&$\\tau$&ESD&$\\overline{\\rm SE}$&ESD&$\\overline{\\rm SE}$&ESD&$\\overline{\\rm SE}$&ESD&$\\overline{\\rm SE}$\\\\\\midrule")
  prev <- NULL
  for (i in seq_len(nrow(cc))) { x <- cc[i,]; cur <- paste(x$scenario,x$n); L <- blk(L,prev,cur); prev <- cur
    L <- c(L, paste(paste(c(x$scenario,x$n,f2(x$tau),
      f3(x$esd_acv_cc),f3(x$se_acv_cc),f3(x$esd_cf_cc),f3(x$se_cf_cc),
      f3(x$esd_acv_ip),f3(x$se_acv_ip),f3(x$esd_cf_ip),f3(x$se_cf_ip)),collapse=" & "),"\\\\")) }
  L <- c(L, "\\bottomrule\\end{tabular}\\end{table}"); writeLines(L, "paper/table_estfull.tex")
}

# ---------- 12. tuning sensitivity (a_n x df_g) ----------------------------
if (file.exists("results/tuning.csv")) {
  tu <- read.csv("results/tuning.csv")
  cells <- unique(tu[, c("n", "tau")]); cells <- cells[order(cells$n, cells$tau), ]
  L <- c("\\begin{table}[t]\\centering\\small",
   "\\caption{Tuning sensitivity to the sparsity bandwidth $a_n\\in\\{.03,.05,.08\\}$ and projection dimension $df_g\\in\\{3,4,5\\}$ (defaults $.05,4$), now over the FULL grid $n\\in\\{500,1000,2000\\}\\times\\tau\\in\\{.1,.5,.9\\}$ (Sc.\\ 2, $R=400$). Per $(n,\\tau)$ we report the \\emph{range} across the $9$ combinations of: B-OS-CF $\\beta_{\\tau,1}$ bias and ESD; the plug-in sandwich SER and its coverage; the worst $\\widehat f$-cap rate; and the median bread condition number. \\textbf{Point estimation is robust to tuning everywhere} (bias range $\\le0.03$, ESD range $\\le0.03$, even at $n{=}500$ and the tails). \\textbf{Inference is not}: at $n{=}500$ the sandwich SER ranges $0.57$--$1.03$ (coverage $0.69$--$0.93$) across the grid---confirming that the $a_n$-dependent $\\widehat f$ destabilizes the analytic SE at small $n$ (the bootstrap is used for inference instead). This tuning-sensitivity of the analytic sandwich \\emph{persists at the tails even at} $n{=}2000$ (for $\\tau{=}0.1,0.9$, SER ranges $0.89$--$1.24$, coverage $0.91$--$0.98$ across the grid); only the central quantile is comparatively stable there. Note this varies the \\emph{sandwich}, whereas the recommended procedure uses the bootstrap---so it bounds the analytic SE's tuning-robustness, not the bootstrap's. $\\widehat f$-cap rates are $\\le2.2\\%$ and $\\widehat S$ is well-conditioned ($\\kappa\\le7$) throughout.}",
   "\\label{tab:tuning}\\begin{tabular}{rr rr rr rr}\\toprule",
   "&&\\multicolumn{2}{c}{bias / ESD range}&\\multicolumn{2}{c}{SER / cov range}&&\\\\",
   "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}",
   "$n$&$\\tau$&bias rng&ESD rng&SER rng&cov rng&$\\widehat f$-cap&$\\kappa_{\\max}$\\\\\\midrule")
  rng <- function(v) paste0("[", f3(min(v)), ",", f3(max(v)), "]")
  rng2 <- function(v) paste0("[", f2(min(v)), ",", f2(max(v)), "]")
  for (i in seq_len(nrow(cells))) {
    z <- tu[tu$n == cells$n[i] & tu$tau == cells$tau[i], ]
    L <- c(L, paste(paste(c(cells$n[i], f2(cells$tau[i]), rng(z$bias_cf), rng(z$esd_cf),
      rng2(z$ser_acv), rng2(z$cov_acv), f3(max(z$fcap)), sprintf("%.0f", max(z$cond_med))), collapse = " & "), "\\\\")) }
  L <- c(L, "\\bottomrule\\end{tabular}\\end{table}"); writeLines(L, "paper/table_tuning.tex")
}

# ---------- 13. censoring orthogonality: base + HARD regime ----
if (file.exists("results/compare.csv")) {
  cp <- read.csv("results/compare.csv")
  base <- cp[cp$regime == "base" & cp$tau == 0.5, ]; base <- base[order(base$scenario, base$n), ]
  L <- c("\\begin{table}[t]\\centering\\small",
   "\\caption{Identification/efficiency and the $\\chi^C=0$ diagnostic (A-CV $\\beta_{\\tau,1}$, $\\tau=0.5$, $R=400$; base-censoring regime). CC = unweighted complete-event; IPCW = estimated $\\widehat G_C$; ORC = oracle $G_{C,0}$. All three are unbiased; ESD$_{\\rm IPCW}\\approx$ESD$_{\\rm ORC}$ (no detectable additional variance from estimating $G_C$) and $\\ge$ESD$_{\\rm CC}$. The gap $\\mathrm{sd}(\\widehat\\beta_{\\rm IPCW}-\\widehat\\beta_{\\rm ORC})/$ESD$_{\\rm IPCW}$ shrinks with $n$ ($n^{-1}$-like, $\\ll n^{-1/2}$).}",
   "\\label{tab:compare}\\begin{tabular}{rr rrr rrr r}\\toprule",
   "Sc.&$n$&bias$_{\\rm CC}$&bias$_{\\rm IPCW}$&bias$_{\\rm ORC}$&ESD$_{\\rm CC}$&ESD$_{\\rm IPCW}$&ESD$_{\\rm ORC}$&gap/ESD\\\\\\midrule")
  prev <- NULL
  for (i in seq_len(nrow(base))) { x <- base[i,]; cur <- x$scenario; L <- blk(L,prev,cur); prev <- cur
    L <- c(L, paste(paste(c(x$scenario,x$n,f4(x$bias_cc),f4(x$bias_ipcw),f4(x$bias_orc),
      f4(x$esd_cc),f4(x$esd_ipcw),f4(x$esd_orc),f4(x$sd_gap/x$esd_ipcw)),collapse=" & "),"\\\\")) }
  L <- c(L, "\\bottomrule\\end{tabular}\\end{table}"); writeLines(L, "paper/table_compare.tex")

  # hard-vs-base weight distribution / positivity diagnostics (Sc 2, tau=0.5)
  h <- cp[cp$scenario == 2 & cp$tau == 0.5, ]; h <- h[order(h$regime, h$n), ]
  L <- c("\\begin{table}[t]\\centering\\small",
   "\\caption{Censoring-weight stress test: the $\\chi^C=0$ gap under the base regime vs.\\ a \\emph{hard} regime (later/wider entry $+$ heavier, more covariate-dependent loss; Sc.\\ 2, $\\tau=0.5$, $R=400$). The hard regime lifts censoring to $\\approx43\\%$ with markedly heavier IPCW-weight tails ($99$th pct $\\approx5$, max $\\approx12$) and lower weight ESS$/n_T$ ($\\approx0.79$ vs.\\ $0.98$), i.e.\\ closer to the positivity boundary; the Cox $\\widehat G_C$ never failed. Even so, the estimated-vs-oracle gap/ESD stays small relative to the estimator ESD and decreases overall ($0.20,0.19,0.13,0.13$ at $n=500,1000,2000,4000$); we do \\emph{not} claim a clean $n^{-1}$ rate---the final two sample sizes are nearly flat---only that the ratio stays small and decreases overall, which is \\emph{consistent with} (but does not establish) the second-order, $\\ll n^{-1/2}$ behaviour. ($\\tau=0.5$ only; tail behaviour under hard weights is not probed.)}",
   "\\label{tab:censhard}\\begin{tabular}{lr rr rrrr r}\\toprule",
   "regime&$n$&\\%cens&ESS/$n_T$&$\\widetilde w$&$w_{95}$&$w_{99}$&$w_{\\max}$&gap/ESD\\\\\\midrule")
  prev <- NULL
  for (i in seq_len(nrow(h))) { x <- h[i,]; cur <- x$regime; L <- blk(L,prev,cur); prev <- cur
    L <- c(L, paste(paste(c(x$regime,x$n,f2(x$p_cens),f2(x$ess_over_nT),f2(x$w_med),
      f2(x$w_p95),f2(x$w_p99),f2(x$w_max),f4(x$sd_gap/x$esd_ipcw)),collapse=" & "),"\\\\")) }
  L <- c(L, "\\bottomrule\\end{tabular}\\end{table}"); writeLines(L, "paper/table_censhard.tex")
}

# ---------- 14. computational reproducibility audit ------------
if (file.exists("results/audit.csv")) {
  au <- read.csv("results/audit.csv"); au <- au[order(au$scenario, au$n, au$tau), ]
  L <- c("\\begin{table}[t]\\centering\\small",
   "\\caption{Computational reproducibility audit, a \\emph{dedicated reduced-$B$ run} ($R=300$ per cell, inner $B{=}60$; \\emph{not} a summary of the headline $B{=}500$ log), $\\tau\\in\\{.25,.5,.75\\}$. Failure/fallback rates that would bias the headline numbers if silently dropped, for \\emph{both} the CC and IPCW pipelines: the QR (QR-solver) and CF (B-OS-CF) failure columns are shown \\emph{separately for both pipelines}; the remaining three columns---$\\widehat f$-cap, $\\kappa(\\widehat S)$, and the inner subject-bootstrap replicate-failure rate (``boot'')---are \\emph{CC-pipeline} diagnostics (the recommended estimator). \\emph{Failed replicates are counted, not discarded}. We also report the bread condition-number distribution (median, $95$th pct) rather than only a permissive (and loose) $\\kappa>10^8$ binary flag. \\emph{Denominators}: QR and CF are over the $R{=}300$ datasets; $\\widehat f$-cap over all fold--visit density estimates; boot over the $R{\\times}B$ inner resamples. \\textbf{The QR and B-OS-CF structural-fit failure rates are $0$ for both the CC and IPCW pipelines}, and the (CC) bootstrap-replicate failure rate is $0$ (we tabulate the CC bootstrap-failure column only, so we do \\emph{not} claim an IPCW bootstrap-failure audit here); a zero-event rate at $R{=}300$ has a $95\\%$ binomial upper bound of $\\approx1\\%$, not an exact $0$. $\\kappa(\\widehat S)$ stays small (median $\\approx4$, $95$th pct $\\le9.1$), so ill-conditioning is not a concern; the only non-zero mode is the $\\widehat f$-cap ($\\le1.0\\%$). \\emph{Scope}: base CC/IPCW pipelines at the central quartile quantiles only; the degenerate small-$n$ growing-sieve bootstrap (Table~\\ref{tab:sieve}) is \\emph{not} covered here and needs a separate distributional audit.}",
   "\\label{tab:audit}\\begin{tabular}{rrr rr rr rr rr}\\toprule",
   "&&&\\multicolumn{2}{c}{CC}&\\multicolumn{2}{c}{IPCW}&\\multicolumn{2}{c}{$\\kappa(\\widehat S)$}&&\\\\",
   "\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}\\cmidrule(lr){8-9}",
   "Sc.&$n$&$\\tau$&QR&CF&QR&CF&med&$q_{95}$&$\\widehat f$-cap&boot\\\\\\midrule")
  prev <- NULL
  for (i in seq_len(nrow(au))) { x <- au[i,]; cur <- paste(x$scenario,x$n); L <- blk(L,prev,cur); prev <- cur
    L <- c(L, paste(paste(c(x$scenario,x$n,f2(x$tau),f3(x$acv_fail_rate),f3(x$cf_fail_rate),
      f3(x$acv_ip_fail_rate),f3(x$cf_ip_fail_rate),f1(x$cond_med),f1(x$cond_q95),f3(x$fcap_rate),f3(x$boot_fail_rate)),collapse=" & "),"\\\\")) }
  L <- c(L, "\\bottomrule\\end{tabular}\\end{table}"); writeLines(L, "paper/table_audit.tex")
}

# ---------- 15. DW vs ORD: cellwise + interval geometry/score (r11/r13 #2/#3) --
if (file.exists("results/ordcenter2.csv")) {
  oc <- read.csv("results/ordcenter2.csv"); oc <- oc[oc$scenario == 2, ]; oc <- oc[order(oc$n, oc$tau), ]
  g  <- read.csv("results/ordcenter_global.csv")
  cap <- paste0(
   "Density-weighted (DW) vs.\\ ordinary (ORD) centering (Scenario~2; \\textbf{all three $\\tau$ fit on the same dataset per replicate}, $R=600$ per $n$; IPCW weight, \\emph{oracle} $f$, B-OS-\\emph{NCF}; only the centering $g$ differs---Sc.~1 has $g_{\\rm DW}{=}g_{\\rm ORD}$ and is omitted). DW covers better in \\emph{every} cell ($\\Delta_{\\rm cov}=+0.01$ to $+0.03$; all discordances favour DW, $c{=}0$): the paired exact McNemar is significant in $\\mathbf{8/9}$ cells unadjusted, $\\mathbf{4/9}$ after Holm, $\\mathbf{3/9}$ after Bonferroni (we report only these cellwise/multiplicity-adjusted results; a single global $p$ is omitted as hard to audit and unnecessary). \\textbf{The effect is bread calibration, not a different estimator}: DW and ORD point estimates coincide to $\\sim10^{-4}$, the DW interval is only $\\approx6\\%$ longer (length/SE ratio $\\approx1.06$) and \\emph{nests} the ORD interval in $\\ge99\\%$ of (rep,$\\tau$). Its proper $95\\%$ interval score (IS; lower$=$better; width $+$ miss penalty) is \\emph{statistically indistinguishable} from ORD's: the per-replicate-clustered mean $\\overline{\\Delta\\mathrm{IS}}=",
   f4(g$dis_clustered_mean), "$ ($95\\%$ CI $[", f4(g$dis_ci_lo), ",", f4(g$dis_ci_hi),
   "]$, includes $0$), and cellwise $\\Delta\\mathrm{IS}={}$IS$_{\\rm DW}-{}$IS$_{\\rm ORD}$ is within $\\approx1.6$ MCSE of $0$ \\emph{except} the smallest cell ($n{=}500,\\tau{=}0.25$: $+0.011$, $2.4$ MCSE, where the wider DW interval is slightly \\emph{worse}-scored). \\textbf{So DW improves coverage through a modestly larger sandwich SE, with \\emph{no detectable change} in the proper interval score} (the cellwise $\\Delta$IS are small and not uniformly favourable). The clustered $\\overline{\\Delta\\mathrm{IS}}$ averages $\\Delta$IS over the three $\\tau$ \\emph{within} each dataset and treats the independent datasets (across $n$) as the clusters; the resulting CI targets the \\emph{equally-weighted} average $\\Delta$IS over the nine evaluated $(n,\\tau)$ design cells (a simulation-design summary, not a single-distribution parameter). This is an oracle-$f$, non-cross-fitted diagnostic and does not by itself carry to the feasible B-OS-CF. Coverage $=$ plug-in sandwich.")
  L <- c("\\begin{table}[t]\\centering\\small", paste0("\\caption{", cap, "}"),
   "\\label{tab:ordcenter}\\begin{tabular}{rr rr rr rr r}\\toprule",
   "&&\\multicolumn{2}{c}{coverage}&\\multicolumn{2}{c}{int.\\ length}&\\multicolumn{2}{c}{int.\\ score}&\\\\",
   "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\\cmidrule(lr){7-8}",
   "$n$&$\\tau$&DW&ORD&DW&ORD&DW&ORD&$\\Delta$IS\\\\\\midrule")
  prev <- NULL
  for (i in seq_len(nrow(oc))) { x <- oc[i,]; cur <- x$n; L <- blk(L,prev,cur); prev <- cur
    L <- c(L, paste(paste(c(x$n,f2(x$tau),f3(x$cov_dw),f3(x$cov_ord),
      f3(x$ilen_dw),f3(x$ilen_ord),f3(x$is_dw),f3(x$is_ord),f4(x$dis_mean)),collapse=" & "),"\\\\")) }
  L <- c(L, "\\bottomrule\\end{tabular}\\end{table}"); writeLines(L, "paper/table_ordcenter.tex")
}

# ---------- 16. studentization diagnostic (round-10 #3) --------------------
if (file.exists("results/studentized.csv")) {
  su <- read.csv("results/studentized.csv"); su <- su[su$scenario == 2, ]
  su <- su[order(su$n, su$tau), ]
  L <- c("\\begin{table}[t]\\centering\\small",
   "\\caption{Studentization diagnostic for the headline CC--A-CV estimator (Scenario~2, where coverage dips concentrate; $R=400$, $B=250$ per replicate). For the standardized statistic $Z_r=(\\widehat\\beta_r-\\beta_0)/\\widehat{\\mathrm{SE}}_r$ we report its skewness, excess kurtosis, and $2.5/97.5\\%$ quantiles; the normal-CI miss rates by side (tr$<$: truth below the CI, $Z>1.96$; tr$>$: truth above, $Z<-1.96$); the cross-replicate correlation $\\rho(|\\widehat\\beta-\\beta_0|,\\widehat{\\mathrm{SE}})$ (large negative $\\Rightarrow$ big errors paired with small SE, which would hide undercoverage); and the coverage of the \\emph{normal}-SE, \\emph{percentile}, and \\emph{basic} bootstrap CIs. \\textbf{At the lowest-coverage cell ($n{=}500$, $\\tau{=}0.25$) $Z$ has near-Gaussian central moments} (skew $-0.01$, excess kurtosis $0.02$) with $\\rho\\approx0$, and normal ($0.93$) and percentile ($0.93$) CIs agree---so the dip is \\emph{not} explained by skew/kurtosis, SE--error dependence, or the choice of CI form. The results \\emph{suggest} a finite-sample upper-tail deviation ($z_{.975}{=}2.34>1.96$); the side-specific miss-rate difference (tr$<$ $0.043$ vs tr$>$ $0.028$, i.e.\\ $0.015\\pm0.013$) is imprecise and not clearly distinguishable from $0$, and the deviation shrinks with $n$. Scenario~1 (omitted) is cleaner. The diagnostic gives \\emph{no} evidence that percentile intervals improve on the normal-SE interval (the basic interval is worse); both show a similar mild finite-sample shortfall, \\emph{consistent with} a finite-sample upper-tail deviation of the studentized statistic. The $z_{.975}{=}2.34$ rests on only $\\sim10$ tail replicates ($R{=}400$), so read it qualitatively; this is a supporting re-run (cov $\\approx0.93$), not a per-replicate decomposition of the headline cell.}",
   "\\label{tab:studentized}\\begin{tabular}{rr rr r@{,}l rr r rrr}\\toprule",
   "&&\\multicolumn{2}{c}{$Z$ shape}&\\multicolumn{2}{c}{$Z$ quant.}&\\multicolumn{2}{c}{CI miss}&&\\multicolumn{3}{c}{coverage}\\\\",
   "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\\cmidrule(lr){7-8}\\cmidrule(lr){10-12}",
   "$n$&$\\tau$&skew&ex.kurt&\\multicolumn{2}{c}{$z_{.025/.975}$}&tr$<$&tr$>$&$\\rho_{e,\\rm SE}$&norm&pct&basic\\\\\\midrule")
  prev <- NULL
  for (i in seq_len(nrow(su))) { x <- su[i,]; cur <- x$n; L <- blk(L,prev,cur); prev <- cur
    L <- c(L, paste(paste(c(x$n,f2(x$tau),f2(x$z_skew),f2(x$z_exkurt),f2(x$z_q025),f2(x$z_q975),
      f3(x$miss_lo),f3(x$miss_hi),f2(x$corr_err_se),f3(x$cov_normal),f3(x$cov_pct),f3(x$cov_basic)),collapse=" & "),"\\\\")) }
  L <- c(L, "\\bottomrule\\end{tabular}\\end{table}"); writeLines(L, "paper/table_studentized.tex")
}

# ---------- 17. fold-randomization sensitivity (round-10 #4/#7) -------------
if (file.exists("results/foldsens.csv")) {
  fs <- read.csv("results/foldsens.csv"); fs <- fs[order(fs$n, fs$tau), ]
  L <- c("\\begin{table}[t]\\centering\\small",
   "\\caption{Fold-randomization sensitivity of the CC \\textbf{B-OS-CF} bootstrap (Scenario~2; $R=250$, $B=120$, $M=8$ split draws). The point estimator is the single \\emph{realized} fold split (as in the headline) or the \\emph{average} over $M$ splits; the bootstrap SE either redraws folds each resample (headline default) or keeps each original subject's realized fold (\\emph{fixed}). SER is relative to the matching point estimator's ESD. \\textbf{Redrawing folds inflates the SER by only $\\approx0.01$--$0.015$ over fixing them} ($\\Delta$SER $=$ SER$_{\\rm redr}-$SER$_{\\rm fix}$; paired MCSE in parentheses $=$ the sd of the per-replicate $\\widehat{\\mathrm{SE}}_{\\rm redr}-\\widehat{\\mathrm{SE}}_{\\rm fix}$ divided by the common ESD and $\\sqrt R$, treating ESD as fixed); what matters is the \\emph{effect size}---fold redrawing adds only $\\approx0.01$--$0.015$ to the SER, a small contribution---so the bulk of the mild B-OS-CF conservativeness is the duplicate-ID fold-leakage correction (Table~\\ref{tab:leakfold}), not fold randomization. The headline bootstrap \\emph{redraws} folds each resample and so marginalizes over fold randomness---this is the recommended inference; the fixed-fold and repeated-split rows are \\emph{sensitivity checks}. Fixed folds are marginally tighter and empirically well calibrated here (coverage $0.948$--$0.968$); averaging the point over splits lowers its ESD but, with the redrawn SE, \\emph{over}-covers ($\\le0.98$). Per-cell coverage MCSE $\\approx0.014$ ($R{=}250$), so cell-level coverage gaps should not be over-read.}",
   "\\label{tab:foldsens}\\begin{tabular}{rr rr rr r rrr}\\toprule",
   "&&\\multicolumn{2}{c}{ESD}&\\multicolumn{2}{c}{SER}&&\\multicolumn{3}{c}{coverage}\\\\",
   "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\\cmidrule(lr){8-10}",
   "$n$&$\\tau$&real&avg&fix&redr&$\\Delta$SER&$\\substack{\\rm real\\\\ \\rm fix}$&$\\substack{\\rm real\\\\ \\rm redr}$&$\\substack{\\rm avg\\\\ \\rm redr}$\\\\\\midrule")
  prev <- NULL
  for (i in seq_len(nrow(fs))) { x <- fs[i,]; cur <- x$n; L <- blk(L,prev,cur); prev <- cur
    L <- c(L, paste(paste(c(x$n,f2(x$tau),f3(x$esd_real),f3(x$esd_avg),f2(x$ser_fixed),f2(x$ser_redrawn),
      paste0(f3(x$ser_diff),"(",f3(x$ser_diff_mcse),")"),
      f3(x$cov_real_fixed),f3(x$cov_real_redrawn),f3(x$cov_avg_redrawn)),collapse=" & "),"\\\\")) }
  L <- c(L, "\\bottomrule\\end{tabular}\\end{table}"); writeLines(L, "paper/table_foldsens.tex")
}

# ---------- 18. growing-sieve bootstrap-SE audit + CAUSAL diagnosis (r10/r11) -
if (file.exists("results/sievedist2.csv")) {
  sd0 <- read.csv("results/sievedist2.csv"); sd0 <- sd0[order(sd0$mode, sd0$n), ]
  L <- c("\\begin{table}[t]\\centering\\small",
   "\\caption{Distributional audit of the deterministic-sieve subject bootstrap at $\\tau{=}0.5$ (same DGP/seed as Table~\\ref{tab:sieve}; $R=250$, $B=150$), with a \\emph{scale-invariant} numerical-stability diagnostic. \\textbf{The median bootstrap SE equals the exact sieve's at every $n$}; only the $n{=}500$ growing sieve blows up---SE$_{\\max}=96$, and the SE$>3\\times$med fraction is $0.8\\%$ ($2/250$ datasets). To test \\emph{why}, we replace the raw condition number (which is huge even where the SE is stable) with scale-invariant degeneracy of the standardized design Gram of $[X_1,X_2,\\text{surface basis}]$: the \\emph{effective rank} $r_{\\rm eff}=(\\sum_j\\lambda_j)^2/\\sum_j\\lambda_j^2$ (a participation ratio over the Gram eigenvalues, \\emph{not} a numerical rank). Near-zero-support surface columns are removed by a deterministic pre-fit screen; \\emph{all} columns surviving that screen are kept in the QR fit (no further numerical-rank truncation, and QR never fails). $N_{\\rm low}=$ the count of surface columns carrying $<3$ rows of mass. The $r_{\\rm eff}$ denominator is the design width $=p{+}K_{\\rm ret}$, the $p{=}2$ covariates plus the retained surface tensor columns (near-zero columns dropped): $27{=}2{+}25$ at $J{=}4$, $63{=}2{+}61$ of $64$ at $J{=}7$. The basis \\emph{is} over-rich at $n{=}500$ ($r_{\\rm eff}\\approx21$ of $63$ columns, $\\approx8$ low-support columns)---\\textbf{but this does \\emph{not} explain the blow-up}: $\\mathrm{corr}(\\log\\widehat{\\mathrm{SE}},-\\log\\lambda_{\\min})\\approx0$, and the $2$ extreme-SE datasets are \\emph{not} among the low-$\\lambda_{\\min}$ (weakest-Gram) ones (coincidence $0/2$). We also audited the blow-up at the bootstrap-\\emph{resample} level: the extreme-$|\\widehat\\beta^*|$ resamples (top $0.5\\%$) are \\emph{indistinguishable} from the bulk in lost-support columns ($7.7$ vs $7.6$), standardized-Gram $\\lambda_{\\min}$ ($1.5$ vs $1.6{\\times}10^{-3}$), and onset diversity ($191$ vs $195$), with all correlations $\\approx0$. So the instability is \\emph{not} a Gram-rank or support failure at either the dataset or the resample level; it appears to be a rare numerical/optimization instability of the over-rich small-$n$ sieve rather than a failure of the measured design-rank diagnostics, and we leave the precise trigger \\emph{unidentified} (it does not affect the recommended fixed sieve). We therefore present it descriptively (a small-$n$ over-rich-sieve artifact, $2/250$ datasets) and \\emph{withdraw} any condition-number causal attribution. The recommended fixed sieve showed no analogous instability in the audited runs. (Column big/$\\cap$rd: the number of datasets with bootstrap SE $>3\\times$ the median, and how many of those also fall in the weakest-Gram decile ($\\lambda_{\\min}$ in the bottom $10\\%$; ``weakest-Gram'', not necessarily rank-deficient).)}",
   "\\label{tab:sievedist}\\begin{tabular}{l r r r r r r r r}\\toprule",
   "&&&\\multicolumn{2}{c}{bootstrap SE}&\\multicolumn{2}{c}{degeneracy}&&\\\\",
   "\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}",
   "sieve&$n$&$J$&med&max&$r_{\\rm eff}$&$N_{\\rm low}$&$\\rho_{\\log{\\rm SE},-\\log\\lambda}$&big/$\\cap$rd\\\\\\midrule")
  prev <- NULL
  for (i in seq_len(nrow(sd0))) { x <- sd0[i,]; cur <- x$mode; L <- blk(L,prev,cur); prev <- cur
    L <- c(L, paste(paste(c(x$mode,x$n,x$J,f3(x$se_med),f2(x$se_max),
      sprintf("%.0f/%.0f",x$eff_rank_med,x$ncol), sprintf("%.0f",x$n_unsupp_med),
      f2(x$corr_logSE_neglogLam), sprintf("%.0f/%.0f",x$n_bigSE,x$n_bigSE_and_rankdef)),collapse=" & "),"\\\\")) }
  L <- c(L, "\\bottomrule\\end{tabular}\\end{table}"); writeLines(L, "paper/table_sievedist.tex")
}

# ---------- 19. fixed-J=3 vs A-CV comparator + CV-loss gaps (round-11 #5) ----
if (file.exists("results/fixedj.csv")) {
  fj <- read.csv("results/fixedj.csv"); fj <- fj[order(fj$n, fj$tau), ]
  L <- c("\\begin{table}[t]\\centering\\small",
   "\\caption{Is the data-adaptive A-CV materially different from a \\emph{fixed} small sieve? (Scenario~2, CC estimator, same replicates; $R=400$, $B=200$.) Since CV selects $\\widehat J{=}3$ in $80$--$86\\%$ of replicates (Table~\\ref{tab:cvfreq}), we put A-CV (CV-selected $\\widehat J$) head-to-head with fixed $J{=}3$. \\textbf{CV genuinely prefers $J{=}3$---it is not tie-breaking noise}: the normalized CV-loss gaps (gap$_J=(\\mathrm{CV}(J){-}\\mathrm{CV}(3))/|\\mathrm{CV}(3)|$, the relative increase in held-out \\emph{complete-event} (CC, $\\omega{=}1$) check loss---the CC member tunes its df by CC-CV, so \\emph{no} $\\widehat G_C$ enters the CC pipeline) are positive (gap$_4\\approx0.004$--$0.029$, $P(\\text{gap}_4{>}0)\\approx0.85$), shrinking with $n$. \\textbf{And A-CV $\\approx$ fixed $J{=}3$}: the point estimates differ by only $\\overline{|\\Delta|}\\approx0.002$--$0.007$, with near-identical bias, ESD, and bootstrap coverage. On the \\emph{surface} IMSE, however, fixed $J{=}3$ is \\emph{consistently lower}---by $\\approx9$--$13\\%$ in the displayed cells; the paired $\\overline{\\Delta\\mathrm{IMSE}}={}$IMSE$_{\\rm ACV}-{}$IMSE$_{J3}$ (last column, MCSE in parentheses) is positive throughout, consistent with occasional overfitting when A-CV selects $J{\\in}\\{4,5\\}$. So in these DGPs \\textbf{CV mostly reproduces the fixed $J{=}3$ fit and offers little measurable gain on $\\beta$}---the headline simulation therefore essentially validates a small fixed sieve. We do not claim CV equals the per-replicate oracle complexity, nor that it ``protects'' against richer fits (the average CV criterion and the majority of replicates favour the smallest candidate, but $\\approx14$--$20\\%$ select $J{\\in}\\{4,5\\}$); only that its preference for $J{=}3$ is systematic, not tie-breaking noise.}",
   "\\label{tab:fixedj}\\begin{tabular}{rr r rr r rr rr r}\\toprule",
   "&&&\\multicolumn{2}{c}{CV gap}&&\\multicolumn{2}{c}{coverage}&\\multicolumn{2}{c}{IMSE}&\\\\",
   "\\cmidrule(lr){4-5}\\cmidrule(lr){7-8}\\cmidrule(lr){9-10}",
   "$n$&$\\tau$&$P(\\widehat J{=}3)$&$\\substack{\\rm CV4\\\\-\\rm CV3}$&$\\substack{\\rm CV5\\\\-\\rm CV3}$&$\\overline{|\\Delta\\beta|}$&A-CV&$J{=}3$&A-CV&$J{=}3$&$\\Delta$IMSE\\\\\\midrule")
  prev <- NULL
  for (i in seq_len(nrow(fj))) { x <- fj[i,]; cur <- x$n; L <- blk(L,prev,cur); prev <- cur
    L <- c(L, paste(paste(c(x$n,f2(x$tau),f2(x$p_df3),f3(x$gap4_med),f3(x$gap5_med),
      f4(x$mean_abs_diff),f3(x$cov_acv),f3(x$cov_j3),f4(x$imse_acv),f4(x$imse_j3),
      paste0(f4(x$dimse_mean),"(",f4(x$dimse_mcse),")")),collapse=" & "),"\\\\")) }
  L <- c(L, "\\bottomrule\\end{tabular}\\end{table}"); writeLines(L, "paper/table_fixedj.tex")
}

# ---------- 20. leaked-vs-clean fold bootstrap: CAUSAL test (round-12 #2) ----
if (file.exists("results/leakfold.csv")) {
  lf <- read.csv("results/leakfold.csv"); lf <- lf[order(lf$n, lf$tau), ]
  L <- c("\\begin{table}[t]\\centering\\small",
   "\\caption{\\textbf{Controlled implementation comparison} of the duplicate-ID fold-leakage correction (CC B-OS-CF, Scenario~2; $R=250$, $B=150$). On the \\emph{same} datasets and resamples the bootstrap SE is computed two ways: \\emph{leaked} (the previous implementation---fold by the resampled id, so duplicate copies of a subject split across folds and one trains a nuisance the other is scored against) vs.\\ \\emph{clean} (fold by orig\\_id, duplicates kept together). The last column is the \\emph{paired} ratio $\\widehat{\\mathrm{SE}}_{\\rm clean}/\\widehat{\\mathrm{SE}}_{\\rm leaked}$ computed per replicate then averaged (MCSE in parentheses $=$ sd of the per-replicate ratio $/\\sqrt R$). \\textbf{The leak systematically tightens the bootstrap SE by $6$--$11\\%$}: SER$_{\\rm leaked}<1$ in nearly every cell with coverage dipping to $\\approx0.93$, whereas enforcing fold discipline restores SER$_{\\rm clean}\\ge1$ and near-nominal coverage ($0.94$--$0.98$). This attributes most of the mild B-OS-CF conservativeness of Table~\\ref{tab:estsim} to the correction; combined with the $\\le0.015$ fold-\\emph{randomization} contribution of Table~\\ref{tab:foldsens}, it explains the bulk of the shift from underestimation to mild conservativeness, though residual overcoverage in several cells (clean SER up to $1.17$) is \\emph{not} uniquely attributable to it. Per-cell coverage MCSE $\\approx0.014$ ($R{=}250$), so cell-level coverage gaps ($\\approx0.016$--$0.020$) should not be over-read individually.}",
   "\\label{tab:leakfold}\\begin{tabular}{rr rr r rr}\\toprule",
   "&&\\multicolumn{2}{c}{SER}&&\\multicolumn{2}{c}{coverage}\\\\",
   "\\cmidrule(lr){3-4}\\cmidrule(lr){6-7}",
   "$n$&$\\tau$&leaked&clean&$\\tfrac{\\rm clean}{\\rm leaked}$&leaked&clean\\\\\\midrule")
  prev <- NULL
  for (i in seq_len(nrow(lf))) { x <- lf[i,]; cur <- x$n; L <- blk(L,prev,cur); prev <- cur
    L <- c(L, paste(paste(c(x$n,f2(x$tau),f3(x$ser_leaked),f3(x$ser_clean),
      paste0(f3(x$se_ratio),"(",f3(x$se_ratio_mcse),")"),
      f3(x$cov_leaked),f3(x$cov_clean)),collapse=" & "),"\\\\")) }
  L <- c(L, "\\bottomrule\\end{tabular}\\end{table}"); writeLines(L, "paper/table_leakfold.tex")
}

cat("wrote tables incl. surface, crossing, corr, sievetail, estse, tuning, compare, censhard, audit, ordcenter, studentized, foldsens, sievedist, fixedj, leakfold\n")
