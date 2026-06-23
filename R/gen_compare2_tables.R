# gen_compare2_tables.R — Sun-vs-quantile comparison tables (5 DGPs).
f3 <- function(x) { x <- ifelse(!is.na(x) & abs(x) < 5e-4, 0, x); ifelse(is.na(x), "--", formatC(x, format = "f", digits = 3)) }
f2 <- function(x) { x <- ifelse(!is.na(x) & abs(x) < 5e-3, 0, x); ifelse(is.na(x), "--", formatC(x, format = "f", digits = 2)) }
cells <- list.files("results", pattern = "^compare2_cell_[0-9]+\\.csv$", full.names = TRUE)
d <- do.call(rbind, lapply(cells, read.csv)); d <- d[order(d$dgp, d$n), ]
dgps <- sort(unique(d$dgp))
dn <- c("1"="\\textbf{DGP1} proportional, homoscedastic mean",
        "2"="\\textbf{DGP2} proportional mean, dispersion varies",
        "3"="\\textbf{DGP3} non-proportional mean (headline)",
        "4"="\\textbf{DGP4} tail-specific effect",
        "5"="\\textbf{DGP5} heavy-tailed error (mean undefined)")
# generic block: label spanning row, then one row per n
block <- function(L, dg, ncol, body) {
  rows <- d[d$dgp == dg, ]
  L <- c(L, sprintf("\\multicolumn{%d}{l}{%s}\\\\", ncol, dn[as.character(dg)]))
  for (i in seq_len(nrow(rows))) L <- c(L, body(rows[i, ]))
  c(L, "\\addlinespace")
}
emit <- function(L, ncol, body) { for (dg in dgps) L <- block(L, dg, ncol, body)
  c(L[-length(L)], "\\bottomrule\\end{tabular}\\end{table}") }

# ---- Table 1: beta_tau recovery + efficiency (Sun ESD vs our ESD) ------------
L <- c("\\begin{table}[t]\\centering\\small",
 "\\caption{\\textbf{Distributional effect recovery and scalar-slope variability (Sun's real estimator vs.\\ ours).} Same Monte-Carlo datasets under Sun et al.'s survival/truncation/fixed-visit design ($R=500$; common interior domain). Sun's profile partial-likelihood slope $\\widehat\\beta_1$ (their compiled \\texttt{marker3.cpp}) and its Monte-Carlo ESD vs.\\ the ESD of our median slope $\\widehat\\beta_{.5,1}$, plus our $\\beta_{\\tau,1}$ bias at $\\tau=0.1,0.5,0.9$ and the recovered tail spread $\\widehat\\beta_{.9}-\\widehat\\beta_{.1}$. \\emph{Both estimators are unbiased when their target is well-defined.} In DGP1 the direct median slope has smaller Monte-Carlo variability than the feasible Sun slope for the shared scalar target ($\\sim$$4\\times$ smaller ESD); this is a point-estimation comparison, not an inferential-efficiency claim. An oracle-$\\alpha$ decomposition shows Sun's feasible $\\widehat\\beta=\\widehat\\theta-\\widehat\\alpha$ (Thm.~3.1) is \\emph{not} less stable because it estimates $\\alpha$: the oracle-$\\alpha$ version $\\widehat\\theta-\\alpha_0$ is in fact \\emph{more} variable, so the dominant variability appears to come from the profile-kernel marker score $\\widehat\\theta$, with $\\widehat\\theta-\\widehat\\alpha$ partly benefiting from covariance cancellation. \\emph{The Sun-ESD vs.\\ Q-ESD comparison is like-for-like only in DGP1, where the mean and median slopes coincide ($\\beta_{\\rm mean}=\\beta_{.5}$); in DGP2--DGP4 the scalar targets differ, and in DGP5 the conditional mean is undefined, so Sun's slope variability is reported only as a fragility diagnostic, not an efficiency comparison.} DGP2 only the quantile slope sees the covariate-driven dispersion; DGP3 no proportional slope exists; DGP4 the effect is confined to the upper tail ($\\beta_{.9}\\!\\gg\\!\\beta_{.5}$) so Sun's scalar hides it; in DGP5 (heavy $t_5$ log-tail, mean ill-posed) Sun's raw ESD is unstable and non-monotone, driven by rare explosive observations, while robust spread summaries shrink with $n$ (Table~\\ref{tab:cmp-heavy}) and our median slope is robust.}",
 "\\label{tab:cmp-beta}\\begin{tabular}{r rr rrr r}\\toprule",
 "&&&\\multicolumn{3}{c}{Q $\\beta_{\\tau,1}$ bias}&spread $\\beta_{.9}{-}\\beta_{.1}$\\\\",
 "\\cmidrule(lr){4-6}",
 "$n$&Sun $\\widehat\\beta_1$ (ESD)&Q ESD$_{.5}$&$.10$&$.50$&$.90$&Q / true\\\\\\midrule")
L <- emit(L, 7, function(x) paste0(paste(c(x$n,
  paste0(f2(x$beta_sun1),"\\,(",f2(x$esd_sun1),")"), f3(x$esd_bq50),
  f3(x$bias_bq1), f3(x$bias_bq5), f3(x$bias_bq9),
  paste0(f2(x$spread90_q),"/",f2(x$spread90_true))), collapse=" & "), "\\\\"))

# ---- Table 2: mean-contrast IMSE (headline) -- both reconstructions, x10^3 ----
k3 <- function(x) f3(1000 * x)   # report IMSE x10^3 to avoid spurious 0.000
L <- c(L, "", "\\begin{table}[t]\\centering\\small",
 "\\caption{\\textbf{Conditional-mean contrast --- the headline} ($\\mathrm{IMSE}\\times10^{3}$). Both methods target the conditional mean. Sun: the \\emph{basic} proportional model ($E[Y|u,t,x]=\\mu_0(u,t)e^{\\beta'x}$, log-contrast $=\\widehat\\beta_1$, constant) and his Sec.~3.5 \\emph{stratified-baseline} extension (Sun-S: his genuine estimator fit separately per $X_1$ group, giving a nonparametric $(u,t)$-varying contrast). Ours: the model-free quantile integral Q-INT ($\\int_0^1 e^{Q_\\tau}d\\tau$) and the tail-stable retransformation Q-LN ($\\widehat Q_{.5}+\\tfrac12\\widehat\\sigma^2$). IMSE of the mean log-contrast $\\Delta(u,t)$. Under non-proportionality (DGP3) $\\Delta$ truly varies: Sun-basic collapses it to a pseudo-average, and the fully stratified Sun-S recovers it only at a large \\emph{variance} cost (two nonparametric surfaces), so its IMSE exceeds Sun-basic and far exceeds Q (both Q reconstructions several-fold lower; the paired $\\text{Q-LN}-\\text{Sun-basic}$ difference is significantly negative, MCSE in parentheses). Sun-S is a fully stratified implementation and is \\emph{not} intended to represent all prespecified low-dimensional interaction bases. In DGP1 the lower Q IMSE reflects the smaller Monte-Carlo variability of the direct shared scalar slope (Table~\\ref{tab:cmp-beta}); in DGP2 Q-LN benefits from the log-location-scale structure, whereas the model-free Q-INT is more tail-sensitive and is \\emph{not} uniformly better than Sun ($\\text{Q-INT}-\\text{Sun}=+2.5\\times10^{-3}$ at $n{=}2000$). In \\textbf{DGP4 Sun-basic is \\emph{best}} (well matched to the proportional mean target). DGP5's mean is undefined ($-$).}",
 "\\label{tab:cmp-contrast}\\begin{tabular}{r rrrr rr}\\toprule",
 "&\\multicolumn{4}{c}{contrast IMSE $\\times10^{3}$}&\\multicolumn{2}{c}{paired $-$Sun $\\times10^{3}$ (MCSE)}\\\\",
 "\\cmidrule(lr){2-5}\\cmidrule(lr){6-7}",
 "$n$&Sun&Sun-S&Q-LN&Q-INT&Q-LN&Q-INT\\\\\\midrule")
L <- emit(L, 7, function(x) paste0(paste(c(x$n, k3(x$ctr_sun), k3(x$ctr_sstrat), k3(x$ctr_q), k3(x$ctr_qint),
  paste0(k3(x$dctr_mean),"\\,(",k3(x$dctr_mcse),")"),
  paste0(k3(x$dctr_int_mean),"\\,(",k3(x$dctr_int_mcse),")")), collapse=" & "), "\\\\"))

# ---- Table 3: mean-surface log-IMSE + Sun bandwidth sensitivity ($x10^3$) -----
L <- c(L, "", "\\begin{table}[t]\\centering\\small",
 "\\caption{\\textbf{Conditional-mean surface and Sun bandwidth sensitivity} (log-IMSE $\\times10^{3}$, $x_1{=}0$). Sun's kernel baseline $\\widehat\\mu_0$ at his density-reference bandwidth $h$ ($2.34\\,\\mathrm{sd}(A)N^{-1/6}$) and at the regression bandwidth $0.5h$, vs.\\ our Q-LN and Q-INT reconstructions. Sun's surface is bandwidth-sensitive---the density rule over-smooths the steep baseline ($h$ column), and $0.5h$ roughly halves the IMSE---so the Sun-vs-Q surface gap is partly a smoothing-bandwidth effect, not a method deficiency; the bandwidth-free slope and contrast (Tables~\\ref{tab:cmp-beta}--\\ref{tab:cmp-contrast}) are the substantive comparisons. DGP5's mean is undefined ($-$).}",
 "\\label{tab:cmp-surface}\\begin{tabular}{rr rr rr}\\toprule",
 "&&\\multicolumn{2}{c}{Sun $\\widehat\\mu_0$}&\\multicolumn{2}{c}{ours}\\\\",
 "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}",
 "$n$&$R$&$h$&$0.5h$&Q-LN&Q-INT\\\\\\midrule")
L <- emit(L, 6, function(x) paste0(paste(c(x$n, x$R, k3(x$lim_sun0_h1), k3(x$lim_sun0),
  k3(x$lim_q0), k3(x$lim_qint0)), collapse=" & "), "\\\\"))

# ---- Table 4: DGP5 heavy-tail robustness diagnostics -------------------------
d5 <- d[d$dgp == 5, ]
L <- c(L, "", "\\begin{table}[t]\\centering\\small",
 "\\caption{\\textbf{DGP5 heavy-tail robustness diagnostics} (Student-$t_5$ log-error; true location shift $\\beta_1=0.5$). Sun's mean-slope sampling distribution over $R=500$ replicates vs.\\ our median slope. The ordinary ESD is large and non-monotone in $n$, but the median, IQR, $95$th percentile, and $3$-IQR-trimmed ESD (tESD) are small and \\emph{shrink} with $n$: the inflated ESD is driven by a few replicates with an extreme maximum marker (the $t_5$ log-tail produces markers exceeding $10^4$ in rare replicates), not a coding artifact. Our median slope is unbiased with a small, shrinking ESD throughout. This is an estimand-validity illustration: the conditional mean (Sun's target) is undefined under a $t_5$ log-tail, so it is not a fair mean-estimation contest.}",
 "\\label{tab:cmp-heavy}\\begin{tabular}{r rr rrrr r}\\toprule",
 "&\\multicolumn{2}{c}{Sun $\\widehat\\beta_1$}&\\multicolumn{4}{c}{Sun $\\widehat\\beta_1$ robust spread}&\\\\",
 "\\cmidrule(lr){2-3}\\cmidrule(lr){4-7}",
 "$n$&mean&ESD&med&IQR&$q_{.95}$&tESD&Q ESD$_{.5}$\\\\\\midrule")
for (i in seq_len(nrow(d5))) { x <- d5[i, ]
  L <- c(L, paste0(paste(c(x$n, f2(x$beta_sun1), f2(x$esd_sun1), f2(x$med_sun1), f2(x$iqr_sun1),
    f2(x$q95_sun1), f2(x$tesd_sun1), f3(x$esd_bq50)), collapse=" & "), "\\\\")) }
L <- c(L, "\\bottomrule\\end{tabular}\\end{table}")

writeLines(L, "paper/table_compare2.tex")
cat("wrote paper/table_compare2.tex (4 tables, 5 DGPs)\n")
print(d[, c("dgp","n","beta_sun1","esd_sun1","ctr_sun","ctr_q","ctr_qint","lim_sun0_h1","lim_sun0","lim_q0","tesd_sun1","maxY")])
