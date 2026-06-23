# gen_compare_tables.R — Phase-1 Sun-vs-quantile tables from results/compare_sun.csv
f3 <- function(x) ifelse(is.na(x), "--", formatC(x, format="f", digits=3))
f2 <- function(x) ifelse(is.na(x), "--", formatC(x, format="f", digits=2))
d <- read.csv("results/compare_sun.csv"); d <- d[order(d$dgp, d$n), ]
g <- function(dg) d[d$dgp == dg, ]

# ---- Table A: DGP1 (mean-correct, Sun-favorable) --------------------------
a <- g(1)
L <- c("\\begin{table}[t]\\centering\\small",
 "\\caption{\\textbf{DGP1 (mean-correct, Sun-favorable).} Conditional-mean model $E[Y|u,t,x]=\\mu_0(u,t)e^{\\beta'x}$ holds exactly and $\\beta_\\tau$ is flat. Sun's profile-kernel mean estimator vs.\\ our quantile estimator with its conditional mean recovered by integrating the fitted quantile process (Q-to-Mean), on the \\emph{same} Monte-Carlo datasets ($R=500$, regime~A, common interior domain). Sun's slope $\\widehat\\beta_1$ recovers the truth ($0.5$); the integrated-quantile $\\beta_\\tau$ is correctly flat (spread $\\approx0$). As expected, \\textbf{Sun is more efficient for the mean surface when its model is correct} ($\\log$-mean IMSE $\\approx2\\times$ smaller), the honest cost of recovering the mean indirectly through quantiles; the mean log-contrast is essentially tied.}",
 "\\label{tab:cmpA}\\begin{tabular}{rr rr rr}\\toprule",
 "&&\\multicolumn{2}{c}{$\\log$-mean IMSE}&&\\\\",
 "\\cmidrule(lr){3-4}",
 "$n$&Sun $\\widehat\\beta_1$&Sun&Q&Sun ctr&Q ctr\\\\\\midrule")
for (i in seq_len(nrow(a))) { x <- a[i,]
  L <- c(L, paste(paste(c(x$n, f2(x$beta_sun1), f3(x$logimse_sun0), f3(x$logimse_q0),
    f3(x$ctr_sun), f3(x$ctr_q)), collapse=" & "), "\\\\")) }
L <- c(L, "\\bottomrule\\end{tabular}\\end{table}")

# ---- Table B: DGP2 (mean true, quantile heterogeneity) --------------------
b <- g(2)
L <- c(L, "", "\\begin{table}[t]\\centering\\small",
 "\\caption{\\textbf{DGP2 (mean true, quantile heterogeneity).} The conditional mean is still proportional, but the covariate's effect on the marker \\emph{distribution} varies with $\\tau$ ($\\beta_{\\tau,1}=\\beta_{{\\rm mean},1}-\\tfrac12[(\\sigma_0{+}\\sigma_1)^2-\\sigma_0^2]+\\sigma_1 z_\\tau$). Sun returns a single mean slope and is \\emph{blind} to this; our estimator recovers the full $\\beta_\\tau$ slope (bias at $\\tau{=}.25/.5/.75$ and the recovered spread $\\widehat\\beta_{.75}-\\widehat\\beta_{.25}$ vs.\\ the truth). \\textbf{Only the quantile method reveals the covariate-driven dispersion.}}",
 "\\label{tab:cmpB}\\begin{tabular}{rr rrr rr}\\toprule",
 "&&\\multicolumn{3}{c}{Q $\\beta_\\tau$ bias}&\\multicolumn{2}{c}{$\\beta_\\tau$ spread}\\\\",
 "\\cmidrule(lr){3-5}\\cmidrule(lr){6-7}",
 "$n$&Sun $\\widehat\\beta_1$&$.25$&$.50$&$.75$&Q&true\\\\\\midrule")
for (i in seq_len(nrow(b))) { x <- b[i,]
  L <- c(L, paste(paste(c(x$n, f2(x$beta_sun1), f3(x$bias_bq25), f3(x$bias_bq50), f3(x$bias_bq75),
    f2(x$spread_q), f2(x$spread_true)), collapse=" & "), "\\\\")) }
L <- c(L, "\\bottomrule\\end{tabular}\\end{table}")

# ---- Table C: DGP3 (quantile true, proportional mean false) ---------------
cc <- g(3)
L <- c(L, "", "\\begin{table}[t]\\centering\\small",
 "\\caption{\\textbf{DGP3 (quantile true, proportional mean false) --- the headline.} Here $\\beta_{\\tau,1}=\\beta_0+\\sigma_1 z_\\tau$ and the induced mean log-contrast $\\Delta(u,t)=\\log\\frac{E[Y|x_1=1]}{E[Y|x_1=0]}=\\beta_0+a(u,t)\\sigma_1+\\tfrac12\\sigma_1^2$ \\emph{varies with $(u,t)$} (true range $\\approx[0.78,1.32]$), so no proportional mean slope exists. Sun's scalar $\\widehat\\beta_1$ is a \\emph{pseudo-average} ($\\approx1.1$--$1.2$). Our quantile estimator recovers $\\beta_\\tau$ (spread $0.66$ vs.\\ true $0.67$) and its induced mean contrast tracks the true surface: \\textbf{the mean-contrast IMSE is significantly lower for Q} ($\\Delta=$ paired Q$-$Sun, MCSE in parentheses). The mean \\emph{level} IMSE still favors Sun (its constant slope captures the average level), but the $(u,t)$-\\emph{structure} of the contrast is what Q recovers and Sun cannot.}",
 "\\label{tab:cmpC}\\begin{tabular}{rr rr rr r}\\toprule",
 "&&\\multicolumn{2}{c}{$\\beta_\\tau$ spread}&\\multicolumn{2}{c}{contrast IMSE}&\\\\",
 "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}",
 "$n$&Sun $\\widehat\\beta_1$&Q&true&Sun&Q&$\\Delta$(MCSE)\\\\\\midrule")
for (i in seq_len(nrow(cc))) { x <- cc[i,]
  L <- c(L, paste(paste(c(x$n, f2(x$beta_sun1), f2(x$spread_q), f2(x$spread_true),
    f3(x$ctr_sun), f3(x$ctr_q), paste0(f3(x$dctr_mean),"(",f3(x$dctr_mcse),")")), collapse=" & "), "\\\\")) }
L <- c(L, "\\bottomrule\\end{tabular}\\end{table}")
writeLines(L, "paper/table_compare_sun.tex")
cat("wrote paper/table_compare_sun.tex (Tables A/B/C)\n")
