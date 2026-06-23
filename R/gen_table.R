# Generate LaTeX Table 1 from a pilot CSV. Usage: Rscript R/gen_table.R <csv> <out.tex>
args <- commandArgs(trailingOnly = TRUE)
infile  <- if (length(args) >= 1) args[1] else "results/pilot_v1_expanded.csv"
outfile <- if (length(args) >= 2) args[2] else "paper/table1.tex"
r <- read.csv(infile)
r <- r[order(r$scenario, r$n, r$tau), ]
has_froz <- "cov_bos_froz" %in% names(r)

f3 <- function(x) ifelse(is.na(x), "--", formatC(x, format = "f", digits = 3))
f2 <- function(x) ifelse(is.na(x), "--", formatC(x, format = "f", digits = 2))

lines <- c(
  "\\begin{table}[t]\\centering\\small",
  "\\caption{Pilot results for the quantile slope $\\beta_{\\tau,1}$, \\textbf{IPCW}",
  "pipeline (IPCW--A-CV, IPCW--A-Rich, IPCW--B-OS-NCF). Bias is",
  "Monte-Carlo over $R=200$ replicates; coverages are over a $150$-replicate",
  "bootstrap subset (Monte-Carlo SE $\\approx0.018$).",
  "``boot'' = full-refit subject bootstrap; ``froz'' = frozen-$G_C$ bootstrap.",
  "A separate paired CC-vs-IPCW study is in Table~\\ref{tab:cc} (distinct",
  "experiment; the IPCW values here are not reused there).}",
  "\\label{tab:pilot}",
  paste0("\\begin{tabular}{rrr r rr r rrr", if (has_froz) "r" else "", "}"),
  "\\toprule",
  paste("Sc. & $n$ & $\\tau$ & $\\beta_{\\tau,1}$ &",
        "bias$_{\\text{ACV}}$ & bias$_{\\text{BOS}}$ & ESD$_{\\text{BOS}}$ &",
        "cov$_{\\text{ACV}}^{\\text{boot}}$ & cov$_{\\text{ARich}}^{\\text{boot}}$ &",
        "cov$_{\\text{BOS}}^{\\text{boot}}$",
        if (has_froz) "& cov$_{\\text{BOS}}^{\\text{froz}}$" else "", "\\\\"),
  "\\midrule"
)
prev <- NULL
for (i in seq_len(nrow(r))) {
  x <- r[i, ]
  blk <- paste(x$scenario, x$n)
  if (!is.null(prev) && prev != blk) lines <- c(lines, "\\addlinespace")
  prev <- blk
  cells <- c(x$scenario, x$n, f2(x$tau), f3(x$true_b1),
             f3(x$bias_acv), f3(x$bias_bos), f3(x$esd_bos),
             f3(x$cov_acv_boot), f3(x$cov_aus_boot), f3(x$cov_bos_boot))
  if (has_froz) cells <- c(cells, f3(x$cov_bos_froz))
  lines <- c(lines, paste(paste(cells, collapse = " & "), "\\\\"))
}
lines <- c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(lines, outfile)
cat("wrote", outfile, "(has_froz =", has_froz, ")\n")
