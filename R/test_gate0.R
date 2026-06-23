# Gate 0a: with KNOWN mu(u,t) offset, plain weighted QR on the observed
# event-subject markers must recover analytic beta_tau. Isolates estimand+QR
# from spline approximation, censoring, and orthogonalization.
suppressMessages(library(quantreg))
source("R/dgp.R")

set.seed(20260613)
taus <- c(0.10, 0.25, 0.50, 0.75, 0.90)

run_gate0a <- function(scenario, n = 4000) {
  d <- generate_data(n, scenario = scenario)
  cat(sprintf("\n== Scenario %d ==  ", scenario))
  print(d$summary)
  mk <- d$mk
  mk <- mk[mk$dT == 1 & mk$U > 0 & mk$U < mk$T, ]   # target: observed onset, pre-onset
  mk$resid <- mk$L - mu_surface(mk$U, mk$T)         # known-mu offset
  p <- d$params
  cat(sprintf("  event-subject marker rows used: %d\n", nrow(mk)))
  cat(sprintf("  %-5s %9s %9s | %9s %9s | %9s %9s\n",
              "tau", "b1.hat", "b1.true", "b2.hat", "b2.true", "int.hat", "int.true"))
  for (tau in taus) {
    fit <- rq(resid ~ X1 + X2, tau = tau, data = mk)
    cf  <- coef(fit)
    tb  <- true_beta(tau, p)
    cat(sprintf("  %-5.2f %9.4f %9.4f | %9.4f %9.4f | %9.4f %9.4f\n",
                tau, cf["X1"], tb["beta1"], cf["X2"], tb["beta2"],
                cf["(Intercept)"], true_intercept_shift(tau, p)))
  }
}

run_gate0a(1)
run_gate0a(2)
