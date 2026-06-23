# ============================================================================
# basis.R — tensor-product B-spline sieve basis for the (u,t) surface.
#
# Pilot choice: REGRESSION-spline sieve (fixed knots, no roughness penalty).
# Complexity is controlled by df per margin; "undersmoothing" (A-US) = more df.
# This is faithful to a sieve estimator and lets rq() return beta_tau directly
# with IPCW weights. A quadratic-penalty smoothing-spline version is deferred.
# ============================================================================
suppressMessages(library(splines))

# Build a tensor-product B-spline design for points (u, t).
# Knots are placed on quantiles of a reference sample (ref_u, ref_t) so the
# basis is stable across bootstrap/CV folds. Returns a matrix WITHOUT intercept
# (the intercept lives in the regression). Degenerate (near-constant / all-zero)
# columns are dropped and the kept-column spec is returned for reuse.
make_tensor_basis <- function(u, t, df_u = 5, df_t = 5,
                              ref_u = u, ref_t = t, degree = 3,
                              keep = NULL) {
  # interior knots from reference quantiles
  ik_u <- if (df_u > degree)
    quantile(ref_u, probs = seq_len(df_u - degree) / (df_u - degree + 1),
             names = FALSE) else numeric(0)
  ik_t <- if (df_t > degree)
    quantile(ref_t, probs = seq_len(df_t - degree) / (df_t - degree + 1),
             names = FALSE) else numeric(0)
  bu_range <- range(ref_u); bt_range <- range(ref_t)

  Bu <- bs(u, knots = ik_u, degree = degree, intercept = TRUE,
           Boundary.knots = bu_range)
  Bt <- bs(t, knots = ik_t, degree = degree, intercept = TRUE,
           Boundary.knots = bt_range)
  # row-wise tensor (Khatri-Rao): column (a,b) = Bu[,a]*Bt[,b]
  nu <- ncol(Bu); nt <- ncol(Bt)
  Z <- matrix(0, nrow = length(u), ncol = nu * nt)
  cc <- 1
  for (a in seq_len(nu)) for (b in seq_len(nt)) {
    Z[, cc] <- Bu[, a] * Bt[, b]; cc <- cc + 1
  }
  if (is.null(keep)) {
    # drop columns with negligible spread on the design points
    sds <- apply(Z, 2, sd)
    keep <- which(sds > 1e-8)
  }
  Z <- Z[, keep, drop = FALSE]
  colnames(Z) <- paste0("B", seq_len(ncol(Z)))
  attr(Z, "keep") <- keep
  attr(Z, "spec") <- list(df_u = df_u, df_t = df_t, degree = degree,
                          ik_u = ik_u, ik_t = ik_t,
                          bu_range = bu_range, bt_range = bt_range, keep = keep)
  Z
}

# Re-evaluate the SAME basis (frozen knots & kept columns) at new points.
eval_tensor_basis <- function(u, t, spec) {
  Bu <- bs(u, knots = spec$ik_u, degree = spec$degree, intercept = TRUE,
           Boundary.knots = spec$bu_range)
  Bt <- bs(t, knots = spec$ik_t, degree = spec$degree, intercept = TRUE,
           Boundary.knots = spec$bt_range)
  nu <- ncol(Bu); nt <- ncol(Bt)
  Z <- matrix(0, nrow = length(u), ncol = nu * nt)
  cc <- 1
  for (a in seq_len(nu)) for (b in seq_len(nt)) {
    Z[, cc] <- Bu[, a] * Bt[, b]; cc <- cc + 1
  }
  Z <- Z[, spec$keep, drop = FALSE]
  colnames(Z) <- paste0("B", seq_len(ncol(Z)))
  Z
}
