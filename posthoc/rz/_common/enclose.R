# Smallest enclosing ellipse per group, for the PCA outlines. Sourced, not run.
#
# ggforce::geom_mark_ellipse can leave points outside (seen on fig 1D), so the
# ellipse is computed here: Khachiyan's minimum-volume enclosing ellipse of
# each group's points, each point first widened to a small circle so every
# dot sits inside with a margin (and three points on a line still get a
# proper ellipse). The fit runs in inches on the panel, so the margin is the
# same on screen whatever the axis ranges; the result maps back to data units.
#
#   enclosing_ellipses(df, "PC1", "PC2", "group", x_per_in, y_per_in, pad_in)
# returns a data frame of polygon vertices (x, y, group) for geom_polygon.

mvee <- function(P, tol = 1e-7, max_iter = 5000) {
  # P: n x 2 points. Returns centre c and matrix A with (p - c)' A (p - c) <= 1.
  n <- nrow(P); Q <- rbind(t(P), rep(1, n)); u <- rep(1 / n, n)
  for (i in seq_len(max_iter)) {
    X <- Q %*% (u * t(Q))
    M <- colSums(Q * solve(X, Q))
    j <- which.max(M); step <- (M[j] - 3) / (3 * (M[j] - 1))
    new_u <- (1 - step) * u; new_u[j] <- new_u[j] + step
    done <- sqrt(sum((new_u - u)^2)) < tol
    u <- new_u
    if (done) break
  }
  c <- as.vector(t(P) %*% u)
  A <- solve(t(P) %*% (u * P) - tcrossprod(c)) / 2
  list(c = c, A = A)
}

enclosing_ellipses <- function(df, x, y, group, x_per_in, y_per_in, pad_in, n = 120) {
  ring <- seq(0, 2 * pi, length.out = 13)[-13]
  out <- lapply(split(df, df[[group]], drop = TRUE), function(g) {
    p <- cbind(g[[x]] / x_per_in, g[[y]] / y_per_in)
    # each point widened to a circle of radius pad_in
    P <- do.call(rbind, lapply(seq_len(nrow(p)), function(i)
      cbind(p[i, 1] + pad_in * cos(ring), p[i, 2] + pad_in * sin(ring))))
    e <- mvee(P)
    # boundary: c + L^-T (cos t, sin t), with A = L L' (Cholesky)
    t <- seq(0, 2 * pi, length.out = n)
    B <- backsolve(chol(e$A), rbind(cos(t), sin(t)))
    data.frame(x = (e$c[1] + B[1, ]) * x_per_in, y = (e$c[2] + B[2, ]) * y_per_in,
               group = g[[group]][1])
  })
  out <- do.call(rbind, out)
  names(out)[names(out) == "group"] <- group
  out
}

# A grey outline layer for a PCA panel. The data-per-inch factors come from
# the point ranges and an estimate of the panel size on the canvas (canvas
# minus axis text and titles at this font scale, minus a bottom legend);
# they only set the shape's aspect and the margin, so an estimate is enough:
# every point is inside the ellipse whatever the factors.
enclosing_outline <- function(df, x, y, group, width, height, font_scale,
                              pad_in = 0.18, legend_in = 0, ...) {
  panel_w <- width  - 0.6  * font_scale
  panel_h <- height - 0.55 * font_scale - legend_in
  xr <- diff(range(df[[x]])) * 1.15; yr <- diff(range(df[[y]])) * 1.3
  ell <- enclosing_ellipses(df, x, y, group, xr / panel_w, yr / panel_h, pad_in * font_scale / 2)
  names(ell)[1:2] <- c(x, y)
  ggplot2::geom_polygon(data = ell,
                        mapping = ggplot2::aes(x = .data[[x]], y = .data[[y]], group = .data[[group]]),
                        colour = "grey55", fill = NA, linewidth = 0.4, inherit.aes = FALSE, ...)
}
