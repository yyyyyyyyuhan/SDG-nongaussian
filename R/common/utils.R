# utility functions shared across binary / poisson / multinomial pipelines

make_pd <- function(S, ridge = 1e-4, eps = 1e-8) {
  S <- (S + t(S)) / 2
  n <- nrow(S)

  S2 <- S + diag(ridge, n)

  ok <- TRUE
  tryCatch(chol(S2), error = function(e) ok <<- FALSE)
  if (ok) return(S2)

  # fall back to eigen-floor when Cholesky fails
  eig <- eigen(S2, symmetric = TRUE)
  vals <- pmax(eig$values, eps)
  S3 <- eig$vectors %*% diag(vals, n) %*% t(eig$vectors)
  (S3 + t(S3)) / 2
}

solve_pd <- function(S, ridge = 1e-6) {
  chol2inv(chol(make_pd(S, ridge = ridge)))
}

scad_prox <- function(z, lambda, a = 3.7) {
  if (lambda <= 1e-10) return(z)

  abs_z <- abs(z)
  res <- z

  idx1 <- abs_z <= 2 * lambda
  if (any(idx1)) {
    res[idx1] <- sign(z[idx1]) * pmax(0, abs_z[idx1] - lambda)
  }

  idx2 <- (abs_z > 2 * lambda) & (abs_z <= a * lambda)
  if (any(idx2)) {
    res[idx2] <- ((a - 1) * z[idx2] - sign(z[idx2]) * a * lambda) / (a - 2)
  }

  res
}

run_cgd_scad <- function(XtX, XtY, Beta_init, lambda, n) {
  p <- ncol(XtX)
  Beta <- Beta_init
  H <- diag(XtX) / n + 1e-5

  for (j in 1:ncol(Beta)) {
    beta_j <- Beta[, j]
    XtX_beta_j <- as.vector(XtX %*% beta_j)
    xty_j <- XtY[, j]

    for (iter in 1:10) {
      max_diff <- 0

      for (k in 1:p) {
        grad_k <- (XtX_beta_j[k] - xty_j[k]) / n
        z <- beta_j[k] - grad_k / H[k]
        beta_new <- scad_prox(z, lambda / H[k])
        diff <- beta_new - beta_j[k]

        if (abs(diff) > 1e-8) {
          beta_j[k] <- beta_new
          XtX_beta_j <- XtX_beta_j + diff * XtX[, k]
          max_diff <- max(max_diff, abs(diff))
        }
      }

      if (max_diff < 1e-5) break
    }

    Beta[, j] <- beta_j
  }

  Beta
}
