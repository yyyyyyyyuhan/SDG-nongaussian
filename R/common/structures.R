# Theta / Omega / Gamma generation for the 6 graph structures used across all outcomes.
# Theta is the precision matrix, Omega = Theta^{-1} the covariance matrix.

build_precision <- function(p_dim, structure_type) {
  Theta <- matrix(0, p_dim, p_dim)
  Omega <- matrix(0, p_dim, p_dim)

  if (structure_type %in% c("AR1", "Block")) {
    if (structure_type == "AR1") {
      for (i in 1:p_dim) {
        for (j in 1:p_dim) {
          Omega[i, j] <- 0.7^abs(i - j)
        }
      }

    } else if (structure_type == "Block") {
      split <- floor(p_dim / 2)
      diag(Omega) <- 1
      for (i in 1:split) {
        for (j in 1:split) {
          if (i != j) Omega[i, j] <- 0.5
        }
      }
      for (i in (split + 1):p_dim) {
        for (j in (split + 1):p_dim) {
          if (i != j) Omega[i, j] <- 0.5
        }
      }
    }

    Theta <- solve(Omega)

  } else {
    if (structure_type == "AR2") {
      diag(Theta) <- 1
      val1 <- 0.5
      val2 <- 0.25

      for (i in 1:(p_dim - 1)) Theta[i, i + 1] <- Theta[i + 1, i] <- val1
      for (i in 1:(p_dim - 2)) Theta[i, i + 2] <- Theta[i + 2, i] <- val2

    } else if (structure_type == "Star") {
      diag(Theta) <- 1
      rho <- 0.5
      Theta[1, 2:p_dim] <- Theta[2:p_dim, 1] <- rho
      Theta[1, 1] <- 1 + (p_dim * 0.2)

    } else if (structure_type == "Circle") {
      diag(Theta) <- 2
      val <- 0.9

      for (i in 1:(p_dim - 1)) Theta[i, i + 1] <- Theta[i + 1, i] <- 1
      Theta[1, p_dim] <- Theta[p_dim, 1] <- val

    } else if (structure_type == "Dense") {
      Theta[, ] <- 1
      diag(Theta) <- 2

    } else {
      stop("structure_type must be one of: AR1, AR2, Block, Star, Circle, Dense")
    }

    min_ev <- min(eigen(Theta, symmetric = TRUE)$values)
    if (min_ev <= 0.01) {
      shift <- abs(min_ev) + 0.1
      diag(Theta) <- diag(Theta) + shift
    }

    Omega <- solve(Theta)
  }

  list(Theta = Theta, Omega = Omega)
}

# sparse Gamma 
build_gamma <- function(p_dim, density = 0.2, target_spec_rad = 0.9) {
  Gamma <- matrix(0, p_dim, p_dim)
  n_nz <- round(p_dim * p_dim * density)
  idx <- sample(1:(p_dim * p_dim), n_nz)
  Gamma[idx] <- rnorm(n_nz, 0, 1) * sample(c(1, -1), n_nz, replace = TRUE)

  ev <- max(abs(eigen(Gamma)$values))
  if (ev > 0) Gamma <- Gamma / ev * target_spec_rad

  Gamma
}

# simulate the latent VAR(1) process with Gaussian noise 
simulate_latent_var1 <- function(T_len, p_dim, Gamma, Omega, burn_in = 50) {
  total <- T_len + burn_in
  Z <- matrix(0, total, p_dim)
  Z[1, ] <- rnorm(p_dim)

  noise <- mvtnorm::rmvnorm(total, mean = rep(0, p_dim), sigma = Omega)

  for (t in 2:total) {
    Z[t, ] <- as.numeric(Gamma %*% Z[t - 1, ]) + as.numeric(noise[t, ])
  }

  Z[(burn_in + 1):total, , drop = FALSE]
}
