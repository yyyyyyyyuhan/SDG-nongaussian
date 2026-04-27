# simulate binary outcomes from the latent VAR(1) process via probit-style thresholding.
# requires R/common/structures.R to be sourced.

library(mvtnorm)

set.seed(123)

generate_data <- function(n_subjects = 50, T_len = 50, p_dim = 10, structure_type = "AR1") {

  pr <- build_precision(p_dim, structure_type)
  Theta <- pr$Theta
  Omega <- pr$Omega

  Gamma <- build_gamma(p_dim)

  # binary case: original code used burn_in = 50, embedded inline rather than via the helper,
  # so we keep the same numbers via simulate_latent_var1(burn_in = 50).
  Y_list <- list()
  for (i in 1:n_subjects) {
    Z_true <- simulate_latent_var1(T_len, p_dim, Gamma, Omega, burn_in = 50)
    # additive measurement noise then thresholding at 0
    W <- Z_true + matrix(rnorm(T_len * p_dim), T_len, p_dim)
    Y_list[[i]] <- (W > 0) * 1
  }

  list(
    Y_list = Y_list,
    True_Omega = Omega,
    True_Theta = Theta,
    True_Gamma = Gamma
  )
}
