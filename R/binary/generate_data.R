# simulate binary outcomes from the latent VAR(1) process
#
# Requires:
# R/common/structures.R

library(mvtnorm)

set.seed(123)

generate_data <- function(n_subjects = 50, T_len = 50, p_dim = 10, structure_type = "AR1") {

  pr <- build_precision(p_dim, structure_type)
  Theta <- pr$Theta
  Omega <- pr$Omega

  Gamma <- build_gamma(p_dim)

  Y_list <- list()
  for (i in 1:n_subjects) {
    Z_true <- simulate_latent_var1(T_len, p_dim, Gamma, Omega, burn_in = 50)
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
