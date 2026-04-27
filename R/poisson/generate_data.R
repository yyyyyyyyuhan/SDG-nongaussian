# simulate Poisson outcomes
#
# Requires:
# R/common/structures.R


library(mvtnorm)

set.seed(123)

generate_data_poisson <- function(n_subjects = 50,
                                  T_len = 50,
                                  p_dim = 10,
                                  structure_type = "AR1",
                                  intercept = 0) {

  pr <- build_precision(p_dim, structure_type)
  Theta <- pr$Theta
  Omega <- pr$Omega

  Gamma <- build_gamma(p_dim)

  Y_list <- list()
  Z_list <- list()
  Lambda_list <- list()

  for (i in 1:n_subjects) {
    Z_true <- simulate_latent_var1(T_len, p_dim, Gamma, Omega, burn_in = 50)

    log_lambda <- Z_true + intercept
    log_lambda[log_lambda > 10] <- 10
    log_lambda[log_lambda < -10] <- -10

    lambda_mat <- exp(log_lambda)
    Y_mat <- matrix(rpois(length(lambda_mat), lambda = as.vector(lambda_mat)),
                    nrow = T_len, ncol = p_dim)

    Y_list[[i]] <- Y_mat
    Z_list[[i]] <- Z_true
    Lambda_list[[i]] <- lambda_mat
  }

  list(Y_list = Y_list,Z_list = Z_list,
       Lambda_list = Lambda_list,
       True_Omega = Omega,True_Theta = Theta,
       True_Gamma = Gamma,Intercept = intercept)
}
