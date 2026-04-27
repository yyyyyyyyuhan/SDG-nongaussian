# simulate multinomial outcomes 
# number of categories K = p_dim + 1.
#
# Requires:
# R/common/structures.R

library(mvtnorm)

set.seed(123)

softmax_baseline <- function(z, eps = 1e-8) {
  a <- c(z, 0)
  a <- a - max(a)
  p <- exp(a)
  p <- p / sum(p)
  p <- pmax(p, eps)
  p / sum(p)
}

generate_data_multinom <- function(n_subjects = 50,T_len = 50,p_dim = 10,
                                   structure_type = "AR1",total_count = 20,
                                   burn_in = 50,intercept = 0,
                                   eta = 10) {

  pr <- build_precision(p_dim, structure_type)
  Theta <- pr$Theta
  Omega <- pr$Omega

  Gamma <- build_gamma(p_dim)

  K <- p_dim + 1

  if (length(total_count) == 1) total_count <- rep(total_count, T_len)
  if (length(total_count) != T_len) stop("total_count must be scalar or length T_len.")

  if (length(intercept) == 1) intercept <- rep(intercept, p_dim)
  if (length(intercept) != p_dim) stop("intercept must be scalar or length p_dim.")

  Y_list <- vector("list", n_subjects)
  Z_list <- vector("list", n_subjects)
  P_list <- vector("list", n_subjects)

  for (subj in 1:n_subjects) {
    Z_true <- simulate_latent_var1(T_len, p_dim, Gamma, Omega, burn_in = burn_in)
    Eta <- sweep(Z_true, 2, intercept, FUN = "+")

    if (!is.null(eta)) {
      Eta[Eta >  eta] <-  eta
      Eta[Eta < -eta] <- -eta
    }

    P_mat <- t(apply(Eta, 1, softmax_baseline))

    Y_mat <- matrix(0, T_len, K)
    for (t in 1:T_len) {
      Y_mat[t, ] <- as.numeric(rmultinom(1, size = total_count[t], prob = P_mat[t, ]))
    }

    Y_list[[subj]] <- Y_mat
    Z_list[[subj]] <- Z_true
    P_list[[subj]] <- P_mat
  }

  list(Y_list = Y_list,Z_list = Z_list,P_list = P_list,
       True_Omega = Omega,True_Theta = Theta,
       True_Gamma = Gamma,
       Total_Count = total_count,Intercept = intercept,
       K = K,p_dim = p_dim)
}
