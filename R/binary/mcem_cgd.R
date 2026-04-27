# MCEM with CGD for binary outcomes.
# E-step: latent normals truncated by sign of Y, sampled via Kalman simulation smoother.
# M-step: SCAD-penalized CGD for Gamma; glasso for Theta.
#
# Requires:
#   R/common/utils.R       (make_pd, scad_prox, run_cgd_scad)
#   R/common/evaluation.R  (calc_metrics, only used downstream)

library(glasso)
library(truncnorm)
library(KFAS)

mcem_cgd_binary <- function(Y_list,
                            lambda_gamma, lambda_theta,
                            max_iter = 30, n_samples = 50, burn_in = 10,
                            ebic_gamma = 0.5, tol = 1e-2,
                            gamma_init = NULL, theta_init = NULL) {

  # warm-start regularization: small lambdas in early iters
  early_lambda_gamma <- 0.05
  early_lambda_theta <- 0.05
  early_iter_gamma <- 5
  early_iter_theta <- 5
  p1_scale <- 20

  stopifnot(is.list(Y_list), length(Y_list) >= 1)

  n_subjects <- length(Y_list)
  T_len <- nrow(Y_list[[1]])
  p_dim <- ncol(Y_list[[1]])
  N_eff <- n_subjects * (T_len - 1)

  Gamma <- if (is.null(gamma_init)) diag(0.5, p_dim) else as.matrix(gamma_init)
  Theta <- if (is.null(theta_init)) diag(1, p_dim) else make_pd(as.matrix(theta_init), ridge = 1e-6)

  Gamma_prev <- Gamma
  Theta_prev <- Theta

  Z_list <- lapply(seq_len(n_subjects), function(i) matrix(0, T_len, p_dim))
  Y_dummy <- matrix(0, T_len, p_dim)

  model <- KFAS::SSModel(Y_dummy ~ -1 + SSMcustom(Z = diag(1, p_dim),
                                                  T = Gamma,
                                                  R = diag(1, p_dim),
                                                  Q = diag(1, p_dim),
                                                  a1 = rep(0, p_dim),
                                                  P1 = diag(p1_scale, p_dim)),
                         H = diag(1, p_dim))

  conv_hist_gamma <- rep(NA, max_iter)
  conv_hist_theta <- rep(NA, max_iter)
  q_hist <- rep(NA, max_iter)
  final_S_resid <- diag(1, p_dim)

  total_draws <- burn_in + n_samples

  for (iter in seq_len(max_iter)) {
    Q_mat <- tryCatch(solve(Theta), error = function(e) diag(1, p_dim))
    Q_mat <- make_pd(Q_mat, ridge = 1e-6)

    model$T[, , 1] <- Gamma
    model$Q[, , 1] <- Q_mat

    s_xtx <- matrix(0, p_dim, p_dim)
    s_xty <- matrix(0, p_dim, p_dim)
    s_yty <- matrix(0, p_dim, p_dim)

    for (subj in seq_len(n_subjects)) {
      Y_curr <- as.matrix(Y_list[[subj]])
      Z_curr <- Z_list[[subj]]

      sub_s_xtx <- matrix(0, p_dim, p_dim)
      sub_s_xty <- matrix(0, p_dim, p_dim)
      sub_s_yty <- matrix(0, p_dim, p_dim)
      kept <- 0

      for (m in seq_len(total_draws)) {
        # truncated-normal step: latent Z* given current Z_curr and Y
        Z_star <- matrix(0, T_len, p_dim)

        for (j in seq_len(p_dim)) {
          idx1 <- which(Y_curr[, j] == 1)
          idx0 <- which(Y_curr[, j] == 0)

          if (length(idx1) > 0) {
            Z_star[idx1, j] <- truncnorm::rtruncnorm(length(idx1), a = 0, b = Inf,
                                                     mean = Z_curr[idx1, j], sd = 1)
          }

          if (length(idx0) > 0) {
            Z_star[idx0, j] <- truncnorm::rtruncnorm(length(idx0), a = -Inf, b = 0,
                                                     mean = Z_curr[idx0, j], sd = 1)
          }
        }

        model_curr <- model
        model_curr$y <- Z_star

        sim <- KFAS::simulateSSM(model_curr, type = "states", nsim = 1, conditional = TRUE)

        if (length(dim(sim)) == 3) {
          Z_draw <- sim[, , 1]
        } else {
          Z_draw <- sim
        }

        Z_curr <- Z_draw

        # accumulate sufficient statistics post burn-in
        if (m > burn_in) {
          Z_t   <- Z_draw[2:T_len, , drop = FALSE]
          Z_tm1 <- Z_draw[1:(T_len - 1), , drop = FALSE]

          sub_s_xtx <- sub_s_xtx + crossprod(Z_tm1)
          sub_s_xty <- sub_s_xty + crossprod(Z_t, Z_tm1)
          sub_s_yty <- sub_s_yty + crossprod(Z_t)
          kept <- kept + 1
        }
      }

      Z_list[[subj]] <- Z_curr
      s_xtx <- s_xtx + sub_s_xtx / max(kept, 1)
      s_xty <- s_xty + sub_s_xty / max(kept, 1)
      s_yty <- s_yty + sub_s_yty / max(kept, 1)
    }

    s_xtx <- (s_xtx + t(s_xtx)) / 2
    s_yty <- (s_yty + t(s_yty)) / 2

    # Gamma update
    eff_lam_gamma <- if (iter <= early_iter_gamma) early_lambda_gamma else lambda_gamma

    Beta_new <- run_cgd_scad(XtX = s_xtx, XtY = t(s_xty),
                             Beta_init = t(Gamma), lambda = eff_lam_gamma, n = N_eff)
    Gamma <- t(Beta_new)

    # Theta update
    S_resid <- s_yty - s_xty %*% t(Gamma) - Gamma %*% t(s_xty) + Gamma %*% s_xtx %*% t(Gamma)
    S_resid <- make_pd(S_resid, ridge = 1e-6)
    final_S_resid <- S_resid
    S_cov <- make_pd(S_resid / N_eff, ridge = 1e-6)

    eff_lam_theta <- if (iter <= early_iter_theta) early_lambda_theta else lambda_theta

    g_fit <- tryCatch(glasso::glasso(S_cov, rho = eff_lam_theta),
                      error = function(e) {
                        glasso::glasso(make_pd(S_cov, ridge = 1e-6 * 10), rho = eff_lam_theta)
                      })

    Theta <- make_pd(g_fit$wi, ridge = 1e-6)

    # convergence checks
    diff_gamma <- norm(Gamma - Gamma_prev, type = "F") / (norm(Gamma_prev, type = "F") + 1e-10)
    diff_theta <- norm(Theta - Theta_prev, type = "F") / (norm(Theta_prev, type = "F") + 1e-10)

    conv_hist_gamma[iter] <- diff_gamma
    conv_hist_theta[iter] <- diff_theta

    ld <- determinant(Theta, logarithm = TRUE)
    logdet_theta <- if (ld$sign <= 0) NA else as.numeric(ld$modulus)

    q_hist[iter] <- 0.5 * sum(diag(Gamma %*% s_xtx %*% t(Gamma))) -
      sum(diag(Gamma %*% t(s_xty))) +
      N_eff * (sum(diag(S_cov %*% Theta)) - ifelse(is.na(logdet_theta), 0, logdet_theta))

    if (iter >= 5 && diff_gamma < tol && diff_theta < tol) {
      conv_hist_gamma <- conv_hist_gamma[seq_len(iter)]
      conv_hist_theta <- conv_hist_theta[seq_len(iter)]
      q_hist <- q_hist[seq_len(iter)]
      break
    }

    Gamma_prev <- Gamma
    Theta_prev <- Theta
  }

  conv_hist_gamma <- stats::na.omit(conv_hist_gamma)
  conv_hist_theta <- stats::na.omit(conv_hist_theta)
  q_hist <- stats::na.omit(q_hist)

  df_gamma <- sum(abs(Gamma) > 1e-4)
  df_theta <- sum(abs(Theta[upper.tri(Theta, diag = FALSE)]) > 1e-4)
  df_total <- df_gamma + df_theta

  ld_final <- determinant(Theta, logarithm = TRUE)
  logdet_theta <- if (ld_final$sign <= 0) NA else as.numeric(ld_final$modulus)

  bic_score <- Inf
  if (!is.na(logdet_theta)) {
    neg_2_logL_proxy <- N_eff * (sum(diag((final_S_resid / N_eff) %*% Theta)) - logdet_theta)
    bic_score <- neg_2_logL_proxy + df_total * log(N_eff) + 2 * ebic_gamma * df_total * log(p_dim)
  }

  list(Gamma = Gamma, Theta = Theta, BIC = bic_score,
       conv_gamma = conv_hist_gamma, conv_theta = conv_hist_theta, q_proxy = q_hist)
}
