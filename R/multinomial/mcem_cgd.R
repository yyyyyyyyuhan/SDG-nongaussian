# MCEM with CGD for multinomial outcomes 
# E-step: multinomial outcomes
# M-step: SCAD-penalized CGD for Gamma; glasso for Theta.
#
# Requires:
# R/common/helper.R   
# R/common/evaluation.R  

library(glasso)
library(KFAS)
library(mvtnorm)

softmax_baseline <- function(z, eps = 1e-8) {
  a <- c(z, 0)
  a <- a - max(a)
  p <- exp(a)
  p <- p / sum(p)
  p <- pmax(p, eps)
  p / sum(p)
}

# init for the latent state
init_Z_from_counts <- function(Y, add = 0.5) {
  T_len <- nrow(Y)
  K <- ncol(Y)
  d <- K - 1
  Z0 <- matrix(0, T_len, d)

  for (t in 1:T_len) {
    base_ct <- Y[t, K] + add
    Z0[t, ] <- log((Y[t, 1:d] + add) / base_ct)
  }
  Z0
}

estep_multinom_mcem <- function(Y, Gamma, Theta, n_samples,
                                Z_init = NULL,m0 = NULL, P0 = NULL,
                                max_inner = 30,tol_inner = 1e-5,
                                ridge = 1e-6) {
  T_len <- nrow(Y)
  K <- ncol(Y)
  d <- K - 1

  if (is.null(Z_init)) Z_ref <- init_Z_from_counts(Y) else Z_ref <- Z_init
  if (is.null(m0)) m0 <- rep(0, d)
  if (is.null(P0)) P0 <- diag(10, d)

  Q <- solve_pd(Theta, ridge = ridge)

  step_size <- 0.3;z_cap <- 8
  divergence_cap <- 1e4;info_floor <- 1e-2

  last_diff <- NA_real_
  model_last <- NULL
  Y_tilde_last <- NULL
  R_list_last <- NULL

  # Laplace to construct gaussian proposal
  for (iter in 1:max_inner) {
    Y_tilde <- matrix(0, T_len, d)
    R_list <- vector("list", T_len)

    for (t in 1:T_len) {
      z0 <- Z_ref[t, ]
      p <- softmax_baseline(z0)
      N_t <- sum(Y[t, ])
      p_sub <- p[1:d]
      g_t <- Y[t, 1:d] - N_t * p_sub
      W_t <- N_t * (diag(p_sub, d) - tcrossprod(p_sub))
      W_t <- 0.5 * (W_t + t(W_t))

      ee <- eigen(W_t, symmetric = TRUE)
      vals <- pmax(ee$values, info_floor)
      W_t_stable <- ee$vectors %*% diag(vals, d) %*% t(ee$vectors)
      W_t_stable <- make_pd(W_t_stable, ridge = ridge)
      W_inv <- solve(W_t_stable)

      Y_tilde[t, ] <- as.numeric(z0 + W_inv %*% g_t)
      R_list[[t]] <- W_inv
    }

    H_array <- array(0, c(d, d, T_len))
    for (t in 1:T_len) {
      H_array[, , t] <- make_pd(R_list[[t]], ridge = ridge)
    }

    model <- KFAS::SSModel(
      Y_tilde ~ -1 + SSMcustom(
        Z = diag(1, d),
        T = Gamma,
        R = diag(1, d),
        Q = make_pd(Q, ridge = ridge),
        a1 = m0,
        P1 = make_pd(P0, ridge = ridge),
        P1inf = matrix(0, d, d),
        n = T_len
      ),
      H = H_array
    )

    out <- KFAS::KFS(model,filtering = "state",smoothing = "state",simplify = FALSE,transform = "ldl")

    alphahat <- as.matrix(out$alphahat)
    if (all(dim(alphahat) == c(T_len, d))) {
      Z_new <- alphahat
    } else if (all(dim(alphahat) == c(d, T_len))) {
      Z_new <- t(alphahat)
    } else {
      stop(sprintf("Dimension is incorrect for alphahat: %s.",
                   paste(dim(alphahat), collapse = " x ")))
    }

    # prevent blowup
    Z_new <- step_size * Z_new + (1 - step_size) * Z_ref
    Z_new[Z_new >  z_cap] <-  z_cap
    Z_new[Z_new < -z_cap] <- -z_cap
    diff_now <- max(abs(Z_new - Z_ref))
    if (!is.finite(diff_now) || diff_now > divergence_cap) {
      Z_new <- pmin(pmax(Z_ref, -z_cap), z_cap)
      diff_now <- max(abs(Z_new - Z_ref))
    }

    Z_ref <- Z_new
    last_diff <- diff_now
    model_last <- model
    Y_tilde_last <- Y_tilde
    R_list_last <- R_list

    if (diff_now < tol_inner) break
  }

  # draw from gaussian posterior
  sims <- KFAS::simulateSSM(object = model_last,type = "states",nsim = n_samples,conditional = TRUE)
  dm <- dim(sims)
  if (is.null(dm) || length(dm) != 3) {
    stop("simulateSSM did not return a 3D array.")
  }

  if (all(dm == c(T_len, d, n_samples))) {
    Z_draws <- sims
  } else if (all(dm == c(d, T_len, n_samples))) {
    Z_draws <- array(0, c(T_len, d, n_samples))
    for (s in 1:n_samples) Z_draws[, , s] <- t(sims[, , s])
  } else {
    stop(sprintf("Dimension is incorrect from simulateSSM: %s.",
                 paste(dm, collapse = " x ")))
  }

  Z_draws[Z_draws >  z_cap] <-  z_cap
  Z_draws[Z_draws < -z_cap] <- -z_cap

  # importance weights
  logw <- numeric(n_samples)
  for (s in 1:n_samples) {
    Zs <- Z_draws[, , s]
    ll_mult <- 0
    ll_pseudo <- 0

    for (t in 1:T_len) {
      p_t <- softmax_baseline(Zs[t, ])
      ll_mult <- ll_mult + sum(Y[t, ] * log(p_t))
      ll_pseudo <- ll_pseudo + mvtnorm::dmvnorm(x = Y_tilde_last[t, ],mean = Zs[t, ],sigma = R_list_last[[t]],log = TRUE)
    }

    logw[s] <- ll_mult - ll_pseudo
  }

  if (all(!is.finite(logw))) {
    weights <- rep(1 / n_samples, n_samples)
  } else {
    m <- max(logw[is.finite(logw)])
    weights <- exp(logw - m)
    weights[!is.finite(weights)] <- 0
    sw <- sum(weights)
    if (!is.finite(sw) || sw <= 0) {
      weights <- rep(1 / n_samples, n_samples)
    } else {
      weights <- weights / sw
    }
  }

  ess <- 1 / sum(weights^2)

  #if ESS collapses
  if (!is.finite(ess) || ess < 5) {
    for (temp in c(2, 5, 10, 20)) {
      if (all(!is.finite(logw))) {
        w_try <- rep(1 / n_samples, n_samples)
      } else {
        logw_temp <- logw / temp
        m <- max(logw_temp[is.finite(logw_temp)])
        w_try <- exp(logw_temp - m)
        w_try[!is.finite(w_try)] <- 0
        sw <- sum(w_try)
        if (!is.finite(sw) || sw <= 0) {
          w_try <- rep(1 / n_samples, n_samples)
        } else {
          w_try <- w_try / sw
        }
      }

      ess_try <- 1 / sum(w_try^2)
      if (is.finite(ess_try) && ess_try > ess) {
        weights <- w_try
        ess <- ess_try
      }
      if (ess >= 5) break
    }
  }

  # weighted posterior mean
  Z_mean <- matrix(0, T_len, d)
  for (s in 1:n_samples) {
    Z_mean <- Z_mean + weights[s] * Z_draws[, , s]
  }

  # weighted sufficient statistics
  s_xtx <- matrix(0, d, d)
  s_xty <- matrix(0, d, d)
  s_yty <- matrix(0, d, d)

  for (s in 1:n_samples) {
    Zs <- Z_draws[, , s]
    ws <- weights[s]

    for (t in 2:T_len) {
      z_prev <- Zs[t - 1, ]
      z_t <- Zs[t, ]
      s_xtx <- s_xtx + ws * tcrossprod(z_prev)
      s_xty <- s_xty + ws * tcrossprod(z_t, z_prev)
      s_yty <- s_yty + ws * tcrossprod(z_t)
    }
  }

  list(Z_mean = Z_mean, Z_mode = Z_ref,
       weights = weights, ess = ess,
       s_xtx = s_xtx, s_xty = s_xty, s_yty = s_yty,
       N_eff = T_len - 1,inner_iter = iter, inner_diff = last_diff)
}

mcem_cgd_multinom <- function(Y_list,lambda_gamma, lambda_theta,
                              n_samples = 100, max_iter = 30, ebic_gamma = 0.5,
                              tol = 1e-2, gamma_init = NULL, theta_init = NULL) {
  Y_list <- lapply(Y_list, as.matrix)
  K <- ncol(Y_list[[1]])
  if (K < 2) stop("Need at least 2 categories.")

  for (i in seq_along(Y_list)) {
    Yi <- Y_list[[i]]
    if (ncol(Yi) != K) stop("All subjects must have the same number of categories.")
    if (nrow(Yi) < 2) stop("Subject has T < 2.")
    if (any(Yi < 0)) stop("Subject has negative counts.")
  }

  #warm start
  early_lambda_gamma <- 0.05;early_lambda_theta <- 0.05
  early_iter_gamma <- 5;early_iter_theta <- 5

  n_subjects <- length(Y_list)
  d <- K - 1

  Gamma <- if (is.null(gamma_init)) diag(0.2, d) else as.matrix(gamma_init)
  Theta <- if (is.null(theta_init)) diag(1, d) else make_pd(as.matrix(theta_init), ridge = ridge)

  Z_init_list <- lapply(Y_list, init_Z_from_counts)
  m0 <- rep(0, d)
  P0 <- diag(10, d)

  Gamma_prev <- Gamma
  Theta_prev <- Theta

  conv_hist_gamma <- rep(NA_real_, max_iter)
  conv_hist_theta <- rep(NA_real_, max_iter)
  q_hist <- rep(NA_real_, max_iter)
  ess_hist <- rep(NA_real_, max_iter)

  final_S_resid <- diag(1, d)
  last_estep_list <- NULL
  N_eff_last <- NA_integer_

  for (iter in seq_len(max_iter)) {

    estep_list <- vector("list", n_subjects)
    s_xtx <- matrix(0, d, d)
    s_xty <- matrix(0, d, d)
    s_yty <- matrix(0, d, d)
    N_eff <- 0
    ess_vec <- numeric(n_subjects)

    # E-step
    for (i in seq_len(n_subjects)) {
      estep_i <- estep_multinom_mcem(Y = Y_list[[i]],
                                     Gamma = Gamma, Theta = Theta,
                                     n_samples = n_samples,Z_init = Z_init_list[[i]],
                                     m0 = m0, P0 = P0)

      estep_list[[i]] <- estep_i
      ess_vec[i] <- estep_i$ess
      s_xtx <- s_xtx + estep_i$s_xtx
      s_xty <- s_xty + estep_i$s_xty
      s_yty <- s_yty + estep_i$s_yty
      N_eff <- N_eff + estep_i$N_eff
    }

    s_xtx <- make_pd(s_xtx, ridge = ridge)
    s_yty <- make_pd(s_yty, ridge = ridge)
    N_eff_last <- N_eff
    ess_hist[iter] <- mean(ess_vec)

    # Gamma Update
    eff_lam_gamma <- if (iter <= early_iter_gamma) early_lambda_gamma else lambda_gamma

    Beta_new <- run_cgd_scad(XtX = s_xtx,XtY = t(s_xty),Beta_init = t(Gamma),lambda = eff_lam_gamma,n = N_eff)
    Gamma <- t(Beta_new)

    # Theta Update
    S_resid <- s_yty -s_xty %*% t(Gamma) -Gamma %*% t(s_xty) + Gamma %*% s_xtx %*% t(Gamma)

    S_resid <- 0.5 * (S_resid + t(S_resid))
    S_cov <- make_pd(S_resid / N_eff, ridge = ridge)

    eff_lam_theta <- if (iter <= early_iter_theta) early_lambda_theta else lambda_theta

    g_fit <- tryCatch(glasso::glasso(S_cov, rho = eff_lam_theta),
      error = function(e) glasso::glasso(make_pd(S_cov, ridge = 1e-5), rho = eff_lam_theta))

    Theta <- 0.5 * (g_fit$wi + t(g_fit$wi))
    diag(Theta) <- diag(Theta) + ridge
    final_S_resid <- S_resid

    ld <- determinant(Theta, logarithm = TRUE)
    logdet_theta <- if (ld$sign <= 0) NA_real_ else as.numeric(ld$modulus)

    q_hist[iter] <- 0.5 * sum(diag(Gamma %*% s_xtx %*% t(Gamma))) -
      sum(diag(Gamma %*% t(s_xty))) + N_eff * (sum(diag(S_cov %*% Theta)) -ifelse(is.na(logdet_theta), 0, logdet_theta))
    
    #convergence check
    diff_gamma <- norm(Gamma - Gamma_prev, type = "F") / (norm(Gamma_prev, type = "F") + 1e-10)
    diff_theta <- norm(Theta - Theta_prev, type = "F") / (norm(Theta_prev, type = "F") + 1e-10)

    conv_hist_gamma[iter] <- diff_gamma
    conv_hist_theta[iter] <- diff_theta

    last_estep_list <- estep_list
    Z_init_list <- lapply(estep_list, function(x) x$Z_mean)

    Gamma_prev <- Gamma
    Theta_prev <- Theta
  }

  conv_hist_gamma <- stats::na.omit(conv_hist_gamma)
  conv_hist_theta <- stats::na.omit(conv_hist_theta)
  q_hist <- stats::na.omit(q_hist)
  ess_hist <- stats::na.omit(ess_hist)

  df_gamma <- sum(abs(Gamma) > 1e-4)
  df_theta <- sum(abs(Theta[upper.tri(Theta, diag = FALSE)]) > 1e-4)
  df_total <- df_gamma + df_theta

  ld_final <- determinant(Theta, logarithm = TRUE)
  logdet_theta <- if (ld_final$sign <= 0) NA_real_ else as.numeric(ld_final$modulus)

  #Compute BIC via ebic
  bic_score <- Inf
  if (!is.na(logdet_theta) && !is.na(N_eff_last) && N_eff_last > 0) {
    neg_2_logL_proxy <- N_eff_last *
      (sum(diag((final_S_resid / N_eff_last) %*% Theta)) - logdet_theta)
    bic_score <- neg_2_logL_proxy +df_total * log(N_eff_last) + 2 * ebic_gamma * df_total * log(d)
  }

  subject_fits <- vector("list", n_subjects)
  for (i in seq_len(n_subjects)) {
    Zi <- last_estep_list[[i]]$Z_mean
    Pi <- t(apply(Zi, 1, softmax_baseline))

    subject_fits[[i]] <- list(Z_mean = Zi,fitted_prob = Pi,ess = last_estep_list[[i]]$ess)
  }

  list(Gamma = Gamma, Theta = Theta,BIC = bic_score,
       conv_gamma = conv_hist_gamma, conv_theta = conv_hist_theta)
}
