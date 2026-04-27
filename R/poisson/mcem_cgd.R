# MCEM with CGD for Poisson outcomes.
# E-step: importance sampling via KFAS::importanceSSM (Poisson observation family).
# M-step: SCAD-penalized CGD for Gamma; glasso for Theta.
#
# Requires:
#   R/common/utils.R       (make_pd, scad_prox, run_cgd_scad)
#   R/common/evaluation.R  (calc_metrics, only used downstream)

library(glasso)
library(KFAS)

# weight cleanup for importanceSSM output
normalize_weights <- function(w, nsim) {
  w <- as.numeric(w)
  w[!is.finite(w)] <- 0
  w[w < 0] <- 0

  if (length(w) != nsim) {
    if (length(w) > nsim) {
      w <- w[1:nsim]
    } else {
      w <- c(w, rep(0, nsim - length(w)))
    }
  }

  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) {
    w <- rep(1 / nsim, nsim)
  } else {
    w <- w / sw
  }

  w
}

poisson_estep_subject <- function(model_subj, T_len, p_dim, nsim, max_try = 3) {
  imp <- NULL
  last_err <- NULL

  for (k in 1:max_try) {
    imp <- tryCatch(
      KFAS::importanceSSM(model_subj, nsim = nsim, antithetics = TRUE),
      error = function(e) {
        last_err <<- e
        NULL
      }
    )
    if (!is.null(imp)) break
  }

  if (is.null(imp)) {
    stop(sprintf("importanceSSM failed after %d tries: %s",
                 max_try,
                 if (is.null(last_err)) "unknown error" else last_err$message))
  }

  samples <- imp$samples
  if (is.null(samples)) stop("importanceSSM returned NULL samples.")

  d <- dim(samples)
  if (length(d) != 3) stop("importanceSSM samples do not have 3 dimensions.")

  # KFAS sometimes returns p x T x nsim instead of T x p x nsim
  if (d[1] == T_len && d[2] == p_dim) {
    samples_use <- samples
  } else if (d[1] == p_dim && d[2] == T_len) {
    nsim_actual <- d[3]
    samples_use <- array(0, dim = c(T_len, p_dim, nsim_actual))
    for (m in 1:nsim_actual) {
      samples_use[, , m] <- t(samples[, , m])
    }
  } else {
    stop("Unexpected samples dimension from importanceSSM.")
  }

  nsim_actual <- dim(samples_use)[3]
  w <- normalize_weights(imp$weights, nsim_actual)
  ess <- 1 / sum(w^2)

  sub_s_xtx <- matrix(0, p_dim, p_dim)
  sub_s_xty <- matrix(0, p_dim, p_dim)
  sub_s_yty <- matrix(0, p_dim, p_dim)

  for (m in 1:nsim_actual) {
    wm <- w[m]
    if (wm <= 0) next

    Z_draw <- samples_use[, , m]
    Z_t   <- Z_draw[2:T_len, , drop = FALSE]
    Z_tm1 <- Z_draw[1:(T_len - 1), , drop = FALSE]

    sub_s_xtx <- sub_s_xtx + wm * crossprod(Z_tm1)
    sub_s_xty <- sub_s_xty + wm * crossprod(Z_t, Z_tm1)
    sub_s_yty <- sub_s_yty + wm * crossprod(Z_t)
  }

  list(s_xtx = sub_s_xtx, s_xty = sub_s_xty, s_yty = sub_s_yty, ess = ess)
}


mcem_cgd_poisson <- function(Y_list,
                             lambda_gamma,
                             lambda_theta,
                             max_iter = 30,
                             n_samples = 50,
                             ebic_gamma = 0.5,
                             tol = 1e-2,
                             gamma_init = NULL,
                             theta_init = NULL,
                             p1_scale = 20) {

  early_lambda_gamma <- 0.05
  early_lambda_theta <- 0.05
  early_iter_gamma <- 5
  early_iter_theta <- 5

  stopifnot(is.list(Y_list), length(Y_list) >= 1)

  n_subjects <- length(Y_list)
  T_len <- nrow(Y_list[[1]])
  p_dim <- ncol(Y_list[[1]])
  N_eff <- n_subjects * (T_len - 1)

  for (i in seq_len(n_subjects)) {
    Y_list[[i]] <- as.matrix(Y_list[[i]])
    stopifnot(nrow(Y_list[[i]]) == T_len, ncol(Y_list[[i]]) == p_dim)
    if (any(Y_list[[i]] < 0)) stop("Poisson outcomes must be nonnegative.")
  }

  Gamma <- if (is.null(gamma_init)) diag(0.5, p_dim) else as.matrix(gamma_init)
  Theta <- if (is.null(theta_init)) diag(1, p_dim) else make_pd(as.matrix(theta_init), ridge = 1e-6)

  Gamma_prev <- Gamma
  Theta_prev <- Theta

  model_list <- lapply(seq_len(n_subjects), function(i) {
    y_subj <- Y_list[[i]]

    KFAS::SSModel(y_subj ~ -1 + SSMcustom(Z = diag(1, p_dim),
                                          T = Gamma,
                                          R = diag(1, p_dim),
                                          Q = diag(1, p_dim),
                                          a1 = rep(0, p_dim),
                                          P1 = diag(p1_scale, p_dim)),
                  distribution = rep("poisson", p_dim))
  })

  conv_hist_gamma <- rep(NA, max_iter)
  conv_hist_theta <- rep(NA, max_iter)
  q_hist <- rep(NA, max_iter)
  ess_hist <- rep(NA, max_iter)

  for (iter in seq_len(max_iter)) {
    Q_mat <- tryCatch(solve(Theta), error = function(e) diag(1, p_dim))
    Q_mat <- make_pd(Q_mat, ridge = 1e-6)

    s_xtx <- matrix(0, p_dim, p_dim)
    s_xty <- matrix(0, p_dim, p_dim)
    s_yty <- matrix(0, p_dim, p_dim)
    ess_vec <- numeric(n_subjects)

    # E-step
    for (subj in seq_len(n_subjects)) {
      model_list[[subj]]$T[, , 1] <- Gamma
      model_list[[subj]]$Q[, , 1] <- Q_mat

      estep_out <- poisson_estep_subject(
        model_subj = model_list[[subj]],
        T_len = T_len,
        p_dim = p_dim,
        nsim = n_samples
      )

      s_xtx <- s_xtx + estep_out$s_xtx
      s_xty <- s_xty + estep_out$s_xty
      s_yty <- s_yty + estep_out$s_yty
      ess_vec[subj] <- estep_out$ess
    }

    s_xtx <- (s_xtx + t(s_xtx)) / 2
    s_yty <- (s_yty + t(s_yty)) / 2
    ess_hist[iter] <- mean(ess_vec)

    # Gamma update
    eff_lam_gamma <- if (iter <= early_iter_gamma) early_lambda_gamma else lambda_gamma

    Beta_new <- run_cgd_scad(
      XtX = s_xtx,
      XtY = t(s_xty),
      Beta_init = t(Gamma),
      lambda = eff_lam_gamma,
      n = N_eff
    )
    Gamma <- t(Beta_new)

    # Theta update
    S_resid <- s_yty - s_xty %*% t(Gamma) - Gamma %*% t(s_xty) + Gamma %*% s_xtx %*% t(Gamma)
    S_resid <- make_pd(S_resid, ridge = 1e-6)
    S_cov <- make_pd(S_resid / N_eff, ridge = 1e-6)

    eff_lam_theta <- if (iter <= early_iter_theta) early_lambda_theta else lambda_theta

    g_fit <- tryCatch(
      glasso::glasso(S_cov, rho = eff_lam_theta),
      error = function(e) {
        glasso::glasso(make_pd(S_cov, ridge = 1e-5), rho = eff_lam_theta)
      }
    )

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

    Gamma_prev <- Gamma
    Theta_prev <- Theta
  }

  conv_hist_gamma <- na.omit(conv_hist_gamma)
  conv_hist_theta <- na.omit(conv_hist_theta)
  q_hist <- na.omit(q_hist)
  ess_hist <- na.omit(ess_hist)

  # approximate BIC using each subject's KFAS log-likelihood at the final params
  df_gamma <- sum(abs(Gamma) > 1e-4)
  df_theta <- sum(abs(Theta[upper.tri(Theta, diag = FALSE)]) > 1e-4)
  df_total <- df_gamma + df_theta

  total_logLik <- 0
  Q_mat_final <- tryCatch(solve(Theta), error = function(e) diag(1, p_dim))
  Q_mat_final <- make_pd(Q_mat_final, ridge = 1e-6)

  for (subj in seq_len(n_subjects)) {
    model_list[[subj]]$T[, , 1] <- Gamma
    model_list[[subj]]$Q[, , 1] <- Q_mat_final

    subj_ll <- tryCatch(
      as.numeric(logLik(model_list[[subj]], nsim = 0)),
      error = function(e) NA_real_
    )

    if (is.na(subj_ll) || !is.finite(subj_ll)) {
      total_logLik <- NA_real_
      break
    } else {
      total_logLik <- total_logLik + subj_ll
    }
  }

  bic_score <- Inf
  if (!is.na(total_logLik) && is.finite(total_logLik)) {
    N_total <- n_subjects * T_len
    bic_score <- -2 * total_logLik +
      df_total * log(N_total) +
      2 * ebic_gamma * df_total * log(p_dim)
  }

  list(Gamma = Gamma, Theta = Theta,
       BIC = bic_score,
       conv_gamma = conv_hist_gamma, conv_theta = conv_hist_theta,
       q_proxy = q_hist, ess_hist = ess_hist)
}
