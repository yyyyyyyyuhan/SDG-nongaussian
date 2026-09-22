# run_one_job_real.R
args <- commandArgs(trailingOnly = TRUE)
lam_g <- as.numeric(args[1])
lam_t <- as.numeric(args[2])

library(Matrix)
source("R/multinomial/mcem_cgd.R")
load("data/processed/hmp2_baseline.RData") 

fit <- mcem_cgd_multinom(Y_list = Y_list,lambda_gamma = lam_g,
                         lambda_theta = lam_t,n_samples = 100,
                         max_iter = 10)

out <- data.frame(
  Lam_G = lam_g, Lam_T = lam_t,
  BIC = fit$BIC,
  N_iter = length(fit$conv_gamma),
  Final_Gamma_diff = tail(fit$conv_gamma, 1),
  Final_Theta_diff = tail(fit$conv_theta, 1),
  Mean_ESS = mean(fit$ess)
)
tag <- sprintf("real_lg_%0.3f_lt_%0.3f", lam_g, lam_t)
csv_fn <- sprintf("res_%s.csv", tag)
write.csv(out, csv_fn, row.names = FALSE)

saveRDS(fit, sprintf("fit_real_lg_%0.4f_lt_%0.4f.rds", lam_g, lam_t))
