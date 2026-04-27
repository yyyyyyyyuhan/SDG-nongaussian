# SLURM job entry: multinomial outcome, one (structure, lam_g, lam_t) combination.
# usage: Rscript cluster/run_one_job_multinomial.R <structure> <lam_g> <lam_t>

args <- commandArgs(trailingOnly = TRUE)
structure <- as.character(args[1])
lam_g <- as.numeric(args[2])
lam_t <- as.numeric(args[3])

suppressPackageStartupMessages({
  library(Matrix)
  library(here)
})

source(here::here("R", "common", "utils.R"))
source(here::here("R", "common", "evaluation.R"))
source(here::here("R", "common", "structures.R"))
source(here::here("R", "multinomial", "mcem_cgd.R"))
source(here::here("R", "multinomial", "generate_data.R"))

data <- generate_data_multinom(n_subjects = 100, T_len = 10, p_dim = 10,
                               structure_type = structure, total_count = 100)

Y_list <- data$Y_list
true_G <- data$True_Gamma
true_T <- data$True_Theta

fit <- mcem_cgd_multinom(Y_list = data$Y_list,
                         lambda_gamma = lam_g,
                         lambda_theta = lam_t,
                         n_samples = 100,
                         max_iter = 100)

mT <- calc_metrics(fit$Theta, true_T, type = "Theta")
mG <- calc_metrics(fit$Gamma, true_G, type = "Gamma")

out <- data.frame(
  Structure = structure,
  Lam_G = lam_g,
  Lam_T = lam_t,
  BIC = fit$BIC,

  T_SEN = as.numeric(mT["SEN"]),
  T_SPE = as.numeric(mT["SPE"]),
  T_MCC = as.numeric(mT["MCC"]),
  T_F1  = as.numeric(mT["F1"]),

  G_SEN = as.numeric(mG["SEN"]),
  G_SPE = as.numeric(mG["SPE"]),
  G_MCC = as.numeric(mG["MCC"]),
  G_F1  = as.numeric(mG["F1"])
)

tag <- sprintf("%s_lg_%0.2f_lt_%0.2f", structure, lam_g, lam_t)
csv_fn <- sprintf("res_%s.csv", tag)
write.csv(out, csv_fn, row.names = FALSE)
cat("Wrote Summary CSV:", csv_fn, "\n")

plot_obj <- list(
  structure = structure,
  lam_g = lam_g,
  lam_t = lam_t,
  conv_gamma = fit$conv_gamma,
  conv_theta = fit$conv_theta,
  final_gamma = fit$Gamma,
  final_theta = fit$Theta,
  true_gamma = true_G,
  true_theta = true_T
)

rds_fn <- sprintf("plot_data_%s.rds", tag)
saveRDS(plot_obj, rds_fn)
cat("Wrote Plotting Data:", rds_fn, "\n")
