# metrics: SEN, SPE, MCC, F1
# off-diagonal upper triangle for symmetric Theta, full matrix for Gamma

calc_metrics <- function(est_mat, true_mat, type = c("Theta", "Gamma"), threshold = 0.01) {
  type <- match.arg(type)

  est_mat  <- as.matrix(est_mat)
  true_mat <- as.matrix(true_mat)

  if (!all(dim(est_mat) == dim(true_mat))) {
    stop("Dimension mismatch")
  }

  if (type == "Theta") {
    est_vec  <- est_mat[upper.tri(est_mat, diag = FALSE)]
    true_vec <- true_mat[upper.tri(true_mat, diag = FALSE)]
  } else {
    est_vec  <- as.vector(est_mat)
    true_vec <- as.vector(true_mat)
  }

  est_bin  <- abs(est_vec) > threshold
  true_bin <- abs(true_vec) > 1e-4

  TP <- as.numeric(sum(est_bin & true_bin))
  TN <- as.numeric(sum(!est_bin & !true_bin))
  FP <- as.numeric(sum(est_bin & !true_bin))
  FN <- as.numeric(sum(!est_bin & true_bin))

  SEN <- ifelse((TP + FN) > 0, TP / (TP + FN), 0)
  SPE <- ifelse((TN + FP) > 0, TN / (TN + FP), 0)
  PRE <- ifelse((TP + FP) > 0, TP / (TP + FP), 0)
  F1  <- ifelse((PRE + SEN) > 0, 2 * PRE * SEN / (PRE + SEN), 0)

  denom <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  MCC <- ifelse(denom > 0, (TP * TN - FP * FN) / denom, 0)

  c(SEN = SEN, SPE = SPE, MCC = MCC, F1 = F1)
}
