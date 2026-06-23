depth <- 100  

N <- dim(Y_array)[1]
T_len <- dim(Y_array)[2]
K <- dim(Y_array)[3]   


Y_prop <- Y_array
for (i in 1:N) for (t in 1:T_len) {
  rs <- sum(Y_array[i, t, ], na.rm = TRUE)
  Y_prop[i, t, ] <- if (rs > 0) Y_array[i, t, ] / rs else 1 / K
}

Y_counts <- round(Y_prop * depth)

for (i in 1:N) for (t in 1:T_len) {
  diff <- depth - sum(Y_counts[i, t, ])
  if (diff != 0) {
    max_idx <- which.max(Y_counts[i, t, ])
    Y_counts[i, t, max_idx] <- Y_counts[i, t, max_idx] + diff
  }
}
storage.mode(Y_counts) <- "integer"
stopifnot(all(apply(Y_counts, c(1, 2), sum) == depth))

Y_list <- lapply(1:N, function(i) Y_counts[i, , ])
names(Y_list) <- dimnames(Y_array)[[1]]

save(Y_list, Y_counts, meta_T, N, T_len, K, depth, file = "hmp2_baseline.RData")