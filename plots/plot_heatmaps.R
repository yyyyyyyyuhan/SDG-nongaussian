get_precision <- function(p_dim, structure_type) {
  Theta <- matrix(0, p_dim, p_dim)
  Omega <- matrix(0, p_dim, p_dim)

  if (structure_type %in% c("AR1", "Block")) {
    if (structure_type == "AR1") {
      for (i in 1:p_dim) {
        for (j in 1:p_dim) {
          Omega[i, j] <- 0.7^abs(i - j)
        }
      }
    } else if (structure_type == "Block") {
      split <- floor(p_dim / 2)
      diag(Omega) <- 1
      for (i in 1:split) {
        for (j in 1:split) {
          if (i != j) Omega[i, j] <- 0.5
        }
      }
      for (i in (split + 1):p_dim) {
        for (j in (split + 1):p_dim) {
          if (i != j) Omega[i, j] <- 0.5
        }
      }
    }
    Theta <- solve(Omega)
  } else {
    if (structure_type == "AR2") {
      diag(Theta) <- 1
      val1 <- 0.5; val2 <- 0.25
      for (i in 1:(p_dim - 1)) Theta[i, i + 1] <- Theta[i + 1, i] <- val1
      for (i in 1:(p_dim - 2)) Theta[i, i + 2] <- Theta[i + 2, i] <- val2
    } else if (structure_type == "Star") {
      diag(Theta) <- 1
      rho <- 0.5
      Theta[1, 2:p_dim] <- Theta[2:p_dim, 1] <- rho
      Theta[1, 1] <- 1 + (p_dim * 0.2)
    } else if (structure_type == "Circle") {
      diag(Theta) <- 2
      val <- 0.9
      for (i in 1:(p_dim - 1)) Theta[i, i + 1] <- Theta[i + 1, i] <- 1
      Theta[1, p_dim] <- Theta[p_dim, 1] <- val
    } else if (structure_type == "Dense") {
      Theta[, ] <- 1
      diag(Theta) <- 2
    } else {
      stop("structure_type must be one of: AR1, AR2, Block, Star, Circle, Dense")
    }
    min_ev <- min(eigen(Theta, symmetric = TRUE)$values)
    if (min_ev <= 0.01) {
      shift <- abs(min_ev) + 0.1
      diag(Theta) <- diag(Theta) + shift
    }
  }
  Theta
}

## Plot precision heatmaps
plot_precision_heatmaps <- function(p_dim = 10,file = NULL,
                                    width = 17,height = 11) {
  structures <- c("AR1", "AR2", "Block", "Star", "Circle", "Dense")
  mats <- lapply(structures, function(s) get_precision(p_dim, s))
  names(mats) <- structures

  global_max <- max(sapply(mats, function(M) max(abs(M))))
  zlim <- c(-global_max, global_max)
  n_cols <- 200
  cols <- colorRampPalette(c("#3B4CC0", "white", "#B40426"))(n_cols)

  axis_cex <- if (p_dim <= 10) 2.35 else if (p_dim <= 20) 1.75 else 1.25
  title_cex <- 3.8
  colorbar_cex <- if (p_dim <= 20) 2.35 else 1.7
  tick_len <- if (p_dim <= 20) -0.025 else -0.015

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  # 2x3 heatmaps 
  layout(matrix(c(1, 2, 3, 7,
                   4, 5, 6, 7), nrow = 2, byrow = TRUE),
         widths  = c(1, 1, 1, 0.25),
         heights = c(1, 1))

  for (s in structures) {
    Theta <- mats[[s]]
    p <- nrow(Theta)

    par(mar = c(4.2, 4.2, 4.8, 1.1))
    image(x = 1:p, y = 1:p,z = t(Theta[p:1, ]),
          col  = cols,zlim = zlim,axes = FALSE,
          xlab = "", ylab = "", main = "")

    if (p_dim <= 20) {
      at_x <- 1:p; lab_x <- 1:p
      at_y <- 1:p; lab_y <- p:1
    } else {
      at_x <- pretty(1:p, n = 6); lab_x <- at_x
      at_y <- pretty(1:p, n = 6); lab_y <- (p + 1) - at_y
    }
    axis(1, at = at_x, labels = lab_x, cex.axis = axis_cex, tck = tick_len,
         lwd = 1.2, lwd.ticks = 1.2, gap.axis = -1)
    axis(2, at = at_y, labels = lab_y, cex.axis = axis_cex, tck = tick_len,
         lwd = 1.2, lwd.ticks = 1.2, gap.axis = -1)
    box(col = "grey40", lwd = 1.2)
    title(main = s, cex.main = title_cex, font.main = 1.5, line = 0.35)
  }

  # --- shared colorbar ---
  par(mar = c(4.2, 1.1, 4.8, 4.3))
  zseq <- seq(zlim[1], zlim[2], length.out = n_cols + 1)
  image(x = 1, y = zseq,z = matrix(zseq[-1], nrow = 1),
        col  = cols,axes = FALSE,xlab = "", ylab = "", main = "")
  ticks <- pretty(zlim, n = 5)
  axis(4, at = ticks, labels = round(ticks, 2), las = 1, cex.axis = colorbar_cex,
       lwd = 1.2, lwd.ticks = 1.2)
  box(col = "grey40", lwd = 1.2)
}

plot_precision_heatmaps(p_dim = 10,  file = "heatmap_p10.pdf")
