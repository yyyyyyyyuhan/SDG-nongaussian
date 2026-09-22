library(igraph)

get_omega_structure <- function(p_dim, structure_type) {
  Omega <- matrix(0, p_dim, p_dim)
  Sigma <- matrix(0, p_dim, p_dim) 
  
  if (structure_type %in% c("AR1", "Block")) {
    if (structure_type == "AR1") {
      for (i in 1:p_dim) {
        for (j in 1:p_dim) {
          Sigma[i, j] <- 0.7^abs(i - j) 
        }
      }
      Omega <- solve(Sigma)
    } else if (structure_type == "Block") {
      split <- floor(p_dim / 2)
      diag(Sigma) <- 1
      Sigma[1:split, 1:split] <- 0.5; diag(Sigma[1:split, 1:split]) <- 1
      Sigma[(split + 1):p_dim, (split + 1):p_dim] <- 0.5; diag(Sigma[(split + 1):p_dim, (split + 1):p_dim]) <- 1
      Omega <- solve(Sigma)
    }
  } else { 
    if (structure_type == "AR2") {
      diag(Omega) <- 1
      val1 <- 0.5; val2 <- 0.25 
      for (i in 1:(p_dim - 1)) Omega[i, i + 1] <- Omega[i + 1, i] <- val1
      for (i in 1:(p_dim - 2)) Omega[i, i + 2] <- Omega[i + 2, i] <- val2
    } else if (structure_type == "Star") {
      diag(Omega) <- 1
      rho <- 0.1 
      Omega[1, 2:p_dim] <- Omega[2:p_dim, 1] <- rho # 节点1连向所有人
      Omega[1, 1] <- 1 
    } else if (structure_type == "Circle") {
      diag(Omega) <- 2
      val <- 0.5 
      for (i in 1:(p_dim - 1)) Omega[i, i + 1] <- Omega[i + 1, i] <- val
      Omega[1, p_dim] <- Omega[p_dim, 1] <- val 
    } else if (structure_type == "Dense") {
      Omega[,] <- 0.1 
      diag(Omega) <- 2   
    }
  }
  return(Omega)
}

plot_structure_fixed <- function(p_dim) {
  types <- c("AR1", "AR2", "Block", "Star", "Circle", "Dense")
  par(mfrow = c(2, 3), mar = c(1, 1, 3, 1))
  
  for (type in types) {
    Theta <- get_omega_structure(p_dim, type)
    diag(Theta) <- 0
    adj_mat <- abs(Theta) > 1e-5
    g <- graph_from_adjacency_matrix(adj_mat, mode = "undirected", diag = FALSE)
    
    l <- matrix(0, nrow = p_dim, ncol = 2) 
    
    if (type == "Star") {
      l[1, ] <- c(0, 0) 
      n_outer <- p_dim - 1
      angles <- seq(0, 2*pi, length.out = n_outer + 1)[1:n_outer]
      
      l[2:p_dim, 1] <- cos(angles) 
      l[2:p_dim, 2] <- sin(angles) 
      
    } else if (type == "Block") {
      split <- floor(p_dim / 2)
      
      angles1 <- seq(0, 2*pi, length.out = split + 1)[1:split]
      l[1:split, 1] <- cos(angles1) - 1.5 
      l[1:split, 2] <- sin(angles1)
      
      n_rest <- p_dim - split
      angles2 <- seq(0, 2*pi, length.out = n_rest + 1)[1:n_rest]
      l[(split + 1):p_dim, 1] <- cos(angles2) + 1.5 
      l[(split + 1):p_dim, 2] <- sin(angles2)
      
    } else {
      l <- layout_in_circle(g)
    }
    v_color <- rep("white", p_dim)

    plot(g,layout = l,
         main = paste0(type),
         vertex.size = 15,vertex.label = 1:p_dim,      
         vertex.label.cex = 0.7,vertex.label.color = "black",
         vertex.color = v_color,vertex.frame.color = "steelblue",
         edge.color = "gray70",edge.width = 1.2)           
  }
}

plot_structure_fixed(10)
