# Plot a Gamma network 
fit <- readRDS("data/processed/fit_real_lg_0.3000_lt_0.0100.rds")
sp10 <- c("F.prau", "B.vulg", "F.plau", "B.unif", "B.ovat","E.rect", "P.dist", "A.hadr", "F.sacc", "B.thet")
G <- fit$Gamma

rownames(G) <- sp10
colnames(G) <- sp10

G0 <- G
diag(G0) <- 0

cutoff <- 0.001
p <- length(sp10)

## Colors by phylum
firm <- c("F.prau", "F.plau", "F.sacc", "A.hadr", "E.rect")
phy <- rep("Bacteroidetes", p)
phy[sp10 %in% firm] <- "Firmicutes"

phy_col <- c("Firmicutes" = "#E67E22","Bacteroidetes" = "#3498DB")
node_col <- unname(phy_col[phy])

get_layout <- function(W, phy, seed = 7, niter = 400) {
  set.seed(seed)
  m <- nrow(W)
  xy <- matrix(runif(2 * m, -0.5, 0.5), ncol = 2)
  xy[phy == "Firmicutes", 2] <- xy[phy == "Firmicutes", 2] + 0.6
  xy[phy == "Bacteroidetes", 2] <- xy[phy == "Bacteroidetes", 2] - 0.6
  
  area <- 4
  k <- sqrt(area / m)
  temp <- 0.25
  
  for (it in seq_len(niter)) {
    move <- matrix(0, m, 2)
        for (i in seq_len(m)) {
      for (j in seq_len(m)) {
        if (i == j) next
        
        v <- xy[i, ] - xy[j, ]
        d <- max(sqrt(sum(v^2)), 1e-6)
        move[i, ] <- move[i, ] + v / d * k^2 / d
      }
    }
    
    for (i in seq_len(m - 1)) {
      for (j in (i + 1):m) {
        v <- xy[i, ] - xy[j, ]
        d <- max(sqrt(sum(v^2)), 1e-6)
        
        a <- W[i, j]
        if (phy[i] == phy[j]) a <- a + 0.18
        if (a == 0) next
        
        f <- d^2 / k * a
        move[i, ] <- move[i, ] - v / d * f
        move[j, ] <- move[j, ] + v / d * f
      }
    }
    
    for (i in seq_len(m)) {
      dd <- max(sqrt(sum(move[i, ]^2)), 1e-6)
      xy[i, ] <- xy[i, ] + move[i, ] / dd * min(dd, temp)
    }
    
    temp <- temp * 0.96
  }
  
  xy[, 1] <-(xy[, 1] - mean(xy[, 1])) /max(abs(xy[, 1] - mean(xy[, 1]))) * 0.9
  xy[, 2] <-(xy[, 2] - mean(xy[, 2])) /max(abs(xy[, 2] - mean(xy[, 2]))) * 0.9
  
  xy
}

push_nodes <- function(xy, min.dist = 0.54, niter = 1000) {
  m <- nrow(xy)
  
  for (it in seq_len(niter)) {
    biggest <- 0
    
    for (i in seq_len(m - 1)) {
      for (j in (i + 1):m) {
        v <- xy[i, ] - xy[j, ]
        d <- sqrt(sum(v^2))
        
        if (d >= min.dist) next
        
        if (d < 1e-8) {
          ang <- 2 * pi * (i + j) / (2 * m)
          u <- c(cos(ang), sin(ang))
          d <- 0
        } else {
          u <- v / d
        }
        
        gap <- min.dist - d
        xy[i, ] <- xy[i, ] + u * (gap / 2 + 1e-4)
        xy[j, ] <- xy[j, ] - u * (gap / 2 + 1e-4)
        biggest <- max(biggest, gap)
      }
    }
    
    xy[, 1] <- pmin(1.08, pmax(-1.08, xy[, 1]))
    xy[, 2] <- pmin(1.08, pmax(-1.08, xy[, 2]))
    
    if (biggest < 1e-4) break
  }
  
  xy
}

push_pair <- function(xy, nm, a, b, d0 = 0.68) {
  i <- match(a, nm)
  j <- match(b, nm)
  
  v <- xy[i, ] - xy[j, ]
  d <- sqrt(sum(v^2))
  
  if (d >= d0) return(xy)
  
  if (d < 1e-8) {
    u <- c(1, 0)
    d <- 0
  } else {
    u <- v / d
  }
  
  xy[i, ] <- xy[i, ] + u * (d0 - d) / 2
  xy[j, ] <- xy[j, ] - u * (d0 - d) / 2
  
  xy
}

W <- (abs(G0) + t(abs(G0))) / 2
W[W < cutoff] <- 0

if (max(W) > 0) {
  W <- W / max(W)
}

xy <- get_layout(W, phy, seed = 7)
xy <- push_nodes(xy, min.dist = 0.54)
xy <- push_pair(xy, sp10, "F.prau", "B.ovat", d0 = 0.74)
xy <- push_pair(xy, sp10, "E.rect", "F.prau", d0 = 0.74)
xy <- push_nodes(xy, min.dist = 0.54)

x <- xy[, 1]
y <- xy[, 2]
rad <- rep(0.22, p)

max_off <- max(abs(G0))

edge_lwd <- function(z) {
  if (max_off == 0) return(1.2)
  pmax(1.2, abs(z) / max_off * 6)
}

draw_gamma_network <- function() {
  old_par <- par(mar = c(0.8, 0.8, 0.8, 0.8),bg = "white",xpd = NA)
  on.exit(par(old_par), add = TRUE)
  
  plot.new()
  plot.window(xlim = c(-1.40, 1.40),ylim = c(-1.40, 1.40),asp = 1)
  
  ## Directed off-diagonal edges
  for (r in seq_len(p)) {
    for (c in seq_len(p)) {
      if (r == c) next
      
      val <- G0[r, c]
      if (abs(val) <= cutoff) next
      
      x0 <- x[c]
      y0 <- y[c]
      x1 <- x[r]
      y1 <- y[r]
      
      dx <- x1 - x0
      dy <- y1 - y0
      dd <- sqrt(dx^2 + dy^2)
      
      if (dd < 1e-8) next
      
      ux <- dx / dd
      uy <- dy / dd
      
      xs <- x0 + ux * rad[c]
      ys <- y0 + uy * rad[c]
      xe <- x1 - ux * rad[r]
      ye <- y1 - uy * rad[r]
      
      this_lwd <- edge_lwd(val)
      this_col <- if (val > 0) "#2C3E50" else "#C0392B"
      reverse_edge <- abs(G0[c, r]) > cutoff
      
      if (reverse_edge) {
        off <- 0.07
        mx <- (xs + xe) / 2 - uy * off
        my <- (ys + ye) / 2 + ux * off
        
        xspline(
          c(xs, mx, xe), c(ys, my, ye),
          shape = 1,
          open = TRUE,
          border = this_col,
          lwd = this_lwd
        )
        
        tx <- xe - mx
        ty <- ye - my
        td <- sqrt(tx^2 + ty^2)
        
        arrows(
          xe - tx / td * 0.04,
          ye - ty / td * 0.04,
          xe,
          ye,
          length = 0.11,
          angle = 25,
          lwd = this_lwd,
          col = this_col
        )
      } else {
        arrows(xs, ys, xe, ye,length = 0.11,angle = 25,lwd = this_lwd,col = this_col)
      }
    }
  }
  
  ## Nodes
  symbols(x, y,circles = rad,add = TRUE,bg = node_col,fg = "white",inches = FALSE,lwd = 2)
  
  text(x, y, sp10, cex = 1.55, font = 2)
  
  legend(x = -1.33,y = 1.33,xjust = 0,yjust = 1,
         legend = names(phy_col),pch = 21,
         pt.bg = unname(phy_col),col = "black",
         pt.cex = 3,pt.lwd = 1,bty = "n",cex = 1.45,
         y.intersp = 1.20,title = "Phylum",title.adj = 0,
         xpd = NA)
}

## Save the network 
pdf(file = "Gamma_network.pdf",width = 8,height = 8,paper = "special",useDingbats = FALSE,bg = "white")
draw_gamma_network()
dev.off()
