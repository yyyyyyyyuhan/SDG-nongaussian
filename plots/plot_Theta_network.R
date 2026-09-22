# Plot the contemporaneous network implied by Theta.

library(corrplot)
library(igraph)


#fit_file <- "fit_real_lg_0.3000_lt_0.0100.rds"
fit <- readRDS("data/processed/fit_real_lg_0.3000_lt_0.0100.rds")
sp10 <- c("F.prau", "B.vulg", "F.plau", "B.unif", "B.ovat","E.rect", "P.dist", "A.hadr", "F.sacc", "B.thet")


Theta <- as.matrix(fit$Theta)
p <- length(sp10)
dimnames(Theta) <- list(sp10, sp10)

## Basic checks
ev <- eigen(Theta, symmetric = TRUE, only.values = TRUE)$values

theta_cutoff <- 0.01

Theta_off <- Theta
diag(Theta_off) <- 0

print(round(Theta, 3))

Sigma <- solve(Theta)
Corr <- cov2cor(Sigma)

dimnames(Sigma) <- list(sp10, sp10)
dimnames(Corr) <- list(sp10, sp10)

print(round(Sigma, 3))
print(round(Corr, 3))

## Strong marginal correlations
strong_pairs <- data.frame(species_1 = character(0),species_2 = character(0),
                           correlation = numeric(0),theta_value = numeric(0),
                           stringsAsFactors = FALSE)

for (i in seq_len(p - 1)) {
  for (j in (i + 1):p) {
    if (abs(Corr[i, j]) > 0.30) {
      strong_pairs[nrow(strong_pairs) + 1, ] <- list(sp10[i],sp10[j],round(Corr[i, j], 3),round(Theta[i, j], 3))
    }
  }
}

if (nrow(strong_pairs) > 0) {
  strong_pairs <- strong_pairs[order(-abs(strong_pairs$correlation)),,drop = FALSE]
  rownames(strong_pairs) <- NULL
}

## Phylum information and colors
phylum <- c("F.prau" = "Firmicutes","B.vulg" = "Bacteroidetes",
            "F.plau" = "Firmicutes","B.unif" = "Bacteroidetes",
            "B.ovat" = "Bacteroidetes","E.rect" = "Firmicutes",
            "P.dist" = "Bacteroidetes","A.hadr" = "Firmicutes",
            "F.sacc" = "Firmicutes","B.thet" = "Bacteroidetes")

## Edge colors
corr_col <- c("Positive correlation" = "#009E49","Negative correlation" = "#D62728")

strong_pairs$phylum_1 <- unname(phylum[strong_pairs$species_1])
strong_pairs$phylum_2 <- unname(phylum[strong_pairs$species_2])
strong_pairs$same_phylum <-strong_pairs$phylum_1 == strong_pairs$phylum_2

print(strong_pairs)

## Construct the undirected network.
edges <- data.frame(from = character(0),to = character(0),
                    theta = numeric(0),corr = numeric(0),
                    stringsAsFactors = FALSE)

for (i in seq_len(p - 1)) {
  for (j in (i + 1):p) {
    if (abs(Theta_off[i, j]) > theta_cutoff) {
      edges[nrow(edges) + 1, ] <- list(sp10[i],sp10[j],Theta_off[i, j],Corr[i, j])
    }
  }
}


gT <- graph_from_data_frame(edges[, c("from", "to")],directed = FALSE,
  vertices = data.frame(name = sp10))

## Use igraph to calculate the layout. 
set.seed(1)
layout_weights <- pmax(abs(edges$corr), 1e-6)
lay <- layout_with_fr(gT, weights = layout_weights)
lay <- norm_coords(lay,xmin = -1.08,xmax = 1.08,ymin = -1.08,ymax = 0.48)
rownames(lay) <- V(gT)$name

push_apart <- function(xy,d_min = 0.54,n_iter = 2000,xmin = -1.08,
                       xmax = 1.08,ymin = -1.08,ymax = 0.48) {
  n <- nrow(xy)
  
  for (it in seq_len(n_iter)) {
    biggest_move <- 0
    for (i in seq_len(n - 1)) {
      for (j in (i + 1):n) {
        z <- xy[i, ] - xy[j, ]
        d <- sqrt(sum(z^2))
        if (d >= d_min) next
        if (d < 1e-8) {
          ang <- 2 * pi * (i + j) / (2 * n)
          u <- c(cos(ang), sin(ang))
          d <- 0
        } else {
          u <- z / d
        }
        
        gap <- d_min - d
        step <- u * (gap / 2 + 1e-4)
        xy[i, ] <- xy[i, ] + step
        xy[j, ] <- xy[j, ] - step
        biggest_move <- max(biggest_move, gap)
      }
    }
    
    xy[, 1] <- pmin(xmax, pmax(xmin, xy[, 1]))
    xy[, 2] <- pmin(ymax, pmax(ymin, xy[, 2]))
    
    if (biggest_move < 1e-4) break
  }
  
  xy
}

lay <- push_apart(lay)
lay <- lay[sp10, , drop = FALSE]

x <- lay[, 1]
y <- lay[, 2]

rad <- rep(0.22, p)
node_col <- unname(phy_col[phylum[sp10]])

max_abs_corr <- max(abs(edges$corr))

edge_lwd <- function(z) {
  if (max_abs_corr == 0) return(1.2)
  pmax(1.2, abs(z) / max_abs_corr * 6)
}

## Draw the correlation heatmap
draw_correlation_heatmap <- function() {
  old_par <- par(mar = c(2, 2, 4, 2), bg = "white")
  on.exit(par(old_par), add = TRUE)
  
  corrplot(Corr,method = "color",type = "upper",
           tl.col = "black",tl.cex = 0.95,
           tl.srt = 45,col = colorRampPalette(c("#009E49", "white", "#D62728"))(100),
           addCoef.col = "black",number.cex = 0.65,
           diag = FALSE,title = "Contemporaneous Correlation (Theta)",
           mar = c(0, 0, 2, 0))
}

## Draw the Theta network
draw_theta_network <- function() {
  old_par <- par(mar = c(0.8, 0.8, 0.8, 0.8),bg = "white",xpd = NA)
  on.exit(par(old_par), add = TRUE)
  
  plot.new()
  plot.window(xlim = c(-1.40, 1.40),ylim = c(-1.40, 1.40),asp = 1)
  
  ## Draw undirected edges 
  for (e in seq_len(nrow(edges))) {
    i <- match(edges$from[e], sp10)
    j <- match(edges$to[e], sp10)
    
    dx <- x[j] - x[i]
    dy <- y[j] - y[i]
    dd <- sqrt(dx^2 + dy^2)
    
    if (dd < 1e-8) next
    
    ux <- dx / dd
    uy <- dy / dd
    
    xs <- x[i] + ux * rad[i]
    ys <- y[i] + uy * rad[i]
    xe <- x[j] - ux * rad[j]
    ye <- y[j] - uy * rad[j]
    
    this_col <- if (edges$corr[e] > 0) {
      unname(corr_col["Positive correlation"])
    } else {
      unname(corr_col["Negative correlation"])
    }
    
    segments(xs, ys, xe, ye,col = this_col,lwd = edge_lwd(edges$corr[e]),lend = "round")
  }
  
  symbols(x, y,circles = rad,add = TRUE,bg = node_col,
         fg = "white",inches = FALSE,lwd = 2)
  
  text(x, y, sp10,cex = 1.55,font = 2,family = "sans")
  
  ## Correlation legend
  legend(x = 1.33,y = 1.33,xjust = 1,yjust = 1,legend = names(corr_col),
         lty = 1,col = unname(corr_col),lwd = 4.5,bty = "n",cex = 1.45,
         y.intersp = 1.20,title.cex = 1.50,title.adj = 0,xpd = NA)
}

## Save the network 
pdf(file = "Theta_network.pdf",width = 8,height = 8,paper = "special",useDingbats = FALSE,bg = "white")
draw_theta_network()
dev.off()


