library(corrplot)
library(igraph)

fit <- readRDS("fit_real_lg_0.3000_lt_0.0100.rds")

Theta <- fit$Theta
dimnames(Theta) <- list(sp10, sp10)

## basic check
ev <- eigen(Theta)$values

dim(Theta); min(ev)
max(ev); max(ev) / min(ev)

Theta_off <- Theta
diag(Theta_off) <- 0

sum(abs(Theta_off[upper.tri(Theta_off)]) > 0.01)
choose(10, 2)

round(Theta, 3)

## covariance and correlation implied by Theta
Sigma <- solve(Theta)
Corr <- cov2cor(Sigma)

dimnames(Sigma) <- list(sp10, sp10)
dimnames(Corr) <- list(sp10, sp10)

round(Sigma, 3)
round(Corr, 3)

## correlation heatmap
par(mar = c(2, 2, 4, 2))

corrplot(Corr, method = "color",type = "upper",
         tl.col = "black",tl.cex = 0.95,tl.srt = 45,
         col = colorRampPalette(c("#3498DB", "white", "#E74C3C"))(100),
         addCoef.col = "black",number.cex = 0.65, diag = FALSE,
         title = "Contemporaneous Correlation (Theta)",
         mar = c(0, 0, 2, 0))

## large marginal correlations
strong_pairs <- NULL

for (i in 1:9) {
  for (j in (i + 1):10) {
    if (abs(Corr[i, j]) > 0.30) {
      tmp <- data.frame(
        species_1 = sp10[i],
        species_2 = sp10[j],
        correlation = round(Corr[i, j], 3),
        theta_value = round(Theta[i, j], 3)
      )
      strong_pairs <- rbind(strong_pairs, tmp)
    }
  }
}

strong_pairs <- strong_pairs[order(-abs(strong_pairs$correlation)), ]
rownames(strong_pairs) <- NULL
strong_pairs

## phylum information
phylum <- c("F.prau" = "Firmicutes", "B.vulg" = "Bacteroidetes",
            "F.plau" = "Firmicutes", "B.unif" = "Bacteroidetes",
            "B.ovat" = "Bacteroidetes", "E.rect" = "Firmicutes",
            "P.dist" = "Bacteroidetes", "A.hadr" = "Firmicutes",
            "F.sacc" = "Firmicutes", "B.thet" = "Bacteroidetes")

phy_col <- c(Firmicutes = "#E67E22",Bacteroidetes = "#3498DB")
corr_col <- c("Positive correlation" = "#009E49","Negative correlation" = "#D62728")

strong_pairs$phylum_1 <- phylum[strong_pairs$species_1]
strong_pairs$phylum_2 <- phylum[strong_pairs$species_2]
strong_pairs$same_phylum <- strong_pairs$phylum_1 == strong_pairs$phylum_2

strong_pairs

## network from nonzero Theta entries
edges <- NULL

for (i in 1:9) {
  for (j in (i + 1):10) {
    if (abs(Theta_off[i, j]) > 0.01) {
      edges <- rbind(edges,data.frame(from = sp10[i],to = sp10[j],
                                      theta = Theta_off[i, j],corr = Corr[i, j]))
    }
  }
}

gT <- graph_from_data_frame(edges[, c("from", "to")],
                            directed = FALSE,vertices = data.frame(name = sp10))

V(gT)$phylum <- phylum[V(gT)$name]
V(gT)$color <- unname(phy_col[V(gT)$phylum])
V(gT)$size <- 38
V(gT)$degree <- degree(gT)

E(gT)$theta <- edges$theta
E(gT)$corr <- edges$corr
E(gT)$color <- ifelse(E(gT)$corr > 0, corr_col[1], corr_col[2])
E(gT)$width <- pmax(2.6, abs(E(gT)$corr) * 7.5)

push_apart <- function(xy, d_min = 0.52, n_iter = 2000) {
  n <- nrow(xy)
  for (it in 1:n_iter) {
    moved <- 0
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        z <- xy[i, ] - xy[j, ]
        d <- sqrt(sum(z^2))
        if (d >= d_min) next
  
        if (d < 1e-8) {
          ang <- 2 * pi * (i + j) / (2 * n)
          z <- c(cos(ang), sin(ang))
          d <- 1
        } else {
          z <- z / d
        }
        
        step <- z * ((d_min - d) / 2 + 1e-4)
        xy[i, ] <- xy[i, ] + step
        xy[j, ] <- xy[j, ] - step
        
        moved <- max(moved, d_min - d)
      }
    }
    
    xy[, 1] <- pmin(1.04, pmax(-1.04, xy[, 1]))
    xy[, 2] <- pmin(0.25, pmax(-1.05, xy[, 2]))
    
    if (moved < 1e-4) break
  }
  xy
}

set.seed(1)

lay <- layout_with_fr(gT, weights = abs(E(gT)$corr))
lay <- norm_coords(lay, xmin = -1, xmax = 1, ymin = -1, ymax = 0.25)
lay <- push_apart(lay)

par(mar = c(2, 2, 2, 2), bg = "white")

plot(gT,layout = lay,rescale = FALSE,
     xlim = c(-1.2, 1.2),ylim = c(-1.2, 1.2),
     vertex.size = V(gT)$size,vertex.color = V(gT)$color,
     vertex.frame.color = "white",vertex.frame.width = 2,
     vertex.label = V(gT)$name,vertex.label.color = "black",
     vertex.label.cex = 1.43,vertex.label.font = 2,
     vertex.label.family = "sans",edge.width = E(gT)$width,
     edge.color = E(gT)$color,edge.curved = 0.1)

legend("topleft",inset = c(0.05, 0.13),legend = names(phy_col),
       pch = 21,pt.bg = phy_col,col = "black",pt.cex = 3.1,
       bty = "n",cex = 1.55,title = "Phylum",title.cex = 1.55)

legend("topright",inset = c(0.05, 0.13),legend = names(corr_col),
       lty = 1,col = corr_col,lwd = 5,bty = "n",
       cex = 1.55,title = "Correlation",title.cex = 1.55)

