fit <- readRDS("fit_real_lg_0.3000_lt_0.0100.rds")
G <- fit$Gamma

sp10 <- c("F.prau", "B.vulg", "F.plau", "B.unif", "B.ovat",
          "E.rect", "P.dist", "A.hadr", "F.sacc", "B.thet")
rownames(G) <- sp10; colnames(G) <- sp10

g_diag <- diag(G)
G0 <- G - diag(g_diag)
cutoff <- 0.001
p <- length(sp10)

firm <- c("F.prau", "F.plau", "F.sacc", "A.hadr", "E.rect")
phy <- ifelse(sp10 %in% firm, "Firmicutes", "Bacteroidetes")
phy_col <- c(Firmicutes = "#E67E22", Bacteroidetes = "#3498DB")
node_col <- phy_col[phy]

max_off  <- max(abs(G0))
max_diag <- max(abs(g_diag))
rad <- rep(0.20, p)


# --- layout ---
W <- (abs(G0) + t(abs(G0))) / 2
W[W < cutoff] <- 0
if (max(W) > 0) W <- W / max(W)

set.seed(7)
m <- p
xy <- matrix(runif(2*m, -0.5, 0.5), ncol = 2)
xy[phy == "Firmicutes",    2] <- xy[phy == "Firmicutes",    2] + 0.6
xy[phy == "Bacteroidetes", 2] <- xy[phy == "Bacteroidetes", 2] - 0.6

k <- sqrt(4/m); temp <- 0.25
for (it in 1:400) {
  move <- matrix(0, m, 2)
  for (i in 1:m) for (j in 1:m) {
    if (i == j) next
    v <- xy[i,] - xy[j,]
    d <- max(sqrt(sum(v^2)), 1e-6)
    move[i,] <- move[i,] + v/d * k^2/d
  }
  for (i in 1:(m-1)) for (j in (i+1):m) {
    v <- xy[i,] - xy[j,]
    d <- max(sqrt(sum(v^2)), 1e-6)
    a <- W[i,j]
    if (phy[i] == phy[j]) a <- a + 0.18
    if (a == 0) next
    f <- d^2/k * a
    move[i,] <- move[i,] - v/d * f
    move[j,] <- move[j,] + v/d * f
  }
  for (i in 1:m) {
    dd <- max(sqrt(sum(move[i,]^2)), 1e-6)
    xy[i,] <- xy[i,] + move[i,]/dd * min(dd, temp)
  }
  temp <- temp * 0.96
}
xy[,1] <- (xy[,1] - mean(xy[,1])) / max(abs(xy[,1] - mean(xy[,1]))) * 0.9
xy[,2] <- (xy[,2] - mean(xy[,2])) / max(abs(xy[,2] - mean(xy[,2]))) * 0.9


# --- push nodes apart ---
for (it in 1:1000) {
  biggest <- 0
  for (i in 1:(m-1)) for (j in (i+1):m) {
    v <- xy[i,] - xy[j,]
    d <- sqrt(sum(v^2))
    if (d >= 0.46) next
    if (d < 1e-8) {
      ang <- 2*pi*(i+j)/(2*m)
      u <- c(cos(ang), sin(ang)); d <- 0
    } else { u <- v/d }
    gap <- 0.46 - d
    xy[i,] <- xy[i,] + u*(gap/2 + 1e-4)
    xy[j,] <- xy[j,] - u*(gap/2 + 1e-4)
    biggest <- max(biggest, gap)
  }
  xy[,1] <- pmin(1.02, pmax(-1.02, xy[,1]))
  xy[,2] <- pmin(1.02, pmax(-1.02, xy[,2]))
  if (biggest < 1e-4) break
}


# --- nudge specific pairs ---
for (pair in list(c("F.prau","B.ovat"), c("E.rect","F.prau"))) {
  i <- match(pair[1], sp10); j <- match(pair[2], sp10)
  v <- xy[i,] - xy[j,]; d <- sqrt(sum(v^2))
  if (d < 0.68) {
    u <- if (d < 1e-8) c(1,0) else v/d
    xy[i,] <- xy[i,] + u*(0.68-d)/2
    xy[j,] <- xy[j,] - u*(0.68-d)/2
  }
}

# push nodes apart again after nudging
for (it in 1:1000) {
  biggest <- 0
  for (i in 1:(m-1)) for (j in (i+1):m) {
    v <- xy[i,] - xy[j,]
    d <- sqrt(sum(v^2))
    if (d >= 0.46) next
    if (d < 1e-8) {
      ang <- 2*pi*(i+j)/(2*m)
      u <- c(cos(ang), sin(ang)); d <- 0
    } else { u <- v/d }
    gap <- 0.46 - d
    xy[i,] <- xy[i,] + u*(gap/2 + 1e-4)
    xy[j,] <- xy[j,] - u*(gap/2 + 1e-4)
    biggest <- max(biggest, gap)
  }
  xy[,1] <- pmin(1.02, pmax(-1.02, xy[,1]))
  xy[,2] <- pmin(1.02, pmax(-1.02, xy[,2]))
  if (biggest < 1e-4) break
}

x <- xy[,1]; y <- xy[,2]
cx <- mean(x); cy <- mean(y)


# --- plot ---
par(mar = c(2,2,2,2), bg = "white")
plot.new()
plot.window(xlim = c(-1.35,1.35), ylim = c(-1.25,1.25), asp = 1)

# directed edges
for (r in 1:p) for (c in 1:p) {
  if (r == c) next
  val <- G0[r, c]
  if (abs(val) <= cutoff) next
  
  dx <- x[r]-x[c]; dy <- y[r]-y[c]
  dd <- sqrt(dx^2+dy^2)
  ux <- dx/dd; uy <- dy/dd
  xs <- x[c]+ux*rad[c]; ys <- y[c]+uy*rad[c]
  xe <- x[r]-ux*rad[r]; ye <- y[r]-uy*rad[r]
  
  this_lwd <- pmax(1.2, abs(val)/max_off * 6)
  this_col <- ifelse(val > 0, "#2C3E50", "#C0392B")
  
  if (abs(G0[c,r]) > cutoff) {
    off <- 0.07
    mx <- (xs+xe)/2 - uy*off
    my <- (ys+ye)/2 + ux*off
    xspline(c(xs,mx,xe), c(ys,my,ye), shape=1, open=TRUE,
            border=this_col, lwd=this_lwd)
    tx <- xe-mx; ty <- ye-my; td <- sqrt(tx^2+ty^2)
    arrows(xe-tx/td*0.04, ye-ty/td*0.04, xe, ye,
           length=0.11, angle=25, lwd=this_lwd, col=this_col)
  } else {
    arrows(xs, ys, xe, ye, length=0.11, angle=25,
           lwd=this_lwd, col=this_col)
  }
}

# self-loops
for (i in 1:p) {
  val <- g_diag[i]
  if (abs(val) <= cutoff) next
  
  dx <- if (sp10[i] == "F.prau") 1 else x[i]-cx
  dy <- if (sp10[i] == "F.prau") 1 else y[i]-cy
  dd <- sqrt(dx^2+dy^2)
  if (dd < 1e-6) { dx <- 1; dy <- 0; dd <- 1 }
  ux <- dx/dd; uy <- dy/dd
  
  rr   <- 0.07
  lx   <- x[i] + ux*(rad[i] + rr*0.45)
  ly   <- y[i] + uy*(rad[i] + rr*0.45)
  base <- atan2(-uy, -ux); gap <- 0.55
  th   <- seq(base+gap, base+2*pi-gap, length.out = 80)
  this_lwd <- pmax(0.8, abs(val)/max_diag * 1.8)
  this_col <- ifelse(val > 0, "#2C3E50", "#E67E22")
  lines(lx+rr*cos(th), ly+rr*sin(th), lwd=this_lwd, col=this_col, lend="round")
  
  th2 <- base+2*pi-gap
  xe2 <- lx+rr*cos(th2); ye2 <- ly+rr*sin(th2)
  tx  <- -sin(th2);       ty  <-  cos(th2)
  arrows(xe2-tx*0.02, ye2-ty*0.02, xe2, ye2,
         length=0.07, angle=25, lwd=this_lwd, col=this_col)
}

# nodes + labels
symbols(x, y, circles=rad, add=TRUE,
        bg=node_col, fg="white", inches=FALSE, lwd=2)
text(x, y, sp10, cex=1.43, font=2)

legend("topleft", inset=c(0.03,0.13),
       legend=names(phy_col), pch=21,
       pt.bg=phy_col, col="black",
       pt.cex=3.1, bty="n", cex=1.55,
       title="Phylum", title.cex=1.55)