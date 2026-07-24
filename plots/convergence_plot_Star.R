library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

# Load data
bin <- readRDS("plots/res/binary/(N,T,p)-(100,10,10)-Star.rds")
poi <- readRDS("plots/res/poisson/(N,T,p)-(100,10,10)-Star.rds")
mul <- readRDS("plots/res/multinomial/(N,T,p)-(100,10,10)-Star.rds")

# Build data frame
n_iter <- 100

df <- data.frame(
  Iteration = rep(1:n_iter, 6),
  Outcome = rep(rep(c("Binary", "Poisson", "Multinomial"), each = n_iter), 2),
  Matrix = rep(c("Gamma", "Theta"), each = n_iter * 3),
  Change = c(
    # Gamma
    bin$convergence$gamma,
    poi$conv_gamma,
    mul$conv_gamma,
    # Theta
    bin$convergence$theta,
    poi$conv_theta,
    mul$conv_theta
  )
)

df$Outcome <- factor(df$Outcome, levels = c("Binary", "Poisson", "Multinomial"))
df$Matrix <- factor(df$Matrix, levels = c("Gamma", "Theta"))

positive_change <- df$Change[is.finite(df$Change) & df$Change > 0]
y_limits <- range(positive_change, na.rm = TRUE)

# Panel A: Gamma
pA <- ggplot(df %>% filter(Matrix == "Gamma"),
             aes(x = Iteration, y = Change, color = Outcome)) +
  geom_line(linewidth = 3) +
  scale_y_log10(limits = y_limits) +
  scale_color_manual(values = c("Binary" = "#1f77b4", "Poisson" = "#ff7f0e", "Multinomial" = "#2ca02c")) +
  labs(
    title = expression(bold("convergence of " * Gamma)),
    x = "Iteration",
    y = "Parameter Change (log scale)",
    color = NULL
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 46, hjust = 0.5),
    axis.title = element_text(size = 46),
    axis.text = element_text(size = 44),
    legend.title = element_text(size = 44),
    legend.text = element_text(size = 42),
    legend.key.size = unit(3.6, "cm"),
    legend.key.width = unit(4.8, "cm"),
    legend.position = "bottom"
  )

# Panel B: Theta
pB <- ggplot(df %>% filter(Matrix == "Theta"),
             aes(x = Iteration, y = Change, color = Outcome)) +
  geom_line(linewidth = 3) +
  scale_y_log10(limits = y_limits) +
  scale_color_manual(values = c("Binary" = "#1f77b4", "Poisson" = "#ff7f0e", "Multinomial" = "#2ca02c")) +
  labs(
    title = expression(bold("convergence of " * Theta)),
    x = "Iteration",
    y = NULL,
    color = NULL
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 46, hjust = 0.5),
    axis.title = element_text(size = 46),
    axis.text = element_text(size = 44),
    legend.title = element_text(size = 44),
    legend.text = element_text(size = 42),
    legend.key.size = unit(3.6, "cm"),
    legend.key.width = unit(4.8, "cm"),
    legend.position = "bottom"
  )

# Combine with shared legend
combined <- plot_grid(
  pA + theme(legend.position = "none"),
  pB + theme(legend.position = "none"),
  ncol = 2, labels = NULL
)

legend <- get_legend(pA + theme(legend.position = "bottom"))
final <- plot_grid(combined, legend, ncol = 1, rel_heights = c(1, 0.13))
