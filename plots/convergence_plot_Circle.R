library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

# Load data
bin <- readRDS("plots/res/binary/(N,T,p)-(100,10,10)-Circle.rds")
poi <- readRDS("plots/res/poisson/(N,T,p)-(100,10,10)-Circle.rds")
mul <- readRDS("plots/res/multinomial/(N,T,p)-(100,10,10)-Circle.rds")


# Build data frame — handle different vector lengths
make_df <- function(gamma_vec, theta_vec, outcome_label) {
  n <- length(gamma_vec)
  rbind(
    data.frame(Iteration = 1:n, Outcome = outcome_label, Matrix = "Gamma", Change = gamma_vec),
    data.frame(Iteration = 1:length(theta_vec), Outcome = outcome_label, Matrix = "Theta", Change = theta_vec)
  )
}

df <- rbind(
  make_df(bin$convergence$gamma, bin$convergence$theta, "Binary"),
  make_df(poi$conv_gamma, poi$conv_theta, "Poisson"),
  make_df(mul$conv_gamma, mul$conv_theta, "Multinomial")
)

# Remove zero
df <- df %>% filter(!is.na(Change) & Change > 0)

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
