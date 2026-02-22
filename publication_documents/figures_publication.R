library(dplyr)
library(tidyr)
library(ggplot2)
# Fig 1

pdf(file="Lakens_figure_1.pdf",width=8,height=5)

## Figure illustrating critical effect size

low_x <- -1
high_x <- 2
y_max <- 2

# Add Type 1 error rate function
add_type1_error <- function(N,
                            side = "right",
                            ncp,
                            col = "#00978d") {
  mult <- ifelse(side == "right", 1, -1)
  crit_d <- mult * abs(qt(0.05 / 2, (N * 2) - 2)) / sqrt(N / 2)
  
  if (side == "right") {
    y <- seq(crit_d, 10, length = 10000)
  } else {
    y <- seq(-10, crit_d, length = 10000)
  }
  
  # determine upperbounds polygon
  suppressWarnings({
    z <- (dt(y * sqrt(N / 2), df = (N * 2) - 2, ncp = ncp) * sqrt(N / 2))
  })
  
  if (side == "right") {
    polygon(c(crit_d, y, 10), c(0, z, 0), col = col)
  } else {
    polygon(c(y, crit_d, crit_d), c(z, 0, 0), col = col)
  }
}

# calculate distribution of d based on t-distribution
calc_d_dist <- function(x, N, ncp = 0) {
  suppressWarnings({
    # generates a lot of warnings sometimes
    dt(x * sqrt(N / 2), df = (N * 2) - 2, ncp = ncp) * sqrt(N / 2)
  })
}


# Set sample size per group and effect size d (assumes equal sample sizes per group)
N <- 40 # sample size per group for independent t-test
d <- 0.5 # please enter positive d only
# Calculate non-centrality parameter - equals t-value from sample
ncp <- d * sqrt(N / 2)

# # or Cumming, page 305
# ncp <- d / (sqrt((1 / N) + (1 / N)))

# calc d-distribution
x <- seq(low_x, high_x, length = 10000) # create x values
d_dist <- calc_d_dist(x, N, ncp)

# Set max Y
y_max <- max(d_dist) + 0.5

crit_d <- abs(qt(0.05 / 2, (N * 2) - 2)) / sqrt(N / 2)

# Create base plot with no x-axis drawn
plot(-10, xlim = c(low_x, high_x), ylim = c(0, y_max),
     xlab = expression("Cohen's " * italic(d)),
     ylab = "",
     main = bquote(
       paste("Distribution of observed effects sizes, N = 40, ",
             italic(d), " = ", .(d), ", ",
             d[crit], " = ", .(round(crit_d, 2)))
     ),
     xaxt = "n",   # disable default x-axis
     yaxt = "n"   # disable default x-axis
)

# Add major ticks at each 0.5
axis(side = 1, at = seq(low_x, high_x, by = 0.5), labels = TRUE, lwd = 1.5)

# Add minor ticks at each 0.1 (no labels)
axis(side = 1, at = seq(low_x, high_x, by = 0.1), labels = FALSE, tcl = -0.25)  # tcl controls tick length

d_dist <- dt(x * sqrt(N / 2), df = (N * 2) - 2, ncp = ncp) * sqrt(N / 2)
lines(x, d_dist, col = "black", type = "l", lwd = 2)

# Add type 1 error rate
add_type1_error(N, "right", ncp)
add_type1_error(N, "left", ncp)

abline(v=d, lty = 2)

dev.off()


pdf(file="Lakens_figure_2.pdf",width=8,height=5)

# # Parameters
# n_sim <- 100000
# show_progress <- TRUE
# 
# # Initialize progress bar if needed
# if (show_progress) {
#   pb <- progress_bar$new(
#     format = "  Simulating [:bar] :percent in :elapsed",
#     total = n_sim, clear = FALSE, width = 60
#   )
# }

# # Run simulations
# sim_results <- lapply(seq_len(n_sim), function(i) {
#   if (show_progress) pb$tick()
#   compare_effect_size_corrections(
#     n = 40,
#     mu1 = 1,
#     mu2 = 0,
#     sd = 2,
#     sided = 1,
#     exact = FALSE,
#     retrodesign_d = NULL
#   )
# })

# saveRDS(sim_results, "n=40d=0.5sim=100000.rds")
# to save time we stored the simulation results, and read them back in.
sim_results <- readRDS("../n=40d=0.5sim=100000.rds")
# Combine into a single data frame
final_sim_table <- bind_rows(sim_results)

# Store dataframe for only significant results, so extreme bias based on p < .05
significant_results <- final_sim_table %>%
  filter(p_value < 0.05)

# Below we first have the code for the plot without bias. Not used in the manuscript
# Reshape the data to long format
long_data <- final_sim_table %>%
  select(observed_d, puniform, type_m) %>%
  pivot_longer(cols = everything(), names_to = "variable", values_to = "value")
long_data <- long_data %>%
  mutate(value = ifelse(variable == "puniform" & value < 0, NA, value)) %>%
  mutate(variable_label = recode(variable,
                                 "puniform" = "p-uniform",
                                 "type_m" = "Type M",
                                 "observed_d" = "Observed d"))

long_data$variable_label <- factor(long_data$variable_label,
                                   levels = c("Observed d", "Type M", "p-uniform"))

ggplot(long_data, aes(x = value, fill = variable_label, color = variable_label)) +
  
  # Histogram for Type M and P-uniform
  geom_histogram(
    data = dplyr::filter(long_data, variable_label != "Observed d"),
    aes(y = after_stat(density)),
    position = "identity", alpha = 0.5,
    binwidth = 0.05, boundary = -0.5
  ) +
  
  # Histogram for Observed d, conditional coloring
  geom_histogram(
    data = dplyr::filter(long_data, variable_label == "Observed d"),
    aes(
      y = after_stat(density),
      fill = after_stat(ifelse(x < 0.445, "Observed d (grey)", "Observed d")),
      color = after_stat(ifelse(x < 0.445, "Observed d (grey)", "Observed d"))
    ),
    position = "identity", alpha = 0.5,
    binwidth = 0.05, boundary = -0.5,
    show.legend = FALSE
  ) +
  
  facet_wrap(~ variable_label, scales = "free") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  xlim(-0.5, 1.5) +
  ylim(0, 1.8) +
  labs(
    
    x = "Cohen's d effect size",
    y = NULL
  ) +
  scale_fill_manual(
    values = c(
      "p-uniform" = "#bb7125",
      "Type M" = "#ffdd00",
      "Observed d" = "#00978d",
      "Observed d (grey)" = "grey70"
    ),
    breaks = c("Type M", "p-uniform", "Observed d")
  ) +
  scale_color_manual(
    values = c(
      "p-uniform" = "#bb7125",
      "Type M" = "#ffdd00",
      "Observed d" = "#00978d",
      "Observed d (grey)" = "grey70"
    ),
    breaks = c("Type M", "p-uniform", "Observed d")
  )

dev.off()
