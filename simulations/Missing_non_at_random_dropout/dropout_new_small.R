# Load required packages
library(haven)
library(fitdistrplus)
library(dplyr)
library(ggplot2)
library(parallel)
library(doParallel)
library(foreach)
library(patchwork)

# Set seed for gamma fitting
set.seed(897)

# Define gamma distribution parameters
shape_est <- 50.16049
rate_est <- 6.29346
scale_est <- 1/6.29346

# Simulation parameters
simulate <- 1000
dropout_levels <- seq(0, 0.5, by = 0.1)
n_control <- 20
n_treatment <- 20

# Function to perform one simulation
run_simulation <- function(sample_id, dropout, shape_est, scale_est, n_control, n_treatment) {
  set.seed(sample_id) # Unique seed for each simulation
  
  # Initialize sim_data with n_control + n_treatment rows
  sim_data <- data.frame(
    group = c(rep(0, n_control), rep(1, n_treatment)),
    sample = sample_id,
    dropout = dropout,
    r_a1c = rgamma(n_control + n_treatment, shape = shape_est, scale = scale_est)
  )
  
  if (dropout > 0 & sample_id > 0) {
    control_idx <- sim_data$group == 0
    n_drop <- round(dropout * n_control)
    if (n_drop > 0) {
      control_data <- sim_data[control_idx, ]
      control_order <- order(control_data$r_a1c, decreasing = TRUE, 
                             runif(nrow(control_data)))
      drop_idx <- control_order[1:n_drop]
      control_drop_idx <- which(control_idx)[drop_idx]
      sim_data <- sim_data[-control_drop_idx, ]
    }
  }
  return(sim_data)
}

# Original dataset (sample=0)
master_data <- data.frame(
  group = c(rep(0, n_control), rep(1, n_treatment)),
  sample = 0,
  dropout = 0,
  r_a1c = rgamma(n_control + n_treatment, shape = shape_est, scale = scale_est)
)

# Set up parallel processing
cores <- min(detectCores() - 1, 8)
registerDoParallel(cores)

# Run simulations in parallel
sim_results <- foreach(
  i = 1:(simulate * length(dropout_levels)),
  .combine = rbind,
  .packages = c("dplyr")
) %dopar% {
  dropout_idx <- ceiling(i / simulate) - 1
  dropout <- dropout_levels[dropout_idx + 1]
  sample_id <- i + (dropout_idx * simulate)
  run_simulation(sample_id, dropout, shape_est, scale_est, n_control, n_treatment)
}

# Stop parallel processing
stopImplicitCluster()

# Combine original and simulated data
Master <- rbind(master_data, sim_results)

# Function to calculate summary statistics for a sample
calculate_summary <- function(sample_data) {
  treat_data <- sample_data %>% filter(group == 1)
  xT <- mean(treat_data$r_a1c, na.rm = TRUE)
  sT <- sd(treat_data$r_a1c, na.rm = TRUE)
  nT <- sum(!is.na(treat_data$r_a1c))
  
  control_data <- sample_data %>% filter(group == 0)
  xC <- mean(control_data$r_a1c, na.rm = TRUE)
  sC <- sd(control_data$r_a1c, na.rm = TRUE)
  nC <- sum(!is.na(control_data$r_a1c))
  
  CVT <- sT / xT
  CVC <- sC / xC
  lnCVR_naive <- log(CVT / CVC)
  correction1 <- 0.5 * (1 / (nT - 1) - 1 / (nC - 1))
  correction2 <- 0.5 * ((sC^2 / (nC * xC^2)) - (sT^2 / (nT * xT^2)))
  lnCVR <- lnCVR_naive + correction1 + correction2
  
  var_lnCVR <- (sT^2 / (nT * xT^2)) + (sC^2 / (nC * xC^2)) +
    0.5 * (nC / (nC - 1)^2) + 0.5 * (nT / (nT - 1)^2) +
    0.5 * (sT^4 / (xT^4 * nT^2)) +
    0.5 * (sC^4 / (xC^4 * nC^2))
  
  se_LnCVR <- sqrt(var_lnCVR)
  lower_95_ci <- lnCVR - 1.96 * se_LnCVR
  upper_95_ci <- lnCVR + 1.96 * se_LnCVR
  z <- lnCVR / se_LnCVR
  p_value <- 2 * (1 - pnorm(abs(z)))
  
  data.frame(
    xT = xT, sT = sT, nT = nT,
    xC = xC, sC = sC, nC = nC,
    lnCVR_naive = lnCVR_naive, lnCVR = lnCVR,
    var_lnCVR = var_lnCVR, se_LnCVR = se_LnCVR,
    lower_95_ci = lower_95_ci, upper_95_ci = upper_95_ci,
    z = z, p_value = p_value,
    dropout = unique(sample_data$dropout),
    sample = unique(sample_data$sample),
    significance = ifelse(p_value <= 0.05, "Significant", "Not Significant")
  )
}

# Calculate Master_summary using lapply
Master_summary <- do.call(rbind, lapply(split(Master, Master$sample), calculate_summary))

# Set PDF path
# Replace with desired path
current_wd <- getwd()

# for each setting: quantiles, proportion significant, and proportion of significant outliers  
write.table(as.data.frame(Master_summary %>%
                            filter(sample != 0) %>%
                            group_by(dropout) %>%
                            summarise(
                              range = paste0(round(quantile(lnCVR),digits=2),collapse = ";"),
                              prop_sig = sum(p_value <= 0.05, na.rm = TRUE) / n()*100,
                              prop_3sigma = sum(p_value <= 0.05 & (lnCVR < -0.51390 | lnCVR > 0.51390), na.rm = TRUE) *100/ n(),
                              .groups = "drop"
                            ) %>%
                            mutate(prop_3sigma = round(prop_3sigma, 3))),row.names=F,sep="\t",quote = F,
            file=paste0(current_wd,"/dropout_small.txt"))


# Calculate proportion of significant p-values for each dropout level
prop_significant <- Master_summary %>%
  filter(sample != 0) %>%
  group_by(dropout) %>%
  summarise(
    n_significant = sum(p_value <= 0.05, na.rm = TRUE),
    n_total = n(),
    prop_significant = n_significant / n_total
  ) %>%
  mutate(prop_significant = round(prop_significant, 3))

# Print the proportion of significant p-values
# print("Proportion of significant LnCVRs by dropout level:")
# print(prop_significant)

# Visualize the proportion of significant p-values
p_prop <- ggplot(prop_significant, aes(x = dropout, y = prop_significant)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(color = "blue", size = 2) +
  labs(
    title = "%LnCVRs (p <= 0.05) by Dropout Level",
    x = "Fraction poorest responders dropping out of control group",
    y = "%LnCVRs (p<= 0.05)"
  ) +
  scale_x_continuous(breaks = dropout_levels) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  theme_minimal(base_size = 12) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    plot.title = element_text(size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5)
  )

# Save the proportion plot
prop_plot_path <- file.path(current_wd,"prop_significant_dropout_small.pdf")
prop_dir <- dirname(prop_plot_path)
if (!dir.exists(prop_dir)) {
  dir.create(prop_dir, recursive = TRUE)
}
ggsave(prop_plot_path, plot = p_prop, width = 8, height = 6, units = "in", dpi = 300)

# Display the proportion plot
# print(p_prop)

# Calculate proportion of significant lnCVR values outside boundaries
prop_3sigma <- Master_summary %>%
  filter(sample != 0) %>%
  group_by(dropout) %>%
  summarise(
    n_3sigma = sum(p_value <= 0.05 & (lnCVR < -0.51390 | lnCVR > 0.51390), na.rm = TRUE),
    n_total = n(),
    prop_3sigma = n_3sigma / n_total
  ) %>%
  mutate(prop_3sigma = round(prop_3sigma, 3))

# Print the proportion of significant lnCVR values outside boundaries
# print("%LnCVRs > 3-sigma & p <= 0.05 by Dropout Level")
# print(prop_3sigma)

# Visualize the proportion of significant lnCVR values outside boundaries
p_3sigma <- ggplot(prop_3sigma, aes(x = dropout, y = prop_3sigma)) +
  geom_line(color = "purple", linewidth = 1) +
  geom_point(color = "purple", size = 2) +
  labs(
    title = "%LnCVRs (p <= 0.05 & > 3-sigma) by Dropout Level",
    x = "Fraction poorest responders dropping out of control group",
    y = "%LnCVRs (p <= 0.05 & > 3-sigma)"
  ) +
  scale_x_continuous(breaks = dropout_levels) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  theme_minimal(base_size = 12) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    plot.title = element_text(size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5)
  )

# Save the 3-sigma proportion plot
sigma_plot_path <- file.path(current_wd,"3_sigma_dropout_small.pdf")
sigma_dir <- dirname(sigma_plot_path)
if (!dir.exists(sigma_dir)) {
  dir.create(sigma_dir, recursive = TRUE)
}
ggsave(sigma_plot_path, plot = p_3sigma, width = 8, height = 6, units = "in", dpi = 300)

# Display the 3-sigma proportion plot
# print(p_3sigma)

# Calculate mean lnCVR for each dropout level
mean_lnCVR <- Master_summary %>%
  filter(sample != 0) %>%
  group_by(dropout) %>%
  summarise(mean_lnCVR = mean(lnCVR, na.rm = TRUE))

# Create the main plot
p <- ggplot(Master_summary %>% filter(sample != 0), aes(x = dropout, y = lnCVR)) +
  geom_point(aes(color = significance)) +
  geom_line(data = mean_lnCVR, aes(y = mean_lnCVR), color = "black", linewidth = 1) +
  scale_color_manual(name = "Significance", values = c("Significant" = "red", "Not Significant" = "green")) +
  labs(
    title = "LnCVR Values by Dropout Level",
    x = "Fraction poorest responders dropping out of control group",
    y = "LnCVR Value"
  ) +
  scale_y_continuous(limits = c(-2.5, 2.5), breaks = seq(-2.5, 2.5, by = 0.5)) +
  scale_x_continuous(breaks = dropout_levels) +
  theme_minimal(base_size = 12) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    plot.title = element_text(size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5)
  )

# Save the main plot
plot_path <- file.path(current_wd,"dropout_small.pdf")
plot_dir <- dirname(plot_path)
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}
ggsave(plot_path, plot = p, width = 8, height = 11, units = "in", dpi = 300)

# Display the main plot
# print(p)

# Combine the three plots into a single PDF with a title
combined_plot <- p / p_prop / p_3sigma +
  plot_layout(heights = c(1, 1, 1), ncol = 1, guides = "collect") +
  plot_annotation(
    title = "Impact of MNAR on LnCVR in Small Samples (20/trial arm): A Simulation Study",
    theme = theme(
      plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      plot.margin = margin(t = 0.5, b = 0.5, l = 0.5, r = 0.5, unit = "in")
    )
  )

# Save the combined plot
combined_plot_path <- file.path(current_wd,"combined_plots_small.pdf")
combined_dir <- dirname(combined_plot_path)
if (!dir.exists(combined_dir)) {
  dir.create(combined_dir, recursive = TRUE)
}
ggsave(combined_plot_path, plot = combined_plot, width = 8, height = 11, units = "in", dpi = 300)

# Display the combined plot
# print(combined_plot)