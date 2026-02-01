# Load required packages
library(haven)
library(fitdistrplus)
library(dplyr)
library(ggplot2)
library(parallel)
library(doParallel)
library(foreach)
library(patchwork)
library(magick)       # For flattening PDF
library(pdftools)     # Required for image_read_pdf

# Set seed for gamma fitting
set.seed(897)

# Define gamma distribution parameters
shape_est <- 50.16049
rate_est <- 6.29346
scale_est <- 1/6.29346

# Simulation parameters
simulate <- 1000
row_duplication_levels <- seq(0, 4, by = 1)
n_control <- 20
n_treatment <- 20

# Function to perform one simulation
run_simulation <- function(sample_id, row_dup, shape_est, scale_est, n_control, n_treatment) {
  sim_data <- data.frame(
    group = c(rep(0, n_control), rep(1, n_treatment)),
    sample = sample_id,
    row_duplication = row_dup
  )
  sim_data$r_a1c <- rgamma(nrow(sim_data), shape = shape_est, scale = scale_est)
  if (row_dup > 0 & sample_id > 0) {
    treatment_idx <- sim_data$group == 1
    n_dup <- row_dup
    if (n_dup > 0 & n_dup <= n_treatment) {
      treatment_data <- sim_data[treatment_idx, ]
      smallest_idx <- which.min(treatment_data$r_a1c)
      smallest_participant <- treatment_data[smallest_idx, ]
      duplicated_participant <- bind_rows(replicate(n_dup, smallest_participant, simplify = FALSE))
      largest_indices <- order(treatment_data$r_a1c, decreasing = TRUE)[1:n_dup]
      keep_indices <- setdiff(1:nrow(treatment_data), largest_indices)
      remaining_participants <- treatment_data[keep_indices, ]
      new_treatment_data <- bind_rows(remaining_participants, duplicated_participant)
      sim_data[treatment_idx, ] <- new_treatment_data[1:sum(treatment_idx), ]
    }
  }
  return(sim_data)
}

# Original dataset (sample=0)
master_data <- data.frame(
  group = c(rep(0, n_control), rep(1, n_treatment)),
  sample = 0,
  row_duplication = 0,
  r_a1c = rgamma(n_control + n_treatment, shape = shape_est, scale = scale_est)
)

# Set up parallel processing
cores <- min(detectCores() - 1, 8)
registerDoParallel(cores)

# Run simulations in parallel
sim_results <- foreach(
  i = 1:(simulate * length(row_duplication_levels)),
  .combine = rbind,
  .packages = c("dplyr")
) %dopar% {
  dup_idx <- ceiling(i / simulate) - 1
  row_dup <- row_duplication_levels[dup_idx + 1]
  sample_id <- i + (dup_idx * simulate)
  run_simulation(sample_id, row_dup, shape_est, scale_est, n_control, n_treatment)
}

# Stop parallel processing
stopImplicitCluster()

# Combine original and simulated data
Master <- rbind(master_data, sim_results)

# Function to calculate summary statistics for a sample
calculate_summary <- function(sample_data) {
  treat_data <- sample_data %>% filter(group == 1)
  control_data <- sample_data %>% filter(group == 0)
  if (nrow(treat_data) == 0 || nrow(control_data) == 0) {
    warning("Empty treatment or control group for sample: ", unique(sample_data$sample))
    return(NULL)
  }
  xT <- mean(treat_data$r_a1c, na.rm = TRUE)
  sT <- sd(treat_data$r_a1c, na.rm = TRUE)
  nT <- sum(!is.na(treat_data$r_a1c))
  xC <- mean(control_data$r_a1c, na.rm = TRUE)
  sC <- sd(control_data$r_a1c, na.rm = TRUE)
  nC <- sum(!is.na(control_data$r_a1c))
  if (xT == 0 || xC == 0) {
    warning("Zero mean in treatment or control group for sample: ", unique(sample_data$sample))
    return(NULL)
  }
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
    row_duplication = unique(sample_data$row_duplication),
    sample = unique(sample_data$sample),
    significance = ifelse(p_value <= 0.05, "Significant", "Not Significant")
  )
}

# Calculate Master_summary using lapply
Master_summary <- do.call(rbind, lapply(split(Master, Master$sample), calculate_summary))

# Ensure row_duplication is a factor with consistent levels
Master_summary <- Master_summary %>%
  mutate(row_duplication = factor(row_duplication, levels = row_duplication_levels))

# Calculate proportion of significant p-values for each row duplication level
prop_significant <- Master_summary %>%
  filter(sample != 0) %>%
  group_by(row_duplication) %>%
  summarise(
    n_significant = sum(p_value <= 0.05, na.rm = TRUE),
    n_total = n(),
    prop_significant = n_significant / n_total
  ) %>%
  mutate(
    prop_significant = round(prop_significant, 3),
    row_duplication = factor(row_duplication, levels = row_duplication_levels)
  )

# Print the proportion of significant p-values
# print("Proportion of significant LnCVRs by row duplication level:")
# print(prop_significant)

# Calculate the means manually for the line plot
mean_summary <- Master_summary %>%
  filter(sample != 0) %>%
  group_by(row_duplication) %>%
  summarise(mean_lnCVR = mean(lnCVR, na.rm = TRUE),
            .groups = 'drop') %>%
  mutate(row_duplication = factor(row_duplication, levels = row_duplication_levels))

# Calculate proportion of significant lnCVR values outside boundaries
prop_3sigma <- Master_summary %>%
  filter(sample != 0) %>%
  group_by(row_duplication) %>%
  summarise(
    n_3sigma = sum(p_value <= 0.05 & (lnCVR < -0.51390 | lnCVR > 0.51390), na.rm = TRUE),
    n_total = n(),
    prop_3sigma = n_3sigma / n_total
  ) %>%
  mutate(
    prop_3sigma = round(prop_3sigma, 3),
    row_duplication = factor(row_duplication, levels = row_duplication_levels)
  )

# Print the proportion of significant lnCVR values outside boundaries
# print("%LnCVRs > 3-sigma & p ≤ 0.05 by Row Duplication Level")
# print(prop_3sigma)

# Create the main plot (p)
p <- ggplot(Master_summary, aes(x = row_duplication, y = lnCVR)) +
  geom_point(aes(color = significance), 
             position = position_jitter(width = 0.2, height = 0), 
             size = 1.5, alpha = 0.7) +
  geom_line(data = mean_summary, 
            aes(x = row_duplication, y = mean_lnCVR, group = 1), 
            color = "black", linewidth = 1) +
  geom_point(data = mean_summary, 
             aes(x = row_duplication, y = mean_lnCVR), 
             color = "black", size = 2) +
  labs(
    title = "Best responders (smallest HbA1c post-treatment) duplicated",
    x = "Number of best responder duplications in treatment group",
    y = "LnCVR",
    color = "Significance"
  ) +
  scale_y_continuous(limits = c(-2.5, 2.5), breaks = seq(-2.5, 2.5, by = 0.5)) +
  scale_x_discrete(breaks = row_duplication_levels) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "green")) +
  theme_minimal(base_family = "sans") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    legend.position = "bottom"
  )

# Visualize the proportion of significant p-values (p_prop)
p_prop <- ggplot(prop_significant, aes(x = row_duplication, y = prop_significant)) +
  geom_line(aes(group = 1), color = "blue", linewidth = 1) +
  geom_point(color = "blue", size = 2) +
  labs(
    title = "%LnCVRs (p ≤ 0.05) by Row Duplication Level",
    x = "Number of best responder duplications in treatment group",
    y = "%LnCVRs (p ≤ 0.05)"
  ) +
  scale_x_discrete(breaks = row_duplication_levels) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  theme_minimal(base_family = "sans") +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 8),
    plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
    plot.margin = margin(t = 0.2, b = 0.3, l = 0.3, r = 0.3, unit = "in")
  )

# Visualize the proportion of significant lnCVR values outside boundaries (p_3sigma)
p_3sigma <- ggplot(prop_3sigma, aes(x = row_duplication, y = prop_3sigma)) +
  geom_line(aes(group = 1), color = "purple", linewidth = 1) +
  geom_point(color = "purple", size = 2) +
  labs(
    title = "%LnCVRs (p ≤ 0.05 & > 3-sigma) by Row Duplication Level",
    x = "Number of best responder duplications in treatment group",
    y = "%LnCVRs (p ≤ 0.05 & > 3-sigma)"
  ) +
  scale_x_discrete(breaks = row_duplication_levels) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  theme_minimal(base_family = "sans") +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 8),
    plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
    plot.margin = margin(t = 0.2, b = 0.3, l = 0.3, r = 0.3, unit = "in")
  )

# Function to flatten PDF with error handling
flatten_pdf <- function(input_path, output_path) {
  tryCatch({
    image_write(image_read_pdf(input_path, density = 600), 
                path = output_path, 
                format = "pdf", 
                density = 600)
    cat("Flattened PDF saved to:", output_path, "\n")
  }, error = function(e) {
    cat("Error flattening PDF", input_path, ":", e$message, "\n")
    cat("Saving original PDF without flattening.\n")
  })
}


# Combine the three plots into a single PDF with proper dimensions
combined_plot <- p / p_prop / p_3sigma +
  plot_layout(
    heights = c(2, 2, 2),  # Each plot 2 inches high
    widths = c(7),         # Each plot 7 inches wide
    ncol = 1               # Vertical stacking
  ) +
  plot_annotation(
    title = "Impact of Row Duplication on LnCVR: Small Sample (20/trial arm)",
    subtitle = "Impact on LnCVR, % significant LnCVR, and % significant LnCVR outliers",
    theme = theme(
      plot.title = element_text(family = "sans", size = 14, hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(family = "sans", size = 10, hjust = 0.5),
      plot.margin = margin(t = 0.5, b = 0.3, l = 0.5, r = 0.5, unit = "in")
    )
  )

# Create plots
# Set PDF path 
# Replace with desired path
current_wd <- getwd()
combined_plot_path <- file.path(current_wd,"/combined_row_dup_plots_small.pdf")
combined_dir <- dirname(combined_plot_path)
if (!dir.exists(combined_dir)) {
  dir.create(combined_dir, recursive = TRUE)
}

# Save with exact 8x11 dimensions
ggsave(combined_plot_path, 
       plot = combined_plot, 
       width = 8, 
       height = 11, 
       units = "in", 
       dpi = 300,
       device = cairo_pdf)

# Flatten the PDF
flatten_pdf(combined_plot_path, paste0(combined_plot_path, "_flattened.pdf"))

# Display the combined plot
# print(combined_plot)