# Load required packages
library(sn)           # For skewed normal distribution
library(dplyr)        # For data manipulation
library(fitdistrplus) # For fitting gamma distribution (optional, can be removed if preferred)
library(ggplot2)      # For plotting
library(doParallel)   # For parallel processing
library(foreach)      # For parallel loops
library(grid)         # For plot layout
library(magick)       # For flattening PDF
library(pdftools)     # Required for image_read_pdf

# Define gamma distribution parameters
shape_est <- 50.16049
rate_est <- 6.29346
scale_est <- 1/6.29346

# Define simulation parameters
subgroup <- seq(0.1, 0.5, by = 0.1) # Subgroup proportions
effect <- seq(0, -1.4, by = -0.2)     # Effect sizes
simulate <- 1000                    # Number of simulations per combination
n_per_group <- 250                   # Sample size per group (treatment and control)
n_subgroups <- length(subgroup)
n_effects <- length(effect)
n_sim <- simulate
n_samples <- n_subgroups * n_effects * n_sim # Total simulations (excluding sample=0)

# Function to perform one simulation
run_simulation <- function(sample_id, subgroup, effect, shape_est, scale_est, n_per_group) {
  set.seed(sample_id) # Unique seed for each simulation
  
  # Initialize sim_data with 2 * n_per_group rows (n_per_group for control, n_per_group for treatment)
  sim_data <- data.frame(
    group = c(rep(0, n_per_group), rep(1, n_per_group)),
    sample = sample_id,
    subgroup = subgroup,
    effect = effect,
    r_a1c = rgamma(2 * n_per_group, shape = shape_est, scale = scale_est)
  )
  
  # Handle group=1: select subgroup and apply effect
  treatment_rows <- sim_data$group == 1
  n_select <- round(n_per_group * subgroup) # Number of treatment group participants to select
  if (n_select > 0) {
    treatment_indices <- which(treatment_rows)
    selected_indices <- sample(treatment_indices, n_select) # Random selection
    sim_data$r_a1c[selected_indices] <- sim_data$r_a1c[selected_indices] + effect
  }
  
  # group=0 retains unmodified rgamma output
  return(sim_data)
}

# Prepare original dataset for sample=0
original_data <- data.frame(
  group = c(rep(0, n_per_group), rep(1, n_per_group)),
  sample = 0,
  subgroup = NA,
  effect = NA,
  r_a1c = rgamma(2 * n_per_group, shape = shape_est, scale = scale_est)
)

# Set up parallel processing
cores <- detectCores() - 1 # Use all but one core
registerDoParallel(cores)
set.seed(123) # Set global seed for reproducibility

# Run simulations in parallel
Master <- foreach(sg = seq_along(subgroup), .combine = rbind, .packages = c("dplyr")) %:%
  foreach(ef = seq_along(effect), .combine = rbind) %:%
  foreach(sim = 1:n_sim, .combine = rbind) %dopar% {
    sample_id <- (sg - 1) * n_effects * n_sim + (ef - 1) * n_sim + sim
    run_simulation(sample_id, subgroup[sg], effect[ef], shape_est, scale_est, n_per_group)
  }

# Append original dataset (sample=0) to Master
Master <- rbind(original_data, Master)

# Create Master_summary dataset
Master_summary <- Master %>%
  group_by(sample) %>%
  summarise(
    # Treatment group (group=1)
    xT = mean(r_a1c[group == 1], na.rm = TRUE),
    sT = sd(r_a1c[group == 1], na.rm = TRUE),
    nT = sum(!is.na(r_a1c[group == 1])),
    # Control group (group=0)
    xC = mean(r_a1c[group == 0], na.rm = TRUE),
    sC = sd(r_a1c[group == 0], na.rm = TRUE),
    nC = sum(!is.na(r_a1c[group == 0])),
    # Calculate LnCVR
    CVT = sT / xT,
    CVC = sC / xC,
    LnCVR_naive = log(CVT / CVC),
    correction1 = 0.5 * (1 / (nT - 1) - 1 / (nC - 1)),
    correction2 = 0.5 * ((sC^2 / (nC * xC^2)) - (sT^2 / (nT * xT^2))),
    LnCVR = LnCVR_naive + correction1 + correction2,
    # Calculate variance of LnCVR
    var_LnCVR = (sT^2 / (nT * xT^2)) + (sC^2 / (nC * xC^2)) +
      0.5 * (nC / (nC - 1)^2) + 0.5 * (nT / (nT - 1)^2) +
      0.5 * (sT^4 / (xT^4 * nT^2)) +
      0.5 * (sC^4 / (xC^4 * nC^2)),
    # Calculate standard error
    se_LnCVR = sqrt(var_LnCVR),
    # Calculate 95% confidence interval
    lower_95_ci = LnCVR - 1.96 * se_LnCVR,
    upper_95_ci = LnCVR + 1.96 * se_LnCVR,
    # Calculate z-test and p-value
    z = LnCVR / se_LnCVR,
    p_value = 2 * (1 - pnorm(abs(z))),
    # Retain simulation parameters
    subgroup = first(subgroup),
    effect = first(effect),
    .groups = "drop"
  )

# Define y-axis limits for plots
# For page 1: LnCVR values (widened to accommodate data range)
y_limits_page1 <- c(-2.0, 2.0) # Adjusted to cover potential LnCVR range
y_breaks_page1 <- seq(-2.0, 2.0, by = 0.5)

# For page 2: Proportion of significant p-values
y_limits_page2 <- c(0, 1)
y_breaks_page2 <- seq(0, 1, by = 0.2)

# For page 3: Proportion of significant LnCVR outside 3-sigma
y_limits_page3 <- c(0, 0.50)
y_breaks_page3 <- seq(0, 0.50, by = 0.10)

# Define effect levels globally to avoid duplication
effect_levels <- as.character(effect) # Use global effect vector: [0, -0.4, -0.8, -1.2, -1.6, -2.0]

# Create plots
# Set PDF path
username <- Sys.getenv("USERNAME")
plot_path <- file.path("C:", "Users", username, "Dropbox", "Carlisle", "R-diabetes",
                       "Submission_manuscripts", "LnCVR-SAS", "LnCVR-publication_SAS", "Revision", "Reviewer_request1_Simulate heterogeneity",
                       "heterogeneity_plots_large.pdf")
# Open PDF device with cairo_pdf to avoid font issues
cairo_pdf(plot_path, width = 8, height = 11)

# Page 1: Original plots (LnCVR values)
# Calculate layout parameters
plot_height <- 1.8 # inches
n_plots <- 5
gap_height <- 0.21875 # inches (0.875 / 4)
title_height <- 0.125 # inches
total_content_height <- title_height + n_plots * plot_height + (n_plots - 1) * gap_height
margin_height <- (11 - total_content_height) / 2 # Equal top and bottom margins (~0.5625 inches)

# Set up plot layout with viewport
grid.newpage()
pushViewport(viewport(y = 0.5, height = total_content_height / 11))

# Add title for page 1
grid.text("Impact of HTE on LnCVR in large Samples (250/trial arm): A Simulation Study",
          y = unit(1, "npc") - unit(0.00, "inches"),
          gp = gpar(fontsize = 12, fontface = "bold"))

# Create plots for each subgroup
plots <- lapply(1:n_plots, function(i) {
  sg <- subgroup[i]
  # Filter data for the specific subgroup and exclude effect = NA
  data_plot <- Master_summary %>% filter(subgroup == sg, !is.na(effect))
  
  # Calculate mean LnCVR for each effect level
  mean_data <- data_plot %>%
    group_by(effect) %>%
    summarise(mean_LnCVR = mean(LnCVR, na.rm = TRUE), .groups = "drop")
  
  ggplot(data_plot, aes(x = factor(effect, levels = effect_levels), y = LnCVR)) +
    geom_point(aes(color = p_value <= 0.05, alpha = p_value <= 0.05), size = 1) +
    geom_line(data = mean_data, aes(x = factor(effect, levels = effect_levels), y = mean_LnCVR, group = 1)) +
    scale_color_manual(values = c("TRUE" = "tomato1", "FALSE" = "forestgreen"),
                       labels = c("TRUE" = "p <= 0.05", "FALSE" = "p > 0.05"),
                       name = "Significance") +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.5),
                       labels = c("TRUE" = "p <= 0.05", "FALSE" = "p > 0.05"),
                       name = "Significance") +
    scale_y_continuous(limits = y_limits_page1, breaks = y_breaks_page1) +
    labs(x = sprintf("HbA1c treatment effect in %.0f percent subgroup", sg * 100),
         y = "LnCVR value") +
    theme_minimal(base_size = 12) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      axis.title.x = element_text(size = 8),
      axis.title.y = element_text(size = 8),
      plot.title = element_text(size = 12, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5)
    )
})

# Arrange plots with gaps on page 1
for (i in 1:n_plots) {
  y_pos <- 1 - (title_height + (i - 1) * (plot_height + gap_height) + plot_height/2) / total_content_height
  pushViewport(viewport(y = unit(y_pos, "npc"),
                        height = unit(plot_height, "inches"),
                        width = unit(7, "inches")))
  print(plots[[i]], vp = viewport())
  popViewport()
}

# Page 2: Proportion of significant LnCVR p-values
# Calculate proportion of significant p-values
prop_significant <- Master_summary %>%
  filter(sample != 0, !is.na(effect)) %>%
  group_by(effect, subgroup) %>%
  summarise(
    n_significant = sum(p_value <= 0.05, na.rm = TRUE),
    n_total = n(),
    prop_significant = n_significant / n_total,
    .groups = "drop"
  ) %>%
  mutate(prop_significant = round(prop_significant, 3))

# New page for proportion plots
grid.newpage()
pushViewport(viewport(y = 0.5, height = total_content_height / 11))

# Add title for page 2
grid.text("Impact of HTE on %LnCVR (p<=0.05) in large Samples: A Simulation Study",
          y = unit(1, "npc") - unit(0.00, "inches"),
          gp = gpar(fontsize = 12, fontface = "bold"))

# Create proportion plots for each subgroup
prop_plots <- lapply(1:n_plots, function(i) {
  sg <- subgroup[i]
  data_plot <- prop_significant %>% filter(subgroup == sg)
  
  ggplot(data_plot, aes(x = factor(effect, levels = effect_levels), y = prop_significant)) +
    geom_point(color = "black", size = 1) +
    geom_line(aes(group = 1), color = "black") +
    scale_y_continuous(limits = y_limits_page2, breaks = y_breaks_page2) +
    labs(x = sprintf("HbA1c treatment effect in %.0f percent subgroup", sg * 100),
         y = "%LnCVR values (p<=0.05)") +
    theme_minimal(base_size = 12) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      axis.title.x = element_text(size = 8),
      axis.title.y = element_text(size = 8),
      plot.title = element_text(size = 12, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      legend.position = "none"
    )
})

# Arrange proportion plots with gaps on page 2
for (i in 1:n_plots) {
  y_pos <- 1 - (title_height + (i - 1) * (plot_height + gap_height) + plot_height/2) / total_content_height
  pushViewport(viewport(y = unit(y_pos, "npc"),
                        height = unit(plot_height, "inches"),
                        width = unit(7, "inches")))
  print(prop_plots[[i]], vp = viewport())
  popViewport()
}

# Page 3: Proportion of significant LnCVR values outside boundaries
# Calculate proportion of significant LnCVR values outside boundaries
prop_3sigma <- Master_summary %>%
  filter(sample != 0, !is.na(effect)) %>%
  group_by(effect, subgroup) %>%
  summarise(
    n_3sigma = sum(p_value <= 0.05 & (LnCVR < -0.51390 | LnCVR > 0.51390), na.rm = TRUE),
    n_total = n(),
    prop_3sigma = n_3sigma / n_total,
    .groups = "drop"
  ) %>%
  mutate(prop_3sigma = round(prop_3sigma, 3))

# New page for 3-sigma proportion plots
grid.newpage()
pushViewport(viewport(y = 0.5, height = total_content_height / 11))

# Add title for page 3
grid.text("Impact of HTE on %LnCVR (p<=0.05 &> 3-sigma) in large Samples: A Simulation Study",
          y = unit(1, "npc") - unit(0.00, "inches"),
          gp = gpar(fontsize = 12, fontface = "bold"))

# Create 3-sigma proportion plots for each subgroup
sigma_plots <- lapply(1:n_plots, function(i) {
  sg <- subgroup[i]
  data_plot <- prop_3sigma %>% filter(subgroup == sg)
  
  ggplot(data_plot, aes(x = factor(effect, levels = effect_levels), y = prop_3sigma)) +
    geom_point(color = "black", size = 1) +
    geom_line(aes(group = 1), color = "black") +
    scale_y_continuous(limits = y_limits_page3, breaks = y_breaks_page3) +
    labs(x = sprintf("HbA1c treatment effect in %.0f percent subgroup", sg * 100),
         y = "%LnCVR values (p<=0.05 & > 3-sigma)") +
    theme_minimal(base_size = 12) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      axis.title.x = element_text(size = 8),
      axis.title.y = element_text(size = 8),
      plot.title = element_text(size = 12, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      legend.position = "none"
    )
})

# Arrange 3-sigma proportion plots with gaps on page 3
for (i in 1:n_plots) {
  y_pos <- 1 - (title_height + (i - 1) * (plot_height + gap_height) + plot_height/2) / total_content_height
  pushViewport(viewport(y = unit(y_pos, "npc"),
                        height = unit(plot_height, "inches"),
                        width = unit(7, "inches")))
  print(sigma_plots[[i]], vp = viewport())
  popViewport()
}

# Close PDF device
dev.off()

# Flatten the PDF using magick
image_write(image_read_pdf(plot_path, density = 600), path = paste0(plot_path, "_flattened.pdf"), format = "pdf", density = 600)

# Stop parallel processing
stopImplicitCluster()

# Calculate mean of LnCVR stratified by subgroup and effect in Master_summary
stratified_mean <- Master_summary %>%
  group_by(subgroup, effect) %>%
  summarise(mean_LnCVR = mean(LnCVR, na.rm = TRUE), .groups = "drop")

# Print the result
print(stratified_mean)
View(stratified_mean)

# Calculate proportion of significant p-values
result <- Master_summary %>%
  group_by(subgroup, effect) %>%
  summarise(
    total_count = n(),
    count_p_0_05 = sum(p_value >= 0 & p_value <= 0.05, na.rm = TRUE),
    proportion = count_p_0_05 / total_count,
    .groups = "drop"
  )

# View the result
print(result, n = 60)