# Load required packages
library(sn)           # For skewed normal distribution
library(dplyr)        # For data manipulation
library(fitdistrplus) # For fitting gamma distribution (optional, can be removed if preferred)
library(doParallel)   # For parallel processing
library(foreach)      # For parallel loops

# Define gamma distribution parameters
shape_est <- 50.16049
rate_est <- 6.29346
scale_est <- 1/6.29346

# Define simulation parameters
subgroup <- seq(0.1, 0.5, by = 0.1) # Subgroup proportions
effect <- seq(0, -2, by = -0.2)     # Effect sizes
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
    # Calculate t-test and p-value
    # z = LnCVR / se_LnCVR,
    # p_value = 2 * (1 - pnorm(abs(z))),
    t = LnCVR / se_LnCVR,
    p_value = 2*pt(abs(t),df= 2*n_per_group -2,lower.tail = F),
    # Welch–Satterthwaite degrees of freedom 
    # note the approximate variance of lnCV = ln(sT/xT) = sT^2 / (nT*xT^2)
    var_T = sT^2 / (nT * xT^2), var_C = sC^2 / (nC * xC^2),
    df_welch = (var_T + var_C)^2 / (var_T^2/(nT - 1) +  var_C^2/(nC - 1)),
    p_value_welch = 2*pt(abs(t),df= df_welch,lower.tail = F),
    # Retain simulation parameters
    subgroup = first(subgroup),
    effect = first(effect),
    .groups = "drop"
  )

# Replace with desired path
current_wd <- getwd()

write.table(as.data.frame(Master_summary %>%
                            filter(sample != 0, !is.na(effect)) %>%
                            group_by(effect, subgroup) %>%
                            summarise(
                              specificity = 100*(1-(sum(p_value_welch <= 0.05 & (LnCVR < -0.5 | LnCVR > 0.5), na.rm = TRUE)/ n())),
                              .groups = "drop"
                            ) %>% arrange(desc(effect)) %>% 
                            mutate(effect=as.character(effect),subgroup = subgroup * 100) %>%
                            tidyr::pivot_wider(names_from = subgroup,values_from=specificity)),row.names=F,sep="\t",quote = F,
            file=paste0(current_wd,"/HTE_large_specificity.txt"))
