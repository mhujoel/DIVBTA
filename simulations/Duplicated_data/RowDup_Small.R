# Load required packages
library(haven)
library(fitdistrplus)
library(dplyr)
library(ggplot2)
library(parallel)
library(doParallel)
library(foreach)

# Set seed for gamma fitting
set.seed(897)

# Define gamma distribution parameters
shape_est <- 50.16049
rate_est <- 6.29346
scale_est <- 1/6.29346

# Simulation parameters
simulate <- 1000
n_control <- 20
n_treatment <- 20
row_duplication_levels <- c(10, 15,16, 18)  # only these three duplication levels

# Function to perform one simulation
run_simulation <- function(sample_id, row_dup, shape_est, scale_est, n_control, n_treatment) {
  set.seed(sample_id) # Unique seed for each simulation
  
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
  level_idx <- ceiling(i / simulate)         
  row_dup   <- row_duplication_levels[level_idx]  
  sample_id <- (level_idx - 1) * simulate + (i %% simulate) + 1
  
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
  t <- lnCVR / se_LnCVR
  p_value <- 2*pt(abs(t),df= n_treatment+n_control -2,lower.tail = F)
  var_T = sT^2 / (nT * xT^2)
  var_C = sC^2 / (nC * xC^2)
  df_welch = (var_T + var_C)^2 / (var_T^2/(nT - 1) +  var_C^2/(nC - 1))
  p_value_welch = 2*pt(abs(t),df= df_welch,lower.tail = F)
  data.frame(
    xT = xT, sT = sT, nT = nT,
    xC = xC, sC = sC, nC = nC,
    lnCVR_naive = lnCVR_naive, lnCVR = lnCVR,
    var_lnCVR = var_lnCVR, se_LnCVR = se_LnCVR,
    t = t, p_value = p_value,p_value_welch=p_value_welch,
    row_duplication = unique(sample_data$row_duplication),
    sample = unique(sample_data$sample)
  )
}

# Calculate Master_summary using lapply
Master_summary <- do.call(rbind, lapply(split(Master, Master$sample), calculate_summary))

# Ensure row_duplication is a factor with consistent levels
Master_summary <- Master_summary %>%
  mutate(row_duplication = factor(row_duplication, levels = row_duplication_levels))

# Replace with desired path
current_wd <- getwd()
write.table(as.data.frame(Master_summary %>%
                            filter(sample != 0) %>%
                            group_by(row_duplication) %>%
                            summarise(
                              sensitivity = 100*(sum(p_value_welch <= 0.05 & (lnCVR < -0.5 | lnCVR > 0.5), na.rm = TRUE)/ n()),
                              .groups = "drop"
                            ) %>% arrange(row_duplication) %>% 
                            tidyr::pivot_wider(names_from = row_duplication,values_from=sensitivity)),row.names=F,sep="\t",quote = F,
            file=paste0(current_wd,"/duplicated_rows_small_sensitivity.txt"))
