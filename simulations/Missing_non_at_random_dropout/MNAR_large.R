# Load required packages
library(haven)
library(fitdistrplus)
library(dplyr)
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
dropout_levels <- seq(0, 0.5, by = 0.1)
n_control <- 250
n_treatment <- 250

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
    dropout = unique(sample_data$dropout),
    sample = unique(sample_data$sample)
  )
}

# Calculate Master_summary using lapply
Master_summary <- do.call(rbind, lapply(split(Master, Master$sample), calculate_summary))

# Set PDF path
# Replace with desired path
current_wd <- getwd()
write.table(as.data.frame(Master_summary %>%
                            filter(sample != 0) %>%
                            group_by(dropout) %>%
                            summarise(
                              specificity = 100*(1-sum(p_value_welch <= 0.05 & (lnCVR < -0.5 | lnCVR > 0.5), na.rm = TRUE)/ n()),
                              .groups = "drop"
                            ) %>%
                            mutate(specificity = round(specificity, 3))),row.names=F,sep="\t",quote = F,
            file=paste0(current_wd,"/specificity_MNAR_large.txt"))

