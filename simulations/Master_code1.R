# Load packages
library(dplyr)
library(tidyr)
library(doParallel)
library(foreach)
# -----------------------------
# Fixed Parameters
# -----------------------------
shape_est <- 50.16049
rate_est  <- 6.29346
scale_est <- 1 / rate_est

# Meta-analysis outputs for PI

# RVE: Correlated Effects Model with Small-Sample Corrections 
# 
# Model: lnCVR ~ 1 
# 
# Number of studies = 175 
# Number of outcomes = 229 (min = 1 , mean = 1.31 , median = 1 , max = 6 )
# Rho = 0.8 
# I.sq = 67.41011 
# Tau.sq = 0.02674155 
# 
# Estimate StdErr t-value dfs P(|t|>) 95% CI.L 95% CI.U Sig
# 1 X.Intercept.  -0.0358 0.0172   -2.08 143  0.0392  -0.0698  -0.0018  **
#   ---
#   Signif. codes: < .01 *** < .05 ** < .10 *
#   ---
#   Note: If df < 4, do not trust the results

M_star    <- -0.0358
Tau_sq    <- 0.02674155  # T²
StdErr    <- 0.0172
VM_star   <- StdErr^2
sd_pi     <- sqrt(Tau_sq + VM_star)
for(multiplier in c(3,4)){
  # multiple=3: 3-sigma (normal). For same coverage under t: qt(1 - (1 - (pnorm(3)-pnorm(-3)))/2, df=143) ≈ 2.98

lower_pi  <- M_star - multiplier * sd_pi
upper_pi  <- M_star + multiplier * sd_pi

# Reference CV under homogeneity (your original)
CV_ref <- 1 / sqrt(shape_est)
mu_ref <- shape_est * scale_est   # ≈ 7.97

# Simulation settings
n_sim     <- 1000  # Per parameter combination
cores     <- detectCores() - 1  # For parallel

# -----------------------------
# Modified Simulation Function
# -----------------------------
run_simulation <- function(sample_id, n_per_arm, subgroup_prop, effect_mag,
                           dropout_prop, fraud_prop,
                           shape_est, scale_est, 
                           M_star, Tau_sq, CV_ref, mu_ref) {
  
  set.seed(sample_id)

  # Always apply realistic between-trial heterogeneity in CV ratio
  true_lnCVR <- rnorm(1, mean = M_star, sd = sqrt(Tau_sq))
  
  CV_C <- CV_ref                  # control arm CV stays at reference level
  CV_T <- CV_C * exp(true_lnCVR)  # treatment arm CV adjusted according to sampled ratio
  
  shape_C <- 1 / CV_C^2
  shape_T <- 1 / CV_T^2
  
  scale_C <- mu_ref / shape_C
  scale_T <- mu_ref / shape_T
  # ── Generate base data ─────────────────────────────────────────────────────
  sim_data <- data.frame(
    original_row = 1:(2 * n_per_arm),
    group        = c(rep(0, n_per_arm), rep(1, n_per_arm)),
    sample       = sample_id,
    r_a1c        = NA_real_
  )
  
  # Generate HbA1c values using possibly different parameters per arm
  sim_data$r_a1c[sim_data$group == 0] <- rgamma(n_per_arm, shape = shape_C, scale = scale_C)
  sim_data$r_a1c[sim_data$group == 1] <- rgamma(n_per_arm, shape = shape_T, scale = scale_T)
  
  # ── Apply heterogeneous treatment effect (HTE) if present ──────────────────
  n_select <- round(n_per_arm * subgroup_prop)
  if (n_select > 0 && effect_mag != 0) {
    treat_indices <- which(sim_data$group == 1)
    selected_indices <- sample(treat_indices, n_select)
    sim_data$r_a1c[selected_indices] <- sim_data$r_a1c[selected_indices] + effect_mag
    sim_data$r_a1c <- pmax(sim_data$r_a1c, 0.01)  # prevent negative/zero values
  }
  
  # ── Apply MNAR dropout (control arm, highest values dropped) ───────────────
  if (dropout_prop > 0) {
    n_dropout <- round(n_per_arm * dropout_prop)
    n_dropout <- min(n_dropout, n_per_arm - 2)  # keep at least 2 observations
    
    if (n_dropout > 0) {
      control_rows <- sim_data[sim_data$group == 0, ]
      sorted_idx <- order(-control_rows$r_a1c)  # highest first
      rows_to_drop <- control_rows$original_row[sorted_idx[1:n_dropout]]
      sim_data <- sim_data[!sim_data$original_row %in% rows_to_drop, ]
    }
  }
  
  # ── Apply fraud (duplicate best/lowest value in treatment arm) ─────────────
  if (fraud_prop > 0) {
    n_fraud <- round(n_per_arm * fraud_prop)
    n_fraud <- min(n_fraud, n_per_arm - 1)
    
    if (n_fraud > 0) {
      treat_data <- sim_data[sim_data$group == 1, ]
      best_value <- min(treat_data$r_a1c)  # assuming lower HbA1c = better
      sorted_treat <- treat_data[order(-treat_data$r_a1c), ]  # highest first? Wait — adjust logic if needed
      rows_to_fraud <- sorted_treat$original_row[1:n_fraud]
      sim_data$r_a1c[sim_data$original_row %in% rows_to_fraud] <- best_value
    }
  }
  
  # ── Compute summary statistics and LnCVR ───────────────────────────────────
  control <- sim_data$r_a1c[sim_data$group == 0]
  treat   <- sim_data$r_a1c[sim_data$group == 1]
  
  nC <- length(control)
  nT <- length(treat)
  
  if (nC < 2 || nT < 2) {
    return(list(data = sim_data, LnCVR = NA_real_, is_significant = FALSE))
  }
  
  xC <- mean(control)
  xT <- mean(treat)
  sC <- sd(control)
  sT <- sd(treat)
  
  cv_control <- if (xC > 0) sC / xC else NA_real_
  cv_treat   <- if (xT > 0) sT / xT else NA_real_
  
  lnCVR_biasedcorrected <- log(cv_treat / cv_control)+(1/2)*(1/(nT-1)-1/(nC-1))+(1/2)*((sC^2 /(nC*xC^2))-(sT^2 /(nT*xT^2)))
  
  LnCVR <- if (!is.na(cv_control) && cv_control > 0) lnCVR_biasedcorrected else NA_real_
  
  if (is.na(LnCVR)) {
    return(list(data = sim_data, LnCVR = NA_real_, is_significant = FALSE))
  }
  
  var_LnCVR <- sC^2/(nC*xC^2) + sC^4/(2*nC^2*xC^4)+0.5*nC/((nC-1)^2)+
    sT^2/(nT*xT^2) + sT^4/(2*nT^2*xT^4) + 0.5*nT/((nT-1)^2)
  se_LnCVR  <- sqrt(var_LnCVR)
  
  if (se_LnCVR == 0 || is.na(se_LnCVR)) {
    return(list(data = sim_data, LnCVR = LnCVR, is_significant = FALSE))
  }
  
  t_stat    <- LnCVR / se_LnCVR
  
  V_T <- sT^2 / (nT * xT^2) + sT^4 / (2 * nT^2 * xT^4) + 0.5*nT/((nT-1)^2)
  V_C <- sC^2 / (nC * xC^2) + sC^4 / (2 * nC^2 * xC^4) + 0.5*nC/((nC-1)^2)
   
  
  # Welch-Satterthwaite df for parallel-arm lnCVR
  df_welch <- (V_T + V_C)^2 / ( (V_T^2 / (nT - 1)) + (V_C^2 / (nC - 1)) )
  
  if (is.na(df_welch) || df_welch <= 0) {
    return(list(data = sim_data, LnCVR = LnCVR, is_significant = FALSE))
  }
  
  p_value   <- 2 * pt(-abs(t_stat), df = df_welch)
  is_significant <- p_value < 0.05
  
  list(data = sim_data, LnCVR = LnCVR, is_significant = is_significant)
}

# -----------------------------
# Parameter Grids
# -----------------------------
# Non-fraud grid for specificity (fraud_prop=0, vary HTE and dropout)
non_fraud_grid <- expand_grid(
  n_per_arm = c(20, 250),
  subgroup_prop = seq(0, 0.5, by = 0.1),
  effect_mag = seq(0, -2.0, by = -0.2),
  dropout_prop = seq(0, 0.5, by = 0.1),
  fraud_prop = 0
)

# Fraud grid for sensitivity (set HTE/dropout to 0, vary fraud)
fraud_grid <- expand_grid(
  n_per_arm = c(20, 250),
  subgroup_prop = 0,
  effect_mag = 0,
  dropout_prop = 0,
  fraud_prop = seq(0.2, 0.9, by = 0.1)
)

# Combine grids
full_grid <- bind_rows(non_fraud_grid, fraud_grid) %>% mutate(scenario = ifelse(fraud_prop == 0, "non_fraud", "fraud"))

# -----------------------------
# Run Simulations in Parallel
# -----------------------------
cl <- makeCluster(cores)
registerDoParallel(cl)

sim_results <- foreach(i = 1:nrow(full_grid), .combine = 'rbind', .packages = c('dplyr')) %dopar% {
  params <- full_grid[i, ]
  
  outside_count <- 0
  for (sim_id in 1:n_sim) {
    result <- run_simulation(sample_id = sim_id + i * n_sim,
                             n_per_arm = params$n_per_arm,
                             subgroup_prop = params$subgroup_prop,
                             effect_mag = params$effect_mag,
                             dropout_prop = params$dropout_prop,
                             fraud_prop = params$fraud_prop,
                             shape_est = shape_est,
                             scale_est = scale_est,
                             M_star    = M_star,
                             Tau_sq    = Tau_sq,
                             CV_ref    = CV_ref,
                             mu_ref    = mu_ref)
    LnCVR <- result$LnCVR
    if (!is.na(LnCVR) && (LnCVR < lower_pi || LnCVR > upper_pi) && result$is_significant) {
      outside_count <- outside_count + 1
    }
  }
  
  proportion_outside <- outside_count / n_sim
  params %>% mutate(proportion_outside = proportion_outside)
}

stopCluster(cl)
# -----------------------------
# Output to CSV (with approximate counts: FP, TP, FN, TN per combination & Sensitivity & Specificity)
# -----------------------------
results_for_csv <- sim_results %>%
  mutate(
    # Raw counts (approximate, rounded)
    FP = ifelse(scenario == "non_fraud", round(proportion_outside * n_sim), NA_real_),
    TN = ifelse(scenario == "non_fraud", round((1 - proportion_outside) * n_sim), NA_real_),
    TP = ifelse(scenario == "fraud",     round(proportion_outside * n_sim), NA_real_),
    FN = ifelse(scenario == "fraud",     round((1 - proportion_outside) * n_sim), NA_real_),
    
    # Calculated performance metrics
    Specificity = ifelse(scenario == "non_fraud", 1 - proportion_outside, NA_real_),
    Sensitivity = ifelse(scenario == "fraud",     proportion_outside,     NA_real_)
  ) %>%
  dplyr::select(
    scenario,
    n_per_arm, subgroup_prop, effect_mag, dropout_prop, fraud_prop,
    FP, TN, TP, FN,
    Specificity, Sensitivity
  )
write.csv(results_for_csv[(results_for_csv$effect_mag == 0 & results_for_csv$subgroup_prop ==0) |
                            (results_for_csv$subgroup_prop > 0 & results_for_csv$dropout_prop == 0), ],
          file = paste0("simulation_results_with_metrics_",multiplier,"sigma.csv"),row.names = FALSE)
}
