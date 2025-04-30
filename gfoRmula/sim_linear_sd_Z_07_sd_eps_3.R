# data simulation
library(gfoRmula)
library(data.table)
library(ggplot2)
library(tidyr)
library(survey)


convert_to_long <- function(dat) {
  # Create an id column
  dat$id <- 1:nrow(dat)
  
  # Reshape the data to long format
  dat_long <- data.frame(
    id = rep(dat$id, times = 3),
    t0 = rep(0:2, each = nrow(dat)),
    treatment = unlist(dat[, c("treatment_1", "treatment_2", "treatment_3")],
                       use.names = FALSE)
  )
  
  # Assign covariates dynamically according to the correct time point
  dat_long$systolic_bp <- c(dat$systolic_bp_1, dat$systolic_bp_2, dat$systolic_bp_3)
  dat_long$diastolic_bp <- c(dat$diastolic_bp_1, dat$diastolic_bp_2, dat$diastolic_bp_3)
  dat_long$electrolyte_balance <- c(dat$electrolyte_balance_1, dat$electrolyte_balance_2,
                                    dat$electrolyte_balance_3)
  dat_long$weight <- c(dat$weight_1, dat$weight_2, dat$weight_3)
  
  # Assign Y_glomerular_filtration only at t0 = 2, otherwise NA
  dat_long$Y_glomerular_filtration <- ifelse(dat_long$t0 == 2,
                                             dat$Y_glomerular_filtration, NA)
  
  dat_long <- dat_long[order(dat_long$id), ]
  
  # Convert to data.table format 
  setDT(dat_long)
  
  return(dat_long)
}

# Function to estimate NICE ATE for given interventions
NICE_estimate_ATE <- function(data, intervention1, intervention2, nsimul = 100,
                              seed = 1234) {
  set.seed(seed)  # Ensure reproducibility
  
  data <- convert_to_long(data)
  
  # Define model parameters
  id <- 'id'
  time_name <- 't0'
  covnames <- c('systolic_bp', 'diastolic_bp', 'electrolyte_balance', 'weight',
                'treatment')
  outcome_name <- 'Y_glomerular_filtration'
  outcome_type <- 'continuous_eof'
  covtypes <- c('normal', 'normal', 'normal', 'normal', 'binary')
  histories <- c(lagged)
  histvars <- list(c('systolic_bp', 'diastolic_bp', 'electrolyte_balance', 'weight',
                     'treatment'))
  
  # Define models
  covparams <- list(
    covmodels = list(
      systolic_bp ~ lag1_systolic_bp + lag1_diastolic_bp + lag1_electrolyte_balance +
        lag1_weight + lag1_treatment + t0,
      diastolic_bp ~ lag1_systolic_bp + lag1_diastolic_bp + lag1_electrolyte_balance +
        lag1_weight + lag1_treatment + t0,
      electrolyte_balance ~ lag1_systolic_bp + lag1_diastolic_bp + lag1_electrolyte_balance +
        lag1_weight + lag1_treatment + t0,
      weight ~ lag1_systolic_bp + lag1_diastolic_bp + lag1_electrolyte_balance +
        lag1_weight + lag1_treatment + t0,
      treatment ~ lag1_treatment + systolic_bp + diastolic_bp + electrolyte_balance + 
        weight + lag1_systolic_bp + lag1_diastolic_bp + lag1_electrolyte_balance +
        lag1_weight + t0
    )
  )
  
  ymodel <- Y_glomerular_filtration ~ treatment + systolic_bp + diastolic_bp +
    electrolyte_balance + weight + lag1_treatment + lag1_systolic_bp +
    lag1_diastolic_bp + lag1_electrolyte_balance + lag1_weight 
  
  int_descript <- c('Intervention 1', 'Intervention 2')
  
  # Run gformula for the given interventions
  gform_cont_eof <- gformula(
    obs_data = data,
    id = id, time_name = time_name,
    covnames = covnames, outcome_name = outcome_name,
    outcome_type = outcome_type, covtypes = covtypes,
    covparams = covparams, ymodel = ymodel,
    intervention1.treatment = intervention1,
    intervention2.treatment = intervention2,
    int_descript = int_descript,
    histories = histories, histvars = histvars,
    nsimul = nsimul, seed = seed
  )
  
  # Extract mean outcomes for the two interventions
  mean_outcome_1 <- gform_cont_eof$result$`g-form mean`[gform_cont_eof$result$Interv. == 1]
  mean_outcome_2 <- gform_cont_eof$result$`g-form mean`[gform_cont_eof$result$Interv. == 2]
  
  # Compute ATE
  ATE <- mean_outcome_2 - mean_outcome_1
  
  # Return results
  return(list(
    ATE = ATE,
    mean_outcome_1 = mean_outcome_1,
    mean_outcome_2 = mean_outcome_2
  ))
}


# Initialize vectors to store results
result_point_estimate_111_000_vec <- numeric(1000)
result_point_estimate_011_000_vec <- numeric(1000)

# Loop over seed numbers
for (i in 1:1000) {
  seed_num <- 20250000 + i
  file_name <- paste0("simulated_data_2_seed_", seed_num, "_sdZ_07_sdEps_3.csv")
  file_path <- file.path("/Users/ruyi/Desktop/452 Final Project/simulation_data/parallel_linear_sd_Z_07_sd_eps_3", file_name)
  
  # Read data
  dat <- read.csv(file_path, header = TRUE)
  
  # Define interventions
  intervention1 <- list(static, rep(0, 3))  # (0,0,0)
  intervention2 <- list(static, rep(1, 3))  # (1,1,1)
  
  # Estimate ATE for (1,1,1) vs (0,0,0)
  result <- NICE_estimate_ATE(dat, intervention1, intervention2, nsimul = 5000, seed = 1234)
  result_point_estimate_111_000_vec[i] <- result$ATE
  
  # Estimate ATE for (0,1,1) vs (0,0,0)
  intervention2 <- list(static, c(0, 1, 1))
  result <- NICE_estimate_ATE(dat, intervention1, intervention2, nsimul = 5000, seed = 1234)
  result_point_estimate_011_000_vec[i] <- result$ATE
  
  print(seed_num)
}


# True values
tau_111_000 <- 5.041
tau_011_000 <- 3.021

H <- 1000

# --- (1,1,1) vs (0,0,0) ---
estimates_111_000 <- result_point_estimate_111_000_vec

# Relative Bias
relative_bias_111_000 <- (mean(estimates_111_000) - tau_111_000) / tau_111_000

# MCSD
mcsd_111_000 <- sqrt(sum((estimates_111_000 - mean(estimates_111_000))^2) / (H - 1))

# RMSE
rmse_111_000 <- sqrt(mean((estimates_111_000 - tau_111_000)^2))


# --- (0,1,1) vs (0,0,0) ---
estimates_011_000 <- result_point_estimate_011_000_vec

# Relative Bias
relative_bias_011_000 <- (mean(estimates_011_000) - tau_011_000) / tau_011_000

# MCSD
mcsd_011_000 <- sqrt(sum((estimates_011_000 - mean(estimates_011_000))^2) / (H - 1))

# RMSE
rmse_011_000 <- sqrt(mean((estimates_011_000 - tau_011_000)^2))


result_point_estimate_111_000_vec
result_point_estimate_011_000_vec


# Print results
cat("=== ATE (1,1,1) vs (0,0,0) ===\n")
cat("Relative Bias:", relative_bias_111_000, "\n")
cat("MCSD:", mcsd_111_000, "\n")
cat("RMSE:", rmse_111_000, "\n\n")

cat("=== ATE (0,1,1) vs (0,0,0) ===\n")
cat("Relative Bias:", relative_bias_011_000, "\n")
cat("MCSD:", mcsd_011_000, "\n")
cat("RMSE:", rmse_011_000, "\n")


