# =============================================================================
#
# Global Sensitivity Analysis using LHS-PRCC
#
# Description: This script performs a global sensitivity analysis on a 16-state
#              epidemiological model. It uses Latin Hypercube Sampling (LHS)
#              from the 'lhs' package to generate parameter sets, runs the model
#              in parallel, and then calculates Partial Rank Correlation
#              Coefficients (PRCC) to quantify parameter impact on key outcomes.
#
# =============================================================================


# --- 1. Load Necessary Libraries ---
# -----------------------------------------------------------------------------
library(deSolve)      # For solving Ordinary Differential Equations (ODEs)
library(lhs)          # For Latin Hypercube Sampling (LHS)
library(sensitivity)  # For PRCC (pcc) function
library(dplyr)        # For data manipulation and pipelines
library(ggplot2)      # For creating high-quality plots
library(doParallel)   # For enabling parallel processing
library(purrr)        # For functional programming, useful with pmap

# --- Helper function to fit Beta distribution from quantiles ---
# This function is crucial for sampling from Beta distributions.
# It finds the shape parameters (alpha and beta) that correspond to
# given quantiles of a Beta distribution.
get.beta.par <- function(p = c(0.025, 0.975), q = c(0.1, 0.9), 
                         show.output = FALSE, plot = FALSE) {
  
  # Ensure quantiles are within (0, 1)
  q <- pmax(pmin(q, 0.9999), 0.0001)
  
  # Objective function to minimize: sum of squared differences between
  # target quantiles and quantiles from a Beta distribution with given shape parameters.
  objective_function <- function(parms, target_q, target_p) {
    alpha <- parms[1]
    beta <- parms[2]
    
    # Ensure shape parameters are positive
    if (alpha <= 0 || beta <= 0) {
      return(Inf)
    }
    
    # Calculate quantiles for the current Beta distribution
    calculated_q <- qbeta(target_p, shape1 = alpha, shape2 = beta)
    
    # Sum of squared errors
    sum_sq_err <- sum((calculated_q - target_q)^2)
    return(sum_sq_err)
  }
  
  # Initial guess for shape parameters
  initial_guess <- c(alpha = 2, beta = 2)
  
  # Optimize to find the best shape parameters
  opt_result <- optim(
    par = initial_guess,
    fn = objective_function,
    target_q = q,
    target_p = p,
    method = "L-BFGS-B",
    lower = c(0.001, 0.001) # Shape parameters must be positive
  )
  
  fitted_parms <- opt_result$par
  
  if (show.output) {
    print(paste("Fitted shape1 (alpha):", fitted_parms[1]))
    print(paste("Fitted shape2 (beta):", fitted_parms[2]))
  }
  
  if (plot) {
    curve(dbeta(x, shape1 = fitted_parms[1], shape2 = fitted_parms[2]), 
          from = 0, to = 1, 
          ylab = "Density", xlab = "Parameter Value", main = "Fitted Beta Distribution")
    points(q, dbeta(q, shape1 = fitted_parms[1], shape2 = fitted_parms[2]), 
           col = "red", pch = 19)
    abline(v = q, lty = 2, col = "blue")
  }
  
  return(fitted_parms)
}


# --- 2. Setup Parallel Processing ---
# -----------------------------------------------------------------------------
# Use a reasonable number of cores, leaving some for system stability
num_cores <- max(1, parallel::detectCores() - 2)
cl <- makeCluster(num_cores)
registerDoParallel(cl)
cat("Parallel processing enabled on", num_cores, "cores.\n")


# --- 3. Define Model, Initial Conditions, and Time ---
# -----------------------------------------------------------------------------
# Set a seed for reproducibility of the random sampling process
set.seed(123)

# --- Initial Conditions (Latest Version) ---
N <- 5.18e6
Nh <- 37588 * 0.85
Nc <- N - Nh
RCh <- round(0.124 * Nh); RIh <- round(0.119 * 0.314 * 0.377 * Nh); SIh <- round(0.119 * 0.314 * 0.56 * Nh);
Sh <- Nh - RCh - RIh - SIh;
RCc <- round(0.257 * Nc); RIc <- round(0.001 * Nc); SIc <- round(0.01 * Nc);
Sc <- Nc - RCc - RIc - SIc;

# State vector for the ODE solver
y_init <- c(Sh = Sh, RCh = RCh, RIh = RIh, SIh = SIh, Sc = Sc, RCc = RCc, RIc = RIc, SIc = SIc,
            Drh = 0, Dsh = 0, Drc = 0, Dsc = 0, 
            C_RIh = 0, C_SIh = 0, C_RIc = 0, C_SIc = 0)

# --- Simulation Time ---
# We will analyze sensitivity at the end of a 10-year period
t_start <- 0
t_end <- 365 * 10
times <- seq(t_start, t_end, by = 365) # Annual output is sufficient
final_time <- max(times)

# --- ODE Model (Latest Version with Constrained Flows) ---
constrain_flow <- function(flow_rate, population, dt = 1) {
  max_flow <- population / dt
  return(min(flow_rate, max_flow))
}

infection_model_ode <- function(t, y, params) {
  with(as.list(c(y, params)), {
    # Ensure state variables are non-negative
    Sh <- max(0, y["Sh"]); RCh <- max(0, y["RCh"]); RIh <- max(0, y["RIh"]); SIh <- max(0, y["SIh"]);
    Sc <- max(0, y["Sc"]); RCc <- max(0, y["RCc"]); RIc <- max(0, y["RIc"]); SIc <- max(0, y["SIc"]);
    
    # Calculate current population sizes for density-dependent terms
    Nh_current <- max(1, Sh + RCh + RIh + SIh)
    Nc_current <- max(1, Sc + RCc + RIc + SIc)
    N_total_living <- Nh_current + Nc_current
    
    # Calculate all flows using the constrain_flow helper function for robustness
    flow_beta_eh_Sh <- constrain_flow(beta_eh * Sh, Sh); flow_beta_hh_Sh <- constrain_flow(beta_hh * Sh * RCh / Nh_current, Sh);
    flow_beta_sh_sih_Sh <- constrain_flow(beta_sh_sih * Sh * SIh / Nh_current, Sh); flow_alpha_dis_Sh <- constrain_flow(alpha_dis * Sh, Sh); flow_mu_d_Sh <- constrain_flow(mu_d * Sh, Sh);
    flow_delta_rch_rih_RCh <- constrain_flow(delta_rch_rih * RCh, RCh); flow_gamma_rch_sh_RCh <- constrain_flow(gamma_rch_sh * RCh, RCh);
    flow_alpha_dis_RCh <- constrain_flow(alpha_dis * RCh, RCh); flow_delta_rch_sih_RCh <- constrain_flow(delta_rch_sih * RCh, RCh); flow_mu_d_RCh <- constrain_flow(mu_d * RCh, RCh);
    flow_gamma_rih_rch_RIh <- constrain_flow(gamma_rih_rch * RIh, RIh); flow_mu_rih_RIh <- constrain_flow(mu_rih * RIh, RIh);
    flow_gamma_sih_sh_SIh <- constrain_flow(gamma_sih_sh * SIh, SIh); flow_mu_sih_SIh <- constrain_flow(mu_sih * SIh, SIh);
    flow_beta_sc_ec_Sc <- constrain_flow(beta_sc_ec * Sc, Sc); flow_beta_sc_ac_Sc <- constrain_flow(beta_sc_ac * Sc * 0.33, Sc);
    flow_beta_hc_Sc <- constrain_flow(beta_hc * Sc * RCc / Nc_current, Sc); flow_beta_sc_sic_Sc <- constrain_flow(beta_sc_sic * Sc * SIc / Nc_current, Sc);
    flow_alpha_adm_Sc <- constrain_flow(alpha_adm * Sc, Sc); flow_mu_d_Sc <- constrain_flow(mu_d * Sc, Sc);
    flow_delta_rcc_ric_RCc <- constrain_flow(delta_rcc_ric * RCc, RCc); flow_delta_rcc_sic_RCc <- constrain_flow(delta_rcc_sic * RCc, RCc);
    flow_gamma_rcc_sc_RCc <- constrain_flow(gamma_rcc_sc * RCc, RCc); flow_alpha_adm_RCc <- constrain_flow(alpha_adm * RCc, RCc); flow_mu_d_RCc <- constrain_flow(mu_d * RCc, RCc);
    flow_gamma_ric_rcc_RIc <- constrain_flow(gamma_ric_rcc * RIc, RIc); flow_mu_ric_RIc <- constrain_flow(mu_ric * RIc, RIc); flow_alpha_adm_RIc <- constrain_flow(alpha_adm * RIc, RIc);
    flow_gamma_sic_sc_SIc <- constrain_flow(gamma_sic_sc * SIc, SIc); flow_mu_sic_SIc <- constrain_flow(mu_sic * SIc, SIc); flow_alpha_adm_SIc <- constrain_flow(alpha_adm * SIc, SIc);
    
    # Define the differential equations
    dSh_dt <- -flow_beta_eh_Sh - flow_beta_hh_Sh - flow_beta_sh_sih_Sh - flow_alpha_dis_Sh - flow_mu_d_Sh + flow_gamma_sih_sh_SIh + flow_gamma_rch_sh_RCh + flow_alpha_adm_Sc;
    dRCh_dt <- -flow_delta_rch_rih_RCh - flow_gamma_rch_sh_RCh - flow_alpha_dis_RCh - flow_delta_rch_sih_RCh - flow_mu_d_RCh + flow_beta_eh_Sh + flow_beta_hh_Sh + flow_gamma_rih_rch_RIh + flow_alpha_adm_RCc;
    dRIh_dt <- -flow_gamma_rih_rch_RIh - flow_mu_rih_RIh + flow_delta_rch_rih_RCh + flow_alpha_adm_RIc;
    dSIh_dt <- -flow_gamma_sih_sh_SIh - flow_mu_sih_SIh + flow_beta_sh_sih_Sh + flow_delta_rch_sih_RCh + flow_alpha_adm_SIc;
    dSc_dt <- -flow_beta_sc_ec_Sc - flow_beta_sc_ac_Sc - flow_beta_hc_Sc - flow_beta_sc_sic_Sc - flow_alpha_adm_Sc - flow_mu_d_Sc + flow_gamma_sic_sc_SIc + flow_gamma_rcc_sc_RCc + flow_alpha_dis_Sh + mu_b * N_total_living + mu_m * N_total_living;
    dRCc_dt <- -flow_delta_rcc_ric_RCc - flow_delta_rcc_sic_RCc - flow_gamma_rcc_sc_RCc - flow_alpha_adm_RCc - flow_mu_d_RCc + flow_beta_sc_ec_Sc + flow_beta_sc_ac_Sc + flow_beta_hc_Sc + flow_gamma_ric_rcc_RIc + flow_alpha_dis_RCh;
    dRIc_dt <- -flow_gamma_ric_rcc_RIc - flow_mu_ric_RIc - flow_alpha_adm_RIc + flow_delta_rcc_ric_RCc;
    dSIc_dt <- -flow_gamma_sic_sc_SIc - flow_mu_sic_SIc - flow_alpha_adm_SIc + flow_beta_sc_sic_Sc + flow_delta_rcc_sic_RCc;
    dDrh_dt <- flow_mu_rih_RIh; dDsh_dt <- flow_mu_sih_SIh; dDrc_dt <- flow_mu_ric_RIc; dDsc_dt <- flow_mu_sic_SIc;
    dC_RIh_dt <- flow_delta_rch_rih_RCh + flow_alpha_adm_RIc; dC_SIh_dt <- flow_beta_sh_sih_Sh + flow_delta_rch_sih_RCh + flow_alpha_adm_SIc;
    dC_RIc_dt <- delta_rcc_ric * RCc; dC_SIc_dt <- flow_beta_sc_sic_Sc + flow_delta_rcc_sic_RCc;
    
    # Return the derivatives as a list
    return(list(c(dSh_dt, dRCh_dt, dRIh_dt, dSIh_dt, dSc_dt, dRCc_dt, dRIc_dt, dSIc_dt, dDrh_dt, dDsh_dt, dDrc_dt, dDsc_dt, dC_RIh_dt, dC_SIh_dt, dC_RIc_dt, dC_SIc_dt)))
  })
}


# --- 4. Define Parameter Ranges and Generate Samples with LHS ---
# -----------------------------------------------------------------------------
cat("\nGenerating parameter samples using Latin Hypercube Sampling...\n")

# --- NEW: Define Parameter Distributions based on provided data ---
# This uses the provided quantiles to fit either Beta or Uniform distributions.
# It also calculates the mean for each parameter based on the fitted distribution.

param_data <- data.frame(
  param = c("beta_hh", "beta_eh", "beta_hc", "beta_sc_ac", "beta_sc_ec", "beta_sh_sih", "delta_rcc_ric",
            "delta_rch_rih", "mu_sih", "mu_rih", "gamma_sih_sh", "gamma_rch_sh", "gamma_rih_rch",
            "delta_rch_sih", "beta_sc_sic", "gamma_sic_sc", "gamma_rcc_sc", "gamma_ric_rcc",
            "mu_sic", "mu_ric", "delta_rcc_sic", "alpha_adm", "alpha_dis", "mu_b", "mu_m", "mu_d"),
  dist = c(rep("beta", 6), rep("unif", 20)), # First 6 are Beta, rest are Uniform
  q_025 = c(0.0024, 1e-6, 0.00174, 0.000024, 1.89e-07, 0.14, 0.00001,
            0.00010513, 0.00233333, 0.00254333, 0.08, 0.00156438, 0.014, 0.000592603,
            0.00000707, 0.1, 0.00179726, 0.014, 0.00182000,   0.00050, 0.000592603,
            0.000316 * 0.7, 0.235 * 0.7, 0.00002518 * 0.7, 0.000002 * 0.7, 0.00001353 * 0.7),
  q_975 = c(0.013597, 0.017507, 0.010412, 0.00355, 0.00105, 0.8917, 0.00016,
            0.06779661, 0.00390000, 0.00554667, 0.148571429, 0.01023333, 0.03466667,
            0.001100548, 0.00000877, 0.185714286, 0.02350000, 0.02600000, 0.00338000,
            0.00124, 0.001100548, 0.000316 * 1.3, 0.235 * 1.3, 0.00002518 * 1.3,
            0.000002 * 1.3, 0.00001353 * 1.3),
  provided_mean = c(0.007909485, 0.002670832, 0.005249949, 0.001469924, 0.00029738,
                    0.227154355, rep(NA, 20)), # Provided means for Beta distributions
  stringsAsFactors = FALSE
)

# Process Beta distributions
beta_params <- param_data %>% filter(dist == "beta")
beta_fits <- beta_params %>%
  rowwise() %>%
  mutate(
    fit = list(get.beta.par(p = c(0.025, 0.975), q = c(q_025, q_975), show.output = FALSE, plot = FALSE)),
    shape1 = fit[1],
    shape2 = fit[2],
    # Calculate mean from fitted parameters if provided_mean is NA
    calc_mean = shape1 / (shape1 + shape2), 
    mean = if_else(is.na(provided_mean), calc_mean, provided_mean) # Use provided_mean if available
  ) %>%
  ungroup() %>%
  select(param, dist, q_025, q_975, provided_mean, shape1, shape2, calc_mean, mean)

# Process Uniform distributions
non_beta_params <- param_data %>% 
  filter(dist == "unif") %>%
  mutate(
    # Calculate mean for uniform as (min+max)/2
    calc_mean = (q_025 + q_975) / 2, 
    mean = if_else(is.na(provided_mean), calc_mean, provided_mean) # Use provided_mean if available
  ) %>%
  select(param, dist, q_025, q_975, provided_mean, calc_mean, mean)

# Combine and prepare the final parameter data structure for sampling
final_param_data <- bind_rows(beta_fits, non_beta_params) %>%
  mutate(
    # Ensure shape1/shape2 are NA for uniform distributions
    shape1 = if_else(dist == "beta", shape1, NA_real_),
    shape2 = if_else(dist == "beta", shape2, NA_real_)
  ) %>%
  select(param, dist, q_025, q_975, provided_mean, shape1, shape2, calc_mean, mean)

# Create the epi_param_distributions list for sampling
epi_param_distributions <- setNames(
  purrr::pmap(list(final_param_data$param, final_param_data$dist, final_param_data$q_025, final_param_data$q_975, 
                   final_param_data$shape1, final_param_data$shape2), 
              function(p, d, min_val, max_val, s1, s2) {
                if (d == "beta") {
                  # For Beta distribution, we store shape1 and shape2
                  list(dist = "beta", shape1 = s1, shape2 = s2, min = min_val, max = max_val) 
                } else {
                  # For Uniform, we store min and max
                  list(dist = "unif", min = min_val, max = max_val)
                }
              }),
  final_param_data$param
)

# Verify the structure of the distributions
print(epi_param_distributions)

# Extract parameter names and set number of samples
param_names_to_vary <- names(epi_param_distributions)
n_params <- length(param_names_to_vary)
n_samples <- 1000 # Number of simulations (increase for higher precision, e.g., 1000)

# --- Use randomLHS from the 'lhs' package ---
# This generates a matrix of probabilities (0-1) for each parameter
lhs_matrix <- randomLHS(n = n_samples, k = n_params)

# Create an empty data frame for parameter values
param_samples <- data.frame(matrix(NA, nrow = n_samples, ncol = n_params))
colnames(param_samples) <- param_names_to_vary

# Convert LHS probabilities to parameter values using their respective distributions
for (i in 1:n_params) {
  param_name <- param_names_to_vary[i]
  dist_info <- epi_param_distributions[[param_name]]
  
  if (dist_info$dist == "unif") {
    # Use the quantile function for the uniform distribution (qunif)
    param_samples[, i] <- qunif(lhs_matrix[, i], min = dist_info$min, max = dist_info$max)
  } else if (dist_info$dist == "beta") {
    # Use the quantile function for the Beta distribution (qbeta)
    param_samples[, i] <- qbeta(lhs_matrix[, i], shape1 = dist_info$shape1, shape2 = dist_info$shape2)
  }
  # Note: Add 'else if' blocks here for other distributions if needed.
}
cat("Parameter sample generation complete.\n")
print(param_samples) 


# --- 5. Run Simulations in Parallel ---
# -----------------------------------------------------------------------------
cat("\nRunning", n_samples, "simulations in parallel...\n")

results_list <- foreach(i = 1:n_samples, .packages = c("deSolve", "dplyr")) %dopar% {
  current_params <- as.list(param_samples[i, ])
  
  # Run the ODE solver
  out <- try(ode(y = y_init, times = times, func = infection_model_ode, parms = current_params, method = "lsoda"), silent = TRUE)
  
  # Return the final row of results if successful
  if (!inherits(out, "try-error")) {
    final_row <- as.data.frame(out) %>% tail(1)
    final_row$sample_id <- i # Add sample ID for tracking
    return(final_row)
  } else {
    # print(paste("Simulation failed for sample:", i)) # Optional: uncomment for debugging
    return(NULL) # Mark failed simulations
  }
}

# Stop the parallel cluster
stopCluster(cl)
cat("Simulations complete.\n")

# Combine results and remove failed runs
all_results <- bind_rows(results_list)
cat(nrow(all_results), "out of", n_samples, "simulations completed successfully.\n")

# If no simulations succeeded, stop here
if (nrow(all_results) == 0) {
  stop("All simulations failed. Please check model and parameter ranges.")
}

# Filter the parameter samples (X) to include only those from successful runs
X <- param_samples[all_results$sample_id, ]


# --- 6. Calculate PRCC for Key Outcomes ---
# -----------------------------------------------------------------------------
cat("\nCalculating PRCC for key epidemiological outcomes...\n")

# Define the outcomes of interest (focused on resistant infections)
outcome_vars <- c("C_RIh", "C_RIc", "Drh", "Drc")
outcome_titles <- c("Cumulative Hospital RI", 
                    "Cumulative Community RI",
                    "Deaths from Hospital RI", 
                    "Deaths from Community RI")

# List to store PRCC results
prcc_results_list <- list()

for (i in seq_along(outcome_vars)) {
  outcome_var <- outcome_vars[i]
  outcome_title <- outcome_titles[i]
  
  # Extract the outcome vector (y) for the current outcome
  y <- all_results[[outcome_var]]
  
  # Compute PRCC with bootstrapping for confidence intervals
  # nboot=100 is fast for testing; use nboot=1000 for final analysis
  # Rank=TRUE is important for PRCC calculation
  prcc_result <- pcc(X = X, y = y, rank = TRUE, nboot = 100) 
  prcc_results_list[[outcome_title]] <- prcc_result
}
cat("PRCC calculations finished.\n")


# --- 7. Visualize PRCC Results as a Heatmap ---
# -----------------------------------------------------------------------------
cat("Generating PRCC heatmap...\n")

# Combine all PRCC results into a single data frame for plotting
all_prcc_df <- lapply(names(prcc_results_list), function(outcome_name) {
  prcc_df <- prcc_results_list[[outcome_name]]$PRCC
  prcc_df$Parameter <- rownames(prcc_df)
  prcc_df$Outcome <- outcome_name
  return(prcc_df)
}) %>% bind_rows()

# Ensure parameters are ordered for better visualization (e.g., alphabetical or by significance)
# You might want to reorder parameters based on their average absolute PRCC for better clarity.
# For now, we'll use the order as they appear.

# Create the heatmap
prcc_heatmap <- ggplot(all_prcc_df, aes(x = Parameter, y = Outcome, fill = original)) +
  geom_tile(color = "white", linewidth = 0.5) +
  # Use scale_fill_gradient2 for a divergent color scale
  scale_fill_gradient2(
    name = "PRCC",
    low = "steelblue", # Color for negative PRCC
    mid = "white",     # Color for PRCC near zero
    high = "firebrick",# Color for positive PRCC
    midpoint = 0,
    limits = c(-1, 1)  # Ensure the scale is from -1 to 1
  ) +
  # Add the PRCC values as text inside each tile
  geom_text(aes(label = sprintf("%.2f", original)), color = "white", size = 3, fontface="bold") + # Changed text color to black for better visibility on all fills
  labs(
    title = "PRCC Global Sensitivity Analysis",
    subtitle = paste("Based on", nrow(all_results), "successful model runs from", n_samples, "LHS samples"),
    x = "Input Parameter",
    y = "Model Outcome"
  ) +
  theme_minimal(base_size = 10) + # Reduced base size slightly for more parameters
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8), # Reduced text size for x-axis
    axis.text.y = element_text(size = 8), # Reduced text size for y-axis
    legend.position = "right",
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  ) +
  # Ensure parameter labels fit if there are many
  guides(fill = guide_colorbar(barwidth = 1.5, barheight = 8))

# Print the heatmap
print(prcc_heatmap)

# --- Optional: Save the heatmap ---
# ggsave("prcc_heatmap.png", plot = prcc_heatmap, width = 12, height = 8, dpi = 300)

cat("\nAnalysis and visualization complete.\n")