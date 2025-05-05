# -------------------------------------------------------------------------
# Causal Simulation Framework with Variable Cluster Sizes
# and Correctly Specified Outcome Model
# -------------------------------------------------------------------------

# Load required packages
library(mvtnorm)
library(SoftBart)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(gridExtra)
library(knitr)
library(parallel)
library(survival)  # for Cox models
library(coxme)     # for Cox mixed models with frailty
library(frailtypack) # alternative for frailty models

# -------------------------------------------------------------------------
# 1. Data Simulation Function with Variable Cluster Sizes
# -------------------------------------------------------------------------

simulate_causal_data_variable_clusters <- function(seed = 123, 
                                                   n = 2000, 
                                                   p = 2, 
                                                   num_groups = 5,
                                                   min_cluster_size = 50,
                                                   max_cluster_size = 500,
                                                   rho_car = 0.7, 
                                                   sigma_r2_true = 0.5,
                                                   true_sigma2 = 0.25, 
                                                   lambda_c = 0.005) {
  
  set.seed(seed)
  
  # Generate variable cluster sizes
  if (n < num_groups * min_cluster_size) {
    stop("Total sample size too small to accommodate minimum cluster size.")
  }
  
  # First assign the minimum number to each cluster
  base_sizes <- rep(min_cluster_size, num_groups)
  remaining <- n - sum(base_sizes)
  
  # Distribute remaining observations randomly but respecting max size
  while (remaining > 0) {
    # Find clusters that haven't reached max size
    available <- which(base_sizes < max_cluster_size)
    
    if (length(available) == 0) break  # All clusters at max size
    
    # Randomly select a cluster
    selected <- sample(available, 1)
    base_sizes[selected] <- base_sizes[selected] + 1
    remaining <- remaining - 1
  }
  
  # If we still have remaining observations, distribute them anyway
  if (remaining > 0) {
    warning(paste("Some clusters will exceed max_cluster_size. Remaining:", remaining))
    # Distribute remaining observations randomly
    for (i in 1:remaining) {
      selected <- sample(1:num_groups, 1)
      base_sizes[selected] <- base_sizes[selected] + 1
    }
  }
  
  # Create group assignments
  group <- rep(1:num_groups, base_sizes)
  
  # Actual sample size might be slightly different due to rounding
  n_actual <- length(group)
  
  # Simulate covariates
  X <- matrix(runif(n_actual * p), n_actual, p)
  colnames(X) <- c("X1", "X2")
  
  # Define the expit (logistic) function
  expit <- function(x) { 1 / (1 + exp(-x)) }
  
  # Generate treatment indicator Z as a function of the confounders
  p_treat <- expit(0.3 * X[, "X1"] - 0.2 * X[, "X2"] + 0.5 * X[, "X1"] * X[, "X2"])
  Z <- rbinom(n_actual, 1, p_treat)
  
  # Create adjacency matrix for the clusters (chain structure)
  W <- matrix(0, nrow = num_groups, ncol = num_groups)
  for(i in 1:(num_groups - 1)) {
    W[i, i+1] <- 1
    W[i+1, i] <- 1
  }
  
  # Generate CAR random effects for clusters
  D_mat <- diag(rowSums(W))
  true_r_vec <- as.vector(rmvnorm(1, mean = rep(0, num_groups),
                                  sigma = sigma_r2_true * solve(D_mat - rho_car * W)))
  
  # Define the true function - with treatment effect that varies with covariates
  f_true <- function(X, Z) {
    sin(pi * X[, "X1"]) + log(1 + X[, "X2"]^2) +
      2 * Z * (X[, "X1"] * X[, "X2"]) + (X[, "X1"]^2) * Z
  }
  
  # Generate log survival times WITHOUT random effect interaction
  # This is the key change from the original code - random effect is simply additive
  sigma_true <- sqrt(true_sigma2)
  
  # Generate true potential outcomes for both treatment options
  # T(0): potential outcome under control
  logT0 <- f_true(X, rep(0, n_actual)) + true_r_vec[group] + rnorm(n_actual, 0, sigma_true)
  T0 <- exp(logT0)
  
  # T(1): potential outcome under treatment
  logT1 <- f_true(X, rep(1, n_actual)) + true_r_vec[group] + rnorm(n_actual, 0, sigma_true)
  T1 <- exp(logT1)
  
  # Generate factual outcomes based on treatment assignment
  logT <- Z * logT1 + (1-Z) * logT0
  T_true <- exp(logT)
  
  # Generate censoring times
  C <- rexp(n_actual, rate = lambda_c)
  
  # Observed survival time and censoring indicator
  Y_obs <- pmin(T_true, C)
  delta <- as.integer(T_true <= C)
  
  # Calculate true marginal causal effects
  true_ACE <- mean(T1 - T0)
  true_RACE <- sapply(seq(1, 20, by = 1), function(t_star) {
    mean(pmin(T1, t_star) - pmin(T0, t_star))
  })
  true_SPCE <- sapply(seq(1, 20, by = 1), function(t) {
    mean(T1 > t) - mean(T0 > t)
  })
  
  # Build a data frame
  sim_data <- data.frame(
    cluster = group,
    X1 = X[, "X1"],
    X2 = X[, "X2"],
    Z  = Z,
    Y  = Y_obs,
    delta = delta,
    T0 = T0,  # Store potential outcomes for later causal effect calculation
    T1 = T1
  )
  
  # Return all simulation objects
  return(list(
    sim_data = sim_data,
    X = X,
    Z = Z,
    group = group,
    W = W,
    cluster_sizes = table(group),  # Return the actual cluster sizes
    true_r_vec = true_r_vec,
    f_true = f_true,
    true_sigma2 = true_sigma2,
    true_ACE = true_ACE,
    true_RACE = true_RACE,
    true_SPCE = true_SPCE,
    true_times = T_true
  ))
}

# -------------------------------------------------------------------------
# 2. Function to Estimate Marginal Causal Effects from MCMC Samples
# -------------------------------------------------------------------------

# Function to estimate the marginal SPCE
get_marginal_SPCE <- function(t, X, group, fit) {
  n <- nrow(X)
  M <- nrow(fit$r)  # Number of MCMC samples
  
  # Initialize storage for survival probabilities
  S1_samples <- matrix(0, nrow = M, ncol = n)
  S0_samples <- matrix(0, nrow = M, ncol = n)
  
  # Obtain propensity scores for all observations
  ps_raw <- fit$propensity_forest$do_predict(X)
  e_hat <- 1 / (1 + exp(-ps_raw))
  
  # For each MCMC sample
  for (m in 1:M) {
    # For each observation
    for (i in 1:n) {
      # Get the cluster for this observation
      g <- group[i]
      
      # Create augmented covariate matrices for treatment = 1 and 0
      X_aug1 <- matrix(c(X[i,], 1, e_hat[i]), nrow = 1)
      X_aug0 <- matrix(c(X[i,], 0, e_hat[i]), nrow = 1)
      
      # Predict f(X) for treatment = 1 and 0
      f1 <- fit$forest$do_predict(X_aug1)
      f0 <- fit$forest$do_predict(X_aug0)
      
      # Get the random effect for this cluster
      r_g <- fit$r[m, g]
      
      # Calculate survival probabilities
      mu1 <- f1 + r_g
      mu0 <- f0 + r_g
      sigma <- sqrt(fit$sigma2[m])
      
      S1_samples[m, i] <- 1 - pnorm((log(t) - mu1) / sigma)
      S0_samples[m, i] <- 1 - pnorm((log(t) - mu0) / sigma)
    }
  }
  
  # Calculate marginal survival for each MCMC sample
  S1_marginal <- rowMeans(S1_samples)
  S0_marginal <- rowMeans(S0_samples)
  
  # Calculate SPCE for each MCMC sample
  SPCE_samples <- S1_marginal - S0_marginal
  
  # Return posterior mean and credible interval
  return(list(
    estimate = mean(SPCE_samples),
    CI = quantile(SPCE_samples, probs = c(0.025, 0.975))
  ))
}

# Function to estimate the marginal ACE
get_marginal_ACE <- function(X, group, fit) {
  n <- nrow(X)
  M <- nrow(fit$r)  # Number of MCMC samples
  
  # Initialize storage for expected survival times
  ET1_samples <- matrix(0, nrow = M, ncol = n)
  ET0_samples <- matrix(0, nrow = M, ncol = n)
  
  # Obtain propensity scores for all observations
  ps_raw <- fit$propensity_forest$do_predict(X)
  e_hat <- 1 / (1 + exp(-ps_raw))
  
  # For each MCMC sample
  for (m in 1:M) {
    # For each observation
    for (i in 1:n) {
      # Get the cluster for this observation
      g <- group[i]
      
      # Create augmented covariate matrices for treatment = 1 and 0
      X_aug1 <- matrix(c(X[i,], 1, e_hat[i]), nrow = 1)
      X_aug0 <- matrix(c(X[i,], 0, e_hat[i]), nrow = 1)
      
      # Predict f(X) for treatment = 1 and 0
      f1 <- fit$forest$do_predict(X_aug1)
      f0 <- fit$forest$do_predict(X_aug0)
      
      # Get the random effect for this cluster
      r_g <- fit$r[m, g]
      
      # Calculate expected survival times (formula for lognormal)
      mu1 <- f1 + r_g
      mu0 <- f0 + r_g
      sigma2 <- fit$sigma2[m]
      
      ET1_samples[m, i] <- exp(mu1 + sigma2/2)
      ET0_samples[m, i] <- exp(mu0 + sigma2/2)
    }
  }
  
  # Calculate marginal expected survival for each MCMC sample
  ET1_marginal <- rowMeans(ET1_samples)
  ET0_marginal <- rowMeans(ET0_samples)
  
  # Calculate ACE for each MCMC sample
  ACE_samples <- ET1_marginal - ET0_marginal
  
  # Return posterior mean and credible interval
  return(list(
    estimate = mean(ACE_samples),
    CI = quantile(ACE_samples, probs = c(0.025, 0.975))
  ))
}

# Function to estimate the marginal RACE
get_marginal_RACE <- function(t_star, X, group, fit) {
  n <- nrow(X)
  M <- nrow(fit$r)  # Number of MCMC samples
  
  # Initialize storage for restricted expected survival times
  ERT1_samples <- matrix(0, nrow = M, ncol = n)
  ERT0_samples <- matrix(0, nrow = M, ncol = n)
  
  # Obtain propensity scores for all observations
  ps_raw <- fit$propensity_forest$do_predict(X)
  e_hat <- 1 / (1 + exp(-ps_raw))
  
  # For each MCMC sample
  for (m in 1:M) {
    # For each observation
    for (i in 1:n) {
      # Get the cluster for this observation
      g <- group[i]
      
      # Create augmented covariate matrices for treatment = 1 and 0
      X_aug1 <- matrix(c(X[i,], 1, e_hat[i]), nrow = 1)
      X_aug0 <- matrix(c(X[i,], 0, e_hat[i]), nrow = 1)
      
      # Predict f(X) for treatment = 1 and 0
      f1 <- fit$forest$do_predict(X_aug1)
      f0 <- fit$forest$do_predict(X_aug0)
      
      # Get the random effect for this cluster
      r_g <- fit$r[m, g]
      
      # Calculate parameters
      mu1 <- f1 + r_g
      mu0 <- f0 + r_g
      sigma <- sqrt(fit$sigma2[m])
      sigma2 <- fit$sigma2[m]
      
      # For lognormal distribution:
      # E[min(T,t*)] = exp(mu + sigma2/2) * Φ((log(t*) - mu - sigma2)/sigma) + 
      #                t* * (1 - Φ((log(t*) - mu)/sigma))
      ERT1_samples[m, i] <- exp(mu1 + sigma2/2) * 
        pnorm((log(t_star) - mu1 - sigma2) / sigma) +
        t_star * (1 - pnorm((log(t_star) - mu1) / sigma))
      
      ERT0_samples[m, i] <- exp(mu0 + sigma2/2) * 
        pnorm((log(t_star) - mu0 - sigma2) / sigma) +
        t_star * (1 - pnorm((log(t_star) - mu0) / sigma))
    }
  }
  
  # Calculate marginal restricted expected survival for each MCMC sample
  ERT1_marginal <- rowMeans(ERT1_samples)
  ERT0_marginal <- rowMeans(ERT0_samples)
  
  # Calculate RACE for each MCMC sample
  RACE_samples <- ERT1_marginal - ERT0_marginal
  
  # Return posterior mean and credible interval
  return(list(
    estimate = mean(RACE_samples),
    CI = quantile(RACE_samples, probs = c(0.025, 0.975))
  ))
}

# -------------------------------------------------------------------------
# NEW: Cox Frailty Model Causal Inference Functions
# -------------------------------------------------------------------------

# Function to fit a Cox frailty model



# Function to fit a Cox frailty model - simplified approach
fit_cox_frailty <- function(data, formula = NULL) {
  # Try different model formulations from most to least complex
  
  # First, try a model with main effects and just treatment interactions
  # (without the three-way interaction which may be causing issues)
  model1_formula <- Surv(Y, delta) ~ Z + X1 + X2 + X1:Z + X2:Z + frailty(cluster)
  model1 <- try(coxph(model1_formula, data = data), silent = TRUE)
  if (!inherits(model1, "try-error")) {
    # Check if model converged properly
    if (!is.null(model1$coefficients) && !any(is.na(model1$coefficients))) {
      return(list(model = model1, type = "coxph", formula = model1_formula))
    }
  }
  
  # If that fails, try an even simpler model with just main effects and treatment
  model2_formula <- Surv(Y, delta) ~ Z + X1 + X2 + frailty(cluster)
  model2 <- try(coxph(model2_formula, data = data), silent = TRUE)
  if (!inherits(model2, "try-error")) {
    if (!is.null(model2$coefficients) && !any(is.na(model2$coefficients))) {
      return(list(model = model2, type = "coxph", formula = model2_formula))
    }
  }
  
  # If frailty still causes issues, try a robust but simplified model without frailty
  # This at least captures the cluster structure even if not through frailty
  model3_formula <- Surv(Y, delta) ~ Z + X1 + X2 + X1:Z + X2:Z + cluster(cluster)
  model3 <- try(coxph(model3_formula, data = data), silent = TRUE)
  if (!inherits(model3, "try-error")) {
    if (!is.null(model3$coefficients) && !any(is.na(model3$coefficients))) {
      return(list(model = model3, type = "coxph_cluster", formula = model3_formula))
    }
  }
  
  # Last resort: standard Cox model with main effects and key interactions
  model4_formula <- Surv(Y, delta) ~ Z + X1 + X2 + X1:Z + X2:Z
  model4 <- try(coxph(model4_formula, data = data), silent = TRUE)
  if (!inherits(model4, "try-error")) {
    if (!is.null(model4$coefficients) && !any(is.na(model4$coefficients))) {
      return(list(model = model4, type = "standard_cox", formula = model4_formula))
    }
  }
  
  # If everything else fails, use the most basic model that should converge
  model5_formula <- Surv(Y, delta) ~ Z
  model5 <- coxph(model5_formula, data = data)
  return(list(model = model5, type = "minimal_cox", formula = model5_formula))
}

# Function to estimate survival probability from Cox model
get_cox_survival_prob <- function(fit_obj, data, newdata_z0, newdata_z1, time_points) {
  fit <- fit_obj$model
  model_type <- fit_obj$type
  
  # Extract survival probabilities at specific time points
  surv_z0 <- try(survfit(fit, newdata = newdata_z0), silent = TRUE)
  surv_z1 <- try(survfit(fit, newdata = newdata_z1), silent = TRUE)
  
  # If prediction fails, return NA values
  if (inherits(surv_z0, "try-error") || inherits(surv_z1, "try-error")) {
    return(list(
      time = time_points,
      s0 = rep(NA, length(time_points)),
      s1 = rep(NA, length(time_points)),
      spce = rep(NA, length(time_points))
    ))
  }
  
  # Get actual times from the survival objects
  times_z0 <- surv_z0$time
  times_z1 <- surv_z1$time
  
  # Get survival probabilities at the closest available times
  s0 <- numeric(length(time_points))
  s1 <- numeric(length(time_points))
  
  for (i in 1:length(time_points)) {
    t <- time_points[i]
    
    # For control group
    if (t <= min(times_z0)) {
      s0[i] <- 1  # Before first event, survival is 1
    } else if (t > max(times_z0)) {
      s0[i] <- min(surv_z0$surv, na.rm = TRUE)  # After last event, use last probability
    } else {
      # Find closest time point at or before t
      idx <- max(which(times_z0 <= t))
      s0[i] <- surv_z0$surv[idx]
    }
    
    # For treatment group
    if (t <= min(times_z1)) {
      s1[i] <- 1
    } else if (t > max(times_z1)) {
      s1[i] <- min(surv_z1$surv, na.rm = TRUE)
    } else {
      idx <- max(which(times_z1 <= t))
      s1[i] <- surv_z1$surv[idx]
    }
  }
  
  # Calculate SPCE
  spce <- s1 - s0
  
  return(list(
    time = time_points,
    s0 = s0,
    s1 = s1,
    spce = spce
  ))
}

# Function to estimate average survival time (ACE) from Cox model
get_cox_ace <- function(fit_obj, data, newdata_z0, newdata_z1, max_time = 15) {
  fit <- fit_obj$model
  model_type <- fit_obj$type
  
  # Get survival curves
  surv_z0 <- try(survfit(fit, newdata = newdata_z0), silent = TRUE)
  surv_z1 <- try(survfit(fit, newdata = newdata_z1), silent = TRUE)
  
  if (inherits(surv_z0, "try-error") || inherits(surv_z1, "try-error")) {
    return(NA)
  }
  
  # Find the max observed time in the predictions
  max_obs_time <- min(max(surv_z0$time), max(surv_z1$time), max_time)
  
  # Calculate RMST using numerical integration (trapezoidal rule)
  # This is more robust than using summary(survfit, rmean=) which can sometimes fail
  
  # For control group
  time_points_z0 <- c(0, surv_z0$time[surv_z0$time <= max_obs_time])
  surv_probs_z0 <- c(1, surv_z0$surv[surv_z0$time <= max_obs_time])
  if (length(time_points_z0) < 2) {
    # If insufficient time points, approximate
    rmst_z0 <- max_obs_time  # Assume survival = 1 for all time
  } else {
    # Trapezoidal integration
    rmst_z0 <- 0
    for (i in 1:(length(time_points_z0)-1)) {
      dt <- time_points_z0[i+1] - time_points_z0[i]
      avg_surv <- (surv_probs_z0[i] + surv_probs_z0[i+1]) / 2
      rmst_z0 <- rmst_z0 + dt * avg_surv
    }
  }
  
  # For treatment group
  time_points_z1 <- c(0, surv_z1$time[surv_z1$time <= max_obs_time])
  surv_probs_z1 <- c(1, surv_z1$surv[surv_z1$time <= max_obs_time])
  if (length(time_points_z1) < 2) {
    rmst_z1 <- max_obs_time
  } else {
    rmst_z1 <- 0
    for (i in 1:(length(time_points_z1)-1)) {
      dt <- time_points_z1[i+1] - time_points_z1[i]
      avg_surv <- (surv_probs_z1[i] + surv_probs_z1[i+1]) / 2
      rmst_z1 <- rmst_z1 + dt * avg_surv
    }
  }
  
  # Calculate ACE
  ace <- rmst_z1 - rmst_z0
  
  return(ace)
}

# Function to estimate RACE from Cox model
get_cox_race <- function(fit_obj, data, newdata_z0, newdata_z1, time_points) {
  fit <- fit_obj$model
  model_type <- fit_obj$type
  
  # Get survival curves once
  surv_z0 <- try(survfit(fit, newdata = newdata_z0), silent = TRUE)
  surv_z1 <- try(survfit(fit, newdata = newdata_z1), silent = TRUE)
  
  if (inherits(surv_z0, "try-error") || inherits(surv_z1, "try-error")) {
    return(list(
      time = time_points,
      race = rep(NA, length(time_points))
    ))
  }
  
  # Initialize storage for RACE estimates
  race_estimates <- numeric(length(time_points))
  
  # For each time point, calculate RMST up to that point
  for (i in 1:length(time_points)) {
    t_star <- time_points[i]
    
    # For control group
    time_points_z0 <- c(0, surv_z0$time[surv_z0$time <= t_star])
    surv_probs_z0 <- c(1, surv_z0$surv[surv_z0$time <= t_star])
    if (length(time_points_z0) < 2) {
      # If no events before t_star, rmst = t_star
      rmst_z0 <- t_star
    } else {
      # Trapezoidal integration
      rmst_z0 <- 0
      for (j in 1:(length(time_points_z0)-1)) {
        dt <- time_points_z0[j+1] - time_points_z0[j]
        avg_surv <- (surv_probs_z0[j] + surv_probs_z0[j+1]) / 2
        rmst_z0 <- rmst_z0 + dt * avg_surv
      }
    }
    
    # For treatment group
    time_points_z1 <- c(0, surv_z1$time[surv_z1$time <= t_star])
    surv_probs_z1 <- c(1, surv_z1$surv[surv_z1$time <= t_star])
    if (length(time_points_z1) < 2) {
      rmst_z1 <- t_star
    } else {
      rmst_z1 <- 0
      for (j in 1:(length(time_points_z1)-1)) {
        dt <- time_points_z1[j+1] - time_points_z1[j]
        avg_surv <- (surv_probs_z1[j] + surv_probs_z1[j+1]) / 2
        rmst_z1 <- rmst_z1 + dt * avg_surv
      }
    }
    
    race_estimates[i] <- rmst_z1 - rmst_z0
  }
  
  return(list(
    time = time_points,
    race = race_estimates
  ))
}

# Improved bootstrap function with more robust error handling
bootstrap_cox_ci <- function(data, time_points, n_bootstrap = 200) {
  # Initialize storage
  boot_spce <- matrix(NA, nrow = n_bootstrap, ncol = length(time_points))
  boot_race <- matrix(NA, nrow = n_bootstrap, ncol = length(time_points))
  boot_ace <- numeric(n_bootstrap)
  
  # Create counterfactual datasets for prediction
  data_all_z0 <- data
  data_all_z0$Z <- 0
  data_all_z1 <- data
  data_all_z1$Z <- 1
  
  # Get the unique clusters
  clusters <- unique(data$cluster)
  
  # Bootstrap by resampling clusters
  successful_boots <- 0
  cat("Starting bootstrap iterations...\n")
  
  for (b in 1:n_bootstrap) {
    if (b %% 10 == 0) {
      cat("Bootstrap iteration:", b, "\n")
    }
    
    # Sample clusters with replacement
    sampled_clusters <- sample(clusters, length(clusters), replace = TRUE)
    
    # Create bootstrap sample
    boot_indices <- unlist(lapply(sampled_clusters, function(c) which(data$cluster == c)))
    boot_data <- data[boot_indices, ]
    
    # Fit Cox model and calculate estimates
    tryCatch({
      # Fit Cox model
      boot_fit <- fit_cox_frailty(boot_data)
      
      # Calculate SPCE
      spce_result <- get_cox_survival_prob(boot_fit, boot_data, data_all_z0, data_all_z1, time_points)
      boot_spce[b, ] <- spce_result$spce
      
      # Calculate RACE
      race_result <- get_cox_race(boot_fit, boot_data, data_all_z0, data_all_z1, time_points)
      boot_race[b, ] <- race_result$race
      
      # Calculate ACE
      boot_ace[b] <- get_cox_ace(boot_fit, boot_data, data_all_z0, data_all_z1, max(time_points))
      
      successful_boots <- successful_boots + 1
    }, error = function(e) {
      if (b %% 10 == 0) {
        cat("Error in bootstrap iteration", b, ":", e$message, "\n")
      }
    })
  }
  
  cat("Successfully completed", successful_boots, "out of", n_bootstrap, "bootstrap iterations\n")
  
  # Calculate percentile confidence intervals
  spce_lower <- apply(boot_spce, 2, function(x) quantile(x, 0.025, na.rm = TRUE))
  spce_upper <- apply(boot_spce, 2, function(x) quantile(x, 0.975, na.rm = TRUE))
  
  race_lower <- apply(boot_race, 2, function(x) quantile(x, 0.025, na.rm = TRUE))
  race_upper <- apply(boot_race, 2, function(x) quantile(x, 0.975, na.rm = TRUE))
  
  ace_ci <- quantile(boot_ace, c(0.025, 0.975), na.rm = TRUE)
  
  return(list(
    spce_ci = data.frame(time = time_points, lower = spce_lower, upper = spce_upper),
    race_ci = data.frame(time = time_points, lower = race_lower, upper = race_upper),
    ace_ci = ace_ci,
    boot_samples = list(spce = boot_spce, race = boot_race, ace = boot_ace),
    successful_boots = successful_boots
  ))
}
# -------------------------------------------------------------------------
# 3. Single Replicate Analysis Function with Cox Frailty Comparison
# -------------------------------------------------------------------------

run_single_replicate_with_cox <- function(seed = 123, 
                                          n = 2000, 
                                          n_iter = 2000, 
                                          burn_in = 500, 
                                          thin = 5,
                                          time_grid = seq(1, 15, by = 1),
                                          min_cluster_size = 50,
                                          max_cluster_size = 500,
                                          num_groups = 5,
                                          n_bootstrap = 100) {
  
  # 1. Simulate data with variable cluster sizes
  sim_data <- simulate_causal_data_variable_clusters(
    seed = seed, 
    n = n,
    num_groups = num_groups,
    min_cluster_size = min_cluster_size,
    max_cluster_size = max_cluster_size
  )
  
  # 2. SOFTBART Method: Fit the model
  cat("Fitting SoftBART model...\n")
  fit <- AFT_mixed_CAR_softbart_causal(
    time = sim_data$sim_data$Y,
    status = sim_data$sim_data$delta,
    X = sim_data$X,
    Z = sim_data$Z,
    group = sim_data$group,
    W = sim_data$W,
    X_test = sim_data$X,
    Z_test = sim_data$Z,
    group_test = sim_data$group,
    n_iter = n_iter, 
    burn_in = burn_in, 
    thin = thin
  )
  
  # 3. Compute SoftBART marginal causal effects
  
  # Marginal SPCE at each time point
  marginal_SPCE <- lapply(time_grid, function(t) {
    get_marginal_SPCE(t, sim_data$X, sim_data$group, fit)
  })
  
  # Marginal ACE
  marginal_ACE <- get_marginal_ACE(sim_data$X, sim_data$group, fit)
  
  # Marginal RACE at each time point
  marginal_RACE <- lapply(time_grid, function(t_star) {
    get_marginal_RACE(t_star, sim_data$X, sim_data$group, fit)
  })
  
  # 4. COX FRAILTY Method: Fit Cox model
  cat("Fitting Cox frailty model...\n")
  cox_fit <- fit_cox_frailty(sim_data$sim_data)
  
  # Create counterfactual datasets for Cox predictions
  data_all_z0 <- sim_data$sim_data
  data_all_z0$Z <- 0
  data_all_z1 <- sim_data$sim_data
  data_all_z1$Z <- 1
  
  # Get Cox SPCE estimates
  cox_spce <- get_cox_survival_prob(cox_fit, sim_data$sim_data, data_all_z0, data_all_z1, time_grid)
  
  # Get Cox ACE estimate
  cox_ace <- get_cox_ace(cox_fit, sim_data$sim_data, data_all_z0, data_all_z1, max(time_grid))
  
  # Get Cox RACE estimates
  cox_race <- get_cox_race(cox_fit, sim_data$sim_data, data_all_z0, data_all_z1, time_grid)
  
  # 5. Bootstrap CI for Cox estimates
  cat("Bootstrapping Cox confidence intervals...\n")
  cox_bootstrap <- bootstrap_cox_ci(sim_data$sim_data, time_grid, n_bootstrap)
  
  # 6. Compute conditional causal effects (by cluster) - SoftBART
  
  # For each cluster, construct a representative subject (using mean values)
  clusters <- sort(unique(sim_data$group))
  cluster_results <- list()
  
  for (i in clusters) {
    # Get data for this cluster
    cluster_data <- sim_data$sim_data[sim_data$sim_data$cluster == i, ]
    
    # Compute mean covariates for the cluster
    X_i <- matrix(colMeans(cluster_data[, c("X1", "X2")]), nrow = 1)
    colnames(X_i) <- c("X1", "X2")
    
    # Get true and estimated cluster-specific random effects
    r_true_i <- sim_data$true_r_vec[i]
    r_est_i <- mean(fit$r[, i])
    
    # CSPCE at each time point
    cspce_results <- lapply(time_grid, function(t) {
      getCSPCE(t = t, X = X_i, r = r_est_i, fit = fit)
    })
    
    # CACE 
    cace_result <- getCACE(X = X_i, r = r_est_i, fit = fit)
    
    # CRACE at each time point
    crace_results <- lapply(time_grid, function(t_star) {
      getCRACE(t_star = t_star, X = X_i, r = r_est_i, fit = fit)
    })
    
    # Store results for this cluster
    cluster_results[[as.character(i)]] <- list(
      X_i = X_i,
      r_true = r_true_i,
      r_est = r_est_i,
      CSPCE = cspce_results,
      CACE = cace_result,
      CRACE = crace_results,
      cluster_size = nrow(cluster_data)
    )
  }
  
  # 7. Compile comparison results
  
  # Create data frame for SPCE comparison
  spce_comparison <- data.frame(
    Time = time_grid,
    True_SPCE = sapply(time_grid, function(t) {
      sim_data$true_SPCE[which(seq(1,20) == t)]
    }),
    SoftBART_SPCE = sapply(marginal_SPCE, function(x) x$estimate),
    SoftBART_SPCE_Lower = sapply(marginal_SPCE, function(x) x$CI[1]),
    SoftBART_SPCE_Upper = sapply(marginal_SPCE, function(x) x$CI[2]),
    Cox_SPCE = cox_spce$spce,
    Cox_SPCE_Lower = cox_bootstrap$spce_ci$lower,
    Cox_SPCE_Upper = cox_bootstrap$spce_ci$upper
  )
  
  # Calculate coverage
  spce_comparison$SoftBART_SPCE_Coverage <- 
    (spce_comparison$True_SPCE >= spce_comparison$SoftBART_SPCE_Lower & 
       spce_comparison$True_SPCE <= spce_comparison$SoftBART_SPCE_Upper)
  
  spce_comparison$Cox_SPCE_Coverage <- 
    (spce_comparison$True_SPCE >= spce_comparison$Cox_SPCE_Lower & 
       spce_comparison$True_SPCE <= spce_comparison$Cox_SPCE_Upper)
  
  # Create data frame for RACE comparison
  race_comparison <- data.frame(
    Time = time_grid,
    True_RACE = sapply(time_grid, function(t) {
      sim_data$true_RACE[which(seq(1,20) == t)]
    }),
    SoftBART_RACE = sapply(marginal_RACE, function(x) x$estimate),
    SoftBART_RACE_Lower = sapply(marginal_RACE, function(x) x$CI[1]),
    SoftBART_RACE_Upper = sapply(marginal_RACE, function(x) x$CI[2]),
    Cox_RACE = cox_race$race,
    Cox_RACE_Lower = cox_bootstrap$race_ci$lower,
    Cox_RACE_Upper = cox_bootstrap$race_ci$upper
  )
  
  # Calculate coverage
  race_comparison$SoftBART_RACE_Coverage <- 
    (race_comparison$True_RACE >= race_comparison$SoftBART_RACE_Lower & 
       race_comparison$True_RACE <= race_comparison$SoftBART_RACE_Upper)
  
  race_comparison$Cox_RACE_Coverage <- 
    (race_comparison$True_RACE >= race_comparison$Cox_RACE_Lower & 
       race_comparison$True_RACE <= race_comparison$Cox_RACE_Upper)
  
  # Create data frame for ACE comparison
  ace_comparison <- data.frame(
    True_ACE = sim_data$true_ACE,
    SoftBART_ACE = marginal_ACE$estimate,
    SoftBART_ACE_Lower = marginal_ACE$CI[1],
    SoftBART_ACE_Upper = marginal_ACE$CI[2],
    Cox_ACE = cox_ace,
    Cox_ACE_Lower = cox_bootstrap$ace_ci[1],
    Cox_ACE_Upper = cox_bootstrap$ace_ci[2]
  )
  
  # Calculate coverage
  ace_comparison$SoftBART_ACE_Coverage <- 
    (ace_comparison$True_ACE >= ace_comparison$SoftBART_ACE_Lower & 
       ace_comparison$True_ACE <= ace_comparison$SoftBART_ACE_Upper)
  
  ace_comparison$Cox_ACE_Coverage <- 
    (ace_comparison$True_ACE >= ace_comparison$Cox_ACE_Lower & 
       ace_comparison$True_ACE <= ace_comparison$Cox_ACE_Upper)
  
  # Return all results including cluster sizes
  return(list(
    # Simulation data
    seed = seed,
    sim_data = sim_data,
    
    # SoftBART results
    softbart_fit = fit,
    softbart_marginal_SPCE = marginal_SPCE,
    softbart_marginal_ACE = marginal_ACE,
    softbart_marginal_RACE = marginal_RACE,
    softbart_cluster_results = cluster_results,
    
    # Cox results
    cox_fit = cox_fit,
    cox_spce = cox_spce,
    cox_ace = cox_ace,
    cox_race = cox_race,
    cox_bootstrap = cox_bootstrap,
    
    # Comparison results
    spce_comparison = spce_comparison,
    race_comparison = race_comparison,
    ace_comparison = ace_comparison,
    
    # Misc
    true_values = list(
      ACE = sim_data$true_ACE,
      RACE = sim_data$true_RACE,
      SPCE = sim_data$true_SPCE
    ),
    time_grid = time_grid,
    cluster_sizes = sim_data$cluster_sizes
  ))
}

# -------------------------------------------------------------------------
# 4. Multiple Replicate Simulation Function with Comparison
# -------------------------------------------------------------------------

run_simulation_study_with_comparison <- function(n_replicates = 10,
                                                 n = 1000, 
                                                 n_iter = 2000,
                                                 burn_in = 500,
                                                 thin = 5,
                                                 time_grid = seq(1, 15, by = 1),
                                                 min_cluster_size = 50,
                                                 max_cluster_size = 500,
                                                 num_groups = 5,
                                                 n_bootstrap = 100,
                                                 use_parallel = FALSE,
                                                 n_cores = NULL) {
  
  if (use_parallel) {
    if (is.null(n_cores)) {
      n_cores <- min(parallel::detectCores() - 1, n_replicates)
    }
    
    # Create a cluster
    cl <- makeCluster(n_cores)
    
    # Export necessary functions and variables
    clusterExport(cl, c("run_single_replicate_with_cox", 
                        "simulate_causal_data_variable_clusters", 
                        "get_marginal_SPCE", "get_marginal_ACE", "get_marginal_RACE",
                        "getCSPCE", "getCACE", "getCRACE", "AFT_mixed_CAR_softbart_causal",
                        "fit_cox_frailty", "get_cox_survival_prob", "get_cox_ace", 
                        "get_cox_race", "bootstrap_cox_ci",
                        "n", "n_iter", "burn_in", "thin", "time_grid", 
                        "min_cluster_size", "max_cluster_size", "num_groups", "n_bootstrap"))
    
    # Load libraries on each cluster
    clusterEvalQ(cl, {
      library(mvtnorm)
      library(SoftBart)
      library(dplyr)
      library(survival)
      library(coxme)
    })
    
    # Run replicates in parallel
    replicate_seeds <- 1:n_replicates
    
    all_replicates <- parLapply(cl, replicate_seeds, function(seed) {
      run_single_replicate_with_cox(
        seed = seed, 
        n = n, 
        n_iter = n_iter, 
        burn_in = burn_in, 
        thin = thin, 
        time_grid = time_grid,
        min_cluster_size = min_cluster_size, 
        max_cluster_size = max_cluster_size,
        num_groups = num_groups,
        n_bootstrap = n_bootstrap
      )
    })
    
    # Stop the cluster
    stopCluster(cl)
    
  } else {
    # Run replicates sequentially
    all_replicates <- list()
    
    for (rep_idx in 1:n_replicates) {
      message(paste("Running replicate", rep_idx, "of", n_replicates))
      all_replicates[[rep_idx]] <- run_single_replicate_with_cox(
        seed = rep_idx, 
        n = n, 
        n_iter = n_iter,
        burn_in = burn_in,
        thin = thin,
        time_grid = time_grid,
        min_cluster_size = min_cluster_size, 
        max_cluster_size = max_cluster_size,
        num_groups = num_groups,
        n_bootstrap = n_bootstrap
      )
    }
  }
  
  # Aggregate results across replicates
  
  # 1. Initialize aggregate data frames
  spce_comparison_by_time <- data.frame(
    Time = rep(time_grid, n_replicates),
    Replicate = rep(1:n_replicates, each = length(time_grid)),
    True_SPCE = NA,
    SoftBART_SPCE = NA,
    SoftBART_SPCE_Coverage = NA,
    Cox_SPCE = NA,
    Cox_SPCE_Coverage = NA
  )
  
  race_comparison_by_time <- data.frame(
    Time = rep(time_grid, n_replicates),
    Replicate = rep(1:n_replicates, each = length(time_grid)),
    True_RACE = NA,
    SoftBART_RACE = NA,
    SoftBART_RACE_Coverage = NA,
    Cox_RACE = NA,
    Cox_RACE_Coverage = NA
  )
  
  ace_comparison_df <- data.frame(
    Replicate = 1:n_replicates,
    True_ACE = NA,
    SoftBART_ACE = NA,
    SoftBART_ACE_Coverage = NA,
    Cox_ACE = NA,
    Cox_ACE_Coverage = NA
  )
  
  # 2. Fill in the data frames
  for (rep_idx in 1:n_replicates) {
    rep_results <- all_replicates[[rep_idx]]
    
    # Fill SPCE comparison
    for (t_idx in 1:length(time_grid)) {
      row_idx <- (rep_idx - 1) * length(time_grid) + t_idx
      spce_comparison_by_time[row_idx, "True_SPCE"] <- 
        rep_results$spce_comparison$True_SPCE[t_idx]
      spce_comparison_by_time[row_idx, "SoftBART_SPCE"] <- 
        rep_results$spce_comparison$SoftBART_SPCE[t_idx]
      spce_comparison_by_time[row_idx, "SoftBART_SPCE_Coverage"] <- 
        rep_results$spce_comparison$SoftBART_SPCE_Coverage[t_idx]
      spce_comparison_by_time[row_idx, "Cox_SPCE"] <- 
        rep_results$spce_comparison$Cox_SPCE[t_idx]
      spce_comparison_by_time[row_idx, "Cox_SPCE_Coverage"] <- 
        rep_results$spce_comparison$Cox_SPCE_Coverage[t_idx]
    }
    
    # Fill RACE comparison
    for (t_idx in 1:length(time_grid)) {
      row_idx <- (rep_idx - 1) * length(time_grid) + t_idx
      race_comparison_by_time[row_idx, "True_RACE"] <- 
        rep_results$race_comparison$True_RACE[t_idx]
      race_comparison_by_time[row_idx, "SoftBART_RACE"] <- 
        rep_results$race_comparison$SoftBART_RACE[t_idx]
      race_comparison_by_time[row_idx, "SoftBART_RACE_Coverage"] <- 
        rep_results$race_comparison$SoftBART_RACE_Coverage[t_idx]
      race_comparison_by_time[row_idx, "Cox_RACE"] <- 
        rep_results$race_comparison$Cox_RACE[t_idx]
      race_comparison_by_time[row_idx, "Cox_RACE_Coverage"] <- 
        rep_results$race_comparison$Cox_RACE_Coverage[t_idx]
    }
    
    # Fill ACE comparison
    ace_comparison_df[rep_idx, "True_ACE"] <- 
      rep_results$ace_comparison$True_ACE
    ace_comparison_df[rep_idx, "SoftBART_ACE"] <- 
      rep_results$ace_comparison$SoftBART_ACE
    ace_comparison_df[rep_idx, "SoftBART_ACE_Coverage"] <- 
      rep_results$ace_comparison$SoftBART_ACE_Coverage
    ace_comparison_df[rep_idx, "Cox_ACE"] <- 
      rep_results$ace_comparison$Cox_ACE
    ace_comparison_df[rep_idx, "Cox_ACE_Coverage"] <- 
      rep_results$ace_comparison$Cox_ACE_Coverage
  }
  
  # 3. Create summary tables
  
  # SPCE summary by time
  spce_summary <- spce_comparison_by_time %>%
    group_by(Time) %>%
    summarize(
      Mean_True_SPCE = mean(True_SPCE),
      Mean_SoftBART_SPCE = mean(SoftBART_SPCE),
      SoftBART_SPCE_Coverage = mean(SoftBART_SPCE_Coverage) * 100,
      SoftBART_SPCE_Bias = mean(SoftBART_SPCE - True_SPCE),
      SoftBART_SPCE_RMSE = sqrt(mean((SoftBART_SPCE - True_SPCE)^2)),
      Mean_Cox_SPCE = mean(Cox_SPCE),
      Cox_SPCE_Coverage = mean(Cox_SPCE_Coverage) * 100,
      Cox_SPCE_Bias = mean(Cox_SPCE - True_SPCE),
      Cox_SPCE_RMSE = sqrt(mean((Cox_SPCE - True_SPCE)^2))
    )
  
  # RACE summary by time
  race_summary <- race_comparison_by_time %>%
    group_by(Time) %>%
    summarize(
      Mean_True_RACE = mean(True_RACE),
      Mean_SoftBART_RACE = mean(SoftBART_RACE),
      SoftBART_RACE_Coverage = mean(SoftBART_RACE_Coverage) * 100,
      SoftBART_RACE_Bias = mean(SoftBART_RACE - True_RACE),
      SoftBART_RACE_RMSE = sqrt(mean((SoftBART_RACE - True_RACE)^2)),
      Mean_Cox_RACE = mean(Cox_RACE),
      Cox_RACE_Coverage = mean(Cox_RACE_Coverage) * 100,
      Cox_RACE_Bias = mean(Cox_RACE - True_RACE),
      Cox_RACE_RMSE = sqrt(mean((Cox_RACE - True_RACE)^2))
    )
  
  # ACE summary
  ace_summary <- ace_comparison_df %>%
    summarize(
      Mean_True_ACE = mean(True_ACE),
      Mean_SoftBART_ACE = mean(SoftBART_ACE),
      SoftBART_ACE_Coverage = mean(SoftBART_ACE_Coverage) * 100,
      SoftBART_ACE_Bias = mean(SoftBART_ACE - True_ACE),
      SoftBART_ACE_RMSE = sqrt(mean((SoftBART_ACE - True_ACE)^2)),
      Mean_Cox_ACE = mean(Cox_ACE),
      Cox_ACE_Coverage = mean(Cox_ACE_Coverage) * 100,
      Cox_ACE_Bias = mean(Cox_ACE - True_ACE),
      Cox_ACE_RMSE = sqrt(mean((Cox_ACE - True_ACE)^2))
    )
  
  # 4. Generate comparison plots
  
  # Coverage probability plot for SPCE
  spce_coverage_plot <- ggplot(spce_summary, aes(x = Time)) +
    geom_line(aes(y = SoftBART_SPCE_Coverage, color = "SoftBART"), size = 1) +
    geom_line(aes(y = Cox_SPCE_Coverage, color = "Cox Frailty"), size = 1) +
    geom_hline(yintercept = 95, linetype = "dashed", color = "red") +
    labs(title = "Coverage Probability for SPCE by Method",
         x = "Time", y = "Coverage (%)", color = "Method") +
    ylim(0, 100) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Coverage probability plot for RACE
  race_coverage_plot <- ggplot(race_summary, aes(x = Time)) +
    geom_line(aes(y = SoftBART_RACE_Coverage, color = "SoftBART"), size = 1) +
    geom_line(aes(y = Cox_RACE_Coverage, color = "Cox Frailty"), size = 1) +
    geom_hline(yintercept = 95, linetype = "dashed", color = "red") +
    labs(title = "Coverage Probability for RACE by Method",
         x = "Time", y = "Coverage (%)", color = "Method") +
    ylim(0, 100) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Bias plot for SPCE
  spce_bias_plot <- ggplot(spce_summary, aes(x = Time)) +
    geom_line(aes(y = SoftBART_SPCE_Bias, color = "SoftBART"), size = 1) +
    geom_line(aes(y = Cox_SPCE_Bias, color = "Cox Frailty"), size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    labs(title = "Bias in SPCE Estimation by Method",
         x = "Time", y = "Bias", color = "Method") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Bias plot for RACE
  race_bias_plot <- ggplot(race_summary, aes(x = Time)) +
    geom_line(aes(y = SoftBART_RACE_Bias, color = "SoftBART"), size = 1) +
    geom_line(aes(y = Cox_RACE_Bias, color = "Cox Frailty"), size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    labs(title = "Bias in RACE Estimation by Method",
         x = "Time", y = "Bias", color = "Method") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Point estimate comparison plot for SPCE
  spce_est_plot <- ggplot(spce_summary, aes(x = Time)) +
    geom_line(aes(y = Mean_True_SPCE, color = "True"), size = 1, linetype = "dashed") +
    geom_line(aes(y = Mean_SoftBART_SPCE, color = "SoftBART"), size = 1) +
    geom_line(aes(y = Mean_Cox_SPCE, color = "Cox Frailty"), size = 1) +
    labs(title = "Mean Estimated SPCE by Method",
         x = "Time", y = "SPCE", color = "Method") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Point estimate comparison plot for RACE
  race_est_plot <- ggplot(race_summary, aes(x = Time)) +
    geom_line(aes(y = Mean_True_RACE, color = "True"), size = 1, linetype = "dashed") +
    geom_line(aes(y = Mean_SoftBART_RACE, color = "SoftBART"), size = 1) +
    geom_line(aes(y = Mean_Cox_RACE, color = "Cox Frailty"), size = 1) +
    labs(title = "Mean Estimated RACE by Method",
         x = "Time", y = "RACE", color = "Method") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Bar plot for ACE comparison
  ace_bar_data <- data.frame(
    Method = c("True", "SoftBART", "Cox Frailty"),
    ACE = c(ace_summary$Mean_True_ACE, ace_summary$Mean_SoftBART_ACE, ace_summary$Mean_Cox_ACE),
    Coverage = c(NA, ace_summary$SoftBART_ACE_Coverage, ace_summary$Cox_ACE_Coverage)
  )
  
  ace_bar_plot <- ggplot(ace_bar_data, aes(x = Method, y = ACE, fill = Method)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
    geom_text(aes(label = paste0("Coverage: ", round(Coverage,1), "%"), y = ACE + 0.1),
              position = position_dodge(width = 0.9), na.rm = TRUE, size = 3.5) +
    labs(title = "Average Causal Effect (ACE) by Method",
         x = "", y = "ACE") +
    theme_minimal() +
    theme(legend.position = "none")
  
  # 5. Create a comprehensive summary table for all metrics
  summary_table <- data.frame(
    Metric = c("Average Causal Effect (ACE)", 
               "Restricted Average Causal Effect (RACE)", 
               "Survival Probability Causal Effect (SPCE)"),
    True_Value = c(ace_summary$Mean_True_ACE, 
                   mean(race_summary$Mean_True_RACE), 
                   mean(spce_summary$Mean_True_SPCE)),
    SoftBART_Est = c(ace_summary$Mean_SoftBART_ACE, 
                     mean(race_summary$Mean_SoftBART_RACE), 
                     mean(spce_summary$Mean_SoftBART_SPCE)),
    SoftBART_Bias = c(ace_summary$SoftBART_ACE_Bias, 
                      mean(race_summary$SoftBART_RACE_Bias), 
                      mean(spce_summary$SoftBART_SPCE_Bias)),
    SoftBART_RMSE = c(ace_summary$SoftBART_ACE_RMSE, 
                      mean(race_summary$SoftBART_RACE_RMSE), 
                      mean(spce_summary$SoftBART_SPCE_RMSE)),
    SoftBART_Coverage = c(ace_summary$SoftBART_ACE_Coverage, 
                          mean(race_summary$SoftBART_RACE_Coverage), 
                          mean(spce_summary$SoftBART_SPCE_Coverage)),
    Cox_Est = c(ace_summary$Mean_Cox_ACE, 
                mean(race_summary$Mean_Cox_RACE), 
                mean(spce_summary$Mean_Cox_SPCE)),
    Cox_Bias = c(ace_summary$Cox_ACE_Bias, 
                 mean(race_summary$Cox_RACE_Bias), 
                 mean(spce_summary$Cox_SPCE_Bias)),
    Cox_RMSE = c(ace_summary$Cox_ACE_RMSE, 
                 mean(race_summary$Cox_RACE_RMSE), 
                 mean(spce_summary$Cox_SPCE_RMSE)),
    Cox_Coverage = c(ace_summary$Cox_ACE_Coverage, 
                     mean(race_summary$Cox_RACE_Coverage), 
                     mean(spce_summary$Cox_SPCE_Coverage))
  )
  
  # Return all results and summaries
  return(list(
    # Raw data
    all_replicates = all_replicates,
    spce_comparison_by_time = spce_comparison_by_time,
    race_comparison_by_time = race_comparison_by_time,
    ace_comparison_df = ace_comparison_df,
    
    # Summary tables
    spce_summary = spce_summary,
    race_summary = race_summary,
    ace_summary = ace_summary,
    summary_table = summary_table,
    
    # Plots
    spce_coverage_plot = spce_coverage_plot,
    race_coverage_plot = race_coverage_plot,
    spce_bias_plot = spce_bias_plot,
    race_bias_plot = race_bias_plot,
    spce_est_plot = spce_est_plot,
    race_est_plot = race_est_plot,
    ace_bar_plot = ace_bar_plot,
    
    # Combined plots
    combined_coverage_plot = spce_coverage_plot + race_coverage_plot + 
      plot_layout(guides = "collect") & theme(legend.position = "bottom"),
    combined_bias_plot = spce_bias_plot + race_bias_plot + 
      plot_layout(guides = "collect") & theme(legend.position = "bottom"),
    
    # Study settings
    n_replicates = n_replicates,
    n = n,
    time_grid = time_grid,
    min_cluster_size = min_cluster_size,
    max_cluster_size = max_cluster_size,
    num_groups = num_groups
  ))
}

# -------------------------------------------------------------------------
# 5. Example Run with Comparison
# -------------------------------------------------------------------------

# Set simulation parameters
n_replicates <- 3   # Small number for testing
sample_size <- 2000    # Reduced sample size for quicker run
n_iter <- 2000         # Reduced for quicker run
burn_in <- 500       # Reduced for quicker run
thin <- 5
time_grid <- seq(1, 10, by = 1)  # Fewer time points for quicker bootstrap
num_groups <- 15
min_cluster_size <- 30
max_cluster_size <- 300
n_bootstrap <-  100   # Reduced for quicker run

# Run simulation with comparison
set.seed(123)
cat("Starting comparison simulation with SoftBART and Cox frailty models...\n")
start_time <- Sys.time()

# Single replicate for testing
test_comparison <- run_single_replicate_with_cox(
  seed = 1,
  n = sample_size, 
  n_iter = n_iter,
  burn_in = burn_in,
  thin = thin,
  time_grid = time_grid,
  min_cluster_size = min_cluster_size,
  max_cluster_size = max_cluster_size,
  num_groups = num_groups,
  n_bootstrap = n_bootstrap
)

# Check coverage comparison
print(test_comparison$spce_comparison[, c("Time", "True_SPCE", "SoftBART_SPCE_Coverage", "Cox_SPCE_Coverage")])

# Full simulation with comparison (commented out - uncomment to run)
 comparison_results <- run_simulation_study_with_comparison(
   n_replicates = n_replicates,
   n = sample_size, 
   n_iter = n_iter,
   burn_in = burn_in,
   thin = thin,
   time_grid = time_grid,
   min_cluster_size = min_cluster_size,
   max_cluster_size = max_cluster_size,
   num_groups = num_groups,
   n_bootstrap = n_bootstrap,
   use_parallel = FALSE  # Set to TRUE for parallel processing
 )

# Display execution time
end_time <- Sys.time()
cat("Comparison simulation completed in:", format(end_time - start_time), "\n")

# Save results
# saveRDS(comparison_results, "softbart_cox_comparison_results.rds")

# -------------------------------------------------------------------------
# 6. Create Summary Plots and Tables for the Comparison
# -------------------------------------------------------------------------

# Function to create a publication-quality table for model comparison
create_comparison_summary_table <- function(results) {
  # Extract summary data
  summary_df <- results$summary_table
  
  # Format the table
  kable(summary_df, format = "markdown", digits = 3,
        col.names = c("Causal Effect", "True Value", 
                      "SoftBART Est", "SoftBART Bias", "SoftBART RMSE", "SoftBART Coverage",
                      "Cox Est", "Cox Bias", "Cox RMSE", "Cox Coverage"),
        caption = "Comparison of SoftBART and Cox Frailty Models for Causal Inference")
}

# Function to create comparison plots
create_comparison_plots <- function(results) {
  # Create a combined coverage plot
  combined_coverage <- (results$spce_coverage_plot + ggtitle("SPCE Coverage")) +
    (results$race_coverage_plot + ggtitle("RACE Coverage")) +
    plot_layout(guides = "collect") & 
    theme(legend.position = "bottom") &
    plot_annotation(
      title = "Coverage Probability Comparison between SoftBART and Cox Frailty Models",
      theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
    )
  
  # Create a combined bias plot
  combined_bias <- (results$spce_bias_plot + ggtitle("SPCE Bias")) +
    (results$race_bias_plot + ggtitle("RACE Bias")) +
    plot_layout(guides = "collect") & 
    theme(legend.position = "bottom") &
    plot_annotation(
      title = "Bias Comparison between SoftBART and Cox Frailty Models",
      theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
    )
  
  # Create a combined estimates plot
  combined_estimates <- (results$spce_est_plot + ggtitle("SPCE Estimates")) +
    (results$race_est_plot + ggtitle("RACE Estimates")) +
    plot_layout(guides = "collect") & 
    theme(legend.position = "bottom") &
    plot_annotation(
      title = "Point Estimates Comparison between SoftBART and Cox Frailty Models",
      theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
    )
  
  # Return all plots
  return(list(
    combined_coverage = combined_coverage,
    combined_bias = combined_bias,
    combined_estimates = combined_estimates,
    ace_plot = results$ace_bar_plot
  ))
}

# Uncomment to run these on your complete simulation results
comp_table <- create_comparison_summary_table(comparison_results)
 comp_plots <- create_comparison_plots(comparison_results)
# 
print(comp_table)
# 
# # Save plots
# ggsave("coverage_comparison.pdf", comp_plots$combined_coverage, width = 10, height = 6)
# ggsave("bias_comparison.pdf", comp_plots$combined_bias, width = 10, height = 6) 
# ggsave("estimates_comparison.pdf", comp_plots$combined_estimates, width = 10, height = 6)
# ggsave("ace_comparison.pdf", comp_plots$ace_plot, width = 8, height = 6)

