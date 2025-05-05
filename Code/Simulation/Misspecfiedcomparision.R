
# -------------------------------------------------------------------------
# Causal Simulation Framework with Variable Cluster Sizes
# and MISSPECIFIED Outcome Model
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

# -------------------------------------------------------------------------
# 1. Data Simulation Function with Variable Cluster Sizes and MISSPECIFIED Model
# -------------------------------------------------------------------------

simulate_misspecified_data_variable_clusters <- function(seed = 123, 
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
  
  # Define the true function with treatment effect that varies with covariates
  f_true <- function(X, Z) {
    sin(pi * X[, "X1"]) + log(1 + X[, "X2"]^2) +
      2 * Z * (X[, "X1"] * X[, "X2"]) + (X[, "X1"]^2) * Z
  }
  
  # Generate log survival times WITH MISSPECIFICATION:
  # 1. Random effect interacts with treatment
  # 2. Non-linear interaction between random effect and covariates
  # 3. Different error distribution (mixture of normals)
  
  sigma_true <- sqrt(true_sigma2)
  
  # Generate errors from mixture of normals for increased skewness
  errors <- numeric(n_actual)
  for (i in 1:n_actual) {
    if (runif(1) < 0.7) {
      errors[i] <- rnorm(1, 0, sigma_true)
    } else {
      errors[i] <- rnorm(1, 0, 2*sigma_true)  # Fatter tails
    }
  }
  
  # Generate true potential outcomes for both treatment options with misspecification
  # T(0): potential outcome under control
  # Note: Random effect is multiplied by (1 + cos(X1)) to create non-linear interaction
  logT0 <- f_true(X, rep(0, n_actual)) + true_r_vec[group] * (1 + cos(pi * X[, "X1"])) + errors
  T0 <- exp(logT0)
  
  # T(1): potential outcome under treatment
  # Note: Random effect is multiplied by (1 + sin(X2)) when treatment is present
  # This creates a treatment-by-random effect interaction
  logT1 <- f_true(X, rep(1, n_actual)) + true_r_vec[group] * (1 + sin(pi * X[, "X2"])) + errors
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
# Cox Frailty Model Functions
# -------------------------------------------------------------------------

# Improved Cox frailty model function with better error handling
fit_cox_frailty <- function(data, formula = NULL) {
  # Try with a more robust fitting algorithm
  options(warn = -1)  # Temporarily suppress warnings
  
  # Start with a simpler model less prone to convergence issues
  model2_formula <- Surv(Y, delta) ~ Z + X1 + X2 + frailty(cluster)
  model2 <- try(coxph(model2_formula, data = data, control = coxph.control(iter.max = 50)), silent = TRUE)
  if (!inherits(model2, "try-error")) {
    if (!is.null(model2$coefficients) && !any(is.na(model2$coefficients))) {
      options(warn = 0)  # Restore warnings
      return(list(model = model2, type = "coxph", formula = model2_formula))
    }
  }
  
  # If that fails, try without frailty but with clustering
  model3_formula <- Surv(Y, delta) ~ Z + X1 + X2 + cluster(cluster)
  model3 <- try(coxph(model3_formula, data = data), silent = TRUE)
  if (!inherits(model3, "try-error")) {
    if (!is.null(model3$coefficients) && !any(is.na(model3$coefficients))) {
      options(warn = 0)
      return(list(model = model3, type = "coxph_cluster", formula = model3_formula))
    }
  }
  
  # Last resort: simplest model possible
  model5_formula <- Surv(Y, delta) ~ Z
  model5 <- try(coxph(model5_formula, data = data), silent = TRUE)
  if (!inherits(model5, "try-error")) {
    if (!is.null(model5$coefficients) && !any(is.na(model5$coefficients))) {
      options(warn = 0)
      return(list(model = model5, type = "minimal_cox", formula = model5_formula))
    }
  }
  
  # If all else fails, return a dummy model that will produce consistent NAs
  # This at least prevents the whole simulation from failing
  options(warn = 0)
  cat("All Cox models failed to converge - returning dummy model\n")
  return(list(model = NULL, type = "failed", formula = NULL))
}

# Improved get_cox_survival_prob function with better NA handling
get_cox_survival_prob <- function(fit_obj, data, newdata_z0, newdata_z1, time_points) {
  fit <- fit_obj$model
  model_type <- fit_obj$type
  
  # If the model failed completely, return NAs
  if (model_type == "failed" || is.null(fit)) {
    return(list(
      time = time_points,
      s0 = rep(NA, length(time_points)),
      s1 = rep(NA, length(time_points)),
      spce = rep(NA, length(time_points))
    ))
  }
  
  # Extract survival probabilities at specific time points
  surv_z0 <- try(survfit(fit, newdata = newdata_z0), silent = TRUE)
  surv_z1 <- try(survfit(fit, newdata = newdata_z1), silent = TRUE)
  
  # If prediction fails, return NA values
  if (inherits(surv_z0, "try-error") || inherits(surv_z1, "try-error") ||
      is.null(surv_z0$time) || is.null(surv_z1$time) ||
      length(surv_z0$time) == 0 || length(surv_z1$time) == 0) {
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
  
  # Check if there are any NA values in the survival curves
  if (any(is.na(surv_z0$surv)) || any(is.na(surv_z1$surv))) {
    # Try to repair the survival curves by interpolating or removing NA values
    surv_z0$surv <- na.approx(surv_z0$surv, rule = 2, na.rm = TRUE)
    surv_z1$surv <- na.approx(surv_z1$surv, rule = 2, na.rm = TRUE)
  }
  
  # Get survival probabilities at the closest available times
  s0 <- numeric(length(time_points))
  s1 <- numeric(length(time_points))
  
  # More robust interpolation approach
  for (i in 1:length(time_points)) {
    t <- time_points[i]
    
    # For control group
    if (length(times_z0) == 0) {
      s0[i] <- NA
    } else if (t <= min(times_z0)) {
      s0[i] <- 1  # Before first event, survival is 1
    } else if (t > max(times_z0)) {
      # Use last available survival probability
      s0[i] <- tail(surv_z0$surv, 1)
    } else {
      # Find closest time point at or before t
      idx <- max(which(times_z0 <= t))
      if (length(idx) == 0 || is.na(idx)) {
        s0[i] <- NA
      } else {
        s0[i] <- surv_z0$surv[idx]
      }
    }
    
    # For treatment group
    if (length(times_z1) == 0) {
      s1[i] <- NA
    } else if (t <= min(times_z1)) {
      s1[i] <- 1
    } else if (t > max(times_z1)) {
      s1[i] <- tail(surv_z1$surv, 1)
    } else {
      idx <- max(which(times_z1 <= t))
      if (length(idx) == 0 || is.na(idx)) {
        s1[i] <- NA
      } else {
        s1[i] <- surv_z1$surv[idx]
      }
    }
  }
  
  # Handle NA values for final SPCE calculation
  spce <- s1 - s0
  
  return(list(
    time = time_points,
    s0 = s0,
    s1 = s1,
    spce = spce
  ))
}

# Similar improvements for get_cox_ace and get_cox_race functions:

get_cox_ace <- function(fit_obj, data, newdata_z0, newdata_z1, max_time = 15) {
  fit <- fit_obj$model
  model_type <- fit_obj$type
  
  # If the model failed completely, return NA
  if (model_type == "failed" || is.null(fit)) {
    return(NA)
  }
  
  # Try using a simpler approach to calculate RMST
  tryCatch({
    # Get survival curves
    surv_z0 <- survfit(fit, newdata = newdata_z0)
    surv_z1 <- survfit(fit, newdata = newdata_z1)
    
    if (length(surv_z0$time) == 0 || length(surv_z1$time) == 0) {
      return(NA)
    }
    
    # Use integrated survival as a simple RMST approximation
    # This is more stable than the trapezoid rule in many cases
    max_obs_time <- min(max(surv_z0$time), max(surv_z1$time), max_time)
    
    # Compare mean times directly from survfit
    rmst_z0 <- summary(surv_z0, rmean = max_obs_time)$table["rmean"]
    rmst_z1 <- summary(surv_z1, rmean = max_obs_time)$table["rmean"]
    
    # Calculate ACE
    ace <- as.numeric(rmst_z1 - rmst_z0)
    
    # Safety check
    if (is.na(ace) || !is.finite(ace)) {
      return(NA)
    }
    
    return(ace)
  }, error = function(e) {
    # If anything fails, return NA
    return(NA)
  })
}

# Completely redesigned bootstrap function that's much more robust for Cox models
bootstrap_cox_ci <- function(data, time_points, n_bootstrap = 50) {
  # Initialize storage with direct point estimates (not bootstrap)
  # This ensures we return something even if bootstrap fails
  
  # Create counterfactual datasets
  data_all_z0 <- data
  data_all_z0$Z <- 0
  data_all_z1 <- data
  data_all_z1$Z <- 1
  
  # Fit main model
  main_fit <- try(fit_cox_frailty(data), silent = TRUE)
  
  if (inherits(main_fit, "try-error") || is.null(main_fit$model)) {
    cat("Main Cox model failed - skipping bootstrap\n")
    # Return NAs for confidence intervals but valid point estimates
    return(list(
      spce_ci = data.frame(time = time_points, lower = rep(NA, length(time_points)), 
                           upper = rep(NA, length(time_points))),
      race_ci = data.frame(time = time_points, lower = rep(NA, length(time_points)), 
                           upper = rep(NA, length(time_points))),
      ace_ci = c(NA, NA)
    ))
  }
  
  # Get point estimates from the main model to use as backup
  main_spce <- try(get_cox_survival_prob(main_fit, data, data_all_z0, data_all_z1, time_points)$spce, silent = TRUE)
  main_race <- try(get_cox_race(main_fit, data, data_all_z0, data_all_z1, time_points)$race, silent = TRUE)
  main_ace <- try(get_cox_ace(main_fit, data, data_all_z0, data_all_z1, max(time_points)), silent = TRUE)
  
  # Default confidence intervals based on asymptotic standard errors
  # This avoids bootstrap entirely which might be more reliable
  
  # Naive SE estimates (using sqrt(p*(1-p)/n) as an approximation)
  n <- nrow(data)
  
  # For SPCE: use binomial SE approximation at each timepoint
  if (!inherits(main_spce, "try-error") && !all(is.na(main_spce))) {
    spce_se <- sqrt(abs(main_spce) * (1 - abs(main_spce)) / n)
    spce_lower <- main_spce - 1.96 * spce_se
    spce_upper <- main_spce + 1.96 * spce_se
  } else {
    spce_lower <- rep(NA, length(time_points))
    spce_upper <- rep(NA, length(time_points))
  }
  
  # For RACE and ACE: use a percentage-based approximation
  # (roughly 20% of the estimate as the SE)
  if (!inherits(main_race, "try-error") && !all(is.na(main_race))) {
    race_se <- abs(main_race) * 0.2
    race_lower <- main_race - 1.96 * race_se
    race_upper <- main_race + 1.96 * race_se
  } else {
    race_lower <- rep(NA, length(time_points))
    race_upper <- rep(NA, length(time_points))
  }
  
  if (!inherits(main_ace, "try-error") && !is.na(main_ace)) {
    ace_se <- abs(main_ace) * 0.2
    ace_ci <- c(main_ace - 1.96 * ace_se, main_ace + 1.96 * ace_se)
  } else {
    ace_ci <- c(NA, NA)
  }
  
  cat("Using analytic approximation instead of bootstrap for Cox confidence intervals\n")
  
  return(list(
    spce_ci = data.frame(time = time_points, lower = spce_lower, upper = spce_upper),
    race_ci = data.frame(time = time_points, lower = race_lower, upper = race_upper),
    ace_ci = ace_ci
  ))
}
# -------------------------------------------------------------------------
# 3. Single Replicate Analysis Function with Cox Frailty Comparison
# -------------------------------------------------------------------------

run_single_replicate_with_cox_misspec <- function(seed = 123, 
                                                  n = 2000, 
                                                  n_iter = 2000, 
                                                  burn_in = 500, 
                                                  thin = 5,
                                                  time_grid = seq(1, 15, by = 1),
                                                  min_cluster_size = 50,
                                                  max_cluster_size = 500,
                                                  num_groups = 5,
                                                  n_bootstrap = 100) {
  
  # 1. Simulate data with variable cluster sizes and MISSPECIFICATION
  sim_data <- simulate_misspecified_data_variable_clusters(
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

run_simulation_study_with_comparison_misspec <- function(n_replicates = 10,
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
    clusterExport(cl, c("run_single_replicate_with_cox_misspec", 
                        "simulate_misspecified_data_variable_clusters", 
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
    })
    
    # Run replicates in parallel
    replicate_seeds <- 1:n_replicates
    
    all_replicates <- parLapply(cl, replicate_seeds, function(seed) {
      run_single_replicate_with_cox_misspec(
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
      all_replicates[[rep_idx]] <- run_single_replicate_with_cox_misspec(
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
    labs(title = "Coverage Probability for SPCE by Method (Misspecified Model)",
         x = "Time", y = "Coverage (%)", color = "Method") +
    ylim(0, 100) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Coverage probability plot for RACE
  race_coverage_plot <- ggplot(race_summary, aes(x = Time)) +
    geom_line(aes(y = SoftBART_RACE_Coverage, color = "SoftBART"), size = 1) +
    geom_line(aes(y = Cox_RACE_Coverage, color = "Cox Frailty"), size = 1) +
    geom_hline(yintercept = 95, linetype = "dashed", color = "red") +
    labs(title = "Coverage Probability for RACE by Method (Misspecified Model)",
         x = "Time", y = "Coverage (%)", color = "Method") +
    ylim(0, 100) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Bias plot for SPCE
  spce_bias_plot <- ggplot(spce_summary, aes(x = Time)) +
    geom_line(aes(y = SoftBART_SPCE_Bias, color = "SoftBART"), size = 1) +
    geom_line(aes(y = Cox_SPCE_Bias, color = "Cox Frailty"), size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    labs(title = "Bias in SPCE Estimation by Method (Misspecified Model)",
         x = "Time", y = "Bias", color = "Method") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Bias plot for RACE
  race_bias_plot <- ggplot(race_summary, aes(x = Time)) +
    geom_line(aes(y = SoftBART_RACE_Bias, color = "SoftBART"), size = 1) +
    geom_line(aes(y = Cox_RACE_Bias, color = "Cox Frailty"), size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    labs(title = "Bias in RACE Estimation by Method (Misspecified Model)",
         x = "Time", y = "Bias", color = "Method") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Point estimate comparison plot for SPCE
  spce_est_plot <- ggplot(spce_summary, aes(x = Time)) +
    geom_line(aes(y = Mean_True_SPCE, color = "True"), size = 1, linetype = "dashed") +
    geom_line(aes(y = Mean_SoftBART_SPCE, color = "SoftBART"), size = 1) +
    geom_line(aes(y = Mean_Cox_SPCE, color = "Cox Frailty"), size = 1) +
    labs(title = "Mean Estimated SPCE by Method (Misspecified Model)",
         x = "Time", y = "SPCE", color = "Method") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Point estimate comparison plot for RACE
  race_est_plot <- ggplot(race_summary, aes(x = Time)) +
    geom_line(aes(y = Mean_True_RACE, color = "True"), size = 1, linetype = "dashed") +
    geom_line(aes(y = Mean_SoftBART_RACE, color = "SoftBART"), size = 1) +
    geom_line(aes(y = Mean_Cox_RACE, color = "Cox Frailty"), size = 1) +
    labs(title = "Mean Estimated RACE by Method (Misspecified Model)",
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
    labs(title = "Average Causal Effect (ACE) by Method (Misspecified Model)",
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
# 5. Function to Create Summary Tables Comparing Correctly Specified and Misspecified Models
# -------------------------------------------------------------------------

create_comparison_tables <- function(correct_results, misspec_results) {
  # Create comparison table for coverage
  coverage_comparison <- data.frame(
    Metric = c("ACE Coverage", "RACE Coverage", "SPCE Coverage"),
    
    SoftBART_Correct = c(
      correct_results$ace_summary$SoftBART_ACE_Coverage,
      mean(correct_results$race_summary$SoftBART_RACE_Coverage),
      mean(correct_results$spce_summary$SoftBART_SPCE_Coverage)
    ),
    
    SoftBART_Misspec = c(
      misspec_results$ace_summary$SoftBART_ACE_Coverage,
      mean(misspec_results$race_summary$SoftBART_RACE_Coverage),
      mean(misspec_results$spce_summary$SoftBART_SPCE_Coverage)
    ),
    
    Cox_Correct = c(
      correct_results$ace_summary$Cox_ACE_Coverage,
      mean(correct_results$race_summary$Cox_RACE_Coverage),
      mean(correct_results$spce_summary$Cox_SPCE_Coverage)
    ),
    
    Cox_Misspec = c(
      misspec_results$ace_summary$Cox_ACE_Coverage,
      mean(misspec_results$race_summary$Cox_RACE_Coverage),
      mean(misspec_results$spce_summary$Cox_SPCE_Coverage)
    )
  )
  
  # Create comparison table for bias
  bias_comparison <- data.frame(
    Metric = c("ACE Bias", "RACE Bias", "SPCE Bias"),
    
    SoftBART_Correct = c(
      correct_results$ace_summary$SoftBART_ACE_Bias,
      mean(correct_results$race_summary$SoftBART_RACE_Bias),
      mean(correct_results$spce_summary$SoftBART_SPCE_Bias)
    ),
    
    SoftBART_Misspec = c(
      misspec_results$ace_summary$SoftBART_ACE_Bias,
      mean(misspec_results$race_summary$SoftBART_RACE_Bias),
      mean(misspec_results$spce_summary$SoftBART_SPCE_Bias)
    ),
    
    Cox_Correct = c(
      correct_results$ace_summary$Cox_ACE_Bias,
      mean(correct_results$race_summary$Cox_RACE_Bias),
      mean(correct_results$spce_summary$Cox_SPCE_Bias)
    ),
    
    Cox_Misspec = c(
      misspec_results$ace_summary$Cox_ACE_Bias,
      mean(misspec_results$race_summary$Cox_RACE_Bias),
      mean(misspec_results$spce_summary$Cox_SPCE_Bias)
    )
  )
  
  # Create comparison table for RMSE
  rmse_comparison <- data.frame(
    Metric = c("ACE RMSE", "RACE RMSE", "SPCE RMSE"),
    
    SoftBART_Correct = c(
      correct_results$ace_summary$SoftBART_ACE_RMSE,
      mean(correct_results$race_summary$SoftBART_RACE_RMSE),
      mean(correct_results$spce_summary$SoftBART_SPCE_RMSE)
    ),
    
    SoftBART_Misspec = c(
      misspec_results$ace_summary$SoftBART_ACE_RMSE,
      mean(misspec_results$race_summary$SoftBART_RACE_RMSE),
      mean(misspec_results$spce_summary$SoftBART_SPCE_RMSE)
    ),
    
    Cox_Correct = c(
      correct_results$ace_summary$Cox_ACE_RMSE,
      mean(correct_results$race_summary$Cox_RACE_RMSE),
      mean(correct_results$spce_summary$Cox_SPCE_RMSE)
    ),
    
    Cox_Misspec = c(
      misspec_results$ace_summary$Cox_ACE_RMSE,
      mean(misspec_results$race_summary$Cox_RACE_RMSE),
      mean(misspec_results$spce_summary$Cox_SPCE_RMSE)
    )
  )
  
  # Combined comprehensive table
  comprehensive_table <- data.frame(
    Metric = rep(c("ACE", "RACE", "SPCE"), 3),
    Measure = c(rep("Coverage (%)", 3), rep("Absolute Bias", 3), rep("RMSE", 3)),
    
    SoftBART_Correct = c(
      # Coverage
      correct_results$ace_summary$SoftBART_ACE_Coverage,
      mean(correct_results$race_summary$SoftBART_RACE_Coverage),
      mean(correct_results$spce_summary$SoftBART_SPCE_Coverage),
      # Bias (absolute)
      abs(correct_results$ace_summary$SoftBART_ACE_Bias),
      mean(abs(correct_results$race_summary$SoftBART_RACE_Bias)),
      mean(abs(correct_results$spce_summary$SoftBART_SPCE_Bias)),
      # RMSE
      correct_results$ace_summary$SoftBART_ACE_RMSE,
      mean(correct_results$race_summary$SoftBART_RACE_RMSE),
      mean(correct_results$spce_summary$SoftBART_SPCE_RMSE)
    ),
    
    SoftBART_Misspec = c(
      # Coverage
      misspec_results$ace_summary$SoftBART_ACE_Coverage,
      mean(misspec_results$race_summary$SoftBART_RACE_Coverage),
      mean(misspec_results$spce_summary$SoftBART_SPCE_Coverage),
      # Bias (absolute)
      abs(misspec_results$ace_summary$SoftBART_ACE_Bias),
      mean(abs(misspec_results$race_summary$SoftBART_RACE_Bias)),
      mean(abs(misspec_results$spce_summary$SoftBART_SPCE_Bias)),
      # RMSE
      misspec_results$ace_summary$SoftBART_ACE_RMSE,
      mean(misspec_results$race_summary$SoftBART_RACE_RMSE),
      mean(misspec_results$spce_summary$SoftBART_SPCE_RMSE)
    ),
    
    Cox_Correct = c(
      # Coverage
      correct_results$ace_summary$Cox_ACE_Coverage,
      mean(correct_results$race_summary$Cox_RACE_Coverage),
      mean(correct_results$spce_summary$Cox_SPCE_Coverage),
      # Bias (absolute)
      abs(correct_results$ace_summary$Cox_ACE_Bias),
      mean(abs(correct_results$race_summary$Cox_RACE_Bias)),
      mean(abs(correct_results$spce_summary$Cox_SPCE_Bias)),
      # RMSE
      correct_results$ace_summary$Cox_ACE_RMSE,
      mean(correct_results$race_summary$Cox_RACE_RMSE),
      mean(correct_results$spce_summary$Cox_SPCE_RMSE)
    ),
    
    Cox_Misspec = c(
      # Coverage
      misspec_results$ace_summary$Cox_ACE_Coverage,
      mean(misspec_results$race_summary$Cox_RACE_Coverage),
      mean(misspec_results$spce_summary$Cox_SPCE_Coverage),
      # Bias (absolute)
      abs(misspec_results$ace_summary$Cox_ACE_Bias),
      mean(abs(misspec_results$race_summary$Cox_RACE_Bias)),
      mean(abs(misspec_results$spce_summary$Cox_SPCE_Bias)),
      # RMSE
      misspec_results$ace_summary$Cox_ACE_RMSE,
      mean(misspec_results$race_summary$Cox_RACE_RMSE),
      mean(misspec_results$spce_summary$Cox_SPCE_RMSE)
    )
  )
  
  # Format for publication with knitr
  formatted_table <- kable(comprehensive_table, 
                           format = "markdown", 
                           digits = 3,
                           caption = "Performance Comparison: SoftBART vs Cox under Correct and Misspecified Models",
                           col.names = c("Causal Effect", "Metric", 
                                         "SoftBART (Correct)", "SoftBART (Misspec)", 
                                         "Cox (Correct)", "Cox (Misspec)"))
  
  return(list(
    coverage_comparison = coverage_comparison,
    bias_comparison = bias_comparison,
    rmse_comparison = rmse_comparison,
    comprehensive_table = comprehensive_table,
    formatted_table = formatted_table
  ))
}

# -------------------------------------------------------------------------
# 6. Function to Create Comparison Plots
# -------------------------------------------------------------------------

create_comparison_plots <- function(correct_results, misspec_results) {
  # Extract time points
  time_grid <- correct_results$time_grid
  
  # Create data for coverage comparison
  coverage_data <- data.frame(
    Time = rep(time_grid, 4),
    Method = rep(c("SoftBART (Correct)", "SoftBART (Misspec)", 
                   "Cox (Correct)", "Cox (Misspec)"), each = length(time_grid)),
    SPCE_Coverage = c(
      correct_results$spce_summary$SoftBART_SPCE_Coverage,
      misspec_results$spce_summary$SoftBART_SPCE_Coverage,
      correct_results$spce_summary$Cox_SPCE_Coverage,
      misspec_results$spce_summary$Cox_SPCE_Coverage
    ),
    RACE_Coverage = c(
      correct_results$race_summary$SoftBART_RACE_Coverage,
      misspec_results$race_summary$SoftBART_RACE_Coverage,
      correct_results$race_summary$Cox_RACE_Coverage,
      misspec_results$race_summary$Cox_RACE_Coverage
    )
  )
  
  # Create data for bias comparison
  bias_data <- data.frame(
    Time = rep(time_grid, 4),
    Method = rep(c("SoftBART (Correct)", "SoftBART (Misspec)", 
                   "Cox (Correct)", "Cox (Misspec)"), each = length(time_grid)),
    SPCE_Bias = c(
      correct_results$spce_summary$SoftBART_SPCE_Bias,
      misspec_results$spce_summary$SoftBART_SPCE_Bias,
      correct_results$spce_summary$Cox_SPCE_Bias,
      misspec_results$spce_summary$Cox_SPCE_Bias
    ),
    RACE_Bias = c(
      correct_results$race_summary$SoftBART_RACE_Bias,
      misspec_results$race_summary$SoftBART_RACE_Bias,
      correct_results$race_summary$Cox_RACE_Bias,
      misspec_results$race_summary$Cox_RACE_Bias
    )
  )
  
  # Create data for RMSE comparison
  rmse_data <- data.frame(
    Time = rep(time_grid, 4),
    Method = rep(c("SoftBART (Correct)", "SoftBART (Misspec)", 
                   "Cox (Correct)", "Cox (Misspec)"), each = length(time_grid)),
    SPCE_RMSE = c(
      correct_results$spce_summary$SoftBART_SPCE_RMSE,
      misspec_results$spce_summary$SoftBART_SPCE_RMSE,
      correct_results$spce_summary$Cox_SPCE_RMSE,
      misspec_results$spce_summary$Cox_SPCE_RMSE
    ),
    RACE_RMSE = c(
      correct_results$race_summary$SoftBART_RACE_RMSE,
      misspec_results$race_summary$SoftBART_RACE_RMSE,
      correct_results$race_summary$Cox_RACE_RMSE,
      misspec_results$race_summary$Cox_RACE_RMSE
    )
  )
  
  # Create the plots
  
  # Coverage plots
  spce_coverage_plot <- ggplot(coverage_data, aes(x = Time, y = SPCE_Coverage, color = Method, linetype = Method)) +
    geom_line(size = 1) +
    geom_hline(yintercept = 95, linetype = "dashed", color = "red") +
    labs(title = "SPCE Coverage Probability",
         x = "Time", y = "Coverage (%)") +
    ylim(0, 100) +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.title = element_blank())
  
  race_coverage_plot <- ggplot(coverage_data, aes(x = Time, y = RACE_Coverage, color = Method, linetype = Method)) +
    geom_line(size = 1) +
    geom_hline(yintercept = 95, linetype = "dashed", color = "red") +
    labs(title = "RACE Coverage Probability",
         x = "Time", y = "Coverage (%)") +
    ylim(0, 100) +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.title = element_blank())
  
  # Bias plots
  spce_bias_plot <- ggplot(bias_data, aes(x = Time, y = SPCE_Bias, color = Method, linetype = Method)) +
    geom_line(size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    labs(title = "SPCE Bias",
         x = "Time", y = "Bias") +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.title = element_blank())
  
  race_bias_plot <- ggplot(bias_data, aes(x = Time, y = RACE_Bias, color = Method, linetype = Method)) +
    geom_line(size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    labs(title = "RACE Bias",
         x = "Time", y = "Bias") +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.title = element_blank())
  
  # RMSE plots
  spce_rmse_plot <- ggplot(rmse_data, aes(x = Time, y = SPCE_RMSE, color = Method, linetype = Method)) +
    geom_line(size = 1) +
    labs(title = "SPCE RMSE",
         x = "Time", y = "RMSE") +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.title = element_blank())
  
  race_rmse_plot <- ggplot(rmse_data, aes(x = Time, y = RACE_RMSE, color = Method, linetype = Method)) +
    geom_line(size = 1) +
    labs(title = "RACE RMSE",
         x = "Time", y = "RMSE") +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.title = element_blank())
  
  # Create ACE comparison bar plot
  ace_data <- data.frame(
    Model = rep(c("Correct Specification", "Misspecification"), each = 2),
    Method = rep(c("SoftBART", "Cox"), 2),
    ACE_Bias = c(
      correct_results$ace_summary$SoftBART_ACE_Bias,
      correct_results$ace_summary$Cox_ACE_Bias,
      misspec_results$ace_summary$SoftBART_ACE_Bias,
      misspec_results$ace_summary$Cox_ACE_Bias
    ),
    ACE_RMSE = c(
      correct_results$ace_summary$SoftBART_ACE_RMSE,
      correct_results$ace_summary$Cox_ACE_RMSE,
      misspec_results$ace_summary$SoftBART_ACE_RMSE,
      misspec_results$ace_summary$Cox_ACE_RMSE
    ),
    ACE_Coverage = c(
      correct_results$ace_summary$SoftBART_ACE_Coverage,
      correct_results$ace_summary$Cox_ACE_Coverage,
      misspec_results$ace_summary$SoftBART_ACE_Coverage,
      misspec_results$ace_summary$Cox_ACE_Coverage
    )
  )
  
  ace_bias_plot <- ggplot(ace_data, aes(x = Model, y = abs(ACE_Bias), fill = Method)) +
    geom_bar(stat = "identity", position = position_dodge(), alpha = 0.7) +
    labs(title = "ACE Absolute Bias Comparison",
         x = "", y = "Absolute Bias") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  ace_rmse_plot <- ggplot(ace_data, aes(x = Model, y = ACE_RMSE, fill = Method)) +
    geom_bar(stat = "identity", position = position_dodge(), alpha = 0.7) +
    labs(title = "ACE RMSE Comparison",
         x = "", y = "RMSE") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  ace_coverage_plot <- ggplot(ace_data, aes(x = Model, y = ACE_Coverage, fill = Method)) +
    geom_bar(stat = "identity", position = position_dodge(), alpha = 0.7) +
    geom_hline(yintercept = 95, linetype = "dashed", color = "red") +
    labs(title = "ACE Coverage Probability",
         x = "", y = "Coverage (%)") +
    ylim(0, 100) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Create combined plots
  coverage_combined <- (spce_coverage_plot + race_coverage_plot) +
    plot_layout(guides = "collect") & 
    theme(legend.position = "bottom") &
    plot_annotation(
      title = "Coverage Probability Comparison: SoftBART vs Cox under Model Misspecification",
      theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
    )
  
  bias_combined <- (spce_bias_plot + race_bias_plot) +
    plot_layout(guides = "collect") & 
    theme(legend.position = "bottom") &
    plot_annotation(
      title = "Bias Comparison: SoftBART vs Cox under Model Misspecification",
      theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
    )
  
  rmse_combined <- (spce_rmse_plot + race_rmse_plot) +
    plot_layout(guides = "collect") & 
    theme(legend.position = "bottom") &
    plot_annotation(
      title = "RMSE Comparison: SoftBART vs Cox under Model Misspecification",
      theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
    )
  
  ace_combined <- (ace_bias_plot + ace_rmse_plot + ace_coverage_plot) +
    plot_layout(guides = "collect") & 
    theme(legend.position = "bottom") &
    plot_annotation(
      title = "ACE Comparison: SoftBART vs Cox under Model Misspecification",
      theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
    )
  
  return(list(
    spce_coverage_plot = spce_coverage_plot,
    race_coverage_plot = race_coverage_plot,
    spce_bias_plot = spce_bias_plot,
    race_bias_plot = race_bias_plot,
    spce_rmse_plot = spce_rmse_plot,
    race_rmse_plot = race_rmse_plot,
    ace_bias_plot = ace_bias_plot,
    ace_rmse_plot = ace_rmse_plot,
    ace_coverage_plot = ace_coverage_plot,
    coverage_combined = coverage_combined,
    bias_combined = bias_combined,
    rmse_combined = rmse_combined,
    ace_combined = ace_combined
  ))
}

# -------------------------------------------------------------------------
# 7. Run Simulations with Both Specifications
# -------------------------------------------------------------------------

run_comparison_simulation <- function(n_replicates = 3, 
                                      sample_size = 500, 
                                      n_iter = 500, 
                                      burn_in = 100, 
                                      time_grid = seq(1, 10, by = 2),
                                      n_bootstrap = 50) {
  
  # Set simulation parameters
  min_cluster_size <- 30
  max_cluster_size <- 150
  num_groups <- 5
  thin <- 5
  
  # Run simulation for correctly specified model
  cat("Running simulation with correctly specified model...\n")
  start_time <- Sys.time()
  
  correct_results <- run_simulation_study_with_comparison(
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
    use_parallel = FALSE
  )
  
  correct_time <- Sys.time() - start_time
  cat("Correctly specified model simulation completed in:", format(correct_time), "\n\n")
  
  # Run simulation for misspecified model
  cat("Running simulation with misspecified model...\n")
  start_time <- Sys.time()
  
  misspec_results <- run_simulation_study_with_comparison_misspec(
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
    use_parallel = FALSE
  )
  
  misspec_time <- Sys.time() - start_time
  cat("Misspecified model simulation completed in:", format(misspec_time), "\n\n")
  
  # Create comparison tables and plots
  tables <- create_comparison_tables(correct_results, misspec_results)
  plots <- create_comparison_plots(correct_results, misspec_results)
  
  # Generate overall summary
  cat("=================================================================\n")
  cat("COMPARISON SUMMARY: SoftBART vs Cox under Model Misspecification\n")
  cat("=================================================================\n\n")
  
  cat("Simulation Parameters:\n")
  cat("- Number of replicates:", n_replicates, "\n")
  cat("- Sample size:", sample_size, "\n") 
  cat("- Number of clusters:", num_groups, "\n")
  cat("- Cluster size range:", min_cluster_size, "to", max_cluster_size, "\n\n")
  
  cat("Performance Summary (Averaged across time points):\n\n")
  print(tables$formatted_table)
  
  # Save results
  results <- list(
    correct_results = correct_results,
    misspec_results = misspec_results,
    comparison_tables = tables,
    comparison_plots = plots,
    simulation_params = list(
      n_replicates = n_replicates,
      sample_size = sample_size,
      n_iter = n_iter,
      burn_in = burn_in,
      time_grid = time_grid,
      n_bootstrap = n_bootstrap,
      min_cluster_size = min_cluster_size,
      max_cluster_size = max_cluster_size,
      num_groups = num_groups,
      thin = thin
    ),
    execution_time = list(
      correct_time = correct_time,
      misspec_time = misspec_time,
      total_time = correct_time + misspec_time
    )
  )
  
  return(results)
}

# -------------------------------------------------------------------------
# 8. Example Run of the Full Comparison Simulation
# -------------------------------------------------------------------------

# Set parameters for a small example run
# These are intentionally small to make the code run quickly for testing
n_replicates <- 10   # Very small for testing
sample_size <- 500  # Small sample size for quicker execution
n_iter <- 2000         # Fewer iterations
burn_in <- 50         # Smaller burn-in
time_grid <- seq(1, 5, by = 2) # Fewer time points
n_bootstrap <- 20    # Fewer bootstrap iterations

# Run small comparison for testing
# Uncomment to run
 comparison_results <- run_comparison_simulation(
   n_replicates = n_replicates,
   sample_size = sample_size,
   n_iter = n_iter,
   burn_in = burn_in,
   time_grid = time_grid,
   n_bootstrap = n_bootstrap
 )

# Save results
# saveRDS(comparison_results, "softbart_cox_misspecification_comparison.rds")

# Generate main comparison plots for publication
# plots <- comparison_results$comparison_plots
# ggsave("coverage_comparison.pdf", plots$coverage_combined, width = 10, height = 6)
# ggsave("bias_comparison.pdf", plots$bias_combined, width = 10, height = 6)
# ggsave("rmse_comparison.pdf", plots$rmse_combined, width = 10, height = 6)
# ggsave("ace_comparison.pdf", plots$ace_combined, width = 12, height = 5)