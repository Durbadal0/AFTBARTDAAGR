# Install and load necessary packages
if (!requireNamespace("mvtnorm", quietly = TRUE)) {
  install.packages("mvtnorm")
}
if (!requireNamespace("SoftBart", quietly = TRUE)) {
  install.packages("SoftBart")
}
library(mvtnorm)
library(SoftBart)

# Function: Bayesian mixed model with spatial CAR random effects and nonparametric f() via SoftBart,
# now treating rho as unknown and updating it via a Metropolis-Hastings step.
AFT_mixed_CAR_softbart <- function(time,status, X, group, W, X_test, group_test = group,
                                        n_iter = 6000, burn_in = 1000, thin = 5,
                                        a_sigma = 1, b_sigma = 1,
                                        a_r = 1, b_r = 1,
                                        rho_r = 0.8,        # initial value for rho
                                        rho_prop_sd = 0.05,   # proposal sd for rho update
                                        hypers = NULL, opts = NULL) {
  # y: response vector (length n)
  # X: training covariate matrix for f() (n x p)
  # group: vector (or factor) indicating spatial group (area) membership for each observation
  # W: I x I spatial adjacency matrix (assumed symmetric, zeros on the diagonal)
  # X_test: test set covariate matrix for f() (n_test x p)
  # group_test: group indicator for X_test (default: same as training groups)
  # n_iter, burn_in, thin: MCMC settings
  # a_sigma, b_sigma: hyperparameters for the inverse-gamma prior on sigma^2
  # a_r, b_r: hyperparameters for the inverse-gamma prior on sigma_r^2 (CAR variance)
  # rho_r: initial value for the spatial autocorrelation parameter
  # rho_prop_sd: proposal standard deviation for the Metropolis-Hastings update of rho
  # hypers, opts: optional SoftBart hyperparameters and options. If NULL, defaults are used.
  y<-log(time)
  n <- length(y)
  p <- ncol(X)
  
  # Unique groups and indices:
  groups <- sort(unique(group))
  I <- length(groups)
  group_indices <- lapply(groups, function(g) which(group == g))
  n_i <- sapply(group_indices, length)
  
  # Precompute D (diagonal matrix of row sums of W)
  D_mat <- diag(rowSums(W))
  
  # Set up SoftBart forest for f() update:
  if (is.null(hypers)) {
    hypers <- Hypers(X = X, Y = y)
  }
  if (is.null(opts)) {
    opts <- Opts(num_burn = 0, num_save = 1, num_thin = 1)
  }
  forest <- MakeForest(hypers, opts)
  
  # Initialize f() using the current forest predictions:
  f_current <- forest$do_predict(X)
  f_test_current <- forest$do_predict(X_test)
  
  # Initialize spatial random effects (one per group) and variance parameters:
  r_current <- rep(0, I)
  sigma2_current <- 1
  sigma_r2_current <- 1
  rho_current <- rho_r  # starting value for rho
  
  # Storage for posterior samples:
  n_save <- floor((n_iter - burn_in) / thin)
  f_samples <- matrix(0, nrow = n_save, ncol = n)         # f(x) for training
  f_test_samples <- matrix(0, nrow = n_save, ncol = nrow(X_test))  # f(x) for test set
  r_samples <- matrix(0, nrow = n_save, ncol = I)           # spatial random effects
  sigma2_samples <- numeric(n_save)
  sigma_r2_samples <- numeric(n_save)
  rho_samples <- numeric(n_save)
  
  sample_index <- 1
  
  for (iter in 1:n_iter) {
    
    
    
    
    
    ## 1. Update f(·) via SoftBart:
    # For each observation, subtract the corresponding spatial effect:
    r_obs <- sapply(group, function(g) r_current[which(groups == g)])
    resid_for_f <- y - r_obs  # f(x) + e
    # Update forest by one Gibbs iteration:
    forest$do_gibbs(X, resid_for_f, X_test, 1)
    # Get updated predictions:
    f_current <- forest$do_predict(X)
    f_test_current <- forest$do_predict(X_test)
    
    ## 2. Update spatial CAR random effects r (as a block):
    # Compute b_i = sum_{j in group i} (y_j - f_current[j]) / sigma2_current.
    b_vec <- sapply(group_indices, function(idx) sum(y[idx] - f_current[idx])) / sigma2_current
    # Likelihood precision from data: diag(n_i)/sigma2_current.
    Lik_prec <- diag(n_i) / sigma2_current
    # CAR prior precision: Q = D - rho*W, where rho = rho_current.
    Q_current_mat <- D_mat - rho_current * W
    # Full conditional precision and covariance:
    Prec_r <- Lik_prec + Q_current_mat / sigma_r2_current
    V_r <- solve(Prec_r)
    mu_r <- V_r %*% b_vec
    r_current <- as.vector(rmvnorm(1, mean = as.vector(mu_r), sigma = V_r))
    
    ## 3. Update sigma^2 (residual variance):
    r_obs <- sapply(group, function(g) r_current[which(groups == g)])
    resid <- y - f_current - r_obs
    shape_sigma <- a_sigma + n/2
    rate_sigma <- b_sigma + 0.5 * sum(resid^2)
    sigma2_current <- 1 / rgamma(1, shape = shape_sigma, rate = rate_sigma)
    
    ## 4. Update sigma_r^2 (CAR variance):
    shape_r <- a_r + I/2
    # Note: use Q based on current rho.
    rate_r <- b_r + 0.5 * as.numeric(t(r_current) %*% (D_mat - rho_current * W) %*% r_current)
    sigma_r2_current <- 1 / rgamma(1, shape = shape_r, rate = rate_r)
    
    ## 5. Metropolis-Hastings update for rho:
    # Propose a new rho via a random-walk proposal:
    rho_prop <- rho_current + rnorm(1, mean = 0, sd = rho_prop_sd)
    # Check that the proposal is within (0,1). (You can adjust the boundaries as needed.)
    if (rho_prop > 0 && rho_prop < 1) {
      # Compute Q for the proposed and current rho:
      Q_prop <- D_mat - rho_prop * W
      log_det_prop <- as.numeric(determinant(Q_prop, logarithm = TRUE)$modulus)
      log_det_current <- as.numeric(determinant(D_mat - rho_current * W, logarithm = TRUE)$modulus)
      # Compute the log density for r under the CAR prior:
      # r ~ MVN(0, sigma_r2 * Q^{-1}) so density ∝ |Q|^(1/2) exp{-0.5/sigma_r2 * r' Q r}
      log_target_prop <- 0.5 * log_det_prop - 0.5/sigma_r2_current * as.numeric(t(r_current) %*% Q_prop %*% r_current)
      log_target_current <- 0.5 * log_det_current - 0.5/sigma_r2_current * as.numeric(t(r_current) %*% (D_mat - rho_current * W) %*% r_current)
      
      # Since the proposal is symmetric and we assume a uniform prior on rho, acceptance probability is:
      log_accept_ratio <- log_target_prop - log_target_current
      if (log(runif(1)) < log_accept_ratio) {
        rho_current <- rho_prop
      }
    }
    # (If the proposal is out of bounds, we simply keep the current rho.)
    
    ## Save samples after burn-in and thinning:
    if (iter > burn_in && ((iter - burn_in) %% thin == 0)) {
      f_samples[sample_index, ] <- f_current
      f_test_samples[sample_index, ] <- f_test_current
      r_samples[sample_index, ] <- r_current
      sigma2_samples[sample_index] <- sigma2_current
      sigma_r2_samples[sample_index] <- sigma_r2_current
      rho_samples[sample_index] <- rho_current
      sample_index <- sample_index + 1
    }
  }
  
  return(list(f_train = f_samples,
              f_test = f_test_samples,
              r = r_samples,
              sigma2 = sigma2_samples,
              sigma_r2 = sigma_r2_samples,
              rho = rho_samples,
              forest = forest))
}





# --- Define the New Function for Bayesian Causal Inference ---
AFT_mixed_CAR_softbart_causal <- function(time, status, X, Z, group, W, 
                                          X_test, Z_test = Z, group_test = group,
                                          n_iter = 6000, burn_in = 1000, thin = 5,
                                          a_sigma = 1, b_sigma = 1,
                                          a_r = 1, b_r = 1,
                                          rho_r = 0.8,        # initial value for rho
                                          rho_prop_sd = 0.05,   # proposal sd for rho update
                                          hypers = NULL, opts = NULL,
                                          n_iter_ps=2000) {
  # time, status: survival outcome and censoring indicator (here we assume observed times)
  # X: training covariate matrix (n x p)
  # Z: binary treatment indicator for training data
  # group: spatial group indicator for each observation
  # W: spatial adjacency matrix for groups (symmetric, zeros on diagonal)
  # X_test: test set covariate matrix
  # Z_test: test set treatment indicator (default same as training if not provided)
  # group_test: test set group indicator (default same as training)
  #
  # The model is:
  #   log T_{ij} = f(X_{ij}, Z_{ij}, e_hat(X_{ij})) + r_i + epsilon,   epsilon ~ N(0,sigma^2)
  # where e_hat(.) is the estimated propensity score (using SBART) for P(Z=1|X).
  
  n <- length(time)
  p <- ncol(X)
  
  # --- 0. Estimate the Propensity Score Using SBART ---
  # We use SoftBart to fit a regression on the binary outcome Z.
  # Since SoftBart is a regression routine, we obtain raw predictions and then transform them via a logistic (expit) function.
  expit <- function(x) 1/(1+exp(-x))
  
  hypers_ps <- Hypers(X = X, Y = Z)  # Note: even though Z is binary, we use SBART to obtain a latent score.
  opts_ps <- Opts(num_burn = 0, num_save = 1, num_thin = 1)
  forest_ps <- MakeForest(hypers_ps, opts_ps)
  forest_ps$do_gibbs(X,Z,X,n_iter_ps)
  # Obtain raw predictions on training and test sets
  ps_raw_train <- forest_ps$do_predict(X)
  ps_raw_test  <- forest_ps$do_predict(X_test)
  
  # Transform to [0,1] probability scale via the logistic function
  e_hat_train <- expit(ps_raw_train)
  e_hat_test  <- expit(ps_raw_test)
  
  # --- 1. Create Augmented Covariate Matrices ---
  # Augment X with treatment indicator Z and the estimated propensity score
  X_aug      <- cbind(X, Z, e_hat_train)
  X_test_aug <- cbind(X_test, Z_test, e_hat_test)
  
  # --- 2. Prepare Data for the AFT Model ---
  # Transform survival times to the log scale
  y <- log(time)
  
  # Identify unique groups and indices:
  groups <- sort(unique(group))
  I <- length(groups)
  group_indices <- lapply(groups, function(g) which(group == g))
  n_i <- sapply(group_indices, length)
  
  # Precompute D (diagonal matrix of row sums of W)
  D_mat <- diag(rowSums(W))
  
  # --- 3. Set Up SoftBart for f(.) ---
  if (is.null(hypers)) {
    hypers <- Hypers(X = X_aug, Y = y)
  }
  if (is.null(opts)) {
    opts <- Opts(num_burn = 0, num_save = 1, num_thin = 1)
  }
  forest <- MakeForest(hypers, opts)
  
  # Initialize f() using the current forest predictions:
  f_current      <- forest$do_predict(X_aug)
  f_test_current <- forest$do_predict(X_test_aug)
  
  # Initialize spatial random effects (one per group) and variance parameters:
  r_current <- rep(0, I)
  sigma2_current <- 1
  sigma_r2_current <- 1
  rho_current <- rho_r  # starting value for rho
  
  # Storage for posterior samples:
  n_save <- floor((n_iter - burn_in) / thin)
  f_samples <- matrix(0, nrow = n_save, ncol = n)               # f(x) for training
  f_test_samples <- matrix(0, nrow = n_save, ncol = nrow(X_test_aug))  # f(x) for test set
  r_samples <- matrix(0, nrow = n_save, ncol = I)                 # spatial random effects
  sigma2_samples <- numeric(n_save)
  sigma_r2_samples <- numeric(n_save)
  rho_samples <- numeric(n_save)
  
  sample_index <- 1
  
  # --- 4. MCMC Sampling ---
  for (iter in 1:n_iter) {
    
    ## 4.1 Update f(·) via SoftBart:
    # For each observation, subtract the corresponding spatial effect:
    r_obs <- sapply(group, function(g) r_current[which(groups == g)])
    resid_for_f <- y - r_obs  # f(x) + error
    
    # Update the forest by one Gibbs iteration:
    forest$do_gibbs(X_aug, resid_for_f, X_test_aug, 1)
    
    # Update predictions:
    f_current      <- forest$do_predict(X_aug)
    f_test_current <- forest$do_predict(X_test_aug)
    
    ## 4.2 Update Spatial CAR Random Effects r (as a Block):
    # Compute b_i = sum_{j in group i} (y_j - f_current[j]) / sigma2_current.
    b_vec <- sapply(group_indices, function(idx) sum(y[idx] - f_current[idx])) / sigma2_current
    # Likelihood precision from data: diag(n_i)/sigma2_current.
    Lik_prec <- diag(n_i) / sigma2_current
    # CAR prior precision: Q = D - rho * W.
    Q_current_mat <- D_mat - rho_current * W
    # Full conditional precision and covariance:
    Prec_r <- Lik_prec + Q_current_mat / sigma_r2_current
    V_r <- solve(Prec_r)
    mu_r <- V_r %*% b_vec
    r_current <- as.vector(rmvnorm(1, mean = as.vector(mu_r), sigma = V_r))
    
    ## 4.3 Update Residual Variance sigma^2:
    r_obs <- sapply(group, function(g) r_current[which(groups == g)])
    resid <- y - f_current - r_obs
    shape_sigma <- a_sigma + n/2
    rate_sigma <- b_sigma + 0.5 * sum(resid^2)
    sigma2_current <- 1 / rgamma(1, shape = shape_sigma, rate = rate_sigma)
    
    ## 4.4 Update CAR Variance sigma_r^2:
    shape_r <- a_r + I/2
    rate_r <- b_r + 0.5 * as.numeric(t(r_current) %*% (D_mat - rho_current * W) %*% r_current)
    sigma_r2_current <- 1 / rgamma(1, shape = shape_r, rate = rate_r)
    
    ## 4.5 Metropolis-Hastings Update for rho:
    rho_prop <- rho_current + rnorm(1, mean = 0, sd = rho_prop_sd)
    if (rho_prop > 0 && rho_prop < 1) {
      Q_prop <- D_mat - rho_prop * W
      log_det_prop <- as.numeric(determinant(Q_prop, logarithm = TRUE)$modulus)
      log_det_current <- as.numeric(determinant(Q_current_mat, logarithm = TRUE)$modulus)
      log_target_prop <- 0.5 * log_det_prop - 0.5/sigma_r2_current * as.numeric(t(r_current) %*% Q_prop %*% r_current)
      log_target_current <- 0.5 * log_det_current - 0.5/sigma_r2_current * as.numeric(t(r_current) %*% Q_current_mat %*% r_current)
      log_accept_ratio <- log_target_prop - log_target_current
      if (log(runif(1)) < log_accept_ratio) {
        rho_current <- rho_prop
      }
    }
    
    ## 4.6 Save Samples After Burn-in and Thinning:
    if (iter > burn_in && ((iter - burn_in) %% thin == 0)) {
      f_samples[sample_index, ] <- f_current
      f_test_samples[sample_index, ] <- f_test_current
      r_samples[sample_index, ] <- r_current
      sigma2_samples[sample_index] <- sigma2_current
      sigma_r2_samples[sample_index] <- sigma_r2_current
      rho_samples[sample_index] <- rho_current
      sample_index <- sample_index + 1
    }
  }
  
  return(list(f_train = f_samples,
              f_test = f_test_samples,
              r = r_samples,
              sigma2 = sigma2_samples,
              sigma_r2 = sigma_r2_samples,
              rho = rho_samples,
              forest = forest,
              propensity_forest = forest_ps))
}








# Helper: logistic (expit) function
expit <- function(x) {
  1 / (1 + exp(-x))
}

# Function to obtain CSPCE: 
# ΔCSPCE(t, X, r) = S(1)(t|X, r) - S(0)(t|X, r)
getCSPCE <- function(t, X, r, fit) {
  # Ensure X is a one-row matrix
  X <- as.matrix(X)
  
  # Obtain estimated propensity score for new X using the fitted propensity forest:
  ps_raw <- fit$propensity_forest$do_predict(X)
  e_hat <- expit(ps_raw)
  
  # Create augmented covariate matrices for treatment = 1 and 0:
  X_aug1 <- cbind(X, 1, e_hat)
  X_aug0 <- cbind(X, 0, e_hat)
  
  # Get MCMC draws for f(·) under treatment 1 and 0 from the outcome forest:
  # (Assuming do_predict returns a vector of predictions, one for each MCMC draw.)
  f1_samples <- as.vector(fit$forest$do_predict(X_aug1))
  f0_samples <- as.vector(fit$forest$do_predict(X_aug0))
  
  # Extract the posterior draws for the residual variance
  sigma2_samples <- fit$sigma2
  
  # Compute survival probabilities for each draw:
  logt <- log(t)
  S1 <- 1 - pnorm((logt - (f1_samples + r)) / sqrt(sigma2_samples))
  S0 <- 1 - pnorm((logt - (f0_samples + r)) / sqrt(sigma2_samples))
  
  # Compute the posterior draws of ΔCSPCE and summarize:
  delta_CSPCE <- S1 - S0
  list(estimate = mean(delta_CSPCE),
       CI = quantile(delta_CSPCE, probs = c(0.025, 0.975)))
}

# Function to obtain CACE: 
# ΔCACE(X, r) = E[T(1)|X, r] - E[T(0)|X, r]
getCACE <- function(X, r, fit) {
  # Ensure X is a one-row matrix
  X <- as.matrix(X)
  
  # Obtain estimated propensity score for new X:
  ps_raw <- fit$propensity_forest$do_predict(X)
  e_hat <- expit(ps_raw)
  
  # Create augmented covariate matrices for treatment = 1 and 0:
  X_aug1 <- cbind(X, 1, e_hat)
  X_aug0 <- cbind(X, 0, e_hat)
  
  # Get MCMC draws for f(·) under treatment 1 and 0:
  f1_samples <- as.vector(fit$forest$do_predict(X_aug1))
  f0_samples <- as.vector(fit$forest$do_predict(X_aug0))
  
  # Extract σ² draws:
  sigma2_samples <- fit$sigma2
  
  # For a lognormal distribution: E[T] = exp(mu + σ²/2),
  # where mu = f(·) + r.
  ET1 <- exp(f1_samples + r + sigma2_samples/2)
  ET0 <- exp(f0_samples + r + sigma2_samples/2)
  
  delta_CACE <- ET1 - ET0
  list(estimate = mean(delta_CACE),
       CI = quantile(delta_CACE, probs = c(0.025, 0.975)))
}

# Function to obtain CRACE: 
# ΔCRACE(X, r) = E[min(T(1), t*)|X, r] - E[min(T(0), t*)|X, r]
getCRACE <- function(t_star, X, r, fit) {
  # Ensure X is a one-row matrix
  X <- as.matrix(X)
  
  # Obtain estimated propensity score for new X:
  ps_raw <- fit$propensity_forest$do_predict(X)
  e_hat <- expit(ps_raw)
  
  # Create augmented covariate matrices for treatment = 1 and 0:
  X_aug1 <- cbind(X, 1, e_hat)
  X_aug0 <- cbind(X, 0, e_hat)
  
  # Get MCMC draws for f(·) under treatment 1 and 0:
  f1_samples <- as.vector(fit$forest$do_predict(X_aug1))
  f0_samples <- as.vector(fit$forest$do_predict(X_aug0))
  
  # Extract σ² draws:
  sigma2_samples <- fit$sigma2
  
  # For each MCMC draw, let mu = f(·) + r and σ = sqrt(σ²).
  # Then, one can show that:
  # E[min(T, t_star)] = exp(mu + σ²/2)*Φ((log(t_star) - mu - σ²)/σ) + 
  #                      t_star * [1 - Φ((log(t_star) - mu)/σ)]
  mu1 <- f1_samples + r
  mu0 <- f0_samples + r
  sigma1 <- sqrt(sigma2_samples)
  sigma0 <- sqrt(sigma2_samples)
  
  EM1 <- exp(mu1 + sigma2_samples/2) * pnorm((log(t_star) - mu1 - sigma2_samples) / sigma1) +
    t_star * (1 - pnorm((log(t_star) - mu1) / sigma1))
  EM0 <- exp(mu0 + sigma2_samples/2) * pnorm((log(t_star) - mu0 - sigma2_samples) / sigma0) +
    t_star * (1 - pnorm((log(t_star) - mu0) / sigma0))
  
  delta_CRACE <- EM1 - EM0
  list(estimate = mean(delta_CRACE),
       CI = quantile(delta_CRACE, probs = c(0.025, 0.975)))
}



#' Compare True vs Estimated Functions in SoftBART Models
#'
#' This function compares the true function values with the estimated function from 
#' a SoftBART model for causal inference in survival data. It can produce comparative 
#' plots and calculate metrics for assessment.
#'
#' @param fit The fitted AFT_mixed_CAR_softbart_causal model object
#' @param X_test Matrix of test covariates to evaluate the functions on
#' @param Z_test Vector of treatment indicators for test data
#' @param true_f Function that takes X and Z and returns the true function values
#' @param group_test Vector of group/cluster indicators for test data
#' @param true_r_vec Vector of true spatial random effects for each cluster
#' @param plot Logical; whether to produce plots (default: TRUE)
#' @param metrics Logical; whether to calculate comparison metrics (default: TRUE)
#' @param grid_size Number of points to use if creating a grid for visualization
#' @return A list containing comparison metrics and optionally plots
#'
#' @examples
#' # Assuming fit is your fitted model and other parameters are defined:
#' results <- compare_softbart_functions(
#'   fit = fit, 
#'   X_test = X, 
#'   Z_test = Z, 
#'   true_f = function(X, Z) { 0.5 * X[, "X1"] - 0.3 * X[, "X2"] + 1.0 * Z },
#'   group_test = group,
#'   true_r_vec = true_r_vec
#' )
#' Compare True vs Estimated Functions in SoftBART Models (Without Random Effects)
#'
#' This function compares the true function values with the estimated function from 
#' a SoftBART model for causal inference, excluding the random effects.
#'
#' @param fit The fitted AFT_mixed_CAR_softbart_causal model object
#' @param X_test Matrix of test covariates to evaluate the functions on
#' @param Z_test Vector of treatment indicators for test data
#' @param true_f Function that takes X and Z and returns the true function values
#' @param plot Logical; whether to produce plots (default: TRUE)
#' @param metrics Logical; whether to calculate comparison metrics (default: TRUE)
#' @param grid_size Number of points to use if creating a grid for visualization
#' @return A list containing comparison metrics and optionally plots
#'
#' @examples
#' # Assuming fit is your fitted model and other parameters are defined:
#' results <- compare_softbart_functions(
#'   fit = fit, 
#'   X_test = X, 
#'   Z_test = Z, 
#'   true_f = function(X, Z) { 0.5 * X[, "X1"] - 0.3 * X[, "X2"] + 1.0 * Z }
#' )
# First, use the original compare_softbart_functions
compare_softbart_functions <- function(fit, X_test, Z_test, true_f, 
                                       plot = TRUE, metrics = TRUE,
                                       grid_size = 100) {
  
  # Ensure X_test is a matrix
  X_test <- as.matrix(X_test)
  
  # Get estimated propensity scores for test data
  ps_raw <- fit$propensity_forest$do_predict(X_test)
  e_hat <- 1 / (1 + exp(-ps_raw))
  
  # Create augmented test data for predictions
  X_aug <- cbind(X_test, Z_test, e_hat)
  
  # Get posterior mean of the function estimates (without random effects)
  f_est <- colMeans(fit$f_test)
  
  # Calculate true function values (without random effects)
  f_true <- true_f(X_test, Z_test)
  
  results <- list()
  
  # Calculate metrics if requested
  if (metrics) {
    # Root Mean Squared Error
    rmse <- sqrt(mean((f_est - f_true)^2))
    # Mean Absolute Error
    mae <- mean(abs(f_est - f_true))
    # Correlation
    corr <- cor(f_est, f_true)
    
    results$metrics <- list(
      rmse = rmse,
      mae = mae,
      correlation = corr
    )
  }
  
  # Create plots if requested
  if (plot) {
    # Prepare data for plotting
    plot_data <- data.frame(
      True = f_true,
      Estimated = f_est,
      Z = Z_test
    )
    
    # Scatter plot of true vs estimated
    scatter_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = True, y = Estimated, color = factor(Z))) +
      ggplot2::geom_point(alpha = 0.7) +
      ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
      ggplot2::labs(
        title = "True vs Estimated Function Values ",
        x = "True f(X,Z)",
        y = "Estimated f(X,Z)",
        color = "Treatment"
      ) +
      ggplot2::theme_minimal()
    
    # Create treatment effect plots if both treatment groups are present
    if (length(unique(Z_test)) > 1) {
      # For treatment effect visualization, we need paired observations
      # This works best with a grid of values rather than original data
      
      # Create a grid of covariate values
      if (ncol(X_test) > 1) {
        # For 2D, create a grid
        x1_vals <- seq(min(X_test[,1]), max(X_test[,1]), length.out = sqrt(grid_size))
        x2_vals <- seq(min(X_test[,2]), max(X_test[,2]), length.out = sqrt(grid_size))
        grid_df <- expand.grid(X1 = x1_vals, X2 = x2_vals)
        
        # Calculate true treatment effect
        X_grid <- as.matrix(grid_df)
        grid_df$TrueTrtEffect <- true_f(X_grid, rep(1, nrow(X_grid))) - 
          true_f(X_grid, rep(0, nrow(X_grid)))
        
        # Calculate estimated treatment effect
        # Need to create augmented data for Z=0 and Z=1
        ps_grid <- fit$propensity_forest$do_predict(X_grid)
        e_hat_grid <- 1 / (1 + exp(-ps_grid))
        
        X_aug1 <- cbind(X_grid, 1, e_hat_grid)
        X_aug0 <- cbind(X_grid, 0, e_hat_grid)
        
        f1_est <- fit$forest$do_predict(X_aug1)
        f0_est <- fit$forest$do_predict(X_aug0)
        
        grid_df$EstTrtEffect <- f1_est - f0_est
        
        # Heatmap plot of treatment effect differences
        te_heatmap <- ggplot2::ggplot(grid_df, ggplot2::aes(x = X1, y = X2)) +
          ggplot2::geom_tile(ggplot2::aes(fill = EstTrtEffect - TrueTrtEffect)) +
          ggplot2::scale_fill_gradient2(
            low = "blue", mid = "white", high = "red",
            midpoint = 0, 
            name = "Est. - True"
          ) +
          ggplot2::labs(
            title = "Treatment Effect Estimation Error",
            x = "X1", 
            y = "X2"
          ) +
          ggplot2::theme_minimal()
        
        # 3D surface plot of true treatment effect
        true_te_plot <- ggplot2::ggplot(grid_df, ggplot2::aes(x = X1, y = X2, fill = TrueTrtEffect)) +
          ggplot2::geom_tile() +
          ggplot2::scale_fill_viridis_c(name = "True Effect") +
          ggplot2::labs(
            title = "True Treatment Effect",
            x = "X1", 
            y = "X2"
          ) +
          ggplot2::theme_minimal()
        
        # 3D surface plot of estimated treatment effect
        est_te_plot <- ggplot2::ggplot(grid_df, ggplot2::aes(x = X1, y = X2, fill = EstTrtEffect)) +
          ggplot2::geom_tile() +
          ggplot2::scale_fill_viridis_c(name = "Est. Effect") +
          ggplot2::labs(
            title = "Estimated Treatment Effect",
            x = "X1", 
            y = "X2"
          ) +
          ggplot2::theme_minimal()
        
        results$plots <- list(
          scatter_plot = scatter_plot,
          true_te_plot = true_te_plot,
          est_te_plot = est_te_plot,
          te_heatmap = te_heatmap
        )
      } else {
        # For 1D, create a simple line plot
        x_vals <- seq(min(X_test), max(X_test), length.out = grid_size)
        X_line <- matrix(x_vals, ncol = 1)
        
        # Calculate true treatment effects
        true_te <- true_f(X_line, rep(1, grid_size)) - 
          true_f(X_line, rep(0, grid_size))
        
        # Calculate estimated treatment effects
        ps_line <- fit$propensity_forest$do_predict(X_line)
        e_hat_line <- 1 / (1 + exp(-ps_line))
        
        X_aug1 <- cbind(X_line, 1, e_hat_line)
        X_aug0 <- cbind(X_line, 0, e_hat_line)
        
        f1_est <- fit$forest$do_predict(X_aug1)
        f0_est <- fit$forest$do_predict(X_aug0)
        
        est_te <- f1_est - f0_est
        
        line_data <- data.frame(
          X = x_vals,
          TrueTrtEffect = true_te,
          EstTrtEffect = est_te
        )
        
        te_line <- ggplot2::ggplot(line_data, ggplot2::aes(x = X)) +
          ggplot2::geom_line(ggplot2::aes(y = TrueTrtEffect, color = "True")) +
          ggplot2::geom_line(ggplot2::aes(y = EstTrtEffect, color = "Estimated")) +
          ggplot2::labs(
            title = "Treatment Effect Comparison",
            x = "X",
            y = "Treatment Effect",
            color = "Source"
          ) +
          ggplot2::scale_color_manual(values = c("True" = "black", "Estimated" = "blue")) +
          ggplot2::theme_minimal()
        
        results$plots <- list(
          scatter_plot = scatter_plot,
          te_line = te_line
        )
      }
    } else {
      results$plots <- list(scatter_plot = scatter_plot)
    }
  }
  
  return(results)
}

# Now let's create a function to display the results in a single combined figure
display_combined_results <- function(results) {
  # Check if required packages are installed
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    stop("Package 'gridExtra' is needed. Please install it with install.packages('gridExtra')")
  }
  if (!requireNamespace("grid", quietly = TRUE)) {
    stop("Package 'grid' is needed. Please install it with install.packages('grid')")
  }
  
  # Create metrics text
  metrics_text <- paste(
    "Performance Metrics:",
    sprintf("RMSE: %.4f", results$metrics$rmse),
    sprintf("MAE: %.4f", results$metrics$mae),
    sprintf("Correlation: %.4f", results$metrics$correlation),
    sep = "\n"
  )
  
  # Create a text grob for metrics
  metrics_grob <- grid::textGrob(
    metrics_text,
    gp = grid::gpar(fontface = "bold", fontsize = 11),
    just = "left"
  )
  
  # Combine scatter plot, te_heatmap, and metrics
  if (!is.null(results$plots$te_heatmap)) {
    combined_plot <- gridExtra::grid.arrange(
      results$plots$scatter_plot,
      results$plots$te_heatmap,
      metrics_grob,
      layout_matrix = rbind(c(1, 2), c(3, 3)),
      heights = c(4, 1)
    )
  } else {
    # If no te_heatmap is available, just show scatter plot and metrics
    combined_plot <- gridExtra::grid.arrange(
      results$plots$scatter_plot,
      metrics_grob,
      heights = c(4, 1)
    )
  }
  
  return(combined_plot)
}

