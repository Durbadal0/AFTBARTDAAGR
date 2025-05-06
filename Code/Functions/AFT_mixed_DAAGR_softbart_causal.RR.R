# Install and load necessary packages
if (!requireNamespace("mvtnorm", quietly = TRUE)) {
  install.packages("mvtnorm")
}
if (!requireNamespace("SoftBart", quietly = TRUE)) {
  install.packages("SoftBart")
}
library(mvtnorm)
library(SoftBart)




# Function: Bayesian mixed model with DAGAR (Directed Acyclic Graph Auto-Regressive) 
# spatial random effects and nonparametric f() via SoftBart
AFT_mixed_DAGAR_softbart_causal <- function(time, status, X, Z, group, W, 
                                            X_test, Z_test = Z, group_test = group,
                                            n_iter = 6000, burn_in = 1000, thin = 5,
                                            a_sigma = 1, b_sigma = 1,
                                            a_r = 1, b_r = 1,
                                            rho = 0.8,        # initial value for rho
                                            rho_prop_sd = 0.05,   # proposal sd for rho update
                                            ordering = NULL,  # ordering for DAGAR construction
                                            hypers = NULL, opts = NULL,
                                            n_iter_ps = 2000) {
  # time, status: survival outcome and censoring indicator
  # X: training covariate matrix (n x p)
  # Z: binary treatment indicator for training data
  # group: spatial group indicator for each observation
  # W: spatial adjacency matrix for groups (symmetric, zeros on diagonal)
  # X_test: test set covariate matrix
  # Z_test: test set treatment indicator (default same as training if not provided)
  # group_test: test set group indicator (default same as training)
  # ordering: specify ordering for DAGAR construction (default: sum of coordinates)
  
  n <- length(time)
  p <- ncol(X)
  
  # Transform survival times to the log scale
  y <- log(time)
  
  # --- 0. Estimate the Propensity Score Using SBART ---
  # We use SoftBart to fit a regression on the binary outcome Z.
  # Since SoftBart is a regression routine, we obtain raw predictions and then transform them via a logistic (expit) function.
  expit <- function(x) 1/(1+exp(-x))
  
  hypers_ps <- Hypers(X = X, Y = Z)  # Note: even though Z is binary, we use SBART to obtain a latent score.
  opts_ps <- Opts(num_burn = 0, num_save = 1, num_thin = 1)
  forest_ps <- MakeForest(hypers_ps, opts_ps)
  forest_ps$do_gibbs(X, Z, X, n_iter_ps)
  
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
  
  # --- 2. Prepare for DAGAR Model Construction ---
  
  # Identify unique groups
  groups <- sort(unique(group))
  k <- length(groups)
  
  # Create mapping from original groups to sequential indices
  group_map <- match(group, groups)
  group_test_map <- match(group_test, groups)
  
  # If ordering is not provided, create one based on a reasonable default
  if (is.null(ordering)) {
    # Create a default ordering based on sum of coordinates
    # The code assumes that there is a reasonable way to order the regions
    # For example, if we have coordinates for each region
    ordering <- 1:k
  }
  
  # Create a permutation vector based on the ordering
  perm <- ordering
  inv_perm <- match(1:k, perm)
  
  # --- 3. Construct DAGAR precision matrix components ---
  
  # Initialize the lower triangular matrix L and diagonal matrix F
  L <- diag(k)
  F_diag <- numeric(k)
  
  # For each node (except the first in ordering)
  for (i in 2:k) {
    # Get the position in ordered sequence
    pos <- perm[i]
    
    # Find directed neighbors (neighbors that come before in the ordering)
    N_i <- which(W[pos, ] == 1 & inv_perm < inv_perm[pos])
    n_i <- length(N_i)
    
    if (n_i > 0) {
      # DAGAR coefficients based on equation (2.4) in the paper
      b_i <- rho / (1 + (n_i - 1) * rho^2)
      
      # Fill in the lower triangular matrix
      L[pos, N_i] <- -b_i
    }
    
    # Fill in diagonal matrix F based on equation (2.4)
    F_diag[pos] <- (1 + (n_i - 1) * rho^2) / (1 - rho^2)
  }
  
  # First node has no directed neighbors
  F_diag[perm[1]] <- 1 / (1 - rho^2)
  
  # --- 4. Set Up SoftBart for f(.) ---
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
  
  # Initialize spatial random effects and variance parameters:
  r_current <- rep(0, k)
  sigma2_current <- 1
  rho_current <- rho  # starting value for rho
  
  # Storage for posterior samples:
  n_save <- floor((n_iter - burn_in) / thin)
  f_samples <- matrix(0, nrow = n_save, ncol = n)
  f_test_samples <- matrix(0, nrow = n_save, ncol = nrow(X_test_aug))
  r_samples <- matrix(0, nrow = n_save, ncol = k)
  sigma2_samples <- numeric(n_save)
  rho_samples <- numeric(n_save)
  
  sample_index <- 1
  
  # --- 5. MCMC Sampling ---
  for (iter in 1:n_iter) {
    
    # 5.1 Update DAGAR precision matrix if rho has changed
    if (iter > 1) {
      # Update L and F with current rho
      for (i in 2:k) {
        pos <- perm[i]
        N_i <- which(W[pos, ] == 1 & inv_perm < inv_perm[pos])
        n_i <- length(N_i)
        
        if (n_i > 0) {
          b_i <- rho_current / (1 + (n_i - 1) * rho_current^2)
          L[pos, N_i] <- -b_i
        }
        
        F_diag[pos] <- (1 + (n_i - 1) * rho_current^2) / (1 - rho_current^2)
      }
      
      F_diag[perm[1]] <- 1 / (1 - rho_current^2)
    }
    
    # 5.2 Update f(·) via SoftBart:
    # For each observation, subtract the corresponding spatial effect:
    r_obs <- r_current[group_map]
    resid_for_f <- y - r_obs
    
    # Update the forest by one Gibbs iteration:
    forest$do_gibbs(X_aug, resid_for_f, X_test_aug, 1)
    
    # Update predictions:
    f_current <- forest$do_predict(X_aug)
    f_test_current <- forest$do_predict(X_test_aug)
    
    # 5.3 Update Spatial Random Effects r
    # Construct precision matrix Q = L'FL using current rho
    F_matrix <- diag(F_diag)
    Q <- t(L) %*% F_matrix %*% L
    
    # Calculate the mean for conditional distribution of r
    resid_y <- y - f_current
    
    # Compute precision and mean for the full conditional of r
    Sigma_r_inv <- Q / sigma2_current
    
    # Add data contribution to diagonal of precision matrix
    for (i in 1:k) {
      indices <- which(group_map == i)
      if (length(indices) > 0) {
        Sigma_r_inv[i, i] <- Sigma_r_inv[i, i] + length(indices) / sigma2_current
      }
    }
    
    # Calculate the mean vector
    mu_r <- numeric(k)
    for (i in 1:k) {
      indices <- which(group_map == i)
      if (length(indices) > 0) {
        mu_r[i] <- sum(resid_y[indices]) / sigma2_current
      }
    }
    
    # Solve for mean
    Sigma_r <- solve(Sigma_r_inv)
    mu_r <- Sigma_r %*% mu_r
    
    # Draw sample from multivariate normal
    r_current <- as.vector(rmvnorm(1, mean = mu_r, sigma = Sigma_r))
    
    # 5.4 Update Residual Variance sigma^2:
    r_obs <- r_current[group_map]
    resid <- y - f_current - r_obs
    shape_sigma <- a_sigma + n/2
    rate_sigma <- b_sigma + 0.5 * sum(resid^2)
    sigma2_current <- 1 / rgamma(1, shape = shape_sigma, rate = rate_sigma)
    
    # 5.5 Metropolis-Hastings Update for rho:
    rho_prop <- rho_current + rnorm(1, mean = 0, sd = rho_prop_sd)
    
    if (rho_prop > 0 && rho_prop < 1) {
      # Calculate log-likelihood for proposed rho
      L_prop <- L  # Create copy of L
      F_diag_prop <- F_diag  # Create copy of F_diag
      
      # Update with proposed rho
      for (i in 2:k) {
        pos <- perm[i]
        N_i <- which(W[pos, ] == 1 & inv_perm < inv_perm[pos])
        n_i <- length(N_i)
        
        if (n_i > 0) {
          b_i <- rho_prop / (1 + (n_i - 1) * rho_prop^2)
          L_prop[pos, N_i] <- -b_i
        }
        
        F_diag_prop[pos] <- (1 + (n_i - 1) * rho_prop^2) / (1 - rho_prop^2)
      }
      
      F_diag_prop[perm[1]] <- 1 / (1 - rho_prop^2)
      
      # Compute proposed precision matrix
      F_matrix_prop <- diag(F_diag_prop)
      Q_prop <- t(L_prop) %*% F_matrix_prop %*% L_prop
      
      # Calculate log determinant (using DAGAR properties)
      log_det_prop <- -k * log(1 - rho_prop^2)
      for (i in 2:k) {
        pos <- perm[i]
        N_i <- which(W[pos, ] == 1 & inv_perm < inv_perm[pos])
        n_i <- length(N_i)
        if (n_i > 0) {
          log_det_prop <- log_det_prop - log(1 + (n_i - 1) * rho_prop^2)
        }
      }
      
      log_det_curr <- -k * log(1 - rho_current^2)
      for (i in 2:k) {
        pos <- perm[i]
        N_i <- which(W[pos, ] == 1 & inv_perm < inv_perm[pos])
        n_i <- length(N_i)
        if (n_i > 0) {
          log_det_curr <- log_det_curr - log(1 + (n_i - 1) * rho_current^2)
        }
      }
      
      # Calculate log density for current and proposed
      quad_form_prop <- t(r_current) %*% Q_prop %*% r_current
      quad_form_curr <- t(r_current) %*% Q %*% r_current
      
      log_target_prop <- 0.5 * log_det_prop - 0.5 * quad_form_prop
      log_target_curr <- 0.5 * log_det_curr - 0.5 * quad_form_curr
      
      # Calculate acceptance ratio (uniform prior on rho)
      log_accept_ratio <- log_target_prop - log_target_curr
      
      if (log(runif(1)) < log_accept_ratio) {
        rho_current <- rho_prop
      }
    }
    
    # 5.6 Save Samples After Burn-in and Thinning:
    if (iter > burn_in && ((iter - burn_in) %% thin == 0)) {
      f_samples[sample_index, ] <- f_current
      f_test_samples[sample_index, ] <- f_test_current
      r_samples[sample_index, ] <- r_current
      sigma2_samples[sample_index] <- sigma2_current
      rho_samples[sample_index] <- rho_current
      sample_index <- sample_index + 1
    }
    
    # Display progress
    if (iter %% 500 == 0) {
      cat("Completed iteration", iter, "out of", n_iter, "\n")
      cat("Current rho:", rho_current, "\n")
      cat("Current sigma2:", sigma2_current, "\n")
    }
  }
  
  # --- 6. Return results ---
  return(list(
    f_train = f_samples,
    f_test = f_test_samples,
    r = r_samples,
    sigma2 = sigma2_samples,
    rho = rho_samples,
    forest = forest,
    propensity_forest = forest_ps,
    ordering = perm
  ))
}