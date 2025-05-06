# Install and load necessary packages
if (!requireNamespace("mvtnorm", quietly = TRUE)) {
  install.packages("mvtnorm")
}
if (!requireNamespace("SoftBart", quietly = TRUE)) {
  install.packages("SoftBart")
}
library(mvtnorm)
library(SoftBart)





AFT_mixed_DAGAR_softbart_improved <- function(time, status, X, group, W, X_test, group_test = group,
                                              n_iter = 6000, burn_in = 1000, thin = 5,
                                              a_sigma = 2, b_sigma = 1,  # More informative priors
                                              a_r = 2, b_r = 1,
                                              rho = 0.5,  # Start at middle of range
                                              rho_prop_sd = 0.1,  # Wider proposal initially
                                              ordering = NULL,
                                              adaptive_mh = TRUE,  # Enable adaptive MH for rho
                                              hypers = NULL, opts = NULL) {
  # Data preparation
  y <- log(time)
  n <- length(y)
  p <- ncol(X)
  
  # Group processing
  groups <- sort(unique(group))
  k <- length(groups)
  group_map <- match(group, groups)
  group_test_map <- match(group_test, groups)
  
  # Ordering setup
  if (is.null(ordering)) {
    ordering <- 1:k
  }
  perm <- ordering
  inv_perm <- match(1:k, perm)
  
  # Function to update DAGAR precision components
  updateDagarComponents <- function(rho_val) {
    L <- diag(k)
    F_diag <- numeric(k)
    
    # First node has no directed neighbors
    F_diag[perm[1]] <- 1 / (1 - rho_val^2)
    
    # For each node (except the first in ordering)
    for (i in 2:k) {
      pos <- perm[i]
      
      # Find directed neighbors properly
      neighbors <- which(W[pos, ] == 1)  # All neighbors
      directed_neighbors <- neighbors[inv_perm[neighbors] < inv_perm[pos]]  # Only previous ones
      n_i <- length(directed_neighbors)
      
      if (n_i > 0) {
        # DAGAR coefficients based on equation (2.4) from the paper
        b_i <- rho_val / (1 + (n_i - 1) * rho_val^2)
        L[pos, directed_neighbors] <- -b_i
      }
      
      # Fill in diagonal matrix F based on equation (2.4)
      F_diag[pos] <- (1 + (n_i - 1) * rho_val^2) / (1 - rho_val^2)
    }
    
    # Calculate precision matrix
    F_matrix <- diag(F_diag)
    Q <- t(L) %*% F_matrix %*% L
    
    return(list(L = L, F_diag = F_diag, Q = Q))
  }
  
  # Initialize DAGAR components with starting rho
  dagar <- updateDagarComponents(rho)
  L <- dagar$L
  F_diag <- dagar$F_diag
  Q <- dagar$Q
  
  # Set up SoftBart forest
  if (is.null(hypers)) {
    hypers <- Hypers(X = X, Y = y)
  }
  if (is.null(opts)) {
    opts <- Opts(num_burn = 0, num_save = 1, num_thin = 1)
  }
  forest <- MakeForest(hypers, opts)
  
  # Initialize parameters
  f_current <- forest$do_predict(X)
  f_test_current <- forest$do_predict(X_test)
  r_current <- rep(0, k)
  sigma2_current <- var(y) / 2  # Better initialization
  sigma_r2_current <- var(y) / 2
  rho_current <- rho
  
  # For adaptive MH for rho
  accept_count <- 0
  adapt_count <- 0
  target_rate <- 0.44  # Roberts and Rosenthal (2001) optimal rate
  
  # Storage for posterior samples
  n_save <- floor((n_iter - burn_in) / thin)
  f_samples <- matrix(0, nrow = n_save, ncol = n)
  f_test_samples <- matrix(0, nrow = n_save, ncol = nrow(X_test))
  r_samples <- matrix(0, nrow = n_save, ncol = k)
  sigma2_samples <- numeric(n_save)
  sigma_r2_samples <- numeric(n_save)
  rho_samples <- numeric(n_save)
  
  sample_index <- 1
  
  # MCMC Sampling
  for (iter in 1:n_iter) {
    # Adjust proposal SD for rho adaptively
    if (adaptive_mh && iter > 100 && iter %% 50 == 0) {
      adapt_count <- adapt_count + 1
      accept_rate <- accept_count / 50
      
      # Adjust proposal standard deviation
      if (accept_rate < target_rate) {
        rho_prop_sd <- max(0.01, rho_prop_sd * 0.9)  # Decrease if acceptance rate is too low
      } else {
        rho_prop_sd <- min(0.5, rho_prop_sd * 1.1)   # Increase if acceptance rate is too high
      }
      
      # Reset counter
      accept_count <- 0
    }
    
    # Update f(·) via SoftBart
    r_obs <- r_current[group_map]
    resid_for_f <- y - r_obs
    
    # Update forest by one Gibbs iteration
    forest$do_gibbs(X, resid_for_f, X_test, 1)
    
    # Get updated predictions
    f_current <- forest$do_predict(X)
    f_test_current <- forest$do_predict(X_test)
    
    # Calculate Q if rho has changed
    if (iter > 1 && abs(rho_current - rho) > 1e-10) {
      dagar <- updateDagarComponents(rho_current)
      L <- dagar$L
      F_diag <- dagar$F_diag
      Q <- dagar$Q
      rho <- rho_current  # Update stored rho value
    }
    
    # Update spatial random effects r
    F_matrix <- diag(F_diag)
    
    # Calculate precision and mean for full conditional of r
    Sigma_r_inv <- Q / sigma_r2_current
    
    # Add data contribution to precision diagonal (clustered observations)
    for (i in 1:k) {
      indices <- which(group_map == i)
      n_i <- length(indices)
      if (n_i > 0) {
        Sigma_r_inv[i, i] <- Sigma_r_inv[i, i] + n_i / sigma2_current
      }
    }
    
    # Calculate mean vector
    mu_r <- numeric(k)
    for (i in 1:k) {
      indices <- which(group_map == i)
      if (length(indices) > 0) {
        mu_r[i] <- sum(y[indices] - f_current[indices]) / sigma2_current
      }
    }
    
    # Use Cholesky decomposition for stable and efficient random effect updates
    tryCatch({
      chol_prec <- chol(Sigma_r_inv)
      z <- rnorm(k)
      # Solve system using Cholesky
      r_current <- backsolve(chol_prec, forwardsolve(t(chol_prec), mu_r) + z)
    }, error = function(e) {
      # Fallback to direct method if Cholesky fails
      Sigma_r <- solve(Sigma_r_inv)
      mu_r <- Sigma_r %*% mu_r
      r_current <<- as.vector(rmvnorm(1, mean = mu_r, sigma = Sigma_r))
    })
    
    # Update residual variance sigma^2
    r_obs <- r_current[group_map]
    resid <- y - f_current - r_obs
    shape_sigma <- a_sigma + n/2
    rate_sigma <- b_sigma + 0.5 * sum(resid^2)
    sigma2_current <- 1 / rgamma(1, shape = shape_sigma, rate = rate_sigma)
    
    # Update spatial variance parameter sigma_r^2
    shape_r <- a_r + k/2
    rate_r <- b_r + 0.5 * as.numeric(t(r_current) %*% Q %*% r_current)
    sigma_r2_current <- 1 / rgamma(1, shape = shape_r, rate = rate_r)
    
    # Metropolis-Hastings update for rho with logit transform
    # logit_rho = log(rho/(1-rho)) transforms (0,1) to (-Inf,Inf)
    logit_rho_current <- log(rho_current / (1 - rho_current))
    logit_rho_prop <- logit_rho_current + rnorm(1, mean = 0, sd = rho_prop_sd)
    rho_prop <- exp(logit_rho_prop) / (1 + exp(logit_rho_prop))  # inverse logit
    
    # Build proposed L and F
    dagar_prop <- updateDagarComponents(rho_prop)
    Q_prop <- dagar_prop$Q
    
    # Calculate log determinants efficiently
    log_det_prop <- 0
    log_det_curr <- 0
    
    # Try using determinant function with logarithm option
    tryCatch({
      log_det_prop <- as.numeric(determinant(Q_prop, logarithm = TRUE)$modulus)
      log_det_curr <- as.numeric(determinant(Q, logarithm = TRUE)$modulus)
    }, error = function(e) {
      # Fallback to determinant calculation from DAGAR paper
      log_det_prop <- -k * log(1 - rho_prop^2)
      log_det_curr <- -k * log(1 - rho_current^2)
      
      # Adjust for directed edges
      for (i in 2:k) {
        pos <- perm[i]
        neighbors <- which(W[pos, ] == 1)
        directed_neighbors <- neighbors[inv_perm[neighbors] < inv_perm[pos]]
        n_i <- length(directed_neighbors)
        
        if (n_i > 0) {
          log_det_prop <- log_det_prop - log(1 + (n_i - 1) * rho_prop^2)
          log_det_curr <- log_det_curr - log(1 + (n_i - 1) * rho_current^2)
        }
      }
    })
    
    # Calculate log density for current and proposed
    quad_form_prop <- t(r_current) %*% Q_prop %*% r_current
    quad_form_curr <- t(r_current) %*% Q %*% r_current
    
    log_target_prop <- 0.5 * log_det_prop - 0.5 * quad_form_prop / sigma_r2_current
    log_target_curr <- 0.5 * log_det_curr - 0.5 * quad_form_curr / sigma_r2_current
    
    # Include Jacobian adjustment for logit transform
    log_jacobian_adj <- log(rho_prop * (1 - rho_prop)) - log(rho_current * (1 - rho_current))
    
    # Calculate acceptance ratio
    log_accept_ratio <- log_target_prop - log_target_curr + log_jacobian_adj
    
    # Accept/reject step
    if (log(runif(1)) < log_accept_ratio) {
      rho_current <- rho_prop
      Q <- Q_prop
      accept_count <- accept_count + 1
    }
    
    # Save samples after burn-in and thinning
    if (iter > burn_in && ((iter - burn_in) %% thin == 0)) {
      f_samples[sample_index, ] <- f_current
      f_test_samples[sample_index, ] <- f_test_current
      r_samples[sample_index, ] <- r_current
      sigma2_samples[sample_index] <- sigma2_current
      sigma_r2_samples[sample_index] <- sigma_r2_current
      rho_samples[sample_index] <- rho_current
      sample_index <- sample_index + 1
    }
    
    # Display progress
    if (iter %% 1000 == 0) {
      cat("Completed iteration", iter, "of", n_iter, "\n")
      cat("Current rho:", rho_current, " sigma2:", sigma2_current, 
          " sigma_r2:", sigma_r2_current, "\n")
      if (adaptive_mh) {
        cat("Current rho proposal SD:", rho_prop_sd, 
            " Recent acceptance rate:", accept_count/(min(iter, 50)), "\n")
      }
    }
  }
  
  # Return results
  return(list(
    f_train = f_samples,
    f_test = f_test_samples,
    r = r_samples,
    sigma2 = sigma2_samples,
    sigma_r2 = sigma_r2_samples,
    rho = rho_samples,
    forest = forest,
    ordering = perm
  ))
}




