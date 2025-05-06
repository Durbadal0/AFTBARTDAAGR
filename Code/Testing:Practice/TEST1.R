# Simulation study for AFT_mixed_DAGAR_softbart function with simpler true function

library(mvtnorm)
library(SoftBart)
library(ggplot2)
library(dplyr)

# Revised DAGAR function
AFT_mixed_DAGAR_softbart <- function(time, status, X, group, W, X_test, group_test = group,
                                     n_iter = 6000, burn_in = 1000, thin = 5,
                                     a_sigma = 1, b_sigma = 1,
                                     a_r = 1, b_r = 1,
                                     rho = 0.8,        # initial value for rho
                                     rho_prop_sd = 0.05,   # proposal sd for rho update
                                     ordering = NULL,  # ordering for DAGAR construction
                                     hypers = NULL, opts = NULL) {
  # time, status: survival outcome and censoring indicator
  # X: covariate matrix (n x p)
  # group: spatial group indicator for each observation
  # W: spatial adjacency matrix (symmetric, zeros on diagonal)
  
  y <- log(time)
  n <- length(y)
  p <- ncol(X)
  
  # Unique groups and indices
  groups <- sort(unique(group))
  k <- length(groups)
  
  # Create group mapping
  group_map <- match(group, groups)
  group_test_map <- match(group_test, groups)
  
  # If ordering is not provided, create a simple one
  if (is.null(ordering)) {
    ordering <- 1:k
  }
  
  # Create permutation vectors
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
      
      # Find directed neighbors properly - neighbors that come before in the ordering
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
    
    return(list(L = L, F_diag = F_diag))
  }
  
  # Initialize DAGAR components with starting rho
  dagar <- updateDagarComponents(rho)
  L <- dagar$L
  F_diag <- dagar$F_diag
  
  # Set up SoftBart forest for f() update
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
  sigma2_current <- 1
  sigma_r2_current <- 1 # Initialize spatial variance
  rho_current <- rho
  
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
    # Update L and F matrices if rho has changed
    if (iter > 1 && abs(rho_current - rho) > 1e-10) {
      dagar <- updateDagarComponents(rho_current)
      L <- dagar$L
      F_diag <- dagar$F_diag
      rho <- rho_current  # Update stored rho value
    }
    
    # Update f(·) via SoftBart
    r_obs <- r_current[group_map]
    resid_for_f <- y - r_obs
    
    # Update forest by one Gibbs iteration
    forest$do_gibbs(X, resid_for_f, X_test, 1)
    
    # Get updated predictions
    f_current <- forest$do_predict(X)
    f_test_current <- forest$do_predict(X_test)
    
    # Update spatial random effects r
    F_matrix <- diag(F_diag)
    Q <- t(L) %*% F_matrix %*% L
    
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
    
    # Solve for mean
    Sigma_r <- solve(Sigma_r_inv)
    mu_r <- Sigma_r %*% mu_r
    
    # Draw sample from multivariate normal
    r_current <- as.vector(rmvnorm(1, mean = mu_r, sigma = Sigma_r))
    
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
    
    # Metropolis-Hastings update for rho
    rho_prop <- rho_current + rnorm(1, mean = 0, sd = rho_prop_sd)
    
    if (rho_prop > 0 && rho_prop < 1) {
      # Build proposed L and F
      dagar_prop <- updateDagarComponents(rho_prop)
      L_prop <- dagar_prop$L
      F_diag_prop <- dagar_prop$F_diag
      
      F_matrix_prop <- diag(F_diag_prop)
      Q_prop <- t(L_prop) %*% F_matrix_prop %*% L_prop
      
      # Calculate log determinants
      log_det_prop <- -k * log(1 - rho_prop^2)
      log_det_curr <- -k * log(1 - rho_current^2)
      
      # Adjust log determinant for the directed graph structure
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
      
      # Calculate log density for current and proposed
      quad_form_prop <- t(r_current) %*% Q_prop %*% r_current
      quad_form_curr <- t(r_current) %*% Q %*% r_current
      
      log_target_prop <- 0.5 * log_det_prop - 0.5 * quad_form_prop / sigma_r2_current
      log_target_curr <- 0.5 * log_det_curr - 0.5 * quad_form_curr / sigma_r2_current
      
      # Calculate acceptance ratio
      log_accept_ratio <- log_target_prop - log_target_curr
      
      if (log(runif(1)) < log_accept_ratio) {
        rho_current <- rho_prop
      }
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

set.seed(12345)

# Simulation parameters
n_clusters <- 36  # Simpler graph structure (6x6 grid)
cluster_sizes <- rpois(n_clusters, lambda = 20)  # Smaller cluster sizes
n_total <- sum(cluster_sizes)  # Total number of observations
n_iter <- 3000    # Fewer iterations for quicker testing
burn_in <- 600
thin <- 5
tau_w <- 4  # Precision for spatial random effects
tau_e <- 2  # Precision for individual-level error
rho_true <- 0.7  # True spatial correlation

cat("Setting up simulation with:\n")
cat("- Number of clusters:", n_clusters, "\n")
cat("- Average cluster size:", mean(cluster_sizes), "\n")
cat("- Total observations:", n_total, "\n")
cat("- MCMC iterations:", n_iter, "\n")

# Create a grid graph structure (6x6 grid)
create_grid_graph <- function(n_row, n_col) {
  n <- n_row * n_col
  W <- matrix(0, n, n)
  
  for (i in 1:n_row) {
    for (j in 1:n_col) {
      node <- (i-1)*n_col + j
      
      # Add east neighbor
      if (j < n_col) W[node, node+1] <- W[node+1, node] <- 1
      
      # Add south neighbor
      if (i < n_row) W[node, node+n_col] <- W[node+n_col, node] <- 1
    }
  }
  return(W)
}

W <- create_grid_graph(6, 6)

# Generate data from DAGAR model with multiple observations per cluster
generate_clustered_data_simple <- function(cluster_sizes, W, rho, tau_w, tau_e) {
  n_clusters <- length(cluster_sizes)
  n_total <- sum(cluster_sizes)
  
  # Create a valid ordering for the DAGAR model
  ordering <- 1:n_clusters
  
  # Construct DAGAR precision matrix based on this ordering
  L <- diag(n_clusters)
  F_diag <- numeric(n_clusters)
  
  # First node has no directed neighbors
  F_diag[1] <- 1 / (1 - rho^2)
  
  # For each subsequent node
  for (i in 2:n_clusters) {
    # Find directed neighbors properly
    neighbors <- which(W[i, ] == 1)
    directed_neighbors <- neighbors[neighbors < i]
    n_i <- length(directed_neighbors)
    
    if (n_i > 0) {
      b_i <- rho / (1 + (n_i - 1) * rho^2)
      L[i, directed_neighbors] <- -b_i
    }
    
    F_diag[i] <- (1 + (n_i - 1) * rho^2) / (1 - rho^2)
  }
  
  # Construct precision matrix
  F_matrix <- diag(F_diag)
  Q <- t(L) %*% F_matrix %*% L
  
  # Generate spatial random effects at the cluster level
  w <- as.vector(rmvnorm(1, rep(0, n_clusters), solve(tau_w * Q)))
  
  # Create group membership vector
  group <- rep(1:n_clusters, cluster_sizes)
  
  # Generate covariates - just 2 covariates for simpler function
  X <- matrix(rnorm(n_total * 2), n_total, 2)
  colnames(X) <- c("X1", "X2")
  
  # Simple true function - linear function
  f_true <- function(X) {
    2 * X[, 1] + X[, 2]
  }
  
  # Generate response with cluster-level random effects
  f_values <- f_true(X)
  y <- f_values + w[group] + rnorm(n_total, 0, 1/sqrt(tau_e))
  
  # Convert to survival time
  time <- exp(y)
  
  # All observations uncensored
  status <- rep(1, n_total)
  
  return(list(
    time = time,
    status = status,
    X = X,
    group = group,
    W = W,
    cluster_sizes = cluster_sizes,
    w_true = w,
    f_true = f_values,
    f_function = f_true,
    rho_true = rho,
    tau_w_true = tau_w,
    tau_e_true = tau_e,
    n_clusters = n_clusters,
    n_total = n_total
  ))
}

# Generate simulated data with simpler function
cat("Generating clustered data with simple linear function...\n")
sim_data <- generate_clustered_data_simple(cluster_sizes, W, rho_true, tau_w, tau_e)

cat("Generated clustered data:\n")
cat("- Number of clusters:", sim_data$n_clusters, "\n")
cat("- Total observations:", sim_data$n_total, "\n")
cat("- Average cluster size:", mean(sim_data$cluster_sizes), "\n")
cat("- Min/Max cluster size:", min(sim_data$cluster_sizes), "/", max(sim_data$cluster_sizes), "\n")

# Fit the model
cat("\nFitting DAGAR model...\n")
fit_dagar <- AFT_mixed_DAGAR_softbart(
  time = sim_data$time,
  status = sim_data$status,
  X = sim_data$X,
  group = sim_data$group,
  W = sim_data$W,
  X_test = sim_data$X,
  group_test = sim_data$group,
  n_iter = n_iter,
  burn_in = burn_in,
  thin = thin,
  rho = 0.5  # Starting value different from true
)

# Evaluate results with detailed diagnostics
evaluate_clustered_results <- function(fit, sim_data) {
  # Calculate posterior means
  f_post_mean <- colMeans(fit$f_train)
  r_post_mean <- colMeans(fit$r)
  rho_post_mean <- mean(fit$rho)
  sigma2_post_mean <- mean(fit$sigma2)
  sigma_r2_post_mean <- mean(fit$sigma_r2)
  
  # Calculate 95% credible intervals
  rho_ci <- quantile(fit$rho, c(0.025, 0.975))
  sigma2_ci <- quantile(fit$sigma2, c(0.025, 0.975))
  sigma_r2_ci <- quantile(fit$sigma_r2, c(0.025, 0.975))
  
  # Recover the fitted times
  fitted_time <- exp(f_post_mean + r_post_mean[sim_data$group])
  
  # Calculate error metrics
  rmse_time <- sqrt(mean((fitted_time - sim_data$time)^2))
  mae_time <- mean(abs(fitted_time - sim_data$time))
  relative_error <- mean(abs((fitted_time - sim_data$time)/sim_data$time))
  
  # Random effects evaluation
  r_rmse <- sqrt(mean((r_post_mean - sim_data$w_true)^2))
  r_mae <- mean(abs(r_post_mean - sim_data$w_true))
  r_corr <- cor(r_post_mean, sim_data$w_true)
  
  # Function evaluation
  f_rmse <- sqrt(mean((f_post_mean - sim_data$f_true)^2))
  f_mae <- mean(abs(f_post_mean - sim_data$f_true))
  f_corr <- cor(f_post_mean, sim_data$f_true)
  
  # Parameter estimates
  rho_error <- abs(rho_post_mean - sim_data$rho_true)
  sigma2_error <- abs(sigma2_post_mean - 1/sim_data$tau_e_true)
  sigma_r2_error <- abs(sigma_r2_post_mean - 1/sim_data$tau_w_true)
  
  # Calculate effective sample size for key parameters
  ess_rho <- coda::effectiveSize(fit$rho)
  ess_sigma2 <- coda::effectiveSize(fit$sigma2)
  ess_sigma_r2 <- coda::effectiveSize(fit$sigma_r2)
  
  # Check if true values are within credible intervals
  rho_coverage <- sim_data$rho_true >= rho_ci[1] && sim_data$rho_true <= rho_ci[2]
  sigma2_coverage <- 1/sim_data$tau_e_true >= sigma2_ci[1] && 1/sim_data$tau_e_true <= sigma2_ci[2]
  sigma_r2_coverage <- 1/sim_data$tau_w_true >= sigma_r2_ci[1] && 1/sim_data$tau_w_true <= sigma_r2_ci[2]
  
  # For random effects, calculate percentage of true values within credible intervals
  r_ci_lower <- apply(fit$r, 2, function(x) quantile(x, 0.025))
  r_ci_upper <- apply(fit$r, 2, function(x) quantile(x, 0.975))
  r_coverage <- mean(sim_data$w_true >= r_ci_lower & sim_data$w_true <= r_ci_upper)
  
  # Print results
  cat("\nModel Evaluation Results:\n")
  cat("\nParameter Estimates:\n")
  cat("  Parameter | True Value | Estimate | 95% CI | Error | Coverage\n")
  cat("  ----------|------------|----------|--------|-------|--------\n")
  cat("  rho       |", sprintf("%.4f", sim_data$rho_true), "|", 
      sprintf("%.4f", rho_post_mean), "|", 
      sprintf("(%.4f, %.4f)", rho_ci[1], rho_ci[2]), "|",
      sprintf("%.4f", rho_error), "|",
      rho_coverage, "\n")
  
  cat("  sigma2    |", sprintf("%.4f", 1/sim_data$tau_e_true), "|", 
      sprintf("%.4f", sigma2_post_mean), "|", 
      sprintf("(%.4f, %.4f)", sigma2_ci[1], sigma2_ci[2]), "|",
      sprintf("%.4f", sigma2_error), "|",
      sigma2_coverage, "\n")
  
  cat("  sigma_r2  |", sprintf("%.4f", 1/sim_data$tau_w_true), "|", 
      sprintf("%.4f", sigma_r2_post_mean), "|", 
      sprintf("(%.4f, %.4f)", sigma_r2_ci[1], sigma_r2_ci[2]), "|",
      sprintf("%.4f", sigma_r2_error), "|",
      sigma_r2_coverage, "\n")
  
  cat("\nEffective Sample Size:\n")
  cat("  rho:", round(ess_rho, 1), "\n")
  cat("  sigma2:", round(ess_sigma2, 1), "\n")
  cat("  sigma_r2:", round(ess_sigma_r2, 1), "\n")
  
  cat("\nRandom Effects Recovery:\n")
  cat("  RMSE:", r_rmse, "\n")
  cat("  MAE:", r_mae, "\n")
  cat("  Correlation:", r_corr, "\n")
  cat("  Coverage Rate:", r_coverage * 100, "%\n")
  
  cat("\nFunction Recovery:\n")
  cat("  RMSE:", f_rmse, "\n")
  cat("  MAE:", f_mae, "\n")
  cat("  Correlation:", f_corr, "\n")
  
  cat("\nTime Prediction:\n")
  cat("  RMSE:", rmse_time, "\n")
  cat("  MAE:", mae_time, "\n")
  cat("  Mean Relative Error:", relative_error * 100, "%\n")
  
  # Create diagnostic plots
  par(mfrow = c(2, 3))
  
  # Plot 1: Trace of rho
  plot(fit$rho, type = "l", main = "Trace Plot for rho", 
       xlab = "Iteration", ylab = "rho")
  abline(h = sim_data$rho_true, col = "red", lty = 2)
  
  # Plot 2: Trace of sigma2
  plot(fit$sigma2, type = "l", main = "Trace Plot for sigma2", 
       xlab = "Iteration", ylab = "sigma2")
  abline(h = 1/sim_data$tau_e_true, col = "red", lty = 2)
  
  # Plot 3: Trace of sigma_r2
  plot(fit$sigma_r2, type = "l", main = "Trace Plot for sigma_r2", 
       xlab = "Iteration", ylab = "sigma_r2")
  abline(h = 1/sim_data$tau_w_true, col = "red", lty = 2)
  
  # Plot 4: True vs. estimated spatial random effects
  plot(sim_data$w_true, r_post_mean, main = "True vs. Estimated Random Effects",
       xlab = "True Random Effects", ylab = "Estimated Random Effects",
       pch = 16, col = "darkblue")
  abline(0, 1, col = "red", lty = 2)
  
  # Plot 5: True vs. estimated function values
  plot(sim_data$f_true, f_post_mean, main = "True vs. Estimated Function",
       xlab = "True Function Values", ylab = "Estimated Function Values",
       pch = 16, col = "darkgreen")
  abline(0, 1, col = "red", lty = 2)
  
  # Plot 6: True vs. fitted survival times
  plot(sim_data$time, fitted_time, main = "True vs. Fitted Times",
       xlab = "True Time", ylab = "Fitted Time",
       pch = 16, col = "purple", cex = 0.5)
  abline(0, 1, col = "red", lty = 2)
  
  par(mfrow = c(1, 1))
  
  # Create distribution plots
  par(mfrow = c(2, 2))
  
  # Plot 1: Posterior distribution of rho
  hist(fit$rho, main = "Posterior Distribution of rho",
       xlab = "rho", freq = FALSE, col = "skyblue")
  abline(v = sim_data$rho_true, col = "red", lty = 2, lwd = 2)
  abline(v = rho_post_mean, col = "blue", lwd = 2)
  legend("topright", legend = c("True", "Estimate"), 
         col = c("red", "blue"), lty = c(2, 1), lwd = 2)
  
  # Plot 2: QQ plot of random effects
  qqplot(sim_data$w_true, r_post_mean, main = "QQ Plot: Random Effects",
         xlab = "True Quantiles", ylab = "Estimated Quantiles")
  abline(0, 1, col = "red", lty = 2)
  
  # Plot 3: Histogram of residuals
  residuals <- fitted_time - sim_data$time
  hist(residuals, main = "Histogram of Time Residuals",
       xlab = "Residuals", freq = FALSE, col = "salmon")
  curve(dnorm(x, mean = mean(residuals), sd = sd(residuals)), 
        add = TRUE, col = "blue", lwd = 2)
  
  # Plot 4: Display cluster sizes vs random effects accuracy
  cluster_accuracy <- data.frame(
    cluster = 1:sim_data$n_clusters,
    size = sim_data$cluster_sizes,
    error = abs(r_post_mean - sim_data$w_true)
  )
  plot(cluster_accuracy$size, cluster_accuracy$error,
       main = "Cluster Size vs. Random Effect Error",
       xlab = "Cluster Size", ylab = "Absolute Error",
       pch = 16, col = "darkblue")
  abline(lm(error ~ size, data = cluster_accuracy), col = "red", lwd = 2)
  
  par(mfrow = c(1, 1))
  
  # Create spatial visualization of random effects
  if (sqrt(n_clusters) == floor(sqrt(n_clusters))) {
    # We have a square grid
    grid_dim <- sqrt(n_clusters)
    
    # Arrange true and estimated random effects on grid
    true_grid <- matrix(sim_data$w_true, nrow = grid_dim, ncol = grid_dim)
    est_grid <- matrix(r_post_mean, nrow = grid_dim, ncol = grid_dim)
    diff_grid <- matrix(r_post_mean - sim_data$w_true, nrow = grid_dim, ncol = grid_dim)
    
    # Create heatmaps
    par(mfrow = c(1, 3))
    
    # Plot 1: True random effects
    image(1:grid_dim, 1:grid_dim, true_grid, 
          main = "True Random Effects",
          xlab = "", ylab = "", axes = FALSE,
          col = heat.colors(100, rev = TRUE))
    box()
    
    # Plot 2: Estimated random effects
    image(1:grid_dim, 1:grid_dim, est_grid, 
          main = "Estimated Random Effects",
          xlab = "", ylab = "", axes = FALSE,
          col = heat.colors(100, rev = TRUE))
    box()
    
    # Plot 3: Differences
    image(1:grid_dim, 1:grid_dim, diff_grid, 
          main = "Differences (Est - True)",
          xlab = "", ylab = "", axes = FALSE,
          col = colorRampPalette(c("blue", "white", "red"))(100))
    box()
    
    par(mfrow = c(1, 1))
  }
  
  # Return evaluation metrics
  return(list(
    parameters = list(
      rho = list(true = sim_data$rho_true, est = rho_post_mean, 
                 ci = rho_ci, error = rho_error, coverage = rho_coverage),
      sigma2 = list(true = 1/sim_data$tau_e_true, est = sigma2_post_mean, 
                    ci = sigma2_ci, error = sigma2_error, coverage = sigma2_coverage),
      sigma_r2 = list(true = 1/sim_data$tau_w_true, est = sigma_r2_post_mean, 
                      ci = sigma_r2_ci, error = sigma_r2_error, coverage = sigma_r2_coverage)
    ),
    ess = list(rho = ess_rho, sigma2 = ess_sigma2, sigma_r2 = ess_sigma_r2),
    random_effects = list(rmse = r_rmse, mae = r_mae, corr = r_corr, coverage = r_coverage),
    function_metrics = list(rmse = f_rmse, mae = f_mae, corr = f_corr),
    time = list(rmse = rmse_time, mae = mae_time, rel_error = relative_error)
  ))
}

# Need coda package for diagnostics
if (!requireNamespace("coda", quietly = TRUE)) {
  install.packages("coda")
}

# Run evaluation
cat("\nEvaluating model fit...\n")
eval_results <- evaluate_clustered_results(fit_dagar, sim_data)

# Summarize validation results
cat("\n=== OVERALL VALIDATION SUMMARY ===\n")
cat("Parameters:\n")
cat("- rho: True =", eval_results$parameters$rho$true, 
    "Estimated =", round(eval_results$parameters$rho$est, 4),
    "95% CI = (", round(eval_results$parameters$rho$ci[1], 4), ",", 
    round(eval_results$parameters$rho$ci[2], 4), ")\n")

cat("- sigma2: True =", eval_results$parameters$sigma2$true, 
    "Estimated =", round(eval_results$parameters$sigma2$est, 4),
    "95% CI = (", round(eval_results$parameters$sigma2$ci[1], 4), ",", 
    round(eval_results$parameters$sigma2$ci[2], 4), ")\n")

cat("- sigma_r2: True =", eval_results$parameters$sigma_r2$true, 
    "Estimated =", round(eval_results$parameters$sigma_r2$est, 4),
    "95% CI = (", round(eval_results$parameters$sigma_r2$ci[1], 4), ",", 
    round(eval_results$parameters$sigma_r2$ci[2], 4), ")\n")

cat("\nRandom Effects Recovery:\n")
cat("- RMSE:", round(eval_results$random_effects$rmse, 4), "\n")
cat("- Correlation:", round(eval_results$random_effects$corr, 4), "\n")
cat("- Coverage Rate:", round(eval_results$random_effects$coverage * 100, 1), "%\n")

cat("\nFunction Recovery:\n")
cat("- RMSE:", round(eval_results$function_metrics$rmse, 4), "\n")
cat("- Correlation:", round(eval_results$function_metrics$corr, 4), "\n")

cat("\nTime Prediction:\n")
cat("- RMSE:", round(eval_results$time$rmse, 4), "\n")
cat("- Mean Relative Error:", round(eval_results$time$rel_error * 100, 1), "%\n")

cat("\nMCMC Diagnostics:\n")
cat("- Effective Sample Size - rho:", round(eval_results$ess$rho, 1), "\n")
cat("- Effective Sample Size - sigma2:", round(eval_results$ess$sigma2, 1), "\n")
cat("- Effective Sample Size - sigma_r2:", round(eval_results$ess$sigma_r2, 1), "\n")

# Final validation assessment
validation_passed <- eval_results$parameters$rho$coverage && 
  eval_results$parameters$sigma2$coverage && 
  eval_results$parameters$sigma_r2$coverage &&
  eval_results$random_effects$corr > 0.8 &&
  eval_results$function_metrics$corr > 0.8 &&
  eval_results$ess$rho > 100 &&
  eval_results$ess$sigma2 > 100 &&
  eval_results$ess$sigma_r2 > 100

if (validation_passed) {
  cat("\nVALIDATION RESULT: PASSED\n")
  cat("The function successfully recovers the true parameters, random effects,\n")
  cat("and function values with good accuracy and chain convergence.\n")
} else {
  cat("\nVALIDATION RESULT: PARTIAL\n")
  cat("Some criteria were not fully met. Check the detailed results above.\n")
}