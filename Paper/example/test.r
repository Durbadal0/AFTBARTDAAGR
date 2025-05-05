# Load necessary libraries
library(truncnorm)
library(MASS)
library(BART)   # for wbart()

# A helper function that updates f() via a BART call.
# Here, we use wbart() to obtain a single posterior draw of f(x) evaluated at X.
update_BART <- function(X, y, sigma2) {
  # Ensure X is a properly formatted matrix with column names
  if (!is.matrix(X)) X <- as.matrix(X)
  if (is.null(colnames(X))) colnames(X) <- paste0("X", 1:ncol(X))
  
  # Now use the properly formatted matrix
  bart_fit <- wbart(x.train = X, y.train = y, ntree = 50, 
                   ndpost = 1, nskip = 10, sigma = sqrt(sigma2))
  
  # f_hat is the posterior draw for f(x) on the training set
  f_hat <- bart_fit$yhat.train.mean
  return(f_hat)
}

# Gibbs sampler function with BART for f() and joint random effects update.
gibbs_sampler <- function(data, X, clusters, n_iter, hyperparams, W) {
  # data: data frame with columns: y (observed time), Delta (censoring indicator)
  # X: design matrix for the BART model (including any covariates needed for f())
  # clusters: vector of cluster IDs (1-indexed) for each observation.
  # n_iter: number of Gibbs iterations.
  # hyperparams: list with hyperparameters (a_sigma, b_sigma, a_tau, b_tau, 
  #              rho_lower, rho_upper).
  # W: spatial weight matrix for clusters (I x I matrix).
  
  # Unique clusters and their counts.
  cluster_ids <- sort(unique(clusters))
  I <- length(cluster_ids)
  n <- nrow(data)
  
  # Diagonal: row sums of W.
  D <- rowSums(W)
  
  # Pre-calculate indices for each cluster.
  cluster_index_list <- lapply(cluster_ids, function(cl) which(clusters == cl))
  n_cluster <- sapply(cluster_index_list, length)  # number of observations in each cluster
  
  # Initialize parameters.
  sigma2 <- 1       # error variance
  tau2 <- 1         # spatial variance for random effects
  rho <- 0.5        # spatial correlation parameter (must lie in [rho_lower, rho_upper])
  r <- rep(0, I)    # random effects (one per cluster)
  
  # Initialize f() as a vector evaluated at each observation.
  # For instance, start with f(x) = 0 for all observations.
  f_hat <- rep(0, n)
  
  # Initialize latent log survival times.
  # For uncensored observations (Delta==1): set to log(y).
  # For censored ones: initialize to log(y) + 0.5.
  latent_y <- log(data$y)
  latent_y[data$Delta == 0] <- latent_y[data$Delta == 0] + 0.5
  
  # Extract hyperparameters.
  a_sigma <- hyperparams$a_sigma
  b_sigma <- hyperparams$b_sigma
  a_tau   <- hyperparams$a_tau
  b_tau   <- hyperparams$b_tau
  
  # Storage for posterior samples.
  samples <- list(sigma2 = numeric(n_iter),
                  tau2   = numeric(n_iter),
                  rho    = numeric(n_iter),
                  r      = matrix(NA, nrow = n_iter, ncol = I),
                  f_hat  = matrix(NA, nrow = n_iter, ncol = n))
  
  for (iter in 1:n_iter) {
    
    ## Step 1: Impute latent log survival times for censored observations.
    # For each observation, compute mu_i = f_hat[i] + r[cluster].
    for (i in 1:n) {
      cl <- clusters[i]
      mu_i <- f_hat[i] + r[cl]
      if (data$Delta[i] == 0) {
        latent_y[i] <- rtruncnorm(1, a = log(data$y[i]), b = Inf,
                                  mean = mu_i, sd = sqrt(sigma2))
      } else {
        latent_y[i] <- log(data$y[i])  # observed events remain fixed.
      }
    }
    
    ## Step 2: Update f() using BART.
    # The "partial response" for f() is given by: y_star = latent_y - r[cluster].
    y_star <- latent_y - sapply(clusters, function(cl) r[cl])
    f_hat <- update_BART(X, y_star, sigma2)
    
    ## Step 3: Joint update of random effects r.
    # For each cluster, compute s_i = sum_{j in cluster i} (latent_y[j] - f_hat[j]).
    s <- sapply(cluster_index_list, function(idx) sum(latent_y[idx] - f_hat[idx]))
    # The full conditional for r is multivariate normal:
    # Precision matrix: P_r = diag(n_i/sigma2) + (diag(D) - rho * W) / tau2.
    precision_r <- diag(n_cluster / sigma2, I) + (diag(D) - rho * W) / tau2
    Sigma_r <- solve(precision_r)
    mu_r <- Sigma_r %*% (s / sigma2)
    r <- as.vector(mvrnorm(1, mu = mu_r, Sigma = Sigma_r))
    
    ## Step 4: Update sigma2 (error variance) via its inverse-gamma full conditional.
    # Residual: latent_y - (f_hat + r[cluster]).
    resid <- numeric(n)
    for (i in 1:n) {
      resid[i] <- latent_y[i] - (f_hat[i] + r[clusters[i]])
    }
    RSS <- sum(resid^2)
    sigma2 <- 1 / rgamma(1, shape = a_sigma + n/2, rate = b_sigma + RSS/2)
    
    ## Step 5: Update tau2 (spatial variance) via its inverse-gamma full conditional.
    quad_form <- as.numeric(t(r) %*% (diag(D) - rho * W) %*% r)
    tau2 <- 1 / rgamma(1, shape = a_tau + I/2, rate = b_tau + quad_form/2)
    
    ## Step 6: Update rho (spatial correlation) via a Metropolis–Hastings step.
    proposal_sd <- 0.1
    rho_prop <- rnorm(1, mean = rho, sd = proposal_sd)
    rho_lower <- hyperparams$rho_lower
    rho_upper <- hyperparams$rho_upper
    if (rho_prop < rho_lower || rho_prop > rho_upper) {
      rho_new <- rho  # reject if out-of-range.
    } else {
      quad_new <- as.numeric(t(r) %*% (diag(D) - rho_prop * W) %*% r)
      quad_old <- quad_form
      log_accept_ratio <- - (quad_new - quad_old) / (2 * tau2)
      if (log(runif(1)) < log_accept_ratio) {
        rho_new <- rho_prop
      } else {
        rho_new <- rho
      }
    }
    rho <- rho_new
    
    ## Store samples.
    samples$sigma2[iter] <- sigma2
    samples$tau2[iter] <- tau2
    samples$rho[iter] <- rho
    samples$r[iter, ] <- r
    samples$f_hat[iter, ] <- f_hat
  }
  
  return(samples)
}

#############################
# Example: Testing the sampler
#############################

set.seed(123)

# Simulate synthetic data.
# Assume 3 clusters, each with 20 subjects.
n_clusters <- 3
n_per_cluster <- 20
n_total <- n_clusters * n_per_cluster

# True underlying function f(x).
# For illustration, we simulate f(x) using a nonlinear function.
# In practice, f(x) is unknown and estimated via BART.
true_f <- function(x) {
  0.5 * sin(x[,2]) + 0.2 * x[,1]
}

# True parameters.
true_sigma2 <- 0.5
true_tau2 <- 0.2
true_r <- rnorm(n_clusters, mean = 0, sd = sqrt(true_tau2))

# Create cluster assignments.
clusters <- rep(1:n_clusters, each = n_per_cluster)

# Generate covariates X.
# Let X have an intercept and one covariate.
X <- as.matrix(cbind(intercept = 1, x = rnorm(n_total)))

# Generate latent log survival times.
f_true <- true_f(X)
latent_true <- f_true + true_r[clusters] + rnorm(n_total, sd = sqrt(true_sigma2))
# Generate survival times on the original scale.
T_true <- exp(latent_true)

# Introduce censoring: set censoring time at the 70th percentile of T_true (about 30% events).
censoring_time <- quantile(T_true, 0.7)
Delta <- as.numeric(T_true <= censoring_time)
y_obs <- pmin(T_true, censoring_time)

# Create data frame.
data <- data.frame(y = y_obs, Delta = Delta)

# Spatial weight matrix for clusters.
# For simplicity, assume clusters are arranged in a chain.
W <- matrix(0, nrow = n_clusters, ncol = n_clusters)
for (i in 1:n_clusters) {
  if (i > 1) W[i, i-1] <- 1
  if (i < n_clusters) W[i, i+1] <- 1
}
# Diagonal remains zero.

# Set hyperparameters.
hyperparams <- list(
  a_sigma = 2, b_sigma = 1,
  a_tau   = 2, b_tau   = 1,
  rho_lower = 0.0, rho_upper = 1.0
)

# Run the Gibbs sampler.
n_iter <- 2000
samples <- gibbs_sampler(data = data, X = X, clusters = clusters,
                         n_iter = n_iter, hyperparams = hyperparams, W = W)

# Diagnostics: summaries of posterior samples.
cat("Posterior mean of sigma^2:", mean(samples$sigma2), "\n")
cat("True sigma^2:", true_sigma2, "\n\n")
cat("Posterior mean of rho:", mean(samples$rho), "\n")
cat("Posterior mean of tau^2:", mean(samples$tau2), "\n")

# (Optional) Compare the estimated f() to the true f() for the training data.
# Use the last draw of f_hat.
f_est <- samples$f_hat[n_iter, ]
plot(f_true, f_est, main = "Estimated f(x) vs. True f(x)",
     xlab = "True f(x)", ylab = "Estimated f(x)")
abline(a = 0, b = 1, col = "red")