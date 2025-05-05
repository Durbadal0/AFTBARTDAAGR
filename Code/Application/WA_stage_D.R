# Comprehensive Analysis of Causal Effects with Spatial SBART Model
# ----------------------------------------------------------------------

# Load required packages
library(mvtnorm)
library(SoftBart)
library(dplyr)
library(readr)
library(ggplot2)
library(gridExtra)
library(viridis)
library(knitr)
library(scales)
library(reshape2)
library(stringr)
library(grid)  # For textGrob function

# Source the AFT_mixed_CAR_softbart_causal function
source("Code/Functions/AFT_mixed_CAR_softbart_causal.R")

# Set seed for reproducibility
set.seed(42)

# ----------------------------------------------------------------------
# 1. Data Loading and Preprocessing
# ----------------------------------------------------------------------

# Load the FCR data
fcr_data <- read_csv("Code/Application/Surv_data2.csv")

# Create modified dataset with time+1
fcr_modified <- data.frame(
  time = fcr_data$`as.numeric.date_diff.` + 1,  # Add 1 to avoid zeros
  status = fcr_data$death,
  county = fcr_data$county,
  Age = fcr_data$Age,
  BX_Delay = fcr_data$BX_Delay,
  HR_p = fcr_data$HR_p,
  Tgrade = fcr_data$Tgrade,
  Race_cat = ifelse(fcr_data$Race == 2, "AA", "WA"),
  Stage_cat = factor(fcr_data$Stage, levels = c(1, 2, 3), labels = c("Local", "Regional", "Distant")),
  Z = ifelse(fcr_data$TX_Delay == 1, 0, 1)
)

# Load the adjacency matrix
W_mat <- read_csv("Code/Application/W.mat.csv", col_types = cols(...1 = col_skip()))
W_mat <- as.matrix(W_mat)
storage.mode(W_mat) <- "numeric"

# Select stratum for analysis
race_val <- "WA"
stage_val <- "Distant"
stratum_key <- paste(race_val, stage_val, sep = "_")

# Filter data for this stratum
data_stratum <- fcr_modified %>%
  filter(Race_cat == race_val, Stage_cat == stage_val) %>%
  filter(!is.na(Age), !is.na(BX_Delay), !is.na(HR_p), !is.na(Tgrade), !is.na(Z))

cat("Sample size:", nrow(data_stratum), "observations\n")

# Get unique counties in this stratum
stratum_counties <- sort(unique(data_stratum$county))
cat("Number of counties in stratum:", length(stratum_counties), "\n")

# Extract the subset of the adjacency matrix for these counties
W_subset <- matrix(0, nrow = length(stratum_counties), ncol = length(stratum_counties))
for(i in 1:length(stratum_counties)) {
  for(j in 1:length(stratum_counties)) {
    county_i <- stratum_counties[i]
    county_j <- stratum_counties[j]
    
    # Check if counties are within bounds of W_mat
    if(county_i <= ncol(W_mat) && county_j <= nrow(W_mat)) {
      W_subset[i, j] <- W_mat[county_i, county_j]
    }
  }
}

# Make sure W_subset is symmetric with zeros on diagonal
W_subset <- (W_subset + t(W_subset))/2
diag(W_subset) <- 0

# Function to scale variables
scale01 <- function(x) {
  if(length(unique(x)) <= 1) return(x)
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

# Prepare covariates
X <- as.matrix(data.frame(
  Age = scale01(data_stratum$Age),
  BX_Delay = scale01(data_stratum$BX_Delay),
  HR_p = scale01(data_stratum$HR_p),
  Tgrade = scale01(data_stratum$Tgrade)
))

# Create group indices
group <- match(data_stratum$county, stratum_counties)

# Extract other variables
Z <- data_stratum$Z
time <- data_stratum$time
status <- data_stratum$status

# ----------------------------------------------------------------------
# 2. Model Fitting
# ----------------------------------------------------------------------

# Define MCMC parameters
n_iter <- 5000   # Increased for better convergence
burn_in <- 2000
thin <- 5

cat("Fitting model with spatial component...\n")

# Always fit the model (don't save/load RDS)
cat("Fitting new model...\n")
fit <- AFT_mixed_CAR_softbart_causal(
  time = time,
  status = status,
  X = X,
  Z = Z,
  group = group,
  W = W_subset,
  X_test = X,
  Z_test = Z,
  group_test = group,
  n_iter = n_iter,
  burn_in = burn_in,
  thin = thin,
  n_iter_ps = 5000  # Increased for better propensity score estimation
)

cat("Model fitting complete\n")

# --------------------------------------------------------------------------
# 3. Calculate County-specific Average Treatment Effects with Credible Intervals
# --------------------------------------------------------------------------

#' Calculate county-specific average treatment effects with credible intervals
#'
#' @param fit The fitted AFT_mixed_CAR_softbart_causal model
#' @param X Covariate matrix used in the model
#' @param Z Treatment vector (1=treated, 0=control)
#' @param group County group indices for each observation
#' @param counties Vector of county IDs
#'
#' @return A data frame with county treatment effects and credible intervals
#'
calculate_county_effects_with_ci <- function(fit, X, Z, group, counties) {
  # Helper function for logistic transformation
  expit <- function(x) 1/(1+exp(-x))
  
  # Extract number of MCMC samples
  n_samples <- nrow(fit$r)
  
  # Calculate sample sizes for each county
  county_sizes <- sapply(1:length(counties), function(i) sum(group == i))
  
  # Create results data frame
  results <- data.frame(
    county = counties,
    sample_size = county_sizes,
    mu_avg = NA_real_,
    mu_avg_lower = NA_real_,
    mu_avg_upper = NA_real_,
    mu_sum = NA_real_,
    mu_sum_lower = NA_real_,
    mu_sum_upper = NA_real_
  )
  
  # Initialize matrices to store MCMC samples of effects
  county_avg_samples <- matrix(NA, nrow = n_samples, ncol = length(counties))
  county_sum_samples <- matrix(NA, nrow = n_samples, ncol = length(counties))
  
  # For each MCMC sample
  cat("Processing MCMC samples...\n")
  for(m in 1:n_samples) {
    if(m %% 100 == 0) cat("  Sample", m, "of", n_samples, "\n")
    
    # Extract parameters for this sample
    r_m <- fit$r[m,]
    sigma2_m <- fit$sigma2[m]
    
    # For each county
    for(i in 1:length(counties)) {
      # Skip counties with no patients
      if(county_sizes[i] == 0) {
        next
      }
      
      # Get spatial effect for this county
      r_i_m <- r_m[i]
      
      # Get indices of patients in this county
      county_indices <- which(group == i)
      n_i <- length(county_indices)
      
      # Initialize vector to store individual treatment effects
      mu_ij <- numeric(n_i)
      
      # For each patient in the county
      for(j_idx in 1:n_i) {
        j <- county_indices[j_idx]
        
        # Extract patient covariates
        X_j <- X[j,, drop = FALSE]
        
        # Get estimated propensity score
        ps_raw <- fit$propensity_forest$do_predict(X_j)
        e_hat_j <- expit(ps_raw)
        
        # Create augmented matrices for treatment and control
        X_aug1 <- cbind(X_j, 1, e_hat_j)  # Treatment
        X_aug0 <- cbind(X_j, 0, e_hat_j)  # Control
        
        # Get function predictions
        f1 <- fit$forest$do_predict(X_aug1)
        f0 <- fit$forest$do_predict(X_aug0)
        
        # Calculate potential outcomes
        # E[T] = exp(mu + sigma^2/2) for lognormal
        mu1 <- f1 + r_i_m
        mu0 <- f0 + r_i_m
        
        E_T1 <- exp(mu1 + 0.5*sigma2_m)
        E_T0 <- exp(mu0 + 0.5*sigma2_m)
        
        # Calculate individual treatment effect
        mu_ij[j_idx] <- E_T1 - E_T0
      }
      
      # Calculate sum and average
      mu_i_sum <- sum(mu_ij)
      mu_i_avg <- mu_i_sum / n_i
      
      # Store results for this MCMC sample
      county_sum_samples[m, i] <- mu_i_sum
      county_avg_samples[m, i] <- mu_i_avg
    }
  }
  
  # Calculate posterior means and credible intervals
  for(i in 1:length(counties)) {
    # Skip counties with no patients
    if(county_sizes[i] == 0) {
      next
    }
    
    # Extract samples for this county
    avg_samples <- county_avg_samples[, i]
    sum_samples <- county_sum_samples[, i]
    
    # Calculate mean and credible intervals
    results$mu_avg[i] <- mean(avg_samples, na.rm = TRUE)
    ci_avg <- quantile(avg_samples, probs = c(0.025, 0.975), na.rm = TRUE)
    results$mu_avg_lower[i] <- ci_avg[1]
    results$mu_avg_upper[i] <- ci_avg[2]
    
    results$mu_sum[i] <- mean(sum_samples, na.rm = TRUE)
    ci_sum <- quantile(sum_samples, probs = c(0.025, 0.975), na.rm = TRUE)
    results$mu_sum_lower[i] <- ci_sum[1]
    results$mu_sum_upper[i] <- ci_sum[2]
  }
  
  # Calculate overall ACE and credible interval
  overall_ace_samples <- colMeans(county_avg_samples, na.rm = TRUE)
  overall_ace <- mean(overall_ace_samples, na.rm = TRUE)
  overall_ace_ci <- quantile(overall_ace_samples, probs = c(0.025, 0.975), na.rm = TRUE)
  
  # Add ACE results to the output
  attr(results, "ACE") <- overall_ace
  attr(results, "ACE_lower") <- overall_ace_ci[1]
  attr(results, "ACE_upper") <- overall_ace_ci[2]
  
  # Return results
  return(results)
}

#' Create a plot of county-specific treatment effects against county size
#'
#' @param county_effects The output from calculate_county_effects_with_ci
#' @param title Plot title
#' @param effect_type Type of effect to plot (either "avg" or "sum")
#' @param highlight_counties Optional vector of county IDs to highlight
#' @param highlight_labels Whether to add labels to highlighted counties
#'
#' @return A ggplot object
#'
plot_county_effects <- function(county_effects, 
                                title = "County-Specific Treatment Effects vs. County Size",
                                effect_type = "avg",
                                highlight_counties = NULL,
                                highlight_labels = TRUE) {
  require(ggplot2)
  
  # Create a copy of the data
  plot_data <- county_effects
  
  # Remove counties with missing values or zero sample size
  plot_data <- plot_data[!is.na(plot_data$mu_avg) & plot_data$sample_size > 0, ]
  
  # Choose columns based on effect type
  if(effect_type == "avg") {
    y_col <- "mu_avg"
    y_lower_col <- "mu_avg_lower"
    y_upper_col <- "mu_avg_upper"
    y_label <- "Average Treatment Effect (days)"
  } else if(effect_type == "sum") {
    y_col <- "mu_sum"
    y_lower_col <- "mu_sum_lower"
    y_upper_col <- "mu_sum_upper"
    y_label <- "Sum of Treatment Effects (days)"
  } else {
    stop("effect_type must be either 'avg' or 'sum'")
  }
  
  # Get overall ACE
  ace <- attr(county_effects, "ACE")
  ace_lower <- attr(county_effects, "ACE_lower")
  ace_upper <- attr(county_effects, "ACE_upper")
  
  # Create base plot
  p <- ggplot(plot_data, aes_string(x = "sample_size", y = y_col)) +
    geom_point(alpha = 0.7, size = 3) +
    geom_errorbar(aes_string(ymin = y_lower_col, ymax = y_upper_col), width = 0) +
    geom_hline(yintercept = ace, linetype = "dashed", color = "red") +
    geom_hline(yintercept = ace_lower, linetype = "dotted", color = "red") +
    geom_hline(yintercept = ace_upper, linetype = "dotted", color = "red") +
    labs(
      title = title,
      subtitle = paste0("Overall ACE: ", round(ace, 2), " days (95% CI: ", 
                        round(ace_lower, 2), " to ", round(ace_upper, 2), ")"),
      x = "County Sample Size",
      y = y_label
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold")
    )
  
  # Highlight specific counties if requested
  if(!is.null(highlight_counties)) {
    # Create highlight data
    highlight_data <- plot_data[plot_data$county %in% highlight_counties, ]
    
    if(nrow(highlight_data) > 0) {
      # Add highlighted points
      p <- p + geom_point(data = highlight_data, 
                          aes_string(x = "sample_size", y = y_col),
                          color = "red", size = 4)
      
      # Add labels if requested
      if(highlight_labels) {
        # Check if ggrepel is available
        if(requireNamespace("ggrepel", quietly = TRUE)) {
          p <- p + ggrepel::geom_text_repel(
            data = highlight_data,
            aes_string(x = "sample_size", y = y_col, label = "county"),
            box.padding = 0.5,
            point.padding = 0.3
          )
        } else {
          # Fall back to regular text if ggrepel not available
          p <- p + geom_text(
            data = highlight_data,
            aes_string(x = "sample_size", y = y_col, label = "county"),
            nudge_x = 5,
            nudge_y = 5
          )
        }
      }
    }
  }
  
  return(p)
}

# First, classify Florida counties as rural or urban
classify_florida_counties <- function() {
  # Based on Florida Department of Health criteria
  rural_counties <- c(
    "Baker", "Bradford", "Calhoun", "Columbia", "DeSoto", "Dixie", "Franklin", 
    "Gadsden", "Gilchrist", "Glades", "Gulf", "Hamilton", "Hardee", "Hendry", 
    "Holmes", "Jackson", "Jefferson", "Lafayette", "Levy", "Liberty", "Madison", 
    "Monroe", "Okeechobee", "Putnam", "Suwannee", "Taylor", "Union", "Wakulla", 
    "Walton", "Washington"
  )
  
  all_counties <- c(
    "Alachua", "Baker", "Bay", "Bradford", "Brevard", "Broward", "Calhoun", 
    "Charlotte", "Citrus", "Clay", "Collier", "Columbia", "DeSoto", "Dixie", 
    "Duval", "Escambia", "Flagler", "Franklin", "Gadsden", "Gilchrist", 
    "Glades", "Gulf", "Hamilton", "Hardee", "Hendry", "Hernando", "Highlands", 
    "Hillsborough", "Holmes", "Indian River", "Jackson", "Jefferson", 
    "Lafayette", "Lake", "Lee", "Leon", "Levy", "Liberty", "Madison", 
    "Manatee", "Marion", "Martin", "Miami-Dade", "Monroe", "Nassau", 
    "Okaloosa", "Okeechobee", "Orange", "Osceola", "Palm Beach", "Pasco", 
    "Pinellas", "Polk", "Putnam", "St. Johns", "St. Lucie", "Santa Rosa", 
    "Sarasota", "Seminole", "Sumter", "Suwannee", "Taylor", "Union", 
    "Volusia", "Wakulla", "Walton", "Washington"
  )
  
  county_classification <- data.frame(
    county_name = all_counties,
    county_id = 1:67,
    rural_urban = ifelse(all_counties %in% rural_counties, "Rural", "Urban")
  )
  
  return(county_classification)
}

# Function to classify counties based on their numeric IDs
classify_counties_by_id <- function(county_ids) {
  florida_counties <- classify_florida_counties()
  
  county_classification <- data.frame(
    county = county_ids,
    rural_urban = "Unknown"
  )
  
  for(i in 1:length(county_ids)) {
    idx <- which(florida_counties$county_id == county_ids[i])
    if(length(idx) > 0) {
      county_classification$rural_urban[i] <- florida_counties$rural_urban[idx]
    }
  }
  
  return(county_classification)
}

# Calculate county effects with credible intervals
county_effects <- calculate_county_effects_with_ci(
  fit = fit,                # Your fitted model
  X = X,                    # Covariate matrix used in the model
  Z = Z,                    # Treatment vector
  group = group,            # County group indices
  counties = stratum_counties  # Vector of county IDs
)

# Print a summary of the results
cat("Average Causal Effect (ACE):", attr(county_effects, "ACE"), 
    "days (95% CI:", attr(county_effects, "ACE_lower"), "to", 
    attr(county_effects, "ACE_upper"), ")\n")

# Identify top counties by effect size
top_counties <- head(county_effects[order(-county_effects$mu_avg), "county"], 3)
print(paste("Top counties by effect size:", paste(top_counties, collapse=", ")))

# Get rural/urban classification for counties
county_class <- classify_counties_by_id(stratum_counties)

# Add rural/urban classification to county_effects
county_effects$rural_urban <- county_class$rural_urban

# Create modified plot function for rural/urban coloring
plot_county_effects_rural_urban <- function(county_effects, 
                                            title = "County-Specific Treatment Effects vs. County Size",
                                            effect_type = "avg",
                                            highlight_counties = NULL,
                                            highlight_labels = TRUE) {
  require(ggplot2)
  
  # Create a copy of the data
  plot_data <- county_effects
  
  # Remove counties with missing values or zero sample size
  plot_data <- plot_data[!is.na(plot_data$mu_avg) & plot_data$sample_size > 0, ]
  
  # Choose columns based on effect type
  if(effect_type == "avg") {
    y_col <- "mu_avg"
    y_lower_col <- "mu_avg_lower"
    y_upper_col <- "mu_avg_upper"
    y_label <- "Days"
  } else if(effect_type == "sum") {
    y_col <- "mu_sum"
    y_lower_col <- "mu_sum_lower"
    y_upper_col <- "mu_sum_upper"
    y_label <- "Sum of Treatment Effects (days)"
  } else {
    stop("effect_type must be either 'avg' or 'sum'")
  }
  
  # Get overall ACE
  ace <- attr(county_effects, "ACE")
  ace_lower <- attr(county_effects, "ACE_lower")
  ace_upper <- attr(county_effects, "ACE_upper")
  
  # Calculate quartiles for rural and urban counties separately
  rural_data <- plot_data[plot_data$rural_urban == "Rural", ]
  urban_data <- plot_data[plot_data$rural_urban == "Urban", ]
  
  # Rural quartiles
  rural_q1 <- quantile(rural_data[[y_col]], 0.25, na.rm = TRUE)
  rural_q2 <- quantile(rural_data[[y_col]], 0.50, na.rm = TRUE) # median
  rural_q3 <- quantile(rural_data[[y_col]], 0.75, na.rm = TRUE)
  
  # Urban quartiles
  urban_q1 <- quantile(urban_data[[y_col]], 0.25, na.rm = TRUE)
  urban_q2 <- quantile(urban_data[[y_col]], 0.50, na.rm = TRUE) # median
  urban_q3 <- quantile(urban_data[[y_col]], 0.75, na.rm = TRUE)
  
  # Create base plot with rural/urban coloring
  p <- ggplot(plot_data, aes_string(x = "sample_size", y = y_col, color = "rural_urban")) +
    # Zero line (bold black dashed)
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1) +
    
    # Rural quartile lines
    geom_hline(yintercept = rural_q1, linetype = "dotted", color = "forestgreen", size = 0.7) +
    geom_hline(yintercept = rural_q2, linetype = "solid", color = "forestgreen", size = 0.7) +
    geom_hline(yintercept = rural_q3, linetype = "dotted", color = "forestgreen", size = 0.7) +
    
    # Urban quartile lines
    geom_hline(yintercept = urban_q1, linetype = "dotted", color = "purple", size = 0.7) +
    geom_hline(yintercept = urban_q2, linetype = "solid", color = "purple", size = 0.7) +
    geom_hline(yintercept = urban_q3, linetype = "dotted", color = "purple", size = 0.7) +
    
    # Overall ACE lines
    geom_hline(yintercept = ace, linetype = "dashed", color = "red") +
    geom_hline(yintercept = ace_lower, linetype = "dotted", color = "red") +
    geom_hline(yintercept = ace_upper, linetype = "dotted", color = "red") +
    
    # Points and error bars
    geom_point(alpha = 0.8, size = 3) +
    geom_errorbar(aes_string(ymin = y_lower_col, ymax = y_upper_col), width = 0, alpha = 0.7) +
    
    # Labels and aesthetics
    labs(
      title = title,
      subtitle = paste0("Overall ACE: ", round(ace, 2), " days (95% CI: ", 
                        round(ace_lower, 2), " to ", round(ace_upper, 2), ")"),
      x = "County Sample Size",
      y = y_label,
      color = "County Type"
    ) +
    scale_color_manual(
      values = c("Rural" = "forestgreen", "Urban" = "purple", "Unknown" = "gray50")
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      legend.position = "bottom"
    ) +
    
    # Add annotations for quartile lines
    annotate("text", x = max(plot_data$sample_size) * 0.95, y = rural_q1, 
             label = paste0("Rural Q1: ", round(rural_q1, 1)), 
             color = "forestgreen", hjust = 1, size = 3) +
    annotate("text", x = max(plot_data$sample_size) * 0.95, y = rural_q2, 
             label = paste0("Rural Q2: ", round(rural_q2, 1)), 
             color = "forestgreen", hjust = 1, size = 3) +
    annotate("text", x = max(plot_data$sample_size) * 0.95, y = rural_q3, 
             label = paste0("Rural Q3: ", round(rural_q3, 1)), 
             color = "forestgreen", hjust = 1, size = 3) +
    
    annotate("text", x = max(plot_data$sample_size) * 0.95, y = urban_q1, 
             label = paste0("Urban Q1: ", round(urban_q1, 1)), 
             color = "purple", hjust = 1, size = 3) +
    annotate("text", x = max(plot_data$sample_size) * 0.95, y = urban_q2, 
             label = paste0("Urban Q2: ", round(urban_q2, 1)), 
             color = "purple", hjust = 1, size = 3) +
    annotate("text", x = max(plot_data$sample_size) * 0.95, y = urban_q3, 
             label = paste0("Urban Q3: ", round(urban_q3, 1)), 
             color = "purple", hjust = 1, size = 3) +
    
    annotate("text", x = max(plot_data$sample_size) * 0.95, y = 0, 
             label = "Zero Effect", 
             color = "black", hjust = 1, size = 3)
  
  return(p)
}

# Create plot of average treatment effects vs. county size with rural/urban coloring
p_avg <- plot_county_effects_rural_urban(
  county_effects = county_effects,
  title = "County-Specific CACE",
  effect_type = "avg",
  highlight_counties = NULL
)

# Display the plot
print(p_avg)

# Save the plot
ggsave("county_effects_plot.png", p_avg, width = 10, height = 8)

# Create Florida county map showing Life Years Gained
# Add this code after calculating county_effects

# Load required packages for mapping
library(maps)
library(sf)

# Get Florida county map data properly
florida_map <- map_data("county", "florida")

# Create a lookup table for county IDs
# The county IDs are in alphabetical order
fl_counties <- c(
  "alachua", "baker", "bay", "bradford", "brevard", "broward", "calhoun", 
  "charlotte", "citrus", "clay", "collier", "columbia", "desoto", "dixie", 
  "duval", "escambia", "flagler", "franklin", "gadsden", "gilchrist", 
  "glades", "gulf", "hamilton", "hardee", "hendry", "hernando", "highlands", 
  "hillsborough", "holmes", "indian river", "jackson", "jefferson", 
  "lafayette", "lake", "lee", "leon", "levy", "liberty", "madison", 
  "manatee", "marion", "martin", "miami-dade", "monroe", "nassau", 
  "okaloosa", "okeechobee", "orange", "osceola", "palm beach", "pasco", 
  "pinellas", "polk", "putnam", "santa rosa", "sarasota", "seminole", 
  "st. johns", "st. lucie", "sumter", "suwannee", "taylor", "union", 
  "volusia", "wakulla", "walton", "washington"
)

# Create county ID lookup
county_lookup <- data.frame(
  county = 1:67,
  county_name = fl_counties,
  stringsAsFactors = FALSE
)

# Prepare data for mapping
# Create a full dataset with all counties
all_counties <- data.frame(
  county = 1:67
)

# Join our county effects with all counties
map_data <- left_join(all_counties, county_effects, by = "county")

# Handle missing values
map_data$mu_avg[is.na(map_data$mu_avg)] <- 0
map_data$mu_avg_lower[is.na(map_data$mu_avg_lower)] <- 0
map_data$mu_avg_upper[is.na(map_data$mu_avg_upper)] <- 0

# Join county names
map_data <- left_join(map_data, county_lookup, by = "county")

# Convert life years to days (ACE values are in days)
map_data$lyg_years <- map_data$mu_avg / 365.25

# Join to Florida map data
florida_map$county_name <- florida_map$subregion
florida_map_data <- left_join(florida_map, 
                              map_data,
                              by = "county_name")

# Create a map showing Life Years Gained
lyg_map <- ggplot() +
  geom_polygon(data = florida_map_data, 
               aes(x = long, y = lat, group = group, fill = lyg_years),
               color = "white", size = 0.2) +
  coord_fixed(1.3) +
  scale_fill_viridis_c(
    option = "plasma", 
    name = "Years",
    na.value = "gray80"
  ) +
  labs(
    title = "County-Specific ACE",
    
    caption = "Counties in gray have no data"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12),
    legend.position = "right"
  )

print(lyg_map)

# Save the map
ggsave("florida_lyg_map.png", lyg_map, width = 10, height = 8)

# Create a map showing the statistical significance
# Counties where the confidence interval doesn't include 0 are significant
florida_map_data$significant <- with(florida_map_data, 
                                     (mu_avg_lower > 0) | (mu_avg_upper < 0))

# Create significance map
sig_map <- ggplot() +
  geom_polygon(data = florida_map_data, 
               aes(x = long, y = lat, group = group, 
                   fill = significant),
               color = "white", size = 0.2) +
  coord_fixed(1.3) +
  scale_fill_manual(
    values = c("FALSE" = "lightgray", "TRUE" = "darkred"),
    name = "Statistically\nSignificant",
    labels = c("No", "Yes"),
    na.value = "gray80"
  ) +
  labs(
    title = "Counties with Significant Treatment Effect",
    subtitle = paste("WA", stage_val, "- 95% Credible Interval excludes zero"),
    caption = "Counties in gray have non-significant effect or no data"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12),
    legend.position = "right"
  )

print(sig_map)

# Save the significance map
ggsave("florida_significance_map.png", sig_map, width = 10, height = 8)

# Create a combined plot showing both maps side by side
combined_maps <- grid.arrange(lyg_map, sig_map, ncol = 2,
                              top = textGrob("Life Years Gained from Treatment by County",
                                             gp = gpar(fontsize = 20, font = 2)))

# Save the combined plot
ggsave("florida_combined_maps.png", combined_maps, width = 16, height = 8)

# Print a message
cat("Florida county maps created and saved.\n")

# --------------------------------------------------------------------------
# 4. Calculate CSPCE (Conditional Survival Probability Causal Effect) for each county
# --------------------------------------------------------------------------

# Create time grid from 0 to 14 years (in days)
time_grid <- seq(1, 14*365.25, length.out = 100)
time_years <- time_grid/365.25  # Convert to years for plotting

#' Calculate Conditional Survival Probability Causal Effect (CSPCE) for a county using all patients
#' @param county_idx Index of the county
#' @param X Covariate matrix for all patients in the county
#' @param Z Treatment vector for all patients in the county
#' @param fit The fitted model
#' @param time_grid Time points at which to evaluate CSPCE
#' @return List with CSPCE values and credible intervals
#'
calculate_county_spce_all_patients <- function(county_idx, X_county, fit, time_grid) {
  # Get number of patients in this county
  n_patients <- nrow(X_county)
  
  # Helper function for logistic transformation
  expit <- function(x) 1/(1+exp(-x))
  
  # Get MCMC samples for interval calculation
  n_samples <- nrow(fit$r)
  n_times <- length(time_grid)
  
  # Initialize array to store SPCE for each patient, time point, and MCMC sample
  spce_samples_patients <- array(NA, dim = c(n_samples, n_times, n_patients))
  
  # For each patient in the county
  for (p in 1:n_patients) {
    # Get patient covariates
    X_p <- X_county[p, , drop = FALSE]
    
    # Get estimated propensity score
    ps_raw <- fit$propensity_forest$do_predict(X_p)
    e_hat_p <- expit(ps_raw)
    
    # Create augmented matrices for treatment and control
    X_aug1 <- cbind(X_p, 1, e_hat_p)  # Treatment
    X_aug0 <- cbind(X_p, 0, e_hat_p)  # Control
    
    # For each MCMC sample
    for (m in 1:n_samples) {
      # Get function predictions for this MCMC sample
      f1 <- fit$forest$do_predict(X_aug1)
      f0 <- fit$forest$do_predict(X_aug0)
      
      # Get spatial effect for this sample
      r_m <- fit$r[m, county_idx]
      
      # Get sigma for this sample
      sigma_m <- sqrt(fit$sigma2[m])
      
      # Calculate survival probabilities for each time point
      for (t in 1:n_times) {
        log_t <- log(time_grid[t])
        
        # S(t) = 1 - Φ((log(t) - (f + r))/sigma)
        S1_mt <- 1 - pnorm((log_t - (f1 + r_m)) / sigma_m)
        S0_mt <- 1 - pnorm((log_t - (f0 + r_m)) / sigma_m)
        
        # Store SPCE sample for this patient
        spce_samples_patients[m, t, p] <- S1_mt - S0_mt
      }
    }
  }
  
  # Average over patients for each MCMC sample and time point
  spce_samples <- apply(spce_samples_patients, c(1, 2), mean)
  
  # Calculate mean and credible intervals for SPCE
  spce_mean <- colMeans(spce_samples)
  spce_lower <- apply(spce_samples, 2, function(x) quantile(x, 0.025))
  spce_upper <- apply(spce_samples, 2, function(x) quantile(x, 0.975))
  
  # Return results
  return(list(
    time = time_grid,
    spce = spce_mean,
    spce_lower = spce_lower,
    spce_upper = spce_upper
  ))
}

# Calculate CSPCE for each county using all patients in the county
county_spce_results <- list()
valid_counties <- c()

for (i in 1:length(stratum_counties)) {
  # Get indices of patients in this county
  county_indices <- which(group == i)
  
  # Skip if county has too few patients
  if (length(county_indices) <= 1) {
    next
  }
  
  cat("Calculating SPCE for county", stratum_counties[i], "with", length(county_indices), "patients...\n")
  
  # Get covariates for all patients in this county
  X_county <- X[county_indices, , drop = FALSE]
  
  # Calculate SPCE for this county using all patients
  county_spce <- calculate_county_spce_all_patients(i, X_county, fit, time_grid)
  
  # Store results
  county_spce_results[[as.character(stratum_counties[i])]] <- county_spce
  valid_counties <- c(valid_counties, stratum_counties[i])
}

# Create a data frame for plotting county-specific SPCE
county_spce_df <- data.frame()

for (county in valid_counties) {
  county_data <- county_spce_results[[as.character(county)]]
  
  # Add to the data frame
  county_df <- data.frame(
    county = as.character(county),
    time = county_data$time / 365.25,  # Convert to years
    spce = county_data$spce,
    spce_lower = county_data$spce_lower,
    spce_upper = county_data$spce_upper
  )
  
  county_spce_df <- rbind(county_spce_df, county_df)
}

# Calculate Marginal SPCE by averaging over all counties
marginal_spce <- data.frame(
  time = time_years,
  spce = numeric(length(time_years)),
  spce_lower = numeric(length(time_years)),
  spce_upper = numeric(length(time_years))
)

# For each time point
for (t in 1:length(time_years)) {
  # Get all SPCE values at this time point
  spce_values <- county_spce_df$spce[county_spce_df$time == time_years[t]]
  spce_lower_values <- county_spce_df$spce_lower[county_spce_df$time == time_years[t]]
  spce_upper_values <- county_spce_df$spce_upper[county_spce_df$time == time_years[t]]
  
  # Calculate mean and confidence intervals
  marginal_spce$spce[t] <- mean(spce_values, na.rm = TRUE)
  marginal_spce$spce_lower[t] <- mean(spce_lower_values, na.rm = TRUE)
  marginal_spce$spce_upper[t] <- mean(spce_upper_values, na.rm = TRUE)
}

# Create a plot for county-specific SPCE curves
# Select a few notable counties for better visualization
notable_counties <- c(
  as.character(which.max(county_effects$mu_avg)),  # County with highest effect
  as.character(which.min(county_effects$mu_avg)),  # County with lowest effect
  "6"                                       # Broward county (assuming ID is 6)
)

# Filter for notable counties
notable_counties_df <- county_spce_df %>%
  filter(county %in% notable_counties)

# Create a color palette for counties
county_colors <- rainbow(length(notable_counties))

# Plot county-specific SPCE
p_county_spce <- ggplot() +
  # Add confidence intervals for each county
  geom_ribbon(data = notable_counties_df, 
              aes(x = time, ymin = spce_lower, ymax = spce_upper, 
                  fill = county), alpha = 0.2) +
  # Add SPCE lines
  geom_line(data = notable_counties_df, 
            aes(x = time, y = spce, color = county), size = 1) +
  # Add a reference line at zero
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  # Customize appearance
  labs(
    title = "County-Specific Survival Probability Causal Effect (SPCE) Over Time",
    subtitle = paste("WA", stage_val, "- Selected Counties"),
    x = "Time (years)",
    y = "SPCE (Treatment - Control)",
    color = "County",
    fill = "County"
  ) +
  scale_color_manual(
    values = county_colors,
    labels = paste("County", notable_counties)
  ) +
  scale_fill_manual(
    values = county_colors,
    labels = paste("County", notable_counties)
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

print(p_county_spce)

# Save county-specific SPCE plot
ggsave("county_specific_spce.png", p_county_spce, width = 10, height = 8)

# Plot Marginal SPCE
p_marginal_spce <- ggplot(marginal_spce, aes(x = time)) +
  # Add confidence interval ribbon
  geom_ribbon(aes(ymin = spce_lower, ymax = spce_upper), 
              fill = "steelblue", alpha = 0.3) +
  # Add SPCE line
  geom_line(aes(y = spce), color = "steelblue", size = 1.2) +
  # Add a reference line at zero
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  # Customize appearance
  labs(
    title = "Marginal Survival Probability Causal Effect (SPCE) Over Time",
    subtitle = paste("WA", stage_val, "- Average across all counties"),
    x = "Time (years)",
    y = "SPCE (Treatment - Control)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold")
  )

print(p_marginal_spce)

# Save marginal SPCE plot
ggsave("marginal_spce.png", p_marginal_spce, width = 10, height = 8)

# Create a combined plot with both county-specific and marginal SPCE
p_combined_spce <- p_county_spce +
  # Add marginal SPCE
  geom_ribbon(data = marginal_spce, 
              aes(x = time, ymin = spce_lower, ymax = spce_upper),
              fill = "black", alpha = 0.2) +
  geom_line(data = marginal_spce, 
            aes(x = time, y = spce),
            color = "black", size = 1.5) +
  # Update legend to include marginal SPCE
  scale_color_manual(
    values = c(county_colors, "black"),
    labels = c(paste("County", notable_counties), "Marginal SPCE")
  ) +
  scale_fill_manual(
    values = c(county_colors, "black"),
    labels = c(paste("County", notable_counties), "Marginal SPCE")
  )

print(p_combined_spce)

# Save combined SPCE plot
ggsave("combined_spce.png", p_combined_spce, width = 12, height = 8)

# Extract CSPCE for specific counties with rural/urban classification

# First, classify Florida counties as rural or urban
classify_florida_counties <- function() {
  # Based on 2020 Census and Florida Department of Health criteria
  # This is a specific list of Florida rural counties according to the Florida Department of Health
  rural_counties <- c(
    "Baker", "Bradford", "Calhoun", "Columbia", "DeSoto", "Dixie", "Franklin", 
    "Gadsden", "Gilchrist", "Glades", "Gulf", "Hamilton", "Hardee", "Hendry", 
    "Holmes", "Jackson", "Jefferson", "Lafayette", "Levy", "Liberty", "Madison", 
    "Monroe", "Okeechobee", "Putnam", "Suwannee", "Taylor", "Union", "Wakulla", 
    "Walton", "Washington"
  )
  
  # All 67 Florida counties
  all_counties <- c(
    "Alachua", "Baker", "Bay", "Bradford", "Brevard", "Broward", "Calhoun", 
    "Charlotte", "Citrus", "Clay", "Collier", "Columbia", "DeSoto", "Dixie", 
    "Duval", "Escambia", "Flagler", "Franklin", "Gadsden", "Gilchrist", 
    "Glades", "Gulf", "Hamilton", "Hardee", "Hendry", "Hernando", "Highlands", 
    "Hillsborough", "Holmes", "Indian River", "Jackson", "Jefferson", 
    "Lafayette", "Lake", "Lee", "Leon", "Levy", "Liberty", "Madison", 
    "Manatee", "Marion", "Martin", "Miami-Dade", "Monroe", "Nassau", 
    "Okaloosa", "Okeechobee", "Orange", "Osceola", "Palm Beach", "Pasco", 
    "Pinellas", "Polk", "Putnam", "St. Johns", "St. Lucie", "Santa Rosa", 
    "Sarasota", "Seminole", "Sumter", "Suwannee", "Taylor", "Union", 
    "Volusia", "Wakulla", "Walton", "Washington"
  )
  
  # Create dataframe with classification
  county_classification <- data.frame(
    county_name = all_counties,
    county_id = 1:67,
    rural_urban = ifelse(all_counties %in% rural_counties, "Rural", "Urban")
  )
  
  return(county_classification)
}

# Function to classify counties based on their numeric IDs
classify_counties_by_id <- function(county_ids) {
  # Get the full classification
  florida_counties <- classify_florida_counties()
  
  # Map county IDs to classification
  county_classification <- data.frame(
    county = county_ids,
    rural_urban = "Unknown"
  )
  
  # Match county IDs to the classification
  for(i in 1:length(county_ids)) {
    idx <- which(florida_counties$county_id == county_ids[i])
    if(length(idx) > 0) {
      county_classification$rural_urban[i] <- florida_counties$rural_urban[idx]
    }
  }
  
  return(county_classification)
}

# Classify the counties
county_class <- classify_counties_by_id(stratum_counties)

# Convert county to character for joining
county_class$county <- as.character(county_class$county)

# Join classification to CSPCE data
county_spce_with_class <- merge(
  county_spce_df,
  county_class,
  by = "county",
  all.x = TRUE
)

# Fill in any missing classifications
county_spce_with_class$rural_urban[is.na(county_spce_with_class$rural_urban)] <- "Unknown"

# Print summary of classified counties
rural_count <- length(unique(county_spce_with_class$county[county_spce_with_class$rural_urban == "Rural"]))
urban_count <- length(unique(county_spce_with_class$county[county_spce_with_class$rural_urban == "Urban"]))
unknown_count <- length(unique(county_spce_with_class$county[county_spce_with_class$rural_urban == "Unknown"]))

cat("Classification summary:\n")
cat("Rural counties:", rural_count, "\n")
cat("Urban counties:", urban_count, "\n")
if (unknown_count > 0) cat("Unknown counties:", unknown_count, "\n")

# Create the plot of all counties
p_all_counties <- ggplot(county_spce_with_class, 
                         aes(x = time, y = spce, 
                             group = county, 
                             color = rural_urban)) +
  # Add individual county lines
  geom_line(alpha = 0.7, size = 0.5) +
  # Add a reference line at zero
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  # Customize appearance
  labs(
    title = "Conditional Survival Probability Causal Effect (CSPCE) Over Time",
    subtitle = "Florida Counties by Rural/Urban Classification",
    x = "Time (years)",
    y = "CSPCE (Treatment - Control)",
    color = "County Type"
  ) +
  scale_color_manual(
    values = c("Rural" = "green", "Urban" = "red", "Unknown" = "gray")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

print(p_all_counties)

# Calculate average CSPCE by time and rural/urban group
avg_cspce <- county_spce_with_class %>%
  filter(rural_urban != "Unknown") %>%
  group_by(time, rural_urban) %>%
  summarize(
    spce_mean = mean(spce, na.rm = TRUE),
    spce_lower = mean(spce_lower, na.rm = TRUE),
    spce_upper = mean(spce_upper, na.rm = TRUE),
    .groups = "drop"
  )

# Create the plot of average CSPCE
p_avg_cspce <- ggplot(avg_cspce, aes(x = time, y = spce_mean, color = rural_urban)) +
  # Add confidence intervals
  geom_ribbon(aes(ymin = spce_lower, ymax = spce_upper, fill = rural_urban), 
              alpha = 0.2, color = NA) +
  # Add mean lines
  geom_line(size = 1.2) +
  # Add a reference line at zero
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  # Customize appearance
  labs(
    title = "Average Conditional Survival Probability Causal Effect (CSPCE)",
    subtitle = "By Rural/Urban Classification",
    x = "Time (years)",
    y = "Average CSPCE (Treatment - Control)",
    color = "County Type",
    fill = "County Type"
  ) +
  scale_color_manual(
    values = c("Rural" = "green", "Urban" = "red")
  ) +
  scale_fill_manual(
    values = c("Rural" = "green", "Urban" = "red")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

print(p_avg_cspce)

# Create a combined plot
combined_plot <- grid.arrange(
  p_all_counties, p_avg_cspce, 
  ncol = 1,
  heights = c(1, 1),
  top = textGrob("CSPCE by Rural/Urban Classification in Florida Counties", 
                 gp = gpar(fontsize = 16, font = 2))
)

print(combined_plot)

# Save the plots
ggsave("rural_urban_all_counties_cspce.png", p_all_counties, width = 12, height = 8)
ggsave("rural_urban_avg_cspce.png", p_avg_cspce, width = 12, height = 8)
ggsave("rural_urban_combined_cspce.png", combined_plot, width = 12, height = 14)

# Analyze the results
summary_stats <- avg_cspce %>%
  group_by(rural_urban) %>%
  summarize(
    max_cspce = max(spce_mean, na.rm = TRUE),
    min_cspce = min(spce_mean, na.rm = TRUE),
    time_at_max = time[which.max(spce_mean)],
    time_at_min = time[which.min(spce_mean)]
  )
print(summary_stats)

# Define counties to extract for specific examples
# Urban, large: Miami-Dade (county 43)
# Urban, small: Sarasota (county 58)
# Rural, large: Marion (county 42)
# Rural, small: Franklin (county 18)

counties_to_plot <- c("43", "58", "42", "18")
county_labels <- c("Miami-Dade (Urban Large)", "Sarasota (Urban Small)", 
                   "Marion (Rural Large)", "Franklin (Rural Small)")

# Create a subset of the CSPCE data for these counties
selected_counties_df <- county_spce_df %>%
  filter(county %in% counties_to_plot)

# Create a color palette for counties
county_colors <- c("red", "orange", "blue", "green")

# Plot county-specific CSPCE
p_selected_counties <- ggplot() +
  # Add confidence intervals for each county
  geom_ribbon(data = selected_counties_df, 
              aes(x = time, ymin = spce_lower, ymax = spce_upper, 
                  fill = county), alpha = 0.2) +
  # Add SPCE lines
  geom_line(data = selected_counties_df, 
            aes(x = time, y = spce, color = county), size = 1) +
  # Add a reference line at zero
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  # Customize appearance
  labs(
    title = "County-Specific Survival Probability Causal Effect (SPCE)",
    subtitle = "Comparison of Rural vs Urban Counties",
    x = "Time (years)",
    y = "SPCE (Treatment - Control)",
    color = "County",
    fill = "County"
  ) +
  scale_color_manual(
    values = county_colors,
    labels = county_labels
  ) +
  scale_fill_manual(
    values = county_colors,
    labels = county_labels
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    legend.box = "vertical"
  )

# View the plot
print(p_selected_counties)

# Save the plot
ggsave("rural_urban_counties_spce.png", p_selected_counties, width = 10, height = 8)

cat("CSPCE plot for selected rural and urban counties created and saved.\n")
cat("Analysis complete!\n")










# Age-specific Causal Effects Analysis with Different Characteristics
# ----------------------------------------------------------------------

# This code creates three plots of age-specific CACE (averaged over all patients of a specific age)
# Colored by: 1) Tumor Grade, 2) Biopsy Delay, and 3) HR status

# Load required packages (assuming these are already loaded from previous code)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(viridis)

#' Calculate age-specific causal effects with confidence intervals, grouped by clinical characteristics
#'
#' @param fit The fitted AFT_mixed_CAR_softbart_causal model
#' @param X Patient covariate matrix
#' @param Z Treatment vector
#' @param group County group indices
#' @param original_data Original data frame with unscaled covariates
#'
#' @return A data frame with age-specific causal effects by clinical characteristics
#'
calculate_age_specific_effects <- function(fit, X, Z, group, original_data) {
  # Helper function for logistic transformation
  expit <- function(x) 1/(1+exp(-x))
  
  # Extract number of MCMC samples
  n_samples <- nrow(fit$r)
  n_patients <- nrow(X)
  
  # Get unique ages
  unique_ages <- sort(unique(original_data$Age))
  
  # Initialize data frames to store results
  # For tumor grade
  results_tgrade <- data.frame()
  
  # For biopsy delay
  results_bxdelay <- data.frame()
  
  # For HR status
  results_hr <- data.frame()
  
  # Define functions to get category for each characteristic
  get_tgrade_category <- function(tgrade) {
    return(factor(tgrade, levels = c(1, 2, 3), labels = c("Grade 1", "Grade 2", "Grade 3")))
  }
  
  get_bxdelay_category <- function(bx_delay) {
    median_delay <- median(original_data$BX_Delay, na.rm = TRUE)
    return(factor(bx_delay > median_delay, levels = c(FALSE, TRUE), 
                  labels = c("Low Delay", "High Delay")))
  }
  
  get_hr_category <- function(hr) {
    # If HR is binary 0/1
    if(all(original_data$HR_p %in% c(0, 1))) {
      return(factor(hr, levels = c(0, 1), labels = c("HR Negative", "HR Positive")))
    } else {
      # If HR is continuous, use a threshold
      threshold <- 0.1  # Commonly used threshold
      return(factor(hr > threshold, levels = c(FALSE, TRUE), 
                    labels = c("HR Negative", "HR Positive")))
    }
  }
  
  # Process each age
  cat("Processing age-specific effects...\n")
  for(age in unique_ages) {
    cat("  Processing age:", age, "\n")
    
    # Get patients of this age
    age_indices <- which(original_data$Age == age)
    
    if(length(age_indices) == 0) {
      next
    }
    
    # Get unique tumor grades for this age
    unique_tgrades <- sort(unique(original_data$Tgrade[age_indices]))
    
    # Process each tumor grade
    for(tgrade in unique_tgrades) {
      # Get patients of this age and tumor grade
      patient_indices <- intersect(age_indices, which(original_data$Tgrade == tgrade))
      
      if(length(patient_indices) == 0) {
        next
      }
      
      # Initialize matrix to store effect samples for these patients
      patient_effect_samples <- matrix(NA, nrow = n_samples, ncol = length(patient_indices))
      
      # For each MCMC sample
      for(m in 1:n_samples) {
        # Extract parameters for this sample
        r_m <- fit$r[m,]
        sigma2_m <- fit$sigma2[m]
        
        # For each patient in this group
        for(j in 1:length(patient_indices)) {
          i <- patient_indices[j]
          
          # Get county for this patient
          county_idx <- group[i]
          
          # Get spatial effect for this county
          r_i_m <- r_m[county_idx]
          
          # Extract patient covariates
          X_i <- X[i,, drop = FALSE]
          
          # Get estimated propensity score
          ps_raw <- fit$propensity_forest$do_predict(X_i)
          e_hat_i <- expit(ps_raw)
          
          # Create augmented matrices for treatment and control
          X_aug1 <- cbind(X_i, 1, e_hat_i)  # Treatment
          X_aug0 <- cbind(X_i, 0, e_hat_i)  # Control
          
          # Get function predictions
          f1 <- fit$forest$do_predict(X_aug1)
          f0 <- fit$forest$do_predict(X_aug0)
          
          # Calculate potential outcomes
          # E[T] = exp(mu + sigma^2/2) for lognormal
          mu1 <- f1 + r_i_m
          mu0 <- f0 + r_i_m
          
          E_T1 <- exp(mu1 + 0.5*sigma2_m)
          E_T0 <- exp(mu0 + 0.5*sigma2_m)
          
          # Calculate individual treatment effect
          patient_effect_samples[m, j] <- E_T1 - E_T0
        }
      }
      
      # Calculate mean effect across MCMC samples for each patient
      mean_effect_by_patient <- colMeans(patient_effect_samples)
      
      # Calculate mean effect and CI across all patients in this group
      mean_effect <- mean(mean_effect_by_patient)
      
      # For CI, use the mean effect across patients for each MCMC sample
      mean_effect_by_sample <- rowMeans(patient_effect_samples)
      ci_lower <- quantile(mean_effect_by_sample, 0.025)
      ci_upper <- quantile(mean_effect_by_sample, 0.975)
      
      # Store results
      tgrade_result <- data.frame(
        Age = age,
        Category = get_tgrade_category(tgrade),
        CharacteristicType = "Tumor Grade",
        CACE = mean_effect,
        CACE_lower = ci_lower,
        CACE_upper = ci_upper,
        n_patients = length(patient_indices)
      )
      
      results_tgrade <- rbind(results_tgrade, tgrade_result)
    }
    
    # Process biopsy delay (split at median)
    median_delay <- median(original_data$BX_Delay, na.rm = TRUE)
    delay_categories <- c(FALSE, TRUE)  # Low delay, High delay
    
    for(delay_cat in delay_categories) {
      if(delay_cat) {
        # High delay
        patient_indices <- intersect(age_indices, which(original_data$BX_Delay > median_delay))
      } else {
        # Low delay
        patient_indices <- intersect(age_indices, which(original_data$BX_Delay <= median_delay))
      }
      
      if(length(patient_indices) == 0) {
        next
      }
      
      # Initialize matrix to store effect samples for these patients
      patient_effect_samples <- matrix(NA, nrow = n_samples, ncol = length(patient_indices))
      
      # For each MCMC sample
      for(m in 1:n_samples) {
        # Extract parameters for this sample
        r_m <- fit$r[m,]
        sigma2_m <- fit$sigma2[m]
        
        # For each patient in this group
        for(j in 1:length(patient_indices)) {
          i <- patient_indices[j]
          
          # Get county for this patient
          county_idx <- group[i]
          
          # Get spatial effect for this county
          r_i_m <- r_m[county_idx]
          
          # Extract patient covariates
          X_i <- X[i,, drop = FALSE]
          
          # Get estimated propensity score
          ps_raw <- fit$propensity_forest$do_predict(X_i)
          e_hat_i <- expit(ps_raw)
          
          # Create augmented matrices for treatment and control
          X_aug1 <- cbind(X_i, 1, e_hat_i)  # Treatment
          X_aug0 <- cbind(X_i, 0, e_hat_i)  # Control
          
          # Get function predictions
          f1 <- fit$forest$do_predict(X_aug1)
          f0 <- fit$forest$do_predict(X_aug0)
          
          # Calculate potential outcomes
          # E[T] = exp(mu + sigma^2/2) for lognormal
          mu1 <- f1 + r_i_m
          mu0 <- f0 + r_i_m
          
          E_T1 <- exp(mu1 + 0.5*sigma2_m)
          E_T0 <- exp(mu0 + 0.5*sigma2_m)
          
          # Calculate individual treatment effect
          patient_effect_samples[m, j] <- E_T1 - E_T0
        }
      }
      
      # Calculate mean effect across MCMC samples for each patient
      mean_effect_by_patient <- colMeans(patient_effect_samples)
      
      # Calculate mean effect and CI across all patients in this group
      mean_effect <- mean(mean_effect_by_patient)
      
      # For CI, use the mean effect across patients for each MCMC sample
      mean_effect_by_sample <- rowMeans(patient_effect_samples)
      ci_lower <- quantile(mean_effect_by_sample, 0.025)
      ci_upper <- quantile(mean_effect_by_sample, 0.975)
      
      # Store results
      delay_result <- data.frame(
        Age = age,
        Category = get_bxdelay_category(delay_cat),
        CharacteristicType = "Biopsy Delay",
        CACE = mean_effect,
        CACE_lower = ci_lower,
        CACE_upper = ci_upper,
        n_patients = length(patient_indices)
      )
      
      results_bxdelay <- rbind(results_bxdelay, delay_result)
    }
    
    # Process HR status
    if(all(original_data$HR_p %in% c(0, 1))) {
      # Binary HR status
      hr_categories <- c(0, 1)
      
      for(hr_cat in hr_categories) {
        patient_indices <- intersect(age_indices, which(original_data$HR_p == hr_cat))
        
        if(length(patient_indices) == 0) {
          next
        }
        
        # Initialize matrix to store effect samples for these patients
        patient_effect_samples <- matrix(NA, nrow = n_samples, ncol = length(patient_indices))
        
        # For each MCMC sample
        for(m in 1:n_samples) {
          # Extract parameters for this sample
          r_m <- fit$r[m,]
          sigma2_m <- fit$sigma2[m]
          
          # For each patient in this group
          for(j in 1:length(patient_indices)) {
            i <- patient_indices[j]
            
            # Get county for this patient
            county_idx <- group[i]
            
            # Get spatial effect for this county
            r_i_m <- r_m[county_idx]
            
            # Extract patient covariates
            X_i <- X[i,, drop = FALSE]
            
            # Get estimated propensity score
            ps_raw <- fit$propensity_forest$do_predict(X_i)
            e_hat_i <- expit(ps_raw)
            
            # Create augmented matrices for treatment and control
            X_aug1 <- cbind(X_i, 1, e_hat_i)  # Treatment
            X_aug0 <- cbind(X_i, 0, e_hat_i)  # Control
            
            # Get function predictions
            f1 <- fit$forest$do_predict(X_aug1)
            f0 <- fit$forest$do_predict(X_aug0)
            
            # Calculate potential outcomes
            # E[T] = exp(mu + sigma^2/2) for lognormal
            mu1 <- f1 + r_i_m
            mu0 <- f0 + r_i_m
            
            E_T1 <- exp(mu1 + 0.5*sigma2_m)
            E_T0 <- exp(mu0 + 0.5*sigma2_m)
            
            # Calculate individual treatment effect
            patient_effect_samples[m, j] <- E_T1 - E_T0
          }
        }
        
        # Calculate mean effect across MCMC samples for each patient
        mean_effect_by_patient <- colMeans(patient_effect_samples)
        
        # Calculate mean effect and CI across all patients in this group
        mean_effect <- mean(mean_effect_by_patient)
        
        # For CI, use the mean effect across patients for each MCMC sample
        mean_effect_by_sample <- rowMeans(patient_effect_samples)
        ci_lower <- quantile(mean_effect_by_sample, 0.025)
        ci_upper <- quantile(mean_effect_by_sample, 0.975)
        
        # Store results
        hr_result <- data.frame(
          Age = age,
          Category = get_hr_category(hr_cat),
          CharacteristicType = "HR Status",
          CACE = mean_effect,
          CACE_lower = ci_lower,
          CACE_upper = ci_upper,
          n_patients = length(patient_indices)
        )
        
        results_hr <- rbind(results_hr, hr_result)
      }
    } else {
      # Continuous HR status (use threshold)
      threshold <- 0.1
      hr_categories <- c(FALSE, TRUE)  # Below threshold, Above threshold
      
      for(hr_cat in hr_categories) {
        if(hr_cat) {
          # HR positive
          patient_indices <- intersect(age_indices, which(original_data$HR_p > threshold))
        } else {
          # HR negative
          patient_indices <- intersect(age_indices, which(original_data$HR_p <= threshold))
        }
        
        if(length(patient_indices) == 0) {
          next
        }
        
        # Initialize matrix to store effect samples for these patients
        patient_effect_samples <- matrix(NA, nrow = n_samples, ncol = length(patient_indices))
        
        # For each MCMC sample
        for(m in 1:n_samples) {
          # Extract parameters for this sample
          r_m <- fit$r[m,]
          sigma2_m <- fit$sigma2[m]
          
          # For each patient in this group
          for(j in 1:length(patient_indices)) {
            i <- patient_indices[j]
            
            # Get county for this patient
            county_idx <- group[i]
            
            # Get spatial effect for this county
            r_i_m <- r_m[county_idx]
            
            # Extract patient covariates
            X_i <- X[i,, drop = FALSE]
            
            # Get estimated propensity score
            ps_raw <- fit$propensity_forest$do_predict(X_i)
            e_hat_i <- expit(ps_raw)
            
            # Create augmented matrices for treatment and control
            X_aug1 <- cbind(X_i, 1, e_hat_i)  # Treatment
            X_aug0 <- cbind(X_i, 0, e_hat_i)  # Control
            
            # Get function predictions
            f1 <- fit$forest$do_predict(X_aug1)
            f0 <- fit$forest$do_predict(X_aug0)
            
            # Calculate potential outcomes
            # E[T] = exp(mu + sigma^2/2) for lognormal
            mu1 <- f1 + r_i_m
            mu0 <- f0 + r_i_m
            
            E_T1 <- exp(mu1 + 0.5*sigma2_m)
            E_T0 <- exp(mu0 + 0.5*sigma2_m)
            
            # Calculate individual treatment effect
            patient_effect_samples[m, j] <- E_T1 - E_T0
          }
        }
        
        # Calculate mean effect across MCMC samples for each patient
        mean_effect_by_patient <- colMeans(patient_effect_samples)
        
        # Calculate mean effect and CI across all patients in this group
        mean_effect <- mean(mean_effect_by_patient)
        
        # For CI, use the mean effect across patients for each MCMC sample
        mean_effect_by_sample <- rowMeans(patient_effect_samples)
        ci_lower <- quantile(mean_effect_by_sample, 0.025)
        ci_upper <- quantile(mean_effect_by_sample, 0.975)
        
        # Store results
        hr_result <- data.frame(
          Age = age,
          Category = get_hr_category(hr_cat),
          CharacteristicType = "HR Status",
          CACE = mean_effect,
          CACE_lower = ci_lower,
          CACE_upper = ci_upper,
          n_patients = length(patient_indices)
        )
        
        results_hr <- rbind(results_hr, hr_result)
      }
    }
  }
  
  # Return all results
  return(list(
    tgrade = results_tgrade,
    bxdelay = results_bxdelay,
    hr = results_hr
  ))
}

# Calculate age-specific effects
age_effects <- calculate_age_specific_effects(
  fit = fit,
  X = X,
  Z = Z,
  group = group,
  original_data = data_stratum
)

# Create a simplified function that doesn't use geom_ribbon at all
create_cace_age_plot <- function(data, characteristic_name, color_values) {
  # Create a basic plot without facets or ribbons
  p <- ggplot(data, aes(x = Age, y = CACE)) +
    # Add points with color and size
    geom_point(aes(color = Category, size = n_patients), alpha = 0.7) +
    # Add lines with color
    geom_line(aes(color = Category, group = Category), alpha = 0.7) +
    # Add error bars with color
    geom_errorbar(aes(ymin = CACE_lower, ymax = CACE_upper, color = Category, group = Category), 
                  width = 0.5, alpha = 0.5) +
    # Add a reference line at zero
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1) +
    # Apply color scale
    scale_color_manual(values = color_values) +
    # Scale size
    scale_size_continuous(name = "Number of\nPatients", range = c(1, 5)) +
    # Customize appearance
    labs(
      title = paste("Age-Specific Conditional Average Causal Effect"),
      subtitle = paste("By", characteristic_name),
      x = "Age (years)",
      y = "CACE (days)",
      color = characteristic_name
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      legend.position = "bottom",
      legend.box = "vertical"
    )
  
  return(p)
}

# Create the three plots using the new function
# 1. Tumor Grade (already works)
plot_cace_by_tgrade <- create_cace_age_plot(
  data = age_effects$tgrade,
  characteristic_name = "Tumor Grade",
  color_values = c("Grade 1" = viridis::viridis(3, option = "plasma", end = 0.8)[1],
                   "Grade 2" = viridis::viridis(3, option = "plasma", end = 0.8)[2],
                   "Grade 3" = viridis::viridis(3, option = "plasma", end = 0.8)[3])
)

# 2. Biopsy Delay - fix by ensuring color values match the exact category labels
plot_cace_by_bxdelay <- create_cace_age_plot(
  data = age_effects$bxdelay,
  characteristic_name = "Biopsy Delay",
  color_values = c("Low Delay" = "blue", "High Delay" = "red")
)

# 3. HR Status - fix by ensuring color values match the exact category labels
plot_cace_by_hr <- create_cace_age_plot(
  data = age_effects$hr,
  characteristic_name = "HR Status",
  color_values = c("HR Negative" = "purple", "HR Positive" = "green")
)

# 2. Create plot of CACE vs Age colored by Biopsy Delay
plot_cace_by_bxdelay <- ggplot(age_effects$bxdelay, 
                               aes(x = Age, y = CACE, color = Category, 
                                   group = Category,
                                   size = n_patients)) +
  # Add confidence intervals
  geom_ribbon(aes(ymin = CACE_lower, ymax = CACE_upper, fill = Category, group = Category), 
              alpha = 0.2, color = NA) +
  # Add points
  geom_point(alpha = 0.7) +
  # Add lines
  geom_line(alpha = 0.7) +
  # Add a reference line at zero
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1) +
  # Use contrasting colors for two categories
  scale_color_manual(values = c("Low Delay" = "blue", "High Delay" = "red")) +
  scale_fill_manual(values = c("Low Delay" = "blue", "High Delay" = "red")) +
  # Scale size
  scale_size_continuous(name = "Number of\nPatients", range = c(1, 5)) +
  # Customize appearance
  labs(
    title = "Age-Specific Conditional Average Causal Effect",
    subtitle = "By Biopsy Delay",
    x = "Age (years)",
    y = "CACE (days)",
    color = "Biopsy Delay",
    fill = "Biopsy Delay"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom",
    legend.box = "vertical"
  )

# 3. Create plot of CACE vs Age colored by HR status
plot_cace_by_hr <- ggplot(age_effects$hr, 
                          aes(x = Age, y = CACE, color = Category, 
                              group = Category,
                              size = n_patients)) +
  # Add confidence intervals
  geom_ribbon(aes(ymin = CACE_lower, ymax = CACE_upper, fill = Category, group = Category), 
              alpha = 0.2, color = NA) +
  # Add points
  geom_point(alpha = 0.7) +
  # Add lines
  geom_line(alpha = 0.7) +
  # Add a reference line at zero
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1) +
  # Use contrasting colors for two categories
  scale_color_manual(values = c("HR Negative" = "purple", "HR Positive" = "green")) +
  scale_fill_manual(values = c("HR Negative" = "purple", "HR Positive" = "green")) +
  # Scale size
  scale_size_continuous(name = "Number of\nPatients", range = c(1, 5)) +
  # Customize appearance
  labs(
    title = "Age-Specific Conditional Average Causal Effect",
    subtitle = "By Hormone Receptor Status",
    x = "Age (years)",
    y = "CACE (days)",
    color = "HR Status",
    fill = "HR Status"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom",
    legend.box = "vertical"
  )

# Display the plots
print(plot_cace_by_tgrade)
print(plot_cace_by_bxdelay)
print(plot_cace_by_hr)

# Save the plots
ggsave("cace_by_age_tumor_grade.png", plot_cace_by_tgrade, width = 10, height = 8)
ggsave("cace_by_age_biopsy_delay.png", plot_cace_by_bxdelay, width = 10, height = 8)
ggsave("cace_by_age_hr_status.png", plot_cace_by_hr, width = 10, height = 8)

# Create a combined plot
combined_plots <- grid.arrange(
  plot_cace_by_tgrade, 
  plot_cace_by_bxdelay, 
  plot_cace_by_hr, 
  ncol = 1,
  heights = c(1, 1, 1)
)

# Save the combined plot
ggsave("combined_cace_by_age_plots.png", combined_plots, width = 12, height = 24)

cat("Age-specific CACE analysis complete!\n")



# Create ggplot2 plots for age-specific CACE
# ----------------------------------------------------------------------

library(ggplot2)
library(dplyr)

# Creating a combined plot with all categories in one plot
create_combined_cace_plot <- function(data, characteristic_name, colors, category_labels) {
  # Create the plot with all categories together
  p <- ggplot(data, aes(x = Age, y = CACE, color = Category, fill = Category, group = Category)) +
    # Add confidence interval as ribbons
    geom_ribbon(aes(ymin = CACE_lower, ymax = CACE_upper), alpha = 0.2, color = NA) +
    # Add lines
    geom_line(size = 1) +
    # Add points without size variation
    geom_point(size = 3, alpha = 0.7) +
    # Add horizontal line at zero
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1) +
    # Set colors
    scale_color_manual(values = colors, labels = category_labels) +
    scale_fill_manual(values = colors, labels = category_labels) +
    # Set title and labels
    labs(
      title = paste("Age-Specific Conditional Average Causal Effect by", characteristic_name),
      x = "Age (years)",
      y = "CACE (days)",
      color = characteristic_name,
      fill = characteristic_name
    ) +
    # Set theme
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      axis.title = element_text(face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold")
    )
  
  return(p)
}

# Create separate function for binary variables to ensure proper coding
create_binary_cace_plot <- function(data, characteristic_name, color_0, color_1, label_0, label_1) {
  # Ensure category is properly coded as factor with correct labels
  data$Category <- factor(data$Category, levels = unique(data$Category))
  
  # Create plot
  p <- ggplot(data, aes(x = Age, y = CACE, color = Category, fill = Category, group = Category)) +
    # Add confidence interval as ribbons
    geom_ribbon(aes(ymin = CACE_lower, ymax = CACE_upper), alpha = 0.2, color = NA) +
    # Add lines
    geom_line(size = 1) +
    # Add points without size variation
    geom_point(size = 3, alpha = 0.7) +
    # Add horizontal line at zero
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1) +
    # Set colors for binary variable (0 and 1)
    scale_color_manual(values = c(color_0, color_1), labels = c(label_0, label_1)) +
    scale_fill_manual(values = c(color_0, color_1), labels = c(label_0, label_1)) +
    # Set title and labels
    labs(
      title = paste("Age-Specific Conditional Average Causal Effect by", characteristic_name),
      x = "Age (years)",
      y = "CACE (days)",
      color = characteristic_name,
      fill = characteristic_name
    ) +
    # Set theme
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      axis.title = element_text(face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold")
    )
  
  return(p)
}

# 1. Tumor Grade plot
# Define colors and labels for tumor grades
tgrade_colors <- c("darkblue", "darkgreen", "darkred")
tgrade_labels <- c("Grade 1", "Grade 2", "Grade 3")

# Create the plot
plot_tgrade <- create_combined_cace_plot(
  data = age_effects$tgrade,
  characteristic_name = "",
  colors = tgrade_colors,
  category_labels = tgrade_labels
)

# 2. Biopsy Delay plot
# Not binary but categorized as 1 and 3
plot_bxdelay <- create_binary_cace_plot(
  data = age_effects$bxdelay,
  characteristic_name = "Biopsy Delay",
  color_0 = "darkblue",  # Color for BD=1 (Short Delay)
  color_1 = "darkred",   # Color for BD=3 (Long Delay)
  label_0 = "Short Delay (BD=1)",
  label_1 = "Long Delay (BD=3)"
)

# 3. HR Status plot
# Binary variable with 0/1 coding
plot_hr <- create_binary_cace_plot(
  data = age_effects$hr,
  characteristic_name = "",
  color_0 = "purple",    # Color for HR=0 (Negative)
  color_1 = "darkgreen", # Color for HR=1 (Positive)
  label_0 = "HR Negative ",
  label_1 = "HR Positive "
)

# Save the plots
ggsave("cace_by_age_tumor_grade.png", plot_tgrade, width = 10, height = 7)
ggsave("cace_by_age_biopsy_delay.png", plot_bxdelay, width = 10, height = 7)
ggsave("cace_by_age_hr_status.png", plot_hr, width = 10, height = 7)

# Create a combined plot with tumor grade and HR status together
library(gridExtra)
library(grid)

# Modify the tumor grade plot to have a more concise title
plot_tgrade_concise <- plot_tgrade + 
  labs(
    title = "Tumor Grade",
    x = "Age (years)",
    y = "CACE (days)"
  )

# Modify the HR status plot to have a more concise title
plot_hr_concise <- plot_hr + 
  labs(
    title = "HR Status",
    x = "Age (years)",
    y = "CACE (days)"
  )

# Combine the plots
combined_plot <- grid.arrange(
  plot_tgrade_concise, plot_hr_concise,
  ncol = 2,
  top = textGrob(
    "Age-Specfic Effects",
    gp = gpar(fontface = "bold", fontsize = 10)
  )
)

# Save the combined plot
ggsave("combined_tgrade_hr_plot.png", combined_plot, width = 14, height = 7)

# Print the combined plot
print(combined_plot)

cat("Combined plot of tumor grade and HR status created successfully!\n")




# Plot county-specific CRACE (Restricted Average Causal Effect) for 5 years vs county size

# Create modified plot function for rural/urban coloring
plot_county_crace_rural_urban <- function(county_crace_df, 
                                          title = "County-Specific 5-Year CRACE") {
  require(ggplot2)
  
  # Create a copy of the data
  plot_data <- county_crace_df
  
  # Remove counties with missing values or zero sample size
  plot_data <- plot_data[!is.na(plot_data$crace) & plot_data$sample_size > 0, ]
  
  # Get overall RACE (mean of county-specific CRACE values)
  race <- mean(plot_data$crace, na.rm = TRUE)
  race_lower <- mean(plot_data$crace_lower, na.rm = TRUE)
  race_upper <- mean(plot_data$crace_upper, na.rm = TRUE)
  
  # Calculate quartiles for rural and urban counties separately
  rural_data <- plot_data[plot_data$rural_urban == "Rural", ]
  urban_data <- plot_data[plot_data$rural_urban == "Urban", ]
  
  # Rural quartiles
  rural_q1 <- quantile(rural_data$crace, 0.25, na.rm = TRUE)
  rural_q2 <- quantile(rural_data$crace, 0.50, na.rm = TRUE) # median
  rural_q3 <- quantile(rural_data$crace, 0.75, na.rm = TRUE)
  
  # Urban quartiles
  urban_q1 <- quantile(urban_data$crace, 0.25, na.rm = TRUE)
  urban_q2 <- quantile(urban_data$crace, 0.50, na.rm = TRUE) # median
  urban_q3 <- quantile(urban_data$crace, 0.75, na.rm = TRUE)
  
  # Create base plot with rural/urban coloring
  p <- ggplot(plot_data, aes(x = sample_size, y = crace, color = rural_urban)) +
    # Zero line (bold black dashed)
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 1) +
    
    # Rural quartile lines
    geom_hline(yintercept = rural_q1, linetype = "dotted", color = "forestgreen", size = 0.7) +
    geom_hline(yintercept = rural_q2, linetype = "solid", color = "forestgreen", size = 0.7) +
    geom_hline(yintercept = rural_q3, linetype = "dotted", color = "forestgreen", size = 0.7) +
    
    # Urban quartile lines
    geom_hline(yintercept = urban_q1, linetype = "dotted", color = "purple", size = 0.7) +
    geom_hline(yintercept = urban_q2, linetype = "solid", color = "purple", size = 0.7) +
    geom_hline(yintercept = urban_q3, linetype = "dotted", color = "purple", size = 0.7) +
    
    # Overall RACE lines
    geom_hline(yintercept = race, linetype = "dashed", color = "red") +
    geom_hline(yintercept = race_lower, linetype = "dotted", color = "red") +
    geom_hline(yintercept = race_upper, linetype = "dotted", color = "red") +
    
    # Points and error bars
    geom_point(alpha = 0.8, size = 3) +
    geom_errorbar(aes(ymin = crace_lower, ymax = crace_upper), width = 0, alpha = 0.7) +
    
    # Labels and aesthetics
    labs(
      title = title,
      subtitle = paste0("Overall RACE: ", round(race, 2), " days (95% CI: ", 
                        round(race_lower, 2), " to ", round(race_upper, 2), ")"),
      x = "County Sample Size",
      y = "Days",
      color = "County Type"
    ) +
    scale_color_manual(
      values = c("Rural" = "forestgreen", "Urban" = "purple", "Unknown" = "gray50")
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      legend.position = "bottom"
    ) +
    
    # Add annotations for quartile lines
    annotate("text", x = max(plot_data$sample_size) * 0.95, y = rural_q1, 
             label = paste0("Rural Q1: ", round(rural_q1, 1)), 
             color = "forestgreen", hjust = 1, size = 3) +
    annotate("text", x = max(plot_data$sample_size) * 0.95, y = rural_q2, 
             label = paste0("Rural Q2: ", round(rural_q2, 1)), 
             color = "forestgreen", hjust = 1, size = 3) +
    annotate("text", x = max(plot_data$sample_size) * 0.95, y = rural_q3, 
             label = paste0("Rural Q3: ", round(rural_q3, 1)), 
             color = "forestgreen", hjust = 1, size = 3) +
    
    annotate("text", x = max(plot_data$sample_size) * 0.95, y = urban_q1, 
             label = paste0("Urban Q1: ", round(urban_q1, 1)), 
             color = "purple", hjust = 1, size = 3) +
    annotate("text", x = max(plot_data$sample_size) * 0.95, y = urban_q2, 
             label = paste0("Urban Q2: ", round(urban_q2, 1)), 
             color = "purple", hjust = 1, size = 3) +
    annotate("text", x = max(plot_data$sample_size) * 0.95, y = urban_q3, 
             label = paste0("Urban Q3: ", round(urban_q3, 1)), 
             color = "purple", hjust = 1, size = 3) +
    
    annotate("text", x = max(plot_data$sample_size) * 0.95, y = 0, 
             label = "Zero Effect", 
             color = "black", hjust = 1, size = 3)
  
  return(p)
}

# Function to calculate CRACE from CSPCE
calculate_crace <- function(time_array, spce_array, max_years = 5) {
  # Convert max_years to days (consistent with time scale in data)
  max_time <- max_years * 365.25
  
  # Filter points within our time range
  valid_indices <- which(time_array <= max_time)
  
  # If no valid points, return 0
  if(length(valid_indices) == 0) {
    return(0)
  }
  
  # Calculate the integral using the trapezoidal rule
  integral <- 0
  for(i in 1:(length(valid_indices)-1)) {
    idx1 <- valid_indices[i]
    idx2 <- valid_indices[i+1]
    
    time1 <- time_array[idx1]
    time2 <- time_array[idx2]
    spce1 <- spce_array[idx1]
    spce2 <- spce_array[idx2]
    
    # Trapezoidal rule: area = (width) * (average height)
    integral <- integral + (time2 - time1) * (spce1 + spce2) / 2
  }
  
  return(integral)
}

# Calculate CRACE for each county
county_crace_results <- list()
for(county in names(county_spce_results)) {
  county_data <- county_spce_results[[county]]
  
  crace <- calculate_crace(county_data$time, county_data$spce, max_years = 5)
  crace_lower <- calculate_crace(county_data$time, county_data$spce_lower, max_years = 5)
  crace_upper <- calculate_crace(county_data$time, county_data$spce_upper, max_years = 5)
  
  county_crace_results[[county]] <- list(
    crace = crace,
    crace_lower = crace_lower,
    crace_upper = crace_upper
  )
}

# Create a data frame for easier analysis and visualization
county_crace_df <- data.frame(
  county = numeric(length(county_crace_results)),
  crace = numeric(length(county_crace_results)),
  crace_lower = numeric(length(county_crace_results)),
  crace_upper = numeric(length(county_crace_results))
)

# Fill the data frame
for(i in 1:length(county_crace_results)) {
  county_id <- as.numeric(names(county_crace_results)[i])
  county_crace_df$county[i] <- county_id
  county_crace_df$crace[i] <- county_crace_results[[i]]$crace
  county_crace_df$crace_lower[i] <- county_crace_results[[i]]$crace_lower
  county_crace_df$crace_upper[i] <- county_crace_results[[i]]$crace_upper
}

# Add county sample size
county_crace_df <- merge(county_crace_df, data.frame(
  county = county_effects$county,
  sample_size = county_effects$sample_size
), by = "county", all.x = TRUE)

# Join with rural/urban classification
county_crace_df <- merge(county_crace_df, county_class, by = "county", all.x = TRUE)

# Create plot of average treatment effects vs. county size with rural/urban coloring
p_crace_by_size <- plot_county_crace_rural_urban(
  county_crace_df = county_crace_df,
  title = "County-Specific 5-Year CRACE"
)

# Display the plot
print(p_crace_by_size)

# Save the plot
ggsave("county_crace_by_size_plot.png", p_crace_by_size, width = 10, height = 8)