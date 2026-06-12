# This estimator puts a prior on the hazard difference, leading to a more stable estimate with higher right-censoring.
# Here we implement this using a Piecewise Exponential Model (PEM), assuming baseline hazard and hazard difference to 
# be constant in pre-defined intervals.

library(dplyr)
library(tidyr)
library(survival)
library(brms)
library(HDInterval)

estimateStablizedBayesianRiskDifference <- function(population, maxTime = 365) {
  population <- population |>
    mutate(id = row_number())
  
  # Convert to Piecewise Exponential Format (Counting Process)
  interval_width <- 30
  cut_points <- unique(c(seq(0, maxTime, by = interval_width), maxTime))
  
  # survSplit breaks each patient's follow-up into the defined intervals
  pem_dat <- survSplit(Surv(survivalTime, y) ~ a + id, 
                       data = population, 
                       cut = cut_points, 
                       episode = "interval_id")
  pem_dat <- pem_dat |>
    filter(tstart != maxTime)

  # Calculate exact exposure survivalTime within each interval and format factors
  pem_dat <- pem_dat |>
    mutate(
      exposure = survivalTime - tstart,
      log_exposure = log(exposure),
      # Create a factor for each survivalTime bin to estimate survivalTime-varying baseline hazards
      interval_fac = as.factor(interval_id),
      a = as.integer(a)
    )

  # Define the model formula
  # 0 + interval_fac estimates a distinct baseline hazard for every survivalTime bin
  # interval_fac:a estimates the treatment effect specifically within each survivalTime bin
  pem_formula <- bf(
    y ~ 0 + interval_fac + interval_fac:a + offset(log_exposure),
    family = poisson(link = "log")
  )
  
  pem_priors <- c(
    set_prior("normal(-5, 2)", class = "b") 
  )
  
  # Dynamically append the strong shrinkage prior for each treatment interaction term
  intervals <- levels(pem_dat$interval_fac)
  
  for (i in intervals) {
    # Construct the exact coefficient name as brms will generate it
    coef_name <- paste0("interval_fac", i, ":a")
    
    # Append to the prior object
    pem_priors <- c(
      pem_priors,
      set_prior("normal(0, 0.2)", class = "b", coef = coef_name)
    )
  }
  
  # (Optional) Verify the prior object looks correct
  # print(pem_priors)
  
  # Fit the model as before
  fit_pem <- brm(
    formula = pem_formula,
    data = pem_dat,
    prior = pem_priors,
    chains = 4,
    iter = 2000,
    cores = 4,
    seed = 123
  )
  
  # Generate results from samples
  
  # Dynamically extract the exact width and end survivalTime for every interval
  interval_info <- pem_dat |>
    group_by(interval_fac) |>
    summarize(
      t_start = min(tstart),
      t_end = max(survivalTime),
      .groups = "drop"
    ) |>
    mutate(width = t_end - t_start)
  
  interval_widths <- interval_info$width
  interval_ends <- interval_info$t_end
  
  # Extract posterior draws
  post_draws <- as_draws_matrix(fit_pem)
  base_cols <- grep("^b_interval_fac[0-9]+$", colnames(post_draws), value = TRUE)
  trt_cols  <- grep("^b_interval_fac[0-9]+:a$", colnames(post_draws), value = TRUE)
  
  n_intervals <- length(base_cols)
  n_draws <- nrow(post_draws)
  
  rd_matrix <- matrix(NA, nrow = n_draws, ncol = n_intervals)
  
  # 2. Calculate the cumulative metrics using the dynamic interval widths
  for (i in 1:n_draws) {
    
    log_haz_0 <- post_draws[i, base_cols]
    log_haz_1 <- log_haz_0 + post_draws[i, trt_cols]
    
    haz_0 <- exp(log_haz_0)
    haz_1 <- exp(log_haz_1)
    
    # Multiply by the exact width of each interval for the Riemann sum
    H_0 <- cumsum(haz_0 * interval_widths)
    H_1 <- cumsum(haz_1 * interval_widths)
    
    F_0 <- 1 - exp(-H_0)
    F_1 <- 1 - exp(-H_1)
    
    rd_matrix[i, ] <- F_1 - F_0
  }

  # 3. Build the summary data frame using the dynamic interval end times
  cis <- apply(rd_matrix, 2, hdi)
  rd_summary <- data.frame(
    timePoint = interval_ends, 
    estimate = apply(rd_matrix, 2, median),
    lb = cis[1, ],
    ub = cis[2, ]
  )
  return(rd_summary)
}