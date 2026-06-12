# This Bayesian risk difference estimator uses a prior on the hazards. (By default Jeffreys prior, beta distribution 
# with alpha = beta = 0.5)

library(survival)
library(HDInterval)

samplePosterior <- function(t, time0, event0, time1, event1, sampleSize = 10000, priorAlpha = 0.5, priorBeta = 0.5) {
  
  computeSurvivalPosterior <- function(times, events) {
    fit <- survfit(Surv(times, events) ~ 1)
    
    # Filter only to times up to t
    idx <- fit$time <= t
    nAtRisk <- fit$n.risk[idx]
    nEvent <- fit$n.event[idx]
    
    # If no events occurred before t, survival is always 1
    if (length(nEvent) == 0) {
      return(rep(1, sampleSize))
    }
    
    # Generate sampleSize draws for the hazard at each time point from Beta posteriors
    # Matrix dimensions: sampleSize x length(nEvent) columns (time points)
    posteriorHazards <- sapply(seq_along(nEvent), function(j) {
      rbeta(sampleSize, 
            shape1 = priorAlpha + nEvent[j], 
            shape2 = priorBeta + nAtRisk[j] - nEvent[j])
    })
    
    # Calculate S(t) = product of (1 - hazard) across all time points up to t
    # apply(..., 1, prod) multiplies across the columns for each row (simulation)
    sample <- apply(1 - posteriorHazards, 1, prod)
    
    return(sample)
  }
  
  posteriorS0 <- computeSurvivalPosterior(time0, event0)
  posteriorS1 <- computeSurvivalPosterior(time1, event1)
  posteriorRiskDifference <- posteriorS0 - posteriorS1 
  return(posteriorRiskDifference)
}

.estimateBayesianRiskDifference <- function(t, population) {
  target <- population |>
    filter(a == 1)
  comparator <- population |>
    filter(a == 0)
  sample <- samplePosterior(t = t, 
                            time0 = comparator$survivalTime,
                            event0 = comparator$y, 
                            time1 = target$survivalTime, 
                            event1 = target$y)
  ci <- hdi(sample)
  rd <- tibble(
    estimate = median(sample),
    lb = ci[1],
    ub = ci[2],
    se = sd(sample)
  )
  return(rd)
}

estimateBayesianRiskDifference <- function(population, 
                                           timePoints = c(180, 365, 730, 1095, 1460)) {
  estimates <- lapply(timePoints, .estimateBayesianRiskDifference, population = population)
  estimates <- bind_rows(estimates) |>
    mutate(timePoint = timePoints)
  return(estimates)
}
