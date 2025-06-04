library(dplyr)

createSimulationSettings = function(
    n = 2500,                                                  # Number of persons to simulate
    pA = 0.5,                                                  # Probability of being exposed
    baselineHazardFunction = function(t) {dweibull(t, 1, 50)}, # Baseline hazard as a function of time
    censorHazard = 0.01,                                       # Daily probability of being censored
    pSusceptible = 0.20,                                       # Probability of being susceptible to exposure effect   
    logHrFunction = function(t) {25 * dweibull(t, 1.5, 10)},   # Log(hazard ratio) as a function of time 
    rdFunction = NULL                                          # Risk difference as a function of time
) {
  if (!is.null(logHrFunction) && !is.null(rdFunction)) {
    stop("Cannot specify both logHrFunction and rdFunction. One of them needs to be NULL")
  }
  settings <- list(
    n = n,
    pA = pA,
    baselineHazardFunction = baselineHazardFunction,
    censorHazard = censorHazard,
    pSusceptible = pSusceptible,
    logHrFunction = logHrFunction,
    rdFunction = rdFunction
  )  
  return(settings)
}
# settings = createSimulationSettings()

simulatePopulation <- function(settings, seed = NULL) {
  set.seed(seed)
  a <- rbinom(settings$n, 1, settings$pA)
  susceptible <- rbinom(settings$n, 1, settings$pSusceptible)
  
  atRisk <- rep(TRUE, settings$n)
  survivalTime <- rep(1001, settings$n)
  y <- rep(0, settings$n)
  averageLogHr <- 0
  denominator <- 0
  targetOverTime <- rep(NA, 1000)
  targetSusceptiblesOverTime <- rep(NA, 1000)
  for (t in seq_len(1000)) {
    nAtRisk <- sum(atRisk)
    if (nAtRisk == 0) {
      break
    }
    targetOverTime[t] <- sum(atRisk & a)
    targetSusceptiblesOverTime[t] <- sum(atRisk & a & susceptible)
    
    baselineHazard <- settings$baselineHazardFunction(t)
    if (is.null(settings$logHrFunction)) {
      logHr <- 0
    } else {
      logHr <- settings$logHrFunction(t)
      averageLogHr <- averageLogHr + logHr * nAtRisk
    }
    if (is.null(settings$rdFunction)) {
      rd <- 0
    } else {
      rd <- settings$rdFunction(t)
    }
    denominator <- denominator + nAtRisk
    hazard <- ifelse(a[atRisk] == 1 & susceptible[atRisk] == 1, exp(logHr) * baselineHazard + rd, baselineHazard)
    outcome <- runif(nAtRisk) < hazard
    censored <- runif(nAtRisk) < settings$censorHazard
    noLongerAtRisk <- outcome | censored
    survivalTime[atRisk][noLongerAtRisk] <- t
    y[atRisk] <- outcome
    atRisk[atRisk] <- !noLongerAtRisk
  }
  averageLogHr <- averageLogHr / denominator
  
  population <- tibble(
    a = a,
    survivalTime = survivalTime,
    y = y,
    susceptible = susceptible
  )
  attr(population, "targetOverTime") <- targetOverTime
  attr(population, "targetSusceptiblesOverTime") <- targetSusceptiblesOverTime
  if (!is.null(settings$logHrFunction)) {
    attr(population, "averageLogHr") <- averageLogHr
  }
  return(population)
}
