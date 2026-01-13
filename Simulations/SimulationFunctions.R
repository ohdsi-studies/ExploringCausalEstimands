library(dplyr)

createSimulationSettings = function(
    n = 2500,                                                  # Number of persons to simulate
    pA = 0.5,                                                  # Probability of being exposed
    baselineHazardFunction = function(t) {dweibull(t, 1, 50)}, # Baseline hazard as a function of time
    baselineHazardMultiplier = runif(n, 0.5, 2),               # For each subject, their multiplier for the hazard function
    censorHazard = 0.01,                                       # Daily probability of being censored
    pEffectSusceptible = 0.20,                                 # Probability of being susceptible to exposure effect   
    pOutcomeSusceptible = 1,                                   # Probability of being susceptible to the outcome
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
    baselineHazardMultiplier = baselineHazardMultiplier,
    censorHazard = censorHazard,
    pEffectSusceptible = pEffectSusceptible,
    pOutcomeSusceptible = pOutcomeSusceptible,
    logHrFunction = logHrFunction,
    rdFunction = rdFunction
  )  
  return(settings)
}
# settings = createSimulationSettings()

simulatePopulation <- function(settings, seed = NULL) {
  set.seed(seed)
  a <- rbinom(settings$n, 1, settings$pA)
  effectSusceptible <- rbinom(settings$n, 1, settings$pEffectSusceptible)
  outcomeSusceptible <- rbinom(settings$n, 1, settings$pOutcomeSusceptible)

  atRisk <- rep(TRUE, settings$n)
  survivalTime <- rep(1001, settings$n)
  y <- rep(0, settings$n)
  trueHazardRatioOverTime <- rep(NA, 1000)
  trueCumulativeHazardRatioOverTime <- rep(NA, 1000)
  trueCumulativeHazardTarget <- 0
  trueCumulativeHazardComparator <- 0
  targetOverTime <- rep(NA, 1000)
  comparatorOverTime <- rep(NA, 1000)
  targetEffectSusceptiblesOverTime <- rep(NA, 1000)
  targetOutcomeSusceptiblesOverTime <- rep(NA, 1000)
  comparatorOutcomeSusceptiblesOverTime <- rep(NA, 1000)
  for (t in seq_len(1000)) {
    nAtRisk <- sum(atRisk)
    if (nAtRisk == 0) {
      break
    }
    targetOverTime[t] <- sum(atRisk & a)
    comparatorOverTime[t] <- sum(atRisk & !a)
    targetEffectSusceptiblesOverTime[t] <- sum(atRisk & a & effectSusceptible)
    targetOutcomeSusceptiblesOverTime[t] <- sum(atRisk & a & outcomeSusceptible)
    comparatorOutcomeSusceptiblesOverTime[t] <- sum(atRisk & !a & outcomeSusceptible)
    
    baselineHazards <- settings$baselineHazardFunction(t) * settings$baselineHazardMultiplier
    if (is.null(settings$logHrFunction)) {
      logHr <- 0
    } else {
      logHr <- settings$logHrFunction(t)
    }
    if (is.null(settings$rdFunction)) {
      rd <- 0
    } else {
      rd <- settings$rdFunction(t)
    }
    hazards <- if_else(a[atRisk] == 1 & effectSusceptible[atRisk] == 1, exp(logHr) * baselineHazards[atRisk] + rd, baselineHazards[atRisk])
    hazards <- if_else(outcomeSusceptible[atRisk] == 1, hazards, 0)
    meanHazardTarget <- mean(hazards[a[atRisk] == 1])
    meanHazardComparator <- mean(hazards[a[atRisk] != 1])
    trueHazardRatioOverTime[t] <- meanHazardTarget / meanHazardComparator
    trueCumulativeHazardTarget <- trueCumulativeHazardTarget + meanHazardTarget
    trueCumulativeHazardComparator <- trueCumulativeHazardComparator + meanHazardComparator
    trueCumulativeHazardRatioOverTime[t] <- trueCumulativeHazardTarget / trueCumulativeHazardComparator
    outcome <- runif(nAtRisk) < hazards
    censored <- runif(nAtRisk) < settings$censorHazard
    noLongerAtRisk <- outcome | censored
    survivalTime[atRisk][noLongerAtRisk] <- t
    y[atRisk] <- outcome
    atRisk[atRisk] <- !noLongerAtRisk
  }
  
  population <- tibble(
    a = a,
    survivalTime = survivalTime,
    y = y,
    effectSusceptible = effectSusceptible,
    outcomeSusceptible = outcomeSusceptible
  )
  attr(population, "trueHazardRatioOverTime") <- trueHazardRatioOverTime
  attr(population, "trueCumulativeHazardRatioOverTime") <- trueCumulativeHazardRatioOverTime
  attr(population, "targetOverTime") <- targetOverTime
  attr(population, "comparatorOverTime") <- comparatorOverTime
  attr(population, "targetEffectSusceptiblesOverTime") <- targetEffectSusceptiblesOverTime
  attr(population, "targetOutcomeSusceptiblesOverTime") <- targetOutcomeSusceptiblesOverTime
  attr(population, "comparatorOutcomeSusceptiblesOverTime") <- comparatorOutcomeSusceptiblesOverTime
  return(population)
}
