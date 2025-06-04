source("Simulations/SimulationFunctions.R")
source("Common/FunctionsForNonHrEstimands.R")

library(survival)
library(ggplot2)
library(dplyr)
library(Cyclops)

maxCores <- 10
sampleSize <- 1000

cluster <- ParallelLogger::makeCluster(maxCores)
ParallelLogger::clusterRequire(cluster, "dplyr")
ParallelLogger::clusterRequire(cluster, "Cyclops")
sourceFun <- function(x) {
  source("Simulations/SimulationFunctions.R")
  source("Common/FunctionsForNonHrEstimands.R")
  return(NULL)
}
invisible(ParallelLogger::clusterApply(cluster, seq_len(maxCores), sourceFun))

simulateOne <- function(seed, settings) {
  population <- simulatePopulation(settings, seed)
  
  cyclopsData <- createCyclopsData(Surv(survivalTime, y) ~ a, modelType = "cox", data = population)
  fit <- fitCyclopsModel(cyclopsData)
  hr <- exp(coef(fit))
  ci <- exp(confint(fit, parm = "a")[c(2, 3)])
  se <- (log(ci[2]) - log(ci[1])) / (qnorm(0.975) - qnorm(0.025))
  
  timePoints <- c(180, 365, 730, 1095, 1460)
  kmEstimates <-computeEstimands(population, timePoints = timePoints, bootstrapSize = 200)
  
  rows <-
    kmEstimates |>
    mutate(
      seed = seed,
      hr = hr,
      lbHr = ci[1],
      ubHr = ci[2],
      seLogHr = se
    )
  return(rows)
}

# Simulate under the null --------------------------------------------------------------------------
settings <- createSimulationSettings(
  logHrFunction = function(t) {0}
)
estimates <- ParallelLogger::clusterApply(cluster, seq_len(sampleSize), simulateOne, settings = settings)
estimates <- bind_rows(estimates)

estimates |>
  group_by(timePoint) |>
  summarise(type1ErrorHr = mean(lbHr > 1 | ubHr < 1),
            type1ErrorRr = mean(lbRrAsymptotics > 1 | ubRrAsymptotics < 1),
            type1ErrorRd = mean(lbRdAsymptotics > 0 | ubRdAsymptotics < 0))
# # A tibble: 5 × 4
# timePoint type1ErrorHr type1ErrorRr type1ErrorRd
# <dbl>        <dbl>        <dbl>        <dbl>
# 1       180        0.055        0.064        0.064
# 2       365        0.055        0.06         0.056
# 3       730        0.055        0.06         0.056
# 4      1095        0.055        0.06         0.056
# 5      1460        0.055        0.06         0.056

# Simulate when true effect one the (non-constant) hazard ratio scale ------------------------------
settings <- createSimulationSettings()
estimates <- ParallelLogger::clusterApply(cluster, seq_len(sampleSize), simulateOne, settings = settings)
estimates <- bind_rows(estimates)

estimates |>
  group_by(timePoint) |>
  summarise(type2ErrorHr = mean(lbHr < 1 & ubHr > 1),
            type2ErrorRr = mean(lbRrAsymptotics < 1 & ubRrAsymptotics > 1),
            type2ErrorRd = mean(lbRdAsymptotics < 0 & ubRdAsymptotics > 0),
            type2ErrorRrPerc = mean(lbRrPercentile < 1 & ubRrPercentile > 1),
            type2ErrorRdPerc = mean(lbRdPercentile < 0 & ubRdPercentile > 0))
# # A tibble: 5 × 6
# timePoint type2ErrorHr type2ErrorRr type2ErrorRd type2ErrorRrPerc type2ErrorRdPerc
# <dbl>        <dbl>        <dbl>        <dbl>            <dbl>            <dbl>
# 1       180        0.024        0.549        0.543            0.534            0.534
# 2       365        0.024        0.614        0.611            0.601            0.601
# 3       730        0.024        0.615        0.612            0.602            0.602
# 4      1095        0.024        0.615        0.612            0.602            0.602
# 5      1460        0.024        0.615        0.612            0.602            0.602

# Simulate when true effect one the (non-constant) risk difference scale ---------------------------
settings <- createSimulationSettings(
  rdFunction = function(t) {dweibull(t, 1, 10)},
  logHrFunction = NULL
)
estimates <- ParallelLogger::clusterApply(cluster, seq_len(sampleSize), simulateOne, settings = settings)
estimates <- bind_rows(estimates)

estimates |>
  group_by(timePoint) |>
  summarise(type2ErrorHr = mean(lbHr < 1 & ubHr > 1),
            type2ErrorRr = mean(lbRrAsymptotics < 1 & ubRrAsymptotics > 1),
            type2ErrorRd = mean(lbRdAsymptotics < 0 & ubRdAsymptotics > 0),
            type2ErrorRrPerc = mean(lbRrPercentile < 1 & ubRrPercentile > 1),
            type2ErrorRdPerc = mean(lbRdPercentile < 0 & ubRdPercentile > 0))
# # A tibble: 5 × 6
# timePoint type2ErrorHr type2ErrorRr type2ErrorRd type2ErrorRrPerc type2ErrorRdPerc
# <dbl>        <dbl>        <dbl>        <dbl>            <dbl>            <dbl>
#   1       180        0.024        0.542        0.537            0.526            0.526
# 2       365        0.024        0.612        0.609            0.582            0.582
# 3       730        0.024        0.618        0.616            0.585            0.585
# 4      1095        0.024        0.618        0.616            0.585            0.585
# 5      1460        0.024        0.618        0.616            0.585            0.585

ParallelLogger::stopCluster(cluster)



# population <- simulatePopulation(settings, seed)
