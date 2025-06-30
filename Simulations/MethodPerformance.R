source("Simulations/SimulationFunctions.R")
source("Common/FunctionsForEstimands.R")

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
  source("Common/FunctionsForEstimands.R")
  return(NULL)
}
invisible(ParallelLogger::clusterApply(cluster, seq_len(maxCores), sourceFun))

simulateOne <- function(seed, settings) {
  population <- simulatePopulation(settings, seed)
  timePoints <- c(2, 7, 30, 90, 180, 365)
  estimates <- computeEstimands(population, timePoints = timePoints, bootstrapSize = 200) |>
    mutate(seed = seed)
  return(estimates)
}

# Simulate under the null --------------------------------------------------------------------------
settings <- createSimulationSettings(
  logHrFunction = function(t) {0}
)
estimates <- ParallelLogger::clusterApply(cluster, seq_len(sampleSize), simulateOne, settings = settings)
estimates <- bind_rows(estimates)

type1Error <- estimates |>
  group_by(timePoint, estimand, ciMethod) |>
  mutate(h0 = if_else(contrast == "ratio", 1, 0)) |>
  summarise(type1Error = mean(lb > h0 | ub < h0), .groups = "drop")
print(type1Error, n = 50)
# timePoint estimand                 ciMethod   type1Error
# <dbl> <chr>                    <chr>           <dbl>
# 1       180 Acceleration Coefficient asymptotic      0.051
# 2       180 Hazard Ratio             asymptotic      0.043
# 3       180 RMST Difference          asymptotic      0.051
# 4       180 RMST Difference          percentile      0.058
# 5       180 RMST Ratio               asymptotic      0.049
# 6       180 RMST Ratio               percentile      0.058
# 7       180 Risk Difference          asymptotic      0.056
# 8       180 Risk Difference          percentile      0.064
# 9       180 Risk Ratio               asymptotic      0.055
# 10       180 Risk Ratio               percentile      0.064
# 11       365 Acceleration Coefficient asymptotic      0.056
# 12       365 Hazard Ratio             asymptotic      0.044
# 13       365 RMST Difference          asymptotic      0.054
# 14       365 RMST Difference          percentile      0.068
# 15       365 RMST Ratio               asymptotic      0.054
# 16       365 RMST Ratio               percentile      0.068
# 17       365 Risk Difference          asymptotic      0.056
# 18       365 Risk Difference          percentile      0.068
# 19       365 Risk Ratio               asymptotic      0.056
# 20       365 Risk Ratio               percentile      0.068
# 21       730 Acceleration Coefficient asymptotic      0.059
# 22       730 Hazard Ratio             asymptotic      0.044
# 23       730 RMST Difference          asymptotic      0.073
# 24       730 RMST Difference          percentile      0.069
# 25       730 RMST Ratio               asymptotic      0.059
# 26       730 RMST Ratio               percentile      0.069
# 27       730 Risk Difference          asymptotic      0.056
# 28       730 Risk Difference          percentile      0.068
# 29       730 Risk Ratio               asymptotic      0.056
# 30       730 Risk Ratio               percentile      0.068
# 31      1095 Acceleration Coefficient asymptotic      0.059
# 32      1095 Hazard Ratio             asymptotic      0.044
# 33      1095 RMST Difference          asymptotic      0.077
# 34      1095 RMST Difference          percentile      0.069
# 35      1095 RMST Ratio               asymptotic      0.059
# 36      1095 RMST Ratio               percentile      0.069
# 37      1095 Risk Difference          asymptotic      0.056
# 38      1095 Risk Difference          percentile      0.068
# 39      1095 Risk Ratio               asymptotic      0.056
# 40      1095 Risk Ratio               percentile      0.068
# 41      1460 Acceleration Coefficient asymptotic      0.059
# 42      1460 Hazard Ratio             asymptotic      0.044
# 43      1460 RMST Difference          asymptotic      0.077
# 44      1460 RMST Difference          percentile      0.069
# 45      1460 RMST Ratio               asymptotic      0.059
# 46      1460 RMST Ratio               percentile      0.069
# 47      1460 Risk Difference          asymptotic      0.056
# 48      1460 Risk Difference          percentile      0.068
# 49      1460 Risk Ratio               asymptotic      0.056
# 50      1460 Risk Ratio               percentile      0.068

# Simulate when true effect one the (non-constant) hazard ratio scale ------------------------------
settings <- createSimulationSettings()
estimates <- ParallelLogger::clusterApply(cluster, seq_len(sampleSize), simulateOne, settings = settings)
estimates <- bind_rows(estimates)
type2Error <- estimates |>
  group_by(timePoint, estimand, ciMethod) |>
  mutate(h0 = if_else(contrast == "ratio", 1, 0)) |>
  summarise(type2Error = mean(lb <= h0 & ub >= h0), .groups = "drop")
print(type2Error, n = 50)
# timePoint estimand                 ciMethod   type2Error
# <dbl> <chr>                    <chr>           <dbl>
# 1       180 Acceleration Coefficient asymptotic      0.117
# 2       180 Hazard Ratio             asymptotic      0.098
# 3       180 RMST Difference          asymptotic      0.237
# 4       180 RMST Difference          percentile      0.224
# 5       180 RMST Ratio               asymptotic      0.24 
# 6       180 RMST Ratio               percentile      0.224
# 7       180 Risk Difference          asymptotic      0.543
# 8       180 Risk Difference          percentile      0.53 
# 9       180 Risk Ratio               asymptotic      0.549
# 10       180 Risk Ratio               percentile      0.53 
# 11       365 Acceleration Coefficient asymptotic      0.129
# 12       365 Hazard Ratio             asymptotic      0.103
# 13       365 RMST Difference          asymptotic      0.407
# 14       365 RMST Difference          percentile      0.393
# 15       365 RMST Ratio               asymptotic      0.409
# 16       365 RMST Ratio               percentile      0.393
# 17       365 Risk Difference          asymptotic      0.603
# 18       365 Risk Difference          percentile      0.59 
# 19       365 Risk Ratio               asymptotic      0.604
# 20       365 Risk Ratio               percentile      0.59 
# 21       730 Acceleration Coefficient asymptotic      0.129
# 22       730 Hazard Ratio             asymptotic      0.103
# 23       730 RMST Difference          asymptotic      0.45 
# 24       730 RMST Difference          percentile      0.468
# 25       730 RMST Ratio               asymptotic      0.481
# 26       730 RMST Ratio               percentile      0.468
# 27       730 Risk Difference          asymptotic      0.602
# 28       730 Risk Difference          percentile      0.59 
# 29       730 Risk Ratio               asymptotic      0.602
# 30       730 Risk Ratio               percentile      0.59 
# 31      1095 Acceleration Coefficient asymptotic      0.13 
# 32      1095 Hazard Ratio             asymptotic      0.103
# 33      1095 RMST Difference          asymptotic      0.446
# 34      1095 RMST Difference          percentile      0.468
# 35      1095 RMST Ratio               asymptotic      0.481
# 36      1095 RMST Ratio               percentile      0.468
# 37      1095 Risk Difference          asymptotic      0.602
# 38      1095 Risk Difference          percentile      0.59 
# 39      1095 Risk Ratio               asymptotic      0.602
# 40      1095 Risk Ratio               percentile      0.59 
# 41      1460 Acceleration Coefficient asymptotic      0.13 
# 42      1460 Hazard Ratio             asymptotic      0.103
# 43      1460 RMST Difference          asymptotic      0.446
# 44      1460 RMST Difference          percentile      0.468
# 45      1460 RMST Ratio               asymptotic      0.481
# 46      1460 RMST Ratio               percentile      0.468
# 47      1460 Risk Difference          asymptotic      0.602
# 48      1460 Risk Difference          percentile      0.59 
# 49      1460 Risk Ratio               asymptotic      0.602
# 50      1460 Risk Ratio               percentile      0.59 


# Simulate when true effect on the (non-constant) risk difference scale ----------------------------
settings <- createSimulationSettings(
  rdFunction = function(t) {dweibull(t, 1, 10)},
  logHrFunction = NULL
)
estimates <- ParallelLogger::clusterApply(cluster, seq_len(sampleSize), simulateOne, settings = settings)
estimates <- bind_rows(estimates)
type2Error <- estimates |>
  group_by(timePoint, estimand, ciMethod) |>
  mutate(h0 = if_else(contrast == "ratio", 1, 0)) |>
  summarise(type2Error = mean(lb <= h0 & ub >= h0), .groups = "drop")
print(type2Error, n = 50)
# # A tibble: 50 × 4
# timePoint estimand                 ciMethod   type2Error
# <dbl> <chr>                    <chr>           <dbl>
# 1       180 Acceleration Coefficient asymptotic      0.001
# 2       180 Hazard Ratio             asymptotic      0.001
# 3       180 RMST Difference          asymptotic      0.01 
# 4       180 RMST Difference          percentile      0.012
# 5       180 RMST Ratio               asymptotic      0.01 
# 6       180 RMST Ratio               percentile      0.012
# 7       180 Risk Difference          asymptotic      0.19 
# 8       180 Risk Difference          percentile      0.174
# 9       180 Risk Ratio               asymptotic      0.194
# 10       180 Risk Ratio               percentile      0.174
# 11       365 Acceleration Coefficient asymptotic      0.002
# 12       365 Hazard Ratio             asymptotic      0.001
# 13       365 RMST Difference          asymptotic      0.068
# 14       365 RMST Difference          percentile      0.063
# 15       365 RMST Ratio               asymptotic      0.068
# 16       365 RMST Ratio               percentile      0.063
# 17       365 Risk Difference          asymptotic      0.25 
# 18       365 Risk Difference          percentile      0.236
# 19       365 Risk Ratio               asymptotic      0.257
# 20       365 Risk Ratio               percentile      0.236
# 21       730 Acceleration Coefficient asymptotic      0.002
# 22       730 Hazard Ratio             asymptotic      0.001
# 23       730 RMST Difference          asymptotic      0.111
# 24       730 RMST Difference          percentile      0.113
# 25       730 RMST Ratio               asymptotic      0.126
# 26       730 RMST Ratio               percentile      0.113
# 27       730 Risk Difference          asymptotic      0.256
# 28       730 Risk Difference          percentile      0.243
# 29       730 Risk Ratio               asymptotic      0.263
# 30       730 Risk Ratio               percentile      0.243
# 31      1095 Acceleration Coefficient asymptotic      0.002
# 32      1095 Hazard Ratio             asymptotic      0.001
# 33      1095 RMST Difference          asymptotic      0.11 
# 34      1095 RMST Difference          percentile      0.114
# 35      1095 RMST Ratio               asymptotic      0.126
# 36      1095 RMST Ratio               percentile      0.114
# 37      1095 Risk Difference          asymptotic      0.256
# 38      1095 Risk Difference          percentile      0.243
# 39      1095 Risk Ratio               asymptotic      0.263
# 40      1095 Risk Ratio               percentile      0.243
# 41      1460 Acceleration Coefficient asymptotic      0.002
# 42      1460 Hazard Ratio             asymptotic      0.001
# 43      1460 RMST Difference          asymptotic      0.11 
# 44      1460 RMST Difference          percentile      0.114
# 45      1460 RMST Ratio               asymptotic      0.126
# 46      1460 RMST Ratio               percentile      0.114
# 47      1460 Risk Difference          asymptotic      0.256
# 48      1460 Risk Difference          percentile      0.243
# 49      1460 Risk Ratio               asymptotic      0.263
# 50      1460 Risk Ratio               percentile      0.243


ParallelLogger::stopCluster(cluster)
