library(ggplot2)
source("Simulations/SimulationFunctions.R")
source("Common/FunctionsForEstimands.R") 
source("ExperimentalEstimands/RiskDifference.R") 
source("ExperimentalEstimands/BayesianRiskDifference.R") 
source("ExperimentalEstimands/StabilizedBayesianRiskDifference.R") 

cluster <- ParallelLogger::makeCluster(10)

settings <- createSimulationSettings()
settings <- createSimulationSettings(
  n = 10000,
  baselineHazardFunction = function(t) {0.0001},
  logHrFunction = NULL,
  rdFunction = function(t) {0.25 * dweibull(t, 2, 120)}
)
plot(1:300, settings$rdFunction(1:300))


population <- simulatePopulation(settings)

# Plot KM:
popForKm <- population |>
  mutate(outcomeCount = y,
         treatment = a)
CohortMethod::plotKaplanMeier(popForKm, dataCutoff = 0.99, fileName = "ExperimentalEstimands/Km.png")


timePoints <- c(2, 30, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360)

estimates <- computeEstimands(population = population, timePoints = timePoints, cluster = cluster)
estimatesRd <- estimateRiskDifference(population = population, timePoints = timePoints)
estimatesBrd <- estimateBayesianRiskDifference(population = population, timePoints = timePoints)
estimatesSbrd <- estimateStablizedBayesianRiskDifference(population, max(timePoints))
vizData <- bind_rows(
  estimates |>
    filter(estimand == "Hazard Ratio") |>
    mutate(type = "Cox"),
  estimates |>
    filter(estimand == "Risk Difference", ciMethod == "percentile") |>
    mutate(type = "Risk Difference (bootstrap)"),
  estimatesRd |>
    mutate(type = "Risk Difference (EL)"),
  estimatesBrd |>
    mutate(type = "Risk Difference (Bayesian)"),
  estimatesSbrd |>
    mutate(type = "Risk Difference (Stabilized Bayesian)")
)
reference <- tibble(type = c("Cox", "Risk Difference (bootstrap)", "Risk Difference (EL)", "Risk Difference (Bayesian)", "Risk Difference (Stabilized Bayesian)"),
                    h0 = c(1, 0, 0, 0, 0))
ggplot(vizData, aes(x = timePoint, y = estimate)) +
  geom_hline(aes(yintercept = h0), data = reference) +
  geom_errorbar(aes(ymin = lb, ymax = ub)) +
  geom_point() +
  scale_y_continuous("Time point") +
  scale_y_continuous("Effect estimate") +
  facet_grid(type ~ ., scales = "free_y")
ggsave("ExperimentalEstimands/ExperimentalEstimands.png")

 ParallelLogger::stopCluster(cluster)




