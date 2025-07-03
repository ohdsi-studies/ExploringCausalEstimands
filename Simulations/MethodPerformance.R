source("Simulations/SimulationFunctions.R")
source("Common/FunctionsForEstimands.R")

library(survival)
library(ggplot2)
library(dplyr)
library(Cyclops)

maxCores <- 16
replications <- 1000
tempFolder <- "e:/temp/simTemp"


dir.create(tempFolder, showWarnings = FALSE)
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

allEstimates <- list()
for (trueEffectType in c("null", "multiplicative", "additive")) {
  for (depletionOfSusceptibles in c(TRUE, FALSE)) {
    for (sampleSize in c("big", "small")) {
      message("Simulating using trueEffectType = ",
              trueEffectType, 
              ", depletion of susceptibles = ", 
              depletionOfSusceptibles,
              ", sample size = ",
              sampleSize)
      
      fileName <- file.path(tempFolder, sprintf("Sim_%s_d%s%s.rds", 
                                                trueEffectType, 
                                                depletionOfSusceptibles,
                                                if (sampleSize == "big") "" else "_small"))
      if (file.exists(fileName)) {
        estimates <- readRDS(fileName)
        # estimates <- estimates |>
        #   mutate(sampleSize = !!sampleSize)
        # saveRDS(estimates, fileName)
      } else {
        settings <- createSimulationSettings()
        if (depletionOfSusceptibles) {
          settings$pSusceptible <- 0.2
        } else {
          settings$pSusceptible <- 1
        }
        if (trueEffectType == "null") {
          settings$logHrFunction <- function(t) {0}
          settings$rdFunction <- NULL
        } else if (trueEffectType == "multiplicative") {
          if (depletionOfSusceptibles) {
            settings$logHrFunction <- function(t) {25 * dweibull(t, 1.5, 10)}
          } else {
            settings$logHrFunction <- function(t) {5 * dweibull(t, 1.5, 10)}
          }
          settings$rdFunction <- NULL
        } else if (trueEffectType == "additive") {
          settings$logHrFunction <- NULL
          if (depletionOfSusceptibles) {
            settings$rdFunction <- function(t) {dweibull(t, 1, 10)}
          } else {
            settings$rdFunction <- function(t) {0.2 * dweibull(t, 1, 10)}
          }
        } else {
          stop("Unknown true effect type: ", trueEffectType)
        }
        if (sampleSize == "big") {
          settings$n <- 2500
        } else if (sampleSize == "small"){
          settings$n <- 1250
        } else {
          stop("Unknown sample size: ", sampleSize)
        }
        
        estimates <- ParallelLogger::clusterApply(cluster, seq_len(replications), simulateOne, settings = settings)
        estimates <- bind_rows(estimates)
        estimates <- estimates |>
          mutate(trueEffectType = !!trueEffectType,
                 depletionOfSusceptibles = !!depletionOfSusceptibles,
                 sampleSize = !!sampleSize)
        saveRDS(estimates, fileName)
      }
      allEstimates[[length(allEstimates) + 1]] <- estimates
    }
  }
}
allEstimates <- bind_rows(allEstimates)
ParallelLogger::stopCluster(cluster)



errorStats <- allEstimates |>
  mutate(h0 = if_else(contrast == "ratio", 1, 0)) |>
  mutate(error = if_else(trueEffectType == "null",
                         lb > h0 | ub < h0,
                         lb <= h0 & ub >= h0)) |>
  group_by(timePoint,
           estimand,
           model,
           contrast,
           ciMethod,
           trueEffectType,
           depletionOfSusceptibles,
           sampleSize) |>
  summarise(error = mean(error, na.rm = TRUE), .groups = "drop") |>
  mutate(errorType = if_else(trueEffectType == "null", "type 1", "type 2"))


library(ggplot2)
library(ggh4x)
colors <- c(
  "#ff7921", # orange
  "#94ad73", # matcha
  "#73655d", # gray
  "#d99b77" # brown
)

vizData <- errorStats |>
  filter(ciMethod == "asymptotic", ) |>
  mutate(Contrast = SqlRender::camelCaseToTitleCase(contrast),
         Model = model,
         trueEffectType = paste("True effect:", SqlRender::camelCaseToTitleCase(trueEffectType)),
         depletionOfSusceptibles = if_else(depletionOfSusceptibles, "Depletion of susceptibles", "No depletion of susceptibles"),
         sampleSize = if_else(sampleSize == "big", "2,500 patients", "1,250 patients"))
optimal <- tibble(
  trueEffectType = c("True effect: Null", "True effect: Multiplicative", "True effect: Additive"),
  y = c(0.05, 0, 0)
)
ggplot(vizData, aes(x = timePoint, y = error, group = estimand, color = Model)) +
  geom_hline(aes(yintercept = y), data = optimal) +
  geom_line(aes(linetype = Contrast), size = 1.25, alpha = 0.75) +
  scale_y_continuous("Error (type 1 or 2)") +
  scale_x_log10("Cutoff time (days)", breaks = unique(vizData$timePoint)) +
  # facet_grid(trueEffectType~depletionOfSusceptibles, scales = "free_y") +
  facet_nested(trueEffectType~depletionOfSusceptibles + sampleSize, scales = "free_y") +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "top"
  )
ggsave("Simulations/Type1And2Error.png", width = 9, height = 6, dpi = 300)

