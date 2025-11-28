source("Simulations/SimulationFunctions.R")
source("Common/FunctionsForEstimands.R")

library(survival)
library(dplyr)
library(Cyclops)

maxCores <- 13
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
  for (depletionOfSusceptibles in c("effect", "outcome", "none")) {
    for (sampleSize in c("big", "small")) {
      message("Simulating using trueEffectType = ",
              trueEffectType, 
              ", depletion of susceptibles = ", 
              depletionOfSusceptibles,
              ", sample size = ",
              sampleSize)
      
      fileName <- file.path(tempFolder, sprintf("Sim_%s_%s_%s.rds", 
                                                trueEffectType, 
                                                depletionOfSusceptibles,
                                                if (sampleSize == "big") "" else "_small"))
      if (file.exists(fileName)) {
        estimates <- readRDS(fileName)
      } else {
        settings <- createSimulationSettings()
        if (depletionOfSusceptibles == "effect") {
          settings$pEffectSusceptible <- 0.2
          settings$pOutcomeSusceptible <- 1
        } else if (depletionOfSusceptibles == "outcome") {
          settings$pEffectSusceptible <- 1
          settings$pOutcomeSusceptible <- 0.2
        } else if (depletionOfSusceptibles == "none") {
          settings$pEffectSusceptible <- 1
          settings$pOutcomeSusceptible <- 1
        } else {
          stop("Unknown depletion of susceptibles: ", depletionOfSusceptibles)
        }
        if (trueEffectType == "null") {
          settings$logHrFunction <- function(t) {0}
          settings$rdFunction <- NULL
        } else if (trueEffectType == "multiplicative") {
          if (depletionOfSusceptibles != "none") {
            settings$logHrFunction <- function(t) {25 * dweibull(t, 1.5, 10)}
          } else {
            settings$logHrFunction <- function(t) {5 * dweibull(t, 1.5, 10)}
          }
          settings$rdFunction <- NULL
        } else if (trueEffectType == "additive") {
          settings$logHrFunction <- NULL
          if (depletionOfSusceptibles != "none") {
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
        settings$baselineHazardMultiplier <- runif(settings$n, 0.5, 2)
        
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
saveRDS(allEstimates, "Simulations/allEstimates.rds")

# Plot type 1 and 2 error ------------------------------------------------------
library(ggplot2)
library(ggh4x)
library(RColorBrewer)
library(dplyr)

allEstimates <- readRDS("Simulations/allEstimates.rds")

errorStats <- allEstimates |>
  filter((ciMethod == "percentile" | model == "Cox" | model == "AFT")) |>
  mutate(model = if_else(estimand == "Cumulative Hazard Ratio", "Cumulative Hazard (KM)", model)) |>
  mutate(model = if_else(model == "Kaplan Meier", "Risk (KM)", model)) |>
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

vizData <- errorStats |>
  mutate(Contrast = SqlRender::camelCaseToTitleCase(contrast),
         Model = model,
         trueEffectType = paste("True effect:", SqlRender::camelCaseToTitleCase(trueEffectType), sep = "\n"),
         depletionOfSusceptibles = case_when(depletionOfSusceptibles == "outcome" ~ "Depletion of outcome susceptibles", 
                                             depletionOfSusceptibles == "effect" ~ "Depletion of effect susceptibles", 
                                             .default = "No depletion of susceptibles"),
         sampleSize = if_else(sampleSize == "big", "2,500 patients", "1,250 patients"))
optimal <- tibble(
  trueEffectType = c("True effect:\nNull", "True effect:\nMultiplicative", "True effect:\nAdditive"),
  y = c(0.05, 0, 0)
)

library(wesanderson)
ggplot(vizData, aes(x = timePoint, y = error, group = estimand, color = Model)) +
  geom_hline(aes(yintercept = y), data = optimal) +
  geom_line(aes(linetype = Contrast), linewidth = 0.75, alpha = 0.65) +
  scale_y_continuous("Error (type 1 or 2)") +
  scale_x_log10("Cutoff time (days)", breaks = unique(vizData$timePoint)) +
  scale_color_manual(values = wes_palette("Darjeeling1", 5, type = "discrete")) +
  # facet_nested(depletionOfSusceptibles + sampleSize ~ trueEffectType, scales = "free_y") +
  facet_nested(trueEffectType ~ depletionOfSusceptibles + sampleSize) +
  theme(
    axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.5, color = "white"),
    legend.position = "top"
  )
ggsave("Simulations/Type1And2Error.png", width = 8, height = 4.5, dpi = 300)
ggsave("Simulations/Type1And2Error.svg", width = 8, height = 4.5)

# ggsave("Simulations/Type1And2Error.png", width = 7, height = 7, dpi = 300)
# ggsave("Simulations/Type1And2Error.svg", width = 7, height = 7, dpi = 300)

e90 <- allEstimates |>
  filter(trueEffectType == "multiplicative",
         depletionOfSusceptibles == "none",
         sampleSize == "big",
         ciMethod == "asymptotic",
         timePoint == 90,
         model == "Kaplan Meier",
         contrast == "difference")

e30 <- allEstimates |>
  filter(trueEffectType == "multiplicative",
         depletionOfSusceptibles == "none",
         sampleSize == "big",
         ciMethod == "asymptotic",
         timePoint == 30,
         model == "Kaplan Meier",
         contrast == "difference")
e90
e30

merged <- inner_join(
  e30 |>
    select(seed, estimate30 = estimate, lb30 = lb, ub30 = ub, se30 = se),
  e90 |>
    select(seed, estimate90 = estimate, lb90 = lb, ub90 = ub, se90 = se),
  by = join_by(seed)
)
ggplot(merged, aes(x = se30, y = se90)) +
  geom_abline(slope = 1) +
  geom_point(alpha = 0.4)

ggplot(merged, aes(x = estimate30, y = estimate90)) +
  geom_abline(slope = 1) +
  geom_point(alpha = 0.4)


e <- allEstimates |>
  filter(trueEffectType == "multiplicative",
         depletionOfSusceptibles == "none",
         sampleSize == "big",
         ciMethod == "asymptotic",
         model %in% c("Kaplan Meier", "Cox"))

ggplot(e, aes(x = timePoint, y = se, group = timePoint)) +
  geom_boxplot() +
  scale_x_log10() +
  scale_y_continuous("Standard error") +
  facet_nested(contrast + model~., scale = "free_y")

ggplot(e, aes(x = timePoint, y = estimate, group = timePoint)) +
  geom_boxplot() +
  scale_x_log10() +
  scale_y_continuous("Point estimate") +
  facet_nested(contrast + model~., scale = "free_y")


e |>
  filter(contrast == "ratio") |>
  mutate(newLb = exp(log(estimate) + qnorm(0.025) * se)) |>
  group_by(model, timePoint) |>
  summarise(mean(newLb < 1))

vizData <- e |>
  filter(seed %in% 1:4, contrast == "ratio") |>
  mutate(seed = paste("Seed", seed),
         model = if_else(estimand == "Cumulative Hazard Ratio", "Cumulative Hazard", model))


ggplot(vizData, aes(x = timePoint, y = estimate)) +
  geom_hline(yintercept = 1) +
  geom_errorbar(aes(ymin = lb, ymax = ub)) +
  geom_point() +
  scale_x_log10("Time point") +
  scale_y_continuous("Ratio") +
  coord_cartesian(ylim = c(0.75, 1.75)) +
  facet_grid(seed ~ model)
ggsave("Simulations/CoxVsKmEstimates.png", width = 6, height = 6)
