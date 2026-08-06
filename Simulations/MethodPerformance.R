# Code for running the simulations and analyzing the simulation results

# Run simulations ----------------------------------------------------------------------------------
source("Simulations/SimulationFunctions.R")
source("Common/FunctionsForEstimands.R")

library(survival)
library(dplyr)
library(Cyclops)

maxCores <- 10
replicationsMain <- 1000
replicationsAdditional <- 200
tempFolder <- "e:/temp/simTemp"


dir.create(tempFolder, showWarnings = FALSE)
ParallelLogger::addDefaultFileLogger(file.path(tempFolder, "log.txt"), name = "FILE_LOGGER")
cluster <- ParallelLogger::makeCluster(maxCores)
ParallelLogger::clusterRequire(cluster, "dplyr")
ParallelLogger::clusterRequire(cluster, "Cyclops")
sourceFun <- function(x) {
  source("Simulations/SimulationFunctions.R")
  source("Common/FunctionsForEstimands.R")
  return(NULL)
}
invisible(ParallelLogger::clusterApply(cluster, seq_len(maxCores), sourceFun))

simulateOne <- function(seed, settings, offset = 0) {
  population <- simulatePopulation(settings, seed)
  timePoints <- c(2, 7, 30, 90, 180, 365)
  if (offset != 0) {
    population <- population |>
      mutate(survivalTime = survivalTime - offset) |>
      filter(survivalTime > 0)
    timePoints <-timePoints - offset
    timePoints <- timePoints[timePoints > 0]
  } 
  if (length(unique(population$a)) < 2) {
    return(NULL)
  }
  estimates <- computeEstimands(population, timePoints = timePoints, bootstrapSize = 200) |>
    mutate(seed = seed)
  
  if (offset != 0) {
    estimates <- estimates |>
      mutate(timePoint = timePoint + offset)
  }
  return(estimates)
}

allEstimates <- list()
kmPlots <- list()
# trueEffectType = "null"
# depletionOfSusceptibles = "effect"
# sampleSize = "big"
# incidence = "high"
# frailtyHeterogeneity = "low"
# offset = 0
for (trueEffectType in c("null", "multiplicative", "additive", "multiplicativeLarge", "additiveLarge")) {
  for (depletionOfSusceptibles in c("effect", "outcome", "none")) {
    for (sampleSize in c("big", "small")) {
      for (incidence in c("high", "low")) {
        for (frailtyHeterogeneity in c("low", "high")) {
          for (offset in c(0, 14)) {
            message("Simulating using trueEffectType = ",
                    trueEffectType, 
                    ", depletion of susceptibles = ", 
                    depletionOfSusceptibles,
                    ", sample size = ",
                    sampleSize,
                    ", incidence = ",
                    incidence,
                    ", frailtyHeterogeneity = ",
                    frailtyHeterogeneity,
                    ", offset = ",
                    offset)
            
            isMain <- trueEffectType %in% c("null", "multiplicative", "additive") &&
              incidence == "high" &&
              frailtyHeterogeneity == "low" &&
              offset == 0
            if (isMain) {
              replications <- replicationsMain
            } else {
              replications <- replicationsAdditional
            }
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
            } else if (trueEffectType == "multiplicativeLarge") {
              if (depletionOfSusceptibles != "none") {
                settings$logHrFunction <- function(t) {250 * dweibull(t, 1.5, 10)}
              } else {
                settings$logHrFunction <- function(t) {50 * dweibull(t, 1.5, 10)}
              }
              settings$rdFunction <- NULL
            } else if (trueEffectType == "additiveLarge") {
              settings$logHrFunction <- NULL
              if (depletionOfSusceptibles != "none") {
                settings$rdFunction <- function(t) {10 * dweibull(t, 1, 10)}
              } else {
                settings$rdFunction <- function(t) {2 * dweibull(t, 1, 10)}
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
            if (incidence == "high") {
              settings$baselineHazardFunction <- function(t) {dweibull(t, 1, 50)}
            } else if (incidence == "low") {
              settings$baselineHazardFunction <- function(t) {0.1 * dweibull(t, 1, 50)}
            } else {
              stop("Unknown incidence: ", incidence)
            }
            if (frailtyHeterogeneity == "low") {
              settings$baselineHazardMultiplier = runif(settings$n, 0.5, 2)
            } else if (frailtyHeterogeneity == "high") {
              settings$baselineHazardMultiplier = runif(settings$n, 0.1, 10)
            } else {
              stop("Unknown frailtyHeterogeneity: ", frailtyHeterogeneity)
            }
            key <- tibble(
              trueEffectType = trueEffectType,
              depletionOfSusceptibles = depletionOfSusceptibles,
              sampleSize = sampleSize,
              incidence = incidence,
              frailtyHeterogeneity = frailtyHeterogeneity,
              offset = offset
            )
            
            # Monte-Carlo simulations
            fileName <- file.path(tempFolder, sprintf("Sim_%s_%s%s%s%s%s.rds", 
                                                      trueEffectType, 
                                                      depletionOfSusceptibles,
                                                      if (sampleSize == "big") "" else "_smallSample",
                                                      if (incidence == "high") "" else "_lowIncidence",
                                                      if (frailtyHeterogeneity == "low") "" else  "_highFh",
                                                      if (offset == 0) "" else "_offset14"))
            if (file.exists(fileName)) {
              estimates <- readRDS(fileName)
              # estimates$offset = offset
              # estimates$incidence = incidence
              # estimates$frailtyHeterogeneity = frailtyHeterogeneity
              # saveRDS(estimates, fileName)
            } else {
              message(sprintf("- Using %d replications", replications))
              estimates <- ParallelLogger::clusterApply(cluster, seq_len(replications), simulateOne, settings = settings, offset = offset)
              estimates <- bind_rows(estimates)
              estimates <- key |>
                bind_cols(estimates)
              saveRDS(estimates, fileName)
            }
            allEstimates[[length(allEstimates) + 1]] <- estimates
            
            # KM plot
            if (offset == 0) {
              fileName <- file.path(tempFolder, sprintf("KM_%s_%s%s%s%s.pdf", 
                                                        trueEffectType, 
                                                        depletionOfSusceptibles,
                                                        if (sampleSize == "big") "" else "_smallSample",
                                                        if (incidence == "high") "" else "_lowIncidence",
                                                        if (frailtyHeterogeneity == "low") "" else  "_highFh"))
              label <- sprintf("## %s, %s, %s, %s, %s", 
                               paste(gsub("Large", "(Large)", SqlRender::camelCaseToTitleCase(trueEffectType)), "true effect"), 
                               paste(gsub("none", "no", depletionOfSusceptibles), "depletion"),
                               if_else(sampleSize == "big", "2,500 patients", "1,250 patients"),
                               paste(incidence, "incidence"),
                               paste(frailtyHeterogeneity, "frailty heterogeneity"))
              if (!file.exists(fileName)) {
                population <- simulatePopulation(settings, seed = 1) |>
                  rename(treatment = a,
                         outcomeCount = y)
                CohortMethod::plotKaplanMeier(population, fileName = fileName)
              }
              kmPlots[[length(kmPlots) + 1]] <- key |>
                mutate(label = !!label,
                       fileName = !!fileName)
            }
          }
        }
      }
    }
  }
}

allEstimates <- bind_rows(allEstimates)
kmPlots <- bind_rows(kmPlots)
ParallelLogger::stopCluster(cluster)
ParallelLogger::unregisterLogger("FILE_LOGGER")
saveRDS(allEstimates, "Simulations/allEstimates.rds")
saveRDS(kmPlots, file.path(tempFolder, "kmPlots.rds"))

# Plot type 1 and 2 error ------------------------------------------------------
library(ggplot2)
library(ggh4x)
library(dplyr)
library(RColorBrewer)

allEstimates <- readRDS("Simulations/allEstimates.rds")

# Note: type 2 error computation looks at 1 side only. We don't give credit for rejection the null on the wrong side:
errorStats <- allEstimates |>
  filter((ciMethod == "percentile" | model == "Cox" | model == "AFT")) |>
  mutate(model = if_else(estimand == "Cumulative Hazard Ratio", "Cumulative Hazard (KM)", model)) |>
  mutate(model = if_else(model == "Kaplan Meier", "Risk (KM)", model),
         temp = lb) |>
  mutate(h0 = if_else(contrast == "ratio", 1, 0),
         lb = if_else(model == "AFT", 1/ub, lb),
         ub = if_else(model == "AFT", 1/temp, ub)) |>
  select(-"temp") |>
  mutate(error = if_else(trueEffectType == "null",
                         lb > h0 | ub < h0,
                         lb <= h0),
         flipError = ub < h0) |>
  group_by(timePoint,
           estimand,
           model,
           contrast,
           ciMethod,
           trueEffectType,
           depletionOfSusceptibles,
           sampleSize,
           incidence,
           frailtyHeterogeneity,
           offset) |>
  summarise(error = mean(error, na.rm = TRUE), 
            flipError = mean(flipError, na.rm = TRUE),
            .groups = "drop") |>
  mutate(errorType = if_else(trueEffectType == "null", "type 1", "type 2"))

# Plot subset for paper:
vizData <- errorStats |>
  filter(trueEffectType %in% c("null", "multiplicative", "additive"),
         offset == 0,
         incidence == "high",
         frailtyHeterogeneity == "low") |>
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

ggplot(vizData, aes(x = timePoint, y = error, group = estimand, color = Model)) +
  geom_hline(aes(yintercept = y), data = optimal) +
  geom_line(aes(linetype = Contrast), linewidth = 0.75, alpha = 0.65) +
  scale_y_continuous("Error (type 1 or 2)") +
  scale_x_log10("Cutoff time (days)", breaks = unique(vizData$timePoint)) +
  scale_color_manual(values = brewer.pal(5, "Set1")) +
  facet_nested(trueEffectType ~ depletionOfSusceptibles + sampleSize) +
  theme(
    axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.5, color = "white"),
    legend.position = "top"
  )
ggsave("Simulations/plots/Type1And2Error.png", width = 8, height = 4.5, dpi = 300)
ggsave("Simulations/plots/Type1And2Error.svg", width = 8, height = 4.5)
ggsave("Simulations/plots/Type1And2Error.pdf", width = 8, height = 4.5)

# Plots for supplement:
# frailtyHeterogeneity = "low"
# offset = 14
for (frailtyHeterogeneity in c("low", "high")) {
  for (offset in c(0)) {
    vizData <- errorStats |>
      filter(offset == !!offset,
             frailtyHeterogeneity == !!frailtyHeterogeneity) |>
      mutate(Contrast = SqlRender::camelCaseToTitleCase(contrast),
             Model = model,
             trueEffectType = paste("True effect:", gsub("Large", " (Large)", SqlRender::camelCaseToTitleCase(trueEffectType)), sep = "\n"),
             depletionOfSusceptibles = case_when(depletionOfSusceptibles == "outcome" ~ "Depletion of outcome susceptibles", 
                                                 depletionOfSusceptibles == "effect" ~ "Depletion of effect susceptibles", 
                                                 .default = "No depletion of susceptibles"),
             incidence = if_else(incidence == "high", "High\nincidence", "Low\nincidence"),
             sampleSize = if_else(sampleSize == "big", "2,500 patients", "1,250 patients"))
    optimal <- vizData |>
      distinct(trueEffectType, incidence) |>
      mutate(y = if_else(grepl("Null", trueEffectType), 0.05, 0))
    
    
    ggplot(vizData, aes(x = timePoint, y = error, group = estimand, color = Model)) +
      geom_hline(aes(yintercept = y), data = optimal) +
      geom_line(aes(linetype = Contrast), linewidth = 0.75, alpha = 0.65) +
      scale_y_continuous("Error (type 1 or 2)") +
      scale_x_log10("Cutoff time (days)", breaks = unique(vizData$timePoint)) +
      scale_color_manual(values = brewer.pal(5, "Set1")) +
      facet_nested(trueEffectType + incidence ~ depletionOfSusceptibles + sampleSize) +
      theme(
        axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.5, color = "white"),
        legend.position = "top"
      )
    baseName <- sprintf("Simulations/plots/Type1And2Error%s%s", frailtyHeterogeneity, offset)
    ggsave(paste0(baseName, ".png"), width = 8, height = 10, dpi = 300)
    ggsave(paste0(baseName, ".svg"), width = 8, height = 10)
    ggsave(paste0(baseName, ".pdf"), width = 8, height = 10)
  }
}
