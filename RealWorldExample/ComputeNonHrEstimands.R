source("RealWorldExample/SetConnectionDetails.R")
source("Common/FunctionsForNonHrEstimands.R")
library(CohortMethod)
library(survival)
library(dplyr)

timePoints <- c(180, 365, 730, 1095, 1460)
bootstrapSize <- 1000

cluster <- ParallelLogger::makeCluster(maxCores)

ref <- getFileReference(outputFolder)

# Using weighting --------------------------------------------------------------
allEstimates <- list()
for (i in 1:nrow(ref)) {
  message(sprintf("Computing estimates for %d of %d", i , nrow(ref)))
  population <- readRDS(file.path(outputFolder, ref$psFile[i]))
  population$y <- population$outcomeCount > 0
  population$weight <- computeWeights(population)
  estimates <- computeEstimands(
    population = population,
    timePoints = timePoints,
    bootstrapSize = bootstrapSize,
    cluster = cluster
  )
  estimate <- estimates |>
    mutate(ncoId = i)
  allEstimates[[i]] <- estimates
}
allEstimates <- do.call(rbind, allEstimates)
saveRDS(estimates, "RealWorldExample/rrEstimatesWeighted.rds")

# Using matching ---------------------------------------------------------------
allEstimates <- list()
for (i in 1:nrow(ref)) {
  message(sprintf("Computing estimates for %d of %d", i , nrow(ref)))
  population <- readRDS(file.path(outputFolder, ref$strataFile[i]))
  population$y <- population$outcomeCount > 0
  estimates <- computeEstimands(
    population = population,
    timePoints = timePoints,
    bootstrapSize = bootstrapSize,
    cluster = cluster
  )
  estimates <- estimates |>
    mutate(ncoId = i)
  allEstimates[[i]] <- estimates
}
allEstimates <- do.call(rbind, allEstimates)
saveRDS(estimates,  "RealWorldExample/rrEstimatesMatching.rds")


ParallelLogger::stopCluster(cluster)
