source("RealWorldExample/SetConnectionDetails.R")
source("Common/FunctionsForEstimands.R")
library(CohortMethod)
library(survival)
library(dplyr)

timePoints <- c(2, 7, 30, 90, 180, 365)
sampleSizes <- c(200000, 100000, 50000, 25000, 12500)
bootstrapSize <- 200

cluster <- ParallelLogger::makeCluster(maxCores)

ref <- getFileReference(outputFolder)

tempFolder <- file.path(outputFolder, "estimatesTemp")
if (!dir.exists(tempFolder)) {
  dir.create(tempFolder)
}

# Using weighting --------------------------------------------------------------
# allEstimates <- list()
# for (i in 1:nrow(ref)) {
#   message(sprintf("Computing estimates for %d of %d", i , nrow(ref)))
#   fileName <- file.path(tempFolder, sprintf("estimatesWeighted_t%d_c%d_o%d.rds", 
#                                             ref$targetId[i],
#                                             ref$comparatorId[i],
#                                             ref$outcomeId[i]))
#   if (file.exists(fileName)) {
#     estimates <- readRDS(fileName)
#   } else {
#     population <- readRDS(file.path(outputFolder, ref$psFile[i]))
#     population$y <- population$outcomeCount > 0
#     population$weight <- computeWeights(population)
#     estimates <- computeEstimands(
#       population = population,
#       timePoints = timePoints,
#       bootstrapSize = bootstrapSize,
#       cluster = cluster
#     )
#     estimates <- estimates |>
#       mutate(targetId = ref$targetId[i],
#              comparatorId = ref$comparatorId[i],
#              outcomeId = ref$outcomeId[i])
#     saveRDS(estimates, fileName)
#   }
#   allEstimates[[i]] <- estimates
# }
# allEstimates <- do.call(rbind, allEstimates)
# saveRDS(allEstimates, "RealWorldExample/estimatesWeighted.rds")

# Using matching ---------------------------------------------------------------
allEstimates <- list()
for (i in 1:nrow(ref)) {
  # sampleSize = sampleSizes[1]
  for (sampleSize in sampleSizes) {
    
    message(sprintf("Computing estimates for %d of %d at sample size %d", i , nrow(ref), sampleSize))
    
    
    fileName <- file.path(tempFolder, sprintf("estimatesMatched_t%d_c%d_o%d_s%d.rds", 
                                              ref$targetId[i],
                                              ref$comparatorId[i],
                                              ref$outcomeId[i],
                                              sampleSize))
    if (file.exists(fileName)) {
      estimates <- readRDS(fileName)
    } else {
      population <- readRDS(file.path(outputFolder, ref$strataFile[i]))
      sampledStratumIds <- sample(unique(population$stratumId), sampleSize / 2)
      population <- population |>
        filter(stratumId %in% sampledStratumIds)
      population$y <- population$outcomeCount > 0
      estimates <- computeEstimands(
        population = population,
        timePoints = timePoints,
        bootstrapSize = bootstrapSize,
        cluster = cluster
      )
      estimates <- estimates |>
        mutate(targetId = ref$targetId[i],
               comparatorId = ref$comparatorId[i],
               outcomeId = ref$outcomeId[i],
               sampleSize = !!sampleSize,
               nOutcomes = sum(population$y))
      saveRDS(estimates, fileName)
    }
    allEstimates[[i]] <- estimates
  }
}
allEstimates <- do.call(rbind, allEstimates)
saveRDS(allEstimates,  "RealWorldExample/estimatesMatched.rds")


ParallelLogger::stopCluster(cluster)
