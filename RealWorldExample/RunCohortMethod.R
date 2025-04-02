source("SetConnectionDetails.R")
library(CohortMethod)
library(survival)

# specify analyses -------------------------------------------------------------

negativeControlConceptIds <- c(72748,73241,73560,75911,76786,77965,78619,81151,81378,81634,133655,134438,136368,137951,139099,140641,140648,140842,141932,194083,195873,196168,199192,201606,259995,373478,374375,376707,377572,378427,380706,432303,432593,433111,433527,433577,434165,434203,434327,436409,437264,438130,438329,439790,440193,440329,441788,443172,444132,4012570,4012934,4083487,4088290,4091513,4092879,4092896,4103640,4103703,4115367,4115402,4149084,4166231,4170770,4180978,4201390,4201717,4202045,4209423,4213540,4231770,4344500,36713918,40480893,40481632,44783954,45757370,46269889,46286594)
# Setting outcomes to be outcome of interest so all data objects are generated for the negative controls:
negativeControlOutcomes <- lapply(negativeControlConceptIds,
                                  function(outcomeId) createOutcome(outcomeId = outcomeId,
                                                                    outcomeOfInterest = TRUE))
tcos <- createTargetComparatorOutcomes(targetId = 12676, # Thiazides
                                       comparatorId = 12672, # ACEis
                                       outcomes = negativeControlOutcomes)
targetComparatorOutcomesList <- list(tcos)

covarSettings <- createDefaultCovariateSettings(excludedCovariateConceptIds = c(1308216, 1310756, 1331235, 1334456, 1335471, 1340128, 1341927, 1342439, 1363749, 1373225, 907013, 974166, 978555, 1395058),
                                                addDescendantsToExclude = TRUE)

getDbCmDataArgs <- createGetDbCohortMethodDataArgs(washoutPeriod = 0,
                                                   restrictToCommonPeriod = FALSE,
                                                   firstExposureOnly = FALSE,
                                                   removeDuplicateSubjects = "keep first",
                                                   covariateSettings = covarSettings)

createStudyPopArgs <- createCreateStudyPopulationArgs(removeSubjectsWithPriorOutcome = TRUE,
                                                      minDaysAtRisk = 1,
                                                      riskWindowStart = 1,
                                                      startAnchor = "cohort start",
                                                      riskWindowEnd = 0,
                                                      endAnchor = "cohort end")

createPsArgs <- createCreatePsArgs(maxCohortSizeForFitting = 250000,
                                   control = createControl(cvType = "auto",
                                                           startingVariance = 0.01,
                                                           tolerance = 2e-07,
                                                           noiseLevel = "quiet",
                                                           cvRepetitions = 1),
                                   estimator = "att")

matchOnPsArgs <- createMatchOnPsArgs(maxRatio = 1)

computeSharedCovBalArgs <- createComputeCovariateBalanceArgs()

computeCovBalArgs <- createComputeCovariateBalanceArgs(covariateFilter = CohortMethod::getDefaultCmTable1Specifications())

fitOutcomeModelArgs <- createFitOutcomeModelArgs(modelType = "cox")

cmAnalysis <- createCmAnalysis(analysisId = 1,
                               description = "1-on-1 Matching",
                               getDbCohortMethodDataArgs = getDbCmDataArgs,
                               createStudyPopArgs = createStudyPopArgs,
                               createPsArgs = createPsArgs,
                               matchOnPsArgs = matchOnPsArgs,
                               computeSharedCovariateBalanceArgs = computeSharedCovBalArgs,
                               computeCovariateBalanceArgs = computeCovBalArgs,
                               fitOutcomeModelArgs = fitOutcomeModelArgs)


cmAnalysisList <- list(cmAnalysis)


# Run analyses -----------------------------------------------------------------
multiThreadingSettings <- createDefaultMultiThreadingSettings(parallel::detectCores())

dir.create(outputFolder)
result <- runCmAnalyses(
  connectionDetails = connectionDetails,
  cdmDatabaseSchema = cdmDatabaseSchema,
  exposureDatabaseSchema = cohortDatabaseSchema,
  exposureTable = cohortTable,
  outcomeDatabaseSchema = cohortDatabaseSchema,
  outcomeTable = cohortTable,
  outputFolder = outputFolder,
  cmAnalysisList = cmAnalysisList,
  targetComparatorOutcomesList = targetComparatorOutcomesList,
  multiThreadingSettings = multiThreadingSettings
)

# Compute risk ratios using weighted KM estimator ------------------------------
computeWeights <- function(data) {
  weights <- ifelse(data$treatment == 1,
                    mean(data$treatment == 1),
                    mean(data$treatment == 0) * data$propensityScore / (1 - data$propensityScore)
  )
  return(weights)
}

calc_rate_ratio <- function(data, indices) {
  require(survival)
  sampled_data <- data[indices, ]
  
  # Fit Kaplan-Meier for treated group
  idx <- sampled_data$treatment == 1
  kmTreated <- survfit(Surv(survivalTime, y) ~ 1,
                       data = sampled_data[idx, ],
                       weights = sampled_data$iptw[idx])
  
  # Fit Kaplan-Meier for control group
  idx <- sampled_data$treatment == 0
  kmControl <- survfit(Surv(survivalTime, y) ~ 1,
                       data = sampled_data[idx, ],
                       weights = sampled_data$iptw[idx])
  
  # Calculate risks at the specified time point
  timePoint <- 365
  riskTreated <- 1 - summary(kmTreated, time = timePoint)$surv
  riskControl <- 1 - summary(kmControl, time = timePoint)$surv
  
  # Calculate and return the rate ratio
  return(riskTreated / riskControl)
}
cluster <- snow::makeCluster(10, type = "SOCK")

ref <- getFileReference(outputFolder)
estimates <- list()
for (i in 1:nrow(ref)) {
  print(i)
  population <- readRDS(file.path(outputFolder, ref$psFile[i]))
  population$y <- population$outcomeCount > 0
  population$weight <- computeWeights(population)
  if (sum(population$y) == 0) {
    estimates[[length(estimates) + 1]] <- data.frame(
      outcomeId = ref$outcomeId[i],
      rr = as.numeric(NA),
      lb = as.numeric(NA),
      ub = as.numeric(NA),
      logRr = as.numeric(NA),
      seLogRr = as.numeric(NA)
    )
  } else {
    population$treatment <- as.factor(population$treatment)  # Ensure treatment is a factor
    population$iptw <- computeWeights(population)
    
    boot_results <- boot::boot(data = population,
                               statistic = calc_rate_ratio,
                               R = 1000,
                               parallel = "snow",
                               ncpus = 10,
                               cl = cluster)
    ci <- tryCatch({
      boot::boot.ci(boot_results, type = "perc")
    },
    error = function(e) {
      warning("Error computing confidence interval: ", e$message)
      return(list(percent = rep(NA, 5)))
    }
    )
    ci <- ci$percent[c(4, 5)]
    rateRatio <- boot_results$t0
    estimates[[length(estimates) + 1]] <- data.frame(
      outcomeId = ref$outcomeId[i],
      rr = rateRatio,
      lb = ci[1],
      ub = ci[2],
      logRr = log(rateRatio),
      seLogRr = (ci[2] - ci[1]) / (2 * qnorm(0.975))
    )
  }
}
estimates <- do.call(rbind, estimates)
saveRDS(estimates, file.path(outputFolder, "rrEstimates.rds"))
snow::stopCluster(cluster)

# Compute risk ratios using unweighted KM estimator in matched population ------
calc_rate_ratio <- function(data, indices) {
  require(survival)
  sampled_data <- data[indices, ]
  if (sum(sampled_data$y[sampled_data$treatment == 1]) == 0 |
      sum(sampled_data$y[sampled_data$treatment == 0]) == 0) {
    return(1)
  }
  
  # Fit Kaplan-Meier for treated group
  idx <- sampled_data$treatment == 1
  kmTreated <- survfit(Surv(survivalTime, y) ~ 1,
                       data = sampled_data[idx, ])
  
  # Fit Kaplan-Meier for control group
  idx <- sampled_data$treatment == 0
  kmControl <- survfit(Surv(survivalTime, y) ~ 1,
                       data = sampled_data[idx, ])
  
  timePoint <- 365
  riskTreated <- 1 - summary(kmTreated, time = timePoint)$surv
  riskControl <- 1 - summary(kmControl, time = timePoint)$surv
  rr <- riskTreated / riskControl
  return(rr)
}

cluster <- snow::makeCluster(10, type = "SOCK")

ref <- getFileReference(outputFolder)
estimates <- list()
for (i in 56:nrow(ref)) {
  print(i)
  population <- readRDS(file.path(outputFolder, ref$strataFile[i]))
  population$y <- population$outcomeCount > 0
  if (sum(population$y[population$treatment == 1]) == 0 |
      sum(population$y[population$treatment == 0]) == 0) {
    estimates[[length(estimates) + 1]] <- data.frame(
      outcomeId = ref$outcomeId[i],
      rr = as.numeric(NA),
      lb = as.numeric(NA),
      ub = as.numeric(NA),
      logRr = as.numeric(NA),
      seLogRr = as.numeric(NA)
    )
  } else {
    population$treatment <- as.factor(population$treatment)  # Ensure treatment is a factor
    
    boot_results <- boot::boot(data = population,
                               statistic = calc_rate_ratio,
                               R = 1000,
                               parallel = "snow",
                               ncpus = 10,
                               cl = cluster)
    ci <- boot::boot.ci(boot_results, type = "perc")
    ci <- ci$percent[c(4, 5)]
    rateRatio <- boot_results$t0
    estimates[[length(estimates) + 1]] <- data.frame(
      outcomeId = ref$outcomeId[i],
      rr = rateRatio,
      lb = ci[1],
      ub = ci[2],
      logRr = log(rateRatio),
      seLogRr = (ci[2] - ci[1]) / (2 * qnorm(0.975))
    )
  }
}
estimates <- do.call(rbind, estimates)
saveRDS(estimates, file.path(outputFolder, "rrEstimatesMatched.rds"))
snow::stopCluster(cluster)
