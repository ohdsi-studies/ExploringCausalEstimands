source("RealWorldExample/SetConnectionDetails.R")
library(CohortMethod)
library(survival)

# specify analyses -------------------------------------------------------------

tcs <- readr::read_csv("RealWorldExample/TCs.csv", show_col_types = FALSE)
negativeControls <- readr::read_csv("RealWorldExample/NegativeControls.csv", show_col_types = FALSE)

targetComparatorOutcomesList <- list()
for (i in seq_len(nrow(tcs))) {
  negativeControlConceptIds <- negativeControls |>
    filter(indicationId == tcs$indicationId[i]) |>
    pull(conceptId)
  negativeControlOutcomes <- lapply(negativeControlConceptIds,
                                    function(outcomeId) createOutcome(outcomeId = outcomeId,
                                                                      outcomeOfInterest = TRUE))
  targetComparatorOutcomesList[[i]] <- createTargetComparatorOutcomes(
    targetId = tcs$targetId[i],
    comparatorId = tcs$comparatorId[i],
    outcomes = negativeControlOutcomes,
    excludedCovariateConceptIds = c(tcs$targetId[i], tcs$comparatorId[i])
  )
}

covarSettings <- createDefaultCovariateSettings(addDescendantsToExclude = TRUE)

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

# Save estimates ---------------------------------------------------------------
hrEstimates <- getResultsSummary(outputFolder) |>
  select(targetId, comparatorId, outcomeId, logRr, seLogRr, ci95Lb, ci95Ub, p, mdrr)
saveRDS(hrEstimates, "RealWorldExample/hrEstimates.rds")

# Compute follup-time ----------------------------------------------------------
ref <- getFileReference(outputFolder)
studyPop <- readRDS(file.path(outputFolder, ref$strataFile[1]))
quantile(studyPop$timeAtRisk, c(0.25, 0.5, 0.75, 0.9, 0.95, 1))
# 25%  50%  75% 100% 
#   52  132  386 6622 
