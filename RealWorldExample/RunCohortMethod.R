source("RealWorldExample/SetConnectionDetails.R")
library(CohortMethod)
library(survival)

# specify analyses -------------------------------------------------------------

tcos <- readr::read_csv("RealWorldExample/TCOs.csv", show_col_types = FALSE)
negativeControls <- readr::read_csv("RealWorldExample/NegativeControls.csv", show_col_types = FALSE)

targetComparatorOutcomesList <- list()
for (i in seq_len(nrow(tcos))) {
  negativeControlConceptIds <- negativeControls |>
    filter(indicationId == tcos$indicationId[i]) |>
    pull(conceptId)
  outcomes <- lapply(c(negativeControlConceptIds, tcos$outcomeId[i]),
                     function(outcomeId) createOutcome(outcomeId = outcomeId,
                                                       outcomeOfInterest = TRUE))
  targetComparatorOutcomesList[[i]] <- createTargetComparatorOutcomes(
    targetId = tcos$targetId[i],
    comparatorId = tcos$comparatorId[i],
    outcomes = outcomes,
    excludedCovariateConceptIds = c(tcos$targetId[i], tcos$comparatorId[i])
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

