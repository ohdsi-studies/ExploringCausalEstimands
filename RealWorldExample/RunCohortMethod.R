# Script for running the CohortMethod package. Assumes CreateCohorts.R has been executed.

source("RealWorldExample/SetConnectionDetails.R")
library(CohortMethod)
library(survival)
library(dplyr)

# Specify analyses -------------------------------------------------------------
tcos <- readr::read_csv("RealWorldExample/TCOs.csv", show_col_types = FALSE)
negativeControls <- readr::read_csv("RealWorldExample/NegativeControls.csv", show_col_types = FALSE)

tcs <- tcos |>
  group_by(targetId, comparatorId) |>
  group_split()
targetComparatorOutcomesList <- list()
for (i in seq_along(tcs)) {
  tc <- tcs[[i]]
  negativeControlConceptIds <- negativeControls |>
    filter(indicationId == tc$indicationId[1]) |>
    pull(conceptId)
  outcomes <- lapply(c(negativeControlConceptIds, tc$outcomeId),
                     function(outcomeId) createOutcome(outcomeId = outcomeId,
                                                       outcomeOfInterest = TRUE))
  targetComparatorOutcomesList[[i]] <- createTargetComparatorOutcomes(
    targetId = tc$targetId[1],
    comparatorId = tc$comparatorId[1],
    outcomes = outcomes,
    excludedCovariateConceptIds = c(tc$targetId[1], tc$comparatorId[1])
  )
}

covarSettings <- createDefaultCovariateSettings(addDescendantsToExclude = TRUE)

getDbCmDataArgs <- createGetDbCohortMethodDataArgs(washoutPeriod = 365,
                                                   restrictToCommonPeriod = TRUE,
                                                   firstExposureOnly = TRUE,
                                                   removeDuplicateSubjects = "keep first",
                                                   covariateSettings = covarSettings)

createStudyPopArgs <- createCreateStudyPopulationArgs(removeSubjectsWithPriorOutcome = TRUE,
                                                      minDaysAtRisk = 1,
                                                      riskWindowStart = 0,
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

# Compute counts ---------------------------------------------------------------
ref <- getFileReference(outputFolder)
refRows <- ref |>
  filter(outcomeId %in% c(1, 2))
for (i in 1:2) {
  if (refRows$outcomeId[i] == 1) {
    print("Suicide and suicidal ideation")
    targetLabel <- "Sertraline"
    comparatorLabel <- "Bupropion"
  } else {
    print("Angioedema")
    targetLabel <- "Lisinopril"
    comparatorLabel <- "Hydrochlorothiazide"
    
  }
  strataPop <- readRDS(file.path(outputFolder, refRows$strataFile[i]))
  attrition <- getAttritionTable(strataPop)
  print(attrition[2, ])
  print(attrition[5, ])
  print(sum(strataPop$outcomeCount > 0))
  print(mean(strataPop$outcomeCount > 0) * 100)
  
  plotKaplanMeier(population = strataPop,
                  targetLabel = targetLabel, 
                  comparatorLabel = comparatorLabel,
                  fileName = file.path("RealWorldExample", sprintf("KM_o%d.png", refRows$outcomeId[i])))
  plotKaplanMeier(population = strataPop,
                  targetLabel = targetLabel, 
                  comparatorLabel = comparatorLabel,
                  fileName = file.path("RealWorldExample", sprintf("KM_o%d.svg", refRows$outcomeId[i])))
  strataPop <- strataPop |>
    mutate(y = outcomeCount > 0)
  fit <- coxph(Surv(survivalTime, y) ~ treatment, data = strataPop)
  phTest <- cox.zph(fit)
  print(phTest)
  
  balance <- readRDS(file.path(outputFolder, refRows$sharedBalanceFile[i]))
  writeLines(sprintf("Number of covariates: %d, max ASDM: %0.2f", nrow(balance), max(abs(balance$afterMatchingStdDiff), na.rm = TRUE)))
}
