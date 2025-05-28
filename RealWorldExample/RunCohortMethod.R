source("RealWorldExample/SetConnectionDetails.R")
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

# Compute follup-time ----------------------------------------------------------
ref <- getFileReference(outputFolder)
studyPop <- readRDS(file.path(outputFolder, ref$strataFile[1]))
quantile(studyPop$timeAtRisk, c(0.25, 0.5, 0.75, 0.9, 0.95, 1))
# 25%  50%  75% 100% 
#   52  132  386 6622 
