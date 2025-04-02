source("SetConnectionDetails.R")

# Cohort definitions have already been saved, so commenting out fetching them
# from WebAPI:
ROhdsiWebApi::authorizeWebApi(
  baseUrl = Sys.getenv("baseUrl"),
  authMethod = "windows"
)
cohortDefinitionSet <- ROhdsiWebApi::exportCohortDefinitionSet(
  baseUrl = Sys.getenv("baseUrl"),
  generateStats = FALSE,
  cohortIds = c(12676, 12672) # Thiazides, ACEi
)
saveRDS(cohortDefinitionSet, "cohortDefinitionSet.rds")
cohortDefinitionSet = readRDS("cohortDefinitionSet.rds")

negativeControlConceptIds <- c(72748,73241,73560,75911,76786,77965,78619,81151,81378,81634,133655,134438,136368,137951,139099,140641,140648,140842,141932,194083,195873,196168,199192,201606,259995,373478,374375,376707,377572,378427,380706,432303,432593,433111,433527,433577,434165,434203,434327,436409,437264,438130,438329,439790,440193,440329,441788,443172,444132,4012570,4012934,4083487,4088290,4091513,4092879,4092896,4103640,4103703,4115367,4115402,4149084,4166231,4170770,4180978,4201390,4201717,4202045,4209423,4213540,4231770,4344500,36713918,40480893,40481632,44783954,45757370,46269889,46286594)

connection <- connect(connectionDetails)

cohortTableNames <- CohortGenerator::getCohortTableNames(cohortTable = cohortTable)
CohortGenerator::createCohortTables(
  connection = connection,
  cohortDatabaseSchema = cohortDatabaseSchema,
  cohortTableNames = cohortTableNames,
)
CohortGenerator::generateCohortSet(
  connection = connection,
  cdmDatabaseSchema = cdmDatabaseSchema,
  cohortDatabaseSchema = cohortDatabaseSchema,
  cohortTableNames = cohortTableNames,
  cohortDefinitionSet = cohortDefinitionSet
)
CohortGenerator::generateNegativeControlOutcomeCohorts(
  connection = connection,
  cdmDatabaseSchema = cdmDatabaseSchema,
  cohortDatabaseSchema = cohortDatabaseSchema,
  cohortTable = cohortTable,
  occurrenceType = "first",
  detectOnDescendants = TRUE,
  negativeControlOutcomeCohortSet = data.frame(
    cohortId = negativeControlConceptIds,
    cohortName = negativeControlConceptIds,
    outcomeConceptId = negativeControlConceptIds
  )
)
disconnect(connection)
