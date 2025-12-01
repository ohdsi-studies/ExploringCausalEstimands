# Script for creating the exposure and outcome cohorts

source("RealWorldExample/SetConnectionDetails.R")
library(Capr)
library(dplyr)
library(CirceR)

tcos <- readr::read_csv("RealWorldExample/TCOs.csv", show_col_types = FALSE)
negativeControls <- readr::read_csv("RealWorldExample/NegativeControls.csv", show_col_types = FALSE)

connection <- connect(connectionDetails)

# Exposure cohorts -------------------------------------------------------------
indicationConceptSets <- list()
mdd <- cs(
  descendants(4191716, 4212469, 4175329, 440383, 40546087), 
  descendants(exclude(377527, 379784, 433440, 435520, 436665, 438727, 442306, 443864, 4224940, 4239471, 36684319, 40481798)),
  name = "Major depressive disorder"
)
mdd <- getConceptSetDetails(mdd, connection, cdmDatabaseSchema)
indicationConceptSets[["Depression"]] <- mdd

hypertensiveDisorder <- cs(
  descendants(316866),
  name = "Hypertensive disorder"
)
hypertensiveDisorder <- getConceptSetDetails(hypertensiveDisorder, connection, cdmDatabaseSchema)
indicationConceptSets[["Hypertension"]] <- hypertensiveDisorder

exposures <- bind_rows(
  tcos |>
    distinct(conceptId = targetId, name = targetName, indicationId),
  tcos |>
    distinct(conceptId = comparatorId, name = comparatorName, indicationId)
)

cohortDefinitionSet  <- list()
for (i in seq_len(nrow(exposures))) {
  drugConceptSet <- cs(
    descendants(exposures$conceptId[i]),
    name = exposures$name[i]
  )
  drugConceptSet <- getConceptSetDetails(drugConceptSet, connection, cdmDatabaseSchema)
  
  indicationConceptSet <- indicationConceptSets[[exposures$indicationId[i]]]
  
  newUserCohort <- cohort(
    entry = entry(
      drugExposure(drugConceptSet, firstOccurrence()),
      observationWindow = continuousObservation(priorDays = 365)
    ),
    attrition = attrition(
      "Has indication" = withAll(
        atLeast(1, conditionOccurrence(indicationConceptSet), duringInterval(eventStarts(-365, 0)))
      )
    ),
    exit = exit(endStrategy = drugExit(drugConceptSet, persistenceWindow = 30, surveillanceWindow = 0))
  )
  json <- as.json(newUserCohort)
  cohortDefinitionSet [[i]] <- tibble(
    cohortId = exposures$conceptId[i],
    cohortName = exposures$name[i],
    json = json,
    sql = buildCohortQuery(json, createGenerateOptions(generateStats = FALSE))
  )
}
cohortDefinitionSet <- bind_rows(cohortDefinitionSet)

# Outcome cohorts --------------------------------------------------------------
suicideAndCSuicidalIdeation <- cs(descendants(435446,439235,440925,444362,4009713,4021339,4092411,4181216,4216115,4219484,4273391,4303690),
                                  name = "Suicide and suicidal ideation")
suicideAndCSuicidalIdeation <- getConceptSetDetails(suicideAndCSuicidalIdeation, connection, cdmDatabaseSchema)
suicideAndCSuicidalIdeationCohort <- cohort(
  entry = entry(
    conditionOccurrence(suicideAndCSuicidalIdeation),
    observation(suicideAndCSuicidalIdeation),
    primaryCriteriaLimit = "First"
  )
)

hyponatremia <- cs(descendants(435515,4232311),
                   name = "Hyponatremia")
hyponatremia <- getConceptSetDetails(hyponatremia, connection, cdmDatabaseSchema)
hyponatremiaCohort <- cohort(
  entry = entry(
    conditionOccurrence(hyponatremia),
    primaryCriteriaLimit = "First"
  )
)

vomitingOrNaussea <- cs(descendants(4012500, 31967, 4101344, 201218),
                        exclude(descendants(26727, 30284, 436166, 440785, 4153517, 4156946, 4274029, 4323686, 40480291, 45757414)),
                                  name = "Vomiting or naussea")
vomitingOrNaussea <- getConceptSetDetails(vomitingOrNaussea, connection, cdmDatabaseSchema)
vomitingOrNausseaCohort <- cohort(
  entry = entry(
    conditionOccurrence(vomitingOrNaussea),
    observation(vomitingOrNaussea),
    primaryCriteriaLimit = "First"
  )
)

orgasmDisorder <- cs(descendants(4221288),
                   name = "Orgasm disorder")
orgasmDisorder <- getConceptSetDetails(orgasmDisorder, connection, cdmDatabaseSchema)
orgasmDisorderCohort <- cohort(
  entry = entry(
    conditionOccurrence(orgasmDisorder),
    primaryCriteriaLimit = "First"
  )
)


angioedema <- cs(descendants(432791,4296370),
                 name = "Angioedema")
angioedema <- getConceptSetDetails(angioedema, connection, cdmDatabaseSchema)
angioedemaOrUrticaria <- cs(descendants(432791,4270861,4296370,4297361),
                           name = "Angioedema or urticaria") 
angioedemaOrUrticaria <- getConceptSetDetails(angioedemaOrUrticaria, connection, cdmDatabaseSchema)
angioedemaCohort <- cohort(
  entry = entry(
    conditionOccurrence(angioedemaOrUrticaria),
    primaryCriteriaLimit = "All",
    qualifiedLimit = "First"
  ),
  attrition = attrition(
    "Angioedema within 3 days" = withAll(
      atLeast(1, conditionOccurrence(angioedema), duringInterval(eventStarts(0,3)))
    )
  )
)
json <- c(as.json(suicideAndCSuicidalIdeationCohort), 
          as.json(angioedemaCohort), 
          as.json(hyponatremiaCohort), 
          as.json(vomitingOrNausseaCohort), 
          as.json(orgasmDisorderCohort))
outcomeCohortDefinitionsSet <- tibble(
  cohortId = c(1, 2, 3, 4, 5),
  cohortName = c("Suicide and suicidal ideation", 
                 "Angioedema",
                 "Hyponatremia",
                 "Vomiting or naussea",
                 "Orgams disorder"),
  json = json,
  sql = c(buildCohortQuery(json[1], createGenerateOptions(generateStats = FALSE)),
          buildCohortQuery(json[2], createGenerateOptions(generateStats = FALSE)),
          buildCohortQuery(json[3], createGenerateOptions(generateStats = FALSE)),
          buildCohortQuery(json[4], createGenerateOptions(generateStats = FALSE)),
          buildCohortQuery(json[5], createGenerateOptions(generateStats = FALSE)))
)
cohortDefinitionSet <- bind_rows(cohortDefinitionSet, outcomeCohortDefinitionsSet)
# json <- as.json(angioedemaCohort)
# writeLines(CirceR::cohortPrintFriendly(json))


# Generate cohorts -------------------------------------------------------------
cohortTableNames <- CohortGenerator::getCohortTableNames(cohortTable)
CohortGenerator::createCohortTables(
  connection = connection, 
  cohortDatabaseSchema = cohortDatabaseSchema,
  cohortTableNames = cohortTableNames
)

CohortGenerator::generateCohortSet(
  connection = connection, 
  cohortDatabaseSchema = cohortDatabaseSchema,
  cohortTableNames = cohortTableNames,
  cdmDatabaseSchema = cdmDatabaseSchema,
  cohortDefinitionSet = cohortDefinitionSet
)

negativeControlOutcomeCohortSet <- negativeControls |>
  mutate(cohortId = conceptId) |>
  select("cohortId", cohortName = "outcomeName", outcomeConceptId = "conceptId") |>
  filter(!duplicated(cohortId))
CohortGenerator::generateNegativeControlOutcomeCohorts(
  connection = connection, 
  cohortDatabaseSchema = cohortDatabaseSchema,
  cohortTable = cohortTable,
  cdmDatabaseSchema = cdmDatabaseSchema,
  negativeControlOutcomeCohortSet = negativeControlOutcomeCohortSet,
  occurrenceType = "first",
  detectOnDescendants = TRUE
)

# Count cohorts ----------------------------------------------------------------
cohortCounts <- CohortGenerator::getCohortCounts(
  connection = connection, 
  cohortDatabaseSchema = cohortDatabaseSchema,
  cohortTable = cohortTable
)
cohortCounts <- cohortCounts |>
  inner_join(
    bind_rows(
      exposures |>
        select(cohortId = "conceptId", cohortName = "name") |>
        mutate(type = "exposure"),
      tcos |>
        select(cohortId = "outcomeId", cohortName = "outcomeName") |>
        mutate(type = "positive control outcome"),
      negativeControlOutcomeCohortSet |>
        select("cohortId", "cohortName") |>
        mutate(type = "negative control outcome")
    ),
    by = join_by("cohortId")
  )
if (!dir.exists(outputFolder)) {
  dir.create(outputFolder, recursive = TRUE)
}
readr::write_csv(cohortCounts, file.path(outputFolder, "cohortCounts.csv"))

disconnect(connection)
