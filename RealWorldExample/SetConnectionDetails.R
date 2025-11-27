library(DatabaseConnector)
options(andromedaTempFolder = "e:/andromedaTemp")

connectionDetails <- createConnectionDetails(
  dbms = "spark",
  connectionString = keyring::key_get("databricksConnectionString"),
  user = "token",
  password = keyring::key_get("databricksToken")
)
cdmDatabaseSchema <- "optum_extended_dod.cdm_optum_extended_dod_v3492"
cohortDatabaseSchema <- "scratch.scratch_mschuemi"
cohortTable  <- "hrs_and_causality"
options(sqlRenderTempEmulationSchema = "scratch.scratch_mschuemi")
outputFolder <- "e:/HRsAndCausality"
maxCores <- 13
