library(DatabaseConnector)
options(andromedaTempFolder = "e:/andromedaTemp")

connectionDetails <- createConnectionDetails(
  dbms = "spark",
  connectionString = keyring::key_get("databricksConnectionString"),
  user = "token",
  password = keyring::key_get("databricksToken")
)
cdmDatabaseSchema <- "merative_mdcr.cdm_merative_mdcr_v3045"
cohortDatabaseSchema <- "scratch.scratch_mschuemi"
cohortTable  <- "hrs_and_causality"
options(sqlRenderTempEmulationSchema = "scratch.scratch_mschuemi")
outputFolder <- "e:/HRsAndCausality"
