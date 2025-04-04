source("RealWorldExample/SetConnectionDetails.R")
library(dplyr)


computeMeanPrecision <- function(seLogRr) {
  # idx <- is.na(seLogRr) | is.infinite(seLogRr) | seLogRr == 0
  # seLogRr[idx] <- 999
  precision <- 1 / seLogRr^2
  return(exp(mean(log(precision))))
}

hrEstimates <- CohortMethod::getResultsSummary(outputFolder) |>
  arrange(outcomeId) |>
  mutate(valid = !is.na(seLogRr) & !is.infinite(seLogRr) & seLogRr != 0) |>
  mutate(seLogRr = if_else(valid, seLogRr, 999))
rrEstimates <- readRDS(file.path(outputFolder, "rrEstimatesMatched.rds")) |>
  arrange(outcomeId) |>
  mutate(valid = !is.na(seLogRr) & !is.infinite(seLogRr) & seLogRr != 0) |>
  mutate(seLogRr = if_else(valid, seLogRr, 999))

mean(hrEstimates$valid)
mean(rrEstimates$valid)
sum(hrEstimates$valid)
sum(rrEstimates$valid)

idxBothValid <- hrEstimates$valid & rrEstimates$valid
mean(idxBothValid)

hrEstimates <- hrEstimates[idxBothValid, ]
rrEstimates <- rrEstimates[idxBothValid, ]

null <- EmpiricalCalibration::fitMcmcNull(hrEstimates$logRr, hrEstimates$seLogRr)
EmpiricalCalibration::plotCalibrationEffect(hrEstimates$logRr,
                                            hrEstimates$seLogRr,
                                            showCis = TRUE,
                                            showExpectedAbsoluteSystematicError = TRUE,
                                            null = null,
                                            xLabel = "Hazard Ratio",
                                            fileName = "hrEstimates.png")
ease <- EmpiricalCalibration::computeExpectedAbsoluteSystematicError(null)
writeLines(sprintf("Type 1 error: %0.1f%%. Mean precision: %0.2f. EASE: %0.2f (%0.2f - %0.2f)", 
                   100*mean(hrEstimates$p < 0.05, na.rm = TRUE),
                   computeMeanPrecision(hrEstimates$seLogRr),
                   ease$ease, 
                   ease$ciLb, 
                   ease$ciUb))
# Type 1 error: 15.1%. Mean precision: 41.91. EASE: 0.07 (0.05 - 0.11)


null <- EmpiricalCalibration::fitMcmcNull(rrEstimates$logRr, hrEstimates$seLogRr)
EmpiricalCalibration::plotCalibrationEffect(rrEstimates$logRr,
                                            rrEstimates$seLogRr,
                                            showCis = TRUE,
                                            showExpectedAbsoluteSystematicError = TRUE,
                                            null = null,
                                            xLabel = "Relative Risk (first 365 days)",
                                            fileName = "rrEstimatesMatched.png")
ease <- EmpiricalCalibration::computeExpectedAbsoluteSystematicError(null)
writeLines(sprintf("Type 1 error: %0.1f%%. Mean precision: %0.2f. EASE: %0.2f (%0.2f - %0.2f)", 
                   100*mean(rrEstimates$lb > 1 | rrEstimates$ub < 1, na.rm = TRUE),
                   computeMeanPrecision(rrEstimates$seLogRr),
                   ease$ease, 
                   ease$ciLb, 
                   ease$ciUb))
# Type 1 error: 15.1%. Mean precision: 16.89. EASE: 0.07 (0.05 - 0.11)


# rrEstimates <- readRDS(file.path(outputFolder, "rrEstimates.rds"))
# EmpiricalCalibration::plotCalibrationEffect(rrEstimates$logRr,
#                                             rrEstimates$seLogRr,
#                                             showCis = TRUE,
#                                             showExpectedAbsoluteSystematicError = TRUE,
#                                             xLabel = "Relative Risk (first 365 days)",
#                                             fileName = "rrEstimates.png")