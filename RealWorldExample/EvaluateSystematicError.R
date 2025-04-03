source("RealWorldExample/SetConnectionDetails.R")


hrEstimates <- CohortMethod::getResultsSummary(outputFolder)
EmpiricalCalibration::plotCalibrationEffect(hrEstimates$logRr,
                                            hrEstimates$seLogRr,
                                            showCis = TRUE,
                                            showExpectedAbsoluteSystematicError = TRUE,
                                            xLabel = "Hazard Ratio",
                                            fileName = "hrEstimates.png")

writeLines(sprintf("Type 1 error: %0.1f%%", 100*mean(hrEstimates$p < 0.05, na.rm = TRUE)))
# Type 1 error: 14.0%

rrEstimates <- readRDS(file.path(outputFolder, "rrEstimates.rds"))
EmpiricalCalibration::plotCalibrationEffect(rrEstimates$logRr,
                                            rrEstimates$seLogRr,
                                            showCis = TRUE,
                                            showExpectedAbsoluteSystematicError = TRUE,
                                            xLabel = "Relative Risk (first 365 days)",
                                            fileName = "rrEstimates.png")


rrEstimates <- readRDS(file.path(outputFolder, "rrEstimatesMatched.rds"))
EmpiricalCalibration::plotCalibrationEffect(rrEstimates$logRr,
                                            rrEstimates$seLogRr,
                                            showCis = TRUE,
                                            showExpectedAbsoluteSystematicError = TRUE,
                                            xLabel = "Relative Risk (first 365 days)",
                                            fileName = "rrEstimatesMatched.png")
writeLines(sprintf("Type 1 error: %0.1f%%", 100*mean(rrEstimates$lb > 1 | rrEstimates$ub < 1, na.rm = TRUE)))
# Type 1 error: 16.1%
