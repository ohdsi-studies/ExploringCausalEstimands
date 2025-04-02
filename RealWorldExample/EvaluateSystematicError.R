source("SetConnectionDetails.R")


hrEstimates <- CohortMethod::getResultsSummary(outputFolder)
EmpiricalCalibration::plotCalibrationEffect(hrEstimates$logRr,
                                            hrEstimates$seLogRr,
                                            showCis = TRUE,
                                            showExpectedAbsoluteSystematicError = TRUE,
                                            xLabel = "Hazard Ratio",
                                            fileName = "hrEstimates.png")



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

