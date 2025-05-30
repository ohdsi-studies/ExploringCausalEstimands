library(dplyr)

mdrrTheshold <- 2

hrEstimates <- readRDS("RealWorldExample/hrEstimates.rds")
nonHrEstimates <- readRDS("RealWorldExample/nonHrEstimatesMatching.rds")

# Restrict to TCOs having sufficient power -------------------------------------
valid <- hrEstimates |>
  filter(mdrr < mdrrTheshold) |>
  select("targetId", "comparatorId", "outcomeId")

hrEstimates <- hrEstimates |>
  inner_join(valid, by = join_by("targetId", "comparatorId", "outcomeId"))
nonHrEstimates <- nonHrEstimates |>
  inner_join(valid, by = join_by("targetId", "comparatorId", "outcomeId"))

# Pivot:
nonHrEstimates <- bind_rows(
    nonHrEstimates |>
      select(targetId, comparatorId, outcomeId, timePoint, estimate = rr, lb = lbRrPercentile, ub = ubRrPercentile, se = seLogRrPercentile) |>
      mutate(estimand = "rr", method = "percentile"),
    nonHrEstimates |>
      select(targetId, comparatorId, outcomeId, timePoint, estimate = rd, lb = lbRdPercentile, ub = ubRdPercentile, se = seRdPercentile) |>
      mutate(estimand = "rd", method = "percentile"),
    nonHrEstimates |>
      select(targetId, comparatorId, outcomeId, timePoint, estimate = rr, lb = lbRrAsymptotics, ub = ubRrAsymptotics, se = seLogRrNormal) |>
      mutate(estimand = "rr", method = "asymptotics"),
    nonHrEstimates |>
      select(targetId, comparatorId, outcomeId, timePoint, estimate = rd, lb = lbRdAsymptotics, ub = ubRdAsymptotics, se = seRdNormal) |>
      mutate(estimand = "rd", method = "asymptotics")
  )

# Compute metrics --------------------------------------------------------------

computeType1Error <- function(lb, ub, h0 = 1) {
  rejectNull <- (!is.na(lb) & h0 < lb) | (!is.na(ub) & h0 > ub)
  return(mean(rejectNull))
}

computeMeanPrecision <- function(seLogRr) {
  seLogRr[is.na(seLogRr)] <- 999
  precision <- 1 / seLogRr^2
  return(exp(mean(log(precision))))
}

computeEase <- function(logRr, seLogRr) {
  null <- EmpiricalCalibration::fitMcmcNull(logRr, seLogRr)
  ease <- EmpiricalCalibration::computeExpectedAbsoluteSystematicError(null)
  return(ease$ease)
}

hrEstimates |>
  summarise(
    count = n(),
    type1Error = computeType1Error(ci95Lb, ci95Ub, h0 = 1),
    meanPrecision = computeMeanPrecision(seLogRr),
    ease = computeEase(logRr, seLogRr)
  )

nonHrEstimates |>
  group_by(timePoint, estimand, method) |>
  summarise(
    count = n(),
    type1Error = computeType1Error(lb, ub, h0 = if_else(estimand == "rr", 1, 0)),
    meanPrecision = computeMeanPrecision(se),
    ease = computeEase(if_else(estimand == "rr", log(estimate), estimate), se),
    .groups = "drop"
  )


