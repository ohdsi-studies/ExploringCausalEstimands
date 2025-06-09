library(dplyr)

mdrrTheshold <- 2

hrEstimates <- readRDS("RealWorldExample/hrEstimates.rds")
nonHrEstimates <- readRDS("RealWorldExample/nonHrEstimatesMatching.rds")
nonHrEstimatesWeighted <- readRDS("RealWorldExample/rrEstimatesWeighted.rds")

# Restrict to TCOs having sufficient power -------------------------------------
valid <- hrEstimates |>
  filter(mdrr < mdrrTheshold) |>
  select("targetId", "comparatorId", "outcomeId")

hrEstimates <- hrEstimates |>
  inner_join(valid, by = join_by("targetId", "comparatorId", "outcomeId"))
nonHrEstimates <- nonHrEstimates |>
  inner_join(valid, by = join_by("targetId", "comparatorId", "outcomeId"))
nonHrEstimatesWeighted <- nonHrEstimatesWeighted |>
  inner_join(valid, by = join_by("targetId", "comparatorId", "outcomeId"))

# Compute metrics --------------------------------------------------------------

pivotNonHr <- function(nonHrEstimates) {
  nonHrEstimates <- bind_rows(
    nonHrEstimates |>
      select(targetId, comparatorId, outcomeId, timePoint, estimate = rr, lb = lbRrPercentile, ub = ubRrPercentile, se = seLogRrPercentile) |>
      mutate(estimand = "rr", method = "percentile"),
    nonHrEstimates |>
      select(targetId, comparatorId, outcomeId, timePoint, estimate = rd, lb = lbRdPercentile, ub = ubRdPercentile, se = seRdPercentile) |>
      mutate(estimand = "rd", method = "percentile"),
    nonHrEstimates |>
      select(targetId, comparatorId, outcomeId, timePoint, estimate = rr, lb = lbRrAsymptotics, ub = ubRrAsymptotics, se = seLogRrAsymptotics) |>
      mutate(estimand = "rr", method = "asymptotics"),
    nonHrEstimates |>
      select(targetId, comparatorId, outcomeId, timePoint, estimate = rd, lb = lbRdAsymptotics, ub = ubRdAsymptotics, se = seRdAsymptotics) |>
      mutate(estimand = "rd", method = "asymptotics")
  )
  return(nonHrEstimates)
}

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
  group_by(targetId, comparatorId) |>
  summarise(
    count = n(),
    type1Error = computeType1Error(ci95Lb, ci95Ub, h0 = 1),
    meanPrecision = computeMeanPrecision(seLogRr),
    ease = computeEase(logRr, seLogRr)
  )
# # A tibble: 3 × 6
# # Groups:   targetId [3]
# targetId comparatorId count type1Error meanPrecision   ease
# <dbl>        <dbl> <int>      <dbl>         <dbl>  <dbl>
# 1   715259       743670    27      0.111          84.3 0.0681
# 2  1334456     40226742    37      0.162         110.  0.0673
# 3  1559684      1560171    37      0.135         106.  0.0499

nonHrEstimates <- nonHrEstimates |>
  pivotNonHr() 

results <- nonHrEstimates |>
  group_by(targetId, comparatorId, timePoint, estimand, method) |>
  summarise(
    count = n(),
    type1Error = computeType1Error(lb, ub, h0 = if_else(estimand == "rr", 1, 0)),
    meanPrecision = computeMeanPrecision(se),
    ease = computeEase(if_else(estimand == "rr", log(estimate), estimate), se),
    .groups = "drop"
  ) 
print(results, n = 60)
# # A tibble: 60 × 9
# targetId comparatorId timePoint estimand method      count type1Error meanPrecision    ease
# <dbl>        <dbl>     <dbl> <chr>    <chr>       <int>      <dbl>         <dbl>   <dbl>
# 1   715259       743670       180 rd       asymptotics    27     0.148       990866.  0.00281
# 2   715259       743670       180 rd       percentile     27     0.148       994606.  0.00281
# 3   715259       743670       180 rr       asymptotics    27     0.148           37.8 0.0857 
# 4   715259       743670       180 rr       percentile     27     0.148           37.8 0.0862 
# 5   715259       743670       365 rd       asymptotics    27     0.0370      343309.  0.00304
# 6   715259       743670       365 rd       percentile     27     0.0741      347694.  0.00304
# 7   715259       743670       365 rr       asymptotics    27     0.0370          49.3 0.0467 
# 8   715259       743670       365 rr       percentile     27     0.0741          49.7 0.0432 
# 9   715259       743670       730 rd       asymptotics    27     0.0741       97233.  0.00371
# 10   715259       743670       730 rd       percentile     27     0.0741       97771.  0.00371
# 11   715259       743670       730 rr       asymptotics    27     0.0741          47.6 0.0530 
# 12   715259       743670       730 rr       percentile     27     0.0741          48.2 0.0544 
# 13   715259       743670      1095 rd       asymptotics    27     0            36118.  0.00425
# 14   715259       743670      1095 rd       percentile     27     0            36433.  0.00426
# 15   715259       743670      1095 rr       asymptotics    27     0               36.8 0.0342 
# 16   715259       743670      1095 rr       percentile     27     0               37.2 0.0354 
# 17   715259       743670      1460 rd       asymptotics    27     0            15307.  0.00485
# 18   715259       743670      1460 rd       percentile     27     0            15357.  0.00487
# 19   715259       743670      1460 rr       asymptotics    27     0               26.0 0.0418 
# 20   715259       743670      1460 rr       percentile     27     0               26.2 0.0413 
# 21  1334456     40226742       180 rd       asymptotics    37     0.0541     3144357.  0.00208
# 22  1334456     40226742       180 rd       percentile     37     0.0541     3162352.  0.00208
# 23  1334456     40226742       180 rr       asymptotics    37     0.0541          39.1 0.0422 
# 24  1334456     40226742       180 rr       percentile     37     0.0541          39.1 0.0420 
# 25  1334456     40226742       365 rd       asymptotics    37     0.108      1161704.  0.00225
# 26  1334456     40226742       365 rd       percentile     37     0.108      1180981.  0.00224
# 27  1334456     40226742       365 rr       asymptotics    37     0.0811          56.3 0.0317 
# 28  1334456     40226742       365 rr       percentile     37     0.108           57.6 0.0338 
# 29  1334456     40226742       730 rd       asymptotics    37     0.216       355395.  0.00273
# 30  1334456     40226742       730 rd       percentile     37     0.162       360027.  0.00273
# 31  1334456     40226742       730 rr       asymptotics    37     0.216           66.0 0.0476 
# 32  1334456     40226742       730 rr       percentile     37     0.162           67.3 0.0520 
# 33  1334456     40226742      1095 rd       asymptotics    37     0.135       147225.  0.00333
# 34  1334456     40226742      1095 rd       percentile     37     0.135       149821.  0.00333
# 35  1334456     40226742      1095 rr       asymptotics    37     0.135           60.2 0.0675 
# 36  1334456     40226742      1095 rr       percentile     37     0.135           60.8 0.0693 
# 37  1334456     40226742      1460 rd       asymptotics    37     0.135        72489.  0.00378
# 38  1334456     40226742      1460 rd       percentile     37     0.135        72947.  0.00378
# 39  1334456     40226742      1460 rr       asymptotics    37     0.135           50.9 0.0630 
# 40  1334456     40226742      1460 rr       percentile     37     0.135           51.0 0.0625 
# 41  1559684      1560171       180 rd       asymptotics    37     0.0270     8834968.  0.00197
# 42  1559684      1560171       180 rd       percentile     37     0.0270     8899301.  0.00197
# 43  1559684      1560171       180 rr       asymptotics    37     0.0270          38.3 0.0367 
# 44  1559684      1560171       180 rr       percentile     37     0.0270          38.3 0.0369 
# 45  1559684      1560171       365 rd       asymptotics    37     0.135      3318916.  0.00210
# 46  1559684      1560171       365 rd       percentile     37     0.135      3366554.  0.00210
# 47  1559684      1560171       365 rr       asymptotics    37     0.135           53.9 0.0437 
# 48  1559684      1560171       365 rr       percentile     37     0.135           54.7 0.0437 
# 49  1559684      1560171       730 rd       asymptotics    37     0.135      1024562.  0.00240
# 50  1559684      1560171       730 rd       percentile     37     0.135      1054871.  0.00239
# 51  1559684      1560171       730 rr       asymptotics    37     0.135           61.0 0.0758 
# 52  1559684      1560171       730 rr       percentile     37     0.135           62.7 0.0759 
# 53  1559684      1560171      1095 rd       asymptotics    37     0.0270      423642.  0.00265
# 54  1559684      1560171      1095 rd       percentile     37     0.0270      427352.  0.00264
# 55  1559684      1560171      1095 rr       asymptotics    37     0.0270          54.4 0.0607 
# 56  1559684      1560171      1095 rr       percentile     37     0.0270          54.5 0.0599 
# 57  1559684      1560171      1460 rd       asymptotics    37     0.135       199984.  0.00321
# 58  1559684      1560171      1460 rd       percentile     37     0.135       202688.  0.00322
# 59  1559684      1560171      1460 rr       asymptotics    37     0.135           45.7 0.0866 
# 60  1559684      1560171      1460 rr       percentile     37     0.135           46.3 0.0860 

nonHrEstimatesWeighted <- nonHrEstimatesWeighted |>
  pivotNonHr()
results <- nonHrEstimatesWeighted |>
  group_by(targetId, comparatorId, timePoint, estimand, method) |>
  summarise(
    count = n(),
    type1Error = computeType1Error(lb, ub, h0 = if_else(estimand == "rr", 1, 0)),
    meanPrecision = computeMeanPrecision(se),
    ease = computeEase(if_else(estimand == "rr", log(estimate), estimate), se),
    .groups = "drop"
  )
print(results, n = 60)

# Forest plots -----------------------------------------------------------------
library(ggplot2)
cohortNames <- readr::read_csv(file.path(outputFolder, "cohortCounts.csv"), show_col_types = FALSE) |>
  filter(type == "exposure")

y <- hrEstimates |>
  inner_join(cohortNames |>
               select(targetId = cohortId, targetName = cohortName), by = join_by(targetId)) |>
  inner_join(cohortNames |>
               select(comparatorId = cohortId, comparatorName = cohortName), by = join_by(comparatorId)) |>
  mutate(label = paste(targetName, "vs", comparatorName)) |>
  group_by(label, targetId, comparatorId) |>
  arrange(outcomeId) |>
  mutate(y = row_number()) |>
  select(targetId, comparatorId, label, outcomeId, y)
vizData <- bind_rows(
  hrEstimates |>
    inner_join(y, by = join_by(targetId, comparatorId, outcomeId)) |>
    mutate(rr = exp(logRr)) |>
    select(label, y, rr, lb = ci95Lb, ub = ci95Ub, logRr, seLogRr) |>
    mutate(significant = if_else(lb > 1 | ub < 1, "Yes", "No"),
           estimand = "Hazard Ratio"),
  nonHrEstimates |>
    filter(timePoint == 1095, estimand == "rr", method == "asymptotics") |>
    inner_join(y, by = join_by(targetId, comparatorId, outcomeId)) |>
    mutate(logRr = log(estimate)) |>
    select(label, y, rr = estimate, lb, ub, logRr, seLogRr = se) |>
    mutate(significant = if_else(lb > 1 | ub < 1, "Yes", "No"),
           estimand = "Risk Ratio"),
  nonHrEstimates |>
    filter(timePoint == 1095, estimand == "rd", method == "asymptotics") |>
    inner_join(y, by = join_by(targetId, comparatorId, outcomeId)) |>
    select(label, y, rr = estimate, lb, ub, logRr = estimate, seLogRr = se) |>
    mutate(significant = if_else(lb > 0 | ub < 0, "Yes", "No"),
           estimand = "Risk Difference")
) 
null <- tibble(
  x = c(1, 1, 0),
  estimand = c("Hazard Ratio", "Risk Ratio", "Risk Difference")
)
ggplot(vizData, aes(x = rr, y = y, color = significant)) +
  geom_vline(aes(xintercept = x), data = null) +
  geom_point(alpha = 0.8) +
  geom_errorbarh(aes(xmin = lb, xmax = ub), alpha = 0.8) +
  scale_color_manual("Reject null hypothesis", values = c("#336B91", "#EB6622")) +
  scale_y_continuous("Negative control outcome") +
  facet_grid(label~estimand, scales = "free") +
  theme(axis.title.x = element_blank(),
        axis.text.y.left = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "bottom")
ggsave("RealWorldExample/EstimatesMatched.png", width = 8, height = 10, dpi = 300)

# Empirical calibration plots
groups <- vizData |>
  group_by(label, estimand) |>
  group_split()

for (group in groups) {
  fileName <- sprintf("RealWorldExample/NCs_%s_%s.png", gsub(" ", "_", group$label[1]), gsub(" ", "_", group$estimand[1]))
  EmpiricalCalibration::plotCalibrationEffect(
    logRrNegatives = group$logRr, 
    seLogRrNegatives = group$seLogRr,
    showExpectedAbsoluteSystematicError = TRUE,
    showCis = TRUE,
    title = group$label[1],
    xLabel = group$estimand[1],
    fileName = fileName
  )
}
