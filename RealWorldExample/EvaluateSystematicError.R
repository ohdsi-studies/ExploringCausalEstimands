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
  summarise(
    count = n(),
    type1Error = computeType1Error(ci95Lb, ci95Ub, h0 = 1),
    meanPrecision = computeMeanPrecision(seLogRr),
    ease = computeEase(logRr, seLogRr)
  )
# # A tibble: 1 × 4
# count type1Error meanPrecision   ease
# <int>      <dbl>         <dbl>  <dbl>
#   1    37      0.189          137. 0.0705

nonHrEstimates |>
  pivotNonHr()
group_by(timePoint, estimand, method) |>
  summarise(
    count = n(),
    type1Error = computeType1Error(lb, ub, h0 = if_else(estimand == "rr", 1, 0)),
    meanPrecision = computeMeanPrecision(se),
    ease = computeEase(if_else(estimand == "rr", log(estimate), estimate), se),
    .groups = "drop"
  )
# # A tibble: 20 × 7
# timePoint estimand method      count type1Error meanPrecision    ease
# <dbl> <chr>    <chr>       <int>      <dbl>         <dbl>   <dbl>
# 1       180 rd       asymptotics    37     0.108      3090103.  0.00212
# 2       180 rd       percentile     37     0.108      3095688.  0.00212
# 3       180 rr       asymptotics    37     0.0811          51.2 0.0726 
# 4       180 rr       percentile     37     0.108           50.9 0.0721 
# 5       365 rd       asymptotics    37     0.189      1164009.  0.00244
# 6       365 rd       percentile     37     0.189      1162677.  0.00244
# 7       365 rr       asymptotics    37     0.189           68.8 0.0847 
# 8       365 rr       percentile     37     0.189           68.0 0.0839 
# 9       730 rd       asymptotics    37     0.216       361188.  0.00277
# 10       730 rd       percentile     37     0.189       368155.  0.00277
# 11       730 rr       asymptotics    37     0.216           78.6 0.0579 
# 12       730 rr       percentile     37     0.189           80.1 0.0589 
# 13      1095 rd       asymptotics    37     0.189       153486.  0.00351
# 14      1095 rd       percentile     37     0.189       154635.  0.00351
# 15      1095 rr       asymptotics    37     0.189           71.8 0.0802 
# 16      1095 rr       percentile     37     0.189           71.9 0.0801 
# 17      1460 rd       asymptotics    37     0.189        76434.  0.00461
# 18      1460 rd       percentile     37     0.162        76982.  0.00462
# 19      1460 rr       asymptotics    37     0.189           62.2 0.0996 
# 20      1460 rr       percentile     37     0.162           62.7 0.100  


nonHrEstimatesWeighted |>
  pivotNonHr() |>
  group_by(timePoint, estimand, method) |>
  summarise(
    count = n(),
    type1Error = computeType1Error(lb, ub, h0 = if_else(estimand == "rr", 1, 0)),
    meanPrecision = computeMeanPrecision(se),
    ease = computeEase(if_else(estimand == "rr", log(estimate), estimate), se),
    .groups = "drop"
  )
# # A tibble: 20 × 7
# timePoint estimand method      count type1Error meanPrecision    ease
# <dbl> <chr>    <chr>       <int>      <dbl>         <dbl>   <dbl>
# 1       180 rd       asymptotics    37      0.108     3730245.  0.00209
# 2       180 rd       percentile     37      0.135     3817377.  0.00209
# 3       180 rr       asymptotics    37      0.108          64.9 0.0674 
# 4       180 rr       percentile     37      0.135          65.9 0.0673 
# 5       365 rd       asymptotics    37      0.162     1478556.  0.00237
# 6       365 rd       percentile     37      0.162     1493983.  0.00237
# 7       365 rr       asymptotics    37      0.162          91.1 0.0802 
# 8       365 rr       percentile     37      0.162          91.9 0.0799 
# 9       730 rd       asymptotics    37      0.108      440635.  0.00264
# 10       730 rd       percentile     37      0.108      442782.  0.00264
# 11       730 rr       asymptotics    37      0.108         101.  0.0512 
# 12       730 rr       percentile     37      0.108         101.  0.0506 
# 13      1095 rd       asymptotics    37      0.162      201158.  0.00310
# 14      1095 rd       percentile     37      0.189      202959.  0.00310
# 15      1095 rr       asymptotics    37      0.162          98.0 0.0558 
# 16      1095 rr       percentile     37      0.189          98.7 0.0550 
# 17      1460 rd       asymptotics    37      0.162       97913.  0.00378
# 18      1460 rd       percentile     37      0.162       98917.  0.00377
# 19      1460 rr       asymptotics    37      0.162          82.2 0.0673 
# 20      1460 rr       percentile     37      0.162          83.1 0.0671 

# Forest plots -----------------------------------------------------------------

library(ggplot2)
y <- hrEstimates |>
  arrange(targetId, comparatorId, outcomeId) |>
  mutate(y = row_number()) |>
  select(targetId, comparatorId, outcomeId, y)
vizData <- bind_rows(
  hrEstimates |>
    inner_join(y, by = join_by(targetId, comparatorId, outcomeId)) |>
    mutate(rr = exp(logRr)) |>
    select(y, rr, lb = ci95Lb, ub = ci95Ub) |>
    mutate(significant = if_else(lb > 1 | ub < 1, "Yes", "No"),
           estimand = "Hazard Ratio"),
  nonHrEstimates |>
    filter(timePoint == 1095, estimand == "rr", method == "asymptotics") |>
    inner_join(y, by = join_by(targetId, comparatorId, outcomeId)) |>
    select(y, rr = estimate, lb, ub) |>
    mutate(significant = if_else(lb > 1 | ub < 1, "Yes", "No"),
           estimand = "Risk Ratio"),
  nonHrEstimates |>
    filter(timePoint == 1095, estimand == "rd", method == "asymptotics") |>
    inner_join(y, by = join_by(targetId, comparatorId, outcomeId)) |>
    select(y, rr = estimate, lb, ub) |>
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
  facet_grid(~estimand, scales = "free_x") +
  theme(axis.title.x = element_blank(),
        axis.text.y.left = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "bottom")
ggsave("RealWorldExample/EstimatesMatched.png", width = 8, height = 6, dpi = 300)
