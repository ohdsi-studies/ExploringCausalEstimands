library(timeEL)
library(dplyr)

.estimateRd <- function(t, population) {
  result <- TwoSampleKaplanMeier(
    time = population$survivalTime,
    status = population$y,
    group = population$a,
    t = t,
    contr = list(tol = 1e-05, algo = 2, k = 3, Trace = TRUE, method = "EL")
  )
  estimate <- tibble(
    timePoint = t,
    estimate = result$table.Diff["EL", "est."],
    lb = result$table.Diff["EL", "lower"],
    ub = result$table.Diff["EL", "upper"]
  )
  return(estimate)
}

estimateRiskDifference <- function(population, 
                                   timePoints = c(180, 365, 730, 1095, 1460)) {
  estimates <- lapply(timePoints, .estimateRd, population = population)
  estimates <- bind_rows(estimates) 
  return(estimates)
}