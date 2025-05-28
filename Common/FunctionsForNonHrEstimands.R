library(CohortMethod)
library(survival)
library(dplyr)

computeWeights <- function(data) {
  weights <- ifelse(data$treatment == 1,
                    mean(data$treatment == 1),
                    mean(data$treatment == 0) * data$propensityScore / (1 - data$propensityScore)
  )
  return(weights)
}

.calculateEstimands <- function(dummy, data, timePoints, sample = FALSE) {
  if (sample) {
    indices <- sample.int(nrow(data), nrow(data), replace = TRUE)
    sampledData <- data[indices, ]
  } else {
    sampledData <- data
  }
  if (sum(sampledData$y[sampledData$treatment == 1]) == 0 |
      sum(sampledData$y[sampledData$treatment == 0]) == 0) {
    results <- list()
    for (i in seq_along(timePoints)) {
      timePoint <- timePoints[i]
      results[[i]] <- tibble(
        timePoint = timePoint,
        rr = 1,
        rd = 0
      )
    }
    results <- bind_rows(results)
    return(results)
  }
  sampledData$treatment <- as.factor(sampledData$treatment)  
  if ("weight" %in% colnames(sampledData)) {
    fit <- survfit(Surv(survivalTime, y) ~ treatment,
                   robust = FALSE,
                   data = sampledData,
                   weights = sampledData$weight)
  } else {
    fit <- survfit(Surv(survivalTime, y) ~ treatment,
                   robust = FALSE,
                   data = sampledData)
    
  }
  
  results <- list()
  for (i in seq_along(timePoints)) {
    timePoint <- timePoints[i]
    riskControl <- 1 - summary(fit,times = timePoint, extend=TRUE)$surv[1]
    riskTreated <- 1 - summary(fit,times = timePoint, extend=TRUE)$surv[2]
    results[[i]] <- tibble(
      timePoint = timePoint,
      rr = riskTreated / riskControl,
      rd = riskTreated - riskControl
    )
  }
  results <- bind_rows(results)
  return(results)
}

computeEstimands <- function(population, 
                             timePoints = c(180, 365, 730, 1095, 1460),
                             bootstrapSize = 1000, 
                             cluster = NULL) {
  if (is.null(cluster)) {
    cluster <- ParallelLogger::makeCluster(1)
    on.exit(ParallelLogger::stopCluster(cluster))
  }
  ParallelLogger::clusterRequire(cluster, "survival")
  ParallelLogger::clusterRequire(cluster, "dplyr")
  mainEstimates <- .calculateEstimands(NA, population, timePoints)
  bootStrap <- ParallelLogger::clusterApply(cluster, 
                                            seq_len(bootstrapSize), 
                                            .calculateEstimands, 
                                            data = population, 
                                            timePoints = timePoints,
                                            sample = TRUE)
  bootStrap <- bind_rows(bootStrap) 
  estimates <- bootStrap |>
    group_by(timePoint) |>
    summarise(
      seLogRrNormal = sqrt(var(log(rr))),
      lbRrPercentile = quantile(rr, 0.025),
      ubRrPercentile = quantile(rr, 0.975),
      seRdNormal = sqrt(var(rd)),
      lbRdPercentile = quantile(rd, 0.025),
      ubRdPercentile = quantile(rd, 0.975),
    ) |>
    inner_join(mainEstimates, by = join_by(timePoint)) |>
    mutate(
      lbRrAsymptotics = exp(log(rr) + qnorm(0.025) * seLogRrNormal),
      ubRrAsymptotics = exp(log(rr) + qnorm(0.975) * seLogRrNormal),
      lbRdAsymptotics = exp(rd + qnorm(0.025) * seRdNormal),
      ubRdAsymptotics = exp(rd + qnorm(0.975) * seRdNormal),
      seLogRrPercentile = (log(ubRrPercentile) - log(lbRrPercentile)) / (2 * qnorm(0.975)),
      seRdPercentile = (ubRdPercentile - lbRdPercentile) / (2 * qnorm(0.975))
    )
  return(estimates)
}
