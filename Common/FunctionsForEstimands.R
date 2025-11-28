library(CohortMethod)
library(survival)
library(dplyr)
library(adjustedCurves)

computeWeights <- function(data) {
  weights <- ifelse(data$treatment == 1,
                    mean(data$treatment == 1),
                    mean(data$treatment == 0) * data$propensityScore / (1 - data$propensityScore)
  )
  return(weights)
}

.calculateKmEstimands <- function(dummy, data, timePoints, sample = FALSE) {
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
        rd = 0,
        chr = 1
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
    cumHazardControl <- summary(fit,times = timePoint, extend=TRUE)$cumhaz[1]
    cumHazardTreated <- summary(fit,times = timePoint, extend=TRUE)$cumhaz[2]
    results[[i]] <- tibble(
      timePoint = timePoint,
      rr = riskTreated / riskControl,
      rd = riskTreated - riskControl,
      chr = cumHazardTreated / cumHazardControl
    )
  }
  results <- bind_rows(results)
  return(results)
}

.truncateData <- function(data, timePoint) {
  data <- data |>
    mutate(y = if_else(survivalTime > timePoint, 0, y),
           survivalTime = if_else(survivalTime > timePoint, timePoint, survivalTime))
  return(data)
}

.calculateAcceleratedFailureTime <- function(dummy, data, sample = FALSE) {
  if (sample) {
    indices <- sample.int(nrow(data), nrow(data), replace = TRUE)
    sampledData <- data[indices, ]
  } else {
    sampledData <- data
  }
  if ("weight" %in% colnames(sampledData)) {
    weights <- sampledData$weight
  } else {
    weights <- NULL
  }
  fit <- tryCatch(survreg(Surv(survivalTime, y) ~ treatment,
                          data = sampledData,
                          control = survreg.control(maxiter = 300),
                          dist = "weibull",
                          weights = weights),
                  error = function(e) {NA}
  )
  if (isTRUE(is.na(fit))) {
    fit <- survreg(Surv(survivalTime, y) ~ treatment,
                   data = sampledData,
                   control = survreg.control(maxiter = 300),
                   dist = "weibull",
                   scale = 1,
                   weights = weights)
  }
  logAf <- coef(fit)["treatment1"]
  af <- exp(logAf)
  logCi <- confint(fit, parm = "treatment1")
  ci <- exp(logCi)
  se <- (logCi[2] - logCi[1]) / (2 * qnorm(0.975))
  estimate <- tibble(
    estimate = af,
    lb = ci[1],
    ub = ci[2],
    se = se
  )
  return(estimate)
}

.calculateCoxEstimate <- function(dummy, data, sample = FALSE) {
  if (sample) {
    indices <- sample.int(nrow(data), nrow(data), replace = TRUE)
    sampledData <- data[indices, ]
  } else {
    sampledData <- data
  }
  if (sum(sampledData$y[sampledData$treatment == 1]) == 0 |
      sum(sampledData$y[sampledData$treatment == 0]) == 0) {
    if (sample) {
      estimate <- tibble(
        logHr = NA
      )
    } else {
      estimate <- tibble(
        estimate = NA,
        lb = NA,
        ub = NA,
        se = NA
      )
    }
  } else {
    cyclopsData <- Cyclops::createCyclopsData(Surv(survivalTime, y) ~ treatment, modelType = "cox", data = sampledData)
    fit <- Cyclops::fitCyclopsModel(cyclopsData)
    if (fit$return_flag != "SUCCESS") {
      estimate <- tibble(
        estimate = NA,
        lb = NA,
        ub = NA,
        se = NA
      )
    } else {
      logHr <- coef(fit)
      if (sample) {
        estimate <- tibble(
          logHr = logHr
        )
      } else {
        hr <- exp(logHr)
        ci <- exp(confint(fit, parm = "treatment1")[c(2, 3)])
        se <- (log(ci[2]) - log(ci[1])) / (qnorm(0.975) - qnorm(0.025))
        estimate <- tibble(
          estimate = hr,
          lb = ci[1],
          ub = ci[2],
          se = se
        )
      }
    }
  }
  return(estimate)
}

.calculateRmst <- function(dummy, data, timePoints, sample = FALSE) {
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
        rmstDiff = NA,
        rmstRatio = NA
      )
    }
    results <- bind_rows(results)
    return(results)
  }
  lastDay <- sampledData |>
    filter(y == 1) |>
    group_by(treatment) |>
    summarise(time = max(survivalTime)) |>
    summarise(min(time)) |>
    pull()
  timePointTruncation <- tibble(timePoint = timePoints) |>
    mutate(truncatedTimePoint = pmin(timePoint, lastDay))
  
  if ("weight" %in% colnames(sampledData)) {
    survObject <- adjustedCurves::adjustedsurv(
      data = sampledData,
      variable = "treatment",
      ev_time = "survivalTime",
      event = "y",
      method = "iptw_km",
      treatment_model = sampledData$weight
    )
  } else {
    survObject <- adjustedCurves::adjustedsurv(
      data = sampledData,
      variable = "treatment",
      ev_time = "survivalTime",
      event = "y",
      method = "km"
    )
  }
  rmstDiff <- adjustedCurves::adjusted_rmst(
    adjsurv = survObject,
    to = unique(timePointTruncation$truncatedTimePoint),
    contrast = "diff"
  )
  rmstRatio <- adjustedCurves::adjusted_rmst(
    adjsurv = survObject,
    to = unique(timePointTruncation$truncatedTimePoint),
    contrast = "ratio"
  )
  results <- timePointTruncation |>
    inner_join(rmstDiff |>
                 rename(truncatedTimePoint = "to",
                        rmstDiff = "diff"),
               by = join_by(truncatedTimePoint)
    ) |>
    inner_join(rmstRatio |>
                 rename(truncatedTimePoint = "to",
                        rmstRatio = "ratio"),
               by = join_by(truncatedTimePoint)
    ) |>
    select(-truncatedTimePoint)
  return(results)
}

computeEstimands <- function(population, 
                             timePoints = c(180, 365, 730, 1095, 1460),
                             bootstrapSize = 200, 
                             cluster = NULL) {
  if (is.null(cluster)) {
    cluster <- ParallelLogger::makeCluster(1)
    on.exit(ParallelLogger::stopCluster(cluster))
  }
  ParallelLogger::clusterRequire(cluster, "survival")
  ParallelLogger::clusterRequire(cluster, "dplyr")
  if ("a" %in% colnames(population)) {
    population <- population |>
      rename(treatment = a)
  }
  population$treatment <- as.factor(population$treatment)  
  
  # Kaplan-Meier-based estimates:
  mainKmEstimates <- .calculateKmEstimands(NA, population, timePoints)
  kmBootStrap <- ParallelLogger::clusterApply(cluster, 
                                              seq_len(bootstrapSize), 
                                              .calculateKmEstimands, 
                                              data = population, 
                                              timePoints = timePoints,
                                              sample = TRUE)
  kmBootStrap <- bind_rows(kmBootStrap) 
  kmRrAsymptotic <- kmBootStrap |>
    group_by(timePoint) |>
    summarise(
      se = sqrt(var(log(rr), na.rm = TRUE))
    ) |>
    inner_join(mainKmEstimates, by = join_by(timePoint)) |>
    transmute(
      timePoint = timePoint,
      estimate = rr,
      lb = exp(log(rr) + qnorm(0.025) * se),
      ub = exp(log(rr) + qnorm(0.975) * se),
      se = se
    ) |>
    mutate(estimand = "Risk Ratio", 
           model = "Kaplan Meier",
           contrast = "ratio",
           ciMethod = "asymptotic")
  kmRdAsymptotic <- kmBootStrap |>
    group_by(timePoint) |>
    summarise(
      se = sqrt(var(rd, na.rm = TRUE))
    ) |>
    inner_join(mainKmEstimates, by = join_by(timePoint)) |>
    transmute(
      timePoint = timePoint,
      estimate = rd,
      lb = rd + qnorm(0.025) * se,
      ub = rd + qnorm(0.975) * se,
      se = se
    ) |>
    mutate(estimand = "Risk Difference", 
           model = "Kaplan Meier",
           contrast = "difference",
           ciMethod = "asymptotic")
  kmChrAsymptotic <- kmBootStrap |>
    group_by(timePoint) |>
    summarise(
      se = sqrt(var(log(chr), na.rm = TRUE))
    ) |>
    inner_join(mainKmEstimates, by = join_by(timePoint)) |>
    transmute(
      timePoint = timePoint,
      estimate = chr,
      lb = exp(log(chr) + qnorm(0.025) * se),
      ub = exp(log(chr) + qnorm(0.975) * se),
      se = se
    ) |>
    mutate(estimand = "Cumulative Hazard Ratio", 
           model = "Kaplan Meier",
           contrast = "ratio",
           ciMethod = "asymptotic")
  kmRrPercentile <- kmBootStrap |>
    group_by(timePoint) |>
    summarise(
      lb = quantile(rr, 0.025, na.rm = TRUE),
      ub = quantile(rr, 0.975, na.rm = TRUE),
    ) |>
    inner_join(mainKmEstimates, by = join_by(timePoint)) |>
    transmute(
      timePoint = timePoint,
      estimate = rr,
      lb = lb,
      ub = ub,
      se = (log(ub) - log(lb)) / (2 * qnorm(0.975)),
    ) |>
    mutate(estimand = "Risk Ratio", 
           model = "Kaplan Meier",
           contrast = "ratio",
           ciMethod = "percentile")
  kmRdPercentile <- kmBootStrap |>
    group_by(timePoint) |>
    summarise(
      lb = quantile(rd, 0.025, na.rm = TRUE),
      ub = quantile(rd, 0.975, na.rm = TRUE),
    ) |>
    inner_join(mainKmEstimates, by = join_by(timePoint)) |>
    transmute(
      timePoint = timePoint,
      estimate = rd,
      lb = lb,
      ub = ub,
      se = (ub - lb) / (2 * qnorm(0.975)),
    ) |>
    mutate(estimand = "Risk Difference", 
           model = "Kaplan Meier",
           contrast = "difference",
           ciMethod = "percentile")
  kmChrPercentile <- kmBootStrap |>
    group_by(timePoint) |>
    summarise(
      lb = quantile(chr, 0.025, na.rm = TRUE),
      ub = quantile(chr, 0.975, na.rm = TRUE),
    ) |>
    inner_join(mainKmEstimates, by = join_by(timePoint)) |>
    transmute(
      timePoint = timePoint,
      estimate = chr,
      lb = lb,
      ub = ub,
      se = (log(ub) - log(lb)) / (2 * qnorm(0.975)),
    ) |>
    mutate(estimand = "Cumulative Hazard Ratio", 
           model = "Kaplan Meier",
           contrast = "ratio",
           ciMethod = "percentile")
  
  kmEstimates <- bind_rows(
    kmRrAsymptotic,
    kmRrPercentile,
    kmRdAsymptotic,
    kmRdPercentile,
    kmChrAsymptotic,
    kmChrPercentile
  )
  
  # Cox and accelerated failure time estimates:
  afHrEstimates <- list()
  for (i in seq_along(timePoints)) {
    timePoint <- timePoints[i]
    truncatedData <- .truncateData(population, timePoint)
    aftMainEstimate <- .calculateAcceleratedFailureTime(NA, truncatedData)
    aftBootStrap <- ParallelLogger::clusterApply(cluster, 
                                                seq_len(bootstrapSize), 
                                                .calculateAcceleratedFailureTime, 
                                                data = truncatedData, 
                                                sample = TRUE)
    aftBootStrap <- bind_rows(aftBootStrap) 
    aftAsymptotic <- aftMainEstimate |>
      mutate(estimand = "Acceleration Coefficient", 
             model = "AFT",
             contrast = "ratio",
             ciMethod = "asymptotic")
    aftPercentile <- tibble(
      estimate = aftMainEstimate$estimate,
      lb = quantile(aftBootStrap$estimate, 0.025, na.rm = TRUE),
      ub = quantile(aftBootStrap$estimate, 0.975, na.rm = TRUE),
      se = sd(aftBootStrap$estimate, na.rm = TRUE),
      estimand = "Acceleration Coefficient", 
      model = "AFT",
      contrast = "ratio",
      ciMethod = "percentile"
    )
    if ("weight" %in% colnames(truncatedData)) {
      hrMainEstimate <- .calculateCoxEstimate(NA, truncatedData)
      hrBootStrap <- ParallelLogger::clusterApply(cluster, 
                                                  seq_len(bootstrapSize), 
                                                  .calculateCoxEstimate, 
                                                  data = population, 
                                                  sample = TRUE)
      hrBootStrap <- bind_rows(hrBootStrap) 
      seLogHr <- sqrt(var(hrBootStrap$logHr))
      ci <- exp(log(hrMainEstimate$estimate) + qnorm(c(0.025, 0.975)) * seLogHr)
      hrEstimate <- tibble(
        estimate = hrMainEstimate$estimate,
        lb = ci[1],
        ub = ci[2],
        se = seLogHr,
        estimand = "Hazard Ratio", 
        model = "Cox",
        contrast = "ratio",
        ciMethod = "asymptotic"
      )
    } else {
      hrEstimate <- .calculateCoxEstimate(NA, truncatedData) |>
        mutate(estimand = "Hazard Ratio", 
               model = "Cox",
               contrast = "ratio",
               ciMethod = "asymptotic")
    }
    afHrEstimates[[i]] <- bind_rows(
      aftAsymptotic,
      aftPercentile,
      hrEstimate
    ) |>
      mutate(timePoint = timePoint)
  }
  afHrEstimates <- bind_rows(afHrEstimates)
  
  # RMST
  mainMrstEstimates <- .calculateRmst(NA, population, timePoints)
  mrstBootStrap <- ParallelLogger::clusterApply(cluster, 
                                                seq_len(bootstrapSize), 
                                                .calculateRmst, 
                                                data = population, 
                                                timePoints = timePoints,
                                                sample = TRUE)
  mrstBootStrap <- bind_rows(mrstBootStrap) 
  mrstDiffAsymptotic <- mrstBootStrap |>
    group_by(timePoint) |>
    summarise(
      se = sqrt(var(rmstDiff, na.rm = TRUE))
    ) |>
    inner_join(mainMrstEstimates, by = join_by(timePoint)) |>
    transmute(
      timePoint = timePoint,
      estimate = rmstDiff,
      lb = rmstDiff + qnorm(0.025) * se,
      ub = rmstDiff + qnorm(0.975) * se,
      se = se
    ) |>
    mutate(estimand = "RMST Difference", 
           model = "RMST",
           contrast = "difference",
           ciMethod = "asymptotic")
  mrstRatioAsymptotic <- mrstBootStrap |>
    group_by(timePoint) |>
    summarise(
      se = sqrt(var(log(rmstRatio), na.rm = TRUE))
    ) |>
    inner_join(mainMrstEstimates, by = join_by(timePoint)) |>
    transmute(
      timePoint = timePoint,
      estimate = rmstRatio,
      lb = exp(log(rmstRatio) + qnorm(0.025) * se),
      ub = exp(log(rmstRatio) + qnorm(0.975) * se),
      se = se
    ) |>
    mutate(estimand = "RMST Ratio", 
           model = "RMST",
           contrast = "ratio",
           ciMethod = "asymptotic")
  mrstDiffPercentile <- mrstBootStrap |>
    group_by(timePoint) |>
    summarise(
      lb = quantile(rmstDiff, 0.025, na.rm = TRUE),
      ub = quantile(rmstDiff, 0.975, na.rm = TRUE)
    ) |>
    inner_join(mainMrstEstimates, by = join_by(timePoint)) |>
    transmute(
      timePoint = timePoint,
      estimate = rmstDiff,
      lb = lb,
      ub = ub,
      se = (ub - lb) / (2 * qnorm(0.975))
    ) |>
    mutate(estimand = "RMST Difference", 
           model = "RMST",
           contrast = "difference",
           ciMethod = "percentile")
  mrstRatioPercentile <- mrstBootStrap |>
    group_by(timePoint) |>
    summarise(
      lb = quantile(rmstRatio, 0.025, na.rm = TRUE),
      ub = quantile(rmstRatio, 0.975, na.rm = TRUE)
    ) |>
    inner_join(mainMrstEstimates, by = join_by(timePoint)) |>
    transmute(
      timePoint = timePoint,
      estimate = rmstRatio,
      lb = lb,
      ub = ub,
      se = (log(ub) - log(lb)) / (2 * qnorm(0.975))
    ) |>
    mutate(estimand = "RMST Ratio", 
           model = "RMST",
           contrast = "ratio",
           ciMethod = "percentile")
  mrstEstimates <- bind_rows(
    mrstDiffAsymptotic,
    mrstDiffPercentile,
    mrstRatioAsymptotic,
    mrstRatioPercentile
  )
  estimates <- bind_rows(
    kmEstimates,
    afHrEstimates,
    mrstEstimates
  )
  return(estimates)
}
