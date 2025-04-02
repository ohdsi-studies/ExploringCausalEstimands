# Some code to better understand causal estimands. Please ignore

library(survival)

# Simplest model: no confounding, constant baseline hazard, proportional hazards assumption holds ----
n <- 1000
baselineHazard <- 0.001
censorHazard <- 0.01
hr <- 2
iter <- 10000
pA <- 0.5

coverage <- 0
for (i in seq_len(iter)) {
  a <- rbinom(n, 1, pA)
  tY <- rexp(n, ifelse(a == 1, hr * baselineHazard, baselineHazard))  
  tC <- rexp(n, censorHazard)
  t <- pmin(tY, tC)
  y <- t == tY
  fit <- coxph(Surv(t, y) ~ a)
  ci <- exp(confint(fit))
  coverage <- coverage + (ci[1] <= hr & ci[2] >= hr)  
}
writeLines(sprintf("Coverage: %0.1f%%", 100 * coverage / iter))
# Coverage: 94.9%


# Non-constant baseline hazard (Weibull) -----------------------------------------------
n <- 1000
lambdaBaselineHazard <- 50 # Scale parameter for the baseline hazard
kBaselineHazard <- 1       # Shape parameter for the baseline hazard
plot(1:100, dweibull(1:100, kBaselineHazard, lambdaBaselineHazard))

censorHazard <- 0.01
hr <- 2
iter <- 10000
pA <- 0.5

coverage <- 0
for (i in seq_len(iter)) {
  a <- rbinom(n, 1, pA)
  delta <- 0.1
  # Simpler to use discrete-time:
  atRisk <- rep(TRUE, n)
  survivalTime <- rep(1001, n)
  for (t in seq_len(1000)) {
    nAtRisk <- sum(atRisk)
    if (nAtRisk == 0) {
      break
    }
    baselineHazard <- dweibull(t, kBaselineHazard, lambdaBaselineHazard)  
    hazard <- ifelse(a[atRisk] == 1, hr * baselineHazard, baselineHazard)
    outcome <- runif(nAtRisk) < hazard
    censored <- runif(nAtRisk) < censorHazard
    noLongerAtRisk <- outcome | censored
    survivalTime[atRisk][noLongerAtRisk] <- t
    y[atRisk] <- outcome
    atRisk[atRisk] <- !noLongerAtRisk
  }
  fit <- coxph(Surv(survivalTime, y) ~ a)
  ci <- exp(confint(fit))
  coverage <- coverage + (ci[1] <= hr & ci[2] >= hr)  
}
writeLines(sprintf("Coverage: %0.1f%%", 100 * coverage / iter))
# Coverage: 95.1%


# Non-proportionality (HR function over time, also Weibull) ----------------------------------------
n <- 1000
lambdaBaselineHazard <- 50 # Scale parameter for the baseline hazard
kBaselineHazard <- 1       # Shape parameter for the baseline hazard
censorHazard <- 0.01
lambdaLogHr <- 20          # Scale parameter for the hazard ratio
kLogHr <- 1.5              # Shape parameter for the hazard ratio
plot(1:100, dweibull(1:100, kLogHr, lambdaLogHr))

iter <- 10000
pA <- 0.5

coverage <- 0
for (i in seq_len(iter)) {
  a <- rbinom(n, 1, pA)

  atRisk <- rep(TRUE, n)
  survivalTime <- rep(1001, n)
  y <- rep(0, n)
  averageLogHr <- 0
  denominator <- 0
  for (t in seq_len(1000)) {
    nAtRisk <- sum(atRisk)
    if (nAtRisk == 0) {
      break
    }
    baselineHazard <- dweibull(t, kBaselineHazard, lambdaBaselineHazard)  
    logHr <- dweibull(t, kLogHr, lambdaLogHr)  
    averageLogHr <- averageLogHr + logHr * nAtRisk
    denominator <- denominator + nAtRisk
    hazard <- ifelse(a[atRisk] == 1, exp(logHr) * baselineHazard, baselineHazard)
    outcome <- runif(nAtRisk) < hazard
    censored <- runif(nAtRisk) < censorHazard
    noLongerAtRisk <- outcome | censored
    survivalTime[atRisk][noLongerAtRisk] <- t
    y[atRisk] <- outcome
    atRisk[atRisk] <- !noLongerAtRisk
  }
  averageLogHr <- averageLogHr / denominator
  fit <- coxph(Surv(survivalTime, y) ~ a)
  ci <- exp(confint(fit))
  coverage <- coverage + (ci[1] <= exp(averageLogHr) & ci[2] >= exp(averageLogHr))  
}
writeLines(sprintf("Coverage: %0.1f%%", 100 * coverage / iter))
# Coverage: 94.9%


# Depletion of susceptibles (susceptibe to treatment effecft) --------------------------------------
n <- 1000
lambdaBaselineHazard <- 50 # Scale parameter for the baseline hazard
kBaselineHazard <- 1       # Shape parameter for the baseline hazard
censorHazard <- 0.01
pSusceptible <- 0.05
lambdaLogHr <- 10          # Scale parameter for the hazard ratio
kLogHr <- 1.5              # Shape parameter for the hazard ratio
multiplierLogHr <- 100
plot(1:100, multiplierLogHr * dweibull(1:100, kLogHr, lambdaLogHr))

iter <- 1000
pA <- 0.5

coverage <- 0
for (i in seq_len(iter)) {
  a <- rbinom(n, 1, pA)
  susceptible <- rbinom(n, 1, pSusceptible)

  atRisk <- rep(TRUE, n)
  survivalTime <- rep(1001, n)
  y <- rep(0, n)
  averageLogHr <- 0
  denominator <- 0
  for (t in seq_len(1000)) {
    nAtRisk <- sum(atRisk)
    if (nAtRisk == 0) {
      break
    }
    baselineHazard <- dweibull(t, kBaselineHazard, lambdaBaselineHazard)  
    logHr <- multiplierLogHr * dweibull(t, kLogHr, lambdaLogHr)  
    averageLogHr <- averageLogHr + logHr * nAtRisk
    denominator <- denominator + nAtRisk
    hazard <- ifelse(a[atRisk] == 1 & susceptible[atRisk] == 1, exp(logHr) * baselineHazard, baselineHazard)
    outcome <- runif(nAtRisk) < hazard
    censored <- runif(nAtRisk) < censorHazard
    noLongerAtRisk <- outcome | censored
    survivalTime[atRisk][noLongerAtRisk] <- t
    y[atRisk] <- outcome
    atRisk[atRisk] <- !noLongerAtRisk
  }
  averageLogHr <- averageLogHr / denominator
  fit <- coxph(Surv(survivalTime, y) ~ a)
  ci <- exp(confint(fit))
  coverage <- coverage + (ci[1] <= exp(averageLogHr) & ci[2] >= exp(averageLogHr))  
}
writeLines(sprintf("Coverage: %0.1f%%", 100 * coverage / iter))
# Coverage: 94.8%

# Depletion of susceptibles (susceptibe to outcomwe) --------------------------------------
n <- 1000
lambdaBaselineHazard <- 50 # Scale parameter for the baseline hazard
kBaselineHazard <- 1       # Shape parameter for the baseline hazard
censorHazard <- 0.01
pSusceptible <- 0.05
lambdaLogHr <- 10          # Scale parameter for the hazard ratio
kLogHr <- 1.5              # Shape parameter for the hazard ratio
multiplierLogHr <- 100
plot(1:100, multiplierLogHr * dweibull(1:100, kLogHr, lambdaLogHr))

iter <- 10000
pA <- 0.5

coverage <- 0
for (i in seq_len(iter)) {
  a <- rbinom(n, 1, pA)
  susceptible <- rbinom(n, 1, pSusceptible)
  nSusceptible <- sum(susceptible)
  
  atRisk <- rep(TRUE, n)
  survivalTime <- rep(1001, n)
  y <- rep(0, n)
  averageLogHr <- 0
  denominator <- 0
  for (t in seq_len(1000)) {
    nAtRisk <- sum(atRisk)
    if (nAtRisk == nSusceptible) {
      break
    }
    baselineHazard <- dweibull(t, kBaselineHazard, lambdaBaselineHazard)  
    logHr <- multiplierLogHr * dweibull(t, kLogHr, lambdaLogHr)  
    averageLogHr <- averageLogHr + logHr * nAtRisk
    denominator <- denominator + nAtRisk
    hazard <- ifelse(a[atRisk] == 1, exp(logHr) * baselineHazard, baselineHazard)
    hazard <- ifelse(susceptible[atRisk], hazard, 0)
    outcome <- runif(nAtRisk) < hazard
    censored <- runif(nAtRisk) < censorHazard
    noLongerAtRisk <- outcome | censored
    survivalTime[atRisk][noLongerAtRisk] <- t
    y[atRisk] <- outcome
    atRisk[atRisk] <- !noLongerAtRisk
  }
  averageLogHr <- averageLogHr / denominator
  fit <- coxph(Surv(survivalTime, y) ~ a)
  ci <- exp(confint(fit))
  coverage <- coverage + (ci[1] <= exp(averageLogHr) & ci[2] >= exp(averageLogHr))  
}
writeLines(sprintf("Coverage: %0.1f%%", 100 * coverage / iter))
# Coverage: 92.3%

