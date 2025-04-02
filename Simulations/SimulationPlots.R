# Code to generate plots for comparing different estimands

library(survival)
library(ggplot2)
library(dplyr)
library(Cyclops)

n <- 2500                  # Number of persons to simulate
pA <- 0.5                  # Probability of being exposed
lambdaBaselineHazard <- 50 # Scale parameter for the baseline hazard
kBaselineHazard <- 1       # Shape parameter for the baseline hazard
censorHazard <- 0.01       # Daily probability of being censored
pSusceptible <- 0.20       # Probability of being susceptible to exposure effect     
lambdaLogHr <- 10          # Scale parameter for the hazard ratio
kLogHr <- 1.5              # Shape parameter for the hazard ratio
multiplierLogHr <- 25      # Multiplier for the hazard ratio

baselineHazardFunction <- function(t) {
  dweibull(t, kBaselineHazard, lambdaBaselineHazard)
}

logHrFunction <- function(t) {
  multiplierLogHr * dweibull(t, kLogHr, lambdaLogHr)
}

# Plots of the generative process -------------------------------
x <- 0:100
vizData <- tibble(
  x = x,
  y = baselineHazardFunction(x)
)
ggplot(vizData, aes(x = x, y = y)) +
  geom_line(color = "#336B91", linewidth = 1) +
  geom_hline(yintercept = 0) +
  scale_x_continuous("Time (days)") +
  scale_y_continuous("Hazard (daily risk of outcome)") +
  theme(
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = "bottom"
  )
ggsave(filename = "Simulations/BaselineHazard.png", width = 5, height = 3.5, dpi = 300)

vizData <- bind_rows(
  tibble(
    x = x,
    y = logHrFunction(x),
    label = "Within susceptibles"
  )
)
ggplot(vizData, aes(x = x, y = y, color = label, group = label)) +
  geom_hline(yintercept = 0) +
  geom_line(linewidth = 1, alpha = 0.7) +
  scale_x_continuous("Time (days)") +
  scale_y_continuous("Hazard Ratio", breaks = log(c(1, 2, 3, 4, 5, 6, 7)), labels = c(1, 2, 3, 4, 5, 6,7 )) +
  scale_color_manual(values = c("#336B91")) +
  theme(
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = "right"
  )
ggsave(filename = "Simulations/HazardRatio.png", width = 7, height = 3.5, dpi = 300)

# Simulation ---------------------------------------------------------
set.seed(123)
a <- rbinom(n, 1, pA)
susceptible <- rbinom(n, 1, pSusceptible)

atRisk <- rep(TRUE, n)
survivalTime <- rep(1001, n)
y <- rep(0, n)
averageLogHr <- 0
denominator <- 0
targetOverTimme <- rep(NA, 1000)
targetSusceptiblesOverTime <- rep(NA, 1000)
for (t in seq_len(1000)) {
  nAtRisk <- sum(atRisk)
  if (nAtRisk == 0) {
    break
  }
  targetOverTime[t] <- sum(atRisk & a)
  targetSusceptiblesOverTime[t] <- sum(atRisk & a & susceptible)
  
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
cyclopsData <- createCyclopsData(Surv(survivalTime, y) ~ a, modelType = "cox")
fit <- fitCyclopsModel(cyclopsData)
hr <- exp(coef(fit))
ci <- exp(confint(fit, parm = "a")[c(2, 3)])
fit <- coxph(Surv(survivalTime, y) ~ a)
ci <- exp(confint(fit))
print(sprintf("Hazard ratio = %0.2f (95%% CI: %0.2f - %0.2f)",
              hr,
              ci[1],
              ci[2]))
# [1] "Hazard ratio = 1.22 (95% CI: 1.15 - 1.29)"

# Plots of estimates -------------------------------------------------------------------------------

# Depletion of susceptibles
vizData <- tibble(
  x = x,
  y = c(pSusceptible, targetSusceptiblesOverTime[1:100] / targetOverTime[1:100])
)
ggplot(vizData, aes(x = x, y = y)) +
  geom_hline(yintercept = 0) +
  geom_line(linewidth = 1, color = "#336B91", alpha = 0.7) +
  scale_x_continuous("Time (days)") +
  scale_y_continuous("Fraction of target susceptible") +
  theme(
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = "bottom"
  )
ggsave(filename = "Simulations/FractionSusceptible.png", width = 5, height = 3.5, dpi = 300)

# True effect given depletion of susceptibles
vizData <- bind_rows(
  tibble(
    x = x,
    y = logHrFunction(x),
    label = "Within susceptibles"
  ),
  tibble(
    x = x,
    y = logHrFunction(x) * c(pSusceptible, targetSusceptiblesOverTime[1:100] / targetOverTime[1:100]),
    label = "Average over target"
  ),
  tibble(
    x = x,
    y = logHrFunction(x) * pSusceptible,
    label = "Average without depletion"
  )
)
ggplot(vizData, aes(x = x, y = y, color = label, group = label)) +
  geom_hline(yintercept = 0) +
  geom_line(linewidth = 1, alpha = 0.7) +
  scale_x_continuous("Time (days)") +
  scale_y_continuous("Hazard Ratio", breaks = log(c(1, 2, 3, 4, 5, 6, 7)), labels = c(1, 2, 3, 4, 5, 6,7 )) +
  scale_color_manual(values = c("#EB6622", "#11A08A", "#336B91")) +
  theme(
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = "right"
  )
ggsave(filename = "Simulations/HazardRatioWithSusceptibles.png", width = 7, height = 3.5, dpi = 300)

# Show HR estimate from Cox model
estimate <- tibble(
  x = c(0, 100),
  y = rep(log(hr), 2),
  ymin = rep(log(ci[1]), 2),
  ymax = rep(log(ci[2]), 2),
  label = "Estimate from Cox model"
)
ggplot(vizData, aes(x = x, y = y, color = label, group = label)) +
  geom_hline(yintercept = 0) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.35, fill = "#FBC511", size = 0, data = estimate) +
  geom_line(data = estimate) +
  geom_line(linewidth = 1, alpha = 0.7) +
  scale_x_continuous("Time (days)") +
  scale_y_continuous("Hazard Ratio", breaks = log(c(1, 2, 3, 4, 5, 6, 7)), labels = c(1, 2, 3, 4, 5, 6,7 )) +
  scale_color_manual(values = c("#EB6622", "#11A08A", "#FBC511", "#336B91")) +
  theme(
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = "right"
  )
ggsave(filename = "Simulations/HazardRatioWithEstimate.png", width = 7, height = 3.5, dpi = 300)

# KM curves
km <- survfit(Surv(survivalTime, y) ~ a)
kms <- summary(km)
targetIdx <- kms$strata == "a=1"
timeIdx <- kms$time <= 100
vizData <- bind_rows(
  tibble(
    x = kms$time[targetIdx & timeIdx],
    y = kms$surv[targetIdx & timeIdx],
    ymin = kms$lower[targetIdx & timeIdx],
    ymax = kms$upper[targetIdx & timeIdx],
    label = "Target"
  ),
  tibble(
    x = kms$time[!targetIdx & timeIdx],
    y = kms$surv[!targetIdx & timeIdx],
    ymin = kms$lower[!targetIdx & timeIdx],
    ymax = kms$upper[!targetIdx & timeIdx],
    label = "Comparator"
  )
)

ggplot(vizData, aes(x = x, y = y, color = label, fill = label, group = label)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.35,  size = 0) +
  geom_line(linewidth = 1, alpha = 0.7) +
  scale_x_continuous("Time (days)") +
  scale_y_continuous("Survival probability") +
  scale_color_manual(values = c("#336B91", "#EB6622")) +
  scale_fill_manual(values = c("#336B91", "#EB6622")) +
  theme(
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = "bottom"
  )
ggsave(filename = "Simulations/KaplanMeier.png", width = 5, height = 3.5, dpi = 300)

# Risk ratio over time
vizData <- bind_rows(
  tibble(
    x = x,
    y = logHrFunction(x),
    label = "Within susceptibles"
  ),
  tibble(
    x = x,
    y = logHrFunction(x) * c(pSusceptible, targetSusceptiblesOverTime[1:100] / targetOverTime[1:100]),
    label = "Average over target"
  ),
  tibble(
    x = x,
    y = logHrFunction(x) * pSusceptible,
    label = "Average without depletion"
  )
)
calculate_rr_bootstrap <- function(data, sample = FALSE) {
  if (sample) {
    boot_data <- data[sample.int(nrow(data), nrow(data), replace = TRUE), ]
  } else {
    boot_data <- data
  }
  
  boot_km <- survfit(Surv(survivalTime, y) ~ a, data = boot_data)
  boot_kms <- summary(boot_km)
  boot_kmsDf <- as_tibble(boot_kms[c("time", "surv", "strata")])
  
  boot_kmEstimate <- inner_join(
    boot_kmsDf |>
      filter(strata == "a=1", time <= 100) |>
      transmute(pTarget = 1 - surv, x = time),
    boot_kmsDf |>
      filter(strata == "a=0", time <= 100) |>
      transmute(pComparator = 1 - surv, x = time),
    by = join_by(x)
  ) |>
    mutate(riskRatio = pTarget / pComparator) |>
    right_join(tibble(x = 1:100), by = join_by(x)) |>
    arrange(x) |>
    select(x, riskRatio) |>
    mutate(riskRatio = if_else(is.na(riskRatio) & x == 1, 1, riskRatio)) |> 
    tidyr::fill(riskRatio, .direction = "down") 
  
  return(boot_kmEstimate)
}
data = data.frame(survivalTime, y, a)
kmEstimate <- calculate_rr_bootstrap(data)
bootstrap <- lapply(1:1000, function(x) calculate_rr_bootstrap(data, sample = TRUE))
bootstrap <- bind_rows(bootstrap)
ci <- bootstrap |>
  group_by(x) |>
  summarise(lower = quantile(riskRatio, 0.025),
            upper = quantile(riskRatio, 0.975))
kmEstimate <- kmEstimate |>
  inner_join(ci, by = join_by(x)) |>
  mutate(label = "KM risk ratio")

ggplot(vizData, aes(x = x, y = y, color = label, group = label)) +
  geom_hline(yintercept = 0) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.35, fill = "#FBC511", size = 0, data = estimate) +
  geom_line(linewidth = 1, alpha = 0.7, data = estimate) +
  geom_ribbon(aes(ymin = log(lower), ymax = log(upper), y = log(riskRatio)), alpha = 0.35, fill = "#69AED5", size = 0, data = kmEstimate) +
  geom_line(aes(y = log(riskRatio)), linewidth = 1, alpha = 0.7, data = kmEstimate) +
  geom_line(linewidth = 1, alpha = 0.7) +
  scale_x_continuous("Time (days)") +
  scale_y_continuous("Hazard Ratio", 
                     breaks = log(c(1, 2, 3, 4, 5, 6, 7)), 
                     labels = c(1, 2, 3, 4, 5, 6, 7),
                     sec.axis = sec_axis(transform = ~., 
                                         breaks = log(c(1, 2, 3, 4, 5, 6, 7)), 
                                         labels = c(1, 2, 3, 4, 5, 6, 7),
                                         name = "Risk ratio")) +
  scale_color_manual(values = c("#EB6622", "#11A08A", "#FBC511", "#69AED5", "#336B91")) +
  theme(
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = "right"
  )

ggsave(filename = "Simulations/HrsAndRrs.png", width = 7, height = 3.5, dpi = 300)
