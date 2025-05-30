# Code to generate plots for comparing different estimands
source("Simulations/SimulationFunctions.R")

library(survival)
library(ggplot2)
library(dplyr)
library(Cyclops)

settings <- createSimulationSettings()

# Plots of the generative process -------------------------------
x <- 0:100
vizData <- tibble(
  x = x,
  y = settings$baselineHazardFunction(x)
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
    y = settings$logHrFunction(x),
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
population <- simulatePopulation(settings, seed = 123)

cyclopsData <- createCyclopsData(Surv(survivalTime, y) ~ a, modelType = "cox", data = population)
fit <- fitCyclopsModel(cyclopsData)
hr <- exp(coef(fit))
ci <- exp(confint(fit, parm = "a")[c(2, 3)])
# fit <- coxph(Surv(survivalTime, y) ~ a, data = population)
# ci <- exp(confint(fit))
print(sprintf("Hazard ratio = %0.2f (95%% CI: %0.2f - %0.2f)",
              hr,
              ci[1],
              ci[2]))
# [1] "Hazard ratio = 1.27 (95% CI: 1.13 - 1.42)"

# Plots of estimates -------------------------------------------------------------------------------

# Depletion of susceptibles
targetSusceptiblesOverTime <- attr(population, "targetSusceptiblesOverTime")
targetOverTime <- attr(population, "targetOverTime")
vizData <- tibble(
  x = x,
  y = c(settings$pSusceptible, targetSusceptiblesOverTime[1:100] / targetOverTime[1:100])
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
    y = settings$logHrFunction(x),
    label = "Within susceptibles"
  ),
  tibble(
    x = x,
    y = settings$logHrFunction(x) * c(settings$pSusceptible, targetSusceptiblesOverTime[1:100] / targetOverTime[1:100]),
    label = "Average over target"
  ),
  tibble(
    x = x,
    y = settings$logHrFunction(x) * settings$pSusceptible,
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
km <- survfit(Surv(survivalTime, y) ~ a, data = population)
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
source("Common/FunctionsForNonHrEstimands.R")
cluster <- ParallelLogger::makeCluster(10)
kmEstimate <-computeEstimands(population, timePoints = 1:100, bootstrapSize = 200, cluster = cluster)
ParallelLogger::stopCluster(cluster)

vizData <- bind_rows(
  tibble(
    x = x,
    y = settings$logHrFunction(x),
    label = "Within susceptibles"
  ),
  tibble(
    x = x,
    y = settings$logHrFunction(x) * c(settings$pSusceptible, targetSusceptiblesOverTime[1:100] / targetOverTime[1:100]),
    label = "Average over target"
  ),
  tibble(
    x = x,
    y = settings$logHrFunction(x) * settings$pSusceptible,
    label = "Average without depletion"
  )
)

kmEstimate <- kmEstimate |>
  select(x = timePoint, riskRatio = rr, lower = lbRrAsymptotics, upper = ubRrAsymptotics) |>
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
