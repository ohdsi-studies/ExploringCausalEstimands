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
print(sprintf("Hazard ratio = %0.2f (95%% CI: %0.2f - %0.2f)",
              hr,
              ci[1],
              ci[2]))

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
vizData$label <- factor(vizData$label, levels = c("Within susceptibles", 
                                                  "Average without depletion", 
                                                  "Average over target"))
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
source("Common/FunctionsForEstimands.R")
cluster <- ParallelLogger::makeCluster(10)
estimates <- computeEstimands(population, timePoints = 1:100, bootstrapSize = 200, cluster = cluster)
ParallelLogger::stopCluster(cluster)

estimates <- estimates |>
  filter(ciMethod == "asymptotic",
         contrast == "ratio",
         model != "Cox") |>
  transmute(x = timePoint,
            y = log(estimate),
            ymin = log(lb),
            ymax = log(ub),
            label = model) |>
  mutate(y = if_else(label == "AFT", -y, y),
         ymin = if_else(label == "AFT", -ymin, ymin),
         ymax = if_else(label == "AFT", -ymax, ymax))

coxEstimate <- tibble(
  x = c(1, 100),
  y = rep(log(hr), 2),
  ymin = rep(log(ci[1]), 2),
  ymax = rep(log(ci[2]), 2),
  label = "Cox"
)

vizData <- bind_rows(
  tibble(
    x = x,
    y = settings$logHrFunction(x),
    label = "Within susceptibles",
    lineType = "solid"
  ),
  tibble(
    x = x,
    y = settings$logHrFunction(x) * c(settings$pSusceptible, targetSusceptiblesOverTime[1:100] / targetOverTime[1:100]),
    label = "Average over target",
    lineType = "dashed"
  ),
  tibble(
    x = x,
    y = settings$logHrFunction(x) * settings$pSusceptible,
    label = "Average without depletion",
    lineType = "dotted"
  ),
  estimates,
  coxEstimate
)

vizData$label <- factor(vizData$label, levels = c("Within susceptibles", 
                                                  "Average without depletion", 
                                                  "Average over target",
                                                  "AFT",
                                                  "Cox",
                                                  "Kaplan Meier",
                                                  "RMST"))
colors <- c("#000000", "#000000", "#000000", "#EB6622", "#11A08A", "#FBC511", "#69AED5", "#336B91")
ggplot(vizData, aes(x = x, y = y, color = label, fill = label, group = label)) +
  geom_hline(yintercept = 0) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2, size = 0) +
  geom_line(aes(linetype = label), linewidth = 1, alpha = 0.7) +
  scale_x_continuous("Time (days)") +
  scale_y_continuous("Hazard Ratio", 
                     breaks = log(c(1, 2, 3, 4, 5, 6, 7)), 
                     labels = c(1, 2, 3, 4, 5, 6, 7),
                     sec.axis = sec_axis(transform = ~., 
                                         breaks = log(c(1, 2, 3, 4, 5, 6, 7)), 
                                         labels = c(1, 2, 3, 4, 5, 6, 7),
                                         name = "Risk ratio")) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "solid", "solid", "solid", "solid")) +
  coord_cartesian(xlim = c(0, 100), ylim = c(log(1), log(7))) +
  theme(
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = "right"
  )

ggsave(filename = "Simulations/HrsAndRrs.png", width = 7, height = 3.5, dpi = 300)
