# Code to generate plots for comparing different estimands
source("Simulations/SimulationFunctions.R")

library(survival)
library(ggplot2)
library(dplyr)
library(Cyclops)
library(RColorBrewer)

settings <- createSimulationSettings()

# Plots of the generative process -------------------------------
x <- 0:100
vizData <- tibble(
  x = x,
  y = settings$baselineHazardFunction(x)
)
ggplot(vizData, aes(x = x, y = y)) +
  geom_line(color = "black", linewidth = 1) +
  geom_hline(yintercept = 0, color = "darkgray") +
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
  geom_hline(yintercept = 0, color = "darkgray") +
  geom_line(linewidth = 1, alpha = 0.7) +
  scale_x_continuous("Time (days)") +
  scale_y_continuous("Hazard Ratio", breaks = log(c(1, 2, 3, 4, 5, 6, 7)), labels = c(1, 2, 3, 4, 5, 6,7 )) +
  scale_color_manual(values = "black") +
  theme(
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = "right"
  )
ggsave(filename = "Simulations/HazardRatio.png", width = 7, height = 3.5, dpi = 300)


rdFunction <- function(t) {dweibull(t, 1, 10)}
vizData <- bind_rows(
  tibble(
    x = x,
    y = rdFunction(x),
    label = "Within susceptibles"
  )
)
ggplot(vizData, aes(x = x, y = y, color = label, group = label)) +
  geom_hline(yintercept = 0, color = "darkgray") +
  geom_line(linewidth = 1, alpha = 0.7) +
  scale_x_continuous("Time (days)") +
  scale_y_continuous("Risk Difference") + 
  scale_color_manual(values = "black") +
  theme(
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = "right"
  )
ggsave(filename = "Simulations/RiskDifference.png", width = 7, height = 3.5, dpi = 300)


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
settingsCounterFactual <- settings
settingsCounterFactual$logHrFunction <- NULL
populationCounterFactual <- simulatePopulation(settingsCounterFactual, seed = 123)
ratio <- c()
for (t in 1:100) {
  exposedOutcomes <- population |>
    filter(a == 1, survivalTime <= t) |>
    summarise(sum(y)) |>
    pull()
  exposedOutcomesCounterFactual <- populationCounterFactual |>
    filter(a == 1, survivalTime <= t) |>
    summarise(sum(y)) |>
    pull()
  ratio[t] <- exposedOutcomes / exposedOutcomesCounterFactual
}

# Plots of estimates -------------------------------------------------------------------------------

# Depletion of susceptibles
targetEffectSusceptiblesOverTime <- attr(population, "targetEffectSusceptiblesOverTime")
targetOverTime <- attr(population, "targetOverTime")
vizData <- tibble(
  x = x,
  y = c(settings$pEffectSusceptible, targetEffectSusceptiblesOverTime[1:100] / targetOverTime[1:100])
)
ggplot(vizData, aes(x = x, y = y)) +
  geom_hline(yintercept = 0, color = "darkgray") +
  geom_line(linewidth = 1, color = "black", alpha = 0.7) +
  scale_x_continuous("Time (days)") +
  scale_y_continuous("Fraction of target susceptible") +
  theme(
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = "bottom"
  )
ggsave(filename = "Simulations/FractionSusceptible.png", width = 5, height = 3.5, dpi = 300)

# True effect given depletion of susceptibles
trueHazardRatioOverTime <- attr(population, "trueHazardRatioOverTime")

vizData <- bind_rows(
  tibble(
    x = x,
    y = settings$logHrFunction(x),
    label = "Within susceptibles"
  ),
  tibble(
    x = x,
    y = c(settings$logHrFunction(0), log(trueHazardRatioOverTime[1:100])),
    label = "Average over population at risk"
  ),
  tibble(
    x = x,
    y = log(exp(settings$logHrFunction(x)) * settings$pEffectSusceptible + (1 - settings$pEffectSusceptible)),
    label = "Average without depletion"
  ),
  tibble(
    x = x,
    y = c(0, log(ratio)),
    label = "Observed outcomes / counterfactual"
  )
)

vizData$label <- factor(vizData$label, levels = c("Within susceptibles", 
                                                  "Average without depletion", 
                                                  "Average over population at risk",
                                                  "Observed outcomes / counterfactual"))
ggplot(vizData, aes(x = x, y = y, linetype = label, group = label)) +
  geom_hline(yintercept = 0, color = "darkgray") +
  geom_line(linewidth = 1, alpha = 0.7) +
  scale_x_continuous("Time (days)") +
  scale_y_continuous("Hazard Ratio", breaks = log(c(1, 2, 3, 4, 5, 6, 7)), labels = c(1, 2, 3, 4, 5, 6,7 )) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "solid", "solid", "solid", "solid")) +
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
         contrast == "ratio") |>
         #model != "Cox") |>
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
    label = "Within susceptibles"
  ),
  tibble(
    x = x,
    y = c(settings$logHrFunction(0), log(trueHazardRatioOverTime[1:100])),
    label = "Average over population at risk"
  ),
  tibble(
    x = x,
    y = log(exp(settings$logHrFunction(x)) * settings$pEffectSusceptible + (1 - settings$pEffectSusceptible)),
    label = "Average without depletion"
  ),
  tibble(
    x = x,
    y = c(0, log(ratio)),
    label = "Observed outcomes / counterfactual"
  ),
  estimates
  # coxEstimate
)

vizData$label <- factor(vizData$label, levels = c("Within susceptibles", 
                                                  "Average without depletion", 
                                                  "Average over population at risk",
                                                  "Observed outcomes / counterfactual",
                                                  "AFT",
                                                  "Cox",
                                                  "Kaplan Meier",
                                                  "RMST"))
# colors <- c("#000000", "#000000", "#000000", "#EB6622", "#11A08A", "#FBC511", "#69AED5", "#336B91")
colors <- c("#000000", "#000000", "#000000", "#000000", brewer.pal(4, "Dark2"))
ggplot(vizData, aes(x = x, y = y, color = label, fill = label, group = label)) +
  geom_hline(yintercept = 0, color = "darkgray") +
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
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash", "solid", "solid", "solid", "solid")) +
  coord_cartesian(xlim = c(0, 100), ylim = c(log(1), log(7))) +
  theme(
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = "right"
  )

ggsave(filename = "Simulations/HrsAndRrs.png", width = 7, height = 3.5, dpi = 300)

# Some builds:
subset <- vizData |>
  filter(label %in% c("Within susceptibles", 
                      "Average without depletion", 
                      "Average over population at risk",
                      "Cox"))
subsetColors <- colors[c(1, 2, 3, 5)]
ggplot(subset, aes(x = x, y = y, color = label, fill = label, group = label)) +
  geom_hline(yintercept = 0, color = "darkgray") +
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
  scale_color_manual(values = subsetColors) +
  scale_fill_manual(values = subsetColors) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "solid", "solid", "solid", "solid")) +
  coord_cartesian(xlim = c(0, 100), ylim = c(log(1), log(7))) +
  theme(
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = "right"
  )

ggsave(filename = "Simulations/HrsAndRrsCoxOnly.png", width = 7, height = 3.5, dpi = 300)

subset <- vizData |>
  filter(label %in% c("Within susceptibles", 
                      "Average without depletion", 
                      "Average over population at risk",
                      "Cox",
                      "Kaplan Meier"))
subsetColors <- colors[c(1, 2, 3, 5, 6)]
ggplot(subset, aes(x = x, y = y, color = label, fill = label, group = label)) +
  geom_hline(yintercept = 0, color = "darkgray") +
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
  scale_color_manual(values = subsetColors) +
  scale_fill_manual(values = subsetColors) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "solid", "solid", "solid", "solid")) +
  coord_cartesian(xlim = c(0, 100), ylim = c(log(1), log(7))) +
  theme(
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = "right"
  )

ggsave(filename = "Simulations/HrsAndRrsCoxAndKmOnly.png", width = 7, height = 3.5, dpi = 300)

