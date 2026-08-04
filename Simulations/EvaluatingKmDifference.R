library(ggplot2)
library(ggh4x)
library(dplyr)
library(RColorBrewer)

# Type 1 and 2 error ----------------

allEstimates <- readRDS("Simulations/allEstimates.rds")

errorStats <- allEstimates |>
  filter(model == "Kaplan Meier", contrast == "difference") |>
  mutate(h0 = 0) |>
  mutate(error = if_else(trueEffectType == "null",
                         lb > h0 | ub < h0,
                         lb <= h0 & ub >= h0)) |>
  group_by(timePoint,
           estimand,
           model,
           contrast,
           ciMethod,
           trueEffectType,
           depletionOfSusceptibles,
           sampleSize) |>
  summarise(error = mean(error, na.rm = TRUE), .groups = "drop") |>
  mutate(errorType = if_else(trueEffectType == "null", "type 1", "type 2"))

vizData <- errorStats |>
  mutate(`CI Method` = ciMethod,
         trueEffectType = paste("True effect:", SqlRender::camelCaseToTitleCase(trueEffectType), sep = "\n"),
         depletionOfSusceptibles = case_when(depletionOfSusceptibles == "outcome" ~ "Depletion of outcome susceptibles", 
                                             depletionOfSusceptibles == "effect" ~ "Depletion of effect susceptibles", 
                                             .default = "No depletion of susceptibles"),
         sampleSize = if_else(sampleSize == "big", "2,500 patients", "1,250 patients"))
optimal <- tibble(
  trueEffectType = c("True effect:\nNull", "True effect:\nMultiplicative", "True effect:\nAdditive"),
  y = c(0.05, 0, 0)
)

ggplot(vizData, aes(x = timePoint, y = error, group = `CI Method`, color = `CI Method`)) +
  geom_hline(aes(yintercept = y), data = optimal) +
  geom_line(linewidth = 0.75, alpha = 0.65) +
  scale_y_continuous("Error (type 1 or 2)") +
  scale_x_log10("Cutoff time (days)", breaks = unique(vizData$timePoint)) +
  scale_color_manual(values = brewer.pal(5, "Set1")) +
  facet_nested(trueEffectType ~ depletionOfSusceptibles + sampleSize) +
  theme(
    axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.5, color = "white"),
    legend.position = "top"
  )
ggsave("Simulations/plots/KmDifferenceType1And2Error.png", width = 8, height = 4.5, dpi = 300)
ggsave("Simulations/plots/KmDifferenceType1And2Error.svg", width = 8, height = 4.5)

# Estimates over time -------------------------------------------
source("Simulations/SimulationFunctions.R")

settings <- createSimulationSettings(
  pEffectSusceptible = 1,
  pOutcomeSusceptible = 0.2,
  logHrFunction = function(t) {25 * dweibull(t, 1.5, 10)}
)

population <- simulatePopulation(settings, seed = 123)

# Compute true counterfactual:
settingsCounterFactual <- settings
settingsCounterFactual$logHrFunction <- NULL
reference <- list()
for (i in seq_len(100)) {
  populationFactual <- simulatePopulation(settings, seed = i)
  populationCounterFactual <- simulatePopulation(settingsCounterFactual, seed = i)
  difference <- tibble(x = 1:365, y = NA, iteration = i)
  for (t in 1:365) {
    exposedOutcomes <- populationFactual |>
      filter(a == 1, survivalTime <= t) |>
      summarise(sum(y)) |>
      pull()
    exposedOutcomesCounterFactual <- populationCounterFactual |>
      filter(a == 1, survivalTime <= t) |>
      summarise(sum(y)) |>
      pull()
    difference$y[t] <- (exposedOutcomes - exposedOutcomesCounterFactual) / settings$n
  }
  reference[[i]] <- difference
}
reference <- bind_rows(reference)

# Risk ratio over time
source("Common/FunctionsForEstimands.R")
x <- 0:365
timePoints <- seq(from = 1, to = 365, by = 7)

cluster <- ParallelLogger::makeCluster(10)
estimates <- computeEstimands(population, timePoints = timePoints, bootstrapSize = 200, cluster = cluster)
ParallelLogger::stopCluster(cluster)

vizData <- estimates |>
  filter(model == "Kaplan Meier", contrast == "difference") |>
  transmute(x = timePoint,
            y = estimate,
            ymin = lb,
            ymax = ub,
            label = ciMethod) 

# vizData <- bind_rows(
#   tibble(
#     x = x,
#     y = c(0, difference),
#     label = "Observed outcomes -\nsimulated counterfactual"
#   ),
#   vizData
# )
# vizData$label <- factor(vizData$label, levels = c(
#   "Observed outcomes -\nsimulated counterfactual",
#   "asymptotic",
#   "percentile"))

colors <- brewer.pal(3, "Set1")
ggplot(vizData, aes(x = x, y = y)) +
  geom_hline(yintercept = 0, color = "darkgray") +
  geom_line(aes(group = iteration), color = "black", linewidth = 0.75, alpha = 0.1, data = reference) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax, color = label, fill = label, group = label), alpha = 0.1, size = 0) +
  geom_line(aes( color = label, group = label), linewidth = 0.75, alpha = 0.7) +
  scale_x_continuous("Time (days)") +
  scale_y_continuous("Risk difference") +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  coord_cartesian(xlim = c(0, 365)) +
  theme(
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = "right",
    legend.key.size = unit(2,"line")
  )

ggsave(filename = "Simulations/plots/KmEstimateOverTime.png", width = 7, height = 4, dpi = 300)
ggsave(filename = "Simulations/plots/KmEstimateOverTime.svg", width = 7, height = 4)

# Explore why KM performs worse for longer time cutoffs --------------------------------------------
library(ggplot2)
library(dplyr)
allEstimates <- readRDS("Simulations/allEstimates.rds")
vizData <- allEstimates |>
  filter(trueEffectType == "multiplicative",
         depletionOfSusceptibles == "none",
         sampleSize == "big",
         (ciMethod == "asymptotic" & model == "Cox") & contrast == "ratio" | (ciMethod == "asymptotic" & estimand == "Risk Difference")) |>
  filter(seed %in% 1:5) |>
  mutate(seed = paste("Seed", seed))

reference <- tibble(estimand = c("Hazard Ratio", "Risk Difference"),
                    h0 = c(1, 0))
ggplot(vizData, aes(x = timePoint, y = estimate)) +
  geom_hline(aes(yintercept = h0), data = reference) +
  geom_errorbar(aes(ymin = lb, ymax = ub)) +
  geom_point() +
  scale_x_log10("Time point") +
  scale_y_continuous("Effect estimate") +
  facet_grid(estimand ~ seed, scale = "free_y")
ggsave("Simulations/plots/KmVsCoxEstimates.png", width = 7, height = 4)
ggsave("Simulations/plots/KmVsCoxEstimates.svg", width = 7, height = 4)
