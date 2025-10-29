library(survival)
library(ggplot2)
library(dplyr)
library(Cyclops)
source("Simulations/SimulationFunctions.R")

settings <- createSimulationSettings()

# Hazard ratio -------------------------------------------------------------------------------------
x <- 0:100
vizData <- tibble(
  y = c(settings$baselineHazardFunction(x), settings$baselineHazardFunction(x) * 2),
  label = c(rep("Comparator", length(x)), rep("Target", length(x))),
  x = c(x, x)
)
ggplot(vizData, aes(x = x, y = y, color = label, group = label)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, color = "darkgray") +
  scale_x_continuous("Time (days)") +
  scale_y_continuous("Hazard") +
  scale_color_manual(values = c("#336B91", "#EB6622")) +
  theme(
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = "bottom"
  )
ggsave(filename = "Simulations/IllustrationHr.svg", width = 5, height = 3.5, dpi = 300)


fit <- survfit(Surv(survivalTime, y) ~ treatment,
               robust = FALSE,
               data = sampledData)