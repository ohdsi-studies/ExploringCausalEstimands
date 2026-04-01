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
vizData1 <- tibble(
  x = x,
  y = settings$baselineHazardFunction(x),
  type = "Baseline hazard"
)
hrFunction <- function(t) {exp(5 * dweibull(t, 1.5, 10))}
vizData2 <- bind_rows(
  tibble(
    x = x,
    y = hrFunction(x),
    type = "Hazard ratio\n(Multiplicative effect)",
    label = "Within susceptibles"
  )
)
rdFunction <- function(t) {0.2 * dweibull(t, 1, 10)}
vizData3 <- bind_rows(
  tibble(
    x = x,
    y = rdFunction(x),
    type = "Risk difference\n(Additive effect)",
    label = "Within susceptibles"
  )
)
vizData <- bind_rows(vizData1, vizData2, vizData3)
reference <- tibble(
  y = c(0, 1, 0),
  type = c("Baseline hazard", 
           "Hazard ratio\n(Multiplicative effect)",
           "Risk difference\n(Additive effect)")
)
vizData$type <- factor(vizData$type, levels = c("Baseline hazard", 
                                                "Hazard ratio\n(Multiplicative effect)",
                                                "Risk difference\n(Additive effect)"))
reference$type <- factor(reference$type, levels = c("Baseline hazard", 
                                                "Hazard ratio\n(Multiplicative effect)",
                                                "Risk difference\n(Additive effect)"))
ggplot(vizData, aes(x = x, y = y)) +
  geom_hline(aes(yintercept = y), color = "darkgray", data = reference) +
  geom_line(linewidth = 1, alpha = 0.7) +
  scale_x_continuous("Time (days)") +
  facet_grid(type ~ ., scales = "free_y", switch = "y") +
  theme(
    panel.grid.minor = element_blank(),
    axis.title.y = element_blank(),
    strip.placement = "outside",
    strip.background = element_blank()
  )
ggsave(filename = "Simulations/plots/GenerativeFunctions.png", width = 4, height = 4, dpi = 300)
ggsave(filename = "Simulations/plots/GenerativeFunctions.svg", width = 4, height = 4)
ggsave(filename = "Simulations/plots/GenerativeFunctions.pdf", width = 4, height = 4)
