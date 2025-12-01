# Code for analyzing the real-world data results. Assumes ComputeEstimands.R has been executed.

library(dplyr)
library(ggplot2)
library(patchwork)
library(wesanderson)

estimatesMatched <- readRDS("RealWorldExample/estimatesMatched.rds")
negativeControls <- readr::read_csv("RealWorldExample/NegativeControls.csv", show_col_types = FALSE)

estimates <- estimatesMatched |>
  filter((ciMethod == "percentile" | model == "Cox" | model == "AFT")) |>
  mutate(model = if_else(estimand == "Cumulative Hazard Ratio", "Cumulative Hazard (KM)", model)) |>
  mutate(model = if_else(model == "Kaplan Meier", "Risk (KM)", model))

# Plot type 1 error --------------------------------------------------------------------------------
computeType1Error <- function(lb, ub, h0 = 1) {
  rejectNull <- (!is.na(lb) & h0 < lb) | (!is.na(ub) & h0 > ub)
  return(mean(rejectNull))
}

type1Error <- estimates |>
  mutate(h0 = if_else(contrast == "ratio", 1, 0)) |>
  filter(outcomeId %in% negativeControls$conceptId) |>
  group_by(targetId, comparatorId, timePoint, estimand, model, contrast, sampleSize) |>
  summarise(
    count = n(),
    type1Error = computeType1Error(lb, ub, h0),
    .groups = "drop"
  )

vizData <- type1Error |>
  mutate(Contrast = SqlRender::camelCaseToTitleCase(contrast),
         Model = model,
         sampleSize = sprintf("N=%s", format(sampleSize, big.mark = ",", scientific = FALSE)),
         example = if_else(targetId == 739138, 
                           "Sertraline vs bupropion for suicide attempt or ideation", 
                           "Lisinopril vs hydrochlorothiazide for angioedema"))


myPalette <- wes_palette("Darjeeling1", 5, type = "discrete")
plot1 <- ggplot(vizData, aes(x = timePoint, y = type1Error, group = estimand, color = Model)) +
  geom_hline(yintercept = 0.05) +
  geom_line(aes(linetype = Contrast), linewidth = 0.75, alpha = 0.75) +
  scale_y_continuous("Type 1 error") +
  scale_x_log10("Cutoff time (days)", breaks = unique(vizData$timePoint)) +
  scale_color_manual(values = myPalette) +
  facet_grid(sampleSize~example) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.5, color = "white"),
    legend.position = "top",
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()
  )
plot1
# ggsave("RealWorldExample/Type1Error.png", width = 8, height = 5, dpi = 300)
# ggsave("RealWorldExample/Type1Error.svg", width = 8, height = 5)

# Plot positive control p-value --------------------------------------------------------------------
computeLogP <- function(estimate, lb, ub, contrast, model) {
  estimate <- if_else(contrast == "ratio", log(estimate), estimate)
  lb <- if_else(contrast == "ratio", log(lb), lb)
  ub <- if_else(contrast == "ratio", log(ub), ub)
  se <- (ub - lb) / (2 * qnorm(0.975))
  estimate <- if_else(model == "AFT", -estimate, estimate)
  z <- estimate / se
  logP <- pnorm(z, log.p = TRUE, lower.tail = FALSE)
  return(logP)
}

estimate <- estimates |>
  filter(targetId == 739138, outcomeId == 1, sampleSize == 50000, model == "AFT")

vizData <- estimates |>
  filter((targetId == 739138 & outcomeId == 1) | (targetId != 739138 & outcomeId == 2)) |>
  mutate(logP = computeLogP(estimate, lb, ub, contrast, model)) |>
  mutate(Contrast = SqlRender::camelCaseToTitleCase(contrast),
         Model = model,
         sampleSize = sprintf("N=%s", format(sampleSize, big.mark = ",", scientific = FALSE)),
         logP = pmax(if_else(is.na(logP), 0, logP), -200),
         example = if_else(targetId == 739138, "Sertraline vs bupropion for suicide attempt or ideation", "Lisinopril vs hydrochlorothiazide for angioedema"))

myPalette <- wes_palette("Darjeeling1", 5, type = "discrete")
breaks <- c(0.5, 0.05, 0.01, 0.001)
labels <- breaks

plot2 <- ggplot(vizData, aes(x = timePoint, y = logP, group = estimand, color = Model)) +
  # geom_line(aes(group = paste(estimand, outcomeId)), color = "#666666", alpha = 0.1, data = vizDataNcs) +
  geom_hline(yintercept = log(0.05)) +
  geom_line(aes(linetype = Contrast), linewidth = 0.75, alpha = 0.75) +
  scale_y_continuous("One-sided P-value", breaks = log(breaks), labels = labels) +
  scale_x_log10("Cutoff time (days)", breaks = unique(vizData$timePoint)) +
  scale_color_manual(values = myPalette) +
  # coord_cartesian(ylim = c(min(vizData$logP), max(vizData$logP))) +
  coord_cartesian(ylim = c(log(0.0001), log(1))) +
  facet_grid(sampleSize~example) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.5, color = "white"),
    legend.position = "none",
    strip.text.x = element_blank()
  )
plot2
# ggsave("RealWorldExample/PvalueOutcomeOfInterest.png", width =8, height = 5, dpi = 300)

# Combine plots ------------------------------------------------------------------------------------
plot1 / plot2
ggsave("RealWorldExample/plots/ErrorAndP.png", width = 8, height = 8, dpi = 300)
ggsave("RealWorldExample/plots/ErrorAndP.svg", width = 8, height = 8)
