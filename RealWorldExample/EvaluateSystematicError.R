library(dplyr)
library(ggplot2)
library(RColorBrewer)

mdrrTheshold <- 2

mdrr <- readRDS("RealWorldExample/mdrr.rds")
estimatesMatched <- readRDS("RealWorldExample/estimatesMatched.rds")
negativeControls <- readr::read_csv("RealWorldExample/NegativeControls.csv", show_col_types = FALSE)


# Restrict to TCOs having sufficient power -------------------------------------
valid <- mdrr |>
  filter(mdrr < mdrrTheshold) |>
  select("targetId", "comparatorId", "outcomeId")

estimates <- estimatesMatched |>
  inner_join(valid, by = join_by("targetId", "comparatorId", "outcomeId")) |>
  mutate(adjustment = "PS matching")

# Plot type 1 error --------------------------------------------------------------------------------
computeType1Error <- function(lb, ub, h0 = 1) {
  rejectNull <- (!is.na(lb) & h0 < lb) | (!is.na(ub) & h0 > ub)
  return(mean(rejectNull))
}

type1Error <- estimates |>
  mutate(h0 = if_else(contrast == "ratio", 1, 0)) |>
  filter(outcomeId %in% negativeControls$conceptId) |>
  group_by(targetId, comparatorId, timePoint, estimand, model, contrast, ciMethod, adjustment) |>
  summarise(
    count = n(),
    type1Error = computeType1Error(lb, ub, h0),
    .groups = "drop"
  )

vizData <- type1Error |>
  filter(ciMethod == "asymptotic") |>
  mutate(Contrast = SqlRender::camelCaseToTitleCase(contrast),
         Model = model,
         example = if_else(targetId == 739138, 
                           "Sertraline vs bupropion for suicide attempt or ideation", 
                           "Lisinopril vs hydrochlorothiazide for angioedema"))

myPalette = brewer.pal(4, "Dark2")
plot1 <- ggplot(vizData, aes(x = timePoint, y = type1Error, group = estimand, color = Model)) +
  geom_hline(yintercept = 0.05) +
  geom_line(aes(linetype = Contrast), linewidth = 0.75, alpha = 0.75) +
  scale_y_continuous("Type 1 error") +
  scale_x_log10("Cutoff time (days)", breaks = unique(vizData$timePoint)) +
  scale_color_manual(values = myPalette) +
  facet_grid(~example) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.5, color = "white"),
    legend.position = "top",
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()
  )
# ggsave("RealWorldExample/Type1Error.png", width = 7, height = 3, dpi = 300)
# ggsave("RealWorldExample/Type1Error.svg", width = 7, height = 3, dpi = 300)

# Plot positive control p-value --------------------------------------------------------------------
# estimate = pValues$estimate; lb = pValues$lb; ub = pValues$ub; contrast = pValues$contrast
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

vizData <- estimates |>
  filter(ciMethod == "asymptotic", (targetId == 739138 & outcomeId == 1) | (targetId != 739138 & outcomeId == 2)) |>
  mutate(logP = computeLogP(estimate, lb, ub, contrast, model)) |>
  mutate(Contrast = SqlRender::camelCaseToTitleCase(contrast),
         Model = model,
         logP = pmax(if_else(is.na(logP), 0, logP), -200),
         example = if_else(targetId == 739138, "Sertraline vs bupropion for suicide attempt or ideation", "Lisinopril vs hydrochlorothiazide for angioedema"))

# Remove AFT from negative controls because it looks awful:
vizDataNcs <- estimates |>
  filter(ciMethod == "asymptotic", outcomeId > 10, model != "AFT") |>
  mutate(logP = computeLogP(estimate, lb, ub, contrast, model)) |>
  mutate(Contrast = SqlRender::camelCaseToTitleCase(contrast),
         Model = model,
         logP = pmax(if_else(is.na(logP), 0, logP), -200),
         example = if_else(targetId == 739138, "Sertraline vs bupropion for suicide attempt or ideation", "Lisinopril vs hydrochlorothiazide for angioedema"))


myPalette = brewer.pal(4, "Dark2")
breaks <- c(0.05, 1e-8, 1e-16, 1e-32, 1e-64)
labels <- c(0.5, parse(text = "10^-8"), parse(text = "10^-16"), parse(text = "10^-32"), parse(text = "10^-64"))

plot2 <- ggplot(vizData, aes(x = timePoint, y = logP, group = estimand, color = Model)) +
  geom_line(aes(group = paste(estimand, outcomeId)), color = "#666666", alpha = 0.1, data = vizDataNcs) +
  geom_hline(yintercept = log(0.05)) +
  geom_line(aes(linetype = Contrast), linewidth = 0.75, alpha = 0.75) +
  scale_y_continuous("One-sided P-value", breaks = log(breaks), labels = labels) +
  scale_x_log10("Cutoff time (days)", breaks = unique(vizData$timePoint)) +
  scale_color_manual(values = myPalette) +
  coord_cartesian(ylim = c(min(vizData$logP), max(vizData$logP))) +
  facet_grid(~example) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.5, color = "white"),
    legend.position = "none",
    strip.text = element_blank()
  )


# ggsave("RealWorldExample/PvalueOutcomeOfInterest.png", width = 7, height = 3, dpi = 300)


library(patchwork)
plot1 / plot2
ggsave("RealWorldExample/ErrorAndP.png", width = 7, height = 6, dpi = 300)
ggsave("RealWorldExample/ErrorAndP.svg", width = 7, height = 6)

# Plot estimates ---------------------------------------------------------------
plotEstimatesPerTimePoint <- function(estimatesSubset, title, fileName) {
  timePoints <- estimatesSubset |>
    pull(timePoint) |>
    unique() |>
    sort()
  vizData <- estimatesSubset |>
    filter(se < 10) |>
    mutate(estimand = gsub("Risk", "KM Risk", estimand)) |>
    mutate(h0 = if_else(contrast == "ratio", 1, 0),
           control = if_else(outcomeId %in% negativeControls$conceptId, "Negative control", "Positive control"),
           logEstimate = if_else(contrast == "ratio", log(estimate), estimate),
           estimand = if_else(contrast == "ratio", paste("Log", estimand, sep = "\n"), estimand),
           timePoint = factor(paste(timePoint, "days"), levels = paste(timePoints, "days")),
           logEstimate = if_else(model == "AFT", -logEstimate, logEstimate)) |>
    filter(!is.infinite(logEstimate))
  seScaleFactors <- vizData |>
    # group_by(timePoint, estimand) |>
    # summarise(meanSe = median(se, na.rm = TRUE), .groups ="drop") |>
    group_by(estimand) |>
    summarise(meanSe = sd(logEstimate, na.rm = TRUE), .groups ="drop") |>
    mutate(scaleFactor = 1 / meanSe)
  vizData <- vizData |>
    # inner_join(seScaleFactors, by = join_by(estimand, timePoint)) |>
    inner_join(seScaleFactors, by = join_by(estimand)) |>
    mutate(scaledSe = se * scaleFactor,
           scaledEstimate = logEstimate * scaleFactor) |>
    arrange(control)
  # estimateLabels <- vizData |>
  #   filter(outcomeId < 10) |>
  #   mutate(text = sprintf("%0.2f (%0.2f - %0.2f)", estimate, lb, ub))
  plot <- ggplot(vizData, aes(x = scaledEstimate, y = scaledSe, color = control, fill = control, shape = control, size = control, alpha = control)) +
    geom_hline(yintercept = 0) +
    geom_abline(aes(intercept = 0, slope = 1/(qnorm(0.025))), linetype = "dashed", data = seScaleFactors) +
    geom_abline(aes(intercept = 0, slope = 1/(qnorm(0.975))), linetype = "dashed", data = seScaleFactors) +
    geom_vline(xintercept = 0) +
    geom_point() +
    ggtitle(title) +
    # geom_label(aes(label = text), x = 0, y = 4, fill = "white", alpha = 0.8, hjust = 0.5, data = estimateLabels) +
    scale_x_continuous("Effect Size Estimate / SD(estimates)", limits = c(-6, 6)) +
    scale_y_continuous("Standard Error / SD(estimates)", limits = c(0, 4)) +
    scale_shape_manual(values = c(16, 23)) +
    scale_size_manual(values = c(2, 4)) +
    scale_color_manual(values = c(rgb(0, 0, 0.8),rgb(0, 0, 0))) +
    scale_fill_manual(values = c(rgb(0, 0, 0.8), rgb(1, 1, 0))) +
    scale_alpha_manual(values = c(0.5, 0.8)) +
    facet_grid(timePoint ~ estimand) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
    )
  ggsave(fileName, plot, width = 9, height = 7, dpi = 300)
  return(plot)
}

estimatesSubset <- estimates |>
  filter(ciMethod == "asymptotic", targetId != 739138, outcomeId > 10 | outcomeId == 2)
plotEstimatesPerTimePoint(estimatesSubset, 
                          title = "Lisinopril vs hydrochlorothiazide for angioedema",
                          fileName = "RealWorldExample/Example1.png")
estimatesSubset <- estimates |>
  filter(ciMethod == "asymptotic", targetId == 739138, outcomeId > 10 | outcomeId == 1)
plotEstimatesPerTimePoint(estimatesSubset, 
                          title = "Sertraline vs bupropion for suicide attempt or ideation",
                          fileName = "RealWorldExample/Example2.png")
