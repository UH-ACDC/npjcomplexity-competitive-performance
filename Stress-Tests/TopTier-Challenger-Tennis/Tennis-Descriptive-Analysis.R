library(readr)
library(dplyr)
library(e1071)  # for skewness()

if (interactive() &&
    requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

# 1. Read the consolidated CSV
tennis_runs <- read_csv(
  "tennis_player_event_runs_knockout_2016_2017_2018_2019_2021_2022_2023_fromRAW.csv",
  show_col_types = FALSE
)

# Sanity check: only "Top" and "Challenger"
table(tennis_runs$stratum)

# 1b. Restrict to events with at least 2 wins
tennis_runs_ge2 <- tennis_runs %>%
  filter(wins_in_event >= 2)

# Optional sanity check
table(tennis_runs_ge2$stratum)

# 2. Compute stats by stratum on wins_in_event >= 2 only
stats_tennis <- tennis_runs_ge2 %>%
  group_by(stratum) %>%
  summarise(
    n        = n(),
    Min      = min(wins_in_event, na.rm = TRUE),
    Q1       = as.numeric(quantile(wins_in_event, 0.25, na.rm = TRUE)),
    Median   = median(wins_in_event, na.rm = TRUE),
    Mean     = mean(wins_in_event, na.rm = TRUE),
    Q3       = as.numeric(quantile(wins_in_event, 0.75, na.rm = TRUE)),
    Max      = max(wins_in_event, na.rm = TRUE),
    SD       = sd(wins_in_event, na.rm = TRUE),
    Skewness = skewness(wins_in_event, type = 3, na.rm = TRUE),
    .groups  = "drop"
  )

stats_tennis

# 3. Build LaTeX rows (Top first)
stats_tennis_latex <- stats_tennis %>%
  arrange(factor(stratum, levels = c("Top", "Challenger"))) %>%
  mutate(
    Dataset = if_else(stratum == "Top",
                      "Tennis Top-tier runs",
                      "Tennis Challenger runs"),
    latex_row = sprintf(
      "%s & %d & %d & %.2f & %d & %.2f & %d & %.2f & %.2f & %.2f \\\\",
      Dataset,
      n,
      Min,
      Q1,
      Median,
      Mean,
      Q3,
      Max,
      SD,
      Skewness
    )
  )

cat(stats_tennis_latex$latex_row, sep = "\n")