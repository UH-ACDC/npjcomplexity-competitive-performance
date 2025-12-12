library(readr)
library(dplyr)
library(e1071)  # for skewness()

if (interactive() &&
    requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

# 1. Read the GBR swimmers totals CSV
gbr_swim <- read_csv("ATHLETES_GBR_SWIMMERS_total.csv",
                     show_col_types = FALSE)

# If you ever wanted to restrict to athletes with >= 2 total medals, you could:
# gbr_swim <- gbr_swim %>% filter(Total_Score >= 2)

# 2. Compute descriptive stats on Total_Score
stats_gbr <- gbr_swim %>%
  summarise(
    n        = n(),
    Min      = min(Total_Score, na.rm = TRUE),
    Q1       = as.numeric(quantile(Total_Score, 0.25, na.rm = TRUE)),
    Median   = median(Total_Score, na.rm = TRUE),
    Mean     = mean(Total_Score, na.rm = TRUE),
    Q3       = as.numeric(quantile(Total_Score, 0.75, na.rm = TRUE)),
    Max      = max(Total_Score, na.rm = TRUE),
    SD       = sd(Total_Score, na.rm = TRUE),
    Skewness = skewness(Total_Score, type = 3, na.rm = TRUE)
  )

stats_gbr

# 3. Build LaTeX row for your table
latex_row <- with(stats_gbr, sprintf(
  "GBR Olympic swimmers & %d & %d & %.2f & %d & %.2f & %d & %.2f & %.2f & %.2f \\\\",
  n,
  Min,
  Q1,
  Median,
  Mean,
  Q3,
  Max,
  SD,
  Skewness
))

cat(latex_row, "\n")