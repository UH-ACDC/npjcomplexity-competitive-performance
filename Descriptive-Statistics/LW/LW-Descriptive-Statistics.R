library(readr)
library(dplyr)
library(e1071)  # for skewness() and kurtosis()

if (interactive() &&
    requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

# 1. Read pilots file
lw <- read_csv("GER_pilots.csv", show_col_types = FALSE)

# 2. Compute entry year and total victories
lw <- lw %>%
  mutate(
    entry_year = case_when(
      y1939 > 0 ~ 1939L,
      y1940 > 0 ~ 1940L,
      y1941 > 0 ~ 1941L,
      y1942 > 0 ~ 1942L,
      y1943 > 0 ~ 1943L,
      y1944 > 0 ~ 1944L,
      y1945 > 0 ~ 1945L,
      TRUE       ~ NA_integer_
    ),
    total_victories = y1939 + y1940 + y1941 + y1942 + y1943 + y1944 + y1945
  )

# 3. Restrict to 1940–1942 entry cohorts
lw_404142 <- lw %>%
  filter(entry_year %in% c(1940L, 1941L, 1942L))

# 4. Descriptive stats in the order of your table header
stats_LW <- lw_404142 %>%
  summarise(
    n       = n(),
    min_val = min(total_victories),
    q1      = as.numeric(quantile(total_victories, 0.25, type = 7)),
    median  = median(total_victories),
    mean    = mean(total_victories),
    q3      = as.numeric(quantile(total_victories, 0.75, type = 7)),
    max_val = max(total_victories),
    sd_val  = sd(total_victories),
    skew    = e1071::skewness(total_victories, type = 2)
  )

# 5. LaTeX row matching: Dataset, n, Min, Q1, Median, Mean, Q3, Max, SD, Skewness
latex_row_LW <- with(
  stats_LW,
  sprintf(
    "LW & %d & %d & %d & %d & %.2f & %d & %d & %.2f & %.2f \\\\",
    n,
    min_val,
    q1,
    median,
    mean,
    q3,
    max_val,
    sd_val,
    skew
  )
)

cat(latex_row_LW, "\n")