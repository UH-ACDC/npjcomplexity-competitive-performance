library(readr)
library(dplyr)
library(e1071)  # for skewness()

# Optional: set working directory to script location in RStudio
if (interactive() &&
    requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

# Helper: compute stats + LaTeX row for one breakdown file
compute_totals_row <- function(df, dataset_label, min_hits = 2) {
  # Detect medal columns of the form XYYYY.Summer / XYYYY.Winter
  year_cols <- grep("^X[0-9]{4}\\.", names(df), value = TRUE)
  if (length(year_cols) == 0) {
    stop(sprintf("No year columns of the form XYYYY.<season> found for %s", dataset_label))
  }
  
  df_totals <- df %>%
    mutate(
      total_hits = rowSums(across(all_of(year_cols)), na.rm = TRUE)
    ) %>%
    # Keep only athletes with at least `min_hits` medals
    filter(total_hits >= min_hits)
  
  stats <- df_totals %>%
    summarise(
      n       = n(),
      min_val = min(total_hits),
      q1      = as.numeric(quantile(total_hits, 0.25, type = 7)),
      median  = median(total_hits),
      mean    = mean(total_hits),
      q3      = as.numeric(quantile(total_hits, 0.75, type = 7)),
      max_val = max(total_hits),
      sd_val  = sd(total_hits),
      skew    = e1071::skewness(total_hits, type = 2)
    )
  
  # Use %.0f for integer-like quantities, %.2f for means/SD/skewness
  latex_row <- with(
    stats,
    sprintf(
      "%s & %.0f & %.0f & %.0f & %.0f & %.2f & %.0f & %.0f & %.2f & %.2f \\\\",
      dataset_label,
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
  
  cat(latex_row, "\n")
}

# Read the breakdown files
fra <- read_csv("ATHLETES_FRA_FENCING_breakdown.csv", show_col_types = FALSE)
usa <- read_csv("ATHLETES_USA_SWIMMERS_breakdown.csv", show_col_types = FALSE)

# FRA fencing (FRAOT_F) — athletes with 2+ medals
compute_totals_row(fra, "FRAOT\\_F", min_hits = 2)

# USA swimmers (USOT_S) — athletes with 2+ medals
compute_totals_row(usa, "USOT\\_S", min_hits = 2)