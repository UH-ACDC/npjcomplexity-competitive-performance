library(readr)
library(dplyr)
library(e1071)  # for skewness()

if (interactive() &&
    requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

compute_scholar_row <- function(file, dataset_label,
                                first_year_max = 2010,
                                min_awards    = 2) {
  df <- read_csv(file, show_col_types = FALSE)
  
  # year columns of the form y1996, y1997, ...
  year_cols <- grep("^y[0-9]{4}$", names(df), value = TRUE)
  if (length(year_cols) == 0) {
    stop(sprintf("No yYYYY columns found in %s", file))
  }
  
  years <- as.integer(sub("^y", "", year_cols))
  
  # matrix of award counts per year
  year_mat <- as.matrix(df[, year_cols])
  
  # first year with > 0 awards for each scholar
  first_year <- apply(year_mat, 1, function(row) {
    idx <- which(row > 0)
    if (length(idx) == 0) return(NA_integer_)
    years[min(idx)]
  })
  
  df_totals <- df %>%
    mutate(
      total_awards = rowSums(across(all_of(year_cols)), na.rm = TRUE),
      first_year   = first_year
    ) %>%
    filter(
      !is.na(first_year),
      first_year <= first_year_max,
      total_awards >= min_awards
    )
  
  stats <- df_totals %>%
    summarise(
      n       = n(),
      min_val = min(total_awards),
      q1      = as.numeric(quantile(total_awards, 0.25, type = 7)),
      median  = median(total_awards),
      mean    = mean(total_awards),
      q3      = as.numeric(quantile(total_awards, 0.75, type = 7)),
      max_val = max(total_awards),
      sd_val  = sd(total_awards),
      skew    = e1071::skewness(total_awards, type = 2)
    )
  
  # Matches header:
  # Dataset, n, Min, Q1, Median, Mean, Q3, Max, SD, Skewness
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

# BIO scholars: first grant <= 2010, total awards >= 2
compute_scholar_row("SCHOLARS_BIO_breakdown.csv", "BIO")

# CS scholars: first grant <= 2010, total awards >= 2
compute_scholar_row("SCHOLARS_CS_breakdown.csv", "CS")