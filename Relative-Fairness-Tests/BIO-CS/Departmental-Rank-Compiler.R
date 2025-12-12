library(readr)
library(dplyr)

if (interactive() &&
    requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

## 1. Read the faculty file (only to grab SchoolRank)
faculty <- read_csv("Faculty_GoogleScholar_Funding_Data_N4190.csv",
                    show_col_types = FALSE)

rank_lookup <- faculty %>%
  select(google_id, dept, SchoolRank)

## 2. Read your trusted breakdown files
bio_break <- read_csv("SCHOLARS_BIO_breakdown.csv",
                      show_col_types = FALSE)
cs_break  <- read_csv("SCHOLARS_CS_breakdown.csv",
                      show_col_types = FALSE)

## 3. Attach SchoolRank WITHOUT changing any existing columns
bio_break_ranked <- bio_break %>%
  left_join(rank_lookup, by = c("google_id", "dept"))

cs_break_ranked <- cs_break %>%
  left_join(rank_lookup, by = c("google_id", "dept"))

## Optional sanity checks
stopifnot(nrow(bio_break_ranked) == nrow(bio_break))
stopifnot(nrow(cs_break_ranked)  == nrow(cs_break))

# How many rows didn’t find a rank? (these are the #NAME? placeholders)
sum(is.na(bio_break_ranked$SchoolRank))  # 34
sum(is.na(cs_break_ranked$SchoolRank))   # 22

## 4. Write the new breakdown files (same data + 1 extra column)
write_csv(bio_break_ranked, "SCHOLARS_BIO_breakdown_withRank.csv")
write_csv(cs_break_ranked,  "SCHOLARS_CS_breakdown_withRank.csv")