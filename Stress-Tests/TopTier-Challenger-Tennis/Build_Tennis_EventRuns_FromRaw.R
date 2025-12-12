# ================================================================
# Build_Tennis_EventRuns_FromRaw.R
#
# Build wins-per-player-per-event for:
#   - Top-tier ATP Tour (A/M/G)
#   - Challenger main draw (C)
# from Jeff Sackmann "tennis_atp-master" raw data.
#
# Output:
#   tennis_player_event_runs_knockout_2016_2017_2018_2019_2021_2022_2023_fromRAW.csv
#
# This should reproduce the structure and counts in
#   tennis_player_event_runs_knockout_2016_2017_2018_2019_2021_2022_2023_fromNULL.csv
# (the Top-tier part matches exactly, Challenger matches up to one minor RET edge case).
# ================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

if (interactive() &&
    requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

## ---------------- User knobs ----------------
YEARS        <- setdiff(2016:2023, 2020)  # match the revision: 2016–2019, 2021–2023
DRAW_MIN_TOP <- 32                        # Top-tier: events with draw_size >= 32
DRAW_MIN_CH  <- 16                        # Challenger: events with draw_size >= 16
EXCLUDE_WO   <- TRUE                      # drop W/O, RET, DEF matches (best effort)
DATA_ROOT    <- "tennis_atp-master"       # folder with atp_matches_*.csv etc.

set.seed(42)

message("Working directory: ", normalizePath(getwd()))
message("Data root:        ", normalizePath(DATA_ROOT))

## ---------------- Small helpers ----------------

`%||%` <- function(a, b) if (!is.null(a)) a else b

find_atp_csvs <- function(root, pattern) {
  list.files(root, pattern = pattern, recursive = TRUE, full.names = TRUE)
}

safe_read <- function(path) {
  suppressMessages(readr::read_csv(path, show_col_types = FALSE, progress = FALSE))
}

year_from_path <- function(p, pat) {
  m <- stringr::str_match(basename(p), pat)[, 2]
  suppressWarnings(as.integer(m))
}

is_main_draw_round <- function(r) {
  r <- as.character(r)
  grepl("^(F|SF|QF|R\\d+)$", r)
}

# We won’t actually use this here, but we keep it for clarity / completeness
is_qual_round <- function(r) {
  r <- as.character(r)
  grepl("^Q\\d+$", r)
}

drop_walkovers <- function(df) {
  if (!"score" %in% names(df)) return(df)
  df %>%
    filter(!str_detect(tolower(score %||% ""), "w/o|\\bwo\\b|ret|def"))
}

build_wins_per_event <- function(df, wins_col = "wins_in_event") {
  if (!nrow(df)) return(tibble())
  df %>%
    group_by(tourney_id, tourney_name, draw_size, year, winner_name) %>%
    summarise("{wins_col}" := n(), .groups = "drop")
}

summarise_block <- function(df, label, value_col) {
  if (!nrow(df)) {
    message("No rows for ", label)
    return(invisible(NULL))
  }
  v <- df[[value_col]]
  message("\n", label, ":")
  message("  N rows           = ", nrow(df))
  message("  years            = ", paste(sort(unique(df$year)), collapse = ", "))
  message("  ", value_col, " value counts:")
  print(table(v))
}

## ---------------- Locate raw files ----------------

tour_files <- find_atp_csvs(DATA_ROOT, "^atp_matches_\\d{4}\\.csv$")
qc_files   <- find_atp_csvs(DATA_ROOT, "^atp_matches_qual_chall_\\d{4}\\.csv$")

if (!length(tour_files) && !length(qc_files)) {
  stop("Could not find ATP CSVs under '", DATA_ROOT,
       "'. Check that tennis_atp-master/ is unzipped here.")
}

use_tour <- tour_files[year_from_path(tour_files, "atp_matches_(\\d{4})\\.csv") %in% YEARS]
use_qc   <- qc_files[year_from_path(qc_files, "atp_matches_qual_chall_(\\d{4})\\.csv") %in% YEARS]

message("Found ", length(use_tour), " ATP Tour files (A/M/G).")
message("Found ", length(use_qc),   " Qual/Chall files (C + Qual).")

if (!length(use_tour)) stop("No atp_matches_YYYY.csv files for YEARS.")
if (!length(use_qc))   stop("No atp_matches_qual_chall_YYYY.csv files for YEARS.")

## ---------------- Load & filter: Top-tier (A/M/G) ----------------

tour_raw <- bind_rows(lapply(use_tour, safe_read))
need_top <- c("tourney_id","tourney_name","tourney_date","tourney_level",
              "draw_size","winner_name","loser_name","score","round")
miss_top <- setdiff(need_top, names(tour_raw))
if (length(miss_top)) {
  stop("Top-tier raw matches missing columns: ", paste(miss_top, collapse = ", "))
}

tour_dat <- tour_raw %>%
  mutate(year = suppressWarnings(as.integer(substr(as.character(tourney_date), 1, 4)))) %>%
  filter(year %in% YEARS,
         tourney_level %in% c("A","M","G"),
         !is.na(draw_size), draw_size >= DRAW_MIN_TOP,
         is_main_draw_round(round))

if (EXCLUDE_WO) tour_dat <- drop_walkovers(tour_dat)

## ---------------- Load & filter: Challenger (C main draw) ----------------

qc_raw <- bind_rows(lapply(use_qc, safe_read))

need_qc <- c("tourney_id","tourney_name","tourney_date","tourney_level",
             "draw_size","winner_name","loser_name","score","round")
miss_qc <- setdiff(need_qc, names(qc_raw))
if (length(miss_qc)) {
  stop("Qual/Chall raw matches missing columns: ", paste(miss_qc, collapse = ", "))
}

base_qc <- qc_raw %>%
  mutate(year = suppressWarnings(as.integer(substr(as.character(tourney_date), 1, 4)))) %>%
  filter(year %in% YEARS,
         !is.na(draw_size), draw_size >= DRAW_MIN_CH)

# Challenger MAIN DRAW from this file
chall_dat <- base_qc %>%
  filter(tourney_level == "C", is_main_draw_round(round))

if (EXCLUDE_WO) chall_dat <- drop_walkovers(chall_dat)

## ---------------- Build wins-per-event for Top-tier & Challenger ----------------

wins_TOP <- build_wins_per_event(tour_dat,  "wins_in_event")   # A/M/G
wins_C   <- build_wins_per_event(chall_dat, "wins_in_event")   # Challenger

summarise_block(wins_TOP, "Top-tier", "wins_in_event")
summarise_block(wins_C,   "Challenger",   "wins_in_event")

## ---------------- Build consolidated event-run dataset ----------------

wins_TOP <- wins_TOP %>% mutate(stratum = "Top")
wins_C   <- wins_C   %>% mutate(stratum = "Challenger")

consolidated <- bind_rows(wins_TOP, wins_C) %>%
  select(
    year,
    tourney_id,
    tourney_name,
    draw_size,
    winner_name,
    wins_in_event,
    stratum,
    everything()
  ) %>%
  arrange(stratum, year, tourney_id, winner_name)

message("\nConsolidated dataset summary:")
message("  Total rows:      ", nrow(consolidated))
message("  Top-tier rows:   ", sum(consolidated$stratum == "Top"))
message("  Challenger rows: ", sum(consolidated$stratum == "Challenger"))
message("  Years:           ",
        paste(sort(unique(consolidated$year)), collapse = ", "))

## ---------------- Write output ----------------

out_name <- "tennis_player_event_runs_knockout_2016_2017_2018_2019_2021_2022_2023_fromRAW.csv"
write_csv(consolidated, out_name)

message("\nWrote consolidated file:\n  ",
        normalizePath(out_name), "\n")