library(readr)
library(dplyr)
library(ggplot2)
library(ineq)
library(scales)

# Optional: set working directory to script location in RStudio
if (interactive() &&
    requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

## 1. Read filled breakdown files
bio_raw <- read_csv("SCHOLARS_BIO_breakdown_withRank_filled.csv",
                    show_col_types = FALSE)
cs_raw  <- read_csv("SCHOLARS_CS_breakdown_withRank_filled.csv",
                    show_col_types = FALSE)

## 2. Helper to add total_awards and first_award_year, then apply modeling filters
prepare_faculty <- function(df) {
  # columns y1996, y1997, ..., y2016
  year_cols <- grep("^y[0-9]{4}$", names(df), value = TRUE)
  years     <- as.integer(sub("y", "", year_cols))
  Y         <- as.matrix(df[, year_cols])
  
  # total awards over the whole period
  total_awards <- rowSums(Y, na.rm = TRUE)
  
  # index of first year with > 0 awards
  first_idx <- max.col(Y > 0, ties.method = "first")
  first_idx[total_awards == 0] <- NA_integer_
  first_year <- ifelse(is.na(first_idx), NA_integer_, years[first_idx])
  
  df %>%
    mutate(
      total_awards     = total_awards,
      first_award_year = first_year
    ) %>%
    filter(
      total_awards >= 2,              # at least 2 awards
      !is.na(first_award_year),
      first_award_year <= 2010,       # first award no later than 2010
      !is.na(SchoolRank)              # must have a department rank
    )
}

bio <- prepare_faculty(bio_raw)
cs  <- prepare_faculty(cs_raw)

## 3. Build cumulative curves by SchoolRank (drop ranks with 0 awards)
make_cum_curve <- function(df, domain) {
  by_rank <- df %>%
    group_by(SchoolRank) %>%
    summarise(
      awards = sum(total_awards, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(awards > 0) %>%      # only departments that actually compete
    arrange(SchoolRank) %>%
    mutate(
      cum_awards = cumsum(awards),
      cum_share  = cum_awards / sum(awards),
      domain     = domain
    )
  
  by_rank
}

curve_bio <- make_cum_curve(bio, "BIO")
curve_cs  <- make_cum_curve(cs,  "CS")

## 4. Ginis on department-level totals
gini_bio <- Gini(curve_bio$awards)
gini_cs  <- Gini(curve_cs$awards)

label_bio <- sprintf("BIO (Gini=%.2f)", gini_bio)
label_cs  <- sprintf("CS (Gini=%.2f)", gini_cs)

curve_bio$label <- label_bio
curve_cs$label  <- label_cs

lorenz_df <- bind_rows(curve_bio, curve_cs)

## 5. Equality line over observed rank range
n_eq <- 200  # resolution of the diagonal
eq_line <- data.frame(
  SchoolRank = seq(min(lorenz_df$SchoolRank),
                   max(lorenz_df$SchoolRank),
                   length.out = n_eq),
  cum_share  = seq(0, 1, length.out = n_eq)
)

## 6. Interpolation helper for annotations
interp_at <- function(curve, r) {
  x <- curve$SchoolRank
  y <- curve$cum_share
  
  if (r <= min(x)) return(y[1])
  if (r >= max(x)) return(y[length(y)])
  
  idx <- max(which(x <= r))
  x0 <- x[idx]; x1 <- x[idx + 1]
  y0 <- y[idx]; y1 <- y[idx + 1]
  
  y0 + (y1 - y0) * (r - x0) / (x1 - x0)
}

mark_ranks <- c(3, 10, 20)

get_annots <- function(curve, domain_label, domain_code) {
  data.frame(
    SchoolRank = mark_ranks,
    cum_share  = vapply(mark_ranks, interp_at, numeric(1), curve = curve),
    label      = domain_label,
    domain     = domain_code
  ) %>%
    mutate(text = paste0(domain, "@", SchoolRank, " = ",
                         percent(cum_share, accuracy = 0.1)))
}

annot_pts <- bind_rows(
  get_annots(curve_bio, label_bio, "BIO"),
  get_annots(curve_cs,  label_cs,  "CS")
) %>%
  mutate(
    # smaller domain-specific offsets so labels sit closer to their points
    dx = ifelse(domain == "BIO",  0.7, -0.7),
    dy = ifelse(domain == "BIO",  0.04, -0.04),
    x_text = SchoolRank + dx,
    y_text = pmin(pmax(cum_share + dy, 0.05), 0.95)
  )

## 7. Colors
cols <- c("#E69F00", "#0072B2")  # orange (BIO), blue (CS)
names(cols) <- c(label_bio, label_cs)

## 8. Plot
gg <- ggplot(lorenz_df,
             aes(x = SchoolRank, y = cum_share, colour = label)) +
  # BIO & CS curves
  geom_line(linewidth = 1.1) +
  # Equality line – black, dashed
  geom_line(data = eq_line,
            aes(x = SchoolRank, y = cum_share),
            inherit.aes = FALSE,
            linetype = "dashed",
            colour  = "black",
            linewidth = 0.7) +
  # Vertical dotted lines at 3, 10, 20
  geom_vline(xintercept = mark_ranks,
             linetype = "dotted", colour = "orange", alpha = 0.6) +
  # Points at the true positions for BIO/CS at the marked ranks
  geom_point(
    data = annot_pts,
    aes(x = SchoolRank, y = cum_share, colour = label),
    inherit.aes = FALSE,
    size = 2
  ) +
  # Labels nudged away, on white boxes (no extra geom_text)
  geom_label(
    data = annot_pts,
    aes(x = x_text, y = y_text, label = text, colour = label),
    inherit.aes = FALSE,
    hjust = 0, vjust = 0.5, size = 3,
    fill = "white",
    linewidth = 0,
    show.legend = FALSE   # <- add this
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  scale_x_continuous(
    breaks = seq(0, max(lorenz_df$SchoolRank), by = 10),
    limits = c(0, max(lorenz_df$SchoolRank))
  ) +
  scale_colour_manual(values = cols, name = NULL) +
  labs(
    title = "Cumulative Share of Awards by Department Rank",
    x     = "Department rank (1 = highest)",
    y     = "Cumulative share of awards [%]"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position   = c(0.8, 0.2),
    legend.background = element_rect(fill = "white", colour = "grey80"),
    panel.grid.minor  = element_blank()
  )

gg

## 9. Save as PDF
ggsave(
  filename = "Cumulative_Awards_by_DepartmentRank_filtered.pdf",
  plot     = gg,
  width    = 7,
  height   = 4.5,
  units    = "in",
  device   = cairo_pdf  # or "pdf" if cairo_pdf is not available
)