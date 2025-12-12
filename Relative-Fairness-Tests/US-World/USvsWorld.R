library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(scales)
library(grid)   # for unit()

if (interactive() &&
    requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

## 1. Read Kaggle data ---------------------------------------------

ath <- read_csv("athlete_events.csv", show_col_types = FALSE)

## 2. Filter to Summer swimming, individual events, with medals ----

swim_indiv <- ath %>%
  filter(
    Season == "Summer",
    Sport  == "Swimming",
    !is.na(Medal),
    Year != 1906,               # drop Intercalated Games
    !str_detect(Event, "Relay"),# exclude relay (team) events
    !str_detect(Event, "Team")  # exclude any other team events
  )

## 3. Mark USA vs Rest of World ------------------------------------

swim_indiv <- swim_indiv %>%
  mutate(
    is_usa = if_else(NOC == "USA" | Team == "United States", 1L, 0L)
  )

## 4. Medal counts and shares by year & medal ----------------------

shares <- swim_indiv %>%
  group_by(Year, Medal) %>%
  summarise(
    total_medals = n(),
    usa_medals   = sum(is_usa),
    .groups      = "drop"
  ) %>%
  mutate(
    usa_share   = 100 * usa_medals / total_medals,
    world_share = 100 - usa_share
  ) %>%
  select(Year, Medal, usa_share, world_share) %>%
  pivot_longer(
    cols      = c(usa_share, world_share),
    names_to  = "region",
    values_to = "share"
  ) %>%
  mutate(
    region = recode(region,
                    usa_share   = "United States",
                    world_share = "Rest of World"),
    # facets top→bottom: Gold, Silver, Bronze
    Medal  = factor(Medal, levels = c("Gold", "Silver", "Bronze"))
  ) %>%
  # Treat 1980 as NA for both US and Rest of World (US boycott)
  mutate(
    share = if_else(Year == 1980, NA_real_, share)
  )

## 5. Create segment IDs to break lines across long gaps (WWI, WWII, 1980) ----

shares <- shares %>%
  arrange(Medal, region, Year) %>%
  group_by(Medal, region) %>%
  mutate(
    gap_years = c(0, diff(Year)),
    seg_id    = cumsum(gap_years > 4)   # new segment if gap > 4 years
  ) %>%
  ungroup()

## 6. Plot ----------------------------------------------------------

cols <- c("Rest of World" = "#0072B2",   # blue
          "United States" = "#D55E00")   # red

p <- ggplot(shares,
            aes(x = Year,
                y = share,
                colour = region,
                group = interaction(region, seg_id))) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_grid(Medal ~ ., scales = "fixed") +
  scale_colour_manual(values = cols, name = NULL) +
  scale_y_continuous(labels = percent_format(scale = 1),
                     limits = c(0, 100)) +
  scale_x_continuous(breaks = seq(1900, max(shares$Year, na.rm = TRUE), by = 20)) +
  labs(
    title = "Swimming Medal Shares by Olympiad (Individual Events)",
    x     = NULL,
    y     = "Share of medals [%]"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position   = "bottom",
    legend.background = element_rect(fill = "white", colour = "grey80"),
    panel.grid.minor  = element_blank(),
    panel.spacing     = unit(1.1, "cm")   # more vertical space between panels
  )

p

## 7. Save to PDF (optional) ---------------------------------------

ggsave(
  filename = "Swimming_Medal_Shares_by_Olympiad.pdf",
  plot     = p,
  width    = 7,
  height   = 5,
  units    = "in"
)