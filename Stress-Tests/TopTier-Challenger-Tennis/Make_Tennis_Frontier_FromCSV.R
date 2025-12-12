# ================================================================
# Make_Tennis_Frontier_FromConsolidated_COMPACT.R
# - Horizontal layout with A. / B. tags outside panels (patchwork)
# ================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

if (interactive() &&
    requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

IN_CSV      <- "tennis_player_event_runs_knockout_2016_2017_2018_2019_2021_2022_2023_fromRAW.csv"
TOOLS_FILE  <- "FrontierTools_min.R"
OUT_FIG_DIR <- "Figures"

if (!file.exists(IN_CSV)) {
  stop("Input CSV not found: ", IN_CSV)
}
if (!file.exists(TOOLS_FILE)) {
  stop("Frontier tools file not found: ", TOOLS_FILE)
}
if (!dir.exists(OUT_FIG_DIR)) dir.create(OUT_FIG_DIR, recursive = TRUE)

# ---- load frontier tools (defines scan_one_multimodel, choose_frontier_dln) ----
source(TOOLS_FILE, local = TRUE)

# ---- read consolidated runs ----
tennis <- read_csv(IN_CSV, show_col_types = FALSE)

# basic sanity checks
req_cols <- c("wins_in_event", "stratum")
missing  <- setdiff(req_cols, names(tennis))
if (length(missing)) {
  stop("Input CSV is missing columns: ", paste(missing, collapse = ", "))
}

message("\nStratum x year counts:")
print(table(tennis$stratum, tennis$year))

message("\nTop-tier wins_in_event distribution:")
print(table(tennis$wins_in_event[tennis$stratum == "Top"]))

message("\nChallenger wins_in_event distribution:")
print(table(tennis$wins_in_event[tennis$stratum == "Challenger"]))

# ---- build scan vectors (wins >= 2) ----
wins_TOP <- tennis %>% filter(stratum == "Top")
wins_C   <- tennis %>% filter(stratum == "Challenger")

vec_TOP_ge2 <- wins_TOP$wins_in_event[wins_TOP$wins_in_event >= 2]
vec_C_ge2   <- wins_C$wins_in_event[wins_C$wins_in_event >= 2]

message("\nVector sizes (wins >= 2):")
message("  Top-tier:     N = ", length(vec_TOP_ge2))
message("  Challenger:   N = ", length(vec_C_ge2))

# ---- run scans ----
scan_top  <- scan_one_multimodel(vec_TOP_ge2)
scan_chal <- scan_one_multimodel(vec_C_ge2)

fstar_top  <- choose_frontier_dln(scan_top)
fstar_chal <- choose_frontier_dln(scan_chal)

# helper to print a compact summary at f*
summarise_frontier <- function(sc, fstar, label) {
  cat("\n---", label, "---\n")
  if (is.na(fstar)) {
    cat("  No valid frontier (choose_frontier_dln returned NA)\n")
    return(invisible(NULL))
  }
  row <- sc %>% filter(f == fstar)
  if (!nrow(row)) {
    cat("  f* =", fstar, "but no matching row in scan table.\n")
    return(invisible(NULL))
  }
  cat(sprintf("  f* = %.2f; kmin* = %s; n_kept = %s\n",
              fstar, row$kmin[1], row$n_kept[1]))
  cat(sprintf("  p_dLN = %.3g; p_Zipf = %.3g; p_Geom = %.3g\n",
              row$p_dLN[1], row$p_Zipf[1], row$p_Geom[1]))
}

summarise_frontier(scan_top,  fstar_top,  "ATP Top-tier")
summarise_frontier(scan_chal, fstar_chal, "ATP Challenger")

# ------------------------------------------------
# Plotting helper (no panel_label; tags via patchwork)
# ------------------------------------------------
plot_scan <- function(df, title_str, fstar, show_xlab = TRUE) {
  long <- df %>%
    select(f, n_kept, p_dLN, p_Zipf, p_Geom) %>%
    pivot_longer(starts_with("p_"),
                 names_to = "model", values_to = "pval") %>%
    mutate(
      model = recode(model,
                     p_dLN  = "dLN",
                     p_Zipf = "Zipf",
                     p_Geom = "Geom"),
      model = factor(model, levels = c("dLN","Zipf","Geom"))
    )
  
  # Values for f* and N_kept at the frontier (if any)
  if (!is.na(fstar)) {
    row        <- df %>% filter(f == fstar)
    f_label    <- round(fstar, 2)
    n_kept_val <- ifelse(nrow(row), row$n_kept[1], NA)
  } else {
    f_label    <- NA
    n_kept_val <- NA
  }
  
  p <- ggplot(long, aes(x = f, y = pval, color = model)) +
    geom_hline(yintercept = 0.05, linetype = "dashed",
               linewidth = 0.25) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.4) +
    scale_color_manual(values = c("dLN"  = "#d62728",
                                  "Zipf" = "#1f77b4",
                                  "Geom" = "#2ca02c"),
                       name = NULL) +
    scale_x_continuous(limits = c(0.20, 0.80),
                       breaks = seq(0.20, 0.80, by = 0.10)) +
    scale_y_continuous(limits = c(0, 1),
                       breaks = seq(0, 1, by = 0.2)) +
    labs(
      title = NULL,  # labels go inside panel
      x     = if (show_xlab) "Retained fraction f" else NULL,
      y     = "KS bootstrap p-value"
    ) +
    theme_minimal(base_size = 9) +
    theme(
      legend.position   = "top",
      axis.title.x      = element_text(size = 9),
      axis.title.y      = element_text(size = 9),
      panel.grid.minor  = element_blank(),
      plot.margin       = margin(2, 4, 2, 4, "pt")
    )
  
  if (!is.na(fstar)) {
    p <- p + geom_vline(xintercept = fstar,
                        linetype = "dotted",
                        linewidth = 0.5)
  }
  
  # Title inside: plain text at top-right
  p <- p + annotate(
    "text",
    x        = 0.80,
    y        = 0.98,
    label    = title_str,
    hjust    = 1,
    vjust    = 1,
    size     = 3,
    fontface = "bold"
  )
  
  # f* and N_kept on one line, plain text, just below the title
  if (!is.na(f_label) && !is.na(n_kept_val)) {
    label_char <- sprintf("f* = %.2f; N_kept = %d", f_label, n_kept_val)
    
    p <- p + annotate(
      "text",
      x      = 0.80,
      y      = 0.90,
      label  = label_char,
      hjust  = 1,
      vjust  = 1,
      size   = 3
    )
  }
  
  p
}

# Left panel: Top-tier
p_left <- plot_scan(
  scan_top,
  "ATP Top-tier — frontier scan",
  fstar_top,
  show_xlab = TRUE
)

# Right panel: Challenger
p_right <- plot_scan(
  scan_chal,
  "ATP Challenger — frontier scan",
  fstar_chal,
  show_xlab = TRUE
)

# Combine panels side-by-side with shared legend (top-right) and A./B. tags
combo <- (p_left | p_right) +
  plot_layout(widths = c(1, 1), guides = "collect") +
  plot_annotation(
    tag_levels = "A",
    tag_prefix = "",
    tag_suffix = "."
  ) &
  theme(
    legend.position   = "top",          # legend stays on top
    legend.justification = "right",     # right-align legend
    legend.box.just   = "right",        # keep box right-justified
    plot.margin       = margin(2, 6, 2, 6, "pt"),
    plot.tag          = element_text(face = "bold", size = 10),
    plot.tag.position = c(0.02, 0.96)   # tags closer to panels (tweak if needed)
  )

pdf_path <- file.path(OUT_FIG_DIR,
                      "Tennis_TopTier_vs_Challenger_FrontierScans.pdf")
png_path <- file.path(OUT_FIG_DIR,
                      "Tennis_TopTier_vs_Challenger_FrontierScans.png")

ggsave(pdf_path, combo, width = 7.6, height = 3.6, device = cairo_pdf)
ggsave(png_path, combo, width = 7.6, height = 3.6, dpi = 300)
message("\nSaved figure:\n  ", normalizePath(pdf_path),
        "\n  ", normalizePath(png_path), "\n")