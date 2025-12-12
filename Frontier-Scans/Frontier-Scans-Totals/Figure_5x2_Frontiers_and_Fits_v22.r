# Figure_5x2_Frontiers_and_Fits_TOTALS_v21_fixed.R
# VISUAL-TWEAKS ONLY vs v20/v21 (computations unchanged):
# - Legends only in first row panels (scan: dLN/Zipf/Geom; fit: Empirical/dLN) and nudged up.
# - Slightly more white space between rows.
# - Common x-axis labels at the bottom ("Retained fraction [%]" and "k").
# - Right-side row labels placed per-row in consistent tall gray boxes.
# - No undefined objects: final assembly uses header_row + body_with_spacers + bottom_row.

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(cowplot); library(stringr); library(grid)
})

if (interactive() &&
    requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

# ---- Font setup (legibility) ----
BASE_FAMILY <- "sans"  # fallback

if (requireNamespace("showtext", quietly = TRUE) &&
    requireNamespace("sysfonts", quietly = TRUE)) {
  
  library(showtext)
  # Pick ONE font; comment the others:
  sysfonts::font_add_google("Fira Sans", "fira")
  # sysfonts::font_add_google("Source Sans 3", "ss3")
  # sysfonts::font_add_google("Lato", "lato")
  
  showtext::showtext_auto()
  BASE_FAMILY <- "fira"   # or "ss3"/"lato" if you chose those
}

# Make it the default everywhere
theme_set(theme_minimal(base_family = BASE_FAMILY, base_size = 10))


# ----------------------------- Settings (unchanged) -----------------------------
F_GRID     <- seq(0.20, 0.80, by = 0.05)
MIN_N_KEPT <- 40
B_BOOT     <- 200
BASE_SEED  <- 202509

FIG_W <- 6.5; FIG_H <- 9

# --- Scan column styling ---
SCAN_LINE_WIDTH <- 0.45   # thinner model lines (was 0.7)
SCAN_POINT_SIZE <- 1.6    # dot size (was 1.9)
SCAN_LINE_ALPHA <- 0.95   # a touch of transparency
REF_LINE_WIDTH  <- 0.30   # width of the y=0.05 reference line

# small horizontal gap between scan & fit columns
H_GAP   <- 0.09                       # try 0.04–0.10 to taste
spacerH <- ggplot() + theme_void()    # blank spacer 

# vertical gap between rows
ROW_GAP  <- 0.09                   # try 0.04–0.10 to taste
spacerV  <- ggplot() + theme_void()

COL_EMP  <- "#4D4D4D"  # grey
COL_DLN  <- "#D62728"  # red
COL_ZIPF <- "#1F77B4"  # blue
COL_GEOM <- "#2CA02C"  # green

# --- Side-label styling (right-hand grey boxes per row) ---
SIDE_FILL   <- "#E9E9E9"  # box color
SIDE_TEXT   <- 12.5       # label font size (pt)
BOX_LEN     <- 1.0        # fraction of row height covered by the box (↑ = longer)
BOX_THICK   <- 0.45       # fraction of the per-row side area width used as thickness (↓ = slimmer)

# ----------------------------- Helpers (unchanged logic) ------------------------
sanitize_names <- function(df) {
  nm <- names(df); nm[is.na(nm) | nm==""] <- paste0("V", which(is.na(nm) | nm==""))
  names(df) <- make.names(nm, unique = TRUE); df
}
extract_numeric_column <- function(df, prefer=c("score","x","total","points","medal_points","total_score",
                                                "victories","kills","hits")) {
  nms <- tolower(names(df))
  hit <- intersect(prefer, nms)
  if (length(hit)) {
    col <- names(df)[match(hit[1], nms)]
    x <- suppressWarnings(as.integer(df[[col]]))
    return(x[!is.na(x)])
  }
  num <- which(sapply(df, is.numeric))
  if (!length(num)) stop("No numeric column found to use as score.")
  rng <- sapply(df[num], function(v) diff(range(v, na.rm = TRUE)))
  col <- names(df)[num[which.max(rng)]]
  x <- suppressWarnings(as.integer(df[[col]]))
  x[!is.na(x)]
}
sum_year_columns <- function(df, year_regex="^y[0-9]{4}$|^y(199[6-9]|20(0[0-9]|1[0-6]))$") {
  ycols <- grep(year_regex, names(df), value = TRUE)
  if (!length(ycols)) stop("No year columns matched for totals.")
  rowSums(df[ycols], na.rm = TRUE)
}
kmin_by_fraction <- function(k, f) {
  k <- sort(k)
  idx <- ceiling((1 - f) * length(k))
  idx <- max(1, min(idx, length(k)))
  as.integer(k[idx])
}

# --------------------- Discrete Lognormal (unchanged) ---------------------------
Phi <- function(z) pnorm(z)
dln_loglik <- function(par, k, kmin) {
  mu <- par[1]; sigma <- max(1e-3, abs(par[2]))
  a  <- (log(kmin - 0.5) - mu)/sigma
  Zd <- max(1e-16, 1 - Phi(a))
  bH <- (log(k + 0.5) - mu)/sigma
  bL <- (log(k - 0.5) - mu)/sigma
  pm <- pmax(Phi(bH) - Phi(bL), 0) / Zd
  sum(log(pmax(pm, 1e-16)))
}
fit_dln <- function(k, kmin) {
  kk <- as.integer(k[k >= kmin]); if (length(kk) < 2) stop("dLN: n<2")
  logs <- log(pmax(as.numeric(kk), 2))
  mu0 <- mean(logs); if (!is.finite(mu0)) mu0 <- median(logs)
  s0  <- sd(logs);   if (!is.finite(s0) || s0 < 1e-3) s0 <- max(0.5, stats::mad(logs))
  opt <- optim(c(mu0, s0),
               fn = function(p) -dln_loglik(p, kk, kmin),
               method = "L-BFGS-B",
               lower = c(mu0 - 5, 0.05), upper = c(mu0 + 5, 5))
  list(mu = opt$par[1], sigma = max(0.05, abs(opt$par[2])))
}
dln_cdf <- function(kq, mu, sigma, kmin) {
  kq <- as.integer(kq); J <- seq.int(kmin, max(kq))
  a  <- (log(kmin - 0.5) - mu)/sigma
  Zd <- max(1e-16, 1 - Phi(a))
  bH <- (log(J + 0.5) - mu)/sigma
  bL <- (log(J - 0.5) - mu)/sigma
  pm <- pmax(Phi(bH) - Phi(bL), 0) / Zd
  Fj <- cumsum(pm)
  Fj[pmax(1, findInterval(kq, J))]
}
dln_ccdf_smooth <- function(x, mu, sigma, kmin) {
  den <- max(1e-12, 1 - Phi((log(kmin - 0.5) - mu)/sigma))
  pmax((1 - Phi((log(x) - mu)/sigma))/den, 1e-12)
}
rdln <- function(n, mu, sigma, kmin) {
  out <- integer(0); thr <- kmin - 0.5
  while (length(out) < n) {
    m <- ceiling((n - length(out)) * 1.6)
    x <- rlnorm(m, meanlog = mu, sdlog = sigma); x <- x[x >= thr]
    if (length(x)) out <- c(out, as.integer(round(x)))
  }
  out[seq_len(n)]
}

# --------------------- Zipf & Geom (unchanged) -----------------------------------
fit_zipf_finite <- function(k, kmin) {
  kk <- as.integer(k[k >= kmin]); ks <- seq.int(kmin, max(kk))
  logs <- sum(log(kk)); n <- length(kk)
  nll <- function(a) { if (a <= 1.0001) return(1e12); w <- ks^(-a); a*logs + n*log(sum(w)) }
  opt <- optimize(nll, interval = c(1.001, 10), maximum = FALSE)
  a <- opt$minimum; w <- ks^(-a); pm <- w/sum(w); cdf <- cumsum(pm)
  list(alpha = a, ks = ks, pmf = pm, cdf = cdf)
}
zipf_cdf <- function(fit, kq) {
  idx <- pmax(1, findInterval(as.integer(kq), fit$ks))
  fit$cdf[pmin(idx, length(fit$cdf))]
}
rzipf_finite <- function(n, fit) sample(fit$ks, size = n, replace = TRUE, prob = fit$pmf)

fit_geom <- function(k, kmin) {
  kk <- as.integer(k[k >= kmin]); y <- kk - kmin; ybar <- mean(y)
  list(theta = if (!is.na(ybar) && ybar > 0) ybar/(1 + ybar) else 1e-4)
}
geom_cdf <- function(kq, theta, kmin) {
  y <- (as.integer(kq) - kmin) + 1L
  pmin(pmax(1 - (theta^y), 0), 1)
}
rgeom_shift <- function(n, theta, kmin) {
  p <- max(1e-12, 1 - theta)
  kmin + rgeom(n, prob = p)
}

# --------------------- KS & bootstrap (unchanged) --------------------------------
ks_stat_discrete <- function(k, F_model) {
  k <- sort(as.integer(k)); u <- sort(unique(k))
  counts <- tabulate(match(k, u), nbins = length(u))
  ecdf_u <- cumsum(counts)/sum(counts)
  modF   <- F_model(u)
  max(abs(ecdf_u - modF))
}
ks_bootstrap_p <- function(kept, cdf_fn, model, params, B = 200L, seed_offset = 0L) {
  set.seed(BASE_SEED + switch(model, dLN=1L, Zipf=2L, Geom=3L) + as.integer(seed_offset))
  n <- length(kept)
  obs <- ks_stat_discrete(kept, cdf_fn)
  geq <- 0L
  for (b in seq_len(B)) {
    sim <- switch(model,
                  dLN  = rdln(n, params$mu, params$sigma, params$kmin),
                  Zipf = rzipf_finite(n, params$zf_fit),
                  Geom = rgeom_shift(n, params$theta, params$kmin)
    )
    stat <- ks_stat_discrete(sim, cdf_fn)
    if (stat >= obs) geq <- geq + 1L
  }
  (geq + 1)/(B + 1)
}

# --------------------- Scans & frontier choice (unchanged) -----------------------
scan_models <- function(x, f_grid = F_GRID, min_n = MIN_N_KEPT) {
  out <- vector("list", length(f_grid))
  for (i in seq_along(f_grid)) {
    f    <- f_grid[i]
    kmin <- kmin_by_fraction(x, f)
    kept <- x[x >= kmin]; n_kept <- length(kept)
    if (n_kept < min_n || length(unique(kept)) < 2L) {
      out[[i]] <- data.frame(f=f, kmin=kmin, n_kept=n_kept,
                             p_dLN=NA_real_, p_Zipf=NA_real_, p_Geom=NA_real_)
      next
    }
    dln <- fit_dln(kept, kmin)
    zf  <- fit_zipf_finite(kept, kmin)
    gm  <- fit_geom(kept, kmin)
    p_dln <- ks_bootstrap_p(kept, function(t) dln_cdf(sort(unique(t)), dln$mu, dln$sigma, kmin),
                            "dLN",  list(mu=dln$mu, sigma=dln$sigma, kmin=kmin), B_BOOT, kmin)
    p_zip <- ks_bootstrap_p(kept, function(t) zipf_cdf(zf, sort(unique(t))),
                            "Zipf", list(zf_fit=zf), B_BOOT, kmin)
    p_geo <- ks_bootstrap_p(kept, function(t) geom_cdf(sort(unique(t)), gm$theta, kmin),
                            "Geom", list(theta=gm$theta, kmin=kmin), B_BOOT, kmin)
    out[[i]] <- data.frame(f=f, kmin=kmin, n_kept=n_kept,
                           p_dLN=p_dln, p_Zipf=p_zip, p_Geom=p_geo)
  }
  bind_rows(out)
}
choose_frontier_dln <- function(scan_df) {
  ok <- scan_df %>% filter(!is.na(p_dLN), p_dLN >= 0.05, n_kept >= MIN_N_KEPT)
  if (nrow(ok) > 0) return(max(ok$f))
  avail <- scan_df %>% filter(n_kept >= MIN_N_KEPT)
  if (nrow(avail) == 0) return(NA_real_)
  max(avail$f, na.rm = TRUE)
}

# --------------------- Plot builders (legend tweak only) -------------------------
scan_plot <- function(sc, show_legend=FALSE) {
  df <- sc %>%
    transmute(f=100*f, dLN=p_dLN, Zipf=p_Zipf, Geom=p_Geom) %>%
    tidyr::pivot_longer(cols=c("dLN","Zipf","Geom"), names_to="model", values_to="p") %>%
    mutate(model=factor(model, levels=c("dLN","Zipf","Geom")))
  
  p <- ggplot(df, aes(f, p, color=model)) +
    geom_hline(yintercept=0.05, linetype="dashed", color="gray60",
               linewidth=REF_LINE_WIDTH) +
    geom_line(na.rm=TRUE, linewidth=SCAN_LINE_WIDTH, alpha=SCAN_LINE_ALPHA) +
    geom_point(na.rm=TRUE, size=SCAN_POINT_SIZE) +
    scale_color_manual(values=c(dLN=COL_DLN, Zipf=COL_ZIPF, Geom=COL_GEOM), name=NULL) +
    coord_cartesian(ylim=c(0,1)) +
    scale_x_continuous(breaks=seq(20,80,20)) +
    labs(x=NULL, y="KS bootstrap p") +
    theme_minimal(base_size=10) +
    theme(panel.grid.minor=element_blank(),
          plot.margin=margin(2,3,2,2))
  
  if (show_legend) {
    p + theme(legend.position=c(0.98,0.98),
              legend.justification=c(1,1),
              legend.background=element_rect(fill=scales::alpha("white",0.75), color=NA),
              legend.direction="vertical")
  } else p + theme(legend.position="none")
}

fit_plot_dln <- function(x, f_star, show_legend=FALSE) {
  if (is.na(f_star)) {
    return(ggplot() + theme_void() +
             annotate("text", x=0.5, y=0.5, label="Insufficient data",
                      size=3.5, colour="gray30"))
  }
  kmin <- kmin_by_fraction(x, f_star)
  kept <- x[x >= kmin]
  if (length(kept) < 2L) {
    return(ggplot() + theme_void() +
             annotate("text", x=0.5, y=0.5, label="Insufficient data",
                      size=3.5, colour="gray30"))
  }
  
  # empirical CCDF
  emp_vals   <- sort(unique(kept))
  emp_counts <- tabulate(match(kept, emp_vals), nbins=length(emp_vals))
  emp_ccdf   <- 1 - cumsum(emp_counts)/sum(emp_counts) + emp_counts/sum(emp_counts)
  df_emp     <- data.frame(x=rep(emp_vals, each=2)[-1],
                           y=rep(emp_ccdf, each=2)[-1],
                           model="Empirical")
  
  # dLN fit
  fit   <- fit_dln(kept, kmin)
  xgrid <- exp(seq(log(kmin-0.5), log(max(emp_vals)+0.5), length.out=600))
  df_d  <- data.frame(x=xgrid, y=dln_ccdf_smooth(xgrid, fit$mu, fit$sigma, kmin), model="dLN")
  
  # zoom
  N_keep <- length(kept); y_floor <- max(1e-4, (1/N_keep)/3)
  x_max  <- {
    idx <- which(df_emp$y >= y_floor)
    if (length(idx)) min(max(df_emp$x[idx]) * 1.15, max(emp_vals)) else max(emp_vals)
  }
  
  p <- ggplot() +
    # draw dLN FIRST (goes underneath), with slight transparency
    geom_line(
      data = df_d, aes(x, y, color = model),
      linewidth = 1.0, alpha = 0.65, lineend = "round"
    ) +
    # draw Empirical SECOND (on top), a hair thicker for visibility
    geom_step(
      data = df_emp, aes(x, y, color = model),
      direction = "hv", linewidth = 0.9
    ) +
    scale_color_manual(values = c(Empirical = COL_EMP, dLN = COL_DLN), name = NULL) +
    scale_x_log10(limits = c(kmin - 0.5, x_max), expand = c(0, 0)) +
    scale_y_log10(limits = c(y_floor, 1.0), expand = c(0, 0)) +
    labs(x = NULL, y = expression(P(K >= k))) +
    theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          plot.margin = margin(2, 3, 2, 2))
  
  if (show_legend) {
    p + theme(
      legend.position      = c(0.965, 0.999),     # a hair higher
      legend.justification = c(1, 1),
      legend.background    = element_blank(),     # <- no white box
      legend.box.background= element_blank(),
      legend.key           = element_blank(),
      legend.key.height    = unit(0.28, "lines"),
      legend.key.width     = unit(0.95, "lines"),
      legend.text          = element_text(size = 9)
    )
  } else {
    p + theme(legend.position = "none")
  }
}

# --------------------- Side labels: tall gray boxes per row ----------------------
make_side_label <- function(txt) {
  ggplot() +
    geom_rect(
      aes(xmin = 0.5 - BOX_THICK/2, xmax = 0.5 + BOX_THICK/2,
          ymin = 0.5 - BOX_LEN/2,   ymax = 0.5 + BOX_LEN/2),
      fill = SIDE_FILL, colour = NA
    ) +
    annotate(
      "text", x = 0.5, y = 0.5, label = txt,
      angle = 270, fontface = "bold",
      size = SIDE_TEXT / ggplot2::.pt, colour = "black"
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_void()
}

# --------------------- Domain loaders (totals; UPDATED) ------------------------

# LW: only pilots whose *first hit* is in 1940, 1941, or 1942
# This should give N = 1148, matching your Descriptive Statistics.
load_LW_totals <- function() {
  df <- sanitize_names(read.csv("GER_pilots.csv", check.names = TRUE, strip.white = TRUE))
  
  ycols <- grep("^y[0-9]{4}$", names(df), value = TRUE)
  if (!length(ycols)) stop("LW: no year columns found.")
  
  years <- as.integer(sub("^y", "", ycols))
  mat   <- as.matrix(df[, ycols, drop = FALSE])
  
  # Total confirmed victories across all years
  totals <- rowSums(mat, na.rm = TRUE)
  
  # First hit year for each pilot
  first_year <- apply(mat, 1, function(x) {
    idx <- which(x > 0)
    if (length(idx) == 0L) NA_integer_ else years[min(idx)]
  })
  
  keep <- !is.na(first_year) & first_year %in% c(1940L, 1941L, 1942L)
  # (All of these have totals >= 2 anyway, so no extra threshold needed.)
  as.integer(totals[keep])
}

# BIO: investigators whose first grant is by 2010, and total grants >= 2
# Should yield N = 524, matching Descriptive Statistics.
load_BIO_totals <- function() {
  df <- sanitize_names(read.csv("SCHOLARS_BIO_breakdown.csv", check.names = TRUE))
  
  ycols <- grep("^y(199[6-9]|20(0[0-9]|1[0-6]))$", names(df), value = TRUE)
  if (!length(ycols)) stop("BIO: no year columns found.")
  
  years <- as.integer(sub("^y", "", ycols))
  mat   <- as.matrix(df[, ycols, drop = FALSE])
  
  totals <- rowSums(mat, na.rm = TRUE)
  
  first_year <- apply(mat, 1, function(x) {
    idx <- which(x > 0)
    if (length(idx) == 0L) NA_integer_ else years[min(idx)]
  })
  
  keep <- !is.na(first_year) & first_year <= 2010L & totals >= 2
  as.integer(totals[keep])
}

# CS: investigators whose first grant is by 2010, and total grants >= 2
# Should yield N = 1323, matching Descriptive Statistics.
load_CS_totals <- function() {
  df <- sanitize_names(read.csv("SCHOLARS_CS_breakdown.csv", check.names = TRUE))
  
  ycols <- grep("^y(199[6-9]|20(0[0-9]|1[0-6]))$", names(df), value = TRUE)
  if (!length(ycols)) stop("CS: no year columns found.")
  
  years <- as.integer(sub("^y", "", ycols))
  mat   <- as.matrix(df[, ycols, drop = FALSE])
  
  totals <- rowSums(mat, na.rm = TRUE)
  
  first_year <- apply(mat, 1, function(x) {
    idx <- which(x > 0)
    if (length(idx) == 0L) NA_integer_ else years[min(idx)]
  })
  
  keep <- !is.na(first_year) & first_year <= 2010L & totals >= 2
  as.integer(totals[keep])
}

# US swimmers: athletes with total score (medals) >= 2
# This already gives N = 446 with the current breakdown file.
load_US_totals <- function() {
  df <- sanitize_names(read.csv("ATHLETES_USA_SWIMMERS_breakdown.csv", check.names = TRUE))
  
  # Sum across all Games columns: X1896.Summer, X1900.Summer, ..., X2016.Summer
  totals <- sum_year_columns(df, "^X[0-9]{4}\\.")
  
  # Keep only swimmers with 2+ medals
  as.integer(totals[totals >= 2])
}

# FRA fencers: athletes with total score (medals) >= 2
# With the raw breakdown file and only "score >= 2", this yields N = 137.
# Your Descriptive Statistics N = 56 must involve an additional temporal
# window or other filter; that isn’t encoded here, so this loader sticks
# to the literal “score >= 2” rule.
load_FRA_totals <- function() {
  df <- sanitize_names(read.csv("ATHLETES_FRA_FENCING_breakdown.csv", check.names = TRUE))
  
  # Sum across all Games columns: X1896.Summer, X1900.Summer, ..., X2016.Summer
  totals <- sum_year_columns(df, "^X[0-9]{4}\\.")
  
  # Keep only fencers with 2+ medals
  as.integer(totals[totals >= 2])
}

# --------------------- Build computations (unchanged) ---------------------------
domains <- list(
  LW      = load_LW_totals(),
  BIO     = load_BIO_totals(),
  CS      = load_CS_totals(),
  USOT_S  = load_US_totals(),
  FRAOT_F = load_FRA_totals()
)

set.seed(BASE_SEED)
scans  <- lapply(domains, scan_models)
fstars <- lapply(scans,  choose_frontier_dln)

# --------------------- Assemble rows with a bit more spacing ---------------------
rows <- list(); nm_vec <- names(domains)
for (i in seq_along(nm_vec)) {
  nm <- nm_vec[i]
  left  <- scan_plot(scans[[nm]], show_legend = (i == 1))          # legend only in first scan
  right <- fit_plot_dln(domains[[nm]], fstars[[nm]], show_legend = (i == 1))  # legend only in first fit
  side  <- make_side_label(nm)
  rows[[i]] <- plot_grid(
    left, spacerH, right, side,
    ncol = 4,
    rel_widths = c(1, H_GAP, 1, 0.18),
    align = "h"
  )
}

# Column headers (tighter vertical space)
header_left  <- ggplot() + theme_void() +
  annotate("text", x=0.5, y=0.5, label="Frontier scan", fontface="bold", size=3.6)
header_right <- ggplot() + theme_void() +
  annotate("text", x=0.5, y=0.5, label="dLN fit at frontier", fontface="bold", size=3.6)
# headers
header_row <- plot_grid(
  header_left, spacerH, header_right, ggplot() + theme_void(),
  ncol = 4,
  rel_widths = c(1, H_GAP, 1, 0.18)
)

# Interleave small vertical spacers between rows to add a bit more white space
spacer <- ggplot() + theme_void()
body_with_spacers <- plot_grid(
  rows[[1]], spacerV,
  rows[[2]], spacerV,
  rows[[3]], spacerV,
  rows[[4]], spacerV,
  rows[[5]],
  ncol = 1,
  rel_heights = c(1, ROW_GAP, 1, ROW_GAP, 1, ROW_GAP, 1, ROW_GAP, 1)
)

# Bottom common x-axis labels for the two main columns
xlab_scan <- ggplot() + theme_void() +
  annotate("text", x=0.5, y=0.5, label="Retained fraction [%]", size=3.4)
xlab_fit  <- ggplot() + theme_void() +
  annotate("text", x=0.5, y=0.5, label="k", size=3.4)
# bottom x-axis labels
bottom_row <- plot_grid(
  xlab_scan, spacerH, xlab_fit, ggplot() + theme_void(),
  ncol = 4,
  rel_widths = c(1, H_GAP, 1, 0.18)
)

# -------- Final figure (no undefined objects; single assembled panel) -----------
final <- plot_grid(
  header_row,
  body_with_spacers,
  bottom_row,
  ncol = 1,
  rel_heights = c(0.065, 0.87, 0.065)
)

ggsave("Figure_5x2_Frontiers_and_Fits_TOTALS_v22.pdf",
       final, width = FIG_W, height = FIG_H, units = "in",
       device = cairo_pdf)   # <- better kerning/embedding than useDingbats=FALSE

# ================== Frontier-scan summary tables ==================

# Helper: vectorized winner picker (handles NAs and ties)
winner_model_vec <- function(pd, pz, pg) {
  vapply(seq_along(pd), function(i) {
    v <- c(dLN = pd[i], Zipf = pz[i], Geom = pg[i])
    if (all(is.na(v))) return(NA_character_)
    mx <- max(v, na.rm = TRUE)
    w  <- names(v)[which(!is.na(v) & v == mx)]
    if (length(w) == 1) w else paste(w, collapse = "=")
  }, character(1))
}

# ---------- 1) Detailed table: every scan point for all domains ----------
detailed_scan_tbl <- do.call(
  rbind,
  lapply(names(scans), function(nm) {
    sc <- scans[[nm]]
    if (!nrow(sc)) return(NULL)
    best_p <- do.call(pmax, c(sc[c("p_dLN","p_Zipf","p_Geom")], list(na.rm = TRUE)))
    best_p[is.infinite(best_p)] <- NA_real_
    transform(
      sc,
      domain      = nm,
      retained_pct = round(100 * f, 1),
      winner       = winner_model_vec(p_dLN, p_Zipf, p_Geom),
      adequate_dLN = !is.na(p_dLN) & p_dLN >= 0.05,
      adequate_any = !is.na(best_p) & best_p >= 0.05
    )[ , c("domain","retained_pct","kmin","n_kept","p_dLN","p_Zipf","p_Geom","winner","adequate_dLN","adequate_any")]
  })
)

outfile_detailed <- "FrontierScan_ALL_POINTS.csv"
write.csv(detailed_scan_tbl, outfile_detailed, row.names = FALSE)
message("Wrote: ", outfile_detailed)

# ---------- 2) Summary table: one row per domain at the chosen frontier f* ----------
summary_scan_tbl <- do.call(
  rbind,
  lapply(names(scans), function(nm) {
    sc <- scans[[nm]]
    f_star <- as.numeric(fstars[[nm]])
    N_total <- length(domains[[nm]])
    if (is.na(f_star) || !nrow(sc)) {
      return(data.frame(
        domain = nm, N_total = N_total,
        f_star = NA_real_, retained_pct = NA_real_,
        kmin_star = NA_integer_, n_kept_star = NA_integer_,
        p_dLN = NA_real_, p_Zipf = NA_real_, p_Geom = NA_real_,
        winner_at_f_star = NA_character_
      ))
    }
    # Find the scan row closest to f*
    i_star <- which.min(abs(sc$f - f_star))
    scs <- sc[i_star, , drop = FALSE]
    
    # Recompute kmin/n_kept from the raw data to be explicit
    kmin_star   <- kmin_by_fraction(domains[[nm]], scs$f)
    kept_star   <- domains[[nm]][domains[[nm]] >= kmin_star]
    n_kept_star <- length(kept_star)
    
    data.frame(
      domain          = nm,
      N_total         = N_total,
      f_star          = scs$f,
      retained_pct    = round(100 * scs$f, 1),
      kmin_star       = kmin_star,
      n_kept_star     = n_kept_star,
      p_dLN           = scs$p_dLN,
      p_Zipf          = scs$p_Zipf,
      p_Geom          = scs$p_Geom,
      winner_at_f_star = winner_model_vec(scs$p_dLN, scs$p_Zipf, scs$p_Geom)
    )
  })
)

outfile_summary <- "FrontierScan_SUMMARY.csv"
write.csv(summary_scan_tbl, outfile_summary, row.names = FALSE)
message("Wrote: ", outfile_summary)

# -------------------------------------------------------------------
# Also write LaTeX table rows for the frontier summary
# -------------------------------------------------------------------

# Nice labels for the domains as they should appear in the table
domain_labels <- c(
  LW      = "LW",
  BIO     = "BIO",
  CS      = "CS",
  USOT_S  = "US swimmers",
  FRAOT_F = "FRA fencers"
)

# Desired row order in the LaTeX table
row_order <- c("LW", "BIO", "CS", "USOT_S", "FRAOT_F")

summary_for_tex <- summary_scan_tbl |>
  dplyr::mutate(
    domain_label = domain_labels[domain]
  ) |>
  dplyr::filter(domain %in% row_order) |>
  dplyr::arrange(factor(domain, levels = row_order))

# Create LaTeX lines.
# Adjust the number of decimals in sprintf() to taste.
latex_rows <- with(summary_for_tex, sprintf(
  "%s & %d & %.2f & %.1f & %d & %d & %.3f & %.3f & %.3f & %s \\\\",
  domain_label,      # pretty name in first column
  N_total,           # N_total
  f_star,            # f^star
  retained_pct,      # retained fraction in percent
  kmin_star,         # k_min^star
  n_kept_star,       # N_kept(f^star)
  p_dLN,             # p_dLN
  p_Zipf,            # p_Zipf
  p_Geom,            # p_Geom
  winner_at_f_star   # winner model at f^star
))

outfile_tex <- "FrontierScan_TABLE_ROWS.tex"
writeLines(latex_rows, outfile_tex)
message("Wrote LaTeX rows to: ", outfile_tex)
# ===================================================================