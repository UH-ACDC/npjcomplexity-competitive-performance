# Make_Period_ScanFit_Figures_ALL_unified_layout.R
# LW (2x3), BIO+CS (4x3 with scan-scan-fit-fit), Athletes 1980–2016 (1x4 scan-scan-fit-fit)

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(cowplot)
})

if (interactive() &&
    requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

# ---- Font setup (same approach as TOTALS) ----
BASE_FAMILY <- "sans"   # fallback

# --- Header band sizing for ATHLETES ---
HEADER_BAND_H_ATH <- 0.15   # was ~0.06; try 0.08–0.10 for taller grey boxes
HEADER_TXT_SZ      <- 3.9    # header text size inside grey boxes

# Side label column width for ATHLETES (tweak as needed)
SIDE_COL_W_ATH <- 0.06   # try 0.05–0.09

have_showtext <- requireNamespace("showtext", quietly = TRUE)
have_sysfonts <- requireNamespace("sysfonts", quietly = TRUE)

if (have_showtext && have_sysfonts) {
  library(showtext); library(sysfonts)
  # Pick one you like; Inter and Source Sans 3 are both airy/readable:
  # sysfonts::font_add_google("Source Sans 3", "Source Sans 3")
  sysfonts::font_add_google("Inter", "Inter")
  showtext_auto()                 # use the added font in all devices
  BASE_FAMILY <- "Inter"
}


# Make the font apply everywhere by default
theme_set(theme_minimal(base_size = 11, base_family = BASE_FAMILY))
theme_update(text = element_text(family = BASE_FAMILY, lineheight = 1.05))

# Make annotate("text", ...) and geom_label() also use the base font automatically
update_geom_defaults("text",  list(family = BASE_FAMILY))
update_geom_defaults("label", list(family = BASE_FAMILY))

# ---------------- Global settings ----------------
F_GRID     <- seq(0.20, 0.80, by = 0.05)    # retained-fraction scan grid
MIN_N_KEPT <- 40                             # minimum N after truncation
B_BOOT     <- 200                            # KS bootstrap reps (bump for camera-ready)
BASE_SEED  <- 202508

# One legend per figure
SHOW_LEGEND_LW      <- TRUE
SHOW_LEGEND_BIO_CS  <- TRUE
SHOW_LEGEND_ATH     <- TRUE

# Colors
COL_EMP  <- "#4D4D4D"  # empirical (grey)
COL_DLN  <- "#D62728"  # dLN (red)
COL_ZIPF <- "#1F77B4"  # Zipf (blue)
COL_GEOM <- "#2CA02C"  # Geom (green)

# ---------------- Utilities ----------------
sanitize_names <- function(df){
  nm <- names(df); nm[is.na(nm) | nm==""] <- paste0("V", which(is.na(nm) | nm==""))
  names(df) <- make.names(nm, unique = TRUE); df
}
first_nonzero_year <- function(row_vec, years_vec){
  ix <- which(row_vec > 0); if(!length(ix)) return(NA_integer_); years_vec[ix[1]]
}
kmin_by_fraction <- function(k, f){
  k <- sort(as.integer(k))
  idx <- ceiling((1 - f) * length(k))
  idx <- max(1, min(idx, length(k)))
  as.integer(k[idx])
}

# ---- models: discrete lognormal, finite Zipf, shifted geometric
Phi <- function(z) pnorm(z)

dln_loglik <- function(par, k, kmin){
  mu <- par[1]; sigma <- max(1e-3, abs(par[2]))
  a  <- (log(kmin-0.5)-mu)/sigma; Zd <- max(1e-16, 1 - Phi(a))
  bH <- (log(k+0.5)-mu)/sigma; bL <- (log(k-0.5)-mu)/sigma
  pm <- pmax(Phi(bH)-Phi(bL), 0) / Zd
  sum(log(pmax(pm, 1e-16)))
}
fit_dln <- function(k, kmin){
  kk <- as.integer(k[k>=kmin]); if(length(kk)<2) stop("dLN: n<2 at kmin.")
  logs <- log(as.numeric(kk))
  mu0 <- mean(logs); s0 <- sd(logs); if(!is.finite(s0) || s0<1e-3) s0 <- 0.5
  opt <- optim(c(mu0,s0), fn=function(p) -dln_loglik(p, kk, kmin),
               method="L-BFGS-B", lower=c(mu0-5,0.05), upper=c(mu0+5,5))
  list(mu=opt$par[1], sigma=max(0.05, abs(opt$par[2])))
}
dln_cdf <- function(kq, mu, sigma, kmin){
  kq <- as.integer(kq); J <- seq.int(kmin, max(kq))
  a  <- (log(kmin-0.5)-mu)/sigma; Zd <- max(1e-16, 1 - Phi(a))
  bH <- (log(J+0.5)-mu)/sigma; bL <- (log(J-0.5)-mu)/sigma
  pm <- pmax(Phi(bH)-Phi(bL), 0) / Zd; Fj <- cumsum(pm)
  Fj[pmax(1, findInterval(kq, J))]
}
dln_ccdf_smooth <- function(x, mu, sigma, kmin){
  den <- max(1e-12, 1 - Phi((log(kmin-0.5)-mu)/sigma))
  pmax((1 - Phi((log(x)-mu)/sigma))/den, 1e-12)
}
rdln <- function(n, mu, sigma, kmin){
  out <- integer(0); thr <- kmin-0.5
  while(length(out) < n){
    m <- ceiling((n-length(out))*1.6)
    x <- rlnorm(m, meanlog=mu, sdlog=sigma); x <- x[x>=thr]
    if(length(x)) out <- c(out, as.integer(round(x)))
  }
  out[seq_len(n)]
}

fit_zipf_finite <- function(k, kmin){
  kk <- as.integer(k[k>=kmin]); ks <- seq.int(kmin, max(kk))
  logs <- sum(log(kk)); n <- length(kk)
  nll <- function(a){ if(a<=1.0001) return(1e12); w <- ks^(-a); a*logs + n*log(sum(w)) }
  opt <- optimize(nll, interval=c(1.001,10), maximum=FALSE); a <- opt$minimum
  w <- ks^(-a); pm <- w/sum(w); cdf <- cumsum(pm)
  list(alpha=a, ks=ks, pmf=pm, cdf=cdf)
}
zipf_cdf <- function(fit, kq){
  idx <- pmax(1, findInterval(as.integer(kq), fit$ks))
  fit$cdf[pmin(idx, length(fit$cdf))]
}
rzipf_finite <- function(n, fit){ sample(fit$ks, size=n, replace=TRUE, prob=fit$pmf) }

fit_geom <- function(k, kmin){
  kk <- as.integer(k[k>=kmin]); y <- kk-kmin; ybar <- mean(y)
  list(theta = if(!is.na(ybar) && ybar>0) ybar/(1+ybar) else 1e-4)
}
geom_cdf <- function(kq, theta, kmin){
  y <- (as.integer(kq)-kmin)+1L
  pmin(pmax(1 - (theta^y), 0), 1)
}
rgeom_shift <- function(n, theta, kmin){
  p <- max(1e-12, 1-theta)
  kmin + rgeom(n, prob = p)
}

# ---- KS + bootstrap
ks_stat_discrete <- function(k, F_model){
  k <- sort(as.integer(k)); u <- sort(unique(k))
  counts <- tabulate(match(k, u), nbins=length(u))
  ecdf_u <- cumsum(counts)/sum(counts); modF <- F_model(u)
  max(abs(ecdf_u - modF))
}
ks_bootstrap_p <- function(kept, cdf_fn, model, params, B=200L, seed_offset=0L){
  set.seed(BASE_SEED + switch(model, dLN=1L, Zipf=2L, Geom=3L) + as.integer(seed_offset))
  n <- length(kept); obs <- ks_stat_discrete(kept, cdf_fn); geq <- 0L
  for(b in seq_len(B)){
    sim <- switch(model,
                  dLN  = rdln(n, params$mu, params$sigma, params$kmin),
                  Zipf = rzipf_finite(n, params$zf_fit),
                  Geom = rgeom_shift(n, params$theta, params$kmin))
    stat <- ks_stat_discrete(sim, cdf_fn); if(stat>=obs) geq <- geq+1L
  }
  (geq + 1)/(B + 1)
}

# ---- scan (multi-model) + frontier choice (by dLN)
scan_one_multimodel <- function(x, f_grid=F_GRID, min_n=MIN_N_KEPT, B=B_BOOT){
  out <- vector("list", length(f_grid))
  for(i in seq_along(f_grid)){
    f <- f_grid[i]; kmin <- kmin_by_fraction(x, f); kept <- x[x>=kmin]; n_kept <- length(kept)
    if(n_kept < min_n || length(unique(kept)) < 2L){
      out[[i]] <- data.frame(f=f, kmin=kmin, n_kept=n_kept,
                             p_dLN=NA_real_, p_Zipf=NA_real_, p_Geom=NA_real_)
      next
    }
    dln <- fit_dln(kept, kmin); zf <- fit_zipf_finite(kept, kmin); gm <- fit_geom(kept, kmin)
    p_dln  <- ks_bootstrap_p(kept, function(t) dln_cdf(sort(unique(t)), dln$mu, dln$sigma, kmin),
                             "dLN",  list(mu=dln$mu, sigma=dln$sigma, kmin=kmin), B, kmin)
    p_zipf <- ks_bootstrap_p(kept, function(t) zipf_cdf(zf, sort(unique(t))),
                             "Zipf", list(zf_fit=zf), B, kmin)
    p_geom <- ks_bootstrap_p(kept, function(t) geom_cdf(sort(unique(t)), gm$theta, kmin),
                             "Geom", list(theta=gm$theta, kmin=kmin), B, kmin)
    out[[i]] <- data.frame(f=f, kmin=kmin, n_kept=n_kept,
                           p_dLN=p_dln, p_Zipf=p_zipf, p_Geom=p_geom)
  }
  bind_rows(out)
}
choose_frontier_dln <- function(scan_df){
  ok <- scan_df %>% filter(!is.na(p_dLN), p_dLN >= 0.05, n_kept >= MIN_N_KEPT)
  if(nrow(ok) > 0) max(ok$f) else {
    avail <- scan_df %>% filter(n_kept >= MIN_N_KEPT)
    if(nrow(avail) == 0) NA_real_ else max(avail$f, na.rm = TRUE)
  }
}

# ---- plotting helpers
make_scan_plot_multimodel <- function(scan_df, title="", show_legend=FALSE){
  df <- scan_df %>% mutate(f_pct=100*f) %>%
    pivot_longer(cols=c("p_dLN","p_Zipf","p_Geom"),
                 names_to="model", values_to="pval") %>%
    mutate(model = dplyr::recode(model, p_dLN="dLN", p_Zipf="Zipf", p_Geom="Geom")) %>%
    filter(is.finite(pval))
  p <- ggplot(df, aes(f_pct, pval, color=model)) +
    geom_hline(yintercept=0.05, linetype="dashed", color="gray60", linewidth=0.4) +
    geom_line(linewidth=0.6, na.rm=TRUE) +
    geom_point(shape=16, size=1.7, na.rm=TRUE) +
    scale_color_manual(values=c(dLN=COL_DLN, Zipf=COL_ZIPF, Geom=COL_GEOM)) +
    coord_cartesian(ylim=c(0,1)) +
    labs(x="Retained fraction [%]", y="KS bootstrap p", title=title) +
    theme_minimal(base_size=11) +
    theme(panel.grid.minor=element_blank(),
          plot.title=element_text(face="bold", hjust=0))
  if(show_legend){
    p + theme(legend.position=c(0.98,0.98), legend.justification=c(1,1),
              legend.direction="vertical", legend.title=element_blank())
  } else {
    p + theme(legend.position="none")
  }
}

make_fit_plot_dln <- function(x, f_star, title = "",
                              show_legend = FALSE,
                              leg_x = 0.985, leg_y = 0.985, leg_bg = 0.80) {
  if (is.na(f_star)) {
    return(ggplot() + theme_void() +
             annotate("text", x=0.5, y=0.5, label="Insufficient data",
                      size=3.5, colour="gray30"))
  }
  kmin <- kmin_by_fraction(x, f_star)
  kept <- x[x >= kmin]
  
  # empirical CCDF (step)
  emp_vals   <- sort(unique(kept))
  emp_counts <- tabulate(match(kept, emp_vals), nbins = length(emp_vals))
  emp_ccdf   <- 1 - cumsum(emp_counts)/sum(emp_counts) + emp_counts/sum(emp_counts)
  df_emp     <- data.frame(
    x = rep(emp_vals, each = 2)[-1],
    y = rep(emp_ccdf, each = 2)[-1],
    model = "Empirical"
  )
  
  # dLN curve
  fit   <- fit_dln(kept, kmin)
  xgrid <- exp(seq(log(kmin - 0.5), log(max(emp_vals) + 0.5), length.out = 600))
  df_d  <- data.frame(x = xgrid,
                      y = dln_ccdf_smooth(xgrid, fit$mu, fit$sigma, kmin),
                      model = "dLN")
  
  # zoom to supported region
  N_keep <- length(kept)
  y_floor <- max(1e-4, (1/N_keep)/3)
  x_max  <- {
    idx <- which(df_emp$y >= y_floor)
    if (length(idx)) min(max(df_emp$x[idx]) * 1.15, max(emp_vals)) else max(emp_vals)
  }
  
  p <- ggplot() +
    geom_step(data = df_emp, aes(x, y, color = model), direction = "hv", linewidth = 0.8) +
    geom_line(data = df_d,   aes(x, y, color = model), linewidth = 1.0) +
    scale_color_manual(values = c(Empirical = COL_EMP, dLN = COL_DLN), name = NULL) +
    scale_x_log10(limits = c(kmin - 0.5, x_max), expand = c(0, 0)) +
    scale_y_log10(limits = c(y_floor, 1.0),        expand = c(0, 0)) +
    labs(x = "k", y = expression(P(K >= k)), title = title) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0),
          plot.margin = margin(2, 3, 2, 2))
  
  if (show_legend) {
    p + theme(legend.position = c(leg_x, leg_y),
              legend.justification = c(1, 1),
              legend.background = element_rect(fill = scales::alpha("white", leg_bg), color = NA),
              legend.direction = "vertical")
  } else {
    p + theme(legend.position = "none")
  }
}

# Like make_fit_plot_dln(), but maps color to 'model' so a legend can be shown.
fit_plot_dln_with_legend <- function(x, f_star, title = "") {
  if (is.na(f_star)) {
    return(ggplot() + theme_void() +
             annotate("text", x = 0.5, y = 0.5, label = "Insufficient data",
                      size = 3.5, colour = "gray30"))
  }
  kmin <- kmin_by_fraction(x, f_star)
  kept <- x[x >= kmin]
  
  # empirical CCDF (step)
  emp_vals   <- sort(unique(kept))
  emp_counts <- tabulate(match(kept, emp_vals), nbins = length(emp_vals))
  emp_ccdf   <- 1 - cumsum(emp_counts)/sum(emp_counts) + emp_counts/sum(emp_counts)
  df_emp     <- data.frame(
    x = rep(emp_vals, each = 2)[-1],
    y = rep(emp_ccdf, each = 2)[-1],
    model = "Empirical"
  )
  
  # dLN fit
  fit   <- fit_dln(kept, kmin)
  xgrid <- exp(seq(log(kmin - 0.5), log(max(emp_vals) + 0.5), length.out = 600))
  df_d  <- data.frame(
    x = xgrid,
    y = dln_ccdf_smooth(xgrid, fit$mu, fit$sigma, kmin),
    model = "dLN"
  )
  
  # zoom to empirically supported region
  N_keep <- length(kept); y_floor <- max(1e-4, (1 / N_keep) / 3)
  x_max  <- {
    idx <- which(df_emp$y >= y_floor)
    if (length(idx)) min(max(df_emp$x[idx]) * 1.15, max(emp_vals)) else max(emp_vals)
  }
  
  ggplot() +
    geom_step(data = df_emp, aes(x, y, color = model), direction = "hv", linewidth = 0.8) +
    geom_line(data = df_d,   aes(x, y, color = model), linewidth = 1.0) +
    scale_color_manual(values = c(Empirical = COL_EMP, dLN = COL_DLN), name = NULL) +
    scale_x_log10(limits = c(kmin - 0.5, x_max), expand = c(0,0)) +
    scale_y_log10(limits = c(y_floor, 1.0), expand = c(0,0)) +
    labs(x = NULL, y = expression(P(K>=k)), title = title) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          plot.margin = margin(2, 3, 2, 2),
          legend.position = c(0.98, 0.985),     # top-right inside
          legend.justification = c(1, 1),
          legend.background = element_rect(fill = scales::alpha("white", 0.80), color = NA),
          legend.direction = "vertical")
}

# --- header & row strip using ggplot2 rectangles
header_cell <- function(label){
  ggplot() +
    geom_rect(aes(xmin=0, xmax=1, ymin=0, ymax=1), fill="#F3F3F3", color=NA) +
    annotate("text", x=0.5, y=0.5, label=label, fontface="bold", size=HEADER_TXT_SZ) +
    theme_void() +
    theme(plot.margin = margin(0,0,0,0))
}
row_label_strip <- function(label){
  ggplot() + theme_void() +
    geom_rect(aes(xmin=0, xmax=1, ymin=0, ymax=1), fill="#E9E9E9", color=NA) +
    annotate("text", x=0.5, y=0.5, label=label, angle=270, fontface="bold")
}
spacer <- function(){ ggplot() + theme_void() }

# ========================= DATA: LW =========================
LW <- sanitize_names(read.csv("GER_pilots.csv", check.names = TRUE, strip.white = TRUE))
year_cols_LW <- grep("^y[0-9]{4}$", names(LW), value = TRUE)
if (!length(year_cols_LW)) stop("LW: need columns like y1939..y1945 in GER_pilots.csv")
years_LW <- as.integer(sub("^y", "", year_cols_LW))
LW$total <- rowSums(LW[year_cols_LW], na.rm = TRUE)
LW$entry_year <- apply(LW[year_cols_LW], 1, first_nonzero_year, years_vec = years_LW)
LW_targets <- c(1940L, 1941L, 1942L)
LW_sets <- lapply(LW_targets, function(y) as.integer(subset(LW, entry_year==y, select=total)[,1]))
names(LW_sets) <- paste0("", LW_targets)

# ========================= DATA: BIO & CS ====================
derive_first_year <- function(df, regex="^y(199[6-9]|20(0[0-9]|1[0-6]))$"){
  ycols <- grep(regex, names(df), value = TRUE)
  years <- as.integer(sub("^y", "", ycols))
  fy <- apply(df[ycols], 1, first_nonzero_year, years_vec = years)
  total <- rowSums(df[ycols], na.rm = TRUE)
  data.frame(first_year = fy, score = as.integer(total))
}
BIO <- sanitize_names(read.csv("SCHOLARS_BIO_breakdown.csv", check.names = TRUE))
CS  <- sanitize_names(read.csv("SCHOLARS_CS_breakdown.csv",  check.names = TRUE))
BIO_fy <- derive_first_year(BIO) %>% filter(!is.na(first_year), first_year<=2010, score>=2)
CS_fy  <- derive_first_year(CS)  %>% filter(!is.na(first_year),  first_year<=2010, score>=2)

assign_period_3 <- function(y){
  cut(y, breaks=c(1995,2000,2005,2010),
      labels = c("1996–2000","2001–2005","2006–2010"),
      right=TRUE, include.lowest=TRUE)
}
BIO_fy$Period <- assign_period_3(BIO_fy$first_year)
CS_fy$Period  <- assign_period_3(CS_fy$first_year)
BIO_periods <- c("1996–2000","2001–2005","2006–2010")
CS_periods  <- c("1996–2000","2001–2005","2006–2010")
BIO_sets <- lapply(BIO_periods, function(p) BIO_fy %>% filter(Period==p) %>% pull(score) %>% as.integer())
CS_sets  <- lapply(CS_periods,  function(p) CS_fy  %>% filter(Period==p)  %>% pull(score) %>% as.integer())
names(BIO_sets) <- BIO_periods
names(CS_sets)  <- CS_periods

# ========================= DATA: ATHLETES (breakdown) ====================
is_summer_col <- function(nm) grepl("^(X)?\\d{4}\\.(Summer)$", nm, ignore.case = TRUE)
year_from_col <- function(nm) as.integer(sub("^X?(\\d{4})\\..*$", "\\1", nm, perl=TRUE))
derive_period_score_from_breakdown <- function(path){
  df <- sanitize_names(read.csv(path, check.names = FALSE))
  nms <- names(df); summer_cols <- nms[is_summer_col(nms)]
  if(!length(summer_cols)) stop("No 'XYYYY.Summer' columns in: ", path)
  yrs <- year_from_col(summer_cols); o <- order(yrs); summer_cols <- summer_cols[o]; yrs <- yrs[o]
  df[summer_cols] <- lapply(df[summer_cols], function(v) suppressWarnings(as.numeric(v)))
  FirstYear <- apply(df[summer_cols], 1, first_nonzero_year, years_vec = yrs)
  Score <- rowSums(df[summer_cols], na.rm = TRUE)
  out <- data.frame(FirstYear=FirstYear, Score=as.integer(Score))
  out <- out %>% filter(!is.na(FirstYear), FirstYear>=1896, FirstYear<=2016, Score>=2)
  out$Period <- cut(out$FirstYear,
                    breaks = c(1895,1912,1936,1976,2016),
                    labels = c("1896–1912","1920–1936","1948–1976","1980–2016"),
                    right = TRUE, include.lowest = TRUE)
  out
}
US_df <- derive_period_score_from_breakdown("ATHLETES_USA_SWIMMERS_breakdown.csv")
FR_df <- derive_period_score_from_breakdown("ATHLETES_FRA_FENCING_breakdown.csv")
US_1980 <- US_df %>% filter(Period=="1980–2016") %>% pull(Score) %>% as.integer()
FR_1980 <- FR_df %>% filter(Period=="1980–2016") %>% pull(Score) %>% as.integer()
US_1980 <- US_1980[US_1980>=2]; FR_1980 <- FR_1980[FR_1980>=2]

# ========================= Figure builders ============================

pair_for_subset <- function(x, want_legend=FALSE, left_title="", right_title=""){
  x <- as.integer(x); x <- x[!is.na(x) & x>=2]
  scn <- scan_one_multimodel(x)
  f_star <- choose_frontier_dln(scn)
  pL <- make_scan_plot_multimodel(scn, title = left_title, show_legend = want_legend)
  pR <- make_fit_plot_dln(x, f_star, title = right_title)
  list(left = pL, right = pR)
}

# ----- 1) LW 2×3 (scan | fit) with row strips -----
lw_rows <- list()
for (i in seq_along(LW_sets)) {
  lab <- names(LW_sets)[i]
  x_i <- LW_sets[[i]]
  
  # build scan + choose frontier
  scn_i <- scan_one_multimodel(x_i)
  f_i   <- choose_frontier_dln(scn_i)
  
  # SCAN: show the dLN/Zipf/Geom legend only on the first row's first panel
  left <- make_scan_plot_multimodel(
    scn_i,
    title = "",
    show_legend = (i == 1)   # << legend on row 1, col 1
  )
  
  # FIT: keep your existing behavior (legend in row 1 fit if desired)
  right <- if (i == 1) {
    fit_plot_dln_with_legend(x_i, f_i, title = "")  # Empirical / dLN legend
  } else {
    make_fit_plot_dln(x_i, f_i, title = "")
  }
  
  row_main       <- plot_grid(left, right, ncol = 2, align = "hv")
  row_with_strip <- plot_grid(row_main, row_label_strip(lab), ncol = 2, rel_widths = c(1, 0.05))
  lw_rows[[i]]   <- row_with_strip
}
lw_body <- plot_grid(plotlist = lw_rows, ncol = 1, align = "v")
lw_header <- plot_grid(header_cell("LW-Frontier scan"),
                       header_cell("LW-dLN fit at frontier"),
                       ncol = 2, rel_widths = c(1,1))
lw_fig <- plot_grid(lw_header, lw_body, ncol = 1, rel_heights = c(0.06, 1))
ggsave("Figure_LW_periods_2x3.pdf", lw_fig, width = 14, height = 10, units = "in", device = cairo_pdf)
message("Wrote: Figure_LW_periods_2x3.pdf")

# ----- 2) BIO + CS 4×3 with scan–scan–fit–fit and row strips -----
bio_cs_rows <- list()
for (i in seq_along(BIO_sets)) {
  per_lab <- names(BIO_sets)[i]
  
  # BIO pair (no legend)
  bp <- pair_for_subset(BIO_sets[[i]], want_legend = FALSE,
                        left_title = "", right_title = "")
  
  # CS pair (scan legend only on first row)
  want_leg_cs <- SHOW_LEGEND_BIO_CS && i == 1
  cp <- pair_for_subset(CS_sets[[i]], want_legend = want_leg_cs,
                        left_title = "", right_title = "")
  
  ## NEW: put Empirical/dLN legend in the top-right (CS fit) panel on the first row
  if (i == 1) {
    scn_cs <- scan_one_multimodel(CS_sets[[i]])
    f_cs   <- choose_frontier_dln(scn_cs)
    cp$right <- make_fit_plot_dln(
      CS_sets[[i]], f_cs, title = "",
      show_legend = TRUE, leg_x = 0.985, leg_y = 0.985, leg_bg = 0.80
    )
    # (Alternatively: cp$right <- fit_plot_dln_with_legend(CS_sets[[i]], f_cs, title = ""))
  }
  
  row_main <- plot_grid(
    bp$left,  cp$left,
    spacer(),
    bp$right, cp$right,
    ncol = 5, rel_widths = c(1, 1, 0.04, 1, 1), align = "hv"
  )
  row_with_strip <- plot_grid(row_main, row_label_strip(per_lab),
                              ncol = 2, rel_widths = c(1, 0.05))
  bio_cs_rows[[i]] <- row_with_strip
}

bio_cs_body <- plot_grid(plotlist = bio_cs_rows, ncol = 1, align = "v")

# New 4-column headers (plus the spacer column)
bio_cs_header <- plot_grid(
  header_cell("BIO - Frontier scan"),
  header_cell("CS - Frontier scan"),
  spacer(),
  header_cell("BIO - dLN fit at frontier"),
  header_cell("CS - dLN fit at frontier"),
  ncol = 5, rel_widths = c(1, 1, 0.04, 1, 1)
)

bio_cs_fig <- plot_grid(bio_cs_header, bio_cs_body, ncol = 1, rel_heights = c(0.06, 1))
ggsave("Figure_BIO_CS_periods_4x3_scans_scans_fits_fits.pdf",
       bio_cs_fig, width = 16, height = 12, units = "in")

# ----- 3) Athletes 1×4 (1980–2016 only) scan–scan–fit–fit -----
scn_US <- scan_one_multimodel(US_1980)
scn_FR <- scan_one_multimodel(FR_1980)
f_US   <- choose_frontier_dln(scn_US)
f_FR   <- choose_frontier_dln(scn_FR)

p_scan_US <- make_scan_plot_multimodel(scn_US, title = NULL,
                                       show_legend = FALSE)
p_scan_FR <- make_scan_plot_multimodel(scn_FR, title = NULL,
                                       show_legend = TRUE)  # keep scan legend off here

# FITS: legend only in the rightmost panel (FRAOT_F)
p_fit_US  <- make_fit_plot_dln(US_1980, f_US, title = NULL,
                               show_legend = FALSE)
p_fit_FR  <- make_fit_plot_dln(FR_1980, f_FR, title = NULL,
                               show_legend = TRUE,  # turn legend on here
                               leg_x = 0.985, leg_y = 0.985, leg_bg = 0.80)

# --- Boxed top titles (gray) for each column
ath_header <- plot_grid(
  header_cell("USOT_S — Frontier scan"),
  header_cell("FRAOT_F — Frontier scan"),
  header_cell("USOT_S — dLN fit at frontier"),
  header_cell("FRAOT_F — dLN fit at frontier"),
  ncol = 4, rel_widths = c(1, 1, 1, 1)
)

# --- Panels row (4 panels) + a gray side strip with the period (right side)
ath_panels <- plot_grid(p_scan_US, p_scan_FR, p_fit_US, p_fit_FR, ncol = 4, align = "hv")

ath_row <- plot_grid(
  ath_panels,
  row_label_strip("1980–2016"),
  ncol = 2, rel_widths = c(1, SIDE_COL_W_ATH), align = "h"
)

# --- Stack header + row; HEADER_BAND_H_ATH controls the height of the gray header band
ath_fig <- plot_grid(
  ath_header,
  ath_row,
  ncol = 1,
  rel_heights = c(HEADER_BAND_H_ATH, 1)
)

ggsave("Figure_Athletes_1980_2016_1x4_scans_scans_fits_fits.pdf",
       ath_fig, width = 16, height = 4.8, units = "in")

# ========================= PERIOD MODEL ADEQUACY TABLE =========================
# Helper: recompute p-values for all 3 models at a given frontier fraction
pvals_at_frontier <- function(x, f_star){
  if (is.na(f_star)) return(list(kmin=NA_integer_, n_kept=0,
                                 p_dLN=NA_real_, p_Zipf=NA_real_, p_Geom=NA_real_))
  kmin <- kmin_by_fraction(x, f_star)
  kept <- x[x >= kmin]
  if (length(kept) < MIN_N_KEPT || length(unique(kept)) < 2L) {
    return(list(kmin=kmin, n_kept=length(kept),
                p_dLN=NA_real_, p_Zipf=NA_real_, p_Geom=NA_real_))
  }
  # fit params
  dln <- fit_dln(kept, kmin)
  zf  <- fit_zipf_finite(kept, kmin)
  gm  <- fit_geom(kept, kmin)
  # KS bootstrap p (same B as scans)
  p_dln <- ks_bootstrap_p(
    kept, function(t) dln_cdf(sort(unique(t)), dln$mu, dln$sigma, kmin),
    "dLN",  list(mu=dln$mu, sigma=dln$sigma, kmin=kmin), B_BOOT, kmin
  )
  p_zip <- ks_bootstrap_p(
    kept, function(t) zipf_cdf(zf, sort(unique(t))),
    "Zipf", list(zf_fit=zf), B_BOOT, kmin
  )
  p_geo <- ks_bootstrap_p(
    kept, function(t) geom_cdf(sort(unique(t)), gm$theta, kmin),
    "Geom", list(theta=gm$theta, kmin=kmin), B_BOOT, kmin
  )
  list(kmin=kmin, n_kept=length(kept), p_dLN=p_dln, p_Zipf=p_zip, p_Geom=p_geo)
}

summ_row <- function(domain, period, x_vec){
  x <- as.integer(x_vec); x <- x[!is.na(x) & x>=2]
  Ntot <- length(x)
  if (Ntot == 0) {
    return(data.frame(
      domain=domain, period=period, N_total=0,
      f_star=NA_real_, kmin=NA_integer_, n_kept=0,
      p_dLN=NA_real_, p_Zipf=NA_real_, p_Geom=NA_real_,
      adequate_dLN=FALSE, adequate_Zipf=FALSE, adequate_Geom=FALSE,
      winner_by_p="None", stringsAsFactors = FALSE
    ))
  }
  scn <- scan_one_multimodel(x)
  f_star <- choose_frontier_dln(scn)  # same rule used for plots
  pv <- pvals_at_frontier(x, f_star)
  pvec <- c(dLN=pv$p_dLN, Zipf=pv$p_Zipf, Geom=pv$p_Geom)
  adequate <- pvec >= 0.05
  # choose winner by max p (ties -> first by order dLN, Zipf, Geom; all NA -> None)
  if (all(!is.finite(pvec))) {
    win <- "None"
  } else {
    mx <- which.max(replace(pvec, !is.finite(pvec), -Inf))
    win <- names(pvec)[mx]
  }
  data.frame(
    domain=domain, period=period, N_total=Ntot,
    f_star=round(f_star, 3), kmin=pv$kmin, n_kept=pv$n_kept,
    p_dLN=round(pv$p_dLN, 4), p_Zipf=round(pv$p_Zipf, 4), p_Geom=round(pv$p_Geom, 4),
    adequate_dLN = isTRUE(adequate["dLN"]),
    adequate_Zipf= isTRUE(adequate["Zipf"]),
    adequate_Geom= isTRUE(adequate["Geom"]),
    winner_by_p = win,
    stringsAsFactors = FALSE
  )
}

tbl_rows <- list()

# LW periods
for(i in seq_along(LW_sets)){
  tbl_rows[[length(tbl_rows)+1]] <- summ_row("LW", names(LW_sets)[i], LW_sets[[i]])
}
# BIO periods
for(i in seq_along(BIO_sets)){
  tbl_rows[[length(tbl_rows)+1]] <- summ_row("BIO", names(BIO_sets)[i], BIO_sets[[i]])
}
# CS periods
for(i in seq_along(CS_sets)){
  tbl_rows[[length(tbl_rows)+1]] <- summ_row("CS", names(CS_sets)[i], CS_sets[[i]])
}
# Athletes (all four periods; many early will be small/NA and that’s OK to show)
US_periods_all  <- c("1896–1912","1920–1936","1948–1976","1980–2016")
FR_periods_all  <- c("1896–1912","1920–1936","1948–1976","1980–2016")

US_sets_all <- lapply(US_periods_all, function(p)
  US_df %>% dplyr::filter(Period==p) %>% dplyr::pull(Score) %>% as.integer())
FR_sets_all <- lapply(FR_periods_all, function(p)
  FR_df %>% dplyr::filter(Period==p) %>% dplyr::pull(Score) %>% as.integer())

for(i in seq_along(US_sets_all)){
  tbl_rows[[length(tbl_rows)+1]] <- summ_row("USOT_S", US_periods_all[i], US_sets_all[[i]])
}
for(i in seq_along(FR_sets_all)){
  tbl_rows[[length(tbl_rows)+1]] <- summ_row("FRAOT_F", FR_periods_all[i], FR_sets_all[[i]])
}

Model_Adequacy_Period <- dplyr::bind_rows(tbl_rows)

# Optional: nice ordering
Model_Adequacy_Period <- Model_Adequacy_Period %>%
  dplyr::mutate(domain = factor(domain, levels=c("LW","BIO","CS","USOT_S","FRAOT_F"))) %>%
  dplyr::arrange(domain, period)

write.csv(Model_Adequacy_Period, "Period_Model_Adequacy_Table.csv", row.names = FALSE)
message("Wrote: Period_Model_Adequacy_Table.csv (", nrow(Model_Adequacy_Period), " rows)")