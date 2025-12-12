# GBROT_S_Composite_Figure_ALT_regular_grid.R
# Composite for GBROT_S swimmers vs USOT_S.
# Panel A: scan on a REGULAR fraction grid (30..100% by 5), evaluated at nearest feasible kmin (f_actual >= f_target).
# Panels B–D: fits and ΔAIC at locked ~80% frontier (actual shown).
# Uses DISCRETE models + KS bootstrap p, Unicode labels, cairo_pdf device.

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(cowplot); library(scales)
})

# ---------------- Files ----------------
FILE_GBR_TOTAL      <- "ATHLETES_GBR_SWIMMERS_total.csv"
FILE_USA_BREAKDOWN  <- "ATHLETES_USA_SWIMMERS_breakdown.csv"

FORCE_COL_GBR <- "Total_Score"   # explicit, avoids guessing

# ---------------- Settings -------------
# Regular fraction grid for the scan
FRACTION_GRID <- seq(0.30, 0.80, by = 0.05)

# Minimum n to attempt KS p at a scan point (after truncation)
SCAN_MINN    <- 8

# AIC / bootstrap controls
AIC_MINN     <- 10
B_BOOT       <- 120      # ↑ for camera-ready
R_DS         <- 300
BASE_SEED    <- 202508

# Frontier used in Panels B–D (locked)
FORCE_F_COMMON <- 0.80   # 80%

# Fonts (set to "" to use system default if DejaVu Sans is unavailable)
BASE_FAMILY <- "DejaVu Sans"

# Colors
COL_EMP  <- "#4D4D4D"
COL_DLN  <- "#D62728"
COL_ZIPF <- "#1F77B4"
COL_GEOM <- "#2CA02C"

# -------------- Utilities -------------
extract_scores <- function(path, force_col=NULL, min_k=2L){
  df <- read.csv(path, check.names=FALSE)
  if(!is.null(force_col)){
    if(!force_col %in% names(df)) stop(sprintf("force_col='%s' not in %s", force_col, path))
    pick <- force_col
  } else {
    nms <- tolower(names(df))
    cand <- c("score","x","total","points","medal_points","total_score","victories","kills","hits")
    hit <- which(nms %in% cand)
    if(length(hit)==0){
      num <- which(sapply(df, is.numeric))
      if(!length(num)) stop("No numeric columns in: ", path)
      rng <- sapply(df[num], function(v) diff(range(v, na.rm=TRUE)))
      pick <- names(df)[num[which.max(rng)]]
    } else pick <- names(df)[hit[1]]
    message(sprintf("Using '%s' from %s as score column.", pick, path))
  }
  x <- suppressWarnings(as.integer(df[[pick]]))
  x <- x[!is.na(x) & x >= min_k]
  if(!length(x)) stop(sprintf("%s has no observations ≥ %d in '%s'", path, min_k, pick))
  x
}

Phi <- function(z) pnorm(z)

# ---- discrete lognormal (truncated at kmin)
dln_loglik <- function(par, k, kmin){
  mu <- par[1]; sigma <- max(1e-3, abs(par[2]))
  a  <- (log(kmin - 0.5) - mu)/sigma
  Zd <- max(1e-16, 1 - Phi(a))
  bH <- (log(k + 0.5) - mu)/sigma
  bL <- (log(k - 0.5) - mu)/sigma
  pm <- pmax(Phi(bH) - Phi(bL), 0) / Zd
  sum(log(pmax(pm, 1e-16)))
}
fit_dln <- function(k, kmin){
  kk <- as.integer(k[k >= kmin]); if(length(kk) < 2) stop("dLN: n<2 at kmin.")
  logs <- log(pmax(as.numeric(kk), 2))
  mu0 <- mean(logs); if(!is.finite(mu0)) mu0 <- median(logs)
  s0  <- sd(logs);   if(!is.finite(s0) || s0<1e-3) s0 <- max(0.5, stats::mad(logs)*1.0)
  mu_lo <- mu0 - 5; mu_hi <- mu0 + 5
  nll <- function(p){ if(p[2]<=0) return(1e-12 + 1e12); val <- -dln_loglik(p, kk, kmin); if(!is.finite(val)) 1e12 else val }
  opt <- optim(c(mu0,s0), fn=nll, method="L-BFGS-B", lower=c(mu_lo,0.05), upper=c(mu_hi,5))
  list(mu=opt$par[1], sigma=max(0.05, abs(opt$par[2])))
}
dln_cdf <- function(kq, mu, sigma, kmin){
  kq <- as.integer(kq); J <- seq.int(kmin, max(kq))
  a <- (log(kmin-0.5)-mu)/sigma; Zd <- max(1e-16, 1-Phi(a))
  bH <- (log(J+0.5)-mu)/sigma; bL <- (log(J-0.5)-mu)/sigma
  pm <- pmax(Phi(bH)-Phi(bL),0)/Zd; Fj <- cumsum(pm)
  Fj[pmax(1, findInterval(kq, J))]
}
dln_ccdf_smooth <- function(x, mu, sigma, kmin){
  den <- max(1e-12, 1 - Phi((log(kmin-0.5)-mu)/sigma))
  pmax((1 - Phi((log(x)-mu)/sigma))/den, 1e-12)
}
rdln <- function(n, mu, sigma, kmin){
  out <- integer(0); thr <- kmin-0.5
  while(length(out)<n){
    m <- ceiling((n-length(out))*1.6)
    x <- rlnorm(m, meanlog=mu, sdlog=sigma); x <- x[x>=thr]
    if(length(x)) out <- c(out, as.integer(round(x)))
  }
  out[seq_len(n)]
}

# ---- finite Zipf (discrete power law)
fit_zipf_finite <- function(k, kmin){
  kk <- as.integer(k[k>=kmin]); ks <- seq.int(kmin, max(kk))
  logs <- sum(log(kk)); n <- length(kk)
  nll <- function(a){ if(a<=1.0001) return(1e12); w <- ks^(-a); a*logs + n*log(sum(w)) }
  opt <- optimize(nll, interval=c(1.001,10), maximum=FALSE); a <- opt$minimum
  w <- ks^(-a); pm <- w/sum(w); cdf <- cumsum(pm)
  list(alpha=a, ks=ks, pmf=pm, cdf=cdf)
}
zipf_cdf <- function(fit, kq){ idx <- pmax(1, findInterval(as.integer(kq), fit$ks)); fit$cdf[pmin(idx, length(fit$cdf))] }
rzipf_finite <- function(n, fit){ sample(fit$ks, size=n, replace=TRUE, prob=fit$pmf) }

# ---- shifted geometric
fit_geom <- function(k, kmin){
  kk <- as.integer(k[k>=kmin]); y <- kk-kmin; ybar <- mean(y)
  list(theta = if(!is.na(ybar) && ybar>0) ybar/(1+ybar) else 1e-4)
}
geom_cdf <- function(kq, theta, kmin){ y <- (as.integer(kq)-kmin)+1L; pmin(pmax(1 - (theta^y),0),1) }
rgeom_shift <- function(n, theta, kmin){ p <- max(1e-12, 1-theta); kmin + rgeom(n, prob=p) }

# ---- discrete KS + bootstrap
ks_stat_discrete <- function(k, F_model){
  k <- sort(as.integer(k)); u <- sort(unique(k))
  counts <- tabulate(match(k, u), nbins=length(u))
  ecdf_u <- cumsum(counts)/sum(counts); modF <- F_model(u)
  max(abs(ecdf_u - modF))
}
ks_bootstrap_p <- function(kept, cdf_fn, model, params, B=120L, seed_offset=0L){
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

# ---- AIC
aic_dln  <- function(kept, kmin, mu, sigma){ ll <- dln_loglik(c(mu,sigma), kept, kmin); 2*2 - 2*ll }
aic_zipf <- function(kept, zf){ idx <- pmax(1, findInterval(kept, zf$ks)); ll <- sum(log(pmax(zf$pmf[pmin(idx, length(zf$pmf))],1e-16))); 2*1 - 2*ll }
aic_geom <- function(kept, kmin, th){ y <- kept-kmin; pm <- (1-th)*(th^y); ll <- sum(log(pmax(pm,1e-16))); 2*1 - 2*ll }

smooth_ccdf_from_discrete <- function(x, ccdf, n=800L){
  lx <- log(as.numeric(x)); ly <- log(pmax(ccdf,1e-16))
  xg <- seq(min(lx), max(lx), length.out=n); yg <- approx(lx, ly, xout=xg, rule=2)$y
  data.frame(x=exp(xg), y=exp(yg))
}

make_totals_from_breakdown <- function(path, min_k=2L){
  df <- read.csv(path, check.names = FALSE)
  num_cols <- sapply(df, is.numeric)
  # drop ID-like columns if needed
  if("name_UID" %in% names(df)) num_cols[which(names(df) == "name_UID")] <- FALSE
  totals <- rowSums(df[, num_cols, drop = FALSE])
  x <- as.integer(totals)
  x <- x[!is.na(x) & x >= min_k]
  if(!length(x)) stop(sprintf("%s has no totals ≥ %d", path, min_k))
  x
}

# ---------------- Load data ----------------
GB <- extract_scores(FILE_GBR_TOTAL, FORCE_COL_GBR, min_k = 2L)
US <- make_totals_from_breakdown(FILE_USA_BREAKDOWN, min_k = 2L)

message(sprintf("GB N = %d, US N = %d (totals ≥ 2)", length(GB), length(US)))

# Threshold table over all unique kmin (feasible cut points)
vals <- sort(unique(GB))
cand_all <- data.frame(
  kmin   = vals,
  n_kept = sapply(vals, function(k) sum(GB >= k))
)
cand_all$f <- cand_all$n_kept / length(GB)

# Map target f to the closest feasible kmin with f >= target; else fallback to max feasible f below target
kmin_for_f <- function(target_f){
  i <- which(cand_all$f >= target_f)
  if(length(i)==0){ j <- which.max(cand_all$f); return(cand_all$kmin[j]) }
  cand_all$kmin[max(i)]
}

# ---------------- Panel A: REGULAR GRID scan ----------------------
scan <- data.frame(f_target = numeric(0), f_actual = numeric(0),
                   p_dLN = numeric(0), p_Zipf = numeric(0), p_Geom = numeric(0),
                   kmin = integer(0), n_kept = integer(0))

for (f_tgt in FRACTION_GRID) {
  kmin  <- kmin_for_f(f_tgt)
  kept  <- GB[GB >= kmin]
  f_act <- length(kept) / length(GB)
  n_kept <- length(kept)
  
  p_dln <- p_zipf <- p_geom <- NA_real_
  if (n_kept >= SCAN_MINN && length(unique(kept)) >= 2L) {
    dln <- fit_dln(kept, kmin)
    zf  <- fit_zipf_finite(kept, kmin)
    gm  <- fit_geom(kept, kmin)
    
    p_dln <- ks_bootstrap_p(kept, function(t) dln_cdf(sort(unique(t)), dln$mu, dln$sigma, kmin),
                            "dLN", list(mu=dln$mu, sigma=dln$sigma, kmin=kmin), B_BOOT, round(100*f_tgt))
    p_zipf <- ks_bootstrap_p(kept, function(t) zipf_cdf(zf, sort(unique(t))),
                             "Zipf", list(zf_fit=zf), B_BOOT, round(100*f_tgt))
    p_geom <- ks_bootstrap_p(kept, function(t) geom_cdf(sort(unique(t)), gm$theta, kmin),
                             "Geom", list(theta=gm$theta, kmin=kmin), B_BOOT, round(100*f_tgt))
  }
  
  scan <- rbind(scan, data.frame(
    f_target = f_tgt, f_actual = f_act, p_dLN = p_dln, p_Zipf = p_zipf, p_Geom = p_geom,
    kmin = kmin, n_kept = n_kept
  ))
}

# ---------------- Locked frontier at 80% ---------------------------
kmin_common <- kmin_for_f(FORCE_F_COMMON)
f_common    <- sum(GB >= kmin_common)/length(GB)
kept_common <- GB[GB >= kmin_common]
if(length(unique(kept_common)) < 2L) stop("Locked frontier leaves <2 unique values.")

# ---------------- Fits at frontier -------------------------------
dln_c  <- fit_dln(kept_common, kmin_common)
zipf_c <- fit_zipf_finite(kept_common, kmin_common)
geom_c <- fit_geom(kept_common, kmin_common)

emp_vals   <- sort(unique(kept_common))
emp_counts <- tabulate(match(kept_common, emp_vals), nbins=length(emp_vals))
emp_ccdf   <- 1 - cumsum(emp_counts)/sum(emp_counts) + emp_counts/sum(emp_counts)

df_emp  <- data.frame(x=emp_vals, y=emp_ccdf, model="Empirical")
df_dln  <- { xg <- exp(seq(log(kmin_common-0.5), log(max(emp_vals)+0.5), length.out=600))
data.frame(x=xg, y=dln_ccdf_smooth(xg, dln_c$mu, dln_c$sigma, kmin_common), model="dLN") }
df_zipf <- { d <- smooth_ccdf_from_discrete(zipf_c$ks, 1 - zipf_c$cdf); d$model <- "Zipf"; d }
geom_x    <- seq.int(kmin_common, max(emp_vals))
geom_ccdf <- 1 - geom_cdf(geom_x, geom_c$theta, kmin_common)
df_geom   <- { d <- smooth_ccdf_from_discrete(geom_x, geom_ccdf); d$model <- "Geom"; d }

N_keep  <- length(kept_common); y_floor <- max(1e-4, (1/N_keep)/3)
x_max_zoom <- { idx <- which(emp_ccdf >= y_floor)
if(length(idx)) min(max(emp_vals[idx])*1.15, max(emp_vals)) else max(emp_vals) }

# ---------------- Downsample ΔAIC --------------------------
set.seed(BASE_SEED)
n_target <- length(kept_common)
gaps_zipf_minus_dln <- gaps_geom_minus_dln <- numeric(0)

for(r in seq_len(R_DS)){
  samp <- sample(US, size=n_target, replace=(n_target > length(US)))
  kept <- samp[samp >= kmin_common]
  if(length(kept) < AIC_MINN) next
  dln_r <- fit_dln(kept, kmin_common); zf_r <- fit_zipf_finite(kept, kmin_common); gm_r <- fit_geom(kept, kmin_common)
  a_dln <- aic_dln(kept, kmin_common, dln_r$mu, dln_r$sigma)
  a_zip <- aic_zipf(kept, zf_r); a_geo <- aic_geom(kept, kmin_common, gm_r$theta)
  gaps_zipf_minus_dln <- c(gaps_zipf_minus_dln, a_zip - a_dln)   # >0 ⇒ dLN better
  gaps_geom_minus_dln <- c(gaps_geom_minus_dln, a_geo - a_dln)   # >0 ⇒ dLN better
}

# ---------------- Plots ------------------------------
# Use target fraction on x for a dense/regular-looking scan
dfA <- scan %>%
  transmute(f_pct = 100*f_target, dLN=p_dLN, Zipf=p_Zipf, Geom=p_Geom) %>%
  pivot_longer(cols=c("dLN","Zipf","Geom"), names_to="model", values_to="pval")

pA <- ggplot(dfA, aes(f_pct, pval, color=model)) +
  geom_hline(yintercept=0.05, linetype="dashed", color="gray60", linewidth=0.4) +
  #geom_vline(xintercept=80, linetype="dotted", color="gray50") +  # show locked frontier
  geom_line(na.rm=TRUE, linewidth=0.6) +
  geom_point(shape=16, size=1.9, na.rm=TRUE) +
  scale_color_manual(values=c(dLN=COL_DLN, Zipf=COL_ZIPF, Geom=COL_GEOM)) +
  scale_x_continuous(limits=c(30, 80), breaks=seq(30,80,10)) +
  coord_cartesian(ylim=c(0,0.9)) +
  labs(x="Retained fraction [%]", y="KS bootstrap p", title="A.  Frontier scan (GBROT_S)") +
  theme_minimal(base_size=11, base_family = BASE_FAMILY) +
  theme(legend.position="none", panel.grid.minor=element_blank(),
        plot.title=element_text(face="bold", hjust=0))

pB <- ggplot() +
  geom_step(data=df_emp,  aes(x,y,color=model), direction="hv", linewidth=0.8) +
  geom_line(data=df_dln,  aes(x,y,color=model), linewidth=1.0) +
  geom_line(data=df_zipf, aes(x,y,color=model), linewidth=0.9) +
  geom_line(data=df_geom, aes(x,y,color=model), linewidth=0.9) +
  scale_color_manual(values=c(Empirical=COL_EMP, dLN=COL_DLN, Zipf=COL_ZIPF, Geom=COL_GEOM)) +
  scale_x_log10(limits=c(kmin_common-0.5, x_max_zoom), expand=c(0,0)) +
  scale_y_log10(limits=c(y_floor, 1.0), expand=c(0,0)) +
  labs(x="k", y="P(K ≥ k)",
       title=sprintf("B.  Model fits at common frontier (GBROT_S)", 100*f_common)) +
  theme_minimal(base_size=11, base_family = BASE_FAMILY) +
  theme(panel.grid.minor=element_blank(),
        legend.position=c(0.98,0.98), legend.justification=c(1,1),
        legend.direction="vertical", legend.title=element_blank(),
        plot.title=element_text(face="bold", hjust=0))

dfC <- data.frame(gap=gaps_zipf_minus_dln)
pC <- ggplot(dfC, aes(gap)) +
  geom_histogram(bins=24, fill=alpha(COL_ZIPF,0.25), color=COL_ZIPF, linewidth=0.6, na.rm=TRUE) +
  geom_vline(xintercept=0, color="black") +
  { if(nrow(dfC)>0) geom_vline(xintercept=mean(dfC$gap, na.rm=TRUE), color=COL_DLN, linetype="dashed") } +
  scale_y_continuous(limits=c(0,50), breaks=seq(0,50,10), labels=seq(0,50,10)) +
  labs(x = "ΔAIC (Zipf − dLN) [>0 → dLN better]", y = "Count",
       title = "C.  USOT_S downsampled to GBROT_S size") +
  theme_minimal(base_size=11, base_family = BASE_FAMILY) +
  theme(panel.grid.minor=element_blank(), plot.title=element_text(face="bold", hjust=0))

dfD <- data.frame(gap=gaps_geom_minus_dln)
pD <- ggplot(dfD, aes(gap)) +
  geom_histogram(bins=24, fill=alpha(COL_GEOM,0.25), color=COL_GEOM, linewidth=0.6, na.rm=TRUE) +
  geom_vline(xintercept=0, color="black") +
  { if(nrow(dfD)>0) geom_vline(xintercept=mean(dfD$gap, na.rm=TRUE), color=COL_DLN, linetype="dashed") } +
  scale_y_continuous(limits=c(0,50)) +
  labs(x = "ΔAIC (Geom − dLN) [>0 → dLN better]", y = NULL,
       title = "D.  USOT_S downsampled to GBROT_S size") +
  theme_minimal(base_size=11, base_family = BASE_FAMILY) +
  theme(panel.grid.minor=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        plot.title=element_text(face="bold", hjust=0))

final <- plot_grid(pA, pB, pC, pD, ncol=2, rel_heights=c(1,1))

# Unicode-safe export (Δ, →, ≥), avoids Type-3 fonts
ggsave("GBROT_S_Swimmers_Composite_FROM_CODE_ALT_regular_grid.pdf",
       final, width = 12.5, height = 9.3, units = "in", device = cairo_pdf)
ggsave("GBROT_S_Swimmers_Composite_FROM_CODE_ALT_regular_grid.png",
       final, width = 12.5, height = 9.3, units = "in")

message(sprintf("Saved ALT figure. Locked frontier target 80%%; actual ~%.0f%% (kmin=%d, n_kept=%d).",
                100*f_common, kmin_common, length(kept_common)))

# ===== Save dLN "win rates" and gap summaries to CSV =========================
summarize_gaps <- function(gaps, label){
  v <- gaps[is.finite(gaps)]
  if (length(v) == 0) {
    return(data.frame(
      comparison = label, n_reps = 0, n_valid = 0,
      pct_dLN_wins = NA_real_, pct_other_wins = NA_real_, pct_ties = NA_real_,
      mean_gap = NA_real_, median_gap = NA_real_, sd_gap = NA_real_,
      q05 = NA_real_, q25 = NA_real_, q75 = NA_real_, q95 = NA_real_
    ))
  }
  qs <- quantile(v, probs = c(0.05, 0.25, 0.75, 0.95), names = FALSE, type = 7)
  data.frame(
    comparison     = label,
    n_reps         = length(gaps),
    n_valid        = length(v),
    pct_dLN_wins   = round(mean(v > 0) * 100, 2),
    pct_other_wins = round(mean(v < 0) * 100, 2),
    pct_ties       = round(mean(v == 0) * 100, 2),
    mean_gap       = round(mean(v), 3),
    median_gap     = round(median(v), 3),
    sd_gap         = round(sd(v), 3),
    q05            = round(qs[1], 3),
    q25            = round(qs[2], 3),
    q75            = round(qs[3], 3),
    q95            = round(qs[4], 3)
  )
}

summary_tbl <- rbind(
  summarize_gaps(gaps_zipf_minus_dln, "Zipf − dLN   (>0 ⇒ dLN better)"),
  summarize_gaps(gaps_geom_minus_dln, "Geom − dLN   (>0 ⇒ dLN better)")
)

out_csv <- "GBROT_S_vs_USOT_S_downsample_dLN_winrates_ALT_regular_grid.csv"
write.csv(summary_tbl, out_csv, row.names = FALSE)
message("Saved: ", out_csv)
# =============================================================================