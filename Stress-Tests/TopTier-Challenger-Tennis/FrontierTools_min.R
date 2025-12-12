## ---------- FrontierTools_min.R (no plotting; just core scan/fit/KS) ----------
## Constants
F_GRID     <- seq(0.20, 0.80, by = 0.05)
MIN_N_KEPT <- 40L
B_BOOT     <- 200L
BASE_SEED  <- 202508L
Phi <- function(z) pnorm(z)

## Utilities
kmin_by_fraction <- function(k, f){
  k <- sort(as.integer(k))
  idx <- ceiling((1 - f) * length(k))
  idx <- max(1L, min(idx, length(k)))
  as.integer(k[idx])
}

## ---- dLN (discrete lognormal) ----
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
  kk <- as.integer(k[k>=kmin]); if(length(kk) < 2L) stop("dLN: n<2 at kmin.")
  logs <- log(as.numeric(kk))
  mu0 <- mean(logs); s0 <- sd(logs); if(!is.finite(s0) || s0 < 1e-3) s0 <- 0.5
  opt <- optim(c(mu0, s0), fn=function(p) -dln_loglik(p, kk, kmin),
               method="L-BFGS-B", lower=c(mu0-5, 0.05), upper=c(mu0+5, 5))
  list(mu=opt$par[1], sigma=max(0.05, abs(opt$par[2])))
}
dln_cdf <- function(kq, mu, sigma, kmin){
  kq <- as.integer(kq); J <- seq.int(kmin, max(kq))
  a  <- (log(kmin-0.5)-mu)/sigma; Zd <- max(1e-16, 1 - Phi(a))
  bH <- (log(J+0.5)-mu)/sigma; bL <- (log(J-0.5)-mu)/sigma
  pm <- pmax(Phi(bH)-Phi(bL), 0) / Zd; Fj <- cumsum(pm)
  Fj[pmax(1, findInterval(kq, J))]
}
rdln <- function(n, mu, sigma, kmin){
  out <- integer(0); thr <- kmin - 0.5
  while (length(out) < n) {
    m <- ceiling((n - length(out))*1.6)
    x <- rlnorm(m, meanlog = mu, sdlog = sigma); x <- x[x >= thr]
    if (length(x)) out <- c(out, as.integer(round(x)))
  }
  out[seq_len(n)]
}

## ---- Zipf (finite discrete power law) ----
fit_zipf_finite <- function(k, kmin){
  kk <- as.integer(k[k>=kmin]); ks <- seq.int(kmin, max(kk))
  logs <- sum(log(kk)); n <- length(kk)
  nll <- function(a){ if(a <= 1.0001) return(1e12); w <- ks^(-a); a*logs + n*log(sum(w)) }
  opt <- optimize(nll, interval=c(1.001, 10), maximum=FALSE); a <- opt$minimum
  w <- ks^(-a); pm <- w/sum(w); cdf <- cumsum(pm)
  list(alpha=a, ks=ks, pmf=pm, cdf=cdf)
}
zipf_cdf <- function(fit, kq){
  idx <- pmax(1, findInterval(as.integer(kq), fit$ks))
  fit$cdf[pmin(idx, length(fit$cdf))]
}
rzipf_finite <- function(n, fit){ sample(fit$ks, size=n, replace=TRUE, prob=fit$pmf) }

## ---- Shifted geometric ----
fit_geom <- function(k, kmin){
  kk <- as.integer(k[k>=kmin]); y <- kk - kmin; ybar <- mean(y)
  list(theta = if(!is.na(ybar) && ybar>0) ybar/(1+ybar) else 1e-4, kmin=kmin)
}
geom_cdf <- function(kq, theta, kmin){
  y <- (as.integer(kq) - kmin) + 1L
  pmin(pmax(1 - (theta^y), 0), 1)
}
rgeom_shift <- function(n, theta, kmin){
  p <- max(1e-12, 1-theta); kmin + rgeom(n, prob=p)
}

## ---- KS and bootstrap ----
ks_stat_discrete <- function(k, F_model){
  k <- sort(as.integer(k)); u <- sort(unique(k))
  counts <- tabulate(match(k, u), nbins = length(u))
  ecdf_u <- cumsum(counts)/sum(counts); modF <- F_model(u)
  max(abs(ecdf_u - modF))
}
ks_bootstrap_p <- function(kept, cdf_fn, sim_draw, B=200L, seed_offset=0L){
  set.seed(BASE_SEED + as.integer(seed_offset))
  n <- length(kept); obs <- ks_stat_discrete(kept, cdf_fn); geq <- 0L
  for(b in seq_len(B)){
    sim <- sim_draw(n)
    stat <- ks_stat_discrete(sim, cdf_fn); if(stat >= obs) geq <- geq + 1L
  }
  (geq + 1)/(B + 1)
}

## ---- Scan all three models over f, then pick frontier for dLN ----
scan_one_multimodel <- function(x, f_grid=F_GRID, min_n=MIN_N_KEPT, B=B_BOOT){
  out <- vector("list", length(f_grid))
  for(i in seq_along(f_grid)){
    f <- f_grid[i]; kmin <- kmin_by_fraction(x, f); kept <- x[x>=kmin]
    n_kept <- length(kept)
    if(n_kept < min_n || length(unique(kept)) < 2L){
      out[[i]] <- data.frame(f=f, kmin=kmin, n_kept=n_kept,
                             p_dLN=NA_real_, p_Zipf=NA_real_, p_Geom=NA_real_)
      next
    }
    # dLN
    dfit <- fit_dln(kept, kmin)
    cdf_d <- function(t) dln_cdf(sort(unique(t)), dfit$mu, dfit$sigma, kmin)
    p_dln <- ks_bootstrap_p(kept, cdf_d,
                            sim_draw=function(n) rdln(n, dfit$mu, dfit$sigma, kmin),
                            B=B, seed_offset=kmin+1L)
    # Zipf
    zfit <- fit_zipf_finite(kept, kmin)
    cdf_z <- function(t) zipf_cdf(zfit, sort(unique(t)))
    p_zip <- ks_bootstrap_p(kept, cdf_z,
                            sim_draw=function(n) rzipf_finite(n, zfit),
                            B=B, seed_offset=kmin+2L)
    # Geom
    gfit <- fit_geom(kept, kmin)
    cdf_g <- function(t) geom_cdf(sort(unique(t)), gfit$theta, kmin)
    p_geo <- ks_bootstrap_p(kept, cdf_g,
                            sim_draw=function(n) rgeom_shift(n, gfit$theta, kmin),
                            B=B, seed_offset=kmin+3L)
    
    out[[i]] <- data.frame(f=f, kmin=kmin, n_kept=n_kept,
                           p_dLN=p_dln, p_Zipf=p_zip, p_Geom=p_geo)
  }
  do.call(rbind, out)
}

choose_frontier_dln <- function(scan_df){
  ok <- scan_df[is.finite(scan_df$p_dLN) & scan_df$p_dLN >= 0.05 & scan_df$n_kept >= MIN_N_KEPT, , drop=FALSE]
  if(nrow(ok) > 0) {
    max(ok$f)
  } else {
    avail <- scan_df[scan_df$n_kept >= MIN_N_KEPT, , drop=FALSE]
    if(nrow(avail) == 0) NA_real_ else max(avail$f, na.rm=TRUE)
  }
}
## ---------- /FrontierTools_min.R ----------