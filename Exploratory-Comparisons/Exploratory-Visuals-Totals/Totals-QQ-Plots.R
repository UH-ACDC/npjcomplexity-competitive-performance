############################################################
## Figure 1: Cross-domain QQ comparisons (system totals)
## Layout: QQ in lower-left, r in upper-right, n on diagonal
############################################################

## Base graphics packages
library(stats)      # cor, lm, density, qqplot
library(graphics)   # hist, plot, rect, text, box
library(grDevices)  # pdf, dev.off

## Set working directory to script location when in RStudio
if (interactive() &&
    requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

############################################################
## 1. Build system-total vectors from the breakdown CSVs
############################################################

## 1A. Luftwaffe pilots: entry cohorts 1940–1942 (n = 1148)
pilots <- read.csv("GER_pilots.csv", stringsAsFactors = FALSE)
year_cols_pilots <- grep("^y[0-9]+", names(pilots), value = TRUE)

pilots$total_victories <- rowSums(pilots[, year_cols_pilots], na.rm = TRUE)

pilots$first_year <- apply(pilots[, year_cols_pilots], 1, function(z) {
  yrs <- as.integer(sub("y", "", year_cols_pilots)[z > 0])
  if (length(yrs) == 0L) NA_integer_ else min(yrs)
})

WW2_GER_df <- subset(pilots, first_year %in% 1940:1942)
WW2_GER    <- WW2_GER_df$total_victories   # n = 1148 expected


## 1B. NSF CS grantees: first award <= 2010, score >= 2 (n = 1323)
cs <- read.csv("SCHOLARS_CS_breakdown.csv", stringsAsFactors = FALSE)
year_cols_cs <- grep("^y[0-9]+", names(cs), value = TRUE)

cs$first_year <- apply(cs[, year_cols_cs], 1, function(z) {
  yrs <- as.integer(sub("y", "", year_cols_cs)[z > 0])
  if (length(yrs) == 0L) NA_integer_ else min(yrs)
})

NSF_CS_df <- subset(cs,
                    !is.na(first_year) &
                      first_year <= 2010 &
                      score >= 2)

NSF_CS <- NSF_CS_df$score    # n = 1323


## 1C. NIH BIO grantees: first award <= 2010, score >= 2 (n = 524)
bio <- read.csv("SCHOLARS_BIO_breakdown.csv", stringsAsFactors = FALSE)
year_cols_bio <- grep("^y[0-9]+", names(bio), value = TRUE)

bio$first_year <- apply(bio[, year_cols_bio], 1, function(z) {
  yrs <- as.integer(sub("y", "", year_cols_bio)[z > 0])
  if (length(yrs) == 0L) NA_integer_ else min(yrs)
})

NIH_BIO_df <- subset(bio,
                     !is.na(first_year) &
                       first_year <= 2010 &
                       score >= 2)

NIH_BIO <- NIH_BIO_df$score  # n = 524


## 1D. U.S. Olympic swimmers: total 3/2/1 score >= 2 (n = 446)
swim_us <- read.csv("ATHLETES_USA_SWIMMERS_breakdown.csv",
                    stringsAsFactors = FALSE)

games_cols_us <- grep("^X[0-9]{4}", names(swim_us), value = TRUE)
swim_us$total_score <- rowSums(swim_us[, games_cols_us], na.rm = TRUE)

OMS_USA_SW_df <- subset(swim_us, total_score >= 2)
OMS_USA_SW    <- OMS_USA_SW_df$total_score   # n = 446


## 1E. French Olympic fencers: total 3/2/1 score >= 2 (n = 137)
fence_fr <- read.csv("ATHLETES_FRA_FENCING_breakdown.csv",
                     stringsAsFactors = FALSE)

games_cols_fr <- grep("^X[0-9]{4}", names(fence_fr), value = TRUE)
fence_fr$total_score <- rowSums(fence_fr[, games_cols_fr], na.rm = TRUE)

OMS_FRA_F_df <- subset(fence_fr, total_score >= 2)
OMS_FRA_F    <- OMS_FRA_F_df$total_score     # n = 137


## Bundle into a list and names for plotting
Data <- list(
  WW2_GER,     # LW
  NSF_CS,      # CS
  NIH_BIO,     # BIO
  OMS_USA_SW,  # USOT_S
  OMS_FRA_F    # FRAOT_F
)

VarNames <- c("LW", "CS", "BIO", "USOT_S", "FRAOT_F")
n <- length(VarNames)

############################################################
## 2. PDF device
############################################################

pdf("Figure_1_QQ_from_breakdown.pdf",
    width  = 8.5,
    height = 7.0)

############################################################
## 3. Plot matrix
##   - mfrow: row-wise filling so (i,j) is row/column
##   - Lower-left (j < i): QQ plots
##   - Upper-right (j > i): r values in boxes
##   - Diagonal: hist + density + "n = ..."
############################################################

par(mfrow = c(n + 1, n + 1),
    mar   = c(1, 1, 1, 1),
    oma   = c(1, 1, 1, 1))
par(mgp = c(0.5, 0.2, 0))

get_header_col <- function(idx) {
  if (idx == 1) {
    "lightskyblue1"         # LW
  } else if (idx %in% 2:3) {
    "darkseagreen1"         # CS, BIO
  } else {
    "gold1"                 # USOT_S, FRAOT_F
  }
}

## Top row: blank + column headers
plot(0, 0, type = "n", xaxt = "n", yaxt = "n", ann = FALSE, bty = "n")
for (j in 1:n) {
  plot(0, 0, type = "n", xaxt = "n", yaxt = "n", ann = FALSE, bty = "n")
  rect(par("usr")[1], par("usr")[3],
       par("usr")[2], par("usr")[4],
       col = get_header_col(j), border = NA)
  text(0, 0, VarNames[j], cex = 1.7)
}

## Remaining rows: row header + cells
for (i in 1:n) {
  
  ## Row header (leftmost column)
  plot(0, 0, type = "n", xaxt = "n", yaxt = "n", ann = FALSE, bty = "n")
  rect(par("usr")[1], par("usr")[3],
       par("usr")[2], par("usr")[4],
       col = get_header_col(i), border = NA)
  text(0, 0, VarNames[i], cex = 1.7)
  
  for (j in 1:n) {
    
    if (i == j) {
      ## Diagonal: histogram + density + N + x-axis ticks
      x <- Data[[i]]
      hist(x,
           probability = TRUE,
           main = "",
           xlab = "",
           ylab = "",
           axes = FALSE,
           col  = "grey85",
           border = "grey60")
      
      ## x-axis tick marks and labels
      xticks <- pretty(range(x))
      axis(1, at = xticks, labels = xticks, cex.axis = 0.7)
      
      lines(density(x), lwd = 2)
      box()
      
      usr <- par("usr")
      text(x = usr[2] - 0.02 * diff(usr[1:2]),
           y = usr[4] - 0.05 * diff(usr[3:4]),
           labels = bquote(italic(N) == .(length(x))),
           adj = c(1, 1),
           cex = 1.0)
    } else if (j < i) {
      ## LOWER-LEFT triangle: QQ plot
      y <- Data[[i]]
      x <- Data[[j]]
      qq <- qqplot(x, y,
                   xlab = "",
                   ylab = "",
                   axes = FALSE)
      abline(lm(qq$y ~ qq$x), col = "red", lwd = 2)
      box()
      
    } else {  ## j > i
      ## UPPER-RIGHT triangle: r value in boxed panel
      y <- Data[[i]]
      x <- Data[[j]]
      qq <- qqplot(x, y, plot = FALSE)
      r  <- cor(qq$x, qq$y)
      
      plot(0, 0,
           type = "n",
           xaxt = "n",
           yaxt = "n",
           xlab = "",
           ylab = "",
           bty  = "n")
      box()
      text(0, 0,
           labels = paste0("r = ", sprintf("%.3f", r)),
           cex = 1.5)
    }
  }
}

dev.off()
############################################################
## End of script – QQ lower-left, r upper-right, n on diagonal
############################################################