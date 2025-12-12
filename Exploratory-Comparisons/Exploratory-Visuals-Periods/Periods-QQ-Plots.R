############################################################
## Generational / era-bounded QQ matrix (17×17)
## Layout: QQ lower-left, r upper-right, N on diagonal
## Output: Figure-Generational-QQ-from-breakdown.pdf
############################################################

## Base graphics / I/O
library(stats)      # cor, lm, density, qqplot
library(graphics)   # hist, plot, rect, text, box, axis
library(grDevices)  # pdf, dev.off

## Optional: set working directory to script location in RStudio
if (interactive() &&
    requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

############################################################
## 1. Build the 17 domain–period vectors
############################################################

## ---------------------------------------------------------
## 1A. Luftwaffe pilots – entry cohorts 1940, 1941, 1942
## ---------------------------------------------------------
pilots <- read.csv("GER_pilots.csv", stringsAsFactors = FALSE)
year_cols_pilots <- grep("^y[0-9]+", names(pilots), value = TRUE)

## Total victories per pilot
pilots$total_victories <- rowSums(pilots[, year_cols_pilots], na.rm = TRUE)

## First year with >0 victories
pilots$first_year <- apply(pilots[, year_cols_pilots], 1, function(z) {
  yrs <- as.integer(sub("y", "", year_cols_pilots)[z > 0])
  if (length(yrs) == 0L) NA_integer_ else min(yrs)
})

## Restrict to aces (>=6) with first victory 1940–1942
LW_df <- subset(pilots,
                !is.na(first_year) &
                  first_year %in% 1940:1942 &
                  total_victories >= 6)

LW_1940 <- subset(LW_df, first_year == 1940)$total_victories
LW_1941 <- subset(LW_df, first_year == 1941)$total_victories
LW_1942 <- subset(LW_df, first_year == 1942)$total_victories

## ---------------------------------------------------------
## 1B. NSF CS grantees – first award cohorts
## ---------------------------------------------------------
cs <- read.csv("SCHOLARS_CS_breakdown.csv", stringsAsFactors = FALSE)
year_cols_cs <- grep("^y[0-9]+", names(cs), value = TRUE)

cs$first_year <- apply(cs[, year_cols_cs], 1, function(z) {
  yrs <- as.integer(sub("y", "", year_cols_cs)[z > 0])
  if (length(yrs) == 0L) NA_integer_ else min(yrs)
})

CS_df <- subset(cs,
                !is.na(first_year) &
                  first_year <= 2010 &
                  score >= 2)

CS_96_00 <- subset(CS_df,
                   first_year >= 1996 & first_year <= 2000)$score
CS_01_05 <- subset(CS_df,
                   first_year >= 2001 & first_year <= 2005)$score
CS_06_10 <- subset(CS_df,
                   first_year >= 2006 & first_year <= 2010)$score

## ---------------------------------------------------------
## 1C. NIH BIO grantees – first award cohorts
## ---------------------------------------------------------
bio <- read.csv("SCHOLARS_BIO_breakdown.csv", stringsAsFactors = FALSE)
year_cols_bio <- grep("^y[0-9]+", names(bio), value = TRUE)

bio$first_year <- apply(bio[, year_cols_bio], 1, function(z) {
  yrs <- as.integer(sub("y", "", year_cols_bio)[z > 0])
  if (length(yrs) == 0L) NA_integer_ else min(yrs)
})

BIO_df <- subset(bio,
                 !is.na(first_year) &
                   first_year <= 2010 &
                   score >= 2)

BIO_96_00 <- subset(BIO_df,
                    first_year >= 1996 & first_year <= 2000)$score
BIO_01_05 <- subset(BIO_df,
                    first_year >= 2001 & first_year <= 2005)$score
BIO_06_10 <- subset(BIO_df,
                    first_year >= 2006 & first_year <= 2010)$score

## ---------------------------------------------------------
## 1D. U.S. Olympic swimmers – era windows
## ---------------------------------------------------------
swim_us <- read.csv("ATHLETES_USA_SWIMMERS_breakdown.csv",
                    stringsAsFactors = FALSE)

games_cols_us <- grep("^X[0-9]{4}", names(swim_us), value = TRUE)
years_us      <- as.integer(sub("^X(\\d+).*", "\\1", games_cols_us))

cols_1896_1912 <- games_cols_us[years_us >= 1896 & years_us <= 1912]
cols_1920_1936 <- games_cols_us[years_us >= 1920 & years_us <= 1936]
cols_1948_1976 <- games_cols_us[years_us >= 1948 & years_us <= 1976]
cols_1980_2016 <- games_cols_us[years_us >= 1980 & years_us <= 2016]

USA_1896_1912_score <- rowSums(swim_us[, cols_1896_1912, drop = FALSE],
                               na.rm = TRUE)
USA_1920_1936_score <- rowSums(swim_us[, cols_1920_1936, drop = FALSE],
                               na.rm = TRUE)
USA_1948_1976_score <- rowSums(swim_us[, cols_1948_1976, drop = FALSE],
                               na.rm = TRUE)
USA_1980_2016_score <- rowSums(swim_us[, cols_1980_2016, drop = FALSE],
                               na.rm = TRUE)

USA_1896_1912 <- USA_1896_1912_score[USA_1896_1912_score >= 2]
USA_1920_1936 <- USA_1920_1936_score[USA_1920_1936_score >= 2]
USA_1948_1976 <- USA_1948_1976_score[USA_1948_1976_score >= 2]
USA_1980_2016 <- USA_1980_2016_score[USA_1980_2016_score >= 2]

## ---------------------------------------------------------
## 1E. French Olympic fencers – era windows
## ---------------------------------------------------------
fence_fr <- read.csv("ATHLETES_FRA_FENCING_breakdown.csv",
                     stringsAsFactors = FALSE)

games_cols_fr <- grep("^X[0-9]{4}", names(fence_fr), value = TRUE)
years_fr      <- as.integer(sub("^X(\\d+).*", "\\1", games_cols_fr))

colsF_1896_1912 <- games_cols_fr[years_fr >= 1896 & years_fr <= 1912]
colsF_1920_1936 <- games_cols_fr[years_fr >= 1920 & years_fr <= 1936]
colsF_1948_1976 <- games_cols_fr[years_fr >= 1948 & years_fr <= 1976]
colsF_1980_2016 <- games_cols_fr[years_fr >= 1980 & years_fr <= 2016]

FRA_1896_1912_score <- rowSums(fence_fr[, colsF_1896_1912, drop = FALSE],
                               na.rm = TRUE)
FRA_1920_1936_score <- rowSums(fence_fr[, colsF_1920_1936, drop = FALSE],
                               na.rm = TRUE)
FRA_1948_1976_score <- rowSums(fence_fr[, colsF_1948_1976, drop = FALSE],
                               na.rm = TRUE)
FRA_1980_2016_score <- rowSums(fence_fr[, colsF_1980_2016, drop = FALSE],
                               na.rm = TRUE)

FRA_1896_1912 <- FRA_1896_1912_score[FRA_1896_1912_score >= 2]
FRA_1920_1936 <- FRA_1920_1936_score[FRA_1920_1936_score >= 2]
FRA_1948_1976 <- FRA_1948_1976_score[FRA_1948_1976_score >= 2]
FRA_1980_2016 <- FRA_1980_2016_score[FRA_1980_2016_score >= 2]

## ---------------------------------------------------------
## 1F. Bundle into list in desired order (17 groups)
## ---------------------------------------------------------
Data <- list(
  LW_1940,
  LW_1941,
  LW_1942,
  CS_96_00,
  CS_01_05,
  CS_06_10,
  BIO_96_00,
  BIO_01_05,
  BIO_06_10,
  USA_1896_1912,
  USA_1920_1936,
  USA_1948_1976,
  USA_1980_2016,
  FRA_1896_1912,
  FRA_1920_1936,
  FRA_1948_1976,
  FRA_1980_2016
)

VarNames <- c(
  "LW 1940", "LW 1941", "LW 1942",
  "CS 96–00", "CS 01–05", "CS 06–10",
  "BIO 96–00", "BIO 01–05", "BIO 06–10",
  "US 1896–1912", "US 1920–1936",
  "US 1948–1976", "US 1980–2016",
  "FRA 1896–1912", "FRA 1920–1936",
  "FRA 1948–1976", "FRA 1980–2016"
)

## Plot labels with line breaks for long names
HeaderLabels <- c(
  "LW 1940", "LW 1941", "LW 1942",
  "CS\n1996–2000", "CS\n2001–2005", "CS\n2006–2010",
  "BIO\n1996–2000", "BIO\n2001–2005", "BIO\n2006–2010",
  "US\n1896–1912", "US\n1920–1936",
  "US\n1948–1976", "US\n1980–2016",
  "FRA\n1896–1912", "FRA\n1920–1936",
  "FRA\n1948–1976", "FRA\n1980–2016"
)

n <- length(VarNames)  # 17

############################################################
## 2. PDF device
############################################################

pdf("Figure-Generational-QQ-from-breakdown.pdf",
    width  = 11,
    height = 11)

############################################################
## 3. Plot matrix (same logic as Fig-1 totals)
############################################################

par(mfrow = c(n + 1, n + 1),
    mar   = c(0.1, 0.1, 0.1, 0.1),
    oma   = c(0.5, 0.5, 0.5, 0.5))
par(mgp = c(0.5, 0.2, 0))

## Header colors by block:
##  1–3  = LW (blue)
##  4–6  = CS (light green)
##  7–9  = BIO (darker green)
## 10–13 = US swimmers (light gold)
## 14–17 = FRA fencers (darker gold)
get_header_col <- function(idx) {
  if (idx %in% 1:3) {
    "lightskyblue1"
  } else if (idx %in% 4:6) {
    "darkseagreen1"
  } else if (idx %in% 7:9) {
    "palegreen2"
  } else if (idx %in% 10:13) {
    "gold1"
  } else {
    "goldenrod1"
  }
}

## Top-left corner (blank)
plot(0, 0, type = "n", xaxt = "n", yaxt = "n", ann = FALSE, bty = "n")

## Top row: column headers
for (j in 1:n) {
  plot(0, 0, type = "n", xaxt = "n", yaxt = "n", ann = FALSE, bty = "n")
  usr <- par("usr")
  rect(usr[1], usr[3], usr[2], usr[4],
       col = get_header_col(j), border = NA)
  text(0, 0, HeaderLabels[j], cex = 1.0)
}

## Remaining rows: row header + n×n cells
for (i in 1:n) {
  
  ## Row header at left
  plot(0, 0, type = "n", xaxt = "n", yaxt = "n", ann = FALSE, bty = "n")
  usr <- par("usr")
  rect(usr[1], usr[3], usr[2], usr[4],
       col = get_header_col(i), border = NA)
  text(0, 0, HeaderLabels[i], cex = 1.0)
  
  for (j in 1:n) {
    
    if (i == j) {
      ## Diagonal: histogram + density + N (italic) with x-axis ticks
      x <- Data[[i]]
      hist(x,
           probability = TRUE,
           main = "",
           xlab = "",
           ylab = "",
           axes = FALSE,
           col  = "grey85",
           border = "grey60")
      xticks <- pretty(range(x))
      axis(1, at = xticks, labels = xticks, cex.axis = 0.6)
      lines(density(x), lwd = 1.2)
      box()
      
      usr_inner <- par("usr")
      text(x = usr_inner[2] - 0.02 * diff(usr_inner[1:2]),
           y = usr_inner[4] - 0.05 * diff(usr_inner[3:4]),
           labels = bquote(italic(N) == .(length(x))),
           adj = c(1, 1),
           cex = 0.85)
      
    } else if (j < i) {
      ## LOWER-LEFT triangle: QQ plot with red regression line
      y <- Data[[i]]
      x <- Data[[j]]
      qq <- qqplot(x, y,
                   xlab = "",
                   ylab = "",
                   axes = FALSE)
      abline(lm(qq$y ~ qq$x), col = "red", lwd = 1.2)
      box()
      
    } else {  ## j > i
      ## UPPER-RIGHT triangle: correlation r in a single rectangle
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
      usr_inner <- par("usr")
      
      bg_col <- if (r < 0.80) "khaki1" else "white"
      rect(usr_inner[1], usr_inner[3], usr_inner[2], usr_inner[4],
           col = bg_col, border = "grey60")
      
      text(0, 0,
           labels = sprintf("%.3f", r),
           cex = 0.95)
    }
  }
}

dev.off()
############################################################
## End of script
############################################################