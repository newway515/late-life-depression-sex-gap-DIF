# =====================================================================
# 09_missingness.R  —  Reviewer point (Stat #5):
#   (a) Describe item-missingness in the ELIGIBLE CHARLS w4 sample
#       (age >= 60, valid sex) BEFORE complete-case filtering.
#   (b) Test whether "any item missing" is related to sex / education /
#       symptom level (logistic regression).
#   (c) Compare the complete-case Δ_latent with a multiple-imputation
#       (MICE) estimate as a sensitivity analysis.
#
# INPUT: ./data/charls_w4_raw.csv  = the eligible sample WITH missing
#   item values retained. Columns needed:
#     depres_d ... fear_d  (10 items, 0-3, NA allowed), female, raeducl, agey
#
#   If you do not yet have this file, create it once from the analysis DB,
#   e.g. (edit table/column names to match your DB):
#     library(DBI); library(RSQLite)
#     con <- dbConnect(SQLite(), "D:/clinicdatabase/SQLitedatabase/cesd_analysis.db")
#     raw <- dbGetQuery(con, "SELECT * FROM charls_cesd_items_long WHERE wave=4")
#     # keep age>=60 & valid sex; map item columns to depres_d..fear_d (0-3, NA kept)
#     write.csv(raw_mapped, "phaseB_scripts/data/charls_w4_raw.csv", row.names=FALSE)
# =====================================================================
suppressMessages({library(dplyr)})
ITEMS <- c("depres","effort","sleepr","whappy","flone","going","bother","mindts","fhope","fear")

fp <- "./data/charls_w4_raw.csv"
if(!file.exists(fp)){
  cat("!! ./data/charls_w4_raw.csv not found.\n",
      "   This script needs the eligible (age>=60, valid sex) sample WITH missing items retained.\n",
      "   See the header of this file for a one-off export template from your DB.\n",
      "   (The complete-case CSVs in ./data/ cannot show the missingness pattern.)\n")
  quit(save="no")
}
d <- read.csv(fp)
X <- d[, paste0(ITEMS,"_d")]
n <- nrow(d)

## (a) description
per_item <- round(100*colMeans(is.na(X)), 2); names(per_item) <- ITEMS
nmiss    <- rowSums(is.na(X))
cat(sprintf("Eligible N = %d\n", n))
cat(sprintf("Any CES-D item missing: %d (%.1f%%); complete cases: %d (%.1f%%)\n",
            sum(nmiss>0), 100*mean(nmiss>0), sum(nmiss==0), 100*mean(nmiss==0)))
cat("Per-item missing %:\n"); print(per_item)
cat(sprintf("Mean # items missing per respondent = %.3f (max %d)\n", mean(nmiss), max(nmiss)))

## (b) predictors of missingness
d$any_missing <- as.integer(nmiss > 0)
d$obs_mean    <- rowMeans(X, na.rm=TRUE)          # severity proxy from observed items
fit <- glm(any_missing ~ female + raeducl + obs_mean, data=d, family=binomial)
cat("\n== Logistic regression: P(any item missing) ~ female + education + observed-symptom mean ==\n")
print(round(summary(fit)$coef, 4))
cat("OR [95% CI]:\n"); print(round(exp(cbind(OR=coef(fit), confint.default(fit))), 3))

## (c) MICE sensitivity vs complete-case  (optional; needs 'mice' + 'mirt')
ok <- requireNamespace("mice", quietly=TRUE) && requireNamespace("mirt", quietly=TRUE)
if(ok){
  library(mice); library(mirt)
  prev <- readRDS("output/out_charls_w4_dif_gap.rds"); anchor <- prev$anchor
  dl_from <- function(Xi, fem){
    Xi <- as.data.frame(lapply(Xi, function(z) as.integer(round(z)))); colnames(Xi) <- ITEMS
    g <- factor(ifelse(fem==1,"F","M"), levels=c("M","F"))
    m <- multipleGroup(Xi,1,group=g,itemtype="graded",
          invariance=c(anchor,"free_means","free_var"), verbose=FALSE)
    coef(m, simplify=TRUE)$F$means[1]
  }
  # complete-case Δ_latent
  cc <- d[nmiss==0, ]; dl_cc <- dl_from(cc[,paste0(ITEMS,"_d")], cc$female)
  # MICE (m=5) on the item matrix + female/educ/age as predictors
  imp_df <- cbind(X, female=d$female, raeducl=d$raeducl, agey=d$agey)
  set.seed(20260709)
  imp <- mice(imp_df, m=5, method="pmm", printFlag=FALSE)
  dls <- sapply(1:5, function(k){ ci <- complete(imp,k); dl_from(ci[,paste0(ITEMS,"_d")], ci$female) })
  cat(sprintf("\n== Δ_latent: complete-case = %.3f | MICE (m=5) mean = %.3f [range %.3f, %.3f] ==\n",
              dl_cc, mean(dls), min(dls), max(dls)))
  saveRDS(list(per_item=per_item, fit=fit, dl_cc=dl_cc, dl_mi=dls), "output/out_missingness.rds")
} else {
  cat("\n(mice/mirt not available -> skipped the CC-vs-MI Δ_latent comparison; install.packages('mice') to enable.)\n")
  saveRDS(list(per_item=per_item, fit=fit), "output/out_missingness.rds")
}
cat("Saved output/out_missingness.rds\n")
