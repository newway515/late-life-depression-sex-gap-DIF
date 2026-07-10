# =====================================================================
# 11_missingness.R  —  Reviewer point (Stat #5), self-contained.
#   Rebuilds the ELIGIBLE CHARLS w4 sample (age>=60, valid sex) WITH
#   missing CES-D items retained (mirrors 00_export_from_db.R filters,
#   minus the complete-item requirement), then:
#     (a) describes item missingness;
#     (b) tests predictors of "any item missing" (base-R glm);
#     (c) compares complete-case vs multiple-imputation Delta_latent.
#   Also writes CSVs into ./export/ so the figure/table are done in Python.
#
# Run in phaseB_scripts/.  Packages: DBI, RSQLite, dplyr, mirt (+ mice for MI).
# =====================================================================
suppressMessages({library(DBI); library(RSQLite); library(dplyr); library(mirt)})
DB <- Sys.getenv("CESD_DB", "D:/clinicdatabase/SQLitedatabase/cesd_analysis.db")
dir.create("export", showWarnings = FALSE)
ITEMS <- c("depres","effort","sleepr","whappy","flone","going","bother","mindts","fhope","fear")
ID <- paste0(ITEMS,"_d")

con <- dbConnect(SQLite(), DB); charls <- dbReadTable(con, "charls_cesd_items_long"); dbDisconnect(con)
charls <- charls %>% mutate(across(all_of(c(ID,"agey","ragender","raeducl")), as.numeric))

## ELIGIBLE w4 (items may be missing) -- same filters as 00 minus item-completeness
elig <- charls %>% filter(wave==4, agey>=60, ragender %in% c(1,2)) %>%
  mutate(female = as.integer(ragender==2))
X <- elig[, ID]
nmiss <- rowSums(is.na(X))
elig$n_missing <- nmiss
elig$any_missing <- as.integer(nmiss > 0)
elig$obs_mean <- rowMeans(X, na.rm=TRUE)          # NaN if all 10 missing

n <- nrow(elig)
n_cc <- sum(nmiss==0 & !is.na(elig$raeducl))
cat(sprintf("Eligible (w4, age>=60, valid sex) N = %d\n", n))
cat(sprintf("Complete CES-D + education (= analytic sample) N = %d  (should match charls_w4.csv)\n", n_cc))
cat(sprintf("Any CES-D item missing: %d (%.1f%%);  education missing: %d (%.1f%%)\n",
            sum(nmiss>0), 100*mean(nmiss>0), sum(is.na(elig$raeducl)), 100*mean(is.na(elig$raeducl))))
per_item <- round(100*colMeans(is.na(X)),2); names(per_item)<-ITEMS
cat("Per-item missing %:\n"); print(per_item)
cat(sprintf("Mean # items missing = %.3f (max %d)\n", mean(nmiss), max(nmiss)))

## (b) predictors of missingness (base R, no extra pkg)
dd <- elig %>% filter(!is.na(raeducl), !is.nan(obs_mean))
fit <- glm(any_missing ~ female + raeducl + obs_mean, data=dd, family=binomial)
cat("\n== glm: P(any item missing) ~ female + education + observed-symptom mean ==\n")
print(round(summary(fit)$coef,4))
cat("OR [95% CI]:\n"); print(round(exp(cbind(OR=coef(fit), confint.default(fit))),3))

## export rows for Python figure/table
write.csv(elig[,c("female","raeducl","agey","n_missing","obs_mean","any_missing")],
          "export/export_missingness_rows.csv", row.names=FALSE)

## (c) complete-case vs MICE Delta_latent (needs mirt + mice)
prev <- readRDS("output/out_charls_w4_dif_gap.rds"); anchor <- prev$anchor
dl_from <- function(Xi, fem){
  Xi <- as.data.frame(lapply(Xi, function(z) as.integer(round(z)))); colnames(Xi)<-ITEMS
  g <- factor(ifelse(fem==1,"F","M"), levels=c("M","F"))
  m <- multipleGroup(Xi,1,group=g,itemtype="graded",
        invariance=c(anchor,"free_means","free_var"), verbose=FALSE)
  coef(m,simplify=TRUE)$F$means[1]
}
cc <- elig[nmiss==0 & !is.na(elig$raeducl), ]
dl_cc <- dl_from(cc[,ID], cc$female)
if(requireNamespace("mice", quietly=TRUE)){
  library(mice)
  base <- elig %>% filter(!is.na(raeducl))          # impute item NAs; keep rows with known education
  imp_df <- cbind(base[,ID], female=base$female, raeducl=base$raeducl, agey=base$agey)
  set.seed(20260709)
  imp <- mice(imp_df, m=5, method="pmm", printFlag=FALSE)
  dls <- sapply(1:5, function(k){ ci<-complete(imp,k); dl_from(ci[,ID], ci$female) })
  out <- data.frame(method=c("complete_case", paste0("MI_",1:5)), delta_latent=c(dl_cc, dls))
  cat(sprintf("\n== Delta_latent: complete-case = %.3f | MICE(m=5) mean = %.3f [range %.3f, %.3f] ==\n",
              dl_cc, mean(dls), min(dls), max(dls)))
} else {
  out <- data.frame(method="complete_case", delta_latent=dl_cc)
  cat("\n(mice not installed -> only complete-case Delta_latent; install.packages('mice') to add MI.)\n")
}
write.csv(out, "export/export_missingness_delta.csv", row.names=FALSE)
cat("\nWrote export/export_missingness_rows.csv and export/export_missingness_delta.csv\n")
cat("Send me those two CSVs; I will make Table S19 + the missingness figure in Python.\n")
