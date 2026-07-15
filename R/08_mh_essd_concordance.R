# =====================================================================
# 08_mh_essd_concordance.R — concordance of MH and IRT-ESSD screens.
#   Agreement between the nonparametric MH (ETS A/B/C) screen and the
#   IRT ESSD effect-size gate. Produces a per-item table + 2x2 crosstab
#   (MH class B/C vs ESSD >= 0.20) + Cohen's kappa + Spearman correlation
#   of the two continuous magnitudes.
# =====================================================================
suppressMessages({library(dplyr)})
ITEMS <- c("depres","effort","sleepr","whappy","flone","going","bother","mindts","fhope","fear")

prev <- readRDS("output/out_charls_w4_dif_gap.rds")
ess  <- prev$essd

# MH ETS delta/class from the exploratory screen file. Try common locations.
paths <- c("./output/charls_dif_screen_w4.csv", "./data/charls_dif_screen_w4.csv",
           "charls_dif_screen_w4.csv", "../charls_dif_screen_w4.csv")
scrp  <- paths[file.exists(paths)][1]
if(is.na(scrp)) stop("Cannot find charls_dif_screen_w4.csv. Run python/charls_analysis_phaseA.py first.")
scr <- read.csv(scrp)

tab <- data.frame(
  item      = ITEMS,
  MH_delta  = scr$sex_ETSdelta[match(ITEMS, scr$item)],
  MH_class  = scr$sex_class[match(ITEMS, scr$item)],
  ESSD      = round(as.numeric(ess[ITEMS]), 3)
)
tab$MH_substantive   <- tab$MH_class %in% c("B","C")
tab$ESSD_substantive <- tab$ESSD >= 0.20
cat("== Per-item MH vs ESSD ==\n"); print(tab, row.names=FALSE)

ct <- table(MH_BorC   = factor(tab$MH_substantive,   c(FALSE,TRUE), c("A","B/C")),
            ESSD_ge020 = factor(tab$ESSD_substantive, c(FALSE,TRUE), c("<0.20",">=0.20")))
cat("\n== 2x2 concordance ==\n"); print(ct)

po <- sum(diag(ct))/sum(ct)
pe <- sum(rowSums(ct)*colSums(ct))/sum(ct)^2
kap <- (po-pe)/(1-pe)
rho <- suppressWarnings(cor(abs(tab$MH_delta), tab$ESSD, method="spearman", use="complete.obs"))
cat(sprintf("\nObserved agreement = %.2f | Cohen's kappa = %.2f | Spearman(|MH Δ|, ESSD) = %.2f\n",
            po, kap, rho))
cat("Note: sleepr is B/C on MH AND ESSD>=0.20; fear & flone reach ESSD>=0.20 but only A on the coarser MH\n",
    "screen at w4 (fear is B in w1-w3), illustrating the greater sensitivity of the IRT effect-size gate.\n")
write.csv(tab, "output/out_mh_essd_concordance.csv", row.names=FALSE)
cat("Saved output/out_mh_essd_concordance.csv\n")
