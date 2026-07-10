# =====================================================================
# 07_essd_threshold_sensitivity.R  —  Reviewer point (Stat #1):
#   Is Δ_latent sensitive to the ESSD effect-size threshold used to
#   define "substantive" DIF? Re-derive the flagged set and Δ_latent at
#   ESSD >= 0.15 / 0.20 (main) / 0.25, holding adj p < 0.01.
#   Also print the full per-item ESSD so coherence is visible.
# Adds an extra row-family to the sensitivity matrix (§3.8 / Table S15).
# =====================================================================
suppressMessages({library(mirt); library(dplyr)})
ITEMS <- c("depres","effort","sleepr","whappy","flone","going","bother","mindts","fhope","fear")

d    <- read.csv("./data/charls_w4.csv")
prev <- readRDS("output/out_charls_w4_dif_gap.rds")
ess  <- prev$essd                                   # named ESSD per item
adjp <- prev$dif$adj_p[match(ITEMS, prev$dif$item)]; names(adjp) <- ITEMS

cat("== Per-item ESSD (sorted) ==\n")
print(round(sort(ess[ITEMS], decreasing=TRUE), 3))

est_delta <- function(substantive){
  anc <- setdiff(ITEMS, substantive); if(length(anc) < 2) anc <- ITEMS
  X <- d[,paste0(ITEMS,"_d")]; colnames(X) <- ITEMS
  X[] <- lapply(X, function(z) as.integer(round(z)))
  g <- factor(ifelse(d$female==1,"F","M"), levels=c("M","F"))
  m <- multipleGroup(X,1,group=g,itemtype="graded",
        invariance=c(anc,"free_means","free_var"), verbose=FALSE)
  coef(m, simplify=TRUE)$F$means[1]
}

cat("\n== Δ_latent vs. ESSD threshold (adj p < 0.01 required) ==\n")
res <- lapply(c(0.15,0.20,0.25), function(cut){
  sub <- ITEMS[ adjp < 0.01 & ess[ITEMS] >= cut ]
  data.frame(ESSD_cut=cut, n_flagged=length(sub),
             flagged=paste(sub, collapse=","), delta_latent=round(est_delta(sub),3))
})
res <- do.call(rbind, res); print(res, row.names=FALSE)
rng <- diff(range(res$delta_latent))
cat(sprintf("\nΔ_latent range across thresholds = %.3f (robust if <= ~0.05)\n", rng))
write.csv(res, "output/out_essd_threshold_sensitivity.csv", row.names=FALSE)
cat("Saved output/out_essd_threshold_sensitivity.csv\n")
