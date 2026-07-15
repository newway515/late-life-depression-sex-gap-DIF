# =====================================================================
# 06_dd_bootstrap.R — joint bootstrap for the decomposition components.
#   Are ΔΔ_meas (=Δ_raw−Δ_latent) and ΔΔ_struct (=Δ_latent−Δ_adj)
#   individually different from 0, and different from EACH OTHER?
#   -> joint bootstrap of Δ_raw / Δ_latent / Δ_adj with the anchor set
#      HELD FIXED (the same anchors used in the main analysis).
# Output: point estimates + 95% percentile CIs for the three estimands,
#   for ΔΔ_meas, ΔΔ_struct, and for (ΔΔ_meas − ΔΔ_struct) with a p value.
# =====================================================================
suppressMessages({library(mirt); library(boot); library(dplyr)})
set.seed(20260709)
N_BOOT <- 2000
ITEMS  <- c("depres","effort","sleepr","whappy","flone","going","bother","mindts","fhope","fear")

d    <- read.csv("./data/charls_w4.csv")
prev <- readRDS("output/out_charls_w4_dif_gap.rds")
anchor <- prev$anchor                       # FIXED anchor set from the main analysis
stopifnot(length(anchor) >= 2)

cohend <- function(a,b){
  s <- sqrt(((length(a)-1)*var(a)+(length(b)-1)*var(b))/(length(a)+length(b)-2)); (mean(a)-mean(b))/s
}

# returns c(Delta_raw, Delta_latent, Delta_adj, DDmeas, DDstruct) for a data slice
stat <- function(dd){
  Xi <- dd[,paste0(ITEMS,"_d")]; colnames(Xi) <- ITEMS
  Xi[] <- lapply(Xi, function(z) as.integer(round(z)))
  g <- factor(ifelse(dd$female==1,"F","M"), levels=c("M","F"))   # 男=参照
  draw <- cohend(dd$cesd[dd$female==1], dd$cesd[dd$female==0])
  m <- tryCatch(multipleGroup(Xi,1,group=g,itemtype="graded",
        invariance=c(anchor,"free_means","free_var"), verbose=FALSE), error=function(e) NULL)
  if(is.null(m)) return(c(NA,NA,NA,NA,NA))
  dlat <- coef(m, simplify=TRUE)$F$means[1]
  fs   <- fscores(m, method="EAP")[,1]
  dadj <- tryCatch(as.numeric(coef(lm(fs ~ dd$female + dd$raeducl + dd$agey))["dd$female"]),
                   error=function(e) NA_real_)
  c(draw=draw, dlat=dlat, dadj=dadj, ddmeas=draw-dlat, ddstruct=dlat-dadj)
}

obs <- stat(d)
cat(sprintf("Fixed anchors: %s\nRunning %d bootstrap resamples (may take ~15-25 min)...\n",
            paste(anchor, collapse=", "), N_BOOT))
bt <- boot(d, function(dat, idx) stat(dat[idx,]), R=N_BOOT)

pct <- function(x) quantile(x, c(.025,.975), na.rm=TRUE)
labs <- c("Delta_raw","Delta_latent","Delta_adj","DDmeas (raw-latent)","DDstruct (latent-adj)")
cat("\n== Estimand / component: point [95% CI] ==\n")
for(i in 1:5){ ci <- pct(bt$t[,i]); cat(sprintf("  %-22s %.3f [%.3f, %.3f]\n", labs[i], obs[i], ci[1], ci[2])) }

diffcol <- bt$t[,4] - bt$t[,5]                          # ΔΔ_meas − ΔΔ_struct
ciD <- pct(diffcol)
pD  <- 2*min(mean(diffcol<=0, na.rm=TRUE), mean(diffcol>=0, na.rm=TRUE))
excl0 <- function(j){ ci<-pct(bt$t[,j]); !(ci[1]<=0 & ci[2]>=0) }
cat(sprintf("\n  DDmeas 95%% CI excludes 0:  %s", excl0(4)))
cat(sprintf("\n  DDstruct 95%% CI excludes 0: %s", excl0(5)))
cat(sprintf("\n  DDmeas − DDstruct = %.3f [95%% CI %.3f, %.3f], two-sided p = %.3f\n",
            obs[4]-obs[5], ciD[1], ciD[2], pD))

saveRDS(list(obs=obs, t=bt$t, diff=diffcol, anchor=anchor), "output/out_dd_bootstrap.rds")
cat("\nSaved output/out_dd_bootstrap.rds\n")
