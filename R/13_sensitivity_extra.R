# =====================================================================
# 13_sensitivity_extra.R  —  Reviewer-2 (round 2) points #2 and #3.
#   (#2) Contemporary-HRS-wave sensitivity: recompute the cross-cohort
#        gradient using a recent HRS wave (2016 and 2018) instead of the
#        most-complete wave (1998), to show China-largest / all-positive
#        does not depend on HRS wave year.
#   (#3) Bifactor-general-factor sensitivity: re-estimate the CHARLS w4
#        latent sex gap scoring the general factor of a bifactor model,
#        to show Delta_latent stays ~0.32 (the positive-affect method
#        factor does not absorb the sex DIF).
# Run in phaseB_scripts/.  Packages: DBI, RSQLite, dplyr, mirt.  ~3-6 min.
# =====================================================================
suppressMessages({library(DBI); library(RSQLite); library(dplyr); library(mirt)})
set.seed(20260709); dir.create("export", showWarnings = FALSE)
DB <- Sys.getenv("CESD_DB", "D:/clinicdatabase/SQLitedatabase/cesd_analysis.db")
COMMON6 <- c("depres","effort","sleepr","whappy","flone","going")
ITEMS   <- c("depres","effort","sleepr","whappy","flone","going","bother","mindts","fhope","fear")

## ---------- (#2) contemporary HRS wave: within-cohort corrected sex gap ----------
## We recompute the HRS common-6 (binary 2PL) DIF-corrected latent sex gap at
## HRS wave 4 (1998), wave 13 (2016) and wave 14 (2018); compare to CHARLS 0.319.
con <- dbConnect(SQLite(), DB); hrs <- dbReadTable(con,"hrs_cesd_items_long"); dbDisconnect(con)
b6 <- paste0(COMMON6,"_d")
hrs <- hrs %>% mutate(across(all_of(c(b6,"agey","ragender")), as.numeric))
corr_gap_binary <- function(d){
  d <- d %>% filter(agey>=60, ragender %in% c(1,2), if_all(all_of(b6), ~ !is.na(.))) %>%
       mutate(female=as.integer(ragender==2))
  if(nrow(d) < 300) return(c(n=nrow(d), gap=NA))
  X <- sapply(COMMON6, function(b) as.integer(d[[paste0(b,"_d")]] >= 1)); X <- as.data.frame(X); colnames(X)<-COMMON6
  g <- factor(ifelse(d$female==1,"F","M"), levels=c("M","F"))
  full <- multipleGroup(X,1,group=g,itemtype="2PL",invariance=c(COMMON6,"free_means","free_var"),verbose=FALSE)
  dif <- DIF(full, which.par=c("a1","d"), scheme="drop", items2test=1:ncol(X), p.adjust="BH")
  adjp <- dif$adj_p[match(COMMON6, rownames(dif))]
  mgC <- multipleGroup(X,1,group=g,itemtype="2PL",verbose=FALSE)
  Th <- matrix(seq(-4,4,length.out=61)); w<-dnorm(Th[,1]); w<-w/sum(w)
  essd <- sapply(seq_along(COMMON6), function(i){
    eR<-expected.item(extract.item(mgC,i,group="M"),Th,min=0); eF<-expected.item(extract.item(mgC,i,group="F"),Th,min=0)
    sdi<-sd(X[[i]],na.rm=TRUE); as.numeric(sum(abs(eF-eR)*w)/ifelse(is.finite(sdi)&&sdi>0,sdi,1))})
  flg <- COMMON6[adjp<0.01 & essd>=0.20]; anc <- setdiff(COMMON6,flg); if(length(anc)<2) anc<-COMMON6
  part <- multipleGroup(X,1,group=g,itemtype="2PL",invariance=c(anc,"free_means","free_var"),verbose=FALSE)
  c(n=nrow(d), gap=as.numeric(coef(part,simplify=TRUE)$F$means[1]))
}
cat("== (#2) HRS common-6 DIF-corrected latent sex gap by wave (year) ==\n")
res2 <- lapply(c(4,13,14), function(w){ r<-corr_gap_binary(hrs %>% filter(wave==w))
  data.frame(hrs_wave=w, approx_year=c(`4`=1998,`13`=2016,`14`=2018)[as.character(w)],
             n=r["n"], corrected_gap=round(r["gap"],3)) })
res2 <- do.call(rbind,res2); rownames(res2)<-NULL; print(res2)
cat("  (compare: CHARLS w4 single-cohort GRM Delta_latent = 0.319; gradient conclusion = China largest, all positive)\n")
write.csv(res2,"export/export_hrs_wave_sensitivity.csv",row.names=FALSE)

## ---------- (#3) bifactor general-factor scoring: CHARLS w4 sex gap ----------
d <- read.csv("./data/charls_w4.csv")
X <- d[,paste0(ITEMS,"_d")]; colnames(X)<-ITEMS; X[]<-lapply(X,function(z)as.integer(round(z)))
g <- factor(ifelse(d$female==1,"F","M"), levels=c("M","F"))
SOMA <- c("effort","sleepr","going"); sfac <- ifelse(ITEMS %in% SOMA, 1L, 2L)
## multiple-group bifactor with all items constrained except group means/vars:
bf <- tryCatch(bfactor(X, sfac, group=g, itemtype="graded",
        invariance=c(ITEMS,"free_means","free_var"),
        technical=list(NCYCLES=3000), verbose=FALSE), error=function(e) NULL)
if(!is.null(bf)){
  gap_bf <- as.numeric(coef(bf, simplify=TRUE)$F$means[1])   # general-factor mean diff F-M
  cat(sprintf("\n== (#3) CHARLS w4 sex gap, bifactor GENERAL factor = %.3f  (vs single-factor Delta_latent 0.319) ==\n", gap_bf))
  write.csv(data.frame(model=c("single_factor_partial_inv","bifactor_general_factor"),
                       gap=c(0.319, round(gap_bf,3))), "export/export_bifactor_gap.csv", row.names=FALSE)
} else cat("\n(#3) bifactor multiple-group did not converge; try technical NCYCLES higher or report single-factor only.\n")

cat("\nWrote export/export_hrs_wave_sensitivity.csv and export/export_bifactor_gap.csv. Send me both.\n")
