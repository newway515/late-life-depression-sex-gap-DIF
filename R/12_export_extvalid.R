# =====================================================================
# 12_export_extvalid.R — export inputs for pooled person-clustered models.
#   The external-validity files are person-wave POOLED (age>=60, all waves):
#     HRS-CIDI  169k rows / 31,572 persons; ELSA-grip 59k rows / 14,161 persons.
#   This script ONLY refits the corrected latent score and exports, per row,
#   the person id + wave + sex + age + raw score + corrected score + outcome,
#   so the cluster-robust (clustered on person) regressions can be done in Python.
#   CHARLS grip (single wave w3) needs no clustering and is unchanged.
# Run from the repository root. Packages: mirt, dplyr. Runtime: ~2-4 min.
# =====================================================================
suppressMessages({library(mirt); library(dplyr)})
set.seed(20260709); dir.create("export", showWarnings = FALSE)
EN <- c("depres","effort","sleepr","whappy","flone","fsad","going","enlife")  # ELSA/HRS 8-item
ESSD_CUT <- 0.20; ALPHA_DIF <- 0.01

corrected_fscore <- function(X, female, itemtype="2PL"){
  X <- as.data.frame(lapply(X, function(v) as.integer(round(v))))
  g <- factor(ifelse(female==1,"F","M"), levels=c("M","F"))
  ncat <- max(as.matrix(X), na.rm=TRUE); dpar <- if(itemtype=="graded") paste0("d",1:ncat) else "d"
  full <- multipleGroup(X,1,group=g,itemtype=itemtype,
            invariance=c(colnames(X),"free_means","free_var"), verbose=FALSE)
  dif <- DIF(full, which.par=c("a1",dpar), scheme="drop", items2test=1:ncol(X), p.adjust="BH")
  adjp <- dif$adj_p[match(colnames(X), rownames(dif))]
  mgC <- multipleGroup(X,1,group=g,itemtype=itemtype, verbose=FALSE)
  Th <- matrix(seq(-4,4,length.out=61)); w <- dnorm(Th[,1]); w <- w/sum(w)
  essd <- sapply(seq_len(ncol(X)), function(i){
    eR <- expected.item(extract.item(mgC,i,group="M"),Th,min=0)
    eF <- expected.item(extract.item(mgC,i,group="F"),Th,min=0)
    sdi <- sd(X[[i]],na.rm=TRUE); as.numeric(sum(abs(eF-eR)*w)/ifelse(is.finite(sdi)&&sdi>0,sdi,1)) })
  flagged <- colnames(X)[ which(adjp < ALPHA_DIF & essd >= ESSD_CUT) ]
  anc <- setdiff(colnames(X), flagged); if(length(anc)<2) anc <- colnames(X)
  part <- multipleGroup(X,1,group=g,itemtype=itemtype,
            invariance=c(anc,"free_means","free_var"), verbose=FALSE)
  list(fs=fscores(part,method="EAP")[,1], flagged=flagged)
}

## ---- HRS-CIDI (pooled) ----
h <- read.csv("./data/hrs_with_cidi.csv") %>% filter(!is.na(cidi_mde3))
Xh <- h[,paste0(EN,"_d")]; colnames(Xh)<-EN
cf <- corrected_fscore(Xh, h$female, "2PL")
cat("HRS flagged:", paste(cf$flagged, collapse=","), "| n=", nrow(h), "| persons=", length(unique(h$hhidpn)), "\n")
write.csv(data.frame(id=h$hhidpn, wave=h$wave, female=h$female, agey=h$agey,
                     cesd=h$cesd, corr=cf$fs, cidi_mde3=h$cidi_mde3),
          "export/export_hrs_cidi.csv", row.names=FALSE)

## ---- ELSA grip (pooled) ----
e <- read.csv("./data/elsa_with_anchors.csv") %>% filter(!is.na(grip_max))
Xe <- e[,paste0(EN,"_d")]; colnames(Xe)<-EN
ce <- corrected_fscore(Xe, e$female, "2PL")
cat("ELSA flagged:", paste(ce$flagged, collapse=","), "| n=", nrow(e), "| persons=", length(unique(e$idauniq)), "\n")
write.csv(data.frame(id=e$idauniq, wave=e$wave, female=e$female, agey=e$agey,
                     cesd=e$cesd, corr=ce$fs, grip_max=e$grip_max),
          "export/export_elsa_grip.csv", row.names=FALSE)

cat("\nWrote export/export_hrs_cidi.csv and export/export_elsa_grip.csv.\n")
cat("Run python/clustered_external_validity.py to fit the person-clustered models.\n")
