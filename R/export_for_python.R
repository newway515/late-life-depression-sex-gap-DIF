# =====================================================================
# export_for_python.R — fit IRT models and write intermediates to ./export/.
# Run from the repository root so ./data and ./output resolve.
# Packages: mirt, boot, dplyr.  Runtime: ~15-25 min (the bootstrap block).
#
# Produces in ./export/ :
#   export_item_dif.csv          (Stat #2: per-item X2/adj_p/ESSD + MH class)
#   export_essd_thresholds.csv   (Stat #1: Delta_latent at ESSD 0.15/0.20/0.25)
#   export_scores_elsa.csv       (Clinical #2: corrected+raw score+covariates, grip rows)
#   export_scores_charls_grip.csv(Clinical #2: same, CHARLS best grip wave)
#   export_dd_draws.csv          (Major #4: 2000x5 bootstrap draws, fixed anchors)
#   export_missingness_delta.csv (Stat #5: CC vs MI Delta_latent) -- IF raw file present
# =====================================================================
suppressMessages({library(mirt); library(boot); library(dplyr)})
set.seed(20260709)
dir.create("export", showWarnings = FALSE)
ok <- function(tag) cat(sprintf("[done] %s\n", tag))

ITEMS <- c("depres","effort","sleepr","whappy","flone","going","bother","mindts","fhope","fear")
EN    <- c("depres","effort","sleepr","whappy","flone","fsad","going","enlife")   # ELSA/HRS 8-item
CH    <- ITEMS
ESSD_CUT <- 0.20; ALPHA_DIF <- 0.01
z <- function(x) as.numeric(scale(x))

d    <- read.csv("./data/charls_w4.csv")
prev <- readRDS("output/out_charls_w4_dif_gap.rds"); anchor <- prev$anchor

corrected_fscore <- function(X, female, itemtype){
  X <- as.data.frame(lapply(X, function(v) as.integer(round(v))))
  g <- factor(ifelse(female==1,"F","M"), levels=c("M","F"))
  ncat <- max(as.matrix(X), na.rm=TRUE)
  dpar <- if(itemtype=="graded") paste0("d",1:ncat) else "d"
  full <- multipleGroup(X,1,group=g,itemtype=itemtype,
            invariance=c(colnames(X),"free_means","free_var"), verbose=FALSE)
  dif <- DIF(full, which.par=c("a1",dpar), scheme="drop", items2test=1:ncol(X), p.adjust="BH")
  adjp <- dif$adj_p[match(colnames(X), rownames(dif))]
  mgC <- multipleGroup(X,1,group=g,itemtype=itemtype, verbose=FALSE)
  Th <- matrix(seq(-4,4,length.out=61)); w <- dnorm(Th[,1]); w <- w/sum(w)
  essd <- sapply(seq_len(ncol(X)), function(i){
    eR <- expected.item(extract.item(mgC,i,group="M"),Th,min=0)
    eF <- expected.item(extract.item(mgC,i,group="F"),Th,min=0)
    sdi <- sd(X[[i]],na.rm=TRUE); as.numeric(sum(abs(eF-eR)*w)/ifelse(is.finite(sdi)&&sdi>0,sdi,1))
  })
  flagged <- colnames(X)[ which(adjp < ALPHA_DIF & essd >= ESSD_CUT) ]
  anc <- setdiff(colnames(X), flagged); if(length(anc) < 2) anc <- colnames(X)
  part <- multipleGroup(X,1,group=g,itemtype=itemtype,
            invariance=c(anc,"free_means","free_var"), verbose=FALSE)
  fscores(part, method="EAP")[,1]
}

## ---- C) per-item DIF + ESSD + MH class (fast) ----
tryCatch({
  mhp <- c("./output/charls_dif_screen_w4.csv","./data/charls_dif_screen_w4.csv",
           "charls_dif_screen_w4.csv","../charls_dif_screen_w4.csv")
  mhp <- mhp[file.exists(mhp)][1]
  it <- data.frame(item=ITEMS,
                   X2=prev$dif$X2[match(ITEMS,prev$dif$item)],
                   adj_p=prev$dif$adj_p[match(ITEMS,prev$dif$item)],
                   ESSD=as.numeric(prev$essd[ITEMS]))
  if(!is.na(mhp)){ mh<-read.csv(mhp)
    it$MH_delta<-mh$sex_ETSdelta[match(ITEMS,mh$item)]; it$MH_class<-mh$sex_class[match(ITEMS,mh$item)] }
  write.csv(it, "export/export_item_dif.csv", row.names=FALSE); ok("export_item_dif.csv")
}, error=function(e) cat("[skip] item_dif:", conditionMessage(e), "\n"))

## ---- B) ESSD-threshold Delta_latent (fast: 3 refits) ----
tryCatch({
  adjp <- prev$dif$adj_p[match(ITEMS,prev$dif$item)]; names(adjp)<-ITEMS; ess<-prev$essd
  est_delta <- function(sub){
    anc<-setdiff(ITEMS,sub); if(length(anc)<2) anc<-ITEMS
    X<-d[,paste0(ITEMS,"_d")]; colnames(X)<-ITEMS; X[]<-lapply(X,function(z)as.integer(round(z)))
    g<-factor(ifelse(d$female==1,"F","M"),levels=c("M","F"))
    m<-multipleGroup(X,1,group=g,itemtype="graded",invariance=c(anc,"free_means","free_var"),verbose=FALSE)
    coef(m,simplify=TRUE)$F$means[1] }
  thr<-do.call(rbind, lapply(c(0.15,0.20,0.25), function(cut){
    sub<-ITEMS[adjp<0.01 & ess[ITEMS]>=cut]
    data.frame(ESSD_cut=cut, n_flagged=length(sub), flagged=paste(sub,collapse=","),
               delta_latent=round(est_delta(sub),4)) }))
  write.csv(thr,"export/export_essd_thresholds.csv",row.names=FALSE); ok("export_essd_thresholds.csv")
}, error=function(e) cat("[skip] essd_thresholds:", conditionMessage(e), "\n"))

## ---- D) corrected+raw scores + covariates for grip models (fast) ----
CANDS <- c("bmi","bmi_c","physact","mets","exercise","vigact","modact","chronic","chronic_n",
           "n_chronic","comorbid","ncond","smoke","smoken","drink","agey")
export_scores <- function(dat, items, itemtype, outfile){
  dat <- dat %>% filter(!is.na(grip_max))
  if(nrow(dat)<200){ cat("[skip]",outfile,"(n<200)\n"); return(invisible()) }
  Xi <- dat[,paste0(items,"_d")]; colnames(Xi)<-items
  corr <- corrected_fscore(Xi, dat$female, itemtype)
  keep <- c("female","cesd","grip_max", intersect(CANDS, names(dat)))
  out <- cbind(dat[,keep,drop=FALSE], corr=corr)
  write.csv(out, outfile, row.names=FALSE); ok(basename(outfile))
}
tryCatch(export_scores(read.csv("./data/elsa_with_anchors.csv"), EN, "2PL",
                       "export/export_scores_elsa.csv"), error=function(e) cat("[skip] elsa scores:",conditionMessage(e),"\n"))
tryCatch({
  chl<-read.csv("./data/charls_w1to4.csv")
  if("grip_max"%in%names(chl)){
    wg<-chl%>%filter(!is.na(grip_max))%>%count(wave)%>%slice_max(n,n=1)%>%pull(wave)
    export_scores(chl%>%filter(wave==wg), CH, "graded", "export/export_scores_charls_grip.csv")
  }
}, error=function(e) cat("[skip] charls scores:",conditionMessage(e),"\n"))

## ---- E) missingness CC vs MI (optional; needs ./data/charls_w4_raw.csv) ----
tryCatch({
  if(file.exists("./data/charls_w4_raw.csv") && requireNamespace("mice",quietly=TRUE)){
    library(mice)
    dr<-read.csv("./data/charls_w4_raw.csv"); Xr<-dr[,paste0(ITEMS,"_d")]
    dl_from<-function(Xi,fem){ Xi<-as.data.frame(lapply(Xi,function(z)as.integer(round(z)))); colnames(Xi)<-ITEMS
      g<-factor(ifelse(fem==1,"F","M"),levels=c("M","F"))
      m<-multipleGroup(Xi,1,group=g,itemtype="graded",invariance=c(anchor,"free_means","free_var"),verbose=FALSE)
      coef(m,simplify=TRUE)$F$means[1] }
    nm<-rowSums(is.na(Xr)); cc<-dr[nm==0,]; dl_cc<-dl_from(cc[,paste0(ITEMS,"_d")],cc$female)
    imp<-mice(cbind(Xr,female=dr$female,raeducl=dr$raeducl,agey=dr$agey),m=5,method="pmm",printFlag=FALSE)
    dls<-sapply(1:5,function(k){ci<-complete(imp,k); dl_from(ci[,paste0(ITEMS,"_d")],ci$female)})
    write.csv(data.frame(method=c("complete_case",paste0("MI_",1:5)), delta_latent=c(dl_cc,dls)),
              "export/export_missingness_delta.csv", row.names=FALSE); ok("export_missingness_delta.csv")
    # also export the raw missingness matrix indicators for Python description
    md<-data.frame(female=dr$female, raeducl=dr$raeducl, agey=dr$agey,
                   n_missing=nm, obs_mean=rowMeans(Xr,na.rm=TRUE))
    write.csv(md, "export/export_missingness_rows.csv", row.names=FALSE); ok("export_missingness_rows.csv")
  } else cat("[skip] missingness: need ./data/charls_w4_raw.csv (+ mice). See 09_missingness.R header for the export template.\n")
}, error=function(e) cat("[skip] missingness:",conditionMessage(e),"\n"))

## ---- A) heavy: joint bootstrap draws with FIXED anchors (LAST) ----
tryCatch({
  cohend<-function(a,b){s<-sqrt(((length(a)-1)*var(a)+(length(b)-1)*var(b))/(length(a)+length(b)-2));(mean(a)-mean(b))/s}
  stat<-function(dd){
    Xi<-dd[,paste0(ITEMS,"_d")]; colnames(Xi)<-ITEMS; Xi[]<-lapply(Xi,function(z)as.integer(round(z)))
    g<-factor(ifelse(dd$female==1,"F","M"),levels=c("M","F"))
    draw<-cohend(dd$cesd[dd$female==1],dd$cesd[dd$female==0])
    m<-tryCatch(multipleGroup(Xi,1,group=g,itemtype="graded",invariance=c(anchor,"free_means","free_var"),verbose=FALSE),error=function(e)NULL)
    if(is.null(m)) return(c(NA,NA,NA,NA,NA))
    dlat<-coef(m,simplify=TRUE)$F$means[1]; fs<-fscores(m,method="EAP")[,1]
    dadj<-tryCatch(as.numeric(coef(lm(fs~dd$female+dd$raeducl+dd$agey))["dd$female"]),error=function(e)NA_real_)
    c(Delta_raw=draw,Delta_latent=dlat,Delta_adj=dadj,DDmeas=draw-dlat,DDstruct=dlat-dadj)
  }
  cat(sprintf("Running 2000 bootstrap resamples (fixed anchors: %s) ... ~15-25 min\n",paste(anchor,collapse=",")))
  bt<-boot(d,function(dat,idx) stat(dat[idx,]),R=2000)
  D<-as.data.frame(bt$t); names(D)<-c("Delta_raw","Delta_latent","Delta_adj","DDmeas","DDstruct")
  write.csv(D,"export/export_dd_draws.csv",row.names=FALSE); ok("export_dd_draws.csv")
}, error=function(e) cat("[skip] dd_draws:",conditionMessage(e),"\n"))

cat("\nAll requested exports were attempted; inspect the ./export/ folder.\n")
