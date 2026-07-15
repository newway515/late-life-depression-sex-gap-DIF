# =====================================================================
# 10_grip_multicov.R — handgrip models with available additional covariates.
#   The grip-strength external-validity model adjusted only for age.
#   Re-fit the depression-score x sex interaction on grip with ADDITIONAL
#   covariates (whatever is available: BMI, physical activity, chronic-
#   disease count, smoking, ...), for the raw and DIF-corrected scores,
#   to show the correction-related attenuation is not an artefact of
#   residual confounding.
#   Self-contained: reuses the corrected_fscore() logic from 04.
# =====================================================================
suppressMessages({library(mirt); library(dplyr)})
set.seed(20260709)
z <- function(x) as.numeric(scale(x))
ESSD_CUT <- 0.20; ALPHA_DIF <- 0.01
EN <- c("depres","effort","sleepr","whappy","flone","fsad","going","enlife")   # ELSA/HRS 8-item
CH <- c("depres","effort","sleepr","whappy","flone","going","bother","mindts","fhope","fear")

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
  anchor <- setdiff(colnames(X), flagged); if(length(anchor) < 2) anchor <- colnames(X)
  part <- multipleGroup(X,1,group=g,itemtype=itemtype,
            invariance=c(anchor,"free_means","free_var"), verbose=FALSE)
  fscores(part, method="EAP")[,1]
}

# candidate extra covariates; only those present in the data are used
CANDS <- c("bmi","bmi_c","physact","mets","exercise","vigact","modact",
           "chronic","chronic_n","n_chronic","comorbid","ncond","smoke","smoken","drink")

grip_multi <- function(d, items, itemtype, tag){
  d <- d %>% filter(!is.na(grip_max), !is.na(agey))
  if(nrow(d) < 200){ cat(sprintf("\n== %s: n<200, skip ==\n", tag)); return(invisible()) }
  Xi <- d[,paste0(items,"_d")]; colnames(Xi) <- items
  d$corr <- z(corrected_fscore(Xi, d$female, itemtype)); d$raw <- z(d$cesd)
  extra <- intersect(CANDS, names(d))
  extra <- extra[sapply(extra, function(c) sum(!is.na(d[[c]])) > 0.8*nrow(d))]  # need decent coverage
  covs  <- paste(c("agey", extra), collapse=" + ")
  cat(sprintf("\n== %s grip (n=%d) | extra covariates used: %s ==\n",
              tag, nrow(d), if(length(extra)) paste(extra,collapse=", ") else "(none available - age only)"))
  fR <- lm(as.formula(paste0("grip_max ~ raw  * female + ", covs)), d)
  fC <- lm(as.formula(paste0("grip_max ~ corr * female + ", covs)), d)
  cat(sprintf("  raw x sex  interaction = %.3f\n", coef(fR)["raw:female"]))
  cat(sprintf("  corr x sex interaction = %.3f  (attenuation vs raw: %.1f%%)\n",
              coef(fC)["corr:female"],
              100*(1 - abs(coef(fC)["corr:female"])/abs(coef(fR)["raw:female"]))))
}

grip_multi(read.csv("./data/elsa_with_anchors.csv"), EN, "2PL", "ELSA")
chl <- read.csv("./data/charls_w1to4.csv")
if("grip_max" %in% names(chl)){
  wg <- chl %>% filter(!is.na(grip_max)) %>% count(wave) %>% slice_max(n,n=1) %>% pull(wave)
  grip_multi(chl %>% filter(wave==wg), CH, "graded", sprintf("CHARLS(w%s)", wg))
}
cat("\nInterpretation: if the corrected-vs-raw attenuation persists after adding available\n",
    "confounders, the sex-alignment of the depression-grip relation is not merely residual confounding.\n")
