# =====================================================================
# 04_external_validity.R — 外部效度确证检验(预注册 §6.7)  [mirt-only 稳健版 v2]
# (A) HRS-CIDI: 原始 CES-D 分 vs DIF 校正潜分 对 cidi_mde3 的预测;
#     主检验 = 分数×性别交互是否因校正而"减弱"(减弱=条目含测量偏倚)。
# (B) 三库握力: grip_max ~ 抑郁分(原始 vs 校正) × 性别 + 年龄。
# 校正潜分 = 部分不变性(全约束→DIF筛[显著 且 ESSD>=0.20]→锚定非DIF→EAP)。
# =====================================================================
suppressMessages({library(mirt); library(pROC); library(dplyr)})
set.seed(20260709); dir.create("output", showWarnings=FALSE)
z <- function(x) as.numeric(scale(x))
ESSD_CUT <- 0.20; ALPHA_DIF <- 0.01

corrected_fscore <- function(X, female, itemtype){
  X <- as.data.frame(lapply(X, function(v) as.integer(round(v))))
  g <- factor(ifelse(female==1,"F","M"), levels=c("M","F"))
  ncat <- max(as.matrix(X), na.rm=TRUE)
  dpar <- if(itemtype=="graded") paste0("d",1:ncat) else "d"
  full <- multipleGroup(X,1,group=g,itemtype=itemtype,
            invariance=c(colnames(X),"free_means","free_var"), verbose=FALSE)  # 全约束(可识别)
  dif <- DIF(full, which.par=c("a1",dpar), scheme="drop",
             items2test=1:ncol(X), p.adjust="BH")
  adjp <- dif$adj_p[match(colnames(X), rownames(dif))]
  ## ESSD 效应量(configural 模型)
  mgC <- multipleGroup(X,1,group=g,itemtype=itemtype, verbose=FALSE)
  Th <- matrix(seq(-4,4,length.out=61)); w <- dnorm(Th[,1]); w <- w/sum(w)
  essd <- sapply(seq_len(ncol(X)), function(i){
    eR <- expected.item(extract.item(mgC,i,group="M"),Th,min=0)
    eF <- expected.item(extract.item(mgC,i,group="F"),Th,min=0)
    sdi <- sd(X[[i]],na.rm=TRUE); as.numeric(sum(abs(eF-eR)*w)/ifelse(is.finite(sdi)&&sdi>0,sdi,1))
  })
  flagged <- colnames(X)[ which(adjp < ALPHA_DIF & essd >= ESSD_CUT) ]
  anchor <- setdiff(colnames(X), flagged)
  if(length(anchor) < 2) anchor <- colnames(X)               # 安全兜底
  part <- multipleGroup(X,1,group=g,itemtype=itemtype,
            invariance=c(anchor,"free_means","free_var"), verbose=FALSE)
  list(fs=fscores(part,method="EAP")[,1], flagged=flagged, essd=setNames(round(essd,3),colnames(X)))
}

EN <- c("depres","effort","sleepr","whappy","flone","fsad","going","enlife")
CH <- c("depres","effort","sleepr","whappy","flone","going","bother","mindts","fhope","fear")

## ---------- (A) HRS-CIDI ----------
h <- read.csv("./data/hrs_with_cidi.csv") %>% filter(!is.na(cidi_mde3))
Xh <- h[,paste0(EN,"_d")]; colnames(Xh) <- EN
cf <- corrected_fscore(Xh, h$female, "2PL")
cat("HRS 实质性 DIF 条目:", paste(cf$flagged, collapse=", "), "\n")
cat("HRS ESSD:\n"); print(cf$essd)
h$corr <- z(cf$fs); h$raw <- z(h$cesd)
mA <- glm(cidi_mde3 ~ raw  * female, data=h, family=binomial)
mB <- glm(cidi_mde3 ~ corr * female, data=h, family=binomial)
cat("\n== HRS 原始分×性别 交互 ==\n"); print(round(summary(mA)$coef["raw:female",],4))
cat("== HRS 校正分×性别 交互 ==\n"); print(round(summary(mB)$coef["corr:female",],4))
cat("主检验 |交互(校正)| < |交互(原始)| ? ->",
    abs(coef(mB)["corr:female"]) < abs(coef(mA)["raw:female"]), "\n")
for(s in c(1,0)){ hs <- h[h$female==s,]
  aR <- as.numeric(auc(hs$cidi_mde3, predict(glm(cidi_mde3~raw ,hs,family=binomial),type="response"), quiet=TRUE))
  aC <- as.numeric(auc(hs$cidi_mde3, predict(glm(cidi_mde3~corr,hs,family=binomial),type="response"), quiet=TRUE))
  cat(sprintf("  sex=%s AUC 原始=%.3f 校正=%.3f\n", ifelse(s==1,"F","M"), aR, aC)) }

## ---------- (B) 握力 ----------
grip_one <- function(d, items, itemtype, tag){
  d <- d %>% filter(!is.na(grip_max))
  if(nrow(d) < 200){ cat(sprintf("\n== %s 握力: 有效样本不足(%d), 跳过 ==\n", tag, nrow(d))); return() }
  Xi <- d[,paste0(items,"_d")]; colnames(Xi) <- items
  d$corr <- z(corrected_fscore(Xi, d$female, itemtype)$fs); d$raw <- z(d$cesd)
  rA <- lm(grip_max ~ raw  * female + agey, d)
  rB <- lm(grip_max ~ corr * female + agey, d)
  cat(sprintf("\n== %s 握力: 原始分×性别 交互=%.3f | 校正分×性别 交互=%.3f (n=%d) ==\n",
      tag, coef(rA)["raw:female"], coef(rB)["corr:female"], nrow(d)))
}
grip_one(read.csv("./data/elsa_with_anchors.csv"), EN, "2PL", "ELSA")
# CHARLS: 握力仅 w1-w3, 主波 w4 无 -> 用 w1-4 长表里握力最全之波
chl <- read.csv("./data/charls_w1to4.csv")
if("grip_max" %in% names(chl)){
  wg <- chl %>% filter(!is.na(grip_max)) %>% count(wave) %>% slice_max(n,n=1) %>% pull(wave)
  grip_one(chl %>% filter(wave==wg), CH, "graded", sprintf("CHARLS(w%s)", wg))
} else cat("\nCHARLS 长表无 grip_max, 跳过握力\n")

cat("\n解释: 校正后交互减弱/两性一致性改善 -> 支持被标记条目含测量偏倚; 否则更倾向真实表型。\n")
saveRDS(list(hrs_mA=coef(mA), hrs_mB=coef(mB), hrs_flagged=cf$flagged), "output/out_external_validity.rds")
