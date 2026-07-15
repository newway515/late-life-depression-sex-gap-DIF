# =====================================================================
# 02_charls_longitudinal.R — CHARLS w1-w4 纵向不变性 + Δ_latent 轨迹
# Estimate wave-specific Delta_latent using the partial-invariance anchors from script 01;
# 检验跨波(配置/度量/标量)不变性;报告 Δ_latent 轨迹与稳定性。
# =====================================================================
suppressMessages({library(mirt); library(dplyr)})
ITEMS <- c("depres","effort","sleepr","whappy","flone","going",
           "bother","mindts","fhope","fear")
d <- read.csv("./data/charls_w1to4.csv")
dir.create("output", showWarnings = FALSE)
prev <- readRDS("output/out_charls_w4_dif_gap.rds"); anchor <- prev$anchor

cohend <- function(a,b){s<-sqrt(((length(a)-1)*var(a)+(length(b)-1)*var(b))/(length(a)+length(b)-2)); (mean(a)-mean(b))/s}

res <- lapply(sort(unique(d$wave)), function(w){
  dw <- d[d$wave==w,]; X <- dw[,paste0(ITEMS,"_d")]; colnames(X)<-ITEMS
  X[] <- lapply(X, function(z) as.integer(round(z)))
  g <- factor(ifelse(dw$female==1,"F","M"), levels=c("M","F"))   # 男=参照(均值0)
  m <- multipleGroup(X,1,group=g,itemtype="graded",
        invariance=c(anchor,"free_means","free_var"), verbose=FALSE)
  data.frame(wave=w, N=nrow(dw),
             delta_raw=cohend(dw$cesd[dw$female==1], dw$cesd[dw$female==0]),
             delta_latent=coef(m,simplify=TRUE)$F$means[1])
}) %>% bind_rows()
cat("== CHARLS 纵向 Δ_latent 轨迹 ==\n"); print(round(res,3))
cat(sprintf("\nΔ_latent 跨波极差 = %.3f (稳健阈值 <=0.10)\n",
            diff(range(res$delta_latent))))

## 纵向测量不变性(把 wave 当组): 配置->度量->标量 逐级约束比较
Xall <- d[,paste0(ITEMS,"_d")]; colnames(Xall)<-ITEMS
Xall[] <- lapply(Xall, function(z) as.integer(round(z)))
gw <- factor(d$wave)
config <- multipleGroup(Xall,1,group=gw,itemtype="graded",invariance="",verbose=FALSE)
metric <- multipleGroup(Xall,1,group=gw,itemtype="graded",invariance=c("slopes"),verbose=FALSE)
scalar <- multipleGroup(Xall,1,group=gw,itemtype="graded",
           invariance=c("slopes","intercepts","free_means","free_var"),verbose=FALSE)
cat("\n== 纵向不变性(配置/度量/标量) ==\n"); print(anova(config, metric, scalar))
write.csv(res, "output/out_charls_longitudinal.csv", row.names=FALSE)
