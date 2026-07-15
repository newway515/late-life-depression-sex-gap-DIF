# =====================================================================
# 01_charls_dif_gap.R  (mirt-only 版, 不依赖 lordif/rms/dmacs)
# CHARLS 主分析(w4): 测量模型 + 性别 DIF(mirt::DIF) + 效应量(期望分数标准化差)
#   + 部分不变性 Δ_latent + Δ_adj + H2 bootstrap。
# Analysis specification: substantive DIF requires adj p<0.01 (BH) and ESSD>=0.20;
# the directional test uses one-sided alpha=0.025 and the bootstrap uses >=2000 draws.
# =====================================================================
suppressMessages({library(mirt); library(boot); library(dplyr)})
set.seed(20260709)
dir.create("output", showWarnings = FALSE)
ALPHA_DIF <- 0.01; ESSD_CUT <- 0.20; N_BOOT <- 2000
ITEMS <- c("depres","effort","sleepr","whappy","flone","going","bother","mindts","fhope","fear")

d <- read.csv("./data/charls_w4.csv")
X <- d[, paste0(ITEMS,"_d")]; colnames(X) <- ITEMS
X[] <- lapply(X, function(z) as.integer(round(z)))          # 0-3 有序, 抑郁方向
# 参照组 = 男(M), 故 factor 先 M 后 F => coef$F$means = Δ_latent(F-M) 直读
grp <- factor(ifelse(d$female==1,"F","M"), levels=c("M","F"))
NCAT <- max(X, na.rm=TRUE)                                    # 一般=3 (0..3)
dpar <- paste0("d", 1:NCAT)                                   # 阈值参数名 d1..d3

cohend <- function(a,b){s<-sqrt(((length(a)-1)*var(a)+(length(b)-1)*var(b))/(length(a)+length(b)-2)); (mean(a)-mean(b))/s}
delta_raw <- cohend(d$cesd[d$female==1], d$cesd[d$female==0])

## ---- 1) 维度检验: 单因子 vs 双因子(躯体 vs 情感) ----
grm1 <- mirt(X, 1, itemtype="graded", verbose=FALSE)
soma <- c("effort","sleepr","going")                          # prespecified somatic factor
spec <- mirt.model(sprintf("affect = %s\n soma = %s",
        paste(which(!(ITEMS %in% soma)), collapse=","),
        paste(which(ITEMS %in% soma), collapse=",")))
grm2 <- try(mirt(X, spec, itemtype="graded", verbose=FALSE), silent=TRUE)
cat("== 维度比较(单因子 vs 双因子) ==\n")
if(!inherits(grm2,"try-error")) print(anova(grm1, grm2)) else cat("双因子未收敛, 采用单因子\n")
cat("\n单因子拟合 M2:\n"); print(M2(grm1, type="C2"))

## ---- 2) 全不变基线 + mirt::DIF(逐条目 drop) ----
mg_full <- multipleGroup(X, 1, group=grp, itemtype="graded",
             invariance=c(ITEMS, "free_means","free_var"), verbose=FALSE)
dif <- DIF(mg_full, which.par=c("a1", dpar), scheme="drop",
           items2test=1:length(ITEMS), p.adjust="BH")
dif$item <- rownames(dif)
cat("\n== mirt::DIF(性别) ==\n"); print(dif[,c("item","X2","df","p","adj_p")])

## ---- 3) 效应量: 期望分数标准化差(ESSD) 每条目 ----
##   关键: 必须用 configural(自由)模型, 两组条目参数各自估计; 全约束模型会给全 0。
mgC <- multipleGroup(X, 1, group=grp, itemtype="graded", verbose=FALSE)   # configural, 各组 N(0,1)
Theta <- matrix(seq(-4,4,length.out=61)); w <- dnorm(Theta[,1]); w <- w/sum(w)
essd <- sapply(seq_along(ITEMS), function(i){
  eR <- expected.item(extract.item(mgC, i, group="M"), Theta, min=0)  # 男(参照)期望分曲线
  eF <- expected.item(extract.item(mgC, i, group="F"), Theta, min=0)  # 女(焦点)期望分曲线
  sdi <- sd(X[[ITEMS[i]]], na.rm=TRUE)                                # 观测条目 SD 做标准化(dMACS 式)
  as.numeric(sum(abs(eF - eR) * w) / ifelse(is.finite(sdi) && sdi > 0, sdi, 1))
})
names(essd) <- ITEMS
cat("\n== ESSD(期望分数标准化差; 越大=DIF 越强) ==\n"); print(round(essd,3))

## ---- 4) 实质性 DIF = 显著(adj_p<.01) 且 ESSD>=0.20 ----
sig <- dif$adj_p < ALPHA_DIF; names(sig) <- dif$item
substantive <- ITEMS[ ITEMS %in% names(which(sig)) & essd[ITEMS] >= ESSD_CUT ]
cat("\n实质性性别 DIF 条目:", paste(substantive, collapse=", "),
    " [H1 定向预测: sleepr, fear]\n")

## ---- 5) Δ_latent: 部分不变性(锚=非实质DIF条目) ----
anchor <- setdiff(ITEMS, substantive)
mg_part <- multipleGroup(X, 1, group=grp, itemtype="graded",
             invariance=c(anchor, "free_means","free_var"), verbose=FALSE)
delta_latent <- coef(mg_part, simplify=TRUE)$F$means[1]       # F 相对 M(=0)
cat(sprintf("\nΔ_raw=%.3f | Δ_latent=%.3f | ΔΔ_meas=%.3f\n",
            delta_raw, delta_latent, delta_raw-delta_latent))

## ---- 6) Δ_adj: 潜分 ~ 性别 + 教育 + 年龄 ----
fs <- fscores(mg_part, method="EAP")[,1]
adj <- lm(fs ~ d$female + d$raeducl + d$agey)
cat(sprintf("Δ_adj(校正教育/年龄后 female 系数)=%.3f\n", coef(adj)["d$female"]))

## ---- 7) H2 单侧检验(α=0.025): bootstrap Δ_latent ----
boot_dl <- function(dat, idx){
  dd <- dat[idx,]; Xi <- dd[,paste0(ITEMS,"_d")]; colnames(Xi)<-ITEMS
  Xi[] <- lapply(Xi, function(z) as.integer(round(z)))
  g <- factor(ifelse(dd$female==1,"F","M"), levels=c("M","F"))
  m <- try(multipleGroup(Xi,1,group=g,itemtype="graded",
        invariance=c(anchor,"free_means","free_var"), verbose=FALSE), silent=TRUE)
  if(inherits(m,"try-error")) return(NA_real_)
  coef(m, simplify=TRUE)$F$means[1]
}
cat(sprintf("\n运行 %d 次 bootstrap(可能十几分钟)...\n", N_BOOT))
bt <- boot(d, boot_dl, R=N_BOOT)
ci <- boot.ci(bt, type="perc", conf=0.95)$percent[4:5]
p_one <- mean(bt$t <= 0, na.rm=TRUE)
cat(sprintf("H2: Δ_latent=%.3f, 95%%CI[%.3f,%.3f], 单侧 p=%.4f (α=0.025)\n",
            delta_latent, ci[1], ci[2], p_one))
cat("H2 判定:", ifelse(p_one<0.025 & delta_latent>0 & abs(delta_latent)<abs(delta_raw),
                       "支持(缩小但不消失、方向不变)","不支持/看敏感性"), "\n")

saveRDS(list(delta_raw=delta_raw, delta_latent=delta_latent, delta_adj=coef(adj)["d$female"],
             anchor=anchor, substantive=substantive, dif=dif, essd=essd, boot=bt),
        "output/out_charls_w4_dif_gap.rds")
cat("\n已保存 output/out_charls_w4_dif_gap.rds\n")
