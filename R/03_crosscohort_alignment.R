# =====================================================================
# 03_crosscohort_alignment.R — 三库共同 6 条目 alignment(cohort x sex)
# Binary 2PL model; six groups = CHARLS/ELSA/HRS x female/male;
# 输出各组潜均值 -> 分离(同 cohort 内 F-M = 性别 DIF 后差) 与(同性别跨 cohort = 文化 DIF);
# 对照 Δ_raw 梯度 China>England>US。sirt::invariance.alignment 与 mirt 互证。
# =====================================================================
suppressMessages({library(mirt); library(sirt); library(dplyr); library(tidyr)})
COMMON6 <- c("depres","effort","sleepr","whappy","flone","going")
d <- read.csv("./data/pooled_common6_binary.csv")
d$grp <- factor(paste(d$cohort, ifelse(d$female==1,"F","M"), sep="_"))
X <- d[, COMMON6]

## A) mirt 多组 2PL + DIF(相对锚), 估各组潜均值
## configural 模型: 每组各自 N(0,1)(可识别); alignment 从各组自由估计的条目参数出发
mg <- multipleGroup(X, 1, group=d$grp, itemtype="2PL", verbose=FALSE)
means <- sapply(coef(mg, simplify=TRUE), function(z) z$means[1])
cat("== configural 各组潜均值(均固定为 0, 正常) ==\n"); print(round(means,3))

## B) sirt alignment: 需要各组 item 参数(lambda, nu)。用多组 2PL 抽取后对齐。
##    简化流程: 直接用 sirt::invariance.alignment 处理各组 loading/intercept 矩阵。
pars <- coef(mg, simplify=TRUE)
lambda <- t(sapply(pars, function(z) z$items[,"a1"]))
nu     <- t(sapply(pars, function(z) z$items[,"d"]))
al <- invariance.alignment(lambda=lambda, nu=nu)
cat("\n== alignment 后各组 factor mean/SD ==\n"); print(round(al$pars,3))

## C) 分离性别 DIF 与文化 DIF: 从对齐后均值构造对比
am <- al$pars$alpha0; names(am) <- rownames(lambda)   # 对齐后各组均值
tab <- data.frame(grp=names(am), aligned_mean=as.numeric(am)) %>%
  separate(grp, c("cohort","sex"), sep="_")
sex_gap <- tab %>% group_by(cohort) %>%
  summarise(delta_latent_FM = aligned_mean[sex=="F"] - aligned_mean[sex=="M"])
cat("\n== 各 cohort alignment 后 Δ_latent(F-M) ==\n"); print(sex_gap)
cat("Reference raw-gap gradient: CHARLS 0.35 > ELSA 0.28 > HRS 0.22\n")

## D) 条目层面: 哪些条目非不变归因于 sex, 哪些归因于 cohort(文化/翻译)
cat("\n== alignment R2(越高=越不变) 与非不变条目见 al$itempars / al$es ==\n")
print(round(al$es, 3))
dir.create("output", showWarnings = FALSE)
saveRDS(list(means=means, alignment=al, sex_gap=sex_gap), "output/out_crosscohort.rds")
