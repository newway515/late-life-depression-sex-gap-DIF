# =====================================================================
# 01b_dimensionality.R — CES-D-10 维度检验(修复版, 用 bfactor 双因子)
# 目的: 补 01 里那个未收敛的双因子 spec。
# 关键修正:
#   旧写法 mirt(X, mirt.model("affect=.. soma=..")) 建的是"两个正交因子"
#   (因子间协方差默认=0), 躯体与情感高度相关时反而比单因子差 -> anova 异常。
#   正解: bfactor() 建"一般因子 + 特异因子"双因子结构(降维积分, 稳定收敛)。
# 判据: 双因子即便 LRT 显著更优, 若 ECV>=0.70 且 ωH>=0.70 -> 一般因子主导,
#       主分析用单因子有据(Reise 2013; Rodriguez et al. 2016)。
# =====================================================================
suppressMessages({library(mirt); library(dplyr)})
set.seed(20260709); dir.create("output", showWarnings = FALSE)
ITEMS <- c("depres","effort","sleepr","whappy","flone","going","bother","mindts","fhope","fear")
SOMA  <- c("effort","sleepr","going")          # 躯体/迟滞特异因子(预注册固定)

d <- read.csv("./data/charls_w4.csv")
X <- d[, paste0(ITEMS,"_d")]; colnames(X) <- ITEMS
X[] <- lapply(X, function(z) as.integer(round(z)))

## ---- 1) 单因子基线 ----
grm1 <- mirt(X, 1, itemtype="graded", verbose=FALSE)

## ---- 2) 双因子(bfactor): 一般因子(全条目) + 躯体/情感两个特异因子 ----
##   sfac: 每条目所属的"特异"因子(1=躯体, 2=情感); 一般因子由 bfactor 自动全载。
sfac <- ifelse(ITEMS %in% SOMA, 1L, 2L)
grm2 <- try(bfactor(X, sfac, itemtype="graded", verbose=FALSE,
                    technical=list(NCYCLES=3000)), silent=TRUE)

cat("========== CES-D-10 维度检验 ==========\n")
if(inherits(grm2,"try-error")){
  cat("bfactor 仍未收敛:\n"); print(grm2); quit(save="no")
}

## ---- 3) LRT + 拟合指数比较 ----
cat("\n== 单因子 vs 双因子 LRT(anova) ==\n"); print(anova(grm1, grm2))
cat("\n== 单因子 M2(C2) ==\n"); m2a <- M2(grm1, type="C2"); print(m2a)
cat("\n== 双因子 M2(C2) ==\n"); m2b <- try(M2(grm2, type="C2"), silent=TRUE)
if(!inherits(m2b,"try-error")) print(m2b) else cat("双因子 M2 跳过(可选)\n")

## ---- 4) 一般因子主导度: ECV + ωH(自标准化载荷) ----
Lmat <- summary(grm2, verbose=FALSE)$rotF          # 标准化载荷矩阵(含 NA=未载)
Lmat[is.na(Lmat)] <- 0
gcol <- which.max(colSums(abs(Lmat) > 1e-6))       # 稳健识别"一般因子"列(载荷最多)
lg   <- Lmat[, gcol]
spec <- setdiff(seq_len(ncol(Lmat)), gcol)

ecv <- sum(lg^2) / sum(Lmat^2)                     # Explained Common Variance
h2    <- rowSums(Lmat^2); theta <- 1 - h2           # 共同度 / 唯一性(标准化)
ss_sp <- sum(sapply(spec, function(j) sum(Lmat[,j])^2))
omega_h   <- sum(lg)^2 / (sum(lg)^2 + ss_sp + sum(theta))
omega_tot <- (sum(lg)^2 + ss_sp) / (sum(lg)^2 + ss_sp + sum(theta))

cat(sprintf("\n== 一般因子主导度 ==\nECV=%.3f | ωH=%.3f | ω_total=%.3f\n",
            ecv, omega_h, omega_tot))
cat("载荷矩阵(标准化):\n"); print(round(Lmat,3))

dominant <- ecv >= 0.70 && omega_h >= 0.70
cat(sprintf("\n判定: 一般因子%s主导 (ECV>=0.70 且 ωH>=0.70: %s)\n  -> 主分析%s采用单因子结构。\n",
    ifelse(dominant,"","未"), dominant, ifelse(dominant,"可","需重审是否")))

saveRDS(list(anova=anova(grm1,grm2), m2_1f=m2a,
             loadings=Lmat, ecv=ecv, omega_h=omega_h, omega_tot=omega_tot,
             dominant=dominant, soma=SOMA),
        "output/out_dimensionality.rds")
cat("\n已保存 output/out_dimensionality.rds\n")
