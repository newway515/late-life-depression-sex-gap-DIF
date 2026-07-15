# =====================================================================
# 05_sensitivity_matrix.R — sensitivity and robustness matrix.
# 逐个扰动设定重估 Δ_latent(CHARLS w4 为主),汇总方向与量级,画热图,给稳健判定。
# 稳健 = 方向全部为正(女>男) 且 跨设定极差 <= 0.10 个 d。
# =====================================================================
suppressMessages({library(mirt); library(dplyr); library(tidyr); library(ggplot2)})
ITEMS <- c("depres","effort","sleepr","whappy","flone","going","bother","mindts","fhope","fear")
COMMON6 <- c("depres","effort","sleepr","whappy","flone","going")
dir.create("output", showWarnings = FALSE)
d0 <- read.csv("./data/charls_w4.csv")
prev <- readRDS("output/out_charls_w4_dif_gap.rds")

est_delta <- function(d, items, anchor, itemtype="graded"){
  X <- d[,paste0(items,"_d")]; colnames(X)<-items
  X[] <- lapply(X, function(z) if(itemtype=="2PL") as.integer(z>=1) else as.integer(round(z)))
  g <- factor(ifelse(d$female==1,"F","M"), levels=c("M","F"))   # 男=参照(均值0)
  anchor <- intersect(anchor, items); if(length(anchor)<2) anchor <- items  # 锚题合法化
  m <- tryCatch(multipleGroup(X,1,group=g,itemtype=itemtype,
        invariance=c(anchor,"free_means","free_var"),verbose=FALSE), error=function(e) NULL)
  if(is.null(m)) return(NA); coef(m,simplify=TRUE)$F$means[1]
}

settings <- list()
settings[["Main (data-driven anchors, 10 items 0-3)"]] <- est_delta(d0, ITEMS, prev$anchor)
settings[["Full-invariance baseline"]]                 <- est_delta(d0, ITEMS, ITEMS)
settings[["Affective-only anchors"]]                   <- est_delta(d0, ITEMS, c("depres","flone","whappy","fhope"))
settings[["Leave-one-out (drop sleepr)"]]              <- est_delta(d0, ITEMS, setdiff(prev$anchor,"sleepr"))
settings[["Common 6 items (binary)"]]                  <- est_delta(d0, COMMON6, intersect(prev$anchor,COMMON6), "2PL")
settings[["w4 recompute (age>=60)"]]                   <- est_delta(read.csv("./data/charls_w1to4.csv") %>%
                                                 filter(wave==4), ITEMS, prev$anchor)  # exported data are age>=60 only
# wave perturbation (w1-w3)
for(w in 1:3){
  dw <- read.csv("./data/charls_w1to4.csv") %>% filter(wave==w)
  settings[[paste0("Wave ",w)]] <- est_delta(dw, ITEMS, prev$anchor)
}

S <- tibble(setting=names(settings), delta_latent=unlist(settings)) %>%
  mutate(delta_raw = prev$delta_raw)
cat("== Sensitivity matrix: latent sex gap (d) across settings ==\n"); print(as.data.frame(S))
rng <- diff(range(S$delta_latent, na.rm=TRUE))
robust <- all(S$delta_latent > 0, na.rm=TRUE) && rng <= 0.10
cat(sprintf("\nAll directions positive: %s | cross-setting range=%.3f (<=0.10) | verdict: %s\n",
    all(S$delta_latent>0,na.rm=TRUE), rng, ifelse(robust,"ROBUST","assumption-sensitive")))

## Heatmap
ggplot(S, aes(x="Latent sex gap", y=reorder(setting, delta_latent), fill=delta_latent)) +
  geom_tile() + geom_text(aes(label=sprintf("%.2f",delta_latent))) +
  scale_fill_gradient2(midpoint=0, low="#2166AC", mid="white", high="#B2182B",
                       name="Latent sex\ngap (d)") +
  labs(title="Robustness of the female-male depression gap (CHARLS w4)",
       subtitle="Latent-mean difference (female - male) across analytic settings",
       y=NULL, x=NULL) + theme_minimal(base_size=11)
ggsave("output/out_sensitivity_heatmap.png", width=7, height=5, dpi=150)
write.csv(S, "output/out_sensitivity_matrix.csv", row.names=FALSE)
